#! /usr/bin/env python
####################################################################################################
# Script to post-process the nifti files output by the export_niftis.m script in this same
# directory. This script creates new volume files <dset>-angle.nii.gz and <dset>-eccen.nii.gz for
# each dataset found in the subject's Outputs directory; also resamples the relevant data to the
# subject's cortical surface.
# Author: Noah C. Benson <nben@nyu.edu>

import argparse, sys, os, six
import neuropythy as ny, numpy as np, nibabel as nib

def main(args):
    # Parse the arguments...
    parser = argparse.ArgumentParser()
    parser.add_argument('sub', metavar='subject', nargs=1,
                        help=('The FreeSurfer subject ID or directory.'))
    parser.add_argument('outdir', metavar='directory', type=str, nargs=1,
                        help='The directory containing the output nifti files.')
    parser.add_argument('-x', '--no-surf', required=False, default=True,
                        dest='surface', action='store_false',
                        help=('If provided, instructs the script not to produce files of the '
                              'pRF parameters resampled onto the cortical surface.'))
    parser.add_argument('-m', '--method', required=False, default='linear', dest='method',
                        help=('The method to use for volume-to-surface interpolation; this may'
                              ' be nearest or linear; the default is linear.'))
    parser.add_argument('-l', '--layer', required=False, default='midgray', dest='layer',
                        help=('Specifies the cortical layer to user in interpolation from volume'
                              ' to surface. By default, uses midgray. May be set to a value'
                              ' between 0 (white) and 1 (pial) to specify an intermediate surface'
                              ' or may be simply white, pial, or midgray.'))
    parser.add_argument('-d', '--subjects-dir', required=False, default=None, dest='sdir',
                        help=('Specifies the subjects directory to use; by default uses the'
                              ' environment variable SUBJECTS_DIR.'))
    parser.add_argument('-y', '--no-invert-y', required=False, default=True, dest='invert_y',
                        action='store_false',
                        help=('If provided, does not invert the y-coordinate exported from'
                              ' VistaSoft; by default this inversion is performed.'))
    parser.add_argument('-v', '--verbose', required=False, default=False, action='store_true',
                        dest='verbose', help='Print verbose output')
    if args[0].startswith('python'): args = args[2:]
    else: args = args[1:]
    args = parser.parse_args(args)
    # Check some of the arguments...
    sub = ny.freesurfer_subject(args.sub[0])
    dosurf = args.surface
    outdir = args.outdir[0]
    if not os.path.isdir(outdir):
        raise RuntimeError('Directory %s does not exist' % outdir)
    else:
        os.chdir(outdir)
    if args.verbose:
        def note(*args):
            six.print_(*args, flush=True)
            return True
    else:
        def note(*args):
            return False
    try: args.layer = float(args.layer)
    except: pass
    # figure out what datasets are here...
    dsets = [fl[:-13] for fl in os.listdir('.') if fl.endswith('-xcrds.nii.gz')]
    # loop over the given EPIs
    for ds in dsets:
        note('Processing Dataset %s...' % ds)
        # import the files...
        note('  - Importing parameters...')
        (x,y,s,v) = [ny.load('%s-%s.nii.gz' % (ds, suff), to='image')
                     for suff in ['xcrds','ycrds','sigma','vexpl']]
        # fix polar angle/eccen
        note('  - Creating polar angle/eccentricity images...')
        ang = np.arctan2(y.get_data(), x.get_data())
        ang = np.mod((90.0 - 180.0/np.pi * ang) + 180, 360) - 180
        a = nib.Nifti1Image(ang, x.affine, x.header)
        e = nib.Nifti1Image(np.sqrt(y.get_data()**2 + x.get_data()**2), x.affine, x.header)
        a.to_filename('%s-angle.nii.gz' % ds)
        e.to_filename('%s-eccen.nii.gz' % ds)
        # surfaces...
        if not dosurf: continue
        note('  - Projecting to surface...')
        to_output = {'xcrds':x, 'ycrds':y, 'angle':a, 'eccen':e, 'sigma': s, 'vexpl': v}
        for (dname,im) in six.iteritems(to_output):
            (ldat, rdat) = sub.image_to_cortex(im, args.layer,
                                               method=args.method, dtype=np.float32,
                                               weight=v)
            # we need to fix the dimensions...
            for (d,h) in zip([ldat,rdat], ['lh','rh']):
                ny.save('%s.%s-%s.mgz' % (h, ds, dname), d)
    # That's it!
    return 0

main(sys.argv)
