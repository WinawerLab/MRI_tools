#! /usr/bin/env python

# Script to reorient volumes into anatomical orientation and to create surface time-series.
# This script is made for use with a registration file (in FreeSurfer's tkreg format) and
# a series of unwarped time-series volumes: i.e., the output of Serra's preprocessing
# script.
# Author: Noah C. Benson <nben@nyu.edu>

import argparse, sys, os, six
import neuropythy as ny, numpy as np, nibabel as nib

def main(args):
    # Parse the arguments...
    parser = argparse.ArgumentParser()
    parser.add_argument('reg', metavar='registration_file', nargs=1,
                        help=('The distort2anat_tkreg.dat or similar file: the registration'
                              ' file, in FreeSurfer\'s tkreg format, to apply to the EPIs.'))
    parser.add_argument('epis', metavar='EPI', type=str, nargs='+',
                        help='The EPI files to be converted to anatomical orientation')
    parser.add_argument('-t', '--tag', required=False, default='_anat', dest='tag', nargs=1,
                        help=('A tag to append to the output filenames; if given as - then'
                              ' overwrites original files. By default, this is "_anat".'))
    parser.add_argument('-s', '--surf', required=False, default=False,
                        dest='surface', action='store_true',
                        help=('If provided, instructs the script to also produce files of the '
                              'time-series resampled on the cortical surface.'))
    parser.add_argument('-o', '--out', required=False, default='.', dest='outdir',
                        help=('The output directory to which the files should be written; by'
                              ' default this is the current directory (.); note that if this'
                              ' directory also contains the EPI files and there is no tag given,'
                              ' then the EPIs will be overwritten.'))
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
    parser.add_argument('-v', '--verbose', required=False, default=False, action='store_true',
                        dest='verbose', help='Print verbose output')
    if args[0].startswith('python'): args = args[2:]
    else: args = args[1:]
    args = parser.parse_args(args)
    # Check some of the arguments...
    epis = args.epis
    if len(epis) < 1: raise RuntimeError('No EPIs given')
    tag = args.tag
    if tag == '-': tag = ''
    dosurf = args.surface
    outdir = args.outdir
    if not os.path.isdir(outdir):
        raise RuntimeError('Directory %s does not exist' % outdir)
    if args.verbose:
        def note(*args):
            six.print_(*args, flush=True)
            return True
    else:
        def note(*args):
            return False
    # Read in the registration file
    args.reg = args.reg[0]
    if not os.path.isfile(args.reg):
        raise RuntimeError('Given registration file not found: %s' % args.reg)
    with open(args.reg, 'r') as f:
        lines = []
        while True:
            s = f.readline()
            if s is None or s == '': break
            lines.append(s)
    # This tells us some info...
    sub = lines[0].strip()
    if args.sdir is not None:
        ny.add_subject_path(args.sdir)
    try: sub = ny.freesurfer_subject(sub)
    except: raise ValueError('No subject %s; you may need to set your SUBJECTS_DIR' % sub)
    affine = np.asarray([[float(ss) for ss in s.split()] for s in lines[4:8]])
    affinv = np.linalg.inv(affine)
    displm = np.dot(sub.voxel_to_native_matrix, np.linalg.inv(sub.voxel_to_vertex_matrix))
    # loop over the given EPIs
    for epi in epis:
        note('Processing EPI %s...' % epi)
        # import the epi file..
        img = ny.load(epi)
        # edit the header...
        note('   - Correcting volume orientation...')
        new_affine = np.dot(displm, np.dot(affinv, ny.freesurfer.tkr_vox2ras(img)))
        newimg = nib.Nifti1Image(img.dataobj, new_affine, img.header)
        (epi_dir,epi_flnm) = os.path.split(epi)
        srf_flnm = (epi_flnm[:-3] if epi_flnm[:-4] in ['.mgz', '.mgh', '.nii'] else epi_flnm[:-6])
        srf_flnm += 'mgz'
        newimg.to_filename(os.path.join(args.outdir, epi_flnm))
        # okay, now project to the surface
        if args.surface:
            note('   - Projecting to surface...')
            (ldat, rdat) = sub.image_to_cortex(newimg, args.layer, method=args.method)
            # we need to fix the dimensions...
            for (d,h) in zip([ldat,rdat], ['lh','rh']):
                im = nib.freesurfer.mghformat.MGHImage(
                    np.transpose(reduce(np.expand_dims, [-1,-1], d), (0,2,3,1)),
                    np.eye(4))
                im.to_filename(os.path.join(args.outdir, h + '.' + srf_flnm))
    # That's it!
    return 0

main(sys.argv)
