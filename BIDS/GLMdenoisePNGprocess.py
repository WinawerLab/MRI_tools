#! /usr/bin/env python
# Script to crawl a directory of GLMdenoise outputs and convert the surface-png files into plots of
# the results on the cortical surface.

import sys, six, warnings, os, PIL

def die(*args):
    if len(args) > 1:
        sys.stderr.write(args[0] % tuple(args[1:]))
        sys.stderr.write('\n')
        sys.stderr.flush()
    elif len(args) > 0:
        sys.stderr.write(args[0])
        sys.stderr.write('\n')
        sys.stderr.flush()
    sys.exit(1)

try: import numpy as np
except: die('Could not import numpy library')
try: import pimms
except: die('Could not import pimms library')
try: import neuropythy as ny
except: die('Could not import neuropythy library')
try: import matplotlib, matplotlib.pyplot as plt
except: die('Could not import matplotlib library')


# okay, go to the directory requested
if len(sys.argv) < 2: die('SYNTAX: GLMdenoisePNGprocess.py <directory>')
pth = sys.argv[1]
if not os.path.isdir(pth): die('Directory %s not found', pth)
pth = os.path.abspath(pth)
os.chdir(pth)

# Deduce the subject and derivatives folder
(pp, f) = os.path.split(pth)
(sub, derpth) = (None, None)
while derpth is None or sub is None:
    (pp,f) = os.path.split(pp)
    if len(f) == 0: break
    elif f.startswith('sub-'): sub = f[4:]
    elif f == 'derivatives': derpth = os.path.join(pp, f)
if derpth is None: die('Could not find derivatives path')
if sub is None: die('Could not deduce subject from BIDS path')
if os.path.isdir(os.path.join(derpth, 'freesurfer')):
    sub = os.path.join(derpth, 'freesurfer', sub)
elif 'SUBJECTS_DIR' in os.environ:
    tmp = os.path.join(os.environ['SUBJECTS_DIR'], sub)
    if os.path.isdir(tmp): sub = tmp
else: die('No derivatives/freesurfer or $SUBJECTS_DIR/%s directory found', sub)

try:
    sub = ny.freesurfer_subject(sub)
    lh = sub.lh
    rh = sub.rh
except: die('Failed to load freesurfer subject')

print('FreeSurfer Path: %s' % sub)
print('Image Path:      %s' % pth)

sys.stdout.write('\nPreparing flat-maps...')
sys.stdout.flush()
try:
    mps_post = {h: ny.map_projection('occipital_pole', h, radius=np.pi/2)
                for h in ['lh','rh']}
    mps_ante = {h: mp.copy(center=-mp.center, center_right=-mp.center_right)
                for (h,mp) in six.iteritems(mps_post)}
    fms_post = {h: mp(sub.hemis[h]) for (h,mp) in six.iteritems(mps_post)}
    fms_ante = {h: mp(sub.hemis[h]) for (h,mp) in six.iteritems(mps_ante)}
    sys.stdout.write(' Done.')
finally:
    sys.stdout.write('\n')
    sys.stdout.flush()

# find all the png files
fls = [fl for fl in os.listdir('.') if fl.lower().endswith('.png')]
# for each, process them:
(nl, nr) = (lh.vertex_count, rh.vertex_count)
nn = nl + nr + 1 # for some reason these are padded...
print(nl, nr, nn)
for fl in fls:
    outfl = fl[:-4] + '_maps.png'
    sys.stdout.write('Processing image %s into %s...' % (fl, outfl))
    sys.stdout.flush()
    try:
        im = PIL.Image.open(fl)
        im = np.asarray(im.convert('F'))
        if im.shape[0] != 2 or im.shape[1] != nn:
            sys.stdout.write('Wrong size: (%d x %d)' % im.shape[:2])
            continue
        mx = np.max(im[0])
        im = im[0] / mx
        (l,r) = (im[:nl], im[nl:])
        # make the plots:
        (fig,axs) = plt.subplots(2, 2, figsize=(7.5,7.5), dpi=72*4)
        for (u,h,axcol) in zip([l,r], ['lh','rh'], axs.T):
            for (fm,ax) in zip([fms_post[h], fms_ante[h]], axcol):
                ny.cortex_plot(fm, color=u[fm.labels], cmap='afmhot',
                               alpha=1, mask=(u[fm.labels] > 0),
                               axes=ax)
                ax.axis('off')
        plt.subplots_adjust(0,0,1,1,0,0)
        plt.savefig(outfl)
        plt.close(fig)
    except: raise
    else: sys.stdout.write(' Done.')
    finally:
        sys.stdout.write('\n')
        sys.stdout.flush()

print('%d PNG files sucessfully processed.' % len(fls))
sys.exit(0)





