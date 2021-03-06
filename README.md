# MRI_tools
Matlab and Python tools for MRI analyses

# Requirements

Python 3 (at least 3.6, but probably others) is required.

Python requirements are listed in the requirements.txt file,
additionally [pydeface](https://github.com/poldracklab/pydeface)
version 2.0 is required to deface anatomical data (in
`prisma_to_BIDS.py`; if either pydeface or FSL is unavailable, we will
simply copy over the anatomical data), which also requires FSL, and
the preprocessing script requires FSL (version 5.0.10, 6.0 doesn't
seem to work) and Freesurfer.

For `bidsGLM.m` and its associated scripts, you will need
[GLMdenoise](https://github.com/kendrickkay/GLMdenoise) to run the GLM
and, optionally, [vistasoft](https://github.com/vistalab/vistasoft) if
you want to run GLMdenoise on `.mgz` files (freesurfer surface files);
for niftis, we use the built-in `niftiread` function.

# How to use

## Preprocessing

For an example of how to run, see the [Winawer lab
wiki](https://wikis.nyu.edu/pages/viewpage.action?pageId=86054639)

## Prisma to BIDS

These functions can be used to transfer data from the default way it
comes off NYU CBI's prisma scanner to a BIDS-compliant structure. This
is only a temporary way of handling things, since eventually CBI will
add BIDS as a possible export option. For an example of how to use
this,
see
[this script](https://github.com/billbrod/spatial-frequency-preferences/blob/master/sfp/transfer_to_BIDS.py) from
the spatial frequency preferences experiment.
