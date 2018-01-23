# MRI_tools
Matlab and Python tools for MRI analyses

# Requirements

Python requirements are listed in the requirements.txt file,
additionally [pydeface](https://github.com/poldracklab/pydeface)
version 2.0 is required to deface anatomical data (in
`prisma_to_BIDS.py`; if either pydeface or FSL is unavailable, we will
simply copy over the anatomical data), which also requires FSL, and
the preprocessing script requires FSL and Freesurfer.

# How to use

For an example of how to run, see the [Winawer lab
wiki](https://wikis.nyu.edu/pages/viewpage.action?pageId=86054639)
