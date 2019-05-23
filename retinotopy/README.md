## Retinotopy Utilities

This directory contains utility scripts for solving the pRF models
from a retinotopy experiment using VistaSoft. The intended usage is to
copy all three of these scripts to your VistaSoft session's `Code`
directory (e.g.,
`Projects/Retinotopy/wl_subj042/20170810_ColorRetinotopy/Code`) then
to invoke them from there. The script may need to be edited if being
used outside the Winawer lab or if you want to analyze subsets of your
retinotopy data instead of just an average of all retinotopy scans.

For a clearer example of how to use these scripts, see [this demo
page](https://noahbenson.github.io/Retinotopy-Tutorial/).

### Author

Noah C. Benson [&lt;nben@nyu.edu&gt](mailto:nben@nyu.edu)

## Snakefile

A `Snakefile` is also presented here (for use with
[snakemake](https://snakemake.readthedocs.io/en/stable/)), which
should help automate retinotopy analyses for BIDS-formatted
directories. The following steps walks you through how to use it

### Setup

1. Download and install [miniconda](https://conda.io/miniconda.html)
   (this just contains python and conda, a very nifty package manager;
   choose python 3.7). Conda is separate from python: it's a package
   manager, which makes it easy to install a variety of python
   libraries. If you've ever used apt on Ubuntu or brew on Macs, then
   you've used a package manager before. Python has its own package
   manager, `pip`, which generally works very well, but in my
   experience I've found conda tends to work with fewer issues. [See
   here](https://stackoverflow.com/questions/20994716/what-is-the-difference-between-pip-and-conda)
   for some more details, but the gist seems to be: conda can handle
   external (non-python) dependencies, whereas pip cannot, and conda
   can create virtual environments (see item 3 in this list), whereas
   `pip` cannot (the standard python way is to use `virtualenv`, which
   also works well). [See
   here](https://jakevdp.github.io/blog/2016/08/25/conda-myths-and-misconceptions/)
   for a blog post from Jake VanderPlas with more details on conda.
2. Create the virtual environment for your retinotopy analyses by
   navigating to this directory and running
   
   ```
   conda env create -f environment.yml
   ```

	If you've already run this and just want to update the
    environment, run:
	
   ```
   conda env update -f environment.yml
   ```
3. Activate the virtual environment:

   ```
   conda activate retinotopy
   ```
4. Make sure you have matlab installed and on your path. Any version
   at or after 2016b should probably work, but I haven't tested this
   with many different versions.
5. Make sure you have
   [vistasoft](https://github.com/vistalab/vistasoft) downloaded and
   know where it is
6. Make sure you have the `winawerlab` and `Tesla` servers mounted and
   you know where they are (if you've already moved your retinotopy
   and freesurfer data off `Tesla`, you don't need it mounted).
7. Open `config.yml` from within this directory in your favorite text
   editor ([see here](https://gettaurus.org/docs/YAMLTutorial/) for a
   quick tutorial for yml/YAML, but it's pretty straightforward). Edit
   the paths here to be correct for your setup (the comments should
   explain what is what). Then edit the second half for the subjects
   and sessions you want to analyze. You can also edit `TASKS` and
   `NRUNS`, but if a `subject:session` pair doesn't show up there, we
   will assume it has `task-prf` and 6 runs, respectively; you only
   need to add `subject:session` if that's not the case.

### Run the analysis

Navigate to the directory containing the `Snakefile` and make sure the
`retinotopy` conda environment is active (type `conda activate
retinotopy` if you're not sure).

To see what steps are necessary to run the full analysis on any given
subject/session, type:

```
snakemake -n -rpk {BIDS_dir}/derivatives/vistasoft/{subject}/{session}/Outputs/lh.inferred_eccen.mgz
```

and replace `{BIDS_dir}` with the location of your BIDS directory
(e.g., `/mnt/winawerlab/Projects/Retinotopy/BIDS`), `{subject}` with
your subject number (e.g., `sub-wlsubj081`) and `{session}` with the
session number (e.g., `ses-01`).

(The `-n` option is a dry-run, `-r` will print out the reason for
running each step, `-p` will print out the command `snakemake` will
run, and `-k` will keep running independent jobs if one fails).

This will print out all the steps necessary to get from where you are
now to where you want to be and why it's running the steps
again. Double-check the steps to make sure everything makes sense
(e.g., it shouldn't be running the `move_off_tesla` step if you've
already moved your data). As a note, if the input to a step has been
updated, `snakemake` will want to run that step again.

If everything looks good, rerun that command without the `-n`:

```
snakemake -rpk {BIDS_dir}/derivatives/vistasoft/{subject}/{session}/Outputs/lh.inferred_eccen.mgz
```

You can also add `-j N`, where `N` is some number, in order to run N
jobs in parallel. By default, however, `snakemake` doesn't do anything
clever and so you may run use up all your memory or cores, in which
case it might throw an error. 

### On the cluster

If you wish to run this on the `prince` cluster, clone and setup the
[snakemake-slurm](https://github.com/billbrod/snakemake-slurm)
profile, then invoke the command:


```
snakemake --profile snakemake-slurm -u cluster.json -j 100 -rpk {BIDS_dir}/derivatives/vistasoft/{subject}/{session}/Outputs/lh.inferred_eccen.mgz
```

or some similarly large number for the number of jobs to run. This
will put logs for each of the jobs inside `{BIDS_dir}/code/{rule}`,
where `{rule}` is the name of the snakemake rule (and bot this naming
and the name of the individual log files should hopefully be pretty
intuitive).

### Run everything!

If you want to do a dry-run for all subjects and sessions at once, do
the following:

```
snakemake -n -rpk bayesian_retinotopy_all
```

(and similarly remove `-n` to do the actual thing, add `-j N` for
multiple jobs in parallel, make the necessary adjustments to run it on
prince, etc.)
