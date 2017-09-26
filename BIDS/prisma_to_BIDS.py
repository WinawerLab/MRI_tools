#!/usr/bin/python
"""script to reorganize prisma data to BIDS format.

The main thing is a simple rearrangement of the data, changing the pathnames, as well as some
extraction of metadata into the required json.

this follows the general principle of "never alter data", so files will be copied instead of moved

For anatomical data, you'll need to check freesurfer recon-all.log (probably at
Acadia/Freesurfer_subjects/subj_num/scripts/recon-all.log) to see which anatomical images were used
in the recon-all call, and pass them as arguments. this will move them and deface them (using
Poldrack's lab pydeface, https://github.com/poldracklab/pydeface).

This is meant to be run one subject at a time, similar to the preprocessing script, also found in
this repo.
"""
import os
import glob
import shutil
import re
import subprocess

# to see if we have pydeface.py available. if this runs without
# hitting the except, it succeeded.
# try:
#     subprocess.call(["pydeface.py", 'in.nii', 'out.nii'])
# except OSError:
#     raise Exception("No pydeface.py! Download and install from https://github.com/poldracklab/pydeface")

# sbref should be in functional with _sbref instead of _func. PE are fmap, last version on BIDS
# specification.

# I assume xyz maps to ijk, which is what BIDS wants
def copy_func(data_dir, output_dir, epis, sbrefs, task_label, session_label=None, acq_label=None,
              rec_label=None, echo_index=None):
    """copy functional data to BIDS format

    looks for both .nii and .nii.gz files

    because we assume the data came off of prisma, we make pretty strong assumptions about its
    layout. We assume that datadir points to a directory that contains many folders, each of which
    is named "##+{str}", where ## is a two digit number, and {str} is some string that describes
    what's in the folder. We will not make use of {str}, except to find the subject name, which is
    expected to be at the beginning of the filename and of the format "wl_subj###". epis and sbrefs
    are lists of numbers, which indicate which of these directories contain the EPI and single-band
    reference images, respectively.

    you must have the same number of epis and sbrefs and it's assumed that they are paired. That
    is, if epis=[2,4,6] and sbrefs=[1,3,5], then directory 02+ contains the epi for run-01 and
    directory 01+ contains the sbref for run-01. The assumption does not depend on the order of
    number in the epis and sbrefs arguments, but instead in their numerical order.

    subject id will be inferred from the filenames. run index will be inferred from directory
    structure.

    output_dir: str, output directory. this should be the super directory. we will create the
    subject directory (if it does not already exist) here, and then the various sub-directories. if
    the subject directory already exists, session_label is required

    task_label: str, task label. str to identify the task, required by BIDS.

    session_label: str, session label, optional. if there is not subject directory in output_dir
    for this subject, then this is unnecessary. if there *is*, then a session_label is needed.

    acq_label, rec_label: strs, optional. further ways to identify this particular scanning
    session, added to functional data filename. see BIDS specification for details.
    
    echo_index: int, optional. another way to identify this particular scanning session, added to
    functional data filename. see BIDS specification for details.
    """
    nifti_str = [os.path.join(data_dir, "%02d+*", "*.nii"),
                 os.path.join(data_dir, "%02d+*", "*.nii.gz")]
    if len(epis) != len(sbrefs):
        raise Exception("You must have the same number of epis and sbrefs, because I assume "
                        "they're paired!")
    epi_files = [glob.glob(n % i) for n in nifti_str for i in epis]
    sbref_files = [glob.glob(n % i) for n in nifti_str for i in sbrefs]
    # this flattens out the lists and then sorts them, in case epis or sbrefs weren't sorted
    epi_files = sorted([i for sublist in epi_files for i in sublist])
    sbref_files = sorted([i for sublist in sbref_files for i in sublist])
    subj_id = re.search("(wl_subj[0-9]+)", epi_files[0]).groups()
    assert len(subj_id)==1, "Found more than one subject_id in path name! %s" % subj_id
    # we can't have any underscores in our subj_id, because BIDS uses them to separate fields.
    subj_id = subj_id[0].replace("_", "")
    subj_dir = os.path.join(output_dir, "sub-" + subj_id)
    filename_kw = {"subj_id": subj_id, "task_label": task_label}
    filename_kw.update(dict((k, "") for k in ["session", "acq", "rec", "echo"]))
    if os.path.exists(subj_dir) and session_label is None:
        raise Exception("subj_dir already exists but you have no session_label set! You need a "
                        "session_label if subject has multiple sessions")
    elif session_label is not None:
        subj_dir = os.path.join(subj_dir, "ses-%s" % session_label)
        filename_kw["session"] = "_ses-%s" % session_label
        if os.path.exists(subj_dir):
            raise Exception("The session_label %s has already been used, choose another!" % session_label)
    if acq_label is not None:
        filename_kw["acq"] = "_acq-%s" % acq_label
    if rec_label is not None:
        filename_kw["rec"] = "_rec-%s" % rec_label
    if echo_index is not None:
        filename_kw["echo"] = "_echo-%02d" % echo_index
    base_filename = "sub-{subj_id}{session}_task-{task_label}{acq}{rec}_run-%02d{echo}"
    base_filename = base_filename.format(**filename_kw)
    epi_filename = base_filename + "_bold.%s"
    sbref_filename = base_filename + "_sbref.%s"
    os.makedirs(subj_dir)
    if len(epi_files) + len(sbref_files) > 1:
        os.makedirs(os.path.join(subj_dir, "func"))
    for i, epi in enumerate(epi_files):
        if "gz" in epi:
            ext = "nii.gz"
        else:
            ext = "nii"
        # we want runs to be 1-indexed
        shutil.copy(epi, os.path.join(subj_dir, "func", epi_filename % (i+1, ext)))
    for i, sbref in enumerate(sbref_files):
        if "gz" in sbref:
            ext = "nii.gz"
        else:
            ext = "nii"
        # we want runs to be 1-indexed
        shutil.copy(sbref, os.path.join(subj_dir, "func", sbref_filename % (i+1, ext)))

def copy_fmap(data_dir, output_dir, distortPE, distortrevPE, PEdim='j'):
    """copy fieldmap data to BIDS format
    
    looks for both .nii and .nii.gz files

    because we assume the data came off of prisma, we make pretty strong assumptions about its
    layout. We assume that datadir points to a directory that contains many folders, each of which
    is named "##+{str}", where ## is a two digit number, and {str} is some string that describes
    what's in the folder. We will not make use of {str}, except to find the subject name, which is
    expected to be at the beginning of the filename and of the format "wl_subj###". distortPE and
    distortrevPE are single integers, which indicate which of these directories contain the phase
    encoding images, in the same and reverse directions, respectively, as the EPIs.

    PEdim: Phase encoding dimension. BIDS expects direction to use i, j, and k, but FSL uses x, y,
    and z. This function accepts both and will map between the two like so: x->i, y->j, z->k.
"""
