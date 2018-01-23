#!/usr/bin/python
"""script to reorganize prisma data to BIDS format.

The main thing is a simple rearrangement of the data, changing the pathnames, as well as some
extraction of metadata into the required json.

This cannot create the events.tsv or participants.tsv files, so you'll have to do this yourself.

This follows the general principle of "never alter data", so files will be copied instead of moved

For anatomical data, you'll need to check freesurfer recon-all.log (probably at
Acadia/Freesurfer_subjects/{subj_num}/scripts/recon-all.log) to see which anatomical images were
used in the recon-all call, and pass them as arguments. this will move them and, if possible,
deface them (using Poldrack's lab pydeface, https://github.com/poldracklab/pydeface).

This is meant to be run one subject at a time, similar to the preprocessing script, also found in
this repo.

For field map directions, I assume xyz maps to ijk, which is what BIDS wants
"""

import os
import glob
import shutil
import re
import subprocess
import json
from collections import Counter
import nibabel as nib
import warnings


def _construct_base_filename(example_file, output_dir, data_type, run_flag=False, task_label=None,
                             session_label=None, acq_label=None, rec_label=None, echo_index=None,
                             dir_flag=False, ce_label=None):
    """helper function to construct base filename, using various elements

    data_type: "func", "fmap", "anat", etc. Will be part of the subj_dir that gets returned

    run_flag: boolean, default False. Whether to add _run-%02d to the end of the base filename or
    not.
    """
    subj_id = re.search("(wl_subj[0-9]+)", example_file).groups()
    assert len(subj_id) == 1, "Found more than one subject_id in path name! %s" % subj_id
    # we can't have any underscores in our subj_id, because BIDS uses them to separate fields.
    subj_id = subj_id[0].replace("_", "")
    subj_dir = os.path.join(output_dir, "sub-" + subj_id)
    filename_kw = {"subj_id": subj_id}
    filename_kw.update(dict((k, "") for k in ["task", "session", "acq", "rec", "echo", "run",
                                              "dir", "ce"]))
    if session_label is None:
        subj_dir = os.path.join(subj_dir, data_type)
        if os.path.exists(os.path.join(subj_dir)):
            raise IOError("subj_dir already exists but you have no session_label set! You need a "
                          "session_label if subject has multiple sessions")
    else:
        subj_dir = os.path.join(subj_dir, "ses-%s" % session_label, data_type)
        filename_kw["session"] = "_ses-%s" % session_label
        if os.path.exists(subj_dir):
            raise IOError("The session_label %s has already been used, choose another!" % session_label)
    if task_label is not None:
        filename_kw['task'] = "_task-%s" % task_label
    if acq_label is not None:
        filename_kw["acq"] = "_acq-%s" % acq_label
    if rec_label is not None:
        filename_kw["rec"] = "_rec-%s" % rec_label
    if echo_index is not None:
        filename_kw["echo"] = "_echo-%02d" % echo_index
    if run_flag:
        filename_kw['run'] = "_run-%02d"
    if dir_flag:
        filename_kw['dir'] = "_dir-%s"
    if ce_label is not None:
        filename_kw['ce'] = '_ce-%s' % ce_label
    base_filename = "sub-{subj_id}{session}{task}{acq}{ce}{rec}{dir}{run}{echo}"
    return base_filename.format(**filename_kw), subj_dir


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
    if not hasattr(epis, '__iter__'):
        epis = [epis]
    if not hasattr(sbrefs, '__iter__'):
        sbrefs = [sbrefs]
    if len(epis) != len(sbrefs):
        raise Exception("You must have the same number of epis and sbrefs, because I assume "
                        "they're paired!")
    epi_files = [glob.glob(n % i) for n in nifti_str for i in epis]
    sbref_files = [glob.glob(n % i) for n in nifti_str for i in sbrefs]
    # this flattens out the lists and then sorts them, in case epis or sbrefs weren't sorted
    epi_files = sorted([i for sublist in epi_files for i in sublist])
    sbref_files = sorted([i for sublist in sbref_files for i in sublist])
    base_filename, subj_dir = _construct_base_filename(
        epi_files[0], output_dir, "func", True, re.sub('[\W_]+', '', task_label), session_label,
        acq_label, rec_label, echo_index)
    epi_filename = base_filename + "_bold.%s"
    sbref_filename = base_filename + "_sbref.%s"
    os.makedirs(subj_dir)
    for i, (epi, sbref) in enumerate(zip(epi_files, sbref_files)):
        if "gz" in epi:
            epi_ext = "nii.gz"
        else:
            epi_ext = "nii"
        # we want runs to be 1-indexed
        shutil.copy(epi, os.path.join(subj_dir, epi_filename % (i+1, epi_ext)))
        if "gz" in sbref:
            sbref_ext = "nii.gz"
        else:
            sbref_ext = "nii"
        # we want runs to be 1-indexed
        shutil.copy(sbref, os.path.join(subj_dir, sbref_filename % (i+1, sbref_ext)))
        nii = nib.load(epi)
        tr = nii.header['pixdim'][4]
        t_units = nii.header.get_xyzt_units()[1]
        if t_units == 'msec':
            tr /= 1000
        elif t_units != 'sec':
            raise Exception("Don't know how to handle units %s for TR" % t_units)
        bold_dict = {"TaskName": task_label, 'RepetitionTime': float(tr)}
        with open(os.path.join(subj_dir, epi_filename % (i+1, 'json')), 'w') as f:
            json.dump(bold_dict, f)


def copy_fmap(data_dir, output_dir, distortPE, distortrevPE, PEdim='j', session_label=None,
              acq_label=None, run_label=None):
    """copy fieldmap data to BIDS format

    looks for both .nii and .nii.gz files

    because we assume the data came off of prisma, we make pretty strong assumptions about its
    layout. We assume that datadir points to a directory that contains many folders, each of which
    is named "##+{str}", where ## is a two digit number, and {str} is some string that describes
    what's in the folder. We will not make use of {str}, except to find the subject name, which is
    expected to be at the beginning of the filename and of the format "wl_subj###". distortPE and
    distortrevPE are single integers, which indicate which of these directories contain the phase
    encoding images, in the same and reverse directions, respectively, as the EPIs.

    we also make the following other assumptions:

    - that this is the BIDS fieldmap case 4 (section 8.9.4 in the specifications pdf): multiple
      phase encoded directions (topup).

    - (following the lead of the prisma_preproc script,) that the readout times are 1 second for
      each volume in each distortion scan and that there are 3 volumes per scan, for a total
      readout time of 3 seconds.

    - that these are inteded for all scans gathered in the same session.

    - the direction label can be inferred from the filename, which will have DISTORTION_%%, where
      %% are two alphanumeric characters that we'll use for the direction label (note that these
      will not be used to infer the phase encoding direction, they're only for readability)

    PEdim: Phase encoding dimension. BIDS expects direction to use i, j, and k, but FSL uses x, y,
    and z. This function accepts both and will map between the two like so: x->i, y->j, z->k.

    session_label: str, session label, optional. if there is not subject directory in output_dir
    for this subject, then this is unnecessary. if there *is*, then a session_label is needed.

    acq_label, run_label: strs, optional. further ways to identify this particular
    scanning session, added to functional data filename. see BIDS specification for details.
    """
    nifti_str = [os.path.join(data_dir, "%02d+*", "*.nii"),
                 os.path.join(data_dir, "%02d+*", "*.nii.gz")]
    distortScan = [glob.glob(n % distortPE) for n in nifti_str]
    # this flattens these out, so it's a list instead of a list of lists
    distortScan = [i for sublist in distortScan for i in sublist]
    distortRevScan = [glob.glob(n % distortrevPE) for n in nifti_str]
    # this flattens these out, so it's a list instead of a list of lists
    distortRevScan = [i for sublist in distortRevScan for i in sublist]
    assert len(distortScan) == 1, "Found more than one forward distortion scan!"
    assert len(distortRevScan) == 1, "Found more than one reverse distortion scan!"
    distortScan = distortScan[0]
    distortRevScan = distortRevScan[0]
    base_filename, subj_dir = _construct_base_filename(distortScan, output_dir, "fmap",
                                                       run_label is not None,
                                                       session_label=session_label,
                                                       acq_label=acq_label, dir_flag=True)
    if run_label is not None:
        # we want to get the run label in there, but leave the direction label for filling in later
        base_filename = base_filename % (run_label, '%s')
    base_filename += "_epi.%s"
    os.makedirs(subj_dir)
    for i, distort_file in enumerate([distortScan, distortRevScan]):
        if "gz" in distort_file:
            ext = "nii.gz"
        else:
            ext = "nii"
        direction = re.search("DISTORTION_([A-Za-z0-9]+)", distort_file).groups()[0]
        shutil.copy(distort_file, os.path.join(subj_dir, base_filename % (direction, ext)))
        phase_dir = {'x': 'i', 'y': 'j', 'z': 'k'}.get(PEdim, PEdim)
        # the reversed distortion scan should have the dash after it.
        if i == 1:
            phase_dir += '-'
        fmap_dict = {'PhaseEncodingDirection': phase_dir, 'TotalReadoutTime': 3}
        with open(os.path.join(subj_dir, base_filename % (direction, "json")), 'w') as f:
            json.dump(fmap_dict, f)


def copy_anat(data_dir, output_dir, anat_nums, modality_label, session_label=None, acq_label=None,
              ce_label=None):
    """copy anatomical data

    we make much less assumptions about the layout of this data.

    anat_paths: string or list of strings. If many strings, then they're all assumed to have the
    same of every label supplied and will be differentiated using run index.
    """
    if modality_label not in ['T1w', 'T2w', 'T1rho', 'T1map', 'T2map', 'T2map', 'T2star', 'FLAIR',
                              'FLASH', 'PD', 'PDmap', 'PDT2', 'inplaneT1', 'inplaneT2', 'angio',
                              'defacemask', 'SWImageandphase']:
        raise Exception("modality_label %s is not a valid BIDS modality!" % modality_label)
    nifti_str = [os.path.join(data_dir, "*%02d+*", "*.nii"),
                 os.path.join(data_dir, "*%02d+*", "*.nii.gz")]
    if not hasattr(anat_nums, '__iter__'):
        anat_nums = [anat_nums]
        run_label = None
    elif len(anat_nums) > 1:
        run_label = range(1, len(anat_nums)+1)
    anat_files = [glob.glob(n % a) for n in nifti_str for a in anat_nums]
    anat_files = [i for sublist in anat_files for i in sublist]
    base_filename, subj_dir = _construct_base_filename(
        anat_files[0], output_dir, "anat", run_label is not None, session_label=session_label,
        acq_label=acq_label)
    base_filename += '_' + modality_label + '.%s'
    os.makedirs(subj_dir)
    for i, anat_file in enumerate(anat_files):
        if 'gz' in anat_file:
            ext = 'nii.gz'
        else:
            ext = 'nii'
        if run_label is None:
            output_path = os.path.join(subj_dir, base_filename % ext)
        else:
            output_path = os.path.join(subj_dir, base_filename % (run_label[i], ext))
        out_code = 0
        try:
            out_code = subprocess.call(['pydeface', '--outfile', output_path, anat_file])
        except OSError:
            warnings.warn("No pydeface (version 2.0 is required)! Download and install from "
                          "https://github.com/poldracklab/pydeface We are simply copying over the "
                          "anatomical data.")
            shutil.copy(anat_file, output_path)
        if out_code == 1:
            raise IOError("Cannot find anatomical file at %s!" % anat_file)
        elif out_code == 2:
            warnings.warn("pydeface encountered an error (see above). We'll simply copy over the"
                          " anatomical data.")
            shutil.copy(anat_file, output_path)


def json_check(dir_to_check):
    """check whether we can consolidate the jsons in a directory

    this goes through all the jsons in the specified directory and determines what keys we can
    consolidate. following the BIDS inheritance principle, this means we create a json in the
    directory further up containing those keys and its assumed that everything in the
    sub-directories has those same values. for example, if all the differnt run*.json files in our
    func directory have the same RepetitionTime (which is likely), this will remove that key from
    all those files and put it in a json file without the run indicator on directory up
    """
    json_files = dict((f, None) for f in glob.glob(os.path.join(dir_to_check, "*.json")))
    for jf in json_files.iterkeys():
        with open(jf) as f:
            json_files[jf] = json.load(f)
    shared_dict = {}
    all_keys = {k for jd in json_files.itervalues() for k in jd.iterkeys()}
    for k in all_keys:
        try:
            all_vals = [jd[k] for jd in json_files.itervalues()]
            most_common_val = Counter(all_vals).most_common()[0][0]
            shared_dict[k] = most_common_val
        except KeyError:
            pass
    for jf, jd in json_files.iteritems():
        for k in shared_dict.iterkeys():
            if shared_dict[k] == jd[k]:
                jd.pop(k)
        os.remove(jf)
        if any(jd.keys()):
            with open(jf, 'w') as f:
                json.dump(jd, f)
    sup_dir = os.path.dirname(os.path.dirname(json_files.keys()[0]))
    name_parts = os.path.split(json_files.keys()[0])[-1].split("_")
    for k in json_files.keys()[1:]:
        name_parts = [n for n in name_parts if n in os.path.split(k)[-1].split("_")]
    with open(os.path.join(sup_dir, "_".join(name_parts)), 'w') as f:
        json.dump(shared_dict, f)
