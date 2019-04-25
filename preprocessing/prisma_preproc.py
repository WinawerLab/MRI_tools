import argparse
import sys
import os
import os.path as op
from glob import glob
import numpy as np
from nipype import Workflow, Node, MapNode, DataSink
from nipype.interfaces import fsl, freesurfer as fs
import json
import warnings
try:
    from bids import BIDSLayout
except IOError:
    warnings.warn("Can't find pybbids, so won't be able to pre-process BIDS data")


def _get_BIDS_name(layout, value, datadir, target=None):
    if value == 'sub':
        found_value_tmp = layout.get_subjects()
    elif value == 'ses':
        found_value_tmp = layout.get_sessions()
    elif value == 'task':
        found_value_tmp = layout.get_tasks()
    found_value = []
    for v in found_value_tmp:
        if "%s-" % value not in v:
            found_value.append("%s-%s" % (value, v))
        else:
            found_value.append(v)
    if len(found_value) != 1:
        if len(found_value) == 0:
            raise Exception("Found no %s name from data directory %s!" % (value, datadir))
        if len(found_value) > 1:
            if target is not None and target in found_value:
                found_value = [target]
            elif value != 'task':
                # we allow multiple tasks, but not multiple subjects or sessions
                raise Exception("Found more than 1 %s name from data directory %s! %s" %
                                (value, datadir, found_value))
    return found_value


def main(arglist):
    """Preprocess NYU CBI Prisma data"""

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-subject',
                        help=('Freesurfer subject id. Note that we use the $SUBJECTS_DIR '
                              'environmental variable to find the required data. Must be set if '
                              'dir_structure==\'prisma\'. if dir_structure==\'bids\', will use '
                              'the inferred BIDS subject id for the Freesurfer subject id *unless*'
                              ' this is set, in which case we will use this instead'))
    parser.add_argument('-datadir', required=True, help='Raw MR data path')
    parser.add_argument('-outdir', help='Output directory path. ignored if dir_structure==\'bids\'')
    parser.add_argument('-epis', nargs='+', type=int,
                        help=('EPI scan numbers, optional. If set (and you\'re using bids data), '
                              'you must only have one task or use the bids_task input argument.'))
    parser.add_argument('-sbref', type=int,
                        help=('Single band reference scan number. Only required if dir_structure is'
                              ' prisma (if bids, we determine it from filenames)'))
    parser.add_argument('-distortPE', type=int,
                        help=('Distortion scan number with same PE as epis. Only required if '
                              'dir_structure is prisma (if bids, we determine it from the '
                              'filenames).'))
    parser.add_argument('-distortrevPE', type=int,
                        help=('Distortion scan number with reverse PE as epis. Only required if '
                              'dir_structure is prisma (if bids, we determine it from the '
                              'filenames).'))
    parser.add_argument('-PEdim', type=str, default='y',
                        help=('PE dimension (x, y, or z). Only necessary if dir_structure is '
                              'prisma (if bids, we determine it from metadata)'))
    parser.add_argument("-plugin", type=str, default="MultiProc",
                        help=("Nipype plugin to use for running this. MultiProc (default) is "
                              "normally fine for when running locally, though it may use up all "
                              "your  computer's resources. Linear is slower, but won't do that."
                              "SLURM should be used on NYU HPC prince cluster."))
    parser.add_argument('-working_dir', default=None,
                        help=("Path to your working directory. By default, this will be within "
                              "your output directory, but you may want to place it elsewhere. For "
                              "example, on the HPC cluster, you may run out of space if this is in"
                              "your /home directory, so you probably want this in /scratch"))
    parser.add_argument('-dir_structure', default='prisma',
                        help=("{prisma, bids}. Is your data directory structured like it just came"
                              " off the prisma scanner ('prisma') or is it BIDS structured "
                              "('bids')? This determines how we look for the various scans. If "
                              "your data is BIDS-structured, then datadir should point to the "
                              "particular session you want to preprocess. Outputs will then be: "
                              "{BIDS_dir}/derivatives/{bids_derivative_name}/{BIDS_subject_name}/"
                              "{BIDS_session_name}/{BIDS_subject_name}_{BIDS_session_name}_"
                              "{BIDS_task_name}_{run}_{bids_suffix}.nii.gz, where BIDS_subject_"
                              "name, BIDS_session_name, and BIDS_task_name are all inferred from "
                              "input data and BIDS_dir is two directories above datadir (since "
                              "datadir corresponds to one BIDS session). We also preprocess tasks"
                              " separately; if you have more than one task, use the -bids_task flag"
                              " to specify which one you want to preprocess"))
    parser.add_argument('-bids_derivative_name', default='preprocessed',
                        help=("the name of the derivatives directory in which to place the "
                              "outputs. ignored if dir_structure=='prisma'"))
    parser.add_argument('-bids_suffix', default='preproc',
                        help=("the suffix to place at the end of the output filenames. ignored if"
                              " dir_structure=='prisma'"))
    parser.add_argument('-bids_task', default=None,
                        help=("Which bids task to preprocess. Only required if you have more than"
                              " one task in your session directory (we only preprocess one task at"
                              " a time). If you only have one task, we can grab its name from the"
                              " filenames."))
    parser.add_argument('-plugin_args', default=None,
                        help=("Any additional arguments to pass to nipype's workflow.run as plugin"
                              "_args. A single entry in the resulting dictionary should be of the"
                              " format arg:val (e.g., n_procs:2) with multiple args separated by a"
                              " comma with no spaces (e.g., n_procs:2,memory_gb:5). see nipype's "
                              "plugin documentation for more details on possible values: "
                              "http://nipype.readthedocs.io/en/latest/users/plugins.html. "))
    args = vars(parser.parse_args(arglist))

    try:
        subj_dir = os.environ['SUBJECTS_DIR']
        if not os.path.isdir(subj_dir):
            raise Exception("Unable to find your Freesurfer SUBJECTS_DIR! Is it on an external "
                            "drive that you have not mounted?")
    except KeyError:
        raise Exception("You do not have your Freesurfer SUBJECTS_DIR environmental variable set!")

    # Session paths and files
    session = dict()
    session['data'] = args['datadir']
    if args['dir_structure'] == 'prisma':
        session['Freesurfer_subject_name'] = args['subject']
        session['nii_temp'] = op.join(session['data'], '%02d+*', '*.nii')
        session['epis'] = [glob(session['nii_temp'] % r)[0] for r in args['epis']]
        session['sbref'] = glob(session['nii_temp'] % args['sbref'])[0]
        session['distort_PE'] = glob(session['nii_temp'] % args['distortPE'])[0]
        session['distort_revPE'] = glob(session['nii_temp'] % args['distortrevPE'])[0]
        session['PE_dim'] = args['PEdim']
        session['out'] = args['outdir']
        session['out_names'] = ["timeseries_corrected_run-%02d.nii.gz" % i
                                for i in range(1, len(session['epis']) + 1)]
    elif args['dir_structure'] == 'bids':
        layout = BIDSLayout(session['data'], validate=False)
        # _gets_BIDS_name returns a list, but we've already checked (for sub and ses) that it only
        # contains one value, so we just grab that
        session['BIDS_subject_name'] = _get_BIDS_name(layout, 'sub', session['data'])[0]
        session['BIDS_session_name'] = _get_BIDS_name(layout, 'ses', session['data'])[0]
        # this will be a list of task names (though in most cases, it will only contain one value)
        session['BIDS_task_names'] = _get_BIDS_name(layout, 'task', session['data'],
                                                    args['bids_task'])
        if args['subject'] is None:
            # we do NOT want "sub-" in the Freesurfer subject name
            session['Freesurfer_subject_name'] = session['BIDS_subject_name'].replace('sub-', '')
        else:
            session['Freesurfer_subject_name'] = args['subject']
            
        test_files = layout.get('tuple', extensions=['nii', 'nii.gz'], suffix='bold',
                                task=[i.replace('task-', '') for i in session["BIDS_task_names"]])
        if args['epis'] is not None:
            # then we assume that args['epis'] gives us the run numbers we want
            run_nums = args['epis']
            assert len(session['BIDS_task_names']) == 1, "If epis argument is set, can only analyze one task!"
            task_names = [session['BIDS_task_names'][0].replace('task-', '') for i in run_nums]
        else:
            run_nums = [i.run for i in test_files]
            task_names = [i.task for i in test_files]
        for i, n in zip(run_nums, task_names):
            if len(layout.get('file', extensions=['nii', 'nii.gz'], suffix='bold', run=i,
                              task=n)) != 1:
                raise Exception("Found zero or multiple bold nifti files with run %s, task %s! We "
                                "require that there only be 1." % i, n)
        session['epis'] = layout.get('file', extensions=['nii', 'nii.gz'], suffix='bold', run=run_nums,
                                     task=task_names)
        if len(session['epis']) == 0:
            raise Exception("Unable to find any epis! This probably means your pybids version isn'"
                            "t new enough; it must be version >= 0.7!")
        session['sbref'] = layout.get('file', extensions=['nii', 'nii.gz'], suffix='sbref')[0]
        distortion_scans = layout.get('file', extensions=['nii', 'nii.gz'], suffix='epi')
        distortion_PEdirections = {}
        for scan in distortion_scans:
            distortion_PEdirections[layout.get_metadata(scan)['PhaseEncodingDirection']] = scan
        epi_PEdirections = []
        for scan in session['epis']:
            epi_PEdirections.append(layout.get_metadata(scan)['PhaseEncodingDirection'])
        if len(set(epi_PEdirections)) != 1:
            raise Exception("Cannot find unique phase encoding direction for your functional data"
                            " in data directory %s!" % session['data'])
        # we want PEdim to be x, y, or z, but coming from BIDS jsons it will be one of i, j, k
        PE_conversion_dict = {'i': 'x', 'j': 'y', 'k': 'z', 'i-': 'x-', 'j-': 'y-', 'k-': 'z-'}
        session['PE_dim'] = PE_conversion_dict[epi_PEdirections[0]]
        session['distort_PE'] = distortion_PEdirections[epi_PEdirections[0]]
        if '-' not in epi_PEdirections[0]:
            session['distort_revPE'] = distortion_PEdirections["%s-" % epi_PEdirections[0]]
        else:
            session['distort_revPE'] = distortion_PEdirections[epi_PEdirections[0].replace('-', '')]
        session['bids_derivative_name'] = args['bids_derivative_name']
        session['bids_suffix'] = args['bids_suffix']
        out_dir = op.abspath(op.join(session['data'], "../..",
                                    ("derivatives/{bids_derivative_name}/{BIDS_subject_name}/"
                                     "{BIDS_session_name}/")))
        session['out'] = out_dir.format(**session)
        # we replace the bids suffix _bold with _preproc (or whatever we set bids_suffix to
        # be). this way, 'bold' can occur elsewhere in the filename without causing any issues.
        session['out_names'] = [op.split(i)[1].replace('_bold', '_'+session['bids_suffix'])
                                for i in session['epis']]
        # MAKE SURE THIS IS GZIPPED
        session['out_names'] = [i + '.gz' if not i.endswith('gz') else i
                                for i in session['out_names']]
    else:
        raise Exception("Don't know what to do with dir_structure %s!" % args['dir_structure'])

    if not os.path.isdir(os.path.join(subj_dir, session['Freesurfer_subject_name'])):
        raise Exception("Unable to find Freesurfer subject dir for subject %s. Maybe you should "
                        "use the -subject flag to specify an alternate name?" %
                        (session['Freesurfer_subject_name']))

    if not op.exists(session['out']):
        os.makedirs(session['out'])

    if args['working_dir'] is not None:
        session['working_dir'] = args['working_dir']
    else:
        session['working_dir'] = session['out']
    if not op.exists(session["working_dir"]):
        os.makedirs(session['working_dir'])

    session['plugin_args'] = {}
    if args['plugin_args'] is not None:
        for arg in args['plugin_args'].split(','):
            if len(arg.split(':')) != 2:
                raise Exception("Your plugin_args is incorrectly formatted, each should contain one colon!")
            k, v = arg.split(':')
            try:
                session['plugin_args'][k] = int(v)
            except ValueError:
                try:
                    session['plugin_args'][k] = float(v)
                except ValueError:
                    session['plugin_args'][k] = v
    session['plugin'] = args['plugin']

    # Dump session info to json
    with open(op.join(session['out'], 'session.json'), 'w') as sess_file:
        json.dump(session, sess_file, sort_keys=True, indent=4)

    # Define preprocessing worklow
    preproc_wf = create_preproc_workflow(session)

    # Execute workflow in parallel
    preproc_wf.run(session["plugin"], plugin_args=session['plugin_args'])


def create_preproc_workflow(session):
    """
    Defines simple functional preprocessing workflow, including motion
    correction, registration to distortion scans, and unwarping. Assumes
    recon-all has been performed on T1, and computes but does not apply
    registration to anatomy.
    """

    #---Create workflow---
    wf = Workflow(name='workflow', base_dir=session['working_dir'])


    #---EPI Realignment---

    # Realign every TR in each functional run to the sbref image using mcflirt
    realign = MapNode(fsl.MCFLIRT(ref_file=session['sbref'],
                                  save_mats=True,
                                  save_plots=True),
                      iterfield=['in_file'], name='realign')
    realign.inputs.in_file = session['epis']
    wf.add_nodes([realign])


    #---Registration to distortion scan---

    # Register the sbref scan to the distortion scan with the same PE using flirt
    reg2dist = Node(fsl.FLIRT(in_file=session['sbref'],
                              reference=session['distort_PE'],
                              out_file='sbref_reg.nii.gz',
                              out_matrix_file='sbref2dist.mat',
                              dof=6),
                    name='reg2distort')
    wf.add_nodes([reg2dist])


    #---Distortion correction---

    # Merge the two distortion scans for unwarping
    distort_scans = [session['distort_PE'], session['distort_revPE']]
    merge_dist = Node(fsl.Merge(in_files=distort_scans,
                                dimension='t',
                                merged_file='distortion_merged.nii.gz'),
                      name='merge_distort')
    wf.add_nodes([merge_dist])

    # Run topup to estimate warpfield and create unwarped distortion scans
    if '-' not in session['PE_dim']:
        PEs = np.repeat([session['PE_dim'], session['PE_dim'] + '-'], 3)
    else:
        PEs = np.repeat([session['PE_dim'], session['PE_dim'].replace('-', '')], 3)
    unwarp_dist = Node(fsl.TOPUP(encoding_direction=list(PEs),
                                 readout_times=[1, 1, 1, 1, 1, 1],
                                 config='b02b0.cnf', fwhm=0),
                       name='unwarp_distort')
    wf.connect(merge_dist, 'merged_file', unwarp_dist, 'in_file')

    # Unwarp sbref image in case it's useful
    unwarp_sbref = Node(fsl.ApplyTOPUP(in_index=[1], method='jac'),
                        name='unwarp_sbref')
    wf.connect([(reg2dist, unwarp_sbref,
                 [('out_file', 'in_files')]),
                (unwarp_dist, unwarp_sbref,
                    [('out_enc_file', 'encoding_file'),
                     ('out_fieldcoef', 'in_topup_fieldcoef'),
                     ('out_movpar', 'in_topup_movpar')])])


    #---Registration to anatomy---

    # Create mean unwarped distortion scan
    mean_unwarped_dist = Node(fsl.MeanImage(dimension='T'),
                              name='mean_unwarped_distort')
    wf.connect(unwarp_dist, 'out_corrected', mean_unwarped_dist, 'in_file')

    # Register mean unwarped distortion scan to anatomy using bbregister
    reg2anat = Node(fs.BBRegister(subject_id=session['Freesurfer_subject_name'],
                                  contrast_type='t2',
                                  init='fsl',
                                  out_reg_file='distort2anat_tkreg.dat',
                                  out_fsl_file='distort2anat_flirt.mat'),
                    name='reg2anat')
    wf.connect(mean_unwarped_dist, 'out_file', reg2anat, 'source_file')


    #---Combine and apply transforms to EPIs---

    # Split EPI runs into 3D files
    split_epis = MapNode(fsl.Split(dimension='t'),
                         iterfield=['in_file'], name='split_epis')
    split_epis.inputs.in_file = session['epis']
    wf.add_nodes([split_epis])

    # Combine the rigid transforms to be applied to each EPI volume
    concat_rigids = MapNode(fsl.ConvertXFM(concat_xfm=True),
                            iterfield=['in_file'],
                            nested=True,
                            name='concat_rigids')
    wf.connect([(realign, concat_rigids,
                 [('mat_file', 'in_file')]),
                (reg2dist, concat_rigids,
                    [('out_matrix_file', 'in_file2')])])

    # Apply rigid transforms and warpfield to each EPI volume
    correct_epis = MapNode(fsl.ApplyWarp(interp='spline', relwarp=True),
                           iterfield=['in_file', 'ref_file', 'premat'],
                           nested=True,
                           name='correct_epis')

    get_warp = lambda warpfields: warpfields[0]
    wf.connect([(split_epis, correct_epis,
                    [('out_files', 'in_file'),
                     ('out_files', 'ref_file')]),
                (concat_rigids, correct_epis,
                    [('out_file', 'premat')]),
                (unwarp_dist, correct_epis,
                    [(('out_warps', get_warp), 'field_file')])])

    # Merge processed files back into 4D nifti
    merge_epis = MapNode(fsl.Merge(dimension='t',
                                   merged_file='timeseries_corrected.nii.gz'),
                         iterfield='in_files',
                         name='merge_epis')
    wf.connect([(correct_epis, merge_epis, [('out_file', 'in_files')])])


    #---Copy important files to main directory---
    substitutions = [('_merge_epis%d/timeseries_corrected.nii.gz' % i, n)
                     for i, n in enumerate(session['out_names'])]
    ds = Node(DataSink(base_directory=os.path.abspath(session['out']),
                       substitutions=substitutions),
              name='outfiles')
    wf.connect(unwarp_dist, 'out_corrected', ds, '@unwarp_dist')
    wf.connect(mean_unwarped_dist, 'out_file', ds, '@mean_unwarped_dist')
    wf.connect(unwarp_sbref, 'out_corrected', ds, '@unwarp_sbref')
    wf.connect(reg2anat, 'out_reg_file', ds, '@reg2anat')
    wf.connect(merge_epis, 'merged_file', ds, '@merge_epis')

    return wf


if __name__ == '__main__':
    main(sys.argv[1:])
