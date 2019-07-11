#! /bin/bash

# all necessary symbolic links to anatomy folder (there should be 7); these are the links:
# General (new project)

function die {
	echo "$*"
	exit 1
}

function dump_params {
	output_directory="$1"
	cat > "$output_directory"/parameters.json <<EOF
{
  "extract": [
    {
      "step": "classification",
      "repository": "cbinyu/heudiconv",
      "version": "v2.6",
      "options": "-c none --heuristic /data/CBI/heuristics/cbi_heuristic.py --bids --overwrite"
    },
    {
      "step": "conversion",
      "repository": "cbinyu/heudiconv",
      "version": "v2.6",
      "options": "-c dcm2niix --heuristic /data/CBI/heuristics/cbi_heuristic.py --bids --overwrite"
    },
    {
      "step": "validation",
      "repository": "bids/validator",
      "version": "1.2.1",
      "options": ""
    }
  ],
  "deface": [
    {
      "step": "deface",
      "repository": "cbinyu/bids_pydeface",
      "version": "v1.2",
      "options": "--n_cpus 5"
    }
  ],
  "fMRIPrep": [
    {
      "step": "participant",
      "repository": "poldracklab/fmriprep",
      "version": "1.3.1",
      "options": "participant --output-space T1w fsnative template --template-resampling-grid native --t2s-coreg --no-submm-recon"
    }
  ],
  "mriqc": [
    {
      "step": "participant",
      "repository": "cbinyu/mriqc",
      "version": "0.15.0",
      "options": "participant --verbose-reports --ica --fft-spikes-detector"
    },
    {
      "step": "group",
      "repository": "cbinyu/mriqc",
      "version": "0.15.0",
      "options": "group"
    }
  ]
}
EOF
}

project_path="$1"
[ -z "$project_path" ] && die "SYNTAX: bidsNewProject.sh <new project path>
Example:
   > bidsNewProject.sh /Volumes/CBIUserData/winawerlab/my_new_project
"

[ -d "$project_path" ] && die "Project path already exists!"

mkdir -p "$project_path" || die "Could not make project path: $project_path"
mkdir -p "$project_path"/derivatives/freesurfer || die "Could not make derivatives/freesurfer directory!"

dump_params "$project_path"/derivatives || die "Failed to create parameters.json file"





