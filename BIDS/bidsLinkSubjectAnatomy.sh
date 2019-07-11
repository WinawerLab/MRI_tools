#! /bin/bash

function die {
	echo "$*"
	exit 1
}
function syntax {
	die "SYNTAX: bidsLinkSubjectAnatomy.sh <new project path>
Examples:
   > bidsLinkSubjectAnatomy.sh bert /Volumes/CBIUserData/winawerlab/my_new_project
   > bidsLinkSubjectAnatomy.sh bert #(uses pwd as project path)
"
}

subject="$1"
[ -z "$subject" ] && syntax


project_path="$2"
[ -z "$project_path" ] && project_path="`pwd`"

[ -d "$project_path"/"sub-${subject}" ] && die "Project directory already exists!"

[ -d "$project_path"/../Anatomy/derivatives/fmriprep/"sub-${subject}" ] || die "Subject not found in Anatomy directory!"
[ -d "$project_path"/../Anatomy/derivatives/freesurfer/"sub-${subject}" ] || die "Subject not found in Anatomy directory!"
[ -d "$project_path"/../Anatomy/"sub-${subject}" ] || die "Subject not found in Anatomy directory!"

if ! [ -a "$project_path"/derivatives/freesurfer/fsaverage ]
then ln -s ../../../Anatomy/derivatives/freesurfer/fsaverage "$project_path"/derivatives/freesurfer/
fi

# make fmriprep directory
[ -d "$project_path"/derivatives/fmriprep/"sub-${subject}" ] || {
	mkdir -p "$project_path"/derivatives/fmriprep/"sub-${subject}" || dir "Could not make fmriprep/sub-<subject> directory!"
}

[ -d "$project_path"/"sub-${subject}" ] || {
	mkdir -p "$project_path"/"sub-${subject}" || dir "Could not make sub-<subject> directory!"
}


# make links
ln -s ../../../../Anatomy/derivatives/fmriprep/"sub-${subject}"/anat "$project_path"/derivatives/fmriprep/"sub-${subject}"/
ln -s ../../../../Anatomy/derivatives/fmriprep/"sub-${subject}"/ses-anat "$project_path"/derivatives/fmriprep/"sub-${subject}"/
ln -s ../../Anatomy/"sub-${subject}"/ses-anat "$project_path"/"sub-${subject}"/
ln -s ../../../Anatomy/derivatives/freesurfer/"sub-${subject}" "$project_path"/derivatives/freesurfer/


