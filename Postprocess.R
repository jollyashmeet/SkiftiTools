#!/usr/bin/env Rscript
# Copyright 2023-2024 Harri Merisaari <haanme@utu.fi>
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Program and also corresponding docker version
version<-"0.2.3"

suppressPackageStartupMessages(library("argparse"))
library(argparse)
library(RNifti)
library(stringr)
source('Nifti2Skifti.R')
source('writeSkifti.R')
source('writeCSV.R')
# Excluded in docker
#source('Skifti_view.R')
#packageVersion("rgl")

# Function to write screen output 
write_print<-function(entry, verbose=0) {
  if(verbose>0) {
    print(paste(as.POSIXlt(Sys.time()), " ", entry, sep=""))
  }
}

# Function to write a log entry
write_log<-function(odir, entry, verbose=0) {
  if(verbose>0) {
    print(paste(as.POSIXlt(Sys.time()), " ", entry, sep=""))
  }
  write(paste(as.POSIXlt(Sys.time()), " ", entry, sep=""), file=paste(odir, "log.txt", sep="/"), append = TRUE)
}

# Function to write selected indexes of CSV to output
write_csv<-function(selected_names, selected_indexes, output_filename) {
  file.create(output_filename, showWarnings = FALSE)
  write(paste("Name", "Selection_index", sep='\t'), file=output_filename, append = TRUE)
  # write and resolve subjectnames at the same time from the file
  for(i in 1:length(selected_indexes)) {
    write(paste(selected_names[i], selected_indexes[i], sep='\t'), file=output_filename, append = TRUE)
  }
} 

parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=1, help="Print extra output 0=None, 1=Verbose, 2=Debug [default: 1]")
parser$add_argument("--path", default="", help = "Path to apply postprocessing [default WORKING DIRECTORY]")
parser$add_argument("--outputpath", default="", help = "Path to save postprocessing output [default WORKING DIRECTORY/postprocess_output]")
parser$add_argument("--TBSSsubfolder", default="out", help = "TBSS subfolder with stats relative to path [default out]")
parser$add_argument("--maskfile", default="/enigmaDTI/ENIGMA_DTI_FA_skeleton_mask.nii.gz", help = "Nifti mask file for obtained data [default /enigmaDTI/ENIGMA_DTI_FA_skeleton_mask.nii.gz]")
parser$add_argument("--name", default="UNKNOWN", help = "Name identifying the data [default UNKNOWN]")
parser$add_argument("--writecsv", default="No", help = "Yes/No Write output as simple CSV file [default:No]")
parser$add_argument("--writepng", default="No", help = "Yes/No Write mean intensity as png [default:No] (not available with docker)")
parser$add_argument("--writemaskcoordinates", default="No", help = "Yes/No Write mask coordinates in ASCII [default:No]")
parser$add_argument("--filetype", default="ASCII", help = "Output file type ASCII tab separated or rDs R archive file [default:ASCII]")
parser$add_argument("--compression", default="none", help = "Compression for Skifti data bz2/zip/none [default:none]")
parser$add_argument("--subjectsfile", default="", help = "Text file for subjectnames, one name in each line [default: use running numbers]")
parser$add_argument("--selectionfile", default="", help = "Text file for selected subjectnames, one name in each line [default: use all]")
parser$add_argument("--labelfile", default="", help = "Nifti file for skifti subregions, only masked region is considered [default:empty (no label file)]")
parser$add_argument("--scalars", default="FA,MD", help = "Comma separated list of scalars from FA,MD,AD,RD [default FA,MD]")

args <- parser$parse_args()

## These are for debugging with RStudio
#args$verbose<-1
#args$path<-'/Users/haanme/LIFESPAN_DTI/lifespan_test_data/src-LIFESPAN-2024-02-DTI-tools-V2-FSL-single-shell-ols'
#args$path<-'/Users/haanme/LIFESPAN_DTI/tbss_model_data/tbss_tryout_healthy_adults'
#args$outputpath<-'/Users/haanme/LIFESPAN_DTI/lifespan_test_data/src-LIFESPAN-2024-02-DTI-tools-V2-FSL-single-shell-ols/postprocess_output'
#args$outputpath<-'/Users/haanme/LIFESPAN_DTI/tbss_model_data/tbss_tryout_healthy_adults/out_ants_tbss_postprocess_output'
#args$TBSSsubfolder<-'out_ants_tbss_enigma_ss'
#args$TBSSsubfolder<-'out'
#args$maskfile<-'/Users/haanme/LIFESPAN_DTI/lifespan_dti/docker/enigmaDTI/ENIGMA_DTI_FA_skeleton_mask.nii.gz'
#args$name<-'TESTNAME'
#args$writecsv<-'Yes'
#args$writepng<-'No'
#args$writemaskcoordinates<-'No'
#args$filetype<-'ASCII'
#args$compression<-'zip'
#args$subjectsfile<-'/Users/haanme/LIFESPAN_DTI/lifespan_test_data/src-LIFESPAN-2024-02-DTI-tools-V2-FSL-single-shell-ols/CASELIST.txt'
#args$subjectsfile<-'/Users/haanme/LIFESPAN_DTI/tbss_model_data/tbss_tryout_FA_FB_infant/CASELIST.txt'
#args$selectionfile<-'/Users/haanme/LIFESPAN_DTI/tbss_model_data/tbss_tryout_FA_FB_infant/CASELIST.txt'
#args$labelfile<-''
#args$scalars<-'FA'

# Start log entries
start_time<-as.POSIXlt(Sys.time())
if ( nchar(args$outputpath) > 0) {
  odir<-args$outputpath
} else {
  odir<-paste(getwd(), '/postprocess_output', sep='')
}
dir.create(odir, showWarnings = FALSE)
invisible(file.create(paste(odir, "log.txt", sep="/"), showWarnings = FALSE))
write_log(odir, paste("Starting LIFESPAN Postprocess version:", version, sep=""), args$verbose)
write_log(odir, paste("Using ", odir, " as base directory for output", sep=''), args$verbose-1)

# List of zipped files
files_to_zip<-c()
# List of files to be removed after processing
files_to_be_removed<-c()
files_to_be_removed<-c(files_to_be_removed, paste(odir, "log.txt", sep="/"))
known_args<-data.frame(parser$parse_known_args()[[1]])
write_log(odir, "Command line arguments:", args$verbose)
write_log(odir, paste("     verbose:[", args$verbose, "]", sep=""), args$verbose)
write_log(odir, paste("     path:[", args$path, "]", sep=""), args$verbose)
write_log(odir, paste("     outputpath:[", args$outputpath, "]", sep=""), args$verbose)
write_log(odir, paste("     TBSSsubfolder:[", args$TBSSsubfolder, "]", sep=""), args$verbose)
write_log(odir, paste("     maskfile:[", args$maskfile, "]", sep=""), args$verbose)
write_log(odir, paste("     name:[", args$name, "]", sep=""), args$verbose)
write_log(odir, paste("     writecsv:[", args$writecsv, "]", sep=""), args$verbose)
write_log(odir, paste("     writepng:[", args$writepng, "]", sep=""), args$verbose)
write_log(odir, paste("     writemaskcoordinates:[", args$writemaskcoordinates, "]", sep=""), args$verbose)
write_log(odir, paste("     filetype:[", args$filetype, "]", sep=""), args$verbose)
write_log(odir, paste("     compression:[", args$compression, "]", sep=""), args$verbose)
write_log(odir, paste("     subjectsfile:[", args$subjectsfile, "]", sep=""), args$verbose)
write_log(odir, paste("     selectionfile:[", args$selectionfile, "]", sep=""), args$verbose)
write_log(odir, paste("     labelfile:[", args$labelfile, "]", sep=""), args$verbose)
write_log(odir, paste("     scalars:[", args$scalars, "]", sep=""), args$verbose)

# Resolve input arguments #
if ( args$verbose > 0) {
  write_print("Writing verbose output to standard error...\n", args$verbose-1) 
}
if ( nchar(args$path) > 0) {
  idir<-args$path
} else {
  idir<-getwd()
}
write_log(odir, paste("Using", idir, " as base directory for input"), args$verbose-1)
if ( nchar(args$name) > 0) {
  name<-args$name
} else {
  name<-paste('UNKNOWN', sub('(-: )','',Sys.time()), sep='')
}
write_print(paste("Using ", name, " as name for output", sep='')-1)

# Subject names
if ( nchar(args$subjectsfile) > 0) {
    subjectsfile<-args$subjectsfile
    write_print(paste("Using ", subjectsfile, " for subject names", sep=""), args$verbose-1)
    if(!file.exists(subjectsfile)) {
        write_log(paste("Subject name file ", subjectsfile, " not found", sep=''), args$verbose)
        stop("Subject name file ", subjectsfile, " not found", sep='')
    }
    # Read csv file
    if(str_detect(subjectsfile, ".csv")) {
        lines<-readLines(subjectsfile)
        subjectnames<-c()
        for(i in 2:length(lines)) {
            subjectname<-str_split(lines[i], ",")[[1]][1]
            subjectnames<-c(subjectnames, subjectname)
        }
    } else {
        subjectnames<-readLines(subjectsfile)
    }
} else {
    subjectnames_selected<-c()
    write_print(paste("Using running numbers as subject names."), args$verbose-1)
}

# Subject selection from all subjects
if ( nchar(args$selectionfile) > 0) {
    selectionfile<-args$selectionfile
    write_print(paste("Using ", selectionfile, " for subject selections", sep=""), args$verbose-1)
    if(!file.exists(selectionfile)) {
        write_log(paste("Subject name file ", selectionfile, " not found", sep=''), args$verbose)
        stop("Subject name file ", selectionfile, " not found", sep='')
    }
    subjectnames_selected<-readLines(selectionfile)
} else {
    subjectnames_selected<-c()
    write_print(paste("Using all subject names as selected."), args$verbose-1)
}

writecsv<-FALSE
if ( nchar(args$writecsv) > 0) {
    if (str_detect(args$writecsv, "Yes")) {
        write_print("Writing CSV output", args$verbose-1)
        writecsv<-TRUE
    }
}
writepng<-FALSE
if ( nchar(args$writepng) > 0) {
    if (str_detect(args$writepng, "Yes")) {
        write_print("Writing skeleton as PNG", args$verbose-1)
        writepng<-TRUE
    }
}
writemaskcoordinates<-FALSE
if ( nchar(args$writemaskcoordinates) > 0) {
    if (str_detect(args$writemaskcoordinates, "Yes")) {
        write_print("Writing mask coordinates in ASCII, this may take more time ", args$verbose-1)
        writemaskcoordinates<-TRUE
    }
}
if(args$filetype=="ASCII") {
  write_print("Using plain ASCII for Skifti data", args$verbose-1)
  datatype="volume-per-row-ASCII"
} else if(args$filetype=="rDs") {
  write_print("Using rDs for Skifti data", args$verbose-1)
  datatype="binary"
} else {
  write_log(paste("Unrecognized output filetype:", args$filetype, sep=''), args$verbose)
  stop("Unrecognized output filetype:", args$filetype, sep='')  
}
if(args$compression=="none") {
  write_print("No compression")
} else if(args$compression=="bz2") {
  write_print("bz2 compression")
} else if(args$compression=="zip") {
  write_print("zip compression")
} else {
  write_log(paste("Unrecognized output compression:", args$compression, sep=''), args$verbose)
  stop("Unrecognized output compression:", args$compression, sep='')  
}

# Scalars to be used
scalars<-str_split(args$scalars,',')[[1]]

# Verify that input files exist
if ( nchar(args$maskfile) > 0) {
  f2<-args$maskfile
} else {
  f2<-'/enigmaDTI/ENIGMA_DTI_FA_skeleton_mask.nii.gz'
}
write_log(odir, paste("Skeleton mask file:", f2, sep=''), args$verbose)
if(!file.exists(f2)) {
  write_log(odir, paste("Skeleton mask file ", f2, " not found", sep=''), args$verbose)
  stop("Skeleton mask file ", f2, " not found", sep='')
}
for(scalar in scalars) {
  f1<-paste(idir, args$TBSSsubfolder, 'stats', paste('all_', scalar, '_skeletonized.nii.gz', sep=''), sep='/')
  if(!file.exists(f1)) {
    write_log(odir, paste("Skeleton data file ", f1, " not found", sep=''), args$verbose)
    stop("Skeleton data file ", f1, " not found", sep='')
  }
}

if ( nchar(args$labelfile) > 0) {
    write_log(odir, paste("Label mask file:", args$labelfile, sep=''), args$labelfile)
} else {
    args$labelfile=NULL
}

# Handle one scalar at a time
for(scalar in scalars) {
  write_log(odir, paste("Starting to process ", scalar, sep=''), args$verbose)
  f1<-paste(idir, args$TBSSsubfolder, 'stats', paste('all_', scalar, '_skeletonized.nii.gz', sep=''), sep='/')
  # write and resolve subject names at the same time from the file, use list of selections, or select all if selections are not given
  selected_indexes<-c()
  selected_names<-c()
  for(i in 1:length(subjectnames)) {
    subjectname<-subjectnames[i]
    if(length(subjectnames_selected)>0) {
      if(subjectname %in% subjectnames_selected) {
        selected_indexes<-c(selected_indexes, i)
        selected_names<-c(selected_names, subjectname)
      }
    } else {
      selected_indexes<-c(selected_indexes, i)
      selected_names<-c(selected_names, subjectname)
    }
  }
  write_log(odir, paste("Subjectnames resolved, total ", length(selected_indexes) , " selected", sep=''), args$verbose)
  if(length(selected_indexes)==0){
    write_log(odir, paste("No subjects to process", sep=''), args$verbose)
  }
  
  # Create skifti object
  sk<-Nifti2Skifti(f1, f2, selected_indexes, write_coordinates=writemaskcoordinates, Nifti_labels=args$labelfile, verbose=(args$verbose > 0))
  if(length(rownames(sk$data)) != length(subjectnames)) {
    write_log(odir, paste("Rows in data ", length(rownames(sk$data)), " and number of subjectsnames ", length(subjectnames), " do not match", sep=''), args$verbose)
    stop("Rows in data ", length(rownames(sk$data)), " and number of subjectsnames ", length(subjectnames), " do not match", sep='')
  }
  # Set the selected subject names
  rownames(sk$data)<-subjectnames
  write_log(odir, paste("Skifti data created", sep=''), args$verbose)
  
  # Scalar avg
  output_avgfile<-paste(odir, '/', scalar, '_subject_list.csv', sep='')
  write_csv(selected_names, selected_indexes, output_avgfile)
  files_to_zip<-c(files_to_zip, output_avgfile)
  files_to_be_removed<-c(files_to_be_removed, output_avgfile)
  
  # Create list of used subjects
  if ( nchar(args$subjectsfile) > 0) {
    output_subjectsfile<-paste(odir, '/', 'casefile.txt', sep='')
    file.create(output_subjectsfile, showWarnings = FALSE)
    for(i in 1:length(subjectnames)) {
      write(subjectnames[i], file=output_subjectsfile, append = TRUE)
    }
    files_to_zip<-c(files_to_zip, output_subjectsfile)
    files_to_be_removed<-c(files_to_be_removed, output_subjectsfile)
  }
  write_log(odir, paste("Text files created"), args$verbose)

  # Set output datatype and write skifti data
  sk$datatype=datatype
  skfilename<-writeSkifti(sk, basename=paste(odir, '/', name, '_', scalar, '_Skiftidata', sep=''), overwrite = TRUE, compress=args$compression, verbose=(args$verbose > 0))
  files_to_zip<-c(files_to_zip, skfilename)
  files_to_be_removed<-c(files_to_be_removed, skfilename)
  write_log(odir, paste("Skifti data written to ", skfilename, sep=''), args$verbose)
  
  # Write CSV
  if(args$writecsv=="Yes") {
    filename<-writeCSV(sk, file=paste(odir, '/', name, ".csv", sep=''), overwrite = TRUE)
    files_to_zip<-c(files_to_zip, filename)
    files_to_be_removed<-c(files_to_be_removed, paste(odir, '/', name, ".csv", sep=''))
    write_log(odir, paste("CSV data written to ", name, ".csv", sep=''), args$verbose)
  }
  
  # Write graphical image of scalar
  if(args$writepng=="Yes") {
    write_log(odir, paste('writepng not supported in docker'), args$verbose)
    error(paste('writepng not supported in docker'))
    
    write_log(odir, paste("Creating graphical image", sep=''), args$verbose)
    img_hdr<-RNifti::niftiHeader(f1)
    total_volumes_to_read<-length(selected_indexes)
    for(i in length(total_volumes_to_read)) {
      if(is.null(selected_indexes)) {
        ii<-i
      } else {
        ii<-selected_indexes[i]
        if(ii < 1 | ii > img_hdr$dim[5]) {
          write_log(odir, aste('Selected volume index ', ii, ' out of bounds [1..', img_hdr$dim[5], ']', sep=' '), args$verbose)
          error(paste('Selected volume index ', ii, ' out of bounds [1..', img_hdr$dim[5], ']', sep=' '))
          return(NULL)
        }
      }     
      cat(paste("\r Reading ", scalar, " data ", i, "/", total_volumes_to_read, sep=""))      
      img<-RNifti::readNifti(f1, internal = FALSE, volumes = c(ii))
      if(i == 1) {
        mean_scalar<-img/total_volumes_to_read
      } else {
        mean_scalar<-mean_scalar+img/total_volumes_to_read
      }
    }
    write_print("")
    img_hdr<-niftiHeader(img)
    mask<-mean_scalar
    mask[mask>0]<-1
    mean_scalar[is.nan(mean_scalar)]<-0
    mean_scalar[mean_scalar>1]<-1
    if(str_detect(scalar, "FA")) {
      longname<-"Fractional Anisotropy"
    } else if(str_detect(scalar, "MD")) {
      longname<-"Mean Diffusivity"
    } else if(str_detect(scalar, "AD")) {
      longname<-"Axial Diffusivity"
    } else if(str_detect(scalar, "RD")) {
      longname<-"Radial Diffusivity"
    }
    save_skeleton(mask, mean_scalar, img_hdr, paste(odir, '/', name, "_mean_", scalar, ".png", sep=''), longname, 1000, keep_temp=TRUE)
    write_log(odir, paste("Skeleton rendering saved to ", odir, '/', name, "_mean_", scalar, ".png", sep=''), args$verbose)
    files_to_zip<-c(files_to_zip, paste(odir, '/', name, "_mean_", scalar, ".png", sep=''))
    files_to_be_removed<-c(files_to_be_removed, paste(odir, '/', name, "_mean_", scalar, ".png", sep=''))
  }
}

write_print("Files which are included in the resulting zip file", args$verbose)
for(file_i in 1:length(files_to_zip)) {
  write_print(paste("    ", files_to_zip[file_i], sep=""), args$verbose)
}

# Add log file with execution time to list of zipped files
end_time<-as.POSIXlt(Sys.time())
write_log(odir, paste("Execution time:", (end_time-start_time), ' seconds', sep=""), args$verbose)
files_to_zip<-c(files_to_zip, paste(odir, "log.txt", sep="/"))

zip(zipfile = paste(odir, '/LIFESPAN_', name, '.zip',sep=''), files = files_to_zip, flags="-j")
# Remove temporary files
for(file_i in 1:length(files_to_be_removed)) {
  invisible(file.remove(files_to_be_removed[file_i]))
}
write_print("Process ended successfully", 1)
