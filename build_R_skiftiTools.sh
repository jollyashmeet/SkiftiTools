# Download all necessary R scripts from haanme's repo
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/Nifti2Skifti.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/writeSkifti.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/writeCSV.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/Skifti_view.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/Skifti2Nifti.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/Skifti_read.R

## Download all necessary R scripts from jollyashmeet's repo
curl -O https://raw.githubusercontent.com/jollyashmeet/SkiftiTools/main/Postprocess.R

# Make enigmaDTI folder and download the skeleton mask
mkdir -p enigmaDTI
curl -o enigmaDTI/ENIGMA_DTI_FA_skeleton_mask.nii.gz https://raw.githubusercontent.com/jollyashmeet/SkiftiTools/main/enigmaDTI/ENIGMA_DTI_FA_skeleton_mask.nii.gz

#Build Docker image
echo docker build . --no-cache -f Dockerfile.txt -t ashjoll/skiftitools:0.2.0 &> build_R.log
docker build . --no-cache -f Dockerfile.txt -t ashjoll/skiftitools:0.2.0
