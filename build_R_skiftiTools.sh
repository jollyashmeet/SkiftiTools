curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/Nifti2Skifti.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/writeSkifti.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/writeCSV.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/Skifti_view.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/Skifti2Nifti.R
curl -O https://raw.githubusercontent.com/haanme/skiftiTools/main/R/Skifti_read.R
curl -O https://raw.githubusercontent.com/jollyashmeet/SkiftiTools/main/Postprocess.R
echo docker build . --no-cache -f Dockerfile.txt -t ashjoll/skiftitools:0.2.0 &> build_R.log
docker build . -f Dockerfile.txt -t ashjoll/skiftitools:0.2.0
