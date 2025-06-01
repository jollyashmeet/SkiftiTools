cp ../git/skiftiTools/R/Nifti2Skifti.R .
cp ../git/skiftiTools/R/writeSkifti.R .
cp ../git/skiftiTools/R/writeCSV.R .
cp ../git/skiftiTools/R/Skifti_view.R .
cp ../git/skiftiTools/R/Skifti2Nifti.R .
cp ../git/skiftiTools/R/Skifti_read.R .
cp ../R/Postprocess.R .
echo docker build . --no-cache -f Dockerfile.txt -t ashjoll/skiftitools:0.2.0 &> build_R.log
docker build . -f Dockerfile.txt -t ashjoll/skiftitools:0.2.0
