#pragma once

class _3DMultimodalRegistration
{
private:

	/* Static variables */
	static const unsigned int                          Dimension = 3;

public:
	typedef  float                              PixelType;
	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;
	typedef itk::ImageMaskSpatialObject<Dimension>::ImageType ImageMaskType;
	//typedef itk::VersorRigid3DTransform< double > TransformType;
	typedef itk::Euler3DTransform< double> TransformType;
	FixedImageType::Pointer fixedImage;
	MovingImageType::Pointer movingImage;

	//typedef itk::RegularStepGradientDescentOptimizerv4<double>   OptimizerType;
	typedef itk::PowellOptimizerv4<double>    OptimizerType;
	//typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double> OptimizerType;

	//typedef itk::MeanSquaresImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;
	typedef itk::MattesMutualInformationImageToImageMetricv4< FixedImageType, MovingImageType > MetricType;

	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType,
		TransformType >           RegistrationType;

	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();
	/*Private functions*/

private:
	/* Timers */
	itk::TimeProbesCollectorBase timer;

	/* flags */
	bool verbose = false;
	bool ok = false;
	bool customized_iso = false;
	bool moving_mask = false;
	bool fixed_mask = false;
	bool debug = false;
	bool resample = false;
	bool resolution = false;
	bool shrinking = false;

	/* operators */
	const double dtr = (atan(1.0) * 4.0) / 180.0;

	double cx = 0;
	double cy = -234;
	double cz = 70;

	/* Registration parameters */
	const unsigned int numberOfLevels = 1;
	double shrinkFactor = 1;
	FixedImageType::SpacingType ResampleSpacing;

	/* filenames */
	std::string Outputfilename = "CBCT_registered.mha";
	std::string OutputTransformfilename = "3DRigidRegistrationTransform.txt";
	//std::string Iterationfilename = "3DRigidRegistrationIt.txt";

	char* fixedImagefilename = NULL;
	char* movingImagefilename = NULL;
	char* movingMaskfilename = NULL;
	char* fixedMaskfilename = NULL;
	char* Orientation = "RAI";

	/* Types */

	TransformType::Pointer  initialTransform = TransformType::New();

	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
	typedef itk::ImageFileReader< ImageMaskType > MaskReaderType;
	typedef itk::ImageMaskSpatialObject < Dimension > MaskType;

	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	MaskReaderType::Pointer movingMaskReader = MaskReaderType::New();
	MaskReaderType::Pointer fixedMaskReader = MaskReaderType::New();

	MaskType::Pointer  spatialObjectMovingMask = MaskType::New();
	MaskType::Pointer  spatialObjectFixedMask = MaskType::New();

	typedef itk::CenteredTransformInitializer<
		TransformType,
		FixedImageType,
		MovingImageType >  TransformInitializerType;
	TransformInitializerType::Pointer initializer =
		TransformInitializerType::New();

	//typedef TransformType::VersorType  VersorType;
	//typedef VersorType::VectorType     VectorType;
	//VersorType     rotation;
	//VectorType     axis;

	typedef OptimizerType::ScalesType       OptimizerScalesType;
	const double translationScale = 1.0;
	const double rotationScale = 1.0 / dtr;

	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;

	typedef  unsigned char                                          OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension >                OutputImageType;
	//typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;
	typedef itk::ImageFileWriter< FixedImageType >                 WriterType;

	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;

	private:
		//bool Initialize(FixedImageType::Pointer &fixedImage, MovingImageType::Pointer &movingImage);
		bool Initialize();
		//bool Resample(FixedImageType::Pointer InputImage, double shrinkfactor, FixedImageType::Pointer &OutputImage);
		bool Resample(FixedImageType::Pointer InputImage, FixedImageType::SpacingType OutputSpacing, FixedImageType::Pointer &OutputImage);
	public:
		_3DMultimodalRegistration(int argc, char * argv[]);
		//RegistrationClass();
		virtual ~_3DMultimodalRegistration();
		bool StartRegistration();
		bool SetLevels();


};

