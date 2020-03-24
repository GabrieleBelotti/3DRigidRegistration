#pragma once

#include "itkTimeProbesCollectorBase.h"
//#include "itkFileOutputWindow.h"

#include "itkImageRegistrationMethodv4.h"
#include "itkImageToImageMetricv4.h"//metric base we include for the pointer
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
//#include "itkCorrelationImageToImageMetricv4.h"

#include "itkVersorRigid3DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkCenteredTransformInitializer.h"

#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkPowellOptimizerv4.h"
#include "itkLBFGSOptimizerv4.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"

#include "itkImageMaskSpatialObject.h"

#include "DicomImageReadPrintTags.cxx"


class _3DRegistration
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

	typedef itk::ObjectToObjectOptimizerBaseTemplate<double> OptimizerBaseType;
	//typedef itk::LBFGSOptimizerv4 LBFGSOptimizerType;
	typedef itk::LBFGSOptimizerv4 OptimizerType;
	//typedef itk::RegularStepGradientDescentOptimizerv4<double>   OptimizerType;
	//typedef itk::PowellOptimizerv4<double>    OptimizerType;
	typedef itk::PowellOptimizerv4<double>    PowellOptimizerType;
	//typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double> OptimizerType;

	typedef itk::ImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;
	typedef itk::MeanSquaresImageToImageMetricv4< FixedImageType, MovingImageType >   MeanSquaresMetricType;
	typedef itk::MattesMutualInformationImageToImageMetricv4< FixedImageType, MovingImageType > MIMetricType;

	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType,
		TransformType >           RegistrationType;

	MetricType::Pointer         metric;
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
	bool RTplan = false;

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
	char* RTplanFilename = NULL;

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
		//bool Resample(FixedImageType::Pointer InputImage, double shrinkfactor, FixedImageType::Pointer &OutputImage);
		bool Resample(FixedImageType::Pointer InputImage, FixedImageType::SpacingType OutputSpacing, FixedImageType::Pointer &OutputImage);
	public:
		_3DRegistration(int argc, char * argv[]);
		//RegistrationClass();
		virtual ~_3DRegistration();
		
		template<typename TMetricType>
		bool Initialize();

		bool StartRegistration();
		bool SetLevels();


};

