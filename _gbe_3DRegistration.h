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
#include "itkBSplineResampleImageFunction.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkCropImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkSpatialOrientationAdapter.h"

#include "itkImageMaskSpatialObject.h"

//#include "DicomImageReadPrintTags.h"
#include "gdcmReader.h"
#include "gdcmGlobal.h"
#include "gdcmDicts.h"
#include "gdcmDict.h"
#include "gdcmAttribute.h"
#include "gdcmStringFilter.h"
#include "gdcmItem.h"
#include "gdcmSequenceOfItems.h"
/***************************/

#define Reg33GetMacro(name, type)                                       \
  virtual type Get##name ()                                         \
    {                                                                 \
    return this->##name;                                          \
    }

#define Reg33SetMacro(name, type)                     \
  virtual void Set##name (const type _arg)          \
    {                                               \
    if ( this->##name != _arg )                   \
      {                                             \
      this->##name = _arg;                        \
      }                                             \
    }

#define Reg33SetCharMacro(name, type)                     \
  virtual void Set##name (type _arg)          \
    {                                               \
    if ( this->##name != _arg )                   \
      {                                             \
      this->##name = _arg;                        \
      }                                             \
    }

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

	typedef FixedImageType::SpacingType    SpacingType;
	typedef FixedImageType::PointType      OriginType;
	typedef FixedImageType::RegionType     RegionType;
	typedef FixedImageType::SizeType       SizeType;
	TransformType::TranslationType translation;

	/*Output Dimensions*/
	SpacingType FinalMovingSpacing;
	OriginType FinalMovingOrigin;
	SizeType FinalMovingSize;

	SpacingType OriginalMovingSpacing;
	OriginType OriginalMovingOrigin;
	SizeType OriginalMovingSize;

	SpacingType OriginalFixedSpacing;
	OriginType OriginalFixedOrigin;
	SizeType OriginalFixedSize;

	/*Enum types*/
	enum MetricSelection
	{
		MMI,
		MSE,
	};


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
	bool fixed_resample = false;
	bool moving_resample = false;
	bool resolution = false;
	bool shrinking = false;
	bool RTplan = false;
	bool like = false;
	bool autocrop = false; //enhances speed (is overriden by RTPlan)

	/* enumerators */
	MetricSelection RegistrationMetric = MMI;

	/* operators */
	const double dtr = (atan(1.0) * 4.0) / 180.0;

	double cx = 0;
	double cy = -234;
	double cz = 70;

	TransformType::InputPointType Isocenter;

	/* Registration parameters */
	const unsigned int numberOfLevels = 1;
	double shrinkFactor = 1;
	//FixedImageType::SpacingType ResampleSpacing;
	FixedImageType::SpacingType FixedImageResampleSpacing;
	FixedImageType::SpacingType MovingImageResampleSpacing;
	float DefaultPixelValue = -1024;

	/* filenames */
	std::string Outputfilename = "CBCT_registered.mha";
	std::string OutputTransformfilename = "3DRigidRegistrationTransform.txt";
	//std::string Iterationfilename = "3DRigidRegistrationIt.txt";

	//char* fixedImagefilename = NULL;
	//char* movingImagefilename = NULL;
	//char* movingMaskfilename = NULL;
	//char* fixedMaskfilename = NULL;
	//char* Orientation = "RAI";
	//char* RTplanFilename = NULL;

	std::string fixedImagefilename = std::string();
	std::string movingImagefilename = std::string();
	std::string movingMaskfilename = std::string();
	std::string fixedMaskfilename = std::string();
	std::string Orientation = "RAI";
	std::string RTplanFilename = std::string();

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

	ResampleFilterType::Pointer IsoAlignment = ResampleFilterType::New();

	//typedef itk::ExtractImageFilter<FixedImageType, FixedImageType> CropFixedFilterType;
	typedef itk::CropImageFilter<FixedImageType, FixedImageType> CropFixedFilterType;
	typedef itk::RegionOfInterestImageFilter<FixedImageType, FixedImageType> ROIFilterType;

	private:
		//bool Initialize(FixedImageType::Pointer &fixedImage, MovingImageType::Pointer &movingImage);
		//bool Resample(FixedImageType::Pointer InputImage, double shrinkfactor, FixedImageType::Pointer &OutputImage);
		bool IsocenterSearch(std::string argv, unsigned int &count, double Isocenter[][3]);
		bool Resample(FixedImageType::Pointer InputImage, FixedImageType::SpacingType OutputSpacing, FixedImageType::Pointer &OutputImage);
		bool Crop(FixedImageType::Pointer Image2Crop, FixedImageType::Pointer ReferenceImage, FixedImageType::Pointer &OutputImage);
		bool ROICrop(FixedImageType::Pointer Image2Crop, FixedImageType::Pointer ReferenceImage, FixedImageType::Pointer & OutputImage);

	public:
		_3DRegistration();
		_3DRegistration(int argc, char * argv[]);
		virtual ~_3DRegistration();
		
		template<typename TMetricType>
		bool Initialize();

		bool StartRegistration();
		bool SetLevels();
		bool ReadIsocenter();

		/*Set/Get Methods*/
		Reg33GetMacro(verbose, bool);
		Reg33SetMacro(verbose, bool);
		Reg33GetMacro(debug, bool);
		Reg33SetMacro(debug, bool);
		Reg33GetMacro(moving_mask, bool);
		Reg33SetMacro(moving_mask, bool);
		Reg33GetMacro(fixed_mask, bool);
		Reg33SetMacro(fixed_mask, bool);
		Reg33GetMacro(RTplan, bool);
		Reg33SetMacro(RTplan, bool);
		Reg33GetMacro(moving_resample, bool);
		Reg33SetMacro(moving_resample, bool);
		Reg33GetMacro(fixed_resample, bool);
		Reg33SetMacro(fixed_resample, bool);
		//Reg33GetMacro(resolution, bool);
		//Reg33SetMacro(resolution, bool);
		Reg33GetMacro(autocrop, bool);
		Reg33SetMacro(autocrop, bool);

		Reg33GetMacro(RegistrationMetric, MetricSelection);
		Reg33SetMacro(RegistrationMetric, MetricSelection);

		Reg33GetMacro(fixedImagefilename, std::string);
		Reg33SetMacro(fixedImagefilename, std::string);
		Reg33GetMacro(movingImagefilename, std::string);
		Reg33SetMacro(movingImagefilename, std::string);
		Reg33GetMacro(fixedMaskfilename, std::string);
		Reg33SetMacro(fixedMaskfilename, std::string);
		Reg33GetMacro(movingMaskfilename, std::string);
		Reg33SetMacro(movingMaskfilename, std::string);
		Reg33GetMacro(RTplanFilename, std::string);
		Reg33SetMacro(RTplanFilename, std::string);		
		Reg33GetMacro(Outputfilename, std::string);
		Reg33SetMacro(Outputfilename, std::string);

		Reg33GetMacro(DefaultPixelValue, float);
		Reg33SetMacro(DefaultPixelValue, float);

		/*char* fixedImagefilename = NULL;
		char* movingImagefilename = NULL;
		char* movingMaskfilename = NULL;
		char* fixedMaskfilename = NULL;
		char* Orientation = "RAI";
		char* RTplanFilename = NULL;
*/
		itkBooleanMacro(verbose);
		itkBooleanMacro(debug);
		itkBooleanMacro(moving_mask);
		itkBooleanMacro(fixed_mask);
		itkBooleanMacro(RTplan);
		itkBooleanMacro(autocrop);
		itkBooleanMacro(moving_resample);
		itkBooleanMacro(fixed_resample);
		//itkBooleanMacro(resolution);

};

