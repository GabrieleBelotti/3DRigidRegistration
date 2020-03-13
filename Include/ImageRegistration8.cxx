/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkTimeProbesCollectorBase.h"
//#include "itkFileOutputWindow.h"

#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
//#include "itkCorrelationImageToImageMetricv4.h"

#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"

#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkPowellOptimizerv4.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"

#include "itkImageMaskSpatialObject.h"

/* GLOBAL VARIABLES */
std::string OutputTransformfilename = "3DRigidRegistrationTransform.txt";
std::string Iterationfilename = "3DRigidRegistrationIt.txt";

/*FUNCTION PROTOTYPES*/
void exe_usage();

/*OBSERVER CLASSES*/
#include "itkCommand.h"
template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{

public:
	using Self = RegistrationInterfaceCommand;
	using Superclass = itk::Command;
	using Pointer = itk::SmartPointer<Self>;
	itkNewMacro(Self);

protected:
	RegistrationInterfaceCommand() = default;
public:
	using RegistrationType = TRegistration;
	using RegistrationPointer = RegistrationType *;
	//using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	typedef itk::PowellOptimizerv4<double> OptimizerType;
	using OptimizerPointer = OptimizerType *;

	void
		Execute(itk::Object* object, const itk::EventObject& event) override
	{
		// Software Guide : EndCodeSnippet

		// Software Guide : BeginLatex
		//
		// First we verify that the event invoked is of the right type,
		// \code{itk::MultiResolutionIterationEvent()}.
		// If not, we return without any further action.
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		if (!(itk::MultiResolutionIterationEvent().CheckEvent(&event)))
		{
			return;
		}
		// Software Guide : EndCodeSnippet

		// Software Guide : BeginLatex
		//
		// We then convert the input object pointer to a RegistrationPointer.
		// Note that no error checking is done here to verify the
		// \code{dynamic\_cast} was successful since we know the actual object
		// is a registration method. Then we ask for the optimizer object
		// from the registration method.
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		auto registration = static_cast<RegistrationPointer>(object);
		auto optimizer =
			static_cast<OptimizerPointer>(registration->GetModifiableOptimizer());
		// Software Guide : EndCodeSnippet

		unsigned int currentLevel = registration->GetCurrentLevel();
		typename RegistrationType::ShrinkFactorsPerDimensionContainerType shrinkFactors =
			registration->GetShrinkFactorsPerDimension(currentLevel);
		typename RegistrationType::SmoothingSigmasArrayType smoothingSigmas =
			registration->GetSmoothingSigmasPerLevel();

		std::cout << "-------------------------------------" << std::endl;
		std::cout << " Current level = " << currentLevel << std::endl;
		std::cout << "    shrink factor = " << shrinkFactors << std::endl;
		std::cout << "    smoothing sigma = ";
		std::cout << smoothingSigmas[currentLevel] << std::endl;
		std::cout << std::endl;

		// Software Guide : BeginLatex
		//
		// If this is the first resolution level we set the learning rate
		// (representing the first step size) and the minimum step length (representing
		// the convergence criterion) to large values.  At each subsequent resolution
		// level, we will reduce the minimum step length by a factor of 5 in order to
		// allow the optimizer to focus on progressively smaller regions. The learning
		// rate is set up to the current step length. In this way, when the
		// optimizer is reinitialized at the beginning of the registration process for
		// the next level, the step length will simply start with the last value used
		// for the previous level. This will guarantee the continuity of the path
		// taken by the optimizer through the parameter space.
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		//if (registration->GetCurrentLevel() == 0)
		//{
		//	optimizer->SetLearningRate(16.00);
		//	optimizer->SetMinimumStepLength(2.5);
		//}
		//else
		//{
		//	optimizer->SetLearningRate(optimizer->GetCurrentStepLength());
		//	optimizer->SetMinimumStepLength(optimizer->GetMinimumStepLength() * 0.2);
		//}
		// Software Guide : EndCodeSnippet
	}

	// Software Guide : BeginLatex
	//
	// Another version of the \code{Execute()} method accepting a \code{const}
	// input object is also required since this method is defined as pure virtual
	// in the base class.  This version simply returns without taking any action.
	//
	// Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet
	void
		Execute(const itk::Object*, const itk::EventObject&) override
	{
		return;
	}
};

class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:
  //typedef itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
  //typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double> OptimizerType;
  typedef itk::PowellOptimizerv4<double>  OptimizerType;
  typedef   const OptimizerType *                             OptimizerPointer;
  std::ofstream windowObs;

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
	  Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
    OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
	if (!windowObs.is_open())
		windowObs.open(Iterationfilename);
	windowObs << optimizer->GetCurrentIteration() << "   ";
	windowObs << optimizer->GetValue() << "   ";
	windowObs << optimizer->GetCurrentPosition() << "\n\n";
	
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};

int MultimodalRegistration(int argc, char *argv[])
{
	std::cout << "Not yet implemented\n";
	return EXIT_SUCCESS;
}

void exe_usage()
{
	std::cerr << "\n";
	std::cerr << "Usage: ImageRegistration3DClass <unsigned int registration_routine> <options> CTname CBCTname -o outputfilename(default:CBCT_ax_aligned_iso.mha)\n";
	std::cerr << "   where <options> is one or more of the following:\n\n";
	//std::cerr << "       <-A>                     Use GPU hardware (to be implemented) [default: cpu]\n";
	std::cerr << "       <-h>							Display (this) usage information\n";
	std::cerr << "       <-v>							Verbose output [default: no]\n";
	std::cerr << "       <-mm filename>					Moving Image Mask [default: no]\n";
	std::cerr << "       <-mf filename>					Fixed Image Mask [default: no] - under construction\n";
	std::cerr << "       <-par filename>				Output parameters filename [default: yes]\n";
	//std::cerr << "       <-p>							Permute image convention order (coronal to axial first) [default: no]\n";
	//std::cerr << "       <-iso double double double>    Align CBCT to isocenter [default: yes]\n";
	//std::cerr << "       <-s>							Scale image intensity [default: no] - calibration work in progress \n";
	//std::cerr << "       <-dbg>                   Debugging output [default: no]\n";

	std::cerr << "       <-o filename>					Output image filename\n\n";
	std::cerr << "										by  Gabriele Belotti\n";
	std::cerr << "										gabriele.belotti@mail.polimi.it\n";
	std::cerr << "										(Politecnico di Milano)\n\n";
	exit(EXIT_FAILURE);
}

/*BACKUP UNTEMPLATED CLASS*/
class RegistrationClass
{
private:

	/* Static variables */
	static const unsigned int                          Dimension = 3;

public:
	RegistrationClass(int argc, char * argv[]);
	//RegistrationClass();
	~RegistrationClass();
	bool StartRegistration();
	bool SetLevels();

	typedef  float                              PixelType;
	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;
	typedef itk::ImageMaskSpatialObject<Dimension>::ImageType ImageMaskType;
	typedef itk::VersorRigid3DTransform< double > TransformType;

	//typedef itk::RegularStepGradientDescentOptimizerv4<double>   OptimizerType;
	typedef itk::PowellOptimizerv4<double>    OptimizerType;
	//typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double> OptimizerType;

	typedef itk::MeanSquaresImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;
	//typedef itk::MattesMutualInformationImageToImageMetricv4< FixedImageType, MovingImageType > MetricType;

	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType,
		TransformType >           RegistrationType;

	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();
	/*Private functions*/
private:
	bool Initialize();

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

	/* operators */
	const double dtr = (atan(1.0) * 4.0) / 180.0;

	double cx = 0;
	double cy = -234;
	double cz = 70;

	/* Registration parameters */
	const unsigned int numberOfLevels = 1;

	/* filenames */
	std::string Outputfilename = "CBCT_registered.mha";
	std::string OutputTransformfilename = "3DRigidRegistrationTransform.txt";
	std::string Iterationfilename = "3DRigidRegistrationIt.txt";

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

	typedef TransformType::VersorType  VersorType;
	typedef VersorType::VectorType     VectorType;
	VersorType     rotation;
	VectorType     axis;
	typedef OptimizerType::ScalesType       OptimizerScalesType;
	const double translationScale = 1.0;
	const double rotationScale = 1.0 / dtr;

	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;

};

RegistrationClass::RegistrationClass(int argc, char *argv[])
{
	if (argc <= 3)
	{
		std::cout << "Not enough input arguments\n";
		exe_usage();
	}
	while (argc > 1)
	{
		this->ok = false;

		if ((this->ok == false) && (strcmp(argv[1], "-iso") == 0))
		{
			argc--; argv++;
			this->cx = atof(argv[1]);
			argc--; argv++;
			this->cy = atof(argv[1]);
			argc--; argv++;
			this->cz = atof(argv[1]);
			argc--; argv++;
			this->customized_iso = true;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-h") == 0))
		{
			argc--; argv++;
			this->ok = true;
			exe_usage();
		}

		if ((this->ok == false) && (strcmp(argv[1], "-v") == 0))
		{
			argc--; argv++;
			this->ok = true;
			this->verbose = true;
			continue;
			std::cout << "Verbose execution selected\n";
		}

		if ((ok == false) && (strcmp(argv[1], "-dbg") == 0))
		{
			argc--; argv++;
			ok = true;
			debug = true;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-o") == 0))
		{
			argc--; argv++;
			this->ok = true;
			this->Outputfilename = argv[1];
			argc--; argv++;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-par") == 0))
		{
			argc--; argv++;
			this->ok = true;
			this->OutputTransformfilename = argv[1];
			argc--; argv++;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-mm") == 0))
		{
			argc--; argv++;
			this->ok = true;
			this->movingMaskfilename = argv[1];
			this->moving_mask = true;
			argc--; argv++;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-mf") == 0))
		{
			argc--; argv++;
			this->ok = true;
			this->fixedMaskfilename = argv[1];
			this->fixed_mask = true;
			argc--; argv++;
			continue;
		}

		if (ok == false)
		{
			if (!strcmp(argv[1], "1"))
			{
				std::cout << "Multimodality selected\n";
				argc--;
				argv++;
			}
			else if (!strcmp(argv[1], "2"))
			{
				std::cout << "Monomodality selected\n";
				argc--;
				argv++;
			}
			else if (this->fixedImagefilename == NULL)
			{
				this->fixedImagefilename = argv[1];
				argc--;
				argv++;
			}

			else if (this->movingImagefilename == NULL)
			{
				this->movingImagefilename = argv[1];
				argc--;
				argv++;
			}
			else
			{
				std::cerr << "ERROR: Cannot parse argument " << argv[1] << std::endl;
				exe_usage();
			}
		}
	}
}

RegistrationClass::~RegistrationClass()
{
}

bool RegistrationClass::StartRegistration()
{
	this->Initialize();


	try
	{
		if (debug)
		{
			metric->Print(std::cout);
		}
		std::cout << "Start Registration \n";
		timer.Start("Registration");
		registration->Update();
		std::cout << "Optimizer stop condition: "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
		timer.Stop("Registration");
	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	std::ofstream parametersFile;
	//parametersFile.open(OutputTransformfilename, std::ofstream::out | std::ofstream::app);
	parametersFile.open(OutputTransformfilename);

	const TransformType::ParametersType finalParameters =
		registration->GetOutput()->Get()->GetParameters();

	const double versorX = finalParameters[0];
	const double versorY = finalParameters[1];
	const double versorZ = finalParameters[2];
	const double finalTranslationX = finalParameters[3];
	const double finalTranslationY = finalParameters[4];
	const double finalTranslationZ = finalParameters[5];
	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();

	// Print out results
	//
	std::cout << std::endl << std::endl;
	std::cout << "Result = " << std::endl;
	std::cout << " versor X      = " << versorX << std::endl;
	std::cout << " versor Y      = " << versorY << std::endl;
	std::cout << " versor Z      = " << versorZ << std::endl;
	std::cout << " Translation X = " << finalTranslationX << std::endl;
	std::cout << " Translation Y = " << finalTranslationY << std::endl;
	std::cout << " Translation Z = " << finalTranslationZ << std::endl;
	std::cout << " Iterations    = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue << std::endl;

	TransformType::Pointer finalTransform = TransformType::New();

	finalTransform->SetFixedParameters(registration->GetOutput()->Get()->GetFixedParameters());
	finalTransform->SetParameters(finalParameters);

	// Software Guide : BeginCodeSnippet
	TransformType::MatrixType matrix = finalTransform->GetMatrix();
	TransformType::OffsetType offset = finalTransform->GetOffset();
	std::cout << "Matrix = " << std::endl << matrix << std::endl;
	std::cout << "Offset = " << std::endl << offset << std::endl;

	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();

	resampler->SetTransform(finalTransform);
	resampler->SetInput(movingImageReader->GetOutput());

	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

	resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resampler->SetOutputOrigin(fixedImage->GetOrigin());
	resampler->SetOutputSpacing(fixedImage->GetSpacing());
	resampler->SetOutputDirection(fixedImage->GetDirection());
	//resampler->SetDefaultPixelValue( 100 );

	typedef  unsigned char                                          OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension >                OutputImageType;
	//typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;
	typedef itk::ImageFileWriter< FixedImageType >                 WriterType;

	WriterType::Pointer      writer = WriterType::New();
	//CastFilterType::Pointer  caster =  CastFilterType::New();

	writer->SetFileName(Outputfilename);

	//caster->SetInput( resampler->GetOutput() );
	writer->SetInput(resampler->GetOutput());
	try
	{
		std::cout << "Writing transformed CBCT\n";
		timer.Start("WriteCBCT");
		writer->Update();
		timer.Stop("WriteCBCT");
	}
	catch (itk::ExceptionObject& error)
	{
		return EXIT_FAILURE;
	}

	//std::ofstream parametersFile;
	//parametersFile.open(OutputTransformfilename);
	parametersFile << std::endl << std::endl;
	parametersFile << "Result = " << std::endl;
	parametersFile << " versor X      = " << versorX << std::endl;
	parametersFile << " versor Y      = " << versorY << std::endl;
	parametersFile << " versor Z      = " << versorZ << std::endl;
	parametersFile << " Translation X = " << finalTranslationX << std::endl;
	parametersFile << " Translation Y = " << finalTranslationY << std::endl;
	parametersFile << " Translation Z = " << finalTranslationZ << std::endl;
	parametersFile << " Iterations    = " << numberOfIterations << std::endl;
	parametersFile << " Metric value  = " << bestValue << std::endl;
	parametersFile << "\nParameters :\n";
	parametersFile << finalParameters << std::endl;
	//parametersFile << "Timing :\n";
	//parametersFile << timer.
	//parametersFile << registration->GetOutput()->Get()->get

	parametersFile.close();
	//timer.ExpandedReport();
	timer.Report();
	return EXIT_SUCCESS;

}

bool RegistrationClass::Initialize()
{
	fixedImageReader->SetFileName(fixedImagefilename);
	movingImageReader->SetFileName(movingImagefilename);
	try
	{
		std::cout << "Reading CT" << std::endl;
		timer.Start("ReadCT");
		fixedImageReader->Update();
		std::cout << "Done\n";
		timer.Stop("ReadCT");
	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "Exception caught" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	try
	{
		std::cout << "Reading CBCT" << std::endl;
		timer.Start("ReadCBCT");
		movingImageReader->Update();
		std::cout << "Done\n";
		timer.Stop("ReadCBCT");

	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "Exception caught" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	if (fixed_mask)
	{
		fixedMaskReader->SetFileName(movingMaskfilename);
		try
		{
			std::cout << "Reading Fixed Image Mask" << std::endl;
			timer.Start("ReadCBCTMask");
			fixedMaskReader->Update();
			std::cout << "Done\n";
			timer.Stop("ReadCBCTMask");

		}
		catch (itk::ExceptionObject& err)
		{
			std::cerr << "Exception caught" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
		spatialObjectFixedMask->SetImage(fixedMaskReader->GetOutput());
		try
		{
			//std::cout << "Reading Moving Image Mask" << std::endl;
			spatialObjectFixedMask->Update();
			//std::cout << "Done\n";

		}
		catch (itk::ExceptionObject& err)
		{
			std::cerr << "Exception caught" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
		if (debug)
			spatialObjectFixedMask->Print(std::cout);
		metric->SetFixedImageMask(spatialObjectFixedMask);
	}

	if (moving_mask)
	{
		movingMaskReader->SetFileName(movingMaskfilename);
		try
		{
			std::cout << "Reading Moving Image Mask" << std::endl;
			timer.Start("ReadCBCTMask");
			movingMaskReader->Update();
			std::cout << "Done\n";
			timer.Stop("ReadCBCTMask");

		}
		catch (itk::ExceptionObject& err)
		{
			std::cerr << "Exception caught" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
		spatialObjectMovingMask->SetImage(movingMaskReader->GetOutput());
		try
		{
			//std::cout << "Reading Moving Image Mask" << std::endl;
			spatialObjectMovingMask->Update();
			//std::cout << "Done\n";

		}
		catch (itk::ExceptionObject& err)
		{
			std::cerr << "Exception caught" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
		if (debug)
			spatialObjectMovingMask->Print(std::cout);
		metric->SetMovingImageMask(spatialObjectMovingMask);
	}


	if (verbose)
		std::cout << "Set Images to Registration method\n";
	registration->SetFixedImage(fixedImageReader->GetOutput());
	registration->SetMovingImage(movingImageReader->GetOutput());
	if (verbose)
		std::cout << "Done\n";

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);

	if (verbose)
		std::cout << "Create Initial transform\n";
	initializer->SetTransform(initialTransform);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());
	initializer->GeometryOn();
	if (verbose)
		std::cout << "Done\n";
	try
	{
		if (verbose)
			std::cout << "Update Initializer\n";
		initializer->InitializeTransform();
		if (verbose)
			std::cout << "Done\n";
	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "Exception caught" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	if (verbose)
		std::cout << "Set Initial rotation \n";
	axis[0] = 0.0;
	axis[1] = 0.0;
	axis[2] = 1.0;
	const double angle = 0;
	rotation.Set(axis, angle);
	initialTransform->SetRotation(rotation);
	if (verbose)
		std::cout << "Done \n";

	registration->SetInitialTransform(initialTransform);
	OptimizerScalesType optimizerScales(initialTransform->GetNumberOfParameters());

	optimizerScales[0] = rotationScale;
	optimizerScales[1] = rotationScale;
	optimizerScales[2] = rotationScale;
	optimizerScales[3] = translationScale;
	optimizerScales[4] = translationScale;
	optimizerScales[5] = translationScale;
	optimizer->SetScales(optimizerScales);
	optimizer->SetNumberOfIterations(200);

	/* POWELL */

	optimizer->SetMaximumIteration(10);
	optimizer->SetMaximumLineIteration(4); // for Powell's method
	optimizer->SetStepLength(4.0);
	optimizer->SetStepTolerance(0.02);
	optimizer->SetValueTolerance(0.001);
	optimizer->DebugOn();

	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);

	this->SetLevels();
	return EXIT_SUCCESS;

}

bool RegistrationClass::SetLevels()
{
	this->shrinkFactorsPerLevel.SetSize(this->numberOfLevels);
	this->smoothingSigmasPerLevel.SetSize(this->numberOfLevels);

	for (unsigned int i = 0; i < this->numberOfLevels; i++)
	{
		this->shrinkFactorsPerLevel[i] = 1;
		this->smoothingSigmasPerLevel[i] = 0;
	}
	registration->SetNumberOfLevels(this->numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
	return EXIT_SUCCESS;
}

//<class T>
//class RegistrationClass
//{
//public:
//	RegistrationClass(int argc, char * argv[]);
//	//RegistrationClass();
//	~RegistrationClass();
//	bool StartRegistration();
//	bool SetLevels();
//	/*Private functions*/
//private:
//	bool Initialize();
//
//private:
//	/* Timers */
//	itk::TimeProbesCollectorBase timer;
//
//	/* flags */
//	bool verbose = false;
//	bool ok = false;
//	bool customized_iso = false;
//	bool moving_mask = false;
//	bool fixed_mask = false;
//	bool debug = false;
//
//	/* operators */
//	const double dtr = (atan(1.0) * 4.0) / 180.0;
//
//	double cx = 0;
//	double cy = -234;
//	double cz = 70;
//
//	/* Registration parameters */
//	const unsigned int numberOfLevels = 1;
//
//	/* Static vaariables */
//	static const unsigned int                          Dimension = 3;
//
//	/* filenames */
//	std::string Outputfilename = "CBCT_registered.mha";
//	std::string OutputTransformfilename = "3DRigidRegistrationTransform.txt";
//	std::string Iterationfilename = "3DRigidRegistrationIt.txt";
//
//	char* fixedImagefilename = NULL;
//	char* movingImagefilename = NULL;
//	char* movingMaskfilename = NULL;
//	char* fixedMaskfilename = NULL;
//	char* Orientation = "RAI";
//
//	/* Types */
//
//	typedef  float                              PixelType;
//	typedef itk::Image< PixelType, Dimension >  FixedImageType;
//	typedef itk::Image< PixelType, Dimension >  MovingImageType;
//	typedef itk::ImageMaskSpatialObject<Dimension>::ImageType ImageMaskType;
//	typedef itk::VersorRigid3DTransform< double > TransformType;
//
//	//typedef itk::RegularStepGradientDescentOptimizerv4<double>   OptimizerType;
//	typedef itk::PowellOptimizerv4<double>    OptimizerType;
//	//typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double> OptimizerType;
//
//	typedef itk::MeanSquaresImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;
//	//typedef itk::MeanSquaresImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;
//
//	//typedef itk::MattesMutualInformationImageToImageMetricv4< FixedImageType, MovingImageType > MetricType;
//
//	typedef itk::ImageRegistrationMethodv4<
//		FixedImageType,
//		MovingImageType,
//		TransformType >           RegistrationType;
//
//	MetricType::Pointer         metric = MetricType::New();
//	OptimizerType::Pointer      optimizer = OptimizerType::New();
//	RegistrationType::Pointer   registration = RegistrationType::New();
//
//	TransformType::Pointer  initialTransform = TransformType::New();
//
//	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
//	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
//	typedef itk::ImageFileReader< ImageMaskType > MaskReaderType;
//	typedef itk::ImageMaskSpatialObject < Dimension > MaskType;
//
//	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
//	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
//	MaskReaderType::Pointer movingMaskReader = MaskReaderType::New();
//	MaskReaderType::Pointer fixedMaskReader = MaskReaderType::New();
//
//	MaskType::Pointer  spatialObjectMovingMask = MaskType::New();
//	MaskType::Pointer  spatialObjectFixedMask = MaskType::New();
//
//	typedef itk::CenteredTransformInitializer<
//		TransformType,
//		FixedImageType,
//		MovingImageType >  TransformInitializerType;
//	TransformInitializerType::Pointer initializer =
//		TransformInitializerType::New();
//
//	typedef TransformType::VersorType  VersorType;
//	typedef VersorType::VectorType     VectorType;
//	VersorType     rotation;
//	VectorType     axis;
//	typedef OptimizerType::ScalesType       OptimizerScalesType;
//	const double translationScale = 1.0;
//	const double rotationScale = 1.0 / dtr;
//
//	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
//	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
//
//};
//
//RegistrationClass<class T>::RegistrationClass(int argc, char *argv[])
//{
//	if (argc <= 3)
//	{
//		std::cout << "Not enough input arguments\n";
//		exe_usage();
//	}
//	while (argc > 1)
//	{
//		this->ok = false;
//
//		if ((this->ok == false) && (strcmp(argv[1], "-iso") == 0))
//		{
//			argc--; argv++;
//			this->cx = atof(argv[1]);
//			argc--; argv++;
//			this->cy = atof(argv[1]);
//			argc--; argv++;
//			this->cz = atof(argv[1]);
//			argc--; argv++;
//			this->customized_iso = true;
//			continue;
//		}
//
//		if ((this->ok == false) && (strcmp(argv[1], "-h") == 0))
//		{
//			argc--; argv++;
//			this->ok = true;
//			exe_usage();
//		}
//
//		if ((this->ok == false) && (strcmp(argv[1], "-v") == 0))
//		{
//			argc--; argv++;
//			this->ok = true;
//			this->verbose = true;
//			continue;
//			std::cout << "Verbose execution selected\n";
//		}
//
//		if ((ok == false) && (strcmp(argv[1], "-dbg") == 0))
//		{
//			argc--; argv++;
//			ok = true;
//			debug = true;
//			continue;
//		}
//
//		if ((this->ok == false) && (strcmp(argv[1], "-o") == 0))
//		{
//			argc--; argv++;
//			this->ok = true;
//			this->Outputfilename = argv[1];
//			argc--; argv++;
//			continue;
//		}
//
//		if ((this->ok == false) && (strcmp(argv[1], "-par") == 0))
//		{
//			argc--; argv++;
//			this->ok = true;
//			this->OutputTransformfilename = argv[1];
//			argc--; argv++;
//			continue;
//		}
//
//		if ((this->ok == false) && (strcmp(argv[1], "-mm") == 0))
//		{
//			argc--; argv++;
//			this->ok = true;
//			this->movingMaskfilename = argv[1];
//			this->moving_mask = true;
//			argc--; argv++;
//			continue;
//		}
//
//		if ((this->ok == false) && (strcmp(argv[1], "-mf") == 0))
//		{
//			argc--; argv++;
//			this->ok = true;
//			this->fixedMaskfilename = argv[1];
//			this->fixed_mask = true;
//			argc--; argv++;
//			continue;
//		}
//
//		if (ok == false)
//		{
//			if (!strcmp(argv[1], "1"))
//			{
//				std::cout << "Multimodality selected\n";
//				argc--;
//				argv++;
//			}
//			else if (!strcmp(argv[1], "2"))
//			{
//				std::cout << "Monomodality selected\n";
//				argc--;
//				argv++;
//			}
//			else if (this->fixedImagefilename == NULL)
//			{
//				this->fixedImagefilename = argv[1];
//				argc--;
//				argv++;
//			}
//
//			else if (this->movingImagefilename == NULL)
//			{
//				this->movingImagefilename = argv[1];
//				argc--;
//				argv++;
//			}
//			else
//			{
//				std::cerr << "ERROR: Cannot parse argument " << argv[1] << std::endl;
//				exe_usage();
//			}
//		}
//	}
//}
//
//RegistrationClass<class T>::~RegistrationClass()
//{
//}
//
//bool RegistrationClass<class T>::StartRegistration()
//{
//	this->Initialize();
//
//
//	try
//	{
//		if (debug)
//		{
//			metric->Print(std::cout);
//		}
//		std::cout << "Start Registration \n";
//		timer.Start("Registration");
//		registration->Update();
//		std::cout << "Optimizer stop condition: "
//			<< registration->GetOptimizer()->GetStopConditionDescription()
//			<< std::endl;
//		timer.Stop("Registration");
//	}
//	catch (itk::ExceptionObject& err)
//	{
//		std::cerr << "ExceptionObject caught !" << std::endl;
//		std::cerr << err << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	std::ofstream parametersFile;
//	//parametersFile.open(OutputTransformfilename, std::ofstream::out | std::ofstream::app);
//	parametersFile.open(OutputTransformfilename);
//
//	const TransformType::ParametersType finalParameters =
//		registration->GetOutput()->Get()->GetParameters();
//
//	const double versorX = finalParameters[0];
//	const double versorY = finalParameters[1];
//	const double versorZ = finalParameters[2];
//	const double finalTranslationX = finalParameters[3];
//	const double finalTranslationY = finalParameters[4];
//	const double finalTranslationZ = finalParameters[5];
//	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
//	const double bestValue = optimizer->GetValue();
//
//	// Print out results
//	//
//	std::cout << std::endl << std::endl;
//	std::cout << "Result = " << std::endl;
//	std::cout << " versor X      = " << versorX << std::endl;
//	std::cout << " versor Y      = " << versorY << std::endl;
//	std::cout << " versor Z      = " << versorZ << std::endl;
//	std::cout << " Translation X = " << finalTranslationX << std::endl;
//	std::cout << " Translation Y = " << finalTranslationY << std::endl;
//	std::cout << " Translation Z = " << finalTranslationZ << std::endl;
//	std::cout << " Iterations    = " << numberOfIterations << std::endl;
//	std::cout << " Metric value  = " << bestValue << std::endl;
//
//	TransformType::Pointer finalTransform = TransformType::New();
//
//	finalTransform->SetFixedParameters(registration->GetOutput()->Get()->GetFixedParameters());
//	finalTransform->SetParameters(finalParameters);
//
//	// Software Guide : BeginCodeSnippet
//	TransformType::MatrixType matrix = finalTransform->GetMatrix();
//	TransformType::OffsetType offset = finalTransform->GetOffset();
//	std::cout << "Matrix = " << std::endl << matrix << std::endl;
//	std::cout << "Offset = " << std::endl << offset << std::endl;
//
//	typedef itk::ResampleImageFilter<
//		MovingImageType,
//		FixedImageType >    ResampleFilterType;
//
//	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
//
//	resampler->SetTransform(finalTransform);
//	resampler->SetInput(movingImageReader->GetOutput());
//
//	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
//
//	resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
//	resampler->SetOutputOrigin(fixedImage->GetOrigin());
//	resampler->SetOutputSpacing(fixedImage->GetSpacing());
//	resampler->SetOutputDirection(fixedImage->GetDirection());
//	//resampler->SetDefaultPixelValue( 100 );
//
//	typedef  unsigned char                                          OutputPixelType;
//	typedef itk::Image< OutputPixelType, Dimension >                OutputImageType;
//	//typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;
//	typedef itk::ImageFileWriter< FixedImageType >                 WriterType;
//
//	WriterType::Pointer      writer = WriterType::New();
//	//CastFilterType::Pointer  caster =  CastFilterType::New();
//
//	writer->SetFileName(Outputfilename);
//
//	//caster->SetInput( resampler->GetOutput() );
//	writer->SetInput(resampler->GetOutput());
//	try
//	{
//		std::cout << "Writing transformed CBCT\n";
//		timer.Start("WriteCBCT");
//		writer->Update();
//		timer.Stop("WriteCBCT");
//	}
//	catch (itk::ExceptionObject& error)
//	{
//		return EXIT_FAILURE;
//	}
//
//	//std::ofstream parametersFile;
//	//parametersFile.open(OutputTransformfilename);
//	parametersFile << std::endl << std::endl;
//	parametersFile << "Result = " << std::endl;
//	parametersFile << " versor X      = " << versorX << std::endl;
//	parametersFile << " versor Y      = " << versorY << std::endl;
//	parametersFile << " versor Z      = " << versorZ << std::endl;
//	parametersFile << " Translation X = " << finalTranslationX << std::endl;
//	parametersFile << " Translation Y = " << finalTranslationY << std::endl;
//	parametersFile << " Translation Z = " << finalTranslationZ << std::endl;
//	parametersFile << " Iterations    = " << numberOfIterations << std::endl;
//	parametersFile << " Metric value  = " << bestValue << std::endl;
//	parametersFile << "\nParameters :\n";
//	parametersFile << finalParameters << std::endl;
//	//parametersFile << "Timing :\n";
//	//parametersFile << timer.
//	//parametersFile << registration->GetOutput()->Get()->get
//
//	parametersFile.close();
//	//timer.ExpandedReport();
//	timer.Report();
//	return EXIT_SUCCESS;
//
//}
//
//bool RegistrationClass<class T>::Initialize()
//{
//	fixedImageReader->SetFileName(fixedImagefilename);
//	movingImageReader->SetFileName(movingImagefilename);
//	try
//	{
//		std::cout << "Reading CT" << std::endl;
//		timer.Start("ReadCT");
//		fixedImageReader->Update();
//		std::cout << "Done\n";
//		timer.Stop("ReadCT");
//	}
//	catch (itk::ExceptionObject& err)
//	{
//		std::cerr << "Exception caught" << std::endl;
//		std::cerr << err << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	try
//	{
//		std::cout << "Reading CBCT" << std::endl;
//		timer.Start("ReadCBCT");
//		movingImageReader->Update();
//		std::cout << "Done\n";
//		timer.Stop("ReadCBCT");
//
//	}
//	catch (itk::ExceptionObject& err)
//	{
//		std::cerr << "Exception caught" << std::endl;
//		std::cerr << err << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	if (fixed_mask)
//	{
//		fixedMaskReader->SetFileName(movingMaskfilename);
//		try
//		{
//			std::cout << "Reading Fixed Image Mask" << std::endl;
//			timer.Start("ReadCBCTMask");
//			fixedMaskReader->Update();
//			std::cout << "Done\n";
//			timer.Stop("ReadCBCTMask");
//
//		}
//		catch (itk::ExceptionObject& err)
//		{
//			std::cerr << "Exception caught" << std::endl;
//			std::cerr << err << std::endl;
//			return EXIT_FAILURE;
//		}
//		spatialObjectFixedMask->SetImage(fixedMaskReader->GetOutput());
//		try
//		{
//			//std::cout << "Reading Moving Image Mask" << std::endl;
//			spatialObjectFixedMask->Update();
//			//std::cout << "Done\n";
//
//		}
//		catch (itk::ExceptionObject& err)
//		{
//			std::cerr << "Exception caught" << std::endl;
//			std::cerr << err << std::endl;
//			return EXIT_FAILURE;
//		}
//		if (debug)
//			spatialObjectFixedMask->Print(std::cout);
//		metric->SetFixedImageMask(spatialObjectFixedMask);
//	}
//
//	if (moving_mask)
//	{
//		movingMaskReader->SetFileName(movingMaskfilename);
//		try
//		{
//			std::cout << "Reading Moving Image Mask" << std::endl;
//			timer.Start("ReadCBCTMask");
//			movingMaskReader->Update();
//			std::cout << "Done\n";
//			timer.Stop("ReadCBCTMask");
//
//		}
//		catch (itk::ExceptionObject& err)
//		{
//			std::cerr << "Exception caught" << std::endl;
//			std::cerr << err << std::endl;
//			return EXIT_FAILURE;
//		}
//		spatialObjectMovingMask->SetImage(movingMaskReader->GetOutput());
//		try
//		{
//			//std::cout << "Reading Moving Image Mask" << std::endl;
//			spatialObjectMovingMask->Update();
//			//std::cout << "Done\n";
//
//		}
//		catch (itk::ExceptionObject& err)
//		{
//			std::cerr << "Exception caught" << std::endl;
//			std::cerr << err << std::endl;
//			return EXIT_FAILURE;
//		}
//		if (debug)
//			spatialObjectMovingMask->Print(std::cout);
//		metric->SetMovingImageMask(spatialObjectMovingMask);
//	}
//
//
//	if (verbose)
//		std::cout << "Set Images to Registration method\n";
//	registration->SetFixedImage(fixedImageReader->GetOutput());
//	registration->SetMovingImage(movingImageReader->GetOutput());
//	if (verbose)
//		std::cout << "Done\n";
//
//	registration->SetMetric(metric);
//	registration->SetOptimizer(optimizer);
//
//	if (verbose)
//		std::cout << "Create Initial transform\n";
//	initializer->SetTransform(initialTransform);
//	initializer->SetFixedImage(fixedImageReader->GetOutput());
//	initializer->SetMovingImage(movingImageReader->GetOutput());
//	initializer->GeometryOn();
//	if (verbose)
//		std::cout << "Done\n";
//	try
//	{
//		if (verbose)
//			std::cout << "Update Initializer\n";
//		initializer->InitializeTransform();
//		if (verbose)
//			std::cout << "Done\n";
//	}
//	catch (itk::ExceptionObject& err)
//	{
//		std::cerr << "Exception caught" << std::endl;
//		std::cerr << err << std::endl;
//		return EXIT_FAILURE;
//	}
//	if (verbose)
//		std::cout << "Set Initial rotation \n";
//	axis[0] = 0.0;
//	axis[1] = 0.0;
//	axis[2] = 1.0;
//	const double angle = 0;
//	rotation.Set(axis, angle);
//	initialTransform->SetRotation(rotation);
//	if (verbose)
//		std::cout << "Done \n";
//
//	registration->SetInitialTransform(initialTransform);
//	OptimizerScalesType optimizerScales(initialTransform->GetNumberOfParameters());
//
//	optimizerScales[0] = rotationScale;
//	optimizerScales[1] = rotationScale;
//	optimizerScales[2] = rotationScale;
//	optimizerScales[3] = translationScale;
//	optimizerScales[4] = translationScale;
//	optimizerScales[5] = translationScale;
//	optimizer->SetScales(optimizerScales);
//	optimizer->SetNumberOfIterations(200);
//
//	/* POWELL */
//
//	optimizer->SetMaximumIteration(10);
//	optimizer->SetMaximumLineIteration(4); // for Powell's method
//	optimizer->SetStepLength(4.0);
//	optimizer->SetStepTolerance(0.02);
//	optimizer->SetValueTolerance(0.001);
//	optimizer->DebugOn();
//
//	// Create the Command observer and register it with the optimizer.
//	//
//	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
//	optimizer->AddObserver(itk::IterationEvent(), observer);
//
//	this->SetLevels();
//	return EXIT_SUCCESS;
//
//}
//
//bool RegistrationClass<class T>::SetLevels()
//{
//	this->shrinkFactorsPerLevel.SetSize(this->numberOfLevels);
//	this->smoothingSigmasPerLevel.SetSize(this->numberOfLevels);
//
//	for (unsigned int i = 0; i < this->numberOfLevels; i++)
//	{
//		this->shrinkFactorsPerLevel[i] = 1;
//		this->smoothingSigmasPerLevel[i] = 0;
//	}
//	registration->SetNumberOfLevels(this->numberOfLevels);
//	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
//	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
//	return EXIT_SUCCESS;
//}
