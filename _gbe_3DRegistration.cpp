#pragma once
#include "_gbe_3DRegistration.h"

/*Global Variables*/
std::string OutputTransformfilename = "3DRigidRegistrationTransform.txt";
std::string Iterationfilename = "3DRigidRegistrationIt.txt";

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
	using RegistrationPointer = RegistrationType * ;
	//using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	//typedef itk::PowellOptimizerv4<double> OptimizerType;
	typedef itk::LBFGSOptimizerv4 OptimizerType;
	using OptimizerPointer = OptimizerType * ;

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
	itkNewMacro(Self);

protected:
	CommandIterationUpdate() {};

public:
	//typedef itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
	//typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double> OptimizerType;
	//typedef itk::PowellOptimizerv4<double>  OptimizerType;
	typedef itk::LBFGSOptimizerv4 OptimizerType;
	typedef   const OptimizerType *                             OptimizerPointer;
	std::ofstream windowObs;

	void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
	{
		Execute((const itk::Object *)caller, event);
	}

	void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		OptimizerPointer optimizer = static_cast<OptimizerPointer>(object);
		if (!itk::IterationEvent().CheckEvent(&event))
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

//int MultimodalRegistration(int argc, char *argv[])
//{
//	std::cout << "Not yet implemented\n";
//	return EXIT_SUCCESS;
//}

void exe_usage()
{
	std::cerr << "\n";
	std::cerr << "Usage: ImageRegistration3DClass <unsigned int registration_routine> <options> CTname CBCTname -o outputfilename(default:CBCT_ax_aligned_iso.mha)\n";
	std::cerr << "   where <options> is one or more of the following:\n\n";
	//std::cerr << "       <-A>                     Use GPU hardware (to be implemented) [default: cpu]\n";
	std::cerr << "       <-h>								Display (this) usage information\n";
	std::cerr << "       <-v>								Verbose output [default: no]\n";
	//std::cerr << "       <-resample <double>>				Downsample by a factor <double>\n";
	std::cerr << "       <--like>  							Match CBCT to CT dimensions\n";
	std::cerr << "       <-c>  								Crop fixed to moving dimensions during registration -- speeds up the process and avoids Mutual Information issues\n";
	std::cerr << "       <-res <double, double, double>>	Resample both images to <double, double, double>\n";
	std::cerr << "       <-fres <double, double, double>>	Resample fixed image to <double, double, double>\n";
	std::cerr << "       <-mres <double, double, double>>	Resample moving image to <double, double, double>\n";
	std::cerr << "       <-dpix <float>>	Set Default Pixel Value to <float>\n";
	std::cerr << "       <-mm filename>						Moving Image Mask [default: no]\n";
	std::cerr << "       <-mf filename>						Fixed Image Mask [default: no] - under construction\n";
	std::cerr << "       <-par filename>					Output parameters filename [default: yes]\n";
	std::cerr << "       <-rt filename>  					Beta version of Isocenter alignment -- not currently working on LPI images\n";
	//std::cerr << "       <-p>							Permute image convention order (coronal to axial first) [default: no]\n";
	//std::cerr << "       <-iso double double double>    Align CBCT to isocenter [default: yes]\n";
	//std::cerr << "       <-s>							Scale image intensity [default: no] - calibration work in progress \n";
	//std::cerr << "       <-dbg>                   Debugging output [default: no]\n";

	std::cerr << "       <-o filename>						Output image filename\n\n";
	std::cerr << "											by  Gabriele Belotti\n";
	std::cerr << "											gabriele.belotti@mail.polimi.it\n";
	std::cerr << "											(Politecnico di Milano)\n\n";
	exit(EXIT_FAILURE);
}

_3DRegistration::_3DRegistration()
{
}

_3DRegistration::_3DRegistration(int argc, char *argv[])
{
	std::cout << "Parsing arguments\n";
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

		if ((ok == false) && (strcmp(argv[1], "-rt") == 0) && (customized_iso == false))
		{
			argc--; argv++;
			std::cout << "RTplan reading\n";
			ok = true;
			this->RTplanFilename = argv[1];
			this->ReadIsocenter();
			//unsigned int count_iso = 0;
			//double IsocenterBase[10][3];
			//if (this->IsocenterSearch(this->RTplanFilename, count_iso, IsocenterBase))
			//{
			//	this->RTplan = false;
			//	std::cerr << "FAILED\n";
			//}
			//else
			//{
			//	this->RTplan = true;
			//	Isocenter[0] = IsocenterBase[0][0];
			//	Isocenter[1] = IsocenterBase[0][1];
			//	Isocenter[2] = IsocenterBase[0][2];
			//}
			argc--; argv++;
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

		if ((this->ok == false) && (strcmp(argv[1], "-c") == 0))
		{
			argc--; argv++;
			this->ok = true;
			this->autocrop = true;
			continue;
			std::cout << "Autocropping selected\n";
		}

		if ((ok == false) && (strcmp(argv[1], "--like") == 0))
		{
			argc--; argv++;
			ok = true;
			this->like = true;
			std::cout << "Selected fixed dimensions for output\n";
		}

		if ((this->ok == false) && (strcmp(argv[1], "-dbg") == 0))
		{
			argc--; argv++;
			ok = true;
			debug = true;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-dpix") == 0))
		{
			argc--; argv++;
			ok = true;
			//this->DefaultPixelValue = atof(argv[1]);
			this->SetDefaultPixelValue(atof(argv[1]));
			argc--; argv++;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-res") == 0))
		{
			argc--; argv++;
			ok = true;
			this->FixedImageResampleSpacing[0] = atof(argv[1]);
			this->MovingImageResampleSpacing[0] = atof(argv[1]);
			argc--; argv++;
			this->FixedImageResampleSpacing[1] = atof(argv[1]);
			this->MovingImageResampleSpacing[1] = atof(argv[1]);
			argc--; argv++;
			this->FixedImageResampleSpacing[2] = atof(argv[1]);
			this->MovingImageResampleSpacing[2] = atof(argv[1]);
			argc--; argv++;
			this->fixed_resample = true;
			this->moving_resample = true;
			this->resolution = true;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-fres") == 0))
		{
			argc--; argv++;
			ok = true;
			this->FixedImageResampleSpacing[0] = atof(argv[1]);
			argc--; argv++;
			this->FixedImageResampleSpacing[1] = atof(argv[1]);
			argc--; argv++;
			this->FixedImageResampleSpacing[2] = atof(argv[1]);
			argc--; argv++;
			this->fixed_resample = true;
			//this->resolution = true;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-mres") == 0))
		{
			argc--; argv++;
			ok = true;
			this->MovingImageResampleSpacing[0] = atof(argv[1]);
			argc--; argv++;
			this->MovingImageResampleSpacing[1] = atof(argv[1]);
			argc--; argv++;
			this->MovingImageResampleSpacing[2] = atof(argv[1]);
			argc--; argv++;
			this->moving_resample = true;
			//this->resolution = true;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-resample") == 0))
		{
			argc--; argv++;
			this->ok = true;
			this->shrinkFactor = atof(argv[1]);
			this->shrinking = true;
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
				this->RegistrationMetric = MMI;
				argc--;
				argv++;
			}
			else if (!strcmp(argv[1], "2"))
			{
				std::cout << "Monomodality selected\n";
				this->RegistrationMetric = MSE;
				argc--;
				argv++;
			}
			else if (this->fixedImagefilename.empty())
			{
				this->fixedImagefilename = argv[1];
				argc--;
				argv++;
			}

			else if (this->movingImagefilename.empty())
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
	if (this->RegistrationMetric == 0)
		this->Initialize<_3DRegistration::MIMetricType>();
	else
		this->Initialize<_3DRegistration::MeanSquaresMetricType>();
}


_3DRegistration::~_3DRegistration()
{
	//registration->Delete();
	//metric->Delete();
	//optimizer->Delete();
}

bool _3DRegistration::ReadIsocenter()
{
	bool error_isocenter = EXIT_FAILURE;
	unsigned int count_iso = 0;
	double IsocenterBase[10][3];
	if (error_isocenter = this->IsocenterSearch(this->RTplanFilename, count_iso, IsocenterBase))
	{
		this->RTplan = false;
		std::cerr << "FAILED\n";
	}
	else
	{
		this->RTplan = true;
		this->Isocenter[0] = IsocenterBase[0][0];
		this->Isocenter[1] = IsocenterBase[0][1];
		this->Isocenter[2] = IsocenterBase[0][2];
	}
	return error_isocenter;
}
bool _3DRegistration::StartRegistration()
{
	
	//if (!this->Initialize())
	//{
	//	if (verbose)
	//		std::cout << "Initialized correctly\n";
	//}
	//this->Initialize<MIMetricType>();
	if (debug)
	{
		fixedImage->Print(std::cout);
		movingImage->Print(std::cout);
	}
	try
	{
		fixedImage->Update();
	}
	catch (itk::ExceptionObject& error)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << error << std::endl;
		return EXIT_FAILURE;
	}

	try
	{
		movingImage->Update();
	}
	catch (itk::ExceptionObject& error)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << error << std::endl;
		return EXIT_FAILURE;
	}
	
	std::cout << "Updated\n";

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
		std::cerr << "ExceptionObject caught - Registration !" << std::endl;
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
	const double finalTranslation[3] = { finalParameters[3], finalParameters[4], finalParameters[5] };
	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();

	// Print out results
	//
	//TransformType::VersorType finalVersor;
	//finalVersor.SetRotationAroundX(versorX);
	//finalVersor.SetRotationAroundX(versorY);
	//finalVersor.SetRotationAroundX(versorZ);
	//finalVersor.GetAngle();

	std::cout << std::endl << std::endl;
	std::cout << "Result = " << std::endl;
	//std::cout << " versor X         = " << versorX << std::endl;
	//std::cout << " versor Y         = " << versorY << std::endl;
	//std::cout << " versor Z         = " << versorZ << std::endl;
	std::cout << " angle around X   = " << versorX/dtr << std::endl;
	std::cout << " angle around Y   = " << versorY/dtr << std::endl;
	std::cout << " angle around Z   = " << versorZ/dtr << std::endl;
	std::cout << " Translation X    = " << finalTranslation[0] << std::endl;
	std::cout << " Translation Y    = " << finalTranslation[1] << std::endl;
	std::cout << " Translation Z    = " << finalTranslation[2] << std::endl;
	std::cout << " Iterations       = " << numberOfIterations << std::endl;
	std::cout << " Metric value     = " << bestValue << std::endl;
	
	TransformType::Pointer finalTransform = TransformType::New();

	finalTransform->SetFixedParameters(registration->GetOutput()->Get()->GetFixedParameters());
	finalTransform->SetParameters(finalParameters);

	// Software Guide : BeginCodeSnippet
	TransformType::MatrixType matrix = finalTransform->GetMatrix();
	TransformType::OffsetType offset = finalTransform->GetOffset();
	std::cout << "Matrix = " << std::endl << matrix << std::endl;
	std::cout << "Offset = " << std::endl << offset << std::endl;

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();

	resampler->SetTransform(finalTransform);
	resampler->SetInput(movingImage);

	/*MovingImageType::PointType*/ FinalMovingOrigin = movingImage->GetOrigin();
	FinalMovingOrigin[0] = FinalMovingOrigin[0] - finalTranslation[0];
	FinalMovingOrigin[1] = FinalMovingOrigin[1] - finalTranslation[1];
	FinalMovingOrigin[2] = FinalMovingOrigin[2] - finalTranslation[2];

	if (this->like)
	{
		//AGGIUNGI OFFSET PER RIPORTARE L'ORIGINE AL POSTO GIUSTO: TIPO ORIGINALFIXEDORIGIN - FINALMOVINGORIGIN
		this->FinalMovingSpacing = this->OriginalFixedSpacing;
		this->FinalMovingSize = this->OriginalFixedSize;
		resampler->SetInput(movingImageReader->GetOutput());
	}
	else
	{
		this->FinalMovingSpacing = fixedImage->GetSpacing();
		this->FinalMovingSize = fixedImage->GetLargestPossibleRegion().GetSize();
	}
	resampler->SetOutputOrigin(FinalMovingOrigin);
	resampler->SetDefaultPixelValue(this->DefaultPixelValue);
	resampler->SetSize(this->FinalMovingSize);
	resampler->SetOutputSpacing(this->FinalMovingSpacing);
	resampler->SetOutputDirection(fixedImage->GetDirection());


	WriterType::Pointer      writer = WriterType::New();
	//CastFilterType::Pointer  caster =  CastFilterType::New();

	writer->SetFileName(Outputfilename);

	//caster->SetInput( resampler->GetOutput() );
	writer->SetInput(resampler->GetOutput());
	try
	{
		std::cout << "Writing registered Moving Image\n";
		timer.Start("WriteCBCT");
		writer->Update();
		timer.Stop("WriteCBCT");
	}
	catch (itk::ExceptionObject& error)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << error << std::endl;
		return EXIT_FAILURE;
	}

	//std::ofstream parametersFile;
	//parametersFile.open(OutputTransformfilename);
	parametersFile << std::endl << std::endl;
	parametersFile << "Result = " << std::endl;
	parametersFile << " versor X      = " << versorX << std::endl;
	parametersFile << " versor Y      = " << versorY << std::endl;
	parametersFile << " versor Z      = " << versorZ << std::endl;
	parametersFile << " Translation X = " << finalTranslation[0] << std::endl;
	parametersFile << " Translation Y = " << finalTranslation[1] << std::endl;
	parametersFile << " Translation Z = " << finalTranslation[2] << std::endl;
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

//bool _3DRegistration::Resample(FixedImageType::Pointer InputImage, double shrinkfactor, FixedImageType::Pointer &OutputImage)
//{
//	ResampleFilterType::Pointer downsampler = ResampleFilterType::New();
//	TransformType::Pointer IdentityTransform = TransformType::New();
//	IdentityTransform->SetIdentity();
//	downsampler->SetTransform(IdentityTransform);
//	downsampler->SetInput(InputImage);
//	downsampler->SetOutputOrigin(InputImage->GetOrigin());
//	downsampler->SetOutputSpacing(InputImage->GetSpacing() / shrinkfactor);
//	downsampler->SetOutputDirection(InputImage->GetDirection());
//	
//	try
//	{
//		if (verbose)
//			std::cout << "Resampling\n";
//		downsampler->Update();
//		OutputImage = downsampler->GetOutput();
//		OutputImage->DisconnectPipeline();
//		if (verbose)
//			std::cout << "Resampling done\n";
//	}
//	catch (itk::ExceptionObject& error)
//	{
//		std::cerr << "ExceptionObject caught !" << std::endl;
//		std::cerr << error << std::endl;
//		return EXIT_FAILURE;
//	}
//	return EXIT_SUCCESS;
//}

bool _3DRegistration::Resample(FixedImageType::Pointer InputImage, FixedImageType::SpacingType OutputSpacing , FixedImageType::Pointer &OutputImage)
{
	ResampleFilterType::Pointer downsampler = ResampleFilterType::New();

	using InterpolatorType = itk::BSplineResampleImageFunction<MovingImageType, double>;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	TransformType::Pointer IdentityTransform = TransformType::New();
	IdentityTransform->SetIdentity();
	downsampler->SetTransform(IdentityTransform);
	FixedImageType::SizeType size = InputImage->GetLargestPossibleRegion().GetSize();
	FixedImageType::SpacingType InputSpacing = InputImage->GetSpacing();

	size[0] = size[0] / (OutputSpacing[0] / InputSpacing[0]);
	size[1] = size[1] / (OutputSpacing[1] / InputSpacing[1]);
	size[2] = size[2] / (OutputSpacing[2] / InputSpacing[2]);
	downsampler->SetSize(size);
	downsampler->SetInterpolator(interpolator);
	//downsampler->SetSize(InputImage->GetLargestPossibleRegion().GetSize());
	//downsampler->SetInterpolator()
	downsampler->SetInput(InputImage);
	downsampler->SetOutputOrigin(InputImage->GetOrigin());
	downsampler->SetOutputSpacing(OutputSpacing);
	downsampler->SetOutputDirection(InputImage->GetDirection());
	downsampler->SetDefaultPixelValue(this->DefaultPixelValue);
	OutputImage = downsampler->GetOutput();
	try
	{
		if (verbose)
			std::cout << "Resampling\n";
		OutputImage->Update();
		if (verbose)
			std::cout << "Resampling done\n";
	}
	catch (itk::ExceptionObject& error)
	{
		std::cerr << "ExceptionObject caught during resampling !" << std::endl;
		std::cerr << error << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

bool _3DRegistration::Crop(FixedImageType::Pointer Image2Crop, FixedImageType::Pointer ReferenceImage, FixedImageType::Pointer & OutputImage)
{
	itk::ImageRegion<this->Dimension> FixedRegion = Image2Crop->GetLargestPossibleRegion();
	itk::ImageRegion<this->Dimension> MovingRegion = ReferenceImage->GetLargestPossibleRegion();
	
	FixedImageType::SizeType InputSize = FixedRegion.GetSize();
	MovingImageType::SizeType ReferenceSize = MovingRegion.GetSize();

	FixedImageType::SizeType OutputSize;

	FixedImageType::SpacingType InputSpacing = Image2Crop->GetSpacing();
	MovingImageType::SpacingType ReferenceSpacing = ReferenceImage->GetSpacing();

	FixedImageType::SizeType LowerCrop{ 0 };
	FixedImageType::SizeType UpperCrop{ 0 };

	FixedImageType::PointType InputOrigin = Image2Crop->GetOrigin();
	FixedImageType::PointType OutputOrigin = ReferenceImage->GetOrigin();
	
	unsigned int padding_value[this->Dimension] = { 10, 10, 5 };

	//FixedRegion.Crop(MovingRegion);
	if (FixedRegion.IsInside(MovingRegion))
	{
		for (int kk = 0; kk < this->Dimension; kk++)
		{
			//OutputSize[kk] = ReferenceSize[kk] / (InputSpacing[kk] / ReferenceSpacing[kk]);
			LowerCrop[kk] = (OutputOrigin[kk] - InputOrigin[kk]) / InputSpacing[kk]; // Check for positivity --> we need to make sure we're superimposing a subregion to the fixed image
			UpperCrop[kk] = (InputSize[kk] - (LowerCrop[kk] + ReferenceSize[kk] / (InputSpacing[kk]/ReferenceSpacing[kk]))); // Check for positivity --> we need to make sure we're superimposing a subregion to the fixed image
			OutputSize[kk] = InputSize[kk] - (LowerCrop[kk] + UpperCrop[kk]);
			/* add some space around the cropped area to avoid losing info */
			//if (LowerCrop[kk] <= 10)
			//	LowerCrop[kk] = 0;
			//else
			//	LowerCrop[kk] -= padding_value[kk];

			//if (UpperCrop[kk] <= 10)
			//	UpperCrop[kk] = 0;
			//else
			//	UpperCrop[kk] -= padding_value[kk];
		}
	}
	else
	{
		std::cerr << "Reference Image is not completely inside the Fixed one --> cropping is not supported\n";
		return EXIT_FAILURE;
	}
	std::cout << "Reference Size " << ReferenceSize << std::endl;
	std::cout << "Input Size " << InputSize << std::endl;
	std::cout << "Output Size " << OutputSize << std::endl;

	std::cout << "Lower crop " << LowerCrop << std::endl;
	std::cout << "Upper crop " << UpperCrop << std::endl;

	CropFixedFilterType::Pointer CropFixedFilter = CropFixedFilterType::New();
	CropFixedFilter->SetInput(Image2Crop);
	CropFixedFilter->SetLowerBoundaryCropSize(LowerCrop);
	CropFixedFilter->SetUpperBoundaryCropSize(UpperCrop);
	//CropFixedFilter->SetExtractionRegion(FixedRegion); //valid only using extractimagefilter
	CropFixedFilter->SetDirectionCollapseToIdentity();

	try
	{
		if (verbose)
			std::cout << "Cropping Fixed Image to Moving image size\n";
		CropFixedFilter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "Exception caught while cropping \n" 
			<< err << std::endl;
		return EXIT_FAILURE;
	}

	OutputImage = CropFixedFilter->GetOutput();

	return EXIT_SUCCESS;
}

bool _3DRegistration::ROICrop(FixedImageType::Pointer Image2Crop, FixedImageType::Pointer ReferenceImage, FixedImageType::Pointer & OutputImage)
{
	double temp; //utility for calculation

	itk::ImageRegion<this->Dimension> FixedRegion = Image2Crop->GetLargestPossibleRegion();
	itk::ImageRegion<this->Dimension> MovingRegion = ReferenceImage->GetLargestPossibleRegion();
	itk::ImageRegion<this->Dimension> OutputRegion;

	const FixedImageType::DirectionType& FixedDirection = fixedImage->GetDirection();
	MovingImageType::DirectionType& OutputDirection = MovingImageType::DirectionType::Matrix();
	OutputDirection.SetIdentity(); // why setting RAI orientation?

	using OrientationAdapterType = itk::SpatialOrientationAdapter;
	itk::SpatialOrientation::ValidCoordinateOrientationFlags FixedOrientationFlag;
	
	const MovingImageType::DirectionType& MovingDirection = movingImage->GetDirection();

	OrientationAdapterType OrientationAdapter;
	FixedOrientationFlag = OrientationAdapter.FromDirectionCosines(FixedDirection);

	FixedImageType::SizeType InputSize = FixedRegion.GetSize();
	MovingImageType::SizeType ReferenceSize = MovingRegion.GetSize();

	FixedImageType::SpacingType InputSpacing = Image2Crop->GetSpacing();
	MovingImageType::SpacingType ReferenceSpacing = ReferenceImage->GetSpacing();
	//if (FixedOrientationFlag == itk::SpatialOrientation::ValidCoordinateOrientationFlags::ITK_COORDINATE_ORIENTATION_RAI)
	//else if (OrientationFlag == itk::SpatialOrientation::ValidCoordinateOrientationFlags::ITK_COORDINATE_ORIENTATION_LPI)
	MovingImageType::SizeType OutputSize;

	FixedImageType::IndexType StartIndex{ 0 };

	FixedImageType::SizeType LowerCrop{ 0 };
	FixedImageType::SizeType UpperCrop{ 0 };

	FixedImageType::PointType InputOrigin = Image2Crop->GetOrigin();
	FixedImageType::PointType OutputOrigin = ReferenceImage->GetOrigin();

	unsigned int padding_value[this->Dimension] = { 10, 10, 5 };

	//if (FixedRegion.IsInside(MovingRegion))
	//{
		if (FixedOrientationFlag == itk::SpatialOrientation::ValidCoordinateOrientationFlags::ITK_COORDINATE_ORIENTATION_RAI)
		{

			for (int kk = 0; kk < this->Dimension; kk++)
			{
				temp = (OutputOrigin[kk] - InputOrigin[kk]) / InputSpacing[kk];
				if (temp < 0) temp = 0;
				StartIndex[kk] = temp;
				LowerCrop[kk] = temp; // Check for positivity --> we need to make sure we're superimposing a subregion to the fixed image
				temp = (InputSize[kk] - (LowerCrop[kk] + ReferenceSize[kk] / (InputSpacing[kk] / ReferenceSpacing[kk])));
				if (temp < 0) temp = 0;
				UpperCrop[kk] = temp; // Check for positivity --> we need to make sure we're superimposing a subregion to the fixed image
				OutputSize[kk] = InputSize[kk] - (LowerCrop[kk] + UpperCrop[kk]);
				//OutputSize[kk] = ReferenceSize[kk] * (ReferenceSpacing[kk] / InputSpacing[kk]);
				/* add some space around the cropped area to avoid losing info */
				//if (LowerCrop[kk] <= 10)
				//	LowerCrop[kk] = 0;
				//else
				//	LowerCrop[kk] -= padding_value[kk];

				//if (UpperCrop[kk] <= 10)
				//	UpperCrop[kk] = 0;
				//else
				//	UpperCrop[kk] -= padding_value[kk];
			}
		}
		else
		{
			temp = (-OutputOrigin[0] + InputOrigin[0]) / InputSpacing[0];
			if (temp < 0) temp = 0;
			StartIndex[0] = temp;
			temp = (-OutputOrigin[1] + InputOrigin[1]) / InputSpacing[1];
			if (temp < 0) temp = 0;
			StartIndex[1] = temp;
			temp = (OutputOrigin[2] - InputOrigin[2]) / InputSpacing[2];
			if (temp < 0) temp = 0;
			StartIndex[2] = temp;
			for (int kk = 0; kk < this->Dimension; kk++)
			{
				LowerCrop[kk] = StartIndex[kk];
				temp = (InputSize[kk] - (LowerCrop[kk] + ReferenceSize[kk] / (InputSpacing[kk] / ReferenceSpacing[kk]))); // Check for positivity --> we need to make sure we're superimposing a subregion to the fixed image
				if (temp < 0) temp = 0;
				UpperCrop[kk] = temp;
				OutputSize[kk] = InputSize[kk] - (LowerCrop[kk] + UpperCrop[kk]);
			}
		}
	//}

	//else
	//{
	//	std::cerr << "Reference Image is not completely inside the Fixed one --> cropping is not supported\n";
	//	return EXIT_FAILURE;
	//}

	std::cout << "Input Size " << InputSize << std::endl;
	std::cout << "Reference Size " << ReferenceSize << std::endl;
	std::cout << "Output Size " << OutputSize << std::endl;

	std::cout << "Start Index " << StartIndex << std::endl;

	std::cout << "Lower crop " << LowerCrop << std::endl;
	std::cout << "Upper crop " << UpperCrop << std::endl;

	//CropFixedFilterType::Pointer CropFixedFilter = CropFixedFilterType::New();
	ROIFilterType::Pointer ROIFilter = ROIFilterType::New();

	//CropFixedFilter->SetInput(Image2Crop);
	//CropFixedFilter->SetLowerBoundaryCropSize(LowerCrop);
	//CropFixedFilter->SetUpperBoundaryCropSize(UpperCrop);
	//CropFixedFilter->SetExtractionRegion(FixedRegion); //valid only using extractimagefilter
	//CropFixedFilter->SetDirectionCollapseToIdentity();

	//try
	//{
	//	if (verbose)
	//		std::cout << "Cropping Fixed Image to Moving image size\n";
	//	CropFixedFilter->Update();
	//}
	//catch(itk::ExceptionObject &err)
	//{
	//	std::cerr << "Exception caught while cropping \n" 
	//		<< err << std::endl;
	//	return EXIT_FAILURE;
	//}

	OutputRegion.SetIndex(StartIndex);
	OutputRegion.SetSize(OutputSize);
	//FixedRegion.Crop(OutputRegion);

	ROIFilter->SetRegionOfInterest(OutputRegion);
	ROIFilter->SetInput(Image2Crop);

	try
	{
		if (verbose)
			std::cout << "Extracting Fixed Image region to Moving image size\n";
		ROIFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "Exception caught while extracting \n"
			<< err << std::endl;
		return EXIT_FAILURE;
	}

	//OutputImage = CropFixedFilter->GetOutput();
	OutputImage = ROIFilter->GetOutput();

	return EXIT_SUCCESS;
}

//bool _3DRegistration::Initialize(FixedImageType::Pointer &fixedImage, MovingImageType::Pointer &movingImage)
template<typename TMetricType>
bool _3DRegistration::Initialize()
{
	/* Set Metric type to caller template */
	metric = TMetricType::New();
	/*Reading CT and CBCT images*/
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
	/*End of reading CT and CBCT images*/

	this->OriginalMovingSpacing = movingImageReader->GetOutput()->GetSpacing();
	this->OriginalMovingOrigin = movingImageReader->GetOutput()->GetOrigin();
	this->OriginalMovingSize = movingImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

	this->OriginalFixedSpacing = fixedImageReader->GetOutput()->GetSpacing();
	this->OriginalFixedOrigin = fixedImageReader->GetOutput()->GetOrigin();
	this->OriginalFixedSize = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();


	/*Reading masks if requested/available*/
	if (fixed_mask)
	{
		fixedMaskReader->SetFileName(fixedMaskfilename);
		try
		{
			std::cout << "Reading Fixed Image Mask" << std::endl;
			timer.Start("ReadCTMask");
			fixedMaskReader->Update();
			std::cout << "Done\n";
			timer.Stop("ReadCTMask");

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
	/*End of reading masks if requested/available*/

	fixedImage = fixedImageReader->GetOutput();
	movingImage = movingImageReader->GetOutput();

	if (this->moving_resample || this->fixed_resample)
	{
		
		if (verbose)
			std::cout << "Resampling by resolution given\n";
		if (this->fixed_resample)
		{

			FixedImageType::Pointer fixedResampledImage;
			
			if (!Resample(fixedImageReader->GetOutput(), FixedImageResampleSpacing, fixedResampledImage))
			{
				if (verbose)
				{
					fixedResampledImage->Print(std::cout);
					std::cout << "Assign resampled image - fixed\n";
				}
				fixedImage = fixedResampledImage;
				if (verbose)
				{
					std::cout << "Done\n";
					//fixedImage->Print(std::cout);
				}
			}
			else
			{
				std::cout << "CT resampling unsuccessful\n";
			}
			if (debug)
			{
				WriterType::Pointer      writer_dbg = WriterType::New();
				//CastFilterType::Pointer  caster =  CastFilterType::New();

				writer_dbg->SetFileName("Debug_Resampling.mha");

				//caster->SetInput( resampler->GetOutput() );
				writer_dbg->SetInput(fixedResampledImage);
				//writer_dbg->SetInput(fixedResampleImage);
				try
				{
					std::cout << "Writing CT Resampled\n";
					timer.Start("WriteCTRes");
					writer_dbg->Update();
					timer.Stop("WriteCTRes");
				}
				catch (itk::ExceptionObject& error)
				{
					std::cerr << "ExceptionObject caught !" << std::endl;
					std::cerr << error << std::endl;
					return EXIT_FAILURE;
				}
			}
		}

		if (moving_resample)
		{

			MovingImageType::Pointer movingResampledImage;

			if (!Resample(movingImageReader->GetOutput(), MovingImageResampleSpacing, movingResampledImage))
			{
				if (verbose)
					std::cout << "Assign resampled image - moving\n";
				movingImage = movingResampledImage;
				if (verbose)
					std::cout << "Done\n";
			}
			else
			{
				std::cout << "CBCT resampling unsuccessful\n";
			}
		}
		// DELETE SHRINKING ROUTINE
		//else if (shrinking)
		//{
		//	if (verbose)
		//		std::cout << "Resampling by shrinking given\n";
		//	if (!Resample(fixedImageReader->GetOutput(), shrinkFactor, fixedResampledImage))
		//	{
		//		if (verbose)
		//			std::cout << "Assign shrinked image - fixed\n";
		//		//fixedImage = fixedResampledImage;
		//	}
		//	if (!Resample(movingImageReader->GetOutput(), shrinkFactor, movingResampledImage))
		//	{
		//		if (verbose)
		//			std::cout << "Assign shrinked image - moving\n";
		//		//movingImage = movingResampledImage;
		//	}
		//}
		
	}

	// Gathering information about the images
 	const SpacingType fixedSpacing = fixedImage->GetSpacing();
	const OriginType  fixedOrigin = fixedImage->GetOrigin();
	const RegionType  fixedRegion = fixedImage->GetLargestPossibleRegion();
	const SizeType    fixedSize = fixedRegion.GetSize();

	const SpacingType movingSpacing = movingImage->GetSpacing();
	const OriginType  movingOrigin = movingImage->GetOrigin();
	const RegionType  movingRegion = movingImage->GetLargestPossibleRegion();
	const SizeType    movingSize = movingRegion.GetSize();

	if (verbose)
		std::cout << "Create Initial transform\n";
	// Aligning the moving image to isocenter if RTPlan is provided --> FIX FOR PRONE PATIENTS AKA LPI CONFIG
	if (RTplan)
	{
		TransformType::InputPointType centerMoving;

		centerMoving[0] = movingOrigin[0] + movingSpacing[0] * movingSize[0] / 2.0;
		centerMoving[1] = movingOrigin[1] + movingSpacing[1] * movingSize[1] / 2.0;
		centerMoving[2] = movingOrigin[2] + movingSpacing[2] * movingSize[2] / 2.0;

		translation = centerMoving - Isocenter;
		std::cout << "Translation to isocenter: " << translation << std::endl;
		initialTransform->SetCenter(Isocenter);

		initialTransform->SetTranslation(translation);

		IsoAlignment->SetTransform(initialTransform);
		IsoAlignment->SetInput(movingImage);
		IsoAlignment->SetOutputSpacing(movingSpacing);
		//IsoAlignment->SetOutputOrigin(fixedOrigin);
		IsoAlignment->SetOutputOrigin(movingOrigin - translation);
		IsoAlignment->SetSize(movingSize);

		try
		{
			if (verbose)
				std::cout << "Update Initializer\n";
			IsoAlignment->Update();
			if (verbose)
				std::cout << "Done\n";
		}
		catch (itk::ExceptionObject& err)
		{
			std::cerr << "Exception caught on initializer" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
		movingImage = IsoAlignment->GetOutput();
	}
	if (this->RTplan || this->autocrop)
	{
		if (this->debug)
		{
			WriterType::Pointer      writer_dbg_iso = WriterType::New();
			//CastFilterType::Pointer  caster =  CastFilterType::New();

			writer_dbg_iso->SetFileName("Debug_isocenter_alignment.mha");

			//caster->SetInput( resampler->GetOutput() );
			writer_dbg_iso->SetInput(movingImage);
			//writer_dbg->SetInput(fixedResampleImage);
			try
			{
				std::cout << "Writing isocenter aligned CBCT\n";
				timer.Start("WriteCBCTiso");
				writer_dbg_iso->Update();
				timer.Stop("WriteCBCTiso");
			}
			catch (itk::ExceptionObject& error)
			{
				std::cerr << "ExceptionObject caught !" << std::endl;
				std::cerr << error << std::endl;
				return EXIT_FAILURE;
			}
		}
		if (this->RegistrationMetric == 0 || this->RegistrationMetric == 1)
		{
			FixedImageType::Pointer CroppedFixed;
			//if (this->Crop(fixedImage, movingImage, CroppedFixed))
			if (this->ROICrop(fixedImage, movingImage, CroppedFixed))
				return EXIT_FAILURE;
			if (this->debug)
			{
				WriterType::Pointer writer_dbg_crop = WriterType::New();
				//CastFilterType::Pointer  caster =  CastFilterType::New();

				writer_dbg_crop->SetFileName("Debug_fixed_crop.mha");

				//writer_dbg_crop->SetInput(fixedImage);
				writer_dbg_crop->SetInput(CroppedFixed);
				try
				{
					std::cout << "Writing cropped CT\n";
					timer.Start("WriteCTcrop");
					writer_dbg_crop->Update();
					timer.Stop("WriteCTcrop");
				}
				catch (itk::ExceptionObject& error)
				{
					std::cerr << "ExceptionObject caught !" << std::endl;
					std::cerr << error << std::endl;
					return EXIT_FAILURE;
				}

				FixedImageReaderType::Pointer reader_dbg_crop = FixedImageReaderType::New();
				//CastFilterType::Pointer  caster =  CastFilterType::New();

				reader_dbg_crop->SetFileName("Debug_fixed_crop.mha");

				try
				{
					std::cout << "Reading cropped CT\n";
					timer.Start("ReadCTcrop");
					reader_dbg_crop->Update();
					timer.Stop("ReadCTcrop");
				}
				catch (itk::ExceptionObject& error)
				{
					std::cerr << "ExceptionObject caught while reading written cropped CT!" << std::endl;
					std::cerr << error << std::endl;
					return EXIT_FAILURE;
				}
				CroppedFixed = reader_dbg_crop->GetOutput();
			}
			fixedImage = CroppedFixed;
			if (debug)
			{
				fixedImage->Print(std::cout);
			}
		}
	}
	

	if (verbose)
		std::cout << "Set Images to Registration method\n";
	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImage);
	if (verbose)
		std::cout << "Done\n";
	
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);

	// INITIALIZER SETS THE INITIAL TRANSFORM AND CAN OVERRIDE PREVIOUS ACTIONS
	//else
	//{
	//	initialTransform->SetIdentity();

	//	initializer->SetTransform(initialTransform);
	//	initializer->SetFixedImage(fixedImage);
	//	initializer->SetMovingImage(movingImage);
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
	//		std::cerr << "Exception caught on initializer" << std::endl;
	//		std::cerr << err << std::endl;
	//		return EXIT_FAILURE;
	//	}
	//	if (verbose)
	//		std::cout << "Set Initial rotation \n";
	//	//axis[0] = 0.0;
	//	//axis[1] = 0.0;
	//	//axis[2] = 1.0;
	//	//const double angle = 0;
	//	//rotation.Set(axis, angle);

	//	if (verbose)
	//		std::cout << "Done \n";

	//	registration->SetInitialTransform(initialTransform);
	//}
	OptimizerScalesType optimizerScales(initialTransform->GetNumberOfParameters());

	optimizerScales[0] = rotationScale;
	optimizerScales[1] = rotationScale;
	optimizerScales[2] = rotationScale;
	optimizerScales[3] = translationScale;
	optimizerScales[4] = translationScale;
	optimizerScales[5] = translationScale;
	optimizer->SetScales(optimizerScales);
	optimizer->SetNumberOfIterations(200);
	optimizer->SetGradientConvergenceTolerance(1e-4);
	//optimizer->SetLineSearchAccuracy(0.5); //This makes for a more accurate line search the lower the value (default is 0.9)
	//optimizer->DoEstimateScalesOff();
	
	/* POWELL */

	//optimizer->SetMaximumIteration(10);
	//optimizer->SetMaximumLineIteration(4); // for Powell's method
	//optimizer->SetStepLength(4.0);
	//optimizer->SetStepTolerance(0.02);
	//optimizer->SetValueTolerance(0.001);
	//optimizer->DebugOn();

	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();  // Check if persistant
	optimizer->AddObserver(itk::IterationEvent(), observer);

	this->SetLevels();

	return EXIT_SUCCESS;

}

bool _3DRegistration::SetLevels()
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

bool _3DRegistration::IsocenterSearch(std::string argv, unsigned int &count, double Isocenter[][3])
{
	const char *filename = argv.c_str();

	const gdcm::Global &g = gdcm::Global::GetInstance();
	const gdcm::Dicts &dicts = g.GetDicts();
	const gdcm::Dict &pubdict = dicts.GetPublicDict();

	gdcm::Reader reader;
	reader.SetFileName(filename);
	if (!reader.Read())
	{
		std::cerr << "Could not read: " << filename << std::endl;
	}
	const char * tag = "Isocenter Position";
	gdcm::File &file = reader.GetFile();
	gdcm::DataSet &ds = file.GetDataSet();
	gdcm::Tag tTag;
	pubdict.GetDictEntryByName(tag, tTag);
	//const uint16_t Element = tTag.GetElement();
	const uint16_t Group = 0x300a;
	const uint16_t Element = 0x012c;
	gdcm::Attribute<Group, Element> IsocenterPosition;

	gdcm::Tag tIonBeamSequence(0x300a, 0x03a2);
	const gdcm::DataElement &NestedAttributesSequence_1 = ds.GetDataElement(tIonBeamSequence);

	pubdict.GetDictEntryByName("Ion Beam Sequence", tIonBeamSequence);
	std::cout << "Found: " << tIonBeamSequence << std::endl;
	std::cout << "Did we find its position? " << ds.FindDataElement(tIonBeamSequence) << std::endl;
	if (!ds.FindDataElement(tIonBeamSequence))
		return EXIT_FAILURE;
	gdcm::SequenceOfItems *sqi_1 = NestedAttributesSequence_1.GetValueAsSQ();
	const unsigned int &nItems_1 = sqi_1->GetNumberOfItems();
	const gdcm::Tag tIonControlPointSequence(0x300a, 0x03a8);
	for (unsigned int i = 1; i < nItems_1 + 1; i++)
	{
		const gdcm::Item &it_1 = sqi_1->GetItem(i);
		const gdcm::DataSet &ds_2 = it_1.GetNestedDataSet();
		if (ds_2.FindDataElement(tIonControlPointSequence))
		{
			const gdcm::DataElement&  NestedAttributesSequence_2 = ds_2.GetDataElement(tIonControlPointSequence);
			gdcm::SequenceOfItems *sqi_2 = NestedAttributesSequence_2.GetValueAsSQ();
			//sqi_2->Print(std::cout);
			//std::cout << "\nNumber of items, second layer: " << sqi_2->GetNumberOfItems() << std::endl;
			const gdcm::Item &it_2 = sqi_2->GetItem(1);
			gdcm::Tag tIsocenter(Group, Element);
			const gdcm::DataElement de_Isocenter = it_2.GetDataElement(tIsocenter);
			IsocenterPosition.SetFromDataElement(de_Isocenter);
			Isocenter[count][0] = IsocenterPosition.GetValues()[0];
			Isocenter[count][1] = IsocenterPosition.GetValues()[1];
			Isocenter[count][2] = IsocenterPosition.GetValues()[2];
			std::cout << "Isocenter: " << Isocenter[count][0] << " " << Isocenter[count][1] << " " << Isocenter[count][2] << std::endl;
			count++;
		}
	}
	std::cout << "Number of Isocenters found : " << count << std::endl;
	return EXIT_SUCCESS;
}