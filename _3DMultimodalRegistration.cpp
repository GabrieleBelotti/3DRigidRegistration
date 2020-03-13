#include "_3DMultimodalRegistration.h"

_3DMultimodalRegistration::_3DMultimodalRegistration(int argc, char *argv[])
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

		if ((this->ok == false) && (strcmp(argv[1], "-dbg") == 0))
		{
			argc--; argv++;
			ok = true;
			debug = true;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-res") == 0))
		{
			argc--; argv++;
			ok = true;
			this->ResampleSpacing[0] = atof(argv[1]);
			argc--; argv++;
			this->ResampleSpacing[1] = atof(argv[1]);
			argc--; argv++;
			this->ResampleSpacing[2] = atof(argv[1]);
			argc--; argv++;
			resample = true;
			resolution = true;
			continue;
		}

		if ((this->ok == false) && (strcmp(argv[1], "-resample") == 0))
		{
			argc--; argv++;
			ok = true;
			shrinkFactor = atof(argv[1]);
			resample = true;
			shrinking = true;
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


_3DMultimodalRegistration::~_3DMultimodalRegistration()
{
	//registration->Delete();
	//metric->Delete();
	//optimizer->Delete();
}

bool _3DMultimodalRegistration::StartRegistration()
{
	
	if (!this->Initialize())
	{
		if (verbose)
			std::cout << "Initialized correctly\n";
	}
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
	const double finalTranslationX = finalParameters[3];
	const double finalTranslationY = finalParameters[4];
	const double finalTranslationZ = finalParameters[5];
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
	std::cout << " Translation X    = " << finalTranslationX << std::endl;
	std::cout << " Translation Y    = " << finalTranslationY << std::endl;
	std::cout << " Translation Z    = " << finalTranslationZ << std::endl;
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

	//FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

	resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resampler->SetOutputOrigin(fixedImage->GetOrigin());
	resampler->SetOutputSpacing(fixedImage->GetSpacing());
	resampler->SetOutputDirection(fixedImage->GetDirection());
	resampler->SetDefaultPixelValue( -1000 );

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

//bool _3DMonomodalRegistration::Resample(FixedImageType::Pointer InputImage, double shrinkfactor, FixedImageType::Pointer &OutputImage)
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

bool _3DMultimodalRegistration::Resample(FixedImageType::Pointer InputImage, FixedImageType::SpacingType OutputSpacing , FixedImageType::Pointer &OutputImage)
{
	ResampleFilterType::Pointer downsampler = ResampleFilterType::New();
	TransformType::Pointer IdentityTransform = TransformType::New();
	IdentityTransform->SetIdentity();
	downsampler->SetTransform(IdentityTransform);
	FixedImageType::SizeType size = InputImage->GetLargestPossibleRegion().GetSize();
	FixedImageType::SpacingType InputSpacing = InputImage->GetSpacing();

	size[0] = size[0] / (OutputSpacing[0] / InputSpacing[0]);
	size[1] = size[1] / (OutputSpacing[1] / InputSpacing[1]);
	size[2] = size[2] / (OutputSpacing[2] / InputSpacing[2]);
	downsampler->SetSize(size);
	//downsampler->SetSize(InputImage->GetLargestPossibleRegion().GetSize());
	downsampler->SetInput(InputImage);
	downsampler->SetOutputOrigin(InputImage->GetOrigin());
	downsampler->SetOutputSpacing(OutputSpacing);
	downsampler->SetOutputDirection(InputImage->GetDirection());
	downsampler->SetDefaultPixelValue(-1000);
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
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << error << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

//bool _3DMonomodalRegistration::Initialize(FixedImageType::Pointer &fixedImage, MovingImageType::Pointer &movingImage)
bool _3DMultimodalRegistration::Initialize()
{
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


	if (resample)
	{
		FixedImageType::Pointer fixedResampledImage;
		MovingImageType::Pointer movingResampledImage;

		if (resolution)
		{
			if (verbose)
				std::cout << "Resampling by resolution given\n";
			if (!Resample(fixedImageReader->GetOutput(), ResampleSpacing, fixedResampledImage))
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
			if (!Resample(movingImageReader->GetOutput(), ResampleSpacing, movingResampledImage))
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
		else if (shrinking)
		{
			if (verbose)
				std::cout << "Resampling by shrinking given\n";
			if (!Resample(fixedImageReader->GetOutput(), shrinkFactor, fixedResampledImage))
			{
				if (verbose)
					std::cout << "Assign shrinked image - fixed\n";
				//fixedImage = fixedResampledImage;
			}
			if (!Resample(movingImageReader->GetOutput(), shrinkFactor, movingResampledImage))
			{
				if (verbose)
					std::cout << "Assign shrinked image - moving\n";
				//movingImage = movingResampledImage;
			}
		}

		WriterType::Pointer      writer_dbg = WriterType::New();
		//CastFilterType::Pointer  caster =  CastFilterType::New();

		writer_dbg->SetFileName("Debug Resampling.mha");

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
	else
	{
		fixedImage = fixedImageReader->GetOutput();
		movingImage = movingImageReader->GetOutput();
	}

	if (verbose)
		std::cout << "Set Images to Registration method\n";
	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImage);
	if (verbose)
		std::cout << "Done\n";

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	
	if (verbose)
		std::cout << "Create Initial transform\n";
	initializer->SetTransform(initialTransform);
	initializer->SetFixedImage(fixedImage);
	initializer->SetMovingImage(movingImage);
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
		std::cerr << "Exception caught on initializer" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	if (verbose)
		std::cout << "Set Initial rotation \n";
	//axis[0] = 0.0;
	//axis[1] = 0.0;
	//axis[2] = 1.0;
	//const double angle = 0;
	//rotation.Set(axis, angle);
	initialTransform->SetIdentity();
	
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
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();  // Check if persistant
	optimizer->AddObserver(itk::IterationEvent(), observer);

	this->SetLevels();

	return EXIT_SUCCESS;

}

bool _3DMultimodalRegistration::SetLevels()
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