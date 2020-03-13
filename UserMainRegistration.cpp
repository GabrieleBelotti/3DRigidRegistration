#include "_gbe_3DRegistration.cpp"

int main(int argc, char *argv[])
{
	std::cout << "Class Declaration\n";
	std::cout << argv[0] << " - " << argv[1] << " - Parameters " << argc << std::endl;
	_3DRegistration *Registration = new _3DRegistration(argc, argv);
	std::cout << "Class Declared\n";
	if (strcmp(argv[1], "1") == 0)
	{
		std::cout << "Mutual Information selected as metric\n";
		Registration->Initialize< _3DRegistration::MIMetricType>();
	}
	else
	{
		std::cout << "Mean Squares selected as metric\n";
		Registration->Initialize< _3DRegistration::MeanSquaresMetricType>();
	}

	std::cout << "Run registration function\n";
	Registration->StartRegistration();
	std::cout << "All done\n";

	if (Registration) delete(Registration);
}
