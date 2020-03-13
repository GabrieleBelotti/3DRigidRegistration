
#include "gdcmReader.h"
#include "gdcmGlobal.h"
#include "gdcmDicts.h"
#include "gdcmDict.h"
#include "gdcmAttribute.h"
#include "gdcmStringFilter.h"
#include "gdcmItem.h"
#include "gdcmSequenceOfItems.h"

void IsocenterSearch(char *argv, unsigned int &count, double Isocenter[][3]);

//int main(int argc, char *argv[])
//{
//	unsigned int count = 0;
//	double Isocenter[5][3];
//	IsocenterSearch(argc, argv, count,Isocenter);
//	std::cout << "Isocenter: " << Isocenter[0][1] << " " << Isocenter[0][1] << " " << Isocenter[0][2] << std::endl;
//	std::cout << "Number of Isocenters found : " << count << std::endl;
//	return EXIT_SUCCESS;
//}

void IsocenterSearch(char *argv, unsigned int &count, double Isocenter[][3])
{
	const char *filename = argv;
	
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
	return;

	//gdcm::Reader reader;
	//reader.SetFileName(filename);
	//if (!reader.Read())
	//{
	//	std::cerr << "Could not read: " << filename << std::endl;
	//	return 1;
	//}

	//// The output of gdcm::Reader is a gdcm::File
	//gdcm::File &file = reader.GetFile();

	//// the dataset is the the set of element we are interested in:
	//gdcm::DataSet &ds = file.GetDataSet();


	//
	//using namespace gdcm;

	//Attribute<0x300a, 0x012c> IsocenterPosition;
	//IsocenterPosition.Set(ds);
	//std::cout <<"Isocenter " << IsocenterPosition.GetValue()<< std::endl;

	//Tag tPatientsName;
	////const DictEntry &de2 =
	//pubdict.GetDictEntryByName("Patient's Name", tPatientsName);

	//std::cout << "Found: " << tPatientsName << std::endl;

	//if (ds.FindDataElement(tPatientsName))
	//{
	//	gdcm::StringFilter sf;
	//	sf.SetFile(file);
	//	std::cout << "Name attribute Value as String: " << sf.ToString(tPatientsName) << std::endl;

	//	std::pair<std::string, std::string> pss
	//		= sf.ToStringPair(tPatientsName);
	//	std::cout << "Attribute Name Checked: " << pss.first << std::endl;
	//	std::cout << "Attribute Value (string): " << pss.second << std::endl;

	//}

	//Tag tIsocenterPosition;
	//pubdict.GetDictEntryByName("Isocenter Position", tIsocenterPosition);
	//
	//std::cout << "Found: " << tIsocenterPosition << std::endl;
	//std::cout << "Did we find its position? " << ds.FindDataElement(tIsocenterPosition) << std::endl;

	//Tag tIonBeamSequence(0x300a, 0x03a2);

	//const DataElement &NestedAttributesSequence_1 = ds.GetDataElement(tIonBeamSequence);

	//pubdict.GetDictEntryByName("Ion Beam Sequence", tIonBeamSequence);
	//std::cout << "Found: " << tIonBeamSequence << std::endl;
	//std::cout << "Did we find its position? " << ds.FindDataElement(tIonBeamSequence) << std::endl;

	//SequenceOfItems *sqi_1 = NestedAttributesSequence_1.GetValueAsSQ();
	//const unsigned int &nItems_1 = sqi_1->GetNumberOfItems();
	//std::cout << "Number of items, first layer: " << nItems_1 <<std::endl;
	//
	////for (int i = 1; i < 4; i++)
	////{
	////	const Item &it = sqi->GetItem(i);
	////	std::cout << it << std::endl;
	////}
	//const Tag tIonControlPointSequence(0x300a,0x03a8);
	//const Item &it_1 = sqi_1->GetItem(1);
	//const DataSet &ds_2 = it_1.GetNestedDataSet();
	//if (ds_2.FindDataElement(tIonControlPointSequence))
	//{
	//	std::cout << "Found second layer\n";
	//	const DataElement&  NestedAttributesSequence_2 = ds_2.GetDataElement(tIonControlPointSequence);
	//	SequenceOfItems *sqi_2 = NestedAttributesSequence_2.GetValueAsSQ();
	//	sqi_2->Print(std::cout);
	//	std::cout << "\nNumber of items, second layer: " <<sqi_2->GetNumberOfItems() <<std::endl;
	//	const Item &it_2 = sqi_2->GetItem(1);
	//	const DataElement de_Isocenter = it_2.GetDataElement(tIsocenterPosition);
	//	IsocenterPosition.SetFromDataElement(de_Isocenter);
	//	const double * Isocenter = IsocenterPosition.GetValues();
	//	//de_Isocenter.GetValueAsSQ
	//	std::cout << "Isocenter: " << de_Isocenter.GetValue() << std::endl;
	//	std::cout << "Isocenter: " << Isocenter[0] << " " << Isocenter[1] << " " << Isocenter[2] << std::endl;
	//}


	//return EXIT_SUCCESS;
}


//double * TagSearch(gdcm::Dict &pubd, const char * tag, const char *filename)
//{
//	double Isocenter[3][5];
//	unsigned int count=0;
//	gdcm::Reader reader;
//	reader.SetFileName(filename);
//	if (!reader.Read())
//	{
//		std::cerr << "Could not read: " << filename << std::endl;
//		return;
//	}
//
//	gdcm::File &file = reader.GetFile();
//	gdcm::DataSet &ds = file.GetDataSet();
//	gdcm::Tag tTag;
//	pubd.GetDictEntryByName(tag, tTag);
//	//const uint16_t Element = tTag.GetElement();
//	const uint16_t Group = 0x300a;
//	const uint16_t Element = 0x012c;
//	gdcm::Attribute<Group, Element> IsocenterPosition;
//
//	gdcm::Tag tIonBeamSequence(0x300a, 0x03a2);
//	const gdcm::DataElement &NestedAttributesSequence_1 = ds.GetDataElement(tIonBeamSequence);
//
//	pubd.GetDictEntryByName("Ion Beam Sequence", tIonBeamSequence);
//	std::cout << "Found: " << tIonBeamSequence << std::endl;
//	std::cout << "Did we find its position? " << ds.FindDataElement(tIonBeamSequence) << std::endl;
//	gdcm::SequenceOfItems *sqi_1 = NestedAttributesSequence_1.GetValueAsSQ();
//	const unsigned int &nItems_1 = sqi_1->GetNumberOfItems();
//	const gdcm::Tag tIonControlPointSequence(0x300a, 0x03a8);
//	for (int i = 1; i < nItems_1 + 1; i++)
//	{
//		const gdcm::Item &it_1 = sqi_1->GetItem(i);
//		const gdcm::DataSet &ds_2 = it_1.GetNestedDataSet();
//		if (ds_2.FindDataElement(tIonControlPointSequence))
//		{
//			const gdcm::DataElement&  NestedAttributesSequence_2 = ds_2.GetDataElement(tIonControlPointSequence);
//			gdcm::SequenceOfItems *sqi_2 = NestedAttributesSequence_2.GetValueAsSQ();
//			sqi_2->Print(std::cout);
//			std::cout << "\nNumber of items, second layer: " << sqi_2->GetNumberOfItems() << std::endl;
//			const gdcm::Item &it_2 = sqi_2->GetItem(1);
//			gdcm::Tag tIsocenter(Group, Element);
//			const gdcm::DataElement de_Isocenter = it_2.GetDataElement(tIsocenter);
//			IsocenterPosition.SetFromDataElement(de_Isocenter);
//			Isocenter[0][count] = IsocenterPosition.GetValues()[0];
//			Isocenter[1][count] = IsocenterPosition.GetValues()[1];
//			Isocenter[2][count] = IsocenterPosition.GetValues()[2];
//		}
//	}
//	return *Isocenter;
//}
