#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cytolib/MemCytoFrame.hpp>
//#include "test_header.hpp"
//float gTol = 0.05;

//unsigned short myTestPolymorphism(){
//	gate * g= NULL;
//
//	rectGate rectg =rectGate();
//	g=&rectg;
//
//	rectGate * newG = dynamic_cast<rectGate*>(g);
//	return newG->getType();
//
//}
//BOOST_AUTO_TEST_SUITE(Polymorph)
//BOOST_AUTO_TEST_CASE(gateCastDown)
//{
//	BOOST_CHECK(myTestPolymorphism() == RECTGATE);
//	BOOST_CHECK(myTestPolymorphism() != POLYGONGATE);
//
//}
//BOOST_AUTO_TEST_SUITE_END()
//
//BOOST_AUTO_TEST_SUITE(RegExp)
//BOOST_AUTO_TEST_CASE(flowChannelMatch)
//{
//	string strPattern = "[FS]SC-[AWH]";
//	boost::regex ex(strPattern);
//
//	BOOST_CHECK(boost::regex_match( "FSC",ex) == false);
//	BOOST_CHECK(boost::regex_match( "FSC-A",ex) == true);
//	BOOST_CHECK(boost::regex_match( "FSC-AD",ex) == false);
//}
//
//BOOST_AUTO_TEST_SUITE_END()
struct globalFixture{
	globalFixture(){
//		cout << "Enter tolerance (e.g. 0.08):" << endl;
//		cin >> gTol;
//		gTol = 0.08;

	};
	~globalFixture(){};

};
BOOST_GLOBAL_FIXTURE(globalFixture);
struct parseFCSFixture{
	parseFCSFixture(): argc(boost::unit_test_framework::framework::master_test_suite().argc),
	           argv(boost::unit_test_framework::framework::master_test_suite().argv)
	{
		/*
		 * parse argv
		 */
		map<string, string> arg_map;
		for(int i = 1; i < argc; i++){
			string thisArg = argv[i];
			vector<string> strSplit;
			boost::split(strSplit,thisArg, boost::is_any_of("="));
			if(strSplit.size() != 2)
				throw(domain_error("invalid arguments!"));
			else{
				string argName = strSplit.at(0);
				boost::replace_first(argName, "--", "");
				string argValue = strSplit.at(1);
				arg_map[argName] = argValue;
			}
		}
		map<string, string>::iterator it;

	};



	~parseFCSFixture(){

	};
   int argc;
   char **argv;
//	testCase myTest;

};

BOOST_FIXTURE_TEST_SUITE(parseFCS,parseFCSFixture)
BOOST_AUTO_TEST_CASE(sample_1071)
{
	double start = clock();

	string filename="../flowCore/misc/sample_1071.001";
	FCS_READ_PARAM config;
	config.data.num_threads = 10;
	MemCytoFrame cytofrm(filename.c_str(), config,false);
	double runtime = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
		cout << runtime << endl;
	BOOST_CHECK_EQUAL(cytofrm.nCol(), 8);
	BOOST_CHECK_EQUAL(cytofrm.nRow(), 23981);
//	BOOST_CHECK_EQUAL_COLLECTIONS(myTest.isEqual.begin(), myTest.isEqual.end(),isTrue.begin(), isTrue.end());

}
BOOST_AUTO_TEST_CASE(double_precision)
{
	double start = clock();

	string filename="../flowCore/misc/double_precision/wishbone_thymus_panel1_rep1.fcs";
	FCS_READ_PARAM config;
	MemCytoFrame cytofrm(filename.c_str(), config,false);
	double runtime = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
	cout << runtime << endl;
	BOOST_CHECK_EQUAL(cytofrm.nCol(), 35);
	BOOST_CHECK_EQUAL(cytofrm.nRow(), 250170);
//	BOOST_CHECK_EQUAL_COLLECTIONS(myTest.isEqual.begin(), myTest.isEqual.end(),isTrue.begin(), isTrue.end());

}
BOOST_AUTO_TEST_CASE(multidata)
{
	double start = clock();

	string filename="../flowCore/misc/multi-datasegment.fcs";
	FCS_READ_PARAM config;
//	MemCytoFrame cytofrm(filename.c_str(), config,false);
//	double runtime = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
//	cout << runtime << endl;
//	BOOST_CHECK_EQUAL(cytofrm.nCol(), 10);
//	BOOST_CHECK_EQUAL(cytofrm.nRow(), 1244);

	config.header.nDataset = 10;
	MemCytoFrame cytofrm(filename.c_str(), config,false);
	BOOST_CHECK_EQUAL(cytofrm.nRow(), 955);
}

BOOST_AUTO_TEST_SUITE_END()
