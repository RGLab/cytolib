#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cytolib/MemCytoFrame.hpp>
#include <cytolib/H5CytoFrame.hpp>
#include "fixture.hpp"
#include <cytolib/global.hpp>
using namespace cytolib;

BOOST_FIXTURE_TEST_SUITE(parseFCS,parseFCSFixture)
BOOST_AUTO_TEST_CASE(sample_1071)
{
//	double start = gettime();//clock();

	string filename="../flowCore/misc/sample_1071.001";
	FCS_READ_PARAM config;
	config.data.num_threads = 4;
	//subset
	config.data.which_lines = {10};
	MemCytoFrame cf1(filename.c_str(), config);
	cf1.read_fcs();
	BOOST_CHECK_EQUAL(cf1.n_rows(), 10);
	config.data.which_lines = {10, 12};
	MemCytoFrame cf2(filename.c_str(), config);
	cf2.read_fcs();
	BOOST_CHECK_EQUAL(cf2.n_rows(), 2);

	config.data.which_lines = {};
	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();

	//	double runtime = (gettime() - start);// / (double)(CLOCKS_PER_SEC / 1000);
//		cout << runtime << endl;
	BOOST_CHECK_EQUAL(cytofrm.n_cols(), 8);
	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 23981);
//	BOOST_CHECK_EQUAL_COLLECTIONS(myTest.isEqual.begin(), myTest.isEqual.end(),isTrue.begin(), isTrue.end());

	string h5file = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
	cytofrm.write_h5(h5file);
	H5CytoFrame h5fr(h5file);
	BOOST_CHECK_EQUAL(h5fr.n_cols(), 8);
	BOOST_CHECK_EQUAL(h5fr.n_rows(), 23981);
//	for(auto c: h5fr.getChannels())
//		cout << c << endl;
//	for(auto c: h5fr.getKeywords())
//			cout << c.first << ";" << c.second << endl;

}
BOOST_AUTO_TEST_CASE(double_precision)
{
	double start = gettime();

	string filename="../flowCore/misc/double_precision/wishbone_thymus_panel1_rep1.fcs";
	FCS_READ_PARAM config;
	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();
	double runtime = (gettime() - start);
	cout << runtime << endl;
	BOOST_CHECK_EQUAL(cytofrm.n_cols(), 35);
	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 250170);
//	BOOST_CHECK_EQUAL_COLLECTIONS(myTest.isEqual.begin(), myTest.isEqual.end(),isTrue.begin(), isTrue.end());

}
BOOST_AUTO_TEST_CASE(multidata1)
{
//	double start = gettime();

	string filename="../flowCore/misc/multi-datasegment.fcs";
	FCS_READ_PARAM config;
//	MemCytoFrame cytofrm(filename.c_str(), config,false);
//	double runtime = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
//	cout << runtime << endl;
//	BOOST_CHECK_EQUAL(cytofrm.nCol(), 10);
//	BOOST_CHECK_EQUAL(cytofrm.nRow(), 1244);

	config.header.nDataset = 10;
	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();

	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 955);

	//pd
	cytofrm.set_pheno_data("colA", "a1");
	BOOST_CHECK_EQUAL(cytofrm.get_pheno_data("colB"), "");
	BOOST_CHECK_EQUAL(cytofrm.get_pheno_data("colA"), "a1");

	string h5file = "/loc/no-backup/mike/shared/test.h5";
	cytofrm.write_h5(h5file);
	H5CytoFrame fr1(h5file);
	BOOST_CHECK_EQUAL(fr1.get_pheno_data("colA"), "a1");
}
BOOST_AUTO_TEST_CASE(multidata2)
{
//	double start = gettime();

	string filename="../flowCore/misc/multi_data_segment.LMD";
	FCS_READ_PARAM config;
//	MemCytoFrame cytofrm(filename.c_str(), config,false);
//	double runtime = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
//	cout << runtime << endl;
//	BOOST_CHECK_EQUAL(cytofrm.nCol(), 10);
//	BOOST_CHECK_EQUAL(cytofrm.nRow(), 1244);

//	config.header.nDataset = 10;
	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();

	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 53691);
}

BOOST_AUTO_TEST_CASE(oddbitwidth)
{


	string filename="../flowCore/misc/Sample 2.fcs";
	FCS_READ_PARAM config;

	config.header.nDataset = 3;
	double start = gettime();//clock();
	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();

	cout << gettime() - start << endl;

	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 10000);
	BOOST_CHECK_EQUAL(cytofrm.get_data()[1], 194);

}
BOOST_AUTO_TEST_CASE(mixedEndian)
{


	string filename="../flowCore/misc/mixedEndian.fcs";
	FCS_READ_PARAM config;
	double start = gettime();//clock();
	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();

	cout << gettime() - start << endl;

	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 30041);
	BOOST_CHECK_EQUAL(cytofrm.get_data()[1], 7447226);

}

BOOST_AUTO_TEST_CASE(specialDelimiter)
{


	string filename="../flowCore/misc/specialDelimiter.fcs";
	FCS_READ_PARAM config;
	double start = gettime();//clock();
	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();

	cout << gettime() - start << endl;

	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 50146);
	BOOST_CHECK_EQUAL(cytofrm.get_data()[1], 9999998);

}

BOOST_AUTO_TEST_CASE(gigantic_file)
{


	string filename="../flowCore/misc/gigantic_file.fcs";
	FCS_READ_PARAM config;
	vector<long> which_lines(1e3);
	for(auto i = 0; i < 1e3; i++)
		which_lines[i] = i;
	config.data.which_lines = which_lines;
	double start = gettime();//clock();
	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();

	cout << gettime() - start << endl;

	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 1e3);
	BOOST_CHECK_CLOSE(cytofrm.get_data()[1], 63418.3242, 1e-6);

	config.data.seed = 2;
	config.data.which_lines = {1000};
	MemCytoFrame cytofrm1(filename.c_str(), config);
	cytofrm1.read_fcs();
	BOOST_CHECK_EQUAL(cytofrm1.n_rows(), 1e3);
	BOOST_CHECK_CLOSE(cytofrm1.get_data()[1], 38579.7617, 1e-6);

}

BOOST_AUTO_TEST_CASE(newline_in_txt)
{


	string filename="../flowWorkspace/wsTestSuite/curlyQuad/example1/A1001.001.fcs";
	FCS_READ_PARAM config;

	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs_header();
	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 0);
	cytofrm.read_fcs_data();
	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 10045);
	BOOST_CHECK_CLOSE(cytofrm.get_data()[1], 9220, 1e-6);

}
BOOST_AUTO_TEST_CASE(samples_F1)
{


	string filename="../flowWorkspace/wsTestSuite/McGill/Treg/samples_F1.fcs";
	FCS_READ_PARAM config;

	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();
	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 1000000);
	BOOST_CHECK_CLOSE(cytofrm.get_data()[1], 60981.75, 1e-6);

}
BOOST_AUTO_TEST_CASE(truncated_data_section)
{


	string filename="../flowCore/misc/truncated_data_section.fcs";
	FCS_READ_PARAM config;

	MemCytoFrame cytofrm(filename.c_str(), config);

	BOOST_CHECK_EXCEPTION(cytofrm.read_fcs(), domain_error, [](const exception & ex){return string(ex.what()).find("corrupted") != string::npos;});
}
BOOST_AUTO_TEST_CASE(skip_log_transformed)
{


	string filename="../flowCore/misc/skip_non_zero_PnE_value.fcs";
	FCS_READ_PARAM config;

	MemCytoFrame cytofrm(filename.c_str(), config);
	cytofrm.read_fcs();
	BOOST_CHECK_EQUAL(cytofrm.n_rows(), 183699);
	BOOST_CHECK_CLOSE(cytofrm.get_range("Dy161Di", ColType::channel, RangeType::instrument).second, 766, 1e-6);

}
BOOST_AUTO_TEST_SUITE_END()
