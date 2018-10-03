#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cytolib/CytoSet.hpp>

#include "fixture.hpp"
using namespace cytolib;
struct CFFixture{
	CFFixture()
	{

		file_path = "../flowWorkspace/wsTestSuite/curlyQuad/example1/A1001.001.fcs";
		fr = MemCytoFrame(file_path, config);
		fr.read_fcs();

		string tmp = generate_temp_filename();
		fr.write_h5(tmp);
		fr_h5.reset(new H5CytoFrame(tmp));
	};

	~CFFixture(){

	};
	MemCytoFrame fr;
	unique_ptr<H5CytoFrame> fr_h5;
   string file_path;
   FCS_READ_PARAM config;
};

BOOST_FIXTURE_TEST_SUITE(CytoFrame_test,CFFixture)

BOOST_AUTO_TEST_CASE(subset_by_cols)
{
	vector<string> channels = fr.get_channels();
	vector<string> markers = fr.get_markers();
	vector<string> sub_channels = {channels[3], channels[1]};
	vector<string> sub_markers = {markers[3], markers[1]};
	CytoFrameView cr_new (CytoFramePtr(new MemCytoFrame(fr)));
	cr_new.cols_(sub_channels, ColType::channel);

	BOOST_CHECK_EQUAL(cr_new.n_cols(), sub_channels.size());

	vector<string> markers_new = cr_new.get_markers();
	BOOST_CHECK_EQUAL_COLLECTIONS(markers_new.begin(), markers_new.end(), sub_markers.begin(), sub_markers.end());

	EVENT_DATA_VEC mat = cr_new.get_data();
	BOOST_CHECK_EQUAL(mat.n_cols, sub_channels.size());//check the size of data matrix
	BOOST_CHECK_EQUAL(fr.n_cols(), channels.size());//original object is not modified

}
BOOST_AUTO_TEST_CASE(subset_by_rows)
{
	unsigned nEvent = fr.n_rows();

	CytoFrameView cr_new (CytoFramePtr(new MemCytoFrame(fr)));
	vector<unsigned> row_idx = {1,3,7};
	cr_new.rows_(row_idx);
	BOOST_CHECK_EQUAL(cr_new.n_rows(), 3);
	BOOST_CHECK_EQUAL(fr.n_rows(), nEvent);

	EVENT_DATA_VEC mat = cr_new.get_data();
	BOOST_CHECK_EQUAL(mat.n_rows, 3);

	CytoFrameView fr_copy = cr_new.copy();
	BOOST_CHECK_EQUAL(fr_copy.n_rows(), 3);

}
BOOST_AUTO_TEST_CASE(set_channel)
{
	vector<string> channels = fr.get_channels();

	MemCytoFrame fr1 = fr;
	string newname = "test";
	fr1.set_channel(channels[2], newname);
	string key = fr1.get_keyword("$P3N");
	BOOST_CHECK_EQUAL(fr1.get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);

}
BOOST_AUTO_TEST_CASE(shallow_copy)
{
	CytoFramePtr fr_orig = fr_h5->copy();
	H5CytoFrame fr1 = *(dynamic_cast<H5CytoFrame*>(fr_orig.get()));

	//update meta data
	string oldname = fr1.get_channels()[2];
	string newname = "test";
	fr1.set_channel(oldname, newname);
	string key = fr1.get_keyword("$P3N");
	BOOST_CHECK_EQUAL(fr1.get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);
	//meta data is not changed for original cp
	BOOST_CHECK_EQUAL(fr_orig->get_channels()[2], oldname);
	BOOST_CHECK_EQUAL(fr_orig->get_keyword("$P3N"), oldname);

	//update data
	EVENT_DATA_VEC dat = fr1.get_data();
	dat[100] = 100;
	fr1.set_data(dat);
	BOOST_CHECK_CLOSE(fr1.get_data()[100], dat[100], 1e-6);//fr1 is updated
	BOOST_CHECK_CLOSE(fr_orig->get_data()[100], dat[100], 1e-6);//original cp is also updated
}

BOOST_AUTO_TEST_CASE(deep_copy)
{
	//deep cp
	CytoFramePtr fr1 = fr_h5->copy();

	//update meta data
	string oldname = fr1->get_channels()[2];
	string newname = "test";
	fr1->set_channel(oldname, newname);
	string key = fr1->get_keyword("$P3N");
	//update data
	EVENT_DATA_VEC dat = fr1->get_data();
	float old_val = dat[100];
	dat[100] = 100;
	fr1->set_data(dat);
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);

	//original cp is not affected by any change in meta and event data
	BOOST_CHECK_EQUAL(fr_h5->get_channels()[2], oldname);
	BOOST_CHECK_EQUAL(fr_h5->get_keyword("$P3N"), oldname);

	BOOST_CHECK_CLOSE(fr1->get_data()[100], dat[100], 1e-6);
	BOOST_CHECK_CLOSE(fr_h5->get_data()[100], old_val, 1e-6);

	//copy_realized

}

BOOST_AUTO_TEST_CASE(CytoFrameView_copy)
{
	//
	CytoFrameView fr1(fr_h5->copy());

	//subset data
	fr1.rows_(vector<unsigned> ({1, 2, 3}));
	fr1.cols_(vector<unsigned> ({2, 4}));

	CytoFrameView fr3 = fr1.copy();// deep view cp
	CytoFrameView fr4 = fr1.copy_realized();// realized view cp
	BOOST_CHECK_EQUAL(fr3.n_rows(), 3);
	BOOST_CHECK_EQUAL(fr3.n_cols(), 2);
	BOOST_CHECK_EQUAL(fr4.n_rows(), 3);
	BOOST_CHECK_EQUAL(fr4.n_cols(), 2);

	//bypass view and directlycheck the underlying data
	BOOST_CHECK_EQUAL(fr3.get_cytoframe_ptr()->n_rows(), 10045);
	BOOST_CHECK_EQUAL(fr3.get_cytoframe_ptr()->n_cols(), 9);
	BOOST_CHECK_EQUAL(fr4.get_cytoframe_ptr()->n_rows(), 3);
	BOOST_CHECK_EQUAL(fr4.get_cytoframe_ptr()->n_cols(), 2);

}
BOOST_AUTO_TEST_SUITE_END()
