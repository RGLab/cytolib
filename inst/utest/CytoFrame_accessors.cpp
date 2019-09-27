#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/lexical_cast.hpp>
#include <cytolib/GatingSet.hpp>
#include <cytolib/H5CytoFrame.hpp>
#include <cytolib/MemCytoFrame.hpp>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#include "fixture.hpp"
using namespace cytolib;
struct CFFixture{
	CFFixture()
	{
		//
		file_path = "../flowWorkspace/wsTestSuite/curlyQuad/example1/A1001.001.fcs";
		fr = MemCytoFrame(file_path, config);
		fr.read_fcs();
		//create h5 version
		string tmp = generate_unique_filename(fs::temp_directory_path(), "", ".h5");
//		cout << tmp << endl;
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

BOOST_AUTO_TEST_CASE(get_time_step)
{
	auto fr1 = MemCytoFrame("../flowWorkspace/output/s5a01.fcs", config);
	fr1.read_fcs();
	BOOST_CHECK_CLOSE(fr1.get_time_step("Time"), 0.0557, 0.01);

}
BOOST_AUTO_TEST_CASE(h5_vs_mem)
{
	//check if h5 version is consistent with mem
	BOOST_CHECK_EQUAL(fr.get_params().size(), fr_h5->get_params().size());
	BOOST_CHECK_EQUAL(fr.get_params().begin()->max, fr_h5->get_params().begin()->max);

}
BOOST_AUTO_TEST_CASE(flags)
{
	//get a safe cp
	string h5file = fr_h5->copy()->get_h5_file_path();
	//load it by default readonly flag
	unique_ptr<H5CytoFrame> fr1(new H5CytoFrame(h5file));
	//update meta data
	string oldname = fr1->get_channels()[2];
	string newname = "test";
	BOOST_CHECK_EXCEPTION(fr1->set_channel(oldname, newname), domain_error,
				[](const domain_error & ex) {return string(ex.what()).find("read-only") != string::npos;});

	//save error handler
	H5E_auto2_t func;
	void* client_data;
	H5::Exception::getAutoPrint(func, &client_data);
	//turn off auto print
	H5::Exception::dontPrint();

	//try to flush to disk
//	BOOST_CHECK_EXCEPTION(fr1->flush_meta(), H5::DataSetIException,
//				[](const H5::DataSetIException & ex) {return ex.getDetailMsg().find("H5Dwrite failed") != string::npos;});
	//load it again to check if meta data is intact
	fr1->load_meta();
	//no change to disk
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], oldname);
	string key = fr1->get_keyword("$P3N");
	BOOST_CHECK_EQUAL(key, oldname);

	//update data
	EVENT_DATA_VEC dat = fr1->get_data();
	float oldval = dat[100];
	float newval = 100;
	dat[100] = newval;

//	BOOST_CHECK_EXCEPTION(fr1->set_data(dat);, H5::DataSetIException,
//				[](const H5::DataSetIException & ex) {return ex.getDetailMsg().find("read-only") != string::npos;});
	BOOST_CHECK_EXCEPTION(fr1->set_data(dat), domain_error,
			[](const exception & ex) {return string(ex.what()).find("read-only") != string::npos;});

	BOOST_CHECK_CLOSE(fr1->get_data()[100], oldval, 1e-6);


	/**
	 * try another open with different flag
	 */
	//destroy it to expect succeed
	fr1.reset();
	//turn on error report
//	H5::Exception::setAutoPrint(func, client_data);
	//clear previous errors from stacks resulted from write attempts
//	H5::Exception::clearErrorStack();
	//try open again
	H5CytoFrame fr3(h5file, false);
	fr3.set_data(dat);
	BOOST_CHECK_CLOSE(fr3.get_data()[100], newval, 1e-6);//fr1 is updated

//	try{
//		fr3.write_h5(h5file);
//	}catch(H5::FileIException & ex){
//		cout << ex.getDetailMsg() << endl;
//	}
//	BOOST_CHECK_EXCEPTION(fr3.write_h5(h5file);, H5::FileIException,
//						[](const H5::FileIException & ex) {return ex.getDetailMsg().find("H5Fcreate failed") != string::npos;});

}

BOOST_AUTO_TEST_CASE(set_range)
{
	MemCytoFrame fr1 = *fr.copy();

	string channel = fr1.get_channels()[1];
	pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> p = make_pair(-1.1, 1.1);
	fr1.set_range(channel, ColType::channel, p);
	auto p1 = fr1.get_range(channel, ColType::channel, RangeType::instrument);
	BOOST_CHECK_EQUAL(p1.first, p.first);
	BOOST_CHECK_EQUAL(p1.second, p.second);
	int idx = fr1.get_col_idx(channel, ColType::channel);
	string key = "flowCore_$P" + to_string(idx + 1) + "Rmax";
	BOOST_CHECK_EQUAL(fr1.get_keyword(key), boost::lexical_cast<string>(p.second));
}
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

	MemCytoFrame fr1 = *fr.copy();
	string oldname = channels[2];
	//query by upper case
	BOOST_CHECK_GT(fr1.get_col_idx(boost::to_lower_copy(oldname), ColType::channel), 0);

	BOOST_CHECK_EXCEPTION(fr1.set_channel(oldname, channels[1]), domain_error,
			[](const exception & ex) {return string(ex.what()).find("already exists") != string::npos;});
	string newname = "test";
	fr1.set_channel(oldname, newname);
	string key = fr1.get_keyword("$P3N");
	BOOST_CHECK_EQUAL(fr1.get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);

	fr1 = *fr.copy_realized({}, {});
	newname = "test1";
	fr1.set_channel(oldname, newname);
	key = fr1.get_keyword("$P3N");
	BOOST_CHECK_EQUAL(fr1.get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);

	//h5
	CytoFramePtr fr2 = fr_h5->copy();
	fr2->set_channel(oldname, newname);
	key = fr2->get_keyword("$P3N");
	BOOST_CHECK_EQUAL(fr2->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);

	string tmp = generate_unique_filename(fs::temp_directory_path(), "", ".h5");
	fr2 = fr_h5->copy_realized({1,2,3}, {1,2}, tmp);//channel order may change
	fr2->set_channel(oldname, newname);
	key = fr2->get_keyword("$P3N");
	BOOST_CHECK_GT(fr2->get_col_idx(newname, ColType::channel), 0);
	BOOST_CHECK_EQUAL(key, newname);

	//h5 is not synced by setter immediately
	H5CytoFrame fr3(tmp);
	BOOST_CHECK_EQUAL(fr3.get_col_idx(newname, ColType::channel), -1);
	BOOST_CHECK_EQUAL(fr3.get_keyword("$P3N"), oldname);
	//cached meta data is NoT flushed to h5 even when the object is destroyed
	fr2.reset();
	fr3 = H5CytoFrame(tmp);
	BOOST_CHECK_GT(fr3.get_col_idx(oldname, ColType::channel), 0);
	BOOST_CHECK_EQUAL(fr3.get_keyword("$P3N"), oldname);

}
BOOST_AUTO_TEST_CASE(shallow_copy)
{
	CytoFramePtr fr_orig = fr_h5->copy();//create a safe copy to test with by deep copying
	//perform shallow copy
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
	BOOST_CHECK_EQUAL(fr_orig->get_h5_file_path(), fr1.get_h5_file_path());
	//update data
	EVENT_DATA_VEC dat = fr1.get_data();
	float newval = 100;
	dat[100] = newval;
	fr1.set_data(dat);
	BOOST_CHECK_CLOSE(fr1.get_data()[100], newval, 1e-6);//fr1 is updated
	BOOST_CHECK_CLOSE(fr_orig->get_data()[100], newval, 1e-6);//original cp is also updated

	//delete fr_orig to ensure its copy still has valid access to h5 data
	fr_orig.reset();
	BOOST_CHECK_CLOSE(fr1.get_data()[100], newval, 1e-6);
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
	float newval = 100;
	dat[100] = newval;
	fr1->set_data(dat);
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);
	BOOST_CHECK_NE(fr_h5->get_h5_file_path(), fr1->get_h5_file_path());
	//original cp is not affected by any change in meta and event data
	BOOST_CHECK_EQUAL(fr_h5->get_channels()[2], oldname);
	BOOST_CHECK_EQUAL(fr_h5->get_keyword("$P3N"), oldname);

	BOOST_CHECK_CLOSE(fr1->get_data()[100], newval, 1e-6);
	BOOST_CHECK_CLOSE(fr_h5->get_data()[100], old_val, 1e-6);



}
BOOST_AUTO_TEST_CASE(flush_meta)
{
	//deep cp
	CytoFramePtr fr1 = fr_h5->copy();
	string h5file = fr1->get_h5_file_path();
	//update meta data
	string oldname = fr1->get_channels()[2];
	string newname = "test";
	fr1->set_channel(oldname, newname);
	//discard change
	fr1->load_meta();
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], oldname);
	BOOST_CHECK_EQUAL(fr1->get_keyword("$P3N"), oldname);

	fr1->set_channel(oldname, newname);
	//flush the change
	fr1->flush_meta();
	fr1->load_meta();
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(fr1->get_keyword("$P3N"), newname);
	//change it back and see if the destructor does the flushing
	fr1->set_channel(newname, oldname);
	fr1.reset();
	//load it back and see the change has NOT taken effect
	fr1.reset(new H5CytoFrame(h5file));
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(fr1->get_keyword("$P3N"), newname);

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
