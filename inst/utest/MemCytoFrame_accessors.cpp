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

	};

	~CFFixture(){

	};
	MemCytoFrame fr;
   string file_path;
   FCS_READ_PARAM config;
};

BOOST_FIXTURE_TEST_SUITE(MemCytoFrame_test,CFFixture)

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

	CytoFramePtr fr_copy = cr_new.copy();
	BOOST_CHECK_EQUAL(fr_copy->n_rows(), 3);

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
BOOST_AUTO_TEST_SUITE_END()
