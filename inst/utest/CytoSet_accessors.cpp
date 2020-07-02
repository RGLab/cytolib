#include <boost/test/unit_test.hpp>
#include <cytolib/GatingSet.hpp>

#include "fixture.hpp"
using namespace cytolib;
struct CSFixture {
	CSFixture() {

		file_paths = {"../flowWorkspace/wsTestSuite/curlyQuad/example1/A1001.001.fcs", "../flowWorkspace/wsTestSuite/curlyQuad/example1/A2002.001.fcs"};
		cs = GatingSet(file_paths, config, fmt);

	};

	~CSFixture() {

	};
	GatingSet cs;
	FileFormat fmt = FileFormat::TILE;
	vector<string> file_paths;
	FCS_READ_PARAM config;
};

BOOST_FIXTURE_TEST_SUITE(cs_test,CSFixture)
BOOST_AUTO_TEST_CASE(test) {
//	class A{
//	public:
//		int & r_;
//		A(int & r):r_(r){};
//	};
//	int * a = new int(3);
//	A b(*a);
//	b.r_ = 4;
//	cout << *a << endl;
//	delete a;
//	b.r_ = 5;
//	cout << b.r_ << endl;

}
BOOST_AUTO_TEST_CASE(constructor) {
	file_paths[1] = file_paths[0];

	BOOST_CHECK_EXCEPTION(GatingSet(file_paths, config, fmt, "");, domain_error,
			[](const exception & ex) {return string(ex.what()).find("already exists") != string::npos;});
}
BOOST_AUTO_TEST_CASE(copy) {
	GatingSet cs1 = cs.copy();
	vector<string> samples = cs.get_sample_uids();
	CytoFrameView fv = cs.get_cytoframe_view_ref(samples[1]);
	CytoFrameView fv1 = cs1.get_cytoframe_view_ref(samples[1]);
	BOOST_CHECK_CLOSE(fv.get_data()[1], fv1.get_data()[1], 1e-6);
	BOOST_CHECK_CLOSE(fv.get_data()[7e4], fv1.get_data()[7e4], 1e-6);
}
BOOST_AUTO_TEST_CASE(subset_by_cols) {
	vector<string> channels = cs.get_channels();
	vector<string> markers = cs.get_markers();
	vector<string> sub_channels = { channels[3], channels[1] };
	vector<string> sub_markers = { markers[3], markers[1] };
	GatingSet cs_new = cs.copy();
	cs_new.cols_(sub_channels, ColType::channel);

	BOOST_CHECK_EQUAL(cs_new.n_cols(), sub_channels.size());

	vector<string> markers_new = cs_new.get_markers();
	BOOST_CHECK_EQUAL_COLLECTIONS(markers_new.begin(), markers_new.end(),
			sub_markers.begin(), sub_markers.end());

//reorder
	vector<string> chnls_sort = channels;
	sort(chnls_sort.begin(), chnls_sort.end());
	cs_new = cs.copy();
	cs_new.cols_(chnls_sort, ColType::channel);
	vector<string> chnls_new = cs_new.get_channels();
	BOOST_CHECK_EQUAL_COLLECTIONS(chnls_new.begin(), chnls_new.end(),
			chnls_sort.begin(), chnls_sort.end());
	//restore
	cs_new.cols_(channels, ColType::channel);
	chnls_new = cs_new.get_channels();
	BOOST_CHECK_EQUAL_COLLECTIONS(chnls_new.begin(), chnls_new.end(),
			channels.begin(), channels.end());


}
BOOST_AUTO_TEST_CASE(subset_by_sample) {
	//check get_sample_uids
	vector<string> samples = cs.get_sample_uids();
	sort(samples.begin(), samples.end());
	vector<string> filenames(file_paths.size());
	transform(file_paths.begin(), file_paths.end(), filenames.begin(),
			[](const string &i) {return path_base_name(i);});
	BOOST_CHECK_EQUAL_COLLECTIONS(samples.begin(), samples.end(),
			filenames.begin(), filenames.end());

	//subset
	vector<string> select = { samples[0] };
	GatingSet cs_new = cs.sub_samples(select);
	BOOST_CHECK_EQUAL(cs_new.size(), 1);
	vector<string> new_samples = cs_new.get_sample_uids();
	sort(new_samples.begin(), new_samples.end());
	sort(select.begin(), select.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(select.begin(), select.end(),
			new_samples.begin(), new_samples.end());
}
BOOST_AUTO_TEST_CASE(cytoframe) {
	vector<string> samples = cs.get_sample_uids();
	sort(samples.begin(), samples.end());
	CytoFrameView fv = cs.get_cytoframe_view_ref(samples[1]);
	CytoFramePtr fr = fv.get_cytoframe_ptr();

	//get_cytoframe from subseted cs
	vector<string>
	select = {samples[0]};
	GatingSet cs_new = cs.sub_samples(select);
	BOOST_CHECK_EXCEPTION(cs_new.get_cytoframe_view_ref(samples[1]);;, domain_error,
			[](const exception & ex) {return string(ex.what()).find("not found") != string::npos;});

	//add frame
	BOOST_CHECK_EXCEPTION(cs_new.add_cytoframe_view(samples[0], CytoFrameView(fr)), domain_error,
			[](const exception & ex) {return string(ex.what()).find("already exists") != string::npos;});
	BOOST_CHECK_EXCEPTION(cs_new.update_cytoframe_view(samples[1], CytoFrameView(fr)), domain_error,
			[](const exception & ex) {return string(ex.what()).find("doesn't exists") != string::npos;});

	cs_new.add_cytoframe_view(samples[1], CytoFrameView(fr));
	BOOST_CHECK_EQUAL(cs_new.size(), 2);
	vector<string> new_samples = cs_new.get_sample_uids();
	sort(new_samples.begin(), new_samples.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(samples.begin(), samples.end(),
			new_samples.begin(), new_samples.end());

}
BOOST_AUTO_TEST_SUITE_END()
