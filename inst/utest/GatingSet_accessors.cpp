#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cytolib/GatingSet.hpp>

#include "fixture.hpp"
using namespace cytolib;
struct GSFixture {
	GSFixture() {

		path = "../flowWorkspace/output/NHLBI/gs/gs";
		gs = GatingSet(path);

	};

	~GSFixture() {

	};
	GatingSet gs;
	string path;
};

BOOST_FIXTURE_TEST_SUITE(GatingSet_test,GSFixture)
BOOST_AUTO_TEST_CASE(get_cytoset) {
	GatingSet cs = gs.get_cytoset();
	BOOST_CHECK_EQUAL(cs.get_sample_uids().size(), gs.size());
	GatingSet cs1(cs);
	BOOST_CHECK_EQUAL(cs1.get_sample_uids().size(), gs.size());

}
BOOST_AUTO_TEST_CASE(subset_cs_by_node) {

	GatingSet gs1 = gs.sub_samples(gs.get_sample_uids());
//	cout << gs1.get_uid() << endl;
	BOOST_CHECK_NE(gs.get_uid(), gs1.get_uid());
	GatingSet cs = gs.get_cytoset("CD8");
	GatingHierarchy frv = *cs.begin()->second;
	BOOST_CHECK_EQUAL(frv.n_rows(), 14564);

	cs = gs.get_cytoset();
	BOOST_CHECK_EQUAL(cs.get_channels().size(), 12);
	frv = *cs.begin()->second;
	BOOST_CHECK_EQUAL(frv.n_rows(), 119531);
}
BOOST_AUTO_TEST_CASE(copy) {
	GatingSet gs1 = gs.copy();
	vector<string> samples = gs.get_sample_uids();
	const GatingSet & cs = gs.get_cytoset();
	const GatingSet & cs1 = gs1.get_cytoset();
	CytoFrameView fv = cs.get_cytoframe_view(samples[0]);
	CytoFrameView fv1 = cs1.get_cytoframe_view(samples[0]);
	BOOST_CHECK_NE(fv.get_h5_file_path(), fv1.get_h5_file_path());
	BOOST_CHECK_CLOSE(fv.get_data()[1], fv1.get_data()[1], 1e-6);
	BOOST_CHECK_CLOSE(fv.get_data()[7e4], fv1.get_data()[7e4], 1e-6);
}
BOOST_AUTO_TEST_CASE(legacy_gs) {
	GatingSet gs1 = GatingSet("../flowWorkspaceData/inst/extdata/legacy_gs/gs_manual/jzgmkCmwZR.pb",GatingSet());
	vector<string> samples = gs1.get_sample_uids();
	BOOST_CHECK_EQUAL(samples[0], "CytoTrol_CytoTrol_1.fcs");

	GatingHierarchyPtr gh = gs1.getGatingHierarchy(samples[0]);
	VertexID_vec vid = gh->getVertices();
	BOOST_CHECK_EQUAL(vid.size(), 24);
	BOOST_CHECK_EQUAL(gh->getNodePath(vid[16]), "/not debris/singlets/CD3+/CD8/38+ DR-");
//	BOOST_CHECK_EQUAL(gh->get_cytoframe_view().get_keyword("$BEGINDATA"), "3264");

	//save legacy to new format
	string tmp = generate_temp_filename();
	bool is_skip_data = true;
	gs1.serialize_pb(tmp, H5Option::skip, is_skip_data);
	gs1 = GatingSet(tmp, is_skip_data);
	gh = gs1.getGatingHierarchy(samples[0]);
	vid = gh->getVertices();
	BOOST_CHECK_EQUAL(vid.size(), 24);
	BOOST_CHECK_EQUAL(gh->getNodePath(vid[16]), "/not debris/singlets/CD3+/CD8/38+ DR-");

	//save new to new format
	tmp = generate_temp_filename();
	gs.serialize_pb(tmp, H5Option::copy);
	gs1 = GatingSet(tmp);
	gh = gs1.getGatingHierarchy(samples[0]);
	vid = gh->getVertices();
	BOOST_CHECK_EQUAL(vid.size(), 24);
	BOOST_CHECK_EQUAL(gh->getNodePath(vid[16]), "/not debris/singlets/CD3+/CD8/38+ DR-");


}
BOOST_AUTO_TEST_CASE(template_constructor) {
	GatingHierarchy gh=*gs.getGatingHierarchy(gs.get_sample_uids()[0]);

	/*
	 * used gh as the template to clone multiple ghs in the new gs
	 */
	vector<pair<string, string>> id_vs_path;
	id_vs_path.push_back(make_pair("aa", "../flowWorkspaceData/inst/extdata/CytoTrol_CytoTrol_2.fcs"));
	GatingSet cs(id_vs_path);
	GatingSet gs1 = GatingSet(gh, cs);

	BOOST_CHECK_EQUAL(gs1.get_sample_uids()[0], "aa");
	GatingHierarchyPtr gh1 = gs1.getGatingHierarchy("aa");
	VertexID_vec vid = gh1->getVertices();
	BOOST_CHECK_EQUAL(vid.size(), 24);
	BOOST_CHECK_EQUAL(gh1->getNodePath(vid[16]), "/not debris/singlets/CD3+/CD8/38+ DR-");
}
//BOOST_AUTO_TEST_CASE(subset_by_sample) {
//	//check get_sample_uids
//	vector<string> samples = gs.get_sample_uids();
//	sort(samples.begin(), samples.end());
//	vector<string> filenames(file_paths.size());
//	transform(file_paths.begin(), file_paths.end(), filenames.begin(),
//			[](const string &i) {return path_base_name(i);});
//	BOOST_CHECK_EQUAL_COLLECTIONS(samples.begin(), samples.end(),
//			filenames.begin(), filenames.end());
//
//	//subset
//	vector<string> select = { samples[0] };
//	GatingSet cs_new = cs.sub_samples(select);
//	BOOST_CHECK_EQUAL(cs_new.size(), 1);
//	vector<string> new_samples = cs_new.get_sample_uids();
//	sort(new_samples.begin(), new_samples.end());
//	sort(select.begin(), select.end());
//	BOOST_CHECK_EQUAL_COLLECTIONS(select.begin(), select.end(),
//			new_samples.begin(), new_samples.end());
//}
BOOST_AUTO_TEST_SUITE_END()
