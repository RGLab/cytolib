#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cytolib/GatingSet.hpp>
#include <experimental/filesystem>
#include <regex>

#include "fixture.hpp"
using namespace cytolib;
struct GSFixture {
	GSFixture() {
		path = "../flowWorkspace/output/NHLBI/gs/gs";
		gs = GatingSet(path,false,true,{},false);

	};

	~GSFixture() {

	};
	GatingSet gs;
	string path;
};

BOOST_FIXTURE_TEST_SUITE(GatingSet_test,GSFixture)
#ifdef HAVE_TILEDB

BOOST_AUTO_TEST_CASE(s3_gs)
{


	CytoCtx ctx(string(std::getenv("AWS_ACCESS_KEY_ID"))
				, string(std::getenv("AWS_SECRET_ACCESS_KEY"))
				, "us-west-1"
				);
	auto remote = "s3://mike-h5/test";

		//convert h5 to tile
	auto gs1 = gs.copy();
	auto gh = gs1.begin()->second;
	auto cfv = gh->get_cytoframe_view();
	string tmp = generate_unique_dir(fs::temp_directory_path().c_str(), "") + ".tile";

	cfv.write_to_disk(tmp, FileFormat::TILE, ctx);
	gh->set_cytoframe_view(CytoFramePtr(new TileCytoFrame(tmp)));
	gs1.serialize_pb(remote, CytoFileOption::copy, false, ctx);
//	tiledb::VFS vfs(*ctx);
//	for(auto p : vfs.ls(remote))
//	{
//		cout << p << endl;
//	}
	auto gs2 = GatingSet(remote,false,true,{},false, ctx);

	auto cf = gs2.begin()->second->get_cytoframe_view();
	auto ch = cf.get_channels();
	BOOST_CHECK_EQUAL(ch.size(), 12);
//

}
#endif
BOOST_AUTO_TEST_CASE(remove_node)
{
	auto gs1 = gs.copy();
	auto gh = gs1.begin()->second;
	gh->removeNode(0);
	BOOST_CHECK_EQUAL(gh->getTree().m_vertices.size(), 23);

	 gs1 = gs.copy();
	 gh = gs1.begin()->second;
	gh->removeNode("root");
	BOOST_CHECK_EQUAL(gh->getTree().m_vertices.size(), 0);

}
BOOST_AUTO_TEST_CASE(quadgate) {

	GatingSet gs1({"../flowWorkspace/output/s5a01.fcs"}, FCS_READ_PARAM());
	auto gh = gs1.begin()->second;
	//test regular rectgate
	shared_ptr<rectGate> g(new rectGate());
	paramPoly p;
	p.setName({"FSC-H", "SSC-H"});
	p.setVertices({coordinate(300,0), coordinate(500,400)});
	g->setParam(p);
	auto id = gh->addGate(g, 0, "test");
	auto cf = MemCytoFrame(*(gh->get_cytoframe_view().get_cytoframe_ptr()));
	gh->gating(cf, 0, true, true);
	BOOST_CHECK_EQUAL(gh->getNodeProperty(gh->getNodeID("test")).getCounts(), 649);
	//quadgate
//	gh->removeNode(id);
	p.setVertices({coordinate(500,600)});
	gh->addGate(gatePtr(new quadGate(p, "123", Q1)), id, "A");
	gh->addGate(gatePtr(new quadGate(p, "123", Q2)), id, "B");
	gh->addGate(gatePtr(new quadGate(p, "123", Q3)), id, "C");
	gh->addGate(gatePtr(new quadGate(p, "123", Q4)), id, "D");
	gh->gating(cf, 0, true, true);
	BOOST_CHECK_EQUAL(gh->getNodeProperty(gh->getNodeID("test")).getCounts()
					, gh->getNodeProperty(gh->getNodeID("A")).getCounts()+
					gh->getNodeProperty(gh->getNodeID("B")).getCounts()+
					gh->getNodeProperty(gh->getNodeID("C")).getCounts()+
					gh->getNodeProperty(gh->getNodeID("D")).getCounts());
	//test pb
	string tmp = generate_unique_dir(fs::temp_directory_path().c_str(), "gs");
	gs1.serialize_pb(tmp, CytoFileOption::copy);
	gs1 = GatingSet(tmp);
	gh = gs1.begin()->second;
	gh->gating(cf, 0, true, true);
	BOOST_CHECK_EQUAL(gh->getNodeProperty(gh->getNodeID("test")).getCounts()
					, gh->getNodeProperty(gh->getNodeID("A")).getCounts()+
					gh->getNodeProperty(gh->getNodeID("B")).getCounts()+
					gh->getNodeProperty(gh->getNodeID("C")).getCounts()+
					gh->getNodeProperty(gh->getNodeID("D")).getCounts());

}
BOOST_AUTO_TEST_CASE(serialize) {
	GatingSet gs1 = gs.copy();
	/*
	 * make changes to meta data
	 */
	auto cf = gs1.begin()->second->get_cytoframe_view_ref();
	string oldc = cf.get_channels()[0];
	string newc = "new_channel";
	cf.set_channel(oldc, newc);
	string kn = cf.get_keywords().begin()->first;
	string kv = "new_key";
	cf.set_keyword(kn, kv);
	string pdn = "pdn";
	string pdv = "pdv";
	cf.set_pheno_data(pdn, pdv);
	//archive
//	string tmp = "/tmp/gsEGnOlP";
	string tmp = generate_unique_dir(fs::temp_directory_path().c_str(), "gs");
	gs1.serialize_pb(tmp, CytoFileOption::copy);
	//load it back
	GatingSet gs2(tmp,false,true,{},true);
	auto cf2 = gs2.begin()->second->get_cytoframe_view_ref();
	BOOST_CHECK_EQUAL(cf2.get_readonly(), true);


	//save gs and see if readonly flag changes
	tmp = generate_unique_dir(fs::temp_directory_path().c_str(), "gs");
	gs2.serialize_pb(tmp, CytoFileOption::copy);
	BOOST_CHECK_EQUAL(cf2.get_readonly(), true);


	BOOST_CHECK_EQUAL(cf2.get_channels()[0], newc);
	BOOST_CHECK_EQUAL(cf2.get_keyword(kn), kv);
	BOOST_CHECK_EQUAL(cf2.get_pheno_data().find(pdn)->second, pdv);

}
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
	vector<string> samples = gs.get_sample_uids();
	const GatingSet & cs = gs.get_cytoset();
	CytoFrameView fv = cs.get_cytoframe_view(samples[0]);
	BOOST_CHECK_EQUAL(fv.get_readonly(), true);

	GatingSet gs1 = gs.copy();
	//flag remains unchanged for original gs
	BOOST_CHECK_EQUAL(fv.get_readonly(), true);
	const GatingSet & cs1 = gs1.get_cytoset();
	CytoFrameView fv1 = cs1.get_cytoframe_view(samples[0]);
	//new cp is writable by default
	BOOST_CHECK_EQUAL(fv1.get_readonly(), false);
	fv1.set_readonly(true);
	BOOST_CHECK_EQUAL(fv1.get_readonly(), true);
	fv1.set_channel("SSC-A", "test");
	BOOST_CHECK_EXCEPTION(fv1.flush_meta(), domain_error,
				[](const domain_error & ex) {return string(ex.what()).find("read-only") != string::npos;});

	BOOST_CHECK_NE(fv.get_uri(), fv1.get_uri());
	BOOST_CHECK_CLOSE(fv.get_data()[1], fv1.get_data()[1], 1e-6);
	BOOST_CHECK_CLOSE(fv.get_data()[7e4], fv1.get_data()[7e4], 1e-6);
}
BOOST_AUTO_TEST_CASE(legacy_gs) {
	auto cs = GatingSet();
	auto sn = "CytoTrol_CytoTrol_1.fcs";
	cs.add_cytoframe_view(sn, CytoFrameView(CytoFramePtr(new MemCytoFrame())));//add dummy frame
	GatingSet gs1 = GatingSet("../flowWorkspaceData/inst/extdata/legacy_gs/gs_manual/jzgmkCmwZR.pb", cs);
	vector<string> samples = gs1.get_sample_uids();
	BOOST_CHECK_EQUAL(samples[0], sn);

	GatingHierarchyPtr gh = gs1.getGatingHierarchy(samples[0]);
	VertexID_vec vid = gh->getVertices();
	BOOST_CHECK_EQUAL(vid.size(), 24);
	BOOST_CHECK_EQUAL(gh->getNodePath(vid[16]), "/not debris/singlets/CD3+/CD8/38+ DR-");
//	BOOST_CHECK_EQUAL(gh->get_cytoframe_view().get_keyword("$BEGINDATA"), "3264");

	//save legacy to new format
	string tmp = generate_unique_dir(fs::temp_directory_path().c_str(), "gs");
	bool is_skip_data = true;
	gs1.serialize_pb(tmp, CytoFileOption::skip, is_skip_data);
	gs1 = GatingSet(tmp, is_skip_data);
	gh = gs1.getGatingHierarchy(samples[0]);
	vid = gh->getVertices();
	BOOST_CHECK_EQUAL(vid.size(), 24);
	BOOST_CHECK_EQUAL(gh->getNodePath(vid[16]), "/not debris/singlets/CD3+/CD8/38+ DR-");

	//save new to new format
	tmp = generate_unique_dir(fs::temp_directory_path().c_str(), "gs");
	gs.serialize_pb(tmp, CytoFileOption::copy);
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
