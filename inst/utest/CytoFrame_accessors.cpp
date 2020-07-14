#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/lexical_cast.hpp>
#include <cytolib/GatingSet.hpp>
#include <cytolib/TileCytoFrame.hpp>
#include <cytolib/H5CytoFrame.hpp>
#include <cytolib/MemCytoFrame.hpp>

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
		string tmp = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
//		cout << tmp << endl;
		fr.write_to_disk(tmp, file_format);
		cf_disk = load_cytoframe(tmp);

	};

	~CFFixture(){

	};
	MemCytoFrame fr;
	FileFormat file_format = FileFormat::TILE;
//	FileFormat file_format = FileFormat::H5;
//	unique_ptr<H5CytoFrame> cf_disk;
	CytoFramePtr cf_disk;
   string file_path;
   FCS_READ_PARAM config;
};

BOOST_FIXTURE_TEST_SUITE(CytoFrame_test,CFFixture)
#ifdef HAVE_TILEDB
BOOST_AUTO_TEST_CASE(TileCytoFrameconstructor)
{
//	TileCytoFrame cf = *(dynamic_cast<TileCytoFrame*>(cf_disk.get()));
//	cf_disk.reset();
//	auto dat = cf.get_data();
//	cout << dat[10] << endl;
	string tmp = generate_unique_filename(fs::temp_directory_path().string(), "", ".tile");
//	auto path = file_path;
//	auto cf = TileCytoFrame(path, config, tmp);
	auto path = "/tmp/Rtmprb2X3k/file1be46cf4a372.tile";
	auto cf = load_cytoframe(path, true);
	auto dat = cf->get_data();
}
BOOST_AUTO_TEST_CASE(tile_write_block_test)
{
	/*
	 * create and wirte array
	 */
	tiledb::Context ctx;
	tiledb::Domain domain(ctx);
	domain.add_dimension(tiledb::Dimension::create<int>(ctx, "cell", {1, 4}, 2));

	tiledb::ArraySchema schema(ctx, TILEDB_DENSE);
	schema.set_domain(domain);
	schema.add_attribute(tiledb::Attribute::create<int>(ctx, "a1"));
	schema.set_tile_order(TILEDB_COL_MAJOR).set_cell_order(TILEDB_COL_MAJOR);

	auto array_uri = "/tmp/test.tile";
	tiledb::VFS vfs(ctx);
	if(vfs.is_dir(array_uri))
		vfs.remove_dir(array_uri);
	tiledb::Array::create(array_uri, schema);
	tiledb::Array array(ctx, array_uri, TILEDB_WRITE);
	tiledb::Query query(ctx, array);
	query.set_layout(TILEDB_GLOBAL_ORDER);
	query.set_subarray({1,4});
	vector<int> buf = {1, 2, 3, 4};

	query.set_buffer("a1", buf);
	query.submit();
	query.finalize();
	string v1 = "v1";
	array.put_metadata("k1", TILEDB_CHAR, v1.size(), v1.c_str());
	v1 = "";
	array.put_metadata("k2", TILEDB_CHAR, v1.size(), v1.c_str());

	/*
	 * open the array and read it
	 */
	array.close();
	tiledb::Array array1(ctx, array_uri, TILEDB_READ);
	vector<int> buf1(4);
	tiledb::Query query1(ctx, array1);
	query1.set_subarray({1,4});
	query1.set_layout(TILEDB_GLOBAL_ORDER);
	query1.set_buffer("a1", buf1);
	query1.submit();
	query1.finalize();

	uint64_t nkw = array1.metadata_num();



	uint32_t v_num;
	tiledb_datatype_t v_type;
	for (uint64_t i = 0; i < nkw; ++i) {
		string key,val;
		const void* v;
		array1.get_metadata_from_index(i, &key, &v_type, &v_num, &v);
		if(v)
		{
			val = string(static_cast<const char *>(v));
			val.resize(v_num);
		}
		cout << key << ":" << val << endl;
	}
//	array.close();
	/*
	 * open it with another array object
	 */
	tiledb::Array array2(ctx, array_uri, TILEDB_WRITE);
	vector<int> buf2= {5,6,7,8};
	tiledb::Query query2(ctx, array2);
	query2.set_subarray({1,4});
	query2.set_layout(TILEDB_GLOBAL_ORDER);
	query2.set_buffer("a1", buf2);
	query2.submit();
//	query2.submit_async([]() { std::cout << "Callback: Query completed\n"; });
	query2.finalize();

	array.close();
	array.open(TILEDB_READ);
	tiledb::Query query3(ctx, array);
	query3.set_subarray({1,4});
	query3.set_layout(TILEDB_GLOBAL_ORDER);
	query3.set_buffer("a1", buf1);
	query3.submit();
	query3.finalize();

	for(auto i : buf1)
		cout << i << endl;
}
BOOST_AUTO_TEST_CASE(tile)
{
//	auto uri = "s3://mike-h5/file24ad1cad7ca7.tile";
	auto uri = "/tmp/test.tile";
	CytoCtx ctx(string(std::getenv("AWS_ACCESS_KEY_ID"))
					, string(std::getenv("AWS_SECRET_ACCESS_KEY"))
					, "us-west-1"
					);
	CytoVFS vfs(ctx);

	if(vfs.is_dir(uri))
		vfs.remove_dir(uri);
	auto h5 = "/tmp/test.h5";
	fr.write_h5(h5);
	H5CytoFrame fr_h5(h5);
	fr.write_tile(uri, ctx);
	auto cf_tile = TileCytoFrame(uri, true, true, ctx);

	auto ch = cf_tile.get_channels();
//	BOOST_CHECK_EQUAL(ch.size(), 9);

	//read all
	double start = gettime();
	auto mat1 = fr_h5.get_data();
	double runtime = (gettime() - start);
	cout << "fr_h5.get_data(): " << runtime << endl;

	start = gettime();
	auto mat2 = cf_tile.get_data();
	runtime = (gettime() - start);
	cout << "cf_tile->get_data(): " << runtime << endl;

	for(auto i : {100,1000,6000})
		BOOST_CHECK_CLOSE(mat1.mem[i], mat2.mem[i], 1);

	//idx by row and col
	auto nrow = 200;
	uvec ridx(nrow);
	srand (1);//seed
	for(int i = 0; i < nrow; i++)
		ridx[i] = rand()%fr.n_rows();

	start = gettime();
	mat1 = fr_h5.get_data(ridx, {3,8});
	runtime = (gettime() - start);
	cout << "fr_h5.get_data(ridx,cidx): " << runtime << endl;

	start = gettime();
	mat2 = cf_tile.get_data(ridx, {3,8});
	runtime = (gettime() - start);
	cout << "cf_tile->get_data(ridx, cidx): " << runtime << endl;
//	for(auto i : {0,1,2,3})
//		BOOST_CHECK_CLOSE(mat1.mem[i], mat2.mem[i], 1);
//	//idx by row
//	mat1 = fr_h5.get_data({1,100}, false);
//	mat2 = cf_tile.get_data({1,100}, false);
//	for(auto i : {0,10,15})
//		BOOST_CHECK_CLOSE(mat1.mem[i], mat2.mem[i], 1);
//	//idx by col
//	mat1 = fr_h5.get_data({1,3}, true);
//	mat2 = cf_tile.get_data({1,3}, true);
//	for(auto i : {0,1000,5500})
//		BOOST_CHECK_CLOSE(mat1.mem[i], mat2.mem[i], 1);
}
BOOST_AUTO_TEST_CASE(s3)
{
//	H5RCytoFrame cf = H5RCytoFrame("https://mike-h5.s3.amazonaws.com/bcell.h5", true);
//	auto ch = cf.get_channels();
//	BOOST_CHECK_EQUAL(ch.size(), 10);


}
#endif
BOOST_AUTO_TEST_CASE(spillover)
{
	auto comp = fr.get_compensation();
	auto markers = fr.get_channels();
	BOOST_CHECK_EQUAL_COLLECTIONS(comp.marker.begin(), comp.marker.end(), markers.begin()+3, markers.end());
	auto txt = comp.to_string();
	auto comp1 = compensation(txt);
	BOOST_CHECK_EQUAL_COLLECTIONS(comp.marker.begin(), comp.marker.end(), comp1.marker.begin(), comp1.marker.end());
	for(unsigned i = 0; i < comp.spillOver.size(); i++)
		BOOST_CHECK_CLOSE(comp.spillOver[i], comp1.spillOver[i], 1);

}
BOOST_AUTO_TEST_CASE(profile_get_data)
{
	auto fr1 = MemCytoFrame("../flowWorkspace/wsTestSuite/profile_get_data.fcs", config);
	fr1.read_fcs();
	double start = gettime();
	auto dat = fr1.get_data();
	double runtime = (gettime() - start);
	cout << "get a copy: " << runtime << endl;

	start = gettime();
	auto &dat_ref = fr1.get_data_ref();
	runtime = (gettime() - start);
	cout << "get a reference: " << runtime << endl;
}
BOOST_AUTO_TEST_CASE(get_time_step)
{
	auto fr1 = MemCytoFrame("../flowWorkspace/output/s5a01.fcs", config);
	fr1.read_fcs();
	string ts = fr1.get_keyword("$BTIM");
	TM_ext btime = parse_time_with_fractional_seconds(ts);
	BOOST_CHECK_EQUAL(btime._time.tm_hour, 9);
	BOOST_CHECK_EQUAL(btime._time.tm_min, 51);
	BOOST_CHECK_EQUAL(btime._time.tm_sec, 34);
	BOOST_CHECK_CLOSE(fr1.get_time_step("Time"), 0.0557, 0.01);

}
BOOST_AUTO_TEST_CASE(h5_vs_mem)
{
	//check if h5 version is consistent with mem
	BOOST_CHECK_EQUAL(fr.get_params().size(), cf_disk->get_params().size());
	BOOST_CHECK_EQUAL(fr.get_params().begin()->max, cf_disk->get_params().begin()->max);

}
BOOST_AUTO_TEST_CASE(flags)
{
	if(file_format == FileFormat::H5)
	{
	//get a safe cp
	string h5file = cf_disk->copy()->get_uri();
	//load it by default readonly flag
	unique_ptr<H5CytoFrame> fr1(new H5CytoFrame(h5file));
	//update meta data
	string oldname = fr1->get_channels()[2];
	string newname = "test";
	fr1->set_channel(oldname, newname);
	BOOST_CHECK_EXCEPTION(fr1->flush_params(), domain_error,
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
	}
//	try{
//		fr3.write_to_disk(h5file);
//	}catch(H5::FileIException & ex){
//		cout << ex.getDetailMsg() << endl;
//	}
//	BOOST_CHECK_EXCEPTION(fr3.write_to_disk(h5file);, H5::FileIException,
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

	string tmp = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
	shared_ptr<CytoFrame> cf;;
	cr_new.write_to_disk(tmp, file_format);
	cf = load_cytoframe(tmp);

	BOOST_CHECK_EQUAL(cf->n_cols(), sub_channels.size());
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

	string tmp = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
	cr_new.write_to_disk(tmp);
	CytoFramePtr cf;
	cf = load_cytoframe(tmp);
	BOOST_CHECK_EQUAL(cf->n_rows(), 3);
	uvec idx = {};
	auto cf1 = cf->copy(idx, idx, "", false);
	BOOST_CHECK_EQUAL(cf1->n_rows(), 0);
	BOOST_CHECK_EQUAL(cf1->n_cols(), 0);
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

	fr1 = *fr.copy({}, true, "");
	newname = "test1";
	fr1.set_channel(oldname, newname);
	key = fr1.get_keyword("$P3N");
	BOOST_CHECK_EQUAL(fr1.get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);

	//h5
	CytoFramePtr fr2 = cf_disk->copy();
	fr2->set_channel(oldname, newname);
	key = fr2->get_keyword("$P3N");
	BOOST_CHECK_EQUAL(fr2->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);

	string tmp = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
	fr2 = cf_disk->copy({1,2,3}, {1,2}, tmp);//channel order may change
	fr2->set_channel(oldname, newname);
	key = fr2->get_keyword("$P3N");
	BOOST_CHECK_GT(fr2->get_col_idx(newname, ColType::channel), 0);
	BOOST_CHECK_EQUAL(key, newname);

	//h5 is not synced by setter immediately
	shared_ptr<CytoFrame> fr3;
	fr3 = load_cytoframe(tmp);
	BOOST_CHECK_EQUAL(fr3->get_col_idx(newname, ColType::channel), -1);
	BOOST_CHECK_EQUAL(fr3->get_keyword("$P3N"), oldname);
	//cached meta data is NoT flushed to h5 even when the object is destroyed
	fr2.reset();
	fr3 = load_cytoframe(tmp);
	BOOST_CHECK_GT(fr3->get_col_idx(oldname, ColType::channel), 0);
	BOOST_CHECK_EQUAL(fr3->get_keyword("$P3N"), oldname);

}

BOOST_AUTO_TEST_CASE(append_columns)
{
  MemCytoFrame fr1 = *fr.copy();
  uvec copy_idx;
  copy_idx << 6 << 7 << 8;
  EVENT_DATA_VEC new_cols = fr1.get_data(copy_idx, true);
  
  // Test all guards
  vector<string> new_names;
  // Empty colname vector
  BOOST_CHECK_EXCEPTION(fr1.append_columns(new_names, new_cols);;, domain_error,
                        [](const exception & ex) {return string(ex.what()).find("Must have equal (nonzero)") != string::npos;});
  // Length mismatch (only 2 names)
  new_names = {"new_channel_1", "new_channel_2"};
  BOOST_CHECK_EXCEPTION(fr1.append_columns(new_names, new_cols);;, domain_error,
                        [](const exception & ex) {return string(ex.what()).find("Must have equal (nonzero)") != string::npos;});
  // Duplicates an existing channel
  new_names = {"new_channel_1", "new_channel_2", "PE-A"};
  BOOST_CHECK_EXCEPTION(fr1.append_columns(new_names, new_cols);;, domain_error,
                        [](const exception & ex) {return string(ex.what()).find("already contains") != string::npos;});
  // Duplicates within the new channels
  new_names = {"new_channel_1", "new_channel_2", "new_channel_2"};
  BOOST_CHECK_EXCEPTION(fr1.append_columns(new_names, new_cols);;, domain_error,
                        [](const exception & ex) {return string(ex.what()).find("Duplicate new channel names detected") != string::npos;});
  // Empty channel name
  new_names = {"new_channel_1", "", "new_channel_2"};
  BOOST_CHECK_EXCEPTION(fr1.append_columns(new_names, new_cols);;, domain_error,
                        [](const exception & ex) {return string(ex.what()).find("must be non-empty strings") != string::npos;});
  
  new_names = {"new_channel_1", "new_channel_2", "new_channel_3"};
  EVENT_DATA_VEC shortened_columns = new_cols.rows(0,42);
  // Columns too short
  BOOST_CHECK_EXCEPTION(fr1.append_columns(new_names, shortened_columns);;, domain_error,
                        [](const exception & ex) {return string(ex.what()).find("same number of rows") != string::npos;});
  
  //  This should succeed
  fr1.append_columns(new_names, new_cols);
  // Check dims
  BOOST_CHECK_EQUAL(fr.n_cols() + new_names.size(), fr1.n_cols());
  BOOST_CHECK_EQUAL(fr.n_rows(), fr1.n_rows());
  vector<cytoParam> old_params = fr.get_params();
  vector<cytoParam> new_params = fr1.get_params();
  BOOST_CHECK_EQUAL(old_params.size() + new_names.size(), new_params.size());
  // Spot check params
  BOOST_CHECK_EQUAL(new_cols.col(1).min(), new_params[fr1.n_cols()-2].min);
  BOOST_CHECK_EQUAL(new_cols.col(2).max(), new_params[fr1.n_cols()-1].max);
  BOOST_CHECK_EQUAL(new_params[fr1.n_cols()-3].PnG, 1);
  BOOST_CHECK_EQUAL(new_params[fr1.n_cols()-2].PnB, 32);
  BOOST_CHECK_EQUAL(new_params[fr1.n_cols()-1].PnE[0], 0.0);
  BOOST_CHECK_EQUAL(new_params[fr1.n_cols()-1].PnE[1], 0.0);
  // Spot check keywords
  BOOST_CHECK_EQUAL(fr1.get_keyword("$P" + to_string(fr1.n_cols()) + "N"), new_names[2]);
  BOOST_CHECK_EQUAL(fr1.get_keyword("$P" + to_string(fr1.n_cols()-1) + "B"), "32");
  BOOST_CHECK_EQUAL(fr1.get_keyword("$P" + to_string(fr1.n_cols()-2) + "E"), "0,0");
  BOOST_CHECK_EQUAL(fr1.get_keyword("$P" + to_string(fr1.n_cols()) + "G"), "");
  BOOST_CHECK_EQUAL(fr1.get_keyword("$P" + to_string(fr1.n_cols()-1) + "R"), to_string(new_params[fr1.n_cols()-2].max + 1));
}

BOOST_AUTO_TEST_CASE(shallow_copy)
{
	CytoFramePtr fr_orig = cf_disk->copy();//create a safe copy to test with by deep copying

	//perform shallow copy
	shared_ptr<CytoFrame> fr1;
	if(file_format == FileFormat::H5)
		fr1.reset(new H5CytoFrame(*dynamic_cast<H5CytoFrame*>(fr_orig.get())));
	else
#ifdef HAVE_TILEDB
		fr1.reset(new TileCytoFrame(*dynamic_cast<TileCytoFrame*>(fr_orig.get())));
#else
	throw(domain_error("unsupported format: " + fmt_to_str(file_format)));

#endif


	//update meta data
	string oldname = fr1->get_channels()[2];
	string newname = "test";
	fr1->set_channel(oldname, newname);
	string key = fr1->get_keyword("$P3N");
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(key, newname);
	//meta data is not changed for original cp
	BOOST_CHECK_EQUAL(fr_orig->get_channels()[2], oldname);
	BOOST_CHECK_EQUAL(fr_orig->get_keyword("$P3N"), oldname);
	BOOST_CHECK_EQUAL(fr_orig->get_uri(), fr1->get_uri());
	//update data
	EVENT_DATA_VEC dat = fr1->get_data();
	float newval = 100;
	dat[100] = newval;
	fr1->set_data(dat);
	BOOST_CHECK_CLOSE(fr1->get_data()[100], newval, 1e-6);//fr1 is updated
	BOOST_CHECK_CLOSE(fr_orig->get_data()[100], newval, 1e-6);//original cp is also updated

	//delete fr_orig to ensure its copy still has valid access to h5 data
	fr_orig.reset();
	BOOST_CHECK_CLOSE(fr1->get_data()[100], newval, 1e-6);
}

BOOST_AUTO_TEST_CASE(deep_copy)
{
	//deep cp
	CytoFramePtr fr1 = cf_disk->copy();
//	fr1->copy(fr1->get_h5_file_path());
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
	BOOST_CHECK_NE(cf_disk->get_uri(), fr1->get_uri());
	//original cp is not affected by any change in meta and event data
	BOOST_CHECK_EQUAL(cf_disk->get_channels()[2], oldname);
	BOOST_CHECK_EQUAL(cf_disk->get_keyword("$P3N"), oldname);

	BOOST_CHECK_CLOSE(fr1->get_data()[100], newval, 1e-6);
	BOOST_CHECK_CLOSE(cf_disk->get_data()[100], old_val, 1e-6);



}
BOOST_AUTO_TEST_CASE(flush_meta)
{
	//deep cp
	CytoFramePtr fr1 = cf_disk->copy();
	for(auto k : cf_disk->get_keywords())
	{
		BOOST_CHECK_EQUAL(fr1->get_keyword(k.first), k.second);

	}
	string h5file = fr1->get_uri();
	//update meta data
	string oldname = fr1->get_channels()[2];
	string newname = "test";
	fr1->set_channel(oldname, newname);
	//discard change
	fr1->load_meta();
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], oldname);
	BOOST_CHECK_EQUAL(fr1->get_keyword("$P3N"), oldname);
	for(auto k : cf_disk->get_keywords())
	{
		BOOST_CHECK_EQUAL(fr1->get_keyword(k.first), k.second);

	}

	fr1->set_channel(oldname, newname);
	//flush the change
	fr1->flush_meta();
	fr1->load_meta();
	BOOST_CHECK_EQUAL(fr1->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(fr1->get_keyword("$P3N"), newname);

	for(auto k : cf_disk->get_keywords())
	{
		auto key = k.first;
		auto v1 = k.second;
		auto v2 = fr1->get_keyword(key);
		if(key!="$P3N"&&key!="$P3S")
			BOOST_CHECK_EQUAL(v1, v2);


	}

	//change it back and see if the destructor does the flushing
	fr1->set_channel(newname, oldname);
	fr1.reset();
	//load it back and see the change has NOT taken effect
	//perform shallow copy
	fr1 = load_cytoframe(h5file);


	BOOST_CHECK_EQUAL(fr1->get_channels()[2], newname);
	BOOST_CHECK_EQUAL(fr1->get_keyword("$P3N"), newname);

}

BOOST_AUTO_TEST_CASE(CytoFrameView_copy)
{
	//
	CytoFrameView fr1(cf_disk->copy());

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
