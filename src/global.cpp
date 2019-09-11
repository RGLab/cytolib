#include <cytolib/global.hpp>
#ifdef ROUT
#include <R_ext/Print.h>
#endif

#include <boost/algorithm/string.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

namespace cytolib
{
	bool my_throw_on_error = true;
	unsigned short g_loglevel = 0;

	void PRINT(string a){
	#ifdef ROUT
	 Rprintf(a.c_str());
	#else
	 cout << a;
	#endif

	}
	void PRINT(const char * a){
	#ifdef ROUT
	 Rprintf(a);
	#else
	 cout << a;
	#endif

	}
	string fs_tmp_path()
	{
		return fs::temp_directory_path();
	}
	string generate_unique_filename(const string & dir, const string & prefix, const string & suffix)
	{

		string tmp = dir + "/" + prefix + "XXXXXX" + suffix;
		int fid = mkstemps(&tmp[0], suffix.size());
		if(fid == -1)
			throw(domain_error("Can't create the unique file: " + tmp));

		close(fid);
		return tmp;
	}

	string generate_unique_dir(const string & dir, const string & prefix)
	{

		string tmp = dir + "/" + prefix + "XXXXXX";
		char * res = mkdtemp(&tmp[0]);
		if(!res)
			throw(domain_error("Can't create the unique dir: " + tmp));

		return tmp;
	}
	/**
	 * Generate time stamp as string
	 * @return
	 */
	string generate_timestamp()
	{
		time_t rawtime;
		time(&rawtime);
		char out[13];
		strftime(out, 13, "%y%m%d%H%M%S", localtime(&rawtime));
		out[12] ='\0';
		return string(out);
	}
	/**
	 * Generate uniquely identifiable id
	 * @return
	 */
	string generate_uid()
	{
		return to_string(boost::uuids::random_generator()());
	}
	string generate_uid_old(int len)//(pseudo guid)
	{
//		int t = time(NULL);//this only returns second-wise precision, not sufficient for distinguishing multiple gatingsets that are generated within short period

		chrono::milliseconds ms = chrono::duration_cast< chrono::milliseconds >(
													chrono::system_clock::now().time_since_epoch()
													);
		int t = ms.count();
//		cout << "random seed value: " << t << endl;
		srand (t);
		char s[len+1];
		static const char alphanum[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
		s[0] = alphanum[rand() % (sizeof(alphanum) - 11) + 10];
		for (int i = 1; i < len; ++i) {
			s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
		}
		s[len] = '\0';
		return string(s);
	}
	string path_dir_name(const string & full_path)
	{
//		//TODO:handle windows path
//		size_t i_last_slash = full_path.find_last_of('/');
//		return full_path.substr(0,i_last_slash);
		return fs::path(full_path).parent_path().string();
	}
	string path_base_name(const string & full_path)
	{
//		//TODO:handle windows path
//		size_t i_last_slash = full_path.find_last_of('/');
//		return full_path.substr(i_last_slash);
		return fs::path(full_path).filename().string();
	}
/**
	 * Parse the time string with fractional seconds
	 * std lib doesn't handle and boost::posix_time is not header-only
	 * @param s_time time string "H:M:S.ss"
	 * @return
	 */
	TM_ext parse_time_with_fractional_seconds(const string s_time){
		TM_ext res;
		vector<string> time_vec;
		//split the H:M:S.ms by .
		boost::split(time_vec, s_time, boost::is_any_of("."));
		//using std lib to parse the first half
		strptime(time_vec[0].c_str(), "%H:%M:%S", &(res._time));
		 //parse the second half as fractional seconds
		if(time_vec.size()==2)
		{
			res.fractional_secs = boost::lexical_cast<EVENT_DATA_TYPE>(time_vec[1]);
		}
		else
			res.fractional_secs = 0;
		return res;
	}
};

