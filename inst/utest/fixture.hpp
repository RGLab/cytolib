/*
 * fixture.hpp
 *
 *  Created on: Feb 28, 2018
 *      Author: wjiang2
 */

#ifndef INST_UTEST_FIXTURE_HPP_
#define INST_UTEST_FIXTURE_HPP_
#include <boost/test/floating_point_comparison.hpp>

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




#endif /* INST_UTEST_FIXTURE_HPP_ */
