// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/ellipse2points.hpp>
#include <stdexcept>
#include <algorithm>
#include <cmath>
using namespace std;
namespace cytolib
{

 ellipse_parsed parseEllipse(vector<float> x, vector<float> y){


// get center
	if(x.size()!=y.size())
	throw runtime_error("invalid antipodal coordinates!");

	int n = x.size();
	float mu_x = 0, mu_y = 0;
	for(auto & i : x)
	mu_x+=i;
	mu_x/=n;
	for(auto & i : y)
	mu_y+=i;
	mu_y/=n;


  //center the antipods
	for(auto & i :x)
	i-=mu_x;
	for(auto & i :y)
	i-=mu_y;

	//compute  a, b

	/*
	 * yet to be proved that L and R are always a pair of opposite antipod points (either b or a)
	 * since top vs bottom has already been proved to be
	 * not necessarily the opposite points. see the example in https://github.com/RGLab/flowWorkspace/issues/142#issuecomment-280752179
	 */
	//far left point
	int L = min_element(x.begin(), x.end()) - x.begin();
	//far right point
	int R = max_element(x.begin(), x.end()) - x.begin();

	 // use rest of two points
	vector<int> ind;
	for(int i = 0; i < n; ++i)
		if(i !=L &&i!=R)
			ind.push_back(i);




	//the lengths of two axis
	int a = sqrt(pow(x[L]-x[R],2) + pow(y[L]-y[R],2))/2;
	int b = sqrt(pow(x[ind[0]]-x[ind[1]],2) + pow(y[ind[0]]-y[ind[1]],2))/2;

	//correct a b if needed
	if(a < b)
	{
		swap(a,b);
		L = ind[0];
		R = ind[1];
	}

	//computing alpha (rotated angle)
	float alpha = atan2(y[R]-y[L],x[R]-x[L]);
	ellipse_parsed res;
	//record two antipods
//	res.x.push_back(x[L]);
//	res.y.push_back(y[L]);
//	res.x.push_back(x[B]);
//	res.y.push_back(y[B]);
	//record mu and a, b, alpha
	res.mu_x = mu_x;
	res.mu_y = mu_y;
	res.a = a;
	res.b = b;
	res.alpha = alpha;
	return res;
}
/**
 * translated from flowClust:::.ellipsePoints R code
 * @param res
 * @param n
 * @return
 */
 matrix toPoly(ellipse_parsed res, int n){


	float a = res.a;
	float b = res.b;
	float alpha = res.alpha;


	float B = min(a,b);
	float A = max(a,b);

	float d2 = (A-B)*(A+B);

	vector<float>x1(n),y1(n);
	for(int i = 0; i <n; i++)
	{
	  float phi = 2*M_PI*i/n;

	  float sp = sin(phi);
	  float cp = cos(phi);
	  float r = a*b / sqrt(B * B + d2 * sp * sp);
	  float x = r * cp;
	  float y = r * sp;
	  // xy are the ellipse points for alpha = 0 and loc = (0,0)

	  //rotate and shift
	  float ca = cos(alpha);
	  float sa = sin(alpha);
	  x1[i] = x * ca - y *sa + res.mu_x;
	  y1[i] = x * sa + y *ca + res.mu_y;


	}

	matrix m;
	m.x = x1;
	m.y = y1;

	return m;

}

};


