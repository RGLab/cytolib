// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/in_polygon.hpp>
namespace cytolib
{

void in_polygon(EVENT_DATA_TYPE * xdata, EVENT_DATA_TYPE * ydata, const vector<cytolib::CYTO_POINT> & vertices, INDICE_TYPE & parentInd, bool is_negated, INDICE_TYPE &res)
{
	unsigned counter;
	EVENT_DATA_TYPE xinters;
	 //find max py
	double p_y_max = max_element(vertices.begin(), vertices.end(),[](const cytolib::CYTO_POINT & v1, const cytolib::CYTO_POINT & v2){return v1.y < v2.y;})->y;
	//TODO: potentially we can speed up by caching pre-calculated p_bottom,top,left,right here to avoid repeated computation within the loop
	for(auto i : parentInd)
	{//iterate over points
		vector<cytolib::CYTO_POINT>::const_iterator p1;
		vector<cytolib::CYTO_POINT>::const_iterator p2;
		counter=0;
		for(p1  = vertices.begin(), p2 = p1 + 1; p1 < vertices.end(); p1++,p2++)
		{// iterate over vertices

			if (p2 == vertices.end())
			{//the last vertice must "loop around"
				p2 = vertices.begin();
			}

		  /*if horizontal ray is in y range of vertex find the x coordinate where
			ray and vertex intersect*/
			const cytolib::CYTO_POINT *  p_bottom = &(*p1);
			const cytolib::CYTO_POINT *  p_top = &(*p2);
			if(p_bottom->y > p_top->y)
				swap(p_bottom, p_top);

			const cytolib::CYTO_POINT *  p_left = p_bottom;
			const cytolib::CYTO_POINT *  p_right = p_top;
			if(p_left->x > p_right->x)
				swap(p_left, p_right);

			if(ydata[i] >= p_bottom->y && ydata[i] < p_top->y &&xdata[i] <= p_right->x && p2->y != p1->y)
			{
				xinters = (ydata[i]-p1->y)*(p2->x-p1->x)/(p2->y-p1->y)+p1->x;
			  /*if intersection x coordinate == point x coordinate it lies on the
				  boundary of the polygon, which means "in"*/
				if(xinters==xdata[i])
				{
				  counter=1;
				  break;
				}
				/*count how many vertices are passed by the ray*/
				if (xinters > xdata[i])counter++;
			 }
			 else if(ydata[i] == p_y_max)//handle cell that is at the same y-level as the top vertex/edge
			 {
				if(p_top->y == ydata[i])//one end of edge reach the same y as top
				{
				  if(p_bottom->y == ydata[i])//horizontal top edge
				  {
					  counter = xdata[i] >= p_left->x && xdata[i] <= p_right->x;//whether on the edge

				  }
				  else//check if on the top vertex
					  counter = xdata[i] == p_top->x;
				}
				break;
			  }

			}
			/*uneven number of vertices passed means "in"*/

			 bool isIn =((counter % 2) != 0);
			 if(isIn != is_negated)
				res.push_back(i);
	}


}

}

