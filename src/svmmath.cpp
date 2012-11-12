/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #1:
 * svmmath.cpp
 *		a routine for intersecting >2 lines (for vanishing point
 *		computation);
 *		routines for computing the homography for the reference
 *		plane and arbitrary polygons
 **************************************************************/

#pragma warning(disable : 4996)

#include "svmmath.h"
#include "jacob.h"
#include "vec.h"
#include <cstring>
#include <cstdio>
#include <assert.h>
#include <iostream>

#include "Eigen/Core"
#include "MinEig.h"

using namespace Eigen;
using namespace std;

//
// TODO 1: BestFitIntersect()
//		Given lines, the list of 3 or more lines to be intersected,
//		find the best fit intersection point.
//		See http://www-2.cs.cmu.edu/~ph/869/www/notes/vanishing.txt.
//	
SVMPoint BestFitIntersect(const std::list<SVMLine> &lines, int imgWidth, int imgHeight)
{
	// check
	if (lines.size() < 2)
	{
		fprintf(stderr, "Not enough lines to compute the best fit.");
		abort();
	}

	SVMPoint bestfit;
	list<SVMLine>::const_iterator iter;

	// To accumulate stuff
	typedef Matrix<double, Dynamic, 3, RowMajor> Matrix3;
	//mineig defs
	double **eigmatrix, *eigenvec, mineig, *mineigvec;
	int rot;
    int minid; //id for mineig
	
	int numLines = (int) lines.size();
	Matrix3 A = Matrix3::Zero(numLines, 3);	
	int globalW=(imgWidth+imgHeight)/4; //get average of width and height, then divide by 2 to get w
	// Transformation for numerical stability

	// Note: iterate through the lines list as follows:
	//		for (iter = lines.begin(); iter != lines.end(); iter++) {
	//			...iter is the pointer to the current line...
	//		}
	// Note: Function to find eigenvector with smallest eigenvalue is MinEig(A, eval, evec)
	//
	/******** BEGIN TODO ********/

	//define blank matrix M for best fit procedure
	double **M = nrmatrix(1,3,1,3);
	int i, j;
	for (i=1; i<=3; i++)
	{
		for (j=1; j<=3; j++)
		{
			M[i][j] = 0;
		}
	}
	Vec3d holder=Vec3d(0,0,0);
	int first=0;
	CTransform3x3 Tshift = CTransform3x3::Translation((float) -imgWidth/2, (float) -imgHeight/2);
	for (iter = lines.begin(); iter != lines.end(); iter++) {
	//1) specify each line's endpoints e1 and e2 in homogeneous coordinates
		SVMPoint *e1,*e2;
		//TODO: translate  by (-imageX/2, -imageY/2)
		e1=iter->pnt1;
		e2=iter->pnt2;
		
		Vec3d p1=Vec3d(e1->u,e1->v,globalW);
		
		Vec3d p2=Vec3d(e2->u,e2->v,globalW);
	//2) compute a homogenous coordinate vector representing the line
   // as the cross product of its two endpoints
		Vec3d c_result=cross(p1,p2);
		
		if(numLines==2){
			//check if first line
			if(first==0){
				//yes, so we keep current cross of I1
				holder=c_result;
				first=1;
			}else{
				//now we have both lines
				Vec3d twoline=cross(c_result,holder);
				//divide by w, convert to svmpoint. bestfit is defined
				bestfit=SVMPoint(twoline[0]/twoline[2],twoline[1]/twoline[2]);
			}
			
		}else{
		//accumulate sums if more than two lines
		M[1][1]+=c_result[0]*c_result[0];
		M[1][2]+=c_result[0]*c_result[1];
		M[1][3]+=c_result[0]*c_result[2];
		M[2][1]+=c_result[0]*c_result[1];
		M[2][2]+=c_result[1]*c_result[1];
		M[2][3]+=c_result[1]*c_result[2];
		M[3][1]+=c_result[0]*c_result[2];
		M[3][2]+=c_result[1]*c_result[2];
		M[3][3]+=c_result[2]*c_result[2];
		}

	}
	
	if(numLines!=2){
		//jacobi decomposition
		jacobi(M, 3, eigenvec, eigmatrix, &rot);
		//find smallest eigenvalues
		mineig=eigenvec[1];
		minid=1;
		if(mineig>eigenvec[2]){
			mineig=eigenvec[2];
			minid=2;
		}
		if(mineig>eigenvec[3]){
			mineig=eigenvec[3];
			minid=3;
		}
		for (i=0; i<3; i++) {
		mineigvec[i] = eigmatrix[i+1][minid];
		}

		bestfit=SVMPoint(mineigvec[0]/mineigvec[2],mineigvec[1]/mineigvec[2]);
	}

//printf("TODO: svmmath.cpp:61\n"); 
//fl_message("TODO: svmmath.cpp:61\n");

	/******** END TODO ********/
		free_nrmatrix(M,1,3,1,3);
	return bestfit;
}


//
// TODO 2: ConvertToPlaneCoordinate()
//		Given a plane defined by points, converts their coordinates into
//		plane coordinates. See the following document for more detail.
//		http://www.cs.cornell.edu/courses/cs4670/2012fa/projects/p4/homography.pdf.
//      The final divisors you apply to the u and v coordinates should be saved uScale and vScale
//
void ConvertToPlaneCoordinate(const vector<SVMPoint>& points, vector<Vec3d>& basisPts, double &uScale, double &vScale)
{
	int numPoints = points.size();

	/******** BEGIN TODO ********/
printf("TODO: svmmath.cpp:101\n"); 
fl_message("TODO: svmmath.cpp:101\n");

	/******** END TODO ********/
}



//
// TODO 3: ComputeHomography()
//		Computes the homography H from the plane specified by "points" to the image plane,
//		and its inverse Hinv.
//		If the plane is the reference plane (isRefPlane == true), don't convert the
//		coordinate system to the plane. Only do this for polygon patches where
//		texture mapping is necessary.
//		Coordinate system conversion is to be implemented in a separate routine
//		ConvertToPlaneCoordinate.
//		For more detailed explaination, see
//		http://www.cs.cornell.edu/courses/cs4670/2012fa/projects/p4/homography.pdf.
//
void ComputeHomography(CTransform3x3 &H, CTransform3x3 &Hinv, const vector<SVMPoint> &points, vector<Vec3d> &basisPts, bool isRefPlane)
{
	int i;
	int numPoints = (int) points.size();
	assert( numPoints >= 4 );

	basisPts.clear();
	if (isRefPlane) // reference plane
	{
		for (i=0; i < numPoints; i++)
		{
			Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
			basisPts.push_back(tmp);
		}
	} 
	else // arbitrary polygon
	{
        double uScale, vScale; // unused in this function
		ConvertToPlaneCoordinate(points, basisPts, uScale, vScale);
	}

	// A: 2n x 9 matrix where n is the number of points on the plane
	//    as discussed in lecture
	int numRows = 2 * numPoints;
	const int numCols = 9;

	typedef Matrix<double, Dynamic, 9, RowMajor> MatrixType;
	MatrixType A = MatrixType::Zero(numRows, numCols);

	/******** BEGIN TODO ********/
	// fill in the entries of A 
printf("TODO: svmmath.cpp:187\n"); 
fl_message("TODO: svmmath.cpp:187\n");


	/******** END TODO ********/

	double eval, h[9];
	MinEig(A, eval, h);

	H[0][0] = h[0];
	H[0][1] = h[1];
	H[0][2] = h[2];

	H[1][0] = h[3];
	H[1][1] = h[4];
	H[1][2] = h[5];

	H[2][0] = h[6];
	H[2][1] = h[7];
	H[2][2] = h[8];

	// compute inverse of H
	if (H.Determinant() == 0)
		fl_alert("Computed homography matrix is uninvertible \n");
	else
		Hinv = H.Inverse();

	int ii;
	printf("\nH=[\n");
	for (ii=0; ii<3; ii++)
		printf("%e\t%e\t%e;\n", H[ii][0]/H[2][2], H[ii][1]/H[2][2], H[ii][2]/H[2][2]);
	printf("]\nHinv=[\n");

	for (ii=0; ii<3; ii++)
		printf("%e\t%e\t%e;\n", Hinv[ii][0]/Hinv[2][2], Hinv[ii][1]/Hinv[2][2], Hinv[ii][2]/Hinv[2][2]);

	printf("]\n\n");
}
