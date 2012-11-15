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
#include <math.h>
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
	int rot,i;
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
	//double eval;
	//double **M = nrmatrix( 1,numLines,1,3);
	//int linecount=1;
	//eigmatrix = nrmatrix(1,numLines,1,3);
	//eigenvec = nrvector(1,3);
	//	mineigvec = nrvector(1,3);
	//for (iter = lines.begin(); iter != lines.end(); iter++) {
	//	SVMPoint *e1,*e2;

	//	e1=iter->pnt1;
	//	e2=iter->pnt2;
	//	Vec3d p1=Vec3d(e1->u,e1->v,globalW);
	//	Vec3d p2=Vec3d(e2->u,e2->v,globalW);
	//	Vec3d line=cross(p1,p2);
	//	//A(linecount,0)=line[0];
	////	A(linecount,1)=line[1];
	////	A(linecount,2)=line[2];
	//	M[linecount][0]=line[0];
	//	M[linecount][1]=line[1];
	//	M[linecount][2]=line[2];
	//	linecount++;
	//}
	////find smallest eigenvalues

	//jacobi(M,numLines, eigenvec, eigmatrix, &rot);
	//	mineig=eigenvec[1];
	//	minid=1;
	//	if(mineig>eigenvec[2]){
	//		mineig=eigenvec[2];
	//		minid=2;
	//	}
	//	if(mineig>eigenvec[3]){
	//		mineig=eigenvec[3];
	//		minid=3;
	//	}
	//	printf("mineig:%d\n",mineig);
	//	printf("eigvec: ");
	//	for (i=0; i<3; i++) {

	//	mineigvec[i] = eigmatrix[i+1][minid];
	//	printf("%f\n",mineigvec[i]);
	//	}
	//MinEig(A,eval,evec);
	//define blank matrix M for best fit procedure
	double **M = nrmatrix(1,3,1,3);
	int  j;
	for (i=1; i<=3; i++)
	{
		for (j=1; j<=3; j++)
		{
			M[i][j] = 0;
		}
	}
	Vec3d holder=Vec3d(0,0,0);
	int first=0;
	int offsetW=imgWidth/2;
	int offsetH=imgHeight/2;
	//CTransform3x3 Tshift = CTransform3x3::Translation((float) -imgWidth/2, (float) -imgHeight/2);
	for (iter = lines.begin(); iter != lines.end(); iter++) {
	//1) specify each line's endpoints e1 and e2 in homogeneous coordinates
		SVMPoint *e1,*e2;
		
		e1=iter->pnt1;
		e2=iter->pnt2;
		//printf("p1:[%f,%f]\n",e1->v,e1->u);
		Vec3d p1=Vec3d(e1->u,e1->v,globalW);
		printf("p1:[%f,%f]\n",p1[0],p1[1]);
		Vec3d p2=Vec3d(e2->u,e2->v,globalW);
		//printf("p2:[%d,%d]\n",p2[0],p2[2]);
	//2) compute a homogenous coordinate vector representing the line
   // as the cross product of its two endpoints
		Vec3d c_result=cross(p1,p2);
		
		if(numLines==2){
			//printf("two lines \n");
			//check if first line
			if(first==0){
				//yes, so we keep current cross of I1
				holder=c_result;
				first=1;
			}else{
				//now we have both lines
				Vec3d twoline=cross(c_result,holder);
				//divide by w, convert to svmpoint. bestfit is defined
				
				bestfit=SVMPoint(twoline[0]/twoline[2]*globalW,twoline[1]/twoline[2]*globalW);
				printf("[%f,%f]\n",twoline[0],twoline[1]);
				//printf("w:%d \n",twoline[1]);
			}
			
		}else{
		//	printf("not two lines \n");
		//accumulate sums if more than two lines
			//M matrix extra credit
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
		eigmatrix = nrmatrix(1,3,1,3);
		eigenvec = nrvector(1,3);
		mineigvec = nrvector(1,3);
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
		//printf("mineig:%d\n",mineig);
	//	printf("eigvec: ");
		for (i=0; i<3; i++) {

		mineigvec[i] = eigmatrix[i+1][minid];
		//printf("%f\n",mineigvec[i]);
		}

		bestfit=SVMPoint(mineigvec[0]/mineigvec[2]*globalW,mineigvec[1]/mineigvec[2]*globalW);
	}
//printf("w:%d, h:%d\n",imgWidth,imgHeight); 
//printf("TODO: svmmath.cpp:61\n"); 
////fl_message("TODO: svmmath.cpp:61\n");
//
	/******** END TODO ********/
		free_nrmatrix(M, 1,numLines,1,3);
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
	int  i,j,k;
	int first=0;
	double dotsum;
	double lastsum=0;
	Vec4d final_q;
	Vec4d p,q,r;
	printf("numpts:%d\n",numPoints);
	/******** BEGIN TODO ********/
//printf("TODO: svmmath.cpp:101\n"); 
//fl_message("TODO: svmmath.cpp:101\n");
	//printf("point1:%f\n",points[0].u,points[0].v,points[0].w);
		 //choose p,q,r from the points, will get 4 or more pts to work with
		 //select first point as r
	
		p=Vec4d(points[1].X,points[1].Y,points[1].Z,points[1].W);
	//  printf("p:%f,%f,%f\n",p[0],p[1],p[2]);
			 r=Vec4d(points[0].X,points[0].Y,points[0].Z,points[0].W);
			//  printf("r:%f,%f,%f\n",r[0],r[1],r[2]);
			// line pr
			 Vec4d pr=p-r;
			 //now loop through left over points to make qr
			
			 for(k=2; k<numPoints; k++){
				 q=Vec4d(points[k].X,points[k].Y,points[k].Z,points[k].W);
				// printf("q:%f,%f,%f\n",q[0],q[1],q[2]);
				 Vec4d qr=q-r;
				  if(first==0){
						final_q=q;
						first=1;
					 }
				 //now find dot product, as close as 0
				// dotsum= pr[0]*qr[0]+ pr[1]*qr[1]+ pr[2]*qr[2];
				 int pr_m=sqrt(pr[0]*pr[0]+pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);
				 int qr_m=sqrt(qr[0]*qr[0]+qr[1]*qr[1]+qr[2]*qr[2]+qr[3]*qr[3]);
				 printf("pr:%f,%f,%f qr:%f,%f,%f\n",pr[0],pr[1],pr[2],qr[0],qr[1],qr[2]);
				// printf("result:%f",acos(pr[0]*qr[0]+pr[1]*qr[1]+pr[2]*qr[2]));
				 dotsum=(pr[0]*qr[0]+pr[1]*qr[1]+pr[2]*qr[2]+pr[3]*qr[3])/(pr_m*qr_m);
				 if(dotsum>1){
					dotsum=1;
				 }
				 dotsum=acos(dotsum)*180/3.14;
				 printf("angle: %f\n",dotsum);
					 if(dotsum>lastsum && dotsum<=90){
						lastsum=dotsum;
						final_q=q;
					 }
					
			 }
		
		 q=final_q;
		
		  Vec4d qr=q-r;
	//	printf("pr: [%f,%f], qr:[%f, %f]\n", final_p[0],final_p[1],final_q[0],final_q[1]);
		Vec4d ex=Vec4d(pr[0]*pr[0],pr[1]*pr[1],pr[2]*pr[2],pr[3]*pr[3]);// dot product of itself gives |p-r|^2
		double ex_m=sqrt(ex[0]+ex[1]+ex[2]+ex[3]); //|p-r|
		 ex=Vec4d(pr[0]/ex_m,pr[1]/ex_m,pr[2]/ex_m,pr[3]/ex_m);

		 double sum=ex[0]*qr[0]+ex[1]*qr[1]+ex[2]*qr[2]+ex[3]*qr[3];
		 Vec4d s=Vec4d(ex[0]*sum,ex[1]*sum,ex[2]*sum,ex[3]*sum);
	
		Vec4d t=qr-s;
		Vec4d ey=Vec4d(t[0]*t[0],t[1]*t[1],t[2]*t[2],t[3]*t[3]);

		double t_m=sqrt(ey[0]+ey[1]+ey[2]+ey[3]);
		ey=Vec4d(t[0]/t_m,t[1]/t_m,t[2]/t_m,t[3]/t_m);

		double min_u=FLT_MAX;
		double min_v=FLT_MAX;
		double max_u=0;
			double max_v=0;
		////loop through all points

		for(i=0; i<numPoints; i++){
			Vec4d a=Vec4d(points[i].X,points[i].Y,points[i].Z,points[i].W);
			Vec4d ar=a-r;
			
			Vec3d tmp = Vec3d(ar[0]*ex[0]+ar[1]*ex[1]+ar[2]*ex[2]+ar[3]*ex[3],ar[0]*ey[0]+ar[1]*ey[1]+ar[2]*ey[2]+ar[3]*ey[3],1);
			basisPts.push_back(tmp);
			if(min_u>basisPts[i][0]){
			
				min_u=basisPts[i][0];
			}			
			if(max_u<basisPts[i][0]){
				max_u=basisPts[i][0];
			}
			//basisPts[i][1]=ar[0]*ey[0]+ar[1]*ey[1]+ar[2]*ey[2];

			if(min_v>basisPts[i][1]){
				min_v=basisPts[i][1];
			}
			if(max_u<basisPts[i][1]){
				max_v=basisPts[i][1];
			}

		}
		uScale=max_u-min_u;
		vScale=max_v-min_v;
	
		
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
	int i,j;
	int numPoints = (int) points.size();
	printf("pts:%d",numPoints);
	assert( numPoints >= 4 );

	basisPts.clear();
	if (isRefPlane) // reference plane
	{
		printf("ref plane\n");
		for (i=0; i < numPoints; i++)
		{
			Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
			basisPts.push_back(tmp);
		}
	} 
	else // arbitrary polygon
	{
			printf("polygon\n");
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
//printf("TODO: svmmath.cpp:187\n"); 
//fl_message("TODO: svmmath.cpp:187\n");
	int n=0;
	for(j=0;j<numRows; j+=2){
		
		A(j,0)=basisPts[n][0]; //x1
		A(j,1)=basisPts[n][1]; //y1
		A(j,2)=1;
		A(j,3)=0;
		A(j,4)=0;
		A(j,5)=0;

		A(j,6)=-points[n].u*basisPts[n][0]; //-x1'*x1
		A(j,7)=-points[n].u*basisPts[n][1]; //-x1'*y1
		A(j,8)=-points[n].u; //-x1'

		//next row
		A(j+1,0)=0;
		A(j+1,1)=0;
		A(j+1,2)=0;
		A(j+1,3)=basisPts[n][0];
		A(j+1,4)=basisPts[n][1];
		A(j+1,5)=1;
		A(j+1,6)=-points[n].v*basisPts[n][0];
		A(j+1,7)=-points[n].v*basisPts[n][1];
		A(j+1,8)=-points[n].v;
		n++;

	}
	

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
