/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #2:
 * ImgView.inl (included from ImgView.cpp)
 *		contains routines for computing the 3D position of points
 ***************************************************************/

//
// TODO 4: sameXY()
//		Computes the 3D position of newPoint using knownPoint
//		that has the same X and Y coordinate, i.e. is directly
//		below or above newPoint.
//		See lecture slide on measuring heights.
//
// HINT1: make sure to dehomogenize points when necessary
// HINT2: there is a degeneracy that you should look out for involving points already in line with the reference
// HINT3: make sure to get the sign of the result right, i.e. whether it is above or below ground
void ImgView::sameXY()
{
	if (pntSelStack.size() < 2)
	{
		fl_alert("Not enough points on the stack.");
		return;
	}

	SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
	SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];

	if( !knownPoint.known() )
	{
		fl_alert("Can't compute relative values for unknown point.");
		return;
	}

	if( refPointOffPlane == NULL )
	{
		fl_alert("Need to specify the reference height first.");
		return;
	}

	/******** BEGIN TODO ********/

	// See the lecture note on measuring heights
	// using a known point directly below the new point.
	Vec3d nP = Vec3d(newPoint.u/newPoint.w, newPoint.v/newPoint.w, 1);
	Vec3d kP = Vec3d(knownPoint.u/knownPoint.w, knownPoint.v/knownPoint.w, 1);
	Vec3d rP = Vec3d(refPointOffPlane->u/refPointOffPlane->w,refPointOffPlane->u/refPointOffPlane->w, 1);

	Vec3d zV = Vec3d(zVanish.u/zVanish.w, zVanish.v/zVanish.w, 1);

	double t_b = (nP-kP)*(nP-kP)-1;
	double r_b = (kP-rP)*(kP-rP)-1;
	double v_r = (zV-rP)*(zV-rP)-1;
	double v_t = (zV-nP)*(zV-nP)-1;
	double h = sqrt(t_b/r_b * v_r/ v_t) * referenceHeight;

	newPoint.X = knownPoint.X;
	newPoint.Y = knownPoint.Y;
	newPoint.Z = h;
	newPoint.W = 1;
	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}



//
// TODO 5: sameZPlane()
//		Compute the 3D position of newPoint using knownPoint
//		that lies on the same plane and whose 3D position is known.
//		See the man on the box lecture slide.
//		If newPoint is on the reference plane (Z==0), use homography (this->H, or simply H) directly.
//
// HINT: For this function, you will only need to use the three vanishing points and the reference homography 
//       (in addition to the known 3D location of knownPoint, and the 2D location of newPoint)
void ImgView::sameZPlane()
{
	if (pntSelStack.size() < 2)
	{
		fl_alert("Not enough points on the stack.");
		return;
	}

	SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
	SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];

	if( !knownPoint.known() )
	{
		fl_alert("Can't compute relative values for unknown point.");
		return;
	}

	/******** BEGIN TODO ********/
	Mat3d mHinv = Mat3d(Hinv[0][0],Hinv[0][1],Hinv[0][2],
							Hinv[1][0],Hinv[1][1],Hinv[1][2],
							Hinv[2][0],Hinv[2][1],Hinv[2][2]);

	Mat3d mH = Mat3d(H[0][0],H[0][1],H[0][2],
					 H[1][0],H[1][1],H[1][2],
					 H[2][0],H[2][1],H[2][2]);


	Vec3d nP = Vec3d(newPoint.u, newPoint.v, newPoint.w);
	Vec3d kP = Vec3d(knownPoint.u, knownPoint.v, knownPoint.w);

	Vec3d b1 = Vec3d(newPoint.u, newPoint.v, newPoint.w);

	if (knownPoint.Z != 0){
		Vec3d xV = Vec3d(xVanish.u, xVanish.v, xVanish.w);
		Vec3d yV = Vec3d(yVanish.u, yVanish.v, yVanish.w);
		Vec3d zV = Vec3d(zVanish.u, zVanish.v, zVanish.w);

		// Vertical Line passing though newpoint	
		Vec3d vLine = cross(nP, zV);
		// Horizon
		Vec3d horizonL = cross(xV, yV);
		// Vanishing point
		Vec3d vL = cross(cross(nP, kP), horizonL);
		// 
		Vec3d b0 = mH * Vec3d(knownPoint.X, knownPoint.Y, 1);
		// Intersection between vertical line and vanishing line on reference plane
		b1 = cross(cross(b0, vL), vLine);
	}
	// Use homography
	
	nP = mHinv * b1;
	
	newPoint.X = nP[0]/nP[2];
	newPoint.Y = nP[1]/nP[2];
	newPoint.Z = knownPoint.Z;
	newPoint.W = 1;
	
	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}

