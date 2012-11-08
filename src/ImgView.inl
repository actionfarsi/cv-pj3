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
	double t_b = (knownPoint.u-newPoint.u)*(knownPoint.u-newPoint.u) +
				 (knownPoint.v-newPoint.v)*(knownPoint.v-newPoint.v);
	double r_b = (knownPoint.u-refPointOffPlane->u)*(knownPoint.u-refPointOffPlane->u) +
				 (knownPoint.v-refPointOffPlane->v)*(knownPoint.v-refPointOffPlane->v);
	double v_r = (zVanish.u-refPointOffPlane->u)*(zVanish.u-refPointOffPlane->u) +
				 (zVanish.v-refPointOffPlane->v)*(zVanish.v-refPointOffPlane->v);
	double v_t = (zVanish.u-newPoint.u)*(zVanish.u-newPoint.u) +
				 (zVanish.v-newPoint.v)*(zVanish.v-newPoint.v);
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
	printf("sameZPlane() to be implemented!\n");
	fl_message("sameZPlane() to be implemented!\n");
	
	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}

