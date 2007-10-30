#include "TKDSpline.h"


ClassImp(TKDSpline)


//_________________________________________________________________
TKDSpline::TKDSpline() :
	TKDInterpolator()
{
}

//_________________________________________________________________
TKDSpline::TKDSpline(Int_t npoints, Int_t ndim) :
	TKDInterpolator(npoints, ndim)
{
	Build();
}


//_________________________________________________________________
void TKDSpline::Build(Int_t)
{
}

