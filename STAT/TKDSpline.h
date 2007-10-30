#ifndef ROOT_TKDSpline
#define ROOT_TKDSpline

#ifndef ROOT_TKDInterpolator
#include "TKDInterpolator.h"
#endif

class TKDSpline : public TKDInterpolator
{
public:
	TKDSpline();
	TKDSpline(Int_t npoints, Int_t ndim);

private:
	void		Build(Int_t ndim = 0);
	
protected:

private:

	ClassDef(TKDSpline, 1)   // spline fitter based on KD tree
};


#endif

