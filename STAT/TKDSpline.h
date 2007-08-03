#ifndef ROOT_TKDSpline
#define ROOT_TKDSpline

#ifndef ROOT_TKDInterpolator
#include "TKDInterpolator.h"
#endif

class TKDSpline : public TKDInterpolator
{
public:
	TKDSpline();
	TKDSpline(TTree *t, const Char_t *var, const Char_t *cut = 0, UInt_t bsize=100);
	TKDSpline(Int_t npoints, Int_t ndim, UInt_t bsize, Float_t **data);

private:
	void		Build();
	
protected:

private:

	ClassDef(TKDSpline, 1)   // spline fitter based on KD tree
};


#endif

