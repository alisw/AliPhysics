#include "TKDSpline.h"



ClassImp(TKDSpline)

/////////////////////////////////////////////////////////////////////
// Memory setup of protected data memebers
// fRefPoints : evaluation point of PDF for each terminal node of underlying KD Tree.
// | 1st terminal node (fNDim point coordinates) | 2nd terminal node (fNDim point coordinates) | ...
//
// fRefValues : evaluation value/error of PDF for each terminal node of underlying KD Tree.
// | 1st terminal node (value) | 2nd terminal node (value) | ... | 1st terminal node (error) | 2nd terminal node (error) | ...
/////////////////////////////////////////////////////////////////////

//_________________________________________________________________
TKDSpline::TKDSpline() : TKDInterpolator()
{
}

//_________________________________________________________________
TKDSpline::TKDSpline(Int_t npoints, Int_t ndim, UInt_t bsize, Float_t **data) : TKDInterpolator(npoints, ndim, bsize, data)
{
	Build();
}


//_________________________________________________________________
TKDSpline::TKDSpline(TTree *t, const Char_t *var, const Char_t *cut, UInt_t bsize) : TKDInterpolator(t, var, cut, bsize)
{
	Build();
}


//_________________________________________________________________
void TKDSpline::Build()
{
}

