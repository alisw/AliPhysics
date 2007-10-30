#ifndef ROOT_TKDInterpolator
#define ROOT_TKDInterpolator

#ifndef ROOT_TKDInterpolatorBase
#include "TKDInterpolatorBase.h"
#endif

class TKDInterpolator : public TKDInterpolatorBase
{
public:
	TKDInterpolator();
	TKDInterpolator(Int_t ndim, Int_t npoints=0);
	~TKDInterpolator();
	void       AddNode(const TKDNodeInfo &ref);
	void       Build(Int_t ndim) {TKDInterpolatorBase::Build(ndim);}
	void       Build(Int_t npoints, Int_t ndim);
	Int_t      GetNodeIndex(const Float_t *p);
	Bool_t     SetNode(Int_t i, const TKDNodeInfo &ref);

private:
	TKDInterpolator(const TKDInterpolator &);
	TKDInterpolator& operator=(const TKDInterpolator &);	

private:
	
	ClassDef(TKDInterpolator, 1)   // LOWESS data interpolator
};


#endif

