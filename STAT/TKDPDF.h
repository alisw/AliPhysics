#ifndef ROOT_TKDPDF
#define ROOT_TKDPDF

#ifndef ROOT_TKDInterpolatorBase
#include "TKDInterpolatorBase.h"
#endif

#ifndef ROOT_TKDTree
#include "TKDTree.h"
#endif

// Non parametric interpolation class based on local polinomial regression.
// The class can be used to approximate PDF together with TKDTree or for
// general regression when the data points are given.

template <typename Index, typename Value> class TKDTree;
typedef class TKDTree<Int_t, Float_t> TKDTreeIF;
class TTree;
class TLinearFitter;
class TKDPDF : public TKDTreeIF, public TKDInterpolatorBase
{
public:
	TKDPDF();
	TKDPDF(TTree *t, const Char_t *var, const Char_t *cut = 0, UInt_t bsize = 100, Long64_t nentries = 1000000000, Long64_t firstentry = 0);
	TKDPDF(Int_t npoints, Int_t ndim, UInt_t bsize, Float_t **data);
	~TKDPDF();

	inline Bool_t     GetDataPoint(Int_t n, Float_t *p) const;
	inline Int_t      GetNodeIndex(const Float_t *p);
					void       DrawNode(Int_t tnode, UInt_t ax1=0, UInt_t ax2=1);

private:
	TKDPDF(const TKDPDF &);
	TKDPDF& operator=(const TKDPDF &);
	        void       Build(Int_t ndim = 0);

					
	ClassDef(TKDPDF, 1)   // data interpolator based on KD tree
};


//__________________________________________________________________
Bool_t	TKDPDF::GetDataPoint(Int_t n, Float_t *p) const
{
	if(n < 0 || n >= fNPoints) return kFALSE;
	if(!fData) return kFALSE;
		
	for(int i=0; i<fNDim; i++) p[i] = fData[i][n];
	return kTRUE;
}

//__________________________________________________________________
Int_t	TKDPDF::GetNodeIndex(const Float_t *p)
{
	return FindNode(p) - fNNodes;
}

#endif

