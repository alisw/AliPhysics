#ifndef ROOT_TKDInterpolator
#define ROOT_TKDInterpolator

#ifndef ROOT_TKDTree
#include "TKDTree.h"
#endif

class TTree;
class TLinearFitter;
class TKDInterpolator : public TKDTreeIF
{
public:
	TKDInterpolator();
	TKDInterpolator(TTree *t, const Char_t *var, const Char_t *cut = 0, UInt_t bsize = 100);
	TKDInterpolator(Int_t npoints, Int_t ndim, UInt_t bsize, Float_t **data);
	~TKDInterpolator();

	        TKDTreeIF* GetHelper() {return fKDhelper;}
	inline Bool_t     GetCOGPoint(Int_t node, Float_t *coord, Float_t &val, Float_t &error) const ;
	inline Bool_t     GetDataPoint(Int_t n, Float_t *p) const;
	        Double_t   Eval(const Double_t *point, Int_t npoints = 12);
	        void       DrawNodes(UInt_t ax1 = 0, UInt_t ax2 = 1, Int_t depth = -1);
	        void       DrawNode(Int_t tnode, UInt_t ax1 = 0, UInt_t ax2 = 1);

private:
	TKDInterpolator(const TKDInterpolator &);
	TKDInterpolator& operator=(const TKDInterpolator &);	
	        void       Build();
	
protected:
	Int_t     fNTNodes;        //Number of evaluation data points
	Float_t   **fRefPoints;    //[fNDim][fNTNodes]
	Float_t   *fRefValues;     //[fNTNodes]

private:
	Int_t			fDepth;        //! depth of the KD Tree structure used
	Double_t	*fTmpPoint;    //! temporary storage for one data point
	TKDTreeIF *fKDhelper;    //! kNN finder
	TLinearFitter *fFitter;  //! linear fitter	

	ClassDef(TKDInterpolator, 1)   // data interpolator based on KD tree
};

//__________________________________________________________________
Bool_t	TKDInterpolator::GetCOGPoint(Int_t node, Float_t *coord, Float_t &val, Float_t &error) const
{
	if(node < 0 || node > fNTNodes) return kFALSE;

	for(int idim=0; idim<fNDim; idim++) coord[idim] = fRefPoints[idim][node];
	val   = fRefValues[node];
	error = fRefValues[node]; // to be implemented
	return kTRUE;
}

//__________________________________________________________________
Bool_t	TKDInterpolator::GetDataPoint(Int_t n, Float_t *p) const
{
	if(n < 0 || n >= fNpoints) return kFALSE;
	
	for(int i=0; i<fNDim; i++) p[i] = fData[i][n];
	return kTRUE;
}


#endif

