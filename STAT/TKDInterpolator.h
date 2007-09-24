#ifndef ROOT_TKDInterpolator
#define ROOT_TKDInterpolator

#ifndef ROOT_TKDTree
#include "TKDTree.h"
#endif
#ifndef ROOT_TVectorD
#include "TVectorD.h"
#endif
#ifndef ROOT_TMatrixD
#include "TMatrixD.h"
#endif


class TTree;
class TLinearFitter;
class TKDInterpolator : public TKDTreeIF
{
public:
	struct TKDNodeInfo
	{
		TKDNodeInfo(const Int_t ndim = 0);
		~TKDNodeInfo();
		void			Build(const Int_t ndim);

		Int_t     fNDim;         // data dimension
		Float_t   *fRefPoint;    //[fNDim] node's COG
		Float_t   fRefValue;     // measured value for node 
		TMatrixD  fCov;          // interpolator covariance matrix
		TVectorD  fPar;          // interpolator parameters
		Bool_t    fPDFstatus;    // status bit for node's PDF

		ClassDef(TKDNodeInfo, 1) // node info for interpolator
	};
public:
	TKDInterpolator();
	TKDInterpolator(TTree *t, const Char_t *var, const Char_t *cut = 0, UInt_t bsize = 100, Long64_t nentries = 1000000000, Long64_t firstentry = 0);
	TKDInterpolator(Int_t npoints, Int_t ndim, UInt_t bsize, Float_t **data);
	~TKDInterpolator();

	        Double_t   Eval(const Double_t *point, Double_t &result, Double_t &error);
	inline Bool_t     GetCOGPoint(Int_t node, Float_t *coord, Float_t &val, Float_t &error) const ;
	inline Bool_t     GetDataPoint(Int_t n, Float_t *p) const;
	        TKDTreeIF* GetHelper() {return fKDhelper;}
	inline Bool_t     GetIntInterpolation() const {return fStatus&1;}
	inline Bool_t     GetSetStore() const {return fStatus&2;}
	inline Bool_t     GetUseWeights() const {return fStatus&4;}
					
					void       DrawNodes(UInt_t ax1 = 0, UInt_t ax2 = 1, Int_t depth = -1);
	        void       DrawNode(Int_t tnode, UInt_t ax1 = 0, UInt_t ax2 = 1);
	        void       GetStatus();
	        void       SetIntInterpolation(const Bool_t on = kTRUE);
	        void       SetSetStore(const Bool_t on = kTRUE);
	        void       SetUseWeights(const Bool_t on = kTRUE);
	
private:
	TKDInterpolator(const TKDInterpolator &);
	TKDInterpolator& operator=(const TKDInterpolator &);	
	        void       Build();
	        Double_t   CookPDF(const Double_t *point, const Int_t node, Double_t &result, Double_t &error);
					
protected:
	Int_t     fNTNodes;        //Number of evaluation data points
	TKDNodeInfo *fTNodes;      //[fNTNodes] interpolation node
// 	Float_t   *fRefValues;     //[fNTNodes]
// 	TMatrixD  *fCov;           //[fNTNodes] cov matrix array for nodes
// 	TVectorD  *fPar;           //[fNTNodes] parameters array for nodes
// 	Bool_t    *fPDFstatus;     //[fNTNodes] status bit for node's PDF

private:
	UChar_t   fStatus;         // status of the interpolator
	Int_t     fLambda;         // number of parameters in polynom
	Int_t			fDepth;          //! depth of the KD Tree structure used
	Float_t   **fRefPoints;    //! temporary storage of COG data
	Double_t	*fBuffer;        //! working space [2*fLambda]
	TKDTreeIF *fKDhelper;      //! kNN finder
	TLinearFitter *fFitter;    //! linear fitter	

	ClassDef(TKDInterpolator, 1)   // data interpolator based on KD tree
};

//__________________________________________________________________
Bool_t	TKDInterpolator::GetCOGPoint(Int_t node, Float_t *coord, Float_t &val, Float_t &err) const
{
	if(node < 0 || node > fNTNodes) return kFALSE;

	for(int idim=0; idim<fNDim; idim++) coord[idim] = fTNodes[node].fRefPoint[idim];
	val = fTNodes[node].fRefValue;
	err = fTNodes[node].fRefValue/TMath::Sqrt(fBucketSize);
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

