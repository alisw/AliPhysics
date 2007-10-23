#ifndef ROOT_TKDInterpolator
#define ROOT_TKDInterpolator

#ifndef ROOT_TKDTree
#include "TKDTree.h"
#endif

// Non parametric interpolation class based on local polinomial regression.
// The class can be used to approximate PDF together with TKDTree or for
// general regression when the data points are given.

template <typename Value> class TVectorT;
typedef struct TVectorT<Double_t> TVectorD;
template <typename Value> class TMatrixT;
typedef class TMatrixT<Double_t> TMatrixD;
template <typename Index, typename Value> class TKDTree;
typedef class TKDTreeIF<Int_t, Float_t> TKDTreeIF;
class TTree;
class TLinearFitter;
class TKDInterpolator : public TKDTreeIF
{
public:
	TKDInterpolator();
	TKDInterpolator(TTree *t, const Char_t *var, const Char_t *cut = 0, UInt_t bsize = 100, Long64_t nentries = 1000000000, Long64_t firstentry = 0);
	TKDInterpolator(Int_t npoints, Int_t ndim, UInt_t bsize, Float_t **data);
	~TKDInterpolator();

	        Double_t   Eval(const Double_t *point, Double_t &result, Double_t &error, Bool_t force = kFALSE);
					Float_t    GetAlpha() const {return fAlpha;}
	inline Bool_t     GetCOGPoint(Int_t node, Float_t *coord, Float_t &val, Float_t &error) const;
	inline Bool_t     GetDataPoint(Int_t n, Float_t *p) const;
	        TKDTreeIF* GetHelper() {return fKDhelper;}
	        Bool_t     GetInterpolationMethod() const {return fStatus&1;}
	        Int_t      GetNTNodes() const {return fNTNodes;}
	        Bool_t     GetStore() const {return fStatus&2;}
	        Bool_t     GetWeights() const {return fStatus&4;}
					
					void       DrawNodes(UInt_t ax1 = 0, UInt_t ax2 = 1, Int_t depth = -1);
	        void       DrawNode(Int_t tnode, UInt_t ax1 = 0, UInt_t ax2 = 1);
	        void       GetStatus();
					void       SetAlpha(const Float_t a){if(a>0.) fAlpha = a;}
	        void       SetInterpolationMethod(const Bool_t on = kTRUE);
	        void       SetStore(const Bool_t on = kTRUE);
	        void       SetWeights(const Bool_t on = kTRUE);
	
private:
	TKDInterpolator(const TKDInterpolator &);
	TKDInterpolator& operator=(const TKDInterpolator &);	
	        void       Build();
					
public:
	class TKDNodeInfo
	{
	public:
		TKDNodeInfo(const Int_t ndim = 0);
		virtual  ~TKDNodeInfo();
		void      Build(const Int_t ndim);
	  Double_t  CookPDF(const Double_t *point, Double_t &result, Double_t &error);
		void      Store(const TVectorD &par, const TMatrixD &cov);
	
		Int_t     fNDim;         // data dimension
		Float_t   *fRefPoint;    //[fNDim] node's COG
		Float_t   fRefValue;     // measured value for node 
		TMatrixD  *fCov;         // interpolator covariance matrix
		TVectorD  *fPar;         // interpolator parameters

	private:
		TKDNodeInfo(const TKDNodeInfo &);
		TKDNodeInfo& operator=(const TKDNodeInfo &);
		
		ClassDef(TKDNodeInfo, 1) // node info for interpolator
	};

protected:
	Int_t     fNTNodes;        //Number of evaluation data points
	TKDNodeInfo *fTNodes;      //[fNTNodes] interpolation node

private:
	UChar_t   fStatus;         // status of the interpolator
	UChar_t   fLambda;         // number of parameters in polynom
	Short_t		fDepth;          //! depth of the KD Tree structure used
	Float_t   fAlpha;          // alpha parameter
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
	if(!fData) return kFALSE;
		
	for(int i=0; i<fNDim; i++) p[i] = fData[i][n];
	return kTRUE;
}


#endif

