#ifndef ROOT_TKDInterpolatorBase
#define ROOT_TKDInterpolatorBase

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

///////////////////////////////////////////////////////////////
//
// Base non parametric interpolation algorithm.
// The class implements local polynomial regression (LOWESS).
// The user will work with daughter classes which implements
// particular data configurations.
//
///////////////////////////////////////////////////////////////

template <typename Value> class TVectorT;
typedef struct TVectorT<Double_t> TVectorD;
template <typename Value> class TMatrixT;
typedef class TMatrixT<Double_t> TMatrixD;
template <typename Index, typename Value> class TKDTree;
typedef class TKDTree<Int_t, Float_t> TKDTreeIF;
class TLinearFitter;
class TClonesArray;
class TKDNodeInfo;
class TKDInterpolatorBase
{
public:	
  TKDInterpolatorBase(Int_t size = 0);
  virtual ~TKDInterpolatorBase();

  Double_t   Eval(const Double_t *point, Double_t &result, Double_t &error, Bool_t force = kFALSE);
  virtual Int_t GetNodeIndex(const Float_t *p) = 0;
  Float_t    GetAlpha() const {return fAlpha;}
  Bool_t     GetCOGPoint(Int_t node, Float_t *&coord, Float_t &val, Float_t &error) const;
  Bool_t     GetInterpolationMethod() const {return fStatus&1;}
  TKDNodeInfo* GetNodeInfo(Int_t inode) const;
  Int_t      GetNTNodes() const {return fNTNodes;}
  void       GetStatus();
  Bool_t     GetStore() const {return fStatus&2;}
  Bool_t     GetWeights() const {return fStatus&4;}

  void       DrawBins(UInt_t ax1 = 0, UInt_t ax2 = 1, Float_t ax1min=-1., Float_t ax1max=1., Float_t ax2min=-1., Float_t ax2max=1.);
  void       SetAlpha(Float_t a){if(a>0.) fAlpha = a;}
  void       SetInterpolationMethod(Bool_t on = kTRUE);
  void       SetStore(Bool_t on = kTRUE);
  void       SetWeights(Bool_t on = kTRUE);

protected:
  virtual void      Build(Int_t nnodes);

private:
  TKDInterpolatorBase(const TKDInterpolatorBase &);
  TKDInterpolatorBase& operator=(const TKDInterpolatorBase &);

protected:
  Int_t         fNSize;       // data dimension
  Int_t         fNTNodes;     //Number of evaluation data points
  TClonesArray  *fTNodes;     //interpolation nodes

//private:
  UChar_t       fStatus;      // status of the interpolator
  UChar_t       fLambda;      // number of parameters in polynom
  Short_t		    fDepth;       //! depth of the KD Tree structure used
  Float_t       fAlpha;       // alpha parameter
  Float_t       **fRefPoints; //! temporary storage of COG data
  Double_t	    *fBuffer;     //! working space [2*fLambda]
  TKDTree<Int_t, Float_t> *fKDhelper;      //! kNN finder
  TLinearFitter *fFitter;     //! linear fitter	

  ClassDef(TKDInterpolatorBase, 1)   // data interpolator based on KD tree
};


#endif

