#ifndef ROOT_TKDInterpolatorBase
#define ROOT_TKDInterpolatorBase

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#ifndef ROOT_TKDNodeInfo
#include "TKDNodeInfo.h"
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
class TKDInterpolatorBase
{
public:	
  enum EKDInterpolatorBase {
    kdN = 4       // increase in the number of PDF if fit failled
   ,kNhelper = 30 // bucket size in helper kdTree
  };
  enum EKDInterpolatorBaseBits {
    kCOG   = 0  // COG interpolation method
   ,kSTORE = 1  // Store interpolation results
   ,kWEIGHTS = 2 // use weights
  };
  TKDInterpolatorBase(Int_t size = 0);
  virtual ~TKDInterpolatorBase();

  Bool_t        Bootstrap();
  Double_t      Eval(const Double_t *point, Double_t &result, Double_t &error, Bool_t force = kFALSE);
  virtual Int_t GetNodeIndex(const Float_t *p) = 0;
  Float_t       GetAlpha() const {return fAlpha;}
  Int_t         GetLambda() const {return fLambda;}
  Int_t         GetSize() const {return fNSize;}
  Bool_t        GetCOGPoint(Int_t node, Float_t *&coord, Float_t &val, Float_t &error) const;
  TKDNodeInfo*  GetNodeInfo(Int_t inode) const;
  Int_t         GetNTNodes() const;
  Bool_t        GetRange(Int_t ax, Float_t &min, Float_t &max) const;
  void          GetStatus(Option_t *opt="");

  Bool_t        HasStore() const {return TESTBIT(fStatus, kSTORE);}
  Bool_t        UseCOG() const {return TESTBIT(fStatus, kCOG);}
  Bool_t        UseWeights() const {return TESTBIT(fStatus, kWEIGHTS);}

  void          DrawProjection(UInt_t ax1 = 0, UInt_t ax2 = 1);
  void          SetAlpha(Float_t a);
  void          SetCOG(Bool_t on = kTRUE) {on ? SETBIT(fStatus, kCOG) : CLRBIT(fStatus, kCOG);}
  void          SetStore(Bool_t on = kTRUE) {on ? SETBIT(fStatus, kSTORE) : CLRBIT(fStatus, kSTORE);}
  void          SetWeights(Bool_t on = kTRUE) {on ? SETBIT(fStatus, kWEIGHTS) : CLRBIT(fStatus, kWEIGHTS);}

protected:
  virtual Bool_t    Build(Int_t nnodes);


  Int_t         fNSize;       //!data dimension
  TClonesArray  *fNodes;     //interpolation nodes
  TKDNodeInfo::TKDNodeDraw  *fNodesDraw; //!graphical representation of interpolation nodes

//private:
  UChar_t       fStatus;      // status of the interpolator
  UChar_t       fLambda;      //! number of parameters in polynom
  Short_t		    fDepth;       //! depth of the KD Tree structure used
  Float_t       fAlpha;       // parameter controlling the size of the region to interpolate n = (1+alpha)*lambda
  Float_t       **fRefPoints; //! temporary storage of COG data
  Double_t	    *fBuffer;     //! working space [2*fLambda]
  TKDTree<Int_t, Float_t> *fKDhelper;      //! kNN finder
  TLinearFitter *fFitter;     //! linear fitter	

private:
  TKDInterpolatorBase(const TKDInterpolatorBase &);
  TKDInterpolatorBase& operator=(const TKDInterpolatorBase &);

  ClassDef(TKDInterpolatorBase, 3)   // data interpolator based on KD tree
};


#endif

