#ifndef ROOT_TKDNodeInfo
#define ROOT_TKDNodeInfo

////////////////////////////////////////////////////////
//
// Bucket representation for TKDInterpolator(Base)
//
// Author Alexandru Bercuci <A.Bercuci@gsi.de>
//
////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TBox
#include "TBox.h"
#endif

#ifndef ROOT_TMarker
#include "TMarker.h"
#endif

template <typename Value> class TVectorT;
typedef class TVectorT<Double_t> TVectorD;
template <typename Value> class TMatrixT;
typedef class TMatrixT<Double_t> TMatrixD;
class TKDNodeInfo : public TObject
{
public:
  friend class TKDPDF;
  friend class TKDInterpolatorBase;
  class TKDNodeDraw : public TBox {
  public:
    TKDNodeDraw();
    ~TKDNodeDraw() {;}
    void  Draw(Option_t* option = "");
    void  Print(const Option_t* option = "") const; // *MENU*
    void  SetNode(TKDNodeInfo*, UChar_t s, UChar_t ax1=0, UChar_t ax2=1);
  private:
    TKDNodeDraw(const TKDNodeDraw & ref);
    TKDNodeDraw& operator=(const TKDNodeDraw & ref);

    TMarker     fCOG;    // COG of the node
    TKDNodeInfo *fNode;  //! node data
    ClassDef(TKDNodeDraw, 1)   // graphical representation of TKDNodeInfo
  };

  TKDNodeInfo(Int_t ndim = 0);
  TKDNodeInfo(const TKDNodeInfo & ref);
  TKDNodeInfo& operator=(const TKDNodeInfo & ref);
  virtual ~TKDNodeInfo();
  Bool_t        CookPDF(const Double_t *point, Double_t &result, Double_t &error) const;
  inline void   GetBoundary(Int_t ax, Float_t &min, Float_t &max) const;
  inline void   GetCOG(Float_t* &p) const;
  Int_t         GetDimension() const { return fNDim/3; }
  void          GetPDF(Float_t &pdf, Float_t &error) const {pdf=fVal[0]; error=fVal[1];}
  Int_t         GetSize() const { return fNDim; }
  Int_t         GetNcov() const { return fNcov; }
  Int_t         GetNpar() const { return fNpar; }
  inline Bool_t Has(const Float_t *p) const;
  void          Print(const Option_t * = "") const;
  void          Store(TVectorD const *par, TMatrixD const *cov=NULL);

  Double_t*      Cov() const  { return fCov; }
  Double_t*      Par() const  { return fPar; }

  void          SetNode(Int_t ndim, Float_t *data, Float_t *pdf);
protected:
  void          Bootstrap();
  void          Build(Int_t ndim);
  void          SetNcov() { fNcov=Int_t(.5*fNpar*(fNpar+1.));}
  void          SetNpar() { Int_t dim=Int_t(fNDim/3); fNpar=Int_t(1 + dim + .5*dim*(dim+1));}
  Float_t*      Data() { return fData;}
  Float_t*      Val()  { return &fVal[0]; }


private:
  Int_t     fNDim;          //! 3 times data dimension
  Float_t   *fData;         //![fNDim] node's data
  Float_t   fVal[2];        //!measured value for node
  Int_t     fNpar;          //number of parameters
  Int_t     fNcov;          //number of covarince elements
  Double_t  *fPar;          //[fNpar] interpolator parameters
  Double_t  *fCov;          //[fNcov] interpolator covariance matrix


  ClassDef(TKDNodeInfo, 2)  // node info for interpolator
};

//_____________________________________________________________________
inline Bool_t TKDNodeInfo::Has(const Float_t *p) const
{
  Int_t ndim = fNDim/3;
  
  Float_t *it = &fData[ndim]; Int_t n(0);
  for(int id=0; id<ndim; id++, it+=2)
    if(p[id]>=it[0] && p[id]<it[1]) n++;
  if(n==ndim) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________
inline void TKDNodeInfo::GetBoundary(Int_t ax, Float_t &min, Float_t &max) const
{
  Int_t ndim = fNDim/3;
  if(ax<0 || ax>=ndim){
    min=0.; max=0.;
    return;
  }
  Float_t *it = &fData[ndim+(ax<<1)];
  min = it[0]; max = it[1];
}

//_____________________________________________________________________
inline void TKDNodeInfo::GetCOG(Float_t* &p) const
{
  Int_t ndim = fNDim/3;
  memcpy(p, fData, ndim*sizeof(Float_t));
}

#endif

