#ifndef ALIITSUSEED_H
#define ALIITSUSEED_H

#include "AliExternalTrackParam.h"
#include "AliITSUAux.h"
using namespace AliITSUAux;


class AliITSUSeed: public AliExternalTrackParam
{
 public:
  enum {kKilled=BIT(14)};
  enum {kF02,kF04,kF12,kF13,kF14,kF24, kF44,kNFElem}; // non-trivial elems of propagation matrix
  enum {kB00,kB01,kB02,kB03,kB04,kB10,kB11,kB12,kB13,kB14, kNBElem}; // non-trivial elems of B matrix (I - K*H)
  //
  AliITSUSeed();
  AliITSUSeed(const AliITSUSeed& src);
  AliITSUSeed &operator=(const AliITSUSeed &src);
  virtual ~AliITSUSeed();
  virtual void    Print(Option_t* option = "") const;
  //
  void            SetLrClusterID(Int_t lr, Int_t cl);
  void            SetLr(Int_t lr)                        {SetLrClusterID(lr,-1);} // lr w/o cluster
  void            SetLrClusterID(UInt_t id)              {fClID = id;}
  void            SetParent(TObject* par)                {fParent = par;}
  void            SetChi2Cl(Double_t v)                  {fChi2Glo += fChi2Cl= v;}
  void            Kill(Bool_t v=kTRUE)                   {SetBit(kKilled, v);}
  //
  UInt_t          GetLrClusterID()                 const {return fClID;}
  Int_t           GetLrCluster(Int_t &lr)          const {return UnpackCluster(fClID,lr);}
  Int_t           GetLayerID()                     const {return UnpackLayer(fClID);}
  Int_t           GetClusterID()                   const {return UnpackCluster(fClID);}
  Bool_t          HasClusterOnLayer(Int_t lr)      const {return fHitsPattern&(0x1<<lr);}
  Int_t           GetNLayersHit()                  const {return NumberOfBitsSet(fHitsPattern);}
  UShort_t        GetHitsPattern()                 const {return fHitsPattern;}
  Float_t         GetChi2Cl()                      const {return fChi2Cl;}
  Float_t         GetChi2Glo()                     const {return fChi2Glo;}
  Float_t         GetChi2GloNrm()                  const;
  Bool_t          IsKilled()                       const {return TestBit(kKilled);}
  //
  TObject*        GetParent()                      const {return fParent;}
  //
  virtual Bool_t  IsSortable()                     const {return kTRUE;}
  virtual Bool_t  IsEqual(const TObject* obj)      const;
  virtual Int_t	  Compare(const TObject* obj)      const;
  //
  // test
  void            ResetFMatrix();
  void            ApplyELoss2FMatrix(Double_t frac, Bool_t beforeProp);
  Bool_t          ApplyMaterialCorrection(Double_t xOverX0, Double_t xTimesRho, Double_t mass, Bool_t beforeProp);
  Bool_t          PropagateToX(Double_t xk, Double_t b);
  //
 protected:
  //
  Double_t              fFMatrix[kNFElem];  // matxif of propagation from prev layer (non-trivial elements)
  Double_t              fResid[2];          // residuals vector
  Double_t              fCombErrI[3];       // inverse combined error matrix
  Double_t              fBMatix[kNBElem];   // I - K*H matix non-trivial elements
  UShort_t              fHitsPattern;       // bit pattern of hits
  UInt_t                fClID;              // packed cluster info (see AliITSUAux::PackCluster)
  Float_t               fChi2Glo;           // current chi2 global
  Float_t               fChi2Cl;            // track-cluster chi2
  TObject*              fParent;            // parent track (in higher tree hierarchy)
  
  ClassDef(AliITSUSeed,1)
};

//_________________________________________________________________________
inline void AliITSUSeed::SetLrClusterID(Int_t lr, Int_t cl)
{
  // assign layer, cluster (if -1 - no hit on this layer)
  fClID = PackCluster(lr,cl);
  if (cl>=0) fHitsPattern |= 0x1<<lr;
}

//_________________________________________________________________________
inline void AliITSUSeed::ResetFMatrix()
{
  // reset transport matrix
  fFMatrix[kF02] = fFMatrix[kF04] = fFMatrix[kF12] = fFMatrix[kF13] = fFMatrix[kF14] = fFMatrix[kF24] = 0;
  fFMatrix[kF44] = 1.0;  // this element accumulates eloss 
}

//_________________________________________________________________________
inline Bool_t AliITSUSeed::ApplyMaterialCorrection(Double_t xOverX0, Double_t xTimesRho, Double_t mass, Bool_t beforeProp)
{
  // apply material correction and modify transport matrix
  double pold = Get1P();
  if (!CorrectForMeanMaterial(xOverX0,xTimesRho,mass)) return kFALSE;
  ApplyELoss2FMatrix( Get1P()/pold, beforeProp);
  return kTRUE;
}


//_________________________________________________________________________
inline void AliITSUSeed::ApplyELoss2FMatrix(Double_t frac, Bool_t beforeProp)
{
  // Accounts for the energy loss in the transport matrix
  // equivalent to multiplying Fmatix by E=diag{1,1,1,1,P4new/P4old}, where P4 is the 1/pt param.
  // If beforeProp is true, then it is assumed that the eloss was applied before the transport,
  // i.e. F' = F * E, otherwise, after transport, F' = E * F
  fFMatrix[kF44] *= frac;
  if (beforeProp) {
    fFMatrix[kF04] *= frac;
    fFMatrix[kF14] *= frac;
    fFMatrix[kF24] *= frac;
  }
}


#endif
