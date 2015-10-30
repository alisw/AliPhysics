#ifndef ALIITSUSEED_H
#define ALIITSUSEED_H

#include "AliExternalTrackParam.h"
#include "AliITSMFTAux.h"
class AliESDtrack;

using namespace AliITSMFTAux;


class AliITSUSeed: public AliExternalTrackParam
{
 public:
  enum {kKilled=BIT(14),kSave=BIT(15),kMark=BIT(16),kFake=BIT(17)};
  enum {kF02,kF04,kF12,kF13,kF14,kF24, kF44,kNFElem}; // non-trivial elems of propagation matrix
  enum {kK00,kK01,kK10,kK11,kK20,kK21,kK30,kK31,kK40,kK41, kNKElem}; // non-trivial elems of gain matrix
  enum {kS00,kS10,kS11,kS20,kS21,kS22,kS30,kS31,kS32,kS33,kS40,kS41,kS42,kS43,kS44,kNSElem}; // elements of 5x5 sym matrix
  enum {kR00,kR22,kNRElem};                        // non trivial elements of rotation matrix
  //
  AliITSUSeed();
  AliITSUSeed(const AliITSUSeed& src);
  AliITSUSeed &operator=(const AliITSUSeed &src);
  virtual ~AliITSUSeed();
  virtual void    Print(Option_t* option = "") const;
  void            Reset()                                {ResetBit(0xffffffff);}
  //
  void            SetLrClusterID(Int_t lr, Int_t cl);
  void            SetLr(Int_t lr)                        {SetLrClusterID(lr,-1);} // lr w/o cluster
  void            SetLrClusterID(UInt_t id)              {fClID = id;}
  void            SetParent(TObject* par)                {fParent = par;}
  void            SetChi2Cl(Double_t v)                  {fChi2Cl= v; v>0 ? fChi2Glo+=v : fChi2Penalty -= v;}
  void            Kill(Bool_t v=kTRUE)                   {SetBit(kKilled, v);}
  void            SetFake(Bool_t v=kTRUE)                {SetBit(kFake, v);}
  void            Save(Bool_t v=kTRUE)                   {SetBit(kSave,v);}
  void            SetMark(Bool_t v=kTRUE)                {SetBit(kMark,v);}
  void            FlagTree(UInt_t bits, Bool_t v=kTRUE);
  void            SetChi2ITSTPC(Float_t v)               {fChi2Match = v;}
  void            SetChi2ITSSA(Float_t v)                {fChi2ITSSA = v;}
  //
  UInt_t          GetLrClusterID()                 const {return fClID;}
  Int_t           GetLrCluster(Int_t &lr)          const {return UnpackCluster(fClID,lr);}
  Int_t           GetLayerID()                     const {return UnpackLayer(fClID);}
  Int_t           GetClusterID()                   const {return UnpackCluster(fClID);}
  Bool_t          HasClusterOnLayer(Int_t lr)      const {return fHitsPattern&(0x1<<lr);}
  Bool_t          HasCluster()                     const {return IsCluster(fClID);}
  Int_t           GetNLayersHit()                  const {return NumberOfBitsSet(fHitsPattern);}
  Int_t           GetNClusters()                   const;
  Int_t           GetClusterIndex(Int_t ind)       const;
  UShort_t        GetHitsPattern()                 const {return fHitsPattern;}
  Float_t         GetChi2Cl()                      const {return fChi2Cl;}
  Float_t         GetChi2Glo()                     const {return fChi2Glo;}
  Float_t         GetChi2Penalty()                 const {return fChi2Penalty;}
  Float_t         GetChi2GloNrm()                  const;
  Float_t         GetChi2ITSTPC()                  const {return fChi2Match;}
  Float_t         GetChi2ITSSA()                   const {return fChi2ITSSA;}
  Float_t         GetQualityVar()                  const;
  Bool_t          IsKilled()                       const {return TestBit(kKilled);}
  Bool_t          IsFake()                         const {return fClID>0&&TestBit(kFake);}
  Bool_t          IsSaved()                        const {return TestBit(kSave);}
  Bool_t          IsMarked()                       const {return TestBit(kMark);}
  Bool_t          ContainsFake()                   const;
  Int_t           FetchClusterInfo(Int_t *clIDarr) const;
  //
  Int_t           GetPoolID()                      const {return int(GetUniqueID())-1;}
  void            SetPoolID(Int_t id)                    {SetUniqueID(id+1);}
  //
  TObject*        GetParent()                      const {return fParent;}
  const AliITSUSeed* GetParent(Int_t lr)              const;
  //
  virtual Bool_t  IsSortable()                     const {return kTRUE;}
  virtual Bool_t  IsEqual(const TObject* obj)      const;
  virtual Int_t	  Compare(const TObject* obj)      const;
  //
  // test
  void            InitFromSeed(const AliExternalTrackParam* seed);
  void            ResetFMatrix();
  void            ApplyELoss2FMatrix(Double_t frac, Bool_t beforeProp);
  Bool_t          ApplyMaterialCorrection(Double_t xOverX0, Double_t xTimesRho, Double_t mass, Bool_t beforeProp);
  Bool_t          PropagateToX(Double_t xk, Double_t b);
  Bool_t          RotateToAlpha(Double_t alpha);
  Bool_t          GetTrackingXAtXAlpha(double xOther,double alpOther,double bz, double &x);
  Double_t        GetPredictedChi2(Double_t p[2],Double_t cov[3]);
  Bool_t          Update();
  Bool_t          Smooth(Double_t vecL[5],Double_t matL[15]);
  Double_t*       ProdABA(const double a[15],const double b[15]) const;
  //
  UInt_t          GetNChildren()                                 const {return fNChildren;}
  Int_t           IncChildren()                                        {return ++fNChildren;}
  Int_t           DecChildren()                                        {return --fNChildren;}
  Double_t        GetResidY()                                    const {return fResid[0];}
  Double_t        GetResidZ()                                    const {return fResid[1];}
  //
#ifdef _ITSU_TUNING_MODE_ // this is used only for tuning histo filling
  void                  SetOrdBranch(UShort_t v) {fOrdBranch=v;}
  void                  SetOrdCand(UShort_t v)   {fOrdCand=v;}
  Int_t                 GetOrdBranch()           const {return fOrdBranch;}
  Int_t                 GetOrdCand()             const {return fOrdCand;}
  Double_t              GetPullY()               const {return fResid[0]*Sqrt(fCovIYZ[0]);}
  Double_t              GetPullZ()               const {return fResid[1]*Sqrt(fCovIYZ[2]);}
#endif
  //
 protected:
  //
  UShort_t              fHitsPattern;       // bit pattern of hits
  UShort_t              fNChildren;         // number of children (prolongations)
  UInt_t                fClID;              // packed cluster info (see AliITSMFTAux::PackCluster)
  Float_t               fChi2Glo;           // current chi2 global (sum of track-cluster chi2's on layers with hit)
  Float_t               fChi2Cl;            // track-cluster chi2 (if >0) or penalty for missing cluster (if < 0)
  Float_t               fChi2Penalty;       // total penalty (e.g. for missing clusters)
  Float_t               fChi2Match;         // ITS/TPC matching chi2 (per NDF)  // RS: to move to separate object of final seed
  Float_t               fChi2ITSSA;         // ITSSA backward fit chi2 (per NDF) // RS: to move to separate object of final seed
  Double_t              fFMatrix[kNFElem];  // matrix of propagation from prev layer (non-trivial elements)
  Double_t              fKMatrix[kNKElem];  // Gain matrix non-trivial elements (note: standard MBF formula uses I-K*H)
  Double_t              fRMatrix[kNRElem];  // rotation matrix non-trivial elements
  Double_t              fCovIYZ[3];         // inverted matrix of propagation + meas errors = [Hi * Pi|i-1 * Hi^T + Ri]^-1
  Double_t              fResid[2];          // residuals vector
  TObject*              fParent;            // parent track (in higher tree hierarchy)
  //
#ifdef _ITSU_TUNING_MODE_ // this is used only for tuning histo filling
  UShort_t              fOrdBranch;         //! order in the seed branching
  UShort_t              fOrdCand;           //! order in candidates tree
#endif
  //
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

//_________________________________________________________________________
inline const AliITSUSeed* AliITSUSeed::GetParent(Int_t lr) const
{
  // get parent at given layer
  const AliITSUSeed* par=this;
  int lrt;
  while( par && (lrt=par->GetLayerID())>=0 ) {
    if (lrt==lr) break;
    par = dynamic_cast<const AliITSUSeed*>(par->GetParent());
  }
  return par;
}

//__________________________________________________________________
inline Int_t AliITSUSeed::GetNClusters() const
{
  // count number of clusters (some layers may have >1 hit)
  int ncl = 0;
  const AliITSUSeed* seed = this;
  while(seed) {
    if (seed->HasCluster()) ncl++;
    seed = (AliITSUSeed*)seed->GetParent();
  }
  return ncl;
  //
}

//__________________________________________________________________
inline void AliITSUSeed::FlagTree(UInt_t bits, Bool_t v)
{
  // set bits on all tree levels
  AliITSUSeed* seed = this;
  do {seed->SetBit(bits,v);} while ( (seed=(AliITSUSeed*)seed->GetParent()) );
  //
}

//__________________________________________________________________
inline Float_t AliITSUSeed::GetQualityVar() const
{
  // return variable for quality check (the smaller - the better)
  //  return GetChi2ITSSA()/Pt(); // this may change! RS
  return GetChi2GloNrm()/Pt();
}


#endif
