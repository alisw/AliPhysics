#ifndef ALITRDTRACKV1_H
#define ALITRDTRACKV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Represents a reconstructed TRD track                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ALIKALMANTRACK_H
#include "AliKalmanTrack.h"
#endif

#ifndef ALITRDSEEDV1_H
#include "AliTRDseedV1.h"
#endif

class AliESDtrack;
class AliTRDcluster;
class AliTRDReconstructor;
class AliTRDtrackV1 : public AliKalmanTrack
{
  friend class AliHLTTRDTrack; // allow HLT special access
public:
  enum ETRDtrackSize { 
    kNdet      = AliTRDgeometry::kNdet
   ,kNstacks   = AliTRDgeometry::kNstack*AliTRDgeometry::kNsector
   ,kNplane    = AliTRDgeometry::kNlayer
   ,kNcham     = AliTRDgeometry::kNstack
   ,kNsect     = AliTRDgeometry::kNsector
   ,kNslice    =   3
   ,kNMLPslice =   8 
   ,kMAXCLUSTERSPERTRACK = 210
  };
  
  // bits from 0-13 are reserved by ROOT (see TObject.h)
  enum ETRDtrackStatus {
    kOwner     = BIT(14)
   ,kStopped   = BIT(15) 
   ,kKink      = BIT(16) 
   ,kPrimary   = BIT(17) 
  };

  // propagation/update error codes (up to 4 bits)
  enum ETRDtrackError {
    kProlongation = 1
   ,kPropagation
   ,kAdjustSector
   ,kSnp
   ,kTrackletInit
   ,kUpdate
   ,kUnknown      = 0xff
  };

  // data/clusters/tracklet error codes (up to 4 bits/layer)
  enum ETRDlayerError {
    kGeometry = 1
   ,kBoundary
   ,kNoClusters
   ,kNoAttach
   ,kNoClustersTracklet
   ,kNoFit
   ,kChi2
  };

  AliTRDtrackV1();
  AliTRDtrackV1(AliTRDseedV1 * const trklts, const Double_t p[5], const Double_t cov[15], Double_t x, Double_t alpha);
  AliTRDtrackV1(const AliESDtrack &ref);
  AliTRDtrackV1(const AliTRDtrackV1 &ref);
  virtual ~AliTRDtrackV1();
  AliTRDtrackV1 &operator=(const AliTRDtrackV1 &ref);
  virtual void   Copy(TObject &ref) const;
 
  Bool_t         CookPID();
  Bool_t         CookLabel(Float_t wrong);
  AliTRDtrackV1* GetBackupTrack() const {return fBackupTrack;}
  Double_t       GetBudget(Int_t i) const { return fBudget[i];}
  AliTRDcluster* GetCluster(Int_t id);
  Int_t          GetClusterIndex(Int_t id) const;
  Float_t        GetEdep() const {return fDE;}
  Int_t          GetESDid() const {return fESDid;}
  inline Float_t GetMomentum(Int_t plane=-1) const;
  inline Int_t   GetNCross();
  inline Int_t   GetNumberOfTracklets() const;
  Double_t       GetPIDsignal() const   { return 0.;}
  Double_t       GetPID(Int_t is) const { return (is >=0 && is < AliPID::kSPECIES) ? fPID[is] : -1.;}
  UChar_t        GetNumberOfTrackletsPID() const;
  Double_t       GetPredictedChi2(const AliTRDseedV1 *tracklet, Double_t *cov) const;
  Double_t       GetPredictedChi2(const AliCluster* /*c*/) const                   { return 0.0; }
  Int_t          GetProlongation(Double_t xk, Double_t &y, Double_t &z) const;
  inline UChar_t GetStatusTRD(Int_t ly=-1) const;
  Int_t          GetSector() const;
  AliTRDseedV1*  GetTracklet(Int_t plane) const {return plane >=0 && plane <kNplane ? fTracklet[plane] : NULL;}
  Int_t          GetTrackletIndex(Int_t plane) const          { return (plane>=0 && plane<kNplane) ? fTrackletIndex[plane] : -1;}
  AliExternalTrackParam*
                 GetTrackIn() const  { return fTrackLow;} 
  AliExternalTrackParam*
                 GetTrackOut() const  { return fTrackHigh;} 
  const Int_t* GetTrackletIndexes() const { return &fTrackletIndex[0];}
  
  Bool_t         IsEqual(const TObject *inTrack) const;
  Bool_t         IsKink() const    { return TestBit(kKink);}
  Bool_t         IsOwner() const   { return TestBit(kOwner);};
  Bool_t         IsPrimary() const   { return TestBit(kPrimary);};
  Bool_t         IsStopped() const { return TestBit(kStopped);};
  Bool_t         IsElectron() const;
  inline static Bool_t IsTrackError(ETRDtrackError error, UInt_t status);
  inline static Bool_t IsLayerError(ETRDlayerError error, Int_t layer, UInt_t status);

  Int_t          MakeBackupTrack();
  void           Print(Option_t *o="") const;

  Bool_t         PropagateTo(Double_t xr, Double_t x0 = 8.72, Double_t rho = 5.86e-3);
  Int_t          PropagateToR(Double_t xr, Double_t step);
  Bool_t         Rotate(Double_t angle, Bool_t absolute = kFALSE);
  void           SetBudget(Int_t i, Double_t b) {if(i>=0 && i<3) fBudget[i] = b;}
  void           SetEdep(Double32_t inDE){fDE = inDE;};
  void           SetESDid(Int_t id) {fESDid = id;}
  void           SetKink(Bool_t k)        { SetBit(kKink, k);}
  void           SetPrimary(Bool_t k)     { SetBit(kPrimary, k);}
  void           SetNumberOfClusters();
  void           SetOwner();
  void           SetPID(Short_t is, Double_t inPID){if (is >=0 && is < AliPID::kSPECIES) fPID[is]=inPID;};
  void           SetPIDquality(UChar_t /*inPIDquality*/) const {/*fPIDquality = inPIDquality*/;};
  inline void    SetStatus(UChar_t stat, Int_t ly=-1);
  void           SetStopped(Bool_t stop) {SetBit(kStopped, stop);}
  void           SetTracklet(AliTRDseedV1 *const trklt,  Int_t index);
  void           SetTrackIn();
  void           SetTrackOut(const AliExternalTrackParam *op=NULL);
  inline void    SetReconstructor(const AliTRDReconstructor *rec);
  inline Float_t StatusForTOF();
  void           UnsetTracklet(Int_t plane);
  Bool_t         Update(const AliCluster *, Double_t, Int_t) { return kFALSE; };
  void           UpdateChi2(Float_t chi2);
  void           UpdateESDtrack(AliESDtrack *t);

private:
  UInt_t       fStatus;                //  Bit map for the status of propagation
  Int_t        fTrackletIndex[kNplane];//  Tracklets index in the tracker list
  Int_t        fESDid;                 //  ESD track id 
  Double32_t   fPID[AliPID::kSPECIES]; //  PID probabilities
  Double32_t   fBudget[3];             //  Integrated material budget
  Double32_t   fDE;                    //  Integrated delta energy
  const AliTRDReconstructor *fkReconstructor;//! reconstructor link 
  AliTRDtrackV1 *fBackupTrack;         //! Backup track
  AliTRDseedV1  *fTracklet[kNplane];   //  Tracklets array defining the track
  AliExternalTrackParam *fTrackLow;    // parameters of the track which enter TRD from below (TPC) 
  AliExternalTrackParam *fTrackHigh;   // parameters of the track which enter TRD from above (HMPID, PHOS) 

  ClassDef(AliTRDtrackV1, 7)          // TRD track - tracklet based
};

//____________________________________________________
inline Float_t AliTRDtrackV1::GetMomentum(Int_t plane) const
{
// Return ESD momentum stored in the tracklet reconstructed in layer = "plane". 
// By default returns the ESD momentum in first tracklet attached to track
  if(plane==-1){
    for(Int_t i(0); i<kNplane; i++){
      if(fTracklet[i]) return fTracklet[i]->GetMomentum();
    }
  } else if( plane >=0 && plane < kNplane){
    if(fTracklet[plane]) return fTracklet[plane]->GetMomentum();
  }
  return -1.;
}

//____________________________________________________
inline Int_t AliTRDtrackV1::GetNCross()
{
  Int_t ncross = 0;
  for(Int_t ip=0; ip<kNplane; ip++){
    if(!fTracklet[ip]) continue;
    ncross += fTracklet[ip]->IsRowCross();
  }
  return ncross;
}

//____________________________________________________
inline Int_t AliTRDtrackV1::GetNumberOfTracklets() const
{
  Int_t n = 0;
  for(Int_t ip=0; ip<kNplane; ip++){
    if(!fTracklet[ip]) continue;
    n++;
  }
  return n;
}

//____________________________________________________
inline UChar_t AliTRDtrackV1::GetStatusTRD(Int_t ly) const
{
  if(ly>=-1 && ly<kNplane) return (fStatus>>((ly+1)*4))&0xf;
  return kUnknown;
}

//____________________________________________________
inline Bool_t AliTRDtrackV1::IsTrackError(ETRDtrackError error, UInt_t status)
{
  return (status&0xf)==UChar_t(error);
}

//____________________________________________________
inline Bool_t AliTRDtrackV1::IsLayerError(ETRDlayerError error, Int_t ly, UInt_t status)
{
  if(ly>=kNplane || ly<0) return kFALSE;
  return ((status>>((ly+1)*4))&0xf) == UChar_t(error);
}

//____________________________________________________
inline void AliTRDtrackV1::SetReconstructor(const AliTRDReconstructor *rec)
{
  for(Int_t ip=0; ip<kNplane; ip++){
    if(!fTracklet[ip]) continue;
    fTracklet[ip]->SetReconstructor(rec);
  }
  fkReconstructor = rec;
}

//____________________________________________________
inline void AliTRDtrackV1::SetStatus(UChar_t status, Int_t ly)
{
  if(ly<kNplane) fStatus|=((status&0xf)<<((ly+1)*4));
  return;
}


//____________________________________________________________________________
inline Float_t AliTRDtrackV1::StatusForTOF()
{
  // OBSOLETE
  // Defines the status of the TOF extrapolation
  //

  if(!fTracklet[5]) return 0.;

  // Definition of res ????
  Float_t res = /*(0.2 + 0.8 * (fN / (fNExpected + 5.0))) **/ (0.4 + 0.6 * fTracklet[5]->GetN() / 20.0);
  res *= (0.25 + 0.8 * 40.0 / (40.0 + fBudget[2]));
  return res;
}

#endif



