#ifndef AliAODTrack_H
#define AliAODTrack_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD track implementation of AliVTrack
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TRef.h>
#include <TBits.h>

#include "AliVTrack.h"
#include "AliAODVertex.h"
#include "AliAODRedCov.h"
#include "AliAODPid.h"
#include "AliExternalTrackParam.h"
 

class AliVVertex;
class AliDetectorPID;
class AliTPCdEdxInfo;
class AliAODEvent;
class AliTOFHeader;

class AliAODTrack : public AliVTrack {

 public:
  
  enum AODTrk_t {kUndef = -1, 
		 kPrimary, 
		 kFromDecayVtx, 
		 kOrphan}; // Please note that this flag does not guarantee that the particle is a Physical Primary, it simply identifies the algorithm which was used to filter the track. In general, the following associations are used (check the filter macro to be sure, as this comment may be outdated): 
                           //kPrimary: TPC only tracks, global constrained tracks, primary tracks, kink mothers; 
                           //kFromDecayVtx: bachelor tracks from cascades, tracks from V0, kink daughters; 
                           //kUndef:TRD matched tracks

  enum AODTrkBits_t {
    kIsDCA=BIT(14),   // set if fPosition is the DCA and not the position of the first point
    kUsedForVtxFit=BIT(15), // set if this track was used to fit the vertex it is attached to
    kUsedForPrimVtxFit=BIT(16), // set if this track was used to fit the primary vertex
    kIsTPCConstrained=BIT(17), // set if this track is a SA TPC track constrained to the SPD vertex, needs to be skipped in any track loop to avoid double counting
    kIsHybridTPCCG=BIT(18), // set if this track can be used as a hybrid track i.e. Gbobal tracks with certain slecetion plus the TPC constrained tracks that did not pass the selection
    kIsGlobalConstrained=BIT(19), // set if this track is a global track constrained to the vertex, needs to be skipped in any track loop to avoid double counting
    kIsHybridGCG=BIT(20)// set if this track can be used as a hybrid track i.e. tracks with certain slecetion plus the global constraint tracks that did not pass the selection
  };


  enum AODTrkFilterBits_t {
    kTrkTPCOnly            = BIT(0), // Standard TPC only tracks
    kTrkITSsa              = BIT(1), // ITS standalone
    kTrkITSConstrained     = BIT(2), // Pixel OR necessary for the electrons
    kTrkElectronsPID       = BIT(3),    // PID for the electrons
    kTrkGlobalNoDCA        = BIT(4), // standard cuts with very loose DCA
    kTrkGlobal             = BIT(5),  // standard cuts with tight DCA cut
    kTrkGlobalSDD          = BIT(6), // standard cuts with tight DCA but with requiring the first SDD cluster instead of an SPD cluster tracks selected by this cut are exclusive to those selected by the previous cut
    kTrkTPCOnlyConstrained = BIT(7) // TPC only tracks: TPConly information constrained to SPD vertex in the filter below
  };
  

  enum AODTrkPID_t {
    kElectron     =  0,
    kMuon         =  1,
    kPion         =  2,
    kKaon         =  3,
    kProton       =  4,
    kDeuteron     =  5,
    kTriton       =  6,
    kHelium3      =  7,
    kAlpha        =  8,
    kUnknown      =  9,
    kMostProbable = -1
  };

  AliAODTrack();
  AliAODTrack(Short_t id,
	      Int_t label,
	      Double_t p[3],
	      Bool_t cartesian,
	      Double_t x[3],
	      Bool_t dca,
	      Double_t covMatrix[21],
	      Short_t q,
	      UChar_t itsClusMap,
	      AliAODVertex *prodVertex,
	      Bool_t usedForVtxFit,
	      Bool_t usedForPrimVtxFit,
	      AODTrk_t ttype=kUndef,
	      UInt_t selectInfo=0,
	      Float_t chi2perNDF = -999.);


  AliAODTrack(Short_t id,
	      Int_t label,
	      Float_t p[3],
	      Bool_t cartesian,
	      Float_t x[3],
	      Bool_t dca,
	      Float_t covMatrix[21],
	      Short_t q,
	      UChar_t itsClusMap,
	      AliAODVertex *prodVertex,
	      Bool_t usedForVtxFit,
	      Bool_t usedForPrimVtxFit,
	      AODTrk_t ttype=kUndef,
	      UInt_t selectInfo=0,
	      Float_t chi2perNDF = -999.);

  virtual ~AliAODTrack();
  AliAODTrack(const AliAODTrack& trk); 
  AliAODTrack& operator=(const AliAODTrack& trk);

  // kinematics
  virtual Double_t OneOverPt() const { return (fMomentum[0] != 0.) ? 1./fMomentum[0] : -999.; }
  virtual Double_t Phi()       const { return fMomentum[1]; }
  virtual Double_t Theta()     const { return fMomentum[2]; }
  
  virtual Double_t Px() const { return fMomentum[0] * TMath::Cos(fMomentum[1]); }
  virtual Double_t Py() const { return fMomentum[0] * TMath::Sin(fMomentum[1]); }
  virtual Double_t Pz() const { return fMomentum[0] / TMath::Tan(fMomentum[2]); }
  virtual Double_t Pt() const { return fMomentum[0]; }
  virtual Double_t P()  const { return TMath::Sqrt(Pt()*Pt()+Pz()*Pz()); }
  virtual Bool_t   PxPyPz(Double_t p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE; }

  virtual Double_t Xv() const { return GetProdVertex() ? GetProdVertex()->GetX() : -999.; }
  virtual Double_t Yv() const { return GetProdVertex() ? GetProdVertex()->GetY() : -999.; }
  virtual Double_t Zv() const { return GetProdVertex() ? GetProdVertex()->GetZ() : -999.; }
  virtual Bool_t   XvYvZv(Double_t x[3]) const { x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return kTRUE; }

  Double_t Chi2perNDF()  const { return fChi2perNDF; }

  UShort_t GetTPCnclsS(Int_t i0=0,Int_t i1=159)  const { 
    UShort_t cl = fTPCSharedMap.CountBits(i0)-fTPCSharedMap.CountBits(i1);
    return cl;
  }
  
  UShort_t GetTPCncls(Int_t i0=0,Int_t i1=159)  const { 
    UShort_t cl = fTPCFitMap.CountBits(i0)-fTPCFitMap.CountBits(i1);
    if(cl==0)cl = fTPCClusterMap.CountBits(i0)-fTPCClusterMap.CountBits(i1);// backward compatibility
    return cl;
  }
  
  UShort_t GetTPCNcls()  const { return GetTPCncls(); }

  Int_t GetNcls(Int_t idet) const;

  virtual Double_t M() const { return M(GetMostProbablePID()); }
  Double_t M(AODTrkPID_t pid) const;
  virtual Double_t E() const { return E(GetMostProbablePID()); }
  Double_t E(AODTrkPID_t pid) const;
  Double_t E(Double_t m) const { return TMath::Sqrt(P()*P() + m*m); }
  virtual Double_t Y() const { return Y(GetMostProbablePID()); }
  Double_t Y(AODTrkPID_t pid) const;
  Double_t Y(Double_t m) const;
  
  virtual Double_t Eta() const { return -TMath::Log(TMath::Tan(0.5 * fMomentum[2])); }

  virtual Short_t  Charge() const {return fCharge; }

  virtual Bool_t   PropagateToDCA(const AliVVertex *vtx, 
	  Double_t b, Double_t maxd, Double_t dz[2], Double_t covar[3]);

  // PID
  virtual const Double_t *PID() const { return fPID; }
  AODTrkPID_t GetMostProbablePID() const;
  void ConvertAliPIDtoAODPID();
  void SetDetPID(AliAODPid *aodpid) {fDetPid = aodpid;}

  void     SetPIDForTracking(Int_t pid) {fPIDForTracking = pid;}
  Int_t    GetPIDForTracking()  const   {return fPIDForTracking;}
  Double_t GetMassForTracking() const;

  template <typename T> void GetPID(T *pid) const {
    for(Int_t i=0; i<10; ++i) pid[i] = fPID ? fPID[i]:0;}
 
  template <typename T> void SetPID(const T *pid) {
    if (pid) {
      if (!fPID) fPID = new Double32_t[10];
      for(Int_t i=0; i<10; ++i) fPID[i]=pid[i];
    }
    else {delete[] fPID; fPID = 0;}
  }
  
  Bool_t IsOn(Int_t mask) const {return (fFlags&mask)>0;}
  ULong_t GetStatus() const { return GetFlags(); }
  ULong_t GetFlags() const { return fFlags; }

  Int_t   GetID() const { return (Int_t)fID; }
  Int_t   GetLabel() const { return fLabel; } 
  void    GetTOFLabel(Int_t *p) const;


  Char_t  GetType() const { return fType;}
  Bool_t  IsPrimaryCandidate() const;
  Bool_t  GetUsedForVtxFit() const { return TestBit(kUsedForVtxFit); }
  Bool_t  GetUsedForPrimVtxFit() const { return TestBit(kUsedForPrimVtxFit); }

  Bool_t  IsHybridGlobalConstrainedGlobal() const { return TestBit(kIsHybridGCG); }
  Bool_t  IsHybridTPCConstrainedGlobal() const { return TestBit(kIsHybridTPCCG); }
  Bool_t  IsTPCOnly() const { return IsTPCConstrained(); } // obsolete bad naming
  Bool_t  IsTPCConstrained() const { return TestBit(kIsTPCConstrained); }
  Bool_t  IsGlobalConstrained() const { return TestBit(kIsGlobalConstrained); }
  //
  Int_t   GetTOFBunchCrossing(Double_t b=0, Bool_t tpcPIDonly=kFALSE) const;
  //
  using AliVTrack::GetP;
  template <typename T> void GetP(T *p) const {
    p[0]=fMomentum[0]; p[1]=fMomentum[1]; p[2]=fMomentum[2];}

//  template <typename T> void GetPxPyPz(T *p) const {
//    p[0] = Px(); p[1] = Py(); p[2] = Pz();}
  Bool_t GetPxPyPz(Double_t *p) const;

  template <typename T> Bool_t GetPosition(T *x) const {
    x[0]=fPosition[0]; x[1]=fPosition[1]; x[2]=fPosition[2];
    return TestBit(kIsDCA);}

  template <typename T> void SetCovMatrix(const T *covMatrix) {
    if(!fCovMatrix) fCovMatrix=new AliAODRedCov<6>();
    fCovMatrix->SetCovMatrix(covMatrix);}

  template <typename T> Bool_t GetCovMatrix(T *covMatrix) const {
    if(!fCovMatrix) return kFALSE;
    fCovMatrix->GetCovMatrix(covMatrix); return kTRUE;}

  Bool_t GetXYZ(Double_t *p) const {
    return GetPosition(p); }  
  
  Bool_t GetXYZAt(Double_t x, Double_t b, Double_t *r) const;
  Bool_t GetXYZatR(Double_t xr,Double_t bz, Double_t *xyz=0, Double_t* alpSect=0) const;  

  Bool_t GetCovarianceXYZPxPyPz(Double_t cv[21]) const {
    return GetCovMatrix(cv);}

  void RemoveCovMatrix() {delete fCovMatrix; fCovMatrix=NULL;}

  Double_t XAtDCA() const { return fPositionAtDCA[0]; } //makes sense only for constrained tracks, returns dummy values for all other tracks
  Double_t YAtDCA() const { return fPositionAtDCA[1]; } //makes sense only for constrained tracks, returns dummy values for all other tracks
  Double_t ZAtDCA() const {
    if (IsMuonTrack()) return fPosition[2];
    else if (TestBit(kIsDCA)) return fPosition[1];
    else return -999.; }                                //makes sense only for constrained tracks, returns dummy values for all other tracks
  Bool_t   XYZAtDCA(Double_t x[3]) const { x[0] = XAtDCA(); x[1] = YAtDCA(); x[2] = ZAtDCA(); return kTRUE; }
  
  Double_t DCA() const {
    if (IsMuonTrack()) return TMath::Sqrt(XAtDCA()*XAtDCA() + YAtDCA()*YAtDCA());
    else if (TestBit(kIsDCA)) return fPosition[0];
    else return -999.; }
  
  Double_t PxAtDCA() const { return fMomentumAtDCA[0]; } //makes sense only for constrained tracks, returns dummy values for all other tracks
  Double_t PyAtDCA() const { return fMomentumAtDCA[1]; } //makes sense only for constrained tracks, returns dummy values for all other tracks
  Double_t PzAtDCA() const { return fMomentumAtDCA[2]; } //makes sense only for constrained tracks, returns dummy values for all other tracks
  Double_t PAtDCA() const { return TMath::Sqrt(PxAtDCA()*PxAtDCA() + PyAtDCA()*PyAtDCA() + PzAtDCA()*PzAtDCA()); }
  Bool_t   PxPyPzAtDCA(Double_t p[3]) const { p[0] = PxAtDCA(); p[1] = PyAtDCA(); p[2] = PzAtDCA(); return kTRUE; }
  
  Double_t GetRAtAbsorberEnd() const { return fRAtAbsorberEnd; }
  
  Double_t GetITSchi2()       const       {return fITSchi2;}
  UChar_t  GetITSClusterMap() const       { return (UChar_t)(fITSMuonClusterMap&0xff); }
  UChar_t  GetITSSharedClusterMap() const { return (UChar_t)((fITSMuonClusterMap&0xff00)>>8); }
  Int_t    GetITSNcls() const; 
  Bool_t   HasPointOnITSLayer(Int_t i) const { return TESTBIT(GetITSClusterMap(),i); }
  Bool_t   HasSharedPointOnITSLayer(Int_t i) const { return TESTBIT(GetITSSharedClusterMap(),i); }
  UInt_t   GetMUONClusterMap() const      { return (fITSMuonClusterMap&0x3ff0000)>>16; }
  UInt_t   GetITSMUONClusterMap() const   { return fITSMuonClusterMap; }
  
  Bool_t  TestFilterBit(UInt_t filterBit) const {return (Bool_t) ((filterBit & fFilterMap) != 0);}
  Bool_t  TestFilterMask(UInt_t filterMask) const {return (Bool_t) ((filterMask & fFilterMap) == filterMask);}
  void    SetFilterMap(UInt_t i){fFilterMap = i;}
  UInt_t  GetFilterMap() const {return fFilterMap;}

  const TBits& GetTPCClusterMap() const {return fTPCClusterMap;}
  const TBits* GetTPCClusterMapPtr() const {return &fTPCClusterMap;}
  const TBits& GetTPCFitMap() const {return fTPCFitMap;}
  const TBits* GetTPCFitMapPtr() const {return &fTPCFitMap;}
  Float_t GetTPCClusterInfo(Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159, Int_t /*type*/=0) const;
  
  const TBits& GetTPCSharedMap() const {return fTPCSharedMap;}
  const TBits* GetTPCSharedMapPtr() const {return &fTPCSharedMap;}
  void    SetTPCClusterMap(const TBits amap) {fTPCClusterMap = amap;}
  void    SetTPCSharedMap(const TBits amap) {fTPCSharedMap = amap;}
  void    SetTPCFitMap(const TBits amap) {fTPCFitMap = amap;}
  void    SetTPCPointsF(UShort_t  findable){fTPCnclsF = findable;}
  void    SetTPCNCrossedRows(UInt_t n)     {fTPCNCrossedRows = n;}
  
  virtual const    AliExternalTrackParam * GetInnerParam() const { return NULL; }
  virtual const    AliExternalTrackParam * GetOuterParam() const { return NULL; }

  UShort_t GetTPCNclsF() const { return fTPCnclsF;}
  UShort_t GetTPCNCrossedRows()  const { return fTPCNCrossedRows;}
  Float_t GetTPCCrossedRows() const {return (Float_t) GetTPCNCrossedRows();}
  Float_t  GetTPCFoundFraction() const { return fTPCNCrossedRows>0 ? float(GetTPCNcls())/fTPCNCrossedRows : 0;}

  // Calorimeter Cluster
  Int_t GetEMCALcluster() const {return fCaloIndex;}
  void SetEMCALcluster(Int_t index) {fCaloIndex=index;}
  Bool_t IsEMCAL() const {return fFlags&kEMCALmatch;}

  Double_t GetTrackPhiOnEMCal() const {return fTrackPhiOnEMCal;}
  Double_t GetTrackEtaOnEMCal() const {return fTrackEtaOnEMCal;}
  Double_t GetTrackPtOnEMCal() const {return fTrackPtOnEMCal;}
  Double_t GetTrackPOnEMCal() const {return TMath::Abs(fTrackEtaOnEMCal) < 1 ? fTrackPtOnEMCal*TMath::CosH(fTrackEtaOnEMCal) : -999;}
  void SetTrackPhiEtaPtOnEMCal(Double_t phi,Double_t eta,Double_t pt) {fTrackPhiOnEMCal=phi;fTrackEtaOnEMCal=eta;fTrackPtOnEMCal=pt;}

  Int_t GetPHOScluster() const {return fCaloIndex;}
  void SetPHOScluster(Int_t index) {fCaloIndex=index;}
  Bool_t IsPHOS() const {return fFlags&kPHOSmatch;}

  //pid signal interface
  Double_t  GetITSsignal()       const { return fDetPid?fDetPid->GetITSsignal():0.;    }
  void      GetITSdEdxSamples(Double_t s[4]) const;
  Double_t  GetITSsignalTunedOnData() const {return fITSsignalTuned ;}
  void      SetITSsignalTunedOnData(Double_t signal) {fITSsignalTuned = signal;}
  Double_t  GetTPCsignal()       const { return fDetPid?fDetPid->GetTPCsignal():0.;    }
  Double_t  GetTPCsignalTunedOnData() const { return fTPCsignalTuned;}
  void      SetTPCsignalTunedOnData(Double_t signal) {fTPCsignalTuned = signal;}
  UShort_t  GetTPCsignalN()      const { return fDetPid?fDetPid->GetTPCsignalN():0;    }
  virtual AliTPCdEdxInfo* GetTPCdEdxInfo() const {return fDetPid?fDetPid->GetTPCdEdxInfo():0;}
  Double_t  GetTPCmomentum()     const { return fDetPid?fDetPid->GetTPCmomentum():0.;  }
  Double_t  GetTPCTgl()          const { return fDetPid?fDetPid->GetTPCTgl():0.;  }
  Double_t  GetTOFsignal()       const { return fDetPid?fDetPid->GetTOFsignal():0.;    }
  Double_t  GetIntegratedLength() const { return fTrackLength;}
  void      SetIntegratedLength(Double_t l) {fTrackLength = l;}
  Double_t  GetTOFsignalTunedOnData() const { return fTOFsignalTuned;}
  void      SetTOFsignalTunedOnData(Double_t signal) {fTOFsignalTuned = signal;}
  Double_t  GetHMPIDsignal()     const; 
  Double_t  GetHMPIDoccupancy()  const;

  Int_t     GetHMPIDcluIdx()     const;
    
  void GetHMPIDtrk(Float_t &x, Float_t &y, Float_t &th, Float_t &ph) const;  
  void GetHMPIDmip(Float_t &x,Float_t &y,Int_t &q,Int_t &nph) const;
  
  Bool_t GetOuterHmpPxPyPz(Double_t *p) const;
  
  void      GetIntegratedTimes(Double_t *times, Int_t nspec=AliPID::kSPECIES) const {if (fDetPid) fDetPid->GetIntegratedTimes(times, nspec);}
  Double_t  GetTRDslice(Int_t plane, Int_t slice) const;
  Double_t  GetTRDsignal()                        const {return fDetPid ? fDetPid->GetTRDsignal() : 0;}
  Double_t  GetTRDmomentum(Int_t plane, Double_t */*sp*/=0x0) const;
  Double_t  GetTRDchi2()                 const {return fDetPid ? fDetPid->GetTRDChi2() : -1;}
  UChar_t   GetTRDncls(Int_t layer)      const;
  UChar_t   GetTRDncls()                 const {return GetTRDncls(-1);}
  UChar_t   GetTRDntrackletsPID() const;
  Int_t     GetNumberOfTRDslices() const { return fDetPid?fDetPid->GetTRDnSlices():0; }
  void      GetHMPIDpid(Double_t */*p*/) const { return; } // TODO: To be implemented properly with the new HMPID object

  void SetMFTClusterPattern(ULong_t mftClusterPattern) { fMFTClusterPattern = mftClusterPattern; }   // AU
  ULong_t GetMFTClusterPattern() { return fMFTClusterPattern; }                                      // AU

  const AliAODEvent* GetAODEvent() const {return fAODEvent;}
  virtual const AliVEvent* GetEvent() const {return (AliVEvent*)fAODEvent;}
  void SetAODEvent(const AliAODEvent* ptr){fAODEvent = ptr;}
  const AliTOFHeader* GetTOFHeader() const;

  AliAODPid    *GetDetPid() const { return fDetPid; }
  AliAODVertex *GetProdVertex() const { return (AliAODVertex*)fProdVertex.GetObject(); }
  
  // print
  void  Print(const Option_t *opt = "") const;

  // setters
  void SetFlags(ULong_t flags) { fFlags = flags; }
  void SetStatus(ULong_t flags) { fFlags|=flags; }
  void ResetStatus(ULong_t flags) { fFlags&=~flags; }

  void SetID(Short_t id) { fID = id; }
  void SetLabel(Int_t label) { fLabel = label; }
  void SetTOFLabel(const Int_t* p);
  template <typename T> void SetPosition(const T *x, Bool_t isDCA = kFALSE);
  template <typename T> void SetP(const T *p, const Bool_t cartesian);
  void SetDCA(Double_t d, Double_t z);
  void SetUsedForVtxFit(Bool_t used = kTRUE) { used ? SetBit(kUsedForVtxFit) : ResetBit(kUsedForVtxFit); }
  void SetUsedForPrimVtxFit(Bool_t used = kTRUE) { used ? SetBit(kUsedForPrimVtxFit) : ResetBit(kUsedForPrimVtxFit); }

  void SetIsTPCOnly(Bool_t b = kTRUE) { SetIsTPCConstrained(b); }// obsolete bad naming

  void SetIsTPCConstrained(Bool_t b = kTRUE) { b ? SetBit(kIsTPCConstrained) : ResetBit(kIsTPCConstrained); }
  void SetIsHybridTPCConstrainedGlobal(Bool_t hybrid = kTRUE) { hybrid ? SetBit(kIsHybridTPCCG) : ResetBit(kIsHybridTPCCG); }

  void SetIsGlobalConstrained(Bool_t b = kTRUE) { b ? SetBit(kIsGlobalConstrained) : ResetBit(kIsGlobalConstrained); }
  void SetIsHybridGlobalConstrainedGlobal(Bool_t hybrid = kTRUE) { hybrid ? SetBit(kIsHybridGCG) : ResetBit(kIsHybridGCG); }



  void SetOneOverPt(Double_t oneOverPt) { fMomentum[0] = 1. / oneOverPt; }
  void SetPt(Double_t pt) { fMomentum[0] = pt; };
  void SetPhi(Double_t phi) { fMomentum[1] = phi; }
  void SetTheta(Double_t theta) { fMomentum[2] = theta; }
  void SetP() {fMomentum[0]=fMomentum[1]=fMomentum[2]=-999.;}

  void SetXYAtDCA(Double_t x, Double_t y) {fPositionAtDCA[0] = x; fPositionAtDCA[1] = y;}
  void SetPxPyPzAtDCA(Double_t pX, Double_t pY, Double_t pZ) {fMomentumAtDCA[0] = pX; fMomentumAtDCA[1] = pY; fMomentumAtDCA[2] = pZ;}
  
  void SetRAtAbsorberEnd(Double_t r) { fRAtAbsorberEnd = r; }
  
  void SetCharge(Short_t q) { fCharge = q; }
  void SetChi2perNDF(Double_t chi2perNDF) { fChi2perNDF = chi2perNDF; }

  void SetITSchi2(Double_t ITSchi2)                         {fITSchi2 = ITSchi2;}
  void SetITSClusterMap(UChar_t itsClusMap)                 { fITSMuonClusterMap = (fITSMuonClusterMap&0xffffff00)|(((UInt_t)itsClusMap)&0xff); }
  void SetITSSharedMap(UChar_t map)                         { fITSMuonClusterMap = (fITSMuonClusterMap&0xffff00ff)|((((UInt_t)map)&0xff)<<8); }
  void SetMuonClusterMap(UInt_t muonClusMap)                { fITSMuonClusterMap = (fITSMuonClusterMap&0xfc00ffff)|((muonClusMap&0x3ff)<<16); }
  void SetITSMuonClusterMap(UInt_t itsMuonClusMap)          { fITSMuonClusterMap = itsMuonClusMap; }
  void SetMUONtrigHitsMapTrg(UInt_t muonTrigHitsMap) { fMUONtrigHitsMapTrg = muonTrigHitsMap; }
  UInt_t GetMUONTrigHitsMapTrg() const { return fMUONtrigHitsMapTrg; }
  void SetMUONtrigHitsMapTrk(UInt_t muonTrigHitsMap) { fMUONtrigHitsMapTrk = muonTrigHitsMap; }
  UInt_t GetMUONTrigHitsMapTrk() const { return fMUONtrigHitsMapTrk; }
  Int_t GetMuonTrigDevSign() const;

  Int_t GetMatchTrigger() const {return fITSMuonClusterMap>>30;}
					//  0 Muon track does not match trigger
					//  1 Muon track match but does not pass pt cut
					//  2 Muon track match Low pt cut
					//  3 Muon track match High pt cut
  void     SetMatchTrigger(Int_t MatchTrigger);
  Bool_t   MatchTrigger() const { return (GetMatchTrigger()>0); }	  //  Muon track matches trigger track
  Bool_t   MatchTriggerLowPt()   const  { return (GetMatchTrigger()>1); } //  Muon track matches trigger track and passes Low pt cut
  Bool_t   MatchTriggerHighPt()  const  { return (GetMatchTrigger()>2); } //  Muon track matches trigger track and passes High pt cut
  Bool_t   MatchTriggerDigits()  const;                                   //  Muon track matches trigger digits
  Double_t GetChi2MatchTrigger() const  { return fChi2MatchTrigger;}
  void     SetChi2MatchTrigger(Double_t Chi2MatchTrigger) {fChi2MatchTrigger = Chi2MatchTrigger; }
  Bool_t   HitsMuonChamber(Int_t MuonChamber, Int_t cathode = -1) const;  // Check if track hits Muon chambers
  Bool_t   IsMuonTrack() const { return ( (GetMUONClusterMap()>0) && !fIsMuonGlobalTrack ) ? kTRUE : kFALSE; }
  
  Bool_t   IsMuonGlobalTrack() const { return fIsMuonGlobalTrack; }                                     // AU
  void     SetIsMuonGlobalTrack(Bool_t isMuonGlobalTrack) { fIsMuonGlobalTrack = isMuonGlobalTrack; }   // AU

  void     Connected(Bool_t flag) {flag ? SETBIT(fITSMuonClusterMap,26) : CLRBIT(fITSMuonClusterMap,26);}
  Bool_t   IsConnected() const {return TESTBIT(fITSMuonClusterMap,26);}

  void     SetProdVertex(TObject *vertex) { fProdVertex = vertex; }
  void     SetType(AODTrk_t ttype) { fType=ttype; }

  // Trasient PID object, is owned by the track
  virtual void  SetDetectorPID(const AliDetectorPID *pid);
  virtual const AliDetectorPID* GetDetectorPID() const { return fDetectorPID; }

  // Dummy
  Int_t    PdgCode() const {return 0;}
  
 private :

  // Momentum & position
  Double32_t    fMomentum[3];       // momemtum stored in pt, phi, theta
  Double32_t    fPosition[3];       // position of first point on track or dca
  
  Double32_t    fMomentumAtDCA[3];  // momentum (px,py,pz) at DCA
  Double32_t    fPositionAtDCA[2];  // trasverse position (x,y) at DCA
  
  Double32_t    fRAtAbsorberEnd;    // transverse position r at the end of the muon absorber
  
  Double32_t    fChi2perNDF;        // chi2/NDF of momentum fit
  Double32_t    fChi2MatchTrigger;  // chi2 of trigger/track matching
  Double32_t*   fPID;               //! [0.,1.,8] pointer to PID object

  Double32_t    fITSchi2;           // ITS chi2

  ULong_t       fFlags;             // reconstruction status flags 
  Int_t         fLabel;             // track label, points back to MC track
  Int_t         fTOFLabel[3];       // TOF label
  Double32_t    fTrackLength;       // Track length
  UInt_t        fITSMuonClusterMap; // map of ITS and muon clusters, one bit per layer
                                    // (ITS: bit 1-8, muon trigger or ITS shared: bit 9-16, muon tracker: bit 17-26, muon match trigger: bit 31-32) 
  UInt_t        fMUONtrigHitsMapTrg; // Muon trigger hits map from trigger
  UInt_t        fMUONtrigHitsMapTrk; // Muon trigger hits map from tracker track extrapolation
  UInt_t        fFilterMap;         // filter information, one bit per set of cuts

  TBits         fTPCFitMap;      // Map of clusters, one bit per padrow; if has a cluster on given padrow which is used in the fit   
  TBits         fTPCClusterMap;     // Map of clusters, one bit per padrow; 1 if has a cluster on given padrow
  TBits         fTPCSharedMap;      // Map of clusters, one bit per padrow; 1 if has a shared cluster on given padrow

  UShort_t      fTPCnclsF;          // findable clusters
  UShort_t      fTPCNCrossedRows;   // n crossed rows

  Short_t       fID;                // unique track ID, points back to the ESD track

  Char_t        fCharge;            // particle charge
  Char_t        fType;              // Track Type, explanation close to the enum AODTrk_t

  Char_t        fPIDForTracking;    // pid using for tracking of ESD track

  Int_t         fCaloIndex;         // index of associated EMCAL/PHOS cluster (AliAODCaloCluster)

  
  AliAODRedCov<6> *fCovMatrix;      // covariance matrix (x, y, z, px, py, pz)
  AliAODPid    *fDetPid;            // more detailed or detector specific raw pid information
  mutable const AliDetectorPID* fDetectorPID; //! transient object to cache calibrated PID information
  TRef          fProdVertex;        // vertex of origin

  Double32_t    fTrackPhiOnEMCal;   // phi of track after being propagated to the EMCal surface (default r = 440 cm)
  Double32_t    fTrackEtaOnEMCal;   // eta of track after being propagated to the EMCal surface (default r = 440 cm)
  Double32_t    fTrackPtOnEMCal;    // pt of track after being propagated to the EMCal surface (default r = 440 cm)

  Bool_t fIsMuonGlobalTrack;        // True if the track is built from the combination of MUON and MFT clusters     // AU

  Double32_t    fITSsignalTuned;    //! ITS signal tuned on data when using MC
  Double32_t    fTPCsignalTuned;    //! TPC signal tuned on data when using MC
  Double32_t    fTOFsignalTuned;    //! TOF signal tuned on data when using MC

  ULong_t fMFTClusterPattern;       // Tells us which MFT clusters are contained in the track, and which one is a good one (if MC)  // AU

  const AliAODEvent* fAODEvent;     //! pointer back to the event the track belongs to

  //---------------------------------------------------------------------------
  //--the calibration interface--
  //--to be used in online calibration/QA
  //--should also be implemented in ESD so it works offline as well
  //-----------
  virtual Int_t GetTrackParam         ( AliExternalTrackParam &p ) const;
  virtual Int_t GetTrackParamRefitted ( AliExternalTrackParam &p ) const;
  virtual Int_t GetTrackParamIp       ( AliExternalTrackParam &p ) const;
  virtual Int_t GetTrackParamTPCInner ( AliExternalTrackParam &p ) const;
  virtual Int_t GetTrackParamOp       ( AliExternalTrackParam &p ) const;
  virtual Int_t GetTrackParamCp       ( AliExternalTrackParam &p ) const;
  virtual Int_t GetTrackParamITSOut   ( AliExternalTrackParam &p ) const;
  Int_t GetNumberOfITSClusters() const { return GetITSNcls();}
  Int_t GetNumberOfTPCClusters() const { return GetTPCncls();}  
  Int_t GetNumberOfTRDClusters() const { return GetTRDncls();}  

  ClassDef(AliAODTrack, 25);
};

inline Bool_t  AliAODTrack::IsPrimaryCandidate() const
{
    // True of track passes primary particle selection (independent of type) 
    // 
    if (fFilterMap) {
	return kTRUE;
    } else {
	return kFALSE;
    }
}

inline Int_t AliAODTrack::GetITSNcls() const 
{
  // Number of points in ITS
  Int_t n=0;
  for(Int_t i=0;i<6;i++) if(HasPointOnITSLayer(i)) n++;
  return n;
}

//______________________________________________________________________________
template <typename T> 
void AliAODTrack::SetPosition(const T *x, const Bool_t dca) 
{
  // set the position

  if (x) {
    if (!dca) {
      ResetBit(kIsDCA);

      fPosition[0] = x[0];
      fPosition[1] = x[1];
      fPosition[2] = x[2];
    } else {
      SetBit(kIsDCA);
      // don't know any better yet
      fPosition[0] = -999.;
      fPosition[1] = -999.;
      fPosition[2] = -999.;
    }
  } else {
    ResetBit(kIsDCA);

    fPosition[0] = -999.;
    fPosition[1] = -999.;
    fPosition[2] = -999.;
  }
}


//______________________________________________________________________________
template <typename T> void AliAODTrack::SetP(const T *p, const Bool_t cartesian) 
{
  // Set the momentum

  if (p) {
    if (cartesian) {
      Double_t pt2 = p[0]*p[0] + p[1]*p[1];
      Double_t pp  = TMath::Sqrt(pt2 + p[2]*p[2]);
      
      fMomentum[0] = TMath::Sqrt(pt2); // pt
      fMomentum[1] = (pt2 != 0.) ? TMath::Pi()+TMath::ATan2(-p[1], -p[0]) : -999; // phi
      fMomentum[2] = (pp != 0.) ? TMath::ACos(p[2] / pp) : -999.; // theta
    } else {
      fMomentum[0] = p[0];  // pt
      fMomentum[1] = p[1];  // phi
      fMomentum[2] = p[2];  // theta
    }
  } else {
    fMomentum[0] = -999.;
    fMomentum[1] = -999.;
    fMomentum[2] = -999.;
  }
}


//template<> void AliAODTrack::SetPosition(const double *, Bool_t);

#endif
