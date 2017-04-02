#ifndef AliNanoAODTrack_H
#define AliNanoAODTrack_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//     AOD NanoAOD track implementation
//     NanoAOD tracks are lightweight tracks which only contain a subset
//     of and ESD/AOD tracks. They are meant to be stored in specialized
//     nanoAOD for specific analysis.
//     The standard constructor takes either a AOD or an ESD track and a 
//     list of comma separated variables.
//     Only those variables are actually allocated and saved with the track.
//     Attempts to use any other variable produces an AliFatal
//    
//     Allowed kin var:
//      pt, theta, phi, chi2perNDF, posx, posy, posz, covmat, posDCAx,
//      posDCAy, pDCAx, pDCAy, pDCAz, RAtAbsorberEnd, TPCncls, TPCnclsF,
//      TPCNCrossedRows, TrackPhiOnEMCal, TrackEtaOnEMCal,
//      TrackPtOnEMCal, ITSsignal, TPCsignal, TPCsignalTuned,
//      TPCsignalN, TPCmomentum, TPCTgl, TOFsignal, integratedLenght,
//      TOFsignalTuned, HMPIDsignal, HMPIDoccupancy, TRDsignal,
//      TRDChi2, TRDnSlices

//    Custom vars
//      if the var name begins by "cst", it allows to add a track
//      property not initially foreseen in the Vtrack or in the
//      original AOD. For instance, you can add the value of the
//      bayesian probability to be a kaon by using cstKBayes. Custom
//      variables can be set at once using
//      AliNanoAODTrack::SetCustomVariables(Double_t *vars), by
//      providing them in the same order as they are defined. // FIXME: to be implemented
//
//
//     TODO
//      - Go through class again after implementation
//      - Decide if you want to create a short or int array
//      - in the constructor, set the covariant matrix only if requested
//
//     INFO/TO BE DECIDED
//      - I think this should not support muons, since muons have already their specialized AOD class
//     Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include <TRef.h>
#include <TBits.h>
#include <iostream>
#include "AliVTrack.h"
#include "AliAODVertex.h"
#include "AliAODRedCov.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "TObjString.h"
#include "TMap.h"
#include "AliNanoAODTrackMapping.h"
#include "AliNanoAODStorage.h"


#include <vector>

class AliVVertex;
class AliDetectorPID;
class AliTPCdEdxInfo;
class AliAODEvent;
class AliAODTrack;
class AliESDTrack;

class AliNanoAODTrack : public AliVTrack, public AliNanoAODStorage {

public:
  
  using TObject::ClassName;
  
  AliNanoAODTrack();
  AliNanoAODTrack(AliAODTrack * aodTrack, const char * vars);
  AliNanoAODTrack(AliESDTrack * esdTrack, const char * vars);
  AliNanoAODTrack(const char * vars);

  virtual ~AliNanoAODTrack();
  AliNanoAODTrack(const AliNanoAODTrack& trk); 
  AliNanoAODTrack& operator=(const AliNanoAODTrack& trk);


  virtual void Clear(Option_t * opt) ;
  
  // kinematics
  virtual Double_t OneOverPt() const { return (Pt() != 0.) ? 1./Pt() : -999.; }
  virtual Double_t Phi()       const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPhi());   }
  virtual Double_t Theta()     const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTheta()); }
  
  virtual Double_t Px() const { return Pt() * TMath::Cos(Phi()); }
  virtual Double_t Py() const { return Pt() * TMath::Sin(Phi()); }
  virtual Double_t Pz() const { return Pt() / TMath::Tan(Theta()); }
  virtual Double_t Pt() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPt()); }
  virtual Double_t P()  const { return TMath::Sqrt(Pt()*Pt()+Pz()*Pz()); }
  virtual Bool_t   PxPyPz(Double_t p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE; }

  virtual Double_t Xv() const { return GetProdVertex() ? GetProdVertex()->GetX() : -999.; }
  virtual Double_t Yv() const { return GetProdVertex() ? GetProdVertex()->GetY() : -999.; }
  virtual Double_t Zv() const { return GetProdVertex() ? GetProdVertex()->GetZ() : -999.; }
  virtual Bool_t   XvYvZv(Double_t x[3]) const { x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return kTRUE; }

  Double_t Chi2perNDF()  const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetChi2PerNDF()); }  
  UShort_t GetTPCNcls()  const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCncls()); } // FIXME: should this be short?

  virtual Double_t M() const { AliFatal("Not Implemented"); return -1; }
  Double_t M(AliAODTrack::AODTrkPID_t pid) const;
  virtual Double_t E() const { AliFatal("Not Implemented"); return -1; }
  Double_t E(AliAODTrack::AODTrkPID_t pid) const;
  Double_t E(Double_t m) const { return TMath::Sqrt(P()*P() + m*m); }
  virtual Double_t Y() const { AliFatal("Not Implemented"); return  -1; }
  Double_t Y(AliAODTrack::AODTrkPID_t pid) const;
  Double_t Y(Double_t m) const;
  
  virtual Double_t Eta() const { return -TMath::Log(TMath::Tan(0.5 * Theta())); }
  virtual Short_t  Charge() const {return fCharge; } // FIXME: leave like this? Create shorts array?
  virtual Double_t GetSign() const {return fCharge; }
  virtual Bool_t   PropagateToDCA(const AliVVertex *vtx, 
				  Double_t b, Double_t maxd, Double_t dz[2], Double_t covar[3]);


  // Bool_t IsOn(Int_t mask) const {return (fFlags&mask)>0;}
  ULong_t GetStatus() const { AliFatal("Not implemented"); return 0; }
  // ULong_t GetFlags() const { return fFlags; }

  //  Int_t   GetID() const { return (Int_t)fID; } // FIXME another int (short)
  Int_t   GetID() const { AliFatal("Not Implemented"); return 0; } // FIXME another int (short)
  Int_t   GetLabel() const { return fLabel; }  // 
  // void    GetTOFLabel(Int_t *p) const;


  //  Bool_t  IsPrimaryCandidate() const;
  // FIXME: do we need any of the following (note that I already removed the bit enums)
  // Bool_t  GetUsedForVtxFit() const { return TestBit(kUsedForVtxFit); }
  // Bool_t  GetUsedForPrimVtxFit() const { return TestBit(kUsedForPrimVtxFit); }

  // Bool_t  IsHybridGlobalConstrainedGlobal() const { return TestBit(kIsHybridGCG); }
  // Bool_t  IsHybridTPCConstrainedGlobal() const { return TestBit(kIsHybridTPCCG); }
  // Bool_t  IsTPCOnly() const { return IsTPCConstrained(); } // obsolete bad naming
  // Bool_t  IsTPCConstrained() const { return TestBit(kIsTPCConstrained); }
  // Bool_t  IsGlobalConstrained() const { return TestBit(kIsGlobalConstrained); }
  //
  //
  template <typename T> void GetP(T *p) const {
    p[0]=Pt(); p[1]=Phi(); p[2]=Theta();}
  //using AliVVtrack::GetP;
  //  template <typename T> void GetPxPyPz(T *p) const {
  //    p[0] = Px(); p[1] = Py(); p[2] = Pz();}
  Bool_t GetPxPyPz(Double_t *p) const;

  
  template <typename T> Bool_t GetPosition(T *x) const {
    x[0]=GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX()); x[1]=GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY()); x[2]=GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosZ());
    return TestBit(AliAODTrack::kIsDCA);}

  // FIXME: only allocate if listed?
  // template <typename T> void SetCovMatrix(const T *covMatrix) {
  //   if(!fCovMatrix) fCovMatrix=new AliAODRedCov<6>();
  //   fCovMatrix->SetCovMatrix(covMatrix);}

  // template <typename T> Bool_t GetCovMatrix(T *covMatrix) const {
  //   if(!fCovMatrix) return kFALSE;
  //   fCovMatrix->GetCovMatrix(covMatrix); return kTRUE;}

  Bool_t GetXYZ(Double_t *p) const {
    return GetPosition(p); }
  
  //using AliVVtrack::GetXYZ;
  Bool_t GetXYZAt(Double_t x, Double_t b, Double_t *r) const;
  
  Bool_t GetCovarianceXYZPxPyPz(Double_t /*cv*/[21]) const {AliFatal("Not implemented"); return 0;}
  //   return GetCovMatrix(cv);}

  // void RemoveCovMatrix() {delete fCovMatrix; fCovMatrix=NULL;}

  Double_t XAtDCA() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosDCAx()); }
  Double_t YAtDCA() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosDCAy()); }
  Double_t ZAtDCA() const { // FIXME: not sure about this one
    if (TestBit(AliAODTrack::kIsDCA)) return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY());
    else return -999.; }
  Bool_t   XYZAtDCA(Double_t x[3]) const { x[0] = XAtDCA(); x[1] = YAtDCA(); x[2] = ZAtDCA(); return kTRUE; }
  
  Double_t DCA() const { // FIXME: not sure about this one
    if (TestBit(AliAODTrack::kIsDCA)) return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX()); // FIXME: Why does this return posX?
    else return -999.; }
  
  Double_t PxAtDCA() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPDCAX()); }
  Double_t PyAtDCA() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPDCAY()); }
  Double_t PzAtDCA() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetPDCAZ()); }
  Double_t PAtDCA() const { return TMath::Sqrt(PxAtDCA()*PxAtDCA() + PyAtDCA()*PyAtDCA() + PzAtDCA()*PzAtDCA()); }
  Bool_t   PxPyPzAtDCA(Double_t p[3]) const { p[0] = PxAtDCA(); p[1] = PyAtDCA(); p[2] = PzAtDCA(); return kTRUE; }
  
  Double_t GetRAtAbsorberEnd() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetRAtAbsorberEnd()); }
  
  // For this whole block of cluster maps I could simply define a cluster map in the int array. For the moment comment all maps. Maybe not neede 
  UChar_t  GetITSClusterMap() const       { AliFatal("Not Implemented"); return 0;};
  // Int_t    GetITSNcls() const; 
  // Bool_t   HasPointOnITSLayer(Int_t i) const { return TESTBIT(GetITSClusterMap(),i); }
  // UShort_t GetHitsPatternInTrigCh() const { return (UShort_t)((fITSMuonClusterMap&0xff00)>>8); }// Fixme: array of uchars or of int?
  // UInt_t   GetMUONClusterMap() const      { return (fITSMuonClusterMap&0x3ff0000)>>16; } // 
  // UInt_t   GetITSMUONClusterMap() const   { return fITSMuonClusterMap; }
  
  // Bool_t  TestFilterBit(UInt_t filterBit) const {return (Bool_t) ((filterBit & fFilterMap) != 0);}
  // Bool_t  TestFilterMask(UInt_t filterMask) const {return (Bool_t) ((filterMask & fFilterMap) == filterMask);}
  // void    SetFilterMap(UInt_t i){fFilterMap = i;}
  // UInt_t  GetFilterMap() const {return fFilterMap;}

  // const TBits& GetTPCClusterMap() const {return fTPCClusterMap;}
  // const TBits* GetTPCClusterMapPtr() const {return &fTPCClusterMap;}
  // const TBits& GetTPCFitMap() const {return fTPCFitMap;}
  // const TBits* GetTPCFitMapPtr() const {return &fTPCFitMap;}
 
  Float_t GetTPCClusterInfo(Int_t /*nNeighbours=3*/, Int_t /*type=0*/, Int_t /*row0=0*/, Int_t /*row1=159*/, Int_t /*type*/=0) const { AliFatal("Not Implemented"); return 0;};  // FIXME What is this? 
  
  // const TBits& GetTPCSharedMap() const {return fTPCSharedMap;}
  // const TBits* GetTPCSharedMapPtr() const {return &fTPCSharedMap;}
  // void    SetTPCClusterMap(const TBits amap) {fTPCClusterMap = amap;}
  // void    SetTPCSharedMap(const TBits amap) {fTPCSharedMap = amap;}
  // void    SetTPCFitMap(const TBits amap) {fTPCFitMap = amap;}
  // 
  void    SetTPCPointsF(UShort_t  findable){fVars[AliNanoAODTrackMapping::GetInstance()->GetTPCnclsF()] = findable;}
  void    SetTPCNCrossedRows(UInt_t n)     {fVars[AliNanoAODTrackMapping::GetInstance()->GetTPCNCrossedRows()] = n;}

  UShort_t GetTPCNclsF() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCnclsF());}
  UShort_t GetTPCNCrossedRows()  const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCNCrossedRows());}
  Float_t  GetTPCFoundFraction() const { return GetTPCNCrossedRows()>0 ? float(GetTPCNcls())/GetTPCNCrossedRows() : 0;}

  // Calorimeter Cluster
  // FIXME: to be implemented
  // Int_t GetEMCALcluster() const {return fCaloIndex;}
  // void SetEMCALcluster(Int_t index) {fCaloIndex=index;}
  // Bool_t IsEMCAL() const {return fFlags&kEMCALmatch;}

  Double_t GetTrackPhiOnEMCal() const {return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTrackPhiOnEMCal());}
  Double_t GetTrackEtaOnEMCal() const {return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTrackEtaOnEMCal());}
  Double_t GetTrackPtOnEMCal() const  {return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTrackPtOnEMCal());}
  Double_t GetTrackPOnEMCal() const {return TMath::Abs(GetTrackEtaOnEMCal()) < 1 ? GetTrackPtOnEMCal()*TMath::CosH(GetTrackEtaOnEMCal()) : -999;}
  void SetTrackPhiEtaPtOnEMCal(Double_t phi,Double_t eta,Double_t pt) {fVars[AliNanoAODTrackMapping::GetInstance()->GetTrackPhiOnEMCal()]=phi;fVars[AliNanoAODTrackMapping::GetInstance()->GetTrackEtaOnEMCal()]=eta;fVars[AliNanoAODTrackMapping::GetInstance()->GetTrackPtOnEMCal()]=pt;}

  //  Int_t GetPHOScluster() const {return fCaloIndex;} // TODO: int array
  //  void SetPHOScluster(Int_t index) {fCaloIndex=index;}
  //  Bool_t IsPHOS() const {return fFlags&kPHOSmatch;}//TODO: fFlags?

  //pid signal interface
  //TODO you can remove the PID object
  Double_t  GetITSsignal()       const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetITSsignal());}
  Double_t  GetTPCsignal()       const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCsignal());}
  Double_t  GetTPCsignalTunedOnData() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCsignalTuned());}
  void      SetTPCsignalTunedOnData(Double_t signal) {fVars[AliNanoAODTrackMapping::GetInstance()->GetTPCsignalTuned()] = signal;}
  UShort_t  GetTPCsignalN()      const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCsignalN());}// FIXME: what is this?
  //  virtual AliTPCdEdxInfo* GetTPCdEdxInfo() const {return fDetPid?fDetPid->GetTPCdEdxInfo():0;} // FIXME: is this needed?
  Double_t  GetTPCmomentum()     const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCmomentum()); }
  Double_t  GetTPCTgl()          const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCTgl());      } // FIXME: what is this?
  Double_t  GetTOFsignal()       const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTOFsignal());   } 
  Double_t  GetIntegratedLength() const { AliFatal("Not implemented"); return 0;} // TODO: implement track lenght
  void      SetIntegratedLength(Double_t/* l*/) {AliFatal("Not implemented");}
  Double_t  GetTOFsignalTunedOnData() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTOFsignalTuned());}
  void      SetTOFsignalTunedOnData(Double_t signal) {fVars[AliNanoAODTrackMapping::GetInstance()->GetTOFsignalTuned()] = signal;}
  Double_t  GetHMPIDsignal()      const {return GetVar(AliNanoAODTrackMapping::GetInstance()->GetHMPIDsignal());}; 
  Double_t  GetHMPIDoccupancy()  const {return GetVar(AliNanoAODTrackMapping::GetInstance()->GetHMPIDoccupancy());}; 
  
      
  
  virtual void GetIntegratedTimes(Double_t */*times*/, Int_t) const { AliFatal("Not implemented"); return;} // FIXME check implementation

  // ________________________________________________________________________________________________________________
  // *** FIXME Those methods used fPID in the AOD track. Reimplement as an array of kin vars? ***
  // It's also a bit ugly that I used a mixed approach here: some methods are commented, some other alifatal
  // void GetHMPIDtrk(Float_t &x, Float_t &y, Float_t &th, Float_t &ph) const;  
  // void GetHMPIDmip(Float_t &x,Float_t &y,Int_t &q,Int_t &nph) const;
  //  Bool_t GetOuterHmpPxPyPz(Double_t *p) const;
  //  Int_t     GetHMPIDcluIdx()     const;// FIXME: array of ints?
  //   void      GetITSdEdxSamples(Double_t s[4]) const; // FIXME: To be reimplemented. Use one kin var for each sample
  Int_t   GetTOFBunchCrossing(Double_t /*b=0*/, Bool_t /*tpcPIDonly=kFALSE*/) const {AliFatal("Not Implemented"); return 0;};
  UChar_t   GetTRDncls(Int_t /*layer*/)                           const {AliFatal("Not Implemented"); return 0;}; 
  Double_t  GetTRDslice(Int_t /*plane*/, Int_t /*slice*/)         const {AliFatal("Not Implemented"); return 0;};
  Double_t  GetTRDmomentum(Int_t /*plane*/, Double_t */*sp*/=0x0) const {AliFatal("Not Implemented"); return 0;};
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  Double_t  GetTRDsignal()         const {return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTRDsignal());}
  Double_t  GetTRDchi2()           const {return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTRDChi2());}
  UChar_t   GetTRDncls()           const {return GetTRDncls(-1);}
  Int_t     GetNumberOfTRDslices() const { return GetVar(AliNanoAODTrackMapping::GetInstance()->GetTRDnSlices()); }

  const AliAODEvent* GetAODEvent() const {return fAODEvent;}// FIXME: change to special event type
  void SetAODEvent(const AliAODEvent* ptr){fAODEvent = ptr;}

  AliAODVertex *GetProdVertex() const { return (AliAODVertex*)fProdVertex.GetObject(); } // FIXME: change to special vertex type
  
  // print
  void  Print(const Option_t *opt = "") const;

  // setters
  // void SetFlags(ULong_t flags) { fFlags = flags; }
  // void SetStatus(ULong_t flags) { fFlags|=flags; }
  // void ResetStatus(ULong_t flags) { fFlags&=~flags; }

  //  void SetID(Short_t id) { fID = id; }
  void SetLabel(Int_t label) { fLabel = label; }
  // void SetTOFLabel(const Int_t* p);
  template <typename T> void SetPosition(const T *x, Bool_t isDCA = kFALSE);
  void SetDCA(Double_t d, Double_t z);
  void SetUsedForVtxFit(Bool_t used = kTRUE) { used ? SetBit(AliAODTrack::kUsedForVtxFit) : ResetBit(AliAODTrack::kUsedForVtxFit); }
  void SetUsedForPrimVtxFit(Bool_t used = kTRUE) { used ? SetBit(AliAODTrack::kUsedForPrimVtxFit) : ResetBit(AliAODTrack::kUsedForPrimVtxFit); }

  void SetIsTPCOnly(Bool_t b = kTRUE) { SetIsTPCConstrained(b); }// obsolete bad naming

  void SetIsTPCConstrained(Bool_t b = kTRUE) { b ? SetBit(AliAODTrack::kIsTPCConstrained) : ResetBit(AliAODTrack::kIsTPCConstrained); }
  void SetIsHybridTPCConstrainedGlobal(Bool_t hybrid = kTRUE) { hybrid ? SetBit(AliAODTrack::kIsHybridTPCCG) : ResetBit(AliAODTrack::kIsHybridTPCCG); }

  void SetIsGlobalConstrained(Bool_t b = kTRUE) { b ? SetBit(AliAODTrack::kIsGlobalConstrained) : ResetBit(AliAODTrack::kIsGlobalConstrained); }
  void SetIsHybridGlobalConstrainedGlobal(Bool_t hybrid = kTRUE) { hybrid ? SetBit(AliAODTrack::kIsHybridGCG) : ResetBit(AliAODTrack::kIsHybridGCG); }



  void SetOneOverPt(Double_t oneOverPt) { fVars[AliNanoAODTrackMapping::GetInstance()->GetPt()] = 1. / oneOverPt; }
  void SetPt(Double_t pt) { fVars[AliNanoAODTrackMapping::GetInstance()->GetPt()] = pt; };
  void SetPhi(Double_t phi) { fVars[AliNanoAODTrackMapping::GetInstance()->GetPhi()] = phi; }
  void SetTheta(Double_t theta) { fVars[AliNanoAODTrackMapping::GetInstance()->GetTheta()] = theta; }
  template <typename T> void SetP(const T *p, Bool_t cartesian = kTRUE);// TODO: WHAT IS THIS FOR?
  void SetP() {AliFatal("Not Implemented");}

  void SetXYAtDCA(Double_t x, Double_t y) {fVars[AliNanoAODTrackMapping::GetInstance()->GetPosDCAx()] = x;  fVars[AliNanoAODTrackMapping::GetInstance()->GetPosDCAy()]= y;}
  void SetPxPyPzAtDCA(Double_t pX, Double_t pY, Double_t pZ) {fVars[AliNanoAODTrackMapping::GetInstance()->GetPDCAX()] = pX; fVars[AliNanoAODTrackMapping::GetInstance()->GetPDCAY()] = pY; fVars[AliNanoAODTrackMapping::GetInstance()->GetPDCAZ()] = pZ;}
  
void SetRAtAbsorberEnd(Double_t r) { fVars[AliNanoAODTrackMapping::GetInstance()->GetRAtAbsorberEnd()] = r; }
  
  void SetCharge(Short_t q) { fCharge = q; }
void SetChi2perNDF(Double_t chi2perNDF) { fVars[AliNanoAODTrackMapping::GetInstance()->GetChi2PerNDF()] = chi2perNDF; }

  // void SetITSClusterMap(UChar_t itsClusMap)                 { fITSMuonClusterMap = (fITSMuonClusterMap&0xffffff00)|(((UInt_t)itsClusMap)&0xff); }
  // void SetHitsPatternInTrigCh(UShort_t hitsPatternInTrigCh) { fITSMuonClusterMap = (fITSMuonClusterMap&0xffff00ff)|((((UInt_t)hitsPatternInTrigCh)&0xff)<<8); }
  // void SetMuonClusterMap(UInt_t muonClusMap)                { fITSMuonClusterMap = (fITSMuonClusterMap&0xfc00ffff)|((muonClusMap&0x3ff)<<16); }
  // void SetITSMuonClusterMap(UInt_t itsMuonClusMap)          { fITSMuonClusterMap = itsMuonClusMap; }
  // void SetMUONtrigHitsMapTrg(UInt_t muonTrigHitsMap) { fMUONtrigHitsMapTrg = muonTrigHitsMap; }
  // UInt_t GetMUONTrigHitsMapTrg() { return fMUONtrigHitsMapTrg; }
  // void SetMUONtrigHitsMapTrk(UInt_t muonTrigHitsMap) { fMUONtrigHitsMapTrk = muonTrigHitsMap; }
  // UInt_t GetMUONTrigHitsMapTrk() { return fMUONtrigHitsMapTrk; }

  void     SetProdVertex(TObject *vertex) { fProdVertex = vertex; }
  // void     SetType(AliAODTrack::AODTrk_t ttype) { fType=ttype; }// FIXME: what is this?

  // 


  // Dummy: FIXME why is this dummy?
  Int_t    PdgCode() const {return 0;}

  //  needed  to inherit from VTrack, but not implemented
  virtual void  SetDetectorPID(const AliDetectorPID */*pid*/)  {AliFatal("Not Implemented"); return ;}; 
  virtual const AliDetectorPID* GetDetectorPID() const {AliFatal("Not Implemented"); return 0;}; 
  virtual UChar_t  GetTRDntrackletsPID() const  {AliFatal("Not Implemented"); return 0;}; 
  virtual void      GetHMPIDpid(Double_t */*p*/) const  {AliFatal("Not Implemented"); return;}; 
  virtual Double_t GetBz() const  {AliFatal("Not Implemented"); return 0;}; 
  virtual void     GetBxByBz(Double_t [3]/*b[3]*/) const  {AliFatal("Not Implemented"); return;}; 
  virtual const    AliExternalTrackParam * GetOuterParam() const {AliFatal("Not Implemented"); return 0;}; 
  virtual const    AliExternalTrackParam * GetInnerParam() const {AliFatal("Not Implemented"); return 0;}; 
  virtual Int_t    GetNcls(Int_t /*idet*/) const {AliFatal("Not Implemented"); return 0;}; 
  virtual const Double_t *PID() const {AliFatal("Not Implemented"); return 0;}; 



private :



  // Momentum & position
  // FIXME: the following was replaced by posx, posy, posz. Check if the names make sense
  //  Double32_t    fPosition[3];       // position of first point on track or dca
  Int_t         fLabel;             // track label, points back to MC track
  TRef          fProdVertex;        // vertex of origin
  Short_t       fCharge; // track charge
  const AliAODEvent* fAODEvent;     //! 

  ClassDef(AliNanoAODTrack, 1);
};

// inline Bool_t  AliNanoAODTrack::IsPrimaryCandidate() const
// {
//   // FIXME: reimplement
//   // True of track passes primary particle selection (independent of type) 
//   // 
//   if (fFilterMap) {
//     return kTRUE;
//   } else {
//     return kFALSE;
//   }
// }

// FIXME: is this needed?
// inline Int_t AliNanoAODTrack::GetITSNcls() const 
// {
//   // Number of points in ITS
//   Int_t n=0;
//   for(Int_t i=0;i<6;i++) if(HasPointOnITSLayer(i)) n++;
//   return n;
// }

//______________________________________________________________________________
template <typename T> 
void AliNanoAODTrack::SetPosition(const T *x, const Bool_t dca) 
{
  // set the position

  if (x) {
    if (!dca) {
      ResetBit(AliAODTrack::kIsDCA);

      fVars[AliNanoAODTrackMapping::GetInstance()->GetPosX()] = x[0];
      fVars[AliNanoAODTrackMapping::GetInstance()->GetPosY()] = x[1];
      fVars[AliNanoAODTrackMapping::GetInstance()->GetPosZ()] = x[2];
    } else {
      SetBit(AliAODTrack::kIsDCA);
      // don't know any better yet
      fVars[AliNanoAODTrackMapping::GetInstance()->GetPosX()] = -999.;
      fVars[AliNanoAODTrackMapping::GetInstance()->GetPosY()] = -999.;
      fVars[AliNanoAODTrackMapping::GetInstance()->GetPosZ()] = -999.;
    }
  } else {
    ResetBit(AliAODTrack::kIsDCA);

    fVars[AliNanoAODTrackMapping::GetInstance()->GetPosX()] = -999.;
    fVars[AliNanoAODTrackMapping::GetInstance()->GetPosY()] = -999.;
    fVars[AliNanoAODTrackMapping::GetInstance()->GetPosZ()] = -999.;
  }
}

//template<> void AliNanoAODTrack::SetPosition(const double *, Bool_t);




#endif
