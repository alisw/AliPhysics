/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Debug tree to look at the distribution of the variable we are cutting on
//
//
#ifndef ALIHFEREDUCEDTRACK_H
#define ALIHFEREDUCEDTRACK_H

#include <TObject.h>
#include <TMath.h>
#include <TBits.h>

class AliHFEreducedTrack : public TObject{
 public:
  typedef enum{
    kV0electron = 0,
    kV0proton = 1,
    kV0pion = 2,
    kV0undef = 3
  } EV0PID_t;
  AliHFEreducedTrack();
  AliHFEreducedTrack(const AliHFEreducedTrack &ref);
  AliHFEreducedTrack &operator=(const AliHFEreducedTrack &ref);
  ~AliHFEreducedTrack() {}
  
  // -------------- Getters ------------------------
  Double_t Pt() const { return TMath::Abs(fSignedPt); }
  Double_t P() const { return fP; }
  Double_t TPCmomentum() const { return fTPCmomentum; }
  Double_t Eta() const { return fEta; }
  Double_t Phi() const { return fPhi; }
  Int_t Charge() const {
    if(fSignedPt >= 0.) return 1.;
    return -1;
  }
  
  Double_t MCPt() const { return TMath::Abs(fMCSignedPt); }
  Double_t MCP() const { return fMCP; }
  Double_t MCEta() const { return fMCEta; }
  Double_t MCPhi() const { return fMCPhi; }
  Int_t MCCharge() const {
    if(fMCSignedPt >= 0.) return 1.;
    return -1;
  }
  Int_t MCPDG() const { return fMCPDG; }
  Int_t MCMotherPdg() const { return fMCMotherPdg; }
  Bool_t MCSignal() const { return fMCSignal; }
  Int_t MCSource() const { return fMCSource; }
  Double_t MCProdVtxX() const { return fMCProdVtx[0]; }
  Double_t MCProdVtxY() const { return fMCProdVtx[1]; }
  Double_t MCProdVtxZ() const { return fMCProdVtx[2]; }
  Double_t MCProdRadius() const { return TMath::Sqrt(fMCProdVtx[0]*fMCProdVtx[0]+fMCProdVtx[1]*fMCProdVtx[1]); }
  
  Double_t DCAr() const { return fDCA[0]; }
  Double_t DCAz() const { return fDCA[1]; }
  
  Bool_t TestFilterBit(Int_t ibit) const { return fFilterBit.TestBitNumber(ibit); }
  Int_t GetTrackID() const { return fTrackID; }
  Bool_t HasITSrefit() const { return fTrackStatus.TestBitNumber(kITSrefit); }
  Bool_t HasTPCrefit() const { return fTrackStatus.TestBitNumber(kTPCrefit); }
  Bool_t HasTOFpid() const { return fTrackStatus.TestBitNumber(kTOFpid); }
  Bool_t IsTOFmismatch() const { return fTrackStatus.TestBitNumber(kTOFmismatch); }
  Bool_t HasEMCALpid() const { return fTrackStatus.TestBitNumber(kEMCALpid); }
  Bool_t IsDoubleCounted() const { return fTrackStatus.TestBitNumber(kDoubleCounted); }
  
  Int_t GetITSnclusters() const { return static_cast<Int_t>(fNclustersITS); }
  Int_t GetTPCnclusters() const { return static_cast<Int_t>(fNclustersTPC); }
  Int_t GetTRDnclusters() const { return static_cast<Int_t>(fNclustersTRD); }
  Bool_t HasITScluster(int nly) const { 
    if(nly > 5) return kFALSE;
    return fITSclusterMap.TestBitNumber(nly);
  }
  Bool_t GetITSlayerStatus(int nly) const { 
    if(nly > 5) return kFALSE;
    return fITSstatusMap.TestBitNumber(nly);
  }
  Int_t GetTPCnclusterPID() const { return static_cast<Int_t>(fNclustersTPCPID); }
  Int_t GetTPCnclustersAll() const { return static_cast<Int_t>(fNclustersTPCAll); }
  Int_t GetTPCcrossedRows() const { return static_cast<Int_t>(fTPCcrossedRows); }
  Int_t GetTPCsharedClusters() const { return static_cast<Int_t>(fTPCsharedClusters); }
  Float_t GetTPCclusterRatio() const { return fTPCclusterRatio; }
  Float_t GetTPCclusterRatioAll() const { return fTPCclusterRatioAll; }
  Int_t GetTRDntrackletsPID() const { return static_cast<UChar_t>(fTRDtrackletsPID); }
  Int_t GetTRDnslices() const { return static_cast<Int_t>(fTRDnslices); }
  Bool_t GetTRDstatus(Int_t layer) const { 
    if(layer > 5) return kFALSE;
    return fTRDlayer.TestBitNumber(layer); 
  }
  Float_t GetTRDchi2() const { return fTRDchi2; }
  
  Double_t GetTPCdEdx() const { return fTPCdEdx; }
  Double_t GetTPCdEdxCorrected() const { return fTPCdEdxCorrected; }
  Double_t GetTPCsigmaEl() const { return fTPCsigmaEl; }
  Double_t GetTPCsigmaElCorrected() const { return fTPCsigmaElCorrected; }
  Double_t GetTOFsigmaEl() const { return fTOFsigmaEl; }
  Float_t GetTOFmismatchProb() const { return fTOFmismatchProb; }
  Double_t GetEMCALEoverP() const { return fEoverP; }
  Double_t GetEMCALSigmaEl() const { return fEMCALsigmaEl; }
  void GetEMCALShowerShape(Double_t showershape[4]) const{
    for(Int_t is = 0; is < 4; is++) showershape[is] = fShowerShape[is]; 
  }
  
  Bool_t IsV0electron() const { return fV0PID == kV0electron; }
  Bool_t IsV0pion() const { return fV0PID == kV0pion; }
  Bool_t IsV0proton() const { return fV0PID == kV0proton; }
  
  
  // -------------- Setters ------------------------
  void SetSignedPt(Double_t abspt, Bool_t positivecharge) { 
    Double_t charge = positivecharge ? 1. : -1;
    fSignedPt = abspt * charge;
  }
  void SetP(Double_t p) { fP = p; }
  void SetTPCmomentum(Double_t p) { fTPCmomentum = p; }
  void SetEta(Double_t eta) { fEta = eta; }
  void SetPhi(Double_t phi) { fPhi = phi; }
  
  void SetFilterBit(Int_t ibit) { fFilterBit.SetBitNumber(ibit); }
  void SetTrackID(Int_t trackID) { fTrackID = trackID; }
  void SetITSrefit() { fTrackStatus.SetBitNumber(kITSrefit); }
  void SetTPCrefit() { fTrackStatus.SetBitNumber(kTPCrefit); }
  void SetTOFpid() { fTrackStatus.SetBitNumber(kTOFpid); }
  void SetTOFmismatch() { fTrackStatus.SetBitNumber(kTOFmismatch); }
  void SetEMCALpid() { fTrackStatus.SetBitNumber(kEMCALpid); }
  void SetDoubleCounted() { fTrackStatus.SetBitNumber(kDoubleCounted); }
  
  void SetMCSignedPt(Double_t abspt, Bool_t positivecharge){
    Double_t charge = positivecharge ? 1. : -1;
    fSignedPt = abspt * charge;
  }
  void SetMCP(Double_t mcp) { fMCP = mcp; }
  void SetMCEta(Double_t mceta) { fMCEta = mceta; }
  void SetMCPhi(Double_t mcphi) { fMCPhi = mcphi; }
  void SetMCPDG(Int_t mcpdg) {fMCPDG = mcpdg; }
  void SetMCMotherPdg(Int_t pdg) { fMCMotherPdg = pdg; }
  void SetMCSignal() { fMCSignal = kTRUE; }
  void SetMCSource(Int_t mcsource) { fMCSource = mcsource; }
  void SetMCProdVtx(Double_t vx, Double_t vy, Double_t vz){
    fMCProdVtx[0] = vx;
    fMCProdVtx[1] = vy;
    fMCProdVtx[2] = vz;
  }
  
  void SetDCA(Float_t dcaR, Float_t dcaZ){
    fDCA[0] = dcaR;
    fDCA[1] = dcaZ;
  }
  
  void SetITSnclusters(int ncls) { fNclustersITS = ncls; }
  void SetTPCnclusters(int ncls) { fNclustersTPC = ncls; }
  void SetTRDnclusters(int ncls) { fNclustersTRD = ncls; }
  void SetITScluster(UInt_t ly){
    if(ly > 5) return;
    fITSclusterMap.SetBitNumber(ly); 
  }
  void SetITSstatus(UInt_t ly){
    if(ly > 5) return;
    fITSstatusMap.SetBitNumber(ly); 
  }
  void SetTPCnclustersPID(Int_t ncls) { fNclustersTPCPID = static_cast<UChar_t>(ncls); }
  void SetTPCnclustersAll(Int_t ncls) { fNclustersTPCAll = static_cast<UChar_t>(ncls); }
  void SetTPCcrossedRows(int nrows) { fTPCcrossedRows = static_cast<UChar_t>(nrows); }
  void SetTPCsharedClusters(Int_t ncls) { fTPCsharedClusters = static_cast<UChar_t>(ncls); }
  void SetTPCclusterRatio(Float_t ratio) { fTPCclusterRatio = ratio; }
  void SetTPCclusterRatioAll(Float_t ratio) { fTPCclusterRatioAll = ratio; }
  void SetTRDntrackletsPID(Int_t ntracklets) { fTRDtrackletsPID = static_cast<UChar_t>(ntracklets); }
  void SetTRDnslices(Int_t nslices) { fTRDnslices = static_cast<UChar_t>(nslices); }
  void SetTRDstatus(Int_t layer) { 
    if(layer > 5) return;
    fTRDlayer.SetBitNumber(static_cast<UInt_t>(layer)); 
  }
  void SetTRDchi2(Float_t chi2) { fTRDchi2 = chi2; }
  
  void SetTPCdEdx(Double_t dEdx) { fTPCdEdx = dEdx; }
  void SetTPCdEdxCorrected(Double_t dEdx) { fTPCdEdxCorrected = dEdx; }
  void SetTPCsigmaEl(Double_t sigma) { fTPCsigmaEl = sigma; }
  void SetTPCsigmaElCorrected(Double_t sigma) { fTPCsigmaElCorrected = sigma; }
  void SetTOFsigmaEl(Double_t sigma) { fTOFsigmaEl = sigma; }
  void SetTOFmismatchProbability(Float_t mismatchProb) { fTOFmismatchProb = mismatchProb; }
  void SetEMCALEoverP(Double_t eop) { fEoverP = eop; }
  void SetEMCALSigmaEl(Double_t sigma) { fEMCALsigmaEl = sigma; }
  void SetEMCALShowerShape(Double_t showershape[4]){
    for(Int_t is = 0; is < 4; is++) fShowerShape[is] = showershape[is]; 
  }
  void SetV0PID(AliHFEreducedTrack::EV0PID_t v0pid) { fV0PID = v0pid; }
  
 private:
  typedef enum{
    kITSrefit = 0,
    kTPCrefit = 1,
    kTOFpid   = 2,
    kTOFmismatch = 3,
    kEMCALpid =4,
    kDoubleCounted = 5,
    kKink = 6
  } TrackStatus_t;
  Double_t fSignedPt;                     // signed pt
  Double_t fP;                            // p
  Double_t fEta;                          // eta
  Double_t fPhi;                          // phi
  Double_t fTPCmomentum;                  // TPC p
  TBits    fFilterBit;                    // filterbit
  Int_t    fTrackID;                      // trackID
  Double_t fMCSignedPt;                   // MCSignedPt
  Double_t fMCP;                          // MCP
  Double_t fMCEta;                        // MCEta
  Double_t fMCPhi;                        // MCPhi
  Int_t    fMCPDG;                        // MCPDG
  Int_t    fMCMotherPdg;                  // MCMP
  Bool_t   fMCSignal;                     // MCSignal
  Int_t    fMCSource;                     // MCSource
  Double_t fMCProdVtx[3];                 // MC prod Vtx
  TBits    fTrackStatus;                  // Track Status
  UChar_t  fNclustersITS;                 // ITS nb cls
  UChar_t  fNclustersTPC;                 // TPC nb cls
  UChar_t  fNclustersTRD;                 // TRD nb cls
  TBits    fITSclusterMap;                // ITS maps
  TBits    fITSstatusMap;                 // ITS status map
  UChar_t  fNclustersTPCPID;              // TPC PID nb cls
  UChar_t  fNclustersTPCAll;              // TPC all nb cls
  UChar_t  fTPCcrossedRows;               // TPC crossed rows
  UChar_t  fTPCsharedClusters;            // TPC shared clusters
  Float_t  fTPCclusterRatio;              // TPC cls ratio
  Float_t  fTPCclusterRatioAll;           // TPC cls ratio all
  UChar_t  fTRDtrackletsPID;              // TRD tracklet PID
  UChar_t  fTRDnslices;                   // TRD nslices
  TBits    fTRDlayer;                     // TRD layer
  Float_t  fTRDchi2;                      // TRD chi2
  Double_t fTPCdEdx;                      // TPC dedx
  Double_t fTPCdEdxCorrected;             // TPC dedx corrected
  Double_t fTPCsigmaEl;                   // TPC sigma el
  Double_t fTPCsigmaElCorrected;          // TPC sigma el corrected
  Double_t fTOFsigmaEl;                   // TOF sigma el
  Float_t  fTOFmismatchProb;              // TOF mismatch prob
  Double_t fEoverP;                       // Eoverp
  Double_t fEMCALsigmaEl;                 // EMCAl sigmal el
  Double_t fShowerShape[4];               // showershape
  Float_t fDCA[2];                        // dca
  EV0PID_t fV0PID;                        // V0pid
  
  ClassDef(AliHFEreducedTrack, 1)
};
#endif
