#ifndef AliAnalysisTaskEfficiency_cxx
#define AliAnalysisTaskEfficiency_cxx
 
class TH1F;
class TH2F;
class TList;
class TProfile;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEfficiency : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEfficiency(const char *name = "AliAnalysisTaskEfficiency");
  virtual ~AliAnalysisTaskEfficiency() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t SelectJouri(AliESDtrack* track);
  virtual void   SetTrackType(Int_t type)      {fTrackType = type;} 
  virtual void   SetCuts(AliESDtrackCuts* cuts){fCuts = cuts;}
  virtual void   SetFieldOn(Bool_t b = kTRUE)  {fFieldOn = b;} 
  
 private:

  // For crude division into A/C side and positive/negative
  enum { kNegA = 0, kPosA, kNegC, kPosC, kTotalAC};
  enum { kTPC = 0, kSPD, kTotalVtx};
  Int_t GetIndexAC(AliESDtrack *track);
  Int_t ConvertePDG(Int_t pdg);

  Bool_t      fFieldOn;

  TList       *fHists;          // List of histos
  TH1F        *fHistRECpt;      // pt of reconstructed tracks
  TH1F        *fHistMCpt;       // pt of MC primaries
  TH1F        *fHistMCNRpt;     // pt of not reconstructed primaries
  TH1F        *fHistFAKEpt;     // pt of fake track
  TH1F        *fHistMULTpt;     // pt of multiply reconstructed tracks

  TH1F        *fHistRECeta;     // eta of reconstructed tracks
  TH1F        *fHistMCeta;      // eta of MC primaries
  TH1F        *fHistMCNReta;    // eta of not reconstructed primaries
  TH1F        *fHistFAKEeta;    // eta of fake track
  TH1F        *fHistMULTeta;    // eta of multiply reconstructed tracks

  TH1F        *fHistRECphi;     // phi of reconstructed tracks
  TH1F        *fHistMCphi;      // phi of MC primaries
  TH1F        *fHistMCNRphi;    // phi of not reconstructed primaries
  TH1F        *fHistFAKEphi;    // phi of fake track
  TH1F        *fHistMULTphi;    // phi of multiply reconstructed tracks
  
  TH1F        *fHistRECHPTeta;  // eta of reconstructed high-pT tracks
  TH1F        *fHistMCHPTeta;   // eta of MC primary high-pT tracks 
  TH1F        *fHistMCNRHPTeta; // eta of not reconstructed primaries
  TH1F        *fHistFAKEHPTeta; // eta of fake high-pT tracks
  
  TH1F        *fHistRECHPTphi;  // phi of reconstructed high-pT tracks
  TH1F        *fHistMCHPTphi;   // phi of MC primary high-pT tracks 
  TH1F        *fHistMCNRHPTphi; // phi of not reconstructed primaries
  TH1F        *fHistFAKEHPTphi; // phi of fake high-pT tracks

  TH1F        *fHistRecMult;    // Reconstruction multiplicity
  TH1F        *fHistNCluster;   // Number of clusters of suspicious tracks

  // Checks on ESD data only
  TH1F        *fh1VtxEff;              // Vtx Reconstruction eff for TPC and SPD
  TH2F        *fh2PhiPadRow[kTotalAC]; // TPC pad row vs phi
  TH2F        *fh2PhiLayer[kTotalAC];  // ITS layer vs phi
  TH2F        *fh2MultSpdChips[2];     // Multiplicity vs fired chips
  TH2F        *fh2VtxTpcSpd;           // Vtx correlation TPC <-> SPD

  // Correlation histos 
  TH2F* fh2VertexCorrelation[kTotalVtx];                    // ESD z-vtx vs MC z-vtx
  TH2F* fh2VertexCorrelationShift[kTotalVtx];               // (MC z-vtx - ESD z-vtx) vs MC z-vtx
  TH1F* fh1VertexShift[kTotalVtx];                          // (MC z-vtx - ESD z-vtx) 
  TH1F* fh1VertexShiftNorm[kTotalVtx];                      // (MC z-vtx - ESD z-vtx) / (sigma_ESD-z-vtx) histogrammed

  TH2F* fh2EtaCorrelation;                       // ESD eta vs MC eta
  TH2F* fh2EtaCorrelationShift;                  // (MC eta - ESD eta) vs MC eta
  TH2F* fh2PhiCorrelation;                       // ESD phi vs MC phi
  TH2F* fh2PhiCorrelationShift;                  // (MC phi - ESD phi) vs MC phi
  TH2F* fh2PtCorrelation;                        // ESD pt vs MC pt
  TH2F* fh2PtCorrelationShift;                   // (MC pt - ESD pt) vs MC pt

 //Tracked Deltaphi histos
  TH1F* fHistDeltaphiprimaries;                // Tracklet dPhi for primaries
  TH1F* fHistDeltaphisecondaries;              // Tracklet dPhi for secondaries
  TH1F* fHistDeltaphireject;                   // Tracklet dPhi for rejected tracks

  //ITS-TPC matching
  TH1F* fHistTPCITSdeltax;                     // TPC-ITS matching Delta_x
  TH1F* fHistTPCITSdeltay;                     // TPC-ITS matching Delta_y
  TH1F* fHistTPCITSdeltaz;                     // TPC-ITS matching Delta_z
  TH1F* fHistTPCITSdeltar;                     // TPC-ITS matching Delta_r
  //TRD-TPC matching
  TH1F* fHistTPCTRDdeltax;                     // TPC-TRD matching Delta_x
  TH1F* fHistTPCTRDdeltay;                     // TPC-TRD matching Delta_y
  TH1F* fHistTPCTRDdeltaz;                     // TPC-TRD matching Delta_z
  TH1F* fHistTPCTRDdeltar;                     // TPC-TRD matching Delta_r
  TH1F* fHistTPCTRDdeltaphi;                   // TPC-TRD matching Delta_phi
  // Pulls
  TProfile* fPtPullsPos;                          // Pt pulls
  TProfile* fPtPullsNeg;                          // Pt pulls
  TProfile* fPtShiftPos;                          // Pt shift
  TProfile* fPtShiftNeg;                          // Pt shift
  // Others
  Int_t fTrackType;
  AliESDtrackCuts* fCuts;                         // List of cuts
  // SPD Vertex
  TH2F* fVertexRvsZ;
  TH2F* fVertexRvsZC;
  // Multi Fluctuations
  TProfile* fEtaMultiFluc;                    // multiplicity fluctuations
  TH1F*     fEtaMulti;                        // multiplicity fluctuations
  TH1F*     fEtaMultiH;                       // multiplicity fluctuations
  TProfile* fPhiMultiFluc;                    // multiplicity fluctuations
  TH1F*     fPhiMulti;                        // multiplicity fluctuations
  TH1F*     fPhiMultiH;                       // multiplicity fluctuations

  //for each nuclei type one histo (d, t, 3He, 4He)
  TH1F        *fHistRECptCharge[8];       // pt of reconstructed tracks
  TH1F        *fHistMCptCharge[8];        // pt of MC primaries
  TH1F        *fHistMCNRptCharge[8];      // pt of not reconstructed primaries
  TH1F        *fHistFAKEptCharge[8];      // pt of fake track

  TH1F        *fHistRECetaCharge[8];      // eta of reconstructed tracks
  TH1F        *fHistMCetaCharge[8];       // eta of MC primaries
  TH1F        *fHistMCNRetaCharge[8];     // eta of not reconstructed primaries
  TH1F        *fHistFAKEetaCharge[8];     // eta of fake track

  TH1F        *fHistRECphiCharge[8];      // phi of reconstructed tracks
  TH1F        *fHistMCphiCharge[8];       // phi of MC primaries
  TH1F        *fHistMCNRphiCharge[8];     // phi of not reconstructed primaries
  TH1F        *fHistFAKEphiCharge[8];     // phi of fake track


  AliAnalysisTaskEfficiency(const AliAnalysisTaskEfficiency&); // not implemented
  AliAnalysisTaskEfficiency& operator=(const AliAnalysisTaskEfficiency&); // not implemented
  
  ClassDef(AliAnalysisTaskEfficiency, 1); // example of analysis
};

#endif
