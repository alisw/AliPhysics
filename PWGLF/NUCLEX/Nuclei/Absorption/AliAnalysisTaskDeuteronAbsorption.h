/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDeuteronAbsorption_H
#define AliAnalysisTaskDeuteronAbsorption_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"
#include "AliPID.h"
#include <string>

class AliPIDResponse;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliESDtrackCuts;
class AliMultSelection;

#define kNabsSpecies 5

class AliAnalysisTaskDeuteronAbsorption : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskDeuteronAbsorption(const char *name = "AliAnalysisTaskDeuteronAbsorption");
  virtual ~AliAnalysisTaskDeuteronAbsorption();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option) {}

  void SetMindEdx(double opt) { fMindEdx = opt; }
  void SetTreeFlag(Bool_t tmode) {fTreemode = tmode;};
  void SetMinTPCsignalN(double signalN = 50) { fMinTPCsignalN = signalN; }
  void SetParticleType(AliPID::EParticleType type = AliPID::kHe3) {ParticleType = type; };
  void SetESDtrackCuts(const AliESDtrackCuts& cuts) { fESDtrackCuts = cuts; }
  void SetMaxNsig(double opt = 6.) { fMaxNSigma = opt; }
  void SetMinNsig(double opt = -6.) { fMinNSigma = opt; }
    
  double GetMindEdx() const { return fMindEdx; }
  double GetMinTPCsignalN() const { return fMinTPCsignalN; }
  
  static const AliPID::EParticleType fgkSpecies[kNabsSpecies];
  static const std::string fgkParticleNames[kNabsSpecies];
  static const double fgkPhiParamPos[4][4];
  static const double fgkPhiParamNeg[4][4];

  bool fUseTRDboundariesCut;
  float fNtpcSigmas;
  AliEventCuts fEventCuts;
  bool fUseTrackCuts;
  bool fUseTOFallClustersInfo;
  
private:
  double fMindEdx = 100.0; /// Cut on the minimum dE/dx in TPC
  int    fMinTPCsignalN = 50; /// Minimum number of PID clusters in the TPC
  Bool_t fTreemode = kFALSE;    // Flag for filling the tree mode
  AliPID::EParticleType ParticleType = AliPID::kHe3;    // to select He3 or triton
  double fMaxNSigma = 6.;
  double fMinNSigma = -6.;
    
  AliPIDResponse *fPIDResponse;   //! pid response
  AliESDtrackCuts fESDtrackCuts;  // input track cuts
                                  //
  TList *fOutputList;             //! output list

  TTree *fTreeTrack;  //! tree for some track parameters
  
  // Variables for the tree
  //Double_t tP;
  Float_t  tCentrality;       // centrality                                                                                      //!
  Int_t    tNglobalTracks;   // number of good tracks per event
  Float_t  tPt;              // pt of the track (at inner wall of the TPC)
  Float_t  tEta;             // eta of the track (at inner wall of the TPC)
  Float_t  tPhi;             // phi of the track (at inner wall of the TPC)
  Float_t  tSign;            // 
  Float_t  tdEdx;            // 
  Float_t  tdEdxExp;         // dEdx expected for 3He (or t)
  Float_t  tdEdxExpSigma;    // sigma expected for 3He (or t)
  Float_t  tnsigTPC;         // nSigma PID to 3He in the TPC
  Float_t  tnsigTOF;         // nSigma PID to 3He in the TOF
  Float_t  tmass2;           // m^2/z^2 of the track based on the TOF
  Float_t  tITSchi2;         // ITS chi2
  Float_t  tTOFsigDx;        // track-to-hit residual in TOF (x-direction)
  Float_t  tTOFsigDz;        // track-to-hit residual in TOF (z-direction)
  Float_t  tTOFchi2;         // chi2 of the hit in the TOF
  Float_t  tTPCchi2;         // chi2
  Float_t  tTPCxRows;        // TPC crossed rows
  Float_t  tTPCxRowsOverFindable; // TPC crossed rows over findable
  Float_t  tDCAxy;           // DCAxy
  Float_t  tDCAz;            // DCAz
  Float_t  tT0res;           // T0 resolution
  Int_t    tT0mask;          // T0 mask
  Int_t    tTRDclsN;         // number of TRD clusters attached to the track
  UChar_t  tTRDntracklets;   // number of TRD tracklets used for tracking
  UChar_t  tTRDNchamberdEdx; // number of chambers used to calculate the TRD truncated mean
  Int_t    tID;              // identification number of the track
  Int_t    tPdgCodeMc;       // pdg code of the track if MC information is available
  Int_t    tTOFchannel;
  UChar_t  tnPIDclsTPC;      // number of clusters used for PID in the TPC
  UChar_t  tITSclsMap;       // ITS cluster map
  Float_t  tMCpt;            // MC pt
  Float_t  tMCabsMom;        // MC absorption momentum
  Float_t  tMCabsRadius;     // MC absorption radius
  Float_t  tMCabsEta;
  Float_t  tMCabsPhi;
  Bool_t   tMCtofMismatch;   //
  Int_t    tNdaughters; // 
  Bool_t   tIsReconstructed; // False for MC particles 
  Bool_t   thasTOF;          // 
  Int_t    tRunNumber;       //
  UChar_t  tPIDforTracking;  //

  UChar_t  tTOFclsN;         // number of cluster candidates in TOF
  Float_t  tNsigmaTOFarray[64], tmass2array[64], tdxTOFarray[64], tdzTOFarray[64], tTOFarray[64], tLengtharray[64], texpTOFarray[64], tsigmaTOFarray[64];
  Int_t    tTOFchannelarray[64];
  Bool_t   tMCtofMismatcharray[64];
  UChar_t  tNmatchableTracks[64];
  //
  TH1F *fHistZv;      //! Primary vertex z distribution
  TH3F *fHist3TPCpid[kNabsSpecies];  //! QA TPC dE/dx per species
  TH3F *fHist3TPCpidAll;  //! QA TPC dE/dx no species selection
  TH3F *fHist3TOFpid[kNabsSpecies];  //! QA TOF beta per species
  TH3F *fHist3TOFpidAll;  //! QA TOF beta no species selection
  TH3F *fHist3TOFmass[kNabsSpecies]; //! QA TOF mass per species
  TH3F *fHist3TOFnsigma[kNabsSpecies]; //! QA TOF nsigma per species
  TH3F *fHist3TOFmassAll; //! QA TOF mass no species selection

  void FillTOFallMatchableHitsInfo(AliESDtrack *track);

  //
  AliAnalysisTaskDeuteronAbsorption(const AliAnalysisTaskDeuteronAbsorption &);            // not implemented
  AliAnalysisTaskDeuteronAbsorption &operator=(const AliAnalysisTaskDeuteronAbsorption &); // not implemented
  
  ClassDef(AliAnalysisTaskDeuteronAbsorption, 4);
};

#endif
