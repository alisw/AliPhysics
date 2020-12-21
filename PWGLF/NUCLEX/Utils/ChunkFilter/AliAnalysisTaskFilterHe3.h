#ifndef AliAnalysisTaskFilterHe3_H
#define AliAnalysisTaskFilterHe3_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "AliEventCuts.h"

class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class AliPIDResponse;
class AliMultSelection;

class AliAnalysisTaskFilterHe3 : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskFilterHe3();
  AliAnalysisTaskFilterHe3(const char *name);
  virtual ~AliAnalysisTaskFilterHe3();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetUseMultTaskCentrality(Bool_t useMultTaskCentrality = kTRUE) { fUseMultTaskCentrality = useMultTaskCentrality; };
  void SetParticleType(AliPID::EParticleType particleType = AliPID::kHe3) {ParticleType = particleType; };
  
  // Setters for some cuts
  void SetMinNsig(Double_t opt) { fMinNSigma = opt; }
  void SetMaxNsig(Double_t opt) { fMaxNSigma = opt; }
  void SetNclsTPC(Double_t opt) { fMinNclsTPC = opt; }
  void SetPtotPosMin(Double_t opt) { fMinPtotPos = opt; }
  void SetPtotNegMin(Double_t opt) { fMinPtotNeg = opt; }
  void SetPtotPosMax(Double_t opt) { fMaxPtotPos = opt; }
  void SetPtotNegMax(Double_t opt) { fMaxPtotNeg = opt; }
  void SetFillSecifTOF(Bool_t opt) { fillSecifTOF = opt; }
  void SetMinMass(Double_t opt) { fMinMass = opt; }
  void SetMaxMass(Double_t opt) { fMaxMass = opt; }
  void SetMaxDCAxy(Double_t opt) { fMaxDCAxy = opt; }
  void SetMaxDCAz(Double_t opt) { fMaxDCAz = opt; }
  void SetManualBetheBlockParam2018(Double_t opt) { fUseManualBetheBlockParam2018 = opt; }
    
private:
  AliESDEvent *fESD;                                    //! input event
  TList *fOutputList;                                   //! output list
  TTree *fListOfFiles;                                  //! list of files
  AliESDtrackCuts *fESDtrackCuts;                       //! input track cuts (secondary & primary)
  AliESDtrackCuts *fESDtrackCutsPrimary;                //! input track cuts (only primary)
  AliPIDResponse *fPIDResponse;                         //! pid response object
  AliMultSelection *fMultSel;                           //! centrality and multiplicity selection
  AliEventCuts fEventCuts;                              //!
  Bool_t fUseMultTaskCentrality;                        // flag: true if centrality should be taken from AliMultSelectionTask, false if traditional method should be used
  AliPID::EParticleType ParticleType = AliPID::kHe3;    // to select He3 or triton
  //
  Int_t fEventIdFile;                                                                                       //! event id in file
  TString fFileName;                                                                                        //! chunk file name
  Float_t fCentrality;                                                                                      //!
  Bool_t fIsEventAccepted;                                                                                  //!
  Int_t fNfilteredParticles;                                                                                //!
  Float_t fSign[8192], fPtot[8192], fNsigmaTPC[8192], fdEdx[8192], fDCAxy[8192], fDCAz[8192], fMass[8192];  //!
  Int_t fNTPCclusters[8192];                                                                                      //!
  Bool_t fhasTOF[8192];                                                                                     //! 

  // Cut params
  Double_t fMinNSigma = -4.0;
  Double_t fMaxNSigma = 10.0;
  Double_t fMinNclsTPC = 50;
  Double_t fMinPtotPos = 1.5;
  Double_t fMinPtotNeg = 0.5;
  Double_t fMaxPtotPos = 20.0;
  Double_t fMaxPtotNeg = 20.0;
  Bool_t fillSecifTOF = kFALSE;
  Double_t fMinMass = 1.0;            
  Double_t fMaxMass = 2.3;            
  Double_t fMaxDCAxy = 0.5;
  Double_t fMaxDCAz = 2.0;
  Bool_t fUseManualBetheBlockParam2018 = kTRUE;
  //
  // histograms
  //
  TH1F *fHistZv;                   //! z-Vertex distribution
  TH2F *fHistdEdxData;             //! PID histogram dEdx all particles
  TProfile *fHistdEdxExpDeuteron;  //! PID-QA
  TProfile *fHistdEdxExpHe3;       //! PID-QA
  TProfile *fHistdEdxExpTriton;    //! PID-QA
  TH2F *fHistTof;                  //! PID histogram TOF all particles
  TH3F *fHistdEdxDeuteronParam[2]; //! PID-QA histogram for deuteron Bethe-Bloch parameterisation
  TH3F *fHistdEdxHe3Param[2];      //! PID-QA histogram for he3 Bethe-Bloch parameterisation
  TH3F *fHistdEdxTritonParam[2];   //! PID-QA histogram for triton Bethe-Bloch parameterisation
  TH1F *fHistCent;                 //! centrality histogram
  //
  //
  // Bethe-Bloch parameterisations (hard-coded for the time being)
  //
  Double_t fParamDeuteron[5];
  Double_t fParamTriton[5];
  Double_t fParamHe3[5];
  Double_t fParamDeuteronMC[5];
  Double_t fParamTritonMC[5];
  Double_t fParamHe3MC[5];
  //
  //
  AliAnalysisTaskFilterHe3(const AliAnalysisTaskFilterHe3 &);
  AliAnalysisTaskFilterHe3 &operator=(const AliAnalysisTaskFilterHe3 &);
  //
  ClassDef(AliAnalysisTaskFilterHe3, 1);
};

#endif
