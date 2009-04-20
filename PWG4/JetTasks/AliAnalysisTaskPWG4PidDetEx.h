#ifndef ALIANALYSISTASKPWG4PidDetEx_H
#define ALIANALYSISTASKPWG4PidDetEx_H

#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisFilter.h"

class AliAnalysisTaskPWG4PidDetEx : public AliAnalysisTaskSE {
 public:
  enum TriggerMode {kMB1, kMB2, kSPDFASTOR}; 

  AliAnalysisTaskPWG4PidDetEx();
  AliAnalysisTaskPWG4PidDetEx(const char *name);
  virtual ~AliAnalysisTaskPWG4PidDetEx();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *); 

  virtual void  SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
  void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void  SetTriggerMode(TriggerMode triggermode) {fTriggerMode = triggermode;}
  void  SetMC(Bool_t analysisMC){fAnalysisMC = analysisMC;}
  void  SetPtCut(Double_t ptCut){fPtCut = ptCut;}
  void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}

  Bool_t IsEventTriggered(AliVEvent* ev, TriggerMode trigger);

  void AnalyzeESD(AliESDEvent* esd);
  void AnalyzeAOD(AliAODEvent* aod);

 private:
  AliAnalysisTaskPWG4PidDetEx(const AliAnalysisTaskPWG4PidDetEx&);
  AliAnalysisTaskPWG4PidDetEx& operator=(const AliAnalysisTaskPWG4PidDetEx&);

  Double_t IntegratedLength(AliVTrack* track) const;
  Double_t MassSquared     (AliVTrack* track) const;
  Double_t KaonDecay       (AliVTrack* track) const;

  static const Double_t fgkCtau;                //  distance for kaon decay
  static const Double_t fgkPionMass;            //  pion mass
  static const Double_t fgkKaonMass;            //  kaon mass
  static const Double_t fgkProtonMass;          //  proton mass


  AliESDEvent* fESD;                //! ESD object
  AliAODEvent* fAOD;                //! AOD object
  TList*       fListOfHists;        //! Output list of histograms
  Double_t     fEtaCut;             //  Eta cut used to select particles
  Double_t     fPtCut;              //  pT cut used to select particles
  Int_t        fXbins;              //  #bins for Pt histos range
  Double_t     fXmin;               //  min X value for histo range
  Double_t     fTOFCutP;            //  max X value for histo range; also the p cut used in TOF for PID

  AliAnalysisFilter* fTrackFilter;  //  Track Filter
  Bool_t       fAnalysisMC;         //  Flag for MC analysis
  TString      fAnalysisType;       //  "ESD" or "AOD"
  TriggerMode  fTriggerMode;        //  Trigger mode


  //Histograms
  TH1I*         fEvents;             //! #analyzed events           
  TH1F*         fEffTot;             //! pT for all charged particles
  TH1F*         fEffPID;             //! pT for charged particles with TOF signal
  TH1F*         fAccP;               //! pT for charged particles with p < fTOFCutP
  TH1F*         fAccPt;              //! pT for charged particles with pT < fTOFCutP
  TProfile*     fKaonDecayCorr;      //! decay correction for Kaons
  TH1F*         fdNdPt;              //! pT dist (Rec)
  TH2F*         fMassAll;            //! mass calculated from TOF vs p
  TH1F*         fdNdPtPion;          //! pT for pions identified with TOF
  TH2F*         fMassPion;           //! mass for pions identified with TOF
  TH2F*         fdEdxTPCPion;        //! dE/dx vs p (TPC) for pions identified with TOF
  TH2F*         fbgTPCPion;          //! dE/dx vs betagamma (TPC) for pions identified with TOF
  TH1F*         fdNdPtKaon;          //! pT for kaons identified with TOF
  TH2F*         fMassKaon;           //! mass for kaons identified with TOF
  TH2F*         fdEdxTPCKaon;        //! dE/dx vs p (TPC) for kaons identified with TOF
  TH2F*         fbgTPCKaon;          //! dE/dx vs betagamma (TPC) for kaons identified with TOF
  TH1F*         fdNdPtProton;        //! pT for protons identified with TOF
  TH2F*         fMassProton;         //! mass for protons identified with TOF
  TH2F*         fdEdxTPCProton;      //! dE/dx vs p (TPC) for protons identified with TOF
  TH2F*         fbgTPCProton;        //! dE/dx vs betagamma (TPC) for protons identified with TOF
  TH1F*         fdNdPtMC;            //! pT dist (MC)
  TH1F*         fdNdPtMCPion;        //! pT dist for pions (MC)
  TH1F*         fdNdPtMCKaon;        //! pT dist for kaons (MC)
  TH1F*         fdNdPtMCProton;      //! pT dist for protons (MC)

  ClassDef(AliAnalysisTaskPWG4PidDetEx, 1);    //Analysis task for PWG4 PID using detector signals 
};

#endif
