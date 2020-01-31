
#ifndef ALIANALYSISTASKHIGHPTDEDX_H
#define ALIANALYSISTASKHIGHPTDEDX_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
//#include <TTreeStream.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 
#include "DebugClassesMultESA2013.h"
#include "AliEventCuts.h"

 
class AliAnalysisTaskHighPtDeDx : public AliAnalysisTaskSE {
 public:
  enum AnalysisMode { kInvalid = -1, kGlobalTrk = 0x1, kTPCTrk = 0x2 }; 
  AliAnalysisTaskHighPtDeDx();
  AliAnalysisTaskHighPtDeDx(const char *name);

  virtual ~AliAnalysisTaskHighPtDeDx();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  AliEventCuts fEventCuts; /// Event cuts

  Bool_t   GetAnalysisMC() { return fAnalysisMC; }
  Bool_t   GetCentFrameworkAliCen() { return fCentFrameworkAliCen; }   
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     
  Double_t GetEtaCutStack() { return fEtaCutStack; }   
  Double_t GetMinPt() { return fMinPt; }   
  Double_t GetMinPtV0() { return fMinPtV0; }
  Double_t GetCosPACut() { return fCosPACut; }
  Double_t GetDecayRCut() { return fDecayRCut; }
  Int_t    GetContributorsVtxCut() { return fContributorsVtxCut; }
  Int_t    GetContributorsVtxSPDCut() { return fContributorsVtxSPDCut; }
  Double_t GetPileupCut() { return fPileupCut; }
  Double_t GetVtxR2Cut() { return fVtxR2Cut; }
  Double_t GetCrossedRowsCut() { return fCrossedRowsCut; }
  Double_t GetCrossedOverFindableCut() { return fCrossedOverFindableCut; }
  ULong64_t GetNegTrackStatus() { return fNegTrackStatus; }
  ULong64_t GetPosTrackStatus() { return fPosTrackStatus; }
  Float_t GetNegTOFExpTDiff() { return fNegTOFExpTDiff; }
  Float_t GetPosTOFExpTDiff() { return fPosTOFExpTDiff; }
  Bool_t   GetRejectKinks() { return fRejectKinks; }
  Bool_t   GetSigmaDedxCut() { return fSigmaDedxCut; }
  
  virtual void  SetTrigger1(UInt_t ktriggerInt1) {ftrigBit1 = ktriggerInt1;}
  virtual void  SetTrigger2(UInt_t ktriggerInt2) {ftrigBit2 = ktriggerInt2;}
  virtual void  SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
  virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
  virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
  virtual void  SetProduceVZEROBranch(Bool_t prodvzerob) {fVZEROBranch = prodvzerob;}
  virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  virtual void  SetCentDetector(const char* centDetector) {fCentDetector = centDetector;}
  virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
  virtual void  SetCentFrameworkAliCen(Bool_t isAliCen) {fCentFrameworkAliCen = isAliCen;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetEtaCutStack(Double_t etaCutStack){fEtaCutStack = etaCutStack;}
  virtual void  SetMinPt(Double_t value) {fMinPt = value;}   
  virtual void  SetMinPtV0(Double_t value) {fMinPtV0 = value;}
  virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
  virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
  virtual void  SetLowPtFraction(Double_t value) {fLowPtFraction = value;}   
  virtual void  SetMassCut(Double_t massCut){fMassCut = massCut;}
  virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) {fAnalysisPbPb = isanaPbPb;}
  virtual void  SetAnalysisRun2(Bool_t isanaRun2) {fAnalysisRun2 = isanaRun2;}
  virtual void  SetCosPACut(Double_t value) {fCosPACut = value;}   
  virtual void  SetDecayRCut(Double_t value) {fDecayRCut = value;}
  virtual void  SetContributorsVtxCut(Int_t value) {fContributorsVtxCut = value;}
  virtual void  SetContributorsVtxSPDCut(Int_t value) {fContributorsVtxSPDCut =  value;}
  virtual void  SetPileupCut(Double_t value) {fPileupCut = value;}
  virtual void  SetVtxR2Cut(Double_t value) { fVtxR2Cut = value;}
  virtual void  SetCrossedRowsCut(Double_t value) {fCrossedRowsCut = value;}
  virtual void  SetCrossedOverFindableCut(Double_t value) {fCrossedOverFindableCut = value;}
  virtual void  SetNegTrackStatus(Double_t value) {fNegTrackStatus = value;}
  virtual void  SetPosTrackStatus(Double_t value) {fPosTrackStatus = value;}
  virtual void  SetNegTOFExpTDiff(Double_t value) {fNegTOFExpTDiff = value;}
  virtual void  SetPosTOFExpTDiff(Double_t value) {fPosTOFExpTDiff = value;}
  virtual void  SetRejectKinks(Bool_t isRejectKinks) {fRejectKinks = isRejectKinks;}
  virtual void  SetSigmaDedxCut(Bool_t isSigmaDedxCut) {fSigmaDedxCut = isSigmaDedxCut;}

  
  //Task Configuration: trigger selection
  /* void SetSelectedTriggerClass1(AliVEvent::EOfflineTriggerTypes trigType) {fTrigType1 = trigType;} */
  /* void SetSelectedTriggerClass2(AliVEvent::EOfflineTriggerTypes trigType) {fTrigType2 = trigType;}  */
  /* void SetSelectedTriggerClass3(AliVEvent::EOfflineTriggerTypes trigType) {fTrigType3 = trigType;}  */

 private:
  virtual Float_t GetVertex(const AliVEvent* event) const;
  virtual void AnalyzeAOD(AliAODEvent* aod); 
  virtual void ProduceArrayTrksAOD(AliAODEvent* event, AnalysisMode anamode );
  virtual void ProduceArrayV0AOD(AliAODEvent* event, AnalysisMode anamode );
  Short_t   GetPidCode(Int_t pdgCode) const;

  void      ProcessMCTruthAOD(); 
  void      Sort(TClonesArray* array, Bool_t isMC);
  Short_t   GetPythiaEventProcessType(Int_t pythiaType);
  Short_t   GetDPMjetEventProcessType(Int_t dpmJetType);
  ULong64_t GetEventIdAsLong(AliVHeader* header) const;

  TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
  Int_t      FindPrimaryMotherLabel(AliStack* stack, Int_t label);

  AliAODMCParticle* FindPrimaryMotherAOD(AliAODMCParticle* startParticle);

  TParticle* FindPrimaryMotherV0(AliStack* stack, Int_t label);
  Int_t      FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps);

  AliAODMCParticle* FindPrimaryMotherAODV0(AliAODMCParticle* startParticle, Int_t& nSteps);

  /* AliVEvent::EOfflineTriggerTypes fTrigType1; // trigger type */
  /* AliVEvent::EOfflineTriggerTypes fTrigType2; // trigger type */
  /* AliVEvent::EOfflineTriggerTypes fTrigType3; // trigger type */
  

  static const Double_t fgkClight;   // Speed of light (cm/ps)

  AliAODEvent* fAOD;                  //! AOD object
  AliMCEvent*  fMC;                   //! MC object
  TClonesArray* fMCArray;             //! MC array for AOD
  AliAnalysisFilter* fTrackFilter;    //  Track Filter, old cuts 2010
  AliAnalysisFilter* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
  AliAnalysisFilter* fTrackFilterTPC; // track filter for TPC only tracks
  TString       fAnalysisType;        //  "ESD" or "AOD"
  TString       fCentDetector;        //  e.g. "V0M" or "V0A"
  Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
  Bool_t        fCentFrameworkAliCen; //   kTRUE: use AliCentrality, kFALSE: use AliMultSelection
  Bool_t        fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
  Bool_t        fAnalysisRun2;        //  true for LHC run-2 analyses
  Bool_t        fVZEROBranch;         //true if you want to store VZERO cells information
  TRandom*      fRandom;              //! random number generator
  DeDxEvent*    fEvent;               //! event pointer
  TClonesArray* fTrackArrayGlobalPar;          //! track array pointer, global tracks
  TClonesArray* fV0ArrayGlobalPar;             //! V0 array pointer, global tracks
  TClonesArray* fTrackArrayMC;        //! MC track array pointer
  TClonesArray* fVZEROArray;          //! array of the v0 cells.


  //
  // Cuts and options
  //
  UInt_t       ftrigBit1;
  UInt_t       ftrigBit2;
  Double_t     fVtxCut;             // Vtx cut on z position in cm
  Double_t     fEtaCut;             // Eta cut used to select particles
  Double_t     fEtaCutStack;        // Eta cut used to select particles - reduce saved stack size
  Double_t     fMinPt;              // Min pt - for histogram limits
  Double_t     fMinPtV0;            // Min pt - for histogram limits - V0s / strangeness part
  Double_t     fLowPtFraction;      // Fraction of tracks below min pt to keep
  Double_t     fCosPACut;           // Min cosPA - for histogram limits
  Double_t     fDecayRCut;          // Min decay radius
  Double_t     fMassCut;            // Reject all v0 with all dmass > masscut!
  Float_t      fMinCent; //minimum centrality
  Float_t      fMaxCent; //maximum centrality

  Int_t        fContributorsVtxCut; // Number of tracks (+1) used to fit this vertex
 
  Int_t        fContributorsVtxSPDCut; // Number of tracks (+1) used to fit this vertex
  Double_t     fPileupCut;          // |zvtx-zvtxSPD| > fPileupCut is considered bad vtx
  Double_t     fVtxR2Cut;           // r = sqrt(x^2+y^2) which is the distance between PV and the z axis
  Double_t     fCrossedRowsCut;     // CrossedRowsTOC
  Double_t     fCrossedOverFindableCut; // CrossedRowsTPC / findable 
  ULong64_t    fNegTrackStatus;     
  ULong64_t    fPosTrackStatus;
  Float_t      fNegTOFExpTDiff; 
  Float_t      fPosTOFExpTDiff; 
  Bool_t       fRejectKinks;        // reject kink daughters
  Bool_t       fSigmaDedxCut;       // dE/dx cut < 3 sigma on proton daughter candidates with momentum < 1 GeV/c:

  //
  // Help variables
  //
  Short_t      fMcProcessType;      // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
  Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only)
  Short_t      fVtxStatus;          // -1 = no vtx, 0 = outside cut, 1 = inside cut
  Float_t      fZvtx;               // z vertex
  Float_t      fZvtxMC;             // z vertex MC (truth)
  Int_t        fRun;                // run no
  ULong64_t    fEventId;            // unique event id

  Int_t        fTriggerInt;         // 0 = kMB, 1 = kCent, 2 = kSemiCent
  Int_t        fV0Finder;           // 0 = oldFinder, 1 = newFinder
  Int_t        fCentFramework;      // 0 = AliCentrality, 1 = AliMultSelection
    
  //
  // Output objects
  //
  TList*        fListOfObjects;     //! Output list of objects
  TH1I*         fEvents;            //! No of accepted events
  TH1I*         fVtx;               //! Event vertex info
  TH1I*         fVtxMC;             //! Event vertex info for ALL MC events
  TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
  TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts
  TH1F* fn1;
  TH1F* fn2;
  TH1F* fcent;
  TTree*        fTree;              //! Debug tree 

  ClassDef(AliAnalysisTaskHighPtDeDx, 1);    //Analysis task for high pt analysis 
};

#endif
