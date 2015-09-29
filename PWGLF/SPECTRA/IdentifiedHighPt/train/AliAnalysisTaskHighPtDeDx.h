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
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 
#include <AliESDtrackCuts.h>
#include "DebugClassesMultESA2013.h"



class AliAnalysisTaskHighPtDeDx : public AliAnalysisTaskSE {
 public:
  enum AnalysisMode { kInvalid = -1, kGlobalTrk = 0x1, kTPCTrk = 0x2 }; 
  AliAnalysisTaskHighPtDeDx();
  AliAnalysisTaskHighPtDeDx(const char *name);

  virtual ~AliAnalysisTaskHighPtDeDx();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     
  Double_t GetEtaCutStack() { return fEtaCutStack; }   
  Double_t GetMinPt() { return fMinPt; }   
  Double_t GetMinPtV0() { return fMinPtV0; }
  Int_t    GetTreeOption() { return fTreeOption; }

  virtual void  SetTrigger1(UInt_t ktriggerInt1) {ftrigBit1 = ktriggerInt1;}
  virtual void  SetTrigger2(UInt_t ktriggerInt2) {ftrigBit2 = ktriggerInt2;}
  virtual void  SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
  virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
  virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
  virtual void  SetProduceVZEROBranch(Bool_t prodvzerob) {fVZEROBranch = prodvzerob;}
  virtual void  SetProduceTPCBranch(Bool_t prodtpcb) {fTPCBranch = prodtpcb;}
  virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetEtaCutStack(Double_t etaCutStack){fEtaCutStack = etaCutStack;}
  virtual void  SetMinPt(Double_t value) {fMinPt = value;}   
  virtual void  SetMinPtV0(Double_t value) {fMinPtV0 = value;}
  virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
  virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
  virtual void  SetLowPtFraction(Double_t value) {fLowPtFraction = value;}   
  virtual void  SetMassCut(Double_t massCut){fMassCut = massCut;}
  virtual void  SetTreeOption(Int_t value) {fTreeOption = value;}    
  virtual void  SetRequireRecV0(Bool_t value) {fRequireRecV0 = value;}
  virtual void  SetStoreMcIn(Bool_t value) {fStoreMcIn = value;}
  virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) {fAnalysisPbPb = isanaPbPb;}

 private:
  virtual Float_t GetVertex(const AliVEvent* event) const;
  virtual void AnalyzeESD(AliESDEvent* esd); 
  virtual void AnalyzeAOD(AliAODEvent* aod); 
  virtual void ProduceArrayTrksESD(AliESDEvent* event, AnalysisMode anamode );
  virtual void ProduceArrayV0ESD(AliESDEvent* event, AnalysisMode anamode );
  virtual void ProduceArrayTrksAOD(AliAODEvent* event, AnalysisMode anamode );
  virtual void ProduceArrayV0AOD(AliAODEvent* event, AnalysisMode anamode );
  Short_t   GetPidCode(Int_t pdgCode) const;
  /* Float_t   GetSpherocity(AliESDEvent* event, AliAnalysisFilter* cuts, Float_t etacut, Float_t ptcut, Bool_t useTPCtrack); */
  /* Float_t   GetSphericity(AliESDEvent* event, AliAnalysisFilter* cuts, Float_t etacut, Float_t ptcut, Bool_t useTPCtrack); */
  /* Float_t   GetSpherocityTrue(AliStack *Stack, Float_t etacut, Float_t ptcut); */
  /* Float_t   GetSphericityTrue(AliStack *Stack, Float_t etacut, Float_t ptcut); */

  void      ProcessMCTruthESD();
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



  static const Double_t fgkClight;   // Speed of light (cm/ps)

  AliESDEvent* fESD;                  //! ESD object
  AliAODEvent* fAOD;                  //! AOD object
  AliMCEvent*  fMC;                   //! MC object
  AliStack*    fMCStack;              //! MC ESD stack
  TClonesArray* fMCArray;             //! MC array for AOD
  AliAnalysisFilter* fTrackFilter;    //  Track Filter, old cuts 2010
  AliAnalysisFilter* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
  AliAnalysisFilter* fTrackFilterTPC; // track filter for TPC only tracks
  TString       fAnalysisType;        //  "ESD" or "AOD"
  Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
  Bool_t        fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
  Bool_t        fVZEROBranch;         //true if you want to store VZERO cells information
  Bool_t        fTPCBranch;           //tru if you want to produce the TPC branch
  TRandom*      fRandom;              //! random number generator
  DeDxEvent*    fEvent;               //! event pointer
  TClonesArray* fTrackArrayGlobalPar;          //! track array pointer, global tracks
  TClonesArray* fTrackArrayTPCPar;          //! track array pointer, tpc track parameters
  TClonesArray* fV0ArrayGlobalPar;             //! V0 array pointer, global tracks
  TClonesArray* fV0ArrayTPCPar;             //! V0 array pointer, tpc tracks
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
  Double_t     fMassCut;            // Reject all v0 with all dmass > masscut!
  Int_t        fTreeOption;         // 0: no tree, >0: enable debug tree
  Float_t      fMinCent; //minimum centrality
  Float_t      fMaxCent; //maximum centrality
  Bool_t       fRequireRecV0;       // Require a v0 before updating tree
                                    // For a spectra analysis we will need to
                                    // keep track also of the empty events
  Bool_t       fStoreMcIn;          // Store MC input tracks
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
