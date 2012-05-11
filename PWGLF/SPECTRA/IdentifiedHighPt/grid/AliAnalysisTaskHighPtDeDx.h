#ifndef ALIANALYSISTASKHIGHPTDEDX_H
#define ALIANALYSISTASKHIGHPTDEDX_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TTreeStream.h>
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

#include "DebugClasses.C"



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
  Double_t GetMinPt() { return fMinPt; }   
  Int_t    GetTreeOption() { return fTreeOption; }  

  virtual void  SetTrigger1(UInt_t ktriggerInt1) {ftrigBit1 = ktriggerInt1;}
  virtual void  SetTrigger2(UInt_t ktriggerInt2) {ftrigBit2 = ktriggerInt2;}
  virtual void  SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
  virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
  virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
  virtual void  SetProduceTPCBranch(Bool_t prodtpcb) {fTPCBranch = prodtpcb;}
  virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetMinPt(Double_t value) {fMinPt = value;}   
  virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
  virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
  virtual void  SetLowPtFraction(Double_t value) {fLowPtFraction = value;}   
  virtual void  SetTreeOption(Int_t value) {fTreeOption = value;}    
  virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) {fAnalysisPbPb = isanaPbPb;}

 private:
  virtual Float_t GetVertex(const AliVEvent* event) const;
  virtual void AnalyzeESD(AliESDEvent* esd); 
  virtual void AnalyzeAOD(AliAODEvent* aod); 
  virtual void ProduceArrayTrksESD(AliESDEvent* event, AnalysisMode anamode );
  virtual void ProduceArrayTrksAOD(AliAODEvent* event, AnalysisMode anamode );
  Short_t   GetPidCode(Int_t pdgCode) const;
  void      ProcessMCTruthESD();
  void      ProcessMCTruthAOD(); 
  void      Sort(TClonesArray* array, Bool_t isMC);
  Short_t   GetPythiaEventProcessType(Int_t pythiaType);
  Short_t   GetDPMjetEventProcessType(Int_t dpmJetType);
  ULong64_t GetEventIdAsLong(AliVHeader* header) const;

  TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
  Int_t      FindPrimaryMotherLabel(AliStack* stack, Int_t label);

  AliAODMCParticle* FindPrimaryMotherAOD(AliAODMCParticle* startParticle);

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
  Bool_t        fTPCBranch;           //tru if you want to produce the TPC branch
  TRandom*      fRandom;              //! random number generator
  DeDxEvent*    fEvent;               //! event pointer
  TClonesArray* fTrackArrayGlobalPar;          //! track array pointer, global tracks
  TClonesArray* fTrackArrayTPCPar;          //! track array pointer, tpc track parameters
  TClonesArray* fTrackArrayMC;        //! MC track array pointer
  TClonesArray* fVZEROArray;          //! array of the v0 cells.


  //
  // Cuts and options
  //
  UInt_t       ftrigBit1;
  UInt_t       ftrigBit2;
  Double_t     fVtxCut;             // Vtx cut on z position in cm
  Double_t     fEtaCut;             // Eta cut used to select particles
  Double_t     fMinPt;              // Min pt - for histogram limits
  Double_t     fLowPtFraction;      // Fraction of tracks below min pt to keep
  Int_t        fTreeOption;         // 0: no tree, >0: enable debug tree
  Float_t      fMinCent; //minimum centrality
  Float_t      fMaxCent; //maximum centrality
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
