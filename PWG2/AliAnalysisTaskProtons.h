#ifndef AliAnalysisTaskProtons_cxx
#define AliAnalysisTaskProtons_cxx

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou
class TString;
class TList;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliProtonAnalysis;
class TF1;

#include "AliAnalysisTask.h"

class AliAnalysisTaskProtons : public AliAnalysisTask {
 public:
  enum TriggerMode { kMB1 = 0, kMB2, kSPDFASTOR }; 
  enum AnalysisMode { kInvalid = -1, kTPC = 0, kHybrid, kGlobal };
  
  AliAnalysisTaskProtons();
  AliAnalysisTaskProtons(const char *name);
  virtual ~AliAnalysisTaskProtons() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetType(const char* type) {fAnalysisType = type;}
  void SetPriorProbabilityFunctions(TF1 *felectrons,
                                    TF1 *fmuons,
                                    TF1 *fpions,
                                    TF1 *fkaons,
                                    TF1 *fprotons) {
    fFunctionUsed = kTRUE;
    fElectronFunction = felectrons;
    fMuonFunction = fmuons;
    fPionFunction = fpions;
    fKaonFunction = fkaons;
    fProtonFunction = fprotons;
  }

  void SetTriggerMode(TriggerMode triggermode) {fTriggerMode = triggermode;}
  void SetAnalysisMode(AnalysisMode analysismode) {fProtonAnalysisMode = analysismode;}
  void SetAcceptedVertexDiamond(Double_t gVx, Double_t gVy, Double_t gVz) {
    fVxMax = gVx; fVyMax = gVy; fVzMax = gVz;}

  static Bool_t IsEventTriggered(const AliESDEvent *esd, 
				 TriggerMode trigger = kMB2);
  static const  AliESDVertex *GetVertex(AliESDEvent *esd, 
					AnalysisMode mode = kTPC,
					Double_t gVx = 100.,
					Double_t gVy = 100.,
					Double_t gVz = 100.);
  
 private:
  AliESDEvent *fESD;    //ESD object 
  AliAODEvent *fAOD;    //AOD object
  AliMCEvent  *fMC;     //MC object 
  
  TString fAnalysisType;//"ESD", "AOD" or "MC"

  TList  *fList; //TList output object 
  
  AliProtonAnalysis *fProtonAnalysis; //analysis object 
  
  TF1 *fElectronFunction; //TF1 for e
  TF1 *fMuonFunction; //TF1 for mu
  TF1 *fPionFunction; //TF1 for pi
  TF1 *fKaonFunction; //TF1 for K
  TF1 *fProtonFunction; //TF1 for p
  
  Bool_t fFunctionUsed; //kTRUE if Functions are used
  
  TriggerMode   fTriggerMode; //Trigger mode
  AnalysisMode  fProtonAnalysisMode; //Analysis mode
  Double_t      fVxMax, fVyMax, fVzMax; //vertex diamond constrain

  AliAnalysisTaskProtons(const AliAnalysisTaskProtons&); // not implemented
  AliAnalysisTaskProtons& operator=(const AliAnalysisTaskProtons&); // not implemented
  
  ClassDef(AliAnalysisTaskProtons, 1); // example of analysis
};

#endif
