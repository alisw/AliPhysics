#ifndef ALIPROTONQAANALYSIS_H
#define ALIPROTONQAANALYSIS_H

/*  See cxx source for full Copyright notice */


/* $Id: AliProtonQAAnalysis.h 29114 2008-10-03 16:49:02Z pchrist $ */

//-------------------------------------------------------------------------
//                       Class AliProtonQAAnalysis
//   This is the class for the baryon (proton) analysis
//
//    Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TList.h"
#include "TArrayI.h"

class TF1;
class TH1F;
class TH3F;

class AliPID;
class AliESDEvent;
class AliESDtrack;
class AliStack;
class AliGenEventHeader;
class AliESDVertex;
class AliMCEvent;
class AliProtonAnalysisBase;

class AliProtonQAAnalysis : public TObject {
 public:
  AliProtonQAAnalysis();
  virtual ~AliProtonQAAnalysis();

  void SetBaseAnalysis(AliProtonAnalysisBase *const baseAnalysis) {
    fProtonAnalysisBase = baseAnalysis;}
  AliProtonAnalysisBase *GetProtonAnalysisBaseObject() const {
    return fProtonAnalysisBase;}

  //Vertex QA
  void RunVertexQA(AliGenEventHeader *header,
		   AliStack *stack,
		   AliESDEvent *esd);
  TList *GetVertexQAList() const {return fQAVertexList;}

  //QA histograms
  void SetQAYPtBins(Int_t nbinsY, Double_t minY, Double_t maxY,
		    Int_t nbinsPt, Double_t minPt, Double_t maxPt);
  void RunQAAnalysis(AliStack *stack, 
		     AliESDEvent *esd,
		     const AliESDVertex *vertex);
  void SetRunQAAnalysis();
  TList *GetGlobalQAList() const {return fGlobalQAList;}

  //Efficiency plots (reconstruction & PID)
  void RunReconstructionEfficiencyAnalysis(AliMCEvent *mcEvent, 
					   AliESDEvent *esd,
					   const AliESDVertex *vertex);
  void RunPIDEfficiencyAnalysis(AliStack *stack, 
				AliESDEvent *esd);
  void RunEfficiencyAnalysis(AliStack *stack, 
			     AliESDEvent *esd,
			     const AliESDVertex *vertex);
  void SetRunEfficiencyAnalysis(Bool_t gUseCuts) {
    fRunEfficiencyAnalysis = kTRUE;
    fUseCutsInEfficiency = gUseCuts;
  }
  TList *GetEfficiencyQAList() const {return fEfficiencyList;}

  //MC analysis
  void RunMCAnalysis(AliStack* stack);
  void SetRunMCAnalysis() {fRunMCAnalysis = kTRUE;}
  void SetMCProcessId(Int_t id) {
    fMCProcessIdFlag = kTRUE;
    fMCProcessId = id;
  }
  void SetMotherParticlePDGCode(Int_t pdgCode) {
    fMotherParticlePDGCodeFlag = kTRUE;
    fMotherParticlePDGCode = pdgCode;
  }
  TList *GetPDGList() const {return fPDGList;}
  TList *GetMCProcessesList() const {return fMCProcessesList;}

  TList *GetAcceptedCutList() const {return fAcceptedCutList;}
  TList *GetRejectedCutList() const {return fRejectedCutList;}
  TList *GetAcceptedDCAList() const {return fAcceptedDCAList;}
  TList *GetRejectedDCAList() const {return fRejectedDCAList;}

 private:
  AliProtonQAAnalysis(const AliProtonQAAnalysis&); // Not implemented
  AliProtonQAAnalysis& operator=(const AliProtonQAAnalysis&);// Not implemented

  void     InitVertexQA();
  void     InitQA();
  void     InitMCAnalysis();
  void     InitCutLists();
  void     InitEfficiencyAnalysis();
  void     FillQA(AliStack *stack,
		  AliESDEvent *esd,
		  const AliESDVertex *vertex, 
		  AliESDtrack* track);
 
  Bool_t   IsLabelUsed(TArrayI array, Int_t label);
  Int_t    ConvertPDGToInt(Int_t pdgCode) const;
  
  AliProtonAnalysisBase *fProtonAnalysisBase;//base analysis object

  Int_t fNBinsY; //number of bins in eta or y
  Float_t fMinY, fMaxY; //min & max value of eta or y
  Int_t fNBinsPt;  //number of bins in pT
  Float_t fMinPt, fMaxPt; //min & max value of pT
  
  //QA histograms
  //Bool_t fQAHistograms; //Boolean to activate the QA histograms
  TList *fGlobalQAList; //TList storing the directories for the QA histograms
  TList *fQAVertexList; //TList storing the vertex QA plots
  TList *fQA2DList; //TList storing the accepted primary/secondary (anti)protons
  TList *fQAPrimaryProtonsAcceptedList; //list of the QA histos for accepted primary protons
  TList *fQAPrimaryProtonsRejectedList; //list of the QA histos for rejected primary protons
  TList *fQASecondaryProtonsAcceptedList; //list of the QA histos for accepted secondary protons
  TList *fQASecondaryProtonsRejectedList; //list of the QA histos for rejected secondary protons
  TList *fQAPrimaryAntiProtonsAcceptedList; //list of the QA histos for accepted primary antiprotons
  TList *fQAPrimaryAntiProtonsRejectedList; //list of the QA histos for rejected primary antiprotons
  TList *fQASecondaryAntiProtonsAcceptedList; //list of the QA histos for accepted secondary antiprotons
  TList *fQASecondaryAntiProtonsRejectedList; //list of the QA histos for rejected secondary antiprotons

  //MC analysis
  TList *fPDGList; //list with the 3D histograms: y-pt-pdg (anti)protons
  TList *fMCProcessesList; //list with the MC processes for every secondary (anti)proton
  Bool_t fRunMCAnalysis; //run this part or not
  Bool_t fMCProcessIdFlag; //flag to see if we should check the process id
  UInt_t fMCProcessId; //process id based on the TMCProcess
  Bool_t fMotherParticlePDGCodeFlag; //flag to see if we should check the pdg code of the mother particle
  Int_t  fMotherParticlePDGCode; //pdg code of the mother particle

  TList *fAcceptedCutList;// list of the cut parameters' histograms
  TList *fRejectedCutList;// list of the cut parameters' histograms
  TList *fAcceptedDCAList;// list of the DCA histograms
  TList *fRejectedDCAList;// list of the DCA histograms

  //Efficiency (reconstruction & PID)
  Bool_t fRunEfficiencyAnalysis; //run this part or not
  Bool_t fUseCutsInEfficiency;//use the cuts in the reco and pid efficiency

  TList *fEfficiencyList;// list of the efficiency histograms

  ClassDef(AliProtonQAAnalysis,1);
};

#endif
