#ifndef AliHFJetsContainerVertex_H
#define AliHFJetsContainerVertex_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* mailto: svallero@to.infn.it */

/* Class defining containers for the HF b-jets analysis */

#include "AliHFJetsContainer.h"
#include "AliHFJetsTaggingVertex.h"
#include "AliEmcalJet.h"
#include "TClonesArray.h"

class AliCFContainer;
class AliAODJet;

class AliHFJetsContainerVertex : public AliHFJetsContainer
{
 public:
  // 4 types of containers: 1) reco jets 2) B jets 3)primary vertex and 4) secondary vertices in a jet
  static const Int_t fgkContTypes=4;
  enum ContType { kJets=0, kBJets, kQaVtx, kJetVtx, kJetVtxData };
  const char *strContType(ContType e)
  {
    switch(e) {
      case kJets: return "kJets";
      case kBJets: return "kBJets";
      case kQaVtx: return "kQaVtx";
      case kJetVtx: return "kJetVtx";
      case kJetVtxData: return "kJetVtxData";
    }
    return "invalid";
  }
  // Constructors
  AliHFJetsContainerVertex();
  AliHFJetsContainerVertex(const char* name,ContType contType);
  AliHFJetsContainerVertex(const AliHFJetsContainerVertex &c);
  // Destructor
  virtual ~AliHFJetsContainerVertex();
  // Assignment operator
  AliHFJetsContainerVertex& operator=(const AliHFJetsContainerVertex& corr);
  // Useful tools
  virtual void Copy(TObject& c) const;

  // Methods imported from task AliAnalysisTaskSEHFJets (A. Rossi)
  virtual void FillStepJets(AliHFJetsContainer::CFSteps step=kCFStepAll, Double_t mult=0, const AliEmcalJet *jet=0x0,Double_t p[3]=0x0,Double_t contr=0,Double_t pt[3]=0x0);

  virtual void FillStepQaVtx(AliHFJetsContainer::CFSteps step=kCFStepAll, Double_t mult=0, const AliEmcalJet *jet=0x0, const TClonesArray *vertices=0x0, Double_t* disp=0x0,Int_t nvtx=0,const AliAODVertex *primVtx=0x0,const TClonesArray *mcPart=0x0,Double_t p[3]=0x0);

  virtual void FillStepJetVtx(AliHFJetsContainer::CFSteps step=kCFStepAll, Double_t mult=0, const AliEmcalJet *jet=0x0, const TClonesArray *vertices=0x0, Int_t nvtx=0,const AliAODVertex *primVtx=0x0,const TClonesArray *mcPart=0x0,Double_t p[3]=0x0,Double_t* disp=0x0);

  virtual void FillStepBJets(AliHFJetsContainer::CFSteps step=kCFStepAll, Double_t mult=0, const AliEmcalJet *jet=0x0,Int_t nvtx=0,Double_t partonnat[3]=0x0, Double_t contribution=0.,Double_t ptpart=0.);

  virtual void FillStepJetVtxData(AliHFJetsContainer::CFSteps step=kCFStepAll, Double_t mult=0, const AliEmcalJet *jet=0x0, const TClonesArray *vertices=0x0, Int_t nvtx=0,const AliAODVertex *primVtx=0x0,Double_t* disp=0x0);


  //AliHFJetsContainer* GetContainer() { return fContainer; }
  //void SetContainer(AliHFJetsContainer* hist) { fContainer = hist; }

protected:
  ContType fType;                   // container type       
  AliHFJetsTaggingVertex *fTagger;  // to use tagging methods       
  //AliCFContainer* fContainer;    // custom container 
  void CreateContainerVertex(ContType contType=kJets); // create containers belonging to this class
  void GetBinningVertex(TString var, Int_t& nBins,Double_t * bins, const char*& axistitle); // returns array of bin limts for relevant vars
  
  ClassDef(AliHFJetsContainerVertex, 1)    // containers for HF b-jets analysis
};

#endif
