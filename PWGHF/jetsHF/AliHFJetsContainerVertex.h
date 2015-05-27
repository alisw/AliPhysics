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
  static const Int_t fgkContTypes = 3;
  /* enum ContType { kJets=0, kQaVtx, kJetVtx, kJetVtxData }; */
  enum ContType { kJetVtx=0, kJetVtxData, kQaVtx };
  const char *strContType(ContType e)
  {
    switch(e) {
      /* case kJets: return "kJets"; */
      case kJetVtx: return "kJetVtx";
      case kJetVtxData: return "kJetVtxData";
      case kQaVtx: return "kQaVtx";
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

  /* virtual void FillStepJets(AliHFJetsContainer::CFSteps step=kCFStepEventSelected, Double_t mult=0, const AliEmcalJet *jet=0x0, Int_t nvtx=0, Double_t partonnat[2]=0x0,Double_t partpt[2]=0x0); */
  
  virtual void FillStepQaVtx(AliHFJetsContainer::CFSteps step = kCFStepEventSelected, Double_t mult = 0, const AliEmcalJet *jet = 0x0,
                             const TClonesArray *vertices = 0x0, Double_t *disp = 0x0, Int_t nvtx = 0, const AliAODVertex *primVtx = 0x0,
                             const TClonesArray *mcPart = 0x0, Double_t p[2] = 0x0, Double_t jetpt_sub = 0.);
  
  virtual void FillStepJetVtx(AliHFJetsContainer::CFSteps step = kCFStepEventSelected, Double_t mult = 0, const AliEmcalJet *jet = 0x0,
                              const TClonesArray *vertices = 0x0, Int_t nvtx = 0, const AliAODVertex *primVtx = 0x0, const TClonesArray *mcPart = 0x0,
                              Double_t partonnat[2] = 0x0, Double_t partpt[2] = 0x0, Double_t *disp = 0x0, Double_t jetpt_sub = 0.);
  
  virtual void FillStepJetVtxData(AliHFJetsContainer::CFSteps step = kCFStepEventSelected, Double_t mult = 0, const AliEmcalJet *jet = 0x0,
                                  const TClonesArray *vertices = 0x0, Int_t nvtx = 0, const AliAODVertex *primVtx = 0x0, Double_t *disp=0x0, Double_t jetpt_sub = 0.);


 
protected:
  ContType fType;                   // container type       
  AliHFJetsTaggingVertex *fTagger;  // to use tagging methods  
  void CreateContainerVertex(ContType contType = kJetVtx); // create containers belonging to this class
  void GetBinningVertex(TString var, Int_t &nBins, Double_t *bins, const char *&axistitle); // returns array of bin limts for relevant vars
  
  ClassDef(AliHFJetsContainerVertex, 1)    // containers for HF b-jets analysis
};

#endif
