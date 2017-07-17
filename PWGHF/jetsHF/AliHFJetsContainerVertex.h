#ifndef ALIHFJETSCONTAINERVERTEX_H
#define ALIHFJETSCONTAINERVERTEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* mailto: svallero@to.infn.it */

/* Class defining containers for the HF b-jets analysis */

//--Root--
class TClonesArray;

//--AliRoot--
class AliAODVertex;
class AliEmcalJet;

//--AliHFJetsClass--
#include "AliHFJetsUtils.h"
#include "AliHFJetsContainer.h"
#include "AliHFJetsTaggingVertex.h"

//-------------------------------------------------------------------------------------

class AliHFJetsContainerVertex : public AliHFJetsContainer {
  
 public:
  
  static const Int_t fgkContTypes = 3;
  
  enum EContTYPE {
    
    kJetVtxSim  = 0,
    kJetVtxData = 1,
    kQaVtx
  };
  
  const char *fgStrContType(EContTYPE type) {
    
    switch(type) {
     
      case kJetVtxSim :
        return "kJetVtxSim";

      case kJetVtxData:
        return "kJetVtxData";

      case kQaVtx:
        return "kQaVtx";
    }
   
    return "invalid";
  }
  
  // Constructors
  AliHFJetsContainerVertex();
  
  AliHFJetsContainerVertex(const char *name, EContTYPE contType);
  
  AliHFJetsContainerVertex(const AliHFJetsContainerVertex &c);
  
  // Destructor
  virtual ~AliHFJetsContainerVertex();
  
  // Assignment operator
  AliHFJetsContainerVertex &operator=(const AliHFJetsContainerVertex &corr);
  
  virtual void FillStepJetVtxSim(CFSteps                step,
                                 const Int_t            nSVtx,
                                 Double_t               evtMult,          // Event multiplicity (not used yet)
                                 Double_t               jetPt_wBkgRej,    // Jet's pt after background subtraction
                                 vctr_pair_dbl_int     &arrVtxDisp,       // Vector of pair with SV_vtx and SV_sigma
                                 const TClonesArray    *arrVtxHF,
                                 const AliAODVertex    *primVtx,
                                 const AliEmcalJet     *jet,
                                 const TClonesArray    *mcPart,
                                 const Double_t        *partonnat,
                                 const Double_t        *partpt,
                                 Double_t               wght);
  
  virtual void FillStepJetVtxData(CFSteps               step,
                                  const Int_t           nSVtx,
                                  Double_t              evtMult,
                                  Double_t              jetPt_wBkgRej,
                                  vctr_pair_dbl_int    &arrVtxDisp,
                                  const TClonesArray   *arrVtxHF,
                                  const AliAODVertex   *primVtx,
                                  const AliEmcalJet    *jet,
                                  Double_t              wght);
  
  virtual void FillStepQaVtx(CFSteps                    step,
                             const Int_t                nSVtx,
                             Double_t                   evtMult,
                             const AliAODVertex        *primVtx,
                             const AliEmcalJet         *jet,
                             const TClonesArray        *arrVtxHF,
                             const TClonesArray        *mcPart,
                             vctr_pair_dbl_int         &arrVtxDisp,
                             Double_t                  *p,
                             Double_t                   wght);
  
protected:
  
  // Useful tools
  void Copy(TObject &c) const;
  
  void CreateContainerVertex(EContTYPE contType = kJetVtxSim);  // create containers belonging to this class
  
  void GetBinningVertex(TString     var,
                        Int_t      &nBins,
                        Double_t   *bins,
                        const char *&axistitle);                // returns array of bin limts for relevant vars
  
private:
  
  EContTYPE fType;                   // container type
  
  AliHFJetsTaggingVertex *fTagger;    // to use tagging methods

  ClassDef(AliHFJetsContainerVertex, 3)    // containers for HF b-jets analysis
};

#endif
