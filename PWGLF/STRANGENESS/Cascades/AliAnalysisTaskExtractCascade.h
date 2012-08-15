/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef ALIANALYSISTASKEXTRACTCASCADE_H
#define ALIANALYSISTASKEXTRACTCASCADE_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskExtractCascade : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskExtractCascade();
  AliAnalysisTaskExtractCascade(const char *name);
  virtual ~AliAnalysisTaskExtractCascade();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;

  void SetIsNuclear           (Bool_t lIsNuclear   = kTRUE ) { fkIsNuclear   = lIsNuclear;   }
  void SetIsLowEnergyPP       (Bool_t lLowEnergyPP = kTRUE ) { fkLowEnergyPP = lLowEnergyPP; }
  
 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  TList  *fListHist;  //! List of Cascade histograms
  TTree  *fTreeCascade;              //! Output Tree, Cascades

  //Objects that have to be streamed:
  AliPIDResponse *fPIDResponse;     // PID response object
  AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition

  //Objects Controlling Task Behaviour 
  // (have to be streamed too or configuration is lost)
  
  Bool_t fkIsNuclear;   //if true, replace multiplicity est. by centrality (default FALSE) 
  Bool_t fkLowEnergyPP; //if true, skip FASTOnly (default FALSE)

//===========================================================================================
//   Variables for tree, cascades
//===========================================================================================

  Int_t fTreeCascVarCharge;         //! 
  Float_t fTreeCascVarMassAsXi;     //! 
  Float_t fTreeCascVarMassAsOmega;  //! 
  Float_t fTreeCascVarPt;           //!
  Float_t fTreeCascVarPtMC;         //!
  Float_t fTreeCascVarRapMC;        //!
  Float_t fTreeCascVarRapXi;        //!
  Float_t fTreeCascVarRapOmega;     //!
  Float_t fTreeCascVarNegEta;       //!
  Float_t fTreeCascVarPosEta;       //!
  Float_t fTreeCascVarBachEta;      //!
  Float_t fTreeCascVarDCACascDaughters; //!
  Float_t fTreeCascVarDCABachToPrimVtx; //!
  Float_t fTreeCascVarDCAV0Daughters;   //!
  Float_t fTreeCascVarDCAV0ToPrimVtx;   //!
  Float_t fTreeCascVarDCAPosToPrimVtx;  //!
  Float_t fTreeCascVarDCANegToPrimVtx;  //!
  Float_t fTreeCascVarCascCosPointingAngle; //!
  Float_t fTreeCascVarCascRadius;           //!
  Float_t fTreeCascVarV0Mass;               //!
  Float_t fTreeCascVarV0CosPointingAngle;   //!
  Float_t fTreeCascVarV0Radius;             //!
  Int_t   fTreeCascVarLeastNbrClusters;     //!
  Int_t   fTreeCascVarMultiplicity;         //!
  Float_t fTreeCascVarDistOverTotMom;       //!
  Int_t   fTreeCascVarPID;         //!
  Int_t   fTreeCascVarPIDBachelor; //!  
  Int_t   fTreeCascVarPIDNegative; //!
  Int_t   fTreeCascVarPIDPositive; //!
  Float_t fTreeCascVarPosTransMom;   //!
  Float_t fTreeCascVarNegTransMom;   //!
  Float_t fTreeCascVarPosTransMomMC; //!
  Float_t fTreeCascVarNegTransMomMC; //!

  Float_t fTreeCascVarNegNSigmaPion;   //!
  Float_t fTreeCascVarNegNSigmaProton; //!
  Float_t fTreeCascVarPosNSigmaPion;   //! 
  Float_t fTreeCascVarPosNSigmaProton; //! 
  Float_t fTreeCascVarBachNSigmaPion;  //! 
  Float_t fTreeCascVarBachNSigmaKaon;  //! 

//===========================================================================================
//   Histograms
//===========================================================================================

   TH1F      *fHistV0MultiplicityBeforeTrigSel;              //! V0 multiplicity distribution
   TH1F      *fHistV0MultiplicityForTrigEvt;                 //! V0 multiplicity distribution
   TH1F      *fHistV0MultiplicityForSelEvt;                  //! V0 multiplicity distribution
   TH1F      *fHistV0MultiplicityForSelEvtNoTPCOnly;         //! V0 multiplicity distribution
   TH1F      *fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup; //! V0 multiplicity distribution

   TH1F      *fHistMultiplicityBeforeTrigSel;     //! multiplicity distribution      
   TH1F      *fHistMultiplicityForTrigEvt;        //! multiplicity distribution
   TH1F      *fHistMultiplicity;                  //! multiplicity distribution
   TH1F      *fHistMultiplicityNoTPCOnly;         //! multiplicity distribution
   TH1F      *fHistMultiplicityNoTPCOnlyNoPileup; //! multiplicity distribution

//---> Generated Histograms

   TH1F      *fHistPVx;                      //! PVx distrib
   TH1F      *fHistPVy;                      //! PVy distrib
   TH1F      *fHistPVz;                      //! PVz distrib
   TH1F      *fHistPVxAnalysis;                      //! PVx distrib
   TH1F      *fHistPVyAnalysis;                      //! PVy distrib
   TH1F      *fHistPVzAnalysis;                      //! PVz distrib

   AliAnalysisTaskExtractCascade(const AliAnalysisTaskExtractCascade&);            // not implemented
   AliAnalysisTaskExtractCascade& operator=(const AliAnalysisTaskExtractCascade&); // not implemented
   
   ClassDef(AliAnalysisTaskExtractCascade, 11);
};

#endif
