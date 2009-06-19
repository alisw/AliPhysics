#ifndef ALIANALYSISTASKCHARMFRACTION_H
#define ALIANALYSISTASKCHARMFRACTION_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskCharmFraction
// AliAnalysisTask for the extraction of the fraction of prompt charm
// using the charm hadron impact parameter to the primary vertex
// Author: Andrea Rossi andrea.rossi@ts.infn.it
//*************************************************************************

class TH1F;
class TH2F;
class AliAODDEvent;
class AliAODMCHeader;

#include "AliAnalysisTask.h"

class AliAnalysisTaskCharmFraction : public AliAnalysisTask {
 public:
  AliAnalysisTaskCharmFraction(const char *name="AliAnalysisTaskCharmFraction");
  AliAnalysisTaskCharmFraction(const char *name,Int_t nptbins);
  
  virtual ~AliAnalysisTaskCharmFraction() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void SetNPtBins(Int_t nbins=10){fnbins=nbins;}
  void SetCheckMC(Bool_t checkMC){fcheckMC=checkMC;}
  void SetCheckMC_D0(Bool_t check_D0){fcheckMCD0=check_D0;}
  void SetCheckMC_2prongs(Bool_t check2prongs){fcheckMC2prongs=check2prongs;}
  void SetCheckMC_prompt(Bool_t checkprompt){fcheckMCprompt=checkprompt;}
  void SetCheckMC_fromB(Bool_t checkfromB){fcheckMCfromB=checkfromB;}
  void SetCheckMC_fromDstar(Bool_t skipD0star){fSkipD0star=skipD0star;}
  void SetUseCuts(Bool_t usecuts){fD0usecuts=usecuts;}
  void SetSideBands(Double_t sideband){fSideBands=sideband;}
  void SetStudyPureBackground(Bool_t back){fStudyPureBackground=back;}
 private:
  AliAODEvent *fAOD;    //AOD object
  AliAnalysisVertexingHF *fVHF;        // Vertexer heavy flavour
  TClonesArray *fArrayD0toKpi;        // Array of D0->Kpi
  TClonesArray *fArrayMC;             // Array of MC info
  AliAODMCHeader *fAODmcHeader;        //AOD header
  TH1F *fhMass;                         //!Inv Mass
  TH1F *fhMassTrue;                     //!MC inv Mass
  TH2F *fhCPtaVSd0d0;                  //! histo of CosPtAngle Vs d0xd0
  TH1F *fhd0xd0;                         //! histo of d0d0
  TH1F *fhCPta;                        //! histo of cosptangle
  TH1F *fhSecVtxZ;                     //! histo of z coord of sec. vtx 
  TH1F *fhSecVtxX;                      //! histo of x coord of sec. vtx
  TH1F *fhSecVtxY;                      //! histo of x coord of sec. vtx
  TH2F *fhSecVtxXY;                      //! histo of x coord of sec. vtx
  TH1F *fhSecVtxPhi;                      //! histo of x coord of sec. vtx
  TH1F        *fhd0D0;                 //!D0 impact par histo all bins
  TH1F        **fhd0D0pt;                 //!D0 impact par histos per pt bin 
  TH1F *fhd0D0VtxTrue;                    //! histo for D0 impact par w.r.t. true vertex, all bins integrated
  TH1F      **fhd0D0VtxTruept;            //! histos for D0 impact par w.r.t. true vertex
  TH1F   *fhMCd0D0;                      //!D0 MC impact par histo all bins
  TH1F      **fhMCd0D0pt;              //!D0 MC impact par histos per pt bin 
  Int_t        fnbins;                //Number of pt bins
  Bool_t       fD0usecuts;            // Switch in the use of the cuts
  Bool_t       fcheckMC;              //  Switch on MC check: minimum check is same mother
  Bool_t       fcheckMCD0;           //  check the mother is a D0
  Bool_t       fcheckMC2prongs;         //  check the decay is in two prongs
  Bool_t       fcheckMCprompt;       //  check the D0 comes from a c quark
  Bool_t       fcheckMCfromB;        //  check the D0 comes from a b quark
  Bool_t       fSkipD0star;           // skip if D0 comes from a D*  
  Bool_t  fStudyd0fromBTrue;         // Flag for analyze true impact par of D0 from B
  Bool_t  fStudyPureBackground;      // Flag to study the background (reverse the selection on the signal)
  Double_t  fSideBands;                //Side bands selection (see cxx)
  AliAnalysisTaskCharmFraction(const AliAnalysisTaskCharmFraction&); // not implemented
  AliAnalysisTaskCharmFraction& operator=(const AliAnalysisTaskCharmFraction&); // not implemented
  
  ClassDef(AliAnalysisTaskCharmFraction,1); // analysis task for prompt charm fraction
};

#endif
