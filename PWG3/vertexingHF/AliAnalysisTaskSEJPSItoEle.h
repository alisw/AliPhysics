#ifndef ALIANALYSISTASKSEJPSITOELE_H
#define ALIANALYSISTASKSEJPSITOELE_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
//                     Class AliAnalysisTaskSEJPSItoEle
//                     AliAnalysisTaskSE class to select 
//              J/psi -> e+e- candidates only and save them 
//                    into a specific stand-alone AOD file 
//              Author: C.Di Giglio, carmelo.digiglio@ba.infn.it
//*************************************************************************

#include "TClonesArray.h"
#include "TH1F.h" 
#include "TH2F.h"
#include "TList.h"
#include "TChain.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSEJPSItoEle : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEJPSItoEle();
  AliAnalysisTaskSEJPSItoEle(const char *name);
  virtual ~AliAnalysisTaskSEJPSItoEle();

  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetCutsJPSI(const Double_t cutsJPSI[9]);
  void SetPtCuts(const Double_t ptcuts[2]);
  void SetAODMCInfo(Bool_t okMCInfo) { fOkAODMC = okMCInfo;}
  void ReadAODMCInfo(AliAODEvent* aodEv, const TClonesArray* inArray);

 private:

  AliAnalysisTaskSEJPSItoEle(const AliAnalysisTaskSEJPSItoEle &source);
  AliAnalysisTaskSEJPSItoEle& operator=(const AliAnalysisTaskSEJPSItoEle& source); 
  
  TList *fOutput;                            //! list send on output slot 0

  //output histograms
  TH1F *fhDecayTimeMCjpsifromB;              //! Pseudo-proper decay time distribution used as template for JPSIs from B
  TH1F *fhDecayTime;                         //! Pseudo-proper decay time distribution
  TH1F *fhDecayTimeOut;                      //! Pseudo-proper decay time distribution (stand-alone AOD)
  TH1F *fhInvMass;                           //! Invariant mass distribution
  TH1F *fhD0;                                //! Impact parameter distribution
  TH1F *fhD0D0;                              //! Product of impact parameters distributions
  TH1F *fhCosThetaStar;                      //! Cosine of decay angle distribution
  TH1F *fhCosThetaPointing;                  //! Cosine of pointing angle distribution
  TH1F *fhDCA;                               //! Distance of closest approach
  TH2F *fhCtsVsD0D0;                         //! Cos theta star Vs. D0D0 distribution

  //like sign pairs histograms 
  TH1F *fHistMassLS;                         //! Invariant mass distribution
  TH1F *fHistCtsLS;                          //! Cosine of decay angle distribution
  TH1F *fHistCtsLSpos;                       //! Cosine of decay angle distribution (++ pairs)
  TH1F *fHistCtsLSneg;                       //! Cosine of decay angle distribution (-- pairs)
  TH1F *fHistCPtaLS;                         //! Cosine of pointing angle distribution
  TH1F *fHistd0d0LS;                         //! Product of impact parameters distributions
  TH1F *fHistDCALS;                          //! Distance of closest approach
  AliAnalysisVertexingHF *fVHF;              //! Vertexer heavy flavour (used to pass the cuts)

  //Int_t totJPSIin;
  //Int_t totJPSIout;

  //like sign spectrum normalization
  Int_t fTotPosPairs;                        //
  Int_t fTotNegPairs;                        //
  Double_t fLsNormalization;                 // Like sign normalization factor

  //flags for analysis
  Bool_t fOkAODMC;                           // Flag to read AOD monte carlo information
  Bool_t fOkLikeSign;                        // Flag to select like sign candidates analysis

  //momentum cuts   
  Double_t fCuts[9];                         // cuts for N-tuple values selection
  Double_t fPtCuts[2];                       // Max and min pt of the candidates

  TClonesArray    *fVerticesHFTClArr;        // Array of heavy flavor vertices to be replicated in stand-alone AOD
  TClonesArray    *fJpsiToEleTClArr;         // Array of J/psi->e+e- candidates to be replicated in stand-alone AOD
  TClonesArray    *fLikeSignTClArr;          // Array of like sign candidates to be replicated in stand-alone AOD
  TClonesArray    *fTracksTClArr;            // Array of tracks belonging to J/psi->e+e- candidates to be replicated in stand-alone AOD 
  TChain          *fChain; 
  AliAODEvent     *fOrigAOD;                 // original AOD event
  AliAODEvent     *fNewAOD;                  // new AOD event with only JPSItoEle candidates stored

  ClassDef(AliAnalysisTaskSEJPSItoEle,1); // AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates
};
#endif
