#ifndef ALIANALYSISTASKHYPERONSPECTRA_H
#define ALIANALYSISTASKHYPERONSPECTRA_H





//#ifndef SIGMA0TESTSPECTRA_H
//#define SIGMA0TESTSPECTRA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to prepeare lists of photons in calorimeters,
// converted photon etc. and fill few QA histograms
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskSigma0.h"

//class TLorentzVector; class AliCaloParticle ;

class AliAnalysisTaskSigma0Spectra : public AliAnalysisTaskSigma0
{
	
 public:
  AliAnalysisTaskSigma0Spectra();
  AliAnalysisTaskSigma0Spectra(const char* name);
  virtual ~AliAnalysisTaskSigma0Spectra() ;// virtual destructor
		
  // Implementation of interface methods
  virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
		
 protected:
//   void ProcessMC();

  void FillSpectrum() ;
  void FillRealMixed() ;
  //  void FillTPCRealMixed() ;
  // void FillEventSig(AliCaloParticle * trig, const Int_t itype = 100) ;
  
  //  Bool_t ArePi0Partners(Double_t m, Double_t pt, const char * type) ;

  void FillCorr200(TLorentzVector * trig, const Int_t itype = 100 ) ;

  void GetArPod( Double_t pos[3], Double_t neg[3], Double_t moth[3],  Double_t arpod[2] );

  //  void AliAnalysisTaskSigma0::GetArmenterosQtAlfa(AliCaloParticle* positiveKFParticle, 
  //					AliCaloParticle * negativeKFParticle, AliKFParticle * gammaKFCandidate, Double_t armenterosQtAlfa[2] );


 private:
  AliAnalysisTaskSigma0Spectra(const AliAnalysisTaskSigma0Spectra&); // Not implemented
  AliAnalysisTaskSigma0Spectra& operator=(const AliAnalysisTaskSigma0Spectra&); // Not implemented
		
 protected:
   //Containers for storing previous events
  TList * fPHOSPi0Events[10] ;      //Container for PHOS photons
  TList * fEMCALPi0Events[10] ;     //Container for EMCAL photons
  TList * fConvPi0Events[10] ;      //Container for conversion photons


  //  TList * fTPCRho0Events[10] ;      //Container for PHOS photons
  //  TList * fTPCRhoplusEvents[1] ;     //Container for EMCAL photons
  //  TList * fTPCRhominusEvents[1] ;     //Container for EMCAL photons
  //  TList * fTPCPhi0Events[10] ;      //Container for conversion photons

  //Current event
  TClonesArray * fPHOSPi0Event ;   //tracks in the current event
  TClonesArray * fEMCALPi0Event ;   //tracks in the current event
  TClonesArray * fConvPi0Event ;   //tracks in the current event


  //  TClonesArray * fTPCRho0Event ;   //tracks in the current event
  //  TClonesArray * fTPCRhoplusEvent ;   //tracks in the current event
  //  TClonesArray * fTPCRhominusEvent ;   //tracks in the current event
  //  TClonesArray * fTPCPhi0Event ;   //tracks in the current event
   		
  //  Int_t fLeadingRho0TPC ;      //Leading rho0
  //  Double_t fELeadingRho0TPC ;  //Energy of the leading rho0
  //  Int_t fLeadingRhoplusTPC ;      //Leading pi0
  // Double_t fELeadingRhoplusTPC ;  //Energy of the leading pi0
  // Int_t fLeadingRhominusTPC ;      //Leading pi0
  // Double_t fELeadingRhominusTPC ;  //Energy of the leading pi0
  // Int_t fLeadingPhi0TPC ;      //Leading phi0
  //  Double_t fELeadingPhi0TPC ;  //Energy of the leading phi0


   		
  ClassDef(AliAnalysisTaskSigma0Spectra, 1); // Analysis task for conversion + calorimeters
};

#endif //ALIANALYSISTASHYPERON_H

// #endif //SIGMA0TESTSPECTRA_H
