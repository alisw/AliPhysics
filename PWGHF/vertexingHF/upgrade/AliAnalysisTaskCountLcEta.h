#ifndef ALIANALYSISTASKCOUNTLCETA_H
#define ALIANALYSISTASKCOUNTLCETA_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#####################################################
//#                                                   # 
//#          Analysis Task for Lc analysis on ESD     #
//#Authors: C. Bianchin (Utrecht University)	      #
//#         and R. Romita (Univ of Liverpool,         # 
//#         Daresbury Lab),                           #
//#         based on a class                          #
//#         by MinJung Kweon, Universitaet Heidelberg #
//#                                                   #
//#         chiara.bianchin@cern.ch                   #
//#                                                   #
//#####################################################

class TParticle;
class TString;
class TList;
class TH1F;
class TLoretzVector;
class AliESDEvent;
class AliAODEvent;
class AliStack;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCountLcEta : public AliAnalysisTaskSE {
	
 public:
  AliAnalysisTaskCountLcEta(const char *name, Int_t ncuts,Double_t* cuts);
  AliAnalysisTaskCountLcEta();
  virtual ~AliAnalysisTaskCountLcEta();
		
  virtual void UserCreateOutputObjects();
	
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void FillHistosL(TParticle *part, AliStack *stack);
  Bool_t SelectTrack(TParticle *p,Bool_t fillh=kFALSE);
  Bool_t SelectTracksForCandidate(TParticle* pion, TParticle* kaon, TParticle* proton);
  Double_t InvMass(TParticle *p1p,TParticle *pn,TParticle *p2p,TLorentzVector *&candp1ppnp2p);
  void SetFillBkgHistos(Bool_t fill=kTRUE) {fFillBkg=fill;}
  Bool_t GetFillBkgHistos() const {return fFillBkg;}
  void SetDataType(TString type){fAnalysisType=type;}
  TString GetDataType() const {return fAnalysisType;}
  void SetEtaAbs(Float_t eta){fEtaAbs=eta;}
  Float_t GetEtaAbs() const {return fEtaAbs;}
  void SetEtaAbsMax(Float_t eta){fEtaAbsMax=eta;}
  Float_t GetEtaAbsMax() const {return fEtaAbsMax;}
  void SetCuts(const Int_t ncuts, Double_t* cuts){fNcuts=ncuts; fCuts=cuts;}
  Double_t* GetCuts() const {return fCuts;}
  void SetCutNames(const Int_t ncuts, TString* cutnames){if (ncuts!=fNcuts) {Printf("ERROR! %d names, expected %d",ncuts,fNcuts); return;} else fCutNames=cutnames;}
  TString* GetCutNames() const {return fCutNames;}

  void SetInvMassCut(Double_t mass){fInvMassCut=mass;}
  Double_t GetInvMassCut()const {return fInvMassCut;}

 private:
  AliAnalysisTaskCountLcEta(const AliAnalysisTaskCountLcEta &source);
  AliAnalysisTaskCountLcEta& operator=(const AliAnalysisTaskCountLcEta& source); 
 	 
  void DisableMCQA() { fEnableMCQA = kFALSE; }

  AliESDEvent 	*fESD;    	   // ESD object
  AliAODEvent 	*fAOD;    	   // AOD object
  TString 	fAnalysisType;	   // "ESD" or "AOD"	

  Long64_t fEvt;			   // event number
  TList* fOutList;             // two outputs, one for mcqa and the other for global 


  Bool_t fEnableMCQA;	// MC QA

  TH1F* fhNevt;            // histogram for book-keeping
  Float_t fEtaAbs;   // eta limit considered
  Float_t fEtaAbsMax;   // max eta limit considered for PID
  Bool_t fFillBkg;   // fill the histograms concerning background
  Int_t fNcuts; // number of selection cuts for MC candidates
  Double_t* fCuts;    //[fNcuts] cut values
  TString* fCutNames; //[fNcuts] names of the cuts variables
  Double_t fLooserPtTrack; // pt cut
  Double_t fInvMassCut;   // inv mass cut

  ClassDef(AliAnalysisTaskCountLcEta,1); // class to study Lc acceptance vs eta
};

#endif
