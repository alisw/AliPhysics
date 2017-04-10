#ifndef AliAnalysisTaskCountLcEta_H
#define AliAnalysisTaskCountLcEta_H
//#####################################################
//#                                                   # 
//#          Analysis Task for Lc analysis on ESD     #
//#Authors: C. Bianchin (Utrecht University)	      #
//#         and R. Romita (Univ of Liverpool,         # 
//#         Daresbury Lab),                           #
//#         based on a class                          #
//#         by MinJung Kweon, Universitaet Heidelberg #
//#                                                   #
//#####################################################

class TString;
class TLoretzVector;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCountLcEta : public AliAnalysisTaskSE {
   
public:
   AliAnalysisTaskCountLcEta(const char *name, Int_t ncuts,Double_t* cuts);
   AliAnalysisTaskCountLcEta();
   virtual ~AliAnalysisTaskCountLcEta(){};
   
   virtual void UserCreateOutputObjects();
   
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);
   void SetFillBkgHistos(Bool_t fill=kTRUE) {fFillBkg=fill;}
   Bool_t GetFillBkgHistos() const {return fFillBkg;}
   void SetDataType(TString type){fAnalysisType=type;}
   TString GetDataType() const {return fAnalysisType;}
   void SetEtaAbs(Float_t eta){fEtaAbs=eta;}
   Float_t GetEtaAbs() const {return fEtaAbs;}
   void SetEtaAbsMax(Float_t eta){fEtaAbsMax=eta;}
   Float_t GetEtaAbsMax() const {return fEtaAbsMax;}
   void SetCuts(Int_t ncuts, Double_t* cuts){fNcuts=ncuts; fCuts=cuts;}
   Double_t* GetCuts() const {return fCuts;}
   void SetCutNames(Int_t ncuts, TString* cutnames){if (ncuts!=fNcuts) {Printf("ERROR! %d names, expected %d",ncuts,fNcuts); return;} else fCutNames=cutnames;}
   TString* GetCutNames() const {return fCutNames;}
   
   void SetInvMassCut(Double_t mass){fInvMassCut=mass;}
   Double_t GetInvMassCut()const {return fInvMassCut;}
   
private:
   AliAnalysisTaskCountLcEta(const AliAnalysisTaskCountLcEta &source);
   AliAnalysisTaskCountLcEta& operator=(const AliAnalysisTaskCountLcEta& source); 
   
   
   void FillHistosL(TParticle *part, AliMCEvent* mcEvent);
   Bool_t SelectTrack(TParticle *p,Bool_t fillh=kFALSE);
   Bool_t SelectTracksForCandidate(TParticle* pion, TParticle* kaon, TParticle* proton);
   Double_t InvMass(TParticle *p1p,TParticle *pn,TParticle *p2p,TLorentzVector *&candp1ppnp2p);
   void FillHistogramsBackgroundCandidates(TParticle *p1,TParticle *p2,TParticle *p3, Double_t etamax);
   
   void DisableMCQA() { fEnableMCQA = kFALSE; }
   
   AliESDEvent 	*fESD;    	   // ESD object
   AliAODEvent 	*fAOD;    	   // AOD object
   TString 	fAnalysisType;	   // "ESD" or "AOD"	
   
   Long64_t fEvt;			   // event number
   TList* fOutList;             // two outputs, one for mcqa and the other for global 
   
   
   Bool_t fEnableMCQA;	
   
   TH1F* fhNevt;
   Float_t fEtaAbs;   // eta limit considered
   Float_t fEtaAbsMax;   // max eta limit considered for PID and background
   Bool_t fFillBkg;   // fill the histograms concerning background
   Int_t fNcuts; // number of selection cuts for MC candidates
   Double_t* fCuts;    //[fNcuts] cut values
   TString* fCutNames; //[fNcuts] names of the cuts variables
   Double_t fLooserPtTrack;
   Double_t fInvMassCut;
   Double_t *fThreesigmas; //[2] number of sigmas
   
   ClassDef(AliAnalysisTaskCountLcEta, 2);	// adding the class to ROOT
};

#endif
