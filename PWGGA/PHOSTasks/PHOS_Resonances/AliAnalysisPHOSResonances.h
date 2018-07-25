#ifndef AliAnalysisPHOSResonances_cxx
#define AliAnalysisPHOSResonances_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for extracting spectra of resonances in PHOS plus tracks
// Authors: Dmitri Peresunko
// 12-July-2018

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliESDtrackCuts;
class AliTriggerAnalysis;
class AliPIDResponse;
class AliAODTrack ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisPHOSResonances : public AliAnalysisTaskSE {
public:
	AliAnalysisPHOSResonances(const char *name = "AliAnalysisPHOSResonances");
        AliAnalysisPHOSResonances(const AliAnalysisPHOSResonances&); 
	AliAnalysisPHOSResonances& operator=(const AliAnalysisPHOSResonances&); 
	virtual ~AliAnalysisPHOSResonances();

	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
        void SetListOfChannels(Int_t *array, Int_t size=18){fListOfChannels.Set(size,array);}
	
private:

        enum Channels_t{gg,gpp,gpm,gppDCA,gpmDCA,gpi0,gpi0Merg,pi0pi0,pi0pi0M,pi0pp,pi0pm,pi0ppDCA,pi0pmDCA,gLam, pi0Lam,ee,gpippim,pi0pippim} ;
        enum PID_t {All,TOF,Disp,CPV,Both,NoPi0,NoPi0Eta} ;
    
	void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
// 	void FillHistogram(Int_t remi,Channels_t ch,PID_t pid,Double_t x) const ; //Fill 1D histogram witn name key
	
	void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
	void FillHistogram(Int_t remi,Channels_t ch,PID_t pid,Double_t x, Double_t y, Double_t z) const ; //Fill 2D histogram witn name key
	void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
	void SelectGamma() ;
	void SelectElectrons() ;
        void SelectHadrons() ;
	void SelectLambda() ;
        Bool_t AcceptTrack(const AliAODTrack *t);
        Double_t Pi0Mass(Double_t pt) ;
        Double_t Pi0Width(Double_t pt) ;
        Double_t EtaMass(Double_t pt) ;
        Double_t EtaWidth(Double_t pt) ;
        Double_t PionDispCut(Double_t m02, Double_t m20, Double_t E) ;
       
 
private:
	THashList * fOutputContainer;  //! final histogram container
	TArrayI fListOfChannels;     //  list of active channels

	Int_t fnPID;                 //  NUmber of PID cuts
// 	AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation
	AliPIDResponse     *fPIDResponse;     //! PID response 
	
	AliAODEvent  * fEvent ;      //! Current event
	TClonesArray * fGamma;       //! List of selected photons
	TClonesArray * fPi0 ;        //! List of selected pi0s
	TClonesArray * fPi0Merged ;  //! List of selected merged pi0s
	TClonesArray * fTracksElm ;  //! list of selected electrons
	TClonesArray * fTracksElp ;  //! list of selected positrons
	TClonesArray * fTracksPim ;  //! list of selected pi-
	TClonesArray * fTracksPip ;  //! list of selected pi+
	TClonesArray * fTracksKm ;   //! list of selected K+
	TClonesArray * fTracksKp ;   //! list of selected K- 
	TClonesArray * fTracksPm ;   //! list of selected p+ 
	TClonesArray * fTracksPp ;   //! list of selected p-  
	TClonesArray * fLambda ;     //! list of selected Lambdas
	TList * fMixGamma ;          //!
	TList * fMixPi0 ;            //!
	TList * fMixPi0Merged ;      //!
	TList * fMixElm ;            //!
	TList * fMixElp ;            //!
	TList * fMixTracksPim ;      //!
	TList * fMixTracksPip ;      //!
	TList * fMixTracksKm ;       //!
	TList * fMixTracksKp ;       //!
	TList * fMixTracksPm ;       //!
	TList * fMixTracksPp ;       //! 
	TList * fMixLambda ;         //!
        
        TH3F * fhHistos[300];      //!
        TH2F * fhPr;               //!

	ClassDef(AliAnalysisPHOSResonances, 1); 
};
#endif
