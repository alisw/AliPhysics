#ifndef ALIANALYSISTASKDIMUONCFCONTAINERBUILDER_H
#define ALIANALYSISTASKDIMUONCFCONTAINERBUILDER_H

/* $Id$ */ 

#include "AliAnalysisTaskSE.h"
#include "TString.h"

//	Analysis task for the building of a dimuon CF container
//	Also some single-muon variables are stored
//	L. Bianchi - Universita' & INFN Torino

class AliCFContainer;
class AliCFManager;

class AliAnalysisTaskDimuonCFContainerBuilder : public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskDimuonCFContainerBuilder();
  AliAnalysisTaskDimuonCFContainerBuilder(const Char_t* name, Bool_t readaod, Bool_t readMC, Bool_t isaccept, Double_t beamEn);
  AliAnalysisTaskDimuonCFContainerBuilder& operator= (const AliAnalysisTaskDimuonCFContainerBuilder& c);
  AliAnalysisTaskDimuonCFContainerBuilder(const AliAnalysisTaskDimuonCFContainerBuilder& c);
  virtual ~AliAnalysisTaskDimuonCFContainerBuilder();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  void     UserCreateOutputObjects();
  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* const io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager() const {return fCFManager;}           // get corr manager
  void           SetQAList(TList* const list) {fQAHistList = list;}

  // Setters and Getters
  Bool_t IsReadAODData   ()   			const {return fReadAODData;}
  void   SetReadAODData    	(Bool_t flag=kTRUE)	  {fReadAODData=flag;}
  void	 SetReadMCinfo     	(Bool_t flag=kTRUE)	  {fReadMCInfo=flag;}
  void   SetIsAccProd	   	(Bool_t flag=kTRUE)	  {fIsAccProduction=flag;}
  void	 SetBeamEnergy     	(Double_t en)	          {fBeamEnergy=en;}
  void   SetChi2Limits     	(Double_t chi2track[])    {fChi2Track[0]=chi2track[0];fChi2Track[1]=chi2track[1];}
  void   SetChi2MatchLimits	(Double_t chi2match[])    {fChi2MatchTrig[0]=chi2match[0];fChi2MatchTrig[1]=chi2match[1];}
  void   SetPtSingMuLimits	(Double_t PtSingle[])     {fPtSingMuCut[0]=PtSingle[0];fPtSingMuCut[1]=PtSingle[1];}
  void   SetThetaSingMuLimits	(Double_t ThetaSingle[])  {fThetaSingMuCut[0]=ThetaSingle[0];fThetaSingMuCut[1]=ThetaSingle[1];}
  void   SetZprimVertLimits	(Double_t Zprimvtx[])     {fzPrimVertexSPD[0]=Zprimvtx[0];fzPrimVertexSPD[1]=Zprimvtx[1];}
  void	 SetCutonZvtxSPD	(Bool_t   cut=kFALSE)	  {fCutOnzVtxSPD=cut;}
  void   SetNContributorsLimits	(Int_t NContr[])       {fNContributors[0]=NContr[0];fNContributors[1]=NContr[1];}
  void	 SetCutonNContributors	(Bool_t   cut=kFALSE)	  {fCutOnNContributors=cut;}
  void	 SetDistinguishTrigClass(Bool_t   dist=kFALSE)	  {fDistinguishTrigClass=dist;}
  void   SetTrigClassMuonName	(TString name = "CMU")	  {fTrigClassMuon=name;}
  void   SetTrigClassInteracName(TString name = "CINT")	  {fTrigClassInteraction=name;}
  void   SetTrigClassMuonSideName(TString name[4])	  {for(Int_t i=0;i<4;i++) fTrigClassMuonSide[i]=name[i];}
  void   SetTrigClassInteracSideName(TString name[4])	  {for(Int_t i=0;i<4;i++) fTrigClassInteractionSide[i]=name[i];}
 
 protected:
  
  Bool_t          	fReadAODData     	;    // flag for AOD/ESD input files
  Bool_t		fReadMCInfo		;    // flag for reading MC info (ESD->Kinematics, AOD->MCbranch)
  Bool_t		fIsAccProduction	;    // flag to activate in case of acceptance MC production (in this case fReadMCInfo==kTRUE)
  AliCFManager   	*fCFManager      	;    // pointer to the CF manager
  TList          	*fQAHistList     	;    // list of QA histograms
  Double_t           	fNevt            	;    // event counter
  Double_t		fBeamEnergy      	;    // Energy of the beam (required for the CS angle)
  TList 		*fOutput         	;    // list of TH in output
  
						     // CUTS ON TRACKS
  Double_t		fChi2Track[2]	 	;    // Cut on chi2 of the tracks ([0]==chi2min, [1]==chi2max)
  Double_t		fChi2MatchTrig[2]	;    // Cut on chi2matchtrigger of the tracks ([0]==chi2Matchmin, [1]==chi2Matchmax)
  Double_t		fPtSingMuCut[2]	 	;    // Cut on pt of single-mu tracks ([0]==ptmin, [1]==ptmax)
  Double_t		fThetaSingMuCut[2]	;    // Cut on polar angle (wrt beam axis) of single-mu tracks ([0]==thetamin, [1]==thetamax)

						     // CUTS ON EVENT
  Double_t		fzPrimVertexSPD[2]	;    // Cut on the z coordinate of the primary vertex in SPD (full ITS for AODs)
  Bool_t		fCutOnzVtxSPD		;    // flag to activate the cut on the z of the primary vertex
  Int_t			fNContributors[2]	;    // Cut on NContributors in SPD
  Bool_t		fCutOnNContributors	;    // flag to activate the cut on NContributors in SPD

						     // CUTS ON THE FIRED TRIGGER CLASS
  TString		fTrigClassMuon		;    // name of the muon trigger class (CMU by default)
  TString		fTrigClassInteraction	;    // name of the interaction trigger class (CINT by default)
  TString		fTrigClassMuonSide[4]	;    // name of the muon trigger classes containing the side
  TString		fTrigClassInteractionSide[4];// name of the interaction trigger classes containing the side
  Bool_t		fDistinguishTrigClass	;    // flag to activate the cut on the fired trigger class
  
  
  
  Double_t Imass  (Double_t e1, Double_t px1, Double_t py1, Double_t pz1, Double_t e2, Double_t px2, Double_t py2, Double_t pz2) const;
  Double_t Rap	  (Double_t e, Double_t pz) const;
  
  Double_t CostCS (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t CostHE (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t PhiCS  (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t PhiHE  (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  
  ClassDef(AliAnalysisTaskDimuonCFContainerBuilder,1);
};

#endif
