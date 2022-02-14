/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskPolarizationTestJP_H
#define AliAnalysisTaskPolarizationTestJP_H

#include <AliAnalysisTaskSE.h>

class AliMuonTrackCuts; 	// Include class for standard muon tack cuts


class AliAnalysisTaskPolarizationTestJP : public AliAnalysisTaskSE  
{
public:
                            AliAnalysisTaskPolarizationTestJP();
                            AliAnalysisTaskPolarizationTestJP(const char *name);
    virtual                 ~AliAnalysisTaskPolarizationTestJP();

    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    virtual void   			    NotifyRun();								  // Implement the Notify run to search for the new parameters at each new runs
   // virtual void            AODAnalysis(AliVEvent *fAOD,Int_t iSelectionCounter);
    
    
    
	  void 					TwoMuonAna(Int_t *pos, Int_t *neg);			    // Analyses two muons and extracs dimuon information
	  void 					TwoMCMuonAna(Int_t *MCpos, Int_t *MCneg);	  // Analyses two MC muons and extracs MC dimuon information
	  void 					SetMC(Bool_t flag){fIsMC = flag;}	          //used to analyse the MC data
    void          SimAnalysis(Int_t iSelectionCounter);                              // simulated analysis  
    void          AODAnalysis(AliVEvent *fAOD,Int_t iSelectionCounter);
  //  void          AODAnalysis(Int_t iSelectionCounter);
    Bool_t       IsTriggered(AliVEvent *fAOD);
    
   Double_t            CosThetaHelicityFrame( TLorentzVector muonPositive,
                                                       TLorentzVector muonNegative,
                                                       TLorentzVector possibleJPsi);
    Double_t            CosThetaCollinsSoper( TLorentzVector muonPositive,
                                                       TLorentzVector muonNegative,
                                                       TLorentzVector possibleJPsi);
    Double_t            CosPhiHelicityFrame( TLorentzVector muonPositive, 
                                                       TLorentzVector muonNegative,
                                                       TLorentzVector possibleJPsi );
    Double_t            CosPhiCollinsSoper( TLorentzVector muonPositive,
                                                       TLorentzVector muonNegative,
                                                       TLorentzVector possibleJPsi );
                                                       
                                                       
   Double_t             TildePhiCalulator(Double_t phi, Double_t costheta);
 
 
   




	  void 					PostAllData();	

  AliMuonTrackCuts* 		fMuonTrackCuts; 					// Use the class as a data member

private:
    
	  Bool_t 					  fIsMC;
	

    AliAODEvent*      fAOD;       		//! input event
    AliMCEvent*				fMC;				//! input MC event

    TList*            fOutputList; 		//! output list
   

   TTree *fRecTree; 			//! analysis tree
   TTree *fGenTree; 			    //! MC tree
	
  
 /*trigger class infoas integers
 
 
 */ 
  
  TLorentzVector fRecPosDaughter; 
  TLorentzVector fRecNegDaughter; 
  TLorentzVector fRecPair_Parent;
  Float_t     fRecPair_ParentMass;
  
  
  
  TLorentzVector fMCPosDaughter; 
  TLorentzVector fMCNegDaughter; 
  TLorentzVector fMCPair_Parent;
  Float_t     fMCPair_ParenttMass;
  
  
  TLorentzVector fRec_ConnectedMCPosDaughter; 
  TLorentzVector fRec_ConnectedMCNegDaughter; 
  TLorentzVector fRec_ConnectedMCPair_Parent;
  Float_t        fRec_ConnectedMCPair_ParenttMass;
  
  
  //trigger decisions
  
  
  Int_t fCMUPDecision;
  Int_t fCMUP6Decision;
	Int_t fCMUP10Decision;
	Int_t fCMUP11Decision;
  Int_t fCMUP13Decision;
	Int_t fCMUP26Decision;
 //znenergy
 	Float_t fZNCEnergy; 
	Float_t fZNAEnergy;
	// Double_t fZPCEnergy; 
	// Double_t fZPAEnergy;
	Float_t fZNATDC[4];
	Float_t fZNCTDC[4];
	// Double_t fZPATDC[4];
	// Double_t fZPCTDC[4];
  
  
  
  //histograms
  
  
  
  
  TH1F*                     fCounterH; 			//! counter for events passing each cut
  TH1D*                     fHistRunCounter;
  //TH1D*                     fHistTriggers;
  TH1D*                     fHistCMUPTriggers;
  TH1D*                     fHistCMUP6Triggers;
  TH1D*                     fHistCMUP10Triggers;
  TH1D*                     fHistCMUP11Triggers;
  TH1D*                     fHistCMUP13Triggers;
  TH1D*                     fHistCMUP26Triggers;
  //TH1I*                     fHistCounter; 




	TClonesArray *fGenPart; 	//! MC particle object
   
  

	//Float_t fMCMuMuM;
 Int_t fRunNumber;
 Int_t fTrgRunNum;

  
  
  


	TTree *fTrgTree; 			    //! trigger info tree
	Int_t fCMUP;
  Int_t fCMUP6;
	Int_t fCMUP10;
	Int_t fCMUP11;
  Int_t fCMUP13;
  Int_t fCMUP26;

  //ANGULAR DISTRIBUTIONS

  Float_t                 fRecHelicityTheta;
  Float_t                 fRecCollinTheta;
  Float_t                 fRecCollinTildePhi;
  Float_t                 fRecHelicityTildePhi;
  Float_t                 fRecHelicityPhi;
  Float_t                 fRecCollinPhi;
  
  
  
  
  Float_t fSimulated_Reconstructed_HelicityTheta;
  Float_t fSimulated_Reconstructed_CollinTheta;
  Float_t fSimulated_Reconstructed_HelicityPhi;
  Float_t fSimulated_Reconstructed_CollinPhi;
  Float_t fSimulated_Reconstructed_CollinTildePhi;
  Float_t fSimulated_Reconstructed_HelicityTildePhi;
    
  
  Float_t                 fMCHelicityTheta;
  Float_t                 fMCCollinTheta;
  Float_t                 fMCCollinTildePhi;
  Float_t                 fMCHelicityTildePhi;
  Float_t                 fMCHelicityPhi;
  Float_t                 fMCCollinPhi;

 






    AliAnalysisTaskPolarizationTestJP(const AliAnalysisTaskPolarizationTestJP&); // not implemented
    AliAnalysisTaskPolarizationTestJP& operator=(const AliAnalysisTaskPolarizationTestJP&); // not implemented

    ClassDef(AliAnalysisTaskPolarizationTestJP, 1);
};

#endif
