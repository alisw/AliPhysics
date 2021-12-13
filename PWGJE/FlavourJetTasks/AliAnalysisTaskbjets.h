/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskbjets_H
#define AliAnalysisTaskbjets_H

#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include  "THistManager.h"
//#include "AliAnalysisTaskParticleInJet.h"

class AliEmcalJet;
class AliAODVertex;
class AliAODTrack;
class AliAODEvent;
class TList;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliEmcalList;
class AliEmcalJetFinder;
class AliAODMCParticle;
class AliMCEvent;
class AliOADBContainer;
class AliVParticle;
#include "AliFJWrapper.h"


class AliAnalysisTaskbjets : public AliAnalysisTaskEmcalJet
{
    public:
                                AliAnalysisTaskbjets();
                                AliAnalysisTaskbjets(const char *name);
        virtual                 ~AliAnalysisTaskbjets();

        virtual void            UserCreateOutputObjects();
        //virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
	Bool_t IsParticleInCone(const AliVParticle* part, const AliEmcalJet* jet, Double_t dRMax);
	virtual Bool_t			 Run();
	
	
	//TObjArray                   fJetCollArray;
    protected:
    
    	void			 ExecOnce();
   //	Bool_t			 FillHistograms();
	
    //	void			 AllocateJetHistograms();
    //	void			 DoJetLoop();
    	
    //	THistManager		 fHistManager;
    
    private:
    
    	
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TH1F*                   fHistPt;        //! dummy histogram
        TH1F*			 fHistNEvents;   //! histogram for total number of events
        TH1F*			 fHistVtxPos;	  //! histogram for primary vertex positions
        TH2F*			 fHistpteta;	  //! histogram for eta distribution
        TH2F*			 fHistptphi;	  //! histogram for phi distribution
        AliPIDResponse*	 fPIDResponse;	  //! pid response object 
        TH2F*			 fHisttpc;	  //! histogram for tpc
        TH1F*			 fHistjet;	  //! histogram for jet-number
        TH1F*			 fHistjetPt;	  //! histogram for jet-Pt
        AliEmcalJet*		 GetPerpendicularPseudoJet (AliEmcalJet*jet_in , bool rev); //
        AliEmcalJet*		 jetrec;	  //
        AliEmcalJet*		 jetmatched;	  //
        AliVParticle*		 vp;	  	  //
        
        Int_t IsMCJetPartonFast(const AliEmcalJet *jet,  Double_t radius,Bool_t &is_udg); //
        
        
	AliJetContainer*  jetconrec; //
        AliJetContainer*  jetcongen; //
       // AliJetContainer*  GetNJets; 
        
        Float_t fJetRecPt; 	//
        Float_t fJetArea;	//
        Float_t fJetMass;	//
        
        //Booleans for settings
        Bool_t   fDoFlavourMatching; //
        
        
        //Cuts
        Int_t fNoJetConstituents;	//
        
        
        AliAnalysisTaskbjets(const AliAnalysisTaskbjets&); // not implemented
        AliAnalysisTaskbjets& operator=(const AliAnalysisTaskbjets&); // not implemented

        ClassDef(AliAnalysisTaskbjets, 1);
};

#endif
