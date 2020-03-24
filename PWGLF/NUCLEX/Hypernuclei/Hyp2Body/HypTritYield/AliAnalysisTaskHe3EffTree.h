#ifndef AliAnalysisTaskHe3EffTree_H
#define AliAnalysisTaskHe3EffTree_H

#include "AliAnalysisTaskSE.h"

class AliPIDResponse;
class AliAnalysisTaskHe3EffTree : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskHe3EffTree();
                                AliAnalysisTaskHe3EffTree(const char *name);
        virtual                 ~AliAnalysisTaskHe3EffTree();
        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
  		void SetParamsHe(Double_t params[6]) { for(Int_t i=0; i < 6; i++) fBetheParamsHe[i] = params[i];};

    private:
		AliESDInputHandler     *fInputHandler;		//!<! Input handler
		AliMCEvent* 			mcEvent;			//! MC event
		AliPIDResponse* 		fPIDResponse; 		//! pid response object
  		AliStack               *fStack;             //!<! MC stack
  		Bool_t                  fMCtrue;        	//< flag for MC events
        AliESDEvent*            fESDevent;           	//! input event
		AliEventCuts            fEventCuts;         //< event cuts as advised by PDG (AliEventCuts)
        TList*                  fOutputList;    	//! output list
		TH1F*					fHistEvents;		//! number of events and trigger info
		TH2F*					fHistdEdx;			//! TPC dEdx histogram
		Double_t                fBetheParamsHe[6];	//< Bethe Aleph He3 Parameter + TPC sigma
		Double_t                fBetheParamsT[6];	
  		TTree                  *fTree;              //< tree containing He3 information
		Int_t					tRunNumber;
		Int_t					tTrigMB;			// trigger info
		Int_t					tTrigHMV0;
		Int_t					tTrigHMSPD;
		Int_t					tTrigHNU;
		Int_t					tTrigHQU;

		Float_t					tSPDFiredChips0;	// multiplicity triggers
		Float_t					tSPDFiredChips1;
		Float_t 				tSPDTracklets;
		Float_t 				tSPDCluster;
		Float_t					tV0Multiplicity;
		Float_t					tNTracks;

		Float_t 				tMultV0M;			// multiplicity estimators
		Float_t 				tMultOfV0M;			
		Float_t 				tMultSPDTracklet;	
		Float_t 				tMultSPDCluster;	
		Float_t 				tMultRef05;			
		Float_t 				tMultRef08;			 

		Float_t 				tCharge;

		Float_t					tPt;				// He3 track parameter
		Float_t					tY;
		Float_t					tEta;
		Float_t					tPhi;
		Float_t					tPx;				
		Float_t					tPy;
		Float_t					tPz;
		Float_t					tE;
		

		Float_t					tP;					// PID parameter 
		Float_t					tHeDEdx;
		Float_t					tHeSigma;
		Float_t					tTOFSignalHe;

		Float_t					tDcaXY;				// impact parameters
		Float_t					tDcaZ;
		Float_t					tSigmaYX;
		Float_t					tSigmaXYZ;
		Float_t					tSigmaZ;

		Int_t					tMCtrue;			// MC info
		Int_t					tPrimary;
		Int_t					tWeak;
		Int_t					tMaterial;
		Int_t					tHypertriton;
  		
		TTree                  *fTreeGen;           //< tree containing generated He3 information
		Float_t 				tGenCharge;
		Float_t					tGenPt;				// He3 track parameter
		Float_t					tGenY;
		Int_t					tGenPrimary;			
		Int_t					tGenHypertriton;

		void 					ProcessMCParticles();
		Float_t 				GetInvPtDevFromBC(Int_t b, Int_t c);
		Double_t 				Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
		Float_t 				GetTOFSignalHe3(AliESDtrack& trackHe, Float_t tP);
		void SetBetheBlochParams(Int_t runNumber);
        AliAnalysisTaskHe3EffTree(const AliAnalysisTaskHe3EffTree&); // not implemented
        AliAnalysisTaskHe3EffTree& operator=(const AliAnalysisTaskHe3EffTree&); // not implemented

        ClassDef(AliAnalysisTaskHe3EffTree, 1);
};

#endif
