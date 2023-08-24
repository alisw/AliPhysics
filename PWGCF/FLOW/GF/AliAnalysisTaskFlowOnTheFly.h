/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskFlowOnTheFly_H
#define AliAnalysisTaskFlowOnTheFly_H

#include "AliAnalysisTaskSE.h"
#include "AliGFWWeights.h"
#include "CorrelationCalculator.h"
#include "AliEventCuts.h"
#include <TComplex.h>
#include "AliMultSelection.h"

#include <TObject.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TComplex.h>
#include <TBits.h>
#include <TRandom3.h>
// AliRoot includes
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include <THnSparse.h>
// Preload to boost program
class TList;
class TF1;
class TH1;
class TH2;
class TH3F;
class TH1F;
class TH2F;
class TH3D;
class TProfile;
class TProfile2D;
class TComplex;
class AliESDEvent;
class AliAODEvent;
class AliVEvent;
class AliVTrack;
class AliVVertex;
class AliESDtrackCuts;
class AliAODITSsaTrackCuts;
class AliInputEventHandler;


class PhysicsProfileFlowOnTheFly : public TObject {
	public:
        PhysicsProfileFlowOnTheFly();
		// Physics profiles
		TProfile*	 fChsc4242;		             	//! SC(4,2)
		TProfile*	 fChsc4242_Gap0;			//! SC(4,2) |#Delta#eta| > 0.0
		TProfile*	 fChsc4242_Gap2;			//! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap4;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap6;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap8;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap10;                       //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap12;                       //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*	 fChsc4242_3sub;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subMMLRA;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subMMLRB;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subLLMRA;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subLLMRB;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subRRMLA;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subRRMLB;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4224_3sub;			//! SC(4,2)_B 3subevent method
		TProfile*	 fChsc4242_3subGap2;			//! SC(4,2)_A 3subevent method |#Delta#eta| > 0.2
		TProfile*	 fChsc4224_3subGap2;			//! SC(4,2)_B 3subevent method |#Delta#eta| > 0.2
		TProfile*	 fChsc3232;				//! SC(3,2)
		TProfile*	 fChsc3232_Gap0;			//! SC(3,2) |#Delta#eta| > 0.0
		TProfile*	 fChsc3232_Gap2;			//! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap4;                        //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap6;                        //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap8;                        //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap10;                       //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap12;                       //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*     	 fChsc3232_3sub;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subMMLRA;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subMMLRB;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subLLMRA;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subLLMRB;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subRRMLA;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subRRMLB;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3223_3sub;			//! SC(3,2)_B 3subevent method
		TProfile*     	 fChsc3232_3subGap2;			//! SC(3,2)_A 3subevent method |#Delta#eta| > 0.2
		TProfile*     	 fChsc3223_3subGap2;			//! SC(3,2)_B 3subevent method |#Delta#eta| > 0.2

		// Standard correlation profiles for different harmonics
		TProfile*	 fChc422;          //!
		TProfile*        fChc532;          //!
		TProfile*        fChc422_Gap0A;    //!
		TProfile*        fChc422_Gap0B;    //!
		TProfile*        fChc532_Gap0A;    //!
		TProfile*        fChc532_Gap0B;    //!
		TProfile*        fChc422_Gap2A;    //!
		TProfile*        fChc422_Gap2B;    //!
		TProfile*        fChc532_Gap2A;    //!
		TProfile*        fChc532_Gap2B;    //!
		TProfile*        fChc422_Gap4A;    //!
		TProfile*        fChc422_Gap4B;    //!
		TProfile*        fChc532_Gap4A;    //!
		TProfile*        fChc532_Gap4B;    //!
		TProfile*        fChc422_Gap6A;    //!
		TProfile*        fChc422_Gap6B;    //!
		TProfile*        fChc532_Gap6A;    //!
		TProfile*        fChc532_Gap6B;    //!
		TProfile*        fChc422_Gap8A;    //!
		TProfile*        fChc422_Gap8B;    //!
		TProfile*        fChc532_Gap8A;    //!
		TProfile*        fChc532_Gap8B;    //!
		TProfile*        fChc422_Gap10A;    //!
		TProfile*        fChc422_Gap10B;    //!
		TProfile*        fChc532_Gap10A;    //!
		TProfile*        fChc532_Gap10B;    //!
		TProfile*        fChc422_Gap12A;    //!
		TProfile*        fChc422_Gap12B;    //!
		TProfile*        fChc532_Gap12A;    //!
		TProfile*        fChc532_Gap12B;    //!

		TProfile*	 fChcn2[5]; 			//! <<2>> in unit bins of Ntrks
		TProfile*    fChcn2_Gap0[5];  		//! <<2>> |#Delta#eta| > 0.0
		TProfile*	 fChcn2_Gap2[5];  		//! <<2>> |#Delta#eta| > 0.2
		TProfile*	 fChcn2_Gap4[5];  		//! <<2>> |#Delta#eta| > 0.4
		TProfile*    fChcn2_Gap6[5];        //! <<2>> |#Delta#eta| > 0.6
		TProfile*	 fChcn2_Gap8[5];  		//! <<2>> |#Delta#eta| > 0.8
		TProfile*	 fChcn2_Gap10[5];  		//! <<2>> |#Delta#eta| > 1.0
		TProfile*	 fChcn2_Gap12[5];  		//! <<2>> |#Delta#eta| > 1.2
		TProfile*	 fChcn2_Gap14[5];  	        //! <<2>> |#Delta#eta| > 1.4
		TProfile*	 fChcn2_Gap16[5];  	        //! <<2>> |#Delta#eta| > 1.6
		TProfile*	 fChcn2_Gap18[5];  	        //! <<2>> |#Delta#eta| > 1.8

		TProfile*	 fChcn2_3subLM[5];  		//! <<2>> left vs. middle subevent
		TProfile*	 fChcn2_3subRM[5];  		//! <<2>> middle vs. right subevent
		TProfile*	 fChcn2_3subLR[5];  		//! <<2>> middle vs. right subevent
		TProfile*	 fChcn2_3subGap2LM[5];          //! <<2>> left vs. middle subevent |#Delta#eta| > 0.2
		TProfile*	 fChcn2_3subGap2RM[5];          //! <<2>> middle vs. right subevent |#Delta#eta| > 0.2

		TProfile*	 fChcn4[5];  			//! <<4>> in unit bins of Ntrks
		TProfile*	 fChcn4_Gap0[5];  		//! <<4>> |#Delta#eta| > 0.0
		TProfile*	 fChcn4_Gap2[5];  		//! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap4[5];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap6[5];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap8[5];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap10[5];               //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap12[5];               //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_3subMMLR[5];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subLLMR[5];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subRRML[5];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subGap2[5];            //! <<4>> 3subevent method |#Delta#eta| > 0.2

	
		//5,6 particle correlation
		//TProfile* fChc5_A42222;  		//! <<5>> 
		//TProfile* fChc5_A52322;  		//! <<5>> 
		//TProfile* fChc6_222222;  		//! <<6>> 
		//TProfile* fChc6_322322;  		//! <<6>> 

		//Addtional 6 particle correlation
		//TProfile* fChsc6222_Gap0;  		//! <<4>> |#Delta#eta| > 1.0
		//TProfile* fChsc6222_Gap10;  		//! <<4>> |#Delta#eta| > 1.0
		//TProfile* fChsc633_Gap0A;  		//! <<3>> |#Delta#eta| > 1.0
		//TProfile* fChsc633_Gap10A;  		//! <<3>> |#Delta#eta| > 1.0
		TProfile*	fChcn6[6];  			//! <<6>> in unit bins of Ntrks
		TProfile*   fChcn6_Gap10[6];               //! <<6>> |#Delta#eta| > 1.0
		TProfile*   fChcn6_Gap0[6];               //! <<6>> |#Delta#eta| > 0.

		// 8 particle correlation
		TProfile*	fChcn8[6];  			//! <<8>> in unit bins of Ntrks
		TProfile*   fChcn8_Gap0[6];               //! <<8>> |#Delta#eta| > 1.0

		private:
		ClassDef(PhysicsProfileFlowOnTheFly, 1);    //Analysis task
};


class AliAnalysisTaskFlowOnTheFly : public AliAnalysisTaskSE  
{
    public:
        enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kCharUnidentified, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter
                                AliAnalysisTaskFlowOnTheFly();
                                AliAnalysisTaskFlowOnTheFly(const char *name);
        virtual                 ~AliAnalysisTaskFlowOnTheFly();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
		virtual void   			NotifyRun();
        virtual void            Terminate(Option_t* option);

        
        virtual void   SetFilterbit(Int_t bit){fFilterbit = bit;}
		virtual void   SetFilterbitDefault(Int_t bit){fFilterbitDefault = bit;}
		virtual void   SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void   SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void   SetVtxCutDefault(Double_t vtxCut){fVtxCutDefault = vtxCut;} // The vtxCut for NtrksCounter
		virtual void   SetMinPt(Double_t minPt){fMinPt = minPt;}
		virtual void   SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
		virtual void	SetEventWeightSetOne(Bool_t useOne){
			fEventWeightSetToOne=useOne;	//Set Event Weight to 1
		}
		virtual void   SetUseCL1Cent(Bool_t CL1cent){fUseCL1Centrality = CL1cent;}
		virtual void   SetTrigger(Int_t trig){fTrigger = trig;}
		virtual void   SetUseWeigthsRunByRun(Bool_t bRunByRun = kTRUE) { fFlowRunByRunWeights = bRunByRun; }
		virtual void   SetUsePeriodWeigths(Bool_t weight = kTRUE) { fFlowPeriodWeights = weight; }
		virtual void   SetUseWeights3D(Bool_t use = kTRUE) { fFlowUse3Dweights = use; }
		virtual void   SetfUseImpactXaxis(Bool_t use = kTRUE) { fUseImpactXaxis = use; }
        virtual void   SetSystFlag(int flag) { fCurrSystFlag = flag; }
        virtual int    GetSystFlag() { return fCurrSystFlag; }
		void SetAMPTCentralityMap(std::vector<double> b, std::vector<double> cent) { for(size_t i(0); i<b.size(); ++i) centralitymap[b[i]]=cent[i]; }
        
        //===================================================================================

    private:

        AliAnalysisTaskFlowOnTheFly(const AliAnalysisTaskFlowOnTheFly&); // not implemented
        AliAnalysisTaskFlowOnTheFly& operator=(const AliAnalysisTaskFlowOnTheFly&); // not implemented
        

		AliEventCuts	fEventCuts;					// Event cuts

		// Cuts and options
		Int_t			fFilterbit;				// flag for filter bit
		Int_t			fFilterbitDefault;			// flag for filter bit (for NtrksCounter)
		Double_t		fEtaCut;				// Eta cut used to select particles
		Double_t		fVtxCut;				// Vtx cut on z position in cm
		Double_t		fVtxCutDefault;				// Vtx cut on z position in cm (for NtrksCounter)
		Double_t		fMinPt;					// Min pt - for histogram limits
		Double_t		fMaxPt;					// Max pt - for histogram limits
		Bool_t			fUseImpactXaxis;		// Use Impact parameter as X axis
		AliMCEvent 		*fMCEvent; 				//! MC event
		Double_t 		fImpactParameterMC;			//! Impact Parameter
		Bool_t			fEventWeightSetToOne;		// Set Event weight to 1
		Bool_t			fAddTPCPileupCuts;			// Pile up Cuts in TPC
		Double_t		fESDvsTPConlyLinearCut;     // ESDvsTPConlyLinearCut : default = 15000
		Bool_t			fUseCL1Centrality;				// switch to choose whether I want to use DCAz cut or not (for systematic studies, otherwise it is in FB selection by default)
		Int_t			fTrigger;				// flag for trigger
		Int_t			fAliTrigger;				// name for trigger
                Int_t                   fCurrSystFlag;                          // Systematics flag


		Double_t			fCurrCentrality;	//! Centrality of Ongoing Event

		// Output objects
		TList*			fListOfObjects;			//! Output list of objects

		// NUE
		TList*			fTrackEfficiency;		//! file with tracking efficiency
		//TH1D*			hTrackEfficiency;		//! histogram with tracking efficiency
		TH1D*			hTrackEfficiencyRun;            //! histogram with tracking efficiency

		// NUA
		bool fFlowRunByRunWeights;                              // flag of whether get the Run by run weight 
		bool fFlowPeriodWeights;                                // flag of whether to use period weight
		bool fFlowUse3Dweights;                                 // flag of whether to use 3d weight

		//
		TList*                  fFlowWeightsList;               //! flowWightsList
		AliGFWWeights*          fWeightsSystematics;            //! Weights for systematics


		TH3F*			hPhiWeight;			//! 3D weight for all periods except LHC15ijl

		// Event histograms
		TH1D*			hEventCount;			//! counting events passing given event cuts
		TH1F*			fCentralityDis;			//! distribution of centrality percentile using V0M estimator


		// // Track histograms
		TH1D*			hPhi;			//! distribution of Phi
		TH1D*			hPhiBefore;			//! distribution of Phi
		TH1D*			fEtaDis;			//! distribution of Phi
		TH1D*			fPtDis;			//! distribution of Phi
		TH1D*			fIP;			//! distribution of Impact parameter

		// Global variables
		double NtrksCounter = 0;               //!
		double NTracksCorrected = 0;			//!
		double NTracksUncorrected = 0;			//!
		int NtrksAfter = 0;                 //!
		int NtrksAfterGap0M = 0;           //!
		int NtrksAfterGap0P = 0;           //!
		int NtrksAfterGap2M = 0;           //!
		int NtrksAfterGap2P = 0;           //!
		int NtrksAfterGap4M = 0;           //!
		int NtrksAfterGap4P = 0;           //!
		int NtrksAfterGap6M = 0;           //!
		int NtrksAfterGap6P = 0;           //!
		int NtrksAfterGap8M = 0;           //!
		int NtrksAfterGap8P = 0;           //!
		int NtrksAfterGap10M = 0;           //!
		int NtrksAfterGap10P = 0;           //!
		int NtrksAfterGap12M = 0;           //!
		int NtrksAfterGap12P = 0;           //!
		int NtrksAfterGap14M = 0;           //!
		int NtrksAfterGap14P = 0;           //!
		int NtrksAfter3subL = 0;            //!
		int NtrksAfter3subM = 0;            //!
		int NtrksAfter3subR = 0;            //!
		
        
		PhysicsProfileFlowOnTheFly multProfile;          //!
		PhysicsProfileFlowOnTheFly multProfile_bin[30];  //!
        
		TRandom3 rand;                      //!
		Int_t bootstrap_value;              //!

		CorrelationCalculator correlator;   //!

		std::map<double,double> centralitymap;
		double xbins[3000+10];              //!
		int nn;                             //!
		//TH1F* MyEventNumber;					//!
        
		void CalculateProfile(PhysicsProfileFlowOnTheFly& profile, double Ntrks);
		void InitProfile(PhysicsProfileFlowOnTheFly& profile, TString);
		void ProcessOnTheFly();
		AliMCEvent *getMCEvent();
		double getAMPTCentrality();
        
        //===================================================================================

        ClassDef(AliAnalysisTaskFlowOnTheFly, 1);
};

#endif
