#ifndef ALIANALYSISTASKUESPECTRADPHI_H
#define ALIANALYSISTASKUESPECTRADPHI_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH3D.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliESDVertex.h>
#include <AliAnalysisFilter.h>
#include <AliVHeader.h>
#include <AliESDtrackCuts.h>
#include <AliEventCuts.h>
#include <AliMultSelection.h>
#include <AliCentrality.h>
#include <AliMultSelectionTask.h>

class AliAnalysisTaskUeSpectraDphi : public AliAnalysisTaskSE {

 public:
		AliAnalysisTaskUeSpectraDphi();
		AliAnalysisTaskUeSpectraDphi(const char *name);
		virtual ~AliAnalysisTaskUeSpectraDphi();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);

		virtual void     AnalyzeESD(AliESDEvent* esd);
		virtual void     AnalyzeMC(AliMCEvent* fMCEvent);
		virtual void     AnalyzeESDforDCA(AliESDEvent* esd);
		virtual void     SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		virtual void     SetTrackCutsDCA(AliAnalysisFilter* fTrackFilterDCA);
		virtual Double_t DeltaPhi(Double_t phia, Double_t phib,
				 Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
		virtual Double_t DeltaPhiOA(Double_t phia, Double_t phib, Double_t range = TMath::Pi());

		Double_t GetVtxCut() { return fVtxCut; }
		Double_t GetEtaCut() { return fEtaCut; }
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		
		virtual void  SetTrigger(UInt_t ktriggerInt = AliVEvent::kINT7) {ftrigBit = ktriggerInt;}
		virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}

 protected:
		void		 Exit(const char *msg);

 private:
		virtual Bool_t selectVertex2015pp(AliESDEvent *esd,Bool_t checkSPDres, Bool_t requireSPDandTrk,Bool_t checkProximity);
		virtual Bool_t IsGoodSPDvertexRes(const AliESDVertex * spdVertex);
		virtual Bool_t isMCEventTrueINEL0(AliMCEvent* fMCEvent);
		virtual ULong64_t GetEventIdAsLong(AliVHeader* header);

		//
		// Cuts and options
		//
		Double_t     fVtxCut;                  // Vtx cut on z position in cm
		Double_t     fEtaCut;                  // Eta cut used to select particles

		AliMCEvent *        fMCEvent;           //! MC event
		AliStack   *        fMCStack;           //!
		AliESDEvent*        fESD;               //!  ESD object
		TRandom*            fRandom;            //! random number generator
		AliEventCuts        fEventCuts;         //!
		AliAnalysisFilter*  fTrackFilter;       //! 
		AliAnalysisFilter*  fTrackFilterDCA;    //!
		AliMultSelection *  fMultSelection;     //!
	

		UInt_t        ftrigBit;		         //
		Bool_t        fPileUpRej;                // kTRUE is pile-up is rejected
		Bool_t        fAnalysisMC;               //

		//
		// Help variables
		//
		Short_t      fTriggeredEventMB;            // 1 = triggered, 0 = not trigged (MC only) */
		Bool_t       fisPS; 		  	   //
		Bool_t       fisTracklet;		   //
		Bool_t       isINEL0Rec;                   //
		Bool_t       isINEL0True;                  //  
		Bool_t       fisMCvtxInZcut;               //
		Int_t        fRun;                         // run no
		ULong64_t    fEventId;                     // unique event id
		Float_t      fdcaxy;                       // 
		Float_t      fdcaz;                        //

		//
		// Output objects
		//
		TList* fListOfObjects;                      //! Output list of objects
		TH1D * fVtxBeforeCuts;                      //! Vertex z dist before cuts
		TH1D * fVtxAfterCuts;                       //! Vertex z dist after cuts
		TH1D * hSelEv;                              //! No of accepted events
		TH1D * hINEL0;                              //!
		TH1D * hEvMultSel;                          //!
		TH1D * hEvDphiSel;                          //!

		TH1D * hPS;	                	    //!		
		TH1D * hVtxPS;	                	    //!				

		TH3D *ptvstrackletsvsdcaData;	            //!	
		TH3D *ptvstrackletsvsdcaPrim;	            //!
		TH3D *ptvstrackletsvsdcaDecs;	            //!
		TH3D *ptvstrackletsvsdcaMatl;	            //!
		TH3D *ptvstrackletsvsdcacentralData;	    //!
		TH3D *ptvstrackletsvsdcacentralPrim;	    //!
		TH3D *ptvstrackletsvsdcacentralDecs;	    //!
		TH3D *ptvstrackletsvsdcacentralMatl;	    //!

		TH1D * hpT;	                	    //!	
		TH1D * hEta;	                	    //!	
		TH1D * hPhi;	                	    //!	
		TH1D * hPtL;	                	    //!	
		TH1D * hEtaL;	                	    //!	
		TH1D * hPhiL;	                	    //!	
		TH1D * hDphi;	                	    //!	
		TH1D * hRefMult08;	                    //!	
		TH1D * hV0Mmult;	                    //!	
		TH2D * hpTvsDphiOA;	                    //!	
		TH2D * hpTvsDphiSA;	                    //!
		TH2D * hMultvsDphiOA;	                    //!
		TH2D * hMultvsDphiSA;	                    //!
		TH3D * hMultvspTvsDphi;	                    //!
		TH3D * hMultvspTvsDphiWLP;	            //! With leading particle


		TH1D  *hDPhiNchTSGT12;                       //!  NchTS > 12 in both Left and Right
		TH1D  *hDPhiNchTSLT8;                        //!  NchTS <  8 in both Left and Right
		TH1D  *hDPhiNchTSGT12LT12;                   //!  NchTS > 12 in Left or Right  and < 12 in left or Right
		TH1D  *hDPhiNchTSGT12LT8;                    //!  NchTS > 12 in Left or Right  and < 8 in left or Right
		TH1D  *hDPhiNchTSGT12GT0;                    //!  NchTS > 12 in Left or Right  and > 0 in left or Right
		
		TH2D  *hMultTSNchTSGT12;                     //!  NchTS > 12 in both Left and Right
		TH2D  *hMultTSNchTSLT8;                      //!  NchTS <  8 in both Left and Right
		TH2D  *hMultTSNchTSGT12GT0;                  //!  NchTS > 12 in Left or Right  and > 0 in left or Right
		TH2D  *hMultTSNchTSGT12LT12;                 //!  NchTS > 12 in Left or Right  and < 12 in left or Right
		TH2D  *hMultTSDNchTSGT12LT8;                 //!  NchTS > 12 in Left or Right  and < 8 in left or Right

		TH1D * hNchTSLeft[10];                       //!
		TH1D * hNchTSRight[10];                      //!
		TH1D * hSumptTSLeft[10];                     //!
		TH1D * hSumptTSRight[10];                    //!
		TH2D * hNchTSLeftvsNchTSRight[10];           //!
		TH2D * hSumptTSLeftvsSumptTSRight[10];       //!
		TH3D * hNchTSLeftvsNchTSRightvsDphi;        //!
	
		//MC....
		TH1D * hINEL0MCTrig;	                    //!
		TH1D * hINEL0MCTrue;	                    //!
		TH1D * hPS_MC;		                    //!		
		TH1D * hVtxPS_MC;	                    //!		
		TH1D * hpTMCTrue;	                    //!
		TH1D * hEtaMCTrue;	                    //!
		TH1D * hPhiMCTrue;	                    //!
		TH1D * hPtLMCTrue;	                    //!
		TH1D * hEtaLMCTrue;	                    //!
		TH1D * hPhiLMCTrue;	                    //!
		TH1D * hDphiMCTrue;	                    //!
		TH2D * hpTvsDphiOAMCTrue;	            //!
		TH2D * hpTvsDphiSAMCTrue;	            //!
		TH2D * hMultvsDphiOAMCTrue;	            //!
		TH2D * hMultvsDphiSAMCTrue;	            //!
		TH3D * hMultvspTvsDphiMCTrue;	            //!
		TH3D * hMultvspTvsDphiWLPMCTrue;	    //! With leading particle

		TH3D *sigLossTrueINEL0;
		TH3D *sigLossTrigINEL0;
		TH3D *primaries;
		TH3D *secondaries;
		
		Double_t ftrackmult08;                      //
		Double_t fv0mpercentile;                    //

		ClassDef(AliAnalysisTaskUeSpectraDphi, 1);    //Analysis task for high pt analysis

};

#endif
