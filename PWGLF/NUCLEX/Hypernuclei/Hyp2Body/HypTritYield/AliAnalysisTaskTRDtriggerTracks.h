/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskTRDtriggerTracks_H
#define AliAnalysisTaskTRDtriggerTracks_H
#include "AliTRDonlineTrackMatching.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"

class AliPIDResponse;

class AliAnalysisTaskTRDtriggerTracks : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskTRDtriggerTracks();
                                AliAnalysisTaskTRDtriggerTracks(const char *name);
        virtual                 ~AliAnalysisTaskTRDtriggerTracks();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
			AliESDEvent*			fESD;           //! input event
			AliPIDResponse*		fPIDResponse; 		//! pid response object
			TList* 						fOutputList;    //! output list
			AliEventCuts      fEventCuts;         //< event cuts as advised by PDG (AliEventCuts) 		

			AliESDInputHandler     *fInputHandler;        //!<! Input handler 					
			TTree                  *fTree;              //< tree containing TRD information
			Int_t						tTrigMB;			// trigger info
			Int_t						tTrigHMV0;
			Int_t						tTrigHMSPD;
			Int_t						tTrigHNU;
			Int_t						tTrigHQU;
			
			Float_t 				tCharge;
			Float_t					tPt;				// 
			Float_t					tPx;				// 
			Float_t					tPy;				// 
			Float_t					tPz;				// 
			Float_t					tY;
			Float_t					tP;					//
			Float_t					tTPCDEdx;
			Float_t					tTOFSignal;
			Float_t					tDcaXY;				// impact parameters
			Float_t					tDcaZ;			

			Int_t										tTRDtrigHNU;
			Int_t										tTRDtrigHQU;			

			Int_t										tTRDPid;
			Int_t										tTRDnTracklets;
			Int_t										tTRDPt;
			Int_t										tTRDLayerMask;
			Float_t									tTRDSagitta;
			
			Float_t 								tMultV0M;
			Float_t 								tMultSPDCluster;
			Float_t		 							tSPDTracklets;

			TH1F * histMultV0HNU; 
			TH1F * histMultSPDHNU; 
			TH1F * histMultRefHNU;
			TH1F * histMultV0HQU;
			TH1F * histMultSPDHQU;
			TH1F * histMultRefHQU;

        AliAnalysisTaskTRDtriggerTracks(const AliAnalysisTaskTRDtriggerTracks&); // not implemented
        AliAnalysisTaskTRDtriggerTracks& operator=(const AliAnalysisTaskTRDtriggerTracks&); // not implemented

        ClassDef(AliAnalysisTaskTRDtriggerTracks, 1);
};

#endif
