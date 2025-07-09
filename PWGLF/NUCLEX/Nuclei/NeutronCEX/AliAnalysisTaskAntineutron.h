/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskAntineutron_H
#define AliAnalysisTaskAntineutron_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskAntineutron : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskAntineutron();
                                AliAnalysisTaskAntineutron(const char *name);
        virtual                 ~AliAnalysisTaskAntineutron();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
  
    private:
        
    
    	Bool_t 			fSimulation;    // analysis of MC simulation
        
        AliInputEventHandler*   fEventHandler;  // for ESDs or AODs
        AliESDEvent*            fESD;           // input event
        AliPIDResponse*         fPIDResponse;   // pid response object
        AliESDpid*              fESDpid;        // basic TPC object
        TList*                  fOutputList;    // output list
        AliMCEvent* 		    fMCevent;       // monte carlo event
        UChar_t                 fITSmap;        // ITS cluster map
        const UChar_t kSPDL1 = 0x01; // SPD Layer 1
        const UChar_t kSPDL2 = 0x02; // SPD Layer
        const UChar_t kSDDL1 = 0x04; // SDD Layer 1
        const UChar_t kSDDL2 = 0x08; // SDD Layer 2
        const UChar_t kSSDL1 = 0x10; // SSD Layer 1
        const UChar_t kSSDL2 = 0x20; // SSD Layer 2
    
        // Trees and variables
        TTree*          fTree_Antip = nullptr;
        Int_t           fmom_pdg;
        Double_t        fmom_E;
        Double_t        fmom_process;
        Double_t        fAntipP;
        Double_t        fAntip_Px;
        Double_t        fAntip_Py;
        Double_t        fAntip_Pz;
        Double_t        fAntipPt;
        Double_t        fAntipE;
        Double_t        fAntipY;
        Double_t        fAntipEta;
        Double_t        fAntipTPCsignal;
        Double_t        fAntipITSsignal;
        Double_t        fAntipTPCpoints;
        Double_t        fAntipTOFsignal;
        Int_t           fAntipncTPC;
        Int_t           fAntipncITS;
        Double_t        fAntipDCAxy;
        Double_t        fAntipDCAz;
    
        TTree*          fTree_proton = nullptr;
        Int_t           fpmom_pdg;
        Double_t        fpmom_E;
        Double_t        fpmom_process;
        Double_t        fpP;
        Double_t        fp_Px;
        Double_t        fp_Py;
        Double_t        fp_Pz;
        Double_t        fpPt;
        Double_t        fpE;
        Double_t        fpY;
        Double_t        fpEta;
        Double_t        fpTPCsignal;
        Double_t        fpITSsignal;
        Double_t        fpTPCpoints;
        Double_t        fpTOFsignal;
        Int_t           fpncTPC;
        Int_t           fpncITS;
        Double_t        fpDCAxy;
        Double_t        fpDCAz;
    
        TTree*          fTree_pair = nullptr;
        Double_t        fPairPSumMag;
        Double_t        fPairVtxSep;
        Double_t        fPairMinDCA;
        Double_t        fPairSecRadius;
        Double_t        fPairVtxDistance;
        Double_t        fSecVtxX, fSecVtxY, fSecVtxZ;
        Double_t        fPrimVtxX, fPrimVtxY, fPrimVtxZ;
    
        // Histograms
        TH1F*           fHistPt;
        TH1F*           fHistNEvents;
        // MC Antineutron histograms
        TH1F*		    fHistAntinEta;
        TH1F*           fHistAntinEk;
        TH1F*           fHistAntinEk_ITScuts;
        TH1F*		    fHistAntipEk;
        TH1F*           fHistpEk;
        // MC CEX histograms
        TH2F* 		fHistcexPx;
        TH2F*		fHistcexPy;
        TH2F*		fHistcexPz;
        TH2F*		fHistcexEk;
        TH2F*		fHistcexP;
        TH1F*       fHistcexDP;
        TH2F*       fHistcexVrtx;
        TH1F*       fHistmomcexpdg;
        TH1F*       fHistcexEta_antin;
        TH1F*       fHistcex_ITScuts;
        TH1F*       fHistcex_angle;
        // MC CEX histograms normalized
        TH1F*       fHistncexPx;
        TH1F*       fHistncexPy;
        TH1F*       fHistncexPz;
        TH1F*       fHistncexEk;
        TH1F*       fHistncexP;
        TH1F*       fHistncex_ITScuts;
        
        TH1F*        fHistRecoAngle;
        TH1F*        fHistRecoP;
        TH1F*        fHistRecoradio;
        TH1F*        fHistRecovtx;
        TH1F*        fHistDCApair;
        TH1F*        fHistRecoVtxSep;
        TH2F*        fHistRecosvtxXY;
        TH1F*        fHistRecosvtxZ;
    
        AliAnalysisTaskAntineutron(const AliAnalysisTaskAntineutron&); // not implemented
        AliAnalysisTaskAntineutron& operator=(const AliAnalysisTaskAntineutron&); // not implemented

        ClassDef(AliAnalysisTaskAntineutron, 1);
};

#endif
