/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCentralJpsi_DG_H
#define AliAnalysisTaskCentralJpsi_DG_H

class TH1;
class TTree;
class TList;
class TFile;
class TBits;
class AliESDtrackCuts;
class AliPIDResponse;
class AliTOFTriggerMask;

#include "AliTimeRangeCut.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCentralJpsi_DG : public AliAnalysisTaskSE
{
    public:
                AliAnalysisTaskCentralJpsi_DG(); // a constructor
                AliAnalysisTaskCentralJpsi_DG(const char *name);
        virtual ~AliAnalysisTaskCentralJpsi_DG(); // a destructor

        virtual void    UserCreateOutputObjects(); // here one can define output objects
        virtual void    UserExec(Option_t* option);	// called for each single event
        virtual void    Terminate(Option_t* option); // usually empty, called at the end

        void    SetIsMC(Bool_t MC){isMC = MC;}
        void    ReplayTriggersMC(AliVEvent *fEvent);
        void    RunMCGenLevel();  
        void    FillMCGenTree(TLorentzVector v);
        void    SetCrossed(Int_t spd[4], TBits &crossed);
        Int_t   GetChipId(Int_t index, Int_t &chipId2, Bool_t debug = 0);
        Bool_t  IsSTGFired(TBits bits, Int_t dphiMin = 4, Int_t dphiMax = 10, Bool_t tolerance = 1);

    private:
        AliPIDResponse  *fPIDResponse;
        AliTimeRangeCut fTimeRangeCut;
        AliESDtrackCuts *fTrackCutsBit4;
        Bool_t          isMC;

        AliVEvent   *fEvent;
        TList       *fOutputList;   //! output list
        TTree       *fTreeJpsi;     //! analysis tree (if MC, then this tree corresponds to MC rec level)
        TTree       *fTreeJpsiMCGen;//! analysis tree on MC gen level
        Int_t       fRunNumber;
        TString     fTriggerName;
        // Histograms:
        TH1D    *hCounterCuts;      //! to count the number of events passing each of the cuts
        TH1D    *hCounterTrigger;   //! to count the number of events per run passing trigger conditions
        TH1D    *hVertexContrib;
        TH1D    *hVertexZ;
        TH2I    *hADdecision;
        TH2I    *hV0decision;
        TH2D    *hTPCdEdx;
        TH2D    *hTPCdEdxMuon;
        TH2D    *hTPCdEdxElectron;
        TH2D    *hPtRecGen;
        // PID, sigmas:
        Double_t fTrk1dEdx;
        Double_t fTrk2dEdx;
        Double_t fTrk1SigIfMu;
        Double_t fTrk1SigIfEl;
        Double_t fTrk2SigIfMu;
        Double_t fTrk2SigIfEl;
        // Kinematics:
        Double_t fPt;   //! transverse momentum
        Double_t fPhi;  //! azimuthal angle
        Double_t fY;    //! rapidity
        Double_t fM;    //! invariant mass
        // Two tracks:
        Double_t fPt1;  //! transverse momenta
        Double_t fPt2;
        Double_t fEta1; //! pseudorapidities
        Double_t fEta2;
        Double_t fPhi1; //! azimuthal angles
        Double_t fPhi2;
        Double_t fQ1;   //! charges
        Double_t fQ2;
        // Vertex info:
        Double_t fVertexZ;
        Int_t    fVertexContrib;
        // Info from the detectors:
        // ZDC
        Double_t fZNA_energy;
        Double_t fZNC_energy;
        // TDC = Time-to-Digital Converter
        Double_t fZNA_time[4];
        Double_t fZNC_time[4];
        // V0:
        Int_t fV0A_dec;
        Int_t fV0C_dec;
        // AD:
        Int_t fADA_dec;
        Int_t fADC_dec;
        // Matching SPD clusters with FOhits
        Bool_t fMatchingSPD;
        TBits fFOCrossFiredChips;
        // Trigger inputs for MC data
        Bool_t  fTriggerInputsMC[11];
        TFile   *fSPDfile;
        TFile   *fTOFfile;
        Int_t   fLoadedRun;
        TH2F    *hTOFeff;
        TH1D    *hSPDeff;
        AliTOFTriggerMask *fTOFmask;
        // MC kinematics on generated level
        Double_t fPtGen;
        Double_t fYGen;
        Double_t fMGen;
        Double_t fPhiGen;

        AliAnalysisTaskCentralJpsi_DG(const AliAnalysisTaskCentralJpsi_DG&); // not implemented
        AliAnalysisTaskCentralJpsi_DG& operator=(const AliAnalysisTaskCentralJpsi_DG&); // not implemented

        ClassDef(AliAnalysisTaskCentralJpsi_DG, 1);
};

#endif
