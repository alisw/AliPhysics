#ifndef AliAnalysisTaskPPvsRT_TPCTOF_H
#define AliAnalysisTaskPPvsRT_TPCTOF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */


// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH3.h>
#include <TProfile.h>
#include <TTreeStream.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 
#include <AliESDtrackCuts.h>
#include <AliPIDResponse.h>
#include "AliTPCPIDResponse.h"
#include <AliEventCuts.h>
#include "AliVTrack.h"



class AliAnalysisTaskPPvsRT_TPCTOF : public AliAnalysisTaskSE
{
public:
    
    
    AliAnalysisTaskPPvsRT_TPCTOF();
    AliAnalysisTaskPPvsRT_TPCTOF(const char *name);
    virtual ~AliAnalysisTaskPPvsRT_TPCTOF();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    
    Bool_t   GetAnalysisMC() { return fAnalysisMC; }
    Double_t GetEtaCut() { return fEtaCut; }
    
    virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
    virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
    virtual void  SetTrackFilter2015PbPb(AliAnalysisFilter* trackF) {fTrackFilter2015PbPb = trackF;}
    virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
    virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
    virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
    virtual void  SetLeadingCut(const double LeadingCut){fLeadingCut = LeadingCut;}
    virtual void  SetNcl(const Int_t ncl){fNcl = ncl;}
    virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
    virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
    virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
    virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) { fAnalysisPbPb = isanaPbPb; }
    virtual void  SetPeriod(const char* Period) { fPeriod = Period; }
    virtual void  SetMeanCh(const Double_t MeanCh) { fMeanChT = MeanCh; }
    
private:
    
    
    AliESDtrack* GetLeadingTrack();
    TObjArray* SortRegions(AliESDtrack* Ltrk);
    virtual Double_t DeltaPhi(Double_t phi, Double_t lphi,
                              Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
    int GetBinRT(TList* lt);
    virtual void ProduceArrayTrksESD(const int& r, TList* lt, const int& cent);
    short   GetPidCode(Int_t pdgCode) const;
    
    bool selectVertex2015pp(AliESDEvent* esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
    bool IsGoodSPDvertexRes(const AliESDVertex* spdVertex = NULL);
    bool IsGoodZvertexPos(AliESDEvent *esd);
    bool PhiCut(const double& pt, double phi, const double& q, const float& mag, TF1* phiCutLow, TF1* phiCutHigh);
    float GetMaxDCApTDep( TF1 *fcut, Double_t pt );
    virtual void SetTrackCuts(AliAnalysisFilter* fTrackFilter);
    double EtaCalibration(const int &centrality, const double &Eta);
    double EtaCalibrationEl(const int &centrality, const double &Eta);
    int GetIndex();
    bool TOFPID(AliESDtrack* track);
    
    static const Double_t fgkClight;   // Speed of light (cm/ps)
    
    AliESDEvent* fESD;                  //! ESD object
    AliAODEvent* fAOD;                  //! AOD object
    AliEventCuts fEventCuts;
    AliMCEvent*  fMC;                   //! MC object
    AliStack*    fMCStack;              //! MC ESD stack
    TClonesArray* fMCArray;             //! MC array for AOD
    AliPIDResponse* fPIDResponse;       //! Pointer to PIDResponse
    AliAnalysisFilter* fTrackFilter2015PbPb;    //  Track Filter, set 2010 with golden cuts
    AliAnalysisFilter* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
    AliAnalysisFilter* fTrackFilterTPC; // track filter for TPC only tracks
    AliAnalysisFilter* fTrackFilter;
    AliAnalysisUtils* utils;
    TString       fAnalysisType;        //  "ESD" or "AOD"
    Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
    Bool_t        fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
    TRandom*      fRandom;              //! random number generator
    
    //
    // Cuts and options
    //
    
    Double_t     fVtxCut;             // Vtx cut on z position in cm
    double       fLeadingCut;
    Int_t        fNcl;
    Double_t     fEtaCut;             // Eta cut used to select particles
    Int_t        cent; //minimum centrality
    Float_t      fMinCent; //minimum centrality
    Float_t      fMaxCent; //maximum centrality
    const Double_t fDeDxMIPMin;
    const Double_t fDeDxMIPMax;
    const Double_t fdEdxHigh;
    const Double_t fdEdxLow;
    TString  fPeriod;
    Double_t fMeanChT;
    
    //
    // Help variables
    //
    Short_t      fMcProcessType;      // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
    Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only)
    Short_t      fVtxStatus;          // -1 = no vtx, 0 = outside cut, 1 = inside cut
    Float_t      fZvtx;               // z vertex
    Float_t      fZvtxMC;             // z vertex MC (truth)
    Int_t        fRun;                // run no
    ULong64_t    fEventId;            // unique event id
    
    //
    // Output objects
    //
    TList*        fListOfObjects;     //! Output list of objects
    TH1F*         fEvents;            //! No of accepted events
    TH1F*         fVtxMC;             //! Event vertex info for ALL MC events
    Bool_t       fdEdxCalibrated;
    Bool_t       fMakePid;
    Bool_t       fLowPt;
    TH1F* fMultN;
    TH1F* fPtN;
    TH1F* fDphiN;
    TH1F* fMultA;
    TH1F* fPtA;
    TH1F* fDphiA;
    TH1F* fMultT;
    TH1F* fPtT;
    TH1F* fMultNMC;
    TH1F* fMultAMC;
    TH1F* fMultTMC;
    TH1F* fDphiT;
    TH1F* fDphi;
    TH1D* fRT;
    TH1D* fRTMC;
    TH2D* fPtLVsRt; 
    TH2D* fPtLVsRtMC; 
    
    // [rt][region][eta]
    TH2D *hPtVsP[4];
    TH2D *hBetavsPMB[4];
    TH1D *hPtpos_TPC[5][3];
    TH1D *hPtpos_TPC_Eta[5][3][4];
    TH1D *hPtneg_TPC[5][3];
    TH1D *hPtneg_TPC_Eta[5][3][4];
    TH1D *hPtpos_TOF[5][3];
    TH1D *hPtpos_TOF_Eta[5][3][4];
    TH1D *hPpos_TOF_Eta[5][3][4];
    TH1D *hPtneg_TOF[5][3];
    TH1D *hPtneg_TOF_Eta[5][3][4];
    TH1D *hPneg_TOF_Eta[5][3][4];
    TH2D *hnSigPipos[5][3][4];
    TH2D *hnSigkpos[5][3][4];
    TH2D *hnSigppos[5][3][4];
    TH2D *hnSigPineg[5][3][4];
    TH2D *hnSigkneg[5][3][4];
    TH2D *hnSigpneg[5][3][4];
    TH2D *hBetavsPpos[5][3][4];
    TH2D *hBetavsPneg[5][3][4];
    
    TF1* fEtaCalibration;
    TF1* fEtaCalibrationEl;
    TF1* fcutDCAxy;
    TF1* fcutLow;
    TF1* fcutHigh;
    
    
    AliAnalysisTaskPPvsRT_TPCTOF(const AliAnalysisTaskPPvsRT_TPCTOF&);            // not implemented
    AliAnalysisTaskPPvsRT_TPCTOF& operator=(const AliAnalysisTaskPPvsRT_TPCTOF&); // not implemented
    
    //TTree*        fTree;              //! Debug tree
    
    ClassDef(AliAnalysisTaskPPvsRT_TPCTOF, 1);    //Analysis task for high pt analysis
};

#endif

