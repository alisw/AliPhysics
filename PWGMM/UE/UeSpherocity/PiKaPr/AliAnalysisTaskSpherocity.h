#ifndef AliAnalysisTaskSpherocity_H
#define AliAnalysisTaskSpherocity_H
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
//#include <AliSpherocityUtils.h>
#include <AliEventCuts.h>
#include "AliVTrack.h"
#include <vector>
using namespace std;



class AliAnalysisTaskSpherocity : public AliAnalysisTaskSE
{
public:
    
    
    AliAnalysisTaskSpherocity();
    AliAnalysisTaskSpherocity(const char *name);
    virtual ~AliAnalysisTaskSpherocity();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    
    Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
    Double_t GetEtaCut() { return fEtaCut; }     
    
    virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
    virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
    virtual void  SetTrackFilter2015PbPb(AliAnalysisFilter* trackF) {fTrackFilter2015PbPb = trackF;}
    virtual void  SetTrackCutsSpherocity(AliAnalysisFilter* fTrackFilter);
    virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
    virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
    virtual void  SetNcl(const Int_t ncl){fNcl = ncl;}
    virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
    virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
    virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
    virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) { fAnalysisPbPb = isanaPbPb; }
    virtual void  SetAnalysisTask(Bool_t PostCalib) { fdEdxCalibrated = PostCalib; }
    virtual void  SetAnalysisPID(Bool_t makePid) { fMakePid = makePid; }
    //		virtual void  SetAddLowPt(Bool_t addlowpt) { fLowPt = addlowpt; }
    virtual void  SetPeriod(Int_t isLHC16l) { fLHC16l = isLHC16l; }
    virtual void  SetEstimator(const Bool_t isV0M) { fisV0Mestimator = isV0M; }
    virtual void  SetJettyCutOff(const Double_t JettyCutOff) { fJettyCutOff = JettyCutOff; }
    virtual void  SetIsotrCutOff(const Double_t IsotrCutOff) { fIsotrCutOff = IsotrCutOff; }
    
private:
    
    virtual void ProduceArrayTrksESD(AliESDEvent* event, const Int_t cent, const Int_t sperocity);
    virtual void ProduceArrayV0ESD(AliESDEvent* event, const Int_t cent, const Int_t sperocity );
    Int_t   GetV0MIndex(Float_t V0MPercentile);
    Int_t   GetMultiplicityIndex(Int_t Mult);
    ////    void    PtRecVsPtTruth(AliESDEvent* event, const Bool_t isjetty);
    Int_t   GetMultiplicityParticles(Double_t etaCut);
    Short_t GetPidCode(Int_t pdgCode) const;
    void    ProcessMCTruthESD( const double cent, const int so );
    
    TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
    Int_t      FindPrimaryMotherLabel(AliStack* stack, Int_t label);
    
    
    TParticle* FindPrimaryMotherV0(AliStack* stack, Int_t label);
    int      FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps);
    bool selectVertex2015pp(AliESDEvent* esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
    bool IsGoodSPDvertexRes(const AliESDVertex* spdVertex = NULL);
    bool IsGoodZvertexPos(AliESDEvent *esd);
    bool PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh);
    float GetMaxDCApTDep(TF1* fcut, Double_t pt );
    float GetSpherocity( TH1D * hphi, TH1D *heta );
    int ReadESDEvent( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TH1D * hphi, TH1D *heta );
    float AnalyseGetSpherocity( const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi );
    bool TOFPID(AliESDtrack* track);
    double EtaCalibrationNeg(const Int_t centrality, const Double_t Eta);
    double EtaCalibrationPos(const Int_t centrality, const Double_t Eta);
    double EtaCalibrationNegEl(const Int_t centrality, const Double_t Eta);
    double EtaCalibrationPosEl(const Int_t centrality, const Double_t Eta);
    
    
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
    //		AliSpherocityUtils* fSpheroUtils;
    TString       fAnalysisType;        //  "ESD" or "AOD"
    Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
    Bool_t        fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
    Bool_t        fisV0Mestimator; 
    TRandom*      fRandom;              //! random number generator
    
    //
    // Cuts and options
    //
    
    Int_t        fNcl;                
    Double_t     fEtaCut;             // Eta cut used to select particles
    Int_t        fCent; //minimum centrality
    Int_t        fSpherocity;
    Float_t      fMinCent; //minimum centrality
    Float_t      fMaxCent; //maximum centrality
    const Double_t fDeDxMIPMin;
    const Double_t fDeDxMIPMax;
    int fdEdxHigh;
    int fdEdxLow;
    Double_t fJettyCutOff;
    Double_t fIsotrCutOff;
    int fMinMult;
    int fNrec;
    float fSizeStep;
    
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
    TH2F*         fEvents;            //! No of accepted events
    TH1F* 	hMult;
    Bool_t       fdEdxCalibrated;
    Bool_t       fMakePid;
    Int_t  fLHC16l;
    TH1F* fcent;
    TH1F* fcentAfterPrimaries;
    TH1F* fcentAfterV0s;
    
    TH1D *hphiso;
    TH1D *hetaso;
    TH2D *hSOrvsV0M;
    TH2D *hSOrvsTrks;
    TH2D *hRefMultVsRefMultPer;
    
    // Histograms for PreCalibration
    
    TH2D *hMIPVsEta[11][3];
    TProfile *pMIPVsEta[11][3];
    TH2D *hMIPVsEtaV0s[11][3];
    TProfile *pMIPVsEtaV0s[11][3];
    TH2D *hPlateauVsEta[11][3];
    TProfile *pPlateauVsEta[11][3];
    TH2D *hPhi[11];
    
    TH2D *hMIPVsV0M[4];
    TProfile *pMIPVsV0M[4];
    TH2D *hMIPVsNch[4];
    TProfile *pMIPVsNch[4];
    
    TH2D     *hMIPVsPhi[11][4][3];
    TProfile *pMIPVsPhi[11][4][3];
    TH2D     *hPlateauVsPhi[11][4][3];
    TProfile *pPlateauVsPhi[11][4][3];
    
    
    // Histograms for PostCalibration
    
    
    TH2D *hPtAll[3];
    TH3D *hPtVsP[3][4];
    TH3D *hDeDxVsP[3][4];
    
    TH3D *hnSigmaPiPos[3][4];
    TH3D *hnSigmaKPos[3][4];
    TH3D *hnSigmaPPos[3][4];
    TH3D *hnSigmaPiNeg[3][4];
    TH3D *hnSigmaKNeg[3][4];
    TH3D *hnSigmaPNeg[3][4];
    
    TH3D *hBetavsPneg[3][4];
    TH3D *hBetavsPpos[3][4];
    
    TH2D *hPtneg_TPC_Eta[3][4];
    TH2D *hPtpos_TPC_Eta[3][4];
    TH2D *hPtneg_TPC[3];
    TH2D *hPtpos_TPC[3];
    TH2D *hPtneg_TOF_Eta[3][4];
    TH2D *hPtpos_TOF_Eta[3][4];
    TH2D *hPneg_TOF_Eta[3][4];
    TH2D *hPpos_TOF_Eta[3][4];
    TH2D *hPtneg_TOF[3];
    TH2D *hPtpos_TOF[3];
    
    TH2D* histPiV0[11][4][3];
    TH2D* histPV0[11][4][3];
    TH2D* histPiTof[11][4][3];
    TH2D* histEV0[11][4][3];
    
    TH3D* hDCApi_TPC[3][2];
    TH3D* hDCAp_TPC[3][2];
    TH3D* hDCApi2_TPC[3][2];
    TH3D* hDCAp2_TPC[3][2];
    
    TH3D* hDCApi_TOF[3][2];
    TH3D* hDCAp_TOF[3][2];
    TH3D* hDCApi2_TOF[3][2];
    TH3D* hDCAp2_TOF[3][2];
    
    TH2D* hMcIn[3][7];
    TH2D* hMcOut[3][7];
    TH2D* hMcInNeg[3][7];
    TH2D* hMcInPos[3][7];
    TH2D* hMcOutNeg[3][7];
    TH2D* hMcOutPos[3][7];
    
    TH3D* hDCAPrim_pi_TPC[3][2];
    TH3D* hDCAWDec_pi_TPC[3][2];
    TH3D* hDCAMate_pi_TPC[3][2];
    TH3D* hDCAPrim2_pi_TPC[3][2];
    TH3D* hDCAWDec2_pi_TPC[3][2];
    TH3D* hDCAMate2_pi_TPC[3][2];
    
    TH3D* hDCAPrim_p_TPC[3][2];
    TH3D* hDCAWDec_p_TPC[3][2];
    TH3D* hDCAMate_p_TPC[3][2];
    TH3D* hDCAPrim2_p_TPC[3][2];
    TH3D* hDCAWDec2_p_TPC[3][2];
    TH3D* hDCAMate2_p_TPC[3][2];
    
    TH3D* hDCAPrim_pi_TOF[3][2];
    TH3D* hDCAWDec_pi_TOF[3][2];
    TH3D* hDCAMate_pi_TOF[3][2];
    TH3D* hDCAPrim2_pi_TOF[3][2];
    TH3D* hDCAWDec2_pi_TOF[3][2];
    TH3D* hDCAMate2_pi_TOF[3][2];
    
    TH3D* hDCAPrim_p_TOF[3][2];
    TH3D* hDCAWDec_p_TOF[3][2];
    TH3D* hDCAMate_p_TOF[3][2];
    TH3D* hDCAPrim2_p_TOF[3][2];
    TH3D* hDCAWDec2_p_TOF[3][2];
    TH3D* hDCAMate2_p_TOF[3][2];
    
    TF1* fEtaCalibrationNeg;
    TF1* fEtaCalibration;
    TF1* felededxfitPos;
    TF1* felededxfitNeg;
    TF1* fcutDCAxy;
    TF1* fcutLow;
    TF1* fcutHigh;
    
    AliAnalysisTaskSpherocity(const AliAnalysisTaskSpherocity&);            // not implemented
    AliAnalysisTaskSpherocity& operator=(const AliAnalysisTaskSpherocity&); // not implemented
    ClassDef(AliAnalysisTaskSpherocity, 1);    //Analysis task for high pt analysis 
};

#endif


