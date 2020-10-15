/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskMcKnoUeSyst_H
#define AliAnalysisTaskMcKnoUeSyst_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1D;
class TH2D;
class TH3D;
class TH1I;
class TProfile;
class THnSparse;

#include "AliAnalysisTaskSE.h"


#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliGenEventHeader.h"



class AliAnalysisTaskMcKnoUeSyst : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskMcKnoUeSyst();
	AliAnalysisTaskMcKnoUeSyst(const char *name);
	
    virtual                 ~AliAnalysisTaskMcKnoUeSyst();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);
    
	void       GetLeadingObject(Bool_t isMC);
	void       GetDetectorResponse();
	void       GetBinByBinCorrections();
	void       GetUEObservables();
    void       GetUEObservablesData();
	void       GetPtLeadingMisRecCorrection();
	void       GetMeanUEObservables(std::vector<Double_t> &gen, std::vector<Double_t> &rec);
	void       GetMultiplicityDistributions();
	void       GetMultiplicityDistributionsData();
	void       SetPtMin(Double_t val)              {fPtMin = val;}   // use differnet ptcuts
	void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}   // use to analyse MC data
	void       SetMCclosureTest(Bool_t mcc = kFALSE)    {fIsMCclosure = mcc;}
	void       SetParametrizationEfficiency(Bool_t ispy = kTRUE)  {fIsPythia = ispy;}
	bool       HasRecVertex();
    //Systematic ============================
    void       SetTPCclustersLow(Bool_t TPCclustersLow = kFALSE) {fTPCclustersLow = TPCclustersLow;}
        void       SetTPCclustersHigh(Bool_t TPCclustersHigh = kFALSE) {fTPCclustersHigh = TPCclustersHigh;}
        void       SetNcrLow(Bool_t NcrLow = kFALSE) {fNcrLow = NcrLow;}
        void       SetNcrHigh(Bool_t NcrHigh = kFALSE) {fNcrHigh = NcrHigh;}
        void       SetChisqTPCLow(Bool_t ChisqTPCLow = kFALSE) {fChisqTPCLow = ChisqTPCLow;}
        void       SetChisqTPCHigh(Bool_t ChisqTPCHigh = kFALSE) {fChisqTPCHigh = ChisqTPCHigh;}
        void       SetChisqITSLow(Bool_t ChisqITSLow = kFALSE) {fChisqITSLow = ChisqITSLow;}
        void       SetChisqITSHigh(Bool_t ChisqITSHigh = kFALSE) {fChisqITSHigh = ChisqITSHigh;}
        void       SetChisqITSmTPCLow(Bool_t ChisqITSmTPCLow = kFALSE) {fChisqITSmTPCLow = ChisqITSmTPCLow;}
        void       SetChisqITSmTPCHigh(Bool_t ChisqITSmTPCHigh = kFALSE) {fChisqITSmTPCHigh = ChisqITSmTPCHigh;}
        void       SetDcazLow(Bool_t DcazLow = kFALSE) {fDcazLow = DcazLow;}
        void       SetDcazHigh(Bool_t DcazHigh = kFALSE) {fDcazHigh = DcazHigh;}
        void       SetGeoTPCLow1(Bool_t GeoTPCLow1 = kFALSE) {fGeoTPCLow1 = GeoTPCLow1;}
        void       SetGeoTPCLow2(Bool_t GeoTPCLow2 = kFALSE) {fGeoTPCLow2 = GeoTPCLow2;}
        void       SetGeoTPCHigh1(Bool_t GeoTPCHigh1 = kFALSE) {fGeoTPCHigh1 = GeoTPCHigh1;}
        void       SetGeoTPCHigh2(Bool_t GeoTPCHigh2 = kFALSE) {fGeoTPCHigh2 = GeoTPCHigh2;}
    void       SetSPDreqVar1(Bool_t SPDreqVar1 = kFALSE) {fSPDreqVar1 = SPDreqVar1;}
    void       SetVertexZCutLow(Bool_t VertexZCutLow = kFALSE) {fVertexZCutLow = VertexZCutLow;}
    void       SetVertexZCutHigh(Bool_t VertexZCutHigh = kFALSE) {fVertexZCutHigh = VertexZCutHigh;}
        //Systematic ============================
	virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );

protected:



private:
	AliESDEvent* fESD;                                        //! input ESD event
	Bool_t       fIsPythia; 
	AliEventCuts fEventCuts;
	AliStack*    fMCStack;                                                 //! MC stack
	AliMCEvent*  fMC;                                               //! MC Event
	Bool_t       fUseMC;                // analyze MC events
	Bool_t       fIsMCclosure;         
	AliAnalysisFilter*  fLeadingTrackFilter;
	AliAnalysisFilter*  fTrackFilterForDCA;
	AliAnalysisFilter*  fTrackFilter;
	TList*                  fOutputList;                                      //! output list in the root file

	Double_t fEtaCut;
	Double_t fPtMin;
	Double_t fLeadPtCutMin;
	Double_t fLeadPtCutMax;
	Double_t fGenLeadPhi; 
	Double_t fGenLeadPt;
	Int_t    fGenLeadIn;
	Double_t fRecLeadPhi; 
	Double_t fRecLeadPt;
	Int_t    fRecLeadIn;

	Float_t fDCAxy;
	Float_t fDCAz;	
    
    Bool_t  fTPCclustersLow;
    Bool_t  fTPCclustersHigh;
    Bool_t  fNcrLow;
    Bool_t  fNcrHigh;
    Bool_t  fChisqTPCLow;
    Bool_t  fChisqTPCHigh;
    Bool_t  fChisqITSLow;
    Bool_t  fChisqITSHigh;
    Bool_t  fChisqITSmTPCLow;
    Bool_t  fChisqITSmTPCHigh;
    Bool_t  fDcazLow;
    Bool_t  fDcazHigh;
    Bool_t  fGeoTPCLow1;
    Bool_t  fGeoTPCLow2;
    Bool_t  fGeoTPCHigh1;
    Bool_t  fGeoTPCHigh2;
    Bool_t  fSPDreqVar1;
    Bool_t  fVertexZCutLow;
    Bool_t  fVertexZCutHigh;
//     // Corrections
//     
//     TH2D * hMCFractions;    // ! Histos for particle abundances from MC
//     
    
    // DCA 
    TH2D * hPTVsDCAData;
    TH2D * hPtDCAPrimary;
    TH2D * hPtDCAWeak;
    TH2D * hPtDCAMat;
    TH2D * hPtDCAall;    
    
	// KNO
	TH1D * hPhiGen[3];
	TH1D * hNchTSGen;
	TH1D * hNchTSGenTest;
	TH1D * hPhiRec[3];
	TH1D * hNchTSRec;
	TH1D * hNchTSRecTest;
	TH1D * hNchTSData;
	TH2D * hNchResponse;

	// UE 
	TH1D * hPtInPrim;
	TH1D * hPtOut;
	TH1D * hPtOutPrim; 
	TH1D * hPtOutSec; 
	TH1D * hCounter;
	TH2D * hNumDenMC[3];
	TH2D * hSumPtMC[3];
	TH2D * hNumDenMCMatch[3];
	TH2D * hSumPtMCMatch[3];
	TH2D * hNumDenMCDd[3];
	TH2D * hSumPtMCDd[3];
	TH2D * hNumDenMCMatchDd[3];
	TH2D * hSumPtMCMatchDd[3];

	TH1D * hPtLeadingTrue;
	TH1D * hPtLeadingMeasured;
	TH1D * hPtLeadingData;
	TH2D * hPtVsPtLeadingMeasured[3];
	TH2D * hPtVsPtLeadingData[3];
	TH2D * hPtVsPtLeadingTrue[3];
	TProfile * pNumDenMeasured[3];
	TProfile * pNumDenData[3];
	TProfile * pNumDenTrue[3];
	TProfile * pSumPtMeasured[3];
	TProfile * pSumPtData[3];
	TProfile * pSumPtTrue[3];

	TProfile * pNumDenMeasuredAll[3];
	TProfile * pNumDenTrueAll[3];
	TProfile * pSumPtMeasuredAll[3];
	TProfile * pSumPtTrueAll[3];

	TProfile * pNumDenMeasuredPS[3];
	TProfile * pNumDenTruePS[3];
	TProfile * pSumPtMeasuredPS[3];
	TProfile * pSumPtTruePS[3];

	TProfile * pNumDenMeasuredPSV[3];
	TProfile * pNumDenTruePSV[3];
	TProfile * pSumPtMeasuredPSV[3];
	TProfile * pSumPtTruePSV[3];

	TProfile * pNumDenMeasuredGood[3];
	TProfile * pNumDenTrueGood[3];
	TProfile * pSumPtMeasuredGood[3];
	TProfile * pSumPtTrueGood[3];

	TH2D * hPtVsUEGenTest[3];
	TH2D * hPtVsUERecTest[3];
	TH2D * hPtVsUEData[3];

	TH1D * hPtInPrimPart[6];
	TH1D * hPtOutPrimPart[6];
    
	TH1D * hPtLeadingRecPS;
	TH1D * hPtLeadingRecPSV;
	TH1D * hPtLeadingRecGood;
	TH1D * hPtLeadingGenPS;
	TH1D * hPtLeadingGenPSV;
	TH1D * hPtLeadingGenGood;
	TH1D * hPtLeadingRecAll;
	TH1D * hPtLeadingGenAll;

	AliAnalysisTaskMcKnoUeSyst(const AliAnalysisTaskMcKnoUeSyst&);                  // not implemented
	AliAnalysisTaskMcKnoUeSyst& operator=(const AliAnalysisTaskMcKnoUeSyst&);       // not implemented

	ClassDef(AliAnalysisTaskMcKnoUeSyst, 3);
};

#endif
