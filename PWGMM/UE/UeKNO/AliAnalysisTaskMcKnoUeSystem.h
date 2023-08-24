/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskMcKnoUeSystem_H
#define AliAnalysisTaskMcKnoUeSystem_H

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



class AliAnalysisTaskMcKnoUeSystem : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskMcKnoUeSystem();
	AliAnalysisTaskMcKnoUeSystem(const char *name);
	
    virtual                 ~AliAnalysisTaskMcKnoUeSystem();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);
    
    //Getter
	void       GetLeadingObject(Bool_t isMC);
	void       GetDetectorResponse();
    void       GetTrackingEfficiencyTPConly();
	void       GetBinByBinCorrections();
	void       GetUEObservables();
    void       GetUEObservablesData();
	void       GetPtLeadingMisRecCorrection();
	void       GetMeanUEObservables(std::vector<Double_t> &gen, std::vector<Double_t> &rec);
	void       GetMultiplicityDistributions();
	void       GetMultiplicityDistributionsData();
    
    //Setter
	void       SetPtMin(Double_t val)                  {fPtMin = val;}        // use differnet ptcuts
	void       SetUseMC(Bool_t mc = kFALSE)            {fUseMC = mc;}   // use to analyse MC data
	void       SetMCclosureTest(Bool_t mcc = kFALSE)   {fIsMCclosure = mcc;}
    void       SetParametrizationEfficiency(Bool_t ispy = kFALSE)          {fIsPythia = ispy;}
	void       SetParametrizationEfficiencyppdata(Bool_t ispp = kFALSE)    {fIsppData = ispp;}
    void       SetParametrizationEfficiencypPbdata(Bool_t ispPb = kFALSE)  {fIspPbData = ispPb;}
    void       SetLeadingPtMin(Double_t PtLmin)    {fLeadPtCutMin = PtLmin;}   // use differnet ptcuts
    void       SetLeadingPtMax(Double_t PtLmax)    {fLeadPtCutMax = PtLmax;}   // use differnet ptcuts
    
    //systematic uncertainties --- variation
    //event selection for vertex position
    void    SetUpVertexZposition(Float_t val = 10.0) { fCutVertexZposition = val;}
    //track selection for leading particle
    void    SetUpMaxFractionSharedTPCClusters(Float_t max = 1e10)              {fCutMaxFractionSharedTPCClusters = max;}
    void    SetUpMinRatioCrossedRowsOverFindableClustersTPC(Float_t min = -1)  {fCutMinRatioCrossedRowsOverFindableClustersTPC = min;}
    void    SetUpCutGeoNcrNcl(Float_t deadZoneWidth = 3.0, Float_t cutGeoNcrNclLength = 130.0)  {fCutGeoNcrNclZone = deadZoneWidth; fCutGeoNcrNclLength = cutGeoNcrNclLength;}
    void    SetUpClusterRequirementITS(Bool_t req = kFALSE)   {fIsRequirementSPD = req;}
    void    SetUpMaxChi2PerClusterITS(Float_t max = 1e10)        {fCutMaxChi2PerClusterITS = max;}
    void    SetUpMaxChi2PerClusterTPC(Float_t max = 1e10)        {fCutMaxChi2PerClusterTPC = max;}
    void    SetUpMaxChi2TPCConstrainedGlobal(Float_t max = 1e10) {fCutMaxChi2TPCConstrainedVsGlobal = max;}
    void    SetUpMaxDCAToVertexZ(Float_t dist = 1e10)            {fCutMaxDCAToVertexZ = dist;}
    //track selection for Nch
    void    SetUpMinNClustersTPC(Int_t min = -1)          {fCutMinNClusterTPC = min;}
    void    SetUpMaxChi2PerClusterTPC_nch(Float_t max = 1e10) {fCutMaxChi2PerClusterTPC_nch = max;}
    void    SetUpMaxDCAToVertexZ_nch(Float_t dist = 1e10)  {fCutMaxDCAToVertexZ_nch = dist;}
    void    SetUpMaxDCAToVertexXY(Float_t dist = 1e10)    {fCutMaxDCAToVertexXY = dist;}
    void    SetUpRequireITSRefit(Bool_t req = kFALSE)   {fIsITSRefit = req;}
    void    SetUpRequireTPCRefit(Bool_t req = kFALSE)   {fIsTPCRefit = req;}
    //event
    
    
	bool       HasRecVertex();
	virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );

protected:



private:
	AliESDEvent* fESD;                                        //! input ESD event
	AliMCEvent*  fMC;                                               //! MC Event
	AliStack*    fMCStack;                                                 //! MC stack
	Bool_t       fUseMC;                // analyze MC events
	Bool_t       fIsMCclosure;         
	Bool_t       fIsPythia;
        Bool_t       fIsppData;
        Bool_t       fIspPbData;
	AliEventCuts fEventCuts;
	AliAnalysisFilter*  fLeadingTrackFilter;
	AliAnalysisFilter*  fTrackFilterForDCA;
	AliAnalysisFilter*  fTrackFilter;
	TList*              fOutputList;                                      //! output list in the root file

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
        
	AliMultSelection *fMultSelection;
	Double_t fRefmult08std;
        Double_t fpercentileV0M;
    
    
    //systematic uncertainties --- variation
    //event selection for vertex position
    Float_t  fCutVertexZposition; //10(default), 5, 15
    //track selection for leading particle
    Float_t  fCutMaxFractionSharedTPCClusters;                //0.4(default), 0.2, 1.0
    Float_t  fCutMinRatioCrossedRowsOverFindableClustersTPC;  //0.8(default), 0.7, 0.9
    Float_t  fCutGeoNcrNclZone, fCutGeoNcrNclLength;          //3cm(default), 2cm, 4cm;     130(default), 120, 140
    Bool_t   fIsRequirementSPD;         //kTRUE--kAny(default), kFALSE--don't use that cut
    Float_t  fCutMaxChi2PerClusterITS;  //36(default), 25, 49
    Float_t  fCutMaxChi2PerClusterTPC;  //4(default), 3, 5
    Float_t  fCutMaxChi2TPCConstrainedVsGlobal;  //36(default), 25, 49
    Float_t  fCutMaxDCAToVertexZ;                //2(default), 1, 5
    //track selection for Nch
    Int_t    fCutMinNClusterTPC;          //50(default), 30, 70
    Float_t  fCutMaxChi2PerClusterTPC_nch; //4(default),  3, 5
    Float_t  fCutMaxDCAToVertexZ_nch;      //3.2(default), 2.0, 4.0
    Float_t  fCutMaxDCAToVertexXY;        //2.4(default), 1.0, 4.0
    Bool_t   fIsITSRefit; //kTRUE(default), kFALSE    
    Bool_t   fIsTPCRefit; //kTRUE(default), kFALSE
    

    
    //histograms
    TH1I * hCounter;
    //vertex Z position
    TH1D * hZvtxAllMeasured;
    TH1D * hZvtxGoodVtxMeasured;
    TH1D * hZvtxCutAccMeasured;
    TH1D * hZvtxAllGen;
    TH1D * hZvtxCutGen;
    TH1D * hZvtxTrigGen;
    TH1D * hZvtxGoodVtxGen;
    TH1D * hZvtxCutAccGen;
    TH1D * hZvtxCutKnoGen;
    // KNO
    TH1D * hNchTSData;
    TH1D * hNchTSminData;
    TH1D * hNchTSmaxData;
    TH1D * hPhiData_TS1;
    TH1D * hPhiData_TS2;

    TH1D * hPhiGen[3];
    TH1D * hPhiRec[3];
    TH1D * hNchTSGen;
    TH1D * hNchTSRec;
    TH2D * hNchTSResponse;
    TH1D * hNchTSGenTest;
    TH1D * hNchTSRecTest;
         
    TH1D * hPhiGen_TS1;
    TH1D * hPhiGen_TS2;
    TH1D * hPhiRec_TS1;
    TH1D * hPhiRec_TS2;
    TH1D * hPhiGenTest_TS1;
    TH1D * hPhiGenTest_TS2;
    TH1D * hPhiRecTest_TS1;
    TH1D * hPhiRecTest_TS2;
        
    TH1D * hNchTSminGen;
    TH1D * hNchTSminRec;
    TH2D * hNchTSminResponse;
    TH1D * hNchTSminGenTest;
    TH1D * hNchTSminRecTest;
        
    TH1D * hNchTSmaxGen;
    TH1D * hNchTSmaxRec;
    TH2D * hNchTSmaxResponse;
    TH1D * hNchTSmaxGenTest;
    TH1D * hNchTSmaxRecTest;
    
    TH1D * hPtPrimGen;
    TH1D * hPtPrimRec;
    


        // DCA
        TH2D * hPTVsDCAData;
        TH2D * hPtDCAPrimary;
    	TH2D * hPtDCAWeak;
    	TH2D * hPtDCAMat;
    	TH2D * hPtDCAall;
	// UE 
	TH1D * hPtInPrim;
	TH1D * hPtOut;
	TH1D * hPtOutPrim; 
	TH1D * hPtOutSec;
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
        TH1D * hRefMult08std;
        TH1D * hMultV0M;
        TH2D * hRefMultvsMultV0M;

	AliAnalysisTaskMcKnoUeSystem(const AliAnalysisTaskMcKnoUeSystem&);                  // not implemented
	AliAnalysisTaskMcKnoUeSystem& operator=(const AliAnalysisTaskMcKnoUeSystem&);       // not implemented

	ClassDef(AliAnalysisTaskMcKnoUeSystem, 3);
};

#endif
