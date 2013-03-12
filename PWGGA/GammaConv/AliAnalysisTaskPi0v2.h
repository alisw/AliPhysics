#ifndef AliAnalysisTaskPi0v2_cxx
#define AliAnalysisTaskPi0v2_cxx

#include "AliAnalysisTaskSE.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "AliLog.h"
#include "AliConversionSelection.h"
#include "AliConversionMesonCuts.h"
#include "AliV0ReaderV1.h"
#include "AliEventplane.h"
#include "TVector2.h"

#include "TProfile.h"
using namespace std;

class AliAnalysisTaskPi0v2 : public AliAnalysisTaskSE{

public:

    enum EEventPlaneMethod{
	kTPC=0,
        kTPCEtaGap=1,
	kV0A=2,
	kV0C=3,
	knEPMethod=4
    };

    enum EEventPlane{
	kEPTPC=0,
	kEPTPCEtaA,
	kEPTPCEtaC,
	kEPV0A,
	kEPV0C,
        knEP
    };

    enum EPDGCode{
	kPi0=111,
	kEta=221
    };

    static const Int_t knBinsPhi=6;

    AliAnalysisTaskPi0v2(const char *name="pi0v2",Int_t harmonic=2);
    AliAnalysisTaskPi0v2(const AliAnalysisTaskPi0v2&); // not implemented
    AliAnalysisTaskPi0v2& operator=(const AliAnalysisTaskPi0v2&); // not implemented
    virtual ~AliAnalysisTaskPi0v2();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetCentralityBins(Double_t *bins,Int_t nbins);

    void SetMeson(EPDGCode meson){fMesonPDGCode=meson;}

    void SetNBinsPhi(Int_t nbins){fNBinsPhi=nbins;}
    void SetV0Reader(AliV0ReaderV1 *v0Reader){fV0Reader=v0Reader;}
    void SetInvMassRange(Double_t range[2]){fInvMassRange[0]=range[0];fInvMassRange[1]=range[1];};
    void SetEtaGap(Double_t gapsize){fEtaGap=gapsize;};

    void SetMesonCuts(const TString cut);
    void SetCuts(AliConversionSelection **conversionselection,Int_t numberOfCuts);

    void SetFillQA(Bool_t fill){fFillQA=fill;}

    void SetEPSelectionMask(Int_t mask[knEPMethod]){for(Int_t ii=0;ii<knEPMethod;ii++)fEPSelectionMask[ii]=mask[ii];};

private:
    Bool_t InitEvent();

    void ProcessGammas(Int_t iCut,EEventPlaneMethod iEP);
    void ProcessPi0s(Int_t iCut,EEventPlaneMethod iEP);
    void ProcessQA();

    Double_t GetPhiwrtRP(Double_t dPhi);
    Double_t GetPhotonPhiwrtRP(AliAODConversionPhoton *gamma,EEventPlaneMethod iEP,Bool_t bDoFlattening=kTRUE);
    Double_t GetPi0PhiwrtRP(AliAODConversionMother *pi0,EEventPlaneMethod iEP,Bool_t bDoFlattening=kTRUE);
    Double_t GetChargedPhiwrtRP(AliVTrack *charged,EEventPlaneMethod iEP,Bool_t bDoFlattening=kTRUE);
    void GetPhotondNdPhi(Int_t *dNdPhi,Int_t iEP,Int_t iCut=0);
    void GetChargeddNdPhi(Int_t *dNdPhi,Int_t &ntot,Int_t iEP);
    Int_t GetPhiBin(Double_t phiwrt);
    Int_t GetPhotonPhiBin(AliAODConversionPhoton *gamma,Int_t iEP);
    Double_t GetMCPhotonPhiwrtRP(TParticle *gamma,EEventPlaneMethod iEP,Bool_t bDoFlattening=kTRUE);
    TVector2 GetEPContribution(AliAODConversionPhoton *gamma);
    Double_t GetEventPlaneAngle(EEventPlaneMethod EPmethod,Double_t eta=0,AliAODConversionPhoton *gamma0=NULL,AliAODConversionPhoton *gamma1=NULL,Bool_t bDoFlattening=kTRUE);
    Double_t GetTPCSubEPEta(EEventPlane ep);
    Double_t GetCorrectedTPCEPAngle(AliAODConversionPhoton *gamma0=NULL,AliAODConversionPhoton *gamma1=NULL,Bool_t bDoFlattening=kTRUE);
    Bool_t SetCentrality();
    void ProcessEventPlane();
    Int_t GetRadialBin(Double_t radius);
    Int_t GetRunIndex(Int_t run);
    Double_t ApplyFlattening(Double_t phi,EEventPlane ep);
    Bool_t GetTPCEventPlane();

    void GetV0EP(AliVEvent * event,Double_t &rpv0a,Double_t &rpv0c);
    void LoadVZEROCalibration(Int_t run);
    void LoadTPCCalibration(Int_t run);
    
    Double_t GetWeight(TObject* track1);
    Double_t GetPhiWeight(TObject* track1);
    TH1F* SelectPhiDist(AliVTrack *track);

    Double_t GetPsiInRange(Double_t phi);

    TObjArray* GetEventPlaneTracks(Int_t &maxID);
    TVector2 GetContributionEP(AliVTrack *track);
    Int_t GetAODEPTrackFilterBit();
   
    // Constants

    enum Ebinsgamma{
	kGammaPt=0,
	kGammadPhi,
	kGammaCent,
	kGammaEPM,
	knbinsGamma
    };

    enum Ebinspi0{
	kPi0Pt=0,
	kPi0Mass,
	kPi0dPhi,
	kPi0Cent,
	kPi0EPM,
	knbinsPi0
    };

    static const Int_t knbinsGammaMult=3;
  
    static const Int_t kGCnYBinsSpectra = 80;
    static const Double_t kGCfirstYBinSpectra = 0.;
    static const Double_t kGClastYBinSpectra = 8.;

    // Class variables and pointer

    AliV0ReaderV1 *fV0Reader; // V0Reader
    Int_t fNCuts; // Number of Photon Cuts for v2 analysis
    AliConversionSelection **fConversionSelection; //[fNCuts] Selection of Particles for given Cut
    TClonesArray *fConversionGammas; //Reconstructed Photons;
    Int_t fNCentralityBins; // Number of Centrality Bins
    Double_t *fCentralityBins; //[fNCentralityBins] CentralityBins for Analysis
    Float_t fCentrality; //Event Centrality
    Int_t fCentralityBin; // Event Centrality Bin
    Int_t fNBinsPhi; // Number of Phi wrt RP bins
    AliEventplane *fEP; // Event Plane Pointer
    Bool_t fUseTPCOnlyTracks; // Use TPC Only Tracks for EP
    Double_t fEtaMax; // Eta Max for analysis;
    Double_t fEtaGap; // Eta Gap
    Double_t fRPTPCEtaA; // TPCEtaA event plane
    Double_t fRPTPCEtaC; // TPCEtaC event plane
    Double_t fRPV0A; // V0A event plane
    Double_t fRPV0C; // V0C event plane
    Double_t fRPTPC; // TPC event plane
    Double_t fRPTPCEtaABF; // TPCEtaA event plane before flattening
    Double_t fRPTPCEtaCBF; // TPCEtaC event plane before flattening
    Double_t fRPV0ABF;// V0A event plane before flattening
    Double_t fRPV0CBF;// V0C event plane before flattening
    Double_t fRPTPCBF;// TPC event plane before flattening
    AliConversionCuts *fConversionCuts; // Cuts used by the V0Reader
    TRandom3 *fRandomizer; // Randomizer for Event Plane Randomisation
    TList *fOutputList; // List for Output (Histograms etc.)
    EPDGCode fMesonPDGCode; // PDG Code of the processed Meson (for MC truth)
    Double_t fInvMassRange[2]; // Inv Mass Range
    Double_t fDeltaPsiRP; // Difference between subEventPlane angles
    Int_t fRunNumber; // current run number
    Int_t fRunIndex; // current internal run index
    Int_t fNEPMethods; // number of EP methods
    Bool_t fFillQA; // Fill QA Histograms
    Int_t fHarmonic; // Harmonic to be analyzed (v2,v3,..)
    Double_t fPsiMax; //Range for Psi
    TString  fPeriod; //"LHC11h","LHC10h"
    Bool_t fIsAOD; // Is AOD? else ESD
    TH1F*	 fPhiDist[4];			// array of Phi distributions used to calculate phi weights
    THnSparse *fSparseDist;               //! THn for eta-charge phi-weighting
    TH1F *fHruns;                         // information about runwise statistics of phi-weights
    Bool_t fDoEPFlattening; // Do flattening
    Int_t fEPSelectionMask[knEPMethod]; // Which EP methods shall be applied
    // Histograms

    TH1F *hNEvents;

    // RP
    TH2F *hRPTPC;
    TH2F *hRPV0A;
    TH2F *hRPV0C;
    TH2F *hRPTPCAC;
    TH2F *hRPV0ATPC;
    TH2F *hRPV0CTPC;
    TH2F *hRPV0AC;
    TH2F *hCos2TPC;
    TH2F *hCos2V0ATPC;
    TH2F *hCos2V0CTPC;
    TH2F *hCos2V0AC;
    TH2F *hRPTPCEtaA;
    TH2F *hRPTPCEtaC;
    TH2F *hRPTPCEtaAC;
    TH2F *hCos2TPCEta;
    TH2F *hCos2V0ATPCEtaA;
    TH2F *hCos2V0ATPCEtaC;
    TH2F *hCos2V0CTPCEtaA;
    TH2F *hCos2V0CTPCEtaC;
    TH2F *hCos2SumWeights;
    TH2F *hEtaTPCEP;

    // Gamma
    TH2F *hGammaMultCent;
    TH2F **hGammaPhi;
    TH2F *hMultChargedvsNGamma;
    TH2F *hMultChargedvsVZERO;
    TH2F *hMultChargedvsSPD;

    THnSparseF *hGammadNdPhi;
    THnSparseF *hGammaMultdPhiTRUE;
    THnSparseF *hGammaMultdPhiRECOTRUE;
    THnSparseF *hGammaMultTRUE;
    THnSparseF *hGammaMultRECOTRUE;
    THnSparseF **hGammaMultdPhi;
    THnSparseF **hGammaMult;

    THnSparseF **hGamma;
    THnSparseF *hGammaFull;

    THnSparseF *hCharged;

    // Pi0
    THnSparseF **hPi0;
    THnSparseF **hPi0BG;

    //V0 Calibration

    static const Int_t nCentrBinV0 = 9; // # cenrality bins
    TProfile *fMultV0;                  // object containing VZERO calibration information
    Float_t fV0Cpol,fV0Apol;            // loaded by OADB
    Float_t fMeanQ[nCentrBinV0][2][2];    // and recentering
    Float_t fWidthQ[nCentrBinV0][2][2];   // ...

    //Event Plane
    //THnSparse *hEPVertex;
    THnSparse *hEPQA;

    ClassDef(AliAnalysisTaskPi0v2, 4); // example of analysis
};

#endif

