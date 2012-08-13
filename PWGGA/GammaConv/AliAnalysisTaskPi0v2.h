#ifndef AliAnalysisTaskPi0v2_cxx
#define AliAnalysisTaskPi0v2_cxx

#include "AliAnalysisTaskSE.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "AliLog.h"
#include "AliConversionSelection.h"
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

    enum EPDGCode{
	kPi0=111,
	kEta=221
    };

    static const Int_t knBinsPhi=6;

    AliAnalysisTaskPi0v2(const char *name);
    virtual ~AliAnalysisTaskPi0v2();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetCentralityBins(Double_t *bins,Int_t nbins);
    void SetRadialBins(Float_t *bins,Int_t nbins);

    void SetMeson(EPDGCode meson){fMesonPDGCode=meson;}

    void SetNBinsPhi(Int_t nbins){fNBinsPhi=nbins;}
    void SetV0Reader(AliV0ReaderV1 *v0Reader){fV0Reader=v0Reader;}
    void SetInvMassRange(Double_t range[2]){fInvMassRange[0]=range[0];fInvMassRange[1]=range[1];};
    void SetEtaGap(Double_t gapsize){fEtaGap=gapsize;};

    void SetMesonCuts(const TString cut);
    void SetCuts(TString *cutarray,Int_t ncuts);

    void SetFillQA(Bool_t fill){fFillQA=fill;}

    void SetWeightMultiplicity(Bool_t b){fWeightMultiplicity=b;}

private:
    Bool_t InitEvent();

    void ProcessGammas(Int_t iCut,EEventPlaneMethod iEP);
    void ProcessPi0s(Int_t iCut,EEventPlaneMethod iEP);
    void ProcessQA();

    void InitConversionSelection();
    Double_t GetPhiwrtRP(Double_t phi);
    Double_t GetPhotonPhiwrtRP(AliAODConversionPhoton *gamma,EEventPlaneMethod iEP);
    Double_t GetPi0PhiwrtRP(AliAODConversionMother *pi0,EEventPlaneMethod iEP);
    Double_t GetChargedPhiwrtRP(AliVTrack *charged,EEventPlaneMethod iEP);
    void GetPhotondNdPhi(Int_t *dNdPhi,Int_t iEP,Int_t iCut=0);
    void GetChargeddNdPhi(Int_t *dNdPhi,Int_t &ntot,Int_t iEP);
    Int_t GetPhiBin(Double_t phiwrt);
    Int_t GetPhotonPhiBin(AliAODConversionPhoton *gamma,Int_t iEP);
    Double_t GetMCPhotonPhiwrtRP(TParticle *gamma,EEventPlaneMethod iEP);
    TVector2 GetEPContribution(AliAODConversionPhoton *gamma);
    Double_t GetEventPlaneAngle(EEventPlaneMethod EPmethod,Double_t eta=0,AliAODConversionPhoton *gamma0=NULL,AliAODConversionPhoton *gamma1=NULL);
    Double_t GetTPCSubEPEta(Double_t etamin,Double_t etamax);
    Double_t GetCorrectedTPCEPAngle(AliAODConversionPhoton *gamma0=NULL,AliAODConversionPhoton *gamma1=NULL);
    Bool_t SetCentrality();
    void ProcessEventPlane();
    Int_t GetRadialBin(Double_t radius);
    Int_t GetRunIndex(Int_t run);

    // For V0 EP
    void GetV0EP(AliVEvent * event);
    void OpenInfoCalibration(Int_t run);

    Double_t ApplyFlatteningTPC(Double_t phi, Double_t c);
    Double_t ApplyFlatteningV0A(Double_t phi, Double_t c);
    Double_t ApplyFlatteningV0C(Double_t phi, Double_t c);


    // Constants

    static const Int_t knbinsGamma=5;
    static const Int_t knbinsGammaMult=4;
    static const Int_t knbinsPi0=5;

    static const Int_t kGCnYBinsSpectra = 80;
    static const Double_t kGCfirstYBinSpectra = 0.;
    static const Double_t kGClastYBinSpectra = 8.;

    // Class variables and pointer

    AliV0ReaderV1 *fV0Reader; // V0Reader
    AliConversionSelection **fConversionSelection; // Selection of Particles for given Cut
    TClonesArray *fConversionGammas; //Reconstructed Photons;
    Int_t fNCentralityBins; // Number of Centrality Bins
    Double_t *fCentralityBins; // CentralityBins for Analysis
    Float_t fCentrality; //Event Centrality
    Int_t fCentralityBin; // Event Centrality Bin
    Int_t fNRadialBins; // Number of Radial Bins for Photon Conversion Point
    Double_t *fRadialBins; // Radial Bins for Photons Conversion Point
    Int_t fNBinsPhi; // Number of Phi wrt RP bins
    AliEventplane *fEP; // Event Plane Pointer
    Bool_t fWeightMultiplicity; // Use Multiplicity Weight
    Double_t fEtaMax; // Eta Max for analysis;
    Double_t fEtaGap; // Eta Gap
    Double_t fRPTPCEtaA; // TPCEtaA event plane
    Double_t fRPTPCEtaC; // TPCEtaC event plane
    Double_t fRPV0A; // V0A event plane
    Double_t fRPV0C; // V0C event plane
    Int_t fNCuts; // NUmber of Photon Cuts for v2 analysis
    TList *fCutList; // Cuts for Photon v2 analysis
    AliConversionCuts *fConversionCuts; // Cuts used by the V0Reader
    TRandom3 *fRandomizer; // Randomizer for Event Plane Randomisation
    TList *fOutputList; // List for Output (Histograms etc.)
    EPDGCode fMesonPDGCode; // PDG Code of the processed Meson (for MC truth)
    Double_t *fInvMassRange; // Inv Mass Range
    Double_t fDeltaPsiRP; // Difference between subEventPlane angles
    Int_t fRunNumber; // current run number
    Int_t fRunIndex; // current internal run index
    Int_t fNEPMethods; // number of EP methods
    Bool_t fFillQA; // Fill QA Histograms

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


    // Gamma
    TH2F *hGammaMultCent;
    TH3F **hGammaPhi;
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


    AliAnalysisTaskPi0v2(const AliAnalysisTaskPi0v2&); // not implemented
    AliAnalysisTaskPi0v2& operator=(const AliAnalysisTaskPi0v2&); // not implemented
  
    ClassDef(AliAnalysisTaskPi0v2, 3); // example of analysis
};

#endif

