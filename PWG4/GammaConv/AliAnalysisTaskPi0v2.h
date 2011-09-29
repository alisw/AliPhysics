#ifndef AliAnalysisTaskPi0v2_cxx
#define AliAnalysisTaskPi0v2_cxx

#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"
#include "AliESDtrack.h"

#include "AliAnalysisTaskSE.h"
#include "TRandom3.h"
#include "TH1.h"
#include "AliLog.h"
#include "AliV0ReaderV1.h"
#include "AliAnalysisTaskPi0Reconstruction.h"
using namespace std;

class AliAnalysisTaskPi0v2 : public AliAnalysisTaskPi0Reconstruction {

public:
    AliAnalysisTaskPi0v2(const char *name);
    virtual ~AliAnalysisTaskPi0v2();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetNBinsPhi(Int_t nphibins){fNBinsPhi=nphibins;}

private:

    void ProcessGammas();
    void ProcessPi0s();
    void ProcessBGPi0s();

    void ProcessChargedParticles();
    void ProcessEventPlaneResolution();

    Int_t GetPhiBin(Double_t phi);

    Int_t fNBinsPhi;

    // Histograms

    TH2F ***hInvMassmphi;
    TH2F ***hBGmphi;
    TH2F *hInvMassPt;
    TH1F *hInvMass;
    TH2F *hBGPt;
    TH1F *hBG;

    TH1F **hPi0TRUE;
    TH2F **hPi0RECOTRUE;

    // RP
    TH1F *hRP;
    TH2F *hRPCentrality;
    TH2F *hRPSubevents;
    TH1F *hRPTrackCuts;
    TH3F *hRPPhiTracks;
    TH2F *hRPDeltaRP;
    TH2F *hRPQxQy;
    TH2F *hRPCosDeltaRP;

    // Others
    TH2F **hInclv2PtInvMassCentrality;

    TH1F **hChargedPt;
     
    // Gamma
    TH2F **hGammav2;



    AliAnalysisTaskPi0v2(const AliAnalysisTaskPi0v2&); // not implemented
    AliAnalysisTaskPi0v2& operator=(const AliAnalysisTaskPi0v2&); // not implemented
  
    ClassDef(AliAnalysisTaskPi0v2, 2); // example of analysis
};

#endif

