#ifndef AliAnalysisTaskPi0EtaV2_H
#define AliAnalysisTaskPi0EtaV2_H

#include "AliAnalysisTaskSE.h"
#include "TH3.h"
#include "TH3F.h"
#include "TProfile2D.h"

class AliAnalysisTaskPi0EtaV2 : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskPi0EtaV2();
    AliAnalysisTaskPi0EtaV2(const char *name);
    virtual ~AliAnalysisTaskPi0EtaV2();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    void IfVZEROCalibOn(bool bVZEROCalibOn) { this->IsVZEROCalibOn = bVZEROCalibOn; }
    void IfQAVZERO(bool bQAVZERO) { this->IsQAVZERO = bQAVZERO; }
    bool LoadCalibHistForThisRun();
    bool GetVZEROPlane();
    double GetEventPlane(double qx, double qy, double harmonic);
    void SetPeriod(TString period) { this->fPeriod = period; }
    void SetnCuts(Int_t nCuts)
    {
        fnCuts = nCuts;
    }
    void SetCent(Double_t cent)
    {
        centSPD1 = cent;
    }
    void SetCutList(TList *CutArray)
    {
        fEventCutArray = CutArray;
    }
    void SetTrainconfig(Int_t config)
    {
        fTrainConfig = config;
    }

private:
    AliVEvent *fInputEvent; // current event
    Int_t fnCuts;           // total cuts
    Int_t fiCut;            // current cut
    Int_t fTrainConfig;     // fTrainConfig
    TString fPeriod;
    TString fOutputAODBranchName;                     ///<  Name of output clusters AOD branch
    TString fOutputBGBranchName;                      ///<  Name of output of background
    TList *fEventCutArray;                            // List with Event Cuts
    TList *fOutputContainer;                          // Output container
    TList **fCutFolder;                               // Array of lists for containers belonging to cut
    TList **fESDList;                                 // Array of lists with histograms with reconstructed properties
    TList **fQAList;                                  // Array of lists with histograms with V0 plane QA
    TH3F **fHistoMotherInvMassPtPhiV0A;               //! array of histogram with signal + BG for same event photon pairs in deta phi, inv Mass, pt
    TH3F **fHistoMotherInvMassPtPhiV0C;               //! array of histogram with signal + BG for same event photon pairs in deta phi, inv Mass, pt
    TH3F **fHistoMotherBackInvMassPtdPhiV0A;          //! array of histogram with BG for mixed event photon pairs in deta phi, inv Mass, pt
    TH3F **fHistoMotherBackInvMassPtdPhiV0C;          //! array of histogram with BG for mixed event photon pairs in deta phi, inv Mass, pt
    TH1D **fEventCount;                               //! array of histogram of event count in centBins
    TH2F **fHistoMotherInvMassPtV0CInPlane;           // ray of histogram with signal + BG for same event photon pairs In Plane, inv Mass, pt
    TH2F **fHistoMotherInvMassPtV0AInPlane;           // ray of histogram with signal + BG for same event photon pairs In Plane, inv Mass, pt
    TH2F **fHistoMotherInvMassPtV0COutPlane;          // ray of histogram with signal + BG for same event photon pairs Out Plane, inv Mass, pt
    TH2F **fHistoMotherInvMassPtV0AOutPlane;          // ray of histogram with signal + BG for same event photon pairs Out Plane, inv Mass, pt
    TH2F **fHistoMotherBackInvMassPtV0CInPlane;       // ray of histogram with only BG for same event photon pairs In Plane, inv Mass, pt
    TH2F **fHistoMotherBackInvMassPtV0AInPlane;       // ray of histogram with only BG for same event photon pairs In Plane, inv Mass, pt
    TH2F **fHistoMotherBackInvMassPtV0COutPlane;      // ray of histogram with only BG for same event photon pairs Out Plane, inv Mass, pt
    TH2F **fHistoMotherBackInvMassPtV0AOutPlane;      // ray of histogram with only BG for same event photon pairs Out Plane, inv Mass, pt
    TProfile2D **fHistoMotherInvMassPtV0CCos2phi;     // ray of histogram with signal + BG for same event photon pairs cos(2(phi-Psi)) inv Mass, pt
    TProfile2D **fHistoMotherInvMassPtV0ACos2phi;     // ray of histogram with signal + BG for same event photon pairs cos(2(phi-Psi)) inv Mass, pt
    TProfile2D **fHistoMotherBackInvMassPtV0CCos2phi; // ray of histogram with only BG for same event photon pairs cos(2(phi-Psi)), inv Mass, pt
    TProfile2D **fHistoMotherBackInvMassPtV0ACos2phi; // ray of histogram with only BG for same event photon pairs cos(2(phi-Psi)), inv Mass, pt
    std::map<int, int> *runNumList;
    TH1F **fhistPhi; // ray of  hist of candidate in a ray invmass range
    TH1F **fhistPhiBG;
    TH1I *fHistRunNumBin;

    int runNum;
    int oldRunNum;
    double centSPD1;
    bool IsVZEROCalibOn; // switch for VZERO qn calib
    bool IsQAVZERO;      // switch for Fill VZERO calib qa
    bool IsUseInOutPlane;
    TList *fListVZEROCalib; // read list for V0 Calib
    TFile *fVZEROCalibFile;
    double fPsi2V0C;
    double fPsi2V0A;
    TH2D **fHist2DPsi2V0CCent;
    TH2D **fHist2DPsi2V0ACent;
    TH1D *hMultV0; // Dobrin
    AliOADBContainer *contMult;
    AliOADBContainer *contQxncm;
    AliOADBContainer *contQyncm;
    AliOADBContainer *contQxnam;
    AliOADBContainer *contQynam;
    // 18q/r
    TH2F *fHCorrectV0ChWeghts;
    // V0C
    TProfile **fProfileV0CQxCentGE;
    TProfile **fProfileV0CQyCentGE;
    TProfile **fProfileV0CQxVtxGE;
    TProfile **fProfileV0CQyVtxGE;
    TH2D **fHist2CalibPsi2V0CCentGE;
    TProfile **fProfileV0AQxCentGE;
    TProfile **fProfileV0AQyCentGE;
    TProfile **fProfileV0AQxVtxGE;
    TProfile **fProfileV0AQyVtxGE;
    TH2D **fHist2CalibPsi2V0ACentGE;
    TProfile **fProfileV0CQxCentRC;
    TProfile **fProfileV0CQyCentRC;
    TProfile **fProfileV0CQxVtxRC;
    TProfile **fProfileV0CQyVtxRC;
    TH2D **fHist2CalibPsi2V0CCentRC;
    TProfile **fProfileV0AQxCentRC;
    TProfile **fProfileV0AQyCentRC;
    TProfile **fProfileV0AQxVtxRC;
    TProfile **fProfileV0AQyVtxRC;
    TH2D **fHist2CalibPsi2V0ACentRC;
    TProfile **fHist2V0Res;
    TH1D *hQx2mV0[2];
    TH1D *hQy2mV0[2];

    AliAnalysisTaskPi0EtaV2(const AliAnalysisTaskPi0EtaV2 &);            // not implemented
    AliAnalysisTaskPi0EtaV2 &operator=(const AliAnalysisTaskPi0EtaV2 &); // not implemented

    ClassDef(AliAnalysisTaskPi0EtaV2, 2);
};

#endif
