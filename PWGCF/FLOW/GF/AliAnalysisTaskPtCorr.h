#ifndef PTCORRTASK_H
#define PTCORRTASK_H

#include "AliAnalysisTaskSE.h"
#include "complex.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliGFWCuts.h"
#include "TString.h"
#include "TRandom.h"
#include "AliESDtrackCuts.h"
#include "AliProfileBS.h"
#include "AliCkContainer.h"
#include "AliPtContainer.h"
#include "TH3D.h"

class TList;
class TH1;
class TRandom;
class TAxis;
class AliVEvent;

using namespace std;

namespace PtCorrFlags {
    enum {
      noeff = 1,
      consteff = 2,
      gausseff = 4,
      flateff = 8,
      powereff = 16,
      realeffin = 32
    };
}

class AliAnalysisTaskPtCorr : public AliAnalysisTaskSE
{
  public:
    AliAnalysisTaskPtCorr();
    AliAnalysisTaskPtCorr(const char *name, bool IsMC, bool isOnTheFly, unsigned int fl_eff, TString ContSubfix); //pseudoeff: 0 = no eff, 1 = subset of particles, 2 = gaussian distribution
    virtual ~AliAnalysisTaskPtCorr();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void NotifyRun();
    void SetSystFlag(Int_t newval) { if(!fGFWSelection) fGFWSelection = new AliGFWCuts(); fGFWSelection->SetupCuts(newval); };
    void SetCentralityEstimator(TString newval) { if(fCentEst) delete fCentEst; fCentEst = new TString(newval); };
    void SetV0MBins(int nBins, double* bins);
    void SetMultiplicityBins(int nBins, double* bins);
    void SetPtBins(int nBins, double *ptbins);
    void SetContSubfix(TString newval) {if(fContSubfix) delete fContSubfix; fContSubfix = new TString(newval); if(!fContSubfix->IsNull()) fContSubfix->Prepend("_"); };
    void SetMPar(int m) { mpar = m; }
    void OverrideMC(bool ismc) { fIsMC = ismc; }
    void SetTrigger(unsigned int newval) {fTriggerType = newval; };
    void SetEta(double eta) { fEta = eta; }
    void SetEtaGap(double eta) { fEtaGap = eta; }
    void TurnOffPileup(bool off) { fPileupOff = off; }
    void SetPileupCut(double cut) { fPUcut = cut; }
    void SetEventWeight(unsigned int weight) { fEventWeight = weight; }
    void SetUseRecNchForMc(bool userec) { fUseRecNchForMC = userec; }
    void SetUseNch(bool usench) { fUseNch = usench; }
    void SetEtaNch(double etanch) { fEtaNch = etanch; }
    void SetNBootstrap(double nboot) { fNbootstrap = nboot; }
    void SetEffFlags(unsigned int fl) { eff_flags = fl; }
    void UseTPCOnlyTracks(bool usetpc) { fUseTPConly = usetpc; };
    void SetPseudoEffPars(double consteff = 0.8, double sigmaeff = 0.05) { fConstEff = consteff; fSigmaEff = sigmaeff; }
  protected:
    AliEventCuts            fEventCuts;
  private:
    AliAnalysisTaskPtCorr(const AliAnalysisTaskPtCorr&);
    AliAnalysisTaskPtCorr& operator=(const AliAnalysisTaskPtCorr&);
    TString*                fCentEst;
    int                     fRunNo; //!
    int                     fSystFlag;
    TString*                fContSubfix;
    bool                    fIsMC; 
    AliMCEvent*             fMCEvent; //!
    TList*                  fCorrList; //!
    TList*                  fQAList; //!
    TList*                  fEfficiencyList;
    TH1D**                  fEfficiencies;
    TH2D**                  fPowerEfficiencies;     
    TString                 fWeightSubfix;
    AliGFWCuts*             fGFWSelection; 
    AliGFWCuts*             fGFWnTrackSelection;
    TAxis*                  fV0MAxis; 
    TAxis*                  fMultiAxis; 
    double*                 fMultiBins; //!
    int                     fNMultiBins; //!
    TAxis*                  fPtAxis;  
    double*                 fPtBins; //!
    int                     fNPtBins; //!
    double                  fEta;
    double                  fEtaNch;
    double                  fEtaGap;
    double                  fPUcut;
    TRandom*                fRndm;
    int                     fNbootstrap; 
    bool                    fUseWeightsOne;
    bool                    fUseRecNchForMC;
    bool                    fPileupOff;
    bool                    fUseNch;
    bool                    fUseTPConly;
    unsigned int            eff_flags;
    int                     mpar;
    double                  fConstEff;
    double                  fSigmaEff;
    vector<vector<double>>  wp;
    vector<vector<double>>  wpP;
    vector<vector<double>>  wpN;
    unsigned int            fEventWeight;
    TH1D*                   fV0MMulti;    //!
    AliPtContainer*         fpt;      //!
    TH1D*                   fEventCount;
    TH2D*                   fNchTrueVsRec; //!
    TH2D*                   fV0MvsMult; //!
    TH3D*                   fPtMoms;  //!
    TH1D*                   fPtDistB;  //!
    TH1D*                   fPtDistA;  //!
    TH1D*                   fPtDistC;  //!
    TH3D*                   fPtDCA; //!
    TH2D*                   fPtVsNTrk; //!
    unsigned int            fTriggerType;
    bool                    fOnTheFly;
    double                  fImpactParameter; 
    map<double,double>      centralitymap;
    int                     EventNo;
    bool AcceptAODTrack(AliAODTrack *tr, double *ltrackXYZ, const double &ptMin, const double &ptMax, double *vtxp);
    bool AcceptAODTrack(AliAODTrack *mtr, double *ltrackXYZ, const double &ptMin, const double &ptMax, double *vtxp, int &nTot);
    bool AcceptAODEvent(AliAODEvent *ev, double *vtxXYZ);
    AliMCEvent* getMCEvent();
    bool CheckTrigger(Double_t lCent);
    double getCentrality();
    void FillPtCorr(AliAODEvent* fAOD, const double &l_Cent, double *vtxXYZ);
    void FillWPCounter(vector<vector<double>> &inarr, double w, double p);
    void FillWPCounter(vector<vector<double>> &inarr, vector<double> w, double p);
    int GetNTracks(AliAODEvent* fAOD, const Double_t &ptmin, const Double_t &ptmax, Double_t *vtxp);
    double *GetBinsFromAxis(TAxis *inax);
    
  ClassDef(AliAnalysisTaskPtCorr,1);
};

#endif