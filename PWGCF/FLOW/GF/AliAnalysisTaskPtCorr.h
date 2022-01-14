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
#include "AliGFWWeights.h"
#include "AliProfileBS.h"
//#include "AliPtContainer.h"

class TList;
class TH1;
class TRandom;
class TAxis;
class AliVEvent;

using namespace std;

class AliAnalysisTaskPtCorr : public AliAnalysisTaskSE
{
  public:
    AliAnalysisTaskPtCorr();
    AliAnalysisTaskPtCorr(const char *name, bool IsMC, bool dodynamics, TString analysisStage, TString ContSubfix);
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
    void OnTheFly(bool otf) { fOnTheFly = otf; }
    void SetTrigger(unsigned int newval) {fTriggerType = newval; };
    void SetUseWeightsOne(bool use) { fUseWeightsOne = use; }
    void DoDynamicCorrelations(bool dodyn) { fdoDynamicCorr = dodyn; }
    AliEventCuts            fEventCuts;
    
  private:
    static int              fFactorial[9];
    static int              fSign[9];
    TString*                fCentEst;
    int                     fRunNo; //!
    int                     fSystFlag;
    TString*                fContSubfix;
    bool                    fIsMC; 
    TList*                  fCorrList; //!
    TList*                  fDynList; //!
    TList*                  fWeightList; //!
    TList*                  fEfficiencyList; //!
    TList*                  fmptList; //!
    TH2D**                  fEfficiency;     
    TH1D**                  fEfficiencies;     
    AliGFWWeights**         fWeights; //!
    TH1D*                   fmpt;
    TString                 fWeightSubfix;
    AliGFWCuts*             fGFWSelection; 
    TAxis*                  fV0MAxis; //!
    TAxis*                  fMultiAxis; //!
    double*                 fMultiBins; //!
    double                  fNMultiBins; //!
    TAxis*                  fPtAxis;  //!
    double*                 fPtBins; //!
    int                     fNPtBins; //!
    double                  fEta;
    double                  fPtMin;
    double                  fPtMax;
    int                     fAnalysisStage;
    TRandom*                fRndm;
    int                     fNbootstrap; 
    bool                    fUseWeightsOne;
    bool                    fdoDynamicCorr;
    int                     mpar;
    TH1D*                   fV0MMulti;    //!
    //AliPtContainer**        fdyncorrnompt;     //!
    AliProfileBS**          fptcorr;     //!
    AliProfileBS**          fdyncorr;     //!
    AliProfileBS*           pfmpt;      //!

    vector<double>          corr;
    vector<double>          dyncorr;
    vector<double>          sumw;

    unsigned int            fTriggerType;
    bool                    fOnTheFly;
    double                  fImpactParameter; 
    map<double,double>      centralitymap;
    

    bool AcceptAODTrack(AliAODTrack *tr, double *ltrackXYZ, const double &ptMin, const double &ptMax, double *vtxp);
    bool AcceptAODEvent(AliAODEvent *ev, double *vtxXYZ);
    bool AcceptMCEvent(AliVEvent* inev);
    bool CheckTrigger(Double_t lCent);
    double getCentrality();
    void FillWeights(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
    void FillPtCorr(AliVEvent* ev, const double &VtxZ, const double &l_Cent, double *vtxXYZ);
    void FillMPT(AliVEvent* ev, const double &VtxZ, const double &l_cent, double *vtxXYZ);
    void FillCorrelationProfiles(const double &l_cent, double &rn);
    void FillDynamicProfiles(const double &l_cent, double &rn);
    template <typename T> void FillWPCounter(T& inarr, double w, double p);
    void MptCounter(double* inarr, int icent);
    template<typename T> void getMomentumCorrelation(T& wp);
    template<typename T> void getDynamicMomentumCorrelation(T& wp, double* mpt);
    //template<typename T> void getDynCorrSkipMpt(T& wp, const double &l_cent, double &rn);
    //template<typename T> void getTermsOfMpt(int k, int m, T& term, T& wp);
    double OrderedAddition(std::vector<double> vec, int size);
    int binomial(const int n, const int m);
    template<typename T> double DynamicP(int k, T& wp, double* mpt);
    bool LoadWeights(const int &lRunNo);
    int GetAnalysisStage(TString instr);
    double *GetBinsFromAxis(TAxis *inax);
    
  ClassDef(AliAnalysisTaskPtCorr,1);
};

#endif