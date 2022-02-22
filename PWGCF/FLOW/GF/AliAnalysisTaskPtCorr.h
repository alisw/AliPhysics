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
    AliAnalysisTaskPtCorr(const char *name, bool IsMC, TString ContSubfix);
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
    void SetEta(double eta) { fEta = eta; }
    void SetEtaGap(double eta) { fEtaGap = eta; }
    void SetPtRange(double ptmin, double ptmax) { fPtMin = ptmin; fPtMax = ptmax; }
    void SetEventWeight(unsigned int weight) { fEventWeight = weight; }
    void SetUseWeightsOne(bool use) { fUseWeightsOne = use; }
  protected:
    AliEventCuts            fEventCuts;
  private:
    AliAnalysisTaskPtCorr(const AliAnalysisTaskPtCorr&);
    AliAnalysisTaskPtCorr& operator=(const AliAnalysisTaskPtCorr&);
    static double           fFactorial[9];
    static int              fSign[9];
    TString*                fCentEst;
    int                     fRunNo; //!
    int                     fSystFlag;
    TString*                fContSubfix;
    bool                    fIsMC; 
    AliMCEvent*             fMCEvent; //!
    TList*                  fCorrList; //!
    TList*                  fEfficiencyList;
    TH1D**                  fEfficiencies;     
    TString                 fWeightSubfix;
    AliGFWCuts*             fGFWSelection; 
    TAxis*                  fV0MAxis; 
    TAxis*                  fMultiAxis; 
    double*                 fMultiBins; //!
    int                     fNMultiBins; //!
    TAxis*                  fPtAxis;  
    double*                 fPtBins; //!
    int                     fNPtBins; //!
    double                  fEta;
    double                  fEtaGap;
    double                  fPtMin;
    double                  fPtMax;
    TRandom*                fRndm;
    int                     fNbootstrap; 
    bool                    fUseWeightsOne;
    int                     mpar;
    unsigned int            fEventWeight;
    TH1D*                   fV0MMulti;    //!
    AliProfileBS**          fptcorr;     //!
    AliProfileBS**          fptcorrP;     //!
    AliProfileBS**          fptcorrN;     //!
    AliProfileBS*           pfmpt;      //!
    AliPtContainer*         fck;      //!
    AliPtContainer*         fskew;      //!
    AliPtContainer*         fkur;      //!
    AliPtContainer*         f5;      //!
    AliPtContainer*         f6;      //!
    unsigned int            fTriggerType;
    bool                    fOnTheFly;
    double                  fImpactParameter; 
    map<double,double>      centralitymap;
    

    bool AcceptAODTrack(AliAODTrack *tr, double *ltrackXYZ, const double &ptMin, const double &ptMax, double *vtxp);
    bool AcceptAODEvent(AliAODEvent *ev, double *vtxXYZ);
    bool AcceptMCEvent(AliVEvent* inev);
    bool CheckTrigger(Double_t lCent);
    double getCentrality();
    void FillPtCorr(AliVEvent* ev, const double &VtxZ, const double &l_Cent, double *vtxXYZ);
    void FillWPCounter(double inarr[10][10], double w, double p);
    void FillMomentumCorrelation(double wp[10][10], double wpP[10][10], double wpN[10][10], const double &l_cent, double &rn);
    void CalculateSingleEvent(vector<double> &corr, vector<double> &sumw, double wp[10][10]);
    void CalculateSubEvent(vector<double> &corrP, vector<double> &sumwP, vector<double> &corrN, vector<double> &sumwN, double wpP[10][10], double wpN[10][10]);
    void FillCorrelationProfile(const vector<double> &corr, const vector<double> &sumw, const double &l_cent, double &rn);
    void FillSubEventProfile(const vector<double> &corrP,const vector<double> &sumwP,const vector<double> &corrN, const vector<double> &sumwN, const double &l_cent, double &rn);
    double OrderedAddition(std::vector<double> vec, int size);
    double *GetBinsFromAxis(TAxis *inax);
    
  ClassDef(AliAnalysisTaskPtCorr,1);
};

#endif