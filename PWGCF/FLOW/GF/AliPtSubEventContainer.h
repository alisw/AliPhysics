#ifndef ALIPTSUBEVENTCONTAINER__H
#define ALIPTSUBEVENTCONTAINER__H
#include <iostream>
#include "TProfile.h"
#include "TNamed.h"
#include "TList.h"
#include "TCollection.h"

using namespace std;

namespace WeightSpace {
    enum kWeight {
        kOne,
        kWpt,
        kWperms,
        kWmaxperm
    };
    enum kObs {
        kCorr,
        kCum,
        kNorm
    };
};

class AliPtSubEventContainer: public TNamed {
    public:
        AliPtSubEventContainer();
        ~AliPtSubEventContainer();
        AliPtSubEventContainer(const char* name, const char* title, int m);
        void Initialize(int nbinsx, double* xbins);
        void Initialize(int nbinsx, double xlow, double xhigh);
        
        void InitializeTwoSub(int nbinsx, double xlow, double xhigh);
        void InitializeThreeSub(int nbinsx, double xlow, double xhigh);
        void FillTwoSubAnalsysis(const vector<vector<double>> &inarrSubA, const vector<vector<double>> &inarrSubB, const double &lMult);
        void FillThreeSubAnalsysis(const vector<vector<double>> &inarrSubA, const vector<vector<double>> &inarrSubB, const vector<vector<double>> &inarrSubC, const double &lMult);

        void SetNamedPostfix(TString postfix){ fPostfix = postfix; }

        void FillRecursive(const vector<vector<double>> &inarr,const double &lMult, TString sub = "");
        vector<double> getEventCorrelation(const vector<vector<double>> &inarr, int mOrder);
        vector<vector<double>> getEventCorrelation(const vector<vector<double>> &inarr);

        TList* GetCorrList() { return fCorrList; }
        TList* GetCumulantList(bool normalized);

        TList* GetTwoSubAnalysisList(){ return fTwoSubAnalysisList; }
        TList* GetThreeSubAnalysisList(){ return fThreeSubAnalysisList; }

        void SetEventWeight(const unsigned int &lWeight) { fEventWeight = lWeight; };
        void CalculateCorrelation();

        void RebinMulti(Int_t nbins);
        void RebinMulti(Int_t nbins, Double_t *binedges);
        TH1* getRecursiveHist(int ind, int m, unsigned int l_obs, bool sub = false);
        TH1* getSubeventCorrelation(int ind, int m);
        int getMpar() { return mpar;}
    protected:
        TList* fCorrList;
        TList* fCumulantList;
        TList* fNormList;

        // Subevent analysis
        TList* fTwoSubAnalysisList;
        TList* fThreeSubAnalysisList;

        TString fPostfix;


        const int mpar;
        unsigned int fEventWeight;
        double GetCorr(const vector<vector<double>> arr, const int m);
        double GetWeight(const vector<vector<double>> arr, const int m);

        double OrderedAddition(vector<double> vec, int size);
        TH1* getPowerHist(TH1* inh, double p);
        void FillRecursiveProfiles(const vector<double> &corr, const vector<double> &sumw, const double &lMult, TString sub);
 
        void CalculateRecursive(bool normalized);
        void CalculateCumulantHists(vector<TH1*> inh, int ind, bool normalized); 
        int factorial(const int n) { return (n<2)?1:factorial(n - 1)*n; }

      private:
        static double           fFactorial[9];
        static int              fSign[9];
    ClassDef(AliPtSubEventContainer,3);
};
#endif