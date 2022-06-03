#ifndef ALIPTCONTAINER__H
#define ALIPTCONTAINER__H

#include "AliProfileBS.h"
#include "TNamed.h"
#include "TList.h"
#include "TCollection.h"

using namespace std;

namespace PtSpace {
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

class AliPtContainer: public TNamed {
    public:
        AliPtContainer();
        ~AliPtContainer();
        AliPtContainer(const char* name, const char* title, int nbinsx, double* xbins, int m, bool sub);
        AliPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, int m, bool sub);
        void Initialize(int nbinsx, double* xbins);
        void Initialize(int nbinsx, double xlow, double xhigh);
        void InitializeSubsamples(const int &nsub);
        void FillCk(const vector<vector<double>> &inarr, const double &lMult, const double &rn);
        void FillSkew(const vector<vector<double>> &inarr, const double &lMult, const double &rn);
        void FillKurtosis(const vector<vector<double>> &inarr, const double &lMult, const double &rn);
        void FillRecursive(const vector<vector<double>> &inarr,const double &lMult, const double &rn, TString sub = "");
        vector<double> getEventCorrelation(const vector<vector<double>> &inarr, int mOrder);
        TList* GetCkList() { return fCkTermList; }
        TList* GetSkewList() { return fSkewTermList; }
        TList* GetKurtosisList() { return fKurtosisTermList; }
        TList* GetCorrList() { return fCorrList; }
        void SetEventWeight(const unsigned int &lWeight) { fEventWeight = lWeight; };
        void CalculateCorrelation();
        TH1* getCkHist(int ind);
        TH1* getSkewHist(int ind);
        TH1* getKurtosisHist(int ind);
        void RebinMulti(Int_t nbins);
        void RebinMulti(Int_t nbins, Double_t *binedges);
        TH1* getRecursiveHist(int ind, int m, unsigned int l_obs, bool sub = false);
        TH1* getSubeventCorrelation(int ind, int m);
        int getMpar() { return mpar;}
        Long64_t Merge(TCollection *collist);
    protected:
        TList* fCkTermList;
        TList* fSkewTermList;
        TList* fKurtosisTermList;
        TList* fCkList;
        TList* fSkewList;
        TList* fKurtosisList;
        TList* fCorrList;
        TList* fSubList;
        TList* fCumulantList;
        TList* fNormList;
        const int mpar;
        unsigned int fEventWeight;
        bool fSubevent;
        double OrderedAddition(vector<double> vec, int size);
        TH1* getPowerHist(TH1* inh, double p);
        void FillRecursiveProfiles(const vector<double> &corr, const vector<double> &sumw, const double &lMult, const double &rn, TString sub);
        void CalculateCk();
        void CalculateSkew();
        void CalculateKurtosis();
        TH1* RecalculateCkHists(vector<TH1*> inh);
        TH1* RecalculateSkewHists(vector<TH1*> inh);
        TH1* RecalculateKurtosisHists(vector<TH1*> inh);
        void CalculateCumulantHists(vector<TH1*> inh, int ind, bool normalized); 
        void CalculateRecursive(bool normalized);
        int factorial(const int n) { return (n<2)?1:factorial(n - 1)*n; }

      private:
        static double           fFactorial[9];
        static int              fSign[9];
        void MergeBSLists(TList *source, TList *target);
    ClassDef(AliPtContainer,1);
};
#endif