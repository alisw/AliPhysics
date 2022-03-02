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
        void FillObs(double inarr[10][10], const double &lMult, const double &rn);
        void FillRecursive(double inarr[10][10],const double &lMult, const double &rn, TString sub = "");
        TList* GetTermList() { return fTermList; }
        TList* GetCorrList() { return fCorrList; }
        void SetEventWeight(const unsigned int &lWeight) { fEventWeight = lWeight; };
        void CalculateCorrelation();
        TH1* getHist(int ind);
        TH1* getAltHist(int ind);
        void RebinMulti(Int_t nbins);
        void RebinMulti(Int_t nbins, Double_t *binedges);
        TH1* getRercusiveHist(int ind, int m, unsigned int l_obs, bool sub = false);
        Long64_t Merge(TCollection *collist);
    protected:
        TList* fTermList;
        TList* fCorrList;
        TList* fSubList;
        TList* fMomentList;
        TList* fCumulantList;
        TList* fNormList;
        int mpar;
        unsigned int fEventWeight;
        bool fSubevent;
        double OrderedAddition(vector<double> vec, int size);
        TH1* getPowerHist(TH1* inh, double p);
        void FillMpt(double inarr[10][10], const double &lMult, const double &rn);
        void FillCk(double inarr[10][10], const double &lMult, const double &rn);
        void FillSkew(double inarr[10][10], const double &lMult, const double &rn);
        void FillKurtosis(double inarr[10][10], const double &lMult, const double &rn);
        void FillRecursiveProfiles(const vector<double> &corr, const vector<double> &sumw, const double &lMult, const double &rn, TString sub);
        void CalculateObs();
        TH1* RecalculateObsHists(vector<TH1*> inh);
        TH1* RecalculateCkHists(vector<TH1*> inh);
        TH1* RecalculateSkewHists(vector<TH1*> inh);
        TH1* RecalculateKurtosisHists(vector<TH1*> inh);
        TH1* getSubeventCorrelation(int ind, int m);
        void CalculateCumulantHists(vector<TH1*> inh, int ind, bool normalized); 
        void CalculateRecursive(bool normalized);
        int factorial(const int n) { return (n<2)?1:factorial(n - 1)*n; }
        int GetWeightIndex(int m);

      private:
        static double           fFactorial[9];
        static int              fSign[9];
        void MergeBSLists(TList *source, TList *target);
    ClassDef(AliPtContainer,1);
};
#endif