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
        kWperms
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
        AliPtContainer(const char* name, const char* title, Int_t nbinsx, Double_t* xbins, Int_t m, Bool_t sub);
        AliPtContainer(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t m, Bool_t sub);
        void Initialize(Int_t nbinsx, Double_t* xbins);
        void Initialize(Int_t nbinsx, Double_t xlow, Double_t xhigh);
        void InitializeSubsamples(const Int_t &nsub);
        void FillCk(const vector<vector<Double_t>> &inarr, const Double_t &lMult, const Double_t &rn);
        void FillSkew(const vector<vector<Double_t>> &inarr, const Double_t &lMult, const Double_t &rn);
        void FillKurtosis(const vector<vector<Double_t>> &inarr, const Double_t &lMult, const Double_t &rn);
        void FillRecursive(const vector<vector<Double_t>> &inarr, Int_t subIndex = 0);
        void FillRecursiveProfiles(const Double_t &lMult, const Double_t &rn, bool exotic = kFALSE);
        vector<Double_t> getEventCorrelation(Int_t mOrder, Int_t subIndex = 0);
        double getExoticEventCorrelation(int w, int p);
        void FillExotic(int mMax, const vector<vector<double>> &inarr);
        TList* GetCkList() { return fCkTermList; }
        TList* GetSkewList() { return fSkewTermList; }
        TList* GetKurtosisList() { return fKurtosisTermList; }
        TList* GetCorrList() { return fCorrList; }
        TList* GetSubList() { return fSubList; }
        void SetEventWeight(const unsigned int &lWeight) { fEventWeight = lWeight; };
        void CalculateCorrelation();
        TH1* getCkHist(Int_t ind);
        TH1* getSkewHist(Int_t ind);
        TH1* getKurtosisHist(Int_t ind);
        void RebinMulti(Int_t nbins);
        void RebinMulti(Int_t nbins, Double_t *binedges);
        TH1* getRecursiveHist(Int_t ind, Int_t m, unsigned int l_obs, Bool_t sub = false);
        TH1* getSubeventCumulantHist(int ind, int m);
        vector<Double_t> getSubeventCorrelation(Int_t mOrder);
        Int_t getMpar() { return mpar;}
        Long64_t Merge(TCollection *collist);
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
        Bool_t fSubevent;
        vector<vector<Double_t>> fCorr; 
        vector<vector<Double_t>> fSumw; 
        vector<vector<Double_t>> fExoticCorr; 
        Double_t OrderedAddition(vector<Double_t> vec);
        TH1* getPowerHist(TH1* inh, Double_t p);
        void CalculateCk();
        void CalculateSkew();
        void CalculateKurtosis();
        TH1* RecalculateCkHists(vector<TH1*> inh);
        TH1* RecalculateSkewHists(vector<TH1*> inh);
        TH1* RecalculateKurtosisHists(vector<TH1*> inh);
        void CalculateCumulantHists(vector<TH1*> inh, Int_t ind, Bool_t normalized); 
        void CalculateRecursive(Bool_t normalized);
        Int_t getSubIndex(TString sub);
        Int_t factorial(const Int_t n) { return (n<2)?1:factorial(n - 1)*n; }

      private:
        static Double_t           fFactorial[9];
        static Int_t              fSign[9];
        static Double_t           fCoeff[5][5][5][5];
        void MergeBSLists(TList *source, TList *target);
    ClassDef(AliPtContainer,1);
};
#endif