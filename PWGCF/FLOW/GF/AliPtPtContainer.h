#ifndef ALIPTPTCONTAINER__H
#define ALIPTPTCONTAINER__H

#include "AliProfileBS.h"
#include "TNamed.h"
#include "TList.h"
#include "TCollection.h"

using namespace std;

namespace PtPtSpace {
    enum kEventWeight {
        kUnity,
        kTuples
    };
};

class AliPtPtContainer: public TNamed {
    public:
        AliPtPtContainer();
        ~AliPtPtContainer();
        AliPtPtContainer(const char* name, const char* title, Int_t nbinsx, Double_t* xbins, Int_t m);
        AliPtPtContainer(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t m);
        void Initialize(Int_t nbinsx, Double_t* xbins);
        void Initialize(Int_t nbinsx, Double_t xlow, Double_t xhigh);
        void InitializeSubsamples(const Int_t &nsub);
        void CalculateCorrelations(const vector<vector<double>> &inarr);
        void CalculateCMTerms(const vector<vector<double>> &inarr);
        void FillProfiles(const Double_t &centmult, const Double_t &rn);
        void FillCMProfiles(const vector<vector<double>> &inarr, const double &centmult, const double &rn);
        vector<Double_t> getEventCorrelation(Int_t mOrder);
        TList* GetCorrList() { return fCorrList; }
        TList* GetCMTermList() { return fCMTermList; }
        void SetEventWeight(const unsigned int &lWeight) { fEventWeight = lWeight; };
        void RebinMulti(Int_t nbins);
        void RebinMulti(Int_t nbins, Double_t *binedges);
        TH1* getCentralMomentHist(Int_t ind, Int_t m);
        TH1* getCumulantHist(Int_t ind, Int_t m);
        TH1* getCorrHist(Int_t ind, Int_t m);
        double getNumerator(const int m) { return fCorr[m]; }
        double getDenominator(const int m) { return fSumw[m]; }
        double getNumeratorCM(const int m, const int t) { return fcmNum[m * (m - 1) / 2 + (m - t - 1)]; }
        double getDenominatorCM(const int m) { return fcmDen[m-1]; }
        Int_t getMpar() { return mpar;}
        Long64_t Merge(TCollection *collist);
        TList* fCMTermList;
        TList* fCorrList;
        TList* fCumulantList;
        TList* fCentralMomentList;
        const int mpar;
        unsigned int fEventWeight;
        vector<Double_t> fCorr;
        vector<Double_t> fSumw;
        vector<Double_t> fcmNum;
        vector<Double_t> fcmDen;
        Double_t OrderedAddition(vector<Double_t> vec);
        void CreateCentralMomentList();
        void CalculateCentralMomentHists(vector<TH1*> inh, int ind, int m, TH1* hMpt);
        void CreateCumulantList();
        void CalculateCumulantHists(vector<TH1*> inh, Int_t ind);
        Int_t factorial(const Int_t n) { return (n<2)?1:factorial(n - 1)*n; }

      private:
        static Double_t           fFactorial[9];
        static Int_t              fSign[9];
        static Double_t           fCoeff[5][5][5][5];
        void MergeBSLists(TList *source, TList *target);
        TH1* raiseHistToPower(TH1* inh, double p);
    ClassDef(AliPtPtContainer,4);
};
#endif
