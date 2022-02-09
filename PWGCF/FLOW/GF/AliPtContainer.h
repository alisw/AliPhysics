#ifndef ALIPTCONTAINER__H
#define ALIPTCONTAINER__H

#include "AliProfileBS.h"
#include "TNamed.h"
#include "TList.h"
#include "TCollection.h"

using namespace std;
class AliPtContainer: public TNamed {
    public:
        AliPtContainer();
        ~AliPtContainer();
        AliPtContainer(const char* name, const char* title, int nbinsx, double* xbins, int m);
        AliPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, int m);
        void Initialize(int nbinsx, double* xbins);
        void Initialize(int nbinsx, double xlow, double xhigh);
        void InitializeSubsamples(const int &nsub);
        void FillObs(double inarr[15][15], const double &lMult, const double &rn);
        TList* GetTermList() { return fTermList; }
        void CalculateCorrelation();
        TH1* getHist(int ind);
        Long64_t Merge(TCollection *collist);
    protected:
        TList* fTermList;
        TList* fObsList;
        int mpar;
        void FillMpt(double inarr[15][15], const double &lMult, const double &rn);
        void FillCk(double inarr[15][15], const double &lMult, const double &rn);
        void FillSkew(double inarr[15][15], const double &lMult, const double &rn);
        void FillKurtosis(double inarr[15][15], const double &lMult, const double &rn);
        void CalculateObs();
        TH1* RecalculateObsHists(vector<TH1*> inh);
        TH1* RecalculateCkHists(vector<TH1*> inh);
        TH1* RecalculateSkewHists(vector<TH1*> inh);
        TH1* RecalculateKurtosisHists(vector<TH1*> inh);
      private:
        void MergeBSLists(TList *source, TList *target);
    ClassDef(AliPtContainer,1);
};
#endif
