#ifndef ALIPTCONTAINER
#define ALIPTCONTAINER

#include "AliProfileBS.h"
#include "TNamed.h"
#include "TList.h"

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
        template<typename T> void FillObs(T& inarr, const double &lMult, const double &rn);
        TList* GetTermList() { return fTermList; }
        void CalculateCorrelation();
        TH1* getHist(int ind);
    protected:
        TList* fTermList;
        TList* fObsList;
        int mpar;
        template<typename T> void FillMpt(T& inarr, const double &lMult, const double &rn);
        template<typename T> void FillCk(T& inarr, const double &lMult, const double &rn);
        template<typename T> void FillSkew(T& inarr, const double &lMult, const double &rn);
        template<typename T> void FillKurtosis(T& inarr, const double &lMult, const double &rn);
        void CalculateObs();
        TH1* RecalculateObsHists(vector<TH1*> inh); 
        TH1* RecalculateCkHists(vector<TH1*> inh);
        TH1* RecalculateSkewHists(vector<TH1*> inh);
        TH1* RecalculateKurtosisHists(vector<TH1*> inh);
    ClassDef(AliPtContainer,1);
};
#endif
