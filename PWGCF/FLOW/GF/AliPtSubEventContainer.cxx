/*
Author: Frederik Kehn Rømner
Modified from AliCkContainer by Vytautas

Container to store the terms required to calculate higher order moments of the pT spectrum, 
with added sub-events inplementation. Generic m par correlation for 2-sub and up to six particle 
correlation for 3-sub. Stores all terms in expression for offline analysis.

If used, modified, or distributed, please aknowledge the original author of this code.
*/

#include "AliPtSubEventContainer.h"

double AliPtSubEventContainer::fFactorial[9] = {1.,1.,2.,6.,24.,120.,720.,5040.,40320.};
int AliPtSubEventContainer::fSign[9] = {1,-1,1,-1,1,-1,1,-1,1};

using namespace WeightSpace;

AliPtSubEventContainer::AliPtSubEventContainer():
    fTwoSubAnalysisList(0),
    fThreeSubAnalysisList(0),
    fCorrList(0),
    fCumulantList(0),
    fNormList(0),
    mpar(0),
    fEventWeight(kWeight::kWmaxperm),
    fPostfix(0)
{};
AliPtSubEventContainer::~AliPtSubEventContainer()
{
    delete fTwoSubAnalysisList;
    delete fThreeSubAnalysisList;
    delete fCorrList;
};

AliPtSubEventContainer::AliPtSubEventContainer(const char* name, const char* title, int m):
    TNamed(name,title),
    fTwoSubAnalysisList(0),
    fThreeSubAnalysisList(0),
    fCorrList(0),
    fCumulantList(0),
    fNormList(0),
    mpar(m),
    fEventWeight(kWeight::kWmaxperm),
    fPostfix(0)
{
    //Initialize(nbinsx,xlow,xhigh);
};
void AliPtSubEventContainer::Initialize(int nbinsx, double* xbins) {
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    for(int m=0;m<mpar;++m)     {
      fCorrList->Add(new TProfile(Form("corr_%ipar",m+1)+fPostfix,this->GetTitle(),nbinsx,xbins));
    }
    //printf("Container %s initialized with m = %i\n",this->GetName(),mpar);
};

void AliPtSubEventContainer::Initialize(int nbinsx, double xlow, double xhigh) {
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    for(int m=0;m<mpar;++m)
    {
      fCorrList->Add(new TProfile(Form("corr_%ipar",m+1)+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
    }
    //printf("Container %s initialized with m = %i\n",this->GetName(),mpar);
};


void AliPtSubEventContainer::InitializeTwoSub(int nbinsx, double xlow, double xhigh){
  if(fTwoSubAnalysisList) delete fTwoSubAnalysisList;
  fTwoSubAnalysisList = new TList();
  fTwoSubAnalysisList->SetOwner(kTRUE);
  fTwoSubAnalysisList->SetName("TwoSubAnalysis");

  fTwoSubAnalysisList->Add(new TList() ); ((TList*)fTwoSubAnalysisList->At(0))->SetName("Sub_A");
  fTwoSubAnalysisList->Add(new TList() ); ((TList*)fTwoSubAnalysisList->At(1))->SetName("Sub_B");
  for (int m(0); m < mpar; m++){
    ((TList*)fTwoSubAnalysisList->At(0))->Add(new TProfile(Form("twosub_corr_%ipar_sub_A",m+1)+fPostfix, this->GetTitle(),nbinsx,xlow,xhigh));
    ((TList*)fTwoSubAnalysisList->At(1))->Add(new TProfile(Form("twosub_corr_%ipar_sub_B",m+1)+fPostfix, this->GetTitle(),nbinsx,xlow,xhigh));
  }
  for (int m(0); m < mpar-1; m++){
    fTwoSubAnalysisList->Add( new TList() ); ((TList*)fTwoSubAnalysisList->At(m+2))->SetName(Form("corr_terms_%ipar",m+2));
    for(int k(0); k <= m; k++){
      ((TList*)fTwoSubAnalysisList->At(m+2))->Add(new TProfile(Form("twosub_corr_%ipar_%i%i", m+2, m+1-k, k+1) + fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
    }
  }
}

void AliPtSubEventContainer::InitializeThreeSub(int nbinsx, double xlow, double xhigh){
  if(mpar < 3){ printf("Three subevent require 3 <= %i", mpar); return;}
  
  if(fThreeSubAnalysisList) delete fThreeSubAnalysisList;
  fThreeSubAnalysisList = new TList();
  fThreeSubAnalysisList->SetOwner(kTRUE);
  fThreeSubAnalysisList->SetName("ThreeSubAnalysis");

  fThreeSubAnalysisList->Add(new TList() ); ((TList*)fThreeSubAnalysisList->At(0))->SetName("Sub_A");
  fThreeSubAnalysisList->Add(new TList() ); ((TList*)fThreeSubAnalysisList->At(1))->SetName("Sub_B");
  fThreeSubAnalysisList->Add(new TList() ); ((TList*)fThreeSubAnalysisList->At(2))->SetName("Sub_C");

  for (int m(0); m < mpar; m++){
    ((TList*)fThreeSubAnalysisList->At(0))->Add(new TProfile(Form("threesub_corr_%ipar_sub_A",m+1)+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
    ((TList*)fThreeSubAnalysisList->At(1))->Add(new TProfile(Form("threesub_corr_%ipar_sub_B",m+1)+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
    ((TList*)fThreeSubAnalysisList->At(2))->Add(new TProfile(Form("threesub_corr_%ipar_sub_C",m+1)+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  }

  // three-par correlation
  fThreeSubAnalysisList->Add( new TList() ); ((TList*)fThreeSubAnalysisList->At(3))->SetName("corr_terms_3par");
  ((TList*)fThreeSubAnalysisList->At(3))->Add(new TProfile("threesub_corr_3par_111"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  
  // four-par correlation    
  if(mpar < 4) return;
  fThreeSubAnalysisList->Add( new TList() ); ((TList*)fThreeSubAnalysisList->At(4))->SetName("corr_terms_4par");
  ((TList*)fThreeSubAnalysisList->At(4))->Add(new TProfile("threesub_corr_4par_112"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(4))->Add(new TProfile("threesub_corr_4par_121"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(4))->Add(new TProfile("threesub_corr_4par_211"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));

  // five-par correlation    
  if(mpar < 5) return;
  fThreeSubAnalysisList->Add( new TList() ); ((TList*)fThreeSubAnalysisList->At(5))->SetName("corr_terms_5par");
  ((TList*)fThreeSubAnalysisList->At(5))->Add(new TProfile("threesub_corr_5par_113"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(5))->Add(new TProfile("threesub_corr_5par_131"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(5))->Add(new TProfile("threesub_corr_5par_311"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(5))->Add(new TProfile("threesub_corr_5par_122"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(5))->Add(new TProfile("threesub_corr_5par_212"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(5))->Add(new TProfile("threesub_corr_5par_221"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));

  // six-par correlation    
  if(mpar < 6) return;
  fThreeSubAnalysisList->Add( new TList() ); ((TList*)fThreeSubAnalysisList->At(6))->SetName("corr_terms_6par");
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_114"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_141"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_411"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_123"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_132"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_213"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_231"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_312"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_321"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
  ((TList*)fThreeSubAnalysisList->At(6))->Add(new TProfile("threesub_corr_6par_222"+fPostfix,this->GetTitle(),nbinsx,xlow,xhigh));
};



vector<vector<double>> AliPtSubEventContainer::getEventCorrelation(const vector<vector<double>> &inarr) {
  vector<double> corr(mpar+1,0.0); corr[0] = 1.0;
  vector<double> sumw(mpar+1,0.0); sumw[0] = 1.0;
  double sumNum = 0;
  double sumDenum = 0;
  vector<double> valNum;
  vector<double> valDenum;
  for(int m(1); m<=mpar; ++m)
  {
    for(int k(1);k<=m;++k)
    {
      valNum.push_back(fSign[k-1]*corr[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][k]);
      valDenum.push_back(fSign[k-1]*sumw[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][0]);
    }
    sumNum = OrderedAddition(valNum, m);
    sumDenum = OrderedAddition(valDenum, m);
    valNum.clear();
    valDenum.clear();
  
    corr[m] = sumNum;
    sumw[m] = sumDenum;
  }
  vector<vector<double>> outvec = {corr,sumw};
  return outvec;
}

vector<double> AliPtSubEventContainer::getEventCorrelation(const vector<vector<double>> &inarr, int mOrder) {
  vector<double> corr(mpar+1,0.0); corr[0] = 1.0;
  vector<double> sumw(mpar+1,0.0); sumw[0] = 1.0;
  double sumNum = 0;
  double sumDenum = 0;
  vector<double> valNum;
  vector<double> valDenum;
  for(int m(1); m<=mpar; ++m)
  {
    for(int k(1);k<=m;++k)
    {
      valNum.push_back(fSign[k-1]*corr[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][k]);
      valDenum.push_back(fSign[k-1]*sumw[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][0]);
    }
    sumNum = OrderedAddition(valNum, m);
    sumDenum = OrderedAddition(valDenum, m);
    valNum.clear();
    valDenum.clear();
  
    corr[m] = sumNum;
    sumw[m] = sumDenum;
  }
  vector<double> outvec = {corr[mOrder],sumw[mOrder]};
  return outvec;
}
void AliPtSubEventContainer::FillRecursive(const vector<vector<double>> &inarr,const double &lMult, TString sub) {
  vector<double> corr(mpar+1,0.0); corr[0] = 1.0;
  vector<double> sumw(mpar+1,0.0); sumw[0] = 1.0;
  double sumNum = 0;
  double sumDenum = 0;
  vector<double> valNum;
  vector<double> valDenum;

  for(int m(1); m<=mpar; ++m)
  {
    for(int k(1);k<=m;++k)
    {
      valNum.push_back(fSign[k-1]*corr[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][k]);
      valDenum.push_back(fSign[k-1]*sumw[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][0]);
    }
    sumNum = OrderedAddition(valNum, m);
    sumDenum = OrderedAddition(valDenum, m);
    valNum.clear();
    valDenum.clear();
  
    corr[m] = sumNum;
    sumw[m] = sumDenum;
  }
  FillRecursiveProfiles(corr,sumw,lMult,sub);
  return;
}

void AliPtSubEventContainer::FillRecursiveProfiles(const vector<double> &corr, const vector<double> &sumw, const double &lMult, TString sub) {
    for(int m=1;m<=mpar;++m)
    {
      if(sumw[m]==0) continue; 
      ((TProfile*)fCorrList->At(m-1))->Fill(lMult,corr[m]/sumw[m],(fEventWeight==WeightSpace::kOne)?1.0:sumw[m]);
    }
    return;
}

void AliPtSubEventContainer::FillTwoSubAnalsysis(const vector<vector<double>> &inarrSubA, const vector<vector<double>> &inarrSubB, const double &lMult){
  vector<vector<double>> arrA = getEventCorrelation(inarrSubA);
  vector<vector<double>> arrB = getEventCorrelation(inarrSubB);
  for (int m(1); m <= mpar; m++){
    ((TProfile*)((TList*)fTwoSubAnalysisList->At(0))->At(m-1))->Fill(lMult,arrA[0][m]/arrA[1][m],(fEventWeight==WeightSpace::kOne)?1.0:arrA[1][m]);
    ((TProfile*)((TList*)fTwoSubAnalysisList->At(1))->At(m-1))->Fill(lMult,arrB[0][m]/arrB[1][m],(fEventWeight==WeightSpace::kOne)?1.0:arrB[1][m]);
  }
  double prod, prodw;
  for (int m(1); m <= mpar-1; m++){
    for(int k(1); k <= m; k++){
      prod =  (arrA[0][1+m-k]/arrA[1][1+m-k])*(arrB[0][k]/arrB[1][k]);
      prodw = arrA[1][1+m-k]*arrB[1][k];
      ((TProfile*)((TList*)fTwoSubAnalysisList->At(m-1+2))->At(k-1))->Fill(lMult,prod,(fEventWeight==WeightSpace::kOne)?1.0:prodw);
    }
  }
}

void AliPtSubEventContainer::FillThreeSubAnalsysis(const vector<vector<double>> &inarrSubA, const vector<vector<double>> &inarrSubB, const vector<vector<double>> &inarrSubC, const double &lMult){
  // arr[0][m] -> Get correlations      arr[1][m] -> get weights
  vector<vector<double>> arrA = getEventCorrelation(inarrSubA);
  vector<vector<double>> arrB = getEventCorrelation(inarrSubB);
  vector<vector<double>> arrC = getEventCorrelation(inarrSubC);
  // fill <pt> in sub-ranges
  for (int m(1); m <= mpar; m++){
    ((TProfile*)((TList*)fThreeSubAnalysisList->At(0))->At(m-1))->Fill(lMult,arrA[0][m]/arrA[1][m],(fEventWeight==WeightSpace::kOne)?1.0:arrA[1][m]);
    ((TProfile*)((TList*)fThreeSubAnalysisList->At(1))->At(m-1))->Fill(lMult,arrB[0][m]/arrB[1][m],(fEventWeight==WeightSpace::kOne)?1.0:arrB[1][m]);
    ((TProfile*)((TList*)fThreeSubAnalysisList->At(2))->At(m-1))->Fill(lMult,arrC[0][m]/arrC[1][m],(fEventWeight==WeightSpace::kOne)?1.0:arrC[1][m]);
  }
  // fill three-par corr
  if(mpar < 3) return;
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(3))->At(0))->Fill(lMult,(arrA[0][1]/arrA[1][1])*(arrB[0][1]/arrB[1][1])*(arrC[0][1]/arrC[1][1]),(fEventWeight==WeightSpace::kOne)?1.0:(arrA[1][1]*arrB[1][1]*arrC[1][1]));
  // fill four-par corr
  if(mpar < 4) return;
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(4))->At(0))->Fill(lMult,(arrA[0][1]/arrA[1][1])*(arrB[0][1]/arrB[1][1])*(arrC[0][2]/arrC[1][2]),(fEventWeight==WeightSpace::kOne)?1.0:(arrA[1][1]*arrB[1][1]*arrC[1][2]));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(4))->At(1))->Fill(lMult,(arrA[0][1]/arrA[1][1])*(arrB[0][2]/arrB[1][2])*(arrC[0][1]/arrC[1][1]),(fEventWeight==WeightSpace::kOne)?1.0:(arrA[1][1]*arrB[1][2]*arrC[1][1]));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(4))->At(2))->Fill(lMult,(arrA[0][2]/arrA[1][2])*(arrB[0][1]/arrB[1][1])*(arrC[0][1]/arrC[1][1]),(fEventWeight==WeightSpace::kOne)?1.0:(arrA[1][2]*arrB[1][1]*arrC[1][1]));
  // fill five-par corr
  if(mpar < 5) return;
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(5))->At(0))->Fill(lMult, GetCorr(arrA,1)*GetCorr(arrB,1)*GetCorr(arrC,3),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,1)*GetWeight(arrB,1)*GetWeight(arrC,3)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(5))->At(1))->Fill(lMult, GetCorr(arrA,1)*GetCorr(arrB,3)*GetCorr(arrC,1),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,1)*GetWeight(arrB,3)*GetWeight(arrC,1)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(5))->At(2))->Fill(lMult, GetCorr(arrA,3)*GetCorr(arrB,1)*GetCorr(arrC,1),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,3)*GetWeight(arrB,1)*GetWeight(arrC,1)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(5))->At(3))->Fill(lMult, GetCorr(arrA,1)*GetCorr(arrB,2)*GetCorr(arrC,2),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,1)*GetWeight(arrB,2)*GetWeight(arrC,2)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(5))->At(4))->Fill(lMult, GetCorr(arrA,2)*GetCorr(arrB,1)*GetCorr(arrC,2),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,2)*GetWeight(arrB,1)*GetWeight(arrC,2)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(5))->At(5))->Fill(lMult, GetCorr(arrA,2)*GetCorr(arrB,2)*GetCorr(arrC,1),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,2)*GetWeight(arrB,2)*GetWeight(arrC,1)));
  // fill six-par corr
  if(mpar < 6) return;
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(0))->Fill(lMult, GetCorr(arrA,1)*GetCorr(arrB,1)*GetCorr(arrC,4),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,1)*GetWeight(arrB,1)*GetWeight(arrC,4)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(1))->Fill(lMult, GetCorr(arrA,1)*GetCorr(arrB,4)*GetCorr(arrC,1),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,1)*GetWeight(arrB,4)*GetWeight(arrC,1)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(2))->Fill(lMult, GetCorr(arrA,4)*GetCorr(arrB,1)*GetCorr(arrC,1),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,4)*GetWeight(arrB,1)*GetWeight(arrC,1)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(3))->Fill(lMult, GetCorr(arrA,1)*GetCorr(arrB,2)*GetCorr(arrC,3),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,1)*GetWeight(arrB,2)*GetWeight(arrC,3)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(4))->Fill(lMult, GetCorr(arrA,1)*GetCorr(arrB,3)*GetCorr(arrC,2),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,1)*GetWeight(arrB,3)*GetWeight(arrC,2)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(5))->Fill(lMult, GetCorr(arrA,2)*GetCorr(arrB,1)*GetCorr(arrC,3),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,2)*GetWeight(arrB,3)*GetWeight(arrC,1)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(6))->Fill(lMult, GetCorr(arrA,2)*GetCorr(arrB,3)*GetCorr(arrC,1),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,2)*GetWeight(arrB,1)*GetWeight(arrC,3)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(7))->Fill(lMult, GetCorr(arrA,3)*GetCorr(arrB,1)*GetCorr(arrC,2),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,3)*GetWeight(arrB,2)*GetWeight(arrC,1)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(8))->Fill(lMult, GetCorr(arrA,3)*GetCorr(arrB,2)*GetCorr(arrC,1),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,3)*GetWeight(arrB,1)*GetWeight(arrC,2)));
  ((TProfile*)((TList*)fThreeSubAnalysisList->At(6))->At(9))->Fill(lMult, GetCorr(arrA,2)*GetCorr(arrB,2)*GetCorr(arrC,2),(fEventWeight==WeightSpace::kOne)?1.0:(GetWeight(arrA,2)*GetWeight(arrB,2)*GetWeight(arrC,2)));
}


TH1* AliPtSubEventContainer::getRecursiveHist(int ind, int m, unsigned int l_obs, bool sub) {
  if(l_obs==kObs::kCum) 
  {
      if(!fCumulantList) CalculateRecursive(false);
      if(!fCumulantList) return 0;
      if(ind+1<fCumulantList->GetEntries()) return (TH1*)fCumulantList->At((ind+1)*mpar+m-1);
  } 
  if(l_obs==kObs::kNorm)
  {
      if(!fNormList) CalculateRecursive(true);
      if(!fNormList) return 0;
      if(ind+1<fNormList->GetEntries()) return (TH1*)fNormList->At((ind+1)*mpar+m-1);
  }
  return 0;
}

void AliPtSubEventContainer::CalculateRecursive(bool normalized) { 
    if(normalized)
    {
        if(fNormList) delete fNormList;
        fNormList = new TList();
        fNormList->SetOwner(); 
        if(!fCumulantList)
        {
            fCumulantList = new TList();
            fCumulantList->SetOwner(); 
            printf("cumulant list created!\n") ;
        }
    }
    else
    {
        if(fCumulantList) delete fCumulantList;
        fCumulantList = new TList();
        fCumulantList->SetOwner();
    }
    vector<TH1*> hTs;
    for(int j=0;j<mpar;++j) {
      ((TProfile*)fCorrList->At(j))->SetErrorOption("g");
    }
    CalculateCumulantHists(hTs,-1,normalized);
    return;
}

void AliPtSubEventContainer::CalculateCumulantHists(vector<TH1*> inh, int ind, bool normalized) {
    auto binomial = [&](const int n, const int m) { return factorial(n)/(factorial(m)*factorial(n-m)); };
    int lMax = mpar;
    if((normalized && fCumulantList->GetEntries()<lMax) || !normalized)
    {
        for(int m=1;m<=mpar;++m)
        {
            TH1* reth = (TH1*)inh[m-1]->Clone(Form("h%i%i_%i",m,ind,normalized));
            TH1* hWeights = (TH1*)inh[m-1]->Clone(Form("hWeights%i%i_%i",m,ind,normalized));
            for(int k=1;k<m;++k)
            {
                TH1* corrh = (TH1*)inh[m-k-1]->Clone(Form("h%i%i%i_%i",m,k,ind,normalized));
                corrh->Multiply((TH1*)fCumulantList->At((ind+1)*mpar+k-1));
                corrh->Scale(binomial(m-1,k-1));
                reth->Add(corrh,-1);
                delete corrh;
            }
            for(int i=1;i<=hWeights->GetNbinsX();++i) reth->SetBinError(i,hWeights->GetBinError(i));
            delete hWeights;
            fCumulantList->Add(reth->Clone(Form("reth%i%i",m,ind)));
        }
    }
    if(normalized)
    {
        TH1* sigmah = (TH1*)fCumulantList->At((ind+1)*mpar+1)->Clone(Form("sigmah%i",ind));
        for(int m = 1; m<=mpar;++m)
        {
            TH1* normh = (TH1*)fCumulantList->At((ind+1)*mpar+m-1)->Clone(Form("normh%i%i",m,ind));
            normh->Divide(getPowerHist(sigmah,0.5*(m)));
            fNormList->Add(normh);
        }
    }

    return;
}

TH1* AliPtSubEventContainer::getPowerHist(TH1* inh, double p) {
    TH1D* reth = (TH1D*)inh->Clone("reth");
    reth->SetName(Form("power%.2f_%s",p,inh->GetName()));
    for(int i=1;i<=inh->GetNbinsX();i++) {
        if(inh->GetBinContent(i)>=0 || std::floor(p) == p)
        {
            reth->SetBinContent(i,pow(inh->GetBinContent(i),p));
            reth->SetBinError(i,p*pow(reth->GetBinContent(i),p-1)*inh->GetBinError(i));
        }
        else
        {
            reth->SetBinContent(i,-999);
            reth->SetBinError(i,0.000000001);
        }

    }   
    return reth;
}


TList *AliPtSubEventContainer::GetCumulantList(bool normalized)
{
  if(normalized){
    CalculateRecursive(true);
    return fNormList;
  }
  CalculateRecursive(false);
  return fCumulantList;
}

double AliPtSubEventContainer::GetCorr(const vector<vector<double>> arr, const int m){
  if(mpar < m){ printf("Could not get correlation order m=%i", m); return 1; }
  return arr[0][m]/arr[1][m];
}

double AliPtSubEventContainer::GetWeight(const vector<vector<double>> arr, const int m){
  if(mpar < m){ printf("Could not get weight order m=%i", m); return 1; }
  return arr[1][m];
}

double AliPtSubEventContainer::OrderedAddition(vector<double> vec, int size) {
  double sum = 0;
  std::sort(vec.begin(), vec.end());
  for(int i = 0; i < size; i++)
  {
    sum += vec[i];
  }
  return sum;
}