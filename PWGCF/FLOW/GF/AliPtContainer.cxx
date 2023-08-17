/*
Author: Emil Gorm Nielsen
Modified from AliCkContainer by Vytautas
Container to store the terms required to calculate higher order moments of the pT spectrum.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#include "AliPtContainer.h"

Double_t AliPtContainer::fFactorial[9] = {1.,1.,2.,6.,24.,120.,720.,5040.,40320.};
Int_t AliPtContainer::fSign[9] = {1,-1,1,-1,1,-1,1,-1,1};
Double_t AliPtContainer::fCoeff[5][5][5][5] = {{{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}}},
                                          {{{0,0,0,0,0},{1,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{1,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}}},
                                          {{{0,0,0,0,0},{0,1,0,0,0},{1,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}}},
                                          {{{0,0,0,0,0},{0,0,1,0,0},{0,1,0,0,0},{1,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,2,-1,0,0},{1,0,0,0,0},{1,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{0,0,0,0,0}},{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}}},
                                          {{{0,0,0,0,0},{0,0,0,1,0},{0,0,1,0,0},{0,1,0,0,0},{1,0,0,0,0}},{{0,0,0,0,0},{0,0,3,-2,0},{0,2,-1,0,0},{1,0,0,0,0},{1,0,0,0,0}},{{0,0,0,0,0},{0,2,-1,0,0},{2./3,2./3,-1./3,0,0},{1,0,0,0,0},{1,0,0,0,0}},{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0}},{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0}}}};
using namespace PtSpace;
AliPtContainer::AliPtContainer():
    fCkTermList(0),
    fSkewTermList(0),
    fKurtosisTermList(0),
    fCkList(0),
    fSkewList(0),
    fKurtosisList(0),
    fCorrList(0),
    fSubList(0),
    fCumulantList(0),
    fNormList(0),
    mpar(0),
    fEventWeight(kWeight::kOne),
    fSubevent(false),
    fCorr(),
    fSumw(),
    fExoticCorr()
{};
AliPtContainer::~AliPtContainer()
{
    delete fCorrList;
    delete fCkTermList;
    delete fSkewTermList;
    delete fKurtosisTermList;
};
AliPtContainer::AliPtContainer(const char* name, const char* title, int nbinsx, double* xbins, int m, bool sub):
    TNamed(name,title),
    fCkTermList(0),
    fSkewTermList(0),
    fKurtosisTermList(0),
    fCkList(0),
    fSkewList(0),
    fKurtosisList(0),
    fCorrList(0),
    fSubList(0),
    fCumulantList(0),
    fNormList(0),
    mpar(m),
    fEventWeight(kWeight::kWperms),
    fSubevent(sub),
    fCorr(),
    fSumw(),
    fExoticCorr()
{
    Initialize(nbinsx,xbins);
};
AliPtContainer::AliPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, int m, bool sub):
    TNamed(name,title),
    fCkTermList(0),
    fSkewTermList(0),
    fKurtosisTermList(0),
    fCkList(0),
    fSkewList(0),
    fKurtosisList(0),
    fCorrList(0),
    fSubList(0),
    fCumulantList(0),
    fNormList(0),
    mpar(m),
    fEventWeight(kWeight::kWperms),
    fSubevent(sub),
    fCorr(),
    fSumw(),
    fExoticCorr()
{
    Initialize(nbinsx,xlow,xhigh);
};
void AliPtContainer::Initialize(int nbinsx, double* xbins) {
    fCorr.resize(3);
    fSumw.resize(3);
    if(fCkTermList) delete fCkTermList;
    fCkTermList = new TList();
    fCkTermList->SetOwner(kTRUE);
    if(fSkewTermList) delete fSkewTermList;
    fSkewTermList = new TList();
    fSkewTermList->SetOwner(kTRUE);
    if(fKurtosisTermList) delete fKurtosisTermList;
    fKurtosisTermList = new TList();
    fKurtosisTermList->SetOwner(kTRUE);
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    if(fSubevent){
        if(fSubList) delete fSubList;
        fSubList = new TList();
        fSubList->SetOwner(kTRUE);
    }
    for(int i=0;i<=2;++i){
        fCkTermList->Add(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xbins));
    }
    for(int i=0;i<=3;++i){
        fSkewTermList->Add(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xbins));
    }
    for(int i=0;i<=4;++i){
        fKurtosisTermList->Add(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xbins));
    }
  for(int m=0;m<mpar;++m){
    if(m<4){
      for(int i(m);i<4;++i) {
        fCorrList->Add(new AliProfileBS(Form("corr_%ipar_%ipc",m+1,i+1),this->GetTitle(),nbinsx,xbins));
      }
    }
    else fCorrList->Add(new AliProfileBS(Form("corr_%ipar",m+1),this->GetTitle(),nbinsx,xbins));
    if(fSubevent) {
        fSubList->Add(new AliProfileBS(Form("corr_%ipar_subP",m+1),this->GetTitle(),nbinsx,xbins));
        fSubList->Add(new AliProfileBS(Form("corr_%ipar_subN",m+1),this->GetTitle(),nbinsx,xbins));
        fSubList->Add(new AliProfileBS(Form("corr_%ipar_sub",m+1),this->GetTitle(),nbinsx,xbins));
    }
  }
  printf("Container %s initialized with m = %i\n",this->GetName(),mpar);
};
void AliPtContainer::Initialize(int nbinsx, double xlow, double xhigh) {
    fCorr.resize(3);
    fSumw.resize(3);
    if(fCkTermList) delete fCkTermList;
    fCkTermList = new TList();
    fCkTermList->SetOwner(kTRUE);
    if(fSkewTermList) delete fSkewTermList;
    fSkewTermList = new TList();
    fSkewTermList->SetOwner(kTRUE);
    if(fKurtosisTermList) delete fKurtosisTermList;
    fKurtosisTermList = new TList();
    fKurtosisTermList->SetOwner(kTRUE);
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    if(fSubevent)
    {
        if(fSubList) delete fSubList;
        fSubList = new TList();
        fSubList->SetOwner(kTRUE);
    }
    for(int i=0;i<=2;++i)
    {
        fCkTermList->Add(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xlow,xhigh));
    }
    for(int i=0;i<=3;++i)
    {
        fSkewTermList->Add(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xlow,xhigh));
    }
    for(int i=0;i<=4;++i)     {
        fKurtosisTermList->Add(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xlow,xhigh));
    }
    for(int m=0;m<mpar;++m)
    {
        fCorrList->Add(new AliProfileBS(Form("corr_%ipar",m+1),this->GetTitle(),nbinsx,xlow,xhigh));
        if(fSubevent) {
            fSubList->Add(new AliProfileBS(Form("corr_%ipar_subP",m+1),this->GetTitle(),nbinsx,xlow,xhigh));
            fSubList->Add(new AliProfileBS(Form("corr_%ipar_subN",m+1),this->GetTitle(),nbinsx,xlow,xhigh));
            fSubList->Add(new AliProfileBS(Form("corr_%ipar_sub",m+1),this->GetTitle(),nbinsx,xlow,xhigh));
        }
    }
    printf("Container %s initialized with m = %i\n",this->GetName(),mpar);
};
void AliPtContainer::InitializeSubsamples(const int &nsub) {
    if(!fCkTermList || !fSkewTermList || !fCorrList) return;
    for(int i=0; i<fCkTermList->GetEntries();++i)
        ((AliProfileBS*)fCkTermList->At(i))->InitializeSubsamples(nsub);
    for(int i=0; i<fSkewTermList->GetEntries();++i)
        ((AliProfileBS*)fSkewTermList->At(i))->InitializeSubsamples(nsub);
    for(int i=0; i<fKurtosisTermList->GetEntries();++i)
        ((AliProfileBS*)fKurtosisTermList->At(i))->InitializeSubsamples(nsub);
    for(int i=0; i<fCorrList->GetEntries();++i)
        ((AliProfileBS*)fCorrList->At(i))->InitializeSubsamples(nsub);
    if(fSubevent) {
        if(!fSubList) return;
        for(int i=0; i<fSubList->GetEntries();++i)
            ((AliProfileBS*)fSubList->At(i))->InitializeSubsamples(nsub);
    }
    return;
};
vector<double> AliPtContainer::getEventCorrelation(int mOrder, Int_t subIndex) {
  vector<double> outvec = {fCorr[subIndex][mOrder],fSumw[subIndex][mOrder]};
  return outvec;
}
double AliPtContainer::getExoticEventCorrelation(int w, int p) {
  if(w>4) { printf("Exotic corr only supported up to m = 4\n"); return 0; }
  if(p>w) { printf("weight order must be >= pt order\n"); return 0; }
  return fExoticCorr[w][p];
}
void AliPtContainer::FillExotic(int mMax, const vector<vector<double>> &inarr) {
  vector<double> outer_terms;
  vector<double> inner_terms;
  fExoticCorr.resize(mMax+1,vector<double>(mMax+1,0.0)); fExoticCorr[0][0] = 1.0;
  double inner_sum;
  double outer_sum;
  for(int m(1);m<=mMax;++m){
    for(int n(0);n<=m;++n){
      for(int k(1);k<=m;++k){
        inner_terms.clear();
        for(int l(0);l<=m-k;++l){
          if(l>m-n) continue;
          if(m-k-l>n || m-k-l<0) continue;
          inner_terms.push_back(fCoeff[m][n][k][l]*fExoticCorr[m-k][m-k-l]*inarr[k][n-m+k+l]);
        }
        inner_sum = OrderedAddition(inner_terms);
        outer_terms.push_back(fSign[k-1]*(fFactorial[m-1]/fFactorial[m-k])*inner_sum);
      }
      outer_sum = OrderedAddition(outer_terms);
      outer_terms.clear();
      fExoticCorr[m][n] = outer_sum;
    }
  }
  return;
}
void AliPtContainer::FillRecursive(const vector<vector<double>> &inarr, Int_t subIndex) {
  if(subIndex<0||subIndex>2) return;
  fCorr[subIndex].clear(); fCorr[subIndex].resize(mpar+1,0); fCorr[subIndex][0] = 1.0;
  fSumw[subIndex].clear(); fSumw[subIndex].resize(mpar+1,0); fSumw[subIndex][0] = 1.0;
  double sumNum = 0;
  double sumDenum = 0;
  vector<double> valNum;
  vector<double> valDenum;
  for(int m(1); m<=mpar; ++m)
  {
    for(int k(1);k<=m;++k)
    {
      valNum.push_back(fSign[k-1]*fCorr[subIndex][m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][k]);
      valDenum.push_back(fSign[k-1]*fSumw[subIndex][m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][0]);
    }
    sumNum = OrderedAddition(valNum);
    sumDenum = OrderedAddition(valDenum);
    valNum.clear();
    valDenum.clear();
  
    fCorr[subIndex][m] = sumNum;
    fSumw[subIndex][m] = sumDenum;
  }
  return;
}
void AliPtContainer::FillRecursiveProfiles(const double &lMult, const double &rn, bool exotic) {
    int k = 0;
    for(int m=1;m<=mpar;++m)
    {
      if(m<5){
        for(int i(m);i<=4;++i){
          if(fSumw[0][m]!=0 && fSumw[0][i]!=0) ((AliProfileBS*)fCorrList->At(k))->FillProfile(lMult,(exotic)?fExoticCorr[i][m]/fSumw[0][i]:fCorr[0][m]/fSumw[0][m],(fEventWeight==PtSpace::kOne)?1.0:fSumw[0][i],rn);
          ++k;
        }
      }
      else if(fSumw[0][m]!=0) { ((AliProfileBS*)fCorrList->At(k))->FillProfile(lMult,fCorr[0][m]/fSumw[0][m],(fEventWeight==PtSpace::kOne)?1.0:fSumw[0][m],rn); ++k;}
      if(fSubevent){
        if(fSumw[1][m]!=0) ((AliProfileBS*)fSubList->At(3*(m-1)))->FillProfile(lMult,fCorr[1][m]/fSumw[1][m],(fEventWeight==PtSpace::kOne)?1.0:fSumw[1][m],rn);
        if(fSumw[2][m]!=0) ((AliProfileBS*)fSubList->At(3*(m-1)+1))->FillProfile(lMult,fCorr[2][m]/fSumw[2][m],(fEventWeight==PtSpace::kOne)?1.0:fSumw[2][m],rn);
        vector<double> vSub = getSubeventCorrelation(m);
        if(vSub[1]!=0) ((AliProfileBS*)fSubList->At(3*(m-1)+2))->FillProfile(lMult,vSub[0],(fEventWeight==PtSpace::kOne)?1.0:vSub[1],rn);
      }
    }
    return;
}
vector<Double_t> AliPtContainer::getSubeventCorrelation(Int_t mOrder){
  vector<double> outvec;
  double val = 0;
  double sumw = 0;
  switch(mOrder){
    case 1:
      val = (fCorr[1][1]/fSumw[1][1]+fCorr[2][1]/fSumw[2][1])/2;
      sumw = (fSumw[1][1]+fSumw[2][1])/2;
      outvec = {val,sumw};
      return outvec;
    case 2:
      val = fCorr[1][1]/fSumw[1][1]*fCorr[2][1]/fSumw[2][1];
      sumw = fSumw[1][1]*fSumw[2][1];
      outvec = {val,sumw};
      return outvec;
    case 3:
      val = (fCorr[1][2]/fSumw[1][2]*fCorr[2][1]/fSumw[2][1] + fCorr[2][2]/fSumw[2][2]*fCorr[1][1]/fSumw[1][1])/2;
      sumw = (fSumw[1][2]*fSumw[2][1] + fSumw[2][2]*fSumw[1][1])/2;
      outvec = {val,sumw};
      return outvec;
    case 4:
      val = fCorr[1][2]/fSumw[1][2]*fCorr[2][2]/fSumw[2][2];
      sumw = fSumw[1][2]*fSumw[2][2];
      outvec = {val,sumw};
      return outvec;
    case 5:
      val = (fCorr[1][2]/fSumw[1][2]*fCorr[2][3]/fSumw[2][3] + fCorr[1][3]/fSumw[1][3]*fCorr[2][2]/fSumw[2][2])/2;
      sumw = (fSumw[1][2]*fSumw[2][3] + fSumw[1][3]*fSumw[2][2])/2;
      outvec = {val,sumw};
      return outvec;
    case 6:
      val = fCorr[1][3]/fSumw[1][3]*fCorr[2][3]/fSumw[2][3];
      sumw = fSumw[1][3]*fSumw[2][3];
      outvec = {val,sumw};
      return outvec;
    case 7:
      val = (fCorr[1][4]/fSumw[1][4]*fCorr[2][3]/fSumw[2][3] + fCorr[1][3]/fSumw[1][3]*fCorr[2][4]/fSumw[2][4])/2;
      sumw = (fSumw[1][4]*fSumw[2][3] + fSumw[1][3]*fSumw[2][4])/2;
      outvec = {val,sumw};
      return outvec;
    case 8:
      val = fCorr[1][4]/fSumw[1][4]*fCorr[2][4]/fSumw[2][4];
      sumw = fSumw[1][4]*fSumw[2][4];
      outvec = {val,sumw};
      return outvec;
    default:
      outvec = {0.0,0.0};
      return outvec;
  }
  outvec = {0,0};
  return outvec;
}
double AliPtContainer::OrderedAddition(vector<double> vec) {
  double sum = 0;
  std::sort(vec.begin(), vec.end());

  for(int i = 0; i < vec.size(); i++)
  {
    sum += vec[i];
  }
  return sum;
}
void AliPtContainer::FillCk(const vector<vector<double>> &inarr, const double &lMult, const double &rn) {
    Double_t lwck = inarr[1][0]*inarr[1][0] - inarr[2][0];
    Double_t lwpt = inarr[1][0];
    if(lwck==0 || lwpt==0) return;
    ((AliProfileBS*)fCkTermList->At(0))->FillProfile(lMult,(inarr[1][1]*inarr[1][1] - inarr[2][2])/lwck,(fEventWeight==PtSpace::kOne)?1.0:lwck,rn);
    ((AliProfileBS*)fCkTermList->At(1))->FillProfile(lMult,(inarr[2][1]-inarr[1][1]*inarr[1][0])/lwck,(fEventWeight==PtSpace::kOne)?1.0:lwck,rn);
    ((AliProfileBS*)fCkTermList->At(2))->FillProfile(lMult,inarr[1][1]/lwpt,(fEventWeight==PtSpace::kOne)?1.0:lwpt,rn);
    return;
}
void AliPtContainer::FillSkew(const vector<vector<double>> &inarr, const double &lMult, const double &rn) {
    double lwpt = inarr[1][0];
    double lwskew = lwpt*lwpt*lwpt - 3*inarr[2][0]*inarr[1][0]+2*inarr[3][0];
    if(lwskew==0 || lwpt==0) return;
    /* Full expression
    <dptdptdpt> =  <(w1p1)^3-3*w2p2*wp + 2*w3p3>
    -3[pt]*   <(w1p1)^2*w1p0 - 2*w2p1*w1p1 + 2*w3p2 - w2p2*w1p0>
    +3[pt]^2* <w1p1*(w1p0)^2 - 2*w2p1*w1p0 + 2*w3p1 - w1p1*w2p0>
    -[pt]^3*  <(w1p0)^3-3*(w2p0)*(w1p0)+2*w3p0>   ----------> Cancelled by division of the weight
    */
    ((AliProfileBS*)fSkewTermList->At(0))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*inarr[1][1] - 3*inarr[2][2]*inarr[1][1]+2*inarr[3][3])/lwskew,(fEventWeight==PtSpace::kOne)?1.0:lwskew,rn);
    ((AliProfileBS*)fSkewTermList->At(1))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*inarr[1][0]-2*inarr[2][1]*inarr[1][1]+2*inarr[3][2]-inarr[2][2]*inarr[1][0])/lwskew,(fEventWeight==PtSpace::kOne)?1.0:lwskew,rn);
    ((AliProfileBS*)fSkewTermList->At(2))->FillProfile(lMult,(inarr[1][1]*inarr[1][0]*inarr[1][0]-2*inarr[2][1]*inarr[1][0]+2*inarr[3][1]-inarr[1][1]*inarr[2][0])/lwskew,(fEventWeight==PtSpace::kOne)?1.0:lwskew,rn);
    ((AliProfileBS*)fSkewTermList->At(3))->FillProfile(lMult,inarr[1][1]/lwpt,(fEventWeight==PtSpace::kOne)?1.0:lwpt,rn);
    return;
}
void AliPtContainer::FillKurtosis(const vector<vector<double>> &inarr, const double &lMult, const double &rn) {
    double lwpt = inarr[1][0];
    double lwkur = lwpt*lwpt*lwpt*lwpt-6*inarr[2][0]*lwpt*lwpt+8*lwpt*inarr[3][0]+3*inarr[2][0]*inarr[2][0]-6*inarr[4][0];
    if(lwkur==0 || lwpt==0) return;
    /* Full expression
    <dptdptdptdpt> = <(w1p1)^4 - 4*w3p3*w1p1 - 6*( w2p2(w1p1)^2-(w2p2)^2-2w3p3*w1p1 ) - 3(w2p2)^2 - 6*w4p4>
    -4*[pt]* <(w1p1)^3*w1p0 - w3p3*w1p0 - 3*w3p2*w1p1 - 3*( w2p2*w1p1*w1p0 - w3p3*w1p0 - w3p2*w1*p1 - w2p2*w2p1 ) - 3*( w2p1*(w1p1)^2-2*w3p2*w1p1-w2p1*w2p2 ) -3*w2p1*w2p2 - 6*w4p3>
    +6*[pt]^2* <(w1p1)^2*(w1p0)^2 - 2*w3p2*w1p0 - 2*w3p1*w1p1 - ( w2p2*(w1p0)^2 - w2p2*w2p0 - 2*w3p2*w1p0 ) - ( w2p0*(w1p1)^2 - w2p2*w2p0 -2*w3p1*w1p1 ) -4*( w2p1*w1p1*w1p0 - w3p2*w1p0 - w3p1*w1p1 - (w2p1)^2) - w2p2*w2p0 - 2*(w2p1)^2 - 6*w4p2>

    Simplified expression
    <dptdptdptdpt> = <(w1p1)^4 + 8*w3p3*w1p1 - 6*w2p2(w1p1)^2 + 3*(w2p2)^2 - 6*w4p4>  
    -4*[pt]* <(w1p1)^3*w1p0 + 2*w3p3*w1p0 + 6*w3p2*w1p1 - 3*w2p2*w1p1*w1p0 - 3*w2p1*(w1p1)^2 + 3*w2p2*w2p1 - 6*w4p3>
    +6*[pt]^2* <(w1p1)^2*(w1p0)^2 - w2p2*(w1p0)^2 - w2p0*(w1p1)^2 + w2p2*w2p0 - 4*w2p1*w1p1*w1p0 + 4*w3p2*w1p0 + 4*w3p1*w1p1 + 2*(w2p1)^2 - 6*w4p2>
    -4*[pt]^3* <w1p1*(w1p0)^3 - 3*w2p1*(w1p0)^2 - 3*w1p1*w2p0*w1p0 + 3*w2p1*w2p0 + 2*w1p1*w3p0 + 6*w3p1*w1p0 - 6*w4p1>
    +[pt]^4*     <(w1p0)^4 + 8*w3p0*w1p0 - 6*w2p0*(w1p0)^2 + 3*(w2p0)^2 - 6*w4p0>
    */
    ((AliProfileBS*)fKurtosisTermList->At(0))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*inarr[1][1]*inarr[1][1]-6*inarr[2][2]*inarr[1][1]*inarr[1][1]
                                                    +3*inarr[2][2]*inarr[2][2]+8*inarr[1][1]*inarr[3][3]-6*inarr[4][4])/lwkur,(fEventWeight==PtSpace::kOne)?1.0:lwkur,rn);
    ((AliProfileBS*)fKurtosisTermList->At(1))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*inarr[1][1]*lwpt-3*inarr[2][2]*inarr[1][1]*lwpt
                                                    -3*inarr[1][1]*inarr[1][1]*inarr[2][1]+3*inarr[2][2]*inarr[2][1]
                                                    +6*inarr[1][1]*inarr[3][2]+2*inarr[3][3]*lwpt-6*inarr[4][3])/lwkur,(fEventWeight==PtSpace::kOne)?1.0:lwkur,rn);
    ((AliProfileBS*)fKurtosisTermList->At(2))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*lwpt*lwpt-inarr[2][2]*lwpt*lwpt-inarr[2][0]*inarr[1][1]*inarr[1][1] 
                                                    +inarr[2][0]*inarr[2][2]-4*inarr[2][1]*inarr[1][1]*lwpt+4*inarr[3][2]*lwpt 
                                                    +4*inarr[3][1]*inarr[1][1]+2*inarr[2][1]*inarr[2][1]-6*inarr[4][2])/lwkur,(fEventWeight==PtSpace::kOne)?1.0:lwkur,rn);
    ((AliProfileBS*)fKurtosisTermList->At(3))->FillProfile(lMult,(inarr[1][1]*lwpt*lwpt*lwpt-3*inarr[2][1]*lwpt*lwpt
                                                    -3*inarr[1][1]*inarr[2][0]*lwpt+3*inarr[2][1]*inarr[2][0]+2*inarr[1][1]*inarr[3][0]
                                                    +6*inarr[3][1]*lwpt-6*inarr[4][1])/lwkur,(fEventWeight==PtSpace::kOne)?1.0:lwkur,rn);
    ((AliProfileBS*)fKurtosisTermList->At(4))->FillProfile(lMult,inarr[1][1]/lwpt,(fEventWeight==PtSpace::kOne)?1.0:lwpt,rn);
    return;
}
void AliPtContainer::RebinMulti(Int_t nbins) {   
    if(fCkTermList) for(Int_t i=0;i<fCkTermList->GetEntries();i++) ((AliProfileBS*)fCkTermList->At(i))->RebinMulti(nbins); 
    if(fSkewTermList) for(Int_t i=0;i<fSkewTermList->GetEntries();i++) ((AliProfileBS*)fSkewTermList->At(i))->RebinMulti(nbins); 
    if(fKurtosisTermList) for(Int_t i=0;i<fKurtosisTermList->GetEntries();i++) ((AliProfileBS*)fKurtosisTermList->At(i))->RebinMulti(nbins); 
    if(fCorrList) for(Int_t i=0;i<fCorrList->GetEntries();i++) ((AliProfileBS*)fCorrList->At(i))->RebinMulti(nbins); 
    if(fSubList) for(Int_t i=0;i<fSubList->GetEntries();i++) ((AliProfileBS*)fSubList->At(i))->RebinMulti(nbins); 
    return;
}
void AliPtContainer::RebinMulti(Int_t nbins, Double_t *binedges) {
    if(fCkTermList) for(Int_t i=0;i<fCkTermList->GetEntries();i++) ((AliProfileBS*)fCkTermList->At(i))->RebinMulti(nbins,binedges); 
    if(fSkewTermList) for(Int_t i=0;i<fSkewTermList->GetEntries();i++) ((AliProfileBS*)fSkewTermList->At(i))->RebinMulti(nbins,binedges); 
    if(fKurtosisTermList) for(Int_t i=0;i<fKurtosisTermList->GetEntries();i++) ((AliProfileBS*)fKurtosisTermList->At(i))->RebinMulti(nbins,binedges); 
    if(fCorrList) for(Int_t i=0;i<fCorrList->GetEntries();i++) ((AliProfileBS*)fCorrList->At(i))->RebinMulti(nbins,binedges); 
    if(fSubList) for(Int_t i=0;i<fSubList->GetEntries();i++) ((AliProfileBS*)fSubList->At(i))->RebinMulti(nbins,binedges); 
    return;
}
TH1* AliPtContainer::getCkHist(int ind) {
  if(!fCkList) CalculateCk();
  if(!fCkList) return 0;
  if(ind+1<fCkList->GetEntries()) return (TH1*)fCkList->At(ind+1);
  return 0;
}
TH1* AliPtContainer::getSkewHist(int ind) {
  if(!fSkewList) CalculateSkew();
  if(!fSkewList) return 0;
  if(ind+1<fSkewList->GetEntries()) return (TH1*)fSkewList->At(ind+1);
  return 0;
}
TH1* AliPtContainer::getKurtosisHist(int ind) {
  if(!fKurtosisList) CalculateKurtosis();
  if(!fKurtosisList) return 0;
  if(ind+1<fKurtosisList->GetEntries()) return (TH1*)fKurtosisList->At(ind+1);
  return 0;
}
void AliPtContainer::CalculateCk() {
  if(fCkList) delete fCkList;
  fCkList = new TList();
  fCkList->SetOwner(kTRUE);
  //Override mpt temporarily override mpt weights:
  ((AliProfileBS*)fCkTermList->At(2))->PresetWeights((AliProfileBS*)fCkTermList->At(0));
  for(Int_t i=-1;i<((AliProfileBS*)fCkTermList->At(0))->getNSubs();i++) {
    vector<TH1*> hTs;
    for(Int_t j=0;j<3;j++) {
      ((AliProfileBS*)fCkTermList->At(j))->SetErrorOption("g");
      hTs.push_back(((AliProfileBS*)fCkTermList->At(j))->getHist(i));
    }
    TH1 *hRec = RecalculateCkHists(hTs);
    hRec->SetName(Form("%s_Recalculated%i",this->GetName(),i+1));
    fCkList->Add(hRec);
    hTs.clear();
  };
  //reset mpt weights
  ((AliProfileBS*)fCkTermList->At(2))->PresetWeights(0);
}
void AliPtContainer::CalculateSkew() {
  if(fSkewList) delete fSkewList;
  fSkewList = new TList();
  fSkewList->SetOwner(kTRUE);
  //Override mpt temporarily override mpt weights:
  ((AliProfileBS*)fSkewTermList->At(3))->PresetWeights((AliProfileBS*)fSkewTermList->At(0));
  for(Int_t i=-1;i<((AliProfileBS*)fSkewTermList->At(0))->getNSubs();i++) {
    vector<TH1*> hTs;
    for(Int_t j=0;j<4;j++) {
      ((AliProfileBS*)fSkewTermList->At(j))->SetErrorOption("g");
      hTs.push_back(((AliProfileBS*)fSkewTermList->At(j))->getHist(i));
    }
    TH1 *hRec = RecalculateSkewHists(hTs);
    hRec->SetName(Form("%s_Recalculated%i",this->GetName(),i+1));
    fSkewList->Add(hRec);
    hTs.clear();
  };
  //reset mpt weights
  ((AliProfileBS*)fSkewTermList->At(3))->PresetWeights(0);
}
void AliPtContainer::CalculateKurtosis() {
  if(fKurtosisList) delete fKurtosisList;
  fKurtosisList = new TList();
  fKurtosisList->SetOwner(kTRUE);
  //Override mpt temporarily override mpt weights:
  ((AliProfileBS*)fKurtosisTermList->At(4))->PresetWeights((AliProfileBS*)fKurtosisTermList->At(0));
  for(Int_t i=-1;i<((AliProfileBS*)fKurtosisTermList->At(0))->getNSubs();i++) {
    vector<TH1*> hTs;
    for(Int_t j=0;j<5;j++) {
      ((AliProfileBS*)fKurtosisTermList->At(j))->SetErrorOption("g");
      hTs.push_back(((AliProfileBS*)fKurtosisTermList->At(j))->getHist(i));
    }
    TH1 *hRec = RecalculateKurtosisHists(hTs);
    hRec->SetName(Form("%s_Recalculated%i",this->GetName(),i+1));
    fKurtosisList->Add(hRec);
    hTs.clear();
  };
  //reset mpt weights
  ((AliProfileBS*)fKurtosisTermList->At(4))->PresetWeights(0);
}
TH1 *AliPtContainer::RecalculateCkHists(vector<TH1*> inh) {
  inh[1]->Multiply(inh[2]);
  inh[1]->Scale(2);
  TH1 *mptsq = (TH1*)inh[2]->Clone("mptSquared");
  mptsq->Multiply(inh[2]);
  TH1 *hWeights = (TH1*)inh[0]->Clone("ForErrors");
  inh[0]->Add(inh[1]);
  inh[0]->Add(mptsq);
  delete mptsq;
  for(Int_t i=1;i<=hWeights->GetNbinsX();i++) inh[0]->SetBinError(i,hWeights->GetBinError(i));
  delete hWeights;
  return (TH1*)inh[0]->Clone("hRec");
}
TH1 *AliPtContainer::RecalculateSkewHists(vector<TH1*> inh) {
  inh[1]->Multiply(inh[3]);
  inh[1]->Scale(3);
  TH1 *mptsq = (TH1*)inh[3]->Clone("mptSquared");
  TH1 *mptcb = (TH1*)inh[3]->Clone("mptCubed");
  mptsq->Multiply(inh[3]);
  mptcb->Multiply(mptsq);
  inh[2]->Multiply(mptsq);
  inh[2]->Scale(3);
  TH1 *hWeights = (TH1*)inh[0]->Clone("ForErrors");
  inh[0]->Add(inh[1],-1);
  inh[0]->Add(inh[2]);
  inh[0]->Add(mptcb,-1);
  delete mptsq;
  delete mptcb;
  for(Int_t i=1;i<=hWeights->GetNbinsX();i++) inh[0]->SetBinError(i,hWeights->GetBinError(i));
  delete hWeights;
  return (TH1*)inh[0]->Clone("hRec");
}
TH1 *AliPtContainer::RecalculateKurtosisHists(vector<TH1*> inh) {
  inh[1]->Multiply(inh[4]);
  inh[1]->Scale(4);
  TH1 *mptsq = (TH1*)inh[4]->Clone("mptSquared");
  mptsq->Multiply(inh[4]);
  TH1 *mptcb = (TH1*)mptsq->Clone("mptCubed");
  mptcb->Multiply(inh[4]);
  TH1 *mpt4th = (TH1*)mptcb->Clone("mpt4th");
  mpt4th->Multiply(inh[4]);
  inh[2]->Multiply(mptsq);
  inh[2]->Scale(6);
  inh[3]->Multiply(mptcb);
  inh[3]->Scale(4);
  TH1 *hWeights = (TH1*)inh[0]->Clone("ForErrors");
  inh[0]->Add(inh[1],-1);
  inh[0]->Add(inh[2]);
  inh[0]->Add(inh[3],-1);
  inh[0]->Add(mpt4th);
  delete mptsq;
  delete mptcb;
  delete mpt4th;
  for(Int_t i=1;i<=hWeights->GetNbinsX();i++) inh[0]->SetBinError(i,hWeights->GetBinError(i));
  delete hWeights;
  return (TH1*)inh[0]->Clone("hRec");
}
TH1* AliPtContainer::getRecursiveHist(int ind, int m, unsigned int l_obs, bool sub) {
  if(l_obs==kObs::kCorr) {
      if(!sub) return (m<5)?((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar_%ipc",m,m)))->getHist(ind):((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar",m)))->getHist(ind);
      else return getSubeventCumulantHist(ind, m);
  }
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
TH1* AliPtContainer::getSubeventCumulantHist(int ind, int m) {
    TH1* reth;
    if(m==2) {
        reth = ((AliProfileBS*)fSubList->At(2))->getHist(ind);
        TH1* sub1N = ((AliProfileBS*)fSubList->At(1))->getHist(ind);
        TH1* sub1P = ((AliProfileBS*)fSubList->At(0))->getHist(ind);
        sub1N->Multiply(sub1P);
        reth->Add(sub1N,-1);
    }
    else if(m==3) {
        reth = ((AliProfileBS*)fSubList->At(2))->getHist(ind);
        TH1* sub2N = ((AliProfileBS*)fSubList->At(3))->getHist(ind);
        TH1* sub1P = ((AliProfileBS*)fSubList->At(0))->getHist(ind);
        TH1* sub1N = ((AliProfileBS*)fSubList->At(1))->getHist(ind);
        reth->Multiply(sub1N);
        sub2N->Multiply(sub1P);
        reth->Add(sub2N);
        reth->Scale(0.5);
    }
    else if(m==4) {
        reth = ((AliProfileBS*)fSubList->At(2))->getHist(ind);
        TH1* sub2N = ((AliProfileBS*)fSubList->At(3))->getHist(ind);
        reth->Multiply(sub2N);
    }
    else { printf("Subevent only for 2-,3- and 4-particle\n"); return 0; }
    return reth;
}
void AliPtContainer::CalculateRecursive(bool normalized) { 
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
    ((AliProfileBS*)fCorrList->At(0))->PresetWeights((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar%s",mpar,(mpar>4)?"":Form("_%ipc",mpar))));
    for(int i=-1;i<((AliProfileBS*)fCorrList->At(0))->getNSubs();++i) {
        vector<TH1*> hTs;
        for(int j=1;j<=mpar;++j) {
            ((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar%s",j,(j>4)?"":Form("_%ipc",j))))->SetErrorOption("g");
            hTs.push_back(((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar%s",j,(j>4)?"":Form("_%ipc",j))))->getHist(i));
        }
        CalculateCumulantHists(hTs,i,normalized);
    }
    ((AliProfileBS*)fCorrList->At(0))->PresetWeights(0);
    return;
}
void AliPtContainer::CalculateCumulantHists(vector<TH1*> inh, int ind, bool normalized) {
    auto binomial = [&](const int n, const int m) { return factorial(n)/(factorial(m)*factorial(n-m)); };
    int lMax = (((AliProfileBS*)fCorrList->At(0))->getNSubs()+1)*mpar;
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
TH1* AliPtContainer::getPowerHist(TH1* inh, double p) {
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
Long64_t AliPtContainer::Merge(TCollection *collist) {
  Long64_t nmerged=0;
  TIter all_PTC(collist);
  AliPtContainer *l_PTC = 0;
  while ((l_PTC = ((AliPtContainer*) all_PTC()))) {
      TList *t_CkTerm = l_PTC->fCkTermList;
      TList *t_SkewTerm = l_PTC->fSkewTermList;
      TList *t_KurtosisTerm = l_PTC->fKurtosisTermList;
      TList *t_Ck = l_PTC->fCkList;
      TList *t_Skew = l_PTC->fSkewList;
      TList *t_Kurtosis = l_PTC->fKurtosisList;
      TList* t_Corr = l_PTC->fCorrList;
      TList* t_Sub = l_PTC->fCorrList;
      TList* t_Cum = l_PTC->fCumulantList;
      TList* t_Norm = l_PTC->fNormList;
      if(t_CkTerm) {
        if(!fCkTermList) fCkTermList = (TList*)t_CkTerm->Clone();
        else MergeBSLists(fCkTermList,t_CkTerm);
        nmerged++;
      };
      if(t_SkewTerm) {
        if(!fSkewTermList) fSkewTermList = (TList*)t_SkewTerm->Clone();
        else MergeBSLists(fSkewTermList,t_SkewTerm);
        nmerged++;
      };
      if(t_KurtosisTerm) {
        if(!fKurtosisTermList) fKurtosisTermList = (TList*)t_KurtosisTerm->Clone();
        else MergeBSLists(fKurtosisTermList,t_KurtosisTerm);
        nmerged++;
      };
      if(t_Ck) {
        if(!fCkList) fCkList = (TList*)t_Ck->Clone();
        else MergeBSLists(fCkList,t_Ck);
        nmerged++;
      };
      if(t_Skew) {
        if(!fSkewList) fSkewList = (TList*)t_Skew->Clone();
        else MergeBSLists(fSkewList,t_Skew);
        nmerged++;
      };
      if(t_Kurtosis) {
        if(!fKurtosisList) fKurtosisList = (TList*)t_Kurtosis->Clone();
        else MergeBSLists(fKurtosisList,t_Kurtosis);
        nmerged++;
      };
      if(t_Corr) {
          if(!fCorrList) fCorrList = (TList*)t_Corr->Clone();
          else MergeBSLists(fCorrList,t_Corr);
      };
      if(t_Sub) {
          if(!fSubList) fSubList = (TList*)t_Sub->Clone();
          else MergeBSLists(fSubList,t_Sub);
      };
      if(t_Cum) {
          if(!fCumulantList) fCumulantList = (TList*)t_Cum->Clone();
          else MergeBSLists(fCumulantList,t_Cum);
      };
      if(t_Norm) {
          if(!fNormList) fNormList = (TList*)t_Norm->Clone();
          else MergeBSLists(fNormList,t_Norm);
      }
  }
  return nmerged;
}
void AliPtContainer::MergeBSLists(TList *source, TList *target) {
  if(source->GetEntries()!=target->GetEntries()) { printf("Number in lists to be merged are not the same, skipping...\n"); return; };
  for(Int_t i=0;i<source->GetEntries();i++) {
    AliProfileBS *l_obj = (AliProfileBS*)source->At(i);
    AliProfileBS *t_obj = (AliProfileBS*)target->At(i);
    l_obj->MergeBS(t_obj);
  };
}