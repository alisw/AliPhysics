/*
Author: Emil Gorm Nielsen
Modified from AliCkContainer by Vytautas
Container to store the terms required to calculate higher order moments of the pT spectrum.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#include "AliPtPtContainer.h"

Double_t AliPtPtContainer::fFactorial[9] = {1.,1.,2.,6.,24.,120.,720.,5040.,40320.};
Int_t AliPtPtContainer::fSign[9] = {1,-1,1,-1,1,-1,1,-1,1};

using namespace PtPtSpace;
AliPtPtContainer::AliPtPtContainer():
    fCMTermList(0),
    fCorrList(0),
    fCentralMomentList(0),
    fCumulantList(0),
    mpar(0),
    fEventWeight(kEventWeight::kUnity),
    fCorr(),
    fSumw(),
    fcmNum(),
    fcmDen()
{};
AliPtPtContainer::~AliPtPtContainer()
{
    delete fCMTermList;
    delete fCorrList;
};
AliPtPtContainer::AliPtPtContainer(const char* name, const char* title, int nbinsx, double* xbins, int m):
    TNamed(name,title),
    fCMTermList(0),
    fCorrList(0),
    fCentralMomentList(0),
    fCumulantList(0),
    mpar(m),
    fEventWeight(kEventWeight::kUnity),
    fCorr(),
    fSumw(),
    fcmNum(),
    fcmDen()
{
    Initialize(nbinsx,xbins);
};
AliPtPtContainer::AliPtPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, int m):
    TNamed(name,title),
    fCMTermList(0),
    fCorrList(0),
    fCentralMomentList(0),
    fCumulantList(0),
    mpar(m),
    fEventWeight(kEventWeight::kUnity),
    fCorr(),
    fSumw(),
    fcmNum(),
    fcmDen()
{
    Initialize(nbinsx,xlow,xhigh);
};
void AliPtPtContainer::Initialize(int nbinsx, double* xbins) {
    if(fCMTermList) delete fCMTermList;
    fCMTermList = new TList();
    fCMTermList->SetOwner(kTRUE);
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    for(int m=0;m<mpar;++m) fCorrList->Add(new AliProfileBS(Form("corr_%ipar",m+1),this->GetTitle(),nbinsx,xbins));
    for(int m=0;m<4;++m) {
        for(int i=0;i<=m;++i)
            fCMTermList->Add(new AliProfileBS(Form("cm%i_Mpt%i",m+1,i),this->GetTitle(),nbinsx,xbins));
    }
    printf("Container %s initialized with m = %i\n",this->GetName(),mpar);
};
void AliPtPtContainer::Initialize(int nbinsx, double xlow, double xhigh) {
    if(fCMTermList) delete fCMTermList;
    fCMTermList = new TList();
    fCMTermList->SetOwner(kTRUE);
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    for(int m=0;m<mpar;++m) fCorrList->Add(new AliProfileBS(Form("corr_%ipar",m+1),this->GetTitle(),nbinsx,xlow,xhigh));
    for(int m=0;m<4;++m) {
        for(int i=0;i<=m;++i)
            fCMTermList->Add(new AliProfileBS(Form("cm%i_Mpt%i",m+1,i),this->GetTitle(),nbinsx,xlow,xhigh));
    }
    printf("Container %s initialized with m = %i\n",this->GetName(),mpar);
};
void AliPtPtContainer::InitializeSubsamples(const int &nsub) {
    if(!fCMTermList || !fCorrList) return;
    for(int i=0; i<fCMTermList->GetEntries();++i)
            ((AliProfileBS*)fCMTermList->At(i))->InitializeSubsamples(nsub);
    for(int i=0; i<fCorrList->GetEntries();++i)
        ((AliProfileBS*)fCorrList->At(i))->InitializeSubsamples(nsub);
    return;
};
vector<double> AliPtPtContainer::getEventCorrelation(int mOrder) {
  vector<double> outvec = {fCorr[mOrder],fSumw[mOrder]};
  return outvec;
}
void AliPtPtContainer::CalculateCorrelations(const vector<vector<double>> &inarr) {
  fCorr.clear(); fCorr.resize(mpar+1,0); fCorr[0] = 1.0;
  fSumw.clear(); fSumw.resize(mpar+1,0); fSumw[0] = 1.0;
  double sumNum = 0;
  double sumDenum = 0;
  vector<double> valNum;
  vector<double> valDenum;
  for(int m(1); m<=mpar; ++m)
  {
    for(int k(1);k<=m;++k)
    {
      valNum.push_back(fSign[k-1]*fCorr[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][k]);
      valDenum.push_back(fSign[k-1]*fSumw[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][0]);
    }
    sumNum = OrderedAddition(valNum);
    sumDenum = OrderedAddition(valDenum);
    valNum.clear();
    valDenum.clear();

    fCorr[m] = sumNum;
    fSumw[m] = sumDenum;
  }
  return;
}
void AliPtPtContainer::FillProfiles(const double &centmult, const double &rn) {
    for(int m=1;m<=mpar;++m)
    {
        if(fSumw[m]!=0) ((AliProfileBS*)fCorrList->At(m-1))->FillProfile(centmult,fCorr[m]/fSumw[m],(fEventWeight==PtPtSpace::kUnity)?1.0:fSumw[m],rn);
    }
    return;
}
void AliPtPtContainer::FillCMProfiles(const vector<vector<double>> &inarr, const double &centmult, const double &rn) {
  fcmNum.clear();
  fcmDen.clear();
  if (inarr[0][0] == 0)
    return;
  double tau1 = inarr[2][0] / pow(inarr[1][0], 2);
  double tau2 = inarr[3][0] / pow(inarr[1][0], 3);
  double tau3 = inarr[4][0] / pow(inarr[1][0], 4);
  fcmDen.push_back(inarr[1][0]);
  fcmDen.push_back(1 - tau1);
  fcmDen.push_back(1 - 3 * tau1 + 2 * tau2);
  fcmDen.push_back(1 - 6 * tau1 + 3 * tau1 * tau1 + 8 * tau2 - 6 * tau3);
  // double weight4 = 1 - 10*tau1 + 15*tau1*tau1 + 20*tau2 - 20*tau1*tau2 - 30*tau3 + 24*tau4;
  if (mpar < 1 || fcmDen[0] == 0)
    return;
  fcmNum.push_back(inarr[1][1] / fcmDen[0]);
  dynamic_cast<AliProfileBS*>(fCMTermList->At(0))->FillProfile(centmult, fcmNum[0], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[0], rn);
  if (mpar < 2 || inarr[2][0] == 0 || fcmDen[1] == 0)
    return;
  fcmNum.push_back(1 / fcmDen[1] * (inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] - tau1 * inarr[2][2] / inarr[2][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(1))->FillProfile(centmult, fcmNum[1], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[1], rn);
  fcmNum.push_back(1 / fcmDen[1] * (-2 * inarr[1][1] / inarr[1][0] + 2 * tau1 * inarr[2][1] / inarr[2][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(2))->FillProfile(centmult, fcmNum[2], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[1], rn);
  if (mpar < 3 || inarr[3][0] == 0 || fcmDen[2] == 0)
    return;
  fcmNum.push_back(1 / fcmDen[2] * (inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] - 3 * tau1 * inarr[2][2] / inarr[2][0] * inarr[1][1] / inarr[1][0] + 2 * tau2 * inarr[3][3] / inarr[3][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(3))->FillProfile(centmult, fcmNum[3], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[2], rn);
  fcmNum.push_back(1 / fcmDen[2] * (-3 * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] + 3 * tau1 * inarr[2][2] / inarr[2][0] + 6 * tau1 * inarr[2][1] / inarr[2][0] * inarr[1][1] / inarr[1][0] - 6 * tau2 * inarr[3][2] / inarr[3][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(4))->FillProfile(centmult, fcmNum[4], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[2], rn);
  fcmNum.push_back(1 / fcmDen[2] * (3 * inarr[1][1] / inarr[1][0] - 6 * tau1 * inarr[2][1] / inarr[2][0] - 3 * tau1 * inarr[1][1] / inarr[1][0] + 6 * tau2 * inarr[3][1] / inarr[3][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(5))->FillProfile(centmult, fcmNum[5], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[2], rn);
  if (mpar < 4 || inarr[4][0] == 0 || fcmDen[3] == 0)
    return;
  fcmNum.push_back(1 / fcmDen[3] * (inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] - 6 * tau1 * inarr[2][2] / inarr[2][0] * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] + 3 * tau1 * tau1 * inarr[2][2] / inarr[2][0] * inarr[2][2] / inarr[2][0] + 8 * tau2 * inarr[3][3] / inarr[3][0] * inarr[1][1] / inarr[1][0] - 6 * tau3 * inarr[4][4] / inarr[4][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(6))->FillProfile(centmult, fcmNum[6], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[3], rn);
  fcmNum.push_back(1 / fcmDen[3] * (-4 * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] + 12 * tau1 * inarr[2][2] / inarr[2][0] * inarr[1][1] / inarr[1][0] + 12 * tau1 * inarr[2][1] / inarr[2][0] * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] - 12 * tau1 * tau1 * inarr[2][2] / inarr[2][0] * inarr[2][1] / inarr[2][0] - 8 * tau2 * inarr[3][3] / inarr[3][0] - 24 * tau2 * inarr[3][2] / inarr[3][0] * inarr[1][1] / inarr[1][0] + 24 * tau3 * inarr[4][3] / inarr[4][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(7))->FillProfile(centmult, fcmNum[7], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[3], rn);
  fcmNum.push_back(1 / fcmDen[3] * (6 * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] - 6 * tau1 * inarr[2][2] / inarr[2][0] - 24 * tau1 * inarr[2][1] / inarr[2][0] * inarr[1][1] / inarr[1][0] - 6 * tau1 * inarr[1][1] / inarr[1][0] * inarr[1][1] / inarr[1][0] + 6 * tau1 * tau1 * inarr[2][2] / inarr[2][0] + 12 * tau1 * tau1 * inarr[2][1] / inarr[2][0] * inarr[2][1] / inarr[2][0] + 24 * tau2 * inarr[3][2] / inarr[3][0] + 24 * tau2 * inarr[3][1] / inarr[3][0] * inarr[1][1] / inarr[1][0] - 36 * tau3 * inarr[4][2] / inarr[4][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(8))->FillProfile(centmult, fcmNum[8], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[3], rn);
  fcmNum.push_back(1 / fcmDen[3] * (-4 * inarr[1][1] / inarr[1][0] + 12 * tau1 * inarr[2][1] / inarr[2][0] + 12 * tau1 * inarr[1][1] / inarr[1][0] - 12 * tau1 * tau1 * inarr[2][1] / inarr[2][0] - 24 * tau2 * inarr[3][1] / inarr[3][0] - 8 * tau2 * inarr[1][1] / inarr[1][0] + 24 * tau3 * inarr[4][1] / inarr[4][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(9))->FillProfile(centmult, fcmNum[9], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[3], rn);
  return;
}
double AliPtPtContainer::OrderedAddition(vector<double> vec) {
  double sum = 0;
  std::sort(vec.begin(), vec.end());

  for(int i = 0; i < vec.size(); i++)
  {
    sum += vec[i];
  }
  return sum;
}
void AliPtPtContainer::RebinMulti(Int_t nbins) {
    if(fCMTermList) for(Int_t i=0;i<fCMTermList->GetEntries();i++) ((AliProfileBS*)fCMTermList->At(i))->RebinMulti(nbins);
    if(fCorrList) for(Int_t i=0;i<fCorrList->GetEntries();i++) ((AliProfileBS*)fCorrList->At(i))->RebinMulti(nbins);
    return;
}
void AliPtPtContainer::RebinMulti(Int_t nbins, Double_t *binedges) {
    if(fCMTermList) for(Int_t i=0;i<fCMTermList->GetEntries();i++) ((AliProfileBS*)fCMTermList->At(i))->RebinMulti(nbins,binedges);
    if(fCorrList) for(Int_t i=0;i<fCorrList->GetEntries();i++) ((AliProfileBS*)fCorrList->At(i))->RebinMulti(nbins,binedges);
    return;
}
TH1* AliPtPtContainer::getCorrHist(int ind, int m) {
    return ((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar",m)))->getHist(ind);
}
TH1* AliPtPtContainer::getCentralMomentHist(int ind, int m) {
  if(!fCentralMomentList) CreateCentralMomentList();
  if(!fCentralMomentList) return 0;
  if(ind+1<fCentralMomentList->GetEntries()) return (TH1*)fCentralMomentList->FindObject(Form("cm%i_%i",m,ind));
  return 0;
}
void AliPtPtContainer::CreateCentralMomentList() {
    if(fCentralMomentList) delete fCentralMomentList;
    fCentralMomentList = new TList();
    fCentralMomentList->SetOwner();
    for(auto m(1); m<=4;++m){
        for(int i=-1;i<((AliProfileBS*)fCMTermList->At(0))->getNSubs();++i) {
            TH1* hMpt = ((AliProfileBS*)fCMTermList->At(0))->getHist(i);
            vector<TH1*> hTs;
            for(int j=0;j<m;++j) {
                ((AliProfileBS*)fCMTermList->FindObject(Form("cm%i_Mpt%i",m,j)))->SetErrorOption("g");
                hTs.push_back(((AliProfileBS*)fCMTermList->FindObject(Form("cm%i_Mpt%i",m,j)))->getHist(i));
            }
            CalculateCentralMomentHists(hTs,i,m,hMpt);
        }
    }
    return;
}
void AliPtPtContainer::CalculateCentralMomentHists(vector<TH1*> inh, int ind, int m, TH1* hMpt) {
    TH1* reth = (TH1*)inh[0]->Clone(Form("cm%i_%i",m,ind));
    for(auto i(1);i<m;++i){
        TH1* mptPow = raiseHistToPower(hMpt,i);
        inh[i]->Multiply(mptPow);
        reth->Add(inh[i]);
    }
    TH1* mptLast = raiseHistToPower(hMpt,m);
    reth->Add(mptLast,(m%2)?(-1):1);
    fCentralMomentList->Add(reth);
    return;
}
TH1* AliPtPtContainer::getCumulantHist(int ind, int m) {
    if(!fCumulantList) CreateCumulantList();
    if(!fCumulantList) return 0;
    if(ind+1<fCumulantList->GetEntries()) return (TH1*)fCumulantList->At((ind+1)*mpar+m-1);
}
void AliPtPtContainer::CreateCumulantList() {
    if(fCumulantList) delete fCumulantList;
    fCumulantList = new TList();
    fCumulantList->SetOwner();
    //((AliProfileBS*)fCorrList->At(0))->PresetWeights((AliProfileBS*)fCorrList->At(mpar-1));
    for(int i=-1;i<((AliProfileBS*)fCorrList->At(0))->getNSubs();++i) {
        vector<TH1*> hTs;
        for(int j=0;j<mpar;++j) {
            ((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar",j+1)))->SetErrorOption("g");
            hTs.push_back(((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar",j+1)))->getHist(i));
        }
        CalculateCumulantHists(hTs,i);
    }
    //((AliProfileBS*)fCorrList->At(0))->PresetWeights(0);
    return;
}
void AliPtPtContainer::CalculateCumulantHists(vector<TH1*> inh, int ind) {
    auto binomial = [&](const int n, const int m) { return factorial(n)/(factorial(m)*factorial(n-m)); };
    for(int m=1;m<=mpar;++m)
    {
        TH1* reth = (TH1*)inh[m-1]->Clone(Form("reth%i_%i",m,ind));
        //TH1* hWeights = (TH1*)inh[m-1]->Clone(Form("hWeights%i_%i",m,ind));
        for(int k=1;k<m;++k)
        {
            TH1* corrh = (TH1*)inh[m-k-1]->Clone(Form("hcorr%i%i_%i",m,k,ind));
            corrh->Multiply((TH1*)fCumulantList->At((ind+1)*mpar+k-1));
            corrh->Scale(binomial(m-1,k-1));
            reth->Add(corrh,-1);
            delete corrh;
        }
        //for(int i=1;i<=hWeights->GetNbinsX();++i) reth->SetBinError(i,hWeights->GetBinError(i));
        //delete hWeights;
        fCumulantList->Add(reth->Clone(Form("kappa%i_%i",m,ind)));
    }
    return;
}
Long64_t AliPtPtContainer::Merge(TCollection *collist) {
  Long64_t nmerged=0;
  TIter all_PTC(collist);
  AliPtPtContainer *l_PTC = 0;
  while ((l_PTC = ((AliPtPtContainer*) all_PTC()))) {
      TList *t_CMTerm = l_PTC->fCMTermList;
      TList* t_Corr = l_PTC->fCorrList;
      TList* t_Cum = l_PTC->fCumulantList;
      TList* t_CM = l_PTC->fCentralMomentList;
      if(t_CMTerm) {
        if(!fCMTermList) fCMTermList = (TList*)t_CMTerm->Clone();
        else MergeBSLists(fCMTermList,t_CMTerm);
        nmerged++;
      };
      if(t_Corr) {
          if(!fCorrList) fCorrList = (TList*)t_Corr->Clone();
          else MergeBSLists(fCorrList,t_Corr);
      };
      if(t_Cum) {
          if(!fCumulantList) fCumulantList = (TList*)t_Cum->Clone();
          else MergeBSLists(fCumulantList,t_Cum);
      };
      if(t_CM) {
          if(!fCentralMomentList) fCentralMomentList = (TList*)t_CM->Clone();
          else MergeBSLists(fCentralMomentList,t_CM);
      }
  }
  return nmerged;
}
void AliPtPtContainer::MergeBSLists(TList *source, TList *target) {
  if(source->GetEntries()!=target->GetEntries()) { printf("Number in lists to be merged are not the same, skipping...\n"); return; };
  for(Int_t i=0;i<source->GetEntries();i++) {
    AliProfileBS *l_obj = (AliProfileBS*)source->At(i);
    AliProfileBS *t_obj = (AliProfileBS*)target->At(i);
    l_obj->MergeBS(t_obj);
  };
}
TH1* AliPtPtContainer::raiseHistToPower(TH1* inh, double p)
{
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
