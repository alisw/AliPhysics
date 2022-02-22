/*
Author: Emil Gorm Nielsen
Modified from AliCkContainer by Vytautas
Container to store the terms required to calculate higher order moments of the pT spectrum.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#include "AliPtContainer.h"

double AliPtContainer::fFactorial[9] = {1.,1.,2.,6.,24.,120.,720.,5040.,40320.};
int AliPtContainer::fSign[9] = {1,-1,1,-1,1,-1,1,-1,1};

using namespace PtSpace;
AliPtContainer::AliPtContainer():
    fTermList(0),
    fCorrList(0),
    fMomentList(0),
    fCumulantList(0),
    fNormList(0),
    mpar(0),
    fEventWeight(kWeight::kWmaxperm)
{};
AliPtContainer::~AliPtContainer()
{
    delete fCorrList;
    delete fTermList;
};
AliPtContainer::AliPtContainer(const char* name, const char* title, int nbinsx, double* xbins, int m):
    TNamed(name,title),
    fTermList(0),
    fCorrList(0),
    fMomentList(0),
    fCumulantList(0),
    fNormList(0),
    mpar(m),
    fEventWeight(kWeight::kWmaxperm)
{
    Initialize(nbinsx,xbins);
};
AliPtContainer::AliPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, int m):
    TNamed(name,title),
    fTermList(0),
    fCorrList(0),
    fMomentList(0),
    fCumulantList(0),
    fNormList(0),
    mpar(m),
    fEventWeight(kWeight::kWmaxperm)
{
    Initialize(nbinsx,xlow,xhigh);
};
void AliPtContainer::Initialize(int nbinsx, double* xbins)
{
    if(fTermList) delete fTermList;
    fTermList = new TList();
    fTermList->SetOwner(kTRUE);
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    vector<AliProfileBS*> fTerms;
    fTerms.clear();
    for(int i=0;i<=mpar;++i)
    {
        fTerms.push_back(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xbins));
        fTermList->Add(fTerms[i]);
    }
    for(int m=0;m<mpar;++m)
    {
        fCorrList->Add(new AliProfileBS(Form("corr_%ipar",m+1),this->GetTitle(),nbinsx,xbins));
    }

};
void AliPtContainer::Initialize(int nbinsx, double xlow, double xhigh)
{
    if(fTermList) delete fTermList;
    fTermList = new TList();
    fTermList->SetOwner(kTRUE);
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    vector<AliProfileBS*> fTerms;
    fTerms.clear();
    for(int i=0;i<=mpar;++i)
    {
        fTerms.push_back(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xlow,xhigh));
        fTermList->Add(fTerms[i]);
    }
    for(int m=0;m<mpar;++m)
    {
        fCorrList->Add(new AliProfileBS(Form("corr_%ipar",m+1),this->GetTitle(),nbinsx,xlow,xhigh));
    }
};
void AliPtContainer::InitializeSubsamples(const int &nsub)
{
    if(!fTermList || !fCorrList) return;
    for(int i=0; i<fTermList->GetEntries();++i)
        ((AliProfileBS*)fTermList->At(i))->InitializeSubsamples(nsub);
    for(int i=0; i<fCorrList->GetEntries();++i)
        ((AliProfileBS*)fCorrList->At(i))->InitializeSubsamples(nsub);
    return;
};
void AliPtContainer::FillRecursive(double inarr[10][10],const double &lMult, const double &rn)
{
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
    if(sumw[m]==0) return;  //Need all weights to be non-zero, so skips filling if any weight is zero
  }
  FillRecursiveProfiles(corr,sumw,lMult,rn);
  return;
}
void AliPtContainer::FillRecursiveProfiles(const vector<double> &corr, const vector<double> &sumw, const double &lMult, const double &rn)
{
    for(int m=1;m<=mpar;++m)
    {
        ((AliProfileBS*)fCorrList->At(m-1))->FillProfile(lMult,corr[m]/sumw[m],sumw[GetWeightIndex(m)],rn);
    }
    return;
}
double AliPtContainer::OrderedAddition(vector<double> vec, int size)
{
  double sum = 0;
  std::sort(vec.begin(), vec.end());

  for(int i = 0; i < size; i++)
  {
    sum += vec[i];
  }
  return sum;
}
void AliPtContainer::FillObs(double inarr[10][10], const double &lMult, const double &rn)
{
    switch(mpar)
    {
        case 1:
            FillMpt(inarr,lMult,rn);
            break;
        case 2:
            FillCk(inarr,lMult,rn);
            break;
        case 3:
            FillSkew(inarr,lMult,rn);
            break;
        case 4:
            FillKurtosis(inarr,lMult,rn);
            break;
        default:
            break;
    }
    return;
}
void AliPtContainer::FillMpt(double inarr[10][10], const double &lMult, const double &rn)
{
    Double_t lwpt = inarr[1][0];
    if(lwpt==0) return;
    ((AliProfileBS*)fTermList->At(1))->FillProfile(lMult,(inarr[1][1])/lwpt,lwpt,rn);
    return;
}
void AliPtContainer::FillCk(double inarr[10][10], const double &lMult, const double &rn)
{
    Double_t lwck = inarr[1][0]*inarr[1][0] - inarr[2][0];
    Double_t lwpt = inarr[1][0];
    if(lwck==0 || lwpt==0) return;
    ((AliProfileBS*)fTermList->At(0))->FillProfile(lMult,(inarr[1][1]*inarr[1][1] - inarr[2][2])/lwck,lwck,rn);
    ((AliProfileBS*)fTermList->At(1))->FillProfile(lMult,(inarr[1][1]*inarr[1][0]-inarr[2][1])/lwck,lwck,rn);
    ((AliProfileBS*)fTermList->At(2))->FillProfile(lMult,inarr[1][1]/lwpt,lwpt,rn);
    return;
}
void AliPtContainer::FillSkew(double inarr[10][10], const double &lMult, const double &rn)
{
    double lwpt = inarr[1][0];
    double lwskew = lwpt*lwpt*lwpt - 3*inarr[2][0]*inarr[1][0]+2*inarr[3][0];
    if(lwskew==0 || lwpt==0) return;
    ((AliProfileBS*)fTermList->At(0))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*inarr[1][1] - 3*inarr[2][2]*inarr[1][1]+2*inarr[3][3])/lwskew,lwskew,rn);
    ((AliProfileBS*)fTermList->At(1))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*inarr[1][0]-2*inarr[2][1]*inarr[1][1]+2*inarr[3][2]-inarr[2][2]*inarr[1][0])/lwskew,lwskew,rn);
    ((AliProfileBS*)fTermList->At(2))->FillProfile(lMult,(inarr[1][1]*inarr[1][0]*inarr[1][0]-2*inarr[2][1]*inarr[1][0]+2*inarr[3][1]-inarr[1][1]*inarr[2][0])/lwskew,lwskew,rn);
    ((AliProfileBS*)fTermList->At(3))->FillProfile(lMult,inarr[1][1]/lwpt,lwpt,rn);
    return;
}
void AliPtContainer::FillKurtosis(double inarr[10][10], const double &lMult, const double &rn)
{
    double lwpt = inarr[1][0];
    double lwkur = lwpt*lwpt*lwpt*lwpt-6*inarr[2][0]*lwpt*lwpt+8*lwpt*inarr[3][0]+3*inarr[2][0]*inarr[2][0]-6*inarr[4][0];
    if(lwkur==0 || lwpt==0) return;
    ((AliProfileBS*)fTermList->At(0))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*inarr[1][1]*inarr[1][1]-6*inarr[2][2]*inarr[1][1]*inarr[1][1]
                                                    +3*inarr[2][2]*inarr[2][2]+8*inarr[1][1]*inarr[3][3]-6*inarr[4][4])/lwkur,lwkur,rn);
    ((AliProfileBS*)fTermList->At(1))->FillProfile(lMult,(inarr[1][1]*inarr[1][1]*inarr[1][1]*lwpt-3*inarr[2][2]*inarr[1][1]*lwpt
                                                    -3*inarr[1][1]*inarr[1][1]*inarr[2][1]+3*inarr[2][2]*inarr[2][1]
                                                    +6*inarr[1][1]*inarr[3][2]+2*inarr[3][3]*lwpt-6*inarr[4][3])/lwkur,lwkur,rn);
    ((AliProfileBS*)fTermList->At(2))->FillProfile(lMult,(-inarr[2][2]*lwpt*lwpt
                                                    -4*inarr[1][1]*inarr[2][1]*lwpt+inarr[2][2]*inarr[2][0]+2*inarr[2][1]*inarr[2][1]
                                                    +4*inarr[1][1]*inarr[3][1]+4*inarr[3][2]*lwpt-6*inarr[4][2])/lwkur,lwkur,rn);
    ((AliProfileBS*)fTermList->At(3))->FillProfile(lMult,(inarr[1][1]*lwpt*lwpt*lwpt-3*inarr[2][1]*lwpt*lwpt
                                                    -3*inarr[1][1]*inarr[2][0]*lwpt+3*inarr[2][1]*inarr[2][0]+2*inarr[1][1]*inarr[3][0]
                                                    +6*inarr[3][1]*lwpt-6*inarr[4][1])/lwkur,lwkur,rn);
    ((AliProfileBS*)fTermList->At(4))->FillProfile(lMult,inarr[1][1]/lwpt,lwpt,rn);
    return;
}
TH1* AliPtContainer::getHist(int ind)
{
  if(!fMomentList) CalculateObs();
  if(!fMomentList) return 0;
  if(ind+1<fMomentList->GetEntries()) return (TH1*)fMomentList->At(ind+1);
  return 0;
}
void AliPtContainer::CalculateObs()
{
  if(fMomentList) delete fMomentList;
  fMomentList = new TList();
  fMomentList->SetOwner(kTRUE);
  //Override mpt temporarily override mpt weights:
  ((AliProfileBS*)fTermList->At(mpar))->PresetWeights((AliProfileBS*)fTermList->At(0));
  for(Int_t i=-1;i<((AliProfileBS*)fTermList->At(0))->getNSubs();i++) {
    vector<TH1*> hTs;
    for(Int_t j=0;j<=mpar;j++) {
      ((AliProfileBS*)fTermList->At(j))->SetErrorOption("g");
      hTs.push_back(((AliProfileBS*)fTermList->At(j))->getHist(i));
    }
    TH1 *hRec = RecalculateObsHists(hTs);
    hRec->SetName(Form("%s_Recalculated%i",this->GetName(),i+1));
    fMomentList->Add(hRec);
  };
  //reset mpt weights
  ((AliProfileBS*)fTermList->At(mpar))->PresetWeights(0);
}
TH1 *AliPtContainer::RecalculateObsHists(vector<TH1*> inh)
{
    TH1* reth = (TH1*)inh[mpar]->Clone("reth");
    switch(mpar)
    {
        case 1:
            reth = inh[0];
            break;
        case 2:
            reth = RecalculateCkHists(inh);
            break;
        case 3:
            reth = RecalculateSkewHists(inh);
            break;
        case 4:
            reth = RecalculateKurtosisHists(inh);
            break;
        default:
            break;
    }
    return reth;
}
TH1 *AliPtContainer::RecalculateCkHists(vector<TH1*> inh)
{
  inh[1]->Multiply(inh[2]);
  inh[1]->Scale(2);
  TH1 *mptsq = (TH1*)inh[2]->Clone("mptSquared");
  mptsq->Multiply(inh[2]);
  TH1 *hWeights = (TH1*)inh[0]->Clone("ForErrors");
  inh[0]->Add(inh[1],-1);
  inh[0]->Add(mptsq);
  delete mptsq;
  for(Int_t i=1;i<=hWeights->GetNbinsX();i++) inh[0]->SetBinError(i,hWeights->GetBinError(i));
  delete hWeights;
  return (TH1*)inh[0]->Clone("hRec");
}
TH1 *AliPtContainer::RecalculateSkewHists(vector<TH1*> inh)
{
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
TH1 *AliPtContainer::RecalculateKurtosisHists(vector<TH1*> inh)
{
  inh[1]->Multiply(inh[4]);
  inh[1]->Scale(4);
  TH1 *mptsq = (TH1*)inh[4]->Clone("mptSquared");
  TH1 *mptcb = (TH1*)inh[4]->Clone("mptCubed");
  TH1 *mpt4th = (TH1*)inh[4]->Clone("mpt4th");
  mptsq->Multiply(inh[4]);
  inh[2]->Multiply(mptsq);
  inh[2]->Scale(6);
  mptcb->Multiply(mptsq);
  inh[3]->Multiply(mptcb);
  inh[3]->Scale(4);
  mpt4th->Multiply(mptcb);
  TH1 *hWeights = (TH1*)inh[0]->Clone("ForErrors");
  inh[0]->Add(inh[1]-1);
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
TH1* AliPtContainer::getRercusiveHist(int ind, int m, unsigned int l_obs)
{
  if(l_obs==kObs::kCorr) return ((AliProfileBS*)fCorrList->At(m-1))->getHist(ind);
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
void AliPtContainer::CalculateRecursive(bool normalized)
{ 
    if(normalized)
    {
        if(fNormList) delete fNormList;
        fNormList = new TList();
        fNormList->SetOwner(); 
        if(!fCumulantList)
        {
            fCumulantList = new TList();
            fCumulantList->SetOwner();  
        }
    }
    else
    {
        if(fCumulantList) delete fCumulantList;
        fCumulantList = new TList();
        fCumulantList->SetOwner();
    }
    for(int i=0;i<mpar;++i)
        ((AliProfileBS*)fCorrList->At(i))->PresetWeights((AliProfileBS*)fTermList->At(mpar));
    for(int i=-1;i<((AliProfileBS*)fCorrList->At(0))->getNSubs();++i) {
        vector<TH1*> hTs;
        for(int j=0;j<mpar;++j) {
            ((AliProfileBS*)fCorrList->At(j))->SetErrorOption("g");
            hTs.push_back(((AliProfileBS*)fCorrList->At(j))->getHist(i));
        }
        CalculateCumulantHists(hTs,i,normalized);
    }
    ((AliProfileBS*)fTermList->At(0))->PresetWeights(0);
    return;
}
void AliPtContainer::CalculateCumulantHists(vector<TH1*> inh, int ind, bool normalized) 
{
    auto binomial = [&](const int n, const int m) { return factorial(n)/(factorial(m)*factorial(n-m)); };
    if((normalized && fCumulantList->GetEntries()<1) || !normalized)
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
    else
    {
        TH1* sigmah = (TH1*)fCumulantList->At((ind+1)*mpar+1)->Clone(Form("sigmah%i",ind));
        for(int m = 0; m<mpar;++m)
        {
            TH1* normh = (TH1*)fCumulantList->At((ind+1)*mpar+m-1)->Clone(Form("normh%i%i",m,ind));
            normh->Divide(getPowerHist(sigmah,0.5*(m+1)));
            fNormList->Add(normh);
        }
    }

    return;
}
TH1* AliPtContainer::getPowerHist(TH1* inh, double p)
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
int AliPtContainer::GetWeightIndex(int m)
{
    if(fEventWeight==kOne) return 0;
    if(fEventWeight==kWpt) return 1;
    if(fEventWeight==kWperms) return m;
    if(fEventWeight==kWmaxperm) return mpar;
    return m;
}
Long64_t AliPtContainer::Merge(TCollection *collist) {
  Long64_t nmerged=0;
  TIter all_PTC(collist);
  AliPtContainer *l_PTC = 0;
  while ((l_PTC = ((AliPtContainer*) all_PTC()))) {
      TList *t_Term = l_PTC->fTermList;
      TList *t_Mom  = l_PTC->fMomentList;
      TList* t_Corr = l_PTC->fCorrList;
      TList* t_Cum = l_PTC->fCumulantList;
      TList* t_Norm = l_PTC->fNormList;
      if(t_Term) {
        if(!fTermList) fTermList = (TList*)t_Term->Clone();
        else MergeBSLists(fTermList,t_Term);
        nmerged++;
      };
      if(t_Mom) {
        if(!fMomentList) fMomentList = (TList*)t_Mom->Clone();
        else MergeBSLists(fMomentList,t_Mom);
      };
      if(t_Corr) {
          if(!fCorrList) fCorrList = (TList*)t_Corr->Clone();
          else MergeBSLists(fCorrList,t_Corr);
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