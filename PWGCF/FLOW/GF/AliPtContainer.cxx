/*
Author: Emil Gorm Nielsen
Modified from AliCkContainer by Vytautas
Container to store the terms required to calculate higher order moments of the pT spectrum.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#include "AliPtContainer.h"

AliPtContainer::AliPtContainer():
    fTermList(0),
    fObsList(0)
{};
AliPtContainer::~AliPtContainer()
{
    delete fTermList;
};
AliPtContainer::AliPtContainer(const char* name, const char* title, int nbinsx, double* xbins, int m):
    TNamed(name,title),
    fTermList(0),
    fObsList(0),
    mpar(m)
{
    Initialize(nbinsx,xbins);
};
AliPtContainer::AliPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, int m):
    TNamed(name,title),
    fTermList(0),
    fObsList(0),
    mpar(m)
{
    Initialize(nbinsx,xlow,xhigh);
};
void AliPtContainer::Initialize(int nbinsx, double* xbins)
{
    if(fTermList) delete fTermList;
    fTermList = new TList();
    fTermList->SetOwner(kTRUE);
    vector<AliProfileBS*> fTerms;
    fTerms.clear();
    for(int i=0;i<=mpar;++i)
    {
        fTerms.push_back(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xbins));
        fTermList->Add(fTerms[i]);
    }
};
void AliPtContainer::Initialize(int nbinsx, double xlow, double xhigh)
{
    if(fTermList) delete fTermList;
    fTermList = new TList();
    fTermList->SetOwner(kTRUE);
    vector<AliProfileBS*> fTerms;
    fTerms.clear();
    for(int i=0;i<=mpar;++i)
    {
        fTerms.push_back(new AliProfileBS(Form("%s_p%i",this->GetName(),i),this->GetTitle(),nbinsx,xlow,xhigh));
        fTermList->Add(fTerms[i]);
    }
};
void AliPtContainer::InitializeSubsamples(const int &nsub)
{
    if(!fTermList) return;
    for(int i=0; i<fTermList->GetEntries();++i)
        ((AliProfileBS*)fTermList->At(i))->InitializeSubsamples(nsub);
    return;
};
void AliPtContainer::FillObs(double inarr[15][15], const double &lMult, const double &rn)
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
void AliPtContainer::FillMpt(double inarr[15][15], const double &lMult, const double &rn)
{
    Double_t lwpt = inarr[1][0];
    if(lwpt==0) return;
    ((AliProfileBS*)fTermList->At(1))->FillProfile(lMult,(inarr[1][1])/lwpt,lwpt,rn);
    return;
}
void AliPtContainer::FillCk(double inarr[15][15], const double &lMult, const double &rn)
{
    Double_t lwck = inarr[1][0]*inarr[1][0] - inarr[2][0];
    Double_t lwpt = inarr[1][0];
    if(lwck==0 || lwpt==0) return;
    ((AliProfileBS*)fTermList->At(0))->FillProfile(lMult,(inarr[1][1]*inarr[1][1] - inarr[2][2])/lwck,lwck,rn);
    ((AliProfileBS*)fTermList->At(1))->FillProfile(lMult,(-inarr[2][1]+inarr[1][1]*inarr[1][0])/lwck,lwck,rn);
    ((AliProfileBS*)fTermList->At(2))->FillProfile(lMult,inarr[1][1]/lwpt,lwpt,rn);
    return;
}
void AliPtContainer::FillSkew(double inarr[15][15], const double &lMult, const double &rn)
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
void AliPtContainer::FillKurtosis(double inarr[15][15], const double &lMult, const double &rn)
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
  if(!fObsList) CalculateObs();
  if(!fObsList) return 0;
  if(ind+1<fObsList->GetEntries()) return (TH1*)fObsList->At(ind+1);
  return 0;
}
void AliPtContainer::CalculateObs()
{
  if(fObsList) delete fObsList;
  fObsList = new TList();
  fObsList->SetOwner(kTRUE);
  //Override mpt temporarily override mpt weights:
  ((AliProfileBS*)fObsList->At(mpar))->PresetWeights((AliProfileBS*)fObsList->At(0));
  for(Int_t i=-1;i<((AliProfileBS*)fObsList->At(0))->getNSubs();i++) {
    vector<TH1*> hTs;
    for(Int_t j=0;j<=mpar;j++) {
      ((AliProfileBS*)fObsList->At(j))->SetErrorOption("g");
      hTs.push_back(((AliProfileBS*)fObsList->At(j))->getHist(i));
    }
    TH1 *hRec = RecalculateObsHists(hTs);
    hRec->SetName(Form("%s_Recalculated%i",this->GetName(),i+1));
    fObsList->Add(hRec);
  };
  //reset mpt weights
  ((AliProfileBS*)fObsList->At(mpar))->PresetWeights(0);
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
Long64_t AliPtContainer::Merge(TCollection *collist) {
  Long64_t nmerged=0;
  TIter all_PTC(collist);
  AliPtContainer *l_PTC = 0;
  while ((l_PTC = ((AliPtContainer*) all_PTC()))) {
      TList *t_Term = l_PTC->fTermList;
      TList *t_Obs  = l_PTC->fObsList;
      if(t_Term) {
        if(!fTermList) fTermList = (TList*)t_Term->Clone();
        else MergeBSLists(fTermList,t_Term);
        nmerged++;
      };
      if(t_Obs) {
        if(!fObsList) fObsList = (TList*)t_Obs->Clone();
        else MergeBSLists(fObsList,t_Obs);
      };
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
