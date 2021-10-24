/*
Author: Vytautas Vislavicius
Container to store the terms required to calculate higher order moments of the pT spectrum.
Extra layer to calculate the moments. Current implementation supports only the first and second moments.
If used, modified, or distributed, please aknowledge the original author of this code.
*/

#include "AliCkContainer.h"
AliCkContainer::AliCkContainer():
  fObsList(0),
  fCkList(0) {
};
AliCkContainer::~AliCkContainer() {
  delete fObsList;
  delete fCkList;
};
AliCkContainer::AliCkContainer(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins):
  TNamed(name, title),
  fObsList(0),
  fCkList(0)
{
  Initialize(nbinsx,xbins);
}
AliCkContainer::AliCkContainer(const char* name, const char* title, Int_t nbinsx, Double_t xmin, Double_t xmax):
  TNamed(name, title),
  fObsList(0),
  fCkList(0)
{
  Initialize(nbinsx,xmin,xmax);
}
void AliCkContainer::Initialize(Int_t nbinsx, const Double_t* xbins) {
  if(fObsList) delete fObsList;
  fObsList = new TList();
  fObsList->SetOwner(kTRUE);
  AliProfileBS *consterm = new AliProfileBS(Form("%s_consterm",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(consterm);
  AliProfileBS *linterm = new AliProfileBS(Form("%s_linterm",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(linterm);
  AliProfileBS *lmpt = new AliProfileBS(Form("%s_mpt",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(lmpt);
};
void AliCkContainer::Initialize(Int_t nbinsx, Double_t xmin, Double_t xmax) {
  if(fObsList) delete fObsList;
  fObsList = new TList();
  fObsList->SetOwner(kTRUE);
  AliProfileBS *consterm = new AliProfileBS(Form("%s_consterm",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(consterm);
  AliProfileBS *linterm = new AliProfileBS(Form("%s_linterm",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(linterm);
  AliProfileBS *lmpt = new AliProfileBS(Form("%s_mpt",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(lmpt);
}
void AliCkContainer::InitializeSubsamples(const Int_t &nrb) {
  if(!fObsList) return;
  for(Int_t i=0;i<fObsList->GetEntries();i++) {
    ((AliProfileBS*)fObsList->At(i))->InitializeSubsamples(nrb);
  }
}
void AliCkContainer::FillObs(const Double_t inArr[5], const Double_t &lMult, const Double_t &rn) {
  // [0] = w1p0
  // [1] = w1p1
  // [2] = w2p2
  // [3] = w2p1
  // [4] = w2p0
  /* Full expression for ck:
0)  <ck> = <wpsq>            = [1]*[1]
1)         -2 mpt <wp * w>   = [1]*[0] (*-2mpt)
2)         +mpt2 <wsq>       = [0]*[0]     (*mpt2)
3)         -<w2p2>           = [2]     (*-1)
4)         +2 mpt <w2p>      = [3]     (*2mpt)
5)         - mpt2 <w2>       = [4]     (*-mpt2)
  //To reduce number of profiles, collect terms w/ the same power of [pt]:
0) <ck>     = <wpsq - w2p2>            = [1]*[1] - [2]
1) +2[pt]*    <w2p1 - w1p1 * w1p0>     = [3] - [1]*[0]
2) +[pt]^2 *  < wsq - w2p0>            = [0]*[0] - [4] ===> This is not even required, since by construction this is one (after normalizing by the weight, which is exactly the same)
  */
  Double_t lwck = inArr[0]*inArr[0] - inArr[4];
  Double_t lwpt = inArr[0];
  if(lwck==0 || lwpt==0) return;
  ((AliProfileBS*)fObsList->At(0))->FillProfile(lMult,(inArr[1]*inArr[1] - inArr[2])/lwck,lwck,rn);
  ((AliProfileBS*)fObsList->At(1))->FillProfile(lMult,(inArr[3] - inArr[1]*inArr[0])/lwck,lwck,rn);
  ((AliProfileBS*)fObsList->At(2))->FillProfile(lMult,inArr[1]/inArr[0],lwpt,rn);
}
TH1 *AliCkContainer::RecalculateCkHists(TH1 **inh) {
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
void AliCkContainer::CalculateCks() {
  if(fCkList) delete fCkList;
  fCkList = new TList();
  fCkList->SetOwner(kTRUE);
  //Override mpt temporarily override mpt weights:
  ((AliProfileBS*)fObsList->At(2))->PresetWeights((AliProfileBS*)fObsList->At(0));
  for(Int_t i=-1;i<((AliProfileBS*)fObsList->At(0))->getNSubs();i++) {
    TH1 *hTs[3];
    for(Int_t j=0;j<3;j++) {
      ((AliProfileBS*)fObsList->At(j))->SetErrorOption("g");
      hTs[j] = ((AliProfileBS*)fObsList->At(j))->getHist(i);
    }
    TH1 *hCkRec = RecalculateCkHists(hTs);
    hCkRec->SetName(Form("Ck_Recalculated%i",i+1));
    fCkList->Add(hCkRec);
    for(Int_t j=0;j<3;j++) delete hTs[j];
  };
  //reset mpt weights
  ((AliProfileBS*)fObsList->At(2))->PresetWeights(0);
}
TH1 *AliCkContainer::getHist(Int_t ind) {
  if(!fCkList) CalculateCks();
  if(!fCkList) return 0;
  if(ind+1<fCkList->GetEntries()) return (TH1*)fCkList->At(ind+1);
  return 0;
}
Long64_t AliCkContainer::Merge(TCollection *collist) {
  Long64_t nmerged=0;
  AliCkContainer *l_CKC = 0;
  TIter all_CKC(collist);
  while ((l_CKC = ((AliCkContainer*) all_CKC()))) {
    TList *tarL = l_CKC->fObsList;
    if(!tarL) continue;
    if(!fObsList) {
      fObsList = (TList*)tarL->Clone();
    } else
      for(Int_t i=0; i<fObsList->GetEntries(); i++) ((AliProfileBS*)fObsList->At(i))->MergeBS((AliProfileBS*)tarL->At(i));
    nmerged++;
  };
  return nmerged;
};
