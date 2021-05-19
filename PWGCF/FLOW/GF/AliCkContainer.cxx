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
  AliProfileBS *lwpsq = new AliProfileBS(Form("%s_wpsq",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(lwpsq);
  AliProfileBS *lwpw = new AliProfileBS(Form("%s_wpw",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(lwpw);
  AliProfileBS *lwsq = new AliProfileBS(Form("%s_wsq",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(lwsq);
  AliProfileBS *lw2p2 = new AliProfileBS(Form("%s_w2p2",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(lw2p2);
  AliProfileBS *lw2p = new AliProfileBS(Form("%s_w2p",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(lw2p);
  AliProfileBS *lw2 = new AliProfileBS(Form("%s_w2",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(lw2);
  AliProfileBS *lmpt = new AliProfileBS(Form("%s_mpt",this->GetName()),this->GetTitle(),nbinsx,xbins);
  fObsList->Add(lmpt);
};
void AliCkContainer::Initialize(Int_t nbinsx, Double_t xmin, Double_t xmax) {
  if(fObsList) delete fObsList;
  fObsList = new TList();
  fObsList->SetOwner(kTRUE);
  AliProfileBS *lwpsq = new AliProfileBS(Form("%s_wpsq",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(lwpsq);
  AliProfileBS *lwpw = new AliProfileBS(Form("%s_wpw",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(lwpw);
  AliProfileBS *lwsq = new AliProfileBS(Form("%s_wsq",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(lwsq);
  AliProfileBS *lw2p2 = new AliProfileBS(Form("%s_w2p2",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(lw2p2);
  AliProfileBS *lw2p = new AliProfileBS(Form("%s_w2p",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(lw2p);
  AliProfileBS *lw2 = new AliProfileBS(Form("%s_w2",this->GetName()),this->GetTitle(),nbinsx,xmin,xmax);
  fObsList->Add(lw2);
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
  */
  Double_t lwck = inArr[0]*inArr[0] - inArr[4];
  Double_t lwpt = inArr[0];
  if(lwck==0 || lwpt==0) return;
  ((AliProfileBS*)fObsList->At(0))->FillProfile(lMult,(inArr[1]*inArr[1])/lwck,lwck,rn);
  ((AliProfileBS*)fObsList->At(1))->FillProfile(lMult,(inArr[1]*inArr[0])/lwck,lwck,rn);
  ((AliProfileBS*)fObsList->At(2))->FillProfile(lMult,(inArr[0]*inArr[0])/lwck,lwck,rn);
  ((AliProfileBS*)fObsList->At(3))->FillProfile(lMult,(inArr[2])/lwck,lwck,rn);
  ((AliProfileBS*)fObsList->At(4))->FillProfile(lMult,(inArr[3])/lwck,lwck,rn);
  ((AliProfileBS*)fObsList->At(5))->FillProfile(lMult,(inArr[4])/lwck,lwck,rn);
  ((AliProfileBS*)fObsList->At(6))->FillProfile(lMult,inArr[1]/inArr[0],lwpt,rn);
}
TH1 *AliCkContainer::RecalculateCkHists(TH1 **inh) {
  inh[1]->Multiply(inh[6]);
  inh[1]->Scale(-2);
  inh[2]->Multiply(inh[6]);
  inh[2]->Multiply(inh[6]);
  inh[3]->Scale(-1);
  inh[4]->Multiply(inh[6]);
  inh[4]->Scale(2);
  inh[5]->Multiply(inh[6]);
  inh[5]->Multiply(inh[6]);
  inh[5]->Scale(-1);
  for(Int_t i=1;i<6;i++) inh[0]->Add(inh[i]);
  for(Int_t i=1;i<=inh[3]->GetNbinsX();i++) inh[0]->SetBinError(i,inh[3]->GetBinError(i));
  return (TH1*)inh[0]->Clone("hRec");
}
void AliCkContainer::CalculateCks() {
  if(fCkList) delete fCkList;
  fCkList = new TList();
  fCkList->SetOwner(kTRUE);
  for(Int_t i=-1;i<((AliProfileBS*)fObsList->At(0))->getNSubs();i++) {
    TH1 *hTs[7];
    for(Int_t j=0;j<7;j++) {
      ((AliProfileBS*)fObsList->At(j))->SetErrorOption("g");
      hTs[j] = ((AliProfileBS*)fObsList->At(j))->getHist(i);
    }
    TH1 *hCkRec = RecalculateCkHists(hTs);
    hCkRec->SetName(Form("Ck_Recalculated%i",i+1));
    fCkList->Add(hCkRec);
    for(Int_t j=0;j<7;j++) delete hTs[j];
  };
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
