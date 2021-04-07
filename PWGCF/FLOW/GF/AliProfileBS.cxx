#include "AliProfileBS.h"
AliProfileBS::AliProfileBS():
  TProfile(),
  fListOfEntries(0),
  fProfInitialized(kFALSE),
  fNSubs(0)
{
};
AliProfileBS::~AliProfileBS() {
  delete fListOfEntries;
};
AliProfileBS::AliProfileBS(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins):
  TProfile(name,title,nbinsx,xbins),
  fListOfEntries(0),
  fProfInitialized(kTRUE),
  fNSubs(0)
{};
AliProfileBS::AliProfileBS(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup):
  TProfile(name,title,nbinsx,xlow,xup),
  fListOfEntries(0),
  fProfInitialized(kFALSE),
  fNSubs(0)
{};
void AliProfileBS::InitializeSubsamples(Int_t nSub) {
  if(nSub<1) {printf("Number of subprofiles has to be > 0!\n"); return; };
  if(fListOfEntries) delete fListOfEntries;
  fListOfEntries = new TList();
  fListOfEntries->SetOwner(kTRUE);
  TProfile *dummyPF = (TProfile*)this;
  for(Int_t i=0;i<nSub;i++) {
    fListOfEntries->Add((TProfile*)dummyPF->Clone(Form("%s_Subpf%i",dummyPF->GetName(),i)));
    ((TProfile*)fListOfEntries->At(i))->Reset();
  }
  fNSubs = nSub;
}
void AliProfileBS::FillProfile(const Double_t &xv, const Double_t& yv, const Double_t &w, const Double_t &rn) {
  TProfile::Fill(xv,yv,w);
  if(!fNSubs) return;
  Int_t targetInd = rn*fNSubs;
  if(targetInd>=fNSubs) targetInd = 0;
  ((TProfile*)fListOfEntries->At(targetInd))->Fill(xv,yv,w);
}
void AliProfileBS::FillProfile(const Double_t &xv, const Double_t &yv, const Double_t &w) {
  TProfile::Fill(xv,yv,w);
}
void AliProfileBS::RebinMulti(Int_t nbins) {
  this->RebinX(nbins);
  if(!fListOfEntries) return;
  for(Int_t i=0;i<fListOfEntries->GetEntries();i++)
    ((TProfile*)fListOfEntries->At(i))->RebinX(nbins);
}
TH1 *AliProfileBS::getHist(Int_t ind) {
  if(ind<0) {
    if((TProfile*)this) return ((TProfile*)this)->ProjectionX(Form("%s_hist",this->GetName()));
    else { printf("Empty AliProfileBS addressed, cannot get a histogram\n"); return 0; };
  } else {
    if(!fListOfEntries) { printf("No subprofiles exist!\n"); return 0; };
    if(ind<fNSubs) return ((TProfile*)fListOfEntries->At(ind))->ProjectionX(Form("%s_sub%i",((TProfile*)fListOfEntries->At(ind))->GetName(),ind));
    else { printf("Trying to fetch subprofile no %i out of %i, not possible\n",ind,fNSubs); return 0;};
  }
  return 0;
}
Long64_t AliProfileBS::Merge(TCollection *collist) {
  Long64_t nmergedpf = TProfile::Merge(collist);
  Long64_t nmerged=0;
  AliProfileBS *l_PBS = 0;
  TIter all_PBS(collist);
  while ((l_PBS = ((AliProfileBS*) all_PBS()))) {
    TList *tarL = l_PBS->fListOfEntries;
    if(!tarL) continue;
    if(!fListOfEntries) {
      fListOfEntries = (TList*)tarL->Clone();
      for(Int_t i=0; i<fListOfEntries->GetEntries(); i++) ((TProfile*)fListOfEntries->At(i))->Reset();
    };
    for(Int_t i=0; i<fListOfEntries->GetEntries(); i++) ((TProfile*)fListOfEntries->At(i))->Add((TProfile*)tarL->At(i));
    nmerged++;
  };
  return nmergedpf;
};
