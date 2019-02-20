#include "AliProfileSubset.h"
AliProfileSubset::AliProfileSubset() { 
};
AliProfileSubset::AliProfileSubset(TProfile *inpf):
  TProfile(*inpf)
{
};
AliProfileSubset::~AliProfileSubset() {
};
Bool_t AliProfileSubset::CopyFromProfile(Int_t StartI, Int_t StopI, TProfile *source) {
  AliProfileSubset *inpf = new AliProfileSubset(source);
  Int_t nTot = StopI-StartI+1;
  TArrayD fBinEntries_Arr = inpf->Get_fBinEntries();
  TArrayD fBinSumw2_Arr   = inpf->Get_fBinSumw2();
  TArrayD fSumw2_Arr      = inpf->Get_fSumw2();

  Double_t *fBinEntries_D = fBinEntries_Arr.GetArray();
  Double_t *fBinSumw2_D   = fBinSumw2_Arr.GetArray();
  Double_t *fSumw2_D      = fSumw2_Arr.GetArray();
  Double_t *arr           = inpf->Get_fArray();
  Double_t fArrBU=0;
  if(StartI>0) { //Override underflow bin
    StartI--;
    nTot++;
    fBinEntries_D[StartI] = 0;
    if(fBinSumw2_Arr.GetSize())
      fBinSumw2_D[StartI] = 0;
    fSumw2_D[StartI]      = 0;
    fArrBU=arr[StartI];
    arr[StartI]           = 0;
  };

  //Overriding local values:
  if(fBinEntries_Arr.GetSize())
    fBinEntries = TArrayD(nTot, &fBinEntries_D[StartI]);
  if(fBinSumw2_Arr.GetSize())
    fBinSumw2   = TArrayD(nTot, &fBinSumw2_D[StartI]);
  if(fSumw2_Arr.GetSize())
    fSumw2      = TArrayD(nTot, &fSumw2_D[StartI]);
  for(Int_t i=0;i<nTot;i++) fArray[i] = arr[StartI+i];
  fYmin      = inpf->GetYMin();
  fYmax      = inpf->GetYMax();
  fErrorMode = inpf->GetErrorMode();
  fTsumwy    = inpf->GetSumwy();
  fTsumwy2   = inpf->GetSumwy2();
  if(fArrBU!=0) arr[StartI] = fArrBU;
  SetEntries(inpf->GetEntries());
  delete inpf;
  return kTRUE;
  
};
