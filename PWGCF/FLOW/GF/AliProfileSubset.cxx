#include "AliProfileSubset.h"
AliProfileSubset::AliProfileSubset(TProfile2D &inpf):
  TProfile2D(inpf) {
}
TProfile* AliProfileSubset::GetSubset(Bool_t onX, const char *name, Int_t firstbin, Int_t lastbin, Int_t l_nbins, Double_t *l_binarray) {
  TString expectedName = ( onX ? "_pfx" : "_pfy" );
  TString pname(name);
  if (pname.IsNull() || name == expectedName)
     pname = TString(GetName() ) + expectedName;
  const TAxis& outAxis = ( onX ? fXaxis : fYaxis );
  const TArrayD *bins = outAxis.GetXbins();
  Int_t firstOutBin = outAxis.GetFirst();
  Int_t lastOutBin = outAxis.GetLast();
  //printf("firstOutBin = %i, lastOutBin = %i\n",firstOutBin,lastOutBin);
  TProfile  * p1 = 0;
  if(l_nbins) p1 = new TProfile(pname, GetTitle(),l_nbins, l_binarray);
  else p1 = new TProfile(pname,GetTitle(),outAxis.GetNbins(),bins->fArray);
  //printf("p1 has %i bins (l_nbins = %i)\n",p1->GetNbinsX(),l_nbins);
  if (fBinSumw2.fN) p1->Sumw2();
  TH2D * h2dW = ProjectionXY("h2temp-W","W");
  TH2D * h2dN = ProjectionXY("h2temp-N","B");
  h2dW->SetDirectory(0); h2dN->SetDirectory(0);
  if(onX) {
    h2dW->GetXaxis()->SetRange(firstOutBin,lastOutBin);
    h2dN->GetXaxis()->SetRange(firstOutBin,lastOutBin);
  } else {
    h2dW->GetYaxis()->SetRange(firstOutBin,lastOutBin);
    h2dN->GetYaxis()->SetRange(firstOutBin,lastOutBin);
  }
  TH1D * h1W = (onX) ? h2dW->ProjectionX("h1temp-W",firstbin,lastbin) : h2dW->ProjectionY("h1temp-W",firstbin,lastbin);
  TH1D * h1N = (onX) ? h2dN->ProjectionX("h1temp-N",firstbin,lastbin) : h2dN->ProjectionY("h1temp-N",firstbin,lastbin);
  h1W->SetDirectory(0); h1N->SetDirectory(0);
  //printf("Asserting, %i vs. %i\n",h1W->fN, p1->fN);
  R__ASSERT( h1W->fN == p1->fN );
  R__ASSERT( h1N->fN == p1->fN );
  R__ASSERT( h1W->GetSumw2()->fN != 0); // h1W should always be a weighted histogram since h2dW is
  for (int i = 0; i < p1->fN ; ++i) {
     p1->fArray[i] = h1W->GetBinContent(i);   // array of profile is sum of all values
     p1->GetSumw2()->fArray[i]  = h1W->GetSumw2()->fArray[i];   // array of content square of profile is weight square of the W projected histogram
     p1->SetBinEntries(i, h1N->GetBinContent(i) );
     if (fBinSumw2.fN) p1->GetBinSumw2()->fArray[i] = h1N->GetSumw2()->fArray[i];    // sum of weight squares are stored to compute errors in h1N histogram
  }
  delete h2dW;
  delete h2dN;
  delete h1W;
  delete h1N;
  p1->SetEntries( p1->GetEffectiveEntries() );
  return p1;
};
void AliProfileSubset::OverrideBinContent(Double_t x, Double_t y, Double_t x2, Double_t y2, Double_t val) {
  if (!fBinSumw2.fN) Sumw2();
  TH2D * h2dW = ProjectionXY("h2temp-W","W");
  TH2D * h2dN = ProjectionXY("h2temp-N","B");
  Int_t binIndex = FindBin(x,y);
  Int_t binIndex2 = FindBin(x2,y2);
  fArray[binIndex] = h2dW->GetBinContent(binIndex2);
  GetSumw2()->fArray[binIndex] = h2dW->GetSumw2()->fArray[binIndex2];
  SetBinEntries(binIndex,h2dN->GetBinContent(binIndex2));
  if(fBinSumw2.fN) GetBinSumw2()->fArray[binIndex] = h2dN->GetSumw2()->fArray[binIndex2];
}
