

#ifndef ALICHARGEONEPWRTRP_H
#define ALICHARGEONEPWRTRP_H

#include <TMath.h>
#include <THn.h>
#include <TProfile.h>
#include <AliQnCorrectionsQnVector.h>

class AliCMEAnalysisCuts;

#include <iostream>
using std::cout;
using std::endl;

//_________________________________________________________________________
class AliChargeOnePwrtRP : public TNamed {

 public:
  AliChargeOnePwrtRP();
  ~AliChargeOnePwrtRP();

  enum Constants {
      Nharmonics=4,
      Nqvectors=4
  };

  void SetBin(Int_t Ct,Int_t Pt,Int_t mPt,Int_t dPt,Int_t Eta,Int_t mEta,Int_t dEta) {fBins[0]=Ct;fBins[1]=Pt;fBins[2]=mPt;fBins[3]=dPt;fBins[4]=Eta;fBins[5]=mEta;fBins[6]=dEta;}

  //Float_t GetValue(Int_t ih, Int_t comp, Int_t bin) {return fQvector[ih][comp]->At(bin);}
  //Int_t Multiplicity(Int_t bin) {return fMultiplicity->At(bin);}
  Float_t GetValue(Int_t ih, Int_t comp, Int_t bin, Int_t charge) {return fQvector[charge][ih-1][comp][bin];}
  Int_t Multiplicity(Int_t bin,Int_t charge) {return fMultiplicity[charge][bin];}
  Short_t Charge(Int_t ch1) {return (ch1==1 ? 0 : 1);}
  void AddTrack(Int_t charge, Float_t x[Nharmonics], Float_t y[Nharmonics]);
  Int_t GetMinHarmonic() const {return fMinHarmonic;}
  Int_t GetMaxHarmonic() const {return fMaxHarmonic;}
  AliCMEAnalysisCuts* GetTrackCuts() {return fTrackCutsOnePwrtRP;}
  void Clear();

  THnF* GetTHn(Int_t ep , Int_t charge, Int_t ih, Int_t comp )   {return  fCorrelationTHn[ep][charge][ih-1][comp];}
  THnF* GetTHn2(Int_t ep, Int_t charge, Int_t ih, Int_t comp )   {return fCorrelationTHn2[ep][charge][ih-1][comp];}
  THnF* GetTHnMult(Int_t ep, Int_t charge)   {return fCorrelationTHnMult[ep][charge];}

  THnF* GetSumQ(Int_t ep , Int_t charge, Int_t ih, Int_t comp )   {return  fSumQ[ep][charge][ih-1][comp];}
  THnF* GetSumQ2(Int_t ep , Int_t charge, Int_t ih, Int_t comp )   {return  fSumQ2[ep][charge][ih-1][comp];}
  THnF* GetSumQMult(Int_t ep, Int_t charge)   {return fSumQMult[ep][charge];}

  void FillCorrelation(Int_t bin, Int_t thnbin, AliQnCorrectionsQnVector* qvecEP, Int_t ep);
  void SetTHn(TString eventSelection, TString tracking1Quality, TString p1,  TString d1, TString d2="", Char_t* ar[]=0x0);
  void Init(TAxis* x, TAxis* y, TAxis* z);
  void SetVars(Int_t varX, Int_t varY=0, Int_t varZ=0) {fVarX=varX;fVarY=varY;fVarZ=varZ;}
  void SetRangeHarmonics(Int_t min, Int_t max) {fMinHarmonic=min; fMaxHarmonic=max;}
  void SetAxes(TAxis x, TAxis y, TAxis z)  {fBinLimits[0]=TAxis(x);fBinLimits[1]=TAxis(y);fBinLimits[2]=TAxis(z);};
  void SetTrackCuts(AliCMEAnalysisCuts* p) {fTrackCutsOnePwrtRP = p;}

  Int_t GetVarX() {return fVarX;}
  Int_t GetVarY() {return fVarY;}
  Int_t GetVarZ() {return fVarZ;}
  Int_t GetTrackVarX() {return fTrackVarX;}
  Int_t GetTrackVarY() {return fTrackVarY;}
  Int_t GetTrackVarZ() {return fTrackVarZ;}
  Int_t GetTrackVarMap(Int_t i) {return fTrackVarMap[i];}
  Int_t GetNdim() {return fNdim;}
  //Int_t Bin() const {return fBins[fVarX]+fBins[fVarY]*fNbinsX+fBins[fVarZ]*fNbinsX*fNbinsY;}
  Int_t Bin() const {return fBins[fTrackVarX]*fNbinsY*fNbinsZ+fBins[fTrackVarY]*fNbinsZ+fBins[fTrackVarZ];}
  //Int_t Bin(Int_t ix, Int_t iy=0, Int_t iz=0) const {return ix+iy*fNbinsX+iz*fNbinsX*fNbinsY;}
  Int_t GetNbinsX() {return fNbinsX;}
  Int_t GetNbinsY() {return fNbinsY;}
  Int_t GetNbinsZ() {return fNbinsZ;}
  Int_t Nbins(Int_t i) const {return fNbin[i];}





 private:
  Int_t fNbinsX;
  Int_t fNbinsY;
  Int_t fNbinsZ;
  Int_t fNbins;
  Int_t fVarX;
  Int_t fVarY;
  Int_t fVarZ;
  Int_t fTrackVarX;
  Int_t fTrackVarY;
  Int_t fTrackVarZ;
  Int_t fTrackVarMap[3];
  Int_t fNdim;
  Int_t fBin;
  Int_t fBins[7];
  TAxis fBinLimits[3];
  Int_t fNbin[3];
  TString fEventSelection;
  TString fTracking1Quality;
  TString fPID1;
  TString fEventPlanes[Nqvectors];
  TString fXaxisLabel;
  TString fYaxisLabel;
  Int_t fMinHarmonic;
  Int_t fMaxHarmonic;



  Int_t* fMultiplicity[3];
  Float_t* fQvector[3][Nharmonics][2];
  AliCMEAnalysisCuts* fTrackCutsOnePwrtRP;

  THnF* fCorrelationTHn[Nqvectors][3][Nharmonics][2];
  THnF* fCorrelationTHn2[Nqvectors][3][Nharmonics][2];
  THnF* fCorrelationTHnMult[Nqvectors][3];

  THnF* fSumQ[Nqvectors][3][Nharmonics][2];
  THnF* fSumQ2[Nqvectors][3][Nharmonics][2];
  THnF* fSumQMult[Nqvectors][3];


  AliChargeOnePwrtRP(const AliChargeOnePwrtRP &c);
  AliChargeOnePwrtRP& operator= (const AliChargeOnePwrtRP &c);

  ClassDef(AliChargeOnePwrtRP, 1);
};



//_________________________________________________________________
inline void AliChargeOnePwrtRP::FillCorrelation(Int_t bin, Int_t thnbin, AliQnCorrectionsQnVector* qvecEP, Int_t ep) {
  //
  // Correlate Qvectors and fill histogram
  //
  Float_t x,y;
  for(Int_t charge=0; charge<3; charge++){
    if(Multiplicity(bin,charge)==0)continue;


    fCorrelationTHnMult[ep][charge]->AddBinContent(thnbin, Multiplicity(bin,charge));
    fSumQMult[ep][charge]->AddBinContent(thnbin, Multiplicity(bin,charge));
    for(Int_t h=fMinHarmonic; h<=fMaxHarmonic; h++){

      x=fQvector[charge][h-1][0][bin] * qvecEP->QxNorm(h);
      y=fQvector[charge][h-1][1][bin] * qvecEP->QyNorm(h);

      //cout<<bin<<"  "<<thnbin<<"  "<<charge<<"  "<<fCorrelationTHn[ep][charge][h-1][0]->GetNbins()<<endl;

      fCorrelationTHn[ep][charge][h-1][0]->AddBinContent(thnbin, x);
      fCorrelationTHn[ep][charge][h-1][1]->AddBinContent(thnbin, y);
      fCorrelationTHn2[ep][charge][h-1][0]->AddBinContent(thnbin, x*x);
      fCorrelationTHn2[ep][charge][h-1][1]->AddBinContent(thnbin, y*y);

      fSumQ[ep][charge][h-1][0]->AddBinContent(thnbin, fQvector[charge][h-1][0][bin]);
      fSumQ[ep][charge][h-1][1]->AddBinContent(thnbin, fQvector[charge][h-1][1][bin]);
      fSumQ2[ep][charge][h-1][0]->AddBinContent(thnbin, fQvector[charge][h-1][0][bin]*fQvector[charge][h-1][0][bin]);
      fSumQ2[ep][charge][h-1][1]->AddBinContent(thnbin, fQvector[charge][h-1][1][bin]*fQvector[charge][h-1][1][bin]);
    }

  }

  return;
}

//_________________________________________________________________
inline void AliChargeOnePwrtRP::AddTrack( Int_t charge, Float_t x[Nharmonics], Float_t y[Nharmonics]) {
  //
  // Construct Q-vectors for one track
  //

//  fMultiplicity->SetAt(fMultiplicity->At(bin)+1,bin);
//
//    //  1st harmonic
//    fQvector[0][0]->SetAt( fQvector[0][0]->At(bin)+x[0], bin);
//    fQvector[0][1]->SetAt( fQvector[0][1]->At(bin)+y[0], bin);
//    //  2nd harmonic
//    fQvector[1][0]->SetAt( fQvector[1][0]->At(bin)+x[1], bin);
//    fQvector[1][1]->SetAt( fQvector[1][1]->At(bin)+y[1], bin);
//    //  3rd harmonic
//    fQvector[2][0]->SetAt( fQvector[2][0]->At(bin)+x[2], bin);
//    fQvector[2][1]->SetAt( fQvector[2][1]->At(bin)+y[2], bin);
//    //  4th harmonic
//    fQvector[3][0]->SetAt( fQvector[3][0]->At(bin)+x[3], bin);
//    fQvector[3][1]->SetAt( fQvector[3][1]->At(bin)+y[3], bin);
//    //  5th harmonic
//    fQvector[4][0]->SetAt( fQvector[4][0]->At(bin)+x[4], bin);
//    fQvector[4][1]->SetAt( fQvector[4][1]->At(bin)+y[4], bin);
//    //  6th harmonic
//    fQvector[5][0]->SetAt( fQvector[5][0]->At(bin)+x[5], bin);
//    fQvector[5][1]->SetAt( fQvector[5][1]->At(bin)+y[5], bin);



  Int_t bin=Bin();
  fMultiplicity[charge][bin]++;
  fMultiplicity[2][bin]++;

  for(Int_t ih=fMinHarmonic; ih<=fMaxHarmonic; ih++){
    fQvector[charge][ih-1][0][bin] += x[ih-1];
    fQvector[charge][ih-1][1][bin] += y[ih-1];
    fQvector[2][ih-1][0][bin] += x[ih-1];
    fQvector[2][ih-1][1][bin] += y[ih-1];
  }

  //cout<<bin<<"  "<<fMultiplicity[charge][bin]<<"  "<<fQvector[charge][1][0][bin]<<"  "<<fQvector[charge][1][1][bin]<<"  "<<x[1]<<"  "<<y[1]<<endl;

  return;
}


//_________________________________________________________________
inline void AliChargeOnePwrtRP::Clear(){
    for(Int_t g=0; g<3; g++) for(Int_t i=0; i<Nharmonics; i++) for(Int_t j=0; j<2; j++) for(Int_t k=0; k<fNbins; k++) fQvector[g][i][j][k] = 0.0;
    for(Int_t g=0; g<3; g++) for(Int_t k=0; k<fNbins; k++) fMultiplicity[g][k] = 0;
}




#endif
