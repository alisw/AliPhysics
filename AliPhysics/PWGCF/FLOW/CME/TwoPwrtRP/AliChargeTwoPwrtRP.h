

#ifndef ALICHARGETWOPWRTRP_H
#define ALICHARGETWOPWRTRP_H

#include <TMath.h>
#include <THn.h>
#include <TProfile.h>
#include <AliQnCorrectionsQnVector.h>
#include "AliChargeOnePwrtRP.h"
#include "AliCMEAnalysisCuts.h"


//_________________________________________________________________________
class AliChargeTwoPwrtRP : public TObject {

 public:
  AliChargeTwoPwrtRP();
  ~AliChargeTwoPwrtRP();


  enum Constants {
      N2p=8,
      N3p=3
  };


  static UChar_t fCharge;
  static UShort_t fBins[7];
  static Float_t fPair[N2p+N3p*AliChargeOnePwrtRP::Nqvectors][2];

  static void SetCharge(Int_t ch1, Int_t ch2) {fCharge = ( ch1==1 ? (ch2==1 ? 0 : 2) : ( ch2==1 ? 3 : 1 ));} // 0: ++   1: --   2: +-   3: -+
  static void SetBin(Int_t Ct,Int_t Pt,Int_t mPt,Int_t dPt,Int_t Eta,Int_t mEta,Int_t dEta) {fBins[0]=Ct;fBins[1]=Pt;fBins[2]=mPt;fBins[3]=dPt;fBins[4]=Eta;fBins[5]=mEta;fBins[6]=dEta;}
  //static Float_t* GetCorrelationArray() {return fTwoPCorrelationArray;}

  void SetVars(Int_t varX, Int_t varY=0, Int_t varZ=0) {fVarX=varX;fVarY=varY;fVarZ=varZ;}
  void SetEventSelected(Bool_t sel) {fEventSelected=sel;}
  void SetAxes(TAxis x, TAxis y, TAxis z)  {fBinLimits[0]=TAxis(x);fBinLimits[1]=TAxis(y);fBinLimits[2]=TAxis(z);};


  Bool_t IsEventSelected() {return fEventSelected;}
  //Float_t Get3pValue(Int_t cor, Int_t com, Int_t bin) {return f3pCorrelations[cor][com]->At(bin);}
  //Float_t Get2pValue(Int_t cor, Int_t com, Int_t bin) {return f2pCorrelations[cor][com]->At(bin);}
  //Int_t Multiplicity(Int_t bin) {return fMultiplicity->At(bin);}
  //Float_t Get3pValue(Int_t charge, Int_t cor, Int_t com, Int_t bin) {return f3pCorrelation[charge][cor][com][bin];}
  //Float_t Get2pValue(Int_t charge, Int_t cor, Int_t com, Int_t bin) {return f2pCorrelation[charge][cor][com][bin];}
  Int_t Multiplicity(Int_t bin, Int_t charge) {return fMult[charge][bin];}
  //void GetTwoParticleCorrelation(Int_t bin, Int_t charge, Float_t* x1, Float_t* y1, Float_t* x2, Float_t* y2);
  //void GetTwoParticleCorrelation(Int_t bin, Int_t charge, Float_t c[7][4], Float_t c2[7][4]);
  //void GetTwoParticleCorrelation(Int_t bin, Int_t charge, Float_t c[N2p+N3p][2]);
  void GetTwoParticleCorrelation();
  void GetRPcorrelations(Int_t bin, Int_t charge, AliQnCorrectionsQnVector* QvecPsi);
  void GetRPCorrelations(Int_t bin, Int_t charge, Float_t* Qx, Float_t* Qy, Int_t n, Int_t m, Int_t k, Int_t factor1, Int_t factor2, Int_t cor, Int_t holder);
  //Short_t Charge(Int_t ch1, Int_t ch2) {if(ch1==1) {if(ch2==1) return 0; else return 2;} else {if(ch2==1) return 3; else return 1;}} // 0: ++   1: --   2: +-   3: -+
  //static UShort_t Charge(Int_t ch1, Int_t ch2) {return ( ch1==1 ? (ch2==1 ? 0 : 2) : ( ch2==1 ? 3 : 1 ));} // 0: ++   1: --   2: +-   3: -+
  AliCMEAnalysisCuts* GetTrack1Cuts() {return fTrackCutsTwoPwrtRP[0];}
  AliCMEAnalysisCuts* GetTrack2Cuts() {return fTrackCutsTwoPwrtRP[1];}
  Int_t Bin(Int_t ix, Int_t iy=0, Int_t iz=0) const {return ix*fNbinsY*fNbinsZ+iy*fNbinsZ+iz;}
  Int_t Bin() const {return fBins[fTrackVarX]*fNbinsY*fNbinsZ+fBins[fTrackVarY]*fNbinsZ+fBins[fTrackVarZ];}
  //TString GetVarXName() {return (fVarX==0 ? "" : fVarNames[fVarX]);}
  //TString GetVarYName() {return (fVarY==0 ? "" : fVarNames[fVarY]);}
  //TString GetVarZName() {return (fVarZ==0 ? "" : fVarNames[fVarZ]);}
  Int_t GetVarX() {return fVarX;}
  Int_t GetVarY() {return fVarY;}
  Int_t GetVarZ() {return fVarZ;}
  Int_t GetTrackVarX() {return fTrackVarX;}
  Int_t GetTrackVarY() {return fTrackVarY;}
  Int_t GetTrackVarZ() {return fTrackVarZ;}
  Int_t GetTrackVarMap(Int_t i) {return fTrackVarMap[i];}
  //Int_t GetNdim() {return ( fNbinsZ!=0 ? 3 : ( fNbinsY!=0 ? 2 : 1 )) ;}
  Int_t GetNdim() {return fNdim;}
  Int_t GetNbins() {return fNbins;}
  Int_t GetNbinsX() {return fNbinsX;}
  Int_t GetNbinsY() {return fNbinsY;}
  Int_t GetNbinsZ() {return fNbinsZ;}
  void Clear();
//  //void Fill2pCorrelation(Int_t charge, Int_t har, Int_t bin, Float_t value, Float_t value2) {f2pCorrelationTHn[charge][har][0]->FillBin(bin, value);f2pCorrelationTHn[charge][har][1]->FillBin(bin, value2);}
//  //void Fill3pCorrelation(Int_t charge, Int_t har, Int_t bin, Float_t value, Float_t value2) {f3pCorrelationTHn[charge][har][0]->FillBin(bin, value);f3pCorrelationTHn[charge][har][1]->FillBin(bin, value2);}
  void FillCorrelation(Int_t thnbin, Int_t arraybin) {
    for(Int_t charge=0; charge<4; charge++){
      //cout<<"fillcor  "<<thnbin<<"  "<<arraybin<<"  "<<charge<<"  "<<fMult[charge][arraybin]<<endl;
    if(fMult[charge][arraybin]==0) continue;
    for(Int_t icor=0; icor<N2p; icor++){
      //if(arraybin>(fNbins-2)) cout<<"exceed "<<arraybin<<"  "<<fNbins<<endl;
      //cout<<"fill thn  "<<arraybin<<"  "<<"  "<<thnbin<<"  "<<f2pCorrelationShort[charge][icor][arraybin]<<"  "<<f2pCorrelationTHn[charge][icor][0]->GetNbins()<<endl;
      f2pCorrelationTHn[charge][icor][0]->AddBinContent(thnbin, f2pCorrelationShort[charge][icor][arraybin]);
      f2pCorrelationTHn[charge][icor][1]->AddBinContent(thnbin, f2pCorrelationShort2[charge][icor][arraybin]);
      //f2pCorrelationTHn[charge][icor][0]->SetEntries(f2pCorrelationTHn[charge][icor][0]->GetEntries()+fMult[charge][arraybin]);
      //f2pCorrelationTHn[charge][icor][1]->SetEntries(f2pCorrelationTHn[charge][icor][1]->GetEntries()+fMult[charge][arraybin]);
    }
    for(Int_t icor=0; icor<N3p; icor++){//N3p; icor++){
      for(Int_t ep=0; ep<AliChargeOnePwrtRP::Nqvectors; ep++) {
        //cout<<"add "<<f3pCorrelationShort[charge][icor][ep][arraybin]<<endl;
        f3pCorrelationTHn[charge][icor][ep][0]->AddBinContent(thnbin, f3pCorrelationShort[charge][icor][ep][arraybin]);
        f3pCorrelationTHn[charge][icor][ep][1]->AddBinContent(thnbin, f3pCorrelationShort2[charge][icor][ep][arraybin]);
      }
    }
    fCorrelationTHnMult[charge]->AddBinContent(thnbin,fMult[charge][arraybin]);
    //if(charge==0) cout<<thnbin<<"  "<<f2pCorrelationShort[charge][0][arraybin]<<endl;
//    //f3pCorrelationTHn[charge][0][0][0]->FillBin(bin, f3pCorrelationShort[charge][0][0][bin]);f3pCorrelationTHn[charge][0][0][1]->FillBin(bin, f3pCorrelationShort2[charge][0][0][bin]);
//    //f3pCorrelationTHn[charge][0][1][0]->FillBin(bin, f3pCorrelationShort[charge][0][1][bin]);f3pCorrelationTHn[charge][0][1][1]->FillBin(bin, f3pCorrelationShort2[charge][0][1][bin]);
   }
    //if(fMult[0][arraybin]!=0) cout<<arraybin<<"  "<<f3pCorrelationShort[0][1][0][arraybin]<<"  "<<f2pCorrelationShort[0][0][arraybin]<<endl;
  }

  void FillCorrelationMult(Int_t charge, Int_t bin) {fCorrelationTHnMult[charge]->AddBinContent(bin, 1.);}
  THnF* Get2pCorrelationTHn(Int_t charge, Int_t har, Int_t var) {return f2pCorrelationTHn[charge][har][var];}
  THnF* Get3pCorrelationTHn(Int_t charge, Int_t har, Int_t ep, Int_t var) {return f3pCorrelationTHn[charge][har][ep][var];}
  THnF* GetCorrelationTHnMult(Int_t charge) {return fCorrelationTHnMult[charge];}

  void SetTHn(TString eventSelection, TString tracking1Quality, TString tracking2Quality,TString p1, TString p2, TString d1, TString d2="", Char_t* ar[]=0x0);

  void SetTrackCuts(AliCMEAnalysisCuts* p1,AliCMEAnalysisCuts* p2) {fTrackCutsTwoPwrtRP[0] = p1;fTrackCutsTwoPwrtRP[1] = p2; fTracking1Quality = p1->GetName(); fTracking2Quality = p2->GetName();}

  void Init(TAxis* x, TAxis* y, TAxis* z);




 private:
  Int_t* fMult[4];
  Float_t* f2pCorrelationShort[4][N2p];
  Float_t* f2pCorrelationShort2[4][N2p];
  Float_t* f3pCorrelationShort[4][N3p][AliChargeOnePwrtRP::Nqvectors];
  Float_t* f3pCorrelationShort2[4][N3p][AliChargeOnePwrtRP::Nqvectors];
  Int_t fNbinsX;
  Int_t fNbinsY;
  Int_t fNbinsZ;
  Int_t fNbins;
  Int_t fNdim;
  Int_t fVarX;
  Int_t fVarY;
  Int_t fVarZ;
  Int_t fTrackVarX;
  Int_t fTrackVarY;
  Int_t fTrackVarZ;
  Int_t fTrackVarMap[3];
  //TString fVarNames[6];
  TAxis fBinLimits[3];

  Bool_t fEventSelected;
  THnF* f2pCorrelationTHn[4][N2p][2];
  THnF* f3pCorrelationTHn[4][N3p][AliChargeOnePwrtRP::Nqvectors][2];
  THnF* fCorrelationTHnMult[4];

  AliCMEAnalysisCuts* fTrackCutsTwoPwrtRP[2];

  TString fEventSelection;
  TString fTracking1Quality;
  TString fTracking2Quality;
  TString fPID1;
  TString fPID2;
  TString fEventPlanes[AliChargeOnePwrtRP::Nqvectors];
  TString fXaxisLabel;
  TString fYaxisLabel;

  AliChargeTwoPwrtRP(const AliChargeTwoPwrtRP &c);
  AliChargeTwoPwrtRP& operator= (const AliChargeTwoPwrtRP &c);

  ClassDef(AliChargeTwoPwrtRP, 1);
};





//_________________________________________________________________
inline void AliChargeTwoPwrtRP::GetTwoParticleCorrelation(){
  //
  // Construct two particle correlations
  //


  Int_t bin=Bin();
  fMult[fCharge][bin]++;

  //if(bin<0||bin>=fNbins) cout<<"error!!  "<<bin<<"  "<<fNbins<<"  "<<fTrackVarX<<"  "<< fBins[fTrackVarX]<<"  "<<fNbinsY<<"  "<<fNbinsZ<<"  "<<fBins[fTrackVarY]<<"  "<<fNbinsZ<<"  "<<fBins[fTrackVarZ]<<endl;

  //cout<<bin<<"  "<<(Int_t )fCharge<<"  "<<fMult[fCharge][bin]<<endl;

    //f2pCorrelationShort[fCharge][0][bin]     += fPair[0][0];
    //f2pCorrelationShort2[fCharge][0][bin]    += fPair[0][1];

    //f2pCorrelationShort[fCharge][1][bin]     += fPair[1][0];
    //f2pCorrelationShort2[fCharge][1][bin]    += fPair[1][1];

    //f3pCorrelationShort[fCharge][0][0][bin]  += fPair[2][0];
    //f3pCorrelationShort2[fCharge][0][0][bin] += fPair[2][1];

    //f3pCorrelationShort[fCharge][0][1][bin]  += fPair[3][0];
    //f3pCorrelationShort2[fCharge][0][1][bin] += fPair[3][1];

    //f3pCorrelationShort[fCharge][1][0][bin]  += fPair[4][0];
    //f3pCorrelationShort2[fCharge][1][0][bin] += fPair[4][1];

    //f3pCorrelationShort[fCharge][1][1][bin]  += fPair[5][0];
    //f3pCorrelationShort2[fCharge][1][1][bin] += fPair[5][1];

    //f3pCorrelationShort[fCharge][0][2][bin]  += fPair[6][0];
    //f3pCorrelationShort2[fCharge][0][2][bin] += fPair[6][1];

    //f3pCorrelationShort[fCharge][0][3][bin]  += fPair[7][0];
    //f3pCorrelationShort2[fCharge][0][3][bin] += fPair[7][1];

    //f3pCorrelationShort[fCharge][1][2][bin]  += fPair[8][0];
    //f3pCorrelationShort2[fCharge][1][2][bin] += fPair[8][1];

    //f3pCorrelationShort[fCharge][1][3][bin]  += fPair[9][0];
    //f3pCorrelationShort2[fCharge][1][3][bin] += fPair[9][1];

    for(Int_t i=0; i<N2p; i++){
        f2pCorrelationShort[fCharge][i][bin]  += fPair[i][0];
        f2pCorrelationShort2[fCharge][i][bin] += fPair[i][1];
    }
    for(Int_t i=0; i<N3p; i++){
      for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++){
        f3pCorrelationShort[fCharge][i][j][bin]  += fPair[N2p+i*AliChargeOnePwrtRP::Nqvectors+j][0];
        f3pCorrelationShort2[fCharge][i][j][bin] += fPair[N2p+i*AliChargeOnePwrtRP::Nqvectors+j][1];
    }}




}






//_________________________________________________________________
inline void AliChargeTwoPwrtRP::GetRPcorrelations(Int_t bin, Int_t charge, AliQnCorrectionsQnVector* QvecPsi){

  Float_t QxPsi[6],QyPsi[6];

    for(Int_t ih=0; ih<6; ++ih) {
      //QxPsi[ih] = QvecPsi->Qx(ih+1);//TMath::Cos(QvecPsi->EventPlane(ih+1));
      //QyPsi[ih] = QvecPsi->Qy(ih+1);//TMath::Sin(QvecPsi->EventPlane(ih+1));
      QxPsi[ih] = TMath::Cos(QvecPsi->EventPlane(ih+1));
      QyPsi[ih] = TMath::Sin(QvecPsi->EventPlane(ih+1));
    }

      // Qx,Qy correlations

   GetRPCorrelations(bin, charge, QxPsi, QyPsi, 1, -1, 2, 1, 1, 0, 0);
   //GetRPCorrelations(bin,QxPsi, QyPsi, 2, -2, 4, 2, 2, 1, 1);
   GetRPCorrelations(bin, charge, QxPsi, QyPsi, 1,  3, 2, 1, 1, 4, 2);
   //GetRPCorrelations(bin,QxPsi, QyPsi, 1,  4, 3, 2, 1, 5, 3);
   //GetRPCorrelations(bin,QxPsi, QyPsi, 1,  5, 2, 2, 2, 6, 4);
   //GetRPCorrelations(bin,QxPsi, QyPsi, 1,  5, 4, 2, 2, 6, 5);

}


//_________________________________________________________________
inline void AliChargeTwoPwrtRP::GetRPCorrelations(Int_t bin, Int_t charge, Float_t* Qx, Float_t* Qy, Int_t n, Int_t m, Int_t k, Int_t factor1, Int_t factor2, Int_t cor, Int_t holder){

    Float_t kf[6][2];
    Float_t x = Qx[k-1];
    Float_t y = Qy[k-1];
    kf[0][0] = Qx[k-1];
    kf[0][1] = Qy[k-1];

    //  2nd harmonic
    kf[1][0] = (2.0*x*x-1);
    kf[1][1] = (2.0*x*y);
    //  3rd harmonic
    kf[2][0] = (4.0*x*x*x-3.0*x);
    kf[2][1] = (3.0*y-4.0*y*y*y);
    //  4th harmonic
    kf[3][0] = (1.0-8.0*x*x*y*y);
    kf[3][1] = (4.0*x*y-8.0*x*y*y*y);

    //f3pCorrelations[holder][0]->SetAt(f2pCorrelations[cor][0]->At(bin)*kf[factor1-1][0]*kf[factor2-1][0], bin);
    //f3pCorrelations[holder][1]->SetAt(f2pCorrelations[cor][0]->At(bin)*kf[factor1-1][1]*kf[factor2-1][1], bin);
    //f3pCorrelations[holder][2]->SetAt(f2pCorrelations[cor][1]->At(bin)*kf[factor1-1][0]*kf[factor2-1][1], bin);
    //f3pCorrelations[holder][3]->SetAt(f2pCorrelations[cor][1]->At(bin)*kf[factor1-1][1]*kf[factor2-1][0], bin);
    //f3pCorrelations[holder][4]->SetAt(f2pCorrelations[cor][2]->At(bin)*kf[factor1-1][0]*kf[factor2-1][1], bin);
    //f3pCorrelations[holder][5]->SetAt(f2pCorrelations[cor][2]->At(bin)*kf[factor1-1][1]*kf[factor2-1][0], bin);
    //f3pCorrelations[holder][6]->SetAt(f2pCorrelations[cor][3]->At(bin)*kf[factor1-1][0]*kf[factor2-1][0], bin);
    //f3pCorrelations[holder][7]->SetAt(f2pCorrelations[cor][3]->At(bin)*kf[factor1-1][1]*kf[factor2-1][1], bin);

    //f3pCorrelation[charge][holder][0][bin] = f2pCorrelation[charge][cor][0][bin]*kf[factor1-1][0]*kf[factor2-1][0];
    //f3pCorrelation[charge][holder][1][bin] = f2pCorrelation[charge][cor][0][bin]*kf[factor1-1][1]*kf[factor2-1][1];
    //f3pCorrelation[charge][holder][2][bin] = f2pCorrelation[charge][cor][1][bin]*kf[factor1-1][0]*kf[factor2-1][1];
    //f3pCorrelation[charge][holder][3][bin] = f2pCorrelation[charge][cor][1][bin]*kf[factor1-1][1]*kf[factor2-1][0];
    //f3pCorrelation[charge][holder][4][bin] = f2pCorrelation[charge][cor][2][bin]*kf[factor1-1][0]*kf[factor2-1][1];
    //f3pCorrelation[charge][holder][5][bin] = f2pCorrelation[charge][cor][2][bin]*kf[factor1-1][1]*kf[factor2-1][0];
    //f3pCorrelation[charge][holder][6][bin] = f2pCorrelation[charge][cor][3][bin]*kf[factor1-1][0]*kf[factor2-1][0];
    //f3pCorrelation[charge][holder][7][bin] = f2pCorrelation[charge][cor][3][bin]*kf[factor1-1][1]*kf[factor2-1][1];
    //if(holder==2) cout<<f3pCorrelation[holder][0][bin] <<"  "<< f2pCorrelation[cor][0][bin]<<"   "<<kf[factor1-1][0]<<"  "<<kf[factor2-1][0]<<endl;


}


#endif
