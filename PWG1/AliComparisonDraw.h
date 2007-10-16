
#ifndef AliComparisonDraw_h
#define AliComparisonDraw_h

#include <iostream>
#include <fstream>
using namespace std;
#include <TSelector.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "AliGenInfo.h"
#include "AliESDRecInfo.h"
#include "AliESDRecV0Info.h"
#include "AliESDRecKinkInfo.h"

class AliESDEvent; 
class AliESD;
class AliESDfriend;
class TH1I;
class TH3F;
class TH3;
class TProfile;
class TProfile2D;
class TGraph2D;
class TGraph; 

class AliComparisonDraw : public TObject {
public :
  AliComparisonDraw(); 
  virtual Bool_t    IsFolder(){return kTRUE;}
  void      InitHisto();
  void      Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC);
  
  //
  //
  void      ProcessEff(AliMCInfo* infoMC, AliESDRecInfo *infoRC);
  void      ProcessResolConstrained(AliMCInfo* infoMC, AliESDRecInfo *infoRC);
  void      ProcessTPCdedx(AliMCInfo* infoMC, AliESDRecInfo *infoRC);
  void      ProcessDCA(AliMCInfo* infoMC, AliESDRecInfo *infoRC);

  void MakePlots();


  //TH1F            GetPtResol(Float_t pt0, Float_t pt1);
  static TH1F*       MakeResol(TH2F * his, Int_t integ, Bool_t type); 
  static TGraph2D *  MakeStat2D(TH3 * his, Int_t delta0, Int_t delta1, Int_t type);
  static TGraph *  MakeStat1D(TH3 * his, Int_t delta1, Int_t type);


public:
  //
  // efficiency 
  //
  static    Bool_t    fBDraw;         //option draw temporary results
  TProfile* fEffTPCPt;      // TPC efficiency as function of Pt (tan+-1)
  TProfile* fEffTPCPtMC;    // MC -TPC efficiency as function of Pt (tan+-1)
  TProfile* fEffTPCPtF;     // efficiency for findable tracks
  //
  TProfile* fEffTPCTan;   // TPC efficiency as function of Tan (pt>0.15
  TProfile* fEffTPCTanMC; // MC -TPC efficiency as function of Tan (pt>0.15)
  TProfile* fEffTPCTanF;  // efficiency for findable tracks Tan (pt>0.15)
  //
  TProfile2D* fEffTPCPtTan;    // TPC efficiency as function of Pt and tan
  TProfile2D* fEffTPCPtTanMC;  // MC -TPC efficiency as function of Pt and tan
  TProfile2D* fEffTPCPtTanF;  // TPC efficiency as function of Pt and tan
  //
  // dEdx resolution
  //
  TH2F* fTPCSignalNormTan; // tpc signal normalized to the mean signal - MC
  TH2F* fTPCSignalNormSPhi;   // tpc signal normalized to the mean signal - MC
  TH2F* fTPCSignalNormTPhi;   // tpc signal normalized to the mean signal - MC
  //
  TH3F* fTPCSignalNormTanSPhi;   // tpc signal normalized to the mean signal - MC
  TH3F* fTPCSignalNormTanTPhi;   // tpc signal normalized to the mean signal - MC
  TH3F* fTPCSignalNormTanSPt;   // tpc signal normalized to the mean signal - MC


  //
  //
  TH2F* fPtResolLPT;        // pt resolution - low pt
  TH2F* fPtResolHPT;        // pt resolution - high pt 
  TH2F* fPtPullLPT;         // pt resolution - low pt
  TH2F* fPtPullHPT;         // pt resolution - high pt 
  //
  // Resolution constrained param
  //
  TH2F   *fCPhiResolTan;   // angular resolution -  constrained
  TH2F   *fCTanResolTan;   // angular resolution -  constrained
  TH2F   *fCPtResolTan;    // pt resolution      -  constrained
  TH2F   *fCPhiPullTan;   // angular resolution -  constrained
  TH2F   *fCTanPullTan;   // angular resolution -  constrained
  TH2F   *fCPtPullTan;    // pt resolution      -  constrained
  //
  // DCA resolution
  //
  TH3F  *fD0TanSPtB1;   // distance to vertex y  
  TH3F  *fD1TanSPtB1;   // distance to vertex z  
  TH3F  *fD0TanSPtL1;   // distance to vertex y  
  TH3F  *fD1TanSPtL1;   // distance to vertex z  

protected:
   ClassDef(AliComparisonDraw,1);
};















#endif
