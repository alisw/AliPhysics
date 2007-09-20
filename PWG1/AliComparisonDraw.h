
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
#include "AliRecInfo.h"

class AliESDEvent; 
class AliESD;
class AliESDfriend;
class TH1I;

 
class AliComparisonDraw : public TObject {
public :
  AliComparisonDraw(); 
  virtual Bool_t    IsFolder(){return kTRUE;}
  void            InitHisto();
  void            Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC);
  //TH1F            GetPtResol(Float_t pt0, Float_t pt1);
public:
  TH2F* fPtResolLPT;        // pt resolution - low pt
  TH2F* fPtResolHPT;        // pt resolution - high pt 
  TH2F* fPtPoolLPT;         // pt resolution - low pt
  TH2F* fPtPoolHPT;         // pt resolution - high pt 
protected:
   ClassDef(AliComparisonDraw,1);
};















#endif
