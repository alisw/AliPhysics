

#ifndef AliTRDtrackingAnalysis_H
#define AliTRDtrackingAnalysis_H

#include "TObject.h"

class TH1D;
class TH2D;
class TTree;
class TObjArray;
class TGraphErrors;

class AliRunLoader;
class AliTRDgeometry;
class AliESD;
class AliTRDcluster;
class AliTRDtracker;

class AliTRDtrackingAnalysis : public TObject {
  
  const char *fPath;
  
  TObjArray *fRefTPC;
  TObjArray *fRefTRD;
  Int_t fLabels[100000];

  AliRunLoader *fLoader;
  TTree  *fEsdTree;
  AliESD *fESD;

  AliTRDtracker *fTracker;

  // histograms 
  TH1D *fDeltaPt;
  TH1D *fDeltaZ;
  TH1D *fDeltaX;
  TH1D *fDeltaYPos;
  TH1D *fDeltaYNeg;

  TH1D *fNPoints;
  TH1D *fNGood;
  
  TH2D *fRefSpace;

  AliTRDgeometry *fGeo;

  TH1D *fClY2;
  TH1D *fClY3;

  TH1D *fTgPhi;
  TH1D *fClYTgPhi[12];

  TGraphErrors *fGrResTgPhi;
  TGraphErrors *fGrMeanTgPhi;


    //TH1D *fPullY2;
  //TH1D *fPullY3;

  TH1D *fTrklY;
  TH1D *fTrklZ;

  TH1D *fClZ;
  TH2D *fClZZ;
  TH2D *fClYY;
  TH2D *fClYX;
  TH1D *fNLabels;
  TH1D *fBits;
  TH1D *fRefDx;

  TH2D *fClZXref;
  TH2D *fClZXcl;

  TH2D *fClPos;

  void CheckFiles();
  void LoadRecPointsFile();
  void LoadRefs();
  Int_t GetReference(Int_t label); 
  Int_t GetMCPosition(Int_t label, Double_t x, Double_t &Y, Double_t &Z, Double_t &tgphi);

  Int_t GetPhiBin(Double_t phi);
  Double_t GetPhi(Int_t bin);

 public:
  
  AliTRDtrackingAnalysis();
  virtual ~AliTRDtrackingAnalysis() {}
  
  void SetPath(const char *path) {fPath = path;}

  void DrawResolutionPt(int startEvent, int stopEvent);  
  void DrawRecPointResolution(int startEvent, int stopEvent);
  //void DrawTrackletResolution(int startEvent, int stopEvent);
 
  ClassDef(AliTRDtrackingAnalysis,1)            // qa for Digits
};

#endif


