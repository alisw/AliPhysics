#ifndef AliToyMCDrawer_H
#define AliToyMCDrawer_H

#include <TGraph2D.h>
#include <TH3F.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>

#include <AliTPCParam.h>
#include <AliTPCROC.h>

#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"

class TPolyMarker3D;
// class AliTPCROC;
class TClonesArray;

/* Visualization class. To use */

/*  AliToyMCDrawer* draw = new AliToyMCDrawer() */
/*  draw->SetFileName("path/to/toyMC.root") */

/*  draw->FillEventArray(Int_t centerEventNumber) */
/*         or                  */
/*  draw->FillEventArray(Double_t time) */
/*    to display with a certain event in the center or at a certain time  */

/*  draw->DrawEvents(Bool_t both, Bool_t before) */
/*    where "both" will display events before and after the middle event and  */
/*    before will show also events before (after) the middle event if true (false)  */
/*    when "both" is false */

class AliToyMCDrawer : public TObject {


 public:
  AliToyMCDrawer();
  AliToyMCDrawer(const AliToyMCDrawer &drawer);
  AliToyMCDrawer& operator = (const AliToyMCDrawer &drawer);

  virtual ~AliToyMCDrawer();
  
  Int_t FillEventArray(Int_t middleEventNbr, Double_t snapShotTime = -1.);
  Int_t FillEventArray(Double_t snapShotTime);
  Int_t GetNumberOfEvents()     const {return fEventArray->GetEntriesFast(); }
  
  void SetFileName(const Char_t* filename) {fFileName = filename;}
  void DrawEvent(AliToyMCEvent *currentEvent, Double_t centerTime, Int_t color);
  void DrawTrack(const AliToyMCTrack *track,  Double_t centerTime, Double_t currentEventTime, Int_t color);
  void DrawTrack2D(const AliToyMCTrack *track,  Double_t centerTime, Double_t currentEventTime, Int_t color);
  void DrawLaserEvent(Int_t nLaserEvents=1, Int_t side=-1, Int_t rod=-1, Int_t bundle=-1, Int_t beam=-1);
  void DrawGeometry();
  void DrawGeometry2D();
  void DrawEvents(Bool_t both = kFALSE, Bool_t before = kTRUE);
  //  void DrawEvents(Bool_t time = kTRUE, Bool_t both = kTRUE, Bool_t before = kTRUE);

  void SetProjectionType(const char* type) { fProjectionType=type; fProjectionType.ToUpper(); }
  void SetRangeTimeZ  (Float_t min, Float_t max) { fTimeZmin=min; fTimeZmax=max; }
  void SetRangeGlobalX(Float_t min, Float_t max) { fGlobalXmin=min; fGlobalXmax=max; }
  void SetRangeGlobalR(Float_t min, Float_t max) { fGlobalXmin=min; fGlobalXmax=max; }
  void SetRangeGlobalY(Float_t min, Float_t max) { fGlobalYmin=min; fGlobalYmax=max; }
  
  const AliToyMCEvent* GetEvent(Int_t eventnr) const {return static_cast<const AliToyMCEvent*>(fEventArray->At(eventnr));}
// private:

  TTree* fInputTree;
  TFile* fInFile;
  TString fFileName;
  AliToyMCEvent* fEvent;
  TClonesArray* fEventArray;
  TH1* fDispHist;
  
  Double_t fCenterTime;
  Double_t fDriftVel;
  AliTPCParam *fTPCParam;
  Double_t fMaxZ0;
  Double_t fIFCRadius;
  Double_t fOFCRadius;
  Double_t fTimeRange;
  AliTPCROC *fRoc;
  TClonesArray *fPoints;
  TClonesArray *fDistPoints;

  TString         fProjectionType;                  // projection type, x,y,z,r combinations
  Float_t         fTimeZmin;                            // Xmin (time axis)
  Float_t         fTimeZmax;                            // Xmax (time axis)
  Float_t         fGlobalXmin;                            // Ymin (global x)
  Float_t         fGlobalXmax;                            // Ymax (global x)
  Float_t         fGlobalYmin;                            // Zmin (global y)
  Float_t         fGlobalYmax;                            // Zmax (global y)
  
  Bool_t ConnectInputTree();
  
  ClassDef(AliToyMCDrawer, 1);

};










#endif
