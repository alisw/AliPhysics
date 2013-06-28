#ifndef AliToyMCDrawer_H
#define AliToyMCDrawer_H

#include <TGraph2D.h>
#include <TH3F.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <AliTPCParam.h>
#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"
class TPolyMarker3D;
class AliTPCROC;
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
  void DrawGeometry();
  void DrawEvents(Bool_t both = kFALSE, Bool_t before = kTRUE);
  //  void DrawEvents(Bool_t time = kTRUE, Bool_t both = kTRUE, Bool_t before = kTRUE);
 
  const AliToyMCEvent* GetEvent(Int_t eventnr) const {return static_cast<const AliToyMCEvent*>(fEventArray->At(eventnr));}
 private:

  TTree* fInputTree;
  TFile* inFile;
  const char* fFileName;
  AliToyMCEvent* fEvent;
  TClonesArray* fEventArray;
  TH3F* fDispHist;
  
  Double_t fCenterTime;
  Double_t fDriftVel;
  AliTPCParam *fTPCParam;
  Double_t fMaxZ0;
  Double_t fOFCRadius;
  Double_t fIFCRadius;
  Double_t fTimeRange;
  AliTPCROC *fRoc;
  TClonesArray *fPoints;
  TClonesArray *fDistPoints; 
  ClassDef(AliToyMCDrawer, 1);

};










#endif
