#ifndef AliFast_H
#define AliFast_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast                                                              //
//                                                                      //
// Main class to control the AliFast program.                           //
//                                                                      //
// This class :                                                         //
//   - Initialises the run default parameters                           //
//   - Provides API to Set/Get run parameters                           //
//   - Creates the support lists (TClonesArrays) for the Event structure//
//   - Creates the physics objects makers                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliRun.h"
//#ifndef ROOT_TTree
#include <TTree.h>
//#endif
//#ifndef AliFHistBrowser_H
#include "AliFHistBrowser.h"
//#endif
//#ifndef AliFBigBang_H
#include "AliFBigBang.h"
//#endif
//#ifndef AliFMaker_H
#include "AliFMaker.h"
//#endif
//#ifndef AliFDet_H
//#include "AliFDet.h"
//#endif

class TBrowser;
class TChain;
class AliFMCMaker;
class AliFTrackMaker;
class AliFDet;
class AliFVirtualDisplay;

class AliFast : public AliRun {

private:
   Int_t               fVersion;           //AliFast version number
   Int_t               fVersionDate;       //AliFast version date
   Int_t               fMode;              //Run mode
   Int_t               fTestTrack;         //Test mode for TrackMaker
   TTree              *fTree;              //!Pointer to the Root tree
   TList              *fMakers;            //List of Makers
//    pointers to standard Makers
   AliFMCMaker        *fMCMaker;           //!Pointer to MCMaker
   AliFTrackMaker     *fTrackMaker;        //!Pointer to TrackMaker
   AliFDet            *fDet;               //!Pointer to Detector
//    flags and switches
   Int_t               fLuminosity;        //Luminosity option (low=1, high=2)
   Bool_t              fBfield;            //B-field (on=1,off=0)
   Bool_t              fSmearing;          //Smearing on=1, off=0
   Int_t               fSUSYcodeLSP;       //Code for SUSY LSP particle
   Int_t               fTrackFinding;      //Track/finding on=1,off=0
   AliFHistBrowser     fHistBrowser;       //Object to Browse Maker histograms
   AliFBigBang         fBigBang;           //!Object to Browse generated particles
   AliFVirtualDisplay *fDisplay;           //!Pointer to Event display object

public:
                      AliFast();
                      AliFast(const char *name, const char *title="The ALIAS fast MonteCarlo");
   virtual           ~AliFast();
   virtual void       Browse(TBrowser *b);
   virtual void       Draw(Option_t *option="");  // *MENU*
   Int_t              GetVersion() {return fVersion;}
   Int_t              GetVersionDate() {return fVersionDate;}
   virtual void       Clear(Option_t *option="");
   AliFVirtualDisplay *Display() {return fDisplay;}
   virtual void       FillClone();
   virtual void       Finish();
   virtual void       GetTreeEvent(Int_t event);  // *MENU*
   virtual void       Init();
   Bool_t             IsFolder() const {return kTRUE;}
   virtual void       Make(Int_t i=0);
   virtual void       Paint(Option_t *option="");
   virtual void       PrintInfo();
   virtual void       SetDefaultParameters();

   TList             *Makers()    {return fMakers;}
   AliFMaker         *Maker(const char *name) {return (AliFMaker*)fMakers->FindObject(name);}
   AliFMCMaker       *MCMaker()       {return fMCMaker;}
   AliFTrackMaker    *TrackMaker()    {return fTrackMaker;}
   AliFDet           *Detector()    {return fDet;}
   TTree             *Tree() {return fTree;}

   Int_t             Run()   {return fRun;}
   Int_t             Event() {return fEvent;}
   Int_t             Mode()  {return fMode;}
   Int_t             TestTrack()  {return fTestTrack;}

//    Getters for flags and switches
   Int_t             Luminosity()   {return fLuminosity;}
   Bool_t            Bfield()       {return fBfield;}
   Bool_t            Smearing()     {return fSmearing;} 
   Int_t             SUSYcodeLSP()  {return fSUSYcodeLSP;}
   Bool_t            TrackFinding() {return fTrackFinding;}


//    Setter for Event Display
   void           SetDisplay(AliFVirtualDisplay *display) {fDisplay=display;}
//    Setters for flags and switches
   void           SetLuminosity(Int_t lumi=1)   {fLuminosity=lumi;}
   void           SetBfield(Bool_t field=1)     {fBfield=field;}
   void           SetSmearing(Bool_t val=1)     {fSmearing=val;} 
   void           SetSUSYcodeLSP(Int_t val=66)  {fSUSYcodeLSP=val;}
   void           SetTrackFinding(Bool_t val=0) {fTrackFinding=val;}

   virtual void   SetRun(Int_t run=1)     {fRun=run;}
   virtual void   SetEvent(Int_t event=1) {fEvent=event;}
   virtual void   SetMode(Int_t mode=0)   {fMode=mode;}
   virtual void   SetTestTrack(Int_t test=0)   {fTestTrack=test;}

   void           FillTree();
   void           InitChain(TChain *chain);
   void           MakeTree(const char* name="T", const char*title="AliFast tree");

   void           SortDown(Int_t n, Float_t *a, Int_t *index, Bool_t down=kTRUE);

   ClassDef(AliFast, 1)   //AliFast control class
};

R__EXTERN AliFast *gAliFast;

#endif

















