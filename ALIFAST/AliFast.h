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

#include "AliFBigBang.h"
#include "AliFHistBrowser.h"
#include "AliRun.h"

class TBrowser;
class TChain;
class TTree;

class AliFDet;
class AliFMCMaker;
class AliFMaker;
class AliFTrackMaker;
class AliFVirtualDisplay;

class AliFast : public AliRun {

public:
                      AliFast();
		      AliFast(const AliFast &afast);
                      AliFast(const char *name, const char *title="The ALIAS fast MonteCarlo");
   virtual           ~AliFast();
   AliFast& operator = (const AliFast &afast) 
     {afast.Copy(*this); return (*this);}

   virtual void       Browse(TBrowser *b);
   virtual void       Draw(Option_t *option="");  // *MENU*
   Int_t              GetVersion() const {return fVersion;}
   Int_t              GetVersionDate() const {return fVersionDate;}
   virtual void       Clear(Option_t *option="");
   AliFVirtualDisplay *FDisplay() {return fFDisplay;}
   virtual void       FillClone();
   virtual void       Finish();
   virtual void       GetTreeEvent(Int_t event);  // *MENU*
   virtual void       Init(const char * setup)
     {AliRun::Init(setup);}
   virtual void       Init();
   Bool_t             IsFolder() const {return kTRUE;}
   virtual void       Make(Int_t i=0);
   virtual void       Paint(Option_t *option="");
   virtual void       PrintInfo();
   virtual void       SetDefaultParameters();

   TList             *Makers()    const {return fMakers;}
   AliFMaker         *Maker(const char *name) const {return (AliFMaker*)fMakers->FindObject(name);}
   AliFMCMaker       *MCMaker()     const {return fMCMaker;}
   AliFTrackMaker    *TrackMaker()  const {return fTrackMaker;}
   AliFDet           *Detector()    const {return fDet;}
   TTree             *Tree()        const {return fTree;}

   void              Run(Int_t nevent, const char * setup)
     {AliRun::Run(nevent,setup);}
   Int_t             Run()          const {return fRun;}
   Int_t             Event()        const {return fEvent;}
   Int_t             Mode()         const {return fMode;}
   Int_t             TestTrack()    const {return fTestTrack;}

//    Getters for flags and switches
   Int_t             Luminosity()   const {return fLuminosity;}
   Bool_t            Bfield()       const {return fBfield;}
   Bool_t            Smearing()     const {return fSmearing;} 
   Int_t             SUSYcodeLSP()  const {return fSUSYcodeLSP;}
   Bool_t            TrackFinding() const {return fTrackFinding;}


//    Setter for Event Display
   void           SetDisplay(AliDisplay *display) 
     {fFDisplay=(AliFVirtualDisplay*)display;}
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

   void           SortDown(Int_t n, Float_t *a, Int_t *index, Bool_t down=kTRUE) const;

private:

   void Copy(TObject &afast) const;

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
   AliFVirtualDisplay *fFDisplay;           //!Pointer to Event display object

   ClassDef(AliFast, 1)   //AliFast control class
};

R__EXTERN AliFast *gAliFast;

#endif

















