#ifndef ALIRUN_H
#define ALIRUN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TArrayF.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TMCProcess.h>
#include <TStopwatch.h>
#include <TVirtualMCApplication.h>
#include <TVirtualMC.h>
#include <TError.h>

class TBranch;
class TBrowser;
class TDatabasePDG;
class TFile;
class TGeometry;
class TList;
class TParticle;
class TRandom;
class TTree;

#include "AliRunLoader.h"
class AliDetector;
class AliDisplay;
class AliGenEventHeader;
class AliGenerator;
class AliHeader;
class AliLego;
class AliLegoGenerator;
class AliLegoGenerator;
class AliMCQA;
class AliMagF;
class AliModule;
class AliStack;
class AliTrackReference;


enum {kKeepBit=1, kDaughtersBit=2, kDoneBit=4};


class AliRun : public TVirtualMCApplication {
public:
   // Creators - distructors
   AliRun();
   AliRun(const char *name, const char *title);
   AliRun(const AliRun &arun);
   virtual ~AliRun();

   AliRun& operator = (const AliRun &arun) 
     {arun.Copy(*this); return (*this);}
   virtual  void  AddHit(Int_t id, Int_t track, Int_t *vol, Float_t *hits) const;
   virtual  void  AddDigit(Int_t id, Int_t *tracks, Int_t *digits) const;
   virtual  void  AddHitList(TCollection *hitList) {fHitLists->Add(hitList);}
   virtual  void  Browse(TBrowser *b);
   virtual  void  Build();
   virtual  void  BuildSimpleGeometry();
   virtual  void  CleanDetectors();
   TObjArray     *Detectors() const {return fModules;}
   TObjArray     *Modules() const {return fModules;}
   Int_t          GetCurrentTrackNumber() const;
   AliDisplay    *Display() const { return fDisplay;}
   virtual  Int_t DistancetoPrimitive(Int_t px, Int_t py) const;
   virtual  void  DumpPart (Int_t i) const;
   virtual  void  DumpPStack () const;
   virtual AliMagF *Field() const {return fField;}
   virtual  void  FinishRun();
   virtual  void  FlagTrack(Int_t track);
   void           AddEnergyDeposit(Int_t id, Float_t edep) 
                                       {fEventEnergy[id]+=edep;}
   void           AddModule(AliModule* mod);
   Int_t          GetEvNumber() const;
   Int_t          GetRunNumber() const {return fRun;}
   void           SetRunNumber(Int_t run) {fRun=run;}
   void           SetEventNrInRun(Int_t event) {fEventNrInRun=event;}
   Int_t          GetEventNrInRun() const {return fEventNrInRun;}
   Int_t          GetEventsPerRun() const {return fEventsPerRun;}
   Int_t          GetDebug() const {return fDebug;}
   AliModule     *GetModule(const char *name) const;
   TList*         GetHitLists() const {return fHitLists ;}
   AliDetector   *GetDetector(const char *name) const;
   AliMCQA       *GetMCQA() const {return fMCQA;}
   Int_t          GetModuleID(const char *name) const;
   virtual  const char *GetBaseFile() const 
    {return fBaseFileName.Data();}
   virtual  Int_t GetEvent(Int_t event);
   virtual  void  SetEvent(Int_t event) {fEvent=event;}
   virtual  void  SetConfigFunction(const char * config="Config();");
   virtual  const char *GetConfigFunction() const 
    {return fConfigFunction.Data();}
   TGeometry     *GetGeometry();
   virtual  void  SetGenEventHeader(AliGenEventHeader* header);
   Int_t          GetNtrack() const;
   virtual  Int_t GetPrimary(Int_t track) const;
   virtual  void  Hits2Digits(const char *detector=0); 
   virtual  void  Hits2SDigits(const char *detector=0)   {Tree2Tree("S",detector);}
   virtual  void  SDigits2Digits(const char *detector=0) {Tree2Tree("D",detector);}
   virtual  void  Digits2Reco(const char *detector=0)    {Tree2Tree("R",detector);}
   virtual  void  InitMC(const char *setup="Config.C");
   virtual  void  Init(const char *setup="Config.C") {InitMC(setup);}
   Bool_t         IsFolder() const {return kTRUE;}
   virtual AliLego* Lego() const {return fLego;}

   TObjArray     *Particles() const;
   TParticle     *Particle(Int_t i) const;
   virtual  void  ResetDigits();
   virtual  void  ResetSDigits();
   virtual  void  ResetHits();
// Track reference related 
   virtual void        AddTrackReference(Int_t label);
   TClonesArray   *TrackReferences()   const {return fTrackReferences;}
   virtual void   RemapTrackReferencesIDs(Int_t *map); //remaping track references MI
   virtual void   ResetTrackReferences();
   
   virtual  void  ResetPoints();
   virtual  void  SetTransPar(const char *filename="$(ALICE_ROOT)/data/galice.cuts");
   virtual  void  SetBaseFile(const char *filename="galice.root");
   virtual  void  ReadTransPar();
   virtual  void  RunMC(Int_t nevent=1, const char *setup="Config.C");
   virtual  void  Run(Int_t nevent=1, const char *setup="Config.C") {RunMC(nevent,setup);}
   virtual  void  RunLego(const char *setup="Config.C",Int_t nc1=60,Float_t c1min=2,Float_t c1max=178,
                          Int_t nc2=60,Float_t c2min=0,Float_t c2max=360,Float_t rmin=0,
                          Float_t rmax=430,Float_t zmax=10000, AliLegoGenerator* gener=NULL);
   virtual  Bool_t IsLegoRun() const {return (fLego!=0);}
   virtual  void  RunReco(const char *detector=0, Int_t first = 0, Int_t last = 0);
   virtual  void  SetCurrentTrack(Int_t track);                           
   virtual  void  SetDebug(const Int_t level=0) {fDebug = level;}
   virtual  void  SetDisplay(AliDisplay *display) {fDisplay = display;}
   virtual  void  SetField(Int_t type=2, Int_t version=1, Float_t scale=1, Float_t maxField=10, char*filename="$(ALICE_ROOT)/data/field01.dat");
   virtual  void  SetField(AliMagF* magField);
   virtual  void  PushTrack(Int_t done, Int_t parent, Int_t pdg, 
			   Float_t *pmom, Float_t *vpos, Float_t *polar, 
			   Float_t tof, TMCProcess mech, Int_t &ntr,
			   Float_t weight = 1, Int_t is = 0);
   virtual  void  PushTrack(Int_t done, Int_t parent, Int_t pdg,
			   Double_t px, Double_t py, Double_t pz, Double_t e,
			   Double_t vx, Double_t vy, Double_t vz, Double_t tof,
			   Double_t polx, Double_t poly, Double_t polz,
			   TMCProcess mech, Int_t &ntr, Float_t weight=1,
			   Int_t is = 0);
   virtual  void  SetHighWaterMark(const Int_t nt);
   
   virtual  void  KeepTrack(const Int_t itra);
   virtual  void  MediaTable();
   virtual  void    TrackingLimits( Float_t rmax=1.e10, Float_t zmax=1.e10) {fTrRmax=rmax; fTrZmax=zmax;}
   virtual  Int_t   DetFromMate(Int_t i) const { return (*fImedia)[i];}
   virtual  AliGenerator* Generator() const {return fGenerator;}
   virtual  void SetGenerator(AliGenerator *generator);
   virtual  void ResetGenerator(AliGenerator *generator);
   virtual  void EnergySummary();
   virtual  TDatabasePDG* PDGDB() const {return fPDGDB;}
   
   // MC Application
   //
   virtual  void  ConstructGeometry();
   virtual  void  InitGeometry();     
   virtual  void  GeneratePrimaries();
   virtual  void  BeginEvent();
   virtual  void  BeginPrimary();
   virtual  void  PreTrack();
   virtual  void  Stepping();         
   virtual  void  PostTrack();
   virtual  void  FinishPrimary();
   virtual  void  FinishEvent();
   virtual  Double_t  TrackingZmax() const {return fTrZmax;}
   virtual  Double_t  TrackingRmax() const {return fTrRmax;}
   virtual  void Field(const Double_t* x, Double_t* b) const;

   //
   // End of MC Application

   TTree         *TreeE() {return (fRunLoader)?fRunLoader->TreeE():0x0;}
   TTree         *TreeK() {return (fRunLoader)?fRunLoader->TreeK():0x0;}
   AliStack      *Stack() {return (fRunLoader)?fRunLoader->Stack():0x0;}
   AliHeader*     GetHeader() {return (fRunLoader)?fRunLoader->GetHeader():0x0;}

   TTree         *TreeD() {MayNotUse("TreeD"); return 0x0;}
   TTree         *TreeS() {MayNotUse("TreeS"); return 0x0;}
   TTree         *TreeR() {MayNotUse("TreeR"); return 0x0;}

   
   void SetRunLoader(AliRunLoader* rloader);
   AliRunLoader* GetRunLoader() const {return fRunLoader;}
//   void SetEventFolderName(const char* eventfoldername);
  virtual  void Announce() const;
   
  virtual  void  InitLoaders(); //prepares run (i.e. creates getters)
  static void Deprecated(TObject *obj, const char *method,
			 const char *replacement) {
    if (obj)
      ::Warning(Form("%s::%s", obj->ClassName(), method),
		"method is depricated\nPlease use: %s", replacement);
    else
      ::Warning(method, "method is depricated\nPlease use: %s", replacement);
  }
protected:
  virtual  void  Tree2Tree(Option_t *option, const char *detector=0);
  Int_t          fRun;               //! Current run number
  Int_t          fEvent;             //! Current event number (from 1)
  Int_t          fEventNrInRun;      //! Current unique event number in run
  Int_t          fEventsPerRun;      //  Number of events per run
  Int_t          fDebug;             //  Debug flag
  TObjArray     *fModules;           //  List of Detectors
  TGeometry     *fGeometry;          //  Pointer to geometry
  AliDisplay    *fDisplay;           //! Pointer to event display
  TStopwatch     fTimer;             //  Timer object
  AliMagF       *fField;             //  Magnetic Field Map
  TVirtualMC    *fMC;                //! Pointer to MonteCarlo object
  TArrayI       *fImedia;            //! Array of correspondence between media and detectors
  Int_t          fNdets;             //  Number of detectors
  Float_t        fTrRmax;            //  Maximum radius for tracking
  Float_t        fTrZmax;            //  Maximu z for tracking
  AliGenerator  *fGenerator;         //  Generator used in the MC
  Bool_t         fInitDone;          //! True when initialisation done
  AliLego       *fLego;              //! Pointer to aliLego object if it exists
  TDatabasePDG  *fPDGDB;             //  Particle factory object
  TList         *fHitLists;          //! Lists of hits to be remapped by PurifyKine
  TArrayF        fEventEnergy;       //! Energy deposit for current event
  TArrayF        fSummEnergy;        //! Energy per event in each volume
  TArrayF        fSum2Energy;        //! Energy squared per event in each volume
  TString        fConfigFunction;    //  Configuration file to be executed
  TRandom       *fRandom;            //  Pointer to the random number generator
  AliMCQA       *fMCQA;              //  Pointer to MC Quality assurance class
  TString        fTransParName;      //  Name of the transport parameters file
  TString        fBaseFileName;      //  Name of the base root file

  TClonesArray *fTrackReferences;     //!list of track references - for one primary track only -MI
  AliRunLoader  *fRunLoader;         //!run getter - written as a separate object
private:
  void Copy(AliRun &arun) const;

  ClassDef(AliRun,8)      //Supervisor class for all Alice detectors
};
 
R__EXTERN  AliRun *gAlice;
 
#endif
