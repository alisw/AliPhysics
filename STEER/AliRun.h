#ifndef AliRun_H
#define AliRun_H
#include <TROOT.h>
#include <TBrowser.h>
#include <TList.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TGeometry.h>
#include "AliModule.h"
#include "AliHeader.h"
#include "AliMagF.h"
#include "AliMC.h"
#include "AliGenerator.h"
#include "AliLego.h"

enum {kMaxModules = 25, kLenModuleName=7};

class AliDisplay;

class AliRun : public TNamed {

protected:
  Int_t         fRun;          //Current run number
  Int_t         fEvent;        //Current event number (from 1)
  Int_t         fNtrack;       //Number of tracks
  Int_t         fHgwmk;        //Last track purified
  Int_t         fCurrent;      //Last track returned from the stack
  Int_t         fDebug;        //Debug flag
  AliHeader     fHeader;       //Header information
  TTree        *fTreeD;        //Pointer to Tree for Digits
  TTree        *fTreeK;        //Pointer to Tree for Kinematics
  TTree        *fTreeH;        //Pointer to Tree for Hits
  TTree        *fTreeE;        //Pointer to Tree for Header
  TTree        *fTreeR;        //Pointer to Tree for Reconstructed Objects
  TObjArray    *fModules;    //List of Detectors
  TClonesArray *fParticles;    //Pointer to list of particles
  TGeometry    *fGeometry;     //Pointer to geometry
  AliDisplay   *fDisplay;      //Pointer to event display
  TStopwatch    fTimer;        //Timer object
  AliMagF      *fField;        //Magnetic Field Map
  AliMC        *fMC;           //pointer to MonteCarlo object
  char          fDnames[kMaxModules][kLenModuleName];
                               //Array of detector names
  TArrayI      *fImedia;       //Array of correspondence between media and detectors
  Int_t         fNdets;        //Number of detectors
  Float_t       fTrRmax;       //Maximum radius for tracking
  Float_t       fTrZmax;       //Maximu z for tracking
  AliGenerator *fGenerator;    //Generator used in the MC
  Int_t        *fIdtmed;       //Array to contain media numbers
  Bool_t        fInitDone;     //true when initialisation done
  AliLego      *fLego;         //pointer to aliLego object if it exists
  
public:
   // Creators - distructors
   AliRun();
   AliRun(const char *name, const char *title);
   virtual ~AliRun();

   virtual  void  AddHit(Int_t id, Int_t track, Int_t *vol, Float_t *hits) const;
   virtual  void  AddDigit(Int_t id, Int_t *tracks, Int_t *digits) const;
   virtual  void  Browse(TBrowser *b);
   virtual  void  Build();
   virtual  void  BuildSimpleGeometry();
   virtual  void  CleanDetectors();
   virtual  void  CleanParents();
   TObjArray     *Detectors() const {return fModules;}
   TObjArray     *Modules() const {return fModules;}
   Int_t          CurrentTrack() const {return fCurrent;}
   AliDisplay    *Display() { return fDisplay;}
   virtual  Int_t DistancetoPrimitive(Int_t px, Int_t py);
   virtual  void  DumpPart (Int_t i);
   virtual  void  DumpPStack ();
   virtual AliMagF *Field() const {return fField;}
   virtual  void  FillTree();
   virtual  void  FinishPrimary();
   virtual  void  FinishEvent();
   virtual  void  FinishRun();
   virtual  void  FlagTrack(Int_t track);
   Int_t          GetEvNumber() const {return fEvent;}
   Int_t          GetRunNumber() const {return fRun;}
   void           SetRunNumber(Int_t run) {fRun=run;}
   Int_t          GetDebug() const {return fDebug;}
   AliModule     *GetModule(const char *name);
   Int_t          GetModuleID(const char *name);
   virtual  Int_t GetEvent(Int_t event);
   TGeometry     *GetGeometry();
   AliHeader     *GetHeader() {return &fHeader;}
   virtual  void  GetNextTrack(Int_t &mtrack, Int_t &ipart, Float_t *pmom, Float_t &e, Float_t *vpos, Float_t *polar, Float_t &tof);
   Int_t          GetNtrack() {return fNtrack;}
   virtual  Int_t GetPrimary(Int_t track);
   virtual  void  Init(const char *setup="Config.C");
   Bool_t         IsFolder() {return kTRUE;}
   virtual  void  MakeTree(Option_t *option="KH");
   TClonesArray  *Particles() {return fParticles;};
   virtual  void  PurifyKine();
   virtual  Int_t PurifyKine(Int_t lastSavedTrack, Int_t nofTracks);
   virtual  void  Reset(Int_t run, Int_t idevent);
   virtual  void  ResetDigits();
   virtual  void  ResetHits();
   virtual  void  ResetPoints();
   virtual  void  SetTransPar(char *filename="$(ALICE_ROOT)/data/galice.cuts");
   virtual  void  ResetStack() {fCurrent=-1;fHgwmk=0;fNtrack=0;fParticles->Clear();}
   virtual  void  Run(Int_t nevent=1, const char *setup="Config.C");
   virtual  void  RunLego(const char *setup="Config.C",Int_t ntheta=60,Float_t themin=2,Float_t themax=178,
			  Int_t nphi=60,Float_t phimin=0,Float_t phimax=360,Float_t rmin=0,
			  Float_t rmax=570,Float_t zmax=10000);
   virtual  void  SetCurrentTrack(Int_t track);                           
   virtual  void  SetDebug(const Int_t level=1) {fDebug = level;}
   virtual  void  SetDisplay(AliDisplay *display) {fDisplay = display;}
   virtual  void  StepManager(Int_t id) const;
   virtual  void  SetField(Int_t type=2, Int_t version=1, Float_t scale=1, Float_t maxField=10, char*filename="$(ALICE_ROOT)/data/field01.dat");
   virtual  void  SetTrack(Int_t done, Int_t parent, Int_t ipart, 
  			       Float_t *pmom, Float_t *vpos, Float_t *polar, 
                               Float_t tof, const char *mecha, Int_t &ntr,
                               Float_t weight=1);
   virtual  void  KeepTrack(const Int_t);
   virtual  void  MediaTable();
   virtual  Float_t TrackingZmax() const {return fTrZmax;}
   virtual  Float_t TrackingRmax() const {return fTrRmax;}
   virtual  void    TrackingLimits( Float_t rmax=1.e10, Float_t zmax=1.e10) {fTrRmax=rmax; fTrZmax=zmax;}
   virtual  Int_t   DetFromMate(Int_t i) { return (*fImedia)[i];}
   virtual  AliGenerator* Generator() {return fGenerator;}
   virtual  void SetGenerator(AliGenerator *generator);
   virtual  void EnergySummary();
   virtual  Int_t* Idtmed() {return fIdtmed;}

  // Functions from GEOCAD
  //_______________________________________________________________________
  
   virtual void ReadEuclid(const char*, Int_t, const char*);
   virtual void ReadEuclidMedia(const char*, Int_t);

   TTree         *TreeD() {return fTreeD;}
   TTree         *TreeE() {return fTreeE;}
   TTree         *TreeH() {return fTreeH;}
   TTree         *TreeK() {return fTreeK;}
   TTree         *TreeR() {return fTreeR;}

  // --------------------------- commons -------------------------------------

   ClassDef(AliRun,1)      //Supervisor class for all Alice detectors
};
 
EXTERN  AliRun *gAlice;
 
#endif
