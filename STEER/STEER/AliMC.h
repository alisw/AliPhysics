#ifndef ALIMC_H
#define ALIMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This is the ALICE implementation of TVirtualMCApplication
// for simulation with different implementations
// of the Virtual MonteCarlo
//

#include <TArrayF.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TMCProcess.h>
#include <TVirtualMCApplication.h>

class TParticle;
class TFile;
class TTree;
class TVirtualMC;

class AliGenerator;
class AliTrackReference;
class AliTransportMonitor;

class AliMC : public TVirtualMCApplication {
public:
    AliMC();
    AliMC(const char *name, const char *title);
    virtual ~AliMC();
    
//
//  MC Application
//
   virtual  void  ConstructGeometry();
   virtual  Bool_t  MisalignGeometry();
   virtual  void  ConstructOpGeometry();
   virtual  void  AddParticles();
   virtual  void  InitGeometry();     
   virtual  void  SetAllAlignableVolumes();     
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
   virtual  Int_t   DetFromMate(Int_t i) const { return (*fImedia)[i];}
//

   virtual  AliGenerator* Generator() const {return fGenerator;}
   virtual  void SetGenerator(AliGenerator *generator);
   virtual  void ResetGenerator(AliGenerator *generator);

//
   virtual  void ReadTransPar();
   virtual  void MediaTable();
   virtual  void EnergySummary();
   virtual  void FinishRun();
   void          AddEnergyDeposit(Int_t id, Float_t edep) 
                                       {fEventEnergy[id]+=edep;}
   virtual  void  ResetHits();
   virtual  void  ResetDigits();
   virtual  void  ResetSDigits();
   virtual  void  TrackingLimits( Float_t rmax=1.e10, Float_t zmax=1.e10)
       {fTrRmax=rmax; fTrZmax=zmax;}
   virtual  void  DecayLimits( Float_t rmin = -1., Float_t rmax = -1., Int_t pdg = 0)
       {fRDecayMin = rmin; fRDecayMax = rmax; fDecayPdg = pdg;}
   
   virtual  void  Init();
   virtual  void  SetTransPar(const char *filename="$(ALICE_ROOT)/data/galice.cuts");
   //PH
   virtual  void  AddHit(Int_t id, Int_t track, Int_t *vol, Float_t *hits) const;
   virtual  void  AddDigit(Int_t id, Int_t *tracks, Int_t *digits) const;
   virtual  void  AddHitList(TCollection *hitList) {if (hitList) fHitLists->Add(hitList);}
   Int_t          GetCurrentTrackNumber() const;
   virtual  void  DumpPart (Int_t i) const;
   virtual  void  DumpPStack () const;
   TList*         GetHitLists() const {return fHitLists ;}
   Int_t          GetNtrack() const;
   virtual  Int_t GetPrimary(Int_t track) const;
   const TObjArray* Particles() const;
   TParticle     *Particle(Int_t i) const;
   virtual  void  PushTrack(Int_t done, Int_t parent, Int_t pdg, 
			   const Float_t *pmom, const Float_t *vpos, const Float_t *polar, 
			   Float_t tof, TMCProcess mech, Int_t &ntr,
			   Float_t weight = 1, Int_t is = 0) const;
   virtual  void  PushTrack(Int_t done, Int_t parent, Int_t pdg,
			   Double_t px, Double_t py, Double_t pz, Double_t e,
			   Double_t vx, Double_t vy, Double_t vz, Double_t tof,
			   Double_t polx, Double_t poly, Double_t polz,
			   TMCProcess mech, Int_t &ntr, Float_t weight=1,
			   Int_t is = 0) const;
   virtual  void  SetHighWaterMark(Int_t nt) const;
   
   virtual  void  KeepTrack(Int_t itra) const;
   virtual  void  FlagTrack(Int_t track) const;
   virtual  void  SetCurrentTrack(Int_t track) const;                           
// Track reference related
   AliTrackReference* AddTrackReference(Int_t label, Int_t id = -999);
   // const TClonesArray*      TrackReferences()   const {return &fTrackReferences;}
   virtual void       RemapTrackReferencesIDs(const Int_t *map); //remaping track references MI
   virtual void       ResetTrackReferences();
   virtual void       FixParticleDecaytime();
// Geometry related
   void           SetGeometryFromCDB();
   Bool_t         IsGeometryFromCDB() const;
// Monitor transport   
   void           SetUseMonitoring(Bool_t flag=kTRUE)      { fUseMonitoring = flag; }
   AliTransportMonitor *GetTransportMonitor() const        { return fMonitor; }
// Random number generator status
   void           SetSaveRndmStatus(Bool_t value)          { fSaveRndmStatus = value; }  
   void           SetSaveRndmStatusPerEvent(Bool_t value)  { fSaveRndmEventStatus = value; }  
   void           SetReadRndmStatus(Bool_t value)          { fReadRndmStatus = value; };   
   void           SetRndmFileName(const TString& fileName) { fRndmFileName = fileName; };   

 private:
   AliMC(const AliMC&); // Not implemented
   AliMC& operator= (const AliMC&); // Not implemented
   void MakeTmpTrackRefsTree();
   void ReorderAndExpandTreeTR();
   void CacheVMCInstance(); // Cache pointer to VMC object (to avoid _tls_ penalties)

 private:
   void RemapHits();
   TVirtualMC    *fMC;                //! cached MC instance (to avoid calls to thread_local gMC=TVirtualMC::GetMC())
   AliGenerator  *fGenerator;         //  Generator used in the MC
   Bool_t         fSaveRndmStatus;    //! Options to save random engine status
   Bool_t         fSaveRndmEventStatus; //! Options to save random engine status for each event
   Bool_t         fReadRndmStatus;    //! Options to read random engine status
   Bool_t         fUseMonitoring;     //! Activate monitoring
   TString        fRndmFileName;      //! The file name of random engine status to be read in
   TArrayF        fEventEnergy;       //! Energy deposit for current event
   TArrayF        fSummEnergy;        //! Energy per event in each volume
   TArrayF        fSum2Energy;        //! Energy squared per event in each volume
   Float_t        fTrRmax;            //  Maximum radius for tracking
   Float_t        fTrZmax;            //  Maximu z for tracking
   Float_t        fRDecayMax;         //  Maximum radius for decay
   Float_t        fRDecayMin;         //  Minimum radius for decay
   Int_t          fDecayPdg;          //  PDG code of particle with forced decay length
   TArrayI       *fImedia;            //! Array of correspondence between media and detectors
   TString        fTransParName;      //  Name of the transport parameters file
   AliTransportMonitor *fMonitor;     //! Transport monitoring tool
   TList         *fHitLists;          //! Lists of hits to be remapped by PurifyKine
   //Temporary Track Reference tree related
   TTree         *fTmpTreeTR;            //! Temporary track reference tree
   TFile         *fTmpFileTR;            //! Temporary track reference file
   TClonesArray   fTrackReferences;      //! List of track references - for one primary track only
   TClonesArray   fTmpTrackReferences;   //! Temporary list of track references - for one primary track only
   ClassDef(AliMC, 4)
};

 
#endif
