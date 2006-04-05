#ifndef ALITRDTRIGGER_H
#define ALITRDTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD trigger class                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliTRDgtuTrack.h"

class TTree;

class AliTRDdigitsManager;
class AliTRDdataArrayI;

class AliRunLoader;
class AliRawReader;

class AliTRDmcmTracklet;
class AliTRDmcm;
class AliTRDmodule;
class AliTRDtrigParam;

class AliTRDtrigger : public TNamed {

 public:  

  enum { kNMCM = 16, kMaxTrackletsPerMCM = 4, kMcmCol = 21 };

  AliTRDtrigger();
  AliTRDtrigger(const Text_t* name, const Text_t* title);
  AliTRDtrigger(const AliTRDtrigger &p);   
  virtual ~AliTRDtrigger();

  AliTRDtrigger &operator=(const AliTRDtrigger &p); 
  virtual void Copy(TObject &p) const;

  void Init();

  void   SetRunLoader(AliRunLoader *rl) { fRunLoader = rl; };
  Bool_t Open(const Char_t *name, Int_t nEvent = 0);
  Bool_t ReadDigits();
  Bool_t ReadDigits(AliRawReader* rawReader);
  Bool_t MakeTracklets(Bool_t makeTracks = kFALSE);
  Bool_t WriteTracklets(Int_t det);
  Bool_t ReadTracklets(AliRunLoader *rl);

  void   AddTracklet(Int_t det, Int_t row, Int_t seed, Int_t n);
  TObjArray *Tracklets() { 
    if(!fTracklets) fTracklets = new TObjArray(400); return fTracklets; 
  };
  void   ResetTracklets() { if(fTracklets) fTracklets->Delete(); };
  void   SetMCMcoordinates(Int_t imcm);
  void   SetParameter(AliTRDtrigParam *trigp) { fTrigParam = trigp; };
  AliTRDtrigParam *GetParameter() const { return fTrigParam; };

  void   MakeTracks(Int_t det);

  AliTRDgtuTrack *GetTrack(Int_t i) const {
    return (AliTRDgtuTrack *)fTracks.UncheckedAt(i);
  }
  void AddTrack(const AliTRDgtuTrack *t, Int_t det) {
    AliTRDgtuTrack * track = new(fTracks[fTracks.GetEntriesFast()]) AliTRDgtuTrack(*t);
    track->SetDetector(det);
  }
  Int_t GetNumberOfTracks() const {return fTracks.GetEntriesFast();}

  Int_t GetNPrimary() const { return fNPrimary; };

 protected:

  AliTRDtrigParam       *fTrigParam;                   //! Trigger class parameters
  AliRunLoader          *fRunLoader;                   //! Run Loader
  AliTRDdigitsManager   *fDigitsManager;               //! TRD digits manager
  TTree                 *fTrackletTree;                //! Tree with tracklets
  TObjArray             *fTracklets;                   //! Array of tracklets

  Int_t                  fNROB;                        //! Number of ROBs in the current chamber
  AliTRDmcm             *fMCM;                         //! Current MCM
  AliTRDmcmTracklet     *fTrk;                         //! Current tracklet
  AliTRDmodule          *fModule;                      //! Current module
  AliTRDgtuTrack        *fGTUtrk;                      //! Current GTU track

  Int_t                  fNtracklets;                  //! Tracklets counter

  AliTRDdataArrayI *fDigits;                           //! Array with digits
  AliTRDdataArrayI *fTrack0;                           //! Track dictionary 0
  AliTRDdataArrayI *fTrack1;                           //! Track dictionary 1
  AliTRDdataArrayI *fTrack2;                           //! Track dictionary 2

  Int_t fNPrimary;                                     //! Number of primary tracks

  TClonesArray           fTracks;                      //! Array of GTU tracks

  ClassDef(AliTRDtrigger,2)                            //  TRD trigger class

};

#endif
