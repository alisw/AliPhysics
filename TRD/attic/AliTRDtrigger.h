#ifndef ALITRDTRIGGER_H
#define ALITRDTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrigger.h 29960 2008-11-18 18:08:05Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD trigger class                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class TTree;
class TClonesArray;
class TObjArray;

class AliRunLoader;
class AliRawReader;

class AliTRDmcmTracklet;
class AliTRDgtuTrack;
class AliTRDmcm;
class AliTRDmodule;
class AliTRDdigitsManager;
class AliTRDarrayDictionary;
class AliTRDarrayADC;
class AliTRDgeometry;

class AliTRDtrigger : public TNamed {

 public:  

  enum { kNMCM = 16, kMaxTrackletsPerMCM = 4, kMcmCol = 21 };

  AliTRDtrigger();
  AliTRDtrigger(const Text_t* name, const Text_t* title);
  AliTRDtrigger(const AliTRDtrigger &p);   
  virtual         ~AliTRDtrigger();
  AliTRDtrigger   &operator=(const AliTRDtrigger &p); 

  virtual void     Copy(TObject &p) const;

          void     Init();

          Bool_t   Open(const Char_t *name, Int_t nEvent = 0);
          Bool_t   ReadDigits();
          Bool_t   ReadDigits(AliRawReader *rawReader);
          Bool_t   ReadDigits(TTree *digitsTree);
          Bool_t   MakeTracklets(Bool_t makeTracks = kFALSE);
          void     MakeTracks(Int_t det);
          Bool_t   WriteTracklets(Int_t det);
          Bool_t   ReadTracklets(AliRunLoader *rl);

          void     AddTracklet(Int_t det, Int_t row, Int_t seed, Int_t n);
          void     AddTrack(const AliTRDgtuTrack *t, Int_t det);
          Bool_t   TestTracklet(Int_t det, Int_t row, Int_t seed, Int_t n);
          TObjArray *Tracklets();
          void     ResetTracklets();

          Int_t    GetNumberOfTracks() const;
          Int_t    GetNPrimary() const                           { return fNPrimary;   };
          AliTRDgtuTrack  *GetTrack(Int_t i) const;

          void     SetRunLoader(AliRunLoader *rl)                { fRunLoader = rl;    };
          void     SetMCMcoordinates(Int_t imcm);

 protected:

          Float_t                fField;                       //! Magnetic field
          AliTRDgeometry        *fGeo;                         //! TRD geometry

          AliRunLoader          *fRunLoader;                   //! Run Loader
          AliTRDdigitsManager   *fDigitsManager;               //! TRD digits manager
          TTree                 *fTrackletTree;                //! Tree with tracklets
          TObjArray             *fTracklets;                   //! Array of tracklets

          Int_t                  fNROB;                        //! Number of ROBs in the current chamber
          AliTRDmcm             *fMCM;                         //! Current MCM
          AliTRDmcmTracklet     *fTrk;                         //! Current tracklet
          AliTRDmcmTracklet     *fTrkTest;                     //! Test tracklet
          AliTRDmodule          *fModule;                      //! Current module
          AliTRDgtuTrack        *fGTUtrk;                      //! Current GTU track

          Int_t                  fNtracklets;                  //! Tracklets counter

          AliTRDarrayADC        *fDigits;                      //! Array with digits
          AliTRDarrayDictionary *fTrack0;                      //! Track dictionary 0
          AliTRDarrayDictionary *fTrack1;                      //! Track dictionary 1
          AliTRDarrayDictionary *fTrack2;                      //! Track dictionary 2

          Int_t fNPrimary;                                     //! Number of primary tracks

          TClonesArray          *fTracks;                      //! Array of GTU tracks

  ClassDef(AliTRDtrigger,6)                                    //  TRD trigger class

};

#endif
