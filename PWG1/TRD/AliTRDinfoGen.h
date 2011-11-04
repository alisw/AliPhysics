#ifndef ALITRDINFOGEN_H
#define ALITRDINFOGEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDinfoGen.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD Performance tender wagon                                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliESDEvent;
class AliMCEvent;
class AliESDfriend;
class AliTRDtrackInfo;
class AliTRDeventInfo;
class AliTRDv0Info;
class AliTRDeventCuts;
class AliESDtrackCuts;
class TObjArray;
class TString;
class TTreeSRedirector;
class AliTRDReconstructor;
class AliTRDgeometry;
class AliTRDinfoGen : public AliAnalysisTaskSE
{
public:
  enum AliTRDinfoGenSteeringBits{
    kMCdata               = BIT(18)
   ,kUseLocalTrkSelection = BIT(19)
   ,kCollision            = BIT(20)
   ,kOCDB                 = BIT(21)
  };
  enum AliTRDinfoGenObjects{
     kTracksESD =  0
    ,kTracksMC
    ,kV0
    ,kTPC
    ,kTRDin
    ,kTRDout
    ,kBarrel
    ,kBarrelMC
    ,kSA
    ,kSAMC
    ,kKink
    ,kKinkMC
    ,kBarrelFriend
    ,kSAFriend
    ,kNObjects
  };
  enum AliTRDinfoGenClasses{
     kStatTrk = 0
    ,kEvType
    ,kBC
    ,kTrigger
    ,kNclasses
  };
  AliTRDinfoGen();
  AliTRDinfoGen(char* name);
  virtual ~AliTRDinfoGen();
  
  static Float_t GetEndITS() { return fgkITS;}
  static Float_t GetEndTPC() { return fgkTPC;}
  static Float_t GetEndTRD() { return fgkTRD;}
  Int_t   GetNRefFigures() const  { return 1;} 
  const char* GetOCDB() const {return fOCDB.Data();}
  Bool_t  GetRefFigure(Int_t ifig);
  Bool_t  Load(const Char_t *fn="AnalysisResults.root", const Char_t *dir="TRD_Performance", const Char_t *name=NULL);

  Bool_t  HasMCdata() const       { return TestBit(kMCdata);};
  // temporary until check with AliAnalysisTaskSE collision selection mechanism
  Bool_t  IsInitOCDB() const {return TestBit(kOCDB);}
  Bool_t  IsCollision() const {return TestBit(kCollision);}
  static const AliTRDReconstructor* Reconstructor() {return fgReconstructor;}
  static AliTRDgeometry*      Geometry() {return fgGeo;}
  void    SetInitOCDB(Bool_t set=kTRUE) {SetBit(kOCDB, set);}
  void    SetCollision(Bool_t set=kTRUE) {SetBit(kCollision, set);}
  //void    SetLocalEvSelection(const AliTRDeventCuts */*cut*/){;} 
  void    SetLocalEvSelection(const AliTRDeventCuts &ec);
  //void    SetLocalTrkSelection(const AliESDtrackCuts */*cut*/){;} 
  void    SetLocalTrkSelection(Bool_t use=kTRUE) {SetBit(kUseLocalTrkSelection, use);}
  void    SetLocalV0Selection(const AliTRDv0Info &v0);
  void    SetMCdata(Bool_t mc = kTRUE) {SetBit(kMCdata, mc);}
  void    SetOCDB(const char *ocdb) {fOCDB=ocdb;}

  Bool_t  UseLocalEvSelection() const {return Bool_t(fEventCut);}
  Bool_t  UseLocalTrkSelection() const {return TestBit(kUseLocalTrkSelection);}
  void    UserCreateOutputObjects();
  void    UserExec(Option_t *);
  void    Terminate(Option_t* option = "");
private:
  // rough radial limits for TRD
  static const Float_t fgkITS;      // end ITS
  static const Float_t fgkTPC;      // end TPC
  static const Float_t fgkTRD;      // end TRD

  // Track selection
  static const Float_t fgkTrkDCAxy; // cm
  static const Float_t fgkTrkDCAz;  // cm
  static const Int_t   fgkNclTPC;   // N clusters TPC
  static const Float_t fgkPt;       // min. pt
  static const Float_t fgkEta;      // eta range
  static AliTRDReconstructor   *fgReconstructor; // single instance of TRD reconstructor used by all tasks
  static AliTRDgeometry        *fgGeo;           // single instance of TRD geometry used by all tasks

  AliTRDinfoGen(const AliTRDinfoGen&);
  AliTRDinfoGen& operator=(const AliTRDinfoGen&);
  TTreeSRedirector* DebugStream();

  AliESDEvent      *fESDev;          //! ESD event
  AliMCEvent       *fMCev;           //! MC event
  // event/track cuts OO - to be used
  AliTRDeventCuts  *fEventCut;       // event cut
  AliESDtrackCuts  *fTrackCut;       // track cut
  AliTRDv0Info     *fV0Cut;          // v0 cut
  TString           fOCDB;           // OCDB location
  AliTRDtrackInfo  *fTrackInfo;      //! Track info
  AliTRDeventInfo  *fEventInfo;	     //! Event info
  AliTRDv0Info     *fV0Info;         //! V0 info
  TObjArray        *fTracksBarrel;   //! Array of barrel tracks
  TObjArray        *fTracksSA;       //! Array of stand alone tracks
  TObjArray        *fTracksKink;     //! Array of kink tracks
  TObjArray        *fV0List;         //! V0 container
  TObjArray        *fContainer;      //! container to store results
  TObjArray        *fRecos;          //! array of reco params
  TTreeSRedirector *fDebugStream;    //! debug stream

  ClassDef(AliTRDinfoGen, 7)         // entry to TRD analysis train
};
#endif
