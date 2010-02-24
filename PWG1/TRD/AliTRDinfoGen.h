#ifndef ALITRDINFOGEN_H
#define ALITRDINFOGEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDinfoGen.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

#ifndef ROOT_TString
#include "TString.h"
#endif

class AliESDEvent;
class AliMCEvent;
class AliESDfriend;
class AliTRDtrackInfo;
class AliTRDeventInfo;
class AliTRDv0Info;
class TObjArray;
class TTreeSRedirector;
class AliTRDinfoGen : public AliTRDrecoTask{
public:
  enum AliTRDinfoGenSteeringBits{
    kUseLocalEvSelection  = BIT(21)
   ,kUseLocalTrkSelection = BIT(22)
   ,kCollision            = BIT(23)
  };
  AliTRDinfoGen();
  virtual ~AliTRDinfoGen();
  
  void    ConnectInputData(Option_t *);
  void    CreateOutputObjects();
  static Float_t GetTPCx() { return fgkTPC;}
  static Float_t GetTOFx() { return fgkTOF;}

  void    Exec(Option_t *);
  Bool_t  IsCollision() const {return TestBit(kCollision);}

  void    SetCollision(Bool_t set=kTRUE) {SetBit(kCollision, set);}
  void    SetLocalEvSelection(Bool_t use=kTRUE) {SetBit(kUseLocalEvSelection, use);}
  void    SetLocalTrkSelection(Bool_t use=kTRUE) {SetBit(kUseLocalTrkSelection, use);}
  void    SetTrigger(const Char_t *trigger) {fEvTrigger = trigger;}

  Bool_t  UseLocalEvSelection() const {return TestBit(kUseLocalEvSelection);}
  Bool_t  UseLocalTrkSelection() const {return TestBit(kUseLocalTrkSelection);}

private:
  // rough radial limits for TRD
  static const Float_t fgkTPC;      // end TPC
  static const Float_t fgkTOF;      // begin TOF
  // Trigger selection
  TString              fEvTrigger;  // list of trigger classes separated by space
  // Vertex selection
  static const Float_t fgkEvVertexZ;// cm
  static const Int_t   fgkEvVertexN;// cm
  // Track selection
  static const Float_t fgkTrkDCAxy; // cm
  static const Float_t fgkTrkDCAz;  // cm
  static const Int_t   fgkNclTPC;   // N clusters TPC
  static const Float_t fgkPt;       // min. pt
  static const Float_t fgkEta;      // eta range

  AliTRDinfoGen(const AliTRDinfoGen&);
  AliTRDinfoGen& operator=(const AliTRDinfoGen&);

  AliESDEvent      *fESDev;          //! ESD event
  AliMCEvent       *fMCev;           //! MC event
  AliESDfriend     *fESDfriend;      //! ESD friends
  AliTRDtrackInfo  *fTrackInfo;      //! Track info
  AliTRDeventInfo  *fEventInfo;		   //! Event info
  TObjArray        *fV0container;    //! V0 container
  AliTRDv0Info     *fV0Info;		     //! V0 info

  ClassDef(AliTRDinfoGen, 2)         // entry to TRD analysis train
};
#endif
