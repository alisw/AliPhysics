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
class AliTRDinfoGen : public AliAnalysisTaskSE
{
public:
  enum AliTRDinfoGenSteeringBits{
    kMCdata               = BIT(18)
   ,kUseLocalEvSelection  = BIT(19)
   ,kUseLocalTrkSelection = BIT(20)
   ,kCollision            = BIT(21)
  };

  AliTRDinfoGen();
  AliTRDinfoGen(char* name);
  virtual ~AliTRDinfoGen();
  
  void    ConnectInputData(Option_t *opt) {AliAnalysisTaskSE::ConnectInputData(opt);}
  static Float_t GetEndITS() { return fgkITS;}
  static Float_t GetEndTPC() { return fgkTPC;}
  static Float_t GetEndTRD() { return fgkTRD;}

  Bool_t  HasMCdata() const       { return TestBit(kMCdata);};
  // temporary until check with AliAnalysisTaskSE collision selection mechannism
  Bool_t  IsCollision() const {return TestBit(kCollision);}
  void    SetCollision(Bool_t set=kTRUE) {SetBit(kCollision, set);}

  void    SetLocalEvSelection(AliTRDeventCuts */*cut*/){;} 
  void    SetLocalEvSelection(Bool_t use=kTRUE) {SetBit(kUseLocalEvSelection, use);}
  void    SetLocalTrkSelection(AliESDtrackCuts */*cut*/){;} 
  void    SetLocalTrkSelection(Bool_t use=kTRUE) {SetBit(kUseLocalTrkSelection, use);}
  void    SetLocalV0Selection(AliTRDv0Info *v0);
  void    SetMCdata(Bool_t mc = kTRUE) {SetBit(kMCdata, mc);}
  void    SetTrigger(const Char_t *trigger);

  Bool_t  UseLocalEvSelection() const {return TestBit(kUseLocalEvSelection);}
  Bool_t  UseLocalTrkSelection() const {return TestBit(kUseLocalTrkSelection);}
  void    UserCreateOutputObjects();
  void    UserExec(Option_t *);

private:
  // rough radial limits for TRD
  static const Float_t fgkITS;      // end ITS
  static const Float_t fgkTPC;      // end TPC
  static const Float_t fgkTRD;      // end TRD

  // Trigger selection
  TString              *fEvTrigger; // list of trigger classes separated by space
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
  TTreeSRedirector* DebugStream();

  AliESDEvent      *fESDev;          //! ESD event
  AliMCEvent       *fMCev;           //! MC event
  // event/track cuts OO - to be used
  AliTRDeventCuts  *fEventCut;       // event cut
  AliESDtrackCuts  *fTrackCut;       // track cut
  AliTRDv0Info     *fV0Cut;          // v0 cut
  AliTRDtrackInfo  *fTrackInfo;      //! Track info
  AliTRDeventInfo  *fEventInfo;	     //! Event info
  AliTRDv0Info     *fV0Info;         //! V0 info
  TObjArray        *fTracksBarrel;   //! Array of barrel tracks
  TObjArray        *fTracksSA;       //! Array of stand alone tracks
  TObjArray        *fTracksKink;     //! Array of kink tracks
  TObjArray        *fV0List;         //! V0 container
  TTreeSRedirector *fDebugStream;    //! debug stream

  ClassDef(AliTRDinfoGen, 5)         // entry to TRD analysis train
};
#endif
