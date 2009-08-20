#ifndef ALITRDRAWSTREAMBASE_H
#define ALITRDRAWSTREAMBASE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDrawStreamBase.h 23387 2008-01-17 17:25:16Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This base class defines access to TRD digits in raw data.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliLog.h"

class AliRawReader;
class AliTRDdigitsManager;

#define TRDMAXTBINS 63
#define TRDMAXADC   21
#define TRDMAXMCM   4 * 16
#define MAXTRACKLETSPERHC 256

#define TRD_NOIMP() AliFatal("Not Implemented for this class. Use inherited.");

class AliTRDrawStreamBase : public TObject
{ // class def begin

 public:

  AliTRDrawStreamBase();
  AliTRDrawStreamBase(AliRawReader *rawReader);
  virtual ~AliTRDrawStreamBase();

  AliTRDrawStreamBase(const AliTRDrawStreamBase& st);
  AliTRDrawStreamBase &operator=(const AliTRDrawStreamBase &);

  //--------------------------------------------------------

  enum STREAMTYPE
    {
      kTRDsimStream  =  0,
      kTRDrealStream =  1,
      kTRDfastStream =  2
    };

  enum { kDDLOffset = 0x400 };                                // Offset for DDL numbers

  static   AliTRDrawStreamBase *GetRawStream();
  static   AliTRDrawStreamBase *GetRawStream(AliRawReader *reader);

  static  void      SetRawStreamVersion(Int_t iver) { fgRawStreamVersion = iver; }
  static  void      SetRawStreamVersion(const char *opt);

  virtual Bool_t    Next() {TRD_NOIMP(); return 0;}          
  //virtual Int_t     NextChamber(AliTRDdigitsManager */*man*/) {TRD_NOIMP(); return 0;} 
  //virtual Int_t     NextChamber(AliTRDdigitsManager */*man*/, UInt_t **/*trackletContainer*/=NULL) {TRD_NOIMP(); return 0;}
  virtual Int_t     NextChamber(AliTRDdigitsManager */*man*/, UInt_t **/*trackletContainer*/=NULL, UShort_t **/*errorCodeContainer*/=NULL) {TRD_NOIMP(); return 0;}
  virtual Bool_t    Init() {TRD_NOIMP(); return -1;}     

  virtual Bool_t    SetRawVersion(Int_t /*fraw*/) {TRD_NOIMP(); return 0;} 
  
  virtual Bool_t    IsCurrentPadShared() const {TRD_NOIMP(); return 0;}
  virtual void      SetSharedPadReadout(Bool_t /*fv*/) {TRD_NOIMP();} 
  virtual Bool_t    IsDataZeroSuppressed() const {TRD_NOIMP(); return 0;}
  
  virtual Bool_t    SetReader(AliRawReader */*reader*/) {TRD_NOIMP(); return 0;} 
    	   
  virtual Bool_t    IsTrackletEnableBitSet() const {TRD_NOIMP(); return 0;}
  virtual Bool_t    IsStackActive(Int_t /*is*/) const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetNofActiveStacks() const {TRD_NOIMP(); return 0;}
  virtual UInt_t   *GetSMstreamPosition() const {TRD_NOIMP(); return 0;}
	    
  virtual Bool_t    IsSMbufferClean() const {TRD_NOIMP(); return 0;}
    
  virtual Bool_t    IsLinkActiveInStack(Int_t /*is*/, Int_t /*il*/) const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetActiveLinksInStack(Int_t /*is*/) const {TRD_NOIMP(); return 0;}
    	    
  virtual Int_t     GetSpecialRawVersion() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetMajorRawVersion() const {TRD_NOIMP(); return 0;}  
  virtual Int_t     GetRawVersion() const {TRD_NOIMP(); return 0;}       
  virtual Int_t     GetMinorRawVersion() const {TRD_NOIMP(); return 0;}  

  virtual Int_t     GetSM() const {TRD_NOIMP(); return 0;}   
  virtual Int_t     GetLayer() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetStack() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetSide() const {TRD_NOIMP(); return 0;} 
  virtual Int_t     GetDCS() const {TRD_NOIMP(); return 0;}  

  virtual Int_t     GetROC() const {TRD_NOIMP(); return 0;}  
  virtual Int_t     GetNumberOfTimeBins() const {TRD_NOIMP(); return 0;} 
  virtual UInt_t    GetBunchCrossCounter() const {TRD_NOIMP(); return 0;}
  virtual UInt_t    GetPreTriggerCounter() const {TRD_NOIMP(); return 0;}
  virtual UInt_t    GetPreTriggerPhase() const {TRD_NOIMP(); return 0;}          

  virtual Int_t     GetRow() const {TRD_NOIMP(); return 0;}   
  virtual Int_t     GetCol() const {TRD_NOIMP(); return 0;}   
  virtual Int_t     GetRowMax() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetColMax() const {TRD_NOIMP(); return 0;}

  // compatibility
  virtual Int_t     GetMaxRow() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetMaxCol() const {TRD_NOIMP(); return 0;}

  virtual Int_t     GetDET() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetDet() const {TRD_NOIMP(); return 0;}
    	    
  virtual Int_t     GetROB() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetMCM() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetEventNumber() const {TRD_NOIMP(); return 0;}
  virtual Int_t     IsMCMcorrupted() const {TRD_NOIMP(); return 0;}

  virtual Int_t    *GetSignals() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetADC() const {TRD_NOIMP(); return 0;}
  virtual Int_t     GetTimeBin() const {TRD_NOIMP(); return 0;}

  virtual Int_t     GetCommonAdditive() const {TRD_NOIMP(); return 0;}
 
  //----------------------------------------------------------
 
 protected:

 private:

  static Int_t fgRawStreamVersion;           // Raw stream version number

  ClassDef(AliTRDrawStreamBase, 0)           // TRD raw stream base class

}; //clas def end

#endif
