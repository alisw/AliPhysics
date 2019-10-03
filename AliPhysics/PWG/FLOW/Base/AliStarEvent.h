/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

/*****************************************************************
  AliStarEvent: Track container for a star track
                                     
  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#ifndef ALISTAREVENT_H
#define ALISTAREVENT_H

class TObjArray;
class AliStarTrack;

class AliStarEvent : public TObject {

 public:
  AliStarEvent();
  AliStarEvent( Int_t n );
  AliStarEvent( const AliStarEvent& track );  
  AliStarEvent& operator=( const AliStarEvent& track );
  virtual ~AliStarEvent();
  virtual void Print( Option_t* option="" ) const;

  Int_t GetRunID() const {return (Int_t)fParams[0];}
  Int_t GetEventNumber() const {return (Int_t)fParams[1];}
  Float_t GetVtxX() const {return fParams[2];}
  Float_t GetVtxY() const {return fParams[3];}
  Float_t GetVtxZ() const {return fParams[4];}
  Float_t GetBField() const {return fParams[5];}
  Int_t GetRefMult() const {return (Int_t)fParams[6];}
  Int_t GetCentralityID() const {return (Int_t)fParams[7];}
  Int_t GetNumberOfPrimaryTracks() const {return (Int_t)fParams[8];}
  Int_t GetNumberOfTracks() const {return (Int_t)fParams[9];}
  const Float_t* GetParams() const {return fParams; }
  Int_t CalculateCentrality(Int_t refMult) const;

  const AliStarTrack* GetTrack( const Int_t i ) const;
  void AddTrack( AliStarTrack* track );
  void Reset();

  void SetRunID( const Int_t p )  { fParams[0]=(Float_t)p;}
  void SetEventNumber( const Int_t p )  { fParams[1]=(Float_t)p;}
  void SetVtxX( const Float_t p )  { fParams[2]=p;}
  void SetVtxY( const Float_t p )  { fParams[3]=p;}
  void SetVtxZ( const Float_t p )  { fParams[4]=p;}
  void SetBField( const Float_t p )  { fParams[5]=p;}
  void SetRefMult( const Int_t p )  { fParams[6]=(Float_t)p;}
  void SetCentralityID( const Int_t p )  { fParams[7]=(Float_t)p;}
  void SetNumberOfPrimaryTracks( const Int_t p )  { fParams[8]=(Float_t)p;}
  void SetNumberOfTracks( const Int_t p )  { fParams[9]=(Float_t)p;}
  void SetParams( const Float_t* params );

 private:
  static const Int_t fgkNparams = 10; //number of params
  Float_t fParams[fgkNparams]; //params
  TObjArray* fTracks; //track collection

  ClassDef(AliStarEvent,1)         // Base class

};

#endif

