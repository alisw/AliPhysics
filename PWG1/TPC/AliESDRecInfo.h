#ifndef ALIESDRECINFO_H
#define ALIESDRECINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliRecInfo                                //
//   collect together MC info and Rec info for comparison purposes 
//                                           - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "AliESDkink.h"
class AliESDEvent;
class AliESDtrack;
class AliV0;
class AliESDkink;
class AliESDfriendTrack;
class AliITStrackMI;
class AliTRDtrackV1;
class AliTPCParam;
class AliTPCseed;

/////////////////////////////////////////////////////////////////////////
class AliESDRecInfo: public TObject {
  friend class  AliRecInfoMaker;
  friend class  AliESDRecV0Info;
  friend class  AliESDRecKinkInfo;

public:
  AliESDRecInfo();
  AliESDRecInfo(const AliESDRecInfo& recinfo); 
  AliESDRecInfo& operator=(const AliESDRecInfo& info);
  ~AliESDRecInfo();
  void Update(AliMCInfo* info,AliTPCParam * par, Bool_t reconstructed);
  void UpdateStatus(AliMCInfo* info, Bool_t reconstructed);
  void UpdatePoints(AliESDtrack* track);
  void UpdateTPC(AliMCInfo* info);
  void UpdateITS(AliMCInfo* info);
  void UpdateTOF(AliMCInfo* info);
  //
  void Reset();
  //
  void AddESDtrack(const AliESDtrack *track, AliMCInfo* info);
  void SetESDtrack(const AliESDtrack *track);
  AliESDtrack *GetESDtrack() const { return fESDtrack;}
  AliESDfriendTrack *GetTrackF() const  { return fTrackF;}
  AliTPCseed *GetTPCtrack() const { return fTPCtrack;}
  AliITStrackMI *GetITStrack() const { return fITStrack;}
  AliTRDtrackV1   *GetTRDtrack() const { return fTRDtrack;}
  Int_t      GetStatus(Int_t i) { return fStatus[i];}
protected:
  //
  Float_t  fTPCPoints[10]; //start , biggest end points,max density .. density at the last 30 pad-rows
  Double_t fTPCinR0[5];   //generated position of the track at inner tpc - radius [3] and fi [4]
  Double_t fTPCinR1[5];   //reconstructed postion of the track           - radius [3] and fi [
  Double_t fTPCinP0[5];   //generated position of the track at inner tpc
  Double_t fTPCinP1[5];   //reconstructed postion of the track
  Double_t fTPCAngle0[2]; // generated angle 
  Double_t fTPCAngle1[2]; //refconstructed angle 
  Double_t fTPCDelta[5];  // deltas
  Double_t fTPCPools[5];  // pools
  Double_t fITSinR0[5];   //generated position of the track at inner tpc
  Double_t fITSinR1[5];   //reconstructed postion of the track
  Double_t fITSinP0[5];   //generated position of the track at inner tpc
  Double_t fITSinP1[5];   //reconstructed postion of the track
  Double_t fITSAngle0[2]; // generated angle 
  Double_t fITSAngle1[2]; //refconstructed angle
  Double_t fITSDelta[5];  // deltas
  Double_t fITSPools[5];  // pools
  Float_t  fTRLocalCoord[3];       //local coordinates of the track ref.
  Int_t    fStatus[4];        // status -0 not found - 1 -only in - 2 -in-out -3 -in -out-refit
  Int_t    fLabels[2];         // labels

  Bool_t   fITSOn;           // ITS refitted inward
  Bool_t   fTRDOn;           // ITS refitted inward
  Float_t  fDeltaP;          //delta of momenta
  Double_t fSign;           // sign
  Int_t    fReconstructed;         //flag if track was reconstructed
  Int_t    fFake;             // fake track
  Int_t    fMultiple;         // number of reconstructions
  Bool_t   fTPCOn;           // TPC refitted inward
  Float_t  fBestTOFmatch;        //best matching between times

private:
  AliESDtrack   *fESDtrack;        // esd track
  AliESDfriendTrack *fTrackF;      // friend track
  AliTPCseed *fTPCtrack;        // tpc track
  AliITStrackMI *fITStrack;        // its track
  AliTRDtrackV1   *fTRDtrack;        // trd track
  TClonesArray   *fTracks;         // esd tracks array
  ClassDef(AliESDRecInfo,2)  // container for 
};


#endif
