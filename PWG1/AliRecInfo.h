#ifndef ALIRECINFO_H
#define ALIRECINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliRecInfo                                //
//   collect together MC info and Rec info for comparison purposes 
//                                           - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliV0.h"
#include "AliESDkink.h"
#include "AliESDfriendTrack.h"
#include "AliITStrackMI.h"
#include "AliTRDtrack.h"
class AliTPCseed;

/////////////////////////////////////////////////////////////////////////
class AliESDRecInfo: public TObject {
  friend class  AliRecInfoMaker;
  friend class  AliESDRecV0Info;
  friend class  AliESDRecKinkInfo;

public:
  AliESDRecInfo();
  AliESDRecInfo(const AliESDRecInfo& recinfo);
  ~AliESDRecInfo();
  void UpdatePoints(AliESDtrack* track);
  void Update(AliMCInfo* info,AliTPCParam * par, Bool_t reconstructed);
  void Reset();
  //
  void SetESDtrack(const AliESDtrack *track);
  AliESDtrack *GetESDtrack() const { return fESDtrack;}
  AliESDfriendTrack *GetTrackF() const  { return fTrackF;}
  AliTPCseed *GetTPCtrack() const { return fTPCtrack;}
  AliITStrackMI *GetITStrack() const { return fITStrack;}
  AliTRDtrack   *GetTRDtrack() const { return fTRDtrack;}
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
  AliTRDtrack   *fTRDtrack;        // trd track
  
  ClassDef(AliESDRecInfo,2)  // container for 
};



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


class AliESDRecV0Info: public TObject {
  friend class  AliRecInfoMaker;
public:
  void Update(Float_t vertex[3]);
protected:
  AliESDRecInfo  fT1;      //track1
  AliESDRecInfo  fT2;      //track2  
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  Double_t       fInvMass;  //reconstructed invariant mass -
  //
  Double_t       fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t       fXr[3];     //rec. position according helix
  //
  Double_t       fRs[2];     // minimum radius in rphi intersection
  Double_t       fDistMinR; // distance at minimal radius
  Double_t       fPm[3];    //momentum at the vertex mother
  Double_t       fAngle[3]; //three angles
  Double_t       fRr;       // rec position of the vertex 
  Int_t          fLab[2];   //MC label of the partecle
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  Int_t          fV0Status;       // status of the kink
  AliV0*         fV0tpc;           // Vo information from reconsturction according TPC
  AliV0*         fV0its;           // Vo information from reconsturction according ITS
  AliV0*         fV0rec;           // V0 information form the reconstruction
  Int_t          fMultiple;     // how man times V0 was recostructed 
  Int_t          fV0Multiple;   // how man times was V0 reconstucted
  Int_t          fRecStatus;    // status form the reconstuction
  ClassDef(AliESDRecV0Info,2)   // container for  
};



class AliESDRecKinkInfo: public TObject {
friend class  AliRecInfoMaker;
public:
  void Update();
protected:
  AliESDRecInfo  fT1;      //track1
  AliESDRecInfo  fT2;      //track2  
  AliESDkink     fKink;    //kink
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  Double_t       fInvMass;  //reconstructed invariant mass -
  //
  Double_t       fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t       fXr[3];     //rec. position according helix
  //
  Double_t       fPm[3];    //momentum at the vertex mother
  Double_t       fAngle[3]; //three angles
  Double_t       fRr;       // rec position of the vertex 
  Double_t       fMinR;     // minimum radius in rphi intersection
  Double_t       fDistMinR; // distance at minimal radius
  Int_t          fLab[2];   //MC label of the partecle
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  Int_t          fStatus;       //status -tracks 
  Int_t          fRecStatus;    //kink -status- 0 - not found  1-good -  fake
  Int_t          fMultiple;     // how many times was kink reconstructed
  Int_t          fKinkMultiple; // how many times was kink reconstructed
  ClassDef(AliESDRecKinkInfo,1)   // container for  
};

#endif
