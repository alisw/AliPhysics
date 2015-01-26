/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

// $Id: AliJTrack.h,v 1.3 2008/01/21 11:56:39 djkim Exp $
////////////////////////////////////////////////////
/*!
  \file AliJTrack.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.1 $
  \date $Date: 2008/05/02 11:56:39 $
*/
////////////////////////////////////////////////////

#ifndef ALIJTRACK_H
#define ALIJTRACK_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "AliJBaseTrack.h"

class AliJTrack : public AliJBaseTrack {

public:
   enum AliJTrkPID {
    kMostProbableAliJ = -1, 
    kElectronAliJ, 
    kMuonAliJ,
    kPionAliJ, 
    kKaonAliJ, 
    kProtonAliJ, 
    kNAliJTrkPID
   };

   enum AliJTrkPIDmethod {
    //kTOF=0, kTPC=1, kTPCTOF=2, kNAliJTrkPIDmethod
    kTOF=0, kNAliJTrkPIDmethod
   };

  AliJTrack();             // default constructor
  AliJTrack(const AliJTrack& a); // copy constructor
  ~AliJTrack(){;}    //destructor

  Double32_t  GetPID(AliJTrkPID p, AliJTrkPIDmethod m) const { return fTrkPID[p][m]; }
  UInt_t      GetFilterMap() const { return fFilterMap; }
  Bool_t      IsFiltered( int i ) const { return TESTBIT( fFilterMap, i); }
  Bool_t      IsFilteredMask( UInt_t mask ) const { return ((Bool_t)(((fFilterMap) & mask) != 0)); };
  int         GetTPCnClust() const {return fTPCnClust;}
  Double32_t     GetTPCdEdx()  const {return fTPCdEdx; }
  Double32_t GetTOFsignal() const {return fTOFsignal;}
  Double32_t GetExpectedTOFsignal( AliJTrkPID p ) const {return fExpTOFsignal[p];}
  Double32_t GetTPCmomentum() const {return fTPCmom;}
  void GetTrackPos( Double32_t *p ) const { p[0] = fTrackPos[0]; p[1] = fTrackPos[1]; p[2] = fTrackPos[2]; }

  void SetPID(AliJTrkPID p, double pro, AliJTrkPIDmethod m){ fTrkPID[p][m] = pro; }
  void SetFilterMap( UInt_t map ){ fFilterMap = map; }
  void SetFilterMap( int i, bool t ){ t?SETBIT(fFilterMap,i):CLRBIT(fFilterMap,i); } 
  void SetTPCnClust(int ival) {fTPCnClust = ival;}
  void SetTPCdEdx(Double32_t dedx ){ fTPCdEdx = dedx; }
  void SetTOFsignal(Double_t tofsig) {fTOFsignal = tofsig;}
  void SetExpectedTOFsignal(AliJTrkPID p, double extime) {fExpTOFsignal[p] = extime;}
  void SetTPCmomentum(Double_t momTPC) {fTPCmom = momTPC;}
  void SetTPCTrack( Double32_t px,  Double32_t py,  Double32_t pz ){ fTPCTrack[0]=px;fTPCTrack[1]=py;fTPCTrack[2]=pz;}
  void SetGCGTrack( Double32_t px,  Double32_t py,  Double32_t pz ){ fGCGTrack[0]=px;fGCGTrack[1]=py;fGCGTrack[2]=pz;}
  void SetTrackPos( Double32_t *p ){ fTrackPos[0]=p[0];fTrackPos[1]=p[1];fTrackPos[2]=p[2];}
  void SetUseTPCTrack(){ SetPxPyPzE( fTPCTrack[0], fTPCTrack[1], fTPCTrack[2], 0 ); } 
  void SetUseGCGTrack(){ SetPxPyPzE( fGCGTrack[0], fGCGTrack[1], fGCGTrack[2], 0 ); } 
  Double32_t * GetTPCTrack(){ return fTPCTrack; }
  Double32_t * GetGCGTrack(){ return fGCGTrack; }

  Double32_t  PtTPC(){ return TMath::Sqrt( fTPCTrack[0]*fTPCTrack[0] +  fTPCTrack[1]*fTPCTrack[1] ); }
  Double32_t  PtGCG(){ return TMath::Sqrt( fGCGTrack[0]*fGCGTrack[0] +  fGCGTrack[1]*fGCGTrack[1] ); }


  AliJTrack& operator=(const AliJTrack& trk);

  //Double32_t GetDCAtoVertexXY(){ return fDCAtoVertexXY; }
  //void SetDCAtoVertexXY( double dxy ){ fDCAtoVertexXY = dxy; }
  //Double32_t GetDCAtoVertexZ(){ return fDCAtoVertexZ; }
  //void SetDCAtoVertexZ( double dxy ){ fDCAtoVertexZ = dxy; }

private:
  Double32_t    fTrkPID[kNAliJTrkPID][kNAliJTrkPIDmethod];   //[0.,1.,8] Array for PID. 
  UInt_t        fFilterMap;                 // bit serious of cuts
  Short_t       fTPCnClust;                 // track TPC nclusters 
  Double32_t    fTPCdEdx;                   // TPC dEdx
  Double32_t	fTOFsignal;		// TOF time
  Double32_t	fExpTOFsignal[kNAliJTrkPID];	// expected TOF time
  Double32_t	fTPCmom;		// TPC momentum to calculate expected TPC dEdx
  Double32_t    fTPCTrack[3];               // px, py, pz for TPCTrack;
  Double32_t    fGCGTrack[3];               // px, py, pz for GCGTrack;
  Double32_t    fTrackPos[3];               // track position

  //Double32_t    fDCAtoVertexXY;             //!
  //Double32_t    fDCAtoVertexZ;             //!

  ClassDef(AliJTrack,1)
};

#endif
