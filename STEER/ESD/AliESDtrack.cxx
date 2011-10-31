/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-----------------------------------------------------------------
//           Implementation of the ESD track class
//   ESD = Event Summary Data
//   This is the class to deal with during the phisics analysis of data
//      Origin: Iouri Belikov, CERN
//      e-mail: Jouri.Belikov@cern.ch
//
//
//
//  What do you need to know before starting analysis
//  (by Marian Ivanov: marian.ivanov@cern.ch)
//
//
//   AliESDtrack:
//   1.  What is the AliESDtrack
//   2.  What informations do we store
//   3.  How to use the information for analysis
//   
//
//   1.AliESDtrack is the container of the information about the track/particle
//     reconstructed during Barrel Tracking.
//     The track information is propagated from one tracking detector to 
//     other using the functionality of AliESDtrack - Current parameters.  
//
//     No global fit model is used.
//     Barrel tracking use Kalman filtering technique, it gives optimal local 
//     track parameters at given point under certian assumptions.
//             
//     Kalman filter take into account additional effect which are 
//     difficult to handle using global fit.
//     Effects:
//        a.) Multiple scattering
//        b.) Energy loss
//        c.) Non homogenous magnetic field
//
//     In general case, following barrel detectors are contributing to 
//     the Kalman track information:
//         a. TPC
//         b. ITS
//         c. TRD
//
//      In general 3 reconstruction itteration are performed:
//         1. Find tracks   - sequence TPC->ITS
//         2. PropagateBack - sequence ITS->TPC->TRD -> Outer PID detectors
//         3. Refit invward - sequence TRD->TPC->ITS
//      The current tracks are updated after each detector (see bellow).
//      In specical cases a track  sanpshots are stored.   
// 
//
//      For some type of analysis (+visualization) track local parameters at 
//      different position are neccesary. A snapshots during the track 
//      propagation are created.
//      (See AliExternalTrackParam class for desctiption of variables and 
//      functionality)
//      Snapshots:
//      a. Current parameters - class itself (AliExternalTrackParam)
//         Contributors: general case TRD->TPC->ITS
//         Preferable usage:  Decission  - primary or secondary track
//         NOTICE - By default the track parameters are stored at the DCA point
//                  to the primary vertex. optimal for primary tracks, 
//                  far from optimal for secondary tracks.
//      b. Constrained parameters - Kalman information updated with 
//         the Primary vertex information 
//         Contributors: general case TRD->TPC->ITS
//         Preferable usage: Use only for tracks selected as primary
//         NOTICE - not real constrain - taken as additional measurement 
//         with corresponding error
//         Function:  
//       const AliExternalTrackParam *GetConstrainedParam() const {return fCp;}
//      c. Inner parameters -  Track parameters at inner wall of the TPC 
//         Contributors: general case TRD->TPC
//         function:
//           const AliExternalTrackParam *GetInnerParam() const { return fIp;}
//
//      d. TPCinnerparam  - contributors - TPC only
//         Contributors:  TPC
//         Preferable usage: Requested for HBT study 
//                         (smaller correlations as using also ITS information)
//         NOTICE - the track parameters are propagated to the DCA to  
//         to primary vertex
//         Optimal for primary, far from optimal for secondary tracks
//         Function:
//    const AliExternalTrackParam *GetTPCInnerParam() const {return fTPCInner;}
//     
//      e. Outer parameters - 
//           Contributors-  general case - ITS-> TPC -> TRD
//           The last point - Outer parameters radius is determined
//           e.a) Local inclination angle bigger than threshold - 
//                Low momenta tracks 
//           e.a) Catastrofic energy losss in material
//           e.b) Not further improvement (no space points)
//           Usage:             
//              a.) Tracking: Starting parameter for Refit inward 
//              b.) Visualization
//              c.) QA
//         NOTICE: Should be not used for the physic analysis
//         Function:
//            const AliExternalTrackParam *GetOuterParam() const { return fOp;}
//
//-----------------------------------------------------------------

#include <TMath.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TMatrixD.h>

#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliKalmanTrack.h"
#include "AliVTrack.h"
#include "AliLog.h"
#include "AliTrackPointArray.h"
#include "TPolyMarker3D.h"
#include "AliTrackerBase.h"
#include "AliTPCdEdxInfo.h"

ClassImp(AliESDtrack)

void SetPIDValues(Double_t * dest, const Double_t * src, Int_t n) {
  // This function copies "n" PID weights from "scr" to "dest"
  // and normalizes their sum to 1 thus producing conditional probabilities.
  // The negative weights are set to 0.
  // In case all the weights are non-positive they are replaced by
  // uniform probabilities

  if (n<=0) return;

  Float_t uniform = 1./(Float_t)n;

  Float_t sum = 0;
  for (Int_t i=0; i<n; i++) 
    if (src[i]>=0) {
      sum+=src[i];
      dest[i] = src[i];
    }
    else {
      dest[i] = 0;
    }

  if(sum>0)
    for (Int_t i=0; i<n; i++) dest[i] /= sum;
  else
    for (Int_t i=0; i<n; i++) dest[i] = uniform;
}

//_______________________________________________________________________
AliESDtrack::AliESDtrack() : 
  AliExternalTrackParam(),
  fCp(0),
  fIp(0),
  fTPCInner(0),
  fOp(0),
  fHMPIDp(0),  
  fFriendTrack(NULL),
  fTPCFitMap(159),//number of padrows
  fTPCClusterMap(159),//number of padrows
  fTPCSharedMap(159),//number of padrows
  fFlags(0),
  fID(0),
  fLabel(0),
  fITSLabel(0),
  fTPCLabel(0),
  fTRDLabel(0),
  fTOFCalChannel(-1),
  fTOFindex(-1),
  fHMPIDqn(0),
  fHMPIDcluIdx(-1),
  fCaloIndex(kEMCALNoMatch),
  fHMPIDtrkTheta(0),
  fHMPIDtrkPhi(0),
  fHMPIDsignal(0),
  fTrackLength(0),
  fdTPC(0),fzTPC(0),
  fCddTPC(0),fCdzTPC(0),fCzzTPC(0),
  fCchi2TPC(0),
  fD(0),fZ(0),
  fCdd(0),fCdz(0),fCzz(0),
  fCchi2(0),
  fITSchi2(0),
  fTPCchi2(0),
  fTPCchi2Iter1(0),
  fTRDchi2(0),
  fTOFchi2(0),
  fHMPIDchi2(0),
  fGlobalChi2(0),
  fITSsignal(0),
  fTPCsignal(0),
  fTPCsignalS(0),
  fTPCdEdxInfo(0),
  fTRDsignal(0),
  fTRDQuality(0),
  fTRDBudget(0),
  fTOFsignal(99999),
  fTOFsignalToT(99999),
  fTOFsignalRaw(99999),
  fTOFsignalDz(999),
  fTOFsignalDx(999),
  fTOFdeltaBC(999),
  fTOFl0l1(999),
  fCaloDx(0),
  fCaloDz(0),
  fHMPIDtrkX(0),
  fHMPIDtrkY(0),
  fHMPIDmipX(0),
  fHMPIDmipY(0),
  fTPCncls(0),
  fTPCnclsF(0),
  fTPCsignalN(0),
  fTPCnclsIter1(0),
  fTPCnclsFIter1(0),
  fITSncls(0),
  fITSClusterMap(0),
  fITSSharedMap(0),
  fTRDncls(0),
  fTRDncls0(0),
  fTRDntracklets(0),
  fTRDnSlices(0),
  fTRDslices(0x0),
  fVertexID(-2),// -2 means an orphan track 
  fESDEvent(0),
  fCacheNCrossedRows(-10),
  fCacheChi2TPCConstrainedVsGlobal(-10),
  fCacheChi2TPCConstrainedVsGlobalVertex(0)
{
  //
  // The default ESD constructor 
  //
  if (!OnlineMode()) fFriendTrack=new AliESDfriendTrack();

  Int_t i;
  for (i=kNITSchi2Std;i--;) fITSchi2Std[i] = 0;
  for (i=0; i<AliPID::kSPECIES; i++) {
    fTrackTime[i]=0.;
    fR[i]=0.;
    fITSr[i]=0.;
    fTPCr[i]=0.;
    fTRDr[i]=0.;
    fTOFr[i]=0.;
    fHMPIDr[i]=0.;
  }
  
  for (i=0; i<3; i++)   { fKinkIndexes[i]=0;}
  for (i=0; i<3; i++)   { fV0Indexes[i]=0;}
  for (i=0;i<kTRDnPlanes;i++) {
    fTRDTimBin[i]=0;
  }
  for (i=0;i<4;i++) {fITSdEdxSamples[i]=0.;}
  for (i=0;i<4;i++) {fTPCPoints[i]=0;}
  for (i=0;i<3;i++) {fTOFLabel[i]=-1;}
  for (i=0;i<10;i++) {fTOFInfo[i]=0;}
  for (i=0;i<12;i++) {fITSModule[i]=-1;}
}

bool AliESDtrack::fgkOnlineMode=false;

//_______________________________________________________________________
AliESDtrack::AliESDtrack(const AliESDtrack& track):
  AliExternalTrackParam(track),
  fCp(0),
  fIp(0),
  fTPCInner(0),
  fOp(0),
  fHMPIDp(0),  
  fFriendTrack(0),
  fTPCFitMap(track.fTPCFitMap),
  fTPCClusterMap(track.fTPCClusterMap),
  fTPCSharedMap(track.fTPCSharedMap),
  fFlags(track.fFlags),
  fID(track.fID),
  fLabel(track.fLabel),
  fITSLabel(track.fITSLabel),
  fTPCLabel(track.fTPCLabel),
  fTRDLabel(track.fTRDLabel),
  fTOFCalChannel(track.fTOFCalChannel),
  fTOFindex(track.fTOFindex),
  fHMPIDqn(track.fHMPIDqn),
  fHMPIDcluIdx(track.fHMPIDcluIdx),
  fCaloIndex(track.fCaloIndex),
  fHMPIDtrkTheta(track.fHMPIDtrkTheta),
  fHMPIDtrkPhi(track.fHMPIDtrkPhi),
  fHMPIDsignal(track.fHMPIDsignal),
  fTrackLength(track.fTrackLength),
  fdTPC(track.fdTPC),fzTPC(track.fzTPC),
  fCddTPC(track.fCddTPC),fCdzTPC(track.fCdzTPC),fCzzTPC(track.fCzzTPC),
  fCchi2TPC(track.fCchi2TPC),
  fD(track.fD),fZ(track.fZ),
  fCdd(track.fCdd),fCdz(track.fCdz),fCzz(track.fCzz),
  fCchi2(track.fCchi2),
  fITSchi2(track.fITSchi2),
  fTPCchi2(track.fTPCchi2),
  fTPCchi2Iter1(track.fTPCchi2Iter1),
  fTRDchi2(track.fTRDchi2),
  fTOFchi2(track.fTOFchi2),
  fHMPIDchi2(track.fHMPIDchi2),
  fGlobalChi2(track.fGlobalChi2),
  fITSsignal(track.fITSsignal),
  fTPCsignal(track.fTPCsignal),
  fTPCsignalS(track.fTPCsignalS),
  fTPCdEdxInfo(0),
  fTRDsignal(track.fTRDsignal),
  fTRDQuality(track.fTRDQuality),
  fTRDBudget(track.fTRDBudget),
  fTOFsignal(track.fTOFsignal),
  fTOFsignalToT(track.fTOFsignalToT),
  fTOFsignalRaw(track.fTOFsignalRaw),
  fTOFsignalDz(track.fTOFsignalDz),
  fTOFsignalDx(track.fTOFsignalDx),
  fTOFdeltaBC(track.fTOFdeltaBC),
  fTOFl0l1(track.fTOFl0l1),
  fCaloDx(track.fCaloDx),
  fCaloDz(track.fCaloDz),
  fHMPIDtrkX(track.fHMPIDtrkX),
  fHMPIDtrkY(track.fHMPIDtrkY),
  fHMPIDmipX(track.fHMPIDmipX),
  fHMPIDmipY(track.fHMPIDmipY),
  fTPCncls(track.fTPCncls),
  fTPCnclsF(track.fTPCnclsF),
  fTPCsignalN(track.fTPCsignalN),
  fTPCnclsIter1(track.fTPCnclsIter1),
  fTPCnclsFIter1(track.fTPCnclsIter1),
  fITSncls(track.fITSncls),
  fITSClusterMap(track.fITSClusterMap),
  fITSSharedMap(track.fITSSharedMap),
  fTRDncls(track.fTRDncls),
  fTRDncls0(track.fTRDncls0),
  fTRDntracklets(track.fTRDntracklets),
  fTRDnSlices(track.fTRDnSlices),
  fTRDslices(0x0),
  fVertexID(track.fVertexID),
  fESDEvent(track.fESDEvent),
  fCacheNCrossedRows(track.fCacheNCrossedRows),
  fCacheChi2TPCConstrainedVsGlobal(track.fCacheChi2TPCConstrainedVsGlobal),
  fCacheChi2TPCConstrainedVsGlobalVertex(track.fCacheChi2TPCConstrainedVsGlobalVertex)
{
  //
  //copy constructor
  //
  for (Int_t i=kNITSchi2Std;i--;) fITSchi2Std[i] = track.fTrackTime[i];
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i]=track.fTrackTime[i];
  for (Int_t i=0;i<AliPID::kSPECIES;i++)  fR[i]=track.fR[i];
  //
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fITSr[i]=track.fITSr[i]; 
  //
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTPCr[i]=track.fTPCr[i]; 
  for (Int_t i=0;i<4;i++) {fITSdEdxSamples[i]=track.fITSdEdxSamples[i];}
  for (Int_t i=0;i<4;i++) {fTPCPoints[i]=track.fTPCPoints[i];}
  for (Int_t i=0; i<3;i++)   { fKinkIndexes[i]=track.fKinkIndexes[i];}
  for (Int_t i=0; i<3;i++)   { fV0Indexes[i]=track.fV0Indexes[i];}
  //
  for (Int_t i=0;i<kTRDnPlanes;i++) {
    fTRDTimBin[i]=track.fTRDTimBin[i];
  }

  if (fTRDnSlices) {
    fTRDslices=new Double32_t[fTRDnSlices];
    for (Int_t i=0; i<fTRDnSlices; i++) fTRDslices[i]=track.fTRDslices[i];
  }


  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTRDr[i]=track.fTRDr[i]; 
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTOFr[i]=track.fTOFr[i];
  for (Int_t i=0;i<3;i++) fTOFLabel[i]=track.fTOFLabel[i];
  for (Int_t i=0;i<10;i++) fTOFInfo[i]=track.fTOFInfo[i];
  for (Int_t i=0;i<12;i++) fITSModule[i]=track.fITSModule[i];
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fHMPIDr[i]=track.fHMPIDr[i];

  if (track.fCp) fCp=new AliExternalTrackParam(*track.fCp);
  if (track.fIp) fIp=new AliExternalTrackParam(*track.fIp);
  if (track.fTPCInner) fTPCInner=new AliExternalTrackParam(*track.fTPCInner);
  if (track.fOp) fOp=new AliExternalTrackParam(*track.fOp);
  if (track.fHMPIDp) fHMPIDp=new AliExternalTrackParam(*track.fHMPIDp);
  if (track.fTPCdEdxInfo) fTPCdEdxInfo = new AliTPCdEdxInfo(*track.fTPCdEdxInfo);

  
  if (track.fFriendTrack) fFriendTrack=new AliESDfriendTrack(*(track.fFriendTrack));
}

//_______________________________________________________________________
AliESDtrack::AliESDtrack(const AliVTrack *track) : 
  AliExternalTrackParam(track),
  fCp(0),
  fIp(0),
  fTPCInner(0),
  fOp(0),
  fHMPIDp(0),  
  fFriendTrack(0),
  fTPCFitMap(159),//number of padrows
  fTPCClusterMap(159),//number of padrows
  fTPCSharedMap(159),//number of padrows
  fFlags(0),
  fID(),
  fLabel(0),
  fITSLabel(0),
  fTPCLabel(0),
  fTRDLabel(0),
  fTOFCalChannel(-1),
  fTOFindex(-1),
  fHMPIDqn(0),
  fHMPIDcluIdx(-1),
  fCaloIndex(kEMCALNoMatch),
  fHMPIDtrkTheta(0),
  fHMPIDtrkPhi(0),
  fHMPIDsignal(0),
  fTrackLength(0),
  fdTPC(0),fzTPC(0),
  fCddTPC(0),fCdzTPC(0),fCzzTPC(0),
  fCchi2TPC(0),
  fD(0),fZ(0),
  fCdd(0),fCdz(0),fCzz(0),
  fCchi2(0),
  fITSchi2(0),
  fTPCchi2(0),
  fTPCchi2Iter1(0),
  fTRDchi2(0),
  fTOFchi2(0),
  fHMPIDchi2(0),
  fGlobalChi2(0),
  fITSsignal(0),
  fTPCsignal(0),
  fTPCsignalS(0),
  fTPCdEdxInfo(0),
  fTRDsignal(0),
  fTRDQuality(0),
  fTRDBudget(0),
  fTOFsignal(99999),
  fTOFsignalToT(99999),
  fTOFsignalRaw(99999),
  fTOFsignalDz(999),
  fTOFsignalDx(999),
  fTOFdeltaBC(999),
  fTOFl0l1(999),
  fCaloDx(0),
  fCaloDz(0),
  fHMPIDtrkX(0),
  fHMPIDtrkY(0),
  fHMPIDmipX(0),
  fHMPIDmipY(0),
  fTPCncls(0),
  fTPCnclsF(0),
  fTPCsignalN(0),
  fTPCnclsIter1(0),
  fTPCnclsFIter1(0),
  fITSncls(0),
  fITSClusterMap(0),
  fITSSharedMap(0),
  fTRDncls(0),
  fTRDncls0(0),
  fTRDntracklets(0),
  fTRDnSlices(0),
  fTRDslices(0x0),
  fVertexID(-2),  // -2 means an orphan track
  fESDEvent(0),
  fCacheNCrossedRows(-10),
  fCacheChi2TPCConstrainedVsGlobal(-10),
  fCacheChi2TPCConstrainedVsGlobalVertex(0)
{
  //
  // ESD track from AliVTrack.
  // This is not a copy constructor !
  //

  if (track->InheritsFrom("AliExternalTrackParam")) {
     AliError("This is not a copy constructor. Use AliESDtrack(const AliESDtrack &) !");
     AliWarning("Calling the default constructor...");
     AliESDtrack();
     return;
  }

  // Reset all the arrays
  Int_t i;
  for (i=kNITSchi2Std;i--;) fITSchi2Std[i] = 0;
  for (i=0; i<AliPID::kSPECIES; i++) {
    fTrackTime[i]=0.;
    fR[i]=0.;
    fITSr[i]=0.;
    fTPCr[i]=0.;
    fTRDr[i]=0.;
    fTOFr[i]=0.;
    fHMPIDr[i]=0.;
  }
  
  for (i=0; i<3; i++)   { fKinkIndexes[i]=0;}
  for (i=0; i<3; i++)   { fV0Indexes[i]=-1;}
  for (i=0;i<kTRDnPlanes;i++) {
    fTRDTimBin[i]=0;
  }
  for (i=0;i<4;i++) {fITSdEdxSamples[i]=0.;}
  for (i=0;i<4;i++) {fTPCPoints[i]=0;}
  for (i=0;i<3;i++) {fTOFLabel[i]=-1;}
  for (i=0;i<10;i++) {fTOFInfo[i]=0;}
  for (i=0;i<12;i++) {fITSModule[i]=-1;}

  // Set the ID
  SetID(track->GetID());

  // Set ITS cluster map
  fITSClusterMap=track->GetITSClusterMap();
  fITSSharedMap=0;

  fITSncls=0;
  for(i=0; i<6; i++) {
    if(HasPointOnITSLayer(i)) fITSncls++;
  }

  // Set TPC ncls 
  fTPCncls=track->GetTPCNcls();


  // Set the combined PID
  const Double_t *pid = track->PID();
  if(pid){
    for (i=0; i<AliPID::kSPECIES; i++) fR[i]=pid[i];
  }
  // AliESD track label
  SetLabel(track->GetLabel());
  // Set the status
  SetStatus(track->GetStatus());
}

//_______________________________________________________________________
AliESDtrack::AliESDtrack(TParticle * part) : 
  AliExternalTrackParam(),
  fCp(0),
  fIp(0),
  fTPCInner(0),
  fOp(0),
  fHMPIDp(0),  
  fFriendTrack(0),
  fTPCFitMap(159),//number of padrows
  fTPCClusterMap(159),//number of padrows
  fTPCSharedMap(159),//number of padrows
  fFlags(0),
  fID(0),
  fLabel(0),
  fITSLabel(0),
  fTPCLabel(0),
  fTRDLabel(0),
  fTOFCalChannel(-1),
  fTOFindex(-1),
  fHMPIDqn(0),
  fHMPIDcluIdx(-1),
  fCaloIndex(kEMCALNoMatch),
  fHMPIDtrkTheta(0),
  fHMPIDtrkPhi(0),
  fHMPIDsignal(0),
  fTrackLength(0),
  fdTPC(0),fzTPC(0),
  fCddTPC(0),fCdzTPC(0),fCzzTPC(0),
  fCchi2TPC(0),
  fD(0),fZ(0),
  fCdd(0),fCdz(0),fCzz(0),
  fCchi2(0),
  fITSchi2(0),
  fTPCchi2(0),
  fTPCchi2Iter1(0),  
  fTRDchi2(0),
  fTOFchi2(0),
  fHMPIDchi2(0),
  fGlobalChi2(0),
  fITSsignal(0),
  fTPCsignal(0),
  fTPCsignalS(0),
  fTPCdEdxInfo(0),
  fTRDsignal(0),
  fTRDQuality(0),
  fTRDBudget(0),
  fTOFsignal(99999),
  fTOFsignalToT(99999),
  fTOFsignalRaw(99999),
  fTOFsignalDz(999),
  fTOFsignalDx(999),
  fTOFdeltaBC(999),
  fTOFl0l1(999),
  fCaloDx(0),
  fCaloDz(0),
  fHMPIDtrkX(0),
  fHMPIDtrkY(0),
  fHMPIDmipX(0),
  fHMPIDmipY(0),
  fTPCncls(0),
  fTPCnclsF(0),
  fTPCsignalN(0),
  fTPCnclsIter1(0),
  fTPCnclsFIter1(0),
  fITSncls(0),
  fITSClusterMap(0),
  fITSSharedMap(0),
  fTRDncls(0),
  fTRDncls0(0),
  fTRDntracklets(0),
  fTRDnSlices(0),
  fTRDslices(0x0),
  fVertexID(-2),  // -2 means an orphan track
  fESDEvent(0),
  fCacheNCrossedRows(-10),
  fCacheChi2TPCConstrainedVsGlobal(-10),
  fCacheChi2TPCConstrainedVsGlobalVertex(0)  
{
  //
  // ESD track from TParticle
  //

  // Reset all the arrays
  Int_t i;
  for (i=kNITSchi2Std;i--;) fITSchi2Std[i] = 0;
  for (i=0; i<AliPID::kSPECIES; i++) {
    fTrackTime[i]=0.;
    fR[i]=0.;
    fITSr[i]=0.;
    fTPCr[i]=0.;
    fTRDr[i]=0.;
    fTOFr[i]=0.;
    fHMPIDr[i]=0.;
  }
  
  for (i=0; i<3; i++)   { fKinkIndexes[i]=0;}
  for (i=0; i<3; i++)   { fV0Indexes[i]=-1;}
  for (i=0;i<kTRDnPlanes;i++) {
    fTRDTimBin[i]=0;
  }
  for (i=0;i<4;i++) {fITSdEdxSamples[i]=0.;}
  for (i=0;i<4;i++) {fTPCPoints[i]=0;}
  for (i=0;i<3;i++) {fTOFLabel[i]=-1;}
  for (i=0;i<10;i++) {fTOFInfo[i]=0;}
  for (i=0;i<12;i++) {fITSModule[i]=-1;}

  // Calculate the AliExternalTrackParam content

  Double_t xref;
  Double_t alpha;
  Double_t param[5];
  Double_t covar[15];

  // Calculate alpha: the rotation angle of the corresponding local system (TPC sector)
  alpha = part->Phi()*180./TMath::Pi();
  if (alpha<0) alpha+= 360.;
  if (alpha>360) alpha -= 360.;

  Int_t sector = (Int_t)(alpha/20.);
  alpha = 10. + 20.*sector;
  alpha /= 180;
  alpha *= TMath::Pi();

  // Covariance matrix: no errors, the parameters are exact
  for (i=0; i<15; i++) covar[i]=0.;

  // Get the vertex of origin and the momentum
  TVector3 ver(part->Vx(),part->Vy(),part->Vz());
  TVector3 mom(part->Px(),part->Py(),part->Pz());

  // Rotate to the local coordinate system (TPC sector)
  ver.RotateZ(-alpha);
  mom.RotateZ(-alpha);

  // X of the referense plane
  xref = ver.X();

  Int_t pdgCode = part->GetPdgCode();

  Double_t charge = 
    TDatabasePDG::Instance()->GetParticle(pdgCode)->Charge();

  param[0] = ver.Y();
  param[1] = ver.Z();
  param[2] = TMath::Sin(mom.Phi());
  param[3] = mom.Pz()/mom.Pt();
  param[4] = TMath::Sign(1/mom.Pt(),charge);

  // Set AliExternalTrackParam
  Set(xref, alpha, param, covar);

  // Set the PID
  Int_t indexPID = 99;

  switch (TMath::Abs(pdgCode)) {

  case  11: // electron
    indexPID = 0;
    break;

  case 13: // muon
    indexPID = 1;
    break;

  case 211: // pion
    indexPID = 2;
    break;

  case 321: // kaon
    indexPID = 3;
    break;

  case 2212: // proton
    indexPID = 4;
    break;

  default:
    break;
  }

  // If the particle is not e,mu,pi,K or p the PID probabilities are set to 0
  if (indexPID < AliPID::kSPECIES) {
    fR[indexPID]=1.;
    fITSr[indexPID]=1.;
    fTPCr[indexPID]=1.;
    fTRDr[indexPID]=1.;
    fTOFr[indexPID]=1.;
    fHMPIDr[indexPID]=1.;

  }
  // AliESD track label
  SetLabel(part->GetUniqueID());

}

//_______________________________________________________________________
AliESDtrack::~AliESDtrack(){ 
  //
  // This is destructor according Coding Conventrions 
  //
  //printf("Delete track\n");
  delete fIp; 
  delete fTPCInner; 
  delete fOp;
  delete fHMPIDp;
  delete fCp; 
  if (fFriendTrack) delete fFriendTrack;
  if (fTPCdEdxInfo) delete fTPCdEdxInfo;
  fFriendTrack=NULL;
  if(fTRDnSlices)
    delete[] fTRDslices;
}

AliESDtrack &AliESDtrack::operator=(const AliESDtrack &source){
  

  if(&source == this) return *this;
  AliExternalTrackParam::operator=(source);

  
  if(source.fCp){
    // we have the trackparam: assign or copy construct
    if(fCp)*fCp = *source.fCp;
    else fCp = new AliExternalTrackParam(*source.fCp);
  }
  else{
    // no track param delete the old one
    if(fCp)delete fCp;
    fCp = 0;
  }

  if(source.fIp){
    // we have the trackparam: assign or copy construct
    if(fIp)*fIp = *source.fIp;
    else fIp = new AliExternalTrackParam(*source.fIp);
  }
  else{
    // no track param delete the old one
    if(fIp)delete fIp;
    fIp = 0;
  }


  if(source.fTPCInner){
    // we have the trackparam: assign or copy construct
    if(fTPCInner) *fTPCInner = *source.fTPCInner;
    else fTPCInner = new AliExternalTrackParam(*source.fTPCInner);
  }
  else{
    // no track param delete the old one
    if(fTPCInner)delete fTPCInner;
    fTPCInner = 0;
  }

  if(source.fTPCdEdxInfo) {
    if(fTPCdEdxInfo) *fTPCdEdxInfo = *source.fTPCdEdxInfo;
    fTPCdEdxInfo = new AliTPCdEdxInfo(*source.fTPCdEdxInfo);
  }

  if(source.fOp){
    // we have the trackparam: assign or copy construct
    if(fOp) *fOp = *source.fOp;
    else fOp = new AliExternalTrackParam(*source.fOp);
  }
  else{
    // no track param delete the old one
    if(fOp)delete fOp;
    fOp = 0;
  }

  
  if(source.fHMPIDp){
    // we have the trackparam: assign or copy construct
    if(fHMPIDp) *fHMPIDp = *source.fHMPIDp;
    else fHMPIDp = new AliExternalTrackParam(*source.fHMPIDp);
  }
  else{
    // no track param delete the old one
    if(fHMPIDp)delete fHMPIDp;
    fHMPIDp = 0;
  }

  // copy also the friend track 
  // use copy constructor
  if(source.fFriendTrack){
    // we have the trackparam: assign or copy construct
    delete fFriendTrack; fFriendTrack=new AliESDfriendTrack(*source.fFriendTrack);
  }
  else{
    // no track param delete the old one
    delete fFriendTrack; fFriendTrack= 0;
  }

  fTPCFitMap = source.fTPCFitMap; 
  fTPCClusterMap = source.fTPCClusterMap; 
  fTPCSharedMap  = source.fTPCSharedMap;  
  // the simple stuff
  fFlags    = source.fFlags; 
  fID       = source.fID;             
  fLabel    = source.fLabel;
  fITSLabel = source.fITSLabel;
  for(int i = 0; i< 12;++i){
    fITSModule[i] = source.fITSModule[i];
  }
  fTPCLabel = source.fTPCLabel; 
  fTRDLabel = source.fTRDLabel;
  for(int i = 0; i< 3;++i){
    fTOFLabel[i] = source.fTOFLabel[i];    
  }
  fTOFCalChannel = source.fTOFCalChannel;
  fTOFindex      = source.fTOFindex;
  fHMPIDqn       = source.fHMPIDqn;
  fHMPIDcluIdx   = source.fHMPIDcluIdx; 
  fCaloIndex    = source.fCaloIndex;
  for (int i=kNITSchi2Std;i--;) fITSchi2Std[i] = source.fITSchi2Std[i];
  for(int i = 0; i< 3;++i){
    fKinkIndexes[i] = source.fKinkIndexes[i]; 
    fV0Indexes[i]   = source.fV0Indexes[i]; 
  }

  for(int i = 0; i< AliPID::kSPECIES;++i){
    fR[i]     = source.fR[i];
    fITSr[i]  = source.fITSr[i];
    fTPCr[i]  = source.fTPCr[i];
    fTRDr[i]  = source.fTRDr[i];
    fTOFr[i]  = source.fTOFr[i];
    fHMPIDr[i] = source.fHMPIDr[i];
    fTrackTime[i] = source.fTrackTime[i];  
  }

  fHMPIDtrkTheta = source.fHMPIDtrkTheta;
  fHMPIDtrkPhi   = source.fHMPIDtrkPhi;
  fHMPIDsignal   = source.fHMPIDsignal; 

  
  fTrackLength   = source. fTrackLength;
  fdTPC  = source.fdTPC; 
  fzTPC  = source.fzTPC; 
  fCddTPC = source.fCddTPC;
  fCdzTPC = source.fCdzTPC;
  fCzzTPC = source.fCzzTPC;
  fCchi2TPC = source.fCchi2TPC;

  fD  = source.fD; 
  fZ  = source.fZ; 
  fCdd = source.fCdd;
  fCdz = source.fCdz;
  fCzz = source.fCzz;
  fCchi2     = source.fCchi2;

  fITSchi2   = source.fITSchi2;             
  fTPCchi2   = source.fTPCchi2;            
  fTPCchi2Iter1   = source.fTPCchi2Iter1;            
  fTRDchi2   = source.fTRDchi2;      
  fTOFchi2   = source.fTOFchi2;      
  fHMPIDchi2 = source.fHMPIDchi2;      

  fGlobalChi2 = source.fGlobalChi2;      

  fITSsignal  = source.fITSsignal;     
  for (Int_t i=0;i<4;i++) {fITSdEdxSamples[i]=source.fITSdEdxSamples[i];}
  fTPCsignal  = source.fTPCsignal;     
  fTPCsignalS = source.fTPCsignalS;    
  for(int i = 0; i< 4;++i){
    fTPCPoints[i] = source.fTPCPoints[i];  
  }
  fTRDsignal = source.fTRDsignal;

  for(int i = 0;i < kTRDnPlanes;++i){
    fTRDTimBin[i] = source.fTRDTimBin[i];   
  }

  if(fTRDnSlices)
    delete[] fTRDslices;
  fTRDslices=0;
  fTRDnSlices=source.fTRDnSlices;
  if (fTRDnSlices) {
    fTRDslices=new Double32_t[fTRDnSlices];
    for(int j = 0;j < fTRDnSlices;++j) fTRDslices[j] = source.fTRDslices[j];
  }

  fTRDQuality =   source.fTRDQuality;     
  fTRDBudget  =   source.fTRDBudget;      
  fTOFsignal  =   source.fTOFsignal;     
  fTOFsignalToT = source.fTOFsignalToT;   
  fTOFsignalRaw = source.fTOFsignalRaw;  
  fTOFsignalDz  = source.fTOFsignalDz;      
  fTOFsignalDx  = source.fTOFsignalDx;      
  fTOFdeltaBC   = source.fTOFdeltaBC;
  fTOFl0l1      = source.fTOFl0l1;
 
  for(int i = 0;i<10;++i){
    fTOFInfo[i] = source.fTOFInfo[i];    
  }

  fHMPIDtrkX = source.fHMPIDtrkX; 
  fHMPIDtrkY = source.fHMPIDtrkY; 
  fHMPIDmipX = source.fHMPIDmipX;
  fHMPIDmipY = source.fHMPIDmipY; 

  fTPCncls    = source.fTPCncls;      
  fTPCnclsF   = source.fTPCnclsF;     
  fTPCsignalN = source.fTPCsignalN;   
  fTPCnclsIter1    = source.fTPCnclsIter1;      
  fTPCnclsFIter1   = source.fTPCnclsFIter1;     

  fITSncls = source.fITSncls;       
  fITSClusterMap = source.fITSClusterMap; 
  fITSSharedMap = source.fITSSharedMap; 
  fTRDncls   = source.fTRDncls;       
  fTRDncls0  = source.fTRDncls0;      
  fTRDntracklets  = source.fTRDntracklets; 
  fVertexID = source.fVertexID;
  return *this;
}



void AliESDtrack::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDtrack *robj = dynamic_cast<AliESDtrack*>(&obj);
  if(!robj)return; // not an AliESDtrack
  *robj = *this;

}



void AliESDtrack::AddCalibObject(TObject * object){
  //
  // add calib object to the list
  //
  if (!fFriendTrack) fFriendTrack  = new AliESDfriendTrack;
  if (!fFriendTrack) return;
  fFriendTrack->AddCalibObject(object);
}

TObject *  AliESDtrack::GetCalibObject(Int_t index){
  //
  // return calib objct at given position
  //
  if (!fFriendTrack) return 0;
  return fFriendTrack->GetCalibObject(index);
}


Bool_t AliESDtrack::FillTPCOnlyTrack(AliESDtrack &track){
  
  // Fills the information of the TPC-only first reconstruction pass
  // into the passed ESDtrack object. For consistency fTPCInner is also filled
  // again



  // For data produced before r26675
  // RelateToVertexTPC was not properly called during reco
  // so you'll have to call it again, before FillTPCOnlyTrack
  //  Float_t p[2],cov[3];
  // track->GetImpactParametersTPC(p,cov); 
  // if(p[0]==0&&p[1]==0) // <- Default values
  //  track->RelateToVertexTPC(esd->GetPrimaryVertexTPC(),esd->GetMagneticField(),kVeryBig);
  

  if(!fTPCInner)return kFALSE;

  // fill the TPC track params to the global track parameters
  track.Set(fTPCInner->GetX(),fTPCInner->GetAlpha(),fTPCInner->GetParameter(),fTPCInner->GetCovariance());
  track.fD = fdTPC;
  track.fZ = fzTPC;
  track.fCdd = fCddTPC;
  track.fCdz = fCdzTPC;
  track.fCzz = fCzzTPC;

  // copy the inner params
  if(track.fIp) *track.fIp = *fIp;
  else track.fIp = new AliExternalTrackParam(*fIp);

  // copy the TPCinner parameters
  if(track.fTPCInner) *track.fTPCInner = *fTPCInner;
  else track.fTPCInner = new AliExternalTrackParam(*fTPCInner);
  track.fdTPC   = fdTPC;
  track.fzTPC   = fzTPC;
  track.fCddTPC = fCddTPC;
  track.fCdzTPC = fCdzTPC;
  track.fCzzTPC = fCzzTPC;
  track.fCchi2TPC = fCchi2TPC;

  // copy all other TPC specific parameters

  // replace label by TPC label
  track.fLabel    = fTPCLabel;
  track.fTPCLabel = fTPCLabel;

  track.fTPCchi2 = fTPCchi2; 
  track.fTPCchi2Iter1 = fTPCchi2Iter1; 
  track.fTPCsignal = fTPCsignal;
  track.fTPCsignalS = fTPCsignalS;
  for(int i = 0;i<4;++i)track.fTPCPoints[i] = fTPCPoints[i];

  track.fTPCncls    = fTPCncls;     
  track.fTPCnclsF   = fTPCnclsF;     
  track.fTPCsignalN =  fTPCsignalN;
  track.fTPCnclsIter1    = fTPCnclsIter1;     
  track.fTPCnclsFIter1   = fTPCnclsFIter1;     

  // PID 
  for(int i=0;i<AliPID::kSPECIES;++i){
    track.fTPCr[i] = fTPCr[i];
    // combined PID is TPC only!
    track.fR[i] = fTPCr[i];
  }
  track.fTPCFitMap = fTPCFitMap;
  track.fTPCClusterMap = fTPCClusterMap;
  track.fTPCSharedMap = fTPCSharedMap;


  // reset the flags
  track.fFlags = kTPCin;
  track.fID    = fID;

  track.fFlags |= fFlags & kTPCpid; //copy the TPCpid status flag
 
  for (Int_t i=0;i<3;i++) track.fKinkIndexes[i] = fKinkIndexes[i];
  
  return kTRUE;
    
}

//_______________________________________________________________________
void AliESDtrack::MakeMiniESDtrack(){
  // Resets everything except
  // fFlags: Reconstruction status flags 
  // fLabel: Track label
  // fID:  Unique ID of the track
  // Impact parameter information
  // fR[AliPID::kSPECIES]: combined "detector response probability"
  // Running track parameters in the base class (AliExternalTrackParam)
  
  fTrackLength = 0;

  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i] = 0;

  // Reset track parameters constrained to the primary vertex
  delete fCp;fCp = 0;

  // Reset track parameters at the inner wall of TPC
  delete fIp;fIp = 0;
  delete fTPCInner;fTPCInner=0;
  // Reset track parameters at the inner wall of the TRD
  delete fOp;fOp = 0;
  // Reset track parameters at the HMPID
  delete fHMPIDp;fHMPIDp = 0;


  // Reset ITS track related information
  fITSchi2 = 0;
  fITSncls = 0;       
  fITSClusterMap=0;
  fITSSharedMap=0;
  fITSsignal = 0;     
  for (Int_t i=0;i<4;i++) fITSdEdxSamples[i] = 0.;
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fITSr[i]=0; 
  fITSLabel = 0;       

  // Reset TPC related track information
  fTPCchi2 = 0;       
  fTPCchi2Iter1 = 0;       
  fTPCncls = 0;       
  fTPCnclsF = 0;       
  fTPCnclsIter1 = 0;       
  fTPCnclsFIter1 = 0;  
  fTPCFitMap = 0;       
  fTPCClusterMap = 0;  
  fTPCSharedMap = 0;  
  fTPCsignal= 0;      
  fTPCsignalS= 0;      
  fTPCsignalN= 0;      
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTPCr[i]=0; 
  fTPCLabel=0;       
  for (Int_t i=0;i<4;i++) fTPCPoints[i] = 0;
  for (Int_t i=0; i<3;i++)   fKinkIndexes[i] = 0;
  for (Int_t i=0; i<3;i++)   fV0Indexes[i] = 0;

  // Reset TRD related track information
  fTRDchi2 = 0;        
  fTRDncls = 0;       
  fTRDncls0 = 0;       
  fTRDsignal = 0;      
  for (Int_t i=0;i<kTRDnPlanes;i++) {
    fTRDTimBin[i]  = 0;
  }
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTRDr[i] = 0; 
  fTRDLabel = 0;       
  fTRDQuality  = 0;
  fTRDntracklets = 0;
  if(fTRDnSlices)
    delete[] fTRDslices;
  fTRDslices=0x0;
  fTRDnSlices=0;
  fTRDBudget  = 0;

  // Reset TOF related track information
  fTOFchi2 = 0;        
  fTOFindex = -1;       
  fTOFsignal = 99999;      
  fTOFCalChannel = -1;
  fTOFsignalToT = 99999;
  fTOFsignalRaw = 99999;
  fTOFsignalDz = 999;
  fTOFsignalDx = 999;
  fTOFdeltaBC = 999;
  fTOFl0l1 = 999;
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTOFr[i] = 0;
  for (Int_t i=0;i<3;i++) fTOFLabel[i] = -1;
  for (Int_t i=0;i<10;i++) fTOFInfo[i] = 0;

  // Reset HMPID related track information
  fHMPIDchi2 = 0;     
  fHMPIDqn = 0;     
  fHMPIDcluIdx = -1;     
  fHMPIDsignal = 0;     
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fHMPIDr[i] = 0;
  fHMPIDtrkTheta = 0;     
  fHMPIDtrkPhi = 0;      
  fHMPIDtrkX = 0;     
  fHMPIDtrkY = 0;      
  fHMPIDmipX = 0;
  fHMPIDmipY = 0;
  fCaloIndex = kEMCALNoMatch;

  // reset global track chi2
  fGlobalChi2 = 0;

  fVertexID = -2; // an orphan track

  delete fFriendTrack; fFriendTrack = 0;
} 

//_______________________________________________________________________
Int_t AliESDtrack::GetPID() const 
{
  // Returns the particle most probable id
  Int_t i;
  for (i=0; i<AliPID::kSPECIES-1; i++) if (fR[i] != fR[i+1]) break;
  //
  if (i == AliPID::kSPECIES-1) return AliPID::kPion;  // If all the probabilities are equal, return the pion mass
  //
  Float_t max=0.;
  Int_t k=-1;
  for (i=0; i<AliPID::kSPECIES; i++) if (fR[i]>max) {k=i; max=fR[i];}
  //
  if (k==0) { // dE/dx "crossing points" in the TPC
    Double_t p=GetP();
    if ((p>0.38)&&(p<0.48))
      if (fR[0]<fR[3]*10.) return AliPID::kKaon;
    if ((p>0.75)&&(p<0.85))
      if (fR[0]<fR[4]*10.) return AliPID::kProton;
    return AliPID::kElectron;
  }
  if (k==1) return AliPID::kMuon; 
  if (k==2||k==-1) return AliPID::kPion;
  if (k==3) return AliPID::kKaon;
  if (k==4) return AliPID::kProton;
  AliWarning("Undefined PID !");
  return AliPID::kPion;
}

//_______________________________________________________________________
Int_t AliESDtrack::GetTOFBunchCrossing(Double_t b) const 
{
  // Returns the number of bunch crossings after trigger (assuming 25ns spacing)
  const double kSpacing = 25e3; // min interbanch spacing
  const double kShift = 0;
  Int_t bcid = kTOFBCNA; // defualt one
  if (!IsOn(kTOFout) || !IsOn(kESDpid)) return bcid; // no info
  //
  double tdif = fTOFsignal;
  if (IsOn(kTIME)) { // integrated time info is there
    int pid = GetPID();
    tdif -= fTrackTime[pid];
  }
  else { // assume integrated time info from TOF radius and momentum
    const double kRTOF = 385.;
    const double kCSpeed = 3.e-2; // cm/ps
    double p = GetP();
    if (p<0.01) return bcid;
    double m = GetMass();
    double curv = GetC(b);
    double path = TMath::Abs(curv)>kAlmost0 ? // account for curvature
      2./curv*TMath::ASin(kRTOF*curv/2.)*TMath::Sqrt(1.+GetTgl()*GetTgl()) : kRTOF;
    tdif -= path/kCSpeed*TMath::Sqrt(1.+m*m/(p*p));
  }
  bcid = TMath::Nint((tdif - kShift)/kSpacing);
  return bcid;
}

//______________________________________________________________________________
Double_t AliESDtrack::M() const
{
  // Returns the assumed mass
  // (the pion mass, if the particle can't be identified properly).
  static Bool_t printerr=kTRUE;
  if (printerr) {
     AliWarning("WARNING !!! ... THIS WILL BE PRINTED JUST ONCE !!!");
     printerr = kFALSE;
     AliWarning("This is the ESD mass. Use it with care !"); 
  }
  return GetMass(); 
}
  
//______________________________________________________________________________
Double_t AliESDtrack::E() const
{
  // Returns the energy of the particle given its assumed mass.
  // Assumes the pion mass if the particle can't be identified properly.
  
  Double_t m = M();
  Double_t p = P();
  return TMath::Sqrt(p*p + m*m);
}

//______________________________________________________________________________
Double_t AliESDtrack::Y() const
{
  // Returns the rapidity of a particle given its assumed mass.
  // Assumes the pion mass if the particle can't be identified properly.
  
  Double_t e = E();
  Double_t pz = Pz();
  if (e != TMath::Abs(pz)) { // energy was not equal to pz
    return 0.5*TMath::Log((e+pz)/(e-pz));
  } else { // energy was equal to pz
    return -999.;
  }
}

//_______________________________________________________________________
Bool_t AliESDtrack::UpdateTrackParams(const AliKalmanTrack *t, ULong_t flags){
  //
  // This function updates track's running parameters 
  //
  Bool_t rc=kTRUE;

  SetStatus(flags);
  fLabel=t->GetLabel();

  if (t->IsStartedTimeIntegral()) {
    SetStatus(kTIME);
    Double_t times[10];t->GetIntegratedTimes(times); SetIntegratedTimes(times);
    SetIntegratedLength(t->GetIntegratedLength());
  }

  Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
  if (fFriendTrack) {
  if (flags==kITSout) fFriendTrack->SetITSOut(*t);
  if (flags==kTPCout) fFriendTrack->SetTPCOut(*t);
  if (flags==kTRDrefit) fFriendTrack->SetTRDIn(*t);
  }
  
  switch (flags) {
    
  case kITSin: 
    fITSchi2Std[0] = t->GetChi2();
    //
  case kITSout: 
    fITSchi2Std[1] = t->GetChi2();
  case kITSrefit:
    {
    fITSchi2Std[2] = t->GetChi2();
    fITSClusterMap=0;
    fITSncls=t->GetNumberOfClusters();
    if (fFriendTrack) {
    Int_t* indexITS = new Int_t[AliESDfriendTrack::kMaxITScluster];
    for (Int_t i=0;i<AliESDfriendTrack::kMaxITScluster;i++) {
	indexITS[i]=t->GetClusterIndex(i);

	if (i<fITSncls) {
	  Int_t l=(indexITS[i] & 0xf0000000) >> 28;
           SETBIT(fITSClusterMap,l);                 
        }
    }
    fFriendTrack->SetITSIndices(indexITS,AliESDfriendTrack::kMaxITScluster);
    delete [] indexITS;
    }

    fITSchi2=t->GetChi2();
    fITSsignal=t->GetPIDsignal();
    fITSLabel = t->GetLabel();
    // keep in fOp the parameters outside ITS for ITS stand-alone tracks 
    if (flags==kITSout) { 
      if (!fOp) fOp=new AliExternalTrackParam(*t);
      else 
        fOp->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
    }   
    }
    break;
    
  case kTPCin: case kTPCrefit:
    {
    fTPCLabel = t->GetLabel();
    if (flags==kTPCin)  {
      fTPCInner=new AliExternalTrackParam(*t); 
      fTPCnclsIter1=t->GetNumberOfClusters();    
      fTPCchi2Iter1=t->GetChi2();
    }
    if (!fIp) fIp=new AliExternalTrackParam(*t);
    else 
      fIp->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
    }
    // Intentionally no break statement; need to set general TPC variables as well
  case kTPCout:
    {
    if (flags & kTPCout){
      if (!fOp) fOp=new AliExternalTrackParam(*t);
      else 
        fOp->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
    }
    fTPCncls=t->GetNumberOfClusters();    
    fTPCchi2=t->GetChi2();
    
    if (fFriendTrack) {  // Copy cluster indices
      Int_t* indexTPC = new Int_t[AliESDfriendTrack::kMaxTPCcluster];
      for (Int_t i=0;i<AliESDfriendTrack::kMaxTPCcluster;i++)         
	indexTPC[i]=t->GetClusterIndex(i);
      fFriendTrack->SetTPCIndices(indexTPC,AliESDfriendTrack::kMaxTPCcluster);
      delete [] indexTPC;
    }
    fTPCsignal=t->GetPIDsignal();
    }
    break;

  case kTRDin: case kTRDrefit:
    break;
  case kTRDout:
    {
    fTRDLabel = t->GetLabel(); 
    fTRDchi2  = t->GetChi2();
    fTRDncls  = t->GetNumberOfClusters();
    if (fFriendTrack) {
      Int_t* indexTRD = new Int_t[AliESDfriendTrack::kMaxTRDcluster];
      for (Int_t i=0;i<AliESDfriendTrack::kMaxTRDcluster;i++) indexTRD[i]=-2;
      for (Int_t i=0;i<6;i++) indexTRD[i]=t->GetTrackletIndex(i);
      fFriendTrack->SetTRDIndices(indexTRD,AliESDfriendTrack::kMaxTRDcluster);
      delete [] indexTRD;
    }    
    
    fTRDsignal=t->GetPIDsignal();
    }
    break;
  case kTRDbackup:
    if (!fOp) fOp=new AliExternalTrackParam(*t);
    else 
      fOp->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
    fTRDncls0 = t->GetNumberOfClusters(); 
    break;
  case kTOFin: 
    break;
  case kTOFout: 
    break;
  case kTRDStop:
    break;
  case kHMPIDout:
  if (!fHMPIDp) fHMPIDp=new AliExternalTrackParam(*t);
    else 
      fHMPIDp->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
    break;
  default: 
    AliError("Wrong flag !");
    return kFALSE;
  }

  return rc;
}

//_______________________________________________________________________
void AliESDtrack::GetExternalParameters(Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters
  //---------------------------------------------------------------------
  x=GetX();
  for (Int_t i=0; i<5; i++) p[i]=GetParameter()[i];
}

//_______________________________________________________________________
void AliESDtrack::GetExternalCovariance(Double_t cov[15]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the cov. matrix
  //---------------------------------------------------------------------
  for (Int_t i=0; i<15; i++) cov[i]=AliExternalTrackParam::GetCovariance()[i];
}

//_______________________________________________________________________
Bool_t AliESDtrack::GetConstrainedExternalParameters
                 (Double_t &alpha, Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns the constrained external track parameters
  //---------------------------------------------------------------------
  if (!fCp) return kFALSE;
  alpha=fCp->GetAlpha();
  x=fCp->GetX();
  for (Int_t i=0; i<5; i++) p[i]=fCp->GetParameter()[i];
  return kTRUE;
}

//_______________________________________________________________________
Bool_t 
AliESDtrack::GetConstrainedExternalCovariance(Double_t c[15]) const {
  //---------------------------------------------------------------------
  // This function returns the constrained external cov. matrix
  //---------------------------------------------------------------------
  if (!fCp) return kFALSE;
  for (Int_t i=0; i<15; i++) c[i]=fCp->GetCovariance()[i];
  return kTRUE;
}

Bool_t
AliESDtrack::GetInnerExternalParameters
                 (Double_t &alpha, Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters 
  // at the inner layer of TPC
  //---------------------------------------------------------------------
  if (!fIp) return kFALSE;
  alpha=fIp->GetAlpha();
  x=fIp->GetX();
  for (Int_t i=0; i<5; i++) p[i]=fIp->GetParameter()[i];
  return kTRUE;
}

Bool_t 
AliESDtrack::GetInnerExternalCovariance(Double_t cov[15]) const {
 //---------------------------------------------------------------------
 // This function returns external representation of the cov. matrix 
 // at the inner layer of TPC
 //---------------------------------------------------------------------
  if (!fIp) return kFALSE;
  for (Int_t i=0; i<15; i++) cov[i]=fIp->GetCovariance()[i];
  return kTRUE;
}

void 
AliESDtrack::SetOuterParam(const AliExternalTrackParam *p, ULong_t flags) {
  //
  // This is a direct setter for the outer track parameters
  //
  SetStatus(flags);
  if (fOp) delete fOp;
  fOp=new AliExternalTrackParam(*p);
}

void 
AliESDtrack::SetOuterHmpParam(const AliExternalTrackParam *p, ULong_t flags) {
  //
  // This is a direct setter for the outer track parameters
  //
  SetStatus(flags);
  if (fHMPIDp) delete fHMPIDp;
  fHMPIDp=new AliExternalTrackParam(*p);
}

Bool_t 
AliESDtrack::GetOuterExternalParameters
                 (Double_t &alpha, Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters 
  // at the inner layer of TRD
  //---------------------------------------------------------------------
  if (!fOp) return kFALSE;
  alpha=fOp->GetAlpha();
  x=fOp->GetX();
  for (Int_t i=0; i<5; i++) p[i]=fOp->GetParameter()[i];
  return kTRUE;
}

Bool_t 
AliESDtrack::GetOuterHmpExternalParameters
                 (Double_t &alpha, Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters 
  // at the inner layer of TRD
  //---------------------------------------------------------------------
  if (!fHMPIDp) return kFALSE;
  alpha=fHMPIDp->GetAlpha();
  x=fHMPIDp->GetX();
  for (Int_t i=0; i<5; i++) p[i]=fHMPIDp->GetParameter()[i];
  return kTRUE;
}

Bool_t 
AliESDtrack::GetOuterExternalCovariance(Double_t cov[15]) const {
 //---------------------------------------------------------------------
 // This function returns external representation of the cov. matrix 
 // at the inner layer of TRD
 //---------------------------------------------------------------------
  if (!fOp) return kFALSE;
  for (Int_t i=0; i<15; i++) cov[i]=fOp->GetCovariance()[i];
  return kTRUE;
}

Bool_t 
AliESDtrack::GetOuterHmpExternalCovariance(Double_t cov[15]) const {
 //---------------------------------------------------------------------
 // This function returns external representation of the cov. matrix 
 // at the inner layer of TRD
 //---------------------------------------------------------------------
  if (!fHMPIDp) return kFALSE;
  for (Int_t i=0; i<15; i++) cov[i]=fHMPIDp->GetCovariance()[i];
  return kTRUE;
}

Int_t AliESDtrack::GetNcls(Int_t idet) const
{
  // Get number of clusters by subdetector index
  //
  Int_t ncls = 0;
  switch(idet){
  case 0:
    ncls = fITSncls;
    break;
  case 1:
    ncls = fTPCncls;
    break;
  case 2:
    ncls = fTRDncls;
    break;
  case 3:
    if (fTOFindex != -1)
      ncls = 1;
    break;
  case 4: //PHOS
    break;
  case 5: //HMPID
    if ((fHMPIDcluIdx >= 0) && (fHMPIDcluIdx < 7000000)) {
      if ((fHMPIDcluIdx%1000000 != 9999) && (fHMPIDcluIdx%1000000 != 99999)) {
	ncls = 1;
      }
    }    
    break;
  default:
    break;
  }
  return ncls;
}

Int_t AliESDtrack::GetClusters(Int_t idet, Int_t *idx) const
{
  // Get cluster index array by subdetector index
  //
  Int_t ncls = 0;
  switch(idet){
  case 0:
    ncls = GetITSclusters(idx);
    break;
  case 1:
    ncls = GetTPCclusters(idx);
    break;
  case 2:
    ncls = GetTRDclusters(idx);
    break;
  case 3:
    if (fTOFindex != -1) {
      idx[0] = fTOFindex;
      ncls = 1;
    }
    break;
  case 4: //PHOS
    break;
  case 5:
    if ((fHMPIDcluIdx >= 0) && (fHMPIDcluIdx < 7000000)) {
      if ((fHMPIDcluIdx%1000000 != 9999) && (fHMPIDcluIdx%1000000 != 99999)) {
	idx[0] = GetHMPIDcluIdx();
	ncls = 1;
      }
    }    
    break;
  case 6: //EMCAL
    break;
  default:
    break;
  }
  return ncls;
}

//_______________________________________________________________________
void AliESDtrack::GetIntegratedTimes(Double_t *times) const {
  // Returns the array with integrated times for each particle hypothesis
  for (Int_t i=0; i<AliPID::kSPECIES; i++) times[i]=fTrackTime[i];
}

//_______________________________________________________________________
void AliESDtrack::SetIntegratedTimes(const Double_t *times) {
  // Sets the array with integrated times for each particle hypotesis
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fTrackTime[i]=times[i];
}

//_______________________________________________________________________
void AliESDtrack::SetITSpid(const Double_t *p) {
  // Sets values for the probability of each particle type (in ITS)
  SetPIDValues(fITSr,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kITSpid);
}

//_______________________________________________________________________
void AliESDtrack::GetITSpid(Double_t *p) const {
  // Gets the probability of each particle type (in ITS)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fITSr[i];
}

//_______________________________________________________________________
Char_t AliESDtrack::GetITSclusters(Int_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  if (idx && fFriendTrack) {
    Int_t *index=fFriendTrack->GetITSindices();
    for (Int_t i=0; i<AliESDfriendTrack::kMaxITScluster; i++) {
      if ( (i>=fITSncls) && (i<6) ) idx[i]=-1;
      else {
	if (index) {
	  idx[i]=index[i];
	}
	else idx[i]= -2;
      }
    }
  }
  return fITSncls;
}

//_______________________________________________________________________
Bool_t AliESDtrack::GetITSModuleIndexInfo(Int_t ilayer,Int_t &idet,Int_t &status,
					 Float_t &xloc,Float_t &zloc) const {
  //----------------------------------------------------------------------
  // This function encodes in the module number also the status of cluster association
  // "status" can have the following values: 
  // 1 "found" (cluster is associated), 
  // 2 "dead" (module is dead from OCDB), 
  // 3 "skipped" (module or layer forced to be skipped),
  // 4 "outinz" (track out of z acceptance), 
  // 5 "nocls" (no clusters in the road), 
  // 6 "norefit" (cluster rejected during refit), 
  // 7 "deadzspd" (holes in z in SPD)
  // Also given are the coordinates of the crossing point of track and module
  // (in the local module ref. system)
  // WARNING: THIS METHOD HAS TO BE SYNCHRONIZED WITH AliITStrackV2::GetModuleIndexInfo()!
  //----------------------------------------------------------------------

  if(fITSModule[ilayer]==-1) {
    idet = -1;
    status=0;
    xloc=-99.; zloc=-99.;
    return kFALSE;
  }

  Int_t module = fITSModule[ilayer];

  idet = Int_t(module/1000000);

  module -= idet*1000000;

  status = Int_t(module/100000);

  module -= status*100000;

  Int_t signs = Int_t(module/10000);

  module-=signs*10000;

  Int_t xInt = Int_t(module/100);
  module -= xInt*100;

  Int_t zInt = module;

  if(signs==1) { xInt*=1; zInt*=1; }
  if(signs==2) { xInt*=1; zInt*=-1; }
  if(signs==3) { xInt*=-1; zInt*=1; }
  if(signs==4) { xInt*=-1; zInt*=-1; }

  xloc = 0.1*(Float_t)xInt;
  zloc = 0.1*(Float_t)zInt;

  if(status==4) idet = -1;

  return kTRUE;
}

//_______________________________________________________________________
UShort_t AliESDtrack::GetTPCclusters(Int_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  if (idx && fFriendTrack) {
    Int_t *index=fFriendTrack->GetTPCindices();

    if (index){
      for (Int_t i=0; i<AliESDfriendTrack::kMaxTPCcluster; i++) idx[i]=index[i];
    }
    else {
      for (Int_t i=0; i<AliESDfriendTrack::kMaxTPCcluster; i++) idx[i]=-2;
    }
  }
  return fTPCncls;
}

//_______________________________________________________________________
Float_t AliESDtrack::GetTPCCrossedRows() const
{
  // This function calls GetTPCClusterInfo with some default parameters which are used in the track selection and caches the outcome
  // because GetTPCClusterInfo is quite time-consuming
  
  if (fCacheNCrossedRows > -1)
    return fCacheNCrossedRows;
  
  fCacheNCrossedRows = GetTPCClusterInfo(2, 1);
  return fCacheNCrossedRows;
}

//_______________________________________________________________________
Float_t AliESDtrack::GetTPCClusterInfo(Int_t nNeighbours/*=3*/, Int_t type/*=0*/, Int_t row0, Int_t row1) const
{
  //
  // TPC cluster information
  // type 0: get fraction of found/findable clusters with neighbourhood definition
  //      1: findable clusters with neighbourhood definition
  //      2: found clusters
  //
  // definition of findable clusters:
  //            a cluster is defined as findable if there is another cluster
  //           within +- nNeighbours pad rows. The idea is to overcome threshold
  //           effects with a very simple algorithm.
  //

  if (type==2) return fTPCClusterMap.CountBits();
  
  Int_t found=0;
  Int_t findable=0;
  Int_t last=-nNeighbours;
  
  Int_t upperBound=fTPCClusterMap.GetNbits();
  if (upperBound>row1) upperBound=row1;
  for (Int_t i=row0; i<upperBound; ++i){
    //look to current row
    if (fTPCClusterMap[i]) {
      last=i;
      ++found;
      ++findable;
      continue;
    }
    //look to nNeighbours before
    if ((i-last)<=nNeighbours) {
      ++findable;
      continue;
    }
    //look to nNeighbours after
    for (Int_t j=i+1; j<i+1+nNeighbours; ++j){
      if (fTPCClusterMap[j]){
        ++findable;
        break;
      }
    }
  }
  if (type==1) return findable;
  
  if (type==0){
    Float_t fraction=0;
    if (findable>0) 
      fraction=(Float_t)found/(Float_t)findable;
    else 
      fraction=0;
    return fraction;
  }  
  return 0;  // undefined type - default value
}

//_______________________________________________________________________
Double_t AliESDtrack::GetTPCdensity(Int_t row0, Int_t row1) const{
  //
  // GetDensity of the clusters on given region between row0 and row1
  // Dead zone effect takin into acoount
  //
  if (!fFriendTrack) return 0.0;
  Int_t good  = 0;
  Int_t found = 0;
  //  
  Int_t *index=fFriendTrack->GetTPCindices();
  for (Int_t i=row0;i<=row1;i++){     
    Int_t idx = index[i];
    if (idx!=-1)  good++;             // track outside of dead zone
    if (idx>0)    found++;
  }
  Float_t density=0.5;
  if (good>(row1-row0)*0.5) density = Float_t(found)/Float_t(good);
  return density;
}

//_______________________________________________________________________
void AliESDtrack::SetTPCpid(const Double_t *p) {  
  // Sets values for the probability of each particle type (in TPC)
  SetPIDValues(fTPCr,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kTPCpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTPCpid(Double_t *p) const {
  // Gets the probability of each particle type (in TPC)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fTPCr[i];
}

//_______________________________________________________________________
UChar_t AliESDtrack::GetTRDclusters(Int_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined TRD clusters 
  //---------------------------------------------------------------------
  if (idx && fFriendTrack) {
    Int_t *index=fFriendTrack->GetTRDindices();

    if (index) {
      for (Int_t i=0; i<AliESDfriendTrack::kMaxTRDcluster; i++) idx[i]=index[i];
    }
    else {
      for (Int_t i=0; i<AliESDfriendTrack::kMaxTRDcluster; i++) idx[i]=-2;
    }
  }
  return fTRDncls;
}

//_______________________________________________________________________
UChar_t AliESDtrack::GetTRDtracklets(Int_t *idx) const {
//
// This function returns the number of TRD tracklets used in tracking
// and it fills the indices of these tracklets in the array "idx" as they 
// are registered in the TRD track list. 
// 
// Caution :
//   1. The idx array has to be allocated with a size >= AliESDtrack::kTRDnPlanes
//   2. The idx array store not only the index but also the layer of the tracklet. 
//      Therefore tracks with TRD gaps contain default values for indices [-1] 

  if (!fFriendTrack) return 0;
  if (!idx) return GetTRDntracklets();
  Int_t *index=fFriendTrack->GetTRDindices();
  Int_t n = 0;
  for (Int_t i=0; i<kTRDnPlanes; i++){ 
    if (index){
      if(index[i]>=0) n++;
      idx[i]=index[i];
    }
    else idx[i] = -2;
  }
  return n;
}

//_______________________________________________________________________
void AliESDtrack::SetTRDpid(const Double_t *p) {  
  // Sets values for the probability of each particle type (in TRD)
  SetPIDValues(fTRDr,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kTRDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTRDpid(Double_t *p) const {
  // Gets the probability of each particle type (in TRD)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fTRDr[i];
}

//_______________________________________________________________________
void    AliESDtrack::SetTRDpid(Int_t iSpecies, Float_t p)
{
  // Sets the probability of particle type iSpecies to p (in TRD)
  fTRDr[iSpecies] = p;
}

Double_t AliESDtrack::GetTRDpid(Int_t iSpecies) const
{
  // Returns the probability of particle type iSpecies (in TRD)
  return fTRDr[iSpecies];
}

//____________________________________________________
Int_t AliESDtrack::GetNumberOfTRDslices() const 
{
  // built in backward compatibility
  Int_t idx = fTRDnSlices - (kTRDnPlanes<<1);
  return idx<18 ? fTRDnSlices/kTRDnPlanes : idx/kTRDnPlanes;
}

//____________________________________________________
Double_t AliESDtrack::GetTRDmomentum(Int_t plane, Double_t *sp) const
{
//Returns momentum estimation and optional its error (sp)
// in TRD layer "plane".

  if (!fTRDnSlices) {
    AliDebug(2, "No TRD info allocated for this track.");
    return -1.;
  }
  if ((plane<0) || (plane>=kTRDnPlanes)) {
    AliWarning(Form("Request for TRD plane[%d] outside range.", plane)); 
    return -1.;
  }

  Int_t idx = fTRDnSlices-(kTRDnPlanes<<1)+plane;
  // Protection for backward compatibility
  if(idx<(GetNumberOfTRDslices()*kTRDnPlanes)) return -1.;

  if(sp) (*sp) = fTRDslices[idx+kTRDnPlanes];
  return fTRDslices[idx];
}

//____________________________________________________
Double_t  AliESDtrack::GetTRDslice(Int_t plane, Int_t slice) const {
  //Gets the charge from the slice of the plane

  if(!fTRDslices) {
    //AliError("No TRD slices allocated for this track !");
    return -1.;
  }
  if ((plane<0) || (plane>=kTRDnPlanes)) {
    AliError("Info for TRD plane not available !");
    return -1.;
  }
  Int_t ns=GetNumberOfTRDslices();
  if ((slice<-1) || (slice>=ns)) {
    //AliError("Wrong TRD slice !");  
    return -1.;
  }

  if(slice>=0) return fTRDslices[plane*ns + slice];

  // return average of the dEdx measurements
  Double_t q=0.; Double32_t *s = &fTRDslices[plane*ns];
  for (Int_t i=0; i<ns; i++, s++) if((*s)>0.) q+=(*s);
  return q/ns;
}

//____________________________________________________
void  AliESDtrack::SetNumberOfTRDslices(Int_t n) {
  //Sets the number of slices used for PID 
  if (fTRDnSlices) return;

  fTRDnSlices=n;
  fTRDslices=new Double32_t[fTRDnSlices];
  
  // set-up correctly the allocated memory
  memset(fTRDslices, 0, n*sizeof(Double32_t));
  for (Int_t i=GetNumberOfTRDslices(); i--;) fTRDslices[i]=-1.;
}

//____________________________________________________
void  AliESDtrack::SetTRDslice(Double_t q, Int_t plane, Int_t slice) {
  //Sets the charge q in the slice of the plane
  if(!fTRDslices) {
    AliError("No TRD slices allocated for this track !");
    return;
  }
  if ((plane<0) || (plane>=kTRDnPlanes)) {
    AliError("Info for TRD plane not allocated !");
    return;
  }
  Int_t ns=GetNumberOfTRDslices();
  if ((slice<0) || (slice>=ns)) {
    AliError("Wrong TRD slice !");
    return;
  }
  Int_t n=plane*ns + slice;
  fTRDslices[n]=q;
}


//____________________________________________________
void AliESDtrack::SetTRDmomentum(Double_t p, Int_t plane, Double_t *sp)
{
  if(!fTRDslices) {
    AliError("No TRD slices allocated for this track !");
    return;
  }
  if ((plane<0) || (plane>=kTRDnPlanes)) {
    AliError("Info for TRD plane not allocated !");
    return;
  }

  Int_t idx = fTRDnSlices-(kTRDnPlanes<<1)+plane;
  // Protection for backward compatibility
  if(idx<GetNumberOfTRDslices()*kTRDnPlanes) return;

  if(sp) fTRDslices[idx+kTRDnPlanes] = (*sp);
  fTRDslices[idx] = p;
}


//_______________________________________________________________________
void AliESDtrack::SetTOFpid(const Double_t *p) {  
  // Sets the probability of each particle type (in TOF)
  SetPIDValues(fTOFr,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kTOFpid);
}

//_______________________________________________________________________
void AliESDtrack::SetTOFLabel(const Int_t *p) {  
  // Sets  (in TOF)
  for (Int_t i=0; i<3; i++) fTOFLabel[i]=p[i];
}

//_______________________________________________________________________
void AliESDtrack::GetTOFpid(Double_t *p) const {
  // Gets probabilities of each particle type (in TOF)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fTOFr[i];
}

//_______________________________________________________________________
void AliESDtrack::GetTOFLabel(Int_t *p) const {
  // Gets (in TOF)
  for (Int_t i=0; i<3; i++) p[i]=fTOFLabel[i];
}

//_______________________________________________________________________
void AliESDtrack::GetTOFInfo(Float_t *info) const {
  // Gets (in TOF)
  for (Int_t i=0; i<10; i++) info[i]=fTOFInfo[i];
}

//_______________________________________________________________________
void AliESDtrack::SetTOFInfo(Float_t*info) {
  // Gets (in TOF)
  for (Int_t i=0; i<10; i++) fTOFInfo[i]=info[i];
}



//_______________________________________________________________________
void AliESDtrack::SetHMPIDpid(const Double_t *p) {  
  // Sets the probability of each particle type (in HMPID)
  SetPIDValues(fHMPIDr,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kHMPIDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetHMPIDpid(Double_t *p) const {
  // Gets probabilities of each particle type (in HMPID)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fHMPIDr[i];
}



//_______________________________________________________________________
void AliESDtrack::SetESDpid(const Double_t *p) {  
  // Sets the probability of each particle type for the ESD track
  SetPIDValues(fR,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kESDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetESDpid(Double_t *p) const {
  // Gets probability of each particle type for the ESD track
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fR[i];
}

//_______________________________________________________________________
Bool_t AliESDtrack::RelateToVertexTPC(const AliESDVertex *vtx, 
Double_t b, Double_t maxd, AliExternalTrackParam *cParam) {
  //
  // Try to relate the TPC-only track parameters to the vertex "vtx", 
  // if the (rough) transverse impact parameter is not bigger then "maxd". 
  //            Magnetic field is "b" (kG).
  //
  // a) The TPC-only paramters are extapolated to the DCA to the vertex.
  // b) The impact parameters and their covariance matrix are calculated.
  // c) An attempt to constrain the TPC-only params to the vertex is done.
  //    The constrained params are returned via "cParam".
  //
  // In the case of success, the returned value is kTRUE
  // otherwise, it's kFALSE)
  // 

  if (!fTPCInner) return kFALSE;
  if (!vtx) return kFALSE;

  Double_t dz[2],cov[3];
  if (!fTPCInner->PropagateToDCA(vtx, b, maxd, dz, cov)) return kFALSE;

  fdTPC = dz[0];
  fzTPC = dz[1];  
  fCddTPC = cov[0];
  fCdzTPC = cov[1];
  fCzzTPC = cov[2];
  
  Double_t covar[6]; vtx->GetCovMatrix(covar);
  Double_t p[2]={GetParameter()[0]-dz[0],GetParameter()[1]-dz[1]};
  Double_t c[3]={covar[2],0.,covar[5]};

  Double_t chi2=GetPredictedChi2(p,c);
  if (chi2>kVeryBig) return kFALSE;

  fCchi2TPC=chi2;

  if (!cParam) return kTRUE;

  *cParam = *fTPCInner;
  if (!cParam->Update(p,c)) return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
Bool_t AliESDtrack::RelateToVertexTPCBxByBz(const AliESDVertex *vtx, 
Double_t b[3], Double_t maxd, AliExternalTrackParam *cParam) {
  //
  // Try to relate the TPC-only track parameters to the vertex "vtx", 
  // if the (rough) transverse impact parameter is not bigger then "maxd". 
  //
  // All three components of the magnetic field ,"b[3]" (kG), 
  // are taken into account.
  //
  // a) The TPC-only paramters are extapolated to the DCA to the vertex.
  // b) The impact parameters and their covariance matrix are calculated.
  // c) An attempt to constrain the TPC-only params to the vertex is done.
  //    The constrained params are returned via "cParam".
  //
  // In the case of success, the returned value is kTRUE
  // otherwise, it's kFALSE)
  // 

  if (!fTPCInner) return kFALSE;
  if (!vtx) return kFALSE;

  Double_t dz[2],cov[3];
  if (!fTPCInner->PropagateToDCABxByBz(vtx, b, maxd, dz, cov)) return kFALSE;

  fdTPC = dz[0];
  fzTPC = dz[1];  
  fCddTPC = cov[0];
  fCdzTPC = cov[1];
  fCzzTPC = cov[2];
  
  Double_t covar[6]; vtx->GetCovMatrix(covar);
  Double_t p[2]={GetParameter()[0]-dz[0],GetParameter()[1]-dz[1]};
  Double_t c[3]={covar[2],0.,covar[5]};

  Double_t chi2=GetPredictedChi2(p,c);
  if (chi2>kVeryBig) return kFALSE;

  fCchi2TPC=chi2;

  if (!cParam) return kTRUE;

  *cParam = *fTPCInner;
  if (!cParam->Update(p,c)) return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
Bool_t AliESDtrack::RelateToVertex(const AliESDVertex *vtx, 
Double_t b, Double_t maxd, AliExternalTrackParam *cParam) {
  //
  // Try to relate this track to the vertex "vtx", 
  // if the (rough) transverse impact parameter is not bigger then "maxd". 
  //            Magnetic field is "b" (kG).
  //
  // a) The track gets extapolated to the DCA to the vertex.
  // b) The impact parameters and their covariance matrix are calculated.
  // c) An attempt to constrain this track to the vertex is done.
  //    The constrained params are returned via "cParam".
  //
  // In the case of success, the returned value is kTRUE
  // (otherwise, it's kFALSE)
  //  

  if (!vtx) return kFALSE;

  Double_t dz[2],cov[3];
  if (!PropagateToDCA(vtx, b, maxd, dz, cov)) return kFALSE;

  fD = dz[0];
  fZ = dz[1];  
  fCdd = cov[0];
  fCdz = cov[1];
  fCzz = cov[2];
  
  Double_t covar[6]; vtx->GetCovMatrix(covar);
  Double_t p[2]={GetParameter()[0]-dz[0],GetParameter()[1]-dz[1]};
  Double_t c[3]={covar[2],0.,covar[5]};

  Double_t chi2=GetPredictedChi2(p,c);
  if (chi2>kVeryBig) return kFALSE;

  fCchi2=chi2;


  //--- Could now these lines be removed ? ---
  delete fCp;
  fCp=new AliExternalTrackParam(*this);  

  if (!fCp->Update(p,c)) {delete fCp; fCp=0; return kFALSE;}
  //----------------------------------------

  fVertexID = vtx->GetID();

  if (!cParam) return kTRUE;

  *cParam = *this;
  if (!cParam->Update(p,c)) return kFALSE; 

  return kTRUE;
}

//_______________________________________________________________________
Bool_t AliESDtrack::RelateToVertexBxByBz(const AliESDVertex *vtx, 
Double_t b[3], Double_t maxd, AliExternalTrackParam *cParam) {
  //
  // Try to relate this track to the vertex "vtx", 
  // if the (rough) transverse impact parameter is not bigger then "maxd". 
  //            Magnetic field is "b" (kG).
  //
  // a) The track gets extapolated to the DCA to the vertex.
  // b) The impact parameters and their covariance matrix are calculated.
  // c) An attempt to constrain this track to the vertex is done.
  //    The constrained params are returned via "cParam".
  //
  // In the case of success, the returned value is kTRUE
  // (otherwise, it's kFALSE)
  //  

  if (!vtx) return kFALSE;

  Double_t dz[2],cov[3];
  if (!PropagateToDCABxByBz(vtx, b, maxd, dz, cov)) return kFALSE;

  fD = dz[0];
  fZ = dz[1];  
  fCdd = cov[0];
  fCdz = cov[1];
  fCzz = cov[2];
  
  Double_t covar[6]; vtx->GetCovMatrix(covar);
  Double_t p[2]={GetParameter()[0]-dz[0],GetParameter()[1]-dz[1]};
  Double_t c[3]={covar[2],0.,covar[5]};

  Double_t chi2=GetPredictedChi2(p,c);
  if (chi2>kVeryBig) return kFALSE;

  fCchi2=chi2;


  //--- Could now these lines be removed ? ---
  delete fCp;
  fCp=new AliExternalTrackParam(*this);  

  if (!fCp->Update(p,c)) {delete fCp; fCp=0; return kFALSE;}
  //----------------------------------------

  fVertexID = vtx->GetID();

  if (!cParam) return kTRUE;

  *cParam = *this;
  if (!cParam->Update(p,c)) return kFALSE; 

  return kTRUE;
}

//_______________________________________________________________________
void AliESDtrack::Print(Option_t *) const {
  // Prints info on the track
  AliExternalTrackParam::Print();
  printf("ESD track info\n") ; 
  Double_t p[AliPID::kSPECIESN] ; 
  Int_t index = 0 ; 
  if( IsOn(kITSpid) ){
    printf("From ITS: ") ; 
    GetITSpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetITSsignal()) ;
  } 
  if( IsOn(kTPCpid) ){
    printf("From TPC: ") ; 
    GetTPCpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetTPCsignal()) ;
  }
  if( IsOn(kTRDpid) ){
    printf("From TRD: ") ; 
    GetTRDpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
      printf("\n           signal = %f\n", GetTRDsignal()) ;
  }
  if( IsOn(kTOFpid) ){
    printf("From TOF: ") ; 
    GetTOFpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetTOFsignal()) ;
  }
  if( IsOn(kHMPIDpid) ){
    printf("From HMPID: ") ; 
    GetHMPIDpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetHMPIDsignal()) ;
  }
} 


//
// Draw functionality
// Origin: Marian Ivanov, Marian.Ivanov@cern.ch
//
void AliESDtrack::FillPolymarker(TPolyMarker3D *pol, Float_t magF, Float_t minR, Float_t maxR, Float_t stepR){
  //
  // Fill points in the polymarker
  //
  TObjArray arrayRef;
  arrayRef.AddLast(new AliExternalTrackParam(*this));
  if (fIp) arrayRef.AddLast(new AliExternalTrackParam(*fIp));
  if (fOp) arrayRef.AddLast(new AliExternalTrackParam(*fOp));
  if (fHMPIDp) arrayRef.AddLast(new AliExternalTrackParam(*fHMPIDp));
  //
  Double_t mpos[3]={0,0,0};
  Int_t entries=arrayRef.GetEntries();
  for (Int_t i=0;i<entries;i++){
    Double_t pos[3];
    ((AliExternalTrackParam*)arrayRef.At(i))->GetXYZ(pos);
    mpos[0]+=pos[0]/entries;
    mpos[1]+=pos[1]/entries;
    mpos[2]+=pos[2]/entries;    
  }
  // Rotate to the mean position
  //
  Float_t fi= TMath::ATan2(mpos[1],mpos[0]);
  for (Int_t i=0;i<entries;i++){
    Bool_t res = ((AliExternalTrackParam*)arrayRef.At(i))->Rotate(fi);
    if (!res) delete arrayRef.RemoveAt(i);
  }
  Int_t counter=0;
  for (Double_t r=minR; r<maxR; r+=stepR){
    Double_t sweight=0;
    Double_t mlpos[3]={0,0,0};
    for (Int_t i=0;i<entries;i++){
      Double_t point[3]={0,0,0};
      AliExternalTrackParam *param = ((AliExternalTrackParam*)arrayRef.At(i));
      if (!param) continue;
      if (param->GetXYZAt(r,magF,point)){
	Double_t weight = 1./(10.+(r-param->GetX())*(r-param->GetX()));
	sweight+=weight;
	mlpos[0]+=point[0]*weight;
	mlpos[1]+=point[1]*weight;
	mlpos[2]+=point[2]*weight;
      }
    }
    if (sweight>0){
      mlpos[0]/=sweight;
      mlpos[1]/=sweight;
      mlpos[2]/=sweight;      
      pol->SetPoint(counter,mlpos[0],mlpos[1], mlpos[2]);
      //      printf("xyz\t%f\t%f\t%f\n",mlpos[0], mlpos[1],mlpos[2]);
      counter++;
    }
  }
}

//_______________________________________________________________________
void AliESDtrack::SetITSdEdxSamples(const Double_t s[4]) {
  //
  // Store the dE/dx samples measured by the two SSD and two SDD layers.
  // These samples are corrected for the track segment length. 
  //
  for (Int_t i=0; i<4; i++) fITSdEdxSamples[i]=s[i];
}

//_______________________________________________________________________
void AliESDtrack::GetITSdEdxSamples(Double_t *s) const {
  //
  // Get the dE/dx samples measured by the two SSD and two SDD layers.  
  // These samples are corrected for the track segment length.
  //
  for (Int_t i=0; i<4; i++) s[i]=fITSdEdxSamples[i];
}


UShort_t   AliESDtrack::GetTPCnclsS(Int_t i0,Int_t i1) const{
  //
  // get number of shared TPC clusters
  //
  return  fTPCSharedMap.CountBits(i0)-fTPCSharedMap.CountBits(i1);
}

UShort_t   AliESDtrack::GetTPCncls(Int_t i0,Int_t i1) const{
  //
  // get number of TPC clusters
  //
  return  fTPCClusterMap.CountBits(i0)-fTPCClusterMap.CountBits(i1);
}

//____________________________________________________________________
Double_t AliESDtrack::GetChi2TPCConstrainedVsGlobal(const AliESDVertex* vtx) const
{
  // Calculates the chi2 between the TPC track (TPCinner) constrained to the primary vertex and the global track
  //
  // Returns -1 in case the calculation failed
  //
  // Value is cached as a non-persistent member.
  //
  // Code adapted from original code by GSI group (Jacek, Marian, Michael)
  
  // cache, ignoring that a different vertex might be passed
  if (fCacheChi2TPCConstrainedVsGlobalVertex == vtx)
    return fCacheChi2TPCConstrainedVsGlobal;
  
  fCacheChi2TPCConstrainedVsGlobal = -1;
  fCacheChi2TPCConstrainedVsGlobalVertex = vtx;
  
  Double_t x[3];
  GetXYZ(x);
  Double_t b[3];
  AliTrackerBase::GetBxByBz(x,b);

  if (!fTPCInner)  { 
    AliWarning("Could not get TPC Inner Param.");
    return fCacheChi2TPCConstrainedVsGlobal;
  }
  
  // clone for constraining
  AliExternalTrackParam* tpcInnerC = new AliExternalTrackParam(*fTPCInner);
  if (!tpcInnerC) { 
    AliWarning("Clone of TPCInnerParam failed.");
    return fCacheChi2TPCConstrainedVsGlobal;  
  }
  
  // transform to the track reference frame 
  Bool_t isOK = tpcInnerC->Rotate(GetAlpha());
  isOK &= tpcInnerC->PropagateTo(GetX(), b[2]);
  if (!isOK) { 
    delete tpcInnerC;
    tpcInnerC = 0; 
    AliWarning("Rotation/Propagation of track failed.") ; 
    return fCacheChi2TPCConstrainedVsGlobal;    
  }  

  // constrain TPCinner 
  isOK = tpcInnerC->ConstrainToVertex(vtx, b);
  
  // transform to the track reference frame 
  isOK &= tpcInnerC->Rotate(GetAlpha());
  isOK &= tpcInnerC->PropagateTo(GetX(), b[2]);

  if (!isOK) {
    AliWarning("ConstrainTPCInner failed.") ;
    delete tpcInnerC;
    tpcInnerC = 0; 
    return fCacheChi2TPCConstrainedVsGlobal;  
  }
  
  // calculate chi2 between vi and vj vectors
  // with covi and covj covariance matrices
  // chi2ij = (vi-vj)^(T)*(covi+covj)^(-1)*(vi-vj)
  TMatrixD deltaT(5,1);
  TMatrixD delta(1,5);
  TMatrixD covarM(5,5);

  for (Int_t ipar=0; ipar<5; ipar++) {
    deltaT(ipar,0) = tpcInnerC->GetParameter()[ipar] - GetParameter()[ipar];
    delta(0,ipar) = tpcInnerC->GetParameter()[ipar] - GetParameter()[ipar];

    for (Int_t jpar=0; jpar<5; jpar++) {
      Int_t index = GetIndex(ipar,jpar);
      covarM(ipar,jpar) = GetCovariance()[index]+tpcInnerC->GetCovariance()[index];
    }
  }
  // chi2 distance TPC constrained and TPC+ITS
  TMatrixD covarMInv = covarM.Invert();
  TMatrixD mat2 = covarMInv*deltaT;
  TMatrixD chi2 = delta*mat2; 
  
  delete tpcInnerC; 
  tpcInnerC = 0;
  
  fCacheChi2TPCConstrainedVsGlobal = chi2(0,0);
  return fCacheChi2TPCConstrainedVsGlobal;
}
