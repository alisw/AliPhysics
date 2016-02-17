/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD track implementation of AliVTrack
//     Author: Markus Oldenburg, CERN
//     Markus.Oldenburg@cern.ch
//-------------------------------------------------------------------------

#include <TVector3.h>
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliVVertex.h"
#include "AliDetectorPID.h"
#include "AliAODEvent.h"
#include "AliAODHMPIDrings.h"
#include "AliTOFHeader.h"

#include "AliAODTrack.h"

ClassImp(AliAODTrack)

//______________________________________________________________________________
AliAODTrack::AliAODTrack() : 
  AliVTrack(),
  fRAtAbsorberEnd(0.),
  fChi2perNDF(-999.),
  fChi2MatchTrigger(0.),
  fPID(0),
  fITSchi2(0),
  fFlags(0),
  fLabel(-999),
  fTOFLabel(),
  fTrackLength(0),
  fITSMuonClusterMap(0),
  fMUONtrigHitsMapTrg(0),
  fMUONtrigHitsMapTrk(0),
  fFilterMap(0),
  fTPCFitMap(),
  fTPCClusterMap(),
  fTPCSharedMap(),
  fTPCnclsF(0),
  fTPCNCrossedRows(0),
  fID(-999),
  fCharge(-99),
  fType(kUndef),
  fPIDForTracking(AliPID::kPion),
  fCaloIndex(kEMCALNoMatch),
  fCovMatrix(NULL),
  fDetPid(NULL),
  fDetectorPID(NULL),
  fProdVertex(NULL),
  fTrackPhiOnEMCal(-999),
  fTrackEtaOnEMCal(-999),
  fTrackPtOnEMCal(-999),
  fIsMuonGlobalTrack(kFALSE),    // AU
  fITSsignalTuned(0.),
  fTPCsignalTuned(0),
  fTOFsignalTuned(99999),
  fMFTClusterPattern(0),         // AU
  fAODEvent(NULL)
{
  // default constructor

  SetP();
  SetPosition((Float_t*)NULL);
  SetXYAtDCA(-999., -999.);
  SetPxPyPzAtDCA(-999., -999., -999.);
  for (Int_t i = 0; i < 3; i++) {fTOFLabel[i] = -1;}
}

//______________________________________________________________________________
AliAODTrack::AliAODTrack(Short_t id,
			 Int_t label, 
			 Double_t p[3],
			 Bool_t cartesian,
			 Double_t x[3],
			 Bool_t isDCA,
			 Double_t covMatrix[21],
			 Short_t charge,
			 UChar_t itsClusMap,
			 AliAODVertex *prodVertex,
			 Bool_t usedForVtxFit,
			 Bool_t usedForPrimVtxFit,
			 AODTrk_t ttype,
			 UInt_t selectInfo,
			 Float_t chi2perNDF) :
  AliVTrack(),
  fRAtAbsorberEnd(0.),
  fChi2perNDF(chi2perNDF),
  fChi2MatchTrigger(0.),
  fPID(0),
  fITSchi2(0),
  fFlags(0),
  fLabel(label),
  fTOFLabel(),
  fTrackLength(0),
  fITSMuonClusterMap(0),
  fMUONtrigHitsMapTrg(0),
  fMUONtrigHitsMapTrk(0),
  fFilterMap(selectInfo),
  fTPCFitMap(),
  fTPCClusterMap(),
  fTPCSharedMap(),
  fTPCnclsF(0),
  fTPCNCrossedRows(0),
  fID(id),
  fCharge(charge),
  fType(ttype),
  fPIDForTracking(AliPID::kPion),
  fCaloIndex(kEMCALNoMatch),
  fCovMatrix(NULL),
  fDetPid(NULL),
  fDetectorPID(NULL),
  fProdVertex(prodVertex),
  fTrackPhiOnEMCal(-999),
  fTrackEtaOnEMCal(-999),
  fTrackPtOnEMCal(-999),
  fIsMuonGlobalTrack(kFALSE),    // AU
  fITSsignalTuned(0),
  fTPCsignalTuned(0),
  fTOFsignalTuned(99999),
  fMFTClusterPattern(0),         // AU
  fAODEvent(NULL)
{
  // constructor
 
  SetP(p, cartesian);
  SetPosition(x, isDCA);
  SetXYAtDCA(-999., -999.);
  SetPxPyPzAtDCA(-999., -999., -999.);
  SetUsedForVtxFit(usedForVtxFit);
  SetUsedForPrimVtxFit(usedForPrimVtxFit);
  if(covMatrix) SetCovMatrix(covMatrix);
  SetITSClusterMap(itsClusMap);
  for (Int_t i=0;i<3;i++) {fTOFLabel[i]=-1;}
}

//______________________________________________________________________________
AliAODTrack::AliAODTrack(Short_t id,
			 Int_t label, 
			 Float_t p[3],
			 Bool_t cartesian,
			 Float_t x[3],
			 Bool_t isDCA,
			 Float_t covMatrix[21],
			 Short_t charge,
			 UChar_t itsClusMap,
			 AliAODVertex *prodVertex,
			 Bool_t usedForVtxFit,
			 Bool_t usedForPrimVtxFit,
			 AODTrk_t ttype,
			 UInt_t selectInfo,
			 Float_t chi2perNDF ) :
  AliVTrack(),
  fRAtAbsorberEnd(0.),
  fChi2perNDF(chi2perNDF),
  fChi2MatchTrigger(0.),
  fPID(0),
  fITSchi2(0),
  fFlags(0),
  fLabel(label),
  fTOFLabel(),
  fTrackLength(0),
  fITSMuonClusterMap(0),
  fMUONtrigHitsMapTrg(0),
  fMUONtrigHitsMapTrk(0),
  fFilterMap(selectInfo),
  fTPCFitMap(),
  fTPCClusterMap(),
  fTPCSharedMap(),
  fTPCnclsF(0),
  fTPCNCrossedRows(0),
  fID(id),
  fCharge(charge),
  fType(ttype),
  fPIDForTracking(AliPID::kPion),
  fCaloIndex(kEMCALNoMatch),
  fCovMatrix(NULL),
  fDetPid(NULL),
  fDetectorPID(NULL),
  fProdVertex(prodVertex),
  fTrackPhiOnEMCal(-999),
  fTrackEtaOnEMCal(-999),
  fTrackPtOnEMCal(-999),
  fIsMuonGlobalTrack(kFALSE),    // AU
  fITSsignalTuned(0),
  fTPCsignalTuned(0),
  fTOFsignalTuned(99999),
  fMFTClusterPattern(0),         // AU
  fAODEvent(NULL)
{
  // constructor
 
  SetP(p, cartesian);
  SetPosition(x, isDCA);
  SetXYAtDCA(-999., -999.);
  SetPxPyPzAtDCA(-999., -999., -999.);
  SetUsedForVtxFit(usedForVtxFit);
  SetUsedForPrimVtxFit(usedForPrimVtxFit);
  if(covMatrix) SetCovMatrix(covMatrix);
  SetITSClusterMap(itsClusMap);
  for (Int_t i=0;i<3;i++) {fTOFLabel[i]=-1;}
}

//______________________________________________________________________________
AliAODTrack::~AliAODTrack() 
{
  // destructor
  delete fCovMatrix;
  delete fDetPid;
  delete fDetectorPID;
  if (fPID) {delete[] fPID; fPID = 0;}
}


//______________________________________________________________________________
AliAODTrack::AliAODTrack(const AliAODTrack& trk) :
  AliVTrack(trk),
  fRAtAbsorberEnd(trk.fRAtAbsorberEnd),
  fChi2perNDF(trk.fChi2perNDF),
  fChi2MatchTrigger(trk.fChi2MatchTrigger),
  fPID(0),
  fITSchi2(trk.fITSchi2),
  fFlags(trk.fFlags),
  fLabel(trk.fLabel),
  fTOFLabel(),
  fTrackLength(trk.fTrackLength),
  fITSMuonClusterMap(trk.fITSMuonClusterMap),
  fMUONtrigHitsMapTrg(trk.fMUONtrigHitsMapTrg),
  fMUONtrigHitsMapTrk(trk.fMUONtrigHitsMapTrk),
  fFilterMap(trk.fFilterMap),
  fTPCFitMap(trk.fTPCFitMap),
  fTPCClusterMap(trk.fTPCClusterMap),
  fTPCSharedMap(trk.fTPCSharedMap),
  fTPCnclsF(trk.fTPCnclsF),
  fTPCNCrossedRows(trk.fTPCNCrossedRows),
  fID(trk.fID),
  fCharge(trk.fCharge),
  fType(trk.fType),
  fPIDForTracking(trk.fPIDForTracking),
  fCaloIndex(trk.fCaloIndex),
  fCovMatrix(NULL),
  fDetPid(NULL),
  fDetectorPID(NULL),
  fProdVertex(trk.fProdVertex),
  fTrackPhiOnEMCal(trk.fTrackPhiOnEMCal),
  fTrackEtaOnEMCal(trk.fTrackEtaOnEMCal),
  fTrackPtOnEMCal(trk.fTrackPtOnEMCal),
  fIsMuonGlobalTrack(trk.fIsMuonGlobalTrack),    // AU
  fITSsignalTuned(trk.fITSsignalTuned),
  fTPCsignalTuned(trk.fTPCsignalTuned),
  fTOFsignalTuned(trk.fTOFsignalTuned),
  fMFTClusterPattern(trk.fMFTClusterPattern),    // AU
  fAODEvent(trk.fAODEvent)
{
  // Copy constructor

  trk.GetP(fMomentum);
  trk.GetPosition(fPosition);
  SetXYAtDCA(trk.XAtDCA(), trk.YAtDCA());
  SetPxPyPzAtDCA(trk.PxAtDCA(), trk.PyAtDCA(), trk.PzAtDCA());
  SetUsedForVtxFit(trk.GetUsedForVtxFit());
  SetUsedForPrimVtxFit(trk.GetUsedForPrimVtxFit());
  if(trk.fCovMatrix) fCovMatrix=new AliAODRedCov<6>(*trk.fCovMatrix);
  if(trk.fDetPid) fDetPid=new AliAODPid(*trk.fDetPid);
  SetPID(trk.fPID);
  if (trk.fDetectorPID) fDetectorPID = new AliDetectorPID(*trk.fDetectorPID);
  for (Int_t i = 0; i < 3; i++) {fTOFLabel[i] = trk.fTOFLabel[i];}  
}

//______________________________________________________________________________
AliAODTrack& AliAODTrack::operator=(const AliAODTrack& trk)
{
  // Assignment operator
  if(this!=&trk) {

    AliVTrack::operator=(trk);

    trk.GetP(fMomentum);
    trk.GetPosition(fPosition);
    SetXYAtDCA(trk.XAtDCA(), trk.YAtDCA());
    SetPxPyPzAtDCA(trk.PxAtDCA(), trk.PyAtDCA(), trk.PzAtDCA());
    fRAtAbsorberEnd    = trk.fRAtAbsorberEnd;
    fChi2perNDF        = trk.fChi2perNDF;
    fChi2MatchTrigger  = trk.fChi2MatchTrigger;
    SetPID( trk.fPID );
    fITSchi2           = trk.fITSchi2;
    fFlags             = trk.fFlags;
    fLabel             = trk.fLabel;    
    fTrackLength       = trk.fTrackLength;
    fITSMuonClusterMap = trk.fITSMuonClusterMap;
    fMUONtrigHitsMapTrg = trk.fMUONtrigHitsMapTrg;
    fMUONtrigHitsMapTrk = trk.fMUONtrigHitsMapTrk;
    fFilterMap         = trk.fFilterMap;
    fTPCFitMap         = trk.fTPCFitMap;
    fTPCClusterMap     = trk.fTPCClusterMap;
    fTPCSharedMap      = trk.fTPCSharedMap;
    fTPCnclsF          = trk.fTPCnclsF;
    fTPCNCrossedRows   = trk.fTPCNCrossedRows;
    fID                = trk.fID;
    fCharge            = trk.fCharge;
    fType              = trk.fType;
    fPIDForTracking    = trk.fPIDForTracking;
    fCaloIndex         = trk.fCaloIndex;
    fTrackPhiOnEMCal   = trk.fTrackPhiOnEMCal;
    fTrackEtaOnEMCal   = trk.fTrackEtaOnEMCal;
    fTrackPtOnEMCal    = trk.fTrackPtOnEMCal;
    fIsMuonGlobalTrack = trk.fIsMuonGlobalTrack;     // AU
    fITSsignalTuned    = trk.fITSsignalTuned;
    fTPCsignalTuned    = trk.fTPCsignalTuned;
    fTOFsignalTuned    = trk.fTOFsignalTuned;
    fMFTClusterPattern = trk.fMFTClusterPattern;     // AU
    
    delete fCovMatrix;
    if(trk.fCovMatrix) fCovMatrix=new AliAODRedCov<6>(*trk.fCovMatrix);
    else fCovMatrix=NULL;


    fProdVertex        = trk.fProdVertex;
    SetUsedForVtxFit(trk.GetUsedForVtxFit());
    SetUsedForPrimVtxFit(trk.GetUsedForPrimVtxFit());

    //detector raw signals
    delete fDetPid;
    if(trk.fDetPid) fDetPid=new AliAODPid(*trk.fDetPid);
    else fDetPid=NULL;

    //calibrated PID cache
    delete fDetectorPID;
    fDetectorPID=0x0;
    if (trk.fDetectorPID) fDetectorPID = new AliDetectorPID(*trk.fDetectorPID);
    for (Int_t i = 0; i < 3; i++) {fTOFLabel[i] = trk.fTOFLabel[i];}  
  }

  return *this;
}

//______________________________________________________________________________
Double_t AliAODTrack::M(AODTrkPID_t pid) const
{
  // Returns the mass.
  // Masses for nuclei don't exist in the PDG tables, therefore they were put by hand.

  switch (pid) {

  case kElectron :
    return 0.000510999; //TDatabasePDG::Instance()->GetParticle(11/*::kElectron*/)->Mass();
    break;

  case kMuon :
    return 0.1056584; //TDatabasePDG::Instance()->GetParticle(13/*::kMuonMinus*/)->Mass();
    break;

  case kPion :
    return 0.13957; //TDatabasePDG::Instance()->GetParticle(211/*::kPiPlus*/)->Mass();
    break;

  case kKaon :
    return 0.4937; //TDatabasePDG::Instance()->GetParticle(321/*::kKPlus*/)->Mass();
    break;

  case kProton :
    return 0.9382720; //TDatabasePDG::Instance()->GetParticle(2212/*::kProton*/)->Mass();
    break;

  case kDeuteron :
    return 1.8756; //TDatabasePDG::Instance()->GetParticle(1000010020)->Mass();
    break;

  case kTriton :
    return 2.8089; //TDatabasePDG::Instance()->GetParticle(1000010030)->Mass();
    break;

  case kHelium3 :
    return 2.8084; //TDatabasePDG::Instance()->GetParticle(1000020030)->Mass();
    break;

  case kAlpha :
    return 3.7274; //TDatabasePDG::Instance()->GetParticle(1000020040)->Mass();
    break;

  case kUnknown :
    return -999.;
    break;

  default :
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliAODTrack::E(AODTrkPID_t pid) const
{
  // Returns the energy of the particle of a given pid.
  
  if (pid != kUnknown) { // particle was identified
    Double_t m = M(pid);
    return TMath::Sqrt(P()*P() + m*m);
  } else { // pid unknown
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliAODTrack::Y(AODTrkPID_t pid) const
{
  // Returns the rapidity of a particle of a given pid.
  
  if (pid != kUnknown) { // particle was identified
    Double_t e = E(pid);
    Double_t pz = Pz();
    if (e>=0 && e!=pz) { // energy was positive (e.g. not -999.) and not equal to pz
      return 0.5*TMath::Log((e+pz)/(e-pz));
    } else { // energy not known or equal to pz
      return -999.;
    }
  } else { // pid unknown
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliAODTrack::Y(Double_t m) const
{
  // Returns the rapidity of a particle of a given mass.
  
  if (m >= 0.) { // mass makes sense
    Double_t e = E(m);
    Double_t pz = Pz();
    if (e>=0 && e!=pz) { // energy was positive (e.g. not -999.) and not equal to pz
      return 0.5*TMath::Log((e+pz)/(e-pz));
    } else { // energy not known or equal to pz
      return -999.;
    }
  } else { // pid unknown
    return -999.;
  }
}

void AliAODTrack::SetTOFLabel(const Int_t *p) {  
  // Sets  (in TOF)
  for (Int_t i = 0; i < 3; i++) fTOFLabel[i]=p[i];
}

//_______________________________________________________________________
void AliAODTrack::GetTOFLabel(Int_t *p) const {
  // Gets (in TOF)
  for (Int_t i=0; i<3; i++) p[i]=fTOFLabel[i];
}

//______________________________________________________________________________
AliAODTrack::AODTrkPID_t AliAODTrack::GetMostProbablePID() const 
{
  // Returns the most probable PID array element.
  
  Int_t nPID = 10;
  AODTrkPID_t loc = kUnknown;
  Bool_t allTheSame = kTRUE;
  if (fPID) {
    Double_t max = 0.;
    for (Int_t iPID = 0; iPID < nPID; iPID++) {
      if (fPID[iPID] >= max) {
	if (fPID[iPID] > max) {
	  allTheSame = kFALSE;
	  max = fPID[iPID];
	  loc = (AODTrkPID_t)iPID;
	} else {
	  allTheSame = kTRUE;
	}
      }
    }
  }
  return allTheSame ? AODTrkPID_t(GetPIDForTracking()) : loc;
}

//______________________________________________________________________________
void AliAODTrack::ConvertAliPIDtoAODPID()
{
  // Converts AliPID array.
  // The numbering scheme is the same for electrons, muons, pions, kaons, and protons.
  // Everything else has to be set to zero.
  if (fPID) {
    fPID[kDeuteron] = 0.;
    fPID[kTriton]   = 0.;
    fPID[kHelium3]  = 0.;
    fPID[kAlpha]    = 0.;
    fPID[kUnknown]  = 0.;
  }
  return;
}

/*
//______________________________________________________________________________
template <typename T> void AliAODTrack::SetPosition(const T *x, const Bool_t dca) 
{
  // set the position

  if (x) {
    if (!dca) {
      ResetBit(kIsDCA);

      fPosition[0] = x[0];
      fPosition[1] = x[1];
      fPosition[2] = x[2];
    } else {
      SetBit(kIsDCA);
      // don't know any better yet
      fPosition[0] = -999.;
      fPosition[1] = -999.;
      fPosition[2] = -999.;
    }
  } else {
    ResetBit(kIsDCA);

    fPosition[0] = -999.;
    fPosition[1] = -999.;
    fPosition[2] = -999.;
  }
}
*/
//______________________________________________________________________________
void AliAODTrack::SetDCA(Double_t d, Double_t z) 
{
  // set the dca
  fPosition[0] = d;
  fPosition[1] = z;
  fPosition[2] = 0.;
  SetBit(kIsDCA);
}

//______________________________________________________________________________
void AliAODTrack::Print(Option_t* /* option */) const
{
  // prints information about AliAODTrack

  printf("Object name: %s   Track type: %s\n", GetName(), GetTitle()); 
  printf("        px = %f\n", Px());
  printf("        py = %f\n", Py());
  printf("        pz = %f\n", Pz());
  printf("        pt = %f\n", Pt());
  printf("      1/pt = %f\n", OneOverPt());
  printf("     theta = %f\n", Theta());
  printf("       phi = %f\n", Phi());
  printf("  chi2/NDF = %f\n", Chi2perNDF());
  printf("    charge = %d\n", Charge());
}

//______________________________________________________________________________
void AliAODTrack::SetMatchTrigger(Int_t matchTrig)
{
  // Set the MUON trigger information
  switch(matchTrig){
    case 0: // 0 track does not match trigger
      fITSMuonClusterMap=fITSMuonClusterMap&0x3fffffff;
      break;
    case 1: // 1 track match but does not pass pt cut
      fITSMuonClusterMap=(fITSMuonClusterMap&0x3fffffff)|0x40000000;
      break;
    case 2: // 2 track match Low pt cut
      fITSMuonClusterMap=(fITSMuonClusterMap&0x3fffffff)|0x80000000;
      break;
    case 3: // 3 track match High pt cut
      fITSMuonClusterMap=fITSMuonClusterMap|0xc0000000;
      break;
    default:
      fITSMuonClusterMap=fITSMuonClusterMap&0x3fffffff;
      AliWarning(Form("unknown case for matchTrig: %d\n",matchTrig));
  }
}

//______________________________________________________________________________
Bool_t AliAODTrack::HitsMuonChamber(Int_t MuonChamber, Int_t cathode) const
{
  // return kTRUE if the track fires the given tracking or trigger chamber.
  // If the chamber is a trigger one:
  // - if cathode = 0 or 1, the track matches the corresponding cathode
  // - if cathode = -1, the track matches both cathodes
  
  if (MuonChamber < 0) return kFALSE;
  
  if (MuonChamber < 10) return TESTBIT(GetMUONClusterMap(), MuonChamber);
  
  if (MuonChamber < 14) {
    
    if (cathode < 0) return TESTBIT(GetMUONTrigHitsMapTrg(), 13-MuonChamber) &&
                            TESTBIT(GetMUONTrigHitsMapTrg(), 13-MuonChamber+4);
    
    if (cathode < 2) return TESTBIT(GetMUONTrigHitsMapTrg(), 13-MuonChamber+(1-cathode)*4);
    
  }
  
  return kFALSE;
}

//______________________________________________________________________________
Bool_t AliAODTrack::MatchTriggerDigits() const
{
  // return kTRUE if the track matches a digit on both planes of at least 2 trigger chambers
  
  Int_t nMatchedChambers = 0;
  for (Int_t ich=10; ich<14; ich++) if (HitsMuonChamber(ich)) nMatchedChambers++;
  
  return (nMatchedChambers >= 2);
}

//______________________________________________________________________________
Int_t AliAODTrack::GetMuonTrigDevSign() const
{
  /// Return the sign of the  MTR deviation

  Int_t signInfo = (Int_t)((fMUONtrigHitsMapTrg>>30)&0x3);
  // Dummy value for old AODs which do not have the info
  if ( signInfo == 0 ) return -999;
  return signInfo - 2;
}

//______________________________________________________________________________
Bool_t AliAODTrack::PropagateToDCA(const AliVVertex *vtx, 
    Double_t b, Double_t maxd, Double_t dz[2], Double_t covar[3])
{
  // compute impact parameters to the vertex vtx and their covariance matrix
  // b is the Bz, needed to propagate correctly the track to vertex 
  // only the track parameters are update after the propagation (pos and mom),
  // not the covariance matrix. This is OK for propagation over short distance
  // inside the beam pipe.
  // return kFALSE is something went wrong

  // allowed only for tracks inside the beam pipe
  Float_t xstart2 = fPosition[0]*fPosition[0]+fPosition[1]*fPosition[1];
  if(xstart2 > 3.*3.) { // outside beampipe radius
    AliError("This method can be used only for propagation inside the beam pipe");
    return kFALSE; 
  }

  // convert to AliExternalTrackParam
  AliExternalTrackParam etp; etp.CopyFromVTrack(this);  

  // propagate
  if(!etp.PropagateToDCA(vtx,b,maxd,dz,covar)) return kFALSE;

  // update track position and momentum
  Double_t mom[3];
  etp.GetPxPyPz(mom);
  SetP(mom,kTRUE);
  etp.GetXYZ(mom);
  SetPosition(mom,kFALSE);


  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODTrack::GetPxPyPz(Double_t p[3]) const 
{
    //---------------------------------------------------------------------
    // This function returns the global track momentum components
    //---------------------------------------------------------------------
  p[0]=Px(); p[1]=Py(); p[2]=Pz();
  return kTRUE;
}


//_______________________________________________________________________
Float_t AliAODTrack::GetTPCClusterInfo(Int_t nNeighbours/*=3*/, Int_t type/*=0*/, Int_t row0, Int_t row1, Int_t bitType ) const
{
  //
  // TPC cluster information 
  // type 0: get fraction of found/findable clusters with neighbourhood definition
  //      1: findable clusters with neighbourhood definition
  //      2: found clusters
  // bitType:
  //      0 - all cluster used
  //      1 - clusters  used for the kalman update
  // definition of findable clusters:
  //            a cluster is defined as findable if there is another cluster
  //           within +- nNeighbours pad rows. The idea is to overcome threshold
  //           effects with a very simple algorithm.
  //

  
  Int_t found=0;
  Int_t findable=0;
  Int_t last=-nNeighbours;
  const TBits & clusterMap = (bitType%2==0) ? fTPCClusterMap : fTPCFitMap;
  
  Int_t upperBound=clusterMap.GetNbits();
  if (upperBound>row1) upperBound=row1;
  for (Int_t i=row0; i<upperBound; ++i){
    //look to current row
    if (clusterMap[i]) {
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
      if (clusterMap[j]){
        ++findable;
        break;
      }
    }
  }
  if (type==2) return found;
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


//______________________________________________________________________________
Double_t  AliAODTrack::GetTRDslice(Int_t plane, Int_t slice) const {
  //
  // return TRD Pid information
  //
  if (!fDetPid) return -1;
  Double32_t *trdSlices=fDetPid->GetTRDslices();
  if (!trdSlices) return -1;
  if ((plane<0) || (plane>=kTRDnPlanes)) {
    return -1.;
  }

  Int_t ns=fDetPid->GetTRDnSlices();
  if ((slice<-1) || (slice>=ns)) {
    return -1.;
  }

  if(slice>=0) return trdSlices[plane*ns + slice];

  // return average of the dEdx measurements
  Double_t q=0.; Double32_t *s = &trdSlices[plane*ns];
  for (Int_t i=0; i<ns; i++, s++) if((*s)>0.) q+=(*s);
  return q/ns;
}

//______________________________________________________________________________
UChar_t AliAODTrack::GetTRDntrackletsPID() const{
  //
  // return number of tracklets calculated from the slices
  //
  if(!fDetPid) return -1;
  return fDetPid->GetTRDntrackletsPID();
}

//______________________________________________________________________________
UChar_t AliAODTrack::GetTRDncls(Int_t layer) const {
  // 
  // return number of TRD clusters
  //
  if(!fDetPid || layer > 5) return -1;
  if(layer < 0) return fDetPid->GetTRDncls();
  else return fDetPid->GetTRDncls(layer);
}

//______________________________________________________________________________
Double_t AliAODTrack::GetTRDmomentum(Int_t plane, Double_t */*sp*/) const
{
  //Returns momentum estimation
  // in TRD layer "plane".

  if (!fDetPid) return -1;
  const Double_t *trdMomentum=fDetPid->GetTRDmomentum();

  if (!trdMomentum) {
    return -1.;
  }
  if ((plane<0) || (plane>=kTRDnPlanes)) {
    return -1.;
  }

  return trdMomentum[plane];
}

//_______________________________________________________________________
Int_t AliAODTrack::GetTOFBunchCrossing(Double_t b, Bool_t) const 
{
  // Returns the number of bunch crossings after trigger (assuming 25ns spacing)
  const double kSpacing = 25e3; // min interbanch spacing
  const double kShift = 0;
  Int_t bcid = kTOFBCNA; // defualt one
  if (!IsOn(kTOFout) || !IsOn(kESDpid)) return bcid; // no info
  //
  double tdif = GetTOFsignal();
  if (IsOn(kTIME)) { // integrated time info is there
    int pid = (int)GetMostProbablePID();
    double ttimes[10]; 
    GetIntegratedTimes(ttimes, pid>=AliPID::kSPECIES ? AliPID::kSPECIESC : AliPID::kSPECIES);
    tdif -= ttimes[pid];
  }
  else { // assume integrated time info from TOF radius and momentum
    const double kRTOF = 385.;
    const double kCSpeed = 3.e-2; // cm/ps
    double p = P();
    if (p<0.001) p = 1.0;
    double m = M();
    double path =  kRTOF;     // mean TOF radius
    if (TMath::Abs(b)>kAlmost0) {  // account for curvature
      double curv = Pt()/(b*kB2C);
      if (curv>kAlmost0) {
	double tgl = Pz()/Pt();
	path = 2./curv*TMath::ASin(kRTOF*curv/2.)*TMath::Sqrt(1.+tgl*tgl);
      }
    }
    tdif -= path/kCSpeed*TMath::Sqrt(1.+m*m/(p*p));
  }
  bcid = TMath::Nint((tdif - kShift)/kSpacing);
  return bcid;
}

void AliAODTrack::SetDetectorPID(const AliDetectorPID *pid)
{
  //
  // Set the detector PID
  //
  if (fDetectorPID) delete fDetectorPID;
  fDetectorPID=pid;
  
}

//_____________________________________________________________________________
Double_t AliAODTrack::GetHMPIDsignal() const
{
  if(fAODEvent->GetHMPIDringForTrackID(fID)) return fAODEvent->GetHMPIDringForTrackID(fID)->GetHmpSignal();
  else return -999.;
}

//_____________________________________________________________________________
Double_t AliAODTrack::GetHMPIDoccupancy() const
{
  if(fAODEvent->GetHMPIDringForTrackID(fID)) return fAODEvent->GetHMPIDringForTrackID(fID)->GetHmpOccupancy();
  else return -999.;
}

//_____________________________________________________________________________
Int_t AliAODTrack::GetHMPIDcluIdx() const
{
  if(fAODEvent->GetHMPIDringForTrackID(fID)) return fAODEvent->GetHMPIDringForTrackID(fID)->GetHmpCluIdx();
  else return -999;
}

//_____________________________________________________________________________
void AliAODTrack::GetHMPIDtrk(Float_t &x, Float_t &y, Float_t &th, Float_t &ph) const
{
  x = -999; y = -999.; th = -999.; ph = -999.;

  const AliAODHMPIDrings *ring=fAODEvent->GetHMPIDringForTrackID(fID);
  if(ring){
    x  = ring->GetHmpTrackX();
    y  = ring->GetHmpTrackY();
    th = ring->GetHmpTrackTheta();
    ph = ring->GetHmpTrackPhi();
  }
}

//_____________________________________________________________________________
void AliAODTrack::GetHMPIDmip(Float_t &x,Float_t &y,Int_t &q, Int_t &nph) const
{
  x = -999; y = -999.; q = -999; nph = -999;
  
  const AliAODHMPIDrings *ring=fAODEvent->GetHMPIDringForTrackID(fID);
  if(ring){
    x   = ring->GetHmpMipX();
    y   = ring->GetHmpMipY();
    q   = (Int_t)ring->GetHmpMipCharge();
    nph = (Int_t)ring->GetHmpNumOfPhotonClusters();
  }
}

//_____________________________________________________________________________
Bool_t AliAODTrack::GetOuterHmpPxPyPz(Double_t *p) const 
{ 
 if(fAODEvent->GetHMPIDringForTrackID(fID)) {fAODEvent->GetHMPIDringForTrackID(fID)->GetHmpMom(p); return kTRUE;}
 
 else return kFALSE;      
}
//_____________________________________________________________________________
Bool_t AliAODTrack::GetXYZAt(Double_t x, Double_t b, Double_t *r) const
{
  //---------------------------------------------------------------------
  // This function returns the global track position extrapolated to
  // the radial position "x" (cm) in the magnetic field "b" (kG)
  //---------------------------------------------------------------------

  //conversion of track parameter representation is
  //based on the implementation of AliExternalTrackParam::Set(...)
  //maybe some of this code can be moved to AliVTrack to avoid code duplication
  Double_t alpha=0.0;
  Double_t radPos2 = fPosition[0]*fPosition[0]+fPosition[1]*fPosition[1];  
  Double_t radMax  = 45.; // approximately ITS outer radius
  if (radPos2 < radMax*radMax) { // inside the ITS     
    alpha = fMomentum[1]; //TMath::ATan2(fMomentum[1],fMomentum[0]); // fMom is pt,phi,theta!
  } else { // outside the ITS
     Float_t phiPos = TMath::Pi()+TMath::ATan2(-fPosition[1], -fPosition[0]);
     alpha = 
     TMath::DegToRad()*(20*((((Int_t)(phiPos*TMath::RadToDeg()))/20))+10);
  }
  //
  // Get the vertex of origin and the momentum
  TVector3 ver(fPosition[0],fPosition[1],fPosition[2]);
  TVector3 mom(Px(),Py(),Pz());
  //
  // Rotate to the local coordinate system
  ver.RotateZ(-alpha);
  mom.RotateZ(-alpha);

  Double_t param0 = ver.Y();
  Double_t param1 = ver.Z();
  Double_t param2 = TMath::Sin(mom.Phi());
  Double_t param3 = mom.Pz()/mom.Pt();
  Double_t param4 = TMath::Sign(1/mom.Pt(),(Double_t)fCharge);

  //calculate the propagated coordinates
  //this is based on AliExternalTrackParam::GetXYZAt(Double_t x, Double_t b, Double_t *r)
  Double_t dx=x-ver.X();
  if(TMath::Abs(dx)<=kAlmost0) return GetXYZ(r);

  Double_t f1=param2;
  Double_t f2=f1 + dx*param4*b*kB2C;

  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
  
  Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
  r[0] = x;
  r[1] = param0 + dx*(f1+f2)/(r1+r2);
  r[2] = param1 + dx*(r2 + f2*(f1+f2)/(r1+r2))*param3;//Thanks to Andrea & Peter
  return Local2GlobalPosition(r,alpha);
}

//_____________________________________________________________________________
Bool_t AliAODTrack::GetXYZatR(Double_t xr,Double_t bz, Double_t *xyz, Double_t* alpSect) const
{
  // This method has 3 modes of behaviour
  // 1) xyz[3] array is provided but alpSect pointer is 0: calculate the position of track intersection 
  //    with circle of radius xr and fill it in xyz array
  // 2) alpSect pointer is provided: find alpha of the sector where the track reaches local coordinate xr
  //    Note that in this case xr is NOT the radius but the local coordinate.
  //    If the xyz array is provided, it will be filled by track lab coordinates at local X in this sector
  // 3) Neither alpSect nor xyz pointers are provided: just check if the track reaches radius xr
  //
  //
  Double_t alpha=0.0;
  Double_t radPos2 = fPosition[0]*fPosition[0]+fPosition[1]*fPosition[1];  
  Double_t radMax  = 45.; // approximately ITS outer radius
  if (radPos2 < radMax*radMax) { // inside the ITS     
    alpha = fMomentum[1]; //TMath::ATan2(fMomentum[1],fMomentum[0]); // fMom is pt,phi,theta!
  } else { // outside the ITS
     Float_t phiPos = TMath::Pi()+TMath::ATan2(-fPosition[1], -fPosition[0]);
     alpha = 
     TMath::DegToRad()*(20*((((Int_t)(phiPos*TMath::RadToDeg()))/20))+10);
  }
  //  
  // Get the vertex of origin and the momentum
  TVector3 ver(fPosition[0],fPosition[1],fPosition[2]);
  TVector3 mom(Px(),Py(),Pz());
  //
  // Rotate to the local coordinate system
  ver.RotateZ(-alpha);
  mom.RotateZ(-alpha);
  //
  Double_t fx = ver.X();
  Double_t fy = ver.Y();
  Double_t fz = ver.Z();
  Double_t sn = TMath::Sin(mom.Phi());
  Double_t tgl = mom.Pz()/mom.Pt();
  Double_t crv = TMath::Sign(1/mom.Pt(),(Double_t)fCharge)*bz*kB2C;
  //
  if ( (TMath::Abs(bz))<kAlmost0Field ) crv=0.;
  //
  // general circle parameterization:
  // x = (r0+tR)cos(phi0) - tR cos(t+phi0)
  // y = (r0+tR)sin(phi0) - tR sin(t+phi0)
  // where qb is the sign of the curvature, tR is the track's signed radius and r0 
  // is the DCA of helix to origin
  //
  double tR = 1./crv;            // track radius signed
  double cs = TMath::Sqrt((1-sn)*(1+sn));
  double x0 = fx - sn*tR;        // helix center coordinates
  double y0 = fy + cs*tR;
  double phi0 = TMath::ATan2(y0,x0);  // angle of PCA wrt to the origin
  if (tR<0) phi0 += TMath::Pi();
  if      (phi0 > TMath::Pi()) phi0 -= 2.*TMath::Pi();
  else if (phi0 <-TMath::Pi()) phi0 += 2.*TMath::Pi();
  double cs0 = TMath::Cos(phi0);
  double sn0 = TMath::Sin(phi0);
  double r0 = x0*cs0 + y0*sn0 - tR; // DCA to origin
  double r2R = 1.+r0/tR;
  //
  //
  if (r2R<kAlmost0) return kFALSE;  // helix is centered at the origin, no specific intersection with other concetric circle
  if (!xyz && !alpSect) return kTRUE;
  double xr2R = xr/tR;
  double r2Ri = 1./r2R;
  // the intersection cos(t) = [1 + (r0/tR+1)^2 - (r0/tR)^2]/[2(1+r0/tR)]
  double cosT = 0.5*(r2R + (1-xr2R*xr2R)*r2Ri);
  if ( TMath::Abs(cosT)>kAlmost1 ) {
    //    printf("Does not reach : %f %f\n",r0,tR);
    return kFALSE; // track does not reach the radius xr
  }
  //
  double t = TMath::ACos(cosT);
  if (tR<0) t = -t;
  // intersection point
  double xyzi[3];
  xyzi[0] = x0 - tR*TMath::Cos(t+phi0);
  xyzi[1] = y0 - tR*TMath::Sin(t+phi0);
  if (xyz) { // if postition is requested, then z is needed:
    double t0 = TMath::ATan2(cs,-sn) - phi0;
    double z0 = fz - t0*tR*tgl;    
    xyzi[2] = z0 + tR*t*tgl;
  }
  else xyzi[2] = 0;
  //
  Local2GlobalPosition(xyzi,alpha);
  //
  if (xyz) {
    xyz[0] = xyzi[0];
    xyz[1] = xyzi[1];
    xyz[2] = xyzi[2];
  }
  //
  if (alpSect) {
    double &alp = *alpSect;
    // determine the sector of crossing
    double phiPos = TMath::Pi()+TMath::ATan2(-xyzi[1],-xyzi[0]);
    int sect = ((Int_t)(phiPos*TMath::RadToDeg()))/20;
    alp = TMath::DegToRad()*(20*sect+10);
    double x2r,f1,f2,r1,r2,dx,dy2dx,yloc=0, ylocMax = xr*TMath::Tan(TMath::Pi()/18); // min max Y within sector at given X
    //
    while(1) {
      Double_t ca=TMath::Cos(alp-alpha), sa=TMath::Sin(alp-alpha);
      if ((cs*ca+sn*sa)<0) {
	AliDebug(1,Form("Rotation to target sector impossible: local cos(phi) would become %.2f",cs*ca+sn*sa));
	return kFALSE;
      }
      //
      f1 = sn*ca - cs*sa;
      if (TMath::Abs(f1) >= kAlmost1) {
	AliDebug(1,Form("Rotation to target sector impossible: local sin(phi) would become %.2f",f1));
	return kFALSE;
      }
      //
      double tmpX =  fx*ca + fy*sa;
      double tmpY = -fx*sa + fy*ca;
      //
      // estimate Y at X=xr
      dx=xr-tmpX;
      x2r = crv*dx;
      f2=f1 + x2r;
      if (TMath::Abs(f2) >= kAlmost1) {
	AliDebug(1,Form("Propagation in target sector failed ! %.10e",f2));
	return kFALSE;
      }
      r1 = TMath::Sqrt((1.-f1)*(1.+f1));
      r2 = TMath::Sqrt((1.-f2)*(1.+f2));
      dy2dx = (f1+f2)/(r1+r2);
      yloc = tmpY + dx*dy2dx;
      if      (yloc>ylocMax)  {alp += 2*TMath::Pi()/18; sect++;}
      else if (yloc<-ylocMax) {alp -= 2*TMath::Pi()/18; sect--;}
      else break;
      if      (alp >= TMath::Pi()) alp -= 2*TMath::Pi();
      else if (alp < -TMath::Pi()) alp += 2*TMath::Pi();
      //      if (sect>=18) sect = 0;
      //      if (sect<=0) sect = 17;
    }
    //
    // if alpha was requested, then recalculate the position at intersection in sector
    if (xyz) {
      xyz[0] = xr;
      xyz[1] = yloc;
      if (TMath::Abs(x2r)<0.05) xyz[2] = fz + dx*(r2 + f2*dy2dx)*tgl;
      else {
	// for small dx/R the linear apporximation of the arc by the segment is OK,
	// but at large dx/R the error is very large and leads to incorrect Z propagation
	// angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
	// The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
	// Similarly, the rotation angle in linear in dx only for dx<<R
	double chord = dx*TMath::Sqrt(1+dy2dx*dy2dx);   // distance from old position to new one
	double rot = 2*TMath::ASin(0.5*chord*crv); // angular difference seen from the circle center
	xyz[2] = fz + rot/crv*tgl;
      }
      Local2GlobalPosition(xyz,alp);
    }
  }
  return kTRUE;    
  //
}

//_______________________________________________________
void  AliAODTrack::GetITSdEdxSamples(Double_t s[4]) const
{
  // get ITS dedx samples
  if (!fDetPid) for (int i=4;i--;) s[i]=0;
  else          for (int i=4;i--;) s[i] = fDetPid->GetITSdEdxSample(i);
}

//_____________________________________________
Double_t AliAODTrack::GetMassForTracking() const
{
  int pid = fPIDForTracking;
  if (pid<AliPID::kPion) pid = AliPID::kPion;
  double m = AliPID::ParticleMass(fPIDForTracking);
  return (fPIDForTracking==AliPID::kHe3 || fPIDForTracking==AliPID::kAlpha) ? -m : m;
}
//_______________________________________________________
const AliTOFHeader* AliAODTrack::GetTOFHeader() const {
  return fAODEvent->GetTOFHeader();
}
  
//_______________________________________________________
Int_t AliAODTrack::GetNcls(Int_t idet) const
{
  // Get number of clusters by subdetector index
  //
  Int_t ncls = 0;
  switch(idet){
  case 0:
    ncls = GetITSNcls();
    break;
  case 1:
    ncls = (Int_t)GetTPCNcls();
    break;
  case 2:
    ncls = (Int_t)GetTRDncls();
    break;
  case 3:
    break;
    /*if (fTOFindex != -1)
      ncls = 1;*/
    break;
  case 4: //PHOS
    break;
  case 5: //HMPID
    break;
    if ((GetHMPIDcluIdx() >= 0) && (GetHMPIDcluIdx() < 7000000)) {
      if ((GetHMPIDcluIdx()%1000000 != 9999) && (GetHMPIDcluIdx()%1000000 != 99999)) {
	ncls = 1;
	}
    }    
    break;
  default:
    break;
  }
  return ncls;
}

Int_t AliAODTrack::GetTrackParam         ( AliExternalTrackParam & ) const {return 0;} 
Int_t AliAODTrack::GetTrackParamRefitted ( AliExternalTrackParam & ) const {return 0;} 
Int_t AliAODTrack::GetTrackParamIp       ( AliExternalTrackParam & ) const {return 0;} 
Int_t AliAODTrack::GetTrackParamTPCInner ( AliExternalTrackParam & ) const {return 0;} 
Int_t AliAODTrack::GetTrackParamOp       ( AliExternalTrackParam & ) const {return 0;} 
Int_t AliAODTrack::GetTrackParamCp       ( AliExternalTrackParam & ) const {return 0;} 
Int_t AliAODTrack::GetTrackParamITSOut   ( AliExternalTrackParam & ) const {return 0;} 

