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

/* $Id$ */

//-----------------------------------------------------------------
//           Implementation of the ESD class
//   This is the class to deal with during the phisical analysis of data
//   This class is generated directly by the reconstruction methods
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "AliESD.h"
#include "AliESDfriend.h"

ClassImp(AliESD)

//______________________________________________________________________________
AliESD::AliESD():
  fEventNumberInFile(0),
  fBunchCrossNumber(0),
  fOrbitNumber(0),
  fPeriodNumber(0),
  fRunNumber(0),
  fTimeStamp(0),
  fEventType(0),
  fTriggerMask(0),
  fTriggerCluster(0),
  fRecoVersion(0),
  fMagneticField(0),
  fZDCN1Energy(0),
  fZDCP1Energy(0),
  fZDCN2Energy(0),
  fZDCP2Energy(0),
  fZDCEMEnergy(0),
  fZDCParticipants(0),
  fT0zVertex(0),
  fSPDVertex(),
  fPrimaryVertex(),
  fSPDMult(),
  fT0timeStart(0),
  fTracks("AliESDtrack",15000),
  fHLTConfMapTracks("AliESDHLTtrack",25000),
  fHLTHoughTracks("AliESDHLTtrack",15000),
  fMuonTracks("AliESDMuonTrack",30),
  fPmdTracks("AliESDPmdTrack",3000),
  fTrdTracks("AliESDTrdTrack",300),
  fV0s("AliESDv0",200),  
  fCascades("AliESDcascade",20),
  fKinks("AliESDkink",4000),
  fCaloClusters("AliESDCaloCluster",10000),
  fEMCALClusters(0), 
  fFirstEMCALCluster(-1),
  fEMCALTriggerPosition(0x0),
  fEMCALTriggerAmplitudes(0x0),
  fPHOSClusters(0), 
  fFirstPHOSCluster(-1),
  fPHOSTriggerPosition(0x0),
  fPHOSTriggerAmplitudes(0x0),
  fESDFMD(0x0),
  fESDVZERO(0x0),
  fErrorLogs("AliRawDataErrorLog",5)

{
  for (Int_t i=0; i<24; i++) {
    fT0time[i] = 0;
    fT0amplitude[i] = 0;
  }
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=0.;
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=0.;
}

//______________________________________________________________________________
AliESD::AliESD(const AliESD& esd):
  TObject(esd),
  fEventNumberInFile(esd.fEventNumberInFile),
  fBunchCrossNumber(esd.fBunchCrossNumber),
  fOrbitNumber(esd.fOrbitNumber),
  fPeriodNumber(esd.fPeriodNumber),
  fRunNumber(esd.fRunNumber),
  fTimeStamp(esd.fTimeStamp),
  fEventType(esd.fEventType),
  fTriggerMask(esd.fTriggerMask),
  fTriggerCluster(esd.fTriggerCluster),
  fRecoVersion(esd.fRecoVersion),
  fMagneticField(esd.fMagneticField),
  fZDCN1Energy(esd.fZDCN1Energy),
  fZDCP1Energy(esd.fZDCP1Energy),
  fZDCN2Energy(esd.fZDCN2Energy),
  fZDCP2Energy(esd.fZDCP2Energy),
  fZDCEMEnergy(esd.fZDCEMEnergy),
  fZDCParticipants(esd.fZDCParticipants),
  fT0zVertex(esd.fT0zVertex),
  fSPDVertex(esd.fSPDVertex),
  fPrimaryVertex(esd.fPrimaryVertex),
  fSPDMult(esd.fSPDMult),
  fT0timeStart(esd.fT0timeStart),
  fTracks(*((TClonesArray*)esd.fTracks.Clone())),
  fHLTConfMapTracks(*((TClonesArray*)esd.fHLTConfMapTracks.Clone())),
  fHLTHoughTracks(*((TClonesArray*)esd.fHLTHoughTracks.Clone())),
  fMuonTracks(*((TClonesArray*)esd.fMuonTracks.Clone())),
  fPmdTracks(*((TClonesArray*)esd.fPmdTracks.Clone())),
  fTrdTracks(*((TClonesArray*)esd.fTrdTracks.Clone())),
  fV0s(*((TClonesArray*)esd.fV0s.Clone())),  
  fCascades(*((TClonesArray*)esd.fCascades.Clone())),
  fKinks(*((TClonesArray*)esd.fKinks.Clone())),
  fCaloClusters(*((TClonesArray*)esd.fCaloClusters.Clone())),
  fEMCALClusters(esd.fEMCALClusters), 
  fFirstEMCALCluster(esd.fFirstEMCALCluster),
  fEMCALTriggerPosition(esd. fEMCALTriggerPosition),
  fEMCALTriggerAmplitudes(esd.fEMCALTriggerAmplitudes),
  fPHOSClusters(esd.fPHOSClusters), 
  fFirstPHOSCluster(esd.fFirstPHOSCluster),
  fPHOSTriggerPosition(esd.fPHOSTriggerPosition),
  fPHOSTriggerAmplitudes(esd.fPHOSTriggerAmplitudes),
  fESDFMD(esd.fESDFMD),
  fESDVZERO(esd.fESDVZERO),
  fErrorLogs(*((TClonesArray*)esd.fErrorLogs.Clone()))
{
  for (Int_t i=0; i<24; i++) {
    fT0time[i] = esd.fT0time[i];
    fT0amplitude[i] = esd.fT0amplitude[i];
  }
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=esd.fDiamondXY[i];
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=esd.fDiamondCovXY[i];
}

//______________________________________________________________________________
AliESD::~AliESD()
{
  //
  // Standard destructor
  //
  fTracks.Delete();
  fHLTConfMapTracks.Delete();
  fHLTHoughTracks.Delete();
  fMuonTracks.Delete();
  fPmdTracks.Delete();
  fTrdTracks.Delete();
  fV0s.Delete();
  fCascades.Delete();
  fKinks.Delete();
  fCaloClusters.Delete();
  delete fESDFMD;
  delete fESDVZERO;
//   fEMCALTriggerPosition->Delete();
//   fEMCALTriggerAmplitudes->Delete();
//   fPHOSTriggerPosition->Delete();
//   fPHOSTriggerAmplitudes->Delete();
//   delete fEMCALTriggerPosition;
//   delete fEMCALTriggerAmplitudes;
//   delete fPHOSTriggerPosition;
//   delete fPHOSTriggerAmplitudes;
  fErrorLogs.Delete();

}

//______________________________________________________________________________
void AliESD::Reset()
{
  fEventNumberInFile=0;
  fBunchCrossNumber=0;
  fOrbitNumber=0;
  fPeriodNumber=0;
  fRunNumber=0;
  fTimeStamp = 0;
  fEventType = 0;
  fTriggerMask=0;
  fTriggerCluster=0;
  fRecoVersion=0;
  fMagneticField=0;
  fZDCN1Energy=0;
  fZDCP1Energy=0;
  fZDCN2Energy=0;
  fZDCP2Energy=0;
  fZDCEMEnergy=0;
  fZDCParticipants=0;
  fT0zVertex=0;
  fT0timeStart = 0;
  new (&fSPDVertex) AliESDVertex();
  new (&fPrimaryVertex) AliESDVertex();
  new (&fSPDMult) AliMultiplicity();
  fTracks.Clear();
  fHLTConfMapTracks.Clear();
  fHLTHoughTracks.Clear();
  fMuonTracks.Clear();
  fPmdTracks.Clear();
  fTrdTracks.Clear();
  fV0s.Clear();
  fCascades.Clear();
  fCaloClusters.Clear();
  fEMCALClusters=0; 
  fFirstEMCALCluster=-1; 
  fPHOSClusters=0; 
  fFirstPHOSCluster=-1; 
  if (fESDFMD) fESDFMD->Clear();
//   fEMCALTriggerPosition->Clear();
//   fEMCALTriggerAmplitudes->Clear();
//   fPHOSTriggerPosition->Clear();
//   fPHOSTriggerAmplitudes->Clear();
  fErrorLogs.Clear();
}

Bool_t  AliESD::RemoveTrack(Int_t /*i*/) {
  // ---------------------------------------------------------
  // Remove track
  // ---------------------------------------------------------

  // Check if this track comes from a reconstructed decay
  // if (yes) return kFALSE

  // Remap the indices of the daughters of recosntructed decays

  // Remove the track
  // delete fTracks->RemoveAt(i);

  // Compress the array with tracks
  // fTracks->Compress();

  return kTRUE;
}

Int_t AliESD::AddV0(const AliESDv0 *v) {
  //
  // Add V0
  //
    Int_t idx=fV0s.GetEntriesFast();
    new(fV0s[idx]) AliESDv0(*v);
    return idx;
}  

//______________________________________________________________________________
void AliESD::Print(Option_t *) const 
{
  //
  // Print header information of the event
  //
  printf("ESD run information\n");
  printf("Event # in file %d Bunch crossing # %d Orbit # %d Period # %d Run # %d Trigger %lld Magnetic field %f \n",
	 GetEventNumberInFile(),
	 GetBunchCrossNumber(),
	 GetOrbitNumber(),
	 GetPeriodNumber(),
	 GetRunNumber(),
	 GetTriggerMask(),
	 GetMagneticField() );
    printf("Vertex: (%.4f +- %.4f, %.4f +- %.4f, %.4f +- %.4f) cm\n",
	   fPrimaryVertex.GetXv(), fPrimaryVertex.GetXRes(),
	   fPrimaryVertex.GetYv(), fPrimaryVertex.GetYRes(),
	   fPrimaryVertex.GetZv(), fPrimaryVertex.GetZRes());
    printf("Mean vertex in RUN: X=%.4f Y=%.4f cm\n",
	   GetDiamondX(),GetDiamondY());
    printf("SPD Multiplicity. Number of tracklets %d \n",
           fSPDMult.GetNumberOfTracklets());
  printf("Event from reconstruction version %d \n",fRecoVersion);
  printf("Number of tracks: \n");
  printf("                 charged   %d\n", GetNumberOfTracks());
  printf("                 hlt CF    %d\n", GetNumberOfHLTConfMapTracks());
  printf("                 hlt HT    %d\n", GetNumberOfHLTHoughTracks());
  printf("                 muon      %d\n", GetNumberOfMuonTracks());
  printf("                 pmd       %d\n", GetNumberOfPmdTracks());
  printf("                 trd       %d\n", GetNumberOfTrdTracks());
  printf("                 v0        %d\n", GetNumberOfV0s());
  printf("                 cascades  %d\n", GetNumberOfCascades());
  printf("                 kinks     %d\n", GetNumberOfKinks());
  printf("                 CaloClusters %d\n", GetNumberOfCaloClusters());
  printf("                 phos      %d\n", GetNumberOfPHOSClusters());
  printf("                 emcal     %d\n", GetNumberOfEMCALClusters());
  printf("                 FMD       %s\n", (fESDFMD ? "yes" : "no"));
  printf("                 VZERO     %s\n", (fESDVZERO ? "yes" : "no"));
}

void AliESD::SetESDfriend(const AliESDfriend *ev) {
  //
  // Attaches the complementary info to the ESD
  //
  if (!ev) return;

  Int_t ntrk=ev->GetNumberOfTracks();

  for (Int_t i=0; i<ntrk; i++) {
    const AliESDfriendTrack *f=ev->GetTrack(i);
    GetTrack(i)->SetFriendTrack(f);
  }
}

void AliESD::GetESDfriend(AliESDfriend *ev) const {
  //
  // Extracts the complementary info from the ESD
  //
  if (!ev) return;

  Int_t ntrk=GetNumberOfTracks();

  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=GetTrack(i);
    const AliESDfriendTrack *f=t->GetFriendTrack();
    ev->AddTrack(f);

    t->ReleaseESDfriendTrack();// Not to have two copies of "friendTrack"

  }
}
