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
  fEventNumber(0),
  fRunNumber(0),
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
  fPHOSClusters(0), 
  fFirstPHOSCluster(-1),
  fESDFMD(0x0)
{
  for (Int_t i=0; i<24; i++) {
    fT0time[i] = 0;
    fT0amplitude[i] = 0;
  }
}
//______________________________________________________________________________
AliESD::AliESD(const AliESD& esd):
  TObject(esd),
  fEventNumber(esd.fEventNumber),
  fRunNumber(esd.fRunNumber),
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
  fPHOSClusters(esd.fPHOSClusters), 
  fFirstPHOSCluster(esd.fFirstPHOSCluster),
  fESDFMD(esd.fESDFMD)
{
  for (Int_t i=0; i<24; i++) {
    fT0time[i] = esd.fT0time[i];
    fT0amplitude[i] = esd.fT0amplitude[i];
  }
}

//______________________________________________________________________________
AliESD & AliESD::operator=(const AliESD& source) {

  // Assignment operator

  if(&source == this) return *this;

  fEventNumber = source.fEventNumber;
  fRunNumber = source.fRunNumber;
  fTriggerMask = source.fTriggerMask;
  fTriggerCluster = source.fTriggerCluster;
  fRecoVersion = source.fRecoVersion;
  fMagneticField = source.fMagneticField;
  fZDCN1Energy = source.fZDCN1Energy;
  fZDCP1Energy = source.fZDCP1Energy;
  fZDCN2Energy = source.fZDCN2Energy;
  fZDCP2Energy = source.fZDCP2Energy;
  fZDCEMEnergy = source.fZDCEMEnergy;
  fZDCParticipants = source.fZDCParticipants;
  fT0zVertex = source.fT0zVertex;
  fSPDVertex = source.fSPDVertex;
  fPrimaryVertex = source.fPrimaryVertex;
  fSPDMult = source.fSPDMult;
  fT0timeStart = source.fT0timeStart;
  fTracks = *((TClonesArray*)source.fTracks.Clone());
  fHLTConfMapTracks = *((TClonesArray*)source.fHLTConfMapTracks.Clone());
  fHLTHoughTracks = *((TClonesArray*)source.fHLTHoughTracks.Clone());
  fMuonTracks = *((TClonesArray*)source.fMuonTracks.Clone());
  fPmdTracks = *((TClonesArray*)source.fPmdTracks.Clone());
  fTrdTracks = *((TClonesArray*)source.fTrdTracks.Clone());
  fV0s = *((TClonesArray*)source.fV0s.Clone());
  fCascades = *((TClonesArray*)source.fCascades.Clone());
  fKinks = *((TClonesArray*)source.fKinks.Clone());
  fCaloClusters = *((TClonesArray*)source.fCaloClusters.Clone());
  fEMCALClusters = source.fEMCALClusters;
  fFirstEMCALCluster = source.fFirstEMCALCluster;
  fPHOSClusters = source.fPHOSClusters;
  fFirstPHOSCluster = source.fFirstPHOSCluster;
  fESDFMD = source.fESDFMD;

  for (Int_t i=0; i<24; i++) {
    fT0time[i] = source.fT0time[i];
    fT0amplitude[i] = source.fT0amplitude[i];
  }

  return *this;

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
}

//______________________________________________________________________________
void AliESD::Reset()
{
  fEventNumber=0;
  fRunNumber=0;
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
}

Int_t AliESD::AddV0(const AliESDv0 *v) {
  //
  // Add V0
  //
    Int_t idx=fV0s.GetEntriesFast();
    AliESDv0 *v0=new(fV0s[idx]) AliESDv0(*v);
    v0->SetID(idx);
    return idx;
}  

//______________________________________________________________________________
void AliESD::Print(Option_t *) const 
{
  //
  // Print header information of the event
  //
  printf("ESD run information\n");
  printf("Event # %d Run # %d Trigger %lld Magnetic field %f \n",
	 GetEventNumber(),
	 GetRunNumber(),
	 GetTriggerMask(),
	 GetMagneticField() );
    printf("Vertex: (%.4f +- %.4f, %.4f +- %.4f, %.4f +- %.4f) cm\n",
	   fPrimaryVertex.GetXv(), fPrimaryVertex.GetXRes(),
	   fPrimaryVertex.GetYv(), fPrimaryVertex.GetYRes(),
	   fPrimaryVertex.GetZv(), fPrimaryVertex.GetZRes());
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
    const AliESDtrack *t=GetTrack(i);
    const AliESDfriendTrack *f=t->GetFriendTrack();
    ev->AddTrack(f);
  }
}
