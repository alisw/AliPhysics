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

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODCluster.h"


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
AliESD & AliESD::operator=(const AliESD& source) {

  // Assignment operator

  if(&source == this) return *this;

  fEventNumberInFile = source.fEventNumberInFile;
  fBunchCrossNumber = source.fBunchCrossNumber;
  fOrbitNumber = source.fOrbitNumber;
  fPeriodNumber = source.fPeriodNumber;
  fRunNumber = source.fRunNumber;
  fTimeStamp   = source.fTimeStamp;
  fEventType   = source.fEventType;
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
  fESDVZERO = source.fESDVZERO;
  fEMCALTriggerPosition=source. fEMCALTriggerPosition;
  fEMCALTriggerAmplitudes=source.fEMCALTriggerAmplitudes;
  fPHOSTriggerPosition=source.fPHOSTriggerPosition;
  fPHOSTriggerAmplitudes=source.fPHOSTriggerAmplitudes;
  fErrorLogs = *((TClonesArray*)source.fErrorLogs.Clone());

  for (Int_t i=0; i<24; i++) {
    fT0time[i] = source.fT0time[i];
    fT0amplitude[i] = source.fT0amplitude[i];
  }
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=source.fDiamondXY[i];
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=source.fDiamondCovXY[i];

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
    const AliESDtrack *t=GetTrack(i);
    const AliESDfriendTrack *f=t->GetFriendTrack();
    ev->AddTrack(f);
  }
}

AliAODEvent *AliESD::CreateAOD() const {
  //
  // Creates and returns the standard AOD from the ESD.
  // Make sure to delete it outside of this function!
  //
  
  // create an AliAOD object 
  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();

  // set arrays and pointers
  Float_t posF[3];
  Double_t pos[3];
  Double_t p[3];
  Double_t covVtx[6];
  Double_t covTr[21];
  Double_t pid[10];
  
  // Multiplicity information needed by the header (to be revised!)
  Int_t nTracks   = GetNumberOfTracks();
  Int_t nPosTracks = 0;
  for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) 
    if (GetTrack(iTrack)->GetSign()> 0) nPosTracks++;
  
  // create the header
  aod->AddHeader(new AliAODHeader(GetRunNumber(),
				  GetBunchCrossNumber(),
				  GetOrbitNumber(),
				  GetPeriodNumber(),
				  nTracks,
				  nPosTracks,
				  nTracks-nPosTracks,
				  GetMagneticField(),
				  -999., // fill muon magnetic field
				  -999., // centrality; to be filled, still
				  GetZDCN1Energy(),
				  GetZDCP1Energy(),
				  GetZDCN2Energy(),
				  GetZDCP2Energy(),
				  GetZDCEMEnergy(),
				  GetTriggerMask(),
				  GetTriggerCluster(),
				  GetEventType()));
  
  Int_t nV0s      = GetNumberOfV0s();
  Int_t nCascades = GetNumberOfCascades();
  Int_t nKinks    = GetNumberOfKinks();
  Int_t nVertices = nV0s + nCascades + nKinks;
  
  aod->ResetStd(nTracks, nVertices);
  AliAODTrack *aodTrack;
  
  // Array to take into account the tracks already added to the AOD
  Bool_t * usedTrack = NULL;
  if (nTracks>0) {
    usedTrack = new Bool_t[nTracks];
    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) usedTrack[iTrack]=kFALSE;
  }
  // Array to take into account the V0s already added to the AOD
  Bool_t * usedV0 = NULL;
  if (nV0s>0) {
    usedV0 = new Bool_t[nV0s];
    for (Int_t iV0=0; iV0<nV0s; ++iV0) usedV0[iV0]=kFALSE;
  }
  // Array to take into account the kinks already added to the AOD
  Bool_t * usedKink = NULL;
  if (nKinks>0) {
    usedKink = new Bool_t[nKinks];
    for (Int_t iKink=0; iKink<nKinks; ++iKink) usedKink[iKink]=kFALSE;
  }
  
  // Access to the AOD container of vertices
  TClonesArray &vertices = *(aod->GetVertices());
  Int_t jVertices=0;
  
  // Access to the AOD container of tracks
  TClonesArray &tracks = *(aod->GetTracks());
  Int_t jTracks=0; 
  
  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (V0, cascades, kinks)
  const AliESDVertex *vtx = GetPrimaryVertex();
  
  vtx->GetXYZ(pos); // position
  vtx->GetCovMatrix(covVtx); //covariance matrix
  
  AliAODVertex * primary = new(vertices[jVertices++])
    AliAODVertex(pos, covVtx, vtx->GetChi2toNDF(), NULL, AliAODVertex::kPrimary);
  
  // Create vertices starting from the most complex objects
  
  // Cascades
  for (Int_t nCascade = 0; nCascade < nCascades; ++nCascade) {
    AliESDcascade *cascade = GetCascade(nCascade);
    
    cascade->GetXYZ(pos[0], pos[1], pos[2]);
    cascade->GetPosCovXi(covVtx);
    
    // Add the cascade vertex
    AliAODVertex * vcascade = new(vertices[jVertices++]) AliAODVertex(pos,
								      covVtx,
								      cascade->GetChi2Xi(), // = chi2/NDF since NDF = 2*2-3
								      primary,
								      AliAODVertex::kCascade);
    
    primary->AddDaughter(vcascade);
    
    // Add the V0 from the cascade. The ESD class have to be optimized...
    // Now we have to search for the corresponding Vo in the list of V0s
    // using the indeces of the positive and negative tracks
    
    Int_t posFromV0 = cascade->GetPindex();
    Int_t negFromV0 = cascade->GetNindex();
    
    AliESDv0 * v0 = 0x0;
    Int_t indV0 = -1;
    
    for (Int_t iV0=0; iV0<nV0s; ++iV0) {
      
      v0 = GetV0(iV0);
      Int_t posV0 = v0->GetPindex();
      Int_t negV0 = v0->GetNindex();
      
      if (posV0==posFromV0 && negV0==negFromV0) {
	indV0 = iV0;
	break;
      }
    }
    
    AliAODVertex * vV0FromCascade = 0x0;
    
    if (indV0>-1 && !usedV0[indV0] ) {
      
      // the V0 exists in the array of V0s and is not used
      
      usedV0[indV0] = kTRUE;
      
      v0->GetXYZ(pos[0], pos[1], pos[2]);
      v0->GetPosCov(covVtx);
      
      vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(pos,
							       covVtx,
							       v0->GetChi2V0(), // = chi2/NDF since NDF = 2*2-3
							       vcascade,
							       AliAODVertex::kV0);
    } else {
      
      // the V0 doesn't exist in the array of V0s or was used
      printf("Error: cascade %d. The V0 %d doesn't exist in the array of V0s or was used!\n", nCascade, indV0);
      cascade->GetXYZ(pos[0], pos[1], pos[2]);
      cascade->GetPosCov(covVtx);
      
      vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(pos,
							       covVtx,
							       v0->GetChi2V0(), // = chi2/NDF since NDF = 2*2-3
							       vcascade,
							       AliAODVertex::kV0);
      vcascade->AddDaughter(vV0FromCascade);
    }
    
    // Add the positive tracks from the V0
    
    if (! usedTrack[posFromV0]) {
      
      usedTrack[posFromV0] = kTRUE;
      
      AliESDtrack *esdTrack = GetTrack(posFromV0);
      esdTrack->GetPxPyPz(p);
      esdTrack->GetXYZ(pos);
      esdTrack->GetCovarianceXYZPxPyPz(covTr);
      esdTrack->GetESDpid(pid);
      
      vV0FromCascade->AddDaughter(
	   aodTrack = new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
							 esdTrack->GetLabel(), 
							 p, 
							 kTRUE,
							 pos,
							 kFALSE,
							 covTr, 
							 (Short_t)esdTrack->GetSign(),
							 esdTrack->GetITSClusterMap(), 
							 pid,
							 vV0FromCascade,
							 kTRUE,  // check if this is right
							 kFALSE, // check if this is right
							 AliAODTrack::kSecondary)
	   );
      aodTrack->ConvertAliPIDtoAODPID();
    }
    else {
      printf("Error: cascade %d track %d has already been used!\n", nCascade, posFromV0);
    }
    
    // Add the negative tracks from the V0
    
    if (!usedTrack[negFromV0]) {
      
      usedTrack[negFromV0] = kTRUE;
      
      AliESDtrack *esdTrack = GetTrack(negFromV0);
      esdTrack->GetPxPyPz(p);
      esdTrack->GetXYZ(pos);
      esdTrack->GetCovarianceXYZPxPyPz(covTr);
      esdTrack->GetESDpid(pid);
      
      vV0FromCascade->AddDaughter(
	   aodTrack = new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
							 esdTrack->GetLabel(),
							 p,
							 kTRUE,
							 pos,
							 kFALSE,
							 covTr, 
							 (Short_t)esdTrack->GetSign(),
							 esdTrack->GetITSClusterMap(), 
							 pid,
							 vV0FromCascade,
							 kTRUE,  // check if this is right
							 kFALSE, // check if this is right
							 AliAODTrack::kSecondary)
	   );
      aodTrack->ConvertAliPIDtoAODPID();
    }
    else {
      printf("Error: cascade %d track %d has already been used!\n", nCascade,  negFromV0);
    }
    
    // Add the bachelor track from the cascade
    
    Int_t bachelor = cascade->GetBindex();
    
    if(!usedTrack[bachelor]) {
      
      usedTrack[bachelor] = kTRUE;
      
      AliESDtrack *esdTrack = GetTrack(bachelor);
      esdTrack->GetPxPyPz(p);
      esdTrack->GetXYZ(pos);
      esdTrack->GetCovarianceXYZPxPyPz(covTr);
      esdTrack->GetESDpid(pid);
      
      vcascade->AddDaughter(
       aodTrack = new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
						     esdTrack->GetLabel(),
						     p,
						     kTRUE,
						     pos,
						     kFALSE,
						     covTr, 
						     (Short_t)esdTrack->GetSign(),
						     esdTrack->GetITSClusterMap(), 
						     pid,
						     vcascade,
						     kTRUE,  // check if this is right
						     kFALSE, // check if this is right
						     AliAODTrack::kSecondary)
       );
      aodTrack->ConvertAliPIDtoAODPID();
    }
    else {
      printf("Error: cascade %d track %d has already been used!\n", nCascade, bachelor);
    }
    
    // Add the primary track of the cascade (if any)
    
  } // end of the loop on cascades
  
    // V0s
  
  for (Int_t nV0 = 0; nV0 < nV0s; ++nV0) {
    
    if (usedV0[nV0]) continue; // skip if aready added to the AOD
    
    AliESDv0 *v0 = GetV0(nV0);
    
    v0->GetXYZ(pos[0], pos[1], pos[2]);
    v0->GetPosCov(covVtx);
    
    AliAODVertex * vV0 = 
      new(vertices[jVertices++]) AliAODVertex(pos,
					      covVtx,
					      v0->GetChi2V0(), // = chi2/NDF since NDF = 2*2-3
					      primary,
					      AliAODVertex::kV0);
    primary->AddDaughter(vV0);
    
    Int_t posFromV0 = v0->GetPindex();
    Int_t negFromV0 = v0->GetNindex();
    
    // Add the positive tracks from the V0
    
    if (!usedTrack[posFromV0]) {
      
      usedTrack[posFromV0] = kTRUE;
      
      AliESDtrack *esdTrack = GetTrack(posFromV0);
      esdTrack->GetPxPyPz(p);
      esdTrack->GetXYZ(pos);
      esdTrack->GetCovarianceXYZPxPyPz(covTr);
      esdTrack->GetESDpid(pid);
      
      vV0->AddDaughter(
        aodTrack = new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
						      esdTrack->GetLabel(), 
						      p, 
						      kTRUE,
						      pos,
						      kFALSE,
						      covTr, 
						      (Short_t)esdTrack->GetSign(),
						      esdTrack->GetITSClusterMap(), 
						      pid,
						      vV0,
						      kTRUE,  // check if this is right
						      kFALSE, // check if this is right
						      AliAODTrack::kSecondary)
	);
      aodTrack->ConvertAliPIDtoAODPID();
    }
    else {
      printf("Error: V0 %d track %d has already been used!\n", nV0, posFromV0);
    }
    
    // Add the negative tracks from the V0
    
    if (!usedTrack[negFromV0]) {
      
      usedTrack[negFromV0] = kTRUE;
      
      AliESDtrack *esdTrack = GetTrack(negFromV0);
      esdTrack->GetPxPyPz(p);
      esdTrack->GetXYZ(pos);
      esdTrack->GetCovarianceXYZPxPyPz(covTr);
      esdTrack->GetESDpid(pid);
      
      vV0->AddDaughter(
	aodTrack = new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
						      esdTrack->GetLabel(),
						      p,
						      kTRUE,
						      pos,
						      kFALSE,
						      covTr, 
						      (Short_t)esdTrack->GetSign(),
						      esdTrack->GetITSClusterMap(), 
						      pid,
						      vV0,
						      kTRUE,  // check if this is right
						      kFALSE, // check if this is right
						      AliAODTrack::kSecondary)
	);
      aodTrack->ConvertAliPIDtoAODPID();
    }
    else {
      printf("Error: V0 %d track %d has already been used!\n", nV0, negFromV0);
    }
    
  } // end of the loop on V0s
  
  // Kinks: it is a big mess the access to the information in the kinks
  // The loop is on the tracks in order to find the mother and daugther of each kink
  
  for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) {
    
    
    AliESDtrack * esdTrack = GetTrack(iTrack);
    
    Int_t ikink = esdTrack->GetKinkIndex(0);
    
    if (ikink) {
      // Negative kink index: mother, positive: daughter
      
      // Search for the second track of the kink
      
      for (Int_t jTrack = iTrack+1; jTrack<nTracks; ++jTrack) {
	
	AliESDtrack * esdTrack1 = GetTrack(jTrack);
	
	Int_t jkink = esdTrack1->GetKinkIndex(0);
	
	if ( TMath::Abs(ikink)==TMath::Abs(jkink) ) {
	  
	  // The two tracks are from the same kink
	  
	  if (usedKink[TMath::Abs(ikink)-1]) continue; // skip used kinks
	  
	  Int_t imother = -1;
	  Int_t idaughter = -1;
	  
	  if (ikink<0 && jkink>0) {
	    
	    imother = iTrack;
	    idaughter = jTrack;
	  }
	  else if (ikink>0 && jkink<0) {
	    
	    imother = jTrack;
	    idaughter = iTrack;
	  }
	  else {
	    printf("Error: Wrong combination of kink indexes: %d %d\n", ikink, jkink);
	    continue;
	  }
	  
	  // Add the mother track
	  
	  AliAODTrack * mother = NULL;
	  
	  if (!usedTrack[imother]) {
	    
	    usedTrack[imother] = kTRUE;
	    
	    AliESDtrack *esdTrack = GetTrack(imother);
	    esdTrack->GetPxPyPz(p);
	    esdTrack->GetXYZ(pos);
	    esdTrack->GetCovarianceXYZPxPyPz(covTr);
	    esdTrack->GetESDpid(pid);
	    
	    mother = 
	      new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
						 esdTrack->GetLabel(),
						 p,
						 kTRUE,
						 pos,
						 kFALSE,
						 covTr, 
						 (Short_t)esdTrack->GetSign(),
						 esdTrack->GetITSClusterMap(), 
						 pid,
						 primary,
						 kTRUE, // check if this is right
						 kTRUE, // check if this is right
						 AliAODTrack::kPrimary);
	    primary->AddDaughter(mother);
	    mother->ConvertAliPIDtoAODPID();
	  }
	  else {
	    printf("Error: kink %d track %d has already been used!\n", TMath::Abs(ikink)-1, imother);
	  }
	  // Add the kink vertex
	  AliESDkink * kink = GetKink(TMath::Abs(ikink)-1);
	  
	  AliAODVertex * vkink = 
	    new(vertices[jVertices++]) AliAODVertex(kink->GetPosition(),
						    NULL,
						    0.,
						    mother,
						    AliAODVertex::kKink);
	  // Add the daughter track
	  
	  AliAODTrack * daughter = NULL;
	  
	  if (!usedTrack[idaughter]) {
	    
	    usedTrack[idaughter] = kTRUE;
	    
	    AliESDtrack *esdTrack = GetTrack(idaughter);
	    esdTrack->GetPxPyPz(p);
	    esdTrack->GetXYZ(pos);
	    esdTrack->GetCovarianceXYZPxPyPz(covTr);
	    esdTrack->GetESDpid(pid);
	    
	    daughter = 
	      new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
						 esdTrack->GetLabel(),
						 p,
						 kTRUE,
						 pos,
						 kFALSE,
						 covTr, 
						 (Short_t)esdTrack->GetSign(),
						 esdTrack->GetITSClusterMap(), 
						 pid,
						 vkink,
						 kTRUE, // check if this is right
						 kTRUE, // check if this is right
						 AliAODTrack::kPrimary);
	    vkink->AddDaughter(daughter);
	    daughter->ConvertAliPIDtoAODPID();
	  }
	  else {
	    printf("Error: kink %d track %d has already been used!\n", TMath::Abs(ikink)-1, idaughter);
	  }
	  
	}
  
      }
      
    }      
    
  }
  
  // Tracks (primary and orphan)
  
  for (Int_t nTrack = 0; nTrack < nTracks; ++nTrack) {
    
    
    if (usedTrack[nTrack]) continue;
    
    AliESDtrack *esdTrack = GetTrack(nTrack);
    esdTrack->GetPxPyPz(p);
    esdTrack->GetXYZ(pos);
    esdTrack->GetCovarianceXYZPxPyPz(covTr);
    esdTrack->GetESDpid(pid);
    
    Float_t impactXY, impactZ;
    
    esdTrack->GetImpactParameters(impactXY,impactZ);
    
    if (impactXY<3) {
      // track inside the beam pipe
      
      primary->AddDaughter(
        aodTrack = new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
						      esdTrack->GetLabel(),
						      p,
						      kTRUE,
						      pos,
						      kFALSE,
						      covTr, 
						      (Short_t)esdTrack->GetSign(),
						      esdTrack->GetITSClusterMap(), 
						      pid,
						      primary,
						      kTRUE, // check if this is right
						      kTRUE, // check if this is right
						      AliAODTrack::kPrimary)
	);
      aodTrack->ConvertAliPIDtoAODPID();
    }
    else {
      // outside the beam pipe: orphan track
      aodTrack =
	new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p,
					   kTRUE,
					   pos,
					   kFALSE,
					   covTr, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   NULL,
					   kFALSE, // check if this is right
					   kFALSE, // check if this is right
					   AliAODTrack::kOrphan);
      aodTrack->ConvertAliPIDtoAODPID();
    }	
  } // end of loop on tracks
  
  // muon tracks
  Int_t nMuTracks = GetNumberOfMuonTracks();
  for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {
    
    AliESDMuonTrack *esdMuTrack = GetMuonTrack(nMuTrack);     
    p[0] = esdMuTrack->Px(); 
    p[1] = esdMuTrack->Py(); 
    p[2] = esdMuTrack->Pz();
    pos[0] = primary->GetX(); 
    pos[1] = primary->GetY(); 
    pos[2] = primary->GetZ();
    
    // has to be changed once the muon pid is provided by the ESD
    for (Int_t i = 0; i < 10; pid[i++] = 0.); pid[AliAODTrack::kMuon]=1.;
    
    primary->AddDaughter(
	   new(tracks[jTracks++]) AliAODTrack(0, // no ID provided
					      0, // no label provided
					      p,
					      kTRUE,
					      pos,
					      kFALSE,
					      NULL, // no covariance matrix provided
					      (Short_t)-99, // no charge provided
					      0, // no ITSClusterMap
					      pid,
					      primary,
					      kTRUE,  // check if this is right
					      kTRUE,  // not used for vertex fit
					      AliAODTrack::kPrimary)
	   );
  }
  
  // Access to the AOD container of clusters
  TClonesArray &clusters = *(aod->GetClusters());
  Int_t jClusters=0;
  
  // Calo Clusters
  Int_t nClusters    = GetNumberOfCaloClusters();
  
  for (Int_t iClust=0; iClust<nClusters; ++iClust) {
    
    AliESDCaloCluster * cluster = GetCaloCluster(iClust);
    
    Int_t id = cluster->GetID();
    Int_t label = -1;
    Float_t energy = cluster->GetClusterEnergy();
    cluster->GetGlobalPosition(posF);
    AliAODVertex *prodVertex = primary;
    AliAODTrack *primTrack = NULL;
    Char_t ttype=AliAODCluster::kUndef;
    
    if (cluster->IsPHOS()) ttype=AliAODCluster::kPHOSNeutral;
    else if (cluster->IsEMCAL()) {
      
      if (cluster->GetClusterType() == AliESDCaloCluster::kPseudoCluster)
	ttype = AliAODCluster::kEMCALPseudoCluster;
      else
	ttype = AliAODCluster::kEMCALClusterv1;  
    
    }
    
    new(clusters[jClusters++]) AliAODCluster(id,
					     label,
					     energy,
					     pos,
					     NULL, // no covariance matrix provided
					     NULL, // no pid for clusters provided
					     prodVertex,
					     primTrack,
					     ttype);
    
  } // end of loop on calo clusters
  
  delete [] usedTrack;
  delete [] usedV0;
  delete [] usedKink;
  
  return aod;
}
