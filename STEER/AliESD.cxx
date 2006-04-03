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

ClassImp(AliESD)

//______________________________________________________________________________
AliESD::AliESD():
  fEventNumber(0),
  fRunNumber(0),
  fTrigger(0),
  fRecoVersion(0),
  fMagneticField(0),
  fZDCN1Energy(0),
  fZDCP1Energy(0),
  fZDCN2Energy(0),
  fZDCP2Energy(0),
  fZDCEMEnergy(0),
  fZDCParticipants(0),
  fT0zVertex(0),
  fPrimaryVertex(),
  fTracks("AliESDtrack",15000),
  fHLTConfMapTracks("AliESDHLTtrack",25000),
  fHLTHoughTracks("AliESDHLTtrack",15000),
  fMuonTracks("AliESDMuonTrack",30),
  fPmdTracks("AliESDPmdTrack",3000),
  fTrdTracks("AliESDTrdTrack",300),
  fV0s("AliESDv0",200),  
  fCascades("AliESDcascade",20),
  fKinks("AliESDkink",4000),
  fV0MIs("AliESDV0MI",4000),
  fCaloClusters("AliESDCaloCluster",10000),
  fEMCALClusters(0), 
  fFirstEMCALCluster(-1),
  fPHOSClusters(0), 
  fFirstPHOSCluster(-1),
  fESDFMD(0x0)
{

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
  fV0MIs.Delete();
  fCaloClusters.Delete();
  delete fESDFMD;
}

void AliESD::UpdateV0PIDs()
{
  //
  //
  //
  Int_t nV0 = GetNumberOfV0MIs();
  for (Int_t i=0;i<nV0;i++){
    AliESDV0MI * v0 = GetV0MI(i);
    AliESDtrack* tp = GetTrack(v0->GetIndex(0));
    AliESDtrack* tm = GetTrack(v0->GetIndex(1));
    if (!tm || !tp){
      printf("BBBUUUUUUUGGGG\n");
    }
    Double_t pp[5],pm[5];
    tp->GetESDpid(pp);
    tm->GetESDpid(pm);
    v0->UpdatePID(pp,pm);    
  }
}

//______________________________________________________________________________
void AliESD::Reset()
{
  fEventNumber=0;
  fRunNumber=0;
  fTrigger=0;
  fRecoVersion=0;
  fMagneticField=0;
  fZDCN1Energy=0;
  fZDCP1Energy=0;
  fZDCN2Energy=0;
  fZDCP2Energy=0;
  fZDCEMEnergy=0;
  fZDCParticipants=0;
  fT0zVertex=0;
  fPrimaryVertex.Reset();
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

//______________________________________________________________________________
void AliESD::Print(Option_t *) const 
{
  //
  // Print header information of the event
  //
  printf("ESD run information\n");
  printf("Event # %d Run # %d Trigger %ld Magnetic field %f \n",
	 GetEventNumber(),
	 GetRunNumber(),
	 GetTrigger(),
	 GetMagneticField() );
  printf("Vertex: (%.4f +- %.4f, %.4f +- %.4f, %.4f +- %.4f) cm\n",
	 fPrimaryVertex.GetXv(), fPrimaryVertex.GetXRes(),
	 fPrimaryVertex.GetYv(), fPrimaryVertex.GetYRes(),
	 fPrimaryVertex.GetZv(), fPrimaryVertex.GetZRes());
  printf("Event from reconstruction version %d \n",fRecoVersion);
  printf("Number of tracks: \n");
  printf("                 charged   %d\n", GetNumberOfTracks());
  printf("                 hlt CF    %d\n", GetNumberOfHLTConfMapTracks());
  printf("                 hlt HT    %d\n", GetNumberOfHLTHoughTracks());
  printf("                 phos      %d\n", GetNumberOfPHOSClusters());
  printf("                 emcal     %d\n", GetNumberOfEMCALClusters());
  printf("                 muon      %d\n", GetNumberOfMuonTracks());
  printf("                 pmd       %d\n", GetNumberOfPmdTracks());
  printf("                 trd       %d\n", GetNumberOfTrdTracks());
  printf("                 v0        %d\n", GetNumberOfV0s());
  printf("                 cascades  %d\n)", GetNumberOfCascades());
  printf("                 kinks     %d\n)", GetNumberOfKinks());
  printf("                 V0MIs     %d\n)", GetNumberOfV0MIs());
  printf("                 CaloClusters %d\n)", GetNumberOfCaloClusters());
  printf("                 FMD       %s\n)", (fESDFMD ? "yes" : "no"));
}
