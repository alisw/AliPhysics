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

/* $Header$ */

/////////////////////////////////////////////////////////////////////////
//                                                                     //
//       Prototype ESD class                                           //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "Riostream.h"

#include "AliESD.h"

ClassImp(AliESD)

//_______________________________________________________________________
AliESD::AliESD():
  fEventNumber(0),
  fRunNumber(0),
  fTrigger(0),
  fRecoVersion(0),
  fBitDDL(0),
  fNSecVertex(0),
  fNParticipants(0),
  fNPartError(0),
  fNElectron(0),
  fNMuons(0),
  fNPions(0),
  fNKaons(0),
  fNProtons(0),
  fNPHOSPhotons(0),
  fNPHOSNeutrons(0),
  fNPHOSCCluster(0),
  fNEMCALCluster(0),
  fNPMDCluster(0),
  fTMaxClusterEnergy(0),
  fTMaxPCharged(0),
  fNCharged(0),
  fTotTranEnergy(0),
  fESDVertex(),
  fSecVertex(0),
  fNonAssTrack(0),
  fPhoton(0),
  fNeutron(0),
  fEMCALCluster(0),
  fPMDCluster(0)
{
  Info("def ctor","Has been called\n");
}

//_______________________________________________________________________
AliESD::AliESD(const AliESD &esd):
  TObject(esd),
  fEventNumber(0),
  fRunNumber(0),
  fTrigger(0),
  fRecoVersion(0),
  fBitDDL(0),
  fNSecVertex(0),
  fNParticipants(0),
  fNPartError(0),
  fNElectron(0),
  fNMuons(0),
  fNPions(0),
  fNKaons(0),
  fNProtons(0),
  fNPHOSPhotons(0),
  fNPHOSNeutrons(0),
  fNPHOSCCluster(0),
  fNEMCALCluster(0),
  fNPMDCluster(0),
  fTMaxClusterEnergy(0),
  fTMaxPCharged(0),
  fNCharged(0),
  fTotTranEnergy(0),
  fESDVertex(),
  fSecVertex(0),
  fNonAssTrack(0),
  fPhoton(0),
  fNeutron(0),
  fEMCALCluster(0),
  fPMDCluster(0)
{
}

ClassImp(AliESDVertex)

//_______________________________________________________________________
AliESDVertex::AliESDVertex():
  fNPrimary(0),
  fCoordinates(3),
  fErrorMatrix(6),
  fPrimaryTracks(0),
  fEffectiveMass(0),
  fEffectiveMassError(0)
{
  cout << "ESDVertex def ctor" << endl;
}

//_______________________________________________________________________
AliESDVertex::AliESDVertex(const AliESDVertex &esdv):
  TObject(esdv),
  fNPrimary(0),
  fCoordinates(0),
  fErrorMatrix(0),
  fPrimaryTracks(0),
  fEffectiveMass(0),
  fEffectiveMassError(0)
{
}

ClassImp(AliESDTrack)

//_______________________________________________________________________
AliESDTrack::AliESDTrack() :
  fTrackID(0),
  fPVertex(5),
  fPEVertex(15),
  fPFMeasPoint(6),
  fPFMeasPointErr(15),
  fPLMeasPoint(6),
  fPLMeasPointErr(15),
  fTrackLength(0),
  fTrackLengthErr(0),
  fStopVertex(0),
  fNPointsITS(0),
  fNPointsTPC(0),
  fNPointsTRD(0),
  fMeanResITS(0),
  fMeanResTPC(0),
  fMeanResTRD(0),
  fGlobalChi2(0),
  fParticleType(0),
  fPIDprobPi(0),
  fPIDprobK(0),
  fPIDprobP(0),
  fPIDprobE(0)
{
  cout << "ESDTrack def ctor" << endl;
}

//_______________________________________________________________________
AliESDTrack::AliESDTrack(const AliESDTrack &esdt):
  TObject(esdt),
  fTrackID(0),
  fPVertex(0),
  fPEVertex(0),
  fPFMeasPoint(0),
  fPFMeasPointErr(0),
  fPLMeasPoint(0),
  fPLMeasPointErr(0),
  fTrackLength(0),
  fTrackLengthErr(0),
  fStopVertex(0),
  fNPointsITS(0),
  fNPointsTPC(0),
  fNPointsTRD(0),
  fMeanResITS(0),
  fMeanResTPC(0),
  fMeanResTRD(0),
  fGlobalChi2(0),
  fParticleType(0),
  fPIDprobPi(0),
  fPIDprobK(0),
  fPIDprobP(0),
  fPIDprobE(0)
{
}

