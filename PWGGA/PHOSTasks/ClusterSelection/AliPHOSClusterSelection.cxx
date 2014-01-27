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

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultiInputEventHandler.h"
#include "AliVEvent.h"

#include "AliPHOSClusterSelection.h"

AliPHOSClusterSelection::AliPHOSClusterSelection()
  : fMinChargedParticleTrackDistance(-1.), 
    fNotUnfolded(false),
    fMaxDispR2(-1.),
    fMaxDispCoreR2(-1.),
    fMaxTOF(-1.)
{
  // Defaults to the most lenient selection allowable
	return;
}

AliPHOSClusterSelection::~AliPHOSClusterSelection()
{
}

Bool_t AliPHOSClusterSelection::IsSelected(AliVCluster* cluster) const
{
  return IsSelectedCPV(cluster)
    && IsSelectedUnfolded(cluster)
    && IsSelectedDisp(cluster)
    && IsSelectedDispCore(cluster)
    && IsSelectedTOF(cluster);
}

Bool_t IsSelectedCPV(AliVCluster* cluster) const
{
  if( 0 > SetMinChargedParticleTrackDistance )
    return true; 
  //TODO: implement positive case
}

AliPHOSClusterSelection* AliPHOSClusterSelection::SetMinChargedParticleTrackDistance(Float_t distance)
{
  // 'distance' set the minimal allowable distance between the cluster
  // and the nearest extrapolated track.
  // If 'distance' is negative, then all clusters are sellected, the selection
  // being "not applied" or "disabled".
  
  fMinChargedParticleTrackDistance = distance;
}

TString AliPHOSClusterSelection::ToString() const
{
  // returns a string an quasi-unique string for whatever selection 
  // parameters the instance contains. The uniqueness of the string
  // is limited by the precision given in the formatting of the string.
  // Take care that the precision is sufficient for your needs.

  return TString::Format("%f_%i_%f_%f_%f",
			 fMinChargedParticleTrackDistance,
			 fNotUnfolded,
			 fMaxDispR2,
			 fMaxDispCoreR2,
			 fMaxTOF
			 );
}


Float_t AliPHOSClusterSelection::SetMinChargedParticleTrackDistance(const TString& string)
{
  TObjArray * objarray = string.Tokenize("_");
  Float_t flt = objarray->At(0)->Atof();
  delete objarray;
  return flt;
}


AliVEvent* AliPHOSClusterSelection::GetCurrentEvent() const
{
  // Hackish way of getting the current event.
  // Its probably not appropriate to call this function outside of
  // AliAnalysisTaskSE::UserExec
  
  AliAnalysisManager* analysisManager = dynamic_cast<AliAnalysisManager*>(AliAnalysisManager::GetAnalysisManager());
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(analysisManager->GetInputEventHandler());
  AliMultiInputEventHandler *multiInputHandler = dynamic_cast<AliMultiInputEventHandler *>(fInputHandler);
  if (multiInputHandler)
    inputHandler = dynamic_cast<AliInputEventHandler *>(multiInputHandler->GetFirstInputEventHandler());
  
  AliVEvent* inputEvent = dynamic_cast<AliVEvent*>(inputHandler->GetEvent());
  if( ! inputEvent ) 
    AliError("Was not able to retrieve event!");
  
  return inputEvent;
}
