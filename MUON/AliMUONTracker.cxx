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

/// \class AliMUONTracker
/// Interface class for use of global tracking framework;
/// reconstruct tracks from recpoints
///
/// \author Christian Finck, SUBATECH Nantes

#include "AliMUONTracker.h"
#include "AliMUONTrackReconstructorK.h"
#include "AliMUONTrackReconstructor.h"
#include "AliMUONRecData.h"
#include "AliLog.h"

//_____________________________________________________________________________
AliMUONTracker::AliMUONTracker()
  : AliTracker(),
    fTriggerCircuit(0x0),
    fMUONData(0x0),
    fTrackReco(0x0)
{
  /// constructor

}
//_____________________________________________________________________________
AliMUONTracker::~AliMUONTracker()
{
  /// dtr
  delete fTrackReco;
}

//_____________________________________________________________________________
void AliMUONTracker::SetOption(Option_t* option)
{
  /// set reconstructor class

  if (!fMUONData)
    AliError("MUONData not defined");

  if (!fTriggerCircuit)
    AliError("TriggerCircuit not defined");

  if (strstr(option,"Original")) 
    fTrackReco = new AliMUONTrackReconstructor(fMUONData);
  else if (strstr(option,"Combi")) 
    fTrackReco = new AliMUONTrackReconstructorK(fMUONData,"Combi");
  else 
    fTrackReco = new AliMUONTrackReconstructorK(fMUONData,"Kalman");

  fTrackReco->SetTriggerCircuit(fTriggerCircuit);

}

//_____________________________________________________________________________
Int_t AliMUONTracker::Clusters2Tracks(AliESD* /*esd*/)
{

  /// clusters2Tracks method
  /// in general tracking framework
   
  // open TClonesArray for reading
  fMUONData->SetTreeAddress("TC,RC");

  // open for writing
  // trigger branch
  fMUONData->MakeBranch("RL"); //trigger track
  fMUONData->SetTreeAddress("RL");
  fTrackReco->EventReconstructTrigger();
  fMUONData->Fill("RL");

  // tracking branch
  fMUONData->MakeBranch("RT"); //track
  fMUONData->SetTreeAddress("RT");
  fTrackReco->EventReconstruct();
  fMUONData->Fill("RT");

  fMUONData->ResetRecTracks();
  fMUONData->ResetRecTriggerTracks();

 
  return kTRUE;
}
