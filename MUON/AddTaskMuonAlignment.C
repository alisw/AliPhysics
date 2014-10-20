/// \ingroup macros
/// \file AddTaskMuonAlignment.C
/// \brief Macro to add an AliMUONAlignmentTask to an analysis train
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN
/// \author Hugo Pereira Da Costa, CEA/Saclay - Irfu/SPhN

/**
AliMUONAlignmentTask calculates alignment correction for the muon spectrometer
based on reconstructed tracks, either with or without magnetic field. It uses AliMillePede2
internally for alignment, and AliMillePedeRecords for storing the derivatives.

The task has two ways of operation

1/ Full mode:
This corresponds to flags: readRecords = kFALSE, doAlignment = kTRUE, writeRecords = kFALSE
- reads tracks from an ESD,
- calculates all derivatives needed to minimize both the tracks
  and the alignment parameters
- fill and invert the (large) matrix corresponding to the global chisquare minimization
- write the resulting alignment parameters, on top of the initiali geometry, in a new CDB

Note that writeRecords can also be set to kTRUE. This will output all derivatives to a special
branch in an AOD for re-running the alignment with no need to process the tracks (see below)

2/ Split mode:
The task must run twice.
First pass with flags: readRecords = kFALSE, writeRecords = kTRUE, doAlignment = kFALSE
This
- reads tracks from an ESD,
- calculates all derivatives needed to minimize both the tracks
  and the alignment parameters
- stores them in a special branch of the output AOD, named "records",
  and containing a TClonesArray of objects typed "AliMillePedeRecord"

Second pass with flags: readRecords = kTRUE, doAlignment = kTRUE
- reads the AliMillePedeRecords stored in the AOD (or a collection of AODs)
- fill and invert the (large) matrix corresponding to the global chisquare minimization
- write the resulting alignment parameters, on top of the initiali geometry, in a new CDB

The split mode is usefull for running on the grid:
- the first pass can be done with multiple jobs (one per ESD chunk) in parallel.
- the second pass must be done one on all created AODs, in a single shot, possibly locally

When performing the alignment (doAlignment=kTRUE), one can "fix", "group", or "constrain"
detector's alignment parameters in order to ease the inversion.

One can redo the second pass as many times as needed, changing these "fixed/group/constrained" parameters,
without re-processing the ESD tracks.
*/

#include "AliAnalysisManager.h"
#include "AliMUONAlignmentTask.h"
#include "AliVEventHandler.h"

AliMUONAlignmentTask *AddTaskMuonAlignment(
  TString oldAlignmentOCDB,
  TString newAlignmentOCDB,
  Bool_t doAlignment = kTRUE,
  Bool_t writeRecords = kTRUE,
  Bool_t readRecords = kFALSE
 )
{

  /// Creates a Muon Alignment task and adds it to the analysis manager.

  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *analysisManager = AliAnalysisManager::GetAnalysisManager();
  if( !analysisManager )
  {
    ::Error("AddTaskMuonAlignment", "No analysis manager to connect to.");
    return NULL;
  }

  // get input event handler and check type
  TString type = analysisManager->GetInputEventHandler()->GetDataType();
  if( readRecords )
  {

    // when reading records, AOD are required
    if (!type.Contains( "AOD" ) )
    {
      Error("AddTaskMuonRefit", "AOD input handler needed!");
      return NULL;
    }

  } else {

    // ESDs are required otherwise
    if( !type.Contains( "ESD" ) )
    {
      Error("AddTaskMuonRefit", "AOD input handler needed!");
      return NULL;
    }

  }

  // Create the task, add it to the manager and configure it.
  AliMUONAlignmentTask *muonAlign = new AliMUONAlignmentTask( "AliMUONAlignmentTask" );
  muonAlign->SetOldAlignStorage( oldAlignmentOCDB );
  muonAlign->SetNewAlignStorage( newAlignmentOCDB );
  muonAlign->SetLoadOCDBOnce( kTRUE );
  muonAlign->SetReadRecords( readRecords );
  muonAlign->SetDoAlignment( doAlignment );
  muonAlign->SetWriteRecords( writeRecords );
  muonAlign->SetMergeAlignmentCDBs( kTRUE );

  analysisManager->AddTask(muonAlign);

  // connect input
  analysisManager->ConnectInput(muonAlign,  0, analysisManager->GetCommonInputContainer());

  // when writting records, also connect output
  if( writeRecords )
  { analysisManager->ConnectOutput(muonAlign,  0, analysisManager->GetCommonOutputContainer()); }

  // return created task
  return muonAlign;

}
