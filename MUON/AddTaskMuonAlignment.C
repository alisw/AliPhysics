/// \ingroup macros
/// \file AddTaskMuonAlignment.C
/// \brief Macro to add an AliMUONAlignmentTask to an analysis train
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN
/// \author Hugo Pereira Da Costa, CEA/Saclay - Irfu/SPhN

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
  muonAlign->SetUnbias( kFALSE );

  analysisManager->AddTask(muonAlign);

  // connect input
  analysisManager->ConnectInput(muonAlign,  0, analysisManager->GetCommonInputContainer());

  // when writting records, also connect output
  if( writeRecords )
  { analysisManager->ConnectOutput(muonAlign,  0, analysisManager->GetCommonOutputContainer()); }

  // return created task
  return muonAlign;

}
