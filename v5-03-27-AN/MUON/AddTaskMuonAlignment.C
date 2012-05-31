/// \ingroup macros
/// \file AddTaskMuonAlignment.C
/// \brief Macro to add an AliMUONAlignmentTask to an analysis train
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN

AliMUONAlignmentTask4 *AddTaskMuonAlignment(bool readrecs, bool doalign, bool writerecs, const char *name = "AliMUONAlignmentTask", const char *newalignocdb = "local://ReAlignOCDB", const char *oldalignocdb = "alien://Folder=/alice/data/2010/OCDB", const char *defaultocdb = "alien://Folder=/alice/data/2011/OCDB", const char *geofilename = "geometry.root")
{
/// Creates a Muon Alignment task and adds it to the analysis manager.

	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		::Error("AddTaskMuonAlignment", "No analysis manager to connect to.");
		return NULL;
	}   
   
	// This task requires an ESD input handler.
	// Check this using the analysis manager.
	//===============================================================================
	TString type = mgr->GetInputEventHandler()->GetDataType();
	if (readrecs) {		
		if (!type.Contains("AOD")) {
			::Error("AddTaskMuonAlignment", "Alignment task is asked to read records but input data type is not AOD! Game over!");
			return NULL;
		} 
	} else {
		if (!type.Contains("ESD")) {
			::Error("AddTaskMuonAlignment", "Alignment task is not asked to read records but input data type is not ESD! Game over!");
			return NULL;
		} 
	}

	// Create the task, add it to the manager and configure it.
	//===========================================================================   
	// Muons
	AliMUONAlignmentTask *muonalign = new AliMUONAlignmentTask(name, newalignocdb, oldalignocdb,  defaultocdb, geofilename);
	muonalign->SetLoadOCDBOnce(true);
	muonalign->SetReadRecords(readrecs);
  muonalign->SetDoAlignment(doalign);
  muonalign->SetWriteRecords(writerecs);
	mgr->AddTask(muonalign);
	
//    // Cuts on primary tracks
//    AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
//    esdTrackCutsL->SetMinNClustersTPC(50);

//    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
//    trackFilter->AddCuts(esdTrackCutsL);

//    muonalign->SetTrackFilter(trackFilter);


	// Create ONLY the output containers for the data produced by the task.
	// Get and connect other common input/output containers via the manager as below
	//==============================================================================
	mgr->ConnectInput(muonalign,  0, mgr->GetCommonInputContainer());
	if (doalign) {
		AliAnalysisDataContainer *listOut = mgr->CreateContainer("output", TList::Class(), AliAnalysisManager::kOutputContainer, "measShifts.root");
		mgr->ConnectOutput(muonalign,  1, listOut);
	}
	if (writerecs) {
		mgr->ConnectOutput(muonalign,  0, mgr->GetCommonOutputContainer());
	}
	return muonalign;
}   
