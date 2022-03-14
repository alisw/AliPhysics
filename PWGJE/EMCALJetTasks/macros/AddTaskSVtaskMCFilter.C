AliAnalysisTaskSVtaskMCFilter* AddTaskSVtaskMCFilter(const char* trkcontname   = "tracks", const char* outtrk = "mytracks", Bool_t fFilterTracks = kTRUE){

   //Steering macro for filter of detector level MC tracks. Removes EPOS part of the event

	
   return  AliAnalysisTaskSVtaskMCFilter::AddTaskSVtaskMCFilter(trkcontname, outtrk, fFilterTracks);

}
