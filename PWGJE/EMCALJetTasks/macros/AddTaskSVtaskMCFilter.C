AliAnalysisTaskSVtaskMCFilter* AddTaskSVtaskMCFilter(const char* trkcontname   = "tracks", const char* outtrk = "mytracks", Bool_t fFilterTracks = kTRUE, const char* clscontname = "caloClusters", const char* outcls = "myclusters", Bool_t fFilterClusters = kFALSE){

   //Steering macro for filter of detector level MC tracks. Removes EPOS part of the event


   return  AliAnalysisTaskSVtaskMCFilter::AddTaskSVtaskMCFilter(trkcontname, outtrk, fFilterTracks, clscontname, outcls, fFilterClusters);

}
