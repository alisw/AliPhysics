AliAnalysisTaskV0Mnorm* AddTaskV0Mnorm(
       const char* mcpariclearraynamePartMC = "mcparticles", //name of mcparticle TClonesArray array for MC particle level jets
       Bool_t      useVertexCut       = kTRUE,  // vertex cut
       Bool_t      usePileUpCut       = kTRUE, // discard pile up event
       const char* suffix = ""                              //SUBWAGON has to be the l
){

   return AliAnalysisTaskV0Mnorm::AddTaskV0Mnorm(
       mcpariclearraynamePartMC,
       useVertexCut,
       usePileUpCut,
       suffix 
          ); 
}
