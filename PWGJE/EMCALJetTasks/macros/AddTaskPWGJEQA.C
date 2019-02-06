// AddTaskPWGJEQA.C

AliAnalysisTaskPWGJEQA* AddTaskPWGJEQA(
                                       const char* ntracks            = "usedefault",
                                       const char* nclusters          = "usedefault",
                                       const char* ncells             = "usedefault",
                                       const char *nGenLev            = "mcparticles",
                                       Bool_t      doTrackQA          = kTRUE,
                                       Bool_t      doCaloQA           = kTRUE,
                                       Bool_t      doJetQA            = kTRUE,
                                       Bool_t      doEventQA          = kTRUE,
                                       Double_t    trackPtCut         = 0.15,
                                       Double_t    clusECut           = 0.30,
                                       const char* suffix             = ""
                                       )
{
	  return AliAnalysisTaskPWGJEQA::AddTaskPWGJEQA(ntracks,nclusters,ncells,nGenLev,
			                                        doTrackQA,doCaloQA,doJetQA,doEventQA,
			                                        trackPtCut,clusECut,suffix);

}
