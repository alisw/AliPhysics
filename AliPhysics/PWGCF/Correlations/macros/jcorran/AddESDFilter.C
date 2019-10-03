AliAnalysisFilter *AddESDFilter()
{
  // Cuts on primary tracks
  AliESDtrackCuts* esdTrackCutsL = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();

  // ITS stand-alone tracks
  AliESDtrackCuts* esdTrackCutsITSsa = new AliESDtrackCuts("ITS stand-alone Track Cuts", "ESD Track Cuts");
  esdTrackCutsITSsa->SetRequireITSStandAlone(kTRUE);

  // Pixel OR necessary for the electrons
  AliESDtrackCuts *itsStrong = new AliESDtrackCuts("ITSorSPD", "pixel requirement for ITS");
  itsStrong->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);


  // PID for the electrons
  AliESDpidCuts *electronID = new AliESDpidCuts("Electrons", "Electron PID cuts");
  electronID->SetTPCnSigmaCut(AliPID::kElectron, 3.);

  // standard cuts with very loose DCA
  AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE); 
  esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
  esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
  esdTrackCutsH->SetDCAToVertex2D(kTRUE);

  // standard cuts with tight DCA cut
  AliESDtrackCuts* esdTrackCutsH2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

  // standard cuts with tight DCA but with requiring the first SDD cluster instead of an SPD cluster
  // tracks selected by this cut are exclusive to those selected by the previous cut
  AliESDtrackCuts* esdTrackCutsH3 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(); 
  esdTrackCutsH3->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  esdTrackCutsH3->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);

  // TPC only tracks: Optionally enable the writing of TPConly information
  // constrained to SPD vertex in the filter below
  AliESDtrackCuts* esdTrackCutsTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  esdTrackCutsTPCOnly->SetMinNClustersTPC(70);

  // R_AA cut with recent change in AliESDtrackCuts in trunk aliroot (111026)
  AliESDtrackCuts* esdTrackCutsRaa = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1); // 1,2,3,4,5,6,7
  esdTrackCutsRaa->SetMinNCrossedRowsTPC(120); // 3 (70 set in above cut)
  esdTrackCutsRaa->SetMaxChi2PerClusterITS(36.); // 8
  esdTrackCutsRaa->SetMaxChi2TPCConstrainedGlobal(36.); // 9
  esdTrackCutsRaa->SetMaxFractionSharedTPCClusters(0.4); // 10

  // Compose the filter
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  // 1
  trackFilter->AddCuts(esdTrackCutsL);
  // 2
  trackFilter->AddCuts(esdTrackCutsITSsa);
  // 4
  trackFilter->AddCuts(itsStrong);
  itsStrong->SetFilterMask(1);        // AND with Standard track cuts 
  // 8
  trackFilter->AddCuts(electronID);
  electronID->SetFilterMask(4);       // AND with Pixel Cuts
  // 16
  trackFilter->AddCuts(esdTrackCutsH);
  // 32
  trackFilter->AddCuts(esdTrackCutsH2);
  // 64
  trackFilter->AddCuts(esdTrackCutsH3);
  // 128 , 1 << 7
  trackFilter->AddCuts(esdTrackCutsTPCOnly);

  //========================================
  //
  //        Add Additional Cuts here
  //
  //========================================

  trackFilter->AddCuts(esdTrackCutsRaa);

  return trackFilter;
}

