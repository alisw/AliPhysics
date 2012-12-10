AliESDtrackCuts *RsnQualityCut(TString cut="pp_LHC11_p4_120") {


   // For RSN analysis, we select Primaries
   Bool_t selPrimaries = kTRUE;

   Printf("RsnQualityCut : %s",cut.Data());
   AliESDtrackCuts *esdTrackCuts = 0;
   if (cut.Contains("pp_LHC11a_p4")) {
      //esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(Bool_t selPrimaries=kTRUE, Int_t clusterCut=1);
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);

      // std AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1)
      esdTrackCuts->SetMaxChi2PerClusterTPC(4);
      esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
      esdTrackCuts->SetRequireTPCRefit(kTRUE);
      // ITS

      esdTrackCuts->SetRequireITSRefit(kTRUE);
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                             AliESDtrackCuts::kAny);
      if(selPrimaries) {
         // 7*(0.0015+0.0050/pt^1.1)
         esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
         esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      }
      esdTrackCuts->SetMaxDCAToVertexZ(2);
      esdTrackCuts->SetDCAToVertex2D(kFALSE);
      esdTrackCuts->SetRequireSigmaToVertex(kFALSE);

      // additional cuts by FilterBit 10

      if (!cut.CompareTo("pp_LHC11a_p4_120")) esdTrackCuts->SetMinNCrossedRowsTPC(120);
      if (!cut.CompareTo("pp_LHC11a_p4_70")) esdTrackCuts->SetMinNCrossedRowsTPC(70);

      esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      esdTrackCuts->SetMaxChi2PerClusterITS(36);
      esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
      esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      esdTrackCuts->SetEtaRange(-0.9,0.9);
      esdTrackCuts->SetPtRange(0.15, 1e10);
   } else if (cut.Contains("pp_LHC11a_p3")) {
      //AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(Bool_t selPrimaries=kTRUE, Int_t clusterCut=0);
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,0);

      // std AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,0);
      esdTrackCuts->SetMinNClustersTPC(70);
      esdTrackCuts->SetMaxChi2PerClusterTPC(4);
      esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
      esdTrackCuts->SetRequireTPCRefit(kTRUE);
      // ITS
      esdTrackCuts->SetRequireITSRefit(kTRUE);
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                             AliESDtrackCuts::kAny);
      if(selPrimaries) {
         // 7*(0.0026+0.0050/pt^1.01)
         esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
         esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      }
      esdTrackCuts->SetMaxDCAToVertexZ(2);
      esdTrackCuts->SetDCAToVertex2D(kFALSE);
      esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
      esdTrackCuts->SetMaxChi2PerClusterITS(36)
   }

   return esdTrackCuts;
}
