AliESDtrackCuts *RsnQualityCut(TString cut="pp_LHC11_p4_120") {


   // For RSN analysis, we select Primaries
   Bool_t selPrimaries = kTRUE;

   Printf("RsnQualityCut : %s",cut.Data());
   AliESDtrackCuts *esdTrackCuts = 0;
   if (cut.Contains("pp_LHC11a_p4_AOD113")) {
      //esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(Bool_t selPrimaries=kTRUE, Int_t clusterCut=1);
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
      if (!cut.CompareTo("pp_LHC11a_p4_AOD113_120")) esdTrackCuts->SetMinNCrossedRowsTPC(120);
      if (!cut.CompareTo("pp_LHC11a_p4_AOD113_70")) esdTrackCuts->SetMinNCrossedRowsTPC(70);

      esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      esdTrackCuts->SetMaxChi2PerClusterITS(36);
      esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
      esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      esdTrackCuts->SetEtaRange(-0.9,0.9);
      esdTrackCuts->SetPtRange(0.15, 1e10);

   } else if (cut.Contains("pp_LHC11a_p3_AOD067_bit4")) {
      //AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(Bool_t selPrimaries=kTRUE, Int_t clusterCut=0);
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,0);
      // additional to std cuts
      esdTrackCuts->SetMaxDCAToVertexXY(2.4);
      esdTrackCuts->SetMaxDCAToVertexZ(3.2);
      esdTrackCuts->SetDCAToVertex2D(kTRUE);

   } else if (cut.Contains("pp_LHC11a_p3_AOD067_bit5")) {
      //AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(Bool_t selPrimaries=kTRUE, Int_t clusterCut=0);
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,0);
   } else if (cut.Contains("STD2010_PRIMARY")) {
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,0);
   } else if (cut.Contains("STD2010_SECONDARY")) {
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,0);
   } else if (cut.Contains("STD2011_PRIMARY")) {
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
   } else if (cut.Contains("STD2011_SECONDARY")) {
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
   } else if (cut.Contains("STD2011_PRIMARY_NCLSTTPC")) {
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,0);
   } else if (cut.Contains("STD2011_SECONDARY_NCLSTTPC")) {
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,0);
   }

   if (cut.BeginsWith("STD")) {
      // DCAXY: 3.5 - 14 sigma
      if (cut.Contains("DCAXY7S")) esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.035/pt^1.01");
      else if (cut.Contains("DCAXY6S")) esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0156+0.03/pt^1.01");
      else if (cut.Contains("DCAXY5S")) esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.013+0.025/pt^1.01");

      // DCAZ: 0.1 - 2 cm
      if (cut.Contains("DCAZ20")) esdTrackCuts->SetMaxDCAToVertexZ(2);
      else if (cut.Contains("DCAZ01")) esdTrackCuts->SetMaxDCAToVertexZ(0.1);
      // MinNClustersTPC: 50-70
      if (cut.Contains("NCLSTTPC50")) esdTrackCuts->SetMinNClustersTPC(50);
      else if (cut.Contains("NCLSTTPC70")) esdTrackCuts->SetMinNClustersTPC(70);
      else if (cut.Contains("NCLSTTPC80")) esdTrackCuts->SetMinNClustersTPC(80);
      // Chi2 in TPC: 4-6
      if (cut.Contains("CHI2TPC04")) esdTrackCuts->SetMaxChi2PerClusterITS(4);
      else if (cut.Contains("CHI2TPC6")) esdTrackCuts->SetMaxChi2PerClusterITS(6);
      // Chi2 in ITS: 36-100
      if (cut.Contains("CHI2ITS036")) esdTrackCuts->SetMaxChi2PerClusterITS(36);
      else if (cut.Contains("CHI2ITS100")) esdTrackCuts->SetMaxChi2PerClusterITS(100);
   }

   return esdTrackCuts;
}
