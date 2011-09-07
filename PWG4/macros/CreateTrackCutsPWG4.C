AliESDtrackCuts *CreateTrackCutsPWG4(Int_t cutMode) {

  //
  // Macro to create track cuts for PWG4 Jet analysis
  // User can select a specific set by indicating cutMode
  // cutMode has 8 digits: first 4 digits additional cuts, last 4 digits standard cuts
  //                       additional cuts are variations of standard cuts (used for hybrid track selection and QA)
  // Numbering starts from 1000 for standard and additional cut numbers

  AliESDtrackCuts *trackCuts  = new AliESDtrackCuts("AliESDtrackCuts");

  TString tag;

  Int_t mod = 10000;

  Bool_t bStdCutsDefined = kFALSE;


  //_____________________________________________________________________
  //                     STANDARD CUTS

  //Get standard cuts: last 4 digits of cutMode
  Int_t stdCutMode = cutMode%mod;

  if(stdCutMode == 1000) {

    bStdCutsDefined = kTRUE;

    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
    trackCuts->SetMinNCrossedRowsTPC(120);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    trackCuts->SetMaxChi2PerClusterITS(36);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 1e10);

    tag = "Global track RAA analysis QM2011 + Chi2ITS<36";

  }

  if(stdCutMode == 1001) {

    bStdCutsDefined = kTRUE;

    // TPC  
    trackCuts->SetMinNClustersTPC(90);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

 
    tag = "Global tracks jet analysis with ITSrefit and NclsIter1=90, noSPD requirement";

  }

  if(stdCutMode == 1002) {

    bStdCutsDefined = kTRUE;

    // TPC  
    trackCuts->SetMinNClustersTPC(80);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

 
    tag = "Global tracks jet analysis with ITSrefit and Ncls=80, noSPD requirement";

  }

  if(stdCutMode == 1003) {

    bStdCutsDefined = kTRUE;

    // tight global tracks
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,1);
    trackCuts->SetMinNClustersTPC(0);
    trackCuts->SetMinNCrossedRowsTPC(120);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.1);// essentially switches it off
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    trackCuts->SetMaxChi2PerClusterITS(36);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);

    tag = "Global tracks ITSTPC2010 + NCrossedRows + loose ITS";

  }
  
  if(stdCutMode == 1004) {

    bStdCutsDefined = kTRUE;

    // TPC  
    trackCuts->SetMinNClustersTPC(70);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);
 
    tag = "Global tracks jet analysis with ITSrefit and NclsIter1=70, noSPD requirement";

  }
  if(stdCutMode == 1005) {

    bStdCutsDefined = kTRUE;

    // TPC  
    trackCuts->SetMinNClustersTPC(70);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 1E+15.);
 
    tag = "Global tracks jet analysis with ITSrefit and NclsIter1=70, noSPD requirement, no upper pt cut";

  }



  if(stdCutMode == 2000) {

    bStdCutsDefined = kTRUE;

    // TPC  
    trackCuts->SetMinNClustersTPC(90);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

 
    tag = "Global tracks jet analysis, loose cuts, NClsIter1=90, no ITS requirements";

  }

  if(stdCutMode == 2001) {

    bStdCutsDefined = kTRUE;

    trackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); 
    trackCuts->SetMinNClustersTPC(70);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

 
    tag = "TPConly track cuts, loose cuts, NCls=70, no ITS requirements";

  }

  if(stdCutMode == 2002) {

    bStdCutsDefined = kTRUE;

    trackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); 
    trackCuts->SetMinNClustersTPC(0);
    trackCuts->SetMinNCrossedRowsTPC(120);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.1);// essentially switches it off  

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

 
    tag = "TPConly track cuts, loose cuts, NCls=70, no ITS requirements";

  }

  if(!bStdCutsDefined) {
    printf("last 4 digits do not represent a predefined set of standard cuts. Returning 0\n");
    return 0;

  }


  //_____________________________________________________________________
  //                     ADDITIONAL CUTS

  //Get additional cut mode: first 4 digits of cutMode
  Int_t addCutMode = (int)((float)cutMode/(float)mod);

  if(addCutMode == 1000) {

    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
 
    tag += " + additonal: SPD any requirement";

  }

  if(addCutMode == 1001) {

    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
 
    tag += " + additional: w/o hits in SPD";

  }

  if(addCutMode == 1002) {

    trackCuts->SetMaxChi2PerClusterITS(1E10);

    tag += " + additional: maxITSChi2=1e10";

  }

  if(addCutMode == 1003) {

    trackCuts->SetMinNClustersTPC(0);
    trackCuts->SetMinNCrossedRowsTPC(0);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.);

    tag += " + additional: minClusters=0 minCrossedRows=0 minCrossedRowsOverFindable=0";

  }

  if(addCutMode == 1004) {

    trackCuts->SetRequireITSRefit(kFALSE);

    tag += " + additional: ITSrefit=kFALSE";

  }

  Printf("Created track cuts for: %s", tag.Data());

  return trackCuts;

}
