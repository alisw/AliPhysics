//=========================================================================//
//                                                                         //
//           ESD track Cuts the for Particle Ratio Fluctuation Study       //
//              Author: Deepika Rathee  || Satyajit Jenara                 //
//                      drathee@cern.ch || sjena@cern.ch                   //
//                       Thu Dec 19 09:09:38 CET 2013                      //
//                                                                         //
//=========================================================================//
//class AliESDtrackCuts;
// If used for Marta's cuts: no eta not dca accepted
//
AliESDtrackCuts *configureNetChargeTrackCut(const Char_t *tname = "sjenaTracksCf", Int_t imode = 5, Int_t cutMode = 1000, Double_t eta = 1.0, Double_t dcaxy = 5., Double_t dcaz = 5.) {

  Double_t gEtaMin = -1*eta;
  Double_t gEtaMax = eta;

  Double_t gPtMin = 0.10;
  Double_t gPtMax = 100.;

  TString tag = Form(" %d : %d : %4.1f : %4.1f %4.1f =>", imode, cutMode, eta, dcaxy, dcaz);
  AliESDtrackCuts *trackCuts  = new AliESDtrackCuts(Form("sjena_%s",tname),Form("satya_%s",tname));

  if (        imode == 0 ) { // "StandardITSTPCTrackCuts2010"
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    
    tag += " Track Cut Used: StandardITSTPCTrackCuts2010";
  } else if ( imode == 1 ) { // "StandardITSTPCTrackCuts2010no"
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
    
    tag += " Track Cut Used: StandardITSTPCTrackCuts2010 No"; 
  } else if ( imode == 2 ) { // "StandardITSTPCTrackCuts2010trueAndCrossRows"
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, 1);
    
    tag += " Track Cut Used: StandardITSTPCTrackCuts2010trueAndCrossRows"; 

  } else if ( imode == 3 ) { // "StandardTPCOnlyTrackCuts"
    trackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    tag += " Track Cut Used: StandardTPCOnlyTrackCuts";
  } else if ( imode == 4 ) { // "StandardITSSATrackCuts2010"
    trackCuts = AliESDtrackCuts::GetStandardITSSATrackCuts2010();

    tag += " Track Cut Used: StandardITSSATrackCuts2010";

  } else if ( imode == 5 ) { // Initial Cuts used for Net-Charge PRL (Though Results are in AODs) - TPC tracks only in analysis
     trackCuts->SetMinNClustersTPC(80);
     trackCuts->SetMaxChi2PerClusterTPC(4.0);
     trackCuts->SetRequireTPCRefit();
     trackCuts->SetAcceptKinkDaughters(kFALSE);
     trackCuts->SetMaxDCAToVertexXY(3.0);
     trackCuts->SetMaxDCAToVertexZ(3.0);
     trackCuts->SetPtRange(0.3,1.5);
     trackCuts->SetEtaRange(-0.8,0.8);
    
     tag += " Track Cut Used: Initial Cuts used for Net-Charge PRL";

  } else if ( imode == 6) { // wide DCA cut
    trackCuts->SetMinNClustersTPC(80);
    trackCuts->SetMinNClustersITS(2);
    trackCuts->SetMaxChi2PerClusterTPC(4.0);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetRequireITSRefit(kTRUE);
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetMaxDCAToVertexXY(dcaxy);
    trackCuts->SetMaxDCAToVertexZ(dcaz);
    trackCuts->SetPtRange(gPtMin,gPtMax);
    trackCuts->SetEtaRange(gEtaMin, gEtaMax);
    
    tag += " Track Cut Used: user DCA Cuts";

  } else if (imode == 7) { // Modified ITS2010 Cuts - with 

    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
    trackCuts->SetMaxChi2PerClusterITS(36);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    trackCuts->SetMaxChi2TPCConstrainedGlobal(36);
    trackCuts->SetPtRange(gPtMin,gPtMax);
    trackCuts->SetEtaRange(gEtaMin, gEtaMax);

    // Extra
    trackCuts->SetMinNClustersTPC(70);  
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); 
    trackCuts->SetRequireITSRefit(kTRUE);                                                     
    trackCuts->SetDCAToVertex2D(kFALSE);   
    trackCuts->SetMaxDCAToVertexZ(dcaz);
    trackCuts->SetRequireSigmaToVertex(kFALSE);  

    tag += " Track Cut Used: Modified ITS2010 Cuts ";                                             
        
  } else if (imode == 8) { /**** My Standard Cut for TPC : drathee****/
  
    trackCuts->SetMinNCrossedRowsTPC(70);                                             
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);                    
    trackCuts->SetMaxChi2PerClusterTPC(4);                                            
    trackCuts->SetAcceptKinkDaughters(kFALSE);                                        
    trackCuts->SetRequireTPCRefit(kTRUE);                                             
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff); 
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kOff); 
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSSD,AliESDtrackCuts::kOff); 
    trackCuts->SetRequireITSRefit(kTRUE);                                             
    trackCuts->SetMaxChi2PerClusterITS(36);                                           
    
    trackCuts->SetDCAToVertex2D(kFALSE);                                              
    trackCuts->SetRequireSigmaToVertex(kFALSE);                                       
       
    //    trackCuts->SetMaxDCAToVertexXY(dcaxy);
    trackCuts->SetMaxDCAToVertexZ(dcaz);
    trackCuts->SetPtRange(gPtMin,gPtMax);
    trackCuts->SetEtaRange(gEtaMin, gEtaMax);    
    trackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");    // Can use TFormula FixMe 

    tag += " Track Cut Used: My Standard Cut for TPC";                                                                

  } else if (imode == 9) { /*** for Correction jochen ***/
    trackCuts->SetMinNCrossedRowsTPC(70);                                            
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);                   
    trackCuts->SetMaxChi2PerClusterTPC(4);                                           
    trackCuts->SetAcceptKinkDaughters(kFALSE);                                       
    trackCuts->SetRequireTPCRefit(kTRUE);                                            
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kOff);
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSSD,AliESDtrackCuts::kOff);
    trackCuts->SetRequireITSRefit(kTRUE);                                            
    trackCuts->SetMaxChi2PerClusterITS(36);                                          
    
    //  trackCuts->SetDCAToVertex2D(kFALSE);                                             
    // trackCuts->SetRequireSigmaToVertex(kFALSE);                                      

    trackCuts->SetDCAToVertex2D(kTRUE);                                             
    trackCuts->SetRequireSigmaToVertex(kTRUE);                                      
       
    trackCuts->SetMaxDCAToVertexXY(dcaxy);
    trackCuts->SetMaxDCAToVertexZ(dcaz);
    trackCuts->SetPtRange(gPtMin,gPtMax);
    trackCuts->SetEtaRange(gEtaMin, gEtaMax);
    
    tag += " Track Cut Used: Mainly Used for Correction";                                                                
    
  } else if (imode == 10) { /*** LF: Standard Cuts ***/
    trackCuts->SetMinNClustersTPC(70);   
    trackCuts->SetMaxChi2PerClusterTPC(4);                                             
    trackCuts->SetAcceptKinkDaughters(kFALSE);                                         
    trackCuts->SetRequireTPCRefit(kTRUE);   
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); 
    
    trackCuts->SetRequireITSRefit(kTRUE);                                             
    trackCuts->SetMaxChi2PerClusterITS(36);                                           
    
    trackCuts->SetDCAToVertex2D(kFALSE);                                              
    trackCuts->SetRequireSigmaToVertex(kFALSE);                                       
    
    //trackCuts->SetMaxDCAToVertexXY(dcaxy);
    trackCuts->SetMaxDCAToVertexZ(dcaz);
    trackCuts->SetPtRange(gPtMin,gPtMax);
    trackCuts->SetEtaRange(gEtaMin, gEtaMax);
    trackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");  
                                 
    tag += " Track Cut Used: LF Customised Cuts";                                                                
  } else if (imode == 99) { // Marta's Magic
    
    // I am not using this moment check it for future
    // TODO: Impliment the Marta's suggestion of Hybrid Like tracks
    // Taken from JET Configure task : 
    // Magic is done by Marta
    
    // Macro to create track cuts for PWG Jet analysis
    // User can select a specific set by indicating cutMode
    // cutMode has 8 digits: first 4 digits additional cuts, last 4 digits standard cuts
    //                       additional cuts are variations of standard cuts (used for hybrid track selection and QA)
    // Numbering starts from 1000 For standard and additional cut numbers
    
  
    Int_t mod = 10000;
    Bool_t bStdCutsDefined = kFALSE;
    
    //Get standard cuts: last 4 digits of cutMode
    Int_t stdCutMode = cutMode%mod;
    
    if(stdCutMode == 1000) {
      bStdCutsDefined = kTRUE;
      trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
      trackCuts->SetMinNCrossedRowsTPC(120);
      trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      trackCuts->SetMaxChi2PerClusterITS(36);
      trackCuts->SetMaxFractionSharedTPCClusters(0.4);
      trackCuts->SetMaxChi2TPCConstrainedGlobal(36);
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
    
    if(stdCutMode == 1006) {
      
      bStdCutsDefined = kTRUE;
      
      // TPC  
      TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
      trackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
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
      trackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      
      trackCuts->SetRequireSigmaToVertex(kFALSE);
      
      trackCuts->SetEtaRange(-0.9,0.9);
      trackCuts->SetPtRange(0.15, 1E+15.);
      
      tag = "Global tracks jet analysis with ITSrefit and NclsIter1=PtDep, noSPD requirement, no upper pt cut, golden chi2";
      
    }
    
    if(stdCutMode == 1007) {
      
      bStdCutsDefined = kTRUE;
      
      trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
      //accept secondaries
      trackCuts->SetMaxDCAToVertexXY(2.4);
      trackCuts->SetMaxDCAToVertexZ(3.2);
      trackCuts->SetDCAToVertex2D(kTRUE);
      
      //
      trackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      
      trackCuts->SetEtaRange(-0.9,0.9);
      trackCuts->SetPtRange(0.15, 1E+15.);
      
      tag = "Global tracks with AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE)";
      
    }
    
    if(stdCutMode == 1008) {
      
      bStdCutsDefined = kTRUE;
      
      trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
      //accept secondaries
      trackCuts->SetMaxDCAToVertexXY(2.4);
      trackCuts->SetMaxDCAToVertexZ(3.2);
      trackCuts->SetDCAToVertex2D(kTRUE);
      
      //
      trackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
      
      trackCuts->SetMaxFractionSharedTPCClusters(0.4);
      
      tag = "Global tracks 2011 with NCrossedRows cut";
      
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
    // trackCuts->SetRequireTPCRefit(kTRUE);
      trackCuts->SetMinNClustersTPC(70);
      
      trackCuts->SetEtaRange(-0.9,0.9);
      trackCuts->SetPtRange(0.15, 100.);
      
      tag = "TPConly track cuts, loose cuts, NCls=70, no ITS requirements";
    }
    
    if(stdCutMode == 2002) {
      
      bStdCutsDefined = kTRUE;
      
      trackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); 
      //  trackCuts->SetRequireTPCRefit(kTRUE);
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
    
    if(addCutMode == 1005) {
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
      tag += " + additional: no SPD requirement (kOff)";
      
    }
    
    Printf("Created track cuts for: %s", tag.Data());
  }

  else {
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
    Printf(" Default 2010 Cuts ");
  }
    
  

  Printf("%s", tag.Data());
  
  trackCuts->Print();

  return trackCuts;
  
}
