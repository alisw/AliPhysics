/* $Id: $ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon - jet (JETAN) correlation analysis with ESDs
//
// Author : Gustavo Conesa Balbastre (INFN-LNF)
//------------------------------------

AliAnaPartCorrMaker*  ConfigAnalysis()
{
  //
  // Configuration goes here
  // 
  printf("======================== \n");
  printf("ConfigAnalysis() \n");
  printf("======================== \n");
  
  
  //Detector Fidutial Cuts
  AliFidutialCut * fidCut = new AliFidutialCut();
  fidCut->DoCTSFidutialCut(kTRUE) ;
  fidCut->DoEMCALFidutialCut(kTRUE) ;
  fidCut->DoPHOSFidutialCut(kTRUE) ;
  
  fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
  fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  fidCut->SetSimplePHOSFidutialCut(0.13,220.,320.);
 
  fidCut->Print("");
  
  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  AliCaloTrackReader *reader = new AliCaloTrackESDReader();
  reader->SetDebug(-1);

  //Switch on or off the detectors information that you want
  reader->SwitchOnEMCAL();
  reader->SwitchOnCTS();
  reader->SwitchOnPHOS();
  reader->SwitchOffEMCALCells();
  reader->SwitchOffPHOSCells();
 

  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  reader->SetCTSPtMin(0.2);
  
  reader->SetFidutialCut(fidCut);
  reader->Print("");
  
  
  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------
  
  //Detector Fidutial Cuts for analysis part
  AliFidutialCut * fidCut2 = new AliFidutialCut();
  fidCut2->DoCTSFidutialCut(kFALSE) ;
  fidCut2->DoEMCALFidutialCut(kTRUE) ;
  fidCut2->DoPHOSFidutialCut(kFALSE) ;
  
  fidCut2->SetSimpleCTSFidutialCut(0.9,0.,360.);
  fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  fidCut2->SetSimplePHOSFidutialCut(0.13,220.,320.);
  fidCut2->Print("");

  AliCaloPID * pid = new AliCaloPID();
  // use selection with simple weights
  pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
  pid->SetEMCALPhotonWeight(0.7);    pid->SetEMCALPi0Weight(0.7);
  
  pid->Print("");

  AliIsolationCut * ic = new AliIsolationCut();
  ic->SetConeSize(0.4);
  ic->SetPtThreshold(1.);
  ic->SetICMethod(AliIsolationCut::kPtThresIC);
  ic->Print("");

  //Photon Analysis
  AliAnaGammaDirect *ana = new AliAnaGammaDirect();
  ana->SetDebug(-1);
  ana->SetMinPt(5.);
  ana->SetCaloPID(pid);
  ana->SetFidutialCut(fidCut2);
  ana->SetIsolationCut(ic) ;
  ana->SetDetector("EMCAL");
  ana->SwitchOnIsolation();
  ana->SwitchOnCaloPID();
  ana->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  ana->SwitchOffFidutialCut();
  ana->SwitchOffInvariantMass();
  ana->SwitchOffDataMC() ;

  ana->Print("");

  //Photon - JETAN correlation
  AliAnaParticleJetFinderCorrelation *ana2 = new AliAnaParticleJetFinderCorrelation();
  ana2->SetDebug(-1);
  ana2->SetConeSize(1);  
  ana2->SetPtThresholdInCone(0.5);
  ana2->SetDeltaPhiCutRange(0.5,5.5);//Mostly Open Cuts 
  ana2->SetRatioCutRange(0.01,3); //Mostly Open Cuts
  ana2->UseJetRefTracks(kFALSE); //Not working now
  ana2->Print("");
  
  //---------------------------------------------------------------------
  // Set  analysis algorithm and reader
  //---------------------------------------------------------------------
  maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(ana,0);
  maker->AddAnalysis(ana2,1);
  maker->SetAODBranchName("PhotonJetCorrelation");
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  //maker->SwitchOffHistogramsMaker() ;  
  maker->SwitchOnAODsMaker()  ;
  //maker->SwitchOffAODsMaker() ; 
  
  maker->Print("");
  //
  printf("======================== \n");
  printf("END ConfigAnalysis() \n");
  printf("======================== \n");
  return maker ;
}

