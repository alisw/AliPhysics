/* $Id: $ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon - jet leading in cone correlation  analysis 
//
// Author : Gustavo Conesa Balbastre (INFN-LNF)
//------------------------------------

AliAnaMaker*  ConfigAnalysis()
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
  ana->SetDetector("PHOS");
  ana->SwitchOnIsolation();
  ana->SwitchOnCaloPID();
  ana->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
  ana->SwitchOffFidutialCut();
  ana->SwitchOffInvariantMass();
  ana->SwitchOffDataMC() ;

  ana->Print("");

  //Photon hadron correlation
  AliAnaParticleJetLeadingConeCorrelation *ana2 = new AliAnaParticleJetLeadingConeCorrelation();
  ana2->SetDebug(-1);
  ana2->SwitchOnCaloPID();
  ana2->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  ana2->SwitchOffFidutialCut();
  ana2->SwitchOffJetsOnlyInCTS();
  ana2->SwitchOffJetsRecalculation();
  //Analysis cuts for leading particle selection
  ana2->SetDeltaPhiCutRange(1.5,4.5); //Back-Leading particle angular cut
  ana2->SetLeadingRatioCutRange(0.,3);//Cut for the momentum of leading
  //Analysis cuts for jet selection
  ana2->SetppCollisions(); //Jet particles Pt threshold for different collisions
  ana2->SetCone(0.7); //Jet cone size
  ana2->SetJetPtThreshold(0.2); //Jet particle threshold 
  ana2->SetJetRatioCutRange(0.7, 1.3);//Only if SwitchOffJetsOnlyInCTS(); and SetJetSelectionMode(2)
  ana2->SetJetCTSRatioCutRange(0.3,1.3); //Only if SwitchOnJetsOnlyInCTS(); and SetJetSelectionMode(2)

  ana2->Print("");
  //---------------------------------------------------------------------
  // Set  analysis algorithm and reader
  //---------------------------------------------------------------------
  maker = new AliAnaMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(ana,0);
  maker->AddAnalysis(ana2,1);
  maker->SetAODBranchName("PhotonJetLCCorrelation");
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

