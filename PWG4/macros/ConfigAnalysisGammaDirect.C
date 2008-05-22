/* $Id:$ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon analysis with ESDs
// Gamma in PHOS, 
// for EMCAL, PHOS by EMCAL where necessary
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
  //fidCut->DoEMCALFidutialCut(kFALSE) ;
  fidCut->DoPHOSFidutialCut(kTRUE) ;
  
  fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
  //fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  fidCut->SetSimplePHOSFidutialCut(0.13,220.,320.);
   
  fidCut->Print("");
  
  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  AliCaloTrackReader *reader = new AliCaloTrackESDReader();
  reader->SetDebug(-1);

  //Switch on or off the detectors information that you want
  reader->SwitchOffEMCAL();
  reader->SwitchOnCTS();
  reader->SwitchOnPHOS();
  
  //Min particle pT
  //reader->SetEMCALPtMin(0.5); 
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
  //fidCut2->DoEMCALFidutialCut(kTRUE) ;
  fidCut2->DoPHOSFidutialCut(kFALSE) ;
  
  fidCut2->SetSimpleCTSFidutialCut(0.9,0.,360.);
  //fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  fidCut2->SetSimplePHOSFidutialCut(0.13,220.,320.);
  fidCut2->Print("");

  AliCaloPID * pid = new AliCaloPID();
  // use selection with simple weights
  pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
  pid->SetEMCALPhotonWeight(0.7);    pid->SetEMCALPi0Weight(0.7);
  // use more complicated selection, particle weight depending on cluster energy
//   pid->UsePHOSPIDWeightFormula(kTRUE);
//   TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
//   TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
//   pid->SetPHOSPhotonWeightFormula(photonF);
//   pid->SetPHOSPi0WeightFormula(pi0F);
  pid->Print("");

  AliIsolationCut * ic = new AliIsolationCut();
  ic->SetConeSize(0.4);
  ic->SetPtThreshold(1.);
  ic->SetICMethod(AliIsolationCut::kPtThresIC);
  ic->Print("");

  AliAnaGammaDirect *ana = new AliAnaGammaDirect();
  ana->SetDebug(-1);
  ana->SetMinPt(5.);
  ana->SetCaloPID(pid);
  ana->SetFidutialCut(fidCut2);
  ana->SetIsolationCut(ic) ;
  ana->SetDetector("PHOS");
  ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms

  ana->SwitchOnCaloPID();
  ana->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
  ana->SwitchOnFidutialCut();
  //Select clusters with no pair, if both clusters with pi0 mass
  ana->SwitchOffInvariantMass();
  //Do isolation cut
  ana->SwitchOnIsolation();

  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  ana->SwitchOffReIsolation();

  //Multiple IC
  ana->SwitchOffSeveralIsolation() ;
//   ana->SwitchOnSeveralIsolation() ;
//   ana->SetNCones(2) ;     
//   ana->SetNPtThresFrac(2) ;     
//   ana->SetConeSizes(0, 0.3) ; ana->SetConeSizes(1, 0.4) ;    
//   ana->SetPtThresholds(0, 0.5) ;     ana->SetPtThresholds(1, 1.) ;   
//   ana->SetPtFractions(0, 1.) ;   ana->SetPtFractions(1, 1.5) ;  

  ana->Print("");

  //---------------------------------------------------------------------
  // Set  analysis algorithm and reader
  //---------------------------------------------------------------------
  maker = new AliAnaMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(ana,0);
  maker->SetAODBranchName("Photon");
  maker->SetAnaDebug(1)  ;
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
