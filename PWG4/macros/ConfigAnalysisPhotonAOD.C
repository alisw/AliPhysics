/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do photon identification analysis with AODs
// Gamma in EMCAL. For PHOS clusters change 
// EMCAL by PHOS where necessary
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
  
  
  //Detector Fiducial Cuts
  AliFiducialCut * fidCut = new AliFiducialCut();
  fidCut->DoCTSFiducialCut(kFALSE) ;
  fidCut->DoEMCALFiducialCut(kFALSE) ;
  fidCut->DoPHOSFiducialCut(kFALSE) ;
  
  //fidCut->SetSimpleCTSFiducialCut(0.9,0.,360.);
  //fidCut->SetSimpleEMCALFiducialCut(0.7,80.,190.);
  //fidCut->SetSimplePHOSFiducialCut(0.13,220.,320.);
   
  fidCut->Print("");
  
  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  AliCaloTrackAODReader *reader = new AliCaloTrackAODReader();
  reader->SetDebug(-1);

  //Switch on or off the detectors information that you want
  reader->SwitchOnEMCAL();
  reader->SwitchOffCTS();
  reader->SwitchOffPHOS();
  
  //Kine
  reader->SwitchOffStack();          // On  by default, remember to SwitchOnMCData() in analysis classes
  reader->SwitchOnAODMCParticles();  // Off by default, remember to SwitchOnMCData() in analysis classes

  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  //reader->SetPHOSPtMin(0.5);
  //reader->SetCTSPtMin(0.2);
  
  //In case of generating jet events (with PYTHIA), if pt hard bin is small
  //reject events with large difference between ptHard and triggered jet	
  //reader->SetPtHardAndJetPtComparison(kTRUE);
	
  reader->SetFiducialCut(fidCut);

  //Embedding/mixing/merging analysis with other events in another AOD file
  reader->SetSecondInputFileName("./aod.root");
  //Standard event loop can have less events, start mixing at event N
  reader->SetSecondInputFirstEvent(0);   
	
  // Analysis with tracks, select only tracks with
  // following bits
	
  //     //We want tracks fitted in the detectors:
  //     ULong_t status=AliAODTrack::kTPCrefit;
  //     status|=AliAODTrack::kITSrefit;
    
  //     We want tracks whose PID bit is set:
  //     ULong_t status =AliAODTrack::kITSpid;
  //     status|=AliAODTrack::kTPCpid;	

  //	reader->SetTrackStatus(status);
	
  reader->Print("");
  
  
  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------
  
  //Detector Fiducial Cuts for analysis part
  AliFiducialCut * fidCut2 = new AliFiducialCut();
  fidCut2->DoCTSFiducialCut(kFALSE) ;
  fidCut2->DoEMCALFiducialCut(kFALSE) ;
  fidCut2->DoPHOSFiducialCut(kFALSE) ;
  
  //fidCut2->SetSimpleCTSFiducialCut(0.9,0.,360.);
  //fidCut2->SetSimpleEMCALFiducialCut(0.7,80.,190.);
  //fidCut2->SetSimplePHOSFiducialCut(0.13,220.,320.);

  fidCut2->Print("");

  AliCaloPID * pid = new AliCaloPID();
  // use selection with simple weights
  //pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
  // use more complicated selection, particle weight depending on cluster energy
//   pid->UsePHOSPIDWeightFormula(kTRUE);
//   TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
//   TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
//   pid->SetPHOSPhotonWeightFormula(photonF);
//   pid->SetPHOSPi0WeightFormula(pi0F);
//Check these cuts for EMCAL
  pid->SetDispersionCut(1.5);
  pid->SetTOFCut(5.e-9);
  pid->SetDebug(-1);
  pid->Print("");

  AliAnaPhoton *ana = new AliAnaPhoton();
  ana->SetDebug(-1);
  ana->SetMinPt(5.);
  ana->SetMinDistanceToBadChannel(2, 4, 5);
  ana->SetCaloPID(pid);
  ana->SetFiducialCut(fidCut2);
  ana->SetCalorimeter("EMCAL");
  ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  ana->SwitchOffCaloPID();
  ana->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
  ana->SwitchOnTrackMatchRejection(); //Only in use when OnCaloPID
  ana->SwitchOffFiducialCut();
  ana->SetOutputAODName("Photons");
  ana->SetOutputAODClassName("AliAODPWG4Particle");
  //Set Histrograms bins and ranges
//	ana->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  ana->Print("");

  //---------------------------------------------------------------------
  // Set  analysis algorithm and reader
  //---------------------------------------------------------------------
  maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(ana,0);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  
  maker->Print("");
  //
  printf("======================== \n");
  printf("END ConfigAnalysis() \n");
  printf("======================== \n");
  return maker ;
}
