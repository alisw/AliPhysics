/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do photon identification analysis with ESDs
// Gamma in PHOS. For EMCAL clusters change 
// PHOS by EMCAL where necessary
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
  fidCut->DoCTSFidutialCut(kFALSE) ;
  fidCut->DoEMCALFidutialCut(kFALSE) ;
  fidCut->DoPHOSFidutialCut(kFALSE) ;
  
  //fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
  //fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  //fidCut->SetSimplePHOSFidutialCut(0.13,220.,320.);
   
  fidCut->Print("");
  
  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
  reader->SetDebug(-1);

  //Switch on or off the detectors information that you want
  reader->SwitchOffEMCAL();
  reader->SwitchOffCTS();
  reader->SwitchOnPHOS();
  
  //Min particle pT
  //reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  //reader->SetCTSPtMin(0.2);
  
  reader->SetFidutialCut(fidCut);
  reader->Print("");
  
  
  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------
  
  //Detector Fidutial Cuts for analysis part
  AliFidutialCut * fidCut2 = new AliFidutialCut();
  fidCut2->DoCTSFidutialCut(kFALSE) ;
  fidCut2->DoEMCALFidutialCut(kFALSE) ;
  fidCut2->DoPHOSFidutialCut(kFALSE) ;
  
  //fidCut2->SetSimpleCTSFidutialCut(0.9,0.,360.);
  //fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  //fidCut2->SetSimplePHOSFidutialCut(0.13,220.,320.);

  fidCut2->Print("");

  AliCaloPID * pid = new AliCaloPID();
  // use selection with simple weights
  pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
  // use more complicated selection, particle weight depending on cluster energy
//   pid->UsePHOSPIDWeightFormula(kTRUE);
//   TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
//   TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
//   pid->SetPHOSPhotonWeightFormula(photonF);
//   pid->SetPHOSPi0WeightFormula(pi0F);

  pid->SetDispersionCut(1.5);
  pid->SetTOFCut(5.e-9);
  pid->SetDebug(-1);
  pid->Print("");

  AliAnaPhoton *ana = new AliAnaPhoton();
  ana->SetDebug(-1);
  ana->SetMinPt(5.);
  ana->SetMinDistanceToBadChannel(2, 4, 5);
  ana->SetCaloPID(pid);
  ana->SetFidutialCut(fidCut2);
  ana->SetCalorimeter("PHOS");
  ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  ana->SwitchOffCaloPID();
  ana->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
  ana->SwitchOnTrackMatchRejection(); //Only in use when OnCaloPID
  ana->SwitchOffFidutialCut();
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
