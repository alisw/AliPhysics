/* $Id:  $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do example analysis with MC
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
  fidCut->DoPHOSFidutialCut(kTRUE) ;
  
  //fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
  //fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  fidCut->SetSimplePHOSFidutialCut(0.13,220.,320.);
  
  //   //Fidutial cut EMCAL, 5 regions
  //   Float_t etamax[]={0.67,0.51,0.16,-0.21,-0.61};
  //   TArrayF etamaxarr(5,etamax);
  //   fidCut->AddEMCALFidCutMaxEtaArray(etamaxarr);
  //   Float_t etamin[]={0.61,0.21,-0.16,-0.51,-0.67};
  //   TArrayF etaminarr(5,etamin);
  //   fidCut->AddEMCALFidCutMinEtaArray(etaminarr);
  //   Float_t phimax[]={99.*TMath::DegToRad(), 119.*TMath::DegToRad(), 139.*TMath::DegToRad(), 
  // 		    159.*TMath::DegToRad(), 179.*TMath::DegToRad(), 189.*TMath::DegToRad()};
  //   TArrayF phimaxarr(6,phimax);
  //   fidCut->AddEMCALFidCutMaxPhiArray(phimaxarr);
  //   Float_t phimin[]={81.*TMath::DegToRad(), 101.*TMath::DegToRad(), 121.*TMath::DegToRad(), 
  // 		    141.*TMath::DegToRad(), 161.*TMath::DegToRad(), 181.*TMath::DegToRad()};
  //   TArrayF phiminarr(6,phimin);
  //   fidCut->AddEMCALFidCutMinPhiArray(phiminarr);
  
  fidCut->Print("");
  
  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  AliCaloTrackMCReader *reader = new AliCaloTrackMCReader();
  reader->SetDebug(-1);
  
 //Switch on or off the detectors information that you want
  reader->SwitchOffEMCAL();
  reader->SwitchOffCTS();
  reader->SwitchOnPHOS();
  reader->SwitchOffEMCALCells();
  reader->SwitchOffPHOSCells();

  //Min particle pT
  reader->SetEMCALPtMin(0.); 
  reader->SetPHOSPtMin(0.);
  reader->SetCTSPtMin(.2);
 
  reader->SetFidutialCut(fidCut);
  reader->Print("");
  
  
  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------
  
  //Detector Fidutial Cuts for analysis part
  AliFidutialCut * fidCut2 = new AliFidutialCut();
  fidCut2->DoCTSFidutialCut(kFALSE) ;
  fidCut2->DoEMCALFidutialCut(kFALSE) ;
  fidCut2->DoPHOSFidutialCut(kTRUE) ;
  fidCut2->SetSimplePHOSFidutialCut(0.1,240.,280.);

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


  AliAnaExample *ana = new AliAnaExample();
  ana->SetDebug(-1);
  ana->SetCaloPID(pid);
  ana->SetFidutialCut(fidCut2);
  ana->SetDetector("PHOS");
  //ana->SwitchOnDataMC();
  ana->SetMinPt(0.);
  ana->SetOutputAODName("Example");
  ana->SetOutputAODClassName("AliAODPWG4Particle"); //Or AliAODPWG4ParticleCorrelation
  //Set Histrograms bins and ranges
//	ana->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;  
  
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
