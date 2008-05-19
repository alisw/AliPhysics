/* $Id:  $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon analysis with ESDs
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
  //Select particles depending on acceptance
  //kFALSE, open cuts
  fidCut->DoCTSFidutialCut(kTRUE) ;
  fidCut->DoEMCALFidutialCut(kTRUE) ;
  fidCut->DoPHOSFidutialCut(kTRUE) ;
  
  //Select particles in one region of the detectors
  fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
  fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  fidCut->SetSimplePHOSFidutialCut(0.13,220.,320.);
  
  //Select particles in N regions of the detectors
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
  
  //   //Fidutial cut PHOS
  //   TArrayF etamaxarr(1,etamax); 
  //   etamaxarr.SetAt(0.12,0);
  //   fidCut->AddPHOSFidCutMaxEtaArray(etamaxarr);
  //   TArrayF etaminarr(1,etamin); 
  //   etaminarr.SetAt(-0.12,0);
  //   fidCut->AddPHOSFidCutMinEtaArray(etaminarr);
  //   TArrayF phimaxarr(1,phimax);
  //   phimaxarr.SetAt(320.*TMath::DegToRad(),0);   
  //   fidCut->AddPHOSFidCutMaxPhiArray(phimaxarr); 
  //   TArrayF phiminarr(1,phimin);
  //   phiminarr.SetAt(220.*TMath::DegToRad(),0);   
  //   fidCut->AddPHOSFidCutMinPhiArray(phiminarr);
  
  //   //Fidutial cut CTS
  //   TArrayF etamaxarr(1,etamax); 
  //   etamaxarr.SetAt(0.12,0);
  //   fidCut->AddCTSFidCutMaxEtaArray(etamaxarr);
  //   TArrayF etaminarr(1,etamin); 
  //   etaminarr.SetAt(-0.12,0);
  //   fidCut->AddCTSFidCutMinEtaArray(etaminarr);
  //   TArrayF phimaxarr(1,phimax);
  //   phimaxarr.SetAt(320.*TMath::DegToRad(),0);   
  //   fidCut->AddCTSFidCutMaxPhiArray(phimaxarr); 
  //   TArrayF phiminarr(1,phimin);
  //   phiminarr.SetAt(220.*TMath::DegToRad(),0);   
  //   fidCut->AddCTSFidCutMinPhiArray(phiminarr);
  
  fidCut->Print("");
  
  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  AliCaloTrackReader *reader = new AliCaloTrackESDReader();
  reader->SetDebug(-1);
  
  //Switch on or off the detectors information that you want
  //It will fill the corresponding detector arrays
  reader->SwitchOnEMCAL();
  reader->SwitchOnCTS();
  reader->SwitchOnPHOS();
  reader->SwitchOnEMCALCells();
  reader->SwitchOnPHOSCells();
  //reader->SwitchOffEMCAL();
  //reader->SwitchOffCTS();
  //reader->SwitchOffPHOS();
  //reader->SwitchOffEMCALCells();
  //reader->SwitchOffPHOSCells();

  //Min particle pT
  //Selections done while filling the detector arrays
  reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  reader->SetCTSPtMin(0.2);
 
  reader->SetFidutialCut(fidCut);

  reader->Print("");
  
  
  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------
  
  //Some settings for the analysis

  //Detector Fidutial Cuts for analysis part
  //You can define different cuts for your analysis
  AliFidutialCut * fidCut2 = new AliFidutialCut();
  fidCut2->DoCTSFidutialCut(kFALSE) ;
  fidCut2->DoEMCALFidutialCut(kTRUE) ;
  fidCut2->DoPHOSFidutialCut(kTRUE) ;
  
  fidCut2->SetSimpleCTSFidutialCut(0.9,0.,360.);
  fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,160.);
  fidCut2->SetSimplePHOSFidutialCut(0.13,220.,280.);

  AliCaloPID * pid = new AliCaloPID();
  // use selection with simple weights
  pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
  pid->SetEMCALPhotonWeight(0.7);    pid->SetEMCALPi0Weight(0.7);
  // use more complicated selection, particle weight depending on cluster energy,
  // only for PHOS
  //   pid->UsePHOSPIDWeightFormula(kTRUE);
  //   TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
  //   TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
  //   pid->SetPHOSPhotonWeightFormula(photonF);
  //   pid->SetPHOSPi0WeightFormula(pi0F);
  
  AliAnaExample *ana = new AliAnaExample();
  ana->SetDebug(-1);
  //Do PID selection with settings defined up
  ana->SetCaloPID(pid);
  ana->SwitchOnCaloPID();
  //Do Acceptance selection with settings defined up
  ana->SetFidutialCut(fidCut2);

  ana->SetPdg(AliCaloPID::kPhoton); //plot identified photons
  ana->SetDetector("PHOS"); //Detector for the analysis
  ana->SetMinPt(2);// Minimum pt of clusters/tracks
  
  ana->SwitchOnDataMC() ; //Access to the stack and fill 
  // some histograms with MC information

  ana->Print("");

  //AliAnaExample *ana2 = new AliAnaExample();
  
  //---------------------------------------------------------------------
  // Set  analysis algorithm and reader
  //---------------------------------------------------------------------
  maker = new AliAnaMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(ana,0);
  //maker->AddAnalysis(ana2,1);
  maker->SetAODBranchName("Test");
  maker->SetAnaDebug(1)  ;
  maker->SwitchOnAODsMaker()  ;
  //maker->SwitchOffAODsMaker() ; 
  maker->SwitchOnHistogramsMaker()  ;
  //maker->SwitchOffHistogramsMaker() ;  

  
  maker->Print("");
  //
  printf("======================== \n");
  printf("END ConfigAnalysis() \n");
  printf("======================== \n");
  return maker ;
}
