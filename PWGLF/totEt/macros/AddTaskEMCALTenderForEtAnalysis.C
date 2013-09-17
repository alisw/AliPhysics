// EMCal tender task adder
// Author: Jiri Kral

AliTender *AddTaskEMCALTender(const char *geoname="EMCAL_COMPLETEV1", AliEMCALRecParam *pars = 0 , Bool_t withNonlinearity = kTRUE , Bool_t withReclusterizing = kFALSE, Int_t trackmatchcuts==0)
{
  // Parameters: geoname = "EMCAL_FIRSTYEARV1" or "EMCAL_COMPLETEV1" or ""

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALTender", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  AliTender* ana = new  AliTender("TenderTask");
  
  mgr->AddTask(ana);

  // Adding EMCAL supply
  AliEMCALTenderSupply *EMCALSupply=new AliEMCALTenderSupply("EMCALtender");  

  EMCALSupply->SetEMCALGeometryName(geoname);  

  // prepare the reco params ------------------------------------------------
	if( pars == 0 ){
		// you can write your reco params here to avoid loading them automatically
		// from OCDB during execution time
		AliEMCALRecParam *params = new AliEMCALRecParam();
		// reclustering parameters
		// use v1 for pp and v2 for PbPb
		params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2);
		params->SetClusteringThreshold(0.1); // 100 MeV
		params->SetMinECut(0.05); //50 MeV  
		params->SetW0(4.5);
		// you may want to enable the timing cut
		params->SetTimeCut(1e6);// s
		params->SetTimeMin(-1);
		params->SetTimeMax(1e6);//s

		EMCALSupply->SetRecParam(params);
	}
	else{
		cout << "------- TENDER is using supplied reco params -------" << endl;
		pars->Print( "reco" );
		cout << "----------------------------------------------------" << endl;
		EMCALSupply->SetRecParam(pars);
	}

  // prepare tender parameters ----------------------------------------------
  EMCALSupply->SetDebugLevel( 0 );

  // fiducial cut
  EMCALSupply->SetNumberOfCellsFromEMCALBorder( 1 );

  // nonlinearity
  EMCALSupply->SetNonLinearityFunction( AliEMCALTenderSupply::kBeamTestCorrected );

  // track matching parameters
  //EMCALSupply->SetMass(0.139);
  //EMCALSupply->SetStep(5);
  //EMCALSupply->SwitchOnCutEtaPhiSum(); 
  //EMCALSupply->SetRCut(0.025);
  EMCALSupply->SwitchOnCutEtaPhiSeparate();
//   EMCALSupply->SetEtaCut(0.015);
//   EMCALSupply->SetPhiCut(0.03);
  if(trackmatchcuts==0){//default
    EMCALSupply->SetEtaCut(0.025);
    EMCALSupply->SetPhiCut(0.05);
  }
  if(trackmatchcuts==1){//tighter
    EMCALSupply->SetEtaCut(0.015);
    EMCALSupply->SetPhiCut(0.03);
  }
  if(trackmatchcuts==2){//looser
    EMCALSupply->SetEtaCut(0.035);
    EMCALSupply->SetPhiCut(0.07);
  }

  // switches ---------------------------------------------------------------
  EMCALSupply->SwitchOnBadCellRemove();
  EMCALSupply->SwitchOnExoticCellRemove();
  EMCALSupply->SwitchOnCalibrateEnergy();
  EMCALSupply->SwitchOnCalibrateTime();
  EMCALSupply->SwitchOnUpdateCell();
  if(withReclusterizing) EMCALSupply->SwitchOnReclustering();
  EMCALSupply->SwitchOnClusterBadChannelCheck();
  EMCALSupply->SwitchOnClusterExoticChannelCheck();
  EMCALSupply->SwitchOnCellFiducialRegion();
  EMCALSupply->SwitchOnReCalibrateCluster();
  EMCALSupply->SwitchOnRecalculateClusPos();
  EMCALSupply->SwitchOnRecalShowerShape();
  EMCALSupply->SwitchOnRecalDistBadChannel();
 if(withNonlinearity) EMCALSupply->SwitchOnNonLinearityCorrection();
else{cout<<"WARNING:  TURNING OFF NONLINEARITY"<<endl;}
  EMCALSupply->SwitchOnTrackMatch();
  

  ana->AddSupply(EMCALSupply);
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
//  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosEmcalTender", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("emcal_tender_event", AliESDEvent::Class(),
                           AliAnalysisManager::kExchangeContainer,"emcal_tender");
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
 
