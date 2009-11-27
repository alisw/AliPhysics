
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
	fidCut->DoEMCALFiducialCut(kTRUE) ;
	fidCut->DoPHOSFiducialCut(kTRUE) ;
	
	//fidCut->SetSimpleCTSFiducialCut(0.9,0.,360.);
	fidCut->SetSimpleEMCALFiducialCut(0.7,80.,190.);
	fidCut->SetSimplePHOSFiducialCut(0.12,220.,320.);
	
	fidCut->Print("");
	
	//-----------------------------------------------------------  
	// Reader
	//-----------------------------------------------------------
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);//10 for lots of messages
	
	//Switch on or off the detectors information that you want
	reader->SwitchOnPHOS();
	reader->SwitchOnEMCAL();
	reader->SwitchOffCTS();
	
	//Min particle pT
	reader->SetEMCALPtMin(0.2); 
	reader->SetPHOSPtMin(0.2);
	//reader->SetCTSPtMin(0.2);
	
	reader->SetFiducialCut(fidCut);
		
	//Remove the temporal AODs we create.	
	reader->SwitchOnCleanStdAOD();	
	
	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	
	AliCaloPID * pid = new AliCaloPID();
	// use selection with simple weights
	//pid->SetEMCALPhotonWeight(0.7);    pid->SetEMCALPi0Weight(0.7);
	pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7);
	pid->SetDispersionCut(2);
	pid->SetTOFCut(5.e-9);
	pid->SetDebug(-1);
	pid->Print("");

	AliAnaPhoton *anaphoton = new AliAnaPhoton();
        anaphoton->SetDebug(-1); //10 for lots of messages
        anaphoton->SetMinPt(0.2);
        anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
        anaphoton->SetCaloPID(pid);
//        anaphoton->SetFiducialCut(fidCut); //More acceptance selections if needed at this level
        anaphoton->SetCalorimeter("PHOS");
        anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
        anaphoton->SwitchOffCaloPID();
        anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
        anaphoton->SwitchOffFiducialCut();
        anaphoton->SetOutputAODName("PhotonsPHOS");
        anaphoton->SetOutputAODClassName("AliAODPWG4Particle");
        anaphoton->AddToHistogramsName("AnaPhotonPHOS_");
	 //Set Histrograms bins and ranges
        anaphoton->SetHistoPtRangeAndNBins(0, 50, 100) ;
        anaphoton->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
        anaphoton->SetHistoEtaRangeAndNBins(-0.2, 0.2, 100) ;
        anaphoton->Print("");

	fidCut->Print("");
	
	// >>> Second Analysis <<<
	
	AliNeutralMesonSelection *nms = new AliNeutralMesonSelection();
	nms->SetInvMassCutRange(0.1, 0.17)	;
	nms->KeepNeutralMesonSelectionHistos(kTRUE);
	//Set Histrograms bins and ranges
	nms->SetHistoERangeAndNBins(0, 50, 100) ;
	nms->SetHistoPtRangeAndNBins(0, 50, 100) ;
	nms->SetHistoAngleRangeAndNBins(0, 0.3, 100) ;
	nms->SetHistoIMRangeAndNBins(0, 0.4, 100) ;  

	AliAnaPi0EbE *anapi0 = new AliAnaPi0EbE();
	anapi0->SetDebug(-1);//10 for lots of messages
	anapi0->SetAnalysisType(AliAnaPi0EbE::kIMCalo);
	anapi0->SetInputAODName("PhotonsPHOS");
	anapi0->SetOutputAODName("Pi0sPHOS");
	anapi0->SetOutputAODClassName("AliAODPWG4Particle");
	anapi0->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anapi0->SetNeutralMesonSelection(nms);
	anapi0->AddToHistogramsName("AnaPi0EbEPHOS_");
	//Set Histrograms bins and ranges
	anapi0->SetHistoPtRangeAndNBins(0, 50, 100) ;
	anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
	anapi0->SetHistoEtaRangeAndNBins(-0.2, 0.2, 100) ;
	anapi0->Print("");
	
        AliAnaOmegaToPi0Gamma *anaomega = new AliAnaOmegaToPi0Gamma();
        anaomega->SetDebug(-1);//10 for lots of messages
        anaomega->SetInputAODName("Pi0sPHOS");
        anaomega->SetInputAODPhotonName("PhotonsPHOS");
	anaomega->SetNPID(3);
        anaomega->SetNEventsMixed(4);
        anaomega->SetNPtBinsMinMax(200,0,20.);
        anaomega->SetNMassBinsMinMas(200,0,1.);
        anaomega->SetPi0MassPeakWidthCut(0.015);
        anaomega->SwitchOnFiducialCut();
        anaomega->SwitchOnDataMC() ;//Access MC stack and fill more histograms
        anaomega->AddToHistogramsName("AnaNeuPHOS_");
        anaomega->Print("");


//for EMCAL
        AliCaloPID * pid1 = new AliCaloPID();
        // use selection with simple weights
        pid1->SetEMCALPhotonWeight(0.7);    pid1->SetEMCALPi0Weight(0.7);
        pid1->SetDispersionCut(2);
        pid1->SetTOFCut(5.e-9);
        pid1->SetDebug(-1);
        pid1->Print("");

        AliAnaPhoton *anaphoton1 = new AliAnaPhoton();
        anaphoton1->SetDebug(-1); //10 for lots of messages
        anaphoton1->SetMinPt(0.2);
        anaphoton1->SetMinDistanceToBadChannel(2, 4, 5);
        anaphoton1->SetCaloPID(pid);
        anaphoton1->SetCalorimeter("EMCAL");
        anaphoton1->SwitchOffDataMC() ;//Access MC stack and fill more histograms
        anaphoton1->SwitchOffCaloPID();
        anaphoton1->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
        anaphoton1->SwitchOffFiducialCut();
        anaphoton1->SetOutputAODName("PhotonsEMCAL");
        anaphoton1->SetOutputAODClassName("AliAODPWG4Particle");
        anaphoton1->AddToHistogramsName("AnaPhotonEMCAL_");
         //Set Histrograms bins and ranges
        anaphoton1->SetHistoPtRangeAndNBins(0, 50, 100) ;
        anaphoton1->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
        anaphoton1->SetHistoEtaRangeAndNBins(-0.8, 0.8, 100) ;
        anaphoton1->Print("");

	AliNeutralMesonSelection *nms1 = new AliNeutralMesonSelection();
        nms1->SetInvMassCutRange(0.1, 0.17)     ;
        nms1->KeepNeutralMesonSelectionHistos(kTRUE);
        //Set Histrograms bins and ranges
        nms1->SetHistoERangeAndNBins(0, 50, 100) ;
        nms1->SetHistoPtRangeAndNBins(0, 50, 100) ;
        nms1->SetHistoAngleRangeAndNBins(0, 0.3, 100) ;
        nms1->SetHistoIMRangeAndNBins(0, 0.4, 100) ;


        AliAnaPi0EbE *anapi01 = new AliAnaPi0EbE();
        anapi01->SetDebug(-1);//10 for lots of messages
        anapi01->SetAnalysisType(AliAnaPi0EbE::kIMCalo);
        anapi01->SetInputAODName("PhotonsEMCAL");
        anapi01->SetOutputAODName("Pi0sEMCAL");
        anapi01->SetOutputAODClassName("AliAODPWG4Particle");
        anapi01->SwitchOffDataMC() ;//Access MC stack and fill more histograms
        anapi01->SetNeutralMesonSelection(nms1);
        anapi01->AddToHistogramsName("AnaPi0EbEEMCAL_");
        //Set Histrograms bins and ranges
        anapi01->SetHistoPtRangeAndNBins(0, 50, 100) ;
        anapi01->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
        anapi01->SetHistoEtaRangeAndNBins(-0.8, 0.8, 100) ;
        anapi01->Print("");

        AliAnaOmegaToPi0Gamma *anaomega1 = new AliAnaOmegaToPi0Gamma();
        anaomega1->SetDebug(-1);//10 for lots of messages
        anaomega1->SetInputAODName("Pi0sEMCAL");
        anaomega1->SetInputAODPhotonName("PhotonsEMCAL");
        anaomega1->SetNPID(3);
        anaomega1->SetNEventsMixed(4);
        //anaomega1->SetNAsyBinsMinMax(200,0,1.);
        anaomega1->SetNPtBinsMinMax(200,0,20.);
        anaomega1->SetNMassBinsMinMas(200,0,1.);
        anaomega1->SetPi0MassPeakWidthCut(0.015);
        anaomega1->SwitchOnFiducialCut();
        anaomega1->SwitchOnDataMC() ;//Access MC stack and fill more histograms
        anaomega1->AddToHistogramsName("AnaNeuEMCAL_");
        anaomega1->Print("");


	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->SetAnaDebug(0);
	maker->AddAnalysis(anaphoton,0);
	maker->AddAnalysis(anapi0,1);
       	maker->AddAnalysis(anaomega,2);

        maker->AddAnalysis(anaphoton1,3);
        maker->AddAnalysis(anapi01,4);
       maker->AddAnalysis(anaomega1,5);

	maker->SwitchOnHistogramsMaker()  ;
	maker->SwitchOnAODsMaker()  ;
	
	maker->Print("");
	//
	printf("======================== \n");
	printf("END ConfigAnalysis() \n");
	printf("======================== \n");
	return maker ;
}
