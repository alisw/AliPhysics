/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon analysis with ESDs
// Gamma in PHOS, 
// for EMCAL, PHOS by EMCAL where necessary
// First find photons with AliAnaPhoton, then 
// isolate them with AliAnaParticleIsolation
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
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);
	
	//Switch on or off the detectors information that you want
	reader->SwitchOffEMCAL();
	reader->SwitchOnCTS();
	reader->SwitchOnPHOS();
	
	//Min particle pT
	reader->SetEMCALPtMin(0.2); 
	reader->SetPHOSPtMin(0.2);
	reader->SetCTSPtMin(0.2);
	
	reader->SetFiducialCut(fidCut);
	
	//     //We want tracks fitted in the detectors:
	//     ULong_t status=AliAODTrack::kTPCrefit;
	//     status|=AliAODTrack::kITSrefit; //(default settings)
	
	//     We want tracks whose PID bit is set:
	//     ULong_t status =AliAODTrack::kITSpid;
	//     status|=AliAODTrack::kTPCpid;	
	
	//	reader->SetTrackStatus(status);
		
	//Remove the temporal AODs we create.	
	reader->SwitchOnCleanStdAOD();	
	
	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	//>>>> First Analysis <<<<
	//Detector Fiducial Cuts for analysis part
	AliFiducialCut * fidCut2 = new AliFiducialCut();
	fidCut2->DoCTSFiducialCut(kFALSE) ;
	fidCut2->DoEMCALFiducialCut(kFALSE) ;
	fidCut2->DoPHOSFiducialCut(kFALSE) ;
	
	//fidCut2->SetSimpleCTSFiducialCut(0.9,0.,360.);
	//fidCut2->SetSimpleEMCALFiducialCut(0.7,80.,190.);
	fidCut2->SetSimplePHOSFiducialCut(0.12,220.,320.);
	
//	//Select particles in N regions of the detectors
//	Float_t etamax[]={0.67,0.51,0.16,-0.21,-0.61};
//	TArrayF etamaxarr(5,etamax);
//	fidCut2->AddEMCALFidCutMaxEtaArray(etamaxarr);
//	Float_t etamin[]={0.61,0.21,-0.16,-0.51,-0.67};
//	TArrayF etaminarr(5,etamin);
//	fidCut2->AddEMCALFidCutMinEtaArray(etaminarr);
//	Float_t phimax[]={99., 119., 139., 159., 179., 189.};
//	TArrayF phimaxarr(6,phimax);
//	fidCut2->AddEMCALFidCutMaxPhiArray(phimaxarr);
//	Float_t phimin[]={81., 101.  , 121. , 141. , 161. , 181. };
//	TArrayF phiminarr(6,phimin);
//	fidCut2->AddEMCALFidCutMinPhiArray(phiminarr);
//	fidCut2->Print("");
//	
	AliCaloPID * pid = new AliCaloPID();
	// use selection with simple weights
	pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
	pid->SetEMCALPhotonWeight(0.7);    pid->SetEMCALPi0Weight(0.7);
	// use more complicated selection, particle weight depending on cluster energy
	// pid->UsePHOSPIDWeightFormula(kTRUE);
	// TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
	// TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
	// pid->SetPHOSPhotonWeightFormula(photonF);
	// pid->SetPHOSPi0WeightFormula(pi0F);
	pid->Print("");
	
	
	AliAnaPhoton *anaphoton = new AliAnaPhoton();
	anaphoton->SetDebug(-1); //10 for lots of messages
	anaphoton->SetMinPt(0.2);
	anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
	anaphoton->SetCaloPID(pid);
	anaphoton->SetFiducialCut(fidCut2); //More acceptance selections if needed at this level
	anaphoton->SetCalorimeter("PHOS");
	anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anaphoton->SwitchOnCaloPID();
	anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
	anaphoton->SwitchOnFiducialCut();
	anaphoton->SetOutputAODName("Photons");
	anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
	//Set Histrograms bins and ranges
	//	anaphoton->SetHistoPtRangeAndNBins(0, 50, 100) ;
	//	anaphoton->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
	//	anaphoton->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
	anaphoton->Print("");
	
	// >>>> Second Analysis <<<< Isolate the photons
	AliIsolationCut * ic = new AliIsolationCut();
	ic->SetConeSize(0.5);
	ic->SetPtThreshold(1.);
	ic->SetICMethod(AliIsolationCut::kPtThresIC);
	ic->Print("");
	
	AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
	anaisol->SetDebug(-1);
	anaisol->SetMinPt(5);
	anaisol->SetInputAODName("Photons");
	anaisol->SetCalorimeter("PHOS");
	anaisol->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	//Select clusters with no pair, if both clusters with pi0 mass
	anaisol->SwitchOffInvariantMass();
	//anaisol->SetNeutralMesonSelection(nms);
	//Do isolation cut
	anaisol->SetIsolationCut(ic);	
	//Do or not do isolation with previously produced AODs.
	//No effect if use of SwitchOnSeveralIsolation()
	anaisol->SwitchOffReIsolation();
	//Multiple IC
	anaisol->SwitchOffSeveralIsolation() ;
	
	anaisol->Print("");
	
		
	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->AddAnalysis(anaphoton,0);
	maker->AddAnalysis(anaisol,1);
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
