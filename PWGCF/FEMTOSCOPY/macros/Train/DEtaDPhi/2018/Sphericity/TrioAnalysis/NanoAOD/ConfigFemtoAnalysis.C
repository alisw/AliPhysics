enum ESys  { kPIpPIpPIp, kPImPImPIm, kPIpPImPIm, kKpKpKp, kKmKmKm, kKpKmKm, kPPP, kAPAPAP, kPAPAP, kPIpKpKp, kPIpKmKm, kPPIpPIp, kPPImPIm, kPKpKp, kPKmKm, kAPKpKp, kAPKmKm, kPIpKpP, nSys };

const char *sysNames[nSys]      = {"PIpPIpPIp", "PImPImPIm", "PIpPImPIm", "KpKpKp", "KmKmKm", "KpKmKm", "PPP","APAPAP","PAPAP", "PIpKpKp", "PIpKmKm", "PPIpPIp", "PPImPIm", "PKpKp","PKmKm","APKpKp","APKmKm", "PIpKpP"};
const bool runSys[nSys]         = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

const int nMultBins = 1;
const int multBins[nMultBins+1] = {3, 2000};
const int runMult[nMultBins]    = {1};
const int vertZMin = -10;
const int vertZMax = 10;

const  double PionMass = 0.13956995;
const  double KaonMass = 0.493677;
const  double ProtonMass = 0.938272013;
const  double LambdaMass = 1.115683;
const  double XiMass = 1.32171;


AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReaderPP(bool mcAnalysis);
AliFemtoTrioAnalysis* GetAnalysis(bool doEventMixing);
AliFemtoSphericityEventCut* GetEventCut(int multMin=3,int multMax=100000,int zVertMin=-8.0,int zVertMax=8.0, double stMin=0, double stMax=1);
AliFemtoESDTrackCut* GetTrackCut(AliFemtoTrio::EPart particle);
void GetParticlesForSystem(ESys system, AliFemtoTrio::EPart &firstParticle, AliFemtoTrio::EPart &secondParticle, AliFemtoTrio::EPart &thirdParticle);

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, bool sepCuts=false, int year=2015, bool ppAnalysis=false, bool eventMixing=false, bool doMonitors=true)
{

  

  // create analysis managers
  AliFemtoManager* Manager = new AliFemtoManager();
  AliFemtoModelManager *modelMgr = new AliFemtoModelManager();
  
  // add event reader
  if(ppAnalysis){
    AliFemtoEventReaderAODChain* ReaderPP = GetReaderPP(mcAnalysis);
    Manager->SetEventReader(ReaderPP);
  }
  else if(year==2015){
    AliFemtoEventReaderAODMultSelection* Reader2015 = GetReader2015(mcAnalysis);
    Manager->SetEventReader(Reader2015);
  }
  else if(year==2011){
    AliFemtoEventReaderAODChain* Reader2011 = GetReader2011(mcAnalysis);
    Manager->SetEventReader(Reader2011);
  }
  
  // declare necessary objects
  AliFemtoTrioAnalysis         *trioAnalysis[500];
  AliFemtoSphericityEventCut     *eventCut[500];

  AliFemtoCutMonitorEventMult   *monitorEventMultPass[500];
  AliFemtoCutMonitorEventMult   *monitorEventMultFail[500];
  AliFemtoCutMonitorCollections  *monitorCollPass[500];
  AliFemtoCutMonitorCollections  *monitorCollFail[500];
  
  AliFemtoCutMonitorParticlePID   *monitorPIDPass[500];
  AliFemtoCutMonitorParticlePID   *monitorPIDFail[500];
  AliFemtoCutMonitorParticleYPt   *monitorYPtPass[500];
  AliFemtoCutMonitorParticleYPt   *monitorYPtFail[500];
  	
  AliFemtoTrioDEtaDPhiFctn         *detadphifctn[5000];
  
  
  // setup analysis
  int anIter = 0;
  
  AliFemtoTrio::EPart firstParticle, secondParticle, thirdParticle;
  
  for (int iSys=0; iSys<nSys; iSys++)
    {
      if (!runSys[iSys]) continue;
    
      // get particle cuts
      GetParticlesForSystem((ESys)iSys,firstParticle,secondParticle, thirdParticle);
      AliFemtoMJTrackCut *firstTrackCut  = GetTrackCut(firstParticle);
      AliFemtoMJTrackCut *secondTrackCut = GetTrackCut(secondParticle);
      AliFemtoMJTrackCut *thirdTrackCut  = GetTrackCut(thirdParticle);

 
      for(int iSphericity=0;iSphericity<3;iSphericity++)
	{
    
	  for(int iMult=0;iMult<nMultBins;iMult++)
	    {
	      if(!runMult[iMult]) continue;
     
	      // create new analysis
	      trioAnalysis[anIter] = GetAnalysis(eventMixing);
        
	      // get event cut
	      if(iSphericity==0)
		eventCut[anIter] = GetEventCut(multBins[iMult],multBins[iMult+1],vertZMin,vertZMax);
	      if(iSphericity==1)
		eventCut[anIter] = GetEventCut(multBins[iMult],multBins[iMult+1],vertZMin,vertZMax,0,0.3);
	      if(iSphericity==2)
		eventCut[anIter] = GetEventCut(multBins[iMult],multBins[iMult+1],vertZMin,vertZMax,0.7,1);
	      
	      // event monitors
	      if(doMonitors)
		{
		  monitorEventMultPass[anIter] = new AliFemtoCutMonitorEventMult(Form("monitorEventMultPass%sSt%iM%i", sysNames[iSys], iSphericity, iMult));
		  monitorEventMultFail[anIter] = new AliFemtoCutMonitorEventMult(Form("monitorEventMultFail%sSt%iM%i", sysNames[iSys], iSphericity, iMult));
		  eventCut[anIter]->AddCutMonitor(monitorEventMultPass[anIter], monitorEventMultFail[anIter]);
	      
		  monitorCollPass[anIter] = new AliFemtoCutMonitorCollections(Form("monitorCollPass%sSt%iM%i", sysNames[iSys],iSphericity, iMult));
		  monitorCollFail[anIter] = new AliFemtoCutMonitorCollections(Form("monitorCollFail%sSt%iM%i", sysNames[iSys],iSphericity, iMult));
		  eventCut[anIter]->AddCutMonitor(monitorCollPass[anIter], monitorCollFail[anIter]);
	    }
	  
	  // setup anallysis cuts
	  trioAnalysis[anIter]->SetEventCut(eventCut[anIter]);
	  trioAnalysis[anIter]->SetV0SharedDaughterCut(true);



	  if(iSys == kPIpPIpPIp || iSys == kPImPImPIm || iSys == kKpKpKp || iSys == kKmKmKm || iSys == kPPP || iSys == kAPAPAP)
	    {
	      trioAnalysis[anIter]->SetFirstParticleCut(firstTrackCut);
	  
	      trioAnalysis[anIter]->SetCollection1type(firstParticle);
	      
	      trioAnalysis[anIter]->SetMinSizePart1Collection(3);
	    }
	  else if(iSys == kPIpPImPIm || iSys == kKpKmKm || iSys == kPAPAP || iSys == kPKpKp || iSys == kPKmKm || iSys == kAPKpKp || iSys == kAPKmKm || iSys == kPIpKpKp || iSys == kPIpKmKm || iSys == kPPIpPIp || iSys == kPPImPIm)
	    {
	      trioAnalysis[anIter]->SetFirstParticleCut(firstTrackCut);
	      trioAnalysis[anIter]->SetSecondParticleCut(secondTrackCut);
	  
	      trioAnalysis[anIter]->SetCollection1type(firstParticle);
	      trioAnalysis[anIter]->SetCollection2type(secondParticle);

	      trioAnalysis[anIter]->SetMinSizePart1Collection(1);
	      trioAnalysis[anIter]->SetMinSizePart2Collection(2);
	    }
	  else if(kPIpKpP)
	    {
	      trioAnalysis[anIter]->SetFirstParticleCut(firstTrackCut);
	      trioAnalysis[anIter]->SetSecondParticleCut(secondTrackCut);
	      trioAnalysis[anIter]->SetThirdParticleCut(thirdTrackCut);

	      trioAnalysis[anIter]->SetCollection1type(firstParticle);
	      trioAnalysis[anIter]->SetCollection2type(secondParticle);
	      trioAnalysis[anIter]->SetCollection3type(thirdParticle);

	      trioAnalysis[anIter]->SetMinSizePart1Collection(1);
	      trioAnalysis[anIter]->SetMinSizePart2Collection(1);
	      trioAnalysis[anIter]->SetMinSizePart3Collection(1);
	      //collection min size by default is 1
	    }



	  if(doMonitors)
	    {
	      if(iSys == kPIpPIpPIp || iSys == kPImPImPIm) //PI
		{
		  monitorPIDPass[anIter] = new AliFemtoCutMonitorParticlePID(Form("monitorPIDPass%sSt%iM%i", sysNames[iSys], iSphericity, iMult),0);//0-pion,1-kaon,2-proton
		  monitorPIDFail[anIter] = new AliFemtoCutMonitorParticlePID(Form("monitorPIDFail%sSt%iM%i", sysNames[iSys], iSphericity, iMult),0);//0-pion,1-kaon,2-proton
		  firstTrackCut->AddCutMonitor(monitorPIDPass[anIter],monitorPIDFail[anIter]);

		  monitorYPtPass[anIter] = new AliFemtoCutMonitorParticleYPt(Form("monitorYPtPass%sSt%iM%i", sysNames[iSys], iSphericity, iMult),PionMass);
		  monitorYPtFail[anIter] = new AliFemtoCutMonitorParticleYPt(Form("monitorYPtFail%sSt%iM%i", sysNames[iSys], iSphericity, iMult),PionMass);
		  firstTrackCut->AddCutMonitor(monitorYPtPass[anIter],monitorYPtFail[anIter]);
		}
	      if(iSys == kKpKpKp || iSys == kKmKmKm) //K
		{
		  monitorPIDPass[anIter] = new AliFemtoCutMonitorParticlePID(Form("monitorPIDPass%sSt%iM%i", sysNames[iSys], iSphericity, iMult),1);//0-pion,1-kaon,2-proton
		  monitorPIDFail[anIter] = new AliFemtoCutMonitorParticlePID(Form("monitorPIDFail%sSt%iM%i", sysNames[iSys], iSphericity, iMult),1);//0-pion,1-kaon,2-proton
		  firstTrackCut->AddCutMonitor(monitorPIDPass[anIter],monitorPIDFail[anIter]);

		  monitorYPtPass[anIter] = new AliFemtoCutMonitorParticleYPt(Form("monitorYPtPass%sSt%iM%i", sysNames[iSys], iSphericity, iMult),KaonMass);
		  monitorYPtFail[anIter] = new AliFemtoCutMonitorParticleYPt(Form("monitorYPtFail%sSt%iM%i", sysNames[iSys], iSphericity, iMult),KaonMass);
		  firstTrackCut->AddCutMonitor(monitorYPtPass[anIter],monitorYPtFail[anIter]);
		}
	      if(iSys == kPPP || iSys == kAPAPAP) //P
		{
		  monitorPIDPass[anIter] = new AliFemtoCutMonitorParticlePID(Form("monitorPIDPass%sSt%iM%i", sysNames[iSys], iSphericity, iMult),2);//0-pion,1-kaon,2-proton
		  monitorPIDFail[anIter] = new AliFemtoCutMonitorParticlePID(Form("monitorPIDFail%sSt%iM%i", sysNames[iSys], iSphericity, iMult),2);//0-pion,1-kaon,2-proton
		  firstTrackCut->AddCutMonitor(monitorPIDPass[anIter],monitorPIDFail[anIter]);

		  monitorYPtPass[anIter] = new AliFemtoCutMonitorParticleYPt(Form("monitorYPtPass%sSt%iM%i", sysNames[iSys], iSphericity, iMult),ProtonMass);
		  monitorYPtFail[anIter] = new AliFemtoCutMonitorParticleYPt(Form("monitorYPtFail%sSt%iM%i", sysNames[iSys], iSphericity, iMult),ProtonMass);
		  firstTrackCut->AddCutMonitor(monitorYPtPass[anIter],monitorYPtFail[anIter]);
		}
	    }
        

        
	  bool doMinv = true;
        
	  // create m_inv distribution and add to the analysis
        
	  // get trio cut
	  AliFemtoTrioCut *trioCut = GetTrioCut((ESys)iSys);
	  detadphifctn[anIter] = new AliFemtoTrioDEtaDPhiFctn(Form("cdedpnocorr%stpcSt%iM%i", sysNames[iSys],iSphericity, iMult),29, 29);
	  detadphifctn[anIter]->SetTrioCut(trioCut);
	  trioAnalysis[anIter]->AddTrioFctn(detadphifctn[anIter]);
                        
	  // add analysis to the manager
	  Manager->AddAnalysis(trioAnalysis[anIter]);
	  anIter++;

      
	    }
	}
    }
  return Manager;
}

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis)
{
  AliFemtoEventReaderAODMultSelection* Reader = new AliFemtoEventReaderAODMultSelection();
  Reader->SetFilterBit(7);
  Reader->SetReadV0(1);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  Reader->SetReadMC(mcAnalysis);


  
  
  return Reader;
}

AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis)
{
  AliFemtoEventReaderAODChain* Reader = new AliFemtoEventReaderAODChain();
  Reader->SetFilterBit(7);
  Reader->SetReadV0(1);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  Reader->SetReadMC(mcAnalysis);
  
  return Reader;
}

AliFemtoEventReaderAODChain* GetReaderPP(bool mcAnalysis)
{
  /*
  AliFemtoEventReaderAODChain* Reader = new AliFemtoEventReaderAODChain();
  Reader->SetFilterMask(96);
  Reader->SetReadV0(true);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);

  Reader->SetUseAliEventCuts(true);
  Reader->SetIsPileUpEvent(true);
  Reader->SetUseMVPlpSelection(true);
  Reader->SetTrackPileUpRemoval(true);
  
  //Reader->SetMinPlpContribSPD(3);
  
  Reader->SetReadMC(mcAnalysis);
  */
  
  AliFemtoEventReaderNanoAODChain *Reader = new AliFemtoEventReaderNanoAODChain();
  Reader->SetFilterMask(filterbit);
  Reader->SetCovMatPresent(false);
  Reader->SetDCAglobalTrack(1); //false for FB7, true for the rest //we do not use DCA at all
  Reader->SetUseMultiplicity("MultSelection.RefMult08.Value");

  Reader->SetReadV0(kTRUE);
  Reader->SetReadCascade(kTRUE);

  return Reader;
}

AliFemtoTrioAnalysis* GetAnalysis(bool doEventMixing)
{
  AliFemtoTrioAnalysis *analysis = new AliFemtoTrioAnalysis();
  analysis->SetDoEventMixing(doEventMixing);
  // here one can put some additional analysis settings in the future
  
  return analysis;
}

AliFemtoSphericityEventCut* GetEventCut(int multMin,int multMax,int zVertMin,int zVertMax, double stMin, double stMax)
{
  AliFemtoSphericityEventCut *eventCut = new AliFemtoSphericityEventCut();
  eventCut->SetEventMult(multMin,multMax);
  eventCut->SetVertZPos(zVertMin,zVertMax);
  eventCut->SetStMin(stMin);
  eventCut->SetStMax(stMax);
  //eventCut->SetEPVZERO(-TMath::Pi()/2.,TMath::Pi()/2.);
  
  return eventCut;
}

AliFemtoMJTrackCut* GetTrackCut(AliFemtoTrio::EPart particle)
{

  AliFemtoMJTrackCut *particleCut = new AliFemtoMJTrackCut();
  particleCut->SetPt(0.15,2.5);

  if(particle == AliFemtoTrio::kPionPlus || particle==AliFemtoTrio::kPionMinus){
    particleCut->SetMostProbable(19);
    particleCut->SetMass(PionMass);
    particleCut->SetPt(0.2,2.5);
  }
  else  if(particle == AliFemtoTrio::kKaonPlus || particle==AliFemtoTrio::kKaonMinus){
    particleCut->SetMostProbable(20);
    particleCut->SetMass(KaonMass);
    particleCut->SetPt(0.3,2.5);
  }
  else if(particle == AliFemtoTrio::kProton || particle==AliFemtoTrio::kAntiProton)
    {
      particleCut->SetMostProbable(21);
      particleCut->SetMass(ProtonMass);
      particleCut->SetPt(0.5,2.5);
    }
  
  if(particle == AliFemtoTrio::kKaonPlus || particle == AliFemtoTrio::kPionPlus || particle == AliFemtoTrio::kProton)
    { particleCut->SetCharge( 1.0); }
  else
    { particleCut->SetCharge(-1.0); }


  double ProtonMass = 0.938272013;
  
  //particleCut->SetMostProbable(18);
  //particleCut->SetMass(ProtonMass);
  particleCut->SetEta(-0.8, 0.8);
  //particleCut->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  //particleCut->SetminTPCncls(80);
  //particleCut->SetRemoveKinks(kTRUE);
  //particleCut->SetLabel(kFALSE);
  //particleCut->SetMaxTPCChiNdof(4.0);
  //particleCut->SetMaxImpactXY(2.8);
  //particleCut->SetMaxImpactZ(3.2);
  particleCut->SetNsigma(3.0);
  particleCut->SetNsigma2(3.0);
  particleCut->SetNsigmaTPCTOF(kTRUE);
  particleCut->SetElectronRejection(kTRUE);
  //particleCut->SetCharge(particle == kEProton ? 1.0 : -1.0);
  


  return particleCut;
}


AliFemtoTrioCut* GetTrioCut(ESys system)
{
  AliFemtoTrioCut *trioCut = new AliFemtoTrioCut();
  
  
  return trioCut;
}


void GetParticlesForSystem(ESys system, AliFemtoTrio::EPart &firstParticle, AliFemtoTrio::EPart &secondParticle, AliFemtoTrio::EPart &thirdParticle)
{
  if(system == kPIpPIpPIp){
    firstParticle  = AliFemtoTrio::kPionPlus;
    secondParticle = AliFemtoTrio::kPionPlus;
    thirdParticle  = AliFemtoTrio::kPionPlus;
  }
  if(system == kPImPImPIm){
    firstParticle  = AliFemtoTrio::kPionMinus;
    secondParticle = AliFemtoTrio::kPionMinus;
    thirdParticle  = AliFemtoTrio::kPionMinus;
  }
  if(system == kPIpPImPIm){
    firstParticle  = AliFemtoTrio::kPionPlus;
    secondParticle = AliFemtoTrio::kPionMinus;
    thirdParticle  = AliFemtoTrio::kPionMinus;
  }
  if(system == kKpKpKp){
    firstParticle  = AliFemtoTrio::kKaonPlus;
    secondParticle = AliFemtoTrio::kKaonPlus;
    thirdParticle  = AliFemtoTrio::kKaonPlus;
  }
  if(system == kKmKmKm){
    firstParticle  = AliFemtoTrio::kKaonMinus;
    secondParticle = AliFemtoTrio::kKaonMinus;
    thirdParticle  = AliFemtoTrio::kKaonMinus;
  }
  if(system == kKpKmKm){
    firstParticle  = AliFemtoTrio::kKaonPlus;
    secondParticle = AliFemtoTrio::kKaonMinus;
    thirdParticle  = AliFemtoTrio::kKaonMinus;
  }
  if(system == kPPP){
    firstParticle  = AliFemtoTrio::kProton;
    secondParticle  = AliFemtoTrio::kProton;
    thirdParticle  = AliFemtoTrio::kProton;
  }
  if(system == kAPAPAP){
    firstParticle  = AliFemtoTrio::kAntiProton;
    secondParticle  = AliFemtoTrio::kAntiProton;
    thirdParticle  = AliFemtoTrio::kAntiProton;
  }
  if(system == kPAPAP){
    firstParticle  = AliFemtoTrio::kProton;
    secondParticle  = AliFemtoTrio::kAntiProton;
    thirdParticle  = AliFemtoTrio::kAntiProton;
  }
  if(system == kPKpKp)
    {
      firstParticle  = AliFemtoTrio::kProton;
      secondParticle  = AliFemtoTrio::kKaonPlus;
      thirdParticle  = AliFemtoTrio::kKaonPlus;
    }
  if(system == kPKmKm)
    {
      firstParticle  = AliFemtoTrio::kProton;
      secondParticle  = AliFemtoTrio::kKaonMinus;
      thirdParticle  = AliFemtoTrio::kKaonMinus;
    }
  if(system == kAPKpKp)
    {
      firstParticle  = AliFemtoTrio::kAntiProton;
      secondParticle  = AliFemtoTrio::kKaonPlus;
      thirdParticle  = AliFemtoTrio::kKaonPlus;
    }
  if(system == kAPKmKm)
    {
      firstParticle  = AliFemtoTrio::kAntiProton;
      secondParticle  = AliFemtoTrio::kKaonMinus;
      thirdParticle  = AliFemtoTrio::kKaonMinus;
    }
if(system == kPIpKpKp)
    {
      firstParticle  = AliFemtoTrio::kPionPlus;
      secondParticle  = AliFemtoTrio::kKaonPlus;
      thirdParticle  = AliFemtoTrio::kKaonPlus;
    }
 if(system == kPIpKmKm)
    {
      firstParticle  = AliFemtoTrio::kPionPlus;
      secondParticle  = AliFemtoTrio::kKaonMinus;
      thirdParticle  = AliFemtoTrio::kKaonMinus;
    }
 if(system == kPPIpPIp)
    {
      firstParticle  = AliFemtoTrio::kProton;
      secondParticle  = AliFemtoTrio::kPionPlus;
      thirdParticle  = AliFemtoTrio::kPionPlus;
    }
 if(system == kPPImPIm)
    {
      firstParticle  = AliFemtoTrio::kProton;
      secondParticle  = AliFemtoTrio::kPionMinus;
      thirdParticle  = AliFemtoTrio::kPionMinus;
    }
 if(system == kPIpKpP)
    {
      firstParticle  = AliFemtoTrio::kPionPlus;
      secondParticle  = AliFemtoTrio::kKaonPlus;
      thirdParticle  = AliFemtoTrio::kProton;
    }

  
}
