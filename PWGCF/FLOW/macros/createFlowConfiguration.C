void createFlowConfiguration() {
//Macro that configures the flow cut objects and returns them. 
//This is used by the AddTaskPIDFlowSP macro, and can be used 
//for local and grid analysis, even on the train
//In the case of a grid analysis, this macro needs to be placed 
//on the file catalogue and then the relevant location needs 
//to be provided to the AddTask as an argument
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWGCFflowBase.so");
  gSystem->Load("libPWGCFflowTasks.so");
}

//_________________________________________________________//
void defineCentralities(Int_t &nCentralities, 
			Double_t *gCentrality) {
  //Part of the code that defines the centrality percentiles 
  //to be used for Pb-Pb
  //static const Int_t gCentralities = 24;
  //Double_t fCentrality[gCentralities] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,40,50};  
  static const Int_t gCentralities = 2;
  Double_t fCentrality[gCentralities] = {0,50};
  for(Int_t iCentrality = 0; iCentrality < gCentralities; iCentrality++) 
    gCentrality[iCentrality] = fCentrality[iCentrality];

  nCentralities = gCentralities;
}

//_________________________________________________________//
AliFlowEventCuts *createFlowEventCutObject(Int_t gCentralityMin = -1,
					   Int_t gCentralityMax = -1,
					   Bool_t isPbPb = kTRUE,
					   Bool_t doQA = kFALSE) {
  //Part of the code that creates the event cut objects
  Bool_t is2011 = kTRUE;
  AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO;
  Double_t gVertexZmin = -10., gVertexZmax = 10.;
  
  //Create the event cut objects
  AliFlowEventCuts *cutsEvent = new AliFlowEventCuts(Form("eventcutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
  if(isPbPb) {
    cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax);
    cutsEvent->SetCentralityPercentileMethod(gCentralityEstimator);
    cutsEvent->SetLHC11h(is2011);
    cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
  }
  else 
    cutsEvent->SetCheckPileup(kFALSE);
  
  cutsEvent->SetPrimaryVertexZrange(gVertexZmin,gVertexZmax);
  cutsEvent->SetQA(doQA);

  //return the object
  return cutsEvent;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowRPCutObject(Int_t gCentralityMin = -1,
					Int_t gCentralityMax = -1,
					Bool_t isVZERO = kFALSE,
					Bool_t doQA = kFALSE) {
  //Part of the code that creates the RP cut objects
  Double_t gEtaMin = -0.8, gEtaMax = 0.8;
  Bool_t is2011 = kTRUE;
  Int_t gMinTPCdedx = 10;
  Int_t gAODfilterBit = 768;

  //Create the event cut objects
  AliFlowTrackCuts *cutsRP = new AliFlowTrackCuts(Form("rpCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
  if(!isVZERO) {
    cutsRP->SetPtRange(0.2,5.);
    cutsRP->SetEtaRange(gEtaMin,gEtaMax);
    cutsRP->SetMinimalTPCdedx(gMinTPCdedx);
    cutsRP->SetAODfilterBit(gAODfilterBit);
  }
  else if(isVZERO) { // use vzero sub analysis
    if(!is2011) 
      cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2010(); 
    if(is2011)  
      cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2011(); 
  }//VZERO
  
  cutsRP->SetQA(doQA);
  
  //return the object
  return cutsRP;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowPOICutObject(Int_t gCentralityMin = -1,
					 Int_t gCentralityMax = -1,
					 TString particleSpecies = "Pion",
					 TString gQvector = "Qa",
					 Double_t gEtaGap = 0.0,
					 Bool_t isVZERO = kFALSE,
					 Bool_t isPbPb = kTRUE,
					 Bool_t doQA = kFALSE) {
  //Part of the code that creates the POI cut objects
  Double_t gEtaMin = -0.8, gEtaMax = 0.8;
  Bool_t is2011 = kTRUE;
  Int_t gMinTPCdedx = 10;
  Int_t gAODfilterBit = 768;
  Bool_t isPID = kTRUE;
  Int_t gCharge = 0;  
  Double_t gDCAvtxXY = -1.;
  Double_t gDCAvtxZ = -1.;
  AliFlowTrackCuts::PIDsource sourcePID=AliFlowTrackCuts::kTOFbayesian;
  AliPID::EParticleType particleType = AliPID::kPion;

  if(particleSpecies.Contains("Pion")) {
    //cout<<"(Panos) Pions"<<endl;
    particleType = AliPID::kPion;
  }
  else if(particleSpecies.Contains("Kaon")) {
    //cout<<"(Panos) Kaons"<<endl;
    particleType = AliPID::kKaon;
  }
  else if(particleSpecies.Contains("Proton")) {
    //cout<<"(Panos) Protons"<<endl;
    particleType = AliPID::kProton;
  }
    
  AliFlowTrackCuts *cutsPOI = new AliFlowTrackCuts(Form("poiCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));

  //for 2010 data to use old TPC PID Response instead of the official one
  if(!is2011) 
    cutsPOI->GetBayesianResponse()->ForceOldDedx(); 
  
  if(!isVZERO) {
    if(gQvector.Contains("Qa")) {
      //cout<<"(Panos): "<<gQvector.Data()<<endl;
      cutsPOI->SetEtaRange(+0.5*gEtaGap,gEtaMax );
    }
    else if(gQvector.Contains("Qb")) {
      //cout<<"(Panos): "<<gQvector.Data()<<endl;
      cutsPOI->SetEtaRange(gEtaMin,-0.5*gEtaGap);
    }
  }
  else if(isVZERO)
    cutsPOI->SetEtaRange(gEtaMin,gEtaMax);
  
  cutsPOI->SetAcceptKinkDaughters(kFALSE);
  
  if(isPID) {
    cutsPOI->SetPID(particleType, sourcePID);
    
    if(isPbPb) 
      cutsPOI->SetPriors((gCentralityMin + gCentralityMax)*0.5);
  }
  
  if (gCharge!=0) cutsPOI->SetCharge(gCharge);

  cutsPOI->SetMinNClustersTPC(70);
  cutsPOI->SetPtRange(0.2,6.);
  cutsPOI->SetAODfilterBit(gAODfilterBit);
  if(gDCAvtxXY > 0) 
    cutsPOI->SetMaxDCAToVertexXY(gDCAvtxXY);
  if(gDCAvtxZ > 0)
    cutsPOI->SetMaxDCAToVertexZ(gDCAvtxZ);
  cutsPOI->SetMinimalTPCdedx(gMinTPCdedx);
  cutsPOI->SetQA(doQA);

  //return the object
  return cutsPOI;
}
