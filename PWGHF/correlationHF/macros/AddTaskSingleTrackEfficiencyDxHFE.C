//DEFINITION OF A FEW CONSTANTS
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


/*__|______________________________________________________________________________|
 |                              -----Info(i)-----                                  |
 |                                                                                 |
 |  AddTask for single track  efficiecy                                            |
 |                                                                                 |
 |   ESDs<-->AODs (ON/OFF)                                                         |
 |   filterbit<-->(Options)                                                        |
 |                                                      Authors:                   |
 |_____________________________________________________________________________|___*/



const Int_t    mintrackrefsTPC = 5 ;
const Int_t    mintrackrefsITS = 4 ;
const Int_t    mintrackrefsTOF = 0;
const Int_t    mintrackrefsMUON = 0 ;
const Int_t    minclustersTPC = 70 ;
const Int_t    minclustersITS = 2 ;
const Bool_t   TPCRefit = kTRUE;
const Bool_t   ITSRefit = kFALSE;
const Bool_t   charge  = kTRUE ;
const Int_t    fBit;


// cuts
const Double_t etamin  = -0.9 ;
const Double_t etamax  =  0.9 ;
const Double_t ptmin =  0.0 ;
const Double_t ptmax =  24.0 ;
const Double_t phimin = -2*TMath::Pi();
const Double_t phimax = 2*TMath::Pi();
const Double_t thetamin = 0;
const Double_t thetamax = TMath::Pi();
const Double_t zvtxmin =  -10.0 ;
const Double_t zvtxmax =  10.0 ;
const Double_t dcamin = 0;  // micron
const Double_t dcamax = 600;  // micron


// Mutliplicity
const Float_t multmin_0_20 = 0;
const Float_t multmax_0_20 = 20;
const Float_t multmin_20_50 = 20;
const Float_t multmax_20_50 = 50;
const Float_t multmin_50_102 = 50;
const Float_t multmax_50_102 = 102;


//  Pt Range
Double_t ptmin_0_6   = 0.0;
Double_t ptmax_0_6   = 6.0;
Double_t ptmin_6_8   = 6.0;
Double_t ptmax_6_8   = 8.0;
Double_t ptmin_8_16  = 8.0;
Double_t ptmax_8_16  = 16.0;
Double_t ptmin_16_24 = 16.0;
Double_t ptmax_16_24 = 24.0;




AliCFSingleTrackEfficiencyTask *AddSingleTrackEfficiencyTask(const Bool_t readAOD = 0,   TString TrackCutsfilename =  "")
{ 

  Info("AliCFSingleTrackEfficiencyTask","SETUP CONTAINER");

  const Int_t nvar   = 5 ; //number of variables on the grid:pt,y,phi
  UInt_t nstep = 7 ; //number of selection steps MC

  const UInt_t ipt = 0;
  const UInt_t iy  = 1;
  const UInt_t iphi = 2;
  const UInt_t itheta  = 3;
  const UInt_t izvtx  = 4;
  // const UInt_t imult = 5;
  
  

  //A1. Bins variation by hand for pt
  const Int_t nbinpt_0_6  = 24 ; //bins in pt from 0 to 6 GeV
  const Int_t nbinpt_6_8  = 4 ; //bins in pt from 6 to 8 GeV
  const Int_t nbinpt_8_16  = 8; //bins in pt from 8 to 16 GeV
  const Int_t nbinpt_16_24  = 1 ; //bins in pt from 16 to 24 GeV


  //A2. Bins variation by hand for other variables
  const Int_t nbin2  = 18 ; //bins in eta
  const Int_t nbin3  = 18 ; //bins in phi
  const Int_t nbin4  = 18 ; //bins in theta
  const Int_t nbin5  = 20 ; //bins in zvtx

  //A3. Bins for multiplicity
  const Int_t nbinmult = 48;  //bins in multiplicity (total number)	
  const Int_t nbinmult_0_20 = 20; //bins in multiplicity between 0 and 20
  const Int_t nbinmult_20_50 = 15; //bins in multiplicity between 20 and 50
  const Int_t nbinmult_50_102 = 13; //bins in multiplicity between 50 and 102
  

  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbinpt_0_6+nbinpt_6_8+nbinpt_8_16+nbinpt_16_24;
  iBin[1]=nbin2;
  iBin[2]=nbin3;
  iBin[3]=nbin4;
  iBin[4]=nbin5;
  //  iBin[5]=nbinmult_0_20 +  nbinmult_20_50 +  nbinmult_50_102  ;

  
  //arrays for lower bounds :
  Double_t *binLimpT=new Double_t[iBin[0]+1];
  Double_t *binLim2=new Double_t[iBin[1]+1];
  Double_t *binLim3=new Double_t[iBin[2]+1];
  Double_t *binLim4=new Double_t[iBin[3]+1];
  Double_t *binLim5=new Double_t[iBin[4]+1];
  //  Double_t *binLimmult=new Double_t[iBin[5]+1];

  // pt
  for(Int_t i=0; i<=nbinpt_0_6; i++) binLimpT[i]=(Double_t)ptmin_0_6 + (ptmax_0_6-ptmin_0_6)/nbinpt_0_6*(Double_t)i ; 
  for(Int_t i=0; i<=nbinpt_6_8; i++) binLimpT[i+nbinpt_0_6]=(Double_t)ptmin_6_8 + (ptmax_6_8-ptmin_6_8)/nbinpt_6_8*(Double_t)i ; 
  for(Int_t i=0; i<=nbinpt_8_16; i++) binLimpT[i+nbinpt_0_6+nbinpt_6_8]=(Double_t)ptmin_8_16 + (ptmax_8_16-ptmin_8_16)/nbinpt_8_16*(Double_t)i ; 
  for(Int_t i=0; i<=nbinpt_16_24; i++) binLimpT[i+nbinpt_0_6+nbinpt_6_8+nbinpt_8_16]=(Double_t)ptmin_16_24 + (ptmax_16_24-ptmin_16_24)/nbinpt_16_24*(Double_t)i ; 


  // Other Variables
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)etamin  + (etamax-etamin)  /nbin2*(Double_t)i ;
  for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)phimin + (phimax-phimin)/nbin3*(Double_t)i ;
  for(Int_t i=0; i<=nbin4; i++) binLim4[i]=(Double_t)thetamin  + (thetamax-thetamin)  /nbin4*(Double_t)i ;
  for(Int_t i=0; i<=nbin5; i++) binLim5[i]=(Double_t)zvtxmin  + (zvtxmax-zvtxmin)  /nbin5*(Double_t)i ;

  
  /*
  // multiplicity bining..
  for(Int_t i=0; i<=nbinmult_0_20; i++) binLimmult[i]=(Double_t)multmin_0_20 + (multmax_0_20-multmin_0_20)/nbinmult_0_20*(Double_t)i ; 
  for(Int_t i=0; i<=nbinmult_20_50; i++) binLimmult[i+nbinmult_0_20]=(Double_t)multmin_20_50 + (multmax_20_50-multmin_20_50)/nbinmult_20_50*(Double_t)i ; 
  for(Int_t i=0; i<=nbinmult_50_102; i++) binLimmult[i+nbinmult_0_20+nbinmult_20_50]=(Double_t)multmin_50_102 + (multmax_50_102-multmin_50_102)/nbinmult_50_102*(Double_t)i ; 
  */

 
  //Container  
  AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
  container -> SetBinLimits(ipt,binLimpT);//pt
  container -> SetBinLimits(iy,binLim2);//eta
  container -> SetBinLimits(iphi,binLim3);//phi
  container -> SetBinLimits(itheta,binLim4);//theta
  container -> SetBinLimits(izvtx,binLim5);//Zvtx
  //  container -> SetBinLimits(imult,binLimmult);

  //Variable Titles 
  container -> SetVarTitle(ipt,"pt");
  container -> SetVarTitle(iy, "y");
  container -> SetVarTitle(iphi,"phi");
  container -> SetVarTitle(itheta, "theata");
  container -> SetVarTitle(izvtx, "Zvtx");
  // container -> SetVarTitle(imult, "Multiplicity");

  //Variable Titles 
  container -> SetStepTitle(0, " MC Particle with Generated Cuts");
  container -> SetStepTitle(1, " MC Particle with Kine Acceptance Cuts");
  container -> SetStepTitle(2, " MC Particle with Track Ref Acceptance Cuts");
  container -> SetStepTitle(3, " Total Reconstructed  Particle ");
  container -> SetStepTitle(4, " Reco Particle With Kine Acceptance Cuts");
  container -> SetStepTitle(5, " Reco Particle to MC True pt particles ");
  container -> SetStepTitle(6, " Reco Particle With Quality Cuts");


  // SET TLIST FOR QA HISTOS
  TList* qaList = new TList();
  TObjArray* emptyList = new TObjArray(0);
    
  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* man = new AliCFManager();
  
  man->SetNStepEvent(2);
  man->SetEventContainer(container);
  man->SetEventCutsList(0,emptyList);//evtmcList);
  man->SetEventCutsList(1,emptyList);//evtrecoList);
  
  man->SetParticleContainer(container);
  man->SetParticleCutsList(0,emptyList);//mcGenList);
  man->SetParticleCutsList(1,emptyList);//mcKineList);
  man->SetParticleCutsList(2,emptyList);//mcaccList);
  man->SetParticleCutsList(3,emptyList);//evtrecoPureList);
  man->SetParticleCutsList(4,emptyList);//recKineList);
  man->SetParticleCutsList(5,emptyList);//fPIDCutList);
  man->SetParticleCutsList(6,emptyList);//fPIDCutList);

  
  // Simulated particle & event cuts
  AliSingleTrackEffCuts* cuts = new AliSingleTrackEffCuts();
  cuts->SetPtRange(ptmin,ptmax);
  cuts->SetEtaRange(etamin,etamax);
  //cuts->SetYRange(Ymin,Ymax);
  cuts->SetIsCharged(charge);
  cuts->SetMinVtxContr(-1);
  cuts->SetMaxVtxZ(zvtxmax);
  cuts->SetNumberOfClusters(mintrackrefsITS,mintrackrefsTPC,mintrackrefsTOF,mintrackrefsMUON);
  cuts->SetTriggerMask(AliVEvent::kAnyINT);
  cuts->SetIsAOD(readAOD);
  
  


 //  Track Quality cuts from HF Correlations Class (Associated Track Cuts)
  AliESDtrackCuts* QualityCuts=new AliESDtrackCuts();
  TString loadLibraries="LoadLibraries.C"; 
  gROOT->LoadMacro(loadLibraries.Data());
  LoadLibraries();
  
  
  Bool_t AssTrackCuts=kFALSE;
  TFile* filecuts=new TFile(TrackCutsfilename.Data());
  
  if( TrackCutsfilename.EqualTo("") ) {
    AssTrackCuts=kFALSE; 
  } else {
    AssTrackCuts=kTRUE; 
    if(!filecuts->IsOpen()){
      cout<<"Track cut object file not found: exit"<<endl;
      return;
    }	 
  }
  
  
  if(AssTrackCuts) {
    AliHFAssociatedTrackCuts* HFAssTrackCuts=new AliHFAssociatedTrackCuts();
    HFAssTrackCuts = (AliHFAssociatedTrackCuts*)filecuts->Get("AssociatedCuts");
    HFAssTrackCuts->SetName("AssociatedCuts");
    HFAssTrackCuts->PrintAll();
    if(!HFAssTrackCuts){
      cout<<"Specific associated track as Quality cuts is not found"<<endl;
      return;
    } 
    //QualityCuts = HFAssTrackCuts->GetESDTrackCuts(); //Invalid option
  }
  else{
    // Track Quality cuts for general methods
    QualityCuts->SetRequireSigmaToVertex(kFALSE);
    QualityCuts->SetDCAToVertex2D(kFALSE);
    QualityCuts->SetMinNClustersTPC(minclustersTPC);
    QualityCuts->SetMinNClustersITS(minclustersITS);
    QualityCuts->SetMaxChi2PerClusterTPC(4);
    QualityCuts->SetMaxDCAToVertexZ(1);
    QualityCuts->SetMaxDCAToVertexXY(1);
    QualityCuts->SetRequireTPCRefit(TPCRefit);
    QualityCuts->SetRequireITSRefit(ITSRefit);
  }
  
  
  //CREATE THE TASK
  printf("CREATE CF Single track task\n");
  
  AliCFSingleTrackEfficiencyTask *task = new AliCFSingleTrackEfficiencyTask("AliCFSingleTrackEfficiencyTask",QualityCuts,cuts);
  task->SetFilterBit(kTRUE);
  task->SetFilterType(0);
  task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  task->SetCFManager(man); //here is set the CF manager
  
  
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }
  
  // This task requires an ESD or AOD input handler and an AOD output handler.
  // Check this using the analysis manager.
  //===============================================================================
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddSingleTrackEfficiencyTask", "AliCFSingleTrackEfficiency task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }   
  
  
  printf(" Create the output container\n");

  // Create and connect containers for input/output
    
  
  // ----- output data -----
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString input1name="cchain0";
  TString output1name="ctree0", output2name="chist0", output3name="ccontainer0",output4name="clist0";
  TString suffix="";
  outputfile += ":PWG_PP_CFtaskSingleTrack";
  outputfile += suffix;
  output1name += suffix;
  output2name += suffix;
  output3name += suffix;
  output4name += suffix;


  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();


  // ----- output data -----
  
  //slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
  //  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer(output1name, TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  //now comes user's output objects :

  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output2name, TH1I::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output3name, AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // output QA histograms
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output4name, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  
  return task;
}
