//DEFINITION OF A FEW CONSTANTS
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

const Float_t phimin = 0.;
const Float_t phimax = 2.*TMath::Pi();
const Float_t etamin = -0.9;
const Float_t etamax = 0.9;

const Int_t   mintrackrefsTPC = 0;
const Int_t   mintrackrefsITS = 0;

/* AddTask setting up the containers and adding the tasks for computing the efficiency of the Hybrid Tracks. Tasks works on ESDs and AODs (datatype established automatically).
Input:
hybridTrackDef: "2010" or "2011", indicates whether the 2010 definition or the 2011 definition of the hybrid tracks should be used. Only important in ESD.
prodType: Data or MC period, 
beamType: "pp", "pPb" or "PbPb"
triggerMask: Event Selection Trigger Mask
bSelHijingParticles: kTRUE, select only particles from Hijing event. kFALSE, select all particles of the cocktail
usePythiaxsec: kFALSE, do not use the Pythia cross section information file. This might reduce the number of errors and file opening issues. 
*/

void AddTaskHybridTrackEfficiency(TString hybridTrackDef = "2011", char *prodType = "LHC11h", TString beamType = "PbPb", UInt_t triggerMask = AliVEvent::kMB, Bool_t bSelHijingParticles = kFALSE, Bool_t usePythiaxsec = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHybridTrackEfficiency", "No analysis manager to connect to.");
    return NULL;
  }  
  TString dataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  if(dataType=="AOD") {
    AddTaskHybridTrackEfficiencyAOD(prodType, beamType, triggerMask, bSelHijingParticles, usePythiaxsec);
  }
  if(dataType=="ESD") {
    if(hybridTrackDef.Contains("2010"))
      AddTaskHybridTrackEfficiencyESD2010(prodType, beamType, triggerMask, bSelHijingParticles, usePythiaxsec);
    if(hybridTrackDef.Contains("2011"))
      AddTaskHybridTrackEfficiencyESD2011(prodType, beamType, triggerMask, bSelHijingParticles, usePythiaxsec);
  }
}

void AddTaskHybridTrackEfficiencyAOD(char *prodType = "LHC11h", TString beamType = "PbPb", UInt_t triggerMask = AliVEvent::kMB, Bool_t bSelHijingParticles = kFALSE, Bool_t usePythiaxsec = kFALSE)
{
  Int_t filterMask1; //standard global tracks
  Int_t filterMask2; //complementary tracks
  Int_t filterMask; //the sum: hybrid tracks

  TString strRunPeriod = TString(prodType);
  strRunPeriod.ToLower();

  if (strRunPeriod == "lhc10h" || strRunPeriod == "lhc11h" ||
      strRunPeriod == "lhc12a" || strRunPeriod == "lhc12b" || strRunPeriod == "lhc12c" || strRunPeriod == "lhc12d" ||
      strRunPeriod == "lhc12e" || strRunPeriod == "lhc12f" || strRunPeriod == "lhc12g" || strRunPeriod == "lhc12g" ||
      strRunPeriod == "lhc12h" || strRunPeriod == "lhc12i" ||
      strRunPeriod == "lhc13b" || strRunPeriod == "lhc13c" || strRunPeriod == "lhc13d" || strRunPeriod == "lhc13e" ||
      strRunPeriod == "lhc13f" || strRunPeriod == "lhc13g" ||
      strRunPeriod == "lhc12a15e" || strRunPeriod == "lhc13b4" || strRunPeriod == "lhc13b4_fix" ||
      strRunPeriod == "lhc13b4_plus" || strRunPeriod == "lhc12a15f" || strRunPeriod.Contains("lhc12a17") || strRunPeriod.Contains("lhc14a1")) {
    filterMask  = 768;
    filterMask1 = 256;
    filterMask2 = 512;
    bIncludeNoITS = kFALSE;
    if(strRunPeriod == "lhc10h") bIncludeNoITS = kTRUE;
  }
  else if (strRunPeriod == "lhc11a" || strRunPeriod == "lhc10hold" || strRunPeriod == "lhc12a15a" || strRunPeriod.Contains("lhc11a2")) {
    filterMask  = 272;
    filterMask1 = 16;
    filterMask2 = 256;
    bIncludeNoITS = kTRUE;
  }
  else {
    ::Error("AddTaskHybridTrackEfficiency","Period string not of predefined type. Add it to the list in this macro.");
    return NULL;
  }

  AliPWG4HighPtSpectra   *taskSpectraSUM   = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,0,1,triggerMask,"CENT_10_SUM",bSelHijingParticles,usePythiaxsec,filterMask ); //Sum: the hybrid tracks.
  AliPWG4HighPtSpectra   *taskSpectraRESTR = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,0,5,triggerMask,"CENT_10_RESTR",bSelHijingParticles,usePythiaxsec,filterMask1); //High quality tracks
  AliPWG4HighPtSpectra   *taskSpectraNOSPD = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,7,5,triggerMask,"CENT_10_COMPL",bSelHijingParticles,usePythiaxsec,filterMask2); //Complementary tracks
  if(beamType=="PbPb" || beamType=="pPb") { //also vary the centrality
    for(int cent = 0; cent<4; cent++) {
      AliPWG4HighPtSpectra *taskSpectraCent0 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,0,1,triggerMask,Form("CENT_%d_SUM",cent),bSelHijingParticles,usePythiaxsec,filterMask );
      AliPWG4HighPtSpectra *taskSpectraCent0 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,0,5,triggerMask,Form("CENT_%d_RESTR",cent),bSelHijingParticles,usePythiaxsec,filterMask1 );
      AliPWG4HighPtSpectra *taskSpectraCent0 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,7,5,triggerMask,Form("CENT_%d_COMPL",cent),bSelHijingParticles,usePythiaxsec,filterMask2 );
    }
  }
}

void AddTaskHybridTrackEfficiencyESD2010(char *prodType = "LHC10h", TString beamType = "PbPb", UInt_t triggerMask = AliVEvent::kMB, Bool_t bSelHijingParticles = kFALSE, Bool_t usePythiaxsec = kFALSE)
{
  AliPWG4HighPtSpectra *taskSpectra00cent10 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,0,0,triggerMask,"CENT_10_RAA"  ,bSelHijingParticles,usePythiaxsec); // RAA track cuts
  AliPWG4HighPtSpectra *taskSpectra01cent10 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,0,1,triggerMask,"CENT_10_RESTR",bSelHijingParticles,usePythiaxsec); // High quality tracks
  AliPWG4HighPtSpectra *taskSpectra70cent10 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,7,0,triggerMask,"CENT_10_NOITS",bSelHijingParticles,usePythiaxsec); // Complementary tracks. Subtype: no ITS refit
  AliPWG4HighPtSpectra *taskSpectra71cent10 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,7,1,triggerMask,"CENT_10_NOSPD",bSelHijingParticles,usePythiaxsec); // Complementary tracks. Subtype: no SPD requirement

  if(beamType=="PbPb" || beamType=="pPb") { //also vary the centrality
    for(int cent=0; cent<4; cent++) {
      AliPWG4HighPtSpectra *taskSpectra00 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,0,0,triggerMask,Form("CENT_%d_RAA"  ,cent),bSelHijingParticles,usePythiaxsec);
      AliPWG4HighPtSpectra *taskSpectra01 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,0,1,triggerMask,Form("CENT_%d_RESTR",cent),bSelHijingParticles,usePythiaxsec);
      AliPWG4HighPtSpectra *taskSpectra70 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,7,0,triggerMask,Form("CENT_%d_NOITS",cent),bSelHijingParticles,usePythiaxsec);
      AliPWG4HighPtSpectra *taskSpectra71 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,7,1,triggerMask,Form("CENT_%d_NOSPD",cent),bSelHijingParticles,usePythiaxsec);
    }
  }
}

void AddTaskHybridTrackEfficiencyESD2011(char *prodType = "LHC10h", TString beamType = "PbPb", UInt_t triggerMask = AliVEvent::kMB, Bool_t bSelHijingParticles = kFALSE, Bool_t usePythiaxsec = kFALSE)
{
  AliPWG4HighPtSpectra *taskSpectra00cent10 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,0,0,triggerMask,"CENT_10_RAA"  ,bSelHijingParticles,usePythiaxsec); // RAA track cuts
  AliPWG4HighPtSpectra *taskSpectra01cent10 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,0,5,triggerMask,"CENT_10_RESTR",bSelHijingParticles,usePythiaxsec); // High quality tracks
  AliPWG4HighPtSpectra *taskSpectra71cent10 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,7,5,triggerMask,"CENT_10_NOSPD",bSelHijingParticles,usePythiaxsec); // Complementary tracks. Subtype: no SPD requirement

  if(beamType=="PbPb" || beamType=="pPb") { //also vary the centrality
    for(int cent=0; cent<4; cent++) {
      AliPWG4HighPtSpectra *taskSpectra00 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,0,0,triggerMask,Form("CENT_%d_RAA"  ,cent),bSelHijingParticles,usePythiaxsec);
      AliPWG4HighPtSpectra *taskSpectra01 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,0,5,triggerMask,Form("CENT_%d_RESTR",cent),bSelHijingParticles,usePythiaxsec);
      AliPWG4HighPtSpectra *taskSpectra71 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,cent,7,5,triggerMask,Form("CENT_%d_NOSPD",cent),bSelHijingParticles,usePythiaxsec);
    }
  }
}

void AddTaskHybridTrackEfficiencyQA_AOD_train(char *prodType = "LHC11h", TString beamType = "PbPb", UInt_t triggerMask = AliVEvent::kMB, Bool_t bSelHijingParticles = kFALSE, Bool_t usePythiaxsec = kFALSE)
{
  Int_t filterMask1; //standard global tracks
  Int_t filterMask2; //complementary tracks
  Int_t filterMask; //the sum: hybrid tracks

  TString strRunPeriod = TString(prodType);
  strRunPeriod.ToLower();

  if (strRunPeriod == "lhc10h" || strRunPeriod == "lhc11h" ||
      strRunPeriod == "lhc12a" || strRunPeriod == "lhc12b" || strRunPeriod == "lhc12c" || strRunPeriod == "lhc12d" ||
      strRunPeriod == "lhc12e" || strRunPeriod == "lhc12f" || strRunPeriod == "lhc12g" || strRunPeriod == "lhc12g" ||
      strRunPeriod == "lhc12h" || strRunPeriod == "lhc12i" ||
      strRunPeriod == "lhc13b" || strRunPeriod == "lhc13c" || strRunPeriod == "lhc13d" || strRunPeriod == "lhc13e" ||
      strRunPeriod == "lhc13f" || strRunPeriod == "lhc13g" ||
      strRunPeriod == "lhc12a15e" || strRunPeriod == "lhc13b4" || strRunPeriod == "lhc13b4_fix" ||
      strRunPeriod == "lhc13b4_plus" || strRunPeriod == "lhc12a15f" || strRunPeriod.Contains("lhc12a17") || strRunPeriod.Contains("lhc14a1")) {
    filterMask  = 768;
    filterMask1 = 256;
    filterMask2 = 512;
  }
  else if (strRunPeriod == "lhc11a" || strRunPeriod == "lhc10hold" || strRunPeriod == "lhc12a15a" || strRunPeriod.Contains("lhc11a2")) {
    filterMask  = 272;
    filterMask1 = 16;
    filterMask2 = 256;
  }
  else {
    ::Error("AddTaskHybridTrackEfficiency","Period string not of predefined type. Add it to the list in this macro.");
    return NULL;
  }

  AliPWG4HighPtSpectra   *taskSpectraSUM   = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,0,1,triggerMask,"SUM"  ,bSelHijingParticles,usePythiaxsec,filterMask ); //Sum: the hybrid tracks.
  AliPWG4HighPtSpectra   *taskSpectraRESTR = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,0,5,triggerMask,"RESTR",bSelHijingParticles,usePythiaxsec,filterMask1); //With SPD and ITS refit.
  AliPWG4HighPtSpectra   *taskSpectraNOSPD = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,10,7,5,triggerMask,"NOSPD",bSelHijingParticles,usePythiaxsec,filterMask2); //Only ITS refit, not SPD.
  if(beamType=="PbPb" || beamType=="pPb") { //also vary the centrality
    AliPWG4HighPtSpectra *taskSpectraCent0 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,0 ,0,1,triggerMask,"CENT0",bSelHijingParticles,usePythiaxsec,filterMask );
    AliPWG4HighPtSpectra *taskSpectraCent1 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,1 ,0,1,triggerMask,"CENT1",bSelHijingParticles,usePythiaxsec,filterMask );
    AliPWG4HighPtSpectra *taskSpectraCent2 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,2 ,0,1,triggerMask,"CENT2",bSelHijingParticles,usePythiaxsec,filterMask );
    AliPWG4HighPtSpectra *taskSpectraCent3 = ConfigureTaskPWG4HighPtSpectra(prodType,beamType,3 ,0,1,triggerMask,"CENT3",bSelHijingParticles,usePythiaxsec,filterMask );
  }
}

AliPWG4HighPtSpectra* ConfigureTaskPWG4HighPtSpectra(char *prodType = "LHC10e14", TString beamType = "PbPb", Int_t centClass = 0, 
						     Int_t trackType = 0, Int_t cuts = 0, UInt_t triggerMask = AliVEvent::kMB,
						     TString taskName, Bool_t bSelHijingParticles = kFALSE,
                                                     Bool_t usePythiaxsec = kTRUE, Int_t filterMask = 0)
{

  /*
    trackType: 0 = global
               1 = TPC stand alone
               2 = TPC stand alone constrained to SPD vertex
               4 = TPC stand alone constrained to SPD vertex with QA track selection on global tracks
               5 = Hybrid tracks: constrained TPConly for which no tight ITS is available
               6 = Hybrid tracks: constrained loose global for which no tight ITS is available
    cuts:      0 (global) = standard ITSTPC2010 a la RAA analysis
               1 (global) = ITSrefit, no SPD requirements -> standard for jet analysis
               2 (global) = ITSrefit + no hits in SPD
               3 (global) = standard ITS tight cuts with nCrossed rows cut for hybrid tracks
               0 (TPC)    = standard TPC + NClusters>70
               1 (TPC)    = standard TPC + NClusters>0 --> to study new TPC QA recommendations
               0 (hybrid 5) = constrained TPConly for which no tight ITS is available
               0 (hybrid 6) = constrained loose global for which no tight ITS is available
   */

  // Creates HighPtSpectra analysis task and adds it to the analysis manager.
  
  //Load common track cut class
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");

  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHybridTrackEfficiency", "No analysis manager to connect to.");
    return NULL;
  }  
  TString dataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHybridTrackEfficiency", "This task requires an input event handler");
    return NULL;
  }  

  // C. Create the task, add it to manager.
  //===========================================================================
 //CONTAINER DEFINITION
  Info("AliPWG4HighPtSpectra","SETUP CONTAINER");
  //the sensitive variables, their indices
  UInt_t ipt   = 0;
  UInt_t iphi  = 1;
  UInt_t ieta  = 2;
  UInt_t incls = 3;

  //Setting up the container grid... 
  UInt_t nstep = 4; //Steps/Modes for containers
  Int_t kStepReconstructed          = 0;
  Int_t kStepSecondaries            = 1;
  Int_t kStepReconstructedMC        = 2;
  Int_t kStepMCAcceptance           = 3;

  //redefine pt ranges in case of Jet-Jet production
  Float_t ptBinEdges[4][4];
  Float_t ptmin =  0.0 ;
  Float_t ptmax =  100.0 ;

  ptBinEdges[0][0] = 2.;
  ptBinEdges[0][1] = 0.2;
  ptBinEdges[1][0] = 6.;
  ptBinEdges[1][1] = 0.5;
  ptBinEdges[2][0] = 20.;
  ptBinEdges[2][1] = 2.;
  ptBinEdges[3][0] = 100.;
  ptBinEdges[3][1] = 5.;
  
  const Int_t nvar   = 4; //number of variables on the grid: pt:phi:eta:NClsIter1

  const Int_t nbin11 = round((ptBinEdges[0][0]-ptmin)/ptBinEdges[0][1]);
  const Int_t nbin12 = round((ptBinEdges[1][0]-ptBinEdges[0][0])/ptBinEdges[1][1])+nbin11;
  const Int_t nbin13 = round((ptBinEdges[2][0]-ptBinEdges[1][0])/ptBinEdges[2][1])+nbin12;
  const Int_t nbin14 = round((ptBinEdges[3][0]-ptBinEdges[2][0])/ptBinEdges[3][1])+nbin13;

  const Int_t nbin1  = nbin14; //bins in pt 
  const Int_t nbin2  =  36;    //bins in phi
  const Int_t nbin3  =  10;     //bins in eta
  const Int_t nbin4  =  1;//6;     //bins in NClsIter1: 0 70 80 90 100 120


  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin1;
  iBin[1]=nbin2;
  iBin[2]=nbin3;
  iBin[3]=nbin4;
   
  //arrays for lower bounds :
  Double_t *binLim1=new Double_t[nbin1+1];
  Double_t *binLim2=new Double_t[nbin2+1];
  Double_t *binLim3=new Double_t[nbin3+1];
  Double_t *binLim4=new Double_t[nbin4+1];
  
  //values for bin lower bounds 
  for(Int_t i=0; i<=nbin1; i++) {
    if(i<=nbin11) binLim1[i]=(Double_t)ptmin + (ptBinEdges[0][0]-ptmin)/nbin11*(Double_t)i ;  
    if(i<=nbin12 && i>nbin11) binLim1[i]=(Double_t)ptBinEdges[0][0] + (ptBinEdges[1][0]-ptBinEdges[0][0])/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ; 
    if(i<=nbin13 && i>nbin12) binLim1[i]=(Double_t)ptBinEdges[1][0] + (ptBinEdges[2][0]-ptBinEdges[1][0])/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ; 
    if(i<=nbin14 && i>nbin13) binLim1[i]=(Double_t)ptBinEdges[2][0] + (ptBinEdges[3][0]-ptBinEdges[2][0])/(nbin14-nbin13)*((Double_t)i-(Double_t)nbin13) ; 
  }
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)phimin + (phimax-phimin)/nbin2*(Double_t)i ;
  for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)etamin + (etamax-etamin)/nbin3*(Double_t)i ;  
  binLim4[0] = 70.;//0.;
  binLim4[1] = 160.;//70.;
  // binLim4[2] = 80.;
  // binLim4[3] = 90.;
  // binLim4[4] = 100.;
  // binLim4[5] = 120.;
  // binLim4[6] = 160.;


  AliCFContainer* containerPos = new AliCFContainer("containerPos","container for positive tracks",nstep,nvar,iBin);
  //setting the bin limits
  containerPos -> SetBinLimits(ipt,binLim1);
  containerPos -> SetBinLimits(iphi,binLim2);
  containerPos -> SetBinLimits(ieta,binLim3);
  containerPos -> SetBinLimits(incls,binLim4);

  AliCFContainer* containerNeg = new AliCFContainer("containerNeg","container for negative tracks",nstep,nvar,iBin);
  //setting the bin limits
  containerNeg -> SetBinLimits(ipt,binLim1);
  containerNeg -> SetBinLimits(iphi,binLim2);
  containerNeg -> SetBinLimits(ieta,binLim3);
  containerNeg -> SetBinLimits(incls,binLim4);
  
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts, only for ESD analysis
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  AliESDtrackCuts *trackCutsReject = 0x0;

  if(dataType=="ESD") {
    //Standard Cuts
    //Set track cuts for global tracks
    if(trackType==0 && cuts==0) {
      // tight global tracks - RAA analysis
      trackCuts = CreateTrackCutsPWGJE(1000);
    }
    if(trackType==0 && cuts==1) {
      //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis
      trackCuts = CreateTrackCutsPWGJE(10001006);
    }
    if(trackType==0 && cuts==5) {
      //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis + NCrossedRowsCut>120 recommended in 2011
      trackCuts = CreateTrackCutsPWGJE(10001008);
    }
    
    if(trackType==0 && cuts==2) {
      //Cuts global tracks with ITSrefit requirement but without SPD
      if (!strcmp(prodType,"LHC12a15e"))
        trackCuts = CreateTrackCutsPWGJE(10011008);
      else
        trackCuts = CreateTrackCutsPWGJE(10011006);
    }
    if(trackType==7 && cuts==0) {
      // no requirements on SPD and ITSrefit failed
      trackCuts = CreateTrackCutsPWGJE(10041006);   //no ITSrefit requirement
      trackCutsReject = CreateTrackCutsPWGJE(1006); //ITSrefit requirement
      trackCutsReject->SetEtaRange(etamin,etamax);
      trackCutsReject->SetPtRange(0.15, 1e10);
    }
    if(trackType==7 && cuts==4) {
      // tight global tracks +  NCrossedRowsCut>120 recommended in 2011
      trackCuts = CreateTrackCutsPWGJE(10041008);
      trackCutsReject = CreateTrackCutsPWGJE(1008);
      trackCutsReject->SetEtaRange(-0.9,0.9);
      trackCutsReject->SetPtRange(0.15, 1e10);
    }
    if(trackType==7 && cuts==1) {
      //Cuts global tracks with ITSrefit requirement but without SPD
      trackCuts = CreateTrackCutsPWGJE(10011006);
    }
    if(trackType==7 && cuts==5) {
      // tight global tracks  + NCrossedRowsCut>120 recommended in 2011
      trackCuts = CreateTrackCutsPWGJE(10011008);
    }
  
    if(trackType==7 && cuts==2) {
      // no requirements on SPD and ITSrefit failed
      trackCuts = CreateTrackCutsPWGJE(10041006);       //no ITSrefit requirement filter 256
      trackCutsReject = CreateTrackCutsPWGJE(10001006); //ITSrefit requirement filter 16
      trackCutsReject->SetEtaRange(etamin,etamax);
      trackCutsReject->SetPtRange(0.15, 1e10);
    }
    if(trackType==7 && cuts==6) {
      // no requirements on SPD and ITSrefit failed
      trackCuts = CreateTrackCutsPWGJE(10041008);       //no ITSrefit requirement filter 256
      trackCutsReject = CreateTrackCutsPWGJE(10001008); //ITSrefit requirement filter 16
      trackCutsReject->SetEtaRange(-0.9,0.9);
      trackCutsReject->SetPtRange(0.15, 1e10);
    }
  
  
  
    if(trackType==1 && cuts==0) {
      //Set track cuts for TPConly tracks
      trackCuts = CreateTrackCutsPWGJE(2001);
    }
    if(trackType==2 && cuts==0) {
       //	      Set track cuts for TPConly constrained tracks
      trackCuts = CreateTrackCutsPWGJE(2001);
    }
    trackCuts->SetEtaRange(etamin,etamax);
    trackCuts->SetPtRange(0.15, 1e10);
  }

  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange(0.15,1e10);
  mcKineCuts->SetEtaRange(etamin,etamax);//-0.5,0.5);
  mcKineCuts->SetRequireIsCharged(kTRUE);

  //Acceptance Cuts
  //AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts","MC acceptance cuts");
  // mcAccCuts->SetMinNHitITS(mintrackrefsITS);
  //mcAccCuts->SetMinNHitTPC(mintrackrefsTPC);

  TObjArray* recMCList = new TObjArray(0);
  TObjArray* secList = new TObjArray(0) ;

  printf("CREATE MC KINE CUTS\n");
  TObjArray* recList = new TObjArray(0);
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  recList->AddLast(mcKineCuts);
  //mcList->AddLast(mcAccCuts);

  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* manPos = new AliCFManager("manPos","Manager for Positive tracks") ;
  manPos->SetParticleContainer(containerPos);
  manPos->SetParticleCutsList(kStepReconstructed,recList);
  manPos->SetParticleCutsList(kStepSecondaries,secList);
  manPos->SetParticleCutsList(kStepReconstructedMC,recMCList);
  manPos->SetParticleCutsList(kStepMCAcceptance,mcList);

  AliCFManager* manNeg = new AliCFManager("manNeg","Manager for Negative tracks") ;
  manNeg->SetParticleContainer(containerNeg);
  manNeg->SetParticleCutsList(kStepReconstructed,recList);
  manNeg->SetParticleCutsList(kStepSecondaries,secList);
  manNeg->SetParticleCutsList(kStepReconstructedMC,recMCList);
  manNeg->SetParticleCutsList(kStepMCAcceptance,mcList);

  TString trigName = "";
  if (triggerMask == AliVEvent::kAnyINT)
    trigName += "kAnyINT";
  else if (triggerMask == AliVEvent::kAny)
    trigName += "kAny";
  else if(triggerMask == AliVEvent::kINT7)
    trigName += "kINT7";
  else if(triggerMask == AliVEvent::kMB)
    trigName += "kMB";
  else if(triggerMask == AliVEvent::kEMC7)
    trigName += "kEMC7";
  else if(triggerMask == AliVEvent::kEMCEJE)
    trigName += "kEMCEJE";
  else if(triggerMask == AliVEvent::kEMCEGA)
    trigName += "kEMCEGA";


  AliPWG4HighPtSpectra *taskPWG4HighPtSpectra = new AliPWG4HighPtSpectra(Form("AliPWG4HighPtSpectra%s_%s",taskName.Data(),trigName.Data()));
  taskPWG4HighPtSpectra->SetTrackType(trackType);
  if(dataType=="AOD")
    taskPWG4HighPtSpectra->SetFilterMask(filterMask);
  else {
    taskPWG4HighPtSpectra->SetCuts(trackCuts);
    taskPWG4HighPtSpectra->SetCutsReject(trackCutsReject);
  }
  taskPWG4HighPtSpectra->SetCFManagerPos(manPos); //here is set the CF manager +
  taskPWG4HighPtSpectra->SetCFManagerNeg(manNeg); //here is set the CF manager -
  taskPWG4HighPtSpectra->SetTriggerMask(triggerMask);
  taskPWG4HighPtSpectra->SelectHIJINGOnly(bSelHijingParticles);
  taskPWG4HighPtSpectra->SetReadAODData(dataType=="AOD"? kTRUE : kFALSE);
  if(!usePythiaxsec)
    taskPWG4HighPtSpectra->SetNoPythiaInfo();

  if(beamType=="PbPb" || beamType=="pPb") {
    taskPWG4HighPtSpectra->SetIsPbPb(kTRUE);
    taskPWG4HighPtSpectra->SetCentralityClass(centClass);
  }
  //  taskPWG4HighPtSpectra->SetSigmaConstrainedMax(5.);


  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  //------ output containers ------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG4_HighPtSpectra%s_%s",taskName.Data(),trigName.Data());

  AliAnalysisDataContainer *coutput0 = 0x0;
  AliAnalysisDataContainer *coutput1 = 0x0;
  AliAnalysisDataContainer *coutput2 = 0x0;
  AliAnalysisDataContainer *cout_cuts0 = 0x0;

  coutput0 = mgr->CreateContainer(Form("chist0HighPtSpectra%s_%s",taskName.Data(),trigName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  coutput1 = mgr->CreateContainer(Form("ccontainer0HighPtSpectra%s_%s",taskName.Data(),trigName.Data()), AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  coutput2 = mgr->CreateContainer(Form("ccontainer1HighPtSpectra%s_%s",taskName.Data(),trigName.Data()), AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  cout_cuts0 = mgr->CreateContainer(Form("qa_SpectraTrackCuts%s_%s",taskName.Data(),trigName.Data()), AliESDtrackCuts::Class(), AliAnalysisManager::kParamContainer,outputfile);
  
  mgr->AddTask(taskPWG4HighPtSpectra);

  mgr->ConnectInput(taskPWG4HighPtSpectra,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4HighPtSpectra,0,coutput0);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,1,coutput1);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,2,coutput2);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,3,cout_cuts0);

  // Return task pointer at the end
  return taskPWG4HighPtSpectra;
}
