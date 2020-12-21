  ///////////////////////////////////////////////////////////////////////////
  ///\file AddTaskEMCALPhotonIsolation.C
  ///\brief Configuration of AliAnalysisTaskEMCALPhotonIsolation
  ///
  /// Version to be used in lego train for testing on pp@7TeV
  ///
  /// \author Lucile Ronflette <lucile.ronflette@cern.ch>, SUBATECH, Nantes
  /// \author Davide Francesco Lodato <davide.francesco.lodato@cern.ch>, Utrecht University
  /// \author Marco Marquard <marco.marquard@cern.ch>, University Frankfurt am Main
  /// \author Erwann Masson <Erwann.Masson@subatech.in2p3.fr>, SUBATECH, Nantes
  ///////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALPhotonIsolation.h"

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)

#include "PWGGA/EMCALTasks/macros/config_PhotonIsolation.C"

#endif


AliAnalysisTaskEMCALPhotonIsolation* AddTaskEMCALPhotonIsolation(
                                                                 const char*            periodstr                 = "LHC11c",
                                                                 const char*            ntracks                   = "EmcalTracks",
                                                                 const char*            nclusters                 = "EmcCaloClusters",
                                                                 const UInt_t           pSel                      = AliVEvent::kEMC7,
                                                                 const TString          dType                     = "ESD",
                                                                 const Bool_t	        bHisto                    = kTRUE,
                                                                 const Int_t	        iOutput	  	          = 0,
                                                                 const Bool_t	        bIsMC  	                  = kFALSE,
                                                                 const Bool_t           bMCNormalization          = kFALSE,
                                                                 const Bool_t           bNLMCut                   = kFALSE,
                                                                 const Int_t            NLMCut                    = 0,
                                                                 const Double_t         EtIso                     = 2.,
                                                                 const Int_t            iIsoMethod                = 1,
                                                                 const Int_t            iEtIsoMethod              = 0,
                                                                 const Int_t            iUEMethod                 = 1,
                                                                 const Bool_t           bUseofTPC                 = kFALSE,
                                                                 const Double_t         TMdeta                    = 0.02,
                                                                 const Double_t         TMdphi                    = 0.03,
                                                                 const Bool_t           bTMClusterRejection       = kTRUE,
                                                                 const Float_t          iIsoConeRadius            = 0.4,
                                                                 const Bool_t           iSmearingSS               = kFALSE,
                                                                 const Float_t          iWidthSSsmear             = 0.,
                                                                 const Float_t          iMean_SSsmear             = 0.,
                                                                 const Bool_t           iExtraIsoCuts             = kFALSE,
								 const Bool_t		SetListNameOutput 	  = kFALSE, //add the output type to the EmcalList name
                                                                 TString                configBasePath            = "",
                                                                 const Int_t            bWhichToSmear             = 0,
                                                                 const Int_t            minNLM                    = 1,
                                                                 const Double_t         TMdetaIso                 = 0.02,
                                                                 const Double_t         TMdphiIso                 = 0.03,
                                                                 const Bool_t           bmcTruth                  = kTRUE,
                                                                 TString                triggerName               = "",
                                                                 const Bool_t           RejectPileUpEvent         = kFALSE,
                                                                 const Int_t            NContrToPileUp            = 3,
                                                                 const Float_t          iFiducialCut              = 0.4,
                                                                 const Bool_t           bANwithNoSameTcard        = kFALSE
                                                                 )
{
  Printf("Preparing neutral cluster analysis\n");
  
  // #### Define manager and data container names
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if(!manager){
    ::Error("AddTaskEMCALPhotonIsolation", "No analysis manager to connect to.");
    return NULL;
  }

  // Use the input containers naming convention with "usedefault" (already used in several EMCal tasks and new correction framework)
  TString trackName(ntracks);
  TString clusName(nclusters);
  TString clusInfix("");

  if(trackName == "usedefault"){
    if(dType == "ESD"){
      trackName = "Tracks";
    }
    else if(dType == "AOD"){
      trackName = "tracks";
    }
    else{
      trackName = "";
    }
  }

  if(clusName == "usedefault"){
    if(dType == "ESD"){
      clusName = "CaloClusters";
    }
    else if(dType == "AOD"){
      clusName = "caloClusters";
    }
    else{
      clusName = "";
    }
  }
  else{
      clusInfix = Form("_%s",clusName.Data());
  }

  // Set the task output container name
  Printf("Creating container name for cluster analysis\n");

  TString myContName("");
  if(bIsMC){
    myContName = Form("Analysis_Neutrals_MC%s",clusInfix.Data());
  }
  else{
    myContName = Form("Analysis_Neutrals%s",clusInfix.Data());
  }

  if(SetListNameOutput){
	  myContName.Append(Form("_output%d",iOutput));
  }
  
  // For the 2012 EGA/L1 analysis, only events with EGA/L1 recalc patches are considered
  // (This configuration requires a TriggerMaker wagon in the train!!!)
  TString period     = periodstr;
  Bool_t  is2012_EGA = kFALSE;
  if((triggerName.Contains("EGArecalc") || triggerName.Contains("L1recalc")) && period.Contains("12"))
    is2012_EGA = kTRUE;
  
  if(triggerName.Contains("EG1") || triggerName.Contains("EGA1")){
    triggerName = "_Trigger_EG1";
  }
  else if(triggerName.Contains("EG2") || triggerName.Contains("EGA2")){
    triggerName = "_Trigger_EG2";
  }
  else if(triggerName.Contains("MB")){
    triggerName = "_Trigger_MB";
  }
  else{
    triggerName = "";
  }

  TString pileUp;
  if(RejectPileUpEvent){
    pileUp = Form("_PU_ON%d", NContrToPileUp);
  }
  else{
    pileUp = "";
  }
  
  myContName.Append(Form("%s_TM_%s_CPVe%.2lf_CPVp%.2lf_IsoMet%d_EtIsoMet%d_UEMet%d_TPCbound_%s_IsoConeR%.1f_NLMCut_%s_minNLM%d_maxNLM%d_SSsmear_%s_Width%.3f_Mean_%.3f_PureIso_%s_WhichSmear_%d%s%s%s",is2012_EGA ? "_L1recalc":"",bTMClusterRejection? "On" :"Off", TMdeta , TMdphi ,iIsoMethod,iEtIsoMethod,iUEMethod,bUseofTPC ? "Yes" : "No",iIsoConeRadius,bNLMCut ? "On": "Off",minNLM, NLMCut, iSmearingSS ? "On":"Off",iWidthSSsmear,iMean_SSsmear,iExtraIsoCuts?"On":"Off",bWhichToSmear,triggerName.Data(),pileUp.Data(),bANwithNoSameTcard?"_NoSameTcard":""));

  // Switch on/off the p-Pb analysis features
  Bool_t i_pPb = kFALSE;
  if (period.Contains("12g") ||
      period.Contains("13b") || period.Contains("13c") || period.Contains("13d") || period.Contains("13e") || period.Contains("13f") ||
      period.Contains("16q") || period.Contains("16r") || period.Contains("16s") || period.Contains("16t"))
    i_pPb = kTRUE;

  // #### Define analysis task
  AliAnalysisTaskEMCALPhotonIsolation* task = new AliAnalysisTaskEMCALPhotonIsolation("Analysis",bHisto);
  
#if defined(__CINT__)
  TString configFile("config_PhotonIsolation.C"); // Name of config file
  //  if(gSystem->AccessPathName(configFile.Data())){ // Check for exsisting file and delete it
  //    gSystem->Exec(Form("rm %s",configFile.Data()));
  //  }

//  if(configBasePath.IsNull()){ // Check if a specific config should be used and copy appropriate file
 configBasePath="$ALICE_PHYSICS/PWGGA/EMCALTasks/macros";
 gSystem->Exec(Form("cp %s/%s .",configBasePath.Data(),configFile.Data()));
//  }
//  else if(configBasePath.Contains("alien:///")){
//    gSystem->Exec(Form("alien_cp %s/%s .",configBasePath.Data(),configFile.Data()));
//  }
//  else{
//    gSystem->Exec(Form("cp %s/%s .",configBasePath.Data(),configFile.Data()));
//  }

  configBasePath=Form("%s/",gSystem->pwd());
  TMD5* MD5calc = new TMD5();
  TMD5* MD5sum = MD5calc->FileChecksum(configFile.Data());
  TString configMD5 = MD5sum->AsString();
  configMD5.Resize(8); // Short hash value for usable extension
  TString configFileMD5 = configFile;
  TDatime time; // Get timestamp
  Int_t timeStamp = time.GetTime();
  configFileMD5.ReplaceAll(".C",Form("_%s_%i.C",configMD5.Data(),timeStamp));

  if(gSystem->AccessPathName(configFileMD5.Data())){ // Increase timeStamp by 1 if file exists
    gSystem->Exec(Form("mv %s %s",configFile.Data(),configFileMD5.Data()));
  }
  else{
    while(!gSystem->AccessPathName(configFileMD5.Data())){
      timeStamp++;
      configFileMD5 =configFile;
      configFileMD5.ReplaceAll(".C",Form("_%s_%i.C",configMD5.Data(),timeStamp));
    }
    gSystem->Exec(Form("mv %s %s",configFile.Data(),configFileMD5.Data()));
  }

  TString configFilePath(configBasePath+"/"+configFileMD5);
  gROOT->LoadMacro(configFilePath.Data());

  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/EMCALTasks/macros/config_PhotonIsolation.C");
#endif
  config_PhotonIsolation();
//  Printf("Path of config file: %s\n",configFilePath.Data());

  // Set histrogram bins and ranges
  task->GetHistogramRangesAndBinning()->SetHistoEtaRangeAndNBins(-0.7, 0.7, 96);                                                 // 2 SM = 96 cells in eta
  if(period.Contains("10")){
    task->GetHistogramRangesAndBinning()->SetHistoPhiRangeAndNBins(80.*TMath::DegToRad(), 120.*TMath::DegToRad(), 48);           // 2*2 SM = 48 cells in phi
  }
  else if(period.Contains("11") || period.Contains("12") || period.Contains("13")){
    task->GetHistogramRangesAndBinning()->SetHistoPhiRangeAndNBins(80.*TMath::DegToRad(), 180.*TMath::DegToRad(), 120);          // 2*5 SM = 120 cells in phi
  }
  else{
    task->GetHistogramRangesAndBinning()->SetHistoPhiRangeAndNBins(80.*TMath::DegToRad(), (340.-40./3.)*TMath::DegToRad(), 296); // 2*10 SM = 296 cells in phi
  }

  // #### Task preferences
  task->SelectCollisionCandidates(pSel);
  task->SetOutputFormat(iOutput);
  task->SetIsoConeRadius(iIsoConeRadius);
  task->SetEtIsoThreshold(EtIso);
  task->SetTMClusterRejection(bTMClusterRejection);
  task->SetCTMdeltaEta(TMdeta);
  task->SetCTMdeltaPhi(TMdphi);
  task->SetCTMdeltaEtaIso(TMdetaIso);
  task->SetCTMdeltaPhiIso(TMdphiIso);
  task->SetIsoMethod(iIsoMethod);
  task->SetEtIsoMethod(iEtIsoMethod);
  task->SetUEMethod(iUEMethod);
  task->SetUSEofTPC(bUseofTPC);
  task->SetMC(bIsMC);
  task->SetM02Smearing(iSmearingSS);
  task->SetWidth4Smear(iWidthSSsmear);
  task->SetMean4Smear(iMean_SSsmear);
  task->SetExtraIsoCuts(iExtraIsoCuts);
  task->SetAnalysispPb(i_pPb);
  task->SetNLMCut(bNLMCut,NLMCut,minNLM);
  task->SetPtBinning(ptBin);
  task->SetM02Binning(M02Bin);
  task->SetEtisoBinning(EtisoBin);
  task->SetEtueBinning(EtueBin);
  task->SetEtaBinning(EtaBin);
  task->SetPhiBinning(PhiBin);
  task->SetLabelBinning(LabelBin);
  task->SetPDGBinning(PDGBin);
  task->SetMomPDGBinning(MomPDGBin);
  task->SetClustPDGBinning(ClustPDGBin);
  task->SetDxBinning(DxBin);
  task->SetDzBinning(DzBin);
  task->SetDecayBinning(DecayBin);
  task->SetSmearForClusters(bWhichToSmear);
  task->SetNeedEmcalGeom(kTRUE);
  task->SetMCtruth(bmcTruth);
  task->SetPeriod(periodstr);
  task->SetRejectPileUpEvent(RejectPileUpEvent);
  task->SetNcontributorsToPileUp(NContrToPileUp);
  task->SetFiducialCut(iFiducialCut);
  task->Set2012L1Analysis(is2012_EGA);
  task->SetANWithNoSameTcard(bANwithNoSameTcard);

  if(triggerName.Contains("EG1") || triggerName.Contains("EGA1")){
    if(bIsMC)
      task->SetTriggerLevel1(1);
    else
      task->SetTrigClass("EMC7EG1-B-NOPF-CENTNOTRD");
  }
  else if(triggerName.Contains("EG2") || triggerName.Contains("EGA2")){
    if(bIsMC)
      task->SetTriggerLevel1(2);
    else
      task->SetTrigClass("EMC7EG2-B-NOPF-CENTNOTRD");
  }

  if(bIsMC && bMCNormalization)
    task->SetIsPythia(kTRUE);
  
  TString name(Form("PhotonIsolation_%s_%s", trackName.Data(), clusName.Data()));
  cout<<"Name of the container "<<name.Data()<<endl;
  
  // Tracks to be used for the track matching (already used in TM task, TPC only tracks)
  AliTrackContainer *trackCont = task->AddTrackContainer(trackName);
  if(!trackCont) Printf("Error with TPC-Only Tracks!!");
  trackCont->SetName("tpconlyMatch");
  trackCont->SetTrackFilterType(AliEmcalTrackSelection::kTPCOnlyTracks);
  
  // Tracks to be used in the analysis (Hybrid Tracks)
  AliTrackContainer * tracksForAnalysis = task->AddTrackContainer(trackName);
  if(!tracksForAnalysis) Printf("Error with Hybrids Tracks!!");
  tracksForAnalysis->SetName("filterTracksAna");
  tracksForAnalysis->SetFilterHybridTracks(kTRUE);
  if(!bIsMC){
    tracksForAnalysis->SetTrackCutsPeriod(periodstr);
    tracksForAnalysis->SetDefTrackCutsPeriod(periodstr);
  }

  // Clusters to be used in the analysis (already filtered)
  AliClusterContainer *clusterCont = task->AddClusterContainer(clusName);

  Printf("Name of track container for matching: %s \nName of track container for isolation: %s \nName of cluster container: %s",trackCont->GetName(),tracksForAnalysis->GetName(),clusterCont->GetName());
  
  // #### Add analysis task
  Printf("Task for neutral cluster analysis created and configured, pass it to AnalysisManager\n");
  manager->AddTask(task);
  
  AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:NeutralClusters",AliAnalysisManager::GetCommonFileName()));

  AliAnalysisDataContainer *cinput = manager->GetCommonInputContainer();
  manager->ConnectInput(task, 0, cinput);
  manager->ConnectOutput(task, 1, contHistos);
  
  return task;
}
