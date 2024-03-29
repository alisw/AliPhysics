/// \file AddTaskPi0IMGammaCorrQA.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Configuration of the analysis QA wagon for PWG-GA EMCal analysis
///
/// Configuration macro of EMCal related PWG-GA  analysis, although it can be
/// also used for PHOS. It does:
/// * a simple photon cluster selection
/// * invariant mass analysis
/// * cluster-charged track correlation analysis (optional)
/// * detector general QA analysis (optional)
/// * track general QA analysis (optional)
///
/// Wagon responsible: Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)


// Set includes for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

// ROOT
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>

// AliPhysics
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVTrack.h"
#include "AliESDtrackCuts.h"

// CaloTrackCorrelations frame, base
#include "AliAnalysisTaskCaloTrackCorrelation.h"
#include "AliAnaCaloTrackCorrMaker.h"
#include "AliCaloTrackESDReader.h"
#include "AliCaloTrackAODReader.h"
#include "AliCalorimeterUtils.h"

// CaloTrackCorrelations frame, analysis
#include "AliAnaPhoton.h"
#include "AliAnaPi0.h"
#include "AliHistogramRanges.h"
#include "AliAnaParticleIsolation.h"
#include "AliAnaParticleHadronCorrelation.h"
#include "AliAnaChargedParticles.h"
#include "AliAnaCalorimeterQA.h"
#include "AliAnaGeneratorKine.h"

// Loaded macros of CaloTrackCorr and other
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#include "ConfigureEMCALRecoUtils.C"
#include "PWGGA/CaloTrackCorrelations/macros/ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C"
#include "PWGGA/CaloTrackCorrelations/macros/CheckActiveEMCalTriggerPerPeriod.C"
#include "PWGGA/CaloTrackCorrelations/macros/GetAlienGlobalProductionVariables.C"
#include "PWGJE/macros/CreateTrackCutsPWGJE.C"


#endif // no CINT


///
/// Configure histograms ranges and bins
///
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges, TString calorimeter, Bool_t caloType,
                            TString collision,               Int_t year)
{
  histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
  
  if(calorimeter=="EMCAL")
  {
    if ( year == 2010 )
    {
      histoRanges->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 120*TMath::DegToRad(), 2*24) ;
      histoRanges->SetHistoXRangeAndNBins(-230,90,120); // QA
      histoRanges->SetHistoYRangeAndNBins(370,450,40);  // QA
    }
    else if ( year < 2014 )
    {           
      histoRanges->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 180*TMath::DegToRad(), 5*24) ;
      histoRanges->SetHistoXRangeAndNBins(-460,90,200); // QA
      histoRanges->SetHistoYRangeAndNBins(100,450,100); // QA    
    }
    else // Run2
    {
      if      (caloType == 0)
        histoRanges->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 187*TMath::DegToRad(),  5*24+8) ;
      else if (caloType == 1) 
        histoRanges->SetHistoPhiRangeAndNBins(260*TMath::DegToRad(), 327*TMath::DegToRad(), 3*24+8) ;
      else                 
        histoRanges->SetHistoPhiRangeAndNBins(80 *TMath::DegToRad(), 327*TMath::DegToRad(), 247) ;
      
      histoRanges->SetHistoXRangeAndNBins(-460,460,230); // QA
      histoRanges->SetHistoYRangeAndNBins(-450,450,225); // QA
    }
    
    histoRanges->SetHistoEtaRangeAndNBins(-0.70, 0.70, 2*48) ;
  }
  else
  {
    histoRanges->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 320*TMath::DegToRad(), 70) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
  }
  
  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 2.9, 150);
  
  // Invariant mass histoRangeslysis
  histoRanges->SetHistoMassRangeAndNBins(0., 0.8, 160) ;
  histoRanges->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  histoRanges->SetHistoOpeningAngleRangeAndNBins(0,0.7,50);
  
  // check if time calibration is on
  histoRanges->SetHistoTimeRangeAndNBins(-250.,250,250);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-150, 150, 150);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.05,0.05,100);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.05,0.05,100);
  histoRanges->SetHistodRRangeAndNBins(0.,0.05,50);//QA
  
  // QA, electron, charged
  histoRanges->SetHistoEOverPRangeAndNBins(0,  2. ,100);
  histoRanges->SetHistodEdxRangeAndNBins  (0.,200.,100);
  
  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,250);
  histoRanges->SetHistoZRangeAndNBins(-350,350,175);
  histoRanges->SetHistoRRangeAndNBins(430,460,30);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,250);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,250);
  
  // QA, correlation
  if ( collision=="PbPb" )
  {
    histoRanges->SetHistoNClusterCellRangeAndNBins(0,100,100);
    histoRanges->SetHistoNClustersRangeAndNBins(0,500,50);
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,2000,200);
  }
  else
  {
    histoRanges->SetHistoNClusterCellRangeAndNBins(0,50,50);
    histoRanges->SetHistoNClustersRangeAndNBins(0,50,50);
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,200,200);
  }
  
  // correlation
  histoRanges->SetHistoDeltaPhiRangeAndNBins(-1.7, 4.8, 163);
  histoRanges->SetHistoDeltaEtaRangeAndNBins(-2.0,2.0,20);
  
  // xE, zT
  histoRanges->SetHistoRatioRangeAndNBins(0.,1.2,120);
  histoRanges->SetHistoHBPRangeAndNBins  (0.,10.,100);
  
  // Isolation
  histoRanges->SetHistoPtInConeRangeAndNBins(0, 25 ,50);
  histoRanges->SetHistoPtSumRangeAndNBins   (0, 50, 50);
  if ( collision.Contains("pPb") )
    histoRanges->SetHistoPtSumRangeAndNBins   (0, 100, 50);
  else if ( collision.Contains("PbPb") )
    histoRanges->SetHistoPtSumRangeAndNBins   (0, 200, 50);
}


///
/// Configure the class handling the events and cluster/tracks filtering.
///
AliCaloTrackReader * ConfigureReader(TString inputDataType, TString collision, Bool_t calibrate,
                                     Int_t   minTime,       Int_t maxTime,
                                     Int_t   minCen,        Int_t maxCen,
                                     Bool_t  simulation,    Int_t year,        Int_t debugLevel)
{
  AliCaloTrackReader * reader = 0;
  if     (inputDataType=="AOD")
    reader = new AliCaloTrackAODReader();
  else if(inputDataType=="ESD")
    reader = new AliCaloTrackESDReader();
  else 
    printf("AddTaskPi0IMGammaCorrQA::ConfigureReader() - Data combination not known input Data=%s\n",
           inputDataType.Data());
  
  reader->SetDebug(debugLevel);//10 for lots of messages
  
  //------------------------
  // Detector input filling
  //------------------------
  
  //Min cluster/track E
  reader->SetEMCALEMin(0.3); 
  reader->SetEMCALEMax(1000); 
  reader->SetPHOSEMin(0.3);
  reader->SetPHOSEMax(1000);
  reader->SetCTSPtMin(0.2);
  reader->SetCTSPtMax(1000);

  // Time cut 
  reader->SwitchOffUseParametrizedTimeCut();
  
  if(calibrate)
  {
    reader->SwitchOnUseEMCALTimeCut() ;
    reader->SetEMCALTimeCut(minTime,maxTime);
  }
  
  reader->SwitchOffUseTrackTimeCut();
  reader->SetTrackTimeCut(-1e10,1e10);

  reader->SwitchOffFiducialCut();
  
  // Tracks
  reader->SwitchOffCTS();
  reader->SwitchOffRejectNoTrackEvents();
  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();
  reader->SwitchOffUseTrackDCACut();
  //reader->SetTrackDCACut(0,0.0105);
  //reader->SetTrackDCACut(1,0.035);
  //reader->SetTrackDCACut(2,1.1);
  
  if(inputDataType=="ESD")
  {
#if defined(__CINT__)
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
#endif

    if(year > 2010)
    {
      //Hybrids 2011
      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
      reader->SetTrackCuts(esdTrackCuts);
      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
      reader->SetTrackComplementaryCuts(esdTrackCuts2);
    }
    else
    {
      //Hybrids 2010
      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001006);
      reader->SetTrackCuts(esdTrackCuts);
      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10041006);
      reader->SetTrackComplementaryCuts(esdTrackCuts2);
    }
  }
  else if(inputDataType=="AOD")
  {
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SetTrackStatus(AliVTrack::kITSrefit);    
  }
  
  // Calorimeter
  
  //reader->SetEMCALClusterListName("");
  
  if(calibrate && !simulation) reader->SwitchOnClusterRecalculation();
  else                         reader->SwitchOffClusterRecalculation();
  
  reader->SwitchOnEMCALCells();
  reader->SwitchOnEMCAL();

  reader->SwitchOffPHOSCells();
  reader->SwitchOffPHOS();
  
  //-----------------
  // Event selection
  //-----------------
  
  reader->SwitchOnEventTriggerAtSE();
  
  reader->SetZvertexCut(10.);
  reader->SwitchOnPrimaryVertexSelection();  // and besides primary vertex
  reader->SwitchOffPileUpEventRejection();   // remove pileup
  reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
  
  if(collision=="PbPb")
  {
    if(year < 2014) reader->SwitchOnAliCentrality();
    reader->SetCentralityBin(minCen,maxCen); // Accept all events, if not select range
    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
  }
  
  if(debugLevel > 0) reader->Print("");
  
  return reader;
}

///
/// Configure the class handling the calorimeter clusters specific methods
///
AliCalorimeterUtils* ConfigureCaloUtils(TString calorimeter, TString trigger, 
                                        Bool_t  simulation , Bool_t  calibrate,
                                        Int_t   year       , Int_t   debugLevel)
{
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(debugLevel);
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  // Search of local maxima in cluster
  cu->SetLocalMaximaCutE(0.1);
  cu->SetLocalMaximaCutEDiff(0.03);

  //cu->SwitchOffClusterPlot();
  cu->SwitchOffRecalculateClusterTrackMatching();
  
  cu->SwitchOnBadChannelsRemoval() ;
  
  //EMCAL settings

  if(!simulation)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  cu->SwitchOffRecalibration(); // Check the reader if it is taken into account during filtering
  cu->SwitchOffRunDepCorrection();

  cu->SwitchOffCorrectClusterLinearity();

  Bool_t bExotic  = kTRUE;
  Int_t  bNonLin  = 0;
  Bool_t bBadMap  = kTRUE;
  
  Bool_t bEnCalib = kFALSE;
  Bool_t bTiCalib = kFALSE;
  
  if(calibrate && !simulation)
  {
    cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
    cu->SwitchOffRunDepCorrection();    
    cu->SwitchOnRecalculateClusterPosition() ;

    bEnCalib = kTRUE;
    bTiCalib = kTRUE;
  }
  
//  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
//
//  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
//  ConfigureEMCALRecoUtils(recou,
//                          simulation,
//                          bExotic,
//                          bNonLin,
//                          bEnCalib,
//                          bBadMap,
//                          bTiCalib,
//                          debugLevel
//                          );
  
  cu->ConfigureEMCALRecoUtils(simulation,
                              bExotic,
                              bNonLin,
                              bEnCalib,
                              bBadMap,
                              bTiCalib,
                              debugLevel);
  
  //recou->SetExoticCellDiffTimeCut(50.);

  if(calorimeter=="PHOS")
  {
    if(year < 2014) cu->SetNumberOfSuperModulesUsed(3);
    else            cu->SetNumberOfSuperModulesUsed(4);
  }
  else
  {
    Int_t nSM     = 20;
    Int_t lastEMC = 11;
    if      (year == 2010) { nSM =  4; lastEMC = 3; }// EMCAL first year
    else if (year <  2014) { nSM = 10; lastEMC = 9; }// EMCAL active 2011-2013
    
    cu->SetNumberOfSuperModulesUsed(nSM);
    
    if      (trigger.Contains("EMCAL"))
    {
      cu->SetFirstSuperModuleUsed( 0);
      cu->SetLastSuperModuleUsed (lastEMC);
    }
    else if (trigger.Contains("DCAL"))
    {
      cu->SetFirstSuperModuleUsed(12);
      cu->SetLastSuperModuleUsed (19);
    }
    else
    {
      cu->SetFirstSuperModuleUsed(0);
      cu->SetLastSuperModuleUsed (cu->GetNumberOfSuperModulesUsed()-1);
    }
    
    printf("AddTaskPi0IMGammaCorrQA::CalorimeterUtils() - nSM %d, first %d, last %d\n",
           cu->GetNumberOfSuperModulesUsed(),cu->GetFirstSuperModuleUsed(), cu->GetLastSuperModuleUsed());
  }

  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  if(debugLevel > 0) cu->Print("");
  
  return cu;
}

///
/// Configure the task doing the first photon cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection ...
///
AliAnaPhoton* ConfigurePhotonAnalysis(TString calorimeter,   Bool_t caloType, TString collision,
                                      TString containerName, Bool_t simulation, 
                                      Int_t year,            Int_t debugLevel)
{
  AliAnaPhoton *ana = new AliAnaPhoton();
  ana->SetDebug(debugLevel); //10 for lots of messages
  
  // cluster selection cuts
  
  ana->SwitchOnFiducialCut(); 
  if(caloType==0)ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.7,  80, 187) ; // EMC 
  else           ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.7, 260, 327) ; // DMC
  ana->GetFiducialCut()->DoEMCALFiducialCut(kTRUE);

  ana->SetCalorimeter(calorimeter);
  
  if(calorimeter == "PHOS")
  {
    ana->SetNCellCut(2);// At least 3 cells
    ana->SetMinPt(0.5);
    ana->SetMinDistanceToBadChannel(2, 4, 5);
    ana->SetTimeCut(-1e10,1e10); // open cut
  }
  else 
  {
    // EMCAL
    ana->SetConstantTimeShift(615); // for MC and uncalibrated data, whenever there is time > 400 ns
    ana->SetNCellCut(1);// At least 2 cells
    ana->SetMinEnergy(0.5); // avoid mip peak at E = 260 MeV
    ana->SetMaxEnergy(1000); 
    ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off 
    // restrict to less than 100 ns when time calibration is on 
    ana->SetMinDistanceToBadChannel(2, 4, 6); 
    // Not useful if M02 cut is already strong
    ana->SetNLMCut(1, 2) ;
  }
  
  ana->SwitchOnTrackMatchRejection() ;
  ana->SwitchOnTMHistoFill() ;
  
  ana->SwitchOnFillOnlyPtHisto();
  
  ana->SwitchOnAcceptanceHistoPerEBin();
  
  //PID cuts (shower shape)
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
  //Not used in bayesian
  
  //EMCAL
  caloPID->SetEMCALLambda0CutMax(0.4); // Rather open
  caloPID->SetEMCALLambda0CutMin(0.10);
  
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
    
  //PHOS
  caloPID->SetPHOSDispersionCut(2.5);
  caloPID->SetPHOSRCut(2.);
  
  ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection
  //if(!simulation)ana->SwitchOnFillPileUpHistograms();

  if(collision.Contains("Pb"))   ana->SwitchOnFillHighMultiplicityHistograms();
  
  // Input / output delta AOD settings
  ana->SetOutputAODName(Form("Photon%s_Calo%d",containerName.Data(),caloType));
  ana->SetOutputAODClassName("AliCaloTrackParticleCorrelation");
  ana->SetInputAODName (Form("Photon%s_Calo%d",containerName.Data(),caloType));

  // Set Histograms name tag, bins and ranges
  ana->AddToHistogramsName(Form("AnaPhoton_Calo%d_",caloType));
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter,caloType,collision,year); // see method below
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms(9);
  ana->FillNPrimaryHistograms(6);
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0 ) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the 2 cluster invariant mass analysis
///
AliAnaPi0* ConfigurePi0Analysis(TString calorimeter  , Bool_t caloType  , TString collision,
                                TString containerName, Bool_t simulation, Int_t   year,
                                Int_t   debugLevel   , Int_t  minCen)
{
  AliAnaPi0 *ana = new AliAnaPi0();
  
  ana->SetDebug(debugLevel);//10 for lots of messages
  
  // Input delta AOD settings
  ana->SetInputAODName(Form("Photon%s_Calo%d",containerName.Data(),caloType));
  
  // Calorimeter settings
  ana->SetCalorimeter(calorimeter);
  
  ana->SwitchOnFiducialCut(); 
  if(caloType==0)ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.7,  80, 187) ; // EMC 
  else           ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.7, 260, 327) ; // DMC
  ana->GetFiducialCut()->DoEMCALFiducialCut(kTRUE);  

  ana->SwitchOnRealCaloAcceptance();
  
  // Settings for pp collision mixing
  ana->SwitchOnOwnMix(); //Off when mixing done with general mixing frame
  
  // Cuts 
  if(calorimeter=="EMCAL") 
  {
    ana->SetPairTimeCut(100);
    
    // Angle cut, avoid pairs with too large angle
    ana->SwitchOnAngleSelection(); 
    ana->SetAngleMaxCut(TMath::DegToRad()*80.); // EMCal: 4 SM in phi, 2 full SMs in eta
    ana->SetAngleCut(0.016); // Minimum angle open, ~cell size
  }
  
  ana->SetNPIDBits(1);
  ana->SetNAsymCuts(1); // no asymmetry cut, previous studies showed small effect.
  // In EMCAL assymetry cut prevents combination of assymetric decays which is the main source of pi0 at high E.
    
  if     (collision == "pp"  )
  {
    ana->SetNCentrBin(1);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(1);
    ana->SetNMaxEvMix(100);
    ana->SetMinPt(0.5);
  }
  else if(collision =="PbPb")
  {
    ana->SetNCentrBin(10);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(4);
    ana->SetNMaxEvMix(10);  
    if(minCen >= 10) ana->SetNMaxEvMix(50); 
    if(minCen >= 50) ana->SetNMaxEvMix(100);
    ana->SetMinPt(1.5);
    ana->SwitchOnFillHighMultiplicityHistograms();
  }
  else if(collision =="pPb")
  {
    ana->SetNCentrBin(1);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(4);
    ana->SetNMaxEvMix(100);
    ana->SetMinPt(0.5);
    ana->SwitchOnFillHighMultiplicityHistograms();
  }

  ana->SwitchOffMultipleCutAnalysis();
  ana->SwitchOnSMCombinations();
  ana->SwitchOffFillAngleHisto();
  ana->SwitchOffFillOriginHisto();
  
  // Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaPi0_Calo%d_",caloType));
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter,caloType,collision,year); // see method below

  if(simulation) ana->SwitchOnDataMC();

  if(debugLevel > 0) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing charged track selection
///
AliAnaChargedParticles* ConfigureChargedAnalysis(TString collision,TString containerName,
                                                 Bool_t simulation, Int_t year, Int_t debugLevel)
{
  AliAnaChargedParticles *ana = new AliAnaChargedParticles();
  ana->SetDebug(debugLevel); //10 for lots of messages
    
  // selection cuts
  
  ana->SetMinPt(0.5);
  ana->SwitchOnFiducialCut();
  Float_t etacut = 0.8;
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(etacut, 0, 360) ; //more restrictive cut in reader and after in isolation

  // histogram switchs
  
  ana->SwitchOffFillVertexBC0Histograms() ;
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  ana->SwitchOffFillTrackMultiplicityHistograms();
  
  // Input / output delta AOD settings
  
  ana->SetOutputAODName(Form("Hadron%s",containerName.Data()));
  ana->SetOutputAODClassName("AliCaloTrackParticleCorrelation");
  ana->SetInputAODName(Form("Hadron%s",containerName.Data()));

  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaHadrons_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),"",kFALSE,collision,year); // see method below
  
  ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 120) ;
  ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.*etacut, 1.*etacut, etacut*100) ;
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the trigger particle hadron correlation
///
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle  , TString calorimeter  , Bool_t caloType,
                                                    TString collision , TString containerName,
                                                    Bool_t  simulation, Int_t   year         , Int_t  debugLevel)
{
  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
  ana->SetDebug(debugLevel);
    
  //if(collision.Contains("Pb"))   ana->SwitchOnFillHighMultiplicityHistograms();
  
  ana->SetMinPt(5);
  
  ana->SwitchOffStudyTracksInCone() ;
  //ana->SwitchOnUEBandSubtractionHistoFill();

  ana->SwitchOffDecayTaggedHistoFill() ;
  ana->SwitchOffSSHistoFill();
  
  ana->SwitchOffLeadingOnly();
  ana->SwitchOffCheckNeutralClustersForLeading();

  ana->SwitchOffTMHistoFill();
  
  // MC
  ana->SwitchOffPrimariesInConeSelection();
  ana->SwitchOffPrimariesPi0DecayStudy() ;
  
  ana->SwitchOnRealCaloAcceptance();
  ana->SwitchOnFiducialCut();
  
  ana->SwitchOnOnlyTH3HistoFill();
  
  if(calorimeter == "EMCAL" && caloType == 0)
  {
    // Avoid borders of EMCal
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.60, 86, 174) ;
  }
  if(calorimeter == "EMCAL" && caloType == 1)
  {
    // Avoid borders of DCal
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.60, 264, 316) ;
  }
 
  AliCaloPID* caloPID = ana->GetCaloPID();
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
  
  //
  // Do settings for main isolation cut class
  //
  AliIsolationCut * ic =  ana->GetIsolationCut();
  ic->SetDebug(debugLevel);
  ic->SetParticleTypeInCone(AliIsolationCut::kOnlyCharged);
  ic->SetICMethod(AliIsolationCut::kSumBkgSubIC);
  if ( collision == "pp" || collision == "pPb" )
  {
    ic->SetPtThreshold(0.5);
    ic->SetSumPtThreshold(2.0) ;
    ic->SetConeSize(0.4);
  }
  if ( collision == "PbPb" )
  {
    ic->SetPtThreshold(3.);
    ic->SetSumPtThreshold(3.0) ;
    ic->SetConeSize(0.3);
  }

  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s_Calo%d",particle.Data(),containerName.Data(),caloType));
  ana->SetAODObjArrayName(Form("%sIso_%s_Calo%d",particle.Data(),containerName.Data(),caloType));
  
  // Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaIsol%s_Calo%d_",particle.Data(),caloType));
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter,caloType,collision,year); // see method below
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the trigger particle hadron correlation
///
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle,  TString calorimeter, Bool_t caloType,
                                                                    TString collision, TString containerName,
                                                                    Bool_t simulation, Int_t year, Int_t debugLevel, Int_t minCen)
{
  AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
  ana->SetDebug(debugLevel);
  
  ana->SetTriggerPtRange(5,100);
  ana->SetAssociatedPtRange(0.2,100);
  //ana->SetDeltaPhiCutRange( TMath::Pi()/2,3*TMath::Pi()/2 ); //[90 deg, 270 deg]
  ana->SetDeltaPhiCutRange  (TMath::DegToRad()*120.,TMath::DegToRad()*240.);
  ana->SetUeDeltaPhiCutRange(TMath::DegToRad()*60. ,TMath::DegToRad()*120.);
  ana->SwitchOffFillEtaGapHistograms();

  ana->SetNAssocPtBins(4);
  ana->SetAssocPtBinLimit(0, 0.5) ;
  ana->SetAssocPtBinLimit(1, 2) ;
  ana->SetAssocPtBinLimit(2, 5) ;
  ana->SetAssocPtBinLimit(3, 10) ;
  ana->SetAssocPtBinLimit(4, 20) ;
  
  ana->SelectIsolated(kFALSE); // do correlation with isolated photons

  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  ana->SwitchOffAbsoluteLeading(); // Select trigger leading particle of all the selected tracks
  ana->SwitchOffNearSideLeading(); // Select trigger leading particle of all the particles at +-90 degrees, default
  
  //ana->SwitchOnLeadHadronSelection();
  //ana->SetLeadHadronPhiCut(TMath::DegToRad()*100., TMath::DegToRad()*260.);
  //ana->SetLeadHadronPtCut(0.5, 100);
  
  // Mixing with own pool
  ana->SwitchOffOwnMix();
  
  ana->SetNZvertBin(20);
  ana->SwitchOffCorrelationVzBin() ;

  //if(collision.Contains("Pb"))   ana->SwitchOnFillHighMultiplicityHistograms();
  
  if(collision=="pp")
  {
    ana->SetNMaxEvMix(100);    
    ana->SwitchOnTrackMultBins();
    ana->SetNTrackMultBin(10); // same as SetNCentrBin(10);
    ana->SetNRPBin(1);
  }
  else 
  {
    ana->SetNMaxEvMix(10);  
    if(minCen >= 10) ana->SetNMaxEvMix(50); 
    if(minCen >= 50) ana->SetNMaxEvMix(100);
    ana->SwitchOffTrackMultBins(); // centrality bins
    ana->SetNCentrBin(10); 
    ana->SetNRPBin(3);
  }
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s_Calo%d",particle.Data(),containerName.Data(),caloType));
  ana->SetAODObjArrayName(Form("%sHadronCorr_%s_Calo%d",particle.Data(),containerName.Data(),caloType));
  
  ana->SwitchOffPi0TriggerDecayCorr();
  ana->SwitchOffDecayTriggerDecayCorr();
  ana->SwitchOffNeutralCorr(); // Do only correlation with TPC
  ana->SwitchOffHMPIDCorrelation();
  ana->SwitchOffFillBradHistograms();
  
  ana->SwitchOnFillXEHistograms() ;
  ana->SwitchOnFillZTHistograms() ;
  ana->SwitchOnFillHBPHistograms() ;
  
  if ( simulation )
  {
    ana->SwitchOffFillZTHistograms() ;
    ana->SwitchOffFillHBPHistograms() ;
  }
  
  // Underlying event
  ana->SwitchOffSeveralUECalculation();
  ana->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("Ana%sHadronCorr_Calo%d_",particle.Data(),caloType));
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter,caloType,collision,year); // see method below
  
  if(simulation) ana->SwitchOnDataMC();
  ana->SetMCGenType(0,5);
  
  if(debugLevel > 0) ana->Print("");
 
  return ana;
}

///
/// Configure the task doing standard calorimeter QA
///
AliAnaCalorimeterQA* ConfigureQAAnalysis(TString calorimeter, TString collision,
                                         Bool_t  simulation , 
                                         Int_t   year,        Int_t   debugLevel)
{
  AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
  ana->SetDebug(debugLevel); //10 for lots of messages
  ana->SetCalorimeter(calorimeter);
  
  //printf("AddTaskPi0IMGammaCorrQA::CofigureQAAnalysis() - calorimeter %s, caloType %d, collision %s, simulation %d, fillCellTime %d, year %d, debugLevel  %d\n",
  //       calorimeter.Data(),caloType,collision.Data(),simulation,fillCellTime,year,debugLevel);
  
  ana->SetTimeCut(-1e10,1e10); // Open time cut
  ana->SetConstantTimeShift(615); // for MC and uncalibrated data, whenever there is time > 400 ns
  ana->SetEMCALCellAmpMin(0.5);
  
  ana->SwitchOffStudyBadClusters() ;
  ana->SwitchOffFillAllTH3Histogram();
  ana->SwitchOffFillAllPositionHistogram();
  ana->SwitchOffFillAllPositionHistogram2();
  ana->SwitchOffStudyBadClusters();
  ana->SwitchOffFillAllPi0Histogram()  ;
  ana->SwitchOffCorrelation();
  ana->SwitchOffFillAllCellAbsIdHistogram();
  ana->SwitchOffFillAllTrackMatchingHistogram();
  ana->SwitchOffFillAllClusterHistogram() ;
  ana->SwitchOffFillAllCellTimeHisto() ;
  ana->SwitchOffFillLowGainHistogram() ;
  if ( !simulation )
  {
    ana->SwitchOnFillAllCellTimeHisto() ;
    ana->SwitchOnFillLowGainHistogram() ;
  } 
  ana->SwitchOnFillAllCellHistogram();

  ana->AddToHistogramsName("QA_Cell_"); //Begining of histograms name
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter, -1, collision,year); // see method below
  
//  ana->SwitchOnFiducialCut(); 
//  if(caloType==0)ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.7,  80, 187) ; // EMC 
//  else           ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.7, 260, 327) ; // DMC
//  
//  ana->GetFiducialCut()->DoEMCALFiducialCut(kTRUE);
  
  //if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0) ana->Print("");
  
  return ana;
}


///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle.
/// \param simulation : A bool identifying the data as simulation.
/// \param collision: A string with the colliding system.
/// \param period : A string with the data period: LHC11h, LHC15n ... from it we extract the year.
/// \param qaan: execute the detector QA analysis.
/// \param hadronan: execute the track QA (>0) and cluster-track correlation analysis (>1).
/// \param calibrate: if OADB was updated with calibration parameters not used in reconstruction, apply them here.
/// \param minTime: minimum time cut, leave it open by default even if calibration available, ns
/// \param maxTime: maximum time cut, leave it open by default even if calibration available, ns
/// \param minCen : An int to select the minimum centrality, -1 means no selection.
/// \param maxCen : An int to select the maximum centrality, -1 means no selection.
/// \param debugLevel : An int to define the debug level of all the tasks.
/// \param suffix : A string with the type of trigger (default: MB, EMC).
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskPi0IMGammaCorrQA(const TString  calorimeter   = "EMCAL",
                                                             Bool_t   simulation    = kFALSE,
                                                             TString  collision     = "pp",
                                                             TString  period        = "",
                                                             const Bool_t   qaan          = kTRUE,
                                                             const Int_t    hadronan      = 2,
                                                             const Bool_t   calibrate     = kFALSE,
                                                             const Int_t    minTime       = -1000,
                                                             const Int_t    maxTime       =  1000,
                                                             const Int_t    minCen        = -1,
                                                             const Int_t    maxCen        = -1,
                                                             const Int_t    debugLevel    = -1,
                                                             const char *   suffix        = "default"
                                                             )
{
  printf("AddTaskPi0IMGammaCorrQA::Start configuration\n");

#if defined(__CINT__)
  printf("AddTaskPi0IMGammaCorrQA::Load macros\n");
  // Load macros
  //
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/GetAlienGlobalProductionVariables.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/CheckActiveEMCalTriggerPerPeriod.C");
#endif
  
  // Check the global variables, and reset the provided ones if empty.
  //S
  TString trigger  = suffix;
  Int_t   year        = -1;
  Bool_t  printGlobal = kFALSE;
  if ( trigger.Contains("default") || trigger.Contains("INT") || trigger.Contains("MB") ) printGlobal = kTRUE;
  
  GetAlienGlobalProductionVariables(simulation,collision,period,year,printGlobal);
  
  // Get the pointer to the existing analysis manager via the static access method.
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTaskPi0IMGammaCorrQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //
  if (!mgr->GetInputEventHandler()) 
  {
    ::Error("AddTaskPi0IMGammaCorrQA", "This task requires an input event handler");
    return NULL;
  }
  
  //
  // Create task
  //
  
  // Name for containers
  TString containerName = Form("%s_Trig_%s",calorimeter.Data(), trigger.Data());
  
  if(collision!="pp" && maxCen>=0) containerName+=Form("Cen%d_%d",minCen,maxCen);
  
  TString taskName =Form("Pi0IM_GammaTrackCorr_%s",containerName.Data());
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (taskName);
  //task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(debugLevel);
  //task->SetBranches("ESD:AliESDRun.,AliESDHeader");
  //task->SetBranches("AOD:header,tracks,vertices,emcalCells,caloClusters");
  
  //
  // Init main analysis maker and pass it to the task
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  task->SetAnalysisMaker(maker);
  
  //
  // Pass the task to the analysis manager
  mgr->AddTask(task);
  
  //
  // Create containers
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(trigger, TList::Class(),
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s:%s",outputfile.Data(),Form("Pi0IM_GammaTrackCorr_%s",calorimeter.Data())));
  
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",trigger.Data()), TList::Class(),
                                                             AliAnalysisManager::kParamContainer, 
                                                             Form("%s_Parameters.root",Form("Pi0IM_GammaTrackCorr_%s",calorimeter.Data())));
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  //==============================================================================
  
  // Do not configure the wagon for certain analysis combinations
  // But create the task so that the sub-wagon train can run
  //
  Bool_t doAnalysis = CheckActiveEMCalTriggerPerPeriod(simulation,trigger,period,year);
  if(!doAnalysis) 
  {
    maker->SwitchOffProcessEvent();
    return task;
  }
  
  // #### Start analysis configuration ####
  //  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Make sure the B field is enabled for track selection, some cuts need it
  //
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);
  
  // Print settings to check all is as expected
  //
  printf("AddTaskPi0IMGammaCorrQA - Task NAME: %s \n",taskName.Data());
  
  printf("AddTaskPi0IMGammaCorrQA - Settings: data <%s>, calo <%s>, MC <%d>, collision <%s>, trigger <%s>, period <%s>, year <%d>,\n"
         "\t \t \t  CaloQA on <%d>, Track QA on <%d>, Make corrections <%d>, %d < time < %d, %d < cen < %d, debug level <%d> \n", 
         inputDataType.Data(), calorimeter.Data(),simulation, collision.Data(),trigger.Data(), period.Data(), year,
         qaan , hadronan, calibrate, minTime, maxTime, minCen, maxCen, debugLevel);
  //
  
  // General frame setting and configuration
  maker->SetReader   ( ConfigureReader   (inputDataType,collision,calibrate,minTime,maxTime,minCen,maxCen,simulation,year,debugLevel) );
  if(hadronan)maker->GetReader()->SwitchOnCTS();
  
  maker->SetCaloUtils( ConfigureCaloUtils(calorimeter,trigger,simulation,calibrate,year,debugLevel) );
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important
  
  // Cell QA
  if ( qaan )  maker->AddAnalysis(ConfigureQAAnalysis(calorimeter,collision,
                                                      simulation,
                                                      year,debugLevel),n++);
  
  // Analysis with EMCal trigger or MB
  if ( !trigger.Contains("DCAL") )
  {
    // Cluster selection
    maker->AddAnalysis(ConfigurePhotonAnalysis(calorimeter,0,collision,containerName,simulation,year,debugLevel)       ,n++); 
    // Previous cluster invariant mass
    maker->AddAnalysis(ConfigurePi0Analysis   (calorimeter,0,collision,containerName,simulation,year,debugLevel,minCen),n++);     
    if ( hadronan > 1 )
    {
      // Isolation of selected clusters by AliAnaPhoton
      maker->AddAnalysis(ConfigureIsolationAnalysis("Photon",calorimeter,0,collision,containerName,simulation,year,debugLevel), n++);
      // Selected clusters-track correlation
      maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",calorimeter,0,collision,containerName,simulation,year,debugLevel,minCen), n++); 
    }
  }
  
  // Analysis with DCal trigger or MB
  if ( year > 2014 && calorimeter=="EMCAL" && !trigger.Contains("EMCAL") )
  {
    // Cluster selection
    maker->AddAnalysis(ConfigurePhotonAnalysis(calorimeter,1,collision,containerName,simulation,year,debugLevel)       , n++); 
    // Previous cluster invariant mass
    maker->AddAnalysis(ConfigurePi0Analysis   (calorimeter,1,collision,containerName,simulation,year,debugLevel,minCen),n++); 
    if ( hadronan > 1 )
    {
      // Isolation of selected clusters by AliAnaPhoton
      maker->AddAnalysis(ConfigureIsolationAnalysis("Photon",calorimeter,1,collision,containerName,simulation,year,debugLevel), n++);
      // Selected clusters-track correlation
      maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",calorimeter,1,collision,containerName,simulation,year,debugLevel,minCen), n++); 
    }
  }
  
  // Charged tracks plots, any trigger
  if ( hadronan )
    maker->AddAnalysis(ConfigureChargedAnalysis(collision,containerName,simulation,year,debugLevel), n++); 
  
  if ( simulation )
  {
    // Calculate the cross section weights, apply them to all histograms 
    // and fill xsec and trial histo. Sumw2 must be activated.
    //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionCalculation(); 
    //maker->SwitchOnSumw2Histograms();
    
    // For recent productions where the cross sections and trials are not stored in separate file
    //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionFromEventHeader() ;
    
    // Just fill cross section and trials histograms.
    maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionHistoFill(); 
    
    // Add control histogram with pT hard to control aplication of weights 
    maker->SwitchOnPtHardHistogram();
  }
  
  //
  // Select events trigger depending on trigger
  //
  if ( !simulation )
  {    
    TString caloTriggerString = "";
    UInt_t mask = ConfigureAndGetEventTriggerMaskAndCaloTriggerString(trigger, year, caloTriggerString);
    
    task ->SelectCollisionCandidates( mask );
    maker->GetReader()->SetFiredTriggerClassName(caloTriggerString);
    printf("AddTaskPi0IMGammaCorrQA - Trigger Mask %d, caloTriggerString <%s>\n", mask, caloTriggerString.Data());
  }
  
  //
  // Final maker settings
  //
  maker->SetAnaDebug(debugLevel)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker() ;
  maker->SwitchOnDataControlHistograms(); 
  if ( simulation )
    maker->SwitchOffDataControlHistograms(); 

  if(debugLevel > 0) maker->Print("");
  
  printf("AddTaskPi0IMGammaCorrQA::End configuration\n");
  
  return task;
}



