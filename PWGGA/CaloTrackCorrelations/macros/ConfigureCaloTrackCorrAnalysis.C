/// \file ConfigureCaloTrackCorrAnalysis.C
/// \ingroup CaloTrackCorrMacros
/// \brief Configuration of gamma-hadron and pi0-hadron + isolation, correlation analysis, with option on which to activate
///
/// Configuration macro for analysis of gamma-hadron and pi0-hadron correlation analysis
/// both pi0 and gamma isolated or not. Extracted from AddTaskGammaHadronCorrelationSelectAnalysis.C. 
/// The base configuration handled in the separate macro AddTaskCaloTrackCorrBase.C and this configuration
/// should be executed like in AddTaskMultipleTrackCutIsoConeAnalysis.C
///
/// It can do:
///   * Photon selection in calorimeter with AliAnaPhoton: Track matching, mild shower shape, NLM ... cuts
///   * Tagging of photon candidate clusters as decay from pi0 or eta or from their side bands with AliAnaPi0EbE (4 calls, one pi0, one eta, one pi0 side band and one eta side band)
///   * Tagging of photon candidate as isolated with AliAnaParticleIsolation
///   * Correlation of the photon and charged tracks, twice, one without isolation condition, and other with isolation condition
///   * Pi0 selection with the identification of merged clusters with AliAnaPi0EbE
///   * Tagging of pi0 as isolated with AliAnaParticleIsolation
///   * Correlation of the pi0 and charged tracks, twice, one without isolation condition, and other with isolation condition
///   * Correlation of photon and reconstructed charged or full jet
///   * Particle correlation and isolation at generator level with AliAnaGeneratorKine
///   * Optionally, the QA tasks AliAnaCalorimeterQA, AliAnaChargedParticle and AliAnaClusterShapeStudies are executed
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

// Set includes for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>

#include "AliAnaPhoton.h"
#include "AliAnaElectron.h"
#include "AliAnaPi0.h"
#include "AliAnaPi0EbE.h"
#include "AliAnaRandomTrigger.h"
#include "AliAnaParticleIsolation.h"
#include "AliAnaParticleHadronCorrelation.h"
#include "AliAnaChargedParticles.h"
#include "AliAnaCalorimeterQA.h"
#include "AliAnaCaloExotics.h"
#include "AliAnaClusterShapeCorrelStudies.h"
#include "AliAnaGeneratorKine.h"
#include "AliAnalysisTaskCaloTrackCorrelation.h"
#include "AliAnaCaloTrackCorrMaker.h"
#include "AliAnaParticleJetFinderCorrelation.h"
#include "AliHistogramRanges.h"

#endif

/// Global name to be composed of the settings, used to set the AOD branch name
TString kAnaCaloTrackCorr = "";

/// Global name to be composed of the analysis components chain and some internal settings
// Some examples of strings: "Photon_MergedPi0_DecayPi0_Isolation_FixIsoConeExcess_Correlation_Bkg_QA_Charged_HighMult_MultiIso_PerSM_PerTCard",
TString kAnaCutsString = ""; 

///
/// Set common mixing/centrality analysis binning
///
/// \param ana : Analysis task where histograms are created and common settings are needed
///
void SetAnalysisMixingCentralityBins(AliAnaCaloTrackCorrBaseClass* ana, TString col)
{
  ana->SwitchOffTrackMultBins();
  //ana->SetNTrackMultBin(10);
  
  ana->SetNZvertBin(10);
  
  if     (col == "pp"  )
  {
    ana->SetNCentrBin(1);
    ana->SetNRPBin(1);
    ana->SetNMaxEvMix(100);
  }
  else if(col == "PbPb")
  {
    ana->SetNCentrBin(10);
    ana->SetNRPBin(4);
    ana->SetNMaxEvMix(10);
    //    if(kAnaCaloTrackCorr.Contains("60_90"))
    //    {
    //      ana->SetNMaxEvMix(50);
    //      ana->SetNCentrBin(2);
    //    }
  }
  else if(col =="pPb")
  {
    ana->SetNCentrBin(1);
    ana->SetNRPBin(4);
    ana->SetNMaxEvMix(100);
  }
}

///
/// Set common detector fiducial cuts for calorimeter acceptance
///
/// \param fidCut  Fiducial cut class pointer
/// \param calorimeter  A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year  The year the data was taken, used to configure some histograms
///
void CalorimeterFiducialCut(AliFiducialCut* fidCut, TString calorimeter, Int_t year)
{
  if      ( calorimeter.Contains("CAL" ) ) fidCut->DoEMCALFiducialCut(kTRUE);  
  else if ( calorimeter.Contains("PHOS") ) fidCut->DoPHOSFiducialCut (kTRUE);  
  
  if      ( calorimeter == "EMCAL" )
  {
    if      ( year > 2014 ) fidCut->SetSimpleEMCALFiducialCut(0.67,  81.2, 185.8) ; //12 SM
    else if ( year > 2010 ) fidCut->SetSimpleEMCALFiducialCut(0.67,  81.2, 178.8) ; //10 SM
    else                    fidCut->SetSimpleEMCALFiducialCut(0.67,  81.2, 118.8) ; // 4 SM
  }
  else if ( calorimeter == "DCAL"  ) 
  {
//  if ( year == 2016 )     fidCut->SetSimpleEMCALFiducialCut(0.67, 261.2, 318.8) ;
//  else                    fidCut->SetSimpleEMCALFiducialCut(0.67, 261.2, 325.8) ;
    
    if ( year == 2016 ) // DCal 2 regions
    {
      Float_t etamax[] = { 0.67,-0.25 };
      fidCut->AddEMCALFidCutMaxEtaArray(2,etamax);
      Float_t etamin[] = { 0.25,-0.67 };
      fidCut->AddEMCALFidCutMinEtaArray(2,etamin);
      
      Float_t phimax[] = { 318.8, 318.8 };
      fidCut->AddEMCALFidCutMaxPhiArray(2,phimax);
      Float_t phimin[] = { 261.2, 261.2 };
      fidCut->AddEMCALFidCutMinPhiArray(2,phimin);
    }
    else // DCal 3 regions
    {
      Float_t etamax[] = { 0.67,-0.25, 0.67 };
      fidCut->AddEMCALFidCutMaxEtaArray(3,etamax);
      Float_t etamin[] = { 0.25,-0.67,-0.67 };
      fidCut->AddEMCALFidCutMinEtaArray(3,etamin);
      
      Float_t phimax[] = { 318.8, 318.8, 325.8 };
      fidCut->AddEMCALFidCutMaxPhiArray(3,phimax);
      Float_t phimin[] = { 261.2, 261.2, 321.2 };
      fidCut->AddEMCALFidCutMinPhiArray(3,phimin);
    }
  }
  else if ( calorimeter == "PHOS"  )
  {
    if      ( year > 2014 ) fidCut->SetSimplePHOSFiducialCut (0.12, 250.5, 319.5) ; 
    else                    fidCut->SetSimplePHOSFiducialCut (0.12, 260.5, 319.5) ;
  }
  
  // Acceptance of EMCal+DCal
  if ( kAnaCaloTrackCorr.Contains("FullCalo") && year > 2014) 
  {
    //fidCut->SetSimpleEMCALFiducialCut(0.67,  81.2, 325.8) ; //20 SM
    
    // EMCal, 1 region - DCal 2 regions, remove 1/3 DCal
    if ( year == 2016 )
    {
      Float_t etamax[] = { 0.67, 0.67,-0.25 };
      fidCut->AddEMCALFidCutMaxEtaArray(3,etamax);
      Float_t etamin[] = {-0.67, 0.25,-0.67 };
      fidCut->AddEMCALFidCutMinEtaArray(3,etamin);
      
      Float_t phimax[] = { 185.8, 318.8, 318.8 };
      fidCut->AddEMCALFidCutMaxPhiArray(3,phimax);
      Float_t phimin[] = {  81.2, 261.2, 261.2 };
      fidCut->AddEMCALFidCutMinPhiArray(3,phimin);
    }
    else // EMCal, 1 region - DCal 3 regions
    {
      Float_t etamax[] = { 0.67, 0.67,-0.25, 0.67 };
      fidCut->AddEMCALFidCutMaxEtaArray(4,etamax);
      Float_t etamin[] = {-0.67, 0.25,-0.67,-0.67 };
      fidCut->AddEMCALFidCutMinEtaArray(4,etamin);
      
      Float_t phimax[] = { 185.8, 318.8, 318.8, 325.8 };
      fidCut->AddEMCALFidCutMaxPhiArray(4,phimax);
      Float_t phimin[] = {  81.2, 261.2, 261.2, 321.2 };
      fidCut->AddEMCALFidCutMinPhiArray(4,phimin);
    }
  }
  
  //  if ( kAnaCutsString.Contains("TightAcc") )
  //  {
  //    if ( calorimeter == "EMCAL" ) fidCut->SetSimpleEMCALFiducialCut(0.27, 103, 157) ; // EMC 
  //  }
}

///
/// Set common detector fiducial cuts for calorimeter acceptance in isolation analysis.
/// Different settings for neutrals present in the cone or not, less tight phi acceptance.
///
/// \param fidCut  Fiducial cut class pointer
/// \param calorimeter  A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year  The year the data was taken, used to configure some histograms
/// \param coneContent Isolation cone content: kChargedAndNeutral, kOnlyCharged, kOnlyNeutral
///
void CalorimeterFiducialCutForIsolationAnalysis( AliFiducialCut* fidCut, TString calorimeter, 
                                                 Int_t year, Int_t coneContent )
{
  if      ( calorimeter.Contains("CAL" ) ) fidCut->DoEMCALFiducialCut(kTRUE);  
  else if ( calorimeter.Contains("PHOS") ) fidCut->DoPHOSFiducialCut (kTRUE);  
  
  if ( coneContent == AliIsolationCut::kOnlyCharged )
  {
    if      ( calorimeter == "EMCAL" )
    {
      if      ( year > 2014 ) fidCut->SetSimpleEMCALFiducialCut(0.52,  81.2, 185.8) ; //12 SM
      else if ( year > 2010 ) fidCut->SetSimpleEMCALFiducialCut(0.52,  81.2, 178.8) ; //10 SM
      else                    fidCut->SetSimpleEMCALFiducialCut(0.52,  81.2, 118.8) ; // 4 SM
    }
    else if ( calorimeter == "DCAL"  ) 
    {
//    if ( year == 2016 )     fidCut->SetSimpleEMCALFiducialCut(0.52, 261.2, 318.8) ;
//    else                    fidCut->SetSimpleEMCALFiducialCut(0.52, 261.2, 325.8) ;
      
      if ( year == 2016 ) // DCal 2 Regions
      {
        Float_t etamax[] = { 0.52,-0.25 };
        fidCut->AddEMCALFidCutMaxEtaArray(2,etamax);
        Float_t etamin[] = { 0.25,-0.52 };
        fidCut->AddEMCALFidCutMinEtaArray(2,etamin);
        
        Float_t phimax[] = { 318.8, 318.8 };
        fidCut->AddEMCALFidCutMaxPhiArray(2,phimax);
        Float_t phimin[] = { 261.2, 261.2 };
        fidCut->AddEMCALFidCutMinPhiArray(2,phimin);
      }
      else // DCal 3 regions
      {
        Float_t etamax[] = { 0.52,-0.25, 0.52 };
        fidCut->AddEMCALFidCutMaxEtaArray(3,etamax);
        Float_t etamin[] = { 0.25,-0.52,-0.52 };
        fidCut->AddEMCALFidCutMinEtaArray(3,etamin);
        
        Float_t phimax[] = { 318.8, 318.8, 325.8 };
        fidCut->AddEMCALFidCutMaxPhiArray(3,phimax);
        Float_t phimin[] = { 261.2, 261.2, 321.2 };
        fidCut->AddEMCALFidCutMinPhiArray(3,phimin);
      }
    }
    else if ( calorimeter == "PHOS"  )
    {
      if      ( year > 2014 ) fidCut->SetSimplePHOSFiducialCut (0.12, 250.5, 319.5) ; 
      else                    fidCut->SetSimplePHOSFiducialCut (0.12, 260.5, 319.5) ;
    }
    
    // Acceptance of EMCal+DCal
    if ( kAnaCaloTrackCorr.Contains("FullCalo") && year > 2014) 
    {
      //fidCut->SetSimpleEMCALFiducialCut(0.52,  81.2, 325.8) ; //20 SM
      
      // EMCal, 1 region - DCal 2 regions, remove 1/3 DCal
      if ( year == 2016 )
      {
        Float_t etamax[] = { 0.52, 0.52,-0.25 };
        fidCut->AddEMCALFidCutMaxEtaArray(3,etamax);
        Float_t etamin[] = {-0.52, 0.25,-0.52 };
        fidCut->AddEMCALFidCutMinEtaArray(3,etamin);
        
        Float_t phimax[] = { 185.8, 318.8, 318.8 };
        fidCut->AddEMCALFidCutMaxPhiArray(3,phimax);
        Float_t phimin[] = {  81.2, 261.2, 261.2 };
        fidCut->AddEMCALFidCutMinPhiArray(3,phimin);
      }
      else // EMCal, 1 region - DCal 3 regions
      {
        Float_t etamax[] = { 0.52, 0.52,-0.25, 0.52 };
        fidCut->AddEMCALFidCutMaxEtaArray(4,etamax);
        Float_t etamin[] = {-0.52, 0.25,-0.52,-0.52 };
        fidCut->AddEMCALFidCutMinEtaArray(4,etamin);
        
        Float_t phimax[] = { 185.8, 318.8, 318.8, 325.8 };
        fidCut->AddEMCALFidCutMaxPhiArray(4,phimax);
        Float_t phimin[] = {  81.2, 261.2, 261.2, 321.2 };
        fidCut->AddEMCALFidCutMinPhiArray(4,phimin);
      }
    } // full calo
  } // only charged
  else // neutrals in cone
  { 
    if      ( calorimeter == "EMCAL" )
    {
      if      ( year > 2014 ) fidCut->SetSimpleEMCALFiducialCut(0.52, 90.5, 176.5) ;
      else if ( year > 2010 ) fidCut->SetSimpleEMCALFiducialCut(0.52, 90.5, 169.5) ;
      else                    fidCut->SetSimpleEMCALFiducialCut(0.52, 90.5, 109.5) ;
    }
    else if ( calorimeter == "DCAL"  ) 
    {
//    if ( year == 2016 ) fidCut->SetSimpleEMCALFiducialCut(0.52, 270.5, 309.5) ;
//    else                fidCut->SetSimpleEMCALFiducialCut(0.52, 270.5, 316.5) ;
      
      if ( year == 2016 ) // DCal 2 regions, remove 1/3 DCal
      {
        Float_t etamax[] = { 0.52,-0.25 };
        fidCut->AddEMCALFidCutMaxEtaArray(2,etamax);
        Float_t etamin[] = { 0.25,-0.52 };
        fidCut->AddEMCALFidCutMinEtaArray(2,etamin);
        
        Float_t phimax[] = { 309.5, 309.5 };
        fidCut->AddEMCALFidCutMaxPhiArray(2,phimax);
        Float_t phimin[] = { 270.5, 270.5 };
        fidCut->AddEMCALFidCutMinPhiArray(2,phimin);
      }
      else // DCal 2 regions, remove 1/3 for trigger, not for cone content 
      {
        Float_t etamax[] = { 0.52,-0.25 };
        fidCut->AddEMCALFidCutMaxEtaArray(2,etamax);
        Float_t etamin[] = { 0.25,-0.52 };
        fidCut->AddEMCALFidCutMinEtaArray(2,etamin);
        
        Float_t phimax[] = { 316.5, 316.5 };
        fidCut->AddEMCALFidCutMaxPhiArray(2,phimax);
        Float_t phimin[] = { 270.5, 270.5 };
        fidCut->AddEMCALFidCutMinPhiArray(2,phimin);
      }
    }
    else if ( calorimeter == "PHOS"  )
    {
      if    ( year > 2014 ) fidCut->SetSimplePHOSFiducialCut (0.12, 260.5, 309.5) ; 
      else                  fidCut->SetSimplePHOSFiducialCut (0.12, 270.5, 309.5) ; 
    }
    
    // Acceptance of EMCal+DCal
    if ( kAnaCaloTrackCorr.Contains("FullCalo") ) 
    {
      //fidCut->SetSimpleEMCALFiducialCut(0.67, 90.5, 316.5) ; //20 SM
      
      // EMCal, 1 region - DCal 2 regions, remove 1/3 DCal
      if ( year == 2016 )
      {
        Float_t etamax[] = { 0.52, 0.52,-0.25 };
        fidCut->AddEMCALFidCutMaxEtaArray(3,etamax);
        Float_t etamin[] = {-0.52, 0.25,-0.52 };
        fidCut->AddEMCALFidCutMinEtaArray(3,etamin);
        
        Float_t phimax[] = { 176.5, 309.5, 309.5 };
        fidCut->AddEMCALFidCutMaxPhiArray(3,phimax);
        Float_t phimin[] = {  90.5, 270.5, 270.5 };
        fidCut->AddEMCALFidCutMinPhiArray(3,phimin);
      }
      else // EMCal, 1 region - DCal 2 regions but 1/3 in cone
      {
        Float_t etamax[] = { 0.52, 0.52,-0.25 };
        fidCut->AddEMCALFidCutMaxEtaArray(3,etamax);
        Float_t etamin[] = {-0.52, 0.25,-0.52 };
        fidCut->AddEMCALFidCutMinEtaArray(3,etamin);
        
        Float_t phimax[] = { 176.5, 316.5, 316.5 };
        fidCut->AddEMCALFidCutMaxPhiArray(3,phimax);
        Float_t phimin[] = {  90.5, 270.5, 270.5 };
        fidCut->AddEMCALFidCutMinPhiArray(3,phimin);
      }
    }
  } // neutrals in cone
  
  if ( kAnaCutsString.Contains("TightAcc") )
  {
    if ( calorimeter == "EMCAL" ) fidCut->SetSimpleEMCALFiducialCut(0.27, 103, 157) ; // EMC 
  }
}

///
/// Set common histograms binning
/// plus other analysis common settings like TRD covered super modules
/// the activation of the MC dedicated histograms and the activation of
/// the debug mode.
///
/// \param ana : Analysis task where histograms are created and common settings are needed
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
///
void SetAnalysisCommonParameters(AliAnaCaloTrackCorrBaseClass* ana, TString histoString,
                                 TString calorimeter,   Int_t  year,
                                 TString col,           Bool_t simulation,
                                 Bool_t  printSettings, Int_t  debug)
{
  SetAnalysisMixingCentralityBins(ana, col);
  
  if ( kAnaCutsString.Contains("SelectEmbed") ) 
  {
    ana->SwitchOnEmbededSignalSelection() ;
    ana->SwitchOnFillEmbededSignalHistograms() ;
  }
  
  // Track matching parameters
  AliCaloPID* caloPID = ana->GetCaloPID();

  // EMCAL
  // E/p
  if ( kAnaCutsString.Contains("TMEoP10" ) ) caloPID->SetEOverP(0,10);
  if ( kAnaCutsString.Contains("TMEoP5"  ) ) caloPID->SetEOverP(0,5);
  if ( kAnaCutsString.Contains("TMEoP3"  ) ) caloPID->SetEOverP(0,3);
  if ( kAnaCutsString.Contains("TMEoP2"  ) ) caloPID->SetEOverP(0,2);
  if ( kAnaCutsString.Contains("TMEoP1.7") ) caloPID->SetEOverP(0,1.7);
  if ( kAnaCutsString.Contains("TMEoP1.5") ) caloPID->SetEOverP(0,1.5);

  // if tm = 1, fixed cuts
  caloPID->SetEMCALDEtaCut(0.020);
  caloPID->SetEMCALDPhiCut(0.030);

  // PHOS
  caloPID->SetPHOSDispersionCut(2.5);
  caloPID->SetPHOSRCut(2.);

  //
  // Histograms ranges
  //
  AliHistogramRanges* histoRanges = ana->GetHistogramRanges();
  
  if ( histoString != "" ) 
    ana->AddToHistogramsName(Form("%s_%s",  histoString.Data(), (ana->GetAddedHistogramsStringToName()).Data()) );

  histoRanges->SetHistoPtRangeAndNBins(0, 200, 200) ; // Energy and pt histograms
  
  if ( kAnaCutsString.Contains("MultiIso"))
    histoRanges->SetHistoPtRangeAndNBins(0, 100, 100) ; // Energy and pt histograms
  
  if(calorimeter=="EMCAL")
  {
    ana->SetFirstSMCoveredByTRD(-1);
    
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
      
      if     (year == 2011) ana->SetFirstSMCoveredByTRD( 6);
      else if(year == 2012 ||
              year == 2013) ana->SetFirstSMCoveredByTRD( 4);
    }
    else // Run2
    {
      histoRanges->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 187*TMath::DegToRad(), 5*24+8) ;
      histoRanges->SetHistoXRangeAndNBins(-460,460,230); // QA, revise
      histoRanges->SetHistoYRangeAndNBins(-450,450,225); // QA, revise
    }
    
    histoRanges->SetHistoEtaRangeAndNBins(-0.70, 0.70, 2*48) ;
  }
  else if(calorimeter=="DCAL")
  {
    histoRanges->SetHistoPhiRangeAndNBins(260*TMath::DegToRad(), 327*TMath::DegToRad(), 3*24+8) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.70, 0.70, 2*48) ;
    histoRanges->SetHistoXRangeAndNBins(-460,460,230); // QA, revise
    histoRanges->SetHistoYRangeAndNBins(-450,450,225); // QA, revise
  } 
  else if(calorimeter=="PHOS")
  {
    histoRanges->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 320*TMath::DegToRad(), 70) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
  }
  else if(calorimeter=="CTS")
  {
    histoRanges->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    histoRanges->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  if ( kAnaCaloTrackCorr.Contains("FullCalo") ) 
  {
    histoRanges->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 327*TMath::DegToRad(), 247) ;
    histoRanges->SetHistoXRangeAndNBins(-460,460,230); // QA, revise
    histoRanges->SetHistoYRangeAndNBins(-450,450,225); // QA, revise
  }
  
  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 2.9, 300);
  
  // Invariant mass histo
  histoRanges->SetHistoMassRangeAndNBins(0., 1., 200) ;
  histoRanges->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  
  // check if time calibration is on
  //histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  histoRanges->SetHistoTimeRangeAndNBins(-400.,400,400);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-400, 400, 800);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.06,0.06,120);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.06,0.06,120);
  histoRanges->SetHistodRRangeAndNBins(0.,0.06,60);//QA
  
  // QA, electron, charged
  histoRanges->SetHistoEOverPRangeAndNBins(0,10,200);
  histoRanges->SetHistodEdxRangeAndNBins(0.,200.,200);
  
  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  histoRanges->SetHistoZRangeAndNBins(-350,350,175);
  histoRanges->SetHistoRRangeAndNBins(430,460,30);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  
  // QA, correlation
  if(col=="PbPb")
  {
    histoRanges->SetHistoNClusterCellRangeAndNBins(0,100,100);
    histoRanges->SetHistoNClustersRangeAndNBins(0,100,100);
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,2000,200);
  }
  else
  {
    histoRanges->SetHistoNClusterCellRangeAndNBins(0,50,50);
    histoRanges->SetHistoNClustersRangeAndNBins(0,50,50);
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,200,200);
  }
  
  // xE, zT
  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,200);
  histoRanges->SetHistoHBPRangeAndNBins  (0.,10.,200);
  
  // Isolation
  histoRanges->SetHistoPtInConeRangeAndNBins(0, 50 , 250);
  histoRanges->SetHistoPtSumRangeAndNBins   (0, 100, 250);
  histoRanges->SetHistoPtSumSubRangeAndNBins(-100, 100, 200);
  if ( col=="PbPb" ) 
  {
    histoRanges->SetHistoPtSumRangeAndNBins   (   0, 200, 200);
    histoRanges->SetHistoPtSumSubRangeAndNBins(-150, 150, 300);
  }

  //
  // MC histograms?
  //
  if(simulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else           ana->SwitchOffDataMC() ;
  
  //
  // Specialized histograms on multiplicity
  //
  if ( kAnaCutsString.Contains("HighMult") ) 
    ana->SwitchOnFillHighMultiplicityHistograms();
  else
    ana->SwitchOffFillHighMultiplicityHistograms();
  
  // Consider real calo acceptance (borders, bad map holes)
  // at generator level
  if ( kAnaCutsString.Contains("MCRealCaloAcc") )
    ana->SwitchOnRealCaloAcceptance();
  else
    ana->SwitchOffRealCaloAcceptance();
  
  //
  // Debug
  //
  if(printSettings) ana->Print("");
  
  ana->SetDebug(debug); // 10 for lots of messages
}

///
/// Configure the task doing the first photon cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection ...
///
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param tm : Index with track matching option (0- no TM; 1-Fixed residuals cut; 2-Track pT dependent residuals)
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaPhoton* ConfigurePhotonAnalysis(TString col,           Bool_t simulation,
                                      TString calorimeter,   Int_t year, Int_t tm,
                                      Bool_t  printSettings, Int_t   debug, 
                                      TString histoString )
{
  AliAnaPhoton *ana = new AliAnaPhoton();
  
  // cluster selection cuts
  //
  ana->SwitchOnFiducialCut(); 
  CalorimeterFiducialCut(ana->GetFiducialCut(), calorimeter, year);
  
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection
  if ( kAnaCutsString.Contains("MultiIso") ) ana->SwitchOffFillShowerShapeHistograms();
    
  if ( kAnaCutsString.Contains("PerSM") ) 
    ana->SwitchOnFillShowerShapeHistogramsPerSM(); 
  
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  if(tm) ana->SwitchOnTrackMatchRejection() ;
  else   ana->SwitchOffTrackMatchRejection() ;
  
  // Fill track matching histograms if track matching activated
  // and only once in case of multiple analysis
  if ( tm && !kAnaCutsString.Contains("MultiIso") ) 
  {
    ana->SwitchOnTMHistoFill() ;
    if(tm > 1) ana->SwitchOnTMTrackPtHistoFill() ;
  }
  else
  {
    ana->SwitchOffTMHistoFill() ;
    ana->SwitchOffTMTrackPtHistoFill() ;
  }
  
  if(calorimeter == "PHOS")
  {
    ana->SetNCellCut(2);// At least 3 cells, check if already set at reader level
    ana->SetMinEnergy(0.3);
    ana->SetMaxEnergy(200);
    ana->SetMinDistanceToBadChannel(2, 4, 5); // could have been already applied at reader level
    ana->SetTimeCut(-1e10,1e10); // open cut
  }
  else
  {//EMCAL
    ana->SetNCellCut(1);// At least 2 cells, check if already set at reader level
    ana->SetMinEnergy(0.7); 
    ana->SetMaxEnergy(300);
    ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off
                                 // restrict to less than 100 ns when time calibration is on
    ana->SetMinDistanceToBadChannel(2, 4, 6); // could have been already applied at reader level
    
    // NLM cut, used in all, exclude clusters with more than 2 maxima
    // Not needed if M02 cut is already strong or clusterizer V2
    ana->SetNLMCut(1, 2) ;
  }
  
  if ( kAnaCutsString.Contains("ShSh") )
  {
    printf("ConfigurePhotonAnalysis() >>> Recalculate shower shape within NxN \n");
    // Open it although it should not matter for V3-V2 clusterizers
    ana->SetNLMCut(0,100);

    // Make sure isolated cells are not clusterized, only neighbours
    ana->SwitchOnNxNShowerShapeOnlyNeighbours();

    // Set the size of the recalculation window
    if      ( kAnaCutsString.Contains("ShSh5x5") )
      ana->SwitchOnUse5x5ShowerShapeHisto();
    else if ( kAnaCutsString.Contains("ShSh7x7") )
      ana->SwitchOnUse7x7ShowerShapeHisto();
  }

  // PID cuts (shower shape and track matching)
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
  
  // EMCAL
  //caloPID->SetEMCALLambda0CutMax(0.27);
  caloPID->SetEMCALLambda0CutMax(10); // open, full shower shape needed for isolation studies
  caloPID->SetEMCALLambda0CutMin(0.10);
  
  // Track matching, pT track dependent cuts
  // E/p cuts and residual cuts for tm=1 set on SetAnalysisCommonParameters()
  if ( tm > 1 ) caloPID->SwitchOnEMCTrackPtDepResMatching();
  
  // Branch AOD settings
  ana->SetOutputAODName(Form("PhotonTrigger_%s",kAnaCaloTrackCorr.Data()));
  ana->SetOutputAODClassName("AliCaloTrackParticleCorrelation");
  
  //Set Histograms name tag, bins and ranges
  //ana->AddToHistogramsName(Form("AnaPhoton_TM%d_",tm));
  ana->AddToHistogramsName("AnaPhoton_");
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug) ; // see method below
  
  if(ana->GetFirstSMCoveredByTRD() > 0)
    printf("ConfigurePhotonAnalysis() >>> Set first SM covered by TRD, SM=%d <<< year %d \n", ana->GetFirstSMCoveredByTRD(),year);
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms (17); // 18 max
  ana->FillNPrimaryHistograms(6); // 6 max

  return ana;
}

///
/// Configure the task doing electron (or charged hadron) cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection, dEdx, E/p, TPC nSigma.
///
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaElectron* ConfigureElectronAnalysis(TString col,           Bool_t simulation,
                                          TString calorimeter,   Int_t  year,
                                          Bool_t  printSettings, Int_t  debug, 
                                          TString histoString)
{
  AliAnaElectron *ana = new AliAnaElectron();
  ana->SetDebug(debug); //10 for lots of messages
  
  ana->SwitchOnFiducialCut(); 
  CalorimeterFiducialCut(ana->GetFiducialCut(), calorimeter, year);
  
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  // What particles to select:
  //ana->FillAODWithElectrons();
  //ana->FillAODWithHadrons();
  ana->FillAODWithAny(); // Both
  
  if(calorimeter == "PHOS")
  {
    ana->SetNCellCut(2);// At least 2 cells
    ana->SetMinEnergy(0.3);
    ana->SetMaxEnergy(200);
    ana->SetMinDistanceToBadChannel(2);
  }
  else 
  {// EMCAL
    ana->SetNCellCut(1);// At least 2 cells
    ana->SetMinEnergy(0.7); // no effect minium EMCAL cut.
    ana->SetMaxEnergy(300); 
    //ana->SetTimeCut(400,900);// Time window of [400-900] ns
    ana->SetMinDistanceToBadChannel(2);
  }
  
  // Electron selection cuts with tracks
  //ana->SetM02Range(0.05, 0.35);
  //ana->SetM20Range(-10, 0.3);
  //ana->SetEOverP(0.9, 1.2);
  //ana->SetNSigma(-1, 3);
  //ana->SetNSigmaForHadron(-10, -4);
  
//  if     (kRunNumber < 146861) ana->SetdEdxCut(72, 90);
//  else if(kRunNumber < 154000) ana->SetdEdxCut(54, 70);
//  else                         ana->SetdEdxCut(74, 90);
//  
//  if(simulation)  ana->SetdEdxCut(80, 100);
      
  ana->SwitchOnFillShowerShapeHistograms();  
  
  ana->SwitchOffFillWeightHistograms()  ;
    
  ana->SetOutputAODName(Form("Electron_%s",kAnaCaloTrackCorr.Data()));
  ana->SetOutputAODClassName("AliCaloTrackParticleCorrelation");

  // Set Histograms name tag, bins and ranges
  
  //ana->AddToHistogramsName(Form("AnaElectron_TM%d_",tm));
  ana->AddToHistogramsName("AnaElectron_");
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug) ; // see method below
    
  //if(printSettings) 
  ana->Print("");
  
  return ana ;
}

///
/// Configure the task doing the pi0 even by event selection (invariant mass or split)
/// and the cluster tagging as decay in different mass windows.
///
/// \param particle : String with particle type, "Pi0" or "Eta". If string contains "SideBand", not peak but side band inspected    
/// \param analysis : analysis type, merged clusters or cluster selected from invariant mass pair or side band
/// \param useSSIso : bool, select pairs tagged as isolated previously  
/// \param useAsy : bool, for merged cluster analysis apply asymmetry cut or not
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param tm : Index with track matching option (0- no TM; 1-Fixed residuals cut; 2-Track pT dependent residuals)
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaPi0EbE* ConfigurePi0EbEAnalysis(TString particle,      Int_t  analysis,
                                      Bool_t useSSIso,       Bool_t useAsy,
                                      TString col,           Bool_t simulation,
                                      TString calorimeter,   Int_t  year,  Int_t tm,
                                      Bool_t  printSettings, Int_t debug, 
                                      TString histoString )
{
  // Configuration of pi0 event by event selection
  
  AliAnaPi0EbE *ana = new AliAnaPi0EbE();
  
  ana->SetDebug(debug);
  
  ana->SetAnalysisType((AliAnaPi0EbE::anaTypes)analysis);
  TString opt = "";
  if(analysis==AliAnaPi0EbE::kIMCaloTracks) opt = "Conv";
  if(analysis==AliAnaPi0EbE::kSSCalo)       opt = "SS";
  
  ana->SwitchOffAllNLMHistoFill();
  ana->SwitchOffSelectedClusterHistoFill();
  
  ana->SwitchOffFillWeightHistograms();
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  //if(!kTime && !simulation) ana->SwitchOnFillEMCALBCHistograms();
  
  if(tm) ana->SwitchOnTrackMatchRejection() ;
  else   ana->SwitchOffTrackMatchRejection() ;
  
  if ( tm && !kAnaCutsString.Contains("Photon") &&  !kAnaCutsString.Contains("MultiIso") )
    ana->SwitchOnTMHistoFill() ;
  else
    ana->SwitchOffTMHistoFill() ;
    
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  // Branch AOD settings
  ana->SetOutputAODName(Form("%s%sTrigger_%s",particle.Data(), opt.Data(), kAnaCaloTrackCorr.Data()));
  printf("ConfigurePi0EbEAnalysis() *** Out branch %s***\n",ana->GetOutputAODName().Data());
  ana->SetOutputAODClassName("AliCaloTrackParticleCorrelation");
  
  if(analysis == AliAnaPi0EbE::kIMCaloTracks) ana->SetInputAODGammaConvName("PhotonsCTS");
  
  //Set Histograms name tag, bins and ranges
  
  //ana->AddToHistogramsName(Form("Ana%s%sEbE_TM%d_",particle.Data(),opt.Data(),tm));
  ana->AddToHistogramsName(Form("Ana%s%sEbE_",particle.Data(),opt.Data()));
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  ///////////////////////////////////
  if(analysis!=AliAnaPi0EbE::kSSCalo)
  {
    ana->SetInputAODName(Form("PhotonTrigger_%s",kAnaCaloTrackCorr.Data()));
    
    ana->SetM02CutForInvMass(0.1,0.35); // Loose SS cut
    
    ana->SwitchOnSelectPairInIsolationCone();
    ana->SetR(0.4);
    ana->SetIsolationCandidateMinPt(5);
    
    if(useSSIso)
    {
      ana->SwitchOnSelectIsolatedDecay();
      //ana->AddToHistogramsName(Form("Ana%s%sEbEIsoDecay_TM%d_",particle.Data(),opt.Data(),tm));
      ana->AddToHistogramsName(Form("Ana%s%sEbEIsoDecay_",particle.Data(),opt.Data()));
      ana->SetOutputAODName(Form("%s%sIsoDecayTrigger_%s",particle.Data(), opt.Data(), kAnaCaloTrackCorr.Data()));
    }
    
    if ( calorimeter.Contains("CAL") && !simulation ) ana->SetPairTimeCut(100);
    
    AliNeutralMesonSelection *nms = ana->GetNeutralMesonSelection();
    nms->SetParticle(particle);
    
    //****
    nms->SetInvMassCutMaxParameters(0,0,0); // Overrule the setting in SetParticle for Pi0 option
                                            //****
    
    // Tighten a bit mass cut with respect to default window
    if(particle=="Pi0") nms->SetInvMassCutRange(0.110,0.160);
    if(particle=="Eta") nms->SetInvMassCutRange(0.520,0.580);
    
    //if(!particle.Contains("SideBand")) nms->SwitchOnAngleSelection();
    //else nms->SwitchOnAngleSelection();
    
    nms->SwitchOffAngleSelection();
    
    if(particle.Contains("Pi0SideBand")) // For pi0, do not consider left band
      nms->SetSideBandCutRanges(-1,0,0.190,0.240);
    
    if(particle.Contains("EtaSideBand")) // For pi0, do not consider left band
      nms->SetSideBandCutRanges(0.410,0.470,0.620,0.680);
    
    nms->KeepNeutralMesonSelectionHistos(kTRUE);
    //nms->SetAngleMaxParam(2,0.2);
    nms->SetHistoERangeAndNBins(0, 20, 80) ;
    //nms->SetHistoIMRangeAndNBins(0, 1, 400);
  }
  else  ///////////////////////////////////
  {
    ana->SwitchOnFiducialCut(); 
    CalorimeterFiducialCut(ana->GetFiducialCut(), calorimeter, year);
 
    // cluster splitting settings
    ana->SetMinEnergy(6);
    ana->SetMaxEnergy(300.);
    
    ana->SetNLMMinEnergy(0, 10);
    ana->SetNLMMinEnergy(1, 6);
    ana->SetNLMMinEnergy(2, 6);
    
    // NLM cut, used in all, exclude clusters with more than 2 maxima
    ana->SetNLMCut(1, 2) ;
    
    //
    ana->SetMinDistanceToBadChannel(2, 4, 6);
    ana->SwitchOnSplitClusterDistToBad();
    ana->SetTimeCut(-1e10,1e10); // Open time cut
    
    AliCaloPID* caloPID = ana->GetCaloPID();
    // pT track dependent cuts
    if(tm > 1) caloPID->SwitchOnEMCTrackPtDepResMatching();
    
    caloPID->SetSplitWidthSigma(3); // cut at 3 sigma of the mean pi0 peak.
    
    if(!useSSIso)
    {
      printf("ConfigurePi0EbEAnalysis() << Do not apply SS cut on merged pi0 analysis >> \n");
      caloPID->SwitchOffSplitShowerShapeCut() ;
      //ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_TM%d_",particle.Data(),opt.Data(),tm));
      ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_",particle.Data(),opt.Data()));
      ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenSS",particle.Data(), opt.Data(), kAnaCaloTrackCorr.Data()));
      caloPID->SetClusterSplittingM02Cut(0.1,10);
    }
    else
    {
      caloPID->SetClusterSplittingM02Cut(0.3,4); // Do the selection in the analysis class and not in the PID method to fill SS histograms
      caloPID->SwitchOnSplitShowerShapeCut() ;
    }
    
    if(useAsy)
    {
      caloPID->SwitchOnSplitAsymmetryCut() ;
      ana->GetCaloPID()->SetSubClusterEnergyMinimum(0,2);
      ana->GetCaloPID()->SetSubClusterEnergyMinimum(1,0.5);
      ana->GetCaloPID()->SetSubClusterEnergyMinimum(2,0.5);
    }
    else
    {
      caloPID->SwitchOffSplitAsymmetryCut() ;
      if(!useSSIso)
      {
        //ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_OpenAsy_TM%d_",particle.Data(),opt.Data(),tm));
        ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_OpenAsy_",particle.Data(),opt.Data()));
        ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenSS_OpenAsy",particle.Data(), opt.Data(), kAnaCaloTrackCorr.Data()));
      }
      else
      {
        //ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenAsy_TM%d_",particle.Data(),opt.Data(),tm));
        ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenAsy_",particle.Data(),opt.Data()));
        ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenAsy",particle.Data(), opt.Data(), kAnaCaloTrackCorr.Data()));
      }
    }
    
    // For Pi0 only if  SwitchOnSimpleSplitMassCut()
    caloPID->SetPi0MassRange(0.10, 0.18);
    caloPID->SetEtaMassRange(0.50, 0.60);
    caloPID->SetPhotonMassRange(0.00, 0.08);
    
    caloPID->SetClusterSplittingMinNCells(6);
    
    //caloPID->SetSplitEnergyFractionMinimum(0, 0.95);
    //caloPID->SetSplitEnergyFractionMinimum(1, 0.95);
    //caloPID->SetSplitEnergyFractionMinimum(2, 0.8);
    
    if(col=="PbPb" || kAnaCaloTrackCorr.Contains("150"))
    {
      caloPID->SetClusterSplittingMinNCells(4);
      //caloPID->SetPi0MassShiftHighECell(0.005);
    }
  }
  ///////////////////////////////////
  
  return  ana;
}


///
/// Configure the task doing the 2 cluster invariant mass analysis
/// from selected clusters by AliAnaPhoton
///
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param bothCalo : Activate combination of photon pairs from different detectors (PHOS+EMCal)
/// \param year : The year the data was taken, used to configure some histograms
/// \param mixOn : bool to activate mixing analysis
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaPi0* ConfigureInvariantMassAnalysis
(TString col,           Bool_t simulation,
 TString calorimeter,   Bool_t bothCalo,   
 Int_t   year,          Bool_t mixOn,
 Bool_t  printSettings, Int_t  debug, 
 TString histoString )
{
  AliAnaPi0 *ana = new AliAnaPi0();
  
  ana->SetDebug(debug);
  
  // Input delta AOD settings
  ana->SetInputAODName(Form("PhotonTrigger_%s",kAnaCaloTrackCorr.Data()));

  if(bothCalo)
  {
    ana->SwitchOnPairWithOtherDetector();
    
    TString otherDetector = Form("PhotonTrigger_%s",kAnaCaloTrackCorr.Data());
    
    if ( calorimeter == "EMCAL" ) 
      otherDetector.ReplaceAll("PhotonTrigger_EMCAL","PhotonTrigger_PHOS");
    else    
      otherDetector.ReplaceAll("PhotonTrigger_PHOS","PhotonTrigger_EMCAL");
    
    ana->SetOtherDetectorInputName(otherDetector);
  }
  
  // Calorimeter settings
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  // Acceptance plots
  //
  ana->SwitchOnFiducialCut();
  CalorimeterFiducialCut(ana->GetFiducialCut(), calorimeter, year);

  // settings for pp collision mixing
  if(mixOn) ana->SwitchOnOwnMix();
  else      ana->SwitchOffOwnMix();
  
  // Cuts
  if (calorimeter == "EMCAL" ) 
  {
    if(year < 2014) ana->SetPairTimeCut(50);
    else            ana->SetPairTimeCut(200); // REMEMBER to remove this when time calib is on
  }
  
  ana->SetNPIDBits(1);
  ana->SetNAsymCuts(1); // no asymmetry cut, previous studies showed small effect.
                        // In EMCAL assymetry cut prevents combination of assymetric decays which is the main source of pi0 at high E.
  
  if     (col == "pp"  )
    ana->SetMinPt(0.7);
  else if(col == "PbPb")
    ana->SetMinPt(1.5);
  else if(col =="pPb")
    ana->SetMinPt(0.7);
  
  // Angle cut, avoid pairs with too large angle
  ana->SwitchOnAngleSelection(); 
  ana->SetAngleMaxCut(TMath::DegToRad()*80.); // EMCal: 4 SM in phi, 2 full SMs in eta
  ana->SetAngleCut(0.017); // Minimum angle open, cell size
  
  if ( kAnaCutsString.Contains("PerSM") ) //&& !bothCalo ) 
    ana->SwitchOnSMCombinations();
  
  ana->SwitchOffMultipleCutAnalysis();
  ana->SwitchOnFillAngleHisto();
  ana->SwitchOnFillOriginHisto(); // MC
  
  // Set Histograms name tag, bins and ranges
  if(!bothCalo) ana->AddToHistogramsName(Form("AnaPi0_"));
  else          ana->AddToHistogramsName(Form("AnaPi0_2Det_"));
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below
  if(printSettings) ana->Print("");
  
  return ana;
}

///
/// Configure the task generating random ghost particles
///
/// \param detector : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaRandomTrigger* ConfigureRandomTriggerAnalysis
( TString detector     , Int_t  year,
  Bool_t  printSettings, Int_t  debug, 
  TString histoString
 )
{
  AliAnaRandomTrigger *ana = new AliAnaRandomTrigger();
  ana->SetDebug(debug); //10 for lots of messages
  
  ana->SetTriggerDetector(detector);

  // selection cuts
  
  // A Flat pT distribution within these numbers
  ana->SetMinPt(10.); 
  ana->SetMaxPt(20.);   
  
  // A flat eta and phi distribution limited to known detectors
  
  ana->SwitchOnFiducialCut();
  CalorimeterFiducialCut(ana->GetFiducialCut(), detector, year);
  
  // First broad cut, then also rely on predefined fiducial cut set above
  if      ( detector=="EMCAL" )
  {
    ana->SetEtaCut(-0.7,0.7);
    ana->SetPhiCut( 80*TMath::DegToRad(), 187*TMath::DegToRad());
  }
  else if ( detector=="DCAL" )
   {
     ana->SetEtaCut(-0.7,0.7);
     ana->SetPhiCut(260*TMath::DegToRad(), 327*TMath::DegToRad());
   }
  else if ( detector=="PHOS" )
  {
    ana->SetEtaCut(-0.13,0.13);
    ana->SetPhiCut(250*TMath::DegToRad(), 320*TMath::DegToRad());
  }
  else if ( detector=="CTS" )
  {
    ana->SetEtaCut(-0.9,0.9);
    ana->SetPhiCut(0, TMath::TwoPi());
  }
  
  // AOD branch
  ana->SetOutputAODName(Form("RandomTrigger_%s",kAnaCaloTrackCorr.Data()));
  ana->SetOutputAODClassName("AliCaloTrackParticleCorrelation");
  
  //Set Histograms name tag, bins and ranges
  ana->AddToHistogramsName("AnaRandom_");
  
  SetAnalysisCommonParameters(ana,histoString,detector,year,"pp",kFALSE,printSettings,debug) ; // see method below

  // Reset histogram bins, no need for fine bins
  ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 360) ;
  ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1, 1, 200) ;
  ana->GetHistogramRanges()->SetHistoPtRangeAndNBins(0, 100, 1) ;
  
  if ( printSettings ) ana->Print("");
  
  return ana;
}

///
/// Configure the isolation cuts 
/// used in ConfigureIsolationAnalysis() and ConfigureHadronCorrelationAnalysis()
///
/// \param ic : Pointer to task doing the isolation
/// \param partInCone : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
/// \param thresType : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
/// \param cone : A float setting the isolation cone size higher limit
/// \param coneMin : A float setting the isolation cone size lower limit
/// \param pth : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param col : A string with the colliding system
/// \param debug : An int to define the debug level of all the tasks
///
void ConfigureIsolationCut(AliIsolationCut * ic,
                           Int_t partInCone, Int_t thresType, 
                           Float_t cone, Float_t coneMin,
                           Float_t pth, TString col, Int_t debug)
{
  ic->SetDebug(debug);
  ic->SetParticleTypeInCone(partInCone);
  ic->SetICMethod(thresType);
  ic->SetPtFraction(0.1);
  ic->SetPtThreshold(0.5); // default, change in next lines
  ic->SetSumPtThreshold(1.0); // default, change in next lines
  
  if ( cone > 0 && pth > 0 )
  {
    ic->SetConeSize(cone);
    ic->SetMinDistToTrigger(coneMin);    
    ic->SetPtThresholdMax(10000);
    
    if ( thresType == AliIsolationCut::kPtThresIC )
    {
      printf("ConfigureIsolationCuts() *** PtThresMin = %1.1f GeV/c *** R = %1.2f *** R min %1.2f\n",pth,cone,coneMin);
      ic->SetPtThreshold(pth);
    }
    
    if ( thresType == AliIsolationCut::kSumPtIC || 
         thresType >= AliIsolationCut::kSumBkgSubIC )
    {
      printf("ConfigureIsolationCuts() *** SumPtMin = %1.1f GeV/c *** R = %1.1f *** R min %1.2f\n",pth,cone,coneMin);
      ic->SetSumPtThreshold(pth);
    }
  }
  else
  {
    printf("ConfigureIsolationCuts() *** Careful, use old hardcoded values\n");
    if ( col == "pp" || col == "pPb" )
    {
      ic->SetPtThreshold(0.5);
      ic->SetSumPtThreshold(1.0) ;
      ic->SetConeSize(0.4);
      ic->SetMinDistToTrigger(-1);
    }
    if ( col == "PbPb" )
    {
      ic->SetPtThreshold(3.);
      ic->SetSumPtThreshold(3.0) ;
      ic->SetConeSize(0.3);
      ic->SetMinDistToTrigger(-1);  
    }
  }
  
  //ic->SwitchOnFillEtaPhiHistograms();
  
  if ( kAnaCutsString.Contains("HighMult") )
    ic->SwitchOnFillHighMultHistograms (); 
  else
    ic->SwitchOffFillHighMultHistograms ();
      
  if ( kAnaCutsString.Contains("FixIsoConeExcess") ) 
    ic->SwitchOnConeExcessCorrection();
  else                                         
    ic->SwitchOffConeExcessCorrection();
  
  ic->SetConeSizeBandGap(0.0);
  if ( kAnaCutsString.Contains("IsoBandUEGap") && thresType > AliIsolationCut::kSumBkgSubIC)
  {
    // do not count UE particles near the cone limit > R+0.05

    if ( kAnaCutsString.Contains("IsoBandUEGapFix05") )
      ic->SetConeSizeBandGap(0.5-cone); // the UE band excludes a cone size of R=0.5
    else
      ic->SetConeSizeBandGap(0.1);      // the UE band excludes a cone size with R=cone+0.1

    printf("\t Add isolation UE gap %1.2f\n",ic->GetConeSizeBandGap());
  }
  
  // Set ratio of neutral energy over charged energy
  // Needed for UE subtraction using perp cones and using neutral in cone   
  // 
  if ( col == "pp" || col == "pPb" )
  {
    ic->SetNeutralOverChargedRatio(0.363,0.0,0.0,0.0); // Erwann Thesis p-Pb
  }
  if ( col == "PbPb" )
  {
    ic->SetNeutralOverChargedRatio(1.49e-1,-2.54e-3,0.0,7.32e-7); // LHC18qr, random cone, eta-band
  //ic->SetNeutralOverChargedRatio(1.29e-1,-2.47e-3,0.0,7.38e-7); // LHC18qr, random cone, phi-band
  //ic->SetNeutralOverChargedRatio(1.26e-1,-2.27e-3,0.0,6.91e-7); // LHC18qr, random cone
    
  //ic->SetNeutralOverChargedRatio(2.13e-1,-3.27e-3,0.0,8.86e-7); // LHC11h, random cone, eta-band
  //ic->SetNeutralOverChargedRatio(1.94e-1,-3.40e-3,0.0,9.37e-7); // LHC11h, random cone, phi-band
  //ic->SetNeutralOverChargedRatio(1.82e-1,-3.04e-3,0.0,8.40e-7); // LHC11h, random cone
  }
}

///
/// Configure the task doing the trigger particle isolation
///
/// \param particle : Particle trigger name
/// \param leading : An int setting the type of leading particle selection: 0, select all;l 1: absolute  leading of charged; 2: absolute  leading of charged and neutral; 3: near side leading absolute of charged; 4: near side leading absolute of charged and neutral
/// \param partInCone : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
/// \param thresType : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
/// \param cone : A float setting the isolation cone size higher limit
/// \param coneMin : A float setting the isolation cone size lower limit
/// \param pth : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param tm : Index with track matching option (0- no TM; 1-Fixed residuals cut; 2-Track pT dependent residuals)
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle,      Int_t   leading,
                                                    Int_t   partInCone,    Int_t   thresType,
                                                    Float_t cone,          Float_t coneMin,
                                                    Float_t pth,           
                                                    TString col,           Bool_t  simulation,
                                                    TString calorimeter,   Int_t   year,       Int_t tm,
                                                    Bool_t  printSettings, Int_t   debug, 
                                                    TString histoString                 )
{
  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
  
  ana->SetDebug(debug);
  
  ana->SetMinPt(5);
  ana->SetMaxPt(200);
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
    
  ana->SwitchOffLeadingOnly();
  ana->SwitchOffCheckNeutralClustersForLeading();
  if( leading > 0 )   ana->SwitchOnLeadingOnly();
  if( leading == 2 ||
     leading == 4)   ana->SwitchOnCheckNeutralClustersForLeading();
  
  // Do at generation level detector cuts and effects
  ana->SwitchOnPrimariesInConeSelection();
  ana->SwitchOffPrimariesPi0DecayStudy() ;
  
  ana->SwitchOffSSHistoFill();

  if(particle.Contains("Photon"))
  {
    ana->SwitchOnSSHistoFill();

    if(kAnaCutsString.Contains("Decay"))
    {
      ana->SwitchOffDecayTaggedHistoFill() ;
      ana->SetNDecayBits(5);
    }
  }
 
  if ( kAnaCutsString.Contains("PerSM") ) 
    ana->SwitchOnFillHistogramsPerSM(); 

  if ( kAnaCutsString.Contains("PerNCells") ) 
    ana->SwitchOnStudyNCellsCut();
  
  if ( kAnaCutsString.Contains("PerTCard") ) 
    ana->SwitchOnFillHistogramsPerTCardIndex(); 
  
  if(!tm)  ana->SwitchOnTMHistoFill();
  else     ana->SwitchOffTMHistoFill();
  
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  // Avoid borders of calorimeter or TPC as much as possible  
  // Open phi acceptance in case of charged only in isolation cone
  
  ana->SwitchOnFiducialCut(); 
  CalorimeterFiducialCutForIsolationAnalysis(ana->GetFiducialCut(), calorimeter, 
                                             year, partInCone);
   
  // Same Eta as EMCal, cut in phi if EMCAL was triggering
  if(particle=="Hadron"  || particle.Contains("CTS") )
  {
    //if(trigger.Contains("EMC"))
    //  ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
    //else
    ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;
  }
  
  // Branch AOD settings
  
  ana->SetInputAODName(Form("%sTrigger_%s",particle.Data(),kAnaCaloTrackCorr.Data()));
  ana->SetAODObjArrayName(Form("IC%sTrigger_%s_R%1.1f_ThMin%1.1f",particle.Data(),kAnaCaloTrackCorr.Data(),cone,pth));
  
  //
  // Do settings for main isolation cut class
  //
  AliIsolationCut * ic =  ana->GetIsolationCut();
  ConfigureIsolationCut(ic,partInCone,thresType,cone,coneMin,pth,col,debug);
 
  // Track matching 
  //
  // By default apply always a TM cut in cone if charged and neutrals
  // in isolation cone
  if ( partInCone == AliIsolationCut::kNeutralAndCharged )
  {
    AliCaloPID* caloPID = ana->GetCaloPID();

    if ( tm > 1 )
      caloPID->SwitchOnEMCTrackPtDepResMatching();
  }
  
  // Set Histograms name tag, bins and ranges
  TString histoStringTag = Form("AnaIsol%s_",particle.Data());
  if ( kAnaCutsString.Contains("MultiIso") )
  {
    if ( kAnaCutsString.Contains("UESubMethods") )
    {
      if      ( thresType == AliIsolationCut::kSumBkgSubIC)        histoStringTag += Form("PerpCone_");
      else if ( thresType == AliIsolationCut::kSumBkgSubEtaBandIC) histoStringTag += Form("EtaBand_" );
      else if ( thresType == AliIsolationCut::kSumBkgSubPhiBandIC) histoStringTag += Form("PhiBand_" );
      else printf("--- Isolation method for method %d not added\n",thresType); 

//      if ( partInCone == AliIsolationCut::kOnlyCharged )
//      {
//        histoStringTag += Form("OnlyCharged_");
//      }
    }

    if ( kAnaCutsString.Contains("IsoR"))
      histoStringTag += Form("R%0.2f_",cone);
  }

  ana->AddToHistogramsName(histoStringTag);
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  if ( particle=="Hadron"  || particle.Contains("CTS") )
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  if ( particle == "Random" )
  {
    ana->GetHistogramRanges()->SetHistoPtRangeAndNBins(0, 100, 1) ;
    ana->GetHistogramRanges()->SetHistoShowerShapeRangeAndNBins(-0.1, 2.9, 1);

    ic->SwitchOnFillEtaPhiHistograms();
    ana->SwitchOffSSHistoFill();
    ana->SwitchOffTMHistoFill();
    ana->SwitchOnIsolationControlHistoFill();
    ana->SwitchOffNonConstantPtBinHistoArray();
    ana->SwitchOnStudyPtCutInCone();
  }

  if ( printSettings ) ic ->Print("");
  
  return ana;
}

///
/// Configure the task doing the trigger particle hadron correlation
///
/// \param particle : Particle trigger name
/// \param leading : An int setting the type of leading particle selection: 0, select all;l 1: absolute  leading of charged; 2: absolute  leading of charged and neutral; 3: near side leading absolute of charged; 4: near side leading absolute of charged and neutral
/// \param bIsolated : Bool setting the analysis for previously isolated triggers
/// \param partInCone : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
/// \param thresType : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
/// \param cone : A float setting the isolation cone size higher limit
/// \param coneMin : A float setting the isolation cone size lower limit
/// \param pth : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param mixOn : A bool to switch the correlation mixing analysis
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param tm : Index with track matching option (0- no TM; 1-Fixed residuals cut; 2-Track pT dependent residuals). In use in case of mixing and neutrals used, specially in isolation.
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle,      Int_t   leading,
                                                                    Bool_t  bIsolated,     Float_t shshMax,
                                                                    Int_t   partInCone,    Int_t   thresType,
                                                                    Float_t cone,          Float_t coneMin,
                                                                    Float_t pth,           Bool_t  mixOn,
                                                                    TString col,           Bool_t  simulation,
                                                                    TString calorimeter,   Int_t   year,     Int_t tm,
                                                                    Bool_t  printSettings, Int_t   debug, 
                                                                    TString histoString                       )
{
  AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
  
  ana->SetTriggerPtRange(5,100);
  ana->SetAssociatedPtRange(0.2,100);
  ana->SetDeltaPhiCutRange  (TMath::DegToRad()*120.,TMath::DegToRad()*240.);
  
  // Underlying event
  ana->SetUeDeltaPhiCutRange(TMath::DegToRad()*60. ,TMath::DegToRad()*120.);
  ana->SwitchOnSeveralUECalculation();
  
  ana->SwitchOffAbsoluteLeading();  // Select trigger leading particle of all the selected tracks
  ana->SwitchOffNearSideLeading();  // Select trigger leading particle of all the particles at +-90 degrees, default
  ana->SwitchOffCheckNeutralClustersForLeading();
  
  if(leading >  0 && leading <  3 ) ana->SwitchOnAbsoluteLeading();
  if(leading >  2 )                 ana->SwitchOnNearSideLeading();
  if(leading == 2 || leading == 4 ) ana->SwitchOnCheckNeutralClustersForLeading();
  
  ana->SwitchOffFillPtImbalancePerPtABinHistograms();
  ana->SwitchOffCorrelationVzBin() ;
  ana->SwitchOffFillEtaGapHistograms();
    
  ana->SwitchOffPi0TriggerDecayCorr();
  
  if(particle.Contains("Photon"))
  {
    if(kAnaCutsString.Contains("Decay"))
    {
      printf("ConfigureHadronCorrelationAnalysis() *** Activate analysis on Tagged clusters as Decay product >> \n");
      ana->SwitchOnDecayTriggerDecayCorr();
      ana->SetNDecayBits(5);
      ana->SwitchOnInvariantMassHistograms();
      if ( kAnaCutsString.Contains("Bkg") && !bIsolated )
        ana->SwitchOnBackgroundBinsTaggedDecayPtInConeHistograms();
    }
    
    //printf("ConfigureHadronCorrelationAnalysis() *** SET M02 limits in correlation task *** \n");
    ana->SetM02Cut(0.10,shshMax);
  }
  
  if ( kAnaCutsString.Contains("PerSM") ) 
    ana->SwitchOnFillHistogramsPerSM(); 
 
  if ( kAnaCutsString.Contains("PerTCard") ) 
    ana->SwitchOnFillHistogramsPerTCardIndex(); 
  
  ana->SetMCGenType(0,7); // Change to higher to include Tagged decays
  
  ana->SwitchOffLeadHadronSelection(); // Open cuts, just fill histograms
  ana->SwitchOffFillLeadHadronHistograms();
  
  if ( kAnaCutsString.Contains("Bkg") && !bIsolated )
  {
    printf("ConfigureHadronCorrelationAnalysis() *** Activate analysis on PtTrig and Bkg bins,\n"
           " \t make sure isolation runs first even if not used to select trigger >> \n");
    ana->SwitchOnBackgroundBinsPtInConeHistograms();
  }
  
  ana->SetLeadHadronPhiCut(TMath::DegToRad()*130, TMath::DegToRad()*230.);
  ana->SetLeadHadronPtCut(0.5, 1000);
  
  // if triggering on PHOS and EMCAL is on
  ana->SwitchOffNeutralCorr(); // Do only correlation with TPC
                               //ana->SetPi0AODBranchName("Pi0EMCAL_TrigEMC7_Cl_TM1");
  
  ana->SwitchOffHMPIDCorrelation();
  
  ana->SwitchOffFillBradHistograms();
  
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  ana->SwitchOnFillDeltaEtaPhiPtTrigHistograms();
  
  ana->SetNAssocPtBins(16); // set last bin [20,30] GeV/c
  // See AliAnaParticleCorrelation::InitParameters();
  // Default bins{0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,16,20,30,40,50,100} GeV/c
  // If you want to change it:
  //  ana->SetAssocPtBinLimit(0, 1) ;
  //  ana->SetAssocPtBinLimit(1, 2) ;
  //  ana->SetAssocPtBinLimit(2, 3) ;
  //  ana->SetAssocPtBinLimit(3, 4) ;
  //  ana->SetAssocPtBinLimit(4, 5) ;
  //  ana->SetAssocPtBinLimit(5, 8) ;
  //  ana->SetAssocPtBinLimit(6, 10) ;
  //  ana->SetAssocPtBinLimit(7, 100);
  
  ana->SetNTriggerPtBins(8); // set first bin [10,12] GeV/c and last bin [40,50] GeV/c
  // See AliAnaParticleCorrelation::InitParameters();
  // Default bins{10,12,16,20,25,30,40,50,75,100} GeV/c
  // If you want to change it:
  //  ana->SetTriggerPtBinLimit(0,  5) ;
  //  ana->SetTriggerPtBinLimit(1, 10) ;
  
  ana->SelectIsolated(bIsolated); // do correlation with isolated photons
  
  // Mixing with own pool
  if(mixOn)
  {
    ana->SwitchOnOwnMix();
    ana->SwitchOnFillNeutralInMixedEvent();
    
    if(bIsolated)
    {
      //Do settings for main isolation cut class
      AliIsolationCut * ic =  ana->GetIsolationCut();
      ConfigureIsolationCut(ic,partInCone,thresType,cone,coneMin,pth,col,debug);
    }
  }
  else
    ana->SwitchOffOwnMix();
  
  ana->SetNZvertBin(20);
  
  // Acceptance cut on trigger particle
  // Avoid borders of calorimeter, same as for isolation
  //
  ana->SwitchOnFiducialCut(); 
  // In case of comparison with and without isolation set the same acceptance in both cases
  if ( kAnaCutsString.Contains("Isolation") )
    CalorimeterFiducialCutForIsolationAnalysis(ana->GetFiducialCut(), calorimeter, year, partInCone);
  else 
    CalorimeterFiducialCut(ana->GetFiducialCut(), calorimeter, year);
    
  // Track matching, in case of mixing and specially in isolation case
  //
  AliCaloPID* caloPID = ana->GetCaloPID();
  
  // pT track dependent cuts
  if(tm > 1) caloPID->SwitchOnEMCTrackPtDepResMatching();
  
  // Input / output delta AOD settings
  //
  ana->SetInputAODName(Form("%sTrigger_%s",particle.Data(),kAnaCaloTrackCorr.Data()));
  ana->SetAODObjArrayName(Form("%sHadronCorrIso%dTrigger_%s",particle.Data(),bIsolated,kAnaCaloTrackCorr.Data()));
  //ana->SetAODNamepTInConeHisto(Form("IC%s_%s_R%1.1f_ThMin%1.1f"           ,particle.Data(),kAnaCaloTrackCorr.Data(),cone,pth));
  
  //Set Histograms name tag, bins and ranges
  //
  //ana->AddToHistogramsName(Form("Ana%sHadronCorr_Iso%d_TM%d_",particle.Data(),bIsolated,tm));
  ana->AddToHistogramsName(Form("Ana%sHadronCorr_Iso%d_",particle.Data(),bIsolated));
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  if ( particle == "Random")
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, 100, 1) ;
  
  return ana;
}

///
/// Configure the task doing the selected tracks checking
///
/// \param simulation : A bool identifying the data as simulation
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaChargedParticles* ConfigureChargedAnalysis
( Bool_t simulation, Bool_t printSettings, Int_t   debug, TString histoString  )
{
  AliAnaChargedParticles *ana = new AliAnaChargedParticles();
  
  ana->SetDebug(debug);
  
  // selection cuts
  
  ana->SetMinPt(0.2);
  ana->SwitchOnFiducialCut();
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360) ; //more restrictive cut in reader and after in isolation
  
  ana->SwitchOffFillVertexBC0Histograms();
  
  ana->SwitchOffFillPileUpHistograms();
  ana->SwitchOffFillTrackBCHistograms();
  
  // Branch AOD settings
  
  ana->SetOutputAODName(Form("HadronTrigger_%s",kAnaCaloTrackCorr.Data()));
  ana->SetOutputAODClassName("AliCaloTrackParticle"); // use if no correlation done
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaHadrons_");
  
  SetAnalysisCommonParameters(ana,histoString,"CTS",2012,"pp",simulation,printSettings,debug); // see method below
  
  return ana;
}

///
/// Configure the task doing cluster shape studies
///
/// \param tm : Index with track matching option (0- no TM; 1-Fixed residuals cut; 2-Track pT dependent residuals)
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaClusterShapeCorrelStudies* ConfigureClusterShape
(Int_t tm, TString col , Bool_t  simulation, 
 TString calorimeter, Int_t year, 
 Bool_t  printSettings, Int_t   debug, TString histoString )
{
  AliAnaClusterShapeCorrelStudies *ana = new AliAnaClusterShapeCorrelStudies();
  
  ana->SetM02Min(-1);
  ana->SetNCellsPerClusterMin(-1);
  
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  if(simulation) ana->SetConstantTimeShift(615);
  
  ana->SwitchOffFiducialCut();
  
  ana->SwitchOnStudyClusterShape();
  
  ana->SwitchOnStudyEMCalModuleCells();
  
  ana->SwitchOnStudyColRowFromCellMax() ;

  ana->SwitchOffStudyClusterShapeParam();
  
  ana->SwitchOffStudyMatchedPID() ;
  
  ana->SwitchOffStudyWeight();
  
  ana->SetNCellBinLimits(3); // set to -1 for no analysis on predefined bins in nCell
  ana->SetDistToBadMin(2);
  
  ana->SwitchOffStudyTCardCorrelation() ;
  ana->SwitchOffStudyExotic();
  ana->SwitchOffStudyInvariantMass();
  ana->SwitchOffStudyCellTime() ;
  
  // PID cuts (Track-matching)
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
  
  if(tm > 1)
  {
    // track pT dependent cut
    caloPID->SwitchOnEMCTrackPtDepResMatching();
    
    // Begining of histograms name
    ana->AddToHistogramsName("Shape_TMDep_");
  }
  else if ( tm )
  {
    // Fix
    caloPID->SwitchOffEMCTrackPtDepResMatching();
    
    // Begining of histograms name
    ana->AddToHistogramsName("Shape_TMFix_"); 
  }
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  AliHistogramRanges* histoRanges = ana->GetHistogramRanges();
  histoRanges->SetHistoPtRangeAndNBins(0, 50, 100) ; // Reduce a bit size
  
  return ana;
}

///
/// Configure the task doing standard calorimeter QA
///
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaCalorimeterQA* ConfigureQAAnalysis(TString col,           Bool_t  simulation,
                                         TString calorimeter,   Int_t   year,
                                         Bool_t  printSettings, Int_t   debug, 
                                         TString histoString  )
{
  AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
  
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  ana->SetTimeCut(-1e10,1e10); // Open time cut
  
  ana->SwitchOffCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
    
  ana->SwitchOffFiducialCut();
  ana->SwitchOffFillAllTH3Histogram();
  ana->SwitchOffFillAllPositionHistogram();
  ana->SwitchOffFillAllPositionHistogram2();
  ana->SwitchOffStudyBadClusters() ;
  ana->SwitchOffFillAllCellTimeHisto() ;
  ana->SwitchOffFillAllCellAbsIdHistogram() ;
  
  ana->SwitchOnFillAllTrackMatchingHistogram();
  
  if      ( kAnaCutsString.Contains("QACellsOnly"))
  {
    ana->SwitchOnFillAllCellHistogram() ;
    ana->SwitchOffFillAllClusterHistogram();
  }
  else if ( kAnaCutsString.Contains("QAClustersOnly"))
  {
    ana->SwitchOffFillAllCellHistogram() ;
    ana->SwitchOnFillAllClusterHistogram();
  }
  else
  {
    ana->SwitchOnFillAllCellHistogram() ;
    ana->SwitchOnFillAllClusterHistogram();
  }
  
  ana->AddToHistogramsName("QA_"); // Begining of histograms name
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  return ana;
}

///
/// Configure the task doing exotic cluster/cell checks
///
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaCaloExotics* ConfigureExoticAnalysis(TString col,           Bool_t  simulation,
                                           TString calorimeter,   Int_t   year,
                                           Bool_t  printSettings, Int_t   debug, 
                                           TString histoString  )
{
  AliAnaCaloExotics *ana = new AliAnaCaloExotics();
  
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  ana->SwitchOnFillCellHisto(); 
  ana->SwitchOffFill1CellHisto(); 
  ana->SwitchOffFillMatchingHisto(); 
  
  ana->SetEMinForExo(50);

  ana->SetCellAmpMin(0.2); 
    
  ana->SetTimeCut(-20,20);
  
  if ( simulation ) 
  {
    ana->SetConstantTimeShift(615);
    ana->SwitchOffFillOpenTimeHisto();
    ana->SetTimeCut(-10000,10000);
  }
  
  ana->AddToHistogramsName("Exo_"); // Begining of histograms name
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below
  ana->GetHistogramRanges()->SetHistoPtRangeAndNBins(0, 300, 150) ; // Energy and pt histograms
  ana->GetHistogramRanges()->SetHistoTimeRangeAndNBins(-602.,602.,301);
  ana->GetHistogramRanges()->SetHistoDiffTimeRangeAndNBins(-401, 401, 401);
  return ana;
}

///
/// Configure the task filling generated particle kinematics histograms
///
/// \param partInCone : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
/// \param thresType : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
/// \param pth : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param cone : A float setting the isolation cone size higher limit
/// \param coneMin : A float setting the isolation cone size lower limit
/// \param col : A string with the colliding system
/// \param simulation : A bool identifying the data as simulation
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaGeneratorKine* ConfigureGenKineAnalysis
(Int_t   partInCone,
 Int_t   thresType,     Float_t pth,    
 Float_t cone,          Float_t coneMin,
 TString col,           Bool_t  simulation,
 TString calorimeter,   Int_t   year,
 Bool_t  printSettings, Int_t   debug, 
 TString histoString        )
{
  AliAnaGeneratorKine *ana = new AliAnaGeneratorKine();
  
  // Trigger detector, acceptance and pT cut
  //
  ana->SetTriggerDetector(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetTriggerDetector(calo);
  }

  ana->SetMinPt(2); // Trigger photon, pi0 minimum pT
  
  CalorimeterFiducialCutForIsolationAnalysis(ana->GetFiducialCutForTrigger(), calorimeter, year, partInCone);
  
  // Particles associated to trigger or isolation cone acceptance and pT cut
  //
  ana->SetCalorimeter(calorimeter);
  if(calorimeter == "DCAL") 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  ana->SetMinChargedPt(0.2);
  ana->SetMinNeutralPt(0.5);
  
  ana->SwitchOnFiducialCut(); 
  CalorimeterFiducialCut(ana->GetFiducialCut(), calorimeter, year);
   
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360);
  
  // Isolation paramters
  //
  AliIsolationCut * ic = ana->GetIsolationCut();
  ConfigureIsolationCut(ic,partInCone,thresType,cone,coneMin,pth,col,debug);

  ana->SwitchOnPartonAnalysis();
  if ( kAnaCutsString.Contains("Generator_NoParton") ) 
    ana->SwitchOffPartonAnalysis();
  
  ana->AddToHistogramsName("AnaGenKine_");
  
  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  return ana;
}

///
/// Configure the task doing the trigger particle jet correlation
///
/// \param calorimeter : A string with the calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param year : The year the data was taken, used to configure some histograms
/// \param isolationMatters : A bool identifying whether isolation matters
/// \param gammaConeSize : A float setting the isolation cone size of photon higher limit
/// \param simulation : A bool identifying the data as simulation
/// \param col : A string with the colliding system
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
///
AliAnaParticleJetFinderCorrelation* ConfigureGammaJetAnalysis
(TString calorimeter = "EMCAL"   , Int_t   year = 2011,
 Bool_t  isolationMatters = kTRUE, Float_t gammaConeSize = 0.3,
 Bool_t  simulation = kFALSE     , TString col = "pp",
 Bool_t  printSettings = kFALSE  , Int_t   debug = -1, 
 TString histoString = "")
{
  AliAnaParticleJetFinderCorrelation *ana = new AliAnaParticleJetFinderCorrelation();
  ana->SetDebug(debug);
  //ana->SetDebug(3);
  TString particle="Photon";
  Bool_t isolation=kTRUE;
  if(gammaConeSize>0) isolation=kTRUE;//only isolated photons 
  else isolation=kFALSE;//only not isolated photons 

  Double_t jetConeSize=0.3;
  Double_t ptThresholdInCone=0.15;
  Double_t jetMinPt=0.15,deltaPhi=1.5,minPtRatio=0.,maxPtRatio=5.;
  
  // Input / output delta AOD settings
  //
  //ana->SetInputAODName(Form("%s%s",particle.Data(),kGammaJetCorrelationName.Data()));
  ana->SetInputAODName(Form("%sTrigger_%s",particle.Data(),kAnaCaloTrackCorr.Data()));
  ana->SetAODObjArrayName(Form("%sJetCorr_%s",particle.Data(),kAnaCaloTrackCorr.Data()));

  ana->SwitchOffFiducialCut();

  ana->SetMakeCorrelationInHistoMaker(kFALSE);//must be set to false due to differences AliAODJet and AliEmcalJet

  ana->SetIsolationMatters(isolationMatters); //set whether isolation matters or not
  ana->SelectIsolated(isolation); // do correlation with isolated photons or not
  ana->SetConeSize(jetConeSize); //was 1 - cone to calculate FF
  ana->SetPtThresholdInCone(ptThresholdInCone);
  ana->SetDeltaPhiCutRange(TMath::Pi()-deltaPhi,TMath::Pi()+deltaPhi);  // Delta phi cut for correlation
  ana->SetJetConeSize(jetConeSize);//jet cone size / check the reco jet name
  ana->SetJetMinPt(jetMinPt);//min jet pt
  ana->SetJetAreaFraction(0.8);//min area fraction 
  ana->SetMinPt(0.3);//min cluster pt repeated from reader
  ana->SetGammaConeSize(TMath::Abs(gammaConeSize));//isolation cone repeated from isolation ana
  ana->SetRatioCutRange(minPtRatio,maxPtRatio); // Delta pt cut for correlation

  ana->UseJetRefTracks(kTRUE); //Working now
  //Set Histograms bins and ranges
  //SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below 0,100,200
  //ana->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  ana->SwitchOnBackgroundJetFromReader();
  //background subtraction for photons
  ana->SwitchOffBackgroundSubtractionGamma();//recomended

  //ana->SwitchOnSaveGJTree();
  ana->SwitchOffSaveGJTree();
  ana->SwitchOnMostOpposite();
  //ana->SwitchOnMostEnergetic();

  ana->SwitchOnHistogramTracks();
  ana->SwitchOnHistogramJetTracks();
  ana->SwitchOnHistogramJetBkg();
  
  ana->AddToHistogramsName(Form("Ana%sJetCorr_",particle.Data()));

  SetAnalysisCommonParameters(ana,histoString,calorimeter,year,col,simulation,printSettings,debug); // see method below

  return ana;
}
  

///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param anaList : Container with the sub-analysis to be attached
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param simulation : A bool identifying the data as simulation
/// \param year : The year the data was taken, used to configure some histograms
/// \param col : A string with the colliding system
/// \param analysisString : String that contains what analysis to activate, options: Photon, Isolation, Correlation, ... see below
/// \param histoString : String to add to histo name in case multiple configurations are considered. Very important!!!!
/// \param shshMax : A float setting the maximum value of the shower shape of the clusters for the correlation analysis
/// \param isoCone : A float setting the isolation cone size higher limit
/// \param isoConeMin : A float setting the isolation cone size lower limit
/// \param isoPtTh : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param isoMethod : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
/// \param isoContent : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
/// \param leading : An int setting the type of leading particle selection: 0, select all;l 1: absolute  leading of charged; 2: absolute  leading of charged and neutral; 3: near side leading absolute of charged; 4: near side leading absolute of charged and neutral
/// \param tm : track matching options: 0- no matching; 1-fixed residual cuts; 2-pT track dependent cut
/// \param mixOn : A bool to switch the correlation mixing analysis
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param trigger: Trigger string
///   
///   Options for analysisString:
///    * Analysis: "Photon","InvMass","Electron", "DecayPi0", "MergedPi0", "Charged", "QA", "Isolation", "Correlation", "Generator", "Random","ClusterShape","Exo", "GammaJet"
///    * Isolation analysis: "MultiIsoUESubMethods","MutiIsoR", "MultiIsoRUESubMethods","TightAcc", "FixIsoConeExcess","IsoBandUEGap","IsoBandUEGapFix05"
///    * Common: "SelectEmbed","HighMult","MCRealCaloAcc","PerSM","PerTCard","PerNCells","Bkg"
///                * Track Matching E/P cut: "TMEoP10","TMEoP5",""TMEoP3","TMEoP2","TMEoP1.7","TMEoP1.5"
///    * QA: QACellsOnly, QAClustersOnly
///    * Photon: Recalculation of shower shape in NxN window: "ShSh5x5", "ShSh7x7"
///
void ConfigureCaloTrackCorrAnalysis
(
 TList*   anaList       = 0x0, 
 TString  calorimeter   = "EMCAL", // "DCAL", "PHOS"
 Bool_t   simulation    = kFALSE,
 Int_t    year          = 2011,
 TString  col           = "pp",
 TString  analysisString= "Photon_MergedPi0_DecayPi0_Isolation_Correlation_QA_Charged",
 TString  histoString   = "",
 Float_t  shshMax       = 0.27,
 Float_t  isoCone       = 0.4,
 Float_t  isoConeMin    = -1,
 Float_t  isoPtTh       = 2,
 Int_t    isoMethod     = AliIsolationCut::kSumPtIC,
 Int_t    isoContent    = AliIsolationCut::kNeutralAndCharged,
 Int_t    leading       = 0,
 Int_t    tm            = 2,
 Bool_t   mixOn         = kTRUE,
 Bool_t   printSettings = kFALSE,
 Int_t    debug         = 0,
 const char *triggerStr = "EMC7"
)
{
  if ( !anaList ) 
  {
    printf("ConfigureCaloTrackCorrAnalysis() << No analysis list passed. Stop >>\n");
    return;
  }
  
  // Check the passed variables
  //
  printf("ConfigureCaloTrackCorrAnalysis() << Settings: Base string <%s>, Analysis string <%s>, "
         "\n calorimeter <%s>, simulation <%d>, year <%d>, col <%s>, "
         "\n shshMax <%2.2f>, R <%1.2f>, Rmin <%1.2f>, isoPtTh <%2.2f>, isoMethod <%d>,isoContent <%d>,"
         "\n leading <%d>, tm <%d>, mixOn <%d>, printSettings <%d>, debug <%d>\n",
         anaList->GetName(), analysisString.Data(), 
         calorimeter.Data(), simulation, year, col.Data(),
         shshMax,isoCone,isoConeMin,isoPtTh,isoMethod,isoContent,
         leading,tm,mixOn,printSettings,debug);
  
  kAnaCutsString    = analysisString;
  
  kAnaCaloTrackCorr = Form("%s",anaList->GetName());

  if ( analysisString.Contains("Isolation") )
    kAnaCaloTrackCorr+= Form("_Iso_Meth%d_Part%d_Pt%1.2f_R%1.2f",isoMethod,isoContent,isoPtTh,isoCone);
  if ( analysisString.Contains("Corr") && analysisString.Contains("Photon") && shshMax > 0)
    kAnaCaloTrackCorr+= Form("_CorrM02_%1.2f",shshMax);
  
  if ( isoConeMin > 0 ) kAnaCaloTrackCorr+=Form("_Rmin%1.2f",isoConeMin);
  if ( leading    > 0 ) kAnaCaloTrackCorr+=Form("_Lead%d"   ,leading);
  if ( tm         > 0 ) kAnaCaloTrackCorr+=Form("_TM%d"     ,tm);
 
  printf("ConfigureCaloTrackCorrAnalysis() <<<< TMP branch internal NAME: %s >>>>>\n",kAnaCaloTrackCorr.Data());
  
  TString trigger = triggerStr;

  // #### Configure analysis ####
  
  Int_t n = anaList->GetEntries();//Analysis number, order is important
  
  //
  // Photon analysis
  //
  if ( analysisString.Contains("Photon") )
  {
    anaList->AddAt(ConfigurePhotonAnalysis
                   (col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); // Photon cluster selection
    
    if ( analysisString.Contains("InvMass") )
      anaList->AddAt(ConfigureInvariantMassAnalysis
                     (col,simulation,calorimeter,/*bothCalo*/kFALSE,
                      year,/*mix*/kTRUE,printSettings,debug,histoString), n++);
    
    if ( analysisString.Contains("DecayPi0") )
    {
      // Pi0 event by event selection, invariant mass and photon tagging from decay
      anaList->AddAt(ConfigurePi0EbEAnalysis
                     ("Pi0"        , AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
                      col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
      
      // Eta event by event selection, invariant mass and photon tagging from decay
      anaList->AddAt(ConfigurePi0EbEAnalysis
                     ("Eta"        , AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
                      col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
      
      // Pi0 out of peak event by event selection, and photon tagging from decay
      anaList->AddAt(ConfigurePi0EbEAnalysis
                     ("Pi0SideBand", AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
                      col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
      
      // Eta out of peak event by event selection, and photon tagging from decay
      anaList->AddAt(ConfigurePi0EbEAnalysis
                     ("EtaSideBand", AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
                      col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
    }
    
    // Photon isolation
    //
    if ( analysisString.Contains("Isolation") )
    {
      if (analysisString.Contains("MultiIsoUESubMethods"))
      {
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Photon", leading, AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kSumBkgSubIC,
                        isoCone,isoConeMin,isoPtTh,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Photon", leading, AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kSumBkgSubEtaBandIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Photon", leading, AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kSumBkgSubPhiBandIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Photon", leading, AliIsolationCut::kOnlyCharged, AliIsolationCut::kSumBkgSubIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Photon", leading, AliIsolationCut::kOnlyCharged, AliIsolationCut::kSumBkgSubEtaBandIC,
                        isoCone,isoConeMin,isoPtTh,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Photon", leading, AliIsolationCut::kOnlyCharged, AliIsolationCut::kSumBkgSubPhiBandIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
      }
      else if (analysisString.Contains("MultiIsoRUESubMethods"))
      {
        //printf("**** MultiIsoRUESub ****\n");
        Int_t nsizes = 4;
        Float_t conesize[] = {0.15,0.2,0.3,0.4};
        for(Int_t isize = 0; isize < nsizes; isize++)
        {
          anaList->AddAt(ConfigureIsolationAnalysis
                         ("Photon",leading,AliIsolationCut::kOnlyCharged, AliIsolationCut::kSumBkgSubIC,//isoContent,isoMethod,
                          conesize[isize],isoConeMin,isoPtTh,
                          col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);

          Int_t isoCont = AliIsolationCut::kNeutralAndCharged;
          if ( calorimeter == "DCAL" )
            isoCont = AliIsolationCut::kOnlyCharged;

          anaList->AddAt(ConfigureIsolationAnalysis
                         ("Photon",leading,isoCont, AliIsolationCut::kSumBkgSubPhiBandIC,
                          conesize[isize],isoConeMin,isoPtTh,
                          col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);

          anaList->AddAt(ConfigureIsolationAnalysis
                         ("Photon",leading,isoCont, AliIsolationCut::kSumBkgSubEtaBandIC,
                          conesize[isize],isoConeMin,isoPtTh,
                          col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);
        }
      }
      else if (analysisString.Contains("MultiIsoR"))
      {
        Int_t nsizes = 6;
        Float_t conesize[] = {0.15,0.2,0.25,0.3,0.35,0.4};
        for(Int_t isize = 0; isize < nsizes; isize++)
        {
          anaList->AddAt(ConfigureIsolationAnalysis
                         ("Photon",leading,isoContent,isoMethod,
                          conesize[isize],isoConeMin,isoPtTh, 
                          col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        }
      }
      else // normal case
      {
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Photon", leading, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
      }
      
    }
    
    // Gamma-hadron correlation
    //
    if( analysisString.Contains("Correlation") )
    {
      if ( !analysisString.Contains("MultiIso") )
      {
        anaList->AddAt(ConfigureHadronCorrelationAnalysis
                       ("Photon", leading, kFALSE, shshMax, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh, mixOn,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
      }
      
      if ( analysisString.Contains("Isolation")  )
      {
        anaList->AddAt(ConfigureHadronCorrelationAnalysis
                       ("Photon", leading, kTRUE,  shshMax, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh, mixOn,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString) , n++); 
      }
    } // correlation
    
    // Gamma-jet correlation
    //
    if ( analysisString.Contains("GammaJet") ) 
    {
      if( analysisString.Contains("Isolation") )
        anaList->AddAt(ConfigureGammaJetAnalysis(calorimeter,year,kTRUE ,isoCone,simulation,col,printSettings,debug,histoString), n++);
      else
        anaList->AddAt(ConfigureGammaJetAnalysis(calorimeter,year,kFALSE,isoCone,simulation,col,printSettings,debug,histoString), n++);
    }// end of gamma-jet correlation
  }
  
  //
  // Merged pi0 analysis
  //
  if ( analysisString.Contains("MergedPi0") )
  {
    // Pi0 event by event selection, cluster splitting
    //
    anaList->AddAt(ConfigurePi0EbEAnalysis
                   ("Pi0", AliAnaPi0EbE::kSSCalo,kTRUE,kTRUE,
                    col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
    
    // Merged pi0 isolation
    //
    if ( analysisString.Contains("Isolation") )
    {
      anaList->AddAt(ConfigureIsolationAnalysis
                     ("Pi0SS", leading, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh,
                      col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);   
      
    }
    
    // Merged Pi0-hadron correlation
    //
    if( analysisString.Contains("Correlation") )
    {
      if ( !analysisString.Contains("MultiIso") )
      {
        anaList->AddAt(ConfigureHadronCorrelationAnalysis
                       ("Pi0SS", leading, kFALSE, shshMax, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh, mixOn,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);  
      }
      
      if ( analysisString.Contains("Isolation") )
      {
        anaList->AddAt(ConfigureHadronCorrelationAnalysis
                       ("Pi0SS", leading, kTRUE,  shshMax, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh, mixOn,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString) , n++); 
      }
    } // correlation
  }
  
  //
  // Random ghost trigger analysis, only for MB data
  //
  if ( analysisString.Contains("Random") && 
      ( trigger.Contains("MB")      ||  trigger.Contains("INT") || 
        trigger.Contains("default") ||  trigger.Contains("Cent")  ) )
  {
    // Random trigger generation over selected acceptance
    //
    anaList->AddAt(ConfigureRandomTriggerAnalysis
                   (calorimeter,year,printSettings,debug,histoString), n++); 
    
    // Random trigger isolation on data 
    //
    if ( analysisString.Contains("Isolation") )
    {
      if (analysisString.Contains("MultiIsoUESubMethods"))
      {
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Random", leading, AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kSumBkgSubIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Random", leading, AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kSumBkgSubEtaBandIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Random", leading, AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kSumBkgSubPhiBandIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Random", leading, AliIsolationCut::kOnlyCharged, AliIsolationCut::kSumBkgSubIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Random", leading, AliIsolationCut::kOnlyCharged, AliIsolationCut::kSumBkgSubEtaBandIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Random", leading, AliIsolationCut::kOnlyCharged, AliIsolationCut::kSumBkgSubPhiBandIC,
                        isoCone,isoConeMin,isoPtTh, 
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);
      }
      else if (analysisString.Contains("MultiIsoRUESubMethods"))
      {
        //printf("**** MultiIsoRUESub ****\n");
        Int_t nsizes = 4;
        Float_t conesize[] = {0.15,0.2,0.3,0.4};
        for(Int_t isize = 0; isize < nsizes; isize++)
        {
          anaList->AddAt(ConfigureIsolationAnalysis
                         ("Random",leading,AliIsolationCut::kOnlyCharged, AliIsolationCut::kSumBkgSubIC,//isoContent,isoMethod,
                          conesize[isize],isoConeMin,isoPtTh,
                          col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);

          Int_t isoCont = AliIsolationCut::kNeutralAndCharged;
          if ( calorimeter == "DCAL" )
            isoCont = AliIsolationCut::kOnlyCharged;

          anaList->AddAt(ConfigureIsolationAnalysis
                         ("Random",leading,isoCont, AliIsolationCut::kSumBkgSubPhiBandIC,
                          conesize[isize],isoConeMin,isoPtTh,
                          col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);

          anaList->AddAt(ConfigureIsolationAnalysis
                         ("Random",leading,isoCont, AliIsolationCut::kSumBkgSubEtaBandIC,
                          conesize[isize],isoConeMin,isoPtTh,
                          col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++);
        }
      }
      else if (analysisString.Contains("MultiIsoR"))
      {
        Int_t   nsizes = 6;
        Float_t conesize[] = {0.15,0.2,0.25,0.3,0.35,0.4};
        for(Int_t isize = 0; isize < nsizes; isize++)
        {
          anaList->AddAt(ConfigureIsolationAnalysis
                         ("Random", leading,isoContent,isoMethod,
                          conesize[isize],isoConeMin,isoPtTh, 
                          col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
        }
      }
      else
      {
        anaList->AddAt(ConfigureIsolationAnalysis
                       ("Random", leading, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); 
      }
      
    }
    
    // Random trigger correlation with data
    //
    if( analysisString.Contains("Correlation") )
    {
      if ( !analysisString.Contains("MultiIso") )
      {
        anaList->AddAt(ConfigureHadronCorrelationAnalysis
                       ("Random", leading, kFALSE, shshMax, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh, mixOn,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString), n++); // Gamma-hadron correlation
      }
      
      if ( analysisString.Contains("Isolation")  )
      {
        anaList->AddAt(ConfigureHadronCorrelationAnalysis
                       ("Random", leading, kTRUE,  shshMax, isoContent,isoMethod,isoCone,isoConeMin,isoPtTh, mixOn,
                        col,simulation,calorimeter,year,tm,printSettings,debug,histoString) , n++); // Isolated gamma hadron correlation
      }
    } // correlation
    
    // Random-jet correlation
    //
    if ( analysisString.Contains("GammaJet") ) 
    {
      if( analysisString.Contains("Isolation") )
        anaList->AddAt(ConfigureGammaJetAnalysis(calorimeter,year,kTRUE ,isoCone,simulation,col,printSettings,debug,histoString), n++);
      else
        anaList->AddAt(ConfigureGammaJetAnalysis(calorimeter,year,kFALSE,isoCone,simulation,col,printSettings,debug,histoString), n++);
    }// end of gamma-jet correlation
  }
  
  // Check the generated kinematics
  if(simulation && analysisString.Contains("Generator")) 
  {
    anaList->AddAt(ConfigureGenKineAnalysis
                   (isoContent,isoMethod,isoPtTh,isoCone,isoConeMin,
                    col,simulation,calorimeter,year,printSettings,debug,histoString), n++);
  }
  
  // Charged analysis
  if ( analysisString.Contains("Charged") ) 
  { 
    anaList->AddAt(ConfigureChargedAnalysis(simulation,printSettings,debug,histoString), n++); // track selection checks
  }
 
  // Charged analysis
  if ( analysisString.Contains("Electron") ) 
  { 
    anaList->AddAt(ConfigureElectronAnalysis
                   (col,simulation,calorimeter,year,printSettings,debug,histoString), n++); // electron/hadron selection checks
  }
  
  // Cluster Shape studies
  if ( analysisString.Contains("ClusterShape") ) 
  { 
    anaList->AddAt(ConfigureClusterShape(tm,col,simulation,calorimeter,year,printSettings,debug,histoString) , n++);
  }
  
  // Calo QA
  if ( analysisString.Contains("QA") ) 
  { 
    anaList->AddAt(ConfigureQAAnalysis(col,simulation,calorimeter,year,printSettings,debug,histoString) , n++);
  }
  
  // Exotics
  if ( analysisString.Contains("Exo") ) 
  { 
    anaList->AddAt(ConfigureExoticAnalysis(col,simulation,calorimeter,year,printSettings,debug,histoString) , n++);
  }

  printf("ConfigureCaloTrackCorrAnalysis() << End configuration for %s with total analysis %d>>\n",
         kAnaCaloTrackCorr.Data(), n);
}


