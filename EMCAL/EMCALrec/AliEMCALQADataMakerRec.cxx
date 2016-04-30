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

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h> 
#include <TH3I.h>
#include <TLine.h>
#include <TText.h>
#include <TProfile.h> 
#include <TProfile2D.h> 
#include <TStyle.h>

// --- AliRoot header files ---
#include "AliDAQ.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliEMCALQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliEMCALDigit.h" 
#include "AliEMCALRecPoint.h" 
#include "AliEMCALRawUtils.h"
#include "AliEMCALReconstructor.h"
#include "AliEMCALRecParam.h"
#include "AliRawReader.h"
#include "AliCaloRawStreamV3.h"
#include "AliEMCALGeoParams.h"
#include "AliRawEventHeaderBase.h"
#include "AliQAManager.h"
#include "AliCDBEntry.h"

#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliCaloRawAnalyzer.h"
#include "AliCaloRawAnalyzerFactory.h"

#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerSTURawStream.h"

using namespace std;

/// \cond CLASSIMP
ClassImp(AliEMCALQADataMakerRec) ;
/// \endcond

//____________________________________________________________________________ 
/// Constructor.
///
/// \param fitAlgo: an Integer value to specify the fitting Algorithm to be used
///                 see EMCALUtils/AliCaloConstants.h enum fitAlgorithm
/// 
AliEMCALQADataMakerRec::AliEMCALQADataMakerRec(Int_t fitAlgo) :
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kEMCAL), "EMCAL Quality Assurance Data Maker"),
  fFittingAlgorithm(0),
  fRawAnalyzer(0),
  fRawAnalyzerTRU(0),
  fGeom(0),
  fSuperModules(20), // number of SuperModules; updated to 20 for EMCal + DCal 
  fFirstPedestalSample(0),
  fLastPedestalSample(3),
  fFirstPedestalSampleTRU(0),
  fLastPedestalSampleTRU(3),
  fMinSignalLG(0),
  fMaxSignalLG(AliEMCALGeoParams::fgkSampleMax),
  fMinSignalHG(0),
  fMaxSignalHG(AliEMCALGeoParams::fgkSampleMax),
  fMinSignalTRU(0),
  fMaxSignalTRU(AliEMCALGeoParams::fgkSampleMax),
  fMinSignalLGLEDMon(0),
  fMaxSignalLGLEDMon(AliEMCALGeoParams::fgkSampleMax),
  fMinSignalHGLEDMon(0),
  fMaxSignalHGLEDMon(AliEMCALGeoParams::fgkSampleMax),
  fCalibRefHistoPro(NULL),
  fCalibRefHistoH2F(NULL),
  fLEDMonRefHistoPro(NULL),
  fHighEmcHistoH2F(NULL)
//  fTextSM(new TText*[fSuperModules]) ,
//  fLineCol(NULL),
//  fLineRow(NULL)
{
  SetFittingAlgorithm(fitAlgo);
  
  fGeom = AliEMCALGeometry::GetInstance();
  
  if(!fGeom)
  {
    AliCDBManager* man = AliCDBManager::Instance();
    Int_t runNumber = man->GetRun();
    fGeom =  AliEMCALGeometry::GetInstanceFromRunNumber(runNumber);
  }
  
  if(!fGeom) 
  {
    AliWarning(Form("Using default geometry in reconstruction!!!"));
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  }
  
  if ( !fGeom ) AliFatal(Form("Could not get geometry!"));
  else          AliInfo (Form("Geometry name: <<%s>>",fGeom->GetName())); 
  
  //  for (Int_t sm = 0 ; sm < fSuperModules ; sm++){
  //    fTextSM[sm] = NULL ;
  //  }
}

//____________________________________________________________________________
///
/// Copy constructor.
///
/// \param qadm: AliEMCALQADataMakerRec
///
AliEMCALQADataMakerRec::AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) :
  AliQADataMakerRec(), 
  fFittingAlgorithm(0),
  fRawAnalyzer(0),
  fRawAnalyzerTRU(0),
  fGeom(0),
  fSuperModules(qadm.GetSuperModules()), 
  fFirstPedestalSample(qadm.GetFirstPedestalSample()), 
  fLastPedestalSample(qadm.GetLastPedestalSample()),  
  fFirstPedestalSampleTRU(qadm.GetFirstPedestalSampleTRU()), 
  fLastPedestalSampleTRU(qadm.GetLastPedestalSampleTRU()),  
  fMinSignalLG(qadm.GetMinSignalLG()),
  fMaxSignalLG(qadm.GetMaxSignalLG()),
  fMinSignalHG(qadm.GetMinSignalHG()),
  fMaxSignalHG(qadm.GetMaxSignalHG()),
  fMinSignalTRU(qadm.GetMinSignalTRU()),
  fMaxSignalTRU(qadm.GetMaxSignalTRU()),
  fMinSignalLGLEDMon(qadm.GetMinSignalLGLEDMon()),
  fMaxSignalLGLEDMon(qadm.GetMaxSignalLGLEDMon()),
  fMinSignalHGLEDMon(qadm.GetMinSignalHGLEDMon()),
  fMaxSignalHGLEDMon(qadm.GetMaxSignalHGLEDMon()),
  fCalibRefHistoPro(NULL),
  fCalibRefHistoH2F(NULL),
  fLEDMonRefHistoPro(NULL),
  fHighEmcHistoH2F(NULL)
//  fTextSM(new TText*[fSuperModules]) ,
//  fLineCol(NULL),
//  fLineRow(NULL)
{
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
  SetFittingAlgorithm(qadm.GetFittingAlgorithm());
  
  //  for (Int_t sm = 0 ; sm < fSuperModules ; sm++){
  //    fTextSM[sm] = qadm.fTextSM[sm] ;
  //  }  
}

//__________________________________________________________________
///
/// Assignment operator.
///
/// \param qadm: AliEMCALQADataMakerRec
///
AliEMCALQADataMakerRec& AliEMCALQADataMakerRec::operator = (const AliEMCALQADataMakerRec& qadm )
{
  this->~AliEMCALQADataMakerRec();
  new(this) AliEMCALQADataMakerRec(qadm);
  //  fLineCol = NULL;
  //  fLineRow = NULL;
  //  for (Int_t sm = 0 ; sm < fSuperModules ; sm++){
  //    fTextSM[sm] = qadm.fTextSM[sm] ;
  //  }    
  return *this;
}
 
//____________________________________________________________________________ 
///
/// Detector specific actions at end of cycle
///
/// \param task
/// \param list of histograms
///
void AliEMCALQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //TH1I *rmsL0=(TH1I*)list[fEventSpecie]->At(kNL0TRUSamples);
  //  if(fCycleCounter)
  //	  GetRawsData(kNEventsPerTower)->Scale(1./fCycleCounter);
  /*for(int i=0; i<AliEMCALTriggerMappingV2::fSTURegionNEta; i++)
   for(int j=0; j<AliEMCALTriggerMappingV2::fSTURegionNPhi; j++)
			FillRawsData(kNL0TRURMS,i,j,rmsL0->GetBinError(i,j));
   */
  // reset triggers list to select all histos
  ResetEventTrigClasses(); 
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kEMCAL, task, list) ;  
}

//____________________________________________________________________________ 
///
/// Get the reference histogram from OCDB
///
void AliEMCALQADataMakerRec::GetCalibRefFromOCDB()
{
  TString sName1("hHighEmcalRawMaxMinusMin") ;
  TString sName2("hLowLEDMonEmcalRawMaxMinusMin") ;
  sName1.Prepend(Form("%s_", AliRecoParam::GetEventSpecieName(AliRecoParam::kCalib))) ; 
  sName2.Prepend(Form("%s_", AliRecoParam::GetEventSpecieName(AliRecoParam::kCalib))) ; 
  
  TString refStorage(AliQAv1::GetQARefStorage()) ;
  if (!refStorage.Contains(AliQAv1::GetLabLocalOCDB()) && !refStorage.Contains(AliQAv1::GetLabAliEnOCDB())) {
    AliFatal(Form("%s is not a valid location for reference data", refStorage.Data())) ; 
  } else {
    AliQAManager* manQA = AliQAManager::QAManager(AliQAv1::kRAWS) ;    
    AliQAv1::SetQARefDataDirName(AliRecoParam::kCalib) ;
    if ( ! manQA->GetLock() ) { 
      manQA->SetDefaultStorage(AliQAv1::GetQARefStorage()) ; 
      manQA->SetSpecificStorage("*", AliQAv1::GetQARefStorage()) ;
      manQA->SetRun(AliCDBManager::Instance()->GetRun()) ; 
      manQA->SetLock() ; 
    }
  
    char * detOCDBDir = Form("%s/%s/%s", GetName(), AliQAv1::GetRefOCDBDirName(), AliQAv1::GetRefDataDirName()) ; 
    AliCDBEntry * entry = manQA->Get(detOCDBDir, manQA->GetRun()) ;
    if (entry) {
      TList * listDetQAD =static_cast<TList *>(entry->GetObject()) ;
      if ( strcmp(listDetQAD->ClassName(), "TList") != 0 ) {
        AliError(Form("Expected a Tlist and found a %s for detector %s", listDetQAD->ClassName(), GetName())) ; 
        listDetQAD = NULL ; 
      }
  
      TObjArray * dirOCDB= NULL ; 
      if ( listDetQAD )
        dirOCDB = static_cast<TObjArray *>(listDetQAD->FindObject(Form("%s/%s", AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), AliRecoParam::GetEventSpecieName(AliRecoParam::kCalib)))) ;       
      if (dirOCDB){
        fCalibRefHistoPro = dynamic_cast<TProfile *>(dirOCDB->FindObject(sName1.Data())) ; 
        fLEDMonRefHistoPro = dynamic_cast<TProfile *>(dirOCDB->FindObject(sName2.Data())) ; 
      }
    }
  }
  
  if(fCalibRefHistoPro && fLEDMonRefHistoPro)
  {
    // Defining histograms binning, each 2D histogram covers all SMs
    Int_t nSMSectors = fSuperModules / 2; // 2 SMs per sector
    
    AliCDBManager* man = AliCDBManager::Instance();
    Int_t runNumber = man->GetRun();
    Int_t nbinsZ, nbinsPhi;
    if (!(fGeom-> GetEMCGeometry()->GetGeoName().Contains("DCAL"))){
      nbinsZ = AliEMCALGeoParams::fgkEMCALCols;
      nbinsPhi = AliEMCALGeoParams::fgkEMCALRows;
    }
    else{
      nbinsZ = AliEMCALTriggerMappingV2::fSTURegionNEta;
      nbinsPhi = AliEMCALTriggerMappingV2::fSTURegionNPhi;
    } 
    if(!fCalibRefHistoH2F)
      fCalibRefHistoH2F =  new TH2F("hCalibRefHisto", "hCalibRefHisto", nbinsZ, -0.5, nbinsZ - 0.5, nbinsPhi, -0.5, nbinsPhi -0.5);
    ConvertProfile2H(fCalibRefHistoPro,fCalibRefHistoH2F) ; 
  } else {
    AliFatal(Form("No reference object with name %s or %s found", sName1.Data(), sName2.Data())) ; 
  }
  
}

//____________________________________________________________________________ 
///
/// Create histograms to controll ESD.
/// Multiple paragraphs are split on multiple lines.
///
void AliEMCALQADataMakerRec::InitESDs()
{
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F * h1 = new TH1F("hESDCaloClusterE",  "ESDs CaloCluster energy in EMCAL;Energy [GeV];Counts",    200, 0., 100.) ; 
  h1->Sumw2() ;
  Add2ESDsList(h1, kESDCaloClusE, !expert, image)  ;                                                     
  
  TH1I * h2 = new TH1I("hESDCaloClusterM", "ESDs CaloCluster multiplicity in EMCAL;# of Clusters;Entries", 100, 0,  100) ; 
  h2->Sumw2() ;
  Add2ESDsList(h2, kESDCaloClusM, !expert, image)  ;
  
  TH1F * h3 = new TH1F("hESDCaloCellA",  "ESDs CaloCell amplitude in EMCAL;Energy [GeV];Counts",    500, 0., 50.) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, kESDCaloCellA, !expert, image)  ;  
  
  TH1I * h4 = new TH1I("hESDCaloCellM", "ESDs CaloCell multiplicity in EMCAL;# of Clusters;Entries", 200, 0,  1000) ; 
  h4->Sumw2() ;
  Add2ESDsList(h4, kESDCaloCellM, !expert, image) ;
  // 
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line	
}

//____________________________________________________________________________ 
/// Create Digits histograms in Digits subdir. 
///
void AliEMCALQADataMakerRec::InitDigits()
{
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I * h0 = new TH1I("hEmcalDigits",    "Digits amplitude distribution in EMCAL;Amplitude [ADC counts];Counts",    500, 0, 500) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0, !expert, image) ;
  TH1I * h1 = new TH1I("hEmcalDigitsMul", "Digits multiplicity distribution in EMCAL;# of Digits;Entries", 200, 0, 2000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, 1, !expert, image) ;
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 
///
/// Create Reconstructed PoInt_ts histograms in RecPoints subdir.
///
void AliEMCALQADataMakerRec::InitRecPoints()
{
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F* h0 = new TH1F("hEMCALRpE","EMCAL RecPoint energies;Energy [GeV];Counts",200, 0.,20.); //GeV
  h0->Sumw2();
  Add2RecPointsList(h0,kRecPE, !expert, image);

  TH1I* h1 = new TH1I("hEMCALRpM","EMCAL RecPoint multiplicities;# of Clusters;Entries",100,0,100);
  h1->Sumw2();
  Add2RecPointsList(h1,kRecPM, !expert, image);

  TH1I* h2 = new TH1I("hEMCALRpDigM","EMCAL RecPoint Digit Multiplicities;# of Digits;Entries",20,0,20);
  h2->Sumw2();
  Add2RecPointsList(h2,kRecPDigM, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
}

//____________________________________________________________________________ 
///
///  Create Raws histograms in Raws subdir.
///
void AliEMCALQADataMakerRec::InitRaws()
{
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  const Option_t *profileOption = "s";
  
  Int_t nTowersPerSM = 2*AliEMCALTriggerMappingV2::fNEta*2*AliEMCALTriggerMappingV2::fNPhi; // number of towers in a SuperModule; 24x48
  Int_t nTot = fSuperModules * nTowersPerSM; // max number of towers in all SuperModules
  
  //Defining histograms binning, each 2D histogram covers all SMs
  Int_t nTotTRUs = AliEMCALTriggerMappingV2::fNTotalTRU;  
  Int_t nSMSectors = fSuperModules / 2; // 2 SMs per sector
  Int_t nbinsZ = 2*AliEMCALTriggerMappingV2::fSTURegionNEta;
  Int_t nbinsPhi = 2*AliEMCALTriggerMappingV2::fSTURegionNPhi;
  
  Int_t nTRUCols = nbinsZ/2; // total TRU columns for 2D TRU histos
  Int_t nTRURows = nbinsPhi/2; // total TRU rows for 2D TRU histos
                               // counter info: number of channels per event (bins are SM index)
  TProfile * h0 = new TProfile("hLowEmcalSupermodules", "Low Gain EMC: # of towers vs SuperMod;SM Id;# of towers",
                               fSuperModules, -0.5, fSuperModules-0.5, profileOption) ;
  Add2RawsList(h0, kNsmodLG, !expert, !image, !saveCorr) ;
  TProfile * h1 = new TProfile("hHighEmcalSupermodules", "High Gain EMC: # of towers vs SuperMod;SM Id;# of towers",  
                               fSuperModules, -0.5, fSuperModules-0.5, profileOption) ; 
  Add2RawsList(h1, kNsmodHG, !expert, !image, !saveCorr) ;
  
  // where did max sample occur? (bins are towers)
  TProfile * h2 = new TProfile("hLowEmcalRawtime", "Low Gain EMC: Time at Max vs towerId;Tower Id;Time [ticks]", 
                               nTot, -0.5, nTot-0.5, profileOption) ;
  Add2RawsList(h2, kTimeLG, expert, !image, !saveCorr) ;
  TProfile * h3 = new TProfile("hHighEmcalRawtime", "High Gain EMC: Time at Max vs towerId;Tower Id;Time [ticks]", 
                               nTot, -0.5, nTot-0.5, profileOption) ;
  Add2RawsList(h3, kTimeHG, expert, !image, !saveCorr) ;
  
  // how much above pedestal was the max sample?  (bins are towers)
  TProfile * h4 = new TProfile("hLowEmcalRawMaxMinusMin", "Low Gain EMC: Max - Min vs towerId;Tower Id;Max-Min [ADC counts]", 
                               nTot, -0.5, nTot-0.5, profileOption) ;
  Add2RawsList(h4, kSigLG, !expert, !image, !saveCorr) ;
  TProfile * h5 = new TProfile("hHighEmcalRawMaxMinusMin", "High Gain EMC: Max - Min vs towerId;Tower Id;Max-Min [ADC counts]",
                               nTot, -0.5, nTot-0.5, profileOption) ;
  Add2RawsList(h5, kSigHG, !expert, !image, !saveCorr) ;
  
  // total counter: channels per event
  TH1I * h6 = new TH1I("hLowNtot", "Low Gain EMC: Total Number of found towers;# of Towers;Counts", 200, 0, nTot) ;
  h6->Sumw2() ;
  Add2RawsList(h6, kNtotLG, !expert, !image, !saveCorr) ;
  TH1I * h7 = new TH1I("hHighNtot", "High Gain EMC: Total Number of found towers;# of Towers;Counts", 200,0, nTot) ;
  h7->Sumw2() ;
  Add2RawsList(h7, kNtotHG, !expert, !image, !saveCorr) ;
  
  // pedestal (bins are towers)
  TProfile * h8 = new TProfile("hLowEmcalRawPed", "Low Gain EMC: Pedestal vs towerId;Tower Id;Pedestal [ADC counts]", 
                               nTot, -0.5, nTot-0.5, profileOption) ;
  Add2RawsList(h8, kPedLG, expert, !image, !saveCorr) ;
  TProfile * h9 = new TProfile("hHighEmcalRawPed", "High Gain EMC: Pedestal vs towerId;Tower Id;Pedestal [ADC counts]",
                               nTot, -0.5, nTot-0.5, profileOption) ;
  Add2RawsList(h9, kPedHG, expert, !image, !saveCorr) ;
  
  
  // now repeat the same for TRU and LEDMon data
  Int_t nTot2x2 = nTotTRUs*AliEMCALTriggerMappingV2::fNModulesInTRU; // max number of TRU channels for all SuperModules
  
  // counter info: number of channels per event (bins are SM index)
  TProfile * hT0 = new TProfile("hTRUEmcalSupermodules", "TRU EMC: # of TRU channels vs SuperMod;SM Id;# of TRU channels",
                                fSuperModules, -0.5, fSuperModules-0.5, profileOption) ;
  Add2RawsList(hT0, kNsmodTRU, expert, !image, !saveCorr) ;
  
  // how much above pedestal was the max sample?  (bins are TRU channels)
  TProfile * hT1 = new TProfile("hTRUEmcalRawMaxMinusMin", "TRU EMC: Max - Min vs 2x2Id;2x2 Id;Max-Min [ADC counts]", 
                                nTot2x2, -0.5, nTot2x2-0.5, profileOption) ;
  Add2RawsList(hT1, kSigTRU, expert, !image, !saveCorr) ;
  
  // total counter: channels per event
  TH1I * hT2 = new TH1I("hTRUNtot", "TRU EMC: Total Number of found TRU channels;# of TRU Channels;Counts", 200, 0, nTot2x2) ;
  hT2->Sumw2() ;
  Add2RawsList(hT2, kNtotTRU, expert, !image, !saveCorr) ;
  
  // L0 trigger hits: # of hits (bins are TRU channels)
  TH2I * hT3 = new TH2I("hTRUEmcalL0hits", "L0 trigger hits: Total number of 2x2 L0 generated",  nTRUCols, -0.5, nTRUCols - 0.5, nTRURows, -0.5, nTRURows-0.5);
  hT3->SetOption("COLZ");
  //hT3->Sumw2();
  Add2RawsList(hT3, kNL0TRU, !expert, image, !saveCorr);
  
  // L0 trigger hits: average time (bins are TRU channels)
  TProfile2D * hT4 = new TProfile2D("hTRUEmcalL0hitsAvgTime", "L0 trigger hits: average time bin", nTRUCols, -0.5, nTRUCols - 0.5, nTRURows, -0.5, nTRURows-0.5, profileOption);
  hT4->SetOption("COLZ");
  Add2RawsList(hT4, kTimeL0TRU, !expert, image, !saveCorr);
  
  // L0 trigger hits: first in the event (bins are TRU channels)
  TH1I * hT5 = new TH1I("hTRUEmcalL0hitsFirst", "L0 trigger hits: First hit in the event", nTot2x2, -0.5, nTot2x2);
  hT5->Sumw2();
  Add2RawsList(hT5, kNL0FirstTRU, expert, !image, !saveCorr);
  
  // L0 trigger hits: average time of first hit in the event (bins are TRU channels)
  TProfile * hT6 = new TProfile("hTRUEmcalL0hitsFirstAvgTime", "L0 trigger hits: average time of first hit", nTot2x2, -0.5, nTot2x2, profileOption); 
  Add2RawsList(hT6, kTimeL0FirstTRU, expert, !image, !saveCorr);
  
  // L0 number of time sample: # of samples over the threshold out of the 13 time sample (bins are TRU channels)
  TH2I * hT7 = new TH2I("hTRUEmcalL0Samples", "L0 fired trigger within a bunch per Patch", nTRUCols, -0.5, nTRUCols-0.5, nTRURows, -0.5, nTRURows-0.5);
  hT7->SetOption("COLZ");
  //hT3->Sumw2();
  Add2RawsList(hT7, kNL0TRUSamples, !expert, image, !saveCorr);
  
  TH1F * hT8 = new TH1F("hRMSOfnTimes", "Dispersion of Time for Fired Triggers per Patch", nTRUCols*nTRURows, -0.5, nTRUCols*nTRURows-0.5);
  hT8->Sumw2();
  Add2RawsList(hT8, kNL0TRURMS, !expert, image, !saveCorr);
  
  // and also LED Mon..
  // LEDMon has both high and low gain channels, just as regular FEE/towers
  Int_t nTotLEDMon = fSuperModules * AliEMCALGeoParams::fgkEMCALLEDRefs; // max number of LEDMon channels for all SuperModules 
  
  // counter info: number of channels per event (bins are SM index)
  TProfile * hL0 = new TProfile("hLowLEDMonEmcalSupermodules", "LowLEDMon Gain EMC: # of strips vs SuperMod;SM Id;# of strips",
                                fSuperModules, -0.5, fSuperModules-0.5, profileOption) ;
  Add2RawsList(hL0, kNsmodLGLEDMon, expert, !image, !saveCorr) ;
  TProfile * hL1 = new TProfile("hHighLEDMonEmcalSupermodules", "HighLEDMon Gain EMC: # of strips vs SuperMod;SM Id;# of strips",  
                                fSuperModules, -0.5, fSuperModules-0.5, profileOption) ; 
  Add2RawsList(hL1, kNsmodHGLEDMon, expert, !image, !saveCorr) ;
  
  // where did max sample occur? (bins are strips)
  TProfile * hL2 = new TProfile("hLowLEDMonEmcalRawtime", "LowLEDMon Gain EMC: Time at Max vs stripId;Strip Id;Time [ticks]", 
                                nTotLEDMon, -0.5, nTotLEDMon-0.5, profileOption) ;
  Add2RawsList(hL2, kTimeLGLEDMon, expert, !image, !saveCorr) ;
  TProfile * hL3 = new TProfile("hHighLEDMonEmcalRawtime", "HighLEDMon Gain EMC: Time at Max vs stripId;Strip Id;Time [ticks]", 
                                nTotLEDMon, -0.5, nTotLEDMon-0.5, profileOption) ;
  Add2RawsList(hL3, kTimeHGLEDMon, expert, !image, !saveCorr) ;
  
  // how much above pedestal was the max sample?  (bins are strips)
  TProfile * hL4 = new TProfile("hLowLEDMonEmcalRawMaxMinusMin", "LowLEDMon Gain EMC: Max - Min vs stripId;Strip Id;Max-Min [ADC counts]", 
                                nTotLEDMon, -0.5, nTotLEDMon-0.5, profileOption) ;
  Add2RawsList(hL4, kSigLGLEDMon, expert, !image, !saveCorr) ;
  TProfile * hL5 = new TProfile("hHighLEDMonEmcalRawMaxMinusMin", "HighLEDMon Gain EMC: Max - Min vs stripId;Strip Id;Max-Min [ADC counts]",
                                nTotLEDMon, -0.5, nTotLEDMon-0.5, profileOption) ;
  Add2RawsList(hL5, kSigHGLEDMon, expert, !image, !saveCorr) ;
  
  // total counter: channels per event
  TH1I * hL6 = new TH1I("hLowLEDMonNtot", "LowLEDMon Gain EMC: Total Number of found strips;# of Strips;Counts", 200, 0, nTotLEDMon) ;
  hL6->Sumw2() ;
  Add2RawsList(hL6, kNtotLGLEDMon, expert, !image, !saveCorr) ;
  TH1I * hL7 = new TH1I("hHighLEDMonNtot", "HighLEDMon Gain EMC: Total Number of found strips;# of Strips;Counts", 200,0, nTotLEDMon) ;
  hL7->Sumw2() ;
  Add2RawsList(hL7, kNtotHGLEDMon, expert, !image, !saveCorr) ;
  
  // pedestal (bins are strips)
  TProfile * hL8 = new TProfile("hLowLEDMonEmcalRawPed", "LowLEDMon Gain EMC: Pedestal vs stripId;Strip Id;Pedestal [ADC counts]", 
                                nTotLEDMon, -0.5, nTotLEDMon-0.5, profileOption) ;
  Add2RawsList(hL8, kPedLGLEDMon, expert, !image, !saveCorr) ;
  TProfile * hL9 = new TProfile("hHighLEDMonEmcalRawPed", "HighLEDMon Gain EMC: Pedestal vs stripId;Strip Id;Pedestal [ADC counts]",
                                nTotLEDMon, -0.5, nTotLEDMon-0.5, profileOption) ;
  Add2RawsList(hL9, kPedHGLEDMon, expert, !image, !saveCorr) ;
  
  // temp 2D amplitude histogram for the current run
  fHighEmcHistoH2F = new TH2F("h2DHighEC2", "High Gain EMC:Max - Min [ADC counts]", nbinsZ, -0.5 , nbinsZ-0.5, nbinsPhi, -0.5, nbinsPhi-0.5);
  fHighEmcHistoH2F->SetDirectory(0) ; // this histo must be memory resident
                                      // add ratio histograms: to comapre the current run with the reference data 
  TH2F * h15 = new TH2F("h2DRatioAmp", "High Gain Ratio to Reference:Amplitude_{current run}/Amplitude_{reference run}", nbinsZ, -0.5 , nbinsZ-0.5, 
                        nbinsPhi, -0.5, nbinsPhi-0.5);
  // settings for display in amore
  h15->SetTitle("Amplitude_{current run}/Amplitude_{reference run}"); 
  h15->SetMaximum(2.0);
  h15->SetMinimum(0.1);
  h15->SetOption("COLZ");
  gStyle->SetOptStat(0);
  Int_t color[] = {4,3,2} ;
  gStyle->SetPalette(3,color);
  h15->GetZaxis()->SetNdivisions(3);
  h15->UseCurrentStyle();
  h15->SetDirectory(0);
  Add2RawsList(h15, k2DRatioAmp, !expert, image, !saveCorr) ;
  
  TH1F * h16 = new TH1F("hRatioDist", "Amplitude_{current run}/Amplitude_{reference run} ratio distribution", nTot, 0., 10000.);
  // h16->SetMinimum(0.1); 
  // h16->SetMaximum(100.);
  gStyle->SetOptStat(0);
  h16->UseCurrentStyle();
  h16->SetDirectory(0);
  Add2RawsList(h16, kRatioDist, !expert, image, !saveCorr) ;
  
  // add two histograms for shifter from the LED monitor system: comapre LED monitor with the reference run
  // to be used for decision whether we need to change reference data
  TH1F * hL10 = new TH1F("hMaxMinusMinLEDMonRatio", "LEDMon amplitude, Ratio to reference run", nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  // settings for display in amore
  hL10->SetTitle("Amplitude_{LEDMon current}/Amplitude_{LEDMon reference}"); 
  hL10->SetMaximum(2.0);
  hL10->SetMinimum(0.1); 
  gStyle->SetOptStat(0);
  hL10->UseCurrentStyle();
  hL10->SetDirectory(0);
  //  hL10->SetOption("E");
  Add2RawsList(hL10, kLEDMonRatio, !expert, image, !saveCorr) ;
  
  TH1F * hL11 = new TH1F("hMaxMinusMinLEDMonRatioDist", "LEDMon amplitude, Ratio distribution", nTotLEDMon, 0, 2);
  // hL11->SetMinimum(0.1) ;
  gStyle->SetOptStat(0);
  hL11->UseCurrentStyle();
  hL11->SetDirectory(0);
  Add2RawsList(hL11, kLEDMonRatioDist, !expert, image, !saveCorr) ;
  
  GetCalibRefFromOCDB();   
  
  // STU histgrams
  
  // histos
  Int_t nSTUCols = nbinsZ/2;
  Int_t nSTURows = nbinsPhi/2;
  // 		kAmpL1, kGL1, kJL1,
  // 		kGL1V0, kJL1V0, kSTUTRU  
  
  TProfile2D *hS0 = new TProfile2D("hL1Amp", "Mean STU signal per Row and Column", nSTUCols, -0.5, nSTUCols-0.5, nSTURows, -0.5, nSTURows-0.5);
  Add2RawsList(hS0, kAmpL1, expert, !image, !saveCorr) ;
  //+5 for better visible error box
  TH2F *hS1 = new TH2F("hL1Gamma", "L1 Gamma patch position (FastOR top-left)", nSTUCols, -0.50, nSTUCols-0.5, nSTURows + 5, -0.5, nSTURows-0.5 + 5); 
  Add2RawsList(hS1, kGL1, !expert, image, !saveCorr) ;
  
  TH2F *hS2 = new TH2F("hL1Jet", "L1 Jet patch position (FastOR top-left)", 12, -0.5, nSTUCols-0.5, 16, 0, nSTURows-0.5);
  Add2RawsList(hS2, kJL1, !expert, image, !saveCorr) ;
  
  TH2I *hS3 = new TH2I("hL1GV0", "L1 Gamma patch amplitude versus V0 signal", 500, 0, 50000, 1500, 0, 1500);
  Add2RawsList(hS3, kGL1V0, expert, image, !saveCorr) ;
  
  TH2I *hS4 = new TH2I("hL1JV0", "L1 Jet patch amplitude versus V0 signal", 500, 0, 50000, 1000, 0, 1000);
  Add2RawsList(hS4, kJL1V0, expert, !image, !saveCorr) ;
  
  TH1I *hS5 = new TH1I("hFrameR","Link between TRU and STU", AliEMCALTriggerMappingV2::fNTotalTRU, 0, AliEMCALTriggerMappingV2::fNTotalTRU);
  Add2RawsList(hS5, kSTUTRU, !expert, image, !saveCorr) ;
  
  hS0->SetOption("COLZ");
  hS1->SetOption("COLZ");
  hS2->SetOption("COLZ");
  hS3->SetOption("COLZ");
  hS4->SetOption("COLZ");
  
  // 
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}


//____________________________________________________________________________
///
/// Make QA data from ESDs.
///
/// \param esd: AliESDEvent
///
void AliEMCALQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  Int_t nTot = 0 ; 
  for ( Int_t index = 0; index < esd->GetNumberOfCaloClusters() ; index++ ) {
    AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
    if( clu->IsEMCAL() ) {
      FillESDsData(kESDCaloClusE,clu->E()) ;
      nTot++ ;
    } 
  }
  FillESDsData(kESDCaloClusM,nTot) ;
  
  // fill calo cells
  AliESDCaloCells* cells = esd->GetEMCALCells();
  FillESDsData(kESDCaloCellM,cells->GetNumberOfCells()) ;
  
  for ( Int_t index = 0; index < cells->GetNumberOfCells() ; index++ ) {
    FillESDsData(kESDCaloCellA,cells->GetAmplitude(index)) ;
  }
  
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
}

//____________________________________________________________________________
///
/// Make the histograms for Raw data.
///
/// \param rawReader: AliRawReader
///
void AliEMCALQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  // Check that all the reference histograms exist before we try to use them - otherwise call InitRaws
  // RS: Attention: the counters are increments after custom modification of eventSpecie
  if (!fCalibRefHistoPro || !fCalibRefHistoH2F || !fLEDMonRefHistoPro || !fHighEmcHistoH2F) {
    InitRaws();
  }
  
  // make sure EMCal was readout during the event
  Int_t emcID = AliDAQ::DetectorID("EMCAL"); // bit 18..
  const UInt_t *detPattern = rawReader->GetDetectorPattern(); 
  UInt_t emcInReadout = ( ((1 << emcID) & detPattern[0]) >> emcID);
  if (! emcInReadout) return; // no poInt_t in looking at this event, if no EMCal data
  
  // setup
  rawReader->Reset() ;
  AliCaloRawStreamV3 in(rawReader,"EMCAL"); 
  rawReader->Select("EMCAL",0,AliDAQ::GetFirstSTUDDL()-1) ; // select EMCAL DDL's 
  
  AliRecoParam::EventSpecie_t saveSpecie = fEventSpecie ;
  if (rawReader->GetType() == AliRawEventHeaderBase::kCalibrationEvent) { 
    SetEventSpecie(AliRecoParam::kCalib) ;	
  }
  const Int_t nTowColsPerSM= 2 * AliEMCALTriggerMappingV2::fNEta;
  const Int_t nTowRowsPerSM= 2 * AliEMCALTriggerMappingV2::fNPhi; // number of rows per SuperModule
  
  const Int_t nTowersPerSM = 2*AliEMCALTriggerMappingV2::fNEta*2*AliEMCALTriggerMappingV2::fNPhi; // number of towers in a SuperModule; 24x48
  const Int_t nRows        = 2*AliEMCALTriggerMappingV2::fNPhi; // number of rows per SuperModule
  const Int_t nStripsPerSM = AliEMCALGeoParams::fgkEMCALLEDRefs; // number of strips per SuperModule 
  const Int_t n2x2PerSM    = AliEMCALTriggerMappingV2::fNTRU * AliEMCALTriggerMappingV2::fNModulesInTRU; // number of TRU 2x2's per SuperModule
  const Int_t n2x2PerTRU   = AliEMCALTriggerMappingV2::fNModulesInTRU;
  const Int_t nTot2x2	   = fSuperModules * n2x2PerSM; // total TRU channel
  
  Int_t maxNumL0PatchesPerTRU, nTotTRUs;
  
  if(fGeom->GetEMCGeometry()->GetGeoName().Contains("DCAL")){
    maxNumL0PatchesPerTRU = (AliEMCALTriggerMappingV2::fNModulesInTRUPhi-1)*(AliEMCALTriggerMappingV2::fNModulesInTRUEta-1); // when 4x4 PatchMode for L0 trigger is ON
    nTotTRUs = AliEMCALTriggerMappingV2::fNTotalTRU; // 52 TRUs : 32 EMCAL + 20  ( 14 + 6 unused ) DCAL
  }
  else{
    maxNumL0PatchesPerTRU = (AliEMCALGeoParams::fgkEMCALTRURows-1)*(AliEMCALGeoParams::fgkEMCALTRUCols -1);
    nTotTRUs = AliEMCALGeoParams::fgkEMCALTRUsPerSM * 10 + 2 ; //32 TRUs ( 3*10 Full_SM + 2 1/3_SM )
  }
  
  //const Int_t nTot2x2 = nTotTRUs* n2x2PerTRU;// total TRU channel
  
  // Set nTRUs array size to maximum possible 52, later it will be filled up to the real maximum
  vector <unsigned int> startBins [AliEMCALTriggerMappingV2::fNTotalTRU][n2x2PerTRU];
  
  // SM counters; decl. should be safe, assuming we don't get more than expected SuperModules..
  Int_t nTotalSMLG[AliEMCALGeoParams::fgkEMCALModules]       = {0};
  Int_t nTotalSMHG[AliEMCALGeoParams::fgkEMCALModules]       = {0};
  Int_t nTotalSMTRU[AliEMCALGeoParams::fgkEMCALModules]      = {0};
  Int_t nTotalSMLGLEDMon[AliEMCALGeoParams::fgkEMCALModules] = {0};
  Int_t nTotalSMHGLEDMon[AliEMCALGeoParams::fgkEMCALModules] = {0};
  
  Int_t rmsForTRUL0[AliEMCALTriggerMappingV2::fSTURegionN]={0};
  
  const Int_t nTRUL0ChannelBits = 10; // used for L0 trigger bits checks
  int firstL0TimeBin = 999;
  int triggers[nTot2x2][24]; // auxiliary array for L0 trigger - TODO remove hardcoded 24
  memset(triggers, 0, sizeof(int) * 24 * nTot2x2);
  
  // The next 2 array must be up to 96 (n2x2PerTRU) to account for the possibility
  // to have 2x2 (single module) L0 trigger patches instead of 4x4 (2x2 modules)
  // Set nTRUs array size to maximum possible 52, later it will be filled up to the real maximum
  vector<Int_t> timeOfFiredPatches[AliEMCALTriggerMappingV2::fNTotalTRU][n2x2PerTRU];
  
  bool is2x2PatchModeActive = false; //default in Run2 is 4x4 Patch for L0 trigger
  Int_t iSM = 0; // SuperModule index 
                 // start loop over input stream  
  while (in.NextDDL()) {
    Int_t iRCU = in.GetDDLNumber() % 2; // RCU0 or RCU1, within SuperModule
    Int_t iDDL = in.GetDDLNumber();
    fRawAnalyzer->SetIsZeroSuppressed( in.GetZeroSupp() ); 
    
    while (in.NextChannel()) {
      Int_t iBranch = in.GetBranch();
      
      iSM = in.GetModule(); // SuperModule
                            //prInt_tf("iSM %d DDL %d", iSM, in.GetDDLNumber()); 
      if (iSM>=0 && iSM<fSuperModules) { // valid module reading
        
        Int_t nsamples = 0;
        vector<AliCaloBunchInfo> bunchlist; 
        while (in.NextBunch()) {
          nsamples += in.GetBunchLength();
          bunchlist.push_back( AliCaloBunchInfo(in.GetStartTimeBin(), in.GetBunchLength(), in.GetSignals() ) );
        }  //Storing of each bunch per channel
        
        if (nsamples > 0) { // this check is needed for when we have zero-supp. on, but not sparse readout
          Float_t time = 0.; 
          Float_t amp  = 0.; 
          // indices for pedestal calc.
          Int_t firstPedSample = 0;
          Int_t lastPedSample  = 0;
          bool isTRUL0IdData   = false;
          
          if (! in.IsTRUData() ) { // high gain, low gain, LED Mon data - all have the same shaper/sampling 
            AliCaloFitResults fitResults = fRawAnalyzer->Evaluate( bunchlist, in.GetAltroCFG1(), in.GetAltroCFG2()); 
            amp  = fitResults.GetAmp();
            time = fitResults.GetTof();	
            firstPedSample = fFirstPedestalSample;
            lastPedSample  = fLastPedestalSample;
          }
          else { // TRU data is special, needs its own analyzer
            AliCaloFitResults fitResults = fRawAnalyzerTRU->Evaluate( bunchlist, in.GetAltroCFG1(), in.GetAltroCFG2()); 
            amp  = fitResults.GetAmp();
            time = fitResults.GetTof();	
            firstPedSample = fFirstPedestalSampleTRU;
            lastPedSample  = fLastPedestalSampleTRU;
            if (in.GetColumn() >= n2x2PerTRU) {
              isTRUL0IdData = true;
            }
          }
          
          // pedestal samples
          Int_t nPed = 0;
          vector<Int_t> pedSamples; 
          
          // select earliest bunch 
          unsigned int bunchIndex = 0;
          unsigned int startBin = bunchlist.at(0).GetStartBin();
          if (bunchlist.size() > 0) {
            for(unsigned int ui=1; ui < bunchlist.size(); ui++ ) {
              if (startBin > bunchlist.at(ui).GetStartBin() ) {
                startBin = bunchlist.at(ui).GetStartBin();
                bunchIndex = ui;
              }
            }
          }//bunch size and select earliest bunch
          
          // check bunch for entries in the pedestal sample range
          Int_t bunchLength = bunchlist.at(bunchIndex).GetLength(); 
          const UShort_t *sig = bunchlist.at(bunchIndex).GetData();
          Int_t timebin = 0;
          
          if (! isTRUL0IdData) { // regular data, can look at pedestals
            for (Int_t i = 0; i<bunchLength; i++) {
              timebin = startBin--;
              if ( firstPedSample<=timebin && timebin<=lastPedSample ) {
                pedSamples.push_back( sig[i] );
                nPed++;
              }	    
            } // i
          }//PEDESTALS FOR REGULARE DATA!
          
          /***************/
          // Storing TRUL0Id Data
          else { // TRU L0 Id Data
                 // which TRU the channel belongs to?
                 //Int_t iTRUId = in.GetModule()*3 + (iRCU*in.GetBranch() + iRCU);    
            Int_t iHWaddress = in.GetHWAddress();
            Int_t iTRUId = fGeom->GetTRUIndexFromOnlineHwAdd(iHWaddress,iRCU,iSM);
            for (Int_t i = 0; i< bunchLength; i++) {
              for( Int_t j = 0; j < nTRUL0ChannelBits; j++ ){
                // check if the bit j is 1
                if( (sig[i] & ( 1 << j )) > 0 ){
                  
                  Int_t iPatchInTRU = (in.GetColumn() - n2x2PerTRU)*nTRUL0ChannelBits+j;
                  
                  if(iPatchInTRU<n2x2PerTRU){
                  	 //printf("iTRU %d\tiPatchinTRU %d\tstartBin %d",iTRUId,iPatchInTRU,startBin)	;
                    timeOfFiredPatches[iTRUId][iPatchInTRU].push_back(startBin);
                    //printf("\tnumber of times a patch has fired: %d\n",timeOfFiredPatches[iTRUId][iPatchInTRU].size());
                    startBins[iTRUId][iPatchInTRU].push_back(startBin);
                    if((int)startBin < firstL0TimeBin) firstL0TimeBin = startBin;
                  }
                  else{
                    //only possible if in.GetColumn()>nFastORinTRU && in.GetColumn()==105 and the 8th bit of all sig[i] is 1
                    if(iPatchInTRU==97)
                      is2x2PatchModeActive=true;
                    else
                      continue;
                  }
                }//check if bit j of sig i is 1
              }//loop over the 10 bits of each sig[i]
              startBin--;
            }//loop over the 13 time samples sig[i]
          }//TRU L0idData
          /***************/ //Go to line 282 for the rest of the L0 code    
          
          
          // fill histograms
          Int_t icol = in.GetColumn();
          Int_t irow = in.GetRow();
          
          Int_t towerId = -1;
          if ( in.IsLowGain() || in.IsHighGain() ) { // regular towers
            
            if (iSM < 10) //FULL EMCAL MODULES
              towerId = iSM*nTowersPerSM + icol*nTowRowsPerSM + irow;
            else if (iSM > 9 && iSM < 12) //1/3 EMCAL = OFFSET_FULL_EMCAL + iSM,irow,icol ->tower_number (nrows=1/3nrows EMCAL FULL)
              towerId = 10*nTowersPerSM + (iSM-10)*nTowColsPerSM*(nTowRowsPerSM/3) + icol*(nTowRowsPerSM/3) + irow;
            else if (iSM > 11 && iSM < 18)//FULL DCAL = OFFSET_FULL_EMCAL+ OFFSET_1/3_EMCAL + iSM,irow,icol ->tower_number
              towerId = (iSM-2)*nTowersPerSM + 2*nTowColsPerSM*(nTowRowsPerSM/3) + icol*nTowRowsPerSM + irow;
            else //1/3 DCAL = OFFSET_FULL_EMCAL+ OFFSET_1/3_EMCAL+ OFFSET_FULL_DCAL+ iSM,irow,icol ->tower_number (nrows=1/3nrows EM(D)CAL FULL)
              towerId = 16*nTowersPerSM + (iSM-16)*nTowColsPerSM*(nTowRowsPerSM/3) + icol*(nTowRowsPerSM/3) + irow;
            
            if ( in.IsLowGain() ) { 
              nTotalSMLG[iSM]++; 
              if ( (amp > fMinSignalLG) && (amp < fMaxSignalLG) ) { 
                FillRawsData(kSigLG,towerId, amp);
                FillRawsData(kTimeLG,towerId, time);
              }
              if (nPed > 0) {
                for (Int_t i=0; i<nPed; i++) {
                  FillRawsData(kPedLG,towerId, pedSamples[i]);
                }
              }
            } // gain==0
            else if ( in.IsHighGain() ) {       	
              nTotalSMHG[iSM]++; 
              if ( (amp > fMinSignalHG) && (amp < fMaxSignalHG) ) { 
                FillRawsData(kSigHG,towerId, amp);
                FillRawsData(kTimeHG,towerId, time);
              } 
              if (nPed > 0) {
                for (Int_t i=0; i<nPed; i++) {
                  FillRawsData(kPedHG,towerId, pedSamples[i]);
                }
              }
            } // gain==1
          } // low or high gain
            // TRU
          else if ( in.IsTRUData() && in.GetColumn()<AliEMCALGeoParams::fgkEMCAL2x2PerTRU) {
            // for TRU data, the mapping class holds the TRU Int_ternal 2x2 number (0..95) in the Column var..
            // Int_t iTRU = (iRCU*in.GetBranch() + iRCU); // TRU0 is from RCU0, TRU1 from RCU1, TRU2 is from branch B on RCU1
            //Int_t iTRU2x2Id = iSM*n2x2PerSM + iTRU*AliEMCALGeoParams::fgkEMCAL2x2PerTRU + in.GetColumn();
            Int_t iHWaddress = in.GetHWAddress();
            Int_t iTRU = fGeom->GetTRUIndexFromOnlineHwAdd(iHWaddress,iRCU,iSM); 
            nTotalSMTRU[iSM]++; 
            Int_t iTRU2x2Id; 
            Bool_t gotAbsFastORId=fGeom->GetAbsFastORIndexFromTRU(iTRU, in.GetColumn(), iTRU2x2Id);
            if(!gotAbsFastORId) 
              continue;
            else{
              nTotalSMTRU[iSM]++; 
              if ( (amp > fMinSignalTRU) && (amp < fMaxSignalTRU) ) { 
                FillRawsData(kSigTRU,iTRU2x2Id, amp);
                // FillRawsData(kTimeTRU,iTRU2x2Id, time);
              }
            }
            // if (nPed > 0) {
            // for (Int_t i=0; i<nPed; i++) {
            // FillRawsData(kPedTRU,iTRU2x2Id, pedSamples[i]);
            // }
            // }
          }//end TRU
           // LED Mon
          else if ( in.IsLEDMonData() ) {
            // for LED Mon data, the mapping class holds the gain info in the Row variable
            // and the Strip number in the Column..
            Int_t gain = in.GetRow();
            Int_t stripId = iSM*nStripsPerSM + in.GetColumn();
            
            if ( gain == 0 ) {
              nTotalSMLGLEDMon[iSM]++;
              if ( (amp > fMinSignalLGLEDMon) && (amp < fMaxSignalLGLEDMon) ) {
                FillRawsData(kSigLGLEDMon,stripId, amp);
                FillRawsData(kTimeLGLEDMon,stripId, time);
              }
              if (nPed > 0) {
                for (Int_t i=0; i<nPed; i++) {
                  FillRawsData(kPedLGLEDMon,stripId, pedSamples[i]);
                }
              }
            } // gain==0
            else if ( gain == 1 ) {
              nTotalSMHGLEDMon[iSM]++;
              if ( (amp > fMinSignalHGLEDMon) && (amp < fMaxSignalHGLEDMon) ) {
                FillRawsData(kSigHGLEDMon,stripId, amp);
                FillRawsData(kTimeHGLEDMon,stripId, time);
              }
              if (nPed > 0) {
                for (Int_t i=0; i<nPed; i++) {
                  FillRawsData(kPedHGLEDMon,stripId, pedSamples[i]);
                }
              }//nPed>0
            } // low or high gain
          } // LEDMon
        } // SM index OK
      } // nsamples>0 check, some data found for this channel; not only trailer/header
    }// end while over channel
  }// end while over DDL's, of input stream
  
  /***************/
  //Filling histograms for TRUL0IdData
  Int_t AbsFastORId=-1;
  Int_t AbsFastORIndexesIn4x4Patch[4]={-1};
  Int_t AbsFastORFirstIndexIn4x4Patch;
  Int_t globTRUCol, globTRURow;
  
  bool worry=false;
  for(int k=0; k < nTotTRUs; k++){
    for(int h=0; h < maxNumL0PatchesPerTRU; h++){
      //printf("ciao: iTRU %d\tiPatchinTRU %d\t number of fired patches: %d\n",k,h,timeOfFiredPatches[k][h].size());
      if(timeOfFiredPatches[k][h].size()==0) continue;
      else{ //Run 1 or 2x2PatchModeON L0Channels = modules in TRU
        Int_t maxTsample=-1, minTsample=20;
        Float_t rms=0.;
        if(is2x2PatchModeActive || !(fGeom->GetEMCGeometry()->GetGeoName().Contains("DCAL"))){
          if(fGeom->GetAbsFastORIndexFromTRU(k, h, AbsFastORId)){
            fGeom->GetPositionInEMCALFromAbsFastORIndex(AbsFastORId, globTRUCol, globTRURow);
            for(int s=0; s<timeOfFiredPatches[k][h].size(); s++){
              FillRawsData(kNL0TRU, globTRUCol, globTRURow);
              FillRawsData(kTimeL0TRU, globTRUCol, globTRURow, timeOfFiredPatches[k][h].at(s));
              
              rms+=timeOfFiredPatches[k][h].at(s)*timeOfFiredPatches[k][h].at(s);
              
              if(timeOfFiredPatches[k][h].at(s)>maxTsample) maxTsample= timeOfFiredPatches[k][h].at(s);
              if(timeOfFiredPatches[k][h].at(s)<minTsample) minTsample= timeOfFiredPatches[k][h].at(s);
              
              triggers[AbsFastORId][startBins[k][h].at(s)] = 1;
            }
            FillRawsData(kNL0TRUSamples, globTRUCol, globTRURow,timeOfFiredPatches[k][h].size());
            //FillRawsData(kNL0TRURMS,globTRURow*48+globTRUCol,maxTsample-minTsample); 
            //printf("rms calculated with  SQRT(sum of times^2 / Ntimes samples) %f\n",TMath::Sqrt(rms/timeOfFiredPatches[k][h].size()));
            FillRawsData(kNL0TRURMS,globTRURow*48+globTRUCol,TMath::Sqrt(rms/timeOfFiredPatches[k][h].size())); 
          }
          else continue;
        }
        else{ //Run2 default L0Channels = 4x4 towers patches = 2x2 modules in TRU
          if(fGeom->GetFastORIndexFromL0Index(k, h, AbsFastORIndexesIn4x4Patch, 4)){
            for(int i=0;i<4;i++)
              //printf("AbsID = %d of the %d module in the %d Patch of the %d TRU\n",AbsFastORIndexesIn4x4Patch[i],i,h,k);
              AbsFastORFirstIndexIn4x4Patch=AbsFastORIndexesIn4x4Patch[0];
            fGeom->GetPositionInEMCALFromAbsFastORIndex(AbsFastORFirstIndexIn4x4Patch, globTRUCol, globTRURow);
            //printf("TRU: %d \t PatchIndex:%d \t number of times for the fired patch: %d\n",k,h,timeOfFiredPatches[k][h].size());
            //printf("First Index %d of the 4 modules making the patch\n", AbsFastORFirstIndexIn4x4Patch);
            for(int s=0; s<timeOfFiredPatches[k][h].size(); s++){
              //printf("%d time the patch has been fired and time: %d\n",s+1,timeOfFiredPatches[k][h].at(s));
              rms+=timeOfFiredPatches[k][h].at(s)*timeOfFiredPatches[k][h].at(s);
              if(timeOfFiredPatches[k][h].at(s)>maxTsample) maxTsample= timeOfFiredPatches[k][h].at(s);
              if(timeOfFiredPatches[k][h].at(s)<minTsample) minTsample= timeOfFiredPatches[k][h].at(s);
              
              FillRawsData(kNL0TRU, globTRUCol, globTRURow);
              FillRawsData(kTimeL0TRU, globTRUCol, globTRURow, timeOfFiredPatches[k][h].at(s));
              triggers[AbsFastORFirstIndexIn4x4Patch][startBins[k][h].at(s)] = 1;
            }//number of times each module in a patch has fired the trigger
            FillRawsData(kNL0TRUSamples, globTRUCol, globTRURow,timeOfFiredPatches[k][h].size());
            //FillRawsData(kNL0TRURMS,globTRURow*48+globTRUCol,maxTsample-minTsample); 
            //printf("rms calculated with  SQRT(sum of times^2 / Ntimes samples) %f\n",TMath::Sqrt(rms/timeOfFiredPatches[k][h].size()));
            // printf("\tfilling the bin %d with this value \n",globTRURow*48+globTRUCol);
            FillRawsData(kNL0TRURMS,globTRURow*48+globTRUCol,TMath::Sqrt(rms/timeOfFiredPatches[k][h].size())); 
          }//position of the patch FOUND!
          else
            continue;
        }//check which kind of patch has been SELECTED
      }//the patch "h" in TRU "k" has been fired
      
    }//loop over the max number of AVAILABLE patches per TRU
  }//loop on nTotTRUs (Run1 or Run2)
  
  /***************/
  // filling some L0 trigger histos
  if( firstL0TimeBin < 999 ){
    for(Int_t i = 0; i < nTot2x2; i++) {	
      if( triggers[i][firstL0TimeBin] > 0 ) {
        // histo->Fill(i,j);
        FillRawsData(kNL0FirstTRU, i);
        FillRawsData(kTimeL0FirstTRU, i, firstL0TimeBin);
      }
    }
  }
  
  // calculate the ratio of the amplitude and fill the histograms, only if the events type is Calib
  // RS: operation on the group of histos kSigHG,k2DRatioAmp,kRatioDist,kLEDMonRatio,kLEDMonRatio,kSigLGLEDMon
  const int hGrp[] = {kSigHG,k2DRatioAmp,kRatioDist,kLEDMonRatio,kLEDMonRatioDist,kSigLGLEDMon};
  if ( rawReader->GetType() == AliRawEventHeaderBase::kCalibrationEvent &&
      CheckCloningConsistency(fRawsQAList, hGrp, sizeof(hGrp)/sizeof(int)) ) {  // RS converting original code to loop over all matching triggers
    int nTrig =IsClonedPerTrigClass(kSigHG,fRawsQAList) ? GetNEventTrigClasses() : 0; // loop over triggers only if histos were cloned
                                                                                      //
    for (int itr=-1;itr<nTrig;itr++) { // start from -1 to acknowledge original histos if they were kept
      TObjArray* trArr = GetMatchingRawsHistosSet(hGrp, sizeof(hGrp)/sizeof(int) ,itr);
      if (!trArr) continue;  // no histos for current trigger
                             //
      Double_t binContent = 0.;
      TProfile* prSigHG      = (TProfile *)trArr->At(0); // kSigHG
      TH1* th2DRatioAmp      = (TH1*) trArr->At(1); // k2DRatioAmp
      TH1* thRatioDist       = (TH1*) trArr->At(2); // kRatioDist
      TH1* thLEDMonRatio     = (TH1*) trArr->At(3); // kLEDMonRatio
      TH1* thLEDMonRatioDist = (TH1*) trArr->At(4); // kLEDMonRatio
      TH1* hSigLGLEDMon      = (TH1*) trArr->At(5); // kSigLGLEDMon
      th2DRatioAmp->Reset("ICE");
      thRatioDist->Reset("ICE");
      thLEDMonRatio->Reset("ICE");
      thLEDMonRatioDist->Reset("ICE");
      th2DRatioAmp->ResetStats();
      thRatioDist->ResetStats();
      thLEDMonRatio->ResetStats();
      thLEDMonRatioDist->ResetStats();
      
      ConvertProfile2H(prSigHG, fHighEmcHistoH2F);  
      // 
      for(Int_t ix = 1; ix <= fHighEmcHistoH2F->GetNbinsX(); ix++) {
        for(Int_t iy = 1; iy <= fHighEmcHistoH2F->GetNbinsY(); iy++) { 
          
          // if(fCalibRefHistoH2F->GetBinContent(ix, iy)){
          
          binContent = fHighEmcHistoH2F->GetBinContent(ix, iy);// /fCalibRefHistoH2F->GetBinContent(ix, iy);
                                                               // }
          
          th2DRatioAmp->SetBinContent(ix, iy, binContent);
          thRatioDist->Fill(binContent);
        }
      } 
      // 
      // Now for LED monitor system, to calculate the ratio as well
      Double_t binError = 0. ;
      // For the binError, we add the relative errors, squared
      Double_t relativeErrorSqr = 0. ;
      // 
      for(int ib = 1; ib <= fLEDMonRefHistoPro->GetNbinsX(); ib++) {
        // 
        if(fLEDMonRefHistoPro->GetBinContent(ib) != 0) {
          binContent = hSigLGLEDMon->GetBinContent(ib) / fLEDMonRefHistoPro->GetBinContent(ib);
          relativeErrorSqr = TMath::Power( (fLEDMonRefHistoPro->GetBinError(ib) / fLEDMonRefHistoPro->GetBinContent(ib)), 2);
          if( hSigLGLEDMon->GetBinContent(ib) != 0) {
            relativeErrorSqr += TMath::Power( (hSigLGLEDMon->GetBinError(ib)/hSigLGLEDMon->GetBinContent(ib)), 2);
          }
        }
        else { // ref. run info is zero
          binContent = -1;
          relativeErrorSqr = 1;
        }
        thLEDMonRatio->SetBinContent(ib, binContent);
        
        binError = sqrt(relativeErrorSqr) * binContent;
        thLEDMonRatio->SetBinError(ib, binError);
        thLEDMonRatioDist->Fill(thLEDMonRatio->GetBinContent(ib));
      }
    } // loop over eventual trigger clones
  } 
  // let's also fill the SM and event counter histograms
  Int_t nTotalHG = 0;
  Int_t nTotalLG = 0;
  Int_t nTotalTRU = 0;
  Int_t nTotalHGLEDMon = 0;
  Int_t nTotalLGLEDMon = 0;
  for (iSM=0; iSM<fSuperModules; iSM++) {  
    nTotalLG += nTotalSMLG[iSM]; 
    nTotalHG += nTotalSMHG[iSM]; 
    nTotalTRU += nTotalSMTRU[iSM]; 
    nTotalLGLEDMon += nTotalSMLGLEDMon[iSM]; 
    nTotalHGLEDMon += nTotalSMHGLEDMon[iSM]; 
    FillRawsData(kNsmodLG,iSM, nTotalSMLG[iSM]); 
    FillRawsData(kNsmodHG,iSM, nTotalSMHG[iSM]); 
    FillRawsData(kNsmodTRU,iSM, nTotalSMTRU[iSM]); 
    FillRawsData(kNsmodLGLEDMon,iSM, nTotalSMLGLEDMon[iSM]); 
    FillRawsData(kNsmodHGLEDMon,iSM, nTotalSMHGLEDMon[iSM]); 
  }
  
  FillRawsData(kNtotLG,nTotalLG);
  FillRawsData(kNtotHG,nTotalHG);
  FillRawsData(kNtotTRU,nTotalTRU);
  FillRawsData(kNtotLGLEDMon,nTotalLGLEDMon);
  FillRawsData(kNtotHGLEDMon,nTotalHGLEDMon);
  
  //Last TRUL0 histogram to be filled with RMSs per channel, so AFTER all Channels hav been read 
  
  
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  SetEventSpecie(saveSpecie) ; 
  
  MakeRawsSTU(rawReader);
  
  // just in case the next rawreader consumer forgets to reset; let's do it here again..
  rawReader->Reset() ;
  return;
}

//____________________________________________________________________________
///
/// Make data from Digits.
///
void AliEMCALQADataMakerRec::MakeDigits()
{
  FillDigitsData(1,fDigitsArray->GetEntriesFast()) ; 
  TIter next(fDigitsArray) ; 
  AliEMCALDigit * digit ; 
  while ( (digit = dynamic_cast<AliEMCALDigit *>(next())) ) {
    FillDigitsData(0, digit->GetAmplitude()) ;
  }  
  //
}

//____________________________________________________________________________
///
/// Make data from Digit Tree.
///
/// \param digitTree: TTree
///
void AliEMCALQADataMakerRec::MakeDigits(TTree * digitTree)
{
  // makes data from Digit Tree
  // RS: Attention: the counters are increments in the MakeDigits()
  if (fDigitsArray) 
    fDigitsArray->Clear("C") ; 
  else
    fDigitsArray = new TClonesArray("AliEMCALDigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("EMCAL") ;
  if ( ! branch ) { AliWarning("EMCAL branch in Digit Tree not found"); return; }
  //
  branch->SetAddress(&fDigitsArray) ;
  branch->GetEntry(0) ; 
  MakeDigits() ; 
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();  
  //  
}

//____________________________________________________________________________
///
/// Make data from RecPoints.
///
/// \param clustersTree: TTree
///
void AliEMCALQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  TBranch *emcbranch = clustersTree->GetBranch("EMCALECARP");
  if (!emcbranch) { 
    AliError("can't get the branch with the EMCAL clusters !");
    return;
  }
  
  TObjArray * emcRecPoints = new TObjArray(100) ;
  emcbranch->SetAddress(&emcRecPoints);
  emcbranch->GetEntry(0);
  
  FillRecPointsData(kRecPM,emcRecPoints->GetEntriesFast()) ; 
  TIter next(emcRecPoints) ; 
  AliEMCALRecPoint * rp ; 
  while ( (rp = dynamic_cast<AliEMCALRecPoint *>(next())) ) {
    FillRecPointsData(kRecPE,rp->GetEnergy()) ;
    FillRecPointsData(kRecPDigM,rp->GetMultiplicity());
  }
  emcRecPoints->Delete();
  delete emcRecPoints;
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
}

//____________________________________________________________________________ 
///
/// Detector specific actions at start of cycle.
/// 
void AliEMCALQADataMakerRec::StartOfDetectorCycle()
{
  
}

//____________________________________________________________________________ 
///
/// Set fitting algorithm and initialize it if this same algorithm was not set before.
///
/// \param fitAlgo: kind of fitting algorithm
///
void AliEMCALQADataMakerRec::SetFittingAlgorithm(Int_t fitAlgo)              
{
  fFittingAlgorithm = fitAlgo; // Not sure we need this
  
  fRawAnalyzer    =  AliCaloRawAnalyzerFactory::CreateAnalyzer(fitAlgo);
  
  //  Init also here the TRU algo, even if it is fixed type.
  fRawAnalyzerTRU = AliCaloRawAnalyzerFactory::CreateAnalyzer(Algo::kFakeAltro);
  fRawAnalyzerTRU->SetFixTau(kTRUE);
  fRawAnalyzerTRU->SetTau(2.5); //  default for TRU shaper
}

//_____________________________________________________________________________________
///
/// Convert TProfile to TH2  
///
/// \param p: TProfile
/// \param histo: TH2 histogram
///
void AliEMCALQADataMakerRec::ConvertProfile2H(TProfile * p, TH2 * histo)
{ 
  //  reset histogram
  histo->Reset("ICE") ; 
  histo->ResetStats(); 
  
  Int_t nbinsProf = p->GetNbinsX();
  
  //  loop through the TProfile p and fill the TH2F histo 
  Double_t binContent = 0;
  Int_t towerNum = 0; //  global tower Id
                      //   i = 0; //  tower Id within SuperModule
  Int_t iSM=0; 
  Int_t iSMSide = 0; // 0=A, 1=C side
  Int_t iSMSector = 0; // 2 SM's per sector  
                       //  indices for 2D plots
  Int_t col2d = 0;
  Int_t row2d = 0;
  Int_t kEMCALRows = 24; // TOWERs  nRows
  Int_t kEMCALCols = 48; // TOWERs  nCols
  Int_t kEMCALTowersPerSM = kEMCALRows * kEMCALCols; // nTOWERS in SM
  
  for (Int_t ibin = 1; ibin <= nbinsProf; ibin++) 
  {
    towerNum = (Int_t) p->GetBinCenter(ibin);
    binContent = p->GetBinContent(ibin);
    
    if(towerNum<11520){
      iSM = towerNum/(kEMCALTowersPerSM);
      col2d = (towerNum/kEMCALRows) % kEMCALCols;
      row2d = towerNum % kEMCALRows;
    }
    else if(towerNum>11519 && towerNum<12288){
      iSM = 10 + (towerNum-11520)/(kEMCALTowersPerSM/3);//for the 10th SM the second part returns 0 (division between integers) 
      col2d = ((towerNum-11520)/(kEMCALRows/3)) % kEMCALCols;
      row2d = (towerNum-11520) % (kEMCALRows/3);
    }
    else if(towerNum>12287 && towerNum<19200){
      iSM = 12 + (towerNum-12288)/(kEMCALTowersPerSM);//for the 12th SM the second part return 0 (division between integers) 
      col2d = ((towerNum-12288)/kEMCALRows) % kEMCALCols;
      row2d = (towerNum-12288) % kEMCALRows;
    }
    else{ //10 EMCAL_FULL + 2 EMCAL_1/3 + 6 DCAL_FULL 
      iSM = 18 + (towerNum-19200)/(kEMCALTowersPerSM/3); //for the 18th SM the second part return 0 (division between integers)
      col2d = ((towerNum-19200)/(kEMCALRows/3)) % kEMCALCols;
      row2d = (towerNum-19200) % (kEMCALRows/3); 
    }
    
    iSMSector = iSM / 2;
    iSMSide = iSM % 2;
    if (iSMSide == 1) // C side, shown to the right
      col2d += kEMCALCols;
    
    if (towerNum>=0 && towerNum<12288)//EMCAL (FULL and 1/3 SM)
	     row2d += iSMSector *kEMCALRows;
    else if (towerNum>12287 && towerNum<19968)//DCAL =  EMCAL Offset (120+8) + Offset in DCAL Sectors (here they continue from 6, not 9) 
	     row2d += 5*kEMCALRows + (kEMCALRows/3) + (iSMSector-6)*kEMCALRows;
    
    histo->SetBinContent(col2d+1, row2d+1, binContent);
  }
} 

/*
//____________________________________________________________________________
///
/// Transform local to global indices.
///
/// \param globRow = global row index 
/// \param globColumn = global column index
/// \param module = number of the module
/// \param ddl = number of the ddl
/// \param branch = number of the branch
/// \param column = number of module column inside the TRU
///
void AliEMCALQADataMakerRec::GetTruChannelPosition( Int_t &globRow, Int_t &globColumn, Int_t module, Int_t ddl, Int_t branch, Int_t column ) const
{ // I THINK THIS WHOLE METHOD SHOULD BE CHANGED BUT NEED THE TRU SCHEME! 
  Int_t mrow;
  Int_t mcol;
  Int_t trow;
  Int_t tcol;
  Int_t drow;
  Int_t rcu;
  // RCU 0 or 1
  rcu = ddl % 2;

  // 12 rows of 2x2s in a module (3 TRUs by 4 rows)
  mrow = (module/2) * 12;
  // 24 columns per module, odd module numbers increased by 24
  mcol = (module%2) * 24;

  // position within TRU coordinates
  tcol = column / 4;
  trow = column % 4;

  // .combine
  if( module%2 == 0 ){   // A side
    // mirror rows
    trow = 3 - trow;

    // TRU in module row addition
    drow = (rcu*branch+rcu) * 4;

  }
  else{   // C side
    // mirror columns
    tcol = 23 - tcol;

    // TRU in module row addition
    drow = (2 - (rcu*branch+rcu)) * 4;
  }

  // output global row/collumn position (0,0 = SMA0, phi = 0, |eta| = max)
  globRow = mrow + drow + trow;
  globColumn = mcol + tcol;
}
 
// ____________________________________________________________________________ 
/// 
/// Create the Raw STU histograms.
/// 
/// \param  AliRawReader
///
void AliEMCALQADataMakerRec::MakeRawsSTU(AliRawReader* rawReader)
{ //  STU specifics
  AliEMCALTriggerSTURawStream* inSTU = new AliEMCALTriggerSTURawStream(rawReader);
  Int_t iEMCALtrig[AliEMCALTriggerMappingV2::fSTURegionNEta][AliEMCALTriggerMappingV2::fSTURegionNPhi]={{0}};
  
  for(int dete=0; dete<2; dete++){
  	if(dete==0){		
	  rawReader->Reset();
	  rawReader->Select("EMCAL",AliDAQ::GetFirstSTUDDL(),AliDAQ::GetFirstSTUDDL());
  	}
  	else{
  	  rawReader->Reset();
  	  rawReader->Select("DCAL",AliDAQ::GetLastSTUDDL(),AliDAQ::GetLastSTUDDL());
  	}
  	printf("detector = %d \n",dete);
  	Int_t h=0; //helper for TRU-STU conversion
  	// L1 segmentation
  	Int_t sizeL1gsubr = 1;
  	Int_t sizeL1gpatch = 2; 
  	Int_t sizeL1jsubr = 4; 
  	
  	if (inSTU->ReadPayLoad()){
          Int_t fw = inSTU->GetFwVersion();
  	  Int_t sizeL1jpatch = 2+(fw >> 16);
  	  printf("Firmware= 0x%x \t sizeL1jetpatch= %d \n",fw,sizeL1jpatch);	
      
          Long64_t mask = inSTU->GetFrameReceived() ^ inSTU->GetRegionEnable();
	
          //32 STU for EMCAL and 14 for DCAL
          printf("max number of STu for detector %d = %d\n",dete,(dete==0) ? 32 : 14);
    
          for (int i = 0; i < (dete==0) ? 32 : 14; i++)
          {
            if (!((mask >> i) &  0x1)) FillRawsData(kSTUTRU, i);
          }
 
          // V0 signal in STU
          Int_t iV0Sig = inSTU->GetV0A()+inSTU->GetV0C();
          for (int i = 0; i < (dete==0) ? 32 : 14; i++)
          {
            printf("current STU : %d\n",i);
            UInt_t adc[96];
            for (Int_t j = 0; j < AliEMCALTriggerMappingV2::fNModulesInTRU; j++) adc[j] = 0;
     
            inSTU->GetADC(i, adc);
            Int_t iTRU = fGeom->GetTRUIndexFromSTUIndex(h,dete);
             
            for (Int_t j = 0; j < AliEMCALTriggerMappingV2::fNModulesInTRU; j++)
	    {
	      Int_t idx;
	      fGeom->GetAbsFastORIndexFromTRU(iTRU, j, idx);// CHANGE WITH TRIGGERMAPPINGV2 METHOD
				
	      Int_t px, py;
	      fGeom->GetPositionInEMCALFromAbsFastORIndex(idx, px, py); // CHANGE WITH TRIGGERMAPPINGV2 METHOD
					
	      iEMCALtrig[px][py] = adc[j];
	    }
  	}
  	
        // L1 Gamma patches
        Int_t iTRUSTU, x, y, etaG, phiG;
        printf("Number of High Threshold Gamma patches for Detector %d : %d\n", dete, inSTU->GetNL1GammaPatch(0));
        for(Int_t i = 0; i < inSTU->GetNL1GammaPatch(0); i++)
	 {
	    etaG=0, phiG=0;
	    if (inSTU->GetL1GammaPatch(i, 0, iTRUSTU, x, y)) // col (0..7), row (0..11)
	     {
	       Int_t iTRU,id;
	       	
	       iTRU = fGeom->GetTRUIndexFromSTUIndex(iTRUSTU,dete);
	        
		   if(!fGeom->GetAbsFastORIndexFromPositionInTRU(iTRU, x, y, id))continue;

	       if(!(fGeom->GetPositionInEMCALFromAbsFastORIndex( id, etaG, phiG ))) continue;
	       else{
	      	    // Position of patch L1G (bottom-left FastOR of the patch)
        		//etaG = etaG - sizeL1gsubr * sizeL1gpatch + 1;
      			//phiG = phiG - sizeL1gsubr * sizeL1gpatch + 1;
      			
      			FillRawsData(kGL1, etaG, phiG);
      			Int_t iL1GPatchAmp = 0;
	            for(Int_t L1Gx = 0; L1Gx < sizeL1gpatch; L1Gx ++)
		         {
		           for(Int_t L1Gy = 0; L1Gy < sizeL1gpatch; L1Gy ++)
		             {
		               if (((etaG-L1Gx)> 0 && (etaG+L1Gx) < AliEMCALTriggerMappingV2::fSTURegionNEta) && 
		               		((phiG-L1Gy)>0 && (phiG+L1Gy) < AliEMCALTriggerMappingV2::fSTURegionNPhi)) 
		      	           iL1GPatchAmp += iEMCALtrig[etaG+L1Gx][phiG+L1Gy];
       				 }
		         }
		         FillRawsData(kGL1V0, iV0Sig, iL1GPatchAmp);
	      	   }
	    	}
	  }

	  // L1 Jet patches
	  printf("Number of High Threshold Gamma patches for Detector %d : %d\n", dete, inSTU->GetNL1JetPatch(0));  
	  
          for (Int_t i = 0; i < inSTU->GetNL1JetPatch(0); i++)
	  {
	    if(inSTU->GetL1JetPatch(i, 0, x, y)) // / col (0,15), row (0,11)
	    {
	      //GetPositionInEMCALFrom Something that crosses the TRU boundaries.
	      //NEED To know how JET PATCHES ARE COMPUTED AND STORED.
	      
	      Int_t etaJ = x;
	      Int_t phiJ = y;
	      
	      //Int_t etaJ = sizeL1jsubr * (11-y-sizeL1jpatch + 1); // CHECK THIS FOR JETS
	      //64 EMCAL rows offset when reading DCAL patches 
	      //Int_t phiJ = AliEMCALGeoParams::fgkEMCALSTURows*dete + sizeL1jsubr * (15-x-sizeL1jpatch + 1);
	      
		  // position of patch L1J (FOR bottom-left)
	      FillRawsData(kJL1, x,y);
					
	      // loop the sum aplitude of FOR in the jet patch
	      Int_t iL1JPatchAmp = 0;
	      for (Int_t L1Jx = 0; L1Jx < sizeL1jpatch*4; L1Jx ++)
		  {
		  	for (Int_t L1Jy = 0; L1Jy < sizeL1jpatch*4; L1Jy ++)
		  	  {
		  	    if (etaJ+L1Jx < AliEMCALTriggerMappingV2::fSTURegionNEta && phiJ+L1Jy < AliEMCALTriggerMappingV2::fSTURegionNPhi) 
		      	  iL1JPatchAmp += iEMCALtrig[etaJ+L1Jx][phiJ+L1Jy];
		      }
		  }
		  // cout << "L1J amp =" << iL1JPatchAmp << endl;
                  FillRawsData(kJL1V0, iV0Sig, iL1JPatchAmp);
	   	}//end-if
	  }//end patches loop		
	  delete inSTU;
     }//end inSTU->ReadPayload()  
  }//end detector loop 
      
  // Fill FOR amplitude histo
  for (Int_t i = 0; i < AliEMCALTriggerMappingV2::fSTURegionNEta; i++)
  {
    for (Int_t j = 0; j < AliEMCALTriggerMappingV2::fSTURegionNPhi; j++)
	{
	  if (iEMCALtrig[i][j] != 0) FillRawsData(kAmpL1, i, j, iEMCALtrig[i][j]);
	}
  }
}
*/

// ____________________________________________________________________________
///
/// Create the Raw STU histograms.
/// 
/// \param  rawReader: AliRawReaded
///
void AliEMCALQADataMakerRec::MakeRawsSTU(AliRawReader* rawReader)
{ //  STU specifics
  Int_t h=0; //helper for TRU-STU conversion
             // L1 segmentation 
             // !!!! IMPORTANT THESE ARE in units of MODULES (2*2 towers) !!!!
  Int_t sizeL1gsubr = 1;
  Int_t sizeL1gpatch = 2; 
  Int_t sizeL1jsubr = 2; 
  Int_t sizeL1Jpatch, detPatchOffset, detPhiOffset;
  for(int det=0;det<2;det++){
    //printf("\n\n\n%d times inside the loop for detector",det);
    if(det==0){	
      rawReader->Reset();
      rawReader->Select("EMCAL",AliDAQ::GetFirstSTUDDL(),AliDAQ::GetFirstSTUDDL());
      sizeL1Jpatch=2; //EMCAL 1 JetPatch = 4*4 subregions of 4x4 towers. = 4*4*2
    }
    else{
      rawReader->Reset();
      rawReader->Select("EMCAL",AliDAQ::GetLastSTUDDL(),AliDAQ::GetLastSTUDDL());
      sizeL1Jpatch=1; //DCAL 1 JetPatch = 2*2 subregions of 4x4 towers= 2*2*1  * 4x4 towers.
    }
    detPhiOffset = AliEMCALGeoParams::fgkEMCALSTURows*det;
    detPatchOffset = sizeL1jsubr*sizeL1jsubr*sizeL1Jpatch -1 ;
    //Printf("detector : %c ", (det==0 ? 'E': 'D'));
    //Printf("subregion %d   patch %d   detPatchOffset %d   detphiOffset %d ", sizeL1jsubr,sizeL1Jpatch,detPatchOffset,detPhiOffset);
    AliEMCALTriggerSTURawStream inSTU(rawReader);
    // FOR DCAL
    Int_t iEMCALtrig[AliEMCALTriggerMappingV2::fSTURegionNEta][AliEMCALTriggerMappingV2::fSTURegionNPhi]={{0}};
    //memset(iEMCALtrig, 0, sizeof(int) * AliEMCALTriggerMappingV2::fSTURegionNEta * AliEMCALTriggerMappingV2::fSTURegionNPhi);
    // Int_t iEMCALtrig[AliEMCALGeoParams::fgkEMCALSTUCols][AliEMCALGeoParams::fgkEMCALSTURows];
    // memset(iEMCALtrig, 0, sizeof(int) * AliEMCALGeoParams::fgkEMCALSTUCols * AliEMCALGeoParams::fgkEMCALSTURows);
    
    if (inSTU.ReadPayLoad()) 
    {
      // Fw version (use in case of change in L1 jet 
      Int_t fw = inSTU.GetFwVersion();
      //printf("\nFirmware value: 0x%x\n",fw);
      //printf("Firmware masked 0xf000 = 0x%x\n",fw & 0xf000);
      Int_t sizeL1jpatch = 2+(fw >> 16) ;
      //Printf("sizeL1jetpatch : %d\n",sizeL1jpatch);
      //Int_t sizeL1jpatch=sizeL1jpatchEta;
      //Int_t sizeL1jpatchPhi ;
      
      //to be checked !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //printf("Firmware masked 0xf000 = 0x%x\n",fw & 0xf000);
      /*
       if(fw & 0xf000 == 0xb000){//EMCAL,same dimensions
      	sizeL1jpatchPhi=sizeL1jpatchEta;
       }
       else if(fw & 0xf000 == 0xd000){//DCAL,should be 16*10 //Let's see
      	sizeL1jpatchPhi= 10;
       }
       */
      
      // To check link
      Long64_t mask = inSTU.GetFrameReceived() ^ inSTU.GetRegionEnable();
      
      for (int i = 0; i < AliEMCALTriggerMappingV2::fNTotalTRU; i++)
      {
        if (!((mask >> i) &  0x1)) FillRawsData(kSTUTRU, i);
      }
      
      // V0 signal in STU
      Int_t iV0Sig = inSTU.GetV0A()+inSTU.GetV0C();
      // FastOR amplitude receive from TRU
      
      for (Int_t i = 0; i < (det==0 ? 32 : 14); i++)
      {
        //	  if(i==34 || i==35 || i==40 || i==41 || i==46 || i==47)continue;
        UInt_t adc[96];
        for (Int_t j = 0; j < 96; j++) adc[j] = 0;
        
        /*  if(i>31) {
         det=1;  //Detector switch for STU (for TRUs EMCal and DCal are the same Detector 
         h=i-32; //General Offset DCAL STU / EMCAL+DCAL TRU
         } 
         if(i>35) h=i-34; //1st row DCAL offset STU/TRU 
         if(i>41) h=i-36; //2nd row DCAL offset STU/TRU
         if(i>47) h=i-38; //3rd row DCAL offset STU/TRU
         */
        
        inSTU.GetADC(i, adc);	        
        Int_t iTRU = fGeom->GetTRUIndexFromSTUIndex(i,det);
        
        for (Int_t j = 0; j < 96; j++)
        {
          Int_t idx;
          if(fGeom->GetTriggerMappingVersion() == 1){
            fGeom->GetAbsFastORIndexFromTRU(iTRU, j, idx);
          } else {
            Int_t jTRU = 0, jADC = 0;
            fGeom->GetTRUFromSTU(i, j, jTRU, jADC, det);
            fGeom->GetAbsFastORIndexFromTRU(jTRU, jADC, idx);
          }
          
          Int_t px, py;
          fGeom->GetPositionInEMCALFromAbsFastORIndex(idx, px, py);
          
          iEMCALtrig[px][py] = adc[j];
        }
      }
      
      // L1 Gamma patches
      Int_t iTRUSTU, x, y,etaG,phiG;
      //Printf("Gamma patches");
      //for(detector=0; detector<2;detector++){ // 0 EMCAL , 1 DCAL
      // printf("Number of L1GammaPatches high threshold: %d\n",inSTU.GetNL1GammaPatch(0) );
      for(Int_t ithresh = 0; ithresh < 2; ithresh++)
      {
        for(Int_t i = 0; i < inSTU.GetNL1GammaPatch(ithresh); i++)
        {
          etaG=0,phiG=0;
          if (inSTU.GetL1GammaPatch(i, ithresh, iTRUSTU, x, y)) // col (0..7), row (0..11)
          {
            Int_t iTRU, jTRU, id;

            if(fGeom->GetTriggerMappingVersion() == 1){
              iTRU = fGeom->GetTRUIndexFromSTUIndex(iTRUSTU, det);
              etaG = 23 - x;
              phiG = y + 4 * int(iTRU / 2); // Position in EMCal frame
              if (iTRU % 2) etaG += 24; // C side
              etaG = etaG - 2;
            } else {
              if(!fGeom->GetTRUFromSTU(iTRUSTU, x, y, jTRU, etaG, phiG, det)) continue;
            }
          
            // Position of patch L1G (bottom-left FastOR of the patch)
            //Printf("position of the found patch (eta,phi)=([0,47],[0,63]) for EMCAL ([0,15]&&[32,47],[64,103]) for DCAL\t\t(%d,%d)",etaG,phiG);
            // etaG = (etaG - sizeL1gsubr) * sizeL1gpatch + 1;
            // phiG = (phiG - sizeL1gsubr) * sizeL1gpatch + 1;
            //Printf("position of the found patch (eta,phi) AFTER CORRECTION FOR PATCH size\t\t(%d,%d)",etaG,phiG);
            FillRawsData(kGL1, etaG, phiG);

            // New position in CALORIMETER!
            // if (iTRU<30)
            // Int_t phiG = 11 - y, etaG = x + 8 * int(iTRU/2); // position with new EMCAL TRU configuration !!check iTRU/2..
            // else if ((iTRU>29 && iTRU<32) || (iTRU>43))
            // Int_t etaG = 23 - x, phiG = y + 4 * int(iTRU/2); // position in EMCAL 1/3 SMs !!check iTRU/2
            // else
            // Int_t phiG = 11 - y, etaG = x + 8 * int(iTRU/2); // position in DCAL SMs !!!!Still have to check this!!!
            // Int_t etaG = 23-x, phiG = y + 4 * int(iTRU/2); // position in EMCal
            // if (iTRU%2) etaG += 24; // C-side// / NEED THE NEW TRU NUMBERING SCHEME !!!!!
            // etaG = etaG - sizeL1gsubr * sizeL1gpatch + 1;

            // loop to sum amplitude of FOR in the gamma patch
            Int_t iL1GPatchAmp = 0;
            for(Int_t L1Gx = 0; L1Gx < sizeL1gpatch; L1Gx ++)
            {
              for(Int_t L1Gy = 0; L1Gy < sizeL1gpatch; L1Gy ++)
              {
                if (((etaG+L1Gx) < AliEMCALTriggerMappingV2::fSTURegionNEta) && ((phiG+L1Gy) < AliEMCALTriggerMappingV2::fSTURegionNPhi))
                  iL1GPatchAmp += iEMCALtrig[etaG+L1Gx][phiG+L1Gy];
                // cout << iEMCALtrig[etaG+L1Gx][phiG+L1Gy] << endl;
              }
            }
            // if (iL1GPatchAmp > 500) cout << "L1G amp =" << iL1GPatchAmp << endl;
            FillRawsData(kGL1V0, iV0Sig, iL1GPatchAmp);
          }
        }
      }
      // L1 Jet patches
      //for(detector=0; detector<2;detector++){ // 0 EMCAL , 1 DCAL
      //printf("Number of L1JetPatches high threshold: %d\n",inSTU.GetNL1JetPatch(0) );
      for(Int_t ithresh = 0; ithresh < 2; ++ithresh)
      {
        for (Int_t i = 0; i < inSTU.GetNL1JetPatch(ithresh); i++)
        {
          if(inSTU.GetL1JetPatch(i, ithresh, x, y)) // Position in patches units, should be the other way around (row, col)
          {
            //GetPositionInEMCALFrom Something that crosses the TRU boundaries.
            Int_t etaJ, phiJ;
            if (fGeom->GetTriggerMappingVersion() == 1) {
              etaJ = 11 - y - sizeL1jpatch + 1;
              phiJ = 15 - x - sizeL1jpatch + 1;
            }
            else {
              etaJ = y;
              phiJ = x;
            }

            etaJ *= 4.;
            phiJ *= 4.;
          

            //NEED To know how JET PATCHES ARE COMPUTED AND STORED.
            //printf("sizeL1jsubr %d\n",sizeL1jsubr);
            //Printf("%dth patch in eta , %dth patch in phi",y , x);
            //Int_t etaJ = y*sizeL1jsubr*sizeL1jsubr + detPatchOffset;
            //Int_t etaJ = sizeL1jsubr * (10-y); // CHECK THIS FOR JETS
            //64 EMCAL rows offset when reading DCAL patches
            //Int_t phiJ = detPhiOffset + x*sizeL1jsubr*sizeL1jsubr + detPatchOffset;
            //Printf("eta and phi EMCAL ([7,39],[7,55]) DCAL ([3,43],[67,99]) \t\t %d,%d\n\n",etaJ,phiJ);
            //Int_t phiJ = AliEMCALGeoParams::fgkEMCALSTURows*det + sizeL1jsubr * ((det==0 ? 14 : 8) -x) ;
            //printf("det %d JetPatch %d  x and y for jet patch %d\t%d and etaJ and phiJ  %d \t %d\n",det, i, x,y,etaJ, phiJ);
            // position of patch L1J (FOR bottom-left)
            FillRawsData(kJL1, etaJ +2, phiJ);
          
            // loop the sum aplitude of FOR in the jet patch
            Int_t iL1JPatchAmp = 0;
            for (Int_t L1Jx = 0; L1Jx < sizeL1jpatch*4; L1Jx ++)
            {
              for (Int_t L1Jy = 0; L1Jy < sizeL1jpatch*4; L1Jy ++)
              {
                if (etaJ+L1Jx < AliEMCALTriggerMappingV2::fSTURegionNEta && phiJ+L1Jy < AliEMCALTriggerMappingV2::fSTURegionNPhi)
                  iL1JPatchAmp += iEMCALtrig[etaJ+L1Jx][phiJ+L1Jy];
              }
            }
          
            // cout << "L1J amp =" << iL1JPatchAmp << endl;
            FillRawsData(kJL1V0, iV0Sig, iL1JPatchAmp);
          }//end-if
        }//end patches loop
      }//end detector loop
    }//end inSTU.ReadPayload()
    
    // Fill FOR amplitude histo
    for (Int_t i = 0; i < AliEMCALTriggerMappingV2::fSTURegionNEta; i++)
    {
      for (Int_t j = 0; j < AliEMCALTriggerMappingV2::fSTURegionNPhi; j++)
      {
        if (iEMCALtrig[i][j] != 0) FillRawsData(kAmpL1, i, j, iEMCALtrig[i][j]);
      }
    }
    
  }
}


