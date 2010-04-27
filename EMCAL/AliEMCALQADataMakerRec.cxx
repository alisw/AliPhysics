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
/*
Based on the QA code for PHOS written by Yves Schutz July 2007

Authors:  J.Klay (Cal Poly) May 2008
          S. Salur LBL April 2008
 
Created one histogram for QA shifter;
The idea:average counts for all the towers should be flat 
Change all existing histograms as experts
 --By Yaxian Mao 

*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h> 
#include <TProfile.h> 

// --- Standard library ---


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

#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliCaloRawAnalyzerFastFit.h"
#include "AliCaloRawAnalyzerNN.h"
#include "AliCaloRawAnalyzerLMS.h"
#include "AliCaloRawAnalyzerPeakFinder.h"
#include "AliCaloRawAnalyzerCrude.h"

using namespace std;

ClassImp(AliEMCALQADataMakerRec)
           
//____________________________________________________________________________ 
AliEMCALQADataMakerRec::AliEMCALQADataMakerRec(fitAlgorithm fitAlgo) : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kEMCAL), "EMCAL Quality Assurance Data Maker"),
  fFittingAlgorithm(0),
  fRawAnalyzer(0),
  fRawAnalyzerTRU(0),
  fSuperModules(4), // FIXME!!! number of SuperModules; 4 for 2009; update default to 12 for later runs..
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
  fMaxSignalHGLEDMon(AliEMCALGeoParams::fgkSampleMax)
{
  // ctor
  SetFittingAlgorithm(fitAlgo);
  fRawAnalyzerTRU = new AliCaloRawAnalyzerLMS();
  fRawAnalyzerTRU->SetFixTau(kTRUE); 
  fRawAnalyzerTRU->SetTau(2.5); // default for TRU shaper
}

//____________________________________________________________________________ 
AliEMCALQADataMakerRec::AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) :
  AliQADataMakerRec(), 
  fFittingAlgorithm(0),
  fRawAnalyzer(0),
  fRawAnalyzerTRU(0),
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
  fMaxSignalHGLEDMon(qadm.GetMaxSignalHGLEDMon())
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
  SetFittingAlgorithm(qadm.GetFittingAlgorithm());
  fRawAnalyzerTRU = new AliCaloRawAnalyzerLMS();
  fRawAnalyzerTRU->SetFixTau(kTRUE); 
  fRawAnalyzerTRU->SetTau(2.5); // default for TRU shaper
}

//__________________________________________________________________
AliEMCALQADataMakerRec& AliEMCALQADataMakerRec::operator = (const AliEMCALQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliEMCALQADataMakerRec();
  new(this) AliEMCALQADataMakerRec(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
	
//  if(fCycleCounter)
//	  GetRawsData(kNEventsPerTower)->Scale(1./fCycleCounter);

  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kEMCAL, task, list) ;  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
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
	
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I * h0 = new TH1I("hEmcalDigits",    "Digits amplitude distribution in EMCAL;Amplitude [ADC counts];Counts",    500, 0, 500) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0, !expert, image) ;
  TH1I * h1 = new TH1I("hEmcalDigitsMul", "Digits multiplicity distribution in EMCAL;# of Digits;Entries", 200, 0, 2000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, 1, !expert, image) ;
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
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

}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitRaws()
{
  // create Raws histograms in Raws subdir
   const Bool_t expert   = kTRUE ; 
   const Bool_t saveCorr = kTRUE ; 
   const Bool_t image    = kTRUE ; 

  int nTowersPerSM = AliEMCALGeoParams::fgkEMCALRows * AliEMCALGeoParams::fgkEMCALCols; // number of towers in a SuperModule; 24x48
  int nTot = fSuperModules * nTowersPerSM; // max number of towers in all SuperModules

  // counter info: number of channels per event (bins are SM index)
  TProfile * h0 = new TProfile("hLowEmcalSupermodules", "Low Gain EMC: # of towers vs SuperMod;SM Id;# of towers",
			       fSuperModules, -0.5, fSuperModules-0.5) ;
  Add2RawsList(h0, kNsmodLG, expert, image, !saveCorr) ;
  TProfile * h1 = new TProfile("hHighEmcalSupermodules", "High Gain EMC: # of towers vs SuperMod;SM Id;# of towers",  
			       fSuperModules, -0.5, fSuperModules-0.5) ; 
  Add2RawsList(h1, kNsmodHG, expert, image, !saveCorr) ;

  // where did max sample occur? (bins are towers)
  TProfile * h2 = new TProfile("hLowEmcalRawtime", "Low Gain EMC: Time at Max vs towerId;Tower Id;Time [ticks]", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h2, kTimeLG, expert, image, !saveCorr) ;
  TProfile * h3 = new TProfile("hHighEmcalRawtime", "High Gain EMC: Time at Max vs towerId;Tower Id;Time [ticks]", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h3, kTimeHG, expert, image, !saveCorr) ;

  // how much above pedestal was the max sample?  (bins are towers)
  TProfile * h4 = new TProfile("hLowEmcalRawMaxMinusMin", "Low Gain EMC: Max - Min vs towerId;Tower Id;Max-Min [ADC counts]", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h4, kSigLG, expert, image, !saveCorr) ;
  TProfile * h5 = new TProfile("hHighEmcalRawMaxMinusMin", "High Gain EMC: Max - Min vs towerId;Tower Id;Max-Min [ADC counts]",
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h5, kSigHG, expert, image, !saveCorr) ;

  // total counter: channels per event
  TH1I * h6 = new TH1I("hLowNtot", "Low Gain EMC: Total Number of found towers;# of Towers;Counts", 200, 0, nTot) ;
  h6->Sumw2() ;
  Add2RawsList(h6, kNtotLG, expert, image, !saveCorr) ;
  TH1I * h7 = new TH1I("hHighNtot", "High Gain EMC: Total Number of found towers;# of Towers;Counts", 200,0, nTot) ;
  h7->Sumw2() ;
  Add2RawsList(h7, kNtotHG, expert, image, !saveCorr) ;

  // pedestal (bins are towers)
  TProfile * h8 = new TProfile("hLowEmcalRawPed", "Low Gain EMC: Pedestal vs towerId;Tower Id;Pedestal [ADC counts]", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h8, kPedLG, expert, image, !saveCorr) ;
  TProfile * h9 = new TProfile("hHighEmcalRawPed", "High Gain EMC: Pedestal vs towerId;Tower Id;Pedestal [ADC counts]",
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h9, kPedHG, expert, image, !saveCorr) ;
	
 //number of events per tower, for shifter fast check  	
  TH1I * h12 = new TH1I("hTowerHG", "High Gains on the Tower;Tower", nTot,0, nTot) ;
  h12->Sumw2() ;
  Add2RawsList(h12, kTowerHG, !expert, image, !saveCorr) ;
  TH1I * h13 = new TH1I("hTowerLG", "Low Gains on the Tower;Tower", nTot,0, nTot) ;
  h13->Sumw2() ;
  Add2RawsList(h13, kTowerLG, !expert, image, !saveCorr) ;		

  // now repeat the same for TRU and LEDMon data
  int nTot2x2 = fSuperModules * AliEMCALGeoParams::fgkEMCALTRUsPerSM * AliEMCALGeoParams::fgkEMCAL2x2PerTRU; // max number of TRU channels for all SuperModules

  // counter info: number of channels per event (bins are SM index)
  TProfile * hT0 = new TProfile("hTRUEmcalSupermodules", "TRU EMC: # of TRU channels vs SuperMod;SM Id;# of TRU channels",
				fSuperModules, -0.5, fSuperModules-0.5) ;
  Add2RawsList(hT0, kNsmodTRU, expert, image, !saveCorr) ;

  // where did max sample occur? (bins are TRU channels)
  TProfile * hT1 = new TProfile("hTRUEmcalRawtime", "TRU EMC: Time at Max vs 2x2Id;2x2 Id;Time [ticks]", 
				nTot2x2, -0.5, nTot2x2-0.5) ;
  Add2RawsList(hT1, kTimeTRU, expert, image, !saveCorr) ;

  // how much above pedestal was the max sample?  (bins are TRU channels)
  TProfile * hT2 = new TProfile("hTRUEmcalRawMaxMinusMin", "TRU EMC: Max - Min vs 2x2Id;2x2 Id;Max-Min [ADC counts]", 
				nTot2x2, -0.5, nTot2x2-0.5) ;
  Add2RawsList(hT2, kSigTRU, expert, image, !saveCorr) ;

  // total counter: channels per event
  TH1I * hT3 = new TH1I("hTRUNtot", "TRU EMC: Total Number of found TRU channels;# of TRU Channels;Counts", 200, 0, nTot2x2) ;
  hT3->Sumw2() ;
  Add2RawsList(hT3, kNtotTRU, expert, image, !saveCorr) ;

  // pedestal (bins are TRU channels)
  TProfile * hT4 = new TProfile("hTRUEmcalRawPed", "TRU EMC: Pedestal vs 2x2Id;2x2 Id;Pedestal [ADC counts]", 
				nTot2x2, -0.5, nTot2x2-0.5) ;
  Add2RawsList(hT4, kPedTRU, expert, image, !saveCorr) ;

  // L0 trigger hits: # of hits (bins are TRU channels)
  TH1I * hT5 = new TH1I("hTRUEmcalL0hits", "L0 trigger hits: Total number of 2x2 L0 generated", nTot2x2, -0.5, nTot2x2);
  hT5->Sumw2();
  Add2RawsList(hT5, kNL0TRU, expert, image, !saveCorr);

  // L0 trigger hits: average time (bins are TRU channels)
  TProfile * hT6 = new TProfile("hTRUEmcalL0hitsAvgTime", "L0 trigger hits: average time bin", nTot2x2, -0.5, nTot2x2); 
  Add2RawsList(hT6, kTimeL0TRU, expert, image, !saveCorr);

  // and also LED Mon..
  // LEDMon has both high and low gain channels, just as regular FEE/towers
  int nTotLEDMon = fSuperModules * AliEMCALGeoParams::fgkEMCALLEDRefs; // max number of LEDMon channels for all SuperModules

  // counter info: number of channels per event (bins are SM index)
  TProfile * hL0 = new TProfile("hLowLEDMonEmcalSupermodules", "LowLEDMon Gain EMC: # of strips vs SuperMod;SM Id;# of strips",
			       fSuperModules, -0.5, fSuperModules-0.5) ;
  Add2RawsList(hL0, kNsmodLGLEDMon, expert, image, !saveCorr) ;
  TProfile * hL1 = new TProfile("hHighLEDMonEmcalSupermodules", "HighLEDMon Gain EMC: # of strips vs SuperMod;SM Id;# of strips",  
			       fSuperModules, -0.5, fSuperModules-0.5) ; 
  Add2RawsList(hL1, kNsmodHGLEDMon, expert, image, !saveCorr) ;

  // where did max sample occur? (bins are strips)
  TProfile * hL2 = new TProfile("hLowLEDMonEmcalRawtime", "LowLEDMon Gain EMC: Time at Max vs stripId;Strip Id;Time [ticks]", 
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL2, kTimeLGLEDMon, expert, image, !saveCorr) ;
  TProfile * hL3 = new TProfile("hHighLEDMonEmcalRawtime", "HighLEDMon Gain EMC: Time at Max vs stripId;Strip Id;Time [ticks]", 
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL3, kTimeHGLEDMon, expert, image, !saveCorr) ;

  // how much above pedestal was the max sample?  (bins are strips)
  TProfile * hL4 = new TProfile("hLowLEDMonEmcalRawMaxMinusMin", "LowLEDMon Gain EMC: Max - Min vs stripId;Strip Id;Max-Min [ADC counts]", 
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL4, kSigLGLEDMon, expert, image, !saveCorr) ;
  TProfile * hL5 = new TProfile("hHighLEDMonEmcalRawMaxMinusMin", "HighLEDMon Gain EMC: Max - Min vs stripId;Strip Id;Max-Min [ADC counts]",
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL5, kSigHGLEDMon, expert, image, !saveCorr) ;

  // total counter: channels per event
  TH1I * hL6 = new TH1I("hLowLEDMonNtot", "LowLEDMon Gain EMC: Total Number of found strips;# of Strips;Counts", 200, 0, nTotLEDMon) ;
  hL6->Sumw2() ;
  Add2RawsList(hL6, kNtotLGLEDMon, expert, image, !saveCorr) ;
  TH1I * hL7 = new TH1I("hHighLEDMonNtot", "HighLEDMon Gain EMC: Total Number of found strips;# of Strips;Counts", 200,0, nTotLEDMon) ;
  hL7->Sumw2() ;
  Add2RawsList(hL7, kNtotHGLEDMon, expert, image, !saveCorr) ;

  // pedestal (bins are strips)
  TProfile * hL8 = new TProfile("hLowLEDMonEmcalRawPed", "LowLEDMon Gain EMC: Pedestal vs stripId;Strip Id;Pedestal [ADC counts]", 
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL8, kPedLGLEDMon, expert, image, !saveCorr) ;
  TProfile * hL9 = new TProfile("hHighLEDMonEmcalRawPed", "HighLEDMon Gain EMC: Pedestal vs stripId;Strip Id;Pedestal [ADC counts]",
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL9, kPedHGLEDMon, expert, image, !saveCorr) ;

}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs

  Int_t nTot = 0 ; 
  for ( Int_t index = 0; index < esd->GetNumberOfCaloClusters() ; index++ ) {
    AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
    if( clu->IsEMCAL() ) {
      GetESDsData(kESDCaloClusE)->Fill(clu->E()) ;
      nTot++ ;
    } 
  }
  GetESDsData(kESDCaloClusM)->Fill(nTot) ;

  //fill calo cells
  AliESDCaloCells* cells = esd->GetEMCALCells();
  GetESDsData(kESDCaloCellM)->Fill(cells->GetNumberOfCells()) ;

  for ( Int_t index = 0; index < cells->GetNumberOfCells() ; index++ ) {
    GetESDsData(kESDCaloCellA)->Fill(cells->GetAmplitude(index)) ;
  }

}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //Fill prepared histograms with Raw digit properties

  // make sure EMCal was readout during the event
  Int_t emcID = AliDAQ::DetectorID("EMCAL"); // bit 18..
  const UInt_t *detPattern = rawReader->GetDetectorPattern(); 
  UInt_t emcInReadout = ( ((1 << emcID) & detPattern[0]) >> emcID);
  if (! emcInReadout) return; // no point in looking at this event, if no EMCal data

  // setup
  rawReader->Reset() ;
  AliCaloRawStreamV3 in(rawReader,"EMCAL"); 
  rawReader->Select("EMCAL", 0, AliEMCALGeoParams::fgkLastAltroDDL) ; //select EMCAL DDL's 

  AliRecoParam::EventSpecie_t saveSpecie = fEventSpecie ;

  if (rawReader->GetType() == AliRawEventHeaderBase::kCalibrationEvent) { 
    SetEventSpecie(AliRecoParam::kCalib) ;
  }
  
  fRawAnalyzer->SetIsZeroSuppressed(true); // TMP - should use stream->IsZeroSuppressed(), or altro cfg registers later

  int nTowersPerSM = AliEMCALGeoParams::fgkEMCALRows * AliEMCALGeoParams::fgkEMCALCols; // number of towers in a SuperModule; 24x48
  int nRows = AliEMCALGeoParams::fgkEMCALRows; // number of rows per SuperModule
  int nStripsPerSM = AliEMCALGeoParams::fgkEMCALLEDRefs; // number of strips per SuperModule
  int n2x2PerSM = AliEMCALGeoParams::fgkEMCALTRUsPerSM * AliEMCALGeoParams::fgkEMCAL2x2PerTRU; // number of TRU 2x2's per SuperModule
  int n2x2PerTRU = AliEMCALGeoParams::fgkEMCAL2x2PerTRU;

  // SM counters; decl. should be safe, assuming we don't get more than expected SuperModules..
  int nTotalSMLG[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nTotalSMHG[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nTotalSMTRU[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nTotalSMLGLEDMon[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nTotalSMHGLEDMon[AliEMCALGeoParams::fgkEMCALModules] = {0};

  int nTRUL0ChannelBits = 10; // used for L0 trigger bits checks

  // start loop over input stream  
  int iSM = 0;
  while (in.NextDDL()) {
    int iRCU = in.GetDDLNumber() % 2; // RCU0 or RCU1, within SuperModule

    while (in.NextChannel()) {
      iSM = in.GetModule(); // SuperModule
      //printf("iSM %d DDL %d", iSM, in.GetDDLNumber()); 
      if (iSM>=0 && iSM<fSuperModules) { // valid module reading

	int nsamples = 0;
	vector<AliCaloBunchInfo> bunchlist; 
	while (in.NextBunch()) {
	  nsamples += in.GetBunchLength();
	  bunchlist.push_back( AliCaloBunchInfo(in.GetStartTimeBin(), in.GetBunchLength(), in.GetSignals() ) );
	} 
	
	if (nsamples > 0) { // this check is needed for when we have zero-supp. on, but not sparse readout
	  Float_t time = 0; 
	  Float_t amp = 0; 
	  // indices for pedestal calc.
	  int firstPedSample = 0;
	  int lastPedSample = 0;
	  bool isTRUL0IdData = false;

	  if (! in.IsTRUData() ) { // high gain, low gain, LED Mon data - all have the same shaper/sampling 
	    AliCaloFitResults fitResults = fRawAnalyzer->Evaluate( bunchlist, in.GetAltroCFG1(), in.GetAltroCFG2()); 
	    amp = fitResults.GetAmp();
	    time = fitResults.GetTof();	
	    firstPedSample = fFirstPedestalSample;
	    lastPedSample = fLastPedestalSample;
	  }
	  else { // TRU data is special, needs its own analyzer
	    AliCaloFitResults fitResults = fRawAnalyzerTRU->Evaluate( bunchlist, in.GetAltroCFG1(), in.GetAltroCFG2()); 
	    amp = fitResults.GetAmp();
	    time = fitResults.GetTof();	
	    firstPedSample = fFirstPedestalSampleTRU;
	    lastPedSample = fLastPedestalSampleTRU;
	    if (in.GetColumn() > n2x2PerTRU) {
	      isTRUL0IdData = true;
	    }
	  }
  
	  // pedestal samples
	  int nPed = 0;
	  vector<int> pedSamples; 
	
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
	  }

	  // check bunch for entries in the pedestal sample range
	  int bunchLength = bunchlist.at(bunchIndex).GetLength(); 
	  const UShort_t *sig = bunchlist.at(bunchIndex).GetData();
	  int timebin = 0;

	  if (! isTRUL0IdData) { // regular data, can look at pedestals
	    for (int i = 0; i<bunchLength; i++) {
	      timebin = startBin--;
	      if ( firstPedSample<=timebin && timebin<=lastPedSample ) {
		pedSamples.push_back( sig[i] );
		nPed++;
	      }	    
	    } // i
	  //	  printf("nPed %d\n", nPed);
	  }
	  else { // TRU L0 Id Data
	    // which TRU the channel belongs to?
	    int iTRUId = in.GetModule()*3 + (iRCU*in.GetBranch() + iRCU);

	    for (int i = 0; i< bunchLength; i++) {
	      for( int j = 0; j < nTRUL0ChannelBits; j++ ){
		// check if the bit j is 1
		if( (sig[i] & ( 1 << j )) > 0 ){
		  int iTRUIdInSM = (in.GetColumn() - n2x2PerTRU)*nTRUL0ChannelBits+j;
		  if(iTRUIdInSM < n2x2PerTRU) {
		    int iTRUAbsId = iTRUIdInSM + n2x2PerTRU * iTRUId;
		    // Fill the histograms
		    GetRawsData(kNL0TRU)->Fill(iTRUAbsId);
		    GetRawsData(kTimeL0TRU)->Fill(iTRUAbsId, startBin);
		  }
		}
	      }
	      startBin--;
	    } // i	
	  } // TRU L0 Id data			

	  // fill histograms
	  if ( in.IsLowGain() || in.IsHighGain() ) { // regular towers
	    int towerId = iSM*nTowersPerSM + in.GetColumn()*nRows + in.GetRow();
	    if ( in.IsLowGain() ) { 
	      nTotalSMLG[iSM]++; 
	      GetRawsData(kTowerLG)->Fill(towerId);
	      if ( (amp > fMinSignalLG) && (amp < fMaxSignalLG) ) { 
		GetRawsData(kSigLG)->Fill(towerId, amp);
		GetRawsData(kTimeLG)->Fill(towerId, time);
	      }
	      if (nPed > 0) {
		for (int i=0; i<nPed; i++) {
		  GetRawsData(kPedLG)->Fill(towerId, pedSamples[i]);
		}
	      }
	    } // gain==0
	    else if ( in.IsHighGain() ) {       	
	      nTotalSMHG[iSM]++; 
	      GetRawsData(kTowerHG)->Fill(towerId);
	      if ( (amp > fMinSignalHG) && (amp < fMaxSignalHG) ) { 
		GetRawsData(kSigHG)->Fill(towerId, amp);
		GetRawsData(kTimeHG)->Fill(towerId, time);
	      } 
	      if (nPed > 0) {
		for (int i=0; i<nPed; i++) {
		  GetRawsData(kPedHG)->Fill(towerId, pedSamples[i]);
		}
	      }
	    } // gain==1
	  } // low or high gain
	  // TRU
	  else if ( in.IsTRUData() && in.GetColumn()<AliEMCALGeoParams::fgkEMCAL2x2PerTRU) {
	    // for TRU data, the mapping class holds the TRU internal 2x2 number (0..95) in the Column var..
	    int iTRU = (iRCU*in.GetBranch() + iRCU); //TRU0 is from RCU0, TRU1 from RCU1, TRU2 is from branch B on RCU1
	    int iTRU2x2Id = iSM*n2x2PerSM + iTRU*AliEMCALGeoParams::fgkEMCAL2x2PerTRU 
	      + in.GetColumn();
	    nTotalSMTRU[iSM]++; 
	    if ( (amp > fMinSignalTRU) && (amp < fMaxSignalTRU) ) { 
	      GetRawsData(kSigTRU)->Fill(iTRU2x2Id, amp);
	      GetRawsData(kTimeTRU)->Fill(iTRU2x2Id, time);
	    }
	    if (nPed > 0) {
	      for (int i=0; i<nPed; i++) {
		GetRawsData(kPedTRU)->Fill(iTRU2x2Id, pedSamples[i]);
	      }
	    }
	  }
	  // LED Mon
	  else if ( in.IsLEDMonData() ) {
	    // for LED Mon data, the mapping class holds the gain info in the Row variable
	    // and the Strip number in the Column..
	    int gain = in.GetRow(); 
	    int stripId = iSM*nStripsPerSM + in.GetColumn();
	  
	    if ( gain == 0 ) { 
	      nTotalSMLGLEDMon[iSM]++; 
	      if ( (amp > fMinSignalLGLEDMon) && (amp < fMaxSignalLGLEDMon) ) { 
		GetRawsData(kSigLGLEDMon)->Fill(stripId, amp);
		GetRawsData(kTimeLGLEDMon)->Fill(stripId, time);
	      }
	      if (nPed > 0) {
		for (int i=0; i<nPed; i++) {
		  GetRawsData(kPedLGLEDMon)->Fill(stripId, pedSamples[i]);
		}
	      }
	    } // gain==0
	    else if ( gain == 1 ) {       	
	      nTotalSMHGLEDMon[iSM]++; 
	      if ( (amp > fMinSignalHGLEDMon) && (amp < fMaxSignalHGLEDMon) ) { 
		GetRawsData(kSigHGLEDMon)->Fill(stripId, amp);
		GetRawsData(kTimeHGLEDMon)->Fill(stripId, time);
	      }
	      if (nPed > 0) {
		for (int i=0; i<nPed; i++) {
		  GetRawsData(kPedHGLEDMon)->Fill(stripId, pedSamples[i]);
		}
	      }
	    } // low or high gain
	  } // LEDMon

	} // SM index OK

      } // nsamples>0 check, some data found for this channel; not only trailer/header
    }// end while over channel 
   
  }//end while over DDL's, of input stream 

  // let's also fill the SM and event counter histograms
  int nTotalHG = 0;
  int nTotalLG = 0;
  int nTotalTRU = 0;
  int nTotalHGLEDMon = 0;
  int nTotalLGLEDMon = 0;
  for (iSM=0; iSM<fSuperModules; iSM++) {  
    nTotalLG += nTotalSMLG[iSM]; 
    nTotalHG += nTotalSMHG[iSM]; 
    nTotalTRU += nTotalSMTRU[iSM]; 
    nTotalLGLEDMon += nTotalSMLGLEDMon[iSM]; 
    nTotalHGLEDMon += nTotalSMHGLEDMon[iSM]; 
    GetRawsData(kNsmodLG)->Fill(iSM, nTotalSMLG[iSM]); 
    GetRawsData(kNsmodHG)->Fill(iSM, nTotalSMHG[iSM]); 
    GetRawsData(kNsmodTRU)->Fill(iSM, nTotalSMTRU[iSM]); 
    GetRawsData(kNsmodLGLEDMon)->Fill(iSM, nTotalSMLGLEDMon[iSM]); 
    GetRawsData(kNsmodHGLEDMon)->Fill(iSM, nTotalSMHGLEDMon[iSM]); 
  }
  
  GetRawsData(kNtotLG)->Fill(nTotalLG);
  GetRawsData(kNtotHG)->Fill(nTotalHG);
  GetRawsData(kNtotTRU)->Fill(nTotalTRU);
  GetRawsData(kNtotLGLEDMon)->Fill(nTotalLGLEDMon);
  GetRawsData(kNtotHGLEDMon)->Fill(nTotalHGLEDMon);
 

  SetEventSpecie(saveSpecie) ; 
  // just in case the next rawreader consumer forgets to reset; let's do it here again..
  rawReader->Reset() ;

  return;
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeDigits()
{
  // makes data from Digits

  GetDigitsData(1)->Fill(fDigitsArray->GetEntriesFast()) ; 
  TIter next(fDigitsArray) ; 
  AliEMCALDigit * digit ; 
  while ( (digit = dynamic_cast<AliEMCALDigit *>(next())) ) {
    GetDigitsData(0)->Fill( digit->GetAmplitude()) ;
  }  
  
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeDigits(TTree * digitTree)
{
  // makes data from Digit Tree
  if (fDigitsArray) 
    fDigitsArray->Clear() ; 
  else
    fDigitsArray = new TClonesArray("AliEMCALDigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("EMCAL") ;
  if ( ! branch ) {
    AliWarning("EMCAL branch in Digit Tree not found") ; 
  } else {
    branch->SetAddress(&fDigitsArray) ;
    branch->GetEntry(0) ; 
    MakeDigits() ; 
  }
  
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  // makes data from RecPoints
  TBranch *emcbranch = clustersTree->GetBranch("EMCALECARP");
  if (!emcbranch) { 
    AliError("can't get the branch with the EMCAL clusters !");
    return;
  }
  
  TObjArray * emcrecpoints = new TObjArray(100) ;
  emcbranch->SetAddress(&emcrecpoints);
  emcbranch->GetEntry(0);
  
  GetRecPointsData(kRecPM)->Fill(emcrecpoints->GetEntriesFast()) ; 
  TIter next(emcrecpoints) ; 
  AliEMCALRecPoint * rp ; 
  while ( (rp = dynamic_cast<AliEMCALRecPoint *>(next())) ) {
    GetRecPointsData(kRecPE)->Fill( rp->GetEnergy()) ;
    GetRecPointsData(kRecPDigM)->Fill(rp->GetMultiplicity());
  }
  emcrecpoints->Delete();
  delete emcrecpoints;
  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::SetFittingAlgorithm(Int_t fitAlgo)              
{
  //Set fitting algorithm and initialize it if this same algorithm was not set before.
  //printf("**** Set Algorithm , number %d ****\n",fitAlgo);
  
  if(fitAlgo == fFittingAlgorithm && fRawAnalyzer) {
    //Do nothing, this same algorithm already set before.
    //printf("**** Algorithm already set before, number %d, %s ****\n",fitAlgo, fRawAnalyzer->GetName());
    return;
  }
  //Initialize the requested algorithm
  if(fitAlgo != fFittingAlgorithm || !fRawAnalyzer) {
    //printf("**** Init Algorithm , number %d ****\n",fitAlgo);
		
    fFittingAlgorithm = fitAlgo; 
    if (fRawAnalyzer) delete fRawAnalyzer;  // delete prev. analyzer if existed.
		
    if (fitAlgo == kFastFit) {
      fRawAnalyzer = new AliCaloRawAnalyzerFastFit();
    }
    else if (fitAlgo == kNeuralNet) {
      fRawAnalyzer = new AliCaloRawAnalyzerNN();
    }
    else if (fitAlgo == kLMS) {
      fRawAnalyzer = new AliCaloRawAnalyzerLMS();
    }
    else if (fitAlgo == kPeakFinder) {
      fRawAnalyzer = new AliCaloRawAnalyzerPeakFinder();
    }
    else if (fitAlgo == kCrude) {
      fRawAnalyzer = new AliCaloRawAnalyzerCrude();
    }
    else {
      AliWarning("EMCAL QA invalid fit algorithm choice") ; 
    }

  }
  return;
}

