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

ClassImp(AliEMCALQADataMakerRec)
           
//____________________________________________________________________________ 
  AliEMCALQADataMakerRec::AliEMCALQADataMakerRec() : 
    AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kEMCAL), "EMCAL Quality Assurance Data Maker"),
    fSuperModules(4), // FIXME!!! number of SuperModules; 4 for 2009; update default to 12 for later runs..
    fFirstPedestalSample(0),
    fLastPedestalSample(15),
    fMinSignalHG(0),
    fMaxSignalHG(AliEMCALGeoParams::fgkSampleMax)
{
  // ctor
}

//____________________________________________________________________________ 
AliEMCALQADataMakerRec::AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) :
  AliQADataMakerRec(), 
  fSuperModules(qadm.GetSuperModules()), 
  fFirstPedestalSample(qadm.GetFirstPedestalSample()), 
  fLastPedestalSample(qadm.GetLastPedestalSample()),  
  fMinSignalHG(qadm.GetMinSignalHG()),
  fMaxSignalHG(qadm.GetMaxSignalHG())
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
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
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kEMCAL, task, list) ;  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F * h1 = new TH1F("hESDCaloClusterE",  "ESDs CaloCluster energy in EMCAL;Energy [MeV];Counts",    200, 0., 20.) ; 
  h1->Sumw2() ;
  Add2ESDsList(h1, kESDCaloClusE, !expert, image)  ;                                                     

  TH1I * h2 = new TH1I("hESDCaloClusterM", "ESDs CaloCluster multiplicity in EMCAL;# of Clusters;Entries", 100, 0,  100) ; 
  h2->Sumw2() ;
  Add2ESDsList(h2, kESDCaloClusM, !expert, image)  ;

  TH1F * h3 = new TH1F("hESDCaloCellA",  "ESDs CaloCell amplitude in EMCAL;Energy [MeV];Counts",    500, 0., 250.) ; 
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
  
  TH1F* h0 = new TH1F("hEMCALRpE","EMCAL RecPoint energies;Energy [MeV];Counts",200, 0.,20.); //GeV
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
  Add2RawsList(h0, kNsmodLG, !expert, image, !saveCorr) ;
  TProfile * h1 = new TProfile("hHighEmcalSupermodules", "High Gain EMC: # of towers vs SuperMod;SM Id;# of towers",  
			       fSuperModules, -0.5, fSuperModules-0.5) ; 
  Add2RawsList(h1, kNsmodHG, !expert, image, !saveCorr) ;

  // where did max sample occur? (bins are towers)
  TProfile * h2 = new TProfile("hLowEmcalRawtime", "Low Gain EMC: Time at Max vs towerId;Tower Id;Time [ticks]", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h2, kTimeLG, !expert, image, !saveCorr) ;
  TProfile * h3 = new TProfile("hHighEmcalRawtime", "High Gain EMC: Time at Max vs towerId;Tower Id;Time [ticks]", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h3, kTimeHG, !expert, image, !saveCorr) ;

  // how much above pedestal was the max sample?  (bins are towers)
  TProfile * h4 = new TProfile("hLowEmcalRawMaxMinusMin", "Low Gain EMC: Max - Min vs towerId;Tower Id;Max-Min [ADC counts]", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h4, kSigLG, !expert, image, !saveCorr) ;
  TProfile * h5 = new TProfile("hHighEmcalRawMaxMinusMin", "High Gain EMC: Max - Min vs towerId;Tower Id;Max-Min [ADC counts]",
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h5, kSigHG, !expert, image, !saveCorr) ;

  // total counter: channels per event
  TH1I * h6 = new TH1I("hLowNtot", "Low Gain EMC: Total Number of found towers;# of Towers;Counts", 200, 0, nTot) ;
  h6->Sumw2() ;
  Add2RawsList(h6, kNtotLG, !expert, image, !saveCorr) ;
  TH1I * h7 = new TH1I("hHighNtot", "High Gain EMC: Total Number of found towers;# of Towers;Counts", 200,0, nTot) ;
  h7->Sumw2() ;
  Add2RawsList(h7, kNtotHG, !expert, image, !saveCorr) ;

  // pedestal (bins are towers)
  TProfile * h8 = new TProfile("hLowEmcalRawPed", "Low Gain EMC: Pedestal vs towerId;Tower Id;Pedestal [ADC counts]", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h8, kPedLG, !expert, image, !saveCorr) ;
  TProfile * h9 = new TProfile("hHighEmcalRawPed", "High Gain EMC: Pedestal vs towerId;Tower Id;Pedestal [ADC counts]",
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h9, kPedHG, !expert, image, !saveCorr) ;

  // pedestal rms (standard dev = sqrt of variance estimator for pedestal) (bins are towers)
  TProfile * h10 = new TProfile("hLowEmcalRawPedRMS", "Low Gain EMC: Pedestal RMS vs towerId;Tower Id;Width [ADC counts]", 
				nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h10, kPedRMSLG, !expert, image, !saveCorr) ;
  TProfile * h11 = new TProfile("hHighEmcalRawPedRMS", "High Gain EMC: Pedestal RMS vs towerId;Tower Id;Width [ADC counts]",
				nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h11, kPedRMSHG, !expert, image, !saveCorr) ;


  // now repeat the same for TRU and LEDMon data
  int nTot2x2 = fSuperModules * AliEMCALGeoParams::fgkEMCALTRUsPerSM * AliEMCALGeoParams::fgkEMCAL2x2PerTRU; // max number of TRU channels for all SuperModules

  // counter info: number of channels per event (bins are SM index)
  TProfile * hT0 = new TProfile("hTRUEmcalSupermodules", "TRU EMC: # of TRU channels vs SuperMod;SM Id;# of TRU channels",
				fSuperModules, -0.5, fSuperModules-0.5) ;
  Add2RawsList(hT0, kNsmodTRU, !expert, image, !saveCorr) ;

  // where did max sample occur? (bins are TRU channels)
  TProfile * hT1 = new TProfile("hTRUEmcalRawtime", "TRU EMC: Time at Max vs 2x2Id;2x2 Id;Time [ticks]", 
				nTot2x2, -0.5, nTot2x2-0.5) ;
  Add2RawsList(hT1, kTimeTRU, !expert, image, !saveCorr) ;

  // how much above pedestal was the max sample?  (bins are TRU channels)
  TProfile * hT2 = new TProfile("hTRUEmcalRawMaxMinusMin", "TRU EMC: Max - Min vs 2x2Id;2x2 Id;Max-Min [ADC counts]", 
				nTot2x2, -0.5, nTot2x2-0.5) ;
  Add2RawsList(hT2, kSigTRU, !expert, image, !saveCorr) ;

  // total counter: channels per event
  TH1I * hT3 = new TH1I("hTRUNtot", "TRU EMC: Total Number of found TRU channels;# of TRU Channels;Counts", 200, 0, nTot2x2) ;
  hT3->Sumw2() ;
  Add2RawsList(hT3, kNtotTRU, !expert, image, !saveCorr) ;

  // pedestal (bins are TRU channels)
  TProfile * hT4 = new TProfile("hTRUEmcalRawPed", "TRU EMC: Pedestal vs 2x2Id;2x2 Id;Pedestal [ADC counts]", 
				nTot2x2, -0.5, nTot2x2-0.5) ;
  Add2RawsList(hT4, kPedTRU, !expert, image, !saveCorr) ;

  // pedestal rms (standard dev = sqrt of variance estimator for pedestal) (bins are TRU channels)
  TProfile * hT5 = new TProfile("hTRUEmcalRawPedRMS", "TRU EMC: Pedestal RMS vs 2x2Id;2x2 Id;Width [ADC counts]", 
				nTot2x2, -0.5, nTot2x2-0.5) ;
  Add2RawsList(hT5, kPedRMSTRU, !expert, image, !saveCorr) ;

  // and also LED Mon..
  // LEDMon has both high and low gain channels, just as regular FEE/towers
  int nTotLEDMon = fSuperModules * AliEMCALGeoParams::fgkEMCALLEDRefs; // max number of LEDMon channels for all SuperModules

  // counter info: number of channels per event (bins are SM index)
  TProfile * hL0 = new TProfile("hLowLEDMonEmcalSupermodules", "LowLEDMon Gain EMC: # of strips vs SuperMod;SM Id;# of strips",
			       fSuperModules, -0.5, fSuperModules-0.5) ;
  Add2RawsList(hL0, kNsmodLGLEDMon, !expert, image, !saveCorr) ;
  TProfile * hL1 = new TProfile("hHighLEDMonEmcalSupermodules", "HighLEDMon Gain EMC: # of strips vs SuperMod;SM Id;# of strips",  
			       fSuperModules, -0.5, fSuperModules-0.5) ; 
  Add2RawsList(hL1, kNsmodHGLEDMon, !expert, image, !saveCorr) ;

  // where did max sample occur? (bins are strips)
  TProfile * hL2 = new TProfile("hLowLEDMonEmcalRawtime", "LowLEDMon Gain EMC: Time at Max vs stripId;Strip Id;Time [ticks]", 
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL2, kTimeLGLEDMon, !expert, image, !saveCorr) ;
  TProfile * hL3 = new TProfile("hHighLEDMonEmcalRawtime", "HighLEDMon Gain EMC: Time at Max vs stripId;Strip Id;Time [ticks]", 
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL3, kTimeHGLEDMon, !expert, image, !saveCorr) ;

  // how much above pedestal was the max sample?  (bins are strips)
  TProfile * hL4 = new TProfile("hLowLEDMonEmcalRawMaxMinusMin", "LowLEDMon Gain EMC: Max - Min vs stripId;Strip Id;Max-Min [ADC counts]", 
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL4, kSigLGLEDMon, !expert, image, !saveCorr) ;
  TProfile * hL5 = new TProfile("hHighLEDMonEmcalRawMaxMinusMin", "HighLEDMon Gain EMC: Max - Min vs stripId;Strip Id;Max-Min [ADC counts]",
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL5, kSigHGLEDMon, !expert, image, !saveCorr) ;

  // total counter: channels per event
  TH1I * hL6 = new TH1I("hLowLEDMonNtot", "LowLEDMon Gain EMC: Total Number of found strips;# of Strips;Counts", 200, 0, nTotLEDMon) ;
  hL6->Sumw2() ;
  Add2RawsList(hL6, kNtotLGLEDMon, !expert, image, !saveCorr) ;
  TH1I * hL7 = new TH1I("hHighLEDMonNtot", "HighLEDMon Gain EMC: Total Number of found strips;# of Strips;Counts", 200,0, nTotLEDMon) ;
  hL7->Sumw2() ;
  Add2RawsList(hL7, kNtotHGLEDMon, !expert, image, !saveCorr) ;

  // pedestal (bins are strips)
  TProfile * hL8 = new TProfile("hLowLEDMonEmcalRawPed", "LowLEDMon Gain EMC: Pedestal vs stripId;Strip Id;Pedestal [ADC counts]", 
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL8, kPedLGLEDMon, !expert, image, !saveCorr) ;
  TProfile * hL9 = new TProfile("hHighLEDMonEmcalRawPed", "HighLEDMon Gain EMC: Pedestal vs stripId;Strip Id;Pedestal [ADC counts]",
			       nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL9, kPedHGLEDMon, !expert, image, !saveCorr) ;

  // pedestal rms (standard dev = sqrt of variance estimator for pedestal) (bins are strips)
  TProfile * hL10 = new TProfile("hLowLEDMonEmcalRawPedRMS", "LowLEDMon Gain EMC: Pedestal RMS vs stripId;Strip Id;Width [ADC counts]", 
				nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL10, kPedRMSLGLEDMon, !expert, image, !saveCorr) ;
  TProfile * hL11 = new TProfile("hHighLEDMonEmcalRawPedRMS", "HighLEDMon Gain EMC: Pedestal RMS vs stripId;Strip Id;Width [ADC counts]",
				nTotLEDMon, -0.5, nTotLEDMon-0.5) ;
  Add2RawsList(hL11, kPedRMSHGLEDMon, !expert, image, !saveCorr) ;
  
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

  //Raw histogram filling not yet implemented
  //
  //Need to figure out how to get the info we want without having to
  //actually run Raw2Digits twice.
  //I suspect what we actually want is a raw digits method, not a true
  //emcal raw data method, but this doesn't seem to be allowed in
  //AliQADataMakerRec.h

  // For now, to avoid redoing the expensive signal fits we just
  // look at max vs min of the signal spextra, a la online usage in
  // AliCaloCalibPedestal

  rawReader->Reset() ;
  AliCaloRawStreamV3 in(rawReader,"EMCAL"); 

  // setup
  int nTowersPerSM = AliEMCALGeoParams::fgkEMCALRows * AliEMCALGeoParams::fgkEMCALCols; // number of towers in a SuperModule; 24x48
  int nRows = AliEMCALGeoParams::fgkEMCALRows; // number of rows per SuperModule
  int nStripsPerSM = AliEMCALGeoParams::fgkEMCALLEDRefs; // number of strips per SuperModule
  int n2x2PerSM = AliEMCALGeoParams::fgkEMCALTRUsPerSM * AliEMCALGeoParams::fgkEMCAL2x2PerTRU; // number of TRU 2x2's per SuperModule

  int sampleMin = 0; 
  int sampleMax = AliEMCALGeoParams::fgkSampleMax; // 0x3ff = 1023 = 10-bit range

  // for the pedestal calculation
  Bool_t selectPedestalSamples = kTRUE;

  // SM counters; decl. should be safe, assuming we don't get more than expected SuperModules..
  int nTotalSMLG[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nTotalSMHG[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nTotalSMTRU[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nTotalSMLGLEDMon[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nTotalSMHGLEDMon[AliEMCALGeoParams::fgkEMCALModules] = {0};

  // indices for the reading
  int iSM = 0;
  int sample = 0;
  int time = 0;
  // counters, on sample level
  int i = 0; // the sample number in current event.
  int maxTime = 0;
  int startBin = 0;

  // calc. quantities
  double meanPed = 0, squaredMean = 0, rmsPed = 0;

  // start loop over input stream  
  while (in.NextDDL()) {
    int iRCU = in.GetDDLNumber() % 2; // RCU0 or RCU1, within SuperModule
    while (in.NextChannel()) {

      // counters
      int max = sampleMin, min = sampleMax; // min and max sample values
      int nsamples = 0;

      // for the pedestal calculation
      int sampleSum = 0; // sum of samples
      int squaredSampleSum = 0; // sum of samples squared
      int nSum = 0; // number of samples in sum
      
      while (in.NextBunch()) {
	const UShort_t *sig = in.GetSignals();
	startBin = in.GetStartTimeBin();
	nsamples += in.GetBunchLength();
	for (i = 0; i < in.GetBunchLength(); i++) {
	  sample = sig[i];
	  time = startBin--;

	  // check if it's a min or max value
	  if (sample < min) min = sample;
	  if (sample > max) {
	    max = sample;
	    maxTime = time;
	  }

	  // should we add it for the pedestal calculation?
	  if ( (fFirstPedestalSample<=time && time<=fLastPedestalSample) || // sample time in range
	       !selectPedestalSamples ) { // or we don't restrict the sample range.. - then we'll take all 
	    sampleSum += sample;
	    squaredSampleSum += sample*sample;
	    nSum++;
	  }
	  
	} // loop over samples in bunch
      } // loop over bunches
    
      if (nsamples > 0) { // this check is needed for when we have zero-supp. on, but not sparse readout

      // calculate pedesstal estimate: mean of possibly selected samples
      if (nSum > 0) {
	meanPed = sampleSum / (1.0 * nSum);
	squaredMean = squaredSampleSum / (1.0 * nSum);
	// The variance (rms squared) is equal to the mean of the squares minus the square of the mean..
	rmsPed = sqrt(squaredMean - meanPed*meanPed); 
      }
      else {
	meanPed = 0;
	squaredMean = 0;
	rmsPed  = 0;
      }

      // it should be enough to check the SuperModule info for each DDL really, but let's keep it here for now
      iSM = in.GetModule(); //The modules are numbered starting from 0

      if (iSM>=0 && iSM<fSuperModules) { // valid module reading, can go on with filling

	if ( in.IsLowGain() || in.IsHighGain() ) { // regular towers
	  int towerId = iSM*nTowersPerSM + in.GetColumn()*nRows + in.GetRow();

	  if ( in.IsLowGain() ) { 
	    //fill the low gain histograms, and counters
	    nTotalSMLG[iSM]++; // one more channel found
	    GetRawsData(kSigLG)->Fill(towerId, max - min);
	    GetRawsData(kTimeLG)->Fill(towerId, maxTime);
	    if (nSum>0) { // only fill pedestal info in case it could be calculated
	      GetRawsData(kPedLG)->Fill(towerId, meanPed);
	      GetRawsData(kPedRMSLG)->Fill(towerId, rmsPed);
	    }
	  } // gain==0
	  else if ( in.IsHighGain() ) {       	
	    //fill the high gain ones
	    nTotalSMHG[iSM]++; // one more channel found
	    int signal = max - min;
	    // only fill the max-min signal info and maxTime, if the
	    // signal was in the selected range 
	    if ( (signal > fMinSignalHG) && (signal < fMaxSignalHG) ) { 
	      GetRawsData(kSigHG)->Fill(towerId, signal);
	      GetRawsData(kTimeHG)->Fill(towerId, maxTime);
	    } // signal
	    if (nSum>0) { // only fill pedestal info in case it could be calculated
	      GetRawsData(kPedHG)->Fill(towerId, meanPed);
	      GetRawsData(kPedRMSHG)->Fill(towerId, rmsPed);
	    }
	  }
	} // low or high gain
	// TRU
	else if ( in.IsTRUData() ) {
	  // for TRU data, the mapping class holds the TRU internal 2x2 number (0..95) in the Column var..
	  int iTRU = iRCU; //TRU0 is from RCU0, TRU1 from RCU1
	  if (iRCU>0 && in.GetBranch()>0) iTRU=2; // TRU2 is from branch B on RCU1
	  int TRU2x2Id = iSM*n2x2PerSM + iTRU*AliEMCALGeoParams::fgkEMCAL2x2PerTRU 
	    + in.GetColumn();

	  //fill the low gain histograms, and counters
	  nTotalSMTRU[iSM]++; // one more channel found
	  GetRawsData(kSigTRU)->Fill(TRU2x2Id, max - min);
	  GetRawsData(kTimeTRU)->Fill(TRU2x2Id, maxTime);
	  if (nSum>0) { // only fill pedestal info in case it could be calculated
	    GetRawsData(kPedTRU)->Fill(TRU2x2Id, meanPed);
	    GetRawsData(kPedRMSTRU)->Fill(TRU2x2Id, rmsPed);
	  }
	}
	// LED Mon
	else if ( in.IsLEDMonData() ) {
	  // for LED Mon data, the mapping class holds the gain info in the Row variable
	  // and the Strip number in the Column..
	  int gain = in.GetRow(); 
	  int stripId = iSM*nStripsPerSM + in.GetColumn();
	  
	  if ( gain == 0 ) { 
	    //fill the low gain histograms, and counters
	    nTotalSMLGLEDMon[iSM]++; // one more channel found
	    GetRawsData(kSigLGLEDMon)->Fill(stripId, max - min);
	    GetRawsData(kTimeLGLEDMon)->Fill(stripId, maxTime);
	    if (nSum>0) { // only fill pedestal info in case it could be calculated
	      GetRawsData(kPedLGLEDMon)->Fill(stripId, meanPed);
	      GetRawsData(kPedRMSLGLEDMon)->Fill(stripId, rmsPed);
	    }
	  } // gain==0
	  else if ( gain == 1 ) {       	
	    //fill the high gain ones
	    nTotalSMHGLEDMon[iSM]++; // one more channel found
	    GetRawsData(kSigHGLEDMon)->Fill(stripId, max - min);
	    GetRawsData(kTimeHGLEDMon)->Fill(stripId, maxTime);
	    if (nSum>0) { // only fill pedestal info in case it could be calculated
	      GetRawsData(kPedHGLEDMon)->Fill(stripId, meanPed);
	      GetRawsData(kPedRMSHGLEDMon)->Fill(stripId, rmsPed);
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
    nTotalLG += nTotalSMLGLEDMon[iSM]; 
    nTotalHG += nTotalSMHGLEDMon[iSM]; 
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
    GetDigitsData(0)->Fill( digit->GetAmp()) ;
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

