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
#include "AliCaloRawStream.h"

ClassImp(AliEMCALQADataMakerRec)
           
//____________________________________________________________________________ 
  AliEMCALQADataMakerRec::AliEMCALQADataMakerRec() : 
    AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kEMCAL), "EMCAL Quality Assurance Data Maker"),
    fSuperModules(4) // FIXME!!! number of SuperModules; 4 for 2009; update default to 12 for later runs..
{
  // ctor
}

//____________________________________________________________________________ 
AliEMCALQADataMakerRec::AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) :
  AliQADataMakerRec()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
  fSuperModules = qadm.GetSuperModules();
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
  
  TH1F * h1 = new TH1F("hESDCaloClusterE",  "ESDs CaloCluster energy in EMCAL",    200, 0., 20.) ; 
  h1->Sumw2() ;
  Add2ESDsList(h1, kESDCaloClusE, !expert, image)  ;                                                     

  TH1I * h2 = new TH1I("hESDCaloClusterM", "ESDs CaloCluster multiplicity in EMCAL", 100, 0,  100) ; 
  h2->Sumw2() ;
  Add2ESDsList(h2, kESDCaloClusM, !expert, image)  ;

  TH1F * h3 = new TH1F("hESDCaloCellA",  "ESDs CaloCell amplitude in EMCAL",    500, 0., 250.) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, kESDCaloCellA, !expert, image)  ;  
 
  TH1I * h4 = new TH1I("hESDCaloCellM", "ESDs CaloCell multiplicity in EMCAL", 200, 0,  1000) ; 
  h4->Sumw2() ;
  Add2ESDsList(h4, kESDCaloCellM, !expert, image) ;
	
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I * h0 = new TH1I("hEmcalDigits",    "Digits amplitude distribution in EMCAL",    500, 0, 500) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0, !expert, image) ;
  TH1I * h1 = new TH1I("hEmcalDigitsMul", "Digits multiplicity distribution in EMCAL", 200, 0, 2000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, 1, !expert, image) ;
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F* h0 = new TH1F("hEMCALRpE","EMCAL RecPoint energies",200, 0.,20.); //GeV
  h0->Sumw2();
  Add2RecPointsList(h0,kRecPE, !expert, image);

  TH1I* h1 = new TH1I("hEMCALRpM","EMCAL RecPoint multiplicities",100,0,100);
  h1->Sumw2();
  Add2RecPointsList(h1,kRecPM, !expert, image);

  TH1I* h2 = new TH1I("hEMCALRpDigM","EMCAL RecPoint Digit Multiplicities",20,0,20);
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

  int nTowersPerSM = 1152; // number of towers in a SuperModule; 24x48
  int nTot = fSuperModules * nTowersPerSM; // max number of towers in all SuperModules

  // counter info: number of channels per event (bins are SM index)
  TProfile * h0 = new TProfile("hLowEmcalSupermodules", "Low Gain EMC: # of towers vs SuperMod",
			       fSuperModules, -0.5, fSuperModules-0.5) ;
  Add2RawsList(h0, kNsmodLG, !expert, image, !saveCorr) ;
  TProfile * h1 = new TProfile("hHighEmcalSupermodules", "High Gain EMC: # of towers vs SuperMod",  
			       fSuperModules, -0.5, fSuperModules-0.5) ; 
  Add2RawsList(h1, kNsmodHG, !expert, image, !saveCorr) ;

  // where did max sample occur? (bins are towers)
  TProfile * h2 = new TProfile("hLowEmcalRawtime", "Low Gain EMC: Time at Max vs towerId", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h2, kTimeLG, !expert, image, !saveCorr) ;
  TProfile * h3 = new TProfile("hHighEmcalRawtime", "High Gain EMC: Time at Max vs towerId", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h3, kTimeHG, !expert, image, !saveCorr) ;

  // how much above pedestal was the max sample?  (bins are towers)
  TProfile * h4 = new TProfile("hLowEmcalRawMaxMinusMin", "Low Gain EMC: Max - Min vs towerId", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h4, kSigLG, !expert, image, !saveCorr) ;
  TProfile * h5 = new TProfile("hHighEmcalRawMaxMinusMin", "High Gain EMC: Max - Min vs towerId",
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h5, kSigHG, !expert, image, !saveCorr) ;

  // total counter: channels per event
  TH1I * h6 = new TH1I("hLowNtot", "Low Gain EMC: Total Number of found towers", 200, 0, nTot) ;
  h6->Sumw2() ;
  Add2RawsList(h6, kNtotLG, !expert, image, !saveCorr) ;
  TH1I * h7 = new TH1I("hHighNtot", "High Gain EMC: Total Number of found towers", 200,0, nTot) ;
  h7->Sumw2() ;
  Add2RawsList(h7, kNtotHG, !expert, image, !saveCorr) ;

  // pedestal (bins are towers)
  TProfile * h8 = new TProfile("hLowEmcalRawPed", "Low Gain EMC: Pedestal vs towerId", 
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h8, kPedLG, !expert, image, !saveCorr) ;
  TProfile * h9 = new TProfile("hHighEmcalRawPed", "High Gain EMC: Pedestal vs towerId",
			       nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h9, kPedHG, !expert, image, !saveCorr) ;

  // pedestal rms (standard dev = sqrt of variance estimator for pedestal) (bins are towers)
  TProfile * h10 = new TProfile("hLowEmcalRawPedRMS", "Low Gain EMC: Pedestal RMS vs towerId", 
				nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h10, kPedRMSLG, !expert, image, !saveCorr) ;
  TProfile * h11 = new TProfile("hHighEmcalRawPedRMS", "High Gain EMC: Pedestal RMS vs towerId",
				nTot, -0.5, nTot-0.5) ;
  Add2RawsList(h11, kPedRMSHG, !expert, image, !saveCorr) ;
  
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
  AliCaloRawStream in(rawReader,"EMCAL"); 

  // setup
  int nTowersPerSM = 1152; // number of towers in a SuperModule; 24x48
  int nRows = 24; // number of rows per SuperModule
  int sampleMin = 0; 
  int sampleMax = 0x3ff; // 1023 = 10-bit range

  // for the pedestal calculation
  Bool_t selectPedestalSamples = kTRUE;
  int firstPedestalSample = 0;
  int lastPedestalSample = 15;

  // SM counters; decl. should be safe, assuming we don't get more than 12 SuperModules..
  int nTotalSMLG[12] = {0};
  int nTotalSMHG[12] = {0};

  // indices for the reading
  int iSM = 0;
  int sample = 0;
  int gain = 0;
  int time = 0;
  // counters, on sample level
  int i = 0; // the sample number in current event.
  int max = sampleMin, min = sampleMax;//Use these for picking the pedestal
  int maxTime = 0;

  // for the pedestal calculation
  int sampleSum = 0; // sum of samples
  int squaredSampleSum = 0; // sum of samples squared
  int nSum = 0; // number of samples in sum
  // calc. quantities
  double meanPed = 0, squaredMean = 0, rmsPed = 0;
  
  while (in.Next()) { // loop over input stream
    sample = in.GetSignal(); //Get the adc signal
    time = in.GetTime();
    if (sample < min) { min = sample; }
    if (sample > max) { 
      max = sample;
      maxTime = time;
    }
    i++;

    // should we add it for the pedestal calculation?
    if ( (firstPedestalSample<=time && time<=lastPedestalSample) || // sample time in range
	 !selectPedestalSamples ) { // or we don't restrict the sample range.. - then we'll take all 
      sampleSum += sample;
      squaredSampleSum += sample*sample;
      nSum++;
    }

    if ( i >= in.GetTimeLength()) {
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

      //If we're here then we're done with this tower
      gain = -1; // init to not valid value
      if ( in.IsLowGain() ) {
	gain = 0;
      }
      else if ( in.IsHighGain() ) {
	gain = 1;
      }
      
      iSM = in.GetModule(); //The modules are numbered starting from 0

      if (iSM>=0 && iSM<fSuperModules) { // valid module reading, can go on with filling

	int towerId = iSM*nTowersPerSM + in.GetColumn()*nRows + in.GetRow();

	if (gain == 0) { 
	  //fill the low gain histograms, and counters
	  nTotalSMLG[iSM]++; // one more channel found
	  GetRawsData(kSigLG)->Fill(towerId, max - min);
	  GetRawsData(kTimeLG)->Fill(towerId, maxTime);
	  if (nSum>0) { // only fill pedestal info in case it could be calculated
	    GetRawsData(kPedLG)->Fill(towerId, meanPed);
	    GetRawsData(kPedRMSLG)->Fill(towerId, rmsPed);
	  }
	} // gain==0
	else if (gain == 1) {       	
	  //fill the high gain ones
	  nTotalSMHG[iSM]++; // one more channel found
	  GetRawsData(kSigHG)->Fill(towerId, max - min);
	  GetRawsData(kTimeHG)->Fill(towerId, maxTime);
	  if (nSum>0) { // only fill pedestal info in case it could be calculated
	    GetRawsData(kPedHG)->Fill(towerId, meanPed);
	    GetRawsData(kPedRMSHG)->Fill(towerId, rmsPed);
	  }
	}
      } // SM index OK

      // reset counters
      max = sampleMin; min = sampleMax;
      maxTime = 0;
      i = 0;
      // also pedestal calc counters
      sampleSum = 0; // sum of samples
      squaredSampleSum = 0; // sum of samples squared
      nSum = 0; // number of samples in sum
    
    }//End if, of channel
   
  }//end while, of stream

  // let's also fill the SM and event counter histograms
  int nTotalHG = 0;
  int nTotalLG = 0;
  for (iSM=0; iSM<fSuperModules; iSM++) {  
    nTotalLG += nTotalSMLG[iSM]; 
    nTotalHG += nTotalSMHG[iSM]; 
    GetRawsData(kNsmodLG)->Fill(iSM, nTotalSMLG[iSM]); 
    GetRawsData(kNsmodHG)->Fill(iSM, nTotalSMHG[iSM]); 
  }
  GetRawsData(kNtotLG)->Fill(nTotalLG);
  GetRawsData(kNtotHG)->Fill(nTotalHG);

  // just in case the next rawreader consumer forgets to reset; let's do it here again..
  rawReader->Reset() ;

  return;
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeDigits(TClonesArray * digits)
{
  // makes data from Digits
  
  GetDigitsData(1)->Fill(digits->GetEntriesFast()) ; 
  TIter next(digits) ; 
  AliEMCALDigit * digit ; 
  while ( (digit = dynamic_cast<AliEMCALDigit *>(next())) ) {
    GetDigitsData(0)->Fill( digit->GetAmp()) ;
  }  
  
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeDigits(TTree * digitTree)
{
  // makes data from Digit Tree
  TClonesArray * digits = new TClonesArray("AliEMCALDigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("EMCAL") ;
  if ( ! branch ) {
    AliWarning("EMCAL branch in Digit Tree not found") ; 
  } else {
    branch->SetAddress(&digits) ;
    branch->GetEntry(0) ; 
    MakeDigits(digits) ; 
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

