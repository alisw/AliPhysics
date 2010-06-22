/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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


/* $Id: $ */

/*
  Based on AliPHOSQADataMaker
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  P. Christiansen, Lund, January 2008

  Updated June 2010:
  ==================

  The "beautification" of the online DQM histograms have been moved to
  an amore macro.  

  The per event RAW histograms have been modified in AliTPCdataQA and
  the copies have therefore also been modified here.

  The AliTPCdataQA can now be configured a bit from here: time bin
  range (extended default range to 1-1000, event range at start:
  0-100000, 1000 events per bin). (At least the parameters are not
  hardcoded:-)

  Implementation:
  ===============

  For the QA of the RAW data we use the class, AliTPCdataQA, from the
  existing TPC Calibration framework (which is more advanced than the
  standard QA framework) and extract the histograms at the end. The
  Analyse method of the AliTPCdataQA class is called in the method,
  EndOfDetectorCycle, and there also: 1d histogram(s) are projected
  and added to the QA list.
*/

#include "AliTPCQADataMakerRec.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TString.h>
#include <TSystem.h>
#include <TBox.h>
#include <TLine.h>
#include <TAxis.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQAChecker.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "AliTPCClustersRow.h"
#include "AliTPCclusterMI.h"
#include "AliSimDigits.h"


ClassImp(AliTPCQADataMakerRec)

//____________________________________________________________________________ 
AliTPCQADataMakerRec::AliTPCQADataMakerRec() : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kTPC), 
		  "TPC Rec Quality Assurance Data Maker"),
fTPCdataQA(NULL),
fRawMaxEvents(100000),
fRawEventsPerBin(1000),
fRawFirstTimeBin(1),
fRawLastTimeBin(1000)
{
  // ctor
  fTPCdataQA = new AliTPCdataQA*[AliRecoParam::kNSpecies] ;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    fTPCdataQA[specie] = NULL ; 
  
  for(Int_t i = 0; i < 6; i++)
    fMapping[i] = 0;
}

//____________________________________________________________________________ 
AliTPCQADataMakerRec::AliTPCQADataMakerRec(const AliTPCQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fTPCdataQA(NULL),
  fRawMaxEvents(qadm.GetRawMaxEvents()),
  fRawEventsPerBin(qadm.GetRawEventsPerBin()),
  fRawFirstTimeBin(qadm.GetRawFirstTimeBin()),
  fRawLastTimeBin(qadm.GetRawLastTimeBin())
{
  //copy ctor 
  // Does not copy the calibration object, instead InitRaws have to be
  // called again
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 

  fTPCdataQA = new AliTPCdataQA*[AliRecoParam::kNSpecies] ;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    fTPCdataQA[specie] = NULL ; 
  
  for(Int_t i = 0; i < 6; i++)
    fMapping[i] = 0;

  //
  // Associate class histogram objects to the copies in the list
  // Could also be done with the indexes
  //

}

//__________________________________________________________________
AliTPCQADataMakerRec& AliTPCQADataMakerRec::operator = (const AliTPCQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliTPCQADataMakerRec();
  new(this) AliTPCQADataMakerRec(qadm);
  return *this;
}

//__________________________________________________________________
AliTPCQADataMakerRec::~AliTPCQADataMakerRec()
{
  // Destructor
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    if ( fTPCdataQA[specie] != NULL )
      delete fTPCdataQA[specie] ; 
  delete[] fTPCdataQA; 

  for(Int_t i = 0; i < 6; i++) 
    delete fMapping[i];
}
 
//____________________________________________________________________________ 
void AliTPCQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    if(fTPCdataQA[specie] != NULL) { // do the final step of the QA for Raw data

      fTPCdataQA[specie]->Analyse(); // 31/1-08 Analyse is now protected against
                           //         RAW data files with no TPC data

      SetEventSpecie(AliRecoParam::ConvertIndex(specie)) ; 
      TH1F * histRawsOccupancy                 = (TH1F*)GetRawsData(kRawsOccupancy) ;
      TH1F * histRawsOccupancyVsSector         = (TH1F*)GetRawsData(kRawsOccupancyVsSector) ;
      TH1F * histRawsNClustersPerEventVsSector = (TH1F*)GetRawsData(kRawsNClustersPerEventVsSector) ;
      TH1F * histRawsQVsSector                 = (TH1F*)GetRawsData(kRawsQVsSector) ;
      TH1F * histRawsQmaxVsSector              = (TH1F*)GetRawsData(kRawsQmaxVsSector) ;
      TH1F * histRawsOccupancyVsEvent          = (TH1F*)GetRawsData(kRawsOccupancyVsEvent) ;
      TH1F * histRawsNclustersVsEvent          = (TH1F*)GetRawsData(kRawsNclustersVsEvent) ;
      if ( !histRawsOccupancy ||
	   !histRawsOccupancyVsSector ||
	   !histRawsNClustersPerEventVsSector ||
	   !histRawsQVsSector ||
	   !histRawsQmaxVsSector ||
	   !histRawsOccupancyVsEvent ||
	   !histRawsNclustersVsEvent ) {
        AliError("Something very wrong here, corrupted memory ?????. Please check\n") ; 
        continue ; 
      }
        
      //Add2RawsList(fTPCdataQA, 0);
      // get the histograms and add them to the output
      // 31/8-08 Histogram is only added if the Calibration class 
      //         receives TPC data 
      const Int_t eventCounter = fTPCdataQA[specie]->GetEventCounter();
      if(eventCounter>0) { // some TPC data has been processed

        // Reset histograms and refill them 
        histRawsOccupancy->Reset();
        histRawsOccupancyVsSector->Reset();
        histRawsNClustersPerEventVsSector->Reset();
        histRawsQVsSector->Reset();
        histRawsQmaxVsSector->Reset();
      
        TH1F* hNormOcc = new TH1F("hNormOcc", 0, 72, 0, 72);
        hNormOcc->Sumw2();
        TH1F* hNormNclusters = new TH1F("hNormNclusters", 0, 72, 0, 72);
        hNormNclusters->Sumw2();

        for (Int_t iSec = 0; iSec < 72; iSec++) {
	
          AliTPCCalROC* occupancyROC = 
          fTPCdataQA[specie]->GetNoThreshold()->GetCalROC(iSec); 
          AliTPCCalROC* nclusterROC = 
          fTPCdataQA[specie]->GetNLocalMaxima()->GetCalROC(iSec); 
          AliTPCCalROC* qROC = 
          fTPCdataQA[specie]->GetMeanCharge()->GetCalROC(iSec); 
          AliTPCCalROC* qmaxROC = 
          fTPCdataQA[specie]->GetMaxCharge()->GetCalROC(iSec); 

          const Int_t nRows = occupancyROC->GetNrows(); 
          for (Int_t iRow = 0; iRow < nRows; iRow++) {

            const Int_t nPads = occupancyROC->GetNPads(iRow); 
            for (Int_t iPad = 0; iPad < nPads; iPad++) {
	      
              histRawsOccupancy->Fill(occupancyROC->GetValue(iRow, iPad));
              hNormOcc->Fill(iSec);
              histRawsOccupancyVsSector
		->Fill(iSec, occupancyROC->GetValue(iRow, iPad));
	      
              const Int_t nClusters = TMath::Nint(nclusterROC->GetValue(iRow, iPad));
	      
              if(nClusters>0) {
		
		hNormNclusters->Fill(iSec,nClusters);
                histRawsNClustersPerEventVsSector->Fill(iSec, nClusters);
                histRawsQVsSector->Fill(iSec, 
					nClusters*qROC->GetValue(iRow, iPad));
                histRawsQmaxVsSector->Fill(iSec, 
					   nClusters*qmaxROC->GetValue(iRow, iPad));
              }
            }
          }
        } // end loop over sectors
      
	// update event histograms - copy info from TPDdataQA histos
	const TH1F* hQAOccVsEvent = fTPCdataQA[specie]->GetHistOccupancyVsEvent();
	const TH1F* hQANclVsEvent = fTPCdataQA[specie]->GetHistNclustersVsEvent();
	
	// In case the histogram limits have changed we have to update
	// them here
 	if(histRawsOccupancy->GetXaxis()->GetXmax()!=
 	   hQAOccVsEvent->GetXaxis()->GetXmax()) {
	  
	  histRawsOccupancyVsEvent->GetXaxis()->Set(histRawsOccupancyVsEvent->GetXaxis()->GetNbins(), hQAOccVsEvent->GetXaxis()->GetXmin(), hQAOccVsEvent->GetXaxis()->GetXmax()); 
	  
	  histRawsNclustersVsEvent->GetXaxis()->Set(histRawsOccupancyVsEvent->GetXaxis()->GetNbins(), hQANclVsEvent->GetXaxis()->GetXmin(), hQANclVsEvent->GetXaxis()->GetXmax()); 
 	}
	
	// reset the number of entries
	histRawsOccupancyVsEvent->SetEntries(0);
	histRawsNclustersVsEvent->SetEntries(0);

	// the two event histograms should have the same number of bins
	const Int_t nBins = hQAOccVsEvent->GetXaxis()->GetNbins();
	for(Int_t bin = 1; bin <= nBins; bin++) {
	  
	  histRawsOccupancyVsEvent->SetBinContent(bin, hQAOccVsEvent->GetBinContent(bin));
	  histRawsNclustersVsEvent->SetBinContent(bin, hQANclVsEvent->GetBinContent(bin));
	}
	
        // Normalize histograms
        histRawsOccupancyVsSector->Divide(hNormOcc);
        histRawsNClustersPerEventVsSector->Scale(1.0/Float_t(eventCounter));
        histRawsQVsSector->Divide(hNormNclusters);
        histRawsQmaxVsSector->Divide(hNormNclusters);
        delete hNormOcc;
        delete hNormNclusters;
	
      }
    }
  }
  AliQAChecker::Instance()->Run(AliQAv1::kTPC, task, list) ;  
}


//____________________________________________________________________________ 
void AliTPCQADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F * histESDclusters = 
    new TH1F("hESDclusters", "N TPC clusters per track; N clusters; Counts",
	     160, 0, 160);
  histESDclusters->Sumw2();
  Add2ESDsList(histESDclusters, KClusters, !expert, image);

  TH1F * histESDratio = 
    new TH1F("hESDratio", "Ratio: TPC clusters / findable; Ratio: cluster/findable; Counts",
	     100, 0, 1);
  histESDratio->Sumw2();
  Add2ESDsList(histESDratio, kRatio, !expert, image);
  
  TH1F * histESDpt = 
    new TH1F("hESDpt", "P_{T} distribution; p_{T} [GeV/c]; Counts",
	     50, 0, 5);
  histESDpt->Sumw2();
  Add2ESDsList(histESDpt, kPt, !expert, image);
}

//____________________________________________________________________________ 
void AliTPCQADataMakerRec::InitRaws()
{
  //
  // Adding the raw 
  //  

  // Modified: 7/7 - 2008
  // Laurent Aphecetche pointed out that the mapping was read from file
  // for each event, so now we read in the map here and set if for 
  // the raw data qa
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    
    // It might happen that we will be in this method a few times because
    // we create all dataQAs at the first call to this method
    if(fTPCdataQA[specie]!=0) // data QA already created
      continue;

    fTPCdataQA[specie] = 
      new AliTPCdataQA(AliRecoParam::ConvertIndex(specie));
    LoadMaps(); // Load Altro maps
    fTPCdataQA[specie]->SetAltroMapping(fMapping); // set Altro mapping
    fTPCdataQA[specie]->SetRangeTime(fRawFirstTimeBin, fRawLastTimeBin); // set time bin interval 
    fTPCdataQA[specie]->SetMaxEvents(fRawMaxEvents);
    fTPCdataQA[specie]->SetEventsPerBin(fRawEventsPerBin);
//    Add2RawsList(fTPCdataQA, kTPCdataQ, !expert, image, !saveCorrA); // This is used by the AMORE monitoring <------- THIS WILL FAIL (YS)
  }

  TH1F * histRawsOccupancy = 
    new TH1F("hRawsOccupancy", "Occupancy (all pads); Occupancy; Counts",
	     100, 0, 1);
  histRawsOccupancy->Sumw2();
  Add2RawsList(histRawsOccupancy, kRawsOccupancy, expert, !image, !saveCorr);
  
  TH1F * histRawsOccupancyVsSector = 
    new TH1F("hRawsOccupancyVsSector", "Occupancy vs sector; Sector; Occupancy",
	     72, 0, 72);
  histRawsOccupancyVsSector->Sumw2();
  histRawsOccupancyVsSector->SetMarkerStyle(20);
  histRawsOccupancyVsSector->SetOption("P");
  histRawsOccupancyVsSector->SetStats(kFALSE);
  Add2RawsList(histRawsOccupancyVsSector, kRawsOccupancyVsSector, !expert, image, !saveCorr);

  TH1F * histRawsNClustersPerEventVsSector = 
    new TH1F("hRawsNClustersPerEventVsSector", "Nclusters per event vs sector; Sector; Nclusters per event",
	     72, 0, 72);
  histRawsNClustersPerEventVsSector->Sumw2();
  Add2RawsList(histRawsNClustersPerEventVsSector, kRawsNClustersPerEventVsSector, expert, !image, !saveCorr);
  
  TH1F * histRawsQVsSector = 
    new TH1F("hRawsQVsSector", "<Q> vs sector; Sector; <Q>",
	     72, 0, 72);
  histRawsQVsSector->Sumw2();
  Add2RawsList(histRawsQVsSector, kRawsQVsSector, expert, !image, !saveCorr);

  TH1F * histRawsQmaxVsSector = 
    new TH1F("hRawsQmaxVsSector", "<Qmax> vs sector; Sector; <Qmax>",
	     72, 0, 72);
  histRawsQmaxVsSector->Sumw2();
  histRawsQmaxVsSector->SetMarkerStyle(20);
  histRawsQmaxVsSector->SetOption("P");
  histRawsQmaxVsSector->SetStats(kFALSE);
  Add2RawsList(histRawsQmaxVsSector, kRawsQmaxVsSector, !expert, image, !saveCorr);

  // Get histogram information from data QA to build copy
  const TH1F* hOccHelp = fTPCdataQA[0]->GetHistOccupancyVsEvent();
  TH1F * histRawsOccupancyVsEvent = 
    CreateEventsHistCopy(hOccHelp, "hRawsOccupancyVsEvent");
  Add2RawsList(histRawsOccupancyVsEvent, kRawsOccupancyVsEvent, !expert, image, !saveCorr);
  
  // Get histogram information from data QA to build copy
  const TH1F* hNclHelp = fTPCdataQA[0]->GetHistNclustersVsEvent();
  TH1F * histRawsNclustersVsEvent = 
    CreateEventsHistCopy(hNclHelp, "hRawsNclustersVsEvent");
  Add2RawsList(histRawsNclustersVsEvent, kRawsNclustersVsEvent, !expert, image, !saveCorr);
}

//____________________________________________________________________________ 
void AliTPCQADataMakerRec::InitDigits()
{
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  TH1F * histDigitsADC = 
    new TH1F("hDigitsADC", "Digit ADC distribution; ADC; Counts",
             1000, 0, 1000);
  histDigitsADC->Sumw2();
  Add2DigitsList(histDigitsADC, kDigitsADC, !expert, image);
}

//____________________________________________________________________________ 
void AliTPCQADataMakerRec::InitRecPoints()
{
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F * histRecPointsQmaxShort = 
    new TH1F("hRecPointsQmaxShort", "Qmax distrbution (short pads); Qmax; Counts",
	     100, 0, 300);
  histRecPointsQmaxShort->Sumw2();
  Add2RecPointsList(histRecPointsQmaxShort, kQmaxShort, !expert, image);

  TH1F * histRecPointsQmaxMedium = 
    new TH1F("hRecPointsQmaxMedium", "Qmax distrbution (medium pads); Qmax; Counts",
	     100, 0, 300);
  histRecPointsQmaxMedium->Sumw2();
  Add2RecPointsList(histRecPointsQmaxMedium, kQmaxMedium, !expert, image);

  TH1F * histRecPointsQmaxLong = 
    new TH1F("hRecPointsQmaxLong", "Qmax distrbution (long pads); Qmax; Counts",
	     100, 0, 300);
  histRecPointsQmaxLong->Sumw2();
  Add2RecPointsList(histRecPointsQmaxLong, kQmaxLong, !expert, image);

  TH1F * histRecPointsQShort = 
    new TH1F("hRecPointsQShort", "Q distrbution (short pads); Q; Counts",
	     100, 0, 2000);
  histRecPointsQShort->Sumw2();
  Add2RecPointsList(histRecPointsQShort, kQShort, !expert, image);

  TH1F * histRecPointsQMedium = 
    new TH1F("hRecPointsQMedium", "Q distrbution (medium pads); Q; Counts",
	     100, 0, 2000);
  histRecPointsQMedium->Sumw2();
  Add2RecPointsList(histRecPointsQMedium, kQMedium, !expert, image);

  TH1F * histRecPointsQLong = 
    new TH1F("hRecPointsQLong", "Q distrbution (long pads); Q; Counts",
	     100, 0, 2000);
  histRecPointsQLong->Sumw2();
  Add2RecPointsList(histRecPointsQLong, kQLong, !expert, image);

  TH1F * histRecPointsRow = 
    new TH1F("hRecPointsRow", "Clusters per row; Row; Counts",
	     159, 0, 159);
  histRecPointsRow->Sumw2();
  Add2RecPointsList(histRecPointsRow, kRow, !expert, image);
}

//____________________________________________________________________________
void AliTPCQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
 
  const Int_t nESDTracks = esd->GetNumberOfTracks();
  Int_t nTPCtracks = 0; 
  for(Int_t i = 0; i < nESDTracks; i++) {
    
    AliESDtrack * track = esd->GetTrack(i);
    
    if ((track->GetStatus() & AliESDtrack::kTPCrefit)==0)
      continue;
    
    nTPCtracks++;
    
    Int_t nTPCclusters         = track->GetTPCNcls();
    Int_t nTPCclustersFindable = track->GetTPCNclsF();
    if ( nTPCclustersFindable<=0) continue;
    GetESDsData(KClusters)->Fill(nTPCclusters);
    GetESDsData(kRatio)->Fill(Float_t(nTPCclusters)/Float_t(nTPCclustersFindable));
    GetESDsData(kPt)->Fill(track->Pt()); 
  }
}

//____________________________________________________________________________
void AliTPCQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //
  // To make QA for the RAW data we use the TPC Calibration framework 
  // to handle the data and then in the end extract the data
  //
  
  GetRawsData(0); // dummy call to init raw data
  rawReader->Reset() ; 
  if (! fTPCdataQA[AliRecoParam::AConvert(fEventSpecie)] ) {
    AliError("Something unexpected here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!") ; 
  } else {  
    fTPCdataQA[AliRecoParam::AConvert(fEventSpecie)]->ProcessEvent(rawReader);  
  }
}

//____________________________________________________________________________
void AliTPCQADataMakerRec::MakeDigits(TTree* digitTree)
{
 
  TBranch* branch = digitTree->GetBranch("Segment");
  AliSimDigits* digArray = 0;
  branch->SetAddress(&digArray);
  
  Int_t nEntries = Int_t(digitTree->GetEntries());
  
  for (Int_t n = 0; n < nEntries; n++) {
    
    digitTree->GetEvent(n);
    
    if (digArray->First())
      do {
        Float_t dig = digArray->CurrentDigit();
        
        GetDigitsData(kDigitsADC)->Fill(dig);
      } while (digArray->Next());    
  }
}

//____________________________________________________________________________
void AliTPCQADataMakerRec::MakeRecPoints(TTree* recTree)
{

  AliTPCClustersRow clrow;
  clrow.SetClass("AliTPCclusterMI");
  clrow.SetArray(0);
  clrow.GetArray()->ExpandCreateFast(10000);
  AliTPCClustersRow * pclrow = &clrow;
  TBranch* branch = recTree->GetBranch("Segment");
  
  branch->SetAddress(&pclrow);

  const Int_t nEntries = Int_t(recTree->GetEntries());
  for (Int_t i = 0; i < nEntries; i++) {
    
    branch->GetEntry(i);
    
    const Int_t nClusters = clrow.GetArray()->GetEntriesFast();
    for (Int_t icl=0; icl < nClusters; icl++){
      
      AliTPCclusterMI* cluster = 
	(AliTPCclusterMI*)clrow.GetArray()->At(icl);
      
      Float_t Qmax = cluster->GetMax();
      Float_t Q    = cluster->GetQ();
      Int_t   row  = cluster->GetRow();

      if(cluster->GetDetector()<36) { // IROC (short pads)

	GetRecPointsData(kQmaxShort)->Fill(Qmax);
	GetRecPointsData(kQShort)->Fill(Q);
      } else { // OROC (medium and long pads)
	row += 63;
	if(cluster->GetRow()<64) { // medium pads

	  GetRecPointsData(kQmaxMedium)->Fill(Qmax);
	  GetRecPointsData(kQMedium)->Fill(Q);
	} else { // long pads

	  GetRecPointsData(kQmaxLong)->Fill(Qmax);
	  GetRecPointsData(kQLong)->Fill(Q);
	}
      }
      
      GetRecPointsData(kRow)->Fill(row);
    } // end loop over clusters
  } // end loop over tree

}

//____________________________________________________________________________
void AliTPCQADataMakerRec::LoadMaps()
{
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/TPC/mapping/Patch";

  for(Int_t i = 0; i < 6; i++) {

    if(fMapping[i]!=0) // mapping already loaded
      continue;
    TString path2 = path;
    path2 += i;
    path2 += ".data";
    fMapping[i] = new AliTPCAltroMapping(path2.Data());
  }
}

//____________________________________________________________________________
void AliTPCQADataMakerRec::ResetDetector(AliQAv1::TASKINDEX_t task)
{
  // Overwrites general method for RAW data.
  // The AliTPCdataQA elements that does the internal processing are
  // in the case they have processed data deleted and new are created

  if ( task != AliQAv1::kRAWS )
    AliQADataMakerRec::ResetDetector(task);
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    
    if ( fTPCdataQA[specie] != NULL) { // exist
      
      if(fTPCdataQA[specie]->GetEventCounter()>0) { // has processed data
	
	// old configuration
	Int_t  firstTime    = fTPCdataQA[specie]->GetFirstTimeBin();
	Int_t  lastTime     = fTPCdataQA[specie]->GetLastTimeBin();
	Int_t  minADC       = fTPCdataQA[specie]->GetAdcMin();
	Int_t  maxADC       = fTPCdataQA[specie]->GetAdcMax();

	//delete old
	delete fTPCdataQA[specie]; 

	// create new
	fTPCdataQA[specie] = new AliTPCdataQA(AliRecoParam::ConvertIndex(specie));
	// configure new
	LoadMaps(); // Load Altro maps
	fTPCdataQA[specie]->SetAltroMapping(fMapping);
	fTPCdataQA[specie]->SetRangeTime(firstTime, lastTime);
	fTPCdataQA[specie]->SetRangeAdc(minADC, maxADC);
	// Here we want to restore the default configuration because
	// the max events and events are adjusted for the last run
	fTPCdataQA[specie]->SetMaxEvents(fRawMaxEvents);
	fTPCdataQA[specie]->SetEventsPerBin(fRawEventsPerBin);
      }
    }
  }
}

//____________________________________________________________________________
TH1F* AliTPCQADataMakerRec::CreateEventsHistCopy(const TH1F* hist, 
						 const Char_t* copyName)
{
  // This method is used to create a copy of the event histograms
  
  TH1F* histCopy = new TH1F(copyName, hist->GetTitle(),
			    hist->GetXaxis()->GetNbins(),
			    hist->GetXaxis()->GetXmin(), 
			    hist->GetXaxis()->GetXmax());
  histCopy->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  histCopy->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
  histCopy->SetMarkerStyle(20);
  histCopy->SetOption("P");
  histCopy->SetStats(kFALSE);

  return histCopy;
}
