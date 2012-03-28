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

  Updated July 2011:
  ==================

  Major changes to accomodate updates of general DQM/QA changes to have per
  trigger histograms (for a given event specie).

  1) One instance of AliTPCdataQA only. (This also solves some old wishes by
  offline team to use less memory because event the 2d arrays for this object
  is not used). This now has a new flag for only keeping DQM info event by
  event! For this reason there is no need for a special DQM reset any more
  between runs!

  2) Fill the histogram for each event. The histograms are no longer filled
  from the AliTPCdataQA but per event.

  3) Use profiles for the RAW info. By adding the profiles event by event we
  get the correct event averages WITHOUT having to normalize in the end!
  Results should therefore also be directly mergable when that feature will
  come. (none of the other histograms are merged).

  This means that from the DQM/QA point of view the TPC DQM is now fully
  standard and should ease future developments.

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
#include <TH1.h> 
#include <TProfile.h> 
#include <TProfile2D.h> 

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
#include <AliDetectorRecoParam.h>


ClassImp(AliTPCQADataMakerRec)

//____________________________________________________________________________ 
AliTPCQADataMakerRec::AliTPCQADataMakerRec() : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kTPC), 
		  "TPC Rec Quality Assurance Data Maker"),
fTPCdataQA(NULL),
fRawFirstTimeBin(1),
fRawLastTimeBin(1000)
{
  // ctor
  
  for(Int_t i = 0; i < 6; i++)
    fMapping[i] = 0;
}

//____________________________________________________________________________ 
AliTPCQADataMakerRec::AliTPCQADataMakerRec(const AliTPCQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fTPCdataQA(NULL),
  fRawFirstTimeBin(qadm.GetRawFirstTimeBin()),
  fRawLastTimeBin(qadm.GetRawLastTimeBin())
{
  //copy ctor 
  // Does not copy the calibration object, instead InitRaws have to be
  // called again
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 

  for(Int_t i = 0; i < 6; i++)
    fMapping[i] = 0;
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
  delete fTPCdataQA; 

  for(Int_t i = 0; i < 6; i++) 
    delete fMapping[i];
}
 
//____________________________________________________________________________ 
void AliTPCQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  ResetEventTrigClasses();

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
  Add2ESDsList(histESDclusters, kClusters, !expert, image);

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
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line
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
  
  // It might happen that we will be in this method a few times
  // (we create the dataQA at the first call to this method)
  if(!fTPCdataQA) {
    fTPCdataQA = new AliTPCdataQA();
    LoadMaps(); // Load Altro maps
    fTPCdataQA->SetAltroMapping(fMapping); // set Altro mapping
    fTPCdataQA->SetRangeTime(fRawFirstTimeBin, fRawLastTimeBin); // set time bin interval 
    fTPCdataQA->SetIsDQM(kTRUE);
  }

  TProfile * histRawsOccupancyVsSector = 
    new TProfile("hRawsOccupancyVsSector", "Occupancy vs sector; Sector; Occupancy",
	     72, 0, 72);
  histRawsOccupancyVsSector->SetMarkerStyle(20);
  histRawsOccupancyVsSector->SetOption("P");
  histRawsOccupancyVsSector->SetStats(kFALSE);
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib )
    Add2RawsList(histRawsOccupancyVsSector, kRawsOccupancyVsSector, expert, image, !saveCorr);
  else
    Add2RawsList(histRawsOccupancyVsSector, kRawsOccupancyVsSector, !expert, image, !saveCorr);
    
  TProfile * histRawsQVsSector = 
    new TProfile("hRawsQVsSector", "<Q> vs sector; Sector; <Q>",
	     72, 0, 72);
  Add2RawsList(histRawsQVsSector, kRawsQVsSector, expert, !image, !saveCorr);

  TProfile * histRawsQmaxVsSector = 
    new TProfile("hRawsQmaxVsSector", "<Qmax> vs sector; Sector; <Qmax>",
	     72, 0, 72);
  histRawsQmaxVsSector->SetMarkerStyle(20);
  histRawsQmaxVsSector->SetOption("P");
  histRawsQmaxVsSector->SetStats(kFALSE);
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib )
    Add2RawsList(histRawsQmaxVsSector, kRawsQmaxVsSector, expert, image, !saveCorr);
  else
    Add2RawsList(histRawsQmaxVsSector, kRawsQmaxVsSector, !expert, image, !saveCorr);

  TProfile2D * histRawsOccupancy2dVsSector = 
    new TProfile2D("hRawsOccupancy2dVsSector", "Occupancy vs sector; Sector; Patch",
		   72, 0, 36, 6, 0, 6);
  histRawsOccupancy2dVsSector->SetOption("COLZ");
  histRawsOccupancy2dVsSector->SetStats(kFALSE);
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib )
    Add2RawsList(histRawsOccupancy2dVsSector, kRawsOccupancy2dVsSector, expert, image, !saveCorr);
  else
    Add2RawsList(histRawsOccupancy2dVsSector, kRawsOccupancy2dVsSector, !expert, image, !saveCorr);
    
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
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
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
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
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
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
    FillESDsData(kClusters,nTPCclusters);
    FillESDsData(kRatio,Float_t(nTPCclusters)/Float_t(nTPCclustersFindable));
    FillESDsData(kPt,track->Pt()); 
  }
  //
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  //
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
  if (! fTPCdataQA ) {

    AliError("No TPC data QA (no call to InitRaws?)!!!!") ; 
  } else {  

    if(fTPCdataQA->GetIsDQM() == kFALSE)
      AliError("Data QA has to be initialized as DQM!!!!") ; 

    // Fill profile data
    fTPCdataQA->ResetProfiles();
    
    if(fTPCdataQA->ProcessEvent(rawReader)) { // means that TPC data was processed  

      fTPCdataQA->FillOccupancyProfile();
      
      // Fill histograms    
      TObjArray *arrRW = GetMatchingRawsData(kRawsOccupancyVsSector); // all kRawsOccupancyVsSector clones matching to triggers
      for (int ih=arrRW->GetEntriesFast();ih--;) {
	TProfile* hRawsOccupancyVsSector = dynamic_cast<TProfile*>(arrRW->At(ih));
	if (hRawsOccupancyVsSector) hRawsOccupancyVsSector->Add(fTPCdataQA->GetHistOccVsSector());
      }
      arrRW = GetMatchingRawsData(kRawsOccupancy2dVsSector);
      for (int ih=arrRW->GetEntriesFast();ih--;) {
	TProfile2D* hRawsOccupancy2dVsSector = dynamic_cast<TProfile2D*>(arrRW->At(ih));
	if (hRawsOccupancy2dVsSector) hRawsOccupancy2dVsSector->Add(fTPCdataQA->GetHistOcc2dVsSector());
      }
      arrRW = GetMatchingRawsData(kRawsQVsSector);
      for (int ih=arrRW->GetEntriesFast();ih--;) {
	TProfile* hRawsQVsSector = dynamic_cast<TProfile*>(arrRW->At(ih));
	if (hRawsQVsSector) hRawsQVsSector->Add(fTPCdataQA->GetHistQVsSector());
      }
      arrRW = GetMatchingRawsData(kRawsQmaxVsSector);
      for (int ih=arrRW->GetEntriesFast();ih--;) {
	TProfile* hRawsQmaxVsSector = dynamic_cast<TProfile*>(arrRW->At(ih));
	if (hRawsQmaxVsSector) hRawsQmaxVsSector->Add(fTPCdataQA->GetHistQmaxVsSector());
      }
      //
      IncEvCountCycleRaws();
      IncEvCountTotalRaws();
      //
    }
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
        
        FillDigitsData(kDigitsADC,dig);
      } while (digArray->Next());    
  }
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}

//____________________________________________________________________________
void AliTPCQADataMakerRec::MakeRecPoints(TTree* recTree)
{
  
  AliTPCClustersRow* clrow = 0x0;
  TBranch* branch = recTree->GetBranch("Segment");  
  branch->SetAddress(&clrow);
  TClonesArray * clarray = 0x0;

  const Int_t nEntries = Int_t(recTree->GetEntries());
  for (Int_t i = 0; i < nEntries; i++) {
    
    branch->GetEntry(i);

    clarray = clrow->GetArray();

    if (!clarray) continue;

    const Int_t nClusters = clarray->GetEntriesFast();
    for (Int_t icl=0; icl < nClusters; icl++){
      
      AliTPCclusterMI* cluster = 
	(AliTPCclusterMI*)clarray->At(icl);
      
      Float_t Qmax = cluster->GetMax();
      Float_t Q    = cluster->GetQ();
      Int_t   row  = cluster->GetRow();

      if(cluster->GetDetector()<36) { // IROC (short pads)

	FillRecPointsData(kQmaxShort,Qmax);
	FillRecPointsData(kQShort,Q);
      } else { // OROC (medium and long pads)
	row += 63;
	if(cluster->GetRow()<64) { // medium pads

	  FillRecPointsData(kQmaxMedium,Qmax);
	  FillRecPointsData(kQMedium,Q);
	} else { // long pads

	  FillRecPointsData(kQmaxLong,Qmax);
	  FillRecPointsData(kQLong,Q);
	}
      }
      
      FillRecPointsData(kRow,row);
    } // end loop over clusters
  } // end loop over tree
  //
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //
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

