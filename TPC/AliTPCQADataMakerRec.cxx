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
*/

/*
  Implementation:

  We have chosen to have the histograms as non-persistent meber to
  allow better debugging. In the copy constructor we then have to
  assign the pointers to the existing histograms in the copied
  list. This have been implemented but not tested.

  For the QA of the RAW data we use the class, AliTPCdataQA, from the
  existing TPC Calibration framework (which is more advanced than the
  standard QA framework) and extract the histograms at the end. This
  has been tested with zero-suppressed data. The Analyse method of the
  AliTPCdataQA class is called in the method, EndOfDetectorCycle, and
  there also: 1d histogram(s) are projected and added to the QA list.
*/

/*
  TODO:
  Sumw2 for RAW histogram(s)?
  RecPoints and ESD could have many more histograms
*/

#include "AliTPCQADataMakerRec.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TString.h>
#include <TSystem.h>

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

ClassImp(AliTPCQADataMakerRec)

//____________________________________________________________________________ 
AliTPCQADataMakerRec::AliTPCQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kTPC), 
		    "TPC Rec Quality Assurance Data Maker"),
  fTPCdataQA(NULL)
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
  fTPCdataQA(NULL)
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

      SetEventSpecie(specie) ; 
      TH1F * histRawsOccupancy                 = (TH1F*)GetRawsData(kOccupancy) ;
      TH1F * histRawsOccupancyVsSector         = (TH1F*)GetRawsData(kOccupancyVsSector) ;
      TH1F * histRawsNClustersPerEventVsSector = (TH1F*)GetRawsData(kNClustersPerEventVsSector) ;
      TH1F * histRawsQVsSector                 = (TH1F*)GetRawsData(kQVsSector) ;
      TH1F * histRawsQmaxVsSector              = (TH1F*)GetRawsData(kQmaxVsSector) ;
      if ( !histRawsOccupancy ||
          !histRawsOccupancyVsSector ||
          !histRawsNClustersPerEventVsSector ||
          !histRawsQVsSector ||
          !histRawsQmaxVsSector) {
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
      
        TH1F* hNorm72 = new TH1F("hNorm72", "histogram to normalize 72 sectors",
                                 72, 0, 72);
        hNorm72->Sumw2();
        TH1F* hNorm108 = new TH1F("hNorm108", "histogram to normalize 108 sectors (medium and long pads are split up)",
                                  108, 0, 108);
        hNorm108->Sumw2();

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

            Int_t helpSector = iSec;
            if(iRow>=64)
              helpSector += 36; // OROC (long pads)

            const Int_t nPads = occupancyROC->GetNPads(iRow); 
            for (Int_t iPad = 0; iPad < nPads; iPad++) {
	
              histRawsOccupancy->Fill(occupancyROC->GetValue(iRow, iPad));
              hNorm72->Fill(iSec);
              histRawsOccupancyVsSector
              ->Fill(iSec, occupancyROC->GetValue(iRow, iPad));

              const Int_t nClusters = TMath::Nint(nclusterROC->GetValue(iRow, iPad));
	    
              if(nClusters>0) {
	      
                histRawsNClustersPerEventVsSector->Fill(iSec, nClusters);
                hNorm108->Fill(helpSector, nClusters);
                histRawsQVsSector->Fill(helpSector, 
                                       nClusters*qROC->GetValue(iRow, iPad));
                histRawsQmaxVsSector->Fill(helpSector, 
                                          nClusters*qmaxROC->GetValue(iRow, iPad));
              }
            }
          }
        } // end loop over sectors
      
        // Normalize histograms
        histRawsOccupancyVsSector->Divide(hNorm72);
        histRawsNClustersPerEventVsSector->Scale(1.0/Float_t(eventCounter));
        histRawsQVsSector->Divide(hNorm108);
        histRawsQmaxVsSector->Divide(hNorm108);
        delete hNorm72;
        delete hNorm108;
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
    fTPCdataQA[specie] = new AliTPCdataQA(AliRecoParam::Convert(specie));
    LoadMaps(); // Load Altro maps
    fTPCdataQA[specie]->SetAltroMapping(fMapping); // set Altro mapping
    fTPCdataQA[specie]->SetRangeTime(100, 920); // set time bin interval 
//    Add2RawsList(fTPCdataQA, kTPCdataQ, !expert, image, !saveCorrA); // This is used by the AMORE monitoring <------- THIS WILL FAIL (YS)
  }

  TH1F * histRawsOccupancy = 
    new TH1F("hRawsOccupancy", "Occupancy (all pads); Occupancy; Counts",
	     100, 0, 1);
  histRawsOccupancy->Sumw2();
  Add2RawsList(histRawsOccupancy, kOccupancy, !expert, image, !saveCorr);
  
  TH1F * histRawsOccupancyVsSector = 
    new TH1F("hRawsOccupancyVsSector", "Occupancy vs sector; Sector; Occupancy",
	     72, 0, 72);
  histRawsOccupancyVsSector->Sumw2();
  Add2RawsList(histRawsOccupancyVsSector, kOccupancyVsSector, !expert, image, !saveCorr);

  TH1F * histRawsNClustersPerEventVsSector = 
    new TH1F("hRawsNClustersPerEventVsSector", "Nclusters per event vs sector; Sector; Nclusters per event",
	     72, 0, 72);
  histRawsNClustersPerEventVsSector->Sumw2();
  Add2RawsList(histRawsNClustersPerEventVsSector, kNClustersPerEventVsSector, !expert, image, !saveCorr);
  
  TH1F * histRawsQVsSector = 
    new TH1F("hRawsQVsSector", "<Q> vs sector (OROC med: 36-71, long: 72-107); Sector; <Q>",
	     108, 0, 108);
  histRawsQVsSector->Sumw2();
  Add2RawsList(histRawsQVsSector, kQVsSector, !expert, image, !saveCorr);

  TH1F * histRawsQmaxVsSector = 
    new TH1F("hRawsQmaxVsSector", "<Qmax> vs sector (OROC med: 36-71, long: 72-107); Sector; <Qmax>",
	     108, 0, 108);
  histRawsQmaxVsSector->Sumw2();
  Add2RawsList(histRawsQmaxVsSector, kQmaxVsSector, !expert, image, !saveCorr);
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
	rawReader->Reset() ; 
  fTPCdataQA[AliRecoParam::AConvert(fEventSpecie)]->ProcessEvent(rawReader);  
 }

//____________________________________________________________________________
void AliTPCQADataMakerRec::MakeRecPoints(TTree* recTree)
{
  AliTPCClustersRow *clrow = new AliTPCClustersRow();
  clrow->SetClass("AliTPCclusterMI");
  clrow->SetArray(0);
  clrow->GetArray()->ExpandCreateFast(10000);

  TBranch* branch = recTree->GetBranch("Segment");
  branch->SetAddress(&clrow);

  const Int_t nEntries = Int_t(recTree->GetEntries());
  for (Int_t i = 0; i < nEntries; i++) {
    
    branch->GetEntry(i);
    
    const Int_t nClusters = clrow->GetArray()->GetEntriesFast();
    for (Int_t icl=0; icl < nClusters; icl++){
      
      AliTPCclusterMI* cluster = 
	(AliTPCclusterMI*)clrow->GetArray()->At(icl);
      
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

  delete clrow;
}

//____________________________________________________________________________
void AliTPCQADataMakerRec::LoadMaps()
{
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/TPC/mapping/Patch";

  for(Int_t i = 0; i < 6; i++) {
    TString path2 = path;
    path2 += i;
    path2 += ".data";
    fMapping[i] = new AliTPCAltroMapping(path2.Data());
  }
}

