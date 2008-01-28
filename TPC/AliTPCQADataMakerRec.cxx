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

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQAChecker.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliTPCCalPad.h"
#include "AliTPCClustersRow.h"
#include "AliTPCclusterMI.h"

ClassImp(AliTPCQADataMakerRec)

//____________________________________________________________________________ 
AliTPCQADataMakerRec::AliTPCQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kTPC), 
		    "TPC Rec Quality Assurance Data Maker"),
  fTPCdataQA(0),
  fHistESDclusters(0),fHistESDratio(0), fHistESDpt(0),
  fHistRawsOccupancy(0),
  fHistRecPointsQmaxShort(0), fHistRecPointsQmaxMedium(0), 
  fHistRecPointsQmaxLong(0), fHistRecPointsQShort(0), 
  fHistRecPointsQMedium(0), fHistRecPointsQLong(0),
  fHistRecPointsRow(0)
{
  // ctor
}

//____________________________________________________________________________ 
AliTPCQADataMakerRec::AliTPCQADataMakerRec(const AliTPCQADataMakerRec& qadm) :
  AliQADataMakerRec()
{
  //copy ctor 
  // Does not copy the calibration object, instead InitRaws have to be
  // called again
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 

  //
  // Associate class histogram objects to the copies in the list
  // Could also be done with the indexes
  //
  fHistESDclusters   = (TH1F*)fESDsQAList->FindObject("hESDclusters");
  fHistESDratio	     = (TH1F*)fESDsQAList->FindObject("hESDratio");
  fHistESDpt         = (TH1F*)fESDsQAList->FindObject("hESDpt");

  fHistRawsOccupancy = (TH1F*)fRawsQAList->FindObject("hRawsOccupancy");

  fHistRecPointsQmaxShort  = 
    (TH1F*)fRecPointsQAList->FindObject("hRecPointsQmaxShort");
  fHistRecPointsQmaxMedium = 
    (TH1F*)fRecPointsQAList->FindObject("hRecPointsQmaxMedium");
  fHistRecPointsQmaxLong   = 
    (TH1F*)fRecPointsQAList->FindObject("hRecPointsQmaxLong"); 
  fHistRecPointsQShort     = 
    (TH1F*)fRecPointsQAList->FindObject("hRecPointsQShort");
  fHistRecPointsQMedium    = 
    (TH1F*)fRecPointsQAList->FindObject("hRecPointsQMedium");
  fHistRecPointsQLong      = 
    (TH1F*)fRecPointsQAList->FindObject("hRecPointsQLong");
  fHistRecPointsRow      = 
    (TH1F*)fRecPointsQAList->FindObject("hRecPointsRow");
}

//__________________________________________________________________
AliTPCQADataMakerRec& AliTPCQADataMakerRec::operator = (const AliTPCQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliTPCQADataMakerRec();
  new(this) AliTPCQADataMakerRec(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliTPCQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX task, TObjArray * list)
{
  //Detector specific actions at end of cycle

  if(fTPCdataQA) { // do the final step of the QA for Raw data

    fTPCdataQA->Analyse();
    
    // get the histograms and add them to the output
    fHistRawsOccupancy = fTPCdataQA->GetNoThreshold()->MakeHisto1D(0, 1, -1);
    Add2RawsList(fHistRawsOccupancy, 0);
  }

  AliQAChecker::Instance()->Run(AliQA::kTPC, task, list) ;  
}

//____________________________________________________________________________ 
void AliTPCQADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir  
  fHistESDclusters = 
    new TH1F("hESDclusters", "N TPC clusters per track; N clusters; Counts",
	     160, 0, 160);
  fHistESDclusters->Sumw2();
  Add2ESDsList(fHistESDclusters, 0);

  fHistESDratio = 
    new TH1F("hESDratio", "Ratio: TPC clusters / findable; Ratio: cluster/findable; Counts",
	     100, 0, 1);
  fHistESDratio->Sumw2();
  Add2ESDsList(fHistESDratio, 1);
  
  fHistESDpt = 
    new TH1F("hESDpt", "P_{T} distribution; p_{T} [GeV/c]; Counts",
	     50, 0, 5);
  fHistESDpt->Sumw2();
  Add2ESDsList(fHistESDpt, 2);
}

//____________________________________________________________________________ 
void AliTPCQADataMakerRec::InitRaws()
{
  fTPCdataQA = new AliTPCdataQA();
  fTPCdataQA->SetRangeTime(0, 999); // take all 1000 time bins
}

//____________________________________________________________________________ 
void AliTPCQADataMakerRec::InitRecPoints()
{
  fHistRecPointsQmaxShort = 
    new TH1F("hRecPointsQmaxShort", "Qmax distrbution (short pads); Qmax; Counts",
	     200, 0, 1000);
  fHistRecPointsQmaxShort->Sumw2();
  Add2RecPointsList(fHistRecPointsQmaxShort, 0);

  fHistRecPointsQmaxMedium = 
    new TH1F("hRecPointsQmaxMedium", "Qmax distrbution (medium pads); Qmax; Counts",
	     200, 0, 1000);
  fHistRecPointsQmaxMedium->Sumw2();
  Add2RecPointsList(fHistRecPointsQmaxMedium, 1);

  fHistRecPointsQmaxLong = 
    new TH1F("hRecPointsQmaxLong", "Qmax distrbution (long pads); Qmax; Counts",
	     200, 0, 1000);
  fHistRecPointsQmaxLong->Sumw2();
  Add2RecPointsList(fHistRecPointsQmaxLong, 2);

  fHistRecPointsQShort = 
    new TH1F("hRecPointsQShort", "Q distrbution (short pads); Q; Counts",
	     200, 0, 5000);
  fHistRecPointsQShort->Sumw2();
  Add2RecPointsList(fHistRecPointsQShort, 3);

  fHistRecPointsQMedium = 
    new TH1F("hRecPointsQMedium", "Q distrbution (medium pads); Q; Counts",
	     200, 0, 5000);
  fHistRecPointsQMedium->Sumw2();
  Add2RecPointsList(fHistRecPointsQMedium, 4);

  fHistRecPointsQLong = 
    new TH1F("hRecPointsQLong", "Q distrbution (long pads); Q; Counts",
	     200, 0, 5000);
  fHistRecPointsQLong->Sumw2();
  Add2RecPointsList(fHistRecPointsQLong, 5);

  fHistRecPointsRow = 
    new TH1F("hRecPointsRow", "Clusters per row; Row; Counts",
	     159, 0, 159);
  fHistRecPointsRow->Sumw2();
  Add2RecPointsList(fHistRecPointsRow, 6);
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
    
    fHistESDclusters->Fill(nTPCclusters);
    fHistESDratio->Fill(Float_t(nTPCclusters)/Float_t(nTPCclustersFindable));
    fHistESDpt->Fill(track->Pt()); 
  }
}

//____________________________________________________________________________
void AliTPCQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //
  // To make QA for the RAW data we use the TPC Calibration framework 
  // to handle the data and then in the end extract the data
  //
  fTPCdataQA->ProcessEvent(rawReader);  
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

	fHistRecPointsQmaxShort->Fill(Qmax);
	fHistRecPointsQShort->Fill(Q);
      } else { // OROC (medium and long pads)
	row += 63;
	if(cluster->GetRow()<64) { // medium pads

	  fHistRecPointsQmaxMedium->Fill(Qmax);
	  fHistRecPointsQMedium->Fill(Q);
	} else { // long pads

	  fHistRecPointsQmaxLong->Fill(Qmax);
	  fHistRecPointsQLong->Fill(Q);
	}
      }
      
      fHistRecPointsRow->Fill(row);
    } // end loop over clusters
  } // end loop over tree

  delete clrow;
}
