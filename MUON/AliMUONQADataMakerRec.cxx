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
#include <TH3F.h> 
#include <TLorentzVector.h>

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliMUONConstants.h"  
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliQAChecker.h"
#include "AliMpBusPatch.h"
#include "AliMUONCluster.h"  
#include "AliMUONRawStreamTracker.h"
#include "AliMUONRawStreamTrigger.h"

#include "AliMUONVClusterStore.h"
#include "AliMUONVCluster.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"

#include "AliMUONDigitMaker.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONLocalTrigger.h"

#include "AliMUONQADataMakerRec.h"

//-----------------------------------------------------------------------------
/// \class AliMUONQADataMakerRec
///
/// MUON base class for quality assurance data (histo) maker
///
/// \author C. Finck

/// \cond CLASSIMP
ClassImp(AliMUONQADataMakerRec)
/// \endcond
           
//____________________________________________________________________________ 
AliMUONQADataMakerRec::AliMUONQADataMakerRec() : 
    AliQADataMakerRec(AliQA::GetDetName(AliQA::kMUON), "MUON Quality Assurance Data Maker"),
    fDigitStore(0x0),
    fTriggerStore(0x0),
    fDigitMaker(0x0)
{
    /// ctor
  fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV1");
  fDigitMaker = new AliMUONDigitMaker(kTRUE,kFALSE);

}

//____________________________________________________________________________ 
AliMUONQADataMakerRec::AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm) :
    AliQADataMakerRec(qadm),
    fDigitStore(0x0),
    fTriggerStore(0x0),
    fDigitMaker(0x0)
{
    ///copy ctor 
    SetName((const char*)qadm.GetName()) ; 
    SetTitle((const char*)qadm.GetTitle()); 

    // Do not copy the digit store and digit maker, but create its own ones
    fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV1");
    fDigitMaker = new AliMUONDigitMaker(kTRUE,kFALSE);
}

//__________________________________________________________________
AliMUONQADataMakerRec& AliMUONQADataMakerRec::operator = (const AliMUONQADataMakerRec& qadm )
{
  /// Assignment operator

  // check assignment to self
  if (this == &qadm) return *this;

  this->~AliMUONQADataMakerRec();
  new(this) AliMUONQADataMakerRec(qadm);
  return *this;
}

//__________________________________________________________________
AliMUONQADataMakerRec::~AliMUONQADataMakerRec()
{
    /// dtor
    delete fDigitStore;
    delete fTriggerStore;
    delete fDigitMaker;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray* list)
{
    ///Detector specific actions at end of cycle
    // do the QA checking
    AliQAChecker::Instance()->Run(AliQA::kMUON, task, list) ;  
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRaws()
{
    /// create Raws histograms in Raws subdir
    TH1I* h0 = new TH1I("hRawBusPatch", "buspatch distribution",  1932, 1, 1932); 
    Add2RawsList(h0, kRawBusPatch);

    TH1I* h1 = new TH1I("hRawCharge", "Charge distribution in rawdata", 4096, 0, 4095); 
    Add2RawsList(h1, kRawCharge);
		
		for (Int_t iDDL = 0; iDDL < 20; ++iDDL) 
		{
			TH1F* h2 = new TH1F(Form("%s%d", "hRawBusPatchDDL", iDDL), Form("%s %d","RAW Buspatch distribution for DDL", iDDL), 50, 0, 49); 
			Add2RawsList(h2, kRawBuspatchDDL+iDDL);
		}

}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPoints()
{
    /// create Reconstructed Points histograms in RecPoints subdir
    TH3F *h2 = new TH3F("hTriggerDigitsBendPlane", "Trigger digits in bending plane",
			4, -0.5, 4. - 0.5,
			18, -0.5, 18. - 0.5,
			7*64, -0.5, 7.*64. - 0.5);
    Add2RecPointsList(h2, kTriggerDigitsBendPlane);

    TH3F *h3 = new TH3F("hTriggerDigitsNonBendPlane", "Trigger digits in non-bending plane",
			4, -0.5, 4. - 0.5,
			18, -0.5, 18. - 0.5,
			112, -0.5, 112. - 0.5);
    Add2RecPointsList(h3, kTriggerDigitsNonBendPlane);
}


//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitESDs()
{
    ///create ESDs histograms in ESDs subdir
  TH1F* h0 = new TH1F("hESDnTracks", "ESDs track number distribution", 30, 0., 30.);  
  Add2ESDsList(h0, kESDnTracks);

  TH1F* h1 = new TH1F("hESDMomentum", "ESDs P distribution", 300, 0., 300) ; 
  Add2ESDsList(h1, kESDMomentum);

  TH1F* h2 = new TH1F("hESDPt", "ESDs Pt distribution", 200, 0., 50) ; 
  Add2ESDsList(h2, kESDPt);

  TH1F* h3 = new TH1F("hESDRapidity", "ESDs rapidity distribution", 200, -4.5,-2.) ; 
  Add2ESDsList(h3, kESDRapidity);

  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); ++i) 
  {
    TH2F* h4 = new TH2F(Form("%s%d", "hESDClusterHitMap", i), 
		     Form("%s %d", "ESD Clusters hit distribution for chamber", i),
		     100, -1*AliMUONConstants::Rmax(i/2), AliMUONConstants::Rmax(i/2),
		     100, -1*AliMUONConstants::Rmax(i/2), AliMUONConstants::Rmax(i/2)); 
    Add2ESDsList(h4, kESDClusterHitMap+i);
  }
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
    /// make QA for rawdata
    Int_t busPatchId;
    UShort_t manuId;  
    UChar_t channelId;
    UShort_t charge;

    rawReader->Reset();
    AliMUONRawStreamTracker rawStream(rawReader);
    rawStream.First();
    while( rawStream.Next(busPatchId, manuId, channelId, charge) ) {
  
      GetRawsData(kRawBusPatch)->Fill(busPatchId);
      GetRawsData(kRawCharge)->Fill(charge);
			Int_t iDDL = rawStream.GetCurentDDL();
			GetRawsData(kRawBuspatchDDL + iDDL)->Fill(AliMpBusPatch::GetLocalBusID(busPatchId, iDDL));
		
		  
    } // Next digit
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRecPoints(TTree* clustersTree)
{
  
    /// makes data from trigger response
      
    // Fired pads info
    fDigitStore->Clear();

    if (!fTriggerStore) fTriggerStore = AliMUONVTriggerStore::Create(*clustersTree);
    fTriggerStore->Clear();
    fTriggerStore->Connect(*clustersTree, false);
    clustersTree->GetEvent(0);

    AliMUONLocalTrigger* locTrg;
    TIter nextLoc(fTriggerStore->CreateLocalIterator());

    while ( ( locTrg = static_cast<AliMUONLocalTrigger*>(nextLoc()) ) ) 
    {
      if (locTrg->IsNull()) continue;
   
      TArrayS xyPattern[2];
      locTrg->GetXPattern(xyPattern[0]);
      locTrg->GetYPattern(xyPattern[1]);

      Int_t nBoard = locTrg->LoCircuit();
       fDigitMaker->TriggerDigits(nBoard, xyPattern, *fDigitStore);
    }

    TIter nextDigit(fDigitStore->CreateIterator());
    AliMUONVDigit* mDigit;
    while ( ( mDigit = static_cast<AliMUONVDigit*>(nextDigit()) ) )
    {
      Int_t detElemId = mDigit->DetElemId();
      Int_t ch = detElemId/100 - 11;
      Int_t slat = detElemId%100;
      Int_t cathode = mDigit->Cathode();
      Int_t ix = mDigit->PadX();
      Int_t iy = mDigit->PadY();
      Int_t maxY = (cathode==0) ? 64 : 1;
      Int_t currPair = ix*maxY + iy;
      ERecPoints hindex 
        = ( cathode == 0 ) ? kTriggerDigitsBendPlane : kTriggerDigitsNonBendPlane;
      
      ((TH3F*)GetRecPointsData(hindex))->Fill(ch, slat, currPair);
    }
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeESDs(AliESDEvent* esd)
{
    /// make QA data from ESDs

    TLorentzVector v1;

    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ; 
    GetESDsData(0)->Fill(nTracks);

    for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {

      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);
      muonTrack->LorentzP(v1);

      GetESDsData(1)->Fill(v1.P());
      GetESDsData(2)->Fill(v1.Pt());
      GetESDsData(3)->Fill(v1.Rapidity());

      TClonesArray clusters =  muonTrack->GetClusters();

      for (Int_t iCluster = 0; iCluster <clusters.GetEntriesFast(); ++iCluster) {
	AliESDMuonCluster* cluster = (AliESDMuonCluster*)clusters.At(iCluster);
	GetESDsData(kESDClusterHitMap+cluster->GetChamberId())
          ->Fill(cluster->GetX(), cluster->GetY());
      }
    }
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::StartOfDetectorCycle()
{
    /// Detector specific actions at start of cycle
  
}
