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
#include <TLorentzVector.h>

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliQAChecker.h"
#include "AliMUONCluster.h"  
#include "AliMUONRawStreamTracker.h"
#include "AliMUONRawStreamTrigger.h"

#include "AliMUONVClusterStore.h"
#include "AliMUONVCluster.h"
#include "AliESDMuonTrack.h"

#include "AliMUONQADataMaker.h"

//-----------------------------------------------------------------------------
/// \class AliMUONQADataMaker
///
/// MUON base class for quality assurance data (histo) maker
///
/// \author C. Finck

/// \cond CLASSIMP
ClassImp(AliMUONQADataMaker)
/// \endcond
           
//____________________________________________________________________________ 
AliMUONQADataMaker::AliMUONQADataMaker() : 
    AliQADataMaker(AliQA::GetDetName(AliQA::kMUON), "MUON Quality Assurance Data Maker"),
    fClusterStore(0x0)
{
    /// ctor
}

//____________________________________________________________________________ 
AliMUONQADataMaker::AliMUONQADataMaker(const AliMUONQADataMaker& qadm) :
    AliQADataMaker()
{
    ///copy ctor 
    SetName((const char*)qadm.GetName()) ; 
    SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliMUONQADataMaker& AliMUONQADataMaker::operator = (const AliMUONQADataMaker& qadm )
{
    /// Equal operator.
    this->~AliMUONQADataMaker();
    new(this) AliMUONQADataMaker(qadm);
    return *this;
}

//__________________________________________________________________
AliMUONQADataMaker::~AliMUONQADataMaker()
{
    /// dtor
    delete fClusterStore;
}

//____________________________________________________________________________ 
void AliMUONQADataMaker::EndOfDetectorCycle(AliQA::TASKINDEX task, TList* list)
{
    ///Detector specific actions at end of cycle
    // do the QA checking
    AliQAChecker::Instance()->Run(AliQA::kMUON, task, list) ;  
}

//____________________________________________________________________________ 
void AliMUONQADataMaker::InitRaws()
{
    /// create Raws histograms in Raws subdir
    TH1I* h0 = new TH1I("hRawBusPatch", "buspatch distribution",  1932, 1, 1932); 
    Add2RawsList(h0, 0);

    TH1I* h1 = new TH1I("hRawCharge", "Charge distribution in rawdata", 4096, 0, 4095); 
    Add2RawsList(h1, 1);

}

//____________________________________________________________________________ 
void AliMUONQADataMaker::InitRecPoints()
{
    /// create Reconstructed Points histograms in RecPoints subdir
    TH1F* h0 = new TH1F("hClusterCharge", "Clusters Charge distribution", 1000, 0., 4095.); 
    Add2RecPointsList(h0, 0);

    TH1I* h1 = new TH1I("hClusterDetElem", "DetElemId distribution in Clusters ", 1000, 100., 1100.); 
    Add2RecPointsList(h1, 1);
}


//____________________________________________________________________________ 
void AliMUONQADataMaker::InitESDs()
{
    ///create ESDs histograms in ESDs subdir
  TH1F* h0 = new TH1F("hESDnTracks", "ESDs track number distribution", 30, 0., 30.);  
  Add2ESDsList(h0, 0);

  TH1F* h1 = new TH1F("hESDMomentum", "ESDs P distribution", 300, 0., 300) ; 
  Add2ESDsList(h1, 1);

  TH1F* h2 = new TH1F("hESDPt", "ESDs Pt distribution", 200, 0., 50) ; 
  Add2ESDsList(h2, 2) ;

  TH1F* h3 = new TH1F("hESDRapidity", "ESDs rapidity distribution", 200, -4.5,-2.) ; 
  Add2ESDsList(h3, 3) ;
}

//____________________________________________________________________________
void AliMUONQADataMaker::MakeRaws(AliRawReader* rawReader)
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
  
      GetRawsData(0)->Fill(busPatchId);
      GetRawsData(1)->Fill(charge);
		  
    } // Next digit

}

//____________________________________________________________________________
void AliMUONQADataMaker::MakeRecPoints(TTree* clustersTree)
{
  
    /// makes data from RecPoints
    if (!fClusterStore)
	fClusterStore = AliMUONVClusterStore::Create(*clustersTree);
    fClusterStore->Connect(*clustersTree, false);
    clustersTree->GetEvent(0);
    
    TIter next(fClusterStore->CreateIterator());

    AliMUONVCluster* clus = 0x0;

    while ( ( clus = static_cast<AliMUONVCluster*>(next()) ) )
    {
      GetRecPointsData(0)->Fill(clus->GetCharge());
      GetRecPointsData(1)->Fill(clus->GetDetElemId());
    }
}

//____________________________________________________________________________
void AliMUONQADataMaker::MakeESDs(AliESDEvent* esd)
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
    }
}

//____________________________________________________________________________ 
void AliMUONQADataMaker::StartOfDetectorCycle()
{
    /// Detector specific actions at start of cycle
  
}
