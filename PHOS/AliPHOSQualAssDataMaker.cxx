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


/* $Id$ */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliPHOSDigit.h"
#include "AliPHOSHit.h"
#include "AliPHOSQualAssDataMaker.h"
#include "AliPHOSCpvRecPoint.h" 
#include "AliPHOSEmcRecPoint.h" 
#include "AliPHOSRecParticle.h" 
#include "AliPHOSTrackSegment.h" 

ClassImp(AliPHOSQualAssDataMaker)
           
//____________________________________________________________________________ 
  AliPHOSQualAssDataMaker::AliPHOSQualAssDataMaker() : 
  AliQualAssDataMaker(AliQualAss::GetDetName(AliQualAss::kPHOS), "PHOS Quality Assurance Data Maker"),
  fhHits(0x0), 
  fhHitsMul(0x0), 
  fhDigits(0x0),
  fhDigitsMul(0x0),
  fhSDigits(0x0),
  fhSDigitsMul(0x0),
  fhEmcRecPoints(0x0),
  fhEmcRecPointsMul(0x0),
  fhCpvRecPoints(0x0),
  fhCpvRecPointsMul(0x0),
  fhTrackSegments(0x0),
  fhTrackSegmentsMul(0x0),
  fhRecParticles(0x0),
  fhRecParticlesMul(0x0),
  fhESDs(0x0),
  fhESDsMul(0x0) 
{
  // ctor
  fDetectorDir = fOutput->GetDirectory(GetName()) ;  
  if (!fDetectorDir) 
    fDetectorDir = fOutput->mkdir(GetName()) ;  
}

//____________________________________________________________________________ 
AliPHOSQualAssDataMaker::AliPHOSQualAssDataMaker(const AliPHOSQualAssDataMaker& qadm) :
  AliQualAssDataMaker(), 
  fhHits(qadm.fhHits), 
  fhHitsMul(qadm.fhHitsMul), 
  fhDigits(qadm.fhDigits),
  fhDigitsMul(qadm.fhDigitsMul),
  fhSDigits(qadm.fhSDigits),
  fhSDigitsMul(qadm.fhSDigitsMul), 
  fhEmcRecPoints(qadm.fhEmcRecPoints),
  fhEmcRecPointsMul(qadm.fhEmcRecPointsMul), 
  fhCpvRecPoints(qadm.fhCpvRecPoints),
  fhCpvRecPointsMul(qadm.fhCpvRecPointsMul), 
  fhTrackSegments(qadm.fhTrackSegments),
  fhTrackSegmentsMul(qadm.fhTrackSegmentsMul), 
  fhRecParticles(qadm.fhRecParticles),
  fhRecParticlesMul(qadm.fhRecParticlesMul), 
  fhESDs(qadm.fhESDs), 
  fhESDsMul(qadm.fhESDsMul) 
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliPHOSQualAssDataMaker& AliPHOSQualAssDataMaker::operator = (const AliPHOSQualAssDataMaker& qadm )
{
  // Equal operator.
  this->~AliPHOSQualAssDataMaker();
  new(this) AliPHOSQualAssDataMaker(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitESDs()
{
  //create ESDs histograms in ESDs subdir
  fhESDs     = new TH1F("hPhosESDs",    "ESDs energy distribution in PHOS",       100, 0., 100.) ; 
  fhESDs->Sumw2() ; 
  fhESDsMul  = new TH1I("hPhosESDsMul", "ESDs multiplicity distribution in PHOS", 100, 0., 100) ; 
  fhESDsMul->Sumw2() ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitHits()
{
  // create Hits histograms in Hits subdir
  fhHits     = new TH1F("hPhosHits",    "Hits energy distribution in PHOS",       100, 0., 100.) ; 
  fhHits->Sumw2() ;
  fhHitsMul  = new TH1I("hPhosHitsMul", "Hits multiplicity distribution in PHOS", 500, 0., 10000) ; 
  fhHitsMul->Sumw2() ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitDigits()
{
  // create Digits histograms in Digits subdir
  fhDigits     = new TH1I("hPhosDigits",    "Digits amplitude distribution in PHOS",    500, 0, 5000) ; 
  fhDigits->Sumw2() ;
  fhDigitsMul  = new TH1I("hPhosDigitsMul", "Digits multiplicity distribution in PHOS", 500, 0, 1000) ; 
  fhDigitsMul->Sumw2() ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitRecParticles()
{
  // create Reconstructed particles histograms in RecParticles subdir
  fhRecParticles     = new TH1F("hPhosRecParticles",    "RecParticles energy distribution in PHOS",       100, 0., 100.) ; 
  fhRecParticles->Sumw2() ;
  fhRecParticlesMul  = new TH1I("hPhosRecParticlesMul", "RecParticles multiplicity distribution in PHOS", 100, 0,  100) ; 
  fhRecParticlesMul->Sumw2() ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  fhEmcRecPoints     = new TH1F("hEmcPhosRecPoints",    "EMCA RecPoints energy distribution in PHOS",       100, 0., 100.) ; 
  fhEmcRecPoints->Sumw2() ;
  fhEmcRecPointsMul  = new TH1I("hEmcPhosRecPointsMul", "EMCA RecPoints multiplicity distribution in PHOS", 100, 0,  100) ; 
  fhEmcRecPointsMul->Sumw2() ;

  fhCpvRecPoints     = new TH1F("hCpvPhosRecPoints",    "CPV RecPoints energy distribution in PHOS",       100, 0., 100.) ; 
  fhCpvRecPoints->Sumw2() ;
  fhCpvRecPointsMul  = new TH1I("hCpvPhosRecPointsMul", "CPV RecPoints multiplicity distribution in PHOS", 100, 0,  100) ; 
  fhCpvRecPointsMul->Sumw2() ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  fhSDigits     = new TH1F("hPhosSDigits",    "SDigits energy distribution in PHOS",       100, 0., 100.) ; 
  fhSDigits->Sumw2() ;
  fhSDigitsMul  = new TH1I("hPhosSDigitsMul", "SDigits multiplicity distribution in PHOS", 500, 0,  10000) ; 
  fhSDigitsMul->Sumw2() ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitTrackSegments()
{
  // create Track Segments histograms in TrackSegments subdir
  fhTrackSegments     = new TH1F("hPhosTrackSegments",    "TrackSegments EMC-CPV distance in PHOS",       500, 0., 5000.) ; 
  fhTrackSegments->Sumw2() ;
  fhTrackSegmentsMul  = new TH1I("hPhosTrackSegmentsMul", "TrackSegments multiplicity distribution in PHOS", 100, 0,  100) ; 
  fhTrackSegmentsMul->Sumw2() ;
}

//____________________________________________________________________________
void AliPHOSQualAssDataMaker::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  
  Int_t maxClu = esd->GetNumberOfPHOSClusters() ; 
  Int_t index = 0, count = 0 ; 
  for ( index = 0 ; index < maxClu; index++ ) {
    AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
    fhESDs->Fill(clu->E()) ;
    count++ ; 
  }
  fhESDsMul->Fill(count) ;
}

//____________________________________________________________________________
void AliPHOSQualAssDataMaker::MakeHits(TObject * data)
{
  //make QA data from Hits

  TClonesArray * hits = dynamic_cast<TClonesArray *>(data) ; 
  if (!hits) {
    AliError("Wrong type of hits container") ; 
  } else {
    fhHitsMul->Fill(hits->GetEntriesFast()) ; 
    TIter next(hits) ; 
    AliPHOSHit * hit ; 
    while ( (hit = dynamic_cast<AliPHOSHit *>(next())) ) {
      fhHits->Fill( hit->GetEnergy()) ;
    }
  } 
}
 
//____________________________________________________________________________
void AliPHOSQualAssDataMaker::MakeDigits(TObject * data)
{
  // makes data from Digits

  TClonesArray * digits = dynamic_cast<TClonesArray *>(data) ; 
  if (!digits) {
    AliError("Wrong type of digits container") ; 
  } else {
    fhDigitsMul->Fill(digits->GetEntriesFast()) ; 
    TIter next(digits) ; 
    AliPHOSDigit * digit ; 
    while ( (digit = dynamic_cast<AliPHOSDigit *>(next())) ) {
      fhDigits->Fill( digit->GetEnergy()) ;
    }  
  }
}

//____________________________________________________________________________
// void AliPHOSQualAssDataMaker::MakeRecParticles(TTree * recpar)
// {
//   // makes data from RecParticles

//   TClonesArray * recparticles = dynamic_cast<TClonesArray*>(fData) ; 
//   fhRecParticlesMul->Fill(recparticles->GetEntriesFast()) ; 
//   TIter next(recparticles) ; 
//   AliPHOSRecParticle * recparticle ; 
//   while ( (recparticle = dynamic_cast<AliPHOSRecParticle *>(next())) ) {
//     fhRecParticles->Fill( recparticle->Energy()) ;
//   }
// }

//____________________________________________________________________________
void AliPHOSQualAssDataMaker::MakeRecPoints(TTree * clustersTree)
{
  {
    // makes data from RecPoints
    TBranch *emcbranch = clustersTree->GetBranch("PHOSEmcRP");
    if (!emcbranch) { 
      AliError("can't get the branch with the PHOS EMC clusters !");
      return;
    }
    TObjArray * emcrecpoints = new TObjArray(100) ;
    emcbranch->SetAddress(&emcrecpoints);
    emcbranch->GetEntry(0);
    
    fhEmcRecPointsMul->Fill(emcrecpoints->GetEntriesFast()) ; 
    TIter next(emcrecpoints) ; 
    AliPHOSEmcRecPoint * rp ; 
    while ( (rp = dynamic_cast<AliPHOSEmcRecPoint *>(next())) ) {
      fhEmcRecPoints->Fill( rp->GetEnergy()) ;
    }
    emcrecpoints->Delete();
    delete emcrecpoints;
  }
  {
    TBranch *cpvbranch = clustersTree->GetBranch("PHOSCpvRP");
    if (!cpvbranch) { 
      AliError("can't get the branch with the PHOS CPV clusters !");
      return;
    }
    TObjArray *cpvrecpoints = new TObjArray(100) ;
    cpvbranch->SetAddress(&cpvrecpoints);
    cpvbranch->GetEntry(0);
    
    fhCpvRecPointsMul->Fill(cpvrecpoints->GetEntriesFast()) ; 
    TIter next(cpvrecpoints) ; 
    AliPHOSCpvRecPoint * rp ; 
    while ( (rp = dynamic_cast<AliPHOSCpvRecPoint *>(next())) ) {
      fhCpvRecPoints->Fill( rp->GetEnergy()) ;
    }
    cpvrecpoints->Delete();
    delete cpvrecpoints;
  }
}

//____________________________________________________________________________
void AliPHOSQualAssDataMaker::MakeSDigits(TObject * data)
{
  // makes data from SDigits
  
  TClonesArray * sdigits = dynamic_cast<TClonesArray *>(data) ; 
  if (!sdigits) {
    AliError("Wrong type of sdigits container") ; 
  } else {
    fhSDigitsMul->Fill(sdigits->GetEntriesFast()) ; 
    TIter next(sdigits) ; 
    AliPHOSDigit * sdigit ; 
    while ( (sdigit = dynamic_cast<AliPHOSDigit *>(next())) ) {
      fhSDigits->Fill( sdigit->GetEnergy()) ;
    } 
  }
}

//____________________________________________________________________________
// void AliPHOSQualAssDataMaker::MakeTrackSegments(TTree * ts)
// {
//   // makes data from TrackSegments

//   TClonesArray * tracksegments = dynamic_cast<TClonesArray*>(fData) ;

//   fhTrackSegmentsMul->Fill(tracksegments->GetEntriesFast()) ; 
//   TIter next(tracksegments) ; 
//   AliPHOSTrackSegment * ts ; 
//   while ( (ts = dynamic_cast<AliPHOSTrackSegment *>(next())) ) {
//     fhTrackSegments->Fill( ts->GetCpvDistance()) ;
//   } 
// }
