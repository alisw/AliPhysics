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
  AliQualAssDataMaker(AliQualAss::GetDetName(AliQualAss::kPHOS), "PHOS Quality Assurance Data Maker")
{
  // ctor
  fDetectorDir = fOutput->GetDirectory(GetName()) ;  
  if (!fDetectorDir) 
    fDetectorDir = fOutput->mkdir(GetName()) ;  
}

//____________________________________________________________________________ 
AliPHOSQualAssDataMaker::AliPHOSQualAssDataMaker(const AliPHOSQualAssDataMaker& qadm) :
  AliQualAssDataMaker()
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
  TH1F * h0 = new TH1F("hPhosESDs",    "ESDs energy distribution in PHOS",       100, 0., 100.) ;  
  h0->Sumw2() ; 
  Add2ESDsList(h0, 0) ;
  TH1I * h1  = new TH1I("hPhosESDsMul", "ESDs multiplicity distribution in PHOS", 100, 0., 100) ; 
  h1->Sumw2() ;
  Add2ESDsList(h1, 1) ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitHits()
{
  // create Hits histograms in Hits subdir
  TH1F * h0 = new TH1F("hPhosHits",    "Hits energy distribution in PHOS",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2HitsList(h0, 0) ;
  TH1I * h1  = new TH1I("hPhosHitsMul", "Hits multiplicity distribution in PHOS", 500, 0., 10000) ; 
  h1->Sumw2() ;
  Add2HitsList(h1, 1) ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitDigits()
{
  // create Digits histograms in Digits subdir
  TH1I * h0 = new TH1I("hPhosDigits",    "Digits amplitude distribution in PHOS",    500, 0, 5000) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0) ;
  TH1I * h1 = new TH1I("hPhosDigitsMul", "Digits multiplicity distribution in PHOS", 500, 0, 1000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, 0) ;
}

//____________________________________________________________________________ 
//void AliPHOSQualAssDataMaker::InitRecParticles()
//{
//  // create Reconstructed particles histograms in RecParticles subdir
//  fhRecParticles     = new TH1F("hPhosRecParticles",    "RecParticles energy distribution in PHOS",       100, 0., 100.) ; 
//  fhRecParticles->Sumw2() ;
//  fhRecParticlesMul  = new TH1I("hPhosRecParticlesMul", "RecParticles multiplicity distribution in PHOS", 100, 0,  100) ; 
//  fhRecParticlesMul->Sumw2() ;
//}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  TH1F * h0 = new TH1F("hEmcPhosRecPoints",    "EMCA RecPoints energy distribution in PHOS",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2RecPointsList(h0, 0) ;
  TH1I * h1 = new TH1I("hEmcPhosRecPointsMul", "EMCA RecPoints multiplicity distribution in PHOS", 100, 0,  100) ; 
  h1->Sumw2() ;
  Add2RecPointsList(h1, 1) ;

  TH1F * h2 = new TH1F("hCpvPhosRecPoints",    "CPV RecPoints energy distribution in PHOS",       100, 0., 100.) ; 
  h2->Sumw2() ;
  Add2RecPointsList(h2, 2) ;
  TH1I * h3 = new TH1I("hCpvPhosRecPointsMul", "CPV RecPoints multiplicity distribution in PHOS", 100, 0,  100) ; 
  h3->Sumw2() ;
  Add2RecPointsList(h3, 3) ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitRaws()
{
  // create Raws histograms in Raws subdir
  TH1F * h0 = new TH1F("hEmcPhosRaws",    "EMCA Raws in PHOS",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2RawsList(h0, 0) ;
}

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  TH1F * h0 = new TH1F("hPhosSDigits",    "SDigits energy distribution in PHOS",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2SDigitsList(h0, 0) ;
  TH1I * h1 = new TH1I("hPhosSDigitsMul", "SDigits multiplicity distribution in PHOS", 500, 0,  10000) ; 
  h1->Sumw2() ;
  Add2SDigitsList(h1, 1) ;
}

//____________________________________________________________________________ 
//void AliPHOSQualAssDataMaker::InitTrackSegments()
//{
//  // create Track Segments histograms in TrackSegments subdir
//  fhTrackSegments     = new TH1F("hPhosTrackSegments",    "TrackSegments EMC-CPV distance in PHOS",       500, 0., 5000.) ; 
//  fhTrackSegments->Sumw2() ;
//  fhTrackSegmentsMul  = new TH1I("hPhosTrackSegmentsMul", "TrackSegments multiplicity distribution in PHOS", 100, 0,  100) ; 
//  fhTrackSegmentsMul->Sumw2() ;
//}

//____________________________________________________________________________
void AliPHOSQualAssDataMaker::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  
  Int_t maxClu = esd->GetNumberOfPHOSClusters() ; 
  Int_t index = 0, count = 0 ; 
  for ( index = 0 ; index < maxClu; index++ ) {
    AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
    GetESDsData(0)->Fill(clu->E()) ;
    count++ ; 
  }
  GetESDsData(1)->Fill(count) ;
}

//____________________________________________________________________________
void AliPHOSQualAssDataMaker::MakeHits(TObject * data)
{
  //make QA data from Hits

  TClonesArray * hits = dynamic_cast<TClonesArray *>(data) ; 
  if (!hits) {
    AliError("Wrong type of hits container") ; 
  } else {
    GetHitsData(1)->Fill(hits->GetEntriesFast()) ; 
    TIter next(hits) ; 
    AliPHOSHit * hit ; 
    while ( (hit = dynamic_cast<AliPHOSHit *>(next())) ) {
      GetHitsData(0)->Fill( hit->GetEnergy()) ;
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
    GetDigitsData(1)->Fill(digits->GetEntriesFast()) ; 
    TIter next(digits) ; 
    AliPHOSDigit * digit ; 
    while ( (digit = dynamic_cast<AliPHOSDigit *>(next())) ) {
      GetDigitsData(0)->Fill( digit->GetEnergy()) ;
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
void AliPHOSQualAssDataMaker::MakeRaws(TTree * clustersTree)
{
    GetRawsData(1)->Fill(99) ; 
}

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
    
    GetRecPointsData(1)->Fill(emcrecpoints->GetEntriesFast()) ; 
    TIter next(emcrecpoints) ; 
    AliPHOSEmcRecPoint * rp ; 
    while ( (rp = dynamic_cast<AliPHOSEmcRecPoint *>(next())) ) {
      GetRecPointsData(0)->Fill( rp->GetEnergy()) ;
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
    
    GetRecPointsData(1)->Fill(cpvrecpoints->GetEntriesFast()) ; 
    TIter next(cpvrecpoints) ; 
    AliPHOSCpvRecPoint * rp ; 
    while ( (rp = dynamic_cast<AliPHOSCpvRecPoint *>(next())) ) {
      GetRecPointsData(0)->Fill( rp->GetEnergy()) ;
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
    GetSDigitsData(1)->Fill(sdigits->GetEntriesFast()) ; 
    TIter next(sdigits) ; 
    AliPHOSDigit * sdigit ; 
    while ( (sdigit = dynamic_cast<AliPHOSDigit *>(next())) ) {
      GetSDigitsData(0)->Fill( sdigit->GetEnergy()) ;
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
