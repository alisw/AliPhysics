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
#include <TH2F.h> 

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
#include "AliPHOSRawDecoder.h"

ClassImp(AliPHOSQualAssDataMaker)
           
//____________________________________________________________________________ 
  AliPHOSQualAssDataMaker::AliPHOSQualAssDataMaker() : 
  AliQualAssDataMaker(AliQualAss::GetDetName(AliQualAss::kPHOS), "PHOS Quality Assurance Data Maker")
{
  // ctor
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
void AliPHOSQualAssDataMaker::EndOfDetectorCycle()
{
  //Detector specific actions at end of cycle
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
  TH1I * h0 = new TH1I("hPhosModules",    "Hits in EMCA PHOS modules",       6, 0, 6) ; 
  h0->Sumw2() ;
  Add2RawsList(h0, 0) ;
  const Int_t modMax = 5 ; 
  TH2I * h1[modMax] ; 
  for (Int_t mod = 0; mod < modMax; mod++) {
   char name[16] ; 
   sprintf(name, "hPHOSMod%d", mod) ; 
   char title[32] ; 
   sprintf(title, "Raws x Columns for PHOS module %d", mod+1) ;  
   h1[mod] = new TH2I(name, title, 64, 0, 63, 56, 0, 55) ; 
   Add2RawsList(h1[mod], mod+1) ;
  }
  TH1F * h6 = new TH1F("hPhosRawtime", "Time of raw hits in PHOS", 100, 0, 100.) ; 
  h6->Sumw2() ;
  Add2RawsList(h6, 6) ;
  TH1F * h7 = new TH1F("hPhosRawEnergy", "Energy of raw hits in PHOS", 1000, 0., 1200.) ; 
  h7->Sumw2() ;
  Add2RawsList(h7, 7) ;
 
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
void AliPHOSQualAssDataMaker::MakeHits(TClonesArray * hits)
{
  //make QA data from Hits

    GetHitsData(1)->Fill(hits->GetEntriesFast()) ; 
    TIter next(hits) ; 
    AliPHOSHit * hit ; 
    while ( (hit = dynamic_cast<AliPHOSHit *>(next())) ) {
      GetHitsData(0)->Fill( hit->GetEnergy()) ;
    }
}
//____________________________________________________________________________
void AliPHOSQualAssDataMaker::MakeDigits(TClonesArray * digits)
{
  // makes data from Digits

    GetDigitsData(1)->Fill(digits->GetEntriesFast()) ; 
    TIter next(digits) ; 
    AliPHOSDigit * digit ; 
    while ( (digit = dynamic_cast<AliPHOSDigit *>(next())) ) {
      GetDigitsData(0)->Fill( digit->GetEnergy()) ;
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
void AliPHOSQualAssDataMaker::MakeRaws(AliRawReader* rawReader)
{
  rawReader->Reset() ; 
  AliPHOSRawDecoder decoder(rawReader);
  decoder.SetOldRCUFormat(kTRUE);
  Int_t count = 0 ; 
  while (decoder.NextDigit()) {
   Int_t module  = decoder.GetModule() ;
   Int_t row     = decoder.GetRow() ;
   Int_t col     = decoder.GetColumn() ;
   Double_t time = decoder.GetTime() ;
   Double_t energy  = decoder.GetEnergy() ;          
   GetRawsData(0)->Fill(module) ; 
   GetRawsData(module)->Fill(row, col) ; 
   GetRawsData(6)->Fill(time) ; 
   GetRawsData(7)->Fill(energy) ; 
   AliInfo(Form(" %d %d %d %d %f %f\n", count, module, row, col, time, energy)) ;
   count++ ; 
  } 
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
void AliPHOSQualAssDataMaker::MakeSDigits(TClonesArray * sdigits)
{
  // makes data from SDigits
  
	GetSDigitsData(1)->Fill(sdigits->GetEntriesFast()) ; 
    TIter next(sdigits) ; 
    AliPHOSDigit * sdigit ; 
    while ( (sdigit = dynamic_cast<AliPHOSDigit *>(next())) ) {
      GetSDigitsData(0)->Fill( sdigit->GetEnergy()) ;
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

//____________________________________________________________________________ 
void AliPHOSQualAssDataMaker::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
