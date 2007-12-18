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
#include "AliPHOSQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliPHOSCpvRecPoint.h" 
#include "AliPHOSEmcRecPoint.h" 
#include "AliPHOSRecParticle.h" 
#include "AliPHOSTrackSegment.h" 
#include "AliPHOSRawDecoder.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSRecoParam.h"

ClassImp(AliPHOSQADataMakerRec)
           
//____________________________________________________________________________ 
  AliPHOSQADataMakerRec::AliPHOSQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kPHOS), "PHOS Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliPHOSQADataMakerRec::AliPHOSQADataMakerRec(const AliPHOSQADataMakerRec& qadm) :
  AliQADataMakerRec()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliPHOSQADataMakerRec& AliPHOSQADataMakerRec::operator = (const AliPHOSQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliPHOSQADataMakerRec();
  new(this) AliPHOSQADataMakerRec(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kPHOS, task, list) ;  
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitESDs()
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
//void AliPHOSQADataMakerRec::InitRecParticles()
//{
//  // create Reconstructed particles histograms in RecParticles subdir
//  fhRecParticles     = new TH1F("hPhosRecParticles",    "RecParticles energy distribution in PHOS",       100, 0., 100.) ; 
//  fhRecParticles->Sumw2() ;
//  fhRecParticlesMul  = new TH1I("hPhosRecParticlesMul", "RecParticles multiplicity distribution in PHOS", 100, 0,  100) ; 
//  fhRecParticlesMul->Sumw2() ;
//}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitRecPoints()
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
void AliPHOSQADataMakerRec::InitRaws()
{
  // create Raws histograms in Raws subdir
  const Int_t modMax = 5 ; 
  TH2I * h0[modMax*2] ; 
  char name[32] ; 
  char title[32] ; 
  for (Int_t mod = 0; mod < modMax; mod++) {
   sprintf(title, "Low Gain Rows x Columns for PHOS module %d", mod) ;  
   sprintf(name, "hLowPHOSxyMod%d", mod) ; 
   h0[mod] = new TH2I(name, title, 64, 1, 65, 56, 1, 57) ; 
   Add2RawsList(h0[mod], mod) ;
   sprintf(title, "High Gain Rows x Columns for PHOS module %d", mod) ;  
   sprintf(name, "hHighPHOSxyMod%d", mod) ; 
   h0[mod+modMax] = new TH2I(name, title, 64, 1, 65, 56, 1, 57) ; 
   Add2RawsList(h0[mod+modMax], mod+modMax) ;
  }
  TH1I * h10 = new TH1I("hLowPhosModules",    "Low Gain Hits in EMCA PHOS modules",       6, 0, 6) ; 
  h10->Sumw2() ;
  Add2RawsList(h10, 10) ;
  TH1I * h11 = new TH1I("hHighPhosModules",    "High Gain Hits in EMCA PHOS modules",       6, 0, 6) ; 
  h11->Sumw2() ;
  Add2RawsList(h11, 11) ;
  TH1F * h12 = new TH1F("hLowPhosRawtime", "Low Gain Time of raw hits in PHOS", 100, 0, 100.) ; 
  h12->Sumw2() ;
  Add2RawsList(h12, 12) ;
  TH1F * h13 = new TH1F("hHighPhosRawtime", "High Gain Time of raw hits in PHOS", 100, 0, 100.) ; 
  h13->Sumw2() ;
  Add2RawsList(h13, 13) ;
  TH1F * h14 = new TH1F("hLowPhosRawEnergy", "Low Gain Energy of raw hits in PHOS", 100, 0., 100.) ; 
  h14->Sumw2() ;
  Add2RawsList(h14, 14) ;
  TH1F * h15 = new TH1F("hHighPhosRawEnergy", "High Gain Energy of raw hits in PHOS", 100, 0., 100.) ; 
  h15->Sumw2() ;
  Add2RawsList(h15, 15) ;
 
}

//____________________________________________________________________________ 
//void AliPHOSQADataMakerRec::InitTrackSegments()
//{
//  // create Track Segments histograms in TrackSegments subdir
//  fhTrackSegments     = new TH1F("hPhosTrackSegments",    "TrackSegments EMC-CPV distance in PHOS",       500, 0., 5000.) ; 
//  fhTrackSegments->Sumw2() ;
//  fhTrackSegmentsMul  = new TH1I("hPhosTrackSegmentsMul", "TrackSegments multiplicity distribution in PHOS", 100, 0,  100) ; 
//  fhTrackSegmentsMul->Sumw2() ;
//}

//____________________________________________________________________________
void AliPHOSQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs

  Int_t count = 0 ; 
  for ( Int_t index = 0; index < esd->GetNumberOfCaloClusters() ; index++ ) {
	AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
	if ( clu->IsPHOS() ) {
		GetESDsData(0)->Fill(clu->E()) ;
		count++ ;
	} 
  }
  GetESDsData(1)->Fill(count) ;
}

//____________________________________________________________________________
// void AliPHOSQADataMakerRec::MakeRecParticles(TTree * recpar)
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
void AliPHOSQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  const Int_t modMax = 5 ; 
  rawReader->Reset() ; 
  AliPHOSRawDecoder decoder(rawReader);
  decoder.SetOldRCUFormat(kFALSE);
  decoder.SubtractPedestals(AliPHOSReconstructor::GetRecoParamEmc()->SubtractPedestals());
  Int_t count = 0 ; 
  while (decoder.NextDigit()) {
   Int_t module  = decoder.GetModule() ;
   Int_t row     = decoder.GetRow() ;
   Int_t col     = decoder.GetColumn() ;
   Double_t time = decoder.GetTime() ;
   Double_t energy  = decoder.GetEnergy() ;     
   Bool_t lowGain = decoder.IsLowGain();
   if (lowGain) {
     GetRawsData(module)->Fill(row, col) ; 
	 GetRawsData(10)->Fill(module) ; 
     GetRawsData(12)->Fill(time) ; 
     GetRawsData(14)->Fill(energy) ; 
   } else {
	 GetRawsData(module+modMax)->Fill(row, col) ; 
	 GetRawsData(11)->Fill(module) ; 
     GetRawsData(13)->Fill(time) ; 
     GetRawsData(15)->Fill(energy) ; 
   }
   //AliInfo(Form(" %d %d %d %d %f %f\n", count, module, row, col, time, energy)) ;
   count++ ; 
  } 
}

//____________________________________________________________________________
void AliPHOSQADataMakerRec::MakeRecPoints(TTree * clustersTree)
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
// void AliPHOSQADataMakerRec::MakeTrackSegments(TTree * ts)
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
void AliPHOSQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
