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
/*
Based on the QA code for PHOS written by Yves Schutz July 2007

Authors:  J.Klay (Cal Poly) May 2008
          S. Salur LBL April 2008

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
#include "AliESDCaloCells.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliEMCALQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliEMCALRecPoint.h" 
#include "AliEMCALRawUtils.h"
#include "AliEMCALReconstructor.h"
#include "AliEMCALRecParam.h"
#include "AliRawReader.h"

ClassImp(AliEMCALQADataMakerRec)
           
//____________________________________________________________________________ 
  AliEMCALQADataMakerRec::AliEMCALQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kEMCAL), "EMCAL Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliEMCALQADataMakerRec::AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) :
  AliQADataMakerRec()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliEMCALQADataMakerRec& AliEMCALQADataMakerRec::operator = (const AliEMCALQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliEMCALQADataMakerRec();
  new(this) AliEMCALQADataMakerRec(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kEMCAL, task, list) ;  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F * h1 = new TH1F("hESDCaloClusterE",  "ESDs CaloCluster energy in EMCAL",    200, 0., 20.) ; 
  h1->Sumw2() ;
  Add2ESDsList(h1, kESDCaloClusE, !expert, image)  ;                                                     

  TH1I * h2 = new TH1I("hESDCaloClusterM", "ESDs CaloCluster multiplicity in EMCAL", 100, 0,  100) ; 
  h2->Sumw2() ;
  Add2ESDsList(h2, kESDCaloClusM, !expert, image)  ;

  TH1F * h3 = new TH1F("hESDCaloCellA",  "ESDs CaloCell amplitude in EMCAL",    500, 0., 250.) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, kESDCaloCellA, !expert, image)  ;  
 
  TH1I * h4 = new TH1I("hESDCaloCellM", "ESDs CaloCell multiplicity in EMCAL", 200, 0,  1000) ; 
  h4->Sumw2() ;
  Add2ESDsList(h4, kESDCaloCellM, !expert, image) ;
	
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F* h0 = new TH1F("hEMCALRpE","EMCAL RecPoint energies",200, 0.,20.); //GeV
  h0->Sumw2();
  Add2RecPointsList(h0,kRecPE, !expert, image);

  TH1I* h1 = new TH1I("hEMCALRpM","EMCAL RecPoint multiplicities",100,0,100);
  h1->Sumw2();
  Add2RecPointsList(h1,kRecPM, !expert, image);

  TH1I* h2 = new TH1I("hEMCALRpDigM","EMCAL RecPoint Digit Multiplicities",20,0,20);
  h2->Sumw2();
  Add2RecPointsList(h2,kRecPDigM, !expert, image);

}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitRaws()
{
  // create Raws histograms in Raws subdir
  //these need more thought
  /*
   const Bool_t expert   = kTRUE ; 
   const Bool_t saveCorr = kTRUE ; 
   const Bool_t image    = kTRUE ; 
  
  TH1I * h0 = new TH1I("hLowEmcalSupermodules",    "Low Gain digits in EMCAL supermodules",       12, 0, 12) ;
  h0->Sumw2() ;
  Add2RawsList(h0, kNsmodLG, !expert, image, !saveCorr) ;
  TH1I * h1 = new TH1I("hHighEmcalSupermodules",   "High Gain Digits in EMCAL supermodules",       12, 0, 12) ;
  h1->Sumw2() ;
  Add2RawsList(h1, kNsmodHG, !expert, image, !saveCorr) ;

  TH1F * h2 = new TH1F("hLowEmcalRawtime", "Low Gain Time of raw digits in EMCAL", 500, -50., 200.) ;
  h2->Sumw2() ;
  Add2RawsList(h2, kLGtime, !expert, image, !saveCorr) ;
  TH1F * h3 = new TH1F("hHighEmcalRawtime", "High Gain Time of raw digits in EMCAL", 500, -50., 200.) ;
  h3->Sumw2() ;
  Add2RawsList(h3, kHGtime, !expert, image, !saveCorr) ;

  TH1F * h4 = new TH1F("hLowEmcalRawEnergy", "Low Gain Energy of raw digits in EMCAL", 500, 0., 1000.) ;
  h4->Sumw2() ;
  Add2RawsList(h4, kSpecLG, !expert, image, !saveCorr) ;
  TH1F * h5 = new TH1F("hHighEmcalRawEnergy", "High Gain Energy of raw digits in EMCAL",500,0., 1000.) ;
  h5->Sumw2() ;
  Add2RawsList(h5, kSpecHG, !expert, image, !saveCorr) ;

  TH1I * h6 = new TH1I("hLowNtot", "Low Gain Total Number of raw digits in EMCAL", 500, 0, 10000) ;
  h6->Sumw2() ;
  Add2RawsList(h6, kNtotLG, !expert, image, !saveCorr) ;
  TH1I * h7 = new TH1I("hHighNtot", "High Gain Total Number of raw digits in EMCAL",500,0, 10000) ;
  h7->Sumw2() ;
  Add2RawsList(h7, kNtotHG, !expert, image, !saveCorr) ;

  TH1F * h8 = new TH1F("hLowEtot", "Low Gain Total Energy of raw digits in EMCAL", 500, 0., 5000.) ;
  h8->Sumw2() ;
  Add2RawsList(h8, kEtotLG, !expert, image, !saveCorr) ;
  TH1F * h9 = new TH1F("hHighEtot", "High Gain Total Energy of raw digits in EMCAL",500,0., 100000.) ;
  h9->Sumw2() ;
  Add2RawsList(h9, kEtotHG, !expert, image, !saveCorr) ;
  */
  
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs

  Int_t nTot = 0 ; 
  for ( Int_t index = 0; index < esd->GetNumberOfCaloClusters() ; index++ ) {
    AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
    if( clu->IsEMCAL() ) {
      GetESDsData(kESDCaloClusE)->Fill(clu->E()) ;
      nTot++ ;
    } 
  }
  GetESDsData(kESDCaloClusM)->Fill(nTot) ;

  //fill calo cells
  AliESDCaloCells* cells = esd->GetEMCALCells();
  GetESDsData(kESDCaloCellM)->Fill(cells->GetNumberOfCells()) ;

  for ( Int_t index = 0; index < cells->GetNumberOfCells() ; index++ ) {
    GetESDsData(kESDCaloCellA)->Fill(cells->GetAmplitude(index)) ;
  }

}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //Fill prepared histograms with Raw digit properties

  //Raw histogram filling not yet implemented
  //
  //Need to figure out how to get the info we want without having to
  //actually run Raw2Digits twice.
  //I suspect what we actually want is a raw digits method, not a true
  //emcal raw data method, but this doesn't seem to be allowed in
  //AliQADataMakerRec.h

	rawReader->Reset() ;
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  // makes data from RecPoints
  TBranch *emcbranch = clustersTree->GetBranch("EMCALECARP");
  if (!emcbranch) { 
    AliError("can't get the branch with the EMCAL clusters !");
    return;
  }
  TObjArray * emcrecpoints = new TObjArray(100) ;
  emcbranch->SetAddress(&emcrecpoints);
  emcbranch->GetEntry(0);
  
  GetRecPointsData(kRecPM)->Fill(emcrecpoints->GetEntriesFast()) ; 
  TIter next(emcrecpoints) ; 
  AliEMCALRecPoint * rp ; 
  while ( (rp = dynamic_cast<AliEMCALRecPoint *>(next())) ) {
    GetRecPointsData(kRecPE)->Fill( rp->GetEnergy()) ;
    GetRecPointsData(kRecPDigM)->Fill(rp->GetMultiplicity());
  }
  emcrecpoints->Delete();
  delete emcrecpoints;
  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
