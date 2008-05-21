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
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliEMCALQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliEMCALRecPoint.h" 
#include "AliEMCALRawUtils.h"
#include "AliEMCALReconstructor.h"
#include "AliEMCALRecParam.h"

ClassImp(AliEMCALQADataMakerRec)
           
//____________________________________________________________________________ 
  AliEMCALQADataMakerRec::AliEMCALQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kEMCAL), "EMCAL Quality Assurance Data Maker")
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
void AliEMCALQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kEMCAL, task, list) ;  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
 
  TH1F * h1 = new TH1F("hESDEmcalSpectrum",  "ESDs spectrum in EMCAL",    200, 0., 20.) ; 
  h1->Sumw2() ;
  Add2ESDsList(h1, kESDSpec)  ;                                                                                                        
  TH1I * h2 = new TH1I("hESDEmcalMul", "ESDs multiplicity distribution in EMCAL", 100, 0,  100) ; 
  h2->Sumw2() ;
  Add2ESDsList(h2, kESDNtot) ;
 
  TH1I * h3 = new TH1I("hESDEmcalEtot", "ESDs Etot in EMCAL", 100, 0,  1000.) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, kESDEtot) ;
 
  TH1F * h4 = new TH1F("hESDEmcalPid",    "ESDs PID distribution in EMCAL",       100, 0., 1.) ;
  h4->Sumw2() ;
  Add2ESDsList(h4, kESDpid) ;
	
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  TH2I * h0 = new TH2I("hRpEMCALxySMod1","RecPoints Rows x Columns for EMCAL supermodule 1", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h0,kRPsmod1) ;
  TH2I * h1 = new TH2I("hRpEMCALxySMod2","RecPoints Rows x Columns for EMCAL supermodule 2", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h1,kRPsmod2) ;
  TH2I * h2 = new TH2I("hRpEMCALxySMod3","RecPoints Rows x Columns for EMCAL supermodule 3", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h2,kRPsmod3) ;
  TH2I * h3 = new TH2I("hRpEMCALxySMod4","RecPoints Rows x Columns for EMCAL supermodule 4", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h3,kRPsmod4) ;
  TH2I * h4 = new TH2I("hRpEMCALxySMod5","RecPoints Rows x Columns for EMCAL supermodule 5", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h4,kRPsmod5) ;
  TH2I * h5 = new TH2I("hRpEMCALxySMod6","RecPoints Rows x Columns for EMCAL supermodule 6", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h5,kRPsmod6) ;
  TH2I * h6 = new TH2I("hRpEMCALxySMod7","RecPoints Rows x Columns for EMCAL supermodule 7", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h6,kRPsmod7) ;
  TH2I * h7 = new TH2I("hRpEMCALxySMod8","RecPoints Rows x Columns for EMCAL supermodule 8", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h7,kRPsmod8) ;
  TH2I * h8 = new TH2I("hRpEMCALxySMod9","RecPoints Rows x Columns for EMCAL supermodule 9", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h8,kRPsmod9) ;
  TH2I * h9 = new TH2I("hRpEMCALxySMod10","RecPoints Rows x Columns for EMCAL supermodule 10", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h9,kRPsmod10) ;
  TH2I * h10 = new TH2I("hRpEMCALxySMod11","RecPoints Rows x Columns for EMCAL supermodule 11", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h10,kRPsmod11) ;
  TH2I * h11 = new TH2I("hRpEMCALxySMod12","RecPoints Rows x Columns for EMCAL supermodule 12", 64, -72., 72., 56, -63., 63.) ;
  Add2RecPointsList(h11,kRPsmod12) ;
 
  TH1F * h12 = new TH1F("hEmcalRecPointsSpectrum",  "RecPoints spectrum in EMCAL",   2000, 0., 20.) ; 
  h12->Sumw2() ;
  Add2RecPointsList(h12, kRPSpec)  ;

  TH1I * h13 = new TH1I("hEmcalRecPointsMul", "RecPoints multiplicity distribution in EMCAL", 100, 0,  100) ; 
  h13->Sumw2() ;
  Add2RecPointsList(h13, kRPNtot) ;

  TH1I * h14 = new TH1I("hEmcalRecPointsEtot", "RecPoints Etot in EMCAL", 200, 0,  200.) ; 
  h14->Sumw2() ;
  Add2RecPointsList(h14, kRPEtot) ;

}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::InitRaws()
{
  // create Raws histograms in Raws subdir
  TH2I * h0 = new TH2I("hHighEMCALxySMod1","High Gain Rows x Columns for EMCAL supermodule 1", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h0,kHGsmod1) ;
  TH2I * h1 = new TH2I("hHighEMCALxySMod2","High Gain Rows x Columns for EMCAL supermodule 2", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h1,kHGsmod2) ;
  TH2I * h2 = new TH2I("hHighEMCALxySMod3","High Gain Rows x Columns for EMCAL supermodule 3", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h2,kHGsmod3) ;
  TH2I * h3 = new TH2I("hHighEMCALxySMod4","High Gain Rows x Columns for EMCAL supermodule 4", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h3,kHGsmod4) ;
  TH2I * h4 = new TH2I("hHighEMCALxySMod5","High Gain Rows x Columns for EMCAL supermodule 5", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h4,kHGsmod5) ;
  TH2I * h5 = new TH2I("hHighEMCALxySMod6","High Gain Rows x Columns for EMCAL supermodule 6", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h5,kHGsmod6) ;
  TH2I * h6 = new TH2I("hHighEMCALxySMod7","High Gain Rows x Columns for EMCAL supermodule 7", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h6,kHGsmod7) ;
  TH2I * h7 = new TH2I("hHighEMCALxySMod8","High Gain Rows x Columns for EMCAL supermodule 8", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h7,kHGsmod8) ;
  TH2I * h8 = new TH2I("hHighEMCALxySMod9","High Gain Rows x Columns for EMCAL supermodule 9", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h8,kHGsmod9) ;
  TH2I * h9 = new TH2I("hHighEMCALxySMod10","High Gain Rows x Columns for EMCAL supermodule 10", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h9,kHGsmod10) ;
  TH2I * h10 = new TH2I("hHighEMCALxySMod11","High Gain Rows x Columns for EMCAL supermodule 11", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h10,kHGsmod11) ;
  TH2I * h11 = new TH2I("hHighEMCALxySMod12","High Gain Rows x Columns for EMCAL supermodule 12", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h11,kHGsmod12) ;


  TH2I * h12 = new TH2I("hLowEMCALxySMod1","Low Gain Rows x Columns for EMCAL supermodule 1", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h12,kLGsmod1) ;
  TH2I * h13 = new TH2I("hLowEMCALxySMod2","Low Gain Rows x Columns for EMCAL supermodule 2", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h13,kLGsmod2) ;
  TH2I * h14 = new TH2I("hLowEMCALxySMod3","Low Gain Rows x Columns for EMCAL supermodule 3", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h14,kLGsmod3) ;
  TH2I * h15 = new TH2I("hLowEMCALxySMod4","Low Gain Rows x Columns for EMCAL supermodule 4", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h15,kLGsmod4) ;    
  TH2I * h16 = new TH2I("hLowEMCALxySMod5","Low Gain Rows x Columns for EMCAL supermodule 5", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h16,kLGsmod5) ;
  TH2I * h17 = new TH2I("hLowEMCALxySMod6","Low Gain Rows x Columns for EMCAL supermodule 6", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h17,kLGsmod6) ;
  TH2I * h18 = new TH2I("hLowEMCALxySMod7","Low Gain Rows x Columns for EMCAL supermodule 7", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h18,kLGsmod7) ;
  TH2I * h19 = new TH2I("hLowEMCALxySMod8","Low Gain Rows x Columns for EMCAL supermodule 8", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h19,kLGsmod8) ;
  TH2I * h20 = new TH2I("hLowEMCALxySMod9","Low Gain Rows x Columns for EMCAL supermodule 9", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h20,kLGsmod9) ;    
  TH2I * h21 = new TH2I("hLowEMCALxySMod10","Low Gain Rows x Columns for EMCAL supermodule 10", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h21,kLGsmod10) ;
  TH2I * h22 = new TH2I("hLowEMCALxySMod11","Low Gain Rows x Columns for EMCAL supermodule 11", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h22,kLGsmod11) ;
  TH2I * h23 = new TH2I("hLowEMCALxySMod12","Low Gain Rows x Columns for EMCAL supermodule 12", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h23,kLGsmod12) ;


  TH1I * h24 = new TH1I("hLowEmcalSupermodules",    "Low Gain Hits in EMCAL supermodules",       12, 0, 12) ;
  h24->Sumw2() ;
  Add2RawsList(h24, kNsmodLG) ;
  TH1I * h25 = new TH1I("hHighEmcalSupermodules",   "High Gain Hits in EMCAL supermodules",       12, 0, 12) ;
  h25->Sumw2() ;
  Add2RawsList(h25, kNsmodHG) ;

  TH1F * h26 = new TH1F("hLowEmcalRawtime", "Low Gain Time of raw hits in EMCAL", 500, -50., 200.) ;
  h26->Sumw2() ;
  Add2RawsList(h26, kLGtime) ;
  TH1F * h27 = new TH1F("hHighEmcalRawtime", "High Gain Time of raw hits in EMCAL", 500, -50., 200.) ;
  h27->Sumw2() ;
  Add2RawsList(h27, kHGtime) ;

  TH1F * h28 = new TH1F("hLowEmcalRawEnergy", "Low Gain Energy of raw hits in EMCAL", 500, 0., 1000.) ;
  h28->Sumw2() ;
  Add2RawsList(h28, kSpecLG) ;
  TH1F * h29 = new TH1F("hHighEmcalRawEnergy", "High Gain Energy of raw hits in EMCAL",500,0., 1000.) ;
  h29->Sumw2() ;
  Add2RawsList(h29, kSpecHG) ;

  TH1F * h30 = new TH1F("hLowNtot", "Low Gain Total Number of raw hits in EMCAL", 500, 0., 5000.) ;
  h30->Sumw2() ;
  Add2RawsList(h30, kNtotLG) ;
  TH1F * h31 = new TH1F("hHighNtot", "High Gain Total Number of raw hits in EMCAL",500,0., 5000.) ;
  h31->Sumw2() ;
  Add2RawsList(h31, kNtotHG) ;

  TH1F * h32 = new TH1F("hLowEtot", "Low Gain Total Energy of raw hits in EMCAL", 500, 0., 5000.) ;
  h32->Sumw2() ;
  Add2RawsList(h32, kEtotLG) ;
  TH1F * h33 = new TH1F("hHighEtot", "High Gain Total Energy of raw hits in EMCAL",500,0., 100000.) ;
  h33->Sumw2() ;
  Add2RawsList(h33, kEtotHG) ;
  
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs

  Int_t nTot = 0 ; 
  Double_t eTot = 0 ; 
  for ( Int_t index = 0; index < esd->GetNumberOfCaloClusters() ; index++ ) {
    AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
    if( clu->IsEMCAL() ) {
      GetESDsData(kESDSpec)->Fill(clu->E()) ;
      Double_t *pid=clu->GetPid() ;
      GetESDsData(kESDpid)->Fill(pid[AliPID::kPhoton]) ;
      eTot+=clu->E() ;
      nTot++ ;
    } 
  }
  GetESDsData(kESDNtot)->Fill(nTot) ;
  GetESDsData(kESDEtot)->Fill(eTot) ;
}

//____________________________________________________________________________
void AliEMCALQADataMakerRec::MakeRaws(AliRawReader* /* rawReader */)
{
  //Fill prepared histograms with Raw digit properties

  //Raw histogram filling not yet implemented

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
  
  GetRecPointsData(kRPNtot)->Fill(emcrecpoints->GetEntriesFast()) ; 
  TIter next(emcrecpoints) ; 
  AliEMCALRecPoint * rp ; 
  Double_t eTot = 0. ; 
  while ( (rp = dynamic_cast<AliEMCALRecPoint *>(next())) ) {
    GetRecPointsData(kRPSpec)->Fill( rp->GetEnergy()) ;
    Int_t smod = rp->GetSuperModuleNumber() ;
    TVector3 pos ;
    rp->GetLocalPosition(pos) ;
    switch(smod){
      case 0: GetRecPointsData(kRPsmod1)->Fill(pos.X(),pos.Z()) ; break ;
      case 1: GetRecPointsData(kRPsmod2)->Fill(pos.X(),pos.Z()) ; break ;
      case 2: GetRecPointsData(kRPsmod3)->Fill(pos.X(),pos.Z()) ; break ;
      case 3: GetRecPointsData(kRPsmod4)->Fill(pos.X(),pos.Z()) ; break ;
      case 4: GetRecPointsData(kRPsmod5)->Fill(pos.X(),pos.Z()) ; break ;
      case 5: GetRecPointsData(kRPsmod6)->Fill(pos.X(),pos.Z()) ; break ;
      case 6: GetRecPointsData(kRPsmod7)->Fill(pos.X(),pos.Z()) ; break ;
      case 7: GetRecPointsData(kRPsmod8)->Fill(pos.X(),pos.Z()) ; break ;
      case 8: GetRecPointsData(kRPsmod9)->Fill(pos.X(),pos.Z()) ; break ;
      case 9: GetRecPointsData(kRPsmod10)->Fill(pos.X(),pos.Z()) ; break ;
      case 10: GetRecPointsData(kRPsmod11)->Fill(pos.X(),pos.Z()) ; break ;
      case 11: GetRecPointsData(kRPsmod12)->Fill(pos.X(),pos.Z()) ; break ;
    }
    
    eTot+= rp->GetEnergy() ;
  }
  GetRecPointsData(kRPEtot)->Fill(eTot) ;
  emcrecpoints->Delete();
  delete emcrecpoints;
  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
