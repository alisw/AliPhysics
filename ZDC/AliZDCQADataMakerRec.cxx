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

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Produces the data needed to calculate the ZDC quality assurance. //
//  QA objects are 1 & 2 Dimensional histograms.                     //
//  author: C. Oppedisano                                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <TProfile.h>
#include <Riostream.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliZDCQADataMakerRec.h"
#include "AliZDCRawStream.h"
#include "AliESDZDC.h"
#include "AliESDEvent.h"

ClassImp(AliZDCQADataMakerRec)
           
//____________________________________________________________________________ 
  AliZDCQADataMakerRec::AliZDCQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kZDC), "ZDC Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliZDCQADataMakerRec::AliZDCQADataMakerRec(const AliZDCQADataMakerRec& qadm) :
  AliQADataMakerRec() 
{
  //copy ctor 
  SetName((const char*)qadm.GetName()); 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliZDCQADataMakerRec& AliZDCQADataMakerRec::operator = (const AliZDCQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliZDCQADataMakerRec();
  new(this) AliZDCQADataMakerRec(qadm);
  return *this;
}

//____________________________________________________________________________

void AliZDCQADataMakerRec::InitRaws()
{
  // create Digits histograms in Digits subdir
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hRawZNCTot = new TH1F("hRawZNCTot", "Raw signal in ZNC", 100, 0., 6000.);
  TH1F * hRawZNATot = new TH1F("hRawZNATot", "Raw signal in ZNA", 100, 0., 6000.);
  TH1F * hRawZPCTot = new TH1F("hRawZPCTot", "Raw signal in ZPC", 100, 0., 10000.);
  TH1F * hRawZPATot = new TH1F("hRawZPATot", "Raw signal in ZPA", 100, 0., 10000.);
  Add2RawsList(hRawZNCTot, 0);
  Add2RawsList(hRawZNATot, 1);
  Add2RawsList(hRawZPCTot, 2);
  Add2RawsList(hRawZPATot, 3);
  //
  TH1F * hRawSumQZNC = new TH1F("hRawSumQZNC", "Raw summed 4 ZNC quadrants",100, 0., 4000.);
  TH1F * hRawSumQZNA = new TH1F("hRawSumQZNA", "Raw summed 4 ZNA quadrants",100, 0., 4000.);
  TH1F * hRawSumQZPC = new TH1F("hRawSumQZPC", "Raw summed 4 ZPC quadrants",100, 0., 4000.);
  TH1F * hRawSumQZPA = new TH1F("hRawSumQZPA", "Raw summed 4 ZPA quadrants",100, 0., 4000.);
  Add2RawsList(hRawSumQZNC, 4, kTRUE);
  Add2RawsList(hRawSumQZNA, 5, kTRUE);
  Add2RawsList(hRawSumQZPC, 6, kTRUE);
  Add2RawsList(hRawSumQZPA, 7, kTRUE);
  //
  TH1F * hRawPMCZNC = new TH1F("hRawPMCZNC", "Raw common ZNC PMT",100, 0., 4000.);
  TH1F * hRawPMCZNA = new TH1F("hRawPMCZNA", "Raw common ZNA PMT",100, 0., 4000.);
  TH1F * hRawPMCZPC = new TH1F("hRawPMCZPC", "Raw common ZPC PMT",100, 0., 4000.);
  TH1F * hRawPMCZPA = new TH1F("hRawPMCZPA", "Raw common ZPA PMT",100, 0., 4000.);
  Add2RawsList(hRawPMCZNC, 8 , kTRUE);
  Add2RawsList(hRawPMCZNA, 9 , kTRUE);
  Add2RawsList(hRawPMCZPC, 10, kTRUE);
  Add2RawsList(hRawPMCZPA, 11, kTRUE);
  // 
/*  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hRawZNCTotlg = new TH1F("hRawZNCTotlg", "Rawit lg signal in ZNC", 100, 0., 6000.);
  TH1F * hRawZNATotlg = new TH1F("hRawZNATotlg", "Rawit lg signal in ZNA", 100, 0., 6000.);
  TH1F * hRawZPCTotlg = new TH1F("hRawZPCTotlg", "Rawit lg signal in ZPC", 100, 0., 10000.);
  TH1F * hRawZPATotlg = new TH1F("hRawZPATotlg", "Rawit lg signal in ZPA", 100, 0., 10000.);
  Add2RawsList(hRawZNCTotlg, 12);
  Add2RawsList(hRawZNATotlg, 13);
  Add2RawsList(hRawZPCTotlg, 14);
  Add2RawsList(hRawZPATotlg, 15);
  //
  TH1F * hRawSumQZNClg = new TH1F("hRawSumQZNClg", "Raw summed 4 lg ZNC quadrants",100, 0., 4000.);
  TH1F * hRawSumQZNAlg = new TH1F("hRawSumQZNAlg", "Raw summed 4 lg ZNA quadrants",100, 0., 4000.);
  TH1F * hRawSumQZPClg = new TH1F("hRawSumQZPClg", "Raw summed 4 lg ZPC quadrants",100, 0., 4000.);
  TH1F * hRawSumQZPAlg = new TH1F("hRawSumQZPAlg", "Raw summed 4 lg ZPA quadrants",100, 0., 4000.);
  Add2RawsList(hRawSumQZNClg, 16, kTRUE);
  Add2RawsList(hRawSumQZNAlg, 17, kTRUE);
  Add2RawsList(hRawSumQZPClg, 18, kTRUE);
  Add2RawsList(hRawSumQZPAlg, 19, kTRUE);
  //
  TH1F * hRawPMCZNClg = new TH1F("hRawPMCZNClg", "Raw common lg ZNC PMT",100, 0., 4000.);
  TH1F * hRawPMCZNAlg = new TH1F("hRawPMCZNAlg", "Raw common lg ZNA PMT",100, 0., 4000.);
  TH1F * hRawPMCZPClg = new TH1F("hRawPMCZPClg", "Raw common lg ZPC PMT",100, 0., 4000.);
  TH1F * hRawPMCZPAlg = new TH1F("hRawPMCZPAlg", "Raw common lg ZPA PMT",100, 0., 4000.);
  Add2RawsList(hRawPMCZNClg, 20, kTRUE);
  Add2RawsList(hRawPMCZNAlg, 21, kTRUE);
  Add2RawsList(hRawPMCZPClg, 22, kTRUE);
  Add2RawsList(hRawPMCZPAlg, 23, kTRUE);*/
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::InitESDs()
{
  //Booking ESDs histograms
  //
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH2F * hZNC  = new TH2F("hZNC", "Centroid in ZNC", 100, -5.,5.,100,-5.,5.);
  TH2F * hZNA  = new TH2F("hZNA", "Centroid in ZNA", 100, -5.,5.,100,-5.,5.);
  Add2ESDsList(hZNC, 0, !expert, image);
  Add2ESDsList(hZNA, 1, !expert, image);
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hESDZNCTot = new TH1F("hESDZNCTot", "Energy in ZNC", 100, 0., 6000.);
  TH1F * hESDZNATot = new TH1F("hESDZNATot", "Energy in ZNA", 100, 0., 6000.);
  TH1F * hESDZPCTot = new TH1F("hESDZPCTot", "Energy in ZPC", 100, 0., 10000.);
  TH1F * hESDZPATot = new TH1F("hESDZPATot", "Energy in ZPA", 100, 0., 10000.);
  Add2ESDsList(hESDZNCTot, 2, !expert, image);
  Add2ESDsList(hESDZNATot, 3, !expert, image);
  Add2ESDsList(hESDZPCTot, 4, !expert, image);
  Add2ESDsList(hESDZPATot, 5, !expert, image);
  //
  TH1F * hESDSumQZNC = new TH1F("hESDSumQZNC", "Sum of 4 ZNC energy",100, 0., 4000.);
  TH1F * hESDSumQZNA = new TH1F("hESDSumQZNA", "Sum of 4 ZNA energy",100, 0., 4000.);
  TH1F * hESDSumQZPC = new TH1F("hESDSumQZPC", "Sum of 4 ZPC energy",100, 0., 4000.);
  TH1F * hESDSumQZPA = new TH1F("hESDSumQZPA", "Sum of 4 ZPA energy",100, 0., 4000.);
  Add2ESDsList(hESDSumQZNC, 6, expert, !image);
  Add2ESDsList(hESDSumQZNA, 7, expert, !image);
  Add2ESDsList(hESDSumQZPC, 8, expert, !image);
  Add2ESDsList(hESDSumQZPA, 9, expert, !image);
  //
  TH1F * hESDPMCZNC = new TH1F("hESDPMCZNC", "Energy in common ZNC PMT",100, 0., 4000.);
  TH1F * hESDPMCZNA = new TH1F("hESDPMCZNA", "Energy in common ZNA PMT",100, 0., 4000.);
  TH1F * hESDPMCZPC = new TH1F("hESDPMCZPC", "Energy in common ZPC PMT",100, 0., 4000.);
  TH1F * hESDPMCZPA = new TH1F("hESDPMCZPA", "Energy in common ZPA PMT",100, 0., 4000.);
  Add2ESDsList(hESDPMCZNC, 10, expert, !image);
  Add2ESDsList(hESDPMCZNA, 11, expert, !image);
  Add2ESDsList(hESDPMCZPC, 12, expert, !image);
  Add2ESDsList(hESDPMCZPA, 13, expert, !image);
  // 
/*  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hESDZNCTotlg = new TH1F("hESDZNCTotlg", "ESD lg signal in ZNC", 100, 0., 6000.);
  TH1F * hESDZNATotlg = new TH1F("hESDZNATotlg", "ESD lg signal in ZNA", 100, 0., 6000.);
  TH1F * hESDZPCTotlg = new TH1F("hESDZPCTotlg", "ESD lg signal in ZPC", 100, 0., 10000.);
  TH1F * hESDZPATotlg = new TH1F("hESDZPATotlg", "ESD lg signal in ZPA", 100, 0., 10000.);
  Add2ESDsList(hESDZNCTotlg, !expert, image);
  Add2ESDsList(hESDZNATotlg, !expert, image);
  Add2ESDsList(hESDZPCTotlg, !expert, image);
  Add2ESDsList(hESDZPATotlg, !expert, image);
  //
  TH1F * hESDSumQZNClg = new TH1F("hESDSumQZNClg", "Sum of 4 lg ZNC sectors",100, 0., 4000.);
  TH1F * hESDSumQZNAlg = new TH1F("hESDSumQZNAlg", "Sum of 4 lg ZNA sectors",100, 0., 4000.);
  TH1F * hESDSumQZPClg = new TH1F("hESDSumQZPClg", "Sum of 4 lg ZPC sectors",100, 0., 4000.);
  TH1F * hESDSumQZPAlg = new TH1F("hESDSumQZPAlg", "Sum of 4 lg ZPA sectors",100, 0., 4000.);
  Add2ESDsList(hESDSumQZNClg, 18, expert, !image);
  Add2ESDsList(hESDSumQZNAlg, 19, expert, !image);
  Add2ESDsList(hESDSumQZPClg, 20, expert, !image);
  Add2ESDsList(hESDSumQZPAlg, 21, expert, !image);
  //
  TH1F * hESDPMCZNClg = new TH1F("hESDPMCZNClg", "Signal in common ZNC lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZNAlg = new TH1F("hESDPMCZNAlg", "Signal in common ZNA lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZPClg = new TH1F("hESDPMCZPClg", "Signal in common ZPC lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZPAlg = new TH1F("hESDPMCZPAlg", "Signal in common ZPA lg PMT",100, 0., 4000.);
  Add2ESDsList(hESDPMCZNClg, 22, expert, !image);
  Add2ESDsList(hESDPMCZNAlg, 23, expert, !image);
  Add2ESDsList(hESDPMCZPClg, 24, expert, !image);
  Add2ESDsList(hESDPMCZPAlg, 25, expert, !image);*/
}
  
//____________________________________________________________________________

void AliZDCQADataMakerRec::MakeRaws(AliRawReader *rawReader)
{
  // Filling Raws QA histos
  //
	rawReader->Reset() ; 
  Float_t sum_ZNC=0., sum_ZNA=0., sum_ZPC=0., sum_ZPA=0.;
  Float_t sumQ_ZNC=0., sumQ_ZNA=0., sumQ_ZPC=0., sumQ_ZPA=0.;
  //Float_t sum_ZNC_lg=0., sum_ZNA_lg=0., sum_ZPC_lg=0., sum_ZPA_lg=0.;
  //Float_t sumQ_ZNC_lg=0., sumQ_ZNA_lg=0., sumQ_ZPC_lg=0., sumQ_ZPA_lg=0.;
  //
  AliZDCRawStream stream(rawReader);
  while(stream.Next()){
    if(stream.IsADCDataWord() && 
     (stream.GetADCModule()==0 || stream.GetADCModule()==1)){
       if(stream.GetSector(0)==1){
         if(stream.GetADCGain()==0){
	   sum_ZNC += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZNC += stream.GetADCValue();
	   else GetRawsData(8)->Fill(stream.GetADCValue());
	 }
	 /*else{
	   sum_ZNC_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZNC_lg += stream.GetADCValue();
	   else GetRawsData(20)->Fill(stream.GetADCValue());
	 }*/
       }
       else if(stream.GetSector(0)==2){
         if(stream.GetADCGain()==0){
	   sum_ZPC += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZPC += stream.GetADCValue();
	   else GetRawsData(10)->Fill(stream.GetADCValue());
	 }
	 /*else{
	   sum_ZPC_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZPC_lg += stream.GetADCValue();
	   else GetRawsData(22)->Fill(stream.GetADCValue());
	 }*/
       }
       else if(stream.GetSector(0)==4){
         if(stream.GetADCGain()==0){
	   sum_ZNA += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZNA += stream.GetADCValue();
	   else GetRawsData(9)->Fill(stream.GetADCValue());
	 }
	 /*else{
	   sum_ZNA_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZNA_lg += stream.GetADCValue();
	   else GetRawsData(21)->Fill(stream.GetADCValue());
	 }*/
       }
       else if(stream.GetSector(0)==5){
         if(stream.GetADCGain()==0){
	   sum_ZPA += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZPA += stream.GetADCValue();
	   else GetRawsData(11)->Fill(stream.GetADCValue());
	 }
	 /*else{
	   sum_ZPA_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZPA_lg += stream.GetADCValue();
	   else GetRawsData(23)->Fill(stream.GetADCValue());
	 }*/
       }
    }
  }
  //
  GetRawsData(0)->Fill(sum_ZNC);
  GetRawsData(1)->Fill(sum_ZNA);
  GetRawsData(2)->Fill(sum_ZPC);
  GetRawsData(3)->Fill(sum_ZPA);
  //
  GetRawsData(4)->Fill(sumQ_ZNC);
  GetRawsData(5)->Fill(sumQ_ZNA);
  GetRawsData(6)->Fill(sumQ_ZPC);
  GetRawsData(7)->Fill(sumQ_ZPA);
  //
  /*GetRawsData(12)->Fill(sum_ZNC_lg);
  GetRawsData(13)->Fill(sum_ZNA_lg);
  GetRawsData(14)->Fill(sum_ZPC_lg);
  GetRawsData(15)->Fill(sum_ZPA_lg);
  //
  GetRawsData(16)->Fill(sumQ_ZNC_lg);
  GetRawsData(17)->Fill(sumQ_ZNA_lg);
  GetRawsData(18)->Fill(sumQ_ZPC_lg);
  GetRawsData(19)->Fill(sumQ_ZPA_lg);*/
  //
  stream.Delete();
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  //
  AliESDZDC * zdcESD =  esd->GetESDZDC();
  //
  Double32_t * centr_ZNC = zdcESD->GetZNCCentroid();
  GetESDsData(0)->Fill(centr_ZNC[0], centr_ZNC[1]);

  Double32_t * centr_ZNA = zdcESD->GetZNACentroid();
  GetESDsData(1)->Fill(centr_ZNA[0], centr_ZNA[1]);

  //
  GetESDsData(2)->Fill(esd->GetZDCN1Energy());
  GetESDsData(3)->Fill(esd->GetZDCN2Energy());
  GetESDsData(4)->Fill(esd->GetZDCP1Energy());
  GetESDsData(5)->Fill(esd->GetZDCP2Energy());
  //
  Double_t sumQZNC=0., sumQZPC=0., sumQZNA=0., sumQZPA=0.;
  //Double_t sumQZNC_lg=0., sumQZPC_lg=0., sumQZNA_lg=0., sumQZPA_lg=0.;
  //
  const Double_t *towZNC, *towZPC, *towZNA, *towZPA;
  //const Double_t *towZNC_lg, *towZPC_lg, *towZNA_lg, *towZPA_lg;
  //
  towZNC = zdcESD->GetZN1TowerEnergy();
  towZPC = zdcESD->GetZP1TowerEnergy();
  towZNA = zdcESD->GetZN2TowerEnergy();
  towZPA = zdcESD->GetZP2TowerEnergy();
  //
  /*towZNC_lg = zdcESD->GetZN1TowerEnergyLR();
  towZPC_lg = zdcESD->GetZP1TowerEnergyLR();
  towZNA_lg = zdcESD->GetZN2TowerEnergyLR();
  towZPA_lg = zdcESD->GetZP2TowerEnergyLR();*/
  //
  for(Int_t i=0; i<5; i++){
     if(i==0){
       GetESDsData(10)->Fill(towZNC[i]);
       GetESDsData(11)->Fill(towZNA[i]);
       GetESDsData(12)->Fill(towZPC[i]);
       GetESDsData(13)->Fill(towZPA[i]);
       //
       /*GetESDsData(22)->Fill(towZNC_lg[i]);
       GetESDsData(23)->Fill(towZNA_lg[i]);
       GetESDsData(24)->Fill(towZPC_lg[i]);
       GetESDsData(25)->Fill(towZPA_lg[i]);*/
     }
     else{
       sumQZNC += towZNC[i];
       sumQZPC += towZPC[i];
       sumQZNA += towZNA[i];
       sumQZPA += towZPA[i];
       //
       /*sumQZNC_lg += towZNC_lg[i];
       sumQZPC_lg += towZPC_lg[i];
       sumQZNA_lg += towZNA_lg[i];
       sumQZPA_lg += towZPA_lg[i];*/
     }
  }
  GetESDsData(6)->Fill(sumQZNC);
  GetESDsData(7)->Fill(sumQZNA);
  GetESDsData(8)->Fill(sumQZPC);
  GetESDsData(9)->Fill(sumQZPA);
  //
  /*GetESDsData(18)->Fill(sumQZNC_lg);
  GetESDsData(19)->Fill(sumQZNA_lg);
  GetESDsData(20)->Fill(sumQZPC_lg);
  GetESDsData(21)->Fill(sumQZPA_lg);*/
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliZDCQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kZDC, task, list) ;  
}

