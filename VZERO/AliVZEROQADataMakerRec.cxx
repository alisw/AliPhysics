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
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliVZEROQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliVZERORawStream.h"
#include "AliVZEROReconstructor.h"


ClassImp(AliVZEROQADataMakerRec)
           
//____________________________________________________________________________ 
  AliVZEROQADataMakerRec::AliVZEROQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kVZERO), "VZERO Quality Assurance Data Maker")
{
  // constructor
}

//____________________________________________________________________________ 
AliVZEROQADataMakerRec::AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) :
  AliQADataMakerRec()
{
  //copy constructor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliVZEROQADataMakerRec& AliVZEROQADataMakerRec::operator = (const AliVZEROQADataMakerRec& qadm )
{
  // Equal operator
  this->~AliVZEROQADataMakerRec();
  new(this) AliVZEROQADataMakerRec(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kVZERO, task, list) ;  
}

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
 
  TH1I * h1 = new TH1I("hVZERONbPMA", "Number of PMs fired in V0A", 100, 0, 99); 
  h1->Sumw2() ;
  Add2ESDsList(h1, 0)  ;                                                                                                        
  TH1I * h2 = new TH1I("hVZERONbPMC", "Number of PMs fired in V0C", 100, 0, 99); 
  h2->Sumw2() ;
  Add2ESDsList(h2, 1) ;
 
  TH1I * h3 = new TH1I("hVZEROMultA", "Multiplicity in V0A", 50, 0, 49) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, 2) ;
 
  TH1I * h4 = new TH1I("hVZEROMultC", "Multiplicity in V0A", 50, 0, 49) ; 
  h4->Sumw2() ;
  Add2ESDsList(h4, 3) ;
	
}


//____________________________________________________________________________ 
// void AliVZEROQADataMakerRec::InitRaws()
// {
//   // create Raws histograms in Raws subdir
//   TH2I * h0 = new TH2I("hHighVZEROxyMod1","High Gain Rows x Columns for VZERO module 1", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h0,kHGmod1) ;
//   TH2I * h1 = new TH2I("hHighVZEROxyMod2","High Gain Rows x Columns for VZERO module 2", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h1,kHGmod2) ;
//   TH2I * h2 = new TH2I("hHighVZEROxyMod3","High Gain Rows x Columns for VZERO module 3", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h2,kHGmod3) ;
//   TH2I * h3 = new TH2I("hHighVZEROxyMod4","High Gain Rows x Columns for VZERO module 4", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h3,kHGmod4) ;
//   TH2I * h4 = new TH2I("hHighVZEROxyMod5","High Gain Rows x Columns for VZERO module 5", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h4,kHGmod5) ;
//   TH2I * h5 = new TH2I("hLowVZEROxyMod1","Low Gain Rows x Columns for VZERO module 1", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h5,kLGmod1) ;
//   TH2I * h6 = new TH2I("hLowVZEROxyMod2","Low Gain Rows x Columns for VZERO module 2", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h6,kLGmod2) ;
//   TH2I * h7 = new TH2I("hLowVZEROxyMod3","Low Gain Rows x Columns for VZERO module 3", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h7,kLGmod3) ;
//   TH2I * h8 = new TH2I("hLowVZEROxyMod4","Low Gain Rows x Columns for VZERO module 4", 64, 0, 64, 56, 0, 56) ;
//   Add2RawsList(h8,kLGmod4) ;                                                                                                               
//   TH2I * h9 = new TH2I("hLowVZEROxyMod5","Low Gain Rows x Columns for VZERO module 5", 64, 0, 64, 56, 0, 56) ;                               
//   Add2RawsList(h9,kLGmod5) ;                                                                                                               
//                                                                                                                                            
//                                                                                                                                            
//   TH1I * h10 = new TH1I("hLowVZEROModules",    "Low Gain Hits in EMCA VZERO modules",       6, 0, 6) ;                                       
//   h10->Sumw2() ;                                                                                                                           
//   Add2RawsList(h10, kNmodLG) ;                                                                                                             
//   TH1I * h11 = new TH1I("hHighVZEROModules",   "High Gain Hits in EMCA VZERO modules",       6, 0, 6) ;                                      
//   h11->Sumw2() ;                                                                                                                           
//   Add2RawsList(h11, kNmodHG) ;                                                                                                             
//                                                                                                                                            
//   TH1F * h12 = new TH1F("hLowVZERORawtime", "Low Gain Time of raw hits in VZERO", 500, -50., 200.) ;                                            
//   h12->Sumw2() ;                                                                                                                           
//   Add2RawsList(h12, kLGtime) ;                                                                                                             
//   TH1F * h13 = new TH1F("hHighVZERORawtime", "High Gain Time of raw hits in VZERO", 500, -50., 200.) ;                                          
//   h13->Sumw2() ;                                                                                                                           
//   Add2RawsList(h13, kHGtime) ;                                                                                                             
//                                                                                                                                            
//   TH1F * h14 = new TH1F("hLowVZERORawEnergy", "Low Gain Energy of raw hits in VZERO", 500, 0., 1000.) ;                                      
//   h14->Sumw2() ;                                                                                                                           
//   Add2RawsList(h14, kSpecLG) ;                                                                                                             
//   TH1F * h15 = new TH1F("hHighVZERORawEnergy", "High Gain Energy of raw hits in VZERO",500,0., 1000.) ;                                      
//   h15->Sumw2() ;                                                                                                                           
//   Add2RawsList(h15, kSpecHG) ;                                                                                                             
//                                                                                                                                            
//   TH1F * h16 = new TH1F("hLowNtot", "Low Gain Total Number of raw hits in VZERO", 500, 0., 5000.) ;                                         
//   h16->Sumw2() ;                                                                                                                           
//   Add2RawsList(h16, kNtotLG) ;                                                                                                             
//   TH1F * h17 = new TH1F("hHighNtot", "High Gain Total Number of raw hits in VZERO",500,0., 5000.) ;                                         
//   h17->Sumw2() ;                                                                                                                           
//   Add2RawsList(h17, kNtotHG) ;                                                                                                             
//                                                                                                                                            
//   TH1F * h18 = new TH1F("hLowEtot", "Low Gain Total Energy of raw hits in VZERO", 500, 0., 5000.) ;                                       
//   h18->Sumw2() ;                                                                                                                           
//   Add2RawsList(h18, kEtotLG) ;                                                                                                             
//   TH1F * h19 = new TH1F("hHighEtot", "High Gain Total Energy of raw hits in VZERO",500,0., 100000.) ;                                       
//   h19->Sumw2() ;                                                                                                                           
//   Add2RawsList(h19, kEtotHG) ;                                                                                                             
//   
// }

//____________________________________________________________________________
void AliVZEROQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs


  AliESDVZERO *esdVZERO=esd->GetVZEROData();
   
  if (esdVZERO) { 
    GetESDsData(0)->Fill(esdVZERO->GetNbPMV0A());
    GetESDsData(1)->Fill(esdVZERO->GetNbPMV0C());  
    GetESDsData(2)->Fill(esdVZERO->GetMTotV0A());
    GetESDsData(3)->Fill(esdVZERO->GetMTotV0C());  
  }

}

//____________________________________________________________________________
// void AliVZEROQADataMakerRec::MakeRaws(AliRawReader* rawReader)
// {
//   //Fill prepared histograms with Raw digit properties
//   rawReader->Reset() ;
//   AliVZERORawDecoder * decoder ;
//   if(strcmp(AliVZEROReconstructor::GetRecoParamEmc()->DecoderVersion(),"v1")==0)
//     decoder=new AliVZERORawDecoderv1(rawReader);
//   else
//     if(strcmp(AliVZEROReconstructor::GetRecoParamEmc()->DecoderVersion(),"v2")==0)
//       decoder=new AliVZERORawDecoderv2(rawReader);
//     else
//       decoder=new AliVZERORawDecoder(rawReader);
//   decoder->SetOldRCUFormat  (AliVZEROReconstructor::GetRecoParamEmc()->IsOldRCUFormat());
//   decoder->SubtractPedestals(AliVZEROReconstructor::GetRecoParamEmc()->SubtractPedestals());
//   Double_t lgEtot=0. ;
//   Double_t hgEtot=0. ;
//   Int_t lgNtot=0 ;
//   Int_t hgNtot=0 ;
// 
//   while (decoder->NextDigit()) {
//    Int_t module  = decoder->GetModule() ;
//    Int_t row     = decoder->GetRow() ;
//    Int_t col     = decoder->GetColumn() ;
//    Double_t time = decoder->GetTime() ;
//    Double_t energy  = decoder->GetEnergy() ;
//    Bool_t lowGain = decoder->IsLowGain();
//    if (lowGain) {
//      if(energy<2.)
//        continue ;
//      switch(module){
//         case 1: GetRawsData(kLGmod1)->Fill(row-0.5,col-0.5) ; break ;
//         case 2: GetRawsData(kLGmod2)->Fill(row-0.5,col-0.5) ; break ;
//         case 3: GetRawsData(kLGmod3)->Fill(row-0.5,col-0.5) ; break ;
//         case 4: GetRawsData(kLGmod4)->Fill(row-0.5,col-0.5) ; break ;
//         case 5: GetRawsData(kLGmod5)->Fill(row-0.5,col-0.5) ; break ;
//      }                                  
//      GetRawsData(kNmodLG)->Fill(module) ;
//      GetRawsData(kLGtime)->Fill(time) ; 
//      GetRawsData(kSpecLG)->Fill(energy) ;    
//      lgEtot+=energy ;
//      lgNtot++ ;   
//    } else {        
//      if(energy<8.)
//        continue ;
//      switch (module){
//          case 1: GetRawsData(kHGmod1)->Fill(row-0.5,col-0.5) ; break ;
//          case 2: GetRawsData(kHGmod2)->Fill(row-0.5,col-0.5) ; break ;
//          case 3: GetRawsData(kHGmod3)->Fill(row-0.5,col-0.5) ; break ;
//          case 4: GetRawsData(kHGmod4)->Fill(row-0.5,col-0.5) ; break ;
//          case 5: GetRawsData(kHGmod5)->Fill(row-0.5,col-0.5) ; break ;
//      }              
//      GetRawsData(kNmodHG)->Fill(module) ; 
//      GetRawsData(kHGtime)->Fill(time) ;  
//      GetRawsData(kSpecHG)->Fill(energy) ;
//      hgEtot+=energy ; 
//      hgNtot++ ;  
//    }                 
//   }                    
//   GetRawsData(kEtotLG)->Fill(lgEtot) ; 
//   GetRawsData(kEtotHG)->Fill(hgEtot) ;  
//   GetRawsData(kNtotLG)->Fill(lgNtot) ;
//   GetRawsData(kNtotHG)->Fill(hgNtot) ;
// }

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::StartOfDetectorCycle()
{
  // Detector specific actions at start of cycle
  
}
