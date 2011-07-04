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
#include <TParameter.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliCaloRawStreamV3.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliPHOSQADataMakerRec.h"
#include "AliPHOSDigit.h"
#include "AliPHOSCpvRecPoint.h" 
#include "AliPHOSEmcRecPoint.h" 
#include "AliPHOSRawFitterv0.h"
#include "AliPHOSRawFitterv1.h"
#include "AliPHOSRawFitterv2.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSRecParticle.h" 
#include "AliPHOSTrackSegment.h" 
#include "AliQAChecker.h"
#include "AliRawReader.h"

ClassImp(AliPHOSQADataMakerRec)
           
//____________________________________________________________________________ 
  AliPHOSQADataMakerRec::AliPHOSQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kPHOS), "PHOS Quality Assurance Data Maker")
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
void AliPHOSQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  ResetEventTrigClasses();

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! IsValidEventSpecie(specie, list)) continue;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie)); 
    //
    for (int itc=-1;itc<GetNTrigClasses();itc++) {  // RS: loop over all trigger clones
      //
      if(GetRawsData(kHGqualMod0,itc) && GetRawsData(kHGmod0,itc))
	GetRawsData(kHGqualMod0,itc)->Divide( GetRawsData(kHGmod0,itc) );
      if(GetRawsData(kHGqualMod1,itc) && GetRawsData(kHGmod1,itc))
	GetRawsData(kHGqualMod1,itc)->Divide( GetRawsData(kHGmod1,itc) );
      if(GetRawsData(kHGqualMod2,itc) && GetRawsData(kHGmod2,itc))
	GetRawsData(kHGqualMod2,itc)->Divide( GetRawsData(kHGmod2,itc) );
      if(GetRawsData(kHGqualMod3,itc) && GetRawsData(kHGmod3,itc))
	GetRawsData(kHGqualMod3,itc)->Divide( GetRawsData(kHGmod3,itc) );
      if(GetRawsData(kHGqualMod4,itc) && GetRawsData(kHGmod4,itc))
	GetRawsData(kHGqualMod4,itc)->Divide( GetRawsData(kHGmod4,itc) );
    } // RS: loop over all trigger clones
  }
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kPHOS, task, list) ;  
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
 
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1F * h1 = new TH1F("hESDPhosSpectrum",  "ESDs spectrum in PHOS; Energy [MeV];Counts"                ,  200, 0.,   20.) ;
  h1->Sumw2() ;
  Add2ESDsList(h1, kESDSpec, !expert, image)  ;

  TH1I * h2 = new TH1I("hESDPhosMul",       "ESDs multiplicity distribution in PHOS; # of clusters;Counts", 100, 0,   100 ) ; 
  h2->Sumw2() ;
  Add2ESDsList(h2, kESDNtot, !expert, image) ;
 
  TH1F * h3 = new TH1F("hESDPhosEtot",      "ESDs total energy;Energy [MeV];Counts"                     , 2000, 0,  200.) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, kESDEtot, !expert, image) ;  //Expert histo
 
  TH1F * h4 = new TH1F("hESDpid",           "ESDs PID distribution in PHOS;Particle Id;Counts"         , 100, 0.,    1.) ;
  h4->Sumw2() ;
  Add2ESDsList(h4, kESDpid, !expert, image) ; //Expert histo
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line	
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  TH1I * h0 = new TH1I("hPhosDigits",    "Digits amplitude distribution in PHOS;Amplitude [ADC counts];Counts",    500, 0, 1000) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, kDigits, !expert, image) ;
  TH1I * h1 = new TH1I("hPhosDigitsMul", "Digits multiplicity distribution in PHOS;# of Digits;Entries", 2000, 0, 10000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, kDigitsMul, !expert, image) ;
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH2I * h0 = new TH2I("hRpPHOSxyMod1","RecPoints Rows x Columns for PHOS module 1;Row #;Column #", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h0,kRPmod1, expert, !image) ;
  TH2I * h1 = new TH2I("hRpPHOSxyMod2","RecPoints Rows x Columns for PHOS module 2;Row #;Column #", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h1,kRPmod2, expert, !image) ;
  TH2I * h2 = new TH2I("hRpPHOSxyMod3","RecPoints Rows x Columns for PHOS module 3;Row #;Column #", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h2,kRPmod3, expert, !image) ;
  TH2I * h3 = new TH2I("hRpPHOSxyMod4","RecPoints Rows x Columns for PHOS module 4;Row #;Column #", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h3,kRPmod4, expert, !image) ;
  TH2I * h4 = new TH2I("hRpPHOSxyMod5","RecPoints Rows x Columns for PHOS module 5;Row #;Column #", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h4,kRPmod5, expert, !image) ;
 
  TH1F * h5 = new TH1F("hEmcPhosRecPointsSpectrum",  "EMC RecPoints spectrum in PHOS;Energy [MeV];Counts",   2000, 0., 20.) ; 
  h5->Sumw2() ;
  Add2RecPointsList(h5, kRPSpec, !expert, image)  ;

  TH1I * h6 = new TH1I("hEmcPhosRecPointsMul", "EMC RecPoints multiplicity distribution in PHOS;# of EMC Clusters;Entries", 100, 0,  100) ; 
  h6->Sumw2() ;
  Add2RecPointsList(h6, kRPNtot, !expert, image) ;

  TH1I * h7 = new TH1I("hEmcPhosRecPointsEtot", "EMC RecPoints Etot;Energy [MeV];Counts", 200, 0,  200.) ; 
  h7->Sumw2() ;
  Add2RecPointsList(h7, kRPEtot, !expert, image) ;

  TH1I * h8 = new TH1I("hCpvPhosRecPointsMul", "CPV RecPoints multiplicity distribution in PHOS;# of CPV clusters;Counts", 100, 0,  100) ; 
  h8->Sumw2() ;
  Add2RecPointsList(h8, kRPNcpv, !expert, image) ;
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitRaws()
{
  // create Raws histograms in Raws subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH2I * h0 = new TH2I("hHighPHOSxyMod0","High Gain in PHOS module 0", 64, 0, 64, 56, 0, 56) ;
  h0->SetXTitle("x, cells"); h0->SetYTitle("z, cells");
  Add2RawsList(h0,kHGmod0, expert, !image, !saveCorr) ;
  TH2I * h1 = new TH2I("hHighPHOSxyMod1","High Gain in PHOS module 1", 64, 0, 64, 56, 0, 56) ;
  h1->SetXTitle("x, cells"); h1->SetYTitle("z, cells");
  Add2RawsList(h1,kHGmod1, expert, !image, !saveCorr) ;
  TH2I * h2 = new TH2I("hHighPHOSxyMod2","High Gain in PHOS module 2", 64, 0, 64, 56, 0, 56) ;
  h2->SetXTitle("x, cells"); h2->SetYTitle("z, cells");
  Add2RawsList(h2,kHGmod2, expert, !image, !saveCorr) ;
  TH2I * h3 = new TH2I("hHighPHOSxyMod3","High Gain in PHOS module 3", 64, 0, 64, 56, 0, 56) ;
  h3->SetXTitle("x, cells"); h3->SetYTitle("z, cells");
  Add2RawsList(h3,kHGmod3, expert, !image, !saveCorr) ;
  TH2I * h4 = new TH2I("hHighPHOSxyMod4","High Gain in PHOS module 4", 64, 0, 64, 56, 0, 56) ;
  h4->SetXTitle("x, cells"); h4->SetYTitle("z, cells");
  Add2RawsList(h4,kHGmod4, expert, !image, !saveCorr) ;

  TH2I * h5 = new TH2I("hLowPHOSxyMod0","Low Gain in PHOS module 0", 64, 0, 64, 56, 0, 56) ;
  h5->SetXTitle("x, cells"); h5->SetYTitle("z, cells");
  Add2RawsList(h5,kLGmod0, expert, !image, !saveCorr) ;
  TH2I * h6 = new TH2I("hLowPHOSxyMod1","Low Gain in PHOS module 1", 64, 0, 64, 56, 0, 56) ;
  h6->SetXTitle("x, cells"); h6->SetYTitle("z, cells");
  Add2RawsList(h6,kLGmod1, expert, !image, !saveCorr) ;
  TH2I * h7 = new TH2I("hLowPHOSxyMod2","Low Gain in PHOS module 2", 64, 0, 64, 56, 0, 56) ;
  h7->SetXTitle("x, cells"); h7->SetYTitle("z, cells");
  Add2RawsList(h7,kLGmod2, expert, !image, !saveCorr) ;
  TH2I * h8 = new TH2I("hLowPHOSxyMod3","Low Gain in PHOS module 3", 64, 0, 64, 56, 0, 56) ;
  h8->SetXTitle("x, cells"); h8->SetYTitle("z, cells");
  Add2RawsList(h8,kLGmod3, expert, !image, !saveCorr) ;
  TH2I * h9 = new TH2I("hLowPHOSxyMod4","Low Gain in PHOS module 4", 64, 0, 64, 56, 0, 56) ;
  h9->SetXTitle("x, cells"); h9->SetYTitle("z, cells");
  Add2RawsList(h9,kLGmod4, expert, !image, !saveCorr) ;

  TH1I * h10 = new TH1I("hLowPhosModules",    "Low Gain Hits in EMCA PHOS modules", 5, 0, 5) ;
  h10->SetXTitle("Module number");
  Add2RawsList(h10, kNmodLG, !expert, image, !saveCorr) ;
  TH1I * h11 = new TH1I("hHighPhosModules",   "High Gain Hits in EMCA PHOS modules",5, 0, 5) ;
  h11->SetXTitle("Module number");
  Add2RawsList(h11, kNmodHG, !expert, image, !saveCorr) ;

  TH1F * h12 = new TH1F("hLowPhosRawtime" , "Low Gain Time of raw hits in PHOS" , 500, -50., 200.) ;
  h12->SetXTitle("Time [samples]");
  h12->Sumw2() ;
  Add2RawsList(h12, kLGtime, expert, !image, !saveCorr) ;
  TH1F * h13 = new TH1F("hHighPhosRawtime", "High Gain Time of raw hits in PHOS", 500, -50., 200.) ;
  h13->SetXTitle("Time [samples]");
  h13->Sumw2() ;
  Add2RawsList(h13, kHGtime, expert, !image, !saveCorr) ;

  TH1F * h14 = new TH1F("hLowPhosRawEnergy" , "Low Gain Energy of raw hits in PHOS" , 512, 0., 1024.) ;
  h14->SetXTitle("Energy [ADC counts]");
  h14->Sumw2() ;
  Add2RawsList(h14, kSpecLG, !expert, image, !saveCorr) ;
  TH1F * h15 = new TH1F("hHighPhosRawEnergy", "High Gain Energy of raw hits in PHOS", 512, 0., 1024.) ;
  h15->SetXTitle("Energy [ADC counts]");
  h15->Sumw2() ;
  Add2RawsList(h15, kSpecHG, !expert, image, !saveCorr) ;

  TH1F * h16 = new TH1F("hLowNtot" , "Low Gain Total Number of raw hits in PHOS" , 500, 0., 5000.) ;
  h16->SetXTitle("Number of hits");
  h16->Sumw2() ;
  Add2RawsList(h16, kNtotLG, !expert, image, saveCorr) ;
  TH1F * h17 = new TH1F("hHighNtot", "High Gain Total Number of raw hits in PHOS", 500, 0., 5000.) ;
  h17->SetXTitle("Number of hits");
  h17->Sumw2() ;
  Add2RawsList(h17, kNtotHG, !expert, image, saveCorr) ;

  TH1F * h18 = new TH1F("hLowEtot" , "Low Gain Total Energy of raw hits in PHOS" , 500, 0., 100000.) ;
  h18->SetXTitle("Energy [ADC counts]");
  h18->Sumw2() ;
  Add2RawsList(h18, kEtotLG, !expert, image, saveCorr) ;
  TH1F * h19 = new TH1F("hHighEtot", "High Gain Total Energy of raw hits in PHOS", 500, 0., 100000.) ;
  h19->SetXTitle("Energy [ADC counts]");
  h19->Sumw2() ;
  Add2RawsList(h19, kEtotHG, !expert, image, saveCorr) ;

  TH2F * h20 = new TH2F("hQualHGxyMod0","High Gain signal quality in module 0", 64, 0, 64, 56, 0, 56) ;
  h20->SetXTitle("x, cells"); h20->SetYTitle("z, cells");
  Add2RawsList(h20,kHGqualMod0, expert, !image, !saveCorr) ;
  TH2F * h21 = new TH2F("hQualHGxyMod1","High Gain signal quality in module 1", 64, 0, 64, 56, 0, 56) ;
  h21->SetXTitle("x, cells"); h21->SetYTitle("z, cells");
  Add2RawsList(h21,kHGqualMod1, expert, !image, !saveCorr) ;
  TH2F * h22 = new TH2F("hQualHGxyMod2","High Gain signal quality in module 2", 64, 0, 64, 56, 0, 56) ;
  h22->SetXTitle("x, cells"); h22->SetYTitle("z, cells");
  Add2RawsList(h22,kHGqualMod2, expert, !image, !saveCorr) ;
  TH2F * h23 = new TH2F("hQualHGxyMod3","High Gain signal quality in module 3", 64, 0, 64, 56, 0, 56) ;
  h23->SetXTitle("x, cells"); h23->SetYTitle("z, cells");
  Add2RawsList(h23,kHGqualMod3, expert, !image, !saveCorr) ;
  TH2F * h24 = new TH2F("hQualHGxyMod4","High Gain signal quality in module 4", 64, 0, 64, 56, 0, 56) ;
  h24->SetXTitle("x, cells"); h24->SetYTitle("z, cells");
  Add2RawsList(h24,kHGqualMod4, expert, !image, !saveCorr) ;

  TH2F * h25 = new TH2F("hQualLGxyMod0","Low Gain signal quality in module 0", 64, 0, 64, 56, 0, 56) ;
  h25->SetXTitle("x, cells"); h25->SetYTitle("z, cells");
  Add2RawsList(h25,kLGqualMod0, expert, !image, !saveCorr) ;
  TH2F * h26 = new TH2F("hQualLGxyMod1","Low Gain signal quality in module 1", 64, 0, 64, 56, 0, 56) ;
  h26->SetXTitle("x, cells"); h26->SetYTitle("z, cells");
  Add2RawsList(h26,kLGqualMod1, expert, !image, !saveCorr) ;
  TH2F * h27 = new TH2F("hQualLGxyMod2","Low Gain signal quality in module 2", 64, 0, 64, 56, 0, 56) ;
  h27->SetXTitle("x, cells"); h27->SetYTitle("z, cells");
  Add2RawsList(h27,kLGqualMod2, expert, !image, !saveCorr) ;
  TH2F * h28 = new TH2F("hQualLGxyMod3","Low Gain signal quality in module 3", 64, 0, 64, 56, 0, 56) ;
  h28->SetXTitle("x, cells"); h28->SetYTitle("z, cells");
  Add2RawsList(h28,kLGqualMod3, expert, !image, !saveCorr) ;
  TH2F * h29 = new TH2F("hQualLGxyMod4","Low Gain signal quality in module 4", 64, 0, 64, 56, 0, 56) ;
  h29->SetXTitle("x, cells"); h29->SetYTitle("z, cells");
  Add2RawsList(h29,kLGqualMod4, expert, !image, !saveCorr) ;

  TH1F * h30 = new TH1F("hLGpedRMS" ,"Low Gain pedestal RMS" ,200,0.,20.) ;
  h30->SetXTitle("RMS [ADC counts]");
  h30->Sumw2() ;
  Add2RawsList(h30,kLGpedRMS, expert, !image, !saveCorr) ;
  TH1F * h31 = new TH1F("hHGpedRMS" ,"High Gain pedestal RMS",200,0.,20.) ;
  h31->SetXTitle("RMS [ADC counts]");
  h31->Sumw2() ;
  Add2RawsList(h31,kHGpedRMS, expert, !image, !saveCorr) ;

  TH2F * h32 = new TH2F("hpedRMSHGxyMod0","High Gain pedestal RMS in module 0", 64, 0, 64, 56, 0, 56) ;
  h32->SetXTitle("x, cells"); h32->SetYTitle("z, cells");
  Add2RawsList(h32,kHGpedRMSMod0, expert, !image, !saveCorr) ;
  TH2F * h33 = new TH2F("hpedRMSHGxyMod1","High Gain pedestal RMS in module 1", 64, 0, 64, 56, 0, 56) ;
  h33->SetXTitle("x, cells"); h33->SetYTitle("z, cells");
  Add2RawsList(h33,kHGpedRMSMod1, expert, !image, !saveCorr) ;
  TH2F * h34 = new TH2F("hpedRMSHGxyMod2","High Gain pedestal RMS in module 2", 64, 0, 64, 56, 0, 56) ;
  h34->SetXTitle("x, cells"); h34->SetYTitle("z, cells");
  Add2RawsList(h34,kHGpedRMSMod2, expert, !image, !saveCorr) ;
  TH2F * h35 = new TH2F("hpedRMSHGxyMod3","High Gain pedestal RMS in module 3", 64, 0, 64, 56, 0, 56) ;
  h35->SetXTitle("x, cells"); h35->SetYTitle("z, cells");
  Add2RawsList(h35,kHGpedRMSMod3, expert, !image, !saveCorr) ;
  TH2F * h36 = new TH2F("hpedRMSHGxyMod4","High Gain pedestal RMS in module 4", 64, 0, 64, 56, 0, 56) ;
  h36->SetXTitle("x, cells"); h36->SetYTitle("z, cells");
  Add2RawsList(h36,kHGpedRMSMod4, expert, !image, !saveCorr) ;

  TH2F * h37 = new TH2F("hpedRMSLGxyMod0","Low Gain pedestal RMS in module 0", 64, 0, 64, 56, 0, 56) ;
  h37->SetXTitle("x, cells"); h37->SetYTitle("z, cells");
  Add2RawsList(h37,kLGpedRMSMod0, expert, !image, !saveCorr) ;
  TH2F * h38 = new TH2F("hpedRMSLGxyMod1","Low Gain pedestal RMS in module 1", 64, 0, 64, 56, 0, 56) ;
  h38->SetXTitle("x, cells"); h38->SetYTitle("z, cells");
  Add2RawsList(h38,kLGpedRMSMod1, expert, !image, !saveCorr) ;
  TH2F * h39 = new TH2F("hpedRMSLGxyMod2","Low Gain pedestal RMS in module 2", 64, 0, 64, 56, 0, 56) ;
  h39->SetXTitle("x, cells"); h39->SetYTitle("z, cells");
  Add2RawsList(h39,kLGpedRMSMod2, expert, !image, !saveCorr) ;
  TH2F * h40 = new TH2F("hpedRMSLGxyMod3","Low Gain pedestal RMS in module 3", 64, 0, 64, 56, 0, 56) ;
  h40->SetXTitle("x, cells"); h40->SetYTitle("z, cells");
  Add2RawsList(h40,kLGpedRMSMod3, expert, !image, !saveCorr) ;
  TH2F * h41 = new TH2F("hpedRMSLGxyMod4","Low Gain pedestal RMS in module 4", 64, 0, 64, 56, 0, 56) ;
  h41->SetXTitle("x, cells"); h41->SetYTitle("z, cells");
  Add2RawsList(h41,kLGpedRMSMod4, expert, !image, !saveCorr) ;
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}

//____________________________________________________________________________
void AliPHOSQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs 
  Int_t nTot = 0 ; 
  Double_t eTot = 0 ; 
  for ( Int_t index = 0; index < esd->GetNumberOfCaloClusters() ; index++ ) {
    AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
    if( clu->IsPHOS() ) {
      FillESDsData(kESDSpec,clu->E()) ;
      const Double_t * pid = clu->GetPID() ;
      FillESDsData(kESDpid,pid[AliPID::kPhoton]) ;
      eTot+=clu->E() ;
      nTot++ ;
    } 
  }
  FillESDsData(kESDNtot,nTot) ;
  FillESDsData(kESDEtot,eTot) ;
  //
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  //
}

//____________________________________________________________________________
void AliPHOSQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //Fill prepared histograms with Raw digit properties

  rawReader->Reset() ;

  const TObjArray* maps = AliPHOSRecoParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  AliAltroMapping *mapping[20];
  for(Int_t i = 0; i < 20; i++) {
    mapping[i] = (AliAltroMapping*)maps->At(i);
  }

  AliCaloRawStreamV3 fRawStream(rawReader,"PHOS",mapping);

  AliPHOSRawFitterv0 * fitter ;
  if     (strcmp(GetRecoParam()->EMCFitterVersion(),"v1")==0)
    fitter=new AliPHOSRawFitterv1();
  else if(strcmp(GetRecoParam()->EMCFitterVersion(),"v2")==0)
    fitter=new AliPHOSRawFitterv2();
  else
    fitter=new AliPHOSRawFitterv0();
  Double_t lgEtot=0. ;
  Double_t hgEtot=0. ;
  Int_t    lgNtot=0 ;
  Int_t    hgNtot=0 ;

  while (fRawStream.NextDDL()) {
    while (fRawStream.NextChannel()) {
      Int_t module   = fRawStream.GetModule();
      Int_t cellX    = fRawStream.GetCellX();
      Int_t cellZ    = fRawStream.GetCellZ();
      Int_t caloFlag = fRawStream.GetCaloFlag(); // 0=LG, 1=HG, 2=TRU

      if(caloFlag!=0 && caloFlag!=1) continue; //TRU data!

      fitter->SetChannelGeo(module+1,cellX+1,cellZ+1,caloFlag);

      if(fitter->GetAmpOffset()==0 && fitter->GetAmpThreshold()==0) {
	Short_t altroCFG1 = fRawStream.GetAltroCFG1();
	Bool_t ZeroSuppressionEnabled = (altroCFG1 >> 15) & 0x1;
	if(ZeroSuppressionEnabled) {
	  Short_t offset = (altroCFG1 >> 10) & 0xf;
	  Short_t threshold = altroCFG1 & 0x3ff;
	  fitter->SubtractPedestals(kFALSE);
	  fitter->SetAmpOffset(offset);
	  fitter->SetAmpThreshold(threshold);
	}
	else
	  fitter->SubtractPedestals(kTRUE);
      }

      Int_t nBunches = 0;
      while (fRawStream.NextBunch()) {
	nBunches++;
	if (nBunches > 1) continue;
	const UShort_t *sig = fRawStream.GetSignals();
	Int_t sigStart      = fRawStream.GetStartTimeBin();
	Int_t sigLength     = fRawStream.GetBunchLength();
	fitter->Eval(sig,sigStart,sigLength);
      } // End of NextBunch()

      Double_t energy = fitter->GetEnergy() ; 
      Double_t time   = fitter->GetTime() ;

      if (caloFlag == 0) { // LG
	FillRawsData(kLGpedRMS,fitter->GetPedestalRMS()) ;
	switch(module){
        case 0: FillRawsData(kLGmod0,cellX,cellZ) ; break ;
        case 1: FillRawsData(kLGmod1,cellX,cellZ) ; break ;
        case 2: FillRawsData(kLGmod2,cellX,cellZ) ; break ;
        case 3: FillRawsData(kLGmod3,cellX,cellZ) ; break ;
        case 4: FillRawsData(kLGmod4,cellX,cellZ) ; break ;
	}
	switch (module){
        case 0: FillRawsData(kLGpedRMSMod0,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
        case 1: FillRawsData(kLGpedRMSMod1,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
        case 2: FillRawsData(kLGpedRMSMod2,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
        case 3: FillRawsData(kLGpedRMSMod3,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
        case 4: FillRawsData(kLGpedRMSMod4,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
	}
	//if quality was evaluated, fill histo
	if(strcmp(GetRecoParam()->EMCFitterVersion(),"v1")==0){
	  switch (module){
	  case 0: FillRawsData(kLGqualMod0,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  case 1: FillRawsData(kLGqualMod1,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  case 2: FillRawsData(kLGqualMod2,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  case 3: FillRawsData(kLGqualMod3,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  case 4: FillRawsData(kLGqualMod4,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  }
	}                                  
	FillRawsData(kNmodLG,module) ;
	FillRawsData(kLGtime,time) ; 
	FillRawsData(kSpecLG,energy) ;    
	lgEtot+=energy ;
	lgNtot++ ;   
      }
      else if (caloFlag == 1) { // HG        
	FillRawsData(kHGpedRMS,fitter->GetPedestalRMS()) ;
	switch (module){
	case 0: FillRawsData(kHGmod0,cellX,cellZ) ; break ;
	case 1: FillRawsData(kHGmod1,cellX,cellZ) ; break ;
	case 2: FillRawsData(kHGmod2,cellX,cellZ) ; break ;
	case 3: FillRawsData(kHGmod3,cellX,cellZ) ; break ;
	case 4: FillRawsData(kHGmod4,cellX,cellZ) ; break ;
	}
	switch (module){
	case 0: FillRawsData(kHGpedRMSMod0,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
	case 1: FillRawsData(kHGpedRMSMod1,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
	case 2: FillRawsData(kHGpedRMSMod2,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
	case 3: FillRawsData(kHGpedRMSMod3,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
	case 4: FillRawsData(kHGpedRMSMod4,cellX,cellZ,fitter->GetPedestalRMS()) ; break ;
	}               
	//if quality was evaluated, fill histo
	if(strcmp(GetRecoParam()->EMCFitterVersion(),"v1")==0){
	  switch (module){
	  case 0: FillRawsData(kHGqualMod0,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  case 1: FillRawsData(kHGqualMod1,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  case 2: FillRawsData(kHGqualMod2,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  case 3: FillRawsData(kHGqualMod3,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  case 4: FillRawsData(kHGqualMod4,cellX,cellZ,fitter->GetSignalQuality()) ; break ;
	  }	  
	}
	FillRawsData(kNmodHG,module) ; 
	FillRawsData(kHGtime,time) ;  
	FillRawsData(kSpecHG,energy) ;
	hgEtot+=energy ; 
	hgNtot++ ;  
      }
    }  // End of NextChannel
  } // End of NextDDL
  delete fitter;

  FillRawsData(kEtotLG,lgEtot) ; 
  TParameter<double> * p;
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->
					FindObject(Form("%s_%s_%s", GetName(), 
							AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), 
							GetRawsData(kEtotLG)->GetName()))) ; 
  if (p) p->SetVal(lgEtot) ; 
  FillRawsData(kEtotHG,hgEtot) ;  
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->
					FindObject(Form("%s_%s_%s", GetName(), 
							AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), 
							GetRawsData(kEtotHG)->GetName()))) ; 
  if (p) p->SetVal(hgEtot) ; 
  FillRawsData(kNtotLG,lgNtot) ;
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->
					FindObject(Form("%s_%s_%s", GetName(), 
							AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), 
							GetRawsData(kNtotLG)->GetName()))) ; 
  if (p) p->SetVal(lgNtot) ; 
  FillRawsData(kNtotHG,hgNtot) ;
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->
					FindObject(Form("%s_%s_%s", 
							GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), 
							GetRawsData(kNtotHG)->GetName()))) ; 
  if (p) p->SetVal(hgNtot) ; 
  //
  IncEvCountCycleRaws();
  IncEvCountTotalRaws();
  //
}

//____________________________________________________________________________
void AliPHOSQADataMakerRec::MakeDigits()
{
  // makes data from Digits
  
  if ( ! GetDigitsData(kDigitsMul) ) InitDigits() ;
  FillDigitsData(kDigitsMul,fDigitsArray->GetEntriesFast()) ; 
  TIter next(fDigitsArray) ; 
  AliPHOSDigit * digit ; 
  while ( (digit = dynamic_cast<AliPHOSDigit *>(next())) ) {
    FillDigitsData(kDigits, digit->GetEnergy()) ;
  }  
  //
}

//____________________________________________________________________________
void AliPHOSQADataMakerRec::MakeDigits(TTree * digitTree)
{
  // makes data from Digit Tree
  if (fDigitsArray) 
    fDigitsArray->Clear() ; 
  else 
    fDigitsArray = new TClonesArray("AliPHOSDigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("PHOS") ;
  if ( ! branch ) {AliWarning("PHOS branch in Digit Tree not found"); return;} 
  branch->SetAddress(&fDigitsArray) ;
  branch->GetEntry(0) ; 
  MakeDigits() ; 
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
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
    
    FillRecPointsData(kRPNtot,emcrecpoints->GetEntriesFast()) ; 
    TIter next(emcrecpoints) ; 
    AliPHOSEmcRecPoint * rp ; 
    Double_t eTot = 0. ; 
    while ( (rp = static_cast<AliPHOSEmcRecPoint *>(next())) ) {
      FillRecPointsData(kRPSpec, rp->GetEnergy()) ;
      Int_t mod = rp->GetPHOSMod() ;
      TVector3 pos ;
      rp->GetLocalPosition(pos) ;
      switch(mod){
        case 1: FillRecPointsData(kRPmod1,pos.X(),pos.Z()) ; break ;
        case 2: FillRecPointsData(kRPmod2,pos.X(),pos.Z()) ; break ;
        case 3: FillRecPointsData(kRPmod3,pos.X(),pos.Z()) ; break ;
        case 4: FillRecPointsData(kRPmod4,pos.X(),pos.Z()) ; break ;
        case 5: FillRecPointsData(kRPmod5,pos.X(),pos.Z()) ; break ;
      }
      eTot+= rp->GetEnergy() ;
    }
    FillRecPointsData(kRPEtot,eTot) ;
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
    
    FillRecPointsData(kRPNcpv,cpvrecpoints->GetEntriesFast()) ; 
    cpvrecpoints->Delete();
    delete cpvrecpoints;
  }
  //
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
