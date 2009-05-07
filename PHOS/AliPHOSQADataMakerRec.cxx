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
#include "AliPHOSRawDecoderv1.h"
#include "AliPHOSRawDecoderv2.h"
#include "AliPHOSReconstructor.h"

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
  if(GetRawsData(kHGqualMod1) && GetRawsData(kHGmod1))
    GetRawsData(kHGqualMod1)->Divide( GetRawsData(kHGmod1) ) ;
  if(GetRawsData(kHGqualMod2) && GetRawsData(kHGmod2))
    GetRawsData(kHGqualMod2)->Divide( GetRawsData(kHGmod2) ) ;
  if(GetRawsData(kHGqualMod3) && GetRawsData(kHGmod3))
    GetRawsData(kHGqualMod3)->Divide( GetRawsData(kHGmod3) ) ;
  if(GetRawsData(kHGqualMod4) && GetRawsData(kHGmod4))
    GetRawsData(kHGqualMod4)->Divide( GetRawsData(kHGmod4) ) ;
  if(GetRawsData(kHGqualMod5) && GetRawsData(kHGmod5))
    GetRawsData(kHGqualMod5)->Divide( GetRawsData(kHGmod5) ) ;
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kPHOS, task, list) ;  
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
 
  Bool_t expert   = kTRUE ; 

  TH1F * h1 = new TH1F("hESDPhosSpectrum",  "ESDs spectrum in PHOS"                ,  200, 0.,   20.) ;
  h1->Sumw2() ;
  Add2ESDsList(h1, kESDSpec, !expert)  ;

  TH1I * h2 = new TH1I("hESDPhosMul",       "ESDs multiplicity distribution in PHOS", 100, 0,   100 ) ; 
  h2->Sumw2() ;
  Add2ESDsList(h2, kESDNtot, !expert) ;
 
  TH1F * h3 = new TH1F("hESDPhosEtot",      "ESDs total energy"                     , 2000, 0,  200.) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, kESDEtot, !expert) ;  //Expert histo
 
  TH1F * h4 = new TH1F("hESDpid",           "ESDs PID distribution in PHOS"         , 100, 0.,    1.) ;
  h4->Sumw2() ;
  Add2ESDsList(h4, kESDpid, !expert) ; //Expert histo
	
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  Bool_t expert   = kTRUE ; 
  
  TH2I * h0 = new TH2I("hRpPHOSxyMod1","RecPoints Rows x Columns for PHOS module 1", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h0,kRPmod1, expert) ;
  TH2I * h1 = new TH2I("hRpPHOSxyMod2","RecPoints Rows x Columns for PHOS module 2", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h1,kRPmod2, expert) ;
  TH2I * h2 = new TH2I("hRpPHOSxyMod3","RecPoints Rows x Columns for PHOS module 3", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h2,kRPmod3, expert) ;
  TH2I * h3 = new TH2I("hRpPHOSxyMod4","RecPoints Rows x Columns for PHOS module 4", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h3,kRPmod4, expert) ;
  TH2I * h4 = new TH2I("hRpPHOSxyMod5","RecPoints Rows x Columns for PHOS module 5", 64, -72., 72., 56, -63., 63.) ;                             
  Add2RecPointsList(h4,kRPmod5, expert) ;
 
  TH1F * h5 = new TH1F("hEmcPhosRecPointsSpectrum",  "EMC RecPoints spectrum in PHOS",   2000, 0., 20.) ; 
  h5->Sumw2() ;
  Add2RecPointsList(h5, kRPSpec, !expert)  ;

  TH1I * h6 = new TH1I("hEmcPhosRecPointsMul", "EMCA RecPoints multiplicity distribution in PHOS", 100, 0,  100) ; 
  h6->Sumw2() ;
  Add2RecPointsList(h6, kRPNtot, !expert) ;

  TH1I * h7 = new TH1I("hEmcPhosRecPointsEtot", "EMC RecPoints Etot", 200, 0,  200.) ; 
  h7->Sumw2() ;
  Add2RecPointsList(h7, kRPEtot, !expert) ;

  TH1I * h8 = new TH1I("hCpvPhosRecPointsMul", "CPV RecPoints multiplicity distribution in PHOS", 100, 0,  100) ; 
  h8->Sumw2() ;
  Add2RecPointsList(h8, kRPNcpv, !expert) ;
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::InitRaws()
{
  // create Raws histograms in Raws subdir
  Bool_t expert   = kTRUE ; 
  Bool_t saveCorr = kTRUE ; 
  
  TH2I * h0 = new TH2I("hHighPHOSxyMod1","High Gain Rows x Columns for PHOS module 1", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h0,kHGmod1, expert, !saveCorr) ;
  TH2I * h1 = new TH2I("hHighPHOSxyMod2","High Gain Rows x Columns for PHOS module 2", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h1,kHGmod2, expert, !saveCorr) ;
  TH2I * h2 = new TH2I("hHighPHOSxyMod3","High Gain Rows x Columns for PHOS module 3", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h2,kHGmod3, expert, !saveCorr) ;
  TH2I * h3 = new TH2I("hHighPHOSxyMod4","High Gain Rows x Columns for PHOS module 4", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h3,kHGmod4, expert, !saveCorr) ;
  TH2I * h4 = new TH2I("hHighPHOSxyMod5","High Gain Rows x Columns for PHOS module 5", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h4,kHGmod5, expert, !saveCorr) ;
  TH2I * h5 = new TH2I("hLowPHOSxyMod1","Low Gain Rows x Columns for PHOS module 1", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h5,kLGmod1, expert, !saveCorr) ;
  TH2I * h6 = new TH2I("hLowPHOSxyMod2","Low Gain Rows x Columns for PHOS module 2", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h6,kLGmod2, expert, !saveCorr) ;
  TH2I * h7 = new TH2I("hLowPHOSxyMod3","Low Gain Rows x Columns for PHOS module 3", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h7,kLGmod3, expert, !saveCorr) ;
  TH2I * h8 = new TH2I("hLowPHOSxyMod4","Low Gain Rows x Columns for PHOS module 4", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h8,kLGmod4, expert, !saveCorr) ;
  TH2I * h9 = new TH2I("hLowPHOSxyMod5","Low Gain Rows x Columns for PHOS module 5", 64, 0, 64, 56, 0, 56) ;
  Add2RawsList(h9,kLGmod5, expert, !saveCorr) ;

  TH1I * h10 = new TH1I("hLowPhosModules",    "Low Gain Hits in EMCA PHOS modules",       6, 0, 6) ;
  h10->Sumw2() ;
  Add2RawsList(h10, kNmodLG, !expert, !saveCorr) ;
  TH1I * h11 = new TH1I("hHighPhosModules",   "High Gain Hits in EMCA PHOS modules",       6, 0, 6) ;
  h11->Sumw2() ;
  Add2RawsList(h11, kNmodHG, !expert, !saveCorr) ;

  TH1F * h12 = new TH1F("hLowPhosRawtime", "Low Gain Time of raw hits in PHOS", 500, -50., 200.) ;
  h12->Sumw2() ;
  Add2RawsList(h12, kLGtime, !expert, !saveCorr) ;
  TH1F * h13 = new TH1F("hHighPhosRawtime", "High Gain Time of raw hits in PHOS", 500, -50., 200.) ;
  h13->Sumw2() ;
  Add2RawsList(h13, kHGtime, !expert, !saveCorr) ;

  TH1F * h14 = new TH1F("hLowPhosRawEnergy", "Low Gain Energy of raw hits in PHOS", 500, 0., 1000.) ;
  h14->Sumw2() ;
  Add2RawsList(h14, kSpecLG, !expert, !saveCorr) ;
  TH1F * h15 = new TH1F("hHighPhosRawEnergy", "High Gain Energy of raw hits in PHOS",500,0., 1000.) ;
  h15->Sumw2() ;
  Add2RawsList(h15, kSpecHG, !expert, !saveCorr) ;

  TH1F * h16 = new TH1F("hLowNtot", "Low Gain Total Number of raw hits in PHOS", 500, 0., 5000.) ;
  h16->Sumw2() ;
  Add2RawsList(h16, kNtotLG, !expert, saveCorr) ;
  TH1F * h17 = new TH1F("hHighNtot", "High Gain Total Number of raw hits in PHOS",500,0., 5000.) ;
  h17->Sumw2() ;
  Add2RawsList(h17, kNtotHG, !expert, saveCorr) ;

  TH1F * h18 = new TH1F("hLowEtot", "Low Gain Total Energy of raw hits in PHOS", 500, 0., 5000.) ;
  h18->Sumw2() ;
  Add2RawsList(h18, kEtotLG, !expert, saveCorr) ;
  TH1F * h19 = new TH1F("hHighEtot", "High Gain Total Energy of raw hits in PHOS",500,0., 100000.) ;
  h19->Sumw2() ;
  Add2RawsList(h19, kEtotHG, !expert, saveCorr) ;

  TH2F * h20 = new TH2F("hQualHGxyMod1","High Gain signal quality Rows x Columns module 1", 64, 0, 64, 56, 0, 56) ;
  h20->SetOption("colz");
  Add2RawsList(h20,kHGqualMod1, expert, !saveCorr) ;
  TH2F * h21 = new TH2F("hQualHGxyMod2","High Gain signal quality Rows x Columns module 2", 64, 0, 64, 56, 0, 56) ;
  h21->SetOption("colz");
  Add2RawsList(h21,kHGqualMod2, expert, !saveCorr) ;
  TH2F * h22 = new TH2F("hQualHGxyMod3","High Gain signal quality Rows x Columns module 3", 64, 0, 64, 56, 0, 56) ;
  h22->SetOption("colz");
  Add2RawsList(h22,kHGqualMod3, expert, !saveCorr) ;
  TH2F * h23 = new TH2F("hQualHGxyMod4","High Gain signal quality Rows x Columns module 4", 64, 0, 64, 56, 0, 56) ;
  h23->SetOption("colz");
  Add2RawsList(h23,kHGqualMod4, expert, !saveCorr) ;
  TH2F * h24 = new TH2F("hQualHGxyMod5","High Gain signal quality Rows x Columns module 5", 64, 0, 64, 56, 0, 56) ;
  h24->SetOption("colz");
  Add2RawsList(h24,kHGqualMod5, expert, !saveCorr) ;
  TH2F * h25 = new TH2F("hQualLGxyMod1","Low Gain signal quality Rows x Columns module 1", 64, 0, 64, 56, 0, 56) ;
  h25->SetOption("colz");
  Add2RawsList(h25,kLGqualMod1, expert, !saveCorr) ;
  TH2F * h26 = new TH2F("hQualLGxyMod2","Low Gain signal quality Rows x Columns module 2", 64, 0, 64, 56, 0, 56) ;
  h26->SetOption("colz");
  Add2RawsList(h26,kLGqualMod2, expert, !saveCorr) ;
  TH2F * h27 = new TH2F("hQualLGxyMod3","Low Gain signal quality Rows x Columns module 3", 64, 0, 64, 56, 0, 56) ;
  h27->SetOption("colz");
  Add2RawsList(h27,kLGqualMod3, expert, !saveCorr) ;
  TH2F * h28 = new TH2F("hQualLGxyMod4","Low Gain signal quality Rows x Columns module 4", 64, 0, 64, 56, 0, 56) ;
  h28->SetOption("colz");
  Add2RawsList(h28,kLGqualMod4, expert, !saveCorr) ;
  TH2F * h29 = new TH2F("hQualLGxyMod5","Low Gain signal quality Rows x Columns module 5", 64, 0, 64, 56, 0, 56) ;
  h29->SetOption("colz");
  Add2RawsList(h29,kLGqualMod5, expert, !saveCorr) ;

  TH1F * h30 = new TH1F("hLGpedRMS","Low Gain pedestal RMS",200,0.,20.) ;
  h30->Sumw2() ;
  Add2RawsList(h30,kLGpedRMS, !expert, !saveCorr) ;
  TH1F * h31 = new TH1F("hHGpedRMS","High Gain pedestal RMS",200,0.,20.) ;
  h31->Sumw2() ;
  Add2RawsList(h31,kHGpedRMS, !expert, !saveCorr) ;

  TH2F * h32 = new TH2F("hpedRMSHGxyMod1","High Gain pedestal RMS Rows x Columns module 1", 64, 0, 64, 56, 0, 56) ;
  h32->SetOption("colz");
  Add2RawsList(h32,kHGpedRMSMod1, expert, !saveCorr) ;
  TH2F * h33 = new TH2F("hpedRMSHGxyMod2","High Gain pedestal RMS Rows x Columns module 2", 64, 0, 64, 56, 0, 56) ;
  h33->SetOption("colz");
  Add2RawsList(h33,kHGpedRMSMod2, expert, !saveCorr) ;
  TH2F * h34 = new TH2F("hpedRMSHGxyMod3","High Gain pedestal RMS Rows x Columns module 3", 64, 0, 64, 56, 0, 56) ;
  h34->SetOption("colz");
  Add2RawsList(h34,kHGpedRMSMod3, expert, !saveCorr) ;
  TH2F * h35 = new TH2F("hpedRMSHGxyMod4","High Gain pedestal RMS Rows x Columns module 4", 64, 0, 64, 56, 0, 56) ;
  h35->SetOption("colz");
  Add2RawsList(h35,kHGpedRMSMod4, expert, !saveCorr) ;
  TH2F * h36 = new TH2F("hpedRMSHGxyMod5","High Gain pedestal RMS Rows x Columns module 5", 64, 0, 64, 56, 0, 56) ;
  h36->SetOption("colz");
  Add2RawsList(h36,kHGpedRMSMod5, expert, !saveCorr) ;
  TH2F * h37 = new TH2F("hpedRMSLGxyMod1","Low Gain pedestal RMS Rows x Columns module 1", 64, 0, 64, 56, 0, 56) ;
  h37->SetOption("colz");
  Add2RawsList(h37,kLGpedRMSMod1, expert, !saveCorr) ;
  TH2F * h38 = new TH2F("hpedRMSLGxyMod2","Low Gain pedestal RMS Rows x Columns module 2", 64, 0, 64, 56, 0, 56) ;
  h38->SetOption("colz");
  Add2RawsList(h38,kLGpedRMSMod2, expert, !saveCorr) ;
  TH2F * h39 = new TH2F("hpedRMSLGxyMod3","Low Gain pedestal RMS Rows x Columns module 3", 64, 0, 64, 56, 0, 56) ;
  h39->SetOption("colz");
  Add2RawsList(h39,kLGpedRMSMod3, expert, !saveCorr) ;
  TH2F * h40 = new TH2F("hpedRMSLGxyMod4","Low Gain pedestal RMS Rows x Columns module 4", 64, 0, 64, 56, 0, 56) ;
  h40->SetOption("colz");
  Add2RawsList(h40,kLGpedRMSMod4, expert, !saveCorr) ;
  TH2F * h41 = new TH2F("hpedRMSLGxyMod5","Low Gain pedestal RMS Rows x Columns module 5", 64, 0, 64, 56, 0, 56) ;
  h41->SetOption("colz");
  Add2RawsList(h41,kLGpedRMSMod5, expert, !saveCorr) ;

  /*
  TH1F * h42 = new TH1F("hLGpedMean","Low Gain pedestal Mean",200,0.,20.) ;
  h42->Sumw2() ;
  Add2RawsList(h42,kLGpedMean, !expert, !saveCorr) ;
  TH1F * h43 = new TH1F("hHGpedMean","High Gain pedestal Mean",200,0.,20.) ;
  h43->Sumw2() ;
  Add2RawsList(h43,kHGpedMean, !expert, !saveCorr) ;

  TH2F * h44 = new TH2F("hpedMeanHGxyMod1","High Gain pedestal Mean Rows x Columns module 1", 64, 0, 64, 56, 0, 56) ;
  h44->SetOption("colz");
  Add2RawsList(h44,kHGpedMeanMod1, expert, !saveCorr) ;
  TH2F * h45 = new TH2F("hpedMeanHGxyMod2","High Gain pedestal Mean Rows x Columns module 2", 64, 0, 64, 56, 0, 56) ;
  h45->SetOption("colz");
  Add2RawsList(h45,kHGpedMeanMod2, expert, !saveCorr) ;
  TH2F * h46 = new TH2F("hpedMeanHGxyMod3","High Gain pedestal Mean Rows x Columns module 3", 64, 0, 64, 56, 0, 56) ;
  h46->SetOption("colz");
  Add2RawsList(h46,kHGpedMeanMod3, expert, !saveCorr) ;
  TH2F * h47 = new TH2F("hpedMeanHGxyMod4","High Gain pedestal Mean Rows x Columns module 4", 64, 0, 64, 56, 0, 56) ;
  h47->SetOption("colz");
  Add2RawsList(h47,kHGpedMeanMod4, expert, !saveCorr) ;
  TH2F * h48 = new TH2F("hpedMeanHGxyMod5","High Gain pedestal Mean Rows x Columns module 5", 64, 0, 64, 56, 0, 56) ;
  h48->SetOption("colz");
  Add2RawsList(h48,kHGpedMeanMod5, expert, !saveCorr) ;
  TH2F * h49 = new TH2F("hpedMeanLGxyMod1","Low Gain pedestal Mean Rows x Columns module 1", 64, 0, 64, 56, 0, 56) ;
  h49->SetOption("colz");
  Add2RawsList(h49,kLGpedMeanMod1, expert, !saveCorr) ;
  TH2F * h50 = new TH2F("hpedMeanLGxyMod2","Low Gain pedestal Mean Rows x Columns module 2", 64, 0, 64, 56, 0, 56) ;
  h50->SetOption("colz");
  Add2RawsList(h50,kLGpedMeanMod2, expert, !saveCorr) ;
  TH2F * h51 = new TH2F("hpedMeanLGxyMod3","Low Gain pedestal Mean Rows x Columns module 3", 64, 0, 64, 56, 0, 56) ;
  h51->SetOption("colz");
  Add2RawsList(h51,kLGpedMeanMod3, expert, !saveCorr) ;
  TH2F * h52 = new TH2F("hpedMeanLGxyMod4","Low Gain pedestal Mean Rows x Columns module 4", 64, 0, 64, 56, 0, 56) ;
  h52->SetOption("colz");
  Add2RawsList(h52,kLGpedMeanMod4, expert, !saveCorr) ;
  TH2F * h53 = new TH2F("hpedMeanLGxyMod5","Low Gain pedestal Mean Rows x Columns module 5", 64, 0, 64, 56, 0, 56) ;
  h53->SetOption("colz");
  Add2RawsList(h53,kLGpedMeanMod5, expert, !saveCorr) ;
  */
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
void AliPHOSQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //Fill prepared histograms with Raw digit properties
  rawReader->Reset() ;
  AliPHOSRawDecoder * decoder ;
  if(strcmp(GetRecoParam()->EMCDecoderVersion(),"v1")==0)
    decoder=new AliPHOSRawDecoderv1(rawReader);
  else
    if(strcmp(GetRecoParam()->EMCDecoderVersion(),"v2")==0)
      decoder=new AliPHOSRawDecoderv2(rawReader);
    else
      decoder=new AliPHOSRawDecoder(rawReader);
  //decoder->SubtractPedestals(GetRecoParam()->EMCSubtractPedestals());
  decoder->SubtractPedestals(kTRUE);
  Double_t lgEtot=0. ;
  Double_t hgEtot=0. ;
  Int_t lgNtot=0 ;
  Int_t hgNtot=0 ;

  while (decoder->NextDigit()) {
   Int_t module  = decoder->GetModule() ;
   Int_t row     = decoder->GetRow() ;
   Int_t col     = decoder->GetColumn() ;
   Double_t time = decoder->GetTime() ;
   Double_t energy  = decoder->GetEnergy() ;
   Bool_t lowGain = decoder->IsLowGain();
   if (lowGain) {
     //if(GetRecoParam()->EMCSubtractPedestals())
       GetRawsData(kLGpedRMS)->Fill(decoder->GetPedestalRMS()) ;
     if(energy<2.)
       continue ;
     switch(module){
        case 1: GetRawsData(kLGmod1)->Fill(row-0.5,col-0.5) ; break ;
        case 2: GetRawsData(kLGmod2)->Fill(row-0.5,col-0.5) ; break ;
        case 3: GetRawsData(kLGmod3)->Fill(row-0.5,col-0.5) ; break ;
        case 4: GetRawsData(kLGmod4)->Fill(row-0.5,col-0.5) ; break ;
        case 5: GetRawsData(kLGmod5)->Fill(row-0.5,col-0.5) ; break ;
     }
     switch (module){
        case 1: ((TH2F*)GetRawsData(kLGpedRMSMod1))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
        case 2: ((TH2F*)GetRawsData(kLGpedRMSMod2))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
        case 3: ((TH2F*)GetRawsData(kLGpedRMSMod3))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
        case 4: ((TH2F*)GetRawsData(kLGpedRMSMod4))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
        case 5: ((TH2F*)GetRawsData(kLGpedRMSMod5))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
     }
     //if quality was evaluated, fill histo
     if(strcmp(GetRecoParam()->EMCDecoderVersion(),"v1")==0){
       switch (module){
         case 1: ((TH2F*)GetRawsData(kLGqualMod1))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
         case 2: ((TH2F*)GetRawsData(kLGqualMod2))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
         case 3: ((TH2F*)GetRawsData(kLGqualMod3))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
         case 4: ((TH2F*)GetRawsData(kLGqualMod4))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
         case 5: ((TH2F*)GetRawsData(kLGqualMod5))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
       }
     }                                  
     GetRawsData(kNmodLG)->Fill(module) ;
     GetRawsData(kLGtime)->Fill(time) ; 
     GetRawsData(kSpecLG)->Fill(energy) ;    
     lgEtot+=energy ;
     lgNtot++ ;   
   } else {        
     //if this isnon-ZS run - fill pedestal RMS
     //if(GetRecoParam()->EMCSubtractPedestals())
       GetRawsData(kHGpedRMS)->Fill(decoder->GetPedestalRMS()) ;
     if(energy<8.)
       continue ;
     switch (module){
         case 1: GetRawsData(kHGmod1)->Fill(row-0.5,col-0.5) ; break ;
         case 2: GetRawsData(kHGmod2)->Fill(row-0.5,col-0.5) ; break ;
         case 3: GetRawsData(kHGmod3)->Fill(row-0.5,col-0.5) ; break ;
         case 4: GetRawsData(kHGmod4)->Fill(row-0.5,col-0.5) ; break ;
         case 5: GetRawsData(kHGmod5)->Fill(row-0.5,col-0.5) ; break ;
     }
     switch (module){
         case 1: ((TH2F*)GetRawsData(kHGpedRMSMod1))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
         case 2: ((TH2F*)GetRawsData(kHGpedRMSMod2))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
         case 3: ((TH2F*)GetRawsData(kHGpedRMSMod3))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
         case 4: ((TH2F*)GetRawsData(kHGpedRMSMod4))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
         case 5: ((TH2F*)GetRawsData(kHGpedRMSMod5))->Fill(row-0.5,col-0.5,decoder->GetPedestalRMS()) ; break ;
     }               
     //if quality was evaluated, fill histo
     if(strcmp(GetRecoParam()->EMCDecoderVersion(),"v1")==0){
       switch (module){
         case 1: ((TH2F*)GetRawsData(kHGqualMod1))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
         case 2: ((TH2F*)GetRawsData(kHGqualMod2))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
         case 3: ((TH2F*)GetRawsData(kHGqualMod3))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
         case 4: ((TH2F*)GetRawsData(kHGqualMod4))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
         case 5: ((TH2F*)GetRawsData(kHGqualMod5))->Fill(row-0.5,col-0.5,decoder->GetSampleQuality()) ; break ;
       }
     }
     GetRawsData(kNmodHG)->Fill(module) ; 
     GetRawsData(kHGtime)->Fill(time) ;  
     GetRawsData(kSpecHG)->Fill(energy) ;
     hgEtot+=energy ; 
     hgNtot++ ;  
   }                 
  }                    
  delete decoder;
  GetRawsData(kEtotLG)->Fill(lgEtot) ; 
  TParameter<double> * p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kEtotLG)->GetName()))) ; 
  p->SetVal(lgEtot) ; 
  GetRawsData(kEtotHG)->Fill(hgEtot) ;  
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kEtotHG)->GetName()))) ; 
  p->SetVal(hgEtot) ; 
  GetRawsData(kNtotLG)->Fill(lgNtot) ;
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kNtotLG)->GetName()))) ; 
  p->SetVal(lgNtot) ; 
  GetRawsData(kNtotHG)->Fill(hgNtot) ;
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kNtotHG)->GetName()))) ; 
  p->SetVal(hgNtot) ; 
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
    
    GetRecPointsData(kRPNtot)->Fill(emcrecpoints->GetEntriesFast()) ; 
    TIter next(emcrecpoints) ; 
    AliPHOSEmcRecPoint * rp ; 
    Double_t eTot = 0. ; 
    while ( (rp = dynamic_cast<AliPHOSEmcRecPoint *>(next())) ) {
      GetRecPointsData(kRPSpec)->Fill( rp->GetEnergy()) ;
      Int_t mod = rp->GetPHOSMod() ;
      TVector3 pos ;
      rp->GetLocalPosition(pos) ;
      switch(mod){
        case 1: GetRecPointsData(kRPmod1)->Fill(pos.X(),pos.Z()) ; break ;
        case 2: GetRecPointsData(kRPmod2)->Fill(pos.X(),pos.Z()) ; break ;
        case 3: GetRecPointsData(kRPmod3)->Fill(pos.X(),pos.Z()) ; break ;
        case 4: GetRecPointsData(kRPmod4)->Fill(pos.X(),pos.Z()) ; break ;
        case 5: GetRecPointsData(kRPmod5)->Fill(pos.X(),pos.Z()) ; break ;
      }
      eTot+= rp->GetEnergy() ;
    }
    GetRecPointsData(kRPEtot)->Fill(eTot) ;
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
    
    GetRecPointsData(kRPNcpv)->Fill(cpvrecpoints->GetEntriesFast()) ; 
    cpvrecpoints->Delete();
    delete cpvrecpoints;
  }
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
