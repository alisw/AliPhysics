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
//  Produces the data needed to calculate the TOF quality assurance. //
//  QA objects are 1 & 2 Dimensional histograms.                     //
//  author: S.Arcelli                                                //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/*
Modified by fbellini on 18/01/2011
- reduced histo binning to reduce size 
- added decoding errors plot
- added channel maps and options for DQM shifters
- new list of recPoints and ESDs plots
- removed hTrmChannels035 and hTrmChannels3671

 Modified by bvonhall on 03/11/2010
   - modified declaration of hTrmChannels035 and hTrmChannels3671 in EndOfDetectorCycle()
     to prevent memory corruption
   
 Modified by adecaro on 18/10/2010
   - fTOFRawStream object set as private member
 
Modified by fbellini on 13/09/2010
  - Set TLines as private members
  - Set image flag for expert histos

Modified by fbellini on 14/06/2010
  - Updated plots
  - use LoadRawDataBuffersV2()

  Modified by fbellini on 10/05/2010
  - Fixed EndOfDetectorCycle() memory corruption bug

  Modified by fbellini on 22/04/2010
   - Added filter for physics events

   Modified by fbellini on 16/04/2010
   - Added EnableDqmShifterOpt() 
   - Modified EndOfDetectorCycle() with options for DQM		
   - Updated ESDs QA

   Modified by fbellini on 30/03/2010
   - Changed raws time histos range to 610ns
   - Added FilterLTMData() and FilterSpare() methods
   - Added check on enabled channels for raw data 		
   - Updated RecPoints QA

   Modified by fbellini on 02/03/2010
   - Fixed raw data decoding methods (use AliTOFRawStream::LoadRawDataBuffer())
   - Added filter for noisy channels and read map from OCDB
   - Added GetCalibData() method
   - Added CheckVolumeID() and CheckEquipID() methods  
   - Updated Raw QA
*/

#include <TClonesArray.h>
#include <TH1F.h> 
#include <TH2F.h> 
#include <TLine.h>
#include <TPaveText.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliTOFRawStream.h"
#include "AliTOFcluster.h"
#include "AliTOFQADataMakerRec.h"
#include "AliTOFrawData.h"
#include "AliTOFGeometry.h"
#include "AliTOFChannelOnlineStatusArray.h"
#include "AliTOFDecoderSummaryData.h"
#include "AliTOFDecoderV2.h"

ClassImp(AliTOFQADataMakerRec)
           
//____________________________________________________________________________ 
  AliTOFQADataMakerRec::AliTOFQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kTOF), "TOF Quality Assurance Data Maker"),
  fCalibData(0x0),
  fEnableNoiseFiltering(kFALSE),
  fEnableDqmShifterOpt(kFALSE),
  fProcessedRawEventN(0),
  fIsSOC(kFALSE),
  fLineExpTimeMin(new TLine(200., 0., 200., 0.)),
  fLineExpTimeMax(new TLine(250., 0., 250., 0.)),
  fLineExpTotMin(new TLine(5., 0., 5., 0.)),
  fLineExpTotMax(new TLine(20., 0., 20., 0.)),
  fTOFRawStream(AliTOFRawStream()),
  fDecoderSummary(new AliTOFDecoderSummaryData())
{
  //
  // ctor
  //   
  for (Int_t sm=0;sm<17;sm++){
    fLineSMid[sm] = new TLine( sm+1, 0., sm+1, 91.);
  }
  //initialize all TRM counters to -1 i.e. invalid value
  // for (Int_t trm=0;trm<720;trm++){
  //   fTRMNoisyArray[trm]=-1;
  //   fTRMHwOkArray[trm]=-1;
  //   fTRMEnabledArray[trm]=-1;
  // }
  
}

//____________________________________________________________________________ 
AliTOFQADataMakerRec::AliTOFQADataMakerRec(const AliTOFQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fCalibData(qadm.fCalibData),
  fEnableNoiseFiltering(qadm.fEnableNoiseFiltering),
  fEnableDqmShifterOpt(qadm.fEnableDqmShifterOpt),
  fProcessedRawEventN(qadm.fProcessedRawEventN),
  fIsSOC(qadm.fIsSOC),
  fLineExpTimeMin(qadm.fLineExpTimeMin),
  fLineExpTimeMax(qadm.fLineExpTimeMax),
  fLineExpTotMin(qadm.fLineExpTotMin),
  fLineExpTotMax(qadm.fLineExpTotMax),
  fTOFRawStream(qadm.fTOFRawStream),
  fDecoderSummary(qadm.fDecoderSummary)
{
  //
  //copy ctor 
  //
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
   
  for (Int_t sm=0;sm<17;sm++){
    fLineSMid[sm]=qadm.fLineSMid[sm];
  }
  
  // for (Int_t trm=0;trm<10;trm++){
  //   fTRMNoisyArray[trm]=qadm.fTRMNoisyArray[trm];
  //   fTRMHwOkArray[trm]=qadm.fTRMHwOkArray[trm];
  //   fTRMEnabledArray[trm]=qadm.fTRMEnabledArray[trm];
  // }
  
}

//__________________________________________________________________
AliTOFQADataMakerRec& AliTOFQADataMakerRec::operator = (const AliTOFQADataMakerRec& qadm )
{
  //
  // assignment operator.
  //
  this->~AliTOFQADataMakerRec();
  new(this) AliTOFQADataMakerRec(qadm);
  return *this;
}
 
//----------------------------------------------------------------------------
AliTOFQADataMakerRec::~AliTOFQADataMakerRec()
{

  fTOFRawStream.Clear();
  delete fLineExpTimeMin;
  delete fLineExpTimeMax;
  delete fLineExpTotMin;
  delete fLineExpTotMax;
  for (Int_t sm=0;sm<10;sm++){
    delete fLineSMid[sm];
  }
  if (fDecoderSummary)
    delete fDecoderSummary;
}
//----------------------------------------------------------------------------
AliTOFChannelOnlineStatusArray* AliTOFQADataMakerRec::GetCalibData() const
{
  //
  // Retrive TOF calib objects from OCDB
  //
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *cdbe=0;
  
  cdbe = man->Get("TOF/Calib/Status",fRun);
  if(!cdbe){
    AliWarning("Load of calibration data from default storage failed!");
    AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    cdbe = man->Get("TOF/Calib/Status",fRun);
  }
  // Retrieval of data in directory TOF/Calib/Data:
  
  AliTOFChannelOnlineStatusArray * array = 0;
  if (cdbe) array = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();
  if (!array)  AliFatal("No calibration data from calibration database !");
  
  return array;
}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::InitRaws()
{
  //
  // create Raws histograms in Raws subdir
  //
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1I * h0 =  new TH1I("hTOFRaws",      "TOF raw hit multiplicity; TOF raw hits number; Events ",200, 0, 200) ;

  TH1I * h1 =  new TH1I("hTOFRawsIA",    "TOF raw hit multiplicity - I/A side; TOF raw hits number ;Events ",200, 0, 200) ;
  TH1I * h2 =  new TH1I("hTOFRawsOA",    "TOF raw hit multiplicity - O/A side; TOF raw hits number ;Events ",200, 0, 200) ;
  TH1I * h3 =  new TH1I("hTOFRawsIC",    "TOF raw hit multiplicity - I/C side; TOF raw hits number ;Events ",200, 0, 200) ;
  TH1I * h4 =  new TH1I("hTOFRawsOC",    "TOF raw hit multiplicity - O/C side; TOF raw hits number ;Events ",200, 0, 200) ;

  TH1F * h5  = new TH1F("hTOFRawsTime", "TOF Raws - Hit time (ns);Measured Hit time [ns];Hits", 250,0. ,610.) ; 
  TH1F * h6  = new TH1F("hTOFRawsTimeIA", "TOF Raws - Hit time (ns) - I/A side;Measured Hit time [ns];Hits", 250,0. ,610.) ; 
  TH1F * h7  = new TH1F("hTOFRawsTimeOA", "TOF Raws - Hit time (ns) - O/A side;Measured Hit time [ns];Hits", 250,0. ,610.) ; 
  TH1F * h8  = new TH1F("hTOFRawsTimeIC", "TOF Raws - Hit time (ns) - I/C side;Measured Hit time [ns];Hits", 250,0. ,610.) ; 
  TH1F * h9  = new TH1F("hTOFRawsTimeOC", "TOF Raws - Hit time (ns) - O/C side;Measured Hit time [ns];Hits", 250,0. ,610.) ; 
  
  TH1F * h10  = new TH1F("hTOFRawsToT", "TOF Raws - Hit ToT (ns);Measured Hit ToT (ns);Hits", 100, 0., 48.8) ; 
  
  TH1F * h11  = new TH1F("hTOFRawsToTIA", "TOF Raws - Hit ToT (ns) - I/A side;Measured Hit ToT (ns);Hits", 100, 0., 48.8) ; 
  TH1F * h12  = new TH1F("hTOFRawsToTOA", "TOF Raws - Hit ToT (ns) - O/A side;Measured Hit ToT (ns);Hits", 100, 0., 48.8) ; 
  TH1F * h13  = new TH1F("hTOFRawsToTIC", "TOF Raws - Hit ToT (ns) - I/C side;Measured Hit ToT (ns);Hits", 100, 0., 48.8) ; 
  TH1F * h14  = new TH1F("hTOFRawsToTOC", "TOF Raws - Hit ToT (ns) - O/C side;Measured Hit ToT (ns);Hits", 100, 0., 48.8) ; 
  
  TH1F * h15 = new TH1F("hTOFRawsLTMHits", "LTMs OR signals; Crate; Counts",  72, 0., 72.);
  TH2F * h16  = new TH2F("hTOFrefMap", "TOF enabled channel reference map;sector;strip",  72, 0., 18., 91, 0., 91.);
  TH2F * h17  = new TH2F("hTOFRawHitMap","TOF raw hit map;sector;strip", 72, 0., 18., 91, 0., 91.);
  TH2I * h18 = new TH2I("hTOFDecodingErrors","Decoding error monitoring; DDL; Error ", 72, 0, 72, 13,1,14);
  
  h18->GetYaxis()->SetBinLabel(1,"DRM ");
  h18->GetYaxis()->SetBinLabel(2,"LTM ");
  h18->GetYaxis()->SetBinLabel(3,"TRM 3 ");
  h18->GetYaxis()->SetBinLabel(4,"TRM 4");
  h18->GetYaxis()->SetBinLabel(5,"TRM 5");
  h18->GetYaxis()->SetBinLabel(6,"TRM 6");
  h18->GetYaxis()->SetBinLabel(7,"TRM 7");
  h18->GetYaxis()->SetBinLabel(8,"TRM 8");
  h18->GetYaxis()->SetBinLabel(9,"TRM 9");
  h18->GetYaxis()->SetBinLabel(10,"TRM 10");
  h18->GetYaxis()->SetBinLabel(11,"TRM 11");
  h18->GetYaxis()->SetBinLabel(12,"TRM 12");
  h18->GetYaxis()->SetBinLabel(13,"recovered");
  
  TH1F * h19  = new TH1F("hTOFOrphansTime", "TOF Raws - Orphans time (ns);Measured Hit time [ns];Hits", 250, 0. ,610.) ; 
  TH2F * h20 = new TH2F("hTOFRawTimeVsTRM035", "TOF raws - Hit time vs TRM - crates 0 to 35; TRM index = DDL*10+TRM(0-9);TOF raw time [ns]", 361, 0., 361., 250, 0., 610.0) ;
  TH2F * h21 = new TH2F("hTOFRawTimeVsTRM3671", "TOF raws - Hit time vs TRM - crates 36 to 72; TRM index = DDL**10+TRM(0-9);TOF raw time [ns]", 361, 360., 721., 250, 0., 610.0) ;
  TH2F * h22 = new TH2F("hTOFTimeVsStrip","TOF Raws - Hit time vs. strip (theta); Strip index;Raws TOF time (ns) ", 91,0.,91, 250, 0., 610.) ; 
  h0->Sumw2() ;
  h1->Sumw2() ;
  h2->Sumw2() ;
  h3->Sumw2() ;
  h4->Sumw2() ;
  //  h5->Sumw2() ;
  h6->Sumw2() ;
  h7->Sumw2() ;
  h8->Sumw2() ;
  h9->Sumw2() ;
  // h10->Sumw2() ;
  h11->Sumw2() ;
  h12->Sumw2() ;
  h13->Sumw2() ;
  h14->Sumw2() ;
  h15->Sumw2() ;
  h16->Sumw2() ;
  h17->Sumw2() ;
  h18->Sumw2() ;
  h19->Sumw2() ;
  h20->Sumw2() ;
  h21->Sumw2() ;
  h22->Sumw2() ;
  
  //add lines for DQM shifter
  //fLineExpTimeMin = new TLine(200., 0., 200., 0.);
  fLineExpTimeMin->SetLineColor(kGreen);
  fLineExpTimeMin->SetLineWidth(2);
  
  //fLineExpTimeMax = new TLine(250., 0., 250., 0.);
  fLineExpTimeMax->SetLineColor(kGreen);
  fLineExpTimeMax->SetLineWidth(2);
  
  //fLineExpTotMin = new TLine( 5., 0., 5., 0.);
  fLineExpTotMin->SetLineColor(kGreen);
  fLineExpTotMin->SetLineWidth(2);
  
  //fLineExpTotMax = new TLine(20., 0., 20., 0.);
  fLineExpTotMax->SetLineColor(kGreen);
  fLineExpTotMax->SetLineWidth(2);
  
  for (Int_t sm=0;sm<17;sm++){
    //fLineSMid[sm] = new TLine( 1+sm, 0., 1+sm, 91.);
    fLineSMid[sm]->SetLineColor(kMagenta);
    fLineSMid[sm]->SetLineWidth(2);
  }
  
  h5->GetListOfFunctions()->Add(fLineExpTimeMin);
  h5->GetListOfFunctions()->Add(fLineExpTimeMax);
  h10->GetListOfFunctions()->Add(fLineExpTotMin);
  h10->GetListOfFunctions()->Add(fLineExpTotMax);
  
  for (Int_t sm=0;sm<17;sm++){
    h16->GetListOfFunctions()->Add(fLineSMid[sm]);
    h17->GetListOfFunctions()->Add(fLineSMid[sm]);
  }
  
  
  TPaveText *phosHoleBox=new TPaveText(13,38,16,53,"b");	
  phosHoleBox->SetFillStyle(0);
  phosHoleBox->SetFillColor(kWhite);
  phosHoleBox->SetLineColor(kMagenta);
  phosHoleBox->SetLineWidth(2);
  phosHoleBox->AddText("PHOS");	
  h16->GetListOfFunctions()->Add(phosHoleBox);
  h17->GetListOfFunctions()->Add(phosHoleBox);

  // h0->SetDrawOption("logy");
  // h5->SetDrawOption("logy");
  // h10->SetDrawOption("logy");

  Add2RawsList(h0,   0, !expert,  image, !saveCorr) ;
  Add2RawsList(h1,   1,  expert,  !image, !saveCorr) ;
  Add2RawsList(h2,   2,  expert,  !image, !saveCorr) ;
  Add2RawsList(h3,   3,  expert,  !image, !saveCorr) ;
  Add2RawsList(h4,   4,  expert,  !image, !saveCorr) ;
  Add2RawsList(h5,   5, !expert,  image, !saveCorr) ;
  Add2RawsList(h6,   6,  expert,  !image, !saveCorr) ;
  Add2RawsList(h7,   7,  expert,  !image, !saveCorr) ;
  Add2RawsList(h8,   8,  expert,  !image, !saveCorr) ;
  Add2RawsList(h9,   9,  expert,  !image, !saveCorr) ;
  Add2RawsList(h10, 10, !expert,  image, !saveCorr) ;
  Add2RawsList(h11, 11,  expert, !image, !saveCorr) ;
  Add2RawsList(h12, 12,  expert, !image, !saveCorr) ;
  Add2RawsList(h13, 13,  expert, !image, !saveCorr) ;
  Add2RawsList(h14, 14,  expert, !image, !saveCorr) ;
  Add2RawsList(h15, 15, !expert,  image, !saveCorr) ;
  Add2RawsList(h16, 16, !expert,  image, !saveCorr) ;
  Add2RawsList(h17, 17, !expert,  image, !saveCorr) ;
  Add2RawsList(h18, 18,  expert, !image, !saveCorr) ;
  Add2RawsList(h19, 19,  expert, !image, !saveCorr) ;
  Add2RawsList(h20, 20,  expert, !image, !saveCorr) ;
  Add2RawsList(h21, 21,  expert, !image, !saveCorr) ;
  Add2RawsList(h22, 22,  expert, !image, !saveCorr) ;
}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::InitRecPoints()
{
  //
  // create RecPoints histograms in RecPoints subdir
  //
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1I * h0 = new TH1I("hTOFRecPoints",    "TOF RecPoints multiplicity ; TOF RecPoints number;Events",200, 0, 200) ;

  TH1F * h1 = new TH1F("hTOFRecPointsTimeIA", "RecPoints Time Spectrum in TOF (ns)- I/A side; Calibrated TOF time [ns];Events", 250, 0., 610.) ; 
  TH1F * h2 = new TH1F("hTOFRecPointsTimeOA", "RecPoints Time Spectrum in TOF (ns)- O/A side; Calibrated TOF time [ns];Events", 250, 0., 610.) ;
  TH1F * h3 = new TH1F("hTOFRecPointsTimeIC", "RecPoints Time Spectrum in TOF (ns)- I/C side; Calibrated TOF time [ns];Events", 250, 0., 610.) ;
  TH1F * h4 = new TH1F("hTOFRecPointsTimeOC", "RecPoints Time Spectrum in TOF (ns)- O/C side; Calibrated TOF time [ns];Events", 250, 0., 610.) ;
  
  TH1F * h5  = new TH1F("hTOFRecPointsRawTimeIA", "RecPoints raw Time Spectrum in TOF (ns)-I/A side; Measured TOF time [ns];Hits", 250, 0., 610.) ; 
  TH1F * h6  = new TH1F("hTOFRecPointsRawTimeOA", "RecPoints raw Time Spectrum in TOF (ns)-O/A side; Measured TOF time [ns];Hits", 250, 0., 610.) ; 
  TH1F * h7  = new TH1F("hTOFRecPointsRawTimeIC", "RecPoints raw Time Spectrum in TOF (ns)-I/C side; Measured TOF time [ns];Hits", 250, 0., 610.) ; 
  TH1F * h8  = new TH1F("hTOFRecPointsRawTimeOC", "RecPoints raw Time Spectrum in TOF (ns)-O/C side; Measured TOF time [ns];Hits", 250, 0., 610.) ; 
 
  TH1F * h9   = new TH1F("hTOFRecPointsToTIA", "RecPoints ToT Spectrum in TOF (ns)-I/A side; Measured TOT [ns];Hits", 100, 0., 48.8 ) ; 
  TH1F * h10  = new TH1F("hTOFRecPointsToTOA", "RecPoints ToT Spectrum in TOF (ns)-O/A side; Measured TOT [ns];Hits", 100, 0., 48.8 ) ; 
  TH1F * h11  = new TH1F("hTOFRecPointsToTIC", "RecPoints ToT Spectrum in TOF (ns)-I/C side; Measured TOT [ns];Hits", 100, 0., 48.8 ) ; 
  TH1F * h12  = new TH1F("hTOFRecPointsToTOC", "RecPoints ToT Spectrum in TOF (ns)-O/C side; Measured TOT [ns];Hits", 100, 0., 48.8 ) ; 
  
  TH2F * h13 = new TH2F("hTOFRecPointsClusMap","RecPoints map; sector ;strip", 72, 0., 18., 91, 0., 91.) ; 
  TH2F * h14 = new TH2F("hTOFRecPointsTimeVsStrip","RecPoints TOF time vs. strip (theta); Strip index; RecPoints TOF time (ns) ",91, 0., 91., 250, 0., 610.) ;
  TH2F * h15 = new TH2F("hTOFRecPointsTimeVsTRM035","TOF RecPoints time vs TRM - crates 0 to 35; TRM index = DDL*10+TRM(0-9);TOF time [ns]", 361, 0., 361., 250, 0., 610.0) ;
  TH2F * h16 = new TH2F("hTOFRecPointsTimeVsTRM3671","TOF RecPoints time vs TRM - crates 36 to 72; TRM index = DDL**10+TRM(0-9);TOF time [ns]", 361, 360., 721., 250, 0., 610.0) ;

  h0->Sumw2() ;
  h1->Sumw2() ;
  h2->Sumw2() ;
  h3->Sumw2() ;
  h4->Sumw2() ;
  h5->Sumw2() ;
  h6->Sumw2() ;
  h7->Sumw2() ;
  h8->Sumw2() ;
  h9->Sumw2() ;
  h10->Sumw2() ;
  h11->Sumw2() ;
  h12->Sumw2() ;
  h13->Sumw2() ;
  h14->Sumw2() ;
  h15->Sumw2() ;
  h16->Sumw2() ;
   
  Add2RecPointsList(h0, 0,   !expert,  image) ;
  Add2RecPointsList(h1, 1,   !expert,  image) ;
  Add2RecPointsList(h2, 2,   !expert,  image) ;
  Add2RecPointsList(h3, 3,   !expert,  image) ;
  Add2RecPointsList(h4, 4,   !expert,  image) ;
  Add2RecPointsList(h5, 5,    expert,  !image) ;
  Add2RecPointsList(h6, 6,    expert, !image) ;
  Add2RecPointsList(h7, 7,    expert, !image) ;
  Add2RecPointsList(h8, 8,    expert, !image) ;
  Add2RecPointsList(h9, 9,   !expert, !image) ;
  Add2RecPointsList(h10, 10, !expert, !image) ;
  Add2RecPointsList(h11, 11, !expert, !image) ;
  Add2RecPointsList(h12, 12, !expert, !image) ;
  Add2RecPointsList(h13, 13,  expert, !image) ;
  Add2RecPointsList(h14, 14,  expert, image) ;
  Add2RecPointsList(h15, 15,  expert, !image) ;
  Add2RecPointsList(h16, 16,  expert, !image) ;

}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::InitESDs()
{
  //
  //create ESDs histograms in ESDs subdir
  //

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1I * h0  = new TH1I("hTOFESDs", "Number of matched TOF tracks per event;Number of TOF matched ESD tracks;Counts", 200, 0, 200) ;  
  TH1F * h1  = new TH1F("hTOFESDsTime", "Matched  ESDs tracks: TOF Time spectrum; Calibrated TOF time [ns];Counts", 250, 0., 610. ) ; 
  TH1F * h2  = new TH1F("hTOFESDsRawTime", "Matched ESDs tracks: TOF raw Time spectrum;Measured TOF time [ns];Counts", 250, 0., 610.) ; 
  TH1F * h3  = new TH1F("hTOFESDsToT", "Matched ESDs tracks: TOF ToT spectrum; ESDs ToT [ns];Counts",100, 0., 48.8) ; 
  TH1F * h4  = new TH1F("hTOFESDskTOFOUT", "p_{T}  distribution of tracks with kTOFout; p_{T} (GeV/c);Counts", 50, 0.20, 5.00) ;  
  TH1F * h5  = new TH1F("hTOFESDskTIME", "p_{T}  distribution of tracks with kTOFout && kTIME; p_{T} (GeV/c);Counts", 50, 0.20, 5.00) ;  
  TH1F * h6  = new TH1F("hTOFESDsMatched", "p_{T} distribution of tracks with kTOFout && TOFtime>0; p_{T} (GeV/c);Counts", 50, 0.20, 5.00) ;  
  TH1F * h7  = new TH1F("hTOFESDsMatchingProb", "TPC-TOF track-matching probability;TOF matching probability (%)  ;Counts",101, -1.0, 100) ;  
  TH1F * h8  = new TH1F("hTOFESDsDiffTime", "ESDs t_{TOF}-t_{exp,pi} spectrum in TOF (ps); t_{TOF}-t_{exp,pi} [ps];Counts", 200, -2440., 2440.) ; 
  TH1F * h9  = new TH1F("hTOFHitsLength", "Matched ESDs tracks: Length Spectrum; Track length [cm];Counts", 800, 0., 800) ; 

  h0->Sumw2() ;
  h1->Sumw2() ;
  h2->Sumw2() ;
  h3->Sumw2() ;
  h4->Sumw2() ;
  h5->Sumw2() ;
  h6->Sumw2() ;
  h7->Sumw2() ;
  h8->Sumw2() ;
  h9->Sumw2() ;

  Add2ESDsList(h0, 0, !expert,  image) ;
  Add2ESDsList(h1, 1, !expert,  image) ;
  Add2ESDsList(h2, 2,  expert,  !image) ;
  Add2ESDsList(h3, 3, !expert,  !image) ;
  Add2ESDsList(h4, 4,  expert,  image) ;
  Add2ESDsList(h5, 5,  expert,  image) ;
  Add2ESDsList(h6, 6,  expert,  image) ; 
  Add2ESDsList(h7, 7,  expert,  image) ; 
  Add2ESDsList(h8, 8,  expert,  !image) ; 
  Add2ESDsList(h9, 9, !expert,  !image) ;
 
}


//____________________________________________________________________________
void AliTOFQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //
  // makes data from Raws
  //
  if (rawReader->GetType()==7) {
   
    Double_t tdc2ns=AliTOFGeometry::TdcBinWidth()*1E-3;//in ns
    Double_t tot2ns=AliTOFGeometry::ToTBinWidth()*1E-3;
    Int_t ntof[5]; /* 0=tot, 1=IA, 2=OA, 3=IC, 4=OC*/
    for (Int_t j=0;j<5;j++){ ntof[j]=0;}
    Int_t equipmentID[5]; //(ddl, trm, chain,tdc,channel)
    Int_t volumeID[5];   //(sector,plate,strip,padX,padZ)
    Int_t volumeID2[5];   //(sector,plate,strip,padZ,padX) to use AliTOFGeometry::GetIndex()
    Int_t chIndex=-1;
    
    TClonesArray * clonesRawData;
    fTOFRawStream.SetRawReader(rawReader);
    
    //uncomment if needed to apply DeltaBC correction
    //fTOFRawStream.ApplyBCCorrections(kTRUE);
    
    fDecoderSummary->Reset();
    for (Int_t iDDL = 0; iDDL < AliTOFGeometry::NDDL()*AliTOFGeometry::NSectors(); iDDL++){
      rawReader->Reset();
      fTOFRawStream.LoadRawDataBuffersV2(iDDL);
      
      //* get decoding error counters
      fDecoderSummary = ( (AliTOFDecoderV2*) fTOFRawStream.GetDecoderV2() )->GetDecoderSummaryData();
      if ( (fDecoderSummary) && (fDecoderSummary ->GetErrorDetected()) ) {
	Int_t errorSlotID=(Int_t) fDecoderSummary->GetErrorSlotID();
	GetRawsData(18)->Fill(iDDL,errorSlotID);
	if (fDecoderSummary -> GetRecoverError() ) 		
	  GetRawsData(18)->Fill(iDDL,13);
      }     
      
      clonesRawData = (TClonesArray*)fTOFRawStream.GetRawData();
      for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {
	AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);
	
	if (tofRawDatum->GetTOF()){
	  equipmentID[0]=iDDL;
	  equipmentID[1]=tofRawDatum->GetTRM(); 
	  equipmentID[2]=tofRawDatum->GetTRMchain();
	  equipmentID[3]=tofRawDatum->GetTDC();
	  equipmentID[4]=tofRawDatum->GetTDCchannel();
	  
	  if (CheckEquipID(equipmentID)){
	    fTOFRawStream.EquipmentId2VolumeId(iDDL, 
					       tofRawDatum->GetTRM(), 
					       tofRawDatum->GetTRMchain(),
					       tofRawDatum->GetTDC(), 
					       tofRawDatum->GetTDCchannel(), 
					       volumeID);
	    //LTM data
	    if (FilterLTMData(equipmentID)) { //counts LTM hits
	      if (equipmentID[2]==1)  { //crate left, A-side or C-side
		GetRawsData(15)->Fill(equipmentID[0]);
	      } else {
		if (equipmentID[0]<36) { GetRawsData(15)->Fill(equipmentID[0]-1); }
		else  { GetRawsData(15)->Fill(equipmentID[0]+1); }
	      }
	      continue;
	    }
	    
	    //TRM data
	    if (CheckVolumeID(volumeID)){  
	      volumeID2[0]=volumeID[0];
	      volumeID2[1]=volumeID[1];
	      volumeID2[2]=volumeID[2];
	      volumeID2[3]=volumeID[4];
	      volumeID2[4]=volumeID[3];
	      chIndex=AliTOFGeometry::GetIndex(volumeID2);
	      
	      if (tofRawDatum->GetTOT()){	    
	       	if (!(fCalibData->GetNoiseStatus(chIndex)==AliTOFChannelOnlineStatusArray::kTOFNoiseBad)
		    && (fCalibData->GetHWStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFHWOk)) {//noise and enabled filter
		  ntof[0]++; //counter for tof hits
		  
		  //fill global spectra for DQM plots
		  GetRawsData(5)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;//in ns
		  GetRawsData(10)->Fill( tofRawDatum->GetTOT()*tot2ns) ;//in ns
		  
		  //fill side-related spectra for experts plots
		  Int_t ddlACside=iDDL/36; // 0 or 1
		  Int_t ddlPerSm=iDDL%4;
		  
		  if (volumeID2[0]>4 && volumeID2[0]<14){       //O side
		    if (ddlPerSm<2){ //A side
		      ntof[1]++;
		      GetRawsData(6)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;
		      GetRawsData(11)->Fill( tofRawDatum->GetTOT()*tot2ns) ;
		    } else {  //C side
		      ntof[3]++;
		      GetRawsData(8)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;
		      GetRawsData(13)->Fill( tofRawDatum->GetTOT()*tot2ns) ;
		    }
		  } else {                                    
		    if (volumeID2[0]<5 || volumeID2[0]>13){   //I side
		      if (ddlPerSm<2){ //A side
			ntof[2]++;
			GetRawsData(7)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;
			GetRawsData(12)->Fill( tofRawDatum->GetTOT()*tot2ns) ;
		      } else {//C side
			ntof[4]++;
			GetRawsData(9)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;
			GetRawsData(14)->Fill( tofRawDatum->GetTOT()*tot2ns) ;
		      }
		    }	
		  }
		  
		  //compute TRM offset
		  Int_t trm= iDDL*10+(equipmentID[1]-3);
		  GetRawsData(20+ddlACside)->Fill(trm,tofRawDatum->GetTOF()*tdc2ns);
		  GetRawsData(22)->Fill(GetStripIndex(volumeID),tofRawDatum->GetTOF()*tdc2ns) ;
		  Short_t fea = volumeID2[4]/12;
		  Float_t hitmapx = volumeID2[0] + ((Double_t)(3 - fea) + 0.5) *0.25;
		  GetRawsData(17)->Fill(hitmapx,GetStripIndex(volumeID2));
		}//noise filter
	      }//end hit selection
	      else { //orphans
		if (!(fCalibData->GetNoiseStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad)
		    && (fCalibData->GetHWStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFHWOk))
		  GetRawsData(19)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;//in ns
	      }//end orphans
	    }//end volumeID check
	  }//end equipID check
	}//end tof check
      }//loop on raw data
      clonesRawData->Clear();
    } // DDL Loop
    
    for (Int_t j=0;j<5;j++) { GetRawsData(j)->Fill(ntof[j]); }
    fProcessedRawEventN++;
    fTOFRawStream.Clear();
  } else {
    AliDebug(1,Form("Event of type %d found. Skipping non-physics event for QA.\n", rawReader->GetType())); 
  }
  
  //fill reference map for DQM shifter only once in a detector cycle 
  if (fIsSOC) {
    Int_t geoId[5]={-1,-1,-1,-1,-1};// pgeoId=(sm, mod, strip, padZ, padX)
    Int_t detId[5]={-1,-1,-1,-1,-1};//detID=(ddl,trm,tdc, chain,channel)
    Int_t trmIndex=-1;
    for (Int_t ch = 0; ch <  fCalibData->GetSize(); ch++) {
      AliTOFGeometry::GetVolumeIndices(ch,geoId);
      AliTOFRawStream::Geant2EquipmentId(geoId,detId); 
      if ((detId[1]<0)||(detId[0]<0)) continue;
      trmIndex=(detId[1]-3)+detId[0]*10;
      
      //set TRM counters
      // if (fCalibData->GetNoiseStatus(ch)==AliTOFChannelOnlineStatusArray::kTOFNoiseBad)
      // 	fTRMNoisyArray[trmIndex]+=1;
      // if (fCalibData->GetHWStatus(ch) == AliTOFChannelOnlineStatusArray::kTOFHWOk)
      // 	fTRMHwOkArray[trmIndex]+=1;
      
      if ( (!(fCalibData->GetNoiseStatus(ch)==AliTOFChannelOnlineStatusArray::kTOFNoiseBad))
      	   && (fCalibData->GetHWStatus(ch) == AliTOFChannelOnlineStatusArray::kTOFHWOk) ){
	//fTRMEnabledArray[trmIndex]+=1;	
	//fill reference map with info from OCDB
	Short_t fea = geoId[4]/12;
	Float_t hitmapx = geoId[0] + ((Double_t)(3 - fea) + 0.5)*0.25;
	GetRawsData(16)->Fill(hitmapx, GetStripIndex(geoId));
      }
    }
    //printf("Counters for noisy, enabled and good channels in TOF  TRMs read from OCDB.\n");
    fIsSOC=kFALSE;
  }
    
  //enable options for DQM shifter
  EnableDqmShifterOpt(kTRUE);
}

//____________________________________________________________________________
void AliTOFQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
    //
  // Make data from Clusters
  //
 
  Double_t tdc2ns=AliTOFGeometry::TdcBinWidth()*1E-3;
  Double_t tot2ns=AliTOFGeometry::ToTBinWidth()*1E-3;
  
  Int_t volumeID2[5];//(sm, plate,strip, padZ,padX)
  //Int_t volumeID[5];//(sm, plate,strip, padX,padZ) 
    
  TBranch *branch=clustersTree->GetBranch("TOF");
  if (!branch) { 
    AliError("can't get the branch with the TOF clusters !");
    return;
  }

  static TClonesArray dummy("AliTOFcluster",10000);
  dummy.Clear();
  TClonesArray *clusters=&dummy;
  branch->SetAddress(&clusters);
  
  // Import the tree
  clustersTree->GetEvent(0);  
 
  GetRecPointsData(0)->Fill((Int_t)clusters->GetEntriesFast()) ; 
  
  TIter next(clusters) ; 
  AliTOFcluster * c ; 
  while ( (c = dynamic_cast<AliTOFcluster *>(next())) ) {
   
      // volumeID2[0] = c->GetDetInd(0);
      // volumeID2[1] = c->GetDetInd(1);
      // volumeID2[2] = c->GetDetInd(2);
      // volumeID2[3] = c->GetDetInd(4); //padX
      // volumeID2[4] = c->GetDetInd(3); //padZ 
      
      for (Int_t i=0;i<5;i++){
	volumeID2[i]=c->GetDetInd(i); //X and Z indeces inverted in RecPoints
      }
      //Int_t chIndex=AliTOFGeometry::GetIndex(volumeID2);
      Int_t iDDL=AliTOFRawStream::Geant2DDL(volumeID2);
      Int_t iTRM=AliTOFRawStream::Geant2TRM(volumeID2);
      Short_t fea = volumeID2[4]/12;
      Float_t hitmapx = volumeID2[0] + ((Double_t)(3 - fea) + 0.5) *0.25;
      Int_t ddlACside=iDDL/36; // 0 or 1
      Int_t ddlPerSm=iDDL%4;
      
      if ((c->GetTDCRAW()) && (c->GetTDC()) && (c->GetToT())){
	if (volumeID2[0]>4 && volumeID2[0]<14){       //I side
	  if (ddlPerSm<2){ //A side
	    GetRecPointsData(1)->Fill( c->GetTDC()*tdc2ns) ;//in ns
	    GetRecPointsData(5)->Fill( c->GetTDCRAW()*tdc2ns) ;//in ns
	    GetRecPointsData(9)->Fill( c->GetToT()*tot2ns) ;//in ns
	  } else {//C side
	    GetRecPointsData(3)->Fill( c->GetTDC()*tdc2ns) ;//in ns
	    GetRecPointsData(7)->Fill( c->GetTDCRAW()*tdc2ns) ;//in ns
	    GetRecPointsData(11)->Fill( c->GetToT()*tot2ns) ;//in ns
	  }
	} else {
	  if (volumeID2[0]<5 || volumeID2[0]>13){       //O side
	    if (ddlPerSm<2){ //A side
	      GetRecPointsData(2)->Fill( c->GetTDC()*tdc2ns) ;//in ns
	      GetRecPointsData(6)->Fill( c->GetTDCRAW()*tdc2ns) ;//in ns
	      GetRecPointsData(10)->Fill( c->GetToT()*tot2ns) ;//in ns
	    } else { //C side
	      GetRecPointsData(4)->Fill( c->GetTDC()*tdc2ns) ;//in ns
	      GetRecPointsData(8)->Fill( c->GetTDCRAW()*tdc2ns) ;//in ns
	      GetRecPointsData(12)->Fill( c->GetToT()*tot2ns) ;//in ns
	    }
	  }
	}
	GetRecPointsData(13)->Fill(hitmapx,GetStripIndex(volumeID2));
	GetRecPointsData(14)->Fill(GetStripIndex(volumeID2), c->GetTDC()*tdc2ns) ;
	Int_t trm= iDDL*10+(iTRM-3);
	if (ddlACside==0) { //A side
	  GetRecPointsData(15)->Fill(trm,c->GetTDC()*tdc2ns);
	} else {//C side
	  GetRecPointsData(16)->Fill(trm,c->GetTDC()*tdc2ns);
	}
      }//hit selection
  }//end while   
  EnableDqmShifterOpt(kFALSE);
}

//____________________________________________________________________________
void AliTOFQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //
  // make QA data from ESDs
  //  
  const Double_t speedOfLight = TMath::C()*1E2*1E-12; // cm/ps
  const Double_t pionMass = 0.13957018; //GeV/c^2

  Int_t ntrk = esd->GetNumberOfTracks() ; 
  Int_t ntpc=0;
  Int_t ntofout=0;
  
    while (ntrk--) {
      AliESDtrack *track=esd->GetTrack(ntrk);
      Double_t tofTime=track->GetTOFsignal();//in ps
      Double_t tofTimeRaw=track->GetTOFsignalRaw();//in ps
      Double_t tofToT=track->GetTOFsignalToT(); //in ps
      
      UInt_t status=track->GetStatus();
      if (track->IsOn(AliESDtrack::kTPCrefit)) {
	ntpc++;
	Double_t y=track->Eta();
	if (TMath::Abs(y)<0.9) { //select TOF acceptance
	  if ((status&AliESDtrack::kTOFout)!=0)  { //define matching
	    ntofout++;
	    GetESDsData(1)->Fill(tofTime*1E-3);
	    GetESDsData(2)->Fill(tofTimeRaw*1E-3); 
	    GetESDsData(3)->Fill(tofToT*1E-3);
	    GetESDsData(4)->Fill(track->Pt());
	    
	    Double_t length =track->GetIntegratedLength();
	    Double_t mom2=(track->Pt()*track->Pt())+(track->Pz()*track->Pz());
	    Double_t piTexp = TMath::Sqrt(1+(pionMass*pionMass/mom2))*length/speedOfLight; //in ps
	    GetESDsData(8)->Fill(tofTime-piTexp);
	    GetESDsData(9)->Fill(length);
	    
	    if ((status&AliESDtrack::kTIME)!=0) 
	      GetESDsData(5)->Fill(track->Pt());
	    
	    if (tofTime>0)
	      GetESDsData(6)->Fill(track->Pt());
	  } //end check on matched tracks
	} 
      }//end check on TPCrefit
    }
    
    GetESDsData(0)->Fill(ntofout) ;
    if(ntpc>0){
      Float_t ratio = (Float_t)ntofout/(Float_t)ntpc*100.; //matching probability
      GetESDsData(7)->Fill(ratio) ;
    }
    
    if(ntofout>0) {
	Float_t ratio = (Float_t)ntofout/(Float_t)ntpc*100; //matched over propagated to TOF outer radius
	GetESDsData(8)->Fill(ratio) ;
    }
    EnableDqmShifterOpt(kFALSE);
}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::StartOfDetectorCycle()
{
  //
  //Detector specific actions at start of cycle
  // ResetAllTRMcounters();
  fCalibData = GetCalibData();
  fIsSOC=kTRUE;
  return;
}  

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    
    if (fEnableDqmShifterOpt){
      // Help make the raw qa histogram easier to interpret for the DQM shifter
      if (!GetRawsData(0) || !GetRawsData(5) || !GetRawsData(10) 
	  || !GetRawsData(15) || !GetRawsData(16) || !GetRawsData(17)) {
	printf("No histogram for DQM found - Possible memory corruption ???. Please check\n") ; 
	continue;
      }
      printf("=========>Processed %i physics raw events. \n",fProcessedRawEventN);
      
      //Double_t monitorPeriodLength=fProcessedRawEventN*600*1E-9;//in s
      
      if (fCalibData){
	//set minima and maxima to allow log scale
	Double_t yTimeMax = GetRawsData(5)->GetMaximum()*1.05;
	Double_t yTotMax = GetRawsData(10)->GetMaximum()*1.05;
	fLineExpTimeMin->SetY2(yTimeMax);
	fLineExpTimeMax->SetY2(yTimeMax);
	fLineExpTotMin->SetY2(yTotMax);
	fLineExpTotMax->SetY2(yTotMax);
	
	for (Int_t j=0;j<18;j++){
	  if ((j==0)||(j==5)||(j==10)||(j==15)||(j==16)||(j==17)) {
	    GetRawsData(j)->GetXaxis()->SetLabelOffset(0.005);
	    GetRawsData(j)->GetXaxis()->SetLabelSize(0.05);
	    GetRawsData(j)->GetXaxis()->SetTitleOffset(0.8);
	    GetRawsData(j)->GetXaxis()->SetTitleSize(0.05);
	    GetRawsData(j)->GetYaxis()->SetLabelOffset(0.005);
	    GetRawsData(j)->GetYaxis()->SetLabelSize(0.06);
	    GetRawsData(j)->GetYaxis()->SetTitleOffset(0.8);
	    GetRawsData(j)->GetYaxis()->SetTitleSize(0.06);	  
	  }
	}
	//make up for all histos 
	for(Int_t j=0;j<5;j++){
	  GetRawsData(j)->SetMarkerColor(kBlue);
	  GetRawsData(j)->SetMarkerStyle(8);
	  GetRawsData(j)->SetMarkerSize(0.7);
	}
	for(Int_t j=5;j<15;j++){
	  GetRawsData(j)->SetLineColor(kBlue);
	  GetRawsData(j)->SetLineWidth(1);
	  GetRawsData(j)->SetMarkerColor(kBlue);
	  //GetRawsData(j)->SetFillColor(kWhite);
	  //GetRawsData(j)->SetDrawOption("bar");
	}
	
	GetRawsData(15)->SetLineColor(kBlue);
	GetRawsData(15)->SetLineWidth(1);
	GetRawsData(15)->SetMarkerStyle(8);
	GetRawsData(15)->SetMarkerSize(0.7);
	GetRawsData(15)->SetMarkerColor(kBlue);//Option("bar");
	      
	GetRawsData(16)->SetOption("colz");
	GetRawsData(17)->SetOption("colz");
	GetRawsData(18)->SetOption("colz"); 
      }
    }//END ENABLE DQM SHIFTER OPT
  } //end for
  AliQAChecker::Instance()->Run(AliQAv1::kTOF, task, list) ;  
}
//____________________________________________________________________________
void AliTOFQADataMakerRec::GetMapIndeces(const Int_t* const in , Int_t* out)
{
  //
  //return appropriate indeces for the theta-phi map
  //

  Int_t npadX = AliTOFGeometry::NpadX();
  Int_t npadZ = AliTOFGeometry::NpadZ();
  Int_t nStripA = AliTOFGeometry::NStripA();
  Int_t nStripB = AliTOFGeometry::NStripB();
  Int_t nStripC = AliTOFGeometry::NStripC();

  Int_t isector = in[0];
  Int_t iplate = in[1];
  Int_t istrip = in[2];
  Int_t ipadX = in[3]; 
  Int_t ipadZ = in[4]; 
  
  Int_t stripOffset = 0;
  switch (iplate) {
  case 0:
    stripOffset = 0;
      break;
  case 1:
    stripOffset = nStripC;
    break;
  case 2:
    stripOffset = nStripC+nStripB;
    break;
  case 3:
    stripOffset = nStripC+nStripB+nStripA;
    break;
  case 4:
    stripOffset = nStripC+nStripB+nStripA+nStripB;
    break;
  default:
    AliDebug(1,Form("Wrong plate number in TOF (%d) !",iplate));
    break;
  };
  Int_t zindex=npadZ*(istrip+stripOffset)+(ipadZ+1);
  Int_t phiindex=npadX*isector+ipadX+1;
  out[0]=zindex;  
  out[1]=phiindex;  
  
}

//---------------------------------------------------------------
Int_t AliTOFQADataMakerRec::GetStripIndex(const Int_t * const in)
{
    /* return tof strip index between 0 and 91 */

  Int_t nStripA = AliTOFGeometry::NStripA();
  Int_t nStripB = AliTOFGeometry::NStripB();
  Int_t nStripC = AliTOFGeometry::NStripC();

  // Int_t isector = in[0];
  Int_t iplate = in[1];
  Int_t istrip = in[2];
  //Int_t ipadX = in[3]; 
  //Int_t ipadZ = in[4]; 
  
  Int_t stripOffset = 0;
  switch (iplate) {
  case 0:
    stripOffset = 0;
      break;
  case 1:
    stripOffset = nStripC;
    break;
  case 2:
    stripOffset = nStripC+nStripB;
    break;
  case 3:
    stripOffset = nStripC+nStripB+nStripA;
    break;
  case 4:
    stripOffset = nStripC+nStripB+nStripA+nStripB;
    break;
  default:
      AliDebug(1,Form("Wrong plate number in TOF (%d) !",iplate));
      stripOffset=-1;
      break;
  };
  
  if (stripOffset<0 || stripOffset>92) return -1;
  else 
      return (stripOffset+istrip);
  
}
//---------------------------------------------------------------
Bool_t  AliTOFQADataMakerRec::CheckVolumeID(const Int_t * const volumeID)
{
    //
    //Checks volume ID validity
    //
   
    for (Int_t j=0;j<5;j++){
	if (volumeID[j]<0) {
	    AliDebug(1,Form("Invalid detector volume index for volumeID[%i]",j));
	    return kFALSE;
	}
    }
    return kTRUE;
    
}

//---------------------------------------------------------------
Bool_t  AliTOFQADataMakerRec::CheckEquipID(const Int_t * const equipmentID)
{
    //
    //Checks equipment ID validity
    
   for (Int_t j=0;j<5;j++){
	if (equipmentID[j]<0) {
	  AliDebug(1,Form("Invalid equipment volume index for equipmentID[%i]",j));
	  return kFALSE;
	}
   }
   return kTRUE;
}
//---------------------------------------------------------------
Bool_t  AliTOFQADataMakerRec::FilterLTMData(const Int_t * const equipmentID) const
{
  /*It returns kTRUE if data come from LTM.
    It thus filters trigger-related signals  */

  Int_t ddl, trm, tdc;
  //if (!CheckEquipID(equipmentID)) return kFALSE;
  ddl = equipmentID[0];
  trm = equipmentID[1];
  tdc = equipmentID[3];
  
  if ((ddl%2==1) && (trm==3) && (tdc>11 && tdc<15))
    return kTRUE;
  else 
    return kFALSE;
 
}
//---------------------------------------------------------------
Bool_t  AliTOFQADataMakerRec::FilterSpare(const Int_t * const equipmentID) const
{
  /*It returns kTRUE if data come from spare 
    equipment ID. 
    So far only check on TRM 3 crate left is implemented */

  Int_t ddl, trm, tdc;
  //if (!CheckEquipID(equipmentID)) return kFALSE;
  ddl = equipmentID[0];
  trm = equipmentID[1];
  tdc = equipmentID[3];
  
  if ((ddl%2==1) && (trm==3) && (tdc>2 && tdc<12))
    return kTRUE;
  else 
    return kFALSE;
 
}
//----------------------------------------------------------------
/*void  AliTOFQADataMakerRec::ResetAllTRMcounters()
{
  //resets all TRM counters to 0
  for (Int_t trm=0;trm<720;trm++){
    fTRMNoisyArray[trm]=-1;
    fTRMHwOkArray[trm]=-1;
    fTRMEnabledArray[trm]=-1;
  } 
  return;
  
}
*/
