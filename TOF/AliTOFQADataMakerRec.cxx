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
Modified by fbellini on 01/11/2011
- added histograms for LTM monitoring
- fix for coverity

Modified by fbellini on 17/10/2011
- fix for memory leak in constructor
- added methods to read histos ranges from config file in DQM
- added CTTM maps + relative methods to retrieve CTTM numbering
- removed obslete comments

Modified by fbellini & rshanoian on 06/07/2011
- changes for trigger classes implementation
- fRunNumber added as private member
- added time vs BCID plot

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
#include <iostream>
#include <fstream>

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

Int_t AliTOFQADataMakerRec::fgNbinsMultiplicity=200; //number of bins in multiplicity plot
Int_t AliTOFQADataMakerRec::fgRangeMinMultiplicity=0;//min range in multiplicity plot
Int_t AliTOFQADataMakerRec::fgRangeMaxMultiplicity=200;//max range in multiplicity plot
Int_t AliTOFQADataMakerRec::fgNbinsTime=250;//number of bins in time plot
const Float_t AliTOFQADataMakerRec::fgkNbinsWidthTime=2.44;//width of bins in time plot
Float_t AliTOFQADataMakerRec::fgRangeMinTime=0.0;//range min in time plot
Float_t AliTOFQADataMakerRec::fgRangeMaxTime=620.0; //range max in time plot
Int_t AliTOFQADataMakerRec::fgCutNmaxFiredMacropad=5;//cut on number of max fired macropad
const Int_t AliTOFQADataMakerRec::fgkFiredMacropadLimit=50;//cut on number of max fired macropad


//____________________________________________________________________________ 
  AliTOFQADataMakerRec::AliTOFQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kTOF), "TOF Quality Assurance Data Maker"),
  fCalibData(0x0),
  fEnableNoiseFiltering(kFALSE),
  fEnableDqmShifterOpt(kFALSE),
  fIsSOC(kFALSE),
  fLineExpTimeMin(0x0),
  fLineExpTimeMax(0x0),
  fLineExpTotMin(0x0),
  fLineExpTotMax(0x0),
  fTOFRawStream(AliTOFRawStream()),
  fDecoderSummary(0),
  fRunNumber(-1),
  fCalib(AliTOFcalib())
{
  //
  // ctor
  //   
  // fLineExpTimeMin = new TLine(200., 0., 200., 0.);
  // fLineExpTimeMax = new TLine(250., 0., 250., 0.);
  // fLineExpTotMin = new TLine(5., 0., 5., 0.);
  // fLineExpTotMax = new TLine(20., 0., 20., 0.);
  for (Int_t sm=0;sm<17;sm++){
    fLineSMid[sm] = new TLine( sm+1, 0., sm+1, 91.);
  }

  for (Int_t sm=0;sm<71;sm++){
    fLineLTMid[sm] = new TLine( sm+1, 0., sm+1, 23.);
  }

  for (Int_t sm=0;sm<22;sm++){
    fLineLTMbitId[sm] = new TLine( 0., sm+1, 72. ,sm+1);
  }
  
}

//____________________________________________________________________________ 
AliTOFQADataMakerRec::AliTOFQADataMakerRec(const AliTOFQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fCalibData(qadm.fCalibData),
  fEnableNoiseFiltering(qadm.fEnableNoiseFiltering),
  fEnableDqmShifterOpt(qadm.fEnableDqmShifterOpt),
  fIsSOC(qadm.fIsSOC),
  fLineExpTimeMin(qadm.fLineExpTimeMin),
  fLineExpTimeMax(qadm.fLineExpTimeMax),
  fLineExpTotMin(qadm.fLineExpTotMin),
  fLineExpTotMax(qadm.fLineExpTotMax),
  fTOFRawStream(qadm.fTOFRawStream),
  fDecoderSummary(qadm.fDecoderSummary),
  fRunNumber(qadm.fRunNumber),
  fCalib(qadm.fCalib)
{
  //
  //copy ctor 
  //
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
   
  for (Int_t sm=0;sm<17;sm++){
    fLineSMid[sm]=qadm.fLineSMid[sm];
  }

 for (Int_t sm=0;sm<71;sm++){
    fLineLTMid[sm] = qadm.fLineLTMid[sm];
  }

  for (Int_t sm=0;sm<22;sm++){
    fLineLTMbitId[sm] = qadm.fLineLTMbitId[sm];
  }
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
  //destructor
  fTOFRawStream.Clear();
  fCalib.Clear();
  if (fLineExpTimeMin)
    delete fLineExpTimeMin;
  if (fLineExpTimeMax)
    delete fLineExpTimeMax;
  if (fLineExpTotMin)
    delete fLineExpTotMin;
  if (fLineExpTotMax)
    delete fLineExpTotMax;
  for (Int_t sm=0;sm<17;sm++){
    if (fLineSMid[sm])
      delete fLineSMid[sm];
  }
  for (Int_t sm=0;sm<71;sm++){
    if (fLineLTMid[sm])
      delete fLineLTMid[sm];
  }
for (Int_t sm=0;sm<22;sm++){
    if (fLineLTMbitId[sm])
      delete fLineLTMbitId[sm];
  }
}
//----------------------------------------------------------------------------
AliTOFChannelOnlineStatusArray* AliTOFQADataMakerRec::GetCalibData() 
{
  //
  // Retrive TOF calib objects from OCDB
  //
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *cdbe=0;
 
  if (fRun<=0) fRunNumber=145288; //reference run from LHC11a
  else fRunNumber=fRun;
  
  if (man->GetRun()!=fRunNumber){
    fRunNumber=man->GetRun();
    AliWarning(Form("Run number mismatch found: setting it to value from current AliCDBManager instance = %i", fRunNumber));
  }
  cdbe = man->Get("TOF/Calib/Status",fRunNumber);
  
  if(!cdbe){
    // for DQM online
    AliWarning("Load of calibration data from default (alien://) storage failed!");
    printf("Calibration data will be loaded from local storage - ok if on DQM station!");
    man->SetDefaultStorage("local:///local/cdb/");
    cdbe = man->Get("TOF/Calib/Status",fRun);
    
    if(!cdbe){
      AliWarning("Load of calibration data from local DQM machine storage failed!");
      AliWarning("Calibration data will be loaded from local ($ALICE_ROOT) storage ");
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
      cdbe = man->Get("TOF/Calib/Status",fRunNumber);
    }
  }
  // Retrieval of data in directory TOF/Calib/Data:
  AliTOFChannelOnlineStatusArray * array = 0;
  if (cdbe) {
    printf("======= OCDB object for TOF retrieved from run %i in %s\n",fRunNumber,cdbe->GetName());
    array = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();
  }
  if (!array)  AliFatal("No calibration data from calibration database !");
  
  fCalib.Init(fRunNumber);
  return array;
}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::InitRaws()
{
  //
  // create Raws histograms in Raws subdir
  //
  ReadHistogramRangeFromFile(gSystem->Getenv("TOFDQMHISTO_CONFIGFILE"));
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1I * h0 =  new TH1I("hTOFRaws","TOF raw hit multiplicity; TOF raw hits number; Events ",fgNbinsMultiplicity, fgRangeMinMultiplicity, fgRangeMaxMultiplicity);
  TH1I * h1 =  new TH1I("hTOFRawsIA","TOF raw hit multiplicity - I/A side; TOF raw hits number;Events ",fgNbinsMultiplicity, fgRangeMinMultiplicity, fgRangeMaxMultiplicity);
  TH1I * h2 =  new TH1I("hTOFRawsOA","TOF raw hit multiplicity - O/A side; TOF raw hits number;Events ",fgNbinsMultiplicity, fgRangeMinMultiplicity, fgRangeMaxMultiplicity);
  TH1I * h3 =  new TH1I("hTOFRawsIC","TOF raw hit multiplicity - I/C side; TOF raw hits number;Events ",fgNbinsMultiplicity, fgRangeMinMultiplicity, fgRangeMaxMultiplicity);
  TH1I * h4 =  new TH1I("hTOFRawsOC","TOF raw hit multiplicity - O/C side; TOF raw hits number;Events ",fgNbinsMultiplicity, fgRangeMinMultiplicity, fgRangeMaxMultiplicity);

  TH1F * h5  = new TH1F("hTOFRawsTime", "TOF Raws - Hit time (ns);Measured Hit time [ns];Hits", fgNbinsTime,fgRangeMinTime,fgRangeMaxTime); 
  TH1F * h6  = new TH1F("hTOFRawsTimeIA", "TOF Raws - Hit time (ns) - I/A side;Measured Hit time [ns];Hits", fgNbinsTime,fgRangeMinTime,fgRangeMaxTime); 
  TH1F * h7  = new TH1F("hTOFRawsTimeOA", "TOF Raws - Hit time (ns) - O/A side;Measured Hit time [ns];Hits", fgNbinsTime,fgRangeMinTime,fgRangeMaxTime); 
  TH1F * h8  = new TH1F("hTOFRawsTimeIC", "TOF Raws - Hit time (ns) - I/C side;Measured Hit time [ns];Hits", fgNbinsTime,fgRangeMinTime,fgRangeMaxTime); 
  TH1F * h9  = new TH1F("hTOFRawsTimeOC", "TOF Raws - Hit time (ns) - O/C side;Measured Hit time [ns];Hits", fgNbinsTime,fgRangeMinTime,fgRangeMaxTime); 
  
  TH1F * h10  = new TH1F("hTOFRawsToT", "TOF Raws - Hit ToT (ns);Measured Hit ToT (ns);Hits", 100, 0., 48.8); 
  
  TH1F * h11  = new TH1F("hTOFRawsToTIA", "TOF Raws - Hit ToT (ns) - I/A side;Measured Hit ToT (ns);Hits", 100, 0., 48.8); 
  TH1F * h12  = new TH1F("hTOFRawsToTOA", "TOF Raws - Hit ToT (ns) - O/A side;Measured Hit ToT (ns);Hits", 100, 0., 48.8); 
  TH1F * h13  = new TH1F("hTOFRawsToTIC", "TOF Raws - Hit ToT (ns) - I/C side;Measured Hit ToT (ns);Hits", 100, 0., 48.8); 
  TH1F * h14  = new TH1F("hTOFRawsToTOC", "TOF Raws - Hit ToT (ns) - O/C side;Measured Hit ToT (ns);Hits", 100, 0., 48.8); 
  
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
  
  TH1F * h19  = new TH1F("hTOFOrphansTime", "TOF Raws - Orphans time (ns);Measured Hit time [ns];Hits",fgNbinsTime,fgRangeMinTime,fgRangeMaxTime); 
  TH2F * h20 = new TH2F("hTOFRawTimeVsTRM035", "TOF raws - Hit time vs TRM - crates 0 to 35; TRM index = DDL*10+TRM(0-9);TOF raw time [ns]", 361, 0., 361.,fgNbinsTime,fgRangeMinTime,fgRangeMaxTime);
  TH2F * h21 = new TH2F("hTOFRawTimeVsTRM3671", "TOF raws - Hit time vs TRM - crates 36 to 72; TRM index = DDL**10+TRM(0-9);TOF raw time [ns]", 361, 360., 721.,fgNbinsTime, fgRangeMinTime,fgRangeMaxTime);
  TH2F * h22 = new TH2F("hTOFTimeVsStrip","TOF raw hit time vs. MRPC (along z axis); MRPC index along z axis; Raws TOF time (ns) ", 91,0.,91,fgNbinsTime,fgRangeMinTime,fgRangeMaxTime); 
  TH2F * h23 = new TH2F("hTOFtimeVsBCID","TOF time vs BCID; BCID; time (ns) ", 3564, 0., 3564.,fgNbinsTime,fgRangeMinTime,fgRangeMaxTime);
  TH2F * h24 = new TH2F("hTOFchannelEfficiencyMap","TOF channels (HWok && efficient && !noisy && !problematic);sector;strip",  72, 0., 18., 91, 0., 91.);
  TH2F * h25 = new TH2F("hTOFhitsCTTM","Map of hit pads according to CTTM numbering;LTM index;bit index",  72, 0., 72., 23, 0., 23.);
  TH2F * h26 = new TH2F("hTOFmacropadCTTM","Map of hit macropads according to CTTM numbering;LTM index; bit index",  72, 0., 72., 23, 0., 23.);
  TH2F * h27 = new TH2F("hTOFmacropadDeltaPhiTime","#Deltat vs #Delta#Phi of hit macropads;#Delta#Phi (degrees);#DeltaBX",  18, 0., 180., 20, 0., 20.0);
  TH2I *h28 = new TH2I("hBXVsCttmBit","BX ID in TOF matching window vs trg channel; trg channel; BX", 1728, 0, 1728, 24, 0, 24); 
  TH2F *h29 = new TH2F("hTimeVsCttmBit","TOF raw time vs trg channel; trg channel; raw time (ns)", 1728, 0., 1728., fgNbinsTime, fgRangeMinTime, fgRangeMaxTime); 

  h25->GetYaxis()->SetTickLength(-0.02);
  h26->GetYaxis()->SetTickLength(-0.02);
  h25->GetYaxis()->SetNdivisions(210);
  h26->GetYaxis()->SetNdivisions(210);
  h25->GetXaxis()->SetTickLength(-0.02);
  h26->GetXaxis()->SetTickLength(-0.02);
  h25->GetXaxis()->SetLabelOffset(0.015);
  h26->GetXaxis()->SetLabelOffset(0.015);
  h25->GetXaxis()->SetNdivisions(515);
  h26->GetXaxis()->SetNdivisions(515);
  
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
  h23->Sumw2() ;
  h24->Sumw2() ;
  h25->Sumw2() ;
  h26->Sumw2() ;
  h27->Sumw2() ;
  h28->Sumw2() ;
  h29->Sumw2() ;

  //add lines for DQM shifter
  fLineExpTimeMin = new TLine(200., 0., 200., 0.);
  fLineExpTimeMax = new TLine(250., 0., 250., 0.);
  fLineExpTotMin = new TLine(5., 0., 5., 0.);
  fLineExpTotMax = new TLine(20., 0., 20., 0.);

  fLineExpTimeMin->SetLineColor(kGreen);
  fLineExpTimeMin->SetLineWidth(2);
  
  fLineExpTimeMax->SetLineColor(kGreen);
  fLineExpTimeMax->SetLineWidth(2);
  
  fLineExpTotMin->SetLineColor(kGreen);
  fLineExpTotMin->SetLineWidth(2);
  
  fLineExpTotMax->SetLineColor(kGreen);
  fLineExpTotMax->SetLineWidth(2);
  
  for (Int_t sm=0;sm<17;sm++){
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
  
  for (Int_t sm=0;sm<71;sm++){
    fLineLTMid[sm]->SetLineColor(kBlack);
    fLineLTMid[sm]->SetLineWidth(1);
    h26->GetListOfFunctions()->Add(fLineLTMid[sm]);
    h25->GetListOfFunctions()->Add(fLineLTMid[sm]);
  }
  for (Int_t sm=0;sm<22;sm++){
    fLineLTMbitId[sm]->SetLineColor(kBlack);
    fLineLTMbitId[sm]->SetLineWidth(1);
    h26->GetListOfFunctions()->Add(fLineLTMbitId[sm]);
    h25->GetListOfFunctions()->Add(fLineLTMbitId[sm]);
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
  Add2RawsList(h15, 15,  expert, !image, !saveCorr) ;
  Add2RawsList(h16, 16,  !expert,  image, !saveCorr) ;
  Add2RawsList(h17, 17,  !expert,  image, !saveCorr) ;
  Add2RawsList(h18, 18,   expert, !image, !saveCorr) ;
  Add2RawsList(h19, 19,   expert, !image, !saveCorr) ;
  Add2RawsList(h20, 20,   expert, !image, !saveCorr) ;
  Add2RawsList(h21, 21,   expert, !image, !saveCorr) ;
  Add2RawsList(h22, 22,  !expert,  image, !saveCorr) ;
  Add2RawsList(h23, 23,  !expert, !image, !saveCorr) ;
  Add2RawsList(h24, 24,  !expert, !image, !saveCorr) ;
  Add2RawsList(h25, 25,  !expert, !image, !saveCorr) ;
  Add2RawsList(h26, 26,  !expert,  image, !saveCorr) ;
  Add2RawsList(h27, 27,  !expert,  image, !saveCorr) ;
  Add2RawsList(h28, 28,  !expert,  !image, !saveCorr) ;
  Add2RawsList(h29, 29,  !expert,  !image, !saveCorr) ;
  
//
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
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
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
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
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line 
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
    Int_t indexCTTM[2]={-1,-1};	    
    Int_t indexGeo2CTTM[2]={-1,-1};
    Float_t macropadPhiTimeUPC[fgkFiredMacropadLimit][2];
    for (Int_t ii=0;ii<2;ii++){
      for (Int_t jj=0;jj<fgkFiredMacropadLimit;jj++){	
	macropadPhiTimeUPC[jj][ii]=-999.0; 
      }
    }

    TClonesArray * clonesRawData;
    fTOFRawStream.SetRawReader(rawReader);
    Int_t BCID=rawReader->GetBCID();
    
    Int_t nFiredMacropad=0,
      iFiredMacropad=-1;
    nFiredMacropad=GetNumberOfFiredMacropad(rawReader);
    
    //uncomment if needed to apply DeltaBC correction
    //fTOFRawStream.ApplyBCCorrections(kTRUE);
    
    if (fDecoderSummary){
      fDecoderSummary->Reset();
    }
    for (Int_t iDDL = 0; iDDL < AliTOFGeometry::NDDL()*AliTOFGeometry::NSectors(); iDDL++){
      rawReader->Reset();
      fTOFRawStream.LoadRawDataBuffersV2(iDDL);
      
      //* get decoding error counters
      fDecoderSummary = ( (AliTOFDecoderV2*) fTOFRawStream.GetDecoderV2() )->GetDecoderSummaryData();
      if ( (fDecoderSummary) && (fDecoderSummary ->GetErrorDetected()) ) {
	Int_t errorSlotID=(Int_t) fDecoderSummary->GetErrorSlotID();
	FillRawsData(18,iDDL,errorSlotID);
	if (fDecoderSummary -> GetRecoverError() ) 		
	  FillRawsData(18,iDDL,13);
      }     
      
      clonesRawData = (TClonesArray*)fTOFRawStream.GetRawData();
      for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {
	AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);
	Float_t tofRawTime=tofRawDatum->GetTOF()*tdc2ns;
	
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
	      if (equipmentID[2]==1)  //crate left, A-side or C-side
		FillRawsData(15,equipmentID[0]);
	      else 
		FillRawsData(15,equipmentID[0]-1); 
	      
	      //retrieve CTTM index (Ltm, bit)
	      GetCTTMIndex(equipmentID, indexCTTM);
	      
	      //get BX index within TOF-matching window
	      Int_t indexBC=-1;
	      indexBC= TMath::Nint(tofRawTime/24.4);

	      Int_t indexCttmChannel=indexCTTM[0]*24+indexCTTM[1];
	      FillRawsData(28,indexCttmChannel,indexBC);
	      FillRawsData(29,indexCttmChannel,tofRawTime);
	      
	      //fired macropad map (from LTM hits) - only for low multi evts (UPC)
	      if ((nFiredMacropad<=fgCutNmaxFiredMacropad)){	
		iFiredMacropad++;
		AliInfo(Form("Event found with %i fired macropads in BCID = %i!",nFiredMacropad,BCID));
		FillRawsData(26,indexCTTM[0],indexCTTM[1]);
		Float_t halfSMphi=-999.0;
		if (indexCTTM[0]<36)
		  halfSMphi=indexCTTM[0]*10.+5.;
		else  halfSMphi=(indexCTTM[0]-36)*10.+5.;
		macropadPhiTimeUPC[iFiredMacropad][0]=halfSMphi;
		macropadPhiTimeUPC[iFiredMacropad][1]=indexBC;
	      }
	    }
	    
	    //TRM data
	    if (CheckVolumeID(volumeID)){  
	      volumeID2[0]=volumeID[0];
	      volumeID2[1]=volumeID[1];
	      volumeID2[2]=volumeID[2];
	      volumeID2[3]=volumeID[4];
	      volumeID2[4]=volumeID[3];
	      chIndex=AliTOFGeometry::GetIndex(volumeID2);
	      
	      //fill hit map according to CTTM numbering
	      GetGeo2CTTMIndex(volumeID2, indexGeo2CTTM);
	      if ((nFiredMacropad<=fgCutNmaxFiredMacropad)){	
		FillRawsData(25,indexGeo2CTTM[0],indexGeo2CTTM[1]);
	      }
	      //hits selection
	      if (tofRawDatum->GetTOT()){	    
	       	if (!(fCalibData->GetNoiseStatus(chIndex)==AliTOFChannelOnlineStatusArray::kTOFNoiseBad)
		    && (fCalibData->GetHWStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFHWOk)) {//noise and enabled filter
		  ntof[0]++; //counter for tof hits
		  
		  //fill global spectra for DQM plots
		  FillRawsData(5, tofRawTime) ;//in ns
		  FillRawsData(10, tofRawDatum->GetTOT()*tot2ns) ;//in ns
		  FillRawsData(23, BCID, tofRawTime) ;//in ns
		  
		  //fill side-related spectra for experts plots
		  Int_t ddlACside=iDDL/36; // 0 or 1
		  Int_t ddlPerSm=iDDL%4;
		  
		  if (volumeID2[0]>4 && volumeID2[0]<14){       //O side
		    if (ddlPerSm<2){ //A side
		      ntof[1]++;
		      FillRawsData(6, tofRawTime) ;
		      FillRawsData(11, tofRawDatum->GetTOT()*tot2ns) ;
		    } else {  //C side
		      ntof[3]++;
		      FillRawsData(8, tofRawTime) ;
		      FillRawsData(13, tofRawDatum->GetTOT()*tot2ns) ;
		    }
		  } else {                                    
		    if (volumeID2[0]<5 || volumeID2[0]>13){   //I side
		      if (ddlPerSm<2){ //A side
			ntof[2]++;
			FillRawsData(7, tofRawTime) ;
			FillRawsData(12, tofRawDatum->GetTOT()*tot2ns) ;
		      } else {//C side
			ntof[4]++;
			FillRawsData(9, tofRawTime) ;
			FillRawsData(14, tofRawDatum->GetTOT()*tot2ns) ;
		      }
		    }	
		  }
		  
		  //compute TRM offset
		  Int_t trm= iDDL*10+(equipmentID[1]-3);
		  FillRawsData(20+ddlACside,trm,tofRawTime);
		  FillRawsData(22,GetStripIndex(volumeID),tofRawTime) ;
		  Short_t fea = volumeID2[4]/12;
		  Float_t hitmapx = volumeID2[0] + ((Double_t)(3 - fea) + 0.5) *0.25;
		  FillRawsData(17,hitmapx,GetStripIndex(volumeID2));
		  if (fCalib.IsChannelEnabled(chIndex,kTRUE,kTRUE))//checks also if efficient and if problematic
		    FillRawsData(24,hitmapx,GetStripIndex(volumeID2));
		}//noise filter
	      }//end hit selection
	      else { //orphans
		if (!(fCalibData->GetNoiseStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad)
		    && (fCalibData->GetHWStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFHWOk))
		  FillRawsData(19, tofRawTime) ;//in ns
	      }//end orphans
	    }//end volumeID check
	  }//end equipID check
	}//end tof check
      }//loop on raw data
      clonesRawData->Clear();
    } // DDL Loop
    
    for (Int_t j=0;j<5;j++) FillRawsData(j,ntof[j]);
    fTOFRawStream.Clear();
  
    if ((nFiredMacropad<=fgCutNmaxFiredMacropad)){
      Float_t deltaPhiMacropad=-999.;
      Float_t deltaTimeMacropad=-999.;
      for (Int_t j=0;j<fgCutNmaxFiredMacropad+1; j++){
	for (Int_t k=j+1;k<fgCutNmaxFiredMacropad+1; k++){
	  if ((macropadPhiTimeUPC[j][0]>0.0)&&(macropadPhiTimeUPC[k][0]>0.0)){
	    deltaPhiMacropad=TMath::Abs(macropadPhiTimeUPC[j][0]-macropadPhiTimeUPC[k][0]);
	    deltaTimeMacropad=TMath::Abs(macropadPhiTimeUPC[j][1]-macropadPhiTimeUPC[k][1]);
	    if (deltaPhiMacropad<=180.)
	      FillRawsData(27, deltaPhiMacropad,deltaTimeMacropad);
	    else
	      FillRawsData(27, TMath::Abs(360.0-deltaPhiMacropad),deltaTimeMacropad);
	  }
	} 
      }    
    }//end cut on number of fired macropad
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

      if ( (!(fCalibData->GetNoiseStatus(ch)==AliTOFChannelOnlineStatusArray::kTOFNoiseBad))
      	   && (fCalibData->GetHWStatus(ch) == AliTOFChannelOnlineStatusArray::kTOFHWOk) ){
	//fill reference map with info from OCDB
	Short_t fea = geoId[4]/12;
	Float_t hitmapx = geoId[0] + ((Double_t)(3 - fea) + 0.5)*0.25;
	FillRawsData(16,hitmapx, GetStripIndex(geoId));
      }
    }
    fIsSOC=kFALSE;
  }
  
  //enable options for DQM shifter
  EnableDqmShifterOpt(kTRUE);
  //
  IncEvCountCycleRaws();
  IncEvCountTotalRaws();
  //
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
 
  FillRecPointsData(0,(Int_t)clusters->GetEntriesFast()) ; 
  
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
	    FillRecPointsData(1, c->GetTDC()*tdc2ns) ;//in ns
	    FillRecPointsData(5, c->GetTDCRAW()*tdc2ns) ;//in ns
	    FillRecPointsData(9, c->GetToT()*tot2ns) ;//in ns
	  } else {//C side
	    FillRecPointsData(3, c->GetTDC()*tdc2ns) ;//in ns
	    FillRecPointsData(7, c->GetTDCRAW()*tdc2ns) ;//in ns
	    FillRecPointsData(11, c->GetToT()*tot2ns) ;//in ns
	  }
	} else {
	  if (volumeID2[0]<5 || volumeID2[0]>13){       //O side
	    if (ddlPerSm<2){ //A side
	      FillRecPointsData(2, c->GetTDC()*tdc2ns) ;//in ns
	      FillRecPointsData(6, c->GetTDCRAW()*tdc2ns) ;//in ns
	      FillRecPointsData(10, c->GetToT()*tot2ns) ;//in ns
	    } else { //C side
	      FillRecPointsData(4, c->GetTDC()*tdc2ns) ;//in ns
	      FillRecPointsData(8, c->GetTDCRAW()*tdc2ns) ;//in ns
	      FillRecPointsData(12, c->GetToT()*tot2ns) ;//in ns
	    }
	  }
	}
	FillRecPointsData(13,hitmapx,GetStripIndex(volumeID2));
	FillRecPointsData(14,GetStripIndex(volumeID2), c->GetTDC()*tdc2ns) ;
	Int_t trm= iDDL*10+(iTRM-3);
	if (ddlACside==0) { //A side
	  FillRecPointsData(15,trm,c->GetTDC()*tdc2ns);
	} else {//C side
	  FillRecPointsData(16,trm,c->GetTDC()*tdc2ns);
	}
      }//hit selection
  }//end while   
  EnableDqmShifterOpt(kFALSE);
  //
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //
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
	    FillESDsData(1,tofTime*1E-3);
	    FillESDsData(2,tofTimeRaw*1E-3); 
	    FillESDsData(3,tofToT*1E-3);
	    FillESDsData(4,track->Pt());
	    
	    Double_t length =track->GetIntegratedLength();
	    Double_t mom2=(track->Pt()*track->Pt())+(track->Pz()*track->Pz());
	    Double_t piTexp = TMath::Sqrt(1+(pionMass*pionMass/mom2))*length/speedOfLight; //in ps
	    FillESDsData(8,tofTime-piTexp);
	    FillESDsData(9,length);
	    
	    if ((status&AliESDtrack::kTIME)!=0) 
	      FillESDsData(5,track->Pt());
	    
	    if (tofTime>0)
	      FillESDsData(6,track->Pt());
	  } //end check on matched tracks
	} 
      }//end check on TPCrefit
    }
    
    FillESDsData(0,ntofout) ;
    if(ntpc>0){
      Float_t ratio = (Float_t)ntofout/(Float_t)ntpc*100.; //matching probability
      FillESDsData(7,ratio) ;
    }
    
    if(ntofout>0) {
	Float_t ratio = (Float_t)ntofout/(Float_t)ntpc*100; //matched over propagated to TOF outer radius
	FillESDsData(8,ratio) ;
    }
    EnableDqmShifterOpt(kFALSE);
    //
    IncEvCountCycleESDs();
    IncEvCountTotalESDs();
    //
}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::StartOfDetectorCycle()
{
  //
  //Detector specific actions at start of cycle
  fCalibData = GetCalibData();
  fIsSOC=kTRUE;
  return;
}  

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  ResetEventTrigClasses();
  //
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) continue ;     
    SetEventSpecie(AliRecoParam::ConvertIndex(specie));  

    for (int itc=-1;itc<GetNTrigClasses();itc++) { // RS: loop over eventual clones per trigger class

      if (fEnableDqmShifterOpt) {
	// RS: fetch the histograms for given trigger class
	TObjArray& arrRW = *GetRawsDataOfTrigClass(itc);
	
	// Help make the raw qa histogram easier to interpret for the DQM shifter
	if (!arrRW[ 0] || !arrRW[ 5] || !arrRW[10] || !arrRW[15] || !arrRW[16] || !arrRW[17]) continue;
	
	printf("=========>Processed %i physics raw of specie %s with TrigGlass %d\n",
	       GetEvCountCycleRaws(itc),AliRecoParam::GetEventSpecieName(specie), itc);
	
	//Double_t monitorPeriodLength=fProcessedRawEventN*600*1E-9;//in s
      
	if (fCalibData){
	  //set minima and maxima to allow log scale
	  Double_t yTimeMax = ((TH1*)arrRW[5])->GetMaximum()*1.05;
	  Double_t yTotMax = ((TH1*)arrRW[10])->GetMaximum()*1.05;
	  fLineExpTimeMin->SetY2(yTimeMax);
	  fLineExpTimeMax->SetY2(yTimeMax);
	  fLineExpTotMin->SetY2(yTotMax);
	  fLineExpTotMax->SetY2(yTotMax);
	  //
	  for (Int_t j=0;j<18;j++){
	    if ((j==0)||(j==5)||(j==10)||(j==15)||(j==16)||(j==17)) {
	      TH1* htmp = (TH1*)arrRW[j];
	      htmp->GetXaxis()->SetLabelOffset(0.005);
	      htmp->GetXaxis()->SetLabelSize(0.05);
	      htmp->GetXaxis()->SetTitleOffset(0.8);
	      htmp->GetXaxis()->SetTitleSize(0.05);
	      htmp->GetYaxis()->SetLabelOffset(0.005);
	      htmp->GetYaxis()->SetLabelSize(0.06);
	      htmp->GetYaxis()->SetTitleOffset(0.8);
	      htmp->GetYaxis()->SetTitleSize(0.06);	  
	    }
	  }
	  //make up for all histos 
	  for(Int_t j=0;j<5;j++) {
	    TH1* htmp = (TH1*)arrRW[j];  
	    if (!htmp) continue;
	    htmp->SetMarkerColor(kBlue);
	    htmp->SetMarkerStyle(8);
	    htmp->SetMarkerSize(0.7);
	  }
	  for(Int_t j=5;j<15;j++) {
	    TH1* htmp = (TH1*)arrRW[j];
	    if (!htmp) continue;
	    htmp->SetLineColor(kBlue);
	    htmp->SetLineWidth(1);
	    htmp->SetMarkerColor(kBlue);
	  }
	  
	  TH1* htmp =  (TH1*)arrRW[15];  
	  htmp->SetLineColor(kBlue);
	  htmp->SetLineWidth(1);
	  htmp->SetMarkerStyle(8);
	  htmp->SetMarkerSize(0.7);
	  htmp->SetMarkerColor(kBlue);//Option("bar");
	  //
	  TString title25 = Form("Map of hit pads according to CTTM numbering (Max Fired Macropad = %i)",fgCutNmaxFiredMacropad);
	  TString title26 = Form("Map of hit macropads according to CTTM numbering (Max Fired Macropad = %i)",fgCutNmaxFiredMacropad);
	  TString title27 = Form("#Deltat vs #Delta#Phi of hit macropads (Max Fired Macropad = %i)",fgCutNmaxFiredMacropad);
	  
	  if ( (htmp=(TH1*)arrRW[16]) ) htmp->SetOption("colz");
	  if ( (htmp=(TH1*)arrRW[17]) ) htmp->SetOption("colz");
	  if ( (htmp=(TH1*)arrRW[18]) ) htmp->SetOption("colz"); 
	  if ( (htmp=(TH1*)arrRW[22]) ) htmp->SetOption("colz"); 
	  if ( (htmp=(TH1*)arrRW[23]) ) htmp->SetOption("colz"); 
	  if ( (htmp=(TH1*)arrRW[24]) ) htmp->SetOption("colz"); 
	  if ( (htmp=(TH1*)arrRW[28]) ) htmp->SetOption("colz"); 
	  if ( (htmp=(TH1*)arrRW[29]) ) htmp->SetOption("colz"); 

	  if ( (htmp=(TH1*)arrRW[25]) ) {
	    htmp->SetOption("colz"); 
	    htmp->SetTitle(title25.Data());
	  }
	  if ( (htmp=(TH1*)arrRW[26]) ) {
	    htmp->SetOption("colz"); 
	    htmp->SetTitle(title26.Data());
	  }
	  if ( (htmp=(TH1*)arrRW[27]) ){
	    htmp->SetOption("colz"); 
	    htmp->SetTitle(title27.Data());
	  }
	}
      }//END ENABLE DQM SHIFTER OPT
    } // RS: loop over trigger classes
  } //end for
  //
  AliQAChecker::Instance()->Run(AliQAv1::kTOF, task, list) ;  
  //
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

//-----------------------------------------------------------------------------
void AliTOFQADataMakerRec::GetGeo2LTMIndex(const Int_t * const detind, Int_t *indexLTM) {
  //
  // getting LTMmatrix indexes for current digit
  //
  Int_t stripId=GetStripIndex(detind);

  if (detind[1]==0 || detind[1]==1 || (detind[1]==2 && detind[2]<=7)) { //A side
    if (detind[4]<24){ //R
      indexLTM[0] = detind[0]*2;
    } else { //L
      indexLTM[0] = detind[0]*2+1;
    }  
    indexLTM[1]=stripId;

  } else { //C side
    if (detind[4]<24){
      indexLTM[0] = detind[0]*2+36;
    } else {
      indexLTM[0] = (detind[0]*2+1)+36;
    }
    indexLTM[1]=90-stripId; 
  }
  
  // if (indexLTM[0]<36) { //A side
  //   if (detind[1] ==0){
  //     indexLTM[1] = detind[2];
  //   }
  //   else if (detind[1] ==1){
  //     indexLTM[1] = detind[2]+nStripB;
  //   }
  //   else if (detind[1] ==2){
  //     indexLTM[1] = detind[2]+19*2;
  //   }
  //   else{
  //     AliError("Smth Wrong!!!");
  //   }
  // }
  // else { //C side
  //   if (detind[1]==2){
  //     if (detind[4]<24)
  // 	indexLTM[1] = (nStripAL-detind[2])+nStripC+nStripB;
  //     else 
  // 	indexLTM[1] = (nStripAR-detind[2])+nStripC+nStripB;
  //   }
  //   else if (detind[1] ==3){
  //     indexLTM[1] = (nStripB-detind[2])+nStripC;
  //   }
  //   else if (detind[1] ==4){
  //     indexLTM[1] = nStripC-detind[2];
  //   }
  //   else{
  //     AliError("Smth Wrong!!!");
  //   }
  // }  
}

//-----------------------------------------------------------------------------
void AliTOFQADataMakerRec::GetGeo2CTTMIndex(Int_t *detind, Int_t *indexCTTM) {
  //
  // Returns CTTM index corresponding to the detector element detind
  //
  GetGeo2LTMIndex(detind,indexCTTM);
  indexCTTM[1]/=2;
  return;
}

//-------------------------------------------------------------------------
Int_t AliTOFQADataMakerRec::GetNumberOfFiredMacropad(AliRawReader * rawReader){
  
  Int_t nFired=0;
  TClonesArray * clonesRawData;  
  if (!rawReader) return 0;  
  for (Int_t iDDL = 0; iDDL < AliTOFGeometry::NDDL()*AliTOFGeometry::NSectors(); iDDL++){
    rawReader->Reset();
    fTOFRawStream.LoadRawDataBuffersV2(iDDL); 
    clonesRawData = (TClonesArray*)fTOFRawStream.GetRawData();
    for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {
      AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);    
      if (tofRawDatum->GetTOF()){	
	if ( (tofRawDatum->GetTRM()==3)&&
	     (tofRawDatum->GetTDC()>11)&&
	     (tofRawDatum->GetTDC()<15)) {	  
	  nFired+=1;  
	}
      }
    }    
  }//loop over DDLs
  return nFired; 
}

//----------------------------------------------------------------
void AliTOFQADataMakerRec::GetCTTMIndex(Int_t *equipid, Int_t *indexCTTM) {
  //
  // Returns CTTM index corresponding to the equipment id equipid, only for LTM hits
  // equipid = (crate, trm, chain, tdc, channel)

  if ( (equipid[1]!=3)||(equipid[3]<12) ){
    indexCTTM[0]=-1;
    indexCTTM[1]=-1;
    return;
  }  
  Int_t modStrip2LTM[3][8]={ { 0, 1, 2, 3, 4, 5, 6, 7},
			     { 8, 9, 10, 11, 12, 13, 14, 15},
	                     {16, 17, 18, 19, 20, 21, 22, 23}
                            }; 

  Int_t DDL2LTMmatrix[72]={0,1,37,36,2,3,39,38,4,5,41,40,6,7,43,42,8,9,45,44,10,11,47,46,12,13,49,48,14,15,51,50,16,17,53,52,18,19,
			   55,54,20,21,57,56,22,23,59,58,24,25,61,60,26,27,63,62,28,29,65,64,30,31,67,66,32,33,69,68,34,35,71,70};

  Int_t itdc=equipid[3]%12;
  Int_t crate=-1;
  if (equipid[2]==0)
   crate=equipid[0]-1;
  else crate=equipid[0];
  
  indexCTTM[0]=DDL2LTMmatrix[crate];
  indexCTTM[1]=modStrip2LTM[itdc][equipid[4]];      
  return;
}

//_____________________________________________________________________________
void AliTOFQADataMakerRec::ReadHistogramRangeFromFile(const Char_t * filename)
{
  //
  // read histogram ranges from configuration file
  //
  if (!filename) {
    AliInfo("Config file with histograms ranges not found or invalid -> use default values.");
    SetDefaultHistogramRange();
    SetDefaultCutNmaxFiredMacropad();
    return;
  }
  
  std::fstream configFile;
  configFile.open(filename, std::fstream::in);
  if (!configFile.is_open()){
    AliInfo("Cannot open config file with histograms ranges -> use default values.");
    SetDefaultHistogramRange();
    return;
  }
  
  //check file size
  Int_t begin = configFile.tellg();
  configFile.seekg(0, std::fstream::end); /* end */
  Int_t end = configFile.tellg();
  Int_t size = end - begin;
  configFile.seekg(0, std::fstream::beg); /* rewind file */
  if (size <= 0){
    AliInfo(Form("Unexpected EOF of config file with histograms ranges. File size: %d -> use default values", size));
    SetDefaultHistogramRange();
    return;
  }
  
  Int_t minMulti=9999, maxMulti=-9999;
  Int_t nbinsMulti=0,nbinsTime=0;
  Float_t minTime=9999.0, maxTime=-9999.0;
  Int_t cutFiredMacropad=0;
  TString endoflist;
  while (!configFile.eof()) {
    configFile >> cutFiredMacropad >> minMulti >> maxMulti >> minTime >> maxTime;
    configFile >> endoflist;
    if (endoflist.Contains("end")) break;
  }

  //set multiplicity histo ranges
  if (minMulti>maxMulti){
    AliInfo("Invalid range for multiplicity histogram set. Changing to default values.");
    SetDefaultMultiHistogramRange();
  } else {
    nbinsMulti = maxMulti-minMulti;
    SetNbinsMultiplicityHisto(nbinsMulti);
    SetMultiplicityHistoRange(minMulti,maxMulti);
    //AliInfo(Form("Setting multiplicity histogram ranges to: multMin = %i - multMax = %i - nMultBins = %i", fgRangeMinMultiplicity, fgRangeMaxMultiplicity, fgNbinsMultiplicity));
  }

  //set time histo ranges
  if (minTime>maxTime){
    AliInfo("Invalid range for time histogram set. Changing to default values.");
    SetDefaultTimeHistogramRange();
  } else {
    nbinsTime = TMath::Nint((maxTime - minTime)/fgkNbinsWidthTime);//ns
    maxTime=minTime+nbinsTime*fgkNbinsWidthTime;//ns
    SetNbinsTimeHisto(nbinsTime);
    SetTimeHistoRange(minTime,maxTime);
    //AliInfo(Form("Setting time histogram ranges to: timeMin = %5.2f ns - timeMax = %5.2f ns - nTimeBins = %i", fgRangeMinTime, fgRangeMaxTime,fgNbinsTime));
  } 
 
  if ((cutFiredMacropad>0)&&(cutFiredMacropad<fgkFiredMacropadLimit)){
    AliInfo("Invalid value for cut on fired macropad. Changing to default values.");
    SetDefaultCutNmaxFiredMacropad();
  } else {
    SetCutNmaxFiredMacropad(cutFiredMacropad);
    //AliInfo(Form("Setting cut on fired macropad to:  = %i",cutFiredMacropad));
  } 
  AliInfo(Form("Setting: multMin = %i - multMax = %i - nMultBins = %i, timeMin = %5.2f ns - timeMax = %5.2f ns - nTimeBins = %i, cutMaxFiredMacropad = %i", 
	       fgRangeMinMultiplicity, fgRangeMaxMultiplicity, fgNbinsMultiplicity, fgRangeMinTime, fgRangeMaxTime,fgNbinsTime, cutFiredMacropad));
  configFile.close();
  return;
}

//_____________________________________________________________________________
void AliTOFQADataMakerRec::SetDefaultHistogramRange()
{
  //
  // set default histogram ranges (tuned on 2011 pp collisions)
  // 
  AliInfo("Setting all histogram ranges to default values.");
  SetDefaultMultiHistogramRange();
  SetDefaultTimeHistogramRange();
  SetDefaultCutNmaxFiredMacropad();
  return;
}

//_____________________________________________________________________________
void AliTOFQADataMakerRec::SetDefaultMultiHistogramRange()
{
  //
  // set default histogram ranges (tuned on 2011 pp collisions)
  // 
  SetMultiplicityHistoRange (0, 200);
  SetNbinsMultiplicityHisto(200);
  AliInfo("Setting Multiplicity histogram ranges to default values.");
  AliInfo(Form("multMin = %i - multMax = %i - nMultBins = %i",
	       fgRangeMinMultiplicity, fgRangeMaxMultiplicity, fgNbinsMultiplicity));
  return;
}

//_____________________________________________________________________________
void AliTOFQADataMakerRec::SetDefaultTimeHistogramRange()
{
  //
  // set default histogram ranges (tuned on 2011 pp collisions)
  // 
  SetNbinsTimeHisto(250);
  SetTimeHistoRange (0.0,610.);   
  
  AliInfo("Setting Time histogram ranges to default values:");
  AliInfo(Form("timeMin = %5.2f ns - timeMax = %5.2f ns - nTimeBins = %i",
	       fgRangeMinTime, fgRangeMaxTime,fgNbinsTime));
  return;
}

//------------------------------------------------------------------
void AliTOFQADataMakerRec::SetDefaultCutNmaxFiredMacropad()
{
  //
  // set default cut on fired macropad 
  // 
  SetCutNmaxFiredMacropad(5); 
  AliInfo(Form("Setting cut on fired macropad to default values: NfiredMacropad = %i", fgCutNmaxFiredMacropad));
  return;
}
