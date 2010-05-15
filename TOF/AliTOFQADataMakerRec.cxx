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


ClassImp(AliTOFQADataMakerRec)
           
//____________________________________________________________________________ 
  AliTOFQADataMakerRec::AliTOFQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kTOF), "TOF Quality Assurance Data Maker"),
  fCalibData(0x0),
  fEnableNoiseFiltering(kFALSE),
    fEnableDqmShifterOpt(kFALSE),
    fProcessedRawEventN(0)
{
  //
  // ctor
  //
   
}

//____________________________________________________________________________ 
AliTOFQADataMakerRec::AliTOFQADataMakerRec(const AliTOFQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fCalibData(qadm.fCalibData),
  fEnableNoiseFiltering(qadm.fEnableNoiseFiltering),
  fEnableDqmShifterOpt(qadm.fEnableDqmShifterOpt),
  fProcessedRawEventN(qadm.fProcessedRawEventN)
{
  //
  //copy ctor 
  //
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
  
 
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

  TH1I * h0 =  new TH1I("hTOFRaws",      "TOF raw hit multiplicity; TOF raw hits number;Counts ",200, 0, 200) ;

  TH1I * h1 =  new TH1I("hTOFRawsIA",    "TOF raw hit multiplicity - I/A side; TOF raw hits number ;Counts ",200, 0, 200) ;
  TH1I * h2 =  new TH1I("hTOFRawsOA",    "TOF raw hit multiplicity - O/A side; TOF raw hits number ;Counts ",200, 0, 200) ;
  TH1I * h3 =  new TH1I("hTOFRawsIC",    "TOF raw hit multiplicity - I/C side; TOF raw hits number ;Counts ",200, 0, 200) ;
  TH1I * h4 =  new TH1I("hTOFRawsOC",    "TOF raw hit multiplicity - O/C side; TOF raw hits number ;Counts ",200, 0, 200) ;

  TH1F * h5  = new TH1F("hTOFRawsTime", "TOF Raws - Hit time (ns);Measured Hit time [ns];Counts", 25000,0. ,610.) ; 

  TH1F * h6  = new TH1F("hTOFRawsTimeIA", "TOF Raws - Hit time (ns) - I/A side;Measured Hit time [ns];Counts", 25000,0. ,610.) ; 
  TH1F * h7  = new TH1F("hTOFRawsTimeOA", "TOF Raws - Hit time (ns) - O/A side;Measured Hit time [ns];Counts", 25000,0. ,610.) ; 
  TH1F * h8  = new TH1F("hTOFRawsTimeIC", "TOF Raws - Hit time (ns) - I/C side;Measured Hit time [ns];Counts", 25000,0. ,610.) ; 
  TH1F * h9  = new TH1F("hTOFRawsTimeOC", "TOF Raws - Hit time (ns) - O/C side;Measured Hit time [ns];Counts", 25000,0. ,610.) ; 

  TH1F * h10  = new TH1F("hTOFRawsToT", "TOF Raws - Hit ToT (ns);Measured Hit ToT (ns);Counts", 1000, 0., 48.8) ; 

  TH1F * h11  = new TH1F("hTOFRawsToTIA", "TOF Raws - Hit ToT (ns) - I/A side;Measured Hit ToT (ns);Counts", 1000, 0., 48.8) ; 
  TH1F * h12  = new TH1F("hTOFRawsToTOA", "TOF Raws - Hit ToT (ns) - O/A side;Measured Hit ToT (ns);Counts", 1000, 0., 48.8) ; 
  TH1F * h13  = new TH1F("hTOFRawsToTIC", "TOF Raws - Hit ToT (ns) - I/C side;Measured Hit ToT (ns);Counts", 1000, 0., 48.8) ; 
  TH1F * h14  = new TH1F("hTOFRawsToTOC", "TOF Raws - Hit ToT (ns) - O/C side;Measured Hit ToT (ns);Counts", 1000, 0., 48.8) ; 

  TH1I * h15 = new TH1I("hTOFRawsLTMHits", "LTM hits ; Crate; Counts",  72, 0., 72.);
  TH1I * h16  = new TH1I("hTOFRawsTRMHits035", "TRM hits  - crates 0 to 35 ;TRM index = SMid(crate*10)+TRM(0-9);Hits",  361, 0., 361.) ;
  TH1I * h17  = new TH1I("hTOFRawsTRMHits3671","TRM hits  - crates 36 to 71 ;TRM index = SMid(crate*10)+TRM(0-9);Hits", 361, 360., 721.) ;
  
  TH1I * h18 = new TH1I("hTOFRawChannelHits","TOF channel hits count; Channel ID; Hits", 158000, 0., 158000);
  
  TH1F * h19  = new TH1F("hTOFOrphansTime", "TOF Raws - Orphans time (ns);Measured Hit time [ns];Counts", 25000, 0. ,610.) ; 
  TH2F * h20 = new TH2F("hTOFRawTimeVsTRM035", "TOF raws - Hit time vs TRM - crates 0 to 35; TRM index = DDL*10+TRM(0-9);TOF raw time [ns]", 361, 0., 361., 250, 0., 610.0) ;
  TH2F * h21 = new TH2F("hTOFRawTimeVsTRM3671", "TOF raws - Hit time vs TRM - crates 36 to 72; TRM index = DDL**10+TRM(0-9);TOF raw time [ns]", 361, 360., 721., 250, 0., 610.0) ;
  TH2F * h22 = new TH2F("hTOFRawToTVsTRM035",  "TOF raws - Hit ToT vs TRM - crates 0 to 35; TRM index = DDL*10+TRM(0-9);TOF raw ToT [ns] ", 361, 0., 361., 100, 0., 48.8) ;
  TH2F * h23 = new TH2F("hTOFRawToTVsTRM3671",  "TOF raws - Hit ToT vs TRM - crates 36 to 72; TRM index = DDL*10+TRM(0-9);TOF raw ToT [ns] ", 361, 360., 721., 100, 0., 48.8) ;
  TH2F * h24 = new TH2F("hTOFTimeVsStrip","TOF Raws - Hit time vs. strip (theta); Strip index;Raws TOF time (ns) ", 91,0.,91, 250, 0., 610.) ; 
  
// TH2F * h25 = new TH2F("hTOFRawsClusMap","Raws vs TOF eta-phi;eta (2*strip+padz);phi (48*sector+padx)",183, -0.5, 182.5,865,-0.5,864.5) ; 

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
  h17->Sumw2() ;
  h18->Sumw2() ;
  h19->Sumw2() ;
  h20->Sumw2() ;
  h21->Sumw2() ;
  h22->Sumw2() ;
  h23->Sumw2() ;
  h24->Sumw2() ;

  Add2RawsList(h0,   0, !expert,  image, !saveCorr) ;
  Add2RawsList(h1,   1,  expert, !image, !saveCorr) ;
  Add2RawsList(h2,   2,  expert, !image, !saveCorr) ;
  Add2RawsList(h3,   3,  expert, !image, !saveCorr) ;
  Add2RawsList(h4,   4,  expert, !image, !saveCorr) ;
  Add2RawsList(h5,   5, !expert,  image, !saveCorr) ;
  Add2RawsList(h6,   6,  expert, !image, !saveCorr) ;
  Add2RawsList(h7,   7,  expert, !image, !saveCorr) ;
  Add2RawsList(h8,   8,  expert, !image, !saveCorr) ;
  Add2RawsList(h9,   9,  expert, !image, !saveCorr) ;
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
  Add2RawsList(h23, 23,  expert, !image, !saveCorr) ;
  Add2RawsList(h24, 24,  expert, !image, !saveCorr) ;

}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::InitRecPoints()
{
  //
  // create RecPoints histograms in RecPoints subdir
  //

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1I * h0 = new TH1I("hTOFRecPoints",    "TOF RecPoints multiplicity ; TOF RecPoints number;Counts",200, 0, 200) ;

  TH1F * h1 = new TH1F("hTOFRecPointsTimeIA", "RecPoints Time Spectrum in TOF (ns)- I/A side; Calibrated TOF time [ns];Counts", 25000, 0., 610.) ; 
  TH1F * h2 = new TH1F("hTOFRecPointsTimeOA", "RecPoints Time Spectrum in TOF (ns)- O/A side; Calibrated TOF time [ns];Counts", 25000, 0., 610.) ;
  TH1F * h3 = new TH1F("hTOFRecPointsTimeIC", "RecPoints Time Spectrum in TOF (ns)- I/C side; Calibrated TOF time [ns];Counts", 25000, 0., 610.) ;
  TH1F * h4 = new TH1F("hTOFRecPointsTimeOC", "RecPoints Time Spectrum in TOF (ns)- O/C side; Calibrated TOF time [ns];Counts", 25000, 0., 610.) ;
  
  TH1F * h5  = new TH1F("hTOFRecPointsRawTimeIA", "RecPoints raw Time Spectrum in TOF (ns)-I/A side; Measured TOF time [ns];Counts", 25000, 0., 610.) ; 
  TH1F * h6  = new TH1F("hTOFRecPointsRawTimeOA", "RecPoints raw Time Spectrum in TOF (ns)-O/A side; Measured TOF time [ns];Counts", 25000, 0., 610.) ; 
  TH1F * h7  = new TH1F("hTOFRecPointsRawTimeIC", "RecPoints raw Time Spectrum in TOF (ns)-I/C side; Measured TOF time [ns];Counts", 25000, 0., 610.) ; 
  TH1F * h8  = new TH1F("hTOFRecPointsRawTimeOC", "RecPoints raw Time Spectrum in TOF (ns)-O/C side; Measured TOF time [ns];Counts", 25000, 0., 610.) ; 
 
  TH1F * h9   = new TH1F("hTOFRecPointsToTIA", "RecPoints ToT Spectrum in TOF (ns)-I/A side; Measured TOT [ns];Counts", 1000, 0., 48.8 ) ; 
  TH1F * h10  = new TH1F("hTOFRecPointsToTOA", "RecPoints ToT Spectrum in TOF (ns)-O/A side; Measured TOT [ns];Counts", 1000, 0., 48.8 ) ; 
  TH1F * h11  = new TH1F("hTOFRecPointsToTIC", "RecPoints ToT Spectrum in TOF (ns)-I/C side; Measured TOT [ns];Counts", 1000, 0., 48.8 ) ; 
  TH1F * h12  = new TH1F("hTOFRecPointsToTOC", "RecPoints ToT Spectrum in TOF (ns)-O/C side; Measured TOT [ns];Counts", 1000, 0., 48.8 ) ; 
  TH1F * h13 = new TH1F("hTOFRawChannelHits","TOF channel hits count; Channel ID; Hits", 158000, 0., 158000);

  TH2F * h14 = new TH2F("hTOFRecPointsClusMap","RecPoints vs TOF phi-eta; eta (2*strip+padz);phi (48*sector+padx)",183, -0.5, 182.5,865,-0.5,864.5) ; 
  TH2F * h15 = new TH2F("hTOFRecTimeVsStrip","RecPoints TOF time vs. strip (theta); Strip index; RecPoints TOF time (ns) ",92,-1.,91, 250, 0., 610.) ;

  TH2F * h16 = new TH2F("hTOFRecPointsTimeVsTRM035","TOF RecPoints time vs TRM - crates 0 to 35; TRM index = DDL*10+TRM(0-9);TOF time [ns]", 361, 0., 361., 250, 0., 610.0) ;
  TH2F * h17 = new TH2F("hTOFRecPointsTimeVsTRM3671","TOF RecPoints time vs TRM - crates 36 to 72; TRM index = DDL**10+TRM(0-9);TOF time [ns]", 361, 360., 721., 250, 0., 610.0) ;
  TH2F * h18 = new TH2F("hTOFRecPointsToTVsTRM035","TOF RecPoints ToT vs TRM - crates 0 to 35; TRM index = DDL*10+TRM(0-9);TOF ToT [ns] ", 361, 0., 361., 100, 0., 48.8) ;
  TH2F * h19 = new TH2F("hTOFRecPointsToTVsTRM3671","TOF RecPoints ToT vs TRM - crates 36 to 72; TRM index = DDL*10+TRM(0-9);TOF ToT [ns] ", 361, 360., 721., 100, 0., 48.8) ;
 
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
    h17->Sumw2() ;
    h18->Sumw2() ;
    h19->Sumw2() ;

  Add2RecPointsList(h0, 0,   !expert,  image) ;
  Add2RecPointsList(h1, 1,   !expert,  image) ;
  Add2RecPointsList(h2, 2,   !expert,  image) ;
  Add2RecPointsList(h3, 3,   !expert,  image) ;
  Add2RecPointsList(h4, 4,   !expert,  image) ;
  Add2RecPointsList(h5, 5,    expert, !image) ;
  Add2RecPointsList(h6, 6,    expert, !image) ;
  Add2RecPointsList(h7, 7,    expert, !image) ;
  Add2RecPointsList(h8, 8,    expert, !image) ;
  Add2RecPointsList(h9, 9,   !expert,  image) ;
  Add2RecPointsList(h10, 10, !expert,  image) ;
  Add2RecPointsList(h11, 11, !expert,  image) ;
  Add2RecPointsList(h12, 12, !expert,  image) ;
  Add2RecPointsList(h13, 13,  expert, !image) ;
  Add2RecPointsList(h14, 14,  expert, !image) ;
  Add2RecPointsList(h15, 15,  expert, !image) ;
  Add2RecPointsList(h16, 16,  expert, !image) ;
  Add2RecPointsList(h17, 17,  expert, !image) ;
  Add2RecPointsList(h18, 18,  expert, !image) ;
  Add2RecPointsList(h19, 19,  expert, !image) ;
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
  TH1F * h1  = new TH1F("hTOFESDsTime", "Matched  ESDs tracks: TOF Time spectrum; Calibrated TOF time [ns];Counts", 25000, 0., 610. ) ; 
  TH1F * h2  = new TH1F("hTOFESDsRawTime", "Matched ESDs tracks: TOF raw Time spectrum;Measured TOF time [ns];Counts", 25000, 0., 610.) ; 
  TH1F * h3  = new TH1F("hTOFESDsToT", "Matched ESDs tracks: TOF ToT spectrum; ESDs ToT [ns];Counts",1000, 0., 48.8) ; 
  TH1F * h4  = new TH1F("hTOFESDsPIDoverMatched", "Fraction of TOF identified tracks over matched TOF tracks per event; identified by TOF / matched tracks [%];Counts", 101, -1., 100.) ;  
  TH1F * h5  = new TH1F("hTOFESDsPID", "Fraction of TOF identified tracks over ESD identified tracks per event;  identified by TOF / ESD identified [%];Counts", 101, -1., 100.) ;  
  TH1F * h6  = new TH1F("hTOFESDsTPCmatched", "TPC-TOF matched tracks momentum distribution (GeV/c); p (GeV/c);Counts", 50, 0.20, 5.00) ;  
  TH1F * h7  = new TH1F("hTOFESDsMatchingProb", "TPC-TOF track-matching probability per event;TPC-TOF track-matching probability (%)  ;Counts",101, -1.0, 100) ;  
  TH1F * h8  = new TH1F("hTOFESDsHitOverTracked", "Fraction of TOF matching tracks over propagated tracks per event; TOF matched / propagated to TOF [%];Counts",101, -1.0, 100.0) ;  
  TH1F * h9  = new TH1F("hTOFESDsDiffTime", "ESDs t_{TOF}-t_{exp} spectrum in TOF (ps); t_{TOF}-t_{exp} [ps];Counts", 200, -2440., 2440.) ; 
  TH1F * h10  = new TH1F("hTOFHitsLength", "Matched ESDs tracks: Length Spectrum; Track length [cm];Counts", 800, 0., 800) ; 
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

  Add2ESDsList(h0, 0, !expert,  image) ;
  Add2ESDsList(h1, 1, !expert,  image) ;
  Add2ESDsList(h2, 2,  expert,  image) ;
  Add2ESDsList(h3, 3, !expert,  image) ;
  Add2ESDsList(h4, 4,  expert,  image) ;
  Add2ESDsList(h5, 5,  expert,  image) ;
  Add2ESDsList(h6, 6,  expert,  image) ; 
  Add2ESDsList(h7, 7,  expert,  image) ; 
  Add2ESDsList(h8, 8,  expert,  image) ; 
  Add2ESDsList(h9, 9, !expert,  image) ;
  Add2ESDsList(h10, 10, !expert,  image) ;
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
	Int_t out[5]; //   out=(indexZ,indexPhi)   
	Int_t chIndex=-1;
	
	TClonesArray * clonesRawData;
	AliTOFRawStream tofInput(rawReader);
	
	//uncomment if needed to apply DeltaBC correction
	//tofInput.ApplyBCCorrections(kTRUE);
	
	for (Int_t iDDL = 0; iDDL < AliTOFGeometry::NDDL()*AliTOFGeometry::NSectors(); iDDL++){
	    rawReader->Reset();
	    
	    tofInput.LoadRawDataBuffers(iDDL);
	    clonesRawData = (TClonesArray*)tofInput.GetRawData();
	    for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {
		AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);
		
		if (tofRawDatum->GetTOF()){
		    
		    equipmentID[0]=iDDL;
		    equipmentID[1]=tofRawDatum->GetTRM(); 
		    equipmentID[2]=tofRawDatum->GetTRMchain();
		    equipmentID[3]=tofRawDatum->GetTDC();
		    equipmentID[4]=tofRawDatum->GetTDCchannel();
		    
		    if (CheckEquipID(equipmentID)){
			tofInput.EquipmentId2VolumeId(iDDL, 
						      tofRawDatum->GetTRM(), 
						      tofRawDatum->GetTRMchain(),
						      tofRawDatum->GetTDC(), 
						      tofRawDatum->GetTDCchannel(), 
						      volumeID);
			if (FilterSpare(equipmentID)) continue;
			if (FilterLTMData(equipmentID)){ //counts LTM hits
			    if (tofRawDatum->GetTOT()) GetRawsData(15)->Fill(equipmentID[0]);
			} else {
			    if (CheckVolumeID(volumeID)){  
				
				GetMapIndeces(volumeID,out);
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
					if (volumeID2[0]>4 && volumeID2[0]<14){       //I side
					    if ((iDDL%4==0)|| (iDDL%4==1)){ //A side
						ntof[1]++;
						GetRawsData(6)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;
						GetRawsData(11)->Fill( tofRawDatum->GetTOT()*tot2ns) ;
					    } else {
						if ((iDDL%4==2)|| (iDDL%4==3)){//C side
						    ntof[3]++;
						    GetRawsData(8)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;
						    GetRawsData(13)->Fill( tofRawDatum->GetTOT()*tot2ns) ;
						}
					    }
					} else {                                    
					    if (volumeID2[0]<5 || volumeID2[0]>13){   //O side
						if ((iDDL%4==0)|| (iDDL%4==1)){ //A side
						    ntof[2]++;
						    GetRawsData(7)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;
						    GetRawsData(12)->Fill( tofRawDatum->GetTOT()*tot2ns) ;
						} else {
						    if ((iDDL%4==2)|| (iDDL%4==3)){//C side
							ntof[4]++;
							GetRawsData(9)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;
							GetRawsData(14)->Fill( tofRawDatum->GetTOT()*tot2ns) ;
						    }
						}
					    }	
					}
					GetRawsData(18)->Fill(chIndex);
					//compute TRM offset
					Int_t trm= iDDL*10+(equipmentID[1]-3);
					if (iDDL>=0 && iDDL<36) {
					    GetRawsData(16)->Fill(trm) ;
					    GetRawsData(20)->Fill(trm,tofRawDatum->GetTOF()*tdc2ns);
					    GetRawsData(22)->Fill(trm,tofRawDatum->GetTOT()*tot2ns);
					}
					if (iDDL>=36 && iDDL<72) {
					    GetRawsData(17)->Fill(trm) ;//in ns 
					    GetRawsData(21)->Fill(trm,tofRawDatum->GetTOF()*tdc2ns);
					    GetRawsData(23)->Fill(trm,tofRawDatum->GetTOT()*tot2ns);
					}				
					GetRawsData(24)->Fill(GetStripIndex(volumeID),tofRawDatum->GetTOF()*tdc2ns) ;
					//GetRawsData(25)->Fill( out[0],out[1]) ;//raw map
				    }//noise filter
				}//end hit selection
				else { //orphans
				    if (!(fCalibData->GetNoiseStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad)
					&& (fCalibData->GetHWStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFHWOk))
					GetRawsData(19)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;//in ns
				}//end orphans
			    }//end volumeID check
			}//end LTM filter
		    }//end equipID check
		}
	    } //while loop
	    clonesRawData->Clear();
	} // DDL Loop
	
	for (Int_t j=0;j<5;j++){
	    GetRawsData(j)->Fill(ntof[j]);
	}
	fProcessedRawEventN++;
	
    } else {
	AliDebug(1,Form("Event of type %d found. Skipping non-physics event for QA.\n", rawReader->GetType())); 
    }
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

  Int_t volumeID[5];
  Int_t volumeID2[5];
  Int_t out[5];

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
   
      volumeID[0] = c->GetDetInd(0);
      volumeID[1] = c->GetDetInd(1);
      volumeID[2] = c->GetDetInd(2);
      volumeID[3] = c->GetDetInd(4); //X and Z indeces inverted in RecPoints
      volumeID[4] = c->GetDetInd(3); //X and Z indeces inverted in RecPoints
      
      for (Int_t i=0;i<5;i++){
	  volumeID2[i]=c->GetDetInd(i); //X and Z indeces inverted in RecPoints
      }
      
      GetMapIndeces(volumeID,out);
      Int_t chIndex=AliTOFGeometry::GetIndex(volumeID2);
      Int_t iDDL=AliTOFRawStream::Geant2DDL(volumeID2);
      Int_t iTRM=AliTOFRawStream::Geant2TRM(volumeID2);

      if (c->GetTDCRAW() && c->GetTDC() && c->GetToT()){
	  if (volumeID2[0]>4 && volumeID2[0]<14){       //I side
	      if ((iDDL%4==0)|| (iDDL%4==1)){ //A side
		  GetRecPointsData(1)->Fill( c->GetTDC()*tdc2ns) ;//in ns
		  GetRecPointsData(5)->Fill( c->GetTDCRAW()*tdc2ns) ;//in ns
		  GetRecPointsData(9)->Fill( c->GetToT()*tot2ns) ;//in ns
		  
	      } else {
		  if ((iDDL%4==2)|| (iDDL%4==3)){//C side
		      GetRecPointsData(3)->Fill( c->GetTDC()*tdc2ns) ;//in ns
		      GetRecPointsData(7)->Fill( c->GetTDCRAW()*tdc2ns) ;//in ns
		      GetRecPointsData(11)->Fill( c->GetToT()*tot2ns) ;//in ns
		  }
	      }
	  } else {
	      if (volumeID2[0]<5 || volumeID2[0]>13){       //O side
		  if ((iDDL%4==0)|| (iDDL%4==1)){ //A side
		      GetRecPointsData(2)->Fill( c->GetTDC()*tdc2ns) ;//in ns
		      GetRecPointsData(6)->Fill( c->GetTDCRAW()*tdc2ns) ;//in ns
		      GetRecPointsData(10)->Fill( c->GetToT()*tot2ns) ;//in ns
		  } else {
		      if ((iDDL%4==2)|| (iDDL%4==3)){//C side
			  GetRecPointsData(4)->Fill( c->GetTDC()*tdc2ns) ;//in ns
			  GetRecPointsData(8)->Fill( c->GetTDCRAW()*tdc2ns) ;//in ns
			  GetRecPointsData(12)->Fill( c->GetToT()*tot2ns) ;//in ns
		      }
		  }
	      }
	  }
	  GetRecPointsData(14)->Fill(out[0],out[1]);
	  GetRecPointsData(15)->Fill(GetStripIndex(volumeID), c->GetTDC()*tdc2ns) ;
	  Int_t trm= iDDL*10+(iTRM-3);
	  if (iDDL>=0 && iDDL<36) {
	      GetRecPointsData(16)->Fill(trm,c->GetTDC()*tdc2ns);
	      GetRecPointsData(18)->Fill(trm,c->GetToT()*tot2ns);
	  }
	  if (iDDL>=36 && iDDL<72) {
	      GetRecPointsData(17)->Fill(trm,c->GetTDC()*tdc2ns);
	      GetRecPointsData(19)->Fill(trm,c->GetToT()*tot2ns);
	  }
	  GetRecPointsData(13)->Fill(chIndex) ;//in ns
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
    Int_t ntof=0;
    Int_t ntofpid=0;
    Int_t nesdpid=0;
    Int_t ntpc=0;
    Int_t ntpctof=0;
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
	    if ((status&AliESDtrack::kTOFout)!=0) {
		ntofout++;
		if (tofTime>0){
		    ntof++;
		    
		    if (TMath::Abs(y)<0.9) {
			GetESDsData(6)->Fill(track->Pt());
			ntpctof++;
		    }
		    GetESDsData(1)->Fill(tofTime*1E-3);
		    GetESDsData(2)->Fill(tofTimeRaw*1E-3); 
		    GetESDsData(3)->Fill(tofToT*1E-3);
		    //check how many tracks where ESD PID is ok 
		    if ((status&AliESDtrack::kESDpid)!=0) nesdpid++; 
		    if (((status&AliESDtrack::kESDpid)&AliESDtrack::kTOFpid)!=0) ntofpid++;
		    
		    Double_t length =track->GetIntegratedLength();
		    Double_t mom2=(track->Pt()*track->Pt())+(track->Pz()*track->Pz());
		    Double_t eTexp = TMath::Sqrt(1+(pionMass*pionMass/mom2))*length/speedOfLight; //in ps
		    GetESDsData(9)->Fill(tofTime-eTexp);
		    GetESDsData(10)->Fill(length);
		} //end check on matched tracks
	    }
	}//end check on TPCrefit
    }
    
    GetESDsData(0)->Fill(ntof) ;
  
    
    if(ntof>0) {
	Float_t ratio = (Int_t)ntofpid/(Int_t)ntof*100; //identified by TOF over matched
	GetESDsData(4)->Fill(ratio) ;
    }
    
    if(nesdpid>0) {
	Float_t ratio = (Float_t)ntofpid/(Float_t)nesdpid *100; //identified by TOF over identified
	GetESDsData(5)->Fill(ratio) ;
    }
    
    if(ntpc>0){
	Float_t ratio = (Float_t)ntof/(Float_t)ntpc*100.; //matching probability
	GetESDsData(7)->Fill(ratio) ;
    }
    
    if(ntofout>0) {
	Float_t ratio = (Float_t)ntof/(Float_t)ntofout*100; //matched over propagated to TOF outer radius
	GetESDsData(8)->Fill(ratio) ;
    }
    EnableDqmShifterOpt(kFALSE);
}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::StartOfDetectorCycle()
{
  //
  //Detector specific actions at start of cycle
 
  fCalibData = GetCalibData();

}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
     for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
	if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
	    continue ; 
	
	AliInfo(Form("Processed %i physics raw events",fProcessedRawEventN));

	if (fEnableDqmShifterOpt){
	    // Help make the raw qa histogram easier to interpret for the DQM shifter
	    // This is still to be optimized...

	    if (!GetRawsData(0) || !GetRawsData(5) || !GetRawsData(10) 
		|| !GetRawsData(15) || !GetRawsData(16) || !GetRawsData(17)) {
		printf("No histogram for DQM found - Possible memory corruption ???. Please check\n") ; 
		continue;
	    }
	    
	    Double_t yTimeMin = GetRawsData(5)->GetMinimum();
	    Double_t yTimeMax = GetRawsData(5)->GetMaximum();
	    Double_t yTotMin = GetRawsData(10)->GetMinimum();
	    Double_t yTotMax = GetRawsData(10)->GetMaximum();
	    
	    TLine* lineExpTimeMin = new TLine(200., yTimeMin, 200., yTimeMax);
	    lineExpTimeMin->SetLineColor(kGreen);
	    lineExpTimeMin->SetLineWidth(2);
	    
	    TLine* lineExpTimeMax = new TLine(250., yTimeMin, 250., yTimeMax);
	    lineExpTimeMax->SetLineColor(kGreen);
	    lineExpTimeMax->SetLineWidth(2);
	    
	    TLine* lineExpTotMin = new TLine( 5., yTotMin, 5., yTotMax);
	    lineExpTotMin->SetLineColor(kGreen);
	    lineExpTotMin->SetLineWidth(2);
	    
	    TLine* lineExpTotMax = new TLine(20., yTotMin, 20., yTotMax);
	    lineExpTotMax->SetLineColor(kGreen);
	    lineExpTotMax->SetLineWidth(2);
	    
	    ((TH1F*)GetRawsData(5))->GetListOfFunctions()->Add(lineExpTimeMin);
	    ((TH1F*)GetRawsData(5))->GetListOfFunctions()->Add(lineExpTimeMax);
	    ((TH1F*)GetRawsData(10))->GetListOfFunctions()->Add(lineExpTotMin);
	    ((TH1F*)GetRawsData(10))->GetListOfFunctions()->Add(lineExpTotMax);
			    
	    //make up for all bistos 
	    for(Int_t j=0;j<20;j++){
	      if (j<5) {
		GetRawsData(j)->SetMarkerColor(kRed);
		GetRawsData(j)->SetMarkerStyle(7);
	      } else {
		GetRawsData(j)->SetLineColor(kBlue+1);
		GetRawsData(j)->SetLineWidth(1);
		GetRawsData(j)->SetMarkerColor(kBlue+1);
	      }
	    }
	    Float_t ySMmax035=GetRawsData(16)->GetMaximum();
	    TLine* lineSMid035[10];
	    Float_t ySMmax3671=GetRawsData(17)->GetMaximum();
	    TLine* lineSMid3671[10];
	    
	    for (Int_t sm=0;sm<10;sm++){
	      lineSMid035[sm] = new TLine( 40*sm, 0, 40*sm, ySMmax035);
	      lineSMid035[sm]->SetLineColor(kMagenta);
	      lineSMid035[sm]->SetLineWidth(2);
	      GetRawsData(16)->GetListOfFunctions()->Add(lineSMid035[sm]);
	      
	      lineSMid3671[sm] = new TLine( 40*sm+360, 0, 40*sm+360, ySMmax3671);
	      lineSMid3671[sm]->SetLineColor(kMagenta);
	      lineSMid3671[sm]->SetLineWidth(2);
	      GetRawsData(17)->GetListOfFunctions()->Add(lineSMid3671[sm]);
	    }
	    
	    for (Int_t j=15;j<19;j++){
	      GetRawsData(j)->SetFillColor(kGray+1);
	      GetRawsData(j)->SetOption("bar");
	    }
	    
	}
    }
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
