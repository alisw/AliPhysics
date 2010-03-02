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
// last modified by F. Bellini (fbellini@cern.ch) on 02/03/2010      //
///////////////////////////////////////////////////////////////////////

 
#include <TClonesArray.h>
#include <TH1F.h> 
#include <TH2F.h> 

#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"

#include "AliTOFcluster.h"
#include "AliTOFQADataMakerRec.h"
#include "AliTOFRawStream.h"
#include "AliTOFrawData.h"
#include "AliTOFGeometry.h"
#include "AliTOFdigit.h"
#include "AliTOFChannelOnlineStatusArray.h"


ClassImp(AliTOFQADataMakerRec)
           
//____________________________________________________________________________ 
  AliTOFQADataMakerRec::AliTOFQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kTOF), "TOF Quality Assurance Data Maker"),
  fCalibData(0x0),
  fEnableNoiseFiltering(kFALSE)
{
  //
  // ctor
  //
}

//____________________________________________________________________________ 
AliTOFQADataMakerRec::AliTOFQADataMakerRec(const AliTOFQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fCalibData(qadm.fCalibData),
  fEnableNoiseFiltering(qadm.fEnableNoiseFiltering)
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

  TH1F * h0 = new TH1F("hTOFRaws",      "Number of TOF Raw Hits; TOF raw hits number;Counts ",1001, -1., 1000.) ;
  TH1F * h1 = new TH1F("hTOFRawsIA",    "Number of TOF Raw Hits - I/A side; TOF raw hits number [10 power];Counts ",1001, -1., 1000.) ;
  TH1F * h2 = new TH1F("hTOFRawsOA",    "Number of TOF Raw Hits - O/A side; TOF raw hits number [10 power];Counts ",1001, -1., 1000.) ;
  TH1F * h3 = new TH1F("hTOFRawsIC",    "Number of TOF Raw Hits - I/C side; TOF raw hits number [10 power];Counts ",1001, -1., 1000.) ;
  TH1F * h4 = new TH1F("hTOFRawsOC",    "Number of TOF Raw Hits - O/C side; TOF raw hits number [10 power] ;Counts ",1001, -1., 1000.) ;
  TH1F * h5  = new TH1F("hTOFRawsTimeIA", "Raws Time Spectrum in TOF (ns) - I/A side;Measured raw TOF time [ns];Counts", 100000,0. ,2440.) ; 
  TH1F * h6  = new TH1F("hTOFRawsTimeOA", "Raws Time Spectrum in TOF (ns) - O/A side;Measured raw TOF time [ns];Counts", 100000,0. ,2440.) ; 
  TH1F * h7  = new TH1F("hTOFRawsTimeIC", "Raws Time Spectrum in TOF (ns) - I/C side;Measured raw TOF time [ns];Counts", 100000,0. ,2440.) ; 
  TH1F * h8  = new TH1F("hTOFRawsTimeOC", "Raws Time Spectrum in TOF (ns) - O/C side;Measured raw TOF time [ns];Counts", 100000,0. ,2440.) ; 
  TH1F * h9  = new TH1F("hTOFRawsToTIA",  "Raws ToT Spectrum in TOF (ns) - I/A side;Measured raw ToT [ns];Counts", 10000, 0., 244.) ; 
  TH1F * h10  = new TH1F("hTOFRawsToTOA", "Raws ToT Spectrum in TOF (ns) - O/A side;Measured raw ToT [ns];Counts", 10000, 0., 244.) ; 
  TH1F * h11  = new TH1F("hTOFRawsToTIC", "Raws ToT Spectrum in TOF (ns) - I/C side;Measured raw ToT [ns];Counts", 10000, 0., 244.) ; 
  TH1F * h12  = new TH1F("hTOFRawsToTOC", "Raws ToT Spectrum in TOF (ns) - O/C side;Measured raw ToT [ns];Counts", 10000, 0., 244.) ; 
  TH1F * h13  = new TH1F("hTOFRawsTRMHitRate035", "TRM hit rate - crates 0 to 35 ;TRM index = DDL*10+TRM(0-9);Counts",  361, 0., 361.) ; 
  TH1F * h14  = new TH1F("hTOFRawsTRMHitRate3671","TRM hit rate - crates 36 to 71 ;TRM index = DDL*10+TRM(0-9);Counts", 361, 361., 722.) ; 
  TH1F * h15  = new TH1F("hTOFRawsVolumeErrorCounter","Invalid hit/volume association error counter per DDL; DDL; error counter ",73, -1., 72.) ; 
  TH2F * h16  = new TH2F("hTOFRawsClusMap","Raws vs TOF eta-phi;eta (2*strip+padz);phi (48*sector+padx)",183, -0.5, 182.5,865,-0.5,864.5) ; 
  TH1F * h17  = new TH1F("hTOFOrphansTime", "Raws Time Spectrum in TOF (ns) for hits with no ToT measurement;Measured raw TOF time [ns];Counts", 100000,0. ,2440.) ; 
  TH1F * h18  = new TH1F("hTOFRawsDecodingErrorCounter",  "Invalid hit/equipment association error counter per DDL; DDL; error counter ",73, -1., 72.) ; 
  TH1F * h19 = new TH1F("hTOFRawsWords",    "Number of TOF Raw words per event; TOF raw words number;Counts ",1001, -1., 1000.) ;
  TH2F * h20 = new TH2F("hTOFRawTimeVsDDL",    "TOF raw time spectrum vs DDL; TOF raw time [ns];DDL ",1000,0. ,2440., 72, 0., 72.) ;
  TH2F * h21 = new TH2F("hTOFRawToTVsDDL",     "TOF raw ToT spectrum vs DDL; TOF raw ToT [ns];DDL ",   100, 0., 244., 72, 0., 72.) ;
  
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
      
  Add2RawsList(h0,   0, !expert,  image, !saveCorr) ;
  Add2RawsList(h1,   1,  expert, !image, !saveCorr) ;
  Add2RawsList(h2,   2,  expert, !image, !saveCorr) ;
  Add2RawsList(h3,   3,  expert, !image, !saveCorr) ;
  Add2RawsList(h4,   4,  expert, !image, !saveCorr) ;
  Add2RawsList(h5,   5, !expert,  image, !saveCorr) ;
  Add2RawsList(h6,   6, !expert,  image, !saveCorr) ;
  Add2RawsList(h7,   7, !expert,  image, !saveCorr) ;
  Add2RawsList(h8,   8, !expert,  image, !saveCorr) ;
  Add2RawsList(h9,   9, !expert,  image, !saveCorr) ;
  Add2RawsList(h10, 10, !expert,  image, !saveCorr) ;
  Add2RawsList(h11, 11, !expert,  image, !saveCorr) ;
  Add2RawsList(h12, 12, !expert,  image, !saveCorr) ;
  Add2RawsList(h13, 13,  expert, !image, !saveCorr) ;
  Add2RawsList(h14, 14,  expert, !image, !saveCorr) ;
  Add2RawsList(h15, 15,  expert, !image, !saveCorr) ;
  Add2RawsList(h16, 16, !expert,  image, !saveCorr) ;
  Add2RawsList(h17, 17,  expert, !image, !saveCorr) ;
  Add2RawsList(h18, 18,  expert, !image, !saveCorr) ;
  Add2RawsList(h19, 19,  expert, !image, !saveCorr) ;
  Add2RawsList(h20, 20,  expert, !image, !saveCorr) ;
  Add2RawsList(h21, 21,  expert, !image, !saveCorr) ;

}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::InitRecPoints()
{
  //
  // create RecPoints histograms in RecPoints subdir
  //

  Bool_t expert = kFALSE;

  TH1F * h0 = new TH1F("hTOFRecPoints",    "Number of TOF RecPoints ",301, -1.02, 5.) ;   h0->Sumw2() ;
  Add2RecPointsList(h0, 0, expert) ;

  TH1F * h1  = new TH1F("hTOFRecPointsTime", "RecPoints Time Spectrum in TOF (ns)", 2000, 0., 200) ; 
  h1->Sumw2() ;
  Add2RecPointsList(h1, 1, expert) ;

  TH1F * h2  = new TH1F("hTOFRecPointsRawTime", "RecPoints raw Time Spectrum in TOF (ns)", 2000, 0., 200) ; 
  h2->Sumw2() ;
  Add2RecPointsList(h2, 2, expert) ;

  TH1F * h3  = new TH1F("hTOFRecPointsToT", "RecPoints ToT Spectrum in TOF (ns)", 500, 0., 50) ; 
  h3->Sumw2() ;
  Add2RecPointsList(h3, 3, expert) ;

  TH2F * h4  = new TH2F("hTOFRecPointsClusMap","RecPoints vs TOF eta-phi",183, -0.5, 182.5,865,-0.5,864.5) ; 
  h4->Sumw2() ;
  Add2RecPointsList(h4, 4, expert) ;

}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::InitESDs()
{
  //
  //create ESDs histograms in ESDs subdir
  //

  Bool_t expert = kFALSE;

  TH1F * h0 = new TH1F("hTOFESDs",    "Number of matched TOF tracks over ESDs",       250, -1., 4.) ;  
  h0->Sumw2() ; 
  Add2ESDsList(h0, 0, expert) ;

  TH1F * h1  = new TH1F("hTOFESDsTime", "Time Spectrum in TOF (ns)", 2000, 0., 200) ; 
  h1->Sumw2() ;
  Add2ESDsList(h1, 1, expert) ;

  TH1F * h2  = new TH1F("hTOFESDsRawTime", "raw Time Spectrum in TOF (ns)", 2000, 0., 200) ; 
  h2->Sumw2() ;
  Add2ESDsList(h2, 2, expert) ;

  TH1F * h3  = new TH1F("hTOFESDsToT", "ToT Spectrum in TOF (ns)", 500, 0., 50) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, 3, expert) ;

  TH1F * h4 = new TH1F("hTOFESDsPID",    "Fraction of matched TOF tracks with good PID glag", 100, 0., 1.) ;  
  h4->Sumw2() ; 
  Add2ESDsList(h4, 4, expert) ;
}

//____________________________________________________________________________
void AliTOFQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //
  // makes data from Raws
  //
  
  Double_t tdc2ns=AliTOFGeometry::TdcBinWidth()*1E-3;
  Double_t tot2ns=AliTOFGeometry::ToTBinWidth()*1E-3;

  Int_t nentries=0;
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
	nentries++;
	AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

	if (tofRawDatum->GetTOF()){
	
	    equipmentID[0]=iDDL;
	    equipmentID[1]=tofRawDatum->GetTRM(); 
	    equipmentID[2]=tofRawDatum->GetTRMchain();
	    equipmentID[3]=tofRawDatum->GetTDC();
	    equipmentID[4]=tofRawDatum->GetTDCchannel();
	    
	    if (!CheckEquipID(equipmentID)) GetRawsData(18)->Fill(equipmentID[0]);
	    else {
		tofInput.EquipmentId2VolumeId(iDDL, 
					      tofRawDatum->GetTRM(), 
					      tofRawDatum->GetTRMchain(),
					      tofRawDatum->GetTDC(), 
					      tofRawDatum->GetTDCchannel(), 
					      volumeID);
		
		if (!CheckVolumeID(volumeID)) GetRawsData(15)->Fill(equipmentID[0]);  
		else {
		    GetMapIndeces(volumeID,out);
		    
		    volumeID2[0]=volumeID[0];
		    volumeID2[1]=volumeID[1];
		    volumeID2[2]=volumeID[2];
		    volumeID2[3]=volumeID[4];
		    volumeID2[4]=volumeID[3];
		    chIndex=AliTOFGeometry::GetIndex(volumeID2);
		    
		    if (tofRawDatum->GetTOT()){	    
			if (!(fCalibData->GetNoiseStatus(chIndex)==AliTOFChannelOnlineStatusArray::kTOFNoiseBad)) {//noise filter
			    ntof[0]++; //counter for tof hits
			    if (volumeID2[0]>4 && volumeID2[0]<14){       //I side
				if ((iDDL%4==0)|| (iDDL%4==1)){ //A side
				    ntof[1]++;
				    GetRawsData(5)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;//in ns
				    GetRawsData(9)->Fill( tofRawDatum->GetTOT()*tot2ns) ;//in ns
				} else {
				    if ((iDDL%4==2)|| (iDDL%4==3)){//C side
					ntof[3]++;
					GetRawsData(7)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;//in ns
					GetRawsData(11)->Fill( tofRawDatum->GetTOT()*tot2ns) ;//in ns
				    }
				}
			    } else {
				if (volumeID2[0]<5 || volumeID2[0]>13){       //O side
				    if ((iDDL%4==0)|| (iDDL%4==1)){ //A side
					ntof[2]++;
					GetRawsData(6)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;//in ns
					GetRawsData(10)->Fill( tofRawDatum->GetTOT()*tot2ns) ;//in ns
				    } else {
					if ((iDDL%4==2)|| (iDDL%4==3)){//C side
					    ntof[4]++;
					    GetRawsData(8)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;//in ns
					    GetRawsData(12)->Fill( tofRawDatum->GetTOT()*tot2ns) ;//in ns
					}
				    }
				}
			    }
			    //compute TRM offset
			    Int_t trm= iDDL*10+(volumeID[1]-3);
			    if (iDDL>=0 && iDDL<36) GetRawsData(13)->Fill(trm) ;//in ns 
			    if (iDDL>=36 && iDDL<72) GetRawsData(14)->Fill(trm) ;//in ns 
			    GetRawsData(16)->Fill( out[0],out[1]) ;//raw map
			    
			    GetRawsData(20)->Fill(tofRawDatum->GetTOF()*tdc2ns,iDDL);
			    GetRawsData(21)->Fill(tofRawDatum->GetTOT()*tot2ns,iDDL);
			    
			}//noise filter
		    }//end hit selection
		    else { //orphans
			if (!(fCalibData->GetNoiseStatus(chIndex) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad))
			    GetRawsData(17)->Fill( tofRawDatum->GetTOF()*tdc2ns) ;//in ns
		    }//end orphans
		}//end volumeID check
	    }//end equipID check
	}
    } // while loop
    
    clonesRawData->Clear();
  } // DDL Loop
  
  for (Int_t j=0;j<5;j++){
      if(ntof[j]<=0) GetRawsData(j)->Fill(-1.) ; 
      else GetRawsData(j)->Fill(ntof[j]);
  }
  if (nentries<=0) GetRawsData(19)->Fill(-1.) ;
  else GetRawsData(19)->Fill(nentries) ;
  
}


//____________________________________________________________________________
void AliTOFQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
    //
    // Make data from Clusters
    //
    
    Double_t tdc2ns=AliTOFGeometry::TdcBinWidth()*1E-3;
    Double_t tot2ns=AliTOFGeometry::ToTBinWidth()*1E-3;
    
    Int_t in[5];
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
    
    Int_t nentries=clusters->GetEntriesFast();
    if(nentries<=0){
	GetRecPointsData(0)->Fill(-1.) ; 
    }else{
	GetRecPointsData(0)->Fill(TMath::Log10(nentries)) ; 
    } 
    
    TIter next(clusters) ; 
    AliTOFcluster * c ; 
    while ( (c = dynamic_cast<AliTOFcluster *>(next())) ) {
	GetRecPointsData(1)->Fill(c->GetTDC()*tdc2ns);
	GetRecPointsData(2)->Fill(c->GetTDCRAW()*tdc2ns);
	GetRecPointsData(3)->Fill(c->GetToT()*tot2ns);
	
	in[0] = c->GetDetInd(0);
	in[1] = c->GetDetInd(1);
	in[2] = c->GetDetInd(2);
	in[3] = c->GetDetInd(4); //X and Z indeces inverted in RecPoints
	in[4] = c->GetDetInd(3); //X and Z indeces inverted in RecPoints
	
	GetMapIndeces(in,out);
	GetRecPointsData(4)->Fill(out[0],out[1]);
	
    }
}

//____________________________________________________________________________
void AliTOFQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //
  // make QA data from ESDs
  //  
  Int_t ntrk = esd->GetNumberOfTracks() ; 
  Int_t ntof=0;
  Int_t ntofpid=0;
  while (ntrk--) {
    AliESDtrack *track=esd->GetTrack(ntrk);
    Double_t tofTime=track->GetTOFsignal()*1E-3;//in ns
    Double_t tofTimeRaw=track->GetTOFsignalRaw()*1E-3;//in ns
    Double_t tofToT=track->GetTOFsignalToT(); //in ns
    if(!(tofTime>0))continue;
    ntof++;
    GetESDsData(1)->Fill(tofTime);
    GetESDsData(2)->Fill(tofTimeRaw); 
    GetESDsData(3)->Fill(tofToT);
    //check how many tracks where ESD PID is ok 
    UInt_t status=track->GetStatus();
    if (((status&AliESDtrack::kESDpid)==0) || 
	((status&AliESDtrack::kTOFpid)==0)) continue;
    ntofpid++;
  }
  
  Int_t nentries=ntof;
  if(nentries<=0){
    GetESDsData(0)->Fill(-1.) ;
  }else{
    GetESDsData(0)->Fill(TMath::Log10(nentries)) ;
  }

  if(ntof>0)GetESDsData(4)->Fill(ntofpid/ntof) ;

}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::StartOfDetectorCycle()
{
  //
  //Detector specific actions at start of cycle
  //to be implemented  

  fCalibData = GetCalibData();

}

//____________________________________________________________________________ 
void AliTOFQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking

  AliQAChecker::Instance()->Run(AliQAv1::kTOF, task, list) ;  
}
//____________________________________________________________________________
void AliTOFQADataMakerRec::GetMapIndeces(Int_t* in , Int_t* out)
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
Bool_t  AliTOFQADataMakerRec::CheckVolumeID(Int_t *volumeID)
{
   
    for (Int_t j=0;j<5;j++){
	if (volumeID[j]<0) {
	    AliDebug(1,Form("Invalid detector volume index for volumeID[%i]",j));
	    return kFALSE;
	}
    }
    return kTRUE;
    
}

//---------------------------------------------------------------
Bool_t  AliTOFQADataMakerRec::CheckEquipID(Int_t *equipmentID)
{
   for (Int_t j=0;j<5;j++){
	if (equipmentID[j]<0) {
	  if (j==0)equipmentID[j]=-1;
	  AliDebug(1,Form("Invalid equipment volume index for equipmentID[%i]",j));
	  return kFALSE;
	}
   }
   return kTRUE;
}

