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

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TLine.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <Riostream.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliQAChecker.h"
#include "AliLog.h"
#include "AliHMPIDDigit.h"
#include "AliHMPIDHit.h"
#include "AliHMPIDDigit.h"
#include "AliHMPIDCluster.h"
#include "AliHMPIDQADataMakerRec.h"
#include "AliHMPIDQAChecker.h"
#include "AliHMPIDParam.h"
#include "AliHMPIDRawStream.h"
#include "AliLog.h"

//.
// HMPID AliHMPIDQADataMakerRec base class
// for QA of reconstruction
// here also errors are calculated
//.

ClassImp(AliHMPIDQADataMakerRec)
           
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AliHMPIDQADataMakerRec::AliHMPIDQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kHMPID), "HMPID Quality Assurance Data Maker"), fLineDdlDatSizeLow(0x0), fLineDdlDatSizeUp(0x0), fLineDdlPadOCcLow(0x0), fLineDdlPadOCcUp(0x0), fChannel(0) 
{
  // ctor
  for(Int_t i=0; i<6; i++) fModline[i]=0x0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDQADataMakerRec::AliHMPIDQADataMakerRec(const AliHMPIDQADataMakerRec& qadm) :
  AliQADataMakerRec(),fLineDdlDatSizeLow(qadm.fLineDdlDatSizeLow), fLineDdlDatSizeUp(qadm.fLineDdlDatSizeUp), fLineDdlPadOCcLow(qadm.fLineDdlPadOCcLow), fLineDdlPadOCcUp(qadm.fLineDdlPadOCcUp), fChannel(qadm.fChannel)
{
  //copy ctor 
  for(Int_t i=0; i<6; i++) fModline[i]=qadm.fModline[i];
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDQADataMakerRec& AliHMPIDQADataMakerRec::operator = (const AliHMPIDQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliHMPIDQADataMakerRec();
  new(this) AliHMPIDQADataMakerRec(qadm);
  return *this;
}
 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *hDigChEvt = new TH1F("hDigChEvt","Chamber occupancy per event;Occupancy [%];Counts",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);
  TH1F *hDigPcEvt = new TH1F("hDigPcEvt","PC occupancy",156,-1,77);
  TH2F *hDigMap[7];
  TH1F *hDigQ[42];
  for(Int_t iCh =0; iCh < 7; iCh++){
    hDigMap[iCh] = new TH2F(Form("MapCh%i",iCh),Form("Digit Map in Chamber %i",iCh),159,0,159,143,0,143);
    for(Int_t iPc =0; iPc < 6; iPc++ ){
      hDigQ[iCh*6+iPc] = new TH1F(Form("QCh%iPc%i        ",iCh,iPc),Form("Charge of digits (ADC) in Chamber %i and PC %i;Charge [ADC counts];Counts",iCh,iPc),4100,0,4100);
    }
  }
  
  Add2DigitsList(hDigChEvt,0, !expert, image);
  Add2DigitsList(hDigPcEvt,1,expert, !image);
  for(Int_t iMap=0; iMap < 7; iMap++) Add2DigitsList(hDigMap[iMap],2+iMap,expert, !image);
  for(Int_t iH =0; iH < 42 ; iH++) Add2DigitsList(hDigQ[iH]    ,9+iH,expert,!image);
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
 
  TProfile *hCluMult          = new TProfile("CluMult"   ,"Cluster multiplicity per chamber;Chamber Id;# of clusters"    , 16, -1 , 7  , 0, 500);
  Add2RecPointsList(hCluMult    , 0,expert, !image);
  
  TH2F *hCluFlg          = new TH2F("CluFlg"      ,"Cluster flag;??;??"                              ,  56  ,-1.5, 12.5, 70, -0.5, 6.5);
  Add2RecPointsList(hCluFlg    , 1,expert, !image);
  
  TH1F *hCluSizeMip[7], *hCluSizePho[7];
  
  TH1F *hCluQSect[42], *hCluQSectZoom[42];
  
  for(Int_t iCh =0; iCh <7; iCh++){
    hCluSizeMip[iCh] = new TH1F(Form("CluSizeMipCh%i",iCh),Form("Cluster size  MIP  (cluster Q > 100 ADC) in Chamber %i;Size [MIP];Counts",iCh),  50  , 0  , 50  );
    Add2RecPointsList(hCluSizeMip[iCh], iCh+2,expert,!image);

    hCluSizePho[iCh]  = new TH1F(Form("CluSizePho%i",iCh ),Form("Cluster size  Phots(cluster Q < 100 ADC) in Chamber %i;Size [MIP];Counts",iCh),  50  , 0  , 50  );
    Add2RecPointsList(hCluSizePho[iCh], iCh+7+2,expert,!image);
    
    for(Int_t iSect =0; iSect < 6; iSect++){
      hCluQSectZoom[iCh*6+iSect] = new TH1F(Form("QClusCh%iSect%iZoom",iCh,iSect) ,Form("Zoom on Cluster charge (ADC) in Chamber %i and sector %i;Charge [ADC counts];Counts",iCh,iSect),100,0,100);
      Add2RecPointsList(hCluQSectZoom[iCh*6+iSect],2+14+iCh*6+iSect,expert,!image);
      
      hCluQSect[iCh*6+iSect] = new TH1F(Form("QClusCh%iSect%i",iCh,iSect) ,Form("Cluster charge (ADC) in Chamber %i and sector %i;Charge [ADC counts];Counts",iCh,iSect),250,0,5000);
      Add2RecPointsList(hCluQSect[iCh*6+iSect],2+14+42+iCh*6+iSect, !expert, image);
    }  
  }
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::InitRaws()
{
//
// Booking QA histo for Raw data
//
// All histograms implemented in InitRaws are used in AMORE. Any change here should be propagated to the amoreHMP-QA as well!!! (clm)
//  
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  const Int_t kNerr = (Int_t)AliHMPIDRawStream::kSumErr+1;
  TH1F *hSumErr[14];
  TH2F *hDilo[14];
  TH2I *hPadMap[42]; //AMORE monitoring
  TH1I *hPadQ[42]; //AMORE monitoring

  fLineDdlDatSizeLow  = new TLine(0.5,932,14.5,932);   fLineDdlDatSizeLow->SetLineColor(kGreen); fLineDdlDatSizeLow->SetLineWidth(2);
  fLineDdlDatSizeUp   = new TLine(0.5,1500,14.5,1500); fLineDdlDatSizeUp->SetLineColor(kGreen);  fLineDdlDatSizeUp->SetLineWidth(2);
  
  fLineDdlPadOCcLow  = new TLine(0.5,0.086,14.5,0.086);   fLineDdlPadOCcLow->SetLineColor(kGreen); fLineDdlPadOCcLow->SetLineWidth(2);
  fLineDdlPadOCcUp   = new TLine(0.5,0.86,14.5,0.86); fLineDdlPadOCcUp->SetLineColor(kGreen);  fLineDdlPadOCcUp->SetLineWidth(2);       

  for(Int_t modcnt=0; modcnt < 6; modcnt++){ fModline[modcnt] = new TLine(0,(1+modcnt)*144,160,(1+modcnt)*144);  }
      
  for(Int_t iddl =0; iddl<AliHMPIDRawStream::kNDDL; iddl++) {
    
    hSumErr[iddl] = new TH1F(Form("hSumErrDDL%i",iddl), Form("Error summary for DDL %i;??;??",iddl), 2*kNerr,0,2*kNerr);
    for(Int_t ilabel=0; ilabel< kNerr; ilabel++) {
      hSumErr[iddl]->GetXaxis()->CenterLabels(kTRUE);
      hSumErr[iddl]->GetXaxis()->SetBinLabel((2*ilabel+1),Form("%i  %s",ilabel+1,AliHMPIDRawStream::GetErrName(ilabel)));
    }
    
   Add2RawsList(hSumErr[iddl],iddl,expert,!image, !saveCorr);
    
    hDilo[iddl] = new TH2F(Form("hDiloDDL%i",iddl),Form("Dilogic response at DDL %i;Row # ;Dilogic #",iddl),24,1,25,10,1,11);
    Add2RawsList(hDilo[iddl],14+iddl,expert,!image, !saveCorr);
  }//DDL loop
  for(Int_t iCh = AliHMPIDParam::kMinCh; iCh <=AliHMPIDParam::kMaxCh ;iCh++) {
    for(Int_t iPc = AliHMPIDParam::kMinPc; iPc <= AliHMPIDParam::kMaxPc ;iPc++) {
      hPadMap[iPc+6*iCh] = new TH2I(Form("hPadMap_Ch_%i_Pc%i",iCh,iPc),Form("Pad Map of Ch: %i Pc: %i;Pad X;Pad Y;",iCh,iPc),80,0,80,48,0,48);
     Add2RawsList(hPadMap[iPc+6*iCh],28+iPc+6*iCh,expert,!image, !saveCorr); 
      hPadQ[iPc+6*iCh]   = new TH1I(Form("hPadQ_Ch_%i_Pc%i",iCh,iPc),Form("Pad Charge of Ch: %i Pc: %i;Pad Q;Entries;",iCh,iPc),4100,0,4100);
      Add2RawsList(hPadQ[iPc+6*iCh],70+iPc+6*iCh,expert,!image, !saveCorr); 
    }//PC loop
  }//Ch loop  
  
  TH2I *hGeneralErrorSummary = new TH2I("GeneralErrorSummary"," DDL index vs Error type plot", 2*kNerr, 0, 2*kNerr, 2*AliHMPIDRawStream::kNDDL,0,2*AliHMPIDRawStream::kNDDL);
  for(Int_t igenlabel =0 ; igenlabel< kNerr; igenlabel++) hGeneralErrorSummary->GetXaxis()->SetBinLabel((2*igenlabel+1),Form("%i  %s",igenlabel+1,AliHMPIDRawStream::GetErrName(igenlabel)));
  Add2RawsList(hGeneralErrorSummary,14+14+42+42, expert, !image, !saveCorr);
  
  //___ for DQM shifter and eLogBook ___ start
  //___ Double booking of histograms since TProfile cannot be display in summary image
  //___ hence TProfile plots will not be shown in QA and LogBook!
 
  TH1F* hHmpDdlDataSize = new TH1F("hHmpDdlDataSize","HMP Data Size per DDL;;Data Size (Bytes)",14,0.5,14.5);
  hHmpDdlDataSize->Sumw2();
  hHmpDdlDataSize->SetOption("P");
  hHmpDdlDataSize->SetMinimum(0);
  for(Int_t iddl=0;iddl<14;iddl++)  hHmpDdlDataSize->GetXaxis()->SetBinLabel(iddl+1,Form("DDL_%d",1535+iddl+1));
  hHmpDdlDataSize->SetStats(0);hHmpDdlDataSize->SetMinimum(0);hHmpDdlDataSize->SetMarkerStyle(20);
  hHmpDdlDataSize->GetListOfFunctions()->Add(fLineDdlDatSizeLow);
  hHmpDdlDataSize->GetListOfFunctions()->Add(fLineDdlDatSizeUp);    
  Add2RawsList(hHmpDdlDataSize,14+14+42+42+1,!expert,image,saveCorr);   //shifter, image
  
  TH1F *fHmpPadOcc = new TH1F("fHmpPadOcc","HMP Average pad occupancy per DDL;;Pad occupancy (%)",14,0.5,14.5);
  fHmpPadOcc->Sumw2();fHmpPadOcc->SetMinimum(0);
  fHmpPadOcc->SetMinimum(0);
  for(Int_t iddl=0;iddl<14;iddl++)  fHmpPadOcc->GetXaxis()->SetBinLabel(iddl+1,Form("DDL_%d",1535+iddl+1));
  fHmpPadOcc->SetStats(0);fHmpPadOcc->SetMinimum(0);fHmpPadOcc->SetMarkerStyle(20);
  fHmpPadOcc->GetListOfFunctions()->Add(fLineDdlPadOCcLow);
  fHmpPadOcc->GetListOfFunctions()->Add(fLineDdlPadOCcUp);  
  Add2RawsList(fHmpPadOcc,14+14+42+42+2,!expert,image,!saveCorr);       //shifter, image

  TH2F* fHmpBigMap = new TH2F("hHmpBigMap","HMP Sum Q Maps Ch: 0-6;Ch 0-6: pad X;Ch0, Ch1, Ch2, Ch3, Ch4, Ch5, Ch6 pad Y ;Sum Q / Nevt",160,0,160,1008,0,1008);  
  fHmpBigMap->SetStats(0);  fHmpBigMap->SetOption("COLZ");
  for(Int_t modcnt=0; modcnt < 6; modcnt++) fHmpBigMap->GetListOfFunctions()->Add(fModline[modcnt]);  
  Add2RawsList(fHmpBigMap,14+14+42+42+3,!expert,image,!saveCorr);       //shifter, image
   
  TH2F *fHmpHvSectorQ = new TH2F("fHmpHvSectorQ","HMP HV Sector vs Q; Q (ADC);HV Sector (Ch0-Sc0,Ch0-Sc1,...);Entries*Q/Nevt",410,1,4101,42,0,42);
  fHmpHvSectorQ->SetStats(0); fHmpHvSectorQ->SetOption("colz");
  Add2RawsList(fHmpHvSectorQ,14+14+42+42+4,!expert,image,!saveCorr);    //shifter, image
  
  // TProfiles
  TProfile* hHmpDdlDataSizePrf = new TProfile("hHmpDdlDataSizePrf","HMP Data Size per DDL;;Data Size (Bytes)",14,0.5,14.5);
  hHmpDdlDataSizePrf->Sumw2();
  hHmpDdlDataSizePrf->SetOption("P");
  hHmpDdlDataSizePrf->SetMinimum(0);
  for(Int_t iddl=0;iddl<14;iddl++)  hHmpDdlDataSizePrf->GetXaxis()->SetBinLabel(iddl+1,Form("DDL_%d",1535+iddl+1));
  hHmpDdlDataSizePrf->SetStats(0);hHmpDdlDataSizePrf->SetMinimum(0);hHmpDdlDataSizePrf->SetMarkerStyle(20);
  Add2RawsList(hHmpDdlDataSizePrf,14+14+42+42+5,expert,!image,saveCorr);   //expert, no image
  
  TProfile *fHmpPadOccPrf = new TProfile("fHmpPadOccPrf","HMP Average pad occupancy per DDL;;Pad occupancy (%)",14,0.5,14.5);
  fHmpPadOccPrf->Sumw2();fHmpPadOccPrf->SetMinimum(0);
  fHmpPadOccPrf->SetMinimum(0);
  for(Int_t iddl=0;iddl<14;iddl++)  fHmpPadOccPrf->GetXaxis()->SetBinLabel(iddl+1,Form("DDL_%d",1535+iddl+1));
  fHmpPadOccPrf->SetStats(0);fHmpPadOccPrf->SetMinimum(0);fHmpPadOccPrf->SetMarkerStyle(20);
  Add2RawsList(fHmpPadOccPrf,14+14+42+42+6,expert,!image,saveCorr);       //expert, no image
  //___ for DQM shifter and eLogBook ___ stop
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::InitESDs()
{
  //
  //Booking ESDs histograms
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH2F*  hCkovP  = new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]"   , 150,      0,  7  ,100, 0, 1)   ;
  TH2F*  hSigP   = new TH2F("SigP"  ,"#sigma_{#theta_c} [mrad];[GeV]", 150,      0,  7  ,100, 0, 1)   ;
  TH2F*  hDifXY  = new TH2F("DifXY" ,"diff"                          , 200,    -10, 10  ,200,-10,10)  ;
  TH2F*  hMvsP = new TH2F("MvsP","Reconstructed Mass vs P",60,0,6,1000,0,1) ;
  TH1F*  hPid[5];
  hPid[0] = new TH1F("PidE" ,"electron response"              , 101, -0.005,1.005)             ;
  hPid[1] = new TH1F("PidMu","#mu response"                   , 101, -0.005,1.005)             ;
  hPid[2] = new TH1F("PidPi","#pi response"                   , 101, -0.005,1.005)             ;
  hPid[3] = new TH1F("PidK" ,"K response"                     , 101, -0.005,1.005)             ;
  hPid[4] = new TH1F("PidP" ,"p response"                     , 101, -0.005,1.005)             ;
  
  Add2ESDsList(hCkovP,0, !expert, image);
  Add2ESDsList(hSigP ,1, expert, !image);
  Add2ESDsList(hDifXY,2, !expert, image);
  Add2ESDsList(hMvsP,3, expert, !image);
  for(Int_t i=0; i< 5; i++) Add2ESDsList(hPid[i],i+4, expert, !image);
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::MakeRaws(AliRawReader *rawReader)
{
//
// Filling Raws QA histos
//
    rawReader->Reset() ; 
    Int_t hmpDaqId = AliDAQ::DetectorID("HMPID");                               // shoudl be number 6
    const UInt_t *detPattern = rawReader->GetDetectorPattern(); 
    UInt_t isHmpInRawData = ( ((1 << hmpDaqId) & detPattern[0]) >> hmpDaqId);   // check the 6th bit if HMP is there or not
    if (! isHmpInRawData ) return;                                              // if HMP is not in the event then skip it
    
    AliHMPIDRawStream stream(rawReader);
    Int_t ddlOcc[14]={0};  
    Int_t isHMPin=0;
    UInt_t word; Int_t Nddl, r, d, a;
    Int_t numPadsInDdl;
               
    while(stream.Next())
     {
       UInt_t ddl=stream.GetDDLNumber(); //returns 0,1,2 ... 13   
       if(ddl > 13) continue;
 
     //  FillRawsData(14+14+42+42+1,ddl+1,stream.GetDdlDataSize());
       FillRawsData(14+14+42+42+5,ddl+1,stream.GetDdlDataSize());
       if(stream.GetDdlDataSize() > 0) 
	 {
	   isHMPin++;
	   //___ fill error histo
           for(Int_t iErr =1; iErr<(Int_t)AliHMPIDRawStream::kSumErr; iErr++){
	     Int_t numOfErr = stream.GetErrors(ddl,iErr);
	     FillRawsData(ddl,iErr,numOfErr);
	     FillRawsData(14+14+42+42,iErr,ddl,iErr); //
           }
	   
	   numPadsInDdl= stream.GetNPads();
           ddlOcc[ddl] = numPadsInDdl;
           FillRawsData(14+14+42+42+6,ddl+1,numPadsInDdl/11520.0*100.0);
            
	   //___ loop on pads from raw data from a ddl
	   for(Int_t iPad=0;iPad<numPadsInDdl;iPad++) {
	     AliHMPIDDigit dig(stream.GetPadArray()[iPad],stream.GetChargeArray()[iPad]);dig.Raw(word,Nddl,r,d,a);    
	     //for DQM shifter 
	     FillRawsData(14+14+42+42+3,dig.PadChX(), dig.Ch()*144+dig.PadChY(),dig.Q());
	     FillRawsData(14+14+42+42+4,dig.Q(),(ddl/2*6)+dig.PadChY()/24,dig.Q());
            
	     FillRawsData(ddl+14,r,d);
	     FillRawsData(28+stream.Pc(Nddl,r,d,a)+6*AliHMPIDParam::DDL2C(ddl),stream.PadPcX(Nddl,r,d,a),stream.PadPcY(Nddl,r,d,a));
	     FillRawsData(70+stream.Pc(Nddl,r,d,a)+6*AliHMPIDParam::DDL2C(ddl),stream.GetChargeArray()[iPad]);
            // FillRawsData(14+14+42+42+6,ddl+1,1);
	   }//pad loop
         }     
     }//next
    
    
    if(isHMPin > 0) { // RS: instead of former fEvtRaw
      IncEvCountCycleRaws();
      IncEvCountTotalRaws();
    }
     
    
}//MakeRaws


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::MakeDigits()
{
  //
  //filling QA histos for Digits
  //

  Int_t i = fChannel ; 
  FillDigitsData(0,i,fDigitsArray->GetEntriesFast()/(48.*80.*6.));
  TIter next(fDigitsArray); 
  AliHMPIDDigit * digit; 
  while ( (digit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
    FillDigitsData(1,10.*i+digit->Pc(),1./(48.*80.));
    FillDigitsData(2+i,digit->PadChX(),digit->PadChY());
    FillDigitsData(9+i*6+digit->Pc(),digit->Q());
  }  
}  
  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::MakeDigits(TTree * digTree)
{
  //
  //Opening the Digit Tree
  //

  if(fDigitsArray) 
    fDigitsArray->Clear() ; 
  else
    fDigitsArray=new TClonesArray("AliHMPIDDigit");
  
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    fChannel = iCh ; 
    fDigitsArray->Clear() ; 
    TBranch *branch = digTree->GetBranch(Form("HMPID%d",iCh));
    branch->SetAddress(&fDigitsArray);
    branch->GetEntry(0); 
    MakeDigits();
  }
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  //
  //filling QA histos for clusters
  //
  AliHMPIDParam *pPar =AliHMPIDParam::Instance();
 
  if (fRecPointsArray) 
    fRecPointsArray->Clear() ; 
  else 
    fRecPointsArray = new TClonesArray("AliHMPIDCluster");
  
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    TBranch *branch = clustersTree->GetBranch(Form("HMPID%d",iCh));
    branch->SetAddress(&fRecPointsArray);
    branch->GetEntry(0);
    FillRecPointsData(0,iCh,fRecPointsArray->GetEntries());
    TIter next(fRecPointsArray);
    AliHMPIDCluster *clu;
    while ( (clu = dynamic_cast<AliHMPIDCluster *>(next())) ) {
      FillRecPointsData(1,clu->Status(),iCh);
      Int_t sect =  pPar->InHVSector(clu->Y());
      if(clu->Q()>100) FillRecPointsData(2+iCh,clu->Size());
      else {
        FillRecPointsData(2+7+iCh,clu->Size());
        FillRecPointsData(2+14+iCh*6+sect,clu->Q());
      }    
      FillRecPointsData(2+14+42+iCh*6+sect,clu->Q());
    }
  }
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //
  //fills QA histos for ESD
  //
 
  for(Int_t iTrk = 0 ; iTrk < esd->GetNumberOfTracks() ; iTrk++){
    AliESDtrack *pTrk = esd->GetTrack(iTrk) ;
    Float_t thetaCkov = -999.;
    if(pTrk->GetHMPIDsignal()<0.) thetaCkov = pTrk->GetHMPIDsignal();
    else                          thetaCkov = pTrk->GetHMPIDsignal() - (Int_t)pTrk->GetHMPIDsignal();;
    FillESDsData(0,pTrk->GetP(),thetaCkov);
    FillESDsData(1, pTrk->GetP(),TMath::Sqrt(pTrk->GetHMPIDchi2()));
    Float_t xm,ym; Int_t q,np;  
    pTrk->GetHMPIDmip(xm,ym,q,np);                       //mip info
    Float_t xRad,yRad,th,ph;        
    pTrk->GetHMPIDtrk(xRad,yRad,th,ph);              //track info at the middle of the radiator
    Float_t xPc = xRad+9.25*TMath::Tan(th)*TMath::Cos(ph); // temporar: linear extrapol (B=0!)
    Float_t yPc = yRad+9.25*TMath::Tan(th)*TMath::Sin(ph); // temporar:          "
    FillESDsData(2,xm-xPc,ym-yPc); //track info
    if(pTrk->GetHMPIDsignal()>0) {
     Double_t a = 1.292*1.292*TMath::Cos(thetaCkov)*TMath::Cos(thetaCkov)-1.;
     if(a > 0) {
    Double_t mass = pTrk->P()*TMath::Sqrt(1.292*1.292*TMath::Cos(thetaCkov)*TMath::Cos(thetaCkov)-1.);
    FillESDsData(3, pTrk->GetP(),mass);
     }
    }
   Double_t pid[5] ;      pTrk->GetHMPIDpid(pid) ;
    for(Int_t i = 0 ; i < 5 ; i++) FillESDsData(4+i,pid[i]) ;
  }
  //
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  //
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void AliHMPIDQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray **histos)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  ResetEventTrigClasses(); // reset triggers list to select all histos
  TH1* htmp = 0;
  //
  for (int itc=-1;itc<GetNTrigClasses();itc++) { // RS: loop over eventual clones per trigger class
    //
    if(task==AliQAv1::kRAWS) {
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
	if (! IsValidEventSpecie(specie, histos) ) continue;
	SetEventSpecie(AliRecoParam::ConvertIndex(specie));
	Int_t nEvtRaw = GetEvCountCycleRaws(itc);
	if (nEvtRaw>0) {
	  //ddl histos scaled by the number of events 
	  for(Int_t iddl=0;iddl<14;iddl++) if ( (htmp=GetData(histos, 14+iddl,itc)) ) htmp->Scale(1.0/(Float_t)nEvtRaw);
	  if ( (htmp=GetData(histos, 14+14+42+42+3, itc)) ) htmp->Scale(1.0/(Float_t)nEvtRaw);
	  if ( (htmp=GetData(histos, 14+14+42+42+4, itc)) ) htmp->Scale(1.0/(Float_t)nEvtRaw);
	}      
	Double_t binval=0,binerr=0;	
	TH1F     *h4    = (TH1F*)    GetData(histos, 14+14+42+42+1, itc);
	TProfile *h4prf = (TProfile*)GetData(histos, 14+14+42+42+5, itc);
	if (h4 && h4prf) {
	  for(Int_t iddl=1;iddl<=14;iddl++) {
	    binval=h4prf->GetBinContent(iddl);  binerr=h4prf->GetBinError(iddl);
	    h4->SetBinContent(iddl,binval);     h4->SetBinError(iddl,binerr);
	  }
	}
	//if (nEvtRaw>0) h4->Scale(1./(Float_t)nEvtRaw);
	//
	TH1F     *h5    = (TH1F*)    GetData(histos, 14+14+42+42+2, itc);
	TProfile *h5prf = (TProfile*)GetData(histos, 14+14+42+42+6, itc);
	if (h4 && h4prf) {
	  for(Int_t iddl=1;iddl<=14;iddl++) {
	    binval=h5prf->GetBinContent(iddl);  binerr=h5prf->GetBinError(iddl);
	    h5->SetBinContent(iddl,binval);     h5->SetBinError(iddl,binerr);
	  }
	}
	//if(nEvtRaw>0) h5->Scale(1./(Float_t)nEvtRaw/11520.0*100.0);           
      }//specie loop     
    }//RAWS
    //
  } // RS: loop over eventual clones per trigger class
  //
  AliQAChecker::Instance()->Run(AliQAv1::kHMPID, task, histos);
  //
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
