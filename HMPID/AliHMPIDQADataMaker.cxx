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
#include <TH1F.h> 
#include <TH2F.h>
#include <TProfile.h>
#include <Riostream.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliQAChecker.h"
#include "AliLog.h"
#include "AliHMPIDDigit.h"
#include "AliHMPIDHit.h"
#include "AliHMPIDCluster.h"
#include "AliHMPIDQADataMaker.h"
#include "AliHMPIDParam.h"
#include "AliHMPIDRawStream.h"
#include "AliLog.h"
ClassImp(AliHMPIDQADataMaker)
           
//____________________________________________________________________________ 
  AliHMPIDQADataMaker::AliHMPIDQADataMaker() : 
  AliQADataMaker(AliQAv1::GetDetName(AliQAv1::kHMPID), "HMPID Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliHMPIDQADataMaker::AliHMPIDQADataMaker(const AliHMPIDQADataMaker& qadm) :
  AliQADataMaker() 
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliHMPIDQADataMaker& AliHMPIDQADataMaker::operator = (const AliHMPIDQADataMaker& qadm )
{
  // Equal operator.
  this->~AliHMPIDQADataMaker();
  new(this) AliHMPIDQADataMaker(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliHMPIDQADataMaker::InitHits()
{
  // create Hits histograms in Hits subdir
  TH1F *hHitQdc=new TH1F("HitQdc","HMPID Hit Qdc all chamber;QDC",500,0,4000);
  Add2HitsList(hHitQdc,0);
  TH2F *hHitMap[7];
  for(Int_t iCh=0;iCh<7;iCh++) {
    hHitMap[iCh]=new TH2F(Form("HMPID HitMap%i",iCh),Form("Ch%i;x_{Hit};y_{Hit}",iCh),162,-1,161,146,-1,145);   
    Add2HitsList(hHitMap[iCh],iCh+1);
  }
  //
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line  
}

//____________________________________________________________________________ 
void AliHMPIDQADataMaker::InitDigits()
{
  // create Digits histograms in Digits subdir
  TH1F *hDigPcEvt = new TH1F("hDigPcEvt","PC occupancy",156,-1,77);
  TH1F *hDigQ     = new TH1F("Q        ","Charge of digits (ADC)     ",3000,0,3000);
  TH1F *hDigChEvt = new TH1F("hDigChEvt","Chamber occupancy per event",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);
  
  TProfile *tDigHighQ = new TProfile("tDigHighQ","Highest charge in chamber  ",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);
  TProfile *tDigChEvt = new TProfile("tDigChEvt","Chamber occupancy per event (profile)",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);
  
  Add2DigitsList(hDigPcEvt,0);
  Add2DigitsList(hDigQ    ,1);
  Add2DigitsList(hDigChEvt,2);
  Add2DigitsList(tDigHighQ,3);
  Add2DigitsList(tDigChEvt,4);
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliHMPIDQADataMaker::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  TH1F   *hSDigits     = new TH1F("hHmpidSDigits",    "SDigits Q  distribution in HMPID",  500, 0., 5000.) ; 
  Add2SDigitsList(hSDigits,0);
  //
  ClonePerTrigClass(AliQAv1::kSDIGITS); // this should be the last line
}

//____________________________________________________________________________ 

void AliHMPIDQADataMaker::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir

  TH1F *hCluEvt=new TH1F("CluPerEvt","Cluster multiplicity"   ,100,0,100);
  TH1F *hCluChi2  =new TH1F("CluChi2"  ,"Chi2 "               ,1000,0,100);
  TH1F *hCluFlg   =new TH1F("CluFlg"   ,"Cluster flag"        ,14,-1.5,12.5); hCluFlg->SetFillColor(5);
  TH1F *hCluSize  =new TH1F("CluSize"  ,"Cluster size        ",100,0,100);
  TH1F *hCluQ     =new TH1F("CluQ"     ,"Cluster charge (ADC)",1000,0,5000);

  Add2RecPointsList(hCluEvt , 0);
  Add2RecPointsList(hCluChi2, 1);
  Add2RecPointsList(hCluFlg , 2);
  Add2RecPointsList(hCluSize, 3);
  Add2RecPointsList(hCluQ   , 4);
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
}
//____________________________________________________________________________

void AliHMPIDQADataMaker::InitRaws()
{
  //
  // Booking QA histo for Raw data
  //
  TH1F *hqPad[14];
  for(Int_t iddl =0; iddl<14; iddl++) {
  hqPad[iddl] = new TH1F(Form("hqPadDDL%i",iddl), Form("Pad Q Entries at DDL %i",iddl), 500,0,5000);
  Add2RawsList(hqPad[iddl],iddl);
  }

  const Int_t nerr = (Int_t)AliHMPIDRawStream::kSumErr+1;
  const char *hnames[nerr]={"RawDataSize","RawMarkerSize","WrongRow","WrongDilogic","WrongPad","EoEFlag",
                             "EoESize","EoEDILOGIC","EoERow","BadSegWord","WrongSeg","RowMarkerSize","NoErrors","Invalid"};

  TH1F *hSumErr = new TH1F("SumErr","Summary of the returned errors",2*nerr,0,nerr);

  for(Int_t ilabel=0; ilabel< nerr; ilabel++) {
  hSumErr->GetXaxis()->CenterLabels(kTRUE);
  hSumErr->GetXaxis()->SetBinLabel((2*ilabel+1),Form("%i  %s",ilabel+1,hnames[ilabel]));

  }
  Add2RawsList(hSumErr,14);
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}

//____________________________________________________________________________
void AliHMPIDQADataMaker::InitESDs()
{
  //
  //Booking ESDs histograms
  TH2F*  hCkovP  = new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]"   , 150,      0,  7  ,100, 0, 1)   ;
  TH2F*  hSigP   = new TH2F("SigP"  ,"#sigma_{#theta_c} [mrad];[GeV]", 150,      0,  7  ,100, 0, 1)   ;
  TH2F*  hMipXY  = new TH2F("MipXY" ,"mip position"                  , 260,      0,130  ,252, 0,126)  ;
  TH2F*  hDifXY  = new TH2F("DifXY" ,"diff"                          , 200,    -10, 10  ,200,-10,10)  ;
  TH1F*  hPid[5];
  hPid[0] = new TH1F("PidE" ,"electron response"              , 101, -0.005,1.005)             ;
  hPid[1] = new TH1F("PidMu","#mu response"                   , 101, -0.005,1.005)             ;
  hPid[2] = new TH1F("PidPi","#pi response"                   , 101, -0.005,1.005)             ;
  hPid[3] = new TH1F("PidK" ,"K response"                     , 101, -0.005,1.005)             ;
  hPid[4] = new TH1F("PidP" ,"p response"                         ,101, -0.005,1.005)             ;
  
  Add2ESDsList(hCkovP,0);
  Add2ESDsList(hSigP ,1);
  Add2ESDsList(hMipXY,2);
  Add2ESDsList(hDifXY,3);
  for(Int_t i=0; i< 5; i++) Add2ESDsList(hPid[i],i+4);
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line
}
//____________________________________________________________________________

void AliHMPIDQADataMaker::MakeHits(TClonesArray * data)
{
 //
 //filling QA histos for Hits
 //
  TClonesArray * hits = dynamic_cast<TClonesArray *>(data) ; 
  if (!hits){
    AliError("Wrong type of hits container") ; 
  } else {
    TIter next(hits); 
    AliHMPIDHit * hit ; 
    while ( (hit = dynamic_cast<AliHMPIDHit *>(next())) ) {
      if(hit->Pid()<500000) FillHitsData(0,hit->Q()) ;
      if(hit->Pid()<500000) FillHitsData(hit->Ch()+1,hit->LorsX(),hit->LorsY());
    }
  } 
  //
}
//___________________________________________________________________________
void AliHMPIDQADataMaker::MakeHits(TTree * data)
{
//
//Opening of the Hit TTree 
//
 TClonesArray *pHits=new TClonesArray("AliHMPIDHit");  data->SetBranchAddress("HMPID",&pHits);
  for(Int_t iEnt=0;iEnt<data->GetEntriesFast();iEnt++){//entries loop
    data->GetEntry(iEnt);
    MakeHits(pHits);
  }//entries loop
  //
  IncEvCountCycleHits();
  IncEvCountTotalHits();
  //
}

//____________________________________________________________________________
void AliHMPIDQADataMaker::MakeDigits(TClonesArray * data)
{
 //
 //filling QA histos for Digits
 //
  TObjArray *chamber = dynamic_cast<TObjArray*>(data);
  if ( !chamber) {
    AliError("Wrong type of digits container") ; 
  } else {
    for(Int_t i =0; i< chamber->GetEntries(); i++)
      {
	TClonesArray * digits = dynamic_cast<TClonesArray*>(chamber->At(i)); 
	FillDigitsData(2,i,digits->GetEntriesFast()/(48.*80.*6.));
        FillDigitsData(4,i,digits->GetEntriesFast()/(48.*80.*6.));
        Double_t highQ=0;
	TIter next(digits); 
	AliHMPIDDigit * digit; 
	while ( (digit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
	  FillDigitsData(0,10.*i+digit->Pc(),1./(48.*80.));
	  FillDigitsData(1,digit->Q());
          if(digit->Q()>highQ) highQ = digit->Q();
	}  
	FillDigitsData(3,i,highQ);	
      }
  }
}
//___________________________________________________________________________
void AliHMPIDQADataMaker::MakeDigits(TTree * data)
{
//
//Opening the Digit Tree
//
 TObjArray *pObjDig=new TObjArray(AliHMPIDParam::kMaxCh+1);
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    TClonesArray *pCA=new TClonesArray("AliHMPIDDigit");
    pObjDig->AddAt(pCA,iCh);
  }

  pObjDig->SetOwner(kTRUE);

  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    data->SetBranchAddress(Form("HMPID%i",iCh),&(*pObjDig)[iCh]);
  }
  data->GetEntry(0);

  MakeDigits((TClonesArray *)pObjDig);
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}

//____________________________________________________________________________

void AliHMPIDQADataMaker::MakeRaws(AliRawReader *rawReader)
{
//
// Filling Raws QA histos
//
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) {
    AliHMPIDRawStream stream(rawReader);
    while(stream.Next())
      {
      UInt_t ddl=stream.GetDDLNumber(); //returns 0,1,2 ... 13
      if((UInt_t)(2*iCh)==ddl || (UInt_t)(2*iCh+1)==ddl) {
       for(Int_t row = 1; row <= AliHMPIDRawStream::kNRows; row++){
        for(Int_t dil = 1; dil <= AliHMPIDRawStream::kNDILOGICAdd; dil++){
          for(Int_t pad = 0; pad < AliHMPIDRawStream::kNPadAdd; pad++){
            if(stream.GetCharge(ddl,row,dil,pad) < 1) continue;
              FillRawsData(ddl,stream.GetCharge(ddl,row,dil,pad));
//              Printf("charge %i",stream.GetCharge(ddl,row,dil,pad));
            }//pad
          }//dil
        }//row
      }//while
    }
    for(Int_t iErr =1; iErr<(Int_t)AliHMPIDRawStream::kSumErr; iErr++){
      Int_t errflag = stream.GetErrors(iErr);

       if(errflag < 0) FillRawsData(14,(Int_t)AliHMPIDRawStream::kSumErr+0.5);
       else if(errflag == 0) FillRawsData(14,(Int_t)AliHMPIDRawStream::kSumErr-0.5);
       else FillRawsData(14,iErr-0.5);
     }
    stream.Delete();
  }//chamber loop
  //
  IncEvCountCycleRaws();
  IncEvCountTotalRaws();
  //
}

//___________________________________________________________________________

void AliHMPIDQADataMaker::MakeSDigits(TClonesArray * data)
{
 //
 //filling QA histos for SDigits
 //
  TClonesArray * sdigits = dynamic_cast<TClonesArray *>(data) ; 
  if (!sdigits) {
    AliError("Wrong type of sdigits container") ; 
  } else {
    TIter next(sdigits) ; 
    AliHMPIDDigit * sdigit ; 
    while ( (sdigit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
	    FillSDigitsData(0,sdigit->Q());
    } 
  }
}
//___________________________________________________________________________
void AliHMPIDQADataMaker::MakeSDigits(TTree * data)
{
  //
  // Opening the SDigit Tree
  //
  TClonesArray * sdigits = new TClonesArray("AliHMPIDDigit", 1000) ;

  TBranch * branch = data->GetBranch("HMPID") ;
  if ( ! branch ) {
    AliError("HMPID SDigit Tree not found") ;
    return;
  }
  branch->SetAddress(&sdigits) ;
  branch->GetEntry(0) ;
  MakeSDigits(sdigits) ;
  //
  IncEvCountCycleSDigits();
  IncEvCountTotalSDigits();
  //
}
//____________________________________________________________________________
void AliHMPIDQADataMaker::MakeRecPoints(TTree * clustersTree)
{
  //
  //filling QA histos for clusters
  //
  TClonesArray *clusters = new TClonesArray("AliHMPIDCluster");
  for(int i=AliHMPIDParam::kMinCh;i<=AliHMPIDParam::kMaxCh;i++){
    TBranch *branch = clustersTree->GetBranch(Form("HMPID%d",i));
    branch->SetAddress(&clusters);
    branch->GetEntry(0);

    FillRecPointsData(0,i,clusters->GetEntries());
    TIter next(clusters);
    AliHMPIDCluster *clu;
    while ( (clu = dynamic_cast<AliHMPIDCluster *>(next())) ) {
      FillRecPointsData(1,clu->Chi2());
      FillRecPointsData(2,clu->Status());
      FillRecPointsData(3,clu->Size());
      FillRecPointsData(4,clu->Q()); 
    }
  }

  clusters->Delete();
  delete clusters;
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //
}

//____________________________________________________________________________
void AliHMPIDQADataMaker::MakeESDs(AliESDEvent * esd)
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
    FillESDsData(2,xm,ym);
    Float_t xRad,yRad,th,ph;        
    pTrk->GetHMPIDtrk(xRad,yRad,th,ph);              //track info at the middle of the radiator
    Float_t xPc = xRad+9.25*TMath::Tan(th)*TMath::Cos(ph); // temporar: linear extrapol (B=0!)
    Float_t yPc = yRad+9.25*TMath::Tan(th)*TMath::Sin(ph); // temporar:          "
    FillESDsData(3,xm-xPc,ym-yPc); //track info
    Double_t pid[5] ;      pTrk->GetHMPIDpid(pid) ;
    for(Int_t i = 0 ; i < 5 ; i++) FillESDsData(4+i,pid[i]) ;
  }
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  //
}
//____________________________________________________________________________
void AliHMPIDQADataMaker::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

void AliHMPIDQADataMaker::EndOfDetectorCycle(AliQAv1::TASKINDEX task, TObjArray * obj)
{
  //Detector specific actions at end of cycle
  // do the QA checking
//  AliQAChecker::Instance()->Run(AliQAv1::kHMPID, task, obj) ;  
}

