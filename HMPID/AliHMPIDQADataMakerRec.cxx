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
#include "AliHMPIDQADataMakerRec.h"
#include "AliHMPIDParam.h"
#include "AliHMPIDRawStream.h"
#include "AliLog.h"
ClassImp(AliHMPIDQADataMakerRec)
           
//____________________________________________________________________________ 
  AliHMPIDQADataMakerRec::AliHMPIDQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kHMPID), "HMPID Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliHMPIDQADataMakerRec::AliHMPIDQADataMakerRec(const AliHMPIDQADataMakerRec& qadm) :
  AliQADataMakerRec() 
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliHMPIDQADataMakerRec& AliHMPIDQADataMakerRec::operator = (const AliHMPIDQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliHMPIDQADataMakerRec();
  new(this) AliHMPIDQADataMakerRec(qadm);
  return *this;
}
 
//____________________________________________________________________________ 

void AliHMPIDQADataMakerRec::InitRecPoints()
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
}
//____________________________________________________________________________

void AliHMPIDQADataMakerRec::InitRaws()
{
//
// Booking QA histo for Raw data
//
  const Int_t nerr = (Int_t)AliHMPIDRawStream::kSumErr+1;
  TH1F *hqPad[14], *hSumErr[14];

  for(Int_t iddl =0; iddl<AliHMPIDRawStream::kNDDL; iddl++) {
    hqPad[iddl] = new TH1F(Form("hqPadDDL%i",iddl), Form("Pad Q Entries at DDL %i",iddl), 500,0,5000);
    Add2RawsList(hqPad[iddl],iddl);
    hSumErr[iddl] = new TH1F(Form("SumErrDDL%i",iddl), Form("Error summary for ddl %i",iddl), 2*nerr,0,2*nerr);
    hSumErr[iddl]->SetYTitle("%");
    
    for(Int_t ilabel=0; ilabel< nerr; ilabel++) {
      hSumErr[iddl]->GetXaxis()->CenterLabels(kTRUE);
      //hSumErr[iddl]->GetXaxis()->SetBinLabel((2*ilabel+1),Form("%i  %s",ilabel+1,hnames[ilabel]));
      hSumErr[iddl]->GetXaxis()->SetBinLabel((2*ilabel+1),Form("%i  %s",ilabel+1,AliHMPIDRawStream::GetErrName(ilabel)));
      }
      
    Add2RawsList(hSumErr[iddl],iddl+14);
 }
  TH1F *hNevRaws = new TH1F("NevRaws","Events per DDL",15,0,15);
  Add2RawsList(hNevRaws,28);
}
//____________________________________________________________________________
void AliHMPIDQADataMakerRec::InitESDs()
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
}
//____________________________________________________________________________
void AliHMPIDQADataMakerRec::MakeRaws(AliRawReader *rawReader)
{
//
// Filling Raws QA histos
//
    AliHMPIDRawStream stream(rawReader);

    while(stream.Next())
     {
       UInt_t ddl=stream.GetDDLNumber(); //returns 0,1,2 ... 13

       for(Int_t iPad=0;iPad<stream.GetNPads();iPad++) {
       GetRawsData(ddl)->Fill(stream.GetChargeArray()[iPad]);}
                                                           
       GetRawsData(28)->Fill(ddl);

       for(Int_t iErr =1; iErr<(Int_t)AliHMPIDRawStream::kSumErr; iErr++){
 
         Int_t NumOfErr = stream.GetErrors(ddl,iErr);

         GetRawsData(ddl+14)->Fill(iErr,NumOfErr);
       }
     }
   stream.Delete();
}
//___________________________________________________________________________
void AliHMPIDQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  //
  //filling QA histos for clusters
  //
  TClonesArray *clusters = new TClonesArray("AliHMPIDCluster");
  for(int i=AliHMPIDParam::kMinCh;i<=AliHMPIDParam::kMaxCh;i++){
    TBranch *branch = clustersTree->GetBranch(Form("HMPID%d",i));
    branch->SetAddress(&clusters);
    branch->GetEntry(0);

    GetRecPointsData(0)->Fill(i,clusters->GetEntries());
    TIter next(clusters);
    AliHMPIDCluster *clu;
    while ( (clu = dynamic_cast<AliHMPIDCluster *>(next())) ) {
      GetRecPointsData(1)->Fill(clu->Chi2());
      GetRecPointsData(2)->Fill(clu->Status());
      GetRecPointsData(3)->Fill(clu->Size());
      GetRecPointsData(4)->Fill(clu->Q()); 
    }
  }

  clusters->Delete();
  delete clusters;
}

//____________________________________________________________________________
void AliHMPIDQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //
  //fills QA histos for ESD
  //
  for(Int_t iTrk = 0 ; iTrk < esd->GetNumberOfTracks() ; iTrk++){
    AliESDtrack *pTrk = esd->GetTrack(iTrk) ;
    GetESDsData(0)->Fill(pTrk->GetP(),pTrk->GetHMPIDsignal());
    GetESDsData(1)->Fill( pTrk->GetP(),TMath::Sqrt(pTrk->GetHMPIDchi2()));
    Float_t xm,ym; Int_t q,np;  
    pTrk->GetHMPIDmip(xm,ym,q,np);                       //mip info
    GetESDsData(2)->Fill(xm,ym);
    Float_t xRad,yRad,th,ph;        
    pTrk->GetHMPIDtrk(xRad,yRad,th,ph);              //track info at the middle of the radiator
    Float_t xPc = xRad+9.25*TMath::Tan(th)*TMath::Cos(ph); // temporar: linear extrapol (B=0!)
    Float_t yPc = yRad+9.25*TMath::Tan(th)*TMath::Sin(ph); // temporar:          "
    GetESDsData(3)->Fill(xm-xPc,ym-yPc); //track info
    Double_t pid[5] ;      pTrk->GetHMPIDpid(pid) ;
    for(Int_t i = 0 ; i < 5 ; i++) GetESDsData(4+i)->Fill(pid[i]) ;
  }
}
//____________________________________________________________________________
void AliHMPIDQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

void AliHMPIDQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX, TObjArray *)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  //AliQAChecker::Instance()->Run(AliQA::kHMPID, task, obj);
  //if(task==AliQA::kRAWS){  

//   for(Int_t iddl=0; iddl<14; iddl++)
//  {
//  if(GetRawsData(28)->GetBinContent(iddl)!=0) GetRawsData(iddl+14)->Scale(100./GetRawsData(28)->GetBinContent(iddl));
//  }
}

