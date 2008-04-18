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
           
//____________________________________________________________________________ 
  AliHMPIDQADataMakerRec::AliHMPIDQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kHMPID), "HMPID Quality Assurance Data Maker")
{
  fEvtRaw=0;
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

  TProfile *hCluMult          = new TProfile("CluMult"   ,"Cluster multiplicity per chamber"    , 16, -1 , 7  , 0, 500);
  Add2RecPointsList(hCluMult    , 0);
  
  TH2F *hCluFlg          = new TH2F("CluFlg"      ,"Cluster flag  "                              ,  56  ,-1.5, 12.5, 70, -0.5, 6.5);
  Add2RecPointsList(hCluFlg    , 1);
  
  TH1F *hCluSizeMip[7], *hCluSizePho[7];
  
  TH1F *hCluQSect[42], *hCluQSectZoom[42];
  
  for(Int_t iCh =0; iCh <7; iCh++){
   hCluSizeMip[iCh] = new TH1F(Form("CluSizeMipCh%i",iCh),Form("Cluster size  MIP  (cluster Q > 100 ADC) in Chamber %i",iCh),  50  , 0  , 50  );
  Add2RecPointsList(hCluSizeMip[iCh], iCh+2);

   hCluSizePho[iCh]  = new TH1F(Form("CluSizePho%i",iCh ),Form("Cluster size  Phots(cluster Q < 100 ADC) in Chamber %i",iCh),  50  , 0  , 50  );
  Add2RecPointsList(hCluSizePho[iCh], iCh+7+2);

    for(Int_t iSect =0; iSect < 6; iSect++){
      hCluQSectZoom[iCh*6+iSect] = new TH1F(Form("QClusCh%iSect%iZoom",iCh,iSect) ,Form("Zoom on Cluster charge (ADC) in Chamber %i and sector %i",iCh,iSect),100,0,100);
      Add2RecPointsList(hCluQSectZoom[iCh*6+iSect],2+14+iCh*6+iSect);

      hCluQSect[iCh*6+iSect] = new TH1F(Form("QClusCh%iSect%i",iCh,iSect) ,Form("Cluster charge (ADC) in Chamber %i and sector %i",iCh,iSect),250,0,5000);
      Add2RecPointsList(hCluQSect[iCh*6+iSect],2+14+42+iCh*6+iSect);

    }  
  }
}
//____________________________________________________________________________

void AliHMPIDQADataMakerRec::InitRaws()
{
//
// Booking QA histo for Raw data
//
  const Int_t kNerr = (Int_t)AliHMPIDRawStream::kSumErr+1;
  TH1F *hSumErr[14];
  TH2F *hDilo[14];

  for(Int_t iddl =0; iddl<AliHMPIDRawStream::kNDDL; iddl++) {
    
    hSumErr[iddl] = new TH1F(Form("SumErrDDL%i",iddl), Form("Error summary for ddl %i",iddl), 2*kNerr,0,2*kNerr);
    for(Int_t ilabel=0; ilabel< kNerr; ilabel++) {
      hSumErr[iddl]->GetXaxis()->CenterLabels(kTRUE);
      hSumErr[iddl]->GetXaxis()->SetBinLabel((2*ilabel+1),Form("%i  %s",ilabel+1,AliHMPIDRawStream::GetErrName(ilabel)));
      }
      
    Add2RawsList(hSumErr[iddl],iddl);
    
    hDilo[iddl] = new TH2F(Form("DiloDDL%i",iddl),Form("Dilogic response at DDL",iddl),24,0,24,10,0,10);
    
      Add2RawsList(hDilo[iddl],14+iddl);
 }
}
//____________________________________________________________________________
void AliHMPIDQADataMakerRec::InitESDs()
{
  //
  //Booking ESDs histograms
   TH2F*  hCkovP  = new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]"   , 150,      0,  7  ,100, 0, 1)   ;
   TH2F*  hSigP   = new TH2F("SigP"  ,"#sigma_{#theta_c} [mrad];[GeV]", 150,      0,  7  ,100, 0, 1)   ;
   TH2F*  hDifXY  = new TH2F("DifXY" ,"diff"                          , 200,    -10, 10  ,200,-10,10)  ;
   TH2F*  hMvsP = new TH2F("MvsP","Reconstructed Mass vs P",60,0,6,1000,0,1) ;
   TH1F*  hPid[5];
   hPid[0] = new TH1F("PidE" ,"electron response"              , 101, -0.005,1.005)             ;
   hPid[1] = new TH1F("PidMu","#mu response"                   , 101, -0.005,1.005)             ;
   hPid[2] = new TH1F("PidPi","#pi response"                   , 101, -0.005,1.005)             ;
   hPid[3] = new TH1F("PidK" ,"K response"                     , 101, -0.005,1.005)             ;
   hPid[4] = new TH1F("PidP" ,"p response"                         ,101, -0.005,1.005)             ;

Add2ESDsList(hCkovP,0);
Add2ESDsList(hSigP ,1);
Add2ESDsList(hDifXY,2);
Add2ESDsList(hMvsP,3);
for(Int_t i=0; i< 5; i++) Add2ESDsList(hPid[i],i+4);


}
//____________________________________________________________________________
void AliHMPIDQADataMakerRec::MakeRaws(AliRawReader *rawReader)
{
//
// Filling Raws QA histos
//
    AliHMPIDRawStream stream(rawReader);

    fEvtRaw++;

    while(stream.Next())
     {
       UInt_t ddl=stream.GetDDLNumber(); //returns 0,1,2 ... 13
 
       for(Int_t iErr =1; iErr<(Int_t)AliHMPIDRawStream::kSumErr; iErr++){
 
         Int_t numOfErr = stream.GetErrors(ddl,iErr);
         GetRawsData(ddl)->Fill(iErr,numOfErr);
       }
       UInt_t word; Int_t Nddl, r, d, a;      
       for(Int_t iPad=0;iPad<stream.GetNPads();iPad++) {
       AliHMPIDDigit dig(stream.GetPadArray()[iPad],stream.GetChargeArray()[iPad]);
       dig.Raw(word,Nddl,r,d,a);
       GetRawsData(ddl+14)->Fill(r,d);
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
  AliHMPIDParam *pPar =AliHMPIDParam::Instance();
  TClonesArray *clusters = new TClonesArray("AliHMPIDCluster");
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    TBranch *branch = clustersTree->GetBranch(Form("HMPID%d",iCh));
    branch->SetAddress(&clusters);
    branch->GetEntry(0);
    GetRecPointsData(0)->Fill(iCh,clusters->GetEntries());
    TIter next(clusters);
    AliHMPIDCluster *clu;
    while ( (clu = dynamic_cast<AliHMPIDCluster *>(next())) ) {
      GetRecPointsData(1)->Fill(clu->Status(),iCh);
      Int_t sect =  pPar->InHVSector(clu->Y());
      if(clu->Q()>100) GetRecPointsData(2+iCh)->Fill(clu->Size());
      else {
        GetRecPointsData(2+7+iCh)->Fill(clu->Size());
	GetRecPointsData(2+14+iCh*6+sect)->Fill(clu->Q());
      }    
      GetRecPointsData(2+14+42+iCh*6+sect)->Fill(clu->Q());
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
    Float_t xRad,yRad,th,ph;        
    pTrk->GetHMPIDtrk(xRad,yRad,th,ph);              //track info at the middle of the radiator
    Float_t xPc = xRad+9.25*TMath::Tan(th)*TMath::Cos(ph); // temporar: linear extrapol (B=0!)
    Float_t yPc = yRad+9.25*TMath::Tan(th)*TMath::Sin(ph); // temporar:          "
    GetESDsData(2)->Fill(xm-xPc,ym-yPc); //track info
    if(pTrk->GetHMPIDsignal()>0) {
     Double_t a = 1.292*1.292*TMath::Cos(pTrk->GetHMPIDsignal())*TMath::Cos(pTrk->GetHMPIDsignal())-1.;
     if(a > 0) {
    Double_t mass = pTrk->P()*TMath::Sqrt(1.292*1.292*TMath::Cos(pTrk->GetHMPIDsignal())*TMath::Cos(pTrk->GetHMPIDsignal())-1.);
    GetESDsData(3)->Fill( pTrk->GetP(),mass);
     }
    }
   Double_t pid[5] ;      pTrk->GetHMPIDpid(pid) ;
    for(Int_t i = 0 ; i < 5 ; i++) GetESDsData(4+i)->Fill(pid[i]) ;
  }
}
//____________________________________________________________________________
void AliHMPIDQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

void AliHMPIDQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray *histos)
{
  //Detector specific actions at end of cycle
  // do the QA checking
 //  AliQAChecker::Instance()->Run(AliQA::kHMPID, task, obj);

  if(task==AliQA::kRAWS) {
    for(Int_t iddl=0;iddl<14;iddl++) {
     TH1F *h = (TH1F*)histos->At(14+iddl); //ddl histos scaled by the number of events 
     h->Scale(1./(Float_t)fEvtRaw);
    }
  }
}

