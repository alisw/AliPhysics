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

//---
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  A. Mastroserio
//---

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <TH1I.h> 
#include <TDirectory.h>
#include <Riostream.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliHMPIDDigit.h"
#include "AliHMPIDHit.h"
#include "AliHMPIDCluster.h"
#include "AliHMPIDQADataMakerRec.h"

ClassImp(AliHMPIDQADataMakerRec)
           
//____________________________________________________________________________ 
  AliHMPIDQADataMakerRec::AliHMPIDQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kHMPID), "HMPID Quality Assurance Data Maker"),
//  fhHitQdc(0x0), 
//  fhSDigits(0x0),
//  fhDigPcEvt(0x0),
//  fhDigChEvt(0x0),
//  fhDigQ(0x0),
//  fhCluEvt(0x0),
//  fhCluChi2(0x0),
//  fhCluQ(0x0),
  fhCluFlg(0x0), 
  fhCluSize(0x0),  
  fhMipCluSize(0x0),
  fhCkovP(0x0),
  fhSigP(0x0),
  fhMipXY(0x0),
  fhDifXY(0x0)
{
  // ctor
//  for(Int_t i=0; i<7; i++) fhHitMap[i]=0x0;
  for(Int_t j=0; j<5; j++) fhPid[j]=0x0;
//   fDetectorDir = fOutput->GetDirectory(GetName()) ;  
//   if (!fDetectorDir) 
//     fDetectorDir = fOutput->mkdir(GetName()) ;  
}

//____________________________________________________________________________ 
AliHMPIDQADataMakerRec::AliHMPIDQADataMakerRec(const AliHMPIDQADataMakerRec& qadm) :
  AliQADataMakerRec(), 
 // fhHitQdc(qadm.fhHitQdc), 
//  fhSDigits(qadm.fhSDigits),
//  fhDigPcEvt(qadm.fhDigPcEvt),
//  fhDigChEvt(qadm.fhDigChEvt),
//  fhDigQ(qadm.fhDigQ),
  fhCluEvt(qadm.fhCluEvt),
  fhCluChi2(qadm.fhCluChi2),
  fhCluQ(qadm.fhCluQ),
  fhCluFlg(qadm.fhCluFlg),
  fhCluSize(qadm.fhCluSize),
  fhMipCluSize(qadm.fhMipCluSize),
  fhCkovP(qadm.fhCkovP),
  fhSigP(qadm.fhSigP),
  fhMipXY(qadm.fhMipXY),
  fhDifXY(qadm.fhDifXY)
{
  //copy ctor 
//  for(Int_t i=0; i<7; i++) fhHitMap[i]=qadm.fhHitMap[i];
  for(Int_t j=0; j<5; j++) fhPid[j]=qadm.fhPid[j];

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
      fhCluEvt=new TH1F("CluPerEvt","# clusters per chamber",16,-1,7);
      fhCluChi2  =new TH1F("CluChi2"  ,"Chi2 "               ,1000,0,100);
      fhCluQ   =new TH1F("CluQ"   ,"Cluster charge"        ,3000,0,3000);
      fhCluFlg   =new TH1F("CluFlg"   ,"Cluster flag"        ,14,-1.5,12.5);
      fhCluSize  =new TH1F("CluSize"  ,"Raw cluster size    ",100,0,100);
      fhMipCluSize =new TH1F("MipCluSize"  ,"Mip cluster size    ",100,0,100);
}
//____________________________________________________________________________
void AliHMPIDQADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir
     fhCkovP  = new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]"   , 150,   0,  7  ,100, 0, 1)   ;
     fhSigP   = new TH2F("SigP"  ,"#sigma_{#theta_c} [mrad];[GeV]", 150,   0,  7  ,100, 0, 1)   ;
     fhMipXY  = new TH2F("MipXY" ,"mip position"                  , 260,   0,130  ,252, 0,126)  ;
     fhDifXY  = new TH2F("DifXY" ,"diff"                          , 200, -10, 10  ,200,-10,10)  ;
     fhPid[0] = new TH1F("PidE" ,"PID: e yellow #mu magenta"  ,100,0,1)                         ;
     fhPid[1] = new TH1F("PidMu","pid of #mu"                 ,100,0,1)                         ;
     fhPid[2] = new TH1F("PidPi","PID: #pi red K green p blue",100,0,1)                         ;
     fhPid[3] = new TH1F("PidK" ,"pid of K"                   ,100,0,1)                         ;
     fhPid[4] = new TH1F("PidP" ,"pid of p"                   ,100,0,1)                         ;
}

//____________________________________________________________________________
void AliHMPIDQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  //fills QA histos for clusters

  TClonesArray *clusters = new TClonesArray("AliHMPIDCluster");
  for(int i=AliHMPIDParam::kMinCh;i<=AliHMPIDParam::kMaxCh;i++){
    TBranch *branch = clustersTree->GetBranch(Form("HMPID%d",i));
    branch->SetAddress(&clusters);
    branch->GetEntry(0);

    fhCluEvt->Fill(i,clusters->GetEntries());
    TIter next(clusters);
    AliHMPIDCluster *clu;
    while ( (clu = dynamic_cast<AliHMPIDCluster *>(next())) ) {;
      fhCluFlg->Fill(clu->Status());  fhCluChi2->Fill(clu->Chi2());  fhCluSize->Fill(clu->Size());
      fhCluQ->Fill(clu->Q()); 
      Int_t qCut=100;
      if(clu->Q()>qCut) {
	fhMipCluSize->SetTitle(Form("Mip cluster size at a Qcut = %i ADC",qCut));
	fhMipCluSize->Fill(clu->Size());
      }
    }
  }

  clusters->Delete();
  delete clusters;
}

//____________________________________________________________________________
void AliHMPIDQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD
  for(Int_t iTrk = 0 ; iTrk < esd->GetNumberOfTracks() ; iTrk++){
    AliESDtrack *pTrk = esd->GetTrack(iTrk) ;
    fhCkovP->Fill(pTrk->GetP(),pTrk->GetHMPIDsignal());
    fhSigP->Fill( pTrk->GetP(),TMath::Sqrt(pTrk->GetHMPIDchi2()));
    Float_t xm,ym; Int_t q,np;  
    pTrk->GetHMPIDmip(xm,ym,q,np);                       //mip info
    fhMipXY->Fill(xm,ym);
    Float_t xRad,yRad,th,ph;        
    pTrk->GetHMPIDtrk(xRad,yRad,th,ph);              //track info at the middle of the radiator
    Float_t xPc = xRad+9.25*TMath::Tan(th)*TMath::Cos(ph); // temporar: linear extrapol (B=0!)
    Float_t yPc = yRad+9.25*TMath::Tan(th)*TMath::Sin(ph); // temporar:          "
    fhDifXY->Fill(xm-xPc,ym-yPc); //track info
    Double_t pid[5] ;      pTrk->GetHMPIDpid(pid) ;
    for(Int_t i = 0 ; i < 5 ; i++) fhPid[i]->Fill(pid[i]) ;
  }
}

