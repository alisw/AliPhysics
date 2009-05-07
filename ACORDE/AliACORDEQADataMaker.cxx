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



//---
//  Produces the data needed to calculate the quality assurance. 
//  ACORDE QA for Hits, Digits, RAW and ESD's
//
//  Authors:
//  
//  Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//  Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP)
//  Arturo Fernandez Tellez <afernan@mail.cern.ch (FCFM-BUAP)
//
//  Created: June 13th 2008
//---

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h>
#include <TH2F.h> 
#include <TDirectory.h>
#include <TObject.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliACORDEdigit.h" 
#include "AliACORDEhit.h"
#include "AliACORDERecPoint.h"
#include "AliACORDEQADataMaker.h"
#include "AliQAChecker.h"
#include "AliACORDERawReader.h"
#include "AliACORDERawStream.h"

ClassImp(AliACORDEQADataMaker)
           
//____________________________________________________________________________ 
AliACORDEQADataMaker::AliACORDEQADataMaker():AliQADataMaker(AliQAv1::GetDetName(AliQAv1::kACORDE), "ACORDE Quality Assurance Data Maker")
{
	// Acorde QA Data Maker
}

//____________________________________________________________________________ 
AliACORDEQADataMaker::AliACORDEQADataMaker(const AliACORDEQADataMaker& qadm):AliQADataMaker() 
{
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}
//__________________________________________________________________
AliACORDEQADataMaker& AliACORDEQADataMaker::operator = (const AliACORDEQADataMaker& qadm )
{
  // Equal operator.
  this->~AliACORDEQADataMaker();
  new(this) AliACORDEQADataMaker(qadm);
  return *this;
}
//____________________________________________________________________________
void AliACORDEQADataMaker::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

}
 
//____________________________________________________________________________ 
void AliACORDEQADataMaker::InitHits()
{

// create Hits histograms in Hits subdir
        TH1F *fAHitsACORDE[8];

        fAHitsACORDE[0] = new TH1F("hACORDEloss" ,"Energy Loss ",1000,0.,1500.);
        fAHitsACORDE[1] = new TH1F("hACORDEPolar" ," Polar Angle ",90,0.,90.);
        fAHitsACORDE[2] = new TH1F("hACORDEAzimuth" ,"Azimuth Angle  ",360,-180.,180.);
        fAHitsACORDE[3] = new TH1F("hACORDEPx" ,"Px Distribution ",60,-30.,30.);
        fAHitsACORDE[4] = new TH1F("hACORDEPy" ,"Py Distribution ",60,-30.,30.);
        fAHitsACORDE[5] = new TH1F("hACORDEPz" ,"Pz Distribution ",60,-30.,30.);
        fAHitsACORDE[6] = new TH1F("hACORDEPt" ,"Pt Distribution ",60,0.,50.);
        fAHitsACORDE[7] = new TH1F("hACORDEpxpz" ,"Pt Distribution ",100,-50.,50.);

        TH2F *hACORDExy = new TH2F("hACORDExy" ,"Dist. xy",2800,-2400.,1400.,200,-4805.,4825.);
        TH2F *hACORDExz = new TH2F("hACORDExz" ,"Dist.xz ",900,-1500.,2850.,1200,-1000.,4000.);
        TH2F *hACORDEyz = new TH2F("hACORDEyz" ,"Dist.yz ",5,817.,819.,1200,-600.,600.);
        TH2F *hACORDEAzimPol =  new TH2F("hACORDEAzimPol" ,"Azimuth vs Polar ",360,-180.,180.,180,0.,180.);

        for(Int_t i=0; i<8; i++)
           Add2HitsList(fAHitsACORDE[i],i);

         Add2HitsList(hACORDExy,8);
         Add2HitsList(hACORDExz,9);
         Add2HitsList(hACORDEyz,10);
         Add2HitsList(hACORDEAzimPol,11);

}
//____________________________________________________________________________ 
void AliACORDEQADataMaker::InitDigits()
{
// create Digits histograms in Digits subdir

   TH1F *    fhDigitsModule;
   TString   modulename;
   modulename = "hDigitsModule";
   fhDigitsModule = new TH1F(modulename.Data(),"hDigitsModuleSingle",60,0,60);
   Add2DigitsList( fhDigitsModule,0);


}

//____________________________________________________________________________ 
void AliACORDEQADataMaker::InitRaws()
{
  // create Raw histograms in Raw subdir
TH1D *fhACORDEBitPattern[4];
fhACORDEBitPattern[0] = new TH1D("ACORDERawDataSM","ACORDE-SingleMuon",60,1,60);//AcordeSingleMuon BitPattern
fhACORDEBitPattern[1] = new TH1D("ACORDERawDataMM","ACORDE-MultiMuon",60,1,60);//AcordeMultiMuon BitPattern
fhACORDEBitPattern[2] = new TH1D("ACORDERawDataSMM","ACORDE-SingleMuonMultiplicity",60,1,60);//AcordeSingleMuon Multiplicity
fhACORDEBitPattern[3] = new TH1D("ACORDERawDataMMM","ACORDE-MultiMuonMultiplicity",60,1,60);//AcordeMultiMuon Multiplicity
for(Int_t i=0;i<4;i++)
{
        Add2RawsList(fhACORDEBitPattern[i],i);
}
}

//____________________________________________________________________________ 

void AliACORDEQADataMaker::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
  // Not needed for ACORDE by now !!!
}
//____________________________________________________________________________
void AliACORDEQADataMaker::InitESDs()
{
  //create ESDs histograms in ESDs subdir

   TH1F *    fhESDsSingle;
   TH1F *    fhESDsMulti;

   TString   name;

   name = "hESDsSingle";
   fhESDsSingle = new TH1F(name.Data(),"hESDsSingle",60,0,60);
   Add2ESDsList( fhESDsSingle,0);

   name = "hESDsMulti";
   fhESDsMulti = new TH1F(name.Data(),"hESDsMulti",60,0,60);
   Add2ESDsList( fhESDsMulti,1);



}
//____________________________________________________________________________

void AliACORDEQADataMaker::MakeHits(TTree *hitTree)
{

 // Here we fill the QA histos for Hits declared above

        printf("Estamos en make Hits");
        TClonesArray * hits = new TClonesArray("AliACORDEhit",1000);
        TBranch * branch = hitTree->GetBranch("ACORDE");
        if (!branch)
        {
                AliWarning("ACORDE branch in Hit Tree not found");
        }else
        {
                if (branch)
                {
                        branch->SetAddress(&hits);
                }else
                {
                        AliError("Branch ACORDE hit not found");
                        exit(111);
                }
                Int_t ntracks = (Int_t)hitTree->GetEntries();
                if (ntracks<=0) return;
                for(Int_t track=0;track<ntracks;track++)
                {
                        branch->GetEntry(track);
                        Int_t nhits = hits->GetEntriesFast();
                        for(Int_t ihit=0;ihit<nhits;ihit++)
                        {
                                AliACORDEhit *AcoHit = (AliACORDEhit*) hits->UncheckedAt(ihit);
                                if (!AcoHit)
                                {
                                        AliError("The unchecked hit doesn't exist");
                                        break;
                                }
                                GetHitsData(0)->Fill(AcoHit->Eloss());
                                GetHitsData(1)->Fill(AcoHit->PolarAngle());
                                GetHitsData(2)->Fill(AcoHit->AzimuthAngle());
                                GetHitsData(3)->Fill(AcoHit->Px());
                                GetHitsData(4)->Fill(AcoHit->Py());
                                GetHitsData(5)->Fill(AcoHit->Pz());
                                GetHitsData(6)->Fill(TMath::Sqrt( (AcoHit->Px())*(AcoHit->Px())+
                                             (AcoHit->Py())*(AcoHit->Py())));
                                if((AcoHit->Py()) != 0.0 ) GetHitsData(7)->Fill(TMath::ATan(AcoHit->Px()/AcoHit->Py()));
                                GetHitsData(8)->Fill( (Float_t)(AcoHit->X()),(Float_t)(AcoHit->Y()) );
                                GetHitsData(9)->Fill( (Float_t)(AcoHit->X()),(Float_t)(AcoHit->Z()) );
                                GetHitsData(10)->Fill( (Float_t)(AcoHit->Y()),(Float_t)(AcoHit->Z()) );
                                GetHitsData(11)->Fill( (Float_t)(AcoHit->AzimuthAngle()),
                                (Float_t)(AcoHit->PolarAngle()));
                        }
                }
        }



}

//____________________________________________________________________________
void AliACORDEQADataMaker::MakeDigits( TTree *digitsTree)
{
  //fills QA histos for Digits


        TClonesArray * digits = new TClonesArray("AliACORDEdigit",1000);
        TBranch * branch = digitsTree->GetBranch("ACORDEdigit");
        if (!branch)
        {
                AliWarning("ACORDE branch in Digits Tree not found");
        }else
        {
                if (branch)
                {
                        branch->SetAddress(&digits);
                }else
                {
                        AliError("Branch ACORDE digit not found");
                        exit(111);
                }
                Int_t ntracks = (Int_t)digitsTree->GetEntries();
                if (ntracks<=0) return;
                printf("Entries in DigitsTree:%d\n",ntracks);
                for(Int_t track=0;track<ntracks;track++)
                {
                        branch->GetEntry(track);
                        Int_t ndigits = digits->GetEntriesFast();
                        for(Int_t idigit=0;idigit<ndigits;idigit++)
                        {
                                AliACORDEdigit *AcoDigit = (AliACORDEdigit*) digits->UncheckedAt(idigit);
                                if (!AcoDigit)
                                {
                                        AliError("The unchecked digit doesn't exist");
                                        break;
                                }
                                GetDigitsData(0)->Fill(AcoDigit->GetModule()-1);
                        }
                }


        }


}


//____________________________________________________________________________
void AliACORDEQADataMaker::MakeRaws( AliRawReader* rawReader)
{

 //fills QA histos for RAW
  rawReader->Reset();
  AliACORDERawStream rawStream(rawReader);
  size_t contSingle=0;
  size_t contMulti=0;
  UInt_t dy[4];

  bool kroSingle[60],kroMulti[60];
  UInt_t tmpDy;

  for(Int_t m=0;m<60;m++) {kroSingle[m]=0;kroMulti[m]=0;}

if(rawStream.Next())
{
        rawReader->NextEvent();
        rawStream.Reset();
        dy[0]=rawStream.GetWord(0);
        dy[1]=rawStream.GetWord(1);
        dy[2]=rawStream.GetWord(2);
        dy[3]=rawStream.GetWord(3);
        tmpDy=dy[0];
        for(Int_t r=0;r<30;++r)
        {
                kroSingle[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[1];
        for(Int_t r=30;r<60;++r)
        {
                kroSingle[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[2];
        for(Int_t r=0;r<30;++r)
        {
                kroMulti[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[3];
        for(Int_t r=30;r<60;++r)
        {
                kroMulti[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        contSingle=0;
        contMulti=0;
        for(Int_t r=0;r<60;++r)
        {
                if(kroSingle[r]==1)
                {
                        GetRawsData(0)->Fill(r+1);
                        contSingle=contSingle+1;
                }
                if(kroMulti[r]==1)
                {
                        GetRawsData(1)->Fill(r+1);
                        contMulti++;
                }

        }GetRawsData(2)->Fill(contSingle);GetRawsData(3)->Fill(contMulti);
}



}
//____________________________________________________________________________
void AliACORDEQADataMaker::MakeRecPoints(TTree * clustersTree)
{
  //fills QA histos for clusters
  // Not needed for ACORDE by now!!!
}

//____________________________________________________________________________
void AliACORDEQADataMaker::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD

        AliESDACORDE * fESDACORDE= esd->GetACORDEData();
        Int_t *fACORDEMultiMuon =fESDACORDE->GetACORDEMultiMuon();
        Int_t *fACORDESingleMuon=fESDACORDE->GetACORDESingleMuon();

       for(int i=0;i<60;i++){
            if(fACORDESingleMuon[i]==1)
               GetESDsData(0) -> Fill(i);
            if(fACORDEMultiMuon[i]==1)
               GetESDsData(1) -> Fill(i);
        }



}

