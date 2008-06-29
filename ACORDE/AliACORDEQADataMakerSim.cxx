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
//  All data must be mergeable objects.

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
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliACORDEdigit.h" 
#include "AliACORDEhit.h"
#include "AliACORDEQADataMakerSim.h"
#include "AliQAChecker.h"
#include "AliACORDERawReader.h"
ClassImp(AliACORDEQADataMakerSim)
           
//____________________________________________________________________________ 
AliACORDEQADataMakerSim::AliACORDEQADataMakerSim():AliQADataMakerSim(AliQA::GetDetName(AliQA::kACORDE), "ACORDE Quality Assurance Data Maker")
{
}
//____________________________________________________________________________ 
AliACORDEQADataMakerSim::AliACORDEQADataMakerSim(const AliACORDEQADataMakerSim& qadm) :
  AliQADataMakerSim() 
{
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}
//__________________________________________________________________
AliACORDEQADataMakerSim& AliACORDEQADataMakerSim::operator = (const AliACORDEQADataMakerSim& qadm )
{
  // Equal operator.
  this->~AliACORDEQADataMakerSim();
  new(this) AliACORDEQADataMakerSim(qadm);
  return *this;
}
//____________________________________________________________________________
void AliACORDEQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
   printf("ACORDE---->Detector specific actions at END of cycle\n................\n");

  AliQAChecker::Instance()->Run(AliQA::kACORDE, task, list) ;
}
//____________________________________________________________________________
void AliACORDEQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  printf("ACORDE---->Detector specific actions at START of cycle\n................\n");
}
//____________________________________________________________________________ 
void AliACORDEQADataMakerSim::InitHits()
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
void AliACORDEQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir

   TH1F *    fhDigitsModule;
   TString   modulename;
   modulename = "hDigitsModule";
   fhDigitsModule = new TH1F(modulename.Data(),"hDigitsModuleSingle",60,0,60);
   Add2DigitsList( fhDigitsModule,0);

}
//____________________________________________________________________________

void AliACORDEQADataMakerSim::MakeHits(TTree *hitTree)
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
void AliACORDEQADataMakerSim::MakeDigits( TTree *digitsTree)
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
