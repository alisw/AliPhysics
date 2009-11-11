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

////////////////////////////////////////////////////////////////////////////////
//
//  Glauber MC implementation
//
//  origin: PHOBOS experiment
//  alification: Mikolaj Krzewicki, Nikhef, mikolaj.krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>
#include <TEllipse.h>
#include <TRandom.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TNtuple.h>
#include <TFile.h>

#include "AliGlauberNucleon.h"
#include "AliGlauberNucleus.h"
#include "AliGlauberMC.h"

ClassImp(AliGlauberMC)

//______________________________________________________________________________
AliGlauberMC::AliGlauberMC(Option_t* NA, Option_t* NB, Double_t xsect) :
   fANucleus(NA),fBNucleus(NB),
   fXSect(0),fNucleonsA(0),fNucleonsB(0),fnt(0),
   fMeanX2(0),fMeanY2(0),fMeanXY(0),fMeanXParts(0),
   fMeanYParts(0),fMeanXSystem(0),fMeanYSystem(0),  
   fMeanX_A(0),fMeanY_A(0),fMeanX_B(0),fMeanY_B(0),fB_MC(0),
   fEvents(0),fTotalEvents(0),fBMin(0),fBMax(0),fMaxNpartFound(0),
   fNpart(0),fNcoll(0),fSx2(0),fSy2(0),fSxy(0)
{
   fBMin = 0;
   fBMax = 20;
   fXSect = xsect;
   
   TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
   TString title(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
   SetName(name);
   SetTitle(title);
}

//______________________________________________________________________________
Bool_t AliGlauberMC::CalcEvent(Double_t bgen)
{
   // prepare event
   fANucleus.ThrowNucleons(-bgen/2.);
   fNucleonsA = fANucleus.GetNucleons();
   fAN = fANucleus.GetN();
   for (Int_t i = 0; i<fAN; i++) {
      AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->At(i));
      nucleonA->SetInNucleusA();
   }
   fBNucleus.ThrowNucleons(bgen/2.);
   fNucleonsB = fBNucleus.GetNucleons();
   fBN = fBNucleus.GetN();
   for (Int_t i = 0; i<fBN; i++) {
      AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->At(i));
      nucleonB->SetInNucleusB();
   }

   // "ball" diameter = distance at which two balls interact
   Double_t d2 = (Double_t)fXSect/(TMath::Pi()*10); // in fm^2

   // for each of the A nucleons in nucleus B
   for (Int_t i = 0; i<fBN; i++) {
      AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->At(i));
      for (Int_t j = 0 ; j < fAN ;j++) {
	 AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->At(j));
         Double_t dx = nucleonB->GetX()-nucleonA->GetX();
         Double_t dy = nucleonB->GetY()-nucleonA->GetY();
         Double_t dij = dx*dx+dy*dy;
         if (dij < d2) {
            nucleonB->Collide();
            nucleonA->Collide();
         }
      }
  }
   
  return CalcResults(bgen);
}

//______________________________________________________________________________
Bool_t AliGlauberMC::CalcResults(Double_t bgen)
{
   // calc results for the given event
   fNpart=0;
   fNcoll=0;
   fMeanX2=0;
   fMeanY2=0;
   fMeanXY=0;
   fMeanXParts=0;
   fMeanYParts=0;
   fMeanXSystem=0;
   fMeanYSystem=0;
   fMeanX_A=0;
   fMeanY_A=0;
   fMeanX_B=0;
   fMeanY_B=0;
  
   for (Int_t i = 0; i<fAN; i++) {
      AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->At(i));
      Double_t xA=nucleonA->GetX();
      Double_t yA=nucleonA->GetY();
      fMeanXSystem  += xA;
      fMeanYSystem  += yA;
      fMeanX_A  += xA;
      fMeanY_A  += yA;

      if(nucleonA->IsWounded()) {
         fNpart++;
         fMeanXParts  += xA;
         fMeanYParts  += yA;
         fMeanX2 += xA * xA;
         fMeanY2 += yA * yA;
         fMeanXY += xA * yA;
      }
   }

   for (Int_t i = 0; i<fBN; i++) {
      AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->At(i));
      Double_t xB=nucleonB->GetX();
      Double_t yB=nucleonB->GetY();
      fMeanXSystem  += xB;
      fMeanYSystem  += yB;
      fMeanX_B  += xB;
      fMeanY_B  += yB;

      if(nucleonB->IsWounded()) {
         fNpart++;
         fMeanXParts  += xB;
         fMeanYParts  += yB;
         fMeanX2 += xB * xB;
         fMeanY2 += yB * yB;
         fMeanXY += xB * yB;
	 fNcoll += nucleonB->GetNColl();
      }
   }

   if (fNpart>0) {
      fMeanXParts /= fNpart;
      fMeanYParts /= fNpart;
      fMeanX2 /= fNpart;
      fMeanY2 /= fNpart;
      fMeanXY /= fNpart;
   } else {
      fMeanXParts = 0;
      fMeanYParts = 0;
      fMeanX2 = 0;
      fMeanY2 = 0;
      fMeanXY = 0;
   }
   
   if(fAN+fBN>0) {
      fMeanXSystem /= (fAN + fBN);
      fMeanYSystem /= (fAN + fBN);
   } else {
      fMeanXSystem = 0;
      fMeanYSystem = 0;
   }
   if(fAN>0) {
      fMeanX_A /= fAN;
      fMeanY_A /= fAN;
   } else {
      fMeanX_A = 0;
      fMeanY_A = 0;
   }

   if(fBN>0) {
      fMeanX_B /= fBN;
      fMeanY_B /= fBN;
   } else {
      fMeanX_B = 0;
      fMeanY_B = 0;
   }

   fSx2=fMeanX2-(fMeanXParts*fMeanXParts);
   fSy2=fMeanY2-(fMeanYParts*fMeanYParts);
   fSxy=fMeanXY-fMeanXParts*fMeanYParts;

   fB_MC = bgen;

   fTotalEvents++;
   if (fNpart>0) fEvents++;

   if (fNpart==0) return kFALSE;
   if (fNpart > fMaxNpartFound) fMaxNpartFound = fNpart;

   return kTRUE;
}

//______________________________________________________________________________
void AliGlauberMC::Draw(Option_t* /*option*/)
{
   fANucleus.Draw(fXSect, 2);
   fBNucleus.Draw(fXSect, 4);

   TEllipse e;
   e.SetFillColor(0);
   e.SetLineColor(1);
   e.SetLineStyle(2);
   e.SetLineWidth(1);
   e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
   e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetTotXSect() const
{
   return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBMax*fBMax/100;
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetTotXSectErr() const
{
   return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * 
      TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

//______________________________________________________________________________
TObjArray *AliGlauberMC::GetNucleons() 
{
   if(!fNucleonsA || !fNucleonsB) return 0;
   fNucleonsA->SetOwner(0);
   fNucleonsB->SetOwner(0);
   TObjArray *allnucleons=new TObjArray(fAN+fBN);
   allnucleons->SetOwner();
   for (Int_t i = 0; i<fAN; i++) {
      allnucleons->Add(fNucleonsA->At(i));
   }
   for (Int_t i = 0; i<fBN; i++) {
      allnucleons->Add(fNucleonsB->At(i));
   }
   return allnucleons;
}

//______________________________________________________________________________
Bool_t AliGlauberMC::NextEvent(Double_t bgen)
{
   if(bgen<0) 
      bgen = TMath::Sqrt((fBMax*fBMax-fBMin*fBMin)*gRandom->Rndm()+fBMin*fBMin);

   return CalcEvent(bgen);
}

//______________________________________________________________________________
void AliGlauberMC::Run(Int_t nevents)
{
   TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
   TString title(Form("%s + %s (x-sect = %d mb)",fANucleus.GetName(),fBNucleus.GetName(),(Int_t) fXSect));
   if (fnt == 0) {
      fnt = new TNtuple(name,title,
                        "Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB");
      fnt->SetDirectory(0);
   }

   for (int i = 0;i<nevents;i++) {

      while(!NextEvent()) {}

      Float_t v[17];
      v[0]  = GetNpart();
      v[1]  = GetNcoll();
      v[2]  = fB_MC;
      v[3]  = fMeanXParts;
      v[4]  = fMeanYParts;
      v[5]  = fMeanX2;
      v[6]  = fMeanY2;
      v[7]  = fMeanXY;
      v[8]  = fSx2;
      v[9]  = fSy2;
      v[10] = fSxy;
      v[11] = fMeanXSystem;
      v[12] = fMeanYSystem;
      v[13] = fMeanX_A;
      v[14] = fMeanY_A;
      v[16] = fMeanX_B;
      v[17] = fMeanY_B;

      fnt->Fill(v);

      if (!(i%50)) std::cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;  
   }
   std::cout << std::endl << "Done!" << std::endl;
}

//---------------------------------------------------------------------------------
void AliGlauberMC::runAndSaveNtuple( Int_t n,
                                     Option_t *sysA,
                                     Option_t *sysB,
                                     Double_t signn,
                                     Double_t mind,
                                     const char *fname)
{
   AliGlauberMC *mcg=new AliGlauberMC(sysA,sysB,signn);
   mcg->SetMinDistance(mind);
   mcg->Run(n);
   TNtuple  *nt=mcg->GetNtuple();
   TFile out(fname,"recreate",fname,9);
   if(nt) nt->Write();
   out.Close();
}

//---------------------------------------------------------------------------------
void AliGlauberMC::runAndSaveNucleons( Int_t n,                    
                                       Option_t *sysA,           
                                       Option_t *sysB,           
                                       Double_t signn,
                                       Double_t mind,
                                       Bool_t verbose,
                                       const char *fname)
{
   AliGlauberMC *mcg=new AliGlauberMC(sysA,sysB,signn);
   mcg->SetMinDistance(mind);
   TFile *out=0;
   if(fname) 
      out=new TFile(fname,"recreate",fname,9);

   for(Int_t ievent=0;ievent<n;ievent++){

      //get an event with at least one collision
      while(!mcg->NextEvent()) {}

      //access, save and (if wanted) print out nucleons
      TObjArray* nucleons=mcg->GetNucleons();
      if(!nucleons) continue;
      if(out)
         nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey);

      if(verbose) {
         cout<<endl<<endl<<"EVENT NO: "<<ievent<<endl;
         cout<<"B = "<<mcg->GetB()<<"  Npart = "<<mcg->GetNpart()<<endl<<endl;
         printf("Nucleus\t X\t Y\t Z\tNcoll\n");
         Int_t nNucls=nucleons->GetEntries();
         for(Int_t iNucl=0;iNucl<nNucls;iNucl++) {
            AliGlauberNucleon *nucl=(AliGlauberNucleon *)nucleons->At(iNucl);
            Char_t nucleus='A';
            if(nucl->IsInNucleusB()) nucleus='B';
            Double_t x=nucl->GetX();
            Double_t y=nucl->GetY();
            Double_t z=nucl->GetZ();
            Int_t ncoll=nucl->GetNColl();
            printf("   %c\t%2.2f\t%2.2f\t%2.2f\t%3d\n",nucleus,x,y,z,ncoll);
         }
      }
   }
   if(out) delete out;
}

