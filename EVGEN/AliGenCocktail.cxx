#include "AliGenerator.h"
#include "AliGenCocktail.h"
#include "TGeant3.h"
#include "AliRun.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <stdlib.h>
ClassImp(AliGenCocktailEntry)
void AliGenCocktailEntry::PrintInfo()
{
printf("\n Generator: %s Generated Events: %d First: %d Last: %d",
       (const char *) fName, fGenerator->NumberParticles(), fFirst, fLast);
}

ClassImp(AliGenCocktail)

AliGenCocktail::AliGenCocktail()
                 :AliGenerator()
{
    fEntries = new TList;
    fNGenerators=0;
}


AliGenCocktail::~AliGenCocktail()
{
    delete fEntries;
}

void AliGenCocktail::
AddGenerator(AliGenerator *Generator, TString Name, Float_t RateExp)
{
//
//  Forward parameters to the new generator
    Generator->SetPtRange(fPtMin,fPtMax);
    Generator->SetMomentumRange(fPMin,fPMax);
    Generator->SetYRange(fYMin,fYMax);
    Generator->
	SetPhiRange(fPhiMin*180/TMath::Pi(),fPhiMax*180/TMath::Pi());
    Generator->
	SetThetaRange(fThetaMin*180/TMath::Pi(),fThetaMax*180/TMath::Pi());
    Generator->
	SetOrigin(fOrigin[0], fOrigin[1], fOrigin[2]);
    Generator->
	SetSigma(fOsigma[0], fOsigma[1], fOsigma[2]);
    Generator->SetVertexSmear(fVertexSmear);
    Generator->SetTrackingFlag(fTrackIt);    
//
//  Add generator to list   
    AliGenCocktailEntry *Entry = 
	new AliGenCocktailEntry(Generator, Name, RateExp);
     fEntries->Add(Entry);
     fNGenerators++;
 }

  void AliGenCocktail::Init()
  {
      TIter next(fEntries);
      AliGenCocktailEntry *Entry;
     //
     // Loop over generators and initialize
     while((Entry = (AliGenCocktailEntry*)next())) {
 	Entry->Generator()->Init();
     }  
 }

 void AliGenCocktail::Generate()
 {
     TIter next(fEntries);
     AliGenCocktailEntry *Entry;
     AliGenCocktailEntry *e1;
     AliGenCocktailEntry *e2;
     TClonesArray *PartArray = gAlice->Particles();
     //
     // Loop over generators and generate events
     Int_t igen=0;
     while((Entry = (AliGenCocktailEntry*)next())) {
	 igen++;
	 if (igen ==1) {
	     Entry->SetFirst(0);
	 } else {
	     Entry->SetFirst((PartArray->GetEntriesFast())+1);
	 }
	 Entry->Generator()->Generate();
	 Entry->SetLast(PartArray->GetEntriesFast());
     }  
     next.Reset();
     while((Entry = (AliGenCocktailEntry*)next())) {
	 Entry->PrintInfo();
     }
     for (Entry=FirstGenerator();
	  Entry;
	  Entry=NextGenerator()
	 ) {
	 Entry->PrintInfo();
     }
     for (FirstGeneratorPair(e1,e2);
	  (e1&&e2);
	  NextGeneratorPair(e1,e2)
	 ){
	 printf("\n -----------------------------");
	 e1->PrintInfo();
	 e2->PrintInfo();
     }
 }

AliGenCocktailEntry *  AliGenCocktail::FirstGenerator()
{
    flnk1 = fEntries->FirstLink();
    if (flnk1) {
	return (AliGenCocktailEntry*) (flnk1->GetObject());
    } else {
	return 0;
    }
}

AliGenCocktailEntry*  AliGenCocktail::NextGenerator()
{
    flnk1 = flnk1->Next();
    if (flnk1) {
	return (AliGenCocktailEntry*) (flnk1->GetObject());
    } else {
	return 0;
    }
}

void AliGenCocktail::
FirstGeneratorPair(AliGenCocktailEntry*& e1, AliGenCocktailEntry*& e2)
{
    flnk2 = flnk1 = fEntries->FirstLink();
    if (flnk1) {
	e2 = e1 = (AliGenCocktailEntry*) (flnk1->GetObject());
    } else {
	e2= e1 = 0;
    }
}

void AliGenCocktail::
NextGeneratorPair(AliGenCocktailEntry*& e1, AliGenCocktailEntry*& e2)
{
    flnk2 = flnk2->Next();
    if (flnk2) {
	e1 = (AliGenCocktailEntry*) (flnk1->GetObject());
	e2 = (AliGenCocktailEntry*) (flnk2->GetObject());	
    } else {
	flnk2 = flnk1 = flnk1->Next();
	if (flnk1) {
	    e1 = (AliGenCocktailEntry*) (flnk1->GetObject());
	    e2 = (AliGenCocktailEntry*) (flnk2->GetObject());
	} else {
	    e1=0;
	    e2=0;
	}
    }
}


void AliGenCocktail::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliGenCocktail.
     TIter next(fEntries);
     AliGenCocktailEntry *Entry;

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliGenerator::Streamer(R__b);
      R__b >> fNGenerators;
      R__b >> fEntries;
// Stream generation related information
      while((Entry = (AliGenCocktailEntry*)next())) {
	  Entry->Streamer(R__b);
      }  
   } else {
       R__b.WriteVersion(AliGenCocktail::IsA());
       AliGenerator::Streamer(R__b);
       R__b << fNGenerators;
       R__b << fEntries;
// Stream generation related information
      while((Entry = (AliGenCocktailEntry*)next())) {
	  Entry->Streamer(R__b);
      }  
   }
}


