/* $Id$ */

// This class contains a number of histograms for diagnostics of a TPC
// read out chamber from the raw data
//
// TODO:
//  
//
//

#include "AliTPCRawHistograms.h"

#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TNtuple.h>

#include <AliTPCRawStream.h>
#include <AliLog.h>


//____________________________________________________________________
ClassImp(AliTPCRawHistograms)

//____________________________________________________________________
AliTPCRawHistograms::AliTPCRawHistograms() 
  : TNamed(),
    fhDigits(0),
    fhSignal(0)
{
  // default constructor
}

//____________________________________________________________________
AliTPCRawHistograms::AliTPCRawHistograms(Int_t detector, const Char_t* /* comment */, Int_t timeStart, Int_t timeStop)
  : TNamed(),
    fhDigits(0),
    fhSignal(0),
    fDigitTree(0)
{
  // constructor 
  
  // make name and title
  if (detector < 0 || detector >= 72) {
    AliDebug(AliLog::kError, Form("Detector %d does not exist", detector));
    return;
  }
      
  Int_t sector = detector%18;
  TString side;
  TString inout;
  if (detector<18 || ( detector>=36 && detector<54))
    side.Form("A");
  else 
    side.Form("B");
  
  if (detector<36)
    inout.Form("IROC");
  else 
    inout.Form("OROC");

  TString name;
  name.Form("sector_%s%d_%s", side.Data(), sector, inout.Data());
  
  SetName(name);
  SetTitle(Form("%s (detector %d)",name.Data(), detector));

  fTimeStart = timeStart;
  fTimeStop  = timeStop;
  
  Float_t yRange   = 45;
  Int_t nPadRows   = 96;
  
  if (TString(name).Contains("IROC")) {
    yRange = 25;
    nPadRows = 63;
  }
  
  // 1 bin for each 0.5 cm
  Int_t nBinsY = Int_t(4*yRange);

  // TODO do NOT attach to the directory!
  fhDigits = new TH3F("fhDigits", Form("signal distribution;row;pad;time", name.Data()), nPadRows, -0.5, -0.5 + nPadRows, 120, -0.5, 119.5, 100, 0, 1200);
  fhSignal = new TH1F("fhSignal", "fhSignal", 200, 0, 2000);
  
  fDigitTree = new TNtuple("fDigitTree", "fDigitTree", "row:pad:time:signal");
}

//____________________________________________________________________
AliTPCRawHistograms::AliTPCRawHistograms(const AliTPCRawHistograms& c) : TNamed(c)
{
  // copy constructor

  ((AliTPCRawHistograms &)c).Copy(*this);
}

//____________________________________________________________________
AliTPCRawHistograms::~AliTPCRawHistograms()
{
  //
  // destructor
  //

}

//____________________________________________________________________
AliTPCRawHistograms &AliTPCRawHistograms::operator=(const AliTPCRawHistograms &c)
{
  // assigment operator

  if (this != &c)
    ((AliTPCRawHistograms &) c).Copy(*this);

  return *this;
}


//____________________________________________________________________
Long64_t AliTPCRawHistograms::Merge(TCollection* list)
{
  // Merge a list of AliTPCRawHistograms objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

   Int_t count = 0;

/*
  // collections of measured and generated histograms
  TList* collectionQmaxVsRow     = new TList;
  TList* collectionQtotVsRow	 = new TList;
  TList* collectionSigmaYVsRow	 = new TList;
  TList* collectionSigmaZVsRow	 = new TList;
		   			
  TList* collectionQmaxProfileYVsRow    = new TList;
  TList* collectionQtotProfileYVsRow    = new TList;
  TList* collectionSigmaYProfileYVsRow  = new TList;
  TList* collectionSigmaZProfileYVsRow  = new TList;

  TList* collectionQmaxProfileZVsRow    = new TList;
  TList* collectionQtotProfileZVsRow    = new TList;
  TList* collectionSigmaYProfileZVsRow  = new TList;
  TList* collectionSigmaZProfileZVsRow  = new TList;

  TList* collectionQtotVsTime  = new TList;
  TList* collectionQmaxVsTime  = new TList;

   while ((obj = iter->Next())) {
    
     AliTPCRawHistograms* entry = dynamic_cast<AliTPCRawHistograms*> (obj);
     if (entry == 0) 
       continue;

     collectionQmaxVsRow          ->Add(entry->fhQmaxVsRow	   );
     collectionQtotVsRow	  ->Add(entry->fhQtotVsRow	   );
     collectionSigmaYVsRow	  ->Add(entry->fhSigmaYVsRow	   );
     collectionSigmaZVsRow	  ->Add(entry->fhSigmaZVsRow	   );
	       		      		       				   
     collectionQmaxProfileYVsRow  ->Add(entry->fhQmaxProfileYVsRow );
     collectionQtotProfileYVsRow  ->Add(entry->fhQtotProfileYVsRow );
     collectionSigmaYProfileYVsRow->Add(entry->fhSigmaYProfileYVsRow);
     collectionSigmaZProfileYVsRow->Add(entry->fhSigmaZProfileYVsRow);

     collectionQmaxProfileZVsRow  ->Add(entry->fhQmaxProfileZVsRow );
     collectionQtotProfileZVsRow  ->Add(entry->fhQtotProfileZVsRow );
     collectionSigmaYProfileZVsRow->Add(entry->fhSigmaYProfileZVsRow);
     collectionSigmaZProfileZVsRow->Add(entry->fhSigmaZProfileZVsRow);

     collectionQtotVsTime->Add(entry->fhQtotVsTime);
     collectionQmaxVsTime->Add(entry->fhQmaxVsTime);


     count++;
   }

   fhQmaxVsRow          ->Merge(collectionQmaxVsRow       );	   
   fhQtotVsRow          ->Merge(collectionQtotVsRow	  );	   
   fhSigmaYVsRow        ->Merge(collectionSigmaYVsRow	  );	   
   fhSigmaZVsRow        ->Merge(collectionSigmaZVsRow	  );	   
   					       		      	     
   fhQmaxProfileYVsRow  ->Merge(collectionQmaxProfileYVsRow  ); 
   fhQtotProfileYVsRow  ->Merge(collectionQtotProfileYVsRow  );
   fhSigmaYProfileYVsRow->Merge(collectionSigmaYProfileYVsRow);
   fhSigmaZProfileYVsRow->Merge(collectionSigmaZProfileYVsRow);

   fhQmaxProfileZVsRow  ->Merge(collectionQmaxProfileZVsRow  ); 
   fhQtotProfileZVsRow  ->Merge(collectionQtotProfileZVsRow  );
   fhSigmaYProfileZVsRow->Merge(collectionSigmaYProfileZVsRow);
   fhSigmaZProfileZVsRow->Merge(collectionSigmaZProfileZVsRow);

   fhQtotVsTime->Merge(collectionQtotVsTime);
   fhQmaxVsTime->Merge(collectionQmaxVsTime);

   delete collectionQmaxVsRow;          
   delete collectionQtotVsRow;  
   delete collectionSigmaYVsRow;	  
   delete collectionSigmaZVsRow;	  
   	       		      	  
   delete collectionQmaxProfileYVsRow;  
   delete collectionQtotProfileYVsRow;  
   delete collectionSigmaYProfileYVsRow;
   delete collectionSigmaZProfileYVsRow;

   delete collectionQmaxProfileZVsRow;  
   delete collectionQtotProfileZVsRow;  
   delete collectionSigmaYProfileZVsRow;
   delete collectionSigmaZProfileZVsRow;

   delete collectionQtotVsTime;
   delete collectionQmaxVsTime;*/

  return count+1;
}

//____________________________________________________________________
void AliTPCRawHistograms::FillDigit(AliTPCRawStream* rawStream, Int_t time) 
{
  //
  // Fills the different histograms with the information from a raw digit
  //
  
  Int_t signal = rawStream->GetSignal();
  Int_t row = rawStream->GetRow();
  Int_t pad = rawStream->GetPad();
  Int_t timeBin = rawStream->GetTime();

  if (signal > 120)
    fhDigits->Fill(row, pad, timeBin, signal);
    
  fhSignal->Fill(signal);
  
  //fDigitTree->Fill(row, pad, timeBin, signal);
}

//____________________________________________________________________
void AliTPCRawHistograms::SaveHistograms()
{
  //
  // saves the histograms
  //

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());

  fhDigits           ->Write();
  fhSignal           ->Write();
  fDigitTree->Write();
  			
  gDirectory->cd("../");
}

//____________________________________________________________________
TCanvas* AliTPCRawHistograms::DrawHistograms(const Char_t* opt) {
  //
  // Draws some histograms and save the canvas as eps and gif file.
  //  

  TCanvas* c = new TCanvas(Form("plots_%s",fName.Data()), fName.Data(), 1200, 1000);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetPadLeftMargin(0.05);

  c->Divide(3,3);

  c->Draw();  

  c->cd(1);
  
  // this is not really a nice way to do it...
  c->GetPad(1)->Delete();
  
  TLatex* name = new TLatex(0.1,0.8,fName.Data());
  name->SetTextSize(0.02);
  name->DrawClone();

  c->cd(2);
  fhDigits->Draw();

  c->cd(3);
  fhSignal->Draw();         
  			
  return c;
}
