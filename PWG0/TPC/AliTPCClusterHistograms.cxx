/* $Id$ */


#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TProfile2D.h>

#include <../TPC/AliTPCclusterMI.h>

#include <AliLog.h>

#include "AliTPCClusterHistograms.h"

//____________________________________________________________________
ClassImp(AliTPCClusterHistograms)

//____________________________________________________________________
AliTPCClusterHistograms::AliTPCClusterHistograms() 
  : TNamed(),
  fhQmaxVsRow(0),          
  fhQtotVsRow(0),          
  fhSigmaYVsRow(0),        
  fhSigmaZVsRow(0),          			
  fhQmaxProfileYVsRow(0), 
  fhQtotProfileYVsRow(0),
  fhSigmaYProfileYVsRow(0),
  fhSigmaZProfileYVsRow(0)
{
  // default constructor
}

//____________________________________________________________________
AliTPCClusterHistograms::AliTPCClusterHistograms(const Char_t* name, const Char_t* title) 
  : TNamed(name, title),
  fhQmaxVsRow(0),          
  fhQtotVsRow(0),          
  fhSigmaYVsRow(0),        
  fhSigmaZVsRow(0),          			
  fhQmaxProfileYVsRow(0), 
  fhQtotProfileYVsRow(0),
  fhSigmaYProfileYVsRow(0),
  fhSigmaZProfileYVsRow(0)
{
  // constructor initializing tnamed

  #define BINNING_Y 100, -50, 50

  fhQmaxVsRow  = new TH2F("QmaxVsPadRow", "Qmax vs. pad row;Pad row;Qmax", 91, -0.5, 90.5, 301, -0.5, 300.5);
  fhQtotVsRow  = new TH2F("QtotVsPadRow", "Qtot vs. pad row;Pad row;Qtot", 91, -0.5, 90.5, 100,  0,  1000);
  
  fhSigmaYVsRow = new TH2F("SigmaYVsPadRow", "Sigma Y vs. pad row;Pad row;#sigma_{Y}", 91, -0.5, 90.5, 100,  0,  0.5);
  fhSigmaZVsRow = new TH2F("SigmaZVsPadRow", "Sigma Z vs. pad row;Pad row;#sigma_{Z}", 91, -0.5, 90.5, 100,  0,  0.5);
  
  fhQmaxProfileYVsRow = new TProfile2D("QmaxMeanYVsPadRow","Mean Qmax, y vs pad row;Pad row;y",91,-0.5,90.5, BINNING_Y);
  fhQtotProfileYVsRow = new TProfile2D("QtotMeanYVsPadRow","Mean Qtot, y vs pad row;Pad row;y",91,-0.5,90.5, BINNING_Y);
  fhSigmaYProfileYVsRow = new TProfile2D("SigmaYMeanYVsPadRow","Mean Sigma y, x vs pad row;Pad row;y",91,-0.5,90.5, BINNING_Y);
  fhSigmaZProfileYVsRow = new TProfile2D("SigmaZMeanYVsPadRow","Mean Sigma y, x vs pad row;Pad row;y",91,-0.5,90.5, BINNING_Y);

}

//____________________________________________________________________
AliTPCClusterHistograms::AliTPCClusterHistograms(const AliTPCClusterHistograms& c) : TNamed(c)
{
  // copy constructor
  ((AliTPCClusterHistograms &)c).Copy(*this);
}

//____________________________________________________________________
AliTPCClusterHistograms::~AliTPCClusterHistograms()
{
  //
  // destructor
  //

  if (fhQmaxVsRow) {
    delete fhQmaxVsRow;
    fhQmaxVsRow = 0;
  }
  if (fhQtotVsRow) {
    delete fhQtotVsRow;
    fhQtotVsRow = 0; 
  }
  if (fhSigmaYVsRow) {
    delete fhSigmaYVsRow;
    fhSigmaYVsRow = 0;
  } 
  if (fhSigmaZVsRow) {
    delete fhSigmaZVsRow;
    fhSigmaZVsRow = 0; 
  }
  if (fhQmaxProfileYVsRow) {
    delete fhQmaxProfileYVsRow;
    fhQmaxProfileYVsRow = 0;
  }
  if (fhQtotProfileYVsRow) {
    delete fhQtotProfileYVsRow;
    fhQtotProfileYVsRow = 0;
  }
  if (fhSigmaYProfileYVsRow) {
    delete fhSigmaYProfileYVsRow;
    fhSigmaYProfileYVsRow = 0;
  }
  if (fhSigmaZProfileYVsRow) {
    delete fhSigmaZProfileYVsRow;
    fhSigmaZProfileYVsRow = 0;
  }

}

//____________________________________________________________________
AliTPCClusterHistograms &AliTPCClusterHistograms::operator=(const AliTPCClusterHistograms &c)
{
  // assigment operator

  if (this != &c)
    ((AliTPCClusterHistograms &) c).Copy(*this);

  return *this;
}


//____________________________________________________________________
Long64_t AliTPCClusterHistograms::Merge(TCollection* list)
{
  // Merge a list of AliTPCClusterHistograms objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of measured and generated histograms
  TList* collectionQmaxVsRow     = new TList;
  TList* collectionQtotVsRow	 = new TList;
  TList* collectionSigmaYVsRow	 = new TList;
  TList* collectionSigmaZVsRow	 = new TList;
		   			
  TList* collectionQmaxProfileYVsRow    = new TList;
  TList* collectionQtotProfileYVsRow    = new TList;
  TList* collectionSigmaYProfileYVsRow  = new TList;
  TList* collectionSigmaZProfileYVsRow  = new TList;

   Int_t count = 0;
   while ((obj = iter->Next())) {
    
     AliTPCClusterHistograms* entry = dynamic_cast<AliTPCClusterHistograms*> (obj);
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

   delete collectionQmaxVsRow;          
   delete collectionQtotVsRow;  
   delete collectionSigmaYVsRow;	  
   delete collectionSigmaZVsRow;	  
   	       		      	  
   delete collectionQmaxProfileYVsRow;  
   delete collectionQtotProfileYVsRow;  
   delete collectionSigmaYProfileYVsRow;
   delete collectionSigmaZProfileYVsRow;

  return count+1;
}

//____________________________________________________________________
void AliTPCClusterHistograms::FillCluster(AliTPCclusterMI* cluster) {
  //
  //
  //

  Int_t padRow =   cluster->GetRow(); 
  Float_t qMax =   cluster->GetMax();
  Float_t qTot =   cluster->GetQ();
  Float_t sigmaY = cluster->GetSigmaY2();
  Float_t sigmaZ = cluster->GetSigmaZ2();
  Float_t y      = cluster->GetY();

  fhQmaxVsRow           ->Fill(padRow, qMax);
  fhQtotVsRow           ->Fill(padRow, qTot);
  			
  fhSigmaYVsRow         ->Fill(padRow, sigmaY);
  fhSigmaZVsRow         ->Fill(padRow, sigmaZ);
  			
  fhQmaxProfileYVsRow   ->Fill(padRow, y, qMax);
  fhQtotProfileYVsRow   ->Fill(padRow, y, qTot); 
  fhSigmaYProfileYVsRow ->Fill(padRow, y, sigmaY);
  fhSigmaZProfileYVsRow ->Fill(padRow, y, sigmaZ);
}


//____________________________________________________________________
void AliTPCClusterHistograms::SaveHistograms()
{
  //
  // saves the histograms
  //

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());

  fhQmaxVsRow           ->Write();
  fhQtotVsRow           ->Write();
  			
  fhSigmaYVsRow         ->Write();
  fhSigmaZVsRow         ->Write();
  			
  fhQmaxProfileYVsRow   ->Write();
  fhQtotProfileYVsRow   ->Write();
  fhSigmaYProfileYVsRow ->Write();
  fhSigmaZProfileYVsRow ->Write();

  gDirectory->cd("../");

}

