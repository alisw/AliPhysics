/* $Id$ */

// This class contains a number of histograms for diagnostics of a TPC
// read out chamber from the reconstructed clusters.
//
// TODO:
//  
//
//

#include "AliTPCClusterHistograms.h"

#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TProfile2D.h>
#include <TLatex.h>

#include <../TPC/AliTPCclusterMI.h>
#include <AliLog.h>


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
  fhSigmaZProfileYVsRow(0),
  fhQmaxProfileZVsRow(0), 
  fhQtotProfileZVsRow(0),
  fhSigmaYProfileZVsRow(0),
  fhSigmaZProfileZVsRow(0),
  fhQtotVsTime(0),  
  fhQmaxVsTime(0),
  fIsIROC(kFALSE),
  fEdgeSuppression(kFALSE)
{
  // default constructor
}

//____________________________________________________________________
AliTPCClusterHistograms::AliTPCClusterHistograms(Int_t detector, const Char_t* comment, Int_t timeStart, Int_t timeStop, Bool_t edgeSuppression)
  : TNamed(),
  fhQmaxVsRow(0),          
  fhQtotVsRow(0),          
  fhSigmaYVsRow(0),        
  fhSigmaZVsRow(0),          			
  fhQmaxProfileYVsRow(0), 
  fhQtotProfileYVsRow(0),
  fhSigmaYProfileYVsRow(0),
  fhSigmaZProfileYVsRow(0),
  fhQmaxProfileZVsRow(0), 
  fhQtotProfileZVsRow(0),
  fhSigmaYProfileZVsRow(0),
  fhSigmaZProfileZVsRow(0),
  fhQtotVsTime(0),  
  fhQmaxVsTime(0),
  fIsIROC(kFALSE),
  fEdgeSuppression(edgeSuppression)
{
  // constructor 
  
  // make name and title
  if (detector < 0 || detector >= 72) {
    AliDebug(AliLog::kError, Form("Detector %d does not exist", detector));
    return;
  }
      
  TString name(FormDetectorName(detector, edgeSuppression, comment));

  if (detector < 36)
    fIsIROC = kTRUE; 
  
  SetName(name);
  SetTitle(Form("%s (detector %d)",name.Data(), detector));

  fTimeStart = timeStart;
  fTimeStop  = timeStop;
  
  #define BINNING_Z 250, 0, 250
  
  Float_t yRange   = 45;
  Int_t nPadRows   = 96;
  
  if (fIsIROC)
  {
    yRange = 25;
    nPadRows = 63;
  }
  
  // 1 bin for each 0.5 cm
  Int_t nBinsY = Int_t(4*yRange);

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //defining histograms and profile plots
  fhQmaxVsRow  = new TH2F("QmaxVsPadRow", "Qmax vs. pad row;Pad row;Qmax", nPadRows+2, -1.5, nPadRows+0.5, 500,  0,  500);
  fhQtotVsRow  = new TH2F("QtotVsPadRow", "Qtot vs. pad row;Pad row;Qtot", nPadRows+2, -1.5, nPadRows+0.5, 400,  0,  4000);
  
  fhSigmaYVsRow = new TH2F("SigmaYVsPadRow", "Sigma Y vs. pad row;Pad row;#sigma_{Y}", nPadRows+2, -1.5, nPadRows+0.5, 100,  0,  0.5);
  fhSigmaZVsRow = new TH2F("SigmaZVsPadRow", "Sigma Z vs. pad row;Pad row;#sigma_{Z}", nPadRows+2, -1.5, nPadRows+0.5, 100,  0,  0.5);
  
  fhQmaxProfileYVsRow = new TProfile2D("MeanQmaxYVsPadRow","Mean Qmax, y vs pad row;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);
  fhQtotProfileYVsRow = new TProfile2D("MeanQtotYVsPadRow","Mean Qtot, y vs pad row;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);
  fhSigmaYProfileYVsRow = new TProfile2D("MeanSigmaYYVsPadRow","Mean Sigma y, y vs pad row;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);
  fhSigmaZProfileYVsRow = new TProfile2D("MeanSigmaZYVsPadRow","Mean Sigma z, y vs pad row;Pad row;y",nPadRows+2, -1.5, nPadRows+0.5, nBinsY, -yRange, yRange);

  fhQmaxProfileZVsRow = new TProfile2D("MeanQmaxZVsPadRow","Mean Qmax, z vs pad row;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);
  fhQtotProfileZVsRow = new TProfile2D("MeanQtotZVsPadRow","Mean Qtot, z vs pad row;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);
  fhSigmaYProfileZVsRow = new TProfile2D("MeanSigmaYZVsPadRow","Mean Sigma y, z vs pad row;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);
  fhSigmaZProfileZVsRow = new TProfile2D("MeanSigmaZZVsPadRow","Mean Sigma z, z vs pad row;Pad row;z",nPadRows+2, -1.5, nPadRows+0.5, BINNING_Z);

  Int_t nTimeBins  = 100;
  
  fhQtotVsTime = new TProfile("MeanQtotVsTime", "Mean Qtot vs. event time stamp; time; Qtot",nTimeBins, fTimeStart, fTimeStop);
  fhQmaxVsTime = new TProfile("MeanQmaxVsTime", "Mean Qmax vs. event time stamp; time; Qmax",nTimeBins, fTimeStart, fTimeStop);
  
  TH1::AddDirectory(oldStatus);
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
  if (fhQmaxProfileZVsRow) {
    delete fhQmaxProfileZVsRow;
    fhQmaxProfileZVsRow = 0;
  }
  if (fhQtotProfileZVsRow) {
    delete fhQtotProfileZVsRow;
    fhQtotProfileZVsRow = 0;
  }
  if (fhSigmaYProfileZVsRow) {
    delete fhSigmaYProfileZVsRow;
    fhSigmaYProfileZVsRow = 0;
  }
  if (fhSigmaZProfileZVsRow) {
    delete fhSigmaZProfileZVsRow;
    fhSigmaZProfileZVsRow = 0;
  }

  if (fhQtotVsTime) {
    delete fhQtotVsTime;
    fhQtotVsTime = 0;
  }
  if (fhQmaxVsTime) {
    delete fhQmaxVsTime;
    fhQmaxVsTime = 0;
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
const char* AliTPCClusterHistograms::FormDetectorName(Int_t detector, Bool_t edgeSuppression, const char* comment)
{
  //
  // creates a readable name from the detector number
  //   
  
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

  if (edgeSuppression)
    name += "_noedge";
  
  if (comment)
    name += comment;

  return name; 
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

  TList* collectionQmaxProfileZVsRow    = new TList;
  TList* collectionQtotProfileZVsRow    = new TList;
  TList* collectionSigmaYProfileZVsRow  = new TList;
  TList* collectionSigmaZProfileZVsRow  = new TList;

  TList* collectionQtotVsTime  = new TList;
  TList* collectionQmaxVsTime  = new TList;

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
   delete collectionQmaxVsTime;

  return count+1;
}

//____________________________________________________________________
void AliTPCClusterHistograms::FillCluster(AliTPCclusterMI* cluster, Int_t time) {
  //
  // Fills the different histograms with the information from the cluster.
  //

  Int_t padRow =   cluster->GetRow(); 
  Float_t qMax =   cluster->GetMax();
  Float_t qTot =   cluster->GetQ();
  Float_t sigmaY = cluster->GetSigmaY2();
  Float_t sigmaZ = cluster->GetSigmaZ2();
  Float_t y      = cluster->GetY();
  Float_t z      = cluster->GetZ();

  if (qMax<=0) {
    printf(Form("\n WARNING: Hi Marian! How can we have Qmax = %f ??? \n \n", qMax));
    return;
  }
  if (qTot<=0) {
    printf(Form("\n WARNING: Hi Marian! How can we have Qtot = %f ??? \n \n ", qTot));
    return;
  } 
  
  if (fEdgeSuppression)
  {
    Float_t limit = 0;
    if (fIsIROC)
    {
      limit = 12 + padRow * (20.0 - 12.0) / 63; 
    }
    else
      limit = 16 + padRow * (36.0 - 16.0) / 96;
      
    if (TMath::Abs(y) > limit)
      return;  
  }

  fhQmaxVsRow           ->Fill(padRow, qMax);
  fhQtotVsRow           ->Fill(padRow, qTot);
  			
  fhSigmaYVsRow         ->Fill(padRow, sigmaY);
  fhSigmaZVsRow         ->Fill(padRow, sigmaZ);
  			
  fhQmaxProfileYVsRow   ->Fill(padRow, y, qMax);
  fhQtotProfileYVsRow   ->Fill(padRow, y, qTot); 
  fhSigmaYProfileYVsRow ->Fill(padRow, y, sigmaY);
  fhSigmaZProfileYVsRow ->Fill(padRow, y, sigmaZ);

  fhQmaxProfileZVsRow   ->Fill(padRow, z, qMax);
  fhQtotProfileZVsRow   ->Fill(padRow, z, qTot); 
  fhSigmaYProfileZVsRow ->Fill(padRow, z, sigmaY);
  fhSigmaZProfileZVsRow ->Fill(padRow, z, sigmaZ);

  if (time>0 & fTimeStart>0 & fTimeStop>0 & time>fTimeStart) {
    //Float_t timeFraction = (time - fTimeStart)/(fTimeStop-fTimeStart); 

    fhQtotVsTime->Fill(time,qTot);
    fhQmaxVsTime->Fill(time,qMax);
  }

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

  fhQmaxProfileZVsRow   ->Write();
  fhQtotProfileZVsRow   ->Write();
  fhSigmaYProfileZVsRow ->Write();
  fhSigmaZProfileZVsRow ->Write();

  TProfile* profileQmaxVsRow = fhQmaxVsRow->ProfileX("MeanQmaxVsRow");
  TProfile* profileQtotVsRow = fhQtotVsRow->ProfileX("MeanQtotVsRow");

  profileQmaxVsRow->SetTitle("Mean Qmax vs. pad row; Pad row; Mean Qmax");
  profileQtotVsRow->SetTitle("Mean Qtot vs. pad row; Pad row; Mean Qmax");

  profileQmaxVsRow->Write();
  profileQtotVsRow->Write();

  if (fhQtotVsTime->GetEntries()>0)
    fhQtotVsTime->Write();

  if (fhQmaxVsTime->GetEntries()>0)
    fhQmaxVsTime->Write();

  gDirectory->cd("../");

}

//____________________________________________________________________
TCanvas* AliTPCClusterHistograms::DrawHistograms(const Char_t* /*opt*/) {
  //
  // Draws some histograms and save the canvas as eps and gif file.
  //  

  TCanvas* c = new TCanvas(fName.Data(), fName.Data(), 1200, 1000);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetPadLeftMargin(0.1);

  c->Divide(3,3);

  c->Draw();  

  c->cd(1);
  
  // this is not really a nice way to do it...
  c->GetPad(1)->Delete();
  
  TLatex* name = new TLatex(0.1,0.8,fName.Data());
  name->SetTextSize(0.02);
  name->DrawClone();

  c->cd(2);
  fhQmaxVsRow->Draw("colz");

  c->cd(3);
  fhQtotVsRow->Draw("colz");         
  			
  c->cd(4);
  fhQmaxProfileYVsRow   ->Draw("colz");
  c->cd(5);
  fhQtotProfileYVsRow   ->Draw("colz");
  c->cd(6);
  fhQmaxProfileZVsRow   ->Draw("colz");
  c->cd(7);
  fhQtotProfileZVsRow   ->Draw("colz");
  c->cd(8);
  fhSigmaYVsRow         ->Draw("colz");

  c->cd(9);
  fhSigmaZVsRow         ->Draw("colz");
  			
  //fhSigmaYProfileYVsRow ->Draw("colz");
  //fhSigmaZProfileYVsRow ->Draw("colz");

  //fhSigmaYProfileZVsRow ->Draw("colz");
  //fhSigmaZProfileZVsRow ->Draw("colz");
  return c;
}
