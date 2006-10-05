/* $Id$ */

// ------------------------------------------------------
//
// Class to handle corrections.
//
// ------------------------------------------------------
//

#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>

#include <AliLog.h>

#include "AliCorrectionMatrix.h"

//____________________________________________________________________
ClassImp(AliCorrectionMatrix)

//____________________________________________________________________
AliCorrectionMatrix::AliCorrectionMatrix() : TNamed(),
  fhMeas(0),
  fhGene(0),
  fhCorr(0)
{
  // default constructor
}

//____________________________________________________________________
AliCorrectionMatrix::AliCorrectionMatrix(const Char_t* name, const Char_t* title) : TNamed(name, title),
  fhMeas(0),
  fhGene(0),
  fhCorr(0)
{
  // constructor initializing tnamed
}

//____________________________________________________________________
AliCorrectionMatrix::AliCorrectionMatrix(const AliCorrectionMatrix& c) : TNamed(c),
  fhMeas(0),
  fhGene(0),
  fhCorr(0)
{
  // copy constructor
  ((AliCorrectionMatrix &)c).Copy(*this);
}

//____________________________________________________________________
AliCorrectionMatrix::~AliCorrectionMatrix()
{
  //
  // destructor
  //

  if (fhMeas)
  {
    delete fhMeas;
    fhMeas = 0;
  }

  if (fhGene)
  {
    delete fhGene;
    fhGene = 0;
  }

  if (fhCorr)
  {
    delete fhCorr;
    fhCorr = 0;
  }
}

//____________________________________________________________________
AliCorrectionMatrix &AliCorrectionMatrix::operator=(const AliCorrectionMatrix &c)
{
  // assigment operator

  if (this != &c)
    ((AliCorrectionMatrix &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
void AliCorrectionMatrix::Copy(TObject& c) const
{
  // copy function

  AliCorrectionMatrix& target = (AliCorrectionMatrix &) c;

  if (fhMeas)
    target.fhMeas = dynamic_cast<TH1*> (fhMeas->Clone());

  if (fhGene)
    target.fhGene = dynamic_cast<TH1*> (fhGene->Clone());

  if (fhCorr)
    target.fhCorr = dynamic_cast<TH1*> (fhCorr->Clone());
}

//________________________________________________________________________
void AliCorrectionMatrix::SetAxisTitles(const Char_t* titleX, const Char_t* titleY, const Char_t* titleZ)
{
  //
  // method for setting the axis titles of the histograms
  //

  fhMeas ->SetXTitle(titleX);  fhMeas ->SetYTitle(titleY);  fhMeas ->SetZTitle(titleZ);
  fhGene ->SetXTitle(titleX);  fhGene ->SetYTitle(titleY);  fhGene ->SetZTitle(titleZ);
  fhCorr ->SetXTitle(titleX);  fhCorr ->SetYTitle(titleY);  fhCorr ->SetZTitle(titleZ);
}

//____________________________________________________________________
Long64_t AliCorrectionMatrix::Merge(TCollection* list)
{
  // Merge a list of AliCorrectionMatrix objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of measured and generated histograms
  TList* collectionMeas = new TList;
  TList* collectionGene = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliCorrectionMatrix* entry = dynamic_cast<AliCorrectionMatrix*> (obj);
    if (entry == 0) 
      continue;

    collectionMeas->Add(entry->GetMeasuredHistogram());
    collectionGene->Add(entry->GetGeneratedHistogram());

    count++;
  }
  fhMeas->Merge(collectionMeas);
  fhGene->Merge(collectionGene);

  delete collectionMeas;
  delete collectionGene;

  return count+1;
}

//____________________________________________________________________
void AliCorrectionMatrix::Divide()
{
  //
  // divide the histograms to get the correction
  // 

  if (!fhMeas || !fhGene)
    return;

  fhCorr->Divide(fhGene, fhMeas, 1, 1, "B");

  Int_t emptyBins = 0;
  for (Int_t x=1; x<=fhCorr->GetNbinsX(); ++x)
    for (Int_t y=1; y<=fhCorr->GetNbinsY(); ++y)
      for (Int_t z=1; z<=fhCorr->GetNbinsZ(); ++z)
        if (fhCorr->GetBinContent(x, y, z) == 0)
          ++emptyBins;

  if (emptyBins > 0)
    printf("INFO: In %s we have %d empty bins (of %d) in the correction map\n", GetName(), emptyBins, fhCorr->GetNbinsX() * fhCorr->GetNbinsY() * fhCorr->GetNbinsZ());
}

//____________________________________________________________________
Bool_t AliCorrectionMatrix::LoadHistograms(const Char_t* fileName, const Char_t* dir)
{
  //
  // loads the histograms from a file
  //

  TFile* fin = TFile::Open(fileName);

  if(!fin) {
    //Info("LoadHistograms",Form(" %s file does not exist",fileName));
    return kFALSE;
  }
  
  if(fhGene)  {delete fhGene;  fhGene=0;}
  if(fhCorr)  {delete fhCorr;  fhCorr=0;}
  if(fhMeas)  {delete fhMeas;  fhMeas=0;}
  
  fhMeas  = dynamic_cast<TH1*> (fin->Get(Form("%s/meas_%s", dir,GetName())));
  if(!fhMeas)  Info("LoadHistograms","No meas. (%s) hist available",Form("%s/meas_%s", dir,GetName()));

  fhGene  = dynamic_cast<TH1*> (fin->Get(Form("%s/gene_%s",dir, GetName())));
  if(!fhGene)  Info("LoadHistograms","No gene. (%s) hist available",Form("%s/gene_%s",dir, GetName()));

  fhCorr  = dynamic_cast<TH1*> (fin->Get(Form("%s/corr_%s",dir, GetName())));
  if(!fhCorr) {
    Info("LoadHistograms","No corr.(%s) hist available",Form("%s/corr_%s",dir, GetName()));
    return kFALSE;
  }
      
  return kTRUE;
}

//____________________________________________________________________
void AliCorrectionMatrix::SaveHistograms()
{
  //
  // saves the histograms
  //

  if (fhMeas)
    fhMeas ->Write();

  if (fhGene)
    fhGene ->Write();

  if (fhCorr)
    fhCorr->Write();
}

//____________________________________________________________________
void AliCorrectionMatrix::DrawHistograms()
{
  //
  // draws all the four histograms on one TCanvas
  //

  TCanvas* canvas = new TCanvas(Form("correction_%s",fName.Data()), 
				Form("correction_%s",fName.Data()), 800, 800);
  canvas->Divide(2, 2);

  canvas->cd(1);
  if (fhMeas)
    fhMeas->Draw("COLZ");
  
  canvas->cd(2);
  if (fhGene)
    fhGene->Draw("COLZ");

  canvas->cd(3);
  if (fhCorr)
    fhCorr->Draw("COLZ");

  canvas->cd(4);

  // add: draw here the stat. errors of the correction histogram
}

//____________________________________________________________________
void AliCorrectionMatrix::ReduceInformation()
{
  // this function deletes the measured and generated histograms to reduce the amount of data
  // in memory

  if (fhMeas)
  {
    delete fhMeas;
    fhMeas = 0;
  }

  if (fhGene)
  {
    delete fhGene;
    fhGene = 0;
  }
}

//____________________________________________________________________
void AliCorrectionMatrix::Reset(Option_t* option)
{
  // resets the histograms

  if (fhGene)
    fhGene->Reset(option);

  if (fhMeas)
    fhMeas->Reset(option);

  if (fhCorr)
    fhCorr->Reset(option);
}
