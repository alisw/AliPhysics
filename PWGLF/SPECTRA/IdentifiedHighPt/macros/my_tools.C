#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

#include "AliHighPtDeDxBase.h"

using namespace std;



//_____________________________________________________________________________


AliHighPtDeDxBase* GetObject(TFile* file, Int_t filter, Bool_t phiCut, 
			     Int_t run, Bool_t etaAbs, 
			     Int_t etaLow, Int_t etaHigh, const Char_t* baseName="filter",
			     const Char_t* endName=0);

TFile* FindFileFresh(const Char_t* fileName);
TFile* FindFile(const Char_t* fileName);
void CutHistogram(TH1* hist, Double_t xMin, Double_t xMax);
void SetHistError(TH1* hist, Double_t error);
void CreateDir(const Char_t* dirName);

//___________________________________________________________________________
AliHighPtDeDxBase* GetObject(TFile* file, Int_t filter, Bool_t phiCut, Int_t run, Bool_t etaAbs, 
			     Int_t etaLow, Int_t etaHigh,
			     const Char_t* baseName, const Char_t* endName)
{
  TString objectName(baseName);
  if(filter>0)
  objectName += filter;
  if(phiCut)
    objectName += "phicut";
  if(run) {
    objectName += "_";
    objectName += run;
  }
  if(etaAbs) {
    objectName += "etaabs";
    objectName += etaLow;
    objectName += etaHigh;
  }

  if(endName)
    objectName += endName;

  cout << "Getting object: " << objectName.Data() << endl;

  return (AliHighPtDeDxBase*)(file->Get(objectName.Data()));
}

//______________________________________________________________________
TFile* FindFileFresh(const Char_t* fileName)
{
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    file->Close();
    delete file;
  }

  file = TFile::Open(fileName, "READ");

  if(!file)
    cout << "File : " << fileName << " was not found" << endl;

  return file;
}

//______________________________________________________________________
TFile* FindFile(const Char_t* fileName)
{
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    return file;
  }

  file = TFile::Open(fileName, "READ");

  if(!file)
    cout << "File : " << fileName << " was not found" << endl;

  return file;
}

//______________________________________________________________________
void CutHistogram(TH1* hist, Double_t xMin, Double_t xMax)
{
  const Int_t n = hist->GetNbinsX();
  
  for(Int_t bin = 1; bin <= n; bin++) {
    
    Float_t x = hist->GetXaxis()->GetBinCenter(bin);
    if(x < xMin) {
      hist->SetBinContent(bin, 0);
      hist->SetBinError(bin, 0);
    } else if(x > xMax) {
      hist->SetBinContent(bin, 0);
      hist->SetBinError(bin, 0);
    }

  }
}

//______________________________________________________________________
void SetHistError(TH1* hist, Double_t error)
{
  const Int_t n = hist->GetNbinsX();
  
  for(Int_t bin = 1; bin <= n; bin++) {
    
    //    Float_t x = hist->GetXaxis()->GetBinCenter(bin);
    hist->SetBinError(bin, error);
  }
}

//______________________________________________________________________
void CreateDir(const Char_t* dirName)
{
  TString pwd(gSystem->pwd());
  gSystem->cd(pwd.Data());
  
  if(gSystem->cd(dirName)) {
    gSystem->cd(pwd.Data());
  } else {
    gSystem->mkdir(dirName, kTRUE); // kTRUE means recursive
  }
}
