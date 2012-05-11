#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

#include "AliHighPtDeDxBase.h"

using namespace std;


class DeDxFitInfo : public TObject
{
 public:
  
  DeDxFitInfo();
  void Print(Option_t* option="") const;

  Double_t MIP;
  Double_t plateau;

  Int_t optionDeDx;
  Int_t nDeDxPar;
  Double_t parDeDx[5];

  Int_t optionSigma;
  Int_t nSigmaPar;
  Double_t parSigma[5];

  TString calibFileName;

  ClassDef(DeDxFitInfo, 1);    // Help class
};


//_____________________________________________________________________________
ClassImp(DeDxFitInfo)

DeDxFitInfo::DeDxFitInfo():
TObject(),
  MIP(0),
  plateau(0),
  optionDeDx(-1),
  nDeDxPar(-1),
  optionSigma(-1),
  nSigmaPar(-1),
  calibFileName("")
{
  // default constructor
  for(Int_t i = 0; i < 5; i++) {
    parDeDx[i]  = 0;
    parSigma[i] = 0;
  }
}

//_________________________________________________________
void DeDxFitInfo::Print(Option_t* option) const
{
  if(option) 
    cout << "Option: " << option << endl;

  cout << ClassName() << " : " << GetName() << endl  
       << "MIP: " << MIP << endl
       << "Plateau: " << plateau << endl
       << "OptionDeDx: " << optionDeDx << endl
       << "nDeDxPar: " << nDeDxPar << endl;
  for(Int_t i = 0; i < nDeDxPar; i++) {
    
    cout << "parDeDx[" << i << "] = " << parDeDx[i] << endl;
  }
  cout << "OptionSigma: " << optionSigma << endl
       << "nSigmaPar: " << nSigmaPar << endl;
  for(Int_t i = 0; i < nSigmaPar; i++) {
    
    cout << "parSigma[" << i << "] = " << parSigma[i] << endl;
  }
  
  if(calibFileName.IsNull()) {
    cout << "No eta calibration file." << endl; 
  } else {
    cout << "Eta calibration file: " << calibFileName.Data() << endl; 
  }
}

//_____________________________________________________________________________


AliHighPtDeDxBase* GetObject(TFile* file, Int_t filter, Bool_t phiCut, 
			     Int_t run, const Char_t* baseName="filter",
			     const Char_t* endName=0);

TFile* FindFileFresh(const Char_t* fileName);
TFile* FindFile(const Char_t* fileName);
void CutHistogram(TH1* hist, Double_t xMin, Double_t xMax);
void SetHistError(TH1* hist, Double_t error);
void CreateDir(const Char_t* dirName);

//___________________________________________________________________________
AliHighPtDeDxBase* GetObject(TFile* file, Int_t filter, Bool_t phiCut, Int_t run,
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
