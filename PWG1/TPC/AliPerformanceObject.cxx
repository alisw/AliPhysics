/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERF, All rights reserved. *
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

//------------------------------------------------------------------------------
// Implementation of abstract AliPerformanceObject class. It keeps information from 
// comparison of reconstructed and MC particle tracks. 
//
// Author: J.Otwinowski 14/04/2008 
// Changes by M.Knichel 15/10/2010
//------------------------------------------------------------------------------

#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TPostScript.h"
#include "TList.h"
#include "TMath.h"

#include "AliLog.h" 
#include "AliESDVertex.h" 
#include "AliPerformanceObject.h" 

using namespace std;

ClassImp(AliPerformanceObject)

//_____________________________________________________________________________
AliPerformanceObject::AliPerformanceObject():
  TNamed("AliPerformanceObject","AliPerformanceObject"),
  fAnalysisMode(-1),
  fRunNumber(-1),
  fHptGenerator(kFALSE),
  fTriggerClass(0),
  fUseTrackVertex(kFALSE),
  fHighMultiplicity(kFALSE),
  fUseKinkDaughters(kTRUE),
  fUseCentralityBin(0)
{
  // constructor
}

//_____________________________________________________________________________
AliPerformanceObject::AliPerformanceObject(const char* name, const char* title, Int_t run, Bool_t highMult):
  TNamed(name,title),
  fAnalysisMode(-1),
  fRunNumber(run),
  fHptGenerator(kFALSE),
  fTriggerClass(0),
  fUseTrackVertex(kFALSE),
  fHighMultiplicity(highMult),
  fUseKinkDaughters(kTRUE),
  fUseCentralityBin(0)
{
  // constructor
}

//_____________________________________________________________________________
AliPerformanceObject::~AliPerformanceObject(){
  // destructor 
}

//_____________________________________________________________________________
void AliPerformanceObject::PrintHisto(Bool_t logz, Char_t * outFileName) {
  // draw all histograms from the folder 
  // and store them in the output *.ps file
 
  // use this folder
  TFolder *folder = this->GetAnalysisFolder();
  if (!folder) {
     AliDebug(AliLog::kError, "folder not available");
     return;
  } 

  TCanvas *can = new TCanvas("can");
  can->Divide(2,2);

  char fname[256];
  const char* suffix=".ps"; 

  if(outFileName) snprintf(fname,256,"%s",outFileName);
  else snprintf(fname,256,"%s%s",folder->GetName(),suffix);
  
  TPostScript *ps = new TPostScript(fname,112);
  Printf("Histograms are stored in %s", fname); 
  TIter iter(folder->GetListOfFolders());

  TH1 *obj = 0;
  Int_t count = 0;
  Int_t pad_count = 0;
  while ((obj = (TH1*)iter()) !=0) {

    TString name(obj->ClassName());
 

    // 4 figures per page
    if((count%4) == 0) {
      pad_count = 0;
      ps->NewPage();
    }

    pad_count++; 
    can->cd(pad_count);
    
    if(obj->TestBit(TH1::kLogX)) 
       gPad->SetLogx(1);
    else    
       gPad->SetLogx(0);

    if (obj->GetYaxis() && obj->GetZaxis()) {
      if(logz) gPad->SetLogz();
      if ( name.CompareTo("TH3D") )
	obj->Draw("colz");
    }
    else { 
      obj->SetMarkerStyle(24);
      obj->SetMarkerSize(1.0);
      if ( name.CompareTo("TH3D") )
	obj->Draw();
    }

    if ((pad_count%4) == 0)  { 
      can->Update();
    }

  count++;
  }
  ps->Close();
}
 

//_____________________________________________________________________________
Double_t * AliPerformanceObject::CreateLogAxis(Int_t nbins, Double_t xmin, Double_t xmax) {
  // retun pointer to the array with log axis
  // it is user responsibility to delete the array
 
  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  
  Double_t *xbins =  new Double_t[nbins+1];

  xbins[0] = xmin;
  for (Int_t i=1;i<=nbins;i++) {
    xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
  }

return xbins;
}

//_____________________________________________________________________________
void AliPerformanceObject::InitHighMult() {

  fHighMultiplicity = kTRUE;
  Init();

}


//_____________________________________________________________________________
void AliPerformanceObject::AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse* hSparse, Int_t xDim, TString* selString) 
{
  TH1 *h1=0;  
  TString name = "h_tpc_" + nameSparse + '_';
  if (selString) { name += *selString + '_'; }
  name.ToLower();
  name += xDim;
  TString title = hSparse->GetAxis(xDim)->GetTitle();  
  if (selString) { title += " (" + *selString + ")"; }
  h1 = hSparse->Projection(xDim);
  h1->SetName(name.Data());
  h1->GetXaxis()->SetTitle(hSparse->GetAxis(xDim)->GetTitle());
  h1->SetTitle(title.Data());  
  aFolderObj->Add(h1);
}


//_____________________________________________________________________________
void AliPerformanceObject::AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse *hSparse, Int_t yDim, Int_t xDim, TString* selString)
{
  TH2 *h2=0;  
  TString name = "h_tpc_" + nameSparse + '_';
  if (selString) { name += *selString + '_'; }
  name.ToLower();
  name += yDim;
  name += '_';
  name += xDim;
  TString title = hSparse->GetAxis(yDim)->GetTitle();
  title += " vs ";
  title += hSparse->GetAxis(xDim)->GetTitle();
  if (selString) { title += " (" + *selString + ")"; }  
  h2 = hSparse->Projection(yDim,xDim);
  h2->SetName(name.Data());
  h2->GetXaxis()->SetTitle(hSparse->GetAxis(xDim)->GetTitle());
  h2->GetYaxis()->SetTitle(hSparse->GetAxis(yDim)->GetTitle());
  h2->SetTitle(title.Data());  
  aFolderObj->Add(h2);
}


//_____________________________________________________________________________
void AliPerformanceObject::AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse *hSparse, Int_t xDim, Int_t yDim, Int_t zDim, TString* selString)
{
  TH3 *h3=0;  
  TString name = "h_tpc_" + nameSparse + '_';
  if (selString) { name += *selString + '_'; }
  name.ToLower();
  name += xDim;
  name += '_';
  name += yDim;
  name += '_';
  name += zDim;
  TString title = hSparse->GetAxis(xDim)->GetTitle();
  title += " vs ";
  title += hSparse->GetAxis(yDim)->GetTitle();
  title += " vs ";
  title += hSparse->GetAxis(zDim)->GetTitle();
  if (selString) { title += " (" + *selString + ")"; }
  h3 = hSparse->Projection(xDim,yDim,zDim);
  h3->SetName(name.Data());
  h3->GetXaxis()->SetTitle(hSparse->GetAxis(xDim)->GetTitle());
  h3->GetYaxis()->SetTitle(hSparse->GetAxis(yDim)->GetTitle());
  h3->GetZaxis()->SetTitle(hSparse->GetAxis(zDim)->GetTitle());
  h3->SetTitle(title.Data());  
  aFolderObj->Add(h3);
}
