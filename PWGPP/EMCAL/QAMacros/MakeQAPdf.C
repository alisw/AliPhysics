#if !defined(__CINT__) || defined(__MAKECINT__)
#include <string.h>
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TText.h"
#include "TError.h"
#endif

///
/// \file MakeQAPdf.C
/// \ingroup EMCALOfflineMacros
/// \brief QA to PDF files
///
/// Make a pdf with all QA plots
///
/// \author Alexandre Shabetai, <Alexandre.Shabetai@subatech.in2p3.fr>  SUBATECH
///

TFile* sourceFile;
TObject* Myobj;
TCanvas* canvasDefault;
TString Mypath;
const char* title;

void recurseOverKeys( TDirectory *target, TString f);

//---------------------------------------------------------------------------------------
void MakeQAPdf(const char* fileName, int outputWidth  = 600,  int outputHeight = 600) 
{

  sourceFile = TFile::Open(fileName);

  canvasDefault = new TCanvas("canvasDefault","testCanvas",outputWidth,outputHeight);
  TString f = TString(fileName).ReplaceAll(".root","");
  
  canvasDefault->Print(Form("%s%s",f.Data(),".pdf["));
  gErrorIgnoreLevel = kInfo + 1;

  // Now actually find all the directories, canvas, and make a pdf..
  recurseOverKeys(sourceFile,f);  

  gPad->Print(Form("%s%s",f.Data(),".pdf]"),Form("Title:%s/%s",Mypath.Data(),title));
  gErrorIgnoreLevel = -1;
  sourceFile->Close(); 

}

//---------------------------------------------------------------------------------------
void recurseOverKeys( TDirectory *target,TString f ) 
{
 
  TString path = TString((char*)strstr( target->GetPath(), ":" ));
  path.Remove(0, 2);
 
  gErrorIgnoreLevel = kInfo + 1; 

  sourceFile->cd(path.Data());

  Mypath = path;
  Mypath.ReplaceAll("CaloQA_","");
  Mypath.ReplaceAll("default","trigMB");
  Mypath.ReplaceAll("trig","");

  TDirectory *current_sourcedir = gDirectory;

  TKey *key;
  TIter nextkey(current_sourcedir->GetListOfKeys());

  while ((key = (TKey*)nextkey())) 
  {

    Myobj = key->ReadObj();

    if (TString(Myobj->IsA()->GetName()).Contains("TCanvas")) 
      {
	title = Myobj->GetTitle();

	((TCanvas*)Myobj)->Print(Form("%s%s",f.Data(),".pdf"),Form("Title:%s/%s",Mypath.Data(),title));
      }
    else if ( Myobj->IsA()->InheritsFrom( "TDirectory" ) ) 
     {

      // Myobj is now the starting point of another iteration
      // obj still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      recurseOverKeys( (TDirectory*)Myobj, f );

    } // end of IF a TDriectory

  } // end of LOOP over keys

}
