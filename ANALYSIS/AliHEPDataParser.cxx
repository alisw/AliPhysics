//-------------------------------------------------------------------------
//           Implementation of Class AliHEPDataParser
//
//  This class is used to save the content of hisograms/graphs in the
//  HEP data format and viceversa. The HEP data format is not strictly
//  defined and there are many variants, the class only support a few
//  of them. More will be added as needed.  The input can be a set of
//  2 TH1, TGraphAsymmErrors or TGraphErrors (one for the stat and one
//  for the syst error). If the second one is a null pointer, only the
//  stat error is printed. The class can also import hepdata ascii
//  file (very preliminary)
// 
//  Author: Michele Floris, CERN
//-------------------------------------------------------------------------


#include "AliHEPDataParser.h"
#include "AliLog.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(AliHEPDataParser)

AliHEPDataParser::AliHEPDataParser() : TObject(), fHistStat(0),  fHistSyst(0),  fGraphStat(0),  fGraphSyst(0),  fHEPDataFileLines(0), fValueName("")
{
  // default ctor

}

AliHEPDataParser::AliHEPDataParser(TH1 * hStat, TH1 * hSyst): TObject(), fHistStat(0),  fHistSyst(0),  fGraphStat(0),  fGraphSyst(0),  fHEPDataFileLines(0), fValueName("y")
{
  //ctor
  fHistStat = hStat;
  fHistSyst = hSyst;
  fHEPDataFileLines = new TObjArray;

}
AliHEPDataParser::AliHEPDataParser(TGraph * grStat, TGraph * grSyst): TObject(), fHistStat(0),  fHistSyst(0),  fGraphStat(0),  fGraphSyst(0),  fHEPDataFileLines(0), fValueName("")
{
  // ctor
  fGraphStat = grStat;
  fGraphSyst = grSyst;
  fHEPDataFileLines = new TObjArray;
}

AliHEPDataParser::AliHEPDataParser(const char * hepfileName): TObject(), fHistStat(0),  fHistSyst(0),  fGraphStat(0),  fGraphSyst(0),  fHEPDataFileLines(0), fValueName("y")
{
  //ctor
  // Put results in graphs
  fGraphSyst = new TGraphAsymmErrors();
  fGraphStat = new TGraphAsymmErrors();
  ifstream infile;
  infile.open(hepfileName);
  TString line;
  Int_t ipoints = 0;
  while (line.ReadLine(infile)) {
    TObjArray * tokens = line.Tokenize(" ");
    if(tokens->GetEntries() < 1) {
      delete tokens;
      AliError("not enough columns");
      return;      
    }
    // TODO: Assumes the format
    // binmin binmax y +-stat +-syst. Try to make it smarter...
    TObjString * binMin = (TObjString*) tokens->At(0);
    TObjString * binMax = (TObjString*) tokens->At(1);
    TObjString * value  = (TObjString*) tokens->At(2);
    TObjString * stat   = (TObjString*) tokens->At(3);
    TObjString * syst   = (TObjString*) tokens->At(4);
    stat->String().ReplaceAll("+-","");
    if(syst) syst->String().ReplaceAll("+-","");
    if (!binMin->String().Atof()) continue; // skip headers
    Float_t binCenter = (binMax->String().Atof() + binMin->String().Atof())/2;
    Float_t binWidth  =  (binMax->String().Atof() - binMin->String().Atof())/2;
    cout << line.Data() << endl;//<< " " << binMin->String().Atof() <<" " <<  binCenter << " " << binWidth << endl;


    fGraphStat->SetPoint(ipoints, binCenter, value->String().Atof());
    fGraphSyst->SetPoint(ipoints, binCenter, value->String().Atof());
    ((TGraphAsymmErrors*)fGraphStat)->SetPointError(ipoints, 
						    binWidth,
						    binWidth,
						    stat->String().Atof(),
						    stat->String().Atof());
    if(syst) ((TGraphAsymmErrors*)fGraphSyst)->SetPointError(ipoints, 
							     binWidth,
							     binWidth,
							     syst->String().Atof(),
							     syst->String().Atof());
    ipoints++;
  }
  infile.close();
    

}

AliHEPDataParser::~AliHEPDataParser(){
  // dtor
  if(fHistStat) delete fHistStat;
  if(fHistSyst) delete fHistSyst;
  if(fGraphStat) delete fGraphStat;
  if(fGraphSyst) delete fGraphSyst;
  if(fHEPDataFileLines) delete fHEPDataFileLines;
}
  
void AliHEPDataParser::SaveHEPDataFile(const char * hepfileName, Bool_t trueUseGraphFalesUseHisto) {
  // Fills fHEPDataFileLines and saves its content to a file
  if(!fHEPDataFileLines)   fHEPDataFileLines = new TObjArray;
  if(trueUseGraphFalesUseHisto) {
    if(!fGraphStat) {
      AliError("Graph not set");
      return;
    }
    Bool_t asym = kFALSE; // check if this has asymmetric errors
    if       (!strcmp(fGraphStat->ClassName(), "TGraphErrors")) asym = kFALSE;
    else if  (!strcmp(fGraphStat->ClassName(), "TGraphAsymmErrors")) asym = kTRUE;
    else     {AliError("Unsupported graph type"); return;}
    Int_t npoint = fGraphStat->GetN();
    if(asym) AliInfo("Assymmetric errors");
    for(Int_t ipoint = 0; ipoint < npoint; ipoint++){            
      if(ipoint == 0) {
	if(fGraphSyst) {
	  if (asym)
	    fHEPDataFileLines->Add(new TObjString(Form("BinCenter %s +stat -stat +syst -syst", fValueName.Data()))); 
	  else 
	    fHEPDataFileLines->Add(new TObjString(Form("BinCenter %s +-stat +-syst", fValueName.Data()))); 
	}
	else {
	  if(asym)
	    fHEPDataFileLines->Add(new TObjString(Form("BinCenter %s +stat -stat", fValueName.Data()))); 
	  else 
	    fHEPDataFileLines->Add(new TObjString(Form("BinCenter %s +-stat", fValueName.Data()))); 
	}
      }
      // Skip empty bins
      if(!fGraphStat->GetY()[ipoint]) continue;
      TObjString * line = new TObjString;      
      if(fGraphSyst) {
	if (asym)
	  line->String().Form("%f %f +%f -%f +%f -%f", 
			      fGraphStat->GetX()[ipoint], fGraphStat->GetY()[ipoint],
			      ((TGraphAsymmErrors*)fGraphStat)->GetEYhigh()[ipoint], 
			      ((TGraphAsymmErrors*)fGraphStat)->GetEYlow()[ipoint], 
			      ((TGraphAsymmErrors*)fGraphSyst)->GetEYhigh()[ipoint], 
			      ((TGraphAsymmErrors*)fGraphSyst)->GetEYlow()[ipoint]);
	else 
	  line->String().Form("%f %f +-%f +-%f", 
			      fGraphStat->GetX()[ipoint], fGraphStat->GetY()[ipoint],
			      ((TGraphErrors*)fGraphStat)->GetEY()[ipoint], 
			      ((TGraphErrors*)fGraphSyst)->GetEY()[ipoint]);

	fHEPDataFileLines->Add(line);
      } else {
	if (asym)
	  line->String().Form("%f %f +%f -%f", 
			      fGraphStat->GetX()[ipoint], fGraphStat->GetY()[ipoint],
			      ((TGraphAsymmErrors*)fGraphStat)->GetEYhigh()[ipoint], ((TGraphAsymmErrors*)fGraphStat)->GetEYlow()[ipoint]);
	else { 
	  line->String().Form("%f %f +-%f", 
			      fGraphStat->GetX()[ipoint], fGraphStat->GetY()[ipoint],
			      ((TGraphErrors*)fGraphStat)->GetEY()[ipoint]);
	}

	fHEPDataFileLines->Add(line);
      }
    }    
  }
  else {
    if(!fHistStat) {
      AliError("Hist not set");
      return;
    }
    Int_t nbin = fHistStat->GetNbinsX();
    
    for(Int_t ibin = 1; ibin <= nbin; ibin++){
      if(ibin == 1) {
	if(fHistSyst)
	  fHEPDataFileLines->Add(new TObjString(Form("BinLow BinHigh %s +-stat +-syst", fValueName.Data()))); 
	else 
	  fHEPDataFileLines->Add(new TObjString(Form("BinLow BinHigh %s +-stat", fValueName.Data())));       	
      }
      // Skip empty bins
      if(!fHistStat->GetBinContent(ibin)) continue;
      TObjString * line = new TObjString;      
      if(fHistSyst) {
	line->String().Form("%f %f %f +-%f +-%f", 
			    fHistStat->GetBinLowEdge(ibin), 
			    fHistStat->GetBinLowEdge(ibin)+fHistStat->GetBinWidth(ibin), 
			    fHistStat->GetBinContent(ibin), fHistStat->GetBinError(ibin), fHistSyst->GetBinError(ibin));
	fHEPDataFileLines->Add(line);
      } else {
	line->String().Form("%f %f %f +-%f", 
			    fHistStat->GetBinLowEdge(ibin), 
			    fHistStat->GetBinLowEdge(ibin)+fHistStat->GetBinWidth(ibin), 
			    fHistStat->GetBinContent(ibin), fHistStat->GetBinError(ibin));
	fHEPDataFileLines->Add(line);
      }
      //      delete line;      
    }           
  }

  TIterator * lineIter = fHEPDataFileLines->MakeIterator();
  TObjString * obj = 0;
  ofstream outfile;
  outfile.open (hepfileName);
  cout << "Saving HEP File " << hepfileName << endl;
  
  while ((obj = (TObjString*) lineIter->Next())) {    
    cout << obj->String().Data() << endl;    
    outfile << obj->String().Data() << endl;
  }
  outfile.close();
}


