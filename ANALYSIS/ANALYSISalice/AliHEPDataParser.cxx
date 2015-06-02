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
#include "TMath.h"
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(AliHEPDataParser)

AliHEPDataParser::AliHEPDataParser() : TObject(), fHistStat(0),  fHistSyst(0),  fGraphStat(0),  fGraphSyst(0),  fHEPDataFileLines(0), fValueName(""), fXaxisName(""), fTitle(""), fReaction(""), fEnergy(""), fRapidityRange(""), fPrecision(5)
{
  // default ctor

}

AliHEPDataParser::AliHEPDataParser(TH1 * hStat, TH1 * hSyst): TObject(), fHistStat(0),  fHistSyst(0),  fGraphStat(0),  fGraphSyst(0),  fHEPDataFileLines(0), fValueName("y"), fXaxisName(""), fTitle(""), fReaction(""), fEnergy(""), fRapidityRange(""), fPrecision(5)
{
  //ctor
  fHistStat = hStat;
  fHistSyst = hSyst;
  fHEPDataFileLines = new TObjArray;

}
AliHEPDataParser::AliHEPDataParser(TGraph * grStat, TGraph * grSyst): TObject(), fHistStat(0),  fHistSyst(0),  fGraphStat(0),  fGraphSyst(0),  fHEPDataFileLines(0), fValueName(""), fXaxisName(""), fTitle(""), fReaction(""), fEnergy(""), fRapidityRange(""), fPrecision(5)
{
  // ctor
  fGraphStat = grStat;
  fGraphSyst = grSyst;
  fHEPDataFileLines = new TObjArray;
}

AliHEPDataParser::AliHEPDataParser(const char * hepfileName): TObject(), fHistStat(0),  fHistSyst(0),  fGraphStat(0),  fGraphSyst(0),  fHEPDataFileLines(0), fValueName("y"), fXaxisName(""), fTitle(""), fReaction(""), fEnergy(""), fRapidityRange(""), fPrecision(5)
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
    TObjArray * tokens = line.Tokenize(" \t");
    if(!tokens) continue;
    if(! ((TObjString*) tokens->At(0))->String().Atof()) {
      //The first column is not a number: those are the headers: skip!
      delete tokens;
      continue;
    }
    if(tokens->GetEntries() != 8) {
      // this should never happen!
      delete tokens;
      AliError(Form("Wrong number of columns %d! Assuming [x xlow xhigh y dystat+ dystat- dysyst+ dysyst-]", tokens->GetEntries()));
      cout << line.Data() << endl;
      return;      
      //continue;
    }
    // FIXME: Assumes the format
    // x xlow xhigh y dystat+ dystat- dysyst+ dysyst-
    TObjString * xbin   = (TObjString*) tokens->At(0);
    TObjString * xlow   = (TObjString*) tokens->At(1);
    TObjString * xhigh  = (TObjString*) tokens->At(2);
    TObjString * value  = (TObjString*) tokens->At(3);
    TObjString * statp  = (TObjString*) tokens->At(4);
    TObjString * statm  = (TObjString*) tokens->At(5);
    TObjString * systp  = (TObjString*) tokens->At(6);
    TObjString * systm  = (TObjString*) tokens->At(7);
    statm->String().ReplaceAll("+-","");
    statp->String().ReplaceAll("+-","");
    if(systp) systp->String().ReplaceAll("+-","");
    if(systm) systm->String().ReplaceAll("+-","");
    //    if (!binMin->String().Atof()) {delete tokens; continue;} // skip headers
    Float_t binCenter = xbin->String().Atof();
    Float_t binWidth  =  (xlow->String().Atof() - xhigh->String().Atof())/2;


    fGraphStat->SetPoint(ipoints, binCenter, value->String().Atof());
    fGraphSyst->SetPoint(ipoints, binCenter, value->String().Atof());
    ((TGraphAsymmErrors*)fGraphStat)->SetPointError(ipoints, 
						    binWidth,
						    binWidth,
						    TMath::Abs(statp->String().Atof()),
						    TMath::Abs(statm->String().Atof()));
    if(systp && systm) ((TGraphAsymmErrors*)fGraphSyst)->SetPointError(ipoints, 
                                                                       binWidth,
                                                                       binWidth,
                                                                       TMath::Abs(systp->String().Atof()),
                                                                       TMath::Abs(systm->String().Atof()));
    ipoints++;
    delete tokens;
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
  // Write headers if relevant
  if(fTitle.Length())         fHEPDataFileLines->Add(new TObjString(fTitle));
  if(fReaction.Length())      fHEPDataFileLines->Add(new TObjString(fReaction));
  if(fEnergy.Length())        fHEPDataFileLines->Add(new TObjString(fEnergy));
  if(fRapidityRange.Length()) fHEPDataFileLines->Add(new TObjString(fRapidityRange));
  if(!fValueName.Length() || !fXaxisName.Length()) AliFatal("At least x and y titles should be given!");
  fHEPDataFileLines->Add(new TObjString(Form("x: %s", fXaxisName.Data())));
  fHEPDataFileLines->Add(new TObjString(Form("y: %s", fValueName.Data())));


  if(trueUseGraphFalesUseHisto) {
    AliWarning("Graph saving not thoroughly tested!!");
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
	  fHEPDataFileLines->Add(new TObjString("x\txlow\txhigh\t+stat\t-stat\t+syst\t-syst"));
	}
	else {
	  fHEPDataFileLines->Add(new TObjString("x\txlow\txhigh\t+stat\t-stat"));
	}
      }
      // Skip empty bins
      if(!fGraphStat->GetY()[ipoint]) continue;
      TObjString * line = new TObjString;    
      if(fGraphSyst) {
	if (asym)
	  line->String().Form("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			      //	line->String().Form("%10s %10s %10s %10s %10s %10s %10s %10s", 
			      GetFixWidthCol(fGraphStat->GetX()[ipoint],                                                                   10).Data(), 			    
			      GetFixWidthCol(RoundToSignificantFigures(fGraphStat->GetX()[ipoint]-((TGraphAsymmErrors*)fGraphStat)->GetEXlow()[ipoint],  fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(fGraphStat->GetX()[ipoint]+((TGraphAsymmErrors*)fGraphStat)->GetEXhigh()[ipoint], fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(fGraphStat->GetY()[ipoint],                            fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphAsymmErrors*)fGraphStat)->GetEYhigh()[ipoint], fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphAsymmErrors*)fGraphStat)->GetEYlow()[ipoint] , fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphAsymmErrors*)fGraphSyst)->GetEYhigh()[ipoint], fPrecision), 10).Data(),
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphAsymmErrors*)fGraphSyst)->GetEYlow()[ipoint] , fPrecision), 10).Data());
	else 
	  line->String().Form("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			      //	line->String().Form("%10s %10s %10s %10s %10s %10s %10s %10s", 
			      GetFixWidthCol(fGraphStat->GetX()[ipoint],                                                                   10).Data(), 			    
			      GetFixWidthCol(fGraphStat->GetX()[ipoint]-fGraphStat->GetEX()[ipoint],                                    10).Data(), 
			      GetFixWidthCol(fGraphStat->GetX()[ipoint]+fGraphStat->GetEX()[ipoint],                                10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(fGraphStat->GetY()[ipoint],                            fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphErrors*)fGraphStat)->GetEY()[ipoint], fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphErrors*)fGraphStat)->GetEY()[ipoint], fPrecision), 10).Data(),
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphErrors*)fGraphSyst)->GetEY()[ipoint], fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphErrors*)fGraphSyst)->GetEY()[ipoint], fPrecision), 10).Data());
	  // line->String().Form("%f %f +-%f +-%f", 
	  // 		      fGraphStat->GetX()[ipoint], RoundToSignificantFigures(fGraphStat->GetY()[ipoint],fPrecision),
	  // 		      RoundToSignificantFigures(((TGraphErrors*)fGraphStat)->GetEY()[ipoint],fPrecision), 
	  // 		      RoundToSignificantFigures(((TGraphErrors*)fGraphSyst)->GetEY()[ipoint],fPrecision));

	fHEPDataFileLines->Add(line);
      } else {
	if (asym)
	  line->String().Form("%s\t%s\t%s\t%s\t%s\t%s",
			      //	line->String().Form("%10s %10s %10s %10s %10s %10s %10s %10s", 
			      GetFixWidthCol(fGraphStat->GetX()[ipoint],                                                                   10).Data(), 			    
			      GetFixWidthCol(fGraphStat->GetX()[ipoint]-fGraphStat->GetEXlow()[ipoint],                                    10).Data(), 
			      GetFixWidthCol(fGraphStat->GetX()[ipoint]+fGraphStat->GetEXhigh()[ipoint],                                   10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(fGraphStat->GetY()[ipoint],                            fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphAsymmErrors*)fGraphStat)->GetEYhigh()[ipoint], fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphAsymmErrors*)fGraphStat)->GetEYlow()[ipoint] , fPrecision), 10).Data()); 
	else 
	  line->String().Form("%s\t%s\t%s\t%s\t%s\t%s",
			      //	line->String().Form("%10s %10s %10s %10s %10s %10s %10s %10s", 
			      GetFixWidthCol(fGraphStat->GetX()[ipoint],                                                                   10).Data(), 			    
			      GetFixWidthCol(fGraphStat->GetX()[ipoint]-fGraphStat->GetEX()[ipoint],                                    10).Data(), 
			      GetFixWidthCol(fGraphStat->GetX()[ipoint]+fGraphStat->GetEX()[ipoint],                                10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(fGraphStat->GetY()[ipoint],                            fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphErrors*)fGraphStat)->GetEY()[ipoint], fPrecision), 10).Data(), 
			      GetFixWidthCol(RoundToSignificantFigures(((TGraphErrors*)fGraphStat)->GetEY()[ipoint], fPrecision), 10).Data());

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
	  fHEPDataFileLines->Add(new TObjString("x\t\txlow\t\txhigh\t\ty\t\tdystat+\t\tdystat-\t\tdysyst+\t\tdysyst-")); 
	else 
	  fHEPDataFileLines->Add(new TObjString("x\t\txlow\t\txhigh\t\ty\t\tdystat+\t\tdystat-"));       	
      }
      // Skip empty bins
      if(!fHistStat->GetBinContent(ibin)) continue;
      TObjString * line = new TObjString;      
      if(fHistSyst) {
	
	line->String().Form("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			    //	line->String().Form("%10s %10s %10s %10s %10s %10s %10s %10s", 
			    GetFixWidthCol(fHistStat->GetBinLowEdge(ibin)+fHistStat->GetBinWidth(ibin)/2,         12).Data(), 			    
			    GetFixWidthCol(fHistStat->GetBinLowEdge(ibin),                                        12).Data(), 
			    GetFixWidthCol(fHistStat->GetBinLowEdge(ibin)+fHistStat->GetBinWidth(ibin),           12).Data(), 
			    GetFixWidthCol(RoundToSignificantFigures(fHistStat->GetBinContent(ibin),fPrecision),  12).Data(), 
			    GetFixWidthCol(RoundToSignificantFigures(fHistStat->GetBinError(ibin),  fPrecision),  12).Data(), 
			    GetFixWidthCol(RoundToSignificantFigures(fHistStat->GetBinError(ibin),  fPrecision),  12).Data(), 
			    GetFixWidthCol(RoundToSignificantFigures(fHistSyst->GetBinError(ibin),  fPrecision),  12).Data(),
			    GetFixWidthCol(RoundToSignificantFigures(fHistSyst->GetBinError(ibin),  fPrecision),  12).Data());

	fHEPDataFileLines->Add(line);
      } else {
	line->String().Form("%s\t%s\t%s\t%s\t%s\t%s", 
			    GetFixWidthCol(fHistStat->GetBinLowEdge(ibin)+fHistStat->GetBinWidth(ibin)/2,         10).Data(), 			    
			    GetFixWidthCol(fHistStat->GetBinLowEdge(ibin),                                        10).Data(), 
			    GetFixWidthCol(fHistStat->GetBinLowEdge(ibin)+fHistStat->GetBinWidth(ibin),           10).Data(), 
			    GetFixWidthCol(RoundToSignificantFigures(fHistStat->GetBinContent(ibin),fPrecision),  10).Data(), 
			    GetFixWidthCol(RoundToSignificantFigures(fHistStat->GetBinError(ibin),  fPrecision),  10).Data(), 
			    GetFixWidthCol(RoundToSignificantFigures(fHistStat->GetBinError(ibin),  fPrecision),  10).Data()); 
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


Double_t AliHEPDataParser::RoundToSignificantFigures(double num, int n) {
  // Rounds num to n significant digits.
  // Recipe from :http://stackoverflow.com/questions/202302/rounding-to-an-arbitrary-number-of-significant-digits
  // Basically the log is used to determine the number of leading 0s, than convert to an integer by multipliing by the expo, 
  // round the integer and shift back.
  if(num == 0) {
    return 0;
  }

  Double_t d = TMath::Ceil(TMath::Log10(num < 0 ? -num: num));
  Int_t power = n - (int) d;

  Double_t magnitude = TMath::Power(10, power);
  Long_t shifted = TMath::Nint(num*magnitude);
  return shifted/magnitude;

}

TString AliHEPDataParser::GetFixWidthCol(Double_t number, Int_t width) {

  // Formats a column to fixed width
  TString col;
  char format[100];
  snprintf(format,100,"%%%d#g", fPrecision);
  col.Form(format, number);
  if(col.Length()>width) AliError("larger than width, cannot align!");

  if(col.Contains("e"))
    while (col.Length() < width) col.Append(" ");
  else
    while (col.Length() < width) col.Append("0");
  
  return col;
}
