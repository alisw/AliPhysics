#include <iostream>
#include <stdio.h>
#include "TString.h"
#include "TPRegexp.h"
using namespace std;
void runcorrectCF(const char* options = "")
{
 
  TString option = TString("");
  bool clean = false;
  bool Merge = false;
  bool Draw = false;
  int drawint =0;
  bool collectMVbins = false;
  bool collectMVbinsfirst = false;
  bool correctforME = false;
  bool correctforMEscan = false;
  bool makeyield = false;
  bool drawgen = false;
  bool result = false;
  bool checkmixed = false;
  bool periods=false;
  bool centralities=false;
  bool CountRawNs = false;
  TString delimiter(" ");
  TStringToken token(options, delimiter);
  bool madeachoice = false;
  cout << options<<endl;
  while (token.NextToken()&&!madeachoice) {
    const char* key=0;
    TString argument=token;
    cout << argument.Data()<<endl;    
    key="pp";
    if (argument.CompareTo(key)==0) {
      option.Clear();
      option.Append("pp");
      cout << "Using 'pp' as option."<<endl;
      continue;
    }
    key="PP";
    if (argument.CompareTo(key)==0) {
      option.Clear();
      option.Append("pp");
      cout << "Using 'pp' as option."<<endl;
      continue;
    }
    key="PbPb";
    if (argument.CompareTo(key)==0) {
      option.Clear();
      option.Append("PbPb");
      cout << "Using 'PbPb' as option."<<endl;
      continue;
    }
    key="pbpb";
    if (argument.CompareTo(key)==0) {
      option.Clear();
      option.Append("PbPb");
      cout << "Using 'PbPb' as option."<<endl;      
      continue;
    }
    key="48";
    if (argument.CompareTo(key)==0) {
      option.Append(" 48");
      continue;
    }
    key="816";
    if (argument.CompareTo(key)==0) {
      option.Append(" 816");
      continue;
    }
    key="gen";
    if (argument.CompareTo(key)==0) {
      drawgen = true;
      cout << "Setting to generated data."<<endl;      
      continue;
    }
    key="All";
    if (argument.CompareTo(key)==0) {
      Draw = true;
      collectMVbins =true;
      correctforME = true;
      makeyield    = true;
      madeachoice = true;
      cout << "Drawing, then Collecting multiplicity bins, then correcting for META and METrigger then calculating the yield." <<endl;
      continue;
    }    
    key="Draw1";
    if (argument.CompareTo(key)==0) {
      Draw = true;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      drawint = 1;
      cout << "Drawing the histograms and creating results.root (must use aliroot to use this)."<<endl;
      continue;
    }    
    key="Draw2";
    if (argument.CompareTo(key)==0) {
      Draw = true;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      drawint = 2;
      cout << "Drawing the histograms and creating results.root (must use aliroot to use this)."<<endl;
      continue;
    }  
    key="Draw3";
    if (argument.CompareTo(key)==0) {
      Draw = true;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      drawint = 3;
      cout << "Drawing the histograms and creating results.root (must use aliroot to use this)."<<endl;
      continue;
    }  
    key="Collect";
    if (argument.CompareTo(key)==0) {
      Draw = false;
      collectMVbins =true;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      cout << "Collecting multiplicity bins in results.root (must have drawn first)."<<endl;
      continue;
    }    
    key="CollectFirst";
    if (argument.CompareTo(key)==0) {
      Draw = false;
      collectMVbinsfirst =true;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      cout << "Collecting multiplicity bins in results.root (must have drawn first)."<<endl;
      continue;
    }    
    key="Correct";
    if (argument.CompareTo(key)==0) {
      Draw = false;
      collectMVbins =false;
      correctforME = true;
      makeyield    = false;
      madeachoice = true;
      cout << "Correcting for META and METrigger in results.root (must have drawn first)."<<endl;
      continue;
    }        
    key="Scan";
    if (argument.CompareTo(key)==0) {
      Draw = false;
      collectMVbins =false;
      correctforMEscan = true;
      makeyield    = false;
      madeachoice = true;
      cout << "Correcting for META and METrigger in results.root (must have drawn first)."<<endl;
      continue;
    }        
    key="Yield";
    if (argument.CompareTo(key)==0) {
      Draw = false;
      collectMVbins =false;
      correctforME = false;
      makeyield    = true;
      madeachoice = true;
      cout << "Calculating Yield in results.root (must have drawn & corrected first)."<<endl;
      continue;
    }        
    key="DrawCollect";
    if (argument.CompareTo(key)==0) {
      Draw = true;
      collectMVbins =true;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      cout << "Drawing the histograms, creating results.root and collecting multiplicity bins."<<endl;
      continue;
    }        
    key="CorrectYield";
    if (argument.CompareTo(key)==0) {
      Draw = true;
      collectMVbins =true;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      cout << "Correcting for META and METrigger in results.root and calculating the Yield (must have drawn first)."<<endl;
      continue;
    }        
    key="Merge";
    if (argument.CompareTo(key)==0) {
      Merge = true;
      Draw = false;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      cout << "Merging (only works in aliroot, not in root standalone)."<<endl;
      continue;
    }          
    key="Result";
    if (argument.CompareTo(key)==0) {
      result = true;
      Draw = false;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      cout << "Draws some canvases to compare results."<<endl;
      continue;
    }          
    key="clean";
    if (argument.CompareTo(key)==0) {
      clean = true;
      Draw = false;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      cout << "Cleaning the directory for compiled objects and result files. You will have to rerun all commands."<<endl;
      continue;
    }
    key="CheckMixed";
    if (argument.CompareTo(key)==0) {
      checkmixed = true;
      Draw = false;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      continue;
    }
    key="Periods";
    if (argument.CompareTo(key)==0) {
      checkmixed = false;
      periods = true;
      Draw = false;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      option.Append("pp");
      
      continue;
    }    
    key="Centralities";
    if (argument.CompareTo(key)==0) {
      checkmixed = false;
      centralities = true;
      Draw = false;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      option.Append("pp");
      
      cout << "Cleaning the directory for compiled objects and result files. You will have to rerun all commands."<<endl;
      continue;
    }    
    key="Count";
    if (argument.CompareTo(key)==0) {
      checkmixed = false;
      CountRawNs = true;
      Draw = false;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
      option.Append("pp");
      
      cout << "Cleaning the directory for compiled objects and result files. You will have to rerun all commands."<<endl;
      continue;
    }    
    
  }
  TString basedir=TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  if(collectMVbins||collectMVbinsfirst||correctforME||makeyield||result||correctforMEscan||checkmixed||periods||centralities||CountRawNs) gROOT->LoadMacro(Form("%scorrectCF.cxx+g",basedir.Data()));
  if(Merge||Draw){	gROOT->LoadMacro(Form("%srunan.C",basedir.Data()));runan("compile");}
  if(Merge)		gROOT->LoadMacro(Form("%sMergeSet.C+g",basedir.Data()));
  if(Draw) 		gROOT->LoadMacro(Form("%sDrawThreeParticleCorrelations.C",basedir.Data()));
  if(!madeachoice) cout << "Please choose an action to perform. This macro can:                  "<<endl
			<< "    option: 	task performed:		                         "<<endl
			<< "    clean     -     cleans the directory                             "<<endl
			<< "    Merge     -     Merges the Runs. Please run with Aliroot.        "<<endl
			<< "    Draw      -     Extracts the histograms and creates results.root."<<endl
			<< "                    Please run with Aliroot.                         "<<endl
			<< "    Collect   -     Collects the Multiplicity bins.                  "<<endl
			<< "                    Draw must be performed first                     "<<endl
			<< "    Correct   -     Corrects the Histograms using the different types"<<endl
			<< "			of Mixed event.                                  "<<endl
			<< "                    Draw must be performed first.                    "<<endl
			<< "    Yield     -     Calculates the Yield in each DPHI bin.           "<<endl
			<< "                    Draw and Correct must be performed first.        "<<endl
			<< "    All       -     Performs Draw, Collect Correct and Yield.        "<<endl
			<< "                    Must be run with Aliroot.                        "<<endl<<endl
			<< "To choose the collision type, please use pp or PbPb.                 "<<endl<<endl;
  if(option.CompareTo("")==0){ cout<< "Defaulting to option 'pp'." <<endl; option.Append("pp");}
  if(Draw && option.CompareTo("pp")==0){option.Append(" 48 816");}
  if(clean){
    cout << "Cleaning builds, removing local libraries. If you want to clean builds, please answer yes. ";
    gROOT->ProcessLine(".rm *.so *.d");
    cout << "Removing results.root and corrected.root. If you want to remove them, please answer yes. ";
    gROOT->ProcessLine(".rm results.root corrected.root");
    cout << "Removing AnalysisResults.root . If you want to remove it, please answer yes. ";
    gROOT->ProcessLine(".rm AnalysisResults.root ");    
  }
  if(Merge){MergeSets();return;}
  if(Draw&&!drawgen&&drawint==1){DrawThreeParticleCorrelations(Form("%s makefile",option.Data())); }//DrawThreeParticleCorrelations(Form("%s makefilemixedall",option.Data()));}
  if(Draw&&!drawgen&&drawint==2){DrawThreeParticleCorrelations(Form("%s makefilemixed",option.Data()));}//DrawThreeParticleCorrelations(Form("%s makefilemixedall",option.Data()));}
  if(Draw&&!drawgen&&drawint==3){DrawThreeParticleCorrelations(Form("%s makefiletriggermixed",option.Data()));}

  if(Draw&&drawgen){ 	DrawThreeParticleCorrelations("makefilegen");DrawThreeParticleCorrelations("makefilegen makefilemixed");DrawThreeParticleCorrelations("makefilegen makefiletriggermixed");}
  if(collectMVbins)	CollectMVbins(true);
  if(collectMVbinsfirst)	CollectMVbinsFirst(option);

  if(correctforME)	correct(option);
  if(correctforMEscan)	ScanCorrections(option);
  if(makeyield)		yield(option);
  if(result)		draw(option);
  if(checkmixed)        Checkmixed(option);
  if(periods)		Periods();
  if(centralities)	Centralities();
  if(CountRawNs)	CountRawNumbers();
  
}
