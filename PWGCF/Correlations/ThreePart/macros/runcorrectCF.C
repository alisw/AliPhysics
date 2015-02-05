#include <iostream>
#include "TString.h"
#include "TPRegexp.h"
using namespace std;
void runcorrectCF(const char* options = "")
{
 
  TString option = TString("");
  bool clean = false;
  bool Merge = false;
  bool Draw = false;
  bool collectMVbins = false;
  bool correctforME = false;
  bool makeyield = false;
  bool drawgen = false;
  TString delimiter(" ");
  TStringToken token(options, delimiter);
  bool madeachoice = false;
  while (token.NextToken()&&!madeachoice) {
    const char* key=0;
    TString argument=token;
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
    key="Draw";
    if (argument.CompareTo(key)==0) {
      Draw = true;
      collectMVbins =false;
      correctforME = false;
      makeyield    = false;
      madeachoice = true;
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
  }
  TString basedir=TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  if(collectMVbins||correctforME||makeyield) gROOT->LoadMacro("correctCF.cxx+g");
  if(Merge||Draw){	gROOT->LoadMacro(Form("%srunan.C",basedir.Data()));runan("compile");}

  if(Merge)		gROOT->LoadMacro("MergeSet.C+g");
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

 
  if(clean){
    cout << "Cleaning builds, removing local libraries. If you want to clean builds, please answer yes. ";
    gROOT->ProcessLine(".rm *.so *.d");
    cout << "Removing results.root and corrected.root. If you want to remove them, please answer yes. ";
    gROOT->ProcessLine(".rm results.root corrected.root");
    cout << "Removing AnalysisResults.root . If you want to remove it, please answer yes. ";
    gROOT->ProcessLine(".rm AnalysisResults.root ");    
  }
  if(Merge){MergeSets();return;}
  if(Draw&&!drawgen){	DrawThreeParticleCorrelations("makefile"); DrawThreeParticleCorrelations("makefilemixed");DrawThreeParticleCorrelations("makefiletriggermixed");DrawThreeParticleCorrelations("makefilemixedall");}
  if(Draw&&drawgen){ 	DrawThreeParticleCorrelations("makefilegen");DrawThreeParticleCorrelations("makefilegen makefilemixed");DrawThreeParticleCorrelations("makefilegen makefiletriggermixed");}
  if(collectMVbins)	CollectMVbins();
  if(correctforME)	correct(option);
  if(makeyield)		yield();
}
