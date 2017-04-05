#include <iostream>
#include <stdio.h>
#include "TString.h"
using namespace std;
void runeffmaker(const char* options = "")
{
  TString basedir=TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  gROOT->LoadMacro(Form("%seffmaker.cxx+g",basedir.Data()));
//   MakeTestHists();
  if(TString(options).Contains("pp")){
    MakeEffHistspp();
  }
  if(TString(options).Contains("PbPb")){
    MakeEffHistsPbPb();
  }
}