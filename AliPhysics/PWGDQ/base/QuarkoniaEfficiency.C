#if !defined(__CINT__) || defined(__MAKECINT__)
//C++ includes
#include <Riostream.h>
#include <stdio.h>
#include <string.h>

//PWG3
#include "AliQuarkoniaEfficiency.h"
#endif


void QuarkoniaEfficiency(double Rapidity, double Pt){
  
  //Read Acceptance
  Double_t eff,error;
  AliQuarkoniaEfficiency * JPsiEff = new AliQuarkoniaEfficiency();
  JPsiEff->Init();
  JPsiEff->GetEfficiency(Rapidity,Pt,eff,error);
  printf(" Efficiency calculations give us: eff = %2.2e  error %2.2e \n",eff,error);

}
