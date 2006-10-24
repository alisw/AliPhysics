#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <TCanvas.h>

#include "AliSTARTCalc.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TFile.h"
#include "AliLog.h"
#include "TObjString.h"

#include "TAxis.h"
#include "TH2F.h"


ClassImp(AliSTARTCalc)

AliSTARTCalc::AliSTARTCalc() 
{
 //
 //	fGraphs.SetOwner(1);
}

AliSTARTCalc::AliSTARTCalc(const char* name)
{
	TString namst = "Calib_";
	namst += name;
	SetName(namst.Data());
	SetTitle(namst.Data());
//	fGraphs.SetOwner(1);
	Reset();
			  
}

//________________________________________________________________
AliSTARTCalc::AliSTARTCalc(const AliSTARTCalc& calibdata) : TNamed(calibdata)
 
{ 
// copy constructor
    SetName(calibdata.GetName());
    SetTitle(calibdata.GetName());
        

}

//________________________________________________________________
AliSTARTCalc &AliSTARTCalc::operator =(const AliSTARTCalc& calibdata)
{
// assignment operator
     SetName(calibdata.GetName());
     SetTitle(calibdata.GetName());

     return *this;
}

//________________________________________________________________
AliSTARTCalc::~AliSTARTCalc()
{
//
}

void AliSTARTCalc::Reset()
{
    memset(fTime,1,24*sizeof(Float_t));
   
}

void AliSTARTCalc::SetTime(Float_t* daqtime, Float_t* time_shift)
{ 
	for(Int_t i=0;i<24;i++){
		if (time_shift[i] != 0.)
		  fTime[i] = daqtime[i]-time_shift[i];
		else 
		  fTime[i] = daqtime[i];
		}
}


void AliSTARTCalc::Print(const Option_t*) const
{
	for(Int_t i=0;i<24;i++){
		printf("Total time %d = %.2f\n",i,fTime[i]);
	}
}



