#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <TCanvas.h>

#include "AliT0Calc.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TFile.h"
#include "AliLog.h"
#include "TObjString.h"

#include "TAxis.h"
#include "TH2F.h"


ClassImp(AliT0Calc)

AliT0Calc::AliT0Calc() 
{
 //
 //	fGraphs.SetOwner(1);
}

AliT0Calc::AliT0Calc(const char* name)
{
	TString namst = "Calib_";
	namst += name;
	SetName(namst.Data());
	SetTitle(namst.Data());
//	fGraphs.SetOwner(1);
	Reset();
			  
}

//________________________________________________________________
AliT0Calc::AliT0Calc(const AliT0Calc& calibdata) : TNamed(calibdata)
 
{ 
// copy constructor
    SetName(calibdata.GetName());
    SetTitle(calibdata.GetName());
        

}

//________________________________________________________________
AliT0Calc &AliT0Calc::operator =(const AliT0Calc& calibdata)
{
// assignment operator
     SetName(calibdata.GetName());
     SetTitle(calibdata.GetName());

     return *this;
}

//________________________________________________________________
AliT0Calc::~AliT0Calc()
{
//
}

void AliT0Calc::Reset()
{
    memset(fTime,1,24*sizeof(Float_t));
   
}

void AliT0Calc::SetTime(Float_t* daqtime, Float_t* time_shift)
{ 
	for(Int_t i=0;i<24;i++){
		if (time_shift[i] != 0.)
		  fTime[i] = daqtime[i]-time_shift[i];
		else 
		  fTime[i] = daqtime[i];
		}
}


void AliT0Calc::Print(const Option_t*) const
{
	for(Int_t i=0;i<24;i++){
		printf("Total time %d = %.2f\n",i,fTime[i]);
	}
}



