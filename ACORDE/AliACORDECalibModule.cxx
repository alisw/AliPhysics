/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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




///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ACORDE calibration                                              //
// Pedro Gonzalez Zamora    pedro.gonzalez@fcfm.buap.mx                      //
// Irais Bautista Guzman    irais@fcfm.buap.mx                               //
// Arturo Fernandez Tellez afernan@cern.ch                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliACORDECalibModule.h"
#include "AliACORDEDataModule.h"
#include "AliLog.h"
#include <TCanvas.h>
                                                                                                                                                             

ClassImp(AliACORDECalibModule)


//________________________________________________________________

AliACORDECalibModule::AliACORDECalibModule():
TObject(),
fHRate(0x0)
{


       // Default constructor
       for(int i=0;i<kNmodules;i++)
       fModule[i]=0x0;


}
//________________________________________________________________
AliACORDECalibModule::AliACORDECalibModule(const AliACORDECalibModule &calibdata):
TObject(),
fHRate(calibdata.fHRate)
{
 //copy constructor
 
  for (Int_t i=0;i<kNmodules;i++)
 {
   fModule[i]=calibdata.fModule[i];
 
 }
 

}
//_______________________________________________________________
AliACORDECalibModule& AliACORDECalibModule:: operator=(const AliACORDECalibModule & calibdata)
{
//assignment operator
 	for(Int_t i=0;i<kNmodules;i++)
	{
 	this->fModule[i]=calibdata.GetModule(i);

	}
    
  this->fHRate=calibdata.fHRate;
  
  return *this;
}



//________________________________________________________________
AliACORDECalibModule::~AliACORDECalibModule()
{

  	// Destructor

}


//________________________________________________________________
void AliACORDECalibModule::SetModule(Int_t module,Float_t value,Bool_t status,const char* name)
{ 
 
     	fModule[module] = new AliACORDEDataModule(value,status,name); 
     	//Set values for each module
}
//________________________________________________________________
void AliACORDECalibModule::SetModuleRate(Int_t module, Float_t value)
{
  	// Set RATE  value
  	
           fModule[module]->SetRate(value);
}
//________________________________________________________________

Float_t AliACORDECalibModule::GetModuleRate(Int_t module)
{
           
          return fModule[module]->GetRate();
}
//________________________________________________________________
         
void AliACORDECalibModule::Create_Histo()
{
           //Create histogram of modules rates actual values   
         
         
           fHRate = new TH1F("fHRate","RATES PER MODULE",60,0,60); 
           

           for(int i=0;i<kNmodules;i++)
           fHRate->Fill(i,fModule[i]->GetRate());

}
//_____________________________________________________________________       
void AliACORDECalibModule::Draw(const Option_t* /*option*/) //Draw the Histogram RATE_MODULE
{
        TCanvas *c1 = new TCanvas
          ("RATES","RATES PER MODULE",200,10,700,500);
           c1->cd();
           fHRate->Draw();
}

//______________________________________________________________________
void AliACORDECalibModule::Print_Module() // Print the status and rates for each module
{


       	for(int module=0;module<kNmodules;module++)
	AliInfo(Form("%s con rate %f  Status %d \n",fModule[module]->GetName(),fModule[module]->GetRate(),fModule[module]->GetStatus()));

}
