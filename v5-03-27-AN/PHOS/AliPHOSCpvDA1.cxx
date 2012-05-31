/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
// Class AliPHOSCpvDA1 accumulates histograms with amplitudes per CPV channel.
// It is intended to run at DAQ or HLT computers.
// Author: Boris Polishchuk, 25 January 2008.
///////////////////////////////////////////////////////////////////////////////

#include "AliPHOSCpvDA1.h"
#include "TString.h"

ClassImp(AliPHOSCpvDA1)

//----------------------------------------------------------------
AliPHOSCpvDA1::AliPHOSCpvDA1(int module) : TNamed(),
 fHistoFile(0),fMod(module)

{
  // Create DA1 ("Calibration DA") object.
  // module is the CPV module number (0..4).
  // Checks existence of histograms which might have been left
  // from the previous runs to continue their filling.
  // Histogram names: module_iX_iZ.
  // Root file name: CPV_ModuleX_Calib.root, where X - module number.
  
  char name[128];
  TString sname="CPV_Module%d_Calib";
  snprintf(name,sname.Length(),sname.Data(),fMod);
  SetName(name);

  char title[128];
  TString stitle="Calibration Detector Algorithm for CPV module %d";
  snprintf(title,stitle.Length(),stitle.Data(),fMod);
  SetTitle(title);

  char rootname[128];
  TString srootname="%s.root";
  snprintf(rootname,srootname.Length(),srootname.Data(),GetName());

  fHistoFile =  new TFile(rootname,"update");

  char hname[128];
  TH1F* hist1=0;
  TString shname="%d_%d_%d";

  for(Int_t iX=0; iX<128; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {
      snprintf(hname,shname.Length(),shname.Data(),fMod,iX,iZ);
      hist1 = (TH1F*)fHistoFile->Get(hname);
      if(hist1) fCharge[iX][iZ] = hist1;
      else
	fCharge[iX][iZ] = 0;
    }
  }
 
}

//-------------------------------------------------------------------
AliPHOSCpvDA1::AliPHOSCpvDA1(const AliPHOSCpvDA1& da) : TNamed(da),
  fHistoFile(0),fMod(da.fMod)
{
  // Copy constructor.

  fHistoFile = new TFile(da.GetName(),"update");

  char hname[128];
  TH1F* hist1=0;

  for(Int_t iX=0; iX<128; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {

      snprintf(hname,128,"%d_%d_%d",fMod,iX,iZ);
      hist1 = (TH1F*)da.fHistoFile->Get(hname);
      if(hist1) fCharge[iX][iZ] = new TH1F(*hist1);
      else
	fCharge[iX][iZ] = 0;
    }
  }
  
}

//-------------------------------------------------------------------
AliPHOSCpvDA1& AliPHOSCpvDA1::operator= (const AliPHOSCpvDA1& da)
{
  //Assignment operator.

  if(this != &da) {

    TString oldname(fHistoFile->GetName());
    TString newname(da.fHistoFile->GetName());

    if(oldname != newname) {
      delete fHistoFile;
      fHistoFile = new TFile(da.fHistoFile->GetName(),"update");
    }

    fMod = da.fMod;

    SetName(da.GetName());
    SetTitle(da.GetTitle());

    for(Int_t iX=0; iX<128; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {

	if(fCharge[iX][iZ]) delete fCharge[iX][iZ];
	fCharge[iX][iZ] = da.fCharge[iX][iZ];
      }
    }
    
  }

  return *this;
}


//-------------------------------------------------------------------
AliPHOSCpvDA1::~AliPHOSCpvDA1()
{
  // Destructor
  
  UpdateHistoFile();
  if(fHistoFile) delete fHistoFile;
  
}

//-------------------------------------------------------------------
void AliPHOSCpvDA1::FillHistograms(Float_t e[128][56]) 
{
  // Fill charge deposit histograms of one event.
  // Charge data is encoded as e[X][Z], 
  // where X(0..127) and Z(0..55) - pad position in the CPV module, 
  // Charge in ADC counts.
  // If no charge read for particular channel, 
  // the correspondent array entry should be filled by zero.
  // WARNING: this function should be called once per event!

  char hname[128];
  char htitl[128];
  
  for(Int_t iX=0; iX<128; iX++) {
    for (Int_t iZ=0; iZ<56; iZ++) {
      
      if(!e[iX][iZ]) continue;
      
      if(fCharge[iX][iZ]) 
	fCharge[iX][iZ]->Fill(e[iX][iZ]);
      else {
	snprintf(hname,128,"%d_%d_%d",fMod,iX,iZ);
	snprintf(htitl,128,"Charge deposited on the pad %d_%d_%d",fMod,iX,iZ);
	fCharge[iX][iZ] = new TH1F(hname,htitl,1024,0.,1024.);
	fCharge[iX][iZ]->Fill(e[iX][iZ]);
      }
      
    }
  }
  
}

//-------------------------------------------------------------------
void AliPHOSCpvDA1::UpdateHistoFile()
{
  // Write histograms to file

  if(!fHistoFile) return;
  if(!fHistoFile->IsOpen()) return;

  TH1F* hist1=0;

  for(Int_t iX=0; iX<128; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {
      
      hist1 = fCharge[iX][iZ];
      if(hist1) hist1->Write(hist1->GetName(),TObject::kWriteDelete);
    }
  }

}

