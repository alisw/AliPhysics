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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTRDCalibraVdriftLinearFit                                           //
//                                                                        //
// Does the Vdrift an ExB calibration by applying a linear fit            //
//                                                                        //
// Author:                                                                //
//   R. Bailhache (R.Bailhache@gsi.de)                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

//Root includes
#include <TObjArray.h>
#include <TH2F.h>
#include <TString.h>
#include <TVectorD.h>
#include <TAxis.h>
#include <TLinearFitter.h>
#include <TMath.h>
#include <TTreeStream.h>

//header file
#include "AliTRDCalibraVdriftLinearFit.h"

ClassImp(AliTRDCalibraVdriftLinearFit) /*FOLD00*/

//_____________________________________________________________________
AliTRDCalibraVdriftLinearFit::AliTRDCalibraVdriftLinearFit() : /*FOLD00*/
  TObject(),
  fVersion(0),
  fLinearFitterHistoArray(540),
  fLinearFitterPArray(540),
  fLinearFitterEArray(540)
{
    //
    // default constructor
    //
}
//_____________________________________________________________________
AliTRDCalibraVdriftLinearFit::AliTRDCalibraVdriftLinearFit(const AliTRDCalibraVdriftLinearFit &ped) : /*FOLD00*/
  TObject(ped),
  fVersion(ped.fVersion),
  fLinearFitterHistoArray(540),
  fLinearFitterPArray(540),
  fLinearFitterEArray(540)
{
    //
    // copy constructor
    //
  for (Int_t idet = 0; idet < 540; idet++){
    const TVectorD     *vectorE     = (TVectorD*)ped.fLinearFitterEArray.UncheckedAt(idet);
    const TVectorD     *vectorP     = (TVectorD*)ped.fLinearFitterPArray.UncheckedAt(idet);
    const TH2F         *hped        = (TH2F*)ped.fLinearFitterHistoArray.UncheckedAt(idet);
    
    if ( vectorE != 0x0 ) fLinearFitterEArray.AddAt(new TVectorD(*vectorE), idet);
    if ( vectorP != 0x0 ) fLinearFitterPArray.AddAt(new TVectorD(*vectorP), idet);
    if ( hped != 0x0 ){
      TH2F *hNew = new TH2F(*hped);
      hNew->SetDirectory(0);
      fLinearFitterHistoArray.AddAt(hNew,idet);
    }
    
  }
}
//_____________________________________________________________________
AliTRDCalibraVdriftLinearFit::AliTRDCalibraVdriftLinearFit(const TObjArray &obja) : /*FOLD00*/
  TObject(),
  fVersion(0),
  fLinearFitterHistoArray(540),
  fLinearFitterPArray(540),
  fLinearFitterEArray(540)
{
  //
  // constructor from a TObjArray
  //
  for (Int_t idet = 0; idet < 540; idet++){
    const TH2F         *hped        = (TH2F*)obja.UncheckedAt(idet);
    if ( hped != 0x0 ){
      TH2F *hNew = new TH2F(*hped);
      hNew->SetDirectory(0);
      fLinearFitterHistoArray.AddAt(hNew,idet);
    }
  }
}
//_____________________________________________________________________
AliTRDCalibraVdriftLinearFit& AliTRDCalibraVdriftLinearFit::operator = (const  AliTRDCalibraVdriftLinearFit &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTRDCalibraVdriftLinearFit(source);

  return *this;
}
//_____________________________________________________________________
AliTRDCalibraVdriftLinearFit::~AliTRDCalibraVdriftLinearFit() /*FOLD00*/
{
  //
  // destructor
  //
  fLinearFitterHistoArray.Delete();
  fLinearFitterPArray.Delete();
  fLinearFitterEArray.Delete();

}
//_____________________________________________________________________________
void AliTRDCalibraVdriftLinearFit::Copy(TObject &c) const
{
  //
  // Copy function
  //

  AliTRDCalibraVdriftLinearFit& target = (AliTRDCalibraVdriftLinearFit &) c;

  // Copy only the histos
  for (Int_t idet = 0; idet < 540; idet++){
    if(fLinearFitterHistoArray.UncheckedAt(idet)){
      TH2F *hped1 = target.GetLinearFitterHisto(idet,kTRUE);
      hped1->SetDirectory(0);
      hped1->Add((const TH1F *)fLinearFitterHistoArray.UncheckedAt(idet));
    }
  
  }

  TObject::Copy(c);

}
//_____________________________________________________________________________
Long64_t AliTRDCalibraVdriftLinearFit::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;
  
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  
  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliTRDCalibraVdriftLinearFit* entry = dynamic_cast<AliTRDCalibraVdriftLinearFit*>(obj);
    if (entry == 0) continue; 

    // Copy only the histos
    for (Int_t idet = 0; idet < 540; idet++){
      if(entry->GetLinearFitterHisto(idet)){
	TH2F *hped1 = GetLinearFitterHisto(idet,kTRUE);
	hped1->SetDirectory(0);
	hped1->Add(entry->GetLinearFitterHisto(idet));
      }
      
    }
    
    count++;
  }
  
  return count;
}
//_____________________________________________________________________
void AliTRDCalibraVdriftLinearFit::Add(AliTRDCalibraVdriftLinearFit *ped)
{
  //
  // Add histo
  //

  fVersion++;

  for (Int_t idet = 0; idet < 540; idet++){
    const TH2F         *hped        = (TH2F*)ped->GetLinearFitterHisto(idet);
    //printf("idet %d\n",idet);
    if ( hped != 0x0 ){
      //printf("add\n");
      TH2F *hped1 = GetLinearFitterHisto(idet,kTRUE);
      //printf("test2\n");
      hped1->SetDirectory(0);
      //printf("test4\n");
      hped1->Add(hped);
      //printf("test3\n");
    }
  }
}
//______________________________________________________________________________________
TH2F* AliTRDCalibraVdriftLinearFit::GetLinearFitterHisto(Int_t detector, Bool_t force)
{
    //
    // return pointer to TH2F histo 
    // if force is true create a new histo if it doesn't exist allready
    //
    if ( !force || fLinearFitterHistoArray.UncheckedAt(detector) )
	return (TH2F*)fLinearFitterHistoArray.UncheckedAt(detector);

    // if we are forced and TLinearFitter doesn't yes exist create it

    // new TH2F
    TString name("LFDV");
    name += detector;
    name += "version";
    name +=  fVersion;
    
    TH2F *lfdv = new TH2F((const Char_t *)name,(const Char_t *) name
			  ,100,-1.0,1.0,100
			  ,-2.0,2.0);
    lfdv->SetXTitle("tan(phi_{track})");
    lfdv->SetYTitle("dy/dt");
    lfdv->SetZTitle("Number of clusters");
    lfdv->SetStats(0);
    lfdv->SetDirectory(0);

    fLinearFitterHistoArray.AddAt(lfdv,detector);
    return lfdv;
}
//______________________________________________________________________________________
Bool_t AliTRDCalibraVdriftLinearFit::GetParam(Int_t detector, TVectorD *param)
{
    //
    // return param for this detector
    //
  if ( fLinearFitterPArray.UncheckedAt(detector) ){
    const TVectorD     *vectorP     = (TVectorD*)fLinearFitterPArray.UncheckedAt(detector);
    if(!param) param = new TVectorD(2);
    for(Int_t k = 0; k < 2; k++){
      (*param)[k] = (*vectorP)[k];
    }
    return kTRUE;
  }
  else return kFALSE;

}
//______________________________________________________________________________________
Bool_t AliTRDCalibraVdriftLinearFit::GetError(Int_t detector, TVectorD *error)
{
    //
    // return error for this detector 
    //
  if ( fLinearFitterEArray.UncheckedAt(detector) ){
    const TVectorD     *vectorE     = (TVectorD*)fLinearFitterEArray.UncheckedAt(detector);
    if(!error) error = new TVectorD(3);
    for(Int_t k = 0; k < 3; k++){
      (*error)[k] = (*vectorE)[k];
    }
    return kTRUE;
  }
  else return kFALSE;

}
//______________________________________________________________________________________
void AliTRDCalibraVdriftLinearFit::Update(Int_t detector, Float_t tnp, Float_t pars1)
{
    //
    // Fill the 2D histos for debugging
    //

  ((TH2F *) GetLinearFitterHisto(detector,kTRUE))->Fill(tnp,pars1);

}
//____________Functions fit Online CH2d________________________________________
void AliTRDCalibraVdriftLinearFit::FillPEArray()
{
  //
  // Fill fLinearFitterPArray and fLinearFitterEArray from inside
  //

  
  Int_t *arrayI = new Int_t[540];
  for(Int_t k = 0; k< 540; k++){
    arrayI[k] = 0; 
  }

  // Loop over histos 
  for(Int_t cb = 0; cb < 540; cb++){
    const TH2F *linearfitterhisto = (TH2F*)fLinearFitterHistoArray.UncheckedAt(cb);
    //printf("Processing the detector cb %d we find %d\n",cb, (Bool_t) linearfitterhisto);    

    if ( linearfitterhisto != 0 ){
      
      // Fill a linearfitter
      TAxis *xaxis = linearfitterhisto->GetXaxis();
      TAxis *yaxis = linearfitterhisto->GetYaxis();
      TLinearFitter linearfitter = TLinearFitter(2,"pol1");
      //printf("test\n");
      for(Int_t ibinx = 0; ibinx < linearfitterhisto->GetNbinsX(); ibinx++){
	for(Int_t ibiny = 0; ibiny < linearfitterhisto->GetNbinsY(); ibiny++){
	  if(linearfitterhisto->GetBinContent(ibinx+1,ibiny+1)>0){
	    Double_t x = xaxis->GetBinCenter(ibinx+1);
	    Double_t y = yaxis->GetBinCenter(ibiny+1);
	    for(Int_t k = 0; k < (Int_t)linearfitterhisto->GetBinContent(ibinx+1,ibiny+1); k++){
	      linearfitter.AddPoint(&x,y);
	      arrayI[cb]++;
	    }
	  }
	}
      }
      
      //printf("Find %d entries for the detector %d\n",arrayI[cb],cb);

      // Eval the linear fitter
      if(arrayI[cb]>10){
	TVectorD  *par  = new TVectorD(2);
	TVectorD   pare = TVectorD(2);
	TVectorD  *parE = new TVectorD(3);
	linearfitter.Eval();
	linearfitter.GetParameters(*par);
	linearfitter.GetErrors(pare);
	Float_t  ppointError =  TMath::Sqrt(TMath::Abs(linearfitter.GetChisquare())/arrayI[cb]);
	(*parE)[0] = pare[0]*ppointError;
	(*parE)[1] = pare[1]*ppointError;
	(*parE)[2] = (Double_t) arrayI[cb];
	fLinearFitterPArray.AddAt(par,cb);
	fLinearFitterEArray.AddAt(parE,cb);
	//par->Print();
	//parE->Print();
      }    
    }// if something
  }
   
}
