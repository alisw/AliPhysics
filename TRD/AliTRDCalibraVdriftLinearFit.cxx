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
#include <TStyle.h>
#include <TCanvas.h>
#include <TTreeStream.h>
#include <TGraphErrors.h>
#include <TDirectory.h>
#include <TTreeStream.h>
#include <TF1.h>

//header file
#include "AliTRDCalibraVdriftLinearFit.h"

ClassImp(AliTRDCalibraVdriftLinearFit) /*FOLD00*/

//_____________________________________________________________________
AliTRDCalibraVdriftLinearFit::AliTRDCalibraVdriftLinearFit() : /*FOLD00*/
  TObject(),
  fVersion(0),
  fNameCalibUsed(""),
  fLinearFitterHistoArray(540),
  fLinearFitterPArray(540),
  fLinearFitterEArray(540),
  fRobustFit(kTRUE),
  fDebugStreamer(0x0),
  fDebugLevel(0),
  fSeeDetector(0)
{
  //
  // default constructor
  //
}
//_____________________________________________________________________
AliTRDCalibraVdriftLinearFit::AliTRDCalibraVdriftLinearFit(const AliTRDCalibraVdriftLinearFit &ped) : /*FOLD00*/
  TObject(ped),
  fVersion(ped.fVersion),
  fNameCalibUsed(ped.fNameCalibUsed),
  fLinearFitterHistoArray(540),
  fLinearFitterPArray(540),
  fLinearFitterEArray(540),
  fRobustFit(kTRUE),
  fDebugStreamer(0x0),
  fDebugLevel(0),
  fSeeDetector(0)
{
    //
    // copy constructor
    //
  for (Int_t idet = 0; idet < 540; idet++){
   
    const TVectorD     *vectorE     = (TVectorD*)ped.fLinearFitterEArray.UncheckedAt(idet);
    const TVectorD     *vectorP     = (TVectorD*)ped.fLinearFitterPArray.UncheckedAt(idet);
    const TH2S         *hped        = (TH2S*)ped.fLinearFitterHistoArray.UncheckedAt(idet);
    
    if ( vectorE != 0x0 ) fLinearFitterEArray.AddAt(new TVectorD(*vectorE), idet);
    if ( vectorP != 0x0 ) fLinearFitterPArray.AddAt(new TVectorD(*vectorP), idet);
    if ( hped != 0x0 ){
      TH2S *hNew = (TH2S *)hped->Clone();
      //hNew->SetDirectory(0);
      fLinearFitterHistoArray.AddAt(hNew,idet);
    }
  }
}
//_____________________________________________________________________
AliTRDCalibraVdriftLinearFit::AliTRDCalibraVdriftLinearFit(const TObjArray &obja) : /*FOLD00*/
  TObject(),
  fVersion(0),
  fNameCalibUsed(""),
  fLinearFitterHistoArray(540),
  fLinearFitterPArray(540),
  fLinearFitterEArray(540),
  fRobustFit(kTRUE),
  fDebugStreamer(0x0),
  fDebugLevel(0),
  fSeeDetector(0)
{
  //
  // constructor from a TObjArray
  //
  for (Int_t idet = 0; idet < 540; idet++){
    const TH2S         *hped        = (TH2S*)obja.UncheckedAt(idet);
    if ( hped != 0x0 ){
      TH2S *hNew = (TH2S *)hped->Clone();
      //hNew->SetDirectory(0);
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
  fLinearFitterHistoArray.SetOwner();
  fLinearFitterPArray.SetOwner();
  fLinearFitterEArray.SetOwner();

  fLinearFitterHistoArray.Delete();
  fLinearFitterPArray.Delete();
  fLinearFitterEArray.Delete();

  if ( fDebugStreamer ) delete fDebugStreamer;

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
      TH2S *hped1 = (TH2S *)target.GetLinearFitterHisto(idet,kTRUE);
      //hped1->SetDirectory(0);
      hped1->Add((const TH2S *)fLinearFitterHistoArray.UncheckedAt(idet));
    }
  }
  
  TObject::Copy(c);

}
//_____________________________________________________________________________
Long64_t AliTRDCalibraVdriftLinearFit::Merge(const TCollection* list) 
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
	  TH2S *hped1 = (TH2S *)GetLinearFitterHisto(idet,kTRUE);
	  Double_t entriesa = hped1->GetEntries();
	  Double_t entriesb = ((TH2S *)entry->GetLinearFitterHisto(idet))->GetEntries();
	  if((entriesa + entriesb) < 5*32767) hped1->Add(entry->GetLinearFitterHisto(idet));
	}
      }
      
      count++;
    }
  
  return count;
}
//_____________________________________________________________________
void AliTRDCalibraVdriftLinearFit::Add(const AliTRDCalibraVdriftLinearFit *ped)
{
  //
  // Add histo
  //

  fVersion++;

  for (Int_t idet = 0; idet < 540; idet++){
    const TH2S         *hped        = (TH2S*)ped->GetLinearFitterHistoNoForce(idet);
    //printf("idet %d\n",idet);
    if ( hped != 0x0 ){
      //printf("add\n");
      TH2S *hped1 = (TH2S *)GetLinearFitterHisto(idet,kTRUE);
      Double_t entriesa = hped1->GetEntries();
      Double_t entriesb = hped->GetEntries();
      if((entriesa + entriesb) < 5*32767) hped1->Add(hped);
    }
  }
}
//______________________________________________________________________________________
TH2S* AliTRDCalibraVdriftLinearFit::GetLinearFitterHisto(Int_t detector, Bool_t force)
{
    //
    // return pointer to TH2F histo 
    // if force is true create a new histo if it doesn't exist allready
    //
    if ( !force || fLinearFitterHistoArray.UncheckedAt(detector) )
	return (TH2S*)fLinearFitterHistoArray.UncheckedAt(detector);

    return GetLinearFitterHistoForce(detector);

}
//______________________________________________________________________________________
TH2S* AliTRDCalibraVdriftLinearFit::GetLinearFitterHistoForce(Int_t detector)
{
  //
  // return pointer to TH2F histo 
  // if NULL create a new histo if it doesn't exist allready
  //
  if (fLinearFitterHistoArray.UncheckedAt(detector))
    return (TH2S*)fLinearFitterHistoArray.UncheckedAt(detector);
  
  // if we are forced and TLinearFitter doesn't yes exist create it
  
  // new TH2F
  TString name("LFDV");
  name += detector;
  name += "version";
  name +=  fVersion;
  
  TH2S *lfdv = new TH2S((const Char_t *)name,(const Char_t *) name
			,36,-0.9,0.9,48
			,-1.2,1.2);
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
  
  TH2S *h = ((TH2S *) GetLinearFitterHisto(detector,kTRUE));
  Double_t nbentries = h->GetEntries();
  if(nbentries < 5*32767) h->Fill(tnp,pars1);

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
    const TH2S *linearfitterhisto = (TH2S*)fLinearFitterHistoArray.UncheckedAt(cb);
    //printf("Processing the detector cb %d we find %d\n",cb, (Bool_t) linearfitterhisto);    

    if ( linearfitterhisto != 0 ){
      
      // Fill a linearfitter
      TAxis *xaxis = linearfitterhisto->GetXaxis();
      TAxis *yaxis = linearfitterhisto->GetYaxis();
      TLinearFitter linearfitter = TLinearFitter(2,"pol1");
      //printf("test\n");
      Double_t integral = linearfitterhisto->Integral();
      //printf("Integral is %f\n",integral);
      Bool_t securitybreaking = kFALSE;
      if(TMath::Abs(integral-1199) < 0.00001) securitybreaking = kTRUE;
      for(Int_t ibinx = 0; ibinx < linearfitterhisto->GetNbinsX(); ibinx++){
	for(Int_t ibiny = 0; ibiny < linearfitterhisto->GetNbinsY(); ibiny++){
	  if(linearfitterhisto->GetBinContent(ibinx+1,ibiny+1)>0){
	    Double_t x = xaxis->GetBinCenter(ibinx+1);
	    Double_t y = yaxis->GetBinCenter(ibiny+1);
	    
	    for(Int_t k = 0; k < (Int_t)linearfitterhisto->GetBinContent(ibinx+1,ibiny+1); k++){
	      if(!securitybreaking){
		linearfitter.AddPoint(&x,y);
		arrayI[cb]++;
	      }
	      else {
		if(arrayI[cb]< 1198){
		  linearfitter.AddPoint(&x,y);
		  arrayI[cb]++; 
		}
	      }
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
	//printf("Fit\n");
	if(((fRobustFit) && (linearfitter.EvalRobust(0.8)==0)) || ((!fRobustFit) && (linearfitter.Eval()==0))) {
	  //if((linearfitter.Eval()==0)) {
	  //printf("Take the param\n");
	  linearfitter.GetParameters(*par);
	  //printf("Done\n");
	  //linearfitter.GetErrors(pare);
	  //Float_t  ppointError =  TMath::Sqrt(TMath::Abs(linearfitter.GetChisquare())/arrayI[cb]);
	  //(*parE)[0] = pare[0]*ppointError;
	  //(*parE)[1] = pare[1]*ppointError;
	  
	  (*parE)[0] = 0.0;
	  (*parE)[1] = 0.0;
	  (*parE)[2] = (Double_t) arrayI[cb];
	  fLinearFitterPArray.AddAt(par,cb);
	  fLinearFitterEArray.AddAt(parE,cb);
	  
	  //par->Print();
	  //parE->Print();
	}
	//printf("Finish\n");
      }
      
      //delete linearfitterhisto;
      
    }// if something

  }

  delete [] arrayI;
   
}

//____________Functions fit Online CH2d________________________________________
void AliTRDCalibraVdriftLinearFit::FillPEArray2()
{
  //
  // Fill fFitterPArray and fFitterEArray from inside
  //

  // Loop over histos 
  TF1 *f1 = new TF1("f1","[0]+[1]*x",-1,1);
      
  for(Int_t cb = 0; cb < 540; cb++){
    //const TH2S *fitterhisto = (TH2S*)fLinearFitterHistoArray.UncheckedAt(cb);
    TH2S *fitterhisto = (TH2S*)fLinearFitterHistoArray.UncheckedAt(cb);
    //printf("Processing the detector cb %d we find %d\n",cb, (Bool_t) fitterhisto);    
  
    if ( fitterhisto != 0 ){
      
      Int_t nEntries=0;
      TGraphErrors *gg=DrawMS(fitterhisto,nEntries);
      // Number of points of the TGraphErrors
      if(gg->GetN() < 20) {
	if(gg) delete gg;
	continue;	
      }
      //printf("det %d, number of points %d, nEntries %d\n",cb,gg->GetN(),nEntries);
      gg->Fit(f1,"Q0");
      
      TVectorD  *par  = new TVectorD(2);
      TVectorD  *parE = new TVectorD(3);
      (*parE)[0] = 0.0;
      (*parE)[1] = 0.0;
      (*parE)[2] = (Double_t) nEntries;
      (*par)[0] = f1->GetParameter(0);
      (*par)[1] = f1->GetParameter(1);
      fLinearFitterPArray.AddAt(par,cb);
      fLinearFitterEArray.AddAt(parE,cb);
            
      if(fDebugLevel==0) {
	if(gg) delete gg;
      }
      else {
	if(cb==fSeeDetector) {
	  gStyle->SetPalette(1);
	  gStyle->SetOptStat(1111);
	  gStyle->SetPadBorderMode(0);
	  gStyle->SetCanvasColor(10);
	  gStyle->SetPadLeftMargin(0.13);
	  gStyle->SetPadRightMargin(0.01);
	  TCanvas *stat = new TCanvas("stat","",50,50,600,800);
	  stat->cd(1);
	  fitterhisto->Draw("colz");
	  gg->Draw("P");
	  f1->Draw("same");
	  break;
	}
      } 
    }
  }
  if(fDebugLevel==0) delete f1;
}

//_________Helper function__________________________________________________
TGraphErrors* AliTRDCalibraVdriftLinearFit::DrawMS(const TH2 *const h2, Int_t &nEntries)
{
  TF1 fg("fg", "gaus", -10., 30.);
  TGraphErrors *gp = new TGraphErrors();

  TAxis *ax(h2->GetXaxis());
  TAxis *ay(h2->GetYaxis());
  TH1D *h1(NULL);
  for(Int_t ipt(0), jpt(1), ig(0); ipt<ax->GetNbins(); ipt++, jpt++){
    h1 = h2->ProjectionY("py", jpt, jpt);
    fg.SetParameter(1, h1->GetMean());
    //Float_t x(ax->GetBinCenter(jpt)); 
    Int_t n=Int_t(h1->Integral(1, h1->GetNbinsX()));
    nEntries+=n;
    if(n < 15){
      //Warning("drawMS()", Form("reject x[%d]=%f on n=%d", jpt, x, n));
      continue;
    }
    h1->Fit(&fg, "WWQ0");
    if(fg.GetNDF()<2){
      //Warning("drawMS()", Form("reject x[%d]=%f on NDF=%d", jpt, x, fg.GetNDF()));
      continue;
    }
    if(((fg.GetParameter(1)+fg.GetParameter(2)/2)>ay->GetXmax()) || ((fg.GetParameter(1)-fg.GetParameter(2)/2)<ay->GetXmin()) || (TMath::Abs(fg.GetParameter(0))< 0.00001)) continue;
    gp->SetPoint(ig, ax->GetBinCenter(jpt), fg.GetParameter(1));
    gp->SetPointError(ig, 0, TMath::Sqrt(pow(fg.GetParError(1),2) + (1/pow(fg.GetParameter(0),2))));
    ig++;
  }
  delete h1;
  return gp;
}
