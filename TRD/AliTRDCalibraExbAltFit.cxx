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

/* $Id: AliTRDCalibraExbAltFit.cxx 46327 2011-01-10 13:29:56Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTRDCalibraExbAltFit                                                 //
//                                                                        //
// Does the ExB calibration by applying a quadratic fit                   //
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
#include <TDirectory.h>
#include <TTreeStream.h>
#include <TGraphErrors.h>
#include <TF1.h>

//header file
#include "AliTRDCalibraExbAltFit.h"

ClassImp(AliTRDCalibraExbAltFit) /*FOLD00*/

//_____________________________________________________________________
AliTRDCalibraExbAltFit::AliTRDCalibraExbAltFit() : /*FOLD00*/
  TObject(),
  fVersion(0),
  fFitterHistoArray(540),
  fFitterPArray(540),
  fFitterEArray(540),
  fRobustFit(kFALSE),
  fDebugStreamer(0x0),
  fDebugLevel(0)
{
  //
  // default constructor
  //
}
//_____________________________________________________________________
AliTRDCalibraExbAltFit::AliTRDCalibraExbAltFit(const AliTRDCalibraExbAltFit &ped) : /*FOLD00*/
  TObject(ped),
  fVersion(ped.fVersion),
  fFitterHistoArray(540),
  fFitterPArray(540),
  fFitterEArray(540),
  fRobustFit(kFALSE),
  fDebugStreamer(0x0),
  fDebugLevel(0)
{
    //
    // copy constructor
    //
  for (Int_t idet = 0; idet < 540; idet++){
   
    const TVectorD     *vectorE     = (TVectorD*)ped.fFitterEArray.UncheckedAt(idet);
    const TVectorD     *vectorP     = (TVectorD*)ped.fFitterPArray.UncheckedAt(idet);
    const TH2S         *hped        = (TH2S*)ped.fFitterHistoArray.UncheckedAt(idet);
    
    if ( vectorE != 0x0 ) fFitterEArray.AddAt(new TVectorD(*vectorE), idet);
    if ( vectorP != 0x0 ) fFitterPArray.AddAt(new TVectorD(*vectorP), idet);
    if ( hped != 0x0 ){
      TH2S *hNew = (TH2S *)hped->Clone();
      //hNew->SetDirectory(0);
      fFitterHistoArray.AddAt(hNew,idet);
    }
  }
}
//_____________________________________________________________________
AliTRDCalibraExbAltFit::AliTRDCalibraExbAltFit(const TObjArray &obja) : /*FOLD00*/
  TObject(),
  fVersion(0),
  fFitterHistoArray(540),
  fFitterPArray(540),
  fFitterEArray(540),
  fRobustFit(kFALSE),
  fDebugStreamer(0x0),
  fDebugLevel(0)
{
  //
  // constructor from a TObjArray
  //
  for (Int_t idet = 0; idet < 540; idet++){
    const TH2S         *hped        = (TH2S*)obja.UncheckedAt(idet);
    if ( hped != 0x0 ){
      TH2S *hNew = (TH2S *)hped->Clone();
      //hNew->SetDirectory(0);
      fFitterHistoArray.AddAt(hNew,idet);
    }
  }
}
//_____________________________________________________________________
AliTRDCalibraExbAltFit& AliTRDCalibraExbAltFit::operator = (const  AliTRDCalibraExbAltFit &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTRDCalibraExbAltFit(source);

  return *this;
}
//_____________________________________________________________________
AliTRDCalibraExbAltFit::~AliTRDCalibraExbAltFit() /*FOLD00*/
{
  //
  // destructor
  //
  fFitterHistoArray.SetOwner();
  fFitterPArray.SetOwner();
  fFitterEArray.SetOwner();

  fFitterHistoArray.Delete();
  fFitterPArray.Delete();
  fFitterEArray.Delete();

  if ( fDebugStreamer ) delete fDebugStreamer;

}
//_____________________________________________________________________________
void AliTRDCalibraExbAltFit::Copy(TObject &c) const
{
  //
  // Copy function
  //

  AliTRDCalibraExbAltFit& target = (AliTRDCalibraExbAltFit &) c;

  // Copy only the histos
  for (Int_t idet = 0; idet < 540; idet++){
    if(fFitterHistoArray.UncheckedAt(idet)){
      TH2S *hped1 = (TH2S *)target.GetFitterHisto(idet,kTRUE);
      //hped1->SetDirectory(0);
      hped1->Add((const TH2S *)fFitterHistoArray.UncheckedAt(idet));
    }
  }
  
  TObject::Copy(c);

}
//_____________________________________________________________________________
Long64_t AliTRDCalibraExbAltFit::Merge(const TCollection* list) 
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
      AliTRDCalibraExbAltFit* entry = dynamic_cast<AliTRDCalibraExbAltFit*>(obj);
      if (entry == 0) continue; 
      
      // Copy only the histos
      for (Int_t idet = 0; idet < 540; idet++){
	if(entry->GetFitterHisto(idet)){
	  TH2S *hped1 = (TH2S *)GetFitterHisto(idet,kTRUE);
	  Double_t entriesa = hped1->GetEntries();
	  Double_t entriesb = ((TH2S *)entry->GetFitterHisto(idet))->GetEntries();
	  if((entriesa + entriesb) < 5*32767) hped1->Add(entry->GetFitterHisto(idet));
	}
      }
      
      count++;
    }
  
  return count;
}
//_____________________________________________________________________
void AliTRDCalibraExbAltFit::Add(const AliTRDCalibraExbAltFit *ped)
{
  //
  // Add histo
  //

  fVersion++;

  for (Int_t idet = 0; idet < 540; idet++){
    const TH2S         *hped        = (TH2S*)ped->GetFitterHistoNoForce(idet);
    //printf("idet %d\n",idet);
    if ( hped != 0x0 ){
      //printf("add\n");
      TH2S *hped1 = (TH2S *)GetFitterHisto(idet,kTRUE);
      Double_t entriesa = hped1->GetEntries();
      Double_t entriesb = hped->GetEntries();
      if((entriesa + entriesb) < 5*32767) hped1->Add(hped);
    }
  }
}
//______________________________________________________________________________________
TH2S* AliTRDCalibraExbAltFit::GetFitterHisto(Int_t detector, Bool_t force)
{
    //
    // return pointer to TH2F histo 
    // if force is true create a new histo if it doesn't exist allready
    //
    if ( !force || fFitterHistoArray.UncheckedAt(detector) )
	return (TH2S*)fFitterHistoArray.UncheckedAt(detector);

    return GetFitterHistoForce(detector);

}
//______________________________________________________________________________________
TH2S* AliTRDCalibraExbAltFit::GetFitterHistoForce(Int_t detector)
{
  //
  // return pointer to TH2F histo 
  // if NULL create a new histo if it doesn't exist allready
  //
  if (fFitterHistoArray.UncheckedAt(detector))
    return (TH2S*)fFitterHistoArray.UncheckedAt(detector);
  
  // if we are forced and TLinearFitter doesn't yes exist create it
  
  // new TH2F
  TString name("LFEXB");
  name += detector;
  name += "version";
  name +=  fVersion;
  
  TH2S *lfdv = new TH2S((const Char_t *)name,(const Char_t *) name,
			30, -TMath::DegToRad()*45, TMath::DegToRad()*45, 
			30, 0.3, 1.4);
  lfdv->SetXTitle("tan(phi_{track})");
  lfdv->SetYTitle("rms");
  lfdv->SetZTitle("Number of tracklets");
  lfdv->SetStats(0);
  lfdv->SetDirectory(0);
  
  fFitterHistoArray.AddAt(lfdv,detector);
  return lfdv;
}
//______________________________________________________________________________________
Bool_t AliTRDCalibraExbAltFit::GetParam(Int_t detector, TVectorD *param)
{
    //
    // return param for this detector
    //
  if ( fFitterPArray.UncheckedAt(detector) ){
    const TVectorD     *vectorP     = (TVectorD*)fFitterPArray.UncheckedAt(detector);
    if(!param) param = new TVectorD(vectorP->GetNoElements());
    for(Int_t k = 0; k < vectorP->GetNoElements(); k++){
      (*param)[k] = (*vectorP)[k];
    }
    return kTRUE;
  }
  else return kFALSE;

}
//______________________________________________________________________________________
Bool_t AliTRDCalibraExbAltFit::GetError(Int_t detector, TVectorD *error)
{
    //
    // return error for this detector 
    //
  if ( fFitterEArray.UncheckedAt(detector) ){
    const TVectorD     *vectorE     = (TVectorD*)fFitterEArray.UncheckedAt(detector);
    if(!error) error = new TVectorD(vectorE->GetNoElements());
    for(Int_t k = 0; k < vectorE->GetNoElements(); k++){
      (*error)[k] = (*vectorE)[k];
    }
    return kTRUE;
  }
  else return kFALSE;

}
//______________________________________________________________________________________
void AliTRDCalibraExbAltFit::Update(Int_t detector, Float_t tnp, Float_t pars1)
{
    //
    // Fill the 2D histos for debugging
    //
  
  TH2S *h = ((TH2S *) GetFitterHisto(detector,kTRUE));
  Double_t nbentries = h->GetEntries();
  if(nbentries < 5*32767) h->Fill(tnp,pars1);

}
//____________Functions fit Online CH2d________________________________________
void AliTRDCalibraExbAltFit::FillPEArray()
{
  //
  // Fill fFitterPArray and fFitterEArray from inside
  //

  
  Int_t *arrayI = new Int_t[540];
  for(Int_t k = 0; k< 540; k++){
    arrayI[k] = 0; 
  }

  // Loop over histos 
  for(Int_t cb = 0; cb < 540; cb++){
    const TH2S *fitterhisto = (TH2S*)fFitterHistoArray.UncheckedAt(cb);
    //printf("Processing the detector cb %d we find %d\n",cb, (Bool_t) fitterhisto);    

    if ( fitterhisto != 0 ){
      
      // Fill a fitter
      TAxis *xaxis = fitterhisto->GetXaxis();
      TAxis *yaxis = fitterhisto->GetYaxis();
      TLinearFitter fitter = TLinearFitter(3,"pol2");
      //printf("test\n");
      Double_t integral = fitterhisto->Integral();
      //printf("Integral is %f\n",integral);
      Bool_t securitybreaking = kFALSE;
      if(TMath::Abs(integral-1199) < 0.00001) securitybreaking = kTRUE;
      for(Int_t ibinx = 0; ibinx < fitterhisto->GetNbinsX(); ibinx++){
	for(Int_t ibiny = 0; ibiny < fitterhisto->GetNbinsY(); ibiny++){
	  if(fitterhisto->GetBinContent(ibinx+1,ibiny+1)>0){
	    Double_t x = xaxis->GetBinCenter(ibinx+1);
	    Double_t y = yaxis->GetBinCenter(ibiny+1);
	    
	    for(Int_t k = 0; k < (Int_t)fitterhisto->GetBinContent(ibinx+1,ibiny+1); k++){
	      if(!securitybreaking){
		fitter.AddPoint(&x,y);
		arrayI[cb]++;
	      }
	      else {
		if(arrayI[cb]< 1198){
		  fitter.AddPoint(&x,y);
		  arrayI[cb]++; 
		}
	      }
	    }
	    
	  }
	}
      }
      
      //printf("Find %d entries for the detector %d\n",arrayI[cb],cb);

      // Eval the  fitter
      if(arrayI[cb]>15){
	TVectorD  *par  = new TVectorD(3);
	TVectorD  pare = TVectorD(2);
	TVectorD  *parE = new TVectorD(3);
	//printf("Fit\n");
	//if((fitter.EvalRobust(0.8)==0)) {
	//if((fitter.EvalRobust()==0)) {
	if(((fRobustFit) && (fitter.EvalRobust(0.8)==0)) || ((!fRobustFit) && (fitter.Eval()==0))) {
	  //if((fitter.Eval()==0)) {
	  //printf("Take the param\n");
	  fitter.GetParameters(*par);
	  //printf("Done\n");
	  //fitter.GetErrors(*pare);
	  //Float_t  ppointError =  TMath::Sqrt(TMath::Abs(fitter.GetChisquare())/arrayI[cb]);
	  //(*parE)[0] = pare[0]*ppointError;
	  //(*parE)[1] = pare[1]*ppointError;
	  (*parE)[0] = 0.0;
	  (*parE)[1] = 0.0;
	  (*parE)[2] = (Double_t) arrayI[cb];
	  fFitterPArray.AddAt(par,cb);
	  fFitterEArray.AddAt(parE,cb);
	  
	  //par->Print();
	  //parE->Print();
	}
	//printf("Finish\n");
      }
      
      //delete fitterhisto;
      
    }// if something

  }

  delete [] arrayI;
   
}

//____________Functions fit Online CH2d________________________________________
void AliTRDCalibraExbAltFit::FillPEArray2()
{
  //
  // Fill fFitterPArray and fFitterEArray from inside
  //

  // Loop over histos 
      
  for(Int_t cb = 0; cb < 540; cb++){
    TH2S *fitterhisto = (TH2S*)fFitterHistoArray.UncheckedAt(cb);
    //printf("Processing the detector cb %d we find %d\n",cb, (Bool_t) fitterhisto);    

    if ( fitterhisto != 0 ){
      TF1 *f1 = new TF1("f1","[0]*(x-[1])**2+[2]",-1,1);
      
      Int_t nEntries=0;
      TGraphErrors *gg=DrawMS(fitterhisto,nEntries);
      // printf("N: %i\n",gg->GetN());
      // printf("entries: %i\n",nEntries);

      if(gg->GetN() < 20) {
	if(gg) delete gg;
	continue;	
      }
      gg->Fit(f1,"Q0");
      
      TVectorD  *par  = new TVectorD(3);
      TVectorD  *parE = new TVectorD(3);
      (*parE)[0] = 0.0;
      (*parE)[1] = 0.0;
      (*parE)[2] = (Double_t) nEntries;
      (*par)[0] = f1->GetParameter(0)*TMath::Power(f1->GetParameter(1),2)+f1->GetParameter(2);
      (*par)[1] = -2*f1->GetParameter(0)*f1->GetParameter(1);
      (*par)[2] = f1->GetParameter(0);
      fFitterPArray.AddAt(par,cb);
      fFitterEArray.AddAt(parE,cb);

      // if ( !fDebugStreamer ) {
      // 	//debug stream
      // 	TDirectory *backup = gDirectory;
      // 	fDebugStreamer = new TTreeSRedirector("TRDdebugCalibraFill.root");
      // 	if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
      // } 
    
      // if(cb==513){
      // 	fitterhisto->Draw();
      // 	gg->Draw("P");
      // 	f1->Draw("same");
      // }
      // double exb=f1->GetParameter(1);
      // (* fDebugStreamer) << "FindP1TrackPHtrackletV1"<<
      // 	//"snpright="<<snpright<<
      // 	"det="<<cb<<
      // 	"h2="<<(TObject*)fitterhisto<<
      // 	"gg="<<gg<<
      // 	"f1="<<f1<<
      // 	"exb="<<exb<<
      // 	"\n";

      //if(cb!=513){
      delete gg; 
      delete f1;
      //}
    }
  }
  //if(fDebugStreamer) delete fDebugStreamer;
}

//_________Helper function__________________________________________________
TGraphErrors* AliTRDCalibraExbAltFit::DrawMS(const TH2 *const h2, Int_t &nEntries)
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
    if(n < 5){
      //Warning("drawMS()", Form("reject x[%d]=%f on n=%d", jpt, x, n));
      continue;
    }
    h1->Fit(&fg, "WWQ0");
    if(fg.GetNDF()<2){
      //Warning("drawMS()", Form("reject x[%d]=%f on NDF=%d", jpt, x, fg.GetNDF()));
      continue;
    }
    if(fg.GetParameter(1)+fg.GetParameter(2)/2>ay->GetXmax() || fg.GetParameter(1)-fg.GetParameter(2)/2<ay->GetXmin() || !fg.GetParameter(0)) continue;
    gp->SetPoint(ig, ax->GetBinCenter(jpt), fg.GetParameter(1));
    gp->SetPointError(ig, 0, TMath::Sqrt(pow(fg.GetParError(1),2) + (1/pow(fg.GetParameter(0),2))));
    ig++;
  }
  delete h1;
  return gp;
}
