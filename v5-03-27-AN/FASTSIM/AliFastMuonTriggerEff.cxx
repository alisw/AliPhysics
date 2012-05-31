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
// The trigger parametrization is computed for background levels 0., 0.5 and 1.
// In order to set a background level different from 0 it is necessary to 
// explicitly force it with:
// ForceBkgLevel(BkgLevel).
// For intermediate background levels, the trigger response is linearly 
// interpolated between these values.
// There is increased granularity in the pT region below 3 GeV. Although
// it does not seem to be necessary it is also possible to interpolate
// between pT bins using SetInt().
// Author: Pietro Cortese (Universita' del Piemonte Orientale - Alessandria 
// and INFN of Torino)


#include "AliFastMuonTriggerEff.h"
#include "TROOT.h"
#include "TFile.h"
#include "stdlib.h"
#include "TH3.h"
#include "TObjString.h"

#define PLIN printf("%s: %d: ",__FILE__,__LINE__)

ClassImp(AliFastMuonTriggerEff)

AliFastMuonTriggerEff::AliFastMuonTriggerEff():
    AliFastResponse("Efficiency", "Muon Trigger Efficiency"),
    fPtMin(0.),
    fPtMax(0.),
    fDpt(0.),
    fnptb(0),
    fPhiMin(0.),
    fPhiMax(0.),
    fDphi(0.),
    fnphib(0),
    fThetaMin(0.),
    fThetaMax(0.),
    fDtheta(0.),
    fnthetab(0),
    fCut(kLow),
    fZones(0),
    fhEffAPt(0),
    fhEffLPt(0),
    fhEffHPt(0),
    fhLX(0),
    fhLY(0),
    fhLZ(0),
    fBkg(0.),
    fTableTitle(0),
    fDescription(0),
    fInt(0),
    fibx(0),
    fiby(0),
    fibz(0)
{
//
// Default constructor
//
}

AliFastMuonTriggerEff::AliFastMuonTriggerEff(const char* Name, const char* Title):
    AliFastResponse(Name, Title),
    fPtMin(0.),
    fPtMax(0.),
    fDpt(0.),
    fnptb(0),
    fPhiMin(0.),
    fPhiMax(0.),
    fDphi(0.),
    fnphib(0),
    fThetaMin(0.),
    fThetaMax(0.),
    fDtheta(0.),
    fnthetab(0),
    fCut(kLow),
    fZones(0),
    fhEffAPt(0),
    fhEffLPt(0),
    fhEffHPt(0),
    fhLX(0),
    fhLY(0),
    fhLZ(0),
    fBkg(0.),
    fTableTitle(0),
    fDescription(0),
    fInt(0),
    fibx(0),
    fiby(0),
    fibz(0)
{
// Another constructor
}

AliFastMuonTriggerEff::AliFastMuonTriggerEff(const AliFastMuonTriggerEff& eff)
    :AliFastResponse(eff),
    fPtMin(0.),
    fPtMax(0.),
    fDpt(0.),
    fnptb(0),
    fPhiMin(0.),
    fPhiMax(0.),
    fDphi(0.),
    fnphib(0),
    fThetaMin(0.),
    fThetaMax(0.),
    fDtheta(0.),
    fnthetab(0),
    fCut(kLow),
    fZones(0),
    fhEffAPt(0),
    fhEffLPt(0),
    fhEffHPt(0),
    fhLX(0),
    fhLY(0),
    fhLZ(0),
    fBkg(0.),
    fTableTitle(0),
    fDescription(0),
    fInt(0),
    fibx(0),
    fiby(0),
    fibz(0)
{
// Copy constructor
    eff.Copy(*this);
}

void AliFastMuonTriggerEff::SetCut(Int_t cut) 
{  
  //
  // Set the pt cut
  if(cut==kLow){
    printf("Selecting Low Pt cut\n");
  }else if(cut==kHigh){
    printf("Selecting High Pt cut\n");
  }else if(cut==kAny){
    printf("Selecting Lowest Pt cut\n");
  }else{
    printf("Don't know cut %d! Selecting Low Pt cut\n",cut);
    cut=kLow;
  }
  fCut = cut;
  
}

Int_t AliFastMuonTriggerEff::SetBkgLevel(Float_t Bkg)
{
  //
  // Set the background level
  //
    if((Bkg!=0.)) {
       printf("%s: Warning: requested Bkg: %f\n",
       __FILE__,Bkg);
       fBkg=0.;
       printf("A consistent treatement of the trigger probability\n");
       printf("within the framework of the fast simulation requires\n");
       printf("requires background 0\n");
       printf("%s: fBkg: set to %f\n",
       __FILE__,fBkg);
    } else {
      fBkg=Bkg;
    }
    if(fZones!=0.) {
	Init();
    }
    return 0;
}

Int_t AliFastMuonTriggerEff::ForceBkgLevel(Float_t Bkg)
{
  //
  // Check and enforce consistency of the background level 
  // 
    if((Bkg!=0.)) {
       printf("%s: Warning: requested Bkg: %f\n",
       __FILE__,Bkg);
       printf("A consistent treatement of the trigger probability\n");
       printf("within the framework of the fast simulation\n");
       printf("requires background 0");
       printf("%s: Continue with fBkg: %f\n",
       __FILE__,Bkg);
    }
    fBkg=Bkg;
    if(fZones!=0.) {
	Init();
    }
    return 0;
}

Int_t AliFastMuonTriggerEff::LoadTables(const Char_t *namet){
  //
  // Load the trigger tables
  //
    Char_t hNameA[100],hNameL[100],hNameH[100];
    snprintf(hNameA, 100, "hEffAPt%s",namet);
    snprintf(hNameL, 100, "hEffLPt%s",namet);
    snprintf(hNameH, 100, "hEffHPt%s",namet);
    fhEffAPt = (TH3F*)gDirectory->Get(hNameA);
    fhEffLPt = (TH3F*)gDirectory->Get(hNameL);
    fhEffHPt = (TH3F*)gDirectory->Get(hNameH);
    if(!fhEffAPt){
      PLIN; printf("%s: histogram %s not found\n",__FILE__,hNameA);
      return -1;
    }
    if(!fhEffLPt){
      PLIN; printf("%s: histogram %s not found\n",__FILE__,hNameL);
      return -2;
    }
    if(!fhEffHPt){
      PLIN; printf("%s: histogram %s not found\n",__FILE__,hNameH);
      return -3;
    }
    return 0;
}

void AliFastMuonTriggerEff::Init()
{
//
//  Initialization
//
    fZones=0;
    Char_t file[100]="$(ALICE_ROOT)/FASTSIM/data/MUONtriggerLUT_V2.4nvdn.root";
    printf("Initializing %s / %s\n", fName.Data(), fTitle.Data());
    printf("using data from file: %s\n",file);
    printf("AliFastMuonTriggerEff: Initialization with background level: %f\n",fBkg);
    TFile *f = new TFile(file);
    if(f->IsZombie()) {
        PLIN; printf("Cannot open file: %s\n",file);
        return;
    }
    f->ls();
    Int_t intb=0;
    Char_t namet[10];
    if(TMath::Abs(fBkg)<0.00001){
      snprintf(namet, 10, "00");
    }else if(TMath::Abs(fBkg-0.5)<0.00001){
      snprintf(namet, 10, "05");
    }else if(TMath::Abs(fBkg-1.0)<0.00001){
      snprintf(namet, 10, "10");
    }else{
      PLIN; printf("A table for Bkg level: %f does not exists\n",fBkg);
      intb=1;
    }
    if(intb){ // Interpolation between background levels
      PLIN; printf("Interpolating Bkg level: %f\n",fBkg);
      TH3F* ha1,*hl1,*hh1,*ha2,*hl2,*hh2,*ha0,*hl0,*hh0;
      Char_t name1[10],name2[10]; Float_t b1,b2;
      if(fBkg>0&&fBkg<0.5){
        snprintf(name1,10, "00");
        snprintf(name2,10, "05");
	b1=0.;
	b2=0.5;
      }else if(fBkg>0.5){
        snprintf(name1, 10, "05");
        snprintf(name2, 10, "10");
	b1=0.5;
	b2=1.0;
	if(fBkg>1.0){
	  for(Int_t i=0; i<10;i++){
	    PLIN; printf("WARNING!!!! You are extrapolating above background 1.0\n");
	  }
	}
      }else{
        PLIN; printf("Bkg level: %f is not supported\n",fBkg);
        return;
      }
      if(LoadTables(name1)){
        PLIN; printf("Error in loading trigger tables\n");
	return;
      }
      PLIN; printf("We use tables for %f and %f to interpolate %f Bkg level\n",b1,b2,fBkg);
      ha0=(TH3F*)fhEffAPt->Clone("hEffAPtXX"); ha0->Reset();
      hl0=(TH3F*)fhEffLPt->Clone("hEffLPtXX"); hl0->Reset();
      hh0=(TH3F*)fhEffHPt->Clone("hEffHPtXX"); hh0->Reset();
      ha1=fhEffAPt;
      hl1=fhEffLPt;
      hh1=fhEffHPt;
      if(LoadTables(name2)){
        PLIN; printf("Error in loading trigger tables\n");
	return;
      }
      ha2=fhEffAPt;
      hl2=fhEffLPt;
      hh2=fhEffHPt;
      fhEffAPt=ha0;
      fhEffLPt=hl0;
      fhEffHPt=hh0;
      Int_t nnx=ha0->GetNbinsX()+1;
      Int_t nny=ha0->GetNbinsY()+1;
      Int_t nnz=ha0->GetNbinsZ()+1;
      for(Int_t ix=0; ix<=nnx; ix++){
	for(Int_t iy=0; iy<=nny; iy++){
	  for(Int_t iz=0; iz<=nnz; iz++){
	    Double_t y1,y2; Float_t cont;
	    y1=ha1->GetBinContent(ix,iy,iz); y2=ha2->GetBinContent(ix,iy,iz);
	    cont=Float_t(y1+(y2-y1)/(b2-b1)*(fBkg-b1)); if(cont>1)cont=1; if(cont<0)cont=0;
	    fhEffAPt->SetBinContent(ix,iy,iz,cont);
	    y1=hl1->GetBinContent(ix,iy,iz); y2=hl2->GetBinContent(ix,iy,iz);
	    cont=Float_t(y1+(y2-y1)/(b2-b1)*(fBkg-b1)); if(cont>1)cont=1; if(cont<0)cont=0;
	    fhEffLPt->SetBinContent(ix,iy,iz,cont);
	    y1=hh1->GetBinContent(ix,iy,iz); y2=hh2->GetBinContent(ix,iy,iz);
	    cont=Float_t(y1+(y2-y1)/(b2-b1)*(fBkg-b1)); if(cont>1)cont=1; if(cont<0)cont=0;
	    fhEffHPt->SetBinContent(ix,iy,iz,cont);
	  }
	}
      }
    }else{ // Use tables computed for selected backgound levels
      printf("Loading tables for background level: %f\n",fBkg);
      if(LoadTables(namet)){
        PLIN; printf("Error in loading trigger tables\n");
	return;
      }
    }
    fhEffAPt->SetDirectory(0);
    fhEffLPt->SetDirectory(0);
    fhEffHPt->SetDirectory(0);
    fhLX=fhEffLPt->GetXaxis();
    fhLY=fhEffLPt->GetYaxis();
    fhLZ=fhEffLPt->GetZaxis();
//
//
    if(f->Get("Description"))
    {
      fDescription=((TObjString*)f->Get("Description"))->GetString();
      printf("%s\n",fDescription.Data());
    }

    fThetaMin = fhEffLPt->GetXaxis()->GetXmin();
    fThetaMax = fhEffLPt->GetXaxis()->GetXmax();
    fnthetab=fhEffLPt->GetNbinsX();
    fDtheta   = (fThetaMax-fThetaMin)/fnthetab;

    fPhiMin   = fhEffLPt->GetYaxis()->GetXmin();
    fPhiMax   = fhEffLPt->GetYaxis()->GetXmax();
    fnphib=fhEffLPt->GetNbinsY();
    fDphi     = (fPhiMax-fPhiMin)/fnphib;

    fPtMin=fhEffLPt->GetZaxis()->GetXmin();
    fPtMax=fhEffLPt->GetZaxis()->GetXmax();
    fnptb=fhEffLPt->GetNbinsZ();
    fDpt      = (fPtMax-fPtMin)/fnptb;

    printf("***** This version of AliFastMuonTriggerEff can use both *****\n");
    printf("***** new and old ALICE reference frames depending on    *****\n");
    printf("***** which LUT has been loaded. You can find below some *****\n");
    printf("***** information on the current parametrization:        *****\n");
    printf("%4d bins in theta [%f:%f]\n",fnthetab,fThetaMin,fThetaMax);
    printf("%4d bins in phi [%f:%f]\n",fnphib,fPhiMin,fPhiMax);
    printf("%4d bins in pt [%f:%f]\n",fnptb,fPtMin,fPtMax);

    fZones=fnthetab*fnphib;

    f->Close();
    if(fInt==0) {
      printf("Interpolation of trigger efficiencies is off!\n");
    } else {
      printf("Interpolation of trigger efficiencies is on!\n");
    }
}

void AliFastMuonTriggerEff::Evaluate(Float_t charge, Float_t pt,Float_t theta,
              Float_t phi, Float_t& effLow, Float_t& effHigh, Float_t& effAny)
{
    //
    //  Trigger efficiency for pt, theta, phi (low, high and "any" cut)
    //
#ifdef MYTRIGDEBUG
    printf("Evaluate(ch=%2.0f, pt=%10.6f, theta=%7.2f, phi=%8.2f ...)\n",charge,pt,theta,phi);
#endif
    effLow=0.;
    effHigh=0.;
    effAny=0;
    if(fZones==0) {
        printf("Call to uninitialized object of class: AliFastMuonTriggerEff\n");
	return;
    }
    if(pt<0) {
        printf("Warning: pt: %f < 0. GeV/c\n",pt);
	return;	
    }

    Int_t iPt   = fhLZ->FindBin((Double_t)pt);
    if(iPt>fnptb)iPt=fnptb;
    Int_t iPhi  = Int_t((phi-fPhiMin)/fDphi);
    if(phi<fPhiMin)iPhi=iPhi-1;
    Int_t iTheta = fhLX->FindBin((Double_t)theta);
#ifdef MYTRIGDEBUG
    printf("Evaluate(ch=%2.0f, pt=%10.6f, theta=%7.2f, phi=%8.2f ...)\n",charge,pt,theta,phi);
    printf(" 0:%1d iPt iTheta iPhi: %d %d %d\n",fInt,iPt,iTheta,iPhi);
#endif
    iPhi=iPhi-2*fnphib*(iPhi/(2*fnphib));
#ifdef MYTRIGDEBUG
    printf(" 1:%1d iPhi converted to: %d for angle equivalence\n",fInt,iPhi);
#endif
    if(iPhi<0)iPhi=-iPhi-1;
    if(iPhi>(fnphib-1))iPhi=2*fnphib-1-iPhi;
#ifdef MYTRIGDEBUG
    printf(" 2:%1d iPhi converted to: %d for the symmetry of the spectrometer\n",fInt,iPhi);
#endif
    if(charge==1.){
    } else if(charge==-1.) {
    iPhi=fnphib-1-iPhi;
#ifdef MYTRIGDEBUG
    printf(" 3:%1d iPhi converted to: %d for the charge symmetry\n",fInt,iPhi);
#endif
    } else {
        printf("Warning: not understand charge: %f\n",charge);
        return;
    }
    if(iTheta<=0||iTheta>fnthetab) {
        printf("Warning: theta: %f outside acceptance\n",theta);
        return;
    }
    if(iPt<0) {
        printf("Warning: what do you mean with pt: %f <0?\n",pt);
	return;
    }
    iPhi++;
#ifdef MYTRIGDEBUG
    printf(" 4:%1d Getting: iTheta, iPhi, iPt: %d %d %d\n",
              fInt,iTheta,iPhi,iPt);
#endif
    effLow =fhEffLPt->GetBinContent(iTheta,iPhi,iPt);
    effHigh=fhEffHPt->GetBinContent(iTheta,iPhi,iPt);
    effAny =fhEffAPt->GetBinContent(iTheta,iPhi,iPt);
#ifdef MYTRIGDEBUG
    printf(" 4:%1d Result: charge, iTheta, iPhi, iPt: %f %d %d %d effLow: %f, effHigh: %f, effAny: %f\n",
              fInt,charge,iTheta,iPhi,iPt,effLow,effHigh,effAny);
#endif
        
    if(fInt==1) {
      Float_t angl,angh,anga;
      Float_t effLowp,effHighp,effAnyp;
      Float_t ptc=(iPt+0.5)*fDpt;  // The center of current bin
      #ifdef MYTRIGDEBUG
        printf(" 5:1 The center of current bin iPt: %d is: %f\n",iPt,ptc);
      #endif
      if(iPt==fnptb) {
        #ifdef MYTRIGDEBUG
	printf(" 6:1 No more points above! No interpolation is needed!\n");
        #endif
	return;        
      }else if(TMath::Abs(ptc-pt) < 1.e-10){
        #ifdef MYTRIGDEBUG
	  printf(" 6:1 No interpolation is needed!\n");
        #endif
	return;
      }else if(ptc>pt){
	// Looking for previous point
        if(iPt>1) {
          effLowp =fhEffLPt->GetBinContent(iTheta,iPhi,iPt-1);
          effHighp=fhEffHPt->GetBinContent(iTheta,iPhi,iPt-1);
          effAnyp =fhEffAPt->GetBinContent(iTheta,iPhi,iPt-1);
          #ifdef MYTRIGDEBUG
	  printf(" 7:1 A simple look to previous point: %d: %f %f\n",iPt-1,effLowp,effHighp);
          #endif
	} else {
	  effLowp=0.;
	  effHighp=0.;
          effAnyp=0;
          #ifdef MYTRIGDEBUG
          printf(" 8:1 result is: %f %f %f\n",effLowp,effHighp,effAnyp);
	  #endif	  
	}
        angl=(effLow-effLowp)/fDpt;
        angh=(effHigh-effHighp)/fDpt;
        anga=(effAny-effAnyp)/fDpt;
      }else{
	// Looking for next point
	if(iPt<fnptb) {
          effLowp =fhEffLPt->GetBinContent(iTheta,iPhi,iPt+1);
          effHighp=fhEffHPt->GetBinContent(iTheta,iPhi,iPt+1);
          effAnyp =fhEffAPt->GetBinContent(iTheta,iPhi,iPt+1);
          #ifdef MYTRIGDEBUG
	  printf(" 7:1 A simple look to next point: %d: %f %f %f\n",iPt-1,effLowp,effHighp,effAnyp);
          #endif
	} else {
	  effLowp=effLow;
	  effHighp=effHigh;
	  effAnyp=effAny;
          #ifdef MYTRIGDEBUG
          printf(" 8:1 result is: pt: %f %f %f\n",effLowp,effHighp,effAnyp);
	  #endif	  
        }
        angl=(effLowp-effLow)/fDpt;
        angh=(effHighp-effHigh)/fDpt;
        anga=(effAnyp-effAny)/fDpt;
      }
      effLow=effLow+angl*(pt-ptc);
      effHigh=effHigh+angh*(pt-ptc);
      effAny=effAny+anga*(pt-ptc);
      #ifdef MYTRIGDEBUG
      printf(" 9:1 the interpolation coefficients are: %f %f %f\n",angl,angh,anga);
      #endif
  }
  #ifdef MYTRIGDEBUG
  printf("10:%1d effLow, effHigh=%f %f %f\n",fInt,effLow,effHigh,effAny);
  #endif
  return;
}



Float_t AliFastMuonTriggerEff::Evaluate(Float_t charge, Float_t pt,
                   Float_t theta, Float_t phi)
{
    //
    // Trigger efficiency for pt, theta, phi depending of fCut
    // 
    if(fZones==0) {
        printf("Call to uninitialized object of class: AliFastMuonTriggerEff\n");
	return 0.;
    }
    Float_t eff;
    Float_t effLow, effHigh, effAny;
    
    Evaluate(charge,pt,theta,phi,effLow,effHigh,effAny);
    if (fCut == kLow) 
	eff  = effLow;
    else if (fCut == kHigh)
	eff  = effHigh;
    else if (fCut == kAny)
        eff  = effAny;
    else
        eff  = 0;

    return eff;
}

AliFastMuonTriggerEff& AliFastMuonTriggerEff::operator=(const  AliFastMuonTriggerEff& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

