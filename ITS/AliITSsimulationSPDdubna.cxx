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

#include <Riostream.h>
#include <TRandom.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliITS.h"
#include "AliITShit.h"
#include "AliITSdigit.h"
#include "AliITSmodule.h"
#include "AliITSMapA2.h" 
#include "AliITSpList.h"
#include "AliITSsimulationSPDdubna.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSresponseSPDdubna.h"

//#define DEBUG

ClassImp(AliITSsimulationSPDdubna)
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Boris Batyunya
// December 20 1999
//
// AliITSsimulationSPDdubna is the simulation of SPDs
//______________________________________________________________________


AliITSsimulationSPDdubna::AliITSsimulationSPDdubna(){
    // constructor

    fResponse = 0;
    fSegmentation = 0;
    fMapA2 = 0;
    fpList = 0;
    fModule = 0;
    fEvent = 0;
    fHis = 0;
    fNoise = 0.;
    fBaseline = 0.;
    fNPixelsZ = 0;
    fNPixelsX = 0;
}
//______________________________________________________________________
AliITSsimulationSPDdubna::AliITSsimulationSPDdubna(AliITSsegmentation *seg,
						   AliITSresponse *resp){
    // standard constructor
    const Double_t kmictocm = 1.0e-4; // convert microns to cm.

    fHis = 0;
    fResponse = resp;
    fSegmentation = seg;
    fModule = 0;
    fEvent = 0;

    fNPixelsZ=GetSeg()->Npz();
    fNPixelsX=GetSeg()->Npx();

    GetResp()->GetNoiseParam(fNoise,fBaseline);
    GetResp()->SetDistanceOverVoltage(kmictocm*GetSeg()->Dy(),50.0);

//    fMapA2 = new AliITSMapA2(GetSeg());
    fMapA2 = 0;

    fpList = new AliITSpList(fNPixelsZ+1,fNPixelsX+1);

}
//______________________________________________________________________
AliITSsimulationSPDdubna::~AliITSsimulationSPDdubna(){
    // destructor

    if(fMapA2) delete fMapA2;

    if (fHis) {
	fHis->Delete(); 
	delete fHis;     
    } // end if fHis
}
//______________________________________________________________________
AliITSsimulationSPDdubna::AliITSsimulationSPDdubna(const 
						   AliITSsimulationSPDdubna 
						   &source):
    AliITSsimulation(source){
    //     Copy Constructor 
    if(&source == this) return;
    this->fMapA2 = source.fMapA2;
    this->fNoise = source.fNoise;
    this->fBaseline = source.fBaseline;
    this->fNPixelsX = source.fNPixelsX;
    this->fNPixelsZ = source.fNPixelsZ;
    this->fHis = source.fHis;
    return;
}
//______________________________________________________________________
AliITSsimulationSPDdubna&  AliITSsimulationSPDdubna::operator=(const 
                                           AliITSsimulationSPDdubna &source){
    //    Assignment operator
    if(&source == this) return *this;
    this->fMapA2 = source.fMapA2;
    this->fNoise = source.fNoise;
    this->fBaseline = source.fBaseline;
    this->fNPixelsX = source.fNPixelsX;
    this->fNPixelsZ = source.fNPixelsZ;
    this->fHis = source.fHis;
    return *this;
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::InitSimulationModule(Int_t module, Int_t event){
    //  This function creates maps to build the list of tracks for each
    //  summable digit.
    //
    //  Inputs:
    //    Int_t module   // Module number to be simulated
    //    Int_t event    // Event number to be simulated
    //
    //  Outputs:
    //    none
    //
    //  Returns:
    //    none

    fModule = module;
    fEvent  = event;
//    fMapA2->ClearMap();
    fpList->ClearMap();
}
//_____________________________________________________________________
void AliITSsimulationSPDdubna::SDigitiseModule(AliITSmodule *mod, Int_t mask,
					       Int_t event){
    //  This function begins the work of creating S-Digits
    //
    //  Inputs:
    //    AliITSmodule *mod  //  module
    //    Int_t mask         //  mask to be applied to the module
    //
    //  Outputs:
    //    none
    //
    //  Return:
    //    test              //  test returns kTRUE if the module contained hits
    //                      //  test returns kFALSE if it did not contain hits

    Int_t module = 0;

    event = 0; // remove unused variable warning.
    if(!(mod->GetNhits())) return;// if module has no hits don't create Sdigits
    fModule = mod->GetIndex();
    HitToSDigit(mod, module, mask, fpList);
    WriteSDigits(fpList);
//    fMapA2->ClearMap();
    fpList->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::WriteSDigits(AliITSpList *pList){
    //  This function adds each S-Digit to pList
    //
    //  Inputs:
    //    AliITSpList *pList
    //
    //  Outputs:
    //    none
    //
    //  Return:
    //    none
    Int_t ix, nix, iz, niz;
    static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");

    pList->GetMaxMapIndex(niz, nix);
    for(iz=0; iz<niz-1; iz++)for(ix=0; ix<nix-1; ix++){
	if(pList->GetSignalOnly(iz+1,ix+1)>0.0){
	    aliITS->AddSumDigit(*(pList->GetpListItem(iz+1,ix+1)));
#ifdef DEBUG
	    cout <<"SDigits " << iz << "," << ix << "," << 
		*(pList->GetpListItem(iz+1,ix+1)) << endl;
#endif
	} // end if pList
    } // end for iz,ix
    return; 
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::FinishSDigitiseModule(){
    //  This function calls SDigitsToDigits which creates Digits from SDigits
    //
    //  Inputs:
    //    none
    //
    //  Outputs:
    //    none
    //  Return
    //    none

    SDigitsToDigits(fModule, fpList);
    return;
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::SDigitsToDigits(Int_t module,
					       AliITSpList *pList){
    //  This function adds electronic noise to the S-Digits and then adds them
    // to  a new pList
    //
    //  Inputs:
    //    Int_t       module  // module number
    //    AliITSpList *pList  // pList
    //
    //  Outputs:
    //    pList is passed along to the functions ChargeToSignal and GetList
    //
    //  Return:
    //    none

    fModule = module;
    ChargeToSignal(pList); // Charge To Signal both adds noise and
//    fMapA2->ClearMap();
    pList->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::DigitiseModule(AliITSmodule *mod, Int_t module,
					      Int_t dummy){
    //  This function creates Digits straight from the hits and then adds
    //  electronic noise to the digits before adding them to pList
    //
    //  Inputs:
    //    AliITSmodule *mod    // module
    //    Int_t        module  // module number  Dummy.
    //    Int_t        dummy
    //
    //  Outputs:
    //    Each of the input variables is passed along to HitToSDigit
    //
    //  Return:
    //    none

    fModule = module = mod->GetIndex();//This calls the module for HitToSDigit
    HitToSDigit(mod,fModule, dummy, fpList);
    ChargeToSignal(fpList);
//    fMapA2->ClearMap();
    fpList->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::UpdateMapSignal(Int_t iz, Int_t ix, Int_t trk,
					       Int_t ht, Int_t module,
					       Double_t signal,
					       AliITSpList *pList){
    //  This function adds a signal to the pList from the pList class
    //
    //  Inputs:
    //    Int_t       iz     // row number
    //    Int_t       ix     // column number
    //    Int_t       trk    // track number
    //    Int_t       ht     // hit number
    //    Double_t    signal // signal strength
    //    AliITSpList *pList // pList
    //
    //  Outputs:
    //    All of the inputs are passed to AliITSpList::AddSignal
    //    Int_t    ix  // row number
    //    Int_t    iz  // column number
    //    Double_t sig // signal strength
    //          // These three variables are defined to preserve the
    //          // assignments used in the function AliITSMapA2::AddSignal
    //
    //  Return:
    //    none

//    fMapA2->AddSignal(iz, ix, signal);
    module = fModule; // remove unused variable warning.
    pList->AddSignal(iz+1,ix+1, trk, ht, fModule, signal);
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::UpdateMapNoise(Int_t iz,
					      Int_t ix, Int_t fModule,
					      Double_t sig, Float_t noise,
					      AliITSpList *pList){
    //  This function adds noise to data in the MapA2 as well as the pList
    //
    //  Inputs:
    //    Int_t       iz       // row number
    //    Int_t       ix       // column number
    //    Int_t       mod     // module number
    //    Double_t    sig     // signal strength
    //    Double_t    noise   // electronic noise generated by ChargeToSignal
    //    AliITSpList *pList  // pList
    //
    //  Outputs:
    //    All of the inputs are passed to AliITSMapA2::AddSignal or
    //    AliITSpList::AddNoise
    //
    //  Return:
    //    none

//    fMapA2->AddSignal(iz, ix, noise);
    sig = 0.0; // remove unused variable warning.
    pList->AddNoise(iz+1,ix+1, fModule, noise);
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::HitToDigit(AliITSmodule *mod, Int_t module,
					  Int_t dummy){
    DigitiseModule(mod, module, dummy);
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::HitToSDigit(AliITSmodule *mod, Int_t module,
					    Int_t dummy,AliITSpList *pList){
    // Does the charge distributions using Gaussian diffusion charge charing.
    const Double_t kmictocm = 1.0e-4; // convert microns to cm.
    TObjArray *hits = mod->GetHits();
    Int_t nhits = hits->GetEntriesFast();
    Int_t h,ix,iz;
    Int_t idtrack;
    Double_t x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0,de=0.0;
    Double_t x,y,z,t,tp,st,dt=0.2,el,sig;
    Double_t thick = kmictocm*GetSeg()->Dy();

    module = dummy = pList->GetNEnteries(); // remove unused varuable warning.
    if(nhits<=0) return;
    for(h=0;h<nhits;h++){
#ifdef DEBUG
	cout << "Hits=" << h << "," << *(mod->GetHit(h)) << endl;
#endif
	if(mod->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,idtrack)){
	st =TMath::Sqrt(x1*x1+y1*y1+z1*z1);
	if(st>0.0){
	    st = (Double_t)((Int_t)(1.0E+04*st)); // number of microns
	    if(st<=0.0) st = 1.0;
	    dt = 1.0/st;
	    for(t=0;t<1.0;t+=dt){ // Integrate over t
		tp = t+0.5*dt;
		el = GetResp()->GeVToCharge((Float_t)(dt*de));
#ifdef DEBUG
		if(el<=0.0) cout << "el="<<el<<" dt="<<dt<<" de="<<de<<endl;
#endif
		x = x0+x1*tp;
		y = y0+y1*tp;
		z = z0+z1*tp;
		GetSeg()->LocalToDet(x,z,ix,iz);
		sig = GetResp()->SigmaDiffusion1D(thick + y);
		SpreadCharge(x,y,z,ix,iz,el,sig,idtrack,
			     mod->GetHitTrackIndex(h),h,mod->GetIndex());
	    } // end for t
	} else { // st == 0.0 deposit it at this point
	    el = GetResp()->GeVToCharge((Float_t)de);
	    x = x0;
	    y = y0;
	    z = z0;
	    GetSeg()->LocalToDet(x,z,ix,iz);
	    sig = GetResp()->SigmaDiffusion1D(thick + y);
	    SpreadCharge(x,y,z,ix,iz,el,sig,
			 idtrack,mod->GetHitTrackIndex(h),h,mod->GetIndex());
	} // end if st>0.0
    }} // Loop over all hits h
}/*
//______________________________________________________________________
void AliITSsimulationSPDdubna::HitToSDigit(AliITSmodule *mod, Int_t module,
					    Int_t dummy,AliITSpList *pList){
    // Does the charge distributions using Gaussian diffusion charge charing.
    const Double_t kmictocm = 1.0e-4; // convert microns to cm.
    TObjArray *hits = mod->GetHits();
    Int_t nhits = hits->GetEntriesFast();
    Int_t h,ix,iz,i,n;
    Int_t idtrack;
    Double_t x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0,de=0.0;
    Double_t x,y,z,*ta,t,tp,st,dt=0.2,el,sig;
    Double_t thick = kmictocm*GetSeg()->Dy();

    if(nhits<=0) return;
    for(h=0;h<nhits;h++){
#ifdef DEBUG
	cout << "Hits=" << h << "," << *(mod->GetHit(h)) << endl;
#endif
	if(mod->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,idtrack)){
	st =TMath::Sqrt(x1*x1+y1*y1+z1*z1);
	if(st>0.0){
	    st =TMath::Sqrt(x1*x1+y1*y1+z1*z1)*(ta[i+1]-ta[i]);
	    ta = CreateFindCellEdges(x0,x1,z0,z1,n);
	    for(i=0;i<n-1;i++){
		dt = TMath::Min((1.0E-4)/st,);
		for(t=ta[i];t<ta[i+1];t+=dt){ // Integrate over t
		tp = t+0.5*dt;
		el = GetResp()->GeVToCharge((Float_t)(dt*de));
#ifdef DEBUG
		if(el<=0.0) cout << "el="<<el<<" dt="<<dt<<" de="<<de<<endl;
#endif
		x = x0+x1*tp;
		y = y0+y1*tp;
		z = z0+z1*tp;
		GetSeg()->LocalToDet(x,z,ix,iz);
		sig = GetResp()->SigmaDiffusion1D(thick + y);
		SpreadCharge(x,y,z,ix,iz,el,sig,idtrack,
			     mod->GetHitTrackIndex(h),h,mod->GetIndex());
	    } // end for t[i]
	    delete[] t;
	} else { // st == 0.0 deposit it at this point
	    el = GetResp()->GeVToCharge((Float_t)de);
	    x = x0;
	    y = y0;
	    z = z0;
	    GetSeg()->LocalToDet(x,z,ix,iz);
	    sig = GetResp()->SigmaDiffusion1D(thick + y);
	    SpreadCharge(x,y,z,ix,iz,el,sig,
			 idtrack,mod->GetHitTrackIndex(h),h,mod->GetIndex());
	} // end if st>0.0
    }} // Loop over all hits h
    }*/
//______________________________________________________________________
void AliITSsimulationSPDdubna::SpreadCharge(Double_t x0,Double_t y0,
					    Double_t z0,Int_t ix0,Int_t iz0,
					    Double_t el,Double_t sig,Int_t t,
					    Int_t ti,Int_t hi,Int_t mod){
    // Spreads the charge over neighboring cells. Assume charge is distributed
    // as charge(x,z) = (el/2*pi*sig*sig)*exp(-arg)
    // arg=((x-x0)*(x-x0)/2*sig*sig)+((z-z0*z-z0)/2*sig*sig)
    // Defined this way, the integral over all x and z is el.
    const Int_t knx = 3,knz = 2;
    const Double_t kRoot2 = 1.414213562; // Sqrt(2).
    const Double_t kmictocm = 1.0e-4; // convert microns to cm.
    Int_t ix,iz,ixs,ixe,izs,ize;
    Float_t x,z;
    Double_t x1,x2,z1,z2,s,sp;

    y0 = ti; // remove unused variable warning.
    if(sig<=0.0) {
	fpList->AddSignal(iz0+1,ix0+1,t,hi,mod,el);
	return;
    } // end if
    sp = 1.0/(sig*kRoot2);
#ifdef DEBUG
    cout << "sig=" << sig << " sp=" << sp << endl;
#endif
    ixs = TMath::Max(-knx+ix0,0);
    ixe = TMath::Min(knx+ix0,GetSeg()->Npx()-1);
    izs = TMath::Max(-knz+iz0,0);
    ize = TMath::Min(knz+iz0,GetSeg()->Npz()-1);
    for(ix=ixs;ix<=ixe;ix++) for(iz=izs;iz<=ize;iz++){
	GetSeg()->DetToLocal(ix,iz,x,z); // pixel center
	x1 = x;
	z1 = z;
	x2  = x1 + 0.5*kmictocm*GetSeg()->Dpx(ix); // Upper
	x1 -= 0.5*kmictocm*GetSeg()->Dpx(ix);  // Lower
	z2  = z1 + 0.5*kmictocm*GetSeg()->Dpz(iz); // Upper
	z1 -= 0.5*kmictocm*GetSeg()->Dpz(iz);  // Lower
	x1 -= x0; // Distance from where track traveled
	x2 -= x0; // Distance from where track traveled
	z1 -= z0; // Distance from where track traveled
	z2 -= z0; // Distance from where track traveled
	s = 0.25; // Correction based on definision of Erfc
	s *= TMath::Erfc(sp*x1) - TMath::Erfc(sp*x2);
#ifdef DEBUG
	cout << "el=" << el << " ix0=" << ix0 << " ix=" << ix << " x0="<< x <<
	    " iz0=" << iz0 << " iz=" << iz << " z0=" << z  << 
	    " sp*x1=" << sp*x1 <<" sp*x2=" << sp*x2 << " s=" << s;
#endif
	s *= TMath::Erfc(sp*z1) - TMath::Erfc(sp*z2);
#ifdef DEBUG
	cout << " sp*z1=" << sp*z1 <<" sp*z2=" << sp*z2 << " s=" << s << endl;
#endif
	fpList->AddSignal(iz+1,ix+1,t,hi,mod,s*el);
    } // end for ix, iz
}
//______________________________________________________________________
Double_t *AliITSsimulationSPDdubna::CreateFindCellEdges(Double_t x0,Double_t x1,
                                             Double_t z0,Double_t z1,Int_t &n){
    // Note: This function is a potensial source for a memory leak. The memory
    // pointed to in its return, must be deleted.
    // Inputs:
    //    Double_t x0   The starting location of the track step in x
    //    Double_t x1   The distance allong x for the track step
    //    Double_t z0   The starting location of the track step in z
    //    Double_t z1   The distance allong z for the track step
    // Output:
    //    Int)t &n      The size of the array returned. Minimal n=2.
    // Return:
    //    The pointer to the array of track steps.
    Int_t ix0,ix1,ix,iz0,iz1,iz,i;
    Double_t x,z,lx,ux,lz,uz,a,b,c,d;
    Double_t *t;

    GetSeg()->LocalToDet(x0,z0,ix0,iz0);
    GetSeg()->LocalToDet(x1,z1,ix1,iz1);
    n = 2 + TMath::Abs(ix1-ix0) + TMath::Abs(iz1-iz0);
    t = new Double_t[n];
    t[0] = 0.0;
    t[n-1] = 1.0;
    x = x0;
    z = z0;
    for(i=1;i<n-1;i++){
	GetSeg()->LocalToDet(x,z,ix,iz);
	GetSeg()->CellBoundries(ix,iz,lx,ux,lz,uz);
	a = (lx-x0)/x1;
	if(a<=t[i-1]) a = 1.0;
	b = (ux-x0)/x1;
	if(b<=t[i-1]) b = 1.0;
	c = (lz-z0)/z1;
	if(c<=t[i-1]) c = 1.0;
	d = (uz-z0)/z1;
	if(d<=t[i-1]) d = 1.0;
	t[i] = TMath::Min(TMath::Min(TMath::Min(a,b),c),d);
	x = x0+x1*(t[i]*1.00000001);
	z = z0+z1*(t[i]*1.00000001);
	i++;
    } // end for i
    return t;
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::HitToSDigitOld(AliITSmodule *mod, Int_t module,
                                           Int_t dummy, AliITSpList *pList){
    // digitize module 
    const Float_t kEnToEl = 2.778e+8; // GeV->charge in electrons 
                                      // for 3.6 eV/pair 
    const Float_t kconv = 10000.;     // cm -> microns

    Float_t spdLength = GetSeg()->Dz();
    Float_t spdWidth = GetSeg()->Dx();
    Float_t spdThickness = GetSeg()->Dy();
    Float_t difCoef, dum;       
    GetResp()->DiffCoeff(difCoef,dum); 
    if(spdThickness > 290) difCoef = 0.00613;  

    Float_t zPix0 = 1e+6;
    Float_t xPix0 = 1e+6;
    Float_t yPrev = 1e+6;   

    Float_t zPitch = GetSeg()->Dpz(0);
    Float_t xPitch = GetSeg()->Dpx(0);
  
    TObjArray *fHits = mod->GetHits();
    module = dummy = mod->GetIndex();
    Int_t nhits = fHits->GetEntriesFast();
    if (!nhits) return;
#ifdef DEBUG
    cout<<"len,wid,thickness,nx,nz,pitchx,pitchz,difcoef ="<<spdLength<<","
	<<spdWidth<<","<<spdThickness<<","<<fNPixelsX<<","<<fNPixelsZ<<","
	<<xPitch<<","<<zPitch<<","<<difCoef<<endl;
#endif
    //  Array of pointers to the label-signal list
    Int_t indexRange[4] = {0,0,0,0};

    // Fill detector maps with GEANT hits
    // loop over hits in the module
    static Bool_t first;
    Int_t lasttrack=-2;
    Int_t hit, iZi, jz, jx;
    Int_t idhit=-1; //!
#ifdef DEBUG
    cout<<"SPDdubna: module,nhits ="<<module<<","<<nhits<<endl;
#endif
    for (hit=0;hit<nhits;hit++) {
        AliITShit *iHit = (AliITShit*) fHits->At(hit);
#ifdef DEBUG
	cout << "Hits=" << hit << "," << *iHit << endl;
#endif
	//Int_t layer = iHit->GetLayer();
        Float_t yPix0 = -spdThickness/2; 

	// work with the idtrack=entry number in the TreeH
	//Int_t idhit,idtrack; //!
	//mod->GetHitTrackAndHitIndex(hit,idtrack,idhit);  //!    
	//Int_t idtrack=mod->GetHitTrackIndex(hit);  
        // or store straight away the particle position in the array
	// of particles : 
	if(iHit->StatusEntering()) idhit=hit;
        Int_t itrack = iHit->GetTrack();
        Int_t dray = 0;
   
	if (lasttrack != itrack || hit==(nhits-1)) first = kTRUE; 

	//Int_t parent = iHit->GetParticle()->GetFirstMother();
        Int_t partcode = iHit->GetParticle()->GetPdgCode();

	//  partcode (pdgCode): 11 - e-, 13 - mu-, 22 - gamma, 111 - pi0,
	// 211 - pi+,  310 - K0s, 321 - K+, 2112 - n, 2212 - p, 3122 - lambda

	Float_t pmod = iHit->GetParticle()->P(); // total momentum at the
	                                           // vertex
        pmod *= 1000;

        if(partcode == 11 && pmod < 6) dray = 1; // delta ray is e-
                                                 // at p < 6 MeV/c

	//  Get hit z and x(r*phi) cordinates for each module (detector)
	//  in local system.

	Float_t zPix = kconv*iHit->GetZL();
	Float_t xPix = kconv*iHit->GetXL();
	Float_t yPix = kconv*iHit->GetYL();

	// Get track status
	Int_t status = iHit->GetTrackStatus();      

	// Check boundaries
	if(zPix  > spdLength/2) {
#ifdef DEBUG
	    cout<<"!!! SPD: z outside ="<<zPix<<endl;
#endif
	    zPix = spdLength/2 - 10;
	}
	if(zPix  < 0 && zPix < -spdLength/2) {
#ifdef DEBUG
	    cout<<"!!! SPD: z outside ="<<zPix<<endl;
#endif
	    zPix = -spdLength/2 + 10;
	}
	if(xPix  > spdWidth/2) {
#ifdef DEBUG
	    cout<<"!!! SPD: x outside ="<<xPix<<endl;
#endif
	    xPix = spdWidth/2 - 10;
	}
	if(xPix  < 0 && xPix < -spdWidth/2) {
#ifdef DEBUG
	    cout<<"!!! SPD: x outside ="<<xPix<<endl;
#endif
	    xPix = -spdWidth/2 + 10;
	}
	Int_t trdown = 0;

	// enter Si or after event in Si
	if (status == 66 ) {  
	    zPix0 = zPix;
	    xPix0 = xPix;
	    yPrev = yPix; 
	} // end if status == 66

	Float_t depEnergy = iHit->GetIonization();
	// skip if the input point to Si       

	if(depEnergy <= 0.) continue;        

	// if track returns to the opposite direction:
	if (yPix < yPrev) {
            trdown = 1;
	} // end if yPix < yPrev

	// take into account the holes diffusion inside the Silicon
	// the straight line between the entrance and exit points in Si is
	// divided into the several steps; the diffusion is considered 
	// for each end point of step and charge
	// is distributed between the pixels through the diffusion.

	//  ---------- the diffusion in Z (beam) direction -------
	Float_t charge = depEnergy*kEnToEl;         // charge in e-
	Float_t drPath = 0.;   
	Float_t tang = 0.;
	Float_t sigmaDif = 0.; 
	Float_t zdif = zPix - zPix0;
	Float_t xdif = xPix - xPix0;
	Float_t ydif = TMath::Abs(yPix - yPrev);
	Float_t ydif0 = TMath::Abs(yPrev - yPix0);

	if(ydif < 1) continue; // ydif is not zero

	Float_t projDif = sqrt(xdif*xdif + zdif*zdif);

	Int_t ndZ = (Int_t)TMath::Abs(zdif/zPitch) + 1;
	Int_t ndX = (Int_t)TMath::Abs(xdif/xPitch) + 1; 

	// number of the steps along the track:
	Int_t nsteps = ndZ;
	if(ndX > ndZ) nsteps = ndX;
	if(nsteps < 20) nsteps = 20;  // minimum number of the steps 

	if (projDif < 5 ) {
	    drPath = (yPix-yPix0)*1.e-4;  
	    drPath = TMath::Abs(drPath);        // drift path in cm
	    sigmaDif = difCoef*sqrt(drPath);    // sigma diffusion in cm
	    sigmaDif = sigmaDif*kconv;         // sigma diffusion in microns
	    nsteps = 1;
	}  // end if projDif < 5

	if(projDif > 5) tang = ydif/projDif;
	Float_t dCharge = charge/nsteps;       // charge in e- for one step
	Float_t dZ = zdif/nsteps;
	Float_t dX = xdif/nsteps;

	for (iZi = 1; iZi <= nsteps;iZi++) {
	    Float_t dZn = iZi*dZ;
	    Float_t dXn = iZi*dX;
	    Float_t zPixn = zPix0 + dZn;
	    Float_t xPixn = xPix0 + dXn;

	    if(projDif >= 5) {
		Float_t dProjn = sqrt(dZn*dZn+dXn*dXn);
		drPath = dProjn*tang*1.e-4; // drift path for iZi+1 step in cm 
		if(trdown == 0) {
		    drPath = TMath::Abs(drPath) + ydif0*1.e-4;
		}// end if trdow ==0
		if(trdown == 1) {
		    drPath = ydif0*1.e-4 - TMath::Abs(drPath);
		    drPath = TMath::Abs(drPath);
		} // end if trdown == 1
		sigmaDif = difCoef*sqrt(drPath);    
		sigmaDif = sigmaDif*kconv;       // sigma diffusion in microns
	    } // end if projdif >= 5

	    zPixn = (zPixn + spdLength/2.);  
	    xPixn = (xPixn + spdWidth/2.);  
            Int_t nZpix, nXpix;
            GetSeg()->GetPadIxz(xPixn,zPixn,nXpix,nZpix);
	    zPitch = GetSeg()->Dpz(nZpix);
            GetSeg()->GetPadTxz(xPixn,zPixn);
	    // set the window for the integration
	    Int_t jzmin = 1;  
	    Int_t jzmax = 3; 
	    if(nZpix == 1) jzmin =2;
	    if(nZpix == fNPixelsZ) jzmax = 2; 

	    Int_t jxmin = 1;  
	    Int_t jxmax = 3; 
	    if(nXpix == 1) jxmin =2;
	    if(nXpix == fNPixelsX) jxmax = 2; 

	    Float_t zpix = nZpix; 
	    Float_t dZright = zPitch*(zpix - zPixn);
	    Float_t dZleft = zPitch - dZright;

	    Float_t xpix = nXpix; 
	    Float_t dXright = xPitch*(xpix - xPixn);
	    Float_t dXleft = xPitch - dXright;

	    Float_t dZprev = 0.;
	    Float_t dZnext = 0.;
	    Float_t dXprev = 0.;
	    Float_t dXnext = 0.;

	    for(jz=jzmin; jz <=jzmax; jz++) {
	        if(jz == 1) {
		    dZprev = -zPitch - dZleft;
		    dZnext = -dZleft;
		} else if(jz == 2) {
		    dZprev = -dZleft;
		    dZnext = dZright;
		} else if(jz == 3) {
		    dZprev = dZright;
		    dZnext = dZright + zPitch;
		} // end if jz
		// kz changes from 1 to the fNofPixels(270)  
		Int_t kz = nZpix + jz -2; 

		Float_t zArg1 = dZprev/sigmaDif;
		Float_t zArg2 = dZnext/sigmaDif;
		Float_t zProb1 = TMath::Erfc(zArg1);
		Float_t zProb2 = TMath::Erfc(zArg2);
		Float_t dZCharge =0.5*(zProb1-zProb2)*dCharge; 


		// ----------- holes diffusion in X(r*phi) direction  --------

		if(dZCharge > 1.) { 
		    for(jx=jxmin; jx <=jxmax; jx++) {
			if(jx == 1) {
			    dXprev = -xPitch - dXleft;
			    dXnext = -dXleft;
			} else if(jx == 2) {
			    dXprev = -dXleft;
			    dXnext = dXright;
			} else if(jx == 3) {
			    dXprev = dXright;
			    dXnext = dXright + xPitch;
			}  // end if jx
			Int_t kx = nXpix + jx -2;
			Float_t xArg1 = dXprev/sigmaDif;
			Float_t xArg2 = dXnext/sigmaDif;
			Float_t xProb1 = TMath::Erfc(xArg1);
			Float_t xProb2 = TMath::Erfc(xArg2);
			Float_t dXCharge =0.5*(xProb1-xProb2)*dZCharge;

			if(dXCharge > 1.) {
			    if (first) {
				indexRange[0]=indexRange[1]=kz-1;
				indexRange[2]=indexRange[3]=kx-1;
				first=kFALSE;
			    } // end if first
			    indexRange[0]=TMath::Min(indexRange[0],kz-1);
			    indexRange[1]=TMath::Max(indexRange[1],kz-1);
			    indexRange[2]=TMath::Min(indexRange[2],kx-1);
			    indexRange[3]=TMath::Max(indexRange[3],kx-1);
/*
			    // build the list of digits for this module	
			    Double_t signal = fMapA2->GetSignal(kz-1,kx-1);
			    signal+=dXCharge;
			    fMapA2->SetHit(kz-1,kx-1,(double)signal);
*/
			    // The calling sequence for UpdateMapSignal was 
			    // moved into the (dx > 1 e-) loop because it 
			    // needs to call signal which is defined inside 
			    // this loop
			    fModule   = module;//Defined because functions 
			                       // called by UpdateMapSignal 
			                       // expect module to be an 
			                       // integer
			    UpdateMapSignal(kz-1,kx-1,
//					    mod->GetHitTrackIndex(hit),
                             ((AliITShit*)(mod->GetHit(hit)))->GetTrack(),
					    hit,fModule,dXCharge,pList);
			}      // dXCharge > 1 e-
		    }       // jx loop
		}       // dZCharge > 1 e-
	    }        // jz loop
	}         // iZi loop
        if (status == 65) {   // the step is inside of Si
	    zPix0 = zPix;
	    xPix0 = xPix;
	} // end if status == 65
	yPrev = yPix;
    }   // hit loop inside the module
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::ChargeToSignal(AliITSpList *pList){
    // add noise and electronics, perform the zero suppression and add the
    // digit to the list
    static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");
    Float_t threshold = (float)GetResp()->MinVal();
    Int_t j;
//    Int_t    digits[3], tracks[3], hits[3];
//    Float_t  charges[3];
    Float_t  electronics;
//    Float_t  phys; 
    Double_t sig;
    const Int_t    nmaxtrk=AliITSdigitSPD::GetNTracks();
    static AliITSdigitSPD dig;

    for(Int_t iz=0; iz<fNPixelsZ; iz++){
	for(Int_t ix=0; ix<fNPixelsX; ix++){
	    electronics = fBaseline + fNoise*gRandom->Gaus();
	    sig = pList->GetSignalOnly(iz+1,ix+1);
	    UpdateMapNoise(iz,ix,fModule,sig,electronics,pList);
#ifdef DEBUG
//	    cout << sig << "+" << electronics <<">threshold=" << threshold 
//		 << endl;
#endif
	    if (sig+electronics > threshold) {
		dig.SetCoord1(iz);
		dig.SetCoord2(ix);
		dig.SetSignal(1);
		Int_t sigspd = (Int_t) pList->GetSignal(iz+1,ix+1);
		dig.SetSignalSPD(sigspd);
		for(j=0;j<nmaxtrk;j++){
//		    charges[j] = 0.0;
		    if (j<pList->GetNEnteries()) {
			dig.SetTrack(j,pList->GetTrack(iz+1,ix+1,j));
			dig.SetHit(j,pList->GetHit(iz+1,ix+1,j));
		    }else { // Default values
			dig.SetTrack(j,-3);
			dig.SetHit(j,-3);
		    } // end if pList
		} // end for j
//		charges[0] = (Float_t) pList->GetSumSignal(iz+1,ix+1);
/*
		if(tracks[0] == tracks[1] && tracks[0] == tracks[2]) {
		    tracks[1] = -3;
		    hits[1] = -1;
		    tracks[2] = -3;
		    hits[2] = -1;
		} else if(tracks[0] == tracks[1] && tracks[0] != tracks[2]) {
		    tracks[1] = -3;
		    hits[1] = -1;   
		} else if(tracks[0] == tracks[2] && tracks[0] != tracks[1]) {
		    tracks[2] = -3;
		    hits[2] = -1;   
		} else if(tracks[1] == tracks[2] && tracks[0] != tracks[1]) {
		    tracks[2] = -3;
		    hits[2] = -1;   
		} // end if
*/
//		phys = 0.0;
#ifdef DEBUG
		cout << iz << "," << ix << "," << 
		    *(pList->GetpListItem(iz+1,ix+1)) << endl;
#endif
//		aliITS->AddSimDigit(0, phys, digits, tracks, hits, charges);
		aliITS->AddSimDigit(0,&dig);
	    } // 
	} // 
    } //
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::CreateHistograms(){
    // create 1D histograms for tests

    printf("SPD - create histograms\n");

    fHis=new TObjArray(fNPixelsZ);
    TString spdName("spd_");
    for (Int_t i=0;i<fNPixelsZ;i++) {
	Char_t pixelz[4];
	sprintf(pixelz,"%d",i+1);
	spdName.Append(pixelz);
	//PH	   (*fHis)[i] = new TH1F(spdName.Data(),"SPD maps",
	//PH                       fNPixelsX,0.,(Float_t) fNPixelsX);
	fHis->AddAt(new TH1F(spdName.Data(),"SPD maps",
			     fNPixelsX,0.,(Float_t) fNPixelsX), i);
    } // end for i
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::ResetHistograms(){
    //
    // Reset histograms for this detector
    //

    for ( int i=0;i<fNPixelsZ;i++ ) {
	//PH	if ((*fHis)[i])    ((TH1F*)(*fHis)[i])->Reset();
	if (fHis->At(i))    ((TH1F*)fHis->At(i))->Reset();
    } // end for i
}
