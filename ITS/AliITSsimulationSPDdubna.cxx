#include <iostream.h>
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
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"




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

    fHis = 0;
    fResponse = resp;
    fSegmentation = seg;
    fModule = 0;
    fEvent = 0;

    fNPixelsZ=fSegmentation->Npz();
    fNPixelsX=fSegmentation->Npx();

    fResponse->GetNoiseParam(fNoise,fBaseline);

    fMapA2 = new AliITSMapA2(fSegmentation);

    fpList = new AliITSpList(fNPixelsZ+1,fNPixelsX+1);

}
//______________________________________________________________________
AliITSsimulationSPDdubna::~AliITSsimulationSPDdubna(){
    // destructor

    delete fMapA2;

    if (fHis) {
	fHis->Delete(); 
	delete fHis;     
    } // end if fHis
}
//______________________________________________________________________
AliITSsimulationSPDdubna::AliITSsimulationSPDdubna(const 
						   AliITSsimulationSPDdubna 
						   &source){
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
    fMapA2->ClearMap();
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

    if(!(mod->GetNhits())) return;// if module has no hits don't create Sdigits
    fModule = mod->GetIndex();
    HitToSDigit(mod, module, mask, fpList);
    WriteSDigits(fpList);
    fMapA2->ClearMap();
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
    Int_t i, ni, j, nj;
    static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");

    pList->GetMaxMapIndex(ni, nj);
    for(i=0; i<ni; i++)for(j=0; j<nj; j++){
	if(pList->GetSignalOnly(i, j)>0.0){
	    aliITS->AddSumDigit(*(pList->GetpListItem(i, j)));
	} // end if pList
    } // end for i,j
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
    fMapA2->ClearMap();
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

    fModule = mod->GetIndex();  //This calls the module for HitToSDigit
    HitToSDigit(mod,fModule, dummy, fpList);
    ChargeToSignal(fpList);
    fMapA2->ClearMap();
    fpList->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::UpdateMapSignal(Int_t i, Int_t j, Int_t trk,
					       Int_t ht, Int_t module,
					       Double_t signal,
					       AliITSpList *pList){
    //  This function adds a signal to the pList from the pList class
    //
    //  Inputs:
    //    Int_t       i      // row number
    //    Int_t       j      // column number
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
    Int_t    iz = j;
    Int_t    ix = i;
    Double_t sig = signal;

    fMapA2->AddSignal(iz, ix, sig);
    pList->AddSignal(i, j, trk, ht, fModule, signal);
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::UpdateMapNoise(Int_t i, Int_t j, Int_t ix,
					      Int_t iz, Int_t fModule,
					      Double_t sig, Float_t noise,
					      AliITSpList *pList){
    //  This function adds noise to data in the MapA2 as well as the pList
    //
    //  Inputs:
    //    Int_t       i == ix // row number
    //    Int_t       j == iz // column number
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

    fMapA2->AddSignal(iz, ix, sig);
    pList->AddNoise(i, j, fModule, noise);
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::HitToDigit(AliITSmodule *mod, Int_t module,
					  Int_t dummy){
    DigitiseModule(mod, module, dummy);
}
//______________________________________________________________________
void AliITSsimulationSPDdubna::HitToSDigit(AliITSmodule *mod, Int_t module,
                                           Int_t dummy, AliITSpList *pList){
    // digitize module 
    const Float_t kEnToEl = 2.778e+8; // GeV->charge in electrons 
                                      // for 3.6 eV/pair 
    const Float_t kconv = 10000.;     // cm -> microns

    Float_t spdLength = fSegmentation->Dz();
    Float_t spdWidth = fSegmentation->Dx();
    Float_t spdThickness = fSegmentation->Dy();
    Float_t difCoef, dum;       
    fResponse->DiffCoeff(difCoef,dum); 
    if(spdThickness > 290) difCoef = 0.00613;  

    Float_t zPix0 = 1e+6;
    Float_t xPix0 = 1e+6;
    Float_t yPrev = 1e+6;   

    Float_t zPitch = fSegmentation->Dpz(0);
    Float_t xPitch = fSegmentation->Dpx(0);
  
    TObjArray *fHits = mod->GetHits();
    module = mod->GetIndex();
    Int_t nhits = fHits->GetEntriesFast();
    if (!nhits) return;

    cout<<"len,wid,thickness,nx,nz,pitchx,pitchz,difcoef ="<<spdLength<<","
	<<spdWidth<<","<<spdThickness<<","<<fNPixelsX<<","<<fNPixelsZ<<","
	<<xPitch<<","<<zPitch<<","<<difCoef<<endl;
    //  Array of pointers to the label-signal list
    Int_t indexRange[4] = {0,0,0,0};

    // Fill detector maps with GEANT hits
    // loop over hits in the module
    static Bool_t first;
    Int_t lasttrack=-2;
    Int_t hit, iZi, jz, jx;
    Int_t idhit=-1; //!
    cout<<"SPDdubna: module,nhits ="<<module<<","<<nhits<<endl;
    for (hit=0;hit<nhits;hit++) {
        AliITShit *iHit = (AliITShit*) fHits->At(hit);
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
	    //cout<<"!!! SPD: z outside ="<<zPix<<endl;
	    zPix = spdLength/2 - 10;
	}
	if(zPix  < 0 && zPix < -spdLength/2) {
	    //cout<<"!!! SPD: z outside ="<<zPix<<endl;
	    zPix = -spdLength/2 + 10;
	}
	if(xPix  > spdWidth/2) {
	    //cout<<"!!! SPD: x outside ="<<xPix<<endl;
	    xPix = spdWidth/2 - 10;
	}
	if(xPix  < 0 && xPix < -spdWidth/2) {
	    //cout<<"!!! SPD: x outside ="<<xPix<<endl;
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
            fSegmentation->GetPadIxz(xPixn,zPixn,nXpix,nZpix);
	    zPitch = fSegmentation->Dpz(nZpix);
            fSegmentation->GetPadTxz(xPixn,zPixn);
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
			    Int_t index = kz-1;
			    if (first) {
				indexRange[0]=indexRange[1]=index;
				indexRange[2]=indexRange[3]=kx-1;
				first=kFALSE;
			    } // end if first
			    indexRange[0]=TMath::Min(indexRange[0],kz-1);
			    indexRange[1]=TMath::Max(indexRange[1],kz-1);
			    indexRange[2]=TMath::Min(indexRange[2],kx-1);
			    indexRange[3]=TMath::Max(indexRange[3],kx-1);

			    // build the list of digits for this module	
			    Double_t signal=fMapA2->GetSignal(index,kx-1);
			    signal+=dXCharge;
			    fMapA2->SetHit(index,kx-1,(double)signal);

			    // The calling sequence for UpdateMapSignal was 
			    // moved into the (dx > 1 e-) loop because it 
			    // needs to call signal which is defined inside 
			    // this loop
			    Int_t i   = kx-1;
			    Int_t j   = kz-1;
			    Int_t trk = mod->GetHitTrackIndex(hit);
			    Int_t ht  = hit;
			    fModule   = module;//Defined because functions 
			                       // called by UpdateMapSignal 
			                       // expect module to be an 
			                       // integer
			    UpdateMapSignal(j,i,trk,ht,fModule,signal,pList);
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

    AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");  

    Float_t threshold = (float)fResponse->MinVal();

    Int_t    digits[3], tracks[3], hits[3], gi, j1;
    Float_t  charges[3];
    Float_t  electronics;
    Float_t  signal;
    Float_t  phys; 
    Double_t sig;
    Int_t    module = 0;
    for(Int_t iz=0; iz<fNPixelsZ; iz++){
	for(Int_t ix=0; ix<fNPixelsX; ix++){
	    electronics = fBaseline + fNoise*gRandom->Gaus();
	    signal = (float)pList->GetSignalOnly(ix, iz);
	    sig = Double_t (signal);  // sig will be passed along to 
	                              // UpdateMapNoise this is necessary so 
	                              // that a signal without electronic
	                              // noise is passed along
	    signal += electronics;
	    gi =iz*fNPixelsX+ix; // global index
	    if (signal > threshold) {
		digits[0]=iz;
		digits[1]=ix;
		digits[2]=1;
		for(j1=0;j1<3;j1++){
		    if (pList->GetTrack(ix, iz, gi)) {
			//b.b.	     tracks[j1]=-3;
			tracks[j1] = (Int_t)(pList->GetTrack(ix, iz, j1)+j1);
			hits[j1] = (Int_t)(pList->GetHit(ix, iz, j1)+j1+6);
		    }else {
			tracks[j1]=-2; //noise
			hits[j1] = -1;
		    } // end if pList
		    charges[j1] = 0;
		} // end for j1

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

		phys = 0;

		Int_t i       = ix; // These variables are declared so to be
		Int_t j       = iz; // passed along to UpdateMapNoise and
		Float_t noise = electronics; // in that function
		UpdateMapNoise(j, i, ix, iz, fModule, sig, noise, pList);
		aliITS->AddSimDigit(0, phys, digits, tracks, hits, charges);
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
