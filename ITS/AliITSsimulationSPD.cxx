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
#include "AliITSsimulationSPD.h"
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"




ClassImp(AliITSsimulationSPD)
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Boris Batyunya
// December 20 1999
//
// AliITSsimulationSPD is the simulation of SPDs
//________________________________________________________________________


AliITSsimulationSPD::AliITSsimulationSPD()
{
  // constructor
  fResponse = 0;
  fSegmentation = 0;
  fMapA2=0;
  fHis = 0;
  fNoise=0.;
  fBaseline=0.;
  fNPixelsZ=0;
  fNPixelsX=0;
}


//_____________________________________________________________________________

AliITSsimulationSPD::AliITSsimulationSPD(AliITSsegmentation *seg, AliITSresponse *resp) {
  // standard constructor

      fHis = 0;
      fResponse = resp;
      fSegmentation = seg;

      fResponse->GetNoiseParam(fNoise,fBaseline);

      fMapA2 = new AliITSMapA2(fSegmentation);

      //
      fNPixelsZ=fSegmentation->Npz();
      fNPixelsX=fSegmentation->Npx();

}

//_____________________________________________________________________________

AliITSsimulationSPD::~AliITSsimulationSPD() { 
  // destructor

  delete fMapA2;

  if (fHis) {
     fHis->Delete(); 
     delete fHis;     
  }                
}


//__________________________________________________________________________
AliITSsimulationSPD::AliITSsimulationSPD(const AliITSsimulationSPD &source){
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

//_________________________________________________________________________
AliITSsimulationSPD& 
  AliITSsimulationSPD::operator=(const AliITSsimulationSPD &source) {
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
//_____________________________________________________________________________

void AliITSsimulationSPD::DigitiseModule(AliITSmodule *mod, Int_t module, Int_t dummy)
{
  // digitize module

    const Float_t kEnToEl = 2.778e+8; // GeV->charge in electrons 
                                      // for 3.6 eV/pair 
    const Float_t kconv = 10000.;     // cm -> microns

    Float_t spdLength = fSegmentation->Dz();
    Float_t spdWidth = fSegmentation->Dx();

    Float_t difCoef, dum; 
    fResponse->DiffCoeff(difCoef,dum); 

    Float_t zPix0 = 1e+6;
    Float_t xPix0 = 1e+6;
    Float_t yPix0 = 1e+6;
    Float_t yPrev = 1e+6;   
    Float_t zP0 = 100.;
    Float_t xP0 = 100.;

    Float_t zPitch = fSegmentation->Dpz(0);
    Float_t xPitch = fSegmentation->Dpx(0);
  
    //cout << "pitch per z: " << zPitch << endl;
    //cout << "pitch per r*phi: " << xPitch << endl;

    TObjArray *fHits = mod->GetHits();
    Int_t nhits = fHits->GetEntriesFast();
    if (!nhits) return;
    //    cout << "module, nhits ="<<module<<","<<nhits<< endl;

  //  Array of pointers to the label-signal list

    Int_t maxNDigits = fNPixelsX*fNPixelsZ + fNPixelsX ;; 
    Float_t  **pList = new Float_t* [maxNDigits]; 
    memset(pList,0,sizeof(Float_t*)*maxNDigits);
    Int_t indexRange[4] = {0,0,0,0};

    // Fill detector maps with GEANT hits
    // loop over hits in the module
    static Bool_t first;
    Int_t lasttrack=-2,idhit=-1;
    Int_t hit, iZi, jz, jx;
    for (hit=0;hit<nhits;hit++) {
        AliITShit *iHit = (AliITShit*) fHits->At(hit);
	Int_t layer = iHit->GetLayer();

	// work with the idtrack=entry number in the TreeH
	//Int_t idhit,idtrack;
	//mod->GetHitTrackAndHitIndex(hit,idtrack,idhit);    
	//Int_t idtrack=mod->GetHitTrackIndex(hit);  
        // or store straight away the particle position in the array
	// of particles : 

        if(iHit->StatusEntering()) idhit=hit;
        Int_t itrack = iHit->GetTrack();
        Int_t dray = 0;
   
	if (lasttrack != itrack || hit==(nhits-1)) first = kTRUE; 

	//        Int_t parent = iHit->GetParticle()->GetFirstMother();
        Int_t partcode = iHit->GetParticle()->GetPdgCode();

//  partcode (pdgCode): 11 - e-, 13 - mu-, 22 - gamma, 111 - pi0, 211 - pi+
//                      310 - K0s, 321 - K+, 2112 - n, 2212 - p, 3122 - lambda

        Float_t px = iHit->GetPXL();
        Float_t py = iHit->GetPYL();
        Float_t pz = iHit->GetPZL();
        Float_t pmod = 1000*sqrt(px*px+py*py+pz*pz);


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
	if(TMath::Abs(zPix) > spdLength/2.) {
	  printf("!! Zpix outside = %f\n",zPix);
	  if(status == 66) zP0=100;
	  continue;
	} 


	if (TMath::Abs(xPix) > spdWidth/2.) {
	  printf("!! Xpix outside = %f\n",xPix);
	  if (status == 66) xP0=100;
	  continue;
	}  

	Float_t zP = (zPix + spdLength/2.)/1000.;  
	Float_t xP = (xPix + spdWidth/2.)/1000.;  

	Int_t trdown = 0;

	// enter Si or after event in Si
	if (status == 66 ) {  
           zPix0 = zPix;
           xPix0 = xPix;
           yPrev = yPix; 
	}   
	// enter Si only
	if (layer == 1 && status == 66 && yPix > 71.) {     
             yPix0 = yPix;
             zP0 = zP;
             xP0 = xP;
	}
	// enter Si only
	if (layer == 2 && status == 66 && yPix < -71.) {    
             yPix0 = yPix;
             zP0 = zP;
             xP0 = xP;
	}       
	Float_t depEnergy = iHit->GetIonization();
	// skip if the input point to Si       
	if(depEnergy <= 0.) continue;        
	// skip if the input point is outside of Si, but the next
	// point is inside of Si
	if(zP0 > 90 || xP0 > 90) continue; 
	// if track returns to the opposite direction:
	if (layer == 1 && yPix > yPrev) {
	    yPix0 = yPrev;
            trdown = 1;
	}
	if (layer == 2 && yPix < yPrev) {
            yPix0 = yPrev;
            trdown = 1;
	} 

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
	Float_t ydif = yPix - yPix0;

	if(TMath::Abs(ydif) < 0.1) continue; // Ydif is not zero

	Float_t projDif = sqrt(xdif*xdif + zdif*zdif);
	Int_t ndZ = (Int_t)TMath::Abs(zdif/zPitch) + 1;
	Int_t ndX = (Int_t)TMath::Abs(xdif/xPitch) + 1; 

	// number of the steps along the track:
	Int_t nsteps = ndZ;
	if(ndX > ndZ) nsteps = ndX;
	if(nsteps < 6) nsteps = 6;  // minimum number of the steps 

	if(TMath::Abs(projDif) > 5.0) tang = ydif/projDif;
	Float_t dCharge = charge/nsteps;       // charge in e- for one step
	Float_t dZ = zdif/nsteps;
	Float_t dX = xdif/nsteps;

	if (TMath::Abs(projDif) < 5.0 ) {
	   drPath = ydif*1.e-4;  
           drPath = TMath::Abs(drPath);        // drift path in cm
	   sigmaDif = difCoef*sqrt(drPath);    // sigma diffusion in cm        
	}  

	for (iZi = 1;iZi <= nsteps;iZi++) {
            Float_t dZn = iZi*dZ;
	    Float_t dXn = iZi*dX;
	    Float_t zPixn = zPix0 + dZn;
	    Float_t xPixn = xPix0 + dXn;

	    if(TMath::Abs(projDif) >= 5.) {
	      Float_t dProjn = sqrt(dZn*dZn+dXn*dXn);
	      if(trdown == 0) {
		drPath = dProjn*tang*1.e-4; // drift path for iZi step in cm 
		drPath = TMath::Abs(drPath);
	      }
	      if(trdown == 1) {
		Float_t dProjn = projDif/nsteps; 
		drPath = (projDif-(iZi-1)*dProjn)*tang*1.e-4;
		drPath = TMath::Abs(drPath);
	      }
	      sigmaDif = difCoef*sqrt(drPath);    
	      sigmaDif = sigmaDif*kconv;         // sigma diffusion in microns
	    }
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
		} 
		if(jz == 2) {
		  dZprev = -dZleft;
		  dZnext = dZright;
		} 
		if(jz == 3) {
		  dZprev = dZright;
		  dZnext = dZright + zPitch;
		} 
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
		     } 
		     if(jx == 2) {
		       dXprev = -dXleft;
		       dXnext = dXright;
		     } 
		     if(jx == 3) {
		       dXprev = dXright;
		       dXnext = dXright + xPitch;
		     } 
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
		       }

                       indexRange[0]=TMath::Min(indexRange[0],kz-1);
                       indexRange[1]=TMath::Max(indexRange[1],kz-1);
                       indexRange[2]=TMath::Min(indexRange[2],kx-1);
                       indexRange[3]=TMath::Max(indexRange[3],kx-1);

		       // build the list of digits for this module	
                       Double_t signal=fMapA2->GetSignal(index,kx-1);
                       signal+=dXCharge;
                       fMapA2->SetHit(index,kx-1,(double)signal);
		     }      // dXCharge > 1 e-
		  }       // jx loop
		}       // dZCharge > 1 e-
	    }        // jz loop
	}         // iZi loop

        if (status == 65) {   // the step is inside of Si
	   zPix0 = zPix;
	   xPix0 = xPix;
        }
	yPrev = yPix;  

	if(dray == 0) {
            GetList(itrack,idhit,pList,indexRange);
	}

	lasttrack=itrack;
    }   // hit loop inside the module

   
    // introduce the electronics effects and do zero-suppression
    ChargeToSignal(pList); 

    // clean memory

    fMapA2->ClearMap();


} 

//---------------------------------------------
void AliITSsimulationSPD::GetList(Int_t label,Int_t idhit,Float_t **pList,Int_t *indexRange)
{
  // lop over nonzero digits

   
  //set protection
  for(int k=0;k<4;k++) {
     if (indexRange[k] < 0) indexRange[k]=0;
  }

  for(Int_t iz=indexRange[0];iz<indexRange[1]+1;iz++){
    for(Int_t ix=indexRange[2];ix<indexRange[3]+1;ix++){

        Float_t signal=fMapA2->GetSignal(iz,ix);

	if (!signal) continue;

        Int_t globalIndex = iz*fNPixelsX+ix; // GlobalIndex starts from 0!
        if(!pList[globalIndex]){

           // 
	   // Create new list (9 elements - 3 signals and 3 tracks + 3 hits)
	   //

           pList[globalIndex] = new Float_t [9];

	   // set list to -3 

	   *pList[globalIndex] = -3.;
	   *(pList[globalIndex]+1) = -3.;
	   *(pList[globalIndex]+2) = -3.;
	   *(pList[globalIndex]+3) =  0.;
	   *(pList[globalIndex]+4) =  0.;
	   *(pList[globalIndex]+5) =  0.;
	   *(pList[globalIndex]+6) = -1.;
	   *(pList[globalIndex]+7) = -1.;
	   *(pList[globalIndex]+8) = -1.;


	   *pList[globalIndex] = (float)label;
	   *(pList[globalIndex]+3) = signal;
	   *(pList[globalIndex]+6) = (float)idhit;
        }
        else{

	  // check the signal magnitude

          Float_t highest = *(pList[globalIndex]+3);
          Float_t middle = *(pList[globalIndex]+4);
          Float_t lowest = *(pList[globalIndex]+5);

          signal -= (highest+middle+lowest);

	  //
	  //  compare the new signal with already existing list
	  //

          if(signal<lowest) continue; // neglect this track

          if (signal>highest){
            *(pList[globalIndex]+5) = middle;
            *(pList[globalIndex]+4) = highest;
            *(pList[globalIndex]+3) = signal;

            *(pList[globalIndex]+2) = *(pList[globalIndex]+1);
            *(pList[globalIndex]+1) = *pList[globalIndex];
            *pList[globalIndex] = label;

            *(pList[globalIndex]+8) = *(pList[globalIndex]+7);
            *(pList[globalIndex]+7) = *(pList[globalIndex]+6);
            *(pList[globalIndex]+6) = idhit;
	  }
          else if (signal>middle){
            *(pList[globalIndex]+5) = middle;
            *(pList[globalIndex]+4) = signal;

            *(pList[globalIndex]+2) = *(pList[globalIndex]+1);
            *(pList[globalIndex]+1) = label;

            *(pList[globalIndex]+8) = *(pList[globalIndex]+7);
            *(pList[globalIndex]+7) = idhit;
	  }
          else{
            *(pList[globalIndex]+5) = signal;
            *(pList[globalIndex]+2) = label;
            *(pList[globalIndex]+8) = idhit;
	  }
        }
    } // end of loop pixels in x
  } // end of loop over pixels in z


}


//---------------------------------------------
void AliITSsimulationSPD::ChargeToSignal(Float_t **pList)
{
  // add noise and electronics, perform the zero suppression and add the
  // digit to the list

  AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");
  

  TRandom *random = new TRandom(); 
  Float_t threshold = (float)fResponse->MinVal();

  Int_t digits[3], tracks[3], hits[3],gi,j1;
  Float_t charges[3];
  Float_t electronics;
  Float_t signal,phys;
  for(Int_t iz=0;iz<fNPixelsZ;iz++){
    for(Int_t ix=0;ix<fNPixelsX;ix++){
      electronics = fBaseline + fNoise*random->Gaus();
      signal = (float)fMapA2->GetSignal(iz,ix);
      signal += electronics;
      gi =iz*fNPixelsX+ix; // global index
      if (signal > threshold) {
	 digits[0]=iz;
	 digits[1]=ix;
	 digits[2]=1;
	 for(j1=0;j1<3;j1++){
	   if (pList[gi]) {
	     tracks[j1]=-3;
	     tracks[j1] = (Int_t)(*(pList[gi]+j1));
	     hits[j1] = (Int_t)(*(pList[gi]+j1+6));
	   }else {
	     tracks[j1]=-2; //noise
	     hits[j1] = -1;
	   }
	   charges[j1] = 0;
	 }

	 if(tracks[0] == tracks[1] && tracks[0] == tracks[2]) {
	   tracks[1] = -3;
           hits[1] = -1;
	   tracks[2] = -3;
           hits[2] = -1;
         } 
	 if(tracks[0] == tracks[1] && tracks[0] != tracks[2]) {
	   tracks[1] = -3;
           hits[1] = -1;   
         } 
	 if(tracks[0] == tracks[2] && tracks[0] != tracks[1]) {
	   tracks[2] = -3;
           hits[2] = -1;   
         } 
	 if(tracks[1] == tracks[2] && tracks[0] != tracks[1]) {
	   tracks[2] = -3;
           hits[2] = -1;   
         } 

         phys=0;
	 aliITS->AddSimDigit(0,phys,digits,tracks,hits,charges);
      }
      if(pList[gi]) delete [] pList[gi];
    }
  }
  delete [] pList;

}


//____________________________________________

void AliITSsimulationSPD::CreateHistograms()
{
  // create 1D histograms for tests

      printf("SPD - create histograms\n");

      fHis=new TObjArray(fNPixelsZ);
      TString spdName("spd_");
      for (Int_t i=0;i<fNPixelsZ;i++) {
	   Char_t pixelz[4];
	   sprintf(pixelz,"%d",i+1);
	   spdName.Append(pixelz);
	   (*fHis)[i] = new TH1F(spdName.Data(),"SPD maps",
                              fNPixelsX,0.,(Float_t) fNPixelsX);
      }
}

//____________________________________________

void AliITSsimulationSPD::ResetHistograms()
{
    //
    // Reset histograms for this detector
    //
    for ( int i=0;i<fNPixelsZ;i++ ) {
	if ((*fHis)[i])    ((TH1F*)(*fHis)[i])->Reset();
    }

}









