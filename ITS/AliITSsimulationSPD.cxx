#include <iostream.h>
#include <TRandom.h>
#include <TH1.h>

#include "AliRun.h"
#include "AliITSMap.h"    // "AliITSMapA2.h"
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
//
//________________________________________________________________________

AliITSsimulationSPD::AliITSsimulationSPD(){
  // constructor
  fResponse     = 0;
  fSegmentation = 0;
  fHis          = 0;
  fNoise        = 0.;
  fBaseline     = 0.;
}
//_____________________________________________________________________________

AliITSsimulationSPD::AliITSsimulationSPD(AliITSsegmentation *seg, AliITSresponse *resp) {
  // constructor
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

void AliITSsimulationSPD::DigitiseModule(AliITSmodule *mod, Int_t module, Int_t dummy) {
  // digitize module
    const Float_t kEntoEl = 2.778e+8; // GeV->charge in electrons 
                                      // for 3.6 eV/pair 
    const Float_t kconv = 10000.;     // cm -> microns

    Float_t spdLenght = fSegmentation->Dz();
    Float_t spdWidth = fSegmentation->Dx();

    Float_t difCoef = fResponse->DiffCoeff();  

    Float_t zPix0 = 1e+6;
    Float_t xPix0 = 1e+6;
    Float_t yPix0 = 1e+6;
    Float_t yPrev = 1e+6;

    Float_t zPitch = fSegmentation->Dpz(0);
    Float_t xPitch = fSegmentation->Dpx(0);
  
    //cout << "pitch per z: " << zPitch << endl;
    //cout << "pitch per r*phi: " << xPitch << endl;

    TObjArray *fHits = mod->GetHits();
    Int_t nhits = fHits->GetEntriesFast();
    if (!nhits) return;



  //  Array of pointers to the label-signal list

    Int_t maxNdigits = fNPixelsX*fNPixelsZ; 
    Float_t  **pList = new Float_t* [maxNdigits]; 
    memset(pList,0,sizeof(Float_t*)*maxNdigits);
    Int_t indexRange[4] = {0,0,0,0};

    // Fill detector maps with GEANT hits
    // loop over hits in the module

    Int_t layer,idtrack,status,trdown,ndZ,ndX,nsteps,iZi,nZpix,nXpix;
    Int_t jzmin,jzmax,jxmin,jxmax,jx,jz,lz,lx,index;
    Float_t zArg1,zArg2,xArg1,xArg2;
    Float_t zProb1,zProb2,xProb1,xProb2;
    Float_t dZCharge,dXCharge;
    Double_t signal;
    Float_t projDif;  // RMS of xDif and zDif
    Float_t dProjn;  // RMS of dXn and dZn
    Float_t depEnergy; // The deposited energy from this hit
    Float_t zPix; // hit position in microns
    Float_t xPix; // hit position in microns
    Float_t yPix; // hit position in microns
    Float_t zPixn; // hit position in microns
    Float_t xPixn; // hit position in microns
    Float_t charge,dCharge;// charge in e-
    Float_t drPath;   
    Float_t tAng,dZ,dX,dXn,dZn;
    Float_t sigmaDif;
    Float_t dZright,dZleft;
    Float_t dXright,dXleft;
    Float_t dZprev,dZnext,dXprev,dXnext;

    static Bool_t first=kTRUE;
    Int_t lasttrack = -2,hit;
    for(hit=0;hit<nhits;hit++) {
        AliITShit *iHit = (AliITShit*) fHits->At(hit);
	layer = iHit->GetLayer();

	// work with the idtrack=entry number in the TreeH
	// Int_t idtrack=mod->GetHitTrackIndex(ii);  
        // or store straight away the particle position in the array
	// of particles : 
        idtrack = iHit->GetTrack();

	//if(module==11 || module==157) {
	  // Int_t track = iHit->fTrack;
	  // Int_t primary = gAlice->GetPrimary(track);
	  // Int_t parent = iHit->GetParticle()->GetFirstMother();
	  // printf("module,hit,track,primary,parent  %d %d %d %d %d \n",
	  //                                   md,hit,track,primary,parent);
	//}
   
	//  Get hit z and x(r*phi) cordinates for each module (detector)
	//  in local system.

	zPix = kconv*iHit->GetZL(); // Geant cm to microns
	xPix = kconv*iHit->GetXL(); // Geant cm to micron
        yPix = kconv*iHit->GetYL(); // Geant cm to micron

	// Get track status
	status = iHit->GetTrackStatus();      
	if(status != 66 && status != 68) {
	  printf("!!!!! no order status  %d\n",status);
	}

	trdown = 0;

	// enter Si or after event in Si
	if (status == 66 ) {  
           zPix0 = zPix;
           xPix0 = xPix;
           yPrev = yPix; 
	}   
	// enter Si only
	if (layer == 1 && status == 66 && yPix > 71.) {     
             yPix0 = yPix;
	}
	// enter Si only
	if (layer == 2 && status == 66 && yPix < -71.) {    
             yPix0 = yPix;
	}
	depEnergy = iHit->GetIonization();
	// skip if the input point to Si       
	if(depEnergy <= 0.) continue;  
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

	charge = depEnergy*kEntoEl;         // charge in e-
	drPath = 0.;   
	tAng = 0.;
	sigmaDif = 0.; 
	Float_t zDif = zPix - zPix0;
	Float_t xDif = xPix - xPix0;
	Float_t yDif = yPix - yPix0;

	if(TMath::Abs(yDif) < 0.1) continue; // yDif is not zero


	projDif = sqrt(xDif*xDif + zDif*zDif);
	ndZ = (Int_t)TMath::Abs(zDif/zPitch) + 1;
	ndX = (Int_t)TMath::Abs(xDif/xPitch) + 1; 

	// number of the steps along the track:
	nsteps = ndZ;
	if(ndX > ndZ) nsteps = ndX;
	if(nsteps < 6) nsteps = 6;  // minimum number of the steps 

	if(TMath::Abs(projDif) > 5.0) tAng = yDif/projDif;
	dCharge = charge/nsteps;       // charge in e- for one step
	dZ = zDif/nsteps;
	dX = xDif/nsteps;

	if (TMath::Abs(projDif) < 5.0 ) {
	   drPath = yDif*1.e-4;  
           drPath = TMath::Abs(drPath);        // drift path in cm
	   sigmaDif = difCoef*sqrt(drPath);    // sigma diffusion in cm        
	}  

	for(iZi = 1;iZi <= nsteps;iZi++) {
            dZn = iZi*dZ;
	    dXn = iZi*dX;
	    zPixn = zPix0 + dZn;
	    xPixn = xPix0 + dXn;

	    if(TMath::Abs(projDif) >= 5.) {
	      dProjn = sqrt(dZn*dZn+dXn*dXn);
	      if(trdown == 0) {
		drPath = dProjn*tAng*1.e-4; // drift path for iZi step in cm 
		drPath = TMath::Abs(drPath);
	      }
	      if(trdown == 1) {
		dProjn = projDif/nsteps; 
		drPath = (projDif-(iZi-1)*dProjn)*tAng*1.e-4;
		drPath = TMath::Abs(drPath);
	      }
	      sigmaDif = difCoef*sqrt(drPath);    
	      sigmaDif = sigmaDif*kconv;         // sigma diffusion in microns
	    }
	    zPixn = (zPixn + spdLenght/2.);  
	    xPixn = (xPixn + spdWidth/2.);
            fSegmentation->GetCellIxz(xPixn,zPixn,nXpix,nZpix);
	    zPitch = fSegmentation->Dpz(nZpix);
	    // set the window for the integration
	    jzmin = 1;
	    jzmax = 3; 
	    if(nZpix == 1) jzmin =2;
	    if(nZpix == fNPixelsZ) jzmax = 2; 

	    jxmin = 1;  
	    jxmax = 3; 
	    if(nXpix == 1) jxmin =2;
	    if(nXpix == fNPixelsX) jxmax = 2; 

	    zPix    = nZpix;
	    dZright = zPitch*(zPix - zPixn);
	    dZleft  = zPitch - dZright;
	    xPix    = nXpix;
	    dXright = xPitch*(xPix - xPixn);
	    dXleft  = xPitch - dXright;
	    dZprev  = 0.;
	    dZnext  = 0.;
	    dXprev  = 0.;
	    dXnext  = 0.;

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
		// lz changes from 1 to the fNofPixels(270)
		lz = nZpix + jz -2; 

		zArg1 = dZprev/sigmaDif;
		zArg2 = dZnext/sigmaDif;
		zProb1 = TMath::Erfc(zArg1);
		zProb2 = TMath::Erfc(zArg2);
		dZCharge =0.5*(zProb1-zProb2)*dCharge; 

		//printf("dZCharge %f \n",dZCharge);

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
		     lx = nXpix + jx -2;  

		     xArg1    = dXprev/sigmaDif;
		     xArg2    = dXnext/sigmaDif;
		     xProb1   = TMath::Erfc(xArg1);
		     xProb2   = TMath::Erfc(xArg2);
		     dXCharge =0.5*(xProb1-xProb2)*dZCharge; 

		     //printf("dXCharge %f \n",dXCharge);

		     if(dXCharge > 1.) {
		       index = lz-1;
		       if (first) {
                          indexRange[0]=indexRange[1]=index;
                          indexRange[2]=indexRange[3]=lx-1;
                          first=kFALSE;
		       }

                       indexRange[0]=TMath::Min(indexRange[0],lz-1);
                       indexRange[1]=TMath::Max(indexRange[1],lz-1);
                       indexRange[2]=TMath::Min(indexRange[2],lx-1);
                       indexRange[3]=TMath::Max(indexRange[3],lx-1);
		       // build the list of digits for this module	
                       signal=fMapA2->GetSignal(index,lx-1);
                       signal+=dXCharge;
                       fMapA2->SetHit(index,lx-1,(double)signal);
		     }      // dXCharge > 1 e-
		  }       // jx loop
		}       // dZCharge > 1 e-
	    }        // jz loop
	}         // iZi loop

        if (status == 65) {   // the step is inside of Si
	   zPix0 = zPix;
	   xPix0 = xPix;
        }
	yPrev = yPix;  //ch

	if (lasttrack != idtrack || hit==(nhits-1)) {
            GetList(idtrack,pList,indexRange);
            first=kTRUE;
	}
	lasttrack=idtrack;
    }   // hit loop inside the module

   
    // introduce the electronics effects and do zero-suppression
    ChargeToSignal(pList); 

    // clean memory

    fMapA2->ClearMap();


} 

//---------------------------------------------
void AliITSsimulationSPD::GetList(Int_t label,Float_t **pList,Int_t *indexRange) {
  // loop over nonzero digits

  Int_t ix,iz,globalIndex;
  Float_t signal;
  Float_t highest,middle,lowest;

  for(iz=indexRange[0];iz<indexRange[1]+1;iz++){
    for(ix=indexRange[2];ix<indexRange[3]+1;ix++){
      
      signal=fMapA2->GetSignal(iz,ix);

        globalIndex = iz*fNPixelsX+ix; // globalIndex starts from 0!
        if(!pList[globalIndex]){
        
           // 
	   // Create new list (6 elements - 3 signals and 3 tracks + total sig)
	   //

           pList[globalIndex] = new Float_t [6];

	   // set list to -2 

	   *pList[globalIndex] = -2.;
	   *(pList[globalIndex]+1) = -2.;
	   *(pList[globalIndex]+2) = -2.;
	   *(pList[globalIndex]+3) =  0.;
	   *(pList[globalIndex]+4) =  0.;
	   *(pList[globalIndex]+5) =  0.;


	   *pList[globalIndex] = (float)label;
	   *(pList[globalIndex]+3) = signal;
        }
        else{

	  // check the signal magnitude

          highest = *(pList[globalIndex]+3);
          middle = *(pList[globalIndex]+4);
          lowest = *(pList[globalIndex]+5);

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
	  }
          else if (signal>middle){
            *(pList[globalIndex]+5) = middle;
            *(pList[globalIndex]+4) = signal;

            *(pList[globalIndex]+2) = *(pList[globalIndex]+1);
            *(pList[globalIndex]+1) = label;
	  }
          else{
            *(pList[globalIndex]+5) = signal;
            *(pList[globalIndex]+2) = label;
	  }
        }
    } // end of loop pixels in x
  } // end of loop over pixels in z


}


//---------------------------------------------
void AliITSsimulationSPD::ChargeToSignal(Float_t **pList) {
  // charge to signal  

  AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");
  

  TRandom *random = new TRandom(); 
  Float_t threshold = (float)fResponse->MinVal();

  Int_t digits[3], tracks[3],gi,j1;
  Float_t charges[3];
  Float_t electronics;
  Float_t signal,phys;
  Int_t iz,ix;
  for(iz=0;iz<fNPixelsZ;iz++){
    for(ix=0;ix<fNPixelsX;ix++){
      electronics = fBaseline + fNoise*random->Gaus();
      signal = (float)fMapA2->GetSignal(iz,ix);
      signal += electronics;
      if (signal > threshold) {
	 digits[0]=iz;
	 digits[1]=ix;
	 digits[2]=1;
	 gi =iz*fNPixelsX+ix; // global index
	 for(j1=0;j1<3;j1++){
	   tracks[j1] = (Int_t)(*(pList[gi]+j1));
	   charges[j1] = 0;
	 }
         phys=0;
	 aliITS->AddDigit(0,phys,digits,tracks,charges);
         if(pList[gi]) delete [] pList[gi];
      }
    }
  }
  delete [] pList;


}


//____________________________________________

void AliITSsimulationSPD::CreateHistograms() {
  // CreateHistograms

      Int_t i;
      for(i=0;i<fNPixelsZ;i++) {
	   TString *spdname = new TString("spd_");
	   Char_t candnum[4];
	   sprintf(candnum,"%d",i+1);
	   spdname->Append(candnum);
	   (*fHis)[i] = new TH1F(spdname->Data(),"SPD maps",
                              fNPixelsX,0.,(Float_t) fNPixelsX);
	   delete spdname;
      }

}

//____________________________________________

void AliITSsimulationSPD::ResetHistograms() {
    //
    // Reset histograms for this detector
    //
    Int_t i;
    for(i=0;i<fNPixelsZ;i++ ) {
	if ((*fHis)[i])    ((TH1F*)(*fHis)[i])->Reset();
    }

}
