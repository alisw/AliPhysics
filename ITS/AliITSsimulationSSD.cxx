#include <stdio.h>
#include <iostream.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TMath.h>

#include "AliITSmodule.h"
#include "AliITSMapA2.h"   
#include "AliITSsegmentationSSD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsimulationSSD.h"
//#include "AliITSdictSSD.h"
#include "AliITSdcsSSD.h"
#include "AliITS.h"
#include "AliRun.h"


ClassImp(AliITSsimulationSSD);
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Enrico Fragiacomo
// July 2000
//
// AliITSsimulationSSD is the simulation of SSDs.

//------------------------------------------------------------
AliITSsimulationSSD::AliITSsimulationSSD(AliITSsegmentation *seg,
                                         AliITSresponse *resp){
  // Constructor

    fSegmentation = seg;
    fResponse = resp;
  Float_t noise[2] = {0.,0.};
  fResponse->GetNoiseParam(noise[0],noise[1]);   // retrieves noise parameters
  cout<<"nois1,2 ="<<noise[0]<<","<<noise[1]<<endl;    
  fDCS = new AliITSdcsSSD(seg,resp); 

    fNstrips = fSegmentation->Npx();
    fPitch = fSegmentation->Dpx(0);
    cout<<" Dx,Dz ="<<fSegmentation->Dx()<<","<<fSegmentation->Dz()<<endl;
    cout<<"fNstrips="<<fNstrips<<" fPitch="<<fPitch<<endl;

    //fP = new TArrayF(fNstrips+1); 
    //fN = new TArrayF(fNstrips+1);

    fMapA2 = new AliITSMapA2(fSegmentation);
     
    //fTracksP = new AliITSdictSSD[fNstrips+1];
    //fTracksN = new AliITSdictSSD[fNstrips+1];
    
    fSteps  = 100;   // still hard-wired - set in SetDetParam and get it via  
                     // fDCS together with the others eventually    
}

//___________________________________________________________________________
AliITSsimulationSSD& AliITSsimulationSSD::operator=(AliITSsimulationSSD 
                                                                      &source){
// Operator =

    if(this==&source) return *this;

    this->fDCS = new AliITSdcsSSD(*(source.fDCS));
    //this->fN   = new TArrayF(*(source.fN));
    //this->fP   = new TArrayF(*(source.fP));
    this->fMapA2 = source.fMapA2;
    //this->fTracksP = new AliITSdictSSD(*(source.fTracksP));
    //this->fTracksN = new AliITSdictSSD(*(source.fTracksN));
    this->fNstrips = source.fNstrips;
    this->fPitch   = source.fPitch;
    this->fSteps   = source.fSteps;
    return *this;
}
//_____________________________________________________________
AliITSsimulationSSD::AliITSsimulationSSD(AliITSsimulationSSD &source){
  // copy constructor

   *this = source;
}
//____________________________________________________________________________
AliITSsimulationSSD::~AliITSsimulationSSD() {
  // destructor    
  //if(fP) delete fP;
  //if(fN) delete fN;
    delete fMapA2;
    //if(fTracksP) delete [] fTracksP;
    //if(fTracksN) delete [] fTracksN;
    delete fDCS;
} 

//_______________________________________________________________
void AliITSsimulationSSD::DigitiseModule(AliITSmodule *mod,Int_t module,
                                         Int_t dummy) {
  // Digitizes hits for one SSD module 

  TObjArray *hits = mod->GetHits();
  Int_t nhits = hits->GetEntriesFast();
  if (!nhits) return;
  //cout<<"!! module, nhits ="<<module<<","<<nhits<<endl; //b.b.
  
  Double_t x0=0.0, y0=0.0, z0=0.0;
  Double_t x1=0.0, y1=0.0, z1=0.0;
  Double_t de=0.0;
  Int_t maxNdigits = 2*fNstrips; 
  Float_t  **pList = new Float_t* [maxNdigits]; 
  memset(pList,0,sizeof(Float_t*)*maxNdigits);
  Int_t indexRange[4] = {0,0,0,0};
  static Bool_t first=kTRUE;
  Int_t lasttrack = -2;
  Int_t idtrack = -2;
    
  for(Int_t i=0; i<nhits; i++) {    
    // LineSegmentL returns 0 if the hit is entering
    // If hits is exiting returns positions of entering and exiting hits
    // Returns also energy loss
    if (mod->LineSegmentL(i, x0, x1, y0, y1, z0, z1, de, idtrack)) {
      //cout<<"!! befor HitToDigit: hit ="<<i<<endl; //b.b.
      HitToDigit(module, x0, y0, z0, x1, y1, z1, de, indexRange, first);
       
      if (lasttrack != idtrack || i==(nhits-1)) {
	GetList(idtrack,pList,indexRange);
	first=kTRUE;
      }
      lasttrack=idtrack;
    }
  }  // end loop over hits
  
  ApplyNoise();
  ApplyCoupling();

  ChargeToSignal(pList);

  fMapA2->ClearMap();
}

//---------------------------------------------------------------
void AliITSsimulationSSD::HitToDigit(Int_t module, Double_t x0, Double_t y0, 
				     Double_t z0, Double_t x1, Double_t y1, 
				     Double_t z1, Double_t de,
				     Int_t *indexRange, Bool_t first) {
  // Turns hits in SSD module into one or more digits.

  AliITSsegmentationSSD *seg = new AliITSsegmentationSSD();
  // AliITSresponseSSD *res = new AliITSresponseSSD();

  Float_t tang[2] = {0.0,0.0};
  seg->Angles(tang[0], tang[1]); // stereo<<  -> tan(stereo)~=stereo
  //fSegmentation->Angles(tang[0], tang[1]); // stereo<<  -> tan(stereo)~=stereo
  Double_t x, y, z;
  Double_t dex=0.0, dey=0.0, dez=0.0;
  Double_t pairs;
  Double_t  ionE = 3.62E-9;   // ionization energy of Si (GeV)
  //ionE = (Double_t)  res->GetIonE();     
    
  Double_t sigma[2] = {0.,0.};  // standard deviation of the diffusion gaussian

  Double_t D[2] = {11.,30.};                // diffusion constant {h,e} (cm**2/sec)
  //D[0] = (Double_t) res->GetDiffusionConstantP();
  //D[1] = (Double_t) res->GetDiffusionConstantN();

  Double_t tdrift[2] = {0.,0.}; // time of drift
  Double_t vdrift[2] = {0.86E6,2.28E6};          // drift velocity (cm/sec)   
  //vdrift[0] = (Double_t) res->GetDriftVelocityP();
  //vdrift[1] = (Double_t) res->GetDriftVelocityN();

  Double_t w;
  Double_t inf[2], sup[2], par0[2];                 
  // Steps in the module are determined "manually" (i.e. No Geant)
  // NumOfSteps divide path between entering and exiting hits in steps 
  Int_t numOfSteps = NumOfSteps(x1, y1, z1, dex, dey, dez);
  
  // Enery loss is equally distributed among steps
  de = de/numOfSteps;
  pairs = de/ionE;             // e-h pairs generated

  //cout<<"Dy ="<<seg->Dy()<<endl;
  //cout<<"numOfSteps ="<<numOfSteps<<endl;
  //cout<<"dex,dey,dez ="<<dex<<","<<dey<<","<<dez<<endl;
  //cout<<"y0,y1 ="<<y0<<","<<y1<<endl;
  for(Int_t j=0; j<numOfSteps; j++) {     // stepping
    //    cout<<"step number ="<<j<<endl;
    x = x0 + (j+0.5)*dex;
    y = y0 + (j+0.5)*dey;
    if ( y > (seg->Dy()/2 +10)*1.0E-4 ) {
      //if ( y > (seg->Dy()*1.0E-4/2) ) {
      //if ( y > (fSegmentation->Dy()*1.0E-4/2) ) {
      // check if particle is within the detector
      cout<<"AliITSsimulationSSD::HitToDigit: Warning: hit out of detector y0,y,dey,j ="<<y0<<","<<y<<","<<dey<<","<<j<<endl;
      return;
    };
    z = z0 + (j+0.5)*dez;

    // calculate drift time
    tdrift[0] = (y+(seg->Dy()*1.0E-4)/2) / vdrift[0]; // y is the minimum path
    tdrift[1] = ((seg->Dy()*1.0E-4)/2-y) / vdrift[1]; // y is the minimum path
    //tdrift[0] = (y+(fSegmentation->Dy()*1.0E-4)/2) / vdrift[0]; // y is the minimum path
    //tdrift[1] = ((fSegmentation->Dy()*1.0E-4)/2-y) / vdrift[1]; // y is the minimum path

    for(Int_t k=0; k<2; k++) {   // both sides    remember: 0=Pside 1=Nside

      tang[k]=TMath::Tan(tang[k]);

      // w is the coord. perpendicular to the strips
      if(k==0) {
	w = (x+(seg->Dx()*1.0E-4)/2) - (z+(seg->Dz()*1.0E-4)/2)*tang[k]; 
	//w = (x+(fSegmentation->Dx()*1.0E-4)/2) - (z+(fSegmentation->Dz()*1.0E-4)/2)*tang[k]; 
	//cout<<"k,x,z,w ="<<k<<","<<x<<","<<z<<","<<w<<endl;
      }
      else {
	w = (x+(seg->Dx()*1.0E-4)/2) + (z-(seg->Dz()*1.0E-4)/2)*tang[k]; 
	//w = (x+(fSegmentation->Dx()*1.0E-4)/2) + (z-(fSegmentation->Dz()*1.0E-4)/2)*tang[k]; 
	//cout<<"k,x,z,w ="<<k<<","<<x<<","<<z<<","<<w<<endl;
      }
      w = w / (fPitch*1.0E-4); // w is converted in units of pitch

      if((w<(-0.5)) || (w>(fNstrips-0.5))) {
	// this check rejects hits in regions not covered by strips
	// 0.5 takes into account boundaries 
	if(k==0) cout<<"AliITSsimulationSSD::HitToDigit: Warning: no strip in this region of P side"<<endl;
	else cout<<"AliITSsimulationSSD::HitToDigit: Warning: no strip in this region of N side"<<endl;
	return;
      }

      // sigma is the standard deviation of the diffusion gaussian

      if(tdrift[k]<0) return;

      sigma[k] = TMath::Sqrt(2*D[k]*tdrift[k]);
      sigma[k] = sigma[k] /(fPitch*1.0E-4);  //units of Pitch
      if(sigma[k]==0.0) { 	
	cout<<"AliITSsimulationSSD::DigitiseModule: Error: sigma=0"<<endl; 
	exit(0);
      }

      par0[k] = pairs;
      // we integrate the diffusion gaussian from -3sigma to 3sigma 
      inf[k] = w - 3*sigma[k];      // 3 sigma from the gaussian average  
      sup[k] = w + 3*sigma[k];      // 3 sigma from the gaussian average
      // IntegrateGaussian does the actual integration of diffusion gaussian
      IntegrateGaussian(k, par0[k], w, sigma[k], inf[k], sup[k], 
      			indexRange, first);
    }  // end for loop over side (0=Pside, 1=Nside)      
  } // end stepping
  delete seg;
}

//____________________________________________________________________

void AliITSsimulationSSD::ApplyNoise() {
  // Apply Noise.
  
  Float_t signal;
  Float_t noise[2] = {0.,0.};
  fResponse->GetNoiseParam(noise[0],noise[1]);   // retrieves noise parameters

  for(Int_t k=0;k<2;k++){                        // both sides (0=Pside, 1=Nside)
    for(Int_t ix=0;ix<fNstrips;ix++){            // loop over strips
      signal = (Float_t) fMapA2->GetSignal(k,ix);// retrieves signal from map

      signal += gRandom->Gaus(0,noise[k]);       // add noise to signal
      if(signal<0.) signal=0.0;                  // in case noise is negative...

      fMapA2->SetHit(k,ix,(Double_t)signal);     // give back signal to map
    } // loop over strip 
  } // loop over k (P or N side)
}

//_________________________________________________________________________

void AliITSsimulationSSD::ApplyCoupling() {
  // Apply the effect of electronic coupling between channels    
  Float_t signal, signalLeft=0, signalRight=0;

  for(Int_t ix=0;ix<fNstrips;ix++){
    if(ix>0.) signalLeft  = (Float_t) fMapA2->GetSignal(0,ix-1)*fDCS->GetCouplingPL();
    else signalLeft = 0.0;
    if(ix<(fNstrips-1)) signalRight = (Float_t) fMapA2->GetSignal(0,ix+1)*fDCS->GetCouplingPR();
    else signalRight = 0.0;
    signal = (Float_t) fMapA2->GetSignal(0,ix);
    signal += signalLeft + signalRight;
    fMapA2->SetHit(0,ix,(Double_t)signal);
    
    if(ix>0.) signalLeft  = (Float_t) fMapA2->GetSignal(1,ix-1)*fDCS->GetCouplingNL();
    else signalLeft = 0.0;
    if(ix<(fNstrips-1)) signalRight = (Float_t) fMapA2->GetSignal(1,ix+1)*fDCS->GetCouplingNR();
    else signalRight = 0.0;
    signal = (Float_t) fMapA2->GetSignal(1,ix);
    signal += signalLeft + signalRight;
    fMapA2->SetHit(1,ix,(Double_t)signal);
  } // loop over strips 
}

//____________________________________________________________________________
Float_t AliITSsimulationSSD::F(Float_t av, Float_t x, Float_t s) {
  // Computes the integral of a gaussian using Error Function
  Float_t sqrt2 = TMath::Sqrt(2.0);
  Float_t sigm2 = sqrt2*s;
  Float_t integral;

  integral = 0.5 * TMath::Erf( (x - av) / sigm2);
  return integral;
} 

//_________________________________________________________________________
void AliITSsimulationSSD::IntegrateGaussian(Int_t k,Double_t par, Double_t w,
					    Double_t sigma, 
					    Double_t inf, Double_t sup,
					    Int_t *indexRange, Bool_t first) {
  // integrate the diffusion gaussian
  // remind: inf and sup are w-3sigma and w+3sigma
  //         we could define them here instead of passing them
  //         this way we are free to introduce asimmetry

  Double_t a=0.0, b=0.0;
  Double_t signal = 0.0, dXCharge1 = 0.0, dXCharge2 = 0.0;
  // dXCharge1 and 2 are the charge to two neighbouring strips
  // Watch that we only involve at least two strips
  // Numbers greater than 2 of strips in a cluster depend on
  //  geometry of the track and delta rays, not charge diffusion!   
  
  Double_t strip = TMath::Floor(w);         // clostest strip on the left

  if ( TMath::Abs((strip - w)) < 0.5) { 
    // gaussian mean is closer to strip on the left
    a = inf;                         // integration starting point
    if((strip+0.5)<=sup) {
      // this means that the tail of the gaussian goes beyond
      // the middle point between strips ---> part of the signal
      // is given to the strip on the right
      b = strip + 0.5;               // integration stopping point
      dXCharge1 = F( w, b, sigma) - F(w, a, sigma);    
      dXCharge2 = F( w, sup, sigma) - F(w ,b, sigma); 
    }
    else { 
      // this means that all the charge is given to the strip on the left
      b = sup;
      dXCharge1 = 0.9973;   // gaussian integral at 3 sigmas
      dXCharge2 = 0.0;
    }

    dXCharge1 = par * dXCharge1;  // normalize by mean of number of carriers
    dXCharge2 = par * dXCharge2;

    // for the time being, signal is the charge
    // in ChargeToSignal signal is converted in ADC channel
    signal = fMapA2->GetSignal(k,strip);
    signal += dXCharge1;

    fMapA2->SetHit(k,strip,(Double_t)signal);
    if(((Int_t) strip) < (fNstrips-1)) {
      // strip doesn't have to be the last (remind: last=fNstrips-1)
      // otherwise part of the charge is lost
      signal = fMapA2->GetSignal(k,(strip+1));
      signal += dXCharge2;
      fMapA2->SetHit(k,(strip+1),(Double_t)signal);
    }
    
    if(dXCharge1 > 1.) {
      if (first) {
	indexRange[k*2+0]=indexRange[k*2+1]=(Int_t) strip;
	first=kFALSE;
      }

      indexRange[k*2+0]=TMath::Min(indexRange[k*2+0],(Int_t) strip);
      indexRange[k*2+1]=TMath::Max(indexRange[k*2+1],(Int_t) strip);
    }      // dXCharge > 1 e-

  }
  else {
    // gaussian mean is closer to strip on the right
    strip++;     // move to strip on the rigth
    b = sup;     // now you know where to stop integrating
    if((strip-0.5)>=inf) { 
      // tail of diffusion gaussian on the left goes left of
      // middle point between strips
      a = strip - 0.5;        // integration starting point
      dXCharge1 = F(w, b, sigma) - F(w, a, sigma);
      dXCharge2 = F(w, a, sigma) - F(w, inf, sigma);
    }
    else {
      a = inf;
      dXCharge1 = 0.9973;   // gaussian integral at 3 sigmas
      dXCharge2 = 0.0;
    }
    
    dXCharge1 = par * dXCharge1;    // normalize by means of carriers
    dXCharge2 = par * dXCharge2;

    // for the time being, signal is the charge
    // in ChargeToSignal signal is converted in ADC channel
    signal = fMapA2->GetSignal(k,strip);
    signal += dXCharge1;
    fMapA2->SetHit(k,strip,(Double_t)signal);
    if(((Int_t) strip) > 0) {
      // strip doesn't have to be the first
      // otherwise part of the charge is lost
      signal = fMapA2->GetSignal(k,(strip-1));
      signal += dXCharge2;
      fMapA2->SetHit(k,(strip-1),(Double_t)signal);
    }
    
    if(dXCharge1 > 1.) {
      if (first) {
	indexRange[k*2+0]=indexRange[k*2+1]=(Int_t) strip;
	first=kFALSE;
      }

      indexRange[k*2+0]=TMath::Min(indexRange[k*2+0],(Int_t) strip);
      indexRange[k*2+1]=TMath::Max(indexRange[k*2+1],(Int_t) strip);
    }      // dXCharge > 1 e-
  }
}


  //_________________________________________________________________________

Int_t AliITSsimulationSSD::NumOfSteps(Double_t x, Double_t y, Double_t z,
				 Double_t & dex,Double_t & dey,Double_t & dez) {  
  // number of steps
  // it also returns steps for each coord
  //AliITSsegmentationSSD *seg = new AliITSsegmentationSSD();

  Double_t step = 25E-4;
  //step = (Double_t) seg->GetStepSize();  // step size (cm)
  Int_t numOfSteps = (Int_t) (TMath::Sqrt(x*x+y*y+z*z)/step); 

  if (numOfSteps < 1) numOfSteps = 1;       // one step, at least

  // we could condition the stepping depending on the incident angle
  // of the track
  dex = x/numOfSteps;
  dey = y/numOfSteps;
  dez = z/numOfSteps;

  return numOfSteps;

}

//---------------------------------------------
void AliITSsimulationSSD::GetList(Int_t label,Float_t **pList,Int_t *indexRange) {
  // loop over nonzero digits
  Int_t ix,globalIndex;
  Float_t signal=0.;
  Float_t highest,middle,lowest;
  // printf("SPD-GetList: indexRange[0] indexRange[1] indexRange[2] indexRange[3] %d %d %d %d\n",indexRange[0], indexRange[1], indexRange[2], indexRange[3]);

  for(Int_t k=0; k<2; k++) {
    for(ix=indexRange[k*2+0];ix<indexRange[k*2+1]+1;ix++){
      if(indexRange[k*2+0]<indexRange[k*2+1]) 
	signal=fMapA2->GetSignal(k,ix);
      
      globalIndex = k*fNstrips+ix; // globalIndex starts from 0!
      if(!pList[globalIndex]){
	
	// 
	// Create new list (6 elements - 3 signals and 3 tracks + total sig)
	//
	
	pList[globalIndex] = new Float_t [6];
	// set list to -1 
	
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
void AliITSsimulationSSD::ChargeToSignal(Float_t **pList) {
  // charge to signal  

  AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");
  
  Float_t threshold = 0.;

  Int_t digits[3], tracks[3],hits[3],gi,j1;
  Float_t charges[3];
  Float_t signal,phys;
  Float_t noise[2] = {0.,0.};
  fResponse->GetNoiseParam(noise[0],noise[1]);
  
  for(Int_t k=0;k<2;k++){         // both sides (0=Pside, 1=Nside)

    // Threshold for zero-suppression
    // It can be defined in AliITSresponseSSD
    //             threshold = (Float_t)fResponse->MinVal(k);
    // I prefer to think adjusting the threshold "manually", looking
    // at the scope, and considering noise standard deviation
    threshold = 4.0*noise[k];      // 4 times noise is a choice
    //cout<<"SSD: k,thresh ="<<k<<","<<threshold<<endl;
    for(Int_t ix=0;ix<fNstrips;ix++){         // loop over strips

      signal = (Float_t) fMapA2->GetSignal(k,ix);

      gi =k*fNstrips+ix; // global index
      if (signal > threshold) {
	digits[0]=k;
	digits[1]=ix;

	// convert to ADC signal
	// conversion factor are rather arbitrary (need tuning)
	// minimum ionizing particle --> ~30000 pairs --> ADC channel 50
	signal = signal*50.0/30000.0;        
	//cout<<"SSD: 1 signal ="<<signal<<endl;
	if(signal>1000.) signal = 1000.0; // if exceeding, accumulate last one
	//cout<<"SSD: 2 signal ="<<signal<<endl;
	digits[2]=(Int_t) signal;

	//gi =k*fNstrips+ix; // global index
	for(j1=0;j1<3;j1++){
	  if (pList[gi]) {
	    tracks[j1] = (Int_t)(*(pList[gi]+j1));
	  }	  
	  else {
	    tracks[j1]=-2; //noise
	  }
	  charges[j1] = 0;
	}

	phys=0;

	hits[0]=0;
	hits[1]=0;
	hits[2]=0;
	aliITS->AddSimDigit(2,phys,digits,tracks,hits,charges);  // finally add digit

	//if(pList[gi]) delete [] pList[gi];
      }
      if(pList[gi]) delete [] pList[gi];
    }
  }
  delete [] pList;
}














