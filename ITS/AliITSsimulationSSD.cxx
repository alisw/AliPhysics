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

#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TRandom.h>

#include "AliITSmodule.h"
#include "AliITSMapA2.h"
#include "AliITSpList.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSSD.h"
//#include "AliITSdcsSSD.h"
#include "AliITS.h"
#include "AliITShit.h"
#include "AliITSdigitSSD.h"
#include "AliRun.h"
#include "AliITSgeom.h"
#include "AliITSsimulationSSD.h"
#include "AliITSTableSSD.h"
//#include "AliITSresponseSSD.h"

ClassImp(AliITSsimulationSSD)
////////////////////////////////////////////////////////////////////////
//                                                                    //
// Author: Enrico Fragiacomo                                          //
//         enrico.fragiacomo@ts.infn.it                               //
// Last revised: march 2006                                           // 
//                                                                    //
// AliITSsimulationSSD is the simulation of SSD.                     //
////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------
AliITSsimulationSSD::AliITSsimulationSSD():AliITSsimulation(),
					   //fDCS(0),
fMapA2(0),
fIonE(0.0),
fDifConst(),
fDriftVel(){
    //default Constructor
    //Inputs:
    // none.
    // Outputs:
    // none.
    // Return:
    //  A default construction AliITSsimulationSSD class
}
//----------------------------------------------------------------------
AliITSsimulationSSD::AliITSsimulationSSD(AliITSDetTypeSim* dettyp):
AliITSsimulation(dettyp),
//fDCS(0),
fMapA2(0),
fIonE(0.0),
fDifConst(),
fDriftVel(){
    // Constructor 
    // Input:
    //   AliITSDetTypeSim    Pointer to the SSD dettype to be used
    // Outputs:
    //   none.
    // Return
    //   A standard constructed AliITSsimulationSSD class

    Init();
}
//----------------------------------------------------------------------
void AliITSsimulationSSD::Init(){
  // Inilizer, Inilizes all of the variable as needed in a standard place.
  // Input:
  //   AliITSsegmentationSSD *seg  Pointer to the SSD segmentation to be used
  //   AliITSCalibrationSSD   *resp Pointer to the SSD responce class to be used
  // Outputs:
  //   none.
  // Return
  //   none.
  AliITSsegmentationSSD* seg = (AliITSsegmentationSSD*)GetSegmentationModel(2);
  
  SetDriftVelocity(); // use default values in .h file
  SetIonizeE();       // use default values in .h file
  SetDiffConst();     // use default values in .h file
  fpList           = new AliITSpList(2,GetNStrips());
  fMapA2           = new AliITSMapA2(seg);
}
//______________________________________________________________________
AliITSsimulationSSD& AliITSsimulationSSD::operator=(
                                         const AliITSsimulationSSD &s){
  // Operator =

  if(this==&s) return *this;

  //  this->fDCS         = new AliITSdcsSSD(*(s.fDCS));
  this->fMapA2       = s.fMapA2;
  this->fIonE        = s.fIonE;
  this->fDifConst[0] = s.fDifConst[0];
  this->fDifConst[1] = s.fDifConst[1];
  this->fDriftVel[0] = s.fDriftVel[0];
  this->fDriftVel[1] = s.fDriftVel[1];
  return *this;
}
//______________________________________________________________________
AliITSsimulation& AliITSsimulationSSD::operator=(
                                         const AliITSsimulation &s){
  // Operator =

  if(this==&s) return *this;
  Error("AliITSsimulationSSD","Not allowed to make a = with "
	"AliITSsimulationSSD Using default creater instead");
  
  return *this;
}
//______________________________________________________________________
AliITSsimulationSSD::AliITSsimulationSSD(const AliITSsimulationSSD &source):
    AliITSsimulation(source),
fMapA2(source.fMapA2),
fIonE(source.fIonE),
fDifConst(),
fDriftVel(){
  // copy constructor
  fDifConst[0] = source.fDifConst[0];
  fDifConst[1] = source.fDifConst[1];
  fDriftVel[0] = source.fDriftVel[0];
  fDriftVel[1] = source.fDriftVel[1];
}
//______________________________________________________________________
AliITSsimulationSSD::~AliITSsimulationSSD() {
  // destructor
  delete fMapA2;
  //delete fDCS;
}
//______________________________________________________________________
void AliITSsimulationSSD::InitSimulationModule(Int_t module,Int_t event){
    // Creates maps to build the list of tracks for each sumable digit
    // Inputs:
    //   Int_t module    // Module number to be simulated
    //   Int_t event     // Event number to be simulated
    // Outputs:
    //   none.
    // Return
    //    none.

    SetModuleNumber(module);
    SetEventNumber(event);
    fMapA2->ClearMap();
    fpList->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSSD::FinishSDigitiseModule(){
    // Does the Sdigits to Digits work
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

  FillMapFrompList(fpList);  // need to check if needed here or not????
  SDigitToDigit(fModule,fpList);
  fpList->ClearMap();
  fMapA2->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSSD::DigitiseModule(AliITSmodule *mod,Int_t,Int_t) {
  // Digitizes hits for one SSD module
  SetModuleNumber(mod->GetIndex());
  
  HitsToAnalogDigits(mod,fpList);
  SDigitToDigit(GetModuleNumber(),fpList);
  
  fpList->ClearMap();
  fMapA2->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSSD::SDigitiseModule(AliITSmodule *mod,Int_t,Int_t) {
  // Produces Summable/Analog digits and writes them to the SDigit tree. 

    HitsToAnalogDigits(mod,fpList);

    WriteSDigits(fpList);
    
    fpList->ClearMap();
    fMapA2->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSSD::SDigitToDigit(Int_t module,AliITSpList *pList){
  // Takes the pList and finishes the digitization.
  
  ApplyNoise(pList,module);
  ApplyCoupling(pList,module);
  ApplyDeadChannels(module);
  
  ChargeToSignal(module,pList);
}
//______________________________________________________________________
void AliITSsimulationSSD::HitsToAnalogDigits(AliITSmodule *mod,
                                             AliITSpList *pList){
    // Loops over all hits to produce Analog/floating point digits. This
    // is also the first task in producing standard digits.
  Int_t lasttrack     = -2;
  Int_t idtrack       = -2;
  Double_t x0=0.0, y0=0.0, z0=0.0;
  Double_t x1=0.0, y1=0.0, z1=0.0;
  Double_t de=0.0;
  Int_t module = mod->GetIndex();
  
  AliITSsegmentationSSD* seg = (AliITSsegmentationSSD*)GetSegmentationModel(2);
  
    TObjArray *hits = mod->GetHits();
    Int_t nhits     = hits->GetEntriesFast();
    if (nhits<=0) return;
    AliITSTableSSD * tav = new AliITSTableSSD(GetNStrips());
    module = mod->GetIndex();
    if ( mod->GetLayer() == 6 ) seg->SetLayer(6);
    if ( mod->GetLayer() == 5 ) seg->SetLayer(5);
    for(Int_t i=0; i<nhits; i++) {    
      // LineSegmentL returns 0 if the hit is entering
      // If hits is exiting returns positions of entering and exiting hits
      // Returns also energy loss
      if(GetDebug(4)){
	cout << i << " ";
	cout << mod->GetHit(i)->GetXL() << " "<<mod->GetHit(i)->GetYL();
	cout << " " << mod->GetHit(i)->GetZL();
	cout << endl;
      } // end if
      if (mod->LineSegmentL(i, x0, x1, y0, y1, z0, z1, de, idtrack)) {
	HitToDigit(module, x0, y0, z0, x1, y1, z1, de,tav);
	if (lasttrack != idtrack || i==(nhits-1)) {
	  GetList(idtrack,i,module,pList,tav);
	} // end if
	lasttrack=idtrack;
      } // end if
    }  // end loop over hits
    delete tav; tav=0;
    return;
}
//----------------------------------------------------------------------
void AliITSsimulationSSD::HitToDigit(Int_t module, Double_t x0, Double_t y0, 
                                     Double_t z0, Double_t x1, Double_t y1, 
                                     Double_t z1, Double_t de,
                                     AliITSTableSSD *tav) {
  
  // hit to digit conversion
  
  AliITSsegmentationSSD* seg = (AliITSsegmentationSSD*)GetSegmentationModel(2);
  // Turns hits in SSD module into one or more digits.
  Float_t tang[2] = {0.0,0.0};
  seg->Angles(tang[0], tang[1]);//stereo<<->tan(stereo)~=stereo
  Double_t x, y, z;
  Double_t dex=0.0, dey=0.0, dez=0.0; 
  Double_t pairs; // pair generation energy per step.
  Double_t sigma[2] = {0.,0.};// standard deviation of the diffusion gaussian
  Double_t tdrift[2] = {0.,0.}; // time of drift
  Double_t w;
  Double_t inf[2], sup[2], par0[2];                 
  
  // Steps in the module are determined "manually" (i.e. No Geant)
  // NumOfSteps divide path between entering and exiting hits in steps 
  Int_t numOfSteps = NumOfSteps(x1, y1, z1, dex, dey, dez);
  // Enery loss is equally distributed among steps
  de    = de/numOfSteps;
  pairs = de/GetIonizeE(); // e-h pairs generated
  for(Int_t j=0; j<numOfSteps; j++) {     // stepping
    x = x0 + (j+0.5)*dex;
    y = y0 + (j+0.5)*dey;
    if ( y > (seg->Dy()/2+10)*1.0E-4 ) {
      // check if particle is within the detector
      Warning("HitToDigit",
	      "hit out of detector y0=%e,y=%e,dey=%e,j =%d module=%d,  exceed=%e",
	      y0,y,dey,j,module, y-(seg->Dy()/2+10)*1.0E-4);
      return;
    } // end if
    z = z0 + (j+0.5)*dez;
    if(GetDebug(4)) cout <<"HitToDigit "<<x<<" "<<y<<" "<<z<< " "
			 <<dex<<" "<<dey<<" "<<dez<<endl;
    // calculate drift time
    // y is the minimum path
    tdrift[0] = (y+(seg->Dy()*1.0E-4)/2)/GetDriftVelocity(0);
    tdrift[1] = ((seg->Dy()*1.0E-4)/2-y)/GetDriftVelocity(1);
    
    for(Int_t k=0; k<2; k++) {   // both sides    remember: 0=Pside 1=Nside
      
      tang[k]=TMath::Tan(tang[k]);
      
      // w is the coord. perpendicular to the strips
      Float_t xp=x*1.e+4,zp=z*1.e+4; // microns
      seg->GetPadTxz(xp,zp);
      if(k==0) w = xp; // P side strip number
      else w = zp; // N side strip number
      
      if((w<(-0.5)) || (w>(GetNStrips()-0.5))) {
	// this check rejects hits in regions not covered by strips
	// 0.5 takes into account boundaries 
	if(GetDebug(4)) cout << "x,z="<<x<<","<<z<<" w="<<w
			     <<" Nstrips="<<GetNStrips()<<endl;
	return; // There are dead region on the SSD sensitive volume.
      } // end if
      
      // sigma is the standard deviation of the diffusion gaussian
      if(tdrift[k]<0) return;
      sigma[k] = TMath::Sqrt(2*GetDiffConst(k)*tdrift[k]);
      sigma[k] /= (GetStripPitch()*1.0E-4);  //units of Pitch
      if(sigma[k]==0.0) { 	
	Error("HitToDigit"," sigma[%d]=0",k);
	exit(0);
      } // end if
      
      par0[k] = pairs;
      // we integrate the diffusion gaussian from -3sigma to 3sigma 
      inf[k] = w - 3*sigma[k]; // 3 sigma from the gaussian average  
      sup[k] = w + 3*sigma[k]; // 3 sigma from the gaussian average
      // IntegrateGaussian does the actual
      // integration of diffusion gaussian
      IntegrateGaussian(k, par0[k], w, sigma[k], inf[k], sup[k],tav);
    }  // end for loop over side (0=Pside, 1=Nside)      
  } // end stepping
}

//______________________________________________________________________
void AliITSsimulationSSD::ApplyNoise(AliITSpList *pList,Int_t module){
  // Apply Noise.
  Int_t ix;
  Double_t signal,noise;
  AliITSCalibrationSSD* res =(AliITSCalibrationSSD*)GetCalibrationModel(module);
  
  // Pside
  for(ix=0;ix<GetNStrips();ix++){      // loop over strips
    
    // noise is gaussian
    noise  = (Double_t) gRandom->Gaus(0,res->GetNoiseP().At(ix));
    
    // need to calibrate noise 
    // NOTE. noise from the calibration database comes uncalibrated, 
    // it needs to be calibrated in order to be added
    // to the signal. It will be decalibrated later on together with the noise    
    noise *= (Double_t) res->GetGainP(ix); 
    
    // noise comes in ADC channels from the calibration database
    // It needs to be converted back to electronVolts
    noise /= res->GetDEvToADC(1.);
    
    // Finally, noise is added to the signal
    signal = noise + fMapA2->GetSignal(0,ix);//get signal from map
    fMapA2->SetHit(0,ix,signal); // give back signal to map
    if(signal>0.0) pList->AddNoise(0,ix,module,noise);
  } // loop over strip 
  
    // Nside
  for(ix=0;ix<GetNStrips();ix++){      // loop over strips
    noise  = (Double_t) gRandom->Gaus(0,res->GetNoiseN().At(ix));// give noise to signal
    noise *= (Double_t) res->GetGainN(ix); 
    noise /= res->GetDEvToADC(1.);
    signal = noise + fMapA2->GetSignal(1,ix);//get signal from map
    fMapA2->SetHit(1,ix,signal); // give back signal to map
    if(signal>0.0) pList->AddNoise(1,ix,module,noise);
  } // loop over strip 
  
}
//______________________________________________________________________
void AliITSsimulationSSD::ApplyCoupling(AliITSpList *pList,Int_t module) {
  // Apply the effect of electronic coupling between channels
  Int_t ix;
  Double_t signal=0;
  AliITSCalibrationSSD* res =(AliITSCalibrationSSD*)GetCalibrationModel(module);
  
  Double_t *contrLeft  = new Double_t[GetNStrips()];
  Double_t *contrRight = new Double_t[GetNStrips()];
  
  // P side coupling
  for(ix=0;ix<GetNStrips();ix++){
    if(ix>0) contrLeft[ix] = fMapA2->GetSignal(0,ix-1)*res->GetCouplingPL();
    else contrLeft[ix] = 0.0;
    if(ix<(GetNStrips()-1)) contrRight[ix] = fMapA2->GetSignal(0,ix+1)*res->GetCouplingPR();
    else contrRight[ix] = 0.0;
  } // loop over strips 
  
  for(ix=0;ix<GetNStrips();ix++){
    signal = contrLeft[ix] + contrRight[ix] - res->GetCouplingPL() * fMapA2->GetSignal(0,ix)
      - res->GetCouplingPR() * fMapA2->GetSignal(0,ix);
    fMapA2->AddSignal(0,ix,signal);
    if(signal>0.0) pList->AddNoise(0,ix,module,signal);
  } // loop over strips 
  
  // N side coupling
  for(ix=0;ix<GetNStrips();ix++){
    if(ix>0) contrLeft[ix] = fMapA2->GetSignal(1,ix-1)*res->GetCouplingNL();
    else contrLeft[ix] = 0.0;
    if(ix<(GetNStrips()-1)) contrRight[ix] = fMapA2->GetSignal(1,ix+1)*res->GetCouplingNR();
    else contrRight[ix] = 0.0;
  } // loop over strips 
  
  for(ix=0;ix<GetNStrips();ix++){
    signal = contrLeft[ix] + contrRight[ix] - res->GetCouplingNL() * fMapA2->GetSignal(0,ix)
      - res->GetCouplingNR() * fMapA2->GetSignal(0,ix);
    fMapA2->AddSignal(1,ix,signal);
    if(signal>0.0) pList->AddNoise(1,ix,module,signal);
  } // loop over strips 
  

  delete [] contrLeft;
  delete [] contrRight; 
}

//______________________________________________________________________
void AliITSsimulationSSD::ApplyDeadChannels(Int_t module) {
  // Kill dead channels setting gain to zero

  Int_t deadentries;

  AliITSCalibrationSSD* res = (AliITSCalibrationSSD*)GetCalibrationModel(module);

  deadentries = res->GetDeadPChannelsList().GetSize();
  //cout<<module<<" "<<deadentries<<endl;
  for(Int_t i=0; i<deadentries; i++) {
    res->AddGainP(res->GetDeadPChannelsList().At(i),0.0);
  }

  deadentries = res->GetDeadNChannelsList().GetSize();
  for(Int_t i=0; i<deadentries; i++) {
    res->AddGainN(res->GetDeadNChannelsList().At(i),0.0);
  }

}

//______________________________________________________________________
Float_t AliITSsimulationSSD::F(Float_t av, Float_t x, Float_t s) {
    // Computes the integral of a gaussian using Error Function
    Float_t sqrt2 = TMath::Sqrt(2.0);
    Float_t sigm2 = sqrt2*s;
    Float_t integral;

    integral = 0.5 * TMath::Erf( (x - av) / sigm2);
    return integral;
}
//______________________________________________________________________
void AliITSsimulationSSD::IntegrateGaussian(Int_t k,Double_t par, Double_t w,
                                            Double_t sigma, 
                                            Double_t inf, Double_t sup,
                                            AliITSTableSSD *tav) {
    // integrate the diffusion gaussian
    // remind: inf and sup are w-3sigma and w+3sigma
    //         we could define them here instead of passing them
    //         this way we are free to introduce asimmetry

    Double_t a=0.0, b=0.0;
    Double_t dXCharge1 = 0.0, dXCharge2 = 0.0;
    // dXCharge1 and 2 are the charge to two neighbouring strips
    // Watch that we only involve at least two strips
    // Numbers greater than 2 of strips in a cluster depend on
    //  geometry of the track and delta rays, not charge diffusion!   

    Double_t strip = TMath::Floor(w);         // closest strip on the left

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
        }else { 
            // this means that all the charge is given to the strip on the left
            b = sup;
            dXCharge1 = 0.9973;   // gaussian integral at 3 sigmas
            dXCharge2 = 0.0;
        } // end if
        dXCharge1 = par * dXCharge1;// normalize by mean of number of carriers
        dXCharge2 = par * dXCharge2;

        // for the time being, signal is the charge
        // in ChargeToSignal signal is converted in ADC channel
        fMapA2->AddSignal(k,(Int_t)strip,dXCharge1);
        tav->Add(k,(Int_t)strip);
        if(((Int_t) strip) < (GetNStrips()-1)) {
            // strip doesn't have to be the last (remind: last=GetNStrips()-1)
            // otherwise part of the charge is lost
            fMapA2->AddSignal(k,((Int_t)strip+1),dXCharge2);
            tav->Add(k,((Int_t)(strip+1)));
        } // end if
    }else{
        // gaussian mean is closer to strip on the right
        strip++;     // move to strip on the rigth
        b = sup;     // now you know where to stop integrating
        if((strip-0.5)>=inf) { 
            // tail of diffusion gaussian on the left goes left of
            // middle point between strips
            a = strip - 0.5;        // integration starting point
            dXCharge1 = F(w, b, sigma) - F(w, a, sigma);
            dXCharge2 = F(w, a, sigma) - F(w, inf, sigma);
        }else {
            a = inf;
            dXCharge1 = 0.9973;   // gaussian integral at 3 sigmas
            dXCharge2 = 0.0;
        } // end if
        dXCharge1 = par * dXCharge1;    // normalize by means of carriers
        dXCharge2 = par * dXCharge2;
        // for the time being, signal is the charge
        // in ChargeToSignal signal is converted in ADC channel
        fMapA2->AddSignal(k,(Int_t)strip,dXCharge1);
        tav->Add(k,(Int_t)strip);
        if(((Int_t) strip) > 0) {
            // strip doesn't have to be the first
            // otherwise part of the charge is lost
            fMapA2->AddSignal(k,((Int_t)strip-1),dXCharge2);
            tav->Add(k,((Int_t)(strip-1)));
        } // end if
    } // end if
}
//______________________________________________________________________
Int_t AliITSsimulationSSD::NumOfSteps(Double_t x, Double_t y, Double_t z,
                                      Double_t &dex,Double_t &dey,
                                      Double_t &dez){
    // number of steps
    // it also returns steps for each coord
    //AliITSsegmentationSSD *seg = new AliITSsegmentationSSD();

    Double_t step = 25E-4;
    //step = (Double_t) seg->GetStepSize();  // step size (cm)
    Int_t numOfSteps = (Int_t) (TMath::Sqrt(x*x+y*y+z*z)/step); 

    if (numOfSteps < 1) numOfSteps = 1;       // one step, at least
    //numOfSteps=1;

    // we could condition the stepping depending on the incident angle
    // of the track
    dex = x/numOfSteps;
    dey = y/numOfSteps;
    dez = z/numOfSteps;
    
    return numOfSteps;
}
//----------------------------------------------------------------------
void AliITSsimulationSSD::GetList(Int_t label,Int_t hit,Int_t mod,
                                  AliITSpList *pList,AliITSTableSSD *tav) {
    // loop over nonzero digits
    Int_t ix,i;
    Double_t signal=0.;

    for(Int_t k=0; k<2; k++) {
        ix=tav->Use(k);
        while(ix>-1){
            signal = fMapA2->GetSignal(k,ix);
            if(signal==0.0) {
                ix=tav->Use(k);
                continue;
            } // end if signal==0.0
            // check the signal magnitude
            for(i=0;i<pList->GetNSignals(k,ix);i++){
                signal -= pList->GetTSignal(k,ix,i);
            } // end for i
            //  compare the new signal with already existing list
            if(signal>0)pList->AddSignal(k,ix,label,hit,mod,signal);
            ix=tav->Use(k);
        } // end of loop on strips
    } // end of loop on P/N side
    tav->Clear();
}
//----------------------------------------------------------------------
void AliITSsimulationSSD::ChargeToSignal(Int_t module,AliITSpList *pList) {
    // charge to signal
    static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");
    Float_t threshold = 0.;
    Int_t size = AliITSdigitSSD::GetNTracks();
    Int_t * digits = new Int_t[size];
    Int_t * tracks = new Int_t[size];
    Int_t * hits = new Int_t[size];
    Int_t j1;
    Float_t charges[3] = {0.0,0.0,0.0};
    Float_t signal;
    AliITSCalibrationSSD* res =(AliITSCalibrationSSD*)GetCalibrationModel(module);

    for(Int_t k=0;k<2;k++){         // both sides (0=Pside, 1=Nside)
      for(Int_t ix=0;ix<GetNStrips();ix++){     // loop over strips

	// if strip is dead -> gain=0
	if( ((k==0)&&(res->GetGainP(ix)==0)) || ((k==1)&&(res->GetGainN(ix)==0))) continue;
	
	// signal has to be uncalibrated
	// In real life, gains are supposed to be calculated from calibration runs,
	// stored in the calibration DB and used in the reconstruction
	// (see AliITSClusterFinderSSD.cxx)
	if(k==0) signal /= res->GetGainP(ix);
	else signal /= res->GetGainN(ix);

	// signal is converted in unit of ADC
	signal = res->GetDEvToADC(fMapA2->GetSignal(k,ix));
	if(signal>4096.) signal = 4096.;//if exceeding, accumulate last one

	// threshold for zero suppression is set on the basis of the noise
	// A good value is 3*sigma_noise
	if(k==0) threshold = res->GetNoiseP().At(ix);
	else threshold = res->GetNoiseN().At(ix);
	threshold *= res->GetZSThreshold(); // threshold at 3 sigma noise
	if(signal < threshold) continue;

	digits[0] = k;
	digits[1] = ix;
	digits[2] = TMath::Nint(signal);
	for(j1=0;j1<size;j1++)if(j1<pList->GetNEntries()){
	  // only three in digit.
	  tracks[j1]  = pList->GetTrack(k,ix,j1);
	  hits[j1]    = pList->GetHit(k,ix,j1);
	}else{
	  tracks[j1]  = -3;
	  hits[j1]    = -1;
	} // end for j1
	// finally add digit
	aliITS->AddSimDigit(2,0,digits,tracks,hits,charges);
      } // end for ix
    } // end for k
    delete [] digits;
    delete [] tracks;
    delete [] hits;
}
//______________________________________________________________________
void AliITSsimulationSSD::WriteSDigits(AliITSpList *pList){
    // Fills the Summable digits Tree
    Int_t i,ni,j,nj;
    static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");

    pList->GetMaxMapIndex(ni,nj);
    for(i=0;i<ni;i++)for(j=0;j<nj;j++){
        if(pList->GetSignalOnly(i,j)>0.0){
            aliITS->AddSumDigit(*(pList->GetpListItem(i,j)));
            if(GetDebug(4)) cout << "pListSSD: "<<*(pList->GetpListItem(i,j))
                                << endl;
        } // end if
    } // end for i,j
  return;
}
//______________________________________________________________________
void AliITSsimulationSSD::FillMapFrompList(AliITSpList *pList){
    // Fills fMap2A from the pList of Summable digits
    Int_t k,ix;

    for(k=0;k<2;k++)for(ix=0;ix<GetNStrips();ix++) 
        fMapA2->AddSignal(k,ix,pList->GetSignal(k,ix));
    return;
}
//______________________________________________________________________
void AliITSsimulationSSD::Print(ostream *os){
    //Standard output format for this class

    //AliITSsimulation::Print(os);
    *os << fIonE <<",";
    *os << fDifConst[0] <<","<< fDifConst[1] <<",";
    *os << fDriftVel[0] <<","<< fDriftVel[1];
    //*os <<","; fDCS->Print(os);
    //*os <<","; fMapA2->Print(os);
}
//______________________________________________________________________
void AliITSsimulationSSD::Read(istream *is){
    // Standard output streaming function.

    //AliITSsimulation::Read(is);
    *is >> fIonE;
    *is >> fDifConst[0] >> fDifConst[1];
    *is >> fDriftVel[0] >> fDriftVel[1];
    //fDCS->Read(is);
    //fMapA2->Read(is);
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSsimulationSSD &source){
    // Standard output streaming function.

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &os,AliITSsimulationSSD &source){
    // Standard output streaming function.

    source.Read(&os);
    return os;
}
//______________________________________________________________________





