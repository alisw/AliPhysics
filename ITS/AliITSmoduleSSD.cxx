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

/*
$Log$
Revision 1.2  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/


#include "AliITSmoduleSSD.h"

    //Technical parameters of detector
const Float_t   AliITSmoduleSSD::kStereo = 0.0175;  //Stereo Angle 17.5 mrad
const Float_t   AliITSmoduleSSD::kTan = 0.0175;  
const Int_t     AliITSmoduleSSD::kNStrips = 768;    //Number of strips on each side
const Float_t   AliITSmoduleSSD::kPitch = 0.095;    //Distance strip - strip (mm)
const Float_t   AliITSmoduleSSD::kX = 72.96;        //X size (mm)
const Float_t   AliITSmoduleSSD::kY = 0.3;          //Y size (mm)
const Float_t   AliITSmoduleSSD::kZ = 40;           //Thickness (mm)
    
    // <------------------------------
  
    //______________________________________________________________
    //  
    // Parameters for simulation
    //______________________________________________________________
      
const Float_t   AliITSmoduleSSD::kSigmaP = 0.003;     //Gaussian sigm
const Float_t   AliITSmoduleSSD::kSigmaN = 0.002;
const Int_t     AliITSmoduleSSD::kSteps  = 10;        //Number of steps 
const Int_t     AliITSmoduleSSD::kTresholdP = 1500;    
const Int_t     AliITSmoduleSSD::kTresholdN = 2500; 
   


ClassImp(AliITSmoduleSSD)

//____________________________________________________________________
//
//  Constructor
//____________________________________________________________________
//                                     


AliITSmoduleSSD::AliITSmoduleSSD() {
    

    //Invalid Strips parameters
    
    fInvalidP = new TArrayS(0);
    fInvalidN = new TArrayS(0);
    
    fNInvalidP = 0;
    fNInvalidN = 0;
    
    
    //DCS parameters
    
    fGainP = 83;
    fGainN = 83;   
    fSNRatioP = 600;
    fSNRatioP = 500;
    
    fCouplingPR = 0.021;
    fCouplingPL = 0.026;
    fCouplingNR = 0.013;
    fCouplingNL = 0.010;

}

AliITSmoduleSSD::AliITSmoduleSSD(Int_t index) {
                                
    fIndex = index;	
    
    //Invalid Strips parameters
    
    fInvalidP = new TArrayS(0);
    fInvalidN = new TArrayS(0);
    
    fNInvalidP = 0;
    fNInvalidN = 0;
    
    
    //DCS parameters
    
    fGainP = 83;
    fGainN = 83;   
    fSNRatioP = 600;
    fSNRatioP = 500;
    
    fCouplingPR = 0.021;
    fCouplingPL = 0.026;
    fCouplingNR = 0.013;
    fCouplingNL = 0.010;	
}

AliITSmoduleSSD::~AliITSmoduleSSD() {
   
    if (!fInvalidP) delete fInvalidP;
    if (!fInvalidN) delete fInvalidN;				
}


//____________________________________________________________________
//
//  Inalid strips menagement
//____________________________________________________________________
//                           



void AliITSmoduleSSD::SetInvalidP(Int_t strip, Bool_t b) {
    
    Bool_t already = kFALSE;
    Int_t i;
    
    for (i=0; i<fNInvalidP; i++) {
        if ((*fInvalidP)[i] == strip) {
           already = kTRUE;
           break;
        }
    }
    
    if (!already && b) {
       (*fInvalidP)[fNInvalidP++] = strip;
    } 
}

void AliITSmoduleSSD::SetInvalidMC(Float_t m, Float_t s) {

    fNInvalid = m;
    fISigma = s;
    SetInvalidMC();
}

void AliITSmoduleSSD::SetInvalidMC() {

    Int_t m = (Int_t)gRandom->Gaus(fNInvalid, fISigma);
    
    for(int i=0; i<m; i++) {
       SetInvalidP((Int_t)(gRandom->Rndm()*kNStrips), kTRUE);
    }
    
}

Bool_t AliITSmoduleSSD::IsValidP(Int_t n) {
    
    for(Int_t i=0; i<fNInvalidP; i++) 
       if ((*fInvalidP)[i] == n) return kFALSE;
    return kTRUE;    
}

Bool_t AliITSmoduleSSD::IsValidN(Int_t n) {
    
    for(Int_t i=0; i<fNInvalidN; i++) 
       if ((*fInvalidN)[i] == n) return kFALSE;
    return kTRUE;    
}


//____________________________________________________________________
//
//  Add digit
//____________________________________________________________________
// 

/*********************************************************************
* 
* AddDigits 
* sets paramerers: layer, ladder detector
* scan tracks wich produced this digit
* creates new SSD DIGTS
* call ITS to add digit to its Array
* set index frm ITS in its own array
*
* S.Radomski 17.09.1999
*
*********************************************************************/

void AliITSmoduleSSD::AddDigit(Int_t strNo, Int_t s, Bool_t p) {
 
    Int_t tracks[3];
    Int_t digits[4];
    AliITSdigit *t = (AliITSdigit*) (new AliITSdigitSSD(tracks, digits, 
							strNo, s, p));
    
    fIdigits->AddAt(((AliITS *)fITS)->AddDigit(t), fNdigits++);
}


//____________________________________________________________________
//
//  Hits to digits
//____________________________________________________________________
//                          



void AliITSmoduleSSD::HitToDigit() {

    Int_t i;                           //for iteration
    fP = new TArrayI(768);
    fN = new TArrayI(768); 
    
    fPtrack1 = new TArrayI(768);
    fPtrack2 = new TArrayI(768);
    fPtrack3 = new TArrayI(768);
    
    fNtrack1 = new TArrayI(768);
    fNtrack2 = new TArrayI(768);
    fNtrack3 = new TArrayI(768);
    
    for(i=0; i<kNStrips; i++) {
       (*fN)[i] = 0;
       (*fP)[i] = 0;
    } // end for i
     
    for(i=0; i<fNhitsM; i++) HitToDigit(i);

    ApplyCoupling();    
    ApplyNoise();
        
    for(i=0; i<fNInvalidP; i++) (*fP)[(*fInvalidP)[i]] = -20;
    for(i=0; i<fNInvalidN; i++) (*fN)[(*fInvalidN)[i]] = -20;
    
    for(i=0; i<kNStrips; i++) 
       if ((*fP)[i]>kTresholdP) AddDigit(i+1, (*fP)[i], kTRUE);
    
    for(i=0; i<kNStrips; i++)
       if ((*fN)[i]>kTresholdN) AddDigit(i+1, (*fN)[i], kFALSE);
       
    delete fN;
    delete fP;
    
    delete fPtrack1;
    delete fPtrack2;
    delete fPtrack3;
    
    delete fPtrack1;
    delete fPtrack2;
    delete fPtrack3;
 
}



void AliITSmoduleSSD::HitToDigit(Int_t hitNo) {
    
    Int_t stripP, stripN, i;
    Float_t dsP, dsN;
    Float_t sP, sN;
    Float_t EP, EN;
    AliITShit *hit = (AliITShit*)((*fHitsM)[hitNo]);
    Float_t dZ = kZ/kSteps*1000, l;

    if(hit->GetIonization()==0.0) return;

    Float_t x =  hit->GetXG();
    Float_t y =  hit->GetYG();

    Float_t dx = 0.0; //TMath::Tan(hit->fTheta)*kZ/kSteps;
    Float_t dy = 0.0; //TMath::Tan(hit->fPhi)*kZ/kSteps;
    l = sqrt(dZ*dZ + dx*dx *1000000);
    
    x -= (kSteps/2 -1) * dx;
    y -= (kSteps/2 -1) * dy;
    
    for (i=1; i<kSteps; i++) {
    
        stripP = GetStripP(x, y);
        dsP = Get2StripP(x, y);
        
        stripN = GetStripN(x, y);
        dsN = Get2StripN(x, y);
        
        EP = gRandom->Landau(fGainP*l, l*10);
        EN = gRandom->Landau(fGainN*l, l*10);
        
        sP = kSigmaP * sqrt(i);
        sN = kSigmaN * sqrt(kSteps-i);

        sP = (i<3 && dsP>0.3 && dsP<0.7)? 0.02 : sP;
        sN = (i>7 && dsN>0.3 && dsN<0.7)? 0.02 : sN;         

        sP = (i==3 && dsP>0.4 && dsP<0.6)? 0.015 : sP;
        sN = (i==7 && dsN>0.4 && dsN<0.6)? 0.015 : sN;        
        
        
        (*fP)[stripP-1]+=(Int_t)(EP*(F(-0.5-dsP,sP)-F(-1.5-dsP,sP)));  
        (*fP)[stripP]  +=(Int_t)(EP*(F(0.5-dsP,sP)-F(-0.5-dsP,sP)));
        (*fP)[stripP+1]+=(Int_t)(EP*(F(1.5-dsP,sP)-F(0.5-dsP,sP)));
        (*fP)[stripP+2]+=(Int_t)(EP*(F(2.5-dsP,sP)-F(1.5-dsP,sP))); 
        
        (*fN)[stripN-1]+=(Int_t)(EN*(F(-0.5-dsN,sN)-F(-1.5-dsN,sN)));
        (*fN)[stripN]  +=(Int_t)(EN*(F(0.5-dsN,sN)-F(-0.5-dsN,sN)));
        (*fN)[stripN+1]+=(Int_t)(EN*(F(1.5-dsN,sN)-F(0.5-dsN,sN)));
        (*fN)[stripN+2]+=(Int_t)(EN*(F(2.5-dsN,sN)-F(1.5-dsN,sN))); 
        
        x += dx; 
        y += dy; 
    }
}


//____________________________________________________________________
//
//  Private Methods for Simulation
//____________________________________________________________________
//                           



void AliITSmoduleSSD::ApplyNoise() {
    
    for(Int_t i = 0; i<kNStrips; i++) {
        (*fP)[i] += (Int_t)gRandom->Gaus(0,fSNRatioP);
        (*fN)[i] += (Int_t)gRandom->Gaus(0,fSNRatioN);
    }
}

void AliITSmoduleSSD::ApplyCoupling() {
    
    for(Int_t i = 1; i<kNStrips-1; i++) {
        (*fP)[i] += (Int_t)((*fP)[i-1]*fCouplingPL + (*fP)[i+1]*fCouplingPR);
        (*fN)[i] += (Int_t)((*fN)[i-1]*fCouplingNL + (*fN)[i+1]*fCouplingNR);
    }
}


//____________________________________________________________________
//
//  Private methods for geometry
//____________________________________________________________________
//                           



Int_t AliITSmoduleSSD::GetStripP(Float_t x, Float_t y) {
    
    Float_t  X = x - y*kTan;
    Int_t strip = (Int_t)(X/kPitch);
    strip = (strip<0)? -1: strip;
    strip = (strip>kNStrips)? -1: strip;
    return strip;
}

Int_t AliITSmoduleSSD::GetStripN(Float_t x, Float_t y) {
    
    Float_t  X = x - kTan*(kY - y);
    Int_t strip = (Int_t)(X/kPitch);
    strip = (strip<0)? -1: strip;
    strip = (strip>kNStrips)? -1: strip;
    return strip;

}

Float_t AliITSmoduleSSD::Get2StripP(Float_t x, Float_t y) {
    
    Int_t n = GetStripP(x,y);
    return (x - y*kTan) / kPitch - n;
}

Float_t AliITSmoduleSSD::Get2StripN(Float_t x, Float_t y) {
    
    Int_t n = GetStripN(x,y);
    return (x - kTan*(kY - y)) / kPitch - n;
}


Bool_t AliITSmoduleSSD::GetCrossing (Float_t &P, Float_t &N) {   

    P *= kPitch;
    N *= kPitch; 
    
    P = (kY * kTan + N + P)/2.0;         // x coordinate
    N = kY - (P-N)/kTan;                 // y coordinate
    
    if (N<0 || N>kY) return kFALSE;
    if (P<0 || P>kX) return kFALSE;
    return kTRUE;   
}

//____________________________________________________________________



