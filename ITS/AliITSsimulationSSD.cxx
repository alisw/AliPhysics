
#include <stdio.h>
#include <TObjArray.h>

#include "AliITSsimulationSSD.h"
#include "AliITSdictSSD.h"
#include "AliITSdcsSSD.h"
#include "AliITS.h"
#include "AliRun.h"


ClassImp(AliITSsimulationSSD);
//------------------------------------------------------------
AliITSsimulationSSD::AliITSsimulationSSD(AliITSsegmentation *seg,
                                         AliITSresponse *resp){
  // Constructor

    fSegmentation = seg;
    fResponse = resp;
    fDCS = new AliITSdcsSSD(seg,resp); 

    fNstrips = fSegmentation->Npx();
    fPitch = fSegmentation->Dpx(0);
    
    fP = new TArrayF(fNstrips+1); 
    fN = new TArrayF(fNstrips+1);
     
    fTracksP = new AliITSdictSSD[fNstrips+1];
    fTracksN = new AliITSdictSSD[fNstrips+1];

    
    fSteps  = 10;   // still hard-wired - set in SetDetParam and get it via  
                     // fDCS together with the others eventually    

}
//___________________________________________________________________________
AliITSsimulationSSD& AliITSsimulationSSD::operator=(AliITSsimulationSSD 
                                                                      &source){
// Operator =
    if(this==&source) return *this;

    this->fDCS = new AliITSdcsSSD(*(source.fDCS));
    this->fN   = new TArrayF(*(source.fN));
    this->fP   = new TArrayF(*(source.fP));
    this->fTracksP = new AliITSdictSSD(*(source.fTracksP));
    this->fTracksN = new AliITSdictSSD(*(source.fTracksN));
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
  // anihilator    

    if(fP) delete fP;
    if(fN) delete fN;
    
    if(fTracksP) delete fTracksP;
    if(fTracksN) delete fTracksN;

    delete fDCS;
} 
//_______________________________________________________________
//
// Hit to digit
//_______________________________________________________________
//
void AliITSsimulationSSD::DigitiseModule(AliITSmodule *mod,Int_t module,
                                         Int_t dummy) {
  // Digitizes one SSD module of hits.
    
    TObjArray *hits = mod->GetHits();
    Int_t nhits = hits->GetEntriesFast();
    if (!nhits) return;

    Int_t i;
    for(i=0; i<fNstrips; i++) {
       (*fP)[i] = 0;
       (*fN)[i] = 0;
       fTracksP[i].ZeroTracks();
       fTracksN[i].ZeroTracks();
    }
    
    for(i=0; i<nhits; i++) {
       Int_t idtrack=mod->GetHitTrackIndex(i);  
       HitToDigit(i,idtrack,nhits,hits);
    }
 
   

    ApplyNoise();
    ApplyCoupling();    
    ApplyThreshold();
    ApplyDAQ();
      

}

//---------------------------------------------------------------

void AliITSsimulationSSD::HitToDigit(Int_t & hitNo,Int_t idtrack,
                                     Int_t nhits,TObjArray *hits) {
  // Turns one or more hits in an SSD module into one or more digits.
     
    Int_t      stripP, stripN, i;
    Float_t    dsP, dsN;
    Float_t    sP, sN;
    Float_t    eP, eN;
    Float_t    arrayEP[10];         // hard-wired number of steps
    Float_t    arrayEN[10];
    Int_t      track = -1;
       
    Float_t    ionization = 0;
    Float_t    signal;
    
    AliITSdictSSD *dict;  
    
 
    // check if this is the right order !!!!!

    AliITShit *hitI = (AliITShit*)hits->At(hitNo++);
    AliITShit *hitE = (AliITShit*)hits->At(hitNo);


    while (!((hitE->StatusExiting()) || 
             (hitE->StatusDisappeared()) ||
             (hitE->StatusStop()))) {
       
        if (++hitNo<nhits) {
           ionization = hitE->GetIonization(); 
           hitE = (AliITShit*)hits->At(hitNo);
        }
    }   
    
        
    if (hitI->GetTrack() == hitE->GetTrack()) 
       //track = idtrack;
       track = hitI->GetTrack();
    else 
       printf("!!! Emergency !!!\n");
     
            
    ionization += hitE->GetIonization();
       
    const Float_t kconvm=10000.;  // cm -> microns

    Float_t xI, yI, zI;
    hitI->GetPositionL(xI, yI, zI);
    
    xI *= kconvm;
    yI *= kconvm;
    zI *= kconvm;
    
    Float_t xE, yE, zE;
    hitE->GetPositionL(xE, yE, zE);
    
    xE *= kconvm;
    yE *= kconvm;
    zE *= kconvm;

    Float_t dx = (xE - xI);
    Float_t dz = (zE - zI);
              
    
    // Debuging
    /*
    fSegmentation->GetPadIxz(xI,zI,stripP,stripN);
   
       printf("%5d %8.3f %8.3f %8.3f %8.3f %d %d  %d\n", 
             hitNo, xI, zI, dx, dz, 
             stripP, stripN, track);
     printf("%10.5f %10d \n", ionization, hitI->fTrack); 
    */ 
    
    // end of debuging   
    
    
    eP=0;
    eN=0;
    
    for (i=0; i<fSteps; i++) {
        
      //        arrayEP[i] = gRandom->Landau(ionization/fSteps, ionization/(4*fSteps));
      //        arrayEN[i] = gRandom->Landau(ionization/fSteps, ionization/(4*fSteps));
        arrayEP[i] = ionization/fSteps;
        arrayEN[i] = ionization/fSteps;
       
        eP += arrayEP[i];
        eN += arrayEN[i];
    } 
       
    const Float_t kconv = 1.0e9 / 3.6;  // GeV -> e-hole pairs
       
    for(i=0; i<fSteps; i++) {
    
        arrayEP[i] = kconv * arrayEP[i] * (ionization / eP); 
        arrayEN[i] = kconv * arrayEN[i] * (ionization / eN);        
    }    
        
    dx /= fSteps;
    dz /= fSteps;  

    Float_t sigmaP, sigmaN; 
    fResponse->SigmaSpread(sigmaP,sigmaN); 

    //printf("SigmaP SigmaN %f %f\n",sigmaP, sigmaN);

    Float_t noiseP, noiseN;
    fResponse->GetNoiseParam(noiseP,noiseN);

    //printf("NoiseP NoiseN %f %f\n",noiseP, noiseN);

    for(i=0; i<fSteps; i++) {
        
        Int_t j;
       
        fSegmentation->GetPadIxz(xI,zI,stripP,stripN);
        //printf("hitNo %d i xI zI stripP stripN %d %f %f %d %d\n",hitNo,i,xI, zI, stripP, stripN);
        dsP    = Get2Strip(1,stripP,xI, zI); // Between 0-1
        dsN    = Get2Strip(0,stripN,xI, zI); // Between 0-1

        sP = sigmaP * sqrt(300. * i / (fSteps));
        sN = sigmaN * sqrt(300. * i /(fSteps-i));


        sP = (i<2        && dsP>0.3 && dsP<0.7)? 20. : sP;  // square of (microns) 
        sN = (i>fSteps-2 && dsN>0.3 && dsN<0.7)? 20. : sN;  // square of (microns) 

        sP = (i==2 && dsP>0.4 && dsP<0.6)? 15. : sP;  // square of (microns) 
        sN = (i==8 && dsN>0.4 && dsN<0.6)? 15. : sN;  // square of (microns)        
        
        for (j=-1; j<2; j++) {
            if (stripP+j<0 || stripP+j>fNstrips) continue;
            signal = arrayEP[i] * TMath::Abs( (F(j+0.5-dsP,sP)-F(j-0.5-dsP,sP)) );
            //printf("SimSSD::HitsToDigits:%d arrayEP[%d]=%e signal=%e\n",j,i,arrayEP[i],signal);
	    if (signal > noiseP/fSteps) {
               (*fP)[stripP+j] += signal;
               dict = (fTracksP+stripP+j);   
               (*dict).AddTrack(track);
            } 
	}  // end for j loop over neighboring strips

        for (j=-1; j<2; j++) {
            if (stripN+j<0 || stripN+j>fNstrips) continue;
            signal = arrayEN[i] * TMath::Abs( (F(j+0.5-dsN,sN)-F(j-0.5-dsN,sN)) );
            //printf("SimSSD::HitsToDigits:%d arrayEN[%d]=%e signal=%e\n",j,i,arrayEN[i],signal);
            if (signal > noiseN/fSteps) {
               (*fN)[stripN+j] += signal;
               dict = (fTracksN+stripN+j);    //co to jest
               (*dict).AddTrack(track);
            } 
        }  // end for j loop over neighboring strips
                
        xI += dx; 
        zI += dz; 
    }
    
    
}


//____________________________________________________________________
//
//  Private Methods for Simulation
//______________________________________________________________________
//

void AliITSsimulationSSD::ApplyNoise() {
  // Apply Noise.
   Float_t noiseP, noiseN;
   fResponse->GetNoiseParam(noiseP,noiseN);
       
    Int_t i;
    for(i = 0; i<fNstrips; i++) {
       (*fP)[i] += gRandom->Gaus(0,noiseP);
       (*fN)[i] += gRandom->Gaus(0,noiseN);
    }
}

//_________________________________________________________________________

void AliITSsimulationSSD::ApplyCoupling() {
  // Apply the effecto of electronic coupling between channels    
    Int_t i;
    for(i = 1; i<fNstrips-1; i++) {
      (*fP)[i] += (*fP)[i-1]*fDCS->GetCouplingPL() + (*fP)[i+1]*fDCS->GetCouplingPR();
      (*fN)[i] += (*fN)[i-1]*fDCS->GetCouplingNL() + (*fN)[i+1]*fDCS->GetCouplingNR();
    }
}

//__________________________________________________________________________

void AliITSsimulationSSD::ApplyThreshold() {
  // Applies the effect of a threshold on the signals for digitization.
   Float_t noiseP, noiseN;
   fResponse->GetNoiseParam(noiseP,noiseN);

   // or introduce the SetThresholds in fResponse  

    Int_t i;
    for(i=0; i<fNstrips; i++) {
       (*fP)[i] = ((*fP)[i] > noiseP*4) ? (*fP)[i] : 0;
       (*fN)[i] = ((*fN)[i] > noiseN*4) ? (*fN)[i] : 0; 
    }
        
}

//__________________________________________________________________________

void AliITSsimulationSSD::ApplyDAQ() {
  // Converts simulated signals to simulated ADC counts
    AliITS *its=(AliITS*)gAlice->GetModule("ITS");

    Float_t noiseP, noiseN;
    fResponse->GetNoiseParam(noiseP,noiseN);

    char opt[30],dummy[20];
    fResponse->ParamOptions(opt,dummy);

    Int_t i,j;
    if (strstr(opt,"SetInvalid")) {
      // Set signal = 0 if invalid strip
      for(i=0; i<fNstrips; i++) {
         if (!(fDCS->IsValidP(i))) (*fP)[i] = 0;
	 if (!(fDCS->IsValidN(i))) (*fN)[i] = 0;
      }
    }
    
    Int_t digits[3], tracks[3], hits[3];
    Float_t charges[3];
    Float_t phys=0;
    for(i=0;i<3;i++) tracks[i]=0;
    for(i=0; i<fNstrips; i++) { 
       if( (strstr(opt,"SetInvalid") && (*fP)[i] < noiseP*4) || !(*fP)[i]) continue;
          digits[0]=1;
          digits[1]=i;
          digits[2]=(int)(*fP)[i];
          for(j=0; j<(fTracksP+i)->GetNTracks(); j++) {
            if(j>2) continue;
	    if((fTracksP+i)->GetNTracks()) tracks[j]=(fTracksP+i)->GetTrack(j);
	    else tracks[j]=-2;
	    charges[j] = 0;
	    hits[j] = -1;
          }
          its->AddSimDigit(2,phys,digits,tracks,hits,charges);
          
    }
    
    
    for(i=0; i<fNstrips; i++) {
       if( (strstr(opt,"SetInvalid") && (*fN)[i] < noiseN*4)|| !(*fN)[i]) continue;
          digits[0]=0;
          digits[1]=i;
          digits[2]=(int)(*fN)[i];
          for( j=0; j<(fTracksN+i)->GetNTracks(); j++) {
            if(j>2) continue;
            if((fTracksN+i)->GetNTracks()) tracks[j]=(fTracksN+i)->GetTrack(j);
	    else tracks[j]=-2;
            charges[j] = 0;
	    hits[j] = -1;
          }
          its->AddSimDigit(2,phys,digits,tracks,hits,charges);
          
    }
    
}


//____________________________________________________________________________

Float_t AliITSsimulationSSD::F(Float_t x, Float_t s) {
  // Computes the integral of a gaussian at the mean valuse x with sigma s.

    //printf("SDD:F(%e,%e)\n",x,s);
    return 0.5*TMath::Erf(x * fPitch / s) ;
} 

//______________________________________________________________________

Float_t AliITSsimulationSSD::Get2Strip(Int_t flag, Int_t iStrip, Float_t x, Float_t z){
  // Returns the relative space between two strips.

    // flag==1 for Pside, 0 for Nside

    Float_t stereoP, stereoN;
    fSegmentation->Angles(stereoP,stereoN);
    
    Float_t tanP=TMath::Tan(stereoP);
    Float_t tanN=TMath::Tan(stereoN);
 
    Float_t dx = fSegmentation->Dx();
    Float_t dz = fSegmentation->Dz();


     x += dx/2;
     z += dz/2; 
    
     if (flag) return (x - z*tanP) / fPitch - iStrip;       // from 0 to 1
     else  return (x - tanN*(dz - z)) / fPitch - iStrip;
}
//____________________________________________________________________________
