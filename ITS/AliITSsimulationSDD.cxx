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


#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <TSystem.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TH1.h>
#include <TFile.h>
#include <TVector.h>
#include <TArrayI.h>
#include <TArrayF.h>

#include "AliRun.h"
#include "AliITS.h"
#include "AliITShit.h"
#include "AliITSdigit.h"
#include "AliITSmodule.h"
#include "AliITSMapA1.h"
#include "AliITSMapA2.h"
#include "AliITSetfSDD.h"
#include "AliITSRawData.h"
#include "AliITSHuffman.h"
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"
#include "AliITSsimulationSDD.h"

ClassImp(AliITSsimulationSDD)
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Piergiorgio Cerello
// November 23 1999
//
// AliITSsimulationSDD is the simulation of SDDs.
  //
//Begin_Html
/*
<img src="picts/ITS/AliITShit_Class_Diagram.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This show the relasionships between the ITS hit class and the rest of Aliroot.
</font>
<pre>
*/
//End_Html
//_____________________________________________________________________________

Int_t power(Int_t b, Int_t e) {
  // compute b to the e power, where both b and e are Int_ts.
  Int_t power = 1,i;
  for(i=0; i<e; i++) power *= b;
  return power;
}

//_____________________________________________

void FastFourierTransform(AliITSetfSDD *alisddetf,Double_t *real,
                          Double_t *imag,Int_t direction) {
  // Do a Fast Fourier Transform
  //printf("FFT: direction %d\n",direction);

  Int_t samples = alisddetf->GetSamples();
  Int_t l = (Int_t) ((log((Float_t) samples)/log(2.))+0.5);
  Int_t m1 = samples;
  Int_t m  = samples/2;
  Int_t m2 = samples/m1;
  Int_t i,j,k;
  for(i=1; i<=l; i++) {
    for(j=0; j<samples; j += m1) {
      Int_t p = 0;
      for(k=j; k<= j+m-1; k++) {
	Double_t wsr = alisddetf->GetWeightReal(p);
        Double_t wsi = alisddetf->GetWeightImag(p);
	if(direction == -1) wsi = -wsi;
	Double_t xr = *(real+k+m);
	Double_t xi = *(imag+k+m);
	*(real+k+m) = wsr*(*(real+k)-xr) - wsi*(*(imag+k)-xi);
	*(imag+k+m) = wsr*(*(imag+k)-xi) + wsi*(*(real+k)-xr);
	*(real+k) += xr;
	*(imag+k) += xi;
	p += m2;
      }
    }
    m1 = m;
    m /= 2;
    m2 += m2;
  } 
  
  for(j=0; j<samples; j++) {
    Int_t j1 = j;
    Int_t p = 0;
    Int_t i1;
    for(i1=1; i1<=l; i1++) {
      Int_t j2 = j1;
      j1 /= 2;
      p = p + p + j2 - j1 - j1;
    }
    if(p >= j) {
      Double_t xr = *(real+j);
      Double_t xi = *(imag+j);
      *(real+j) = *(real+p);
      *(imag+j) = *(imag+p);
      *(real+p) = xr;
      *(imag+p) = xi;
    }
  }
  if(direction == -1) {
    for(i=0; i<samples; i++) {
      *(real+i) /= samples;
      *(imag+i) /= samples;
    }
  }
  return;
}
//_____________________________________________________________________________

AliITSsimulationSDD::AliITSsimulationSDD(){
  // Default constructor

  fResponse = 0;
  fSegmentation = 0;
  fHis = 0;
  fHitMap1 = 0;
  fHitMap2 = 0;
  fElectronics = 0;
  fStream = 0;
  fD.Set(0);
  fT1.Set(0);
  fT2.Set(0);
  fTol.Set(0);
  fNoise.Set(0);
  fBaseline.Set(0);
  SetScaleFourier();
  SetPerpendTracksFlag();
  SetDoFFT();
  SetCheckNoise();
  fInZR = 0;
  fInZI = 0;
  fOutZR = 0;
  fOutZI = 0;
  fNofMaps = 0;
  fMaxNofSamples = 0;
  fITS = 0;
  fTreeB=0;
}
//_____________________________________________________________________________
AliITSsimulationSDD::AliITSsimulationSDD(AliITSsimulationSDD &source)
{
  // Copy constructor to satify Coding roules only.
  if(this==&source) return;
  printf("Not allowed to make a copy of AliITSsimulationSDD "
         "Using default creater instead\n");
  AliITSsimulationSDD();
}
//_____________________________________________________________________________
AliITSsimulationSDD& AliITSsimulationSDD::operator=(AliITSsimulationSDD &source)
{
  // Assignment operator to satify Coding roules only.
  if(this==&source) return *this;
  printf("Not allowed to make a = with AliITSsimulationSDD "
         "Using default creater instead\n");
  return *this ;
}
//_____________________________________________________________________________

AliITSsimulationSDD::AliITSsimulationSDD(AliITSsegmentation *seg,AliITSresponse *resp) 
{
  // Standard Constructor

      fHis=0;
      fTreeB=0;
      fResponse = resp;
      fSegmentation = seg;
      SetScaleFourier();
      SetPerpendTracksFlag();
      SetDoFFT();
      SetCheckNoise();

      fHitMap2 = new AliITSMapA2(fSegmentation,fScaleSize,1);
      fHitMap1 = new AliITSMapA1(fSegmentation);

      //
      fNofMaps=fSegmentation->Npz();
      fMaxNofSamples=fSegmentation->Npx();

      Float_t sddLength = fSegmentation->Dx();
      Float_t sddWidth = fSegmentation->Dz();

      Int_t dummy=0;
      Float_t anodePitch = fSegmentation->Dpz(dummy);
      Double_t timeStep = (Double_t)fSegmentation->Dpx(dummy);
      Float_t driftSpeed=fResponse->DriftSpeed();    

      if(anodePitch*(fNofMaps/2) > sddWidth) {
         Warning("AliITSsimulationSDD",
           "Too many anodes %d or too big pitch %f \n",fNofMaps/2,anodePitch);
      }

      if(timeStep*fMaxNofSamples < sddLength/driftSpeed) {
         Error("AliITSsimulationSDD",
                             "Time Interval > Allowed Time Interval: exit\n");
         return;
      }

      fElectronics = new AliITSetfSDD(timeStep/fScaleSize);

      char opt1[20], opt2[20];
      fResponse->ParamOptions(opt1,opt2);
      fParam=opt2;
      char *same = strstr(opt1,"same");
      if (same) {
         fNoise.Set(0);
         fBaseline.Set(0);
      } else {
         fNoise.Set(fNofMaps);
         fBaseline.Set(fNofMaps);
      }
      
      //
      const char *kopt=fResponse->ZeroSuppOption();
        if (strstr(fParam,"file") ) {
	  fD.Set(fNofMaps);
	  fT1.Set(fNofMaps);
          if (strstr(kopt,"2D")) {
	    fT2.Set(fNofMaps);
            fTol.Set(0);
            Init2D();       // desactivate if param change module by module
          } else if(strstr(kopt,"1D"))  {
            fT2.Set(2);
            fTol.Set(2);
            Init1D();      // desactivate if param change module by module
	  }
	} else {
          fD.Set(2);
	  fTol.Set(2);
	  fT1.Set(2);
	  fT2.Set(2);
	  SetCompressParam();
	}


	Bool_t write=fResponse->OutputOption();
	if(write && strstr(kopt,"2D")) MakeTreeB();

        // call here if baseline does not change by module
        // ReadBaseline();

        fITS = (AliITS*)gAlice->GetModule("ITS");
        Int_t size=fNofMaps*fMaxNofSamples;
	fStream = new AliITSInStream(size); 
	
	fInZR = new Double_t [fScaleSize*fMaxNofSamples];
	fInZI = new Double_t [fScaleSize*fMaxNofSamples];
	fOutZR = new Double_t [fScaleSize*fMaxNofSamples];
	fOutZI = new Double_t [fScaleSize*fMaxNofSamples];  

}


//_____________________________________________________________________________

AliITSsimulationSDD::~AliITSsimulationSDD() { 
  // destructor

  delete fHitMap1;
  delete fHitMap2;
  delete fStream;
  delete fElectronics;

  fD.Set(0);
  fT1.Set(0);
  fT2.Set(0);
  fTol.Set(0);
  fNoise.Set(0);
  fBaseline.Set(0);
  fITS = 0;

  if (fHis) {
     fHis->Delete(); 
     delete fHis;     
  }     
  if(fTreeB) delete fTreeB;           
  if(fInZR) delete [] fInZR;
  if(fInZI) delete [] fInZI;	
  if(fOutZR) delete [] fOutZR;
  if(fOutZI) delete [] fOutZI;
}
//_____________________________________________________________________________

void AliITSsimulationSDD::DigitiseModule(AliITSmodule *mod,Int_t md,Int_t ev){
  // create maps to build the lists of tracks
  // for each digit

    fModule=md;
    fEvent=ev;

    TObjArray *fHits = mod->GetHits();
    Int_t nhits = fHits->GetEntriesFast();
    if (!nhits && fCheckNoise) {
        ChargeToSignal();
        GetNoise();
	fHitMap2->ClearMap();
        return;
    } else if (!nhits) return;

    //printf("simSDD: module nhits %d %d\n",md,nhits);


    TObjArray *list=new TObjArray;
    static TClonesArray *padr=0;
    if(!padr) padr=new TClonesArray("TVector",1000);
    Int_t arg[6] = {0,0,0,0,0,0}; 
    fHitMap1->SetArray(list);


    Int_t nofAnodes=fNofMaps/2;

    Float_t sddLength = fSegmentation->Dx();
    Float_t sddWidth = fSegmentation->Dz();

    Int_t dummy=0;
    Float_t anodePitch = fSegmentation->Dpz(dummy);
    Float_t timeStep = fSegmentation->Dpx(dummy);

    Float_t driftSpeed=fResponse->DriftSpeed();    

    // Piergiorgio's part (apart for few variables which I made float
    // when i thought that can be done

    // Fill detector maps with GEANT hits
    // loop over hits in the module

    const Float_t kconv=1.0e+6;  // GeV->KeV
    Int_t ii;
    Int_t idhit=-1;
    for(ii=0; ii<nhits; ii++) {
      AliITShit *hit = (AliITShit*) fHits->At(ii);
      Float_t xL[3];
      hit = (AliITShit*) fHits->At(ii);
      hit->GetPositionL(xL[0],xL[1],xL[2]);
      Int_t hitDetector = hit->GetDetector();

      if(hit->StatusEntering()) idhit=ii;
      
      Int_t nOfSplits = 5;
      if(fFlag) nOfSplits = 1;
      // Deposited energy in keV
      Float_t avpath = 0.;
      Float_t avanod = 0.;
      Float_t depEnergy = kconv*hit->GetIonization()/nOfSplits;
      AliITShit *hit1 = 0;
      Float_t xL1[3];
      if(fFlag && depEnergy != 0.) continue;
      if(depEnergy == 0.) {	
	  ii++;
	  hit1 = (AliITShit*) fHits->At(ii);
	  hit1->GetPositionL(xL1[0],xL1[1],xL1[2]);
      } else {
	  xL1[0] = xL[0];
	  xL1[1] = xL[1];
	  xL1[2] = xL[2];
      }

      // scale path to simulate a perpendicular track

      if(depEnergy == 0.) depEnergy = kconv*hit1->GetIonization()/nOfSplits;
      // continue if the particle did not lose energy
      // passing through detector
      if (!depEnergy) {
	  printf("This particle has passed without losing energy!\n");
	  continue;
      }

      if (fFlag) {
	  Float_t pathInSDD = TMath::Sqrt((xL[0]-xL1[0])*(xL[0]-xL1[0])+(xL[1]-xL1[1])*(xL[1]-xL1[1])+(xL[2]-xL1[2])*(xL[2]-xL1[2]));
	  if(pathInSDD) depEnergy *= (0.03/pathInSDD); 
      }
	
      for(Int_t kk=0;kk<nOfSplits;kk++) {
	Float_t avDrft =  
                xL[0]+(xL1[0]-xL[0])*((kk+0.5)/((Float_t) nOfSplits));
	Float_t avAnode = 
                xL[2]+(xL1[2]-xL[2])*((kk+0.5)/((Float_t) nOfSplits));
	Float_t driftPath = 10000.*avDrft;
	  
	avpath = xL1[0];
	avanod = xL1[2];

	Int_t iWing = 2;
	if(driftPath < 0) {
	   iWing = 1;
	   driftPath = -driftPath;
	}
	driftPath = sddLength-driftPath;
	Int_t detector = 2*(hitDetector-1) + iWing;
	if(driftPath < 0) {
	   cout << "Warning: negative drift path " << driftPath << endl;
	   continue;
	}
	 
	//   Drift Time
	Float_t driftTime = driftPath/driftSpeed;
	Int_t timeSample = (Int_t) (fScaleSize*driftTime/timeStep + 1);
	if(timeSample > fScaleSize*fMaxNofSamples) {
	   cout << "Warning: Wrong Time Sample: " << timeSample << endl;
	   continue;
	}

	//   Anode
	Float_t xAnode = 10000.*(avAnode)/anodePitch + nofAnodes/2;  // +1?
	if((xAnode+1)*anodePitch > sddWidth || xAnode*anodePitch < 0.) 
	  { cout << "Warning: Z = " << xAnode*anodePitch << endl; }
	Int_t iAnode = (Int_t) (1.+xAnode); // xAnode?
	if(iAnode < 0 || iAnode > nofAnodes) {
	  cout << "Warning: Wrong iAnode: " << iAnode << endl;
	  continue;
	} 


	// work with the idtrack=entry number in the TreeH for the moment
	//Int_t idhit,idtrack;
	//mod->GetHitTrackAndHitIndex(ii,idtrack,idhit);    
	//Int_t idtrack=mod->GetHitTrackIndex(ii);  
        // or store straight away the particle position in the array
	// of particles and take idhit=ii only when part is entering (this
	// requires FillModules() in the macro for analysis) : 
	Int_t itrack = hit->GetTrack();

	//  Signal 2d Shape
	Float_t diffCoeff, s0;
	fResponse->DiffCoeff(diffCoeff,s0);
    
	// Squared Sigma along the anodes
	Double_t sigma2A = 2.*diffCoeff*driftTime+s0*s0;
	Double_t sigmaA  = TMath::Sqrt(sigma2A);
	Double_t sigmaT  = sigmaA/driftSpeed;
    
	// Peak amplitude in nanoAmpere
	Double_t eVpairs = 3.6;
	Double_t amplitude = fScaleSize*160.*depEnergy/(timeStep*eVpairs*2.*acos(-1.)*sigmaT*sigmaA);
    
	Float_t nsigma=fResponse->NSigmaIntegration();
	// Spread the charge 
	// Pixel index
	Int_t ja = iAnode;
	Int_t jt = timeSample;
	// Sub-pixel index
	Int_t nsplit = 4; // hard-wired
	nsplit = (nsplit+1)/2*2;
	// Sub-pixel size
	Double_t aStep = anodePitch/(nsplit*fScaleSize);
	Double_t tStep = timeStep/(nsplit*fScaleSize);
	// Define SDD window corresponding to the hit
	Int_t anodeWindow = (Int_t) (fScaleSize*nsigma*sigmaA/anodePitch + 1);
	Int_t timeWindow = (Int_t) (fScaleSize*nsigma*sigmaT/timeStep + 1);
	Int_t jamin = (ja - anodeWindow/2 - 1)*fScaleSize*nsplit + 1;
	Int_t jamax = (ja + anodeWindow/2)*fScaleSize*nsplit;
	if(jamin <= 0) jamin = 1;
	if(jamax > fScaleSize*nofAnodes*nsplit) jamax = fScaleSize*nofAnodes*nsplit;
	Int_t jtmin = (jt - timeWindow*3 - 1)*nsplit + 1; //hard-wired
	Int_t jtmax = (jt + timeWindow*3)*nsplit; //hard-wired
	if(jtmin <= 0) jtmin = 1;
	if(jtmax > fScaleSize*fMaxNofSamples*nsplit) jtmax = fScaleSize*fMaxNofSamples*nsplit;

	Double_t rlAnode = log(aStep*amplitude);

	// Spread the charge in the anode-time window
        Int_t ka;
	for(ka=jamin; ka <=jamax; ka++) {
	  Int_t ia = (ka-1)/(fScaleSize*nsplit) + 1;
	  if(ia <= 0) { cout << "Warning: ia < 1: " << endl; continue; }
	  if(ia > nofAnodes) ia = nofAnodes;
	  Double_t aExpo = (aStep*(ka-0.5)-xAnode*anodePitch)/sigmaA;
	  Double_t anodeAmplitude = rlAnode - 0.5*aExpo*aExpo;
	  // Protect against overflows
	  if(anodeAmplitude > -87.3)
	    anodeAmplitude = exp(anodeAmplitude);
	  else
	    anodeAmplitude = 0;
	  Int_t index = ((detector+1)%2)*nofAnodes+ia-1; // index starts from 0
	  if(anodeAmplitude) {
	    Double_t rlTime = log(tStep*anodeAmplitude);
            Int_t kt;
	    for(kt=jtmin; kt<=jtmax; kt++) {
	      Int_t it = (kt-1)/nsplit+1;  // it starts from 1
	      if(it<=0) { cout << "Warning: it < 1: " << endl; continue; } 
	      if(it>fScaleSize*fMaxNofSamples) it = fScaleSize*fMaxNofSamples;
	      Double_t tExpo = (tStep*(kt-0.5)-driftTime)/sigmaT;
	      Double_t timeAmplitude = rlTime - 0.5*tExpo*tExpo;
	      // Protect against overflows
	      if(timeAmplitude > -87.3){
		timeAmplitude = exp(timeAmplitude);
	      } else
		timeAmplitude = 0;

	      // build the list of digits for this module	
	      arg[0]=index;
	      arg[1]=it;
	      arg[2]=itrack;
	      arg[3]=idhit;
	      ListOfFiredCells(arg,timeAmplitude,list,padr);
	} // loop over time in window 
      } // end if anodeAmplitude
    } // loop over anodes in window
   } // end loop over "sub-hits"
  } // end loop over hits

  // introduce the electronics effects and do zero-suppression if required
  Int_t nentries=list->GetEntriesFast();
  if (nentries) {

    //TStopwatch timer;
    ChargeToSignal(); 
    //timer.Stop(); timer.Print();

    const char *kopt=fResponse->ZeroSuppOption();
    ZeroSuppression(kopt);
  } 

  // clean memory
  list->Delete();
  delete list; 
                      
  padr->Delete(); 

  fHitMap1->ClearMap();
  fHitMap2->ClearMap();

  //gObjectTable->Print();
}


//____________________________________________

void AliITSsimulationSDD::ListOfFiredCells(Int_t *arg,Double_t timeAmplitude,
                                           TObjArray *list,TClonesArray *padr){
  // Returns the list of "fired" cells.

                    Int_t index=arg[0];
                    Int_t ik=arg[1];
                    Int_t idtrack=arg[2];
                    Int_t idhit=arg[3];
                    Int_t counter=arg[4];
                    Int_t countadr=arg[5];
                   
                    Double_t charge=timeAmplitude;
		    charge += fHitMap2->GetSignal(index,ik-1);
		    fHitMap2->SetHit(index, ik-1, charge);

                    Int_t digits[3];
		    Int_t it=(Int_t)((ik-1)/fScaleSize);
		    
		    digits[0]=index;
		    digits[1]=it;
		    digits[2]=(Int_t)timeAmplitude;
                    Float_t phys;
		    if (idtrack >= 0) phys=(Float_t)timeAmplitude;
		    else phys=0;
                   
		    Double_t cellcharge=0.;
		    AliITSTransientDigit* pdigit;
		    // build the list of fired cells and update the info
		    if (!fHitMap1->TestHit(index, it)) {
		      
		        new((*padr)[countadr++]) TVector(3);
			TVector &trinfo=*((TVector*) (*padr)[countadr-1]);
			trinfo(0)=(Float_t)idtrack;
			trinfo(1)=(Float_t)idhit;
			trinfo(2)=(Float_t)timeAmplitude;

			list->AddAtAndExpand(
			    new AliITSTransientDigit(phys,digits),counter);
			
			fHitMap1->SetHit(index, it, counter);
			counter++;
			pdigit=(AliITSTransientDigit*)list->
                                                      At(list->GetLast());
			// list of tracks
			TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			trlist->Add(&trinfo);

		    } else {
			pdigit=
                         (AliITSTransientDigit*) fHitMap1->GetHit(index, it);
			for(Int_t kk=0;kk<fScaleSize;kk++) {
			  cellcharge += fHitMap2->GetSignal(index,fScaleSize*it+kk);
			}
			// update charge
			(*pdigit).fSignal=(Int_t)cellcharge;
			(*pdigit).fPhysics+=phys;			
			// update list of tracks
			TObjArray* trlist=(TObjArray*)pdigit->TrackList();
			Int_t lastentry=trlist->GetLast();
			TVector *ptrkp=(TVector*)trlist->At(lastentry);
			TVector &trinfo=*ptrkp;
			Int_t lasttrack=Int_t(trinfo(0));
			//Int_t lasthit=Int_t(trinfo(1));
			Float_t lastcharge=(trinfo(2));
			
 			if (lasttrack==idtrack ) {
			    lastcharge+=(Float_t)timeAmplitude;
			    trlist->RemoveAt(lastentry);
			    trinfo(0)=lasttrack;
			    //trinfo(1)=lasthit; // or idhit
			    trinfo(1)=idhit;
			    trinfo(2)=lastcharge;
			    trlist->AddAt(&trinfo,lastentry);
			} else {
			  
		            new((*padr)[countadr++]) TVector(3);
			    TVector &trinfo=*((TVector*) (*padr)[countadr-1]);
			    trinfo(0)=(Float_t)idtrack;
			    trinfo(1)=(Float_t)idhit;
			    trinfo(2)=(Float_t)timeAmplitude;
			  
			    trlist->Add(&trinfo);
			}

#ifdef print
			// check the track list - debugging
                        Int_t trk[20], htrk[20];
                        Float_t chtrk[20];  
			Int_t nptracks=trlist->GetEntriesFast();
			if (nptracks > 2) {
                            Int_t tr;
			    for (tr=0;tr<nptracks;tr++) {
				TVector *pptrkp=(TVector*)trlist->At(tr);
				TVector &pptrk=*pptrkp;
				trk[tr]=Int_t(pptrk(0));
				htrk[tr]=Int_t(pptrk(1));
				chtrk[tr]=(pptrk(2));
                                printf("nptracks %d \n",nptracks);
				// set printings
			    }
			} // end if nptracks
#endif
		    } //  end if pdigit

                    arg[4]=counter;
                    arg[5]=countadr;


}


//____________________________________________

void AliITSsimulationSDD::AddDigit(Int_t i, Int_t j, Int_t signal){
  // Adds a Digit.
    // tag with -1 signals coming from background tracks
    // tag with -2 signals coming from pure electronic noise

    Int_t digits[3], tracks[3], hits[3];
    Float_t phys, charges[3];

    Int_t trk[20], htrk[20];
    Float_t chtrk[20];  

    Bool_t do10to8=fResponse->Do10to8();

    if(do10to8) signal=Convert8to10(signal); 
    AliITSTransientDigit *obj = (AliITSTransientDigit*)fHitMap1->GetHit(i,j);
    digits[0]=i;
    digits[1]=j;
    digits[2]=signal;
    if (!obj) {
        phys=0;
        Int_t k;
        for (k=0;k<3;k++) {
	  tracks[k]=-2;
          charges[k]=0;
          hits[k]=-1;
	}
        fITS->AddSimDigit(1,phys,digits,tracks,hits,charges); 
    } else {
      phys=obj->fPhysics;
      TObjArray* trlist=(TObjArray*)obj->TrackList();
      Int_t nptracks=trlist->GetEntriesFast();

      if (nptracks > 20) {
	 cout<<"Attention - nptracks > 20 "<<nptracks<<endl;
	 nptracks=20;
      }
      Int_t tr;
      for (tr=0;tr<nptracks;tr++) {
	  TVector &pp  =*((TVector*)trlist->At(tr));
	  trk[tr]=Int_t(pp(0));
	  htrk[tr]=Int_t(pp(1));
	  chtrk[tr]=(pp(2));
      }
      if (nptracks > 1) {
	//printf("nptracks > 2  -- %d\n",nptracks);
	  SortTracks(trk,chtrk,htrk,nptracks);
      }
      Int_t i;
      if (nptracks < 3 ) {
	 for (i=0; i<nptracks; i++) {
	     tracks[i]=trk[i];
	     charges[i]=chtrk[i];
	     hits[i]=htrk[i];
	 }
	 for (i=nptracks; i<3; i++) {
	     tracks[i]=-3;
	     hits[i]=-1;
	     charges[i]=0;
	 }
      } else {
	 for (i=0; i<3; i++) {
	     tracks[i]=trk[i];
	     charges[i]=chtrk[i];
	     hits[i]=htrk[i];
	 }
      }

      fITS->AddSimDigit(1,phys,digits,tracks,hits,charges); 
 
    }

}

//____________________________________________

void AliITSsimulationSDD::SortTracks(Int_t *tracks,Float_t *charges,Int_t *hits,Int_t ntr){
  //
  // Sort the list of tracks contributing to a given digit
  // Only the 3 most significant tracks are acctually sorted
  //
  
  //
  //  Loop over signals, only 3 times
  //

  
  Float_t qmax;
  Int_t jmax;
  Int_t idx[3] = {-3,-3,-3};
  Float_t jch[3] = {-3,-3,-3};
  Int_t jtr[3] = {-3,-3,-3};
  Int_t jhit[3] = {-3,-3,-3};
  Int_t i,j,imax;
  
  if (ntr<3) imax=ntr;
  else imax=3;
  for(i=0;i<imax;i++){
    qmax=0;
    jmax=0;
    
    for(j=0;j<ntr;j++){
      
      if((i == 1 && j == idx[i-1] )
	 ||(i == 2 && (j == idx[i-1] || j == idx[i-2]))) continue;
      
      if(charges[j] > qmax) {
	qmax = charges[j];
	jmax=j;
      }       
    } 
    
    if(qmax > 0) {
      idx[i]=jmax;
      jch[i]=charges[jmax]; 
      jtr[i]=tracks[jmax]; 
      jhit[i]=hits[jmax]; 
    }
    
  } 
  
  for(i=0;i<3;i++){
    if (jtr[i] == -3) {
         charges[i]=0;
         tracks[i]=-3;
         hits[i]=-1;
    } else {
         charges[i]=jch[i];
         tracks[i]=jtr[i];
         hits[i]=jhit[i];
    }
  }

}
//____________________________________________
void AliITSsimulationSDD::ChargeToSignal() {
  // add baseline, noise, electronics and ADC saturation effects


  Float_t maxadc = fResponse->MaxAdc();    
  Float_t topValue = fResponse->MagicValue();
  Float_t norm = maxadc/topValue;

  char opt1[20], opt2[20];
  fResponse->ParamOptions(opt1,opt2);
  char *read = strstr(opt1,"file");

  Float_t baseline, noise; 

  if (read) {
      static Bool_t readfile=kTRUE;
      //read baseline and noise from file
      if (readfile) ReadBaseline();
      readfile=kFALSE;
  } else fResponse->GetNoiseParam(noise,baseline);

  Float_t contrib=0;

  TRandom random; 
  Int_t i,k,kk; 

  if(!fDoFFT) {
    for (i=0;i<fNofMaps;i++) {
        if  (read && i<fNofMaps) GetAnodeBaseline(i,baseline,noise);
	for(k=0; k<fScaleSize*fMaxNofSamples; k++) {
	   fInZR[k] = fHitMap2->GetSignal(i,k);
	   contrib = baseline + noise*random.Gaus();
	   fInZR[k] += contrib;
	}
	for(k=0; k<fMaxNofSamples; k++) {
	   Float_t newcont = 0.;
	   Float_t maxcont = 0.;
	   for(kk=0;kk<fScaleSize;kk++) {
	     newcont = fInZR[fScaleSize*k+kk];
	     if(newcont > maxcont) maxcont = newcont;
	   }
	   newcont = maxcont;
	   Double_t signal = newcont*norm;
	   if (signal >= maxadc) signal = maxadc -1;
	   // back to analog: ?
	   //signal /=norm;
	   fHitMap2->SetHit(i,k,signal);
	}  
    } // loop over anodes
    return;
  } // end if DoFFT

  for (i=0;i<fNofMaps;i++) {
      if  (read && i<fNofMaps) GetAnodeBaseline(i,baseline,noise);
      for(k=0; k<fScaleSize*fMaxNofSamples; k++) {
	fInZR[k] = fHitMap2->GetSignal(i,k);
	contrib = baseline + noise*random.Gaus();
	fInZR[k] += contrib;
	fInZI[k] = 0.;
      }
      FastFourierTransform(fElectronics,&fInZR[0],&fInZI[0],1);
      for(k=0; k<fScaleSize*fMaxNofSamples; k++) {
 	Double_t rw = fElectronics->GetTraFunReal(k);
	Double_t iw = fElectronics->GetTraFunImag(k);
  	fOutZR[k] = fInZR[k]*rw - fInZI[k]*iw;
  	fOutZI[k] = fInZR[k]*iw + fInZI[k]*rw;
      }
      FastFourierTransform(fElectronics,&fOutZR[0],&fOutZI[0],-1);
      for(k=0; k<fMaxNofSamples; k++) {
	Float_t newcont = 0.;
	//Float_t totcont = 0.;
	Float_t maxcont = 0.;
	for(kk=0;kk<fScaleSize;kk++) {
	  newcont = fOutZR[fScaleSize*k+kk];
	  if(newcont > maxcont) maxcont = newcont;
	  //	  totcont += (0.25*Out_ZR[4*k+kk]);
	}
	newcont = maxcont;
	Double_t signal = newcont*norm;
	if (signal >= maxadc) signal = maxadc -1;
	// back to analog: ?
	// comment the line below because you want to keep the signal in ADCs
	// convert back to nA in cluster finder
	//signal /=norm;
	fHitMap2->SetHit(i,k,signal);
      }      
  } // loop over anodes
  return;

}

//____________________________________________
void AliITSsimulationSDD::GetAnodeBaseline(Int_t i,Float_t &baseline,
                                           Float_t &noise){
  // Returns the Baseline for a particular anode.
    baseline=fBaseline[i];
    noise=fNoise[i];

}

//____________________________________________
void AliITSsimulationSDD::CompressionParam(Int_t i,Int_t &db,Int_t &tl,
                                           Int_t &th){
  // Returns the compression alogirthm parameters
   Int_t size = fD.GetSize();
   if (size > 2 ) {
      db=fD[i]; tl=fT1[i]; th=fT2[i];
   } else {
      if (size <= 2 && i>=fNofMaps/2) {
	db=fD[1]; tl=fT1[1]; th=fT2[1];
      } else {
        db=fD[0]; tl=fT1[0]; th=fT2[0];
      }
   }
}
//____________________________________________
void AliITSsimulationSDD::CompressionParam(Int_t i,Int_t &db,Int_t &tl){
  // returns the compression alogirthm parameters
   Int_t size = fD.GetSize();
   if (size > 2 ) {
      db=fD[i]; tl=fT1[i];
   } else {
      if (size <= 2 && i>=fNofMaps/2) {
	db=fD[1]; tl=fT1[1]; 
      } else {
        db=fD[0]; tl=fT1[0]; 
      }
   }

}
//____________________________________________
void AliITSsimulationSDD::SetCompressParam(){
  // Sets the compression alogirthm parameters  
   Int_t cp[8],i;
   
   fResponse->GiveCompressParam(cp);
   for (i=0; i<2; i++) {
       fD[i]  =cp[i];
       fT1[i] =cp[i+2];
       fT2[i] =cp[i+4];
       fTol[i]=cp[i+6];
       printf("\n i, fD, fT1, fT2, fTol %d %d %d %d %d\n",
                                      i,fD[i],fT1[i],fT2[i],fTol[i]);
   }
}

//____________________________________________
void AliITSsimulationSDD::ReadBaseline(){
      // read baseline and noise from file - either a .root file and in this
      // case data should be organised in a tree with one entry for each
      // module => reading should be done accordingly
      // or a classic file and do smth. like this:
  //
  // Read baselines and noise for SDD
  //


    Int_t na,pos;
    Float_t bl,n;
    char input[100], base[100], param[100];
    char *filtmp;

    fResponse->Filenames(input,base,param);
    fFileName=base;
//
    filtmp = gSystem->ExpandPathName(fFileName.Data());
    FILE *bline = fopen(filtmp,"r");
    printf("filtmp %s\n",filtmp);
    na = 0;

    if(bline) {
       while(fscanf(bline,"%d %f %f",&pos, &bl, &n) != EOF) {
          if (pos != na+1) {
             Error("ReadBaseline","Anode number not in increasing order!",
                   filtmp);
             exit(1);
	  }
          fBaseline[na]=bl;
          fNoise[na]=n;
          na++;
       }
    } else {
      Error("ReadBaseline"," THE BASELINE FILE %s DOES NOT EXIST !",
	  filtmp);
      exit(1);
    } // end if(bline)
    
    fclose(bline);
    delete [] filtmp;
} 

//____________________________________________
Int_t AliITSsimulationSDD::Convert10to8(Int_t signal) {
  // To the 10 to 8 bit lossive compression.
  // code from Davide C. and Albert W.

   if (signal < 128)  return signal;
   if (signal < 256)  return (128+((signal-128)>>1));
   if (signal < 512)  return (192+((signal-256)>>3));
   if (signal < 1024) return (224+((signal-512)>>4));
   return 0;

}

//____________________________________________
Int_t AliITSsimulationSDD::Convert8to10(Int_t signal) {
  // Undo the lossive 10 to 8 bit compression.
  // code from Davide C. and Albert W.
   if (signal < 0 || signal > 255) {
       printf("<Convert8to10> out of range %d \n",signal);
       return 0;
   }
   
   if (signal < 128) return signal;
   if (signal < 192) {
     if (TMath::Odd(signal)) return (128+((signal-128)<<1));
     else  return (128+((signal-128)<<1)+1);
   }
   if (signal < 224) {
     if (TMath::Odd(signal)) return (256+((signal-192)<<3)+3);
     else  return (256+((signal-192)<<3)+4);
   }
   if (TMath::Odd(signal)) return (512+((signal-224)<<4)+7);
   else  return (512+((signal-224)<<4)+7);
   return 0;

}

//____________________________________________
AliITSMap*   AliITSsimulationSDD::HitMap(Int_t i){
  //Return the correct map.
    return ((i==0)? fHitMap1 : fHitMap2);
}


//____________________________________________
void AliITSsimulationSDD::ZeroSuppression(const char *option) {
  // perform the zero suppresion
  if (strstr(option,"2D")) {
    //Init2D();              // activate if param change module by module
    Compress2D();
  } else if (strstr(option,"1D")) {
    //Init1D();              // activate if param change module by module
    Compress1D();  
  } else StoreAllDigits();  

}

//____________________________________________
void AliITSsimulationSDD::Init2D(){
     // read in and prepare arrays: fD, fT1, fT2
     //                         savemu[nanodes], savesigma[nanodes] 
      // read baseline and noise from file - either a .root file and in this
      // case data should be organised in a tree with one entry for each
      // module => reading should be done accordingly
      // or a classic file and do smth. like this ( code from Davide C. and
      // Albert W.) :
  //
  // Read 2D zero-suppression parameters for SDD
  //

    if (!strstr(fParam,"file")) return;

    Int_t na,pos,tempTh;
    Float_t mu,sigma;
    Float_t *savemu = new Float_t [fNofMaps];
    Float_t *savesigma = new Float_t [fNofMaps];
    char input[100],basel[100],par[100];
    char *filtmp;


    Int_t minval = fResponse->MinVal();

    fResponse->Filenames(input,basel,par);
    fFileName=par;

//
    filtmp = gSystem->ExpandPathName(fFileName.Data());
    FILE *param = fopen(filtmp,"r");
    na = 0;

    if(param) {
       while(fscanf(param,"%d %f %f",&pos, &mu, &sigma) != EOF) {
          if (pos != na+1) {
             Error("Init2D ","Anode number not in increasing order!",
                   filtmp);
             exit(1);
	  }
          savemu[na]=mu;
          savesigma[na]=sigma;
          if ((2.*sigma) < mu) {
              fD[na] = (Int_t)floor(mu - 2.0*sigma + 0.5);
              mu = 2.0 * sigma;
	  } else fD[na] = 0;
          tempTh = (Int_t)floor(mu+2.25*sigma+0.5) - minval;
          if (tempTh < 0) tempTh=0;
          fT1[na] = tempTh;
          tempTh = (Int_t)floor(mu+3.0*sigma+0.5) - minval;
          if (tempTh < 0) tempTh=0;
          fT2[na] = tempTh;
          na++;
       } // end while

    } else {
      Error("Init2D "," THE FILE %s DOES NOT EXIST !",
	  filtmp);
      exit(1);
    } // end if(param)
    
    fclose(param);
    delete [] filtmp;
    delete [] savemu;
    delete [] savesigma;
} 

//____________________________________________
void AliITSsimulationSDD::Compress2D(){
  //
  // simple ITS cluster finder -- online zero-suppression conditions
  // 
  //

    Int_t db,tl,th;  
    Int_t minval = fResponse->MinVal();
    Bool_t write=fResponse->OutputOption();   
    Bool_t do10to8=fResponse->Do10to8();

    Int_t nz, nl, nh, low, i, j; 

    for (i=0; i<fNofMaps; i++) {
        CompressionParam(i,db,tl,th);
        nz=0; 
        nl=0;
        nh=0;
        low=0;
	for (j=0; j<fMaxNofSamples; j++) {
	    Int_t signal=(Int_t)(fHitMap2->GetSignal(i,j));
	    signal -= db; // if baseline eq. is done here
            if (signal <= 0) {nz++; continue;}
	    if ((signal - tl) < minval) low++;
            if ((signal - th) >= minval) {
	        nh++;
		Bool_t cond=kTRUE;
		FindCluster(i,j,signal,minval,cond);
		if (cond && ((TMath::Abs(fHitMap2->GetSignal(i,j-1))-th)>=minval)) {
		  if(do10to8) signal = Convert10to8(signal);
		  AddDigit(i,j,signal);
		}
	    } else if ((signal - tl) >= minval) nl++;
       } // loop time samples
       if (write) TreeB()->Fill(nz,nl,nh,low,i+1);
    } // loop anodes  

      char hname[30];
      if (write) {
	sprintf(hname,"TNtuple%d_%d",fModule,fEvent);
	TreeB()->Write(hname);
	// reset tree
        TreeB()->Reset();
      }

} 

//_____________________________________________________________________________
void  AliITSsimulationSDD::FindCluster(Int_t i,Int_t j,Int_t signal,
                                       Int_t minval,Bool_t &cond){
//
//  Find clusters according to the online 2D zero-suppression algorithm
//

    Bool_t do10to8=fResponse->Do10to8();

    Bool_t high=kFALSE;

    fHitMap2->FlagHit(i,j);
//
//  check the online zero-suppression conditions
//  
    const Int_t maxNeighbours = 4;

    Int_t nn;
    Int_t dbx,tlx,thx;  
    Int_t xList[maxNeighbours], yList[maxNeighbours];
    fSegmentation->Neighbours(i,j,&nn,xList,yList);
    Int_t in,ix,iy,qns;
    for (in=0; in<nn; in++) {
	ix=xList[in];
        iy=yList[in];
        if (fHitMap2->TestHit(ix,iy)==kUnused) {
	   CompressionParam(ix,dbx,tlx,thx);
           Int_t qn = (Int_t)(fHitMap2->GetSignal(ix,iy));
	   qn -= dbx; // if baseline eq. is done here
	   if ((qn-tlx) < minval) {
	      fHitMap2->FlagHit(ix,iy);
	      continue;
	   } else {
              if ((qn - thx) >= minval) high=kTRUE;
              if (cond) {
                 if(do10to8) signal = Convert10to8(signal);
	 	 AddDigit(i,j,signal);
	      }
	      if(do10to8) qns = Convert10to8(qn);
	      else qns=qn;
	      if (!high) AddDigit(ix,iy,qns);
	      cond=kFALSE;
	      if(!high) fHitMap2->FlagHit(ix,iy);
	   }
	} // TestHit
    } // loop over neighbours

}

//____________________________________________
void AliITSsimulationSDD::Init1D(){
  // this is just a copy-paste of input taken from 2D algo
  // Torino people should give input
  //
  // Read 1D zero-suppression parameters for SDD
  //

    if (!strstr(fParam,"file")) return;

    Int_t na,pos,tempTh;
    Float_t mu,sigma;
    Float_t *savemu = new Float_t [fNofMaps];
    Float_t *savesigma = new Float_t [fNofMaps];
    char input[100],basel[100],par[100];
    char *filtmp;


    Int_t minval = fResponse->MinVal();
    fResponse->Filenames(input,basel,par);
    fFileName=par;

//  set first the disable and tol param
    SetCompressParam();
//
    filtmp = gSystem->ExpandPathName(fFileName.Data());
    FILE *param = fopen(filtmp,"r");
    na = 0;

    if (param) {
          fscanf(param,"%d %d %d %d ", &fT2[0], &fT2[1], &fTol[0], &fTol[1]);
          while(fscanf(param,"%d %f %f",&pos, &mu, &sigma) != EOF) {
	       if (pos != na+1) {
		  Error("Init1D ","Anode number not in increasing order!",
                   filtmp);
		  exit(1);
	       }
	       savemu[na]=mu;
	       savesigma[na]=sigma;
	       if ((2.*sigma) < mu) {
		 fD[na] = (Int_t)floor(mu - 2.0*sigma + 0.5);
		 mu = 2.0 * sigma;
	       } else fD[na] = 0;
	       tempTh = (Int_t)floor(mu+2.25*sigma+0.5) - minval;
	       if (tempTh < 0) tempTh=0;
	       fT1[na] = tempTh;
	       na++;
	  } // end while
    } else {
      Error("Init1D "," THE FILE %s DOES NOT EXIST !",
	  filtmp);
      exit(1);
    } // end if(param)
    
    fclose(param);
    delete [] filtmp;
    delete [] savemu;
    delete [] savesigma;



}
 
//____________________________________________
void AliITSsimulationSDD::Compress1D(){
    // 1D zero-suppression algorithm (from Gianluca A.)

    Int_t dis,tol,thres,decr,diff;  

    UChar_t *str=fStream->Stream();
    Int_t counter=0;

    Bool_t do10to8=fResponse->Do10to8();

    Int_t last=0;
    Int_t k,i,j;
    for (k=0; k<2; k++) {
         tol = Tolerance(k);
         dis = Disable(k);  
	 for (i=0; i<fNofMaps/2; i++) {
	     Bool_t firstSignal=kTRUE;
             Int_t idx=i+k*fNofMaps/2;
	     CompressionParam(idx,decr,thres); 
	     for (j=0; j<fMaxNofSamples; j++) {
		 Int_t signal=(Int_t)(fHitMap2->GetSignal(idx,j));
                 signal -= decr;  // if baseline eq.
		 if(do10to8) signal = Convert10to8(signal);
		 if (signal <= thres) {
                     signal=0;
                     diff=128; 
                     last=0; 
                     // write diff in the buffer for HuffT
                     str[counter]=(UChar_t)diff;
                     counter++;
                     continue;
		 }
                 diff=signal-last;
                 if (diff > 127) diff=127;
                 if (diff < -128) diff=-128;
   
		 if (signal < dis) {
		   // tol has changed to 8 possible cases ? - one can write
		   // this if(TMath::Abs(diff)<tol) ... else ...
		    if(TMath::Abs(diff)<tol) diff=0;
		    // or keep it as it was before
		    /*
                    if (tol==1 && (diff >= -2 && diff <= 1)) diff=0;
                    if (tol==2 && (diff >= -4 && diff <= 3)) diff=0;
                    if (tol==3 && (diff >= -16 && diff <= 15)) diff=0;
		    */
                    AddDigit(idx,j,last+diff);
		 } else {
                   AddDigit(idx,j,signal);
		 }
                 
                 diff += 128;
                 // write diff in the buffer used to compute Huffman tables
                 if (firstSignal) str[counter]=(UChar_t)signal;
		 else str[counter]=(UChar_t)diff;
		 counter++;

                 last=signal;
	         firstSignal=kFALSE;
 	     } // loop time samples
	 } // loop anodes  one half of detector 
    }

    // check
    fStream->CheckCount(counter);

    // open file and write out the stream of diff's
   
    static Bool_t open=kTRUE;
    static TFile *outFile;
    Bool_t write = fResponse->OutputOption();
 
    if (write ) {
	if(open) {
	    SetFileName("stream.root");
	    cout<<"filename "<<fFileName<<endl;
	    outFile=new TFile(fFileName,"recreate");
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
	}	    
	open=kFALSE;
	outFile->cd();
        fStream->Write();
    }  // endif write	

     fStream->ClearStream();

     // back to galice.root file

     TTree *fAli=gAlice->TreeK();
     TFile *file = 0;
	    
     if (fAli) file =fAli->GetCurrentFile();
     file->cd();


} 
//____________________________________________
void AliITSsimulationSDD::StoreAllDigits(){
  // if non-zero-suppressed data 

    Bool_t do10to8=fResponse->Do10to8();

    Int_t i, j, digits[3];
    for (i=0; i<fNofMaps; i++) {
        for (j=0; j<fMaxNofSamples; j++) {
             Int_t signal=(Int_t)(fHitMap2->GetSignal(i,j));
	     if(do10to8) signal = Convert10to8(signal);
	     if(do10to8) signal = Convert8to10(signal); 
             digits[0]=i;
             digits[1]=j;
             digits[2]=signal;
	     fITS->AddRealDigit(1,digits);
	}
    }
} 
//____________________________________________

void AliITSsimulationSDD::CreateHistograms(Int_t scale){
  // Creates histograms of maps for debugging

      Int_t i;

      fHis=new TObjArray(fNofMaps);
      for (i=0;i<fNofMaps;i++) {
	   TString sddName("sdd_");
	   Char_t candNum[4];
	   sprintf(candNum,"%d",i+1);
	   sddName.Append(candNum);
	   (*fHis)[i] = new TH1F(sddName.Data(),"SDD maps",
                scale*fMaxNofSamples,0.,(Float_t) scale*fMaxNofSamples);
      }

}
//____________________________________________
void AliITSsimulationSDD::FillHistograms(){
  // fill 1D histograms from map
  if (!fHis) return; 
  
  for( Int_t i=0; i<fNofMaps; i++) {
    TH1F *hist =(TH1F *)fHis->UncheckedAt(i);
    Int_t nsamples = hist->GetNbinsX();
    for( Int_t j=0; j<nsamples; j++) {
      Double_t signal=fHitMap2->GetSignal(i,j);
      hist->Fill((Float_t)j,signal);
    }
  }
}

//____________________________________________

void AliITSsimulationSDD::ResetHistograms(){
    //
    // Reset histograms for this detector
    //
    Int_t i;
    for (i=0;i<fNofMaps;i++ ) {
	if ((*fHis)[i])    ((TH1F*)(*fHis)[i])->Reset();
    }

}


//____________________________________________

TH1F *AliITSsimulationSDD::GetAnode(Int_t wing, Int_t anode) { 
  // Fills a histogram from a give anode.  
  if (!fHis) return 0;

  if(wing <=0 || wing > 2) {
    cout << "Wrong wing number: " << wing << endl;
    return NULL;
  }
  if(anode <=0 || anode > fNofMaps/2) {
    cout << "Wrong anode number: " << anode << endl;
    return NULL;
  }

  Int_t index = (wing-1)*fNofMaps/2 + anode-1;
  return (TH1F*)((*fHis)[index]); 
}

//____________________________________________

void AliITSsimulationSDD::WriteToFile(TFile *hfile) {
  // Writes the histograms to a file 
  if (!fHis) return;

  hfile->cd();
  Int_t i;
  for(i=0; i<fNofMaps; i++)  (*fHis)[i]->Write(); //fAdcs[i]->Write();
  return;
}
//____________________________________________
Float_t AliITSsimulationSDD::GetNoise() {  
  // Returns the noise value

  //Bool_t do10to8=fResponse->Do10to8();
  //noise will always be in the liniar part of the signal

  Int_t decr;
  Int_t threshold=fT1[0];

  char opt1[20], opt2[20];
  fResponse->ParamOptions(opt1,opt2);
  fParam=opt2;
  char *same = strstr(opt1,"same");
  Float_t noise,baseline;
  if (same) {
    fResponse->GetNoiseParam(noise,baseline);
  } else {
     static Bool_t readfile=kTRUE;
     //read baseline and noise from file
     if (readfile) ReadBaseline();
     readfile=kFALSE;
  }

   TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
   if(c2) delete c2->GetPrimitive("noisehist");
   if(c2) delete c2->GetPrimitive("anode");
   else     c2=new TCanvas("c2");
   c2->cd();
   c2->SetFillColor(0);

   TH1F *noisehist = new TH1F("noisehist","noise",100,0.,(float)2*threshold);
   TH1F *anode = new TH1F("anode","Anode Projection",fMaxNofSamples,0.,(float)fMaxNofSamples);
  Int_t i,k;
  for (i=0;i<fNofMaps;i++) {
    CompressionParam(i,decr,threshold); 
    if  (!same) GetAnodeBaseline(i,baseline,noise);
    anode->Reset();
    for (k=0;k<fMaxNofSamples;k++) {
      Float_t signal=(Float_t)fHitMap2->GetSignal(i,k);
      //if (signal <= (float)threshold) noisehist->Fill(signal-baseline);
      if (signal <= (float)threshold) noisehist->Fill(signal);
      anode->Fill((float)k,signal);
    }
    anode->Draw();
    c2->Update();
  }
  TF1 *gnoise = new TF1("gnoise","gaus",0.,threshold);
  noisehist->Fit("gnoise","RQ");
  noisehist->Draw();
  c2->Update();
  Float_t mnoise = gnoise->GetParameter(1);
  cout << "mnoise : " << mnoise << endl;
  Float_t rnoise = gnoise->GetParameter(2);
  cout << "rnoise : " << rnoise << endl;
  delete noisehist;
  return rnoise;
}
