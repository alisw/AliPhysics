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



#include "AliRun.h"
#include "AliITSetfSDD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSHuffman.h"

const Int_t kMaxNeighbours = 4;

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
  // copute b to the e power, where bothe b and e are Int_ts.
  Int_t power = 1,i;
  for(i=0; i<e; i++) power *= b;
  return power;
}

//_____________________________________________

void FastFourierTransform(AliITSetfSDD *alisddetf,Double_t *real,
                          Double_t *imag,Int_t direction) {
  // Do a Fast Fourier Transform

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
  fD.Set(0);
  fT1.Set(0);
  fT2.Set(0);
  fTol.Set(0);
  fInZR = 0;
  fInZI = 0;
  fOutZR = 0;
  fOutZI = 0;

}
//_____________________________________________________________________________
AliITSsimulationSDD::AliITSsimulationSDD(AliITSsimulationSDD &source){
  // Copy constructor to satify Coding roules only.
  if(this==&source) return;
  printf("Not allowed to make a copy of AliITSsimulationSDD "
         "Using default creater instead\n");
  AliITSsimulationSDD();
}
//_____________________________________________________________________________
AliITSsimulationSDD& AliITSsimulationSDD::operator=(AliITSsimulationSDD 
                                                     &source){
  // Copy constructor to satify Coding roules only.
  if(this==&source) return *this;
  printf("Not allowed to make a = with AliITSsimulationSDD "
         "Using default creater instead\n");
  return *this ;
}
//_____________________________________________________________________________

AliITSsimulationSDD::AliITSsimulationSDD(AliITSsegmentation *seg,
                                         AliITSresponse *resp) {
  // Constructor
      fResponse = resp;
      fSegmentation = seg;

      fHitMap2 = new AliITSMapA2(fSegmentation);
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

      fElectronics = new AliITSetfSDD(timeStep);

      Option_t *opt1, *opt2;
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
      Option_t *opt=fResponse->ZeroSuppOption();
        if (strstr(fParam,"file") ) {
	  fD.Set(fNofMaps);
	  fT1.Set(fNofMaps);
          if (strstr(opt,"2D")) {
	    fT2.Set(fNofMaps);
            fTol.Set(0);
            Init2D();       // desactivate if param change module by module
          } else if(strstr(opt,"1D"))  {
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
	if(write && strstr(opt,"2D")) MakeTreeB();

        // call here if baseline does not change by module
        // ReadBaseline();

        fITS = (AliITS*)gAlice->GetModule("ITS");
        Int_t size=fNofMaps*fMaxNofSamples;
	fStream = new AliITSInStream(size); 
	
   fInZR = new Double_t [fMaxNofSamples];
   fInZI = new Double_t [fMaxNofSamples];
   fOutZR = new Double_t [fMaxNofSamples];
   fOutZI = new Double_t [fMaxNofSamples];  

}


//_____________________________________________________________________________

AliITSsimulationSDD::~AliITSsimulationSDD() { 
  // destructor

  delete fHitMap1;
  delete fHitMap2;
  delete fStream;

  fD.Set(0);
  fT1.Set(0);
  fT2.Set(0);
  fTol.Set(0);
  fNoise.Set(0);
  fBaseline.Set(0);

  if (fHis) {
     fHis->Delete(); 
     delete fHis;     
  }  
  
   delete [] fInZR;
   delete [] fInZI;	
	delete [] fOutZR;
	delete [] fOutZI;
	
	delete  fInZR;
   delete  fInZI;  
	delete  fOutZR;
	delete  fOutZI;				 
}
//_____________________________________________________________________________

void AliITSsimulationSDD::DigitiseModule(AliITSmodule *mod,Int_t md,Int_t ev){
  // create maps to build the lists of tracks
  // for each digit

    fModule=md;
    fEvent=ev;

    TObjArray *fHits = mod->GetHits();
    Int_t nhits = fHits->GetEntriesFast();
    if (!nhits) return;


    TObjArray *list=new TObjArray;
    static TClonesArray *padr=0;
    if(!padr) padr=new TClonesArray("TVector",1000);
    Int_t arg[5] = {0,0,0,0,0}; 
    fHitMap1->SetArray(list);


    Int_t NofAnodes=fNofMaps/2;

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

    const Float_t kconv=1000000.;  // GeV->KeV
    Int_t ii;
    for(ii=0; ii<nhits; ii++) {
	AliITShit *hit = (AliITShit*) fHits->At(ii);
	Int_t hitDetector = hit->GetDetector();
	Float_t xL[3];
	hit->GetPositionL(xL[0],xL[1],xL[2]);
	//  cout << "hit local coordinates: " << xL[0] << "," << xL[1] << "," << xL[2] << endl;
	// Deposited energy in keV
	Float_t avpath = 0.;
	Float_t avanod = 0.;
	Float_t depEnergy = kconv*hit->GetIonization();
	AliITShit *hit1 = 0;
	if(depEnergy == 0.) { 
	  ii++;
	  Float_t xL1[3];
	  hit1 = (AliITShit*) fHits->At(ii);
	  hit1->GetPositionL(xL1[0],xL1[1],xL1[2]);
	  //cout << "hit1 local coordinates: " << xL1[0] << "," << xL1[1] << "," << xL1[2] << endl;
	  //cout << "radius1: " << TMath::Sqrt(xL1[0]*xL1[0]+xL1[1]*xL1[1]) << ", azimuth: " << TMath::ATan2(xL1[0],xL1[1]) << endl;
	  avpath = xL1[0];
	  avanod = xL1[2];
	  depEnergy = kconv*hit1->GetIonization();
	}
	Float_t avDrft = xL[0]+avpath;
	Float_t avAnode = xL[2]+avanod;

	if(avpath != 0.) avDrft /= 2.;
	if(avanod != 0.) avAnode /= 2.;

	Float_t driftPath = 10000.*avDrft;
	//printf("sddLength %f avDrft driftPath %f %f\n",sddLength,avDrft, driftPath);
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
	Int_t timeSample = (Int_t) (driftTime/timeStep + 1);
	if(timeSample > fMaxNofSamples) {
	  cout << "Warning: Wrong Time Sample: " << timeSample << endl;
	  continue;
	}

	//   Anode
	Float_t xAnode = 10000.*(avAnode)/anodePitch + NofAnodes/2;  // +1?
	//    Int_t iAnode = 0.5+xAnode; // xAnode?
	if(xAnode*anodePitch > sddWidth || xAnode*anodePitch < 0.) 
             { cout << "Warning: Z = " << xAnode*anodePitch << endl; }
	Int_t iAnode = (Int_t) (1.+xAnode); // xAnode?
	//    cout << "iAnode " << iAnode << endl;
	if(iAnode < 0 || iAnode > NofAnodes) {
	  cout << "Warning: Wrong iAnode: " << iAnode << endl;
	  continue;
	} 


	// work with the idtrack=entry number in the TreeH
	// Int_t idtrack=mod->GetHitTrackIndex(ii);  
        // or store straight away the particle position in the array
	// of particles : 
        Int_t idtrack = hit->GetTrack();

	//  Signal 2d Shape
	Double_t qRef = (Double_t)fResponse->Qref();
	Double_t diffCoeff = (Double_t)fResponse->DiffCoeff();
    
	Double_t gamma = 1. + 0.155*depEnergy/qRef;
	// Squared Sigma along the anodes
	Double_t sigma2A = 2.*diffCoeff*driftTime*gamma;
	Double_t sigmaT  = TMath::Sqrt(sigma2A)/driftSpeed;
    
	// Peak amplitude in nanoAmpere
	Double_t eVpairs = 3.6;
	Double_t amplitude = 160.*depEnergy/(timeStep*eVpairs*2.*acos(-1.)*sigmaT*TMath::Sqrt(sigma2A));
    
	// Spread the charge 
	// Pixel index
	Int_t ja = iAnode;
	Int_t jt = timeSample;
	// Sub-pixel index
	Int_t nsplit = 8;
	nsplit = (nsplit+1)/2*2;
	// Sub-pixel size
	Double_t aStep = anodePitch/nsplit;
	Double_t tStep = timeStep/nsplit;
	// Define SDD window corresponding to the hit
	Int_t anodeWindow = (Int_t) (4.*TMath::Sqrt(sigma2A)/anodePitch + 1);
	Int_t timeWindow = (Int_t) (4.*sigmaT/timeStep + 1);
	Int_t jamin = (ja - anodeWindow/2 - 1)*nsplit + 1;
	Int_t jamax = (ja + anodeWindow/2)*nsplit;
	if(jamin <= 0) jamin = 1;
	if(jamax > NofAnodes*nsplit) jamax = NofAnodes*nsplit;
	Int_t jtmin = (jt - timeWindow/2 - 1)*nsplit + 1;
	Int_t jtmax = (jt + timeWindow/2)*nsplit;
	if(jtmin <= 0) jtmin = 1;
	if(jtmax > fMaxNofSamples*nsplit) jtmax = fMaxNofSamples*nsplit;
	Double_t rlAnode = log(aStep*amplitude);
	// Spread the charge in the anode-time window
        Int_t ka;
	for(ka=jamin; ka <=jamax; ka++) {
	  Int_t ia = (ka-1)/nsplit + 1;
	  if(ia <= 0) { cout << "Warning: ia < 1: " << endl; continue; }
	  if(ia > NofAnodes) ia = NofAnodes;
	  Double_t aExpo = aStep*(ka)-xAnode*anodePitch;
	  Double_t anodeAmplitude = rlAnode - 0.5*aExpo*aExpo/sigma2A;
	  // Protect against overflows
	  if(anodeAmplitude > -87.3)
	    anodeAmplitude = exp(anodeAmplitude);
	  else
	    anodeAmplitude = 0;
	  if(anodeAmplitude) {
	    Double_t rlTime = log(tStep*anodeAmplitude);
            Int_t kt;
	    for(kt=jtmin; kt<=jtmax; kt++) {
	      Int_t it = (kt-1)/nsplit+1;
	      if(it<=0) { cout << "Warning: it < 1: " << endl; continue; } 
	      if(it>fMaxNofSamples) it = fMaxNofSamples;
	      Double_t tExpo = (tStep*(kt)-driftTime)/sigmaT;
	      Double_t timeAmplitude = rlTime - 0.5*tExpo*tExpo;
	      // Protect against overflows
	      if(timeAmplitude > -87.3)
		timeAmplitude = exp(timeAmplitude);
	      else
		timeAmplitude = 0;

	      Int_t index = ((detector+1)%2)*NofAnodes+ia-1;
	      // build the list of digits for this module	
	      arg[0]=index;
	      arg[1]=it;
	      arg[2]=idtrack;
	      ListOfFiredCells(arg,timeAmplitude,list,padr);
	} // loop over time in window 
      } // end if anodeAmplitude
    } // loop over anodes in window
  } // end loop over hits

  Int_t nentries=list->GetEntriesFast();
  // introduce the electronics effects and do zero-suppression if required
  if (nentries) {
    ChargeToSignal(); 

    Option_t *opt=fResponse->ZeroSuppOption();
    ZeroSuppression(opt);
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
                    Int_t it=arg[1];
                    Int_t idtrack=arg[2];
                    Int_t counter=arg[3];
                    Int_t countadr=arg[4];
                   
                    Int_t digits[3];

		    digits[0]=index;
		    digits[1]=it-1;
		    digits[2]=(Int_t)timeAmplitude;
                    Float_t phys;
		    if (idtrack >= 0) phys=(Float_t)timeAmplitude;
		    else phys=0;
                   
                    Double_t charge=timeAmplitude;
		    AliITSTransientDigit* pdigit;
		    // build the list of fired cells and update the info
		    if (!fHitMap1->TestHit(index, it-1)) {
		      
		        new((*padr)[countadr++]) TVector(2);
			TVector &trinfo=*((TVector*) (*padr)[countadr-1]);
			trinfo(0)=(Float_t)idtrack;
			trinfo(1)=(Float_t)timeAmplitude;

			list->AddAtAndExpand(
			    new AliITSTransientDigit(phys,digits),counter);
			
			fHitMap1->SetHit(index, it-1, counter);
			fHitMap2->SetHit(index, it-1, charge);
			counter++;

			pdigit=(AliITSTransientDigit*)list->
                                                      At(list->GetLast());
			// list of tracks
			TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			trlist->Add(&trinfo);

		    } else {
			pdigit=
                         (AliITSTransientDigit*) fHitMap1->GetHit(index, it-1);
                        charge += fHitMap2->GetSignal(index,it-1);
			fHitMap2->SetHit(index, it-1, charge);
			// update charge
			(*pdigit).fSignal=(Int_t)charge;
			(*pdigit).fPhysics+=phys;			
			// update list of tracks
			TObjArray* trlist=(TObjArray*)pdigit->TrackList();
			Int_t lastentry=trlist->GetLast();
			TVector *ptrkp=(TVector*)trlist->At(lastentry);
			TVector &trinfo=*ptrkp;
			Int_t lasttrack=Int_t(trinfo(0));
			Float_t lastcharge=(trinfo(1));
			
 			if (lasttrack==idtrack ) {
			    lastcharge+=(Float_t)timeAmplitude;
			    trlist->RemoveAt(lastentry);
			    trinfo(0)=lasttrack;
			    trinfo(1)=lastcharge;
			    trlist->AddAt(&trinfo,lastentry);
			} else {
			  
		            new((*padr)[countadr++]) TVector(2);
			    TVector &trinfo=*((TVector*) (*padr)[countadr-1]);
			    trinfo(0)=(Float_t)idtrack;
			    trinfo(1)=(Float_t)timeAmplitude;
			  
			    trlist->Add(&trinfo);
			}

#ifdef print
			// check the track list - debugging
                        Int_t trk[50];
                        Float_t chtrk[50];  
			Int_t nptracks=trlist->GetEntriesFast();
			if (nptracks > 2) {
                            Int_t tr;
			    for(tr=0;tr<nptracks;tr++) {
				TVector *pptrkp=(TVector*)trlist->At(tr);
				TVector &pptrk=*pptrkp;
				trk[tr]=Int_t(pptrk(0));
				chtrk[tr]=(pptrk(1));
                                printf("nptracks %d \n",nptracks);
				// set printings
			    }
			} // end if nptracks
#endif
		    } //  end if pdigit

                    arg[3]=counter;
                    arg[4]=countadr;


}


//____________________________________________

void AliITSsimulationSDD::AddDigit(Int_t i, Int_t j, Int_t signal){
  // Adds a Digit.
    // tag with -1 signals coming from background tracks
    // tag with -2 signals coming from pure electronic noise

    Int_t digits[3], tracks[3];
    Float_t phys, charges[3];

    Int_t trk[20];
    Float_t chtrk[20];  

    signal=Convert8to10(signal); // set a flag in case non-ZS are 10-bit
    AliITSTransientDigit *obj = (AliITSTransientDigit*)fHitMap1->GetHit(i,j);
    digits[0]=i;
    digits[1]=j;
    digits[2]=signal;
    // printf("module anode, time, signal %d %d %d %d\n",fModule,i,j,signal);
    if (!obj) {
        phys=0;
        Int_t k;
        for(k=0;k<3;k++) {
	  tracks[k]=-2;
          charges[k]=0;
	}
        fITS->AddDigit(1,phys,digits,tracks,charges); 
    } else {
      phys=obj->fPhysics;
      //printf("AddDigit - test: fCoord1 fCoord2 fSignal %d %d %d i j signal %d %d %d \n",obj->fCoord1,obj->fCoord2,obj->fSignal,i,j,signal);

      TObjArray* trlist=(TObjArray*)obj->TrackList();
      Int_t nptracks=trlist->GetEntriesFast();

      if (nptracks > 20) {
	 cout<<"Attention - nptracks > 20 "<<nptracks<<endl;
	 nptracks=20;
      }
      Int_t tr;
      for(tr=0;tr<nptracks;tr++) {
	  TVector &pp  =*((TVector*)trlist->At(tr));
	  trk[tr]=Int_t(pp(0));
	  chtrk[tr]=(pp(1));
      }
      if (nptracks > 1) {
	  //printf("AddDigit: nptracks %d\n",nptracks);
	  SortTracks(trk,chtrk,nptracks);
      }
      Int_t i;
      if (nptracks < 3 ) {
	 for(i=0; i<nptracks; i++) {
	     tracks[i]=trk[i];
	     charges[i]=chtrk[i];
	 }
	 for(i=nptracks; i<3; i++) {
	     tracks[i]=0;
	     charges[i]=0;
	 }
      } else {
	 for(i=0; i<3; i++) {
	     tracks[i]=trk[i];
	     charges[i]=chtrk[i];
	 }
      }

      fITS->AddDigit(1,phys,digits,tracks,charges); 
 
    }

}

//____________________________________________

void AliITSsimulationSDD::SortTracks(Int_t *tracks,Float_t *charges,Int_t ntr){
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
    }
    
  } 
  
  for(i=0;i<3;i++){
    if (jtr[i] == -3) {
         charges[i]=0;
         tracks[i]=0;
    } else {
         charges[i]=jch[i];
         tracks[i]=jtr[i];
    }
  }

}
//____________________________________________
void AliITSsimulationSDD::ChargeToSignal() {
  // add baseline, noise, electronics and ADC saturation effects

//  Double_t InZR[fMaxNofSamples];
//  Double_t InZI[fMaxNofSamples];
//  Double_t OutZR[fMaxNofSamples];
//  Double_t OutZI[fMaxNofSamples];
  
   


  Float_t maxadc = fResponse->MaxAdc();    
  Float_t TopValue = fResponse->MagicValue();
  Float_t norm = maxadc/TopValue;


  Option_t *opt1, *opt2;
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
  Bool_t first=kTRUE;

  TRandom *random = new TRandom();
    Int_t i,k; 
    for(i=0;i<=fNofMaps;i++) {
    if  (read) GetAnodeBaseline(i,baseline,noise);
    if (!first) FastFourierTransform(fElectronics,&fInZR[0],&fInZI[0],1);
    for(k=0;k<fMaxNofSamples;k++) {
        if (!first) {
	   // analog to digital ?
	   Double_t signal = fOutZR[k]*norm;
           if (signal > maxadc) signal = maxadc;
	   // back to analog: ?
	   signal /=norm;
	   //printf("ChargeToSignal: signal %f\n",signal);
	   fHitMap2->SetHit(i-1,k,signal);

 	   Double_t rw = fElectronics->GetTraFunReal(k);
	   Double_t iw = fElectronics->GetTraFunImag(k);
	   fOutZR[k] = fInZR[k]*rw - fInZI[k]*iw;
	   fOutZI[k] = fInZR[k]*iw + fInZI[k]*rw;
	   if(i+1 < fNofMaps) fInZR[k] = fHitMap2->GetSignal(i+1,k);
	}

        if (first) {
             fInZR[k] = fHitMap2->GetSignal(i,k);
	 }
        fInZI[k] = 0.;
        // add baseline and noise 
        contrib = baseline + noise*random->Gaus();
        fInZR[k] += contrib;

    } // loop over time
    
    if (first) {
         FastFourierTransform(fElectronics,&fInZR[0],&fInZI[0],1);
	 for(k=0; k<fMaxNofSamples; k++) {
	     Double_t rw = fElectronics->GetTraFunReal(k);
	     Double_t iw = fElectronics->GetTraFunImag(k);
	     fOutZR[k] = fInZR[k]*rw - fInZI[k]*iw;
	     fOutZI[k] = fInZR[k]*iw + fInZI[k]*rw;
	     fInZR[k] = fHitMap2->GetSignal(i+1,k);
	     fInZI[k] = 0.;
	     // add baseline and noise 
	     contrib = baseline + noise*random->Gaus();
	     fInZR[k] += contrib;
	  }
    }
    FastFourierTransform(fElectronics,&fOutZR[0],&fOutZI[0],-1);
    first = kFALSE;
  } // loop over anodes

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
   for(i=0; i<2; i++) {
       fD[i]  =cp[i];
       fT1[i] =cp[i+2];
       fT2[i] =cp[i+4];
       fTol[i]=cp[i+6];
       /*
       printf("\n i, fD, fT1, fT2, fTol %d %d %d %d %d\n",
                                      i,fD[i],fT1[i],fT2[i],fTol[i]);
       */
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
    const char *kinput, *base,*kparam;
    char *filtmp;

    fResponse->Filenames(kinput,base,kparam);
    fFileName=base;
//
    filtmp = gSystem->ExpandPathName(fFileName.Data());
    FILE *bline = fopen(filtmp,"r");
    printf("filtmp %s\n",filtmp);
    na = 0;

    if(bline) {
       while(fscanf(bline,"%d %f %f",&pos, &bl, &n) != EOF) {
	 //printf("na, pos, bl, n %d %d %f %f\n",na, pos, bl, n);
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
void AliITSsimulationSDD::ZeroSuppression(Option_t *option) {
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
    const char *kinput,*kbasel,*kpar;
    char *filtmp;


    Int_t minval = fResponse->MinVal();

    fResponse->Filenames(kinput,kbasel,kpar);
    fFileName=kpar;

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
	 delete savemu;
	 delete savesigma;
	 


} 
//____________________________________________
void AliITSsimulationSDD::Compress2D(){
  //
  // simple ITS cluster finder -- online zero-suppression conditions
  // 
  //

  //printf("Compress2D!\n");

    Int_t db,tl,th;  
    Int_t minval = fResponse->MinVal();
    Bool_t write=fResponse->OutputOption();   

    Int_t nz, nl, nh, low, i, j; 

    for(i=0; i<fNofMaps; i++) {
        CompressionParam(i,db,tl,th);
        nz=0; 
        nl=0;
        nh=0;
        low=0;
	for(j=0; j<fMaxNofSamples; j++) {
	    Int_t signal=(Int_t)(fHitMap2->GetSignal(i,j));
	    signal -= db; // if baseline eq. is done here
            if (signal <= 0) {nz++; continue;}
	    if ((signal - tl) < minval) low++;
            if ((signal - th) >= minval) {
	        nh++;
		Bool_t cond=kTRUE;
		//printf("Compress2D : i j %d %d signal %d\n",i,j,signal);
		FindCluster(i,j,signal,minval,cond);
	    } else if ((signal - tl) >= minval) nl++;
       } // loop time samples
       if (write) TreeB()->Fill(nz,nl,nh,low,i+1);
       //if (nz != 256 && low != 256) printf("i, nz, nl, nh  low %d %d %d %d %d\n",i,nz,nl,nh,low);
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
                                       Int_t minval,Bool_t cond){
//
//  Find clusters according to the online 2D zero-suppression algorithm
//

    Bool_t high=kFALSE;

    fHitMap2->FlagHit(i,j);
//
//  check the online zero-suppression conditions
//  
    Int_t nn;
    Int_t dbx,tlx,thx;  
    Int_t Xlist[kMaxNeighbours], Ylist[kMaxNeighbours];
    fSegmentation->Neighbours(i,j,&nn,Xlist,Ylist);
    Int_t in,ix,iy;
    for(in=0; in<nn; in++) {
	ix=Xlist[in];
        iy=Ylist[in];
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
                 signal = Convert10to8(signal);
		 //printf("FindCl -cond : i j %d %d signal %d\n",i,j,signal);
	 	 AddDigit(i,j,signal);
	      }
	      Int_t qns = Convert10to8(qn);
	      //printf("FindCl : i j %d %d qns %d\n",i,j,qns);
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
    Float_t savemu[fNofMaps], savesigma[fNofMaps];
    const char *kinput,*kbasel,*kpar;
    char *filtmp;


    Int_t minval = fResponse->MinVal();
    fResponse->Filenames(kinput,kbasel,kpar);
    fFileName=kpar;

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


}
 
//____________________________________________
void AliITSsimulationSDD::Compress1D(){
    // 1D zero-suppression algorithm (from Gianluca A.)

    Int_t dis,tol,thres,decr,diff;  
    //char *dfile=strstr(fParam,"file");

    UChar_t *str=fStream->Stream();
    Int_t counter=0;

    Int_t last=0,k,i,j;
    for(k=1; k<=2; k++) {
         tol = Tolerance(k-1);
         dis = Disable(k-1);  
	 for(i=0; i<fNofMaps/2; i++) {
	     Bool_t firstSignal=kTRUE;
	     CompressionParam(k*i,decr,thres); 
	     for(j=0; j<fMaxNofSamples; j++) {
		 Int_t signal=(Int_t)(fHitMap2->GetSignal(k*i,j));
                 signal -= decr;  // if baseline eq.
		 signal = Convert10to8(signal);
		 if (signal < thres) {
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
                    if (tol==1 && (diff >= -2 && diff <= 1)) diff=0;
                    if (tol==2 && (diff >= -4 && diff <= 3)) diff=0;
                    if (tol==3 && (diff >= -16 && diff <= 15)) diff=0;
                    AddDigit(k*i,j,last+diff);
		 } else {
                   AddDigit(k*i,j,signal);
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
    static TFile *OutFile;
    Bool_t write = fResponse->OutputOption();
 
    if (write ) {
	if(open) {
	    SetFileName("stream.root");
	    cout<<"filename "<<fFileName<<endl;
	    OutFile=new TFile(fFileName,"recreate");
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
	}	    
	open=kFALSE;
	OutFile->cd();
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

    Int_t digits[3],i,j;

    for(i=0; i<fNofMaps; i++) {
        for(j=0; j<fMaxNofSamples; j++) {
             Int_t signal=(Int_t)(fHitMap2->GetSignal(i,j));
	     signal = Convert10to8(signal);
	     signal = Convert8to10(signal); // ?
             digits[0]=i;
             digits[1]=j;
             digits[2]=signal;
             fITS->AddRealDigit(1,digits);
	}
    }
} 
//____________________________________________

void AliITSsimulationSDD::CreateHistograms(){
  // Creates histograms of maps for debugging

      Int_t i;
      for(i=0;i<fNofMaps;i++) {
	   TString *sddName = new TString("sdd_");
	   Char_t candNum[4];
	   sprintf(candNum,"%d",i+1);
	   sddName->Append(candNum);
	   (*fHis)[i] = new TH1F(sddName->Data(),"SDD maps",
                              fMaxNofSamples,0.,(Float_t) fMaxNofSamples);
	   delete sddName;
      }

}
//____________________________________________

void AliITSsimulationSDD::ResetHistograms(){
    //
    // Reset histograms for this detector
    //
    Int_t i;
    for(i=0;i<fNofMaps;i++ ) {
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
Float_t AliITSsimulationSDD::GetNoise(Float_t threshold) {  
  // Returns the noise value
  if (!fHis) return 0.;

  TH1F *noisehist = new TH1F("noisehist","noise",100,0.,threshold);
  Int_t i,k;
  for(i=0;i<fNofMaps;i++) {
    Int_t nOfBinsA = ((TH1F*)(*fHis)[i])->GetNbinsX();
    for(k=0;k<nOfBinsA;k++) {
      Float_t content = ((TH1F*)(*fHis)[i])->GetBinContent(k+1);
      if (content < threshold) noisehist->Fill(content);
    }
  }
  TF1 *gnoise = new TF1("gnoise","gaus",0.,threshold);
  noisehist->Fit("gnoise","RQ");
  noisehist->Draw();
  Float_t mnoise = gnoise->GetParameter(1);
  cout << "mnoise : " << mnoise << endl;
  Float_t rnoise = gnoise->GetParameter(2);
  cout << "rnoise : " << rnoise << endl;
  delete noisehist;
  return rnoise;
}
void AliITSsimulationSDD::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliITSsimulationSDD.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliITSsimulation::Streamer(R__b);
      R__b >> fITS;
      R__b >> fHitMap1;
      R__b >> fHitMap2;
      R__b >> fStream;
      R__b >> fElectronics;
      fD.Streamer(R__b);
      fT1.Streamer(R__b);
      fT2.Streamer(R__b);
      fTol.Streamer(R__b);
      fBaseline.Streamer(R__b);
      fNoise.Streamer(R__b);
      R__b >> fTreeB;
      //R__b.ReadArray(fParam); // Not to be printed out?
      fFileName.Streamer(R__b);
      R__b >> fNofMaps;
      R__b >> fMaxNofSamples;
      R__b >> fModule;
      R__b >> fEvent;
      R__b >> fHis;
   } else {
      R__b.WriteVersion(AliITSsimulationSDD::IsA());
      AliITSsimulation::Streamer(R__b);
      R__b << fITS;
      R__b << fHitMap1;
      R__b << fHitMap2;
      R__b << fStream;
      R__b << fElectronics;
      fD.Streamer(R__b);
      fT1.Streamer(R__b);
      fT2.Streamer(R__b);
      fTol.Streamer(R__b);
      fBaseline.Streamer(R__b);
      fNoise.Streamer(R__b);
      R__b << fTreeB;
      //R__b.WriteArray(fParam, __COUNTER__); // Not to be printed out?
      fFileName.Streamer(R__b);
      R__b << fNofMaps;
      R__b << fMaxNofSamples;
      R__b << fModule;
      R__b << fEvent;
      R__b << fHis;
   }
}
