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

/* $Id: AliTOFT0v1.cxx,v 1.8 2010/01/19 16:32:20 noferini Exp $ */

//_________________________________________________________________________
// This is a TTask that made the calculation of the Time zero using TOF.
// Description: The algorithm used to calculate the time zero of interaction
// using TOF detector is the following.
// We select in the ESD some "primary" particles - or tracks in the following - 
// that strike the TOF detector (the larger part are pions, kaons or protons). 
// We choose a set of 10 selected tracks, for each track You have the length
// of the track when the TOF is reached, 
// the momentum and the time of flight
// given by the TOF detector.
// Let consider now only one set of 10 tracks (the algorithm is the same for all sets).
// Assuming the (mass) hypothesis that each track can be AUT a pion, AUT a kaon, AUT a proton,
// we consider all the 3 at 10 possible cases. 
// For each track in each (mass) configuration
// (a configuration can be e.g. pion/pion/kaon/proton/pion/proton/kaon/kaon/pion/pion)
// we calculate the time zero (we know in fact the velocity of the track after 
// the assumption about its mass, the time of flight given by the TOF, and the 
// corresponding path travelled till the TOF detector). Then for each mass configuration we have
// 10 time zero and we can calculate the ChiSquare for the current configuration using the 
// weighted mean over all 10 time zero.
// We call the best assignment the mass configuration that gives the minimum value of the ChiSquare. 
// We plot the weighted mean over all 10 time zero for the best assignment, 
// the ChiSquare for the best assignment and the corresponding confidence level.
// The strong assumption is the MC selection of primary particles. It will be introduced
// in the future also some more realistic simulation about this point. 
// Use case:
// root [0] AliTOFT0v1 * tzero = new AliTOFT0v1("galice.root")
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
// root [1] tzero->ExecuteTask()
// root [2] tzero->ExecuteTask("tim")
//             // available parameters:
//             tim - print benchmarking information
//             all - print usefull informations about the number of misidentified tracks 
//                   and a comparison about the true configuration (known from MC) and the best
//                   assignment
// Different Selections for pp and Pb-Pb: Momentum Range, Max Time, # pions 
//-- Author: F. Pierella
//-- Mod By Silvia Arcelli, Francesco Noferini, Barbara Guerzoni
//////////////////////////////////////////////////////////////////////////////

#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliTOFT0v1.h"

ClassImp(AliTOFT0v1)
           
//____________________________________________________________________________ 
AliTOFT0v1::AliTOFT0v1():
  TObject(),
  fLowerMomBound(0.5),
  fUpperMomBound(3),  
  fTimeResolution(0.80e-10), 
  fTimeCorr(0.), 
  fEvent(0x0)
//   fCalib(0x0)
{
  //
  // default constructor
  //

  Init(NULL);
    
}

           
//____________________________________________________________________________ 
AliTOFT0v1::AliTOFT0v1(AliESDEvent* event): 
  TObject(),
  fLowerMomBound(0.5),
  fUpperMomBound(3.0),  
  fTimeResolution(0.80e-10), 
  fTimeCorr(0.), 
  fEvent(event)
//   fCalib(0x0)
{
  //
  // real constructor
  //
  
  Init(event);

}

/* copy-constructor and operator= suppresed 

//____________________________________________________________________________ 
AliTOFT0v1::AliTOFT0v1(const AliTOFT0v1 & tzero):
  TObject(),
  fLowerMomBound(tzero.fLowerMomBound),
  fUpperMomBound(tzero.fUpperMomBound),  
  fTimeResolution(tzero.fTimeResolution), 
  fTimeCorr(tzero.fTimeCorr), 
  fEvent(tzero.fEvent)
//   fCalib(tzero.fCalib)
{
  //
  // copy constructor
  //
    
  fT0SigmaT0def[0]=tzero.fT0SigmaT0def[0];
  fT0SigmaT0def[1]=tzero.fT0SigmaT0def[1];
  fT0SigmaT0def[2]=tzero.fT0SigmaT0def[2];
  fT0SigmaT0def[3]=tzero.fT0SigmaT0def[3];

}

//____________________________________________________________________________ 
AliTOFT0v1& AliTOFT0v1::operator=(const AliTOFT0v1 &tzero)
{
 //
  // assign. operator
  //

  if (this == &tzero)
    return *this;
  
  fLowerMomBound=tzero.fLowerMomBound;
  fUpperMomBound=tzero.fUpperMomBound;  
  fTimeResolution=tzero.fTimeResolution; 
  fTimeCorr=tzero.fTimeCorr; 
  fEvent=tzero.fEvent;
//   fCalib=tzero.fCalib;
  fT0SigmaT0def[0]=tzero.fT0SigmaT0def[0];
  fT0SigmaT0def[1]=tzero.fT0SigmaT0def[1];
  fT0SigmaT0def[2]=tzero.fT0SigmaT0def[2];
  fT0SigmaT0def[3]=tzero.fT0SigmaT0def[3];

  return *this;
}

*/
//____________________________________________________________________________ 
AliTOFT0v1::~AliTOFT0v1()
{
  // dtor
//   fCalib=NULL;
  fEvent=NULL;

}
//____________________________________________________________________________ 

void
AliTOFT0v1::Init(AliESDEvent *event) 
{

  /* 
   * init
   */

  fEvent = event;
  fT0SigmaT0def[0]=0.;
  fT0SigmaT0def[1]=0.6;
  fT0SigmaT0def[2]=0.;
  fT0SigmaT0def[3]=0.;

}

//____________________________________________________________________________ 
void AliTOFT0v1::SetTimeResolution(Double_t timeresolution){
  // Set the TOF time resolution
  fTimeResolution=timeresolution;
}
//____________________________________________________________________________
//____________________________________________________________________________
Double_t * AliTOFT0v1::DefineT0(Option_t *option) 
{ 
  // Caluclate the Event Time using the ESD TOF time

  fT0SigmaT0def[0]=0.;
  fT0SigmaT0def[1]=0.600;
  fT0SigmaT0def[2]=0.;
  fT0SigmaT0def[3]=0.;

 Float_t timeresolutioninns=fTimeResolution*(1.e+9); // convert in [ns]
  
  const Int_t nmaxtracksinset=10;
//   if(strstr(option,"all")){
//     cout << "Selecting primary tracks with momentum between " << fLowerMomBound << " GeV/c and " << fUpperMomBound << " GeV/c" << endl;
//     cout << "Memorandum: 0 means PION | 1 means KAON | 2 means PROTON" << endl;
//   }
  
  Int_t nsets=0;
  Int_t nUsedTracks=0;
  Int_t ngoodsetsSel= 0;
  Float_t t0bestSel[300];
  Float_t eT0bestSel[300];
  Float_t chiSquarebestSel[300];
  Float_t confLevelbestSel[300];
  Float_t t0bestallSel=0.;
  Float_t eT0bestallSel=0.;
  Float_t sumWt0bestallSel=0.;
  Float_t eMeanTzeroPi=0.;
  Float_t meantzeropi=0.;
  Float_t sumAllweightspi=0.;
  Double_t t0def=-999;
  Double_t deltat0def=999;
  Int_t ngoodtrktrulyused=0;
  Int_t ntracksinsetmyCut = 0;

  Int_t ntrk=fEvent->GetNumberOfTracks();
  
  AliESDtrack **tracks=new AliESDtrack*[ntrk];
  Int_t ngoodtrk=0;
  Int_t ngoodtrkt0 =0;
  Float_t mintime =1E6;
  
  // First Track loop, Selection of good tracks

  for (Int_t itrk=0; itrk<ntrk; itrk++) {
    AliESDtrack *t=fEvent->GetTrack(itrk);
    Double_t momOld=t->GetP();
    Double_t mom=momOld-0.0036*momOld;
    if ((t->GetStatus()&AliESDtrack::kTIME)==0) continue;
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) continue;
    Double_t time=t->GetTOFsignal();
    
    time*=1.E-3; // tof given in nanoseconds	   
    if (!(mom<=fUpperMomBound && mom>=fLowerMomBound))continue;
   
    if (!AcceptTrack(t)) continue;

    if(t->GetP() < fLowerMomBound || t->GetIntegratedLength() < 350 || t->GetTOFsignalToT() < 0.000000001)continue; //skip decays
    if(time <= mintime) mintime=time;
    tracks[ngoodtrk]=t;
    ngoodtrk++;
  }
  
  
//    cout << " N. of ESD tracks                    : " << ntrk << endl;
//    cout << " N. of preselected tracks            : " << ngoodtrk << endl;
//    cout << " Minimum tof time in set (in ns)                 : " << mintime << endl;
  
  AliESDtrack **gtracks=new AliESDtrack*[ngoodtrk];
  
  for (Int_t jtrk=0; jtrk< ngoodtrk; jtrk++) {
    AliESDtrack *t=tracks[jtrk];
    Double_t time=t->GetTOFsignal();

    if((time-mintime*1.E3)<50.E3){ // For pp and per 
      gtracks[ngoodtrkt0]=t;
      ngoodtrkt0++;
    }
  }
  

  Int_t nseteq = (ngoodtrkt0-1)/nmaxtracksinset + 1;
  Int_t nmaxtracksinsetCurrent=ngoodtrkt0/nseteq;
  if(nmaxtracksinsetCurrent*nseteq < ngoodtrkt0) nmaxtracksinsetCurrent++;

  if(ngoodtrkt0<2){
//     cout << "less than 2 tracks, skip event " << endl;
    t0def=-999;
    deltat0def=0.600;
    fT0SigmaT0def[0]=t0def;
    fT0SigmaT0def[1]=deltat0def;
    fT0SigmaT0def[2]=ngoodtrkt0;
    fT0SigmaT0def[3]=ngoodtrkt0;
    //goto finish;
  }
  if(ngoodtrkt0>=2){
  // Decide how many tracks in set 
    Int_t ntracksinset = std::min(ngoodtrkt0,nmaxtracksinsetCurrent);
    Int_t nset=1;

    if(ngoodtrkt0>nmaxtracksinsetCurrent) {nset= (Int_t)(ngoodtrkt0/ntracksinset)+1;} 
        
    // Loop over selected sets
    
    if(nset>=1){
      for (Int_t i=0; i< nset; i++) {   
	
	Float_t t0best=999.;
	Float_t eT0best=999.;
	Float_t chisquarebest=99999.;
	Int_t npionbest=0;
	
	Int_t ntracksinsetmy=0;      
	AliESDtrack **tracksT0=new AliESDtrack*[ntracksinset];
	for (Int_t itrk=0; itrk<ntracksinset; itrk++) {
	  Int_t index = itrk+i*ntracksinset;
	  if(index < ngoodtrkt0){
	    AliESDtrack *t=gtracks[index];
	    tracksT0[itrk]=t;
	    ntracksinsetmy++;
	  }
	}
	
	// Analyse it
	
	Int_t   assparticle[nmaxtracksinset];
	Float_t exptof[nmaxtracksinset][3];
	Float_t timeofflight[nmaxtracksinset];
	Float_t momentum[nmaxtracksinset];
	Float_t timezero[nmaxtracksinset];
	Float_t weightedtimezero[nmaxtracksinset];
	Float_t beta[nmaxtracksinset];
	Float_t texp[nmaxtracksinset];
	Float_t dtexp[nmaxtracksinset];
	Float_t sqMomError[nmaxtracksinset];
	Float_t sqTrackError[nmaxtracksinset];
	Float_t massarray[3]={0.13957,0.493677,0.9382723};
	Float_t tracktoflen[nmaxtracksinset];
	Float_t besttimezero[nmaxtracksinset];
	Float_t besttexp[nmaxtracksinset];
	Float_t besttimeofflight[nmaxtracksinset];
	Float_t bestmomentum[nmaxtracksinset];
	Float_t bestchisquare[nmaxtracksinset];
	Float_t bestweightedtimezero[nmaxtracksinset];
	Float_t bestsqTrackError[nmaxtracksinset];
	Int_t imass[nmaxtracksinset];
	
	for (Int_t j=0; j<ntracksinset; j++) {
	  assparticle[j] = 3;
	  timeofflight[j] = 0;
	  momentum[j] = 0;
	  timezero[j] = 0;
	  weightedtimezero[j] = 0;
	  beta[j] = 0;
	  texp[j] = 0;
	  dtexp[j] = 0;
	  sqMomError[j] = 0;
	  sqTrackError[j] = 0;
	  tracktoflen[j] = 0;
	  besttimezero[j] = 0;
	  besttexp[j] = 0;
	  besttimeofflight[j] = 0;
	  bestmomentum[j] = 0;
	  bestchisquare[j] = 0;
	  bestweightedtimezero[j] = 0;
	  bestsqTrackError[j] = 0;
	  imass[j] = 1;
	}
	
	for (Int_t j=0; j<ntracksinsetmy; j++) {
	  AliESDtrack *t=tracksT0[j];
	  Double_t momOld=t->GetP();
	  Double_t mom=momOld-0.0036*momOld;
	  Double_t time=t->GetTOFsignal();
	  
	  time*=1.E-3; // tof given in nanoseconds	   
	  Double_t exptime[10]; t->GetIntegratedTimes(exptime);
	  Double_t toflen=t->GetIntegratedLength();
	  toflen=toflen/100.; // toflen given in m 
	  
	  timeofflight[j]=time;
	  tracktoflen[j]=toflen;
	  exptof[j][0]=exptime[2]*1.E-3+fTimeCorr;// in ns
	  exptof[j][1]=exptime[3]*1.E-3+fTimeCorr;
	  exptof[j][2]=exptime[4]*1.E-3+fTimeCorr;
	  momentum[j]=mom;
	  assparticle[j]=3;
	  
	} //end  for (Int_t j=0; j<ntracksinsetmy; j++) {
	
	for (Int_t itz=0; itz<ntracksinsetmy;itz++) {
	  beta[itz]=momentum[itz]/sqrt(massarray[0]*massarray[0]
				       +momentum[itz]*momentum[itz]);
	  sqMomError[itz]= ((1.-beta[itz]*beta[itz])*0.01)*((1.-beta[itz]*beta[itz])*0.01)*(tracktoflen[itz]/(0.299792*beta[itz]))*(tracktoflen[itz]/(0.299792*beta[itz])); 
	  sqTrackError[itz]=(timeresolutioninns*timeresolutioninns+sqMomError[itz]); //in ns
	  timezero[itz]=exptof[itz][0]-timeofflight[itz];// in ns
	  weightedtimezero[itz]=timezero[itz]/sqTrackError[itz];
	  sumAllweightspi+=1./sqTrackError[itz];
	  meantzeropi+=weightedtimezero[itz];	
	} // end loop for (Int_t itz=0; itz< ntracksinset;itz++)
	
	
	// Then, Combinatorial Algorithm
	
	if(ntracksinsetmy<2 )break;
	
	for (Int_t j=0; j<ntracksinsetmy; j++) {
	  imass[j] = 3;
	}
	
	Int_t ncombinatorial = Int_t(TMath::Power(3,ntracksinsetmy));
	
	// Loop on mass hypotheses
	for (Int_t k=0; k < ncombinatorial;k++) {
	  for (Int_t j=0; j<ntracksinsetmy; j++) {
	    imass[j] = (k % Int_t(TMath::Power(3,ntracksinsetmy-j)))/Int_t(TMath::Power(3,ntracksinsetmy-j-1));
	    texp[j]=exptof[j][imass[j]];
	    dtexp[j]=GetMomError(imass[j], momentum[j], texp[j]);
	  }
	  Float_t sumAllweights=0.;
	  Float_t meantzero=0.;
	  Float_t eMeanTzero=0.;
	  
	  for (Int_t itz=0; itz<ntracksinsetmy;itz++) {
	    sqTrackError[itz]=
	      (timeresolutioninns*
	       timeresolutioninns
	       +dtexp[itz]*dtexp[itz]*1E-6); //in ns2
	    
	    timezero[itz]=texp[itz]-timeofflight[itz];// in ns		    	  
	    
	    weightedtimezero[itz]=timezero[itz]/sqTrackError[itz];
	    sumAllweights+=1./sqTrackError[itz];
	    meantzero+=weightedtimezero[itz];
	    
	  } // end loop for (Int_t itz=0; itz<15;itz++)
	  
	  meantzero=meantzero/sumAllweights; // it is given in [ns]
	  eMeanTzero=sqrt(1./sumAllweights); // it is given in [ns]
	  
	  // calculate chisquare
	  
	  Float_t chisquare=0.;		
	  for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++) {
	    chisquare+=(timezero[icsq]-meantzero)*(timezero[icsq]-meantzero)/sqTrackError[icsq];
	    
	  } // end loop for (Int_t icsq=0; icsq<15;icsq++) 
	  
	  if(chisquare<=chisquarebest){
	    for(Int_t iqsq = 0; iqsq<ntracksinsetmy; iqsq++) {
	      
	      bestsqTrackError[iqsq]=sqTrackError[iqsq]; 
	      besttimezero[iqsq]=timezero[iqsq]; 
	      bestmomentum[iqsq]=momentum[iqsq]; 
	      besttimeofflight[iqsq]=timeofflight[iqsq]; 
	      besttexp[iqsq]=texp[iqsq]; 
	      bestweightedtimezero[iqsq]=weightedtimezero[iqsq]; 
	      bestchisquare[iqsq]=(timezero[iqsq]-meantzero)*(timezero[iqsq]-meantzero)/sqTrackError[iqsq]; 
	    }
	    
	    Int_t npion=0;
	    for (Int_t j=0; j<ntracksinsetmy; j++) {
	      assparticle[j]=imass[j];
	      if(imass[j] == 0) npion++;
	    }
	    npionbest=npion;
	    chisquarebest=chisquare;	      
	    t0best=meantzero;
	    eT0best=eMeanTzero;
	  } // close if(dummychisquare<=chisquare)
	  
	}
	
	Double_t chi2cut[nmaxtracksinset];
	chi2cut[0] = 0;
	chi2cut[1] = 6.6; // corresponding to a C.L. of 0.01
	for (Int_t j=2; j<ntracksinset; j++) {
	  chi2cut[j] = chi2cut[1] * TMath::Sqrt(j*1.);
	}
	
	Double_t chi2singlecut = chi2cut[ntracksinsetmy-1]/ntracksinsetmy + TMath::Abs(chisquarebest-chi2cut[ntracksinsetmy-1])/ntracksinsetmy;
	
//  	printf("tracks removed with a chi2 > %f (chi2total = %f w.r.t. the limit of %f)\n",chi2singlecut,chisquarebest,chi2cut[ntracksinsetmy-1]);
	
	Bool_t kRedoT0 = kFALSE;
        ntracksinsetmyCut = ntracksinsetmy;
	Bool_t usetrack[nmaxtracksinset];
	for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++) {
	  usetrack[icsq] = kTRUE;
	  if((bestchisquare[icsq] > chisquarebest*0.5 && ntracksinsetmy > 2) || (bestchisquare[icsq] > chi2singlecut)){
	    kRedoT0 = kTRUE;
	    ntracksinsetmyCut--;
	    usetrack[icsq] = kFALSE;
	  }
	} // end loop for (Int_t icsq=0; icsq<15;icsq++) 
	
	//	printf("ntrackinsetmy = %i - %i\n",ntracksinsetmy,ntracksinsetmyCut);
	
	// Loop on mass hypotheses Redo
	if(kRedoT0 && ntracksinsetmyCut > 1){
	  //	  printf("Redo T0\n");
	  for (Int_t k=0; k < ncombinatorial;k++) {
	    for (Int_t j=0; j<ntracksinsetmy; j++) {
	      imass[j] = (k % Int_t(TMath::Power(3,ntracksinsetmy-j))) / Int_t(TMath::Power(3,ntracksinsetmy-j-1));
	      texp[j]=exptof[j][imass[j]];
	      dtexp[j]=GetMomError(imass[j], momentum[j], texp[j]);
	    }
	    
	    Float_t sumAllweights=0.;
	    Float_t meantzero=0.;
	    Float_t eMeanTzero=0.;
	    
	    for (Int_t itz=0; itz<ntracksinsetmy;itz++) {
	      if(! usetrack[itz]) continue;
	      sqTrackError[itz]=
		(timeresolutioninns*
		 timeresolutioninns
		 +dtexp[itz]*dtexp[itz]*1E-6); //in ns2
	      
	      timezero[itz]=texp[itz]-timeofflight[itz];// in ns		    	  
	      
	      weightedtimezero[itz]=timezero[itz]/sqTrackError[itz];
	      sumAllweights+=1./sqTrackError[itz];
	      meantzero+=weightedtimezero[itz];
	      
	    } // end loop for (Int_t itz=0; itz<15;itz++)
	    
	    meantzero=meantzero/sumAllweights; // it is given in [ns]
	    eMeanTzero=sqrt(1./sumAllweights); // it is given in [ns]
	    
	    // calculate chisquare
	    
	    Float_t chisquare=0.;		
	    for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++) {
	      if(! usetrack[icsq]) continue;
	      chisquare+=(timezero[icsq]-meantzero)*(timezero[icsq]-meantzero)/sqTrackError[icsq];
	      
	    } // end loop for (Int_t icsq=0; icsq<15;icsq++) 
	    
	    Int_t npion=0;
	    for (Int_t j=0; j<ntracksinsetmy; j++) {
	      assparticle[j]=imass[j];
	      if(imass[j] == 0) npion++;
	    }
	    
	    if(chisquare<=chisquarebest){
	      for(Int_t iqsq = 0; iqsq<ntracksinsetmy; iqsq++) {
		if(! usetrack[iqsq]) continue;
		bestsqTrackError[iqsq]=sqTrackError[iqsq]; 
		besttimezero[iqsq]=timezero[iqsq]; 
		bestmomentum[iqsq]=momentum[iqsq]; 
		besttimeofflight[iqsq]=timeofflight[iqsq]; 
		besttexp[iqsq]=texp[iqsq]; 
		bestweightedtimezero[iqsq]=weightedtimezero[iqsq]; 
		bestchisquare[iqsq]=(timezero[iqsq]-meantzero)*(timezero[iqsq]-meantzero)/sqTrackError[iqsq]; 
	      }
	      
	      npionbest=npion;
	      chisquarebest=chisquare;	      
	      t0best=meantzero;
	      eT0best=eMeanTzero;
	    } // close if(dummychisquare<=chisquare)
	    
	  }
	}
		
	// filling histos
	Float_t confLevel=999;
	
	// Sets with decent chisquares
	
	if(chisquarebest<999.){
	  Double_t dblechisquare=(Double_t)chisquarebest;
	  confLevel=(Float_t)TMath::Prob(dblechisquare,ntracksinsetmyCut-1); 
//  	  cout << " Set Number " << nsets << endl;	
//  	  cout << "Best Assignment, selection " << assparticle[0] << 
//  	    assparticle[1] << assparticle[2] << 
//  	    assparticle[3] << assparticle[4] << 
//  	    assparticle[5] << endl;
//  	  cout << " Chisquare of the set "<< chisquarebest <<endl;
//  	  cout << " C.L. of the set "<< confLevel <<endl;
//  	  cout << " T0 for this set (in ns)  " << t0best << endl;

	  for(Int_t icsq=0; icsq<ntracksinsetmy;icsq++){

	    if(! usetrack[icsq]) continue;
	    
//  	    cout << "Track # " << icsq  << " T0 offsets = " 
//  		 << besttimezero[icsq]-t0best << 
//  	      " track error = "  << bestsqTrackError[icsq]
//  		 << " Chisquare = " << bestchisquare[icsq] 
//  		 << " Momentum  = " << bestmomentum[icsq] 
//  		 << " TOF   = "     << besttimeofflight[icsq] 
//  		 << " TOF tracking  = " << besttexp[icsq]
//  		 << " is used = " << usetrack[icsq] << endl;
	  }
	  
	  // Pick up only those with C.L. >1%
	  //	  if(confLevel>0.01 && ngoodsetsSel<200){
	  if(confLevel>0.01 && ngoodsetsSel<200){
	    chiSquarebestSel[ngoodsetsSel]=chisquarebest;
	    confLevelbestSel[ngoodsetsSel]=confLevel;
	    t0bestSel[ngoodsetsSel]=t0best/eT0best/eT0best;
	    eT0bestSel[ngoodsetsSel]=1./eT0best/eT0best;
	    t0bestallSel += t0best/eT0best/eT0best;
	    sumWt0bestallSel += 1./eT0best/eT0best;
	    ngoodsetsSel++;
	    ngoodtrktrulyused+=ntracksinsetmyCut;	    
	  }
	  else{
	    //	    printf("conflevel = %f -- ngoodsetsSel = %i -- ntrackset = %i\n",confLevel,ngoodsetsSel,ntracksinsetmy);
	  }
	}	
	delete[] tracksT0;
	nsets++;
	
      } // end for the current set
      
      nUsedTracks =  ngoodtrkt0;  
      if(strstr(option,"all")){
	if(sumAllweightspi>0.){
	  meantzeropi=meantzeropi/sumAllweightspi; // it is given in [ns]
	  eMeanTzeroPi=sqrt(1./sumAllweightspi); // it is given in [ns]
	}      
	
	if(sumWt0bestallSel>0){
	  t0bestallSel  = t0bestallSel/sumWt0bestallSel;
	  eT0bestallSel = sqrt(1./sumWt0bestallSel);
	  
	}// end of if(sumWt0bestallSel>0){
	
// 	cout << "T0 all " << t0bestallSel << " +/- " << eT0bestallSel << "Number of tracks used: "<<ngoodtrktrulyused<<endl;
      }
      
      t0def=t0bestallSel;
      deltat0def=eT0bestallSel;
      if ((TMath::Abs(t0bestallSel) < 0.001)&&(TMath::Abs(eT0bestallSel)<0.001)){
	t0def=-999; deltat0def=0.600;
      }
      
      fT0SigmaT0def[0]=t0def;
      fT0SigmaT0def[1]=TMath::Sqrt(deltat0def*deltat0def*ngoodtrktrulyused/(ngoodtrktrulyused-1));
      fT0SigmaT0def[2]=ngoodtrkt0;
      fT0SigmaT0def[3]=ngoodtrktrulyused;
    }
  }
  
  //   if(strstr(option,"tim") || strstr(option,"all")){
  //     cout << "AliTOFT0v1:" << endl ;
  //}
  
  if(fT0SigmaT0def[1] < 0.01) fT0SigmaT0def[1] = 0.6;

  return fT0SigmaT0def;
  }
//__________________________________________________________________
Double_t * AliTOFT0v1::DefineT0(Option_t *option,Float_t pMinCut,Float_t pMaxCut) 
{ 
  // Caluclate the Event Time using the ESD TOF time

  fT0SigmaT0def[0]=0.;
  fT0SigmaT0def[1]=0.600;
  fT0SigmaT0def[2]=0.;
  fT0SigmaT0def[3]=0.;

 Float_t timeresolutioninns=fTimeResolution*(1.e+9); // convert in [ns]
  
  const Int_t nmaxtracksinset=10;
//   if(strstr(option,"all")){
//     cout << "Selecting primary tracks with momentum between " << fLowerMomBound << " GeV/c and " << fUpperMomBound << " GeV/c" << endl;
//     cout << "Memorandum: 0 means PION | 1 means KAON | 2 means PROTON" << endl;
//   }
  
  
  Int_t nsets=0;
  Int_t nUsedTracks=0;
  Int_t ngoodsetsSel= 0;
  Float_t t0bestSel[300];
  Float_t eT0bestSel[300];
  Float_t chiSquarebestSel[300];
  Float_t confLevelbestSel[300];
  Float_t t0bestallSel=0.;
  Float_t eT0bestallSel=0.;
  Float_t sumWt0bestallSel=0.;
  Float_t eMeanTzeroPi=0.;
  Float_t meantzeropi=0.;
  Float_t sumAllweightspi=0.;
  Double_t t0def=-999;
  Double_t deltat0def=999;
  Int_t ngoodtrktrulyused=0;
  Int_t ntracksinsetmyCut = 0;

  Int_t ntrk=fEvent->GetNumberOfTracks();
  
  AliESDtrack **tracks=new AliESDtrack*[ntrk];
  Int_t ngoodtrk=0;
  Int_t ngoodtrkt0 =0;
  Float_t mintime =1E6;
  
  // First Track loop, Selection of good tracks

  for (Int_t itrk=0; itrk<ntrk; itrk++) {
    AliESDtrack *t=fEvent->GetTrack(itrk);
    Double_t momOld=t->GetP();
    Double_t mom=momOld-0.0036*momOld;
    if ((t->GetStatus()&AliESDtrack::kTIME)==0) continue;
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) continue;
    Double_t time=t->GetTOFsignal();
    
    time*=1.E-3; // tof given in nanoseconds	   
    if (!(mom<=fUpperMomBound && mom>=fLowerMomBound))continue;
   
    if (!AcceptTrack(t)) continue;

    if(t->GetP() < fLowerMomBound || t->GetIntegratedLength() < 350 || t->GetTOFsignalToT() < 0.000000001)continue; //skip decays
    if(t->GetP() > pMinCut && t->GetP() < pMaxCut) continue;
    if(time <= mintime) mintime=time;
    tracks[ngoodtrk]=t;
    ngoodtrk++;
  }
  
  
//    cout << " N. of ESD tracks                    : " << ntrk << endl;
//    cout << " N. of preselected tracks            : " << ngoodtrk << endl;
//    cout << " Minimum tof time in set (in ns)                 : " << mintime << endl;
  
  AliESDtrack **gtracks=new AliESDtrack*[ngoodtrk];
  
  for (Int_t jtrk=0; jtrk< ngoodtrk; jtrk++) {
    AliESDtrack *t=tracks[jtrk];
    Double_t time=t->GetTOFsignal();

    if((time-mintime*1.E3)<50.E3){ // For pp and per 
      gtracks[ngoodtrkt0]=t;
      ngoodtrkt0++;
    }
  }
  

  Int_t nseteq = (ngoodtrkt0-1)/nmaxtracksinset + 1;
  Int_t nmaxtracksinsetCurrent=ngoodtrkt0/nseteq;
  if(nmaxtracksinsetCurrent*nseteq < ngoodtrkt0) nmaxtracksinsetCurrent++;

  if(ngoodtrkt0<2){
//     cout << "less than 2 tracks, skip event " << endl;
    t0def=-999;
    deltat0def=0.600;
    fT0SigmaT0def[0]=t0def;
    fT0SigmaT0def[1]=deltat0def;
    fT0SigmaT0def[2]=ngoodtrkt0;
    fT0SigmaT0def[3]=ngoodtrkt0;
    //goto finish;
  }
  if(ngoodtrkt0>=2){
  // Decide how many tracks in set 
    Int_t ntracksinset = std::min(ngoodtrkt0,nmaxtracksinsetCurrent);
    Int_t nset=1;

    if(ngoodtrkt0>nmaxtracksinsetCurrent) {nset= (Int_t)(ngoodtrkt0/ntracksinset)+1;} 
        
    // Loop over selected sets
    
    if(nset>=1){
      for (Int_t i=0; i< nset; i++) {   
	
	Float_t t0best=999.;
	Float_t eT0best=999.;
	Float_t chisquarebest=99999.;
	Int_t npionbest=0;
	
	Int_t ntracksinsetmy=0;      
	AliESDtrack **tracksT0=new AliESDtrack*[ntracksinset];
	for (Int_t itrk=0; itrk<ntracksinset; itrk++) {
	  Int_t index = itrk+i*ntracksinset;
	  if(index < ngoodtrkt0){
	    AliESDtrack *t=gtracks[index];
	    tracksT0[itrk]=t;
	    ntracksinsetmy++;
	  }
	}
	
	// Analyse it
	
	Int_t   assparticle[nmaxtracksinset];
	Float_t exptof[nmaxtracksinset][3];
	Float_t timeofflight[nmaxtracksinset];
	Float_t momentum[nmaxtracksinset];
	Float_t timezero[nmaxtracksinset];
	Float_t weightedtimezero[nmaxtracksinset];
	Float_t beta[nmaxtracksinset];
	Float_t texp[nmaxtracksinset];
	Float_t dtexp[nmaxtracksinset];
	Float_t sqMomError[nmaxtracksinset];
	Float_t sqTrackError[nmaxtracksinset];
	Float_t massarray[3]={0.13957,0.493677,0.9382723};
	Float_t tracktoflen[nmaxtracksinset];
	Float_t besttimezero[nmaxtracksinset];
	Float_t besttexp[nmaxtracksinset];
	Float_t besttimeofflight[nmaxtracksinset];
	Float_t bestmomentum[nmaxtracksinset];
	Float_t bestchisquare[nmaxtracksinset];
	Float_t bestweightedtimezero[nmaxtracksinset];
	Float_t bestsqTrackError[nmaxtracksinset];
	Int_t imass[nmaxtracksinset];
	
	for (Int_t j=0; j<ntracksinset; j++) {
	  assparticle[j] = 3;
	  timeofflight[j] = 0;
	  momentum[j] = 0;
	  timezero[j] = 0;
	  weightedtimezero[j] = 0;
	  beta[j] = 0;
	  texp[j] = 0;
	  dtexp[j] = 0;
	  sqMomError[j] = 0;
	  sqTrackError[j] = 0;
	  tracktoflen[j] = 0;
	  besttimezero[j] = 0;
	  besttexp[j] = 0;
	  besttimeofflight[j] = 0;
	  bestmomentum[j] = 0;
	  bestchisquare[j] = 0;
	  bestweightedtimezero[j] = 0;
	  bestsqTrackError[j] = 0;
	  imass[j] = 1;
	}
	
	for (Int_t j=0; j<ntracksinsetmy; j++) {
	  AliESDtrack *t=tracksT0[j];
	  Double_t momOld=t->GetP();
	  Double_t mom=momOld-0.0036*momOld;
	  Double_t time=t->GetTOFsignal();
	  
	  time*=1.E-3; // tof given in nanoseconds	   
	  Double_t exptime[10]; t->GetIntegratedTimes(exptime);
	  Double_t toflen=t->GetIntegratedLength();
	  toflen=toflen/100.; // toflen given in m 
	  
	  timeofflight[j]=time;
	  tracktoflen[j]=toflen;
	  exptof[j][0]=exptime[2]*1.E-3+fTimeCorr;// in ns
	  exptof[j][1]=exptime[3]*1.E-3+fTimeCorr;
	  exptof[j][2]=exptime[4]*1.E-3+fTimeCorr;
	  momentum[j]=mom;
	  assparticle[j]=3;
	  
	} //end  for (Int_t j=0; j<ntracksinsetmy; j++) {
	
	for (Int_t itz=0; itz<ntracksinsetmy;itz++) {
	  beta[itz]=momentum[itz]/sqrt(massarray[0]*massarray[0]
				       +momentum[itz]*momentum[itz]);
	  sqMomError[itz]= ((1.-beta[itz]*beta[itz])*0.01)*((1.-beta[itz]*beta[itz])*0.01)*(tracktoflen[itz]/(0.299792*beta[itz]))*(tracktoflen[itz]/(0.299792*beta[itz])); 
	  sqTrackError[itz]=(timeresolutioninns*timeresolutioninns+sqMomError[itz]); //in ns
	  timezero[itz]=exptof[itz][0]-timeofflight[itz];// in ns
	  weightedtimezero[itz]=timezero[itz]/sqTrackError[itz];
	  sumAllweightspi+=1./sqTrackError[itz];
	  meantzeropi+=weightedtimezero[itz];	
	} // end loop for (Int_t itz=0; itz< ntracksinset;itz++)
	
	
	// Then, Combinatorial Algorithm
	
	if(ntracksinsetmy<2 )break;
	
	for (Int_t j=0; j<ntracksinsetmy; j++) {
	  imass[j] = 3;
	}
	
	Int_t ncombinatorial = Int_t(TMath::Power(3,ntracksinsetmy));
	
	// Loop on mass hypotheses
	for (Int_t k=0; k < ncombinatorial;k++) {
	  for (Int_t j=0; j<ntracksinsetmy; j++) {
	    imass[j] = (k % Int_t(TMath::Power(3,ntracksinsetmy-j)))/Int_t(TMath::Power(3,ntracksinsetmy-j-1));
	    texp[j]=exptof[j][imass[j]];
	    dtexp[j]=GetMomError(imass[j], momentum[j], texp[j]);
	  }
	  Float_t sumAllweights=0.;
	  Float_t meantzero=0.;
	  Float_t eMeanTzero=0.;
	  
	  for (Int_t itz=0; itz<ntracksinsetmy;itz++) {
	    sqTrackError[itz]=
	      (timeresolutioninns*
	       timeresolutioninns
	       +dtexp[itz]*dtexp[itz]*1E-6); //in ns2
	    
	    timezero[itz]=texp[itz]-timeofflight[itz];// in ns		    	  
	    
	    weightedtimezero[itz]=timezero[itz]/sqTrackError[itz];
	    sumAllweights+=1./sqTrackError[itz];
	    meantzero+=weightedtimezero[itz];
	    
	  } // end loop for (Int_t itz=0; itz<15;itz++)
	  
	  meantzero=meantzero/sumAllweights; // it is given in [ns]
	  eMeanTzero=sqrt(1./sumAllweights); // it is given in [ns]
	  
	  // calculate chisquare
	  
	  Float_t chisquare=0.;		
	  for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++) {
	    chisquare+=(timezero[icsq]-meantzero)*(timezero[icsq]-meantzero)/sqTrackError[icsq];
	    
	  } // end loop for (Int_t icsq=0; icsq<15;icsq++) 
	  
	  if(chisquare<=chisquarebest){
	    for(Int_t iqsq = 0; iqsq<ntracksinsetmy; iqsq++) {
	      
	      bestsqTrackError[iqsq]=sqTrackError[iqsq]; 
	      besttimezero[iqsq]=timezero[iqsq]; 
	      bestmomentum[iqsq]=momentum[iqsq]; 
	      besttimeofflight[iqsq]=timeofflight[iqsq]; 
	      besttexp[iqsq]=texp[iqsq]; 
	      bestweightedtimezero[iqsq]=weightedtimezero[iqsq]; 
	      bestchisquare[iqsq]=(timezero[iqsq]-meantzero)*(timezero[iqsq]-meantzero)/sqTrackError[iqsq]; 
	    }
	    
	    Int_t npion=0;
	    for (Int_t j=0; j<ntracksinsetmy; j++) {
	      assparticle[j]=imass[j];
	      if(imass[j] == 0) npion++;
	    }
	    npionbest=npion;
	    chisquarebest=chisquare;	      
	    t0best=meantzero;
	    eT0best=eMeanTzero;
	  } // close if(dummychisquare<=chisquare)
	  
	}
	
	Double_t chi2cut[nmaxtracksinset];
	chi2cut[0] = 0;
	chi2cut[1] = 6.6; // corresponding to a C.L. of 0.01
	for (Int_t j=2; j<ntracksinset; j++) {
	  chi2cut[j] = chi2cut[1] * TMath::Sqrt(j*1.);
	}
	
	Double_t chi2singlecut = chi2cut[ntracksinsetmy-1]/ntracksinsetmy + TMath::Abs(chisquarebest-chi2cut[ntracksinsetmy-1])/ntracksinsetmy;
	
//  	printf("tracks removed with a chi2 > %f (chi2total = %f w.r.t. the limit of %f)\n",chi2singlecut,chisquarebest,chi2cut[ntracksinsetmy-1]);
	
	Bool_t kRedoT0 = kFALSE;
        ntracksinsetmyCut = ntracksinsetmy;
	Bool_t usetrack[nmaxtracksinset];
	for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++) {
	  usetrack[icsq] = kTRUE;
	  if((bestchisquare[icsq] > chisquarebest*0.5 && ntracksinsetmy > 2) || (bestchisquare[icsq] > chi2singlecut)){
	    kRedoT0 = kTRUE;
	    ntracksinsetmyCut--;
	    usetrack[icsq] = kFALSE;
	  }
	} // end loop for (Int_t icsq=0; icsq<15;icsq++) 
	
	//	printf("ntrackinsetmy = %i - %i\n",ntracksinsetmy,ntracksinsetmyCut);
	
	// Loop on mass hypotheses Redo
	if(kRedoT0 && ntracksinsetmyCut > 1){
	  //	  printf("Redo T0\n");
	  for (Int_t k=0; k < ncombinatorial;k++) {
	    for (Int_t j=0; j<ntracksinsetmy; j++) {
	      imass[j] = (k % Int_t(TMath::Power(3,ntracksinsetmy-j))) / Int_t(TMath::Power(3,ntracksinsetmy-j-1));
	      texp[j]=exptof[j][imass[j]];
	      dtexp[j]=GetMomError(imass[j], momentum[j], texp[j]);
	    }
	    
	    Float_t sumAllweights=0.;
	    Float_t meantzero=0.;
	    Float_t eMeanTzero=0.;
	    
	    for (Int_t itz=0; itz<ntracksinsetmy;itz++) {
	      if(! usetrack[itz]) continue;
	      sqTrackError[itz]=
		(timeresolutioninns*
		 timeresolutioninns
		 +dtexp[itz]*dtexp[itz]*1E-6); //in ns2
	      
	      timezero[itz]=texp[itz]-timeofflight[itz];// in ns		    	  
	      
	      weightedtimezero[itz]=timezero[itz]/sqTrackError[itz];
	      sumAllweights+=1./sqTrackError[itz];
	      meantzero+=weightedtimezero[itz];
	      
	    } // end loop for (Int_t itz=0; itz<15;itz++)
	    
	    meantzero=meantzero/sumAllweights; // it is given in [ns]
	    eMeanTzero=sqrt(1./sumAllweights); // it is given in [ns]
	    
	    // calculate chisquare
	    
	    Float_t chisquare=0.;		
	    for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++) {
	      if(! usetrack[icsq]) continue;
	      chisquare+=(timezero[icsq]-meantzero)*(timezero[icsq]-meantzero)/sqTrackError[icsq];
	      
	    } // end loop for (Int_t icsq=0; icsq<15;icsq++) 
	    
	    Int_t npion=0;
	    for (Int_t j=0; j<ntracksinsetmy; j++) {
	      assparticle[j]=imass[j];
	      if(imass[j] == 0) npion++;
	    }
	    
	    if(chisquare<=chisquarebest){
	      for(Int_t iqsq = 0; iqsq<ntracksinsetmy; iqsq++) {
		if(! usetrack[iqsq]) continue;
		bestsqTrackError[iqsq]=sqTrackError[iqsq]; 
		besttimezero[iqsq]=timezero[iqsq]; 
		bestmomentum[iqsq]=momentum[iqsq]; 
		besttimeofflight[iqsq]=timeofflight[iqsq]; 
		besttexp[iqsq]=texp[iqsq]; 
		bestweightedtimezero[iqsq]=weightedtimezero[iqsq]; 
		bestchisquare[iqsq]=(timezero[iqsq]-meantzero)*(timezero[iqsq]-meantzero)/sqTrackError[iqsq]; 
	      }
	      
	      npionbest=npion;
	      chisquarebest=chisquare;	      
	      t0best=meantzero;
	      eT0best=eMeanTzero;
	    } // close if(dummychisquare<=chisquare)
	    
	  }
	}
		
	// filling histos
	Float_t confLevel=999;
	
	// Sets with decent chisquares
	
	if(chisquarebest<999.){
	  Double_t dblechisquare=(Double_t)chisquarebest;
	  confLevel=(Float_t)TMath::Prob(dblechisquare,ntracksinsetmyCut-1); 
//  	  cout << " Set Number " << nsets << endl;	
//  	  cout << "Best Assignment, selection " << assparticle[0] << 
//  	    assparticle[1] << assparticle[2] << 
//  	    assparticle[3] << assparticle[4] << 
//  	    assparticle[5] << endl;
//  	  cout << " Chisquare of the set "<< chisquarebest <<endl;
//  	  cout << " C.L. of the set "<< confLevel <<endl;
//  	  cout << " T0 for this set (in ns)  " << t0best << endl;

	  for(Int_t icsq=0; icsq<ntracksinsetmy;icsq++){

	    if(! usetrack[icsq]) continue;
	    
//  	    cout << "Track # " << icsq  << " T0 offsets = " 
//  		 << besttimezero[icsq]-t0best << 
//  	      " track error = "  << bestsqTrackError[icsq]
//  		 << " Chisquare = " << bestchisquare[icsq] 
//  		 << " Momentum  = " << bestmomentum[icsq] 
//  		 << " TOF   = "     << besttimeofflight[icsq] 
//  		 << " TOF tracking  = " << besttexp[icsq]
//  		 << " is used = " << usetrack[icsq] << endl;
	  }
	  
	  // Pick up only those with C.L. >1%
	  //	  if(confLevel>0.01 && ngoodsetsSel<200){
	  if(confLevel>0.01 && ngoodsetsSel<200){
	    chiSquarebestSel[ngoodsetsSel]=chisquarebest;
	    confLevelbestSel[ngoodsetsSel]=confLevel;
	    t0bestSel[ngoodsetsSel]=t0best/eT0best/eT0best;
	    eT0bestSel[ngoodsetsSel]=1./eT0best/eT0best;
	    t0bestallSel += t0best/eT0best/eT0best;
	    sumWt0bestallSel += 1./eT0best/eT0best;
	    ngoodsetsSel++;
	    ngoodtrktrulyused+=ntracksinsetmyCut;	    
	  }
	  else{
	    //	    printf("conflevel = %f -- ngoodsetsSel = %i -- ntrackset = %i\n",confLevel,ngoodsetsSel,ntracksinsetmy);
	  }
	}	
	delete[] tracksT0;
	nsets++;
	
      } // end for the current set
      
      nUsedTracks =  ngoodtrkt0;  
      if(strstr(option,"all")){
	if(sumAllweightspi>0.){
	  meantzeropi=meantzeropi/sumAllweightspi; // it is given in [ns]
	  eMeanTzeroPi=sqrt(1./sumAllweightspi); // it is given in [ns]
	}      
	
	if(sumWt0bestallSel>0){
	  t0bestallSel  = t0bestallSel/sumWt0bestallSel;
	  eT0bestallSel = sqrt(1./sumWt0bestallSel);
	  
	}// end of if(sumWt0bestallSel>0){
	
// 	cout << "T0 all " << t0bestallSel << " +/- " << eT0bestallSel << "Number of tracks used: "<<ngoodtrktrulyused<<endl;
      }
      
      t0def=t0bestallSel;
      deltat0def=eT0bestallSel;
      if ((TMath::Abs(t0bestallSel) < 0.001)&&(TMath::Abs(eT0bestallSel)<0.001)){
	t0def=-999; deltat0def=0.600;
      }
      
      fT0SigmaT0def[0]=t0def;
      fT0SigmaT0def[1]=TMath::Sqrt(deltat0def*deltat0def*ngoodtrktrulyused/(ngoodtrktrulyused-1));
      fT0SigmaT0def[2]=ngoodtrkt0;
      fT0SigmaT0def[3]=ngoodtrktrulyused;
    }
  }
  
  //   if(strstr(option,"tim") || strstr(option,"all")){
  //     cout << "AliTOFT0v1:" << endl ;
  //}

  if(fT0SigmaT0def[1] < 0.01) fT0SigmaT0def[1] = 0.6;

  return fT0SigmaT0def;
  }
//__________________________________________________________________
Float_t AliTOFT0v1::GetMomError(Int_t index, Float_t mom, Float_t texp) const
{
  // Take the error extimate for the TOF time in the track reconstruction

  static const Double_t kMasses[]={
    0.000511, 0.105658, 0.139570, 0.493677, 0.938272, 1.875613
  };

  Double_t mass=kMasses[index+2];
  Double_t dpp=0.02;      //mean relative pt resolution;
  //  if(mom > 1) dpp = 0.02*mom;
  Double_t sigma=dpp*texp*1E3/(1.+ mom*mom/(mass*mass));

  sigma =TMath::Sqrt(sigma*sigma);

  return sigma;
}

//__________________________________________________________________
Bool_t AliTOFT0v1::AcceptTrack(AliESDtrack *track)
{

  /* TPC refit */
  if (!(track->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
  /* do not accept kink daughters */
  if (track->GetKinkIndex(0)>0) return kFALSE;
  /* N clusters TPC */
  if (track->GetTPCclusters(0) < 50) return kFALSE;
  /* chi2 TPC */
  if (track->GetTPCchi2()/Float_t(track->GetTPCclusters(0)) > 3.5) return kFALSE;
  /* sigma to vertex */
  if (GetSigmaToVertex(track) > 4.) return kFALSE;
  
  /* accept track */
  return kTRUE;

}

//____________________________________________________________________
Float_t AliTOFT0v1::GetSigmaToVertex(AliESDtrack* esdTrack) const
{
  // Calculates the number of sigma to the vertex.

  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);
  
  if (bCov[0]<=0 || bCov[2]<=0) {
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // -----------------------------------
  // How to get to a n-sigma cut?
  //
  // The accumulated statistics from 0 to d is
  //
  // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
  // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
  //
  // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-d**2)/2)
  // Can this be expressed in a different way?

  if (bRes[0] == 0 || bRes[1] ==0)
    return -1;

  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  // work around precision problem
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  // 1e-15 corresponds to nsigma ~ 7.7
  if (TMath::Exp(-d * d / 2) < 1e-15)
    return 1000;

  Float_t nSigma = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  return nSigma;
}
