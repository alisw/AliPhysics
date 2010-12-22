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
#include "TBenchmark.h"
#include "AliPID.h"
#include "AliESDpid.h"

ClassImp(AliTOFT0v1)
           
//____________________________________________________________________________ 
AliTOFT0v1::AliTOFT0v1(AliESDpid *extPID):
  TObject(),
  fLowerMomBound(0.5),
  fUpperMomBound(3),  
  fTimeCorr(0.), 
  fEvent(0x0),
  fPIDesd(extPID),
  fTracks(new TObjArray(10)),
  fGTracks(new TObjArray(10)),
  fTracksT0(new TObjArray(10))
{
  //
  // default constructor
  //
  if(AliPID::ParticleMass(0) == 0) new AliPID();

  if(!fPIDesd){
    fPIDesd = new AliESDpid();
  }

  Init(NULL);
    
}

//____________________________________________________________________________ 
AliTOFT0v1::AliTOFT0v1(AliESDEvent* event,AliESDpid *extPID): 
  TObject(),
  fLowerMomBound(0.5),
  fUpperMomBound(3.0),  
  fTimeCorr(0.), 
  fEvent(event),
  fPIDesd(extPID),
  fTracks(new TObjArray(10)),
  fGTracks(new TObjArray(10)),
  fTracksT0(new TObjArray(10))
{
  //
  // real constructor
  //
  if(AliPID::ParticleMass(0) == 0) new AliPID();

  if(!fPIDesd){
    fPIDesd = new AliESDpid();
  }

  Init(event);

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
  fTimeCorr=tzero.fTimeCorr; 
  fEvent=tzero.fEvent;
  fT0SigmaT0def[0]=tzero.fT0SigmaT0def[0];
  fT0SigmaT0def[1]=tzero.fT0SigmaT0def[1];
  fT0SigmaT0def[2]=tzero.fT0SigmaT0def[2];
  fT0SigmaT0def[3]=tzero.fT0SigmaT0def[3];

  fTracks=tzero.fTracks;
  fGTracks=tzero.fGTracks;
  fTracksT0=tzero.fTracksT0;

  for (Int_t ii=0; ii<tzero.fTracks->GetEntries(); ii++)
    fTracks->AddLast(tzero.fTracks->At(ii));

  for (Int_t ii=0; ii<tzero.fGTracks->GetEntries(); ii++)
    fGTracks->AddLast(tzero.fGTracks->At(ii));

  for (Int_t ii=0; ii<tzero.fTracksT0->GetEntries(); ii++)
    fTracksT0->AddLast(tzero.fTracksT0->At(ii));

  return *this;
}

//____________________________________________________________________________ 
AliTOFT0v1::~AliTOFT0v1()
{
  // dtor
  fEvent=NULL;
  
  if (fTracks) {
    fTracks->Clear();
    delete fTracks;
    fTracks=0x0;
  }

  if (fGTracks) {
    fGTracks->Clear();
    delete fGTracks;
    fGTracks=0x0;
  }

  if (fTracksT0) {
    fTracksT0->Clear();
    delete fTracksT0;
    fTracksT0=0x0;
  }

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
Double_t * AliTOFT0v1::DefineT0(Option_t *option,Float_t pMinCut,Float_t pMaxCut) 
{ 
  TBenchmark *bench=new TBenchmark();
  bench->Start("t0computation");

  // Caluclate the Event Time using the ESD TOF time

  fT0SigmaT0def[0]=0.;
  fT0SigmaT0def[1]=0.600;
  fT0SigmaT0def[2]=0.;
  fT0SigmaT0def[3]=0.;
  
  const Int_t nmaxtracksinsetMax=10;
  Int_t nmaxtracksinset=10;
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
  
  Int_t ngoodtrk=0;
  Int_t ngoodtrkt0 =0;
  Float_t meantime =0;
  
  // First Track loop, Selection of good tracks

  fTracks->Clear();
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

    if(t->GetIntegratedLength() < 350)continue; //skip decays
    if(t->GetP() > pMinCut && t->GetP() < pMaxCut) continue;

    meantime+=time;
    fTracks->AddLast(t);
    ngoodtrk++;
  }

  if(ngoodtrk > 1) meantime /= ngoodtrk;

  if(ngoodtrk>22) nmaxtracksinset = 6;

  fGTracks->Clear();
  for (Int_t jtrk=0; jtrk< fTracks->GetEntries(); jtrk++) {
    AliESDtrack *t=(AliESDtrack*)fTracks->At(jtrk);
    //    Double_t time=t->GetTOFsignal();
    //    if((time-meantime*1.E3)<50.E3){ // For pp and per 
    fGTracks->AddLast(t);
    ngoodtrkt0++;
      //    }
  }

  fTracks->Clear();

  Int_t nseteq = (ngoodtrkt0-1)/nmaxtracksinset + 1;
  Int_t nmaxtracksinsetCurrent=ngoodtrkt0/nseteq;
  if(nmaxtracksinsetCurrent*nseteq < ngoodtrkt0) nmaxtracksinsetCurrent++;


  if(ngoodtrkt0<2){
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
	//	printf("Set %i of %i\n",i+1,nset);
	Float_t t0best=999.;
	Float_t eT0best=999.;
	Float_t chisquarebest=99999.;
	Int_t npionbest=0;
	
	fTracksT0->Clear();
	Int_t ntracksinsetmy=0;      
	for (Int_t itrk=0; itrk<ntracksinset; itrk++) {
	  Int_t index = itrk+i*ntracksinset;
	  if(index < fGTracks->GetEntries()){
	    AliESDtrack *t=(AliESDtrack*)fGTracks->At(index);
	    fTracksT0->AddLast(t);
	    ntracksinsetmy++;
	  }
	}

	// Analyse it
	
	Int_t   assparticle[nmaxtracksinsetMax];
	Float_t exptof[nmaxtracksinsetMax][3];
	Float_t timeofflight[nmaxtracksinsetMax];
	Float_t momentum[nmaxtracksinsetMax];
	Float_t timezero[nmaxtracksinsetMax];
	Float_t weightedtimezero[nmaxtracksinsetMax];
	Float_t beta[nmaxtracksinsetMax];
	Float_t texp[nmaxtracksinsetMax];
	Float_t dtexp[nmaxtracksinsetMax];
	Float_t sqMomError[nmaxtracksinsetMax];
	Float_t sqTrackError[nmaxtracksinsetMax];
	Float_t massarray[3]={0.13957,0.493677,0.9382723};
	Float_t tracktoflen[nmaxtracksinsetMax];
	Float_t besttimezero[nmaxtracksinsetMax];
	Float_t besttexp[nmaxtracksinsetMax];
	Float_t besttimeofflight[nmaxtracksinsetMax];
	Float_t bestmomentum[nmaxtracksinsetMax];
	Float_t bestchisquare[nmaxtracksinsetMax];
	Float_t bestweightedtimezero[nmaxtracksinsetMax];
	Float_t bestsqTrackError[nmaxtracksinsetMax];
	Int_t imass[nmaxtracksinsetMax];
	
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
	
	for (Int_t j=0; j<fTracksT0->GetEntries(); j++) {
	  AliESDtrack *t=(AliESDtrack*)fTracksT0->At(j);
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
	  sqTrackError[itz]=sqMomError[itz]; //in ns
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
	
	Int_t ncombinatorial = ToCalculatePower(3,ntracksinsetmy);
	
	// Loop on mass hypotheses
	for (Int_t k=0; k < ncombinatorial;k++) {
	  for (Int_t j=0; j<ntracksinsetmy; j++) {
	    imass[j] = (k % ToCalculatePower(3,ntracksinsetmy-j))/ToCalculatePower(3,ntracksinsetmy-j-1);
	    texp[j]=exptof[j][imass[j]];
	    dtexp[j]=GetMomError(imass[j], momentum[j], texp[j]);
	  }

	  Float_t sumAllweights=0.;
	  Float_t meantzero=0.;
	  Float_t eMeanTzero=0.;
	  
	  for (Int_t itz=0; itz<ntracksinsetmy;itz++) {
	    sqTrackError[itz]=dtexp[itz]*dtexp[itz]*1E-6; //in ns2
	    
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
	
	Double_t chi2cut[nmaxtracksinsetMax];
	chi2cut[0] = 0;
	chi2cut[1] = 6.6; // corresponding to a C.L. of 0.01
	for (Int_t j=2; j<ntracksinset; j++) {
	  chi2cut[j] = chi2cut[1] * TMath::Sqrt(j*1.);
	}
	
	Double_t chi2singlecut = chi2cut[ntracksinsetmy-1]/ntracksinsetmy + TMath::Abs(chisquarebest-chi2cut[ntracksinsetmy-1])/ntracksinsetmy;
	
	//	printf("tracks removed with a chi2 > %f (chi2total = %f w.r.t. the limit of %f)\n",chi2singlecut,chisquarebest,chi2cut[ntracksinsetmy-1]);
	
	Bool_t kRedoT0 = kFALSE;
        ntracksinsetmyCut = ntracksinsetmy;
	Bool_t usetrack[nmaxtracksinsetMax];
	for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++) {
	  usetrack[icsq] = kTRUE;
	  if((bestchisquare[icsq] > chisquarebest*0.5 && ntracksinsetmy > 2) || (bestchisquare[icsq] > chi2singlecut)){
	      kRedoT0 = kTRUE;
	      ntracksinsetmyCut--;
	      usetrack[icsq] = kFALSE;
	      //	      printf("tracks chi2 = %f\n",bestchisquare[icsq]);
	  }
	} // end loop for (Int_t icsq=0; icsq<15;icsq++) 
	
	// Loop on mass hypotheses Redo
	if(kRedoT0 && ntracksinsetmyCut > 1){
	  //	  printf("Redo T0\n");
	  for (Int_t k=0; k < ncombinatorial;k++) {
	    for (Int_t j=0; j<ntracksinsetmy; j++) {
	      imass[j] = (k % ToCalculatePower(3,ntracksinsetmy-j)) / ToCalculatePower(3,ntracksinsetmy-j-1);
	      texp[j]=exptof[j][imass[j]];
	      dtexp[j]=GetMomError(imass[j], momentum[j], texp[j]);
	    }
	    
	    Float_t sumAllweights=0.;
	    Float_t meantzero=0.;
	    Float_t eMeanTzero=0.;
	    
	    for (Int_t itz=0; itz<ntracksinsetmy;itz++) {
	      if(! usetrack[itz]) continue;
	      sqTrackError[itz]=dtexp[itz]*dtexp[itz]*1E-6; //in ns2
	      
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
	    
	    if(chisquare<=chisquarebest && npion>0){
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
	//	printf("Chi2best of the set = %f \n",chisquarebest);
	
	if(chisquarebest<999.){
	  Double_t dblechisquare=(Double_t)chisquarebest;
	  confLevel=(Float_t)TMath::Prob(dblechisquare,ntracksinsetmyCut-1); 

	  Int_t ntrackincurrentsel=0;
	  for(Int_t icsq=0; icsq<ntracksinsetmy;icsq++){

	    if(! usetrack[icsq]) continue;
	    
	    ntrackincurrentsel++;
	  }
	  
	  //	  printf("%i) CL(Chi2) = %f < 0.01\n",ngoodsetsSel,confLevel);

	  // Pick up only those with C.L. >1%
	  if(confLevel>0.01 && ngoodsetsSel<200){
	    chiSquarebestSel[ngoodsetsSel]=chisquarebest;
	    confLevelbestSel[ngoodsetsSel]=confLevel;
	    t0bestSel[ngoodsetsSel]=t0best/eT0best/eT0best;
	    eT0bestSel[ngoodsetsSel]=1./eT0best/eT0best;
	    t0bestallSel += t0best/eT0best/eT0best;
	    sumWt0bestallSel += 1./eT0best/eT0best;
	    ngoodsetsSel++;
	    ngoodtrktrulyused+=ntracksinsetmyCut;
	    //	    printf("T0best = %f +/- %f (%i-%i) -- conflevel = %f\n",t0best,eT0best,ntrackincurrentsel,npionbest,confLevel);
	  }
	  else{
	    //	    printf("conflevel = %f -- ngoodsetsSel = %i -- ntrackset = %i\n",confLevel,ngoodsetsSel,ntracksinsetmy);
	  }
	}	
	fTracksT0->Clear();
	nsets++;
	
      } // end for the current set
      
      //Redo the computation of the best in order to esclude very bad samples
	if(ngoodsetsSel > 1){
	    Double_t t0BestStep1 = t0bestallSel/sumWt0bestallSel;
	    Int_t nsamples=ngoodsetsSel;
	    ngoodsetsSel=0;
	    t0bestallSel=0;
	    sumWt0bestallSel=0;
	    for (Int_t itz=0; itz<nsamples;itz++) {
		if(TMath::Abs(t0bestSel[itz]/eT0bestSel[itz]-t0BestStep1)<1.0){
		    t0bestallSel += t0bestSel[itz];
		    sumWt0bestallSel += eT0bestSel[itz];	      
		    ngoodsetsSel++;   
		    //	      printf("not rejected %f +/- %f\n",t0bestSel[itz]/eT0bestSel[itz],1./TMath::Sqrt(eT0bestSel[itz]));
		}
		else{
		  //	      printf("rejected %f +/- %f\n",t0bestSel[itz]/eT0bestSel[itz],1./TMath::Sqrt(eT0bestSel[itz]));
		}
	    }
	}
	if(ngoodsetsSel < 1){
	    sumWt0bestallSel = 0.0;
	}
      //--------------------------------End recomputation

      nUsedTracks =  ngoodtrkt0;  
      if(strstr(option,"all")){
	if(sumAllweightspi>0.){
	  meantzeropi=meantzeropi/sumAllweightspi; // it is given in [ns]
	  eMeanTzeroPi=sqrt(1./sumAllweightspi); // it is given in [ns]
	}      
	
	//	printf("t0bestallSel = %f -- eT0bestallSel = %f\n",t0bestallSel,sumWt0bestallSel);

	if(sumWt0bestallSel>0){
	  t0bestallSel  = t0bestallSel/sumWt0bestallSel;
	  eT0bestallSel = sqrt(1./sumWt0bestallSel);
	  //	  printf("Final) t0bestallSel = %f -- eT0bestallSel = %f\n",t0bestallSel,eT0bestallSel);	  
	}// end of if(sumWt0bestallSel>0){
	
      }
      
      t0def=t0bestallSel;
      deltat0def=eT0bestallSel;
      
      fT0SigmaT0def[0]=t0def;
      fT0SigmaT0def[1]=TMath::Sqrt(deltat0def*deltat0def);//*ngoodtrktrulyused/(ngoodtrktrulyused-1));
      fT0SigmaT0def[2]=ngoodtrkt0;
      fT0SigmaT0def[3]=ngoodtrktrulyused;
    }
  }

  fGTracks->Clear();

  if(fT0SigmaT0def[1] < 0.00001) fT0SigmaT0def[1] = 0.6;

  bench->Stop("t0computation");

  fT0SigmaT0def[4]=bench->GetRealTime("t0computation");
  fT0SigmaT0def[5]=bench->GetCpuTime("t0computation");

//   bench->Print("t0computation");
//   printf("(%4.1f < p < %4.1f GeV/c) T0-TOF =%9.1f +/- %5.1f ps (n_track = %i)\n\n",pMinCut,pMaxCut,-fT0SigmaT0def[0]*1000,fT0SigmaT0def[1]*1000,Int_t(fT0SigmaT0def[3]));

  delete bench;
  bench=NULL;

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

  Float_t sigma = fPIDesd->GetTOFResponse().GetExpectedSigma(mom,texp,mass);

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

  //Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
  Float_t d = TMath::Sqrt(ToCalculatePower(b[0]/bRes[0],2) + ToCalculatePower(b[1]/bRes[1],2));

  // work around precision problem
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  // 1e-15 corresponds to nsigma ~ 7.7
  if (TMath::Exp(-d * d / 2) < 1e-15)
    return 1000;

  Float_t nSigma = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  return nSigma;
}
//____________________________________________________________________

Bool_t AliTOFT0v1::CheckTPCMatching(AliESDtrack *track,Int_t imass) const{
    Bool_t status = kFALSE;
    
    Double_t exptimes[5];
    track->GetIntegratedTimes(exptimes);
    
    Float_t dedx = track->GetTPCsignal();
    
    Double_t ptpc[3];
    track->GetInnerPxPyPz(ptpc);
    Float_t momtpc=TMath::Sqrt(ptpc[0]*ptpc[0] + ptpc[1]*ptpc[1] + ptpc[2]*ptpc[2]);

    if(imass > 2 || imass < 0) return status;
    Int_t i = imass+2;
    
    AliPID::EParticleType type=AliPID::EParticleType(i);
    
    Float_t dedxExp = fPIDesd->GetTPCResponse().GetExpectedSignal(momtpc,type);
    Float_t resolutionTPC = fPIDesd->GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),type);
	
    if(TMath::Abs(dedx - dedxExp) < 5 * resolutionTPC){
	status = kTRUE;
    }
    
    return status;
}

//____________________________________________________________________
Float_t AliTOFT0v1::ToCalculatePower(Float_t base, Int_t exponent) const
{
  //
  // Returns base^exponent
  //

  Float_t power=1.;

  for (Int_t ii=exponent; ii>0; ii--)
      power=power*base;

  return power;

}
//____________________________________________________________________
Int_t AliTOFT0v1::ToCalculatePower(Int_t base, Int_t exponent) const
{
  //
  // Returns base^exponent
  //

  Int_t power=1;

  for (Int_t ii=exponent; ii>0; ii--)
      power=power*base;

  return power;

}
