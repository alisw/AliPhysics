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

/* $Id: AliTOFT0v2.cxx,v 1.8 2017/02/13 16:32:20 noferini Exp $ */

#include "AliTOFT0v2.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "TBenchmark.h"
#include "AliPID.h"
#include "AliESDpid.h"

ClassImp(AliTOFT0v2)

//____________________________________________________________________________ 
AliTOFT0v2::AliTOFT0v2(AliESDpid *extPID):
  AliTOFT0v1(extPID),
  fCurrentStartTime(0)
{}

//____________________________________________________________________________ 
AliTOFT0v2::AliTOFT0v2(AliESDEvent* event,AliESDpid *extPID): 
  AliTOFT0v1(event,extPID),
  fCurrentStartTime(0)
{}
//____________________________________________________________________________ 
AliTOFT0v2::~AliTOFT0v2()
{}
//____________________________________________________________________________ 
void AliTOFT0v2::Reset(){
  for(Int_t i=0;i < fgkMaxTracks;i++){
    fStartTimeIndex[i] = -1;
    fTrackLabel[i] = -1;
  }
  fCurrentStartTime = 0;
}
//____________________________________________________________________________ 
//____________________________________________________________________________ 
Double_t * AliTOFT0v2::DefineT0(Option_t *option,Float_t pMinCut,Float_t pMaxCut,Int_t isCalibrationMode) 
{ 
  // calibration mode (default=0): 0=no cal, 1=first half sample of tracks, 2=second half sample of tracks

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
  
  Float_t confLevelThr = 0.001;

  Int_t nsets=0;
  Int_t nUsedTracks=0;
  Int_t ngoodsetsSel= 0;
  Float_t t0bestSel[350];
  Float_t eT0bestSel[350];
  Int_t ntracksUsedInSample[350];
  Float_t chiSquarebestSel[350];
  Float_t confLevelbestSel[350];
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
  Int_t ntracksinset=0;

  // First Track loop, Selection of good tracks

  Int_t ntrackmapped = 0;
  Int_t internalMapping[fgkMaxTracks]; // internal mapping
  std::fill_n(internalMapping, fgkMaxTracks, -1);
  Int_t internalMappingSample[350];

  fTracks->Clear();
  for (Int_t itrk=0; itrk<ntrk; itrk++) {
   if (ntrackmapped >= fgkMaxTracks) continue;
    AliESDtrack *t=fEvent->GetTrack(itrk);
    Double_t momOld=t->GetP();
    Double_t mom=momOld-0.0036*momOld;

    Int_t sample = Int_t(mom*1000)%2 + 1;

    if ((t->GetStatus()&AliESDtrack::kTIME)==0) continue;
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) continue;

    Double_t time=t->GetTOFsignal();
    
    time*=1.E-3; // tof given in nanoseconds	   
    if (!(mom<=fUpperMomBound && mom>=fLowerMomBound))continue;

    if (!AcceptTrack(t)) continue;

    if(t->GetIntegratedLength() < 350)continue; //skip decays
    if(t->GetP() > pMinCut && t->GetP() < pMaxCut) continue;

    if(isCalibrationMode && sample!=isCalibrationMode) continue;

    // mapping of TOF tracks
    fTrackLabel[ntrackmapped] = itrk;
    ntrackmapped++;

    //    printf("%i) status=%i\n",ntrackmapped-1,fStartTimeIndex[ntrackmapped-1]);

    if(fStartTimeIndex[ntrackmapped-1] != -1) continue;

    meantime+=time;
    fTracks->AddLast(t);
    internalMapping[ngoodtrk] = ntrackmapped-1;
    ngoodtrk++;
  }

  if(ngoodtrk > 1) meantime /= ngoodtrk;

  if(ngoodtrk>22) nmaxtracksinset = 6;

  fGTracks->Clear();

  for (Int_t jtrk=0; jtrk< fTracks->GetEntries(); jtrk++) {
    AliESDtrack *t=(AliESDtrack*)fTracks->At(jtrk);
    fGTracks->AddLast(t);
    ngoodtrkt0++;
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
    ntracksinset = std::min(ngoodtrkt0,nmaxtracksinsetCurrent);
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
	    fStartTimeIndex[internalMapping[index]] = fCurrentStartTime;
	    ntracksinsetmy++;
	  }
	}

	// Analyse it
	
	Int_t   assparticle[nmaxtracksinsetMax];
	Float_t exptof[nmaxtracksinsetMax][3];
	Float_t momErr[nmaxtracksinsetMax][3];
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
	  Double_t exptime[AliPID::kSPECIESC]; 
	  t->GetIntegratedTimes(exptime,AliPID::kSPECIESC);
	  Double_t toflen=t->GetIntegratedLength();
	  toflen=toflen/100.; // toflen given in m 
	  
	  timeofflight[j]=time;
	  tracktoflen[j]=toflen;
	  exptof[j][0]=exptime[2]*1.E-3+fTimeCorr;// in ns
	  exptof[j][1]=exptime[3]*1.E-3+fTimeCorr;
	  exptof[j][2]=exptime[4]*1.E-3+fTimeCorr;
	  momentum[j]=mom;
	  assparticle[j]=3;

	  // in principle GetMomError only depends on two indices k=imass[j] and j itslef (see blow in the ncombinatorial loop)
	  // so it should be possible to make a lookup in order to speed up the code:
	  if (fOptFlag) {
	    momErr[j][0]=GetMomError(0, momentum[j], exptof[j][0]);
	    momErr[j][1]=GetMomError(1, momentum[j], exptof[j][1]);
	    momErr[j][2]=GetMomError(2, momentum[j], exptof[j][2]);
	  }
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
	
	//	if(ntracksinsetmy<2 )break;
	
	for (Int_t j=0; j<ntracksinsetmy; j++) {
	  imass[j] = 3;
	}

	Int_t ncombinatorial;
	if (fOptFlag) ncombinatorial = fLookupPowerThree[ntracksinsetmy];
	else ncombinatorial = ToCalculatePower(3,ntracksinsetmy);


	// Loop on mass hypotheses
	for (Int_t k=0; k < ncombinatorial;k++) {
	  for (Int_t j=0; j<ntracksinsetmy; j++) {
	    imass[j] = (k % fLookupPowerThree[ntracksinsetmy-j])/fLookupPowerThree[ntracksinsetmy-j-1];
	    texp[j]=exptof[j][imass[j]];
            if (fOptFlag) dtexp[j]=momErr[j][imass[j]]; // see comments above in the initialisation of momErr
 	    else dtexp[j]=GetMomError(imass[j], momentum[j], texp[j]);
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
	
	Double_t chi2singlecut = chi2cut[ntracksinsetmy-1] + TMath::Abs(chisquarebest-chi2cut[ntracksinsetmy-1])/ntracksinsetmy;
	
	//	printf("tracks removed with a chi2 > %f (chi2total = %f w.r.t. the limit of %f)\n",chi2singlecut,chisquarebest,chi2cut[ntracksinsetmy-1]);
	
	Bool_t kRedoT0 = kFALSE;
        ntracksinsetmyCut = ntracksinsetmy;
	Bool_t usetrack[nmaxtracksinsetMax];
	for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++) {
	  usetrack[icsq] = kTRUE;
	  if(TMath::Prob(bestchisquare[icsq],1) < confLevelThr){
	    //  printf("check %f %f %i %f %f- P= %f\n",bestchisquare[icsq],chisquarebest,ntracksinsetmy,chi2singlecut,chi2cut[ntracksinsetmy-1],TMath::Prob(bestchisquare[icsq],ntracksinsetmy-1));
	    //	    getchar();
	      kRedoT0 = kTRUE;
	      ntracksinsetmyCut--;
	      usetrack[icsq] = kFALSE;
	      Int_t index = icsq+i*ntracksinset;
	      if(internalMapping[index] > -1) fStartTimeIndex[internalMapping[index]] = -1;
	      //	      printf("tracks chi2 = %f\n",bestchisquare[icsq]);
	  }
	} // end loop for (Int_t icsq=0; icsq<15;icsq++) 
	
	// Loop on mass hypotheses Redo
	if(kRedoT0 && ntracksinsetmyCut > 1){
	  //	  printf("Redo T0\n");
	  for (Int_t k=0; k < ncombinatorial;k++) {
	    for (Int_t j=0; j<ntracksinsetmy; j++) {
	      imass[j] = (k % fLookupPowerThree[ntracksinsetmy-j]) / fLookupPowerThree[ntracksinsetmy-j-1];
	      texp[j]=exptof[j][imass[j]];
	      if (fOptFlag) dtexp[j]=momErr[j][imass[j]]; // see comments above in the initialisation of momErr
	      else dtexp[j]=GetMomError(imass[j], momentum[j], texp[j]);
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

	Double_t dblechisquare=(Double_t)chisquarebest;
	confLevel=(Float_t)TMath::Prob(dblechisquare,ntracksinsetmyCut-1); 
	
	Int_t ntrackincurrentsel=0;
	for(Int_t icsq=0; icsq<ntracksinsetmy;icsq++){
	  
	  if(! usetrack[icsq]) continue;
	  
	  ntrackincurrentsel++;
	}
	
	//	  printf("%i) CL(Chi2) = %f < 0.01\n",ngoodsetsSel,confLevel);
	
	// Pick up only those with C.L. >confLevelThr
	if(confLevel>confLevelThr && ngoodsetsSel<350){
	  chiSquarebestSel[ngoodsetsSel]=chisquarebest;
	  confLevelbestSel[ngoodsetsSel]=confLevel;
	  ntracksUsedInSample[ngoodsetsSel]=ntracksinsetmyCut;
	  t0bestSel[ngoodsetsSel]=t0best/eT0best/eT0best;
	  eT0bestSel[ngoodsetsSel]=1./eT0best/eT0best;
	  t0bestallSel += t0best/eT0best/eT0best;
	  sumWt0bestallSel += 1./eT0best/eT0best;
	  internalMappingSample[ngoodsetsSel] = i;
	  ngoodsetsSel++;
	  ngoodtrktrulyused+=ntracksinsetmyCut;
	  //	    printf("T0best = %f +/- %f (%i-%i) -- conflevel = %f\n",t0best,eT0best,ntrackincurrentsel,npionbest,confLevel);
	}
	else{
	  // free all tracks in the current sample from start time computation
	  for (Int_t icsq=0; icsq<ntracksinsetmy;icsq++){
	    Int_t index = icsq+i*ntracksinset;
	    if(internalMapping[index] > -1) fStartTimeIndex[internalMapping[index]] = -1;
	  }  
	  //	    printf("conflevel = %f -- ngoodsetsSel = %i -- ntrackset = %i\n",confLevel,ngoodsetsSel,ntracksinsetmy);
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
		  ngoodtrktrulyused -= ntracksUsedInSample[itz];
		  // free all tracks in the current sample from start time computation
		  for (Int_t icsq=0; icsq<ntracksinset;icsq++){
		    Int_t index = icsq+internalMappingSample[itz]*ntracksinset;
		    if(internalMapping[index] > -1) fStartTimeIndex[internalMapping[index]] = -1;
		  }
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

  Int_t nused=0;
  for(Int_t i=0;i < 2000;i++){
    if(fStartTimeIndex[i]==fCurrentStartTime) nused++;
  }
  //  printf("Track flagged = %i\n",nused);
//   bench->Print("t0computation");
//   printf("(%4.1f < p < %4.1f GeV/c) T0-TOF =%9.1f +/- %5.1f ps (n_track = %i)\n\n",pMinCut,pMaxCut,-fT0SigmaT0def[0]*1000,fT0SigmaT0def[1]*1000,Int_t(fT0SigmaT0def[3]));

  delete bench;
  bench=NULL;

  fCurrentStartTime++;

  return fT0SigmaT0def;
  }
