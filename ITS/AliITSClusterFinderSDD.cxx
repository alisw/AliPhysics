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

#include <TFile.h>
#include <TMath.h>
#include <math.h>

#include "AliITSClusterFinderSDD.h"
#include "AliITSMapA1.h"
#include "AliITS.h"
#include "AliITSdigit.h"
#include "AliITSRawCluster.h"
#include "AliITSRecPoint.h"
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"
#include "AliRun.h"



ClassImp(AliITSClusterFinderSDD)

//----------------------------------------------------------
AliITSClusterFinderSDD::AliITSClusterFinderSDD
(AliITSsegmentation *seg, AliITSresponse *response, TClonesArray *digits, TClonesArray *recp)   
{
  // constructor

    fSegmentation=seg;
    fResponse=response;
    fDigits=digits;
    fClusters=recp;
    fNclusters= fClusters->GetEntriesFast();
    SetCutAmplitude();
    SetDAnode();
    SetDTime();
    SetMinPeak();
    SetMinNCells();
    SetMaxNCells();
    SetTimeCorr();
    SetMinCharge();
    fMap=new AliITSMapA1(fSegmentation,fDigits,fCutAmplitude);

}

//_____________________________________________________________________________
AliITSClusterFinderSDD::AliITSClusterFinderSDD()
{
  // constructor
    fSegmentation=0;
    fResponse=0;
    fDigits=0;
    fClusters=0;
    fNclusters=0;
    fMap=0;
    SetCutAmplitude();
    SetDAnode();
    SetDTime();
    SetMinPeak();
    SetMinNCells();
    SetMaxNCells();
    SetTimeCorr();
    SetMinCharge();

}

//_____________________________________________________________________________
AliITSClusterFinderSDD::~AliITSClusterFinderSDD()
{
    // destructor

    if(fMap) delete fMap;

}
//__________________________________________________________________________
AliITSClusterFinderSDD::AliITSClusterFinderSDD(const AliITSClusterFinderSDD &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fMap = source.fMap ;
  this->fCutAmplitude = source.fCutAmplitude ;
  this->fDAnode = source.fDAnode ;
  this->fDTime = source.fDTime ;
  this->fTimeCorr = source.fTimeCorr ;
  this->fMinPeak = source.fMinPeak ;
  this->fMinNCells = source.fMinNCells ;
  this->fMaxNCells = source.fMaxNCells ;
  return;
}

//_________________________________________________________________________
AliITSClusterFinderSDD& 
  AliITSClusterFinderSDD::operator=(const AliITSClusterFinderSDD &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fMap = source.fMap ;
  this->fCutAmplitude = source.fCutAmplitude ;
  this->fDAnode = source.fDAnode ;
  this->fDTime = source.fDTime ;
  this->fTimeCorr = source.fTimeCorr ;
  this->fMinPeak = source.fMinPeak ;
  this->fMinNCells = source.fMinNCells ;
  this->fMaxNCells = source.fMaxNCells ;
  return *this;
}


//_____________________________________________________________________________

void AliITSClusterFinderSDD::Find1DClusters()
{
  // find 1D clusters

  AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
  
  // retrieve the parameters 
  Int_t fNofMaps = fSegmentation->Npz();
  Int_t fMaxNofSamples = fSegmentation->Npx();
  Int_t fNofAnodes = fNofMaps/2;
  Int_t dummy=0;
  Float_t fTimeStep = fSegmentation->Dpx(dummy);
  Float_t fSddLength = fSegmentation->Dx();
  Float_t fDriftSpeed = fResponse->DriftSpeed();
  
  Float_t anodePitch = fSegmentation->Dpz(dummy);
  // map the signal
  fMap->SetThreshold(fCutAmplitude);

  fMap->FillMap();
  
  Float_t noise;
  Float_t baseline;
  fResponse->GetNoiseParam(noise,baseline);
  
  Int_t nofFoundClusters = 0;
  Int_t i;
  Float_t **dfadc = new Float_t*[fNofAnodes];
  for(i=0;i<fNofAnodes;i++) dfadc[i] = new Float_t[fMaxNofSamples];
  Float_t fadc = 0.;
  Float_t fadc1 = 0.;
  Float_t fadc2 = 0.;
  Int_t j,k,idx,l,m;
  for(j=0;j<2;j++) {
    for(k=0;k<fNofAnodes;k++) {
      idx = j*fNofAnodes+k;
      // signal (fadc) & derivative (dfadc)
      dfadc[k][255]=0.;
      for(l=0; l<fMaxNofSamples; l++) {
	fadc2=(Float_t)fMap->GetSignal(idx,l);
	if(l>0) fadc1=(Float_t)fMap->GetSignal(idx,l-1);
	if(l>0) dfadc[k][l-1] = fadc2-fadc1;
      } // samples
    } // anodes
    
    for(k=0;k<fNofAnodes;k++) {
      //cout << "Anode: " << k+1 << ", Wing: " << j+1 << endl;
      idx = j*fNofAnodes+k;
      
      Int_t imax = 0;
      Int_t imaxd = 0;
      Int_t it=0;
      while(it <= fMaxNofSamples-3) {
	
	imax = it;
	imaxd = it;
	// maximum of signal	  
	
	Float_t fadcmax = 0.;
	Float_t dfadcmax = 0.;
	Int_t lthrmina = 1;
	Int_t lthrmint = 3;
	
	Int_t lthra = 1;
	Int_t lthrt = 0;
	
	for(m=0;m<20;m++) {
	  Int_t id = it+m;
	  if(id>=fMaxNofSamples) break;
	  fadc=(float)fMap->GetSignal(idx,id);
	  if(fadc > fadcmax) { fadcmax = fadc; imax = id;}
	  if(fadc > (float)fCutAmplitude) { 
	    lthrt++; 
	  }
	  
	  if(dfadc[k][id] > dfadcmax) {
	    dfadcmax = dfadc[k][id];
	    imaxd = id;
	  }
	}
	it = imaxd;
	
	if(fMap->TestHit(idx,imax) == kEmpty) {it++; continue;}
	
	// cluster charge
	Int_t tstart = it-2;
	if(tstart < 0) tstart = 0;
	
	Bool_t ilcl = 0;
	if(lthrt >= lthrmint && lthra >= lthrmina) ilcl = 1;

	if(ilcl) {
	  nofFoundClusters++;
	  
	  Int_t tstop = tstart;
	  Float_t dfadcmin = 10000.;
	  Int_t ij;
	  for(ij=0; ij<20; ij++) {
	    if(tstart+ij > 255) { tstop = 255; break; }
	    fadc=(float)fMap->GetSignal(idx,tstart+ij);
	    if((dfadc[k][tstart+ij] < dfadcmin) && (fadc > fCutAmplitude)) {
	      tstop = tstart+ij+5;
	      if(tstop > 255) tstop = 255;
	      dfadcmin = dfadc[k][it+ij];
	    }
	  }
	  
	  Float_t clusterCharge = 0.;
	  Float_t clusterAnode = k+0.5;
	  Float_t clusterTime = 0.;
	  Float_t clusterMult = 0.;
	  Float_t clusterPeakAmplitude = 0.;
	  Int_t its,peakpos=-1;
	  Float_t n, baseline;
	  fResponse->GetNoiseParam(n,baseline);
	  for(its=tstart; its<=tstop; its++) {
	    fadc=(float)fMap->GetSignal(idx,its);
	    if(fadc>baseline)
	      fadc-=baseline;
	    else
	      fadc=0.;
	    clusterCharge += fadc;
	    // as a matter of fact we should take the peak pos before FFT
	    // to get the list of tracks !!!
	    if(fadc > clusterPeakAmplitude) {
	      clusterPeakAmplitude = fadc;
	      //peakpos=fMap->GetHitIndex(idx,its);
	      Int_t shift=(int)(fTimeCorr/fTimeStep);
	      if(its>shift && its<(fMaxNofSamples-shift)) peakpos=fMap->GetHitIndex(idx,its+shift);
	      else peakpos=fMap->GetHitIndex(idx,its);
	      if(peakpos<0) peakpos=fMap->GetHitIndex(idx,its);
	    }
	    clusterTime += fadc*its;
	    if(fadc > 0) clusterMult++;
	    if(its == tstop) {
	      clusterTime /= (clusterCharge/fTimeStep);   // ns
	      if(clusterTime > fTimeCorr) clusterTime -= fTimeCorr;   // ns
	    }
	  }
	  
	  Float_t clusteranodePath = (clusterAnode - fNofAnodes/2)*anodePitch;
	  Float_t clusterDriftPath = clusterTime*fDriftSpeed;
	  clusterDriftPath = fSddLength-clusterDriftPath;

	  if(clusterCharge <= 0.) break;
	  AliITSRawClusterSDD clust(j+1,clusterAnode,clusterTime,clusterCharge,clusterPeakAmplitude,peakpos,0.,0.,clusterDriftPath,clusteranodePath,clusterMult,0,0,0,0,0,0,0);
	  iTS->AddCluster(1,&clust);
	  it = tstop;
	} // ilcl
	
	it++;
	
      } // while (samples)
    } // anodes
  } // detectors (2)
  
  //fMap->ClearMap(); 
  
  for(i=0;i<fNofAnodes;i++) delete[] dfadc[i];
  delete [] dfadc;
  
  return;
  
}

//_____________________________________________________________________________

void AliITSClusterFinderSDD::Find1DClustersE()
{

	AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");

	// retrieve the parameters 
	Int_t fNofMaps = fSegmentation->Npz();
	Int_t fMaxNofSamples = fSegmentation->Npx();
	Int_t fNofAnodes = fNofMaps/2;
	Int_t dummy=0;
	Float_t fTimeStep = fSegmentation->Dpx( dummy );
	Float_t fSddLength = fSegmentation->Dx();
	Float_t fDriftSpeed = fResponse->DriftSpeed();
	Float_t anodePitch = fSegmentation->Dpz( dummy );
	Float_t n, baseline;
	fResponse->GetNoiseParam( n, baseline );

	// map the signal
	fMap->SetThreshold( fCutAmplitude );
	fMap->FillMap();
	Int_t nClu = 0;
	
//	cout << "Search  cluster... "<< endl;
 	for( Int_t j=0; j<2; j++ ) 
	{
    	for( Int_t k=0; k<fNofAnodes; k++ ) 
		{
      		Int_t idx = j*fNofAnodes+k;
			
			Bool_t On = kFALSE;
			Int_t start = 0;
			Int_t nTsteps = 0;
			Float_t fmax = 0.;
			Int_t lmax = 0;
			Float_t charge = 0.;
			Float_t time = 0.;
			Float_t anode = k+0.5;
			Int_t peakpos = -1;

      		for( Int_t l=0; l<fMaxNofSamples; l++ ) 
			{
				Float_t fadc = (Float_t)fMap->GetSignal( idx, l );
				if( fadc > 0.0 )
				{
					if( On == kFALSE && l<fMaxNofSamples-4)  // star RawCluster (reset var.)
					{ 
						Float_t fadc1 = (Float_t)fMap->GetSignal( idx, l+1 );
						if( fadc1 < fadc ) continue;
						start = l;
						fmax = 0.;
						lmax = 0;
						time = 0.;
						charge = 0.; 
						On = kTRUE; 
						nTsteps = 0;
					}
					
					nTsteps++ ;
	  				if( fadc > baseline )
	      				fadc -= baseline;
	    			else
	      				fadc=0.;

	    			charge += fadc;
					time += fadc*l;

					if( fadc > fmax ) 
					{ 
						fmax = fadc; 
						lmax = l; 
						Int_t shift = (Int_t)(fTimeCorr/fTimeStep + 0.5);
	                    if( l > shift && l < (fMaxNofSamples-shift) )  
							peakpos = fMap->GetHitIndex( idx, l+shift );
	      				else 
							peakpos = fMap->GetHitIndex( idx, l );
	      				if( peakpos < 0 ) peakpos = fMap->GetHitIndex( idx, l );
					}
				}
				else
				{
					if( On == kTRUE )
					{	
						if( nTsteps > 2 ) //  min # of timesteps for a RawCluster
						{
							// Found a RawCluster...
							Int_t stop = l-1;
							time /= (charge/fTimeStep);   // ns
				//			time = lmax*fTimeStep;   // ns
							if( time > fTimeCorr ) time -= fTimeCorr;   // ns
							Float_t anodePath = (anode - fNofAnodes/2)*anodePitch;
	  						Float_t DriftPath = time*fDriftSpeed;
	  						DriftPath = fSddLength-DriftPath;
					  		AliITSRawClusterSDD clust( j+1, anode, time, charge,
	                            fmax, peakpos, 0., 0., DriftPath, anodePath, nTsteps
				    			, start, stop, start, stop, 1, k, k );
	  						iTS->AddCluster( 1, &clust );
						//	clust.PrintInfo();
							nClu++;
						}
						On = kFALSE;
					}
				}
      		} // samples
    	} // anodes
	} // wings

//	cout << "# Rawclusters " << nClu << endl; 	
	return; 

}

//_____________________________________________________________________________

Int_t AliITSClusterFinderSDD::SearchPeak( Float_t *spect, Int_t xdim, Int_t zdim, 
										  Int_t *peakX, Int_t *peakZ, Float_t *peakAmp, Float_t minpeak )
{
	Int_t npeak = 0;    // # peaks
    Int_t i,j;

	// search peaks
	for( Int_t z=1; z<zdim-1; z++ )
	{
		for( Int_t x=2; x<xdim-3; x++ )
		{
			Float_t Sxz = spect[x*zdim+z];
			Float_t Sxz1 = spect[(x+1)*zdim+z];
			Float_t Sxz2 = spect[(x-1)*zdim+z];

			// search a local max. in s[x,z]
			if( Sxz < minpeak || Sxz1 <= 0 || Sxz2 <= 0 ) continue;
			if( Sxz >= spect[(x+1)*zdim+z  ] && Sxz >= spect[(x-1)*zdim+z  ] &&
				Sxz >= spect[x*zdim    +z+1] && Sxz >= spect[x*zdim    +z-1] &&
				Sxz >= spect[(x+1)*zdim+z+1] && Sxz >= spect[(x+1)*zdim+z-1] &&
				Sxz >= spect[(x-1)*zdim+z+1] && Sxz >= spect[(x-1)*zdim+z-1] )
			{
				// peak found
				peakX[npeak] = x;
				peakZ[npeak] = z;
				peakAmp[npeak] = Sxz;
				npeak++;
			}
		}
	}			

	// search groups of peaks with same amplitude.
	Int_t *Flag = new Int_t[npeak];
	for( i=0; i<npeak; i++ ) Flag[i] = 0;
	for( i=0; i<npeak; i++ )
	{
		for( j=0; j<npeak; j++ )
		{
			if( i==j) continue;
			if( Flag[j] > 0 ) continue;
			if( peakAmp[i] == peakAmp[j] && abs(peakX[i]-peakX[j])<=1 && abs(peakZ[i]-peakZ[j])<=1 )
			{
				if( Flag[i] == 0) Flag[i] = i+1;
				Flag[j] = Flag[i];
			}
		}
	}

	// make average of peak groups	
	for( i=0; i<npeak; i++ )
	{
		Int_t nFlag = 1;
		if( Flag[i] <= 0 ) continue;
		for( j=0; j<npeak; j++ )
		{
			if( i==j ) continue;
			if( Flag[j] != Flag[i] ) continue;
			peakX[i] += peakX[j];
			peakZ[i] += peakZ[j];
			nFlag++;
			npeak--;
			for( Int_t k=j; k<npeak; k++ )
			{
				peakX[k] = peakX[k+1];
				peakZ[k] = peakZ[k+1];
				peakAmp[k] = peakAmp[k+1];
				Flag[k] = Flag[k+1];
			}	
			j--;
		}

		if( nFlag > 1 )
		{
			peakX[i] /= nFlag;
			peakZ[i] /= nFlag;
		}
	}
	
	delete [] Flag;
	return( npeak );
}


void AliITSClusterFinderSDD::PeakFunc( Int_t xdim, Int_t zdim, Float_t *par, Float_t *spe, Float_t *Integral) 
{
  Int_t Electronics = fResponse->Electronics(); // 1 = PASCAL, 2 = OLA
  Int_t param_peak = 5;
  // par -> paramiters..
  // par[0]  number of peaks.
  // for each peak i=1, ..., par[0]
  // 		par[i] = Ampl.
  // 		par[i+1] = xpos
  // 		par[i+2] = zpos
  // 		par[i+3] = tau
  // 		par[i+4] = sigma.
  Int_t npeak = (Int_t)par[0];
  
  memset( spe, 0, sizeof( Float_t )*zdim*xdim );
  
  Int_t k = 1;
  for( Int_t i=0; i<npeak; i++ )
    {
      if( Integral != 0 ) Integral[i] = 0.;
      Float_t sigmaA2 = par[k+4]*par[k+4]*2.;
      Float_t T2 = par[k+3];   //PASCAL
      if(Electronics == 2) { T2 *= T2; T2 *= 2; } // OLA
      for( Int_t z=0; z<zdim; z++ ) {
	for( Int_t x=0; x<xdim; x++ ) {
	  Float_t z2 = (z-par[k+2])*(z-par[k+2])/sigmaA2;
	  Float_t x2 = 0.;
	  if(Electronics == 1) // PASCAL
	    x2 = (x-par[k+1]+T2)/T2;
	  else if(Electronics == 2) //OLA
	    x2 = (x-par[k+1])*(x-par[k+1])/T2;
	  else
	    cout << "Wrong Electronics" << endl;
	  // Float_t signal = (x2 > 0.) ? par[k] * x2*x2 * exp( -2*x2+2. - z2 ) : 0.0; // RCCR
	  Float_t signal = 0.;
          if(Electronics == 1)
	    signal = (x2 > 0.) ? par[k] * x2 * exp( -x2+1. - z2 ) : 0.0;
	  else if(Electronics == 2) //OLA
	    signal = par[k]  * exp( -x2 - z2 );
	  else
	    cout << "Wrong Electronics" << endl;
	  
	  spe[x*zdim+z] += signal;
	  if( Integral != 0 ) Integral[i] += signal; 
	}
      }	
      k += param_peak;
    }
  return;
}


/*
void AliITSClusterFinderSDD::PeakFunc( Int_t xdim, Int_t zdim, Float_t *par, Float_t *spe, Float_t *Integral=0 ) 
{
 Int_t param_peak = 5;
// par -> paramiters..
// par[0]  number of peaks.
// for each peak i=1, ..., par[0]
// 		par[i] = Ampl.
// 		par[i+1] = xpos
// 		par[i+2] = zpos
// 		par[i+3] = tau
// 		par[i+4] = sigma.
	Int_t npeak = (Int_t)par[0];

	memset( spe, 0, sizeof( Float_t )*zdim*xdim );
	
	Int_t k = 1;
	for( Int_t i=0; i<npeak; i++ )
	{
		if( Integral != 0 ) Integral[i] = 0.;
		Float_t sigmaA2 = par[k+4]*par[k+4]*2.;
		Float_t T2 = par[k+3]*par[k+3]*2.; 
		for( Int_t z=0; z<zdim; z++ )
		{
			for( Int_t x=0; x<xdim; x++ )
			{
				Float_t z2 = (z-par[k+2])*(z-par[k+2])/sigmaA2;
				Float_t x2 = (x-par[k+1])*(x-par[k+1])/T2;
				Float_t signal = par[k]  * exp( -x2 - z2 );
				spe[x*zdim+z] += signal;
				if( Integral != 0 ) Integral[i] += signal; 
			}
		}	
		k += param_peak;
	}
	return;
}
*/

Float_t AliITSClusterFinderSDD::chisq( Int_t xdim, Int_t zdim, Float_t *spe, Float_t *speFit )
{
	// EVALUATES UNNORMALIZED CHI-SQUARED
	
	Float_t chi2 = 0.;
	for( Int_t z=0; z<zdim; z++ )
	{
		for( Int_t x=1; x<xdim-1; x++ )
		{
			Int_t index = x*zdim+z;
			Float_t tmp = spe[index] - speFit[index];
			chi2 += tmp*tmp;
		}
	}	
	return( chi2 );
}


void AliITSClusterFinderSDD::minim( Int_t xdim, Int_t zdim, Float_t *param, Float_t *prm0, Float_t *steprm, Float_t *chisqr, 
		Float_t *spe, Float_t *speFit )
{
	Int_t   k, nnn, mmm, i;
	Float_t p1, delta, d1, chisq1, p2, chisq2, t, p3, chisq3, a, b, p0, chisqt;
	
	Int_t param_peak = 5;
	Int_t npeak = (Int_t)param[0];
	for( k=1; k<(npeak*param_peak+1); k++ ) prm0[k] = param[k];

	for( k=1; k<(npeak*param_peak+1); k++ ) 
	{
		p1 = param[k];
		delta = steprm[k];
		d1 = delta;

		// ENSURE THAT STEP SIZE IS SENSIBLY LARGER THAN MACHINE ROUND OFF
		if( fabs( p1 ) > 1.0E-6 ) 
   			if ( fabs( delta/p1 ) < 1.0E-4 ) delta = p1/1000;
		else  delta = (Float_t)1.0E-4;

		//  EVALUATE CHI-SQUARED AT FIRST TWO SEARCH POINTS
		PeakFunc( xdim, zdim, param, speFit );
		chisq1 = chisq( xdim, zdim, spe, speFit );

		p2 = p1+delta;
		param[k] = p2;

		PeakFunc( xdim, zdim, param, speFit );
		chisq2 = chisq( xdim, zdim, spe, speFit );

		if( chisq1 < chisq2 ) 
		{
    	// REVERSE DIRECTION OF SEARCH IF CHI-SQUARED IS INCREASING
    		delta = -delta;
    		t = p1;
    		p1 = p2;
    		p2 = t;
    		t = chisq1;
    		chisq1 = chisq2;
    		chisq2 = t;
    	}

		i = 1; nnn = 0;
		do {   // INCREMENT param(K) UNTIL CHI-SQUARED STARTS TO INCREASE
	    	nnn++;
   	 		p3 = p2 + delta;
    		mmm = nnn - (nnn/5)*5;  // multiplo de 5
    		if( mmm == 0 ) 
			{
    			d1 = delta;
    			// INCREASE STEP SIZE IF STEPPING TOWARDS MINIMUM IS TOO SLOW 
    			delta *= 5;
    		}
    		param[k] = p3;

			// Constrain paramiters
			Int_t kpos = (k-1) % param_peak;
			switch( kpos ) 
			{
				case 0 :
					if( param[k] <= 20 ) param[k] = fMinPeak;
				case 1 :
					if( fabs( param[k] - prm0[k] ) > 1.5 ) param[k] = prm0[k];
				case 2 :
					if( fabs( param[k] - prm0[k] ) > 1. ) param[k] = prm0[k];
				case 3 :
					if( param[k] < .5 ) param[k] = .5;	
				case 4 :
					if( param[k] < .288 ) param[k] = .288;	// 1/sqrt(12) = 0.288
			};
	
			PeakFunc( xdim, zdim, param, speFit );
			chisq3 = chisq( xdim, zdim, spe, speFit );
		
	    	if( chisq3 < chisq2 && nnn < 50 ) 
			{
   		 		p1 = p2;
    			p2 = p3;
    			chisq1 = chisq2;
    			chisq2 = chisq3;
    		}
    		else i=0;

    	} while( i );

		// FIND MINIMUM OF PARABOLA DEFINED BY LAST THREE POINTS
		a = chisq1*(p2-p3)+chisq2*(p3-p1)+chisq3*(p1-p2);
		b = chisq1*(p2*p2-p3*p3)+chisq2*(p3*p3-p1*p1)+chisq3*(p1*p1-p2*p2);
		if( a!=0 ) p0 = (Float_t)(0.5*b/a);
     		  else p0 = 10000;

		//---IN CASE OF NEARLY EQUAL CHI-SQUARED AND TOO SMALL STEP SIZE PREVENT
		//   ERRONEOUS EVALUATION OF PARABOLA MINIMUM
		//---NEXT TWO LINES CAN BE OMITTED FOR HIGHER PRECISION MACHINES

		//dp = (Float_t) max (fabs(p3-p2), fabs(p2-p1));
		//if( fabs( p2-p0 ) > dp ) p0 = p2;
		param[k] = p0;
		
		// Constrain paramiters
		Int_t kpos = (k-1) % param_peak;
		switch( kpos ) 
		{
			case 0 :
				if( param[k] <= 20 ) param[k] = fMinPeak;   
			case 1 :
				if( fabs( param[k] - prm0[k] ) > 1.5 ) param[k] = prm0[k];
			case 2 :
				if( fabs( param[k] - prm0[k] ) > 1. ) param[k] = prm0[k];
			case 3 :
				if( param[k] < .5 ) param[k] = .5;	
			case 4 :
				if( param[k] < .288 ) param[k] = .288;  // 1/sqrt(12) = 0.288	
		};
	
		PeakFunc( xdim, zdim, param, speFit );
		chisqt = chisq( xdim, zdim, spe, speFit );

		// DO NOT ALLOW ERRONEOUS INTERPOLATION
		if( chisqt <= *chisqr ) 
			*chisqr = chisqt;
   		else  
			param[k] = prm0[k];

		// OPTIMIZE SEARCH STEP FOR EVENTUAL NEXT CALL OF MINIM
		steprm[k] = (param[k]-prm0[k])/5;
		if( steprm[k] >= d1 ) steprm[k] = d1/5;
	}

	// EVALUATE FIT AND CHI-SQUARED FOR OPTIMIZED PARAMETERS
	PeakFunc( xdim, zdim, param, speFit );
	*chisqr = chisq( xdim, zdim, spe, speFit );
	return;
}

Int_t AliITSClusterFinderSDD::noLinearFit( Int_t xdim, Int_t zdim, Float_t *param, Float_t *spe, Int_t *niter, Float_t *chir )
{
	const Float_t chilmt = 0.01; 		   //	relative accuracy 	  
	const Int_t   nel = 3;				   //	for parabolic minimization  
	const Int_t   nstop = 50;			   //	Max. iteration number	  
	const Int_t   param_peak = 5;

	Int_t npeak = (Int_t)param[0];
	
	// RETURN IF NUMBER OF DEGREES OF FREEDOM IS NOT POSITIVE 
	if( (xdim*zdim - npeak*param_peak) <= 0 ) return( -1 );
	Float_t deg_free = (xdim*zdim - npeak*param_peak)-1;

	Int_t   n, k, iter_num = 0;
	Float_t *prm0 = new Float_t[npeak*param_peak+1];
	Float_t *step = new Float_t[npeak*param_peak+1];
	Float_t *schi = new Float_t[npeak*param_peak+1]; 
	Float_t *sprm[3];
	sprm[0] = new Float_t[npeak*param_peak+1];
	sprm[1] = new Float_t[npeak*param_peak+1];
	sprm[2] = new Float_t[npeak*param_peak+1];
	
	Float_t  chi0, chi1, reldif, a, b, prmin, dp;
	
	Float_t *speFit = new Float_t[ xdim*zdim ];
	PeakFunc( xdim, zdim, param, speFit );
	chi0 = chisq( xdim, zdim, spe, speFit );
	chi1 = chi0;


	for( k=1; k<(npeak*param_peak+1); k++) prm0[k] = param[k];

	for( k=1 ; k<(npeak*param_peak+1); k+=param_peak ) 
	{
    	step[k] = param[k] / 20.0 ; 		
    	step[k+1] = param[k+1] / 50.0;
		step[k+2] = param[k+2] / 50.0;		 
    	step[k+3] = param[k+3] / 20.0;		 
    	step[k+4] = param[k+4] / 20.0;		 
    }

	Int_t out = 0;
	do 
	{
	    iter_num++;
    	chi0 = chi1;
		
    	minim( xdim, zdim, param, prm0, step, &chi1, spe, speFit );
    	reldif = ( chi1 > 0 ) ? ((Float_t) fabs( chi1-chi0)/chi1 ) : 0;

	    // EXIT conditions
		if( reldif < (float) chilmt ) 	
		{
   		 	*chir  = (chi1>0) ? (float) TMath::Sqrt (chi1/deg_free) :0;
    		*niter = iter_num;
			out = 0;
			break;
    	}

    	if( (reldif < (float)(5*chilmt)) && (iter_num > nstop) ) 
		{
    		*chir = (chi1>0) ?(float) TMath::Sqrt (chi1/deg_free):0;
    		*niter = iter_num;
			out = 0;
			break;
    	}

    	if( iter_num > 5*nstop ) 
		{
    		*chir  = (chi1>0) ?(float) TMath::Sqrt (chi1/deg_free):0;
    		*niter = iter_num;
			out = 1;
			break;
    	}

    	if( iter_num <= nel ) continue;

    	n = iter_num - (iter_num/nel)*nel; // EXTRAPOLATION LIMIT COUNTER N
    	if( n > 3 || n == 0 ) continue;
    	schi[n-1] = chi1;
    	for( k=1; k<(npeak*param_peak+1); k++ ) sprm[n-1][k] = param[k];
    	if( n != 3 ) continue;

 		//   -EVALUATE EXTRAPOLATED VALUE OF EACH PARAMETER BY FINDING MINIMUM OF
 		//    PARABOLA DEFINED BY LAST THREE CALLS OF MINIM

    	for( k=1; k<(npeak*param_peak+1); k++ ) 
		{
    		Float_t tmp0 = sprm[0][k];
    		Float_t tmp1 = sprm[1][k];
    		Float_t tmp2 = sprm[2][k];
    		a  = schi[0]*(tmp1-tmp2) + schi[1]*(tmp2-tmp0);
    		a += (schi[2]*(tmp0-tmp1));
    		b  = schi[0]*(tmp1*tmp1-tmp2*tmp2);
    		b += (schi[1]*(tmp2*tmp2-tmp0*tmp0)+(schi[2]*(tmp0*tmp0-tmp1*tmp1)));
    		if ((double)a < 1.0E-6) prmin = 0;
    			   else prmin = (float) (0.5*b/a);
    		dp = 5*(tmp2-tmp0);

    		if (fabs(prmin-tmp2) > fabs(dp)) prmin = tmp2+dp;
    		param[k] = prmin;
    		step[k]  = dp/10; // OPTIMIZE SEARCH STEP
    	}

	} while( kTRUE );

	delete [] prm0;
	delete [] step;
	delete [] schi; 
	delete [] sprm[0];
	delete [] sprm[1];
	delete [] sprm[2];
	delete [] speFit;

	return( out );
}
//_____________________________________________________________________________
void AliITSClusterFinderSDD::ResolveClustersE()
{
	// The function to resolve clusters if the clusters overlapping exists

  Int_t i;

	AliITS *iTS = (AliITS*)gAlice->GetModule( "ITS" );
	// get number of clusters for this module
	Int_t nofClusters = fClusters->GetEntriesFast();
	nofClusters -= fNclusters;

	Int_t fNofMaps = fSegmentation->Npz();
	Int_t fNofAnodes = fNofMaps/2;
	Int_t fMaxNofSamples = fSegmentation->Npx();
	Int_t dummy=0;
	Double_t fTimeStep = fSegmentation->Dpx( dummy );
	Double_t fSddLength = fSegmentation->Dx();
	Double_t fDriftSpeed = fResponse->DriftSpeed();
	Double_t anodePitch = fSegmentation->Dpz( dummy );
	Float_t n, baseline;
	fResponse->GetNoiseParam( n, baseline );
	Int_t Electronics = fResponse->Electronics(); // 1 = PASCAL, 2 = OLA
	
	// fill Map of signals
	fMap->FillMap(); 

	for( Int_t j=0; j<nofClusters; j++ ) 
	{ 
		// get cluster information
		AliITSRawClusterSDD *clusterJ = (AliITSRawClusterSDD*) fClusters->At( j );
		Int_t astart = clusterJ->Astart();
		Int_t astop = clusterJ->Astop();
		Int_t tstart = clusterJ->Tstartf();
		Int_t tstop = clusterJ->Tstopf();
		Int_t wing = (Int_t)clusterJ->W();
		if( wing == 2 ) 
		{
			astart += fNofAnodes; 
			astop  += fNofAnodes;
		} 
		Int_t xdim = tstop-tstart+3;
		Int_t zdim = astop-astart+3;
		Float_t *sp = new Float_t[ xdim*zdim+1 ];
		memset( sp, 0, sizeof(Float_t)*(xdim*zdim+1) );

		// make a local map from cluster region
		for( Int_t ianode=astart; ianode<=astop; ianode++ )
		{
			for( Int_t itime=tstart; itime<=tstop; itime++ )
			{
				Float_t fadc = fMap->GetSignal( ianode, itime );
				if( fadc > baseline ) fadc -= (Double_t)baseline;
				else fadc = 0.;
				Int_t index = (itime-tstart+1)*zdim+(ianode-astart+1);
				sp[index] = fadc;
			} // time loop
		} // anode loop
		
		// search peaks on cluster
		const Int_t np = 150;
		Int_t peakX1[np];
		Int_t peakZ1[np];
		Float_t peakAmp1[np];
		Int_t npeak = SearchPeak( sp, xdim, zdim, peakX1, peakZ1, peakAmp1, fMinPeak );

		// if multiple peaks, split cluster
		if( npeak >= 1 )
		{
		//	cout << "npeak " << npeak << endl;
		//	clusterJ->PrintInfo();
			
			Float_t *par = new Float_t[npeak*5+1];
			par[0] = (Float_t)npeak;		
			
			// Initial paramiters in cell dimentions
			Int_t k1 = 1;
			for( i=0; i<npeak; i++ ) {
			  par[k1] = peakAmp1[i];
			  par[k1+1] = peakX1[i]; // local time pos. [timebin]
			  par[k1+2] = peakZ1[i]; // local anode pos. [anodepitch]
			  if(Electronics == 1) 
			    par[k1+3] = 2.; // PASCAL
			  else if(Electronics == 2) 
			    par[k1+3] = 0.7; // tau [timebin]  OLA 
			  par[k1+4] = .4;    // sigma	[anodepich]
			  k1+=5;
			}			
			Int_t niter;
			Float_t chir;			
			noLinearFit( xdim, zdim, par, sp, &niter, &chir );

			Float_t peakX[np];
			Float_t peakZ[np];
			Float_t sigma[np];
			Float_t tau[np];
			Float_t peakAmp[np];
			Float_t Integral[np];
			
			//get integrals => charge for each peak
			PeakFunc( xdim, zdim, par, sp, Integral );
			
			k1 = 1;
			for( i=0; i<npeak; i++ ) 
			{
				peakAmp[i] = par[k1];
				peakX[i] = par[k1+1];
				peakZ[i] = par[k1+2];
				tau[i] = par[k1+3];
				sigma[i] = par[k1+4];
				k1+=5;
			}
			
			// calculate paramiter for new clusters
			for( i=0; i<npeak; i++ )
			{
				AliITSRawClusterSDD clusterI( *clusterJ );
				Int_t newAnode = peakZ1[i]-1 + astart;
				Int_t newiTime = peakX1[i]-1 + tstart;

				Int_t shift = (Int_t)(fTimeCorr/fTimeStep + 0.5);
				if( newiTime > shift && newiTime < (fMaxNofSamples-shift) )  shift = 0;
				Int_t peakpos = fMap->GetHitIndex( newAnode, newiTime+shift );
				clusterI.SetPeakPos( peakpos );
				clusterI.SetPeakAmpl( peakAmp1[i] );

				Float_t newAnodef = peakZ[i] - 0.5 + astart;
				Float_t newiTimef = peakX[i] - 1 + tstart;				
				if( wing == 2 ) newAnodef -= fNofAnodes; 
				Float_t AnodePath = (newAnodef - fNofAnodes/2)*anodePitch;
				newiTimef *= fTimeStep;
				if( newiTimef > fTimeCorr ) newiTimef -= fTimeCorr;
				if(Electronics == 1) {
				  newiTimef *= 0.999438;    // PASCAL
				  newiTimef += (6./fDriftSpeed - newiTimef/3000.); 
				}
				else if(Electronics == 2)
				  newiTimef *= 0.99714;    // OLA

				Float_t DriftPath = fSddLength - newiTimef * fDriftSpeed;
				Float_t sign = ( wing == 1 ) ? -1. : 1.;
				clusterI.SetX( DriftPath*sign * 0.0001 );	
				clusterI.SetZ( AnodePath * 0.0001 );
				clusterI.SetAnode( newAnodef );
				clusterI.SetTime( newiTimef );
				clusterI.SetAsigma( sigma[i]*anodePitch );
				clusterI.SetTsigma( tau[i]*fTimeStep );
				clusterI.SetQ( Integral[i] );
				
			//	clusterI.PrintInfo();
				iTS->AddCluster( 1, &clusterI );
			}
			fClusters->RemoveAt( j );
   			delete [] par;
		}
		else cout <<" --- Peak not found!!!!  minpeak=" << fMinPeak<< 
		            " cluster peak=" << clusterJ->PeakAmpl() << endl << endl;
		
		delete [] sp;
	} // cluster loop

	fClusters->Compress();
	fMap->ClearMap(); 
}


//_____________________________________________________________________________
void  AliITSClusterFinderSDD::GroupClusters()
{
  // group clusters
  Int_t dummy=0;
  Float_t fTimeStep = fSegmentation->Dpx(dummy);


  // get number of clusters for this module
  Int_t nofClusters = fClusters->GetEntriesFast();
  nofClusters -= fNclusters;

  AliITSRawClusterSDD *clusterI;
  AliITSRawClusterSDD *clusterJ;

  Int_t *label = new Int_t [nofClusters];
  Int_t i,j;
  for(i=0; i<nofClusters; i++) label[i] = 0;
  for(i=0; i<nofClusters; i++) { 
    if(label[i] != 0) continue;
    for(j=i+1; j<nofClusters; j++) { 
      if(label[j] != 0) continue;
      clusterI = (AliITSRawClusterSDD*) fClusters->At(i);
      clusterJ = (AliITSRawClusterSDD*) fClusters->At(j);
      // 1.3 good
      if(clusterI->T() < fTimeStep*60) fDAnode = 4.2;  // TB 3.2  
      if(clusterI->T() < fTimeStep*10) fDAnode = 1.5;  // TB 1.
      Bool_t pair = clusterI->Brother(clusterJ,fDAnode,fDTime);
      if(!pair) continue;
      //      clusterI->PrintInfo();
      //      clusterJ->PrintInfo();
      clusterI->Add(clusterJ);
      label[j] = 1;
      fClusters->RemoveAt(j);
      j=i; // <- Ernesto
    } // J clusters  
    label[i] = 1;
  } // I clusters
  fClusters->Compress();
  
  delete [] label;
  return;

}

//_____________________________________________________________________________

void AliITSClusterFinderSDD::SelectClusters()
{
  // get number of clusters for this module
  Int_t nofClusters = fClusters->GetEntriesFast();
  nofClusters -= fNclusters;

  Int_t i;
  for(i=0; i<nofClusters; i++) { 
    AliITSRawClusterSDD *clusterI = (AliITSRawClusterSDD*) fClusters->At(i);
    Int_t rmflg = 0;
    Float_t wy = 0.;
    if(clusterI->Anodes() != 0.) {
      wy = ((Float_t) clusterI->Samples())/clusterI->Anodes();
    }
    Int_t amp = (Int_t) clusterI->PeakAmpl();
    Int_t cha = (Int_t) clusterI->Q();
    if(amp < fMinPeak) rmflg = 1;  
    if(cha < fMinCharge) rmflg = 1;
    if(wy < fMinNCells) rmflg = 1;
    //if(wy > fMaxNCells) rmflg = 1;
    if(rmflg) fClusters->RemoveAt(i);
  } // I clusters
  fClusters->Compress();
  return;

}

//_____________________________________________________________________________

void AliITSClusterFinderSDD::ResolveClusters()
{

// The function to resolve clusters if the clusters overlapping exists

   AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");

  // get number of clusters for this module
   Int_t nofClusters = fClusters->GetEntriesFast();
   nofClusters -= fNclusters;
   //   cout<<"Resolve Cl: nofClusters, fNclusters ="<<nofClusters<<","<<fNclusters<<endl;

   Int_t fNofMaps = fSegmentation->Npz();
   Int_t fNofAnodes = fNofMaps/2;
   Int_t dummy=0;
   Double_t fTimeStep = fSegmentation->Dpx(dummy);
   Double_t fSddLength = fSegmentation->Dx();
   Double_t fDriftSpeed = fResponse->DriftSpeed();
   Double_t anodePitch = fSegmentation->Dpz(dummy);
   Float_t n, baseline;
   fResponse->GetNoiseParam(n,baseline);
   Float_t dzz_1A = anodePitch * anodePitch / 12;

  // fill Map of signals
    fMap->FillMap(); 

  Int_t j,i,ii,ianode,anode,itime;
  Int_t wing,astart,astop,tstart,tstop,nanode;
  Double_t fadc,ClusterTime;
  Double_t q[400],x[400],z[400]; // digit charges and coordinates

  for(j=0; j<nofClusters; j++) { 

    AliITSRawClusterSDD *clusterJ = (AliITSRawClusterSDD*) fClusters->At(j);

    Int_t ndigits = 0;
    astart=clusterJ->Astart();
    astop=clusterJ->Astop();
    tstart=clusterJ->Tstartf();
    tstop=clusterJ->Tstopf();
    nanode=clusterJ->Anodes();  // <- Ernesto
    wing=(Int_t)clusterJ->W();
    if(wing == 2) {
	astart += fNofAnodes; 
	astop  += fNofAnodes;
    } 

//	    cout<<"astart,astop,tstart,tstop ="<<astart<<","<<astop<<","<<tstart<<","<<tstop<<endl;

	    // clear the digit arrays
   for(ii=0; ii<400; ii++) { 
       q[ii] = 0.; 
       x[ii] = 0.;
       z[ii] = 0.;
   }
    
   for(ianode=astart; ianode<=astop; ianode++) { 
    for(itime=tstart; itime<=tstop; itime++) { 
	  fadc=fMap->GetSignal(ianode,itime);
	  if(fadc>baseline) {
	    fadc-=(Double_t)baseline;
	      q[ndigits] = fadc*(fTimeStep/160);  // KeV
              anode = ianode;
              if(wing == 2) anode -= fNofAnodes;
	      z[ndigits] = (anode + 0.5 - fNofAnodes/2)*anodePitch;
	      ClusterTime = itime*fTimeStep;
	      if(ClusterTime > fTimeCorr) ClusterTime -= fTimeCorr;   // ns
	      x[ndigits] = fSddLength - ClusterTime*fDriftSpeed;
	      if(wing == 1) x[ndigits] *= (-1);
//              cout<<"ianode,itime,fadc ="<<ianode<<","<<itime<<","<<fadc<<endl;
//	      cout<<"wing,anode,ndigits,charge ="<<wing<<","<<anode<<","<<ndigits<<","<<q[ndigits]<<endl;
	      ndigits++;
	      continue;
	  }
	      fadc=0;
	      //	      cout<<"fadc=0, ndigits ="<<ndigits<<endl;
    } // time loop
   } // anode loop
   //     cout<<"for new cluster ndigits ="<<ndigits<<endl;


   // Fit cluster to resolve for two separate ones --------------------

   Double_t qq=0., xm=0., zm=0., xx=0., zz=0., xz=0.;
   Double_t dxx=0., dzz=0., dxz=0.;
   Double_t scl = 0., tmp, tga, elps = -1.;
   Double_t xfit[2], zfit[2], qfit[2];
   Double_t pitchz = anodePitch*1.e-4;             // cm
   Double_t pitchx = fTimeStep*fDriftSpeed*1.e-4;  // cm
   Double_t sigma2;
   Int_t nfhits;
   Int_t nbins = ndigits;
   Int_t separate = 0;

   // now, all lengths are in microns

   for (ii=0; ii<nbins; ii++) {
       qq += q[ii];
       xm += x[ii]*q[ii];
       zm += z[ii]*q[ii];
       xx += x[ii]*x[ii]*q[ii];
       zz += z[ii]*z[ii]*q[ii];
       xz += x[ii]*z[ii]*q[ii];
   }

   xm /= qq;
   zm /= qq;
   xx /= qq;
   zz /= qq;
   xz /= qq;

   dxx = xx - xm*xm;
   dzz = zz - zm*zm;
   dxz = xz - xm*zm;
   
   // shrink the cluster in the time direction proportionaly to the 
   // dxx/dzz, which lineary depends from the drift path
 
	// new  Ernesto........	 
   if( nanode == 1 )
   {
 	   dzz = dzz_1A; // for one anode cluster dzz = anode**2/12
 	   scl = TMath::Sqrt( 7.2/(-0.57*xm*1.e-3+71.8) );
   }
   if( nanode == 2 )
   {
 	   scl = TMath::Sqrt( (-0.18*xm*1.e-3+21.3)/(-0.57*xm*1.e-3+71.8) );
   }
   
   if( nanode == 3 )	   
   {
 	   scl = TMath::Sqrt( (-0.5*xm*1.e-3+34.5)/(-0.57*xm*1.e-3+71.8) );
   }

   if( nanode > 3 )	   
   {
 	   scl = TMath::Sqrt( (1.3*xm*1.e-3+49.)/(-0.57*xm*1.e-3+71.8) );
   }

   //   cout<<"1 microns: zm,dzz,xm,dxx,dxz,qq ="<<zm<<","<<dzz<<","<<xm<<","<<dxx<<","<<dxz<<","<<qq<<endl;

 //  old Boris.........
 //  tmp=29730. - 585.*fabs(xm/1000.); 
 //  scl=TMath::Sqrt(tmp/130000.);
   
   xm *= scl;
   xx *= scl*scl;
   xz *= scl;

   
   dxx = xx - xm*xm;
//   dzz = zz - zm*zm;
   dxz = xz - xm*zm;

   //   cout<<"microns: zm,dzz,xm,dxx,xz,dxz,qq ="<<zm<<","<<dzz<<","<<xm<<","<<dxx<<","<<xz<<","<<dxz<<","<<qq<<endl;

//   if(dzz < 7200.) dzz = 7200.; // for one anode cluster dzz = anode**2/12
  
   if (dxx < 0.) dxx=0.;

   // the data if no cluster overlapping (the coordunates are in cm) 
   nfhits = 1;
   xfit[0] = xm*1.e-4;
   zfit[0] = zm*1.e-4;
   qfit[0] = qq;

//   if(nbins < 7) cout<<"**** nbins ="<<nbins<<endl;
  
   if (nbins >= 7) {
      if (dxz==0.) tga=0.;
      else {
         tmp=0.5*(dzz-dxx)/dxz;
         tga = (dxz<0.) ? tmp-TMath::Sqrt(tmp*tmp+1) : tmp+TMath::Sqrt(tmp*tmp+1);
      }
      elps=(tga*tga*dxx-2*tga*dxz+dzz)/(dxx+2*tga*dxz+tga*tga*dzz);

      // change from microns to cm
      xm *= 1.e-4; 
      zm *= 1.e-4; 
      zz *= 1.e-8;
      xx *= 1.e-8;
      xz *= 1.e-8;
      dxz *= 1.e-8;
      dxx *= 1.e-8;
      dzz *= 1.e-8;

   //   cout<<"cm: zm,dzz,xm,dxx,xz,dxz,qq ="<<zm<<","<<dzz<<","<<xm<<","<<dxx<<","<<xz<<","<<dxz<<","<<qq<<endl;

     for (i=0; i<nbins; i++) {     
       x[i] = x[i] *= scl;
       x[i] = x[i] *= 1.e-4;
       z[i] = z[i] *= 1.e-4;
     }

     //     cout<<"!!! elps ="<<elps<<endl;

     if (elps < 0.3) { // try to separate hits 
         separate = 1;
         tmp=atan(tga);
         Double_t cosa=cos(tmp),sina=sin(tmp);
         Double_t a1=0., x1=0., xxx=0.;
         for (i=0; i<nbins; i++) {
             tmp=x[i]*cosa + z[i]*sina;
             if (q[i] > a1) {
                a1=q[i];
                x1=tmp;
             }
             xxx += tmp*tmp*tmp*q[i];
         }
         xxx /= qq;
         Double_t z12=-sina*xm + cosa*zm;
         sigma2=(sina*sina*xx-2*cosa*sina*xz+cosa*cosa*zz) - z12*z12;
         xm=cosa*xm + sina*zm;
         xx=cosa*cosa*xx + 2*cosa*sina*xz + sina*sina*zz;
         Double_t x2=(xx - xm*x1 - sigma2)/(xm - x1);
         Double_t r=a1*2*TMath::ACos(-1.)*sigma2/(qq*pitchx*pitchz);
         for (i=0; i<33; i++) { // solve a system of equations
             Double_t x1_old=x1, x2_old=x2, r_old=r;
             Double_t c11=x1-x2;
             Double_t c12=r;
             Double_t c13=1-r;
             Double_t c21=x1*x1 - x2*x2;
             Double_t c22=2*r*x1;
             Double_t c23=2*(1-r)*x2;
             Double_t c31=3*sigma2*(x1-x2) + x1*x1*x1 - x2*x2*x2;
             Double_t c32=3*r*(sigma2 + x1*x1);
             Double_t c33=3*(1-r)*(sigma2 + x2*x2);
             Double_t f1=-(r*x1 + (1-r)*x2 - xm);
             Double_t f2=-(r*(sigma2 + x1*x1) + (1-r)*(sigma2 + x2*x2) - xx);
             Double_t f3=-(r*x1*(3*sigma2+x1*x1) + (1-r)*x2*(3*sigma2+x2*x2)-xxx);
             Double_t d=c11*c22*c33 + c21*c32*c13 + c12*c23*c31 - c31*c22*c13 - c21*c12*c33 - c32*c23*c11;
             if (d==0.) {
                cout<<"*********** d=0 ***********\n";
                break;
             }
             Double_t dr=f1*c22*c33 + f2*c32*c13 + c12*c23*f3 -
                       f3*c22*c13 - f2*c12*c33 - c32*c23*f1;
             Double_t d1=c11*f2*c33 + c21*f3*c13 + f1*c23*c31 -
                       c31*f2*c13 - c21*f1*c33 - f3*c23*c11;
             Double_t d2=c11*c22*f3 + c21*c32*f1 + c12*f2*c31 -
                       c31*c22*f1 - c21*c12*f3 - c32*f2*c11;
             r  += dr/d;
             x1 += d1/d;
             x2 += d2/d;

             if (fabs(x1-x1_old) > 0.0001) continue;
             if (fabs(x2-x2_old) > 0.0001) continue;
             if (fabs(r-r_old)/5 > 0.001) continue;

             a1=r*qq*pitchx*pitchz/(2*TMath::ACos(-1.)*sigma2);
             Double_t a2=a1*(1-r)/r;
             qfit[0]=a1; xfit[0]=x1*cosa - z12*sina; zfit[0]=x1*sina + z12*cosa;
             qfit[1]=a2; xfit[1]=x2*cosa - z12*sina; zfit[1]=x2*sina + z12*cosa;
             nfhits=2;
             break; // Ok !
         }
         if (i==33) cerr<<"No more iterations ! "<<endl;
    } // end of attempt to separate overlapped clusters
   } // end of nbins cut 

   if(elps < 0.) cout<<" elps=-1 ="<<elps<<endl;
   if(elps >0. && elps< 0.3 && nfhits == 1) cout<<" small elps, nfh=1 ="<<elps<<","<<nfhits<<endl;
   if(nfhits == 2) cout<<" nfhits=2 ="<<nfhits<<endl;

   for (i=0; i<nfhits; i++) {
       xfit[i] *= (1.e+4/scl);
       if(wing == 1) xfit[i] *= (-1);
       zfit[i] *= 1.e+4;
       //       cout<<" ---------  i,xfiti,zfiti,qfiti ="<<i<<","<<xfit[i]<<","<<zfit[i]<<","<<qfit[i]<<endl;
   }

    Int_t ncl = nfhits;
    if(nfhits == 1 && separate == 1) {
      cout<<"!!!!! no separate"<<endl;
      ncl = -2;
    } 

   if(nfhits == 2) {
     cout << "Split cluster: " << endl;
     clusterJ->PrintInfo();
     cout << " in: " << endl;
     for (i=0; i<nfhits; i++) {

	// AliITSRawClusterSDD *clust = new AliITSRawClusterSDD(wing,-1,-1,(Float_t)qfit[i],ncl,0,0,(Float_t)xfit[i],(Float_t)zfit[i],0,0,0,0,tstart,tstop,astart,astop);
	//	AliITSRawClusterSDD *clust = new AliITSRawClusterSDD(wing,-1,-1,(Float_t)qfit[i],0,0,0,(Float_t)xfit[i],(Float_t)zfit[i],0,0,0,0,tstart,tstop,astart,astop,ncl);

	// ???????????
       // if(wing == 1) xfit[i] *= (-1);
	      Float_t Anode = (zfit[i]/anodePitch+fNofAnodes/2-0.5);
	      Float_t Time = (fSddLength - xfit[i])/fDriftSpeed;
	Float_t clusterPeakAmplitude = clusterJ->PeakAmpl();
	Float_t peakpos = clusterJ->PeakPos();

	Float_t clusteranodePath = (Anode - fNofAnodes/2)*anodePitch;
	Float_t clusterDriftPath = Time*fDriftSpeed;
	clusterDriftPath = fSddLength-clusterDriftPath;

	AliITSRawClusterSDD *clust = new AliITSRawClusterSDD(wing,Anode,Time,qfit[i],
	                            clusterPeakAmplitude,peakpos,0.,0.,clusterDriftPath,clusteranodePath,clusterJ->Samples()/2
				    ,tstart,tstop,0,0,0,astart,astop);
	clust->PrintInfo();
	iTS->AddCluster(1,clust);
	//	cout<<"new cluster added: tstart,tstop,astart,astop,x,ncl ="<<tstart<<","<<tstop<<","<<astart<<","<<astop<<","<<xfit[i]<<","<<ncl<<endl;
	delete clust;
     }// nfhits loop
     fClusters->RemoveAt(j);

   } // if nfhits = 2
  } // cluster loop

  fClusters->Compress();
  fMap->ClearMap(); 
	    
  return;
}


//_____________________________________________________________________________

void AliITSClusterFinderSDD::GetRecPoints()
{
  // get rec points

  AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");

  // get number of clusters for this module
  Int_t nofClusters = fClusters->GetEntriesFast();
  nofClusters -= fNclusters;

  const Float_t kconvGeV = 1.e-6; // GeV -> KeV
  const Float_t kconv = 1.0e-4; 
  const Float_t kRMSx = 38.0*kconv; // microns->cm ITS TDR Table 1.3
  const Float_t kRMSz = 28.0*kconv; // microns->cm ITS TDR Table 1.3


  Int_t i;
  Int_t ix, iz, idx=-1;
  AliITSdigitSDD *dig=0;
  Int_t ndigits=fDigits->GetEntriesFast();
  for(i=0; i<nofClusters; i++) { 
    AliITSRawClusterSDD *clusterI = (AliITSRawClusterSDD*)fClusters->At(i);
    if(!clusterI) Error("SDD: GetRecPoints","i clusterI ",i,clusterI);
    if(clusterI) idx=clusterI->PeakPos();
    if(idx>ndigits) Error("SDD: GetRecPoints","idx ndigits",idx,ndigits);
    // try peak neighbours - to be done 
    if(idx && idx <= ndigits) dig = (AliITSdigitSDD*)fDigits->UncheckedAt(idx);
    if(!dig) {
	// try cog
	fSegmentation->GetPadIxz(clusterI->X(),clusterI->Z(),ix,iz);
	dig = (AliITSdigitSDD*)fMap->GetHit(iz-1,ix-1);
	// if null try neighbours
	if (!dig) dig = (AliITSdigitSDD*)fMap->GetHit(iz-1,ix); 
	if (!dig) dig = (AliITSdigitSDD*)fMap->GetHit(iz-1,ix+1); 
        if (!dig) printf("SDD: cannot assign the track number!\n");
    }

    AliITSRecPoint rnew;
    rnew.SetX(clusterI->X());
    rnew.SetZ(clusterI->Z());
    rnew.SetQ(clusterI->Q());   // in KeV - should be ADC
    rnew.SetdEdX(kconvGeV*clusterI->Q());
    rnew.SetSigmaX2(kRMSx*kRMSx);
    rnew.SetSigmaZ2(kRMSz*kRMSz);
    if(dig) rnew.fTracks[0]=dig->fTracks[0];
    if(dig) rnew.fTracks[1]=dig->fTracks[1];
    if(dig) rnew.fTracks[2]=dig->fTracks[2];
    //printf("SDD: i %d track1 track2 track3 %d %d %d x y %f %f\n",i,rnew.fTracks[0],rnew.fTracks[1],rnew.fTracks[2],clusterI->X(),clusterI->Z());
    iTS->AddRecPoint(rnew);
  } // I clusters

  fMap->ClearMap();
}

//_____________________________________________________________________________

void AliITSClusterFinderSDD::FindRawClusters()
{
  // find raw clusters
    Find1DClustersE();
    GroupClusters();
    SelectClusters();
    ResolveClustersE();
    GetRecPoints();
}
//_____________________________________________________________________________

void AliITSClusterFinderSDD::Print()
{
  // Print SDD cluster finder Parameters

   cout << "**************************************************" << endl;
   cout << " Silicon Drift Detector Cluster Finder Parameters " << endl;
   cout << "**************************************************" << endl;
   cout << "Number of Clusters: " << fNclusters << endl;
   cout << "Anode Tolerance: " << fDAnode << endl;
   cout << "Time  Tolerance: " << fDTime << endl;
   cout << "Time  correction (electronics): " << fTimeCorr << endl;
   cout << "Cut Amplitude (threshold): " << fCutAmplitude << endl;
   cout << "Minimum Amplitude: " << fMinPeak << endl;
   cout << "Minimum Charge: " << fMinCharge << endl;
   cout << "Minimum number of cells/clusters: " << fMinNCells << endl;
   cout << "Maximum number of cells/clusters: " << fMaxNCells << endl;
   cout << "**************************************************" << endl;
}
