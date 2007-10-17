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
  $Id$
 */
/////////////////////////////////////////////////////////////////////////// 
//  Cluster finder                                                       //
//  for Silicon                                                          //
//  Drift Detector                                                       //
////////////////////////////////////////////////////////////////////////// 


#include "AliITSClusterFinderSDD.h"
#include "AliITSMapA1.h"
#include "AliITSRawClusterSDD.h"
#include "AliITSRecPoint.h"
#include "AliITSdigitSDD.h"
#include "AliITSDetTypeRec.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSgeom.h"
#include "AliLog.h"

ClassImp(AliITSClusterFinderSDD)

//______________________________________________________________________
AliITSClusterFinderSDD::AliITSClusterFinderSDD():
AliITSClusterFinder(),
fNclusters(0),
fDAnode(0.0),
fDTime(0.0),
fTimeCorr(0.0),
fCutAmplitude(0),
fMinPeak(0),
fMinCharge(0),
fMinNCells(0),
fMaxNCells(0){
    // default constructor
}
//______________________________________________________________________
AliITSClusterFinderSDD::AliITSClusterFinderSDD(AliITSDetTypeRec* dettyp,
                                               TClonesArray *digits,
                                               TClonesArray *recp):
AliITSClusterFinder(dettyp),
fNclusters(0),
fDAnode(0.0),
fDTime(0.0),
fTimeCorr(0.0),
fCutAmplitude(0),
fMinPeak(0),
fMinCharge(0),
fMinNCells(0),
fMaxNCells(0){
    // standard constructor

    SetDigits(digits);
    SetClusters(recp);
    SetCutAmplitude(fDetTypeRec->GetITSgeom()->GetStartSDD());
    SetDAnode();
    SetDTime();
    SetMinPeak((Int_t)((AliITSCalibrationSDD*)GetResp(fDetTypeRec->GetITSgeom()->GetStartSDD()))->GetNoiseAfterElectronics(0)*5);
    SetMinNCells();
    SetMaxNCells();
    SetTimeCorr();
    SetMinCharge();
    SetMap(new AliITSMapA1(GetSeg(),Digits(),fCutAmplitude));
}
//______________________________________________________________________
void AliITSClusterFinderSDD::SetCutAmplitude(Int_t mod,Double_t nsigma){
    // set the signal threshold for cluster finder
    Double_t baseline,noiseAfterEl;

    AliITSresponseSDD* res  = (AliITSresponseSDD*)((AliITSCalibrationSDD*)GetResp(mod))->GetResponse();
    const char *option=res->ZeroSuppOption();
    Int_t nanodes = GetResp(mod)->Wings()*GetResp(mod)->Channels()*GetResp(mod)->Chips();
    fCutAmplitude.Set(nanodes);
    for(Int_t ian=0;ian<nanodes;ian++){
      noiseAfterEl = ((AliITSCalibrationSDD*)GetResp(mod))->GetNoiseAfterElectronics(ian);
      if((strstr(option,"1D")) || (strstr(option,"2D"))){ 
	fCutAmplitude[ian] = (Int_t)(nsigma*noiseAfterEl);
      }
      else{
	baseline=GetResp(mod)->GetBaseline(ian);
	fCutAmplitude[ian] = (Int_t)((baseline + nsigma*noiseAfterEl));
      } 
    }
}
//______________________________________________________________________
void AliITSClusterFinderSDD::Find1DClusters(){
    // find 1D clusters
  
    // retrieve the parameters 
    Int_t fNofMaps       = GetSeg()->Npz();
    Int_t fMaxNofSamples = GetSeg()->Npx();
    Int_t fNofAnodes     = fNofMaps/2;
    Int_t dummy          = 0;
    Double_t fTimeStep    = GetSeg()->Dpx(dummy);
    Double_t fSddLength   = GetSeg()->Dx();
    Double_t anodePitch   = GetSeg()->Dpz(dummy);
    AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
    AliITSresponseSDD* res  = (AliITSresponseSDD*)((AliITSCalibrationSDD*)GetResp(fModule))->GetResponse();
    const char *option=res->ZeroSuppOption();

    // map the signal
    Map()->ClearMap();
    Map()->SetThresholdArr(fCutAmplitude);
    Map()->FillMap2();
  
    Int_t nofFoundClusters = 0;
    Int_t i;
    Double_t **dfadc = new Double_t*[fNofAnodes];
    for(i=0;i<fNofAnodes;i++) dfadc[i] = new Double_t[fMaxNofSamples];
    Double_t fadc  = 0.;
    Double_t fadc1 = 0.;
    Double_t fadc2 = 0.;
    Int_t j,k,idx,l,m;
    for(j=0;j<2;j++) {
        for(k=0;k<fNofAnodes;k++) {
            idx = j*fNofAnodes+k;
            // signal (fadc) & derivative (dfadc)
            dfadc[k][255]=0.;
            for(l=0; l<fMaxNofSamples; l++) {
                fadc2=(Double_t)Map()->GetSignal(idx,l);
                if(l>0) fadc1=(Double_t)Map()->GetSignal(idx,l-1);
                if(l>0) dfadc[k][l-1] = fadc2-fadc1;
            } // samples
        } // anodes

        for(k=0;k<fNofAnodes;k++) {
            AliDebug(5,Form("Anode: %d, Wing: %d",k+1,j+1));
            idx = j*fNofAnodes+k;
            Int_t imax  = 0;
            Int_t imaxd = 0;
            Int_t it    = 0;
            while(it <= fMaxNofSamples-3) {
                imax  = it;
                imaxd = it;
                // maximum of signal          
                Double_t fadcmax  = 0.;
                Double_t dfadcmax = 0.;
                Int_t lthrmina   = 1;
                Int_t lthrmint   = 3;
                Int_t lthra      = 1;
                Int_t lthrt      = 0;
                for(m=0;m<20;m++) {
                    Int_t id = it+m;
                    if(id>=fMaxNofSamples) break;
                    fadc=(float)Map()->GetSignal(idx,id);
                    if(fadc > fadcmax) { fadcmax = fadc; imax = id;}
                    if(fadc > (float)fCutAmplitude[idx])lthrt++; 
                    if(dfadc[k][id] > dfadcmax) {
                        dfadcmax = dfadc[k][id];
                        imaxd    = id;
                    } // end if
                } // end for m
                it = imaxd;
                if(Map()->TestHit(idx,imax) == kEmpty) {it++; continue;}
                // cluster charge
                Int_t tstart = it-2;
                if(tstart < 0) tstart = 0;
                Bool_t ilcl = 0;
                if(lthrt >= lthrmint && lthra >= lthrmina) ilcl = 1;
                if(ilcl) {
                    nofFoundClusters++;
                    Int_t tstop      = tstart;
                    Double_t dfadcmin = 10000.;
                    Int_t ij;
                    for(ij=0; ij<20; ij++) {
                        if(tstart+ij > 255) { tstop = 255; break; }
                        fadc=(float)Map()->GetSignal(idx,tstart+ij);
                        if((dfadc[k][tstart+ij] < dfadcmin) && 
                           (fadc > fCutAmplitude[idx])) {
                            tstop = tstart+ij+5;
                            if(tstop > 255) tstop = 255;
                            dfadcmin = dfadc[k][it+ij];
                        } // end if
                    } // end for ij

                    Double_t clusterCharge = 0.;
                    Double_t clusterAnode  = k+0.5;
                    Double_t clusterTime   = 0.;
                    Int_t   clusterMult   = 0;
                    Double_t clusterPeakAmplitude = 0.;
                    Int_t its,peakpos     = -1;
		    		   
		    for(its=tstart; its<=tstop; its++) {
                        fadc=(float)Map()->GetSignal(idx,its);
			if(!((strstr(option,"1D")) || (strstr(option,"2D")))){
			  Double_t baseline=GetResp(fModule)->GetBaseline(idx);
			  if(fadc>baseline) fadc -= baseline;
			  else fadc = 0.;
			}
                        clusterCharge += fadc;
                        // as a matter of fact we should take the peak
                        // pos before FFT
                        // to get the list of tracks !!!
                        if(fadc > clusterPeakAmplitude) {
                            clusterPeakAmplitude = fadc;
                            //peakpos=Map()->GetHitIndex(idx,its);
                            Int_t shift = (int)(fTimeCorr/fTimeStep);
                            if(its>shift && its<(fMaxNofSamples-shift))
                                peakpos  = Map()->GetHitIndex(idx,its+shift);
                            else peakpos = Map()->GetHitIndex(idx,its);
                            if(peakpos<0) peakpos =Map()->GetHitIndex(idx,its);
                        } // end if
                        clusterTime += fadc*its;
                        if(fadc > 0) clusterMult++;
                        if(its == tstop) {
                            clusterTime /= (clusterCharge/fTimeStep);   // ns
                            if(clusterTime>fTimeCorr) clusterTime -=fTimeCorr;
                            //ns
                        } // end if
                    } // end for its

                    Double_t clusteranodePath = (clusterAnode - fNofAnodes/2)*
                                                 anodePitch;
		    Double_t clusterDriftPath = (Double_t)cal->GetDriftPath(clusterTime,clusteranodePath);
                    clusterDriftPath = fSddLength-clusterDriftPath;
                    if(clusterCharge <= 0.) break;
                    AliITSRawClusterSDD clust(j+1,//i
                                              clusterAnode,clusterTime,//ff
                                              clusterCharge, //f
                                              clusterPeakAmplitude, //f
                                              peakpos, //i
                                              0.,0.,clusterDriftPath,//fff
                                              clusteranodePath, //f
                                              clusterMult, //i
                                              0,0,0,0,0,0,0);//7*i
		    fDetTypeRec->AddCluster(1,&clust);
                    it = tstop;
                } // ilcl
                it++;
            } // while (samples)
        } // anodes
    } // detectors (2)

    for(i=0;i<fNofAnodes;i++) delete[] dfadc[i];
    delete [] dfadc;

    return;
}
//______________________________________________________________________
void AliITSClusterFinderSDD::Find1DClustersE(){
    // find 1D clusters
    // retrieve the parameters 
    Int_t fNofMaps = GetSeg()->Npz();
    Int_t fMaxNofSamples = GetSeg()->Npx();
    Int_t fNofAnodes = fNofMaps/2;
    Int_t dummy=0;
    Double_t fTimeStep = GetSeg()->Dpx( dummy );
    Double_t fSddLength = GetSeg()->Dx();
    Double_t anodePitch = GetSeg()->Dpz( dummy );
    AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
    Map()->ClearMap();
    Map()->SetThresholdArr( fCutAmplitude );
    Map()->FillMap2();

    AliITSresponseSDD* res  = (AliITSresponseSDD*)cal->GetResponse();
    const char *option=res->ZeroSuppOption();

    Int_t nClu = 0;
    //        cout << "Search  cluster... "<< endl;
    for( Int_t j=0; j<2; j++ ){
        for( Int_t k=0; k<fNofAnodes; k++ ){
            Int_t idx = j*fNofAnodes+k;	    
            Bool_t on = kFALSE;
            Int_t start = 0;
            Int_t nTsteps = 0;
            Double_t fmax = 0.;
            Int_t lmax = 0;
            Double_t charge = 0.;
            Double_t time = 0.;
            Double_t anode = k+0.5;
            Int_t peakpos = -1;
            for( Int_t l=0; l<fMaxNofSamples; l++ ){
                Double_t fadc = (Double_t)Map()->GetSignal( idx, l );
                if( fadc > 0.0 ){
                    if( on == kFALSE && l<fMaxNofSamples-4){
                        // star RawCluster (reset var.)
                        Double_t fadc1 = (Double_t)Map()->GetSignal( idx, l+1 );
                        if( fadc1 < fadc ) continue;
                        start = l;
                        fmax = 0.;
                        lmax = 0;
                        time = 0.;
                        charge = 0.; 
                        on = kTRUE; 
                        nTsteps = 0;
                    } // end if on...
                    nTsteps++ ;
		    if(!((strstr(option,"1D")) || (strstr(option,"2D")))){
		      Double_t baseline=GetResp(fModule)->GetBaseline(idx);
		      if( fadc > baseline ) fadc -= baseline;
		      else fadc=0.;
		    }
                    charge += fadc;
                    time += fadc*l;
                    if( fadc > fmax ){ 
                        fmax = fadc; 
                        lmax = l; 
                        Int_t shift = (Int_t)(fTimeCorr/fTimeStep + 0.5);
                        if( l > shift && l < (fMaxNofSamples-shift) )  
                            peakpos = Map()->GetHitIndex( idx, l+shift );
                        else
                            peakpos = Map()->GetHitIndex( idx, l );
                        if( peakpos < 0) peakpos = Map()->GetHitIndex(idx,l);
                    } // end if fadc
                }else{ // end fadc>0
                    if( on == kTRUE ){        
                        if( nTsteps > 2 ){
                            //  min # of timesteps for a RawCluster
                            // Found a RawCluster...
                            Int_t stop = l-1;
                            time /= (charge/fTimeStep);   // ns
                                // time = lmax*fTimeStep;   // ns
                            if( time > fTimeCorr ) time -= fTimeCorr;   // ns
                            Double_t anodePath =(anode-fNofAnodes/2)*anodePitch;
                            
			    Double_t driftPath = (Double_t)cal->GetDriftPath(time,anode);
                            driftPath = fSddLength-driftPath;
                            AliITSRawClusterSDD clust(j+1,anode,time,charge,
                                                      fmax, peakpos,0.,0.,
                                                      driftPath,anodePath,
                                                      nTsteps,start,stop,
                                                      start, stop, 1, k, k );
			    fDetTypeRec->AddCluster( 1, &clust );
                            if(AliDebugLevel()>=5) clust.PrintInfo();
                            nClu++;
                        } // end if nTsteps
                        on = kFALSE;
                    } // end if on==kTRUE
                } // end if fadc>0
            } // samples
        } // anodes
    } // wings
    AliDebug(3,Form("# Rawclusters %d",nClu));         
    return; 
}
//_______________________________________________________________________
Int_t AliITSClusterFinderSDD::SearchPeak(Double_t *spect,Int_t xdim,Int_t zdim,
                                         Int_t *peakX, Int_t *peakZ, 
                                         Double_t *peakAmp, Double_t minpeak ){
    // search peaks on a 2D cluster
    Int_t npeak = 0;    // # peaks
    Int_t i,j;
    // search peaks
    for( Int_t z=1; z<zdim-1; z++ ){
        for( Int_t x=1; x<xdim-2; x++ ){
            Double_t sxz = spect[x*zdim+z];
            Double_t sxz1 = spect[(x+1)*zdim+z];
            Double_t sxz2 = spect[(x-1)*zdim+z];
            // search a local max. in s[x,z]
            if( sxz < minpeak || sxz1 <= 0 || sxz2 <= 0 ) continue;
            if( sxz >= spect[(x+1)*zdim+z  ] && sxz >= spect[(x-1)*zdim+z  ] &&
                sxz >= spect[x*zdim    +z+1] && sxz >= spect[x*zdim    +z-1] &&
                sxz >= spect[(x+1)*zdim+z+1] && sxz >= spect[(x+1)*zdim+z-1] &&
                sxz >= spect[(x-1)*zdim+z+1] && sxz >= spect[(x-1)*zdim+z-1] ){
                // peak found
                peakX[npeak] = x;
                peakZ[npeak] = z;
                peakAmp[npeak] = sxz;
                npeak++;
            } // end if ....
        } // end for x
    } // end for z
    // search groups of peaks with same amplitude.
    Int_t *flag = new Int_t[npeak];
    for( i=0; i<npeak; i++ ) flag[i] = 0;
    for( i=0; i<npeak; i++ ){
        for( j=0; j<npeak; j++ ){
            if( i==j) continue;
            if( flag[j] > 0 ) continue;
            if( peakAmp[i] == peakAmp[j] && 
                TMath::Abs(peakX[i]-peakX[j])<=1 && 
                TMath::Abs(peakZ[i]-peakZ[j])<=1 ){
                if( flag[i] == 0) flag[i] = i+1;
                flag[j] = flag[i];
            } // end if ...
        } // end for j
    } // end for i
    // make average of peak groups        
    for( i=0; i<npeak; i++ ){
        Int_t nFlag = 1;
        if( flag[i] <= 0 ) continue;
        for( j=0; j<npeak; j++ ){
            if( i==j ) continue;
            if( flag[j] != flag[i] ) continue;
            peakX[i] += peakX[j];
            peakZ[i] += peakZ[j];
            nFlag++;
            npeak--;
            for( Int_t k=j; k<npeak; k++ ){
                peakX[k] = peakX[k+1];
                peakZ[k] = peakZ[k+1];
                peakAmp[k] = peakAmp[k+1];
                flag[k] = flag[k+1];
            } // end for k        
            j--;
        } // end for j
        if( nFlag > 1 ){
            peakX[i] /= nFlag;
            peakZ[i] /= nFlag;
        } // end fi nFlag
    } // end for i
    delete [] flag;
    return( npeak );
}
//______________________________________________________________________
void AliITSClusterFinderSDD::PeakFunc( Int_t xdim, Int_t zdim, Double_t *par,
                                       Double_t *spe, Double_t *integral){
    // function used to fit the clusters
    // par -> parameters..
    // par[0]  number of peaks.
    // for each peak i=1, ..., par[0]
    //                 par[i] = Ampl.
    //                 par[i+1] = xpos
    //                 par[i+2] = zpos
    //                 par[i+3] = tau
    //                 par[i+4] = sigma.
    Int_t electronics = GetResp(fModule)->GetElectronics(); // 1 = PASCAL, 2 = OLA
    const Int_t knParam = 5;
    Int_t npeak = (Int_t)par[0];

    memset( spe, 0, sizeof( Double_t )*zdim*xdim );

    Int_t k = 1;
    for( Int_t i=0; i<npeak; i++ ){
        if( integral != 0 ) integral[i] = 0.;
        Double_t sigmaA2 = par[k+4]*par[k+4]*2.;
        Double_t t2 = par[k+3];   // PASCAL
        if( electronics == 2 ) { t2 *= t2; t2 *= 2; } // OLA
        for( Int_t z=0; z<zdim; z++ ){
            for( Int_t x=0; x<xdim; x++ ){
                Double_t z2 = (z-par[k+2])*(z-par[k+2])/sigmaA2;
                Double_t x2 = 0.;
                Double_t signal = 0.;
                if( electronics == 1 ){ // PASCAL
                    x2 = (x-par[k+1]+t2)/t2;
                    signal = (x2>0.) ? par[k]*x2*exp(-x2+1.-z2) :0.0; // RCCR2
                //  signal =(x2>0.) ? par[k]*x2*x2*exp(-2*x2+2.-z2 ):0.0;//RCCR
                }else if( electronics == 2 ) { // OLA
                    x2 = (x-par[k+1])*(x-par[k+1])/t2;
                    signal = par[k]  * exp( -x2 - z2 );
                } else {
                    Warning("PeakFunc","Wrong SDD Electronics = %d",
                            electronics);
                    // exit( 1 );
                } // end if electronicx
                spe[x*zdim+z] += signal;
                if( integral != 0 ) integral[i] += signal;
            } // end for x
        } // end for z
        k += knParam;
    } // end for i
    return;
}
//__________________________________________________________________________
Double_t AliITSClusterFinderSDD::ChiSqr( Int_t xdim, Int_t zdim, Double_t *spe,
                                        Double_t *speFit ) const{
    // EVALUATES UNNORMALIZED CHI-SQUARED
    Double_t chi2 = 0.;
    for( Int_t z=0; z<zdim; z++ ){
        for( Int_t x=1; x<xdim-1; x++ ){
            Int_t index = x*zdim+z;
            Double_t tmp = spe[index] - speFit[index];
            chi2 += tmp*tmp;
        } // end for x
    } // end for z
    return( chi2 );
}
//_______________________________________________________________________
void AliITSClusterFinderSDD::Minim( Int_t xdim, Int_t zdim, Double_t *param,
                                    Double_t *prm0,Double_t *steprm,
                                    Double_t *chisqr,Double_t *spe,
                                    Double_t *speFit ){
    // 
    Int_t   k, nnn, mmm, i;
    Double_t p1, delta, d1, chisq1, p2, chisq2, t, p3, chisq3, a, b, p0, chisqt;
    const Int_t knParam = 5;
    Int_t npeak = (Int_t)param[0];
    for( k=1; k<(npeak*knParam+1); k++ ) prm0[k] = param[k];
    for( k=1; k<(npeak*knParam+1); k++ ){
        p1 = param[k];
        delta = steprm[k];
        d1 = delta;
        // ENSURE THAT STEP SIZE IS SENSIBLY LARGER THAN MACHINE ROUND OFF
        if( TMath::Abs( p1 ) > 1.0E-6 ) 
            if ( TMath::Abs( delta/p1 ) < 1.0E-4 ) delta = p1/1000;
            else  delta = (Double_t)1.0E-4;
        //  EVALUATE CHI-SQUARED AT FIRST TWO SEARCH POINTS
        PeakFunc( xdim, zdim, param, speFit );
        chisq1 = ChiSqr( xdim, zdim, spe, speFit );
        p2 = p1+delta;
        param[k] = p2;
        PeakFunc( xdim, zdim, param, speFit );
        chisq2 = ChiSqr( xdim, zdim, spe, speFit );
        if( chisq1 < chisq2 ){
            // REVERSE DIRECTION OF SEARCH IF CHI-SQUARED IS INCREASING
            delta = -delta;
            t = p1;
            p1 = p2;
            p2 = t;
            t = chisq1;
            chisq1 = chisq2;
            chisq2 = t;
        } // end if
        i = 1; nnn = 0;
        do {   // INCREMENT param(K) UNTIL CHI-SQUARED STARTS TO INCREASE
            nnn++;
            p3 = p2 + delta;
            mmm = nnn - (nnn/5)*5;  // multiplo de 5
            if( mmm == 0 ){
                d1 = delta;
                // INCREASE STEP SIZE IF STEPPING TOWARDS MINIMUM IS TOO SLOW 
                delta *= 5;
            } // end if
            param[k] = p3;
            // Constrain paramiters
            Int_t kpos = (k-1) % knParam;
            switch( kpos ){
            case 0 :
                if( param[k] <= 20 ) param[k] = fMinPeak;
                break;
            case 1 :
                if( TMath::Abs( param[k] - prm0[k] ) > 1.5 ) param[k] = prm0[k];
                break;
            case 2 :
                if( TMath::Abs( param[k] - prm0[k] ) > 1. ) param[k] = prm0[k];
                break;
            case 3 :
                if( param[k] < .5 ) param[k] = .5;        
                break;
            case 4 :
                if( param[k] < .288 ) param[k] = .288;// 1/sqrt(12) = 0.288
                if( param[k] > zdim*.5 ) param[k] = zdim*.5;
                break;
            }; // end switch
            PeakFunc( xdim, zdim, param, speFit );
            chisq3 = ChiSqr( xdim, zdim, spe, speFit );
            if( chisq3 < chisq2 && nnn < 50 ){
                p1 = p2;
                p2 = p3;
                chisq1 = chisq2;
                chisq2 = chisq3;
            }else i=0;
        } while( i );
        // FIND MINIMUM OF PARABOLA DEFINED BY LAST THREE POINTS
        a = chisq1*(p2-p3)+chisq2*(p3-p1)+chisq3*(p1-p2);
        b = chisq1*(p2*p2-p3*p3)+chisq2*(p3*p3-p1*p1)+chisq3*(p1*p1-p2*p2);
        if( a!=0 ) p0 = (Double_t)(0.5*b/a);
        else p0 = 10000;
        //--IN CASE OF NEARLY EQUAL CHI-SQUARED AND TOO SMALL STEP SIZE PREVENT
        //   ERRONEOUS EVALUATION OF PARABOLA MINIMUM
        //---NEXT TWO LINES CAN BE OMITTED FOR HIGHER PRECISION MACHINES
        //dp = (Double_t) max (TMath::Abs(p3-p2), TMath::Abs(p2-p1));
        //if( TMath::Abs( p2-p0 ) > dp ) p0 = p2;
        param[k] = p0;
        // Constrain paramiters
        Int_t kpos = (k-1) % knParam;
        switch( kpos ){
        case 0 :
            if( param[k] <= 20 ) param[k] = fMinPeak;   
            break;
        case 1 :
            if( TMath::Abs( param[k] - prm0[k] ) > 1.5 ) param[k] = prm0[k];
            break;
        case 2 :
            if( TMath::Abs( param[k] - prm0[k] ) > 1. ) param[k] = prm0[k];
            break;
        case 3 :
            if( param[k] < .5 ) param[k] = .5;        
            break;
        case 4 :
            if( param[k] < .288 ) param[k] = .288;  // 1/sqrt(12) = 0.288
            if( param[k] > zdim*.5 ) param[k] = zdim*.5;
            break;
        }; // end switch
        PeakFunc( xdim, zdim, param, speFit );
        chisqt = ChiSqr( xdim, zdim, spe, speFit );
        // DO NOT ALLOW ERRONEOUS INTERPOLATION
        if( chisqt <= *chisqr ) *chisqr = chisqt;
        else param[k] = prm0[k];
        // OPTIMIZE SEARCH STEP FOR EVENTUAL NEXT CALL OF MINIM
        steprm[k] = (param[k]-prm0[k])/5;
        if( steprm[k] >= d1 ) steprm[k] = d1/5;
    } // end for k
    // EVALUATE FIT AND CHI-SQUARED FOR OPTIMIZED PARAMETERS
    PeakFunc( xdim, zdim, param, speFit );
    *chisqr = ChiSqr( xdim, zdim, spe, speFit );
    return;
}
//_________________________________________________________________________
Int_t AliITSClusterFinderSDD::NoLinearFit( Int_t xdim, Int_t zdim, 
                                           Double_t *param, Double_t *spe, 
                                           Int_t *niter, Double_t *chir ){
    // fit method from Comput. Phys. Commun 46(1987) 149
    const Double_t kchilmt = 0.01;  //        relative accuracy           
    const Int_t   knel = 3;        //        for parabolic minimization  
    const Int_t   knstop = 50;     //        Max. iteration number          
    const Int_t   knParam = 5;
    Int_t npeak = (Int_t)param[0];
    // RETURN IF NUMBER OF DEGREES OF FREEDOM IS NOT POSITIVE 
    if( (xdim*zdim - npeak*knParam) <= 0 ) return( -1 );
    Double_t degFree = (xdim*zdim - npeak*knParam)-1;
    Int_t   n, k, iterNum = 0;
    Double_t *prm0 = new Double_t[npeak*knParam+1];
    Double_t *step = new Double_t[npeak*knParam+1];
    Double_t *schi = new Double_t[npeak*knParam+1]; 
    Double_t *sprm[3];
    sprm[0] = new Double_t[npeak*knParam+1];
    sprm[1] = new Double_t[npeak*knParam+1];
    sprm[2] = new Double_t[npeak*knParam+1];
    Double_t  chi0, chi1, reldif, a, b, prmin, dp;
    Double_t *speFit = new Double_t[ xdim*zdim ];
    PeakFunc( xdim, zdim, param, speFit );
    chi0 = ChiSqr( xdim, zdim, spe, speFit );
    chi1 = chi0;
    for( k=1; k<(npeak*knParam+1); k++) prm0[k] = param[k];
        for( k=1 ; k<(npeak*knParam+1); k+=knParam ){
            step[k] = param[k] / 20.0 ;
            step[k+1] = param[k+1] / 50.0;
            step[k+2] = param[k+2] / 50.0;                 
            step[k+3] = param[k+3] / 20.0;                 
            step[k+4] = param[k+4] / 20.0;                 
        } // end for k
    Int_t out = 0;
    do{
        iterNum++;
            chi0 = chi1;
            Minim( xdim, zdim, param, prm0, step, &chi1, spe, speFit );
            reldif = ( chi1 > 0 ) ? ((Double_t) TMath::Abs( chi1-chi0)/chi1 ) : 0;
        // EXIT conditions
        if( reldif < (float) kchilmt ){
            *chir  = (chi1>0) ? (float) TMath::Sqrt (chi1/degFree) :0;
            *niter = iterNum;
            out = 0;
            break;
        } // end if
        if( (reldif < (float)(5*kchilmt)) && (iterNum > knstop) ){
            *chir = (chi1>0) ?(float) TMath::Sqrt (chi1/degFree):0;
            *niter = iterNum;
            out = 0;
            break;
        } // end if
        if( iterNum > 5*knstop ){
            *chir  = (chi1>0) ?(float) TMath::Sqrt (chi1/degFree):0;
            *niter = iterNum;
            out = 1;
            break;
        } // end if
        if( iterNum <= knel ) continue;
        n = iterNum - (iterNum/knel)*knel; // EXTRAPOLATION LIMIT COUNTER N
        if( n > 3 || n == 0 ) continue;
        schi[n-1] = chi1;
        for( k=1; k<(npeak*knParam+1); k++ ) sprm[n-1][k] = param[k];
        if( n != 3 ) continue;
        // -EVALUATE EXTRAPOLATED VALUE OF EACH PARAMETER BY FINDING MINIMUM OF
        //    PARABOLA DEFINED BY LAST THREE CALLS OF MINIM
        for( k=1; k<(npeak*knParam+1); k++ ){
            Double_t tmp0 = sprm[0][k];
            Double_t tmp1 = sprm[1][k];
            Double_t tmp2 = sprm[2][k];
            a  = schi[0]*(tmp1-tmp2) + schi[1]*(tmp2-tmp0);
            a += (schi[2]*(tmp0-tmp1));
            b  = schi[0]*(tmp1*tmp1-tmp2*tmp2);
            b += (schi[1]*(tmp2*tmp2-tmp0*tmp0)+(schi[2]*
                                             (tmp0*tmp0-tmp1*tmp1)));
            if ((double)a < 1.0E-6) prmin = 0;
            else prmin = (float) (0.5*b/a);
            dp = 5*(tmp2-tmp0);
            if( TMath::Abs(prmin-tmp2) > TMath::Abs(dp) ) prmin = tmp2+dp;
            param[k] = prmin;
            step[k]  = dp/10; // OPTIMIZE SEARCH STEP
        } // end for k
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

//______________________________________________________________________
void AliITSClusterFinderSDD::ResolveClusters(){
    // The function to resolve clusters if the clusters overlapping exists
    Int_t i;
    // get number of clusters for this module
    Int_t nofClusters = NClusters();
    nofClusters -= fNclusters;
    Int_t fNofMaps = GetSeg()->Npz();
    Int_t fNofAnodes = fNofMaps/2;
    //Int_t fMaxNofSamples = GetSeg()->Npx();
    Int_t dummy=0;
    Double_t fTimeStep = GetSeg()->Dpx( dummy );
    Double_t fSddLength = GetSeg()->Dx();
    Double_t anodePitch = GetSeg()->Dpz( dummy );
    Int_t electronics =GetResp(fModule)->GetElectronics(); // 1 = PASCAL, 2 = OLA
    AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
    AliITSresponseSDD* res  = (AliITSresponseSDD*)cal->GetResponse();
    const char *option=res->ZeroSuppOption();
    

    for( Int_t j=0; j<nofClusters; j++ ){ 
        // get cluster information
        AliITSRawClusterSDD *clusterJ=(AliITSRawClusterSDD*) Cluster(j);
        Int_t astart = clusterJ->Astart();
        Int_t astop = clusterJ->Astop();
        Int_t tstart = clusterJ->Tstartf();
        Int_t tstop = clusterJ->Tstopf();
        Int_t wing = (Int_t)clusterJ->W();
        if( wing == 2 ){
            astart += fNofAnodes; 
            astop  += fNofAnodes;
        } // end if 
        Int_t xdim = tstop-tstart+3;
        Int_t zdim = astop-astart+3;
        if( xdim > 50 || zdim > 30 ) { 
            Warning("ResolveClusters","xdim: %d , zdim: %d ",xdim,zdim);
            continue;
        }
        Double_t *sp = new Double_t[ xdim*zdim+1 ];
        memset( sp, 0, sizeof(Double_t)*(xdim*zdim+1) );
        
        // make a local map from cluster region
        for( Int_t ianode=astart; ianode<=astop; ianode++ ){
            for( Int_t itime=tstart; itime<=tstop; itime++ ){
                Double_t fadc = Map()->GetSignal( ianode, itime );	       
		if(!((strstr(option,"1D")) || (strstr(option,"2D")))){
		  Double_t baseline=GetResp(fModule)->GetBaseline(ianode);
		  if( fadc > baseline ) fadc -= (Double_t)baseline;
		  else fadc = 0.;
		}
                Int_t index = (itime-tstart+1)*zdim+(ianode-astart+1);
                sp[index] = fadc;
            } // time loop
        } // anode loop
        
        // search peaks on cluster
        const Int_t kNp = 150;
        Int_t peakX1[kNp];
        Int_t peakZ1[kNp];
        Double_t peakAmp1[kNp];
        Int_t npeak = SearchPeak(sp,xdim,zdim,peakX1,peakZ1,peakAmp1,fMinPeak);

        // if multiple peaks, split cluster
        if( npeak >= 1 ){
            //        cout << "npeak " << npeak << endl;
            //        clusterJ->PrintInfo();
            Double_t *par = new Double_t[npeak*5+1];
            par[0] = (Double_t)npeak;                
            // Initial parameters in cell dimentions
            Int_t k1 = 1;
            for( i=0; i<npeak; i++ ){
                par[k1] = peakAmp1[i];
                par[k1+1] = peakX1[i]; // local time pos. [timebin]
                par[k1+2] = peakZ1[i]; // local anode pos. [anodepitch]
                if( electronics == 1 ) par[k1+3] = 2.; // PASCAL
                else if(electronics==2) par[k1+3] = 0.7;//tau [timebin] OLA 
                par[k1+4] = .4;    // sigma        [anodepich]
                k1 += 5;
            } // end for i                        
            Int_t niter;
            Double_t chir;                        
            NoLinearFit( xdim, zdim, par, sp, &niter, &chir );
            Double_t peakX[kNp];
            Double_t peakZ[kNp];
            Double_t sigma[kNp];
            Double_t tau[kNp];
            Double_t peakAmp[kNp];
            Double_t integral[kNp];
            //get integrals => charge for each peak
            PeakFunc( xdim, zdim, par, sp, integral );
            k1 = 1;
            for( i=0; i<npeak; i++ ){
                peakAmp[i] = par[k1];
                peakX[i]   = par[k1+1];
                peakZ[i]   = par[k1+2];
                tau[i]     = par[k1+3];
                sigma[i]   = par[k1+4];
                k1+=5;
            } // end for i
            // calculate parameter for new clusters
            for( i=0; i<npeak; i++ ){
                AliITSRawClusterSDD clusterI( *clusterJ );
            
                Int_t newAnode = peakZ1[i]-1 + astart;

            //    Int_t newiTime = peakX1[i]-1 + tstart;
            //    Int_t shift = (Int_t)(fTimeCorr/fTimeStep + 0.5);
            //    if( newiTime > shift && newiTime < (fMaxNofSamples-shift) ) 
            //        shift = 0;
            //    Int_t peakpos = Map()->GetHitIndex(newAnode,newiTime+shift );
            //    clusterI.SetPeakPos( peakpos );
	    
                clusterI.SetPeakAmpl( peakAmp1[i] );
                Double_t newAnodef = peakZ[i] - 0.5 + astart;
                Double_t newiTimef = peakX[i] - 1 + tstart;
                if( wing == 2 ) newAnodef -= fNofAnodes; 
                Double_t anodePath = (newAnodef - fNofAnodes/2)*anodePitch;
                newiTimef *= fTimeStep;
                if( newiTimef > fTimeCorr ) newiTimef -= fTimeCorr;
                if( electronics == 1 ){
                //    newiTimef *= 0.999438;    // PASCAL
                //    newiTimef += (6./fDriftSpeed - newiTimef/3000.);
                }else if( electronics == 2 )
                    newiTimef *= 0.99714;    // OLA
                    
                Int_t timeBin = (Int_t)(newiTimef/fTimeStep+0.5);    
                Int_t peakpos = Map()->GetHitIndex( newAnode, timeBin );
                if( peakpos < 0 ) { 
                    for( Int_t ii=0; ii<3; ii++ ) {
                        peakpos = Map()->GetHitIndex( newAnode, timeBin+ii );
                        if( peakpos > 0 ) break;
                        peakpos = Map()->GetHitIndex( newAnode, timeBin-ii );
                        if( peakpos > 0 ) break;
                    }
                }
                
                if( peakpos < 0 ) { 
                    //Warning("ResolveClusters",
                    //        "Digit not found for cluster");
                    //if(AliDebugLevel()>=3) clusterI.PrintInfo(); 
                   continue;
                }
                clusterI.SetPeakPos( peakpos ); 
		Float_t dp = cal->GetDriftPath(newiTimef,newAnodef);
		Double_t driftPath = fSddLength - (Double_t)dp;
                Double_t sign = ( wing == 1 ) ? -1. : 1.;
		Double_t xcoord = driftPath*sign * 0.0001;
		Double_t zcoord = anodePath * 0.0001;
		CorrectPosition(zcoord,xcoord);
                clusterI.SetX( xcoord );        
                clusterI.SetZ( zcoord );
                clusterI.SetAnode( newAnodef );
                clusterI.SetTime( newiTimef );
                clusterI.SetAsigma( sigma[i]*anodePitch );
                clusterI.SetTsigma( tau[i]*fTimeStep );
                clusterI.SetQ( integral[i] );
                
		fDetTypeRec->AddCluster( 1, &clusterI );
            } // end for i
            Clusters()->RemoveAt( j );
            delete [] par;
        } else {  // something odd
            Warning( "ResolveClusters",
                     "--- Peak not found!!!!  minpeak=%d ,cluster peak= %f"
                     " , module= %d",
                     fMinPeak, clusterJ->PeakAmpl(),GetModule()); 
            clusterJ->PrintInfo();
            Warning( "ResolveClusters"," xdim= %d zdim= %d", xdim-2, zdim-2 );
        }
        delete [] sp;
    } // cluster loop
    Clusters()->Compress();
//    Map()->ClearMap(); 
}
//________________________________________________________________________
void  AliITSClusterFinderSDD::GroupClusters(){
    // group clusters
    Int_t dummy=0;
    Double_t fTimeStep = GetSeg()->Dpx(dummy);
    // get number of clusters for this module
    Int_t nofClusters = NClusters();
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
            clusterI = (AliITSRawClusterSDD*) Cluster(i);
            clusterJ = (AliITSRawClusterSDD*) Cluster(j);
            // 1.3 good
            if(clusterI->T() < fTimeStep*60) fDAnode = 4.2;  // TB 3.2  
            if(clusterI->T() < fTimeStep*10) fDAnode = 1.5;  // TB 1.
            Bool_t pair = clusterI->Brother(clusterJ,fDAnode,fDTime);
            if(!pair) continue;
            if(AliDebugLevel()>=4){
                clusterI->PrintInfo();
                clusterJ->PrintInfo();
            } // end if AliDebugLevel
            clusterI->Add(clusterJ);
            label[j] = 1;
            Clusters()->RemoveAt(j);
            j=i; // <- Ernesto
        } // J clusters  
        label[i] = 1;
    } // I clusters
    Clusters()->Compress();

    delete [] label;
    return;
}
//________________________________________________________________________
void AliITSClusterFinderSDD::SelectClusters(){
    // get number of clusters for this module
    Int_t nofClusters = NClusters();

    nofClusters -= fNclusters;
    Int_t i;
    for(i=0; i<nofClusters; i++) { 
        AliITSRawClusterSDD *clusterI =(AliITSRawClusterSDD*) Cluster(i);
        Int_t rmflg = 0;
        Double_t wy = 0.;
        if(clusterI->Anodes() != 0.) {
            wy = ((Double_t) clusterI->Samples())/clusterI->Anodes();
        } // end if
        Int_t amp = (Int_t) clusterI->PeakAmpl();
        Int_t cha = (Int_t) clusterI->Q();
        if(amp < fMinPeak) rmflg = 1;  
        if(cha < fMinCharge) rmflg = 1;
        if(wy < fMinNCells) rmflg = 1;
        //if(wy > fMaxNCells) rmflg = 1;
        if(rmflg) Clusters()->RemoveAt(i);
    } // I clusters
    Clusters()->Compress();
    return;
}

//______________________________________________________________________
void AliITSClusterFinderSDD::GetRecPoints(AliITSCalibrationSDD* cal){
    // get rec points
  
    // get number of clusters for this module
    Int_t nofClusters = NClusters();
    nofClusters -= fNclusters;
    const Double_t kconvGeV = 1.e-6; // GeV -> KeV
    const Double_t kconv = 1.0e-4; 
    const Double_t kcmToMicrons = 10000.;
    const Double_t kRMSx = 38.0*kconv; // microns->cm ITS TDR Table 1.3
    const Double_t kRMSz = 28.0*kconv; // microns->cm ITS TDR Table 1.3
    Int_t nAnodes=GetSeg()->NpzHalf();
    Int_t i;
    Int_t ix, iz, idx=-1;
    AliITSdigitSDD *dig=0;
    Int_t ndigits=NDigits();

    Int_t lay,lad,det;
    fDetTypeRec->GetITSgeom()->GetModuleId(fModule,lay,lad,det);
    Int_t ind=(lad-1)*fDetTypeRec->GetITSgeom()->GetNdetectors(lay)+(det-1);
    Int_t lyr=(lay-1);


    for(i=0; i<nofClusters; i++) { 
        AliITSRawClusterSDD *clusterI = (AliITSRawClusterSDD*)Cluster(i);
        if(!clusterI) Error("SDD: GetRecPoints","i clusterI ",i,clusterI);
        if(clusterI) idx=clusterI->PeakPos();
        if(idx>ndigits) Error("SDD: GetRecPoints","idx ndigits",idx,ndigits);
        // try peak neighbours - to be done 
        if(idx&&idx<= ndigits) dig =(AliITSdigitSDD*)GetDigit(idx);
        if(!dig) {
            // try cog
	    Float_t xMicrons=clusterI->X()*kcmToMicrons;
	    Float_t zMicrons=clusterI->Z()*kcmToMicrons;
	    Float_t zAnode=zMicrons/GetSeg()->Dpz(0)+nAnodes/2;
	    Float_t driftSpeed=cal->GetDriftSpeedAtAnode(zAnode);	    
	    Float_t driftPath=GetSeg()->Dx()-TMath::Abs(xMicrons);
	    ix=1+(Int_t)(driftPath/driftSpeed/GetSeg()->Dpx(0));
	    iz=1+(Int_t)zAnode;
	    if(xMicrons>0) iz+=nAnodes;
            dig = (AliITSdigitSDD*)Map()->GetHit(iz-1,ix-1);
            // if null try neighbours
            if (!dig) dig = (AliITSdigitSDD*)Map()->GetHit(iz-1,ix); 
            if (!dig) dig = (AliITSdigitSDD*)Map()->GetHit(iz-1,ix+1); 
            if (!dig) printf("SDD: cannot assign the track number!\n");
        } //  end if !dig

	Int_t lab[4] = {-3141593,-3141593,-3141593,ind};
	if (dig) {
	  lab[0] = dig->GetTrack(0);
	  lab[1] = dig->GetTrack(1);
	  lab[2] = dig->GetTrack(2);
	}
	Float_t hit[5] = {clusterI->X(),clusterI->Z(),kRMSx*kRMSx,kRMSz*kRMSz,clusterI->Q()};
	Int_t info[3] = {0,0,lyr};

        AliITSRecPoint rnew(lab,hit,info,kTRUE);
        rnew.SetdEdX(kconvGeV*clusterI->Q());

	fDetTypeRec->AddRecPoint(rnew);        
    } // I clusters
//    Map()->ClearMap();
}
//______________________________________________________________________
void AliITSClusterFinderSDD::FindRawClusters(Int_t mod){
    // find raw clusters
    
    SetModule(mod);
    SetCutAmplitude(mod);
    AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(mod);
    Int_t nanodes=GetSeg()->Npz();
    Int_t noise=0;
    for(Int_t i=0;i<nanodes;i++){
      noise+=(Int_t)cal->GetNoiseAfterElectronics(i);
    }    
    SetMinPeak((noise/nanodes)*5);
    Find1DClustersE();
    GroupClusters();
    SelectClusters();
    ResolveClusters();
    GetRecPoints(cal);
}
//_______________________________________________________________________
void AliITSClusterFinderSDD::PrintStatus() const{
    // Print SDD cluster finder Parameters

    cout << "**************************************************" << endl;
    cout << " Silicon Drift Detector Cluster Finder Parameters " << endl;
    cout << "**************************************************" << endl;
    cout << "Number of Clusters: " << fNclusters << endl;
    cout << "Anode Tolerance: " << fDAnode << endl;
    cout << "Time  Tolerance: " << fDTime << endl;
    cout << "Time  correction (electronics): " << fTimeCorr << endl;
    cout << "Cut Amplitude (threshold): " << fCutAmplitude[0] << endl;
    cout << "Minimum Amplitude: " << fMinPeak << endl;
    cout << "Minimum Charge: " << fMinCharge << endl;
    cout << "Minimum number of cells/clusters: " << fMinNCells << endl;
    cout << "Maximum number of cells/clusters: " << fMaxNCells << endl;
    cout << "**************************************************" << endl;
}

//_________________________________________________________________________
void AliITSClusterFinderSDD::CorrectPosition(Double_t &z, Double_t&y){
  //correction of coordinates using the maps stored in the DB

  AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
  static const Int_t nbint = cal->GetMapTimeNBin();
  static const Int_t nbina = cal->Chips()*cal->Channels();
  Float_t stepa = (GetSeg()->Dpz(0))/10000.; //anode pitch in cm
  Float_t stept = (GetSeg()->Dx()/cal->GetMapTimeNBin()/2.)/10.;

  Int_t bint = TMath::Abs((Int_t)(y/stept));
  if(y>=0) bint+=(Int_t)(nbint/2.);
  if(bint>nbint) AliError("Wrong bin number!");

  Int_t bina = TMath::Abs((Int_t)(z/stepa));
  if(z>=0) bina+=(Int_t)(nbina/2.);
  if(bina>nbina) AliError("Wrong bin number!");

  Double_t devz = (Double_t)cal->GetMapACell(bina,bint)/10000.;
  Double_t devx = (Double_t)cal->GetMapTCell(bina,bint)/10000.;
  z+=devz;
  y+=devx;

}
