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
  $Log$
  Revision 1.36  2004/01/27 16:12:03  masera
  Coding conventions for AliITSdigitXXX classes and AliITSTrackerV1

  Revision 1.35  2003/11/10 16:33:50  masera
  Changes to obey our coding conventions

  Revision 1.34  2003/09/11 13:48:52  masera
  Data members of AliITSdigit classes defined as protected (They were public)

  Revision 1.33  2003/07/21 14:20:51  masera
  Fix to track labes in SDD Rec-points

  Revision 1.31.2.1  2003/07/16 13:18:04  masera
  Proper fix to track labels associated to SDD rec-points

  Revision 1.31  2003/05/19 14:44:41  masera
  Fix to track labels associated to SDD rec-points

  Revision 1.30  2003/03/03 16:34:35  masera
  Corrections to comply with coding conventions

  Revision 1.29  2002/10/25 18:54:22  barbera
  Various improvements and updates from B.S.Nilsen and T. Virgili

  Revision 1.28  2002/10/22 14:45:29  alibrary
  Introducing Riostream.h

  Revision 1.27  2002/10/14 14:57:00  hristov
  Merging the VirtualMC branch to the main development branch (HEAD)

  Revision 1.23.4.2  2002/10/14 13:14:07  hristov
  Updating VirtualMC to v3-09-02

  Revision 1.26  2002/09/09 17:23:28  nilsen
  Minor changes in support of changes to AliITSdigitS?D class'.

  Revision 1.25  2002/05/10 22:29:40  nilsen
  Change my Massimo Masera in the default constructor to bring things into
  compliance.

  Revision 1.24  2002/04/24 22:02:31  nilsen
  New SDigits and Digits routines, and related changes,  (including new
  noise values).

 */
// 
//  Cluster finder 
//  for Silicon
//  Drift Detector
//
#include <Riostream.h>

#include <TMath.h>
#include <math.h>

#include "AliITSClusterFinderSDD.h"
#include "AliITSMapA1.h"
#include "AliITS.h"
#include "AliITSdigitSDD.h"
#include "AliITSRawClusterSDD.h"
#include "AliITSRecPoint.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliRun.h"

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
AliITSClusterFinderSDD::AliITSClusterFinderSDD(AliITSsegmentation *seg,
                                               AliITSresponse *response,
                                               TClonesArray *digits,
                                               TClonesArray *recp):
AliITSClusterFinder(seg,response),
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
    SetCutAmplitude();
    SetDAnode();
    SetDTime();
    SetMinPeak((Int_t)(((AliITSresponseSDD*)GetResp())->
                       GetNoiseAfterElectronics()*5));
    //    SetMinPeak();
    SetMinNCells();
    SetMaxNCells();
    SetTimeCorr();
    SetMinCharge();
    SetMap(new AliITSMapA1(GetSeg(),Digits(),fCutAmplitude));
}
//______________________________________________________________________
void AliITSClusterFinderSDD::SetCutAmplitude(Double_t nsigma){
    // set the signal threshold for cluster finder
    Double_t baseline,noise,noiseAfterEl;

    GetResp()->GetNoiseParam(noise,baseline);
    noiseAfterEl = ((AliITSresponseSDD*)GetResp())->GetNoiseAfterElectronics();
    fCutAmplitude = (Int_t)((baseline + nsigma*noiseAfterEl));
}
//______________________________________________________________________
void AliITSClusterFinderSDD::Find1DClusters(){
    // find 1D clusters
    static AliITS *iTS = (AliITS*)gAlice->GetModule("ITS");
  
    // retrieve the parameters 
    Int_t fNofMaps       = GetSeg()->Npz();
    Int_t fMaxNofSamples = GetSeg()->Npx();
    Int_t fNofAnodes     = fNofMaps/2;
    Int_t dummy          = 0;
    Double_t fTimeStep    = GetSeg()->Dpx(dummy);
    Double_t fSddLength   = GetSeg()->Dx();
    Double_t fDriftSpeed  = GetResp()->DriftSpeed();  
    Double_t anodePitch   = GetSeg()->Dpz(dummy);

    // map the signal
    Map()->ClearMap();
    Map()->SetThreshold(fCutAmplitude);
    Map()->FillMap();
  
    Double_t noise;
    Double_t baseline;
    GetResp()->GetNoiseParam(noise,baseline);
  
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
            if(GetDebug(5)) cout<<"Anode: "<<k+1<<", Wing: "<<j+1<< endl;
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
                    if(fadc > (float)fCutAmplitude)lthrt++; 
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
                           (fadc > fCutAmplitude)) {
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
                    Double_t n, baseline;
                    GetResp()->GetNoiseParam(n,baseline);
                    for(its=tstart; its<=tstop; its++) {
                        fadc=(float)Map()->GetSignal(idx,its);
                        if(fadc>baseline) fadc -= baseline;
                        else fadc = 0.;
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
                    Double_t clusterDriftPath = clusterTime*fDriftSpeed;
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
                    iTS->AddCluster(1,&clust);
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
    static AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
    // retrieve the parameters 
    Int_t fNofMaps = GetSeg()->Npz();
    Int_t fMaxNofSamples = GetSeg()->Npx();
    Int_t fNofAnodes = fNofMaps/2;
    Int_t dummy=0;
    Double_t fTimeStep = GetSeg()->Dpx( dummy );
    Double_t fSddLength = GetSeg()->Dx();
    Double_t fDriftSpeed = GetResp()->DriftSpeed();
    Double_t anodePitch = GetSeg()->Dpz( dummy );
    Double_t n, baseline;
    GetResp()->GetNoiseParam( n, baseline );
    // map the signal
    Map()->ClearMap();
    Map()->SetThreshold( fCutAmplitude );
    Map()->FillMap();
    
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
                    if( fadc > baseline ) fadc -= baseline;
                    else fadc=0.;
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
                            Double_t driftPath = time*fDriftSpeed;
                            driftPath = fSddLength-driftPath;
                            AliITSRawClusterSDD clust(j+1,anode,time,charge,
                                                      fmax, peakpos,0.,0.,
                                                      driftPath,anodePath,
                                                      nTsteps,start,stop,
                                                      start, stop, 1, k, k );
                            iTS->AddCluster( 1, &clust );
                            if(GetDebug(5)) clust.PrintInfo();
                            nClu++;
                        } // end if nTsteps
                        on = kFALSE;
                    } // end if on==kTRUE
                } // end if fadc>0
            } // samples
        } // anodes
    } // wings
    if(GetDebug(3)) cout << "# Rawclusters " << nClu << endl;         
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
    Int_t electronics = GetResp()->Electronics(); // 1 = PASCAL, 2 = OLA
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
        if( fabs( p1 ) > 1.0E-6 ) 
            if ( fabs( delta/p1 ) < 1.0E-4 ) delta = p1/1000;
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
                if( fabs( param[k] - prm0[k] ) > 1.5 ) param[k] = prm0[k];
                break;
            case 2 :
                if( fabs( param[k] - prm0[k] ) > 1. ) param[k] = prm0[k];
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
        //dp = (Double_t) max (fabs(p3-p2), fabs(p2-p1));
        //if( fabs( p2-p0 ) > dp ) p0 = p2;
        param[k] = p0;
        // Constrain paramiters
        Int_t kpos = (k-1) % knParam;
        switch( kpos ){
        case 0 :
            if( param[k] <= 20 ) param[k] = fMinPeak;   
            break;
        case 1 :
            if( fabs( param[k] - prm0[k] ) > 1.5 ) param[k] = prm0[k];
            break;
        case 2 :
            if( fabs( param[k] - prm0[k] ) > 1. ) param[k] = prm0[k];
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
            reldif = ( chi1 > 0 ) ? ((Double_t) fabs( chi1-chi0)/chi1 ) : 0;
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
            if( fabs(prmin-tmp2) > fabs(dp) ) prmin = tmp2+dp;
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
void AliITSClusterFinderSDD::ResolveClustersE(){
    // The function to resolve clusters if the clusters overlapping exists
    Int_t i;
    static AliITS *iTS = (AliITS*)gAlice->GetModule( "ITS" );
    // get number of clusters for this module
    Int_t nofClusters = NClusters();
    nofClusters -= fNclusters;
    Int_t fNofMaps = GetSeg()->Npz();
    Int_t fNofAnodes = fNofMaps/2;
    //Int_t fMaxNofSamples = GetSeg()->Npx();
    Int_t dummy=0;
    Double_t fTimeStep = GetSeg()->Dpx( dummy );
    Double_t fSddLength = GetSeg()->Dx();
    Double_t fDriftSpeed = GetResp()->DriftSpeed();
    Double_t anodePitch = GetSeg()->Dpz( dummy );
    Double_t n, baseline;
    GetResp()->GetNoiseParam( n, baseline );
    Int_t electronics = GetResp()->Electronics(); // 1 = PASCAL, 2 = OLA

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
            Warning("ResolveClustersE","xdim: %d , zdim: %d ",xdim,zdim);
            continue;
        }
        Double_t *sp = new Double_t[ xdim*zdim+1 ];
        memset( sp, 0, sizeof(Double_t)*(xdim*zdim+1) );
        
        // make a local map from cluster region
        for( Int_t ianode=astart; ianode<=astop; ianode++ ){
            for( Int_t itime=tstart; itime<=tstop; itime++ ){
                Double_t fadc = Map()->GetSignal( ianode, itime );
                if( fadc > baseline ) fadc -= (Double_t)baseline;
                else fadc = 0.;
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
                    //Warning("ResolveClustersE",
                    //        "Digit not found for cluster");
                    //if(GetDebug(3)) clusterI.PrintInfo(); 
                   continue;
                }
                clusterI.SetPeakPos( peakpos );    
                Double_t driftPath = fSddLength - newiTimef * fDriftSpeed;
                Double_t sign = ( wing == 1 ) ? -1. : 1.;
                clusterI.SetX( driftPath*sign * 0.0001 );        
                clusterI.SetZ( anodePath * 0.0001 );
                clusterI.SetAnode( newAnodef );
                clusterI.SetTime( newiTimef );
                clusterI.SetAsigma( sigma[i]*anodePitch );
                clusterI.SetTsigma( tau[i]*fTimeStep );
                clusterI.SetQ( integral[i] );
                
                iTS->AddCluster( 1, &clusterI );
            } // end for i
            Clusters()->RemoveAt( j );
            delete [] par;
        } else {  // something odd
            Warning( "ResolveClustersE",
                     "--- Peak not found!!!!  minpeak=%d ,cluster peak= %f"
                     " , module= %d",
                     fMinPeak, clusterJ->PeakAmpl(),GetModule()); 
            clusterJ->PrintInfo();
            Warning( "ResolveClustersE"," xdim= %d zdim= %d", xdim-2, zdim-2 );
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
            if(GetDebug(4)){
                clusterI->PrintInfo();
                clusterJ->PrintInfo();
            } // end if GetDebug
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
//__________________________________________________________________________
void AliITSClusterFinderSDD::ResolveClusters(){
    // The function to resolve clusters if the clusters overlapping exists
/*    AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
    // get number of clusters for this module
    Int_t nofClusters = NClusters();
    nofClusters -= fNclusters;
    //cout<<"Resolve Cl: nofClusters, fNclusters ="<<nofClusters<<","
    // <<fNclusters<<endl;
    Int_t fNofMaps = GetSeg()->Npz();
    Int_t fNofAnodes = fNofMaps/2;
    Int_t dummy=0;
    Double_t fTimeStep = GetSeg()->Dpx(dummy);
    Double_t fSddLength = GetSeg()->Dx();
    Double_t fDriftSpeed = GetResp()->DriftSpeed();
    Double_t anodePitch = GetSeg()->Dpz(dummy);
    Double_t n, baseline;
    GetResp()->GetNoiseParam(n,baseline);
    Double_t dzz_1A = anodePitch * anodePitch / 12;
    // fill Map of signals
    Map()->FillMap(); 
    Int_t j,i,ii,ianode,anode,itime;
    Int_t wing,astart,astop,tstart,tstop,nanode;
    Double_t fadc,ClusterTime;
    Double_t q[400],x[400],z[400]; // digit charges and coordinates
    for(j=0; j<nofClusters; j++) { 
        AliITSRawClusterSDD *clusterJ=(AliITSRawClusterSDD*) Cluster(j);
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
        }  // end if
        // cout<<"astart,astop,tstart,tstop ="<<astart<<","<<astop<<","
        //      <<tstart<<","<<tstop<<endl;
        // clear the digit arrays
        for(ii=0; ii<400; ii++) { 
            q[ii] = 0.; 
            x[ii] = 0.;
            z[ii] = 0.;
        } // end for ii

        for(ianode=astart; ianode<=astop; ianode++) { 
            for(itime=tstart; itime<=tstop; itime++) { 
                fadc=Map()->GetSignal(ianode,itime);
                if(fadc>baseline) {
                    fadc-=(Double_t)baseline;
                    q[ndigits] = fadc*(fTimeStep/160);  // KeV
                    anode = ianode;
                    if(wing == 2) anode -= fNofAnodes;
                    z[ndigits] = (anode + 0.5 - fNofAnodes/2)*anodePitch;
                    ClusterTime = itime*fTimeStep;
                    if(ClusterTime > fTimeCorr) ClusterTime -= fTimeCorr;// ns
                    x[ndigits] = fSddLength - ClusterTime*fDriftSpeed;
                    if(wing == 1) x[ndigits] *= (-1);
                    // cout<<"ianode,itime,fadc ="<<ianode<<","<<itime<<","
                    //     <<fadc<<endl;
                    // cout<<"wing,anode,ndigits,charge ="<<wing<<","
                    //      <<anode<<","<<ndigits<<","<<q[ndigits]<<endl;
                    ndigits++;
                    continue;
                } //  end if
                fadc=0;
                //              cout<<"fadc=0, ndigits ="<<ndigits<<endl;
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
        } // end for ii
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
        if( nanode == 1 ){
            dzz = dzz_1A; // for one anode cluster dzz = anode**2/12
            scl = TMath::Sqrt( 7.2/(-0.57*xm*1.e-3+71.8) );
        } // end if
        if( nanode == 2 ){
            scl = TMath::Sqrt( (-0.18*xm*1.e-3+21.3)/(-0.57*xm*1.e-3+71.8) );
        } // end if
        if( nanode == 3 ){
            scl = TMath::Sqrt( (-0.5*xm*1.e-3+34.5)/(-0.57*xm*1.e-3+71.8) );
        } // end if
        if( nanode > 3 ){
            scl = TMath::Sqrt( (1.3*xm*1.e-3+49.)/(-0.57*xm*1.e-3+71.8) );
        } // end if
        //   cout<<"1 microns: zm,dzz,xm,dxx,dxz,qq ="<<zm<<","<<dzz<<","
        //  <<xm<<","<<dxx<<","<<dxz<<","<<qq<<endl;
        //  old Boris.........
        //  tmp=29730. - 585.*fabs(xm/1000.); 
        //  scl=TMath::Sqrt(tmp/130000.);
   
        xm *= scl;
        xx *= scl*scl;
        xz *= scl;

        dxx = xx - xm*xm;
        //   dzz = zz - zm*zm;
        dxz = xz - xm*zm;
        //   cout<<"microns: zm,dzz,xm,dxx,xz,dxz,qq ="<<zm<<","<<dzz<<","
        // <<xm<<","<<dxx<<","<<xz<<","<<dxz<<","<<qq<<endl;
        // if(dzz < 7200.) dzz=7200.;//for one anode cluster dzz = anode**2/12
  
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
                tga = (dxz<0.) ? tmp-TMath::Sqrt(tmp*tmp+1) : 
                                                   tmp+TMath::Sqrt(tmp*tmp+1);
            } // end if dxz
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
            //   cout<<"cm: zm,dzz,xm,dxx,xz,dxz,qq ="<<zm<<","<<dzz<<","
            //  <<xm<<","<<dxx<<","<<xz<<","<<dxz<<","<<qq<<endl;
            for (i=0; i<nbins; i++) {     
                x[i] = x[i] *= scl;
                x[i] = x[i] *= 1.e-4;
                z[i] = z[i] *= 1.e-4;
            } // end for i
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
                    } // end if
                    xxx += tmp*tmp*tmp*q[i];
                } // end for i
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
                    Double_t f2=-(r*(sigma2+x1*x1)+(1-r)*(sigma2+x2*x2)- xx);
                    Double_t f3=-(r*x1*(3*sigma2+x1*x1)+(1-r)*x2*
                                                         (3*sigma2+x2*x2)-xxx);
                    Double_t d=c11*c22*c33+c21*c32*c13+c12*c23*c31-
                                       c31*c22*c13 - c21*c12*c33 - c32*c23*c11;
                    if (d==0.) {
                        cout<<"*********** d=0 ***********\n";
                        break;
                    } // end if
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
                    qfit[0]=a1; xfit[0]=x1*cosa - z12*sina; zfit[0]=x1*sina + 
                                                                z12*cosa;
                    qfit[1]=a2; xfit[1]=x2*cosa - z12*sina; zfit[1]=x2*sina + 
                                                                z12*cosa;
                    nfhits=2;
                    break; // Ok !
                } // end for i
                if (i==33) cerr<<"No more iterations ! "<<endl;
            } // end of attempt to separate overlapped clusters
        } // end of nbins cut 
        if(elps < 0.) cout<<" elps=-1 ="<<elps<<endl;
        if(elps >0. && elps< 0.3 && nfhits == 1) cout<<" small elps, nfh=1 ="
                                                     <<elps<<","<<nfhits<<endl;
        if(nfhits == 2) cout<<" nfhits=2 ="<<nfhits<<endl;
        for (i=0; i<nfhits; i++) {
            xfit[i] *= (1.e+4/scl);
            if(wing == 1) xfit[i] *= (-1);
            zfit[i] *= 1.e+4;
            //       cout<<" ---------  i,xfiti,zfiti,qfiti ="<<i<<","
            // <<xfit[i]<<","<<zfit[i]<<","<<qfit[i]<<endl;
        } // end for i
        Int_t ncl = nfhits;
        if(nfhits == 1 && separate == 1) {
            cout<<"!!!!! no separate"<<endl;
            ncl = -2;
        }  // end if
        if(nfhits == 2) {
            cout << "Split cluster: " << endl;
            clusterJ->PrintInfo();
            cout << " in: " << endl;
            for (i=0; i<nfhits; i++) {
                // AliITSRawClusterSDD *clust = new AliITSRawClusterSDD(wing,
                                               -1,-1,(Double_t)qfit[i],ncl,0,0,
                                               (Double_t)xfit[i],
                                               (Double_t)zfit[i],0,0,0,0,
                                                tstart,tstop,astart,astop);
            //    AliITSRawClusterSDD *clust = new AliITSRawClusterSDD(wing,-1,
            //                                 -1,(Double_t)qfit[i],0,0,0,
            //                                  (Double_t)xfit[i],
            //                                  (Double_t)zfit[i],0,0,0,0,
            //                                  tstart,tstop,astart,astop,ncl);
            // ???????????
            // if(wing == 1) xfit[i] *= (-1);
            Double_t Anode = (zfit[i]/anodePitch+fNofAnodes/2-0.5);
            Double_t Time = (fSddLength - xfit[i])/fDriftSpeed;
            Double_t clusterPeakAmplitude = clusterJ->PeakAmpl();
            Double_t peakpos = clusterJ->PeakPos();
            Double_t clusteranodePath = (Anode - fNofAnodes/2)*anodePitch;
            Double_t clusterDriftPath = Time*fDriftSpeed;
            clusterDriftPath = fSddLength-clusterDriftPath;
            AliITSRawClusterSDD *clust = new AliITSRawClusterSDD(wing,Anode,
                                                                 Time,qfit[i],
                                               clusterPeakAmplitude,peakpos,
                                               0.,0.,clusterDriftPath,
                                         clusteranodePath,clusterJ->Samples()/2
                                    ,tstart,tstop,0,0,0,astart,astop);
            clust->PrintInfo();
            iTS->AddCluster(1,clust);
            //    cout<<"new cluster added: tstart,tstop,astart,astop,x,ncl ="
            // <<tstart<<","<<tstop<<","<<astart<<","<<astop<<","<<xfit[i]
            // <<","<<ncl<<endl;
            delete clust;
        }// nfhits loop
        Clusters()->RemoveAt(j);
    } // if nfhits = 2
} // cluster loop
Clusters()->Compress();
Map()->ClearMap(); 
*/
    return;
}
//______________________________________________________________________
void AliITSClusterFinderSDD::GetRecPoints(){
    // get rec points
    static AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
    // get number of clusters for this module
    Int_t nofClusters = NClusters();
    nofClusters -= fNclusters;
    const Double_t kconvGeV = 1.e-6; // GeV -> KeV
    const Double_t kconv = 1.0e-4; 
    const Double_t kRMSx = 38.0*kconv; // microns->cm ITS TDR Table 1.3
    const Double_t kRMSz = 28.0*kconv; // microns->cm ITS TDR Table 1.3
    Int_t i;
    Int_t ix, iz, idx=-1;
    AliITSdigitSDD *dig=0;
    Int_t ndigits=NDigits();
    for(i=0; i<nofClusters; i++) { 
        AliITSRawClusterSDD *clusterI = (AliITSRawClusterSDD*)Cluster(i);
        if(!clusterI) Error("SDD: GetRecPoints","i clusterI ",i,clusterI);
        if(clusterI) idx=clusterI->PeakPos();
        if(idx>ndigits) Error("SDD: GetRecPoints","idx ndigits",idx,ndigits);
        // try peak neighbours - to be done 
        if(idx&&idx<= ndigits) dig =(AliITSdigitSDD*)GetDigit(idx);
        if(!dig) {
            // try cog
            GetSeg()->GetPadIxz(clusterI->X(),clusterI->Z(),ix,iz);
            dig = (AliITSdigitSDD*)Map()->GetHit(iz-1,ix-1);
            // if null try neighbours
            if (!dig) dig = (AliITSdigitSDD*)Map()->GetHit(iz-1,ix); 
            if (!dig) dig = (AliITSdigitSDD*)Map()->GetHit(iz-1,ix+1); 
            if (!dig) printf("SDD: cannot assign the track number!\n");
        } //  end if !dig
        AliITSRecPoint rnew;
        rnew.SetX(clusterI->X());
        rnew.SetZ(clusterI->Z());
        rnew.SetQ(clusterI->Q());   // in KeV - should be ADC
        rnew.SetdEdX(kconvGeV*clusterI->Q());
        rnew.SetSigmaX2(kRMSx*kRMSx);
        rnew.SetSigmaZ2(kRMSz*kRMSz);

        if(dig) rnew.fTracks[0]=dig->GetTrack(0);
        if(dig) rnew.fTracks[1]=dig->GetTrack(1);
        if(dig) rnew.fTracks[2]=dig->GetTrack(2);

        iTS->AddRecPoint(rnew);
    } // I clusters
//    Map()->ClearMap();
}
//______________________________________________________________________
void AliITSClusterFinderSDD::FindRawClusters(Int_t mod){
    // find raw clusters
    
    SetModule(mod);
    Find1DClustersE();
    GroupClusters();
    SelectClusters();
    ResolveClustersE();
    GetRecPoints();
}
//_______________________________________________________________________
void AliITSClusterFinderSDD::Print() const{
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
