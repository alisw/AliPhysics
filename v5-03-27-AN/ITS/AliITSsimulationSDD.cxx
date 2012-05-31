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

/* $Id$ */

#include <Riostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TRandom.h>
#include <TROOT.h>
#include "AliITS.h"
#include "AliITSMapA2.h"
#include "AliITSdigitSPD.h"
#include "AliITSetfSDD.h"
#include "AliITSmodule.h"
#include "AliITShit.h"
#include "AliITSpList.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliITSsimulationSDD.h"
#include "AliLog.h"
#include "AliRun.h"

ClassImp(AliITSsimulationSDD)
////////////////////////////////////////////////////////////////////////
// Version: 0                                                         //
// Written by Piergiorgio Cerello                                     //
// November 23 1999                                                   //
//                                                                    //
// AliITSsimulationSDD is the simulation of SDDs.                     //
////////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSsimulationSDD::AliITSsimulationSDD():
AliITSsimulation(),
fITS(0),
fHitMap2(0),
fHitSigMap2(0),
fHitNoiMap2(0),
fElectronics(0),
fInZR(0),
fInZI(0),
fOutZR(0),
fOutZI(0),
fAnodeFire(0),
fHis(0),
fFlag(kFALSE),
fCrosstalkFlag(kFALSE),
fDoFFT(1),
fNofMaps(0),
fMaxNofSamples(0),
fScaleSize(0){
    // Default constructor
    SetPerpendTracksFlag();
    SetCrosstalkFlag();
    SetDoFFT();
}

//______________________________________________________________________
AliITSsimulationSDD::AliITSsimulationSDD(AliITSDetTypeSim* dettyp):
AliITSsimulation(dettyp),
fITS(0),
fHitMap2(0),
fHitSigMap2(0),
fHitNoiMap2(0),
fElectronics(0),
fInZR(0),
fInZI(0),
fOutZR(0),
fOutZI(0),
fAnodeFire(0),
fHis(0),
fFlag(kFALSE),
fCrosstalkFlag(kFALSE),
fDoFFT(1),
fNofMaps(0),
fMaxNofSamples(0),
fScaleSize(0){
    // Default Constructor
  Init();
}
//______________________________________________________________________
void AliITSsimulationSDD::Init(){
    // Standard Constructor

    AliITSsegmentationSDD* seg = (AliITSsegmentationSDD*)GetSegmentationModel(1);
    fScaleSize = ScaleFourier(seg);
    SetPerpendTracksFlag();
    SetCrosstalkFlag();
    SetDoFFT();

    AliITSSimuParam* simpar = fDetType->GetSimuParam();
    fpList = new AliITSpList( seg->Npz(),
                              fScaleSize*seg->Npx() );
    fHitSigMap2 = new AliITSMapA2(seg,fScaleSize,1);
    fHitNoiMap2 = new AliITSMapA2(seg,fScaleSize,1);
    fHitMap2 = fHitSigMap2;

    fNofMaps = seg->Npz();
    fMaxNofSamples = seg->Npx();
    fAnodeFire = new Bool_t [fNofMaps];
    
    Float_t sddWidth  = seg->Dz();
    Float_t anodePitch = seg->Dpz(0);
    Double_t timeStep  = (Double_t)seg->Dpx(0);

    if(anodePitch*(fNofMaps/2) > sddWidth) {
      AliWarning(Form("Too many anodes %d or too big pitch %f ",
                fNofMaps/2,anodePitch));
    } // end if


    fElectronics = new AliITSetfSDD(timeStep/fScaleSize,
                                    simpar->GetSDDElectronics());


    fITS       = (AliITS*)gAlice->GetModule("ITS");
 
    fInZR  = new Double_t [fScaleSize*fMaxNofSamples];
    fInZI  = new Double_t [fScaleSize*fMaxNofSamples];
    fOutZR = new Double_t [fScaleSize*fMaxNofSamples];
    fOutZI = new Double_t [fScaleSize*fMaxNofSamples];  
}
//______________________________________________________________________
AliITSsimulationSDD::~AliITSsimulationSDD() { 
    // destructor

    //    delete fpList;
    delete fHitSigMap2;
    delete fHitNoiMap2;
    delete fElectronics;

    fITS = 0;

    if (fHis) {
        fHis->Delete(); 
        delete fHis;     
    } // end if fHis
    if(fInZR)  delete [] fInZR;
    if(fInZI)  delete [] fInZI;        
    if(fOutZR) delete [] fOutZR;
    if(fOutZI) delete [] fOutZI;
    if(fAnodeFire) delete [] fAnodeFire;
}
//______________________________________________________________________
void AliITSsimulationSDD::InitSimulationModule( Int_t module, Int_t event ) {
    // create maps to build the lists of tracks for each summable digit
    fModule = module;
    fEvent  = event;
    ClearMaps();
    memset(fAnodeFire,0,sizeof(Bool_t)*fNofMaps);    
}
//______________________________________________________________________
void AliITSsimulationSDD::ClearMaps() {
    // clear maps
    fpList->ClearMap();
    fHitSigMap2->ClearMap();
    fHitNoiMap2->ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSDD::FastFourierTransform(Double_t *real,
                          Double_t *imag,Int_t direction) {
    // Do a Fast Fourier Transform

    Int_t samples = fElectronics->GetSamples();
    Int_t l = (Int_t) ((log((Float_t) samples)/log(2.))+0.5);
    Int_t m1 = samples;
    Int_t m  = samples/2;
    Int_t m2 = samples/m1;
    Int_t i,j,k;
    for(i=1; i<=l; i++) {
        for(j=0; j<samples; j += m1) {
            Int_t p = 0;
            for(k=j; k<= j+m-1; k++) {
                Double_t wsr = fElectronics->GetWeightReal(p);
                Double_t wsi = fElectronics->GetWeightImag(p);
                if(direction == -1) wsi = -wsi;
                Double_t xr = *(real+k+m);
                Double_t xi = *(imag+k+m);
                *(real+k+m) = wsr*(*(real+k)-xr) - wsi*(*(imag+k)-xi);
                *(imag+k+m) = wsr*(*(imag+k)-xi) + wsi*(*(real+k)-xr);
                *(real+k) += xr;
                *(imag+k) += xi;
                p += m2;
            } // end for k
        } // end for j
        m1 = m;
        m /= 2;
        m2 += m2;
    } // end for i
    for(j=0; j<samples; j++) {
        Int_t j1 = j;
        Int_t p = 0;
        Int_t i1;
        for(i1=1; i1<=l; i1++) {
            Int_t j2 = j1;
            j1 /= 2;
            p = p + p + j2 - j1 - j1;
        } // end for i1
        if(p >= j) {
            Double_t xr = *(real+j);
            Double_t xi = *(imag+j);
            *(real+j) = *(real+p);
            *(imag+j) = *(imag+p);
            *(real+p) = xr;
            *(imag+p) = xi;
        } // end if p>=j
    } // end for j
    if(direction == -1) {
        for(i=0; i<samples; i++) {
            *(real+i) /= samples;
            *(imag+i) /= samples;
        } // end for i
    } // end if direction == -1
    return;
}

//______________________________________________________________________
void AliITSsimulationSDD::SDigitiseModule(AliITSmodule *mod,Int_t md,Int_t ev){
    // digitize module using the "slow" detector simulator creating
    // summable digits.

    TObjArray *fHits = mod->GetHits();
    Int_t nhits      = fHits->GetEntriesFast();
    if( !nhits ) return;

    InitSimulationModule( md, ev );
    HitsToAnalogDigits( mod );  // fills fHitMap2 which is = fHitSigmap2
    ChargeToSignal( fModule,kFALSE,kTRUE ); // - Process signal adding gain without adding noise
    fHitMap2 = fHitNoiMap2;   // - Swap to noise map
    ChargeToSignal( fModule,kTRUE,kFALSE );  // - Process only noise
    fHitMap2 = fHitSigMap2;   // - Return to signal map
    WriteSDigits();
    ClearMaps();
}
//______________________________________________________________________
Bool_t AliITSsimulationSDD::AddSDigitsToModule(TClonesArray *pItemArray,
                                               Int_t mask ) {
    // Add Summable digits to module maps.
    AliITSSimuParam* simpar = fDetType->GetSimuParam();
    Int_t    nItems = pItemArray->GetEntries();
    Double_t maxadc = simpar->GetSDDMaxAdc();
    Bool_t sig = kFALSE;
    
    // cout << "Adding "<< nItems <<" SDigits to module " << fModule << endl;
    for( Int_t i=0; i<nItems; i++ ) {
        AliITSpListItem * pItem = (AliITSpListItem *)(pItemArray->At( i ));
        if( pItem->GetModule() != fModule ) {
            Error( "AliITSsimulationSDD","Error reading, SDigits module "
                   "%d != current module %d: exit",
                   pItem->GetModule(), fModule );
            return sig;
        } // end if

        if(pItem->GetSignal()>0.0 ) sig = kTRUE;
        
        fpList->AddItemTo( mask, pItem ); // Add SignalAfterElect + noise
        AliITSpListItem * pItem2 = fpList->GetpListItem( pItem->GetIndex() );
        Double_t sigAE = pItem2->GetSignalAfterElect();
        if( sigAE >= maxadc ) sigAE = maxadc-1; // avoid overflow signal
        Int_t ia;
        Int_t it;
        fpList->GetMapIndex( pItem->GetIndex(), ia, it );
        fHitMap2->SetHit( ia, it, sigAE );
        fAnodeFire[ia] = kTRUE;
    }
    return sig;
}
//______________________________________________________________________
void AliITSsimulationSDD::FinishSDigitiseModule() {
    // digitize module using the "slow" detector simulator from
    // the sum of summable digits.
    FinishDigits() ;
    ClearMaps();
}
//______________________________________________________________________
void AliITSsimulationSDD::DigitiseModule(AliITSmodule *mod,Int_t md,Int_t ev){
    // create maps to build the lists of tracks for each digit

    TObjArray *fHits = mod->GetHits();
    Int_t nhits      = fHits->GetEntriesFast();

    InitSimulationModule( md, ev );
    if( !nhits ) return;
        
    HitsToAnalogDigits( mod );
    ChargeToSignal( fModule,kTRUE,kTRUE );  // process signal + noise

    for( Int_t i=0; i<fNofMaps; i++ ) {
        for( Int_t j=0; j<fMaxNofSamples; j++ ) {
            Int_t jdx = j*fScaleSize;
            Int_t index = fpList->GetHitIndex( i, j );
            AliITSpListItem pItemTmp2( fModule, index, 0. );
            // put the fScaleSize analog digits in only one
            for( Int_t ik=0; ik<fScaleSize; ik++ ) {
                AliITSpListItem *pItemTmp = fpList->GetpListItem( i, jdx+ik );
                if( pItemTmp == 0 ) continue;
                pItemTmp2.Add( pItemTmp );
            }
            fpList->DeleteHit( i, j );
            fpList->AddItemTo( 0, &pItemTmp2 );
        }
    }
    FinishDigits();
    ClearMaps();
}
//______________________________________________________________________
void AliITSsimulationSDD::FinishDigits() {
    // introduce the electronics effects and do zero-suppression if required

    if( fCrosstalkFlag ) ApplyCrosstalk(fModule);

    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);
    Bool_t isZeroSupp = res->GetZeroSupp();
    if (isZeroSupp) Compress2D();
    else StoreAllDigits();
}
//______________________________________________________________________
void AliITSsimulationSDD::HitsToAnalogDigits( AliITSmodule *mod ) {
    // create maps to build the lists of tracks for each digit
  AliITSsegmentationSDD* seg = (AliITSsegmentationSDD*)GetSegmentationModel(1);
  AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);
  AliITSSimuParam* simpar = fDetType->GetSimuParam();
  TObjArray *hits     = mod->GetHits();
  Int_t      nhits    = hits->GetEntriesFast();

  //    Int_t      arg[6]   = {0,0,0,0,0,0};
  Int_t     nofAnodes  = fNofMaps/2;
  Double_t  sddLength  = seg->Dx();
  Double_t  anodePitch = seg->Dpz(0);
  Double_t  timeStep   = seg->Dpx(0);
  Double_t  driftSpeed ;  // drift velocity (anode dependent)
  Double_t  nanoampToADC       = simpar->GetSDDMaxAdc()/simpar->GetSDDDynamicRange(); //   maxadc/topValue;
  Double_t  cHloss     = simpar->GetSDDChargeLoss();
  Float_t   dfCoeff, s1; 
  simpar->GetSDDDiffCoeff(dfCoeff,s1); // Signal 2d Shape
  Double_t  eVpairs    = simpar->GetGeVToCharge()*1.0E9; // 3.6 eV by def.
  Double_t  nsigma     = simpar->GetNSigmaIntegration(); //
  Int_t     nlookups   = simpar->GetGausNLookUp();       //
  Float_t   jitter     = simpar->GetSDDJitterError(); // 
  Float_t   mapsmear   = simpar->GetSDDCorrMapPrecision(); // 
  Float_t   trigDelay  = simpar->GetSDDTrigDelay(); // compensation for MC time zero
  if(res->IsAMAt20MHz()) trigDelay+=12.5; // compensation for discretization step

  Float_t   timeZero=fDetType->GetResponseSDD()->GetTimeZero(fModule);
  Float_t   adcscale   = fDetType->GetResponseSDD()->GetADCtokeV(fModule);
  adcscale/=simpar->GetSDDkeVtoADC();

  // Piergiorgio's part (apart for few variables which I made float
  // when i thought that can be done
  // Fill detector maps with GEANT hits
  // loop over hits in the module
  
  const Float_t kconv = 1.0e+6;  // GeV->KeV
  Int_t     itrack      = 0;
  Int_t     iWing;       // which detector wing/side.
  Int_t     ii,kk,ka,kt; // loop indexs
  Int_t     ia,it,index; // sub-pixel integration indexies
  Int_t     iAnode;      // anode number.
  Int_t     timeSample;  // time buckett.
  Int_t     anodeWindow; // anode direction charge integration width
  Int_t     timeWindow;  // time direction charge integration width
  Int_t     jamin,jamax; // anode charge integration window
  Int_t     jtmin,jtmax; // time charge integration window
  Int_t     nsplitAn;    // the number of splits in anode and time windows
  Int_t     nsplitTb;    // the number of splits in anode and time windows
  Int_t     nOfSplits;   // number of times track length is split into
  Float_t   nOfSplitsF;  // Floating point version of nOfSplits.
  Float_t   kkF;         // Floating point version of loop index kk.
  Double_t  pathInSDD; // Track length in SDD.
  Double_t  drPath; // average position of track in detector. in microns
  Double_t  drTime; // Drift time
  Double_t  avDrft;  // x position of path length segment in cm.
  Double_t  avAnode; // Anode for path length segment in Anode number (float)
  Double_t  zAnode;  // Floating point anode number.
  Double_t  driftPath; // avDrft in microns.
  Double_t  width;     // width of signal at anodes.
  Double_t  depEnergy; // Energy deposited in this GEANT step.
  Double_t  xL[3],dxL[3]; // local hit coordinates and diff.
  Double_t  sigA; // sigma of signal at anode.
  Double_t  sigT; // sigma in time/drift direction for track segment
  Double_t  aStep,aConst; // sub-pixel size and offset anode
  Double_t  tStep,tConst; // sub-pixel size and offset time
  Double_t  amplitude; // signal amplitude for track segment in nanoAmpere
  Double_t  chargeloss; // charge loss for track segment.
  Double_t  anodeAmplitude; // signal amplitude in anode direction
  Double_t  aExpo;          // exponent of Gaussian anode direction
  Double_t  timeAmplitude;  // signal amplitude in time direction
  Double_t  tExpo;          // exponent of Gaussian time direction
  Double_t  tof;            // Time of flight in ns of this step.    
  
  for(ii=0; ii<nhits; ii++) {
    if(!mod->LineSegmentL(ii,xL[0],dxL[0],xL[1],dxL[1],xL[2],dxL[2],
			  depEnergy,itrack)) continue;
    Float_t xloc=xL[0];
    if(xloc>0) iWing=0; // left side, carlos channel 0
    else iWing=1; // right side
    
    Float_t zloc=xL[2]+0.5*dxL[2];
    zAnode=seg->GetAnodeFromLocal(xloc,zloc); // anode number in the range 0.-511.
    driftSpeed = res->GetDriftSpeedAtAnode(zAnode);
    driftSpeed+= fDetType->GetResponseSDD()->GetDeltaVDrift(fModule,zAnode>255);

    if(timeStep*fMaxNofSamples < sddLength/driftSpeed) {
      AliWarning("Time Interval > Allowed Time Interval");
    }
    depEnergy  *= kconv;
    if (!depEnergy) {
      AliDebug(1,
	       Form("fTrack = %d hit=%d module=%d This particle has passed without losing energy!",
		    itrack,ii,mod->GetIndex()));
      continue;
      // continue if the particle did not lose energy
      // passing through detector
    } // end if !depEnergy
     
    tof=0.;
    AliITShit* h=(AliITShit*)hits->At(ii);
    if(h){ 
      tof=h->GetTOF()*1E9; 
      AliDebug(1,Form("TOF for hit %d on mod %d (particle %d)=%g",ii,fModule,h->Track(),tof));
    }

    Float_t corrx=0, corrz=0;
    res->GetShiftsForSimulation(xL[2],xL[0],corrz,corrx,seg);
    xL[2]-=corrz;
    xL[0]-=corrx;
    xL[0] += 0.0001*gRandom->Gaus( 0, mapsmear); //
    xL[0] += 0.0001*gRandom->Gaus( 0, jitter ); //

    pathInSDD = TMath::Sqrt(dxL[0]*dxL[0]+dxL[1]*dxL[1]+dxL[2]*dxL[2]);
    
    if (fFlag && pathInSDD) { depEnergy *= (0.03/pathInSDD); }
    drPath = TMath::Abs(10000.*(dxL[0]+2.*xL[0])*0.5);
    drPath = sddLength-drPath;
    if(drPath < 0) {
      AliInfo( // this should be fixed at geometry level
	       Form("negative drift path drPath=%e sddLength=%e dxL[0]=%e xL[0]=%e",
		    drPath,sddLength,dxL[0],xL[0]));
      continue;
    } // end if drPath < 0
    
    // Compute number of segments to brake step path into
    drTime = drPath/driftSpeed;  //   Drift Time
    sigA   = TMath::Sqrt(2.*dfCoeff*drTime+s1*s1);// Sigma along the anodes
    // calcuate the number of time the path length should be split into.
    nOfSplits = (Int_t) (1. + 10000.*pathInSDD/sigA);
    if(fFlag) nOfSplits = 1;
    
    // loop over path segments, init. some variables.
    depEnergy /= nOfSplits;
    nOfSplitsF = (Float_t) nOfSplits;
    Float_t theAverage=0.,theSteps=0.;
    for(kk=0;kk<nOfSplits;kk++) { // loop over path segments
      kkF       = (Float_t) kk + 0.5;
      avDrft    = xL[0]+dxL[0]*kkF/nOfSplitsF;
      avAnode   = xL[2]+dxL[2]*kkF/nOfSplitsF;
      theSteps+=1.;
      theAverage+=avAnode;
      zAnode = seg->GetAnodeFromLocal(avDrft,avAnode);
      driftSpeed = res->GetDriftSpeedAtAnode(zAnode);	
      driftSpeed+= fDetType->GetResponseSDD()->GetDeltaVDrift(fModule,zAnode>255);
      driftPath = TMath::Abs(10000.*avDrft);
      driftPath = sddLength-driftPath;
      if(driftPath < 0) {
	AliDebug(1, // this should be fixed at geometry level
		 Form("negative drift path driftPath=%e sddLength=%e avDrft=%e dxL[0]=%e xL[0]=%e",
		      driftPath,sddLength,avDrft,dxL[0],xL[0]));
	continue;
      } // end if driftPath < 0
      drTime     = driftPath/driftSpeed; // drift time for segment.
      // Sigma along the anodes for track segment.
      sigA       = TMath::Sqrt(2.*dfCoeff*drTime+s1*s1);
      sigT       = sigA/driftSpeed;

      drTime+=tof; // take into account Time Of Flight from production point
      drTime-=trigDelay;
      drTime+=timeZero;
      timeSample = (Int_t) (fScaleSize*drTime/timeStep + 1.001); // time bin in range 1-256 !!!
      if(zAnode>nofAnodes) zAnode-=nofAnodes;  // to have the anode number between 0. and 256.
      iAnode = (Int_t) (1.001+zAnode); // iAnode in range 1-256 !!!!
      
	// Peak amplitude in nanoAmpere
      amplitude  = fScaleSize*160.*depEnergy/
	(timeStep*eVpairs*2.*acos(-1.));
      chargeloss = 1.-cHloss*driftPath/1000.;
      amplitude *= chargeloss;
      amplitude *= adcscale;
      width  = 2.*nsigma/(nlookups-1);
      // Spread the charge 
      nsplitAn = 4; 
      nsplitTb=4;
      aStep  = anodePitch/(nsplitAn*sigA);
      aConst = zAnode*anodePitch/sigA;
      tStep  = timeStep/(nsplitTb*fScaleSize*sigT);
      tConst = drTime/sigT;
      // Define SDD window corresponding to the hit
      anodeWindow = (Int_t)(nsigma*sigA/anodePitch+1);
      timeWindow  = (Int_t) (fScaleSize*nsigma*sigT/timeStep+1.);
      jamin = (iAnode - anodeWindow - 2)*nsplitAn+1;
      if(jamin <= 0) jamin = 1;
      if(jamin > nofAnodes*nsplitAn){ 
	AliDebug(1,Form("Energy deposition completely outside anode acceptance: anode min=%d",jamin));
	continue;
      }
      jamax = (iAnode + anodeWindow + 2)*nsplitAn;
      if(jamax > nofAnodes*nsplitAn) jamax = nofAnodes*nsplitAn;
      if(jamax <=0){ 
	AliDebug(1,Form("Energy deposition completely outside anode acceptance: anode max=%d",jamax));
	continue;
      }
      jtmin = (Int_t)(timeSample-timeWindow-2)*nsplitTb+1;
      if(jtmin <= 0) jtmin = 1;
      if(jtmin > fScaleSize*fMaxNofSamples*nsplitTb){ 
	AliDebug(1,Form("Energy deposition completely outside time acceptance: time sample min=%d  tof=%f",jtmin,tof));
	continue; 
      }
      jtmax = (Int_t)(timeSample+timeWindow+2)*nsplitTb;
      if(jtmax > fScaleSize*fMaxNofSamples*nsplitTb) jtmax = fScaleSize*fMaxNofSamples*nsplitTb;
      if(jtmax <= 0){
	AliDebug(1,Form("Energy deposition completely outside time acceptance: time sample max=%d  tof=%f",jtmax,tof));
	continue; 
      }

      // Spread the charge in the anode-time window
      for(ka=jamin; ka <=jamax; ka++) {	  
	ia = (ka-1)/nsplitAn + 1;
	if(ia <= 0) ia=1; 
	if(ia > nofAnodes) ia = nofAnodes;
	aExpo     = (aStep*(ka-0.5)-aConst);
	if(TMath::Abs(aExpo) > nsigma)  anodeAmplitude = 0.;
	else {
	  Int_t theBin = (Int_t) ((aExpo+nsigma)/width+0.5);
	  anodeAmplitude = amplitude*simpar->GetGausLookUp(theBin);
	}
	// index starts from 0
	index = iWing*nofAnodes+ia-1;
	if(anodeAmplitude){
	  for(kt=jtmin; kt<=jtmax; kt++) {
	    it = (kt-1)/nsplitTb+1;  // it starts from 1
	    if(it<=0) it=1;
	    if(it>fScaleSize*fMaxNofSamples)
	      it = fScaleSize*fMaxNofSamples;
	    tExpo    = (tStep*(kt-0.5)-tConst);
	    if(TMath::Abs(tExpo) > nsigma) timeAmplitude = 0.;
	    else {
	      Int_t theBin = (Int_t) ((tExpo+nsigma)/width+0.5);
	      timeAmplitude = anodeAmplitude*simpar->GetGausLookUp(theBin)*aStep*tStep;
	    }
	    timeAmplitude *= nanoampToADC;
	    //         ListOfFiredCells(arg,timeAmplitude,alst,padr);
	    Double_t charge = timeAmplitude;
	    charge += fHitMap2->GetSignal(index,it-1);
	    fHitMap2->SetHit(index, it-1, charge);
	    fpList->AddSignal(index,it-1,itrack,ii-1,
			      mod->GetIndex(),timeAmplitude);
	    fAnodeFire[index] = kTRUE;
	  }  // end loop over time in window               
	} // end if anodeAmplitude 
      } // loop over anodes in window
    } // end loop over "sub-hits"
  } // end loop over hits
}

//____________________________________________
void AliITSsimulationSDD::AddDigit( Int_t i, Int_t j, Int_t signalc, Int_t signale) {
  // Adds a Digit.
  Int_t size = AliITSdigit::GetNTracks();

  Int_t digits[3];
  Int_t * tracks = new Int_t[size];
  Int_t * hits = new Int_t[size];
  Float_t phys;
  Float_t * charges = new Float_t[size];

  digits[0] = i;
  digits[1] = j;
  digits[2] = signalc;

  AliITSpListItem *pItem = fpList->GetpListItem( i, j );
  if( pItem == 0 ) {
    phys = 0.0;
    for( Int_t l=0; l<size; l++ ) {
      tracks[l]  = 0;
      hits[l]    = 0;
      charges[l] = 0.0;
    }
  } else {
    Int_t idtrack =  pItem->GetTrack( 0 );
    if( idtrack >= 0 ) phys = pItem->GetSignal();  
    else phys = 0.0;

    for( Int_t l=0; l<size; l++ ) if(l<pItem->GetMaxKept()) {
      tracks[l]  = pItem->GetTrack( l );
      hits[l]    = pItem->GetHit( l );
      charges[l] = pItem->GetSignal( l );
    }else{
      tracks[l]  = -3;
      hits[l]    = -1;
      charges[l] = 0.0;
    }// end for if
  }

  fITS->AddSimDigit( 1, phys, digits, tracks, hits, charges, signale ); 
  delete [] tracks;
  delete [] hits;
  delete [] charges;
}
//______________________________________________________________________
void AliITSsimulationSDD::ChargeToSignal(Int_t mod,Bool_t bAddNoise, Bool_t bAddGain) {
  // add baseline, noise, gain, electronics and ADC saturation effects
  // apply dead channels

  AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(mod);
  Double_t baseline=0; 
  Double_t noise=0; 
  Double_t gain=0; 
  Float_t contrib=0;
  Int_t i,k,kk;
  AliITSSimuParam* simpar = fDetType->GetSimuParam();
  Float_t maxadc = simpar->GetSDDMaxAdc();    
  Int_t nGroup=fScaleSize;
  if(res->IsAMAt20MHz()){
    nGroup=fScaleSize/2;
  }

  for (i=0;i<fNofMaps;i++) {
    if( !fAnodeFire[i] ) continue;
    baseline = res->GetBaseline(i);
    noise = res->GetNoise(i);
    gain = res->GetChannelGain(i)/fDetType->GetAverageGainSDD();
    if(res->IsBad()) gain=0.;
    if( res->IsChipBad(res->GetChip(i)) )gain=0.;
    for(k=0; k<fScaleSize*fMaxNofSamples; k++) {
      fInZR[k]  = fHitMap2->GetSignal(i,k);
      if(bAddGain) fInZR[k]*=gain;
      if( bAddNoise ) {
	contrib   = (baseline + noise*gRandom->Gaus());
	fInZR[k] += contrib;
      }
      fInZI[k]  = 0.;
    } // end for k
    if(!fDoFFT) {      
      for(k=0; k<fMaxNofSamples; k++) {
	Double_t newcont = 0.;
	Double_t maxcont = 0.;
	for(kk=0;kk<fScaleSize;kk++) {
	  newcont = fInZR[fScaleSize*k+kk];
	  if(newcont > maxcont) maxcont = newcont;
	} // end for kk
	newcont = maxcont;
	if (newcont >= maxadc) newcont = maxadc -1;
	if(newcont >= baseline){
	  Warning("","newcont=%f>=baseline=%f",newcont,baseline);
	} // end if
	  // back to analog: ?
	fHitMap2->SetHit(i,k,newcont);
      }  // end for k
    }else{
      FastFourierTransform(&fInZR[0],&fInZI[0],1);
      for(k=0; k<fScaleSize*fMaxNofSamples; k++) {
	Double_t rw = fElectronics->GetTraFunReal(k);
	Double_t iw = fElectronics->GetTraFunImag(k);
	fOutZR[k]   = fInZR[k]*rw - fInZI[k]*iw;
	fOutZI[k]   = fInZR[k]*iw + fInZI[k]*rw;
      } // end for k
      FastFourierTransform(&fOutZR[0],&fOutZI[0],-1);
      for(k=0; k<fMaxNofSamples; k++) {
	Double_t newcont1 = 0.;
	Double_t maxcont1 = 0.;
	for(kk=0;kk<nGroup;kk++) {
	  newcont1 = fOutZR[fScaleSize*k+kk];
	  if(newcont1 > maxcont1) maxcont1 = newcont1;
	} // end for kk
	newcont1 = maxcont1;
	if (newcont1 >= maxadc) newcont1 = maxadc -1;
	fHitMap2->SetHit(i,k,newcont1);
      } // end for k
    }
  } // end for i loop over anodes
  return;
}

//______________________________________________________________________
void AliITSsimulationSDD::ApplyCrosstalk(Int_t mod) {
    // function add the crosstalk effect to signal
    // temporal function, should be checked...!!!
  
    // create and inizialice crosstalk map
    Float_t* ctk = new Float_t[fNofMaps*fMaxNofSamples+1];
    memset( ctk, 0, sizeof(Float_t)*(fNofMaps*fMaxNofSamples+1) );
    AliITSCalibrationSDD* calibr = (AliITSCalibrationSDD*)GetCalibrationModel(mod);
    for( Int_t z=0; z<fNofMaps; z++ ) {
      Double_t baseline = calibr->GetBaseline(z);
        Bool_t on = kFALSE;
        Int_t tstart = 0;
        Int_t tstop = 0;
        Int_t nTsteps = 0;
        
        for( Int_t l=0; l<fMaxNofSamples; l++ ) {
            Float_t fadc = (Float_t)fHitMap2->GetSignal( z, l );
            if( fadc > baseline ) {
                if( on == kFALSE && l<fMaxNofSamples-4 ) {
                    Float_t fadc1 = (Float_t)fHitMap2->GetSignal( z, l+1 );
                    if( fadc1 < fadc ) continue;
                    on = kTRUE;
                    nTsteps = 0;
                    tstart = l;
                }
                nTsteps++;
            }
            else { // end fadc > baseline
                if( on == kTRUE ) {        
                    if( nTsteps > 2 ) {
                        tstop = l;
                        // make smooth derivative
                        Float_t* dev = new Float_t[fMaxNofSamples+1];
                        memset( dev, 0, sizeof(Float_t)*(fMaxNofSamples+1) );
                        for( Int_t i=tstart; i<tstop; i++ ) {   
                            if( i > 2 && i < fMaxNofSamples-2 )
                                dev[i] = -0.2*fHitMap2->GetSignal( z,i-2 ) 
                                    -0.1*fHitMap2->GetSignal( z,i-1 ) 
                                    +0.1*fHitMap2->GetSignal( z,i+1 ) 
                                    +0.2*fHitMap2->GetSignal( z,i+2 );
                        }
                        
                        // add crosstalk contribution to neibourg anodes  
                        for( Int_t i=tstart; i<tstop; i++ ) {
                            Int_t anode = z - 1;
                            Int_t i1 = (Int_t)((i-tstart)*.61+tstart+0.5); // 
                            Float_t ctktmp =  -dev[i1] * 0.25;
                            if( anode > 0 ) {
                                ctk[anode*fMaxNofSamples+i] += ctktmp;
                            }
                            anode = z + 1;
                            if( anode < fNofMaps ) {
                                ctk[anode*fMaxNofSamples+i] += ctktmp;
                            }
                        }
                        delete [] dev;
                        
                    } // if( nTsteps > 2 )
                    on = kFALSE;
                }  // if( on == kTRUE )
            }  // else
        }
    }
    
    for( Int_t a=0; a<fNofMaps; a++ )
        for( Int_t t=0; t<fMaxNofSamples; t++ ) {     
            Float_t signal = fHitMap2->GetSignal(a,t)+ctk[a*fMaxNofSamples+t];
            fHitMap2->SetHit( a, t, signal );
        }

    delete [] ctk;
}

//______________________________________________________________________
Int_t AliITSsimulationSDD::Convert10to8(Int_t signal) const {
    // To the 10 to 8 bit lossive compression.
    // code from Davide C. and Albert W.

    if (signal < 128)  return signal;
    if (signal < 256)  return (128+((signal-128)>>1));
    if (signal < 512)  return (192+((signal-256)>>3));
    if (signal < 1024) return (224+((signal-512)>>4));
    return 0;
}
//______________________________________________________________________
Int_t AliITSsimulationSDD::Convert8to10(Int_t signal) const {
  // Decompression from 8 to 10 bit

  if (signal < 0 || signal > 255) {
    AliWarning(Form("Signal value %d out of range",signal));
    return 0;
  } // end if signal <0 || signal >255

  if (signal < 128) return signal;
  if (signal < 192) {
    if (TMath::Odd(signal)) return (128+((signal-128)<<1));
    else  return (128+((signal-128)<<1)+1);
  } // end if signal < 192
  if (signal < 224) {
    if (TMath::Odd(signal)) return (256+((signal-192)<<3)+3);
    else  return (256+((signal-192)<<3)+4);
  } // end if signal < 224
  if (TMath::Odd(signal)) return (512+((signal-224)<<4)+7);
  return (512+((signal-224)<<4)+8);
}
//______________________________________________________________________
void AliITSsimulationSDD::Compress2D(){
  // 2D zero-suppression algorithm as described in ALICE-INT-1999-28 V10
  AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);  
  for (Int_t iWing=0; iWing<2; iWing++) {
    Int_t tL=res->GetZSLowThreshold(iWing);
    Int_t tH=res->GetZSHighThreshold(iWing);
    for (Int_t i=0; i<fNofMaps/2; i++) {  
      Int_t ian=i+iWing*fNofMaps/2;
      if( !fAnodeFire[ian] ) continue;
      for (Int_t itb=0; itb<fMaxNofSamples; itb++) {
	Int_t nLow=0, nHigh=0;      
	Float_t cC=fHitMap2->GetSignal(ian,itb);
	if(cC<=tL) continue;
	nLow++; // cC is greater than tL
	if(cC>tH) nHigh++;
	//                     N
	// Get "quintuple":   WCE
	//                     S
	Float_t wW=0.;
	if(itb>0) wW=fHitMap2->GetSignal(ian,itb-1);
	if(wW>tL) nLow++;
	if(wW>tH) nHigh++;
	Float_t eE=0.;
	if(itb<fMaxNofSamples-1) eE=fHitMap2->GetSignal(ian,itb+1);
	if(eE>tL) nLow++;
	if(eE>tH) nHigh++;
	Float_t nN=0.;
	if(i<(fNofMaps/2-1)) nN=fHitMap2->GetSignal(ian+1,itb);
	if(nN>tL) nLow++;
	if(nN>tH) nHigh++;
	Float_t sS=0.;
	if(i>0) sS=fHitMap2->GetSignal(ian-1,itb);
	if(sS>tL) nLow++;
	if(sS>tH) nHigh++;
 	
	if(nLow>=2 && nHigh>=1){
	  Int_t signal=(Int_t)cC;
	  Int_t signalc = Convert10to8(signal);
	  Int_t signale = Convert8to10(signalc);
	  signalc-=tL; // subtract low threshold after 10 to 8 bit compression
	  if(signalc>=4) AddDigit(ian,itb,signalc,signale);  // store C 
	}
      }
    }
  }
}


//______________________________________________________________________
void AliITSsimulationSDD::StoreAllDigits(){
  // store digits for non-zero-suppressed data
  for (Int_t ian=0; ian<fNofMaps; ian++) {
    for (Int_t itb=0; itb<fMaxNofSamples; itb++){
      Int_t signal=(Int_t)(fHitMap2->GetSignal(ian,itb));
      Int_t signalc = Convert10to8(signal);
      Int_t signale = Convert8to10(signalc);
      AddDigit(ian,itb,signalc,signale);  
    } 
  }
} 
//______________________________________________________________________
void AliITSsimulationSDD::CreateHistograms(Int_t scale){
  // Creates histograms of maps for debugging
  Int_t i;
  
  fHis=new TObjArray(fNofMaps);
  for (i=0;i<fNofMaps;i++) {
    TString sddName;
    sddName.Form("sdd_%d",i+1);
    fHis->AddAt(new TH1F(sddName.Data(),"SDD maps",scale*fMaxNofSamples,
			 0.,(Float_t) scale*fMaxNofSamples), i);
  } // end for i
}
//______________________________________________________________________
void AliITSsimulationSDD::FillHistograms(){
    // fill 1D histograms from map

    if (!fHis) return;

    for( Int_t i=0; i<fNofMaps; i++) {
        TH1F *hist =(TH1F *)fHis->UncheckedAt(i);
        Int_t nsamples = hist->GetNbinsX();
        for( Int_t j=0; j<nsamples; j++) {
            Double_t signal=fHitMap2->GetSignal(i,j);
            hist->Fill((Float_t)j,signal);
        } // end for j
    } // end for i
}
//______________________________________________________________________
void AliITSsimulationSDD::ResetHistograms(){
    // Reset histograms for this detector
    Int_t i;

    for (i=0;i<fNofMaps;i++ ) {
        if (fHis->At(i))    ((TH1F*)fHis->At(i))->Reset();
    } // end for i
}
//______________________________________________________________________
TH1F *AliITSsimulationSDD::GetAnode(Int_t wing, Int_t anode) { 
    // Fills a histogram from a give anode.  

    if (!fHis) return 0;

    if(wing <=0 || wing > 2) {
        Warning("GetAnode","Wrong wing number: %d",wing);
        return NULL;
    } // end if wing <=0 || wing >2
    if(anode <=0 || anode > fNofMaps/2) {
        Warning("GetAnode","Wrong anode number: %d",anode);
        return NULL;
    } // end if ampde <=0 || andoe > fNofMaps/2

    Int_t index = (wing-1)*fNofMaps/2 + anode-1;
    return (TH1F*)(fHis->At(index));
}
//______________________________________________________________________
void AliITSsimulationSDD::WriteToFile(TFile *hfile) {
    // Writes the histograms to a file

    if (!fHis) return;

    hfile->cd();
    Int_t i;
    for(i=0; i<fNofMaps; i++)  fHis->At(i)->Write(); //fAdcs[i]->Write();
    return;
}
//______________________________________________________________________
void AliITSsimulationSDD::WriteSDigits(){
    // Fills the Summable digits Tree
    static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");

    for( Int_t i=0; i<fNofMaps; i++ ) {
        if( !fAnodeFire[i] ) continue;
	for( Int_t j=0; j<fMaxNofSamples; j++ ) {
            Double_t sig = fHitMap2->GetSignal( i, j );
            if( sig > 0.2 ) {
                Int_t jdx = j*fScaleSize;
                Int_t index = fpList->GetHitIndex( i, j );
                AliITSpListItem pItemTmp2( fModule, index, 0. );
                // put the fScaleSize analog digits in only one
                for( Int_t ik=0; ik<fScaleSize; ik++ ) {
                    AliITSpListItem *pItemTmp = fpList->GetpListItem(i,jdx+ik);
                    if( pItemTmp == 0 ) continue;
                    pItemTmp2.Add( pItemTmp );
                }
                pItemTmp2.AddSignalAfterElect( fModule, index, sig );
                pItemTmp2.AddNoise(fModule,index,fHitNoiMap2->GetSignal(i,j));
                aliITS->AddSumDigit( pItemTmp2 );
            } // end if (sig > 0.2)
        }
    }
    return;
}
//______________________________________________________________________
void AliITSsimulationSDD::PrintStatus() const {
    // Print SDD simulation Parameters

    cout << "**************************************************" << endl;
    cout << "   Silicon Drift Detector Simulation Parameters   " << endl;
    cout << "**************************************************" << endl;
    cout << "Flag for Perpendicular tracks: " << (Int_t) fFlag << endl;
    cout << "Flag to switch off electronics: " << (Int_t) fDoFFT << endl;
    cout << "Number of Anodes used: " << fNofMaps << endl;
    cout << "Number of Time Samples: " << fMaxNofSamples << endl;
    cout << "Scale size factor: " << fScaleSize << endl;
    cout << "**************************************************" << endl;
}
