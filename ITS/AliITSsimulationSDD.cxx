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
#include <string.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TRandom.h>

#include "AliITS.h"
#include "AliITSMapA2.h"
#include "AliITSRawData.h"
#include "AliITSdigitSPD.h"
#include "AliITSetfSDD.h"
#include "AliITSmodule.h"
#include "AliITSpList.h"
#include "AliITSresponseSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSsegmentationSDD.h"
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
Int_t power(Int_t b, Int_t e) {
    // compute b to the e power, where both b and e are Int_ts.
    Int_t power = 1,i;

    for(i=0; i<e; i++) power *= b;
    return power;
}
//______________________________________________________________________
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
AliITSsimulationSDD::AliITSsimulationSDD():
AliITSsimulation(),
fITS(0),
fHitMap2(0),
fHitSigMap2(0),
fHitNoiMap2(0),
fStream(0),
fElectronics(0),
fInZR(0),
fInZI(0),
fOutZR(0),
fOutZI(0),
fAnodeFire(0),
fHis(0),
fD(),
fT1(),
fT2(),
fTol(),
fTreeB(0),
fParam(0),
fFileName(),
fFlag(kFALSE),
fCheckNoise(kFALSE),
fCrosstalkFlag(kFALSE),
fDoFFT(1),
fNofMaps(0),
fMaxNofSamples(0),
fScaleSize(0){
    // Default constructor
    SetScaleFourier();
    SetPerpendTracksFlag();
    SetCrosstalkFlag();
    SetDoFFT();
    SetCheckNoise();
}
//______________________________________________________________________
AliITSsimulationSDD::AliITSsimulationSDD(const AliITSsimulationSDD &source) :
    AliITSsimulation(source),
fITS(source.fITS),
fHitMap2(source.fHitMap2),
fHitSigMap2(source.fHitSigMap2),
fHitNoiMap2(source.fHitNoiMap2),
fStream(source.fStream),
fElectronics(source.fElectronics),
fInZR(source.fInZR),
fInZI(source.fInZI),
fOutZR(source.fOutZR),
fOutZI(source.fOutZI),
fAnodeFire(source.fAnodeFire),
fHis(source.fHis),
fD(source.fD),
fT1(source.fT1),
fT2(source.fT2),
fTol(source.fTol),
fTreeB(source.fTreeB),
fParam(source.fParam),
fFileName(source.fFileName),
fFlag(source.fFlag),
fCheckNoise(source.fCheckNoise),
fCrosstalkFlag(source.fCrosstalkFlag),
fDoFFT(source.fDoFFT),
fNofMaps(source.fNofMaps),
fMaxNofSamples(source.fMaxNofSamples),
fScaleSize(source.fScaleSize){
    // Copy constructor to satify Coding roules only.

}
//______________________________________________________________________
AliITSsimulationSDD& AliITSsimulationSDD::operator=(const AliITSsimulationSDD &src){
    // Assignment operator to satify Coding roules only.

    if(this==&src) return *this;
    Error("AliITSsimulationSDD","Not allowed to make a = with "
          "AliITSsimulationSDD Using default creater instead");
    return *this ;
}
//______________________________________________________________________
AliITSsimulation& AliITSsimulationSDD::operator=(const AliITSsimulation &src){
    // Assignment operator to satify Coding roules only.

    if(this==&src) return *this;
    Error("AliITSsimulationSSD","Not allowed to make a = with "
          "AliITSsimulationSDD Using default creater instead");
    return *this ;
}

//______________________________________________________________________
AliITSsimulationSDD::AliITSsimulationSDD(AliITSDetTypeSim* dettyp):
AliITSsimulation(dettyp),
fITS(0),
fHitMap2(0),
fHitSigMap2(0),
fHitNoiMap2(0),
fStream(0),
fElectronics(0),
fInZR(0),
fInZI(0),
fOutZR(0),
fOutZI(0),
fAnodeFire(0),
fHis(0),
fD(),
fT1(),
fT2(),
fTol(),
fTreeB(0),
fParam(),
fFileName(),
fFlag(kFALSE),
fCheckNoise(kFALSE),
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

    SetScaleFourier();
    SetPerpendTracksFlag();
    SetCrosstalkFlag();
    SetDoFFT();
    SetCheckNoise();

    AliITSsegmentationSDD* seg = (AliITSsegmentationSDD*)GetSegmentationModel(1);
    
    AliITSresponseSDD* res = (AliITSresponseSDD*)fDetType->GetResponse(1);
    fpList = new AliITSpList( seg->Npz(),
                              fScaleSize*seg->Npx() );
    fHitSigMap2 = new AliITSMapA2(seg,fScaleSize,1);
    fHitNoiMap2 = new AliITSMapA2(seg,fScaleSize,1);
    fHitMap2 = fHitSigMap2;

    fNofMaps = seg->Npz();
    fMaxNofSamples = seg->Npx();
    fAnodeFire = new Bool_t [fNofMaps];
    
    Float_t sddLength = seg->Dx();
    Float_t sddWidth  = seg->Dz();

    Int_t dummy        = 0;
    Float_t anodePitch = seg->Dpz(dummy);
    Double_t timeStep  = (Double_t)seg->Dpx(dummy);
    Float_t driftSpeed = res->DriftSpeed();

    if(anodePitch*(fNofMaps/2) > sddWidth) {
        Warning("AliITSsimulationSDD",
                "Too many anodes %d or too big pitch %f \n",
                fNofMaps/2,anodePitch);
    } // end if

    if(timeStep*fMaxNofSamples < sddLength/driftSpeed) {
        Error("AliITSsimulationSDD",
              "Time Interval > Allowed Time Interval: exit\n");
        return;
    } // end if

    fElectronics = new AliITSetfSDD(timeStep/fScaleSize,
                                    res->Electronics());

    char opt1[20], opt2[20];
    res->ParamOptions(opt1,opt2);
    fParam = opt2;

    const char *kopt=res->ZeroSuppOption();
    fD.Set(fNofMaps);
    fT1.Set(fNofMaps);
    fT2.Set(fNofMaps);
    fTol.Set(fNofMaps);
 
    Bool_t write = res->OutputOption();
    if(write && strstr(kopt,"2D")) MakeTreeB();
  
    fITS       = (AliITS*)gAlice->GetModule("ITS");
    Int_t size = fNofMaps*fMaxNofSamples;
    fStream    = new AliITSInStream(size);

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
    delete fStream;
    delete fElectronics;

    fITS = 0;

    if (fHis) {
        fHis->Delete(); 
        delete fHis;     
    } // end if fHis
    if(fTreeB) delete fTreeB;           
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
void AliITSsimulationSDD::SDigitiseModule(AliITSmodule *mod,Int_t md,Int_t ev){
    // digitize module using the "slow" detector simulator creating
    // summable digits.

    TObjArray *fHits = mod->GetHits();
    Int_t nhits      = fHits->GetEntriesFast();
    if( !nhits ) return;

    InitSimulationModule( md, ev );
    HitsToAnalogDigits( mod );
    ChargeToSignal( fModule,kFALSE ); // - Process signal without add noise
    fHitMap2 = fHitNoiMap2;   // - Swap to noise map
    ChargeToSignal( fModule,kTRUE );  // - Process only noise
    fHitMap2 = fHitSigMap2;   // - Return to signal map
    WriteSDigits();
    ClearMaps();
}
//______________________________________________________________________
Bool_t AliITSsimulationSDD::AddSDigitsToModule(TClonesArray *pItemArray,
                                               Int_t mask ) {
    // Add Summable digits to module maps.
   AliITSresponseSDD* res = (AliITSresponseSDD*)fDetType->GetResponse(1);
    Int_t    nItems = pItemArray->GetEntries();
    Double_t maxadc = res->MaxAdc();
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

    if( !nhits && fCheckNoise ) {
        ChargeToSignal( fModule,kTRUE );  // process noise
        GetNoise();
        ClearMaps();
        return;
    } else 
        if( !nhits ) return;
        
    HitsToAnalogDigits( mod );
    ChargeToSignal( fModule,kTRUE );  // process signal + noise

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

    ApplyDeadChannels(fModule);
    if( fCrosstalkFlag ) ApplyCrosstalk(fModule);

    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);
    const char *kopt = res->GetZeroSuppOption();
    ZeroSuppression( kopt );
}
//______________________________________________________________________
void AliITSsimulationSDD::HitsToAnalogDigits( AliITSmodule *mod ) {
    // create maps to build the lists of tracks for each digit

  AliITSsegmentationSDD* seg = (AliITSsegmentationSDD*)GetSegmentationModel(1);
  AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);
  TObjArray *hits     = mod->GetHits();
    Int_t      nhits    = hits->GetEntriesFast();

    //    Int_t      arg[6]   = {0,0,0,0,0,0};
    Int_t     dummy      = 0;
    Int_t     nofAnodes  = fNofMaps/2;
    Double_t  sddLength  = seg->Dx();
    Double_t  sddWidth   = seg->Dz();
    Double_t  anodePitch = seg->Dpz(dummy);
    Double_t  timeStep   = seg->Dpx(dummy);
    Double_t  driftSpeed = res->GetDriftSpeed();
    //Float_t   maxadc     = res->GetMaxAdc();    
    //Float_t   topValue   = res->GetDynamicRange();
    Double_t  norm       = res->GetMaxAdc()/res->GetDynamicRange(); //   maxadc/topValue;
    Double_t  cHloss     = res->GetChargeLoss();
    Float_t   dfCoeff, s1; res->DiffCoeff(dfCoeff,s1); // Signal 2d Shape
    Double_t  eVpairs    = res->GetGeVToCharge()*1.0E9; // 3.6 eV by def.
    Double_t  nsigma     = res->GetNSigmaIntegration(); //
    Int_t     nlookups   = res->GetGausNLookUp();       //
    Float_t   jitter     = res->GetJitterError(); // 

    // Piergiorgio's part (apart for few variables which I made float
    // when i thought that can be done
    // Fill detector maps with GEANT hits
    // loop over hits in the module

    const Float_t kconv = 1.0e+6;  // GeV->KeV
    Int_t     itrack      = 0;
    Int_t     hitDetector; // detector number (lay,lad,hitDetector)
    Int_t     iWing;       // which detector wing/side.
    Int_t     detector;    // 2*(detector-1)+iWing
    Int_t     ii,kk,ka,kt; // loop indexs
    Int_t     ia,it,index; // sub-pixel integration indexies
    Int_t     iAnode;      // anode number.
    Int_t     timeSample;  // time buckett.
    Int_t     anodeWindow; // anode direction charge integration width
    Int_t     timeWindow;  // time direction charge integration width
    Int_t     jamin,jamax; // anode charge integration window
    Int_t     jtmin,jtmax; // time charge integration window
    Int_t     ndiv;        // Anode window division factor.
    Int_t     nsplit;      // the number of splits in anode and time windows==1.
    Int_t     nOfSplits;   // number of times track length is split into
    Float_t   nOfSplitsF;  // Floating point version of nOfSplits.
    Float_t   kkF;         // Floating point version of loop index kk.
    Double_t  pathInSDD; // Track length in SDD.
    Double_t  drPath; // average position of track in detector. in microns
    Double_t  drTime; // Drift time
    Double_t  nmul;   // drift time window multiplication factor.
    Double_t  avDrft;  // x position of path length segment in cm.
    Double_t  avAnode; // Anode for path length segment in Anode number (float)
    Double_t  xAnode;  // Floating point anode number.
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
    //  Double_t tof;            // Time of flight in ns of this step.    

    for(ii=0; ii<nhits; ii++) {
        if(!mod->LineSegmentL(ii,xL[0],dxL[0],xL[1],dxL[1],xL[2],dxL[2],
                              depEnergy,itrack)) continue;
        xL[0] += 0.0001*gRandom->Gaus( 0, jitter ); //
        depEnergy  *= kconv;
        hitDetector = mod->GetDet();
        //tof         = 1.E+09*(mod->GetHit(ii)->GetTOF()); // tof in ns.
        //if(tof>sddLength/driftSpeed) continue; // hit happed too late.
        
        // scale path to simulate a perpendicular track
        // continue if the particle did not lose energy
        // passing through detector
        if (!depEnergy) {
            AliDebug(1,
             Form("fTrack = %d hit=%d module=%d This particle has passed without losing energy!",
		  itrack,ii,mod->GetIndex()));
            continue;
        } // end if !depEnergy

        pathInSDD = TMath::Sqrt(dxL[0]*dxL[0]+dxL[1]*dxL[1]+dxL[2]*dxL[2]);

        if (fFlag && pathInSDD) { depEnergy *= (0.03/pathInSDD); }
        drPath = 10000.*(dxL[0]+2.*xL[0])*0.5;
        if(drPath < 0) drPath = -drPath;
        drPath = sddLength-drPath;
        if(drPath < 0) {
            AliDebug(1, // this should be fixed at geometry level
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
        for(kk=0;kk<nOfSplits;kk++) { // loop over path segments
            kkF       = (Float_t) kk + 0.5;
            avDrft    = xL[0]+dxL[0]*kkF/nOfSplitsF;
            avAnode   = xL[2]+dxL[2]*kkF/nOfSplitsF;
            driftPath = 10000.*avDrft;

            iWing = 2;  // Assume wing is 2
            if(driftPath < 0) { // if wing is not 2 it is 1.
                iWing     = 1;
                driftPath = -driftPath;
            } // end if driftPath < 0
            driftPath = sddLength-driftPath;
            detector  = 2*(hitDetector-1) + iWing;
            if(driftPath < 0) {
                AliDebug(1, // this should be fixed at geometry level
                 Form("negative drift path driftPath=%e sddLength=%e avDrft=%e dxL[0]=%e xL[0]=%e",
		      driftPath,sddLength,avDrft,dxL[0],xL[0]));
                continue;
            } // end if driftPath < 0

            //   Drift Time
            drTime     = driftPath/driftSpeed; // drift time for segment.
            timeSample = (Int_t) (fScaleSize*drTime/timeStep + 1);
            // compute time Sample including tof information. The tof only 
            // effects the time of the signal is recoreded and not the
            // the defusion.
            // timeSample = (Int_t) (fScaleSize*(drTime+tof)/timeStep + 1);
            if(timeSample > fScaleSize*fMaxNofSamples) {
                Warning("HitsToAnalogDigits","Wrong Time Sample: %e",
                        timeSample);
                continue;
            } // end if timeSample > fScaleSize*fMaxNoofSamples

            //   Anode
            xAnode = 10000.*(avAnode)/anodePitch + nofAnodes/2;  // +1?
            if(xAnode*anodePitch > sddWidth || xAnode*anodePitch < 0.) 
                Warning("HitsToAnalogDigits",
                        "Exceedubg sddWidth=%e Z = %e",
                        sddWidth,xAnode*anodePitch);
            iAnode = (Int_t) (1.+xAnode); // xAnode?
            if(iAnode < 1 || iAnode > nofAnodes) {
                Warning("HitToAnalogDigits","Wrong iAnode: 1<%d>%d  (xanode=%e)",
                        iAnode,nofAnodes, xAnode);
                continue;
            } // end if iAnode < 1 || iAnode > nofAnodes

            // store straight away the particle position in the array
            // of particles and take idhit=ii only when part is entering (this
            // requires FillModules() in the macro for analysis) :
    
            // Sigma along the anodes for track segment.
            sigA       = TMath::Sqrt(2.*dfCoeff*drTime+s1*s1);
            sigT       = sigA/driftSpeed;
            // Peak amplitude in nanoAmpere
            amplitude  = fScaleSize*160.*depEnergy/
                (timeStep*eVpairs*2.*acos(-1.)*sigT*sigA);
            amplitude *= timeStep/25.; // WARNING!!!!! Amplitude scaling to 
            // account for clock variations 
            // (reference value: 40 MHz)
            chargeloss = 1.-cHloss*driftPath/1000;
            amplitude *= chargeloss;
            width  = 2.*nsigma/(nlookups-1);
            // Spread the charge 
            // Pixel index
            ndiv = 2;
            nmul = 3.; 
            if(drTime > 1200.) { 
                ndiv = 4;
                nmul = 1.5;
            } // end if drTime > 1200.
            // Sub-pixel index
            nsplit = 4; // hard-wired //nsplit=4;nsplit = (nsplit+1)/2*2;
            // Sub-pixel size see computation of aExpo and tExpo.
            aStep  = anodePitch/(nsplit*fScaleSize*sigA);
            aConst = xAnode*anodePitch/sigA;
            tStep  = timeStep/(nsplit*fScaleSize*sigT);
            tConst = drTime/sigT;
            // Define SDD window corresponding to the hit
            anodeWindow = (Int_t)(fScaleSize*nsigma*sigA/anodePitch+1);
            timeWindow  = (Int_t) (fScaleSize*nsigma*sigT/timeStep+1.);
            jamin = (iAnode - anodeWindow/ndiv - 1)*fScaleSize*nsplit +1;
            jamax = (iAnode + anodeWindow/ndiv)*fScaleSize*nsplit;
            if(jamin <= 0) jamin = 1;
            if(jamax > fScaleSize*nofAnodes*nsplit) 
                jamax = fScaleSize*nofAnodes*nsplit;
            // jtmin and jtmax are Hard-wired
            jtmin = (Int_t)(timeSample-timeWindow*nmul-1)*nsplit+1;
            jtmax = (Int_t)(timeSample+timeWindow*nmul)*nsplit;
            if(jtmin <= 0) jtmin = 1;
            if(jtmax > fScaleSize*fMaxNofSamples*nsplit) 
                jtmax = fScaleSize*fMaxNofSamples*nsplit;
            // Spread the charge in the anode-time window
            for(ka=jamin; ka <=jamax; ka++) {
                ia = (ka-1)/(fScaleSize*nsplit) + 1;
                if(ia <= 0) {
                    Warning("HitsToAnalogDigits","ia < 1: ");
                    continue;
                } // end if
                if(ia > nofAnodes) ia = nofAnodes;
                aExpo     = (aStep*(ka-0.5)-aConst);
                if(TMath::Abs(aExpo) > nsigma)  anodeAmplitude = 0.;
                else {
                    dummy          = (Int_t) ((aExpo+nsigma)/width);
                    anodeAmplitude = amplitude*res->GetGausLookUp(dummy);
                } // end if TMath::Abs(aEspo) > nsigma
                // index starts from 0
                index = ((detector+1)%2)*nofAnodes+ia-1;
                if(anodeAmplitude) for(kt=jtmin; kt<=jtmax; kt++) {
                    it = (kt-1)/nsplit+1;  // it starts from 1
                    if(it<=0){
                        Warning("HitsToAnalogDigits","it < 1:");
                        continue;
                    } // end if 
                    if(it>fScaleSize*fMaxNofSamples)
                        it = fScaleSize*fMaxNofSamples;
                    tExpo    = (tStep*(kt-0.5)-tConst);
                    if(TMath::Abs(tExpo) > nsigma) timeAmplitude = 0.;
                    else {
                        dummy         = (Int_t) ((tExpo+nsigma)/width);
                        timeAmplitude = anodeAmplitude*
                            res->GetGausLookUp(dummy);
                    } // end if TMath::Abs(tExpo) > nsigma
                    // build the list of Sdigits for this module        
                    //                    arg[0]     = index;
                    //                    arg[1]     = it;
                    //                    arg[2]     = itrack; // track number
                    //                    arg[3]     = ii-1; // hit number.
                    timeAmplitude *= norm;
                    timeAmplitude *= 10;
                    //         ListOfFiredCells(arg,timeAmplitude,alst,padr);
                    Double_t charge = timeAmplitude;
                    charge += fHitMap2->GetSignal(index,it-1);
                    fHitMap2->SetHit(index, it-1, charge);
                    fpList->AddSignal(index,it-1,itrack,ii-1,
                                      mod->GetIndex(),timeAmplitude);
                    fAnodeFire[index] = kTRUE;                 
                } // end if anodeAmplitude and loop over time in window
            } // loop over anodes in window
        } // end loop over "sub-hits"
    } // end loop over hits
}
/*
//______________________________________________________________________
void AliITSsimulationSDD::ListOfFiredCells(Int_t *arg,Double_t timeAmplitude,
					   TObjArray *alist,TClonesArray *padr){
    // Returns the list of "fired" cells.

    Int_t index     = arg[0];
    Int_t ik        = arg[1];
    Int_t idtrack   = arg[2];
    Int_t idhit     = arg[3];
    Int_t counter   = arg[4];
    Int_t countadr  = arg[5];
    Double_t charge = timeAmplitude;
    charge += fHitMap2->GetSignal(index,ik-1);
    fHitMap2->SetHit(index, ik-1, charge);

    Int_t digits[3];
    Int_t it = (Int_t)((ik-1)/fScaleSize);
    digits[0] = index;
    digits[1] = it;
    digits[2] = (Int_t)timeAmplitude;
    Float_t phys;
    if (idtrack >= 0) phys = (Float_t)timeAmplitude;
    else phys = 0;

    Double_t cellcharge = 0.;
    AliITSTransientDigit* pdigit;
    // build the list of fired cells and update the info
    if (!fHitMap1->TestHit(index, it)) {
        new((*padr)[countadr++]) TVector(3);
        TVector &trinfo=*((TVector*) (*padr)[countadr-1]);
        trinfo(0) = (Float_t)idtrack;
        trinfo(1) = (Float_t)idhit;
        trinfo(2) = (Float_t)timeAmplitude;
      
        alist->AddAtAndExpand(new AliITSTransientDigit(phys,digits),counter);
        fHitMap1->SetHit(index, it, counter);
        counter++;
        pdigit=(AliITSTransientDigit*)alist->At(alist->GetLast());
        // list of tracks
        TObjArray *trlist=(TObjArray*)pdigit->TrackList();
        trlist->Add(&trinfo);
    } else {
        pdigit = (AliITSTransientDigit*) fHitMap1->GetHit(index, it);
        for(Int_t kk=0;kk<fScaleSize;kk++) {
            cellcharge += fHitMap2->GetSignal(index,fScaleSize*it+kk);
        }  // end for kk
        // update charge
        (*pdigit).fSignal = (Int_t)cellcharge;
        (*pdigit).fPhysics += phys;                        
        // update list of tracks
        TObjArray* trlist = (TObjArray*)pdigit->TrackList();
        Int_t lastentry = trlist->GetLast();
        TVector *ptrkp = (TVector*)trlist->At(lastentry);
        TVector &trinfo = *ptrkp;
        Int_t lasttrack = Int_t(trinfo(0));
        Float_t lastcharge=(trinfo(2));
        if (lasttrack==idtrack ) {
            lastcharge += (Float_t)timeAmplitude;
            trlist->RemoveAt(lastentry);
            trinfo(0) = lasttrack;
            trinfo(1) = idhit;
            trinfo(2) = lastcharge;
            trlist->AddAt(&trinfo,lastentry);
        } else {                  
            new((*padr)[countadr++]) TVector(3);
            TVector &trinfo=*((TVector*) (*padr)[countadr-1]);
            trinfo(0) = (Float_t)idtrack;
            trinfo(1) = (Float_t)idhit;
            trinfo(2) = (Float_t)timeAmplitude;
            trlist->Add(&trinfo);
        } // end if lasttrack==idtrack

        if(AliDebugLevel()){
            // check the track list - debugging
            Int_t trk[20], htrk[20];
            Float_t chtrk[20];  
            Int_t nptracks = trlist->GetEntriesFast();
            if (nptracks > 2) {
                Int_t tr;
                for (tr=0;tr<nptracks;tr++) {
                    TVector *pptrkp = (TVector*)trlist->At(tr);
                    TVector &pptrk  = *pptrkp;
                    trk[tr]   = Int_t(pptrk(0));
                    htrk[tr]  = Int_t(pptrk(1));
                    chtrk[tr] = (pptrk(2));
                    cout << "nptracks "<<nptracks << endl;
                } // end for tr
            } // end if nptracks
        } // end if AliDebugLevel()
    } //  end if pdigit

    // update counter and countadr for next call.
    arg[4] = counter;
    arg[5] = countadr;
}
*/
//____________________________________________
void AliITSsimulationSDD::AddDigit( Int_t i, Int_t j, Int_t signal ) {
    // Adds a Digit.
    Int_t size = AliITSdigitSPD::GetNTracks();
    Int_t digits[3];
    Int_t * tracks = new Int_t[size];
    Int_t * hits = new Int_t[size];
    Float_t phys;
    Float_t * charges = new Float_t[size];

    digits[0] = i;
    digits[1] = j;
    digits[2] = signal;

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

    fITS->AddSimDigit( 1, phys, digits, tracks, hits, charges ); 
    delete [] tracks;
    delete [] hits;
    delete [] charges;
}
//______________________________________________________________________
void AliITSsimulationSDD::ChargeToSignal(Int_t mod,Bool_t bAddNoise) {
    // add baseline, noise, electronics and ADC saturation effects

    char opt1[20], opt2[20];
    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(mod);
    res->GetParamOptions(opt1,opt2);
    Double_t baseline=0; 
    Double_t noise=0; 

    Float_t contrib=0;
    Int_t i,k,kk;
    Float_t maxadc = res->GetMaxAdc();    
    if(!fDoFFT) {
        for (i=0;i<fNofMaps;i++) {
            if( !fAnodeFire[i] ) continue;
	    baseline = res->GetBaseline(i);
	    noise = res->GetNoise(i);
	    
            for(k=0; k<fScaleSize*fMaxNofSamples; k++) {
                fInZR[k]  = fHitMap2->GetSignal(i,k);
                if( bAddNoise ) {
                    contrib   = (baseline + noise*gRandom->Gaus());
                    fInZR[k] += contrib;
                }
            } // end for k
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
                    Warning("","newcont=%d>=baseline=%d",newcont,baseline);
                } // end if
                // back to analog: ?
                fHitMap2->SetHit(i,k,newcont);
            }  // end for k
        } // end for i loop over anodes
        return;
    } // end if DoFFT

    for (i=0;i<fNofMaps;i++) {
        if( !fAnodeFire[i] ) continue;
	baseline = res->GetBaseline(i);
	noise = res->GetNoise(i);
	for(k=0; k<fScaleSize*fMaxNofSamples; k++) {
            fInZR[k]  = fHitMap2->GetSignal(i,k);
            if( bAddNoise ) {
                contrib   = (baseline + noise*gRandom->Gaus());
                fInZR[k] += contrib;
            }
            fInZI[k]  = 0.;
        } // end for k
        FastFourierTransform(fElectronics,&fInZR[0],&fInZI[0],1);
        for(k=0; k<fScaleSize*fMaxNofSamples; k++) {
            Double_t rw = fElectronics->GetTraFunReal(k);
            Double_t iw = fElectronics->GetTraFunImag(k);
            fOutZR[k]   = fInZR[k]*rw - fInZI[k]*iw;
            fOutZI[k]   = fInZR[k]*iw + fInZI[k]*rw;
        } // end for k
        FastFourierTransform(fElectronics,&fOutZR[0],&fOutZI[0],-1);
        for(k=0; k<fMaxNofSamples; k++) {
            Double_t newcont1 = 0.;
            Double_t maxcont1 = 0.;
            for(kk=0;kk<fScaleSize;kk++) {
                newcont1 = fOutZR[fScaleSize*k+kk];
                if(newcont1 > maxcont1) maxcont1 = newcont1;
            } // end for kk
            newcont1 = maxcont1;
            if (newcont1 >= maxadc) newcont1 = maxadc -1;
            fHitMap2->SetHit(i,k,newcont1);
        } // end for k
    } // end for i loop over anodes
    return;
}
//____________________________________________________________________
void AliITSsimulationSDD::ApplyDeadChannels(Int_t mod) {    
    // Set dead channel signal to zero
    AliITSCalibrationSDD * calibr = (AliITSCalibrationSDD *)GetCalibrationModel(mod);
    AliITSsegmentationSDD* seg = (AliITSsegmentationSDD*)GetSegmentationModel(1);
    // nothing to do
    if( calibr->IsDead() ||   
        ( calibr->GetDeadChips() == 0 &&
	  calibr->GetDeadChannels() == 0 ) ) return;  
    
    // static AliITS *iTS = (AliITS*)gAlice->GetModule( "ITS" );

    Int_t fMaxNofSamples = seg->Npx();    
    // AliITSgeom *geom = iTS->GetITSgeom();
    // Int_t firstSDDMod = geom->GetStartDet( 1 );
    // loop over wings
    for( Int_t j=0; j<2; j++ ) {
      // Int_t mod = (fModule-firstSDDMod)*2 + j;
      for( Int_t u=0; u<calibr->Chips(); u++ )
	for( Int_t v=0; v<calibr->Channels(); v++ ) {
	  Float_t gain = calibr->Gain(j, u, v );
	  for( Int_t k=0; k<fMaxNofSamples; k++ ) {
	    Int_t i = j*calibr->Chips()*calibr->Channels() +
	      u*calibr->Channels() + 
	      v;
	    Double_t signal =  gain * fHitMap2->GetSignal( i, k );
	    fHitMap2->SetHit( i, k, signal );  ///
	  }
	}
    } 
}
//______________________________________________________________________
void AliITSsimulationSDD::ApplyCrosstalk(Int_t mod) {
    // function add the crosstalk effect to signal
    // temporal function, should be checked...!!!
    AliITSsegmentationSDD* seg = (AliITSsegmentationSDD*)GetSegmentationModel(1);
  
    Int_t fNofMaps = seg->Npz();
    Int_t fMaxNofSamples = seg->Npx();

    // create and inizialice crosstalk map
    Float_t* ctk = new Float_t[fNofMaps*fMaxNofSamples+1];
    if( ctk == NULL ) {
        Error( "ApplyCrosstalk", "no memory for temporal map: exit \n" );
        return;
    }
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
                        if( ctk == NULL ) {
                            Error( "ApplyCrosstalk", 
                                   "no memory for temporal array: exit \n" );
                            return;
                        }
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
void AliITSsimulationSDD::CompressionParam(Int_t i,Int_t &db,Int_t &tl,
                                           Int_t &th) const{
    // Returns the compression alogirthm parameters
  db=fD[i]; 
  tl=fT1[i]; 
  th=fT2[i];
}
//______________________________________________________________________
void AliITSsimulationSDD::CompressionParam(Int_t i,Int_t &db,Int_t &tl) const{
    // returns the compression alogirthm parameters

    db=fD[i]; 
    tl=fT1[i];
 
}
//______________________________________________________________________
void AliITSsimulationSDD::SetCompressParam(){ 
    // Sets the compression alogirthm parameters  
    AliITSCalibrationSDD* calibr = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);
    for(Int_t ian = 0; ian<fNofMaps;ian++){
      fD[ian] = (Int_t)(calibr->GetBaseline(ian));
      fT1[ian] = (Int_t)(2.*calibr->GetNoiseAfterElectronics(ian)+0.5);
      fT2[ian] = 0;   // used by 2D clustering - not defined yet
      fTol[ian] = 0; // used by 2D clustering - not defined yet
    }
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
/*
//______________________________________________________________________
AliITSMap*   AliITSsimulationSDD::HitMap(Int_t i){
    //Return the correct map.

    return ((i==0)? fHitMap1 : fHitMap2);
}
*/
//______________________________________________________________________
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
//______________________________________________________________________
void AliITSsimulationSDD::Init2D(){
    // read in and prepare arrays: fD, fT1, fT2
    //                         savemu[nanodes], savesigma[nanodes] 
    // read baseline and noise from file - either a .root file and in this
    // case data should be organised in a tree with one entry for each
    // module => reading should be done accordingly
    // or a classic file and do smth. like this ( code from Davide C. and
    // Albert W.) :
    // Read 2D zero-suppression parameters for SDD

    if (!strstr(fParam.Data(),"file")) return;

    Int_t na,pos,tempTh;
    Float_t mu,sigma;
    Float_t *savemu    = new Float_t [fNofMaps];
    Float_t *savesigma = new Float_t [fNofMaps];
    char input[100],basel[100],par[100];
    char *filtmp;
    Double_t tmp1,tmp2;
    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);

    res->Thresholds(tmp1,tmp2);
    Int_t minval = static_cast<Int_t>(tmp1);

    res->Filenames(input,basel,par);
    fFileName = par;
    //
    filtmp = gSystem->ExpandPathName(fFileName.Data());
    FILE *param = fopen(filtmp,"r");
    na = 0;

    if(param) {
        while(fscanf(param,"%d %f %f",&pos, &mu, &sigma) != EOF) {
            if (pos != na+1) {
                Error("Init2D","Anode number not in increasing order!",filtmp);
                exit(1);
            } // end if pos != na+1
            savemu[na] = mu;
            savesigma[na] = sigma;
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
        Error("Init2D","THE FILE %s DOES NOT EXIST !",filtmp);
        exit(1);
    } // end if(param)

    fclose(param);
    delete [] filtmp;
    delete [] savemu;
    delete [] savesigma;
}
//______________________________________________________________________
void AliITSsimulationSDD::Compress2D(){
    // simple ITS cluster finder -- online zero-suppression conditions

    Int_t db,tl,th; 
    Double_t tmp1,tmp2;
    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);

    res->Thresholds(tmp1,tmp2); 
    Int_t minval   = static_cast<Int_t>(tmp1);
    Bool_t write   = res->OutputOption();   
    Bool_t do10to8 = res->Do10to8();
    Int_t nz, nl, nh, low, i, j; 
    SetCompressParam();
    for (i=0; i<fNofMaps; i++) {
        CompressionParam(i,db,tl,th);
        nz  = 0; 
        nl  = 0;
        nh  = 0;
        low = 0;
        for (j=0; j<fMaxNofSamples; j++) {
            Int_t signal=(Int_t)(fHitMap2->GetSignal(i,j));
            signal -= db; // if baseline eq. is done here
            if (signal <= 0) {nz++; continue;}
            if ((signal - tl) < minval) low++;
            if ((signal - th) >= minval) {
                nh++;
                Bool_t cond=kTRUE;
                FindCluster(i,j,signal,minval,cond);
                if(cond && j &&
                   ((TMath::Abs(fHitMap2->GetSignal(i,j-1))-th)>=minval)){
                    if(do10to8) signal = Convert10to8(signal);
                    AddDigit(i,j,signal);
                } // end if cond&&j&&()
            } else if ((signal - tl) >= minval) nl++;
        } // end for j loop time samples
        if (write) TreeB()->Fill(nz,nl,nh,low,i+1);
    } //end for i loop anodes

    char hname[30];
    if (write) {
        sprintf(hname,"TNtuple%d_%d",fModule,fEvent);
        TreeB()->Write(hname);
        // reset tree
        TreeB()->Reset();
    } // end if write
}
//______________________________________________________________________
void  AliITSsimulationSDD::FindCluster(Int_t i,Int_t j,Int_t signal,
                                       Int_t minval,Bool_t &cond){
    // Find clusters according to the online 2D zero-suppression algorithm
    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);
    AliITSsegmentationSDD* seg = (AliITSsegmentationSDD*)GetSegmentationModel(1);
  
    Bool_t do10to8 = res->Do10to8();
    Bool_t high    = kFALSE;

    fHitMap2->FlagHit(i,j);
    //
    //  check the online zero-suppression conditions
    //  
    const Int_t kMaxNeighbours = 4;
    Int_t nn;
    Int_t dbx,tlx,thx;  
    Int_t xList[kMaxNeighbours], yList[kMaxNeighbours];
    seg->Neighbours(i,j,&nn,xList,yList);
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
                } // end if cond
                if(do10to8) qns = Convert10to8(qn);
                else qns=qn;
                if (!high) AddDigit(ix,iy,qns);
                cond=kFALSE;
                if(!high) fHitMap2->FlagHit(ix,iy);
            } // end if qn-tlx < minval
        } // end if  TestHit
    } // end for in loop over neighbours
}
//______________________________________________________________________
void AliITSsimulationSDD::Init1D(){
    // this is just a copy-paste of input taken from 2D algo
    // Torino people should give input
    // Read 1D zero-suppression parameters for SDD
    
    if (!strstr(fParam.Data(),"file")) return;

    Int_t na,pos,tempTh;
    Float_t mu,sigma;
    Float_t *savemu    = new Float_t [fNofMaps];
    Float_t *savesigma = new Float_t [fNofMaps];
    char input[100],basel[100],par[100];
    char *filtmp;
    Double_t tmp1,tmp2;
    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);

    res->Thresholds(tmp1,tmp2);
    Int_t minval = static_cast<Int_t>(tmp1);

    res->Filenames(input,basel,par);
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
                Error("Init1D","Anode number not in increasing order!",filtmp);
                exit(1);
            } // end if pos != na+1
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
        Error("Init1D","THE FILE %s DOES NOT EXIST !",filtmp);
        exit(1);
    } // end if(param)

    fclose(param);
    delete [] filtmp;
    delete [] savemu;
    delete [] savesigma;
} 
//______________________________________________________________________
void AliITSsimulationSDD::Compress1D(){
    // 1D zero-suppression algorithm (from Gianluca A.)
    Int_t    dis,tol,thres,decr,diff;
    UChar_t *str=fStream->Stream();
    Int_t    counter=0;
    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);

    Bool_t   do10to8=res->Do10to8();
    Int_t    last=0;
    Int_t    k,i,j;
    SetCompressParam();
    for (k=0; k<2; k++) {
        tol = Tolerance(k);
        dis = Disable(k);  
        for (i=0; i<fNofMaps/2; i++) {
            Bool_t firstSignal=kTRUE;
            Int_t idx=i+k*fNofMaps/2;
            if( !fAnodeFire[idx] ) continue;
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
                } // end if signal <= thres
                diff=signal-last;
                if (diff > 127) diff=127;
                if (diff < -128) diff=-128;
                if (signal < dis) {
                    // tol has changed to 8 possible cases ? - one can write
                    // this if(TMath::Abs(diff)<tol) ... else ...
                    if(TMath::Abs(diff)<tol) diff=0;
                    // or keep it as it was before
                    AddDigit(idx,j,last+diff);
                } else {
                    AddDigit(idx,j,signal);
                } // end if singal < dis
                diff += 128;
                // write diff in the buffer used to compute Huffman tables
                if (firstSignal) str[counter]=(UChar_t)signal;
                else str[counter]=(UChar_t)diff;
                counter++;
                last=signal;
                firstSignal=kFALSE;
            } // end for j loop time samples
        } // end for i loop anodes  one half of detector 
    } //  end for k

    // check
    fStream->CheckCount(counter);

    // open file and write out the stream of diff's
    static Bool_t open=kTRUE;
    static TFile *outFile;
    Bool_t write = res->OutputOption();
    TDirectory *savedir = gDirectory;
 
    if (write ) {
        if(open) {
            SetFileName("stream.root");
            cout<<"filename "<<fFileName<<endl;
            outFile=new TFile(fFileName,"recreate");
            cout<<"I have opened "<<fFileName<<" file "<<endl;
        } // end if open
        open = kFALSE;
        outFile->cd();
        fStream->Write();
    }  // endif write

    fStream->ClearStream();

    // back to galice.root file
    if(savedir) savedir->cd();
}
//______________________________________________________________________
void AliITSsimulationSDD::StoreAllDigits(){
    // if non-zero-suppressed data
    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);

    Bool_t do10to8 = res->Do10to8();
    Int_t i, j, digits[3];

    for (i=0; i<fNofMaps; i++) {
        for (j=0; j<fMaxNofSamples; j++) {
            Int_t signal=(Int_t)(fHitMap2->GetSignal(i,j));
            if(do10to8) signal = Convert10to8(signal);
            digits[0] = i;
            digits[1] = j;
            digits[2] = signal;
            fITS->AddRealDigit(1,digits);
        } // end for j
    } // end for i
} 
//______________________________________________________________________
void AliITSsimulationSDD::CreateHistograms(Int_t scale){
    // Creates histograms of maps for debugging
    Int_t i;

    fHis=new TObjArray(fNofMaps);
    for (i=0;i<fNofMaps;i++) {
        TString sddName("sdd_");
        Char_t candNum[4];
        sprintf(candNum,"%d",i+1);
        sddName.Append(candNum);
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
Float_t AliITSsimulationSDD::GetNoise() {  
    // Returns the noise value
    //Bool_t do10to8=GetResp()->Do10to8();
    //noise will always be in the liniar part of the signal
    Int_t decr;
    Int_t threshold = fT1[0];
    char opt1[20], opt2[20];
    AliITSCalibrationSDD* res = (AliITSCalibrationSDD*)GetCalibrationModel(fModule);
    SetCompressParam();
    res->GetParamOptions(opt1,opt2);
    fParam=opt2;
    Double_t noise,baseline;
    //GetBaseline(fModule);

    TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
    if(c2) delete c2->GetPrimitive("noisehist");
    if(c2) delete c2->GetPrimitive("anode");
    else     c2=new TCanvas("c2");
    c2->cd();
    c2->SetFillColor(0);

    TH1F *noisehist = new TH1F("noisehist","noise",100,0.,(float)2*threshold);
    TH1F *anode = new TH1F("anode","Anode Projection",fMaxNofSamples,0.,
                           (float)fMaxNofSamples);
    Int_t i,k;
    for (i=0;i<fNofMaps;i++) {
        CompressionParam(i,decr,threshold); 
	baseline = res->GetBaseline(i);
	noise = res->GetNoise(i);
        anode->Reset();
        for (k=0;k<fMaxNofSamples;k++) {
            Float_t signal=(Float_t)fHitMap2->GetSignal(i,k);
            //if (signal <= (float)threshold) noisehist->Fill(signal-baseline);
            if (signal <= (float)(threshold+decr)) noisehist->Fill(signal);
            anode->Fill((float)k,signal);
        } // end for k
        anode->Draw();
        c2->Update();
    } // end for i
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
    cout << "Flag for noise checking: " << (Int_t) fCheckNoise << endl;
    cout << "Flag to switch off electronics: " << (Int_t) fDoFFT << endl;
    cout << "Number pf Anodes used: " << fNofMaps << endl;
    cout << "Number of Time Samples: " << fMaxNofSamples << endl;
    cout << "Scale size factor: " << fScaleSize << endl;
    cout << "**************************************************" << endl;
}
