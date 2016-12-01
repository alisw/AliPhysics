
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

#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TString.h>
#include "AliITSMFTDigitPix.h"
#include "AliITSMFTHit.h"
#include "AliITSMFTChip.h"
#include "AliITSMFTSensMap.h"
#include "AliITSMFTCalibrationPix.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliLog.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliMathBase.h"
#include "AliITSMFTSDigit.h"
#include "AliITSMFTParamList.h"
#include "AliITSMFTAux.h"
#include "AliMC.h"
#include "TParticle.h"
#include "AliITSMFTSimulationPix.h"
#include "AliITSMFTGeomTGeo.h"

using namespace TMath;



ClassImp(AliITSMFTSimulationPix)
////////////////////////////////////////////////////////////////////////
//  Version: 1
//  Modified by D. Elia, G.E. Bruno, H. Tydesjo
//  Fast diffusion code by Bjorn S. Nilsen
//  March-April 2006
//  October     2007: GetCalibrationObjects() removed
//
//  Version: 0
//  Written by Boris Batyunya
//  December 20 1999
//
//  Adapted for pixels of ITS upgrade July 2012, ruben.shahoyan@cern.ch
//
//  AliITSMFTSimulationPix is to do the simulation of pixels
//
//  2013 Feb: Added MonoPix response and nois calculation al la MIMOSA32 (levente.molnar@cern.ch)
//
////////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSMFTSimulationPix::AliITSMFTSimulationPix()
:  fTanLorAng(0)
,fGlobalChargeScale(1.0)
,fSpread2DHisto(0)
,fSpreadFun(0)
,fROTimeFun(0)
{
    // Default constructor.
    SetUniqueID(AliITSMFTAux::kChipTypePix);
}

//______________________________________________________________________
AliITSMFTSimulationPix::AliITSMFTSimulationPix(AliITSMFTSimuParam* sim,AliITSMFTSensMap* map)
:AliITSMFTSimulation(sim,map)
,fTanLorAng(0)
,fGlobalChargeScale(1.0)
,fSpread2DHisto(0)
,fSpreadFun(0)
,fROTimeFun(0)
{
    // standard constructor
    SetUniqueID(AliITSMFTAux::kChipTypePix);
    Init();
}

//______________________________________________________________________
AliITSMFTSimulationPix::AliITSMFTSimulationPix(const AliITSMFTSimulationPix &s)
:AliITSMFTSimulation(s)
,fTanLorAng(s.fTanLorAng)
,fGlobalChargeScale(s.fGlobalChargeScale)
,fSpread2DHisto(s.fSpread2DHisto)
,fSpreadFun(s.fSpreadFun)
,fROTimeFun(s.fROTimeFun)
{
    //     Copy Constructor
}


//______________________________________________________________________
AliITSMFTSimulationPix::~AliITSMFTSimulationPix()
{
    // destructor
    // only the sens map is owned and it is deleted by ~AliITSMFTSimulation
}

//______________________________________________________________________
AliITSMFTSimulationPix& AliITSMFTSimulationPix::operator=(const AliITSMFTSimulationPix &s)
{
    //    Assignment operator
    if (&s == this) return *this;
    AliITSMFTSimulation::operator=(s);
    fSpread2DHisto = s.fSpread2DHisto;
    //
    fGlobalChargeScale = s.fGlobalChargeScale;
    fSpreadFun    = s.fSpreadFun;
    fROTimeFun    = s.fROTimeFun;
    //
    return *this;
}

//______________________________________________________________________
void AliITSMFTSimulationPix::Init()
{
    // Initilization
    if (fSimuParam->GetPixLorentzDrift()) SetTanLorAngle(fSimuParam->GetPixLorentzHoleWeight());
}

//______________________________________________________________________
Bool_t AliITSMFTSimulationPix::SetTanLorAngle(Double_t weightHole)
{
    // This function set the Tangent of the Lorentz angle.
    // A weighted average is used for electrons and holes
    // Input: Double_t weightHole: wheight for hole: it should be in the range [0,1]
    // output: Bool_t : kTRUE in case of success
    //
    if (weightHole<0) {
        weightHole=0.;
        AliWarning("You have asked for negative Hole weight");
        AliWarning("I'm going to use only electrons");
    }
    if (weightHole>1) {
        weightHole=1.;
        AliWarning("You have asked for weight > 1");
        AliWarning("I'm going to use only holes");
    }
    Double_t weightEle=1.-weightHole;
    AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
    if (!fld) AliFatal("The field is not initialized");
    Double_t bz = fld->SolenoidField();
    fTanLorAng = Tan(weightHole*fSimuParam->LorentzAngleHole(bz) +
                     weightEle*fSimuParam->LorentzAngleElectron(bz));
    return kTRUE;
}

//_____________________________________________________________________
void AliITSMFTSimulationPix::SDigitiseChip(TClonesArray *sdarray)
{
    //  This function begins the work of creating S-Digits.

    AliDebug(10,Form("In event %d chip %d there are %d hits", fEvent, fChip->GetIndex(),fChip->GetNHits()));
    if (fChip->GetNHits()) {
      if(fResponseParam->GetParameter(AliITSMFTSimuParam::kDigitalSim) == 0 ) Hits2SDigitsFast(); // analogue chip response
      else Hits2SDigitsFastDigital();                                         // digital chip response
    }
    if (!fSensMap->GetEntries()) return;
    WriteSDigits(sdarray);
    ClearMap();
    //
}

//______________________________________________________________________
void AliITSMFTSimulationPix::WriteSDigits(TClonesArray *sdarray)
{
    //  This function adds each S-Digit to pList
    int nsd = fSensMap->GetEntries();


    for (int i=0;i<nsd;i++) {
        AliITSMFTSDigit* sd = (AliITSMFTSDigit*)fSensMap->At(i); // ordered in index
        if (!(sd->GetSumSignal()>0) || fSensMap->IsDisabled(sd)) continue;
        new( (*sdarray)[sdarray->GetEntriesFast()]) AliITSMFTSDigit(*sd);
    }
    return;
}

//______________________________________________________________________
void AliITSMFTSimulationPix::FinishSDigitiseChip(TObjArray *detDigits)
{
    //  This function calls SDigitsToDigits which creates Digits from SDigits
    FrompListToDigits(detDigits);
    ClearMap();
    return;
}

//______________________________________________________________________
void AliITSMFTSimulationPix::DigitiseChip(TObjArray *detDigits)
{
    //  This function creates Digits straight from the hits and then adds
    //  electronic noise to the digits before adding them to pList
    //  Each of the input variables is passed along to Hits2SDigits
    //
    if(fResponseParam->GetParameter(AliITSMFTSimuParam::kDigitalSim) == 0 ) Hits2SDigitsFast(); // analogue chip response
    else  Hits2SDigitsFastDigital();                                         // digital chip response
    FinishSDigitiseChip(detDigits);
}

//______________________________________________________________________
void AliITSMFTSimulationPix::Hits2SDigits()
{
    // Does the charge distributions using Gaussian diffusion charge charing.
    Int_t nhits = fChip->GetNHits();
    if (!nhits) return;
    //
    Int_t h,ix,iz,i;
    Int_t idtrack;
    Float_t x,y,z;  // keep coordinates float (required by AliSegmentation)
    Double_t tof,x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0,de=0.0;
    Double_t t,tp,st,dt=0.2,el;
    Double_t thick = 0.5*fSeg->Dy();  // Half Thickness

    //
    for (h=0;h<nhits;h++) {
        //
        if (!fChip->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,tof,idtrack)) continue;
        st = Sqrt(x1*x1+y1*y1+z1*z1);
        if (st>0.0) {
            st = (Double_t)((Int_t)(st*1e4)); // number of microns
            if (st<=1.0) st = 1.0;
            dt = 1.0/st;               // RS TODO: do we need 1 micron steps?
            double dy = dt*thick;
            y = -0.5*dy;
            for (t=0.0;t<1.0;t+=dt) { // Integrate over t
                tp  = t+0.5*dt;
                x   = x0+x1*tp;
                y  += dy;
                z   = z0+z1*tp;
                if (!(fSeg->LocalToDet(x,z,ix,iz))) continue; // outside
                el  = fGlobalChargeScale * dt * de / fSimuParam->GetGeVToCharge();
                //
                if (fSimuParam->GetPixLorentzDrift()) x += y*fTanLorAng;
                SpreadCharge2D(x,z,y,ix,iz,el,tof,idtrack,h);
            } // end for t
        } else { // st == 0.0 deposit it at this point
            x   = x0;
            y   = y0 + 0.5*thick;
            z   = z0;
            if (!(fSeg->LocalToDet(x,z,ix,iz))) continue; // outside
            el  = fGlobalChargeScale * de / fSimuParam->GetGeVToCharge();
            if (fSimuParam->GetPixLorentzDrift()) x += y*fTanLorAng;
            SpreadCharge2D(x,z,y,ix,iz,el,tof,idtrack,h);
        } // end if st>0.0
    } // Loop over all hits h
    //
    // Coupling
    int nd = fSensMap->GetEntriesUnsorted(); // use unsorted access when possible, since it is faster
    AliITSMFTSDigit* dg = 0;
    switch (fSimuParam->GetPixCouplingOption()) {
        case AliITSMFTSimuParam::kNoCouplingPix :
            break;
        case AliITSMFTSimuParam::kNewCouplingPix :
            for (i=nd;i--;) {
                dg = (AliITSMFTSDigit*)fSensMap->AtUnsorted(i);
                if (fSensMap->IsDisabled(dg)) continue;
                SetCoupling(dg);
            }
            break;
        case AliITSMFTSimuParam::kOldCouplingPix:
            for (i=nd;i--;) {
                dg = (AliITSMFTSDigit*)fSensMap->AtUnsorted(i);
                if (fSensMap->IsDisabled(dg)) continue;
                SetCouplingOld(dg);
            }
            break;
        default:
            break;

    } // end switch
    if (GetDebug(2)) AliInfo(Form("Finished fCoupling=%d",fSimuParam->GetPixCouplingOption()));
}

//______________________________________________________________________
void AliITSMFTSimulationPix::Hits2SDigitsFast()
{
    // Does the charge distributions using Gaussian diffusion charge charing.    // Inputs:
    //    AliITSMFTChip *mod  Pointer to this chip
    //
    TObjArray *hits = fChip->GetHits();
    Int_t nhits = hits->GetEntriesFast();
    if (nhits<=0) return;
    //
    Int_t h,ix,iz,i;
    Int_t idtrack;
    Float_t x,y,z; // keep coordinates float (required by AliSegmentation)
    Double_t tof,x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0;
    Double_t step,st,el,de=0.0;
    Double_t minDim = Min(fSeg->Dpx(1),fSeg->Dpz(1)); // RStmp: smallest pitch
    Double_t thick = fSeg->Dy();
    //
    for (h=0;h<nhits;h++) {
        //
        if (!fChip->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,tof,idtrack)) continue;
        //
        st = Sqrt(x1*x1+y1*y1+z1*z1);
        if (st>0.0) {
            int np = int(1.5*st/minDim);  //RStmp: inject the points in such a way that there is ~1.5 point per cell
            np = TMath::Max(1.0*np,fResponseParam->GetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps));
            AliDebug(10,Form(" Number of charge injection steps is set to %d ",np));
            double dstep = 1./np;
            double dy = dstep*thick;
            y = -0.5*dy;
            step = -0.5*dstep;
            for (i=0;i<np;i++) {          //RStmp Integrate over t
                //      for (i=0;i<kn10;i++) { // Integrate over t
                step  += dstep;  // RStmp kti[i];
                x   = x0+x1*step;
                y  += dy;
                z   = z0+z1*step;
                if (!(fSeg->LocalToDet(x,z,ix,iz))) continue; // outside
                el  = fGlobalChargeScale * dstep * de/fSimuParam->GetGeVToCharge();
                if (fSimuParam->GetPixLorentzDrift()) x += y*fTanLorAng;
                SpreadCharge2D(x,z,y,ix,iz,el,tof,idtrack,h);
            } // end for i // End Integrate over t
        }
        else { // st == 0.0 deposit it at this point
            x   = x0;
            y   = y0+0.5*thick;
            z   = z0;
            if (!(fSeg->LocalToDet(x,z,ix,iz))) continue; // outside
            el  = fGlobalChargeScale * de / fSimuParam->GetGeVToCharge();
            if (fSimuParam->GetPixLorentzDrift()) x += y*fTanLorAng;
            SpreadCharge2D(x,z,y,ix,iz,el,tof,idtrack,h);
        } // end if st>0.0

    } // Loop over all hits h

    // Coupling
    int nd = fSensMap->GetEntriesUnsorted(); // use unsorted access when possible, since it is faster
    AliITSMFTSDigit* dg = 0;
    switch (fSimuParam->GetPixCouplingOption()) {
        case AliITSMFTSimuParam::kNoCouplingPix :
            break;
        case AliITSMFTSimuParam::kNewCouplingPix :
            for (i=nd;i--;) {
                dg = (AliITSMFTSDigit*)fSensMap->AtUnsorted(i);
                if (fSensMap->IsDisabled(dg)) continue;
                SetCoupling(dg);
            }
        case AliITSMFTSimuParam::kOldCouplingPix:
            for (i=nd;i--;) {
                dg = (AliITSMFTSDigit*)fSensMap->AtUnsorted(i);
                if (fSensMap->IsDisabled(dg)) continue;
                SetCouplingOld(dg);
            }
            break;
        default:
            break;
    } // end switch
    if (GetDebug(2)) AliInfo(Form("Finished fCoupling=%d",fSimuParam->GetPixCouplingOption()));
}

//______________________________________________________________________
void AliITSMFTSimulationPix::SpreadCharge2D(Double_t x0,Double_t z0, Double_t dy, Int_t ix0,Int_t iz0,
                                          Double_t el, Double_t tof, Int_t tID, Int_t hID)
{
    // Spreads the charge over neighboring cells. Assume charge is distributed
    // as charge(x,z) = (el/2*pi*sigx*sigz)*exp(-arg)
    // arg=((x-x0)*(x-x0)/2*sigx*sigx)+((z-z0*z-z0)/2*sigz*sigz)
    // Defined this way, the integral over all x and z is el.
    // Inputs:
    //    Double_t x0   x position of point where charge is liberated (local)
    //    Double_t z0   z position of point where charge is liberated (local)
    //    Double_t dy   distance from the entrance surface (diffusion sigma may depend on it)
    //    Int_t    ix0  row of cell corresponding to point x0
    //    Int_t    iz0  columb of cell corresponding to point z0
    //    Double_t el   number of electrons liberated in this step
    //    Double_t sigx Sigma difusion along x for this step (y0 dependent)
    //    Double_t sigz Sigma difusion along z for this step (z0 dependent)
    //    Int_t    tID  track number
    //    Int_t    hID  hit "hit" index number
    //
    Int_t ix,iz,ixs,ixe,izs,ize;
    Float_t x,z;   // keep coordinates float (required by AliSegmentation)
    Float_t xdioshift = 0 , zdioshift = 0 ;
    Double_t s,dtIn[kNDtSpread]; // data transfered to spread function for integral calculation
    //
    if (GetDebug(2)) AliInfo(Form("(x0=%e,z0=%e,dy=%e, ix0=%d,iz0=%d,el=%e,tID=%d,hID=%d)",x0,z0,dy,ix0,iz0,el,tID,hID));
    //
    Double_t &x1 = dtIn[kCellX1];
    Double_t &x2 = dtIn[kCellX2];
    Double_t &z1 = dtIn[kCellZ1];
    Double_t &z2 = dtIn[kCellZ2];
    //
    int nx = GetResponseParam()->GetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs);
    int nz = GetResponseParam()->GetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs);
    //
    dtIn[kCellYDepth]  = dy;
    ixs = Max(-nx+ix0,0);
    ixe = Min( nx+ix0,fSeg->Npx()-1);
    izs = Max(-nz+iz0,0);
    ize = Min( nz+iz0,fSeg->Npz()-1);
    for (ix=ixs;ix<=ixe;ix++)
        for (iz=izs;iz<=ize;iz++) {
            //
            // Check if the hit is inside readout window
            int cycleRO = (((AliITSMFTSimulationPix*)this)->*AliITSMFTSimulationPix::fROTimeFun)(ix,iz,tof);
            if (Abs(cycleRO)>kMaxROCycleAccept) continue;
            //
            fSeg->DetToLocal(ix,iz,x,z); // pixel center
            xdioshift = zdioshift = 0;
            double dxi = fSeg->Dpx(ix);
            double dzi = fSeg->Dpz(iz);
            CalcDiodeShiftInPixel(ix,iz,xdioshift,zdioshift);    // Check and apply diode shift if needed
            xdioshift *= dxi;
            zdioshift *= dzi;
            dxi *= 0.5;
            dzi *= 0.5;
            //      printf("DShift: %d %d -> %.4f %.4f\n",ix,iz,xdioshift,zdioshift);
            x1  = (x + xdioshift) - x0;   // calculate distance of cell boundaries from injection center
            z1  = (z + zdioshift) - z0;
            x2  = x1 + dxi; // Upper
            x1 -= dxi;      // Lower
            z2  = z1 + dzi; // Upper
            z1 -= dzi;      // Lower
            s   = el* (((AliITSMFTSimulationPix*)this)->*AliITSMFTSimulationPix::fSpreadFun)(dtIn);
            if (s>fSimuParam->GetPixMinElToAdd()) UpdateMapSignal(iz,ix,tID,hID,s,cycleRO);
        } // end for ix, iz
    //
}

//______________________________________________________________________
Double_t AliITSMFTSimulationPix::SpreadFunDoubleGauss2D(const Double_t *dtIn)
{
    // calculate integral of charge in the cell with boundaries at X=dtIn[kCellX1]:dtIn[kCellX2]
    // and Z=dtIn[kCellZ1]:dtIn[kCellZ2]
    // The spread function is assumed to be double gaussian in 2D
    // Parameters should be: mean0,sigma0, mean1,sigma1, relative strenght of 2nd gaussian wrt 1st one
    //
    // 1st gaussian
    double intg1 = GausInt2D(fResponseParam->GetParameter(kG2SigX0),  // sigX
                             dtIn[kCellX1]-fResponseParam->GetParameter(kG2MeanX0),      // x1-xmean
                             dtIn[kCellX2]-fResponseParam->GetParameter(kG2MeanX0),      // x2-xmean
                             fResponseParam->GetParameter(kG2SigZ0),  // sigZ
                             dtIn[kCellZ1]-fResponseParam->GetParameter(kG2MeanZ0),    // z1-zmean
                             dtIn[kCellZ2]-fResponseParam->GetParameter(kG2MeanZ0));   // z2-zmean
    // 2nd gaussian
    double intg2 = GausInt2D(fResponseParam->GetParameter(kG2SigX1),  // sigX
                             dtIn[kCellX1]-fResponseParam->GetParameter(kG2MeanX1),      // x1-xmean
                             dtIn[kCellX2]-fResponseParam->GetParameter(kG2MeanX1),      // x2-xmean
                             fResponseParam->GetParameter(kG2SigZ1),  // sigZ
                             dtIn[kCellZ1]-fResponseParam->GetParameter(kG2MeanZ1),    // z1-zmean
                             dtIn[kCellZ2]-fResponseParam->GetParameter(kG2MeanZ1));   // z2-zmean
    double scl = fResponseParam->GetParameter(kG2ScaleG2);
    return (intg1+intg2*scl)/(1+scl);
    //
}


//______________________________________________________________________
Double_t AliITSMFTSimulationPix::SpreadFrom2DHisto(const Double_t *dtIn)
{
    // calculate integral of charge in the cell with boundaries at X=dtIn[kCellX1]:dtIn[kCellX2]
    // and Z=dtIn[kCellZ1]:dtIn[kCellZ2]
    // The spread function integral is taken from fSpread2DHisto extracted from the sensor response parameters
    // list in the method SetResponseParam. The histo must return the fraction of charge integrates in the
    // cell whose center is on the distance X=(dtIn[kCellX1]+dtIn[kCellX2])/2 and Z=(dtIn[kCellZ1]+dtIn[kCellZ2])/2
    // from the injection point.
    //
    Double_t qpixfrac = 0;
    Double_t xintp = 1e4*(dtIn[kCellX1]+dtIn[kCellX2])/2.0;
    Double_t zintp = 1e4*(dtIn[kCellZ1]+dtIn[kCellZ2])/2.0;
    //
    qpixfrac =  fSpread2DHisto->Interpolate(xintp,zintp); //the PSF map is given in um but the dtIn is in cm so we need to convert it
    //
    return qpixfrac;
}

//______________________________________________________________________
Double_t AliITSMFTSimulationPix::SpreadFunGauss2D(const Double_t *dtIn)
{
    // calculate integral of charge in the cell with boundaries at X=dtIn[kCellX1]:dtIn[kCellX2]
    // and Z=dtIn[kCellZ1]:dtIn[kCellZ2]
    // The spread function is assumed to be gaussian in 2D
    // Parameters should be: mean0,sigma0
    return GausInt2D(fResponseParam->GetParameter(kG1SigX),  // sigX
                     dtIn[kCellX1]-fResponseParam->GetParameter(kG1MeanX),    // x1-xmean
                     dtIn[kCellX2]-fResponseParam->GetParameter(kG1MeanX),    // x2-xmean
                     fResponseParam->GetParameter(kG1SigZ),  // sigZ
                     dtIn[kCellZ1]-fResponseParam->GetParameter(kG1MeanZ),    // z1-zmean
                     dtIn[kCellZ2]-fResponseParam->GetParameter(kG1MeanZ));   // z2-zmean
    //
}

//______________________________________________________________________
void AliITSMFTSimulationPix::RemoveDeadPixels()
{
    // Removes dead pixels on each chip (ladder)
    // This should be called before going from sdigits to digits (i.e. from FrompListToDigits)

    AliITSMFTCalibrationPix* calObj = (AliITSMFTCalibrationPix*) GetCalibDead();
    if (!calObj) return;
    //
    if (calObj->IsBad()) {ClearMap(); return;} // whole chip is masked
    //
    // prepare the list of r/o cycles seen
    Char_t cyclesSeen[2*kMaxROCycleAccept+1];
    int ncycles = 0;
    for (int i=(2*kMaxROCycleAccept+1);i--;) if (fCyclesID[i]) cyclesSeen[ncycles++]=i-kMaxROCycleAccept;

    // remove single bad pixels one by one
    int nsingle = calObj->GetNrBadSingle();
    UInt_t col,row;
    Int_t cycle;
    for (int i=nsingle;i--;) {
        calObj->GetBadPixelSingle(i,row,col);
        for (int icl=ncycles;icl--;) fSensMap->DeleteItem(col,row,cyclesSeen[icl]);
    }
    int nsd = fSensMap->GetEntriesUnsorted();
    for (int isd=nsd;isd--;) {
        AliITSMFTSDigit* sd = (AliITSMFTSDigit*)fSensMap->AtUnsorted(isd);
        if (fSensMap->IsDisabled(sd)) continue;
        fSensMap->GetMapIndex(sd->GetUniqueID(),col,row,cycle);
        int chip = fSeg->GetChipFromChannel(0,col);
        //    if (calObj->IsChipMarkedBad(chip)) fSensMap->Disable(sd); // this will simple mark the hit as bad
        if (calObj->IsChipMarkedBad(chip)) fSensMap->DeleteItem(sd); // this will suppress hit in the sorted list
    }
    //
}

//______________________________________________________________________
void AliITSMFTSimulationPix::AddNoisyPixels()
{
    // Adds noisy pixels on each chip (ladder)
    // This should be called before going from sdigits to digits (i.e. FrompListToDigits)
    AliITSMFTCalibrationPix* calObj = (AliITSMFTCalibrationPix*) GetCalibNoisy();
    if (!calObj) { AliDebug(10,Form("  No Calib Object for Noise!!! ")); return;}
    for (Int_t i=calObj->GetNrBad(); i--;)
    {
        if ( fResponseParam->GetParameter(AliITSMFTSimuParam::kDigitalSim) < 1.0 )
            UpdateMapNoise(calObj->GetBadColAt(i), calObj->GetBadRowAt(i),10*fSimuParam->GetPixThreshold(fChip->GetIndex()));
        else
            UpdateMapNoise(calObj->GetBadColAt(i), calObj->GetBadRowAt(i),kNoisyPixOCDB );
    }
    //
}

//______________________________________________________________________
void AliITSMFTSimulationPix::FrompListToDigits(TObjArray *detDigits)
{
    // add noise and electronics, perform the zero suppression and add the digits to the list
    //
    // RS may use for ordered random sample generation dl.acm.org/ft_gateway.cfm?id=356313&type=pdf
    //
    int nsd = fSensMap->GetEntriesUnsorted();  // sdigits added from the signal
    //
    // add different kinds of noise.
    Bool_t addNoisy = fSimuParam->GetPixAddNoisyFlag() && (nsd>0 || fSimuParam->GetPixNoiseInAllMod()); // do we generate noise?
    if (addNoisy) {
        AddNoisyPixels();       // constantly noisy channels
        AddRandomNoisePixels(0.0); // random noise: at the moment generate noise only for instance 0
        nsd = fSensMap->GetEntriesUnsorted();
    }
    //
    if (nsd && fSimuParam->GetPixRemoveDeadFlag()) {
        RemoveDeadPixels();
        // note that here we shall use GetEntries instead of GetEntriesUnsorted since the
        // later operates on the array where the elements are not removed by flagged
        nsd = fSensMap->GetEntries();
    }
    if (!nsd) return; // nothing to digitize
    //
    UInt_t row,col;
    Int_t iCycle,modId = fChip->GetIndex();
    Double_t sig;
    const Int_t    knmaxtrk=AliITSMFTDigitPix::GetNTracks();
    static AliITSMFTDigitPix dig;
    //
    for (int i=0;i<nsd;i++) {
        AliITSMFTSDigit* sd = (AliITSMFTSDigit*)fSensMap->At(i); // ordered in index
        if (fSensMap->IsDisabled(sd)) continue;
        //
	sig=sd->GetSumSignal();
        if ( fResponseParam->GetParameter(AliITSMFTSimuParam::kDigitalSim) < 1.0 &&
	     sig<=fSimuParam->GetPixThreshold(modId)) continue;   //Threshold only applies in analogue simulation
        //
        if (Abs(sig)>2147483647.0) { //RS?
            //PH 2147483647 is the max. integer
            //PH This apparently is a problem which needs investigation
            AliWarning(Form("Too big or too small signal value %f",sig));
        }
        fSensMap->GetMapIndex(sd->GetUniqueID(),col,row,iCycle);
        dig.SetCoord1(col);
        dig.SetCoord2(row);
        dig.SetROCycle(iCycle);
        dig.SetSignal((Int_t)sig);
        dig.SetSignalPix((Int_t)sig);
        int ntr = sd->GetNTracks();
        for (int j=0;j<ntr;j++) {
            dig.SetTrack(j,sd->GetTrack(j));
            dig.SetHit(j,sd->GetHit(j));
        }
        for (int j=ntr;j<knmaxtrk;j++) {
            dig.SetTrack(j,-3);
            dig.SetHit(j,-1);
        }

        Int_t branch=AliITSMFTAux::kChipTypePix;
        TClonesArray &ldigits = *((TClonesArray*)detDigits->At(branch));
        int nd = ldigits.GetEntriesFast();
        switch(branch){
        case AliITSMFTAux::kChipTypePix:
           new(ldigits[nd]) AliITSMFTDigitPix(dig);
           break;
        default:
           AliFatal(Form("Unknown digits branch %d",branch));
        }

    }
    //
}

//______________________________________________________________________
Int_t AliITSMFTSimulationPix::AddRandomNoisePixels(Double_t tof)
{
  // create random noisy sdigits above threshold
  //
  int modId = fChip->GetIndex();
  int npix = fSeg->GetNPads();
  int ncand = gRandom->Poisson( npix*fSimuParam->GetPixFakeRate() );
  if (ncand<1) return 0;
  //
  double probNoisy,noiseSig,noiseMean,thresh;
  //
  UInt_t row,col;
  Int_t iCycle;
  static TArrayI ordSampleInd(100),ordSample(100); //RS!!! static is not thread-safe!!!
  ncand = GenOrderedSample(npix,ncand,ordSample,ordSampleInd);
  int* ordV = ordSample.GetArray();
  int* ordI = ordSampleInd.GetArray();
  //
  if ( fResponseParam->GetParameter(AliITSMFTSimuParam::kDigitalSim) < 1.0 ) {
    thresh = fSimuParam->GetPixThreshold(modId);
    fSimuParam->GetPixNoise(modId, noiseSig, noiseMean);
      probNoisy = AliITSMFTSimuParam::CalcProbNoiseOverThreshold(noiseMean,noiseSig,thresh); // prob. to have noise above threshold
    //
    for (int j=0;j<ncand;j++) {
      fSensMap->GetMapIndex((UInt_t)ordV[ordI[j]],col,row,iCycle);   // create noisy digit
      iCycle = (((AliITSMFTSimulationPix*)this)->*AliITSMFTSimulationPix::fROTimeFun)(row,col,tof);
      UpdateMapNoise(col,row,AliITSMFTSimuParam::GenerateNoiseQFunction(probNoisy,noiseMean,noiseSig),  iCycle);
    }
  }
  else {
    for (int j=0;j<ncand;j++) {
      fSensMap->GetMapIndex((UInt_t)ordV[ordI[j]],col,row,iCycle);   // create noisy digit
      iCycle = (((AliITSMFTSimulationPix*)this)->*AliITSMFTSimulationPix::fROTimeFun)(row,col,tof);
      UpdateMapNoise(col,row,kNoisyPixRnd, iCycle);
    }
  }
  return ncand;
}


//______________________________________________________________________
void AliITSMFTSimulationPix::SetCoupling(AliITSMFTSDigit* old)
{
    //  Take into account the coupling between adiacent pixels.
    //  The parameters probcol and probrow are the probability of the
    //  signal in one pixel shared in the two adjacent pixels along
    //  the column and row direction, respectively.
    //  Note pList is goten via GetMap() and chip is not need any more.
    //  Otherwise it is identical to that coded by Tiziano Virgili (BSN).
    UInt_t col,row;
    Int_t iCycle;
    Double_t pulse1,pulse2;
    Double_t couplR=0.0,couplC=0.0;
    Double_t xr=0.;
    //
    fSensMap->GetMapIndex(old->GetUniqueID(),col,row,iCycle);
    int cycle = iCycle;
    fSimuParam->GetPixCouplingParam(couplC,couplR);
    if (GetDebug(2)) AliInfo(Form("(col=%d,row=%d,couplC=%e couplR=%e",
                                  col,row,couplC,couplR));
    pulse2 = pulse1 = old->GetSignal();
    if (pulse1<fSimuParam->GetPixMinElToAdd()) return; // too small signal
    for (Int_t isign=-1;isign<=1;isign+=2) {
        //
        // loop in col direction
        int j1 = int(col) + isign;
        xr = gRandom->Rndm();
        if ( !((j1<0) || (j1>fSeg->Npz()-1) || (xr>couplC)) ) UpdateMapSignal(UInt_t(j1),row,old->GetTrack(0),old->GetHit(0),pulse1,cycle);
        //
        // loop in row direction
        int j2 = int(row) + isign;
        xr = gRandom->Rndm();
        if ( !((j2<0) || (j2>fSeg->Npx()-1) || (xr>couplR)) ) UpdateMapSignal(col,UInt_t(j2),old->GetTrack(0),old->GetHit(0),pulse2,cycle);
    }
    //
}

//______________________________________________________________________
void AliITSMFTSimulationPix::SetCouplingOld(AliITSMFTSDigit* old)
{
    //  Take into account the coupling between adiacent pixels.
    //  The parameters probcol and probrow are the fractions of the
    //  signal in one pixel shared in the two adjacent pixels along
    //  the column and row direction, respectively.
    // Inputs:
    // old            existing AliITSMFTSDigit
    // ntrack         track incex number
    // idhit          hit index number
    // chip         chip number
    //
    UInt_t col,row;
    Int_t cycle;
    Int_t modId = fChip->GetIndex();
    Double_t pulse1,pulse2;
    Double_t couplR=0.0,couplC=0.0;
    //
    fSensMap->GetMapIndex(old->GetUniqueID(),col,row,cycle);
    fSimuParam->GetPixCouplingParam(couplC,couplR);
    if (GetDebug(3)) AliInfo(Form("(col=%d,row=%d,roCycle=%d)  couplC=%e couplR=%e",col,row,cycle,couplC,couplR));
    //
    if (old->GetSignal()<fSimuParam->GetPixMinElToAdd()) return; // too small signal
    for (Int_t isign=-1;isign<=1;isign+=2) {// loop in col direction
        pulse2 = pulse1 = old->GetSignal();
        //
        int j1 = int(col)+isign;
        pulse1 *= couplC;
        if ((j1<0)||(j1>fSeg->Npz()-1)||(pulse1<fSimuParam->GetPixThreshold(modId))) pulse1 = old->GetSignal();
        else UpdateMapSignal(UInt_t(j1),row,old->GetTrack(0),old->GetHit(0),pulse1,cycle);

        // loop in row direction
        int j2 = int(row) + isign;
        pulse2 *= couplR;
        if ((j2<0)||(j2>(fSeg->Npx()-1))||(pulse2<fSimuParam->GetPixThreshold(modId))) pulse2 = old->GetSignal();
        else UpdateMapSignal(col,UInt_t(j2),old->GetTrack(0),old->GetHit(0),pulse2,cycle);
    } // for isign
}

//______________________________________________________________________
void AliITSMFTSimulationPix::SetResponseParam(AliITSMFTParamList* resp)
{
    // attach response parameterisation data
    fResponseParam = resp;
    //
    int spreadID = Nint(fResponseParam->GetParameter(AliITSMFTSimuParam::kChargeSpreadType));
    const char* hname = 0;
    fSpread2DHisto = 0;
    //
    switch (spreadID) {
            //
        case kSpreadFunHisto:
            fSpreadFun = &AliITSMFTSimulationPix::SpreadFrom2DHisto;
            hname = fResponseParam->GetParName(AliITSMFTSimuParam::kChargeSpreadType);
            if (!(fSpread2DHisto=(TH2*)fResponseParam->GetParamObject(hname)))
                AliFatal(Form("Did not find 2D histo %s for charge spread parameterization",hname));
            break;
            //
        case kSpreadFunDoubleGauss2D:
            fSpreadFun = &AliITSMFTSimulationPix::SpreadFunDoubleGauss2D;
            break;
            //
        case kSpreadFunGauss2D:
            fSpreadFun = &AliITSMFTSimulationPix::SpreadFunGauss2D;
            break;
            //
        default: AliFatal(Form("Did not find requested spread function id=%d",spreadID));
    }
    //
    int readoutType = Nint(fResponseParam->GetParameter(AliITSMFTSimuParam::kReadOutSchemeType));
    switch (readoutType) {
        case AliITSMFTSimuParam::kReadOutStrobe:
            fROTimeFun = &AliITSMFTSimulationPix::GetReadOutCycle;
            break;
        case AliITSMFTSimuParam::kReadOutRollingShutter:
            fROTimeFun = &AliITSMFTSimulationPix::GetReadOutCycleRollingShutter;
            break;
        default: AliFatal(Form("Did not find requested readout time type id=%d",readoutType));
    }

    //___ Set the Rolling Shutter read-out window
    fReadOutCycleLength = fResponseParam->GetParameter(AliITSMFTSimuParam::kReadOutCycleLength);
    //___ Pixel discrimination threshold, and the S/N cut
    fSimuParam->SetPixThreshold(fResponseParam->GetParameter(AliITSMFTSimuParam::kPixNoiseMPV) *fResponseParam->GetParameter(AliITSMFTSimuParam::kPixSNDisrcCut) , fResponseParam->GetParameter(AliITSMFTSimuParam::kPixSNDisrcCut),-1); //for all chips
    //___ Minimum number of electrons to add
    fSimuParam->SetPixMinElToAdd(fResponseParam->GetParameter(AliITSMFTSimuParam::kPixMinElToAdd));
    //___ Set the Pixel Noise MPV and Sigma (the noise distribution is Landau not Gauss due to RTN)
    fSimuParam->SetPixNoise( fResponseParam->GetParameter(AliITSMFTSimuParam::kPixNoiseMPV), fResponseParam->GetParameter(AliITSMFTSimuParam::kPixNoiseSigma), -1); //for all chips
    //___ Pixel fake hit rate
    fSimuParam->SetPixFakeRate( fResponseParam->GetParameter(AliITSMFTSimuParam::kPixFakeRate) );
    //___ To apply the noise or not
    if (  fResponseParam->GetParameter(AliITSMFTSimuParam::kPixNoiseIsOn) > 0.01)  fSimuParam->SetPixAddNoisyFlag(kTRUE);
    else fSimuParam->SetPixAddNoisyFlag(kFALSE);
    //
    if(fResponseParam->GetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod) > 0.01 ) fSimuParam->SetPixNoiseInAllMod(kTRUE);
    else fSimuParam->SetPixNoiseInAllMod(kFALSE);
    //
    //  Double_t vGeVToQ = fSimuParam->GetGeVToCharge();
    fGlobalChargeScale = fResponseParam->GetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale);

    AliDebug(10,Form("=============== Setting the response start ============================"));
    AliDebug(10,Form("=============== Digital (1) / Analogue (0) simu: %f",fResponseParam->GetParameter(AliITSMFTSimuParam::kDigitalSim)));
    AliDebug(10,Form("=============== RO type: %d",readoutType));
    AliDebug(10,Form("=============== RO cycle lenght: %lf",fReadOutCycleLength));
    AliDebug(10,Form("=============== Noise MPV: %lf",fResponseParam->GetParameter(AliITSMFTSimuParam::kPixNoiseMPV)));
    AliDebug(10,Form("=============== Noise Sigma: %lf",fResponseParam->GetParameter(AliITSMFTSimuParam::kPixNoiseSigma)));
    AliDebug(10,Form("=============== Fake rate: %lf",fResponseParam->GetParameter(AliITSMFTSimuParam::kPixFakeRate)));
    AliDebug(10,Form("=============== Noise On/Off: %d",fSimuParam->GetPixAddNoisyFlag()));
    AliDebug(10,Form("=============== Noise in all mod on/off: %d",fSimuParam->GetPixNoiseInAllMod()));
    AliDebug(10,Form("=============== Global Charge scale: %lf",fGlobalChargeScale));
    AliDebug(10,Form("=============== Setting the response done  ============================"));

}

//______________________________________________________________________
Int_t AliITSMFTSimulationPix::GetReadOutCycleRollingShutter(Int_t row, Int_t col, Double_t hitTime)
{
    //
    // Get the read-out cycle of the hit in the given column/row of the sensor.
    // hitTime is the time of the subhit (hit is divided to nstep charge deposit) in seconds
    // globalPhaseShift gives the start of the RO for the cycle in pixel wrt the LHC clock
    // GetRollingShutterWindow give the with of the rolling shutter read out window
    //
    double tmin = fReadOutCycleOffset + fReadOutCycleLength*(double(row)/fSeg->Npx()-1.);
    int cycle = Nint( (hitTime-tmin)/fReadOutCycleLength - 0.5 );
    AliDebug(3,Form("Rolling shutter at row%d/col%d: particle time: %e, tmin: %e : tmax %e -> cycle:%d",row,col,hitTime,tmin,
                    tmin+fReadOutCycleLength,cycle));
    return cycle;
    //
}

//______________________________________________________________________
Int_t AliITSMFTSimulationPix::GetReadOutCycle(Int_t row, Int_t col, Double_t hitTime)
{
    //
    // Check whether the hit is in the read out window of the given column/row of the sensor
    //
    AliDebug(3,Form("Strobe readout: row%d/col%d: particle time: %e, tmin: %e, tmax %e",
                    row,col,hitTime,fReadOutCycleOffset,fReadOutCycleOffset+fReadOutCycleLength));
    hitTime -= fReadOutCycleOffset-0.5*fReadOutCycleLength;
    return (hitTime<0 || hitTime>fReadOutCycleLength) ? kMaxROCycleAccept+1 : 0;
    //
}

//_______________________________________________________________________
void AliITSMFTSimulationPix::CalcDiodeShiftInPixel(Int_t xrow, Int_t zcol, Float_t &x, Float_t &z)
{
    //
    // Calculates the shift of the diode wrt the geometric center of the pixel.
    // It is needed for staggerred pixel layout or double diode pixels with assymetric center
    // The shift can depend on the column or line or both...
    // The x and z are passed in cm
    //
    ((AliITSMFTSegmentationPix*)fSeg)->GetDiodShift(xrow,zcol,x,z);
    //
}

//______________________________________________________________________
void AliITSMFTSimulationPix::Hits2SDigitsFastDigital()
{
    // Does the digital chip response simulation
    //    AliITSMFTChip *mod  Pointer to this chip
    //


    TObjArray *hits = fChip->GetHits();
    Int_t nhits = hits->GetEntriesFast();
    if (nhits<=0) return;
    //
    Int_t h,ix,iz;
    Int_t idtrack;
    Double_t tof,x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0;
    Double_t el,de=0.0;

    //
    for (h=0;h<nhits;h++) {
        //
        if (!fChip->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,tof,idtrack)) continue;
        //
	//double st = Sqrt(x1*x1+y1*y1+z1*z1);

        //___ place hit to the middle of the path segment - CHANGE LATER !
	// keep coordinates float (required by AliSegmentation)
        float x   = (x0+x1)/2.0;
        //float y   = (y0+y1)/2.0;
        float z   = (z0+z1)/2.0;
        //
        if (!(fSeg->LocalToDet(x,z,ix,iz))) continue; // outside
        //
        // Check if the hit is inside readout window
        int cycleRO = (((AliITSMFTSimulationPix*)this)->*AliITSMFTSimulationPix::fROTimeFun)(ix,iz,tof);
        if (Abs(cycleRO)>kMaxROCycleAccept) continue;
        //
        //
        el  = de / fSimuParam->GetGeVToCharge();
        //
        PlaceDigitalPixels(x,z,el,tof,idtrack,h);
        //
    } // Loop over all hits h

}
//______________________________________________________________________
void AliITSMFTSimulationPix::PlaceDigitalPixels(Double_t x0hit,Double_t z0hit, Double_t el, Double_t tof, Int_t tID, Int_t hID)
{
    // Place the digital pixel positions on the sensor
    // Inputs:
    //    Double_t x0hit   x position of point where charge is liberated (local) - hit
    //    Double_t z0hit   z position of point where charge is liberated (local) - hit
    //    Double_t el   number of electrons liberated in this step
    //    Double_t sigx Sigma difusion along x for this step (y0 dependent)
    //    Double_t sigz Sigma difusion along z for this step (z0 dependent)
    //    Int_t    tID  track number
    //    Int_t    hID  hit "hit" index number
    //


    Int_t ix,iz,nx,nz;
    Float_t x,z;   // keep coordinates float (required by AliSegmentation)
    Float_t distX = 0, distZ = 0;

    //___ TEMPORARY - place a fixed pattern cluster COG to a distance d(x,z) away from hit - TEMPORARY
    // Pattern used (not realistic ebye but averaging over events yes):
    //  -+-
    //  +++
    //  -+-
    //
    //___ This list should come from a look up table based on CluTypeID as well as COG coord

    //
    TRandom3 rnd;
    distX = rnd.Uniform(-5.0*1e-4, 5.0*1e-4); //in um
    distZ = rnd.Uniform(-5.0*1e-4, 5.0*1e-4); //in um
    //
    x = x0hit + distX;
    z = z0hit + distZ;
    //
    if(!fSeg->LocalToDet(x,z,ix,iz)) return; // if clu CoG is outside of the chip skipp the cluster -> refine later
    //
    const Int_t nCluPixels = 5;
    Int_t aPixListX[nCluPixels] = { 0, -1, 0, 1,  0};
    Int_t aPixListZ[nCluPixels] = { 1,  0, 0, 0, -1};
    //
    Double_t s = el / 1.0 / nCluPixels;
    //
    int cycleRO;
    //
    for (Int_t ipix = 0 ; ipix < nCluPixels; ipix++)
    {
        nx = ix + aPixListX[ipix];
        nz = iz + aPixListZ[ipix];
        cycleRO = (((AliITSMFTSimulationPix*)this)->*AliITSMFTSimulationPix::fROTimeFun)(ix,iz,tof);
        if ( nx >= 0 && nx <= fSeg -> Npx() && nz >= 0 && nz <= fSeg -> Npz() ) UpdateMapSignal(nz,nx,tID,hID,s,cycleRO); //if the pixel is in the detector
    }

}
