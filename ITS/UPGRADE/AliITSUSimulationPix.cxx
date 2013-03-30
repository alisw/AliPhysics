/*
 questions to experts: why RemoveDeadPixels should be called before FrompListToDigits ? 

 
*/

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
#include <TString.h>
#include "AliITSU.h"
#include "AliITSUDigitPix.h"
#include "AliITSUHit.h"
#include "AliITSUModule.h"
#include "AliITSUSensMap.h"
#include "AliITSUCalibrationPix.h"
#include "AliITSUSegmentationPix.h"
#include "AliITSUSimulationPix.h"
#include "AliLog.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliMathBase.h"
#include "AliITSUSimuParam.h"
#include "AliITSUSDigit.h"
#include "AliITSUParamList.h"

using std::cout;
using std::endl;
using namespace TMath;

ClassImp(AliITSUSimulationPix)
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
//  AliITSUSimulationPix is to do the simulation of pixels
//
//  2013 Feb: Added MonoPix response and nois calculation al la MIMOSA32 (levente.molnar@cern.ch)
//
////////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSUSimulationPix::AliITSUSimulationPix()
:  fTanLorAng(0)
  ,fReadOutCycleLength(25e-6)
  ,fReadOutCycleOffset(0)
  ,fGlobalChargeScale(1.0)
  ,fSpread2DHisto(0)
  ,fSpreadFun(0)
  ,fROTimeFun(0)
{
   // Default constructor.
  SetUniqueID(AliITSUGeomTGeo::kDetTypePix);
}

//______________________________________________________________________
AliITSUSimulationPix::AliITSUSimulationPix(AliITSUSimuParam* sim,AliITSUSensMap* map)
  :AliITSUSimulation(sim,map)
  ,fTanLorAng(0)
  ,fReadOutCycleLength(25e-6)
  ,fReadOutCycleOffset(0)
  ,fGlobalChargeScale(1.0)
  ,fSpread2DHisto(0)
  ,fSpreadFun(0)
  ,fROTimeFun(0)
{
  // standard constructor
  SetUniqueID(AliITSUGeomTGeo::kDetTypePix);
  Init();
}

//______________________________________________________________________
AliITSUSimulationPix::AliITSUSimulationPix(const AliITSUSimulationPix &s) 
  :AliITSUSimulation(s)
  ,fTanLorAng(s.fTanLorAng)
  ,fReadOutCycleLength(s.fReadOutCycleLength)
  ,fReadOutCycleOffset(s.fReadOutCycleOffset)
  ,fGlobalChargeScale(s.fGlobalChargeScale)
  ,fSpread2DHisto(s.fSpread2DHisto)
  ,fSpreadFun(s.fSpreadFun)
  ,fROTimeFun(s.fROTimeFun)
{
  //     Copy Constructor
}


//______________________________________________________________________
AliITSUSimulationPix::~AliITSUSimulationPix()
{
  // destructor
  // only the sens map is owned and it is deleted by ~AliITSUSimulation
}

//______________________________________________________________________
AliITSUSimulationPix& AliITSUSimulationPix::operator=(const AliITSUSimulationPix &s)
{
  //    Assignment operator
  if (&s == this) return *this;
  AliITSUSimulation::operator=(s);
  fReadOutCycleLength = s.fReadOutCycleLength;
  fReadOutCycleOffset = s.fReadOutCycleOffset;
  fSpread2DHisto = s.fSpread2DHisto;
  fGlobalChargeScale = s.fGlobalChargeScale;
  fSpreadFun    = s.fSpreadFun;
  fROTimeFun    = s.fROTimeFun;
  //
  return *this;
}

//______________________________________________________________________
void AliITSUSimulationPix::Init()
{
  // Initilization
  if (fSimuParam->GetPixLorentzDrift()) SetTanLorAngle(fSimuParam->GetPixLorentzHoleWeight());
}

//______________________________________________________________________
Bool_t AliITSUSimulationPix::SetTanLorAngle(Double_t weightHole) 
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
void AliITSUSimulationPix::SDigitiseModule()
{
  //  This function begins the work of creating S-Digits.
    
  AliDebug(10,Form("In event %d module %d there are %d hits", fEvent, fModule->GetIndex(),fModule->GetNHits()));       
  //      
  if (fModule->GetNHits()) Hits2SDigitsFast();
  //
  if (fSimuParam->GetPixAddNoisyFlag())   AddNoisyPixels();
  if (fSimuParam->GetPixRemoveDeadFlag()) RemoveDeadPixels(); 
  if (!fSensMap->GetEntries()) return;
  WriteSDigits();
  ClearMap();
  //
}

//______________________________________________________________________
void AliITSUSimulationPix::WriteSDigits()
{
  //  This function adds each S-Digit to pList
  static AliITSU *aliITS = (AliITSU*)gAlice->GetModule("ITS");
  int nsd = fSensMap->GetEntries();
  for (int i=0;i<nsd;i++) {
    AliITSUSDigit* sd = (AliITSUSDigit*)fSensMap->At(i); // ordered in index
    if (!sd->GetSumSignal()>0 || fSensMap->IsDisabled(sd)) continue;
    aliITS->AddSumDigit(*sd);
  }
  return; 
}

//______________________________________________________________________
void AliITSUSimulationPix::FinishSDigitiseModule()
{
   //  This function calls SDigitsToDigits which creates Digits from SDigits
  FrompListToDigits();
  ClearMap();
  return;
}

//______________________________________________________________________
void AliITSUSimulationPix::DigitiseModule()
{
  //  This function creates Digits straight from the hits and then adds
  //  electronic noise to the digits before adding them to pList
  //  Each of the input variables is passed along to Hits2SDigits
  //
  // pick charge spread function
  Hits2SDigitsFast();
  //
  if (fSimuParam->GetPixAddNoisyFlag())   AddNoisyPixels();
  if (fSimuParam->GetPixRemoveDeadFlag()) RemoveDeadPixels();
  FrompListToDigits();
  ClearMap();
}

//______________________________________________________________________
void AliITSUSimulationPix::Hits2SDigits()
{
  // Does the charge distributions using Gaussian diffusion charge charing.
  Int_t nhits = fModule->GetNHits();
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
    if (!fModule->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,tof,idtrack)) continue;
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
	el  = dt * de / fSimuParam->GetGeVToCharge();
	//
	if (fSimuParam->GetPixLorentzDrift()) x += y*fTanLorAng;
	SpreadCharge2D(x,z,y,ix,iz,el,tof,idtrack,h);
      } // end for t
    } else { // st == 0.0 deposit it at this point
      x   = x0;
      y   = y0 + 0.5*thick;
      z   = z0;
      if (!(fSeg->LocalToDet(x,z,ix,iz))) continue; // outside
      el  = de / fSimuParam->GetGeVToCharge();
      if (fSimuParam->GetPixLorentzDrift()) x += y*fTanLorAng;
      SpreadCharge2D(x,z,y,ix,iz,el,tof,idtrack,h);
    } // end if st>0.0    
  } // Loop over all hits h
  //
  // Coupling
  int nd = fSensMap->GetEntriesUnsorted(); // use unsorted access when possible, since it is faster
  AliITSUSDigit* dg = 0;
  switch (fSimuParam->GetPixCouplingOption()) {
  case AliITSUSimuParam::kNewCouplingPix :
    for (i=nd;i--;) {
      dg = (AliITSUSDigit*)fSensMap->AtUnsorted(i);
      if (fSensMap->IsDisabled(dg)) continue;
      SetCoupling(dg,idtrack,h);
    } 
    break;
  case AliITSUSimuParam::kOldCouplingPix:
    for (i=nd;i--;) {
      dg = (AliITSUSDigit*)fSensMap->AtUnsorted(i);
      if (fSensMap->IsDisabled(dg)) continue;
      SetCouplingOld(dg,idtrack,h);
    } 
    break;
  default:
    break;
    
  } // end switch
  if (GetDebug(2)) AliInfo(Form("Finished fCoupling=%d",fSimuParam->GetPixCouplingOption()));
}

//______________________________________________________________________
void AliITSUSimulationPix::Hits2SDigitsFast()
{
  // Does the charge distributions using Gaussian diffusion charge charing.    // Inputs:
  //    AliITSUModule *mod  Pointer to this module
  //
  TObjArray *hits = fModule->GetHits();
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
    if (!fModule->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,tof,idtrack)) continue;
    //
    st = Sqrt(x1*x1+y1*y1+z1*z1); 
    if (st>0.0) {
      int np = int(1.5*st/minDim);  //RStmp: inject the points in such a way that there is ~1.5 point per cell
      np = TMath::Max(1.0*np,fResponseParam->GetParameter(kSpreadFunMinSteps));   
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
      el  = de / fSimuParam->GetGeVToCharge();
      if (fSimuParam->GetPixLorentzDrift()) x += y*fTanLorAng;
      SpreadCharge2D(x,z,y,ix,iz,el,tof,idtrack,h);
    } // end if st>0.0
    
  } // Loop over all hits h
  
  // Coupling
  int nd = fSensMap->GetEntriesUnsorted(); // use unsorted access when possible, since it is faster
  AliITSUSDigit* dg = 0;
  switch (fSimuParam->GetPixCouplingOption()) {
  case AliITSUSimuParam::kNewCouplingPix :
    for (i=nd;i--;) {
      dg = (AliITSUSDigit*)fSensMap->AtUnsorted(i);
      if (fSensMap->IsDisabled(dg)) continue;
      SetCoupling(dg,idtrack,h);
    } 
  case AliITSUSimuParam::kOldCouplingPix:
    for (i=nd;i--;) {
      dg = (AliITSUSDigit*)fSensMap->AtUnsorted(i);
      if (fSensMap->IsDisabled(dg)) continue;
      SetCouplingOld(dg,idtrack,h);
    } 
    break;
  default:
    break;
  } // end switch
  if (GetDebug(2)) AliInfo(Form("Finished fCoupling=%d",fSimuParam->GetPixCouplingOption()));
}

//______________________________________________________________________
void AliITSUSimulationPix::SpreadCharge2D(Double_t x0,Double_t z0, Double_t dy, Int_t ix0,Int_t iz0,
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
  int nx = GetResponseParam()->GetParameter(kSpreadFunParamNXoffs);
  int nz = GetResponseParam()->GetParameter(kSpreadFunParamNZoffs);
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
      if ( !(((AliITSUSimulationPix*)this)->*AliITSUSimulationPix::fROTimeFun)(ix,iz,tof) ) continue;
      //
      fSeg->DetToLocal(ix,iz,x,z); // pixel center
      xdioshift = zdioshift = 0;
      CalcDiodeShiftInPixel(ix,iz,xdioshift,zdioshift);    // Check and apply diode shift if needed
      double dxi = 0.5*fSeg->Dpx(ix);
      double dzi = 0.5*fSeg->Dpz(iz);
      x1  = (x + xdioshift) - x0;   // calculate distance of cell boundaries from injection center
      z1  = (z + zdioshift) - z0;
      x2  = x1 + dxi; // Upper
      x1 -= dxi;      // Lower
      z2  = z1 + dzi; // Upper
      z1 -= dzi;      // Lower
      s   = el* (((AliITSUSimulationPix*)this)->*AliITSUSimulationPix::fSpreadFun)(dtIn);
      if (s>fSimuParam->GetPixMinElToAdd()) UpdateMapSignal(iz,ix,tID,hID,s);
    } // end for ix, iz
  //
}

//______________________________________________________________________
Double_t AliITSUSimulationPix::SpreadFunDoubleGauss2D(const Double_t *dtIn)
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
Double_t AliITSUSimulationPix::SpreadFrom2DHisto(const Double_t *dtIn)
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
Double_t AliITSUSimulationPix::SpreadFunGauss2D(const Double_t *dtIn)
{
  // calculate integral of charge in the cell with boundaries at X=dtIn[kCellX1]:dtIn[kCellX2] 
  // and Z=dtIn[kCellZ1]:dtIn[kCellZ2] 
  // The spread function is assumed to be gaussian in 2D
  // Parameters should be: mean0,sigma0
  return GausInt2D(fResponseParam->GetParameter(kG1SigX),  // sigX
		   fResponseParam->GetParameter(kG1SigZ),  // sigZ
		   dtIn[kCellX1]-fResponseParam->GetParameter(kG1MeanX),    // x1-xmean
		   dtIn[kCellX2]-fResponseParam->GetParameter(kG1MeanX),    // x2-xmean
		   dtIn[kCellZ1]-fResponseParam->GetParameter(kG1MeanZ),    // z1-zmean
		   dtIn[kCellZ2]-fResponseParam->GetParameter(kG1MeanZ));   // z2-zmean
  //
} 

//______________________________________________________________________
void AliITSUSimulationPix::RemoveDeadPixels() 
{
  // Removes dead pixels on each module (ladder)
  // This should be called before going from sdigits to digits (FrompListToDigits)
  
  AliITSUCalibrationPix* calObj = (AliITSUCalibrationPix*) GetCalibDead();
  if (!calObj) return;
  //
  if (calObj->IsBad()) {ClearMap(); return;} // whole module is masked
  //
  // remove single bad pixels one by one
  int nsingle = calObj->GetNrBadSingle();
  UInt_t col,row;
  for (int i=nsingle;i--;) {
    calObj->GetBadPixelSingle(i,row,col);
    fSensMap->DeleteItem(col,row);
  }
  int nsd = fSensMap->GetEntriesUnsorted();
  for (int isd=nsd;isd--;) {
    AliITSUSDigit* sd = (AliITSUSDigit*)fSensMap->AtUnsorted(isd);
    if (fSensMap->IsDisabled(sd)) continue;
    fSensMap->GetMapIndex(sd->GetUniqueID(),col,row);
    int chip = fSeg->GetChipFromChannel(0,col);
    //    if (calObj->IsChipMarkedBad(chip)) fSensMap->Disable(sd); // this will simple mark the hit as bad
    if (calObj->IsChipMarkedBad(chip)) fSensMap->DeleteItem(sd); // this will suppress hit in the sorted list
  }
  //
}

//______________________________________________________________________
void AliITSUSimulationPix::AddNoisyPixels() 
{
  // Adds noisy pixels on each module (ladder)
  // This should be called before going from sdigits to digits (FrompListToDigits)
  AliITSUCalibrationPix* calObj = (AliITSUCalibrationPix*) GetCalibNoisy();
  if (!calObj) { AliDebug(10,Form("  No Calib Object for Noise!!! ")); return;}
  for (Int_t i=calObj->GetNrBad(); i--;) UpdateMapNoise(calObj->GetBadColAt(i), calObj->GetBadRowAt(i), 
							10*fSimuParam->GetPixThreshold(fModule->GetIndex()));
  //
}

//______________________________________________________________________
void AliITSUSimulationPix::FrompListToDigits() 
{
  // add noise and electronics, perform the zero suppression and add the
  // digit to the list
  static AliITSU *aliITS = (AliITSU*)gAlice->GetModule("ITS");
  UInt_t ix,iz;
  Double_t sig;
  const Int_t    knmaxtrk=AliITSdigit::GetNTracks();
  static AliITSUDigitPix dig;
  // RS: in principle:
  // 1) for every pixel w/o hit we have to generate a noise and activate the pixel if the noise exceeds the threshold. 
  // 2) for every pixel with hit we should add random noise and check if the total signal exceeds the threshold
  // With many channels this will be too time consuming, hence I do the following
  // 1) Precalculate the probability that the nois alone will exceed the threshold (probnoisy). 
  // 2) Chose randomly empty pixels according to fakes rate
  // 3) For pixels having a hits apply the usual noise and compare with threshold
  //
  // RS may use for ordered random sample generation dl.acm.org/ft_gateway.cfm?id=356313&type=pdf
  //
  int maxInd = fSensMap->GetMaxIndex();
  int modId = fModule->GetIndex();
  //
  int nsd = fSensMap->GetEntries();
  Int_t prevID=0,curID=0;
  TArrayI ordSampleInd(100),ordSample(100);
  //
  double probNoisy,noiseSig,noiseMean,thresh = fSimuParam->GetPixThreshold(modId);
  fSimuParam->GetPixNoise(modId, noiseSig, noiseMean);
  probNoisy = AliITSUSimuParam::CalcProbNoiseOverThreshold(noiseMean,noiseSig,thresh); // prob. to have noise above threshold
  //
  for (int i=0;i<nsd;i++) {
    AliITSUSDigit* sd = (AliITSUSDigit*)fSensMap->At(i); // ordered in index
    if (fSensMap->IsDisabled(sd)) continue;
    curID = (int)sd->GetUniqueID();
    //
    if ( fSimuParam->GetPixAddNoisyFlag() ) { 
      CreateNoisyDigits(prevID,curID,probNoisy, noiseSig, noiseMean);
      prevID = curID+1;
      //
      // add noise also to sdigit with signal
      sd->AddNoise(AliITSUSimuParam::GenerateNoiseQFunction(0,noiseMean,noiseSig));
    }
    //
    if ((sig=sd->GetSumSignal())<=fSimuParam->GetPixThreshold(modId)) continue;
    if (Abs(sig)>2147483647.0) { //RS?
      //PH 2147483647 is the max. integer
      //PH This apparently is a problem which needs investigation
      AliWarning(Form("Too big or too small signal value %f",sig));
    }
    fSensMap->GetMapIndex(sd->GetUniqueID(),iz,ix);
    dig.SetCoord1(iz);
    dig.SetCoord2(ix);
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
    aliITS->AddSimDigit(AliITSUGeomTGeo::kDetTypePix, &dig);
  }
  // if needed, add noisy pixels with id from last real hit to maxID
  if (fSimuParam->GetPixAddNoisyFlag() && 
      (nsd>0  || ( nsd==0 && fSimuParam->GetPixNoiseInAllMod()) )) {
    AliDebug(10,Form("CreateNoisyDigits2( %d ,%d) module %d",prevID,maxInd,modId));
    CreateNoisyDigits(prevID,maxInd,probNoisy, noiseSig, noiseMean);
  }
  // 
}

//______________________________________________________________________
Int_t AliITSUSimulationPix::CreateNoisyDigits(Int_t minID,Int_t maxID,double probNoisy, double noise, double base)
{//RS PAYATT2
  // create random noisy digits above threshold within id range [minID,maxID[
  // see FrompListToDigits for details
  //
  static AliITSU *aliITS = (AliITSU*)gAlice->GetModule("ITS");
  UInt_t ix,iz;
  static AliITSUDigitPix dig;
  static TArrayI ordSampleInd(100),ordSample(100); //RS!!! static is not thread-safe!!!
  const Int_t    knmaxtrk=AliITSdigit::GetNTracks();
  //
  Int_t ncand = 0;
  int npix = maxID-minID;
  if (npix<1 || (ncand=gRandom->Poisson(npix * fSimuParam->GetPixFakeRate() ))<1) return ncand; // decide how many noisy pixels will be added
  ncand = GenOrderedSample(npix,ncand,ordSample,ordSampleInd); 
  int* ordV = ordSample.GetArray();
  int* ordI = ordSampleInd.GetArray();
  for (int j=0;j<ncand;j++) {
    fSensMap->GetMapIndex((UInt_t)ordV[ordI[j]],iz,ix);   // create noisy digit
    dig.SetCoord1(iz);
    dig.SetCoord2(ix);
    dig.SetSignal(1);
    dig.SetSignalPix((Int_t)AliITSUSimuParam::GenerateNoiseQFunction(probNoisy,base,noise));
    for (int k=knmaxtrk;k--;) {
      dig.SetTrack(k,-3);
      dig.SetHit(k,-1);
    }
    aliITS->AddSimDigit(AliITSUGeomTGeo::kDetTypePix,&dig);
    AliDebug(10,Form("Add noisy pixel %d(%d/%d) Noise=%d",ordV[ordI[j]],iz,ix,dig.GetSignalPix()));
  }
  return ncand;
}

//______________________________________________________________________
void AliITSUSimulationPix::SetCoupling(AliITSUSDigit* old, Int_t ntrack, Int_t idhit) 
{
  //  Take into account the coupling between adiacent pixels.
  //  The parameters probcol and probrow are the probability of the
  //  signal in one pixel shared in the two adjacent pixels along
  //  the column and row direction, respectively.
  //  Note pList is goten via GetMap() and module is not need any more.
  //  Otherwise it is identical to that coded by Tiziano Virgili (BSN).
  //Begin_Html
  /*
    <img src="picts/ITS/barimodel_3.gif">
     </pre>
     <br clear=left>
     <font size=+2 color=red>
     <a href="mailto:tiziano.virgili@cern.ch"></a>.
     </font>
     <pre>
   */
   //End_Html
   // Inputs:
  // old                  existing AliITSUSDigit
  // Int_t ntrack         track incex number
  // Int_t idhit          hit index number
  UInt_t col,row;
  Double_t pulse1,pulse2;
  Double_t couplR=0.0,couplC=0.0;
  Double_t xr=0.;
  //
  fSensMap->GetMapIndex(old->GetUniqueID(),col,row);
  fSimuParam->GetPixCouplingParam(couplC,couplR);
  if (GetDebug(2)) AliInfo(Form("(col=%d,row=%d,ntrack=%d,idhit=%d)  couplC=%e couplR=%e",
				col,row,ntrack,idhit,couplC,couplR));
  pulse2 = pulse1 = old->GetSignal();
  if (pulse1<fSimuParam->GetPixMinElToAdd()) return; // too small signal
  for (Int_t isign=-1;isign<=1;isign+=2) {
    //
    // loop in col direction
    int j1 = int(col) + isign;
    xr = gRandom->Rndm();
    if ( !((j1<0) || (j1>fSeg->Npz()-1) || (xr>couplC)) ) UpdateMapSignal(UInt_t(j1),row,ntrack,idhit,pulse1);
    //
    // loop in row direction
    int j2 = int(row) + isign;
    xr = gRandom->Rndm();
    if ( !((j2<0) || (j2>fSeg->Npx()-1) || (xr>couplR)) ) UpdateMapSignal(col,UInt_t(j2),ntrack,idhit,pulse2);
  } 
  //
}

//______________________________________________________________________
void AliITSUSimulationPix::SetCouplingOld(AliITSUSDigit* old, Int_t ntrack,Int_t idhit) 
{
  //  Take into account the coupling between adiacent pixels.
  //  The parameters probcol and probrow are the fractions of the
  //  signal in one pixel shared in the two adjacent pixels along
  //  the column and row direction, respectively.
  //Begin_Html
  /*
    <img src="picts/ITS/barimodel_3.gif">
    </pre>
    <br clear=left>
    <font size=+2 color=red>
    <a href="mailto:Rocco.Caliandro@ba.infn.it"></a>.
    </font>
    <pre>
  */
  //End_Html
  // Inputs:
  // old            existing AliITSUSDigit
  // ntrack         track incex number
  // idhit          hit index number
  // module         module number
  //
  UInt_t col,row;
  Int_t modId = fModule->GetIndex();
  Double_t pulse1,pulse2;
  Double_t couplR=0.0,couplC=0.0;
  //
  fSensMap->GetMapIndex(old->GetUniqueID(),col,row);
  fSimuParam->GetPixCouplingParam(couplC,couplR);
  if (GetDebug(3)) AliInfo(Form("(col=%d,row=%d,ntrack=%d,idhit=%d)  couplC=%e couplR=%e",
				col,row,ntrack,idhit,couplC,couplR));
 //
 if (old->GetSignal()<fSimuParam->GetPixMinElToAdd()) return; // too small signal
 for (Int_t isign=-1;isign<=1;isign+=2) {// loop in col direction
   pulse2 = pulse1 = old->GetSignal();
   //
   int j1 = int(col)+isign;
   pulse1 *= couplC;    
   if ((j1<0)||(j1>fSeg->Npz()-1)||(pulse1<fSimuParam->GetPixThreshold(modId))) pulse1 = old->GetSignal();
   else UpdateMapSignal(UInt_t(j1),row,ntrack,idhit,pulse1);
   
   // loop in row direction
   int j2 = int(row) + isign;
   pulse2 *= couplR;
   if ((j2<0)||(j2>(fSeg->Npx()-1))||(pulse2<fSimuParam->GetPixThreshold(modId))) pulse2 = old->GetSignal();
   else UpdateMapSignal(col,UInt_t(j2),ntrack,idhit,pulse2);
 } // for isign
}

//______________________________________________________________________
void AliITSUSimulationPix::GenerateReadOutCycleOffset()
{
  // Generate randomly the strobe
  // phase w.r.t to the LHC clock
  fReadOutCycleOffset = fReadOutCycleLength*gRandom->Rndm();
  // fReadOutCycleOffset = 25e-9*gRandom->Rndm(); // clm: I think this way we shift too much 10-30 us! The global shift should be between the BCs?!
  // RS: 25 ns is too small number, the staggering will not work. Let's at the moment keep fully random shift (still, no particle from correct
  // collision will be lost) untill real number is specified
 //
}

//______________________________________________________________________
void AliITSUSimulationPix::SetResponseParam(AliITSUParamList* resp)
{
  // attach response parameterisation data
  fResponseParam = resp;
  //
  int spreadID = Nint(fResponseParam->GetParameter(AliITSUSimulationPix::kChargeSpreadType));
  const char* hname = 0;
  fSpread2DHisto = 0;
  //
  switch (spreadID) {
    //
  case kSpreadFunHisto: 
    fSpreadFun = &AliITSUSimulationPix::SpreadFrom2DHisto;
    hname = fResponseParam->GetParName(AliITSUSimulationPix::kChargeSpreadType);
    if (!(fSpread2DHisto=(TH2*)fResponseParam->GetParamObject(hname))) 
      AliFatal(Form("Did not find 2D histo %s for charge spread parameterization",hname));
    break;
    //
  case kSpreadFunDoubleGauss2D: 
    fSpreadFun = &AliITSUSimulationPix::SpreadFunDoubleGauss2D; 
    break;
    //
  case kSpreadFunGauss2D: 
    fSpreadFun = &AliITSUSimulationPix::SpreadFunGauss2D;       
    break;
    //
  default: AliFatal(Form("Did not find requested spread function id=%d",spreadID));
  }
  //
  int readoutType = Nint(fResponseParam->GetParameter(kReadOutSchemeType));
  switch (readoutType) {
  case kReadOutStrobe: 
    fROTimeFun = &AliITSUSimulationPix::IsHitInReadOutWindow;
    break;
  case kReadOutRollingShutter: 
    fROTimeFun = &AliITSUSimulationPix::IsHitInReadOutWindowRollingShutter;
    break;
  default: AliFatal(Form("Did not find requested readout time type id=%d",readoutType));
  }
  //___ Set the Rolling Shutter read-out window 
  fReadOutCycleLength = fResponseParam->GetParameter(kReadOutCycleLength);
  //___ Pixel discrimination threshold, and the S/N cut
  fSimuParam->SetPixThreshold(fResponseParam->GetParameter(kPixNoiseMPV) *fResponseParam->GetParameter(kPixSNDisrcCut) , fResponseParam->GetParameter(kPixSNDisrcCut),-1); //for all modules
  //___ Minimum number of electrons to add 
  fSimuParam->SetPixMinElToAdd(fResponseParam->GetParameter(kPixMinElToAdd));
  //___ Set the Pixel Noise MPV and Sigma (the noise distribution is Landau not Gauss due to RTN)
  fSimuParam->SetPixNoise( fResponseParam->GetParameter(kPixNoiseMPV), fResponseParam->GetParameter(kPixNoiseSigma), -1); //for all modules
  //___ Pixel fake hit rate
  fSimuParam->SetPixFakeRate( fResponseParam->GetParameter(kPixFakeRate) );
  //___ To apply the noise or not
  if (  fResponseParam->GetParameter(kPixNoiseIsOn) > 0.01)  fSimuParam->SetPixAddNoisyFlag(kTRUE);
  else fSimuParam->SetPixAddNoisyFlag(kFALSE);
  //
  if(fResponseParam->GetParameter(kPixNoiseInAllMod) > 0.01 ) fSimuParam->SetPixNoiseInAllMod(kTRUE);
  else fSimuParam->SetPixNoiseInAllMod(kFALSE);
  //
  //  Double_t vGeVToQ = fSimuParam->GetGeVToCharge();
  fGlobalChargeScale = fResponseParam->GetParameter(kSpreadFunGlobalQScale);
    
  AliDebug(10,Form("=============== Setting the response start ============================"));
  AliDebug(10,Form("=============== RO type: %d",readoutType));
  AliDebug(10,Form("=============== RO cycle lenght: %lf",fReadOutCycleLength));
  AliDebug(10,Form("=============== Noise MPV: %lf",fResponseParam->GetParameter(kPixNoiseMPV)));
  AliDebug(10,Form("=============== Noise Sigma: %lf",fResponseParam->GetParameter(kPixNoiseSigma)));
  AliDebug(10,Form("=============== Fake rate: %lf",fResponseParam->GetParameter(kPixFakeRate)));
  AliDebug(10,Form("=============== Noise On/Off: %d",fSimuParam->GetPixAddNoisyFlag()));
  AliDebug(10,Form("=============== Noise in all mod on/off: %d",fSimuParam->GetPixNoiseInAllMod()));
  AliDebug(10,Form("=============== Global Charge scale: %lf",fGlobalChargeScale));
  AliDebug(10,Form("=============== Setting the response done  ============================"));
  
  
  
}

//______________________________________________________________________
Bool_t AliITSUSimulationPix::IsHitInReadOutWindowRollingShutter(Int_t row, Int_t col, Double_t hitTime)
{
  //
  // Check whether the hit is in the read out window of the given column/row of the sensor
  // hitTime is the time of the subhit (hit is divided to nstep charge deposit) in seconds
  // globalPhaseShift gives the start of the RO for the cycle in pixel wrt the LHC clock
  // GetRollingShutterWindow give the with of the rolling shutter read out window
  //
  double timePerRow = fReadOutCycleLength / fSeg->Npx();
  double tmax = fReadOutCycleOffset + timePerRow*(row+1);
  double tmin = tmax - fReadOutCycleLength;
  AliDebug(3,Form("Rolling shutter at row%d/col%d: particle time: %e, tmin: %e : %e",row,col,hitTime,tmin,tmax));
  return (hitTime<tmin || hitTime>tmax) ? kFALSE : kTRUE;
  //  
}

//______________________________________________________________________
Bool_t AliITSUSimulationPix::IsHitInReadOutWindow(Int_t row, Int_t col, Double_t hitTime)
{
  //
  // Check whether the hit is in the read out window of the given column/row of the sensor
  //
  AliDebug(3,Form("Strobe readout: row%d/col%d: particle time: %e, tmin: %e, tmax %e",
		  row,col,hitTime,fReadOutCycleOffset,fReadOutCycleOffset+fReadOutCycleLength));
  hitTime -= fReadOutCycleOffset+0.5*fReadOutCycleLength;
  return (hitTime<0 || hitTime>fReadOutCycleLength) ? kFALSE : kTRUE;
  //  
}

//_______________________________________________________________________
void AliITSUSimulationPix::CalcDiodeShiftInPixel(Int_t xlin, Int_t zcol, Float_t &x, Float_t &)
{
  //
  // Calculates the shift of the diode wrt the geometric center of the pixel.
  // It is needed for staggerred pixel layout or double diode pixels with assymetric center
  // The shift can depend on the column or line or both...
  // The x and z are passed in cm 
  //
  
  TString parTitle = fResponseParam->GetTitle();
  
  // M32terP31 is staggered the diode shift within pixel depends on the column
  if ( parTitle.Contains("M32terP31") ) 
  {
    if ( zcol%2 == 0 )    x += 0.30 * fSeg->Dpx(xlin);
    else                  x -= 0.19 * fSeg->Dpx(xlin);
  }
  
 
}
//_______________________________________________________________________
