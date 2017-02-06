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

#include <Riostream.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TString.h>
#include <TGeoMatrix.h>
#include "AliITS.h"
#include "AliITSdigitSPD.h"
#include "AliITShit.h"
#include "AliITSmodule.h"
#include "AliITSpList.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsimulationSPD.h"
#include "AliITSgeomTGeo.h"
#include "AliLog.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliMathBase.h"

//#define DEBUG

using std::endl;
using std::cout;
ClassImp(AliITSsimulationSPD)
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
 //
 // AliITSsimulationSPD is to do the simulation of SPDs.
 //
 ////////////////////////////////////////////////////////////////////////

 //______________________________________________________________________
 AliITSsimulationSPD::AliITSsimulationSPD():
  AliITSsimulation(),
  fHis(0),
  fSPDname(),
  fCoupling(),
  fLorentz(kFALSE),
  fTanLorAng(0),
  fStrobe(kTRUE),
  fStrobeLenght(4),
  fStrobePhase(-12.5e-9){
   // Default constructor.
   // Inputs:
   //    none.
   // Outputs:
   //    none.
   // Return:
   //    A default constructed AliITSsimulationSPD class.

   AliDebug(1,Form("Calling default constructor"));
   //    Init();
  }
//______________________________________________________________________
AliITSsimulationSPD::AliITSsimulationSPD(AliITSDetTypeSim *dettyp):
 AliITSsimulation(dettyp),
 fHis(0),
 fSPDname(),
 fCoupling(),
 fLorentz(kFALSE),
 fTanLorAng(0),
 fStrobe(kTRUE),
 fStrobeLenght(4),
 fStrobePhase(-12.5e-9){
  // standard constructor
  // Inputs:
  //    AliITSsegmentation *seg  A pointer to the segmentation class
  //                             to be used for this simulation
  //    AliITSCalibration     *resp A pointer to the responce class to
  //                             be used for this simulation
  // Outputs:
  //    none.
  // Return:
  //    A default constructed AliITSsimulationSPD class.

  AliDebug(1,Form("Calling standard constructor "));
  Init();
 }
//______________________________________________________________________
void AliITSsimulationSPD::Init(){
 // Initilization
 // Inputs:
 //    none.
 // Outputs:
 //    none.
 // Return:
 //    none.
 const Double_t kmictocm = 1.0e-4; // convert microns to cm.

 SetModuleNumber(0);
 SetEventNumber(0);
 SetMap(new AliITSpList(GetNPixelsZ(),GetNPixelsX()));
 AliITSSimuParam* simpar = fDetType->GetSimuParam();
 AliITSsegmentationSPD* seg = (AliITSsegmentationSPD*)GetSegmentationModel(0);
 Double_t bias = simpar->GetSPDBiasVoltage();
 //    cout << "Bias Voltage --> " << bias << endl; // dom    
 simpar->SetDistanceOverVoltage(kmictocm*seg->Dy(),bias);
 // set kind of coupling ("old" or "new")
 char opt[20];
 simpar->GetSPDCouplingOption(opt);
 char *old = strstr(opt,"old");
 if (old) {
  fCoupling=2;
 } else {
  fCoupling=1;
 } // end if
 SetLorentzDrift(simpar->GetSPDLorentzDrift());
 if (fLorentz) SetTanLorAngle(simpar->GetSPDLorentzHoleWeight());
 //SetStrobeGeneration(kFALSE);
 if (fStrobe) GenerateStrobePhase();
}
//______________________________________________________________________
Bool_t AliITSsimulationSPD::SetTanLorAngle(Double_t WeightHole) {
 // This function set the Tangent of the Lorentz angle. 
 // A weighted average is used for electrons and holes 
 // Input: Double_t WeightHole: wheight for hole: it should be in the range [0,1]
 // output: Bool_t : kTRUE in case of success
 //
 if(!fDetType) {
  AliError("AliITSsimulationSPD::SetTanLorAngle: AliITSDetTypeSim* fDetType not set ");
  return kFALSE;}
  if(WeightHole<0) {
   WeightHole=0.;
   AliWarning("AliITSsimulationSPD::SetTanLorAngle: You have asked for negative Hole weight");
   AliWarning("AliITSsimulationSPD::SetTanLorAngle: I'm going to use only electrons");
  }
  if(WeightHole>1) {
   WeightHole=1.;
   AliWarning("AliITSsimulationSPD::SetTanLorAngle: You have asked for weight > 1");
   AliWarning("AliITSsimulationSPD::SetTanLorAngle: I'm going to use only holes");
  }
  Double_t WeightEle=1.-WeightHole;
  AliITSSimuParam* simpar = fDetType->GetSimuParam();
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!fld) AliFatal("The field is not initialized");
  Double_t bz = fld->SolenoidField();
  // RS: "-" sign to account for the fact that e.g. holes move against local Y axis 
  fTanLorAng = -TMath::Tan(WeightHole*simpar->LorentzAngleHole(bz) +
			   WeightEle*simpar->LorentzAngleElectron(bz));
  return kTRUE;
}
//______________________________________________________________________
AliITSsimulationSPD::~AliITSsimulationSPD(){
 // destructor
 // Inputs:
 //    none.
 // Outputs:
 //    none.
 // Return:
 //     none.

 if (fHis) {
  fHis->Delete(); 
  delete fHis;     
 } // end if fHis
}
//______________________________________________________________________
AliITSsimulationSPD::AliITSsimulationSPD(const 
  AliITSsimulationSPD 
  &s) : AliITSsimulation(s),
 fHis(s.fHis),
 fSPDname(s.fSPDname),
 fCoupling(s.fCoupling),
 fLorentz(s.fLorentz),
 fTanLorAng(s.fTanLorAng),
 fStrobe(s.fStrobe),
 fStrobeLenght(s.fStrobeLenght),
 fStrobePhase(s.fStrobePhase){
  //     Copy Constructor
  // Inputs:
  //    AliITSsimulationSPD &s The original class for which
  //                                this class is a copy of
  // Outputs:
  //    none.
  // Return:

 }
//______________________________________________________________________
AliITSsimulationSPD&  AliITSsimulationSPD::operator=(const 
  AliITSsimulationSPD &s){
 //    Assignment operator
 // Inputs:
 //    AliITSsimulationSPD &s The original class for which
 //                                this class is a copy of
 // Outputs:
 //    none.
 // Return:

 if(&s == this) return *this;
 this->fHis = s.fHis;
 fCoupling  = s.fCoupling;
 fSPDname   = s.fSPDname;
 fLorentz   = s.fLorentz;
 fTanLorAng = s.fTanLorAng;
 fStrobe       = s.fStrobe;
 fStrobeLenght = s.fStrobeLenght;
 fStrobePhase  = s.fStrobePhase;
 return *this;
}
/*
//______________________________________________________________________
AliITSsimulation&  AliITSsimulationSPD::operator=(const 
AliITSsimulation &s){
//    Assignment operator
// Inputs:
//    AliITSsimulationSPD &s The original class for which
//                                this class is a copy of
// Outputs:
//    none.
// Return:

if(&s == this) return *this;
Error("AliITSsimulationSPD","Not allowed to make a = with "
"AliITSsimulationSPD","Using default creater instead");

return *this;
}
*/
//______________________________________________________________________
void AliITSsimulationSPD::InitSimulationModule(Int_t module, Int_t event){
 //  This function creates maps to build the list of tracks for each
 //  summable digit. Inputs defined by base class.
 //  Inputs:
 //    Int_t module   // Module number to be simulated
 //    Int_t event    // Event number to be simulated
 //  Outputs:
 //    none
 //  Returns:
 //    none

 AliDebug(1,Form("(module=%d,event=%d)",module,event));
 const Double_t kmictocm = 1.0e-4; // convert microns to cm.
 AliITSSimuParam* simpar = fDetType->GetSimuParam();
 AliITSsegmentationSPD* seg = (AliITSsegmentationSPD*)GetSegmentationModel(0);
 SetModuleNumber(module);
 SetEventNumber(event);
 simpar->SetDistanceOverVoltage(kmictocm*seg->Dy(),simpar->GetSPDBiasVoltage(module)); 
 ClearMap();
}
//_____________________________________________________________________
void AliITSsimulationSPD::SDigitiseModule(AliITSmodule *mod,Int_t,
  Int_t event){
 //  This function begins the work of creating S-Digits.  Inputs defined
 //  by base class.
 //  Inputs:
 //    AliITSmodule *mod  //  module
 //    Int_t              //  not used
 //    Int_t event        //  Event number
 //  Outputs:
 //    none
 //  Return:
 //    test              //  test returns kTRUE if the module contained hits
 //                      //  test returns kFALSE if it did not contain hits

 AliDebug(1,Form("(mod=%p, ,event=%d)",mod,event));
 if(!(mod->GetNhits())){
  AliDebug(1,Form("In event %d module %d there are %d hits returning.",
     event, mod->GetIndex(),mod->GetNhits()));
  return;// if module has no hits don't create Sdigits
 } // end if
 SetModuleNumber(mod->GetIndex());
 if (fStrobe) if(event != GetEventNumber()) GenerateStrobePhase(); 
 SetEventNumber(event);
 InitSimulationModule( GetModuleNumber() , event );
 // HitToSDigit(mod);
 HitToSDigitFast(mod);
 if (fDetType->GetSimuParam()->GetSPDAddNoisyFlag())   AddNoisyPixels();
 if (fDetType->GetSimuParam()->GetSPDRemoveDeadFlag()) RemoveDeadPixels();

 //    cout << "After Remove in SDigitiseModule !!!!!" << endl; // dom
 //    cout << "Module " << mod->GetIndex() << " Event " << event << endl; // dom
 WriteSDigits();
 ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSPD::WriteSDigits(){
 //  This function adds each S-Digit to pList
 //  Inputs:
 //    none.
 //  Outputs:
 //    none.
 //  Return:
 //    none
 Int_t ix, nix, iz, niz;
 static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");

 AliDebug(1,Form("Writing SDigits for module %d",GetModuleNumber()));
 //    cout << "WriteSDigits for module " << GetModuleNumber() << endl; // dom
 GetMap()->GetMaxMapIndex(niz, nix);
 for(iz=0; iz<niz; iz++)for(ix=0; ix<nix; ix++){
  if(GetMap()->GetSignalOnly(iz,ix)>0.0){
   //            cout << " Signal gt 0  iz ix " << iz << ix << " Module " << GetModuleNumber() << endl; // dom
   aliITS->AddSumDigit(*(GetMap()->GetpListItem(iz,ix)));
   if(AliDebugLevel()>0) {
    AliDebug(1,Form("%d, %d",iz,ix));
    cout << *(GetMap()->GetpListItem(iz,ix)) << endl;
   } // end if GetDebug
  } // end if GetMap()->GetSignalOnly(iz,ix)>0.0
 } // end for iz,ix
 return; 
}
//______________________________________________________________________
void AliITSsimulationSPD::FinishSDigitiseModule(){
 //  This function calls SDigitsToDigits which creates Digits from SDigits
 //  Inputs:
 //    none
 //  Outputs:
 //    none
 //  Return
 //    none

 AliDebug(1,"()");
 //    cout << "FinishSDigitiseModule for module " << GetModuleNumber() << endl; // dom
 FrompListToDigits(); // Charge To Signal both adds noise and
 ClearMap();
 return;
}
//______________________________________________________________________
void AliITSsimulationSPD::DigitiseModule(AliITSmodule *mod,Int_t,
  Int_t event){
 //  This function creates Digits straight from the hits and then adds
 //  electronic noise to the digits before adding them to pList
 //  Each of the input variables is passed along to HitToSDigit
 //  Inputs:
 //    AliITSmodule *mod     module
 //    Int_t                 Dummy.
 //    Int_t                 Dummy
 //  Outputs:
 //     none.
 //  Return:
 //    none.

 if (fStrobe) if(event != GetEventNumber()) GenerateStrobePhase();
 AliDebug(1,Form("(mod=%p,,0)",mod));
 // HitToSDigit(mod);
 InitSimulationModule( mod->GetIndex(), event );
 HitToSDigitFast(mod);

 if (fDetType->GetSimuParam()->GetSPDAddNoisyFlag())   AddNoisyPixels();
 if (fDetType->GetSimuParam()->GetSPDRemoveDeadFlag()) RemoveDeadPixels();
 //    cout << "After Remove in DigitiseModule in module " << mod->GetIndex() << endl; // dom
 FrompListToDigits();
 ClearMap();
}
//______________________________________________________________________
void AliITSsimulationSPD::HitToSDigit(AliITSmodule *mod){
 // Does the charge distributions using Gaussian diffusion charge charing.
 // Inputs:
 //    AliITSmodule *mod  Pointer to this module
 // Output:
 //    none.
 // Return:
 //    none.

 const Double_t kmictocm = 1.0e-4; // convert microns to cm.
 const Double_t kBunchLenght = 25e-9; // LHC clock
 TObjArray *hits = mod->GetHits();
 Int_t nhits = hits->GetEntriesFast();
 Int_t h,ix,iz,i;
 Int_t idtrack;
 Double_t x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0,de=0.0,ld=0.0;
 Double_t x,y,z,t,tp,st,dt=0.2,el,sig,sigx,sigz,fda;
 AliITSsegmentationSPD* seg = (AliITSsegmentationSPD*)GetSegmentationModel(0);
 AliITSSimuParam *simpar = fDetType->GetSimuParam();
 Double_t thick = 0.5*kmictocm*seg->Dy();  // Half Thickness
 simpar->GetSPDSigmaDiffusionAsymmetry(fda);  //

 // At the implementation time of the SPD readout strobe simulation for pixel hits and fastor separately,
 // the readout window is opened in the time range [-100 ns,200 ns] where 0 is the collision time, whereas 
 // the fastor window is opened in the range [0, 100 ns] where 0 is the collision time. 
 const Int_t *hitWin = simpar->GetSPDHitStrobe(); // left and right limit of the window.  Default is [-1,2] = [-100 ns, 200 ns]
 const Int_t *foWin = simpar->GetSPDFoStrobe(); // left and right limit of the window. Default is [0,1] = [0,100 ns]

 // further comment : to account for the 4 possible BC within 100 ns where the collision may occur, the hit and fo windows are 
 // shifted back and forth according to the fStrobePhase value that changes from 12.5 ns up to 87.5 ns. This shift 
 // should simulate the 4 possible bcmod4 BC where the collision occurs within the SPD clock of 100 ns


 AliDebug(1,Form("(mod=%p) fCoupling=%d",mod,fCoupling));
 if(nhits<=0) return;
 Double_t hitTimeMin = hitWin[0]*(Double_t)fStrobeLenght*kBunchLenght + (Double_t)fStrobePhase;
 Double_t hitTimeMax = hitWin[1]*(Double_t)fStrobeLenght*kBunchLenght + (Double_t)fStrobePhase;
 Double_t foTimeMin = foWin[0]*(Double_t)fStrobeLenght*kBunchLenght + (Double_t)fStrobePhase;
 Double_t foTimeMax = foWin[1]*(Double_t)fStrobeLenght*kBunchLenght + (Double_t)fStrobePhase;
 // if module local Z is inverted wrt lab Z, it will see opposite B in its frame
 double tanLorAngMod = AliITSgeomTGeo::GetMatrix(fModule)->GetRotationMatrix()[8]>0 ? fTanLorAng:-fTanLorAng;
 
 for(h=0;h<nhits;h++){
  if(AliDebugLevel()>0) {
   AliDebug(1,Form("Hits, %d", h));
   cout << *(mod->GetHit(h)) << endl;
  } // end if GetDebug

  // Check if the hit is inside readout window
  if (fStrobe){
   if ((mod->GetHit(h)->GetTOF() < hitTimeMin) ||  (mod->GetHit(h)->GetTOF() > hitTimeMax)) {
    continue;
   }
  }
  if(!mod->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,idtrack)) continue;
  st = TMath::Sqrt(x1*x1+y1*y1+z1*z1);
  if(st>0.0){
   st = (Double_t)((Int_t)(st/kmictocm)); // number of microns
   if(st<=1.0) st = 1.0;
   dt = 1.0/st;
   for(t=0.0;t<1.0;t+=dt){ // Integrate over t
    tp  = t+0.5*dt;
    x   = x0+x1*tp;
    y   = y0+y1*tp;
    z   = z0+z1*tp;
    if(!(seg->LocalToDet(x,z,ix,iz))) continue; // outside
    //el  = res->GeVToCharge((Double_t)(dt*de));
    el  = dt * de / simpar->GetGeVToCharge();
    if(GetDebug(1)){
     if(el<=0.0) cout<<"el="<<el<<" dt="<<dt
      <<" de="<<de<<endl;
    } // end if GetDebug
    sig = simpar->SigmaDiffusion1D(TMath::Abs(thick + y)); 
    sigx=sig;
    sigz=sig*fda;
    if (fLorentz) ld=(y+thick)*tanLorAngMod;
    SpreadChargeAsym(x,z,ix,iz,el,sigx,sigz,ld,idtrack,h);
   } // end for t
  } else { // st == 0.0 deposit it at this point
   x   = x0;
   y   = y0;
   z   = z0;
   if(!(seg->LocalToDet(x,z,ix,iz))) continue; // outside
   //el  = res->GeVToCharge((Double_t)de);
   el  = de / simpar->GetGeVToCharge();
   sig = simpar->SigmaDiffusion1D(TMath::Abs(thick + y));
   sigx=sig;
   sigz=sig*fda;
   if (fLorentz) ld=(y+thick)*tanLorAngMod;
   SpreadChargeAsym(x,z,ix,iz,el,sigx,sigz,ld,idtrack,h);
  } // end if st>0.0

 } // Loop over all hits h

 // Coupling
 switch (fCoupling) {
  default:
   break;
  case 1: //case 3:
   for(i=0;i<GetMap()->GetEntries();i++) 
    if(GetMap()->GetpListItem(i)==0) continue;
    else{
     GetMap()->GetMapIndex(GetMap()->GetpListItem(i)->GetIndex(),iz,ix);
     SetCoupling(iz,ix);
    } // end for i
   break;
  case 2: // case 4:
   for(i=0;i<GetMap()->GetEntries();i++) 
    if(GetMap()->GetpListItem(i)==0) continue;
    else{
     GetMap()->GetMapIndex(GetMap()->GetpListItem(i)->GetIndex(),iz,ix);
     SetCouplingOld(iz,ix);
    } // end for i
   break;
 } // end switch

 // complete the FO information

 for(Int_t j=0;j<GetMap()->GetEntries();j++){
  if(GetMap()->GetpListItem(j)==0) continue;
  Int_t iiz, iix;
  Int_t pListIndex = GetMap()->GetpListItem(j)->GetIndex();
  GetMap()->GetMapIndex(pListIndex,iiz,iix);
  for(Int_t k=0; k<10; k++){
   Int_t hitNb = GetMap()->GetHit(iiz,iix,k);
   if(hitNb>-1) {
    if(GetMap()->GetpListItem(iiz,iix)->GetSignal(k)<0.001) {
     // if there is no signal, skip..
     continue;
    }

    if (fStrobe) {
     Float_t tofHit =  mod->GetHit(hitNb)->GetTOF();
     if((tofHit < foTimeMin) ||  (tofHit > foTimeMax)) {
     continue;
     }
    }
    // add flag to account for the fo readout strobe 
    GetMap()->GetpListItem(j)->SetIsInFoStrobe(k);
   }
  } 
 } // end loop over plists

 if(fStrobe && GetDebug(2)) Info("HitToSDigit","Finished Fo Readout Strobe check");

}
//______________________________________________________________________
void AliITSsimulationSPD::HitToSDigitFast(AliITSmodule *mod){
 // Does the charge distributions using Gaussian diffusion charge charing.    // Inputs:
 //    AliITSmodule *mod  Pointer to this module
 // Output:
 //    none.
 // Return:
 //    none.

 const Double_t kmictocm = 1.0e-4; // convert microns to cm.
 const Int_t kn10=10;
 const Double_t kti[kn10]={7.443716945e-3,2.166976971e-1,3.397047841e-1,
  4.325316833e-1,4.869532643e-1,5.130467358e-1,
  5.674683167e-1,6.602952159e-1,7.833023029e-1,
  9.255628306e-1};
 const Double_t kwi[kn10]={1.477621124e-1,1.346333597e-1,1.095431813e-1,
  7.472567455e-2,3.333567215e-2,3.333567215e-2,
  7.472567455e-2,1.095431813e-1,1.346333597e-1,
  1.477621124e-1};
 const Double_t kBunchLenght = 25e-9; // LHC clock
 TObjArray *hits = mod->GetHits();
 Int_t nhits = hits->GetEntriesFast();
 Int_t h,ix,iz,i;
 Int_t idtrack;
 Double_t x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0,de=0.0,ld=0.0;
 Double_t x,y,z,t,st,el,sig,sigx,sigz,fda;
 AliITSsegmentationSPD* seg = (AliITSsegmentationSPD*)GetSegmentationModel(0);
 AliITSSimuParam* simpar = fDetType->GetSimuParam();
 Double_t thick = 0.5*kmictocm*seg->Dy();  // Half thickness
 simpar->GetSPDSigmaDiffusionAsymmetry(fda);
 //    cout << "Half Thickness " << thick << endl;  // dom
 //    cout << "Diffusion asymm " << fda << endl;  // dom

 // At the implementation time of the SPD readout strobe simulation for pixel hits and fastor separately,
 // the readout window is opened in the time range [-100 ns,200 ns] where 0 is the collision time, whereas 
 // the fastor window is opened in the range [0, 100 ns] where 0 is the collision time. 
 const Int_t *hitWin = simpar->GetSPDHitStrobe(); // left and right limit of the window.  Default is [-1,2] = [-100 ns, 200 ns]
 const Int_t *foWin = simpar->GetSPDFoStrobe(); // left and right limit of the window. Default is [0,1] = [0,100 ns]

 // further comment : to account for the 4 possible BC within 100 ns where the collision may occur, the hit and fo windows are 
 // shifted back and forth according to the fStrobePhase value that changes from 12.5 ns up to 87.5 ns. This shift 
 // should simulate the 4 possible bcmod4 BC where the collision occurs within the SPD clock of 100 ns

 AliDebug(1,Form("(mod=%p) fCoupling=%d",mod,fCoupling));
 if(nhits<=0) return;

 Double_t hitTimeMin = hitWin[0]*(Double_t)fStrobeLenght*kBunchLenght + (Double_t)fStrobePhase;
 Double_t hitTimeMax = hitWin[1]*(Double_t)fStrobeLenght*kBunchLenght + (Double_t)fStrobePhase;
 Double_t foTimeMin = foWin[0]*(Double_t)fStrobeLenght*kBunchLenght + (Double_t)fStrobePhase;
 Double_t foTimeMax = foWin[1]*(Double_t)fStrobeLenght*kBunchLenght + (Double_t)fStrobePhase;
 // if module local Z is inverted wrt lab Z, it will see opposite B in its frame
 double tanLorAngMod = AliITSgeomTGeo::GetMatrix(fModule)->GetRotationMatrix()[8]>0 ? fTanLorAng:-fTanLorAng;
 
 for(h=0;h<nhits;h++){
  if(AliDebugLevel()>0) {
   AliDebug(1,Form("Hits, %d", h));
   cout << *(mod->GetHit(h)) << endl;
  } // end if GetDebug
  // Check if the hit is inside readout window
  if (fStrobe){   
   if ((mod->GetHit(h)->GetTOF() < hitTimeMin) ||  (mod->GetHit(h)->GetTOF() > hitTimeMax)) {
    continue;
   }
  }

  if(!mod->LineSegmentL(h,x0,x1,y0,y1,z0,z1,de,idtrack)) {
   continue;
  }
  st = TMath::Sqrt(x1*x1+y1*y1+z1*z1);
  if(st>0.0) for(i=0;i<kn10;i++){ // Integrate over t
   t   = kti[i];
   x   = x0+x1*t;
   y   = y0+y1*t;
   z   = z0+z1*t;
   if(!(seg->LocalToDet(x,z,ix,iz))) {
    continue; // outside
   } 
   el  = kwi[i]*de/simpar->GetGeVToCharge();
   if(GetDebug(1)){
    if(el<=0.0) cout<<"el="<<el<<" kwi["<<i<<"]="<<kwi[i]
     <<" de="<<de<<endl;
   } // end if GetDebug
   sig = simpar->SigmaDiffusion1D(TMath::Abs(thick + y));
   sigx=sig;
   sigz=sig*fda;
   if (fLorentz) ld=(y+thick)*tanLorAngMod;
   SpreadChargeAsym(x,z,ix,iz,el,sigx,sigz,ld,idtrack,h);
   //                cout << "sigx sigz " << sigx << " " << sigz << endl; // dom
  } // end for i // End Integrate over t
  else { // st == 0.0 deposit it at this point
   x   = x0;
   y   = y0;
   z   = z0;
   if(!(seg->LocalToDet(x,z,ix,iz))) {
    continue; // outside
   } 
   el  = de / simpar->GetGeVToCharge();
   sig = simpar->SigmaDiffusion1D(TMath::Abs(thick + y));
   sigx=sig;
   sigz=sig*fda;
   if (fLorentz) ld=(y+thick)*tanLorAngMod;
   SpreadChargeAsym(x,z,ix,iz,el,sigx,sigz,ld,idtrack,h);
  } // end if st>0.0

 } // Loop over all hits h

 // Coupling
 switch (fCoupling) {
  default:
   break;
  case 1: // case 3:
   for(i=0;i<GetMap()->GetEntries();i++)
    if(GetMap()->GetpListItem(i)==0) continue;
    else{
     GetMap()->GetMapIndex(GetMap()->GetpListItem(i)->GetIndex(),iz,ix);
     SetCoupling(iz,ix);
    } // end for i
   break;
  case 2: // case 4:
   for(i=0;i<GetMap()->GetEntries();i++)
    if(GetMap()->GetpListItem(i)==0) continue;
    else{
     GetMap()->GetMapIndex(GetMap()->GetpListItem(i)->GetIndex(),iz,ix);  
     SetCouplingOld(iz,ix);
    } // end for i
   break;
 } // end switch
 if(GetDebug(2))Info("HitToSDigit","Finished fCoupling=%d",fCoupling);

 for(Int_t j=0;j<GetMap()->GetEntries();j++){
  if(GetMap()->GetpListItem(j)==0) continue;
  Int_t iiz, iix;
  Int_t pListIndex = GetMap()->GetpListItem(j)->GetIndex();
  GetMap()->GetMapIndex(pListIndex,iiz,iix);
  for(Int_t k=0; k<10; k++){
   Int_t hitNb = GetMap()->GetHit(iiz,iix,k);
   if(hitNb>-1) {
    if(GetMap()->GetpListItem(iiz,iix)->GetSignal(k)<0.001) {
     // if there is no signal, skip..
     continue;
    }

    if (fStrobe) {
     Float_t tofHit =  mod->GetHit(hitNb)->GetTOF();
     if((tofHit < foTimeMin) ||  (tofHit > foTimeMax)) {
     continue;
     }
    }
    // add flag to account for the fo readout strobe 
    GetMap()->GetpListItem(j)->SetIsInFoStrobe(k);
   }
  } 
 } // end loop over plists

 if(fStrobe && GetDebug(2)) Info("HitToSDigit","Finished Fo Readout Strobe check");
}
//______________________________________________________________________
void AliITSsimulationSPD::SpreadCharge(Double_t x0,Double_t z0,
  Int_t ix0,Int_t iz0,
  Double_t el,Double_t sig,Double_t ld,
  Int_t t,Int_t hi){
 // Spreads the charge over neighboring cells. Assume charge is distributed
 // as charge(x,z) = (el/2*pi*sig*sig)*exp(-arg)
 // arg=((x-x0)*(x-x0)/2*sig*sig)+((z-z0*z-z0)/2*sig*sig)
 // if fLorentz=kTRUE, then x0=x0+ld (Lorentz drift taken into account)
 // Defined this way, the integral over all x and z is el.
 // Inputs:
 //    Double_t x0   x position of point where charge is liberated
 //    Double_t z0   z position of point where charge is liberated
 //    Int_t    ix0  row of cell corresponding to point x0
 //    Int_t    iz0  columb of cell corresponding to point z0
 //    Double_t el   number of electrons liberated in this step
 //    Double_t sig  Sigma difusion for this step (y0 dependent)
 //    Double_t ld   lorentz drift in x for this step (y0 dependent)
 //    Int_t    t    track number
 //    Int_t    ti   hit track index number
 //    Int_t    hi   hit "hit" index number
 // Outputs:
 //     none.
 // Return:
 //     none.
 const Int_t knx = 3,knz = 2;
 const Double_t kRoot2 = 1.414213562; // Sqrt(2).
 const Double_t kmictocm = 1.0e-4; // convert microns to cm.
 Int_t ix,iz,ixs,ixe,izs,ize;
 Float_t x,z;
 Double_t x1,x2,z1,z2,s,sp;
 AliITSsegmentationSPD* seg = (AliITSsegmentationSPD*)GetSegmentationModel(0);


 if(GetDebug(4)) Info("SpreadCharge","(x0=%e,z0=%e,ix0=%d,iz0=%d,el=%e,"
   "sig=%e,t=%d,i=%d)",x0,z0,ix0,iz0,el,sig,t,hi);
 if(sig<=0.0) { // if sig<=0 No diffusion to simulate.
  GetMap()->AddSignal(iz0,ix0,t,hi,GetModuleNumber(),el);
  if(GetDebug(2)){
   cout << "sig<=0.0=" << sig << endl;
  } // end if GetDebug
  return;
 } // end if
 sp = 1.0/(sig*kRoot2);
 if(GetDebug(2)){
  cout << "sig=" << sig << " sp=" << sp << endl;
 } // end if GetDebug
 ixs = TMath::Max(-knx+ix0,0);
 ixe = TMath::Min(knx+ix0,seg->Npx()-1);
 izs = TMath::Max(-knz+iz0,0);
 ize = TMath::Min(knz+iz0,seg->Npz()-1);
 for(ix=ixs;ix<=ixe;ix++) for(iz=izs;iz<=ize;iz++){
  seg->DetToLocal(ix,iz,x,z); // pixel center
  x1  = x;
  z1  = z;
  x2  = x1 + 0.5*kmictocm*seg->Dpx(ix); // Upper
  x1 -= 0.5*kmictocm*seg->Dpx(ix);  // Lower
  z2  = z1 + 0.5*kmictocm*seg->Dpz(iz); // Upper
  z1 -= 0.5*kmictocm*seg->Dpz(iz);  // Lower
  x1 -= x0+ld; // Distance from where track traveled (taking into account the Lorentz drift)
  x2 -= x0+ld; // Distance from where track traveled (taking into account the Lorentz drift)
  z1 -= z0; // Distance from where track traveled
  z2 -= z0; // Distance from where track traveled
  s   = 0.25; // Correction based on definision of Erfc
  s  *= AliMathBase::ErfcFast(sp*x1) - AliMathBase::ErfcFast(sp*x2);
  if(GetDebug(3)){
   cout <<"el="<<el<<" ix0="<<ix0<<" ix="<<ix<<" x0="<<x<<
    " iz0="<<iz0<<" iz="<<iz<<" z0="<<z<< 
    " sp*x1="<<sp*x1<<" sp*x2="<<sp*x2<<" s="<<s;
  } // end if GetDebug
  s  *= AliMathBase::ErfcFast(sp*z1) - AliMathBase::ErfcFast(sp*z2);
  if(GetDebug(3)){
   cout<<" sp*z1="<<sp*z1<<" sp*z2="<<sp*z2<<" s="<<s<< endl;
  } // end if GetDebug
  GetMap()->AddSignal(iz,ix,t,hi,GetModuleNumber(),s*el);
 } // end for ix, iz
}
//______________________________________________________________________
void AliITSsimulationSPD::SpreadChargeAsym(Double_t x0,Double_t z0,
  Int_t ix0,Int_t iz0,
  Double_t el,Double_t sigx,Double_t sigz,
  Double_t ld,Int_t t,Int_t hi){
 // Spreads the charge over neighboring cells. Assume charge is distributed
 // as charge(x,z) = (el/2*pi*sigx*sigz)*exp(-arg)
 // arg=((x-x0)*(x-x0)/2*sigx*sigx)+((z-z0*z-z0)/2*sigz*sigz)
 // if fLorentz=kTRUE, then x0=x0+ld (Lorentz drift taken into account)
 // Defined this way, the integral over all x and z is el.
 // Inputs:
 //    Double_t x0   x position of point where charge is liberated
 //    Double_t z0   z position of point where charge is liberated
 //    Int_t    ix0  row of cell corresponding to point x0
 //    Int_t    iz0  columb of cell corresponding to point z0
 //    Double_t el   number of electrons liberated in this step
 //    Double_t sigx Sigma difusion along x for this step (y0 dependent)
 //    Double_t sigz Sigma difusion along z for this step (y0 dependent)
 //    Double_t ld   lorentz drift in x for this stip (y0 dependent)
 //    Int_t    t    track number
 //    Int_t    ti   hit track index number
 //    Int_t    hi   hit "hit" index number
 // Outputs:
 //     none.
 // Return:
 //     none.
 const Int_t knx = 3,knz = 2;
 const Double_t kRoot2 = 1.414213562; // Sqrt(2).
 const Double_t kmictocm = 1.0e-4; // convert microns to cm.
 Int_t ix,iz,ixs,ixe,izs,ize;
 Float_t x,z;
 Double_t x1,x2,z1,z2,s,spx,spz;
 AliITSsegmentationSPD* seg = (AliITSsegmentationSPD*)GetSegmentationModel(0);


 if(GetDebug(4)) Info("SpreadChargeAsym","(x0=%e,z0=%e,ix0=%d,iz0=%d,el=%e,"
   "sigx=%e, sigz=%e, t=%d,i=%d)",x0,z0,ix0,iz0,el,sigx,sigz,t,hi);
 if(sigx<=0.0 || sigz<=0.0) { // if sig<=0 No diffusion to simulate.
  GetMap()->AddSignal(iz0,ix0,t,hi,GetModuleNumber(),el);
  if(GetDebug(2)){
   cout << "sigx<=0.0=" << sigx << endl;
   cout << "sigz<=0.0=" << sigz << endl;
  } // end if GetDebug
  return;
 } // end if
 spx = 1.0/(sigx*kRoot2);     spz = 1.0/(sigz*kRoot2);
 if(GetDebug(2)){
  cout << "sigx=" << sigx << " spx=" << spx << endl;
  cout << "sigz=" << sigz << " spz=" << spz << endl;
 } // end if GetDebug
 ixs = TMath::Max(-knx+ix0,0);
 ixe = TMath::Min(knx+ix0,seg->Npx()-1);
 izs = TMath::Max(-knz+iz0,0);
 ize = TMath::Min(knz+iz0,seg->Npz()-1);
 for(ix=ixs;ix<=ixe;ix++) for(iz=izs;iz<=ize;iz++){
  seg->DetToLocal(ix,iz,x,z); // pixel center
  x1  = x;
  z1  = z;
  x2  = x1 + 0.5*kmictocm*seg->Dpx(ix); // Upper
  x1 -= 0.5*kmictocm*seg->Dpx(ix);  // Lower
  z2  = z1 + 0.5*kmictocm*seg->Dpz(iz); // Upper
  z1 -= 0.5*kmictocm*seg->Dpz(iz);  // Lower
  x1 -= x0+ld; // Distance from where track traveled (taking into account the Lorentz drift)
  x2 -= x0+ld; // Distance from where track traveled (taking into account the Lorentz drift)
  z1 -= z0; // Distance from where track traveled
  z2 -= z0; // Distance from where track traveled
  s   = 0.25; // Correction based on definision of Erfc
  s  *= AliMathBase::ErfcFast(spx*x1) - AliMathBase::ErfcFast(spx*x2);
  if(GetDebug(3)){
   cout <<"el="<<el<<" ix0="<<ix0<<" ix="<<ix<<" x0="<<x<<
    " iz0="<<iz0<<" iz="<<iz<<" z0="<<z<< 
    " spx*x1="<<spx*x1<<" spx*x2="<<spx*x2<<" s="<<s;
  } // end if GetDebug
  s  *= AliMathBase::ErfcFast(spz*z1) - AliMathBase::ErfcFast(spz*z2);
  if(GetDebug(3)){
   cout<<" spz*z1="<<spz*z1<<" spz*z2="<<spz*z2<<" s="<<s<< endl;
  } // end if GetDebug
  GetMap()->AddSignal(iz,ix,t,hi,GetModuleNumber(),s*el);
 } // end for ix, iz
}
//______________________________________________________________________
void AliITSsimulationSPD::RemoveDeadPixels(){
 // Removes dead pixels on each module (ladder)
 // This should be called before going from sdigits to digits (FrompListToDigits)
 Int_t mod = GetModuleNumber();
 AliITSCalibrationSPD* calObj = (AliITSCalibrationSPD*) fDetType->GetCalibrationModel(mod);

 Int_t nrDead = calObj->GetNrBad();
 for (Int_t i=0; i<nrDead; i++) {
  GetMap()->DeleteHit(calObj->GetBadColAt(i), calObj->GetBadRowAt(i));
 }
}
//______________________________________________________________________
void AliITSsimulationSPD::AddNoisyPixels() {
 // Adds noisy pixels on each module (ladder)
 // This should be called before going from sdigits to digits (FrompListToDigits)
 Int_t mod = GetModuleNumber();
 AliITSCalibrationSPD* calObj = (AliITSCalibrationSPD*) fDetType->GetSPDNoisyModel(mod);

 Int_t nrNoisy = calObj->GetNrBad();
 for (Int_t i=0; i<nrNoisy; i++) {
  // adding 10 times the threshold will for sure make this pixel fire...
  GetMap()->AddNoise(calObj->GetBadColAt(i), calObj->GetBadRowAt(i), mod, 10*GetThreshold());
 }
}
//______________________________________________________________________
void AliITSsimulationSPD::FrompListToDigits(){
 // add noise and electronics, perform the zero suppression and add the
 // digit to the list
 // Inputs:
 //    none.
 // Outputs:
 //    none.
 // Return:
 //    none.
 static AliITS *aliITS = (AliITS*)gAlice->GetModule("ITS");
 Int_t j,ix,iz;
 Double_t  electronics;
 Double_t sig;
 const Int_t    knmaxtrk=AliITSdigit::GetNTracks();
 static AliITSdigitSPD dig;
 AliITSSimuParam *simpar = fDetType->GetSimuParam();
 if(GetDebug(1)) Info("FrompListToDigits","()");
 for(iz=0; iz<GetNPixelsZ(); iz++) for(ix=0; ix<GetNPixelsX(); ix++){
  // NEW (for the moment plugged by hand, in the future possibly read from Data Base)
  // here parametrize the efficiency of the pixel along the row for the test columns (1,9,17,25)
  //        if(iz==1 || iz == 9 || iz == 17 || iz == 25) {
  //        Double_t eff,p1=0.,p2=0.; 
  //        Double_t x=ix;
  //        switch (iz) {
  //          case 1:   p1=0.63460;p2=0.42438E-01;break;  
  //          case 9:   p1=0.41090;p2=0.75914E-01;break;
  //	  case 17:  p1=0.31883;p2=0.91502E-01;break;
  //	  case 25:  p1=0.48828;p2=0.57975E-01;break;
  //         } // end switch
  //          eff=1.-p1*exp(-p2*x);
  //          if (gRandom->Rndm() >= eff) continue;
  //        } // end  if 
  // END parametrize the efficiency
  // 
  electronics = simpar->ApplySPDBaselineAndNoise();
  UpdateMapNoise(ix,iz,electronics);
  //
  // Apply Threshold and write Digits.
  sig = GetMap()->GetSignalOnly(iz,ix);
  FillHistograms(ix,iz,sig+electronics);
  if(GetDebug(3)){
   cout<<sig<<"+"<<electronics<<">threshold("<<ix<<","<<iz
    <<")="<<GetThreshold() <<endl;
  } // end if GetDebug
  // if (sig+electronics <= GetThreshold()) continue;
  if (GetMap()->GetSignal(iz,ix) <= GetThreshold()) continue;
  dig.SetCoord1(iz);
  dig.SetCoord2(ix);
  dig.SetSignal(1);

  //        dig.SetSignalSPD((Int_t) GetMap()->GetSignal(iz,ix));
  Double_t aSignal =  GetMap()->GetSignal(iz,ix);
  if (TMath::Abs(aSignal)>2147483647.0) {
   //PH 2147483647 is the max. integer
   //PH This apparently is a problem which needs investigation
   AliWarning(Form("Too big or too small signal value %f",aSignal));
   aSignal = TMath::Sign((Double_t)2147483647,aSignal);
  }
  dig.SetSignalSPD((Int_t)aSignal);

  for(j=0;j<knmaxtrk;j++){
   if (j<GetMap()->GetNEntries()) {
    dig.SetTrack(j,GetMap()->GetTrack(iz,ix,j));
    dig.SetHit(j,GetMap()->GetHit(iz,ix,j));
   }else { // Default values
    dig.SetTrack(j,-3);
    dig.SetHit(j,-1);
   } // end if GetMap()
  } // end for j
  if(GetDebug(3)){
   cout<<iz<<","<<ix<<","<<*(GetMap()->GetpListItem(iz,ix))<<endl;
  } // end if GetDebug
  aliITS->AddSimDigit(0,&dig);
  // simulate fo signal response for this pixel hit:
  Float_t difference = TMath::Abs(GetMap()->GetSignal(iz,ix)-GetMap()->GetSignalFo(iz,ix));
  // only hit which released a good signal within the FO strobe are used to create the fastor
  if(GetMap()->GetSignalFo(iz,ix) >= GetThreshold()) fDetType->ProcessSPDDigitForFastOr(fModule, dig.GetCoord1(), dig.GetCoord2());
 } //  for ix/iz
}
//______________________________________________________________________
void AliITSsimulationSPD::CreateHistograms(){
 // create 1D histograms for tests
 // Inputs:
 //    none.
 // Outputs:
 //    none.
 // Return:
 //     none.

 if(GetDebug(1)) Info("CreateHistograms","create histograms");

 fHis = new TObjArray(GetNPixelsZ());
 fSPDname="spd_";
 for(Int_t i=0;i<GetNPixelsZ();i++) {
  Char_t pixelz[4];
  snprintf(pixelz,3,"%d",i);
  fSPDname.Append(pixelz);
  fHis->AddAt(new TH1F(fSPDname.Data(),"SPD maps",
     GetNPixelsX(),0.,(Double_t)GetNPixelsX()),i);
 } // end for i
}
//______________________________________________________________________
void AliITSsimulationSPD::FillHistograms(Int_t ix,Int_t iz,Double_t v){
 // Fill the histogram
 // Inputs:
 //    none.
 // Outputs:
 //    none.
 // Return:
 //     none.

 if(!GetHistArray()) return; // Only fill if setup.
 if(GetDebug(2)) Info("FillHistograms","fill histograms");
 GetHistogram(iz)->Fill(ix,v);
}
//______________________________________________________________________
void AliITSsimulationSPD::ResetHistograms(){
 // Reset histograms for this detector
 // Inputs:
 //    none.
 // Outputs:
 //    none.
 // Return:
 //     none.

 if(!GetHistArray()) return; // Only fill if setup.
 if(GetDebug(2)) Info("FillHistograms","fill histograms");
 for ( int i=0;i<GetNPixelsZ();i++ ) {
  if (fHis->At(i))    ((TH1F*)fHis->At(i))->Reset();
 } // end for i
}

//______________________________________________________________________
void AliITSsimulationSPD::SetCoupling(Int_t col, Int_t row) {
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
 //    Int_t col            z cell index
 //    Int_t row            x cell index
 // Outputs:
 //    none.
 // Return:
 //     none.
 Int_t j1,j2,flag=0;
 Double_t pulse1,pulse2;
 Double_t couplR=0.0,couplC=0.0;
 Double_t xr=0.;

 GetCouplings(couplC,couplR);
 if(GetDebug(3)) Info("SetCoupling","(col=%d,row=%d) "
   "Calling SetCoupling couplC=%e couplR=%e",
   col,row,couplC,couplR);
 j1 = col;
 j2 = row;
 pulse1 = GetMap()->GetSignalOnly(col,row);
 pulse2 = pulse1;
 for (Int_t isign=-1;isign<=1;isign+=2){// loop in col direction
  do{
   j1 += isign;
   xr = gRandom->Rndm();
   if ((j1<0) || (j1>GetNPixelsZ()-1) || (xr>couplC)){
    j1 = col;
    flag = 1;
   }else{
    UpdateMapNoise(row,j1,pulse1);
    //  flag = 0;
    flag = 1; // only first next!!
   } // end if
  } while(flag == 0);
  // loop in row direction
  do{
   j2 += isign;
   xr = gRandom->Rndm();
   if ((j2<0) || (j2>GetNPixelsX()-1) || (xr>couplR)){
    j2 = row;
    flag = 1;
   }else{
    UpdateMapNoise(j2,col,pulse2);
    //  flag = 0;
    flag = 1; // only first next!!
   } // end if
  } while(flag == 0);
 } // for isign
}
//______________________________________________________________________
void AliITSsimulationSPD::SetCouplingOld(Int_t col, Int_t row) {
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
 //    Int_t col            z cell index
 //    Int_t row            x cell index
 //    Int_t module         module number
 // Outputs:
 //    none.
 // Return:
 //     none.
 Int_t j1,j2,flag=0;
 Double_t pulse1,pulse2;
 Double_t couplR=0.0,couplC=0.0;

 GetCouplings(couplC,couplR);

 //  Debugging ...
 //    cout << "Threshold --> " << GetThreshold() << endl;  // dom
 //    cout << "Couplings --> " << couplC << " " << couplR << endl;  // dom

 if(GetDebug(3)) Info("SetCouplingOld","(col=%d,row=%d) "
   "Calling SetCoupling couplC=%e couplR=%e",
   col,row,couplC,couplR);
 for (Int_t isign=-1;isign<=1;isign+=2){// loop in col direction
  pulse1 = GetMap()->GetSignalOnly(col,row);
  pulse2 = pulse1;
  j1 = col;
  j2 = row;
  do{
   j1 += isign;
   pulse1 *= couplC;
   if ((j1<0)||(j1>GetNPixelsZ()-1)||(pulse1<GetThreshold())){
    pulse1 = GetMap()->GetSignalOnly(col,row);
    j1 = col;
    flag = 1;
   }else{
    UpdateMapNoise(row,j1,pulse1);
    // flag = 0;
    flag = 1;  // only first next !!
   } // end if
  } while(flag == 0);
  // loop in row direction
  do{
   j2 += isign;
   pulse2 *= couplR;
   if((j2<0)||(j2>(GetNPixelsX()-1))||(pulse2<GetThreshold())){
    pulse2 = GetMap()->GetSignalOnly(col,row);
    j2 = row;
    flag = 1;
   }else{
    UpdateMapNoise(j2,col,pulse2);
    // flag = 0;
    flag = 1; // only first next!!
   } // end if
  } while(flag == 0);
 } // for isign
}
//______________________________________________________________________
void AliITSsimulationSPD::GenerateStrobePhase()
{
 // Generate randomly the strobe
 // phase w.r.t to the LHC clock
 // Done once per event
 const Double_t kBunchLenght = 25e-9; // LHC clock
 fStrobePhase = ((Double_t)gRandom->Integer(fStrobeLenght))*kBunchLenght-
  (Double_t)fStrobeLenght*kBunchLenght+
  kBunchLenght/2;
}

