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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCExBEffectiveSector class                                                   //
// Correct for the rest of ExB effect which are not covered yet by physical models
//
// Motivation:
//   ExB correction: 
//      dr    =  c0* integral(Er/Ez) + c1* integral(Erphi/Ez)
//      drphi = -c1* integral(Er/Ez) + c0* integral(Erphi/Ez)
//   Where:      
//   wt = Bz*(k*vdrift/E)           ~ 0.3 at B=0.5 T 
//   c0 = 1/(1+T2*T2*wt*wt) 
//   c1 = T1*wt/(1+T1*T1*wt*wt)
//   
//  
//  3 correction maps 0 implemented as histogram used
//  R-Phi correction map obtained minimizing residuals betwee the track
//        and space points (AliTPCcalibAlign class). Track is defined using
//        the points from the refernce plain at the middle of the TPC
//        and vertex
//        Corrected primar tracks straight pointing to the primary vertex
//
//  R distortion - obtained using the cluster residuals in the setup with
//                 plus and minus field 
//                 Only high momenta tracks used for this calibration (1 GeV threshold)
//     drphi_plus-drphi_minus=-2*c1 integral(Er/Ez)
//               - Erphi/Ez cancels   
////////////////////////////////////////////////////////////////////////////
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"

#include "TMath.h"
#include "AliTPCROC.h"
#include "TFile.h"
#include "TAxis.h"
#include "TTree.h"
#include "TTreeStream.h"
#include "THnSparse.h"
#include "THnBase.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "AliTPCExBEffectiveSector.h"
ClassImp(AliTPCExBEffectiveSector)

AliTPCExBEffectiveSector::AliTPCExBEffectiveSector()
  : AliTPCCorrection("ExB_effectiveSector","ExB effective sector"),
    fC0(1.),fC1(0.), 
    fCorrectionR(0),        // radial correction
    fCorrectionRPhi(0),     // r-phi correction
    fCorrectionZ(0)        // z correction
{
  //
  // default constructor
  //
}

AliTPCExBEffectiveSector::~AliTPCExBEffectiveSector() {
  //
  // default destructor
  //
  delete fCorrectionR;        // radial correction
  delete fCorrectionRPhi;     // r-phi correction
  delete fCorrectionZ;        // z correction
}



void AliTPCExBEffectiveSector::Init() {
  //
  // Initialization funtion
  //
  
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);
}

void AliTPCExBEffectiveSector::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function 
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);
}



void AliTPCExBEffectiveSector::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction using the lookup table (histogram) of distortion
  // The histogram is created as poscl - postrack
  //   
  dx[0]=0;
  dx[1]=0;
  dx[2]=0;
  if (!fCorrectionRPhi) return;
  Double_t phi      = TMath::ATan2(x[1],x[0]);
  Double_t r        = TMath::Sqrt(x[1]*x[1]+x[0]*x[0]);
  Double_t sector   = 9.*phi/TMath::Pi();
  if (sector<0) sector+=18.;
  Double_t        kZ=x[2]/r;
  //
  if (kZ>1.2)        kZ= 1.2;
  if (kZ<-1.2)       kZ= -1.2;
  if (roc%36<18)  kZ= TMath::Abs(kZ);
  if (roc%36>=18) kZ=-TMath::Abs(kZ);
  if (TMath::Abs(kZ)<0.15){
    kZ = (roc%36<18) ? 0.15:-0.15;
  }  
  //
  Double_t dlR=0;
  Double_t dlRPhi=0;
  Double_t dlZ=0;
  Double_t rr=TMath::Max(r,fCorrectionRPhi->GetYaxis()->GetXmin()+0.01);
  rr=TMath::Min(rr,fCorrectionRPhi->GetYaxis()->GetXmax()-0.01);
  Double_t kZZ=TMath::Max(kZ,fCorrectionRPhi->GetZaxis()->GetXmin()+0.001);
  kZZ=TMath::Min(kZZ,fCorrectionRPhi->GetZaxis()->GetXmax()-0.001);

  if (fCorrectionRPhi) {  
    //    dlRPhi= -fCorrectionRPhi->Interpolate(sector,rr,kZZ);
    dlRPhi= -fCorrectionRPhi->GetBinContent(fCorrectionRPhi->FindBin(sector,rr,kZZ));
  }
  if (fCorrectionR)    {
    //    dlR= -fCorrectionR->Interpolate(sector,rr,kZZ);
    dlR= -fCorrectionR->GetBinContent(fCorrectionR->FindBin(sector,rr,kZZ));
  }
  if (fCorrectionZ)    {
    //    dlZ= -fCorrectionZ->Interpolate(sector,rr,kZZ);
    dlZ= -fCorrectionZ->GetBinContent(fCorrectionZ->FindBin(sector,rr,kZZ));
  }
  Double_t dr    = fC0*dlR  + fC1*dlRPhi;
  Double_t drphi = -fC1*dlR + fC0*dlRPhi;
   // Calculate distorted position
  if ( r > 0.0 ) {
    r   =  r   + dr;
    phi =  phi + drphi/r;
  }
  // Calculate correction in cartesian coordinates
  dx[0] = r * TMath::Cos(phi) - x[0];
  dx[1] = r * TMath::Sin(phi) - x[1];
  dx[2] = dlZ; 

}

void AliTPCExBEffectiveSector::Print(const Option_t* option) const {
  //
  // Print function to check the settings (e.g. the twist in the X direction)
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\t%s\n",GetName(),GetTitle());  
  if (opt.Contains("a")) { // Print all details
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C0: %1.4f, C1: %1.4f \n",fC0,fC1);
  }    
}

