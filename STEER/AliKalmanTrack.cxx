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

//-------------------------------------------------------------------------
//                Implementation of the AliKalmanTrack class
//   that is the base for AliTPCtrack, AliITStrackV2 and AliTRDtrack
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "AliKalmanTrack.h"
#include "AliPDG.h"
#include "TPDGCode.h"
#include "TDatabasePDG.h"

ClassImp(AliKalmanTrack)

Double_t AliKalmanTrack::fgConvConst;

//_______________________________________________________________________
AliKalmanTrack::AliKalmanTrack():
  fLab(-3141593),
  fChi2(0),
  fMass(0.13957),
  fN(0)
{
  //
  // Default constructor
  //
    if (fgConvConst==0) 
      Fatal("AliKalmanTrack()","The magnetic field has not been set !\n"); 
    
    fStartTimeIntegral = kFALSE;
    fIntegratedLength = 0;
    for(Int_t i=0; i<5; i++) fIntegratedTime[i] = 0;
}

//_______________________________________________________________________
AliKalmanTrack::AliKalmanTrack(const AliKalmanTrack &t):
  TObject(t),
  fLab(t.fLab),
  fChi2(t.fChi2),
  fMass(t.fMass),
  fN(t.fN)
{
  //
  // Copy constructor
  //
  if (fgConvConst==0) 
    Fatal("AliKalmanTrack(const AliKalmanTrack&)",
          "The magnetic field has not been set !\n"); 

  fStartTimeIntegral = t.fStartTimeIntegral;
  fIntegratedLength = t.fIntegratedLength;
  
  for (Int_t i=0; i<5; i++) 
    fIntegratedTime[i] = t.fIntegratedTime[i];
}
//_______________________________________________________________________
void AliKalmanTrack::StartTimeIntegral() 
{
  //
  // Start time integration
  // To be called at Vertex by ITS tracker
  //
  
  //if (fStartTimeIntegral) 
  //  Warning("StartTimeIntegral", "Reseting Recorded Time.");

  fStartTimeIntegral = kTRUE;
  for(Int_t i=0; i<fgkTypes; i++) fIntegratedTime[i] = 0;  
  fIntegratedLength = 0;
}
//_______________________________________________________________________
void AliKalmanTrack:: AddTimeStep(Double_t length) 
{
  // 
  // Add step to integrated time
  // this method should be called by a sublasses at the end
  // of the PropagateTo function or by a tracker
  // each time step is made.
  //
  // If integration not started function does nothing
  //
  // Formula
  // dt = dl * sqrt(p^2 + m^2) / p
  // p = pT * (1 + tg^2 (lambda) )
  //
  // pt = 1/external parameter [4]
  // tg lambda = external parameter [3]
  //
  //
  // Sylwester Radomski, GSI
  // S.Radomski@gsi.de
  // 
  
  static const Double_t kcc = 2.99792458e-2;

  if (!fStartTimeIntegral) return;
  
  fIntegratedLength += length;

  static Int_t pdgCode[fgkTypes]  = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
  TDatabasePDG *db = TDatabasePDG::Instance();

  Double_t xr, param[5];
  Double_t pt, tgl;
  
  GetExternalParameters(xr, param);
  pt =  1/param[4] ;
  tgl = param[3];

  Double_t p = TMath::Abs(pt * TMath::Sqrt(1+tgl*tgl));

  if (length > 100) return;

  for (Int_t i=0; i<fgkTypes; i++) {
    
    Double_t mass = db->GetParticle(pdgCode[i])->Mass();
    Double_t correction = TMath::Sqrt( pt*pt * (1 + tgl*tgl) + mass * mass ) / p;
    Double_t time = length * correction / kcc;

    //cout << mass << "\t" << pt << "\t" << p << "\t" 
    //     << correction << endl;

    fIntegratedTime[i] += time;
  }
}

//_______________________________________________________________________

Double_t AliKalmanTrack::GetIntegratedTime(Int_t pdg) const 
{
  //
  // Return integrated time hypothesis for a given particle
  // type assumption.
  //
  // Input parameter:
  // pdg - Pdg code of a particle type
  //


  if (!fStartTimeIntegral) {
    Warning("GetIntegratedTime","Time integration not started");
    return 0.;
  }

  static Int_t pdgCode[fgkTypes] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};

  for (Int_t i=0; i<fgkTypes; i++)
    if (pdgCode[i] == TMath::Abs(pdg)) return fIntegratedTime[i];

  Warning(":GetIntegratedTime","Particle type [%d] not found", pdg);
  return 0;
}
//_______________________________________________________________________

void AliKalmanTrack::PrintTime() const
{
  // For testing
  // Prints time for all hypothesis
  //

  static Int_t pdgCode[fgkTypes] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};

  for (Int_t i=0; i<fgkTypes; i++)
    printf("%d: %.2f  ", pdgCode[i], fIntegratedTime[i]);
  printf("\n");  
}

//_______________________________________________________________________

