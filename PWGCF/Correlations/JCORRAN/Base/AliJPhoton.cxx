/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

// $Id: AliJPhoton.cxx,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliJPhoton.cxx
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.4 $
  \date $Date: 2008/05/08 13:44:45 $
  */
////////////////////////////////////////////////////

#include "AliJBaseTrack.h"
#include "AliJPhoton.h"
#include <TF1.h>

ClassImp(AliJPhoton);
//______________________________________________________________________________
AliJPhoton::AliJPhoton() : 
    AliJBaseTrack(),
    fChi2(-999),              
    fTof(-999),                   
    fX(-999),               
    fY(-999),                  
    fZ(-999),
    fProbPhot(-999), 
    fTrackDx(-1), 
    fTrackDz(-1), 
  fCaloType(kUndef),
  fDistToBadChannel(-999),
  fDispersion(-999),
  fM20(-999),
  fM02(-999),
  fEmcCpvDist(-999),
  fNCells(-999),
  fNTracksMatched(-999),
  fSuperModuleId(-999),
  fCellsAbsId(0x0),
  fCellsAmpFraction(0x0),
  fNEMCLabel(0),
  fEMCLabel(0x0),
  fEMax(-999),
  fECross(-999),
  fECore(-999),
  fNLM(-999),
  fCellsIndex(0x0)
//  fCellsAmp(0x0)

{
  // default constructor
  for(Int_t i=0;i<kUnknownAli+1;i++) fCaloPID[i]=-1;

  SetPID((Double_t*)NULL);
}

//_____________________________________________________________________________
AliJPhoton::AliJPhoton(const AliJPhoton& a) : 
    AliJBaseTrack(a),
    fChi2(a.fChi2),
    fTof(a.fTof),
    fX(a.fX),
    fY(a.fY),
    fZ(a.fZ),
    fProbPhot(a.fProbPhot), 
    fTrackDx(a.fTrackDx), 
    fTrackDz(a.fTrackDz), 
  fCaloType(a.fCaloType),
  fDistToBadChannel(a.fDistToBadChannel),
  fDispersion(a.fDispersion),
  fM20(a.fM20),
  fM02(a.fM02),
  fEmcCpvDist(a.fEmcCpvDist),
  fNCells(a.fNCells),
  fNTracksMatched(a.fNTracksMatched),
  fSuperModuleId(a.fSuperModuleId),
  fCellsAbsId(NULL),
  fCellsAmpFraction(NULL),
  fNEMCLabel(a.fNEMCLabel),
  fEMCLabel(NULL),
  fEMax(a.fEMax),
  fECross(a.fECross),
  fECore(a.fECore),
  fNLM(a.fNLM),
  fCellsIndex(NULL)
//  fCellsAmp(NULL)

{
  //copy constructor
  for(Int_t i=0;i<kUnknownAli+1;i++) fCaloPID[i] = a.fCaloPID[i];
  SetCellsAbsId( a.fCellsAbsId );
  SetCellsAmplitudeFraction( a.fCellsAmpFraction );
  SetCellsIndex( a.fCellsIndex );
	SetEMCLabel( a.fEMCLabel );
//  SetCellsAmplitude( a.fCellsAmp );
}


//_____________________________________________________________________________
AliJPhoton::~AliJPhoton(){
  // destructor
  Clear("");
}

//_____________________________________________________________________________
void AliJPhoton::Clear( Option_t * /* option */ ){
  // cleaner
    if(fCellsAbsId)       delete [] fCellsAbsId;
    if(fCellsAmpFraction) delete [] fCellsAmpFraction;
    if( fEMCLabel )       delete [] fEMCLabel;
    if( fCellsIndex )     delete [] fCellsIndex;
    fCellsAbsId = 0;
    fCellsAmpFraction = 0;
		fEMCLabel = 0;
    fCellsIndex = 0;
//    if(fCellsAmp)         delete [] fCellsAmp;
}

//_____________________________________________________________________________
AliJPhoton& AliJPhoton::operator=(const AliJPhoton& photon){
  //operator=    
  if(this != &photon){
    AliJBaseTrack::operator=(photon);
    fChi2 = photon.fChi2;
    fTof  = photon.fTof;
    fX    = photon.fX;
    fY    = photon.fY;
    fZ    = photon.fZ;
    fProbPhot = photon.fProbPhot;
    fTrackDx  = photon.fTrackDx;
    fTrackDz  = photon.fTrackDz;
    fCaloType = photon.fCaloType;
    for(Int_t i=0; i<kUnknownAli+1; i++){
      fCaloPID[i] = photon.fCaloPID[i];
    }
    fDistToBadChannel = photon.fDistToBadChannel;
    fDispersion    = photon.fDispersion;
    fM20           = photon.fM20;
    fM02           = photon.fM02;
    fEmcCpvDist    = photon.fEmcCpvDist;
    fNCells        = photon.fNCells;
    fNTracksMatched = photon.fNTracksMatched;
    fSuperModuleId =  photon.fSuperModuleId;
    SetCellsAbsId( photon.fCellsAbsId );
    SetCellsAmplitudeFraction( photon.fCellsAmpFraction );
    SetEMCLabel( photon.fEMCLabel );
//    SetCellsAmplitude( photon.fCellsAmp );
    fEMax = photon.fEMax;
    fECross = photon.fECross;
    fECore = photon.fECore;
    fNLM  = photon.fNLM;
    SetCellsIndex( photon.fCellsIndex );

  }

  return *this;
}
//______________________________________________________________________________
void  AliJPhoton::SetCellsAbsId(const UShort_t *array)
{
    //  Set the array of cell absId numbers 
    if (fNCells) {
        if(fCellsAbsId){ delete [] fCellsAbsId; fCellsAbsId = NULL; }
  fCellsAbsId = new  UShort_t[fNCells];
  for (Int_t i = 0; i < fNCells; i++) fCellsAbsId[i] = array[i];
    }
}

//_______________________________________________________________________
void  AliJPhoton::SetCellsAmplitudeFraction(const Double32_t *array)
{
    //  Set the array of cell amplitude fraction
    if (fNCells) {
  if(fCellsAmpFraction){ delete [] fCellsAmpFraction; fCellsAmpFraction = NULL;}
  fCellsAmpFraction = new  Double32_t[fNCells];
  for (Int_t i = 0; i < fNCells; i++) fCellsAmpFraction[i] = array[i];
    }
}

//______________________________________________________________________________
void  AliJPhoton::SetEMCLabel(const Int_t *array)
{
    //  Set the array of cell absId numbers 
    if (fNEMCLabel) {
        if(fEMCLabel){ delete [] fEMCLabel; fEMCLabel = NULL; }
  fEMCLabel = new  Int_t[fNEMCLabel];
  for (Int_t i = 0; i < fNEMCLabel; i++) fEMCLabel[i] = array[i];
    }
}

//______________________________________________________________________________
void  AliJPhoton::SetCellsIndex(const Int_t *array)
{
  //  Set the array of cell absId numbers 
  if( !array ){
    if( fCellsIndex ){
      delete [] fCellsIndex;
      fCellsIndex = NULL;
    }
  }
  else if (fNCells) {
    if( fCellsIndex ){
      delete [] fCellsIndex;
      fCellsIndex = NULL;
    }
    
    fCellsIndex = new  Int_t[fNCells];
    for (Int_t i = 0; i < fNCells; i++)
      fCellsIndex[i] = array[i];
  }
}

//______________________________________________________________________________
void  AliJPhoton::SetCellIndex(const Int_t pos, const Int_t ind)
{
  //  Set the array of cell absId numbers 

  Int_t i;
  if (fNCells) {
    if( ! fCellsIndex ){
      fCellsIndex = new  Int_t[fNCells];
      for( i = 0; i < fNCells; i++ )
        fCellsIndex[i] = -1;
    }
  }
  
  fCellsIndex[pos] = ind;
}

//______________________________________________________________________________
void  AliJPhoton::ClearCellsIndex()
{
  //  Clears the array of cell absId numbers 
  if(fCellsIndex){
    delete [] fCellsIndex;
    fCellsIndex = NULL;
  }
}

// //_______________________________________________________________________
// void  AliJPhoton::SetCellsAmplitude(const Double32_t *array)
// {
//     //  Set the array of cell amplitude fraction
//     if (fNCells) {
//   if(fCellsAmp){ delete [] fCellsAmp; fCellsAmp = NULL;}
//   fCellsAmp = new  Double_t[fNCells];
//   for (Int_t i = 0; i < fNCells; i++) fCellsAmp[i] = array[i];
//     }
// }

//______________________________________________________________________________
void  AliJPhoton::SetPID(const Double32_t *pid) {
   //set pid
    if(pid){
       for(Int_t i=0; i<kUnknownAli; ++i) fCaloPID[i]=pid[i];
       SetProbPhot(fCaloPID[kPhotonAli]);
    }else{
      for(Int_t i=0; i<kUnknownAli; fCaloPID[i++]=0.){} 
      fCaloPID[kUnknownAli]=1.;
    }
}


//______________________________________________________________________________
particleType AliJPhoton::GetParticleType() {

    // return the most problable particle type
    // Note: following the AliCaloPID implementation 
    
    //Init default weights 
    Float_t wPhoton = 0.75 ;
    Float_t wPi0 = 0.8 ;
    Float_t wElectron = 0.5 ;
    Float_t wCharged = 0.5 ;
    Float_t wNeutral = 0.5 ;

    Bool_t usePHOSweightFormula = kTRUE;
    //Formula to set the PID weight threshold for photon or pizero
    // ALICE-INT-2005-016 (2007)
    TFormula* wPhotonPHOSFormula = 
    new TFormula("photonWeight","0.75*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
    TFormula* wPi0PHOSFormula = 
    new TFormula("pi0Weight","0.80*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");

    if(fCaloType == kPHOSCalo && usePHOSweightFormula){
       wPhoton  = wPhotonPHOSFormula->Eval(E()) ;
       wPi0 = wPi0PHOSFormula->Eval(E());
      // cout << "wPhotonF " << wPhoton << " wphoton " << fCaloPID[kPhotonAli] << " wPi0F " << wPi0 << " wpi0 " << fCaloPID[kPi0Ali] << endl;
    }
    
    if(fCaloType == kEMCALCalo){
       wPhoton   =  0.8 ;
       wPi0      =  0.5 ;
       wElectron =  0.8 ;
       wCharged  =  0.5 ;
       wNeutral  =  0.5 ;
    }
    
  
   particleType pdg = kJHadron;

  //Select most probable ID
  if(fCaloType == kPHOSCalo){
    if(fCaloPID[kPhotonAli] > wPhoton) pdg = kJPhoton ;
    else if(fCaloPID[kPi0Ali] > wPi0) pdg = kJPizero ; 
    //else if(fCaloPID[kElectronAli] > wElectron)  pdg = electron ;
    //else if(fCaloPID[kEleConAli] >  wElectron) pdg = electronCon ;
    //else if(chargedHadronWeight > wCharged) pdg = chargedHadron ;  
    //else if(neutralHadronWeight > wNeutral) pdg = neutralHadron ; 
    //else if(allChargedWeight >  allNeutralWeight)
    // pdg = chargedUnknown ; 
    //else 
    //  pdg = neutralUnknown ;
  }
  else{//EMCAL
    //Temporal solution, electrons and photons not differenciated
    if(fCaloPID[kPhotonAli] + fCaloPID[kElectronAli]  > wPhoton) pdg = kJPhoton ;
    else if(fCaloPID[kPi0Ali] > wPi0) pdg = kJPizero ; 
    //else if(chargedHadronWeight + neutralHadronWeight > wCharged) pdg = chargedHadron ;  
    //else if(neutralHadronWeight + chargedHadronWeight > wNeutral) pdg = neutralHadron ; 
    //else pdg =  neutralUnknown ;

  }
  
  delete wPhotonPHOSFormula;
  delete wPi0PHOSFormula;
  
  //if(fCaloType == kPHOSCalo) cout << "pdg " << pdg << endl;
  
  return pdg ;

}


