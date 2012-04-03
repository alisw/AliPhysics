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
#include "TF1.h"

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
  fNCells(0),
  fSuperModuleId(-999),
  fCellsAbsId(0x0),
  fCellsAmpFraction(0x0)

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
  fSuperModuleId(a.fSuperModuleId),
  fCellsAbsId(NULL),
  fCellsAmpFraction(NULL)

{
  //copy constructor
  for(Int_t i=0;i<kUnknownAli+1;i++) fCaloPID[i] = a.fCaloPID[i];
  SetCellsAbsId( a.fCellsAbsId );
  SetCellsAmplitudeFraction( a.fCellsAmpFraction );
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
    fSuperModuleId =  photon.fSuperModuleId;
    SetCellsAbsId( photon.fCellsAbsId );
    SetCellsAmplitudeFraction( photon.fCellsAmpFraction );

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
void  AliJPhoton::SetPID(const Double_t *pid) {
   //set pid
    if(pid){
       for(Int_t i=0; i<kUnknownAli+1; ++i) fCaloPID[i]=pid[i];
       SetProbPhot(fCaloPID[kPhotonAli]);
    }else{
      for(Int_t i=0; i<kUnknownAli+1; fCaloPID[i++]=0.){} 
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
    
  
   particleType pdg = kHadron;

  //Select most probable ID
  if(fCaloType == kPHOSCalo){
    if(fCaloPID[kPhotonAli] > wPhoton) pdg = kPhoton ;
    else if(fCaloPID[kPi0Ali] > wPi0) pdg = kPizero ; 
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
    if(fCaloPID[kPhotonAli] + fCaloPID[kElectronAli]  > wPhoton) pdg = kPhoton ;
    else if(fCaloPID[kPi0Ali] > wPi0) pdg = kPizero ; 
    //else if(chargedHadronWeight + neutralHadronWeight > wCharged) pdg = chargedHadron ;  
    //else if(neutralHadronWeight + chargedHadronWeight > wNeutral) pdg = neutralHadron ; 
    //else pdg =  neutralUnknown ;

  }
  
  delete wPhotonPHOSFormula;
  delete wPi0PHOSFormula;
  
  //if(fCaloType == kPHOSCalo) cout << "pdg " << pdg << endl;
  
  return pdg ;

}


