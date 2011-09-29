#ifndef ALIKFCONVERSIONPHOTON_H
#define ALIKFCONVERSIONPHOTON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class containing the aod information from conversions
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---

#include "TMath.h"
#include "AliConversionPhotonBase.h"
#include "AliKFParticle.h"
#include "AliV0Reader.h"
#include "AliESDEvent.h"
#include "AliPID.h"
#include "AliStack.h"
#include "TParticle.h"

class AliConversionPhotonBase;

class AliKFConversionPhoton : public AliKFParticle, public AliConversionPhotonBase {

 public: 

  //Constructors
  AliKFConversionPhoton();    
  AliKFConversionPhoton(AliKFParticle &kfparticle);
  AliKFConversionPhoton(AliV0Reader *fV0Reader);
  AliKFConversionPhoton(const AliKFParticle &fCurrentNegativeKFParticle,const AliKFParticle &fCurrentPositiveKFParticle);

  //Copy Constructor
  AliKFConversionPhoton(const AliKFConversionPhoton & g);           
  //assignment operator
  AliKFConversionPhoton & operator = (const AliKFConversionPhoton & g);

  //Destructor
  virtual ~AliKFConversionPhoton() {;}

  //
 void SetArmenterosQtAlpha(Double_t armenteros[2],const AliKFParticle &fCurrentNegativeKFParticle,const AliKFParticle &fCurrentPositiveKFParticle);
  void ConstructGamma(const AliKFParticle &fCurrentNegativeKFParticle,const AliKFParticle &fCurrentPositiveKFParticle);


  Double_t Phi() const;

  // GetInvariantMass

  Double_t M() const {return AliKFParticle::GetMass();}
  Double_t Pt() const {return AliKFParticle::GetPt();}
  Double_t P() const {return AliKFParticle::GetP();}
  Double_t Eta() const {return AliKFParticle::GetEta();}

  virtual Double_t GetPhotonMass() const {return AliKFParticle::GetMass();}
  virtual Double_t GetPhotonPt() const {return AliKFParticle::GetPt();}
  virtual Double_t GetPhotonP() const {return AliKFParticle::GetP();}
  virtual Double_t GetPhotonEta() const {return AliKFParticle::GetEta();}

  ClassDef(AliKFConversionPhoton,1)
};

#endif
