#ifndef ALIAODPARTICLE_H
#define ALIAODPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// base class for AOD particles
//
/////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "AliAnalysis.h"

class AliTrackPoints;
class AliClusterMap;

class AliAODParticle : public TObject {
public:
   AliAODParticle(){}
   virtual ~AliAODParticle(){}
  // kinematics
  virtual TLorentzVector   FourMomentum() const = 0;
  virtual TVector3         Momentum() const {return FourMomentum().Vect();};
  virtual Double_t         Mass() const {return FourMomentum().M();};
  virtual Double_t         E() const {return FourMomentum().E();};
  virtual Double_t         P() const {return FourMomentum().P();};
  virtual Double_t         Pt() const {return FourMomentum().Pt();};
  virtual Double_t         Px() const {return FourMomentum().Px();};
  virtual Double_t         Py() const {return FourMomentum().Py();};
  virtual Double_t         Pz() const {return FourMomentum().Pz();};
  virtual Double_t         Phi() const {return FourMomentum().Phi();};
  virtual Double_t         Theta() const {return FourMomentum().Theta();};
  virtual Double_t         Eta() const {return FourMomentum().Eta();};
  virtual Double_t         Y() const {return FourMomentum().Rapidity();};

  // PID
  virtual Double_t         Charge() const = 0;
  virtual Double_t         GetProbability(Int_t pdg) const = 0;
  virtual Int_t            GetMostProbable() const = 0;
  
  virtual Int_t            GetPdgCode() const = 0;//We need to assume some PID (f.e. energy calculation) 
                                                  //sotimes one track can apear in analysis twise (e.g. ones as pion ones as kaon)
  
  // vertices
  virtual TVector3         ProductionVertex() const = 0;
  virtual Double_t         Vx() const {return ProductionVertex().X();};
  virtual Double_t         Vy() const {return ProductionVertex().Y();};
  virtual Double_t         Vz() const {return ProductionVertex().Z();};
  virtual AliAODParticle*  Mother() const {return NULL;};
  virtual Bool_t           HasDecayVertex() const {return kFALSE;};
  virtual TVector3         DecayVertex() const {return TVector3();};
  virtual Int_t            NumberOfDaughters() const {return 0;};
  virtual AliAODParticle*  Daughter(Int_t /*index*/) const {return NULL;};


  // type information
  virtual Bool_t           IsSimulated() {return kFALSE;};
  virtual Bool_t           IsTrack() {return kFALSE;};
  virtual Bool_t           IsCluster() {return kFALSE;};

  //HBT specific 
  virtual AliTrackPoints*  GetTrackPoints() const {return 0x0;}
  virtual AliClusterMap*   GetClusterMap() const {return 0x0;}
  virtual void             Print() const = 0;

  static void    SetDebug(Int_t dbg=1){fgDebug=dbg;}
  static Int_t   GetDebug(){return fgDebug;}

private:
  static Int_t fgDebug;//! debug level for all the analysis package

  ClassDef(AliAODParticle,1)  // base class for AOD particles
};

#endif
