#ifndef ALIVAODPARTICLE_H
#define ALIVAODPARTICLE_H
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

class AliTrackPoints;
class AliClusterMap;

class AliVAODParticle : public TObject {
public:
  AliVAODParticle() : TObject(){}
  virtual ~AliVAODParticle(){}

  AliVAODParticle(const AliVAODParticle& in);
  virtual AliVAODParticle& operator=(const AliVAODParticle& in); 
   

  // kinematics
  virtual Double_t         Mass() const = 0;
  virtual Double_t         E() const = 0;
  virtual Double_t         P() const = 0;
  virtual Double_t         Pt() const = 0;
  virtual Double_t         Px() const = 0;
  virtual Double_t         Py() const = 0;
  virtual Double_t         Pz() const = 0;
  virtual Double_t         Phi() const = 0;
  virtual Double_t         Theta() const = 0;
  virtual Double_t         Eta() const = 0;
  virtual Double_t         Y() const = 0;
  
  virtual void             SetMomentum(Double_t/*px*/,Double_t/*py*/,Double_t/*pz*/,Double_t/*E*/) = 0;
  virtual void             SetProductionVertex(Double_t /*vx*/, Double_t /*vy*/, Double_t /*vz*/, Double_t /*t*/) = 0;

  // PID
  virtual void             SetPdgCode(Int_t pdg, Float_t prob = 1.0) = 0;
  virtual Double_t         GetProbability(Int_t pdg) const = 0;
  virtual Double_t         GetPidProb() const = 0;//returns probability of being particle type defined by GetPdgCode() 
  virtual Int_t            GetMostProbable() const = 0;
  
  virtual Int_t            GetPdgCode() const = 0;//We need to assume some PID (f.e. energy calculation) 
                                                  //sometimes one track can apear in analysis twise (e.g. ones as pion ones as kaon)
  virtual Int_t            GetNumberOfPids() const = 0; //returns number of non zero PID probabilities
  virtual Int_t            GetNthPid         (Int_t idx) const = 0;//These two methods are made to be able to
  virtual Float_t          GetNthPidProb     (Int_t idx) const = 0;//copy pid information i.e. in copy ctors

  virtual Double_t         Charge() const = 0;
  
  // vertices
  virtual Double_t         Vx() const = 0;
  virtual Double_t         Vy() const = 0;
  virtual Double_t         Vz() const = 0;
  virtual Double_t         T()  const {return 0.0;};

  virtual AliVAODParticle*  Mother() const {return NULL;};
  virtual Bool_t           HasDecayVertex() const {return kFALSE;};
  virtual TVector3         DecayVertex() const {return TVector3();};
  virtual Int_t            NumberOfDaughters() const {return 0;};
  virtual AliVAODParticle* Daughter(Int_t /*index*/) const {return NULL;};

  virtual Int_t            GetUID() const { return 0;}//returns unique ID of this track 
                                                      //(may happen than the same track is selected
                                                      //twise, f.g. as a pion and as a kaon than both have the same UID)
  virtual void             SetUID(Int_t /*id*/){/* *this */}
  // type information
  virtual Bool_t           IsSimulated() {return kFALSE;};
  virtual Bool_t           IsTrack() {return kFALSE;};
  virtual Bool_t           IsCluster() {return kFALSE;};

  //HBT specific 
  virtual AliTrackPoints*  GetTPCTrackPoints() const {return 0x0;}
  virtual AliTrackPoints*  GetITSTrackPoints() const {return 0x0;}
  virtual AliClusterMap*   GetClusterMap() const {return 0x0;}
  virtual void             Print(const Option_t * opt = "") const = 0;

  static void    SetDebug(Int_t dbg=1){fgDebug=dbg;}
  static Int_t   GetDebug(){return fgDebug;}
  virtual void   Clear(Option_t * /*option*/ ="");
private:
  static Int_t fgDebug;//! debug level for all the analysis package

  ClassDef(AliVAODParticle,1)  // base class for AOD particles
};

#endif
