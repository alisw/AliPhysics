#ifndef ALIDIELECTRONPAIR_H
#define ALIDIELECTRONPAIR_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#                  AliDielectronPair                        #
//#               Class to store pair information             #
//#                                                           #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <TMath.h>
#include <TRef.h>
#include <TLorentzVector.h>

#include <AliKFParticle.h>
#include <AliVParticle.h>

class AliVTrack;

//TODO
//TODO: Should we inherit from AliVTrack in order to built another AliDielectronPair of KF with it?
//TODO
class AliDielectronPair : public AliVParticle {
public:
  AliDielectronPair();
  virtual ~AliDielectronPair();
  
  AliDielectronPair(AliVTrack * const particle1, Int_t pid1,
                    AliVTrack * const particle2, Int_t pid2, Char_t type);

  //TODO:  copy constructor + assignment operator
  
  void SetTracks(AliVTrack * const particle1, Int_t pid1,
                 AliVTrack * const particle2, Int_t pid2);

  //AliVParticle interface
  // kinematics
  virtual Double_t Px() const { return fPair.GetPx(); }
  virtual Double_t Py() const { return fPair.GetPy(); }
  virtual Double_t Pz() const { return fPair.GetPz(); }
  virtual Double_t Pt() const { return fPair.GetPt(); }
  virtual Double_t P() const  { return fPair.GetP();  }
  virtual Bool_t   PxPyPz(Double_t p[3]) const { p[0]=Px(); p[1]=Py(); p[2]=Pz(); return kTRUE; }
  
  virtual Double_t Xv() const { return fPair.GetX(); }
  virtual Double_t Yv() const { return fPair.GetY(); }
  virtual Double_t Zv() const { return fPair.GetZ(); }
  virtual Bool_t   XvYvZv(Double_t x[3]) const { x[0]=Xv(); x[1]=Xv(); x[2]=Zv(); return kTRUE; }
  
  virtual Double_t OneOverPt() const { return Pt()>0.?1./Pt():0.; }  //TODO: check
  virtual Double_t Phi()       const { return fPair.GetPhi();}
  virtual Double_t Theta()     const { return Pz()!=0?TMath::ATan(Pt()/Pz()):0.; } //TODO: check
  
  
  virtual Double_t E() const { return fPair.GetE();    }
  virtual Double_t M() const { return fPair.GetMass(); }
  
  virtual Double_t Eta() const { return fPair.GetEta();}
  virtual Double_t Y()  const  { return TLorentzVector(Px(),Py(),Pz(),E()).Rapidity();}
  
  virtual Short_t Charge() const    { return fPair.GetQ();}
  virtual Int_t   GetLabel() const  { return -1; }  //TODO: check
  // PID
  virtual const Double_t *PID() const { return 0;} //TODO: check

  Double_t OpeningAngle() const { return fOpeningAngle; }

  UChar_t GetType() const { return fType; }
  void SetType(Char_t type) { fType=type; }
  // internal KF particle
  const AliKFParticle& GetKFParticle() const { return fPair; }

  // daughter references
  AliVParticle* GetFirstDaughter()   const { return dynamic_cast<AliVParticle*>(fRefD1.GetObject()); }
  AliVParticle* GetSecondDaughter()  const { return dynamic_cast<AliVParticle*>(fRefD2.GetObject()); }
  // Dummy
  Int_t PdgCode() const {return 0;}
  
private:
  Double_t fOpeningAngle; // opening angle of the pair
  Char_t  fType;         // type of the pair e.g. like sign SE, unlike sign SE, ... see AliDielectron

  AliKFParticle fPair;   // KF particle internally used for pair calculation

  TRef fRefD1;           // Reference to first daughter
  TRef fRefD2;           // Reference to second daughter
  
  ClassDef(AliDielectronPair,1)
};

#endif
