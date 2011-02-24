#ifndef ALIDIELECTRONPAIR_H
#define ALIDIELECTRONPAIR_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

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
  virtual Bool_t   XvYvZv(Double_t x[3]) const { x[0]=Xv(); x[1]=Yv(); x[2]=Zv(); return kTRUE; }
  
  virtual Double_t OneOverPt() const { return Pt()>0.?1./Pt():0.; }  //TODO: check
  virtual Double_t Phi()       const { return fPair.GetPhi();}
  virtual Double_t Theta()     const { return Pz()!=0?TMath::ATan(Pt()/Pz()):0.; } //TODO: check
  
  
  virtual Double_t E() const { return fPair.GetE();    }
  virtual Double_t M() const { return fPair.GetMass(); }
  
  virtual Double_t Eta() const { return fPair.GetEta();}
  virtual Double_t Y()  const  { 
    if((E()*E()-Px()*Px()-Py()*Py()-Pz()*Pz())>0.) return TLorentzVector(Px(),Py(),Pz(),E()).Rapidity();
    else return -1111.;
  }
  
  virtual Short_t Charge() const    { return fPair.GetQ();}
  virtual Int_t   GetLabel() const  { return fLabel;      }
  // PID
  virtual const Double_t *PID() const { return 0;} //TODO: check
  // Dummy
  Int_t PdgCode() const {return 0;}


  UChar_t GetType() const { return fType; }
  void SetType(Char_t type) { fType=type; }

  void SetLabel(Int_t label) {fLabel=label;}
  
  //inter leg information
  Double_t OpeningAngle()         const { return fD1.GetAngle(fD2);                             }
  Double_t DistanceDaughters()    const { return fD1.GetDistanceFromParticle(fD2);              }
  Double_t DistanceDaughtersXY()  const { return fD1.GetDistanceFromParticleXY(fD2);            }
  Double_t DeviationDaughters()   const { return fD1.GetDeviationFromParticle(fD2);             }
  Double_t DeviationDaughtersXY() const { return fD1.GetDeviationFromParticleXY(fD2);           }
  Double_t DeltaEta()             const { return TMath::Abs(fD1.GetEta()-fD2.GetEta());         }
  Double_t DeltaPhi()             const { Double_t dphi=TMath::Abs(fD1.GetPhi()-fD2.GetPhi());
                                          return (dphi>TMath::Pi())?dphi-TMath::Pi():dphi;      }
  // calculate cos(theta*) and phi* in HE and CS pictures
  void GetThetaPhiCM(Double_t &thetaHE, Double_t &phiHE, Double_t &thetaCS, Double_t &phiCS) const;
  
  Double_t ThetaPhiCM(Bool_t isHE, Bool_t isTheta) const;
  static Double_t ThetaPhiCM(const AliVParticle* d1, const AliVParticle* d2, 
			                       const Bool_t isHE, const Bool_t isTheta);
  // internal KF particle
  const AliKFParticle& GetKFParticle()       const { return fPair; }
  const AliKFParticle& GetKFFirstDaughter()  const { return fD1;   }
  const AliKFParticle& GetKFSecondDaughter() const { return fD2;   }
  
  // daughter references
  void SetRefFirstDaughter(AliVParticle * const track)  {fRefD1 = track;}
  void SetRefSecondDaughter(AliVParticle * const track) {fRefD2 = track;}
  
  AliVParticle* GetFirstDaughter()   const { return dynamic_cast<AliVParticle*>(fRefD1.GetObject()); }
  AliVParticle* GetSecondDaughter()  const { return dynamic_cast<AliVParticle*>(fRefD2.GetObject()); }

  
private:
  Char_t   fType;         // type of the pair e.g. like sign SE, unlike sign SE, ... see AliDielectron
  Int_t    fLabel;        // MC label
  
  AliKFParticle fPair;   // KF particle internally used for pair calculation
  AliKFParticle fD1;     // KF particle first daughter
  AliKFParticle fD2;     // KF particle1 second daughter
  
  TRef fRefD1;           // Reference to first daughter
  TRef fRefD2;           // Reference to second daughter
  
  ClassDef(AliDielectronPair,3)
};

#endif
