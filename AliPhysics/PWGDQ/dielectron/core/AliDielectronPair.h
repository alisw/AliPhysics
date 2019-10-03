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
#include <TRandom3.h>

#include <AliKFParticle.h>
#include <AliVParticle.h>
#include <AliVEvent.h>

class AliVVertex;
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

  AliDielectronPair(const AliKFParticle * const particle1,
                    const AliKFParticle * const particle2,
                    AliVTrack * const refParticle1,
                    AliVTrack * const refParticle2,
                    Char_t type);

//TODO:  copy constructor + assignment operator

  void SetTracks(AliVTrack * const particle1, Int_t pid1,
                 AliVTrack * const particle2, Int_t pid2);

  void SetGammaTracks(AliVTrack * const particle1, Int_t pid1,
		      AliVTrack * const particle2, Int_t pid2);

  void SetTracks(const AliKFParticle * const particle1,
                 const AliKFParticle * const particle2,
                 AliVTrack * const refParticle1,
                 AliVTrack * const refParticle2);

  static void SetRandomizeDaughters(Bool_t random=kTRUE) { fRandomizeDaughters=random; }
  //static Bool_t GetRandomizeDaughters() { return fRandomizeDaughters; }

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

  //
  // Double_t GetLXY(const AliVVertex * const vtx) const;
  // Double_t GetPseudoProperTime(const AliVVertex * const vtx) const;


  UChar_t GetType() const { return fType; }
  void SetType(Char_t type) { fType=type; }
  static void SetBeamEnergy(AliVEvent *ev, Double_t beamEbyHand=-1.);

  // MC information
  void SetLabel(Int_t label) {fLabel=label;}
  void SetPdgCode(Int_t pdgCode) { fPdgCode=pdgCode; }
  Int_t PdgCode() const {return fPdgCode;}

  void SetProductionVertex(const AliKFParticle &Vtx) { fPair.SetProductionVertex(Vtx); }

  //inter leg information
  Double_t GetKFChi2()            const { return fPair.GetChi2();                               }
  Int_t    GetKFNdf()             const { return fPair.GetNDF();                                }
  Double_t OpeningAngle()         const { return fD1.GetAngle(fD2);                             }
  Double_t OpeningAngleXY()       const { return fD1.GetAngleXY(fD2);                           }
  Double_t OpeningAngleRZ()       const { return fD1.GetAngleRZ(fD2);                           }
  Double_t DistanceDaughters()    const { return fD1.GetDistanceFromParticle(fD2);              }
  Double_t DistanceDaughtersXY()  const { return fD1.GetDistanceFromParticleXY(fD2);            }
  Double_t DeviationDaughters()   const { return fD1.GetDeviationFromParticle(fD2);             }
  Double_t DeviationDaughtersXY() const { return fD1.GetDeviationFromParticleXY(fD2);           }
  Double_t DeltaEta()             const { return TMath::Abs(fD1.GetEta()-fD2.GetEta());         }
//   Double_t DeltaPhi()             const { Double_t dphi=TMath::Abs(fD1.GetPhi()-fD2.GetPhi());
//                                           return (dphi>TMath::Pi())?dphi-TMath::Pi():dphi;      }
  Double_t DeltaPhi()             const { return fD1.GetAngleXY(fD2);     }
  Double_t DeltaCotTheta()        const;

  // calculate cos(theta*) and phi* in HE and CS pictures
  void GetThetaPhiCM(Double_t &thetaHE, Double_t &phiHE, Double_t &thetaCS, Double_t &phiCS) const;


  Double_t ThetaPhiCM(Bool_t isHE, Bool_t isTheta) const;
  static Double_t ThetaPhiCM(const AliVParticle* d1, const AliVParticle* d2,
			                       Bool_t isHE, Bool_t isTheta);

  Double_t PsiPair(Double_t MagField)const; //Angle cut w.r.t. to magnetic field
  Double_t PhivPair(Double_t MagField)const; //Angle of ee plane w.r.t. to magnetic field

  //Calculate the angle between ee decay plane and variables
  Double_t GetPairPlaneAngle(Double_t kv0CrpH2, Int_t VariNum) const;

  Double_t GetCosPointingAngle(const AliVVertex *primVtx) const;
  Double_t GetArmAlpha() const;
  Double_t GetArmPt()    const;
  void GetDCA(const AliVVertex *primVtx, Double_t d0z0[2]) const;

  // Calculate inner product of strong magnetic field and ee plane
  Double_t PairPlaneMagInnerProduct(Double_t ZDCrpH1) const;


  // internal KF particle
  const AliKFParticle& GetKFParticle()       const { return fPair; }
  const AliKFParticle& GetKFFirstDaughter()  const { return fD1;   }
  const AliKFParticle& GetKFSecondDaughter() const { return fD2;   }

  // daughter references
  void SetRefFirstDaughter(AliVParticle * const track)  {fRefD1 = track;}
  void SetRefSecondDaughter(AliVParticle * const track) {fRefD2 = track;}

  AliVParticle* GetFirstDaughterP()   const { return dynamic_cast<AliVParticle*>(fRefD1.GetObject()); }
  AliVParticle* GetSecondDaughterP()  const { return dynamic_cast<AliVParticle*>(fRefD2.GetObject()); }

  void SetKFUsage(Bool_t KFUsage) {fKFUsage = KFUsage;}
  Bool_t GetKFUsage() const {return fKFUsage;}



private:
  Char_t   fType;         // type of the pair e.g. like sign SE, unlike sign SE, ... see AliDielectron
  Int_t    fLabel;        // MC label
  Int_t    fPdgCode;      // pdg code in case it is a MC particle
  static Double_t fBeamEnergy; //!beam energy

  AliKFParticle fPair;   // KF particle internally used for pair calculation
  AliKFParticle fD1;     // KF particle first daughter
  AliKFParticle fD2;     // KF particle1 second daughter

  TRef fRefD1;           // Reference to first daughter
  TRef fRefD2;           // Reference to second daughter

  Bool_t fKFUsage;       // Use KF for vertexing

  static Bool_t   fRandomizeDaughters;
  static TRandom3 fRandom3;

  ClassDef(AliDielectronPair,5)
};

#endif
