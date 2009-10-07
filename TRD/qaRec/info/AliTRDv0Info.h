#ifndef ALITRDV0INFO_H
#define ALITRDV0INFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#ifndef Root_TObject
#include "TObject.h"
#endif

#ifndef ALITRDGEOMETRY_H
#include "AliTRDgeometry.h"
#endif

#ifndef ALIPID_H
#include "AliPID.h"
#endif


class AliESDv0;
class AliESDtrack;
class AliESDEvent;
class AliTRDtrackV1;
class AliTRDtrackInfo;
class AliTRDv0Info : public TObject
{
public:
  enum ETRDv0Info{
    kNV0param = 10
    ,kNlayer   = AliTRDgeometry::kNlayer
    ,kNDetectors = 3//TPC, TOF, ITS (TOF and ITS not implemented yet)
    ,kNDaughters = 2//for positive and negative track
    ,kNDecays = 4//number of decay types considered for reference data (conversions, K0s, Lambda, Anti-Lambda)  
    ,kNMomBins = 2//number of different momentum bins to consider for different cuts; first example: below/above 2.5 GeV -> to be refined!
  };

  enum EDaughter{
    kNeg = 0
    ,kPos = 1
  };

  enum EDecayType{
    kGamma = 0
    ,kK0s = 1
    ,kLambda = 2
    ,kAntiLambda = 3
  };

  enum EDetector{
    kTPC = 0
    ,kTOF = 1
    ,kITS = 2
  };


  AliTRDv0Info();
  virtual ~AliTRDv0Info(){}

  Float_t Pplus[2*kNlayer];
  Float_t Pminus[2*kNlayer];

 
  void Print(Option_t *opt=0x0) const;
 
  Bool_t GetV0PID(Int_t ipart, AliTRDtrackInfo *track);//decides if a track is accepted for one of the reference samples!!

  //Set values of measured/calculated variables:
  void SetQuality(Int_t Quality){fQuality = Quality;}
  void SetPplus(Int_t iLayer, Float_t Pplus){fPplus[iLayer] = Pplus;}
  void SetPminus(Int_t iLayer, Float_t Pminus){fPminus[iLayer] = Pminus;}
  void SetDCA(Float_t DCA){fDCA = DCA;}
  void SetMomentum(Float_t Momentum){fMomentum = Momentum;}
  void SetPointingAngle(Float_t PointingAngle){fPointingAngle = PointingAngle;}
  void SetOpenAngle(Float_t OpenAngle){fOpenAngle = OpenAngle;}
  void SetPsiPair(Float_t PsiPair){fPsiPair = PsiPair;}
  void SetRadius(Float_t Radius){fRadius = Radius;}
  void SetInvMass(Int_t iDecay, Float_t InvMass){fInvMass[iDecay] = InvMass;}
  void SetDetPID(Int_t iDaughter, Int_t iDetector, Int_t iSpecies, Float_t DetPID){fDetPID[iDaughter][iDetector][iSpecies] = DetPID;}

//____________________________________________________________
 //Set cut values:

 void SetUpDCA(Int_t iDecay, Float_t UpDCA){fUpDCA[iDecay] = UpDCA;}
 void SetUpPointingAngle(Int_t iDecay, Float_t UpPointingAngle){fUpPointingAngle[iDecay] = UpPointingAngle;}
 void SetUpOpenAngle(Int_t iDecay, Float_t UpOpenAngle){fUpOpenAngle[iDecay] = UpOpenAngle;}
 void SetDownOpenAngle(Int_t iDecay, Float_t DownOpenAngle){fDownOpenAngle[iDecay] = DownOpenAngle;}
 void SetUpPsiPair(Int_t iDecay, Float_t UpPsiPair){fUpPsiPair[iDecay] = UpPsiPair;}
 void SetDownPsiPair(Int_t iDecay, Float_t DownPsiPair){fDownPsiPair[iDecay] = DownPsiPair;}
 void SetUpRadius(Int_t iDecay, Float_t UpRadius){fUpRadius[iDecay] = UpRadius;}
 void SetDownRadius(Int_t iDecay, Float_t DownRadius){fDownRadius[iDecay] = DownRadius;}
 void SetUpInvMass(Int_t iDecay, Int_t iMomentum, Double_t UpInvMass){fUpInvMass[iDecay][iMomentum] = UpInvMass;}
 void SetDownInvMass(Int_t iDecay, Double_t DownInvMass){fDownInvMass[iDecay] = DownInvMass;}
 void SetDownTPCPIDneg(Int_t iDecay, Double_t DownTPCPIDneg){fDownTPCPIDneg[iDecay] = DownTPCPIDneg;}
 void SetDownTPCPIDpos(Int_t iDecay, Double_t DownTPCPIDpos){fDownTPCPIDpos[iDecay] = DownTPCPIDpos;}

 

private:
  AliTRDv0Info(const AliTRDv0Info&);
  AliTRDv0Info& operator=(const AliTRDv0Info&);

 void GetESDv0Info(AliESDv0 *esdv0);//gets most of the variables below
  void GetDetectorPID();//operating with likelihood values of different detectors
  Int_t Quality(AliESDv0 *esdv0);//checks for track/vertex quality criteria
  Double_t InvMass(Int_t part1, Int_t part2, AliESDv0 *esdv0);//invariant mass of mother
  Float_t PsiPair(AliESDv0 *esdv0);//angle between daughters in plane perpendicular to magnetic field (characteristically around zero for conversions)
  Float_t OpenAngle(AliESDv0 *esdv0);//opening angle between V0 daughters; close to zero for conversions
  Float_t Radius(AliESDv0 *esdv0);//distance of secondary to primary vertex in x-y-plane 
  Float_t DCA() const {return fDCA;}//distance of closest approach between supposed daughter tracks
  Float_t PointingAngle() const {return fPointingAngle;}//pointing angle: between vector from primary to secondary vertex and reconstructed momentum of V0 mother particle
  Float_t V0Momentum(AliESDv0 *esdv0);//reconstructed momentum of V0 mother particle
  void V0fromTrack(AliTRDtrackInfo *track, Int_t ivertex);//checks if a track belongs to a vertex found by V0 finder
  
  AliESDEvent *fESD;


  Bool_t fHasV0; //Does this track belong to a vertex from a V0 finder?
 
  Int_t fQuality;              // track quality status for both V0 daughters; OnFly, TPCrefit, Kinks, TPC clusters
 
  Float_t fPplus[2*kNlayer];    // momentum and variance for the positive daughter  
  Float_t fPminus[2*kNlayer];   // momentum and variance for the negative daughter  
  Double_t fDetPID[kNDaughters][kNDetectors][AliPID::kSPECIES]; // PID provided by TPC, TOF and ITS

  Float_t fMomentum;  // Momentum of track at the vertex

  Float_t fDCA;  // Distance of closest approach of daughter tracks
  
  Float_t fPointingAngle;// = TMath::ACos(esdv0->GetV0CosineOfPointingAngle()); // Cosine of pointing angle
  
  Float_t fOpenAngle;  // opening angle between daughters
  
  Float_t fPsiPair; // /Angle between daughter momentum plane and plane perpendicular to magnetic field
  
  Double_t fInvMass[kNDecays];  // invariant mass for different decay scenarios (conversions, K0s, Lambda->p+pi-, Lambda->p-pi+)

  Double_t fMagField;

  Float_t fRadius; //distance of decay/conversion from primary vertex in x-y plane

  Int_t fTrackID;//track index


  Float_t fV0Momentum; //V0 mother's momentum

  //____________________________________________________________
  //Upper and lower limits for cut variables:

  Float_t fUpDCA[kNDecays];  
  
  Float_t fUpPointingAngle[kNDecays];
  
  Float_t fUpOpenAngle[kNDecays];

  Float_t fDownOpenAngle[kNDecays];
  
  Float_t fUpPsiPair[kNDecays];

  Float_t fDownPsiPair[kNDecays];
  
  Double_t fUpInvMass[kNDecays][kNMomBins]; 

  Double_t fDownInvMass[kNDecays]; 

  Float_t fUpRadius[kNDecays];

  Float_t fDownRadius[kNDecays];

  Float_t fDownTPCPIDneg[AliPID::kSPECIES];

  Float_t fDownTPCPIDpos[AliPID::kSPECIES];

 
  AliESDtrack *fTrackP; //positive daughter
  AliESDtrack *fTrackN; //negative daughter
  AliESDtrack *fTrack; //the current track in the ESDtrack loop (either positive or negative)


  Int_t fNindex;//indices of positive and negative daughter track
  Int_t fPindex;
  
  
  ClassDef(AliTRDv0Info, 0) // extracted V0 MC information
};


#endif

