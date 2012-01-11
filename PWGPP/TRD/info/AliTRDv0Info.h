#ifndef ALITRDV0INFO_H
#define ALITRDV0INFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDv0Info.h 34132 2009-08-06 11:18:32Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

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
class AliKFParticle;
class AliKFVertex;
class AliVEvent;
class AliVTrack;
//class AliTRDtrackV1;
class AliTRDtrackInfo;
class AliTRDv0Info : public TObject
{
public:
  enum ETRDv0Info{
    kNV0param = 10
    ,kNlayer   = AliTRDgeometry::kNlayer
    ,kNMomBins = 2//number of different momentum bins to consider for different cuts; first example: below /above 2.5 GeV-> to be refined!
    ,kNArmenteros = 2//number of Armenteros-Polanski parameters
  };

  enum EDaughter{  
    kNeg = 0
    ,kPos
    ,kNDaughters  //for positive and negative track
  };

  enum EDecayType{
    kGamma = 0
    ,kK0s
    ,kLambda
    ,kAntiLambda
    ,kNDecays     //number of decay types considered for reference data (conversions, K0s, Lambda, Anti-Lambda)
  };

  enum EDetector{
    kTPC = 0
    ,kTOF
    ,kITS
    ,kNDetectors  //TPC, TOF, ITS (TOF and ITS not implemented yet)
  };


  AliTRDv0Info();
  AliTRDv0Info(const AliTRDv0Info &ref);
  virtual ~AliTRDv0Info(){}
  

  Int_t GetQuality() const {return fQuality;}
  Float_t GetDCA() const {return fDCA;}
  Float_t GetPointingAngle() const {return fPointingAngle;}
  Float_t GetOpenAngle() const {return fOpenAngle;}
  Float_t GetPsiPair() const {return fPsiPair;}
  Float_t GetRadius() const {return fRadius;}
  Float_t GetV0Momentum() const {return fV0Momentum;}
  Double_t GetInvMass(Int_t iDecay) const {return fInvMass[iDecay];}
  Float_t GetDetPID(Int_t iDaughter, Int_t iDetector, Int_t iSpecies) const {return fDetPID[iDaughter][iDetector][iSpecies];}
  Float_t GetComPID(Int_t iDaughter, Int_t iSpecies) const {return fComPID[iDaughter][iSpecies];}
  Float_t GetTPCdEdx(Int_t iDaughter) const {return fTPCdEdx[iDaughter];}
  Float_t GetChi2ndf(Int_t decay) const {return fChi2ndf[decay];}
  
  //Get Cut values:
  
  Float_t GetUpDCA(Int_t iDecay) const {return fUpDCA[iDecay];}
  Float_t GetUpPointingAngle(Int_t iDecay) const {return fUpPointingAngle[iDecay];}
  Float_t GetUpOpenAngle(Int_t iDecay) const {return fUpOpenAngle[iDecay];}
  Float_t GetDownOpenAngle(Int_t iDecay) const {return fDownOpenAngle[iDecay];}
  Float_t GetUpPsiPair(Int_t iDecay) const {return fUpPsiPair[iDecay];}
  Float_t GetDownPsiPair(Int_t iDecay) const {return fDownPsiPair[iDecay];}
  Float_t GetUpRadius(Int_t iDecay) const {return fUpRadius[iDecay];}
  Float_t GetDownRadius(Int_t iDecay) const {return fDownRadius[iDecay];}
  Double_t GetUpInvMass(Int_t iDecay, Int_t iMomentum) const {return fUpInvMass[iDecay][iMomentum];}
  Double_t GetDownInvMass(Int_t iDecay) const {return fDownInvMass[iDecay];}
  Float_t GetDownTPCPIDneg(Int_t iPart) const {return fDownTPCPIDneg[iPart];}
  Float_t GetDownTPCPIDpos(Int_t iPart) const {return fDownTPCPIDpos[iPart];}
  Float_t GetDownComPIDneg(Int_t iPart) const {return fDownComPIDneg[iPart];}
  Float_t GetDownComPIDpos(Int_t iPart) const {return fDownComPIDpos[iPart];}
  Float_t GetDownComPIDnegPart(Int_t iPart) const {return fDownComPIDnegPart[iPart];}
  Float_t GetDownComPIDposPart(Int_t iPart) const {return fDownComPIDposPart[iPart];}
  
  
  Int_t  GetPID(Int_t ipart, AliTRDtrackInfo *track);
  Int_t  HasTrack(const AliTRDtrackInfo * const track) const;
  Int_t  HasTrack(Int_t ti) const;
  AliESDtrack *GetV0Daughter(Int_t sign); //Get positive of negative daughter of decay
  
  void   Print(Option_t *opt="") const;
  
  void   SetMagField(Float_t b) {fMagField = b;}
  void   SetV0tracks(AliESDtrack *p, AliESDtrack *n);
  void SetInputEvent(AliVEvent *e)      { fInputEvent = e; };
  void SetPrimaryVertex(AliKFVertex *v) { fPrimaryVertex = v; };
  
  //Set values of measured/calculated variables:
  void SetQuality(Int_t SQuality){fQuality = SQuality;}
  void SetDCA(Float_t SDCA){fDCA = SDCA;}
  void SetPointingAngle(Float_t SPointingAngle){fPointingAngle = SPointingAngle;}
  void SetOpenAngle(Float_t SOpenAngle){fOpenAngle = SOpenAngle;}
  void SetPsiPair(Float_t SPsiPair){fPsiPair = SPsiPair;}
  void SetRadius(Float_t SRadius){fRadius = SRadius;}
  void SetInvMass(Int_t iDecay, Float_t SInvMass){fInvMass[iDecay] = SInvMass;}
  void SetDetPID(Int_t iDaughter, Int_t iDetector, Int_t iSpecies, Float_t SDetPID){fDetPID[iDaughter][iDetector][iSpecies] = SDetPID;}
  void SetComPID(Int_t iDaughter, Int_t iSpecies, Float_t SComPID){fComPID[iDaughter][iSpecies] = SComPID;}
  void SetTPCdEdx(Int_t iDaughter, Float_t STpcdEdx){fTPCdEdx[iDaughter] = STpcdEdx;}
  void SetV0Momentum(Float_t SV0Momentum){fV0Momentum = SV0Momentum;}
  void SetChi2ndf(Int_t decay, Float_t SChi2ndf){fChi2ndf[decay]=SChi2ndf;}

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
  void SetDownTPCPIDneg(Int_t iPart, Double_t DownTPCPIDneg){fDownTPCPIDneg[iPart] = DownTPCPIDneg;}
  void SetDownTPCPIDpos(Int_t iPart, Double_t DownTPCPIDpos){fDownTPCPIDpos[iPart] = DownTPCPIDpos;}
  void SetDownComPIDneg(Int_t iPart, Double_t DownComPIDneg){fDownComPIDneg[iPart] = DownComPIDneg;}
  void SetDownComPIDpos(Int_t iPart, Double_t DownComPIDpos){fDownComPIDpos[iPart] = DownComPIDpos;}
  void SetDownComPIDnegPart(Int_t iPart, Double_t DownComPIDnegPart){fDownComPIDnegPart[iPart] = DownComPIDnegPart;}
  void SetDownComPIDposPart(Int_t iPart, Double_t DownComPIDposPart){fDownComPIDposPart[iPart] = DownComPIDposPart;}
 
  
  void SetV0Info(const AliESDv0 *v0);//gets most of the variables below



private:
  AliTRDv0Info& operator=(const AliTRDv0Info&);

  void GetDetectorPID();//operating with likelihood values of different detectors
  void CombinePID();//Bayesian combination of TPC and TOF likelihoods
  Bool_t TPCdEdxCuts(Int_t part, const AliTRDtrackInfo * const track);//direct cuts on TPC dE/dx

  Bool_t GetTPCdEdx();//TPC dE/dx values from both tracks
  Int_t Quality(const AliESDv0 * const esdv0);//checks for track/vertex quality criteria
  Double_t InvMass(Int_t part1, Int_t part2, const AliESDv0 *esdv0) const;//invariant mass of mother
  Float_t PsiPair(const AliESDv0 *esdv0);//angle between daughters in plane perpendicular to magnetic field (characteristically around zero for conversions)
  Float_t OpenAngle(const AliESDv0 *esdv0);//opening angle between V0 daughters; close to zero for conversions
  Float_t Radius(const AliESDv0 *esdv0);//distance of secondary to primary vertex in x-y-plane 
  Float_t DCA() const {return fDCA;}//distance of closest approach between supposed daughter tracks
  Float_t PointingAngle() const {return fPointingAngle;}//pointing angle: between vector from primary to secondary vertex and reconstructed momentum of V0 mother particle
  Float_t V0Momentum(const AliESDv0 *esdv0) const;//reconstructed momentum of V0 mother particle
  Bool_t V0SignCheck();//checks if daughters have opposite signs
  Bool_t  Armenteros(const AliESDv0 *esdv0, Int_t species);//the famous Armenteros-Polanski cut
  Double_t KFChi2ndf(Int_t part1, Int_t part2,Int_t decay);//Chi2ndf from KF
  AliKFParticle *CreateMotherParticle(const AliESDtrack *pdaughter, const AliESDtrack *ndaughter, Int_t pspec, Int_t nspec);//Mother Particle from KF
 

  //____________________________________________________________
  //Upper and lower limits for cut variables:

  Float_t fUpDCA[kNDecays];                  // DCA, upper limit
  Float_t fUpPointingAngle[kNDecays];        // pointing angle, upper limit
  Float_t fUpOpenAngle[kNDecays];            // opening angle, upper limit
  Float_t fDownOpenAngle[kNDecays];          // opening angle, lower limit
  Float_t fUpPsiPair[kNDecays];              // psi angle, upper limit
  Float_t fDownPsiPair[kNDecays];            // psi angle, lower limit
 
  Double_t fUpChi2ndf[kNDecays];             // upper Chi2/NDF limit
  Float_t fUpRadius[kNDecays];               // radius, upper limit
  Float_t fDownRadius[kNDecays];             // radius, lower limit
  Float_t fDownTPCPIDneg[AliPID::kSPECIES];  // TPC PID negatives, lower limit
  Float_t fDownTPCPIDpos[AliPID::kSPECIES];  // TPC PID positives, lower limit
  Float_t fDownComPIDneg[AliPID::kSPECIES];  // Combined PID negatives, lower limit
  Float_t fDownComPIDpos[AliPID::kSPECIES];  // Combined PID positives, lower limit
  Float_t fDownComPIDnegPart[AliPID::kSPECIES]; // Combined PID positive partner daughters (NOT the daughter track that would go into the reference data; here: pion daughters from Lambda decays; lower limit
  Float_t fDownComPIDposPart[AliPID::kSPECIES]; // Combined PID positive partner daughters (NOT the daughter track that would go into the reference data; here: pion daughters from Lambda decays; lower limit
  Double_t fUpInvMass[kNDecays][kNMomBins];  // invariant mass, upper limit
  Double_t fDownInvMass[kNDecays];           // invariant mass, lower limit
  
 
  Double_t fChi2ndf[kNDecays];//Chi2/NDF from KF
  
  
  Double_t fInvMass[kNDecays];  // invariant mass for different decay scenarios (conversions, K0s, Lambda->p+pi-, Lambda->p-pi+)
  Int_t fQuality;              // track quality status for both V0 daughters; OnFly, TPCrefit, Kinks, TPC clusters
  
  Double_t fDetPID[kNDaughters][kNDetectors][AliPID::kSPECIES]; // PID provided by TPC, TOF and ITS
  Double_t fComPID[kNDaughters][AliPID::kSPECIES];//Combined PID, momentarily from TPC and TOF only
  
  Float_t fDCA;  // Distance of closest approach of daughter tracks
  Float_t fPointingAngle;// = TMath::ACos(esdv0->GetV0CosineOfPointingAngle()); // Cosine of pointing angle
  Float_t fOpenAngle;  // opening angle between daughters
  Float_t fPsiPair; // /Angle between daughter momentum plane and plane perpendicular to magnetic field
  
  Bool_t fArmenteros[kNDecays];// Array for the Armenteros yes/no decision for all decays
  Double_t fMagField; //magnetic field strength
  Float_t fRadius; //distance of decay/conversion from primary vertex in x-y plane
  Float_t fV0Momentum; //V0 mother's momentum

  Int_t       fNindex; //indices of positive and negative daughter track
  Int_t       fPindex; //indices of positive and negative daughter track

  Float_t fTPCdEdx[kNDaughters]; //Energy deposition in the TPC
  
 
  AliVEvent            *fInputEvent;    // Input Event
  AliKFVertex          *fPrimaryVertex; // primary vertex
  
  AliESDtrack *fTrackP; //!positive daughter
  AliESDtrack *fTrackN; //!negative daughter
  
 
  
  ClassDef(AliTRDv0Info, 1) // extracted V0 MC information
};


#endif

