/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliTRDv0Info.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Gathers all information necessary for reference data selection about  //
//  the track and (in case) its corresponding V0.                         //
//  Carries out the selection of electrons (from gamma conversions),      //
//  pions (from K0s decays) and protons (from Lambda and Anti-Lambda      //
//  decays) by cuts specific for the respective decay and particle        //
//  species.                                                              //
//  (M.Heide, 2009/10/06)                                                 //
//                                                                        //
//  Authors:                                                              //
//   Alex Bercuci <A.Bercuci@gsi.de>                                      //
//   Alex Wilk    <wilka@uni-muenster.de>                                 //
//   Markus Heide <mheide@uni-muenster.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TDatabasePDG.h"

#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliLog.h"
#include "TVector3.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliTRDv0Info.h"
#include "AliTRDtrackInfo.h"
#include "AliTRDtrackInfo.h"

ClassImp(AliTRDv0Info)

//_________________________________________________
AliTRDv0Info::AliTRDv0Info()
  : TObject()
  ,fQuality(0)
  ,fDCA(10)
  ,fPointingAngle(10)
  ,fOpenAngle(10)
  ,fPsiPair(99)
  ,fMagField(0)
  ,fRadius(0)
  ,fV0Momentum(0)
  ,fTrackP(NULL)
  ,fTrackN(NULL)
  ,fNindex(0)
  ,fPindex(0)
  ,fInputEvent(NULL)
  ,fPrimaryVertex(NULL)
{
  //
  // Default constructor
  //

  memset(fPplus, 0, 2*kNlayer*sizeof(Float_t));
  memset(fPminus, 0, 2*kNlayer*sizeof(Float_t));
  memset(fDetPID, 0, 2*kNDaughters*kNDetectors*AliPID::kSPECIES*sizeof(Float_t));
  memset(fComPID, 0, 2*kNDaughters*AliPID::kSPECIES*sizeof(Float_t));
  memset(fInvMass, 0, kNMomBins*kNDecays*sizeof(Double_t));
  memset(fArmenteros, 0, kNDecays*sizeof(Bool_t));
  memset(fTPCdEdx, 0, kNDecays*sizeof(Float_t));
  memset(fChi2ndf, 0, kNDecays*sizeof(Double_t));

  /////////////////////////////////////////////////////////////////////////////
  //Set Cut values: First specify decay in brackets, then the actual cut value!
  ///////////////////////////////////////////////////////////////////////////// 

  //Upper limit for distance of closest approach of two daughter tracks :
  fUpDCA[kGamma] = 1000.;
  fUpDCA[kK0s] = 0.08;
  fUpDCA[kLambda] = 0.2;
  fUpDCA[kAntiLambda] = 0.2;

  //Upper limit for pointing angle (= angle between between vector from primary to secondary vertex and reconstructed momentum of V0 mother particle) :
  fUpPointingAngle[kGamma] = 0.03;
  fUpPointingAngle[kK0s] = 0.03;
  fUpPointingAngle[kLambda] = 0.04;
  fUpPointingAngle[kAntiLambda] = 0.04;

  //Upper limit for invariant mass of V0 mother :
  fUpInvMass[kGamma][0] = 0.05;// second pair of brackets is for momentum bin: 0: below mother momentm of 2.5 GeV
  fUpInvMass[kGamma][1] = 0.07;//1: above 2.5 GeV
  fUpInvMass[kK0s][0] = fUpInvMass[kK0s][1] = 0.50265;
  fUpInvMass[kLambda][0] = fUpInvMass[kLambda][1] = 1.1207;
  fUpInvMass[kAntiLambda][0] = fUpInvMass[kAntiLambda][1] = 1.1207;

  //Lower limit for invariant mass of V0 mother :
  fDownInvMass[kGamma] = -1.;
  fDownInvMass[kK0s] = 0.49265;
  fDownInvMass[kLambda] = 1.107;
  fDownInvMass[kAntiLambda] = 1.107;

  //Upper limit for KF Chi2/NDF value;
  fUpChi2ndf[kGamma] = 10000.;//7.;
  fUpChi2ndf[kK0s] = 10000.;//5.;
  fUpChi2ndf[kLambda] = 10000.;//5.;
  fUpChi2ndf[kAntiLambda] = 10000.;//5.;

  //Lower limit for distance from secondary vertex to primary vertex in x-y plane :
  fDownRadius[kGamma] = 6.;
  fDownRadius[kK0s] = 0.;
  fDownRadius[kLambda] = 0.;
  fDownRadius[kAntiLambda] = 0.;

  //Upper limit for distance from secondary vertex to primary vertex in x-y plane :
  fUpRadius[kGamma] = 1000.;
  fUpRadius[kK0s] = 20.;
  fUpRadius[kLambda] = 1000.;
  fUpRadius[kAntiLambda] = 1000.;

  //Upper limit for opening angle between two daughter tracks (characteristically near zero for conversions) :
  fUpOpenAngle[kGamma] = 0.1;
  fUpOpenAngle[kK0s] = 3.15;
  fUpOpenAngle[kLambda] = 3.15;
  fUpOpenAngle[kAntiLambda] = 3.15;

  //Upper limit for angle between daughter momentum plane and plane perpendicular to magnetic field (characteristically around zero for conversions) :
  fUpPsiPair[kGamma] = 0.05;
  fUpPsiPair[kK0s] = 1.6;
  fUpPsiPair[kLambda] = 1.6;
  fUpPsiPair[kAntiLambda] = 1.6;

  //Lower limit for likelihood value of TPC PID :
  fDownTPCPIDneg[AliPID::kElectron] = 0.;
  fDownTPCPIDpos[AliPID::kElectron] = 0.;

  fDownTPCPIDneg[AliPID::kMuon] = 0.;
  fDownTPCPIDpos[AliPID::kMuon] = 0.;

  fDownTPCPIDneg[AliPID::kPion] = 0.;
  fDownTPCPIDpos[AliPID::kPion] = 0.;

  fDownTPCPIDneg[AliPID::kKaon] = 0.;
  fDownTPCPIDpos[AliPID::kKaon] = 0.;

  fDownTPCPIDneg[AliPID::kProton] = 0.;
  fDownTPCPIDpos[AliPID::kProton] = 0.;

 //Lower limit for likelihood value of combined PID :
  fDownComPIDneg[AliPID::kElectron] = 0.;
  fDownComPIDpos[AliPID::kElectron] = 0.;

  fDownComPIDneg[AliPID::kMuon] = 0.;
  fDownComPIDpos[AliPID::kMuon] = 0.;

  fDownComPIDneg[AliPID::kPion] = 0.;
  fDownComPIDpos[AliPID::kPion] = 0.;

  fDownComPIDneg[AliPID::kKaon] = 0.;
  fDownComPIDpos[AliPID::kKaon] = 0.;

  fDownComPIDneg[AliPID::kProton] = 0.;
  fDownComPIDpos[AliPID::kProton] = 0.;

 //Lower limit for likelihood value of combined PID for daughter track which doesn't enter reference data (here: pion daughters from Lambda decays:
  fDownComPIDnegPart[AliPID::kElectron] = 0.;
  fDownComPIDposPart[AliPID::kElectron] = 0.;

  fDownComPIDnegPart[AliPID::kMuon] = 0.;
  fDownComPIDposPart[AliPID::kMuon] = 0.;

  fDownComPIDnegPart[AliPID::kPion] = 0.;
  fDownComPIDposPart[AliPID::kPion] = 0.;

  fDownComPIDnegPart[AliPID::kKaon] = 0.;
  fDownComPIDposPart[AliPID::kKaon] = 0.;

  fDownComPIDnegPart[AliPID::kProton] = 0.;
  fDownComPIDposPart[AliPID::kProton] = 0.;

  //Parameters for data with well-calibrated PID (after usage of tender):
  /* //Lower limit for likelihood value of TPC PID :
  fDownTPCPIDneg[AliPID::kElectron] = 0.21;
  fDownTPCPIDpos[AliPID::kElectron] = 0.21;

  fDownTPCPIDneg[AliPID::kMuon] = 0.21;
  fDownTPCPIDpos[AliPID::kMuon] = 0.21;

  fDownTPCPIDneg[AliPID::kPion] = 0.21;
  fDownTPCPIDpos[AliPID::kPion] = 0.21;

  fDownTPCPIDneg[AliPID::kKaon] = 0.21;
  fDownTPCPIDpos[AliPID::kKaon] = 0.21;

  fDownTPCPIDneg[AliPID::kProton] = 0.21;
  fDownTPCPIDpos[AliPID::kProton] = 0.21;

  //Lower limit for likelihood value of combined PID :
  fDownComPIDneg[AliPID::kElectron] = 0.21;
  fDownComPIDpos[AliPID::kElectron] = 0.21;

  fDownComPIDneg[AliPID::kMuon] = 0.21;
  fDownComPIDpos[AliPID::kMuon] = 0.21;

  fDownComPIDneg[AliPID::kPion] = 0.9;
  fDownComPIDpos[AliPID::kPion] = 0.9;

  fDownComPIDneg[AliPID::kKaon] = 0.21;
  fDownComPIDpos[AliPID::kKaon] = 0.21;

  fDownComPIDneg[AliPID::kProton] = 0.9;
  fDownComPIDpos[AliPID::kProton] = 0.9;

 //Lower limit for likelihood value of combined PID for daughter track which doesn't enter reference data (here: pion daughters from Lambda decays:
  fDownComPIDnegPart[AliPID::kElectron] = 0.05;
  fDownComPIDposPart[AliPID::kElectron] = 0.05;

  fDownComPIDnegPart[AliPID::kMuon] = 0.05;
  fDownComPIDposPart[AliPID::kMuon] = 0.05;

  fDownComPIDnegPart[AliPID::kPion] = 0.05;
  fDownComPIDposPart[AliPID::kPion] = 0.05;

  fDownComPIDnegPart[AliPID::kKaon] = 0.05;
  fDownComPIDposPart[AliPID::kKaon] = 0.05;

  fDownComPIDnegPart[AliPID::kProton] = 0.05;
  fDownComPIDposPart[AliPID::kProton] = 0.05;*/
}

//_________________________________________________
AliTRDv0Info::AliTRDv0Info(const AliTRDv0Info &ref)
  : TObject()
  ,fQuality(ref.fQuality)
  ,fDCA(ref.fDCA)
  ,fPointingAngle(ref.fPointingAngle)
  ,fOpenAngle(ref.fOpenAngle)
  ,fPsiPair(ref.fPsiPair)
  ,fMagField(ref.fMagField)
  ,fRadius(ref.fRadius)
  ,fV0Momentum(ref.fV0Momentum)
  ,fTrackP(ref.fTrackP)
  ,fTrackN(ref.fTrackN)
  ,fNindex(ref.fNindex)
  ,fPindex(ref.fPindex)
  ,fInputEvent(ref.fInputEvent)
  ,fPrimaryVertex(ref.fPrimaryVertex)
{
  //
  // Copy constructor
  //

  memcpy(fPplus, ref.fPplus, 2*kNlayer*sizeof(Float_t));
  memcpy(fPminus, ref.fPminus, 2*kNlayer*sizeof(Float_t));
  memcpy(fDetPID, ref.fDetPID, 2*kNDaughters*kNDetectors*AliPID::kSPECIES*sizeof(Float_t));
  memcpy(fComPID, ref.fComPID, 2*kNDaughters*AliPID::kSPECIES*sizeof(Float_t));
  memcpy(fInvMass, ref.fInvMass, kNMomBins*kNDecays*sizeof(Double_t));
  memcpy(fArmenteros, ref.fArmenteros, kNDecays*sizeof(Bool_t));
  memcpy(fChi2ndf, ref.fChi2ndf, kNDecays*sizeof(Double_t));
  memcpy(fTPCdEdx, ref.fTPCdEdx, kNDaughters*sizeof(Float_t));

  //Upper limit for distance of closest approach of two daughter tracks :
  memcpy(fUpDCA, ref.fUpDCA, kNDecays*sizeof(Float_t));
  memcpy(fUpPointingAngle, ref.fUpPointingAngle, kNDecays*sizeof(Float_t));
  memcpy(fUpOpenAngle, ref.fUpOpenAngle, kNDecays*sizeof(Float_t));
  memcpy(fDownOpenAngle, ref.fDownOpenAngle, kNDecays*sizeof(Float_t));
  memcpy(fUpPsiPair, ref.fUpPsiPair, kNDecays*sizeof(Float_t));
  memcpy(fDownPsiPair, ref.fDownPsiPair, kNDecays*sizeof(Float_t));
  memcpy(fUpInvMass, ref.fUpInvMass, kNDecays*kNMomBins*sizeof(Double_t));
  memcpy(fDownInvMass, ref.fDownInvMass, kNDecays*sizeof(Double_t));
  memcpy(fUpChi2ndf, ref.fUpChi2ndf, kNDecays*sizeof(Double_t));
  memcpy(fUpRadius, ref.fUpRadius, kNDecays*sizeof(Float_t));
  memcpy(fDownRadius, ref.fDownRadius, kNDecays*sizeof(Float_t));
  memcpy(fDownTPCPIDneg, ref.fDownTPCPIDneg, AliPID::kSPECIES*sizeof(Float_t));
  memcpy(fDownTPCPIDpos, ref.fDownTPCPIDpos, AliPID::kSPECIES*sizeof(Float_t));
  memcpy(fDownComPIDneg, ref.fDownComPIDneg, AliPID::kSPECIES*sizeof(Float_t));
  memcpy(fDownComPIDpos, ref.fDownComPIDpos, AliPID::kSPECIES*sizeof(Float_t));
  memcpy(fDownComPIDnegPart, ref.fDownComPIDnegPart, AliPID::kSPECIES*sizeof(Float_t));
  memcpy(fDownComPIDposPart, ref.fDownComPIDposPart, AliPID::kSPECIES*sizeof(Float_t));
}

//_________________________________________________
void AliTRDv0Info::SetV0Info(AliESDv0 *esdv0)
{//Gets values of ESDv0 and daughter track properties
  //See header file for description of variables

  fQuality = Quality(esdv0);//Attributes an Int_t to the V0 due to quality cuts (= 1 if V0 is accepted, other integers depending on cut which excludes the vertex)    

  fRadius = Radius(esdv0);//distance from secondary vertex to primary vertex in x-y plane
      
  fDCA = esdv0->GetDcaV0Daughters();//distance of closest approach of two daughter tracks
      
  fPointingAngle = TMath::ACos(esdv0->GetV0CosineOfPointingAngle());// pointing angle (= angle between between vector from primary to secondary vertex and reconstructed momentum of V0 mother particle)
      
  fOpenAngle = OpenAngle(esdv0);//Opening angle between two daughter tracks
      
  fPsiPair = PsiPair(esdv0);//Angle between daughter momentum plane and plane perpendicular to magnetic field

  fV0Momentum = V0Momentum(esdv0);//Reconstructed momentum of the mother particle
      
  //4 decay types : conversions, K0s, Lambda, Anti-Lambda 
    //five particle types: electrons, muons, pions, kaons, protons (muons and kaons not involved)
  for(Int_t idecay(0), part1(-1), part2(-1); idecay < kNDecays; idecay++){

    fArmenteros[idecay]=Armenteros(esdv0, idecay);//Attribute the Armenteros yes/no decision for every decay type
    if(idecay == kLambda){ //protons and pions from Lambda
      part1 = AliPID::kProton;
      part2 = AliPID::kPion;
    } else if(idecay == kAntiLambda) { //antiprotons and pions from Anti-Lambda
      part1 = AliPID::kPion;
      part2 = AliPID::kProton;
    } else if(idecay == kK0s) {//pions from K0s
      part1 = part2 = AliPID::kPion;
    } else if(idecay == kGamma) {//electrons from conversions
      part1 = part2 = AliPID::kElectron;
    } 
    fInvMass[idecay] = InvMass(part1, part2, esdv0);//Calculate invariant mass for all of our four supposed decays

    // Comment out until bug fix is provided
    // A.Bercuci 14. July 2010
    //fChi2ndf[idecay] = KFChi2ndf(part1, part2,idecay);
   
  }
  //Gets all likelihood values from TPC, TOF and ITS PID for the fDetPID[kNDaughters][kNDetectors][AliPID::kSPECIES] array
  GetDetectorPID();
  //Bayesian combination of likelihoods from TPC and TOF
  CombinePID();
  //TPC dE/dx values for both tracks
  GetTPCdEdx();

}
//_________________________________________________
Float_t  AliTRDv0Info::V0Momentum(AliESDv0 *esdv0) const
{
  //
  // Reconstructed momentum of V0 mother particle
  //

  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};


  esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter; 
  esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;
  
  
  return TMath::Sqrt((mn[0]+mp[0])*(mn[0]+mp[0]) + (mn[1]+mp[1])*(mn[1]+mp[1])+(mn[2]+mp[2])*(mn[2]+mp[2]));
}

//_________________________________________________
Double_t AliTRDv0Info::InvMass(Int_t part1, Int_t part2, AliESDv0 *esdv0) const
{
  //
  // Invariant mass of reconstructed V0 mother
  //

  const Double_t kpmass[5] = {AliPID::ParticleMass(AliPID::kElectron),AliPID::ParticleMass(AliPID::kMuon),AliPID::ParticleMass(AliPID::kPion),AliPID::ParticleMass(AliPID::kKaon),AliPID::ParticleMass(AliPID::kProton)};
  //Masses of electrons, muons, pions, kaons and protons, as implemented in ROOT


  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};  

  esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;
  
  Double_t mass1 = kpmass[part1];//sets supposed rest masses for both daughters: positive
  Double_t mass2 = kpmass[part2];//negative   

  //Calculate daughters' energies :
  Double_t e1    = TMath::Sqrt(mass1*mass1+
            mp[0]*mp[0]+
            mp[1]*mp[1]+
            mp[2]*mp[2]);
  Double_t e2    = TMath::Sqrt(mass2*mass2+
            mn[0]*mn[0]+
            mn[1]*mn[1]+
            mn[2]*mn[2]);  

  //Sum of daughter momenta :   
  Double_t momsum =  
    (mn[0]+mp[0])*(mn[0]+mp[0])+
    (mn[1]+mp[1])*(mn[1]+mp[1])+
    (mn[2]+mp[2])*(mn[2]+mp[2]);

  //invariant mass :	  	     
  Double_t mInv = TMath::Sqrt((e1+e2)*(e1+e2)-momsum);

  return mInv;
  
}
//_________________________________________________
Float_t AliTRDv0Info::OpenAngle(AliESDv0 *esdv0)
{//Opening angle between two daughter tracks
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
    

  esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;

  
  fOpenAngle = TMath::ACos((mp[0]*mn[0] + mp[1]*mn[1] + mp[2]*mn[2])/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1] + mp[2]*mp[2])*TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1] + mn[2]*mn[2])));
  
  return fOpenAngle;
}

//_________________________________________________
Float_t AliTRDv0Info::PsiPair(AliESDv0 *esdv0)
{//Angle between daughter momentum plane and plane perpendicular to magnetic field
  Double_t x, y, z;
  esdv0->GetXYZ(x,y,z);//Reconstructed coordinates of V0; to be replaced by Markus Rammler's method in case of conversions!
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
  

  esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter; 


  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) -  TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis

  Double_t radiussum = TMath::Sqrt(x*x + y*y) + 50;//radius to which tracks shall be propagated

  Double_t momPosProp[3];
  Double_t momNegProp[3];
    
  AliExternalTrackParam nt(*fTrackN), pt(*fTrackP);
    
  fPsiPair = 4.;

  if(nt.PropagateTo(radiussum,fMagField) == 0)//propagate tracks to the outside
    fPsiPair =  -5.;
  if(pt.PropagateTo(radiussum,fMagField) == 0)
    fPsiPair = -5.;
  pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
  nt.GetPxPyPz(momNegProp);
  
  Double_t pEle =
    TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
  Double_t pPos =
    TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter
    
  Double_t scalarproduct =
    momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta
    
  Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

  fPsiPair =  TMath::Abs(TMath::ASin(deltat/chipair));  

  return fPsiPair; 

}
//_________________________________________________
Double_t AliTRDv0Info::KFChi2ndf(Int_t part1, Int_t part2,Int_t decay){
  //Calculates Kalman filter Chi2/NDF
  Int_t mothers[4]={22,310,3122,3122};

  const Double_t partMass=TDatabasePDG::Instance()->GetParticle(mothers[decay])->Mass();
  const Double_t massWidth[4] = {0.001, 0., 0., 0.};
 
  AliKFParticle *kfMother = CreateMotherParticle(fTrackP, fTrackN, part1, part2);
 
  // Lambda
  if(!kfMother) {
  return kFALSE;
  }
  
  // production vertex is set in the 'CreateMotherParticle' function
  kfMother->SetMassConstraint(partMass, massWidth[decay]);
 
  Double_t chi2ndf = (kfMother->GetChi2()/kfMother->GetNDF());
 
  if(kfMother)delete kfMother;
  return chi2ndf; 
}
//________________________________________________________________
AliKFParticle *AliTRDv0Info::CreateMotherParticle(AliESDtrack *pdaughter, AliESDtrack *ndaughter, Int_t pspec, Int_t nspec){
  //
  // Creates a mother particle
  //
  AliKFParticle pkfdaughter(*pdaughter, pspec);
  AliKFParticle nkfdaughter(*ndaughter, nspec);
  
 
  // Create the mother particle 
  AliKFParticle *m = new AliKFParticle(pkfdaughter, nkfdaughter);
 
  AliKFVertex improvedVertex = *fPrimaryVertex;
  improvedVertex += *m;
  m->SetProductionVertex(improvedVertex);
  

  return m;
}
//_________________________________________________
Int_t AliTRDv0Info::HasTrack(AliTRDtrackInfo * const track)
{
//Checks if track is a secondary vertex daughter (due to V0 finder)
  
  Int_t trackID(track->GetTrackId());//index of the track
  return HasTrack(trackID);
}

//_________________________________________________
Int_t AliTRDv0Info::HasTrack(Int_t trackID)
{
  //comparing index of track with indices of pos./neg. V0 daughter :
  if(fNindex==trackID) return -1;
  else if(fPindex==trackID) return 1;
  else return 0;
}

//_________________________________________________
void AliTRDv0Info::GetDetectorPID()
{//PID likelihoods from TPC, TOF, and ITS, for all particle species

  fTrackN->GetTPCpid(fDetPID[kNeg][kTPC]);
  fTrackP->GetTPCpid(fDetPID[kPos][kTPC]);
  fTrackN->GetTOFpid(fDetPID[kNeg][kTOF]);
  fTrackP->GetTOFpid(fDetPID[kPos][kTOF]);
  fTrackN->GetITSpid(fDetPID[kNeg][kITS]);
  fTrackP->GetITSpid(fDetPID[kPos][kITS]);

  Long_t statusN = fTrackN->GetStatus(); 
  Long_t statusP = fTrackP->GetStatus(); 
  
  if(!(statusN & AliESDtrack::kTPCpid)){
         for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	fDetPID[kNeg][kTPC][iPart] = 0.2;
      }    
  }
  if(!(statusN & AliESDtrack::kTOFpid)){
    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      fDetPID[kNeg][kTOF][iPart] = 0.2;
    }    
    
  }
  if(!(statusN & AliESDtrack::kITSpid)){
    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      fDetPID[kNeg][kITS][iPart] = 0.2;
    }    
  }
  if(!(statusP & AliESDtrack::kTPCpid)){
    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      fDetPID[kPos][kTPC][iPart] = 0.2;
    }    
  }
  if(!(statusP & AliESDtrack::kTOFpid)){
    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      fDetPID[kPos][kTOF][iPart] = 0.2;
    }    
    
  }
  if(!(statusP & AliESDtrack::kITSpid)){
    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
      fDetPID[kPos][kITS][iPart] = 0.2;
    }    
  }

}
//____________________________________________________________________________________
void AliTRDv0Info::CombinePID()
{
  Double_t partrat[AliPID::kSPECIES] = {0.208, 0.010, 0.662, 0.019, 0.101};
  
  for(Int_t iSign = 0; iSign < kNDaughters; iSign++)
    {
    for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++)
      {
      fComPID[iSign][iPart] = (partrat[iPart]*fDetPID[iSign][kTPC][iPart]*fDetPID[iSign][kTOF][iPart])/((partrat[0]*fDetPID[iSign][kTPC][0]*fDetPID[iSign][kTOF][0])+(partrat[1]*fDetPID[iSign][kTPC][1]*fDetPID[iSign][kTOF][1])+(partrat[2]*fDetPID[iSign][kTPC][2]*fDetPID[iSign][kTOF][2])+(partrat[3]*fDetPID[iSign][kTPC][3]*fDetPID[iSign][kTOF][3])+(partrat[4]*fDetPID[iSign][kTPC][4]*fDetPID[iSign][kTOF][4]));
      
      }
    }
}
//_________________________________________________
void AliTRDv0Info::GetTPCdEdx()
{
  fTPCdEdx[kNeg] = fTrackN->GetTPCsignal();
  fTPCdEdx[kPos] = fTrackP->GetTPCsignal();

}
//_________________________________________________
Bool_t AliTRDv0Info::TPCdEdxCuts(Int_t part, AliTRDtrackInfo * const track)
{
  //Bethe-Bloch lines
  Double_t alephParameters[5];
  
  // data
  alephParameters[0] = 0.0283086;
  alephParameters[1] = 2.63394e+01;
  alephParameters[2] = 5.04114e-11;
  alephParameters[3] = 2.12543e+00;
  alephParameters[4] = 4.88663e+00;
  

  Double_t deposit = 0;
  Float_t x = 0;
  if(HasTrack(track) == 1){
    x = fTrackP->P();
    deposit = fTPCdEdx[kPos];
  }
  else if(HasTrack(track) == -1){
    x = fTrackN->P();
    deposit = fTPCdEdx[kNeg];
  }
  else{
    printf("No track found");
    return 0;
  }
  if(x < 0.2)return 0;

  Float_t upLimits[5]={85,1000,50*AliExternalTrackParam::BetheBlochAleph(x/0.13957, alephParameters[0], alephParameters[1], alephParameters[2], alephParameters[3],  alephParameters[4])+6,1000,50*AliExternalTrackParam::BetheBlochAleph(x/0.93827, alephParameters[0], alephParameters[1], alephParameters[2], alephParameters[3],  alephParameters[4])+10};
  Float_t downLimits[5]={62,40,50*AliExternalTrackParam::BetheBlochAleph(x/0.13957, alephParameters[0], alephParameters[1], alephParameters[2], alephParameters[3],  alephParameters[4])-6,40,50*AliExternalTrackParam::BetheBlochAleph(x/0.93827, alephParameters[0], alephParameters[1], alephParameters[2], alephParameters[3],  alephParameters[4])-11};
  
  
  if(x < 0.7){
    downLimits[4]=90;
  }
  if(x < 1.25){
    upLimits[0] = 85;
  }
  else{
    downLimits[0] = 64;
  }


  if(deposit < downLimits[part])
    return 0;
  if(deposit > upLimits[part])
    return 0;

 
  return 1;

}
//_________________________________________________
Float_t AliTRDv0Info::Radius(AliESDv0 *esdv0)
{//distance from secondary vertex to primary vertex in x-y plane
  Double_t x, y, z;
  esdv0->GetXYZ(x,y,z); //Reconstructed coordinates of V0
  fRadius = TMath::Sqrt(x*x + y*y);
  return fRadius;

}

//_________________________________________________
Int_t AliTRDv0Info::Quality(AliESDv0 *const esdv0)
{ 
  //
  // Checking track and V0 quality status in order to exclude vertices based on poor information
  //

  Float_t nClsN;
  nClsN = fTrackN->GetTPCNcls();//number of found clusters in TPC for negative track
  Float_t nClsFN;
  nClsFN = fTrackN->GetTPCNclsF();//number of findable clusters in TPC for negative track
  Float_t nClsP;
  nClsP = fTrackP->GetTPCNcls();//number of found clusters in TPC for positive track
  Float_t nClsFP;
  nClsFP = fTrackP->GetTPCNclsF();//number of findable clusters in TPC for positive track
  
  fQuality = 0;


  if (!(esdv0->GetOnFlyStatus()))//accept only vertices from online V0 finder
    return -1;

  Float_t clsRatioN; 
  Float_t clsRatioP;

  if((nClsFN < 80) || (nClsFP < 80)) return -2;//reject all V0s where at least one track has less than 80 TPC clusters

 // Chi2 per TPC cluster
  Int_t nTPCclustersP = fTrackP->GetTPCclusters(0);
  Int_t nTPCclustersN = fTrackN->GetTPCclusters(0);
  Float_t chi2perTPCclusterP = fTrackP->GetTPCchi2()/Float_t(nTPCclustersP);
  Float_t chi2perTPCclusterN = fTrackN->GetTPCchi2()/Float_t(nTPCclustersN);
 
  if((chi2perTPCclusterN > 3.5)||(chi2perTPCclusterP > 3.5)) return -3;//reject all V0s where at least one track has a chi2 above 3.5
    
  clsRatioN = nClsN/nClsFN; //ratios of found to findable clusters in TPC 
  clsRatioP = nClsP/nClsFP;
  
  if((clsRatioN < 0.6)||(clsRatioP < 0.6))//exclude tracks with low ratio of found to findable TPC clusters
    return -4;
 
  if (!((fTrackP->GetStatus() &
  AliESDtrack::kTPCrefit)))//accept only vertices in which both tracks have TPC refit
    return -5;
  if (!((fTrackN->GetStatus() &
  AliESDtrack::kTPCrefit)))
    return -6;	
  if (fTrackP->GetKinkIndex(0)>0  ||
      fTrackN->GetKinkIndex(0)>0 )//exclude tracks with kinks
    return -7;
 
  if(!(V0SignCheck()))
       return -8;
  fQuality = 1;
  return fQuality;
}
//________________________________________________________________
Bool_t AliTRDv0Info::V0SignCheck(){
  //
  // Check if v0 daughters really carry opposite charges
  //
 
  Int_t qP = fTrackP->Charge();
  Int_t qN = fTrackN->Charge();

  if((qP*qN) != -1) return kFALSE;

  return kTRUE;
}
//___________________________________________________________________
Bool_t AliTRDv0Info::Armenteros(AliESDv0 *esdv0, Int_t decay){
  //
  // computes the Armenteros variables for given V0
  //
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};  
  Double_t mm[3] = {0,0,0};  
 
  if(V0SignCheck()){
    esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
    esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
  }
  else{
    esdv0->GetPPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
    esdv0->GetNPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
  }
  esdv0->GetPxPyPz(mm[0],mm[1],mm[2]); //reconstructed cartesian momentum components of mother

  TVector3 vecN(mn[0],mn[1],mn[2]);
  TVector3 vecP(mp[0],mp[1],mp[2]);
  TVector3 vecM(mm[0],mm[1],mm[2]);
  
  Double_t thetaP = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
  Double_t thetaN = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
  
  Double_t alfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/
    ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) ;
  Double_t qt = vecP.Mag()*sin(thetaP);

  Float_t ap[2];
  ap[0] = alfa;
  ap[1] = qt;

  Double_t LcutAP[2];//Lambda/Anti-Lambda cuts
  if(decay == 0){
    // armenteros cuts
    const Double_t cutAlpha[2] = {0.35, 0.45};   // [0.35, 0.45]
    const Double_t cutQT = 0.015;
    if(TMath::Abs(ap[0]) > cutAlpha[0] && TMath::Abs(ap[0]) < cutAlpha[1]) return kFALSE;
   
    if(ap[1] > cutQT) return kFALSE;
  }
  
  else if(decay == 1){
    const Double_t cutQT = 0.1075;
    const Double_t cutAP = 0.22 * TMath::Sqrt( TMath::Abs( (1-ap[0]*ap[0]/(0.92*0.92)) ) );
    if(ap[1] < cutQT) return kFALSE;
    if(ap[1] > cutAP) return kFALSE;
  }
  else if(decay == 2){
    const Double_t cutQT = 0.03;
    const Double_t cutAlpha = 0.7;  // VERY strong - should supress the overlap with K0
    LcutAP[0] = 1.0 - (ap[0]-0.7 * ap[0]-0.7)*1.1 - 0.87;
    if(TMath::Abs(ap[0]) > cutAlpha) return kFALSE;
    if(ap[1] < cutQT) return kFALSE;
    if(ap[1] > LcutAP[0]) return kFALSE;

  }
  else if(decay == 3){
    const Double_t cutQT = 0.03;
    const Double_t cutAlpha = 0.7;  // VERY strong - should supress the overlap with K0
    LcutAP[1] = 1.0 - (ap[0]+0.7 * ap[0]+0.7)*1.1 - 0.87;
    if(TMath::Abs(ap[0]) > cutAlpha) return kFALSE;
    if(ap[1] < cutQT) return kFALSE;
    if(ap[1] > LcutAP[1]) return kFALSE;
  }
  return kTRUE;
}
//_________________________________________________
Int_t AliTRDv0Info::GetPID(Int_t ipart, AliTRDtrackInfo *track)
{
  // Decides if track is accepted for one of the reference data samples
  Int_t cutCode = -99;
  if(!(track)) {
    AliError("No track info");
    return -1;
  }
  if(!HasTrack(track)){
    AliDebug(2, "Track not attached to v0.");
    return -2;
  }

  //translate ipart to decay (Anti-Lambda will be treated separately)
  Int_t iDecay = -1;
  switch(ipart){
  case AliPID::kElectron: iDecay = kGamma; break;
  case AliPID::kPion: iDecay = kK0s; break;
  case AliPID::kProton: iDecay = kLambda; break;
  default:
    AliWarning(Form("Hypothesis \"ipart=%d\" not handled", ipart));
    return -3;
  }

  //... it fulfills our quality criteria
  if(!(fQuality == 1)) return -4;
  //... distance of closest approach between daughters is reasonably small
  if((fDCA > fUpDCA[iDecay])) return -5;
  //... pointing angle between momentum of mother particle and vector from prim. to sec. vertex is small
  if((fPointingAngle > fUpPointingAngle[iDecay])) return -6;
  //... x-y plane distance of decay point to prim. vertex is bigger than a certain minimum value (for conversions)
  if((fRadius < fDownRadius[iDecay])) return -7;
  //...or smaller than a maximum value (for K0s)
  if((fRadius > fUpRadius[iDecay])) return -8;
  //... opening angle is close enough to zero (for conversions)
  if((fOpenAngle > fUpOpenAngle[iDecay])) return -9;
  //... Psi-pair angle is close enough to zero(for conversions)
  if((TMath::Abs(fPsiPair) > fUpPsiPair[iDecay])) return -10;
 


  //Mother momentum slots above/below 2.5 GeV
  Int_t iPSlot(fV0Momentum > 2.5);
  Int_t trackID(track->GetTrackId());

  //specific cut criteria :
  if(ipart == AliPID::kProton) {
    if((fInvMass[kK0s] < fUpInvMass[kK0s][iPSlot]) && (fInvMass[kK0s] > fDownInvMass[kK0s])) return -11;//explicit exclusion of K0s decays

    if(fOpenAngle < (0.3 - 0.2*fV0Momentum))return -9;

 

    //for proton sample: separate treatment of Lamba and Anti-Lambda decays:
    //for Anti-Lambda:
    //Combined PID likelihoods high enough for pi+ and anti-proton ; invariant mass calculated postulating these two particle species...
    //if((fComPID[kNeg][AliPID::kProton] > fDownComPIDneg[AliPID::kProton]) && (fComPID[kPos][AliPID::kPion] > fDownComPIDposPart[AliPID::kPion])) {
    //if((fDetPID[kNeg][kTPC][AliPID::kProton] > fDownTPCPIDneg[AliPID::kProton]) && (fDetPID[kPos][kTPC][AliPID::kPion] > fDownTPCPIDpos[AliPID::kPion])){
    if((TPCdEdxCuts(ipart, track))){//momentary solution: direct cut on TPC dE/dx
      if(fNindex == trackID) {//we're only interested in the anti-proton
	if(fArmenteros[kAntiLambda]){//Armenteros condition has to be fulfilled	  
	  if(fChi2ndf[kAntiLambda] < fUpChi2ndf[kAntiLambda]){//Kalman filter Chi2/NDF not allowed to be too large
	    if((fInvMass[kAntiLambda] < fUpInvMass[kAntiLambda][iPSlot]) && (fInvMass[kAntiLambda] > fDownInvMass[kAntiLambda])){  
        return 1;
	    } else cutCode = -15;
	  }
	  else cutCode =-14;
	}
	else cutCode = -13;
      }
    }
    else cutCode = -12;
    //for Lambda:
    //TPC PID likelihoods high enough for pi- and proton ; invariant mass calculated accordingly
    //if((fComPID[kNeg][AliPID::kPion] > fDownComPIDnegPart[AliPID::kPion]) && (fComPID[kPos][AliPID::kProton] > fDownComPIDpos[AliPID::kProton])) {
    //if((fDetPID[kNeg][kTPC][AliPID::kPion] > fDownTPCPIDneg[AliPID::kPion]) && (fDetPID[kPos][kTPC][AliPID::kProton] > fDownTPCPIDpos[AliPID::kProton])){
    if((TPCdEdxCuts(ipart, track))){//momentary solution: direct TPC dE/dx cuts
      if(fPindex == trackID) {
	if(fArmenteros[kLambda]){
	  if(fChi2ndf[kLambda] < fUpChi2ndf[kLambda]){
	    if((fInvMass[kLambda] < fUpInvMass[kLambda][iPSlot]) && (fInvMass[kLambda] > fDownInvMass[kLambda])){ 
        return 1;
	    } else cutCode = -15;
	  }
	  else cutCode = -14;
	}
	else cutCode = -13;
      }
    }
    else cutCode = -12;
    return cutCode;
  }
 
  //for K0s decays: equal TPC PID likelihood criteria for both daughters ; invariant mass calculated postulating two pions
  if(ipart == AliPID::kPion) {
    
    if(fOpenAngle < (1.0/(fV0Momentum + 0.3) - 0.1))
      return -9;
    
    //explicit exclusion of Lambda decays
    if((fInvMass[kLambda] < fUpInvMass[kLambda][iPSlot]) && (fInvMass[kLambda] > fDownInvMass[kLambda])) return -11;
    //explicit exclusion of Anti-Lambda decays
    if((fInvMass[kAntiLambda] < fUpInvMass[kAntiLambda][iPSlot]) && (fInvMass[kAntiLambda] > fDownInvMass[kAntiLambda])) return -11;
    
    //if((fDetPID[kNeg][kTPC][ipart] < fDownTPCPIDneg[ipart]) || (fDetPID[kPos][kTPC][ipart] < fDownTPCPIDpos[ipart])) return -12;
    if(!(TPCdEdxCuts(ipart, track))){//momentary solution: direct TPC dE/dx cuts
      return -12;
    }
  }
  
  
  //for photon conversions: equal combined PID likelihood criteria for both daughters ; invariant mass calculated postulating two electrons
  //No Lambda/K0s exclusion is provided, since these contributions hardly ever interfere with gamma invariant mass!
  //Float_t momentum(track->GetESDinfo()->GetOuterParam()->P());
  if(ipart == AliPID::kElectron) {
    //if(momentum > 1.75) {//since combined PID performs a little worse in simulations than TPC standalone for higher momenta, ONLY TPC PID is used here
    //if((fDetPID[kNeg][kTPC][ipart] < fDownTPCPIDneg[ipart]) || (fDetPID[kPos][kTPC][ipart] < fDownTPCPIDpos[ipart])) return -12;
    //} else {//for low momenta, combined PID from TOF and TPC is used to get rid of proton contamination
    //if((fComPID[kNeg][ipart] > fDownComPIDneg[ipart]) && (fComPID[kPos][ipart] > fDownComPIDpos[ipart])) return 1;
    //}
    if(!(TPCdEdxCuts(ipart, track))){//momentary solution for direct TPC dE/dx cut
      return -12;
    }
    
  }
  
 
  //Armenteros-Polanski cut
  if(!(fArmenteros[iDecay])) return -13;
  
  //Kalman filter Chi2/NDF cut
  if(fChi2ndf[iDecay] > fUpChi2ndf[iDecay]) return -14;

  //Invariant mass cut for K0s and photons, assuming two pions/two electrons as daughters:
 
  if((fInvMass[iDecay] > fUpInvMass[iDecay][iPSlot]) || (fInvMass[iDecay] < fDownInvMass[iDecay])) {
    return -15;
    
  }
   
  return 1;
}


//_________________________________________________
void AliTRDv0Info::Print(Option_t *opt) const
{
  printf("V0 P[%d] N[%d]\n", fPindex, fNindex);
  printf("  DCA[%5.3f] Radius[%5.3f]\n", fDCA, fRadius);
  printf("  Angles : Pointing[%5.3f] Open[%5.3f] Psi[%5.3f]\n", fPointingAngle, fOpenAngle, fPsiPair);
  if(strcmp(opt, "a")!=0) return;
  printf("  Reconstructed PID\n"
         "  sgn spec   ITS   TPC   TOF   COM\n");
  for(Int_t idt=0; idt<kNDaughters; idt++){
    printf("   %c", idt?'-':'+');
    for(Int_t is(0); is<AliPID::kSPECIES; is++){ 
      printf("%s%s%s", is==0?"   ":"       ", AliPID::ParticleShortName(is), (is==1||is==2)?"  ":"   ");
      for(Int_t id(0); id<kNDetectors; id++){
        printf("%5.1f ", 1.e2*fDetPID[idt][id][is]);
      }
      printf("%5.1f\n", 1.e2*fComPID[idt][is]);
    }
  }
}

//_________________________________________________
void AliTRDv0Info::SetV0tracks(AliESDtrack *p, AliESDtrack *n) 
{
  fTrackP = p; fPindex = p->GetID();
  fTrackN = n; fNindex = n->GetID();
}


