/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                Dielectron TrackRotator                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzRotation.h"
#include "TGenPhaseSpace.h"
#include "TSystem.h"

#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliESDtrack.h"

#include "AliDielectronVarManager.h"
#include "AliDielectronHelper.h"
#include "AliDielectronTrackRotator.h"


ClassImp(AliDielectronTrackRotator)

AliDielectronTrackRotator::AliDielectronTrackRotator() :
  TNamed(),
  fIterations(1),
  fRotationType(kRotateBothRandom),
  fStartAnglePhi(TMath::Pi()),
  fConeAnglePhi(TMath::Pi()/6.),
  fKeepLocalY(kFALSE),
  fRotateAroundMother(kFALSE),
  fTypeOfModernRotation(fRotateULS),
  fkArrTracksP(0x0),
  fkArrTracksN(0x0),
  fkArrTracksPRotation(0x0),
  fkArrTracksNRotation(0x0),
  fCurrentIteration(0),
  fCurrentTackP(0),
  fCurrentTackN(0),
  fCurrentPairPP(0),
  fCurrentPairMM(0),
  fLastPairSent(kFALSE),
  fEvent(0x0),
  fTrack1(),
  fTrack2(),
  fVTrackP(0x0),
  fVTrackN(0x0),
  fChargeTrack1(+1),
  fChargeTrack2(-1),
  fPdgLeg1(-11),
  fPdgLeg2(11),
  fSameTracks(kTRUE),
  fRotatedTracksP(),
  fRotatedTracksN(),
  fRotatedTracksWeightP(),
  fRotatedTracksWeightN(),
  fArrTrackPairs(),
  fMinimalPtCut(0.),
  fMaximalPtCut(999.),
  fMinimalEtaCut(-999),
  fMaximalEtaCut(999.),
  fWeight(0.),
  fRotateTrackCorrectionMap(),
  fUseAccMap(kTRUE)

{
  //
  // Default Constructor
  //
  gRandom->SetSeed();
}

//______________________________________________
AliDielectronTrackRotator::AliDielectronTrackRotator(const char* name, const char* title) :
  TNamed(name, title),
  fIterations(1),
  fRotationType(kRotateBothRandom),
  fStartAnglePhi(TMath::Pi()),
  fConeAnglePhi(TMath::Pi()/6.),
  fKeepLocalY(kFALSE),
  fRotateAroundMother(kFALSE),
  fTypeOfModernRotation(fRotateULS),
  fkArrTracksP(0x0),
  fkArrTracksN(0x0),
  fkArrTracksPRotation(0x0),
  fkArrTracksNRotation(0x0),
  fCurrentIteration(0),
  fCurrentTackP(0),
  fCurrentTackN(0),
  fCurrentPairPP(0),
  fCurrentPairMM(0),
  fLastPairSent(kFALSE),
  fEvent(0x0),
  fTrack1(),
  fTrack2(),
  fVTrackP(0x0),
  fVTrackN(0x0),
  fChargeTrack1(+1),
  fChargeTrack2(-1),
  fPdgLeg1(-11),
  fPdgLeg2(11),
  fSameTracks(kTRUE),
  fRotatedTracksP(),
  fRotatedTracksN(),
  fRotatedTracksWeightP(),
  fRotatedTracksWeightN(),
  fArrTrackPairs(),
  fMinimalPtCut(0.),
  fMaximalPtCut(999.),
  fMinimalEtaCut(-999),
  fMaximalEtaCut(999.),
  fWeight(0.),
  fRotateTrackCorrectionMap(),
  fUseAccMap(kTRUE)
{
  //
  // Named Constructor
  //
  gRandom->SetSeed();
}

//______________________________________________
AliDielectronTrackRotator::~AliDielectronTrackRotator()
{
  //
  // Default Destructor
  //
  
}

//______________________________________________
void AliDielectronTrackRotator::Reset()
{
  //
  // Reset the current iterators
  //
  fCurrentIteration=0;
  fCurrentTackP=0;
  fCurrentTackN=0;
  fCurrentPairPP=0;
  fCurrentPairMM=0;
}

//______________________________________________
Bool_t AliDielectronTrackRotator::NextCombination()
{
  //
  // Perform track rotation of the tracks in the track arrays as long as there are possible combinations
  //

  fSameTracks = kFALSE;
  // if no electron/positron candidates in the event, skip it
  Int_t nPos = fkArrTracksP->GetEntriesFast();
  Int_t nNeg = fkArrTracksN->GetEntriesFast();
  if (!fkArrTracksP || !fkArrTracksN ||
      (nPos < 1 || nNeg < 1)) {
    Reset();
    return kFALSE;
  }
  // Rotate both tracks around the virtual mother particle of the pair. This way the tracks should keep
  // all the kinematic correlations to the other tracks in the event. Then after rotation generate
  // pairs with one of the rotated tracks to all the other tracks in the sample. Rotations are 
  // calculated in CalculatePairsFromRotationAroundMother() and stored in fArrTrackPairs.
  // Residual logic is to send pair after pair until all pairs are sent
  if (fRotateAroundMother == kTRUE){
    if (fLastPairSent == kTRUE){
      fLastPairSent = kFALSE;
      //AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationMultiplicity, 1.);
      //AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracks, 0);
      return kFALSE;
    }

    // calculate rotated pairs, set weight according to the number of LS pairs. This should auto-normalize the rotated pairs
    if (fCurrentIteration == 0){

      const double numberLSPairs = (nPos*(nPos-1) + nNeg*(nNeg-1)) / 2.;
 
      CalculatePairsFromRotationAroundMother();
      //if (numberLSPairs > fArrTrackPairs.size()) CalculateLikeSignPairs(); 
      if (fArrTrackPairsPM.size() == 0) return kFALSE;
      
    }

    if(fCurrentIteration < fArrTrackPairsPM.size()){

      fTrack1  = fArrTrackPairsPM[fCurrentIteration].kf1;
      fTrack2  = fArrTrackPairsPM[fCurrentIteration].kf2;
      fVTrackP = &(fArrTrackPairsPM[fCurrentIteration].vtrack1);
      fVTrackN = &(fArrTrackPairsPM[fCurrentIteration].vtrack2);
      fChargeTrack1 = +1;
      fChargeTrack2 = -1;

      ++fCurrentIteration;
    } 

    else if(fCurrentPairPP < fArrTrackPairsPP.size()){

      fTrack1  = fArrTrackPairsPP[fCurrentPairPP].kf1;
      fTrack2  = fArrTrackPairsPP[fCurrentPairPP].kf2;
      fVTrackP = &(fArrTrackPairsPP[fCurrentPairPP].vtrack1);
      fVTrackN = &(fArrTrackPairsPP[fCurrentPairPP].vtrack2);
      fChargeTrack1 = +1;
      fChargeTrack2 = +1;

      ++fCurrentPairPP;
    } 

    else if(fCurrentPairMM < fArrTrackPairsMM.size()){
                            
      fTrack1  = fArrTrackPairsMM[fCurrentPairMM].kf1;
      fTrack2  = fArrTrackPairsMM[fCurrentPairMM].kf2;
      fVTrackP = &(fArrTrackPairsMM[fCurrentPairMM].vtrack1);
      fVTrackN = &(fArrTrackPairsMM[fCurrentPairMM].vtrack2);
      fChargeTrack1 = -1;
      fChargeTrack2 = -1;
                                                              
      ++fCurrentPairMM;
    } 


    // after sending last pair all variables should be re-initialized
    if (fCurrentIteration == fArrTrackPairsPM.size() && fCurrentPairPP == fArrTrackPairsPP.size() && fCurrentPairMM == fArrTrackPairsMM.size() ){
      Reset();
      fLastPairSent = kTRUE;
    }
    return kTRUE;
  }

  // Do standard rotation method
  else {
    if (!fkArrTracksP || !fkArrTracksP) {
      Reset();
      return kFALSE;
    }

    if (nPos==0||nNeg==0){
      Reset();
      return kFALSE;
    }

    if (fCurrentIteration==fIterations){
      fCurrentIteration=0;
      ++fCurrentTackP;
    }

    if (fCurrentTackP==nPos){
      ++fCurrentTackN;
      fCurrentTackP=0;
    }

    if (fCurrentTackN==nNeg){
      Reset();
      return kFALSE;
    }

    if (!RotateTracks()){
      Reset();
      return kFALSE;
    }

    ++fCurrentIteration;
    return kTRUE;
  }
  
}

//______________________________________________
Bool_t AliDielectronTrackRotator::RotateTracks()
{
  //
  // Actual track rotation
  // Find out particle type and perform the rotation
  //

  AliVTrack *trackP=dynamic_cast<AliVTrack*>(fkArrTracksP->UncheckedAt(fCurrentTackP));
  AliVTrack *trackN=dynamic_cast<AliVTrack*>(fkArrTracksN->UncheckedAt(fCurrentTackN));
  if(trackN == trackP){
    fSameTracks=kTRUE;
    return kTRUE;;
  }
  fSameTracks = kFALSE;

  fTrack1.Initialize();
  fTrack2.Initialize();
  fVTrackP=0x0;
  fVTrackN=0x0;
  if (!trackP||!trackN) return kFALSE;
  fTrack1+=AliKFParticle(*trackP,fPdgLeg1);
  fTrack2+=AliKFParticle(*trackN,fPdgLeg2);

  fVTrackP=trackP;
  fVTrackN=trackN;

  Double_t angle  = fStartAnglePhi+(2*gRandom->Rndm()-1)*fConeAnglePhi;
  Int_t    charge = TMath::Nint(gRandom->Rndm());
  if( fKeepLocalY){// only rotate by multiples of one TPC chamber size
    angle = (Int_t) ( angle /  TMath::Pi() * 9 ) ;
  }

  if (fRotationType==kRotatePositive||(fRotationType==kRotateBothRandom&&charge==0)){
    AliDielectronHelper::RotateKFParticle(&fTrack1, angle, fEvent);
  }

  if (fRotationType==kRotateNegative||(fRotationType==kRotateBothRandom&&charge==1)){
    AliDielectronHelper::RotateKFParticle(&fTrack2, angle, fEvent);
  }

  return kTRUE;
}

//______________________________________________
void AliDielectronTrackRotator::CalculatePairsFromRotationAroundMother()
{
  //
  // Track regeneration
  //

  Int_t nNeg = fkArrTracksN->GetEntriesFast();
  Int_t nPos = fkArrTracksP->GetEntriesFast();

  Int_t nNegPt = nNeg; 
  Int_t nPosPt = nPos;

  const double numberLSPairs = (nPos*(nPos-1) + nNeg*(nNeg-1));
  const double numberLS_PP = nPos*(nPos-1);
  const double numberLS_MM = nNeg*(nNeg-1);

  //const double numberRotatedPairs_PP = (fIterations * (nPos * nNeg)) * (nPos-1);
  //const double numberRotatedPairs_MM = (fIterations * (nPos * nNeg)) * (nNeg-1);

  const double numberRotatedPairs_PP = ((nPos * nNeg)) * (nPos-1);
  const double numberRotatedPairs_MM = ((nPos * nNeg)) * (nNeg-1);

  static const double electron_mass = AliPID::ParticleMass(AliPID::kElectron);
  //static const double electron_mass = 0.;
  //TGenPhaseSpace decayer;
  //const Double_t arr[2] = {electron_mass, electron_mass};

  //if (fTypeOfModernRotation == fRotateWithTGenPhaseSpace || fTypeOfModernRotation == fRotateULS){
  int ULSpair_counter = 0;
  int ULSpair_counter_PP = 0;
  int ULSpair_counter_MM = 0;
  double rotated_pairs = 0;  

    for (Int_t iPos = 0; iPos < nPos; ++iPos){
      // Loop over all positive particles
      AliVTrack *trackP = dynamic_cast<AliVTrack*>(fkArrTracksP->UncheckedAt(iPos));

      AliKFParticle KFpos_orig;
      KFpos_orig.Initialize();
      KFpos_orig += AliKFParticle(*trackP, fPdgLeg1);

      KFpos_orig.SetField(0.2);

      for (Int_t iNeg = 0; iNeg < nNeg; ++iNeg){
        // Loop over all negative particles and rotate them by rotation_angle
        AliVTrack *trackN = dynamic_cast<AliVTrack*>(fkArrTracksN->UncheckedAt(iNeg));

        AliKFParticle KFneg_orig;
        KFneg_orig.Initialize();
        KFneg_orig += AliKFParticle(*trackN, fPdgLeg2);

        KFneg_orig.SetField(0.2);

        // calculate axis to rotate from momentum components of virtual mother
        AliKFParticle KFmother(KFpos_orig,KFneg_orig); 

        if(KFmother.GetMass() < 0.01){
	  KFmother = AliKFParticle(KFpos_orig,KFneg_orig, kTRUE);
        }

        bool rotatedLS_PP = false;
        bool rotatedLS_MM = false; 
        unsigned int rotated_count = 0;
        unsigned int loop_count = 1;
        double tot_weight = 0;
        while (loop_count < fIterations+1){
          //const double rotation_angle = pow(-1.,loop_count) * 180. * TMath::Pi()/180.;
          const double rotation_angle = gRandom->Rndm() * TMath::TwoPi(); //TMath::Pi() / 2.;
          //const double rotation_angle = pow(-1.,loop_count) * gRandom->Rndm()  * TMath::Pi();
        
          double posWeightSum = 0; 
          double negWeightSum = 0;

  
          AliKFParticle KFpos = KFpos_orig;
          AliKFParticle KFneg = KFneg_orig;

          RotateKFParticle(&KFneg, rotation_angle, &KFmother, fEvent);
          RotateKFParticle(&KFpos, rotation_angle, &KFmother, fEvent);

	  Double_t xyz_P[3];
          Double_t pxpypz_P[3];
          Double_t cov_P[21];

          xyz_P[0] = KFpos.X();
          xyz_P[1] = KFpos.Y();
          xyz_P[2] = KFpos.Z();
         
          pxpypz_P[0] = KFpos.Px();
          pxpypz_P[1] = KFpos.Py();
          pxpypz_P[2] = KFpos.Pz();

	  
          Double_t xyz_N[3];
          Double_t pxpypz_N[3];
          Double_t cov_N[21];

          xyz_N[0] = KFneg.X();
          xyz_N[1] = KFneg.Y();
          xyz_N[2] = KFneg.Z();
         
          pxpypz_N[0] = KFneg.Px();
          pxpypz_N[1] = KFneg.Py();
          pxpypz_N[2] = KFneg.Pz();


          TParticle partP(-11,trackP->GetStatus(),trackP->GetMother(),trackN->GetMother(),0,0,pxpypz_P[0],pxpypz_P[1],pxpypz_P[2], trackP->E(), xyz_P[0],xyz_P[1],xyz_P[2],0);
          TParticle partN(11,trackN->GetStatus(),trackN->GetMother(),trackP->GetMother(),0,0,pxpypz_N[0],pxpypz_N[1],pxpypz_N[2], trackN->E(), xyz_N[0],xyz_N[1],xyz_N[2],0);

          AliESDtrack esdTrackN(&partN);
          AliESDtrack esdTrackP(&partP);

	  esdTrackN.SetLabel(trackN->GetLabel());
          esdTrackP.SetLabel(trackP->GetLabel());

  //        TLorentzVector LvecPos;
  //        LvecPos.SetPtEtaPhiM(KFpos_orig.GetPt(), KFpos_orig.GetEta(), KFpos_orig.GetPhi(), electron_mass);
  //        TLorentzVector LvecNeg;
  //        LvecNeg.SetPtEtaPhiM(KFneg_orig.GetPt(), KFneg_orig.GetEta(), KFneg_orig.GetPhi(), electron_mass);

  //        TLorentzVector LvecMom;
  //        LvecMom.SetPtEtaPhiM(KFmother.GetPt(), KFmother.GetEta(), KFmother.GetPhi(), KFmother.GetMass());
  //        TLorentzVector LvecMom2 = LvecPos + LvecNeg;


       //   AliESDtrack trackP_rot2_esd = RotateVTrack(trackP ,rotation_angle, &LvecMom, fEvent);
       //   AliESDtrack trackN_rot2_esd = RotateVTrack(trackN ,rotation_angle, &LvecMom, fEvent);

       //   AliESDtrack trackP_rot3_esd = RotateVTrack(trackP ,rotation_angle, &LvecMom2, fEvent);
       //   AliESDtrack trackN_rot3_esd = RotateVTrack(trackN ,rotation_angle, &LvecMom2, fEvent);


          if ((KFpos.GetPt()  > fMinimalPtCut  && KFpos.GetPt()  < fMaximalPtCut  &&
               KFpos.GetEta() > fMinimalEtaCut && KFpos.GetEta() < fMaximalEtaCut) &&
              (KFneg.GetPt()  > fMinimalPtCut  && KFneg.GetPt()  < fMaximalPtCut  &&
               KFneg.GetEta() > fMinimalEtaCut && KFneg.GetEta() < fMaximalEtaCut) ){
       

            loop_count++;

            Double_t weight_rotAng = 1.;
            //Double_t weight_rotAng = GetWeightFromRotation2(rotation_angle);

	    //has to be moved in the future?
            AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationAngle, weight_rotAng);
            
            //fArrTrackPairs.emplace_back(KFpos, KFneg, esdTrackP, esdTrackN, 0, weight_rotAng);
            fArrTrackPairsPM.emplace_back(KFpos, KFneg, esdTrackP, esdTrackN, 0, weight_rotAng);
            fArrTrackPairsPM.back().rotAng = rotation_angle;


            //if ((KFpos.GetPt()  > fMinimalPtCut  && KFpos.GetPt()  < fMaximalPtCut  &&
            //     KFpos.GetEta() > fMinimalEtaCut && KFpos.GetEta() < fMaximalEtaCut) ){
     
 
              const double weight_pos = GetWeightFromRotation(&KFpos);
              const double weight_neg = GetWeightFromRotation(&KFneg);

              nNegPt = nNeg; 
              nPosPt = nPos;

              //tot_weight += weight_pos * weight_neg;

              //if(weight_pos != 0 && weight_neg != 0){
              

                //double phiv = PhivPair(fEvent->GetMagneticField(),trackP->Charge(),trackN->Charge(),LvecPos.Vect(),LvecNeg.Vect());
                //if(phiv < 2.){
                 
                //fill selected particle into the corresponding track arrays
                fRotatedTracksP.push_back(KFpos);
                Double_t RotWeightP = 1.;

                for (int iPos2 = 0; iPos2 < nPos; ++iPos2){
                  if(iPos2 == iPos) continue;
                  //AliVTrack *trackPToPair = dynamic_cast<AliVTrack*>(fkArrTracksP->UncheckedAt(iPos2));
                  AliESDtrack *esdTrackPToPair = dynamic_cast<AliESDtrack*>(fkArrTracksP->UncheckedAt(iPos2));
                  AliKFParticle KFpos2;
                  KFpos2.Initialize();
                  KFpos2 += AliKFParticle(*esdTrackPToPair, fPdgLeg1);
                  if ((KFpos2.GetPt()  < fMinimalPtCut  || KFpos2.GetPt()  > fMaximalPtCut ||
                       KFpos2.GetEta() < fMinimalEtaCut || KFpos2.GetEta() > fMaximalEtaCut) ){
                    nPosPt--;
                    continue;
                  }
                  else rotatedLS_PP = true;
          
                  double weight = 1.; 
                  //double weight = GetWeightFromOpeningAngle(&KFpos, &KFpos2); 
                  //if(weight > RotWeightP)
                  //  RotWeightP = weight;
                  //weight = weight_pos * weight_neg * (1. /( nPos-1 ));  
                  //weight = weight_pos;// * weight_rotAng * GetWeightFromOpeningAngle(&KFpos, &KFpos2); 
                  weight *= weight_pos; 
                  //weight = weight_pos * GetWeightFromRotation(&KFpos2) * weight_rotAng * GetWeightFromOpeningAngle(&KFpos, &KFpos2) * (numberLS_PP / numberRotatedPairs_PP);  
                  //const double weight = weight_pos * GetWeightFromRotation(&KFpos2)* nPos / (fIterations * (nPos * nNeg));
                  //const double weight = weight_pos * GetWeightFromRotation(&KFpos2) * numberLS_PP / (fIterations * (nPos * nNeg) * (nPos-1));
                  //fArrTrackPairs.emplace_back(KFpos, KFpos2, esdTrackP, *esdTrackPToPair, 1, weight);
                  fArrTrackPairsPP.emplace_back(KFpos, KFpos2, esdTrackP, *esdTrackPToPair, 1, weight);
                  posWeightSum += weight;
                  fArrTrackPairsPP.back().rotAng = rotation_angle;
                }
              //}
            //}

            //if ((KFneg.GetPt()  > fMinimalPtCut  && KFneg.GetPt()  < fMaximalPtCut  &&
            //     KFneg.GetEta() > fMinimalEtaCut && KFneg.GetEta() < fMaximalEtaCut) ){
     
                
                fRotatedTracksN.push_back(KFneg);                       
                Double_t RotWeightN = 1.;


                for (int iNeg2 = 0; iNeg2 < nNeg; ++iNeg2){
                  if(iNeg2 == iNeg) continue;
                  //AliVTrack *trackNToPair = dynamic_cast<AliVTrack*>(fkArrTracksN->UncheckedAt(iNeg2));
                  AliESDtrack *esdTrackNToPair = dynamic_cast<AliESDtrack*>(fkArrTracksN->UncheckedAt(iNeg2));
                  AliKFParticle KFneg2;
                  KFneg2.Initialize();
                  KFneg2 += AliKFParticle(*esdTrackNToPair, fPdgLeg2);


                  if ((KFneg2.GetPt()  < fMinimalPtCut  || KFneg2.GetPt()  > fMaximalPtCut ||
                       KFneg2.GetEta() < fMinimalEtaCut || KFneg2.GetEta() > fMaximalEtaCut) ){
                    nNegPt--;
                    continue; 
                  }
                  else rotatedLS_MM = true;

		  
                  double weight = 1.;
                  //double weight = GetWeightFromOpeningAngle(&KFneg, &KFneg2);
                  // if(weight > RotWeightN)
                  //  RotWeightN = weight;
                  //weight = weight_neg * weight_pos * (1. /( nNeg-1 ));
                  //weight = weight_neg;// * weight_rotAng * GetWeightFromOpeningAngle(&KFneg, &KFneg2);
                  weight *= weight_neg;// * GetWeightFromOpeningAngle(&KFneg, &KFneg2);
                  //weight = weight_neg * GetWeightFromRotation(&KFneg2) * weight_rotAng * GetWeightFromOpeningAngle(&KFneg, &KFneg2) * (numberLS_MM / numberRotatedPairs_MM);
                  //const double weight = weight_neg * GetWeightFromRotation(&KFneg2) * nNeg / (fIterations * (nPos * nNeg));
                  //const double weight = weight_neg * GetWeightFromRotation(&KFneg2) * numberLS_MM / (fIterations * (nPos * nNeg) * (nNeg-1));
                  //fArrTrackPairs.emplace_back(KFneg, KFneg2, esdTrackN, *esdTrackNToPair, 2, weight);
                  fArrTrackPairsMM.emplace_back(KFneg, KFneg2, esdTrackN, *esdTrackNToPair, 2, weight);
                  negWeightSum += weight;
                  fArrTrackPairsMM.back().rotAng = rotation_angle;
		  Double_t dS, dS1;
                  KFneg.GetDStoParticle(KFneg2, dS, dS1);
		  //std::cout << "DS to particle: " << dS << "  " << dS1 << std::endl;

                }

                fRotatedTracksWeightP.push_back(RotWeightP);
                fRotatedTracksWeightN.push_back(RotWeightN);
                //int i1 = 0;
                //if(fArrTrackPairs.size()-1-i1 < 0){ 
                //  while (fArrTrackPairs[fArrTrackPairs.size()-1-i1].charged_tracks != 0){
                //     
                //     if(fArrTrackPairs[fArrTrackPairs.size()-1-i1].charged_tracks == 1){
                //       //fArrTrackPairs[fArrTrackPairs.size()-1-i1].weight /= posWeightSum / (numberLS_PP/ (fIterations*(nNeg*nPos)*(nPos-1)));// * (tot_weight/rotated_count);
                //       fArrTrackPairs[fArrTrackPairs.size()-1-i1].weight *= ((nPos-1)/posWeightSum) * (numberLS_PP/ (fIterations*(nNeg*nPos)*(nPos-1)));// * (tot_weight/rotated_count);
                //       //fArrTrackPairs[fArrTrackPairs.size()-1-i1].weight *= (numberLS_PP/ (fIterations*(nNeg*nPos)*(nPos-1)));
                //       //fArrTrackPairs[fArrTrackPairs.size()-1-i1].weight /= posWeightSum ;// * (tot_weight/rotated_count);

                //     }
                //      if(fArrTrackPairs[fArrTrackPairs.size()-1-i1].charged_tracks == 2){
                //       //fArrTrackPairs[fArrTrackPairs.size()-1-i1].weight /= negWeightSum / (numberLS_MM/ (fIterations*(nNeg*nPos)*(nNeg-1)));// * (tot_weight/rotated_count);
                //       fArrTrackPairs[fArrTrackPairs.size()-1-i1].weight *= ((nNeg-1)/negWeightSum) * (numberLS_MM/ (fIterations*(nNeg*nPos)*(nNeg-1)));// * (tot_weight/rotated_count);

                //       //fArrTrackPairs[fArrTrackPairs.size()-1-i1].weight /= negWeightSum ;// * (tot_weight/rotated_count);

                //     }
                //     i1++; 
                //    
                //     if(fArrTrackPairs.size()-1-i1 < 0){ 
                //       break;
                //     } 
                //   }
                // }


              //}
            //}// non zero acceptance weight
            //else continue; 
          } // if within detector acceptance
          //else loop_count++;
        } // end while loop
        int counter = 0;
        int LSpair_counter = 0;
        if(rotated_count !=0){
          ULSpair_counter++;
          if(rotatedLS_PP)  ULSpair_counter_PP++;
          if(rotatedLS_MM)  ULSpair_counter_MM++;
        }

         //const double numberLSPairsPt = (nPosPt*(nPosPt-1) + nNegPt*(nNegPt-1));
         //const double numberLSpt_PP = nPosPt*(nPosPt-1);
         //const double numberLSpt_MM = nNegPt*(nNegPt-1);
        
         //const double numberRotatedPairsPt_PP = (nPosPt-1);
         //const double numberRotatedPairsPt_MM = (nNegPt-1);

         for(Int_t pair_i = 0; pair_i < fArrTrackPairsPP.size(); pair_i++)
	   fArrTrackPairsPP[pair_i].weight *= numberRotatedPairs_PP/fArrTrackPairsPP.size();

         for(Int_t pair_i = 0; pair_i < fArrTrackPairsMM.size(); pair_i++)
	   fArrTrackPairsMM[pair_i].weight *= numberRotatedPairs_MM/fArrTrackPairsMM.size() ;
//TODO
//        for (size_t i = 0; i < rotated_count;){
//          if(fArrTrackPairs.size()-1-counter < 0){ 
//            break;
//          }
//          if(fArrTrackPairs[fArrTrackPairs.size()-1-counter].charged_tracks == 0){
//            i++;
//            rotated_pairs += (numberRotatedPairs_PP + numberRotatedPairs_MM);
//          }
//          else {
//            //fArrTrackPairs[fArrTrackPairs.size()-1-counter].weight /= (( ((nNeg-1) + (nPos-1)) )/rotated_count);
//              fArrTrackPairs[fArrTrackPairs.size()-1-counter].weight /= rotated_count;// * (tot_weight/rotated_count);
//            //fArrTrackPairs[fArrTrackPairs.size()-1-counter].weight /= tot_weight;
//            LSpair_counter++;
//          }
//          counter++; 
//          
//        }

        //if(rotated_count!=0) rotated_pairs += numberRotatedPairs_PP/rotated_count + numberRotatedPairs_MM/rotated_count;
        
        //rotated_pairs += numberRotatedPairs_PP/rotated_count + numberRotatedPairs_MM/rotated_count;
        //} // end fRotateULS cause
      }
    }
  //}
 

   //const double numberLSPairsPt = (nPosPt*(nPosPt-1) + nNegPt*(nNegPt-1));
   //const double numberLSpt_PP = nPosPt*(nPosPt-1);
   //const double numberLSpt_MM = nNegPt*(nNegPt-1);
   
   //const double numberRotatedPairs_PP = (fIterations * (nPos * nNeg)) * (nPos-1);
   //const double numberRotatedPairs_MM = (fIterations * (nPos * nNeg)) * (nNeg-1);
   
   const double numberRotatedPairsPt_PP = (ULSpair_counter) * (nPosPt-1);
   const double numberRotatedPairsPt_MM = (ULSpair_counter) * (nNegPt-1);


 

  if (fArrTrackPairs.size() == 0) return;
 
  fWeight = 0.;
  double numberRMPairs = 0;
  double numberRMPairs_PP = 0;
  double numberRMPairs_MM = 0;
  for (size_t i = 0; i < fArrTrackPairs.size(); ++i){
    if (fArrTrackPairs[i].charged_tracks == 1){
      //fArrTrackPairs[i].weight /= (ULSpair_counter_PP * ( nPosPt-1 ) );
      //fArrTrackPairs[i].weight /= (ULSpair_counter * ( nPosPt-1 ) )/numberLSpt_PP;
        fArrTrackPairs[i].weight /= ((ULSpair_counter * ( nPos-1 )*fIterations ))/numberLS_PP;
        //fArrTrackPairs[i].weight /= ((ULSpair_counter * (nPosPt-1) *fIterations ));
 //       fArrTrackPairs[i].weight /= ((nNeg *fIterations ));

      numberRMPairs_PP++;
    }
    if (fArrTrackPairs[i].charged_tracks == 2){
      //fArrTrackPairs[i].weight /= (ULSpair_counter_MM * ( nNegPt-1 ) );
      //fArrTrackPairs[i].weight /= (ULSpair_counter * ( nNegPt-1 ) )/numberLSpt_MM;
      fArrTrackPairs[i].weight /= ((ULSpair_counter * ( nNeg-1 )*fIterations ))/numberLS_MM;
      //fArrTrackPairs[i].weight /= ((nPos * fIterations ));
      numberRMPairs_MM++; 
    }
    if (fArrTrackPairs[i].charged_tracks != 0 ){
      fWeight += fArrTrackPairs[i].weight;
      numberRMPairs++;
    }
  }
  for (size_t i = 0; i < fArrTrackPairs.size(); ++i){
    if (fArrTrackPairs[i].charged_tracks == 1)
      fArrTrackPairs[i].weight /= fWeight/numberLSPairs;
    if (fArrTrackPairs[i].charged_tracks == 2)
      fArrTrackPairs[i].weight /= fWeight/numberLSPairs;
  }
  fWeight = 0 ;
  for (size_t i = 0; i < fArrTrackPairs.size(); ++i){
    if (fArrTrackPairs[i].charged_tracks != 0){
      fWeight += fArrTrackPairs[i].weight;
    }
  }

  //TODO
  if (fWeight == 0) fWeight = 10000;
  //for (size_t i = 0; i < RotatedTack_counter; ++i){
  //  //fArrTrackPairs[i].weight = fArrTrackPairs[i].weight * numberLS_PP / fWeight;
  //  fArrTrackPairs[i].weight = fArrTrackPairs[i].weight * numberLS_PP / numberRMPairs;
  //  //std::cout << "weight PP: " << fArrTrackPairs[i].weight << std::endl;
  //}
  //for (size_t i = RotatedTack_counter; i < numberRMPairs; ++i){
  //  //fArrTrackPairs[i].weight = fArrTrackPairs[i].weight * numberLS_MM / fWeight;
  //  fArrTrackPairs[i].weight = fArrTrackPairs[i].weight * numberLS_MM / numberRMPairs;
  //  //std::cout << "weight MM: " << fArrTrackPairs[i].weight << std::endl; 
  //}


  double weight_LS_pairs = 0;
  if(numberRMPairs != 0) weight_LS_pairs = numberLSPairs / numberRMPairs;
  double weight_LS_pairs_PP = 0;
  if(numberRMPairs_PP != 0)  weight_LS_pairs_PP = numberLSPairs / numberRMPairs_PP;
  double weight_LS_pairs_MM = 0;
  if(numberRMPairs_MM != 0)  weight_LS_pairs_MM = numberLSPairs / numberRMPairs_MM;

  AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationMultiplicity   , weight_LS_pairs);
  AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationMultiplicity_PP, weight_LS_pairs_PP);
  AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationMultiplicity_MM, weight_LS_pairs_MM);

 
  //AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationMultiplicity   , numberLSPairs / numberRMPairs);
  //AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationMultiplicity_PP, numberLS_PP / numberRMPairs_PP);
  //AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationMultiplicity_MM, numberLS_MM / numberRMPairs_MM);
  AliDielectronVarManager::SetValue(AliDielectronVarManager::kNumberOfLSPairs, numberLSPairs);
  AliDielectronVarManager::SetValue(AliDielectronVarManager::kNumberOfRotatedPairs, numberRMPairs);

  AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracksForPairSum, fWeight); // should be identical to number of LS pairs

  return;
}


//______________________________________________
void AliDielectronTrackRotator::CalculateLikeSignPairs(){

  Int_t nNeg = fkArrTracksN->GetEntriesFast();
  Int_t nPos = fkArrTracksP->GetEntriesFast();

  for (Int_t iPos = 0; iPos < nPos; ++iPos){
    // Loop over all positive particles
    AliVTrack *trackP = dynamic_cast<AliVTrack*>(fkArrTracksP->UncheckedAt(iPos));
    AliKFParticle KFpos;
    KFpos.Initialize();
    KFpos += AliKFParticle(*trackP, fPdgLeg1);
    const double weight_pos = GetWeightFromRotation(&KFpos);
    for (Int_t iPos2 = iPos+1; iPos2 < nPos; ++iPos2){
      // Loop over all positive particles
      AliVTrack *trackP2 = dynamic_cast<AliVTrack*>(fkArrTracksP->UncheckedAt(iPos2));
      AliKFParticle KFpos2;
      KFpos2.Initialize();
      KFpos2 += AliKFParticle(*trackP2, fPdgLeg1);
      const double weight = weight_pos * GetWeightFromRotation(&KFpos2);
      //fArrTrackPairs.emplace_back(KFpos, KFpos2, trackP, trackP2, 1, weight);
    }
  }

  for (Int_t iNeg = 0; iNeg < nNeg; ++iNeg){
    // Loop over all positive particles
    AliVTrack *trackN = dynamic_cast<AliVTrack*>(fkArrTracksN->UncheckedAt(iNeg));
    AliKFParticle KFneg;
    KFneg.Initialize();
    KFneg += AliKFParticle(*trackN, fPdgLeg1);
    const double weight_neg = GetWeightFromRotation(&KFneg);
    for (Int_t iNeg2 = iNeg+1; iNeg2 < nNeg; ++iNeg2){
      // Loop over all positive particles
      AliVTrack *trackN2 = dynamic_cast<AliVTrack*>(fkArrTracksN->UncheckedAt(iNeg2));
      AliKFParticle KFneg2;
      KFneg2.Initialize();
      KFneg2 += AliKFParticle(*trackN2, fPdgLeg1);
      const double weight = weight_neg * GetWeightFromRotation(&KFneg2);
      //fArrTrackPairs.emplace_back(KFneg, KFneg2, trackN, trackN2, 2, weight);
    }
  }
}


//______________________________________________
void AliDielectronTrackRotator::SetRotatedTrackWeightMap(TString filename, TString histoname){
  TFile* file = TFile::Open(filename.Data(), "READ");
  printf("%p\n", file);
  if (file == 0x0){
    gSystem->Exec(Form("alien_cp alien://%s .",filename.Data()));
    printf("Copy rotated track map from Alien\n");
    TObjArray *arrNames=filename.Tokenize("/");
    TString name = arrNames->Last()->GetName();
    file = TFile::Open(Form("%s", name.Data()));
  }
  else {
    printf("Track Correction Map loaded\n");
  }
  if (file == nullptr){
    AliFatal(Form("Rotated-Track-Weighting file %s not found!", filename.Data()));
  }
  fRotateTrackCorrectionMap = *(dynamic_cast<TH3F*>(file->Get(histoname.Data())));
  if (&fRotateTrackCorrectionMap == nullptr){
    AliFatal(Form("Weighting histogram %s not found!", histoname.Data()));
  }
  fRotateTrackCorrectionMap.SetDirectory(0); 
  file->Close();
}

//______________________________________________
Double_t AliDielectronTrackRotator::GetWeightFromRotation(AliKFParticle* part){
  if(!fUseAccMap){
    return 1;
  }
  else{
    int bin_pt  = fRotateTrackCorrectionMap.GetXaxis()->FindBin(part->GetPt() );
    const int bin_eta = fRotateTrackCorrectionMap.GetYaxis()->FindBin(part->GetEta());
   
    if(bin_pt < fRotWeight_minPtBin){ 
      if(fRotWeight_minPtBin <= 1) bin_pt = 1;
      else bin_pt = fRotWeight_minPtBin;
    }
    if(bin_pt > fRotWeight_maxPtBin){ 
      if(fRotWeight_maxPtBin >= fRotateTrackCorrectionMap.GetXaxis()->GetNbins()) bin_pt = fRotateTrackCorrectionMap.GetXaxis()->GetNbins();
      else bin_pt = fRotWeight_maxPtBin;
    }

    double phi = part->GetPhi();
    if (phi < 0) phi += TMath::TwoPi();
    const int bin_phi = fRotateTrackCorrectionMap.GetZaxis()->FindBin(phi);
    
    Double_t weight = fRotateTrackCorrectionMap.GetBinContent(bin_pt, bin_eta, bin_phi);

    Int_t i = 1;
    while(weight <= 0){
      weight = fRotateTrackCorrectionMap.GetBinContent(bin_pt-i, bin_eta, bin_phi);  
      i++;
      if (bin_pt-i <= 1)
        break;
    }
    
    return weight;
  }
}


//______________________________________________
void AliDielectronTrackRotator::SetRotatedPairWeightMap(TString filename, TString histoname){
  TFile* file = TFile::Open(filename.Data(), "READ");
  printf("%p\n", file);
  if (file == 0x0){
    gSystem->Exec(Form("alien_cp alien://%s .",filename.Data()));
    printf("Copy rotated pair map from Alien\n");
    TObjArray *arrNames=filename.Tokenize("/");
    TString name = arrNames->Last()->GetName();
    file = TFile::Open(Form("%s", name.Data()));
  }
  else {
    printf("Pair Correction Map loaded\n");
  }
  if (file == nullptr){
    AliFatal(Form("Rotated-Pair-Weighting file %s not found!", filename.Data()));
  }
  fRotatePairCorrectionMap = *(dynamic_cast<TH2F*>(file->Get(histoname.Data())));
  if (&fRotatePairCorrectionMap == nullptr){
    AliFatal(Form("Weighting histogram %s not found!", histoname.Data()));
  }
  fRotatePairCorrectionMap.SetDirectory(0); 
  file->Close();
}


//______________________________________________
Double_t AliDielectronTrackRotator::GetWeightFromOpeningAngle(AliKFParticle* KFpos, AliKFParticle* KFneg){
  if(!fUseAccMap){
    return 1;
  }

  //needs to be changed to KF baseline
  else{
    // calculate axis to rotate from momentum components of virtual mother
    static const double electron_mass = AliPID::ParticleMass(AliPID::kElectron);
    TLorentzVector LvecPos;
    LvecPos.SetPtEtaPhiM(KFpos->GetPt(), KFpos->GetEta(), KFpos->GetPhi(), electron_mass);
    TLorentzVector LvecNeg;
    LvecNeg.SetPtEtaPhiM(KFneg->GetPt(), KFneg->GetEta(), KFneg->GetPhi(), electron_mass);
    TLorentzVector LvecMother = LvecPos + LvecNeg;
    Double_t ptee         = LvecMother.Pt();
    Double_t openingAngle = LvecPos.Angle(LvecNeg.Vect());

    if(openingAngle > 0.5) return 1.;

    int bin_ptee  = fRotatePairCorrectionMap.GetYaxis()->FindBin(ptee);
    const int bin_opAng = fRotatePairCorrectionMap.GetXaxis()->FindBin(openingAngle);
   

    double weight = fRotatePairCorrectionMap.GetBinContent(bin_opAng, bin_ptee);
    Int_t i = 1;
    while(weight <= 0){
      weight = fRotatePairCorrectionMap.GetBinContent(bin_opAng, bin_ptee-i);  
      i++;
      if (bin_ptee-i <= 1)
        return 1.;
    } 
    return weight;
  }
}

//______________________________________________
void AliDielectronTrackRotator::SetRotatedPairWeightMap2(TString filename, TString histoname){
  TFile* file = TFile::Open(filename.Data(), "READ");
  printf("%p\n", file);
  if (file == 0x0){
    gSystem->Exec(Form("alien_cp alien://%s .",filename.Data()));
    printf("Copy rotated pair map from Alien\n");
    TObjArray *arrNames=filename.Tokenize("/");
    TString name = arrNames->Last()->GetName();
    file = TFile::Open(Form("%s", name.Data()));
  }
  else {
    printf("Pair Correction Map loaded\n");
  }
  if (file == nullptr){
    AliFatal(Form("Rotated-Pair-Weighting2 file %s not found!", filename.Data()));
  }
  fRotatePairCorrectionMap2 = *(dynamic_cast<TH1F*>(file->Get(histoname.Data())));
  if (&fRotatePairCorrectionMap2 == nullptr){
    AliFatal(Form("Weighting histogram %s not found!", histoname.Data()));
  }
  fRotatePairCorrectionMap2.SetDirectory(0); 
  file->Close();
}


//______________________________________________
Double_t AliDielectronTrackRotator::GetWeightFromRotation2(Double_t rotAng){

    const int bin_rotAng = fRotatePairCorrectionMap2.GetXaxis()->FindBin(rotAng);
    Double_t weight = fRotatePairCorrectionMap2.GetBinContent(bin_rotAng);
    return weight;
}



//______________________________________________
Double_t AliDielectronTrackRotator::PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2) //const
{
  /// Following the idea to use opening of collinear pairs in magnetic field from e.g. PHENIX
  /// to identify conversions. Angle between ee plane and magnetic field is calculated (0 to pi).
  /// Due to tracking to the primary vertex, conversions with no intrinsic opening angle
  /// always end up as pair in "cowboy" configuration. The function as defined here then
  /// returns values close to pi.
  /// Correlated Like Sign pairs (from double conversion / dalitz + conversion) may show up
  /// at pi or at 0 depending on which leg has the higher momentum. (not checked yet)
  /// This expected ambiguity is not seen due to sorting of track arrays in this framework.
  /// To reach the same result as for ULS (~pi), the legs are flipped for LS.
  /// from PWGDQ/dielectron/AliDielectronPair.cxx

  //Define local buffer variables for leg properties
  Double_t px1=-9999.,py1=-9999.,pz1=-9999.;
  Double_t px2=-9999.,py2=-9999.,pz2=-9999.;

  TVector3 fD1=dau1;
  TVector3 fD2=dau2;
  Int_t    d1Q=charge1;
  //Int_t    d2Q=charge2;
  
  if (charge1*charge2 > 0.) { // Like Sign
    if(MagField<0){ // inverted behaviour
      if(d1Q>0){
        px1 = fD1.Px();   py1 = fD1.Py();   pz1 = fD1.Pz();
        px2 = fD2.Px();   py2 = fD2.Py();   pz2 = fD2.Pz();
      }else{
        px1 = fD2.Px();   py1 = fD2.Py();   pz1 = fD2.Pz();
        px2 = fD1.Px();   py2 = fD1.Py();   pz2 = fD1.Pz();
      }
    }else{
      if(d1Q>0){
        px1 = fD2.Px();   py1 = fD2.Py();   pz1 = fD2.Pz();
        px2 = fD1.Px();   py2 = fD1.Py();   pz2 = fD1.Pz();
      }else{
        px1 = fD1.Px();   py1 = fD1.Py();   pz1 = fD1.Pz();
        px2 = fD2.Px();   py2 = fD2.Py();   pz2 = fD2.Pz();
      }
    }
  }
  else { // Unlike Sign
    if(MagField>0){ // regular behaviour
      if(d1Q>0){
        px1 = fD1.Px();
        py1 = fD1.Py();
        pz1 = fD1.Pz();

        px2 = fD2.Px();
        py2 = fD2.Py();
        pz2 = fD2.Pz();
      }else{
        px1 = fD2.Px();
        py1 = fD2.Py();
        pz1 = fD2.Pz();

        px2 = fD1.Px();
        py2 = fD1.Py();
        pz2 = fD1.Pz();
      }
    }else{
      if(d1Q>0){
        px1 = fD2.Px();
        py1 = fD2.Py();
        pz1 = fD2.Pz();

        px2 = fD1.Px();
        py2 = fD1.Py();
        pz2 = fD1.Pz();
      }else{
        px1 = fD1.Px();
        py1 = fD1.Py();
        pz1 = fD1.Pz();

        px2 = fD2.Px();
        py2 = fD2.Py();
        pz2 = fD2.Pz();
      }
    }
  }

  Double_t px = px1+px2;
  Double_t py = py1+py2;
  Double_t pz = pz1+pz2;
  Double_t dppair = TMath::Sqrt(px*px+py*py+pz*pz);

  //unit vector of (pep+pem)
  Double_t pl = dppair;
  Double_t ux = px/pl;
  Double_t uy = py/pl;
  Double_t uz = pz/pl;
  Double_t ax = uy/TMath::Sqrt(ux*ux+uy*uy);
  Double_t ay = -ux/TMath::Sqrt(ux*ux+uy*uy);

  //momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
  //Double_t ptep = iep->Px()*ax + iep->Py()*ay;
  //Double_t ptem = iem->Px()*ax + iem->Py()*ay;

  Double_t pxep = px1;
  Double_t pyep = py1;
  Double_t pzep = pz1;
  Double_t pxem = px2;
  Double_t pyem = py2;
  Double_t pzem = pz2;

  //vector product of pep X pem
  Double_t vpx = pyep*pzem - pzep*pyem;
  Double_t vpy = pzep*pxem - pxep*pzem;
  Double_t vpz = pxep*pyem - pyep*pxem;
  Double_t vp = sqrt(vpx*vpx+vpy*vpy+vpz*vpz);
  //Double_t thev = acos(vpz/vp);

  //unit vector of pep X pem
  Double_t vx = vpx/vp;
  Double_t vy = vpy/vp;
  Double_t vz = vpz/vp;

  //The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
  Double_t wx = uy*vz - uz*vy;
  Double_t wy = uz*vx - ux*vz;
  //Double_t wz = ux*vy - uy*vx;
  //Double_t wl = sqrt(wx*wx+wy*wy+wz*wz);
  // by construction, (wx,wy,wz) must be a unit vector.
  // measure angle between (wx,wy,wz) and (ax,ay,0). The angle between them
  // should be small if the pair is conversion.
  // this function then returns values close to pi!
  Double_t cosPhiV = wx*ax + wy*ay;
  Double_t phiv = TMath::ACos(cosPhiV);

  return phiv;
}

void AliDielectronTrackRotator::RotateKFParticle(AliKFParticle * kfParticle,Double_t angle, AliKFParticle * kfMother, const AliVEvent * const ev){
  // Before rotate needs to be moved to the vertex position; move back after rotation

  Double_t dx = kfMother->GetX();
  Double_t dy = kfMother->GetY();
  Double_t dz = kfMother->GetZ();

  if (ev){
    dx = ev->GetPrimaryVertex()->GetX()-0.;
    dy = ev->GetPrimaryVertex()->GetY()-0.;
    dz = ev->GetPrimaryVertex()->GetZ()-0.;
  }


  kfParticle->X() = kfParticle->GetX() - dx;
  kfParticle->Y() = kfParticle->GetY() - dy;
  kfParticle->Z() = kfParticle->GetZ() - dz;


  kfMother->X() = kfMother->GetX() - dx;
  kfMother->Y() = kfMother->GetY() - dy;
  kfMother->Z() = kfMother->GetZ() - dz;
 

  TLorentzVector LvecMother;
  LvecMother.SetPtEtaPhiM(kfMother->GetPt(), kfMother->GetEta(), kfMother->GetPhi(), kfMother->GetMass());
  TVector3 vec3Axis = (LvecMother).Vect();

  TRotation trans;
  trans.Rotate(angle,vec3Axis);
  

  // Rotate the kf particle
  //Double_t c = cos(angle);
  //Double_t s = sin(angle);

  Double_t mA[8][ 8];
  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<8; j++){
      mA[i][j] = 0;
    }
  }
  //mA[0][0] =  c;  mA[0][1] = s;
  //mA[1][0] = -s;  mA[1][1] = c;
  //mA[3][3] =  c;  mA[3][4] = s;
  //mA[4][3] = -s;  mA[4][4] = c;

  mA[0][0] =  trans.XX();  mA[0][1] = trans.XY();  mA[0][2] = trans.XZ();
  mA[1][0] =  trans.YX();  mA[1][1] = trans.YY();  mA[1][2] = trans.YZ();
  mA[2][0] =  trans.ZX();  mA[2][1] = trans.ZY();  mA[2][2] = trans.ZZ();

  mA[3][3] =  trans.XX();  mA[3][4] = trans.XY();  mA[3][5] = trans.XZ();
  mA[4][3] =  trans.YX();  mA[4][4] = trans.YY();  mA[4][5] = trans.YZ();
  mA[5][3] =  trans.ZX();  mA[5][4] = trans.ZY();  mA[5][5] = trans.ZZ();

  mA[6][6] = 1; mA[7][7] = 1;

  Double_t mAC[8][8];
  Double_t mAp[8];

  for( Int_t i=0; i<8; i++ ){
    mAp[i] = 0;
    for( Int_t k=0; k<8; k++){
      mAp[i]+=mA[i][k] * kfParticle->GetParameter(k);
    }
  }

  for( Int_t i=0; i<8; i++){
    kfParticle->Parameter(i) = mAp[i];
  }

  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<8; j++ ){
      mAC[i][j] = 0;
      for( Int_t k=0; k<8; k++ ){
        mAC[i][j]+= mA[i][k] * kfParticle->GetCovariance(k,j);
      }
    }
  }

  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<=i; j++ ){
      Double_t xx = 0;
      for( Int_t k=0; k<8; k++){
        xx+= mAC[i][k]*mA[j][k];
      }
      kfParticle->Covariance(i,j) = xx;
    }
  }

  kfParticle->X() = kfParticle->GetX() + dx;
  kfParticle->Y() = kfParticle->GetY() + dy;
  kfParticle->Z() = kfParticle->GetZ() + dz;

  kfMother->X() = kfMother->GetX() + dx;
  kfMother->Y() = kfMother->GetY() + dy;
  kfMother->Z() = kfMother->GetZ() + dz;

}


Int_t AliDielectronTrackRotator::IJ( Int_t i, Int_t j ){ 
    return ( j<=i ) ? i*(i+1)/2+j :j*(j+1)/2+i;
}

AliESDtrack AliDielectronTrackRotator::RotateVTrack(AliVTrack * Vtrack,Double_t angle, TLorentzVector* LvecMother, const AliVEvent * const ev){
  // Before rotate needs to be moved to the vertex position; move back after rotation

  Double_t dx = 0;
  Double_t dy = 0;
  Double_t dz = 0;

  if (ev){
    dx = ev->GetPrimaryVertex()->GetX()-0.;
    dy = ev->GetPrimaryVertex()->GetY()-0.;
    dz = ev->GetPrimaryVertex()->GetZ()-0.;
  }


  Double_t cov[21] = {0};
  Vtrack->GetCovarianceXYZPxPyPz(cov);
  const Double_t* param = Vtrack->GetParameter();

  TVector3 vec3Axis = (LvecMother)->Vect();

  TRotation trans;
  trans.Rotate(angle,vec3Axis);
  

  // Rotate the kf particle
  //Double_t c = cos(angle);
  //Double_t s = sin(angle);

  Double_t mA[8][ 8];
  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<8; j++){
      mA[i][j] = 0;
    }
  }
  //mA[0][0] =  c;  mA[0][1] = s;
  //mA[1][0] = -s;  mA[1][1] = c;
  //mA[3][3] =  c;  mA[3][4] = s;
  //mA[4][3] = -s;  mA[4][4] = c;

  mA[0][0] =  trans.XX();  mA[0][1] = trans.XY();  mA[0][2] = trans.XZ();
  mA[1][0] =  trans.YX();  mA[1][1] = trans.YY();  mA[1][2] = trans.YZ();
  mA[2][0] =  trans.ZX();  mA[2][1] = trans.ZY();  mA[2][2] = trans.ZZ();

  mA[3][3] =  trans.XX();  mA[3][4] = trans.XY();  mA[3][5] = trans.XZ();
  mA[4][3] =  trans.YX();  mA[4][4] = trans.YY();  mA[4][5] = trans.YZ();
  mA[5][3] =  trans.ZX();  mA[5][4] = trans.ZY();  mA[5][5] = trans.ZZ();

  mA[6][6] = 1; mA[7][7] = 1;

  Double_t mAC[8][8];
  Double_t mAp[8];

  for( Int_t i=0; i<8; i++ ){
    mAp[i] = 0;
    for( Int_t k=0; k<8; k++){

      if(k==0)      mAp[i]+=mA[i][k] * ( Vtrack->GetX() - dx );
      else if(k==1) mAp[i]+=mA[i][k] * ( Vtrack->GetY() - dy );
      else if(k==2) mAp[i]+=mA[i][k] * ( Vtrack->GetZ() - dz );
      else if(k==3) mAp[i]+=mA[i][k] * ( Vtrack->Px() );
      else if(k==4) mAp[i]+=mA[i][k] * ( Vtrack->Py() );
      else if(k==5) mAp[i]+=mA[i][k] * ( Vtrack->Pz() );
      else if(k==6) mAp[i]+=mA[i][k] * ( Vtrack->E() );
      else if(k==7) mAp[i]+=mA[i][k] * 0.;
      else          mAp[i]+=mA[i][k] * 0.;
    }
  }

  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<8; j++ ){
      mAC[i][j] = 0;
      for( Int_t k=0; k<8; k++ ){
        mAC[i][j]+= mA[i][k] * cov[IJ(k,j)];
      }
    }
  }

  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<=i; j++ ){
      Double_t xx = 0;
      for( Int_t k=0; k<8; k++){
        xx+= mAC[i][k]*mA[j][k];
      }
      cov[IJ(i,j)] = xx;
    }
  }


  TParticle part(11,0,0,0,0,0,mAp[3],mAp[4],mAp[5],mAp[6], dx,dy,dz,0);

  AliESDtrack esdTrack(&part);
  return esdTrack;
//  AliVTrack* Vtrack_rot = dynamic_cast<AliVTrack*>(&esdTrack);

//  return Vtrack_rot;
}


