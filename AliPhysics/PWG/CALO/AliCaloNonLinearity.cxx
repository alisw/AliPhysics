/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Nicolas Schmidt, Friederike Bock                              *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling non linearity correction for
// calorimeter cluster energies
//---------------------------------------------
////////////////////////////////////////////////

#include "AliCaloNonLinearity.h"
#include "TF1.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include <vector>

class iostream;

ClassImp(AliCaloNonLinearity)


AliCaloNonLinearity::AliCaloNonLinearity(const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fCurrentMC(kNoMC)
{
    // Default constructor
}

AliCaloNonLinearity::~AliCaloNonLinearity() 
{
    // default deconstructor
}

// Return non linearity corrected cluster energy or in case of failure the value -100
Float_t AliCaloNonLinearity::GetCorrectedEnergy(Float_t clusterEnergy, Int_t isMC, Int_t switchNonLin, Int_t clusterType, MCSetEnum periodEnum)
{
  // Check if cluster energy is defined
  if (!clusterEnergy) {
    AliInfo("Cluster energy pointer null!\n");
    return -100;
  }

  Float_t energy = clusterEnergy;
  if (energy < 0.05) {
    // Clusters with less than 50 MeV or negative are not possible
    AliInfo(Form("Too Low Cluster energy!, E = %f < 0.05 GeV\n",energy));
    return -100;
  }

  // Obtain enum for period name string
//   fCurrentMC = FindEnumForMCSetString(periodName);
  fCurrentMC = periodEnum;
  AliInfo(Form("AliCaloNonLinearity:Period enum has been set to %o\n",fCurrentMC )) ;

  Bool_t fPeriodNameAvailable = kTRUE;

  switch(switchNonLin){

    // Standard NonLinearity -
    case 1:
      if( fClusterType == 1|| fClusterType == 3){
        // standard kPi0MCv5 for MC and kSDMv5 for data from Jason
        energy *= FunctionNL_kPi0MCv5(energy);
        if(isMC == 0) energy *= FunctionNL_kSDMv5(energy);
      } else if ( fClusterType == 2 ){
          // Nonlin from PHOS group only MC part
          if(isMC != 0) {
              if( fCurrentMC==k14j4 ){
                  energy *= FunctionNL_PHOSOnlyMC(energy, 1.008, 0.015, 0.4);
                  // for LHC13bc
              } else if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c || fCurrentMC == kPPb5T13P2HIJAdd){
                  energy *= FunctionNL_PHOSOnlyMC(energy, 1.0135, 0.018, 1.9);
              } else if(  // pp 5 TeV 2015
                  fCurrentMC == k16h3  || fCurrentMC == k16h8a || fCurrentMC == k16h8b || fCurrentMC == k16k3a  || fCurrentMC == k16k5a ||  fCurrentMC == k16k5b || fCurrentMC == k17e2 ||
                  // PbPb 5 TeV 2015
                  fCurrentMC == k16k3b ||
                  // pPb 5 TeV 2016
                  fCurrentMC == kPPb5T16EPOS || fCurrentMC == kPPb5T16DPMJet ||
                  // pPb 8 TeV 2016
                  fCurrentMC == k17f3a || fCurrentMC == k17f3b || fCurrentMC == k17f4a || fCurrentMC == k17f4b ||
                  // XeXe 5.44 TeV 2017
                  fCurrentMC == kXeXe5T17HIJING
              ){
                  energy *= FunctionNL_PHOSOnlyMC(energy, 1.012, -0.06, 0.7);
              }
          }
      }
      break;

    // kPi0MCv3 for MC and kTestBeamv3 for data
    case 2:
      if (fClusterType == 1|| fClusterType == 3){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
        else energy *= FunctionNL_kPi0MCv3(energy);
      }
      break;
    // kPi0MCv3 for MC and kTestBeamv2 for data
    case 3:
      if (fClusterType == 1|| fClusterType == 3){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
        else energy *= FunctionNL_kPi0MCv3(energy);
      }
      break;

    // kPi0MCv2 for MC and kTestBeamv3 for data
    case 4:
      if (fClusterType == 1|| fClusterType == 3){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
        else energy *= FunctionNL_kPi0MCv2(energy);
      }
      break;
    // kPi0MCv2 for MC and kTestBeamv2 for data
    case 5:
      if (fClusterType == 1|| fClusterType == 3){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
        else energy *= FunctionNL_kPi0MCv2(energy);
      }
      break;

    // kPi0MCv1 for MC and kTestBeamv3 for data
    case 6:
      if (fClusterType == 1|| fClusterType == 3){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
        else energy *= FunctionNL_kPi0MCv1(energy);
      }
      break;
    // kPi0MCv1 for MC and kTestBeamv2 for data
    case 7:
      if (fClusterType == 1|| fClusterType == 3){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
        else energy *= FunctionNL_kPi0MCv1(energy);
      }
      break;

    // kPi0MCv6 for MC and kSDMv6 for data
    case 8:
      if (fClusterType == 1|| fClusterType == 3){
        if(isMC == 0) energy *= FunctionNL_kSDMv6(energy);
        else energy *= FunctionNL_kPi0MCv6(energy);
      }
      break;
    // case 11 of the 8 TeV (LHC15h1 PYTHIA8) nonlinearity as a general case
    case 9:
      if(isMC>0){
        if (fClusterType == 1){
          energy /= FunctionNL_kSDM(energy, 0.96874*0.991*0.9958*0.999, -3.76064, -0.193181);
        }
      }
      break;

//----------------------------------------------------------------------------------------------------------

// *************** 10 + x **** default tender settings - pp

    // NonLinearity pp ConvCalo - only shifting MC - no timing cut
    case 11:
      label_case_11:
      if(isMC>0){
        // 8TeV LHC12x
        //pass1
        if( fCurrentMC==k14e2b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.983251, -3.44339, -1.70998);

        //pass2
        } else if( fCurrentMC == kPP8T12P2Pyt8 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.96874*0.991*0.9958*0.999, -3.76064, -0.193181);

        } else if( fCurrentMC == kPP8T12P2Pho ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.969703*0.989*0.9969*0.9991, -3.80387, -0.200546);

        } else if( fCurrentMC == kPP8T12P2JJ ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.974859*0.987*0.996, -3.85842, -0.405277);

        // 2.76TeV LHC11a/LHC13g
        } else if( fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.984889*0.995*0.9970, -3.65456, -1.12744);

        } else if(fCurrentMC==k12f1b){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.984384*0.995*0.9970, -3.30287, -1.48516);

        } else if( fCurrentMC==kPP2T11P4JJ || fCurrentMC==k15g1b || fCurrentMC==kPP2T13P1JJ || fCurrentMC==k15a3b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.981892*0.995*0.9970, -5.43438, -1.05468);

        // 7 TeV LHC10x
        } else if( fCurrentMC==k14j4 ){ //v3
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.974525*0.986*0.999, -4.00247, -0.453046) ;
            energy /= FunctionNL_kSDM(energy, 0.988038, -4.27667, -0.196969);
            energy /= FunctionNL_kSDM(energy, 0.997544, -4.5662, -0.459687);
          }
        // pp 5.02 TeV LHC15n
        // pass2
        } else if( fCurrentMC==k16h8a ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.969799, -4.11836, -0.293151);

        } else if( fCurrentMC==k16h8b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.969944, -4.02916, -0.366743);

        } else if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||  fCurrentMC==k17e2 ) {
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.972156, -4.10515, -0.381273);
            energy /= FunctionNL_kSDM(energy, 0.979999, -4.39136, -0.102332);
          }
          if(fClusterType==3 && fCurrentMC==k16k5a) {
            energy /= 0.9870110951;
            energy /= FunctionNL_kSDM(energy, 0.992345, -2.33772, -6.1127);
          }
          if(fClusterType==3 && fCurrentMC==k17e2) {
            energy /= FunctionNL_kSDM(energy, 0.986513, 0.430032, -10.99999);
            energy /= 0.9908118231;
          }
        } else if( fCurrentMC==k16k5b ) {
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.974173, -4.07732, -0.570223);
          if(fClusterType==3) {
            energy /= 0.9872826260;
            energy /= 0.9930726691;
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 12:
      label_case_12:
      if(isMC>0){
        // 8TeV LHC12x
        //pass1
        if( fCurrentMC==k14e2b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(2.0*energy, 0.967301, -3.1683, -0.653058);

          //pass2
        } else if( fCurrentMC == kPP8T12P2Pyt8 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.963379*0.9985*0.9992, -3.61217, -0.614043);

        } else if( fCurrentMC == kPP8T12P2Pho ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.96105*0.999*0.9996, -3.62239, -0.556256);

        } else if( fCurrentMC == kPP8T12P2JJ  ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.960596*0.999*0.999, -3.48444, -0.766862);

        // 2.76TeV LHC11a/LHC13g
        } else if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(2.0*energy, 0.966151*0.995*0.9981, -2.97974, -0.29463);

        } else if( fCurrentMC==k12f1b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(2.0*energy, 0.988814*0.995*0.9981, 0.335011, -4.30322);

        } else if( fCurrentMC==kPP2T11P4JJ || fCurrentMC==k15g1b || fCurrentMC==kPP2T13P1JJ || fCurrentMC==k15a3b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(2.0*energy, 0.979994*0.995*0.9981, -3.24431, -0.760205);

        // 7TeV LHC10x
        } else if(  fCurrentMC==k14j4 ){ //v3
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.962095*0.9991*0.9993, -3.63967, -0.747825) ;
            energy /= FunctionNL_kSDM(energy, 0.988922, -4.47811, -0.132757);
            energy /= FunctionNL_kSDM(energy, 0.99738, -4.82724, -0.281305);
          } else if(fClusterType==2){ //const fit
            energy *= 1.0111857903;
          }
        // 5 TeV LHC15n
        //pass2
        } else if( fCurrentMC==k16h8a ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.958994, -4.48233, -0.0314569);

        } else if( fCurrentMC==k16h8b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.960074, -3.31954, -1.14748);

        //pass3
        } else if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||  fCurrentMC==k17e2 ) {
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.965224, -3.04336, -1.85638);
          if(fClusterType==3 && (fCurrentMC==k16k5a ||  fCurrentMC==k17e2)) energy /= 0.9835764493;
        } else if( fCurrentMC==k16k5b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.960074, -3.31954, -1.14748);
          if(fClusterType==3) energy /= FunctionNL_kSDM(energy, 0.981191, -1.93399, -2.60859);


          //pp 13 TeV LHC16
        } else if ( fCurrentMC==kPP13T16P1Pyt8 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.957323, -3.55283, -0.57881);

        } else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    case 13:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_11;// goto previous case for shifting MC
      }
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 14:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_12;// goto previous case for shifting MC
      }
      break;

    // NonLinearity ConvCalo - kPi0MC + kSDM
    case 15:
      if (fClusterType == 1 || fClusterType == 3){
        // 8TeV LHC12x
        if ( fCurrentMC==k14e2b || fCurrentMC == kPP8T12P2Pyt8 || fCurrentMC == kPP8T12P2Pho  || fCurrentMC == k12pp8TeV || fCurrentMC == kPP8T12P2JJ ){
          energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04979, 1.3, 0.0967998, 219.381, 63.1604, 1.011);
          if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9846, -3.319, -2.033);

        // 2.76TeV LHC11a/LHC13g
        } else if ( fCurrentMC == k12f1a || fCurrentMC == k12i3 || fCurrentMC == k15g2 || fCurrentMC == k12f1b ||
                    fCurrentMC == kPP2T11P4JJ || fCurrentMC == k15g1b || fCurrentMC == kPP2T13P1JJ || fCurrentMC == k15a3b ||
                    fCurrentMC == k11pp2760GeV || fCurrentMC == k13pp2760GeV
                  ) {
          energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04123, 1.045, 0.0967998, 219.381, 63.1604, 1.014);
          if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9807*0.995*0.9970, -3.377, -0.8535);
        }
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity Calo - kPi0MC + kSDM
    case 16:
      if (fClusterType == 1 || fClusterType == 3){
        // 8TeV LHC12x
        if ( fCurrentMC==k14e2b  || fCurrentMC == kPP8T12P2Pyt8 || fCurrentMC == kPP8T12P2Pho  || fCurrentMC == k12pp8TeV || fCurrentMC == kPP8T12P2JJ ){
          energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06539, 1.121, 0.0967998, 219.381, 63.1604, 1.011);
          if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9676, -3.216, -0.6828);

        // 2.76TeV LHC11a/LHC13g
        } else if ( fCurrentMC == k12f1a || fCurrentMC == k12i3 || fCurrentMC == k15g2 || fCurrentMC == k12f1b ||
                    fCurrentMC == kPP2T11P4JJ || fCurrentMC == k15g1b || fCurrentMC == kPP2T13P1JJ || fCurrentMC == k15a3b ||
                    fCurrentMC == k11pp2760GeV || fCurrentMC == k13pp2760GeV
                  ) {
          energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06115, 0.9535, 0.0967998, 219.381, 63.1604, 1.013);
          if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9772*0.995*0.9981, -3.256, -0.4449);
        }
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // New PCM-EMC nonlinearity with energy squared
    case 17:
      if(isMC>0){
         if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||fCurrentMC==k16k5b ||  fCurrentMC==k17e2 ) {
          if(fClusterType==1){
            energy /= (FunctionNL_kSDM(energy, 0.945037*1.005, -3.42935, -0.384718));
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;

// *************** 20 + x **** modified tender Settings 1 - pp
    // NonLinearity pp ConvCalo - only shifting MC - no timing cut
    case 21:
      label_case_21:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0443938253, -0.0691830812, -0.1247555443, 1.1673716264, -0.1853095466, -0.0848801702) - 0.0055);
        } else if(fCurrentMC==k15g2){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.1716155406, -0.1962930603, -0.0193959829, 1.0336659741, -0.0467778485, -0.4407662248) - 0.0055);
        } else if(fCurrentMC==k12f1b){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0166321784, -0.0440799552, -0.2611899222, 1.0636538464, -0.0816662488, -0.2173961316) - 0.007);
        } else if( fCurrentMC==kPP2T11P4JJ || fCurrentMC==k15g1b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.1100193881, -0.1389194936, -0.0800000242, 1.1673716264, -0.1853095466, -0.0848801702) - 0.017);
        } else if( fCurrentMC==kPP2T13P1JJ || fCurrentMC==k15a3b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0520183153, -0.0806102847, -0.1450415920, 1.0336724056, -0.0467844121, -0.4406992764) - 0.016);
        // 8TeV
        } else if( fCurrentMC == kPP8T12P2Pyt8 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0654169768, -0.0935785719, -0.1137883054, 1.1814766150, -0.1980098061, -0.0854569214) - 0.0138);
        } else if( fCurrentMC == kPP8T12P2Pho ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0652493513, -0.0929276101, -0.1113762695, 1.1837801885, -0.1999914832, -0.0854569214) - 0.0145);
        } else if( fCurrentMC == kPP8T12P2JJ ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0489259285, -0.0759079646, -0.1239772934, 1.1835846739, -0.1998987993, -0.0854186691) - 0.014);
        // 7 TeV
        } else if( fCurrentMC == k14j4 ){ //v3
          if(fClusterType==1){
            energy /= (FunctionNL_DPOW(energy, 1.1082846035, -0.1369968318, -0.0800000002, 1.1850179319, -0.1999999950, -0.0863054172) - 0.015);
            energy /= FunctionNL_kSDM(energy, 0.988248, -4.26369, -0.208921) ;
            energy /= FunctionNL_kSDM(energy, 0.997359, -4.51031, -0.460041) ;
          }
        // 5 TeV LHC15n
        //pass2
        } else if( fCurrentMC==k16h8a ){
          if(fClusterType==1) energy /= (FunctionNL_DExp(energy, 0.9831956962, 1.2383793944, -3.2676359751, 1.0121710221, 0.6588125132, -3.1578818630));
        } else if( fCurrentMC==k16h8b ){
          if(fClusterType==1) energy /= (FunctionNL_DExp(energy, 0.9912139474, 0.3721971884, -3.6440765835, 1.0141024579, 0.5574244401, -3.1894624833));
        //pass3
        } else if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||  fCurrentMC==k17e2 ) {
          if(fClusterType==1){
            energy /= (FunctionNL_DPOW(energy, 1.1251817141, -0.8420328354, -0.0800008954, 0.9562653194, -0.6683378769, -0.1064755375));
            energy /= FunctionNL_kSDM(energy, 0.977706, -4.21058, -0.0938915) ;
          } else if( fClusterType==3) {
            if( fCurrentMC==k16k5a ) {
              energy /= (FunctionNL_DPOW(energy, 0.9943969544,-0.0181151588,-0.4999998851,1.0288066416,-0.0367913727,-0.4995137932));
              energy /= FunctionNL_DPOW(energy, 1.0055560859, -0.0213391278, -0.4999999991, 1.1047136553, -0.1141567995, -0.1573142879);
              energy /= FunctionNL_DPOW(energy, 1.0275381918, -0.0400165029, -0.4999999995, 1.0703233524, -0.0855441426, -0.2099590700);
            } else if( fCurrentMC==k17e2 ) {
              energy /= (FunctionNL_DPOW(energy, 0.9943969544,-0.0181151588,-0.4999998851,1.0288066416,-0.0367913727,-0.4995137932));
              energy /= FunctionNL_DPOW(energy, 1.0055560859, -0.0213391278, -0.4999999991, 1.1047136553, -0.1141567995, -0.1573142879);
              energy /= FunctionNL_DPOW(energy, 1.0275381918, -0.0400165029, -0.4999999995, 1.0703233524, -0.0855441426, -0.2099590700);
            }
          }

        } else if( fCurrentMC==k16k5b ){
          if(fClusterType==1) energy /= (FunctionNL_DExp(energy, 0.9842689920, 0.9150246921, -3.6796298486, 1.0113148506, 0.6876891951, -3.1672234730));
          if(fClusterType==3) {
            energy /= (FunctionNL_DPOW(energy, 1.1343351836,-0.1571288013,-0.0800000607,1.0288066416,-0.0367913727,-0.4995137932));
            energy /= FunctionNL_DPOW(energy, 1.1105555600, -0.1266067088, -0.0800000497, 1.1047136553, -0.1141567995, -0.1573142879);
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 22:
      label_case_22:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 0.9980625418, -0.0564782662, -0.5, 1.0383412435, -0.0851830429, -0.4999999996) - 0.00175);
        } else if( fCurrentMC==k15g2 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0795372569, -0.1347324732, -0.1630736190, 1.1614181498, -0.199995361, -0.1711378093) - 0.0035);
        } else if( fCurrentMC==k12f1b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0232969083, -0.090409434, -0.3592406513, 1.0383412435, -0.0851830429, -0.4999999996) + 0.0007);
        } else if( fCurrentMC==kPP2T11P4JJ || fCurrentMC==k15g1b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0106037132, -0.0748250591, -0.4999999996, 1.0383412435, -0.0851830429, -0.4999999996) - 0.014);
        } else if( fCurrentMC==kPP2T13P1JJ || fCurrentMC==k15a3b ) {
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0119417393, -0.0755250741, -0.4999999996, 1.1614181498, -0.1999995361, -0.1711378093) - 0.006);
        //8TeV
        } else if( fCurrentMC == kPP8T12P2Pyt8 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.1389201636, -0.1999994717, -0.1622237979, 1.1603460704, -0.1999999989, -0.2194447313) - 0.0025);
        } else if( fCurrentMC == kPP8T12P2Pho ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0105301622, -0.0732424689, -0.5000000000, 1.0689250170, -0.1082682369, -0.4388156470) - 0.001);
        } else if( fCurrentMC == kPP8T12P2JJ ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 0.9922456908, -0.0551212559, -0.5000000000, 1.0513459039, -0.0894163252, -0.5000000000) + 0.002);
        // 7 TeV
        } else if( fCurrentMC == k14j4 ){ //v3
          if(fClusterType==1){
            energy /= (FunctionNL_DPOW(energy, 1.0074002842, -0.0682543971, -0.4509341085, 1.1224162203, -0.1586806096, -0.2458351112) - 0.003) ;
            energy /= FunctionNL_kSDM(energy, 0.99598, -5.03134, -0.269278) ;
            energy /= FunctionNL_kSDM(energy, 0.997738, -4.91921, -0.377381) ;
          }
        // 5 TeV LHC15n
        //pass2
        } else if( fCurrentMC==k16h8a ){
          if(fClusterType==1) energy /= (FunctionNL_DExp(energy, 0.9747084556, 1.3652950049, -1.7832191813, 1.0039014622, 1.3657547071, -1.7852900827));
        } else if( fCurrentMC==k16h8b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0193460981, -0.0851635674, -0.4984580141, 1.0588985795, -0.0957023147, -0.4999999998));
        //pass3
        } else if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||  fCurrentMC==k17e2 ) {
          if(fClusterType==1){
            energy /= (FunctionNL_DExp(energy, 0.9762188425, 0.9120374996, -2.3012968797, 1.0049037083, 1.2643533472, -1.8927172439));
            energy /= (FunctionNL_kSDM(energy, 0.983808, -4.25003, -0.0977335)- 0.003);
          } else if(fClusterType==3){
            if(fCurrentMC==k17e2) energy /= 0.9825370234*0.9993152454;
            if(fCurrentMC==k16k5a) energy /= 0.9825370234*0.9993152454;
          }
        } else if( fCurrentMC==k16k5b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0193460981, -0.0851635674, -0.4984580141, 1.0588985795, -0.0957023147, -0.4999999998));
          if(fClusterType==3) energy /= (FunctionNL_DPOW(energy, 0.9629798154, -0.0178058455, -0.4999999880, 1.1467423891, -0.1999980199, -0.1753999427));

        //pp 13 TeV LHC16
        } else if ( fCurrentMC==kPP13T16P1Pyt8 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0187401756, -0.0857332791, -0.5000000000, 1.1585209386, -0.1999999989, -0.2646540338));
        } else fPeriodNameAvailable = kFALSE;
      }
      break;
    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    case 23:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_21;// goto previous case for shifting MC
      }
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 24:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_22;// goto previous case for shifting MC
      }
      break;

    // New PCM-EMC nonlinearity with energy squared
    case 27:
      if(isMC>0){
         if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||fCurrentMC==k16k5b ||  fCurrentMC==k17e2 ) {
          if(fClusterType==1){
            energy /= (FunctionNL_DPOW(energy, 1.1497456392, -0.1999999732, -0.0839303140, 1.1818406492, -0.1999998957, -0.1434322871) + 0.0055);
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;

// *************** 30 + x **** modified tender Settings 2 - pp

// *************** 40 + x **** default tender Settings - pPb
    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 41:
      label_case_41:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.967546, -3.57657, -0.233837) ; // with TM pt dep
            energy /= FunctionNL_kSDM(energy, 0.987513, -4.34641, -0.522125) ;
            energy /= 0.9935;
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /=  FunctionNL_kSDM(energy, 0.968868, -3.38407, -0.318188) ;
            energy /=  (FunctionNL_kSDM(energy, 0.987931, -4.13218, -0.583746)*0.9953479301) ;//with TM pt dep
          }
        } else if( fCurrentMC==kPPb5T16EPOS ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_kSDM(energy, 0.951944, -3.38177, -0.597868) ; //  2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1682295822, -0.1999999973, -0.2088780018, 1.1653083892, -0.1999999998, -0.2697014136 ); //  2018 03 22
          } else if (fClusterType==2){
            energy /= FunctionNL_kSDM(energy, 0.977985, -2.97438, -0.598613) ;
            energy /= 1.02231;
          }
        } else if( fCurrentMC==kPPb5T16DPMJet ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_kSDM(energy, 0.944223, -3.2442, -0.599841) ;   //  2018 02 20
            energy /= FunctionNL_kSDM(energy, 0.99851, -4.18038, -0.289851) ;   //  2018 03 22
          } else if (fClusterType==2){
            energy /= FunctionNL_kSDM(energy, 0.966329, -2.60954, -0.712271) ;
            energy /= 1.0105;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 42:
      label_case_42:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ) {
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.973301, -3.66136, -1.20116) ; //with TM pt dep
            energy /= (FunctionNL_kSDM(energy, 0.987611, -4.14227, -0.282541) * 1.0036264536 );
            energy /= 0.9935;
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.962047, -3.18433, -0.586904); //with TM pt dep
            energy /= FunctionNL_kSDM(energy, 0.990771, -4.29086, -0.27403);
          }
        } else if( fCurrentMC==kPPb5T16EPOS ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_kSDM(energy, 0.949402, -3.17052, -0.57999) ; // 2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1749094095, -0.1999999946, -0.2097130855, 1.1820209963, -0.1999999999, -0.1811167881 ); // 2018 03 22
          }
        } else if( fCurrentMC==kPPb5T16DPMJet ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_kSDM(energy, 0.950512, -2.68457, -0.989215) ; // 2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1720550395, -0.2, -0.1990753909, 1.1820209963, -0.1999999998, -0.1811167881 ); // 2018 03 22
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb ConvCalo  - kTestBeamv3 + shifting MC
    case 43:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_41;// goto previous case for shifting MC
      }
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 44:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_42;// goto previous case for shifting MC
      }
      break;
    // NonLinearity LHC13 pPb ConvCalo  - applying f^2
    case 45:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.967546, -3.57657, -0.233837) ;
            //apply again the same
            energy /= FunctionNL_kSDM(energy, 0.967546, -3.57657, -0.233837) ;
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /=  FunctionNL_kSDM(energy, 0.968868, -3.38407, -0.318188) ;
            //apply again the same
            energy /=  FunctionNL_kSDM(energy, 0.968868, -3.38407, -0.318188) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

// *************** 50 + x **** modified tender Settings 1 - pPb
    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 51:
      label_case_51:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ) {
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 0.9910691195, 0.4901455923, -3.6647921806, 1.0255088817, 0.3070452373, -2.9149185308); //with TM pt dep
            energy /= FunctionNL_kSDM(energy, 0.989111, -4.26219, -0.819192);
            energy /= 0.9935;
          } else if(fClusterType==2){
            energy /= ( 0.994914734 * 0.9964 ); // additional factors
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 0.9978241421, 0.2054669115, -3.7888984452, 1.0255088817, 0.3070452373, -2.9149185308) ; //with TM pt dep
            energy /= (FunctionNL_kSDM(energy, 0.986673, -4.14594, -0.450765)* 0.9953727823);
          } else if(fClusterType==2){
            energy /= ( 0.993485*0.9971126333 );
          }
        } else if( fCurrentMC==kPPb5T16EPOS ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_DExp(energy, 0.9772393830, 0.9651600903, -2.8485741777, 1.0436698408, 0.4584792411, -2.3634185342); // 2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1771657563, -0.1999999970, -0.1955919316, 1.1820209963, -0.1999999999, -0.1811167881 ) ; // 2018 03 22
          } else if(fClusterType==2){
            energy /= (0.949117*1.02231) ; //first iteration with constant
          }
        } else if( fCurrentMC==kPPb5T16DPMJet ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_DExp(energy, 0.9822883129, 0.4885659420, -2.8847511641, 1.0461791627, 0.4463927828, -2.2994619341); // 2018 02 20
            energy /= FunctionNL_kSDM(energy, 0.992184, -4.22097, -0.561982); // 2018 03 22
          } else if(fClusterType==2){
            energy /= ( 0.997*0.9965200155 ); // additional factors
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 52:
      label_case_52:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ) {
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 0.9795532189, 0.8578583955, -2.3447892540, 1.0165873637, 0.6999387334, -2.1324782465) ;//with TM pt dep
            energy /= (FunctionNL_kSDM(energy, 0.990609, -4.37834, -0.304314) * 1.0040232773) ;
            energy /= 0.9935;
          } else if(fClusterType==2) {
            energy /= (FunctionNL_DExp(energy, 1.0154938040, 0.3062978125, -3.9089772679, 1.0061692542, 513.7621552761, -3566.4426936867 ) * 0.996512);
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 0.9746342307, 0.9576270870, -2.5098585110, 1.0165871862, 0.6999571530, -2.1324658480) ; //with TM pt dep
            energy /= FunctionNL_kSDM(energy, 0.993562, -4.52817, -0.366368) ;
          } else if(fClusterType==2){
            energy /= (FunctionNL_DPOW(energy, 1.0154784875, -0.0161589457, -0.4999999976, 1.0086650887, -0.0010000001, -0.0800000139 ) * 0.9983468115 );
          }
       } else if( fCurrentMC==kPPb5T16EPOS ) {
          if(fClusterType==1|| fClusterType==3 ){
            energy /= FunctionNL_DExp(energy, 0.9707112053, 1.4050445187, -2.0357906356, 1.0241095707, 0.9217457498, -1.9020815528) ;//2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1753551048, -0.1999999981, -0.2060941701, 1.1822845437, -0.2, -0.1793151582 ) ; // 2018 03 22
          } else if(fClusterType==2) {
            energy /= (FunctionNL_DExp(energy, 1.0154938040, 0.3062978125, -3.9089772679, 1.0061692542, 513.7621552761, -3566.4426936867 ) * 0.996512);
          }
       } else if( fCurrentMC==kPPb5T16DPMJet ) {
          if(fClusterType==1|| fClusterType==3 ){
            energy /= FunctionNL_DExp(energy, 0.9706706146, 0.9781531357, -2.5633710383, 1.0355397924, 0.6750800461, -1.9817285526) ;//2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1717517490, -0.1999999942, -0.2126460833, 1.1820209963, -0.1999999999, -0.1811167881 ) ; // 2018 03 22
          } else if(fClusterType==2) {
            energy /= (FunctionNL_DExp(energy, 1.0154938040, 0.3062978125, -3.9089772679, 1.0061692542, 513.7621552761, -3566.4426936867 ) * 0.996512);
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;
    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 53:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_51;// goto previous case for shifting MC
      }
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 54:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_52;// goto previous case for shifting MC
      }
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 /kMCTestBeam+shift from PCM-EMC
    case 55:
      if (fClusterType == 1 || fClusterType == 3){
        if(isMC == 0) {
          energy *= FunctionNL_kTestBeamv3(energy);
        } else{
          energy *= FunctionNL_kPi0MCv3(energy);
          if( fCurrentMC==kPPb5T16EPOS )          energy /= 0.9854013683;
          else if ( fCurrentMC==kPPb5T16DPMJet )  energy /= 0.9818891524;
        }

      }
      break;
    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 /kMCTestBeam+shift from EMC
    case 56:
      if (fClusterType == 1 || fClusterType == 3){
        if(isMC == 0) {
          energy *= FunctionNL_kTestBeamv3(energy);
        } else{
          energy *= FunctionNL_kPi0MCv3(energy);
          if( fCurrentMC==kPPb5T16EPOS )          energy /= 0.9887044419;
          else if ( fCurrentMC==kPPb5T16DPMJet )  energy /= 0.9891917142;
        }
      }
      break;

// *************** 60 + x **** modified tender Settings 2 - pPb


// *************** 70 + x **** default tender Settings - PbPb

    // NonLinearity LHC11h - PbPb 2.76TeV - 0-10% centrality
    case 71:
      if(isMC>0){
        if( fCurrentMC==k14a1 ){
          if (fClusterType == 1 || fClusterType == 3){
            energy /= 0.972607; //
          }
        } else if( fCurrentMC==kXeXe5T17HIJING ){
          if (fClusterType == 1 ){
            energy /= FunctionNL_kSDM(energy, 0.973646, -0.901289, -4.32682) ; // based on peripheral XeXe
          } else if (fClusterType == 2 ){
            energy /= (0.9599764493*0.9873); // based on peripheral XeXe
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
    // NonLinearity LHC11h - PbPb 2.76TeV - 20-50% centrality
    case 72:
      if(isMC>0){
        if( fCurrentMC==k14a1 ){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 1.00926, -2.42107, -1.60995);
          }
        } else if( fCurrentMC==kXeXe5T17HIJING ){
          if (fClusterType == 1 ){
            energy /= 0.9737701536; // based on peripheral XeXe
          } else if (fClusterType == 2 ){
            energy /= (0.9350697962*1.01); // based on semi-central XeXe
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

    // *************** 80 + x **** modified tender Settings 1 - PbPb

    // NonLinearity LHC15o PbPb ConvCalo  - only shifting MC
    case 81:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (FunctionNL_kSDM(energy, 0.95597, -3.09059, -0.702889)*1.008) ;
          }
        } else if( fCurrentMC==kXeXe5T17HIJING ){
          if (fClusterType == 1 ){
            energy /= (FunctionNL_DPOW(energy, 1.0547527663, -0.0927180446, -0.0800012482, 1.0254208020, -0.0345156682, -0.4999999199)); // based on peripheral XeXe
          } else if (fClusterType == 2 ){
            energy /= (0.9764119296*0.9794); // based on peripheral XeXe
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

    // NonLinearity LHC15o PbPb Calo  - only shifting MC
    case 82:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.095, -0.175739, 0.00776757) ;
          }
        } else if( fCurrentMC==kXeXe5T17HIJING ){
          if (fClusterType == 1 ){
            energy /= (FunctionNL_DPOW(energy, 1.1223479533, -0.1999999659, -0.1954398178, 1.0373041075, -0.0913333671, -0.5000000000)); // based on peripheral XeXe
          } else if (fClusterType == 2 ){
            energy /= (FunctionNL_DExp(energy, 0.9840385879, 0.4801589926, -2.8482099501, 1.0214220397, 5.8987542970, -12.5701079799)*1.0148) ; // based on  semi-central XeXe
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - MB
    case 83:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.955001) ;
            energy /= FunctionNL_DExp(energy, 1.0380275426, 0.7534354400, -2.2110408210, 1.0408879042, 0.4399353376, -2.9554918759) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - 0-10%
    case 84:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.954566) ;
            energy /= FunctionNL_DExp(energy, 1.0548582854, 1.5096237243, -1.6079078305, 1.0538380642, 124049.7, -38409.5) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - 10-20%
    case 85:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.9365) ;
            energy /= FunctionNL_DExp(energy, 1.0380275426, 0.7534354400, -2.2110408210, 1.0408879042, 0.4399353376, -2.9554918759) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - 20-50%
    case 86:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.948553) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - 50-90%
    case 87:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.95306) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

// *************** 90 + x **** modified tender Settings 2 - PbPb

      // NonLinearity LHC15o PbPb ConvCalo  - only shifting MC
    case 91:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 1.0026971373, -0.0320283624, -0.4999999953, 1.0750656618, -0.0855019990, -0.4571523301);
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

      // NonLinearity LHC15o PbPb Calo  - only shifting MC
    case 92:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 1.0541217488, -0.1111428177, -0.4999999983, 1.0782958817, -0.0706389211, -0.4999999959);
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;


//----------------------------------------------------------------------------------------------------------

    default:
      AliFatal(Form("NonLinearity correction not defined for cut: '%d' ! Returning...",switchNonLin));
      return -100;

  }

  if(!fPeriodNameAvailable){
    AliFatal(Form("NonLinearity correction not defined for fCurrentMC: '%o'! Please check nonlin switch (%d) as well as function AliCaloNonLinearity::GetCorrectedEnergy. Correction failed, returning...",fCurrentMC,switchNonLin));
    return -100;
  }

  // Return nonlinearity corrected cluster energy
  return energy;
}



//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kPi0MC(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6){
  return ( p6 / ( p0 * ( 1. / ( 1. + p1 * exp( -e / p2 ) ) * 1. / ( 1. + p3 * exp( ( e - p4 ) / p5 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2){
  return ( p0 + exp( p1 + ( p2 * e ) ) );
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_DPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5){
  Float_t ret = 1;
  if ((p3 +  p4 * TMath::Power(e,p5 ) ) != 0)
    ret = ( (p0 +  p1 * TMath::Power(e,p2 ) )/(p3 +  p4 * TMath::Power(e,p5 ) ) );
  if (ret != 0.)
    return ret;
  else
    return 1.;
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_DExp(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5){
  Float_t ret = 1;
  if ( (p3 - TMath::Exp(-p4*e+p5) ) != 0)
    ret = ( (p0 - TMath::Exp(-p1*e+p2) )/(p3 - TMath::Exp(-p4*e+p5) ) );
  if (ret != 0.)
    return ret;
  else
    return 1.;
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_PHOSOnlyMC(Float_t e, Float_t p0, Float_t p1, Float_t p2){
  return p0*(1+p1/(1.+e*e/p2/p2)) ;
}



//************************************************************************
// predefined functions:
//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kPi0MCv1(Float_t e){
  return ( 1.014 * exp( 0.03329 / e ) ) + ( ( -0.3853 / ( 0.5423 * 2. * TMath::Pi() ) * exp( -( e + 0.4335 ) * ( e + 0.4335 ) / (2. * 0.5423 * 0.5423 ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kPi0MCv2(Float_t e){
  return ( 0.311111 / TMath::Power( e - 0.571666, 0.567995 ) + 1 );
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kPi0MCv3(Float_t e){
  return ( 1.0 / ( 0.981039 * ( 1. / ( 1. + 0.113508 * exp( -e / 1.00173 ) ) * 1. / ( 1. + 0.0967998 * exp( ( e - 219.381 ) / 63.1604 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kPi0MCv5(Float_t e){
  return ( 1.01286 / ( 1.0 * ( 1. / ( 1. + 0.0664778 * exp( -e / 1.57 ) ) * 1. / ( 1. + 0.0967998 * exp( ( e - 219.381 ) / 63.1604 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kPi0MCv6(Float_t e){
  return ( 1.00437 / ( 1.0 * ( 1. / ( 1. + 0.0797873 * exp( -e / 1.68322 ) ) * 1. / ( 1. + 0.0806098 * exp( ( e - 244.586 ) / 116.938 ) ) ) ) );
}

// only shifting data, to be used with kPi0MCv5 before
//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kSDMv5(Float_t e){
  return ( 0.964 + exp( -3.132 + ( -0.435 * 2.0 * e ) ) );
}

// be careful: different definition than kSDMv5
//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kSDMv6(Float_t e){
  return ( 0.987054 / ( 1.0 * ( 1. / ( 1. + 0.237767 * exp( -e / 0.651203 ) ) * 1. / ( 1. + 0.183741 * exp( ( e - 155.427 ) / 17.0335 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kTestBeamv2(Float_t e){
  return ( 0.968 / ( 0.983504 *( 1. / ( 1. + 0.210106 * exp( -e / 0.897274 ) ) * 1. / ( 1. + 0.0829064 * exp( ( e - 152.299 ) / 31.5028 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionNL_kTestBeamv3(Float_t e){
  return ( 0.9615 / ( 0.976941 *( 1. / ( 1. + 0.162310 * exp( -e / 1.08689 ) ) * 1. / ( 1. + 0.0819592 * exp( ( e - 152.338 ) / 30.9594 ) ) ) ) );
}

//************************************************************************
//________________________________________________________________________
Float_t AliCaloNonLinearity::FunctionM02(Float_t E, Float_t a, Float_t b, Float_t c, Float_t d, Float_t e){
  return ( exp( a+ b*E ) + c + d*E + e/E);
}

// returns the period enumerator for a given period string 
// can be used to obtain the enum for GetCorrectedEnergy
//________________________________________________________________________
AliCaloNonLinearity::MCSetEnum AliCaloNonLinearity::FindEnumForMCSetString(TString namePeriod){
  // pp 7TeV MB MCs pass 4
  if(       namePeriod.CompareTo("LHC14j4b")==0 ||
            namePeriod.CompareTo("LHC14j4c")==0 ||
            namePeriod.CompareTo("LHC14j4d")==0 ||
            namePeriod.CompareTo("LHC14j4e")==0 ||
            namePeriod.CompareTo("LHC14j4f")==0)        return k14j4;

  // pp 2.76 TeV LHC11a anch MC's pass 4
  else if(  namePeriod.CompareTo("LHC12f1a")==0)        return k12f1a;
  else if(  namePeriod.CompareTo("LHC12f1b")==0)        return k12f1b;
  else if(  namePeriod.CompareTo("LHC12i3")==0)         return k12i3;
  else if(  namePeriod.CompareTo("LHC15g1a")==0)        return kPP2T11P4JJ;
  else if(  namePeriod.CompareTo("LHC15g1b")==0)        return k15g1b;

  // PbPb 2.76TeV 2011 MC add sig
  else if(  namePeriod.CompareTo("LHC14a1a")==0 ||
            namePeriod.CompareTo("LHC14a1b")==0 ||
            namePeriod.CompareTo("LHC14a1c")==0)        return k14a1;

    // pp 8 TeV MC MCs
  // pass 1
  else if(  namePeriod.CompareTo("LHC14e2b")==0)        return k14e2b;
  // pass 2
  else if(  namePeriod.CompareTo("LHC12P2Pyt8")==0 ||
            namePeriod.CompareTo("LHC15h1")==0 ||
            namePeriod.CompareTo("LHC15h1a1")==0 ||
            namePeriod.CompareTo("LHC15h1b")==0 ||
            namePeriod.CompareTo("LHC15h1c")==0 ||
            namePeriod.CompareTo("LHC15h1d")==0 ||
            namePeriod.CompareTo("LHC15h1f")==0 ||
            namePeriod.CompareTo("LHC15h1g")==0 ||
            namePeriod.CompareTo("LHC15h1h")==0 ||
            namePeriod.CompareTo("LHC15h1i")==0)        return kPP8T12P2Pyt8;
  else if(  namePeriod.CompareTo("LHC12P2Pho")==0 ||
            namePeriod.CompareTo("LHC15h2")==0 ||
            namePeriod.CompareTo("LHC15h2a")==0 ||
            namePeriod.CompareTo("LHC15h2b")==0 ||
            namePeriod.CompareTo("LHC15h2c")==0 ||
            namePeriod.CompareTo("LHC15h2d")==0 ||
            namePeriod.CompareTo("LHC15h2f")==0 ||
            namePeriod.CompareTo("LHC15h2g")==0 ||
            namePeriod.CompareTo("LHC15h2h")==0 ||
            namePeriod.CompareTo("LHC15h2i")==0)        return kPP8T12P2Pho;
  // pp 8 TeV JJ MC pass 2
  else if ( namePeriod.CompareTo("LHC12P2JJ") == 0 ||
            namePeriod.CompareTo("LHC16c2") == 0 ||
            namePeriod.CompareTo("LHC16c2_plus") == 0 ) return kPP8T12P2JJ;

  // pPb 5 TeV 2013 MC pass 2
  else if(  namePeriod.Contains("LHC13b2_efix"))        return kPPb5T13P2DPMJet;
  else if(  namePeriod.Contains("LHC13e7"))             return kPPb5T13P2HIJAdd;
  // pPb 5 TeV 2013 MC JJ EMC enhanced
  else if ( namePeriod.CompareTo("LHC16c3a") == 0 ||
            namePeriod.CompareTo("LHC16c3a2") == 0 )    return k16c3a;
  else if ( namePeriod.CompareTo("LHC16c3b") == 0 ||
            namePeriod.CompareTo("LHC16c3b2") == 0 )    return k16c3b;
  // pPb 5 TeV 2013 MC GJ
  else if ( namePeriod.CompareTo("LHC16c3c") == 0 ||
            namePeriod.CompareTo("LHC16c3c2") == 0 )    return k16c3c;

  // pp 2.76 TeV LHC13g anch MC's pass 4
  else if(  namePeriod.CompareTo("LHC15g2")==0)         return k15g2;
  else if(  namePeriod.CompareTo("LHC15a3a")==0 ||
            namePeriod.CompareTo("LHC15a3a_plus")==0)   return kPP2T13P1JJ;
  else if(  namePeriod.CompareTo("LHC15a3b")==0)        return k15a3b;

  // pp 13 TeV 2015 MB pass 2
  else if ( namePeriod.CompareTo("LHC15P2Pyt8") == 0 ||
            namePeriod.CompareTo("LHC17i4") == 0 ||
            namePeriod.CompareTo("LHC17i4_2") == 0 ||
            namePeriod.CompareTo("LHC17g7") == 0 )      return kPP13T15P2Pyt8;
  else if ( namePeriod.CompareTo("LHC15P2EPos") == 0 ||
            namePeriod.CompareTo("LHC16d3") == 0 )      return kPP13T15P2EPOS;
  // pp 13 TeV 2015 HF prod
  else if ( namePeriod.CompareTo("LHC15k5a") == 0 ||
            namePeriod.CompareTo("LHC15k5b") == 0 ||
            namePeriod.CompareTo("LHC15k5c") == 0 ||
            namePeriod.CompareTo("LHC15k5a2") == 0 ||
            namePeriod.CompareTo("LHC15k5b2") == 0 ||
            namePeriod.CompareTo("LHC15k5c2") == 0 )    return k15k5;

  // pp 5 TeV 2015 MB MC pass 2
  else if ( namePeriod.CompareTo("LHC16h8a") == 0 )     return k16h8a;
  else if ( namePeriod.CompareTo("LHC16h8b") == 0 )     return k16h8b;
  else if ( namePeriod.CompareTo("LHC16k3a") == 0 ||                    // special pileup prods
            namePeriod.CompareTo("LHC16k3a2") == 0 )    return k16k3a;
  // pp 5 TeV 2015 MB MC pass 3
  else if ( namePeriod.CompareTo("LHC16k5a") == 0  )    return k16k5a;
  else if ( namePeriod.CompareTo("LHC16k5b") == 0  )    return k16k5b;
  // pp 5 TeV 2015 MB MC pass 4
  else if ( namePeriod.CompareTo("LHC17e2") == 0  )     return k17e2;
  // pp 5 TeV 2015 JJ pass 3
  else if ( namePeriod.CompareTo("LHC16h3") == 0 )      return k16h3;
  // pp 5 TeV 2017 MB MC pass 1
  else if ( namePeriod.CompareTo("LHC17l4b") == 0 ||
            namePeriod.CompareTo("LHC17l4b_fast") == 0 ||
            namePeriod.CompareTo("LHC17l4b_cent") == 0 ||
            namePeriod.CompareTo("LHC17l4b_cent_woSDD") == 0)     return k17l4b;
  else if ( namePeriod.CompareTo("LHC17l3b") == 0 ||
            namePeriod.CompareTo("LHC17l3b_fast") == 0 ||
            namePeriod.CompareTo("LHC17l3b_cent") == 0 ||
            namePeriod.CompareTo("LHC17l3b_cent_woSDD") == 0)     return k17l3b;
  // PbPb 5TeV 2015 MB prods
  else if ( namePeriod.CompareTo("LHC16g1") == 0 ||
            namePeriod.CompareTo("LHC16g1a") == 0 ||
            namePeriod.CompareTo("LHC16g1b") == 0 ||
            namePeriod.CompareTo("LHC16g1c") == 0 ||
            namePeriod.CompareTo("LHC16h4") == 0)       return kPbPb5T15HIJING;
  else if ( namePeriod.CompareTo("LHC16k3b") == 0 ||                    // special pileup prods
            namePeriod.CompareTo("LHC16k3b2") == 0 )    return k16k3b;

  // pp 13 TeV 2016 MB prod
  else if ( namePeriod.CompareTo("LHC16P1Pyt8") == 0 ||
            namePeriod.CompareTo("LHC17f6") == 0 ||
            namePeriod.CompareTo("LHC17d17") == 0 ||
            namePeriod.CompareTo("LHC17f5") == 0 ||
            namePeriod.CompareTo("LHC17d3") == 0 ||
            namePeriod.CompareTo("LHC17e5") == 0 ||
            namePeriod.CompareTo("LHC17d20a1") == 0 ||
            namePeriod.CompareTo("LHC17d20a1_extra") == 0 ||
            namePeriod.CompareTo("LHC17d20a2") == 0 ||
            namePeriod.CompareTo("LHC17d20a2_extra") == 0 ||
            namePeriod.CompareTo("LHC17d16") == 0 ||
            namePeriod.CompareTo("LHC17d18") == 0 ||
            namePeriod.CompareTo("LHC17f9") == 0 ||
            namePeriod.CompareTo("LHC17f9_test") == 0 ) return kPP13T16P1Pyt8;
  else if ( namePeriod.CompareTo("LHC16P1Pyt8LowB") == 0 ||
            namePeriod.CompareTo("LHC17d1") == 0 )      return kPP13T16P1Pyt8LowB;
  else if ( namePeriod.CompareTo("LHC16P1EPOS") == 0 ||
            namePeriod.CompareTo("LHC17d20b1") == 0 ||
            namePeriod.CompareTo("LHC17d20b2") == 0 )   return kPP13T16P1EPOS;
  // pp 13 TeV 2016 JJ prod
  else if ( namePeriod.CompareTo("LHC16P1JJ") == 0 ||
            namePeriod.CompareTo("LHC17f8a") == 0 ||
            namePeriod.CompareTo("LHC17f8c") == 0 ||
            namePeriod.CompareTo("LHC17f8d") == 0 ||
            namePeriod.CompareTo("LHC17f8e") == 0 ||
            namePeriod.CompareTo("LHC17f8f") == 0 ||
            namePeriod.CompareTo("LHC17f8g") == 0 ||
            namePeriod.CompareTo("LHC17f8h") == 0 ||
            namePeriod.CompareTo("LHC17f8i") == 0 ||
            namePeriod.CompareTo("LHC17f8j") == 0 ||
            namePeriod.CompareTo("LHC17f8k") == 0 )     return kPP13T16P1JJ;
  else if ( namePeriod.CompareTo("LHC16P1JJLowB") == 0 ||
            namePeriod.CompareTo("LHC17f8b") == 0  )    return kPP13T16P1JJLowB;
  // pp 13 TeV 2016 HF prods
  else if ( namePeriod.CompareTo("LHC17h8a") == 0 )     return k17h8a;
  else if ( namePeriod.CompareTo("LHC17h8b") == 0 )     return k17h8b;
  else if ( namePeriod.CompareTo("LHC17h8c") == 0 )     return k17h8c;
  else if ( namePeriod.CompareTo("LHC17c3b1") == 0 )    return k17c3b1;
  else if ( namePeriod.CompareTo("LHC17c3a1") == 0 )    return k17c3a1;
  else if ( namePeriod.CompareTo("LHC17c3b2") == 0 )    return k17c3b2;
  else if ( namePeriod.CompareTo("LHC17c3a2") == 0 )    return k17c3a2;

  // pPb 5 TeV 2016 MB MC
  else if ( namePeriod.CompareTo("LHC17f2a") == 0  )    return kPPb5T16EPOS;
  else if ( namePeriod.CompareTo("LHC17f2b") == 0  )    return kPPb5T16DPMJet;
  // pPb 5 TeV 2016 JJ MC
  else if ( namePeriod.CompareTo("LHC17g8a") == 0  )    return k17g8a;
  // pPb 5 TeV 2016 HF prod
  else if ( namePeriod.CompareTo("LHC17d2a") == 0 )     return k17d2a;
  else if ( namePeriod.CompareTo("LHC17d2b") == 0 )     return k17d2b;

  // pPb 8 TeV 2016 MB MC
  else if ( namePeriod.CompareTo("LHC17f3a") == 0  )    return k17f3a;
  else if ( namePeriod.CompareTo("LHC17f3b") == 0  )    return k17f3b;
  else if ( namePeriod.CompareTo("LHC17f4a") == 0  )    return k17f4a;
  else if ( namePeriod.CompareTo("LHC17f4b") == 0  )    return k17f4b;
  // pPb 8 TeV 2016 JJ MC
  else if ( namePeriod.CompareTo("LHC17g8b") == 0  )    return k17g8b;
  else if ( namePeriod.CompareTo("LHC17g8c") == 0  )    return k17g8c;

  //pp 13 TeV LHC17
  else if ( namePeriod.CompareTo("LHC17k1") ==0 )       return k17k1; // HF low B
  else if ( namePeriod.CompareTo("LHC17k4") ==0 ||
            namePeriod.CompareTo("LHC17h11") ==0 ||
            namePeriod.CompareTo("LHC17h1") == 0 ||
            namePeriod.CompareTo("LHC17l5") == 0 )      return kPP13T17P1Pyt8;
  else if ( namePeriod.CompareTo("LHC17h7b") ==0 )      return kPP13T17P1Pho;
  else if ( namePeriod.CompareTo("LHC17h7a") ==0 )      return kPP13T17P1Pyt6;

  else if ( namePeriod.CompareTo("LHC17j5a") ==0 ||
            namePeriod.CompareTo("LHC17j5b") ==0 ||
            namePeriod.CompareTo("LHC17j5c") ==0 ||
            namePeriod.CompareTo("LHC17j5d") ==0 ||
            namePeriod.CompareTo("LHC17j5e") ==0 )      return kPP13T17P1Pyt8Str;
  else if ( namePeriod.CompareTo("LHC17h3") == 0 )      return kPP13T17P1Pyt8LowB;
  // XeXe 5.44 TeV 2017 MB MC
  else if ( namePeriod.CompareTo("LHC17j7") == 0 )      return kXeXe5T17HIJING;


  // data starts here
  else if ( namePeriod.CompareTo("LHC10b") == 0 ||
            namePeriod.CompareTo("LHC10c") == 0 ||
            namePeriod.CompareTo("LHC10d") == 0 ||
            namePeriod.CompareTo("LHC10e") == 0 ||
            namePeriod.CompareTo("LHC10f") == 0 ||
            namePeriod.CompareTo("LHC10g") == 0 )       return k10pp7TeV;
  else if ( namePeriod.CompareTo("LHC10h") == 0 )       return k10PbPb2760GeV;
  else if ( namePeriod.CompareTo("LHC11a") == 0 )       return k11pp2760GeV;
  else if ( namePeriod.CompareTo("LHC11b") == 0 ||
            namePeriod.CompareTo("LHC11c") == 0 ||
            namePeriod.CompareTo("LHC11d") == 0 ||
            namePeriod.CompareTo("LHC11e") == 0 ||
            namePeriod.CompareTo("LHC11f") == 0 ||
            namePeriod.CompareTo("LHC11g") == 0 )       return k11pp7TeV;
  else if ( namePeriod.CompareTo("LHC11h") == 0 )       return k11PbPb2760GeV;
  else if ( namePeriod.CompareTo("LHC12a") == 0 ||
            namePeriod.CompareTo("LHC12b") == 0 ||
            namePeriod.CompareTo("LHC12c") == 0 ||
            namePeriod.CompareTo("LHC12d") == 0 ||
            namePeriod.CompareTo("LHC12e") == 0 ||
            namePeriod.CompareTo("LHC12f") == 0 ||
            namePeriod.CompareTo("LHC12g") == 0 ||
            namePeriod.CompareTo("LHC12h") == 0 ||
            namePeriod.CompareTo("LHC12i") == 0 )       return k12pp8TeV;
  else if ( namePeriod.CompareTo("LHC13b") == 0 ||
            namePeriod.CompareTo("LHC13c") == 0 ||
            namePeriod.CompareTo("LHC13d") == 0 ||
            namePeriod.CompareTo("LHC13e") == 0 ||
            namePeriod.CompareTo("LHC13f") == 0 )       return k13pPb5023GeV;
  else if ( namePeriod.CompareTo("LHC13g") == 0 )       return k13pp2760GeV;
  else if ( namePeriod.CompareTo("LHC15fm") == 0 ||
            namePeriod.CompareTo("LHC15f") == 0 ||
            namePeriod.CompareTo("LHC15g") == 0 ||
            namePeriod.CompareTo("LHC15h") == 0 ||
            namePeriod.CompareTo("LHC15i") == 0 ||
            namePeriod.CompareTo("LHC15j") == 0 ||
            namePeriod.CompareTo("LHC15k") == 0 ||
            namePeriod.CompareTo("LHC15l") == 0 ||
            namePeriod.CompareTo("LHC15m") == 0 )       return k15pp13TeV;
  else if ( namePeriod.CompareTo("LHC15n") == 0 )       return k15pp5TeV;
  else if ( namePeriod.CompareTo("LHC15o") == 0 )       return k15PbPb5TeV;
  else if ( namePeriod.CompareTo("LHC16f") == 0 )       return k16pp13TeVLow;
  else if ( namePeriod.CompareTo("LHC16dp") == 0 ||
            namePeriod.CompareTo("LHC16d") == 0 ||
            namePeriod.CompareTo("LHC16e") == 0 ||
            namePeriod.CompareTo("LHC16g") == 0 ||
            namePeriod.CompareTo("LHC16h") == 0 ||
            namePeriod.CompareTo("LHC16i") == 0 ||
            namePeriod.CompareTo("LHC16j") == 0 ||
            namePeriod.CompareTo("LHC16k") == 0 ||
            namePeriod.CompareTo("LHC16l") == 0 ||
            namePeriod.CompareTo("LHC16m") == 0 ||
            namePeriod.CompareTo("LHC16n") == 0 ||
            namePeriod.CompareTo("LHC16o") == 0 ||
            namePeriod.CompareTo("LHC16p") == 0 )       return k16pp13TeV;
  else if ( namePeriod.CompareTo("LHC16qt") == 0 ||
            namePeriod.CompareTo("LHC16q") == 0 ||
            namePeriod.CompareTo("LHC16t") == 0 )       return k16pPb5023GeV;
  else if ( namePeriod.CompareTo("LHC16rs") == 0 ||
            namePeriod.CompareTo("LHC16r") == 0 ||
            namePeriod.CompareTo("LHC16s") == 0 )       return k16pPb8TeV;
  else if ( namePeriod.CompareTo("LHC17cr") == 0 ||
            namePeriod.CompareTo("LHC17c") == 0 ||
            namePeriod.CompareTo("LHC17d") == 0 ||
            namePeriod.CompareTo("LHC17e") == 0 ||
            namePeriod.CompareTo("LHC17f") == 0 ||
            namePeriod.CompareTo("LHC17h") == 0 ||
            namePeriod.CompareTo("LHC17i") == 0 ||
            namePeriod.CompareTo("LHC17j") == 0 ||
            namePeriod.CompareTo("LHC17k") == 0 ||
            namePeriod.CompareTo("LHC17l") == 0 ||
            namePeriod.CompareTo("LHC17m") == 0 ||
            namePeriod.CompareTo("LHC17o") == 0 ||
            namePeriod.CompareTo("LHC17r") == 0 )       return k17pp13TeV;
  else if ( namePeriod.CompareTo("LHC17g") == 0 )       return k17pp13TeVLow;
  else if ( namePeriod.CompareTo("LHC17pq") == 0 ||
            namePeriod.CompareTo("LHC17p") == 0 ||
            namePeriod.CompareTo("LHC17q") == 0  )      return k17pp5TeV;
  else if ( namePeriod.CompareTo("LHC17n") == 0 )       return k17XeXe5440GeV;
  else return kNoMC;
}
