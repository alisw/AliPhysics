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
    // Standard NonLinearity - standard kPi0MCv5 for MC and kSDMv5 for data from Jason
    case 1:
      if( clusterType == 0 || clusterType == 1|| clusterType == 3){
        energy *= FunctionNL_kPi0MCv5(energy);
        if(isMC == 0) energy *= FunctionNL_kSDMv5(energy);
      }
      else if ( clusterType == 2 ){
        if(isMC>0)
          energy = FunctionNL_PHOS(energy, 1.008, 0.015, 0.4);
      }
      break;

    // kPi0MCv3 for MC and kTestBeamv3 for data
    case 2:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
      else energy *= FunctionNL_kPi0MCv3(energy);
      break;
    // kPi0MCv3 for MC and kTestBeamv2 for data
    case 3:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
      else energy *= FunctionNL_kPi0MCv3(energy);
      break;

    // kPi0MCv2 for MC and kTestBeamv3 for data
    case 4:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
      else energy *= FunctionNL_kPi0MCv2(energy);
      break;
    // kPi0MCv2 for MC and kTestBeamv2 for data
    case 5:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
      else energy *= FunctionNL_kPi0MCv2(energy);
      break;

    // kPi0MCv1 for MC and kTestBeamv3 for data
    case 6:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
      else energy *= FunctionNL_kPi0MCv1(energy);
      break;
    // kPi0MCv1 for MC and kTestBeamv2 for data
    case 7:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
      else energy *= FunctionNL_kPi0MCv1(energy);
      break;

    // kPi0MCv6 for MC and kSDMv6 for data
    case 8:
      if(isMC == 0) energy *= FunctionNL_kSDMv6(energy);
      else energy *= FunctionNL_kPi0MCv6(energy);
      break;

//----------------------------------------------------------------------------------------------------------

// *************** 10 + x **** default tender settings - pp

    // NonLinearity pp ConvCalo - only shifting MC - no timing cut
    case 11:
      label_case_11:
      if(isMC>0){
        // 8TeV LHC12x
        //pass1
        if( fCurrentMC==k14e2a || fCurrentMC==k14e2b )
          energy /= FunctionNL_kSDM(energy, 0.983251, -3.44339, -1.70998);

        else if( fCurrentMC==k14e2c )
          energy /= FunctionNL_kSDM(energy, 0.984462, -3.00363, -2.63773);

        //pass2
        else if( fCurrentMC == k15h1 )
          energy /= FunctionNL_kSDM(energy, 0.96874*0.991*0.9958*0.999, -3.76064, -0.193181);

        else if( fCurrentMC == k15h2 )
          energy /= FunctionNL_kSDM(energy, 0.969703*0.989*0.9969*0.9991, -3.80387, -0.200546);

        else if( fCurrentMC == k16c2 || fCurrentMC == k16c2_plus )
          energy /= FunctionNL_kSDM(energy, 0.974859*0.987*0.996, -3.85842, -0.405277);

        // 2.76TeV LHC11a/LHC13g
        else if( fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 )
          energy /= FunctionNL_kSDM(energy, 0.984889*0.995*0.9970, -3.65456, -1.12744);

        else if(fCurrentMC==k12f1b)
          energy /= FunctionNL_kSDM(energy, 0.984384*0.995*0.9970, -3.30287, -1.48516);

        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b || fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(energy, 0.981892*0.995*0.9970, -5.43438, -1.05468);

        // 7 TeV LHC10x
//        else if( fCurrentMC==k14j4 ) //v1
//          energy /= FunctionNL_kSDM(energy, 0.973866*0.99*0.996*0.999, -4.06436, -0.379);
        else if( fCurrentMC==k14j4 ) //v2
          energy /= FunctionNL_kSDM(energy, 0.974525*0.986*0.999, -4.00247, -0.453046);

        // pp 5.02 TeV LHC15n
        else if( fCurrentMC==k16k5a ) {
          if(clusterType==3) energy /= FunctionNL_kSDM(energy, 0.980211, -4.374598, -0.171988);
          if(clusterType==1) energy /= FunctionNL_kSDM(energy, 0.984876, -9.999609, -4.999891);
        }
        else if( fCurrentMC==k16k5b ) {
          if(clusterType==3) energy /= FunctionNL_kSDM(energy, 0.981417, -2.772002, -0.955616);
          if(clusterType==1) energy /= FunctionNL_kSDM(energy, 0.980275, -3.172374, -0.730326);
        }

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 12:
      label_case_12:
      if(isMC>0){
        // 8TeV LHC12x
        //pass1
        if( fCurrentMC==k14e2a || fCurrentMC==k14e2b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.967301, -3.1683, -0.653058);

        else if( fCurrentMC==k14e2c )
          energy /= FunctionNL_kSDM(2.0*energy, 0.96728, -2.96279, -0.903677);

        //pass2
        else if( fCurrentMC == k15h1 )
          energy /= FunctionNL_kSDM(energy, 0.963379*0.9985*0.9992, -3.61217, -0.614043);

        else if( fCurrentMC == k15h2 )
          energy /= FunctionNL_kSDM(energy, 0.96105*0.999*0.9996, -3.62239, -0.556256);

        else if( fCurrentMC == k16c2 || fCurrentMC == k16c2_plus )
          energy /= FunctionNL_kSDM(energy, 0.960596*0.999*0.999, -3.48444, -0.766862);

        // 2.76TeV LHC11a/LHC13g
        else if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 )
          energy /= FunctionNL_kSDM(2.0*energy, 0.966151*0.995*0.9981, -2.97974, -0.29463);

        else if( fCurrentMC==k12f1b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.988814*0.995*0.9981, 0.335011, -4.30322);

        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b || fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.979994*0.995*0.9981, -3.24431, -0.760205);

        // 7TeV LHC10x
//        else if(  fCurrentMC==k14j4 ) //v1
//          energy /= FunctionNL_kSDM(energy, 0.955095*0.9991, -3.44162, -0.486573);
        else if(  fCurrentMC==k14j4 ) //v2
          energy /= FunctionNL_kSDM(energy, 0.962095*0.9991*0.9993, -3.63967, -0.747825);

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    case 13:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_11;// goto previous case for shifting MC
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 14:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_12;// goto previous case for shifting MC
      break;

    // NonLinearity ConvCalo - kPi0MC + kSDM
    case 15:
      // 8TeV LHC12x
      if ( fCurrentMC==k14e2a || fCurrentMC==k14e2b  || fCurrentMC==k14e2c || fCurrentMC == k15h1 || fCurrentMC == k15h2  || fCurrentMC == k12pp8TeV || fCurrentMC == k16c2 || fCurrentMC == k16c2_plus ){
        energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04979, 1.3, 0.0967998, 219.381, 63.1604, 1.011);
        if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9846, -3.319, -2.033);

      // 2.76TeV LHC11a/LHC13g
      } else if ( fCurrentMC == k12f1a || fCurrentMC == k12i3 || fCurrentMC == k15g2 || fCurrentMC == k12f1b ||
                  fCurrentMC == k15g1a || fCurrentMC == k15g1b || fCurrentMC == k15a3a || fCurrentMC == k15a3a_plus || fCurrentMC == k15a3b ||
                  fCurrentMC == k11pp2760GeV || fCurrentMC == k13pp2760GeV
                ) {
        energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04123, 1.045, 0.0967998, 219.381, 63.1604, 1.014);
        if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9807*0.995*0.9970, -3.377, -0.8535);
      }
      else fPeriodNameAvailable = kFALSE;
      break;

    // NonLinearity Calo - kPi0MC + kSDM
    case 16:
      // 8TeV LHC12x
      if ( fCurrentMC==k14e2a || fCurrentMC==k14e2b  || fCurrentMC==k14e2c || fCurrentMC == k15h1 || fCurrentMC == k15h2  || fCurrentMC == k12pp8TeV || fCurrentMC == k16c2 || fCurrentMC == k16c2_plus ){
        energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06539, 1.121, 0.0967998, 219.381, 63.1604, 1.011);
        if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9676, -3.216, -0.6828);

      // 2.76TeV LHC11a/LHC13g
      } else if ( fCurrentMC == k12f1a || fCurrentMC == k12i3 || fCurrentMC == k15g2 || fCurrentMC == k12f1b ||
                  fCurrentMC == k15g1a || fCurrentMC == k15g1b || fCurrentMC == k15a3a || fCurrentMC == k15a3a_plus || fCurrentMC == k15a3b ||
                  fCurrentMC == k11pp2760GeV || fCurrentMC == k13pp2760GeV
                ) {
        energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06115, 0.9535, 0.0967998, 219.381, 63.1604, 1.013);
        if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9772*0.995*0.9981, -3.256, -0.4449);
      }
      else fPeriodNameAvailable = kFALSE;
      break;

// *************** 20 + x **** modified tender Settings 1 - pp
    // NonLinearity pp ConvCalo - only shifting MC - no timing cut
    case 21:
      label_case_21:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 )
          energy /= (FunctionNL_DPOW(energy, 1.0443938253, -0.0691830812, -0.1247555443, 1.1673716264, -0.1853095466, -0.0848801702) - 0.0055);
        else if(fCurrentMC==k15g2)
          energy /= (FunctionNL_DPOW(energy, 1.1716155406, -0.1962930603, -0.0193959829, 1.0336659741, -0.0467778485, -0.4407662248) - 0.0055);
        else if(fCurrentMC==k12f1b)
          energy /= (FunctionNL_DPOW(energy, 1.0166321784, -0.0440799552, -0.2611899222, 1.0636538464, -0.0816662488, -0.2173961316) - 0.007);
        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b )
          energy /= (FunctionNL_DPOW(energy, 1.1100193881, -0.1389194936, -0.0800000242, 1.1673716264, -0.1853095466, -0.0848801702) - 0.017);
        else if( fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= (FunctionNL_DPOW(energy, 1.0520183153, -0.0806102847, -0.1450415920, 1.0336724056, -0.0467844121, -0.4406992764) - 0.016);
        // 8TeV
        else if( fCurrentMC == k15h1 )
          energy /= (FunctionNL_DPOW(energy, 1.0654169768, -0.0935785719, -0.1137883054, 1.1814766150, -0.1980098061, -0.0854569214) - 0.0138);
        else if( fCurrentMC == k15h2 )
          energy /= (FunctionNL_DPOW(energy, 1.0652493513, -0.0929276101, -0.1113762695, 1.1837801885, -0.1999914832, -0.0854569214) - 0.0145);
        else if( fCurrentMC == k16c2 || fCurrentMC == k16c2_plus )
          energy /= (FunctionNL_DPOW(energy, 1.0489259285, -0.0759079646, -0.1239772934, 1.1835846739, -0.1998987993, -0.0854186691) - 0.014);
        // 7 TeV
//        else if( fCurrentMC == k14j4 ) //v1
//          energy /= (FunctionNL_DPOW(energy, 1.1086453117, -0.1373335557, -0.0800000000, 1.1855482066, -0.1999999504, -0.0830177063) - 0.014);
        else if( fCurrentMC == k14j4 ) //v2
          energy /= (FunctionNL_DPOW(energy, 1.1082846035, -0.1369968318, -0.0800000002, 1.1850179319, -0.1999999950, -0.0863054172) - 0.015);
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 22:
      label_case_22:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 )
          energy /= (FunctionNL_DPOW(energy, 0.9980625418, -0.0564782662, -0.5, 1.0383412435, -0.0851830429, -0.4999999996) - 0.00175);
        else if( fCurrentMC==k15g2 )
          energy /= (FunctionNL_DPOW(energy, 1.0795372569, -0.1347324732, -0.1630736190, 1.1614181498, -0.199995361, -0.1711378093) - 0.0035);
        else if( fCurrentMC==k12f1b )
          energy /= (FunctionNL_DPOW(energy, 1.0232969083, -0.090409434, -0.3592406513, 1.0383412435, -0.0851830429, -0.4999999996) + 0.0007);
        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b )
          energy /= (FunctionNL_DPOW(energy, 1.0106037132, -0.0748250591, -0.4999999996, 1.0383412435, -0.0851830429, -0.4999999996) - 0.014);
        else if( fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= (FunctionNL_DPOW(energy, 1.0119417393, -0.0755250741, -0.4999999996, 1.1614181498, -0.1999995361, -0.1711378093) - 0.006);
        //8TeV
        else if( fCurrentMC == k15h1 )
          energy /= (FunctionNL_DPOW(energy, 1.1389201636, -0.1999994717, -0.1622237979, 1.1603460704, -0.1999999989, -0.2194447313) - 0.0025);
        else if( fCurrentMC == k15h2 )
          energy /= (FunctionNL_DPOW(energy, 1.0105301622, -0.0732424689, -0.5000000000, 1.0689250170, -0.1082682369, -0.4388156470) - 0.001);
        else if( fCurrentMC == k16c2 || fCurrentMC == k16c2_plus )
          energy /= (FunctionNL_DPOW(energy, 0.9922456908, -0.0551212559, -0.5000000000, 1.0513459039, -0.0894163252, -0.5000000000) + 0.002);
        // 7 TeV
//        else if( fCurrentMC == k14j4 ) //v1
//          energy /= (FunctionNL_DPOW(energy, 0.9994789138, -0.0601419399, -0.4999999999, 1.1635744933, -0.1999999978, -0.1985578372) - 0.005);
        else if( fCurrentMC == k14j4 ) //v2
          energy /= (FunctionNL_DPOW(energy, 1.0074002842, -0.0682543971, -0.4509341085, 1.1224162203, -0.1586806096, -0.2458351112) - 0.003);
        else fPeriodNameAvailable = kFALSE;
      }
      break;
    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    case 23:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_21;// goto previous case for shifting MC
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 24:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_22;// goto previous case for shifting MC
      break;


// *************** 30 + x **** modified tender Settings 2 - pp
    case 31:
      label_case_31:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 )
          energy /= FunctionNL_kSDM(energy, 0.983176*0.9945, -3.91107, -0.697613);
        else if(fCurrentMC==k15g2)
          energy /= FunctionNL_kSDM(energy, 0.972574*0.9942, -3.19191, -0.946239);
        else if(fCurrentMC==k12f1b)
          energy /= FunctionNL_kSDM(energy, 0.981893*0.9930, -4.05476, -0.710661);
        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b )
          energy /= FunctionNL_kSDM(energy, 0.983176*0.993*0.99, -1.85546, -3.37696);
        else if( fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(energy, 0.977035*0.9835, -3.82187, -1.04332);
        else fPeriodNameAvailable = kFALSE;
      }
      break;
    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 32:
      label_case_32:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 )
          energy /= FunctionNL_kSDM(energy, 0.974358*0.9987, -2.18037, -1.91622);
        else if( fCurrentMC==k15g2 )
          energy /= FunctionNL_kSDM(energy, 0.963307*0.9962, -3.27998, -0.589806);
        else if( fCurrentMC==k12f1b )
          energy /= FunctionNL_kSDM(energy, 0.97499*0.9995, -0.180148, -4.78066);
        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b )
          energy /= FunctionNL_kSDM(energy, 0.974424*0.998*0.992, -0.533785, -4.06374);
        else if ( fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(energy, 0.963307*0.995, -4.01949, -0.38667);
        else fPeriodNameAvailable = kFALSE;
      }
      break;
    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    case 33:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_31;// goto previous case for shifting MC
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 34:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_32;// goto previous case for shifting MC
      break;


// *************** 40 + x **** default tender Settings - pPb
    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 41:
      label_case_41:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ){
          energy /= FunctionNL_kSDM(energy, 0.995*0.978578, -3.80517, -0.581197);//v4
          energy /= FunctionNL_kSDM(energy, 0.996179, -5.33609, -0.477463);//v5
        }
        else if( fCurrentMC==k13e7 ) energy /= FunctionNL_kSDM(energy, 0.979813, -3.53445, -0.733067);//v0

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 42:
      label_case_42:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c )
          energy /= FunctionNL_kSDM(energy, 1.002*0.970383, -3.65936, -0.721139);//v4
        else if( fCurrentMC==k13e7 )
          energy /= FunctionNL_kSDM(energy, 0.970537, -3.36675, -0.958747);//v0
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb ConvCalo  - kTestBeamv3 + shifting MC
    case 43:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_41;// goto previous case for shifting MC
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 44:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_42;// goto previous case for shifting MC
      break;

    // NonLinearity LHC13 pPb Calo - excluding the two lowest pT points
    case 49:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC==k13e7 || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c )
          energy /= FunctionNL_kSDM(energy, 0.973302, -3.12524, -1.13546);
        else fPeriodNameAvailable = kFALSE;
      }
      break;

// *************** 50 + x **** modified tender Settings 1 - pPb
    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 51:
      label_case_51:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c )
          energy /= (FunctionNL_DPOW(energy, 1.0754004911, -0.0992327361, -0.0802161499, 1.1849304274, -0.1999999986, -0.0828138864) - 0.005);//v4
        else if( fCurrentMC==k13e7 )
          energy /= FunctionNL_DPOW(energy, 1.0546114304, -0.0758513555, -0.0800000002, 1.1849400584, -0.1999999970, -0.0826417756);//v4
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 52:
      label_case_52:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c )
          energy /= (FunctionNL_DPOW(energy, 1.0358207569, -0.0914347267, -0.3683743201, 1.1549558754, -0.1942615277, -0.2216281109) + 0.002);//v4
        else if( fCurrentMC==k13e7 )
          energy /= FunctionNL_DPOW(energy, 1.0149551972, -0.0697288693, -0.4586527438, 1.1549558754, -0.1942615277, -0.2216281109);//v4
        else fPeriodNameAvailable = kFALSE;
      }
      break;


// *************** 60 + x **** modified tender Settings 2 - pPb


// *************** 70 + x **** default tender Settings - PbPb


// *************** 80 + x **** modified tender Settings 1 - PbPb


// *************** 90 + x **** modified tender Settings 2 - PbPb



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
Float_t AliCaloNonLinearity::FunctionNL_PHOS(Float_t e, Float_t p0, Float_t p1, Float_t p2){
  return (0.0241+1.0504*e+0.000249*e*e)*p0*(1+p1/(1.+e*e/p2/p2)) ;
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
  if(       namePeriod.CompareTo("LHC14e2a")==0)        return k14e2a;
  else if(  namePeriod.CompareTo("LHC14e2b")==0)        return k14e2b;
  else if(  namePeriod.CompareTo("LHC14e2c")==0)        return k14e2c;
  else if(  namePeriod.CompareTo("LHC12f1a")==0)        return k12f1a;
  else if(  namePeriod.CompareTo("LHC12f1b")==0)        return k12f1b;
  else if(  namePeriod.CompareTo("LHC12i3")==0)         return k12i3;
  else if(  namePeriod.CompareTo("LHC15g1a")==0)        return k15g1a;
  else if(  namePeriod.CompareTo("LHC15g1b")==0)        return k15g1b;
  else if(  namePeriod.CompareTo("LHC15g2")==0)         return k15g2;
  else if(  namePeriod.CompareTo("LHC15a3a")==0)        return k15a3a;
  else if(  namePeriod.CompareTo("LHC15a3a_plus")==0)   return k15a3a_plus;
  else if(  namePeriod.CompareTo("LHC15a3b")==0)        return k15a3b;
  else if(  namePeriod.Contains("LHC13b2_efix"))        return k13b2_efix;
  else if(  namePeriod.Contains("LHC13e7"))             return k13e7;
  else if(  namePeriod.CompareTo("LHC15h1a1")==0 ||
            namePeriod.CompareTo("LHC15h1b")==0 ||
            namePeriod.CompareTo("LHC15h1c")==0 ||
            namePeriod.CompareTo("LHC15h1d")==0 ||
            namePeriod.CompareTo("LHC15h1f")==0 ||
            namePeriod.CompareTo("LHC15h1g")==0 ||
            namePeriod.CompareTo("LHC15h1h")==0 ||
            namePeriod.CompareTo("LHC15h1i")==0)        return k15h1;
  else if(  namePeriod.CompareTo("LHC15h2a")==0 ||
            namePeriod.CompareTo("LHC15h2b")==0 ||
            namePeriod.CompareTo("LHC15h2c")==0 ||
            namePeriod.CompareTo("LHC15h2d")==0 ||
            namePeriod.CompareTo("LHC15h2f")==0 ||
            namePeriod.CompareTo("LHC15h2g")==0 ||
            namePeriod.CompareTo("LHC15h2h")==0 ||
            namePeriod.CompareTo("LHC15h2i")==0)        return k15h2;
  else if(  namePeriod.CompareTo("LHC14j4b")==0 ||
            namePeriod.CompareTo("LHC14j4c")==0 ||
            namePeriod.CompareTo("LHC14j4d")==0 ||
            namePeriod.CompareTo("LHC14j4e")==0 ||
            namePeriod.CompareTo("LHC14j4f")==0)        return k14j4;
  else if ( namePeriod.CompareTo("LHC16c2") == 0 )      return k16c2;
  else if ( namePeriod.CompareTo("LHC16c2_plus") == 0 ) return k16c2_plus;
  else if ( namePeriod.CompareTo("LHC16c3a") == 0 )     return k16c3a;
  else if ( namePeriod.CompareTo("LHC16c3b") == 0 )     return k16c3b;
  else if ( namePeriod.CompareTo("LHC16c3c") == 0 )     return k16c3c;
  else if ( namePeriod.CompareTo("LHC16h3") == 0 )      return k16h3;
  else if ( namePeriod.CompareTo("LHC16h3b") == 0 )     return k16h3b;
  else if ( namePeriod.CompareTo("LHC16h8a") == 0 )     return k16h8a;
  else if ( namePeriod.CompareTo("LHC16h8b") == 0 )     return k16h8b;
  else if ( namePeriod.CompareTo("LHC16k3a") == 0 ||
            namePeriod.CompareTo("LHC16k3a2") == 0 )     return k16k3a;
  else if ( namePeriod.CompareTo("LHC16k3b") == 0 ||
            namePeriod.CompareTo("LHC16k3b2") == 0 )     return k16k3b;
  else if ( namePeriod.CompareTo("LHC16k5a") == 0  )     return k16k5a;
  else if ( namePeriod.CompareTo("LHC16k5b") == 0  )     return k16k5b;
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
  else if ( namePeriod.CompareTo("LHC15f") == 0 ||
            namePeriod.CompareTo("LHC15g") == 0 ||
            namePeriod.CompareTo("LHC15h") == 0 ||
            namePeriod.CompareTo("LHC15i") == 0 ||
            namePeriod.CompareTo("LHC15j") == 0 ||
            namePeriod.CompareTo("LHC15k") == 0 ||
            namePeriod.CompareTo("LHC15l") == 0 ||
            namePeriod.CompareTo("LHC15m") == 0 )       return k15pp13TeV;
  else if ( namePeriod.CompareTo("LHC15n") == 0 )       return k15pp5TeV;
  else if ( namePeriod.CompareTo("LHC15o") == 0 )       return k15PbPb5TeV;
  else return kNoMC;
}