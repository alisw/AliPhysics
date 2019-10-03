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


#include <TArrayD.h>
#include <TArrayL.h>
#include <TObject.h>
#include <TMath.h>

#include "AliAnalysisEmcalJetHelperEA.h"

#include <stdio.h>
#include <stdlib.h>

using namespace PWGJE::EMCALJetTasks;
using namespace std;

// Run by run mean V0 multiplicities detector level and  particle level in pp 13 TeV 
// Author Filip Krizek   (1.10. 2019)

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisEmcalJetHelperEA);
/// \endcond

//________________________________________________________________________________________

AliAnalysisEmcalJetHelperEA::AliAnalysisEmcalJetHelperEA():
   TObject(),
   fnRun(0),
   fMeanV0A_PartLevel(13.6185),
   fMeanV0C_PartLevel(13.4742),
   fMeanV0M_PartLevel(27.0927),
   fMeanV0A_DetLevel(41.7158),
   fMeanV0C_DetLevel(58.3741),
   fMeanV0M_DetLevel(100.09)
{
   //default constructor

   fMeanV0A.Set(2000);
   fMeanV0C.Set(2000);
   fMeanV0M.Set(2000);
   fRuns.Set(2000);
 

   for(Int_t i=0; i<2000; i++){
      fMeanV0A[i] = 1.;
      fMeanV0C[i] = 1.;
      fMeanV0M[i] = 1.;
      fRuns[i]    = 0;
   }



   Int_t io = 0; 
   // 	FILTER_p-p_208_LHC16d  
   fRuns[io] = 252235;  fMeanV0A[io] = 50.5654;  fMeanV0C[io] = 77.4974;  fMeanV0M[io] = 127.957;  io++;
   fRuns[io] = 252248;  fMeanV0A[io] = 50.2638;  fMeanV0C[io] = 77.2641;  fMeanV0M[io] = 127.430;  io++;
   fRuns[io] = 252271;  fMeanV0A[io] = 50.6799;  fMeanV0C[io] = 77.8131;  fMeanV0M[io] = 128.399;  io++;
   fRuns[io] = 252310;  fMeanV0A[io] = 50.4039;  fMeanV0C[io] = 77.0152;  fMeanV0M[io] = 127.311;  io++;
   fRuns[io] = 252317;  fMeanV0A[io] = 50.4902;  fMeanV0C[io] = 77.4425;  fMeanV0M[io] = 127.834;  io++;
   fRuns[io] = 252319;  fMeanV0A[io] = 50.4825;  fMeanV0C[io] = 77.6853;  fMeanV0M[io] = 128.069;  io++;
   fRuns[io] = 252322;  fMeanV0A[io] = 50.5860;  fMeanV0C[io] = 77.6280;  fMeanV0M[io] = 128.112;  io++;
   fRuns[io] = 252325;  fMeanV0A[io] = 50.5471;  fMeanV0C[io] = 77.7205;  fMeanV0M[io] = 128.171;  io++;
   fRuns[io] = 252326;  fMeanV0A[io] = 50.6368;  fMeanV0C[io] = 77.7966;  fMeanV0M[io] = 128.327;  io++;
   fRuns[io] = 252330;  fMeanV0A[io] = 50.6204;  fMeanV0C[io] = 77.4760;  fMeanV0M[io] = 127.998;  io++;

   //FILTER_p-p_208_LHC16e
   fRuns[io] = 253437;  fMeanV0A[io] = 50.0851;  fMeanV0C[io] = 75.9878;  fMeanV0M[io] = 125.979;  io++;
   fRuns[io] = 253478;  fMeanV0A[io] = 49.9131;  fMeanV0C[io] = 75.9397;  fMeanV0M[io] = 125.745;  io++;
   fRuns[io] = 253481;  fMeanV0A[io] = 49.8021;  fMeanV0C[io] = 75.9260;  fMeanV0M[io] = 125.627;  io++;
   fRuns[io] = 253482;  fMeanV0A[io] = 49.7041;  fMeanV0C[io] = 75.7614;  fMeanV0M[io] = 125.364;  io++;
   fRuns[io] = 253488;  fMeanV0A[io] = 49.6817;  fMeanV0C[io] = 75.7138;  fMeanV0M[io] = 125.290;  io++;
   fRuns[io] = 253517;  fMeanV0A[io] = 50.1114;  fMeanV0C[io] = 75.9579;  fMeanV0M[io] = 125.967;  io++;
   fRuns[io] = 253529;  fMeanV0A[io] = 49.8242;  fMeanV0C[io] = 75.6876;  fMeanV0M[io] = 125.407;  io++;
   fRuns[io] = 253530;  fMeanV0A[io] = 49.7841;  fMeanV0C[io] = 75.6081;  fMeanV0M[io] = 125.285;  io++;
   fRuns[io] = 253563;  fMeanV0A[io] = 50.3303;  fMeanV0C[io] = 75.8570;  fMeanV0M[io] = 126.085;  io++;
   fRuns[io] = 253589;  fMeanV0A[io] = 50.3221;  fMeanV0C[io] = 75.9468;  fMeanV0M[io] = 126.164;  io++;
   fRuns[io] = 253591;  fMeanV0A[io] = 50.1024;  fMeanV0C[io] = 75.5257;  fMeanV0M[io] = 125.533;  io++;

   //   FILTER_p-p_208_LHC16g 
   fRuns[io] = 254128;  fMeanV0A[io] = 45.1506;  fMeanV0C[io] = 70.3466;  fMeanV0M[io] = 115.379;  io++;
   fRuns[io] = 254147;  fMeanV0A[io] = 45.0000;  fMeanV0C[io] = 70.0847;  fMeanV0M[io] = 114.968;  io++;
   fRuns[io] = 254149;  fMeanV0A[io] = 44.8180;  fMeanV0C[io] = 69.9104;  fMeanV0M[io] = 114.602;  io++;
   fRuns[io] = 254174;  fMeanV0A[io] = 44.8628;  fMeanV0C[io] = 69.4357;  fMeanV0M[io] = 114.181;  io++;
   fRuns[io] = 254175;  fMeanV0A[io] = 44.6071;  fMeanV0C[io] = 69.0512;  fMeanV0M[io] = 113.532;  io++;
   fRuns[io] = 254178;  fMeanV0A[io] = 44.4549;  fMeanV0C[io] = 68.9076;  fMeanV0M[io] = 113.240;  io++;
   fRuns[io] = 254193;  fMeanV0A[io] = 44.2683;  fMeanV0C[io] = 68.5877;  fMeanV0M[io] = 112.734;  io++;
   fRuns[io] = 254199;  fMeanV0A[io] = 44.1858;  fMeanV0C[io] = 68.3724;  fMeanV0M[io] = 112.436;  io++;
   fRuns[io] = 254204;  fMeanV0A[io] = 44.0420;  fMeanV0C[io] = 68.3472;  fMeanV0M[io] = 112.259;  io++;
   fRuns[io] = 254205;  fMeanV0A[io] = 43.9808;  fMeanV0C[io] = 68.2510;  fMeanV0M[io] = 112.107;  io++;
   fRuns[io] = 254293;  fMeanV0A[io] = 44.2721;  fMeanV0C[io] = 68.7171;  fMeanV0M[io] = 112.880;  io++;
   fRuns[io] = 254302;  fMeanV0A[io] = 44.1524;  fMeanV0C[io] = 68.4022;  fMeanV0M[io] = 112.434;  io++;
   fRuns[io] = 254303;  fMeanV0A[io] = 44.3168;  fMeanV0C[io] = 68.4798;  fMeanV0M[io] = 112.660;  io++;
   fRuns[io] = 254304;  fMeanV0A[io] = 44.1027;  fMeanV0C[io] = 68.2963;  fMeanV0M[io] = 112.283;  io++;
   fRuns[io] = 254330;  fMeanV0A[io] = 44.2390;  fMeanV0C[io] = 68.3718;  fMeanV0M[io] = 112.496;  io++;
   fRuns[io] = 254331;  fMeanV0A[io] = 44.0142;  fMeanV0C[io] = 67.9854;  fMeanV0M[io] = 111.874;  io++;
   fRuns[io] = 254332;  fMeanV0A[io] = 43.7450;  fMeanV0C[io] = 67.4760;  fMeanV0M[io] = 111.100;  io++;

   //FILTER_p-p_208_LHC16h
   fRuns[io] = 254604;  fMeanV0A[io] = 44.3933;  fMeanV0C[io] = 67.9688;  fMeanV0M[io] = 112.232;  io++;
   fRuns[io] = 254606;  fMeanV0A[io] = 44.1275;  fMeanV0C[io] = 67.6188;  fMeanV0M[io] = 111.618;  io++;
   fRuns[io] = 254621;  fMeanV0A[io] = 43.0362;  fMeanV0C[io] = 66.2734;  fMeanV0M[io] = 109.184;  io++;
   fRuns[io] = 254629;  fMeanV0A[io] = 42.8620;  fMeanV0C[io] = 65.9648;  fMeanV0M[io] = 108.692;  io++;
   fRuns[io] = 254630;  fMeanV0A[io] = 42.8106;  fMeanV0C[io] = 65.9450;  fMeanV0M[io] = 108.621;  io++;
   fRuns[io] = 254632;  fMeanV0A[io] = 42.7535;  fMeanV0C[io] = 65.8630;  fMeanV0M[io] = 108.491;  io++;
   fRuns[io] = 254640;  fMeanV0A[io] = 42.6849;  fMeanV0C[io] = 65.5940;  fMeanV0M[io] = 108.145;  io++;
   fRuns[io] = 254644;  fMeanV0A[io] = 42.5947;  fMeanV0C[io] = 65.4563;  fMeanV0M[io] = 107.920;  io++;
   fRuns[io] = 254646;  fMeanV0A[io] = 42.6034;  fMeanV0C[io] = 65.3453;  fMeanV0M[io] = 107.820;  io++;
   fRuns[io] = 254648;  fMeanV0A[io] = 42.4720;  fMeanV0C[io] = 65.3031;  fMeanV0M[io] = 107.636;  io++;
   fRuns[io] = 254649;  fMeanV0A[io] = 42.3782;  fMeanV0C[io] = 65.1520;  fMeanV0M[io] = 107.396;  io++;
   fRuns[io] = 254651;  fMeanV0A[io] = 42.2944;  fMeanV0C[io] = 65.0897;  fMeanV0M[io] = 107.257;  io++;
   fRuns[io] = 254652;  fMeanV0A[io] = 42.1544;  fMeanV0C[io] = 64.8009;  fMeanV0M[io] = 106.826;  io++;
   fRuns[io] = 254653;  fMeanV0A[io] = 42.1629;  fMeanV0C[io] = 64.8054;  fMeanV0M[io] = 106.839;  io++;
   fRuns[io] = 254654;  fMeanV0A[io] = 41.9686;  fMeanV0C[io] = 64.5126;  fMeanV0M[io] = 106.367;  io++;
   fRuns[io] = 254983;  fMeanV0A[io] = 43.8181;  fMeanV0C[io] = 67.6738;  fMeanV0M[io] = 111.373;  io++;
   fRuns[io] = 254984;  fMeanV0A[io] = 43.7140;  fMeanV0C[io] = 67.3728;  fMeanV0M[io] = 110.968;  io++;
   fRuns[io] = 255079;  fMeanV0A[io] = 42.1432;  fMeanV0C[io] = 64.5147;  fMeanV0M[io] = 106.517;  io++;
   fRuns[io] = 255082;  fMeanV0A[io] = 41.8928;  fMeanV0C[io] = 64.1709;  fMeanV0M[io] = 105.944;  io++;
   fRuns[io] = 255085;  fMeanV0A[io] = 41.9863;  fMeanV0C[io] = 64.2953;  fMeanV0M[io] = 106.139;  io++;
   fRuns[io] = 255086;  fMeanV0A[io] = 41.8729;  fMeanV0C[io] = 64.2697;  fMeanV0M[io] = 106.011;  io++;
   fRuns[io] = 255091;  fMeanV0A[io] = 41.8177;  fMeanV0C[io] = 64.1242;  fMeanV0M[io] = 105.808;  io++;
   fRuns[io] = 255111;  fMeanV0A[io] = 41.9749;  fMeanV0C[io] = 64.3180;  fMeanV0M[io] = 106.151;  io++;
   fRuns[io] = 255154;  fMeanV0A[io] = 42.2102;  fMeanV0C[io] = 64.6172;  fMeanV0M[io] = 106.696;  io++;
   fRuns[io] = 255159;  fMeanV0A[io] = 41.8646;  fMeanV0C[io] = 64.1199;  fMeanV0M[io] = 105.853;  io++;
   fRuns[io] = 255162;  fMeanV0A[io] = 41.8096;  fMeanV0C[io] = 64.2714;  fMeanV0M[io] = 105.947;  io++;
   fRuns[io] = 255167;  fMeanV0A[io] = 41.6892;  fMeanV0C[io] = 64.0895;  fMeanV0M[io] = 105.655;  io++;
   fRuns[io] = 255171;  fMeanV0A[io] = 41.6272;  fMeanV0C[io] = 63.7777;  fMeanV0M[io] = 105.271;  io++;
   fRuns[io] = 255173;  fMeanV0A[io] = 41.5744;  fMeanV0C[io] = 63.8693;  fMeanV0M[io] = 105.299;  io++;
   fRuns[io] = 255174;  fMeanV0A[io] = 41.5868;  fMeanV0C[io] = 64.0504;  fMeanV0M[io] = 105.523;  io++;
   fRuns[io] = 255176;  fMeanV0A[io] = 41.5327;  fMeanV0C[io] = 63.7467;  fMeanV0M[io] = 105.145;  io++;
   fRuns[io] = 255177;  fMeanV0A[io] = 41.3161;  fMeanV0C[io] = 63.6038;  fMeanV0M[io] = 104.784;  io++;
   fRuns[io] = 255180;  fMeanV0A[io] = 41.2691;  fMeanV0C[io] = 63.5439;  fMeanV0M[io] = 104.679;  io++;
   fRuns[io] = 255181;  fMeanV0A[io] = 41.4022;  fMeanV0C[io] = 63.6936;  fMeanV0M[io] = 104.990;  io++;
   fRuns[io] = 255182;  fMeanV0A[io] = 41.2617;  fMeanV0C[io] = 63.5596;  fMeanV0M[io] = 104.681;  io++;
   fRuns[io] = 255240;  fMeanV0A[io] = 41.3592;  fMeanV0C[io] = 63.6545;  fMeanV0M[io] = 104.874;  io++;
   fRuns[io] = 255242;  fMeanV0A[io] = 41.2594;  fMeanV0C[io] = 63.4500;  fMeanV0M[io] = 104.579;  io++;
   fRuns[io] = 255247;  fMeanV0A[io] = 40.8939;  fMeanV0C[io] = 63.0076;  fMeanV0M[io] = 103.778;  io++;
   fRuns[io] = 255248;  fMeanV0A[io] = 40.9299;  fMeanV0C[io] = 63.3164;  fMeanV0M[io] = 104.111;  io++;
   fRuns[io] = 255249;  fMeanV0A[io] = 40.8012;  fMeanV0C[io] = 63.0696;  fMeanV0M[io] = 103.729;  io++;
   fRuns[io] = 255251;  fMeanV0A[io] = 40.6592;  fMeanV0C[io] = 62.7180;  fMeanV0M[io] = 103.244;  io++;
   fRuns[io] = 255252;  fMeanV0A[io] = 40.4342;  fMeanV0C[io] = 62.3911;  fMeanV0M[io] = 102.689;  io++;
   fRuns[io] = 255253;  fMeanV0A[io] = 40.3290;  fMeanV0C[io] = 62.3603;  fMeanV0M[io] = 102.557;  io++;
   fRuns[io] = 255255;  fMeanV0A[io] = 40.2766;  fMeanV0C[io] = 62.2510;  fMeanV0M[io] = 102.392;  io++;
   fRuns[io] = 255256;  fMeanV0A[io] = 40.1143;  fMeanV0C[io] = 62.1446;  fMeanV0M[io] = 102.123;  io++;
   fRuns[io] = 255275;  fMeanV0A[io] = 40.2298;  fMeanV0C[io] = 62.0907;  fMeanV0M[io] = 102.186;  io++;
   fRuns[io] = 255276;  fMeanV0A[io] = 39.8723;  fMeanV0C[io] = 61.6588;  fMeanV0M[io] = 101.395;  io++;
   fRuns[io] = 255280;  fMeanV0A[io] = 39.6103;  fMeanV0C[io] = 61.4261;  fMeanV0M[io] = 100.880;  io++;
   fRuns[io] = 255283;  fMeanV0A[io] = 39.5115;  fMeanV0C[io] = 61.3425;  fMeanV0M[io] = 100.726;  io++;
   fRuns[io] = 255350;  fMeanV0A[io] = 39.9328;  fMeanV0C[io] = 61.8045;  fMeanV0M[io] = 101.613;  io++;
   fRuns[io] = 255351;  fMeanV0A[io] = 39.7017;  fMeanV0C[io] = 61.5193;  fMeanV0M[io] = 101.081;  io++;
   fRuns[io] = 255352;  fMeanV0A[io] = 39.6336;  fMeanV0C[io] = 61.2114;  fMeanV0M[io] = 100.696;  io++;
   fRuns[io] = 255398;  fMeanV0A[io] = 39.6248;  fMeanV0C[io] = 61.6169;  fMeanV0M[io] = 101.103;  io++;
   fRuns[io] = 255402;  fMeanV0A[io] = 39.3788;  fMeanV0C[io] = 61.4122;  fMeanV0M[io] = 100.654;  io++;
   fRuns[io] = 255407;  fMeanV0A[io] = 39.0887;  fMeanV0C[io] = 60.9864;  fMeanV0M[io] = 99.9386;  io++;
   fRuns[io] = 255415;  fMeanV0A[io] = 38.9101;  fMeanV0C[io] = 60.7736;  fMeanV0M[io] = 99.5384;  io++;
   fRuns[io] = 255418;  fMeanV0A[io] = 38.8214;  fMeanV0C[io] = 60.6265;  fMeanV0M[io] = 99.3192;  io++;
   fRuns[io] = 255419;  fMeanV0A[io] = 38.6851;  fMeanV0C[io] = 60.5110;  fMeanV0M[io] = 99.0518;  io++;
   fRuns[io] = 255420;  fMeanV0A[io] = 38.6734;  fMeanV0C[io] = 60.5612;  fMeanV0M[io] = 99.1040;  io++;
   fRuns[io] = 255421;  fMeanV0A[io] = 38.6321;  fMeanV0C[io] = 60.4518;  fMeanV0M[io] = 98.9379;  io++;
   fRuns[io] = 255440;  fMeanV0A[io] = 38.6810;  fMeanV0C[io] = 60.2846;  fMeanV0M[io] = 98.8260;  io++;
   fRuns[io] = 255442;  fMeanV0A[io] = 38.3223;  fMeanV0C[io] = 60.0070;  fMeanV0M[io] = 98.1893;  io++;
   fRuns[io] = 255447;  fMeanV0A[io] = 38.1909;  fMeanV0C[io] = 59.7998;  fMeanV0M[io] = 97.8487;  io++;
   fRuns[io] = 255463;  fMeanV0A[io] = 38.0230;  fMeanV0C[io] = 59.5622;  fMeanV0M[io] = 97.4363;  io++;
   fRuns[io] = 255465;  fMeanV0A[io] = 37.9340;  fMeanV0C[io] = 59.4199;  fMeanV0M[io] = 97.2053;  io++;
   fRuns[io] = 255466;  fMeanV0A[io] = 37.8256;  fMeanV0C[io] = 59.2448;  fMeanV0M[io] = 96.9101;  io++;
   fRuns[io] = 255467;  fMeanV0A[io] = 37.8828;  fMeanV0C[io] = 59.3474;  fMeanV0M[io] = 97.0807;  io++;

   //FILTER_p-p_208_LHC16i
   fRuns[io] = 255539;  fMeanV0A[io] = 37.5186;  fMeanV0C[io] = 58.8381;  fMeanV0M[io] = 96.2088;  io++;
   fRuns[io] = 255540;  fMeanV0A[io] = 37.3519;  fMeanV0C[io] = 58.6419;  fMeanV0M[io] = 95.8565;  io++;
   fRuns[io] = 255541;  fMeanV0A[io] = 37.4124;  fMeanV0C[io] = 58.5785;  fMeanV0M[io] = 95.8520;  io++;
   fRuns[io] = 255542;  fMeanV0A[io] = 37.2186;  fMeanV0C[io] = 58.5283;  fMeanV0M[io] = 95.6009;  io++;
   fRuns[io] = 255543;  fMeanV0A[io] = 37.0002;  fMeanV0C[io] = 58.3689;  fMeanV0M[io] = 95.2210;  io++;
   fRuns[io] = 255577;  fMeanV0A[io] = 37.4128;  fMeanV0C[io] = 58.6536;  fMeanV0M[io] = 95.9215;  io++;
   fRuns[io] = 255582;  fMeanV0A[io] = 37.1460;  fMeanV0C[io] = 58.5270;  fMeanV0M[io] = 95.5242;  io++;
   fRuns[io] = 255583;  fMeanV0A[io] = 37.0478;  fMeanV0C[io] = 58.2216;  fMeanV0M[io] = 95.1242;  io++;
   fRuns[io] = 255591;  fMeanV0A[io] = 36.6901;  fMeanV0C[io] = 57.8495;  fMeanV0M[io] = 94.3830;  io++;
   fRuns[io] = 255614;  fMeanV0A[io] = 36.9364;  fMeanV0C[io] = 57.9284;  fMeanV0M[io] = 94.7144;  io++;
   fRuns[io] = 255615;  fMeanV0A[io] = 36.7757;  fMeanV0C[io] = 57.8926;  fMeanV0M[io] = 94.5190;  io++;
   fRuns[io] = 255616;  fMeanV0A[io] = 36.6153;  fMeanV0C[io] = 57.6735;  fMeanV0M[io] = 94.1351;  io++;
   fRuns[io] = 255617;  fMeanV0A[io] = 36.4980;  fMeanV0C[io] = 57.6068;  fMeanV0M[io] = 93.9495;  io++;
   fRuns[io] = 255618;  fMeanV0A[io] = 36.4774;  fMeanV0C[io] = 57.5393;  fMeanV0M[io] = 93.8620;  io++;

   //FILTER_p-p_208_LHC16j
   fRuns[io] = 256219;  fMeanV0A[io] = 37.1402;  fMeanV0C[io] = 58.7736;  fMeanV0M[io] = 95.7671;  io++;
   fRuns[io] = 256223;  fMeanV0A[io] = 37.0116;  fMeanV0C[io] = 58.4900;  fMeanV0M[io] = 95.3546;  io++;
   fRuns[io] = 256227;  fMeanV0A[io] = 36.7793;  fMeanV0C[io] = 58.1633;  fMeanV0M[io] = 94.7929;  io++;
   fRuns[io] = 256228;  fMeanV0A[io] = 36.7532;  fMeanV0C[io] = 58.0143;  fMeanV0M[io] = 94.6035;  io++;
   fRuns[io] = 256231;  fMeanV0A[io] = 36.7907;  fMeanV0C[io] = 58.0782;  fMeanV0M[io] = 94.7234;  io++;
   fRuns[io] = 256281;  fMeanV0A[io] = 37.0437;  fMeanV0C[io] = 58.2514;  fMeanV0M[io] = 95.1520;  io++;
   fRuns[io] = 256282;  fMeanV0A[io] = 36.9954;  fMeanV0C[io] = 58.2528;  fMeanV0M[io] = 95.0925;  io++;
   fRuns[io] = 256283;  fMeanV0A[io] = 36.7552;  fMeanV0C[io] = 58.0481;  fMeanV0M[io] = 94.6464;  io++;
   fRuns[io] = 256284;  fMeanV0A[io] = 36.5313;  fMeanV0C[io] = 57.8809;  fMeanV0M[io] = 94.2634;  io++;
   fRuns[io] = 256287;  fMeanV0A[io] = 36.4749;  fMeanV0C[io] = 57.8109;  fMeanV0M[io] = 94.1351;  io++;
   fRuns[io] = 256289;  fMeanV0A[io] = 36.4211;  fMeanV0C[io] = 57.6723;  fMeanV0M[io] = 93.9507;  io++;
   fRuns[io] = 256290;  fMeanV0A[io] = 36.3352;  fMeanV0C[io] = 57.5379;  fMeanV0M[io] = 93.7254;  io++;
   fRuns[io] = 256292;  fMeanV0A[io] = 36.4077;  fMeanV0C[io] = 57.7148;  fMeanV0M[io] = 93.9394;  io++;
   fRuns[io] = 256295;  fMeanV0A[io] = 36.3243;  fMeanV0C[io] = 57.5382;  fMeanV0M[io] = 93.7058;  io++;
   fRuns[io] = 256297;  fMeanV0A[io] = 36.2789;  fMeanV0C[io] = 57.3116;  fMeanV0M[io] = 93.4503;  io++;
   fRuns[io] = 256299;  fMeanV0A[io] = 36.2909;  fMeanV0C[io] = 57.4243;  fMeanV0M[io] = 93.5609;  io++;
   fRuns[io] = 256302;  fMeanV0A[io] = 36.1952;  fMeanV0C[io] = 57.3191;  fMeanV0M[io] = 93.3636;  io++;
   fRuns[io] = 256307;  fMeanV0A[io] = 36.0897;  fMeanV0C[io] = 57.1639;  fMeanV0M[io] = 93.0987;  io++;
   fRuns[io] = 256309;  fMeanV0A[io] = 36.1055;  fMeanV0C[io] = 57.2673;  fMeanV0M[io] = 93.2262;  io++;
   fRuns[io] = 256311;  fMeanV0A[io] = 36.0997;  fMeanV0C[io] = 57.0719;  fMeanV0M[io] = 93.0143;  io++;
   fRuns[io] = 256356;  fMeanV0A[io] = 36.3877;  fMeanV0C[io] = 57.2140;  fMeanV0M[io] = 93.4760;  io++;
   fRuns[io] = 256361;  fMeanV0A[io] = 36.1659;  fMeanV0C[io] = 57.0607;  fMeanV0M[io] = 93.0680;  io++;
   fRuns[io] = 256362;  fMeanV0A[io] = 36.0275;  fMeanV0C[io] = 56.9170;  fMeanV0M[io] = 92.7952;  io++;
   fRuns[io] = 256363;  fMeanV0A[io] = 36.0008;  fMeanV0C[io] = 57.0383;  fMeanV0M[io] = 92.8851;  io++;
   fRuns[io] = 256364;  fMeanV0A[io] = 35.8209;  fMeanV0C[io] = 56.6544;  fMeanV0M[io] = 92.3170;  io++;
   fRuns[io] = 256365;  fMeanV0A[io] = 35.8589;  fMeanV0C[io] = 56.8648;  fMeanV0M[io] = 92.5704;  io++;
   fRuns[io] = 256366;  fMeanV0A[io] = 35.8088;  fMeanV0C[io] = 56.8071;  fMeanV0M[io] = 92.4614;  io++;
   fRuns[io] = 256368;  fMeanV0A[io] = 35.7529;  fMeanV0C[io] = 56.7149;  fMeanV0M[io] = 92.3115;  io++;
   fRuns[io] = 256371;  fMeanV0A[io] = 35.6694;  fMeanV0C[io] = 56.5066;  fMeanV0M[io] = 92.0214;  io++;
   fRuns[io] = 256372;  fMeanV0A[io] = 35.7052;  fMeanV0C[io] = 56.5651;  fMeanV0M[io] = 92.1074;  io++;
   fRuns[io] = 256373;  fMeanV0A[io] = 35.5747;  fMeanV0C[io] = 56.4439;  fMeanV0M[io] = 91.8633;  io++;
   fRuns[io] = 256415;  fMeanV0A[io] = 35.6585;  fMeanV0C[io] = 56.4481;  fMeanV0M[io] = 91.9484;  io++;
   fRuns[io] = 256417;  fMeanV0A[io] = 35.5611;  fMeanV0C[io] = 56.3188;  fMeanV0M[io] = 91.7255;  io++;
   fRuns[io] = 256418;  fMeanV0A[io] = 35.5078;  fMeanV0C[io] = 56.2743;  fMeanV0M[io] = 91.6465;  io++;

   // FILTER_p-p_208_LHC16k
   fRuns[io] = 256941;  fMeanV0A[io] = 35.8673;  fMeanV0C[io] = 57.0384;  fMeanV0M[io] = 92.7494;  io++;
   fRuns[io] = 256942;  fMeanV0A[io] = 35.9572;  fMeanV0C[io] = 57.3068;  fMeanV0M[io] = 93.1182;  io++;
   fRuns[io] = 256944;  fMeanV0A[io] = 35.9087;  fMeanV0C[io] = 57.1713;  fMeanV0M[io] = 92.8950;  io++;
   fRuns[io] = 257011;  fMeanV0A[io] = 36.1700;  fMeanV0C[io] = 57.8416;  fMeanV0M[io] = 93.8623;  io++;
   fRuns[io] = 257012;  fMeanV0A[io] = 36.1213;  fMeanV0C[io] = 57.6519;  fMeanV0M[io] = 93.6313;  io++;
   fRuns[io] = 257021;  fMeanV0A[io] = 35.6626;  fMeanV0C[io] = 57.2270;  fMeanV0M[io] = 92.7322;  io++;
   fRuns[io] = 257026;  fMeanV0A[io] = 35.7959;  fMeanV0C[io] = 57.1731;  fMeanV0M[io] = 92.8133;  io++;
   fRuns[io] = 257028;  fMeanV0A[io] = 35.3958;  fMeanV0C[io] = 56.7488;  fMeanV0M[io] = 91.9837;  io++;
   fRuns[io] = 257077;  fMeanV0A[io] = 35.2656;  fMeanV0C[io] = 56.4169;  fMeanV0M[io] = 91.5278;  io++;
   fRuns[io] = 257080;  fMeanV0A[io] = 35.0458;  fMeanV0C[io] = 56.2599;  fMeanV0M[io] = 91.1579;  io++;
   fRuns[io] = 257082;  fMeanV0A[io] = 35.0104;  fMeanV0C[io] = 56.2010;  fMeanV0M[io] = 91.0478;  io++;
   fRuns[io] = 257084;  fMeanV0A[io] = 34.9459;  fMeanV0C[io] = 56.1414;  fMeanV0M[io] = 90.9201;  io++;
   fRuns[io] = 257086;  fMeanV0A[io] = 35.1486;  fMeanV0C[io] = 56.2448;  fMeanV0M[io] = 91.2334;  io++;
   fRuns[io] = 257092;  fMeanV0A[io] = 35.0701;  fMeanV0C[io] = 56.0977;  fMeanV0M[io] = 91.0103;  io++;
   fRuns[io] = 257095;  fMeanV0A[io] = 35.1247;  fMeanV0C[io] = 56.1700;  fMeanV0M[io] = 91.1423;  io++;
   fRuns[io] = 257100;  fMeanV0A[io] = 35.1447;  fMeanV0C[io] = 56.1846;  fMeanV0M[io] = 91.1734;  io++;
   fRuns[io] = 257136;  fMeanV0A[io] = 35.1865;  fMeanV0C[io] = 56.2713;  fMeanV0M[io] = 91.3032;  io++;
   fRuns[io] = 257137;  fMeanV0A[io] = 35.1489;  fMeanV0C[io] = 56.4869;  fMeanV0M[io] = 91.4783;  io++;
   fRuns[io] = 257138;  fMeanV0A[io] = 35.0171;  fMeanV0C[io] = 56.3349;  fMeanV0M[io] = 91.1913;  io++;
   fRuns[io] = 257139;  fMeanV0A[io] = 35.0323;  fMeanV0C[io] = 55.9002;  fMeanV0M[io] = 90.7633;  io++;
   fRuns[io] = 257141;  fMeanV0A[io] = 34.9695;  fMeanV0C[io] = 56.1466;  fMeanV0M[io] = 90.9392;  io++;
   fRuns[io] = 257144;  fMeanV0A[io] = 35.0630;  fMeanV0C[io] = 56.2356;  fMeanV0M[io] = 91.1716;  io++;
   fRuns[io] = 257204;  fMeanV0A[io] = 35.2115;  fMeanV0C[io] = 56.1456;  fMeanV0M[io] = 91.2081;  io++;
   fRuns[io] = 257206;  fMeanV0A[io] = 35.0043;  fMeanV0C[io] = 55.9774;  fMeanV0M[io] = 90.8189;  io++;
   fRuns[io] = 257209;  fMeanV0A[io] = 34.9221;  fMeanV0C[io] = 56.1361;  fMeanV0M[io] = 90.9053;  io++;
   fRuns[io] = 257224;  fMeanV0A[io] = 34.8305;  fMeanV0C[io] = 56.0172;  fMeanV0M[io] = 90.6969;  io++;
   fRuns[io] = 257260;  fMeanV0A[io] = 35.1471;  fMeanV0C[io] = 56.0325;  fMeanV0M[io] = 91.0099;  io++;
   fRuns[io] = 257318;  fMeanV0A[io] = 35.1065;  fMeanV0C[io] = 56.3292;  fMeanV0M[io] = 91.3058;  io++;
   fRuns[io] = 257320;  fMeanV0A[io] = 35.0981;  fMeanV0C[io] = 56.0159;  fMeanV0M[io] = 90.9247;  io++;
   fRuns[io] = 257322;  fMeanV0A[io] = 35.2225;  fMeanV0C[io] = 56.2579;  fMeanV0M[io] = 91.3241;  io++;
   fRuns[io] = 257330;  fMeanV0A[io] = 34.9714;  fMeanV0C[io] = 56.2221;  fMeanV0M[io] = 91.0419;  io++;
   fRuns[io] = 257358;  fMeanV0A[io] = 35.1628;  fMeanV0C[io] = 56.4380;  fMeanV0M[io] = 91.4504;  io++;
   fRuns[io] = 257364;  fMeanV0A[io] = 35.1075;  fMeanV0C[io] = 56.2617;  fMeanV0M[io] = 91.2087;  io++;
   fRuns[io] = 257433;  fMeanV0A[io] = 35.1881;  fMeanV0C[io] = 56.4082;  fMeanV0M[io] = 91.4310;  io++;
   fRuns[io] = 257457;  fMeanV0A[io] = 35.1612;  fMeanV0C[io] = 56.4401;  fMeanV0M[io] = 91.4279;  io++;
   fRuns[io] = 257468;  fMeanV0A[io] = 35.2193;  fMeanV0C[io] = 56.2915;  fMeanV0M[io] = 91.3445;  io++;
   fRuns[io] = 257474;  fMeanV0A[io] = 34.9012;  fMeanV0C[io] = 56.1337;  fMeanV0M[io] = 90.8746;  io++;
   fRuns[io] = 257487;  fMeanV0A[io] = 34.7644;  fMeanV0C[io] = 55.9609;  fMeanV0M[io] = 90.5636;  io++;
   fRuns[io] = 257488;  fMeanV0A[io] = 34.7579;  fMeanV0C[io] = 55.9305;  fMeanV0M[io] = 90.5197;  io++;
   fRuns[io] = 257490;  fMeanV0A[io] = 34.7348;  fMeanV0C[io] = 55.7305;  fMeanV0M[io] = 90.2991;  io++;
   fRuns[io] = 257491;  fMeanV0A[io] = 34.5574;  fMeanV0C[io] = 55.7892;  fMeanV0M[io] = 90.1774;  io++;
   fRuns[io] = 257492;  fMeanV0A[io] = 34.4556;  fMeanV0C[io] = 55.4383;  fMeanV0M[io] = 89.7372;  io++;
   fRuns[io] = 257530;  fMeanV0A[io] = 34.7708;  fMeanV0C[io] = 55.7783;  fMeanV0M[io] = 90.3918;  io++;
   fRuns[io] = 257531;  fMeanV0A[io] = 34.6157;  fMeanV0C[io] = 55.6137;  fMeanV0M[io] = 90.0649;  io++;
   fRuns[io] = 257537;  fMeanV0A[io] = 34.6162;  fMeanV0C[io] = 55.8753;  fMeanV0M[io] = 90.3181;  io++;
   fRuns[io] = 257539;  fMeanV0A[io] = 34.2402;  fMeanV0C[io] = 55.1709;  fMeanV0M[io] = 89.2075;  io++;
   fRuns[io] = 257540;  fMeanV0A[io] = 34.5451;  fMeanV0C[io] = 55.4936;  fMeanV0M[io] = 89.8840;  io++;
   fRuns[io] = 257541;  fMeanV0A[io] = 34.4958;  fMeanV0C[io] = 55.6146;  fMeanV0M[io] = 89.9581;  io++;
   fRuns[io] = 257560;  fMeanV0A[io] = 34.4471;  fMeanV0C[io] = 55.5452;  fMeanV0M[io] = 89.8230;  io++;
   fRuns[io] = 257561;  fMeanV0A[io] = 34.3444;  fMeanV0C[io] = 55.2989;  fMeanV0M[io] = 89.4524;  io++;
   fRuns[io] = 257562;  fMeanV0A[io] = 34.4013;  fMeanV0C[io] = 55.3879;  fMeanV0M[io] = 89.6276;  io++;
   fRuns[io] = 257566;  fMeanV0A[io] = 34.2958;  fMeanV0C[io] = 55.1906;  fMeanV0M[io] = 89.3190;  io++;
   fRuns[io] = 257587;  fMeanV0A[io] = 34.4798;  fMeanV0C[io] = 55.5156;  fMeanV0M[io] = 89.8391;  io++;
   fRuns[io] = 257588;  fMeanV0A[io] = 34.4259;  fMeanV0C[io] = 55.3681;  fMeanV0M[io] = 89.6433;  io++;
   fRuns[io] = 257590;  fMeanV0A[io] = 34.4005;  fMeanV0C[io] = 55.3964;  fMeanV0M[io] = 89.6325;  io++;
   fRuns[io] = 257592;  fMeanV0A[io] = 34.0767;  fMeanV0C[io] = 55.1734;  fMeanV0M[io] = 89.0998;  io++;
   fRuns[io] = 257594;  fMeanV0A[io] = 34.1380;  fMeanV0C[io] = 55.0025;  fMeanV0M[io] = 88.9803;  io++;
   fRuns[io] = 257595;  fMeanV0A[io] = 34.1794;  fMeanV0C[io] = 55.1849;  fMeanV0M[io] = 89.2123;  io++;
   fRuns[io] = 257601;  fMeanV0A[io] = 34.0930;  fMeanV0C[io] = 55.1357;  fMeanV0M[io] = 89.0756;  io++;
   fRuns[io] = 257604;  fMeanV0A[io] = 34.1850;  fMeanV0C[io] = 55.0610;  fMeanV0M[io] = 89.0678;  io++;
   fRuns[io] = 257605;  fMeanV0A[io] = 34.1916;  fMeanV0C[io] = 55.1231;  fMeanV0M[io] = 89.1406;  io++;
   fRuns[io] = 257606;  fMeanV0A[io] = 34.1534;  fMeanV0C[io] = 55.0648;  fMeanV0M[io] = 89.0448;  io++;
   fRuns[io] = 257630;  fMeanV0A[io] = 34.1130;  fMeanV0C[io] = 55.0898;  fMeanV0M[io] = 89.0299;  io++;
   fRuns[io] = 257632;  fMeanV0A[io] = 34.0969;  fMeanV0C[io] = 54.9790;  fMeanV0M[io] = 88.9049;  io++;
   fRuns[io] = 257635;  fMeanV0A[io] = 34.0442;  fMeanV0C[io] = 55.0672;  fMeanV0M[io] = 88.9622;  io++;
   fRuns[io] = 257636;  fMeanV0A[io] = 34.0794;  fMeanV0C[io] = 55.0698;  fMeanV0M[io] = 88.9961;  io++;
   fRuns[io] = 257642;  fMeanV0A[io] = 34.0883;  fMeanV0C[io] = 55.0050;  fMeanV0M[io] = 88.9307;  io++;
   fRuns[io] = 257644;  fMeanV0A[io] = 34.1488;  fMeanV0C[io] = 55.1188;  fMeanV0M[io] = 89.1051;  io++;
   fRuns[io] = 257682;  fMeanV0A[io] = 34.2680;  fMeanV0C[io] = 55.2556;  fMeanV0M[io] = 89.3588;  io++;
   fRuns[io] = 257684;  fMeanV0A[io] = 34.1812;  fMeanV0C[io] = 55.1400;  fMeanV0M[io] = 89.1520;  io++;
   fRuns[io] = 257685;  fMeanV0A[io] = 34.0817;  fMeanV0C[io] = 55.0672;  fMeanV0M[io] = 89.0009;  io++;
   fRuns[io] = 257687;  fMeanV0A[io] = 34.0579;  fMeanV0C[io] = 54.9359;  fMeanV0M[io] = 88.8196;  io++;
   fRuns[io] = 257688;  fMeanV0A[io] = 34.0335;  fMeanV0C[io] = 54.9373;  fMeanV0M[io] = 88.8055;  io++;
   fRuns[io] = 257689;  fMeanV0A[io] = 34.0482;  fMeanV0C[io] = 54.8261;  fMeanV0M[io] = 88.7096;  io++;
   fRuns[io] = 257691;  fMeanV0A[io] = 34.0311;  fMeanV0C[io] = 55.1322;  fMeanV0M[io] = 88.9834;  io++;
   fRuns[io] = 257692;  fMeanV0A[io] = 33.9866;  fMeanV0C[io] = 54.8754;  fMeanV0M[io] = 88.6994;  io++;
   fRuns[io] = 257694;  fMeanV0A[io] = 33.8594;  fMeanV0C[io] = 54.7840;  fMeanV0M[io] = 88.4702;  io++;
   fRuns[io] = 257697;  fMeanV0A[io] = 34.0147;  fMeanV0C[io] = 54.9848;  fMeanV0M[io] = 88.8402;  io++;
   fRuns[io] = 257724;  fMeanV0A[io] = 34.2079;  fMeanV0C[io] = 54.9049;  fMeanV0M[io] = 88.9657;  io++;
   fRuns[io] = 257725;  fMeanV0A[io] = 34.0242;  fMeanV0C[io] = 54.8696;  fMeanV0M[io] = 88.7270;  io++;
   fRuns[io] = 257727;  fMeanV0A[io] = 33.9073;  fMeanV0C[io] = 54.7185;  fMeanV0M[io] = 88.4570;  io++;
   fRuns[io] = 257733;  fMeanV0A[io] = 33.9461;  fMeanV0C[io] = 54.7555;  fMeanV0M[io] = 88.5378;  io++;
   fRuns[io] = 257734;  fMeanV0A[io] = 33.8559;  fMeanV0C[io] = 54.7282;  fMeanV0M[io] = 88.3964;  io++;
   fRuns[io] = 257735;  fMeanV0A[io] = 33.9153;  fMeanV0C[io] = 54.7709;  fMeanV0M[io] = 88.5266;  io++;
   fRuns[io] = 257737;  fMeanV0A[io] = 33.9110;  fMeanV0C[io] = 54.7812;  fMeanV0M[io] = 88.5281;  io++;
   fRuns[io] = 257754;  fMeanV0A[io] = 34.0443;  fMeanV0C[io] = 54.8883;  fMeanV0M[io] = 88.7630;  io++;
   fRuns[io] = 257757;  fMeanV0A[io] = 33.8942;  fMeanV0C[io] = 54.7471;  fMeanV0M[io] = 88.4886;  io++;
   fRuns[io] = 257765;  fMeanV0A[io] = 33.8387;  fMeanV0C[io] = 54.6651;  fMeanV0M[io] = 88.3353;  io++;
   fRuns[io] = 257773;  fMeanV0A[io] = 33.7751;  fMeanV0C[io] = 54.4777;  fMeanV0M[io] = 88.0735;  io++;
   fRuns[io] = 257797;  fMeanV0A[io] = 33.9985;  fMeanV0C[io] = 54.8825;  fMeanV0M[io] = 88.7225;  io++;
   fRuns[io] = 257798;  fMeanV0A[io] = 33.9258;  fMeanV0C[io] = 54.8001;  fMeanV0M[io] = 88.5524;  io++;
   fRuns[io] = 257799;  fMeanV0A[io] = 33.8330;  fMeanV0C[io] = 54.6827;  fMeanV0M[io] = 88.3463;  io++;
   fRuns[io] = 257800;  fMeanV0A[io] = 33.8402;  fMeanV0C[io] = 54.6348;  fMeanV0M[io] = 88.3171;  io++;
   fRuns[io] = 257803;  fMeanV0A[io] = 33.7953;  fMeanV0C[io] = 54.5798;  fMeanV0M[io] = 88.2044;  io++;
   fRuns[io] = 257804;  fMeanV0A[io] = 33.7492;  fMeanV0C[io] = 54.5191;  fMeanV0M[io] = 88.0967;  io++;
   fRuns[io] = 257850;  fMeanV0A[io] = 34.0029;  fMeanV0C[io] = 54.9463;  fMeanV0M[io] = 88.7898;  io++;
   fRuns[io] = 257851;  fMeanV0A[io] = 33.9474;  fMeanV0C[io] = 54.7082;  fMeanV0M[io] = 88.4893;  io++;
   fRuns[io] = 257853;  fMeanV0A[io] = 33.9767;  fMeanV0C[io] = 54.8568;  fMeanV0M[io] = 88.6601;  io++;
   fRuns[io] = 257855;  fMeanV0A[io] = 33.9963;  fMeanV0C[io] = 54.8370;  fMeanV0M[io] = 88.6531;  io++;
   fRuns[io] = 257892;  fMeanV0A[io] = 34.2373;  fMeanV0C[io] = 55.0003;  fMeanV0M[io] = 89.0605;  io++;
   fRuns[io] = 257936;  fMeanV0A[io] = 33.8016;  fMeanV0C[io] = 54.6222;  fMeanV0M[io] = 88.2644;  io++;
   fRuns[io] = 257937;  fMeanV0A[io] = 33.8328;  fMeanV0C[io] = 54.5511;  fMeanV0M[io] = 88.2237;  io++;
   fRuns[io] = 257939;  fMeanV0A[io] = 33.8926;  fMeanV0C[io] = 54.7041;  fMeanV0M[io] = 88.4249;  io++;
   fRuns[io] = 257957;  fMeanV0A[io] = 34.0013;  fMeanV0C[io] = 54.6232;  fMeanV0M[io] = 88.4598;  io++;
   fRuns[io] = 257960;  fMeanV0A[io] = 33.8829;  fMeanV0C[io] = 54.7559;  fMeanV0M[io] = 88.4795;  io++;
   fRuns[io] = 257963;  fMeanV0A[io] = 33.9014;  fMeanV0C[io] = 54.8750;  fMeanV0M[io] = 88.5852;  io++;
   fRuns[io] = 257979;  fMeanV0A[io] = 33.7818;  fMeanV0C[io] = 54.5875;  fMeanV0M[io] = 88.1964;  io++;
   fRuns[io] = 257986;  fMeanV0A[io] = 33.8270;  fMeanV0C[io] = 54.7549;  fMeanV0M[io] = 88.4121;  io++;
   fRuns[io] = 257989;  fMeanV0A[io] = 33.7785;  fMeanV0C[io] = 54.5739;  fMeanV0M[io] = 88.1846;  io++;
   fRuns[io] = 257992;  fMeanV0A[io] = 33.7187;  fMeanV0C[io] = 54.5034;  fMeanV0M[io] = 88.0588;  io++;
   fRuns[io] = 258003;  fMeanV0A[io] = 33.6952;  fMeanV0C[io] = 54.4450;  fMeanV0M[io] = 87.9706;  io++;
   fRuns[io] = 258008;  fMeanV0A[io] = 33.8297;  fMeanV0C[io] = 54.4297;  fMeanV0M[io] = 88.0965;  io++;
   fRuns[io] = 258012;  fMeanV0A[io] = 33.7258;  fMeanV0C[io] = 54.4729;  fMeanV0M[io] = 88.0357;  io++;
   fRuns[io] = 258014;  fMeanV0A[io] = 33.6348;  fMeanV0C[io] = 54.1607;  fMeanV0M[io] = 87.6289;  io++;
   fRuns[io] = 258017;  fMeanV0A[io] = 33.8320;  fMeanV0C[io] = 54.4048;  fMeanV0M[io] = 88.0484;  io++;
   fRuns[io] = 258019;  fMeanV0A[io] = 33.7582;  fMeanV0C[io] = 54.3864;  fMeanV0M[io] = 87.9712;  io++;
   fRuns[io] = 258039;  fMeanV0A[io] = 33.7144;  fMeanV0C[io] = 54.3493;  fMeanV0M[io] = 87.8991;  io++;
   fRuns[io] = 258041;  fMeanV0A[io] = 33.6080;  fMeanV0C[io] = 54.3917;  fMeanV0M[io] = 87.8211;  io++;
   fRuns[io] = 258042;  fMeanV0A[io] = 33.7555;  fMeanV0C[io] = 54.4033;  fMeanV0M[io] = 87.9732;  io++;
   fRuns[io] = 258045;  fMeanV0A[io] = 33.5827;  fMeanV0C[io] = 54.2221;  fMeanV0M[io] = 87.6243;  io++;
   fRuns[io] = 258049;  fMeanV0A[io] = 33.6021;  fMeanV0C[io] = 54.1695;  fMeanV0M[io] = 87.5937;  io++;
   fRuns[io] = 258053;  fMeanV0A[io] = 33.7061;  fMeanV0C[io] = 54.2864;  fMeanV0M[io] = 87.8173;  io++;
   fRuns[io] = 258059;  fMeanV0A[io] = 33.5714;  fMeanV0C[io] = 54.1484;  fMeanV0M[io] = 87.5255;  io++;
   fRuns[io] = 258060;  fMeanV0A[io] = 33.5509;  fMeanV0C[io] = 54.0680;  fMeanV0M[io] = 87.4483;  io++;
   fRuns[io] = 258062;  fMeanV0A[io] = 33.6478;  fMeanV0C[io] = 54.2884;  fMeanV0M[io] = 87.7828;  io++;
   fRuns[io] = 258063;  fMeanV0A[io] = 33.6103;  fMeanV0C[io] = 54.2620;  fMeanV0M[io] = 87.7019;  io++;
   fRuns[io] = 258107;  fMeanV0A[io] = 33.6708;  fMeanV0C[io] = 54.2998;  fMeanV0M[io] = 87.7975;  io++;
   fRuns[io] = 258108;  fMeanV0A[io] = 33.5904;  fMeanV0C[io] = 54.1055;  fMeanV0M[io] = 87.5235;  io++;
   fRuns[io] = 258109;  fMeanV0A[io] = 33.5371;  fMeanV0C[io] = 53.9906;  fMeanV0M[io] = 87.3686;  io++;
   fRuns[io] = 258113;  fMeanV0A[io] = 33.5400;  fMeanV0C[io] = 54.0088;  fMeanV0M[io] = 87.3809;  io++;
   fRuns[io] = 258114;  fMeanV0A[io] = 33.4998;  fMeanV0C[io] = 54.0149;  fMeanV0M[io] = 87.3506;  io++;
   fRuns[io] = 258117;  fMeanV0A[io] = 33.5150;  fMeanV0C[io] = 54.1051;  fMeanV0M[io] = 87.4443;  io++;
   fRuns[io] = 258178;  fMeanV0A[io] = 33.6601;  fMeanV0C[io] = 54.2856;  fMeanV0M[io] = 87.7795;  io++;
   fRuns[io] = 258197;  fMeanV0A[io] = 33.6359;  fMeanV0C[io] = 54.2511;  fMeanV0M[io] = 87.7128;  io++;
   fRuns[io] = 258198;  fMeanV0A[io] = 33.4870;  fMeanV0C[io] = 54.1956;  fMeanV0M[io] = 87.5046;  io++;
   fRuns[io] = 258202;  fMeanV0A[io] = 33.5757;  fMeanV0C[io] = 54.2124;  fMeanV0M[io] = 87.6135;  io++;
   fRuns[io] = 258203;  fMeanV0A[io] = 33.4970;  fMeanV0C[io] = 54.0601;  fMeanV0M[io] = 87.3754;  io++;
   fRuns[io] = 258204;  fMeanV0A[io] = 33.5278;  fMeanV0C[io] = 54.1134;  fMeanV0M[io] = 87.4734;  io++;
   fRuns[io] = 258256;  fMeanV0A[io] = 33.6056;  fMeanV0C[io] = 54.1960;  fMeanV0M[io] = 87.6307;  io++;
   fRuns[io] = 258257;  fMeanV0A[io] = 33.6290;  fMeanV0C[io] = 54.2359;  fMeanV0M[io] = 87.7114;  io++;
   fRuns[io] = 258258;  fMeanV0A[io] = 33.6267;  fMeanV0C[io] = 54.3988;  fMeanV0M[io] = 87.8593;  io++;
   fRuns[io] = 258270;  fMeanV0A[io] = 33.4656;  fMeanV0C[io] = 54.0777;  fMeanV0M[io] = 87.3727;  io++;
   fRuns[io] = 258271;  fMeanV0A[io] = 33.4399;  fMeanV0C[io] = 54.1840;  fMeanV0M[io] = 87.4614;  io++;
   fRuns[io] = 258273;  fMeanV0A[io] = 33.5541;  fMeanV0C[io] = 54.1482;  fMeanV0M[io] = 87.5079;  io++;
   fRuns[io] = 258274;  fMeanV0A[io] = 33.4251;  fMeanV0C[io] = 54.1523;  fMeanV0M[io] = 87.4001;  io++;
   fRuns[io] = 258278;  fMeanV0A[io] = 33.6740;  fMeanV0C[io] = 54.1824;  fMeanV0M[io] = 87.6797;  io++;
   fRuns[io] = 258299;  fMeanV0A[io] = 33.8104;  fMeanV0C[io] = 54.3708;  fMeanV0M[io] = 88.0504;  io++;
   fRuns[io] = 258301;  fMeanV0A[io] = 33.5234;  fMeanV0C[io] = 54.1699;  fMeanV0M[io] = 87.5192;  io++;
   fRuns[io] = 258302;  fMeanV0A[io] = 33.4664;  fMeanV0C[io] = 54.0427;  fMeanV0M[io] = 87.3341;  io++;
   fRuns[io] = 258303;  fMeanV0A[io] = 33.4606;  fMeanV0C[io] = 54.1800;  fMeanV0M[io] = 87.4631;  io++;
   fRuns[io] = 258306;  fMeanV0A[io] = 33.3634;  fMeanV0C[io] = 54.0243;  fMeanV0M[io] = 87.2119;  io++;
   fRuns[io] = 258307;  fMeanV0A[io] = 33.3435;  fMeanV0C[io] = 53.9848;  fMeanV0M[io] = 87.1498;  io++;
   fRuns[io] = 258332;  fMeanV0A[io] = 33.5055;  fMeanV0C[io] = 54.2450;  fMeanV0M[io] = 87.5817;  io++;
   fRuns[io] = 258336;  fMeanV0A[io] = 33.4258;  fMeanV0C[io] = 53.9973;  fMeanV0M[io] = 87.2586;  io++;
   fRuns[io] = 258359;  fMeanV0A[io] = 33.4304;  fMeanV0C[io] = 54.0624;  fMeanV0M[io] = 87.3328;  io++;
   fRuns[io] = 258387;  fMeanV0A[io] = 33.4862;  fMeanV0C[io] = 53.6115;  fMeanV0M[io] = 86.9388;  io++;
   fRuns[io] = 258391;  fMeanV0A[io] = 33.3935;  fMeanV0C[io] = 53.9894;  fMeanV0M[io] = 87.2134;  io++;
   fRuns[io] = 258393;  fMeanV0A[io] = 33.2514;  fMeanV0C[io] = 53.8952;  fMeanV0M[io] = 86.9739;  io++;
   fRuns[io] = 258426;  fMeanV0A[io] = 33.2694;  fMeanV0C[io] = 54.1194;  fMeanV0M[io] = 87.2168;  io++;
   fRuns[io] = 258452;  fMeanV0A[io] = 33.3884;  fMeanV0C[io] = 54.1565;  fMeanV0M[io] = 87.3673;  io++;
   fRuns[io] = 258454;  fMeanV0A[io] = 33.2402;  fMeanV0C[io] = 53.7973;  fMeanV0M[io] = 86.8605;  io++;
   fRuns[io] = 258456;  fMeanV0A[io] = 33.1942;  fMeanV0C[io] = 53.9379;  fMeanV0M[io] = 86.9609;  io++;
   fRuns[io] = 258477;  fMeanV0A[io] = 33.1808;  fMeanV0C[io] = 53.8505;  fMeanV0M[io] = 86.8596;  io++;
   fRuns[io] = 258499;  fMeanV0A[io] = 33.2266;  fMeanV0C[io] = 53.7577;  fMeanV0M[io] = 86.8135;  io++;
   fRuns[io] = 258537;  fMeanV0A[io] = 33.3846;  fMeanV0C[io] = 54.1004;  fMeanV0M[io] = 87.3051;  io++;

   // FILTER_p-p_208_LHC16l
   fRuns[io] = 258919;  fMeanV0A[io] = 34.4914;  fMeanV0C[io] = 56.2004;  fMeanV0M[io] = 90.5395;  io++;
   fRuns[io] = 258923;  fMeanV0A[io] = 34.3178;  fMeanV0C[io] = 55.8782;  fMeanV0M[io] = 90.0409;  io++;
   fRuns[io] = 258962;  fMeanV0A[io] = 34.3707;  fMeanV0C[io] = 55.6343;  fMeanV0M[io] = 89.8506;  io++;
   fRuns[io] = 258964;  fMeanV0A[io] = 34.1137;  fMeanV0C[io] = 55.4591;  fMeanV0M[io] = 89.4036;  io++;
   fRuns[io] = 259088;  fMeanV0A[io] = 33.9508;  fMeanV0C[io] = 55.2812;  fMeanV0M[io] = 89.0623;  io++;
   fRuns[io] = 259090;  fMeanV0A[io] = 34.3084;  fMeanV0C[io] = 55.7677;  fMeanV0M[io] = 89.8728;  io++;
   fRuns[io] = 259091;  fMeanV0A[io] = 34.3212;  fMeanV0C[io] = 55.9455;  fMeanV0M[io] = 90.0905;  io++;
   fRuns[io] = 259096;  fMeanV0A[io] = 34.1366;  fMeanV0C[io] = 55.5434;  fMeanV0M[io] = 89.5387;  io++;
   fRuns[io] = 259099;  fMeanV0A[io] = 34.1129;  fMeanV0C[io] = 55.3375;  fMeanV0M[io] = 89.2860;  io++;
   fRuns[io] = 259117;  fMeanV0A[io] = 34.0526;  fMeanV0C[io] = 55.1900;  fMeanV0M[io] = 89.0938;  io++;
   fRuns[io] = 259118;  fMeanV0A[io] = 33.8573;  fMeanV0C[io] = 54.8806;  fMeanV0M[io] = 88.5843;  io++;
   fRuns[io] = 259162;  fMeanV0A[io] = 34.1587;  fMeanV0C[io] = 55.3975;  fMeanV0M[io] = 89.3769;  io++;
   fRuns[io] = 259204;  fMeanV0A[io] = 33.9156;  fMeanV0C[io] = 55.0741;  fMeanV0M[io] = 88.8159;  io++;
   fRuns[io] = 259257;  fMeanV0A[io] = 33.9806;  fMeanV0C[io] = 55.1582;  fMeanV0M[io] = 88.9619;  io++;
   fRuns[io] = 259261;  fMeanV0A[io] = 33.9649;  fMeanV0C[io] = 55.2691;  fMeanV0M[io] = 89.0716;  io++;
   fRuns[io] = 259263;  fMeanV0A[io] = 33.9173;  fMeanV0C[io] = 55.0013;  fMeanV0M[io] = 88.7473;  io++;
   fRuns[io] = 259264;  fMeanV0A[io] = 33.9868;  fMeanV0C[io] = 55.2807;  fMeanV0M[io] = 89.1018;  io++;
   fRuns[io] = 259269;  fMeanV0A[io] = 33.7863;  fMeanV0C[io] = 55.0210;  fMeanV0M[io] = 88.6334;  io++;
   fRuns[io] = 259270;  fMeanV0A[io] = 33.9441;  fMeanV0C[io] = 54.9004;  fMeanV0M[io] = 88.6752;  io++;
   fRuns[io] = 259271;  fMeanV0A[io] = 33.9592;  fMeanV0C[io] = 55.0255;  fMeanV0M[io] = 88.8104;  io++;
   fRuns[io] = 259272;  fMeanV0A[io] = 33.9610;  fMeanV0C[io] = 55.0699;  fMeanV0M[io] = 88.8690;  io++;
   fRuns[io] = 259273;  fMeanV0A[io] = 33.8020;  fMeanV0C[io] = 54.8484;  fMeanV0M[io] = 88.4843;  io++;
   fRuns[io] = 259274;  fMeanV0A[io] = 33.8162;  fMeanV0C[io] = 54.6482;  fMeanV0M[io] = 88.3091;  io++;
   fRuns[io] = 259302;  fMeanV0A[io] = 34.0061;  fMeanV0C[io] = 55.3003;  fMeanV0M[io] = 89.1493;  io++;
   fRuns[io] = 259303;  fMeanV0A[io] = 33.9418;  fMeanV0C[io] = 54.9118;  fMeanV0M[io] = 88.6828;  io++;
   fRuns[io] = 259305;  fMeanV0A[io] = 33.7853;  fMeanV0C[io] = 54.8316;  fMeanV0M[io] = 88.4500;  io++;
   fRuns[io] = 259307;  fMeanV0A[io] = 33.7006;  fMeanV0C[io] = 54.6174;  fMeanV0M[io] = 88.1527;  io++;
   fRuns[io] = 259334;  fMeanV0A[io] = 33.7397;  fMeanV0C[io] = 54.8103;  fMeanV0M[io] = 88.3796;  io++;
   fRuns[io] = 259336;  fMeanV0A[io] = 33.7724;  fMeanV0C[io] = 54.8777;  fMeanV0M[io] = 88.4803;  io++;
   fRuns[io] = 259339;  fMeanV0A[io] = 33.7440;  fMeanV0C[io] = 54.9377;  fMeanV0M[io] = 88.5151;  io++;
   fRuns[io] = 259340;  fMeanV0A[io] = 33.7349;  fMeanV0C[io] = 54.6236;  fMeanV0M[io] = 88.1999;  io++;
   fRuns[io] = 259341;  fMeanV0A[io] = 33.7142;  fMeanV0C[io] = 54.7088;  fMeanV0M[io] = 88.2560;  io++;
   fRuns[io] = 259342;  fMeanV0A[io] = 33.7241;  fMeanV0C[io] = 54.7186;  fMeanV0M[io] = 88.2702;  io++;
   fRuns[io] = 259378;  fMeanV0A[io] = 33.7283;  fMeanV0C[io] = 54.5212;  fMeanV0M[io] = 88.0738;  io++;
   fRuns[io] = 259382;  fMeanV0A[io] = 33.7264;  fMeanV0C[io] = 54.5205;  fMeanV0M[io] = 88.0848;  io++;
   fRuns[io] = 259388;  fMeanV0A[io] = 33.6508;  fMeanV0C[io] = 54.7131;  fMeanV0M[io] = 88.2044;  io++;
   fRuns[io] = 259389;  fMeanV0A[io] = 33.6353;  fMeanV0C[io] = 54.7209;  fMeanV0M[io] = 88.1889;  io++;
   fRuns[io] = 259394;  fMeanV0A[io] = 33.5886;  fMeanV0C[io] = 54.5773;  fMeanV0M[io] = 87.9949;  io++;
   fRuns[io] = 259395;  fMeanV0A[io] = 33.5998;  fMeanV0C[io] = 54.6285;  fMeanV0M[io] = 88.0468;  io++;
   fRuns[io] = 259396;  fMeanV0A[io] = 33.6836;  fMeanV0C[io] = 54.5338;  fMeanV0M[io] = 88.0419;  io++;
   fRuns[io] = 259473;  fMeanV0A[io] = 33.7134;  fMeanV0C[io] = 54.4865;  fMeanV0M[io] = 88.0632;  io++;
   fRuns[io] = 259477;  fMeanV0A[io] = 33.6195;  fMeanV0C[io] = 54.5453;  fMeanV0M[io] = 87.9935;  io++;
   fRuns[io] = 259747;  fMeanV0A[io] = 33.7514;  fMeanV0C[io] = 54.9986;  fMeanV0M[io] = 88.5649;  io++;
   fRuns[io] = 259748;  fMeanV0A[io] = 33.3619;  fMeanV0C[io] = 54.4581;  fMeanV0M[io] = 87.6605;  io++;
   fRuns[io] = 259750;  fMeanV0A[io] = 33.7655;  fMeanV0C[io] = 54.7701;  fMeanV0M[io] = 88.3708;  io++;
   fRuns[io] = 259751;  fMeanV0A[io] = 33.8231;  fMeanV0C[io] = 54.8894;  fMeanV0M[io] = 88.5492;  io++;
   fRuns[io] = 259752;  fMeanV0A[io] = 33.5726;  fMeanV0C[io] = 54.6453;  fMeanV0M[io] = 88.0599;  io++;
   fRuns[io] = 259756;  fMeanV0A[io] = 33.6019;  fMeanV0C[io] = 54.6708;  fMeanV0M[io] = 88.1054;  io++;
   fRuns[io] = 259781;  fMeanV0A[io] = 33.8090;  fMeanV0C[io] = 54.8454;  fMeanV0M[io] = 88.4894;  io++;
   fRuns[io] = 259788;  fMeanV0A[io] = 33.7429;  fMeanV0C[io] = 54.3738;  fMeanV0M[io] = 87.9595;  io++;
   fRuns[io] = 259789;  fMeanV0A[io] = 33.4697;  fMeanV0C[io] = 54.3261;  fMeanV0M[io] = 87.6266;  io++;
   fRuns[io] = 259822;  fMeanV0A[io] = 33.7797;  fMeanV0C[io] = 54.6678;  fMeanV0M[io] = 88.2951;  io++;
   fRuns[io] = 259841;  fMeanV0A[io] = 33.5933;  fMeanV0C[io] = 54.6251;  fMeanV0M[io] = 88.0475;  io++;
   fRuns[io] = 259842;  fMeanV0A[io] = 33.5941;  fMeanV0C[io] = 54.7403;  fMeanV0M[io] = 88.1638;  io++;
   fRuns[io] = 259860;  fMeanV0A[io] = 33.5704;  fMeanV0C[io] = 54.6549;  fMeanV0M[io] = 88.0547;  io++;
   fRuns[io] = 259866;  fMeanV0A[io] = 33.4203;  fMeanV0C[io] = 54.4165;  fMeanV0M[io] = 87.6564;  io++;
   fRuns[io] = 259867;  fMeanV0A[io] = 33.3623;  fMeanV0C[io] = 54.2432;  fMeanV0M[io] = 87.4240;  io++;
   fRuns[io] = 259868;  fMeanV0A[io] = 33.3224;  fMeanV0C[io] = 54.2497;  fMeanV0M[io] = 87.3971;  io++;
   fRuns[io] = 259888;  fMeanV0A[io] = 33.4895;  fMeanV0C[io] = 54.6938;  fMeanV0M[io] = 88.0095;  io++;

   // FILTER_p-p_208_LHC16o
   fRuns[io] = 262424;  fMeanV0A[io] = 33.1716;  fMeanV0C[io] = 54.5147;  fMeanV0M[io] = 87.5309;  io++;
   fRuns[io] = 262425;  fMeanV0A[io] = 32.8117;  fMeanV0C[io] = 53.9649;  fMeanV0M[io] = 86.6110;  io++;
   fRuns[io] = 262426;  fMeanV0A[io] = 32.8542;  fMeanV0C[io] = 53.9821;  fMeanV0M[io] = 86.6766;  io++;
   fRuns[io] = 262428;  fMeanV0A[io] = 32.7779;  fMeanV0C[io] = 53.9545;  fMeanV0M[io] = 86.5628;  io++;
   fRuns[io] = 262705;  fMeanV0A[io] = 32.3286;  fMeanV0C[io] = 52.7626;  fMeanV0M[io] = 84.9150;  io++;
   fRuns[io] = 262706;  fMeanV0A[io] = 32.2742;  fMeanV0C[io] = 52.7204;  fMeanV0M[io] = 84.8204;  io++;
   fRuns[io] = 262708;  fMeanV0A[io] = 32.3069;  fMeanV0C[io] = 52.7254;  fMeanV0M[io] = 84.8655;  io++;
   fRuns[io] = 262713;  fMeanV0A[io] = 32.2850;  fMeanV0C[io] = 52.6470;  fMeanV0M[io] = 84.7681;  io++;
   fRuns[io] = 262717;  fMeanV0A[io] = 32.3547;  fMeanV0C[io] = 52.6694;  fMeanV0M[io] = 84.8592;  io++;
   fRuns[io] = 262719;  fMeanV0A[io] = 32.3753;  fMeanV0C[io] = 52.4809;  fMeanV0M[io] = 84.6828;  io++;
   fRuns[io] = 262723;  fMeanV0A[io] = 32.3669;  fMeanV0C[io] = 52.6110;  fMeanV0M[io] = 84.8077;  io++;
   fRuns[io] = 262725;  fMeanV0A[io] = 32.2609;  fMeanV0C[io] = 52.5619;  fMeanV0M[io] = 84.6476;  io++;
   fRuns[io] = 262727;  fMeanV0A[io] = 32.3546;  fMeanV0C[io] = 52.5115;  fMeanV0M[io] = 84.6925;  io++;
   fRuns[io] = 262760;  fMeanV0A[io] = 32.5369;  fMeanV0C[io] = 52.7995;  fMeanV0M[io] = 85.1643;  io++;
   fRuns[io] = 262768;  fMeanV0A[io] = 32.4218;  fMeanV0C[io] = 52.8399;  fMeanV0M[io] = 85.0956;  io++;
   fRuns[io] = 262776;  fMeanV0A[io] = 32.3452;  fMeanV0C[io] = 52.5163;  fMeanV0M[io] = 84.6953;  io++;
   fRuns[io] = 262777;  fMeanV0A[io] = 32.2278;  fMeanV0C[io] = 52.3482;  fMeanV0M[io] = 84.3866;  io++;
   fRuns[io] = 262778;  fMeanV0A[io] = 32.2641;  fMeanV0C[io] = 52.5005;  fMeanV0M[io] = 84.5925;  io++;
   fRuns[io] = 262841;  fMeanV0A[io] = 32.4048;  fMeanV0C[io] = 52.5091;  fMeanV0M[io] = 84.6816;  io++;
   fRuns[io] = 262842;  fMeanV0A[io] = 32.2539;  fMeanV0C[io] = 52.6869;  fMeanV0M[io] = 84.7702;  io++;
   fRuns[io] = 262844;  fMeanV0A[io] = 32.3652;  fMeanV0C[io] = 52.6189;  fMeanV0M[io] = 84.8219;  io++;
   fRuns[io] = 262847;  fMeanV0A[io] = 32.2984;  fMeanV0C[io] = 52.7030;  fMeanV0M[io] = 84.8430;  io++;
   fRuns[io] = 262849;  fMeanV0A[io] = 32.3052;  fMeanV0C[io] = 52.7582;  fMeanV0M[io] = 84.8942;  io++;
   fRuns[io] = 262853;  fMeanV0A[io] = 32.2402;  fMeanV0C[io] = 52.6325;  fMeanV0M[io] = 84.6874;  io++;
   fRuns[io] = 262855;  fMeanV0A[io] = 32.3658;  fMeanV0C[io] = 52.7268;  fMeanV0M[io] = 84.9279;  io++;
   fRuns[io] = 262858;  fMeanV0A[io] = 32.3234;  fMeanV0C[io] = 52.6583;  fMeanV0M[io] = 84.8160;  io++;
   fRuns[io] = 263331;  fMeanV0A[io] = 32.9629;  fMeanV0C[io] = 53.8127;  fMeanV0M[io] = 86.6068;  io++;
   fRuns[io] = 263332;  fMeanV0A[io] = 33.1905;  fMeanV0C[io] = 53.6835;  fMeanV0M[io] = 86.7065;  io++;
   fRuns[io] = 263487;  fMeanV0A[io] = 33.0687;  fMeanV0C[io] = 53.4993;  fMeanV0M[io] = 86.3999;  io++;
   fRuns[io] = 263490;  fMeanV0A[io] = 32.9641;  fMeanV0C[io] = 53.4498;  fMeanV0M[io] = 86.2472;  io++;
   fRuns[io] = 263496;  fMeanV0A[io] = 33.0050;  fMeanV0C[io] = 53.6025;  fMeanV0M[io] = 86.4407;  io++;
   fRuns[io] = 263497;  fMeanV0A[io] = 32.9627;  fMeanV0C[io] = 53.5156;  fMeanV0M[io] = 86.3024;  io++;
   fRuns[io] = 263529;  fMeanV0A[io] = 32.8766;  fMeanV0C[io] = 53.2303;  fMeanV0M[io] = 85.9221;  io++;
   fRuns[io] = 263647;  fMeanV0A[io] = 32.9041;  fMeanV0C[io] = 53.2547;  fMeanV0M[io] = 85.9968;  io++;
   fRuns[io] = 263652;  fMeanV0A[io] = 32.8304;  fMeanV0C[io] = 53.1841;  fMeanV0M[io] = 85.8458;  io++;
   fRuns[io] = 263654;  fMeanV0A[io] = 32.7514;  fMeanV0C[io] = 52.9946;  fMeanV0M[io] = 85.5660;  io++;
   fRuns[io] = 263657;  fMeanV0A[io] = 32.7581;  fMeanV0C[io] = 53.2443;  fMeanV0M[io] = 85.8239;  io++;
   fRuns[io] = 263662;  fMeanV0A[io] = 32.7050;  fMeanV0C[io] = 52.9973;  fMeanV0M[io] = 85.5338;  io++;
   fRuns[io] = 263663;  fMeanV0A[io] = 32.6559;  fMeanV0C[io] = 52.8324;  fMeanV0M[io] = 85.3227;  io++;
   fRuns[io] = 263682;  fMeanV0A[io] = 32.6113;  fMeanV0C[io] = 52.9237;  fMeanV0M[io] = 85.3639;  io++;
   fRuns[io] = 263690;  fMeanV0A[io] = 32.5205;  fMeanV0C[io] = 52.7757;  fMeanV0M[io] = 85.1197;  io++;
   fRuns[io] = 263691;  fMeanV0A[io] = 32.4299;  fMeanV0C[io] = 52.6707;  fMeanV0M[io] = 84.9209;  io++;
   fRuns[io] = 263737;  fMeanV0A[io] = 32.6375;  fMeanV0C[io] = 52.7276;  fMeanV0M[io] = 85.1828;  io++;
   fRuns[io] = 263738;  fMeanV0A[io] = 32.8049;  fMeanV0C[io] = 52.9717;  fMeanV0M[io] = 85.5888;  io++;
   fRuns[io] = 263739;  fMeanV0A[io] = 32.4961;  fMeanV0C[io] = 52.7274;  fMeanV0M[io] = 85.0564;  io++;
   fRuns[io] = 263741;  fMeanV0A[io] = 32.4759;  fMeanV0C[io] = 52.6980;  fMeanV0M[io] = 84.9944;  io++;
   fRuns[io] = 263743;  fMeanV0A[io] = 32.3602;  fMeanV0C[io] = 52.6006;  fMeanV0M[io] = 84.7846;  io++;
   fRuns[io] = 263744;  fMeanV0A[io] = 32.5137;  fMeanV0C[io] = 52.8012;  fMeanV0M[io] = 85.1538;  io++;
   fRuns[io] = 263784;  fMeanV0A[io] = 32.6565;  fMeanV0C[io] = 52.9577;  fMeanV0M[io] = 85.4295;  io++;
   fRuns[io] = 263785;  fMeanV0A[io] = 32.4930;  fMeanV0C[io] = 52.7411;  fMeanV0M[io] = 85.0441;  io++;
   fRuns[io] = 263786;  fMeanV0A[io] = 32.4532;  fMeanV0C[io] = 52.7558;  fMeanV0M[io] = 85.0333;  io++;
   fRuns[io] = 263787;  fMeanV0A[io] = 32.3531;  fMeanV0C[io] = 52.6168;  fMeanV0M[io] = 84.7975;  io++;
   fRuns[io] = 263790;  fMeanV0A[io] = 32.2899;  fMeanV0C[io] = 52.6865;  fMeanV0M[io] = 84.7761;  io++;
   fRuns[io] = 263792;  fMeanV0A[io] = 32.3756;  fMeanV0C[io] = 52.5329;  fMeanV0M[io] = 84.7348;  io++;
   fRuns[io] = 263793;  fMeanV0A[io] = 32.3352;  fMeanV0C[io] = 52.7306;  fMeanV0M[io] = 84.8872;  io++;
   fRuns[io] = 263803;  fMeanV0A[io] = 32.3789;  fMeanV0C[io] = 52.5263;  fMeanV0M[io] = 84.7635;  io++;
   fRuns[io] = 263810;  fMeanV0A[io] = 32.4503;  fMeanV0C[io] = 52.6851;  fMeanV0M[io] = 84.9668;  io++;
   fRuns[io] = 263863;  fMeanV0A[io] = 32.2071;  fMeanV0C[io] = 52.1904;  fMeanV0M[io] = 84.2348;  io++;
   fRuns[io] = 263866;  fMeanV0A[io] = 32.1510;  fMeanV0C[io] = 52.2407;  fMeanV0M[io] = 84.1977;  io++;
   fRuns[io] = 263905;  fMeanV0A[io] = 32.2192;  fMeanV0C[io] = 52.3454;  fMeanV0M[io] = 84.3739;  io++;
   fRuns[io] = 263916;  fMeanV0A[io] = 32.1949;  fMeanV0C[io] = 52.4044;  fMeanV0M[io] = 84.4249;  io++;
   fRuns[io] = 263917;  fMeanV0A[io] = 32.2780;  fMeanV0C[io] = 52.7188;  fMeanV0M[io] = 84.8265;  io++;
   fRuns[io] = 263920;  fMeanV0A[io] = 32.0175;  fMeanV0C[io] = 52.3409;  fMeanV0M[io] = 84.1820;  io++;
   fRuns[io] = 263923;  fMeanV0A[io] = 32.3544;  fMeanV0C[io] = 52.2223;  fMeanV0M[io] = 84.3665;  io++;
   fRuns[io] = 263977;  fMeanV0A[io] = 32.2402;  fMeanV0C[io] = 52.4428;  fMeanV0M[io] = 84.4877;  io++;
   fRuns[io] = 263978;  fMeanV0A[io] = 32.0426;  fMeanV0C[io] = 52.2269;  fMeanV0M[io] = 84.0745;  io++;
   fRuns[io] = 263981;  fMeanV0A[io] = 32.1343;  fMeanV0C[io] = 52.2137;  fMeanV0M[io] = 84.1885;  io++;
   fRuns[io] = 263984;  fMeanV0A[io] = 32.0159;  fMeanV0C[io] = 52.2136;  fMeanV0M[io] = 84.0531;  io++;
   fRuns[io] = 263985;  fMeanV0A[io] = 31.9996;  fMeanV0C[io] = 52.2512;  fMeanV0M[io] = 84.0921;  io++;
   fRuns[io] = 264033;  fMeanV0A[io] = 31.8216;  fMeanV0C[io] = 52.0516;  fMeanV0M[io] = 83.6944;  io++;
   fRuns[io] = 264035;  fMeanV0A[io] = 31.7032;  fMeanV0C[io] = 51.7066;  fMeanV0M[io] = 83.1755;  io++;

   // FILTER_p-p_208_LHC16p
   fRuns[io] = 264076;  fMeanV0A[io] = 54.0549;  fMeanV0C[io] = 83.1567;  fMeanV0M[io] = 137.138;  io++;
   fRuns[io] = 264078;  fMeanV0A[io] = 54.1113;  fMeanV0C[io] = 83.3221;  fMeanV0M[io] = 137.356;  io++;
   fRuns[io] = 264082;  fMeanV0A[io] = 54.0991;  fMeanV0C[io] = 83.1547;  fMeanV0M[io] = 137.183;  io++;
   fRuns[io] = 264085;  fMeanV0A[io] = 54.1647;  fMeanV0C[io] = 83.0550;  fMeanV0M[io] = 137.141;  io++;
   fRuns[io] = 264086;  fMeanV0A[io] = 54.3115;  fMeanV0C[io] = 83.5542;  fMeanV0M[io] = 137.803;  io++;
   fRuns[io] = 264109;  fMeanV0A[io] = 54.1899;  fMeanV0C[io] = 83.2636;  fMeanV0M[io] = 137.375;  io++;
   fRuns[io] = 264110;  fMeanV0A[io] = 54.0797;  fMeanV0C[io] = 83.2679;  fMeanV0M[io] = 137.286;  io++;
   fRuns[io] = 264129;  fMeanV0A[io] = 54.1224;  fMeanV0C[io] = 83.1708;  fMeanV0M[io] = 137.216;  io++;
   fRuns[io] = 264137;  fMeanV0A[io] = 54.2878;  fMeanV0C[io] = 83.2537;  fMeanV0M[io] = 137.456;  io++;
   fRuns[io] = 264138;  fMeanV0A[io] = 53.8682;  fMeanV0C[io] = 83.1692;  fMeanV0M[io] = 136.948;  io++;
   fRuns[io] = 264139;  fMeanV0A[io] = 53.9860;  fMeanV0C[io] = 83.0464;  fMeanV0M[io] = 136.956;  io++;
   fRuns[io] = 264164;  fMeanV0A[io] = 53.9944;  fMeanV0C[io] = 83.2151;  fMeanV0M[io] = 137.128;  io++;
   fRuns[io] = 264168;  fMeanV0A[io] = 53.8682;  fMeanV0C[io] = 82.8789;  fMeanV0M[io] = 136.668;  io++;
   fRuns[io] = 264188;  fMeanV0A[io] = 53.8533;  fMeanV0C[io] = 82.8134;  fMeanV0M[io] = 136.591;  io++;
   fRuns[io] = 264190;  fMeanV0A[io] = 53.7856;  fMeanV0C[io] = 82.9129;  fMeanV0M[io] = 136.624;  io++;
   fRuns[io] = 264194;  fMeanV0A[io] = 53.7871;  fMeanV0C[io] = 82.9461;  fMeanV0M[io] = 136.664;  io++;
   fRuns[io] = 264197;  fMeanV0A[io] = 53.9310;  fMeanV0C[io] = 83.0801;  fMeanV0M[io] = 136.924;  io++;
   fRuns[io] = 264198;  fMeanV0A[io] = 53.8011;  fMeanV0C[io] = 82.6978;  fMeanV0M[io] = 136.424;  io++;
   fRuns[io] = 264232;  fMeanV0A[io] = 53.9393;  fMeanV0C[io] = 82.8527;  fMeanV0M[io] = 136.721;  io++;
   fRuns[io] = 264233;  fMeanV0A[io] = 53.8217;  fMeanV0C[io] = 82.7621;  fMeanV0M[io] = 136.502;  io++;
   fRuns[io] = 264235;  fMeanV0A[io] = 53.9278;  fMeanV0C[io] = 82.8783;  fMeanV0M[io] = 136.719;  io++;
   fRuns[io] = 264238;  fMeanV0A[io] = 53.9581;  fMeanV0C[io] = 82.9749;  fMeanV0M[io] = 136.844;  io++;
   fRuns[io] = 264259;  fMeanV0A[io] = 53.9508;  fMeanV0C[io] = 82.8556;  fMeanV0M[io] = 136.708;  io++;
   fRuns[io] = 264260;  fMeanV0A[io] = 53.8180;  fMeanV0C[io] = 82.5797;  fMeanV0M[io] = 136.290;  io++;
   fRuns[io] = 264261;  fMeanV0A[io] = 53.9529;  fMeanV0C[io] = 82.6324;  fMeanV0M[io] = 136.510;  io++;
   fRuns[io] = 264262;  fMeanV0A[io] = 53.9211;  fMeanV0C[io] = 83.0578;  fMeanV0M[io] = 136.886;  io++;
   fRuns[io] = 264264;  fMeanV0A[io] = 53.9748;  fMeanV0C[io] = 82.8723;  fMeanV0M[io] = 136.752;  io++;
   fRuns[io] = 264265;  fMeanV0A[io] = 53.9994;  fMeanV0C[io] = 82.8057;  fMeanV0M[io] = 136.724;  io++;
   fRuns[io] = 264266;  fMeanV0A[io] = 53.9268;  fMeanV0C[io] = 82.8188;  fMeanV0M[io] = 136.657;  io++;
   fRuns[io] = 264267;  fMeanV0A[io] = 53.9356;  fMeanV0C[io] = 82.8281;  fMeanV0M[io] = 136.689;  io++;
   fRuns[io] = 264273;  fMeanV0A[io] = 53.4152;  fMeanV0C[io] = 82.4261;  fMeanV0M[io] = 135.773;  io++;
   fRuns[io] = 264277;  fMeanV0A[io] = 53.7031;  fMeanV0C[io] = 82.9018;  fMeanV0M[io] = 136.534;  io++;
   fRuns[io] = 264279;  fMeanV0A[io] = 53.7084;  fMeanV0C[io] = 82.5741;  fMeanV0M[io] = 136.207;  io++;
   fRuns[io] = 264281;  fMeanV0A[io] = 53.8219;  fMeanV0C[io] = 83.1079;  fMeanV0M[io] = 136.864;  io++;
   fRuns[io] = 264305;  fMeanV0A[io] = 54.5496;  fMeanV0C[io] = 83.6104;  fMeanV0M[io] = 138.109;  io++;
   fRuns[io] = 264306;  fMeanV0A[io] = 54.0384;  fMeanV0C[io] = 83.0286;  fMeanV0M[io] = 136.981;  io++;
   fRuns[io] = 264312;  fMeanV0A[io] = 53.6803;  fMeanV0C[io] = 82.5010;  fMeanV0M[io] = 136.069;  io++;
   fRuns[io] = 264336;  fMeanV0A[io] = 53.7493;  fMeanV0C[io] = 82.6942;  fMeanV0M[io] = 136.364;  io++;
   fRuns[io] = 264341;  fMeanV0A[io] = 53.6894;  fMeanV0C[io] = 82.1597;  fMeanV0M[io] = 135.771;  io++;
   fRuns[io] = 264345;  fMeanV0A[io] = 53.7531;  fMeanV0C[io] = 82.6611;  fMeanV0M[io] = 136.344;  io++;
   fRuns[io] = 264346;  fMeanV0A[io] = 53.6183;  fMeanV0C[io] = 82.7642;  fMeanV0M[io] = 136.299;  io++;
   fRuns[io] = 264347;  fMeanV0A[io] = 53.7273;  fMeanV0C[io] = 82.6561;  fMeanV0M[io] = 136.300;  io++;

   // FILTER_p-p_208_LHC17c
   fRuns[io] = 270581;  fMeanV0A[io] = 53.6741;  fMeanV0C[io] = 79.6333;  fMeanV0M[io] = 133.222;  io++;
   fRuns[io] = 270661;  fMeanV0A[io] = 53.3540;  fMeanV0C[io] = 79.2331;  fMeanV0M[io] = 132.505;  io++;
   fRuns[io] = 270663;  fMeanV0A[io] = 53.0106;  fMeanV0C[io] = 78.8799;  fMeanV0M[io] = 131.810;  io++;
   fRuns[io] = 270665;  fMeanV0A[io] = 52.8295;  fMeanV0C[io] = 78.7320;  fMeanV0M[io] = 131.464;  io++;
   fRuns[io] = 270667;  fMeanV0A[io] = 53.3266;  fMeanV0C[io] = 79.1771;  fMeanV0M[io] = 132.414;  io++;

   // FILTER_p-p_208_LHC17e
   fRuns[io] = 270822;  fMeanV0A[io] = 53.4779;  fMeanV0C[io] = 79.1327;  fMeanV0M[io] = 132.520;  io++;
   fRuns[io] = 270824;  fMeanV0A[io] = 53.9797;  fMeanV0C[io] = 79.6179;  fMeanV0M[io] = 133.516;  io++;
   fRuns[io] = 270827;  fMeanV0A[io] = 54.1372;  fMeanV0C[io] = 79.8064;  fMeanV0M[io] = 133.865;  io++;
   fRuns[io] = 270828;  fMeanV0A[io] = 53.4874;  fMeanV0C[io] = 79.4045;  fMeanV0M[io] = 132.802;  io++;
   fRuns[io] = 270830;  fMeanV0A[io] = 53.3119;  fMeanV0C[io] = 78.7494;  fMeanV0M[io] = 131.978;  io++;

   // FILTER_p-p_208_LHC17f
   fRuns[io] = 270854;  fMeanV0A[io] = 54.2001;  fMeanV0C[io] = 78.8550;  fMeanV0M[io] = 132.966;  io++;
   fRuns[io] = 270855;  fMeanV0A[io] = 53.9221;  fMeanV0C[io] = 78.5313;  fMeanV0M[io] = 132.373;  io++;
   fRuns[io] = 270856;  fMeanV0A[io] = 53.7244;  fMeanV0C[io] = 78.5810;  fMeanV0M[io] = 132.228;  io++;
   fRuns[io] = 270861;  fMeanV0A[io] = 53.6014;  fMeanV0C[io] = 78.4466;  fMeanV0M[io] = 131.965;  io++;
   fRuns[io] = 270865;  fMeanV0A[io] = 52.9113;  fMeanV0C[io] = 77.5194;  fMeanV0M[io] = 130.346;  io++;

   // FILTER_p-p_208_LHC17h
   fRuns[io] = 271870;  fMeanV0A[io] = 53.4253;  fMeanV0C[io] = 77.6433;  fMeanV0M[io] = 130.985;  io++;
   fRuns[io] = 271871;  fMeanV0A[io] = 53.3939;  fMeanV0C[io] = 77.6289;  fMeanV0M[io] = 130.938;  io++;
   fRuns[io] = 271873;  fMeanV0A[io] = 53.4553;  fMeanV0C[io] = 77.8687;  fMeanV0M[io] = 131.239;  io++;
   fRuns[io] = 271874;  fMeanV0A[io] = 53.4622;  fMeanV0C[io] = 77.8199;  fMeanV0M[io] = 131.199;  io++;
   fRuns[io] = 271880;  fMeanV0A[io] = 53.3454;  fMeanV0C[io] = 77.5936;  fMeanV0M[io] = 130.845;  io++;
   fRuns[io] = 271886;  fMeanV0A[io] = 53.2773;  fMeanV0C[io] = 77.5834;  fMeanV0M[io] = 130.777;  io++;
   fRuns[io] = 272018;  fMeanV0A[io] = 53.1846;  fMeanV0C[io] = 77.1358;  fMeanV0M[io] = 130.244;  io++;
   fRuns[io] = 272020;  fMeanV0A[io] = 53.0868;  fMeanV0C[io] = 76.9609;  fMeanV0M[io] = 129.955;  io++;
   fRuns[io] = 272036;  fMeanV0A[io] = 52.8818;  fMeanV0C[io] = 76.2678;  fMeanV0M[io] = 129.064;  io++;
   fRuns[io] = 272038;  fMeanV0A[io] = 52.5681;  fMeanV0C[io] = 76.1117;  fMeanV0M[io] = 128.598;  io++;
   fRuns[io] = 272039;  fMeanV0A[io] = 52.4046;  fMeanV0C[io] = 76.0310;  fMeanV0M[io] = 128.346;  io++;
   fRuns[io] = 272040;  fMeanV0A[io] = 52.4014;  fMeanV0C[io] = 76.0261;  fMeanV0M[io] = 128.340;  io++;
   fRuns[io] = 272042;  fMeanV0A[io] = 52.1931;  fMeanV0C[io] = 75.8957;  fMeanV0M[io] = 128.000;  io++;
   fRuns[io] = 272076;  fMeanV0A[io] = 52.1674;  fMeanV0C[io] = 75.6410;  fMeanV0M[io] = 127.719;  io++;
   fRuns[io] = 272100;  fMeanV0A[io] = 52.1715;  fMeanV0C[io] = 75.6939;  fMeanV0M[io] = 127.778;  io++;
   fRuns[io] = 272101;  fMeanV0A[io] = 52.0067;  fMeanV0C[io] = 75.2956;  fMeanV0M[io] = 127.213;  io++;
   fRuns[io] = 272123;  fMeanV0A[io] = 51.5262;  fMeanV0C[io] = 75.0685;  fMeanV0M[io] = 126.494;  io++;
   fRuns[io] = 272151;  fMeanV0A[io] = 51.9838;  fMeanV0C[io] = 75.5625;  fMeanV0M[io] = 127.447;  io++;
   fRuns[io] = 272152;  fMeanV0A[io] = 51.8182;  fMeanV0C[io] = 75.2100;  fMeanV0M[io] = 126.950;  io++;
   fRuns[io] = 272153;  fMeanV0A[io] = 51.6823;  fMeanV0C[io] = 74.8613;  fMeanV0M[io] = 126.455;  io++;
   fRuns[io] = 272154;  fMeanV0A[io] = 51.4581;  fMeanV0C[io] = 74.6170;  fMeanV0M[io] = 125.979;  io++;
   fRuns[io] = 272155;  fMeanV0A[io] = 51.6622;  fMeanV0C[io] = 74.7943;  fMeanV0M[io] = 126.362;  io++;
   fRuns[io] = 272156;  fMeanV0A[io] = 51.3973;  fMeanV0C[io] = 74.6683;  fMeanV0M[io] = 125.984;  io++;
   fRuns[io] = 272194;  fMeanV0A[io] = 51.4691;  fMeanV0C[io] = 74.8235;  fMeanV0M[io] = 126.198;  io++;
   fRuns[io] = 272335;  fMeanV0A[io] = 51.5982;  fMeanV0C[io] = 74.9777;  fMeanV0M[io] = 126.482;  io++;
   fRuns[io] = 272340;  fMeanV0A[io] = 51.4677;  fMeanV0C[io] = 74.8457;  fMeanV0M[io] = 126.224;  io++;
   fRuns[io] = 272359;  fMeanV0A[io] = 51.3306;  fMeanV0C[io] = 74.6622;  fMeanV0M[io] = 125.909;  io++;
   fRuns[io] = 272360;  fMeanV0A[io] = 51.0299;  fMeanV0C[io] = 74.3888;  fMeanV0M[io] = 125.332;  io++;
   fRuns[io] = 272388;  fMeanV0A[io] = 51.1437;  fMeanV0C[io] = 74.5417;  fMeanV0M[io] = 125.617;  io++;
   fRuns[io] = 272389;  fMeanV0A[io] = 50.9623;  fMeanV0C[io] = 74.3924;  fMeanV0M[io] = 125.252;  io++;
   fRuns[io] = 272394;  fMeanV0A[io] = 50.9290;  fMeanV0C[io] = 74.2519;  fMeanV0M[io] = 125.088;  io++;
   fRuns[io] = 272395;  fMeanV0A[io] = 50.7202;  fMeanV0C[io] = 74.1708;  fMeanV0M[io] = 124.798;  io++;
   fRuns[io] = 272399;  fMeanV0A[io] = 50.7192;  fMeanV0C[io] = 74.4011;  fMeanV0M[io] = 125.030;  io++;
   fRuns[io] = 272400;  fMeanV0A[io] = 50.5527;  fMeanV0C[io] = 73.8907;  fMeanV0M[io] = 124.355;  io++;
   fRuns[io] = 272411;  fMeanV0A[io] = 50.4993;  fMeanV0C[io] = 73.6742;  fMeanV0M[io] = 124.082;  io++;
   fRuns[io] = 272413;  fMeanV0A[io] = 50.3774;  fMeanV0C[io] = 73.4175;  fMeanV0M[io] = 123.714;  io++;
   fRuns[io] = 272461;  fMeanV0A[io] = 50.7607;  fMeanV0C[io] = 73.8036;  fMeanV0M[io] = 124.480;  io++;
   fRuns[io] = 272462;  fMeanV0A[io] = 50.7460;  fMeanV0C[io] = 73.6974;  fMeanV0M[io] = 124.358;  io++;
   fRuns[io] = 272463;  fMeanV0A[io] = 50.4886;  fMeanV0C[io] = 73.5760;  fMeanV0M[io] = 123.968;  io++;
   fRuns[io] = 272466;  fMeanV0A[io] = 50.2826;  fMeanV0C[io] = 73.5741;  fMeanV0M[io] = 123.771;  io++;
   fRuns[io] = 272468;  fMeanV0A[io] = 50.1944;  fMeanV0C[io] = 73.2536;  fMeanV0M[io] = 123.354;  io++;
   fRuns[io] = 272521;  fMeanV0A[io] = 50.5546;  fMeanV0C[io] = 73.5881;  fMeanV0M[io] = 124.053;  io++;
   fRuns[io] = 272574;  fMeanV0A[io] = 50.2838;  fMeanV0C[io] = 73.9828;  fMeanV0M[io] = 124.179;  io++;
   fRuns[io] = 272575;  fMeanV0A[io] = 50.3029;  fMeanV0C[io] = 74.0504;  fMeanV0M[io] = 124.252;  io++;
   fRuns[io] = 272577;  fMeanV0A[io] = 50.1818;  fMeanV0C[io] = 73.9806;  fMeanV0M[io] = 124.069;  io++;
   fRuns[io] = 272585;  fMeanV0A[io] = 50.0835;  fMeanV0C[io] = 73.6951;  fMeanV0M[io] = 123.687;  io++;
   fRuns[io] = 272607;  fMeanV0A[io] = 50.4852;  fMeanV0C[io] = 73.9010;  fMeanV0M[io] = 124.283;  io++;
   fRuns[io] = 272608;  fMeanV0A[io] = 50.2609;  fMeanV0C[io] = 73.6934;  fMeanV0M[io] = 123.865;  io++;
   fRuns[io] = 272610;  fMeanV0A[io] = 50.1322;  fMeanV0C[io] = 73.6460;  fMeanV0M[io] = 123.681;  io++;
   fRuns[io] = 272620;  fMeanV0A[io] = 49.8823;  fMeanV0C[io] = 73.4033;  fMeanV0M[io] = 123.193;  io++;
   fRuns[io] = 272690;  fMeanV0A[io] = 50.0614;  fMeanV0C[io] = 73.3130;  fMeanV0M[io] = 123.316;  io++;
   fRuns[io] = 272691;  fMeanV0A[io] = 50.0595;  fMeanV0C[io] = 73.5772;  fMeanV0M[io] = 123.547;  io++;
   fRuns[io] = 272712;  fMeanV0A[io] = 50.0850;  fMeanV0C[io] = 73.7999;  fMeanV0M[io] = 123.812;  io++;
   fRuns[io] = 272747;  fMeanV0A[io] = 49.7171;  fMeanV0C[io] = 73.0238;  fMeanV0M[io] = 122.645;  io++;
   fRuns[io] = 272749;  fMeanV0A[io] = 49.5864;  fMeanV0C[io] = 72.9406;  fMeanV0M[io] = 122.440;  io++;
   fRuns[io] = 272760;  fMeanV0A[io] = 49.4270;  fMeanV0C[io] = 72.7850;  fMeanV0M[io] = 122.128;  io++;
   fRuns[io] = 272763;  fMeanV0A[io] = 49.5206;  fMeanV0C[io] = 72.7976;  fMeanV0M[io] = 122.212;  io++;
   fRuns[io] = 272764;  fMeanV0A[io] = 49.5446;  fMeanV0C[io] = 72.8155;  fMeanV0M[io] = 122.259;  io++;
   fRuns[io] = 272782;  fMeanV0A[io] = 49.8392;  fMeanV0C[io] = 73.0260;  fMeanV0M[io] = 122.775;  io++;
   fRuns[io] = 272783;  fMeanV0A[io] = 49.4300;  fMeanV0C[io] = 72.8302;  fMeanV0M[io] = 122.174;  io++;
   fRuns[io] = 272784;  fMeanV0A[io] = 49.2696;  fMeanV0C[io] = 72.5920;  fMeanV0M[io] = 121.774;  io++;
   fRuns[io] = 272828;  fMeanV0A[io] = 49.5802;  fMeanV0C[io] = 72.8905;  fMeanV0M[io] = 122.381;  io++;
   fRuns[io] = 272829;  fMeanV0A[io] = 49.4084;  fMeanV0C[io] = 72.7338;  fMeanV0M[io] = 122.049;  io++;
   fRuns[io] = 272833;  fMeanV0A[io] = 49.4586;  fMeanV0C[io] = 72.9376;  fMeanV0M[io] = 122.304;  io++;
   fRuns[io] = 272834;  fMeanV0A[io] = 49.3389;  fMeanV0C[io] = 72.6745;  fMeanV0M[io] = 121.922;  io++;
   fRuns[io] = 272836;  fMeanV0A[io] = 49.1600;  fMeanV0C[io] = 72.5357;  fMeanV0M[io] = 121.615;  io++;
   fRuns[io] = 272870;  fMeanV0A[io] = 49.3026;  fMeanV0C[io] = 72.7585;  fMeanV0M[io] = 121.969;  io++;
   fRuns[io] = 272871;  fMeanV0A[io] = 49.2714;  fMeanV0C[io] = 72.4963;  fMeanV0M[io] = 121.662;  io++;
   fRuns[io] = 272873;  fMeanV0A[io] = 49.1887;  fMeanV0C[io] = 72.4722;  fMeanV0M[io] = 121.564;  io++;
   fRuns[io] = 272880;  fMeanV0A[io] = 49.0320;  fMeanV0C[io] = 72.2156;  fMeanV0M[io] = 121.154;  io++;
   fRuns[io] = 272903;  fMeanV0A[io] = 49.1923;  fMeanV0C[io] = 72.5495;  fMeanV0M[io] = 121.642;  io++;
   fRuns[io] = 272905;  fMeanV0A[io] = 48.9982;  fMeanV0C[io] = 72.4414;  fMeanV0M[io] = 121.351;  io++;
   fRuns[io] = 272932;  fMeanV0A[io] = 49.2980;  fMeanV0C[io] = 72.7207;  fMeanV0M[io] = 121.914;  io++;
   fRuns[io] = 272933;  fMeanV0A[io] = 49.1996;  fMeanV0C[io] = 72.4224;  fMeanV0M[io] = 121.533;  io++;
   fRuns[io] = 272934;  fMeanV0A[io] = 49.1754;  fMeanV0C[io] = 72.2044;  fMeanV0M[io] = 121.304;  io++;
   fRuns[io] = 272935;  fMeanV0A[io] = 49.0940;  fMeanV0C[io] = 72.2362;  fMeanV0M[io] = 121.227;  io++;
   fRuns[io] = 272939;  fMeanV0A[io] = 49.0813;  fMeanV0C[io] = 72.2202;  fMeanV0M[io] = 121.213;  io++;
   fRuns[io] = 272947;  fMeanV0A[io] = 48.9342;  fMeanV0C[io] = 71.9362;  fMeanV0M[io] = 120.769;  io++;
   fRuns[io] = 272949;  fMeanV0A[io] = 48.8327;  fMeanV0C[io] = 72.3800;  fMeanV0M[io] = 121.113;  io++;
   fRuns[io] = 272976;  fMeanV0A[io] = 49.0288;  fMeanV0C[io] = 72.3273;  fMeanV0M[io] = 121.272;  io++;
   fRuns[io] = 272983;  fMeanV0A[io] = 49.0231;  fMeanV0C[io] = 72.3012;  fMeanV0M[io] = 121.245;  io++;
   fRuns[io] = 272985;  fMeanV0A[io] = 49.0029;  fMeanV0C[io] = 72.0154;  fMeanV0M[io] = 120.920;  io++;
   fRuns[io] = 273009;  fMeanV0A[io] = 48.8342;  fMeanV0C[io] = 71.9887;  fMeanV0M[io] = 120.729;  io++;
   fRuns[io] = 273010;  fMeanV0A[io] = 49.0456;  fMeanV0C[io] = 72.0603;  fMeanV0M[io] = 121.013;  io++;
   fRuns[io] = 273077;  fMeanV0A[io] = 49.1776;  fMeanV0C[io] = 72.5627;  fMeanV0M[io] = 121.644;  io++;
   fRuns[io] = 273099;  fMeanV0A[io] = 48.9950;  fMeanV0C[io] = 72.2096;  fMeanV0M[io] = 121.109;  io++;
   fRuns[io] = 273100;  fMeanV0A[io] = 48.8923;  fMeanV0C[io] = 72.2975;  fMeanV0M[io] = 121.099;  io++;
   fRuns[io] = 273103;  fMeanV0A[io] = 48.8253;  fMeanV0C[io] = 72.3443;  fMeanV0M[io] = 121.073;  io++;

   // FILTER_p-p_208_LHC17i
   fRuns[io] = 273591;  fMeanV0A[io] = 50.2660;  fMeanV0C[io] = 74.5786;  fMeanV0M[io] = 124.749;  io++;
   fRuns[io] = 273592;  fMeanV0A[io] = 50.2421;  fMeanV0C[io] = 74.5092;  fMeanV0M[io] = 124.681;  io++;
   fRuns[io] = 273593;  fMeanV0A[io] = 50.0925;  fMeanV0C[io] = 74.3288;  fMeanV0M[io] = 124.331;  io++;
   fRuns[io] = 273653;  fMeanV0A[io] = 50.4190;  fMeanV0C[io] = 74.7928;  fMeanV0M[io] = 125.119;  io++;
   fRuns[io] = 273654;  fMeanV0A[io] = 50.0285;  fMeanV0C[io] = 74.5648;  fMeanV0M[io] = 124.502;  io++;
   fRuns[io] = 273824;  fMeanV0A[io] = 50.0905;  fMeanV0C[io] = 74.3002;  fMeanV0M[io] = 124.318;  io++;
   fRuns[io] = 273825;  fMeanV0A[io] = 50.1727;  fMeanV0C[io] = 73.8634;  fMeanV0M[io] = 123.942;  io++;
   fRuns[io] = 273885;  fMeanV0A[io] = 50.4226;  fMeanV0C[io] = 74.0361;  fMeanV0M[io] = 124.349;  io++;
   fRuns[io] = 273886;  fMeanV0A[io] = 50.3457;  fMeanV0C[io] = 73.9031;  fMeanV0M[io] = 124.159;  io++;
   fRuns[io] = 273887;  fMeanV0A[io] = 50.1684;  fMeanV0C[io] = 73.7543;  fMeanV0M[io] = 123.830;  io++;
   fRuns[io] = 273889;  fMeanV0A[io] = 49.9876;  fMeanV0C[io] = 73.5779;  fMeanV0M[io] = 123.476;  io++;
   fRuns[io] = 273918;  fMeanV0A[io] = 50.1777;  fMeanV0C[io] = 73.6809;  fMeanV0M[io] = 123.769;  io++;
   fRuns[io] = 273942;  fMeanV0A[io] = 50.0613;  fMeanV0C[io] = 73.5322;  fMeanV0M[io] = 123.492;  io++;
   fRuns[io] = 273943;  fMeanV0A[io] = 50.0785;  fMeanV0C[io] = 73.4061;  fMeanV0M[io] = 123.387;  io++;
   fRuns[io] = 273946;  fMeanV0A[io] = 49.9929;  fMeanV0C[io] = 73.4604;  fMeanV0M[io] = 123.362;  io++;
   fRuns[io] = 273985;  fMeanV0A[io] = 49.8681;  fMeanV0C[io] = 73.3381;  fMeanV0M[io] = 123.114;  io++;
   fRuns[io] = 273986;  fMeanV0A[io] = 49.6755;  fMeanV0C[io] = 72.9331;  fMeanV0M[io] = 122.511;  io++;
   fRuns[io] = 274058;  fMeanV0A[io] = 49.8427;  fMeanV0C[io] = 73.3757;  fMeanV0M[io] = 123.131;  io++;
   fRuns[io] = 274092;  fMeanV0A[io] = 49.9666;  fMeanV0C[io] = 73.4847;  fMeanV0M[io] = 123.358;  io++;
   fRuns[io] = 274094;  fMeanV0A[io] = 49.6313;  fMeanV0C[io] = 73.0089;  fMeanV0M[io] = 122.541;  io++;
   fRuns[io] = 274125;  fMeanV0A[io] = 49.6027;  fMeanV0C[io] = 72.7830;  fMeanV0M[io] = 122.284;  io++;
   fRuns[io] = 274147;  fMeanV0A[io] = 49.5850;  fMeanV0C[io] = 72.8873;  fMeanV0M[io] = 122.383;  io++;
   fRuns[io] = 274148;  fMeanV0A[io] = 49.3790;  fMeanV0C[io] = 72.7186;  fMeanV0M[io] = 121.998;  io++;
   fRuns[io] = 274174;  fMeanV0A[io] = 49.4585;  fMeanV0C[io] = 73.1411;  fMeanV0M[io] = 122.508;  io++;
   fRuns[io] = 274212;  fMeanV0A[io] = 49.7309;  fMeanV0C[io] = 73.1030;  fMeanV0M[io] = 122.749;  io++;
   fRuns[io] = 274232;  fMeanV0A[io] = 49.6857;  fMeanV0C[io] = 73.0319;  fMeanV0M[io] = 122.616;  io++;
   fRuns[io] = 274258;  fMeanV0A[io] = 49.9177;  fMeanV0C[io] = 73.9083;  fMeanV0M[io] = 123.732;  io++;
   fRuns[io] = 274259;  fMeanV0A[io] = 49.4281;  fMeanV0C[io] = 73.0707;  fMeanV0M[io] = 122.383;  io++;
   fRuns[io] = 274263;  fMeanV0A[io] = 49.4057;  fMeanV0C[io] = 72.9003;  fMeanV0M[io] = 122.202;  io++;
   fRuns[io] = 274264;  fMeanV0A[io] = 49.3304;  fMeanV0C[io] = 72.8312;  fMeanV0M[io] = 122.060;  io++;
   fRuns[io] = 274266;  fMeanV0A[io] = 49.2764;  fMeanV0C[io] = 72.6248;  fMeanV0M[io] = 121.794;  io++;
   fRuns[io] = 274268;  fMeanV0A[io] = 49.3339;  fMeanV0C[io] = 72.7421;  fMeanV0M[io] = 121.979;  io++;
   fRuns[io] = 274269;  fMeanV0A[io] = 49.3169;  fMeanV0C[io] = 72.6382;  fMeanV0M[io] = 121.860;  io++;
   fRuns[io] = 274270;  fMeanV0A[io] = 49.2947;  fMeanV0C[io] = 72.2474;  fMeanV0M[io] = 121.439;  io++;
   fRuns[io] = 274271;  fMeanV0A[io] = 49.2076;  fMeanV0C[io] = 72.4259;  fMeanV0M[io] = 121.540;  io++;
   fRuns[io] = 274276;  fMeanV0A[io] = 49.2957;  fMeanV0C[io] = 72.3895;  fMeanV0M[io] = 121.592;  io++;
   fRuns[io] = 274278;  fMeanV0A[io] = 49.1525;  fMeanV0C[io] = 72.3075;  fMeanV0M[io] = 121.390;  io++;
   fRuns[io] = 274280;  fMeanV0A[io] = 49.1367;  fMeanV0C[io] = 72.1894;  fMeanV0M[io] = 121.237;  io++;
   fRuns[io] = 274281;  fMeanV0A[io] = 49.1475;  fMeanV0C[io] = 72.3028;  fMeanV0M[io] = 121.361;  io++;
   fRuns[io] = 274283;  fMeanV0A[io] = 49.0898;  fMeanV0C[io] = 72.0661;  fMeanV0M[io] = 121.070;  io++;
   fRuns[io] = 274329;  fMeanV0A[io] = 49.1445;  fMeanV0C[io] = 72.3343;  fMeanV0M[io] = 121.379;  io++;
   fRuns[io] = 274352;  fMeanV0A[io] = 49.1765;  fMeanV0C[io] = 72.4365;  fMeanV0M[io] = 121.508;  io++;
   fRuns[io] = 274360;  fMeanV0A[io] = 48.8969;  fMeanV0C[io] = 72.1477;  fMeanV0M[io] = 120.947;  io++;
   fRuns[io] = 274363;  fMeanV0A[io] = 48.8316;  fMeanV0C[io] = 72.1676;  fMeanV0M[io] = 120.906;  io++;
   fRuns[io] = 274364;  fMeanV0A[io] = 48.9224;  fMeanV0C[io] = 72.1595;  fMeanV0M[io] = 120.980;  io++;
   fRuns[io] = 274385;  fMeanV0A[io] = 49.4778;  fMeanV0C[io] = 72.8796;  fMeanV0M[io] = 122.256;  io++;
   fRuns[io] = 274386;  fMeanV0A[io] = 49.5394;  fMeanV0C[io] = 73.2377;  fMeanV0M[io] = 122.673;  io++;
   fRuns[io] = 274387;  fMeanV0A[io] = 49.1774;  fMeanV0C[io] = 72.6473;  fMeanV0M[io] = 121.738;  io++;
   fRuns[io] = 274388;  fMeanV0A[io] = 48.9727;  fMeanV0C[io] = 72.4297;  fMeanV0M[io] = 121.311;  io++;
   fRuns[io] = 274389;  fMeanV0A[io] = 48.8751;  fMeanV0C[io] = 72.3094;  fMeanV0M[io] = 121.094;  io++;
   fRuns[io] = 274390;  fMeanV0A[io] = 48.8300;  fMeanV0C[io] = 72.1632;  fMeanV0M[io] = 120.899;  io++;
   fRuns[io] = 274442;  fMeanV0A[io] = 49.2726;  fMeanV0C[io] = 72.5233;  fMeanV0M[io] = 121.704;  io++;

   // FILTER_p-p_208_LHC17j
   fRuns[io] = 274593;  fMeanV0A[io] = 49.6974;  fMeanV0C[io] = 73.2608;  fMeanV0M[io] = 122.881;  io++;
   fRuns[io] = 274594;  fMeanV0A[io] = 49.9968;  fMeanV0C[io] = 73.8848;  fMeanV0M[io] = 123.781;  io++;
   fRuns[io] = 274595;  fMeanV0A[io] = 49.3534;  fMeanV0C[io] = 73.0484;  fMeanV0M[io] = 122.308;  io++;
   fRuns[io] = 274596;  fMeanV0A[io] = 49.5712;  fMeanV0C[io] = 73.2036;  fMeanV0M[io] = 122.691;  io++;
   fRuns[io] = 274601;  fMeanV0A[io] = 49.4072;  fMeanV0C[io] = 72.9879;  fMeanV0M[io] = 122.306;  io++;
   fRuns[io] = 274653;  fMeanV0A[io] = 49.5447;  fMeanV0C[io] = 72.8710;  fMeanV0M[io] = 122.324;  io++;
   fRuns[io] = 274657;  fMeanV0A[io] = 49.3259;  fMeanV0C[io] = 72.6970;  fMeanV0M[io] = 121.931;  io++;
   fRuns[io] = 274667;  fMeanV0A[io] = 49.3602;  fMeanV0C[io] = 72.7367;  fMeanV0M[io] = 122.010;  io++;
   fRuns[io] = 274669;  fMeanV0A[io] = 49.2213;  fMeanV0C[io] = 72.6294;  fMeanV0M[io] = 121.756;  io++;
   fRuns[io] = 274671;  fMeanV0A[io] = 49.1885;  fMeanV0C[io] = 72.5030;  fMeanV0M[io] = 121.601;  io++;

   // FILTER_p-p_208_LHC17k
   fRuns[io] = 274690;  fMeanV0A[io] = 49.4361;  fMeanV0C[io] = 73.1383;  fMeanV0M[io] = 122.480;  io++;
   fRuns[io] = 274708;  fMeanV0A[io] = 49.1409;  fMeanV0C[io] = 72.9620;  fMeanV0M[io] = 122.012;  io++;
   fRuns[io] = 274801;  fMeanV0A[io] = 49.6845;  fMeanV0C[io] = 73.4095;  fMeanV0M[io] = 123.005;  io++;
   fRuns[io] = 274802;  fMeanV0A[io] = 49.6791;  fMeanV0C[io] = 73.1944;  fMeanV0M[io] = 122.791;  io++;
   fRuns[io] = 274803;  fMeanV0A[io] = 49.5945;  fMeanV0C[io] = 73.3128;  fMeanV0M[io] = 122.827;  io++;
   fRuns[io] = 274806;  fMeanV0A[io] = 49.4945;  fMeanV0C[io] = 73.1497;  fMeanV0M[io] = 122.548;  io++;
   fRuns[io] = 274815;  fMeanV0A[io] = 49.3836;  fMeanV0C[io] = 72.7867;  fMeanV0M[io] = 122.081;  io++;
   fRuns[io] = 274821;  fMeanV0A[io] = 49.0595;  fMeanV0C[io] = 72.8498;  fMeanV0M[io] = 121.807;  io++;
   fRuns[io] = 274822;  fMeanV0A[io] = 49.2628;  fMeanV0C[io] = 72.7364;  fMeanV0M[io] = 121.900;  io++;
   fRuns[io] = 274877;  fMeanV0A[io] = 49.3088;  fMeanV0C[io] = 72.8673;  fMeanV0M[io] = 122.088;  io++;
   fRuns[io] = 274878;  fMeanV0A[io] = 49.2687;  fMeanV0C[io] = 72.5762;  fMeanV0M[io] = 121.762;  io++;
   fRuns[io] = 274882;  fMeanV0A[io] = 49.1693;  fMeanV0C[io] = 72.4499;  fMeanV0M[io] = 121.526;  io++;
   fRuns[io] = 274886;  fMeanV0A[io] = 48.9891;  fMeanV0C[io] = 72.4220;  fMeanV0M[io] = 121.326;  io++;
   fRuns[io] = 274978;  fMeanV0A[io] = 48.7226;  fMeanV0C[io] = 71.9111;  fMeanV0M[io] = 120.540;  io++;
   fRuns[io] = 274979;  fMeanV0A[io] = 48.5456;  fMeanV0C[io] = 71.6292;  fMeanV0M[io] = 120.082;  io++;
   fRuns[io] = 275067;  fMeanV0A[io] = 48.8254;  fMeanV0C[io] = 71.8224;  fMeanV0M[io] = 120.560;  io++;
   fRuns[io] = 275068;  fMeanV0A[io] = 48.7343;  fMeanV0C[io] = 71.9820;  fMeanV0M[io] = 120.623;  io++;
   fRuns[io] = 275073;  fMeanV0A[io] = 48.5070;  fMeanV0C[io] = 71.7341;  fMeanV0M[io] = 120.147;  io++;
   fRuns[io] = 275075;  fMeanV0A[io] = 48.5411;  fMeanV0C[io] = 71.7163;  fMeanV0M[io] = 120.153;  io++;
   fRuns[io] = 275076;  fMeanV0A[io] = 48.4893;  fMeanV0C[io] = 71.5680;  fMeanV0M[io] = 119.972;  io++;
   fRuns[io] = 275149;  fMeanV0A[io] = 48.5669;  fMeanV0C[io] = 71.9396;  fMeanV0M[io] = 120.401;  io++;
   fRuns[io] = 275150;  fMeanV0A[io] = 48.6348;  fMeanV0C[io] = 71.7380;  fMeanV0M[io] = 120.285;  io++;
   fRuns[io] = 275151;  fMeanV0A[io] = 48.4067;  fMeanV0C[io] = 71.4647;  fMeanV0M[io] = 119.781;  io++;
   fRuns[io] = 275173;  fMeanV0A[io] = 48.5900;  fMeanV0C[io] = 71.6601;  fMeanV0M[io] = 120.159;  io++;
   fRuns[io] = 275174;  fMeanV0A[io] = 48.3523;  fMeanV0C[io] = 71.2877;  fMeanV0M[io] = 119.549;  io++;
   fRuns[io] = 275177;  fMeanV0A[io] = 48.3906;  fMeanV0C[io] = 71.2395;  fMeanV0M[io] = 119.528;  io++;
   fRuns[io] = 275180;  fMeanV0A[io] = 48.2788;  fMeanV0C[io] = 71.2679;  fMeanV0M[io] = 119.464;  io++;
   fRuns[io] = 275184;  fMeanV0A[io] = 48.3695;  fMeanV0C[io] = 71.3580;  fMeanV0M[io] = 119.638;  io++;
   fRuns[io] = 275188;  fMeanV0A[io] = 48.2679;  fMeanV0C[io] = 71.1312;  fMeanV0M[io] = 119.309;  io++;
   fRuns[io] = 275239;  fMeanV0A[io] = 48.6627;  fMeanV0C[io] = 71.9006;  fMeanV0M[io] = 120.463;  io++;
   fRuns[io] = 275245;  fMeanV0A[io] = 48.4415;  fMeanV0C[io] = 71.6733;  fMeanV0M[io] = 120.021;  io++;
   fRuns[io] = 275246;  fMeanV0A[io] = 48.2894;  fMeanV0C[io] = 71.2098;  fMeanV0M[io] = 119.397;  io++;
   fRuns[io] = 275247;  fMeanV0A[io] = 48.3141;  fMeanV0C[io] = 71.1526;  fMeanV0M[io] = 119.371;  io++;
   fRuns[io] = 275283;  fMeanV0A[io] = 48.3930;  fMeanV0C[io] = 71.3861;  fMeanV0M[io] = 119.675;  io++;
   fRuns[io] = 275314;  fMeanV0A[io] = 48.3818;  fMeanV0C[io] = 71.4387;  fMeanV0M[io] = 119.737;  io++;
   fRuns[io] = 275322;  fMeanV0A[io] = 48.1681;  fMeanV0C[io] = 71.2355;  fMeanV0M[io] = 119.308;  io++;
   fRuns[io] = 275324;  fMeanV0A[io] = 48.2073;  fMeanV0C[io] = 71.0885;  fMeanV0M[io] = 119.204;  io++;
   fRuns[io] = 275326;  fMeanV0A[io] = 47.9992;  fMeanV0C[io] = 70.9966;  fMeanV0M[io] = 118.901;  io++;
   fRuns[io] = 275328;  fMeanV0A[io] = 48.0596;  fMeanV0C[io] = 71.0316;  fMeanV0M[io] = 118.997;  io++;
   fRuns[io] = 275332;  fMeanV0A[io] = 48.0221;  fMeanV0C[io] = 70.7904;  fMeanV0M[io] = 118.721;  io++;
   fRuns[io] = 275333;  fMeanV0A[io] = 47.8933;  fMeanV0C[io] = 70.7745;  fMeanV0M[io] = 118.576;  io++;
   fRuns[io] = 275360;  fMeanV0A[io] = 47.9891;  fMeanV0C[io] = 70.9015;  fMeanV0M[io] = 118.793;  io++;
   fRuns[io] = 275361;  fMeanV0A[io] = 47.8616;  fMeanV0C[io] = 70.6260;  fMeanV0M[io] = 118.387;  io++;
   fRuns[io] = 275369;  fMeanV0A[io] = 47.8081;  fMeanV0C[io] = 70.7618;  fMeanV0M[io] = 118.465;  io++;
   fRuns[io] = 275372;  fMeanV0A[io] = 47.7468;  fMeanV0C[io] = 70.5111;  fMeanV0M[io] = 118.162;  io++;
   fRuns[io] = 275401;  fMeanV0A[io] = 47.7460;  fMeanV0C[io] = 70.6796;  fMeanV0M[io] = 118.335;  io++;
   fRuns[io] = 275404;  fMeanV0A[io] = 47.7016;  fMeanV0C[io] = 70.5196;  fMeanV0M[io] = 118.119;  io++;
   fRuns[io] = 275406;  fMeanV0A[io] = 47.6938;  fMeanV0C[io] = 70.3862;  fMeanV0M[io] = 117.981;  io++;
   fRuns[io] = 275443;  fMeanV0A[io] = 47.8046;  fMeanV0C[io] = 70.6607;  fMeanV0M[io] = 118.370;  io++;
   fRuns[io] = 275448;  fMeanV0A[io] = 47.7077;  fMeanV0C[io] = 70.4228;  fMeanV0M[io] = 118.032;  io++;
   fRuns[io] = 275452;  fMeanV0A[io] = 47.6239;  fMeanV0C[io] = 70.1450;  fMeanV0M[io] = 117.683;  io++;
   fRuns[io] = 275453;  fMeanV0A[io] = 47.6725;  fMeanV0C[io] = 70.1925;  fMeanV0M[io] = 117.775;  io++;
   fRuns[io] = 275456;  fMeanV0A[io] = 47.6044;  fMeanV0C[io] = 70.4087;  fMeanV0M[io] = 117.920;  io++;
   fRuns[io] = 275457;  fMeanV0A[io] = 47.6590;  fMeanV0C[io] = 70.1171;  fMeanV0M[io] = 117.682;  io++;
   fRuns[io] = 275459;  fMeanV0A[io] = 47.5447;  fMeanV0C[io] = 70.1885;  fMeanV0M[io] = 117.637;  io++;
   fRuns[io] = 275467;  fMeanV0A[io] = 47.6700;  fMeanV0C[io] = 70.1641;  fMeanV0M[io] = 117.747;  io++;
   fRuns[io] = 275471;  fMeanV0A[io] = 47.7455;  fMeanV0C[io] = 70.2907;  fMeanV0M[io] = 117.926;  io++;
   fRuns[io] = 275472;  fMeanV0A[io] = 47.5815;  fMeanV0C[io] = 70.3797;  fMeanV0M[io] = 117.891;  io++;
   fRuns[io] = 275515;  fMeanV0A[io] = 47.8499;  fMeanV0C[io] = 70.5666;  fMeanV0M[io] = 118.312;  io++;
   fRuns[io] = 275558;  fMeanV0A[io] = 47.7911;  fMeanV0C[io] = 70.5661;  fMeanV0M[io] = 118.264;  io++;
   fRuns[io] = 275559;  fMeanV0A[io] = 47.6404;  fMeanV0C[io] = 70.3156;  fMeanV0M[io] = 117.883;  io++;
   fRuns[io] = 275612;  fMeanV0A[io] = 47.5896;  fMeanV0C[io] = 70.3783;  fMeanV0M[io] = 117.860;  io++;
   fRuns[io] = 275617;  fMeanV0A[io] = 47.5266;  fMeanV0C[io] = 70.4028;  fMeanV0M[io] = 117.830;  io++;
   fRuns[io] = 275621;  fMeanV0A[io] = 47.3731;  fMeanV0C[io] = 70.2453;  fMeanV0M[io] = 117.519;  io++;
   fRuns[io] = 275622;  fMeanV0A[io] = 47.5025;  fMeanV0C[io] = 70.2305;  fMeanV0M[io] = 117.634;  io++;
   fRuns[io] = 275623;  fMeanV0A[io] = 47.4678;  fMeanV0C[io] = 70.0714;  fMeanV0M[io] = 117.424;  io++;
   fRuns[io] = 275624;  fMeanV0A[io] = 47.5096;  fMeanV0C[io] = 70.2897;  fMeanV0M[io] = 117.695;  io++;
   fRuns[io] = 275647;  fMeanV0A[io] = 47.5702;  fMeanV0C[io] = 70.3507;  fMeanV0M[io] = 117.824;  io++;
   fRuns[io] = 275648;  fMeanV0A[io] = 47.5005;  fMeanV0C[io] = 70.2415;  fMeanV0M[io] = 117.640;  io++;
   fRuns[io] = 275650;  fMeanV0A[io] = 47.4471;  fMeanV0C[io] = 70.1572;  fMeanV0M[io] = 117.513;  io++;
   fRuns[io] = 275661;  fMeanV0A[io] = 47.3917;  fMeanV0C[io] = 70.0408;  fMeanV0M[io] = 117.322;  io++;
   fRuns[io] = 275664;  fMeanV0A[io] = 47.3517;  fMeanV0C[io] = 70.1474;  fMeanV0M[io] = 117.408;  io++;
   fRuns[io] = 275847;  fMeanV0A[io] = 48.1717;  fMeanV0C[io] = 71.4433;  fMeanV0M[io] = 119.524;  io++;
   fRuns[io] = 276097;  fMeanV0A[io] = 48.5274;  fMeanV0C[io] = 71.2214;  fMeanV0M[io] = 119.668;  io++;
   fRuns[io] = 276098;  fMeanV0A[io] = 48.1574;  fMeanV0C[io] = 71.0525;  fMeanV0M[io] = 119.113;  io++;
   fRuns[io] = 276099;  fMeanV0A[io] = 47.9949;  fMeanV0C[io] = 71.0837;  fMeanV0M[io] = 118.981;  io++;
   fRuns[io] = 276102;  fMeanV0A[io] = 48.1152;  fMeanV0C[io] = 70.9106;  fMeanV0M[io] = 118.904;  io++;
   fRuns[io] = 276104;  fMeanV0A[io] = 48.1623;  fMeanV0C[io] = 71.0480;  fMeanV0M[io] = 119.106;  io++;
   fRuns[io] = 276135;  fMeanV0A[io] = 48.1831;  fMeanV0C[io] = 71.3942;  fMeanV0M[io] = 119.473;  io++;
   fRuns[io] = 276140;  fMeanV0A[io] = 48.1473;  fMeanV0C[io] = 71.1490;  fMeanV0M[io] = 119.208;  io++;
   fRuns[io] = 276145;  fMeanV0A[io] = 47.5328;  fMeanV0C[io] = 70.1184;  fMeanV0M[io] = 117.585;  io++;
   fRuns[io] = 276166;  fMeanV0A[io] = 48.3912;  fMeanV0C[io] = 71.1991;  fMeanV0M[io] = 119.471;  io++;
   fRuns[io] = 276169;  fMeanV0A[io] = 48.2920;  fMeanV0C[io] = 70.9927;  fMeanV0M[io] = 119.182;  io++;
   fRuns[io] = 276170;  fMeanV0A[io] = 47.9869;  fMeanV0C[io] = 70.8249;  fMeanV0M[io] = 118.717;  io++;
   fRuns[io] = 276177;  fMeanV0A[io] = 47.8960;  fMeanV0C[io] = 70.5601;  fMeanV0M[io] = 118.367;  io++;
   fRuns[io] = 276178;  fMeanV0A[io] = 48.0769;  fMeanV0C[io] = 70.9982;  fMeanV0M[io] = 118.973;  io++;
   fRuns[io] = 276205;  fMeanV0A[io] = 48.0656;  fMeanV0C[io] = 70.9920;  fMeanV0M[io] = 118.964;  io++;
   fRuns[io] = 276230;  fMeanV0A[io] = 48.1770;  fMeanV0C[io] = 71.0347;  fMeanV0M[io] = 119.104;  io++;
   fRuns[io] = 276257;  fMeanV0A[io] = 48.0643;  fMeanV0C[io] = 71.2597;  fMeanV0M[io] = 119.233;  io++;
   fRuns[io] = 276259;  fMeanV0A[io] = 47.8995;  fMeanV0C[io] = 70.5906;  fMeanV0M[io] = 118.394;  io++;
   fRuns[io] = 276290;  fMeanV0A[io] = 48.1672;  fMeanV0C[io] = 70.9525;  fMeanV0M[io] = 119.017;  io++;
   fRuns[io] = 276292;  fMeanV0A[io] = 47.8769;  fMeanV0C[io] = 70.6476;  fMeanV0M[io] = 118.433;  io++;
   fRuns[io] = 276294;  fMeanV0A[io] = 47.9562;  fMeanV0C[io] = 70.8424;  fMeanV0M[io] = 118.710;  io++;
   fRuns[io] = 276297;  fMeanV0A[io] = 47.7896;  fMeanV0C[io] = 70.5442;  fMeanV0M[io] = 118.244;  io++;
   fRuns[io] = 276302;  fMeanV0A[io] = 47.6560;  fMeanV0C[io] = 70.6031;  fMeanV0M[io] = 118.172;  io++;
   fRuns[io] = 276348;  fMeanV0A[io] = 48.0721;  fMeanV0C[io] = 70.8778;  fMeanV0M[io] = 118.875;  io++;
   fRuns[io] = 276351;  fMeanV0A[io] = 47.8088;  fMeanV0C[io] = 70.4937;  fMeanV0M[io] = 118.217;  io++;
   fRuns[io] = 276435;  fMeanV0A[io] = 47.7487;  fMeanV0C[io] = 70.6149;  fMeanV0M[io] = 118.265;  io++;
   fRuns[io] = 276437;  fMeanV0A[io] = 47.7408;  fMeanV0C[io] = 70.5759;  fMeanV0M[io] = 118.235;  io++;
   fRuns[io] = 276438;  fMeanV0A[io] = 47.7326;  fMeanV0C[io] = 70.5045;  fMeanV0M[io] = 118.148;  io++;
   fRuns[io] = 276439;  fMeanV0A[io] = 47.7145;  fMeanV0C[io] = 70.5274;  fMeanV0M[io] = 118.160;  io++;
   fRuns[io] = 276462;  fMeanV0A[io] = 47.8612;  fMeanV0C[io] = 70.6572;  fMeanV0M[io] = 118.424;  io++;
   fRuns[io] = 276506;  fMeanV0A[io] = 47.7181;  fMeanV0C[io] = 70.4317;  fMeanV0M[io] = 118.050;  io++;
   fRuns[io] = 276507;  fMeanV0A[io] = 47.5349;  fMeanV0C[io] = 70.1805;  fMeanV0M[io] = 117.643;  io++;
   fRuns[io] = 276508;  fMeanV0A[io] = 47.7335;  fMeanV0C[io] = 70.4014;  fMeanV0M[io] = 118.041;  io++;

   // FILTER_p-p_208_LHC17l
   fRuns[io] = 276551;  fMeanV0A[io] = 47.6952;  fMeanV0C[io] = 70.4245;  fMeanV0M[io] = 118.023;  io++;
   fRuns[io] = 276552;  fMeanV0A[io] = 47.6101;  fMeanV0C[io] = 70.2790;  fMeanV0M[io] = 117.792;  io++;
   fRuns[io] = 276553;  fMeanV0A[io] = 47.6440;  fMeanV0C[io] = 70.2442;  fMeanV0M[io] = 117.793;  io++;
   fRuns[io] = 276556;  fMeanV0A[io] = 47.5785;  fMeanV0C[io] = 69.8654;  fMeanV0M[io] = 117.342;  io++;
   fRuns[io] = 276557;  fMeanV0A[io] = 47.4947;  fMeanV0C[io] = 70.1151;  fMeanV0M[io] = 117.517;  io++;
   fRuns[io] = 276608;  fMeanV0A[io] = 47.7653;  fMeanV0C[io] = 70.4996;  fMeanV0M[io] = 118.166;  io++;
   fRuns[io] = 276644;  fMeanV0A[io] = 47.6481;  fMeanV0C[io] = 70.4390;  fMeanV0M[io] = 117.988;  io++;
   fRuns[io] = 276670;  fMeanV0A[io] = 47.9609;  fMeanV0C[io] = 71.0836;  fMeanV0M[io] = 118.962;  io++;
   fRuns[io] = 276671;  fMeanV0A[io] = 47.8982;  fMeanV0C[io] = 70.5279;  fMeanV0M[io] = 118.321;  io++;
   fRuns[io] = 276672;  fMeanV0A[io] = 47.7587;  fMeanV0C[io] = 70.4432;  fMeanV0M[io] = 118.081;  io++;
   fRuns[io] = 276674;  fMeanV0A[io] = 47.7458;  fMeanV0C[io] = 70.3684;  fMeanV0M[io] = 118.037;  io++;
   fRuns[io] = 276675;  fMeanV0A[io] = 47.7089;  fMeanV0C[io] = 70.3985;  fMeanV0M[io] = 118.014;  io++;
   fRuns[io] = 276762;  fMeanV0A[io] = 47.8361;  fMeanV0C[io] = 70.3933;  fMeanV0M[io] = 118.126;  io++;
   fRuns[io] = 276916;  fMeanV0A[io] = 48.0012;  fMeanV0C[io] = 70.8739;  fMeanV0M[io] = 118.748;  io++;
   fRuns[io] = 276917;  fMeanV0A[io] = 47.7800;  fMeanV0C[io] = 70.5511;  fMeanV0M[io] = 118.236;  io++;
   fRuns[io] = 276920;  fMeanV0A[io] = 47.7374;  fMeanV0C[io] = 70.5359;  fMeanV0M[io] = 118.180;  io++;
   fRuns[io] = 276967;  fMeanV0A[io] = 48.2902;  fMeanV0C[io] = 71.6620;  fMeanV0M[io] = 119.851;  io++;
   fRuns[io] = 276969;  fMeanV0A[io] = 48.2605;  fMeanV0C[io] = 70.8025;  fMeanV0M[io] = 118.992;  io++;
   fRuns[io] = 276970;  fMeanV0A[io] = 47.7806;  fMeanV0C[io] = 70.7702;  fMeanV0M[io] = 118.446;  io++;
   fRuns[io] = 276971;  fMeanV0A[io] = 47.8429;  fMeanV0C[io] = 70.7176;  fMeanV0M[io] = 118.476;  io++;
   fRuns[io] = 276972;  fMeanV0A[io] = 47.9182;  fMeanV0C[io] = 70.9821;  fMeanV0M[io] = 118.782;  io++;
   fRuns[io] = 277015;  fMeanV0A[io] = 47.8780;  fMeanV0C[io] = 70.7378;  fMeanV0M[io] = 118.527;  io++;
   fRuns[io] = 277016;  fMeanV0A[io] = 47.9278;  fMeanV0C[io] = 70.7856;  fMeanV0M[io] = 118.620;  io++;
   fRuns[io] = 277017;  fMeanV0A[io] = 47.8422;  fMeanV0C[io] = 70.5587;  fMeanV0M[io] = 118.300;  io++;
   fRuns[io] = 277037;  fMeanV0A[io] = 48.0390;  fMeanV0C[io] = 70.6852;  fMeanV0M[io] = 118.628;  io++;
   fRuns[io] = 277073;  fMeanV0A[io] = 48.1501;  fMeanV0C[io] = 70.9917;  fMeanV0M[io] = 119.052;  io++;
   fRuns[io] = 277076;  fMeanV0A[io] = 47.8177;  fMeanV0C[io] = 70.5298;  fMeanV0M[io] = 118.256;  io++;
   fRuns[io] = 277079;  fMeanV0A[io] = 47.8114;  fMeanV0C[io] = 70.3020;  fMeanV0M[io] = 118.017;  io++;
   fRuns[io] = 277082;  fMeanV0A[io] = 47.6292;  fMeanV0C[io] = 70.3810;  fMeanV0M[io] = 117.901;  io++;
   fRuns[io] = 277087;  fMeanV0A[io] = 47.6896;  fMeanV0C[io] = 70.2716;  fMeanV0M[io] = 117.859;  io++;
   fRuns[io] = 277091;  fMeanV0A[io] = 47.7243;  fMeanV0C[io] = 70.1714;  fMeanV0M[io] = 117.789;  io++;
   fRuns[io] = 277117;  fMeanV0A[io] = 47.7280;  fMeanV0C[io] = 70.2837;  fMeanV0M[io] = 117.915;  io++;
   fRuns[io] = 277121;  fMeanV0A[io] = 47.6319;  fMeanV0C[io] = 70.0789;  fMeanV0M[io] = 117.611;  io++;
   fRuns[io] = 277155;  fMeanV0A[io] = 47.9335;  fMeanV0C[io] = 70.3470;  fMeanV0M[io] = 118.206;  io++;
   fRuns[io] = 277180;  fMeanV0A[io] = 47.6803;  fMeanV0C[io] = 70.3170;  fMeanV0M[io] = 117.898;  io++;
   fRuns[io] = 277181;  fMeanV0A[io] = 47.7635;  fMeanV0C[io] = 70.1852;  fMeanV0M[io] = 117.818;  io++;
   fRuns[io] = 277182;  fMeanV0A[io] = 47.7665;  fMeanV0C[io] = 70.1862;  fMeanV0M[io] = 117.859;  io++;
   fRuns[io] = 277183;  fMeanV0A[io] = 47.5358;  fMeanV0C[io] = 69.9863;  fMeanV0M[io] = 117.403;  io++;
   fRuns[io] = 277184;  fMeanV0A[io] = 47.5723;  fMeanV0C[io] = 70.2074;  fMeanV0M[io] = 117.662;  io++;
   fRuns[io] = 277188;  fMeanV0A[io] = 47.5489;  fMeanV0C[io] = 70.0487;  fMeanV0M[io] = 117.513;  io++;
   fRuns[io] = 277189;  fMeanV0A[io] = 47.4521;  fMeanV0C[io] = 69.9192;  fMeanV0M[io] = 117.290;  io++;
   fRuns[io] = 277193;  fMeanV0A[io] = 47.6457;  fMeanV0C[io] = 70.1871;  fMeanV0M[io] = 117.730;  io++;
   fRuns[io] = 277194;  fMeanV0A[io] = 47.4832;  fMeanV0C[io] = 69.7494;  fMeanV0M[io] = 117.141;  io++;
   fRuns[io] = 277196;  fMeanV0A[io] = 47.4906;  fMeanV0C[io] = 70.1156;  fMeanV0M[io] = 117.502;  io++;
   fRuns[io] = 277197;  fMeanV0A[io] = 47.3949;  fMeanV0C[io] = 69.8423;  fMeanV0M[io] = 117.136;  io++;
   fRuns[io] = 277256;  fMeanV0A[io] = 47.4566;  fMeanV0C[io] = 70.2373;  fMeanV0M[io] = 117.605;  io++;
   fRuns[io] = 277257;  fMeanV0A[io] = 47.5216;  fMeanV0C[io] = 70.2457;  fMeanV0M[io] = 117.673;  io++;
   fRuns[io] = 277262;  fMeanV0A[io] = 47.5190;  fMeanV0C[io] = 70.2244;  fMeanV0M[io] = 117.638;  io++;
   fRuns[io] = 277293;  fMeanV0A[io] = 48.0182;  fMeanV0C[io] = 70.7495;  fMeanV0M[io] = 118.667;  io++;
   fRuns[io] = 277310;  fMeanV0A[io] = 47.6425;  fMeanV0C[io] = 70.3100;  fMeanV0M[io] = 117.861;  io++;
   fRuns[io] = 277312;  fMeanV0A[io] = 47.5128;  fMeanV0C[io] = 70.1000;  fMeanV0M[io] = 117.514;  io++;
   fRuns[io] = 277314;  fMeanV0A[io] = 47.5175;  fMeanV0C[io] = 70.2110;  fMeanV0M[io] = 117.633;  io++;
   fRuns[io] = 277360;  fMeanV0A[io] = 47.5793;  fMeanV0C[io] = 70.2940;  fMeanV0M[io] = 117.780;  io++;
   fRuns[io] = 277383;  fMeanV0A[io] = 47.9087;  fMeanV0C[io] = 70.8410;  fMeanV0M[io] = 118.667;  io++;
   fRuns[io] = 277384;  fMeanV0A[io] = 47.5741;  fMeanV0C[io] = 70.0829;  fMeanV0M[io] = 117.564;  io++;
   fRuns[io] = 277385;  fMeanV0A[io] = 47.5265;  fMeanV0C[io] = 70.3253;  fMeanV0M[io] = 117.802;  io++;
   fRuns[io] = 277386;  fMeanV0A[io] = 47.6024;  fMeanV0C[io] = 70.4510;  fMeanV0M[io] = 117.978;  io++;
   fRuns[io] = 277389;  fMeanV0A[io] = 47.4024;  fMeanV0C[io] = 69.8592;  fMeanV0M[io] = 117.172;  io++;
   fRuns[io] = 277416;  fMeanV0A[io] = 47.6157;  fMeanV0C[io] = 70.1602;  fMeanV0M[io] = 117.680;  io++;
   fRuns[io] = 277417;  fMeanV0A[io] = 47.4776;  fMeanV0C[io] = 70.1743;  fMeanV0M[io] = 117.550;  io++;
   fRuns[io] = 277418;  fMeanV0A[io] = 47.3854;  fMeanV0C[io] = 70.2455;  fMeanV0M[io] = 117.525;  io++;
   fRuns[io] = 277472;  fMeanV0A[io] = 47.3709;  fMeanV0C[io] = 70.0799;  fMeanV0M[io] = 117.320;  io++;
   fRuns[io] = 277473;  fMeanV0A[io] = 47.3714;  fMeanV0C[io] = 70.0113;  fMeanV0M[io] = 117.293;  io++;
   fRuns[io] = 277476;  fMeanV0A[io] = 47.5416;  fMeanV0C[io] = 70.1655;  fMeanV0M[io] = 117.601;  io++;
   fRuns[io] = 277477;  fMeanV0A[io] = 47.5455;  fMeanV0C[io] = 70.1511;  fMeanV0M[io] = 117.602;  io++;
   fRuns[io] = 277478;  fMeanV0A[io] = 47.4865;  fMeanV0C[io] = 70.1775;  fMeanV0M[io] = 117.560;  io++;
   fRuns[io] = 277479;  fMeanV0A[io] = 47.2454;  fMeanV0C[io] = 69.7820;  fMeanV0M[io] = 116.927;  io++;
   fRuns[io] = 277530;  fMeanV0A[io] = 47.4426;  fMeanV0C[io] = 70.0379;  fMeanV0M[io] = 117.394;  io++;
   fRuns[io] = 277531;  fMeanV0A[io] = 47.4862;  fMeanV0C[io] = 69.8532;  fMeanV0M[io] = 117.257;  io++;
   fRuns[io] = 277534;  fMeanV0A[io] = 47.5127;  fMeanV0C[io] = 69.9139;  fMeanV0M[io] = 117.320;  io++;
   fRuns[io] = 277536;  fMeanV0A[io] = 47.2971;  fMeanV0C[io] = 70.1481;  fMeanV0M[io] = 117.358;  io++;
   fRuns[io] = 277537;  fMeanV0A[io] = 47.4272;  fMeanV0C[io] = 69.9273;  fMeanV0M[io] = 117.256;  io++;
   fRuns[io] = 277574;  fMeanV0A[io] = 47.6726;  fMeanV0C[io] = 70.2307;  fMeanV0M[io] = 117.801;  io++;
   fRuns[io] = 277575;  fMeanV0A[io] = 47.2681;  fMeanV0C[io] = 69.8407;  fMeanV0M[io] = 117.007;  io++;
   fRuns[io] = 277576;  fMeanV0A[io] = 47.1368;  fMeanV0C[io] = 69.7677;  fMeanV0M[io] = 116.804;  io++;
   fRuns[io] = 277577;  fMeanV0A[io] = 47.1377;  fMeanV0C[io] = 69.7195;  fMeanV0M[io] = 116.770;  io++;
   fRuns[io] = 277721;  fMeanV0A[io] = 47.4002;  fMeanV0C[io] = 69.5548;  fMeanV0M[io] = 116.859;  io++;
   fRuns[io] = 277722;  fMeanV0A[io] = 47.4025;  fMeanV0C[io] = 70.0235;  fMeanV0M[io] = 117.355;  io++;
   fRuns[io] = 277723;  fMeanV0A[io] = 47.2115;  fMeanV0C[io] = 69.6248;  fMeanV0M[io] = 116.749;  io++;
   fRuns[io] = 277725;  fMeanV0A[io] = 47.1834;  fMeanV0C[io] = 69.6720;  fMeanV0M[io] = 116.754;  io++;
   fRuns[io] = 277745;  fMeanV0A[io] = 47.3104;  fMeanV0C[io] = 69.5937;  fMeanV0M[io] = 116.803;  io++;
   fRuns[io] = 277746;  fMeanV0A[io] = 47.1592;  fMeanV0C[io] = 69.6473;  fMeanV0M[io] = 116.698;  io++;
   fRuns[io] = 277747;  fMeanV0A[io] = 47.0830;  fMeanV0C[io] = 69.4083;  fMeanV0M[io] = 116.405;  io++;
   fRuns[io] = 277749;  fMeanV0A[io] = 47.2599;  fMeanV0C[io] = 69.5771;  fMeanV0M[io] = 116.745;  io++;
   fRuns[io] = 277794;  fMeanV0A[io] = 47.4563;  fMeanV0C[io] = 70.0191;  fMeanV0M[io] = 117.368;  io++;
   fRuns[io] = 277795;  fMeanV0A[io] = 47.2753;  fMeanV0C[io] = 69.7699;  fMeanV0M[io] = 116.935;  io++;
   fRuns[io] = 277799;  fMeanV0A[io] = 47.0350;  fMeanV0C[io] = 69.6128;  fMeanV0M[io] = 116.548;  io++;
   fRuns[io] = 277800;  fMeanV0A[io] = 47.1712;  fMeanV0C[io] = 69.7291;  fMeanV0M[io] = 116.806;  io++;
   fRuns[io] = 277801;  fMeanV0A[io] = 47.1671;  fMeanV0C[io] = 69.8360;  fMeanV0M[io] = 116.904;  io++;
   fRuns[io] = 277802;  fMeanV0A[io] = 46.8883;  fMeanV0C[io] = 69.4555;  fMeanV0M[io] = 116.237;  io++;
   fRuns[io] = 277805;  fMeanV0A[io] = 47.2048;  fMeanV0C[io] = 69.6007;  fMeanV0M[io] = 116.696;  io++;
   fRuns[io] = 277834;  fMeanV0A[io] = 47.2185;  fMeanV0C[io] = 69.6999;  fMeanV0M[io] = 116.825;  io++;
   fRuns[io] = 277836;  fMeanV0A[io] = 47.3483;  fMeanV0C[io] = 69.2906;  fMeanV0M[io] = 116.577;  io++;
   fRuns[io] = 277841;  fMeanV0A[io] = 47.0067;  fMeanV0C[io] = 69.4061;  fMeanV0M[io] = 116.293;  io++;
   fRuns[io] = 277842;  fMeanV0A[io] = 47.1481;  fMeanV0C[io] = 69.4535;  fMeanV0M[io] = 116.498;  io++;
   fRuns[io] = 277845;  fMeanV0A[io] = 47.0533;  fMeanV0C[io] = 69.5833;  fMeanV0M[io] = 116.547;  io++;
   fRuns[io] = 277847;  fMeanV0A[io] = 47.2418;  fMeanV0C[io] = 69.5710;  fMeanV0M[io] = 116.698;  io++;
   fRuns[io] = 277848;  fMeanV0A[io] = 46.8961;  fMeanV0C[io] = 69.2096;  fMeanV0M[io] = 116.001;  io++;
   fRuns[io] = 277870;  fMeanV0A[io] = 47.2211;  fMeanV0C[io] = 69.7525;  fMeanV0M[io] = 116.867;  io++;
   fRuns[io] = 277876;  fMeanV0A[io] = 47.1956;  fMeanV0C[io] = 69.5829;  fMeanV0M[io] = 116.673;  io++;
   fRuns[io] = 277897;  fMeanV0A[io] = 47.1185;  fMeanV0C[io] = 69.5328;  fMeanV0M[io] = 116.556;  io++;
   fRuns[io] = 277898;  fMeanV0A[io] = 47.2076;  fMeanV0C[io] = 69.4210;  fMeanV0M[io] = 116.533;  io++;
   fRuns[io] = 277899;  fMeanV0A[io] = 47.1173;  fMeanV0C[io] = 69.4302;  fMeanV0M[io] = 116.462;  io++;
   fRuns[io] = 277900;  fMeanV0A[io] = 46.9724;  fMeanV0C[io] = 69.4127;  fMeanV0M[io] = 116.283;  io++;
   fRuns[io] = 277903;  fMeanV0A[io] = 46.9277;  fMeanV0C[io] = 69.1975;  fMeanV0M[io] = 116.025;  io++;
   fRuns[io] = 277904;  fMeanV0A[io] = 47.0724;  fMeanV0C[io] = 69.2245;  fMeanV0M[io] = 116.180;  io++;
   fRuns[io] = 277907;  fMeanV0A[io] = 46.9989;  fMeanV0C[io] = 69.0393;  fMeanV0M[io] = 115.948;  io++;
   fRuns[io] = 277930;  fMeanV0A[io] = 47.4844;  fMeanV0C[io] = 70.1489;  fMeanV0M[io] = 117.527;  io++;
   fRuns[io] = 277952;  fMeanV0A[io] = 47.4747;  fMeanV0C[io] = 69.4238;  fMeanV0M[io] = 116.815;  io++;
   fRuns[io] = 277987;  fMeanV0A[io] = 47.3402;  fMeanV0C[io] = 69.6010;  fMeanV0M[io] = 116.842;  io++;
   fRuns[io] = 277989;  fMeanV0A[io] = 46.9519;  fMeanV0C[io] = 69.6914;  fMeanV0M[io] = 116.545;  io++;
   fRuns[io] = 277991;  fMeanV0A[io] = 47.1540;  fMeanV0C[io] = 69.5543;  fMeanV0M[io] = 116.628;  io++;
   fRuns[io] = 277996;  fMeanV0A[io] = 47.0729;  fMeanV0C[io] = 69.4706;  fMeanV0M[io] = 116.455;  io++;
   fRuns[io] = 278121;  fMeanV0A[io] = 47.0660;  fMeanV0C[io] = 69.3718;  fMeanV0M[io] = 116.331;  io++;
   fRuns[io] = 278122;  fMeanV0A[io] = 47.0096;  fMeanV0C[io] = 69.1590;  fMeanV0M[io] = 116.064;  io++;
   fRuns[io] = 278123;  fMeanV0A[io] = 46.7641;  fMeanV0C[io] = 69.3338;  fMeanV0M[io] = 115.992;  io++;
   fRuns[io] = 278126;  fMeanV0A[io] = 46.9326;  fMeanV0C[io] = 69.1175;  fMeanV0M[io] = 115.954;  io++;
   fRuns[io] = 278127;  fMeanV0A[io] = 46.9536;  fMeanV0C[io] = 69.1496;  fMeanV0M[io] = 115.996;  io++;
   fRuns[io] = 278158;  fMeanV0A[io] = 47.5851;  fMeanV0C[io] = 70.0317;  fMeanV0M[io] = 117.531;  io++;
   fRuns[io] = 278164;  fMeanV0A[io] = 46.9199;  fMeanV0C[io] = 68.9832;  fMeanV0M[io] = 115.811;  io++;
   fRuns[io] = 278165;  fMeanV0A[io] = 46.8767;  fMeanV0C[io] = 69.1714;  fMeanV0M[io] = 115.947;  io++;
   fRuns[io] = 278166;  fMeanV0A[io] = 46.9689;  fMeanV0C[io] = 69.2249;  fMeanV0M[io] = 116.096;  io++;
   fRuns[io] = 278167;  fMeanV0A[io] = 46.8298;  fMeanV0C[io] = 68.9090;  fMeanV0M[io] = 115.644;  io++;
   fRuns[io] = 278189;  fMeanV0A[io] = 47.0030;  fMeanV0C[io] = 69.1518;  fMeanV0M[io] = 116.059;  io++;
   fRuns[io] = 278191;  fMeanV0A[io] = 46.8210;  fMeanV0C[io] = 68.7961;  fMeanV0M[io] = 115.520;  io++;
   fRuns[io] = 278215;  fMeanV0A[io] = 46.9981;  fMeanV0C[io] = 69.1845;  fMeanV0M[io] = 116.090;  io++;
   fRuns[io] = 278216;  fMeanV0A[io] = 46.8432;  fMeanV0C[io] = 68.8909;  fMeanV0M[io] = 115.637;  io++;

   // FILTER_p-p_208_LHC17m
   fRuns[io] = 278914;  fMeanV0A[io] = 47.5584;  fMeanV0C[io] = 70.8999;  fMeanV0M[io] = 118.366;  io++;
   fRuns[io] = 278915;  fMeanV0A[io] = 47.3598;  fMeanV0C[io] = 70.6420;  fMeanV0M[io] = 117.910;  io++;
   fRuns[io] = 278936;  fMeanV0A[io] = 47.2613;  fMeanV0C[io] = 70.2696;  fMeanV0M[io] = 117.436;  io++;
   fRuns[io] = 278939;  fMeanV0A[io] = 47.1330;  fMeanV0C[io] = 69.9903;  fMeanV0M[io] = 117.035;  io++;
   fRuns[io] = 278941;  fMeanV0A[io] = 47.1021;  fMeanV0C[io] = 69.9612;  fMeanV0M[io] = 116.967;  io++;
   fRuns[io] = 278959;  fMeanV0A[io] = 47.4857;  fMeanV0C[io] = 69.9176;  fMeanV0M[io] = 117.309;  io++;
   fRuns[io] = 278960;  fMeanV0A[io] = 47.4229;  fMeanV0C[io] = 70.1337;  fMeanV0M[io] = 117.467;  io++;
   fRuns[io] = 278963;  fMeanV0A[io] = 47.3773;  fMeanV0C[io] = 70.0271;  fMeanV0M[io] = 117.296;  io++;
   fRuns[io] = 278964;  fMeanV0A[io] = 47.3311;  fMeanV0C[io] = 69.9214;  fMeanV0M[io] = 117.157;  io++;
   fRuns[io] = 278999;  fMeanV0A[io] = 47.5001;  fMeanV0C[io] = 70.1385;  fMeanV0M[io] = 117.528;  io++;
   fRuns[io] = 279000;  fMeanV0A[io] = 47.6205;  fMeanV0C[io] = 70.2633;  fMeanV0M[io] = 117.790;  io++;
   fRuns[io] = 279005;  fMeanV0A[io] = 47.6344;  fMeanV0C[io] = 70.0891;  fMeanV0M[io] = 117.629;  io++;
   fRuns[io] = 279007;  fMeanV0A[io] = 47.5517;  fMeanV0C[io] = 69.9478;  fMeanV0M[io] = 117.409;  io++;
   fRuns[io] = 279008;  fMeanV0A[io] = 47.5525;  fMeanV0C[io] = 69.9198;  fMeanV0M[io] = 117.374;  io++;
   fRuns[io] = 279035;  fMeanV0A[io] = 47.8736;  fMeanV0C[io] = 70.0757;  fMeanV0M[io] = 117.832;  io++;
   fRuns[io] = 279036;  fMeanV0A[io] = 47.6760;  fMeanV0C[io] = 70.0893;  fMeanV0M[io] = 117.663;  io++;
   fRuns[io] = 279041;  fMeanV0A[io] = 47.3668;  fMeanV0C[io] = 69.9352;  fMeanV0M[io] = 117.207;  io++;
   fRuns[io] = 279043;  fMeanV0A[io] = 47.4847;  fMeanV0C[io] = 69.8768;  fMeanV0M[io] = 117.268;  io++;
   fRuns[io] = 279044;  fMeanV0A[io] = 47.2765;  fMeanV0C[io] = 69.6056;  fMeanV0M[io] = 116.790;  io++;
   fRuns[io] = 279068;  fMeanV0A[io] = 47.7310;  fMeanV0C[io] = 69.8889;  fMeanV0M[io] = 117.538;  io++;
   fRuns[io] = 279069;  fMeanV0A[io] = 47.3268;  fMeanV0C[io] = 69.5654;  fMeanV0M[io] = 116.802;  io++;
   fRuns[io] = 279073;  fMeanV0A[io] = 47.2348;  fMeanV0C[io] = 69.3284;  fMeanV0M[io] = 116.466;  io++;
   fRuns[io] = 279074;  fMeanV0A[io] = 47.3272;  fMeanV0C[io] = 69.5974;  fMeanV0M[io] = 116.830;  io++;
   fRuns[io] = 279075;  fMeanV0A[io] = 47.2860;  fMeanV0C[io] = 69.3461;  fMeanV0M[io] = 116.530;  io++;
   fRuns[io] = 279106;  fMeanV0A[io] = 47.4136;  fMeanV0C[io] = 69.6230;  fMeanV0M[io] = 116.955;  io++;
   fRuns[io] = 279107;  fMeanV0A[io] = 47.3107;  fMeanV0C[io] = 69.5521;  fMeanV0M[io] = 116.770;  io++;
   fRuns[io] = 279117;  fMeanV0A[io] = 47.0102;  fMeanV0C[io] = 69.3016;  fMeanV0M[io] = 116.195;  io++;
   fRuns[io] = 279118;  fMeanV0A[io] = 47.2331;  fMeanV0C[io] = 69.4484;  fMeanV0M[io] = 116.597;  io++;
   fRuns[io] = 279122;  fMeanV0A[io] = 47.1588;  fMeanV0C[io] = 69.2411;  fMeanV0M[io] = 116.302;  io++;
   fRuns[io] = 279123;  fMeanV0A[io] = 47.1671;  fMeanV0C[io] = 69.3351;  fMeanV0M[io] = 116.408;  io++;
   fRuns[io] = 279130;  fMeanV0A[io] = 47.1527;  fMeanV0C[io] = 69.2403;  fMeanV0M[io] = 116.292;  io++;
   fRuns[io] = 279155;  fMeanV0A[io] = 47.1385;  fMeanV0C[io] = 69.3710;  fMeanV0M[io] = 116.413;  io++;
   fRuns[io] = 279157;  fMeanV0A[io] = 46.9764;  fMeanV0C[io] = 69.1988;  fMeanV0M[io] = 116.081;  io++;
   fRuns[io] = 279199;  fMeanV0A[io] = 47.1607;  fMeanV0C[io] = 69.1800;  fMeanV0M[io] = 116.250;  io++;
   fRuns[io] = 279201;  fMeanV0A[io] = 46.8602;  fMeanV0C[io] = 69.0006;  fMeanV0M[io] = 115.760;  io++;
   fRuns[io] = 279207;  fMeanV0A[io] = 47.0177;  fMeanV0C[io] = 69.0963;  fMeanV0M[io] = 116.016;  io++;
   fRuns[io] = 279208;  fMeanV0A[io] = 47.1172;  fMeanV0C[io] = 69.1357;  fMeanV0M[io] = 116.143;  io++;
   fRuns[io] = 279232;  fMeanV0A[io] = 47.2713;  fMeanV0C[io] = 69.8804;  fMeanV0M[io] = 117.046;  io++;
   fRuns[io] = 279234;  fMeanV0A[io] = 46.8622;  fMeanV0C[io] = 68.9877;  fMeanV0M[io] = 115.752;  io++;
   fRuns[io] = 279235;  fMeanV0A[io] = 46.9300;  fMeanV0C[io] = 69.0008;  fMeanV0M[io] = 115.835;  io++;
   fRuns[io] = 279238;  fMeanV0A[io] = 46.8384;  fMeanV0C[io] = 68.8022;  fMeanV0M[io] = 115.550;  io++;
   fRuns[io] = 279242;  fMeanV0A[io] = 46.9701;  fMeanV0C[io] = 68.8797;  fMeanV0M[io] = 115.759;  io++;
   fRuns[io] = 279264;  fMeanV0A[io] = 46.8995;  fMeanV0C[io] = 69.0007;  fMeanV0M[io] = 115.809;  io++;
   fRuns[io] = 279265;  fMeanV0A[io] = 46.9632;  fMeanV0C[io] = 69.0384;  fMeanV0M[io] = 115.910;  io++;
   fRuns[io] = 279267;  fMeanV0A[io] = 46.8088;  fMeanV0C[io] = 68.8255;  fMeanV0M[io] = 115.532;  io++;
   fRuns[io] = 279268;  fMeanV0A[io] = 46.6536;  fMeanV0C[io] = 68.7255;  fMeanV0M[io] = 115.292;  io++;
   fRuns[io] = 279270;  fMeanV0A[io] = 46.6660;  fMeanV0C[io] = 68.8739;  fMeanV0M[io] = 115.439;  io++;
   fRuns[io] = 279273;  fMeanV0A[io] = 46.6351;  fMeanV0C[io] = 68.5582;  fMeanV0M[io] = 115.089;  io++;
   fRuns[io] = 279274;  fMeanV0A[io] = 46.6412;  fMeanV0C[io] = 68.5081;  fMeanV0M[io] = 115.047;  io++;
   fRuns[io] = 279309;  fMeanV0A[io] = 46.8474;  fMeanV0C[io] = 68.9172;  fMeanV0M[io] = 115.661;  io++;
   fRuns[io] = 279310;  fMeanV0A[io] = 46.6545;  fMeanV0C[io] = 68.7816;  fMeanV0M[io] = 115.341;  io++;
   fRuns[io] = 279312;  fMeanV0A[io] = 46.6150;  fMeanV0C[io] = 68.5463;  fMeanV0M[io] = 115.057;  io++;
   fRuns[io] = 279342;  fMeanV0A[io] = 46.6159;  fMeanV0C[io] = 68.6330;  fMeanV0M[io] = 115.161;  io++;
   fRuns[io] = 279344;  fMeanV0A[io] = 46.5496;  fMeanV0C[io] = 68.6228;  fMeanV0M[io] = 115.065;  io++;
   fRuns[io] = 279348;  fMeanV0A[io] = 46.6181;  fMeanV0C[io] = 68.5695;  fMeanV0M[io] = 115.086;  io++;
   fRuns[io] = 279349;  fMeanV0A[io] = 46.5086;  fMeanV0C[io] = 68.3710;  fMeanV0M[io] = 114.781;  io++;
   fRuns[io] = 279354;  fMeanV0A[io] = 46.6075;  fMeanV0C[io] = 68.5430;  fMeanV0M[io] = 115.048;  io++;
   fRuns[io] = 279355;  fMeanV0A[io] = 46.6085;  fMeanV0C[io] = 68.4336;  fMeanV0M[io] = 114.955;  io++;
   fRuns[io] = 279391;  fMeanV0A[io] = 46.6563;  fMeanV0C[io] = 68.5193;  fMeanV0M[io] = 115.087;  io++;
   fRuns[io] = 279410;  fMeanV0A[io] = 46.6432;  fMeanV0C[io] = 68.5582;  fMeanV0M[io] = 115.107;  io++;
   fRuns[io] = 279435;  fMeanV0A[io] = 46.8159;  fMeanV0C[io] = 68.5669;  fMeanV0M[io] = 115.286;  io++;
   fRuns[io] = 279439;  fMeanV0A[io] = 46.6805;  fMeanV0C[io] = 68.7043;  fMeanV0M[io] = 115.301;  io++;
   fRuns[io] = 279441;  fMeanV0A[io] = 46.6945;  fMeanV0C[io] = 68.6137;  fMeanV0M[io] = 115.201;  io++;
   fRuns[io] = 279483;  fMeanV0A[io] = 46.5595;  fMeanV0C[io] = 68.5172;  fMeanV0M[io] = 114.967;  io++;
   fRuns[io] = 279487;  fMeanV0A[io] = 46.4140;  fMeanV0C[io] = 68.4282;  fMeanV0M[io] = 114.743;  io++;
   fRuns[io] = 279488;  fMeanV0A[io] = 46.4798;  fMeanV0C[io] = 68.4841;  fMeanV0M[io] = 114.872;  io++;
   fRuns[io] = 279491;  fMeanV0A[io] = 46.4193;  fMeanV0C[io] = 68.6263;  fMeanV0M[io] = 114.941;  io++;
   fRuns[io] = 279550;  fMeanV0A[io] = 46.5750;  fMeanV0C[io] = 68.6536;  fMeanV0M[io] = 115.129;  io++;
   fRuns[io] = 279559;  fMeanV0A[io] = 46.4585;  fMeanV0C[io] = 68.4557;  fMeanV0M[io] = 114.812;  io++;
   fRuns[io] = 279630;  fMeanV0A[io] = 46.0290;  fMeanV0C[io] = 67.6344;  fMeanV0M[io] = 113.556;  io++;
   fRuns[io] = 279632;  fMeanV0A[io] = 45.9176;  fMeanV0C[io] = 67.4828;  fMeanV0M[io] = 113.307;  io++;
   fRuns[io] = 279641;  fMeanV0A[io] = 45.8750;  fMeanV0C[io] = 67.5333;  fMeanV0M[io] = 113.312;  io++;
   fRuns[io] = 279642;  fMeanV0A[io] = 45.7174;  fMeanV0C[io] = 67.3256;  fMeanV0M[io] = 112.929;  io++;
   fRuns[io] = 279676;  fMeanV0A[io] = 46.0945;  fMeanV0C[io] = 67.7100;  fMeanV0M[io] = 113.713;  io++;
   fRuns[io] = 279677;  fMeanV0A[io] = 46.1691;  fMeanV0C[io] = 67.6554;  fMeanV0M[io] = 113.724;  io++;
   fRuns[io] = 279679;  fMeanV0A[io] = 45.9076;  fMeanV0C[io] = 67.4638;  fMeanV0M[io] = 113.283;  io++;
   fRuns[io] = 279682;  fMeanV0A[io] = 45.9569;  fMeanV0C[io] = 67.5730;  fMeanV0M[io] = 113.432;  io++;
   fRuns[io] = 279683;  fMeanV0A[io] = 45.9086;  fMeanV0C[io] = 67.3293;  fMeanV0M[io] = 113.140;  io++;
   fRuns[io] = 279684;  fMeanV0A[io] = 45.8835;  fMeanV0C[io] = 67.4842;  fMeanV0M[io] = 113.266;  io++;
   fRuns[io] = 279687;  fMeanV0A[io] = 46.5700;  fMeanV0C[io] = 67.8653;  fMeanV0M[io] = 114.328;  io++;
   fRuns[io] = 279688;  fMeanV0A[io] = 45.6938;  fMeanV0C[io] = 67.4284;  fMeanV0M[io] = 113.021;  io++;
   fRuns[io] = 279689;  fMeanV0A[io] = 45.8310;  fMeanV0C[io] = 67.5651;  fMeanV0M[io] = 113.286;  io++;
   fRuns[io] = 279715;  fMeanV0A[io] = 46.1017;  fMeanV0C[io] = 67.6245;  fMeanV0M[io] = 113.629;  io++;
   fRuns[io] = 279718;  fMeanV0A[io] = 45.8911;  fMeanV0C[io] = 67.5453;  fMeanV0M[io] = 113.335;  io++;
   fRuns[io] = 279719;  fMeanV0A[io] = 45.8329;  fMeanV0C[io] = 67.5760;  fMeanV0M[io] = 113.287;  io++;
   fRuns[io] = 279747;  fMeanV0A[io] = 46.1922;  fMeanV0C[io] = 67.7411;  fMeanV0M[io] = 113.834;  io++;
   fRuns[io] = 279749;  fMeanV0A[io] = 45.9653;  fMeanV0C[io] = 67.5984;  fMeanV0M[io] = 113.467;  io++;
   fRuns[io] = 279773;  fMeanV0A[io] = 45.8898;  fMeanV0C[io] = 67.5382;  fMeanV0M[io] = 113.323;  io++;
   fRuns[io] = 279826;  fMeanV0A[io] = 46.0809;  fMeanV0C[io] = 67.7898;  fMeanV0M[io] = 113.765;  io++;
   fRuns[io] = 279827;  fMeanV0A[io] = 45.9058;  fMeanV0C[io] = 67.4861;  fMeanV0M[io] = 113.281;  io++;
   fRuns[io] = 279830;  fMeanV0A[io] = 45.9507;  fMeanV0C[io] = 67.6602;  fMeanV0M[io] = 113.509;  io++;
   fRuns[io] = 279853;  fMeanV0A[io] = 46.0248;  fMeanV0C[io] = 67.5692;  fMeanV0M[io] = 113.493;  io++;
   fRuns[io] = 279854;  fMeanV0A[io] = 45.8736;  fMeanV0C[io] = 67.4600;  fMeanV0M[io] = 113.235;  io++;
   fRuns[io] = 279855;  fMeanV0A[io] = 46.0342;  fMeanV0C[io] = 67.8271;  fMeanV0M[io] = 113.766;  io++;
   fRuns[io] = 279879;  fMeanV0A[io] = 46.0531;  fMeanV0C[io] = 67.7138;  fMeanV0M[io] = 113.656;  io++;
   fRuns[io] = 280051;  fMeanV0A[io] = 45.5491;  fMeanV0C[io] = 66.8105;  fMeanV0M[io] = 112.260;  io++;
   fRuns[io] = 280052;  fMeanV0A[io] = 45.4991;  fMeanV0C[io] = 66.7990;  fMeanV0M[io] = 112.191;  io++;
   fRuns[io] = 280066;  fMeanV0A[io] = 45.3800;  fMeanV0C[io] = 66.7512;  fMeanV0M[io] = 112.024;  io++;
   fRuns[io] = 280107;  fMeanV0A[io] = 45.4867;  fMeanV0C[io] = 66.9044;  fMeanV0M[io] = 112.287;  io++;
   fRuns[io] = 280108;  fMeanV0A[io] = 45.4719;  fMeanV0C[io] = 66.9426;  fMeanV0M[io] = 112.313;  io++;
   fRuns[io] = 280111;  fMeanV0A[io] = 45.5357;  fMeanV0C[io] = 66.9535;  fMeanV0M[io] = 112.378;  io++;
   fRuns[io] = 280114;  fMeanV0A[io] = 45.5325;  fMeanV0C[io] = 66.9442;  fMeanV0M[io] = 112.377;  io++;
   fRuns[io] = 280118;  fMeanV0A[io] = 45.4180;  fMeanV0C[io] = 67.0443;  fMeanV0M[io] = 112.337;  io++;
   fRuns[io] = 280126;  fMeanV0A[io] = 45.4677;  fMeanV0C[io] = 66.8080;  fMeanV0M[io] = 112.180;  io++;
   fRuns[io] = 280131;  fMeanV0A[io] = 45.3852;  fMeanV0C[io] = 66.5905;  fMeanV0M[io] = 111.882;  io++;
   fRuns[io] = 280134;  fMeanV0A[io] = 45.4478;  fMeanV0C[io] = 66.9062;  fMeanV0M[io] = 112.261;  io++;
   fRuns[io] = 280135;  fMeanV0A[io] = 45.5855;  fMeanV0C[io] = 66.8675;  fMeanV0M[io] = 112.343;  io++;
   fRuns[io] = 280140;  fMeanV0A[io] = 45.5346;  fMeanV0C[io] = 66.8598;  fMeanV0M[io] = 112.292;  io++;

   // FILTER_p-p_208_LHC17o
   fRuns[io] = 280282;  fMeanV0A[io] = 45.5355;  fMeanV0C[io] = 67.1971;  fMeanV0M[io] = 112.610;  io++;
   fRuns[io] = 280284;  fMeanV0A[io] = 45.3364;  fMeanV0C[io] = 67.4261;  fMeanV0M[io] = 112.635;  io++;
   fRuns[io] = 280285;  fMeanV0A[io] = 45.2837;  fMeanV0C[io] = 66.9026;  fMeanV0M[io] = 112.091;  io++;
   fRuns[io] = 280286;  fMeanV0A[io] = 45.4859;  fMeanV0C[io] = 67.1130;  fMeanV0M[io] = 112.489;  io++;
   fRuns[io] = 280290;  fMeanV0A[io] = 45.4136;  fMeanV0C[io] = 67.0327;  fMeanV0M[io] = 112.334;  io++;
   fRuns[io] = 280310;  fMeanV0A[io] = 45.5346;  fMeanV0C[io] = 67.0857;  fMeanV0M[io] = 112.532;  io++;
   fRuns[io] = 280312;  fMeanV0A[io] = 45.6416;  fMeanV0C[io] = 67.1391;  fMeanV0M[io] = 112.672;  io++;
   fRuns[io] = 280348;  fMeanV0A[io] = 45.8808;  fMeanV0C[io] = 67.5609;  fMeanV0M[io] = 113.345;  io++;
   fRuns[io] = 280349;  fMeanV0A[io] = 45.6028;  fMeanV0C[io] = 67.1107;  fMeanV0M[io] = 112.618;  io++;
   fRuns[io] = 280350;  fMeanV0A[io] = 45.4319;  fMeanV0C[io] = 66.9553;  fMeanV0M[io] = 112.283;  io++;
   fRuns[io] = 280351;  fMeanV0A[io] = 45.5463;  fMeanV0C[io] = 67.0348;  fMeanV0M[io] = 112.487;  io++;
   fRuns[io] = 280374;  fMeanV0A[io] = 45.5042;  fMeanV0C[io] = 67.0734;  fMeanV0M[io] = 112.475;  io++;
   fRuns[io] = 280375;  fMeanV0A[io] = 45.3379;  fMeanV0C[io] = 66.8896;  fMeanV0M[io] = 112.110;  io++;
   fRuns[io] = 280403;  fMeanV0A[io] = 45.5050;  fMeanV0C[io] = 67.1524;  fMeanV0M[io] = 112.554;  io++;
   fRuns[io] = 280405;  fMeanV0A[io] = 45.4915;  fMeanV0C[io] = 67.1147;  fMeanV0M[io] = 112.485;  io++;
   fRuns[io] = 280406;  fMeanV0A[io] = 45.2458;  fMeanV0C[io] = 67.1826;  fMeanV0M[io] = 112.325;  io++;
   fRuns[io] = 280412;  fMeanV0A[io] = 45.5299;  fMeanV0C[io] = 67.1690;  fMeanV0M[io] = 112.593;  io++;
   fRuns[io] = 280415;  fMeanV0A[io] = 45.3890;  fMeanV0C[io] = 66.8723;  fMeanV0M[io] = 112.161;  io++;
   fRuns[io] = 280419;  fMeanV0A[io] = 45.3816;  fMeanV0C[io] = 67.0878;  fMeanV0M[io] = 112.350;  io++;
   fRuns[io] = 280443;  fMeanV0A[io] = 45.6940;  fMeanV0C[io] = 66.9450;  fMeanV0M[io] = 112.537;  io++;
   fRuns[io] = 280445;  fMeanV0A[io] = 45.3995;  fMeanV0C[io] = 66.8869;  fMeanV0M[io] = 112.176;  io++;
   fRuns[io] = 280446;  fMeanV0A[io] = 45.3620;  fMeanV0C[io] = 66.9236;  fMeanV0M[io] = 112.185;  io++;
   fRuns[io] = 280447;  fMeanV0A[io] = 45.2946;  fMeanV0C[io] = 67.0842;  fMeanV0M[io] = 112.274;  io++;
   fRuns[io] = 280448;  fMeanV0A[io] = 45.5020;  fMeanV0C[io] = 67.3603;  fMeanV0M[io] = 112.741;  io++;
   fRuns[io] = 280490;  fMeanV0A[io] = 45.5690;  fMeanV0C[io] = 67.1112;  fMeanV0M[io] = 112.571;  io++;
   fRuns[io] = 280499;  fMeanV0A[io] = 45.2806;  fMeanV0C[io] = 66.9991;  fMeanV0M[io] = 112.173;  io++;
   fRuns[io] = 280518;  fMeanV0A[io] = 45.3653;  fMeanV0C[io] = 66.9941;  fMeanV0M[io] = 112.255;  io++;
   fRuns[io] = 280519;  fMeanV0A[io] = 45.3012;  fMeanV0C[io] = 66.9911;  fMeanV0M[io] = 112.174;  io++;
   fRuns[io] = 280546;  fMeanV0A[io] = 45.5463;  fMeanV0C[io] = 67.2099;  fMeanV0M[io] = 112.652;  io++;
   fRuns[io] = 280547;  fMeanV0A[io] = 45.3663;  fMeanV0C[io] = 66.9790;  fMeanV0M[io] = 112.247;  io++;
   fRuns[io] = 280550;  fMeanV0A[io] = 45.2444;  fMeanV0C[io] = 66.8500;  fMeanV0M[io] = 111.991;  io++;
   fRuns[io] = 280551;  fMeanV0A[io] = 45.1582;  fMeanV0C[io] = 66.5637;  fMeanV0M[io] = 111.612;  io++;
   fRuns[io] = 280574;  fMeanV0A[io] = 45.5773;  fMeanV0C[io] = 67.0881;  fMeanV0M[io] = 112.563;  io++;
   fRuns[io] = 280581;  fMeanV0A[io] = 45.2802;  fMeanV0C[io] = 66.8854;  fMeanV0M[io] = 112.053;  io++;
   fRuns[io] = 280583;  fMeanV0A[io] = 45.1489;  fMeanV0C[io] = 66.8991;  fMeanV0M[io] = 111.944;  io++;
   fRuns[io] = 280613;  fMeanV0A[io] = 45.4793;  fMeanV0C[io] = 66.8521;  fMeanV0M[io] = 112.220;  io++;
   fRuns[io] = 280634;  fMeanV0A[io] = 45.3269;  fMeanV0C[io] = 66.9508;  fMeanV0M[io] = 112.164;  io++;
   fRuns[io] = 280636;  fMeanV0A[io] = 45.2487;  fMeanV0C[io] = 66.4151;  fMeanV0M[io] = 111.561;  io++;
   fRuns[io] = 280637;  fMeanV0A[io] = 45.1821;  fMeanV0C[io] = 66.6735;  fMeanV0M[io] = 111.745;  io++;
   fRuns[io] = 280639;  fMeanV0A[io] = 45.2449;  fMeanV0C[io] = 66.8700;  fMeanV0M[io] = 112.015;  io++;
   fRuns[io] = 280645;  fMeanV0A[io] = 45.1542;  fMeanV0C[io] = 66.8181;  fMeanV0M[io] = 111.872;  io++;
   fRuns[io] = 280647;  fMeanV0A[io] = 45.2772;  fMeanV0C[io] = 66.4556;  fMeanV0M[io] = 111.658;  io++;
   fRuns[io] = 280671;  fMeanV0A[io] = 45.3446;  fMeanV0C[io] = 66.8434;  fMeanV0M[io] = 112.084;  io++;
   fRuns[io] = 280679;  fMeanV0A[io] = 45.1949;  fMeanV0C[io] = 66.7462;  fMeanV0M[io] = 111.834;  io++;
   fRuns[io] = 280681;  fMeanV0A[io] = 45.2247;  fMeanV0C[io] = 66.8675;  fMeanV0M[io] = 111.984;  io++;
   fRuns[io] = 280705;  fMeanV0A[io] = 45.2603;  fMeanV0C[io] = 66.7130;  fMeanV0M[io] = 111.866;  io++;
   fRuns[io] = 280706;  fMeanV0A[io] = 45.2252;  fMeanV0C[io] = 66.5995;  fMeanV0M[io] = 111.716;  io++;
   fRuns[io] = 280729;  fMeanV0A[io] = 45.2763;  fMeanV0C[io] = 66.8903;  fMeanV0M[io] = 112.068;  io++;
   fRuns[io] = 280753;  fMeanV0A[io] = 45.1461;  fMeanV0C[io] = 66.9397;  fMeanV0M[io] = 111.975;  io++;
   fRuns[io] = 280754;  fMeanV0A[io] = 45.2890;  fMeanV0C[io] = 67.0691;  fMeanV0M[io] = 112.260;  io++;
   fRuns[io] = 280755;  fMeanV0A[io] = 45.1680;  fMeanV0C[io] = 66.8859;  fMeanV0M[io] = 111.949;  io++;
   fRuns[io] = 280756;  fMeanV0A[io] = 45.0770;  fMeanV0C[io] = 66.8247;  fMeanV0M[io] = 111.796;  io++;
   fRuns[io] = 280757;  fMeanV0A[io] = 45.1874;  fMeanV0C[io] = 66.8034;  fMeanV0M[io] = 111.878;  io++;
   fRuns[io] = 280761;  fMeanV0A[io] = 45.1414;  fMeanV0C[io] = 66.8556;  fMeanV0M[io] = 111.890;  io++;
   fRuns[io] = 280762;  fMeanV0A[io] = 44.9090;  fMeanV0C[io] = 66.8842;  fMeanV0M[io] = 111.697;  io++;
   fRuns[io] = 280763;  fMeanV0A[io] = 45.2536;  fMeanV0C[io] = 66.7356;  fMeanV0M[io] = 111.895;  io++;
   fRuns[io] = 280764;  fMeanV0A[io] = 45.2283;  fMeanV0C[io] = 66.8701;  fMeanV0M[io] = 111.992;  io++;
   fRuns[io] = 280765;  fMeanV0A[io] = 45.1304;  fMeanV0C[io] = 66.8019;  fMeanV0M[io] = 111.843;  io++;
   fRuns[io] = 280766;  fMeanV0A[io] = 45.0961;  fMeanV0C[io] = 66.6870;  fMeanV0M[io] = 111.676;  io++;
   fRuns[io] = 280767;  fMeanV0A[io] = 45.1334;  fMeanV0C[io] = 66.8340;  fMeanV0M[io] = 111.862;  io++;
   fRuns[io] = 280768;  fMeanV0A[io] = 44.8909;  fMeanV0C[io] = 66.4428;  fMeanV0M[io] = 111.246;  io++;
   fRuns[io] = 280786;  fMeanV0A[io] = 45.1542;  fMeanV0C[io] = 66.8642;  fMeanV0M[io] = 111.913;  io++;
   fRuns[io] = 280787;  fMeanV0A[io] = 45.1296;  fMeanV0C[io] = 66.8469;  fMeanV0M[io] = 111.879;  io++;
   fRuns[io] = 280792;  fMeanV0A[io] = 44.9527;  fMeanV0C[io] = 66.7072;  fMeanV0M[io] = 111.554;  io++;
   fRuns[io] = 280793;  fMeanV0A[io] = 45.0339;  fMeanV0C[io] = 66.7113;  fMeanV0M[io] = 111.617;  io++;
   fRuns[io] = 280842;  fMeanV0A[io] = 45.1862;  fMeanV0C[io] = 66.9600;  fMeanV0M[io] = 112.038;  io++;
   fRuns[io] = 280844;  fMeanV0A[io] = 45.1921;  fMeanV0C[io] = 67.0385;  fMeanV0M[io] = 112.126;  io++;
   fRuns[io] = 280847;  fMeanV0A[io] = 45.1702;  fMeanV0C[io] = 66.9716;  fMeanV0M[io] = 112.038;  io++;
   fRuns[io] = 280848;  fMeanV0A[io] = 45.2566;  fMeanV0C[io] = 67.2557;  fMeanV0M[io] = 112.393;  io++;
   fRuns[io] = 280849;  fMeanV0A[io] = 45.3744;  fMeanV0C[io] = 67.1185;  fMeanV0M[io] = 112.400;  io++;
   fRuns[io] = 280854;  fMeanV0A[io] = 45.1432;  fMeanV0C[io] = 66.8773;  fMeanV0M[io] = 111.924;  io++;
   fRuns[io] = 280856;  fMeanV0A[io] = 45.3142;  fMeanV0C[io] = 66.9513;  fMeanV0M[io] = 112.181;  io++;
   fRuns[io] = 280880;  fMeanV0A[io] = 45.4410;  fMeanV0C[io] = 67.2230;  fMeanV0M[io] = 112.569;  io++;
   fRuns[io] = 280897;  fMeanV0A[io] = 44.9838;  fMeanV0C[io] = 66.6981;  fMeanV0M[io] = 111.560;  io++;
   fRuns[io] = 280936;  fMeanV0A[io] = 45.1734;  fMeanV0C[io] = 67.0298;  fMeanV0M[io] = 112.096;  io++;
   fRuns[io] = 280940;  fMeanV0A[io] = 45.3068;  fMeanV0C[io] = 67.2874;  fMeanV0M[io] = 112.489;  io++;
   fRuns[io] = 280943;  fMeanV0A[io] = 45.2504;  fMeanV0C[io] = 67.1454;  fMeanV0M[io] = 112.297;  io++;
   fRuns[io] = 280947;  fMeanV0A[io] = 45.1479;  fMeanV0C[io] = 66.8313;  fMeanV0M[io] = 111.869;  io++;
   fRuns[io] = 280990;  fMeanV0A[io] = 45.1624;  fMeanV0C[io] = 67.3422;  fMeanV0M[io] = 112.384;  io++;
   fRuns[io] = 280994;  fMeanV0A[io] = 44.8098;  fMeanV0C[io] = 66.8118;  fMeanV0M[io] = 111.534;  io++;
   fRuns[io] = 280996;  fMeanV0A[io] = 44.7443;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.357;  io++;
   fRuns[io] = 280997;  fMeanV0A[io] = 44.7766;  fMeanV0C[io] = 66.8905;  fMeanV0M[io] = 111.558;  io++;
   fRuns[io] = 280998;  fMeanV0A[io] = 44.7738;  fMeanV0C[io] = 66.7456;  fMeanV0M[io] = 111.414;  io++;
   fRuns[io] = 280999;  fMeanV0A[io] = 44.8665;  fMeanV0C[io] = 66.6767;  fMeanV0M[io] = 111.436;  io++;
   fRuns[io] = 281032;  fMeanV0A[io] = 44.8833;  fMeanV0C[io] = 66.7927;  fMeanV0M[io] = 111.564;  io++;
   fRuns[io] = 281033;  fMeanV0A[io] = 44.6934;  fMeanV0C[io] = 66.5216;  fMeanV0M[io] = 111.114;  io++;
   fRuns[io] = 281035;  fMeanV0A[io] = 44.8106;  fMeanV0C[io] = 66.7816;  fMeanV0M[io] = 111.494;  io++;
   fRuns[io] = 281036;  fMeanV0A[io] = 44.7005;  fMeanV0C[io] = 66.6556;  fMeanV0M[io] = 111.244;  io++;
   fRuns[io] = 281060;  fMeanV0A[io] = 45.5495;  fMeanV0C[io] = 67.8417;  fMeanV0M[io] = 113.300;  io++;
   fRuns[io] = 281061;  fMeanV0A[io] = 45.0308;  fMeanV0C[io] = 67.5032;  fMeanV0M[io] = 112.433;  io++;
   fRuns[io] = 281062;  fMeanV0A[io] = 44.9673;  fMeanV0C[io] = 67.2450;  fMeanV0M[io] = 112.112;  io++;
   fRuns[io] = 281080;  fMeanV0A[io] = 44.8555;  fMeanV0C[io] = 66.9124;  fMeanV0M[io] = 111.687;  io++;
   fRuns[io] = 281081;  fMeanV0A[io] = 44.9134;  fMeanV0C[io] = 66.7856;  fMeanV0M[io] = 111.596;  io++;
   fRuns[io] = 281179;  fMeanV0A[io] = 45.1292;  fMeanV0C[io] = 66.8856;  fMeanV0M[io] = 111.915;  io++;
   fRuns[io] = 281180;  fMeanV0A[io] = 45.1596;  fMeanV0C[io] = 66.8875;  fMeanV0M[io] = 111.948;  io++;
   fRuns[io] = 281181;  fMeanV0A[io] = 44.8743;  fMeanV0C[io] = 66.8456;  fMeanV0M[io] = 111.618;  io++;
   fRuns[io] = 281189;  fMeanV0A[io] = 45.0556;  fMeanV0C[io] = 66.8877;  fMeanV0M[io] = 111.839;  io++;
   fRuns[io] = 281190;  fMeanV0A[io] = 44.9675;  fMeanV0C[io] = 67.0752;  fMeanV0M[io] = 111.939;  io++;
   fRuns[io] = 281191;  fMeanV0A[io] = 44.9711;  fMeanV0C[io] = 66.8247;  fMeanV0M[io] = 111.697;  io++;
   fRuns[io] = 281212;  fMeanV0A[io] = 45.0640;  fMeanV0C[io] = 66.9719;  fMeanV0M[io] = 111.936;  io++;
   fRuns[io] = 281213;  fMeanV0A[io] = 45.1604;  fMeanV0C[io] = 66.9935;  fMeanV0M[io] = 112.054;  io++;
   fRuns[io] = 281240;  fMeanV0A[io] = 45.0523;  fMeanV0C[io] = 66.9582;  fMeanV0M[io] = 111.901;  io++;
   fRuns[io] = 281241;  fMeanV0A[io] = 44.9396;  fMeanV0C[io] = 66.8919;  fMeanV0M[io] = 111.727;  io++;
   fRuns[io] = 281242;  fMeanV0A[io] = 44.9527;  fMeanV0C[io] = 66.8370;  fMeanV0M[io] = 111.675;  io++;
   fRuns[io] = 281243;  fMeanV0A[io] = 44.9827;  fMeanV0C[io] = 66.7967;  fMeanV0M[io] = 111.674;  io++;
   fRuns[io] = 281244;  fMeanV0A[io] = 45.0738;  fMeanV0C[io] = 66.7527;  fMeanV0M[io] = 111.722;  io++;
   fRuns[io] = 281271;  fMeanV0A[io] = 45.2303;  fMeanV0C[io] = 67.3252;  fMeanV0M[io] = 112.472;  io++;
   fRuns[io] = 281273;  fMeanV0A[io] = 45.0247;  fMeanV0C[io] = 66.8344;  fMeanV0M[io] = 111.751;  io++;
   fRuns[io] = 281275;  fMeanV0A[io] = 44.9913;  fMeanV0C[io] = 66.7698;  fMeanV0M[io] = 111.680;  io++;
   fRuns[io] = 281277;  fMeanV0A[io] = 44.8939;  fMeanV0C[io] = 66.7281;  fMeanV0M[io] = 111.520;  io++;
   fRuns[io] = 281301;  fMeanV0A[io] = 45.0142;  fMeanV0C[io] = 66.8266;  fMeanV0M[io] = 111.731;  io++;
   fRuns[io] = 281321;  fMeanV0A[io] = 45.0522;  fMeanV0C[io] = 66.7500;  fMeanV0M[io] = 111.694;  io++;
   fRuns[io] = 281415;  fMeanV0A[io] = 44.9017;  fMeanV0C[io] = 67.3152;  fMeanV0M[io] = 112.107;  io++;
   fRuns[io] = 281441;  fMeanV0A[io] = 44.9062;  fMeanV0C[io] = 66.7254;  fMeanV0M[io] = 111.530;  io++;
   fRuns[io] = 281443;  fMeanV0A[io] = 44.8921;  fMeanV0C[io] = 66.7150;  fMeanV0M[io] = 111.495;  io++;
   fRuns[io] = 281444;  fMeanV0A[io] = 44.9984;  fMeanV0C[io] = 66.8996;  fMeanV0M[io] = 111.810;  io++;
   fRuns[io] = 281446;  fMeanV0A[io] = 45.0362;  fMeanV0C[io] = 66.7819;  fMeanV0M[io] = 111.714;  io++;
   fRuns[io] = 281449;  fMeanV0A[io] = 44.9713;  fMeanV0C[io] = 66.8942;  fMeanV0M[io] = 111.779;  io++;
   fRuns[io] = 281450;  fMeanV0A[io] = 44.9411;  fMeanV0C[io] = 66.7430;  fMeanV0M[io] = 111.580;  io++;
   fRuns[io] = 281475;  fMeanV0A[io] = 44.9852;  fMeanV0C[io] = 66.7944;  fMeanV0M[io] = 111.667;  io++;
   fRuns[io] = 281477;  fMeanV0A[io] = 45.0283;  fMeanV0C[io] = 66.6172;  fMeanV0M[io] = 111.545;  io++;
   fRuns[io] = 281509;  fMeanV0A[io] = 44.9428;  fMeanV0C[io] = 66.7099;  fMeanV0M[io] = 111.543;  io++;
   fRuns[io] = 281511;  fMeanV0A[io] = 44.8783;  fMeanV0C[io] = 66.7271;  fMeanV0M[io] = 111.480;  io++;
   fRuns[io] = 281557;  fMeanV0A[io] = 45.0998;  fMeanV0C[io] = 66.6082;  fMeanV0M[io] = 111.595;  io++;
   fRuns[io] = 281562;  fMeanV0A[io] = 44.9329;  fMeanV0C[io] = 66.7540;  fMeanV0M[io] = 111.579;  io++;
   fRuns[io] = 281563;  fMeanV0A[io] = 44.9119;  fMeanV0C[io] = 66.5690;  fMeanV0M[io] = 111.370;  io++;
   fRuns[io] = 281568;  fMeanV0A[io] = 44.9871;  fMeanV0C[io] = 66.6209;  fMeanV0M[io] = 111.511;  io++;
   fRuns[io] = 281569;  fMeanV0A[io] = 44.8148;  fMeanV0C[io] = 66.4740;  fMeanV0M[io] = 111.177;  io++;
   fRuns[io] = 281574;  fMeanV0A[io] = 44.8804;  fMeanV0C[io] = 66.3797;  fMeanV0M[io] = 111.173;  io++;
   fRuns[io] = 281583;  fMeanV0A[io] = 44.8138;  fMeanV0C[io] = 66.6529;  fMeanV0M[io] = 111.358;  io++;
   fRuns[io] = 281592;  fMeanV0A[io] = 44.6784;  fMeanV0C[io] = 66.5507;  fMeanV0M[io] = 111.133;  io++;
   fRuns[io] = 281633;  fMeanV0A[io] = 44.9260;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.540;  io++;
   fRuns[io] = 281892;  fMeanV0A[io] = 44.8142;  fMeanV0C[io] = 66.1647;  fMeanV0M[io] = 110.897;  io++;
   fRuns[io] = 281893;  fMeanV0A[io] = 44.3514;  fMeanV0C[io] = 65.6740;  fMeanV0M[io] = 109.910;  io++;
   fRuns[io] = 281894;  fMeanV0A[io] = 44.3902;  fMeanV0C[io] = 65.7523;  fMeanV0M[io] = 110.032;  io++;
   fRuns[io] = 281895;  fMeanV0A[io] = 44.3739;  fMeanV0C[io] = 65.7561;  fMeanV0M[io] = 110.028;  io++;
   fRuns[io] = 281915;  fMeanV0A[io] = 44.6005;  fMeanV0C[io] = 65.8568;  fMeanV0M[io] = 110.346;  io++;
   fRuns[io] = 281916;  fMeanV0A[io] = 44.4004;  fMeanV0C[io] = 65.7322;  fMeanV0M[io] = 110.021;  io++;
   fRuns[io] = 281918;  fMeanV0A[io] = 44.5497;  fMeanV0C[io] = 65.8153;  fMeanV0M[io] = 110.269;  io++;
   fRuns[io] = 281920;  fMeanV0A[io] = 44.5406;  fMeanV0C[io] = 65.7390;  fMeanV0M[io] = 110.165;  io++;
   fRuns[io] = 281928;  fMeanV0A[io] = 44.4097;  fMeanV0C[io] = 65.5548;  fMeanV0M[io] = 109.858;  io++;
   fRuns[io] = 281931;  fMeanV0A[io] = 44.6852;  fMeanV0C[io] = 66.1994;  fMeanV0M[io] = 110.783;  io++;
   fRuns[io] = 281932;  fMeanV0A[io] = 44.4489;  fMeanV0C[io] = 65.7691;  fMeanV0M[io] = 110.110;  io++;
   fRuns[io] = 281939;  fMeanV0A[io] = 44.5202;  fMeanV0C[io] = 65.7536;  fMeanV0M[io] = 110.159;  io++;
   fRuns[io] = 281940;  fMeanV0A[io] = 44.5177;  fMeanV0C[io] = 65.8935;  fMeanV0M[io] = 110.300;  io++;
   fRuns[io] = 281953;  fMeanV0A[io] = 44.5299;  fMeanV0C[io] = 65.8725;  fMeanV0M[io] = 110.293;  io++;
   fRuns[io] = 281956;  fMeanV0A[io] = 44.3228;  fMeanV0C[io] = 65.7231;  fMeanV0M[io] = 109.936;  io++;
   fRuns[io] = 281961;  fMeanV0A[io] = 44.4486;  fMeanV0C[io] = 65.7437;  fMeanV0M[io] = 110.085;  io++;

   // FILTER_p-p_208_LHC17r
   fRuns[io] = 282528;  fMeanV0A[io] = 44.5751;  fMeanV0C[io] = 65.8345;  fMeanV0M[io] = 110.308;  io++;
   fRuns[io] = 282544;  fMeanV0A[io] = 44.2678;  fMeanV0C[io] = 65.4773;  fMeanV0M[io] = 109.627;  io++;
   fRuns[io] = 282545;  fMeanV0A[io] = 44.2855;  fMeanV0C[io] = 65.4964;  fMeanV0M[io] = 109.669;  io++;
   fRuns[io] = 282546;  fMeanV0A[io] = 44.1408;  fMeanV0C[io] = 65.4330;  fMeanV0M[io] = 109.473;  io++;
   fRuns[io] = 282573;  fMeanV0A[io] = 44.3559;  fMeanV0C[io] = 65.5686;  fMeanV0M[io] = 109.810;  io++;
   fRuns[io] = 282575;  fMeanV0A[io] = 44.3761;  fMeanV0C[io] = 65.6121;  fMeanV0M[io] = 109.880;  io++;
   fRuns[io] = 282579;  fMeanV0A[io] = 44.3710;  fMeanV0C[io] = 65.6803;  fMeanV0M[io] = 109.925;  io++;
   fRuns[io] = 282580;  fMeanV0A[io] = 44.4169;  fMeanV0C[io] = 65.5248;  fMeanV0M[io] = 109.834;  io++;
   fRuns[io] = 282606;  fMeanV0A[io] = 44.5594;  fMeanV0C[io] = 65.6523;  fMeanV0M[io] = 110.111;  io++;
   fRuns[io] = 282607;  fMeanV0A[io] = 44.3277;  fMeanV0C[io] = 65.4114;  fMeanV0M[io] = 109.629;  io++;
   fRuns[io] = 282608;  fMeanV0A[io] = 44.2839;  fMeanV0C[io] = 65.4020;  fMeanV0M[io] = 109.576;  io++;
   fRuns[io] = 282609;  fMeanV0A[io] = 44.3728;  fMeanV0C[io] = 65.5314;  fMeanV0M[io] = 109.797;  io++;
   fRuns[io] = 282618;  fMeanV0A[io] = 44.3684;  fMeanV0C[io] = 65.6726;  fMeanV0M[io] = 109.939;  io++;
   fRuns[io] = 282620;  fMeanV0A[io] = 44.4642;  fMeanV0C[io] = 65.4995;  fMeanV0M[io] = 109.856;  io++;
   fRuns[io] = 282622;  fMeanV0A[io] = 44.3402;  fMeanV0C[io] = 65.5665;  fMeanV0M[io] = 109.810;  io++;
   fRuns[io] = 282629;  fMeanV0A[io] = 44.3156;  fMeanV0C[io] = 65.4220;  fMeanV0M[io] = 109.624;  io++;
   fRuns[io] = 282651;  fMeanV0A[io] = 44.6365;  fMeanV0C[io] = 65.9243;  fMeanV0M[io] = 110.423;  io++;
   fRuns[io] = 282666;  fMeanV0A[io] = 44.3446;  fMeanV0C[io] = 65.6009;  fMeanV0M[io] = 109.835;  io++;
   fRuns[io] = 282667;  fMeanV0A[io] = 44.3031;  fMeanV0C[io] = 65.5014;  fMeanV0M[io] = 109.694;  io++;
   fRuns[io] = 282670;  fMeanV0A[io] = 44.3852;  fMeanV0C[io] = 65.4583;  fMeanV0M[io] = 109.732;  io++;
   fRuns[io] = 282671;  fMeanV0A[io] = 44.2074;  fMeanV0C[io] = 65.5636;  fMeanV0M[io] = 109.665;  io++;
   fRuns[io] = 282673;  fMeanV0A[io] = 44.4285;  fMeanV0C[io] = 65.7165;  fMeanV0M[io] = 110.033;  io++;
   fRuns[io] = 282676;  fMeanV0A[io] = 44.4527;  fMeanV0C[io] = 65.5462;  fMeanV0M[io] = 109.895;  io++;
   fRuns[io] = 282677;  fMeanV0A[io] = 44.2885;  fMeanV0C[io] = 65.5971;  fMeanV0M[io] = 109.776;  io++;
   fRuns[io] = 282700;  fMeanV0A[io] = 44.5101;  fMeanV0C[io] = 65.8797;  fMeanV0M[io] = 110.289;  io++;
   fRuns[io] = 282702;  fMeanV0A[io] = 44.3138;  fMeanV0C[io] = 65.7604;  fMeanV0M[io] = 109.971;  io++;
   fRuns[io] = 282703;  fMeanV0A[io] = 44.2851;  fMeanV0C[io] = 65.5566;  fMeanV0M[io] = 109.735;  io++;
   fRuns[io] = 282704;  fMeanV0A[io] = 44.3626;  fMeanV0C[io] = 65.4834;  fMeanV0M[io] = 109.741;  io++;

   // FILTER_p-p_208_LHC18b
   fRuns[io] = 285009;  fMeanV0A[io] = 55.2521;  fMeanV0C[io] = 81.8254;  fMeanV0M[io] = 136.997;  io++;
   fRuns[io] = 285011;  fMeanV0A[io] = 13.5845;  fMeanV0C[io] = 12.608;   fMeanV0M[io] = 26.0796;  io++;
   fRuns[io] = 285012;  fMeanV0A[io] = 9.16816;  fMeanV0C[io] = 7.7265;   fMeanV0M[io] = 17.1483;  io++;
   fRuns[io] = 285013;  fMeanV0A[io] = 33.9831;  fMeanV0C[io] = 43.7996;  fMeanV0M[io] = 77.6915;  io++;
   fRuns[io] = 285014;  fMeanV0A[io] = 34.0567;  fMeanV0C[io] = 44.1135;  fMeanV0M[io] = 78.0708;  io++;
   fRuns[io] = 285015;  fMeanV0A[io] = 34.2392;  fMeanV0C[io] = 44.4197;  fMeanV0M[io] = 78.5631;  io++;
   fRuns[io] = 285064;  fMeanV0A[io] = 34.0049;  fMeanV0C[io] = 43.6075;  fMeanV0M[io] = 77.4904;  io++;
   fRuns[io] = 285065;  fMeanV0A[io] = 33.5307;  fMeanV0C[io] = 42.9036;  fMeanV0M[io] = 76.3311;  io++;
   fRuns[io] = 285066;  fMeanV0A[io] = 33.486;   fMeanV0C[io] = 42.9486;  fMeanV0M[io] = 76.3393;  io++;
   fRuns[io] = 285106;  fMeanV0A[io] = 33.0943;  fMeanV0C[io] = 42.3562;  fMeanV0M[io] = 75.3466;  io++;
   fRuns[io] = 285108;  fMeanV0A[io] = 33.0744;  fMeanV0C[io] = 42.237;   fMeanV0M[io] = 75.2164;  io++;
   fRuns[io] = 285125;  fMeanV0A[io] = 33.2548;  fMeanV0C[io] = 42.4402;  fMeanV0M[io] = 75.5897;  io++;
   fRuns[io] = 285127;  fMeanV0A[io] = 33.608;   fMeanV0C[io] = 43.0354;  fMeanV0M[io] = 76.5353;  io++;
   fRuns[io] = 285165;  fMeanV0A[io] = 32.8191;  fMeanV0C[io] = 41.9321;  fMeanV0M[io] = 74.6483;  io++;
   fRuns[io] = 285200;  fMeanV0A[io] = 34.0951;  fMeanV0C[io] = 43.9659;  fMeanV0M[io] = 77.9468;  io++;
   fRuns[io] = 285202;  fMeanV0A[io] = 34.2534;  fMeanV0C[io] = 44.0942;  fMeanV0M[io] = 78.2445;  io++;
   fRuns[io] = 285203;  fMeanV0A[io] = 34.137;   fMeanV0C[io] = 43.8929;  fMeanV0M[io] = 77.9142;  io++;
   fRuns[io] = 285222;  fMeanV0A[io] = 54.4815;  fMeanV0C[io] = 80.3744;  fMeanV0M[io] = 134.766;  io++;
   fRuns[io] = 285224;  fMeanV0A[io] = 54.5747;  fMeanV0C[io] = 80.7792;  fMeanV0M[io] = 135.286;  io++;
   fRuns[io] = 285327;  fMeanV0A[io] = 55.1115;  fMeanV0C[io] = 80.9864;  fMeanV0M[io] = 136.026;  io++;
   fRuns[io] = 285328;  fMeanV0A[io] = 54.8856;  fMeanV0C[io] = 80.7743;  fMeanV0M[io] = 135.59;  io++;
   fRuns[io] = 285347;  fMeanV0A[io] = 55.2806;  fMeanV0C[io] = 81.0572;  fMeanV0M[io] = 136.26;  io++;
   fRuns[io] = 285364;  fMeanV0A[io] = 55.2525;  fMeanV0C[io] = 81.1235;  fMeanV0M[io] = 136.306;  io++;
   fRuns[io] = 285365;  fMeanV0A[io] = 54.8061;  fMeanV0C[io] = 80.8271;  fMeanV0M[io] = 135.551;  io++;
   fRuns[io] = 285396;  fMeanV0A[io] = 54.7919;  fMeanV0C[io] = 80.453;   fMeanV0M[io] = 135.162;  io++;

   // FILTER_p-p_208_LHC18d
   fRuns[io] = 285978;  fMeanV0A[io] = 52.6167;  fMeanV0C[io] = 78.1349;  fMeanV0M[io] = 130.667;  io++;
   fRuns[io] = 285979;  fMeanV0A[io] = 52.6559;  fMeanV0C[io] = 78.1399;  fMeanV0M[io] = 130.726;  io++;
   fRuns[io] = 285980;  fMeanV0A[io] = 52.3108;  fMeanV0C[io] = 77.921;   fMeanV0M[io] = 130.154;  io++;
   fRuns[io] = 286014;  fMeanV0A[io] = 52.5196;  fMeanV0C[io] = 78.0187;  fMeanV0M[io] = 130.451;  io++;
   fRuns[io] = 286025;  fMeanV0A[io] = 52.2718;  fMeanV0C[io] = 78.102;   fMeanV0M[io] = 130.307;  io++;
   fRuns[io] = 286064;  fMeanV0A[io] = 52.4505;  fMeanV0C[io] = 77.9472;  fMeanV0M[io] = 130.319;  io++;
   fRuns[io] = 286124;  fMeanV0A[io] = 52.3515;  fMeanV0C[io] = 77.9117;  fMeanV0M[io] = 130.177;  io++;
   fRuns[io] = 286127;  fMeanV0A[io] = 52.4224;  fMeanV0C[io] = 77.78;    fMeanV0M[io] = 130.123;  io++;
   fRuns[io] = 286129;  fMeanV0A[io] = 52.5074;  fMeanV0C[io] = 78.0631;  fMeanV0M[io] = 130.494;  io++;
   fRuns[io] = 286130;  fMeanV0A[io] = 52.6798;  fMeanV0C[io] = 78.1209;  fMeanV0M[io] = 130.714;  io++;
   fRuns[io] = 286159;  fMeanV0A[io] = 52.4128;  fMeanV0C[io] = 77.7735;  fMeanV0M[io] = 130.102;  io++;
   fRuns[io] = 286198;  fMeanV0A[io] = 52.542;   fMeanV0C[io] = 77.896;   fMeanV0M[io] = 130.354;  io++;
   fRuns[io] = 286199;  fMeanV0A[io] = 52.2918;  fMeanV0C[io] = 77.7458;  fMeanV0M[io] = 129.955;  io++;
   fRuns[io] = 286201;  fMeanV0A[io] = 52.3478;  fMeanV0C[io] = 77.7447;  fMeanV0M[io] = 130.014;  io++;
   fRuns[io] = 286202;  fMeanV0A[io] = 52.2706;  fMeanV0C[io] = 77.7824;  fMeanV0M[io] = 129.974;  io++;
   fRuns[io] = 286203;  fMeanV0A[io] = 52.3935;  fMeanV0C[io] = 77.8652;  fMeanV0M[io] = 130.178;  io++;
   fRuns[io] = 286229;  fMeanV0A[io] = 52.3883;  fMeanV0C[io] = 77.6589;  fMeanV0M[io] = 129.968;  io++;
   fRuns[io] = 286230;  fMeanV0A[io] = 52.2572;  fMeanV0C[io] = 77.6396;  fMeanV0M[io] = 129.805;  io++;
   fRuns[io] = 286231;  fMeanV0A[io] = 52.366;   fMeanV0C[io] = 77.7735;  fMeanV0M[io] = 130.053;  io++;
   fRuns[io] = 286254;  fMeanV0A[io] = 52.6608;  fMeanV0C[io] = 77.9607;  fMeanV0M[io] = 130.512;  io++;
   fRuns[io] = 286257;  fMeanV0A[io] = 52.2977;  fMeanV0C[io] = 77.6823;  fMeanV0M[io] = 129.899;  io++;
   fRuns[io] = 286258;  fMeanV0A[io] = 52.1941;  fMeanV0C[io] = 77.5952;  fMeanV0M[io] = 129.712;  io++;
   fRuns[io] = 286261;  fMeanV0A[io] = 52.2966;  fMeanV0C[io] = 77.3074;  fMeanV0M[io] = 129.504;  io++;
   fRuns[io] = 286263;  fMeanV0A[io] = 52.25;    fMeanV0C[io] = 77.8197;  fMeanV0M[io] = 129.993;  io++;
   fRuns[io] = 286282;  fMeanV0A[io] = 52.2815;  fMeanV0C[io] = 77.5496;  fMeanV0M[io] = 129.754;  io++;
   fRuns[io] = 286284;  fMeanV0A[io] = 52.2809;  fMeanV0C[io] = 77.5871;  fMeanV0M[io] = 129.783;  io++;
   fRuns[io] = 286287;  fMeanV0A[io] = 52.2406;  fMeanV0C[io] = 77.5238;  fMeanV0M[io] = 129.691;  io++;
   fRuns[io] = 286288;  fMeanV0A[io] = 52.1208;  fMeanV0C[io] = 77.4366;  fMeanV0M[io] = 129.475;  io++;
   fRuns[io] = 286289;  fMeanV0A[io] = 52.1161;  fMeanV0C[io] = 77.5717;  fMeanV0M[io] = 129.616;  io++;
   fRuns[io] = 286308;  fMeanV0A[io] = 52.3949;  fMeanV0C[io] = 77.788;   fMeanV0M[io] = 130.1;  io++;
   fRuns[io] = 286309;  fMeanV0A[io] = 52.2685;  fMeanV0C[io] = 77.6832;  fMeanV0M[io] = 129.844;  io++;
   fRuns[io] = 286310;  fMeanV0A[io] = 52.2154;  fMeanV0C[io] = 77.4518;  fMeanV0M[io] = 129.576;  io++;
   fRuns[io] = 286311;  fMeanV0A[io] = 52.4999;  fMeanV0C[io] = 77.9524;  fMeanV0M[io] = 130.393;  io++;
   fRuns[io] = 286312;  fMeanV0A[io] = 52.0642;  fMeanV0C[io] = 77.4978;  fMeanV0M[io] = 129.491;  io++;
   fRuns[io] = 286313;  fMeanV0A[io] = 52.0934;  fMeanV0C[io] = 77.4069;  fMeanV0M[io] = 129.416;  io++;
   fRuns[io] = 286314;  fMeanV0A[io] = 52.1274;  fMeanV0C[io] = 77.6984;  fMeanV0M[io] = 129.765;  io++;
   fRuns[io] = 286336;  fMeanV0A[io] = 52.0606;  fMeanV0C[io] = 77.3254;  fMeanV0M[io] = 129.306;  io++;
   fRuns[io] = 286337;  fMeanV0A[io] = 52.0064;  fMeanV0C[io] = 76.9894;  fMeanV0M[io] = 128.92;  io++;
   fRuns[io] = 286340;  fMeanV0A[io] = 51.8017;  fMeanV0C[io] = 77.1237;  fMeanV0M[io] = 128.824;  io++;
   fRuns[io] = 286341;  fMeanV0A[io] = 52.0006;  fMeanV0C[io] = 77.2122;  fMeanV0M[io] = 129.125;  io++;
   fRuns[io] = 286345;  fMeanV0A[io] = 52.0381;  fMeanV0C[io] = 77.1804;  fMeanV0M[io] = 129.14;  io++;
   fRuns[io] = 286348;  fMeanV0A[io] = 51.9719;  fMeanV0C[io] = 77.1224;  fMeanV0M[io] = 129.013;  io++;
   fRuns[io] = 286349;  fMeanV0A[io] = 51.9632;  fMeanV0C[io] = 77.1958;  fMeanV0M[io] = 129.087;  io++;
   fRuns[io] = 286350;  fMeanV0A[io] = 52.0582;  fMeanV0C[io] = 77.2;     fMeanV0M[io] = 129.175;  io++;

   // FILTER_p-p_208_LHC18e
   fRuns[io] = 286380;  fMeanV0A[io] = 52.2194;  fMeanV0C[io] = 77.6518;  fMeanV0M[io] = 129.796;  io++;
   fRuns[io] = 286426;  fMeanV0A[io] = 52.5906;  fMeanV0C[io] = 77.6302;  fMeanV0M[io] = 130.151;  io++;
   fRuns[io] = 286427;  fMeanV0A[io] = 52.2609;  fMeanV0C[io] = 77.3846;  fMeanV0M[io] = 129.559;  io++;
   fRuns[io] = 286428;  fMeanV0A[io] = 52.187;   fMeanV0C[io] = 77.4279;  fMeanV0M[io] = 129.538;  io++;
   fRuns[io] = 286454;  fMeanV0A[io] = 52.3228;  fMeanV0C[io] = 77.6131;  fMeanV0M[io] = 129.852;  io++;
   fRuns[io] = 286455;  fMeanV0A[io] = 52.2174;  fMeanV0C[io] = 77.5624;  fMeanV0M[io] = 129.706;  io++;
   fRuns[io] = 286482;  fMeanV0A[io] = 52.5504;  fMeanV0C[io] = 77.8364;  fMeanV0M[io] = 130.333;  io++;
   fRuns[io] = 286502;  fMeanV0A[io] = 52.2199;  fMeanV0C[io] = 77.3987;  fMeanV0M[io] = 129.538;  io++;
   fRuns[io] = 286508;  fMeanV0A[io] = 52.134;   fMeanV0C[io] = 77.3175;  fMeanV0M[io] = 129.359;  io++;
   fRuns[io] = 286509;  fMeanV0A[io] = 52.0677;  fMeanV0C[io] = 77.2245;  fMeanV0M[io] = 129.206;  io++;
   fRuns[io] = 286511;  fMeanV0A[io] = 51.9751;  fMeanV0C[io] = 77.7522;  fMeanV0M[io] = 129.592;  io++;
   fRuns[io] = 286566;  fMeanV0A[io] = 52.7348;  fMeanV0C[io] = 77.6577;  fMeanV0M[io] = 130.3;  io++;
   fRuns[io] = 286567;  fMeanV0A[io] = 52.2842;  fMeanV0C[io] = 77.4272;  fMeanV0M[io] = 129.619;  io++;
   fRuns[io] = 286568;  fMeanV0A[io] = 52.1079;  fMeanV0C[io] = 77.2863;  fMeanV0M[io] = 129.312;  io++;
   fRuns[io] = 286569;  fMeanV0A[io] = 52.1655;  fMeanV0C[io] = 77.0853;  fMeanV0M[io] = 129.166;  io++;
   fRuns[io] = 286591;  fMeanV0A[io] = 52.2678;  fMeanV0C[io] = 77.8534;  fMeanV0M[io] = 130.004;  io++;
   fRuns[io] = 286592;  fMeanV0A[io] = 52.1537;  fMeanV0C[io] = 77.2772;  fMeanV0M[io] = 129.348;  io++;
   fRuns[io] = 286633;  fMeanV0A[io] = 53.2105;  fMeanV0C[io] = 79.1528;  fMeanV0M[io] = 132.285;  io++;
   fRuns[io] = 286653;  fMeanV0A[io] = 51.9812;  fMeanV0C[io] = 77.0928;  fMeanV0M[io] = 128.984;  io++;
   fRuns[io] = 286661;  fMeanV0A[io] = 51.9539;  fMeanV0C[io] = 77.0664;  fMeanV0M[io] = 128.928;  io++;
   fRuns[io] = 286695;  fMeanV0A[io] = 52.2161;  fMeanV0C[io] = 77.2537;  fMeanV0M[io] = 129.389;  io++;
   fRuns[io] = 286731;  fMeanV0A[io] = 52.288;   fMeanV0C[io] = 77.5656;  fMeanV0M[io] = 129.774;  io++;
   fRuns[io] = 286799;  fMeanV0A[io] = 52.1701;  fMeanV0C[io] = 77.1681;  fMeanV0M[io] = 129.257;  io++;
   fRuns[io] = 286801;  fMeanV0A[io] = 52.0286;  fMeanV0C[io] = 77.058;   fMeanV0M[io] = 129.004;  io++;
   fRuns[io] = 286805;  fMeanV0A[io] = 52.0151;  fMeanV0C[io] = 77.0511;  fMeanV0M[io] = 128.991;  io++;
   fRuns[io] = 286809;  fMeanV0A[io] = 52.2227;  fMeanV0C[io] = 77.6888;  fMeanV0M[io] = 129.815;  io++;
   fRuns[io] = 286846;  fMeanV0A[io] = 52.2838;  fMeanV0C[io] = 77.3642;  fMeanV0M[io] = 129.569;  io++;
   fRuns[io] = 286850;  fMeanV0A[io] = 51.9245;  fMeanV0C[io] = 77.1487;  fMeanV0M[io] = 128.992;  io++;
   fRuns[io] = 286852;  fMeanV0A[io] = 51.9804;  fMeanV0C[io] = 77.154;   fMeanV0M[io] = 129.048;  io++;
   fRuns[io] = 286874;  fMeanV0A[io] = 52.043;   fMeanV0C[io] = 77.1783;  fMeanV0M[io] = 129.135;  io++;
   fRuns[io] = 286876;  fMeanV0A[io] = 51.9762;  fMeanV0C[io] = 77.0462;  fMeanV0M[io] = 128.94;  io++;
   fRuns[io] = 286877;  fMeanV0A[io] = 52.0077;  fMeanV0C[io] = 77.0784;  fMeanV0M[io] = 128.991;  io++;
   fRuns[io] = 286907;  fMeanV0A[io] = 52.0239;  fMeanV0C[io] = 77.0533;  fMeanV0M[io] = 129.002;  io++;
   fRuns[io] = 286910;  fMeanV0A[io] = 52.0222;  fMeanV0C[io] = 76.9398;  fMeanV0M[io] = 128.876;  io++;
   fRuns[io] = 286911;  fMeanV0A[io] = 51.8708;  fMeanV0C[io] = 76.8903;  fMeanV0M[io] = 128.683;  io++;
   fRuns[io] = 286930;  fMeanV0A[io] = 51.9631;  fMeanV0C[io] = 77.0681;  fMeanV0M[io] = 128.943;  io++;
   fRuns[io] = 286931;  fMeanV0A[io] = 51.7828;  fMeanV0C[io] = 76.9136;  fMeanV0M[io] = 128.616;  io++;
   fRuns[io] = 286932;  fMeanV0A[io] = 51.7154;  fMeanV0C[io] = 76.743;   fMeanV0M[io] = 128.381;  io++;
   fRuns[io] = 286933;  fMeanV0A[io] = 51.7479;  fMeanV0C[io] = 76.875;   fMeanV0M[io] = 128.532;  io++;
   fRuns[io] = 286936;  fMeanV0A[io] = 51.7394;  fMeanV0C[io] = 76.5765;  fMeanV0M[io] = 128.229;  io++;
   fRuns[io] = 286937;  fMeanV0A[io] = 51.6711;  fMeanV0C[io] = 76.8347;  fMeanV0M[io] = 128.414;  io++;


   // FILTER_p-p_208_LHC18f
   fRuns[io] = 287000;  fMeanV0A[io] = 51.7496;  fMeanV0C[io] = 76.8361;  fMeanV0M[io] = 128.5;  io++;
   fRuns[io] = 287021;  fMeanV0A[io] = 52.124;   fMeanV0C[io] = 77.1482;  fMeanV0M[io] = 129.201;  io++;
   fRuns[io] = 287063;  fMeanV0A[io] = 52.0098;  fMeanV0C[io] = 77.1583;  fMeanV0M[io] = 129.076;  io++;
   fRuns[io] = 287064;  fMeanV0A[io] = 51.8192;  fMeanV0C[io] = 77.0838;  fMeanV0M[io] = 128.825;  io++;
   fRuns[io] = 287066;  fMeanV0A[io] = 51.8184;  fMeanV0C[io] = 77.0215;  fMeanV0M[io] = 128.76;  io++;
   fRuns[io] = 287071;  fMeanV0A[io] = 51.8689;  fMeanV0C[io] = 77.0099;  fMeanV0M[io] = 128.797;  io++;
   fRuns[io] = 287072;  fMeanV0A[io] = 51.7485;  fMeanV0C[io] = 76.8152;  fMeanV0M[io] = 128.481;  io++;
   fRuns[io] = 287077;  fMeanV0A[io] = 51.713;   fMeanV0C[io] = 76.9792;  fMeanV0M[io] = 128.602;  io++;
   fRuns[io] = 287137;  fMeanV0A[io] = 51.8783;  fMeanV0C[io] = 77.1173;  fMeanV0M[io] = 128.908;  io++;
   fRuns[io] = 287155;  fMeanV0A[io] = 51.7863;  fMeanV0C[io] = 77.0036;  fMeanV0M[io] = 128.703;  io++;
   fRuns[io] = 287185;  fMeanV0A[io] = 52.5877;  fMeanV0C[io] = 77.9283;  fMeanV0M[io] = 130.471;  io++;
   fRuns[io] = 287201;  fMeanV0A[io] = 51.7921;  fMeanV0C[io] = 76.9224;  fMeanV0M[io] = 128.622;  io++;
   fRuns[io] = 287202;  fMeanV0A[io] = 51.7848;  fMeanV0C[io] = 76.9209;  fMeanV0M[io] = 128.621;  io++;
   fRuns[io] = 287203;  fMeanV0A[io] = 51.8239;  fMeanV0C[io] = 76.9434;  fMeanV0M[io] = 128.685;  io++;
   fRuns[io] = 287204;  fMeanV0A[io] = 51.8126;  fMeanV0C[io] = 76.847;   fMeanV0M[io] = 128.571;  io++;
   fRuns[io] = 287208;  fMeanV0A[io] = 51.8501;  fMeanV0C[io] = 77.1585;  fMeanV0M[io] = 128.937;  io++;
   fRuns[io] = 287209;  fMeanV0A[io] = 51.8737;  fMeanV0C[io] = 76.7923;  fMeanV0M[io] = 128.592;  io++;
   fRuns[io] = 287248;  fMeanV0A[io] = 51.9153;  fMeanV0C[io] = 77.0536;  fMeanV0M[io] = 128.894;  io++;
   fRuns[io] = 287249;  fMeanV0A[io] = 51.7987;  fMeanV0C[io] = 76.9959;  fMeanV0M[io] = 128.71;  io++;
   fRuns[io] = 287250;  fMeanV0A[io] = 51.7679;  fMeanV0C[io] = 76.8606;  fMeanV0M[io] = 128.559;  io++;
   fRuns[io] = 287251;  fMeanV0A[io] = 51.8384;  fMeanV0C[io] = 77.0086;  fMeanV0M[io] = 128.77;  io++;
   fRuns[io] = 287254;  fMeanV0A[io] = 51.8829;  fMeanV0C[io] = 76.9351;  fMeanV0M[io] = 128.733;  io++;
   fRuns[io] = 287283;  fMeanV0A[io] = 51.7557;  fMeanV0C[io] = 77.156;   fMeanV0M[io] = 128.826;  io++;
   fRuns[io] = 287323;  fMeanV0A[io] = 52.4941;  fMeanV0C[io] = 77.5554;  fMeanV0M[io] = 129.959;  io++;
   fRuns[io] = 287324;  fMeanV0A[io] = 52.3113;  fMeanV0C[io] = 77.6721;  fMeanV0M[io] = 129.885;  io++;
   fRuns[io] = 287325;  fMeanV0A[io] = 51.9591;  fMeanV0C[io] = 77.063;   fMeanV0M[io] = 128.949;  io++;
   fRuns[io] = 287343;  fMeanV0A[io] = 51.9019;  fMeanV0C[io] = 77.1004;  fMeanV0M[io] = 128.905;  io++;
   fRuns[io] = 287344;  fMeanV0A[io] = 51.665;   fMeanV0C[io] = 76.9331;  fMeanV0M[io] = 128.515;  io++;
   fRuns[io] = 287346;  fMeanV0A[io] = 51.8043;  fMeanV0C[io] = 77.0672;  fMeanV0M[io] = 128.788;  io++;
   fRuns[io] = 287347;  fMeanV0A[io] = 51.7523;  fMeanV0C[io] = 76.9784;  fMeanV0M[io] = 128.641;  io++;
   fRuns[io] = 287349;  fMeanV0A[io] = 51.7251;  fMeanV0C[io] = 76.9843;  fMeanV0M[io] = 128.618;  io++;
   fRuns[io] = 287353;  fMeanV0A[io] = 51.9095;  fMeanV0C[io] = 76.944;   fMeanV0M[io] = 128.765;  io++;
   fRuns[io] = 287355;  fMeanV0A[io] = 51.7266;  fMeanV0C[io] = 76.8758;  fMeanV0M[io] = 128.526;  io++;
   fRuns[io] = 287356;  fMeanV0A[io] = 51.6506;  fMeanV0C[io] = 76.7725;  fMeanV0M[io] = 128.34;  io++;
   fRuns[io] = 287360;  fMeanV0A[io] = 51.7852;  fMeanV0C[io] = 76.8544;  fMeanV0M[io] = 128.566;  io++;
   fRuns[io] = 287380;  fMeanV0A[io] = 52.1157;  fMeanV0C[io] = 77.1391;  fMeanV0M[io] = 129.176;  io++;
   fRuns[io] = 287381;  fMeanV0A[io] = 51.7922;  fMeanV0C[io] = 76.8343;  fMeanV0M[io] = 128.52;  io++;
   fRuns[io] = 287385;  fMeanV0A[io] = 51.6171;  fMeanV0C[io] = 76.5233;  fMeanV0M[io] = 128.059;  io++;
   fRuns[io] = 287387;  fMeanV0A[io] = 51.6648;  fMeanV0C[io] = 76.7807;  fMeanV0M[io] = 128.365;  io++;
   fRuns[io] = 287388;  fMeanV0A[io] = 51.6172;  fMeanV0C[io] = 76.7127;  fMeanV0M[io] = 128.248;  io++;
   fRuns[io] = 287389;  fMeanV0A[io] = 51.5324;  fMeanV0C[io] = 76.6313;  fMeanV0M[io] = 128.078;  io++;
   fRuns[io] = 287413;  fMeanV0A[io] = 51.7286;  fMeanV0C[io] = 76.6695;  fMeanV0M[io] = 128.311;  io++;
   fRuns[io] = 287451;  fMeanV0A[io] = 51.728;   fMeanV0C[io] = 77.0541;  fMeanV0M[io] = 128.704;  io++;
   fRuns[io] = 287480;  fMeanV0A[io] = 51.6599;  fMeanV0C[io] = 76.7909;  fMeanV0M[io] = 128.366;  io++;
   fRuns[io] = 287481;  fMeanV0A[io] = 51.6479;  fMeanV0C[io] = 76.8815;  fMeanV0M[io] = 128.431;  io++;
   fRuns[io] = 287484;  fMeanV0A[io] = 51.6601;  fMeanV0C[io] = 76.7206;  fMeanV0M[io] = 128.292;  io++;
   fRuns[io] = 287486;  fMeanV0A[io] = 51.6674;  fMeanV0C[io] = 76.685;   fMeanV0M[io] = 128.276;  io++;
   fRuns[io] = 287513;  fMeanV0A[io] = 51.8862;  fMeanV0C[io] = 77.0712;  fMeanV0M[io] = 128.867;  io++;
   fRuns[io] = 287516;  fMeanV0A[io] = 51.76;    fMeanV0C[io] = 76.9936;  fMeanV0M[io] = 128.671;  io++;
   fRuns[io] = 287517;  fMeanV0A[io] = 51.9888;  fMeanV0C[io] = 77.0378;  fMeanV0M[io] = 128.945;  io++;
   fRuns[io] = 287518;  fMeanV0A[io] = 51.614;   fMeanV0C[io] = 76.8043;  fMeanV0M[io] = 128.34;  io++;
   fRuns[io] = 287521;  fMeanV0A[io] = 51.6472;  fMeanV0C[io] = 76.9077;  fMeanV0M[io] = 128.462;  io++;
   fRuns[io] = 287524;  fMeanV0A[io] = 51.6282;  fMeanV0C[io] = 76.7354;  fMeanV0M[io] = 128.273;  io++;
   fRuns[io] = 287575;  fMeanV0A[io] = 51.8064;  fMeanV0C[io] = 76.9062;  fMeanV0M[io] = 128.639;  io++;
   fRuns[io] = 287578;  fMeanV0A[io] = 51.5895;  fMeanV0C[io] = 76.7064;  fMeanV0M[io] = 128.221;  io++;
   fRuns[io] = 287654;  fMeanV0A[io] = 51.6042;  fMeanV0C[io] = 76.6736;  fMeanV0M[io] = 128.202;  io++;
   fRuns[io] = 287656;  fMeanV0A[io] = 51.429;   fMeanV0C[io] = 76.5177;  fMeanV0M[io] = 127.855;  io++;
   fRuns[io] = 287657;  fMeanV0A[io] = 51.4698;  fMeanV0C[io] = 76.7682;  fMeanV0M[io] = 128.155;  io++;
   fRuns[io] = 287658;  fMeanV0A[io] = 51.3672;  fMeanV0C[io] = 76.7381;  fMeanV0M[io] = 128.011;  io++;

   // FILTER_p-p_208_LHC18g
   fRuns[io] = 288619;  fMeanV0A[io] = 52.5386;  fMeanV0C[io] = 78.3106;  fMeanV0M[io] = 130.766;  io++;
   fRuns[io] = 288640;  fMeanV0A[io] = 52.3432;  fMeanV0C[io] = 77.852;   fMeanV0M[io] = 130.112;  io++;
   fRuns[io] = 288642;  fMeanV0A[io] = 52.1074;  fMeanV0C[io] = 77.6374;  fMeanV0M[io] = 129.657;  io++;
   fRuns[io] = 288644;  fMeanV0A[io] = 51.9189;  fMeanV0C[io] = 77.5177;  fMeanV0M[io] = 129.359;  io++;
   fRuns[io] = 288650;  fMeanV0A[io] = 51.8638;  fMeanV0C[io] = 77.3969;  fMeanV0M[io] = 129.181;  io++;
   fRuns[io] = 288687;  fMeanV0A[io] = 52.1;     fMeanV0C[io] = 77.2917;  fMeanV0M[io] = 129.306;  io++;
   fRuns[io] = 288689;  fMeanV0A[io] = 52.1288;  fMeanV0C[io] = 77.2345;  fMeanV0M[io] = 129.298;  io++;
   fRuns[io] = 288690;  fMeanV0A[io] = 52.1314;  fMeanV0C[io] = 77.3761;  fMeanV0M[io] = 129.426;  io++;
   fRuns[io] = 288743;  fMeanV0A[io] = 52.2469;  fMeanV0C[io] = 77.7047;  fMeanV0M[io] = 129.857;  io++;
   fRuns[io] = 288748;  fMeanV0A[io] = 52.3296;  fMeanV0C[io] = 78.023;   fMeanV0M[io] = 130.27;  io++;
   fRuns[io] = 288750;  fMeanV0A[io] = 52.2386;  fMeanV0C[io] = 77.9058;  fMeanV0M[io] = 130.071;  io++;


   // FILTER_p-p_208_LHC18h
   fRuns[io] = 288804;  fMeanV0A[io] = 51.8602;  fMeanV0C[io] = 77.2465;  fMeanV0M[io] = 129.014;  io++;
   fRuns[io] = 288806;  fMeanV0A[io] = 51.7138;  fMeanV0C[io] = 77.1213;  fMeanV0M[io] = 128.754;  io++;


   // FILTER_p-p_208_LHC18i
   fRuns[io] = 288861;  fMeanV0A[io] = 53.3457;  fMeanV0C[io] = 78.8432;  fMeanV0M[io] = 132.12;  io++;
   fRuns[io] = 288862;  fMeanV0A[io] = 53.038;   fMeanV0C[io] = 78.5457;  fMeanV0M[io] = 131.506;  io++;
   fRuns[io] = 288863;  fMeanV0A[io] = 52.3928;  fMeanV0C[io] = 77.3812;  fMeanV0M[io] = 129.729;  io++;
   fRuns[io] = 288864;  fMeanV0A[io] = 52.623;   fMeanV0C[io] = 78.1954;  fMeanV0M[io] = 130.734;  io++;
   fRuns[io] = 288868;  fMeanV0A[io] = 52.11;    fMeanV0C[io] = 77.6959;  fMeanV0M[io] = 129.707;  io++;
   fRuns[io] = 288902;  fMeanV0A[io] = 52.3148;  fMeanV0C[io] = 77.9176;  fMeanV0M[io] = 130.154;  io++;
   fRuns[io] = 288903;  fMeanV0A[io] = 52.272;   fMeanV0C[io] = 78.07;    fMeanV0M[io] = 130.259;  io++;
   fRuns[io] = 288908;  fMeanV0A[io] = 52.1213;  fMeanV0C[io] = 77.4056;  fMeanV0M[io] = 129.433;  io++;
   fRuns[io] = 288909;  fMeanV0A[io] = 52.1628;  fMeanV0C[io] = 77.8335;  fMeanV0M[io] = 129.911;  io++;


   // FILTER_p-p_208_LHC18j
   fRuns[io] = 288943;  fMeanV0A[io] = 51.5825;  fMeanV0C[io] = 77.714;   fMeanV0M[io] = 129.203;  io++;

   // FILTER_p-p_208_LHC18k
   fRuns[io] = 289165;  fMeanV0A[io] = 51.7402;  fMeanV0C[io] = 77.182;   fMeanV0M[io] = 128.832;  io++;
   fRuns[io] = 289166;  fMeanV0A[io] = 51.8466;  fMeanV0C[io] = 77.13;    fMeanV0M[io] = 128.884;  io++;
   fRuns[io] = 289167;  fMeanV0A[io] = 51.8459;  fMeanV0C[io] = 77.5012;  fMeanV0M[io] = 129.265;  io++;
   fRuns[io] = 289169;  fMeanV0A[io] = 51.9775;  fMeanV0C[io] = 77.5281;  fMeanV0M[io] = 129.428;  io++;
   fRuns[io] = 289172;  fMeanV0A[io] = 52.1654;  fMeanV0C[io] = 77.7965;  fMeanV0M[io] = 129.877;  io++;
   fRuns[io] = 289175;  fMeanV0A[io] = 51.9446;  fMeanV0C[io] = 77.4985;  fMeanV0M[io] = 129.359;  io++;
   fRuns[io] = 289176;  fMeanV0A[io] = 51.9435;  fMeanV0C[io] = 77.4534;  fMeanV0M[io] = 129.326;  io++;
   fRuns[io] = 289177;  fMeanV0A[io] = 52.0222;  fMeanV0C[io] = 77.5201;  fMeanV0M[io] = 129.457;  io++;
   fRuns[io] = 289198;  fMeanV0A[io] = 52.358;   fMeanV0C[io] = 77.7195;  fMeanV0M[io] = 129.986;  io++;
   fRuns[io] = 289199;  fMeanV0A[io] = 52.0842;  fMeanV0C[io] = 77.8183;  fMeanV0M[io] = 129.818;  io++;
   fRuns[io] = 289200;  fMeanV0A[io] = 51.9986;  fMeanV0C[io] = 77.6431;  fMeanV0M[io] = 129.558;  io++;
   fRuns[io] = 289201;  fMeanV0A[io] = 52.061;   fMeanV0C[io] = 77.5992;  fMeanV0M[io] = 129.585;  io++;


   // FILTER_p-p_208_LHC18l
   fRuns[io] = 289240;  fMeanV0A[io] = 52.2585;  fMeanV0C[io] = 77.5607;  fMeanV0M[io] = 129.716;  io++;
   fRuns[io] = 289241;  fMeanV0A[io] = 52.0142;  fMeanV0C[io] = 77.4558;  fMeanV0M[io] = 129.389;  io++;
   fRuns[io] = 289242;  fMeanV0A[io] = 52.0008;  fMeanV0C[io] = 77.3972;  fMeanV0M[io] = 129.305;  io++;
   fRuns[io] = 289243;  fMeanV0A[io] = 52.0388;  fMeanV0C[io] = 77.3875;  fMeanV0M[io] = 129.347;  io++;
   fRuns[io] = 289247;  fMeanV0A[io] = 52.0264;  fMeanV0C[io] = 77.2262;  fMeanV0M[io] = 129.171;  io++;
   fRuns[io] = 289249;  fMeanV0A[io] = 51.9727;  fMeanV0C[io] = 77.322;   fMeanV0M[io] = 129.215;  io++;
   fRuns[io] = 289253;  fMeanV0A[io] = 52.0342;  fMeanV0C[io] = 77.1517;  fMeanV0M[io] = 129.101;  io++;
   fRuns[io] = 289254;  fMeanV0A[io] = 51.8127;  fMeanV0C[io] = 76.8771;  fMeanV0M[io] = 128.605;  io++;
   fRuns[io] = 289275;  fMeanV0A[io] = 51.9516;  fMeanV0C[io] = 77.1479;  fMeanV0M[io] = 129.025;  io++;
   fRuns[io] = 289276;  fMeanV0A[io] = 51.9228;  fMeanV0C[io] = 76.9921;  fMeanV0M[io] = 128.816;  io++;
   fRuns[io] = 289277;  fMeanV0A[io] = 51.79;    fMeanV0C[io] = 76.8182;  fMeanV0M[io] = 128.523;  io++;
   fRuns[io] = 289278;  fMeanV0A[io] = 51.9554;  fMeanV0C[io] = 76.7349;  fMeanV0M[io] = 128.595;  io++;
   fRuns[io] = 289280;  fMeanV0A[io] = 51.9235;  fMeanV0C[io] = 77.0264;  fMeanV0M[io] = 128.877;  io++;
   fRuns[io] = 289281;  fMeanV0A[io] = 52.1044;  fMeanV0C[io] = 76.9804;  fMeanV0M[io] = 128.995;  io++;
   fRuns[io] = 289300;  fMeanV0A[io] = 52.1588;  fMeanV0C[io] = 77.159;   fMeanV0M[io] = 129.221;  io++;
   fRuns[io] = 289303;  fMeanV0A[io] = 51.8904;  fMeanV0C[io] = 76.9738;  fMeanV0M[io] = 128.778;  io++;
   fRuns[io] = 289306;  fMeanV0A[io] = 51.7409;  fMeanV0C[io] = 76.9081;  fMeanV0M[io] = 128.577;  io++;
   fRuns[io] = 289308;  fMeanV0A[io] = 51.7041;  fMeanV0C[io] = 76.7773;  fMeanV0M[io] = 128.396;  io++;
   fRuns[io] = 289309;  fMeanV0A[io] = 51.7146;  fMeanV0C[io] = 76.5458;  fMeanV0M[io] = 128.173;  io++;
   fRuns[io] = 289353;  fMeanV0A[io] = 51.666;   fMeanV0C[io] = 76.7771;  fMeanV0M[io] = 128.355;  io++;
   fRuns[io] = 289354;  fMeanV0A[io] = 51.7561;  fMeanV0C[io] = 76.7093;  fMeanV0M[io] = 128.389;  io++;
   fRuns[io] = 289355;  fMeanV0A[io] = 51.7672;  fMeanV0C[io] = 76.8363;  fMeanV0M[io] = 128.515;  io++;
   fRuns[io] = 289356;  fMeanV0A[io] = 51.7034;  fMeanV0C[io] = 76.7109;  fMeanV0M[io] = 128.327;  io++;
   fRuns[io] = 289365;  fMeanV0A[io] = 51.6219;  fMeanV0C[io] = 76.6483;  fMeanV0M[io] = 128.166;  io++;
   fRuns[io] = 289366;  fMeanV0A[io] = 51.6518;  fMeanV0C[io] = 76.7704;  fMeanV0M[io] = 128.344;  io++;
   fRuns[io] = 289367;  fMeanV0A[io] = 51.6724;  fMeanV0C[io] = 76.8135;  fMeanV0M[io] = 128.405;  io++;
   fRuns[io] = 289368;  fMeanV0A[io] = 51.7612;  fMeanV0C[io] = 76.7278;  fMeanV0M[io] = 128.398;  io++;
   fRuns[io] = 289369;  fMeanV0A[io] = 51.6089;  fMeanV0C[io] = 76.5883;  fMeanV0M[io] = 128.128;  io++;
   fRuns[io] = 289370;  fMeanV0A[io] = 51.5653;  fMeanV0C[io] = 76.4295;  fMeanV0M[io] = 127.884;  io++;
   fRuns[io] = 289373;  fMeanV0A[io] = 51.5303;  fMeanV0C[io] = 76.4517;  fMeanV0M[io] = 127.897;  io++;
   fRuns[io] = 289374;  fMeanV0A[io] = 51.6396;  fMeanV0C[io] = 76.4563;  fMeanV0M[io] = 128.017;  io++;
   fRuns[io] = 289426;  fMeanV0A[io] = 51.9574;  fMeanV0C[io] = 76.9491;  fMeanV0M[io] = 128.819;  io++;
   fRuns[io] = 289444;  fMeanV0A[io] = 51.6149;  fMeanV0C[io] = 76.6977;  fMeanV0M[io] = 128.225;  io++;
   fRuns[io] = 289462;  fMeanV0A[io] = 51.6037;  fMeanV0C[io] = 76.6411;  fMeanV0M[io] = 128.15;  io++;
   fRuns[io] = 289463;  fMeanV0A[io] = 51.6381;  fMeanV0C[io] = 76.6397;  fMeanV0M[io] = 128.195;  io++;
   fRuns[io] = 289465;  fMeanV0A[io] = 51.5987;  fMeanV0C[io] = 76.746;   fMeanV0M[io] = 128.257;  io++;
   fRuns[io] = 289466;  fMeanV0A[io] = 51.5798;  fMeanV0C[io] = 76.5579;  fMeanV0M[io] = 128.049;  io++;
   fRuns[io] = 289468;  fMeanV0A[io] = 51.4764;  fMeanV0C[io] = 76.846;   fMeanV0M[io] = 128.248;  io++;
   fRuns[io] = 289493;  fMeanV0A[io] = 51.924;   fMeanV0C[io] = 77.0113;  fMeanV0M[io] = 128.855;  io++;
   fRuns[io] = 289494;  fMeanV0A[io] = 51.7669;  fMeanV0C[io] = 76.8876;  fMeanV0M[io] = 128.557;  io++;
   fRuns[io] = 289521;  fMeanV0A[io] = 52.1248;  fMeanV0C[io] = 76.9479;  fMeanV0M[io] = 129.007;  io++;
   fRuns[io] = 289547;  fMeanV0A[io] = 51.7164;  fMeanV0C[io] = 76.8402;  fMeanV0M[io] = 128.471;  io++;
   fRuns[io] = 289721;  fMeanV0A[io] = 51.7513;  fMeanV0C[io] = 76.7275;  fMeanV0M[io] = 128.398;  io++;
   fRuns[io] = 289723;  fMeanV0A[io] = 51.5563;  fMeanV0C[io] = 76.4216;  fMeanV0M[io] = 127.889;  io++;
   fRuns[io] = 289724;  fMeanV0A[io] = 51.4702;  fMeanV0C[io] = 76.4743;  fMeanV0M[io] = 127.856;  io++;
   fRuns[io] = 289729;  fMeanV0A[io] = 51.5103;  fMeanV0C[io] = 76.6061;  fMeanV0M[io] = 128.033;  io++;
   fRuns[io] = 289731;  fMeanV0A[io] = 51.5567;  fMeanV0C[io] = 76.7395;  fMeanV0M[io] = 128.228;  io++;
   fRuns[io] = 289732;  fMeanV0A[io] = 51.467;   fMeanV0C[io] = 76.4716;  fMeanV0M[io] = 127.855;  io++;
   fRuns[io] = 289757;  fMeanV0A[io] = 51.7271;  fMeanV0C[io] = 76.8151;  fMeanV0M[io] = 128.449;  io++;
   fRuns[io] = 289775;  fMeanV0A[io] = 51.7824;  fMeanV0C[io] = 76.8697;  fMeanV0M[io] = 128.566;  io++;
   fRuns[io] = 289808;  fMeanV0A[io] = 51.6899;  fMeanV0C[io] = 77.0409;  fMeanV0M[io] = 128.652;  io++;
   fRuns[io] = 289811;  fMeanV0A[io] = 51.3624;  fMeanV0C[io] = 76.4141;  fMeanV0M[io] = 127.702;  io++;
   fRuns[io] = 289814;  fMeanV0A[io] = 51.3224;  fMeanV0C[io] = 76.4159;  fMeanV0M[io] = 127.64;  io++;
   fRuns[io] = 289815;  fMeanV0A[io] = 51.3708;  fMeanV0C[io] = 76.5202;  fMeanV0M[io] = 127.807;  io++;
   fRuns[io] = 289816;  fMeanV0A[io] = 51.645;   fMeanV0C[io] = 76.3626;  fMeanV0M[io] = 127.887;  io++;
   fRuns[io] = 289817;  fMeanV0A[io] = 51.5962;  fMeanV0C[io] = 76.7998;  fMeanV0M[io] = 128.3;  io++;
   fRuns[io] = 289818;  fMeanV0A[io] = 51.2436;  fMeanV0C[io] = 76.4734;  fMeanV0M[io] = 127.637;  io++;
   fRuns[io] = 289830;  fMeanV0A[io] = 51.3208;  fMeanV0C[io] = 76.5179;  fMeanV0M[io] = 127.748;  io++;
   fRuns[io] = 289849;  fMeanV0A[io] = 51.7393;  fMeanV0C[io] = 76.8408;  fMeanV0M[io] = 128.494;  io++;
   fRuns[io] = 289852;  fMeanV0A[io] = 51.4999;  fMeanV0C[io] = 76.5313;  fMeanV0M[io] = 127.943;  io++;
   fRuns[io] = 289854;  fMeanV0A[io] = 51.4219;  fMeanV0C[io] = 76.6347;  fMeanV0M[io] = 127.985;  io++;
   fRuns[io] = 289855;  fMeanV0A[io] = 51.3322;  fMeanV0C[io] = 76.0443;  fMeanV0M[io] = 127.282;  io++;
   fRuns[io] = 289856;  fMeanV0A[io] = 51.396;   fMeanV0C[io] = 76.5514;  fMeanV0M[io] = 127.868;  io++;
   fRuns[io] = 289857;  fMeanV0A[io] = 51.3949;  fMeanV0C[io] = 76.222;   fMeanV0M[io] = 127.534;  io++;
   fRuns[io] = 289879;  fMeanV0A[io] = 51.17;    fMeanV0C[io] = 76.2705;  fMeanV0M[io] = 127.38;  io++;
   fRuns[io] = 289880;  fMeanV0A[io] = 51.6607;  fMeanV0C[io] = 76.6305;  fMeanV0M[io] = 128.19;  io++;
   fRuns[io] = 289884;  fMeanV0A[io] = 51.5114;  fMeanV0C[io] = 76.3406;  fMeanV0M[io] = 127.761;  io++;
   fRuns[io] = 289928;  fMeanV0A[io] = 51.0937;  fMeanV0C[io] = 76.0376;  fMeanV0M[io] = 127.041;  io++;
   fRuns[io] = 289931;  fMeanV0A[io] = 51.48;    fMeanV0C[io] = 76.5259;  fMeanV0M[io] = 127.917;  io++;
   fRuns[io] = 289935;  fMeanV0A[io] = 51.3228;  fMeanV0C[io] = 76.4237;  fMeanV0M[io] = 127.673;  io++;
   fRuns[io] = 289940;  fMeanV0A[io] = 51.4163;  fMeanV0C[io] = 76.5142;  fMeanV0M[io] = 127.835;  io++;
   fRuns[io] = 289941;  fMeanV0A[io] = 51.3622;  fMeanV0C[io] = 76.1488;  fMeanV0M[io] = 127.417;  io++;
   fRuns[io] = 289943;  fMeanV0A[io] = 51.3023;  fMeanV0C[io] = 76.1292;  fMeanV0M[io] = 127.337;  io++;
   fRuns[io] = 289965;  fMeanV0A[io] = 51.4256;  fMeanV0C[io] = 76.6159;  fMeanV0M[io] = 127.962;  io++;
   fRuns[io] = 289966;  fMeanV0A[io] = 51.3992;  fMeanV0C[io] = 76.205;   fMeanV0M[io] = 127.518;  io++;
   fRuns[io] = 289971;  fMeanV0A[io] = 51.1295;  fMeanV0C[io] = 76.0855;  fMeanV0M[io] = 127.132;  io++;


   // FILTER_p-p_208_LHC18m
   fRuns[io] = 290323;  fMeanV0A[io] = 52.2064;  fMeanV0C[io] = 77.3633;  fMeanV0M[io] = 129.486;  io++;
   fRuns[io] = 290327;  fMeanV0A[io] = 51.9285;  fMeanV0C[io] = 76.7395;  fMeanV0M[io] = 128.605;  io++;
   fRuns[io] = 290350;  fMeanV0A[io] = 52.3693;  fMeanV0C[io] = 77.4215;  fMeanV0M[io] = 129.702;  io++;
   fRuns[io] = 290374;  fMeanV0A[io] = 52.3645;  fMeanV0C[io] = 77.295;   fMeanV0M[io] = 129.576;  io++;
   fRuns[io] = 290375;  fMeanV0A[io] = 52.2277;  fMeanV0C[io] = 77.3294;  fMeanV0M[io] = 129.474;  io++;
   fRuns[io] = 290376;  fMeanV0A[io] = 52.4972;  fMeanV0C[io] = 77.7087;  fMeanV0M[io] = 130.15;  io++;
   fRuns[io] = 290399;  fMeanV0A[io] = 52.7293;  fMeanV0C[io] = 77.8846;  fMeanV0M[io] = 130.535;  io++;
   fRuns[io] = 290401;  fMeanV0A[io] = 52.3672;  fMeanV0C[io] = 77.5119;  fMeanV0M[io] = 129.776;  io++;
   fRuns[io] = 290404;  fMeanV0A[io] = 52.274;   fMeanV0C[io] = 77.4237;  fMeanV0M[io] = 129.617;  io++;
   fRuns[io] = 290411;  fMeanV0A[io] = 52.4419;  fMeanV0C[io] = 77.5173;  fMeanV0M[io] = 129.868;  io++;
   fRuns[io] = 290412;  fMeanV0A[io] = 52.274;   fMeanV0C[io] = 77.4275;  fMeanV0M[io] = 129.622;  io++;
   fRuns[io] = 290423;  fMeanV0A[io] = 52.2335;  fMeanV0C[io] = 77.296;   fMeanV0M[io] = 129.464;  io++;
   fRuns[io] = 290425;  fMeanV0A[io] = 52.1335;  fMeanV0C[io] = 77.194;   fMeanV0M[io] = 129.232;  io++;
   fRuns[io] = 290426;  fMeanV0A[io] = 51.8635;  fMeanV0C[io] = 76.9449;  fMeanV0M[io] = 128.678;  io++;
   fRuns[io] = 290427;  fMeanV0A[io] = 52.1392;  fMeanV0C[io] = 77.2513;  fMeanV0M[io] = 129.306;  io++;
   fRuns[io] = 290456;  fMeanV0A[io] = 51.9943;  fMeanV0C[io] = 76.9586;  fMeanV0M[io] = 128.866;  io++;
   fRuns[io] = 290458;  fMeanV0A[io] = 51.9844;  fMeanV0C[io] = 76.8064;  fMeanV0M[io] = 128.697;  io++;
   fRuns[io] = 290459;  fMeanV0A[io] = 52.0189;  fMeanV0C[io] = 77.0551;  fMeanV0M[io] = 128.989;  io++;
   fRuns[io] = 290467;  fMeanV0A[io] = 51.9952;  fMeanV0C[io] = 76.7881;  fMeanV0M[io] = 128.695;  io++;
   fRuns[io] = 290469;  fMeanV0A[io] = 51.8914;  fMeanV0C[io] = 76.8597;  fMeanV0M[io] = 128.67;  io++;
   fRuns[io] = 290499;  fMeanV0A[io] = 51.9294;  fMeanV0C[io] = 76.9302;  fMeanV0M[io] = 128.781;  io++;
   fRuns[io] = 290500;  fMeanV0A[io] = 51.9248;  fMeanV0C[io] = 77.3379;  fMeanV0M[io] = 129.253;  io++;
   fRuns[io] = 290501;  fMeanV0A[io] = 51.9219;  fMeanV0C[io] = 76.8495;  fMeanV0M[io] = 128.688;  io++;
   fRuns[io] = 290538;  fMeanV0A[io] = 52.0808;  fMeanV0C[io] = 77.0395;  fMeanV0M[io] = 129.032;  io++;
   fRuns[io] = 290539;  fMeanV0A[io] = 51.9219;  fMeanV0C[io] = 77.011;   fMeanV0M[io] = 128.844;  io++;
   fRuns[io] = 290540;  fMeanV0A[io] = 52.0272;  fMeanV0C[io] = 76.8824;  fMeanV0M[io] = 128.834;  io++;
   fRuns[io] = 290544;  fMeanV0A[io] = 51.7968;  fMeanV0C[io] = 77.0168;  fMeanV0M[io] = 128.728;  io++;
   fRuns[io] = 290549;  fMeanV0A[io] = 51.8156;  fMeanV0C[io] = 76.7858;  fMeanV0M[io] = 128.511;  io++;
   fRuns[io] = 290550;  fMeanV0A[io] = 51.7506;  fMeanV0C[io] = 76.5703;  fMeanV0M[io] = 128.233;  io++;
   fRuns[io] = 290553;  fMeanV0A[io] = 51.8604;  fMeanV0C[io] = 76.6603;  fMeanV0M[io] = 128.444;  io++;
   fRuns[io] = 290588;  fMeanV0A[io] = 52.2537;  fMeanV0C[io] = 77.5131;  fMeanV0M[io] = 129.665;  io++;
   fRuns[io] = 290590;  fMeanV0A[io] = 51.8789;  fMeanV0C[io] = 76.7163;  fMeanV0M[io] = 128.516;  io++;
   fRuns[io] = 290612;  fMeanV0A[io] = 52.3618;  fMeanV0C[io] = 77.3577;  fMeanV0M[io] = 129.652;  io++;
   fRuns[io] = 290613;  fMeanV0A[io] = 52.0724;  fMeanV0C[io] = 77.2287;  fMeanV0M[io] = 129.224;  io++;
   fRuns[io] = 290614;  fMeanV0A[io] = 51.9877;  fMeanV0C[io] = 76.6384;  fMeanV0M[io] = 128.537;  io++;
   fRuns[io] = 290615;  fMeanV0A[io] = 52.2762;  fMeanV0C[io] = 77.2987;  fMeanV0M[io] = 129.479;  io++;
   fRuns[io] = 290627;  fMeanV0A[io] = 51.8114;  fMeanV0C[io] = 76.7725;  fMeanV0M[io] = 128.496;  io++;
   fRuns[io] = 290632;  fMeanV0A[io] = 51.8287;  fMeanV0C[io] = 76.8006;  fMeanV0M[io] = 128.539;  io++;
   fRuns[io] = 290645;  fMeanV0A[io] = 51.6914;  fMeanV0C[io] = 76.534;   fMeanV0M[io] = 128.139;  io++;
   fRuns[io] = 290660;  fMeanV0A[io] = 51.6548;  fMeanV0C[io] = 76.4602;  fMeanV0M[io] = 128.032;  io++;
   fRuns[io] = 290665;  fMeanV0A[io] = 51.7612;  fMeanV0C[io] = 76.5901;  fMeanV0M[io] = 128.267;  io++;
   fRuns[io] = 290687;  fMeanV0A[io] = 51.7082;  fMeanV0C[io] = 76.5158;  fMeanV0M[io] = 128.147;  io++;
   fRuns[io] = 290689;  fMeanV0A[io] = 51.6002;  fMeanV0C[io] = 76.9167;  fMeanV0M[io] = 128.466;  io++;
   fRuns[io] = 290766;  fMeanV0A[io] = 51.6419;  fMeanV0C[io] = 76.4301;  fMeanV0M[io] = 127.984;  io++;
   fRuns[io] = 290787;  fMeanV0A[io] = 51.5211;  fMeanV0C[io] = 75.9733;  fMeanV0M[io] = 127.409;  io++;
   fRuns[io] = 290790;  fMeanV0A[io] = 51.6447;  fMeanV0C[io] = 75.9601;  fMeanV0M[io] = 127.514;  io++;
   fRuns[io] = 290841;  fMeanV0A[io] = 51.6479;  fMeanV0C[io] = 76.557;   fMeanV0M[io] = 128.133;  io++;
   fRuns[io] = 290843;  fMeanV0A[io] = 51.2746;  fMeanV0C[io] = 76.0943;  fMeanV0M[io] = 127.265;  io++;
   fRuns[io] = 290846;  fMeanV0A[io] = 51.4794;  fMeanV0C[io] = 76.1097;  fMeanV0M[io] = 127.505;  io++;
   fRuns[io] = 290848;  fMeanV0A[io] = 51.5275;  fMeanV0C[io] = 76.2839;  fMeanV0M[io] = 127.726;  io++;
   fRuns[io] = 290853;  fMeanV0A[io] = 51.468;   fMeanV0C[io] = 76.3142;  fMeanV0M[io] = 127.688;  io++;
   fRuns[io] = 290860;  fMeanV0A[io] = 51.4459;  fMeanV0C[io] = 76.1283;  fMeanV0M[io] = 127.489;  io++;
   fRuns[io] = 290862;  fMeanV0A[io] = 51.3741;  fMeanV0C[io] = 76.1945;  fMeanV0M[io] = 127.465;  io++;
   fRuns[io] = 290886;  fMeanV0A[io] = 51.5112;  fMeanV0C[io] = 76.109;   fMeanV0M[io] = 127.536;  io++;
   fRuns[io] = 290887;  fMeanV0A[io] = 51.2043;  fMeanV0C[io] = 75.9277;  fMeanV0M[io] = 127.052;  io++;
   fRuns[io] = 290888;  fMeanV0A[io] = 51.3418;  fMeanV0C[io] = 76.129;   fMeanV0M[io] = 127.398;  io++;
   fRuns[io] = 290894;  fMeanV0A[io] = 51.0142;  fMeanV0C[io] = 75.6584;  fMeanV0M[io] = 126.585;  io++;
   fRuns[io] = 290895;  fMeanV0A[io] = 51.0335;  fMeanV0C[io] = 75.8954;  fMeanV0M[io] = 126.854;  io++;
   fRuns[io] = 290932;  fMeanV0A[io] = 51.6141;  fMeanV0C[io] = 76.233;   fMeanV0M[io] = 127.752;  io++;
   fRuns[io] = 290935;  fMeanV0A[io] = 51.2697;  fMeanV0C[io] = 76.1677;  fMeanV0M[io] = 127.36;  io++;
   fRuns[io] = 290941;  fMeanV0A[io] = 51.3983;  fMeanV0C[io] = 76.0139;  fMeanV0M[io] = 127.319;  io++;
   fRuns[io] = 290943;  fMeanV0A[io] = 51.2535;  fMeanV0C[io] = 75.8574;  fMeanV0M[io] = 127.028;  io++;
   fRuns[io] = 290944;  fMeanV0A[io] = 51.3412;  fMeanV0C[io] = 75.859;   fMeanV0M[io] = 127.107;  io++;
   fRuns[io] = 290948;  fMeanV0A[io] = 51.2248;  fMeanV0C[io] = 75.6629;  fMeanV0M[io] = 126.808;  io++;
   fRuns[io] = 290974;  fMeanV0A[io] = 51.461;   fMeanV0C[io] = 76.1595;  fMeanV0M[io] = 127.541;  io++;
   fRuns[io] = 290975;  fMeanV0A[io] = 51.4044;  fMeanV0C[io] = 75.8265;  fMeanV0M[io] = 127.15;  io++;
   fRuns[io] = 290976;  fMeanV0A[io] = 51.2275;  fMeanV0C[io] = 75.8258;  fMeanV0M[io] = 126.977;  io++;
   fRuns[io] = 290979;  fMeanV0A[io] = 51.2649;  fMeanV0C[io] = 75.7462;  fMeanV0M[io] = 126.932;  io++;
   fRuns[io] = 290980;  fMeanV0A[io] = 51.0679;  fMeanV0C[io] = 75.6935;  fMeanV0M[io] = 126.663;  io++;
   fRuns[io] = 291002;  fMeanV0A[io] = 51.1839;  fMeanV0C[io] = 75.707;   fMeanV0M[io] = 126.802;  io++;
   fRuns[io] = 291003;  fMeanV0A[io] = 51.0936;  fMeanV0C[io] = 75.6163;  fMeanV0M[io] = 126.622;  io++;
   fRuns[io] = 291004;  fMeanV0A[io] = 51.0627;  fMeanV0C[io] = 75.5463;  fMeanV0M[io] = 126.524;  io++;
   fRuns[io] = 291005;  fMeanV0A[io] = 51.1333;  fMeanV0C[io] = 75.6028;  fMeanV0M[io] = 126.652;  io++;
   fRuns[io] = 291006;  fMeanV0A[io] = 51.2903;  fMeanV0C[io] = 75.8008;  fMeanV0M[io] = 127.002;  io++;
   fRuns[io] = 291035;  fMeanV0A[io] = 51.3027;  fMeanV0C[io] = 75.7608;  fMeanV0M[io] = 126.984;  io++;
   fRuns[io] = 291037;  fMeanV0A[io] = 51.0086;  fMeanV0C[io] = 75.8098;  fMeanV0M[io] = 126.746;  io++;
   fRuns[io] = 291041;  fMeanV0A[io] = 51.0047;  fMeanV0C[io] = 75.496;   fMeanV0M[io] = 126.414;  io++;
   fRuns[io] = 291065;  fMeanV0A[io] = 51.3357;  fMeanV0C[io] = 76.1275;  fMeanV0M[io] = 127.402;  io++;
   fRuns[io] = 291066;  fMeanV0A[io] = 51.0591;  fMeanV0C[io] = 75.5136;  fMeanV0M[io] = 126.487;  io++;
   fRuns[io] = 291069;  fMeanV0A[io] = 51.0437;  fMeanV0C[io] = 75.2797;  fMeanV0M[io] = 126.225;  io++;
   fRuns[io] = 291093;  fMeanV0A[io] = 51.036;   fMeanV0C[io] = 75.7405;  fMeanV0M[io] = 126.699;  io++;
   fRuns[io] = 291100;  fMeanV0A[io] = 51.0369;  fMeanV0C[io] = 75.6177;  fMeanV0M[io] = 126.557;  io++;
   fRuns[io] = 291101;  fMeanV0A[io] = 51.073;   fMeanV0C[io] = 75.6854;  fMeanV0M[io] = 126.675;  io++;
   fRuns[io] = 291110;  fMeanV0A[io] = 50.8681;  fMeanV0C[io] = 75.4631;  fMeanV0M[io] = 126.255;  io++;
   fRuns[io] = 291111;  fMeanV0A[io] = 50.8683;  fMeanV0C[io] = 75.5417;  fMeanV0M[io] = 126.333;  io++;
   fRuns[io] = 291116;  fMeanV0A[io] = 50.9629;  fMeanV0C[io] = 75.5657;  fMeanV0M[io] = 126.445;  io++;
   fRuns[io] = 291143;  fMeanV0A[io] = 51.1769;  fMeanV0C[io] = 75.7062;  fMeanV0M[io] = 126.791;  io++;
   fRuns[io] = 291188;  fMeanV0A[io] = 51.7658;  fMeanV0C[io] = 76.3079;  fMeanV0M[io] = 127.982;  io++;
   fRuns[io] = 291209;  fMeanV0A[io] = 50.9873;  fMeanV0C[io] = 75.4563;  fMeanV0M[io] = 126.362;  io++;
   fRuns[io] = 291240;  fMeanV0A[io] = 50.8475;  fMeanV0C[io] = 75.4179;  fMeanV0M[io] = 126.184;  io++;
   fRuns[io] = 291257;  fMeanV0A[io] = 51.1024;  fMeanV0C[io] = 75.3795;  fMeanV0M[io] = 126.4;  io++;
   fRuns[io] = 291262;  fMeanV0A[io] = 50.5052;  fMeanV0C[io] = 74.1741;  fMeanV0M[io] = 124.559;  io++;
   fRuns[io] = 291263;  fMeanV0A[io] = 50.8428;  fMeanV0C[io] = 75.0773;  fMeanV0M[io] = 125.833;  io++;
   fRuns[io] = 291265;  fMeanV0A[io] = 50.9601;  fMeanV0C[io] = 75.1068;  fMeanV0M[io] = 125.981;  io++;
   fRuns[io] = 291266;  fMeanV0A[io] = 50.8517;  fMeanV0C[io] = 75.0087;  fMeanV0M[io] = 125.77;  io++;
   fRuns[io] = 291282;  fMeanV0A[io] = 50.7607;  fMeanV0C[io] = 75.1423;  fMeanV0M[io] = 125.835;  io++;
   fRuns[io] = 291284;  fMeanV0A[io] = 50.7404;  fMeanV0C[io] = 75.0348;  fMeanV0M[io] = 125.686;  io++;
   fRuns[io] = 291285;  fMeanV0A[io] = 50.9093;  fMeanV0C[io] = 75.182;   fMeanV0M[io] = 126.001;  io++;
   fRuns[io] = 291286;  fMeanV0A[io] = 50.7714;  fMeanV0C[io] = 75.1617;  fMeanV0M[io] = 125.845;  io++;
   fRuns[io] = 291360;  fMeanV0A[io] = 50.8601;  fMeanV0C[io] = 75.2204;  fMeanV0M[io] = 125.991;  io++;
   fRuns[io] = 291361;  fMeanV0A[io] = 50.745;   fMeanV0C[io] = 74.9994;  fMeanV0M[io] = 125.66;  io++;
   fRuns[io] = 291362;  fMeanV0A[io] = 50.7276;  fMeanV0C[io] = 75.0378;  fMeanV0M[io] = 125.681;  io++;
   fRuns[io] = 291363;  fMeanV0A[io] = 50.7185;  fMeanV0C[io] = 75.0882;  fMeanV0M[io] = 125.719;  io++;
   fRuns[io] = 291375;  fMeanV0A[io] = 50.785;   fMeanV0C[io] = 75.0806;  fMeanV0M[io] = 125.78;  io++;
   fRuns[io] = 291377;  fMeanV0A[io] = 50.7999;  fMeanV0C[io] = 75.032;   fMeanV0M[io] = 125.747;  io++;
   fRuns[io] = 291397;  fMeanV0A[io] = 50.616;   fMeanV0C[io] = 74.914;   fMeanV0M[io] = 125.436;  io++;
   fRuns[io] = 291399;  fMeanV0A[io] = 50.7689;  fMeanV0C[io] = 75.0046;  fMeanV0M[io] = 125.685;  io++;
   fRuns[io] = 291400;  fMeanV0A[io] = 50.8546;  fMeanV0C[io] = 75.1863;  fMeanV0M[io] = 125.96;  io++;
   fRuns[io] = 291402;  fMeanV0A[io] = 50.7901;  fMeanV0C[io] = 75.1557;  fMeanV0M[io] = 125.864;  io++;
   fRuns[io] = 291416;  fMeanV0A[io] = 51.0085;  fMeanV0C[io] = 75.337;   fMeanV0M[io] = 126.261;  io++;
   fRuns[io] = 291417;  fMeanV0A[io] = 50.9048;  fMeanV0C[io] = 75.2709;  fMeanV0M[io] = 126.09;  io++;
   fRuns[io] = 291420;  fMeanV0A[io] = 50.7557;  fMeanV0C[io] = 75.0813;  fMeanV0M[io] = 125.757;  io++;
   fRuns[io] = 291424;  fMeanV0A[io] = 50.6108;  fMeanV0C[io] = 74.9915;  fMeanV0M[io] = 125.532;  io++;
   fRuns[io] = 291447;  fMeanV0A[io] = 50.7826;  fMeanV0C[io] = 75.1305;  fMeanV0M[io] = 125.824;  io++;
   fRuns[io] = 291451;  fMeanV0A[io] = 50.7774;  fMeanV0C[io] = 75.0818;  fMeanV0M[io] = 125.776;  io++;
   fRuns[io] = 291453;  fMeanV0A[io] = 50.6252;  fMeanV0C[io] = 74.949;   fMeanV0M[io] = 125.475;  io++;
   fRuns[io] = 291456;  fMeanV0A[io] = 50.5742;  fMeanV0C[io] = 74.874;   fMeanV0M[io] = 125.372;  io++;
   fRuns[io] = 291457;  fMeanV0A[io] = 50.7311;  fMeanV0C[io] = 75.2503;  fMeanV0M[io] = 125.901;  io++;
   fRuns[io] = 291481;  fMeanV0A[io] = 50.6914;  fMeanV0C[io] = 75.0141;  fMeanV0M[io] = 125.621;  io++;
   fRuns[io] = 291482;  fMeanV0A[io] = 50.608;   fMeanV0C[io] = 75.0408;  fMeanV0M[io] = 125.564;  io++;
   fRuns[io] = 291484;  fMeanV0A[io] = 50.5775;  fMeanV0C[io] = 74.916;   fMeanV0M[io] = 125.41;  io++;
   fRuns[io] = 291485;  fMeanV0A[io] = 50.521;   fMeanV0C[io] = 74.8318;  fMeanV0M[io] = 125.287;  io++;
   fRuns[io] = 291590;  fMeanV0A[io] = 50.712;   fMeanV0C[io] = 75.0769;  fMeanV0M[io] = 125.714;  io++;
   fRuns[io] = 291614;  fMeanV0A[io] = 50.7597;  fMeanV0C[io] = 75.2908;  fMeanV0M[io] = 125.968;  io++;
   fRuns[io] = 291615;  fMeanV0A[io] = 50.6086;  fMeanV0C[io] = 75.1972;  fMeanV0M[io] = 125.722;  io++;
   fRuns[io] = 291618;  fMeanV0A[io] = 50.7193;  fMeanV0C[io] = 75.3039;  fMeanV0M[io] = 125.93;  io++;
   fRuns[io] = 291622;  fMeanV0A[io] = 50.6394;  fMeanV0C[io] = 75.2342;  fMeanV0M[io] = 125.783;  io++;
   fRuns[io] = 291624;  fMeanV0A[io] = 50.7174;  fMeanV0C[io] = 75.1785;  fMeanV0M[io] = 125.815;  io++;
   fRuns[io] = 291626;  fMeanV0A[io] = 50.7501;  fMeanV0C[io] = 75.3436;  fMeanV0M[io] = 126.017;  io++;
   fRuns[io] = 291657;  fMeanV0A[io] = 50.6458;  fMeanV0C[io] = 75.0227;  fMeanV0M[io] = 125.591;  io++;
   fRuns[io] = 291661;  fMeanV0A[io] = 50.4051;  fMeanV0C[io] = 74.7819;  fMeanV0M[io] = 125.108;  io++;
   fRuns[io] = 291665;  fMeanV0A[io] = 50.3961;  fMeanV0C[io] = 74.6755;  fMeanV0M[io] = 124.99;  io++;
   fRuns[io] = 291690;  fMeanV0A[io] = 50.5075;  fMeanV0C[io] = 74.8373;  fMeanV0M[io] = 125.256;  io++;
   fRuns[io] = 291697;  fMeanV0A[io] = 50.3674;  fMeanV0C[io] = 74.5537;  fMeanV0M[io] = 124.83;  io++;
   fRuns[io] = 291698;  fMeanV0A[io] = 50.4206;  fMeanV0C[io] = 74.601;   fMeanV0M[io] = 124.929;  io++;
   fRuns[io] = 291706;  fMeanV0A[io] = 50.3676;  fMeanV0C[io] = 74.4799;  fMeanV0M[io] = 124.771;  io++;
   fRuns[io] = 291729;  fMeanV0A[io] = 50.4841;  fMeanV0C[io] = 74.5696;  fMeanV0M[io] = 124.978;  io++;
   fRuns[io] = 291755;  fMeanV0A[io] = 50.421;   fMeanV0C[io] = 74.4562;  fMeanV0M[io] = 124.796;  io++;
   fRuns[io] = 291756;  fMeanV0A[io] = 50.4138;  fMeanV0C[io] = 74.5146;  fMeanV0M[io] = 124.849;  io++;
   fRuns[io] = 291760;  fMeanV0A[io] = 49.9842;  fMeanV0C[io] = 73.9182;  fMeanV0M[io] = 123.837;  io++;
   fRuns[io] = 291762;  fMeanV0A[io] = 50.3519;  fMeanV0C[io] = 74.2138;  fMeanV0M[io] = 124.476;  io++;
   fRuns[io] = 291766;  fMeanV0A[io] = 50.3912;  fMeanV0C[io] = 74.382;   fMeanV0M[io] = 124.693;  io++;
   fRuns[io] = 291768;  fMeanV0A[io] = 50.2169;  fMeanV0C[io] = 74.4569;  fMeanV0M[io] = 124.595;  io++;
   fRuns[io] = 291769;  fMeanV0A[io] = 50.2042;  fMeanV0C[io] = 74.2467;  fMeanV0M[io] = 124.366;  io++;
   fRuns[io] = 291795;  fMeanV0A[io] = 50.2813;  fMeanV0C[io] = 74.194;   fMeanV0M[io] = 124.37;  io++;
   fRuns[io] = 291796;  fMeanV0A[io] = 50.182;   fMeanV0C[io] = 74.2238;  fMeanV0M[io] = 124.315;  io++;
   fRuns[io] = 291803;  fMeanV0A[io] = 50.126;   fMeanV0C[io] = 74.0807;  fMeanV0M[io] = 124.118;  io++;
   fRuns[io] = 291942;  fMeanV0A[io] = 50.8277;  fMeanV0C[io] = 74.9658;  fMeanV0M[io] = 125.705;  io++;
   fRuns[io] = 291943;  fMeanV0A[io] = 50.4185;  fMeanV0C[io] = 74.6935;  fMeanV0M[io] = 125.028;  io++;
   fRuns[io] = 291944;  fMeanV0A[io] = 50.4155;  fMeanV0C[io] = 74.4784;  fMeanV0M[io] = 124.811;  io++;
   fRuns[io] = 291945;  fMeanV0A[io] = 50.2576;  fMeanV0C[io] = 74.5963;  fMeanV0M[io] = 124.77;  io++;
   fRuns[io] = 291946;  fMeanV0A[io] = 50.2668;  fMeanV0C[io] = 74.3811;  fMeanV0M[io] = 124.57;  io++;
   fRuns[io] = 291948;  fMeanV0A[io] = 50.3698;  fMeanV0C[io] = 74.4037;  fMeanV0M[io] = 124.691;  io++;
   fRuns[io] = 291953;  fMeanV0A[io] = 50.2098;  fMeanV0C[io] = 74.4893;  fMeanV0M[io] = 124.619;  io++;
   fRuns[io] = 291976;  fMeanV0A[io] = 50.438;   fMeanV0C[io] = 74.5328;  fMeanV0M[io] = 124.895;  io++;
   fRuns[io] = 291977;  fMeanV0A[io] = 50.3539;  fMeanV0C[io] = 74.5111;  fMeanV0M[io] = 124.773;  io++;
   fRuns[io] = 291982;  fMeanV0A[io] = 50.3484;  fMeanV0C[io] = 74.5091;  fMeanV0M[io] = 124.772;  io++;
   fRuns[io] = 292012;  fMeanV0A[io] = 50.3787;  fMeanV0C[io] = 74.5009;  fMeanV0M[io] = 124.805;  io++;
   fRuns[io] = 292040;  fMeanV0A[io] = 50.2878;  fMeanV0C[io] = 74.4715;  fMeanV0M[io] = 124.68;  io++;
   fRuns[io] = 292060;  fMeanV0A[io] = 50.3762;  fMeanV0C[io] = 74.5415;  fMeanV0M[io] = 124.828;  io++;
   fRuns[io] = 292061;  fMeanV0A[io] = 50.1547;  fMeanV0C[io] = 74.2039;  fMeanV0M[io] = 124.266;  io++;
   fRuns[io] = 292062;  fMeanV0A[io] = 50.0739;  fMeanV0C[io] = 74.1687;  fMeanV0M[io] = 124.151;  io++;
   fRuns[io] = 292067;  fMeanV0A[io] = 50.1455;  fMeanV0C[io] = 74.2078;  fMeanV0M[io] = 124.271;  io++;
   fRuns[io] = 292075;  fMeanV0A[io] = 50.1079;  fMeanV0C[io] = 74.2648;  fMeanV0M[io] = 124.284;  io++;
   fRuns[io] = 292077;  fMeanV0A[io] = 50.0815;  fMeanV0C[io] = 74.2284;  fMeanV0M[io] = 124.227;  io++;
   fRuns[io] = 292080;  fMeanV0A[io] = 50.1456;  fMeanV0C[io] = 74.2549;  fMeanV0M[io] = 124.314;  io++;
   fRuns[io] = 292081;  fMeanV0A[io] = 50.1529;  fMeanV0C[io] = 74.3164;  fMeanV0M[io] = 124.376;  io++;
   fRuns[io] = 292106;  fMeanV0A[io] = 50.2397;  fMeanV0C[io] = 74.3258;  fMeanV0M[io] = 124.468;  io++;
   fRuns[io] = 292107;  fMeanV0A[io] = 50.1127;  fMeanV0C[io] = 74.2079;  fMeanV0M[io] = 124.223;  io++;
   fRuns[io] = 292108;  fMeanV0A[io] = 49.8848;  fMeanV0C[io] = 74.1902;  fMeanV0M[io] = 123.983;  io++;
   fRuns[io] = 292109;  fMeanV0A[io] = 50.0066;  fMeanV0C[io] = 73.9749;  fMeanV0M[io] = 123.895;  io++;
   fRuns[io] = 292114;  fMeanV0A[io] = 50.0039;  fMeanV0C[io] = 74.1549;  fMeanV0M[io] = 124.071;  io++;
   fRuns[io] = 292115;  fMeanV0A[io] = 50.1165;  fMeanV0C[io] = 74.0713;  fMeanV0M[io] = 124.099;  io++;
   fRuns[io] = 292140;  fMeanV0A[io] = 50.612;   fMeanV0C[io] = 74.2423;  fMeanV0M[io] = 124.734;  io++;
   fRuns[io] = 292160;  fMeanV0A[io] = 49.9772;  fMeanV0C[io] = 74.0375;  fMeanV0M[io] = 123.922;  io++;
   fRuns[io] = 292161;  fMeanV0A[io] = 49.8721;  fMeanV0C[io] = 73.9128;  fMeanV0M[io] = 123.697;  io++;
   fRuns[io] = 292162;  fMeanV0A[io] = 50.0978;  fMeanV0C[io] = 74.2378;  fMeanV0M[io] = 124.239;  io++;
   fRuns[io] = 292163;  fMeanV0A[io] = 49.9412;  fMeanV0C[io] = 74.094;   fMeanV0M[io] = 123.943;  io++;
   fRuns[io] = 292164;  fMeanV0A[io] = 49.9112;  fMeanV0C[io] = 73.7908;  fMeanV0M[io] = 123.608;  io++;
   fRuns[io] = 292166;  fMeanV0A[io] = 49.9559;  fMeanV0C[io] = 73.9632;  fMeanV0M[io] = 123.83;  io++;
   fRuns[io] = 292167;  fMeanV0A[io] = 49.8423;  fMeanV0C[io] = 74.1131;  fMeanV0M[io] = 123.876;  io++;
   fRuns[io] = 292168;  fMeanV0A[io] = 49.9679;  fMeanV0C[io] = 73.8901;  fMeanV0M[io] = 123.77;  io++;
   fRuns[io] = 292192;  fMeanV0A[io] = 50.5756;  fMeanV0C[io] = 74.856;   fMeanV0M[io] = 125.331;  io++;
   fRuns[io] = 292218;  fMeanV0A[io] = 50.3083;  fMeanV0C[io] = 75.0763;  fMeanV0M[io] = 125.332;  io++;
   fRuns[io] = 292240;  fMeanV0A[io] = 50.0515;  fMeanV0C[io] = 74.0994;  fMeanV0M[io] = 124.064;  io++;
   fRuns[io] = 292241;  fMeanV0A[io] = 49.8018;  fMeanV0C[io] = 73.6494;  fMeanV0M[io] = 123.351;  io++;
   fRuns[io] = 292242;  fMeanV0A[io] = 50.0911;  fMeanV0C[io] = 74.1197;  fMeanV0M[io] = 124.127;  io++;
   fRuns[io] = 292265;  fMeanV0A[io] = 50.3178;  fMeanV0C[io] = 74.4806;  fMeanV0M[io] = 124.696;  io++;
   fRuns[io] = 292273;  fMeanV0A[io] = 50.0132;  fMeanV0C[io] = 73.9449;  fMeanV0M[io] = 123.872;  io++;
   fRuns[io] = 292298;  fMeanV0A[io] = 49.9694;  fMeanV0C[io] = 74.0659;  fMeanV0M[io] = 123.948;  io++;
   fRuns[io] = 292397;  fMeanV0A[io] = 50.1524;  fMeanV0C[io] = 73.9305;  fMeanV0M[io] = 123.994;  io++;
   fRuns[io] = 292398;  fMeanV0A[io] = 49.8891;  fMeanV0C[io] = 73.6731;  fMeanV0M[io] = 123.476;  io++;
   fRuns[io] = 292405;  fMeanV0A[io] = 49.7805;  fMeanV0C[io] = 73.8791;  fMeanV0M[io] = 123.57;  io++;
   fRuns[io] = 292406;  fMeanV0A[io] = 49.8001;  fMeanV0C[io] = 73.7708;  fMeanV0M[io] = 123.471;  io++;
   fRuns[io] = 292428;  fMeanV0A[io] = 49.9245;  fMeanV0C[io] = 73.7608;  fMeanV0M[io] = 123.599;  io++;
   fRuns[io] = 292429;  fMeanV0A[io] = 49.9626;  fMeanV0C[io] = 73.8308;  fMeanV0M[io] = 123.707;  io++;
   fRuns[io] = 292430;  fMeanV0A[io] = 49.8629;  fMeanV0C[io] = 73.736;   fMeanV0M[io] = 123.507;  io++;
   fRuns[io] = 292432;  fMeanV0A[io] = 49.6795;  fMeanV0C[io] = 73.5804;  fMeanV0M[io] = 123.166;  io++;
   fRuns[io] = 292434;  fMeanV0A[io] = 49.866;   fMeanV0C[io] = 73.8708;  fMeanV0M[io] = 123.641;  io++;
   fRuns[io] = 292456;  fMeanV0A[io] = 49.9206;  fMeanV0C[io] = 73.6988;  fMeanV0M[io] = 123.524;  io++;
   fRuns[io] = 292457;  fMeanV0A[io] = 49.7806;  fMeanV0C[io] = 73.7145;  fMeanV0M[io] = 123.407;  io++;
   fRuns[io] = 292460;  fMeanV0A[io] = 49.7266;  fMeanV0C[io] = 73.7688;  fMeanV0M[io] = 123.413;  io++;
   fRuns[io] = 292461;  fMeanV0A[io] = 49.6857;  fMeanV0C[io] = 73.466;   fMeanV0M[io] = 123.067;  io++;
   fRuns[io] = 292495;  fMeanV0A[io] = 49.9439;  fMeanV0C[io] = 73.7631;  fMeanV0M[io] = 123.614;  io++;
   fRuns[io] = 292496;  fMeanV0A[io] = 49.7652;  fMeanV0C[io] = 73.6415;  fMeanV0M[io] = 123.317;  io++;
   fRuns[io] = 292497;  fMeanV0A[io] = 49.7462;  fMeanV0C[io] = 73.5187;  fMeanV0M[io] = 123.171;  io++;
   fRuns[io] = 292500;  fMeanV0A[io] = 49.6902;  fMeanV0C[io] = 73.6213;  fMeanV0M[io] = 123.221;  io++;
   fRuns[io] = 292521;  fMeanV0A[io] = 49.832;   fMeanV0C[io] = 73.8112;  fMeanV0M[io] = 123.566;  io++;
   fRuns[io] = 292523;  fMeanV0A[io] = 49.7378;  fMeanV0C[io] = 73.6219;  fMeanV0M[io] = 123.267;  io++;
   fRuns[io] = 292524;  fMeanV0A[io] = 49.5803;  fMeanV0C[io] = 73.6883;  fMeanV0M[io] = 123.183;  io++;
   fRuns[io] = 292526;  fMeanV0A[io] = 49.7314;  fMeanV0C[io] = 73.7769;  fMeanV0M[io] = 123.428;  io++;
   fRuns[io] = 292553;  fMeanV0A[io] = 49.9847;  fMeanV0C[io] = 73.7035;  fMeanV0M[io] = 123.598;  io++;
   fRuns[io] = 292554;  fMeanV0A[io] = 49.7367;  fMeanV0C[io] = 73.6218;  fMeanV0M[io] = 123.278;  io++;
   fRuns[io] = 292557;  fMeanV0A[io] = 49.7533;  fMeanV0C[io] = 73.6644;  fMeanV0M[io] = 123.331;  io++;
   fRuns[io] = 292559;  fMeanV0A[io] = 49.6587;  fMeanV0C[io] = 73.5529;  fMeanV0M[io] = 123.126;  io++;
   fRuns[io] = 292560;  fMeanV0A[io] = 49.7035;  fMeanV0C[io] = 73.5564;  fMeanV0M[io] = 123.182;  io++;
   fRuns[io] = 292563;  fMeanV0A[io] = 49.4707;  fMeanV0C[io] = 73.5671;  fMeanV0M[io] = 122.952;  io++;
   fRuns[io] = 292584;  fMeanV0A[io] = 49.9365;  fMeanV0C[io] = 73.8223;  fMeanV0M[io] = 123.674;  io++;
   fRuns[io] = 292586;  fMeanV0A[io] = 49.631;   fMeanV0C[io] = 73.7053;  fMeanV0M[io] = 123.257;  io++;
   fRuns[io] = 292693;  fMeanV0A[io] = 49.5756;  fMeanV0C[io] = 73.4867;  fMeanV0M[io] = 122.98;  io++;
   fRuns[io] = 292695;  fMeanV0A[io] = 49.4186;  fMeanV0C[io] = 73.2884;  fMeanV0M[io] = 122.594;  io++;
   fRuns[io] = 292696;  fMeanV0A[io] = 49.4866;  fMeanV0C[io] = 73.3803;  fMeanV0M[io] = 122.779;  io++;
   fRuns[io] = 292698;  fMeanV0A[io] = 49.739;   fMeanV0C[io] = 73.5405;  fMeanV0M[io] = 123.184;  io++;
   fRuns[io] = 292701;  fMeanV0A[io] = 49.6586;  fMeanV0C[io] = 73.4206;  fMeanV0M[io] = 122.987;  io++;
   fRuns[io] = 292704;  fMeanV0A[io] = 49.5705;  fMeanV0C[io] = 73.2863;  fMeanV0M[io] = 122.763;  io++;
   fRuns[io] = 292737;  fMeanV0A[io] = 49.7659;  fMeanV0C[io] = 73.4917;  fMeanV0M[io] = 123.167;  io++;
   fRuns[io] = 292739;  fMeanV0A[io] = 49.6033;  fMeanV0C[io] = 73.4368;  fMeanV0M[io] = 122.955;  io++;
   fRuns[io] = 292744;  fMeanV0A[io] = 49.59;    fMeanV0C[io] = 73.4275;  fMeanV0M[io] = 122.92;  io++;
   fRuns[io] = 292747;  fMeanV0A[io] = 49.4789;  fMeanV0C[io] = 73.3995;  fMeanV0M[io] = 122.78;  io++;
   fRuns[io] = 292748;  fMeanV0A[io] = 49.5781;  fMeanV0C[io] = 73.5499;  fMeanV0M[io] = 123.036;  io++;
   fRuns[io] = 292750;  fMeanV0A[io] = 49.4233;  fMeanV0C[io] = 73.5836;  fMeanV0M[io] = 122.882;  io++;
   fRuns[io] = 292752;  fMeanV0A[io] = 49.5382;  fMeanV0C[io] = 73.6201;  fMeanV0M[io] = 123.068;  io++;
   fRuns[io] = 292803;  fMeanV0A[io] = 49.8223;  fMeanV0C[io] = 73.8253;  fMeanV0M[io] = 123.566;  io++;
   fRuns[io] = 292804;  fMeanV0A[io] = 49.6846;  fMeanV0C[io] = 73.6316;  fMeanV0M[io] = 123.224;  io++;
   fRuns[io] = 292809;  fMeanV0A[io] = 49.6753;  fMeanV0C[io] = 73.3522;  fMeanV0M[io] = 122.962;  io++;
   fRuns[io] = 292810;  fMeanV0A[io] = 49.5036;  fMeanV0C[io] = 73.4954;  fMeanV0M[io] = 122.931;  io++;
   fRuns[io] = 292811;  fMeanV0A[io] = 49.4969;  fMeanV0C[io] = 73.3655;  fMeanV0M[io] = 122.783;  io++;
   fRuns[io] = 292831;  fMeanV0A[io] = 49.6615;  fMeanV0C[io] = 73.5562;  fMeanV0M[io] = 123.135;  io++;
   fRuns[io] = 292832;  fMeanV0A[io] = 49.6061;  fMeanV0C[io] = 73.3479;  fMeanV0M[io] = 122.875;  io++;
   fRuns[io] = 292834;  fMeanV0A[io] = 49.5875;  fMeanV0C[io] = 73.6374;  fMeanV0M[io] = 123.141;  io++;
   fRuns[io] = 292836;  fMeanV0A[io] = 49.5534;  fMeanV0C[io] = 73.2622;  fMeanV0M[io] = 122.736;  io++;
   fRuns[io] = 292839;  fMeanV0A[io] = 49.5916;  fMeanV0C[io] = 73.5847;  fMeanV0M[io] = 123.087;  io++;


   // FILTER_p-p_208_LHC18n
   fRuns[io] = 293357;  fMeanV0A[io] = 49.8311;  fMeanV0C[io] = 74.5382;  fMeanV0M[io] = 124.285;  io++;
   fRuns[io] = 293359;  fMeanV0A[io] = 49.637;   fMeanV0C[io] = 74.0711;  fMeanV0M[io] = 123.62;  io++;

   // FILTER_p-p_208_LHC18o
   fRuns[io] = 293475;  fMeanV0A[io] = 49.7659;  fMeanV0C[io] = 73.744;   fMeanV0M[io] = 123.424;  io++;
   fRuns[io] = 293570;  fMeanV0A[io] = 49.5091;  fMeanV0C[io] = 73.2193;  fMeanV0M[io] = 122.618;  io++;
   fRuns[io] = 293571;  fMeanV0A[io] = 49.3074;  fMeanV0C[io] = 73.0624;  fMeanV0M[io] = 122.259;  io++;
   fRuns[io] = 293573;  fMeanV0A[io] = 49.2341;  fMeanV0C[io] = 73.2035;  fMeanV0M[io] = 122.36;  io++;
   fRuns[io] = 293578;  fMeanV0A[io] = 49.4036;  fMeanV0C[io] = 73.2991;  fMeanV0M[io] = 122.624;  io++;
   fRuns[io] = 293579;  fMeanV0A[io] = 49.2845;  fMeanV0C[io] = 73.1327;  fMeanV0M[io] = 122.33;  io++;
   fRuns[io] = 293582;  fMeanV0A[io] = 49.3292;  fMeanV0C[io] = 73.3216;  fMeanV0M[io] = 122.568;  io++;
   fRuns[io] = 293583;  fMeanV0A[io] = 49.5044;  fMeanV0C[io] = 73.3352;  fMeanV0M[io] = 122.735;  io++;
   fRuns[io] = 293587;  fMeanV0A[io] = 49.4712;  fMeanV0C[io] = 73.3275;  fMeanV0M[io] = 122.708;  io++;
   fRuns[io] = 293588;  fMeanV0A[io] = 49.5661;  fMeanV0C[io] = 73.3643;  fMeanV0M[io] = 122.841;  io++;
   fRuns[io] = 293691;  fMeanV0A[io] = 49.7579;  fMeanV0C[io] = 73.6636;  fMeanV0M[io] = 123.333;  io++;
   fRuns[io] = 293692;  fMeanV0A[io] = 49.6233;  fMeanV0C[io] = 73.5947;  fMeanV0M[io] = 123.133;  io++;
   fRuns[io] = 293695;  fMeanV0A[io] = 49.6024;  fMeanV0C[io] = 73.6045;  fMeanV0M[io] = 123.122;  io++;
   fRuns[io] = 293696;  fMeanV0A[io] = 49.7104;  fMeanV0C[io] = 73.4656;  fMeanV0M[io] = 123.089;  io++;
   fRuns[io] = 293698;  fMeanV0A[io] = 49.5285;  fMeanV0C[io] = 73.6158;  fMeanV0M[io] = 123.069;  io++;
   fRuns[io] = 293740;  fMeanV0A[io] = 49.7938;  fMeanV0C[io] = 73.6502;  fMeanV0M[io] = 123.357;  io++;
   fRuns[io] = 293741;  fMeanV0A[io] = 49.6388;  fMeanV0C[io] = 73.5383;  fMeanV0M[io] = 123.089;  io++;
   fRuns[io] = 293770;  fMeanV0A[io] = 49.9619;  fMeanV0C[io] = 74.0614;  fMeanV0M[io] = 123.931;  io++;
   fRuns[io] = 293773;  fMeanV0A[io] = 49.7797;  fMeanV0C[io] = 73.6188;  fMeanV0M[io] = 123.31;  io++;
   fRuns[io] = 293774;  fMeanV0A[io] = 49.614;   fMeanV0C[io] = 73.4814;  fMeanV0M[io] = 123.015;  io++;
   fRuns[io] = 293776;  fMeanV0A[io] = 49.5406;  fMeanV0C[io] = 73.6128;  fMeanV0M[io] = 123.056;  io++;
   fRuns[io] = 293802;  fMeanV0A[io] = 49.648;   fMeanV0C[io] = 73.5199;  fMeanV0M[io] = 123.081;  io++;
   fRuns[io] = 293805;  fMeanV0A[io] = 49.5808;  fMeanV0C[io] = 73.5204;  fMeanV0M[io] = 123.016;  io++;
   fRuns[io] = 293806;  fMeanV0A[io] = 49.5885;  fMeanV0C[io] = 74.114;   fMeanV0M[io] = 123.61;  io++;
   fRuns[io] = 293807;  fMeanV0A[io] = 49.4301;  fMeanV0C[io] = 73.4441;  fMeanV0M[io] = 122.787;  io++;
   fRuns[io] = 293809;  fMeanV0A[io] = 49.4911;  fMeanV0C[io] = 73.4561;  fMeanV0M[io] = 122.849;  io++;
   fRuns[io] = 293829;  fMeanV0A[io] = 49.7735;  fMeanV0C[io] = 73.7392;  fMeanV0M[io] = 123.418;  io++;
   fRuns[io] = 293830;  fMeanV0A[io] = 49.6545;  fMeanV0C[io] = 73.6278;  fMeanV0M[io] = 123.194;  io++;
   fRuns[io] = 293831;  fMeanV0A[io] = 49.4882;  fMeanV0C[io] = 73.4254;  fMeanV0M[io] = 122.819;  io++;
   fRuns[io] = 293856;  fMeanV0A[io] = 49.5936;  fMeanV0C[io] = 73.4857;  fMeanV0M[io] = 122.991;  io++;
   fRuns[io] = 293886;  fMeanV0A[io] = 49.6957;  fMeanV0C[io] = 73.5106;  fMeanV0M[io] = 123.113;  io++;
   fRuns[io] = 293891;  fMeanV0A[io] = 49.3805;  fMeanV0C[io] = 73.0865;  fMeanV0M[io] = 122.384;  io++;
   fRuns[io] = 293893;  fMeanV0A[io] = 49.3858;  fMeanV0C[io] = 73.1276;  fMeanV0M[io] = 122.42;  io++;
   fRuns[io] = 293896;  fMeanV0A[io] = 49.4432;  fMeanV0C[io] = 73.2277;  fMeanV0M[io] = 122.584;  io++;
   fRuns[io] = 293898;  fMeanV0A[io] = 49.4526;  fMeanV0C[io] = 73.3552;  fMeanV0M[io] = 122.722;  io++;

   // FILTER_p-p_208_LHC18p
   fRuns[io] = 294009;  fMeanV0A[io] = 49.4815;  fMeanV0C[io] = 73.2223;  fMeanV0M[io] = 122.6;  io++;
   fRuns[io] = 294010;  fMeanV0A[io] = 49.2818;  fMeanV0C[io] = 72.7721;  fMeanV0M[io] = 121.97;  io++;
   fRuns[io] = 294011;  fMeanV0A[io] = 49.2622;  fMeanV0C[io] = 72.7619;  fMeanV0M[io] = 121.943;  io++;
   fRuns[io] = 294012;  fMeanV0A[io] = 49.0446;  fMeanV0C[io] = 72.696;   fMeanV0M[io] = 121.651;  io++;
   fRuns[io] = 294013;  fMeanV0A[io] = 49.1902;  fMeanV0C[io] = 72.8259;  fMeanV0M[io] = 121.924;  io++;
   fRuns[io] = 294131;  fMeanV0A[io] = 49.1402;  fMeanV0C[io] = 72.9002;  fMeanV0M[io] = 121.946;  io++;
   fRuns[io] = 294152;  fMeanV0A[io] = 49.115;   fMeanV0C[io] = 72.7627;  fMeanV0M[io] = 121.788;  io++;
   fRuns[io] = 294154;  fMeanV0A[io] = 49.0092;  fMeanV0C[io] = 72.4962;  fMeanV0M[io] = 121.419;  io++;
   fRuns[io] = 294155;  fMeanV0A[io] = 48.9969;  fMeanV0C[io] = 72.7162;  fMeanV0M[io] = 121.641;  io++;
   fRuns[io] = 294156;  fMeanV0A[io] = 49.0983;  fMeanV0C[io] = 72.6905;  fMeanV0M[io] = 121.694;  io++;
   fRuns[io] = 294199;  fMeanV0A[io] = 49.268;   fMeanV0C[io] = 72.8685;  fMeanV0M[io] = 122.043;  io++;
   fRuns[io] = 294200;  fMeanV0A[io] = 49.1882;  fMeanV0C[io] = 72.8406;  fMeanV0M[io] = 121.934;  io++;
   fRuns[io] = 294201;  fMeanV0A[io] = 49.1848;  fMeanV0C[io] = 72.7845;  fMeanV0M[io] = 121.879;  io++;
   fRuns[io] = 294205;  fMeanV0A[io] = 49.2267;  fMeanV0C[io] = 72.909;   fMeanV0M[io] = 122.043;  io++;
   fRuns[io] = 294208;  fMeanV0A[io] = 49.2588;  fMeanV0C[io] = 72.9157;  fMeanV0M[io] = 122.093;  io++;
   fRuns[io] = 294210;  fMeanV0A[io] = 49.4136;  fMeanV0C[io] = 72.8864;  fMeanV0M[io] = 122.236;  io++;
   fRuns[io] = 294212;  fMeanV0A[io] = 49.1322;  fMeanV0C[io] = 72.7464;  fMeanV0M[io] = 121.792;  io++;
   fRuns[io] = 294241;  fMeanV0A[io] = 49.2994;  fMeanV0C[io] = 72.9824;  fMeanV0M[io] = 122.189;  io++;
   fRuns[io] = 294242;  fMeanV0A[io] = 49.1955;  fMeanV0C[io] = 73.0759;  fMeanV0M[io] = 122.179;  io++;
   fRuns[io] = 294307;  fMeanV0A[io] = 49.0523;  fMeanV0C[io] = 72.8586;  fMeanV0M[io] = 121.825;  io++;
   fRuns[io] = 294308;  fMeanV0A[io] = 49.0749;  fMeanV0C[io] = 72.8544;  fMeanV0M[io] = 121.833;  io++;
   fRuns[io] = 294310;  fMeanV0A[io] = 49.1419;  fMeanV0C[io] = 72.8778;  fMeanV0M[io] = 121.932;  io++;
   fRuns[io] = 294524;  fMeanV0A[io] = 48.8067;  fMeanV0C[io] = 72.9134;  fMeanV0M[io] = 121.652;  io++;
   fRuns[io] = 294525;  fMeanV0A[io] = 49.0437;  fMeanV0C[io] = 72.8546;  fMeanV0M[io] = 121.805;  io++;
   fRuns[io] = 294526;  fMeanV0A[io] = 49.006;   fMeanV0C[io] = 72.7466;  fMeanV0M[io] = 121.657;  io++;
   fRuns[io] = 294527;  fMeanV0A[io] = 49.023;   fMeanV0C[io] = 72.9067;  fMeanV0M[io] = 121.844;  io++;
   fRuns[io] = 294529;  fMeanV0A[io] = 49.1457;  fMeanV0C[io] = 73.0901;  fMeanV0M[io] = 122.153;  io++;
   fRuns[io] = 294530;  fMeanV0A[io] = 49.2182;  fMeanV0C[io] = 73.1218;  fMeanV0M[io] = 122.249;  io++;
   fRuns[io] = 294531;  fMeanV0A[io] = 49.1358;  fMeanV0C[io] = 73.0355;  fMeanV0M[io] = 122.081;  io++;
   fRuns[io] = 294553;  fMeanV0A[io] = 49.3222;  fMeanV0C[io] = 73.2748;  fMeanV0M[io] = 122.508;  io++;
   fRuns[io] = 294556;  fMeanV0A[io] = 49.3066;  fMeanV0C[io] = 72.9958;  fMeanV0M[io] = 122.209;  io++;
   fRuns[io] = 294558;  fMeanV0A[io] = 49.2026;  fMeanV0C[io] = 72.9602;  fMeanV0M[io] = 122.08;  io++;
   fRuns[io] = 294562;  fMeanV0A[io] = 49.4496;  fMeanV0C[io] = 73.4219;  fMeanV0M[io] = 122.791;  io++;
   fRuns[io] = 294563;  fMeanV0A[io] = 49.1275;  fMeanV0C[io] = 72.9041;  fMeanV0M[io] = 121.94;  io++;
   fRuns[io] = 294586;  fMeanV0A[io] = 49.1613;  fMeanV0C[io] = 72.8544;  fMeanV0M[io] = 121.921;  io++;
   fRuns[io] = 294587;  fMeanV0A[io] = 49.0328;  fMeanV0C[io] = 72.7764;  fMeanV0M[io] = 121.724;  io++;
   fRuns[io] = 294590;  fMeanV0A[io] = 49.131;   fMeanV0C[io] = 72.7516;  fMeanV0M[io] = 121.796;  io++;
   fRuns[io] = 294591;  fMeanV0A[io] = 49.0686;  fMeanV0C[io] = 72.7167;  fMeanV0M[io] = 121.701;  io++;
   fRuns[io] = 294593;  fMeanV0A[io] = 49.1431;  fMeanV0C[io] = 72.8436;  fMeanV0M[io] = 121.896;  io++;
   fRuns[io] = 294632;  fMeanV0A[io] = 49.0821;  fMeanV0C[io] = 72.88;    fMeanV0M[io] = 121.869;  io++;
   fRuns[io] = 294633;  fMeanV0A[io] = 49.0115;  fMeanV0C[io] = 72.7179;  fMeanV0M[io] = 121.639;  io++;
   fRuns[io] = 294636;  fMeanV0A[io] = 48.9452;  fMeanV0C[io] = 72.6133;  fMeanV0M[io] = 121.469;  io++;
   fRuns[io] = 294653;  fMeanV0A[io] = 49.3018;  fMeanV0C[io] = 73.0849;  fMeanV0M[io] = 122.296;  io++;
   fRuns[io] = 294703;  fMeanV0A[io] = 49.479;   fMeanV0C[io] = 73.3472;  fMeanV0M[io] = 122.721;  io++;
   fRuns[io] = 294710;  fMeanV0A[io] = 49.058;   fMeanV0C[io] = 72.6571;  fMeanV0M[io] = 121.62;  io++;
   fRuns[io] = 294715;  fMeanV0A[io] = 49.151;   fMeanV0C[io] = 72.8188;  fMeanV0M[io] = 121.868;  io++;
   fRuns[io] = 294718;  fMeanV0A[io] = 48.9351;  fMeanV0C[io] = 72.716;   fMeanV0M[io] = 121.549;  io++;
   fRuns[io] = 294722;  fMeanV0A[io] = 48.7994;  fMeanV0C[io] = 72.5268;  fMeanV0M[io] = 121.239;  io++;
   fRuns[io] = 294741;  fMeanV0A[io] = 49.0705;  fMeanV0C[io] = 72.4611;  fMeanV0M[io] = 121.429;  io++;
   fRuns[io] = 294742;  fMeanV0A[io] = 48.9882;  fMeanV0C[io] = 72.47;    fMeanV0M[io] = 121.356;  io++;
   fRuns[io] = 294744;  fMeanV0A[io] = 48.9959;  fMeanV0C[io] = 72.3949;  fMeanV0M[io] = 121.293;  io++;
   fRuns[io] = 294745;  fMeanV0A[io] = 48.8054;  fMeanV0C[io] = 72.3145;  fMeanV0M[io] = 121.02;  io++;
   fRuns[io] = 294746;  fMeanV0A[io] = 48.6852;  fMeanV0C[io] = 72.2943;  fMeanV0M[io] = 120.897;  io++;
   fRuns[io] = 294747;  fMeanV0A[io] = 48.7699;  fMeanV0C[io] = 72.3131;  fMeanV0M[io] = 120.988;  io++;
   fRuns[io] = 294749;  fMeanV0A[io] = 48.6504;  fMeanV0C[io] = 71.9548;  fMeanV0M[io] = 120.505;  io++;
   fRuns[io] = 294769;  fMeanV0A[io] = 49.021;   fMeanV0C[io] = 72.7553;  fMeanV0M[io] = 121.69;  io++;
   fRuns[io] = 294772;  fMeanV0A[io] = 48.9768;  fMeanV0C[io] = 72.5788;  fMeanV0M[io] = 121.475;  io++;
   fRuns[io] = 294774;  fMeanV0A[io] = 48.8557;  fMeanV0C[io] = 72.349;   fMeanV0M[io] = 121.114;  io++;
   fRuns[io] = 294775;  fMeanV0A[io] = 48.9028;  fMeanV0C[io] = 72.3289;  fMeanV0M[io] = 121.133;  io++;
   fRuns[io] = 294805;  fMeanV0A[io] = 49.2338;  fMeanV0C[io] = 72.564;   fMeanV0M[io] = 121.698;  io++;
   fRuns[io] = 294809;  fMeanV0A[io] = 48.8531;  fMeanV0C[io] = 72.1744;  fMeanV0M[io] = 120.939;  io++;
   fRuns[io] = 294813;  fMeanV0A[io] = 48.8999;  fMeanV0C[io] = 72.3289;  fMeanV0M[io] = 121.132;  io++;
   fRuns[io] = 294815;  fMeanV0A[io] = 25.9875;  fMeanV0C[io] = 11.6106;  fMeanV0M[io] = 37.3443;  io++;
   fRuns[io] = 294816;  fMeanV0A[io] = 48.8221;  fMeanV0C[io] = 72.1895;  fMeanV0M[io] = 120.924;  io++;
   fRuns[io] = 294817;  fMeanV0A[io] = 48.9095;  fMeanV0C[io] = 72.3253;  fMeanV0M[io] = 121.143;  io++;
   fRuns[io] = 294818;  fMeanV0A[io] = 48.8407;  fMeanV0C[io] = 72.3275;  fMeanV0M[io] = 121.099;  io++;
   fRuns[io] = 294852;  fMeanV0A[io] = 49.0492;  fMeanV0C[io] = 72.4037;  fMeanV0M[io] = 121.369;  io++;
   fRuns[io] = 294875;  fMeanV0A[io] = 48.9251;  fMeanV0C[io] = 72.2737;  fMeanV0M[io] = 121.11;  io++;
   fRuns[io] = 294880;  fMeanV0A[io] = 48.8054;  fMeanV0C[io] = 72.0503;  fMeanV0M[io] = 120.767;  io++;
   fRuns[io] = 294883;  fMeanV0A[io] = 48.7705;  fMeanV0C[io] = 72.0967;  fMeanV0M[io] = 120.768;  io++;
   fRuns[io] = 294884;  fMeanV0A[io] = 48.8807;  fMeanV0C[io] = 72.2328;  fMeanV0M[io] = 121.026;  io++;
   fRuns[io] = 294916;  fMeanV0A[io] = 48.7493;  fMeanV0C[io] = 72.0047;  fMeanV0M[io] = 120.662;  io++;
   fRuns[io] = 294925;  fMeanV0A[io] = 48.641;   fMeanV0C[io] = 71.9452;  fMeanV0M[io] = 120.491;  io++;




   fnRun = io;
}
//_____________________________________________________________________________________________

Double_t AliAnalysisEmcalJetHelperEA::GetV0M(Int_t runnumber) const { 

   //get V0M for the given run number
   Long64_t irun = 0;

   if(runnumber < fRuns[0]){
       printf("AliAnalysisEmcalJetHelperEA: RUN NOT FOUND  %d", runnumber); //index will remain 0
   }else{

      irun = TMath::BinarySearch((Long64_t) fnRun, fRuns.GetArray(), runnumber); //index of the given run number

      if(fRuns[irun] != runnumber){

         printf("AliAnalysisEmcalJetHelperEA: RUN NOT FOUND  %d", runnumber);

         if(runnumber > fRuns[fnRun-1])  irun = fnRun-1;  //index will correspond to the last run in the list
      }
   }

   return  fMeanV0M[irun]; 
}

//_____________________________________________________________________________________________

Double_t AliAnalysisEmcalJetHelperEA::GetV0A(Int_t runnumber) const { 

   //get V0A for the given run number
   Long64_t irun = 0;

   if(runnumber < fRuns[0]){
       printf("AliAnalysisEmcalJetHelperEA: RUN NOT FOUND  %d", runnumber); //index will remain 0
   }else{

      irun = TMath::BinarySearch((Long64_t) fnRun, fRuns.GetArray(), runnumber); //index of the given run number

      if(fRuns[irun] != runnumber){

         printf("AliAnalysisEmcalJetHelperEA: RUN NOT FOUND  %d", runnumber);

         if(runnumber > fRuns[fnRun-1])  irun = fnRun-1;  //index will correspond to the last run in the list
      }
   }

   return  fMeanV0A[irun]; 
}

//_____________________________________________________________________________________________

Double_t AliAnalysisEmcalJetHelperEA::GetV0C(Int_t runnumber) const { 

   //get V0C for the given run number
   Long64_t irun = 0;

   if(runnumber < fRuns[0]){
       printf("AliAnalysisEmcalJetHelperEA: RUN NOT FOUND  %d", runnumber); //index will remain 0
   }else{

      irun = TMath::BinarySearch((Long64_t) fnRun, fRuns.GetArray(), runnumber); //index of the given run number

      if(fRuns[irun] != runnumber){

         printf("AliAnalysisEmcalJetHelperEA: RUN NOT FOUND  %d", runnumber);

         if(runnumber > fRuns[fnRun-1])  irun = fnRun-1;  //index will correspond to the last run in the list
      }
   }

   return  fMeanV0C[irun]; 
}







