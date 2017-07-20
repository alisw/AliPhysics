#include <Riostream.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1F.h>
#include <TVector.h>
#include <TRefArray.h>
#include <TArrayS.h>
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraphErrors.h"
#include "TPostScript.h"
#include "TLegend.h"
#include "TH2I.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TFile.h"
#include "TMath.h"
#include "TLeaf.h" 
#include "TBranch.h" 

//#include "/cebaf/cebaf/EMCAL/cosmicsAnalysis/macros/defineMyPalette40All.C"

 char SMP2Name[][100]={"SMA0","SMC0","SMA1","SMC1","SMA2","SMC2","SMA3","SMC3","SMA4","SMC4","SMA5","SMC5","SMA9","SMC9","SMA10","SMC10","SMA11","SMC11","SMA12","SMC12"};
 char SMnumber[][100]={"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"};

 enum detType {kEMCAL,kEMCALthird,kDCAL,kDCALthird};
 int detTypeType[]={kEMCAL,kEMCALthird,kDCAL,kDCALthird};
 char detTypeString[][100]={"EMCAL","EMCALthird","DCAL","DCALthird"};
 int SMdetType[]={kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCALthird,kEMCALthird,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCALthird,kDCALthird};
 const int kNbColEMCAL=48;
 const int kNbRowEMCAL=24;
 const int kNbSMEMCAL=10;
 const int kNbColEMCALthird=kNbColEMCAL;
 const int kNbRowEMCALthird=(int)(kNbRowEMCAL/3);
 const int kNbSMEMCALthird=2;
 const int kNbColDCAL=32;
 const int kNbRowDCAL=kNbRowEMCAL;
 const int kNbSMDCAL=6;
 const int kNbColDCALthird=kNbColEMCALthird;
 const int kNbRowDCALthird=kNbRowEMCALthird;
 const int kNbSMDCALthird=2;
 const int kNbSMtot=kNbSMEMCAL+kNbSMEMCALthird+kNbSMDCAL+kNbSMDCALthird;
 const int kTabNbCol[4]={kNbColEMCAL,kNbColEMCALthird,kNbColDCAL,kNbColDCALthird};
 const int kTabNbRow[4]={kNbRowEMCAL,kNbRowEMCALthird,kNbRowDCAL,kNbRowDCALthird};
 const int kTabNbSM[4]={kNbSMEMCAL,kNbSMEMCALthird,kNbSMDCAL,kNbSMDCALthird};
 const int kNbColMax=kNbColEMCAL;
 const int kNbRowMax=kNbRowEMCAL;
 const int kNbColOffsetDCAL=kNbColEMCAL-kNbColDCAL;
 





 const int kNbFitParams=6;
 const int kNbExtraParams=7;
 //const int kNbExtraParamsToBeRead=4; //Previous version of the code had 10 parameters written in output_calibPi0_parameters.txt, up to the analysis of input/pass3/ (included).
 const int kNbExtraParamsToBeRead=7; //Valid from the analysis of input/pass4/.
 const int kNbTotParams=kNbFitParams+kNbExtraParams+1; //We add the calib coeff.



///
/// \file Pi0CalibInvMassAnalysis3.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief Pi0CalibInvMassAnalysis3
///
/// Compiled execution :
///   int a=0b0101;
///   aliroot -b -q 'macros/pi0CalibInvMassAnalysis3.C++(a)' >& output_calibPi0.txt &
/// or else :
///   aliroot -b -q 'macros/pi0CalibInvMassAnalysis3.C++()' >& output_calibPi0.txt &
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
//__________________________________________________________________
void  Pi0CalibInvMassAnalysis3(int choice=0b0110)
//choice indicates which SM are used :
//  10^0 = EMCAL, 10^1 = EMCAL thirds, 10^2 = DCAL, 10^3 = DCAL thirds,
//  0 = no, 1 = yes.
//Therefore : old EMCAL SMs only -> choice = 0b0001 ;
//            thirds only -> choice = 0b1010 ;
//            DCAL only -> choice = 0b0100.
//Caution, do not write eg '0001', because integers with a leading zero in C are interpreted as octal values.
//Instead, write '0b0001' to state it's binary.
{int m,n,i,j,k,iCol,iRow,iSM,iPt,iStruct,iBin,discardFlag;
 int testChoice,isFirstIteration,shiftCol,shiftRow,cmpt,flag[kNbFitParams+kNbExtraParams],sm2,c2,r2,cmpt2,flag2[kNbFitParams+kNbExtraParams],tmpFlag,tmpFlagTab[23],cmpt3,sm3,c3,r3,cmptEntriesDiscard,cmptAmpDiscard,cmptMasked,cmptEdge,cmptNDF;
 int kNbColMinEff,kNbColMaxEff,cEff;
 int tabChoiceCalos[4],tmpChoice,nbTowersConsideredTot,choiceNoEMCAL,choiceNoDCAL;
 int kCalibrateMidRapBorderCol,kCalibrateLastBorderCol,kCalibrateTouchingBorderRow,kCalibrateLastBorderRow,flagBorderTower;
 double tabMin[kNbFitParams+kNbExtraParams],tabMax[kNbFitParams+kNbExtraParams],nbEntriesMax,cutMin[kNbFitParams+kNbExtraParams],cutMax[kNbFitParams+kNbExtraParams],param2[kNbTotParams],tabMinDiff[9],tabMaxDiff[9],cutEntriesDiscard,cutAmpDiscard;
 double tmpX,tmpY,fitPar3,fitPar4,fitPar5,nbTot,nbSig,intgS,intgSbincounting,intgN,intgKhi2,coeff,coeff3,maxGaussAmpl,maxIntgAmpl;
 double coeffPower,coeffPowerPrev,fitRangeMin,fitRangeMax,maxHisto,fitResultMean,fitResultSigma,fitKhi2SmallRange,fitKhi2PeakLeft,fitKhi2PeakRight,valStatUncert,paramFitStatUncert_a,paramFitStatUncert_b,uncertKSpec,uncertSpec,uncertDist;
 char txtPS[100];
 
 int kTabColOffset[kNbSMtot];
 for (int i=0;i<kNbSMtot;i++)
    {kTabColOffset[i]=0;
     if (SMdetType[i] == kDCAL) kTabColOffset[i]=kNbColOffsetDCAL;
     }
 
 

printf("\n\n\n\n##############\n###############\n#### Change one thing : when, at the last iteration, the calib coeff = 1, tower is untrusted and product of all calib coeffs of previous iterations must be 1, independantly of their actual value. Check mail a Catherine Apr 18th, 2013 w/ subject 'Calib status'.\n##################\n#################!!!!!!!!!!!!!!!!!!!!!\n\n\n");





 
 
 //*****************************************
 // Define parameters for running the code :
 //*****************************************
 
 //CUSTOMIZE customize :
 isFirstIteration=0; //Not the first iteration.
 testChoice=0; //Not a test.
 //isFirstIteration=1; //This is the first iteration.
 //testChoice=1; //Test with all SMs, 2 col.
 //testChoice=2; //Test with 1 SM, 2 col.
 //testChoice=3; //Test with all SMs but don't plot the inv mass distrib for all the towers (lighter ps file)
 
 //CUSTOMIZE customize : (if necessary) :
 //Calib coeff = mass ratio to the power of coeffPower.
 //coeffPower=1.0; //This choice considers that the observed mass shift is due to an energy shift equally shared by both clusters (ultra-safe but slower convergence).
 coeffPower=1.5; //This choice is an intermediate solution. //Chosen from the analysis of input/pass4/ (previous passes : coeff = 1.0).
 //coeffPower=2.0; //This choice considers that the observed mass shift is entirely due to an energy shift of the considered tower (fastest convergence but may fail when many towers around are all shifted the same way).
 coeffPowerPrev=1.0;
 //coeffPowerPrev=1.5;
 
 //CUSTOMIZE customize : (if necessary) :
 //0 = don't calibrate, 1 = do calibrate.
 kCalibrateMidRapBorderCol=0;
 kCalibrateLastBorderCol=0;
 kCalibrateTouchingBorderRow=0;
 kCalibrateLastBorderRow=0;
 
 paramFitStatUncert_a=18.1227/100.; //This is for the relative uncert : uncert/mean=paramFitStatUncert_a/sqrt(N_pi0) (+) paramFitStatUncert_b. A value of 15./100. means a 15%/sqrt(N) uncertainty.
 paramFitStatUncert_b=0.151497/100.;
 uncertKSpec=0.5; //Want the uncertainty on the mean to be lower than uncertKSpec*uncertSpec : ensure we know the mean, before testing that this mean is close to mPDG.
 uncertSpec=0.02; //Tolerate an additional 2% decalibration with respect to the test-beam status.
 uncertDist=2.; //How much times the tolerated decalibration away from mPDG do we consider that a tower is not yet OK.

 choiceNoEMCAL=0;
 choiceNoDCAL=0;
 for (i=4-1;i>=0;i--)
    {tmpChoice=choice>>i;
     tabChoiceCalos[i]=tmpChoice;
     choice-=tmpChoice*pow(2,i);
     }
 nbTowersConsideredTot=0;
 printf("\n\n---------------------------------\n| Running for :");
 for (i=0;i<4;i++)
    {if (tabChoiceCalos[i] == 1)
        {nbTowersConsideredTot+=kTabNbSM[detTypeType[i]]*kTabNbCol[detTypeType[i]]*kTabNbRow[detTypeType[i]];
         printf(" %s",detTypeString[i]);
         if (i<2) choiceNoEMCAL++;
         else choiceNoDCAL++;
         }
     }
 printf(".\n");
 if (choiceNoEMCAL == 0) printf("| Chose to run with no EMCAL supermodules.\n");
 if (choiceNoDCAL == 0) printf("| Chose to run with no DCAL supermodules.\n");
 if (isFirstIteration == 0) printf("| Declared that this is NOT the first iteration.\n");
    else printf("| Declared that this IS the first iteration.\n");
 printf("---------------------------------\n\n");

 double mPDG = 134.9766; //MeV
 double typicalWidth = 13.; //MeV
 
 
 switch (testChoice)
    {case 0 : kNbColMinEff=0;
              kNbColMaxEff=kNbColMax;
              printf("Calib mode, max tower col span %d -> %d.\n",kNbColMinEff,kNbColMaxEff-1);
              break;
     case 1 : kNbColMinEff=0;
              kNbColMaxEff=2;
              printf("\n########### TEST mode, max tower col span %d -> %d.\n",kNbColMinEff,kNbColMaxEff-1);
              printf("#############################################\n\n");
              gSystem->Exec("sleep 2");
              break;
     case 2 : kNbColMinEff=0;
              kNbColMaxEff=2;
              printf("\n########### TEST mode, 1 SM and max tower col span %d -> %d.\n",kNbColMinEff,kNbColMaxEff-1);
              printf("#############################################\n\n");
              gSystem->Exec("sleep 2");
              break;
     case 3 : kNbColMinEff=0;
              kNbColMaxEff=kNbColMax;
              printf("\n########### TEST mode, max tower col span %d -> %d.\n",kNbColMinEff,kNbColMaxEff-1);
              printf("########### Don't plot the inv mass histoes per tower.\n");
              printf("#############################################\n\n");
              gSystem->Exec("sleep 2");
              break;
     default : printf("\nUnknown test mode option, exiting.\n\n");
               return;
     }
 
 printf("Running with target mass %f MeV.\n",mPDG);
 if (mPDG != 134.9766) printf("\n########### CAUTION ! This is not the PDG mass %f MeV !\n",134.9766);
 printf("\n");
 
 Double_t pi0massP1(Double_t *x, Double_t *par);
 Double_t pi0massP2(Double_t *x, Double_t *par);
 Double_t pi0massP3(Double_t *x, Double_t *par);
 
 int rebin=2;

 




 //*****************************************
 // Define input/output files :
 //*****************************************

 //Current iteration :
 //CUSTOMIZE customize :
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults-fMinNCells2.root","read"); //Test
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults-1_3SMstriggerfMinNCells2.root","read"); //Test
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults-tiersSM.root","read"); //Test
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass0/pi0calib.root","read"); //pass0
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults-LHC15imuoncalo_fMinNCells2WithMaskOnTowers.root","read"); //Test
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResultsLHC15jtot.root","read"); //Test
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults-LHC15jMartinRunsFull.root","read"); //Test
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15iBadMap.root","read"); //
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15iTrigEMC.root","read"); //
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15iMaskTowersByHand.root","read"); //
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15jMaskTowersByHand.root","read"); //
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15iMaskTowersByHand26oct2015.root","read"); //
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15jBadMap27oct2015.root","read"); //
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15jOnlyEdgesMaskTowersByHand28oct2015.root","read"); //
 
 //pass 0 :
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass0/AnalysisResults_LHC15sumijMaskTowersByHand.root","read");
 //pass 1 :
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass1/AnalysisResults_LHC15sumijMaskTowersByHandpass1VeryHighTowers.root","read");
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15sumijBadMappass1VeryHighTowersMinCells3.root","read"); //Test cut minCells at 3.
 //pass 2 :
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass2/AnalysisResults_LHC15sumijMaskTowersByHandpass2VeryHighTowers.root","read"); //Chose that one.
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass2/AnalysisResults_LHC15sumijMaskTowersByHandpass2.root","read");
 //pass 3 :
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass3/AnalysisResults_LHC15sumijMaskTowersByHandpass3.root","read");
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass3/AnalysisResults_LHC15sumijMaskTowersByHandpass3VeryHighTowers.root","read");
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass3/AnalysisResults_LHC15sumij4HV.root","read"); //For HV calculation.
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass3/AnalysisResults_LHC15sumijMaskTowersByHandpass3WithTimeCalibAndAndersValues.root","read"); //Chose that one.
 //pass 4 :
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass4/AnalysisResults_LHC15sumijMaskTowersByHand24Nov2016_pass4_NoNegCoeffs.root","read");
 //pass 5 :
 //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass5/AnalysisResults_LHC15sumijMaskTowersByHand27Nov2016_pass5.root","read");

  //Various studies :
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15sumij2012Calib.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15sumijBadMap2012Calib.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResults_LHC15sumijFakeDeadRegions2012Calib.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/TestLowMassPeakTowers/AnalysisResults_LHC15iBadMapNoMask.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/TestLowMassPeakTowers/AnalysisResults_LHC15iBadMapMask.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/TestLowMassPeakTowers/AnalysisResults_LHC15iNoBadMapNoMask.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/TestLowMassPeakTowers/AnalysisResults_LHC15iNoBadMapMask.root","read");
  //Test fitted mean change wrt sample used :
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/testInvMassPeakMeanVsSampleUsed/LHC15i_Sample1.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/testInvMassPeakMeanVsSampleUsed/LHC15i_Sample2.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/testInvMassPeakMeanVsSampleUsed/LHC15i_Sample2.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/AnalysisResultsLHC15ij_pass3TenderOn.root","read");
  //Test factor 2 :
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/LHC15iTenderOnPass2AllOn.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/LHC15iTenderOffPass2L1PhaseOn.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/LHC15iTenderOnPass2L1PhaseOn.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/studiesWith2015ppData/testInvMassPeakMeanVsSampleUsed/inputSample1/LHC15i_Sample1_TenderOn_muoncalopass2.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/LHC15ijTenderOnTimeCutAppliedNoBadMap.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/LHC15jTenderOnTimeCutAppliedNoBadMapMaskAppliedForAnders.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/studiesWith2015ppData/testAlignment/input/LHC15i_DifferentSM.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/studiesWith2015ppData/testSpecialTriggerEffect/input/LHC15j_SpecialTriggerOnly2by4Structure_SpecialMask.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/studiesWith2015ppData/testBkgShapeEffect/input/LHC15i_TenderNoBadMapTightTimeCutAppliedTightM02Plus1to3_Good_Emax0_7_FixOpAngle_0_112.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/studiesWith2015ppData/testBkgShapeEffect/input/LHC15ij_RemoveLeftShoulder_5MeVShift.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/studiesWith2015ppData/testBkgShapeEffect/input/LHC15ij_Ref.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass3/LHC15sumijTenderOnTimeCutAppliedNoBadMapMaskApplied_WithAndersValues.root","read");
  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass3/AnalysisResults_LHC15sumijMaskTowersByHandpass3WithTimeCalib.root","read");

  //TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass3/AnalysisResults_LHC15iMaskTowersByHandpass3WithTimeCalibSameSM.root","read");

 TFile *f05 = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass5/AnalysisResults_LHC15sumijMaskMoreTowersByHand_pass5.root","read");
    
 FILE *txtFileCalibIn = NULL;
 FILE *txtFileParamsIn = NULL;
 FILE *txtFileParamsOut = NULL;
 //Previous iterations :
 //Pointless ; to change at 1st iteration :
 //CUSTOMIZE customize :
 // Pass 8 2012 EMCAL :
 /*txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0/output2012/pass7/output_calibPi0_coeffs_clean.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0/output2012/pass7/output_calibPi0_parameters.txt","r");*/
 // Pass 0 DCAL :
 /*txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass0_DCAL/output_calibPi0_coeffs_clean.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass0_DCAL/output_calibPi0_parameters.txt","r");*/
 // Pass 0 third-SMs :
 /*txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass0_thirdSMs/output_calibPi0_coeffs_clean.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass0_thirdSMs/output_calibPi0_parameters.txt","r");*/
 // Pass 0 DCAL + third-SMs :
 /*txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass0_mergedDCALandThirdSMs/output_calibPi0_coeffs_clean_veryHighTowers_onlyDCALandThirdSMs.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass0_mergedDCALandThirdSMs/output_calibPi0_parameters_onlyDCALandThirdSMs.txt","r");*/
 // Pass 1 DCAL + third-SMs :
 /*txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass1_DCALandThirdSMs/output_calibPi0_coeffs_clean_veryHighTowers_onlyDCALandThirdSMs.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass1_DCALandThirdSMs/output_calibPi0_parameters.txt","r");*/
 // Pass 2 DCAL + third SMs
 /*txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_coeffs_clean_OnlyDCALandThirdSMs.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_parameters.txt","r");*/ //Chose this one.
 /*txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_coeffs_clean_veryHighGains_OnlyDCALandThirdSMs.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_parameters.txt","r");*/
 // Pass 3 DCal + thirds SM
 /*txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass3_DCALandThirdSMsWithAndersCoeffsTowersVeryFarFromPDG/output_calibPi0_coeffs_clean.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass3_DCALandThirdSMsWithAndersCoeffsTowersVeryFarFromPDG/output_calibPi0_parameters.txt","r");*/
 // Pass 4 DCal + thirds SM EMCAL
 txtFileCalibIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass4_DCALandEMCALThirds/output_calibPi0_coeffs_clean.txt","r");
 txtFileParamsIn = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass4_DCALandEMCALThirds/output_calibPi0_parameters.txt","r");


 TString txtFileParamsOutName("output_calibPi0_parameters.txt");
 TString txtFileCalibOutName("output_calibPi0_coeffs.txt");
 txtFileParamsOut = fopen(txtFileParamsOutName.Data(),"w");
 ofstream txtFileCalibOut(txtFileCalibOutName.Data());






 //*****************************************
 // Define graphical tools :
 //*****************************************

  //defineMyPalette40All(30,5);
 
 Int_t colorIndex[]={1,4,6,8};
 Int_t colorIndexStruct[]={kBlack,kGreen+2,kRed};
 Int_t colorIndexZones[]={kGreen+2,kGreen+3,kAzure-5,kBlue+2,kPink+2,kRed,kRed+2};
 Int_t colorIndexRatiosZones[]={kGreen+2,kBlue,kRed,kRed+2,kBlack};

 TLine *linePDGMass = new TLine(mPDG,0.,mPDG,1.);
 TLine *lineFittedMass = new TLine(mPDG,0.,mPDG,1.);
 TLine *lineFittedWidth1 = new TLine(mPDG,0.,mPDG,1.);
 TLine *lineFittedWidth2 = new TLine(mPDG,0.,mPDG,1.);
 TLine *linePreviousFittedMass = new TLine(mPDG,0.,mPDG,1.);
 TLine *linePreviousFittedWidth1 = new TLine(mPDG,0.,mPDG,1.);
 TLine *linePreviousFittedWidth2 = new TLine(mPDG,0.,mPDG,1.);
 TLine *linePreviousDesiredMass = new TLine(mPDG,0.,mPDG,1.);
 TLine *linePreviousUntrustedFit = new TLine(mPDG,0.,mPDG,1.);
 linePDGMass->SetLineColor(kSpring+4);
 lineFittedMass->SetLineColor(kMagenta);
 lineFittedWidth1->SetLineColor(kMagenta-10);
 lineFittedWidth2->SetLineColor(kMagenta-10);
 TLine *lineSMborderVEMCAL = new TLine(kNbColEMCAL-0.5,-0.5,kNbColEMCAL-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
 TLine *lineSMborderVDCALthird = new TLine(kNbColDCALthird-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5,kNbColDCALthird-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
 TLine *lineSMborderVDCAL1 = new TLine(kNbColDCAL-0.5,-0.5,kNbColDCAL-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5);
 TLine *lineSMborderVDCAL2 = new TLine(2*kNbColEMCAL-kNbColDCAL-0.5,-0.5,2*kNbColEMCAL-kNbColDCAL-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5);
 TLine **lineSMborderHEMCAL,**lineSMborderHDCAL;
 lineSMborderHEMCAL = new TLine*[(int)((kNbSMEMCAL+1)/2)];
 lineSMborderHDCAL = new TLine*[(int)((kNbSMDCAL+1)/2)];
 for (i=0;i<(int)((kNbSMEMCAL+1)/2);i++) lineSMborderHEMCAL[i] = new TLine(-0.5,(i+1)*kNbRowEMCAL-0.5,2.*kNbColEMCAL-0.5,(i+1)*kNbRowEMCAL-0.5);
 for (i=0;i<(int)((kNbSMDCAL+1)/2);i++) lineSMborderHDCAL[i] = new TLine(-0.5,(i+1)*kNbRowDCAL-0.5,2.*kNbColDCALthird-0.5,(i+1)*kNbRowDCAL-0.5);
 TLine *lineMin = new TLine(0.,0.,1.,1.);
 TLine *lineMax = new TLine(0.,0.,1.,1.);
 lineMin->SetLineColor(kGreen+2);
 lineMax->SetLineColor(kGreen+2);
 
 gStyle->SetPaintTextFormat("1.0f"); //For histo hSpaceEntriesDiscard.
 gStyle->SetTextSize(0.5);
 
 //gStyle->SetStatFontSize(0.05); //So as to increase the size of the tower numbers in the stat box. But we don't want to do that for all the histoes in the file, so we prefer the following instead.
 TPaveStats *ptrPaveStat;
 TPave *paveUncertStatus = new TPave(0.95,1.0,0.845,0.850,2,"NDC");
 paveUncertStatus->SetFillStyle(1001);
 TPave *paveMassStatus = new TPave(0.95,1.0,0.845,0.850,2,"NDC");
 paveMassStatus->SetFillStyle(1001);
 TPave **paveCutsStatus;
 paveCutsStatus = new TPave*[kNbFitParams];
 for (i=0;i<kNbFitParams;i++)
    {paveCutsStatus[i] = new TPave(0.95,1.0,0.845,0.850,4,"NDC");
     //paveCutsStatus[i]->SetFillStyle(1001);
     }

 char **varName;
 varName = new char*[kNbFitParams+kNbExtraParams];
 char **varNameLong;
 varNameLong = new char*[kNbFitParams+kNbExtraParams];
 varName[0]="A";
 varName[1]="#mu (MeV)";
 varName[2]="#sigma (MeV)";
 varName[3]="c";
 varName[4]="b";
 varName[5]="a";
 varName[6]="I";
 varName[7]="Khi2/Ndf";
 varName[8]="S";
 varName[9]="I-S";
 varName[10]="Khi2/Ndf peak";
 varName[11]="Khi2/Ndf peakLeft";
 varName[12]="Khi2/Ndf peakRight";
 varNameLong[0]="Amplitude of gaussian fit";
 varNameLong[1]="Mean of gaussian fit";
 varNameLong[2]="Width of gaussian fit";
 varNameLong[3]="Polynom coeff c";
 varNameLong[4]="Polynom coeff b";
 varNameLong[5]="Polynom coeff a";
 varNameLong[6]="Histogram integral";
 varNameLong[7]="Khi2/Ndf";
 varNameLong[8]="Integral of gaussian fit";
 varNameLong[9]="Histo-gaussian integrals";
 varNameLong[10]="Khi2/Ndf around peak";
 varNameLong[11]="Khi2/Ndf half-peak left";
 varNameLong[12]="Khi2/Ndf half-peak right";






 //*****************************************
 // Define plots'ranges and cut values :
 //*****************************************
 
 //CUSTOMIZE customize : (if necessary) :
 //Pass0 on DCAL and third-SMs :
 //Pass1 on DCAL and third-SMs :
 //Pass2 on DCAL and third-SMs :
 //Pass3 on DCAL and third-SMs :
 tabMin[0]=0.;
 tabMin[1]=78.;
 tabMin[2]=1.;
 tabMin[3]=-150.;
 tabMin[4]=-3.0;
 tabMin[5]=-0.014;
 tabMin[6]=0.;
 tabMin[7]=0.;
 tabMin[8]=-100.;
 tabMin[9]=-700.;
 tabMin[10]=0.;
 tabMin[11]=0.;
 tabMin[12]=0.;
 tabMax[0]=1200.;
 tabMax[1]=161.;
 tabMax[2]=36.;
 tabMax[3]=600.;
 tabMax[4]=5.0;
 tabMax[5]=0.0100;
 tabMax[6]=100000.;
 tabMax[7]=30.;
 tabMax[8]=20000.;
 tabMax[9]=80000.;
 tabMax[10]=30.;
 tabMax[11]=7.5;
 tabMax[12]=7.5;
 nbEntriesMax=100000.;

 tabMinDiff[0]=0.;
 tabMinDiff[1]=0.8;
 tabMinDiff[2]=0;
 tabMinDiff[3]=-20.;
 tabMinDiff[4]=-0.7;
 tabMinDiff[5]=-0.004;
 tabMinDiff[6]=0.7;
 tabMinDiff[7]=0.7;
 tabMinDiff[8]=0.;
 tabMaxDiff[0]=4.;
 tabMaxDiff[1]=1.2;
 tabMaxDiff[2]=4.;
 tabMaxDiff[3]=25.;
 tabMaxDiff[4]=1.;
 tabMaxDiff[5]=0.004;
 tabMaxDiff[6]=1.8;
 tabMaxDiff[7]=1.9;
 tabMaxDiff[8]=4.;

 cutEntriesDiscard=500.; //It's impossible to have a decent pi0 peak with less than 500 entries in the histo. Shy peaks start to appear around 800 entries, very few are found below 1000.
 cutAmpDiscard=8.001;
 
 //CUSTOMIZE customize :
  
 //Pass0 on DCAL :
 /*cutMin[0]=10.;
 cutMin[0]=70.; //Check made later.
 cutMin[1]=128.;
 cutMin[2]=8.0; //Default = 7.0.
 cutMin[3]=-50.;
 cutMin[4]=-1.00;
 cutMin[5]=-0.008;
 cutMin[6]=500.;
 cutMin[7]=0.;
 cutMin[8]=200.;
 cutMin[8]=1500.; //Check made later.
 cutMin[9]=500.;
 cutMax[0]=950.;
 cutMax[1]=153.;
 cutMax[2]=24.;
 cutMax[3]=400.;
 cutMax[4]=2.3;
 cutMax[5]=0.0022;
 cutMax[6]=43000.;
 cutMax[7]=3.6;
 cutMax[8]=13500.;
 cutMax[9]=40000.;*/

 //Pass0 (w/ time cuts) on DCAL :
// cutMin[0]=40.;
// cutMin[1]=130.;
// cutMin[2]=8.0;
// cutMin[3]=-25.;
// cutMin[4]=-1.20;
// cutMin[5]=-0.006;
// cutMin[6]=1000.;
// cutMin[7]=0.3;
// cutMin[8]=450.;
// cutMin[9]=250.;
// cutMin[10]=0.15;
// cutMin[11]=0.15;
// cutMin[12]=0.15;
// cutMax[0]=1000.;
// cutMax[1]=140.;
// cutMax[2]=20.;
// cutMax[3]=250.;
// cutMax[4]=2.0;
// cutMax[5]=0.0022;
// cutMax[6]=32000.;
// cutMax[7]=3.6;
// cutMax[8]=12000.;
// cutMax[9]=21000.;
// cutMax[10]=4.7;
// cutMax[11]=9.;
// cutMax[12]=4.0;
    
//Pass0 (w/ time cuts) on EMCAL :
 /*cutMin[0]=40.;
 cutMin[1]=124.;
 cutMin[2]=8.0;
 cutMin[3]=-25.;
 cutMin[4]=-1.0;
 cutMin[5]=-0.004;
 cutMin[6]=1000.;
 cutMin[7]=0.3;
 cutMin[8]=400.;
 cutMin[9]=250.;
 cutMin[10]=0.15;
 cutMin[11]=0.15;
 cutMin[12]=0.15;
 cutMax[0]=800.;
 cutMax[1]=142.;
 cutMax[2]=20.;
 cutMax[3]=200.;
 cutMax[4]=1.5;
 cutMax[5]=0.0022;
 cutMax[6]=30000.;
 cutMax[7]=3.6;
 cutMax[8]=10000.;
 cutMax[9]=18000.;
 cutMax[10]=4.7;
 cutMax[11]=9.;
 cutMax[12]=4.0;*/

 //Pass0 on third-SMs :
 /*cutMin[0]=10.;
 cutMin[1]=133.; //Default = 128.
 cutMin[2]=11.0; //Default = 7.0.
 cutMin[3]=-50.;
 cutMin[4]=-1.20;
 cutMin[5]=-0.008;
 cutMin[6]=500.;
 cutMin[7]=0.;
 cutMin[8]=200.;
 cutMin[9]=500.;
 cutMax[0]=500.;
 cutMax[1]=155.;
 cutMax[2]=32.;
 cutMax[3]=220.;
 cutMax[4]=2.3;
 cutMax[5]=0.0022;
 cutMax[6]=30000.;
 cutMax[7]=3.6;
 cutMax[8]=8500.;
 cutMax[9]=35000.;*/

 //Pass1 on DCAL and third-SMs :
 /*cutMin[0]=40.;
 cutMin[1]=128.;
 cutMin[2]=7.0;
 cutMin[3]=-50.;
 cutMin[4]=-1.20;
 cutMin[5]=-0.006;
 cutMin[6]=2000.;
 cutMin[7]=0.3;
 cutMin[8]=900.;
 cutMin[9]=500.;
 cutMax[0]=1030.;
 cutMax[1]=142.;
 cutMax[2]=20.;
 cutMax[3]=300.;
 cutMax[4]=2.3;
 cutMax[5]=0.0022;
 cutMax[6]=38000.;
 cutMax[7]=3.6;
 cutMax[8]=13000.;
 cutMax[9]=25000.;*/

 //Pass2 on DCAL and third-SMs :
 /*cutMin[0]=20.;
 cutMin[1]=128.;
 cutMin[2]=7.0;
 cutMin[3]=-50.;
 cutMin[4]=-1.20;
 cutMin[5]=-0.006;
 cutMin[6]=1000.;
 cutMin[7]=0.3;
 cutMin[8]=450.;
 cutMin[9]=250.;
 cutMax[0]=1030.;
 cutMax[1]=142.;
 cutMax[2]=20.;
 cutMax[3]=300.;
 cutMax[4]=2.3;
 cutMax[5]=0.0022;
 cutMax[6]=38000.;
 cutMax[7]=3.6;
 cutMax[8]=13000.;
 cutMax[9]=25000.;*/

 //Pass3 on DCAL and third-SMs :
 /*cutMin[0]=40.;
 cutMin[1]=128.;
 cutMin[2]=7.0;
 cutMin[3]=-25.;
 cutMin[4]=-1.20;
 cutMin[5]=-0.006;
 cutMin[6]=1000.;
 cutMin[7]=0.3;
 cutMin[8]=450.;
 cutMin[9]=250.;
 cutMin[10]=0.15;
 cutMin[11]=0.15;
 cutMin[12]=0.15;
 cutMax[0]=900.;
 cutMax[1]=142.;
 cutMax[2]=20.;
 cutMax[3]=230.;
 cutMax[4]=1.8;
 cutMax[5]=0.0022;
 cutMax[6]=32000.;
 cutMax[7]=3.6;
 cutMax[8]=10800.;
 cutMax[9]=20000.;
 cutMax[10]=4.7;
 cutMax[11]=9.;
 cutMax[12]=4.0;*/

 //Pass4 on DCAL and EMCAL third-SMs :
 /*cutMin[0]=40.;
 cutMin[1]=128.;
 cutMin[2]=8.0;
 cutMin[3]=-25.;
 cutMin[4]=-1.20;
 cutMin[5]=-0.006;
 cutMin[6]=1000.;
 cutMin[7]=0.3;
 cutMin[8]=450.;
 cutMin[9]=250.;
 cutMin[10]=0.15;
 cutMin[11]=0.15;
 cutMin[12]=0.15;
 cutMax[0]=900.;
 cutMax[1]=142.;
 cutMax[2]=20.;
 cutMax[3]=230.;
 cutMax[4]=1.8;
 cutMax[5]=0.0022;
 cutMax[6]=32000.;
 cutMax[7]=3.6;
 cutMax[8]=11000.;
 cutMax[9]=20000.;
 cutMax[10]=4.7;
 cutMax[11]=9.;
 cutMax[12]=4.0;*/

 //Pass5 (last diagnostic) on DCAL and EMCAL third-SMs :
 cutMin[0]=40.;
 cutMin[1]=130.;
 cutMin[2]=8.0;
 cutMin[3]=-25.;
 cutMin[4]=-1.20;
 cutMin[5]=-0.006;
 cutMin[6]=1000.;
 cutMin[7]=0.3;
 cutMin[8]=450.;
 cutMin[9]=250.;
 cutMin[10]=0.15;
 cutMin[11]=0.15;
 cutMin[12]=0.15;
 cutMax[0]=1000.;
 cutMax[1]=140.;
 cutMax[2]=20.;
 cutMax[3]=250.;
 cutMax[4]=2.0;
 cutMax[5]=0.0022;
 cutMax[6]=32000.;
 cutMax[7]=3.6;
 cutMax[8]=12000.;
 cutMax[9]=21000.;
 cutMax[10]=4.7;
 cutMax[11]=9.;
 cutMax[12]=4.0;

    
 //LHC15i for special trigger study :
 /*cutMin[0]=20.;
 cutMin[1]=128.;
 cutMin[2]=7.0;
 cutMin[3]=-40.;
 cutMin[4]=-0.4;
 cutMin[5]=-0.003;
 cutMin[6]=1000.;
 cutMin[7]=0.3;
 cutMin[8]=450.;
 cutMin[9]=250.;
 cutMin[10]=0.;
 cutMax[0]=400.;
 cutMax[1]=142.;
 cutMax[2]=20.;
 cutMax[3]=100.;
 cutMax[4]=1.;
 cutMax[5]=0.001;
 cutMax[6]=20000.;
 cutMax[7]=3.6;
 cutMax[8]=5000.;
 cutMax[9]=15000.;
 cutMax[10]=3.5;*/
    
 //LHC15j Special trigger for special trigger study :
 /*cutMin[0]=20.;
 cutMin[1]=128.;
 cutMin[2]=7.0;
 cutMin[3]=-40.;
 cutMin[4]=-0.4;
 cutMin[5]=-0.002;
 cutMin[6]=1000.;
 cutMin[7]=0.3;
 cutMin[8]=450.;
 cutMin[9]=250.;
 cutMin[10]=0.;
 cutMax[0]=100.;
 cutMax[1]=142.;
 cutMax[2]=20.;
 cutMax[3]=60.;
 cutMax[4]=0.8;
 cutMax[5]=0.001;
 cutMax[6]=10000.;
 cutMax[7]=3.6;
 cutMax[8]=3000.;
 cutMax[9]=10000.;
 cutMax[10]=3.5;*/

 //LHC15i for left shoulder :
 /*cutMin[0]=20.;
 cutMin[1]=132.;
 cutMin[2]=9.0;
 cutMin[3]=-50.;
 cutMin[4]=-1.20;
 cutMin[5]=-0.006;
 cutMin[6]=1000.;
 cutMin[7]=0.3;
 cutMin[8]=450.;
 cutMin[9]=250.;
 cutMin[10]=0.;
 cutMax[0]=1030.; //changed
 cutMax[1]=138.;
 cutMax[2]=20.;
 cutMax[3]=300.;  //changed
 cutMax[4]=2.3;
 cutMax[5]=0.0022;
 cutMax[6]=38000.;
 cutMax[7]=3.6;
 cutMax[8]=13000.;
 cutMax[9]=25000.;
 cutMax[10]=2.0;*/



 //*****************************************
 // Define pT bins :
 //*****************************************

 //CUSTOMIZE customize : (if necessary) :
 int nbPtBins;
 nbPtBins=21; //Number of bins in pT
 
 //Table with borders between pT bins :
 double *tabPtBins;
 tabPtBins = new double[nbPtBins-1];
 tabPtBins[0]=1.0;//ok pour energy cut = 0.5
 tabPtBins[1]=1.2;
 tabPtBins[2]=1.4;
 tabPtBins[3]=1.6;
 tabPtBins[4]=1.8;
 tabPtBins[5]=2.0;
 tabPtBins[6]=2.2;
 tabPtBins[7]=2.4;
 tabPtBins[8]=2.6;
 tabPtBins[9]=2.8;
 tabPtBins[10]=3.0;
 tabPtBins[11]=3.3;
 tabPtBins[12]=3.6;
 tabPtBins[13]=3.9;
 tabPtBins[14]=4.2;
 tabPtBins[15]=4.6;
 tabPtBins[16]=5.0;
 tabPtBins[17]=5.4;
 tabPtBins[18]=5.8;
 tabPtBins[19]=6.2;
 
 





 //*****************************************
 // Get histoes from input file :
 //*****************************************

 //CUSTOMIZE customize : (if necessary) :
 //TList *l05 = (TList *) f05->Get("Pi0Calibration_TrigCEMC7");
 //TList *l05 = (TList *) f05->Get("Pi0Calibration_TrigCINT");
 //TList *l05 = (TList *) f05->Get("pi0calib");
 //TList *l05 = (TList *) f05->Get("Pi0Calibration_TrigCEMC");
 TList *l05 = (TList *) f05->Get("Pi0Calibration_Trig");
 //TList *l05 = (TList *) f05->Get("Pi0Calibration_TrigEMC");
 if (l05 == NULL)
    {printf("\n\nNull pointer to TList.\n");
     printf("Try with this other name instead :\n");
     return;
     }
 
 const int kNbZones=7;
 TH2F **hAllM_05_SM,**hAllM_05_SM_masked,***hAllM_05_SM_Zones,**hAllM_05_Zones;
 hAllM_05_SM = new TH2F*[kNbSMtot];
 hAllM_05_SM_masked = new TH2F*[kNbSMtot];
 hAllM_05_SM_Zones = new TH2F**[kNbSMtot];
 hAllM_05_Zones = new TH2F*[kNbZones];
 TH2F *hAllM_05      = (TH2F *) l05->FindObject("hmgg");
 TH2F *hAllM_05_masked = (TH2F *)hAllM_05->Clone("hmgg_masked");
 for (j=0;j<kNbZones;j++)
    {hAllM_05_Zones[j] = (TH2F *)hAllM_05->Clone(Form("hmgg_Zone%d",j+1));
     }
 for (i=0;i<kNbSMtot;i++)
    {hAllM_05_SM[i]      = (TH2F *) l05->FindObject(Form("hmgg_SM%d",i));
     hAllM_05_SM_masked[i]      = (TH2F *) l05->FindObject(Form("hmgg_SM%d_MaskFrame",i));
     hAllM_05_SM_Zones[i] = new TH2F*[kNbZones];
     for (j=0;j<kNbZones;j++)
        {hAllM_05_SM_Zones[i][j] = (TH2F *) l05->FindObject(Form("hmgg_SM%d_Zone%d",i,j+1));
         }
     }
 /*hAllM_05_SM[0]      = (TH2F *) l05->FindObject("hmgg_SM0");
 hAllM_05_SM[1]      = (TH2F *) l05->FindObject("hmgg_SM1");
 hAllM_05_SM[2]      = (TH2F *) l05->FindObject("hmgg_SM2");
 hAllM_05_SM[3]      = (TH2F *) l05->FindObject("hmgg_SM3");
 hAllM_05_SM[4]      = (TH2F *) l05->FindObject("hmgg_SM4");
 hAllM_05_SM[5]      = (TH2F *) l05->FindObject("hmgg_SM5");
 hAllM_05_SM[6]      = (TH2F *) l05->FindObject("hmgg_SM6");
 hAllM_05_SM[7]      = (TH2F *) l05->FindObject("hmgg_SM7");
 hAllM_05_SM[8]      = (TH2F *) l05->FindObject("hmgg_SM8");
 hAllM_05_SM[9]      = (TH2F *) l05->FindObject("hmgg_SM9");
 hAllM_05_SM_masked[0]      = (TH2F *) l05->FindObject("hmgg_SM0_MaskFrame");
 hAllM_05_SM_masked[1]      = (TH2F *) l05->FindObject("hmgg_SM1_MaskFrame");
 hAllM_05_SM_masked[2]      = (TH2F *) l05->FindObject("hmgg_SM2_MaskFrame");
 hAllM_05_SM_masked[3]      = (TH2F *) l05->FindObject("hmgg_SM3_MaskFrame");
 hAllM_05_SM_masked[4]      = (TH2F *) l05->FindObject("hmgg_SM4_MaskFrame");
 hAllM_05_SM_masked[5]      = (TH2F *) l05->FindObject("hmgg_SM5_MaskFrame");
 hAllM_05_SM_masked[6]      = (TH2F *) l05->FindObject("hmgg_SM6_MaskFrame");
 hAllM_05_SM_masked[7]      = (TH2F *) l05->FindObject("hmgg_SM7_MaskFrame");
 hAllM_05_SM_masked[8]      = (TH2F *) l05->FindObject("hmgg_SM8_MaskFrame");
 hAllM_05_SM_masked[9]      = (TH2F *) l05->FindObject("hmgg_SM9_MaskFrame");*/

 TH1F *hNevts = (TH1F*)l05->FindObject("hNEvents");
 printf("\n\nNumber of events in histo hNEvents : %d, nb of entries %d (%d Mevts).\n\n",(int)hNevts->GetBinContent(1),(int)hNevts->GetEntries(),(int)(hNevts->GetEntries()/1000000.));

 hAllM_05->Reset();
 hAllM_05_masked->Reset();
 for (j=0;j<kNbZones;j++)
    {hAllM_05_Zones[j]->Reset();
     }
 for (i=0;i<kNbSMtot;i++)
    {if (tabChoiceCalos[SMdetType[i]] == 1)
        {hAllM_05->Add(hAllM_05_SM[i]);
         hAllM_05_masked->Add(hAllM_05_SM_masked[i]);
         for (j=0;j<kNbZones;j++)
            {hAllM_05_Zones[j]->Add(hAllM_05_SM_Zones[i][j]);
             }
         }
     }
 




 

 //*****************************************
 // Instanciate pT histoes and graphs :
 //*****************************************

 TH1F **hAllMPerSM,**hMassDistrPerSMPerPtPerZone;
 hAllMPerSM = new TH1F*[3*(nbPtBins+1)*(kNbSMtot+1)];
 hMassDistrPerSMPerPtPerZone = new TH1F*[kNbZones*(nbPtBins+1)*(kNbSMtot+1)];
 // pT bin "nbPtBins" : pT-integrated.    SM bin "kNbSMtot" : integrated over SMs.
 for (i=0;i<(kNbSMtot+1);i++)
    {for (j=0;j<(nbPtBins+1);j++)
        {for (k=0;k<3;k++)
            {hAllMPerSM[(3*(nbPtBins+1))*i+3*j+k] = new TH1F(Form("hAllMPerSM%d_pT%d_%d",i,j,k),Form("hAllMPerSM%d_pT%d_%d",i,j,k),hAllM_05->GetNbinsX(),hAllM_05->GetXaxis()->GetXmin(),hAllM_05->GetXaxis()->GetXmax());
             hAllMPerSM[(3*(nbPtBins+1))*i+3*j+k]->SetXTitle("#pi^{0} inv mass (MeV)");
             hAllMPerSM[(3*(nbPtBins+1))*i+3*j+k]->SetYTitle("Counts");
             //hAllMPerSM[(3*(nbPtBins+1))*i+3*j+k]->SetStats(0);
             hAllMPerSM[(3*(nbPtBins+1))*i+3*j+k]->SetLineColor(colorIndexStruct[k]);
             }
         for (k=0;k<kNbZones;k++)
            {hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*i+kNbZones*j+k] = new TH1F(Form("hMassDistrPerSMPerPtPerZone_SM%d_pT%d_Zone%d",i,j,k),Form("hMassDistrPerSMPerPtPerZone%d_pT%d_Zone%d",i,j,k),hAllM_05->GetNbinsX(),hAllM_05->GetXaxis()->GetXmin(),hAllM_05->GetXaxis()->GetXmax());
             hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*i+kNbZones*j+k]->SetXTitle("#pi^{0} inv mass (MeV)");
             hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*i+kNbZones*j+k]->SetYTitle("Counts");
             hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*i+kNbZones*j+k]->SetLineColor(colorIndexZones[k]);
             }
         }
     }
 
 double *tabGraphX,*tabGraphYmass,*tabGraphYsig,*tabGraphYnumber,*tabGraphYnumberRatio,*tabErrGraphX,*tabErrGraphYmass,*tabErrGraphYsig,*tabErrGraphYnumber,*tabErrGraphYnumberRatio;
 double *tabTmpGraphYmass,*tabTmpGraphYsig,*tabTmpGraphYnumber,*tabTmpGraphYnumberRatio,*tabTmpErrGraphYmass,*tabTmpErrGraphYsig,*tabTmpErrGraphYnumber,*tabTmpErrGraphYnumberRatio;
 double cmptPtMoyLow,cmptPtMoyHigh;
 double pTmoyLow,pTmoyHigh;
 
 tabGraphX = new double[nbPtBins];
 tabGraphYmass = new double[kNbZones*nbPtBins];
 tabGraphYsig = new double[kNbZones*nbPtBins];
 tabGraphYnumber = new double[kNbZones*nbPtBins];
 tabGraphYnumberRatio = new double[((int)(kNbZones/2))*(kNbZones-1)*nbPtBins];
 tabErrGraphX = new double[nbPtBins];
 tabErrGraphYmass = new double[kNbZones*nbPtBins];
 tabErrGraphYsig = new double[kNbZones*nbPtBins];
 tabErrGraphYnumber = new double[kNbZones*nbPtBins];
 tabErrGraphYnumberRatio = new double[((int)(kNbZones/2))*(kNbZones-1)*nbPtBins];
 tabTmpGraphYmass = new double[nbPtBins];
 tabTmpGraphYsig = new double[nbPtBins];
 tabTmpGraphYnumber = new double[nbPtBins];
 tabTmpGraphYnumberRatio = new double[nbPtBins];
 tabTmpErrGraphYmass = new double[nbPtBins];
 tabTmpErrGraphYsig = new double[nbPtBins];
 tabTmpErrGraphYnumber = new double[nbPtBins];
 tabTmpErrGraphYnumberRatio = new double[nbPtBins];
 
 TGraphErrors **tgeMass,**tgeSig,**tgeNumber,**tgeNumberRatio,**tgeMassZones,**tgeSigZones,**tgeNumberZones,**tgeNumberRatioZones;
 tgeMass = new TGraphErrors*[3*(kNbSMtot+1)];
 tgeSig = new TGraphErrors*[3*(kNbSMtot+1)];
 tgeNumber = new TGraphErrors*[3*(kNbSMtot+1)];
 tgeNumberRatio = new TGraphErrors*[3*(kNbSMtot+1)];
 tgeMassZones = new TGraphErrors*[kNbZones*(kNbSMtot+1)];
 tgeSigZones = new TGraphErrors*[kNbZones*(kNbSMtot+1)];
 tgeNumberZones = new TGraphErrors*[kNbZones*(kNbSMtot+1)];
 tgeNumberRatioZones = new TGraphErrors*[((int)(kNbZones/2))*(kNbZones-1)*(kNbSMtot+1)];




 
 
 //*****************************************
 // Instanciate other histoes :
 //*****************************************

 TH2F **hSpace;
 TH2F **hSpaceIntg;
 TH2F **hSpaceCuts;
 TH2F **hSpaceCoeff;
 TH2F **hSpaceDiff;
 TH2F **hSpaceFitErr;
 TH2F **hSpaceEntriesDiscard;
 TH2F **hSpaceEntriesDiscardHack;
 hSpace = new TH2F*[kNbFitParams*kNbSMtot];
 hSpaceIntg = new TH2F*[kNbExtraParams*kNbSMtot];
 hSpaceCuts = new TH2F*[5*kNbSMtot];
 hSpaceCoeff = new TH2F*[1*kNbSMtot];
 hSpaceDiff = new TH2F*[kNbTotParams*kNbSMtot];
 hSpaceFitErr = new TH2F*[1*kNbSMtot];
 hSpaceEntriesDiscard = new TH2F*[1*kNbSMtot];
 hSpaceEntriesDiscardHack = new TH2F*[1*kNbSMtot];
 for (i=0;i<kNbSMtot;i++)
    {for (j=0;j<kNbFitParams;j++)
        {hSpace[kNbFitParams*i+j] = new TH2F(Form("hSpace_SM%d_%d",i,j),Form("hSpace_SM%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpace[kNbFitParams*i+j]->SetXTitle("Column");
         hSpace[kNbFitParams*i+j]->SetYTitle("Row");
         hSpace[kNbFitParams*i+j]->SetStats(0);
         hSpace[kNbFitParams*i+j]->SetContour(30);
         }
     for (j=0;j<kNbExtraParams;j++)
        {hSpaceIntg[kNbExtraParams*i+j] = new TH2F(Form("hSpaceIntg_SM%d_%d",i,j),Form("hSpaceIntg_SM%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpaceIntg[kNbExtraParams*i+j]->SetXTitle("Column");
         hSpaceIntg[kNbExtraParams*i+j]->SetYTitle("Row");
         hSpaceIntg[kNbExtraParams*i+j]->SetStats(0);
         hSpaceIntg[kNbExtraParams*i+j]->SetContour(30);
         }
     for (j=0;j<5;j++)
        {hSpaceCuts[5*i+j] = new TH2F(Form("hSpaceCuts_SM%d_%d",i,j),Form("hSpaceCuts_SM%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpaceCuts[5*i+j]->SetXTitle("Column");
         hSpaceCuts[5*i+j]->SetYTitle("Row");
         hSpaceCuts[5*i+j]->SetStats(0);
         hSpaceCuts[5*i+j]->SetContour(30);
         }
     for (j=0;j<1;j++)
        {hSpaceCoeff[1*i+j] = new TH2F(Form("hSpaceCoeff_SM%d_%d",i,j),Form("hSpaceCoeff_SM%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpaceCoeff[1*i+j]->SetXTitle("Column");
         hSpaceCoeff[1*i+j]->SetYTitle("Row");
         hSpaceCoeff[1*i+j]->SetStats(0);
         hSpaceCoeff[1*i+j]->SetContour(30);
         hSpaceFitErr[1*i+j] = new TH2F(Form("hSpaceFitErr_SM%d_%d",i,j),Form("hSpaceFitErr_SM%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpaceFitErr[1*i+j]->SetXTitle("Column");
         hSpaceFitErr[1*i+j]->SetYTitle("Row");
         hSpaceFitErr[1*i+j]->SetStats(0);
         hSpaceFitErr[1*i+j]->SetContour(30);
         hSpaceEntriesDiscard[1*i+j] = new TH2F(Form("hSpaceEntriesDiscard_SM%d_%d",i,j),Form("hSpaceEntriesDiscard_SM%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpaceEntriesDiscard[1*i+j]->SetXTitle("Column");
         hSpaceEntriesDiscard[1*i+j]->SetYTitle("Row");
         hSpaceEntriesDiscard[1*i+j]->SetStats(0);
         hSpaceEntriesDiscard[1*i+j]->SetContour(30);
         hSpaceEntriesDiscardHack[1*i+j] = new TH2F(Form("hSpaceEntriesDiscardHack_SM%d_%d",i,j),Form("hSpaceEntriesDiscardHack_SM%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpaceEntriesDiscardHack[1*i+j]->SetStats(0);
         hSpaceEntriesDiscardHack[1*i+j]->Fill(-1,5.); //Plot must have an entry to plot some text, thus to change the text size ! But putting the entry out of the drawable range also works.
         hSpaceEntriesDiscardHack[1*i+j]->SetMarkerSize(1.5);
         }
     for (j=0;j<kNbTotParams;j++)
        {hSpaceDiff[kNbTotParams*i+j] = new TH2F(Form("hSpaceDiff_SM%d_%d",i,j),Form("hSpaceDiff_SM%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpaceDiff[kNbTotParams*i+j]->SetXTitle("Column");
         hSpaceDiff[kNbTotParams*i+j]->SetYTitle("Row");
         hSpaceDiff[kNbTotParams*i+j]->SetStats(0);
         hSpaceDiff[kNbTotParams*i+j]->SetContour(30);
         }
     }
 TH1F *hAllDistribNbEntries;
 hAllDistribNbEntries = new TH1F("hAllDistribNbEntries","hAllDistribNbEntries",100,0.,nbEntriesMax);
 hAllDistribNbEntries->SetXTitle("Nb of entries");
 hAllDistribNbEntries->SetYTitle("Counts");
 hAllDistribNbEntries->SetStats(0);
 TH2F *hAllDistribNbEntriesCorrel;
 hAllDistribNbEntriesCorrel = new TH2F("hAllDistribNbEntriesCorrel","hAllDistribNbEntriesCorrel",100,tabMin[6],tabMax[6]/4.,100,0.,nbEntriesMax/4.);
 hAllDistribNbEntriesCorrel->SetXTitle("I");
 hAllDistribNbEntriesCorrel->SetYTitle("Nb of entries");
 hAllDistribNbEntriesCorrel->SetStats(0);
 TH1F **hAllDistrib;
 TH1F **hAllDistribDistance;
 TH1F **hAllDistribDistance2;
 TH1F *hAllDistribMass;
 TH1F *hAllOldDistribMass;
 TH1F *hMassPerSM;
 TH1F *hMassOldPerSM;
 TH1F **hAllDistribMassPerSM;
 TH1F **hAllOldDistribMassPerSM;
 TH1F **hAllDistribIntg;
 TH2F **hAllSpaceEMCAL;
 TH2F **hAllSpaceEMCALIntg;
 TH2F **hAllSpaceEMCALCoeff;
 TH2F **hAllSpaceEMCALDistance;
 TH2F **hAllSpaceEMCALDiff;
 TH2F **hAllSpaceDCAL;
 TH2F **hAllSpaceDCALIntg;
 TH2F **hAllSpaceDCALCoeff;
 TH2F **hAllSpaceDCALDistance;
 TH2F **hAllSpaceDCALDiff;
 TH1F **hAllOldDistrib;
 TH1F **hAllDiff;
 TH1F **hAllDiffAllTw;
 hAllDistribMassPerSM = new TH1F*[kNbSMtot];
 hAllOldDistribMassPerSM = new TH1F*[kNbSMtot];
 hAllDistrib = new TH1F*[kNbFitParams];
 hAllDistribDistance = new TH1F*[3];
 hAllDistribDistance2 = new TH1F*[2*2];
 hAllDistribIntg = new TH1F*[kNbExtraParams];
 hAllSpaceEMCAL = new TH2F*[kNbFitParams];
 hAllSpaceEMCALIntg = new TH2F*[kNbExtraParams];
 hAllSpaceEMCALCoeff = new TH2F*[1];
 hAllSpaceEMCALDistance = new TH2F*[3];
 hAllSpaceEMCALDiff = new TH2F*[kNbTotParams];
 hAllSpaceDCAL = new TH2F*[kNbFitParams];
 hAllSpaceDCALIntg = new TH2F*[kNbExtraParams];
 hAllSpaceDCALCoeff = new TH2F*[1];
 hAllSpaceDCALDistance = new TH2F*[3];
 hAllSpaceDCALDiff = new TH2F*[kNbTotParams];
 hAllOldDistrib = new TH1F*[kNbFitParams+kNbExtraParams];
 hAllDiff = new TH1F*[kNbTotParams];
 hAllDiffAllTw = new TH1F*[kNbTotParams];

 
 for (j=0;j<kNbSMtot;j++)
    {hAllDistribMassPerSM[j] = new TH1F(Form("hAllDistribMassPerSM_%d",j),Form("hAllDistribMassPerSM_%d",j),150,mPDG-10.,mPDG+10.);
     hAllDistribMassPerSM[j]->SetXTitle(Form("#pi^{0} inv mass (MeV) SM %d",j));
     hAllDistribMassPerSM[j]->SetYTitle("Counts");
     hAllDistribMassPerSM[j]->SetStats(0);
     hAllOldDistribMassPerSM[j] = new TH1F(Form("hAllOldDistribMassPerSM%d",j),Form("hAllOldDistribMassPerSM%d",j),150,mPDG-10.,mPDG+10.);
     hAllOldDistribMassPerSM[j]->SetXTitle(Form("#pi^{0} inv mass (MeV) SM %d",j));
     hAllOldDistribMassPerSM[j]->SetYTitle("Counts");
     hAllOldDistribMassPerSM[j]->SetStats(0);
     }
 for (j=0;j<kNbFitParams;j++)
    {hAllSpaceEMCAL[j] = new TH2F(Form("hAllSpaceEMCAL_%d",j),Form("hAllSpaceEMCAL_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCAL[j]->SetXTitle("Column");
     hAllSpaceEMCAL[j]->SetYTitle("Row");
     hAllSpaceEMCAL[j]->SetStats(0);
     hAllSpaceEMCAL[j]->SetContour(30);
     hAllSpaceDCAL[j] = new TH2F(Form("hAllSpaceDCAL_%d",j),Form("hAllSpaceDCAL_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCAL[j]->SetXTitle("Column");
     hAllSpaceDCAL[j]->SetYTitle("Row");
     hAllSpaceDCAL[j]->SetStats(0);
     hAllSpaceDCAL[j]->SetContour(30);
     hAllDistrib[j] = new TH1F(Form("hAllDistrib_%d",j),Form("hAllDistrib_%d",j),100,tabMin[j],tabMax[j]);
     hAllDistrib[j]->SetXTitle(Form("%s",varName[j]));
     hAllDistrib[j]->SetYTitle("Counts");
     hAllDistrib[j]->SetStats(0);
     }
 for (j=0;j<3;j++)
    {hAllDistribDistance[j] = new TH1F(Form("hAllDistribDistance_%d",j),Form("hAllDistribDistance_%d",j),100,0.,10.);
     hAllDistribDistance[j]->SetXTitle("(#mu_{fit}-#mu_{PDG})/#sigma");
     hAllDistribDistance[j]->SetYTitle("Counts");
     hAllDistribDistance[j]->SetStats(0);
     }
 for (j=0;j<2;j++)
    {hAllDistribDistance2[2*j+0] = new TH1F(Form("hAllDistribDistance2_%d_0",j),Form("hAllDistribDistance2_%d_0",j),100,0.,5.);
     hAllDistribDistance2[2*j+0]->SetXTitle("#sigma_{#mu} (stat) (%)");
     hAllDistribDistance2[2*j+0]->SetYTitle("Counts");
     hAllDistribDistance2[2*j+0]->SetStats(0);
     hAllDistribDistance2[2*j+1] = new TH1F(Form("hAllDistribDistance2_%d_1",j),Form("hAllDistribDistance2_%d_1",j),100,0.,0.2);
     hAllDistribDistance2[2*j+1]->SetXTitle("|#mu_{fit}-#mu_{PDG}|/#mu_{PDG}");
     hAllDistribDistance2[2*j+1]->SetYTitle("Counts");
     hAllDistribDistance2[2*j+1]->SetStats(0);
     }
 hAllDistribMass = new TH1F("hAllDistribMass","hAllDistribMass",200,mPDG-25.,mPDG+25.);
 hAllDistribMass->SetXTitle("#pi^{0} inv mass (MeV)");
 hAllDistribMass->SetYTitle("Counts");
 hAllDistribMass->SetStats(0);
 hAllOldDistribMass = new TH1F("hAllOldDistribMass","hAllOldDistribMass",200,mPDG-25.,mPDG+25.);
 hAllOldDistribMass->SetXTitle("#pi^{0} inv mass (MeV)");
 hAllOldDistribMass->SetYTitle("Counts");
 hAllOldDistribMass->SetStats(0);
 hMassPerSM = new TH1F("hMassPerSM","hMassPerSM",kNbSMtot,-0.5,kNbSMtot-0.5);
 hMassPerSM->SetXTitle("SM number");
 hMassPerSM->SetYTitle("#pi^{0} inv mass (MeV)");
 hMassPerSM->SetStats(0);
 hMassOldPerSM = new TH1F("hMassOldPerSM","hMassOldPerSM",kNbSMtot,-0.42,kNbSMtot-0.42); //Move a bit so as to plot together with hMassPerSM.
 hMassOldPerSM->SetXTitle("SM number");
 hMassOldPerSM->SetYTitle("#pi^{0} inv mass (MeV)");
 hMassOldPerSM->SetStats(0);
 for (j=0;j<kNbExtraParams;j++)
    {hAllSpaceEMCALIntg[j] = new TH2F(Form("hAllSpaceEMCALIntg_%d",j),Form("hAllSpaceEMCALIntg_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALIntg[j]->SetXTitle("Column");
     hAllSpaceEMCALIntg[j]->SetYTitle("Row");
     hAllSpaceEMCALIntg[j]->SetStats(0);
     hAllSpaceEMCALIntg[j]->SetContour(30);
     hAllSpaceDCALIntg[j] = new TH2F(Form("hAllSpaceDCALIntg_%d",j),Form("hAllSpaceDCALIntg_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALIntg[j]->SetXTitle("Column");
     hAllSpaceDCALIntg[j]->SetYTitle("Row");
     hAllSpaceDCALIntg[j]->SetStats(0);
     hAllSpaceDCALIntg[j]->SetContour(30);
     hAllDistribIntg[j] = new TH1F(Form("hAllDistribIntg_%d",j),Form("hAllDistribIntg_%d",j),100,tabMin[j+kNbFitParams],tabMax[j+kNbFitParams]);
     hAllDistribIntg[j]->SetXTitle(Form("%s",varName[j+kNbFitParams]));
     hAllDistribIntg[j]->SetYTitle("Counts");
     hAllDistribIntg[j]->SetStats(0);
     }
 for (j=0;j<1;j++)
    {hAllSpaceEMCALCoeff[j] = new TH2F(Form("hAllSpaceEMCALCoeff%d",j),Form("hAllSpaceEMCALCoeff%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALCoeff[j]->SetXTitle("Column");
     hAllSpaceEMCALCoeff[j]->SetYTitle("Row");
     hAllSpaceEMCALCoeff[j]->SetStats(0);
     hAllSpaceEMCALCoeff[j]->SetContour(30);
     hAllSpaceDCALCoeff[j] = new TH2F(Form("hAllSpaceDCALCoeff%d",j),Form("hAllSpaceDCALCoeff%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALCoeff[j]->SetXTitle("Column");
     hAllSpaceDCALCoeff[j]->SetYTitle("Row");
     hAllSpaceDCALCoeff[j]->SetStats(0);
     hAllSpaceDCALCoeff[j]->SetContour(30);
     }
 for (j=0;j<3;j++)
    {hAllSpaceEMCALDistance[j] = new TH2F(Form("hAllSpaceEMCALDistance%d",j),Form("hAllSpaceEMCALDistance%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALDistance[j]->SetXTitle("Column");
     hAllSpaceEMCALDistance[j]->SetYTitle("Row");
     hAllSpaceEMCALDistance[j]->SetStats(0);
     hAllSpaceEMCALDistance[j]->SetContour(30);
     hAllSpaceDCALDistance[j] = new TH2F(Form("hAllSpaceDCALDistance%d",j),Form("hAllSpaceDCALDistance%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALDistance[j]->SetXTitle("Column");
     hAllSpaceDCALDistance[j]->SetYTitle("Row");
     hAllSpaceDCALDistance[j]->SetStats(0);
     hAllSpaceDCALDistance[j]->SetContour(30);
     }
 for (j=0;j<kNbFitParams+kNbExtraParams;j++)
    {hAllOldDistrib[j] = new TH1F(Form("hAllOldDistrib_%d",j),Form("hAllOldDistrib_%d",j),100,tabMin[j],tabMax[j]);
     hAllOldDistrib[j]->SetXTitle(Form("%s",varName[j]));
     hAllOldDistrib[j]->SetYTitle("Counts");
     hAllOldDistrib[j]->SetStats(0);
     }
 for (j=0;j<kNbTotParams;j++)
    {hAllSpaceEMCALDiff[j] = new TH2F(Form("hAllSpaceEMCALDiff%d",j),Form("hAllSpaceEMCALDiff%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALDiff[j]->SetXTitle("Column");
     hAllSpaceEMCALDiff[j]->SetYTitle("Row");
     hAllSpaceEMCALDiff[j]->SetStats(0);
     hAllSpaceEMCALDiff[j]->SetContour(30);
     hAllSpaceDCALDiff[j] = new TH2F(Form("hAllSpaceDCALDiff%d",j),Form("hAllSpaceDCALDiff%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALDiff[j]->SetXTitle("Column");
     hAllSpaceDCALDiff[j]->SetYTitle("Row");
     hAllSpaceDCALDiff[j]->SetStats(0);
     hAllSpaceDCALDiff[j]->SetContour(30);
     hAllDiff[j] = new TH1F(Form("hAllDiff%d",j),Form("hAllDiff%d",j),100,tabMinDiff[j],tabMaxDiff[j]);
     hAllDiff[j]->SetXTitle(Form("Ratio var_%d",j));
     hAllDiff[j]->SetYTitle("Counts");
     hAllDiff[j]->SetStats(0);
     hAllDiffAllTw[j] = new TH1F(Form("hAllDiffAllTw%d",j),Form("hAllDiffAllTw%d",j),100,tabMinDiff[j],tabMaxDiff[j]);
     hAllDiffAllTw[j]->SetXTitle(Form("Ratio var_%d",j));
     hAllDiffAllTw[j]->SetYTitle("Counts");
     hAllDiffAllTw[j]->SetStats(0);
     }
  
  TH2F *hCorrelMuVsA = new TH2F("hCorrelMuVsA","hCorrelMuVsA",100,tabMin[0],tabMax[0],100,tabMin[1],tabMax[1]);
  hCorrelMuVsA->SetXTitle("A");
  hCorrelMuVsA->SetYTitle("Mu");
  hCorrelMuVsA->SetStats(0);
  hCorrelMuVsA->SetContour(30);
  TH2F *hCorrelSigVsA = new TH2F("hCorrelSigVsA","hCorrelSigVsA",100,tabMin[0],tabMax[0],100,tabMin[2],tabMax[2]);
  hCorrelSigVsA->SetXTitle("A");
  hCorrelSigVsA->SetYTitle("Sig");
  hCorrelSigVsA->SetStats(0);
  hCorrelSigVsA->SetContour(30);
  TH2F *hCorrelSigVsMu = new TH2F("hCorrelSigVsMu","hCorrelSigVsMu",100,tabMin[1],tabMax[1],100,tabMin[2],tabMax[2]);
  hCorrelSigVsMu->SetXTitle("Mu");
  hCorrelSigVsMu->SetYTitle("Sig");
  hCorrelSigVsMu->SetStats(0);
  hCorrelSigVsMu->SetContour(30);
  TH2F *hCorrelBVsA = new TH2F("hCorrelBVsA","hCorrelBVsA",100,tabMin[3],tabMax[3],100,tabMin[4],tabMax[4]);
  hCorrelBVsA->SetXTitle("C");
  hCorrelBVsA->SetYTitle("B");
  hCorrelBVsA->SetStats(0);
  hCorrelBVsA->SetContour(30);
  TH2F *hCorrelCVsA = new TH2F("hCorrelCVsA","hCorrelCVsA",100,tabMin[3],tabMax[3],100,tabMin[5],tabMax[5]);
  hCorrelCVsA->SetXTitle("C");
  hCorrelCVsA->SetYTitle("A");
  hCorrelCVsA->SetStats(0);
  hCorrelCVsA->SetContour(30);
  TH2F *hCorrelCVsB = new TH2F("hCorrelCVsB","hCorrelCVsB",100,tabMin[4],tabMax[4],100,tabMin[5],tabMax[5]);
  hCorrelCVsB->SetXTitle("B");
  hCorrelCVsB->SetYTitle("A");
  hCorrelCVsB->SetStats(0);
  hCorrelCVsB->SetContour(30);
  TH2F *hCorrelISVsI = new TH2F("hCorrelISVsI","hCorrelISVsI",100,tabMin[6],tabMax[6],100,tabMin[9],tabMax[9]);
  hCorrelISVsI->SetXTitle("I");
  hCorrelISVsI->SetYTitle("IS");
  hCorrelISVsI->SetStats(0);
  hCorrelISVsI->SetContour(30);
  TH2F *hCorrelSVsI = new TH2F("hCorrelSVsI","hCorrelSVsI",100,tabMin[6],tabMax[6],100,tabMin[8],tabMax[8]);
  hCorrelSVsI->SetXTitle("I");
  hCorrelSVsI->SetYTitle("S");
  hCorrelSVsI->SetStats(0);
  hCorrelSVsI->SetContour(30);
  TH2F *hCorrelSVsIS = new TH2F("hCorrelSVsIS","hCorrelSVsIS",100,tabMin[9],tabMax[9],100,tabMin[8],tabMax[8]);
  hCorrelSVsIS->SetXTitle("IS");
  hCorrelSVsIS->SetYTitle("S");
  hCorrelSVsIS->SetStats(0);
  hCorrelSVsIS->SetContour(30);
  TH2F *hCorrelSbincountingVsS = new TH2F("hCorrelSbincountingVsS","hCorrelSbincountingVsS",100,tabMin[8],tabMax[8],100,tabMin[8],tabMax[8]);
  hCorrelSbincountingVsS->SetXTitle("S");
  hCorrelSbincountingVsS->SetYTitle("S bin counting");
  hCorrelSbincountingVsS->SetStats(0);
  hCorrelSbincountingVsS->SetContour(30);
  TH2F *hCorrelMuVsI = new TH2F("hCorrelMuVsI","hCorrelMuVsI",100,tabMin[6],tabMax[6],100,tabMin[1],tabMax[1]);
  hCorrelMuVsI->SetXTitle("I");
  hCorrelMuVsI->SetYTitle("Mu");
  hCorrelMuVsI->SetStats(0);
  hCorrelMuVsI->SetContour(30);
  TH2F *hCorrelMuVsS = new TH2F("hCorrelMuVsS","hCorrelMuVsS",100,tabMin[8],tabMax[8],100,tabMin[1],tabMax[1]);
  hCorrelMuVsS->SetXTitle("S");
  hCorrelMuVsS->SetYTitle("Mu");
  hCorrelMuVsS->SetStats(0);
  hCorrelMuVsS->SetContour(30);
  TH2F *hCorrelSigVsI = new TH2F("hCorrelSigVsI","hCorrelSigVsI",100,tabMin[6],tabMax[6],100,tabMin[2],tabMax[2]);
  hCorrelSigVsI->SetXTitle("I");
  hCorrelSigVsI->SetYTitle("Sig");
  hCorrelSigVsI->SetStats(0);
  hCorrelSigVsI->SetContour(30);
  TH2F *hCorrelSigVsS = new TH2F("hCorrelSigVsS","hCorrelSigVsS",100,tabMin[8],tabMax[8],100,tabMin[2],tabMax[2]);
  hCorrelSigVsS->SetXTitle("S");
  hCorrelSigVsS->SetYTitle("Sig");
  hCorrelSigVsS->SetStats(0);
  hCorrelSigVsS->SetContour(30);
  TH2F *hCorrelDvsSig = new TH2F("hCorrelDvsSig","hCorrelDvsSig",100,tabMin[2],tabMax[2],100,0.,5.);
  hCorrelDvsSig->SetXTitle("#sigma");
  hCorrelDvsSig->SetYTitle("d_{fit}");
  hCorrelDvsSig->SetStats(0);
  hCorrelDvsSig->SetContour(30);
  TH2F *hCorrelCmpD = new TH2F("hCorrelCmpD","hCorrelCmpD",100,0.,5.,100,-1.,5.);
  hCorrelCmpD->SetXTitle("d_{fit}");
  hCorrelCmpD->SetYTitle("(d_{N}-d_{fit})/d_{fit}");
  hCorrelCmpD->SetStats(0);
  hCorrelCmpD->SetContour(30);
  
  TH2F *hAllSpaceEMCALSoverN = new TH2F("hAllSpaceEMCALSoverN","hAllSpaceEMCALSoverN",2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
  hAllSpaceEMCALSoverN->SetXTitle("Column");
  hAllSpaceEMCALSoverN->SetYTitle("Row");
  hAllSpaceEMCALSoverN->SetStats(0);
  hAllSpaceEMCALSoverN->SetContour(30);
  TH2F *hAllSpaceDCALSoverN = new TH2F("hAllSpaceDCALSoverN","hAllSpaceDCALSoverN",2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
  hAllSpaceDCALSoverN->SetXTitle("Column");
  hAllSpaceDCALSoverN->SetYTitle("Row");
  hAllSpaceDCALSoverN->SetStats(0);
  hAllSpaceDCALSoverN->SetContour(30);
  TH1F *hCoeff = new TH1F("hCoeff","hCoeff",100,0.9,1.1);
  hCoeff->SetXTitle("Coeff");
  hCoeff->SetYTitle("Counts");
  hCoeff->SetStats(0);
  TH1F *hCoeffOld = new TH1F("hCoeffOld","hCoeffOld",100,0.9,1.1);
  hCoeffOld->SetXTitle("Coeff");
  hCoeffOld->SetYTitle("Counts");
  hCoeffOld->SetStats(0);
  TH1F *hCoeffOldCorr = new TH1F("hCoeffOldCorr","hCoeffOldCorr",100,0.9,1.1);
  hCoeffOldCorr->SetXTitle("Coeff");
  hCoeffOldCorr->SetYTitle("Counts");
  hCoeffOldCorr->SetStats(0);
  TH1F *hCoeffOldLarge = new TH1F("hCoeffOldLarge","hCoeffOldLarge",100,0.4,1.8);
  hCoeffOldLarge->SetXTitle("Coeff");
  hCoeffOldLarge->SetYTitle("Counts");
  hCoeffOldLarge->SetStats(0);
  TH1F *hCoeffOldCorrLarge = new TH1F("hCoeffOldCorrLarge","hCoeffOldCorrLarge",100,0.4,1.8);
  hCoeffOldCorrLarge->SetXTitle("Coeff");
  hCoeffOldCorrLarge->SetYTitle("Counts");
  hCoeffOldCorrLarge->SetStats(0);
  TH1F *hCtrlEvolMu = new TH1F("hCtrlEvolMu","hCtrlEvolMu",100,-3.,3.);
  hCtrlEvolMu->SetXTitle("(#mu_{2}-m_{PDG})/(#mu_{1}-m_{PDG})");
  hCtrlEvolMu->SetYTitle("Counts");
  hCtrlEvolMu->SetStats(0);
  TH2F *h2CtrlEvolMu1 = new TH2F("h2CtrlEvolMu1","h2CtrlEvolMu1",200,mPDG-10.,mPDG+20.,100,-3.,3.);
  h2CtrlEvolMu1->SetXTitle("#mu_{1}");
  h2CtrlEvolMu1->SetYTitle("(#mu_{2}-m_{PDG})/(#mu_{1}-m_{PDG})");
  h2CtrlEvolMu1->SetStats(0);
  h2CtrlEvolMu1->SetContour(30);
  TH2F *h2CtrlEvolMu2 = new TH2F("h2CtrlEvolMu2","h2CtrlEvolMu2",200,-0.1,0.2,200,-0.1,0.2);
  h2CtrlEvolMu2->SetXTitle("(#mu_{1}-m_{PDG})/m_{PDG}");
  h2CtrlEvolMu2->SetYTitle("(#mu_{1}-#mu_{2})/m_{PDG}");
  h2CtrlEvolMu2->SetStats(0);
  h2CtrlEvolMu2->SetContour(30);
  TH1F *hCtrlEvolCoeff = new TH1F("hCtrlEvolCoeff","hCtrlEvolCoeff",100,-3.,3.);
  hCtrlEvolCoeff->SetXTitle("(c_{2}-1)/(c_{1}-1)");
  hCtrlEvolCoeff->SetYTitle("Counts");
  hCtrlEvolCoeff->SetStats(0);
  TH2F *h2CtrlEvolCoeff1 = new TH2F("h2CtrlEvolCoeff1","h2CtrlEvolCoeff1",200,0.8,1.2,100,-3.,3.);
  h2CtrlEvolCoeff1->SetXTitle("c_{1}");
  h2CtrlEvolCoeff1->SetYTitle("(c_{2}-1)/(c_{1}-1)");
  h2CtrlEvolCoeff1->SetStats(0);
  h2CtrlEvolCoeff1->SetContour(30);
  TH2F *h2CtrlEvolCoeff2 = new TH2F("h2CtrlEvolCoeff2","h2CtrlEvolCoeff2",200,-0.1,0.2,200,-0.1,0.2);
  h2CtrlEvolCoeff2->SetXTitle("c_{1}-1");
  h2CtrlEvolCoeff2->SetYTitle("c_{1}-c_{2}");
  h2CtrlEvolCoeff2->SetStats(0);
  h2CtrlEvolCoeff2->SetContour(30);






 
 
 //*****************************************
 // Fill in pT histoes and graphs :
 //*****************************************


 //hAllMPerSM[(3*(nbPtBins+1))*iSM+3*iPt+iStruct]

 for (i=0;i<nbPtBins-2;i++)
    {tabGraphX[i+1]=tabPtBins[i]+(tabPtBins[i+1]-tabPtBins[i])/2.;
     tabErrGraphX[i+1]=(tabPtBins[i+1]-tabPtBins[i])/2.;
     }
 /*TH1F **hPtInvMass;
 hPtInvMass = new TH1F*[nbPtBins];
 for (j=0;j<nbPtBins;j++)
   {hPtInvMass[j] = new TH1F(Form("hPtInvMass_%d",j),Form("hPtInvMass_%d",j),hAllM_05->GetNbinsX(),hAllM_05->GetXaxis()->GetXmin(),hAllM_05->GetXaxis()->GetXmax());
     hPtInvMass[j]->SetXTitle("#pi^{0} inv mass (MeV)");
     hPtInvMass[j]->SetYTitle("Counts");
     hPtInvMass[j]->SetStats(0);
     hPtInvMass[j]->SetEntries(0);
     }*/
 pTmoyLow=0.;
 pTmoyHigh=0.;
 cmptPtMoyLow=0;
 cmptPtMoyHigh=0;
 for (i=0;i<hAllM_05->GetNbinsY();i++)
    {for (k=0;k<nbPtBins-1;k++)
        {if (hAllM_05->GetYaxis()->GetBinCenter(i+1) < tabPtBins[k]) break;
         }
     for (j=0;j<hAllM_05->GetNbinsX();j++)
        {/*hPtInvMass[k]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05->GetBinContent(j+1,i+1));
         hPtInvMass[k]->SetEntries(hPtInvMass[k]->GetEntries()+hAllM_05->GetBinContent(j+1,i+1));*/
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+0]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+0]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+0]->GetEntries()+hAllM_05->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+0]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+0]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+0]->GetEntries()+hAllM_05->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+1]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_masked->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+1]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+1]->GetEntries()+hAllM_05_masked->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+1]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_masked->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+1]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+1]->GetEntries()+hAllM_05_masked->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+2]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05->GetBinContent(j+1,i+1)-hAllM_05_masked->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+2]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*k+2]->GetEntries()+hAllM_05->GetBinContent(j+1,i+1)-hAllM_05_masked->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+2]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05->GetBinContent(j+1,i+1)-hAllM_05_masked->GetBinContent(j+1,i+1));
         hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+2]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*kNbSMtot+3*nbPtBins+2]->GetEntries()+hAllM_05->GetBinContent(j+1,i+1)-hAllM_05_masked->GetBinContent(j+1,i+1));
         for (m=0;m<kNbZones;m++)
            {hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*kNbSMtot+kNbZones*k+m]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_Zones[m]->GetBinContent(j+1,i+1));
             hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*kNbSMtot+kNbZones*k+m]->SetEntries(hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*kNbSMtot+kNbZones*k+0]->GetEntries()+hAllM_05_Zones[m]->GetBinContent(j+1,i+1));
             hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*kNbSMtot+kNbZones*nbPtBins+m]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_Zones[m]->GetBinContent(j+1,i+1));
             hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*kNbSMtot+kNbZones*nbPtBins+m]->SetEntries(hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*kNbSMtot+kNbZones*nbPtBins+0]->GetEntries()+hAllM_05_Zones[m]->GetBinContent(j+1,i+1));
             }
         for (iSM=0;iSM<kNbSMtot;iSM++)
            {hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+0]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_SM[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+0]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+0]->GetEntries()+hAllM_05_SM[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+0]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_SM[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+0]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+0]->GetEntries()+hAllM_05_SM[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+1]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_SM_masked[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+1]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+1]->GetEntries()+hAllM_05_SM_masked[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+1]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_SM_masked[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+1]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+1]->GetEntries()+hAllM_05_SM_masked[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+2]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_SM[iSM]->GetBinContent(j+1,i+1)-hAllM_05_SM_masked[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+2]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*iSM+3*k+2]->GetEntries()+hAllM_05_SM[iSM]->GetBinContent(j+1,i+1)-hAllM_05_SM_masked[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+2]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_SM[iSM]->GetBinContent(j+1,i+1)-hAllM_05_SM_masked[iSM]->GetBinContent(j+1,i+1));
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+2]->SetEntries(hAllMPerSM[(3*(nbPtBins+1))*iSM+3*nbPtBins+2]->GetEntries()+hAllM_05_SM[iSM]->GetBinContent(j+1,i+1)-hAllM_05_SM_masked[iSM]->GetBinContent(j+1,i+1));
             for (m=0;m<kNbZones;m++)
                {hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*k+m]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_SM_Zones[iSM][m]->GetBinContent(j+1,i+1));
                 hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*k+m]->SetEntries(hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*k+0]->GetEntries()+hAllM_05_SM_Zones[iSM][m]->GetBinContent(j+1,i+1));
                 hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*nbPtBins+m]->Fill(hAllM_05->GetXaxis()->GetBinCenter(j+1),hAllM_05_SM_Zones[iSM][m]->GetBinContent(j+1,i+1));
                 hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*nbPtBins+m]->SetEntries(hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*nbPtBins+0]->GetEntries()+hAllM_05_SM_Zones[iSM][m]->GetBinContent(j+1,i+1));
                 }
             }
         if ((hAllM_05->GetXaxis()->GetBinCenter(j+1) > mPDG-10.) && (hAllM_05->GetXaxis()->GetBinCenter(j+1) < mPDG+10.))
            {if (k == 0)
                {cmptPtMoyLow+=hAllM_05->GetBinContent(j+1,i+1);
                 pTmoyLow+=hAllM_05->GetYaxis()->GetBinCenter(i+1)*hAllM_05->GetBinContent(j+1,i+1);
                 }
             if (k == nbPtBins-1)
                {cmptPtMoyHigh+=hAllM_05->GetBinContent(j+1,i+1);
                 pTmoyHigh+=hAllM_05->GetYaxis()->GetBinCenter(i+1)*hAllM_05->GetBinContent(j+1,i+1);
                 }
             }
         }
     }
 if (cmptPtMoyLow != 0.) pTmoyLow=pTmoyLow/cmptPtMoyLow;
 pTmoyHigh=pTmoyHigh/cmptPtMoyHigh;
 printf("\nMean pT low = %f, high = %f\n\n",pTmoyLow,pTmoyHigh);
 tabGraphX[0]=pTmoyLow;
 tabErrGraphX[0]=tabGraphX[1]-tabErrGraphX[1]-pTmoyLow;
 tabGraphX[nbPtBins-1]=pTmoyHigh;
 tabErrGraphX[nbPtBins-1]=pTmoyHigh-tabGraphX[nbPtBins-2]-tabErrGraphX[nbPtBins-2];
  

 

 





 
 
 //*****************************************
 // Calculate some variables to help later fits convergence :
 //*****************************************


  TH1D *hAllM_05Proj = new TH1D();
  hAllM_05Proj = hAllM_05->ProjectionX("hAllM_05Proj",0,-1,"e");
  hAllM_05Proj->SetLineColor(colorIndex[1]);


  //Fitting function
  Int_t polN = 2;
  
  TF1 *fitfun = 0;
  TF1 *fitfun3 = 0;
  TF1 *fitfunBkgPol2 = 0;
  if      (polN == 1)
    fitfun = new TF1("fitfun",pi0massP1,100,250,5);
  else if (polN == 2)
    fitfun = new TF1("fitfun",pi0massP2,100,250,6);
  fitfun->SetLineColor(kRed);
  fitfun->SetLineWidth(2);
  fitfun3 = new TF1("fitfun3",pi0massP3,100,250,7);
  fitfun3->SetLineColor(kRed);
  fitfun3->SetLineWidth(2);
  
  fitfun->SetParName(0,"A");
  fitfun->SetParName(1,"m_{0}");
  fitfun->SetParName(2,"#sigma");
  fitfun->SetParName(3,"a_{0}");
  fitfun->SetParName(4,"a_{1}");
  if (polN == 2)
    fitfun->SetParName(5,"a_{2}");
  
  maxGaussAmpl=3.e+6;
  fitfun->SetParLimits(0,  8.,maxGaussAmpl);
  fitfun->SetParLimits(1,  80.,160.);
  fitfun->SetParLimits(2,  2.,35.);
  fitfun3->SetParLimits(0,  8.,maxGaussAmpl);
  fitfun3->SetParLimits(1,  80.,160.);
  fitfun3->SetParLimits(2,  2.,35.);

  fitfunBkgPol2 = new TF1("fitfunBkgPol2","[0] + [1]*x + [2]*x*x",100.,250.);
  fitfunBkgPol2->SetLineWidth(1);
  fitfunBkgPol2->SetLineColor(kMagenta);
 
  gStyle->SetOptTitle(0);
  hAllM_05Proj->Rebin(rebin);
  
  hAllM_05Proj->SetXTitle("M_{#gamma,#gamma} (MeV/c^{2})");
  
  tmpX=210.;
  tmpY=hAllM_05Proj->GetBinContent(hAllM_05Proj->FindBin(tmpX));
  fitPar4=2.*tmpY/tmpX;
  fitPar5=-tmpY/(tmpX*tmpX);
  fitfun->SetParameter(0,(double)hAllM_05Proj->GetEntries()/100.);
  fitfun->SetParameter(1,mPDG);
  fitfun->SetParameter(2,typicalWidth);
  fitfun->SetParameter(3,0.);
  fitfun->SetParameter(4,fitPar4);
  fitfun->SetParameter(5,fitPar5);
  printf("Fit start parameters : amplitude %f, b %f, a %f\n\n",fitfun->GetParameter(0),fitfun->GetParameter(4),fitfun->GetParameter(5));
  hAllM_05Proj->Fit("fitfun","","",50,250);
  fitPar3=fitfun->GetParameter(3);
  fitPar4=fitfun->GetParameter(4);
  fitPar5=fitfun->GetParameter(5);
  nbTot=hAllM_05Proj->GetEntries();
  nbSig=fitfun->GetParameter(0);


  
   





 
 
 //*****************************************
 // Create root and postscript output files :
 //*****************************************

 TFile *rootFileOut = new TFile("output_calibPi0.root","RECREATE");
 char psfile[100];
 sprintf(psfile,"output_calibPi0.ps");
 const int cWidth=500;
 const int cHeight=(int)(500*(29./21.));
 TCanvas *c1 = new TCanvas("c1","EMCal cosmics analysis",cWidth,cHeight);
 TPostScript *ps = new TPostScript(psfile,111);


 hAllM_05Proj->Write();
 




 
 
 //*****************************************
 // Draw pT histoes and graphs :
 //*****************************************


 //hAllMPerSM[(3*(nbPtBins+1))*iSM+3*iPt+iStruct]
 
 // Inv mass plots per pT bin per SM :
 printf("\nInv mass plots per pT bin per SM.\n");
 linePDGMass->SetLineColor(kBlue);
 for (iSM=0;iSM<(kNbSMtot+1);iSM++)
    {if ((iSM<kNbSMtot) && (tabChoiceCalos[SMdetType[iSM]] == 0))
        {printf("Skip SM %d of type %s.\n",iSM,detTypeString[SMdetType[iSM]]);
         continue;
         }
     for (j=0;j<(int)((double)nbPtBins/8.)+1;j++)
        {ps->NewPage();
         c1->Clear();
         c1->Divide(2,4);
         for (i=0;i<8;i++)
            {if ((j*8+i) > nbPtBins) break;
             c1->cd(i+1);
             if ((j*8+i) != nbPtBins) printf("\nNow fitting pT bin %d",j*8+i);
               else printf("\nNow fitting all pT");
             if (iSM != kNbSMtot) printf(" in SM %d -- nEntries = %f\n",iSM,hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+0]->GetEntries());
                else printf(" for all SMs -- nEntries = %f\n",hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+0]->GetEntries());
             hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+0]->Draw();
             if ((j*8+i) < (nbPtBins-1))
                {fitRangeMin=14.28*tabPtBins[j*8+i]+37.15;
                 if (tabPtBins[j*8+i] > 3.0) fitRangeMin=80.;
                 fitRangeMax=11.54*tabPtBins[j*8+i]+201.5;
                 }
                 else //For the 2 last histoes (last pT bin and pT-integrated)
                {fitRangeMin=80.;
                 fitRangeMax=11.54*tabPtBins[nbPtBins-2]+201.5;
                 }
             for (iStruct=0;iStruct<3;iStruct++)
                {fitfun->SetParameter(0,hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+iStruct]->GetEntries()*(double)nbSig/(double)nbTot);
                 fitfun->SetParameter(1,mPDG);
                 fitfun->SetParameter(2,typicalWidth);
                 fitfun->SetParameter(3,fitPar3*(double)hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+iStruct]->GetEntries()/(double)nbTot);
                 fitfun->SetParameter(4,fitPar4*(double)hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+iStruct]->GetEntries()/(double)nbTot);
                 fitfun->SetParameter(5,fitPar5*(double)hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+iStruct]->GetEntries()/(double)nbTot);
                 //One has to put the histo's drawing options as 3rd argument ; drawing the histoes afterwards doesn't work ! (Only the last one is drawn, because fitting draws it without option SAME).
                 hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+iStruct]->Fit("fitfun","R","SAME",fitRangeMin,fitRangeMax);
                 switch (iStruct)
                    {case 0 : lineFittedMass->SetLineColor(colorIndexStruct[iStruct]);
                              lineFittedMass->SetX1(fitfun->GetParameter(1));
                              lineFittedMass->SetX2(fitfun->GetParameter(1));
                              lineFittedMass->SetY2(1.05*hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+0]->GetMaximum());
                              lineFittedMass->Draw("SAME");
                              break;
                     case 1 : lineFittedWidth1->SetLineColor(colorIndexStruct[iStruct]);
                              lineFittedWidth1->SetX1(fitfun->GetParameter(1));
                              lineFittedWidth1->SetX2(fitfun->GetParameter(1));
                              lineFittedWidth1->SetY2(1.05*hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+0]->GetMaximum());
                              lineFittedWidth1->Draw("SAME");
                              break;
                     case 2 : lineFittedWidth2->SetLineColor(colorIndexStruct[iStruct]);
                              lineFittedWidth2->SetX1(fitfun->GetParameter(1));
                              lineFittedWidth2->SetX2(fitfun->GetParameter(1));
                              lineFittedWidth2->SetY2(1.05*hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+0]->GetMaximum());
                              lineFittedWidth2->Draw("SAME");
                              break;
                     default : ;
                     }
                 if ((j*8+i) < nbPtBins)
                    {tabGraphYmass[3*(j*8+i)+iStruct]=fitfun->GetParameter(1);
                     tabGraphYsig[3*(j*8+i)+iStruct]=fitfun->GetParameter(2);
                     tabGraphYnumber[3*(j*8+i)+iStruct]=fitfun->GetParameter(2)*fitfun->GetParameter(0)*TMath::Sqrt(2.*TMath::Pi())/hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+iStruct]->GetBinWidth(1);
                     tabErrGraphYmass[3*(j*8+i)+iStruct]=fitfun->GetParError(1);
                     tabErrGraphYsig[3*(j*8+i)+iStruct]=fitfun->GetParError(2);
                     if (fitfun->Integral(mPDG-3.*fitfun->GetParameter(2),mPDG+3.*fitfun->GetParameter(2)) >= 0) tabErrGraphYnumber[3*(j*8+i)+iStruct]=TMath::Sqrt(fitfun->Integral(mPDG-3.*fitfun->GetParameter(2),mPDG+3.*fitfun->GetParameter(2)));
                        else tabErrGraphYnumber[3*(j*8+i)+iStruct]=0.; //Fit integral can be <0 (due to fitted bkgnd).
                     }
                 hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+iStruct]->Write();
                 }
             if ((j*8+i) < nbPtBins)
                {if ((tabGraphYnumber[3*(j*8+i)+1] != 0.) && (tabGraphYnumber[3*(j*8+i)+2] != 0.))
                    {tabGraphYnumberRatio[3*(j*8+i)+0]=tabGraphYnumber[3*(j*8+i)+2]/tabGraphYnumber[3*(j*8+i)+1];
                     tabErrGraphYnumberRatio[3*(j*8+i)+0]=tabGraphYnumberRatio[3*(j*8+i)+0]*TMath::Sqrt((tabErrGraphYnumber[3*(j*8+i)+2]/tabGraphYnumber[3*(j*8+i)+2])*(tabErrGraphYnumber[3*(j*8+i)+2]/tabGraphYnumber[3*(j*8+i)+2])+(tabErrGraphYnumber[3*(j*8+i)+1]/tabGraphYnumber[3*(j*8+i)+1])*(tabErrGraphYnumber[3*(j*8+i)+1]/tabGraphYnumber[3*(j*8+i)+1]));
                     }
                    else
                    {tabGraphYnumberRatio[3*(j*8+i)+0]=0.;
                     tabErrGraphYnumberRatio[3*(j*8+i)+0]=0.;
                     }
                 }
             linePDGMass->SetY2(1.05*hAllMPerSM[(3*(nbPtBins+1))*iSM+3*(j*8+i)+0]->GetMaximum());
             linePDGMass->Draw("SAME");
             c1->Update(); //Has to come before the ps->TextNDC, otherwise the text is drawn below the histo (=> invisible).
             if (iSM != kNbSMtot) sprintf(txtPS,"Pi0 inv mass in SM %d,",iSM);
                else sprintf(txtPS,"Pi0 inv mass in all SMs,");
             if ((j*8+i) != nbPtBins) sprintf(txtPS,"%s %4.2f < pT < %4.2f GeV",txtPS,tabGraphX[j*8+i]-tabErrGraphX[j*8+i],tabGraphX[j*8+i]+tabErrGraphX[j*8+i]);
                else sprintf(txtPS,"%s pT-integrated",txtPS);
             ps->TextNDC(0.3,0.93,txtPS);
             }
         c1->Update();
         }
     for (iStruct=0;iStruct<3;iStruct++)
        {for (iPt=0;iPt<nbPtBins;iPt++)
            {tabTmpGraphYmass[iPt]=tabGraphYmass[3*iPt+iStruct];
             tabTmpErrGraphYmass[iPt]=tabErrGraphYmass[3*iPt+iStruct];
             tabTmpGraphYsig[iPt]=tabGraphYsig[3*iPt+iStruct];
             tabTmpErrGraphYsig[iPt]=tabErrGraphYsig[3*iPt+iStruct];
             tabTmpGraphYnumber[iPt]=tabGraphYnumber[3*iPt+iStruct];
             tabTmpErrGraphYnumber[iPt]=tabErrGraphYnumber[3*iPt+iStruct];
             tabTmpGraphYnumberRatio[iPt]=tabGraphYnumberRatio[3*iPt+iStruct];
             tabTmpErrGraphYnumberRatio[iPt]=tabErrGraphYnumberRatio[3*iPt+iStruct];
             }
         tgeMass[3*iSM+iStruct] = new TGraphErrors(nbPtBins,tabGraphX,tabTmpGraphYmass,tabErrGraphX,tabTmpErrGraphYmass);
         tgeSig[3*iSM+iStruct] = new TGraphErrors(nbPtBins,tabGraphX,tabTmpGraphYsig,tabErrGraphX,tabTmpErrGraphYsig);
         tgeNumber[3*iSM+iStruct] = new TGraphErrors(nbPtBins,tabGraphX,tabTmpGraphYnumber,tabErrGraphX,tabTmpErrGraphYnumber);
         tgeNumberRatio[3*iSM+iStruct] = new TGraphErrors(nbPtBins,tabGraphX,tabTmpGraphYnumberRatio,tabErrGraphX,tabTmpErrGraphYnumberRatio);
         tgeMass[3*iSM+iStruct]->SetName(Form("tgeMass_SM%d_%d",iSM,iStruct));
         tgeSig[3*iSM+iStruct]->SetName(Form("tgeSig_SM%d_%d",iSM,iStruct));
         tgeNumber[3*iSM+iStruct]->SetName(Form("tgeNumber_SM%d_%d",iSM,iStruct));
         tgeNumberRatio[3*iSM+iStruct]->SetName(Form("tgeNumberRatio_SM%d_%d",iSM,iStruct));
         tgeMass[3*iSM+iStruct]->SetLineColor(colorIndexStruct[iStruct]);
         tgeSig[3*iSM+iStruct]->SetLineColor(colorIndexStruct[iStruct]);
         tgeNumber[3*iSM+iStruct]->SetLineColor(colorIndexStruct[iStruct]);
         tgeNumberRatio[3*iSM+iStruct]->SetLineColor(colorIndexStruct[iStruct]);
         tgeMass[3*iSM+iStruct]->GetHistogram()->SetXTitle("#pi^{0} p_{#perp} (GeV)");
         tgeMass[3*iSM+iStruct]->GetHistogram()->SetYTitle("Fitted #pi^{0} inv mass (MeV)");
         tgeMass[3*iSM+iStruct]->GetHistogram()->SetMinimum(120.);
         tgeMass[3*iSM+iStruct]->GetHistogram()->SetMaximum(150.);
         tgeSig[3*iSM+iStruct]->GetHistogram()->SetXTitle("#pi^{0} p_{#perp} (GeV)");
         tgeSig[3*iSM+iStruct]->GetHistogram()->SetYTitle("Fitted #pi^{0} width (MeV)");
         tgeSig[3*iSM+iStruct]->GetHistogram()->SetMinimum(7.5);
         tgeSig[3*iSM+iStruct]->GetHistogram()->SetMaximum(22.);
         tgeNumber[3*iSM+iStruct]->GetHistogram()->SetXTitle("#pi^{0} p_{#perp} (GeV)");
         tgeNumber[3*iSM+iStruct]->GetHistogram()->SetYTitle("Number of #pi^{0} (counts)");
         tgeNumberRatio[3*iSM+iStruct]->GetHistogram()->SetXTitle("#pi^{0} p_{#perp} (GeV)");
         tgeNumberRatio[3*iSM+iStruct]->GetHistogram()->SetYTitle("Ratios of numbers of #pi^{0}");
         }
     }

 for (iSM=0;iSM<(kNbSMtot+1);iSM++)
    {if ((iSM<kNbSMtot) && (tabChoiceCalos[SMdetType[iSM]] == 0))
        {printf("Skip SM %d of type %s.\n",iSM,detTypeString[SMdetType[iSM]]);
         continue;
         }
     for (j=0;j<(int)((double)nbPtBins/8.)+1;j++)
        {ps->NewPage();
         c1->Clear();
         c1->Divide(2,4);
         for (i=0;i<8;i++)
            {if ((j*8+i) > nbPtBins) break;
             c1->cd(i+1);
             if ((j*8+i) != nbPtBins) printf("\nNow fitting pT bin %d",j*8+i);
               else printf("\nNow fitting all pT");
             if (iSM != kNbSMtot) printf(" in SM %d -- nEntries = %f\n",iSM,hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+0]->GetEntries());
                else printf(" for all SMs -- nEntries = %f\n",hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+0]->GetEntries());
             hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+3]->Draw(); //3rd of these histoes is the highest.
             if ((j*8+i) < (nbPtBins-1))
                {fitRangeMin=14.28*tabPtBins[j*8+i]+37.15;
                 if (tabPtBins[j*8+i] > 3.0) fitRangeMin=80.;
                 fitRangeMax=11.54*tabPtBins[j*8+i]+201.5;
                 }
                 else //For the 2 last histoes (last pT bin and pT-integrated)
                {fitRangeMin=80.;
                 fitRangeMax=11.54*tabPtBins[nbPtBins-2]+201.5;
                 }
             for (iStruct=0;iStruct<kNbZones;iStruct++)
                {fitfun->SetParameter(0,hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+iStruct]->GetEntries()*(double)nbSig/(double)nbTot);
                 fitfun->SetParameter(1,mPDG);
                 fitfun->SetParameter(2,typicalWidth);
                 fitfun->SetParameter(3,fitPar3*(double)hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+iStruct]->GetEntries()/(double)nbTot);
                 fitfun->SetParameter(4,fitPar4*(double)hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+iStruct]->GetEntries()/(double)nbTot);
                 fitfun->SetParameter(5,fitPar5*(double)hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+iStruct]->GetEntries()/(double)nbTot);
                 //One has to put the histo's drawing options as 3rd argument ; drawing the histoes afterwards doesn't work ! (Only the last one is drawn, because fitting draws it without option SAME).
                 hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+iStruct]->Fit("fitfun","R","SAME",fitRangeMin,fitRangeMax);
                 if ((j*8+i) < nbPtBins)
                    {tabGraphYmass[kNbZones*(j*8+i)+iStruct]=fitfun->GetParameter(1);
                     tabGraphYsig[kNbZones*(j*8+i)+iStruct]=fitfun->GetParameter(2);
                     tabGraphYnumber[kNbZones*(j*8+i)+iStruct]=fitfun->GetParameter(2)*fitfun->GetParameter(0)*TMath::Sqrt(2.*TMath::Pi())/hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+iStruct]->GetBinWidth(1);
                     tabErrGraphYmass[kNbZones*(j*8+i)+iStruct]=fitfun->GetParError(1);
                     tabErrGraphYsig[kNbZones*(j*8+i)+iStruct]=fitfun->GetParError(2);
                     if (fitfun->Integral(mPDG-3.*fitfun->GetParameter(2),mPDG+3.*fitfun->GetParameter(2)) >= 0) tabErrGraphYnumber[kNbZones*(j*8+i)+iStruct]=TMath::Sqrt(fitfun->Integral(mPDG-3.*fitfun->GetParameter(2),mPDG+3.*fitfun->GetParameter(2)));
                        else tabErrGraphYnumber[kNbZones*(j*8+i)+iStruct]=0.; //Fit integral can be <0 (due to fitted bkgnd).
                     }
                 hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+iStruct]->Write();
                 }
             if ((j*8+i) < nbPtBins)
                {if ((tabGraphYnumber[kNbZones*(j*8+i)+0] != 0.) && (tabGraphYnumber[kNbZones*(j*8+i)+1] != 0.))
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+0]=tabGraphYnumber[kNbZones*(j*8+i)+1]/tabGraphYnumber[kNbZones*(j*8+i)+0];
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+0]=tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+0]*TMath::Sqrt((tabErrGraphYnumber[kNbZones*(j*8+i)+1]/tabGraphYnumber[kNbZones*(j*8+i)+1])*(tabErrGraphYnumber[kNbZones*(j*8+i)+1]/tabGraphYnumber[kNbZones*(j*8+i)+1])+(tabErrGraphYnumber[kNbZones*(j*8+i)+0]/tabGraphYnumber[kNbZones*(j*8+i)+0])*(tabErrGraphYnumber[kNbZones*(j*8+i)+0]/tabGraphYnumber[kNbZones*(j*8+i)+0]));
                     }
                    else
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+0]=0.;
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+0]=0.;
                     }
                 if ((tabGraphYnumber[kNbZones*(j*8+i)+2] != 0.) && (tabGraphYnumber[kNbZones*(j*8+i)+3] != 0.))
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+1]=tabGraphYnumber[kNbZones*(j*8+i)+3]/tabGraphYnumber[kNbZones*(j*8+i)+2];
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+1]=tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+2]*TMath::Sqrt((tabErrGraphYnumber[kNbZones*(j*8+i)+3]/tabGraphYnumber[kNbZones*(j*8+i)+3])*(tabErrGraphYnumber[kNbZones*(j*8+i)+3]/tabGraphYnumber[kNbZones*(j*8+i)+3])+(tabErrGraphYnumber[kNbZones*(j*8+i)+2]/tabGraphYnumber[kNbZones*(j*8+i)+2])*(tabErrGraphYnumber[kNbZones*(j*8+i)+2]/tabGraphYnumber[kNbZones*(j*8+i)+2]));
                     }
                    else
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+1]=0.;
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+1]=0.;
                     }
                 if ((tabGraphYnumber[kNbZones*(j*8+i)+4] != 0.) && (tabGraphYnumber[kNbZones*(j*8+i)+5] != 0.))
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+2]=tabGraphYnumber[kNbZones*(j*8+i)+5]/tabGraphYnumber[kNbZones*(j*8+i)+4];
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+2]=tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+4]*TMath::Sqrt((tabErrGraphYnumber[kNbZones*(j*8+i)+5]/tabGraphYnumber[kNbZones*(j*8+i)+5])*(tabErrGraphYnumber[kNbZones*(j*8+i)+5]/tabGraphYnumber[kNbZones*(j*8+i)+5])+(tabErrGraphYnumber[kNbZones*(j*8+i)+4]/tabGraphYnumber[kNbZones*(j*8+i)+4])*(tabErrGraphYnumber[kNbZones*(j*8+i)+4]/tabGraphYnumber[kNbZones*(j*8+i)+4]));
                     }
                    else
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+2]=0.;
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+2]=0.;
                     }
                 if ((tabGraphYnumber[kNbZones*(j*8+i)+5] != 0.) && (tabGraphYnumber[kNbZones*(j*8+i)+6] != 0.))
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+3]=tabGraphYnumber[kNbZones*(j*8+i)+6]/tabGraphYnumber[kNbZones*(j*8+i)+5];
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+3]=tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+5]*TMath::Sqrt((tabErrGraphYnumber[kNbZones*(j*8+i)+6]/tabGraphYnumber[kNbZones*(j*8+i)+6])*(tabErrGraphYnumber[kNbZones*(j*8+i)+6]/tabGraphYnumber[kNbZones*(j*8+i)+6])+(tabErrGraphYnumber[kNbZones*(j*8+i)+5]/tabGraphYnumber[kNbZones*(j*8+i)+5])*(tabErrGraphYnumber[kNbZones*(j*8+i)+5]/tabGraphYnumber[kNbZones*(j*8+i)+5]));
                     }
                    else
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+3]=0.;
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+3]=0.;
                     }
                 if ((tabGraphYnumber[kNbZones*(j*8+i)+4] != 0.) && (tabGraphYnumber[kNbZones*(j*8+i)+6] != 0.))
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+4]=tabGraphYnumber[kNbZones*(j*8+i)+6]/tabGraphYnumber[kNbZones*(j*8+i)+4];
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+4]=tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+4]*TMath::Sqrt((tabErrGraphYnumber[kNbZones*(j*8+i)+6]/tabGraphYnumber[kNbZones*(j*8+i)+6])*(tabErrGraphYnumber[kNbZones*(j*8+i)+6]/tabGraphYnumber[kNbZones*(j*8+i)+6])+(tabErrGraphYnumber[kNbZones*(j*8+i)+4]/tabGraphYnumber[kNbZones*(j*8+i)+4])*(tabErrGraphYnumber[kNbZones*(j*8+i)+4]/tabGraphYnumber[kNbZones*(j*8+i)+4]));
                     }
                    else
                    {tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+4]=0.;
                     tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*(j*8+i)+4]=0.;
                     }
                 }
             linePDGMass->SetY2(1.05*hMassDistrPerSMPerPtPerZone[(kNbZones*(nbPtBins+1))*iSM+kNbZones*(j*8+i)+0]->GetMaximum());
             linePDGMass->Draw("SAME");
             c1->Update(); //Has to come before the ps->TextNDC, otherwise the text is drawn below the histo (=> invisible).
             if (iSM != kNbSMtot) sprintf(txtPS,"Pi0 inv mass in SM %d,",iSM);
                else sprintf(txtPS,"Pi0 inv mass in all SMs,");
             if ((j*8+i) != nbPtBins) sprintf(txtPS,"%s %4.2f < pT < %4.2f GeV",txtPS,tabGraphX[j*8+i]-tabErrGraphX[j*8+i],tabGraphX[j*8+i]+tabErrGraphX[j*8+i]);
                else sprintf(txtPS,"%s pT-integrated",txtPS);
             ps->TextNDC(0.3,0.93,txtPS);
             }
         c1->Update();
         }
     for (iStruct=0;iStruct<kNbZones;iStruct++)
        {for (iPt=0;iPt<nbPtBins;iPt++)
            {tabTmpGraphYmass[iPt]=tabGraphYmass[kNbZones*iPt+iStruct];
             tabTmpErrGraphYmass[iPt]=tabErrGraphYmass[kNbZones*iPt+iStruct];
             tabTmpGraphYsig[iPt]=tabGraphYsig[kNbZones*iPt+iStruct];
             tabTmpErrGraphYsig[iPt]=tabErrGraphYsig[kNbZones*iPt+iStruct];
             tabTmpGraphYnumber[iPt]=tabGraphYnumber[kNbZones*iPt+iStruct];
             tabTmpErrGraphYnumber[iPt]=tabErrGraphYnumber[kNbZones*iPt+iStruct];
             }
         tgeMassZones[kNbZones*iSM+iStruct] = new TGraphErrors(nbPtBins,tabGraphX,tabTmpGraphYmass,tabErrGraphX,tabTmpErrGraphYmass);
         tgeSigZones[kNbZones*iSM+iStruct] = new TGraphErrors(nbPtBins,tabGraphX,tabTmpGraphYsig,tabErrGraphX,tabTmpErrGraphYsig);
         tgeNumberZones[kNbZones*iSM+iStruct] = new TGraphErrors(nbPtBins,tabGraphX,tabTmpGraphYnumber,tabErrGraphX,tabTmpErrGraphYnumber);
         tgeMassZones[kNbZones*iSM+iStruct]->SetName(Form("tgeMassZones_SM%d_%d",iSM,iStruct));
         tgeSigZones[kNbZones*iSM+iStruct]->SetName(Form("tgeSigZones_SM%d_%d",iSM,iStruct));
         tgeNumberZones[kNbZones*iSM+iStruct]->SetName(Form("tgeNumberZones_SM%d_%d",iSM,iStruct));
         tgeMassZones[kNbZones*iSM+iStruct]->SetLineColor(colorIndexZones[iStruct]);
         tgeSigZones[kNbZones*iSM+iStruct]->SetLineColor(colorIndexZones[iStruct]);
         tgeNumberZones[kNbZones*iSM+iStruct]->SetLineColor(colorIndexZones[iStruct]);
         tgeMassZones[kNbZones*iSM+iStruct]->GetHistogram()->SetXTitle("#pi^{0} p_{#perp} (GeV)");
         tgeMassZones[kNbZones*iSM+iStruct]->GetHistogram()->SetYTitle("Fitted #pi^{0} inv mass (MeV)");
         tgeMassZones[kNbZones*iSM+iStruct]->GetHistogram()->SetMinimum(120.);
         tgeMassZones[kNbZones*iSM+iStruct]->GetHistogram()->SetMaximum(150.);
         tgeSigZones[kNbZones*iSM+iStruct]->GetHistogram()->SetXTitle("#pi^{0} p_{#perp} (GeV)");
         tgeSigZones[kNbZones*iSM+iStruct]->GetHistogram()->SetYTitle("Fitted #pi^{0} width (MeV)");
         tgeSigZones[kNbZones*iSM+iStruct]->GetHistogram()->SetMinimum(7.5);
         tgeSigZones[kNbZones*iSM+iStruct]->GetHistogram()->SetMaximum(22.);
         tgeNumberZones[kNbZones*iSM+iStruct]->GetHistogram()->SetXTitle("#pi^{0} p_{#perp} (GeV)");
         tgeNumberZones[kNbZones*iSM+iStruct]->GetHistogram()->SetYTitle("Number of #pi^{0} (counts)");
         }
     for (iStruct=0;iStruct<((int)(kNbZones/2))*(kNbZones-1);iStruct++)
        {for (iPt=0;iPt<nbPtBins;iPt++)
            {tabTmpGraphYnumberRatio[iPt]=tabGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*iPt+iStruct];
             tabTmpErrGraphYnumberRatio[iPt]=tabErrGraphYnumberRatio[((int)(kNbZones/2))*(kNbZones-1)*iPt+iStruct];
             }
         tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+iStruct] = new TGraphErrors(nbPtBins,tabGraphX,tabTmpGraphYnumberRatio,tabErrGraphX,tabTmpErrGraphYnumberRatio);
         tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+iStruct]->SetName(Form("tgeNumberRatioZones_SM%d_%d",iSM,iStruct));
         tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+iStruct]->SetLineColor(colorIndexRatiosZones[iStruct]);
         tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+iStruct]->GetHistogram()->SetXTitle("#pi^{0} p_{#perp} (GeV)");
         tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+iStruct]->GetHistogram()->SetYTitle("Ratios of numbers of #pi^{0}");
         }
     }
  
 //pT vs mInv TH2F graphs for all SMs :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,2);
 c1->cd(1);
 gPad->SetLogz();
 hAllM_05->SetContour(30);
 hAllM_05->Draw("COLZ");
 c1->Update();
 sprintf(txtPS,"Pi0 inv mass vs pT, all SMs");
 ps->TextNDC(0.2,0.94,txtPS);
 hAllM_05->Write();
 c1->cd(2);
 gPad->SetLogz();
 hAllM_05_masked->SetContour(30);
 hAllM_05_masked->Draw("COLZ");
 c1->Update();
 sprintf(txtPS,"Pi0 inv mass vs pT, all SMs");
 ps->TextNDC(0.2,0.94,txtPS);
 hAllM_05_masked->Write();
 c1->Update();
 
 // Peak characteristics vs pT, for all SMs :
 printf("\nPeak characteristics vs pT.\n");
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,2);
 c1->cd(1);
 tgeNumber[3*kNbSMtot+0]->SetName("tgeNumber_ALLSM");
 tgeNumber[3*kNbSMtot+1]->SetName("tgeNumber_ALLSM1");
 tgeNumber[3*kNbSMtot+2]->SetName("tgeNumber_ALLSM2");
 tgeNumber[3*kNbSMtot+0]->Draw("AL");
 tgeNumber[3*kNbSMtot+0]->GetHistogram()->SetMinimum(0.);
 tgeNumber[3*kNbSMtot+0]->Draw("AL");
 tgeNumber[3*kNbSMtot+1]->Draw("SAME");
 tgeNumber[3*kNbSMtot+2]->Draw("SAME");
 c1->Update();
 sprintf(txtPS,"Integral of gaussian fit wrt pT, all SMs");
 ps->TextNDC(0.2,0.94,txtPS);
 for (iStruct=0;iStruct<3;iStruct++) tgeNumber[3*kNbSMtot+iStruct]->Write();
 c1->cd(2);
 tgeNumberRatio[3*kNbSMtot+0]->SetName("tgeNumberRatio_ALLSM");
 tgeNumberRatio[3*kNbSMtot+1]->SetName("tgeNumberRatio_ALLSM1");
 tgeNumberRatio[3*kNbSMtot+2]->SetName("tgeNumberRatio_ALLSM2");
 tgeNumberRatio[3*kNbSMtot+0]->Draw("AL");
 tgeNumberRatio[3*kNbSMtot+0]->GetHistogram()->SetMinimum(0.0);
 tgeNumberRatio[3*kNbSMtot+0]->GetHistogram()->SetMaximum(0.6);
 tgeNumberRatio[3*kNbSMtot+0]->Draw("AL");
 for (iStruct=0;iStruct<3;iStruct++) tgeNumberRatio[3*kNbSMtot+iStruct]->Write();
 c1->cd(3);
 TLine *lineMassH = new TLine(tgeMass[3*kNbSMtot+0]->GetXaxis()->GetXmin(),mPDG,tgeMass[3*kNbSMtot+0]->GetXaxis()->GetXmax(),mPDG);
 lineMassH->SetLineColor(kMagenta);
 tgeMass[3*kNbSMtot+0]->Draw("AL");
 tgeMass[3*kNbSMtot+1]->Draw("L SAME");
 tgeMass[3*kNbSMtot+2]->Draw("L SAME");
 lineMassH->Draw();
 c1->Update();
 sprintf(txtPS,"Mean of gaussian fit wrt pT, all SMs");
 ps->TextNDC(0.2,0.94,txtPS);
 for (iStruct=0;iStruct<3;iStruct++) tgeMass[3*kNbSMtot+iStruct]->Write();
 c1->cd(4);
 tgeSig[3*kNbSMtot+0]->Draw("AL");
 tgeSig[3*kNbSMtot+1]->Draw("L SAME");
 tgeSig[3*kNbSMtot+2]->Draw("L SAME");
 c1->Update();
 sprintf(txtPS,"Width of gaussian fit wrt pT, all SMs");
 ps->TextNDC(0.2,0.94,txtPS);
 for (iStruct=0;iStruct<3;iStruct++) tgeSig[3*kNbSMtot+iStruct]->Write();
 c1->Update();
 
 // Peak characteristics vs pT, per SM :
 for (iSM=0;iSM<kNbSMtot;iSM++)
    {if (tabChoiceCalos[SMdetType[iSM]] == 0)
        {printf("Skip SM %d of type %s.\n",iSM,detTypeString[SMdetType[iSM]]);
         continue;
         }
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,2);
     c1->cd(1);
     tgeNumber[3*iSM+0]->Draw("AL");
     tgeNumber[3*iSM+0]->GetHistogram()->SetMinimum(0.);
     tgeNumber[3*iSM+0]->Draw("AL");
     tgeNumber[3*iSM+1]->Draw("L SAME");
     tgeNumber[3*iSM+2]->Draw("L SAME");
     c1->Update();
     sprintf(txtPS,"Integral of gaussian fit wrt pT, SM %d",iSM);
     ps->TextNDC(0.2,0.94,txtPS);
     for (iStruct=0;iStruct<3;iStruct++) tgeNumber[3*iSM+iStruct]->Write();
     c1->cd(2);
     tgeNumberRatio[3*iSM+0]->Draw("AL");
     tgeNumberRatio[3*iSM+0]->GetHistogram()->SetMinimum(0.0);
     tgeNumberRatio[3*iSM+0]->GetHistogram()->SetMaximum(0.6);
     tgeNumberRatio[3*iSM+0]->Draw("AL");
     for (iStruct=0;iStruct<3;iStruct++) tgeNumberRatio[3*iSM+iStruct]->Write();
     c1->cd(3);
     lineMassH->SetX1(tgeMass[3*iSM+0]->GetXaxis()->GetXmin());
     lineMassH->SetX2(tgeMass[3*iSM+0]->GetXaxis()->GetXmax());
     lineMassH->SetLineColor(kMagenta);
     tgeMass[3*iSM+0]->Draw("AL");
     tgeMass[3*iSM+1]->Draw("L SAME");
     tgeMass[3*iSM+2]->Draw("L SAME");
     lineMassH->Draw();
     c1->Update();
     sprintf(txtPS,"Mean of gaussian fit wrt pT, SM %d",iSM);
     ps->TextNDC(0.2,0.94,txtPS);
     for (iStruct=0;iStruct<3;iStruct++) tgeMass[3*iSM+iStruct]->Write();
     c1->cd(4);
     tgeSig[3*iSM+0]->Draw("AL");
     tgeSig[3*iSM+1]->Draw("L SAME");
     tgeSig[3*iSM+2]->Draw("L SAME");
     c1->Update();
     sprintf(txtPS,"Width of gaussian fit wrt pT, SM %d",iSM);
     ps->TextNDC(0.2,0.94,txtPS);
     for (iStruct=0;iStruct<3;iStruct++) tgeSig[3*iSM+iStruct]->Write();
     c1->Update();
     }
 
 // Peak characteristics in zones vs pT, for all SMs :
 printf("\nPeak characteristics vs pT for specific zones.\n");
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,2);
 c1->cd(1);
 for (iStruct=0;iStruct<kNbZones;iStruct++) tgeNumberZones[kNbZones*kNbSMtot+iStruct]->SetName(Form("tgeNumberZones_ALLSM%d",iStruct));
 tgeNumberZones[kNbZones*kNbSMtot+3]->Draw("AL"); //3rd of these histoes is the highest.
 tgeNumberZones[kNbZones*kNbSMtot+3]->GetHistogram()->SetMinimum(0.);
 tgeNumberZones[kNbZones*kNbSMtot+3]->Draw("AL"); //3rd of these histoes is the highest.
 for (iStruct=0;iStruct<kNbZones;iStruct++) tgeNumberZones[kNbZones*kNbSMtot+iStruct]->Draw("SAME");
 c1->Update();
 sprintf(txtPS,"Integral of gaussian fit wrt pT, all SMs");
 ps->TextNDC(0.2,0.94,txtPS);
 for (iStruct=0;iStruct<kNbZones;iStruct++) tgeNumberZones[kNbZones*kNbSMtot+iStruct]->Write();
 c1->cd(2);
 for (iStruct=0;iStruct<5;iStruct++) tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*kNbSMtot+iStruct]->SetName(Form("tgeNumberRatioZones_ALLSM%d",iStruct));
 tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*kNbSMtot+0]->Draw("AL");
 tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*kNbSMtot+0]->GetHistogram()->SetMinimum(0.0);
 tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*kNbSMtot+0]->GetHistogram()->SetMaximum(4.5);
 tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*kNbSMtot+0]->Draw("AL");
 for (iStruct=1;iStruct<5;iStruct++) tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*kNbSMtot+iStruct]->Draw("L SAME");
 for (iStruct=0;iStruct<5;iStruct++) tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*kNbSMtot+iStruct]->Write();
 c1->cd(3);
 lineMassH->SetLineColor(kMagenta);
 tgeMassZones[kNbZones*kNbSMtot+0]->Draw("AL");
 for (iStruct=1;iStruct<kNbZones;iStruct++) tgeMassZones[kNbZones*kNbSMtot+iStruct]->Draw("L SAME");
 lineMassH->Draw();
 c1->Update();
 sprintf(txtPS,"Mean of gaussian fit wrt pT, all SMs");
 ps->TextNDC(0.2,0.94,txtPS);
 for (iStruct=0;iStruct<kNbZones;iStruct++) tgeMassZones[kNbZones*kNbSMtot+iStruct]->Write();
 c1->cd(4);
 tgeSigZones[kNbZones*kNbSMtot+0]->Draw("AL");
 for (iStruct=1;iStruct<kNbZones;iStruct++) tgeSigZones[kNbZones*kNbSMtot+iStruct]->Draw("L SAME");
 c1->Update();
 sprintf(txtPS,"Width of gaussian fit wrt pT, all SMs");
 ps->TextNDC(0.2,0.94,txtPS);
 for (iStruct=0;iStruct<kNbZones;iStruct++) tgeSigZones[kNbZones*kNbSMtot+iStruct]->Write();
 c1->Update();
 
 // Peak characteristics in zones vs pT, per SM :
 for (iSM=0;iSM<kNbSMtot;iSM++)
    {if (tabChoiceCalos[SMdetType[iSM]] == 0)
        {printf("Skip SM %d of type %s.\n",iSM,detTypeString[SMdetType[iSM]]);
         continue;
         }
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,2);
     c1->cd(1);
     tgeNumberZones[kNbZones*iSM+3]->Draw("AL"); //3rd of these histoes is the highest.
     tgeNumberZones[kNbZones*iSM+3]->GetHistogram()->SetMinimum(0.);
     tgeNumberZones[kNbZones*iSM+3]->Draw("AL"); //3rd of these histoes is the highest.
     for (iStruct=0;iStruct<kNbZones;iStruct++) tgeNumberZones[kNbZones*iSM+iStruct]->Draw("L SAME");
     c1->Update();
     sprintf(txtPS,"Integral of gaussian fit wrt pT, SM %d",iSM);
     ps->TextNDC(0.2,0.94,txtPS);
     for (iStruct=0;iStruct<kNbZones;iStruct++) tgeNumberZones[kNbZones*iSM+iStruct]->Write();
     c1->cd(2);
     tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+0]->Draw("AL");
     tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+0]->GetHistogram()->SetMinimum(0.0);
     tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+0]->GetHistogram()->SetMaximum(4.5);
     tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+0]->Draw("AL");
     for (iStruct=1;iStruct<5;iStruct++) tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+iStruct]->Draw("L SAME");
     for (iStruct=0;iStruct<5;iStruct++) tgeNumberRatioZones[((int)(kNbZones/2))*(kNbZones-1)*iSM+iStruct]->Write();
     c1->cd(3);
     lineMassH->SetX1(tgeMass[3*iSM+0]->GetXaxis()->GetXmin());
     lineMassH->SetX2(tgeMass[3*iSM+0]->GetXaxis()->GetXmax());
     lineMassH->SetLineColor(kMagenta);
     tgeMassZones[kNbZones*iSM+0]->Draw("AL");
     for (iStruct=1;iStruct<kNbZones;iStruct++) tgeMassZones[kNbZones*iSM+iStruct]->Draw("L SAME");
     lineMassH->Draw();
     c1->Update();
     sprintf(txtPS,"Mean of gaussian fit wrt pT, SM %d",iSM);
     ps->TextNDC(0.2,0.94,txtPS);
     for (iStruct=0;iStruct<kNbZones;iStruct++) tgeMassZones[kNbZones*iSM+iStruct]->Write();
     c1->cd(4);
     tgeSigZones[kNbZones*iSM+0]->Draw("AL");
     for (iStruct=1;iStruct<kNbZones;iStruct++) tgeSigZones[kNbZones*iSM+iStruct]->Draw("L SAME");
     c1->Update();
     sprintf(txtPS,"Width of gaussian fit wrt pT, SM %d",iSM);
     ps->TextNDC(0.2,0.94,txtPS);
     for (iStruct=0;iStruct<kNbZones;iStruct++) tgeSigZones[kNbZones*iSM+iStruct]->Write();
     c1->Update();
     }
 
 
  
 
  




 
 
 //*****************************************
 // Loop over all towers :
 // Book histoes and calculate calibration coefficient :
 //*****************************************


 linePDGMass->SetLineColor(kSpring+4);
 lineFittedMass->SetLineColor(kMagenta);
 lineFittedWidth1->SetLineColor(kMagenta-10);
 lineFittedWidth2->SetLineColor(kMagenta-10);
 linePreviousFittedMass->SetLineColor(kMagenta);
 linePreviousFittedWidth1->SetLineColor(kMagenta-10);
 linePreviousFittedWidth2->SetLineColor(kMagenta-10);
 linePreviousDesiredMass->SetLineColor(kOrange+7);
 linePreviousUntrustedFit->SetLineColor(kOrange+7);
 linePreviousDesiredMass->SetLineWidth(2.);
 linePreviousUntrustedFit->SetLineWidth(2.);
 gStyle->SetOptStat(11); //Just name and nb entries.
 
 n=0;
 cmptEntriesDiscard=0;
 cmptAmpDiscard=0;
 cmptMasked=0;
 cmptEdge=0;
 printf("\n\nNow looping over towers.\nConsider %d towers, total gaussian amplitude %f.\nChanging gaussian amplitude par limit from %f to ",nbTowersConsideredTot,nbSig,maxGaussAmpl);
 maxGaussAmpl=50.*nbSig/(double)nbTowersConsideredTot;
 maxIntgAmpl=10.*nbTot/(double)nbTowersConsideredTot;
 printf("%f\n",maxGaussAmpl);
 fitfun->SetParLimits(0,8.,maxGaussAmpl);
 for(Int_t sm = 0; sm < kNbSMtot; sm++) //Loop on SM's.
    {if (tabChoiceCalos[SMdetType[sm]] == 0)
        {printf("Skip SM %d of type %s.\n",sm,detTypeString[SMdetType[sm]]);
         continue;
         }
     if (testChoice == 2)
        {if (!((sm == 0) || (sm == kNbSMEMCAL) || (sm == kNbSMEMCAL+kNbSMEMCALthird) || (sm == kNbSMEMCAL+kNbSMEMCALthird+kNbSMDCAL))) continue;
         }
     printf("\nSM %d : columns ",sm);
     fflush(stdout);
     for(Int_t c = kNbColMinEff ; c < TMath::Min(kNbColMaxEff,kTabNbCol[SMdetType[sm]]); c+=1) //Loop on columns.
        {printf("%d ",c);
         fflush(stdout);
         for(Int_t r = 0 ; r < kTabNbRow[SMdetType[sm]]; r+=1) //Loop on rows.
            {if (((r%12) == 0) && (testChoice != 3))
                {c1->Update();
                 ps->NewPage();
                 c1->Clear();
                 c1->Divide(3,4);
                 n=0;
                 }
             n++;
             cEff=c+kTabColOffset[sm]; //Column number used for (col,row) histoes must be pushed by 16 for DCAL SMs. //But don't use it in the end, because is the single-SM (col,row) plots are shifted by 16 to match the "visual view", then it's more difficult to find a tower which is everywhere in this code referenced by its offline coordinates (not shifted by 16). Thus shift only in hAllSpace* histoes in the end.
             if (testChoice != 3) c1->cd(n);
             //Read the values stored at the previous iteration :
             fscanf(txtFileParamsIn," %d %d %d %d ",&sm2,&c2,&r2,&cmpt2);
             for (k=0;k<kNbFitParams+kNbExtraParamsToBeRead;k++) fscanf(txtFileParamsIn,"%d ",flag2+k);
             //for (k=0;k<kNbTotParams-1;k++) fscanf(txtFileParamsIn,"%le ",param2+k); //Change to the following line, so as to be able to read files with fewer parameters written. NB : file always has 1 more than the sum : the calib coeff.
             for (k=0;k<kNbFitParams+kNbExtraParamsToBeRead;k++) fscanf(txtFileParamsIn,"%le ",param2+k);
             fscanf(txtFileParamsIn,"%le\n",param2+(kNbTotParams-1)); //To be sure to read the \n
             fscanf(txtFileCalibIn," %d %d %d %d %lf\n",&cmpt3,&sm3,&c3,&r3,&coeff3);
             tmpFlag=0;
             if (sm != sm3) tmpFlag++;
             if (c != c3) tmpFlag++;
             if (r != r3) tmpFlag++;
             if ((tmpFlag > 0) && (testChoice != 1)) //Shows up when one doesn't loop over all columns of a SM, and looks at more than the 1st SM.
                {if (testChoice == 0)
                    {printf("$$$ Problem sync input files : (%d,%d,%d) vs (%d,%d,%d) vs (%d,%d,%d)\n",sm,c,r,sm2,c2,r2,sm3,c3,r3);
                     }
                 }
             discardFlag=0;
             if (coeff3 < 0)
                {//Store the value of coeff3 "as a flag", then restore the value it should have as a coefficient :
                 discardFlag=(int)coeff3;
                 coeff3=1.00;
                 }
             //Get the histogram and setup the fit function :
             TH1F *h05      = (TH1F *) l05->FindObject(Form("%d_%d_%d",sm,c,r));
             h05->SetLineColor(colorIndex[1]);
             h05->Rebin(rebin);
             h05->SetXTitle("M_{#gamma,#gamma} (MeV/c^{2})");
             if (h05->GetEntries() != 0) hAllDistribNbEntries->Fill(h05->GetEntries());
             fitfun->SetParameter(0,h05->GetEntries()*(double)nbSig/(double)nbTot);
             if (fitfun->GetParameter(0) < 8.)
                {fitfun->SetParameter(0,8.); //With the new root version (today = Oct 5th 2012), this produces a run-time error "Error in ROOT::Math::ParameterSettings>: Invalid lower/upper bounds - ignoring the bounds" for the towers (80 for the 4 old SMs) with low number of entries. Caused by starting value of param outside the parameter limits set above.
                 //Now 07 Oct 2015, the error message changed to "lower/upper bounds outside current parameter value. The value will be set to (low+up)/2". Since it occurs for only 5 towers and doesn't prevent for a correct analysis, we just ignore.
                 }
             fitfun->SetParameter(1,mPDG);
             fitfun->SetParameter(2,typicalWidth);
             fitfun->SetParameter(3,fitPar3*(double)h05->GetEntries()/(double)nbTot);
             fitfun->SetParameter(4,fitPar4*(double)h05->GetEntries()/(double)nbTot);
             fitfun->SetParameter(5,fitPar5*(double)h05->GetEntries()/(double)nbTot);
             intgN=h05->Integral(1,h05->GetNbinsX());
             //When the histogram has entries, do the fit and calculate various khi2's :
             if (intgN != 0) // if intgN!=0
                {h05->Fit("fitfun","Q","",50,300);
                 fitfunBkgPol2->SetParameter(0,fitfun->GetParameter(3));
                 fitfunBkgPol2->SetParameter(1,fitfun->GetParameter(4));
                 fitfunBkgPol2->SetParameter(2,fitfun->GetParameter(5));
                 fitResultMean=fitfun->GetParameter(1);
                 fitResultSigma=fitfun->GetParameter(2);
                 cmptNDF=0;
                 fitKhi2SmallRange=0.;
                 for (iBin=h05->FindBin(fitResultMean-fitResultSigma);iBin<h05->FindBin(fitResultMean+fitResultSigma)+1;iBin++)
                    {if (h05->GetBinContent(iBin) == 0.) continue;
                     if (h05->GetBinError(iBin) == 0.)
                        {printf("## Weird bin %d (center %f) for tower (%d,%d,%d) : null error but content is %f : skip.\n",iBin,h05->GetBinCenter(iBin),sm,c,r,h05->GetBinContent(iBin));
                         continue;
                         }
                     fitKhi2SmallRange+=pow(h05->GetBinContent(iBin)-fitfun->Eval(h05->GetBinCenter(iBin)),2)/pow(h05->GetBinError(iBin),2);
                     cmptNDF++;
                     }
                 if (cmptNDF != 0) fitKhi2SmallRange=fitKhi2SmallRange/cmptNDF;
                 else fitKhi2SmallRange=0.;
                 cmptNDF=0;
                 fitKhi2PeakLeft=0.;
                 for (iBin=h05->FindBin(fitResultMean-3.*fitResultSigma);iBin<h05->FindBin(fitResultMean)+1;iBin++)
                    {if (h05->GetBinContent(iBin) == 0.) continue;
                     if (h05->GetBinError(iBin) == 0.)
                        {printf("## Weird bin %d (center %f) for tower (%d,%d,%d) : null error but content is %f : skip.\n",iBin,h05->GetBinCenter(iBin),sm,c,r,h05->GetBinContent(iBin));
                         continue;
                         }
                     fitKhi2PeakLeft+=pow(h05->GetBinContent(iBin)-fitfun->Eval(h05->GetBinCenter(iBin)),2)/pow(h05->GetBinError(iBin),2);
                     cmptNDF++;
                     }
                 if (cmptNDF != 0) fitKhi2PeakLeft=fitKhi2PeakLeft/cmptNDF;
                 else fitKhi2PeakLeft=0.;
                 cmptNDF=0;
                 fitKhi2PeakRight=0.;
                 for (iBin=h05->FindBin(fitResultMean);iBin<h05->FindBin(fitResultMean+3.*fitResultSigma)+1;iBin++)
                    {if (h05->GetBinContent(iBin) == 0.) continue;
                     if (h05->GetBinError(iBin) == 0.)
                        {printf("## Weird bin %d (center %f) for tower (%d,%d,%d) : null error but content is %f : skip.\n",iBin,h05->GetBinCenter(iBin),sm,c,r,h05->GetBinContent(iBin));
                         continue;
                         }
                     fitKhi2PeakRight+=pow(h05->GetBinContent(iBin)-fitfun->Eval(h05->GetBinCenter(iBin)),2)/pow(h05->GetBinError(iBin),2);
                     cmptNDF++;
                     }
                 if (cmptNDF != 0) fitKhi2PeakRight=fitKhi2PeakRight/cmptNDF;
                 else fitKhi2PeakRight=0.;
                 } //end if intgN!=0
             else //Do not fit if noisy/masked (empty) tower.
                {fitfun->SetParameter(0,-10.);
                 fitfun->SetParameter(1,0.);
                 fitfun->SetParameter(2,0.);
                 fitResultMean=0.;
                 fitResultSigma=0.;
                 fitfun->SetChisquare(0.);
                 cmptNDF=0;
                 fitKhi2SmallRange=0.;
                 fitKhi2PeakLeft=0.;
                 fitKhi2PeakRight=0.;
                 cmptMasked++;
                 hSpaceEntriesDiscard[1*sm+0]->SetBinContent(c+1,r+1,4);
                 } //end else (i.e. intgN==0).
             //Calculate the number of pi0's, by doing the integral of the gaussian peak, and by subtracting the bkg to the histogram integral (bin-counting) :
             if (h05->GetBinWidth(1) != 0.) intgS=fitResultSigma*fitfun->GetParameter(0)*TMath::Sqrt(2.*3.1415926535)/h05->GetBinWidth(1);
             else {printf("Zero width, put arbitrarily intgS to 0 for tower (%d;%d;%d).\n",sm,c,r); intgS=0.;}
             intgSbincounting=0.;
             for (iBin=h05->FindBin(fitResultMean-3.*fitResultSigma);iBin<h05->FindBin(fitResultMean+3.*fitResultSigma)+1;iBin++) intgSbincounting+=h05->GetBinContent(iBin)-fitfunBkgPol2->Eval(h05->GetBinCenter(iBin));
             //Calculate the calibration coefficient :
             coeff=pow(fitResultMean/mPDG,coeffPower);

          //Draw the inv mass plot together with the fit and several indicators (mass lines, flags,...), book histograms with the fit parameters :
          if (testChoice != 3)
          //if ((isFirstIteration == 0) && (testChoice != 3)) //Choose this if running DCAL with this has the code break.
             {h05->Write();
              h05->SetAxisRange(50.,250.,"X");
              h05->SetMinimum(0.);
              h05->SetStats();
              h05->Draw("HE");
              c1->Update();//Needed to get the pointer to the stat box, otherwise the pointer is null (the stat box is added to the list of object only when the histo is drawn, which, in compiled code, happens only when the canvas is updated).
              ptrPaveStat = (TPaveStats*)h05->FindObject("stats");
              ptrPaveStat->SetTextSize(0.04);
              //printf("%f %f %f %f\n",ptrPaveStat->GetX1NDC(),ptrPaveStat->GetX2NDC(),ptrPaveStat->GetY1NDC(),ptrPaveStat->GetY2NDC());
              //Default values are : 0.780 0.980 0.855 0.935.
              ptrPaveStat->SetX1NDC(0.75);
              ptrPaveStat->SetX2NDC(1.0);
              ptrPaveStat->SetY2NDC(0.95);
              fitfunBkgPol2->Draw("SAME");
              maxHisto=h05->GetMaximum()*(1+gStyle->GetHistTopMargin());
              //This comes from http://roottalk.root.cern.narkive.com/h0eAzqAn/taxis-strange-behaviour
              //h05->GetMaximum() is not satisfactory (the upper edge of the histo is not the max, and is not even always 1.05*max. Besides, since we are plotting in a subrange of the histo, the max of the histo may stand outside of the plotted range, and therefore the upper axis limit is not calculated wrt the histo's GetMaximum() anymore).
              //The Y axis really has only 1 bin, and its UpEdge is equal to GetXmax, even after the histo is plotted.
              //linePDGMass->SetY2(1.05*h05->GetMaximum());
              //printf("%f %f    %f  %d %d %f\n",h05->GetMaximum(),h05->GetYaxis()->GetXmax(),h05->GetMaximum()*(1+gStyle->GetHistTopMargin()),h05->GetYaxis()->GetNbins(),h05->GetYaxis()->GetLast(),h05->GetYaxis()->GetBinUpEdge(h05->GetYaxis()->GetLast()));
              linePDGMass->SetY2(1.05*maxHisto);
              linePDGMass->Draw("SAME");
              lineFittedMass->SetX1(fitResultMean);
              lineFittedMass->SetX2(fitResultMean);
              //lineFittedMass->SetY2(1.05*h05->GetMaximum());
              lineFittedMass->SetY2(0.85*maxHisto);
              lineFittedMass->Draw("SAME");
              lineFittedWidth1->SetX1(fitResultMean-2.*fitResultSigma);
              lineFittedWidth1->SetX2(fitResultMean-2.*fitResultSigma);
              //lineFittedWidth1->SetY2(0.8*h05->GetMaximum());
              lineFittedWidth1->SetY2(0.75*maxHisto);
              lineFittedWidth1->Draw("SAME");
              lineFittedWidth2->SetX1(fitResultMean+2.*fitResultSigma);
              lineFittedWidth2->SetX2(fitResultMean+2.*fitResultSigma);
              //lineFittedWidth2->SetY2(0.8*h05->GetMaximum());
              lineFittedWidth2->SetY2(0.75*maxHisto);
              lineFittedWidth2->Draw("SAME");
              paveUncertStatus->SetX1NDC(0.92); //We have to set them here again, apparently because root needs the plot to be drawn before it knows what NDCs mean.
              paveUncertStatus->SetX2NDC(1.0);
              paveUncertStatus->SetY1NDC(0.83);
              paveUncertStatus->SetY2NDC(0.75);
              paveMassStatus->SetX1NDC(0.94); //We have to set them here again, apparently because root needs the plot to be drawn before it knows what NDCs mean.
              paveMassStatus->SetX2NDC(0.98);
              paveMassStatus->SetY1NDC(0.81);
              paveMassStatus->SetY2NDC(0.77);
              paveUncertStatus->SetFillColor(kWhite);
              paveMassStatus->SetFillColor(kWhite);
              if (intgS > 0)
                 {valStatUncert=TMath::Sqrt(pow(paramFitStatUncert_b,2) + pow(paramFitStatUncert_a/TMath::Sqrt(intgS),2)); //This is the relative statistical uncertainty on the mean, in %.
                  //First check that the uncertainty on the mass is lower than the decalibration specifications :
                  if (valStatUncert > uncertKSpec*uncertSpec)
                     {paveUncertStatus->SetFillColor(kRed);
                      hSpaceCuts[5*sm+2]->Fill(c,r,4.);
                      }
                     else
                     {paveUncertStatus->SetFillColor(kGreen+1);
                      hSpaceCuts[5*sm+2]->Fill(c,r,1.9);
                      }
                  //Then check that the mass difference with the PDG is lower than uncertDist times the tolerated decalibration.
                  if (TMath::Abs(fitResultMean-mPDG)/mPDG > uncertDist*uncertSpec)
                     {paveMassStatus->SetFillColor(kRed);
                      hSpaceCuts[5*sm+4]->Fill(c,r,4.);
                      }
                     else
                     {paveMassStatus->SetFillColor(kGreen+1);
                      hSpaceCuts[5*sm+4]->Fill(c,r,1.9);
                      }
                  hAllDistribDistance[2]->Fill(TMath::Abs(fitResultMean-mPDG)/(fitResultMean*valStatUncert));
                  hSpaceCuts[5*sm+1]->Fill(c,r,100.*valStatUncert);
                  hAllDistribDistance2[2*0+0]->Fill(100.*valStatUncert);
                  if ((discardFlag >= 0) || (discardFlag <= -10)) hAllDistribDistance2[2*1+0]->Fill(100.*valStatUncert);
                  hSpaceCuts[5*sm+3]->Fill(c,r,TMath::Abs(fitResultMean-mPDG)/mPDG);
                  hAllDistribDistance2[2*0+1]->Fill(TMath::Abs(fitResultMean-mPDG)/mPDG);
                  if ((discardFlag >= 0) || (discardFlag <= -10)) hAllDistribDistance2[2*1+1]->Fill(TMath::Abs(fitResultMean-mPDG)/mPDG);
                  }
              paveUncertStatus->Draw("SAME");
              paveMassStatus->Draw("SAME");
              linePreviousFittedMass->SetX1(param2[1]);
              linePreviousFittedMass->SetX2(param2[1]);
              linePreviousFittedMass->SetY1(0.85*maxHisto);
              linePreviousFittedMass->SetY2(1.05*maxHisto);
              linePreviousFittedMass->Draw("SAME");
              linePreviousFittedWidth1->SetX1(param2[1]-2.*param2[2]);
              linePreviousFittedWidth1->SetX2(param2[1]-2.*param2[2]);
              linePreviousFittedWidth1->SetY1(0.85*maxHisto);
              linePreviousFittedWidth1->SetY2(maxHisto);
              linePreviousFittedWidth1->Draw("SAME");
              linePreviousFittedWidth2->SetX1(param2[1]+2.*param2[2]);
              linePreviousFittedWidth2->SetX2(param2[1]+2.*param2[2]);
              linePreviousFittedWidth2->SetY1(0.85*maxHisto);
              linePreviousFittedWidth2->SetY2(maxHisto);
              linePreviousFittedWidth2->Draw("SAME");
              linePreviousDesiredMass->SetX1(mPDG*pow(coeff3,1./coeffPowerPrev));
              linePreviousDesiredMass->SetX2(mPDG*pow(coeff3,1./coeffPowerPrev));
              linePreviousDesiredMass->SetY1(0.95*maxHisto);
              linePreviousDesiredMass->SetY2(1.05*maxHisto);
              linePreviousDesiredMass->Draw("SAME");
              if (coeff3 == 1.)
                 {linePreviousUntrustedFit->SetX1(mPDG-20.);
                  linePreviousUntrustedFit->SetX2(mPDG+20.);
                  linePreviousUntrustedFit->SetY1(maxHisto);
                  linePreviousUntrustedFit->SetY2(maxHisto);
                  linePreviousUntrustedFit->Draw("SAME");
                  }
              }//end if testChoice!=3
    
          //Book histograms with the integrals and khi2's : 
          for (j=0;j<kNbFitParams;j++) hSpace[kNbFitParams*sm+j]->Fill(c,r,fitfun->GetParameter(j));
          hSpaceIntg[kNbExtraParams*sm+0]->Fill(c,r,intgN);
          if (fitfun->GetNDF() > 0)
             {intgKhi2=fitfun->GetChisquare()/(double)fitfun->GetNDF();
              hSpaceIntg[kNbExtraParams*sm+1]->Fill(c,r,intgKhi2);
              }
             else
             {if (intgN > 6) printf("\n#### Weird NDF in tower (%d,%d,%d) : khi2 %f   NDF %f\n",sm,c,r,fitfun->GetChisquare(),(double)fitfun->GetNDF()); //NDF is 0 if more fit parameters (6) than filled bins. Don't want to calculate the nb of filled bins, so just check the number of entries.
              intgKhi2=0.;
              }
          hSpaceIntg[kNbExtraParams*sm+2]->Fill(c,r,intgS);
          hSpaceIntg[kNbExtraParams*sm+3]->Fill(c,r,intgN-intgS);
          hSpaceIntg[kNbExtraParams*sm+4]->Fill(c,r,fitKhi2SmallRange);
          hSpaceIntg[kNbExtraParams*sm+5]->Fill(c,r,fitKhi2PeakLeft);
          hSpaceIntg[kNbExtraParams*sm+6]->Fill(c,r,fitKhi2PeakRight);
          hSpaceFitErr[1*sm+0]->Fill(c,r,fitfun->GetParError(1));
          //Check if tower passes some quality cuts :
          cmpt=0;
          for (k=0;k<kNbFitParams+kNbExtraParams;k++) flag[k]=0;
          for (k=0;k<kNbFitParams;k++)
             {if (fitfun->GetParameter(k) < cutMin[k]) flag[k]++;
              if (fitfun->GetParameter(k) > cutMax[k]) flag[k]++;
              }
          if (intgN < cutMin[6]) flag[6]++;
          if (intgN > cutMax[6]) flag[6]++;
          if (intgKhi2 < cutMin[7]) flag[7]++;
          if (intgKhi2 > cutMax[7]) flag[7]++;
          if (intgS < cutMin[8]) flag[8]++;
          if (intgS > cutMax[8]) flag[8]++;
          if (fitKhi2SmallRange < cutMin[10]) flag[10]++;
          if (fitKhi2SmallRange > cutMax[10]) flag[10]++;
          for (k=0;k<kNbFitParams+kNbExtraParams;k++)
             {if (flag[k] != 0) cmpt++;
              }
          //Add coloured flags on the plot for the fit parameters :
          if (testChoice != 3)
             {for (k=0;k<kNbFitParams;k++)
                 {paveCutsStatus[k]->SetX1NDC(0.92+0.04*(int)(k/3));
                  paveCutsStatus[k]->SetX2NDC(0.96+0.04*(int)(k/3));
                  paveCutsStatus[k]->SetY1NDC(0.63-0.04*(k%3));
                  paveCutsStatus[k]->SetY2NDC(0.67-0.04*(k%3));
                  paveCutsStatus[k]->SetFillColor(19); //kGray is too dark and they haven't created a lighter version.
                  if (flag[k] != 0) paveCutsStatus[k]->SetFillColor(kRed);
                  //paveCutsStatus[k]->SetFillColor(k+2);
                  //paveCutsStatus[k]->SetBorderSize(2);
                  //paveCutsStatus[k]->SetLineColor(1);
                  //paveCutsStatus[k]->SetLineWidth(4); //Don't manage to have this work (draw a square with a border).
                  paveCutsStatus[k]->Draw("SAME");
                  }
              c1->Update(); //For the TLine's and TPave's. Nothing else is drawn beyond this line.
              }
          
          //Decide the output to the files as a function of whether the cuts were passed or not :
          flagBorderTower=1;
          if (c == 0) flagBorderTower=flagBorderTower*kCalibrateLastBorderCol;
          if (c == kTabNbCol[SMdetType[sm]]-1)
             {if (SMdetType[sm] == kDCAL) flagBorderTower=flagBorderTower*kCalibrateLastBorderCol;
                 else flagBorderTower=flagBorderTower*kCalibrateMidRapBorderCol;
              }
          if (r == 0)
             {if ((sm < 2) || ((sm >= kNbSMEMCAL+kNbSMEMCALthird) && (sm < kNbSMEMCAL+kNbSMEMCALthird+2))) flagBorderTower=flagBorderTower*kCalibrateLastBorderRow;
                 else flagBorderTower=flagBorderTower*kCalibrateTouchingBorderRow;
              }
          if (r == kTabNbRow[SMdetType[sm]]-1)
             {if (((sm >= kNbSMEMCAL+kNbSMEMCALthird-2) && (sm < kNbSMEMCAL+kNbSMEMCALthird)) || (sm >= kNbSMtot-2)) flagBorderTower=flagBorderTower*kCalibrateLastBorderRow;
                 else flagBorderTower=flagBorderTower*kCalibrateTouchingBorderRow;
              }
          if (flagBorderTower == 0) //This tower is on one of the edges that we don't want to calibrate.
             {//Use _coeffs flag 4. Put calib coeff -9 to flag them as "don't want to look".
              cmptEdge++;
              txtFileCalibOut << "4 " << sm << " " << c << " " << r << " -9" << endl;
              }//Remaining towers are to be looked at.
             else if (discardFlag < 0) //Decided at the previous iteration that we don't want to look at his tower anymore.
             {//Just copy the same _coeffs flag "cmpt3", and put calib coeff 1.00.
              txtFileCalibOut << cmpt3 << " " << sm << " " << c << " " << r << " " << discardFlag << endl;
              }//Remaining towers are to be looked at.
              else if (cmpt == 0) //This tower is a good tower.
             {//Use _coeffs flag 0. Put _params flag "cmpt" -1.
              cmpt=-1;
              txtFileCalibOut << "0 " << sm << " " << c << " " << r << " " << coeff << endl;
              hSpaceCoeff[1*sm+0]->Fill(c,r,coeff);
              hCoeff->Fill(coeff);
              }//Remaining towers failed at least one cut.
              else if (intgN == 0) //Empty tower (was noisy/masked) (no counts in histo).
             {//Put _params flag "cmpt" 0.
              cmpt=0; //Doesn't show up in hSpaceCuts anyway
              txtFileCalibOut << "0 " << sm << " " << c << " " << r << " 1" << endl;
              }//Remaining towers failed at least a cut and are not empty.
              else if (h05->GetEntries() > cutEntriesDiscard) //Normal number of entries.
             {//Use _coeffs flag 1. Keep _params flag "cmpt" = nb of failed quality cuts.
              txtFileCalibOut << "1 " << sm << " " << c << " " << r << " " << coeff << " //Signif=";
              //txtFileCalibOut << "1 " << sm << " " << c << " " << r << " " << coeff << " //Signif=";
              //txtFileCalibOut << intgS/(4.*10.*(fitfun->GetParameter(3)+mPDG*fitfun->GetParameter(4)+fitfun->GetParameter(5)*(mPDG*mPDG+4.*10.*10./3.))) << " //Dist="; //Ca c'est S/N, now S/sqrt(S+N) :
              if (fitfun->Integral(mPDG-2.*10.,mPDG+2.*10.) > 0)
                 txtFileCalibOut << intgS/TMath::Sqrt(fitfun->Integral(mPDG-2.*10.,mPDG+2.*10.)) << " //Dist=";
                 else txtFileCalibOut << 0.00 << " //Dist=";
              if (fitfun->GetParError(1) != 0.)
                 txtFileCalibOut << TMath::Abs(fitfun->GetParameter(1)-mPDG)/fitfun->GetParError(1) << " //Amp=";
                 else txtFileCalibOut << 9999999999. << " //Amp=";
              if (flag[0] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << fitfun->GetParameter(0) << " //Mu=";
              if (flag[1] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << fitfun->GetParameter(1) << " //Sig=";
              if (flag[2] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << fitfun->GetParameter(2) << " //c=";
              if (flag[3] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << fitfun->GetParameter(3) << " //b=";
              if (flag[4] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << fitfun->GetParameter(4) << " //a=";
              if (flag[5] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << fitfun->GetParameter(5) << " //I=";
              if (flag[6] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << intgN << " //KhiSqNdf=";
              if (flag[7] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << intgKhi2 << " //S=";
              if (flag[8] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << intgS << " //Khi2Peak=";
              if (flag[10] != 0) txtFileCalibOut << "**";
              txtFileCalibOut << fitKhi2SmallRange;
              txtFileCalibOut << endl;
              }//Remaining towers failed the cuts because they have a small number of entries.
              else //Too low number of entries : problem with this tower. Put automatically calib coeff to 1.00 and drop its coords in txt and ps.
             {//Use _coeffs flag 2. Keep _params flag "cmpt" = nb of failed quality cuts.
              txtFileCalibOut << "2 " << sm << " " << c << " " << r << " " << "1.00" << endl;
              cmptEntriesDiscard++;
              hSpaceEntriesDiscard[1*sm+0]->SetBinContent(c+1,r+1,0.8);
              }
          
          //Now we are back in the general loop for all the towers.
          //Write the value of all parameters to the parameter file :
          fprintf(txtFileParamsOut,"%d %d %d %d ",sm,c,r,cmpt);
          for (k=0;k<kNbFitParams+kNbExtraParams;k++) fprintf(txtFileParamsOut,"%d ",flag[k]);
          for (k=0;k<kNbFitParams;k++) fprintf(txtFileParamsOut,"%e ",fitfun->GetParameter(k));
          fprintf(txtFileParamsOut,"%e %e %e %e %e %e %e %e\n",intgN,intgKhi2,intgS,intgN-intgS,fitKhi2SmallRange,fitKhi2PeakLeft,fitKhi2PeakRight,coeff);
          //Fill in histograms :
          for (k=0;k<kNbFitParams+kNbExtraParamsToBeRead;k++) hAllOldDistrib[k]->Fill(param2[k]);
          hAllOldDistribMass->Fill(param2[1]);
          hAllOldDistribMassPerSM[sm]->Fill(param2[1]);
          hCoeffOld->Fill(param2[kNbTotParams-1]);
          hCoeffOldLarge->Fill(param2[kNbTotParams-1]);
          hCoeffOldCorr->Fill(coeff3);
          hCoeffOldCorrLarge->Fill(coeff3);
          if ((isFirstIteration == 0) && (coeff3 == 1.00) && (hSpaceEntriesDiscard[1*sm+0]->GetBinContent(c+1,r+1) == 0.)) hSpaceEntriesDiscard[1*sm+0]->SetBinContent(c+1,r+1,1.9);
          //Difference between both histoes : param2[kNbTotParams-1] is the raw coeff, while coeff3 comes from the "clean" coeff file (hand-corrected).
          hSpaceCuts[5*sm+0]->Fill(c,r,cmpt);
          hCorrelMuVsA->Fill(fitfun->GetParameter(0),fitfun->GetParameter(1));
          hCorrelSigVsA->Fill(fitfun->GetParameter(0),fitfun->GetParameter(2));
          hCorrelSigVsMu->Fill(fitfun->GetParameter(1),fitfun->GetParameter(2));
          hCorrelBVsA->Fill(fitfun->GetParameter(3),fitfun->GetParameter(4));
          hCorrelCVsA->Fill(fitfun->GetParameter(3),fitfun->GetParameter(5));
          hCorrelCVsB->Fill(fitfun->GetParameter(4),fitfun->GetParameter(5));
          hCorrelISVsI->Fill(intgN,intgN-intgS);
          hCorrelSVsI->Fill(intgN,intgS);
          hCorrelSVsIS->Fill(intgN-intgS,intgS);
          hCorrelSbincountingVsS->Fill(intgS,intgSbincounting);
          hCorrelMuVsI->Fill(intgN,fitfun->GetParameter(1));
          hCorrelMuVsS->Fill(intgS,fitfun->GetParameter(1));
          hCorrelSigVsI->Fill(intgN,fitfun->GetParameter(2));
          hCorrelSigVsS->Fill(intgS,fitfun->GetParameter(2));
          if (h05->GetEntries() != 0) hAllDistribNbEntriesCorrel->Fill(intgN,h05->GetEntries());
          for (k=0;k<3;k++)
             {if (param2[k] != 0.) hAllDiffAllTw[k]->Fill(fitfun->GetParameter(k)/param2[k]);
              }
          for (k=3;k<6;k++) hAllDiffAllTw[k]->Fill(fitfun->GetParameter(k)-param2[k]);
          if (param2[6] != 0.) hAllDiffAllTw[6]->Fill(intgN/param2[6]);
          if (param2[7] != 0.) hAllDiffAllTw[7]->Fill(intgKhi2/param2[7]);
          if (param2[8] != 0.) hAllDiffAllTw[8]->Fill(intgS/param2[8]);
          if (param2[9] != 0.) hAllDiffAllTw[9]->Fill((intgN-intgS)/param2[9]);
          if (coeff3 != 1.)
             {hCtrlEvolCoeff->Fill((coeff-1.)/(coeff3-1.));
              h2CtrlEvolCoeff1->Fill(coeff3,(coeff-1.)/(coeff3-1.));
              }
          h2CtrlEvolCoeff2->Fill(coeff3-1.,coeff3-coeff);
          if ((cmpt == -1) && (cmpt2 == -1)) //Book only if we are sure the tower is valid. NB : most of the discarded towers are valid too...
             {if (param2[1] != mPDG)
                 {hCtrlEvolMu->Fill((fitfun->GetParameter(1)-mPDG)/(param2[1]-mPDG));
                  h2CtrlEvolMu1->Fill(param2[1],(fitfun->GetParameter(1)-mPDG)/(param2[1]-mPDG));
                  }
              h2CtrlEvolMu2->Fill((param2[1]-mPDG)/mPDG,(param2[1]-fitfun->GetParameter(1))/mPDG);
              for (k=0;k<3;k++)
                 {if (param2[k] != 0.)
                     {hSpaceDiff[kNbTotParams*sm+k]->Fill(c,r,fitfun->GetParameter(k)/param2[k]);
                      hAllDiff[k]->Fill(fitfun->GetParameter(k)/param2[k]);
                      }
                  }
              for (k=3;k<6;k++)
                 {hSpaceDiff[kNbTotParams*sm+k]->Fill(c,r,fitfun->GetParameter(k)-param2[k]);
                  hAllDiff[k]->Fill(fitfun->GetParameter(k)-param2[k]);
                  }
              if (param2[6] != 0.)
                 {hSpaceDiff[kNbTotParams*sm+6]->Fill(c,r,intgN/param2[6]);
                  hAllDiff[6]->Fill(intgN/param2[6]);
                  }
              if (param2[7] != 0.)
                 {hSpaceDiff[kNbTotParams*sm+7]->Fill(c,r,intgKhi2/param2[7]);
                  hAllDiff[7]->Fill(intgKhi2/param2[7]);
                  }
              if (param2[8] != 0.)
                 {hSpaceDiff[kNbTotParams*sm+8]->Fill(c,r,intgS/param2[8]);
                  hAllDiff[8]->Fill(intgS/param2[8]);
                  }
              if (param2[9] != 0.)
                 {hSpaceDiff[kNbTotParams*sm+9]->Fill(c,r,(intgN-intgS)/param2[9]);
                  hAllDiff[9]->Fill((intgN-intgS)/param2[9]);
                  }
              }
        }//End loop on rows.
      }//End loop on columns.
      
      c1->Update();
    
    printf("\n"); 
      
    }//End loop on SM's.
 gStyle->SetOptStat(1111); //Back to the standard style (name, nb entries, mean, RMS).
 


 








 
 //*****************************************
 // Draw all remaining histoes :
 //*****************************************

 printf("\n\nThere are %d masked/empty/dead towers.\n\nList of the %d towers of which calib coeff has automatically been put to 1.00 because nbEntries < %d :\n",cmptMasked,cmptEntriesDiscard,(int)cutEntriesDiscard);
 fflush(stdout);
 gSystem->Exec(Form("less %s | grep -E '^2 ' | sed 's/^2 //g' | sed 's/ 1.0$//g'",txtFileCalibOutName.Data()));
 
 printf("\n\nList of the %d towers of which calib coeff has automatically been put to 1.0 because %f <Amp < %f :\n",cmptAmpDiscard,8.,cutAmpDiscard);
 fflush(stdout);
 gSystem->Exec(Form("less %s | grep -E '^3 ' | sed 's/^3 //g' | sed 's/ 1.0$//g'",txtFileCalibOutName.Data()));

 printf("\n\nList of the %d towers of which calib coeff has automatically been put to 1.0 because it's on a SM edge :\n",cmptEdge);
 fflush(stdout);
 gSystem->Exec(Form("less %s | grep -E '^4 ' | sed 's/^4 //g' | sed 's/ -9$//g'",txtFileCalibOutName.Data()));

 printf("\n\nStart drawing final plots...\n\n");
 

 //Spatial distribution per SM :
 for (i=0;i<kNbSMtot;i++)
    {if ((i<kNbSMtot) && (tabChoiceCalos[SMdetType[i]] == 0))
        {printf("Skip SM %d of type %s.\n",i,detTypeString[SMdetType[i]]);
         continue;
         }
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,(int)((kNbFitParams+1)/2));
     hSpace[kNbFitParams*i+0]->SetMinimum(0.);
     if (hSpace[kNbFitParams*i+0]->GetMaximum() > 3.0*tabMax[0]) hSpace[kNbFitParams*i+0]->SetMaximum(3.0*tabMax[0]);
     nbTot=hSpace[kNbFitParams*i+1]->GetMaximum();
     for (k=0;k<hSpace[kNbFitParams*i+1]->GetNbinsX();k++)
        {for (j=0;j<hSpace[kNbFitParams*i+1]->GetNbinsY();j++)
            {if ((hSpace[kNbFitParams*i+1]->GetBinContent(k+1,j+1) < nbTot) && (hSpace[kNbFitParams*i+1]->GetBinContent(k+1,j+1) > 0)) nbTot=hSpace[kNbFitParams*i+1]->GetBinContent(k+1,j+1);
             }
         }
     hSpace[kNbFitParams*i+1]->SetMinimum(0.98*nbTot);
     if (hSpace[kNbFitParams*i+1]->GetMaximum() > 190.) hSpace[kNbFitParams*i+1]->SetMaximum(190.);
     hSpace[kNbFitParams*i+2]->SetMinimum(tabMin[2]);
     if (hSpace[kNbFitParams*i+3]->GetMaximum() > 3.0*tabMax[0]) hSpace[kNbFitParams*i+3]->SetMaximum(3.0*tabMax[0]);
     if (hSpace[kNbFitParams*i+4]->GetMaximum() > 3.0*tabMax[0]) hSpace[kNbFitParams*i+4]->SetMaximum(3.0*tabMax[0]);
     if (hSpace[kNbFitParams*i+5]->GetMaximum() > 3.0*tabMax[0]) hSpace[kNbFitParams*i+5]->SetMaximum(3.0*tabMax[0]);
     for (j=0;j<kNbFitParams;j++)
        {c1->cd(j+1);
         k=(j%2)*(int)((j+5)/2)+(1-(j%2))*(int)(j/2);
         hSpace[kNbFitParams*i+k]->Draw("COLZ");
         hSpaceEntriesDiscard[1*i+0]->SetMarkerColor(kWhite); //This changes also the text color for the Draw("TEXT") !!
         //Also changes size and color of test in ps->TextNDC !! => SetTxtColor necessary. But no function to set the font...
         //NB : SetMarkerColor has the same effect if placed when the histo is instanciated.
         hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
         //c1->Update();
         //hSpaceEntriesDiscard[1*i+0]->SetMarkerSize(3);
         hSpaceEntriesDiscardHack[1*i+0]->Draw("TEXT SAME");
         //gStyle->SetTextSize(0.5);
         //hSpaceEntriesDiscard[1*i+0]->UseCurrentStyle();
         c1->Update();
         sprintf(txtPS,"%s for SM %d",varNameLong[k],i);
         //H mettre a la fin, sinon est overwrite par le plot (=> disparait).
         //The histo hSpaceEntriesDiscardHack is needed to have the figures of hSpaceEntriesDiscard plotted in small font, while the ps text (TextNDC) is plotted in readable font size.
         ps->SetTextColor(1);
         ps->TextNDC(0.2,0.94,txtPS);
         hSpace[kNbFitParams*i+k]->Write();
         }
     c1->Update();

     for (j=0;j<kNbExtraParams;j++)
        {hSpaceIntg[kNbExtraParams*i+j]->Write();
         }
     //if (hSpaceIntg[kNbExtraParams*i+0]->GetMaximum() > maxIntgAmpl) hSpaceIntg[kNbExtraParams*i+0]->SetMaximum(maxIntgAmpl);
     //Take the previous line out, to be able to catch the noisy towers.
     //if (hSpaceIntg[kNbExtraParams*i+1]->GetMaximum() > maxIntgAmpl) hSpaceIntg[kNbExtraParams*i+1]->SetMaximum(maxIntgAmpl);
     //Take the previous line out, because now that histo stores the khi2.
     hSpaceIntg[kNbExtraParams*i+0]->SetMinimum(0.);
     hSpaceIntg[kNbExtraParams*i+1]->SetMinimum(0.);
     hSpaceIntg[kNbExtraParams*i+2]->SetMinimum(0.);
     hSpaceIntg[kNbExtraParams*i+3]->SetMinimum(0.);
     if (hSpaceIntg[kNbExtraParams*i+2]->GetMaximum() > 3.0*tabMax[8]) hSpaceIntg[kNbExtraParams*i+2]->SetMaximum(3.0*tabMax[8]);
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hSpaceIntg[kNbExtraParams*i+0]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     hSpaceEntriesDiscardHack[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d",varNameLong[0+kNbFitParams],i);
     ps->SetTextColor(1);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     c1->cd(3);
     hSpaceIntg[kNbExtraParams*i+2]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     hSpaceEntriesDiscardHack[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d",varNameLong[2+kNbFitParams],i);
     ps->SetTextColor(1);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     c1->cd(5);
     hSpaceIntg[kNbExtraParams*i+3]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     hSpaceEntriesDiscardHack[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d",varNameLong[3+kNbFitParams],i);
     ps->SetTextColor(1);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     hSpaceIntg[kNbExtraParams*i+0]->SetMaximum(1.5*tabMax[kNbFitParams+0]);
     c1->cd(2);
     hSpaceIntg[kNbExtraParams*i+0]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d (zoomed)",varNameLong[0+kNbFitParams],i);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
 
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hSpaceIntg[kNbExtraParams*i+1]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     hSpaceEntriesDiscardHack[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d",varNameLong[1+kNbFitParams],i);
     ps->SetTextColor(1);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     c1->cd(3);
     hSpaceIntg[kNbExtraParams*i+4]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     hSpaceEntriesDiscardHack[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d",varNameLong[4+kNbFitParams],i);
     ps->SetTextColor(1);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     hSpaceIntg[kNbExtraParams*i+1]->SetMaximum(5.0);
     hSpaceIntg[kNbExtraParams*i+4]->SetMaximum(5.0);
     hSpaceIntg[kNbExtraParams*i+5]->SetMaximum(5.0);
     hSpaceIntg[kNbExtraParams*i+6]->SetMaximum(5.0);
     c1->cd(2);
     hSpaceIntg[kNbExtraParams*i+1]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d (zoomed)",varNameLong[1+kNbFitParams],i);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(4);
     hSpaceIntg[kNbExtraParams*i+4]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d (zoomed)",varNameLong[4+kNbFitParams],i);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(5);
     hSpaceIntg[kNbExtraParams*i+5]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     hSpaceEntriesDiscardHack[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d",varNameLong[5+kNbFitParams],i);
     ps->SetTextColor(1);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(6);
     hSpaceIntg[kNbExtraParams*i+6]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     hSpaceEntriesDiscardHack[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"%s for SM %d",varNameLong[6+kNbFitParams],i);
     ps->SetTextColor(1);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();

     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hSpaceCuts[5*i+1]->SetMinimum(0.);
     hSpaceCuts[5*i+1]->SetMaximum(2.);
     hSpaceCuts[5*i+1]->Draw("COLZ");
     c1->Update();
     sprintf(txtPS,"Stat uncert on mean for SM %d",i);
     ps->TextNDC(0.2,0.94,txtPS);
     hSpaceCuts[5*i+1]->Write();
     c1->Update();
     c1->cd(2);
     hSpaceCuts[5*i+2]->SetMinimum(0.);
     hSpaceCuts[5*i+2]->SetMaximum(5.05);
     hSpaceCuts[5*i+2]->Draw("COLZ");
     c1->Update();
     sprintf(txtPS,"Stat uncert on mean for SM %d",i);
     ps->TextNDC(0.2,0.94,txtPS);
     hSpaceCuts[5*i+2]->Write();
     c1->Update();
     c1->cd(3);
     hSpaceCuts[5*i+3]->SetMinimum(0.);
     hSpaceCuts[5*i+3]->SetMaximum(0.08);
     hSpaceCuts[5*i+3]->Draw("COLZ");
     c1->Update();
     sprintf(txtPS,"Relative distance to mPDG for SM %d",i);
     ps->TextNDC(0.2,0.94,txtPS);
     hSpaceCuts[5*i+3]->Write();
     c1->Update();
     c1->cd(4);
     hSpaceCuts[5*i+4]->SetMinimum(0.);
     hSpaceCuts[5*i+4]->SetMaximum(5.05);
     hSpaceCuts[5*i+4]->Draw("COLZ");
     c1->Update();
     sprintf(txtPS,"Relative distance to mPDG for SM %d",i);
     ps->TextNDC(0.2,0.94,txtPS);
     hSpaceCuts[5*i+4]->Write();
     c1->Update();
     c1->cd(6);
     hSpaceCuts[5*i+0]->SetMinimum(0.);
     hSpaceCuts[5*i+0]->SetMaximum(kNbFitParams+4);
     hSpaceCuts[5*i+0]->Draw("COLZ");
     hSpaceEntriesDiscard[1*i+0]->Draw("TEXT SAME");
     c1->Update();
     sprintf(txtPS,"Untrusted fits for SM %d",i);
     ps->TextNDC(0.2,0.94,txtPS);
     hSpaceCuts[5*i+0]->Write();
     c1->Update();
     }
 

 //Now put all spatial distribs together to have a full calorimeter real view :
 //The way the loops are imbricated is not easily readable but should be the fastest for processing.
 for (j=0;j<kNbFitParams;j++)
    {for (i=0;i<kNbSMEMCAL+kNbSMEMCALthird;i++)
        {if ((i%2) == 0) shiftCol=0;
            else shiftCol=kNbColMax;
         shiftRow=(int)(i/2)*kNbRowMax;
         for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {hAllSpaceEMCAL[j]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,hSpace[kNbFitParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 //Book in 1-D histoes only if tower not masked/empty ; no way to know it at this stage except that mean and sigma are both 0 :
                 if ((hSpace[kNbFitParams*i+1]->GetBinContent(iCol+1,iRow+1) != 0.) && (hSpace[kNbFitParams*i+2]->GetBinContent(iCol+1,iRow+1) != 0.)) hAllDistrib[j]->Fill(hSpace[kNbFitParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 if (j == 1)
                    {hAllDistribMass->Fill(hSpace[kNbFitParams*i+j]->GetBinContent(iCol+1,iRow+1));
                     hAllDistribMassPerSM[i]->Fill(hSpace[kNbFitParams*i+j]->GetBinContent(iCol+1,iRow+1));
                     }
                 }
             }
         }
     }
 for (j=0;j<kNbFitParams;j++)
    {for (i=kNbSMEMCAL+kNbSMEMCALthird;i<kNbSMtot;i++)
        {if ((i%2) == 0) shiftCol=0;
            else shiftCol=kNbColMax+kTabColOffset[i];
         shiftRow=(int)((i-(kNbSMEMCAL+kNbSMEMCALthird))/2)*kNbRowMax;
         for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {hAllSpaceDCAL[j]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,hSpace[kNbFitParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 //Book in 1-D histoes only if tower not masked/empty ; no way to know it at this stage except that mean and sigma are both 0 :
                 if ((hSpace[kNbFitParams*i+1]->GetBinContent(iCol+1,iRow+1) != 0.) && (hSpace[kNbFitParams*i+2]->GetBinContent(iCol+1,iRow+1) != 0.)) hAllDistrib[j]->Fill(hSpace[kNbFitParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 if (j == 1)
                    {hAllDistribMass->Fill(hSpace[kNbFitParams*i+j]->GetBinContent(iCol+1,iRow+1));
                     hAllDistribMassPerSM[i]->Fill(hSpace[kNbFitParams*i+j]->GetBinContent(iCol+1,iRow+1));
                     }
                 }
             }
         }
     }
 for (j=0;j<kNbExtraParams;j++)
    {for (i=0;i<kNbSMEMCAL+kNbSMEMCALthird;i++)
        {if ((i%2) == 0) shiftCol=0;
            else shiftCol=kNbColMax;
         shiftRow=(int)(i/2)*kNbRowMax;
         for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {hAllSpaceEMCALIntg[j]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,hSpaceIntg[kNbExtraParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 //Book in 1-D histoes only if tower not masked/empty ; no way to know it at this stage except that mean and sigma are both 0 :
                 if ((hSpace[kNbFitParams*i+1]->GetBinContent(iCol+1,iRow+1) != 0.) && (hSpace[kNbFitParams*i+2]->GetBinContent(iCol+1,iRow+1) != 0.)) hAllDistribIntg[j]->Fill(hSpaceIntg[kNbExtraParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 }
             }
         }
     }
 for (j=0;j<kNbExtraParams;j++)
    {for (i=kNbSMEMCAL+kNbSMEMCALthird;i<kNbSMtot;i++)
        {if ((i%2) == 0) shiftCol=0;
            else shiftCol=kNbColMax+kTabColOffset[i];
         shiftRow=(int)((i-(kNbSMEMCAL+kNbSMEMCALthird))/2)*kNbRowMax;
         for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {hAllSpaceDCALIntg[j]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,hSpaceIntg[kNbExtraParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 //Book in 1-D histoes only if tower not masked/empty ; no way to know it at this stage except that mean and sigma are both 0 :
                 if ((hSpace[kNbFitParams*i+1]->GetBinContent(iCol+1,iRow+1) != 0.) && (hSpace[kNbFitParams*i+2]->GetBinContent(iCol+1,iRow+1) != 0.)) hAllDistribIntg[j]->Fill(hSpaceIntg[kNbExtraParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 }
             }
         }
     }
 for (j=0;j<1;j++)
    {for (i=0;i<kNbSMEMCAL+kNbSMEMCALthird;i++)
        {if ((i%2) == 0) shiftCol=0;
            else shiftCol=kNbColMax;
         shiftRow=(int)(i/2)*kNbRowMax;
         for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {hAllSpaceEMCALCoeff[j]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,hSpaceCoeff[1*i+j]->GetBinContent(iCol+1,iRow+1));
                 }
             }
         }
     }
 for (j=0;j<1;j++)
    {for (i=kNbSMEMCAL+kNbSMEMCALthird;i<kNbSMtot;i++)
        {if ((i%2) == 0) shiftCol=0;
            else shiftCol=kNbColMax+kTabColOffset[i];
         shiftRow=(int)((i-(kNbSMEMCAL+kNbSMEMCALthird))/2)*kNbRowMax;
         for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {hAllSpaceDCALCoeff[j]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,hSpaceCoeff[1*i+j]->GetBinContent(iCol+1,iRow+1));
                 }
             }
         }
     }
 for (i=0;i<kNbSMEMCAL+kNbSMEMCALthird;i++)
    {if ((i%2) == 0) shiftCol=0;
        else shiftCol=kNbColMax;
     shiftRow=(int)(i/2)*kNbRowMax;
     for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
        {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
            {//hAllSpaceEMCALDistance[0] = dist_N = (mu_fit - mPDG)/(sig_fit/sqrt(intgS))
             //hAllSpaceEMCALDistance[1] = dist_fit = (mu_fit - mPDG) / (fit uncertainty on mu_fit)
             //hAllSpaceEMCALDistance[2] = (dist_N - dist_fit) / dist_fit
             if ((hAllSpaceEMCAL[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow) == 0.) && (hAllSpaceEMCAL[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow) == 0.)) continue; //A noisy tower, do nothing.
             hAllSpaceEMCALDistance[1]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,TMath::Abs(hAllSpaceEMCAL[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow)-mPDG)/hSpaceFitErr[1*i+0]->GetBinContent(iCol+1,iRow+1));
             if (hAllSpaceEMCALIntg[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow) != 0.) hAllSpaceEMCALDistance[0]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,TMath::Abs(hAllSpaceEMCAL[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow)-mPDG)/(hAllSpaceEMCAL[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow)/TMath::Sqrt(hAllSpaceEMCALIntg[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow))));
                else hAllSpaceEMCALDistance[0]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,49.);
             if (hAllSpaceEMCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow) != 0.) hAllSpaceEMCALDistance[2]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,(hAllSpaceEMCALDistance[0]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow)-hAllSpaceEMCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow))/hAllSpaceEMCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             hAllDistribDistance[0]->Fill(hAllSpaceEMCALDistance[0]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             hAllDistribDistance[1]->Fill(hAllSpaceEMCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             hCorrelDvsSig->Fill(hAllSpaceEMCAL[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow),hAllSpaceEMCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             hCorrelCmpD->Fill(hAllSpaceEMCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow),hAllSpaceEMCALDistance[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             }
         }
     }
 for (i=kNbSMEMCAL+kNbSMEMCALthird;i<kNbSMtot;i++)
    {if ((i%2) == 0) shiftCol=0;
        else shiftCol=kNbColMax+kTabColOffset[i];
     shiftRow=(int)((i-(kNbSMEMCAL+kNbSMEMCALthird))/2)*kNbRowMax;
     for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
        {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
            {//hAllSpaceDCALDistance[0] = dist_N = (mu_fit - mPDG)/(sig_fit/sqrt(intgS))
             //hAllSpaceDCALDistance[1] = dist_fit = (mu_fit - mPDG) / (fit uncertainty on mu_fit)
             //hAllSpaceDCALDistance[2] = (dist_N - dist_fit) / dist_fit
             if ((hAllSpaceDCAL[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow) == 0.) && (hAllSpaceDCAL[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow) == 0.)) continue; //A noisy tower, do nothing.
             hAllSpaceDCALDistance[1]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,TMath::Abs(hAllSpaceDCAL[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow)-mPDG)/hSpaceFitErr[1*i+0]->GetBinContent(iCol+1,iRow+1));
             if (hAllSpaceDCALIntg[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow) != 0.) hAllSpaceDCALDistance[0]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,TMath::Abs(hAllSpaceDCAL[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow)-mPDG)/(hAllSpaceDCAL[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow)/TMath::Sqrt(hAllSpaceDCALIntg[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow))));
                else hAllSpaceDCALDistance[0]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,49.);
             if (hAllSpaceDCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow) != 0.) hAllSpaceDCALDistance[2]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,(hAllSpaceDCALDistance[0]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow)-hAllSpaceDCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow))/hAllSpaceDCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             hAllDistribDistance[0]->Fill(hAllSpaceDCALDistance[0]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             hAllDistribDistance[1]->Fill(hAllSpaceDCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             hCorrelDvsSig->Fill(hAllSpaceDCAL[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow),hAllSpaceDCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             hCorrelCmpD->Fill(hAllSpaceDCALDistance[1]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow),hAllSpaceDCALDistance[2]->GetBinContent(iCol+1+shiftCol,iRow+1+shiftRow));
             }
         }
     }
 for (j=0;j<kNbTotParams;j++)
    {for (i=0;i<kNbSMEMCAL+kNbSMEMCALthird;i++)
        {if ((i%2) == 0) shiftCol=0;
            else shiftCol=kNbColMax;
         shiftRow=(int)(i/2)*kNbRowMax;
         for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {hAllSpaceEMCALDiff[j]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,hSpaceDiff[kNbTotParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 }
             }
         }
     }
 for (j=0;j<kNbTotParams;j++)
    {for (i=kNbSMEMCAL+kNbSMEMCALthird;i<kNbSMtot;i++)
        {if ((i%2) == 0) shiftCol=0;
            else shiftCol=kNbColMax+kTabColOffset[i];
         shiftRow=(int)((i-(kNbSMEMCAL+kNbSMEMCALthird))/2)*kNbRowMax;
         for (iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {hAllSpaceDCALDiff[j]->SetBinContent(iCol+1+shiftCol,iRow+1+shiftRow,hSpaceDiff[kNbTotParams*i+j]->GetBinContent(iCol+1,iRow+1));
                 }
             }
         }
     }

 if (choiceNoEMCAL != 0)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,(int)((kNbFitParams+1)/2));
     hAllSpaceEMCAL[0]->SetMinimum(0.);
     hAllSpaceEMCAL[0]->SetMaximum(1.2*tabMax[0]);
     hAllSpaceEMCAL[1]->SetMinimum(tabMin[1]+30.);
     hAllSpaceEMCAL[1]->SetMaximum(tabMax[1]-10.);
     hAllSpaceEMCAL[2]->SetMinimum(tabMin[2]+5.);
     hAllSpaceEMCAL[2]->SetMaximum(tabMax[2]-5.);
     for (j=0;j<kNbFitParams;j++)
        {c1->cd(j+1);
         k=(j%2)*(int)((j+5)/2)+(1-(j%2))*(int)(j/2);
         hAllSpaceEMCAL[k]->Draw("COLZ");
         lineSMborderVEMCAL->Draw("SAME");
         for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
         c1->Update();
         sprintf(txtPS,"%s for all SMs",varNameLong[k]);
         ps->TextNDC(0.2,0.94,txtPS);
         hAllSpaceEMCAL[k]->Write();
         //if(k==1)hAllSpaceEMCAL[k]->SaveAs("MeanGaussianFitAllSM.root");
        }
     c1->Update();
    
     for (j=0;j<kNbExtraParams;j++)
        {hAllSpaceEMCALIntg[j]->Write();
         }
     //if (hAllSpaceEMCALIntg[0]->GetMaximum() > maxIntgAmpl) hAllSpaceEMCALIntg[0]->SetMaximum(maxIntgAmpl);
     //Take the previous line out, to be able to catch the noisy towers.
     //if (hAllSpaceEMCALIntg[1]->GetMaximum() > maxIntgAmpl) hAllSpaceEMCALIntg[1]->SetMaximum(maxIntgAmpl);
     //Take the previous line out, because now that histo stores the khi2.
     hAllSpaceEMCALIntg[0]->SetMinimum(0.);
     hAllSpaceEMCALIntg[1]->SetMinimum(0.);
     hAllSpaceEMCALIntg[2]->SetMaximum(1.2*tabMax[8]);
     hAllSpaceEMCALIntg[3]->SetMinimum(0.);
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hAllSpaceEMCALIntg[0]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[0+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(3);
     hAllSpaceEMCALIntg[2]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[2+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(5);
     hAllSpaceEMCALIntg[3]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[3+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     hAllSpaceEMCALIntg[0]->SetMaximum(1.2*tabMax[kNbFitParams+0]);
     c1->cd(2);
     hAllSpaceEMCALIntg[0]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs (zoomed)",varNameLong[2+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
 
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hAllSpaceEMCALIntg[1]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[1+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(3);
     hAllSpaceEMCALIntg[4]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[4+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     hAllSpaceEMCALIntg[1]->SetMaximum(4.0);
     hAllSpaceEMCALIntg[4]->SetMaximum(4.0);
     c1->cd(2);
     hAllSpaceEMCALIntg[1]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs (zoomed)",varNameLong[1+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(4);
     hAllSpaceEMCALIntg[4]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs (zoomed)",varNameLong[4+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(5);
     hAllSpaceEMCALSoverN->Divide(hAllSpaceEMCALIntg[2],hAllSpaceEMCALIntg[3]);
     nbTot=hAllSpaceEMCALSoverN->GetMaximum();
     for (i=0;i<hAllSpaceEMCALSoverN->GetNbinsX();i++)
        {for (j=0;j<hAllSpaceEMCALSoverN->GetNbinsY();j++)
            {if ((hAllSpaceEMCALSoverN->GetBinContent(i+1,j+1) < nbTot) && (hAllSpaceEMCALSoverN->GetBinContent(i+1,j+1) > 0)) nbTot=hAllSpaceEMCALSoverN->GetBinContent(i+1,j+1);
             }
         }
     hAllSpaceEMCALSoverN->SetMinimum(0.98*nbTot);
     hAllSpaceEMCALSoverN->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"Signal to noise ratio for all SMs");
     ps->TextNDC(0.2,0.94,txtPS);
     hAllSpaceEMCALSoverN->Write();
     c1->Update();
     c1->cd(6);
     hAllSpaceEMCALSoverN->SetMaximum(2.);
     hAllSpaceEMCALSoverN->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"Signal to noise ratio for all SMs (zoomed)");
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     }
 
 if (choiceNoDCAL != 0)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,(int)((kNbFitParams+1)/2));
     hAllSpaceDCAL[0]->SetMinimum(0.);
     hAllSpaceDCAL[0]->SetMaximum(1.2*tabMax[0]);
     hAllSpaceDCAL[1]->SetMinimum(tabMin[1]+30.);
     hAllSpaceDCAL[1]->SetMaximum(tabMax[1]-10.);
     hAllSpaceDCAL[2]->SetMinimum(tabMin[2]+5.);
     hAllSpaceDCAL[2]->SetMaximum(tabMax[2]-5.);
     for (j=0;j<kNbFitParams;j++)
        {c1->cd(j+1);
         k=(j%2)*(int)((j+5)/2)+(1-(j%2))*(int)(j/2);
         hAllSpaceDCAL[k]->Draw("COLZ");
         lineSMborderVDCALthird->Draw("SAME");
         lineSMborderVDCAL1->Draw("SAME");
         lineSMborderVDCAL2->Draw("SAME");
         for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
         c1->Update();
         sprintf(txtPS,"%s for all SMs",varNameLong[k]);
         ps->TextNDC(0.2,0.94,txtPS);
         hAllSpaceDCAL[k]->Write();
         //if(k==1)hAllSpaceDCAL[k]->SaveAs("MeanGaussianFitAllSM.root");
        }
     c1->Update();
    
     for (j=0;j<kNbExtraParams;j++)
        {hAllSpaceDCALIntg[j]->Write();
         }
     //if (hAllSpaceDCALIntg[0]->GetMaximum() > maxIntgAmpl) hAllSpaceDCALIntg[0]->SetMaximum(maxIntgAmpl);
     //Take the previous line out, to be able to catch the noisy towers.
     //if (hAllSpaceDCALIntg[1]->GetMaximum() > maxIntgAmpl) hAllSpaceDCALIntg[1]->SetMaximum(maxIntgAmpl);
     //Take the previous line out, because now that histo stores the khi2.
     hAllSpaceDCALIntg[0]->SetMinimum(0.);
     hAllSpaceDCALIntg[1]->SetMinimum(0.);
     hAllSpaceDCALIntg[2]->SetMaximum(1.2*tabMax[8]);
     hAllSpaceDCALIntg[3]->SetMinimum(0.);
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hAllSpaceDCALIntg[0]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[0+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     c1->cd(3);
     hAllSpaceDCALIntg[2]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[2+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     c1->cd(5);
     hAllSpaceDCALIntg[3]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[3+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     hAllSpaceDCALIntg[0]->SetMaximum(1.2*tabMax[kNbFitParams+0]);
     c1->cd(2);
     hAllSpaceDCALIntg[0]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs (zoomed)",varNameLong[0+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
 
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hAllSpaceDCALIntg[1]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[1+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(3);
     hAllSpaceDCALIntg[4]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs",varNameLong[4+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->Update();
     hAllSpaceDCALIntg[1]->SetMaximum(4.0);
     hAllSpaceDCALIntg[4]->SetMaximum(4.0);
     c1->cd(2);
     hAllSpaceDCALIntg[1]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs (zoomed)",varNameLong[1+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(4);
     hAllSpaceDCALIntg[4]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"%s for all SMs (zoomed)",varNameLong[4+kNbFitParams]);
     ps->TextNDC(0.2,0.94,txtPS);
     c1->cd(5);
     hAllSpaceDCALSoverN->Divide(hAllSpaceDCALIntg[2],hAllSpaceDCALIntg[3]);
     nbTot=hAllSpaceDCALSoverN->GetMaximum();
     for (i=0;i<hAllSpaceDCALSoverN->GetNbinsX();i++)
        {for (j=0;j<hAllSpaceDCALSoverN->GetNbinsY();j++)
            {if ((hAllSpaceDCALSoverN->GetBinContent(i+1,j+1) < nbTot) && (hAllSpaceDCALSoverN->GetBinContent(i+1,j+1) > 0)) nbTot=hAllSpaceDCALSoverN->GetBinContent(i+1,j+1);
             }
         }
     hAllSpaceDCALSoverN->SetMinimum(0.98*nbTot);
 /*printf("Min %f max %f\n",hAllSpaceDCALSoverN->GetMinimum(),hAllSpaceDCALSoverN->GetMaximum());
 for (int ccc=0;ccc<2*48;ccc++)
 {for (int rrr=0;rrr<3.33*24;rrr++)
   {printf("Col row %d %d  value %e",ccc,rrr,hAllSpaceDCALSoverN->GetBinContent(ccc+1,rrr+1));
    if (hAllSpaceDCALSoverN->GetBinContent(ccc+1,rrr+1) > 10000.) printf("***********");
    printf("\n");
    }
  }*/
     hAllSpaceDCALSoverN->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"Signal to noise ratio for all SMs");
     ps->TextNDC(0.2,0.94,txtPS);
     hAllSpaceDCALSoverN->Write();
     c1->Update();
     c1->cd(6);
     hAllSpaceDCALSoverN->SetMaximum(2.);
     hAllSpaceDCALSoverN->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     sprintf(txtPS,"Signal to noise ratio for all SMs (zoomed)");
     ps->TextNDC(0.2,0.94,txtPS);
     hAllSpaceDCALSoverN->Write();
     c1->Update();
     }
 

 //Inv mass distrib per SM :
 TF1 *gaussMass1 = new TF1("gaussMass1","gaus(0)",130.,140.);
 TF1 *gaussMass2 = new TF1("gaussMass2","gaus(0)",130.,140.);
 gaussMass1->SetParLimits(0,1.,10000.);
 gaussMass1->SetParLimits(1,125.,145.);
 gaussMass1->SetParLimits(2,0.01,20.);
 gaussMass2->SetParLimits(0,1.,10000.);
 gaussMass2->SetParLimits(1,125.,145.);
 gaussMass2->SetParLimits(2,0.01,20.);
 gaussMass1->SetParameter(0,300.);
 gaussMass1->SetParameter(1,135.);
 gaussMass1->SetParameter(2,1.5);
 gaussMass2->SetParameter(0,300.);
 gaussMass2->SetParameter(1,135.);
 gaussMass2->SetParameter(2,1.5);
 gaussMass1->SetLineColor(kViolet+8);
 gaussMass2->SetLineColor(kPink+8);
 gaussMass1->SetLineWidth(1);
 gaussMass2->SetLineWidth(1);
 linePDGMass->SetLineColor(kTeal-5);
 double tabPeakMeanPerSM[kNbSMtot],tabPeakSigmaPerSM[kNbSMtot],tabPeakKhi2PerSM[kNbSMtot];
 double tabPeakOldMeanPerSM[kNbSMtot],tabPeakOldSigmaPerSM[kNbSMtot],tabPeakOldKhi2PerSM[kNbSMtot];

 if (choiceNoEMCAL != 0)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (i=0;i<kNbSMEMCAL;i++)
        {c1->cd(i+1);
         gPad->SetLogy();
         hAllDistribMassPerSM[i]->Fit(gaussMass1,"BRM 0"); //H option "0" (a fortiori "N" marche aussi), sinon screws up le ps file (pages ou canvas vides inseres).
         tabPeakMeanPerSM[i]=gaussMass1->GetParameter(1);
         tabPeakSigmaPerSM[i]=gaussMass1->GetParameter(2);
         if (gaussMass1->GetNDF() != 0.) tabPeakKhi2PerSM[i]=gaussMass1->GetChisquare()/gaussMass1->GetNDF();
         else tabPeakKhi2PerSM[i]=0.;
         hAllDistribMassPerSM[i]->SetLineColor(4);
         hAllDistribMassPerSM[i]->Draw();
         hAllDistribMassPerSM[i]->Write();
         gaussMass1->Draw("SAME");
         linePDGMass->SetY2(1.05*hAllDistribMassPerSM[i]->GetMaximum());
         linePDGMass->Draw("SAME");
         hAllOldDistribMassPerSM[i]->Fit(gaussMass2,"BRM 0"); //H option "0" (a fortiori "N" marche aussi), sinon screws up le ps file (pages ou canvas vides inseres).
         tabPeakOldMeanPerSM[i]=gaussMass2->GetParameter(1);
         tabPeakOldSigmaPerSM[i]=gaussMass2->GetParameter(2);
         if (gaussMass2->GetNDF() != 0.) tabPeakOldKhi2PerSM[i]=gaussMass2->GetChisquare()/gaussMass2->GetNDF();
         else tabPeakOldKhi2PerSM[i]=0.;
         hAllOldDistribMassPerSM[i]->SetLineColor(2);
         hAllOldDistribMassPerSM[i]->Draw("SAME");
         hAllOldDistribMassPerSM[i]->Write();
         gaussMass2->Draw("SAME");
         linePDGMass->SetY2(1.05*hAllOldDistribMassPerSM[i]->GetMaximum());
         linePDGMass->Draw("SAME");
         c1->RedrawAxis();
         c1->Update();
         }
     c1->Update();
     }
 
 if (choiceNoDCAL != 0)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (i=kNbSMEMCAL;i<kNbSMtot;i++)
        {c1->cd(i+1-kNbSMEMCAL);
         gPad->SetLogy();
         hAllDistribMassPerSM[i]->Fit(gaussMass1,"BRM 0"); //H option "0" (a fortiori "N" marche aussi), sinon screws up le ps file (pages ou canvas vides inseres).
         tabPeakMeanPerSM[i]=gaussMass1->GetParameter(1);
         tabPeakSigmaPerSM[i]=gaussMass1->GetParameter(2);
         if (gaussMass1->GetNDF() != 0.) tabPeakKhi2PerSM[i]=gaussMass1->GetChisquare()/gaussMass1->GetNDF();
         else tabPeakKhi2PerSM[i]=0.;
         hAllDistribMassPerSM[i]->SetLineColor(4);
         hAllDistribMassPerSM[i]->Draw();
         hAllDistribMassPerSM[i]->Write();
         gaussMass1->Draw("SAME");
         linePDGMass->SetY2(1.05*hAllDistribMassPerSM[i]->GetMaximum());
         linePDGMass->Draw("SAME");
         hAllOldDistribMassPerSM[i]->Fit(gaussMass2,"BRM 0"); //H option "0" (a fortiori "N" marche aussi), sinon screws up le ps file (pages ou canvas vides inseres).
         tabPeakOldMeanPerSM[i]=gaussMass2->GetParameter(1);
         tabPeakOldSigmaPerSM[i]=gaussMass2->GetParameter(2);
         if (gaussMass2->GetNDF() != 0.) tabPeakOldKhi2PerSM[i]=gaussMass2->GetChisquare()/gaussMass2->GetNDF();
         else tabPeakOldKhi2PerSM[i]=0.;
         hAllOldDistribMassPerSM[i]->SetLineColor(2);
         hAllOldDistribMassPerSM[i]->Draw("SAME");
         hAllOldDistribMassPerSM[i]->Write();
         gaussMass2->Draw("SAME");
         linePDGMass->SetY2(1.05*hAllOldDistribMassPerSM[i]->GetMaximum());
         linePDGMass->Draw("SAME");
         c1->RedrawAxis();
         c1->Update();
         }
     c1->Update();
     }
 
 printf("\n");
 for (i=0;i<kNbSMtot;i++)
    {printf("Mass distr in SM %d : new : %6.2f sig %4.2f MeV, khi2 = %4.2f       vs old : %6.2f sig %4.2f MeV, khi2 = %4.2f\n",i,tabPeakMeanPerSM[i],tabPeakSigmaPerSM[i],tabPeakKhi2PerSM[i],tabPeakOldMeanPerSM[i],tabPeakOldSigmaPerSM[i],tabPeakOldKhi2PerSM[i]);
     hMassPerSM->SetBinContent(i+1,tabPeakMeanPerSM[i]);
     hMassPerSM->SetBinError(i+1,tabPeakOldSigmaPerSM[i]);
     hMassOldPerSM->SetBinContent(i+1,tabPeakOldMeanPerSM[i]);
     hMassOldPerSM->SetBinError(i+1,tabPeakOldSigmaPerSM[i]);
     }
 printf("-----------------------------------------------------\n");

 //Correction factors distribs :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,4);
 c1->cd(1);
 hCoeff->SetLineColor(4);
 hCoeffOldCorr->SetLineColor(2);
 hCoeffOld->SetLineColor(kOrange+7);
 gPad->SetLogy();
 hCoeff->Draw();
 if (isFirstIteration == 0)
    {hCoeffOld->Draw("SAME");
     hCoeffOldCorr->Draw("SAME");
     }
 hCoeff->Write();
 hCoeffOld->Write();
 hCoeffOldCorr->Write();
 if (isFirstIteration == 0)
    {c1->cd(3);
     hCoeffOldCorrLarge->SetLineColor(2);
     hCoeffOldLarge->SetLineColor(kOrange+7);
     gPad->SetLogy();
     hCoeffOldLarge->Draw();
     hCoeffOldCorrLarge->Draw("SAME");
     hCoeffOldLarge->Write();
     hCoeffOldCorrLarge->Write();
     }
 c1->cd(4);
 gaussMass1->SetParameter(0,300.);
 gaussMass1->SetParameter(1,135.);
 gaussMass1->SetParameter(2,1.5);
 gaussMass2->SetParameter(0,300.);
 gaussMass2->SetParameter(1,135.);
 gaussMass2->SetParameter(2,1.5);
 hAllDistribMass->Fit(gaussMass1,"BRM 0"); //H option "0" (a fortiori "N" marche aussi), sinon screws up le ps file (pages ou canvas vides inseres).
 if (isFirstIteration == 0) hAllOldDistribMass->Fit(gaussMass2,"BRM 0"); //Et H option "I", sinon converge vers n'importe-quoi (avec "L" aussi, et dans tous les cas en settant a 10000 l'error des empty bins aussi).
 //Si on met option I, cette version de root breake (probleme de librairie). Mais now ca marche sans l'option I (as expected), => OK.
 printf("\nTower mass distribution : new : %f sig %f MeV, khi2 = %f\n",gaussMass1->GetParameter(1),gaussMass1->GetParameter(2),gaussMass1->GetChisquare()/gaussMass1->GetNDF());
 if (isFirstIteration == 0) printf("                          old : %f sig %f MeV, khi2 = %f\n\n",gaussMass2->GetParameter(1),gaussMass2->GetParameter(2),gaussMass2->GetChisquare()/gaussMass2->GetNDF());
 hAllDistribMass->SetLineColor(4);
 hAllDistribMass->Draw();
 hAllDistribMass->Write();
 gaussMass1->Draw("SAME");
 linePDGMass->SetLineColor(kTeal-5);
 linePDGMass->SetY2(1.05*hAllDistribMass->GetMaximum());
 linePDGMass->Draw("SAME");
 c1->Update();
 if (isFirstIteration == 0)
    {c1->cd(2);
     hAllOldDistribMass->SetLineColor(2);
     hAllOldDistribMass->Draw();
     hAllOldDistribMass->Write();
     gaussMass2->Draw("SAME");
     linePDGMass->SetY2(1.05*hAllOldDistribMass->GetMaximum());
     linePDGMass->Draw("SAME");
     c1->Update();
     }
 c1->cd(6);
 hMassPerSM->SetLineColor(4);
 hMassOldPerSM->SetLineColor(2);
 hMassPerSM->SetMaximum(mPDG+2.5);
 hMassPerSM->SetMinimum(mPDG-2.5);
 hMassPerSM->Draw();
 hMassPerSM->Write();
 if (isFirstIteration == 0) hMassOldPerSM->Draw("SAME");
 hMassOldPerSM->Write();
 lineMassH->SetX1(hMassPerSM->GetXaxis()->GetXmin());
 lineMassH->SetX2(hMassPerSM->GetXaxis()->GetXmax());
 lineMassH->SetLineColor(kTeal-5);
 lineMassH->Draw("SAME");
 c1->Update();
 if (choiceNoEMCAL != 0)
    {c1->cd(7);
     nbTot=hAllSpaceEMCALCoeff[0]->GetMaximum();
     for (i=0;i<hAllSpaceEMCALCoeff[0]->GetNbinsX();i++)
        {for (j=0;j<hAllSpaceEMCALCoeff[0]->GetNbinsY();j++)
            {if ((hAllSpaceEMCALCoeff[0]->GetBinContent(i+1,j+1) < nbTot) && (hAllSpaceEMCALCoeff[0]->GetBinContent(i+1,j+1) != 0))      nbTot=hAllSpaceEMCALCoeff[0]->GetBinContent(i+1,j+1);
             }
         }
     hAllSpaceEMCALCoeff[0]->SetMinimum(0.98*nbTot);
     hAllSpaceEMCALCoeff[0]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     }
 hAllSpaceEMCALCoeff[0]->Write();
 if (choiceNoDCAL != 0)
    {c1->cd(8);
     nbTot=hAllSpaceDCALCoeff[0]->GetMaximum();
     for (i=0;i<hAllSpaceDCALCoeff[0]->GetNbinsX();i++)
        {for (j=0;j<hAllSpaceDCALCoeff[0]->GetNbinsY();j++)
            {if ((hAllSpaceDCALCoeff[0]->GetBinContent(i+1,j+1) < nbTot) && (hAllSpaceDCALCoeff[0]->GetBinContent(i+1,j+1) != 0))      nbTot=hAllSpaceDCALCoeff[0]->GetBinContent(i+1,j+1);
             }
         }
     hAllSpaceDCALCoeff[0]->SetMinimum(0.98*nbTot);
     hAllSpaceDCALCoeff[0]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     }
 hAllSpaceDCALCoeff[0]->Write();
 c1->Update();
 printf("-----------------------------------------------------\n");


 //Add also distribs of distance to mPDG in err bars
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(1);
 hAllSpaceEMCALDistance[0]->SetMinimum(0.);
 hAllSpaceEMCALDistance[0]->SetMaximum(10.);
 hAllSpaceEMCALDistance[0]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALDistance[0]->Write();
 c1->cd(3);
 hAllSpaceEMCALDistance[1]->SetMinimum(0.);
 hAllSpaceEMCALDistance[1]->SetMaximum(10.);
 hAllSpaceEMCALDistance[1]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALDistance[1]->Write();
 c1->cd(5);
 hAllSpaceEMCALDistance[2]->SetMinimum(-1.);
 hAllSpaceEMCALDistance[2]->SetMaximum(4.);
 hAllSpaceEMCALDistance[2]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALDistance[2]->Write();
 c1->cd(2);
 hAllSpaceDCALDistance[0]->SetMinimum(0.);
 hAllSpaceDCALDistance[0]->SetMaximum(10.);
 hAllSpaceDCALDistance[0]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALDistance[0]->Write();
 c1->cd(4);
 hAllSpaceDCALDistance[1]->SetMinimum(0.);
 hAllSpaceDCALDistance[1]->SetMaximum(10.);
 hAllSpaceDCALDistance[1]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALDistance[1]->Write();
 c1->cd(6);
 hAllSpaceDCALDistance[2]->SetMinimum(-1.);
 hAllSpaceDCALDistance[2]->SetMaximum(4.);
 hAllSpaceDCALDistance[2]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALDistance[2]->Write();
 c1->Update();
 
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(1);
 TF1 *fGauss = new TF1("fGauss","gaus(0)",0.,4.);
 TF1 *fGauss2 = new TF1("fGauss2","gaus(0)",0.,4.);
 fGauss->SetParameter(0,100.);
 fGauss->SetParameter(1,0.);
 fGauss->SetParameter(2,1.);
 fGauss->SetLineWidth(1);
 fGauss->SetLineColor(kOrange-6);
 fGauss2->SetParameter(0,100.);
 fGauss2->FixParameter(1,0.);
 fGauss2->SetParameter(2,1.);
 fGauss2->SetLineWidth(1);
 fGauss2->SetLineColor(kTeal-7);
 hAllDistribDistance[0]->SetLineColor(kRed);
 hAllDistribDistance[1]->SetLineColor(kGreen+2);
 hAllDistribDistance[2]->SetLineColor(kBlue);
 hAllDistribDistance[1]->Fit(fGauss,"BR");
 hAllDistribDistance[1]->Fit(fGauss2,"BR");
 printf("\nDistrib tower normalized distances : mu = %5.3f, sig = %5.3f, khi2/ndf = %f\n",fGauss->GetParameter(1),fGauss->GetParameter(2),fGauss->GetChisquare()/(double)TMath::Max((double)fGauss->GetNDF(),0.0000000001));
 printf("                       with 0 mean : mu = %5.3f, sig = %5.3f, khi2/ndf = %f\n\n",fGauss2->GetParameter(1),fGauss2->GetParameter(2),fGauss2->GetChisquare()/(double)TMath::Max((double)fGauss2->GetNDF(),0.0000000001));
 hAllDistribDistance[2]->Draw();
 hAllDistribDistance[1]->Draw("SAME");
 fGauss->Draw("SAME");
 fGauss2->Draw("SAME");
 hAllDistribDistance[0]->Draw("SAME");
 hAllDistribDistance[0]->Write();
 hAllDistribDistance[1]->Write();
 hAllDistribDistance[2]->Write();
 c1->cd(3);
 hCorrelDvsSig->Draw("COLZ");
 hCorrelDvsSig->Write();
 c1->cd(5);
 hCorrelCmpD->Draw("COLZ");
 hCorrelCmpD->Write();
 c1->cd(2);
 gPad->SetLogy();
 hAllDistribDistance2[2*0+0]->Draw();
 hAllDistribDistance2[2*1+0]->SetLineColor(kGreen+2);
 hAllDistribDistance2[2*1+0]->Draw("SAME");
 lineMax->SetX1(100.*uncertKSpec*uncertSpec);
 lineMax->SetX2(100.*uncertKSpec*uncertSpec);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribDistance2[2*0+0]->GetMaximum());
 lineMax->Draw("SAME");
 hAllDistribDistance2[2*0+0]->Write();
 hAllDistribDistance2[2*1+0]->Write();
 c1->Update();
 c1->cd(4);
 gPad->SetLogy();
 TF1 *fExpectedUncertDistrib = new TF1("fExpectedUncertDistrib","gaus",0.,hAllDistribDistance2[2*0+1]->GetXaxis()->GetXmax());
 fExpectedUncertDistrib->SetParameter(0,(2.*nbTowersConsideredTot*hAllDistribDistance2[2*0+1]->GetBinWidth(1))/(uncertSpec*TMath::Sqrt(2.*TMath::Pi()))); //Factor of 2 because of the absolute value.
 fExpectedUncertDistrib->SetParameter(1,0.);
 fExpectedUncertDistrib->SetParameter(2,uncertSpec);
 hAllDistribDistance2[2*0+1]->Draw();
 hAllDistribDistance2[2*1+1]->SetLineColor(kGreen+2);
 hAllDistribDistance2[2*1+1]->Draw("SAME");
 fExpectedUncertDistrib->Draw("SAME");
 lineMax->SetX1(uncertDist*uncertSpec);
 lineMax->SetX2(uncertDist*uncertSpec);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribDistance2[2*0+1]->GetMaximum());
 lineMax->Draw("SAME");
 hAllDistribDistance2[2*0+1]->Write();
 hAllDistribDistance2[2*1+1]->Write();
 fExpectedUncertDistrib->Write();
 c1->Update();


 //Control of convergence
 if (isFirstIteration == 0)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hCtrlEvolMu->Draw();
     hCtrlEvolMu->Write();
     c1->cd(3);
     h2CtrlEvolMu1->Draw("COLZ");
     h2CtrlEvolMu1->Write();
     c1->cd(5);
     TLine *bissec = new TLine(h2CtrlEvolMu2->GetYaxis()->GetXmin(),h2CtrlEvolMu2->GetYaxis()->GetXmin(),h2CtrlEvolMu2->GetYaxis()->GetXmax(),h2CtrlEvolMu2->GetYaxis()->GetXmax());
     TLine *bi2ssec = new TLine(h2CtrlEvolMu2->GetYaxis()->GetXmin()/2.,h2CtrlEvolMu2->GetYaxis()->GetXmin(),h2CtrlEvolMu2->GetYaxis()->GetXmax()/2.,h2CtrlEvolMu2->GetYaxis()->GetXmax());
     TLine *axeX = new TLine(h2CtrlEvolMu2->GetXaxis()->GetXmin(),0.,h2CtrlEvolMu2->GetXaxis()->GetXmax(),0.);
     TLine *axeY = new TLine(0.,h2CtrlEvolMu2->GetYaxis()->GetXmin(),0.,h2CtrlEvolMu2->GetYaxis()->GetXmax());
     h2CtrlEvolMu2->Draw("COLZ");
     bissec->Draw("SAME");
     bi2ssec->Draw("SAME");
     axeX->Draw("SAME");
     axeY->Draw("SAME");
     h2CtrlEvolMu2->Write();
     c1->cd(2);
     hCtrlEvolCoeff->Draw();
     hCtrlEvolCoeff->Write();
     c1->cd(4);
     h2CtrlEvolCoeff1->Draw("COLZ");
     h2CtrlEvolCoeff1->Write();
     c1->cd(6);
     TLine *bissec2 = new TLine(h2CtrlEvolCoeff2->GetYaxis()->GetXmin(),h2CtrlEvolCoeff2->GetYaxis()->GetXmin(),h2CtrlEvolCoeff2->GetYaxis()->GetXmax(),h2CtrlEvolCoeff2->GetYaxis()->GetXmax());
     TLine *bi2ssec2 = new TLine(h2CtrlEvolCoeff2->GetYaxis()->GetXmin()/2.,h2CtrlEvolCoeff2->GetYaxis()->GetXmin(),h2CtrlEvolCoeff2->GetYaxis()->GetXmax()/2.,h2CtrlEvolCoeff2->GetYaxis()->GetXmax());
     TLine *axeX2 = new TLine(h2CtrlEvolCoeff2->GetXaxis()->GetXmin(),0.,h2CtrlEvolCoeff2->GetXaxis()->GetXmax(),0.);
     TLine *axeY2 = new TLine(0.,h2CtrlEvolCoeff2->GetYaxis()->GetXmin(),0.,h2CtrlEvolCoeff2->GetYaxis()->GetXmax());
     h2CtrlEvolCoeff2->Draw("COLZ");
     bissec2->Draw("SAME");
     bi2ssec2->Draw("SAME");
     axeX2->Draw("SAME");
     axeY2->Draw("SAME");
     h2CtrlEvolCoeff2->Write();
     c1->Update();
     }
 

 //1-dimensional distribs :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,(int)((kNbFitParams+1)/2));
 for (j=0;j<kNbFitParams;j++)
    {c1->cd(j+1);
     k=(j%2)*(int)((j+5)/2)+(1-(j%2))*(int)(j/2);
     gPad->SetLogy();
     hAllDistrib[k]->SetLineColor(4);
     hAllOldDistrib[k]->SetLineColor(2);
     hAllDistrib[k]->Draw();
     hAllOldDistrib[k]->Draw("SAME");
     lineMin->SetX1(cutMin[k]);
     lineMin->SetX2(cutMin[k]);
     lineMin->SetY1(0.);
     lineMin->SetY2(hAllDistrib[k]->GetMaximum());
     lineMax->SetX1(cutMax[k]);
     lineMax->SetX2(cutMax[k]);
     lineMax->SetY1(0.);
     lineMax->SetY2(hAllDistrib[k]->GetMaximum());
     lineMin->Draw("SAME");
     lineMax->Draw("SAME");
     if (k == 1)
        {linePDGMass->SetY2(1.5*hAllDistrib[k]->GetMaximum());
         linePDGMass->Draw("SAME");
         }
     c1->RedrawAxis();
     c1->Update(); //Car sinon, update les parametres du TLine de tous les plots avec ceux du dernier plot trac\'e.
     hAllDistrib[k]->Write();
     hAllOldDistrib[k]->Write();
     }
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 for (j=0;j<kNbExtraParams;j++)
    {hAllDistribIntg[j]->SetLineColor(4);
     hAllDistribIntg[j]->Write();
     if (j < kNbExtraParamsToBeRead)
        {hAllOldDistrib[j+kNbFitParams]->SetLineColor(2);
         hAllOldDistrib[j+kNbFitParams]->Write();
         }
     }
 c1->cd(1);
 gPad->SetLogy();
 hAllDistribIntg[0]->Draw();
 hAllOldDistrib[0+kNbFitParams]->Draw("SAME");
 lineMin->SetX1(cutMin[0+kNbFitParams]);
 lineMin->SetX2(cutMin[0+kNbFitParams]);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribIntg[0]->GetMaximum());
 lineMax->SetX1(cutMax[0+kNbFitParams]);
 lineMax->SetX2(cutMax[0+kNbFitParams]);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribIntg[0]->GetMaximum());
 lineMin->Draw("SAME");
 lineMax->Draw("SAME");
 c1->RedrawAxis();
 c1->Update(); //Car sinon, update les parametres du TLine de tous les plots avec ceux du dernier plot trac\'e.
 c1->cd(2);
 //That histo is not filled with the towers with 0 entries.
 hAllDistribNbEntries->SetLineColor(4);
 hAllDistribNbEntries->Draw();
 lineMin->SetX1(cutEntriesDiscard);
 lineMin->SetX2(cutEntriesDiscard);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribNbEntries->GetMaximum());
 lineMin->Draw("SAME");
 hAllDistribNbEntries->Write();
 hAllDistribNbEntriesCorrel->Write();
 c1->RedrawAxis();
 c1->Update();
 c1->cd(3);
 gPad->SetLogy();
 hAllDistribIntg[2]->Draw();
 hAllOldDistrib[2+kNbFitParams]->Draw("SAME");
 lineMin->SetX1(cutMin[2+kNbFitParams]);
 lineMin->SetX2(cutMin[2+kNbFitParams]);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribIntg[2]->GetMaximum());
 lineMax->SetX1(cutMax[2+kNbFitParams]);
 lineMax->SetX2(cutMax[2+kNbFitParams]);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribIntg[2]->GetMaximum());
 lineMin->Draw("SAME");
 lineMax->Draw("SAME");
 c1->RedrawAxis();
 c1->Update(); //Car sinon, update les parametres du TLine de tous les plots avec ceux du dernier plot trac\'e.
 c1->cd(4);
 gPad->SetLogy();
 hAllDistribIntg[3]->Draw();
 hAllOldDistrib[3+kNbFitParams]->Draw("SAME");
 lineMin->SetX1(cutMin[3+kNbFitParams]);
 lineMin->SetX2(cutMin[3+kNbFitParams]);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribIntg[3]->GetMaximum());
 lineMax->SetX1(cutMax[3+kNbFitParams]);
 lineMax->SetX2(cutMax[3+kNbFitParams]);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribIntg[3]->GetMaximum());
 lineMin->Draw("SAME");
 lineMax->Draw("SAME");
 c1->RedrawAxis();
 c1->Update(); //Car sinon, update les parametres du TLine de tous les plots avec ceux du dernier plot trac\'e.

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(1);
 gPad->SetLogy();
 hAllDistribIntg[1]->Draw();
 hAllOldDistrib[1+kNbFitParams]->Draw("SAME");
 lineMin->SetX1(cutMin[1+kNbFitParams]);
 lineMin->SetX2(cutMin[1+kNbFitParams]);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribIntg[1]->GetMaximum());
 lineMax->SetX1(cutMax[1+kNbFitParams]);
 lineMax->SetX2(cutMax[1+kNbFitParams]);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribIntg[1]->GetMaximum());
 lineMin->Draw("SAME");
 lineMax->Draw("SAME");
 c1->RedrawAxis();
 c1->Update(); //Car sinon, update les parametres du TLine de tous les plots avec ceux du dernier plot trac\'e.
 c1->cd(2);
 gPad->SetLogy();
 hAllDistribIntg[4]->Draw();
 hAllOldDistrib[4+kNbFitParams]->Draw("SAME");
 lineMin->SetX1(cutMin[4+kNbFitParams]);
 lineMin->SetX2(cutMin[4+kNbFitParams]);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribIntg[4]->GetMaximum());
 lineMax->SetX1(cutMax[4+kNbFitParams]);
 lineMax->SetX2(cutMax[4+kNbFitParams]);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribIntg[4]->GetMaximum());
 lineMin->Draw("SAME");
 lineMax->Draw("SAME");
 c1->RedrawAxis();
 c1->Update(); //Car sinon, update les parametres du TLine de tous les plots avec ceux du dernier plot trac\'e.
 c1->cd(3);
 gPad->SetLogy();
 hAllDistribIntg[5]->Draw();
 //hAllOldDistrib[5+kNbFitParams]->Draw("SAME");
 lineMin->SetX1(cutMin[5+kNbFitParams]);
 lineMin->SetX2(cutMin[5+kNbFitParams]);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribIntg[5]->GetMaximum());
 lineMax->SetX1(cutMax[5+kNbFitParams]);
 lineMax->SetX2(cutMax[5+kNbFitParams]);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribIntg[5]->GetMaximum());
 lineMin->Draw("SAME");
 lineMax->Draw("SAME");
 c1->RedrawAxis();
 c1->Update(); //Car sinon, update les parametres du TLine de tous les plots avec ceux du dernier plot trac\'e.
 c1->cd(4);
 gPad->SetLogy();
 hAllDistribIntg[6]->Draw();
 //hAllOldDistrib[6+kNbFitParams]->Draw("SAME");
 lineMin->SetX1(cutMin[6+kNbFitParams]);
 lineMin->SetX2(cutMin[6+kNbFitParams]);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribIntg[6]->GetMaximum());
 lineMax->SetX1(cutMax[6+kNbFitParams]);
 lineMax->SetX2(cutMax[6+kNbFitParams]);
 lineMax->SetY1(0.);
 lineMax->SetY2(hAllDistribIntg[6]->GetMaximum());
 lineMin->Draw("SAME");
 lineMax->Draw("SAME");
 c1->RedrawAxis();
 c1->Update(); //Car sinon, update les parametres du TLine de tous les plots avec ceux du dernier plot trac\'e.
 c1->cd(5);
 gPad->SetLogy();
 hAllDistribIntg[5]->SetLineColor(kViolet-5);
 hAllDistribIntg[6]->SetLineColor(kAzure+9);
 hAllDistribIntg[5]->Rebin(2);
 hAllDistribIntg[6]->Rebin(2);
 hAllDistribIntg[4]->Scale(0.5);
 hAllDistribIntg[6]->Draw();
 hAllDistribIntg[5]->Draw("SAME");
 hAllDistribIntg[4]->Draw("SAME");
 lineMin->SetX1(1.0);
 lineMin->SetX2(1.0);
 lineMin->SetY1(0.);
 lineMin->SetY2(hAllDistribIntg[6]->GetMaximum());
 lineMin->Draw("SAME");
 c1->RedrawAxis();
 c1->Update();

 //Masked/empty towers, and towers of which calib coeff has been automatically put to 1.0 for lack of stat :
 if (choiceNoEMCAL != 0)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (i=0;i<kNbSMEMCAL;i++)
        {c1->cd(i+1);
         //Red = masked/empty tower ; blue = lack of stat.
         hSpaceEntriesDiscard[1*i+0]->SetMinimum(0.);
         hSpaceEntriesDiscard[1*i+0]->SetMaximum(5.05);
         hSpaceEntriesDiscard[1*i+0]->Draw("COLZ");
         hSpaceEntriesDiscard[1*i+0]->Write();
         c1->RedrawAxis();
         }
     c1->Update();
     }
 
 if (choiceNoDCAL != 0)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (i=kNbSMEMCAL;i<kNbSMtot;i++)
        {c1->cd(i+1-kNbSMEMCAL);
         //Red = masked/empty tower ; blue = lack of stat.
         hSpaceEntriesDiscard[1*i+0]->SetMinimum(0.);
         hSpaceEntriesDiscard[1*i+0]->SetMaximum(5.05);
         hSpaceEntriesDiscard[1*i+0]->Draw("COLZ");
         hSpaceEntriesDiscard[1*i+0]->Write();
         c1->RedrawAxis();
         }
     c1->Update();
     }
 

 //Correlations by group of 3 :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(1);
 hCorrelMuVsA->Draw("COLZ");
 c1->cd(3);
 hCorrelSigVsA->Draw("COLZ");
 c1->cd(5);
 hCorrelSigVsMu->Draw("COLZ");
 c1->cd(2);
 hCorrelBVsA->Draw("COLZ");
 c1->cd(4);
 hCorrelCVsA->Draw("COLZ");
 c1->cd(6);
 hCorrelCVsB->Draw("COLZ");
 hCorrelMuVsA->Write();
 hCorrelSigVsA->Write();
 hCorrelSigVsMu->Write();
 hCorrelBVsA->Write();
 hCorrelCVsA->Write();
 hCorrelCVsB->Write();
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,4);
 c1->cd(1);
 hCorrelISVsI->Draw("COLZ");
 c1->cd(3);
 hCorrelSVsI->Draw("COLZ");
 c1->cd(5);
 hCorrelSVsIS->Draw("COLZ");
 c1->cd(7);
 hCorrelSbincountingVsS->Draw("COLZ");
 TLine *bissecCorrel = new TLine(hCorrelSbincountingVsS->GetYaxis()->GetXmin(),hCorrelSbincountingVsS->GetYaxis()->GetXmin(),hCorrelSbincountingVsS->GetYaxis()->GetXmax(),hCorrelSbincountingVsS->GetYaxis()->GetXmax());
 bissecCorrel->Draw("SAME");
 c1->cd(2);
 hCorrelMuVsI->Draw("COLZ");
 c1->cd(4);
 hCorrelMuVsS->Draw("COLZ");
 c1->cd(6);
 hCorrelSigVsI->Draw("COLZ");
 c1->cd(8);
 hCorrelSigVsS->Draw("COLZ");
 hCorrelISVsI->Write();
 hCorrelSVsI->Write();
 hCorrelSVsIS->Write();
 hCorrelSbincountingVsS->Write();
 hCorrelMuVsI->Write();
 hCorrelMuVsS->Write();
 hCorrelSigVsI->Write();
 hCorrelSigVsS->Write();
 c1->Update();
 

 //Changes wrt previous pass (space-like) :
 if (isFirstIteration == 0)
    {if (choiceNoEMCAL != 0)
     {ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     for (j=0;j<6;j++)
        {c1->cd(j+1);
         k=(j%2)*(int)((j+5)/2)+(1-(j%2))*(int)(j/2);
         nbTot=hAllSpaceEMCALDiff[k]->GetMaximum();
         for (i=0;i<hAllSpaceEMCALDiff[k]->GetNbinsX();i++)
            {for (int jm=0;jm<hAllSpaceEMCALDiff[k]->GetNbinsY();jm++)
                {if ((hAllSpaceEMCALDiff[k]->GetBinContent(i+1,jm+1) < nbTot) && (hAllSpaceEMCALDiff[k]->GetBinContent(i+1,jm+1) != 0)) nbTot=hAllSpaceEMCALDiff[k]->GetBinContent(i+1,jm+1);
                 }
             }
         if ((j%2) == 0) hAllSpaceEMCALDiff[k]->SetMinimum(0.98*nbTot);
         hAllSpaceEMCALDiff[k]->Draw("COLZ");
         lineSMborderVEMCAL->Draw("SAME");
         for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
         hAllSpaceEMCALDiff[k]->Write();
         }
     c1->Update();
    
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     for (j=0;j<4;j++)
        {c1->cd(j+1);
         if (j == 3) c1->cd(j+2);
         nbTot=hAllSpaceEMCALDiff[j+6]->GetMaximum();
         for (i=0;i<hAllSpaceEMCALDiff[j+6]->GetNbinsX();i++)
            {for (int jm=0;jm<hAllSpaceEMCALDiff[j+6]->GetNbinsY();jm++)
                {if ((hAllSpaceEMCALDiff[j+6]->GetBinContent(i+1,jm+1) < nbTot) && (hAllSpaceEMCALDiff[j+6]->GetBinContent(i+1,jm+1) != 0)) nbTot=hAllSpaceEMCALDiff[j+6]->GetBinContent(i+1,jm+1);
                 }
             }
         hAllSpaceEMCALDiff[j+6]->SetMinimum(0.98*nbTot);
         hAllSpaceEMCALDiff[j+6]->Draw("COLZ");
         lineSMborderVEMCAL->Draw("SAME");
         for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
         hAllSpaceEMCALDiff[j+6]->Write();
         }
     c1->Update();
     }
     
     if (choiceNoDCAL != 0)
     {ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     for (j=0;j<6;j++)
        {c1->cd(j+1);
         k=(j%2)*(int)((j+5)/2)+(1-(j%2))*(int)(j/2);
         nbTot=hAllSpaceDCALDiff[k]->GetMaximum();
         for (i=0;i<hAllSpaceDCALDiff[k]->GetNbinsX();i++)
            {for (int jm=0;jm<hAllSpaceDCALDiff[k]->GetNbinsY();jm++)
                {if ((hAllSpaceDCALDiff[k]->GetBinContent(i+1,jm+1) < nbTot) && (hAllSpaceDCALDiff[k]->GetBinContent(i+1,jm+1) != 0)) nbTot=hAllSpaceDCALDiff[k]->GetBinContent(i+1,jm+1);
                 }
             }
         if ((j%2) == 0) hAllSpaceDCALDiff[k]->SetMinimum(0.98*nbTot);
         hAllSpaceDCALDiff[k]->Draw("COLZ");
         lineSMborderVDCALthird->Draw("SAME");
         lineSMborderVDCAL1->Draw("SAME");
         lineSMborderVDCAL2->Draw("SAME");
         for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
         hAllSpaceDCALDiff[k]->Write();
         }
     c1->Update();
    
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     for (j=0;j<4;j++)
        {c1->cd(j+1);
         if (j == 3) c1->cd(j+2);
         nbTot=hAllSpaceDCALDiff[j+6]->GetMaximum();
         for (i=0;i<hAllSpaceDCALDiff[j+6]->GetNbinsX();i++)
            {for (int jm=0;jm<hAllSpaceDCALDiff[j+6]->GetNbinsY();jm++)
                {if ((hAllSpaceDCALDiff[j+6]->GetBinContent(i+1,jm+1) < nbTot) && (hAllSpaceDCALDiff[j+6]->GetBinContent(i+1,jm+1) != 0)) nbTot=hAllSpaceDCALDiff[j+6]->GetBinContent(i+1,jm+1);
                 }
             }
         hAllSpaceDCALDiff[j+6]->SetMinimum(0.98*nbTot);
         hAllSpaceDCALDiff[j+6]->Draw("COLZ");
         lineSMborderVDCALthird->Draw("SAME");
         lineSMborderVDCAL1->Draw("SAME");
         lineSMborderVDCAL2->Draw("SAME");
         for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
         hAllSpaceDCALDiff[j+6]->Write();
         }
     c1->Update();
     }
    
     //Changes wrt previous pass (1-D) :
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     for (j=0;j<6;j++)
        {c1->cd(j+1);
         k=(j%2)*(int)((j+5)/2)+(1-(j%2))*(int)(j/2);
         gPad->SetLogy();
         hAllDiff[k]->SetLineColor(kAzure-8);
         hAllDiff[k]->SetFillColor(kAzure-9);
         hAllDiffAllTw[k]->SetLineColor(kBlue+1);
         hAllDiffAllTw[k]->Draw();
         hAllDiff[k]->Draw("SAME");
         hAllDiffAllTw[k]->Draw("SAME");
         hAllDiff[k]->Write();
         hAllDiffAllTw[k]->Write();
         c1->RedrawAxis();
         }
     c1->Update();
    
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     for (j=0;j<4;j++)
        {c1->cd(j+1);
         if (j == 3) c1->cd(j+2);
         gPad->SetLogy();
         hAllDiff[j+6]->SetLineColor(kAzure-8);
         hAllDiff[j+6]->SetFillColor(kAzure-9);
         hAllDiffAllTw[j+6]->SetLineColor(kBlue+1);
         hAllDiffAllTw[j+6]->Draw();
         hAllDiff[j+6]->Draw("SAME");
         hAllDiffAllTw[j+6]->Draw("SAME");
         hAllDiff[j+6]->Write();
         hAllDiffAllTw[j+6]->Write();
         c1->RedrawAxis();
         }
     c1->Update();
     }
 

 ps->Close();
 
 txtFileCalibOut.close();
 fclose(txtFileParamsOut);
 fclose(txtFileParamsIn);
 fclose(txtFileCalibIn);
 rootFileOut->Close();

 printf("Fini...\n");


 return;
 }








//-----------------------------------------------------------------------------
Double_t pi0massP3(Double_t *x, Double_t *par)
{
 if (par[2] == 0.) printf("Unvalid (=zero) gaussian width.\n");
 Double_t gaus;
 if (par[2] != 0.) gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
    else gaus = 99999999.; 
 Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
 return gaus+back;
}










//-----------------------------------------------------------------------------
Double_t pi0massP2(Double_t *x, Double_t *par)
{
 if (par[2] == 0.) printf("Unvalid (=zero) gaussian width.\n");
 Double_t gaus;
 if (par[2] != 0.) gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
    else gaus = 99999999.; 
 Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0];
 return gaus+back;
}










//-----------------------------------------------------------------------------
Double_t pi0massP1(Double_t *x, Double_t *par)
{
 Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
 Double_t back = par[3] + par[4]*x[0];
 return gaus+back;
}












