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
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TFile.h"
#include "TMath.h"
#include "TLeaf.h" 
#include "TBranch.h" 

//#include "/cebaf/cebaf/EMCAL/cosmicsAnalysis/macros/defineMyPalette40All.C"

///
/// \file MultiplyPi0CalibrationFactors_TextToHisto.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief Mutlpily old calibration parameters by new iteration
///
/// How to run :
///   aliroot -b -q "macros/MultiplyPi0CalibrationFactors_TextToHisto.C" >& multiplyPi0CalibrationFactors_TextToHisto_forIterXXX.out &
///
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
void MultiplyPi0CalibrationFactors_TextToHisto()
{int a,b,c,d,aTot,bTot,cTot,dTot,i,j,jFile,kNbFiles,iSM;
 float e,eTot;
 double minHist,maxHist,minHistProduct;
 
 
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
 
 
 //CUSTOMIZE customize :
 kNbFiles=5;
 //kNbFiles=2; To create special files like that for Anders.

 
 char psfile[150];
 FILE *txtFileOut = NULL;
 
 //Standard calibration process :
 txtFileOut  =   fopen(Form("multiplyPi0CalibrationFactors_TextToHisto_forIter%d.txt",kNbFiles),"w");
 TFile * f = new TFile(Form("multiplyPi0CalibrationFactors_TextToHisto_forIter%d.root",kNbFiles),"recreate");
 sprintf(psfile,Form("multiplyPi0CalibrationFactors_TextToHisto_forIter%d.ps",kNbFiles));
  
 //calib 2012 coeff for EMCAL, and DCAL + third SMs coeffs at 1.0 :
 /*txtFileOut  =   fopen("multiplyPi0CalibrationFactors_TextToHisto_forEMCALandDCALwith2012coeffs.txt","w");
 TFile * f = new TFile("multiplyPi0CalibrationFactors_TextToHisto_forEMCALandDCALwith2012coeffs.root","recreate");
 sprintf(psfile,"multiplyPi0CalibrationFactors_TextToHisto_forEMCALandDCALwith2012coeffs.ps");*/
 
 //calib 2012 coeff for EMCAL (pi0 only, no E-spectrum calib), and DCAL + third SMs coeffs at 1.0 :
 /*txtFileOut  =   fopen("multiplyPi0CalibrationFactors_TextToHisto_EMCALcoeffs2012pizOnlyNoEspectra_DCALandThirdsAllOne.txt","w");
 TFile * f = new TFile("multiplyPi0CalibrationFactors_TextToHisto_EMCALcoeffs2012pizOnlyNoEspectra_DCALandThirdsAllOne.root","recreate");
 sprintf(psfile,"multiplyPi0CalibrationFactors_TextToHisto_EMCALcoeffs2012pizOnlyNoEspectra_DCALandThirdsAllOne.ps");*/
 
 //calib 2012 coeff for EMCAL (pi0 only, no E-spectrum calib), and DCAL + third SMs coeffs as obtained after pass2 :
 /*txtFileOut  =   fopen("multiplyPi0CalibrationFactors_TextToHisto_forIter3WithAndersValues.txt","w");
 TFile * f = new TFile("multiplyPi0CalibrationFactors_TextToHisto_forIter3WithAndersValues.root","recreate");
 sprintf(psfile,"multiplyPi0CalibrationFactors_TextToHisto_forIter3WithAndersValues.ps");*/
 
  //defineMyPalette40All(30,5);
  
 //Coeff files to be read and multiplied :
 FILE **calibCoeffsInput;
 calibCoeffsInput = new FILE*[kNbFiles];
 //CUSTOMIZE customize :
 //calibCoeffsInput[0] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass0_mergedDCALandThirdSMs/output_calibPi0_coeffs_clean.txt","r");
 calibCoeffsInput[0] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass0_mergedDCALandThirdSMs/output_calibPi0_coeffs_clean_veryHighTowers.txt","r"); //Chosen among various pass0, for pass1.
 //calibCoeffsInput[1] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass1_DCALandThirdSMs/output_calibPi0_coeffs_clean.txt","r");
 calibCoeffsInput[1] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass1_DCALandThirdSMs/output_calibPi0_coeffs_clean_veryHighTowers.txt","r"); //Chosen among various pass1, for pass2.
 //calibCoeffsInput[2] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_coeffs_clean.txt","r"); //Chosen in 2015 among various pass2, for pass3.
 //calibCoeffsInput[2] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_coeffs_clean_veryHighGains.txt","r");
 //calibCoeffsInput[2] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_coeffs_clean_finalFile4HVcalculation.txt","r"); //Chosen in 2015 among various pass2, for the 2015 new HV calculation.
 calibCoeffsInput[2] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_coeffs_clean_AndersCorrectionTowersVeryFar.txt","r"); //Chosen in 2016 among various pass2, for pass3 : includes correction for towers with peak very far from PDG mass.
 calibCoeffsInput[3] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass3_DCALandThirdSMsWithAndersCoeffsTowersVeryFarFromPDG/output_calibPi0_coeffs_cleanForMultiply.txt","r"); //Chosen among various pass3, for pass4.
 calibCoeffsInput[4] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass4_DCALandEMCALThirds/output_calibPi0_coeffs_cleanForMultiply.txt","r"); //From pass4, for pass5.
 
 //calib 2012 coeff for EMCAL, and DCAL + third SMs coeffs at 1.0 :
 //(The EMCAL coeffs come from /cebaf/cebaf/EMCAL/calibPi0/RecalibrationFactors2012_10SM_final.txt and include the E spectra coeffs).
 //calibCoeffsInput[0] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/calib2012onEMCALandDCAL/multiplyPi0CalibrationFactors2012WithFakeValuesForDCALandThirds_TextToHisto_Final.txt","r");
 
 //calib 2012 coeff for EMCAL (pi0 only, no E-spectrum calib), and DCAL + third SMs coeffs at 1.0 :
 //calibCoeffsInput[0] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/calib2012onEMCALandDCAL/RecalibrationFactors2012_10SM_iter8_addDCALandThirdsAllOne.txt","r");
 
 //calib 2012 coeff for EMCAL (pi0 only, no E-spectrum calib), and DCAL + third SMs coeffs as obtained after pass2 :
 //calibCoeffsInput[0] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/calib2012onEMCALandDCAL/RecalibrationFactors2012_10SM_iter8_addDCALandThirdsAfterPass2.txt","r");
 //calibCoeffsInput[1] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/studiesWith2015ppData/testInvMassPeakMeanVsSampleUsed/output_Sample2/output_calibPi0_coeffs_clean.txt","r");
    
 //calib 2012 coeff for EMCAL (pi0 only, no E-spectrum calib), and DCAL + third SMs coeffs as obtained after pass2 + Anders values for towers with peak far from pi0 PDG mass:
 /*calibCoeffsInput[0] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/calib2012onEMCALandDCAL/RecalibrationFactors2012_10SM_iter8_addDCALandThirdsAfterPass2.txt","r");
 calibCoeffsInput[1] = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_coeffs_clean_WithAndersValues.txt","r");*/
 
 
 //Histoes definition :
 TH2F **h2Coeffs;
 h2Coeffs = new TH2F*[kNbSMtot];
 for (i=0;i<kNbSMtot;i++)
    {h2Coeffs[i] = new TH2F(Form("EMCALRecalFactors_SM%d",i),Form("EMCALRecalFactors_SM%d",i),kNbColMax,0,kNbColMax,kNbRowMax,0,kNbRowMax); //All SMs use the same size.
     }
 
 TH2F **hSpace;
 TH1F **hDistr; //Coeff for all SMs, per step
 TH1F **hDistrPerSM; //Coeff for all steps, per SM
 TH1F **hDistrNoone; //Coeff for all SMs, per step, exclude towers at 1.
 TH1F **hDistrNoonePerSM; //Coeff for all steps, per SM, exclude towers at 1.
 hSpace = new TH2F*[(kNbFiles+1)*kNbSMtot];
 hDistr = new TH1F*[(kNbFiles+1)*2];
 hDistrPerSM = new TH1F*[kNbSMtot];
 hDistrNoone = new TH1F*[(kNbFiles+1)*2];
 hDistrNoonePerSM = new TH1F*[kNbSMtot];
 for (i=0;i<kNbSMtot;i++)
    {hDistrPerSM[i] = new TH1F(Form("hDistrPerSM%d",i),Form("hDistrPerSM%d",i),100,0.6,1.4);
     hDistrPerSM[i]->SetXTitle("Coefficient");
     hDistrPerSM[i]->SetYTitle("Counts");
     hDistrPerSM[i]->SetStats(0);
     hDistrNoonePerSM[i] = new TH1F(Form("hDistrNoonePerSM%d",i),Form("hDistrNoonePerSM%d",i),100,0.6,1.4);
     hDistrNoonePerSM[i]->SetXTitle("Coefficient");
     hDistrNoonePerSM[i]->SetYTitle("Counts");
     hDistrNoonePerSM[i]->SetStats(0);
     for (j=0;j<kNbFiles+1;j++)
        {hSpace[(kNbFiles+1)*i+j] = new TH2F(Form("hSpace_%d_%d",i,j),Form("hSpace_%d_%d",i,j),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
         hSpace[(kNbFiles+1)*i+j]->SetXTitle("Column");
         hSpace[(kNbFiles+1)*i+j]->SetYTitle("Row");
         hSpace[(kNbFiles+1)*i+j]->SetStats(0);
         hSpace[(kNbFiles+1)*i+j]->SetContour(30);
         }
     }
 for (j=0;j<kNbFiles+1;j++)
    {hDistr[j] = new TH1F(Form("hDistr_%d",j),Form("hDistr_%d",j),100,0.0,3.0);
     hDistr[j]->SetMaximum(20.);
     hDistr[(kNbFiles+1)+j] = new TH1F(Form("hDistrZm_%d",j),Form("hDistrZm_%d",j),100,0.85,1.15);
     hDistrNoone[j] = new TH1F(Form("hDistrNoone_%d",j),Form("hDistrNoone_%d",j),100,0.0,3.0);
     hDistrNoone[j]->SetMaximum(20.);
     hDistrNoone[(kNbFiles+1)+j] = new TH1F(Form("hDistrNooneZm_%d",j),Form("hDistrNooneZm_%d",j),100,0.85,1.15);
     for (i=0;i<2;i++)
        {hDistr[(kNbFiles+1)*i+j]->SetXTitle("Coefficient");
         hDistr[(kNbFiles+1)*i+j]->SetYTitle("Counts");
         hDistr[(kNbFiles+1)*i+j]->SetStats(0);
         hDistr[(kNbFiles+1)*i+j]->SetLineColor(4);
         hDistrNoone[(kNbFiles+1)*i+j]->SetXTitle("Coefficient");
         hDistrNoone[(kNbFiles+1)*i+j]->SetYTitle("Counts");
         hDistrNoone[(kNbFiles+1)*i+j]->SetStats(0);
         hDistrNoone[(kNbFiles+1)*i+j]->SetLineColor(4);
         }
     }

 //Calculation of coeffs :
 for (iSM=0;iSM<kNbSMtot;iSM++)
    {for(i = 0;i<kTabNbCol[SMdetType[iSM]]*kTabNbRow[SMdetType[iSM]];i++)
        {jFile=0;
         fscanf(calibCoeffsInput[jFile], " %d  %d  %d   %d   %f\n", &a, &b, &c,&d,&e);
         //printf( " %d  %d  %d   %d   %f\n", a, b, c,d,e);
         if (b != iSM) printf("$$$ File reading out of sync : received SM %d while expected %d.\n",b,iSM);
         if (e < 0.) e=1.00; //When coeff is used as a flag for towers not to be touched again.
         hSpace[(kNbFiles+1)*b+jFile]->Fill(c,d,e);
         hDistr[jFile]->Fill(e);
         hDistr[(kNbFiles+1)+jFile]->Fill(e);
         if (e != 1.)
            {hDistrNoone[jFile]->Fill(e);
             hDistrNoone[(kNbFiles+1)+jFile]->Fill(e);
             }
         aTot=a;
         bTot=b;
         cTot=c;
         dTot=d;
         eTot=e;
         for (jFile=1;jFile<kNbFiles;jFile++)
            {fscanf(calibCoeffsInput[jFile], " %d  %d  %d   %d   %f\n", &a, &b, &c,&d,&e);
             if (bTot != b) printf("$$$ Files reading out of sync : (%d,%d,%d) vs (%d,%d,%d)\n",b,c,d,bTot,cTot,dTot);
             if (cTot != c) printf("$$$ Files reading out of sync : (%d,%d,%d) vs (%d,%d,%d)\n",b,c,d,bTot,cTot,dTot);
             if (dTot != d) printf("$$$ Files reading out of sync : (%d,%d,%d) vs (%d,%d,%d)\n",b,c,d,bTot,cTot,dTot);
             if (e < 0.) e=1.00; //When coeff is used as a flag for towers not to be touched again.
             hSpace[(kNbFiles+1)*b+jFile]->Fill(c,d,e);
             hDistr[jFile]->Fill(e);
             hDistr[(kNbFiles+1)+jFile]->Fill(e);
             if (e != 1.)
                {hDistrNoone[jFile]->Fill(e);
                 hDistrNoone[(kNbFiles+1)+jFile]->Fill(e);
                 }
             eTot*=e;
             /*if (jFile == (kNbFiles-1))
                {if (e == 1) eTot=1.0; //If tower not trusted at last iteration, don't trust it for the whole process.
                 }*/
             //We take the previous lines out and use this only at the very last iteration, in a special code similar to this one. In this way, towers can be made converge, then untrusted for a while and re-converged if needed.
             }
         hSpace[(kNbFiles+1)*bTot+kNbFiles]->Fill(cTot,dTot,eTot);
         hDistr[kNbFiles]->Fill(eTot);
         hDistr[(kNbFiles+1)+kNbFiles]->Fill(eTot);
         hDistrPerSM[bTot]->Fill(eTot);
         if (eTot != 1.)
            {hDistrNoone[kNbFiles]->Fill(eTot);
             hDistrNoone[(kNbFiles+1)+kNbFiles]->Fill(eTot);
             hDistrNoonePerSM[bTot]->Fill(eTot);
             }
         fprintf(txtFileOut,"%d %d %d %f\n",bTot,cTot,dTot,eTot);
         if(eTot==0) printf("\n###### Coeff cannot be 0 ! Please check tower : %d %d %d %f #######\n\n",bTot,cTot,dTot,eTot);
         else h2Coeffs[bTot]->SetBinContent(cTot,dTot,1./eTot);
         }
     }


 // Draw plots :

 const int cWidth=500;
 const int cHeight=(int)(500*(29./21.));
 TCanvas *c1 = new TCanvas("c1","EMCal cosmics analysis",cWidth,cHeight);
 TPostScript *ps = new TPostScript(psfile,111);
 
 //Spatial distrib of coeffs, per SM, for each pass and integrated over all passes.
 for (i=0;i<kNbSMtot;i++)
    {minHist=1.;
     maxHist=1.;
     minHistProduct=1.;
     for (j=0;j<kNbFiles;j++)
        {//if (hSpace[(kNbFiles+1)*i+j]->GetMinimum() < minHist) minHist=hSpace[(kNbFiles+1)*i+j]->GetMinimum();
         if (hSpace[(kNbFiles+1)*i+j]->GetMaximum() > maxHist) maxHist=hSpace[(kNbFiles+1)*i+j]->GetMaximum();
         for (int iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
            {for (int iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
                {if (hSpace[(kNbFiles+1)*i+j]->GetBinContent(iCol+1,iRow+1) < minHist) minHist=hSpace[(kNbFiles+1)*i+j]->GetBinContent(iCol+1,iRow+1);
                 }
             }
         }
     for (int iCol=0;iCol<kTabNbCol[SMdetType[i]];iCol++)
        {for (int iRow=0;iRow<kTabNbRow[SMdetType[i]];iRow++)
            {if (hSpace[(kNbFiles+1)*i+kNbFiles]->GetBinContent(iCol+1,iRow+1) < minHistProduct) minHistProduct=hSpace[(kNbFiles+1)*i+kNbFiles]->GetBinContent(iCol+1,iRow+1);
             }
         }
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (j=0;j<kNbFiles;j++)
        {c1->cd(j+1);
         hSpace[(kNbFiles+1)*i+j]->SetMinimum(minHist);
         hSpace[(kNbFiles+1)*i+j]->SetMaximum(maxHist);
         hSpace[(kNbFiles+1)*i+j]->Draw("COLZ");
         hSpace[(kNbFiles+1)*i+j]->Write();
         }
     c1->cd(kNbFiles+1);
     hSpace[(kNbFiles+1)*i+kNbFiles]->SetMinimum(minHistProduct);
     hSpace[(kNbFiles+1)*i+kNbFiles]->Draw("COLZ");
     hSpace[(kNbFiles+1)*i+kNbFiles]->Write();
     c1->Update();
     }
 
 //1-D distribs of all calib coeffs, per pass and all passes integrated.
 //Page "1" : zoomed (in y) ; page "0" : not zoomed.
 for (i=0;i<2;i++)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (j=0;j<kNbFiles;j++)
        {c1->cd(j+1);
         /*printf("Histo coeffs pass %d ",j+1);
         switch (i)
            {case 0 : printf("(large range) ");
                      break;
             case 1 : printf("(reduced range) ");
                      break;
             default : printf("##### unknown histo i=%d ##### \n",i);
             }
         printf(": mean = %f,   RMS = %f\n",hDistr[(kNbFiles+1)*i+j]->GetMean(),hDistr[(kNbFiles+1)*i+j]->GetRMS());*/
         hDistr[(kNbFiles+1)*i+j]->Draw();
         hDistr[(kNbFiles+1)*i+j]->Write();
         }
     c1->cd(kNbFiles+1);
     /*printf("Histo 'product of all coeffs' ");
     switch (i)
        {case 0 : printf("(large range) ");
                  break;
         case 1 : printf("(reduced range) ");
                  break;
         default : printf("##### unknown histo i=%d ##### \n",i);
         }
     printf(": mean = %f,   RMS = %f\n",hDistr[(kNbFiles+1)*i+kNbFiles]->GetMean(),hDistr[(kNbFiles+1)*i+kNbFiles]->GetRMS());*/
     /*if (i == 1)
        {TF1 *func = new TF1("func","gaus(0)",0.85,1.15);
         func->SetParLimits(0,10.,500.);
         func->SetParLimits(1,0.9,1.1);
         func->SetParLimits(2,0.01,0.2);
         func->SetParameter(0,hDistr[(kNbFiles+1)+kNbFiles]->GetMaximum());
         func->SetParameter(1,1.);
         func->SetParameter(2,0.1);
         printf("\n");
         hDistr[(kNbFiles+1)+kNbFiles]->Fit(func,"BRM"); //Option I crashe
         //printf("Histo moy = %f,  RMS = %f\n",hDistr[(kNbFiles+1)+kNbFiles]->GetMean(),hDistr[(kNbFiles+1)+kNbFiles]->GetRMS());
         printf("\nFit over 'product of all coeffs' (reduced range) : mean = %f,   RMS = %f\n",func->GetParameter(1),func->GetParameter(2));
         }*/
     hDistr[(kNbFiles+1)*i+kNbFiles]->Draw();
     hDistr[(kNbFiles+1)*i+kNbFiles]->Write();
     c1->Update();
     }
 
 //Same thing now without the towers with coeff at 1 :
 for (i=0;i<2;i++)
    {ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (j=0;j<kNbFiles;j++)
        {c1->cd(j+1);
         printf("Histo coeffs pass %d ",j+1);
         switch (i)
            {case 0 : printf("(large range) ");
                      break;
             case 1 : printf("(reduced range) ");
                      break;
             default : printf("##### unknown histo i=%d ##### \n",i);
             }
         printf(": mean = %f,   RMS = %f\n",hDistrNoone[(kNbFiles+1)*i+j]->GetMean(),hDistrNoone[(kNbFiles+1)*i+j]->GetRMS());
         hDistrNoone[(kNbFiles+1)*i+j]->Draw();
         hDistrNoone[(kNbFiles+1)*i+j]->Write();
         }
     c1->cd(kNbFiles+1);
     printf("Histo 'product of all coeffs' ");
     switch (i)
        {case 0 : printf("(large range) ");
                  break;
         case 1 : printf("(reduced range) ");
                  break;
         default : printf("##### unknown histo i=%d ##### \n",i);
         }
     printf(": mean = %f,   RMS = %f\n",hDistrNoone[(kNbFiles+1)*i+kNbFiles]->GetMean(),hDistrNoone[(kNbFiles+1)*i+kNbFiles]->GetRMS());
     if (i == 1)
        {TF1 *func = new TF1("func","gaus(0)",0.85,1.15);
         func->SetParLimits(0,10.,500.);
         func->SetParLimits(1,0.9,1.1);
         func->SetParLimits(2,0.01,0.2);
         func->SetParameter(0,hDistrNoone[(kNbFiles+1)+kNbFiles]->GetMaximum());
         func->SetParameter(1,1.);
         func->SetParameter(2,0.1);
         printf("\n");
         hDistrNoone[(kNbFiles+1)+kNbFiles]->Fit(func,"BRM"); //Option I crashe
         //printf("Histo moy = %f,  RMS = %f\n",hDistrNoone[(kNbFiles+1)+kNbFiles]->GetMean(),hDistrNoone[(kNbFiles+1)+kNbFiles]->GetRMS());
         printf("\nFit over 'product of all coeffs' (reduced range) : mean = %f,   RMS = %f\n",func->GetParameter(1),func->GetParameter(2));
         }
     hDistrNoone[(kNbFiles+1)*i+kNbFiles]->Draw();
     hDistrNoone[(kNbFiles+1)*i+kNbFiles]->Write();
     c1->Update();
     }
 
 //Coeff (all iterations multiplied) distrib SM per SM
 printf("\n\nHistoes 'product of all coeffs' for each SM :\n");

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=0;j<kNbSMEMCAL;j++)
    {c1->cd(j+1);
     //printf("  SM %d : mean = %f,   RMS = %f\n",j,hDistrPerSM[j]->GetMean(),hDistrPerSM[j]->GetRMS());
     hDistrPerSM[j]->Draw();
     hDistrPerSM[j]->Write();
     }
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=kNbSMEMCAL;j<kNbSMtot;j++)
    {c1->cd(j-kNbSMEMCAL+1);
     //printf("  SM %d : mean = %f,   RMS = %f\n",j,hDistrPerSM[j]->GetMean(),hDistrPerSM[j]->GetRMS());
     hDistrPerSM[j]->Draw();
     hDistrPerSM[j]->Write();
     }
 c1->Update();

 //Same thing excluding the towers which have coeff equal to 1 :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=0;j<kNbSMEMCAL;j++)
    {c1->cd(j+1);
     printf("  SM %d : mean = %f,   RMS = %f\n",j,hDistrNoonePerSM[j]->GetMean(),hDistrNoonePerSM[j]->GetRMS());
     hDistrNoonePerSM[j]->Draw();
     hDistrNoonePerSM[j]->Write();
     }
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=kNbSMEMCAL;j<kNbSMtot;j++)
    {c1->cd(j-kNbSMEMCAL+1);
     printf("  SM %d : mean = %f,   RMS = %f\n",j,hDistrNoonePerSM[j]->GetMean(),hDistrNoonePerSM[j]->GetRMS());
     hDistrNoonePerSM[j]->Draw();
     hDistrNoonePerSM[j]->Write();
     }
 c1->Update();

 
 ps->Close();
 
 for (i=0;i<kNbSMtot;i++) h2Coeffs[i]->Write();

 f->Close();
 
 fclose(txtFileOut);
 
 for (jFile=kNbFiles-1;jFile>-1;jFile--) fclose(calibCoeffsInput[jFile]);

 return;
 
 }




