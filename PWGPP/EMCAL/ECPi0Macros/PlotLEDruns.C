#include <signal.h>
//#include <sys/time.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iosfwd>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>

#include <map>
#include <utility>
#include <iterator>


#include <TApplication.h>
#include <TChain.h> 
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h> 
#include <TCanvas.h>
#include <TFrame.h>
#include <TPostScript.h>
#include <TLine.h> 
#include <TGaxis.h> 
#include <TStyle.h> 
#include <TGraphErrors.h> 
#include <TMath.h>
#include <TMatrixF.h>
#include <TText.h>


//#include "/cebaf/faivre/recherche/utilities/defineMyPalette2011.C"


namespace std {} using namespace std;


bool wled=true;
bool ledonly=true;
int stripPlot=-1;
char RunlistFilename[200]="";


const double coefFactorWanted=0.0162;





enum detType {kEMCAL,kEMCALthird,kDCAL,kDCALthird};
char detTypeString[][100]={"EMCAL","EMCALthird","DCAL","DCALthird"};

// voir http://dsilverm.web.cern.ch/dsilverm/fee/addrP2.html
char SMP2Name[][100]={"SMA0","SMC0","SMA1","SMC1","SMA2","SMC2","SMA3","SMC3","SMA4","SMC4","SMA5","SMC5","SMA9","SMC9","SMA10","SMC10","SMA11","SMC11","SMA12","SMC12"};
int SMdetType[]={kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCALthird,kEMCALthird,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCALthird,kDCALthird};



///
/// \file Pi0CalibInvMassAnalysis3.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief Pi0CalibInvMassAnalysis3
///
/// How to compile :
///   .L plotLEDruns.C++
///
/// How to run (uncompiled) :
///
/// 1st and 3rd run = same old settings ; 2nd run = new settings to be tested.
///
///   .x macros/plotLEDruns.C(run1,run2,run3,"toto")
///
/// CAUTION : some things at the beginning are "To be checked / changed by user" (search "customize").
///
/// Mapping conventions :
/// Read OCDB files (in Alice-offline mapping) and keep them in Alice-offline mapping.
/// Read amplitude histograms (in Grenoble/electronic mapping) and flip them to Alice-offline mapping.
/// Row/col TH2F plots are in Alice-offline mapping.
/// "Calculations" (i.e. things inside the code) are all done in Alice-offline mapping.
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
void PlotLEDruns(int run1,int run2,int run3,TString ident="")
{int i,j,iRun,iSM,iCol,iRow,iColGlob,iRowGlob,tabCol[3],NumRun[3],cmptIntg[3],cmptAmpl[3],cmptEmismatch[3];
 char name[150],name2[150],fchNameOCDBcoeffsFile[200];
 double maxloc,ampli0,ampli1,ampli2,paramFactor,E0,E1,E2;
 TFile *rootFileIn;
 
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
 
  //defineMyPalette2011(30,5);
 
 double coeffOCDBbefore[kNbSMtot][kNbColMax][kNbRowMax]={{{0.}}};
 double coeffOCDBafter[kNbSMtot][kNbColMax][kNbRowMax]={{{0.}}};

 //CUSTOMIZE customize the runMode and the cut values :
 int runMode=0; // 0 = compare HV changes ; 1 = compare various runs with the same HV.
 double cutWarnERatio=0.05;
 double cutWarnAmpRatio=0.05;
 double cutMinIntg=5.;
 double cutMinAmp=20.;
 double maxIntg=2800.;
 double maxAmp=1100.;
 double maxE=20.;
 double multFactGraphAmp=0.7;
 double multFactGraphE=0.6;
 
 NumRun[0]=run1;
 NumRun[1]=run2;
 NumRun[2]=run3;
 
 tabCol[0]=kGreen+2;
 tabCol[1]=kRed;
 tabCol[2]=kBlue;
 
 TH1F *hampHighGainTow[3][kNbSMtot][kNbRowMax][kNbColMax]; //(SM,row,col), sic (not col,row).
 
 switch (runMode)
    {case 0 : printf("\n--> Run mode = compare HV changes (runs have different HVs).\n\n");
              break;
     case 1 : printf("\n--> Run mode = compare various runs with the same HV.\n\n");
              break;
     default : printf("Unknown run mode %d, exiting.\n",runMode);
               return;
     }
 
 char fchNameBase[200],fchNameRoot[200],fchNamePs[200],fchNameOut[200];
 FILE *inFile,*outFile,*OCDBcoeffsFile;
 if (ident.IsNull()) sprintf(fchNameBase,"output_plotLEDruns_runs%d.%d.%d",run1,run2,run3);
    else sprintf(fchNameBase,"output_plotLEDruns_%s",ident.Data());
 sprintf(fchNameRoot,"%s.root",fchNameBase);
 sprintf(fchNamePs,"%s.ps",fchNameBase);
 sprintf(fchNameOut,"%s.out",fchNameBase);
 outFile=fopen(fchNameOut,"w");

 fprintf(outFile,"\nAnalyzed run numbers :  %d  %d  %d\n",run1,run2,run3);
 fprintf(outFile,"Chosen run mode :  %d\n",runMode);
 
 // Load the OCDB factors :
 //CUSTOMIZE customize every year : old OCDB coeffs :
 //sprintf(fchNameOCDBcoeffsFile,"/cebaf/faivre/recherche/calibPi0/recalculateEMCAL_HV_Apr2013/OCDBparamsAfterCalib2012_allSMs.txt"); //This is WRONG
 sprintf(fchNameOCDBcoeffsFile,"/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/OCDBparamsAfterCalib2015.txt");
 fprintf(outFile,"Load old OCDB coeffs from :  %s\n",fchNameOCDBcoeffsFile);
 OCDBcoeffsFile = fopen(fchNameOCDBcoeffsFile,"r");
 if (!OCDBcoeffsFile) {printf("File %s can not be found\n",fchNameOCDBcoeffsFile);exit(-1);}
 while (fscanf(OCDBcoeffsFile,"%d %d %d %lf",&iSM,&iCol,&iRow,&paramFactor)>0)
    {coeffOCDBbefore[iSM][iCol][iRow]=paramFactor;
     }
 fclose(OCDBcoeffsFile);
 for (iSM=0;iSM<kNbSMtot;iSM++)
    {//CUSTOMIZE customize every HV calculation iteration : new OCDB coeffs :
     //sprintf(fchNameOCDBcoeffsFile,"/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/output_HVrecalculation_pass2/%s/NewParamFactor.txt",SMP2Name[iSM]);
     sprintf(fchNameOCDBcoeffsFile,"/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/output_HVrecalculation_pass2_handMadeFinalCorrections_createSetBiasScripts/%s/NewParamFactor.txt",SMP2Name[iSM]);
     fprintf(outFile,"Load new OCDB coeffs from :  %s\n",fchNameOCDBcoeffsFile);
     OCDBcoeffsFile = fopen(fchNameOCDBcoeffsFile,"r");
     if (!OCDBcoeffsFile)
        {printf("File %s can not be found\nWill use fake coeffs for SM %s !\n",fchNameOCDBcoeffsFile,SMP2Name[iSM]);
         fprintf(outFile,"File %s can not be found\nWill use fake coeffs for SM %s !\n",fchNameOCDBcoeffsFile,SMP2Name[iSM]);
         }
        else
        {while (fscanf(OCDBcoeffsFile,"%d %d %d %lf",&iSM,&iCol,&iRow,&paramFactor)>0) coeffOCDBafter[iSM][iCol][iRow]=paramFactor;
         fclose(OCDBcoeffsFile);
         }
     }
 
 const int cWidth=500;
 const int cHeight=(int)(500*(30./21.));
 TCanvas *c1 = new TCanvas("c1","EMCal cosmics analysis",cWidth,cHeight);
 TPostScript *ps = new TPostScript(fchNamePs,111);
 
 /*From http://root.cern.ch/root/HowtoHistogram.html :
 "Like any other ROOT objects, histograms can be written to a file. Example: hist->Write(). When an histogram is created, a reference to this histogram is added in the list of objects of the current directory in memory. If no file is opened at the time of the creation, the histogram is added to the default list of objects in memory. If a file exists, the histogram is added to the list of memory objects associated to this file."
 That's the reason why, when I close rootFileIn, I encounter problems as soon as I try to access hampHighGainTow[iRun][iSM][i][j] afterwards.
 But : if I close rootFileIn after the loop, also the hAllSpace... histoes are defined while rootFileIn is still opened, and thus associated to it, and thus causes problems when we need to manipulate them afterwards (Draw,...). -- although the segm viol doesn't occur systematically (but is reproducable for a given code config).
 Instanciating the hampHighGainTow[iRun][iSM][i][j]'s before rootFileIn is opened doesn't solve. Cloning works but doesn't solve either.
 So : I must close nothing...
 */
 
 for (iRun=0;iRun<3;iRun++)
    {for (iSM=0;iSM<kNbSMtot;iSM++)
        {sprintf(name,"/cebaf/cebaf/EMCAL/anaLHC/ledana/root4gain/gain_SM%d_%d.root",iSM,NumRun[iRun]);
         printf("Opening %s\n",name);
         rootFileIn=new TFile(name,"READONLY");
         if (rootFileIn->IsZombie()) {exit(-1);}
         
         for (iCol=0;iCol<kTabNbCol[SMdetType[iSM]];iCol++)
            {for (iRow=0;iRow<kTabNbRow[SMdetType[iSM]];iRow++) //Flip from Grenoble/electronic mapping to Alice-offline mapping.
                {if ((iSM%2) == 0) sprintf(name,"hampHighGain_r%i_c%i",iRow,iCol); //Side A.
                    else sprintf(name,"hampHighGain_r%i_c%i",(kTabNbRow[SMdetType[iSM]]-1)-iRow,(kTabNbCol[SMdetType[iSM]]-1)-iCol); //Side C.
                 hampHighGainTow[iRun][iSM][iRow][iCol]=(TH1F*)rootFileIn->Get(name);
                 if( !(hampHighGainTow[iRun][iSM][iRow][iCol]))
                    {printf("hist %s not found in %s. Aborting\n",name,rootFileIn->GetName());
                     exit(-1);
                     }
                 }
             }
         }
     }

 TH1F **hDistrIntg;
 TH1F **hDistrAmp;
 TH1F **hDistrE;
 TH1F **hDistrIntgRatio;
 TH1F **hDistrAmpRatio;
 TH1F **hDistrERatio;
 TH1F **hDistrIntgRatioZm;
 TH1F **hDistrAmpRatioZm;
 TH1F **hDistrERatioZm;
 TH2F **hAllSpaceEMCALIntg;
 TH2F **hAllSpaceEMCALIntgRatio;
 TH2F **hAllSpaceEMCALAmp;
 TH2F **hAllSpaceEMCALAmpRatio;
 TH2F **hAllSpaceEMCALE;
 TH2F **hAllSpaceEMCALERatio;
 TH2F **hAllSpaceDCALIntg;
 TH2F **hAllSpaceDCALIntgRatio;
 TH2F **hAllSpaceDCALAmp;
 TH2F **hAllSpaceDCALAmpRatio;
 TH2F **hAllSpaceDCALE;
 TH2F **hAllSpaceDCALERatio;
 TH2F **h2Amp;
 TH2F **h2E;
 TH2F **h2RatioVsIntg;
 TH2F **h2RatioVsAmp;
 hDistrIntg = new TH1F*[3*(kNbSMtot+1)];
 hDistrAmp = new TH1F*[3*(kNbSMtot+1)];
 hDistrE = new TH1F*[3*(kNbSMtot+1)];
 hDistrIntgRatio = new TH1F*[2*(kNbSMtot+1)];
 hDistrAmpRatio = new TH1F*[2*(kNbSMtot+1)];
 hDistrERatio = new TH1F*[2*(kNbSMtot+1)];
 hDistrIntgRatioZm = new TH1F*[2*(kNbSMtot+1)];
 hDistrAmpRatioZm = new TH1F*[2*(kNbSMtot+1)];
 hDistrERatioZm = new TH1F*[2*(kNbSMtot+1)];
 hAllSpaceEMCALIntg = new TH2F*[3];
 hAllSpaceEMCALIntgRatio = new TH2F*[2];
 hAllSpaceEMCALAmp = new TH2F*[3];
 hAllSpaceEMCALAmpRatio = new TH2F*[2];
 hAllSpaceEMCALE = new TH2F*[3];
 hAllSpaceEMCALERatio = new TH2F*[2];
 hAllSpaceDCALIntg = new TH2F*[3];
 hAllSpaceDCALIntgRatio = new TH2F*[2];
 hAllSpaceDCALAmp = new TH2F*[3];
 hAllSpaceDCALAmpRatio = new TH2F*[2];
 hAllSpaceDCALE = new TH2F*[3];
 hAllSpaceDCALERatio = new TH2F*[2];
 h2Amp = new TH2F*[2*(kNbSMtot+1)];
 h2E = new TH2F*[2*(kNbSMtot+1)];
 h2RatioVsIntg = new TH2F*[6*(kNbSMtot+1)];
 h2RatioVsAmp = new TH2F*[6*(kNbSMtot+1)];
 //Figure "3" here stands for number of analyzed runs ; "2" for the ratios or comparisons, because run2/run1 and run3/run1.
 for (j=0;j<3;j++)
    {hAllSpaceEMCALIntg[j] = new TH2F(Form("hAllSpaceEMCALIntg_%d",j),Form("hAllSpaceEMCALIntg_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALIntg[j]->SetXTitle("Column");
     hAllSpaceEMCALIntg[j]->SetYTitle("Row");
     hAllSpaceEMCALIntg[j]->SetStats(0);
     hAllSpaceEMCALIntg[j]->SetContour(30);
     hAllSpaceEMCALAmp[j] = new TH2F(Form("hAllSpaceEMCALAmp_%d",j),Form("hAllSpaceEMCALAmp_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALAmp[j]->SetXTitle("Column");
     hAllSpaceEMCALAmp[j]->SetYTitle("Row");
     hAllSpaceEMCALAmp[j]->SetStats(0);
     hAllSpaceEMCALAmp[j]->SetContour(30);
     hAllSpaceEMCALE[j] = new TH2F(Form("hAllSpaceEMCALE_%d",j),Form("hAllSpaceEMCALE_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALE[j]->SetXTitle("Column");
     hAllSpaceEMCALE[j]->SetYTitle("Row");
     hAllSpaceEMCALE[j]->SetStats(0);
     hAllSpaceEMCALE[j]->SetContour(30);
     hAllSpaceDCALIntg[j] = new TH2F(Form("hAllSpaceDCALIntg_%d",j),Form("hAllSpaceDCALIntg_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALIntg[j]->SetXTitle("Column");
     hAllSpaceDCALIntg[j]->SetYTitle("Row");
     hAllSpaceDCALIntg[j]->SetStats(0);
     hAllSpaceDCALIntg[j]->SetContour(30);
     hAllSpaceDCALAmp[j] = new TH2F(Form("hAllSpaceDCALAmp_%d",j),Form("hAllSpaceDCALAmp_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALAmp[j]->SetXTitle("Column");
     hAllSpaceDCALAmp[j]->SetYTitle("Row");
     hAllSpaceDCALAmp[j]->SetStats(0);
     hAllSpaceDCALAmp[j]->SetContour(30);
     hAllSpaceDCALE[j] = new TH2F(Form("hAllSpaceDCALE_%d",j),Form("hAllSpaceDCALE_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALE[j]->SetXTitle("Column");
     hAllSpaceDCALE[j]->SetYTitle("Row");
     hAllSpaceDCALE[j]->SetStats(0);
     hAllSpaceDCALE[j]->SetContour(30);
     for (i=0;i<kNbSMtot+1;i++)
        {hDistrIntg[3*i+j] = new TH1F(Form("hDistrIntg_%d_%d",i,j),Form("hDistrIntg_%d_%d",i,j),100,0.,maxIntg);
         hDistrIntg[3*i+j]->SetXTitle("Integral");
         hDistrIntg[3*i+j]->SetYTitle("Counts");
         hDistrIntg[3*i+j]->SetStats(0);
         hDistrIntg[3*i+j]->SetLineColor(tabCol[j]);
         hDistrAmp[3*i+j] = new TH1F(Form("hDistrAmp_%d_%d",i,j),Form("hDistrAmp_%d_%d",i,j),100,0.,maxAmp);
         hDistrAmp[3*i+j]->SetXTitle("Amplitude");
         hDistrAmp[3*i+j]->SetYTitle("Counts");
         hDistrAmp[3*i+j]->SetStats(0);
         hDistrAmp[3*i+j]->SetLineColor(tabCol[j]);
         hDistrE[3*i+j] = new TH1F(Form("hDistrE_%d_%d",i,j),Form("hDistrE_%d_%d",i,j),100,0.,maxE);
         hDistrE[3*i+j]->SetXTitle("Energy");
         hDistrE[3*i+j]->SetYTitle("Counts");
         hDistrE[3*i+j]->SetStats(0);
         hDistrE[3*i+j]->SetLineColor(tabCol[j]);
         }
     }
 for (j=0;j<2;j++)
    {hAllSpaceEMCALIntgRatio[j] = new TH2F(Form("hAllSpaceEMCALIntgRatio_%d",j),Form("hAllSpaceEMCALIntgRatio_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALIntgRatio[j]->SetXTitle("Column");
     hAllSpaceEMCALIntgRatio[j]->SetYTitle("Row");
     hAllSpaceEMCALIntgRatio[j]->SetStats(0);
     hAllSpaceEMCALIntgRatio[j]->SetContour(30);
     hAllSpaceEMCALAmpRatio[j] = new TH2F(Form("hAllSpaceEMCALAmpRatio_%d",j),Form("hAllSpaceEMCALAmpRatio_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALAmpRatio[j]->SetXTitle("Column");
     hAllSpaceEMCALAmpRatio[j]->SetYTitle("Row");
     hAllSpaceEMCALAmpRatio[j]->SetStats(0);
     hAllSpaceEMCALAmpRatio[j]->SetContour(30);
     hAllSpaceEMCALERatio[j] = new TH2F(Form("hAllSpaceEMCALERatio_%d",j),Form("hAllSpaceEMCALERatio_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
     hAllSpaceEMCALERatio[j]->SetXTitle("Column");
     hAllSpaceEMCALERatio[j]->SetYTitle("Row");
     hAllSpaceEMCALERatio[j]->SetStats(0);
     hAllSpaceEMCALERatio[j]->SetContour(30);
     hAllSpaceDCALIntgRatio[j] = new TH2F(Form("hAllSpaceDCALIntgRatio_%d",j),Form("hAllSpaceDCALIntgRatio_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALIntgRatio[j]->SetXTitle("Column");
     hAllSpaceDCALIntgRatio[j]->SetYTitle("Row");
     hAllSpaceDCALIntgRatio[j]->SetStats(0);
     hAllSpaceDCALIntgRatio[j]->SetContour(30);
     hAllSpaceDCALAmpRatio[j] = new TH2F(Form("hAllSpaceDCALAmpRatio_%d",j),Form("hAllSpaceDCALAmpRatio_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALAmpRatio[j]->SetXTitle("Column");
     hAllSpaceDCALAmpRatio[j]->SetYTitle("Row");
     hAllSpaceDCALAmpRatio[j]->SetStats(0);
     hAllSpaceDCALAmpRatio[j]->SetContour(30);
     hAllSpaceDCALERatio[j] = new TH2F(Form("hAllSpaceDCALERatio_%d",j),Form("hAllSpaceDCALERatio_%d",j),2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
     hAllSpaceDCALERatio[j]->SetXTitle("Column");
     hAllSpaceDCALERatio[j]->SetYTitle("Row");
     hAllSpaceDCALERatio[j]->SetStats(0);
     hAllSpaceDCALERatio[j]->SetContour(30);
     for (i=0;i<kNbSMtot+1;i++)
        {hDistrIntgRatio[2*i+j] = new TH1F(Form("hDistrIntgRatio_%d_%d",i,j),Form("hDistrIntgRatio_%d_%d",i,j),100,0.,4.);
         hDistrIntgRatio[2*i+j]->SetXTitle("Integral ratio");
         hDistrIntgRatio[2*i+j]->SetYTitle("Counts");
         hDistrIntgRatio[2*i+j]->SetStats(0);
         hDistrIntgRatio[2*i+j]->SetLineColor(tabCol[2-j]);
         hDistrAmpRatio[2*i+j] = new TH1F(Form("hDistrAmpRatio_%d_%d",i,j),Form("hDistrAmpRatio_%d_%d",i,j),100,0.,3.);
         hDistrAmpRatio[2*i+j]->SetXTitle("Amplitude ratio");
         hDistrAmpRatio[2*i+j]->SetYTitle("Counts");
         hDistrAmpRatio[2*i+j]->SetStats(0);
         hDistrAmpRatio[2*i+j]->SetLineColor(tabCol[2-j]);
         hDistrERatio[2*i+j] = new TH1F(Form("hDistrERatio_%d_%d",i,j),Form("hDistrERatio_%d_%d",i,j),100,0.,3.);
         hDistrERatio[2*i+j]->SetXTitle("Energy ratio");
         hDistrERatio[2*i+j]->SetYTitle("Counts");
         hDistrERatio[2*i+j]->SetStats(0);
         hDistrERatio[2*i+j]->SetLineColor(tabCol[2-j]);
         hDistrIntgRatioZm[2*i+j] = new TH1F(Form("hDistrIntgRatioZm_%d_%d",i,j),Form("hDistrIntgRatioZm_%d_%d",i,j),100,0.,2.0);
         hDistrIntgRatioZm[2*i+j]->SetXTitle("Integral ratio");
         hDistrIntgRatioZm[2*i+j]->SetYTitle("Counts");
         hDistrIntgRatioZm[2*i+j]->SetStats(0);
         hDistrIntgRatioZm[2*i+j]->SetLineColor(tabCol[2-j]);
         hDistrAmpRatioZm[2*i+j] = new TH1F(Form("hDistrAmpRatioZm_%d_%d",i,j),Form("hDistrAmpRatioZm_%d_%d",i,j),100,0.7,1.4);
         hDistrAmpRatioZm[2*i+j]->SetXTitle("Amplitude ratio");
         hDistrAmpRatioZm[2*i+j]->SetYTitle("Counts");
         hDistrAmpRatioZm[2*i+j]->SetStats(0);
         hDistrAmpRatioZm[2*i+j]->SetLineColor(tabCol[2-j]);
         hDistrERatioZm[2*i+j] = new TH1F(Form("hDistrERatioZm_%d_%d",i,j),Form("hDistrERatioZm_%d_%d",i,j),100,0.7,1.4);
         hDistrERatioZm[2*i+j]->SetXTitle("Energy ratio");
         hDistrERatioZm[2*i+j]->SetYTitle("Counts");
         hDistrERatioZm[2*i+j]->SetStats(0);
         hDistrERatioZm[2*i+j]->SetLineColor(tabCol[2-j]);
         h2Amp[2*i+j] = new TH2F(Form("h2Amp%d_%d",i,j),Form("h2Amp%d_%d",i,j),100,0.,multFactGraphAmp*maxAmp,100,0.,multFactGraphAmp*maxAmp);
         h2Amp[2*i+j]->SetXTitle("Amplitude 1st run");
         h2Amp[2*i+j]->SetYTitle("Amplitude other run");
         h2Amp[2*i+j]->SetStats(0);
         h2Amp[2*i+j]->SetContour(30);
         h2E[2*i+j] = new TH2F(Form("h2E%d_%d",i,j),Form("h2E%d_%d",i,j),100,0.,multFactGraphE*maxE,100,0.,multFactGraphE*maxE);
         h2E[2*i+j]->SetXTitle("Energy 1st run");
         h2E[2*i+j]->SetYTitle("Energy other run");
         h2E[2*i+j]->SetStats(0);
         h2E[2*i+j]->SetContour(30);
         }
     }
 for (j=0;j<6;j++)
    {for (i=0;i<kNbSMtot+1;i++)
        {h2RatioVsIntg[6*i+j] = new TH2F(Form("h2RatioVsIntg%d_%d",i,j),Form("h2RatioVsIntg%d_%d",i,j),100,0.,maxIntg,100,0.6,1.5);
         h2RatioVsIntg[6*i+j]->SetXTitle("Integral");
         if (runMode == 0) h2RatioVsIntg[6*i+j]->SetYTitle("Energy ratio");
         if (runMode == 1) h2RatioVsIntg[6*i+j]->SetYTitle("Amplitude ratio");
         h2RatioVsIntg[6*i+j]->SetStats(0);
         h2RatioVsIntg[6*i+j]->SetContour(30);
         h2RatioVsAmp[6*i+j] = new TH2F(Form("h2RatioVsAmp%d_%d",i,j),Form("h2RatioVsAmp%d_%d",i,j),100,0.,multFactGraphAmp*maxAmp,100,0.6,1.5);
         h2RatioVsAmp[6*i+j]->SetXTitle("Amplitude");
         if (runMode == 0) h2RatioVsAmp[6*i+j]->SetYTitle("Energy ratio");
         if (runMode == 1) h2RatioVsAmp[6*i+j]->SetYTitle("Amplitude ratio");
         h2RatioVsAmp[6*i+j]->SetStats(0);
         h2RatioVsAmp[6*i+j]->SetContour(30);
         }
     }
 TH2F *hAllSpaceEMCALOCDBfactorsBefore = new TH2F("hAllSpaceEMCALOCDBfactorsBefore","hAllSpaceEMCALOCDBfactorsBefore",2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
 TH2F *hAllSpaceEMCALRatioAmpliOCDB = new TH2F("hAllSpaceEMCALRatioAmpliOCDB","hAllSpaceEMCALRatioAmpliOCDB",2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
 TH2F *hAllSpaceEMCALRatioEnergyOCDB = new TH2F("hAllSpaceEMCALRatioEnergyOCDB","hAllSpaceEMCALRatioEnergyOCDB",2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird,-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
 TH2F *hAllSpaceDCALOCDBfactorsBefore = new TH2F("hAllSpaceDCALOCDBfactorsBefore","hAllSpaceDCALOCDBfactorsBefore",2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
 TH2F *hAllSpaceDCALRatioAmpliOCDB = new TH2F("hAllSpaceDCALRatioAmpliOCDB","hAllSpaceDCALRatioAmpliOCDB",2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
 TH2F *hAllSpaceDCALRatioEnergyOCDB = new TH2F("hAllSpaceDCALRatioEnergyOCDB","hAllSpaceDCALRatioEnergyOCDB",2*kNbColMax,-0.5,2.*kNbColMax-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird,-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
 TH2F *h2CorrelEnergyOCDB = new TH2F("h2CorrelEnergyOCDB","h2CorrelEnergyOCDB",100,0.3,3.0,100,0.,4.);
 TH2F *h2CorrelEnergyOCDBZm = new TH2F("h2CorrelEnergyOCDBZm","h2CorrelEnergyOCDBZm",100,0.7,1.4,100,0.7,1.5);
 TH2F *h2CorrelEnergyOCDBZmZm = new TH2F("h2CorrelEnergyOCDBZmZm","h2CorrelEnergyOCDBZmZm",100,0.9,1.1,100,0.97,1.03);
 hAllSpaceEMCALOCDBfactorsBefore->SetXTitle("Column");
 hAllSpaceEMCALOCDBfactorsBefore->SetYTitle("Row");
 hAllSpaceEMCALOCDBfactorsBefore->SetStats(0);
 hAllSpaceEMCALOCDBfactorsBefore->SetContour(30);
 hAllSpaceEMCALRatioAmpliOCDB->SetXTitle("Column");
 hAllSpaceEMCALRatioAmpliOCDB->SetYTitle("Row");
 hAllSpaceEMCALRatioAmpliOCDB->SetStats(0);
 hAllSpaceEMCALRatioAmpliOCDB->SetContour(30);
 hAllSpaceEMCALRatioEnergyOCDB->SetXTitle("Column");
 hAllSpaceEMCALRatioEnergyOCDB->SetYTitle("Row");
 hAllSpaceEMCALRatioEnergyOCDB->SetStats(0);
 hAllSpaceEMCALRatioEnergyOCDB->SetContour(30);
 hAllSpaceDCALOCDBfactorsBefore->SetXTitle("Column");
 hAllSpaceDCALOCDBfactorsBefore->SetYTitle("Row");
 hAllSpaceDCALOCDBfactorsBefore->SetStats(0);
 hAllSpaceDCALOCDBfactorsBefore->SetContour(30);
 hAllSpaceDCALRatioAmpliOCDB->SetXTitle("Column");
 hAllSpaceDCALRatioAmpliOCDB->SetYTitle("Row");
 hAllSpaceDCALRatioAmpliOCDB->SetStats(0);
 hAllSpaceDCALRatioAmpliOCDB->SetContour(30);
 hAllSpaceDCALRatioEnergyOCDB->SetXTitle("Column");
 hAllSpaceDCALRatioEnergyOCDB->SetYTitle("Row");
 hAllSpaceDCALRatioEnergyOCDB->SetStats(0);
 hAllSpaceDCALRatioEnergyOCDB->SetContour(30);
 h2CorrelEnergyOCDB->SetXTitle("Old OCDB coeff");
 h2CorrelEnergyOCDB->SetYTitle("Energy");
 h2CorrelEnergyOCDB->SetStats(0);
 h2CorrelEnergyOCDB->SetContour(30);
 h2CorrelEnergyOCDBZm->SetXTitle("Old OCDB coeff");
 h2CorrelEnergyOCDBZm->SetYTitle("Energy");
 h2CorrelEnergyOCDBZm->SetStats(0);
 h2CorrelEnergyOCDBZm->SetContour(30);
 h2CorrelEnergyOCDBZmZm->SetXTitle("Old OCDB coeff");
 h2CorrelEnergyOCDBZmZm->SetYTitle("Energy");
 h2CorrelEnergyOCDBZmZm->SetStats(0);
 h2CorrelEnergyOCDBZmZm->SetContour(30);

 TH1F *h;
 for (i=0;i<3;i++)
    {cmptIntg[i]=0;
     cmptAmpl[i]=0;
     cmptEmismatch[i]=0;
     }
 printf("\nStart reading data. Tower coords given below are in Alice-offline mapping.\n");
 fprintf(outFile,"\nTower coords given below are in Alice-offline mapping.\n\n");
 for(iSM=0;iSM<kNbSMtot;iSM++)
    {printf("Doing SM %d (%s, %s)...\n",iSM,SMP2Name[iSM],detTypeString[SMdetType[iSM]]);
     for (iCol=0;iCol<kTabNbCol[SMdetType[iSM]];iCol++)
        {for (iRow=0;iRow<kTabNbRow[SMdetType[iSM]];iRow++)
            {iColGlob=(iSM%2)*(2*kNbColEMCAL-kTabNbCol[SMdetType[iSM]])+iCol;
             if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                {iRowGlob=(int)(iSM/2)*kNbRowEMCAL+iRow;
                 if (runMode == 0) hAllSpaceEMCALOCDBfactorsBefore->SetBinContent(iColGlob+1,iRowGlob+1,coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted);
                 }
                else
                {iRowGlob=(int)((iSM-kNbSMEMCAL-kNbSMEMCALthird)/2)*kNbRowDCAL+iRow;
                 if (runMode == 0) hAllSpaceDCALOCDBfactorsBefore->SetBinContent(iColGlob+1,iRowGlob+1,coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted);
                 }
             //Look at the integrals :
             for (i=0;i<3;i++)
                {hDistrIntg[3*iSM+i]->Fill(hampHighGainTow[i][iSM][iRow][iCol]->Integral());
                 hDistrIntg[3*kNbSMtot+i]->Fill(hampHighGainTow[i][iSM][iRow][iCol]->Integral());
                 if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                    hAllSpaceEMCALIntg[i]->SetBinContent(iColGlob+1,iRowGlob+1,hampHighGainTow[i][iSM][iRow][iCol]->Integral());
                    else
                    hAllSpaceDCALIntg[i]->SetBinContent(iColGlob+1,iRowGlob+1,hampHighGainTow[i][iSM][iRow][iCol]->Integral());
                 }
             if (hampHighGainTow[0][iSM][iRow][iCol]->Integral() <= cutMinIntg)
                {fprintf(outFile,"Too few counts 1 (%d) in tower (SM,col,row) (%d,%d,%d), skip it.\n",(int)hampHighGainTow[0][iSM][iRow][iCol]->Integral(),iSM,iCol,iRow);
                 cmptIntg[0]++;
                 continue;
                 }
             if (hampHighGainTow[1][iSM][iRow][iCol]->Integral() > cutMinIntg)
                {if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                    hAllSpaceEMCALIntgRatio[1]->SetBinContent(iColGlob+1,iRowGlob+1,hampHighGainTow[1][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                    else
                    hAllSpaceDCALIntgRatio[1]->SetBinContent(iColGlob+1,iRowGlob+1,hampHighGainTow[1][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 hDistrIntgRatio[2*iSM+1]->Fill(hampHighGainTow[1][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 hDistrIntgRatio[2*kNbSMtot+1]->Fill(hampHighGainTow[1][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 hDistrIntgRatioZm[2*iSM+1]->Fill(hampHighGainTow[1][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 hDistrIntgRatioZm[2*kNbSMtot+1]->Fill(hampHighGainTow[1][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 }
                else
                {fprintf(outFile,"Too few counts 2 (%d) in tower (SM,col,row) (%d,%d,%d) for run 2 while run 1 OK.\n",(int)hampHighGainTow[1][iSM][iRow][iCol]->Integral(),iSM,iCol,iRow);
                 cmptIntg[1]++; 
                 }
             if (hampHighGainTow[2][iSM][iRow][iCol]->Integral() > cutMinIntg)
                {if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                    hAllSpaceEMCALIntgRatio[0]->SetBinContent(iColGlob+1,iRowGlob+1,hampHighGainTow[2][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                    else
                    hAllSpaceDCALIntgRatio[0]->SetBinContent(iColGlob+1,iRowGlob+1,hampHighGainTow[2][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 hDistrIntgRatio[2*iSM+0]->Fill(hampHighGainTow[2][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 hDistrIntgRatio[2*kNbSMtot+0]->Fill(hampHighGainTow[2][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 hDistrIntgRatioZm[2*iSM+0]->Fill(hampHighGainTow[2][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 hDistrIntgRatioZm[2*kNbSMtot+0]->Fill(hampHighGainTow[2][iSM][iRow][iCol]->Integral()/hampHighGainTow[0][iSM][iRow][iCol]->Integral());
                 }
                else
                {fprintf(outFile,"Too few counts 3 (%d) in tower (SM,col,row) (%d,%d,%d) for run 3 while run 1 OK.\n",(int)hampHighGainTow[2][iSM][iRow][iCol]->Integral(),iSM,iCol,iRow);
                 cmptIntg[2]++; 
                 }
             //Look at the amplitudes :
             h=hampHighGainTow[0][iSM][iRow][iCol];
             maxloc=h->GetBinCenter(h->GetMaximumBin());
             h->SetAxisRange(maxloc-15.,maxloc+15.);
             h->SetAxisRange(maxloc-2*h->GetRMS(),maxloc+2*h->GetRMS());
             ampli0=h->GetMean();
             hDistrAmp[3*iSM+0]->Fill(ampli0);
             hDistrAmp[3*kNbSMtot+0]->Fill(ampli0);
             if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                hAllSpaceEMCALAmp[0]->SetBinContent(iColGlob+1,iRowGlob+1,ampli0);
                else
                hAllSpaceDCALAmp[0]->SetBinContent(iColGlob+1,iRowGlob+1,ampli0);
             if (runMode == 0) E0=ampli0*coeffOCDBbefore[iSM][iCol][iRow];
             hDistrE[3*iSM+0]->Fill(E0);
             hDistrE[3*kNbSMtot+0]->Fill(E0);
             if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                hAllSpaceEMCALE[0]->SetBinContent(iColGlob+1,iRowGlob+1,E0);
                else
                hAllSpaceDCALE[0]->SetBinContent(iColGlob+1,iRowGlob+1,E0);
             if (ampli0 <= cutMinAmp)
                {fprintf(outFile,"Too low amplitude 1 (%d) in tower (SM,col,row) (%d,%d,%d), skip it.\n",(int)ampli0,iSM,iCol,iRow);
                 cmptAmpl[0]++;
                 continue;
                 }
             //CUSTOMIZE customize
             // I don't remember what this does, leave it like this.
             /////if(h->GetRMS()==0.) continue;
             if (hampHighGainTow[1][iSM][iRow][iCol]->Integral() > cutMinIntg)
                {h=hampHighGainTow[1][iSM][iRow][iCol];
                 maxloc=h->GetBinCenter(h->GetMaximumBin());
                 h->SetAxisRange(maxloc-15.,maxloc+15.);
                 h->SetAxisRange(maxloc-2*h->GetRMS(),maxloc+2*h->GetRMS());
                 ampli1=h->GetMean();
                 hDistrAmp[3*iSM+1]->Fill(ampli1);
                 hDistrAmp[3*kNbSMtot+1]->Fill(ampli1);
                 if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                    hAllSpaceEMCALAmp[1]->SetBinContent(iColGlob+1,iRowGlob+1,ampli1);
                    else
                    hAllSpaceDCALAmp[1]->SetBinContent(iColGlob+1,iRowGlob+1,ampli1);
                 if (runMode == 0) E1=ampli1*coeffOCDBafter[iSM][iCol][iRow];
                 hDistrE[3*iSM+1]->Fill(E1);
                 hDistrE[3*kNbSMtot+1]->Fill(E1);
                 if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                    hAllSpaceEMCALE[1]->SetBinContent(iColGlob+1,iRowGlob+1,E1);
                    else
                    hAllSpaceDCALE[1]->SetBinContent(iColGlob+1,iRowGlob+1,E1);
                 if (ampli1 > cutMinAmp)
                    {if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                        {hAllSpaceEMCALAmpRatio[1]->SetBinContent(iColGlob+1,iRowGlob+1,ampli1/ampli0);
                         if (runMode == 0) hAllSpaceEMCALRatioAmpliOCDB->SetBinContent(iColGlob+1,iRowGlob+1,(ampli1/ampli0)/(coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted));
                         }
                        else
                        {hAllSpaceDCALAmpRatio[1]->SetBinContent(iColGlob+1,iRowGlob+1,ampli1/ampli0);
                         if (runMode == 0) hAllSpaceDCALRatioAmpliOCDB->SetBinContent(iColGlob+1,iRowGlob+1,(ampli1/ampli0)/(coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted));
                         }
                     hDistrAmpRatio[2*iSM+1]->Fill(ampli1/ampli0);
                     hDistrAmpRatio[2*kNbSMtot+1]->Fill(ampli1/ampli0);
                     hDistrAmpRatioZm[2*iSM+1]->Fill(ampli1/ampli0);
                     hDistrAmpRatioZm[2*kNbSMtot+1]->Fill(ampli1/ampli0);
                     if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                        {hAllSpaceEMCALERatio[1]->SetBinContent(iColGlob+1,iRowGlob+1,E1/E0);
                         if (runMode == 0) hAllSpaceEMCALRatioEnergyOCDB->SetBinContent(iColGlob+1,iRowGlob+1,(E1/E0)/(coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted));
                         }
                        else
                        {hAllSpaceDCALERatio[1]->SetBinContent(iColGlob+1,iRowGlob+1,E1/E0);
                         if (runMode == 0) hAllSpaceDCALRatioEnergyOCDB->SetBinContent(iColGlob+1,iRowGlob+1,(E1/E0)/(coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted));
                         }
                     if (runMode == 0)
                        {h2CorrelEnergyOCDB->Fill(coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted,E1/E0);
                         h2CorrelEnergyOCDBZm->Fill(coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted,E1/E0);
                         h2CorrelEnergyOCDBZmZm->Fill(coeffOCDBbefore[iSM][iCol][iRow]/coefFactorWanted,E1/E0);
                         }
                     hDistrERatio[2*iSM+1]->Fill(E1/E0);
                     hDistrERatio[2*kNbSMtot+1]->Fill(E1/E0);
                     hDistrERatioZm[2*iSM+1]->Fill(E1/E0);
                     hDistrERatioZm[2*kNbSMtot+1]->Fill(E1/E0);
                     if (runMode == 0)
                        {if (TMath::Abs(TMath::Log(E1)-TMath::Log(E0)) > TMath::Log(1.+cutWarnERatio))
                            {fprintf(outFile,"== E mismatch for tower (%d,%d,%d), ratio=%f  amp1=%f amp2=%f  OCDB1=%f OCDB2=%f\n",iSM,iCol,iRow,E1/E0,ampli0,ampli1,coeffOCDBbefore[iSM][iCol][iRow],coeffOCDBafter[iSM][iCol][iRow]);
                             cmptEmismatch[0]++;
                             }
                         }
                     if (runMode == 1)
                        {if (TMath::Abs(TMath::Log(ampli1)-TMath::Log(ampli0)) > TMath::Log(1.+cutWarnAmpRatio))
                            {fprintf(outFile,"== Amplitude mismatch for tower (%d,%d,%d), ratio=%f  amp0=%f amp1=%f\n",iSM,iCol,iRow,ampli1/ampli0,ampli0,ampli1);
                             cmptEmismatch[0]++;
                             }
                         }
                     h2Amp[2*iSM+1]->Fill(ampli0,ampli1);
                     h2Amp[2*kNbSMtot+1]->Fill(ampli0,ampli1);
                     h2E[2*iSM+1]->Fill(E0,E1);
                     h2E[2*kNbSMtot+1]->Fill(E0,E1);
                     }
                    else
                    {fprintf(outFile,"Too low amplitude 2 (%d) in tower (SM,col,row) (%d,%d,%d) for run 2 while run 1 OK.\n",(int)ampli1,iSM,iCol,iRow);
                     cmptAmpl[1]++;
                     }
                 }
             if (runMode == 0)
                {h2RatioVsIntg[6*iSM+0]->Fill(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),E1/E0);
                 h2RatioVsIntg[6*kNbSMtot+0]->Fill(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),E1/E0);
                 h2RatioVsIntg[6*iSM+2]->Fill(hampHighGainTow[1][iSM][iRow][iCol]->Integral(),E1/E0);
                 h2RatioVsIntg[6*kNbSMtot+2]->Fill(hampHighGainTow[1][iSM][iRow][iCol]->Integral(),E1/E0);
                 h2RatioVsIntg[6*iSM+4]->Fill(TMath::Min(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),hampHighGainTow[1][iSM][iRow][iCol]->Integral()),E1/E0);
                 h2RatioVsIntg[6*kNbSMtot+4]->Fill(TMath::Min(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),hampHighGainTow[1][iSM][iRow][iCol]->Integral()),E1/E0);
                 h2RatioVsAmp[6*iSM+0]->Fill(ampli0,E1/E0);
                 h2RatioVsAmp[6*kNbSMtot+0]->Fill(ampli0,E1/E0);
                 h2RatioVsAmp[6*iSM+2]->Fill(ampli1,E1/E0);
                 h2RatioVsAmp[6*kNbSMtot+2]->Fill(ampli1,E1/E0);
                 h2RatioVsAmp[6*iSM+4]->Fill(TMath::Min(ampli0,ampli1),E1/E0);
                 h2RatioVsAmp[6*kNbSMtot+4]->Fill(TMath::Min(ampli0,ampli1),E1/E0);
                 }
             if (runMode == 1)
                {h2RatioVsIntg[6*iSM+0]->Fill(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),ampli1/ampli0);
                 h2RatioVsIntg[6*kNbSMtot+0]->Fill(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),ampli1/ampli0);
                 h2RatioVsIntg[6*iSM+2]->Fill(hampHighGainTow[1][iSM][iRow][iCol]->Integral(),ampli1/ampli0);
                 h2RatioVsIntg[6*kNbSMtot+2]->Fill(hampHighGainTow[1][iSM][iRow][iCol]->Integral(),ampli1/ampli0);
                 h2RatioVsIntg[6*iSM+4]->Fill(TMath::Min(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),hampHighGainTow[1][iSM][iRow][iCol]->Integral()),ampli1/ampli0);
                 h2RatioVsIntg[6*kNbSMtot+4]->Fill(TMath::Min(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),hampHighGainTow[1][iSM][iRow][iCol]->Integral()),ampli1/ampli0);
                 h2RatioVsAmp[6*iSM+0]->Fill(ampli0,ampli1/ampli0);
                 h2RatioVsAmp[6*kNbSMtot+0]->Fill(ampli0,ampli1/ampli0);
                 h2RatioVsAmp[6*iSM+2]->Fill(ampli1,ampli1/ampli0);
                 h2RatioVsAmp[6*kNbSMtot+2]->Fill(ampli1,ampli1/ampli0);
                 h2RatioVsAmp[6*iSM+4]->Fill(TMath::Min(ampli0,ampli1),ampli1/ampli0);
                 h2RatioVsAmp[6*kNbSMtot+4]->Fill(TMath::Min(ampli0,ampli1),ampli1/ampli0);
                 }
             if (hampHighGainTow[2][iSM][iRow][iCol]->Integral() > cutMinIntg)
                {h=hampHighGainTow[2][iSM][iRow][iCol];
                 maxloc=h->GetBinCenter(h->GetMaximumBin());
                 h->SetAxisRange(maxloc-15.,maxloc+15.);
                 h->SetAxisRange(maxloc-2*h->GetRMS(),maxloc+2*h->GetRMS());
                 ampli2=h->GetMean();
                 hDistrAmp[3*iSM+2]->Fill(ampli2);
                 hDistrAmp[3*kNbSMtot+2]->Fill(ampli2);
                 if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                    hAllSpaceEMCALAmp[2]->SetBinContent(iColGlob+1,iRowGlob+1,ampli2);
                    else
                    hAllSpaceDCALAmp[2]->SetBinContent(iColGlob+1,iRowGlob+1,ampli2);
                 if (runMode == 0) E2=ampli2*coeffOCDBbefore[iSM][iCol][iRow];
                 hDistrE[3*iSM+2]->Fill(E2);
                 hDistrE[3*kNbSMtot+2]->Fill(E2);
                 if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                    hAllSpaceEMCALE[2]->SetBinContent(iColGlob+1,iRowGlob+1,E2);
                    else
                    hAllSpaceDCALE[2]->SetBinContent(iColGlob+1,iRowGlob+1,E2);
                 if (ampli2 > cutMinAmp)
                    {if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                        hAllSpaceEMCALAmpRatio[0]->SetBinContent(iColGlob+1,iRowGlob+1,ampli2/ampli0);
                        else
                        hAllSpaceDCALAmpRatio[0]->SetBinContent(iColGlob+1,iRowGlob+1,ampli2/ampli0);
                     hDistrAmpRatio[2*iSM+0]->Fill(ampli2/ampli0);
                     hDistrAmpRatio[2*kNbSMtot+0]->Fill(ampli2/ampli0);
                     hDistrAmpRatioZm[2*iSM+0]->Fill(ampli2/ampli0);
                     hDistrAmpRatioZm[2*kNbSMtot+0]->Fill(ampli2/ampli0);
                     if ((SMdetType[iSM] == kEMCAL) || (SMdetType[iSM] == kEMCALthird))
                        hAllSpaceEMCALERatio[0]->SetBinContent(iColGlob+1,iRowGlob+1,E2/E0);
                        else
                        hAllSpaceDCALERatio[0]->SetBinContent(iColGlob+1,iRowGlob+1,E2/E0);
                     hDistrERatio[2*iSM+0]->Fill(E2/E0);
                     hDistrERatio[2*kNbSMtot+0]->Fill(E2/E0);
                     hDistrERatioZm[2*iSM+0]->Fill(E2/E0);
                     hDistrERatioZm[2*kNbSMtot+0]->Fill(E2/E0);
                     if (runMode == 1)
                        {if (TMath::Abs(TMath::Log(ampli2)-TMath::Log(ampli0)) > TMath::Log(1.+cutWarnAmpRatio)) fprintf(outFile,"== Amplitude mismatch for tower (%d,%d,%d), ratio=%f  amp0=%f amp2=%f\n",iSM,iCol,iRow,ampli2/ampli0,ampli0,ampli2);
                         }
                     h2Amp[2*iSM+0]->Fill(ampli0,ampli2);
                     h2Amp[2*kNbSMtot+0]->Fill(ampli0,ampli2);
                     h2E[2*iSM+0]->Fill(E0,E2);
                     h2E[2*kNbSMtot+0]->Fill(E0,E2);
                     }
                    else
                    {fprintf(outFile,"Too low amplitude 3 (%d) in tower (SM,col,row) (%d,%d,%d) for run 3 while run 1 OK.\n",(int)ampli2,iSM,iCol,iRow);
                     cmptAmpl[2]++;
                     }
                 }
             if (runMode == 0)
                {h2RatioVsIntg[6*iSM+1]->Fill(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),E2/E0);
                 h2RatioVsIntg[6*kNbSMtot+1]->Fill(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),E2/E0);
                 h2RatioVsIntg[6*iSM+3]->Fill(hampHighGainTow[2][iSM][iRow][iCol]->Integral(),E2/E0);
                 h2RatioVsIntg[6*kNbSMtot+3]->Fill(hampHighGainTow[2][iSM][iRow][iCol]->Integral(),E2/E0);
                 h2RatioVsIntg[6*iSM+5]->Fill(TMath::Min(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),hampHighGainTow[2][iSM][iRow][iCol]->Integral()),E2/E0);
                 h2RatioVsIntg[6*kNbSMtot+5]->Fill(TMath::Min(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),hampHighGainTow[2][iSM][iRow][iCol]->Integral()),E2/E0);
                 h2RatioVsAmp[6*iSM+1]->Fill(ampli0,E2/E0);
                 h2RatioVsAmp[6*kNbSMtot+1]->Fill(ampli0,E2/E0);
                 //h2RatioVsAmp[6*iSM+3]->Fill(ampli1,E2/E0);
                 //h2RatioVsAmp[6*kNbSMtot+3]->Fill(ampli1,E2/E0);
                 //Julien, bug corrected June 18, 2015.
                 h2RatioVsAmp[6*iSM+3]->Fill(ampli2,E2/E0);
                 h2RatioVsAmp[6*kNbSMtot+3]->Fill(ampli2,E2/E0);
                 h2RatioVsAmp[6*iSM+5]->Fill(TMath::Min(ampli0,ampli2),E2/E0);
                 h2RatioVsAmp[6*kNbSMtot+5]->Fill(TMath::Min(ampli0,ampli2),E2/E0);
                 }
             if (runMode == 1)
                {h2RatioVsIntg[6*iSM+1]->Fill(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),ampli2/ampli0);
                 h2RatioVsIntg[6*kNbSMtot+1]->Fill(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),ampli2/ampli0);
                 h2RatioVsIntg[6*iSM+3]->Fill(hampHighGainTow[2][iSM][iRow][iCol]->Integral(),ampli2/ampli0);
                 h2RatioVsIntg[6*kNbSMtot+3]->Fill(hampHighGainTow[2][iSM][iRow][iCol]->Integral(),ampli2/ampli0);
                 h2RatioVsIntg[6*iSM+5]->Fill(TMath::Min(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),hampHighGainTow[2][iSM][iRow][iCol]->Integral()),ampli2/ampli0);
                 h2RatioVsIntg[6*kNbSMtot+5]->Fill(TMath::Min(hampHighGainTow[0][iSM][iRow][iCol]->Integral(),hampHighGainTow[2][iSM][iRow][iCol]->Integral()),ampli2/ampli0);
                 h2RatioVsAmp[6*iSM+1]->Fill(ampli0,ampli2/ampli0);
                 h2RatioVsAmp[6*kNbSMtot+1]->Fill(ampli0,ampli2/ampli0);
                 //h2RatioVsAmp[6*iSM+3]->Fill(ampli1,ampli2/ampli0);
                 //h2RatioVsAmp[6*kNbSMtot+3]->Fill(ampli1,ampli2/ampli0);
                 //Julien, bug corrected June 18, 2015.
                 h2RatioVsAmp[6*iSM+3]->Fill(ampli2,ampli2/ampli0);
                 h2RatioVsAmp[6*kNbSMtot+3]->Fill(ampli2,ampli2/ampli0);
                 h2RatioVsAmp[6*iSM+5]->Fill(TMath::Min(ampli0,ampli2),ampli2/ampli0);
                 h2RatioVsAmp[6*kNbSMtot+5]->Fill(TMath::Min(ampli0,ampli2),ampli2/ampli0);
                 }
             }
         }
     }

 //rootFileIn->Close();
 TFile *rootFileOut = new TFile(fchNameRoot,"RECREATE");
 //Must be put once rootFileIn is closed, otherwise tries to write in rootFileIn (and fails because not permitted). Not sure why.

 
 fprintf(outFile,"\nTotal towers with too few counts in run 1 : %d (cut at %f)\n",cmptIntg[0],cutMinIntg);
 fprintf(outFile,"   Among towers OK : %d towers nonetheless have too few counts in run 2,\n                     %d towers nonetheless have too few counts in run 3.\n",cmptIntg[1],cmptIntg[2]);
 fprintf(outFile,"Total towers with too low amplitude in run 1 : %d (cut at %f)\n",cmptAmpl[0],cutMinAmp);
 fprintf(outFile,"   Among towers OK : %d towers nonetheless have too low amp in run 2,\n                     %d towers nonetheless have too low amp in run 3.\n",cmptAmpl[1],cmptAmpl[2]);
 if (runMode == 0) fprintf(outFile,"\nTotal towers with energy mismatch : %d (cut at %f %%)\n\n",cmptEmismatch[0],100*cutWarnERatio);
 if (runMode == 1)  fprintf(outFile,"\nTotal towers with amplitude mismatch : %d (cut at %f %%)\n\n",cmptEmismatch[0],100*cutWarnAmpRatio);


 printf("Draw histograms.\n");
 
 TLine *lineSMborderVEMCAL = new TLine(kNbColEMCAL-0.5,-0.5,kNbColEMCAL-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
 TLine *lineSMborderVDCALthird = new TLine(kNbColDCALthird-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5,kNbColDCALthird-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
 TLine *lineSMborderVDCAL1 = new TLine(kNbColDCAL-0.5,-0.5,kNbColDCAL-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5);
 TLine *lineSMborderVDCAL2 = new TLine(2*kNbColEMCAL-kNbColDCAL-0.5,-0.5,2*kNbColEMCAL-kNbColDCAL-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5);
 TLine **lineSMborderHEMCAL,**lineSMborderHDCAL;
 lineSMborderHEMCAL = new TLine*[(int)((kNbSMEMCAL+1)/2)];
 lineSMborderHDCAL = new TLine*[(int)((kNbSMDCAL+1)/2)];
 for (i=0;i<(int)((kNbSMEMCAL+1)/2);i++) lineSMborderHEMCAL[i] = new TLine(-0.5,(i+1)*kNbRowEMCAL-0.5,2.*kNbColEMCAL-0.5,(i+1)*kNbRowEMCAL-0.5);
 for (i=0;i<(int)((kNbSMDCAL+1)/2);i++) lineSMborderHDCAL[i] = new TLine(-0.5,(i+1)*kNbRowDCAL-0.5,2.*kNbColDCALthird-0.5,(i+1)*kNbRowDCAL-0.5);

 //Integrals, all SM and spatial :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,4);
 c1->cd(1);
 hDistrIntg[3*kNbSMtot+0]->Draw();
 hDistrIntg[3*kNbSMtot+1]->Draw("SAME");
 hDistrIntg[3*kNbSMtot+2]->Draw("SAME");
 gPad->SetLogy();
 for (i=0;i<3;i++) hDistrIntg[3*kNbSMtot+i]->Write();
 for (j=0;j<3;j++)
    {c1->cd(2*j+3);
     hAllSpaceEMCALIntg[j]->SetMaximum(maxIntg);
     hAllSpaceEMCALIntg[j]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     hAllSpaceEMCALIntg[j]->Write();
     }
 for (j=0;j<3;j++)
    {c1->cd(2*j+4);
     hAllSpaceDCALIntg[j]->SetMaximum(maxIntg);
     hAllSpaceDCALIntg[j]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     hAllSpaceDCALIntg[j]->Write();
     }
 c1->Update();

 //Integrals, for each SM :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=0;j<kNbSMEMCAL;j++)
    {c1->cd(j+1);
     hDistrIntg[3*j+0]->Draw();
     hDistrIntg[3*j+1]->Draw("SAME");
     hDistrIntg[3*j+2]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<3;i++) hDistrIntg[3*j+i]->Write();
     }
 c1->Update();
 
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=kNbSMEMCAL;j<kNbSMtot;j++)
    {c1->cd(j-kNbSMEMCAL+1);
     hDistrIntg[3*j+0]->Draw();
     hDistrIntg[3*j+1]->Draw("SAME");
     hDistrIntg[3*j+2]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<3;i++) hDistrIntg[3*j+i]->Write();
     }
 c1->Update();

 //Integral ratios, all SM and spatial :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(1);
 hDistrIntgRatio[2*kNbSMtot+0]->SetMinimum(0.6);
 hDistrIntgRatio[2*kNbSMtot+0]->Draw();
 hDistrIntgRatio[2*kNbSMtot+1]->Draw("SAME");
 gPad->SetLogy();
 for (i=0;i<2;i++) hDistrIntgRatio[2*kNbSMtot+i]->Write();
 c1->cd(2);
 hDistrIntgRatioZm[2*kNbSMtot+0]->SetMinimum(0.6);
 hDistrIntgRatioZm[2*kNbSMtot+0]->Draw();
 hDistrIntgRatioZm[2*kNbSMtot+1]->Draw("SAME");
 gPad->SetLogy();
 for (i=0;i<2;i++) hDistrIntgRatioZm[2*kNbSMtot+i]->Write();
 for (i=0;i<2;i++)
    {hAllSpaceEMCALIntgRatio[i]->SetMinimum(0.5);
     hAllSpaceEMCALIntgRatio[i]->SetMaximum(2.0);
     }
 c1->cd(5);
 hAllSpaceEMCALIntgRatio[0]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALIntgRatio[0]->Write();
 c1->cd(3);
 hAllSpaceEMCALIntgRatio[1]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALIntgRatio[1]->Write();
 c1->Update();
 for (i=0;i<2;i++)
    {hAllSpaceEMCALIntgRatio[i]->SetMinimum(0.9);
     hAllSpaceEMCALIntgRatio[i]->SetMaximum(1.1);
     }
 c1->cd(6);
 hAllSpaceEMCALIntgRatio[0]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALIntgRatio[0]->Write();
 c1->cd(4);
 hAllSpaceEMCALIntgRatio[1]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALIntgRatio[1]->Write();
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 for (i=0;i<2;i++)
    {hAllSpaceDCALIntgRatio[i]->SetMinimum(0.5);
     hAllSpaceDCALIntgRatio[i]->SetMaximum(2.0);
     }
 c1->cd(5);
 hAllSpaceDCALIntgRatio[0]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALIntgRatio[0]->Write();
 c1->cd(3);
 hAllSpaceDCALIntgRatio[1]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALIntgRatio[1]->Write();
 c1->Update();
 for (i=0;i<2;i++)
    {hAllSpaceDCALIntgRatio[i]->SetMinimum(0.9);
     hAllSpaceDCALIntgRatio[i]->SetMaximum(1.1);
     }
 c1->cd(6);
 hAllSpaceDCALIntgRatio[0]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALIntgRatio[0]->Write();
 c1->cd(4);
 hAllSpaceDCALIntgRatio[1]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALIntgRatio[1]->Write();
 c1->Update();

 //Integral ratios, for each SM :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=0;j<kNbSMEMCAL;j++)
    {c1->cd(j+1);
     hDistrIntgRatio[2*j+0]->Draw();
     hDistrIntgRatio[2*j+1]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<2;i++)
        {hDistrIntgRatio[2*j+i]->Write();
         hDistrIntgRatioZm[2*j+i]->Write();
         }
     }
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=kNbSMEMCAL;j<kNbSMtot;j++)
    {c1->cd(j-kNbSMEMCAL+1);
     hDistrIntgRatio[2*j+0]->Draw();
     hDistrIntgRatio[2*j+1]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<2;i++)
        {hDistrIntgRatio[2*j+i]->Write();
         hDistrIntgRatioZm[2*j+i]->Write();
         }
     }
 c1->Update();

 //Amplitudes, all SM and spatial :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,4);
 c1->cd(1);
 hDistrAmp[3*kNbSMtot+0]->Draw();
 hDistrAmp[3*kNbSMtot+1]->Draw("SAME");
 hDistrAmp[3*kNbSMtot+2]->Draw("SAME");
 gPad->SetLogy();
 for (i=0;i<3;i++) hDistrAmp[3*kNbSMtot+i]->Write();
 for (j=0;j<3;j++)
    {c1->cd(2*j+3);
     hAllSpaceEMCALAmp[j]->SetMaximum(maxAmp);
     hAllSpaceEMCALAmp[j]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     hAllSpaceEMCALAmp[j]->Write();
     }
 for (j=0;j<3;j++)
    {c1->cd(2*j+4);
     hAllSpaceDCALAmp[j]->SetMaximum(maxAmp);
     hAllSpaceDCALAmp[j]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     hAllSpaceDCALAmp[j]->Write();
     }
 c1->Update();

 //Amplitudes, for each SM :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=0;j<kNbSMEMCAL;j++)
    {c1->cd(j+1);
     hDistrAmp[3*j+0]->Draw();
     hDistrAmp[3*j+1]->Draw("SAME");
     hDistrAmp[3*j+2]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<3;i++) hDistrAmp[3*j+i]->Write();
     }
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=kNbSMEMCAL;j<kNbSMtot;j++)
    {c1->cd(j-kNbSMEMCAL+1);
     hDistrAmp[3*j+0]->Draw();
     hDistrAmp[3*j+1]->Draw("SAME");
     hDistrAmp[3*j+2]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<3;i++) hDistrAmp[3*j+i]->Write();
     }
 c1->Update();

 //Amplitude ratios, all SM and spatial :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(1);
 hDistrAmpRatio[2*kNbSMtot+0]->SetMinimum(0.6);
 hDistrAmpRatio[2*kNbSMtot+0]->Draw();
 hDistrAmpRatio[2*kNbSMtot+1]->Draw("SAME");
 gPad->SetLogy();
 for (i=0;i<2;i++) hDistrAmpRatio[2*kNbSMtot+i]->Write();
 c1->cd(2);
 hDistrAmpRatioZm[2*kNbSMtot+0]->SetMinimum(0.6);
 hDistrAmpRatioZm[2*kNbSMtot+0]->Draw();
 hDistrAmpRatioZm[2*kNbSMtot+1]->Draw("SAME");
 gPad->SetLogy();
 for (i=0;i<2;i++) hDistrAmpRatioZm[2*kNbSMtot+i]->Write();
 for (i=0;i<2;i++)
    {hAllSpaceEMCALAmpRatio[i]->SetMinimum(0.5);
     hAllSpaceEMCALAmpRatio[i]->SetMaximum(2.0);
     }
 c1->cd(5);
 hAllSpaceEMCALAmpRatio[0]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALAmpRatio[0]->Write();
 c1->cd(3);
 hAllSpaceEMCALAmpRatio[1]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALAmpRatio[1]->Write();
 c1->Update();
 for (i=0;i<2;i++)
    {hAllSpaceEMCALAmpRatio[i]->SetMinimum(0.9);
     hAllSpaceEMCALAmpRatio[i]->SetMaximum(1.1);
     }
 c1->cd(6);
 hAllSpaceEMCALAmpRatio[0]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALAmpRatio[0]->Write();
 c1->cd(4);
 hAllSpaceEMCALAmpRatio[1]->Draw("COLZ");
 lineSMborderVEMCAL->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
 hAllSpaceEMCALAmpRatio[1]->Write();
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 for (i=0;i<2;i++)
    {hAllSpaceDCALAmpRatio[i]->SetMinimum(0.5);
     hAllSpaceDCALAmpRatio[i]->SetMaximum(2.0);
     }
 c1->cd(5);
 hAllSpaceDCALAmpRatio[0]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALAmpRatio[0]->Write();
 c1->cd(3);
 hAllSpaceDCALAmpRatio[1]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALAmpRatio[1]->Write();
 c1->Update();
 for (i=0;i<2;i++)
    {hAllSpaceDCALAmpRatio[i]->SetMinimum(0.9);
     hAllSpaceDCALAmpRatio[i]->SetMaximum(1.1);
     }
 c1->cd(6);
 hAllSpaceDCALAmpRatio[0]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALAmpRatio[0]->Write();
 c1->cd(4);
 hAllSpaceDCALAmpRatio[1]->Draw("COLZ");
 lineSMborderVDCALthird->Draw("SAME");
 lineSMborderVDCAL1->Draw("SAME");
 lineSMborderVDCAL2->Draw("SAME");
 for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
 hAllSpaceDCALAmpRatio[1]->Write();
 c1->Update();
 
 //Amplitude ratios, for each SM :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=0;j<kNbSMEMCAL;j++)
    {c1->cd(j+1);
     hDistrAmpRatio[2*j+0]->Draw();
     hDistrAmpRatio[2*j+1]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<3;i++) hDistrAmpRatio[2*j+i]->Write();
     for (i=0;i<2;i++) hDistrAmpRatioZm[2*j+i]->Write();
     }
 c1->Update();

 ps->NewPage();
 c1->Clear();
 c1->Divide(2,5);
 for (j=kNbSMEMCAL;j<kNbSMtot;j++)
    {c1->cd(j-kNbSMEMCAL+1);
     hDistrAmpRatio[2*j+0]->Draw();
     hDistrAmpRatio[2*j+1]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<3;i++) hDistrAmpRatio[2*j+i]->Write();
     for (i=0;i<2;i++) hDistrAmpRatioZm[2*j+i]->Write();
     }
 c1->Update();

 if (runMode == 0)
    {//Energies, all SM and spatial :
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,4);
     c1->cd(1);
     hDistrE[3*kNbSMtot+0]->Draw();
     hDistrE[3*kNbSMtot+1]->Draw("SAME");
     hDistrE[3*kNbSMtot+2]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<3;i++) hDistrE[3*kNbSMtot+i]->Write();
     for (j=0;j<3;j++)
        {c1->cd(2*j+3);
         hAllSpaceEMCALE[j]->SetMaximum(maxE);
         hAllSpaceEMCALE[j]->Draw("COLZ");
         lineSMborderVEMCAL->Draw("SAME");
         for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
         hAllSpaceEMCALE[j]->Write();
         }
     for (j=0;j<3;j++)
        {c1->cd(2*j+4);
         hAllSpaceDCALE[j]->SetMaximum(maxE);
         hAllSpaceDCALE[j]->Draw("COLZ");
         lineSMborderVDCALthird->Draw("SAME");
         lineSMborderVDCAL1->Draw("SAME");
         lineSMborderVDCAL2->Draw("SAME");
         for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
         hAllSpaceDCALE[j]->Write();
         }
     c1->Update();
    
     //Energies, for each SM :
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (j=0;j<kNbSMEMCAL;j++)
        {c1->cd(j+1);
         hDistrE[3*j+0]->Draw();
         hDistrE[3*j+1]->Draw("SAME");
         hDistrE[3*j+2]->Draw("SAME");
         gPad->SetLogy();
         for (i=0;i<3;i++) hDistrE[3*j+i]->Write();
         }
     c1->Update();
    
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (j=kNbSMEMCAL;j<kNbSMtot;j++)
        {c1->cd(j-kNbSMEMCAL+1);
         hDistrE[3*j+0]->Draw();
         hDistrE[3*j+1]->Draw("SAME");
         hDistrE[3*j+2]->Draw("SAME");
         gPad->SetLogy();
         for (i=0;i<3;i++) hDistrE[3*j+i]->Write();
         }
     c1->Update();
    
     //Energy ratios, all SM and spatial :
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     hDistrERatio[2*kNbSMtot+0]->SetMinimum(0.6);
     hDistrERatio[2*kNbSMtot+0]->Draw();
     hDistrERatio[2*kNbSMtot+1]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<2;i++) hDistrERatio[2*kNbSMtot+i]->Write();
     c1->cd(2);
     hDistrERatioZm[2*kNbSMtot+0]->SetMinimum(0.6);
     hDistrERatioZm[2*kNbSMtot+0]->Draw();
     hDistrERatioZm[2*kNbSMtot+1]->Draw("SAME");
     gPad->SetLogy();
     for (i=0;i<2;i++) hDistrERatioZm[2*kNbSMtot+i]->Write();
     for (i=0;i<2;i++)
        {hAllSpaceEMCALERatio[i]->SetMinimum(0.5);
         hAllSpaceEMCALERatio[i]->SetMaximum(2.0);
         }
     c1->cd(5);
     hAllSpaceEMCALERatio[0]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     hAllSpaceEMCALERatio[0]->Write();
     c1->cd(3);
     hAllSpaceEMCALERatio[1]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     hAllSpaceEMCALERatio[1]->Write();
     c1->Update();
     for (i=0;i<2;i++)
        {hAllSpaceEMCALERatio[i]->SetMinimum(0.9);
         hAllSpaceEMCALERatio[i]->SetMaximum(1.1);
         }
     c1->cd(6);
     hAllSpaceEMCALERatio[0]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     hAllSpaceEMCALERatio[0]->Write();
     c1->cd(4);
     hAllSpaceEMCALERatio[1]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     hAllSpaceEMCALERatio[1]->Write();
     c1->Update();
    
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     for (i=0;i<2;i++)
        {hAllSpaceDCALERatio[i]->SetMinimum(0.5);
         hAllSpaceDCALERatio[i]->SetMaximum(2.0);
         }
     c1->cd(5);
     hAllSpaceDCALERatio[0]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     hAllSpaceDCALERatio[0]->Write();
     c1->cd(3);
     hAllSpaceDCALERatio[1]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     hAllSpaceDCALERatio[1]->Write();
     c1->Update();
     for (i=0;i<2;i++)
        {hAllSpaceDCALERatio[i]->SetMinimum(0.9);
         hAllSpaceDCALERatio[i]->SetMaximum(1.1);
         }
     c1->cd(6);
     hAllSpaceDCALERatio[0]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     hAllSpaceDCALERatio[0]->Write();
     c1->cd(4);
     hAllSpaceDCALERatio[1]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     hAllSpaceDCALERatio[1]->Write();
     c1->Update();
    
     //Energy ratios, for each SM :
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (j=0;j<kNbSMEMCAL;j++)
        {c1->cd(j+1);
         hDistrERatio[2*j+0]->Draw();
         hDistrERatio[2*j+1]->Draw("SAME");
         gPad->SetLogy();
         for (i=0;i<3;i++) hDistrERatio[2*j+i]->Write();
         for (i=0;i<2;i++) hDistrERatioZm[2*j+i]->Write();
         }
     c1->Update();
     
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,5);
     for (j=kNbSMEMCAL;j<kNbSMtot;j++)
        {c1->cd(j-kNbSMEMCAL+1);
         hDistrERatio[2*j+0]->Draw();
         hDistrERatio[2*j+1]->Draw("SAME");
         gPad->SetLogy();
         for (i=0;i<3;i++) hDistrERatio[2*j+i]->Write();
         for (i=0;i<2;i++) hDistrERatioZm[2*j+i]->Write();
         }
     c1->Update();
     }

 //Relationship between OCDB factors (requested gain change) and amplitude ratio, energy ratio : 
 if (runMode == 0)
    {TLine *ligne1H = new TLine(0.,1.,1.,1.);
     TF1 *fInv = new TF1("fInv","1/x",0.,4.); //No correction (?)
     TF1 *fInv2 = new TF1("fInv2","1/(x*x)",0.,4.); //Swapped correction
     fInv->SetLineWidth(0.5);
     fInv->SetLineColor(1);
     fInv2->SetLineWidth(0.5);
     fInv2->SetLineColor(1);
     
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,4);
     c1->cd(1);
     hAllSpaceEMCALAmpRatio[1]->SetMinimum(0.8);
     hAllSpaceEMCALAmpRatio[1]->SetMaximum(1.2);
     hAllSpaceEMCALAmpRatio[1]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->cd(3);
     hAllSpaceEMCALERatio[1]->SetMinimum(0.8);
     hAllSpaceEMCALERatio[1]->SetMaximum(1.2);
     hAllSpaceEMCALERatio[1]->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->cd(5);
     hAllSpaceEMCALOCDBfactorsBefore->SetMinimum(0.25);
     hAllSpaceEMCALOCDBfactorsBefore->SetMaximum(1.8);
     hAllSpaceEMCALOCDBfactorsBefore->Draw("COLZ");
     hAllSpaceEMCALOCDBfactorsBefore->Write();
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     c1->cd(6);
     hAllSpaceEMCALOCDBfactorsBefore->SetMinimum(0.8);
     hAllSpaceEMCALOCDBfactorsBefore->SetMaximum(1.2);
     hAllSpaceEMCALOCDBfactorsBefore->Draw("COLZ");
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->cd(2);
     hAllSpaceEMCALRatioAmpliOCDB->SetMinimum(0.7);
     hAllSpaceEMCALRatioAmpliOCDB->SetMaximum(1.3);
     hAllSpaceEMCALRatioAmpliOCDB->Draw("COLZ");
     hAllSpaceEMCALRatioAmpliOCDB->Write();
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->cd(4);
     hAllSpaceEMCALRatioEnergyOCDB->SetMinimum(0.6);
     hAllSpaceEMCALRatioEnergyOCDB->SetMaximum(1.8);
     hAllSpaceEMCALRatioEnergyOCDB->Draw("COLZ");
     hAllSpaceEMCALRatioEnergyOCDB->Write();
     lineSMborderVEMCAL->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMEMCAL+1)/2);jm++) lineSMborderHEMCAL[jm]->Draw("SAME");
     c1->Update();
     
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,4);
     c1->cd(1);
     hAllSpaceDCALAmpRatio[1]->SetMinimum(0.8);
     hAllSpaceDCALAmpRatio[1]->SetMaximum(1.2);
     hAllSpaceDCALAmpRatio[1]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->cd(3);
     hAllSpaceDCALERatio[1]->SetMinimum(0.8);
     hAllSpaceDCALERatio[1]->SetMaximum(1.2);
     hAllSpaceDCALERatio[1]->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->cd(5);
     hAllSpaceDCALOCDBfactorsBefore->SetMinimum(0.25);
     hAllSpaceDCALOCDBfactorsBefore->SetMaximum(1.8);
     hAllSpaceDCALOCDBfactorsBefore->Draw("COLZ");
     hAllSpaceDCALOCDBfactorsBefore->Write();
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     c1->cd(6);
     hAllSpaceDCALOCDBfactorsBefore->SetMinimum(0.8);
     hAllSpaceDCALOCDBfactorsBefore->SetMaximum(1.2);
     hAllSpaceDCALOCDBfactorsBefore->Draw("COLZ");
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->cd(2);
     hAllSpaceDCALRatioAmpliOCDB->SetMinimum(0.7);
     hAllSpaceDCALRatioAmpliOCDB->SetMaximum(1.3);
     hAllSpaceDCALRatioAmpliOCDB->Draw("COLZ");
     hAllSpaceDCALRatioAmpliOCDB->Write();
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->cd(4);
     hAllSpaceDCALRatioEnergyOCDB->SetMinimum(0.6);
     hAllSpaceDCALRatioEnergyOCDB->SetMaximum(1.8);
     hAllSpaceDCALRatioEnergyOCDB->Draw("COLZ");
     hAllSpaceDCALRatioEnergyOCDB->Write();
     lineSMborderVDCALthird->Draw("SAME");
     lineSMborderVDCAL1->Draw("SAME");
     lineSMborderVDCAL2->Draw("SAME");
     for (int jm=0;jm<(int)((kNbSMDCAL+1)/2);jm++) lineSMborderHDCAL[jm]->Draw("SAME");
     c1->Update();
     
     ps->NewPage();
     c1->Clear();
     c1->Divide(2,3);
     c1->cd(1);
     h2CorrelEnergyOCDB->Draw("COLZ");
     h2CorrelEnergyOCDB->Write();
     gPad->SetLogz();
     c1->Update();
     c1->cd(3);
     h2CorrelEnergyOCDBZm->Draw("COLZ");
     h2CorrelEnergyOCDBZm->Write();
     gPad->SetLogz();
     c1->cd(2);
     ligne1H->SetX1(h2CorrelEnergyOCDB->GetXaxis()->GetXmin());
     ligne1H->SetX2(h2CorrelEnergyOCDB->GetXaxis()->GetXmax());
     h2CorrelEnergyOCDB->Draw("COLZ");
     ligne1H->Draw("SAME");
     fInv->Draw("SAME");
     fInv2->Draw("SAME");
     c1->Update();
     c1->cd(4);
     ligne1H->SetX1(h2CorrelEnergyOCDBZm->GetXaxis()->GetXmin());
     ligne1H->SetX2(h2CorrelEnergyOCDBZm->GetXaxis()->GetXmax());
     h2CorrelEnergyOCDBZm->Draw("COLZ");
     ligne1H->Draw("SAME");
     fInv->Draw("SAME");
     fInv2->Draw("SAME");
     c1->Update();
     c1->cd(6);
     ligne1H->SetX1(h2CorrelEnergyOCDBZmZm->GetXaxis()->GetXmin());
     ligne1H->SetX2(h2CorrelEnergyOCDBZmZm->GetXaxis()->GetXmax());
     h2CorrelEnergyOCDBZmZm->Draw("COLZ");
     ligne1H->Draw("SAME");
     fInv->Draw("SAME");
     h2CorrelEnergyOCDBZmZm->Write();
     c1->Update();
     }

 //Correlation plots, all SM :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(5);
 h2Amp[2*kNbSMtot+0]->Draw("COLZ");
 c1->cd(6);
 if (runMode == 0) h2E[2*kNbSMtot+0]->Draw("COLZ");
 c1->cd(3);
 h2Amp[2*kNbSMtot+1]->Draw("COLZ");
 c1->cd(4);
 if (runMode == 0) h2E[2*kNbSMtot+1]->Draw("COLZ");
 for (i=0;i<2;i++)
    {h2Amp[2*kNbSMtot+i]->Write();
     h2E[2*kNbSMtot+i]->Write();
     }
 c1->Update();
 
 for (iSM=0;iSM<kNbSMtot;iSM++)
    {if ((iSM%4) == 0)
        {c1->Update();
         ps->NewPage();
         c1->Clear();
         c1->Divide(2,4);
         }
     c1->cd(2*(iSM%4)+1);
     h2Amp[2*iSM+1]->Draw("COLZ");
     h2Amp[2*iSM+1]->Write();
     c1->cd(2*(iSM%4)+2);
     if (runMode == 0) h2E[2*iSM+1]->Draw("COLZ");
     h2E[2*iSM+1]->Write();
     }

 for (iSM=0;iSM<kNbSMtot;iSM++)
    {if ((iSM%4) == 0)
        {c1->Update();
         ps->NewPage();
         c1->Clear();
         c1->Divide(2,4);
         }
     c1->cd(2*(iSM%4)+1);
     h2Amp[2*iSM+0]->Draw("COLZ");
     h2Amp[2*iSM+0]->Write();
     c1->cd(2*(iSM%4)+2);
     if (runMode == 0) h2E[2*iSM+0]->Draw("COLZ");
     h2E[2*iSM+0]->Write();
     }
 
 c1->Update();

 //Correlation plots ratios vs intgs, all SM :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 for (j=0;j<6;j++)
    {c1->cd(j+1);
     h2RatioVsIntg[6*kNbSMtot+j]->Draw("COLZ");
     h2RatioVsIntg[6*kNbSMtot+j]->Write();
     if ((j%2) == 0) gPad->SetLogz();
     }
 c1->Update();
 
 //Correlation plots ratios vs amplitudes, all SM :
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 for (j=0;j<6;j++)
    {c1->cd(j+1);
     h2RatioVsAmp[6*kNbSMtot+j]->Draw("COLZ");
     h2RatioVsAmp[6*kNbSMtot+j]->Write();
     }
 c1->Update();
 
 
 printf("Drawing done.\n");

 
 ps->Close();
 printf("Closed ps file.\n");
 rootFileOut->Close();
 printf("Closed root file.\n");
 fclose(outFile);
 printf("Closed txt file, quit root.\n");

 return;
 }







