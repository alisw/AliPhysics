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


enum detType {kEMCAL,kEMCALthird,kDCAL,kDCALthird};
char detTypeString[][100]={"EMCAL","EMCALthird","DCAL","DCALthird"};

// see http://dsilverm.web.cern.ch/dsilverm/fee/addrP2.html
char SMP2Name[][100]={"SMA0","SMC0","SMA1","SMC1","SMA2","SMC2","SMA3","SMC3","SMA4","SMC4","SMA5","SMC5","SMA9","SMC9","SMA10","SMC10","SMA11","SMC11","SMA12","SMC12"};
int SMdetType[]={kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCALthird,kEMCALthird,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCALthird,kDCALthird};


///
/// \file ComparePi0CalibResults.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief Compare analysis results ...
///
/// How to compile :
///   .L ComparePi0CalibResults.C++
///
/// How to run (uncompiled) :
///
///   .x macros/ComparePi0CalibResults.C("toto")
///
/// CAUTION : some things at the beginning are "To be checked / changed by user" (search "customize")./
//
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
void ComparePi0CalibResults(TString ident="")
{int i,j,iRun,iSM,iCol,iRow,cmptPlots;
 char name1[150],name2[150];
 TFile *rootFileIn1,*rootFileIn2;
 
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
 
 TH1F *hInvMass1[kNbSMtot][kNbColMax][kNbRowMax];
 TH1F *hInvMass2[kNbSMtot][kNbColMax][kNbRowMax];
 
 FILE *outFile;
 char fchNameBase[200],fchNameRoot[200],fchNamePs[200],fchNameOut[200];
 sprintf(fchNameBase,"output_comparePi0CalibResults%s",ident.Data());
 sprintf(fchNameRoot,"%s.root",fchNameBase);
 sprintf(fchNamePs,"%s.ps",fchNameBase);
 sprintf(fchNameOut,"%s.out",fchNameBase);
 outFile=fopen(fchNameOut,"w");

 
 const int cWidth=500;
 const int cHeight=(int)(500*(30./21.));
 TCanvas *c1 = new TCanvas("c1","EMCal cosmics analysis",cWidth,cHeight);
 TPostScript *ps = new TPostScript(fchNamePs,111);

 sprintf(name1,"/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/TestLowMassPeakTowers/output_BadMapNoMask/output_calibPi0.root");
 sprintf(name2,"/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/TestLowMassPeakTowers/output_BadMapMask/output_calibPi0.root");
 rootFileIn1 = new TFile(name1,"READONLY");
 rootFileIn2 = new TFile(name2,"READONLY");
 
 for (iSM=0;iSM<kNbSMtot;iSM++)
 //for (iSM=0;iSM<2;iSM++)
    {if (iSM == 0) continue;
     if (iSM == 2) continue;
     if (iSM == 4) continue;
     if (iSM == 5) continue;
     if (iSM == 6) continue;
     if (iSM == 7) continue;
     if (iSM == 8) continue;
     if (iSM == 9) continue;
     if (iSM == 10) continue;
     if (iSM == 11) continue;
     if (iSM == 18) continue;
     if (iSM == 19) continue;
     for (iCol=0;iCol<kTabNbCol[SMdetType[iSM]];iCol++)
        {for (iRow=0;iRow<kTabNbRow[SMdetType[iSM]];iRow++)
	    {hInvMass1[iSM][iCol][iRow] = (TH1F*)rootFileIn1->Get(Form("%d_%d_%d",iSM,iCol,iRow));
	     }
	 }
     }

 TFile *rootFileOut = new TFile(fchNameRoot,"RECREATE");
 cmptPlots=0;
 ps->NewPage();
 c1->Clear();
 c1->Divide(3,4);
 for (iSM=0;iSM<kNbSMtot;iSM++)
 //for (iSM=0;iSM<2;iSM++)
    {printf("Beginning SM %d.\n",iSM);
     if (iSM == 0) continue;
     if (iSM == 2) continue;
     if (iSM == 4) continue;
     if (iSM == 5) continue;
     if (iSM == 6) continue;
     if (iSM == 7) continue;
     if (iSM == 8) continue;
     if (iSM == 9) continue;
     if (iSM == 10) continue;
     if (iSM == 11) continue;
     if (iSM == 18) continue;
     if (iSM == 19) continue;
     for (iCol=0;iCol<kTabNbCol[SMdetType[iSM]];iCol++)
        {for (iRow=0;iRow<kTabNbRow[SMdetType[iSM]];iRow++)
	    {hInvMass2[iSM][iCol][iRow] = (TH1F*)rootFileIn2->Get(Form("%d_%d_%d",iSM,iCol,iRow));
	     cmptPlots++;
	     if (!hInvMass2[iSM][iCol][iRow]) continue;
	     hInvMass2[iSM][iCol][iRow]->Divide(hInvMass1[iSM][iCol][iRow]);
	     c1->cd(cmptPlots);
             hInvMass2[iSM][iCol][iRow]->SetMaximum(1.1);
             hInvMass2[iSM][iCol][iRow]->SetMinimum(0.0);
             hInvMass2[iSM][iCol][iRow]->SetAxisRange(50.,250.,"X");
	     hInvMass2[iSM][iCol][iRow]->Draw();
	     hInvMass2[iSM][iCol][iRow]->Write();
	     if (cmptPlots == 3*4)
	        {c1->Update();
		 ps->NewPage();
		 c1->Clear();
                 c1->Divide(3,4);
		 cmptPlots=0;
		 }
	     }
	 }
     }
         


 

 TLine *lineSMborderVEMCAL = new TLine(kNbColEMCAL-0.5,-0.5,kNbColEMCAL-0.5,(int)((kNbSMEMCAL+1)/2)*kNbRowEMCAL+(int)((kNbSMEMCALthird+1)/2)*kNbRowEMCALthird-0.5);
 TLine *lineSMborderVDCALthird = new TLine(kNbColDCALthird-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5,kNbColDCALthird-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL+(int)((kNbSMDCALthird+1)/2)*kNbRowDCALthird-0.5);
 TLine *lineSMborderVDCAL1 = new TLine(kNbColDCAL-0.5,-0.5,kNbColDCAL-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5);
 TLine *lineSMborderVDCAL2 = new TLine(2*kNbColEMCAL-kNbColDCAL-0.5,-0.5,2*kNbColEMCAL-kNbColDCAL-0.5,(int)((kNbSMDCAL+1)/2)*kNbRowDCAL-0.5);
 TLine **lineSMborderHEMCAL,**lineSMborderHDCAL;
 lineSMborderHEMCAL = new TLine*[(int)((kNbSMEMCAL+1)/2)];
 lineSMborderHDCAL = new TLine*[(int)((kNbSMDCAL+1)/2)];
 for (i=0;i<(int)((kNbSMEMCAL+1)/2);i++) lineSMborderHEMCAL[i] = new TLine(-0.5,(i+1)*kNbRowEMCAL-0.5,2.*kNbColEMCAL-0.5,(i+1)*kNbRowEMCAL-0.5);
 for (i=0;i<(int)((kNbSMDCAL+1)/2);i++) lineSMborderHDCAL[i] = new TLine(-0.5,(i+1)*kNbRowDCAL-0.5,2.*kNbColDCALthird-0.5,(i+1)*kNbRowDCAL-0.5);

 
 printf("Drawing done.\n");

 
 ps->Close();
 printf("Closed ps file.\n");
 rootFileOut->Close();
 printf("Closed root file.\n");
 fclose(outFile);
 printf("Closed txt file, quit root.\n");

 return;
 }







