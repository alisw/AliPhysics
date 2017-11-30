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


// voir http://dsilverm.web.cern.ch/dsilverm/fee/addrP2.html
char SMP2Name[][100]={"SMA0","SMC0","SMA1","SMC1","SMA2","SMC2","SMA3","SMC3","SMA4","SMC4"};
char SMcalibName[][100]={"US2","US1","EU2","EU1","US3","US5","US4","EU3","US7","US6"};
char SMnumber[][100]={"0","1","2","3","4","5","6","7","8","9"};

const int kNbCol=48;
const int kNbRow=24;
const int kNbSM=10;



///
/// \file MergeHVscanFilesVariousScans2.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief MergeHVscanFilesVariousScans2
///
/// How to run :
///
/// root -b -q 'macros/mergeHVscanFilesVariousScans2.C' > & output_mergeHVscanFilesVariousScans.out &
///
/// This code reads the 3 parameters p0,p1,p2 of the gain = f(HV) curve for various HV scans. 
/// When the parameters are found not to be valid for a tower, they are replaced by the parameters of 
/// the most recent previous HV scan found to have valid parameters.
///
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
void MergeHVscanFilesVariousScans2(void)
{int ism,icol,irow,iScan,cmpt,cmptEye;
 
 double p0,p1,p2,p1Max,p2Max;
 double p0tab[48][24]={{0.}};
 double p1tab[48][24]={{0.}};
 double p2tab[48][24]={{0.}};
 double yearTab[48][24]={{0.}};
 double monthTab[48][24]={{0.}};
 
  //defineMyPalette2011(30,5);
 
 //CUSTOMIZE customize :
 const int kNbHVscans=3;
 const int kNbMaxDiscardTowers=200;
   
 FILE *inFile,*outParamFile;
 char paramFileName[200],baseName[200];
 TString **HVscansPaths,**discardFiles;
 HVscansPaths = new TString*[kNbHVscans];
 discardFiles = new TString*[kNbHVscans];
 int HVscansYears[kNbHVscans];
 int HVscansMonths[kNbHVscans];
 double discardScans[kNbHVscans][10][48][24]={{{{0.}}}};
 
 //CUSTOMIZE customize :
 sprintf(baseName,"mergedHVscanParamsApr2013");
 //CUSTOMIZE customize :
 discardFiles[0] = new TString("000");
 HVscansPaths[0] = new TString("/cebaf/Web/EMCALpub/biasScanP2/feb2011/Hvscan");
 HVscansYears[0]=2011;
 HVscansMonths[0]=2;
 discardFiles[1] = new TString("000");
 HVscansPaths[1] = new TString("/cebaf/Web/EMCALpub/biasScanP2/Mars2012/Hvscan");
 HVscansYears[1]=2012;
 HVscansMonths[1]=3;
 discardFiles[2] = new TString("/cebaf/faivre/recherche/calibPi0/recalculateEMCAL_HV_Apr2013/HVscanMars2013/listBadScansByEye.txt");
 HVscansPaths[2] = new TString("/cebaf/Web/EMCALpub/biasScanP2/Mars2013/Hvscan"); 
 HVscansYears[2]=2013;
 HVscansMonths[2]=3;

 p1Max=50.0;
 p2Max=0.055;
 TH1F *h0 = new TH1F("h0","h0",200,-150.,150.);
 TH1F *h1 = new TH1F("h1","h1",200,0.,p1Max);
 TH1F *h2 = new TH1F("h2","h2",200,0.,p2Max);
 TH2F *h10 = new TH2F("h10","h10",100,-150.,150.,100,0.,p1Max);
 TH2F *h20 = new TH2F("h20","h20",100,-150.,150.,100,0.,p2Max);
 TH2F *h21 = new TH2F("h21","h21",100,0.,p1Max,100,0.,p2Max);
 h0->SetStats(0);
 h1->SetStats(0);
 h2->SetStats(0);
 h10->SetStats(0);
 h20->SetStats(0);
 h21->SetStats(0);
 h10->SetContour(30);
 h20->SetContour(30);
 h21->SetContour(30);
 
 printf("\nTower numbers given in Grenoble/electronic mapping convention.\n\n");
 
 for (iScan=0;iScan<kNbHVscans;iScan++)
    {inFile=fopen(discardFiles[iScan]->Data(),"r");
     if (inFile)
        {printf("Read discard file scan %d : %s\nFound towers : ",iScan,discardFiles[iScan]->Data());
         while (fscanf(inFile," %d %d %d\n",&ism,&icol,&irow)>0)
            {discardScans[iScan][ism][icol][irow]=1;
             printf("(%d,%d,%d) ",ism,icol,irow);
             }
         fclose(inFile);
         }
         else printf("By-eye discarded LED HV-scan file %s not found (HV scan %d)\n",discardFiles[iScan]->Data(),iScan);
     printf("\n");
     }
 
 gSystem->Exec(Form("mkdir %s",baseName));
 
 for (ism=0;ism<kNbSM;ism++)
    {printf("\n____________________________________\nDoing SM %d...\n",ism);
     gSystem->Exec(Form("mkdir %s/HvscanSM%s_%s_%s",baseName,SMnumber[ism],SMcalibName[ism],SMP2Name[ism]));
     for (icol=0;icol<kNbCol;icol++)
        {for (irow=0;irow<kNbRow;irow++)
            {p0tab[icol][irow]=0.;
             p1tab[icol][irow]=0.;
             p2tab[icol][irow]=0.;
             }
         }
     for (iScan=kNbHVscans-1;iScan>=0;iScan--)
        {if (iScan==kNbHVscans-1) printf("*** LED scan %d (%02d/%d)",iScan+1,HVscansMonths[iScan],HVscansYears[iScan]);
            else printf("*** LED scan %d (%02d/%d)\n    Check scan %d : ",iScan+1,HVscansMonths[iScan],HVscansYears[iScan],iScan+2);
         sprintf(paramFileName,"%sSM%s_%s_%s/parameters.txt",HVscansPaths[iScan]->Data(),SMnumber[ism],SMcalibName[ism],SMP2Name[ism]);
         inFile=fopen(paramFileName,"r");
         cmpt=0;
         cmptEye=0;
         if (inFile)
            {while (fscanf(inFile," %d %d %lf %lf %lf\n",&icol,&irow,&p0,&p1,&p2)>0)
                {if (((p1tab[icol][irow]==0.) && (p2tab[icol][irow]==0.)) || ((p1tab[icol][irow]<=0.) || (p2tab[icol][irow]<=0.))) // Previous LED HV-scan flagged not good, or has bad parameters.
                    {if (iScan != (kNbHVscans-1)) printf("(%d,%d) ",icol,irow);
                     p0tab[icol][irow]=p0;
                     p1tab[icol][irow]=p1;
                     p2tab[icol][irow]=p2;
                     yearTab[icol][irow]=HVscansYears[iScan];
                     monthTab[icol][irow]=HVscansMonths[iScan];
                     cmpt++;
                     }
                    else
                    {if ((iScan < kNbHVscans-1) && (discardScans[iScan+1][ism][icol][irow] == 1)) // Previous LED HV-scan tagged as bad by eye-inspection.
                        {printf("eye(%d,%d) ",icol,irow);
                         p0tab[icol][irow]=p0;
                         p1tab[icol][irow]=p1;
                         p2tab[icol][irow]=p2;
                         yearTab[icol][irow]=HVscansYears[iScan];
                         monthTab[icol][irow]=HVscansMonths[iScan];
                         cmptEye++;
                         }
                     }
                 }
             fclose(inFile);
             }
             else printf("File %s not found\n",paramFileName);
         printf("\n");
         if ((iScan != (kNbHVscans-1)) && (cmpt > 0)) printf("    %d towers found with missing info.\n",cmpt);
         if ((iScan != (kNbHVscans-1)) && (cmptEye > 0)) printf("    %d towers tagged as bad by eye-inspection.\n",cmptEye);
         }

     printf("*** Writing to file : ");
     cmpt=0;
     outParamFile=fopen(Form("%s/HvscanSM%s_%s_%s/parameters.txt",baseName,SMnumber[ism],SMcalibName[ism],SMP2Name[ism]),"w");
     for(icol=0;icol<kNbCol;icol++)
        {for(irow=0;irow<kNbRow;irow++)
            {//fprintf(outParamFile,"%2d %2d %lf %lf %lf %d %d\n",icol,irow,p0tab[icol][irow],p1tab[icol][irow],p2tab[icol][irow],yearTab[icol][irow],monthTab[icol][irow]);
             fprintf(outParamFile,"%2d %2d %lf %lf %lf\n",icol,irow,p0tab[icol][irow],p1tab[icol][irow],p2tab[icol][irow]);
             if (((p1tab[icol][irow]==0.) && (p2tab[icol][irow]==0.)) || ((p1tab[icol][irow]<=0.) || (p2tab[icol][irow]<=0.))) // All LED HV-scans flagged not good, or had bad parameters.
                {printf("(%d,%d) ",icol,irow);
                 cmpt++;
                 }
                else
                {h0->Fill(p0tab[icol][irow]);
                 h1->Fill(p1tab[icol][irow]);
                 h2->Fill(p2tab[icol][irow]);
                 h10->Fill(p0tab[icol][irow],p1tab[icol][irow]);
                 h20->Fill(p0tab[icol][irow],p2tab[icol][irow]);
                 h21->Fill(p1tab[icol][irow],p2tab[icol][irow]);
                 //if (p2tab[icol][irow] < 0.006) printf(">>>>>Check tower (%d %d %d), params (%f %f %f)\n",ism,icol,irow,p0tab[icol][irow],p1tab[icol][irow],p2tab[icol][irow]);
                 //if (p2tab[icol][irow] < 0.01*TMath::Exp(-0.06*p1tab[icol][irow])) printf(">>>>>Check tower (%d %d %d), params (%f %f %f)\n",ism,icol,irow,p0tab[icol][irow],p1tab[icol][irow],p2tab[icol][irow]);
                 }
             }
         }
     fclose(outParamFile);
     printf("\n    %d towers remain with missing info.\n",cmpt);
     }
 
 const int cWidth=500;
 const int cHeight=(int)(500*(29./21.));
 TCanvas *c1 = new TCanvas("c1","EMCal cosmics analysis",cWidth,cHeight);
 
 c1->Divide(2,3);
 c1->cd(1);
 gPad->SetLogy();
 h0->Draw();
 c1->cd(3);
 gPad->SetLogy();
 h1->Draw();
 c1->cd(5);
 gPad->SetLogy();
 h2->Draw();
 c1->cd(2);
 gPad->SetLogz();
 h10->Draw("COLZ");
 c1->cd(4);
 gPad->SetLogz();
 h20->Draw("COLZ");
 c1->cd(6);
 gPad->SetLogz();
 h21->Draw("COLZ");
 c1->Update();

 return;
 }










