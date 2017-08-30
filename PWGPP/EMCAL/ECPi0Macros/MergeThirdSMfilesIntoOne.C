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


#include <TSystem.h> //To have gSystem known at compilation time
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




namespace std {} using namespace std;


const double kMaxHV = 395.; // default max voltage limit; could possibly be relaxed in future
const double kMinHV = 210.; // Hexa encoding limit
const double coefFactorWanted=0.0162;


// voir http://dsilverm.web.cern.ch/dsilverm/fee/addrP2.html
 char SMP2Name[][100]={"SMA0","SMC0","SMA1","SMC1","SMA2","SMC2","SMA3","SMC3","SMA4","SMC4","SMA5","SMC5","SMA9","SMC9","SMA10","SMC10","SMA11","SMC11","SMA12","SMC12"};
 char SMcalibName[][100]={"US2","US1","EU2","EU1","US3","US5","US4","EU3","US7","US6","US8C0","US8C2","DCN1","DJP1","DUS2","DUS1","DJP2","DUS3","CN1A","CN1C"};
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

 int lastSM;


///
/// \file MergeThirdSMfilesIntoOne.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief MergeThirdSMfilesIntoOne
///
/// How to run :
///
/// root -b -q 'MergeThirdSMfilesIntoOne.C++()'
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
void MergeThirdSMfilesIntoOne()
{int iCol,iRow,col,row;
 double value;
 double tabValueSM_0to7[kNbColEMCALthird][kNbRowEMCALthird];
 double tabValueSM_16to23[kNbColEMCALthird][kNbRowEMCALthird];
 
 FILE *fileIn_0to7,*fileIn_16to23,*fileOut;

 const char fchNameIn_0to7[] = "/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/output_HVrecalculation_EMCALthirds_pass2/SMA5/NewBias.txt";
 const char fchNameIn_16to23[] = "/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/output_HVrecalculation_EMCALthirds_pass2/SMC5/NewBias.txt";
 const char fchNameOut[] = "NewFile_merged.txt";
 
 printf("\nRunning with : \n");
 printf("   - file rows  0-7  : %s\n",fchNameIn_0to7);
 printf("   - file rows 16-23 : %s\n",fchNameIn_16to23);
 
 fileIn_0to7 = fopen(fchNameIn_0to7,"r");
 if (!fileIn_0to7) {printf("File %s can not be found\n",fchNameIn_0to7); exit(-1);}
 while (fscanf(fileIn_0to7," %d %d %lf\n",&col,&row,&value) > 0)
    {tabValueSM_0to7[col][row]=value;
     }
 fclose(fileIn_0to7);

 fileIn_16to23 = fopen(fchNameIn_16to23,"r");
 if (!fileIn_16to23) {printf("File %s can not be found\n",fchNameIn_16to23); exit(-1);}
 while (fscanf(fileIn_16to23," %d %d %lf\n",&col,&row,&value) > 0)
    {tabValueSM_16to23[col][row]=value;
     }
 fclose(fileIn_16to23);
 
 fileOut = fopen(fchNameOut,"w");
 for (iCol=0;iCol<kNbColMax;iCol++)
    {for (iRow=0;iRow<kNbRowMax;iRow++)
        {if (iRow < kNbRowEMCALthird) fprintf(fileOut,"%d %d %7.3f\n",iCol,iRow,tabValueSM_0to7[iCol][iRow]);
         if ((iRow >= kNbRowEMCALthird) && (iRow < 2*kNbRowEMCALthird)) fprintf(fileOut,"%d %d 0\n",iCol,iRow);
         if (iRow >= 2*kNbRowEMCALthird) fprintf(fileOut,"%d %d %7.3f\n",iCol,iRow,tabValueSM_16to23[iCol][iRow-2*kNbRowEMCALthird]);
         }
     }
 fclose (fileOut);
 


 return;
 }








