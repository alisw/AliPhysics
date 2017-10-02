

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


///
/// \file PrepareHV_SM1819.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief PrepareHV_SM1819
///
/// How to run :
///
/// root -b -q 'PrepareHV_SM1819.C(<smNb>)' >& output_prepareHV_SM1819_SMxxx.out &
///
/// Mapping convention for printf's : electronic.
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
void PrepareHV_SM1819(int nbSM)
{int iSM,iCol,iRow,sm,col,row,k,offset;
 double tabP2[4],tabCoeffs[kNbSMtot][kNbColMax][kNbRowMax],tabHV[kNbSMtot][kNbColMax][kNbRowMax];
 double hv,coeff;
 tabP2[0]=0.03;
 tabP2[1]=0.025;
 tabP2[2]=0.02;
 tabP2[3]=0.01;

 if ((nbSM != 18) && (nbSM != 19))
    {printf("Wrong SM number.\n");
     return;
     }
 
 FILE *fch = new FILE;
 fch = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/createOCDB_4_with2015data/multiplyPi0CalibrationFactors_TextToHisto_Final.txt","r");
 //Convention is offline mapping. Need to swap to electronic mapping.
 //Coeff convention is pi0 calib in the file read.
 for (iSM=0;iSM<kNbSMtot;iSM++)
    {for (iCol=0;iCol<kTabNbCol[SMdetType[iSM]];iCol++)
        {for (iRow=0;iRow<kTabNbRow[SMdetType[iSM]];iRow++)
            {fscanf(fch," %d %d %d %lf\n",&sm,&col,&row,&coeff);
             if ((sm != iSM) || (col != iCol) || (row != iRow)) printf("=========== Out of sync OCDB (%d %d %d) vs (%d %d %d) ============\n",sm,col,row,iSM,iCol,iRow);
	     if (nbSM == 18) tabCoeffs[sm][col][row]=1./coeff;
	     else tabCoeffs[sm][(kTabNbCol[SMdetType[iSM]]-1)-col][(kTabNbRow[SMdetType[iSM]]-1)-row]=1./coeff;
             }
         }
     }
 
 fclose(fch);

 fch = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/output_HVrecalculation_pass2_createSetBiasScripts/SMA12/NewBias.txt","r");
 //In this file, the voltages actually applied to SM18 during calib data-taking were rows 16-23. Those actually applied to SM19 were rows 0-7.
 //Convention is electronic mapping.
 offset=0;
 if (nbSM == 18) offset=16;
 for (iCol=0;iCol<kTabNbCol[SMdetType[nbSM]];iCol++)
    {for (iRow=0;iRow<kNbRowMax;iRow++)
        {fscanf(fch," %d %d %lf\n",&col,&row,&hv);
         if ((col != iCol) || (row != iRow)) printf("=========== Out of sync HV (%d %d %d) vs (%d %d %d) ============\n",nbSM,col,row,nbSM,iCol,iRow);
         tabHV[nbSM][col][row]=hv;
         }
     }

 fclose(fch);

 for (iCol=0;iCol<kTabNbCol[SMdetType[nbSM]];iCol++)
    {for (iRow=0;iRow<kNbRowDCALthird;iRow++)
        {coeff=tabCoeffs[nbSM][iCol][iRow];
         hv=tabHV[nbSM][iCol][iRow+offset];
         printf("%d %d %d %d coeff %f prevHV %7.3f",nbSM,iCol,iRow,iRow+16,coeff,hv);
         if (coeff == 1.) printf(" no HV change");
         else
            {printf(" HVs");
             for (k=0;k<4;k++)
                {printf(" %6.2f",hv+log(coeff)/tabP2[k]);
                 if (k == 1) printf("**");
                 }
             }
         printf("\n");
         }
     }
 
 
 return;
 }






