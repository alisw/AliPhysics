

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
/// \file ConvertOCDBtoText.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief Text to OCDB conversion
///
/// Parameters stored to a txt file, transfered to an OCDB file
///
/// How to run :
///
/// (CAUTION : customize first what has to be)
///
/// aliroot -b
/// .x ConvertOCDBtoText.C 
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
void ConvertOCDBtoText(void)
{  
  //CUSTOMIZE customize :
  //TFile * fIn = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/createOCDB_4_with2015data/Run0_999999999_v1_s0_OCDBcoeffsLHC15ij_EMCALnotValid.root","READ");
  TFile * fIn = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/createOCDB_4_with2015data/Run0_999999999_v0_s2_OCDBcoeffsLHC15ijAfterCalib_EMCALnotValid.root","READ");
  //TFile * fIn = new TFile("OCDBparamsCalib_allOneForTests.root","READ");

  FILE *txtFileOut = NULL;
  //CUSTOMIZE customize :
  txtFileOut = fopen("OCDBparamsAfterCalib2015.txt","w");
  //txtFileOut = fopen("OCDBparamsCalib_allOneForTests.txt","w");
  
  AliCDBEntry * cdb = (AliCDBEntry*) fIn->Get("AliCDBEntry");    
  AliEMCALCalibData * cparam =  cdb->GetObject();
  
  for(Int_t m = 0; m < kNbSMtot; m++){  
    for (Int_t c = 0; c < kTabNbCol[SMdetType[m]]; c++){
      for (Int_t r = 0; r < kTabNbRow[SMdetType[m]]; r++){
        
        Float_t adc=cparam->GetADCchannel(m,c,r);
        //Explanation of what the params below are : fprintf(txtFileOut,"SM %d, col %d, row %d, parameter  %1.4f\n",m, c, r,adc);
        //fprintf(txtFileOut,"%d %d %d %1.4f\n",m, c, r,adc);
        //Take %1.4f out as introduces an unnecessary rounding at the 1% level.
        fprintf(txtFileOut,"%d %d %d %f\n",m, c, r,adc);
        
      }
    }
  }
 
 printf("Done.\n");
 fclose(txtFileOut);
 
 return; 
 }




