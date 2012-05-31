#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TSystem.h>
#include <TMap.h>
#include <TString.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TGrid.h>
#include "STEER/AliDCSValue.h"
#include "STEER/AliCDBEntry.h"
#include "ZDC/AliZDCDataDCS.h"

#endif

void CheckZDCDCSMap(Int_t nRun=0)
{
   if(nRun==0){
     printf("\n\n YOU MUST PROVIDE A RUN NUMBER != 0!!! \n\n");
     return;
   }
  
   TGrid::Connect("alien:",0,0,"t");
  
   char fName[150];
   sprintf(fName,"alien:///alice/data/2010/Reference/ZDC/DCS/Data/Run%d_%d_v1_s0.root",nRun,nRun);
  
   TString  aliasNames[28];
   aliasNames[0]  = "ZDC_ZNA_POS.actual.position";
   aliasNames[1]  = "ZDC_ZPA_POS.actual.position";
   aliasNames[2]  = "ZDC_ZNC_POS.actual.position";
   aliasNames[3]  = "ZDC_ZPC_POS.actual.position";
   aliasNames[4]  = "ZDC_ZNA_HV0.actual.vMon";
   aliasNames[5]  = "ZDC_ZNA_HV1.actual.vMon";
   aliasNames[6]  = "ZDC_ZNA_HV2.actual.vMon";
   aliasNames[7]  = "ZDC_ZNA_HV3.actual.vMon";
   aliasNames[8]  = "ZDC_ZNA_HV4.actual.vMon";
   aliasNames[9]  = "ZDC_ZPA_HV0.actual.vMon";
   aliasNames[10] = "ZDC_ZPA_HV1.actual.vMon";
   aliasNames[11] = "ZDC_ZPA_HV2.actual.vMon";
   aliasNames[12] = "ZDC_ZPA_HV3.actual.vMon";
   aliasNames[13] = "ZDC_ZPA_HV4.actual.vMon";
   aliasNames[14] = "ZDC_ZNC_HV0.actual.vMon";
   aliasNames[15] = "ZDC_ZNC_HV1.actual.vMon";
   aliasNames[16] = "ZDC_ZNC_HV2.actual.vMon";
   aliasNames[17] = "ZDC_ZNC_HV3.actual.vMon";
   aliasNames[18] = "ZDC_ZNC_HV4.actual.vMon";
   aliasNames[19] = "ZDC_ZPC_HV0.actual.vMon";
   aliasNames[20] = "ZDC_ZPC_HV1.actual.vMon";
   aliasNames[21] = "ZDC_ZPC_HV2.actual.vMon";
   aliasNames[22] = "ZDC_ZPC_HV3.actual.vMon";
   aliasNames[23] = "ZDC_ZPC_HV4.actual.vMon";
   aliasNames[24] = "ZDC_ZEM_HV0.actual.vMon";
   aliasNames[25] = "ZDC_ZEM_HV1.actual.vMon";   
   aliasNames[26] = "ZDC_REFA_HV.actual.vMon";
   aliasNames[27] = "ZDC_REFC_HV.actual.vMon";

  TFile *file = TFile::Open(fName);
  //file->ls();
  
  AliCDBEntry *entry = (AliCDBEntry*)file->Get("AliCDBEntry");
  TMap *dcsAliasMap = dynamic_cast<TMap*> (entry->GetObject());
  dcsAliasMap->SetOwner(1);          
  //dcsAliasMap->Print("");
  TObjArray     *aliasArr;
  AliDCSValue   *value;  
   
  for(int j=0; j<28; j++){
  
     aliasArr = (TObjArray*)  (dcsAliasMap->GetValue(aliasNames[j].Data()));
     if(!aliasArr){
       printf("Alias %s has no DP value stored!\n", aliasNames[j].Data());
       continue;
     }
     
     Int_t nentries = aliasArr->GetEntries();
     printf("************ Alias: %s has  %d DP values collected\n", 
     	aliasNames[j].Data(), nentries);

     TIter iterarray(aliasArr);

     UInt_t  *time = new UInt_t[nentries];
     Float_t *val  = new Float_t[nentries];

     UInt_t ne=0;
     while((value = (AliDCSValue*) iterarray.Next())) {
       time[ne] = value->GetTimeStamp();
       val[ne] = value->GetFloat();
       //printf(" %d  - time %d  value %1.4f\n",ne, 
       //        value->GetTimeStamp(),value->GetFloat());
       printf(" %s\n", value->ToString());
       ne++;
     }
  }
  
                    
}
