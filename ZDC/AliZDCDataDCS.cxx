///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class for ZDC DCS data                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCDataDCS.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TCanvas.h>
#include <TDatime.h>
#include <TString.h>
#include <TObjString.h>
#include <TStyle.h>
#include <TTimeStamp.h>

ClassImp(AliZDCDataDCS)

//---------------------------------------------------------------
AliZDCDataDCS::AliZDCDataDCS():
   TObject(),
   fRun(0),
   fStartTime(0),
   fEndTime(0),
   fStartTimeDCSQuery(0),
   fEndTimeDCSQuery(0),
//   fTimeStamp(0x0), 
//   fHVData(0x0), 
   fIsProcessed(kFALSE)
{
  // Default constructor
  for(Int_t i=0; i<kNAliases; i++) fAliasNames[i] = "";
  for(Int_t i=0; i<kNAlignDet; i++) fAlignData[i] = 0.;
}

//---------------------------------------------------------------
AliZDCDataDCS::AliZDCDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime,
			 UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery):
   TObject(),
   fRun(nRun),
   fStartTime(startTime),
   fEndTime(endTime),
   fStartTimeDCSQuery(startTimeDCSQuery),
   fEndTimeDCSQuery(endTimeDCSQuery),
//   fTimeStamp(0x0), 
//   fHVData(0x0), 
   fIsProcessed(kFALSE)
{
   // Standard constructor
   
   AliDebug(2,Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s \n\tStartTime DCS Query %s \n\tEndTime DCS Query %s", nRun,
   TTimeStamp(startTime).AsString(),
   TTimeStamp(endTime).AsString(), 
   TTimeStamp(startTimeDCSQuery).AsString(), 
   TTimeStamp(endTimeDCSQuery).AsString()));

   Init();

}

//---------------------------------------------------------------
AliZDCDataDCS::AliZDCDataDCS(const AliZDCDataDCS & data):
  TObject(data), 
  fRun(data.fRun),
  fStartTime(data.fStartTime),
  fEndTime(data.fEndTime),
  fStartTimeDCSQuery(data.fStartTimeDCSQuery),
  fEndTimeDCSQuery(data.fEndTimeDCSQuery),
//  fTimeStamp(data.fTimeStamp), 
//  fHVData(data.fHVData), 
  fIsProcessed(data.fIsProcessed)
{

  // copy constructor

  for(int i=0;i<kNAliases;i++) fAliasNames[i] = data.fAliasNames[i];
  for(int i=0;i<kNAlignDet;i++) fAlignData[i] = data.fAlignData[i];
    
}

//---------------------------------------------------------------
AliZDCDataDCS& AliZDCDataDCS:: operator=(const AliZDCDataDCS & data) { 

  // assignment operator

  if (this == &data)
    return *this;

  TObject::operator=(data);
  fRun = data.GetRun();
  fStartTime = data.GetStartTime();
  fEndTime = data.GetEndTime();
  fStartTimeDCSQuery = data.GetStartTimeDCSQuery();
  fEndTimeDCSQuery = data.GetEndTimeDCSQuery();
//  fTimeStamp  = data.GetTimeStamp();
//  fHVData  = data.GetHVData();
  fIsProcessed = data.fIsProcessed; 

  for(int i=0;i<kNAliases;i++) fAliasNames[i] = data.GetAliasName(i);
  for(int i=0;i<kNAlignDet;i++) fAlignData[i] = data.GetAlignData(i);

  return *this;
}

//---------------------------------------------------------------
AliZDCDataDCS::~AliZDCDataDCS() 
{
  // Destructor
}

//---------------------------------------------------------------
Bool_t AliZDCDataDCS::ProcessData(TMap& aliasMap)
{
   // Data processing

   if(!(fAliasNames[0])) Init();
   
   AliInfo(Form(" Start Time = %i",fStartTime));
   AliInfo(Form(" End Time = %i",fEndTime));
   AliInfo(Form(" Start Time DCS Query= %i",fStartTimeDCSQuery));
   AliInfo(Form(" End Time DCS Query= %i",fEndTimeDCSQuery));

   if (fEndTime==fStartTime){
     AliError(Form(" Run with null time length: start time = %i = end time = %i",fStartTime,fEndTime));
     return kFALSE;
   }
   
   TObjArray   *aliasArr;
   AliDCSValue *aValue;
  
   for(int j=0; j<kNAliases; j++){
      //printf(" Processing alias %d  aliasName %s \n", j, fAliasNames[j].Data());
      
      aliasArr = (TObjArray*) (aliasMap.GetValue(fAliasNames[j].Data()));
      if(!aliasArr){
   	AliWarning(Form("Alias %s not found!", fAliasNames[j].Data()));
   	//printf(" AliZDCDataDCS: Alias %s not found!\n", fAliasNames[j].Data());
	continue;
      }

      Introduce(j, aliasArr);

      Int_t nentries = aliasArr->GetEntries();
      if(nentries<=2){
        AliWarning(Form("Alias %s has just %d entries!", fAliasNames[j].Data(), nentries));
//        continue;
      }

      Float_t *time = new Float_t[nentries];
      Float_t *val  = new Float_t[nentries];

      TIter iterarray(aliasArr);
      
      UInt_t ne=0;
      Float_t sum=0.;
      Int_t nMeasures=0;
      while((aValue = (AliDCSValue*) iterarray.Next())){
        val[ne] = aValue->GetFloat();
        time[ne] = (Float_t) (aValue->GetTimeStamp());
        if(j<4){
	  sum += val[ne];
	  nMeasures++;
	}
	else{
	  //fHVData[ne] = val[ne];
          //fTimeStamp[ne] = time[ne];
	}
        ne++;
      }
      //
      if(j<4 && nMeasures!=0) fAlignData[j] = sum/nMeasures;
      
      delete[] val;
      delete[] time;   
   }
  
   fIsProcessed=kTRUE;
   return kTRUE;
   
}

//---------------------------------------------------------------
void AliZDCDataDCS::Init()
{
   // Initialization

   fAliasNames[0] = "ZDC_ZNA_POS.actual.position";
   fAliasNames[1] = "ZDC_ZPA_POS.actual.position";
   fAliasNames[2] = "ZDC_ZNC_POS.actual.position";
   fAliasNames[3] = "ZDC_ZPC_POS.actual.position";
   //
   fAliasNames[4]  = "ZDC_ZNA_HV0.actual.vMon";
   fAliasNames[5]  = "ZDC_ZNA_HV1.actual.vMon";
   fAliasNames[6]  = "ZDC_ZNA_HV2.actual.vMon";
   fAliasNames[7]  = "ZDC_ZNA_HV3.actual.vMon";
   fAliasNames[8]  = "ZDC_ZNA_HV4.actual.vMon";
   //
   fAliasNames[9]   = "ZDC_ZPA_HV0.actual.vMon";
   fAliasNames[10]  = "ZDC_ZPA_HV1.actual.vMon";
   fAliasNames[11]  = "ZDC_ZPA_HV2.actual.vMon";
   fAliasNames[12]  = "ZDC_ZPA_HV3.actual.vMon";
   fAliasNames[13]  = "ZDC_ZPA_HV4.actual.vMon";
   //
   fAliasNames[14]  = "ZDC_ZNC_HV0.actual.vMon";
   fAliasNames[15]  = "ZDC_ZNC_HV1.actual.vMon";
   fAliasNames[16]  = "ZDC_ZNC_HV2.actual.vMon";
   fAliasNames[17]  = "ZDC_ZNC_HV3.actual.vMon";
   fAliasNames[18]  = "ZDC_ZNC_HV4.actual.vMon";
   //
   fAliasNames[19]  = "ZDC_ZPC_HV0.actual.vMon";
   fAliasNames[20]  = "ZDC_ZPC_HV1.actual.vMon";
   fAliasNames[21]  = "ZDC_ZPC_HV2.actual.vMon";
   fAliasNames[22]  = "ZDC_ZPC_HV3.actual.vMon";
   fAliasNames[23]  = "ZDC_ZPC_HV4.actual.vMon";
   //
   fAliasNames[24]  = "ZDC_ZEM_HV0.actual.vMon";
   fAliasNames[25]  = "ZDC_ZEM_HV1.actual.vMon";
   //
   fAliasNames[26]  = "ZDC_REFA_HV.actual.vMon";
   fAliasNames[27]  = "ZDC_REFC_HV.actual.vMon";
 
}

//---------------------------------------------------------------
void AliZDCDataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr)
{
   // Getting array of DCS aliases
   
   int entries = aliasArr->GetEntries();
   printf("************ Alias: %s has  %d DP values collected\n",
   	fAliasNames[numAlias].Data(),entries);

}
