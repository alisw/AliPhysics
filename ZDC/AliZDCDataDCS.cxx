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
#include <TGraph.h>
#include <TH1.h>
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
   fGraphs("TGraph",kNGraphs),
   fIsProcessed(kFALSE)
{
  // Default constructor
  for(Int_t i=0; i<kNAliases; i++){   
     fAliasNames[i] = "";
     fCalibData[i] = 0.; 
     fTimeStamp[i] = 0.; 
  }
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
   fGraphs("TGraph",kNGraphs),
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
  fGraphs(data.fGraphs),
  fIsProcessed(data.fIsProcessed)
{

  // copy constructor

  for(int i=0;i<kNAliases;i++) {
    fAliasNames[i]=data.fAliasNames[i];
  }
  
    
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

  for(int i=0;i<kNAliases;i++) {
    fAliasNames[i]=data.GetAliasName(i);
    fCalibData[i]=data.GetCalibData(i);
    fTimeStamp[i]=data.GetTimeStamp(i);
  }

  return *this;
}

//---------------------------------------------------------------
AliZDCDataDCS::~AliZDCDataDCS() 
{
  // Destructor
  fGraphs.Clear("C");
}

//---------------------------------------------------------------
Bool_t AliZDCDataDCS::ProcessData(TMap& aliasMap)
{
   // Data processing

   if (fEndTime==fStartTime){
     AliError(Form(" Run with null time length: start time = %i = end time = %i",fStartTime,fEndTime));
     return kFALSE;
   }
   
   AliInfo(Form(" Start Time = %i",fStartTime));
   AliInfo(Form(" End Time = %i",fEndTime));
   AliInfo(Form(" Start Time DCS Query= %i",fStartTimeDCSQuery));
   AliInfo(Form(" End Time DCS Query= %i",fEndTimeDCSQuery));
   
   TObjArray *aliasArr;
   AliDCSValue* aValue;
  
   for(int j=0; j<kNAliases; j++){
      aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[j].Data());
      if(!aliasArr){
   	AliWarning(Form("Alias %s not found!", fAliasNames[j].Data()));
   	return kFALSE;
      }
      Introduce(j, aliasArr);

      if(aliasArr->GetEntries()<=2){
   	AliWarning(Form("Alias %s has just %d entries!",
   			fAliasNames[j].Data(),aliasArr->GetEntries()));
   	continue;
      }

      TIter iterarray(aliasArr);

      Double_t *time = new Double_t[aliasArr->GetEntries()];
      Double_t *val = new Double_t[aliasArr->GetEntries()];

      UInt_t ne=0;
      while((aValue = (AliDCSValue*) iterarray.Next())) {
   	val[ne] = aValue->GetFloat();
   	time[ne] = (Float_t) (aValue->GetTimeStamp());
   	fCalibData[ne] = val[ne];
   	fTimeStamp[ne] = time[ne];
   	ne++;
      }
      
//      if(j>=4) CreateGraph(j, aliasArr->GetEntries(), time, val); // fill graphs 
      //
      delete[] val;
      delete[] time;	      
   }
   //
   fIsProcessed=kTRUE;
   
   return kTRUE;
}

//---------------------------------------------------------------
void AliZDCDataDCS::Init()
{
   // Initialization
   
   TH1::AddDirectory(kFALSE);

   fGraphs.SetOwner(1);

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
   AliDebug(2,Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
   AliDebug(2,Form("	   %d DP values collected",entries));

}


//---------------------------------------------------------------
void AliZDCDataDCS::CreateGraph(int i, int dim, const Double_t *x, const Double_t *y)
{

   // Create graphics
   
   TGraph *gr = new(fGraphs[fGraphs.GetEntriesFast()]) TGraph(dim, x, y);

   gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetTitle(fAliasNames[i].Data());
//   AliDebug(2,Form("Array entries: %d",fGraphs.GetEntriesFast()));


}

//---------------------------------------------------------------
void AliZDCDataDCS::Draw(const Option_t* /*option*/)
{
  // Draw graphics

  fIsProcessed=1;
  if(!fIsProcessed) return;
  
  if(fGraphs.GetEntries()==0)  return;
  
  TCanvas *cg1;
  TString canvas1Name="ZN1_HVs";
  cg1=new TCanvas(canvas1Name,canvas1Name,40,40,600,600);
  cg1->Divide(2,2);
  cg1->cd(1);
  ((TGraph*) fGraphs.UncheckedAt(0))->SetMarkerStyle(20);
  ((TGraph*) fGraphs.UncheckedAt(0))->Draw("ALP");
  cg1->cd(2);
  ((TGraph*) fGraphs.UncheckedAt(1))->SetMarkerStyle(20);
  ((TGraph*) fGraphs.UncheckedAt(1))->Draw("ALP");
  cg1->cd(3);
  ((TGraph*) fGraphs.UncheckedAt(2))->SetMarkerStyle(20);
  ((TGraph*) fGraphs.UncheckedAt(2))->Draw("ALP");
  cg1->cd(4);
  ((TGraph*) fGraphs.UncheckedAt(3))->SetMarkerStyle(20);
  ((TGraph*) fGraphs.UncheckedAt(3))->Draw("ALP");
  
  TCanvas *cg2;
  TString canvas2Name="ZP1_HVs";
  cg2=new TCanvas(canvas2Name,canvas2Name,80,80,600,600);
  cg2->Divide(2,2);
  cg2->cd(1);
  ((TGraph*) fGraphs.UncheckedAt(5))->SetMarkerStyle(20);
  ((TGraph*) fGraphs.UncheckedAt(5))->Draw("ALP");
  cg2->cd(2);
  ((TGraph*) fGraphs.UncheckedAt(6))->SetMarkerStyle(20);
  ((TGraph*) fGraphs.UncheckedAt(6))->Draw("ALP");
  cg2->cd(3);
  ((TGraph*) fGraphs.UncheckedAt(7))->SetMarkerStyle(20);
  ((TGraph*) fGraphs.UncheckedAt(7))->Draw("ALP");
  cg2->cd(4);
  ((TGraph*) fGraphs.UncheckedAt(8))->SetMarkerStyle(20);
  ((TGraph*) fGraphs.UncheckedAt(8))->Draw("ALP");
 
}

