#include "AliZDCDataDCS.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TGraph.h>
#include <TDatime.h>
#include <TStyle.h>
#include <TCanvas.h>

ClassImp(AliZDCDataDCS)

//---------------------------------------------------------------
AliZDCDataDCS::AliZDCDataDCS():
	TObject(),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
	fGraphs("TGraph",kNGraphs),
	fIsProcessed(kFALSE)
{}

//---------------------------------------------------------------
AliZDCDataDCS::AliZDCDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fGraphs("TGraph",kNGraphs),
	fIsProcessed(kFALSE)
{
	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", nRun,
	TTimeStamp(startTime).AsString(),
	TTimeStamp(endTime).AsString()));

	Init();

}

//---------------------------------------------------------------
AliZDCDataDCS::~AliZDCDataDCS() {

	fGraphs.Clear("C");
}

//---------------------------------------------------------------
void AliZDCDataDCS::ProcessData(TMap& aliasMap){

	TObjArray *aliasArr;
	AliDCSValue* aValue;
	for(int j=4; j<kNAliases; j++){
		aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[j].Data());
		if(!aliasArr){
			AliError(Form("Alias %s not found!", fAliasNames[j].Data()));
			continue;
		}
		Introduce(j, aliasArr);

		if(aliasArr->GetEntries()<2){
			AliError(Form("Alias %s has just %d entries!",
					fAliasNames[j].Data(),aliasArr->GetEntries()));
			continue;
		}

		TIter iterarray(aliasArr);

		Double_t *time = new Double_t[aliasArr->GetEntries()];
		Double_t *val = new Double_t[aliasArr->GetEntries()];

		UInt_t ne=0;
		while((aValue = (AliDCSValue*) iterarray.Next())) {
		  val[ne] = aValue->GetSimpleValue().GetFloat();
		  time[ne] = (Double_t) (aValue->GetTimeStamp());
		  ne++;
		}

		// fill graphs 
		CreateGraph(j, aliasArr->GetEntries(), time, val);
		delete[] val;
		delete[] time;
	}


	fIsProcessed=kTRUE;


}

//---------------------------------------------------------------
void AliZDCDataDCS::Init(){

	TH1::AddDirectory(kFALSE);

	fGraphs.SetOwner(1);

	for(int i=0;i<kNAliases;i++){
		/*if(i<4){
		  fAliasNames[i] = "ZDC.Position";
		}
		else fAliasNames[i] = "ZDC.HVValue";*/
		fAliasNames[i] = "DCSAlias";
		fAliasNames[i] += i;
	}

}

//---------------------------------------------------------------
void AliZDCDataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr){

	int entries=aliasArr->GetEntries();
	AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
	AliInfo(Form("    	%d DP values collected",entries));

}


//---------------------------------------------------------------
void AliZDCDataDCS::CreateGraph(int i, int dim, const Double_t *x, const Double_t *y)
{

	TGraph *gr = new(fGraphs[fGraphs.GetEntriesFast()]) TGraph(dim, x, y);

	gr->GetXaxis()->SetTimeDisplay(1);
	gr->SetTitle(fAliasNames[i].Data());

	AliInfo(Form("Array entries: %d",fGraphs.GetEntriesFast()));


}

//---------------------------------------------------------------
void AliZDCDataDCS::Draw(const Option_t* /*option*/)
{
// Draw graphs

  fIsProcessed=1;
  if(!fIsProcessed) return;
  
  if(fGraphs.GetEntries()==0)  return;
  
  TCanvas *cg1;
  TString canvas1Name="Graphs1";
  cg1=new TCanvas(canvas1Name,canvas1Name,40,40,600,600);
  cg1->Divide(2,2);
  cg1->cd(1);
  ((TGraph*) fGraphs.UncheckedAt(0))->Draw("AP");
  cg1->cd(2);
  ((TGraph*) fGraphs.UncheckedAt(1))->Draw("AP");
  cg1->cd(3);
  ((TGraph*) fGraphs.UncheckedAt(2))->Draw("AP");
  cg1->cd(4);
  ((TGraph*) fGraphs.UncheckedAt(3))->Draw("AP");
 
}

