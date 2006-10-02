/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


//
// This is an example of a container class
// of the data retreieved from the DCS archive DB.
// It is called by the detector's Preprocessor and
// it is stored in the CDB.
//

#include "AliTPCDataDCS.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TDatime.h>
#include <TStyle.h>
#include <TCanvas.h>

ClassImp(AliTPCDataDCS)

//---------------------------------------------------------------
AliTPCDataDCS::AliTPCDataDCS():
	TObject(),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
	fGraphs("TGraph",kNGraphs),
	fIsProcessed(kFALSE)
{
// default constructor

	for(int i=0;i<kNHistos;i++) fHv[i]=0;
	fFunc=0;
}

//---------------------------------------------------------------
AliTPCDataDCS::AliTPCDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fGraphs("TGraph",kNGraphs),
	fIsProcessed(kFALSE)
{
// constructor

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", nRun,
	TTimeStamp(startTime).AsString(),
	TTimeStamp(endTime).AsString()));

	for(int i=0;i<kNHistos;i++) fHv[i]=0;
	fFunc=0;

}

//---------------------------------------------------------------
AliTPCDataDCS::~AliTPCDataDCS()
{
// destructor

	Clean();
}

//---------------------------------------------------------------
void AliTPCDataDCS::Clean()
{
// destructor

	for(int i=0;i<kNHistos;i++) {
		if(fHv[i]) delete fHv[i];
		fHv[i]=0;
	}
	fGraphs.Clear("C");
	if(fFunc) delete fFunc;
	fFunc=0;
}

//---------------------------------------------------------------
void AliTPCDataDCS::Init()
{
// Init alias names and histos

	TH1::AddDirectory(kFALSE);

	fGraphs.SetOwner(1);

	for(int i=0;i<kNAliases;i++){
		fAliasNames[i] = "TpcHvSect0";
		fAliasNames[i] += i;
		fAliasNames[i] += ".FloatValue";
	}

	for(int i=0;i<kNHistos;i++){
		fHv[i] = new TH1F(fAliasNames[i].Data(),fAliasNames[i].Data(), 20, kHvMin, kHvMax);
		fHv[i]->GetXaxis()->SetTitle("Hv");
	}
}

//---------------------------------------------------------------
void AliTPCDataDCS::ProcessData(TMap& aliasMap)
{
// process TMap of alias values retrieved from the DCS archive DB

	Clean();
	Init();

	TObjArray *aliasArr;
	AliDCSValue* aValue;
	for(int j=0; j<kNAliases; j++){
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
		while ((aValue = (AliDCSValue*) iterarray.Next())) {
		val[ne] = aValue->GetFloat();
		time[ne] = (Double_t) (aValue->GetTimeStamp());
		// fill histos (alias 0-2)
		if(j < 3) fHv[j]->Fill(val[ne]);
		ne++;
		}
		// fill graphs (alias 3-5)
		if(j >= 3) CreateGraph(j, aliasArr->GetEntries(), time, val);
		delete[] val;
		delete[] time;
	}

	// calculate mean and rms of the first two histos
	for(int i=0;i<kNHistos;i++){
		fMean[i] = fHv[i]->GetMean();
		fWidth[i] = fHv[i]->GetRMS();
	}

	// pol1 fit of the first graph
//	if(fGraphs.GetEntries() > 0){
//		((TGraph*) fGraphs.UncheckedAt(0))->Fit("pol1");
//		fFunc = ((TGraph*) fGraphs.UncheckedAt(0))->GetFunction("pol1");
//	}

	fIsProcessed=kTRUE;


}

//---------------------------------------------------------------
void AliTPCDataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr)
{
// Say something about alias value array being processed

	int entries=aliasArr->GetEntries();
	AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
	AliInfo(Form("    	%d DP values collected",entries));

}

//---------------------------------------------------------------
void AliTPCDataDCS::CreateGraph(int i, int dim, const Double_t *x, const Double_t *y)
{
// create TGraph to be filled with DP values

	TGraph *gr = new(fGraphs[fGraphs.GetEntriesFast()]) TGraph(dim, x, y);

	gr->GetXaxis()->SetTimeDisplay(1);
	gr->SetTitle(fAliasNames[i].Data());

	AliInfo(Form("Array entries: %d",fGraphs.GetEntriesFast()));


}

//---------------------------------------------------------------
void AliTPCDataDCS::Draw(const Option_t* /*option*/)
{
// Draw all histos and graphs

  if(!fIsProcessed) return;

  TCanvas *ch;
  TString canvasHistoName="Histos";
  ch=new TCanvas(canvasHistoName,canvasHistoName,20,20,600,600);
  ch->Divide(2,2);
  ch->cd(1);
  fHv[0]->Draw();
  ch->cd(2);
  fHv[1]->Draw();
  ch->cd(3);
  fHv[2]->Draw();


  if(fGraphs.GetEntries() == 0) return;

  TCanvas *cg;
  TString canvasGraphName="Graphs";
  cg=new TCanvas(canvasGraphName,canvasGraphName,40,40,600,600);
  cg->Divide(2,2);
  cg->cd(1);
  ((TGraph*) fGraphs.UncheckedAt(0))->Draw("alp");

  cg->cd(2);
  ((TGraph*) fGraphs.UncheckedAt(1))->Draw("alp");
  cg->cd(3);
  ((TGraph*) fGraphs.UncheckedAt(2))->Draw("alp");

  if(fFunc){
  	cg->cd(4);
  	fFunc->Draw("l");
  }

}

