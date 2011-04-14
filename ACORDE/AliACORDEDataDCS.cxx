/****************************************************
	
  AliACORDEDataDCS class
  create to make a pointer to the
  ACORDE data DCS points

  Author: Pedro Gonzalez (CIEMAT, Madrid)
  ACORDE-DCS creator: Mario Ivan Martinez Hdez
			<mim@fcfm.buap.mx>

  Last update: Fix of coding violations
  Mario Rodriguez C. (FCFM-BUAP)
  <mrodrigu@mail.cern.ch>

*****************************************************/
#include "AliACORDEDataDCS.h"

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

ClassImp(AliACORDEDataDCS)

//---------------------------------------------------------------
AliACORDEDataDCS::AliACORDEDataDCS():
	TObject(),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
	fGraphs("TGraph",kNGraphs),
        fFunc(0),
	fIsProcessed(kFALSE)
{
	for(int i=0;i<kNHistos;i++) 
	{
		fHv[i]=0x0;
		fMean[i] = fWidth[i] = 0.0;
	}
        
}

//---------------------------------------------------------------
AliACORDEDataDCS::AliACORDEDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fGraphs("TGraph",kNGraphs),
        fFunc(0),
	fIsProcessed(kFALSE)
{
// Init of class AliACORDEDataDCS
// Gettin the TimeStamp an put it on a string

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", nRun,
	TTimeStamp(startTime).AsString(),
	TTimeStamp(endTime).AsString()));

       
	Init();

}

//---------------------------------------------------------------
AliACORDEDataDCS::~AliACORDEDataDCS() {

	for(int i=0;i<kNHistos;i++) {delete fHv[i]; fHv[i]=0;}
	fGraphs.Clear("C");
	fFunc=0;
}
//---------------------------------------------------------------

AliACORDEDataDCS::AliACORDEDataDCS(const AliACORDEDataDCS & data):
TObject(),
fRun(0),
fStartTime(0),
fEndTime(0),
fGraphs("TGraph",kNGraphs),
fFunc(0),
fIsProcessed(kFALSE)
{
// Setting the initial values
// fRUn, Start of Run, End of Run, IsProcessed

	fRun=data.fRun;
	fStartTime=data.fStartTime;
	fEndTime=data.fEndTime;
	fFunc=data.fFunc;
	fIsProcessed=data.fIsProcessed;


        for(int i=0;i<kNAliases;i++){fAliasNames[i] = data.fAliasNames[i];}

        for(int i=0;i<kNHistos;i++){fHv[i]=data.fHv[i];}



        
}
//--------------------------------------------------------------
AliACORDEDataDCS& AliACORDEDataDCS:: operator=(const AliACORDEDataDCS & data) { 

	
        this->fRun=data.fRun;
	this->fStartTime=data.fStartTime;
	this->fEndTime=data.fEndTime;
	this->fFunc=data.fFunc;
	this->fIsProcessed=data.fIsProcessed;


        for(int i=0;i<kNAliases;i++){this->fAliasNames[i] = data.fAliasNames[i];}

        for(int i=0;i<kNHistos;i++){this->fHv[i]=data.fHv[i];}

	
         return *this;
 
}
//---------------------------------------------------------------
void AliACORDEDataDCS::ProcessData(TMap& aliasMap)
{
// Process of the data from the aliases DCS-data points

	if(!(fHv[0])) Init();

	TObjArray *aliasArr;
	AliDCSValue* aValue;

	for(int j=0; j<kNAliases; j++)
        {
		aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[j].Data());
		if(!aliasArr)
                {
			AliError(Form("Alias %s not found!", fAliasNames[j].Data()));
			continue;
		}
		Introduce(j, aliasArr);

		if(aliasArr->GetEntries()<2)
                {
		AliError(Form("Alias %s has just %d entries!",
					fAliasNames[j].Data(),aliasArr->GetEntries()));
			continue;
		}

		TIter iterarray(aliasArr);

		Double_t *time = new Double_t[aliasArr->GetEntries()];
		Double_t *val  = new Double_t[aliasArr->GetEntries()];

		UInt_t ne=0;

		while ((aValue = (AliDCSValue*) iterarray.Next())) 
                {
                val[ne] = aValue->GetFloat();
		time[ne] = (Double_t) (aValue->GetTimeStamp());
		fHv[j]->Fill(val[ne]);
		ne++;
		}

               
		
		CreateGraph(j, aliasArr->GetEntries(), time, val);
		delete[] val;
		delete[] time;
	}


	// calculate mean and rms of the first two histos
	for(int i=0;i<kNHistos;i++)
        {
		fMean[i] = fHv[i]->GetMean();
		fWidth[i] = fHv[i]->GetRMS();
	}


	fIsProcessed=kTRUE;


}

//---------------------------------------------------------------
void AliACORDEDataDCS::Init()
{
// Init of AliACORDEDatDCS procedure
// Loop over the aliases

	TH1::AddDirectory(kFALSE);

	fGraphs.SetOwner(1);

        TString aliasName;

	for(int i=0;i<kNAliases;i++){

                aliasName.Form("ACO_HV_MODULE%02d_VMON",i); 
		fAliasNames[i] = aliasName;
	}

	for(int i=0;i<kNHistos;i++)
        {
		fHv[i] = new TH1F(fAliasNames[i].Data(),fAliasNames[i].Data(), 20, kHvMin, kHvMax);
		fHv[i]->GetXaxis()->SetTitle("Hv");
	}
}

//---------------------------------------------------------------
void AliACORDEDataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr)
{

	int entries=aliasArr->GetEntries();
	AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
	AliInfo(Form("    	%d DP values collected",entries));

}

//---------------------------------------------------------------
void AliACORDEDataDCS::CreateGraph(int i, int dim, const Double_t *x, const Double_t *y)
{
// Create the plots for the ACORDE DCS 

	TGraph *gr = new(fGraphs[fGraphs.GetEntriesFast()]) TGraph(dim, x, y);

	gr->GetXaxis()->SetTimeDisplay(1);
	gr->SetTitle(fAliasNames[i].Data());

	AliInfo(Form("Array entries: %d",fGraphs.GetEntriesFast()));


}

//---------------------------------------------------------------
void AliACORDEDataDCS::Draw(const Option_t* /*option*/)
{
// Draw all histos and graphs

  if(!fIsProcessed) return;

  TString canvasHistoName;
  TCanvas *ch[10];
  
  for (int i=0;i<10;i++)
  {
  canvasHistoName.Form("ACO_HV_MODULE");
  ch[i]=new TCanvas(canvasHistoName,canvasHistoName,20,20,600,600);
  ch[i]->Divide(2,3);

    for(int j=0;j<6;j++)
    { 
    ch[i]->cd(j+1);
    ((TGraph*) fGraphs.UncheckedAt(i*6+j))->SetMarkerStyle(20);
    ((TGraph*) fGraphs.UncheckedAt(i*6+j))->Draw("alp");
    }

  }

 
}

