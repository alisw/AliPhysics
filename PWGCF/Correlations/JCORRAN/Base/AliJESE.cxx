#include "AliJESE.h"
#include <TH1.h>
#include <TSpline.h>
#include <TFile.h>
#include <TGrid.h>
#include <AliOADBContainer.h>
#include <AliAODEvent.h>
#include <AliAODVertex.h>
#include <AliAODTrack.h>
#include <AliAODHeader.h>
#include <AliAODVZERO.h>

AliJESE::AliJESE(){
	//
}

AliJESE::~AliJESE(){
	//
}

bool AliJESE::Initialize(){
	//
	TGrid::Connect("alien:");
    poadbf = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibV0Run1.root");
	if(!poadbf){
		printf("Unabled to open V0 calibration file.\n");
		return false;
	}

	static const char *pobjNames[OBJECT_COUNT] = {
		"hMultV0BefCorPfpx",
		"fqxa2m",
		"fqxa2s",
		"fqya2m",
		"fqya2s",
		"fqxc2m",
		"fqxc2s",
		"fqyc2m",
		"fqyc2s"
	};
	for(uint i = 0; i < OBJECT_COUNT; ++i){
		poadb[i] = (AliOADBContainer*)poadbf->Get(pobjNames[i]);
		if(!poadb[i])
			return false;
	}

	psplf = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibSpV0CRun1.root");
	if(!psplf){
		printf("Unable to open spline file.\n");
		return false;
	}
	for(uint i = 0; i < 90; ++i)
		psplineQ2c[i] = (TSpline3*)psplf->Get(Form("hqc2Int_%u",i));
	
	return true;
}

void AliJESE::Destroy(){
	poadbf->Close();
	psplf->Close();
}

double AliJESE::Getqc2Perc(AliAODEvent *pevent, float cent, uint runN){
	//
	TH1D *pmultV0 = (TH1D*)poadb[OBJECT_MULTV0]->GetObject(runN);
	if(!pmultV0){
		printf("No mult V0\n");
		return -1.0;
	}

	AliAODVZERO *paodV0 = pevent->GetVZEROData();

	double Qxc2 = 0.0, Qyc2 = 0.0, Nc = 0.0;
	for(uint i = 0; i < 4; ++i){ //V0C
		for(uint j = 0; j < 8; ++j){
			uint iv0 = 8*i+j;
			double mult = paodV0->GetMultiplicity(iv0);
			double phiV0 = TMath::PiOver4()*(0.5+j);//iv0%8);
			double multcorr = mult/pmultV0->GetBinContent(iv0+1)*pmultV0->GetBinContent(8*i+1);
			if(multcorr < 0.0)
				continue;
			Qxc2 += TMath::Cos(2.0*phiV0)*multcorr;
			Qyc2 += TMath::Sin(2.0*phiV0)*multcorr;
			Nc += multcorr;
		}
	}
	if(Nc <= 0.0){
		printf("Nc <= 0\n");
		return -1.0;
	}

	double Qxa2 = 0.0, Qya2 = 0.0, Na = 0.0;
	for(uint i = 4; i < 8; ++i){ //V0A
		for(uint j = 0; j < 8; ++j){
			uint iv0 = 8*i+j;
			double mult = paodV0->GetMultiplicity(iv0);
			double phiV0 = TMath::PiOver4()*(0.5+j);//iv0%8);
			double multcorr = mult/pmultV0->GetBinContent(iv0+1)*pmultV0->GetBinContent(8*i+1);
			if(multcorr < 0.0)
				continue;
			Qxa2 += TMath::Cos(2.0*phiV0)*multcorr;
			Qya2 += TMath::Sin(2.0*phiV0)*multcorr;
			Na += multcorr;
		}
	}
	if(Na <= 0.0){
		printf("Nc <= 0\n");
		return -1.0;
	}
	
	uint centIndex = (uint)TMath::Floor(cent);
	if(centIndex >= 90){
		printf("centIndex > 90\n");
		return -1.0;
	}
	
	/*double Qxa2corr = (Qxa2-((TH1D*)poadb[OBJECT_QXA2M]->GetObject(runN))->GetBinContent(centIndex+1))/((TH1D*)poadb[OBJECT_QXA2S]->GetObject(runN))->GetBinContent(centIndex+1);
	double Qya2corr = (Qya2-((TH1D*)poadb[OBJECT_QYA2M]->GetObject(runN))->GetBinContent(centIndex+1))/((TH1D*)poadb[OBJECT_QYA2S]->GetObject(runN))->GetBinContent(centIndex+1);
	double Qxc2corr = (Qxc2-((TH1D*)poadb[OBJECT_QXC2M]->GetObject(runN))->GetBinContent(centIndex+1))/((TH1D*)poadb[OBJECT_QXC2S]->GetObject(runN))->GetBinContent(centIndex+1);
	double Qyc2corr = (Qyc2-((TH1D*)poadb[OBJECT_QYC2M]->GetObject(runN))->GetBinContent(centIndex+1))/((TH1D*)poadb[OBJECT_QYC2S]->GetObject(runN))->GetBinContent(centIndex+1);*/

	double Qxc2corrSE = Qxc2-((TH1D*)poadb[OBJECT_QXC2M]->GetObject(runN))->GetBinContent(centIndex+1);
	double Qyc2corrSE = Qyc2-((TH1D*)poadb[OBJECT_QYC2M]->GetObject(runN))->GetBinContent(centIndex+1);

	double qc2 = TMath::Sqrt((Qxc2corrSE*Qxc2corrSE+Qyc2corrSE*Qyc2corrSE)/Nc);
	double qc2f = psplineQ2c[centIndex]->Eval(qc2);

	return qc2f;
}

