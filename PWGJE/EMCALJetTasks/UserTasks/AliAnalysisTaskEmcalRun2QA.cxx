/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <cstring>
#include <Riostream.h>

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <THashList.h>

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliTrackContainer.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliPicoTrack.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"

#include "AliAnalysisTaskEmcalRun2QA.h"

ClassImp(AliAnalysisTaskEmcalRun2QA)

//________________________________________________________________________
AliAnalysisTaskEmcalRun2QA::AliAnalysisTaskEmcalRun2QA() :
AliAnalysisTaskEmcalLight("AliAnalysisTaskEmcalRun2QA", kTRUE),
fCellEnergyCut(0.05),
fParticleLevel(kFALSE),
fIsMC(kFALSE),
fCentMethod2(""),
fCentMethod3(""),
fDoV0QA(0),
fDoEPQA(0),
fDoLeadingObjectPosition(0),
fMaxCellsInCluster(50),
fPtBinWidth(0.5),
fMaxPt(150),
fSeparateEMCalDCal(kTRUE),
fIsEmbedded(kFALSE),
fCent2(0),
fCent3(0),
fVZERO(0),
fV0ATotMult(0),
fV0CTotMult(0),
fLeadingTrack(),
//fGeoName("EMCAL_FIRSTYEARV1"),//fGeomEMCal(0),
fCaloUtils(0),

fHistCellIDvsE(0),fHistCellIDvsELow(0),fHistCellEtaPhi(0),
fHistClusterIDvsE(0),fHistClusterIDvsELow(0),fHistClusterEtaPhi(0),fHistNrCellsInCluster(0),

fOutputList_Event(), fOutputList_Cell(), fOutputList_Cluster(),
fHistManager("AliAnalysisTaskEmcalRun2QA")
{
	InitConstants();
}
//________________________________________________________________________
AliAnalysisTaskEmcalRun2QA::AliAnalysisTaskEmcalRun2QA(const char *name):
AliAnalysisTaskEmcalLight(name, kTRUE),
fCellEnergyCut(0.05),
fParticleLevel(kFALSE),
fIsMC(kFALSE),
fCentMethod2(""),
fCentMethod3(""),
fDoV0QA(0),
fDoEPQA(0),
fDoLeadingObjectPosition(0),
fMaxCellsInCluster(50),
fPtBinWidth(0.5),
fMaxPt(150),
fSeparateEMCalDCal(kTRUE),
fIsEmbedded(kFALSE),
fCent2(0),
fCent3(0),
fVZERO(0),
fV0ATotMult(0),
fV0CTotMult(0),
//fGeoName("EMCAL_FIRSTYEARV1"),//fGeomEMCal(0),//
fCaloUtils(0),

fHistCellIDvsE(0),fHistCellIDvsELow(0),fHistCellEtaPhi(0),
fHistClusterIDvsE(0),fHistClusterIDvsELow(0),fHistClusterEtaPhi(0),fHistNrCellsInCluster(0),

fOutputList_Event(), fOutputList_Cell(), fOutputList_Cluster(),
fHistManager(name)
{
	InitConstants();
}
//________________________________________________________________________
void AliAnalysisTaskEmcalRun2QA::InitConstants()
{
	//..Standard
	fRtoD=180.0/TMath::Pi();
	memset(fNTotClusters, 0, sizeof(Int_t)*2);
	fDebug=0; //printout some debug lines

	SetMakeGeneralHistograms(kTRUE);

	//..super ugly but I can not retrieve the info directly from the
	//..AliEMCALGeometry because it is initialized later than the OutputObjects (where it is needed)
	//number of modules, 12 for EMCal + 8 for DCAL
	//2013 fNumOfSuperMod=12;	//HARD CODED eeeeewwwhh!
	fNumOfSuperMod=20;	//2015 HARD CODED eeeeewwwhh!
	if(!fCaloUtils)fCaloUtils = new AliCalorimeterUtils();
}
//________________________________________________________________________
AliAnalysisTaskEmcalRun2QA::~AliAnalysisTaskEmcalRun2QA()
{
   //..Destructor
   if (fCaloUtils) delete fCaloUtils ;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalRun2QA::UserCreateOutputObjects()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskEmcalRun2QA::UserCreateOutputObjects()"<<endl;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   General Event Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

	AliEmcalContainer* cont = 0;

	TString histname;
	TString title;

	Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);

	fOutputList_Event    = new TList();
	fOutputList_Event    ->SetOwner();
	fOutputList_Event    ->SetName("EventQA");

	fOutputList_Cell    = new TList();
	fOutputList_Cell    ->SetOwner();
	fOutputList_Cell    ->SetName("CellQA");

	fOutputList_Cluster    = new TList();
	fOutputList_Cluster    ->SetOwner();
	fOutputList_Cluster    ->SetName("ClusterQA");

	fOutput->Add(fOutputList_Event);
	fOutput->Add(fOutputList_Cell);
	fOutput->Add(fOutputList_Cluster);


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Settings for Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	//..per module!
	Int_t fNMaxCols = AliEMCALGeoParams::fgkEMCALCols; //Eta direction = 48
	Int_t fNMaxRows = AliEMCALGeoParams::fgkEMCALRows; //Phi direction = 24

	TString Calo="EMCal";

	//cout<<"Number of EMCal modules : "<<AliEMCALGeoParams::fgkEMCALModules<<", fNModules: "<<fNModules<<", fGeom: "<< fNumOfSuperMod<<endl;
	//Int_t NumSupMod = fGeom->GetNumberOfSuperModules();  //GetNumberOfSuperModules()
	//cout<<"Number of Supermodules: "<<NumSupMod<<endl;

	if(Calo=="EMCal")
	{
		fNMaxColsAbs = 2*fNMaxCols;  //two supermodules next to each other
		//..as the modules come in rows of 2 the maximum length is half the number of modules x row in a module
		//..what about smaller modules?????
		fNMaxRowsAbs = Int_t(fNumOfSuperMod/2)*fNMaxRows;
	}
	else
	{
		fNMaxColsAbs=fNumOfSuperMod*fNMaxRows;
	}

	cout<<"Info1: "<<endl;
	//number of rows
	Int_t NRows = AliEMCALGeoParams::fgkEMCALRows;
	cout<<"A"<<endl;
	Int_t NCol  = AliEMCALGeoParams::fgkEMCALCols;
	cout<<"B"<<endl;
	Int_t NSMod = AliEMCALGeoParams::fgkEMCALModules;
//    cout<<" Number of super modules :"<<fGeom->GetNumberOfSuperModules()<<", "<<fCaloUtils->GetNumberOfSuperModulesUsed()<<"(fGeom, fCaloUtils)"<<endl;
    cout<<"number of mo from geoparams:  "<<NSMod<<endl;
    cout<<"number of col from geoparams: "<<NCol<<endl;
    cout<<"number of row from geoparams: "<<NRows<<endl;
	cout<<"Number of rows: "<<NRows<<", rowsMax: "<<fNMaxRowsAbs<<", Number of collumns:"<<NCol<<", colMax: "<<fNMaxColsAbs<<endl;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Cluster Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	TIter nextClusColl(&fClusterCollArray);
	while ((cont = static_cast<AliEmcalContainer*>(nextClusColl())))
	{
		fHistManager.CreateHistoGroup(cont->GetArrayName());
		for (Int_t i = 0; i < fNcentBins; i++)
		{
			histname = TString::Format("%s/fHistRejectionReason_%d", cont->GetArrayName().Data(), i);
			title = histname + ";Rejection reason;#it{E}_{cluster} (GeV);counts";
			TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 40, 0, 100);
			SetRejectionReasonLabels(hist->GetXaxis());

			histname = TString::Format("%s/fHistClusPosition_%d", cont->GetArrayName().Data(), i);
			title = histname + ";#it{x} (cm);#it{y} (cm);#it{z} (cm)";
			fHistManager.CreateTH3(histname.Data(), title.Data(), 50, -500, 500, 50, -500, 500, 50, -500, 500);

			histname = TString::Format("%s/fHistClusPhiEtaEnergy_%d", cont->GetArrayName().Data(), i);
			title = histname + ";#eta;#phi;#it{E}_{cluster} (GeV)";
			fHistManager.CreateTH3(histname.Data(), title.Data(), 50, -1, 1, 50, 0, TMath::TwoPi(), nPtBins, 0, fMaxPt);

			//Three different energies
			histname = TString::Format("%s/fHistClusEnergy_%d", cont->GetArrayName().Data(), i);
			title = histname + ";#it{E}_{cluster} (GeV);counts";
			fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

			histname = TString::Format("%s/fHistClusNonLinCorrEnergy_%d", cont->GetArrayName().Data(), i);
			title = histname + ";#it{E}_{cluster} (GeV);counts";
			fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

			histname = TString::Format("%s/fHistClusHadCorrEnergy_%d", cont->GetArrayName().Data(), i);
			title = histname + ";#it{E}_{cluster} (GeV);counts";
			fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
		}

		cout<<"create histograms for cluster name: "<<cont->GetArrayName()<<endl;

		fHistClusterIDvsE    = new TH2*[fNumOfSuperMod+1];
		fHistClusterIDvsELow = new TH2*[fNumOfSuperMod+1];
		fHistClusterEtaPhi   = new TH2*[fNumOfSuperMod+1];
		fHistNrCellsInCluster= new TH1*[fNumOfSuperMod+1];

		for (Int_t i = 0; i < fNumOfSuperMod+1; i++)
		{
			fHistClusterIDvsE[i] = new TH2F(Form("fHistClusterIDvsE_%d",i),Form("fHistClusterIDvsE_%d",i),18000,0,18000,100,0,50);
			if(i==fNumOfSuperMod)fHistClusterIDvsE[i]->GetXaxis()->SetTitle("Leading cell ID (all SuperModules");
			else                  fHistClusterIDvsE[i]->GetXaxis()->SetTitle(Form("Leading cell ID (SuperMod No. %d)",i));
			fHistClusterIDvsE[i]->GetYaxis()->SetTitle("Energy probably in GeV");
			fOutputList_Cluster->Add(fHistClusterIDvsE[i]);

			fHistClusterIDvsELow[i] = new TH2F(Form("fHistClusterIDvsELow_%d",i),Form("fHistClusterIDvsELow_%d",i),18000,0,18000,200,0,2);
			if(i==fNumOfSuperMod)fHistClusterIDvsELow[i]->GetXaxis()->SetTitle("Leading cell ID (all SuperModules");
			else                  fHistClusterIDvsELow[i]->GetXaxis()->SetTitle(Form("Leading cell ID (SuperMod No. %d)",i));
			fHistClusterIDvsELow[i]->GetYaxis()->SetTitle("Energy probably in GeV");
			fOutputList_Cluster->Add(fHistClusterIDvsELow[i]);

/*			//			fHistClusterEtaPhi[i] = new TH2F(Form("fHistClusterEtaPhi_%d",i),Form("fHistClusterEtaPhi_%d",i),1120,60,130,400,-1,1);  //0.012625, 0.0125, 0.0136
			fHistClusterEtaPhi[i] = new TH2F(Form("fHistClusterEtaPhi_%d",i),Form("fHistClusterEtaPhi_%d",i),5600,60,130,400,-1,1);
			if(i==fNumOfSuperMod)fHistClusterEtaPhi[i]->GetXaxis()->SetTitle("Leading cell #phi (all SuperModules");
			else                  fHistClusterEtaPhi[i]->GetXaxis()->SetTitle(Form("Leading cell #phi (SuperMod No. %d)",i));
			fHistClusterEtaPhi[i]->GetYaxis()->SetTitle("Leading cell #eta");
			fOutputList_Cluster->Add(fHistClusterEtaPhi[i]);
*/
            //For row and collumn
			fHistClusterEtaPhi[i] = new TH2F(Form("fHistClusterEtaPhi_%d",i),Form("fHistClusterEtaPhi_%d",i),fNMaxColsAbs+2,-1.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+2,-1.5,fNMaxRowsAbs+0.5);
			if(i==fNumOfSuperMod)fHistClusterEtaPhi[i]->GetXaxis()->SetTitle("Leading cell column (#eta direction)(all SuperModules)");
			else                 fHistClusterEtaPhi[i]->GetXaxis()->SetTitle(Form("Leading cell column (#eta direction)(SuperMod No. %d)",i));
			fHistClusterEtaPhi[i]->GetYaxis()->SetTitle("Leading cell row (#phi direction)");
			fOutputList_Cluster->Add(fHistClusterEtaPhi[i]);


			fHistNrCellsInCluster[i] = new TH1F(Form("fHistNrCellsInCluster_%d",i),Form("fHistNrCellsInCluster_%d",i),30,0,30);
			if(i==fNumOfSuperMod)fHistNrCellsInCluster[i]->GetXaxis()->SetTitle("Number of cells (all SuperModules)");
			else                  fHistNrCellsInCluster[i]->GetXaxis()->SetTitle(Form("Number of cells (SuperMod No. %d)",i));
			fHistNrCellsInCluster[i]->GetYaxis()->SetTitle("N_{Cluster}");
			fOutputList_Cluster->Add(fHistNrCellsInCluster[i]);
		}
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Cells Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if (!fCaloCellsName.IsNull())
	{

		fHistCellIDvsE    = new TH2*[fNumOfSuperMod+1];
		fHistCellIDvsELow = new TH2*[fNumOfSuperMod+1];
		fHistCellEtaPhi   = new TH2*[fNumOfSuperMod+1];

		for (Int_t i = 0; i < fNumOfSuperMod+1; i++)
		{
			fHistCellIDvsE[i] = new TH2F(Form("fHistCellIDvsE_%d",i),Form("fHistCellIDvsE_%d",i),18000,0,18000,100,0,50);
			if(i==fNumOfSuperMod)  fHistCellIDvsE[i]->GetXaxis()->SetTitle("Cell ID (all SuperModules");
			else                   fHistCellIDvsE[i]->GetXaxis()->SetTitle(Form("Cell ID (SuperMod No. %d)",i));
			fHistCellIDvsE[i]->GetYaxis()->SetTitle("Energy probably in GeV");
			fOutputList_Cell->Add(fHistCellIDvsE[i]);

			fHistCellIDvsELow[i] = new TH2F(Form("fHistCellIDvsELow_%d",i),Form("fHistCellIDvsELow_%d",i),18000,0,18000,200,0,2);
			if(i==fNumOfSuperMod)  fHistCellIDvsELow[i]->GetXaxis()->SetTitle("Cell ID (all SuperModules");
			else                   fHistCellIDvsELow[i]->GetXaxis()->SetTitle(Form("Cell ID (SuperMod No. %d)",i));
			fHistCellIDvsELow[i]->GetYaxis()->SetTitle("Energy probably in GeV");
			fOutputList_Cell->Add(fHistCellIDvsELow[i]);

			/*
			 //for eta phi
			if(i==fNumOfSuperMod)fHistCellEtaPhi[i] = new TH2F(Form("fHistCellEtaPhi_%d",i),Form("fHistCellEtaPhi_%d",i),700,60,130,200,-1,1);  //superfine (not necessary)
			else                  fHistCellEtaPhi[i] = new TH2F(Form("fHistCellEtaPhi_%d",i),Form("fHistCellEtaPhi_%d",i),5600,60,130,400,-1,1);  //superfine (not necessary)
			if(i==fNumOfSuperMod)fHistCellEtaPhi[i]->GetXaxis()->SetTitle("Cell #phi (all SuperModules)");
			else                  fHistCellEtaPhi[i]->GetXaxis()->SetTitle(Form("Cell #phi (SuperMod No. %d)",i));
			fHistCellEtaPhi[i]->GetYaxis()->SetTitle("Cell #eta");
			fOutputList_Cell->Add(fHistCellEtaPhi[i]);
			*/

            //For row and collumn
			fHistCellEtaPhi[i] = new TH2F(Form("fHistCellEtaPhi_%d",i),Form("fHistCellEtaPhi_%d",i),fNMaxColsAbs+2,-1.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+2,-1.5,fNMaxRowsAbs+0.5);
			if(i==fNumOfSuperMod)fHistCellEtaPhi[i]->GetXaxis()->SetTitle("column (#eta direction)(all SuperModules)");
			else                 fHistCellEtaPhi[i]->GetXaxis()->SetTitle(Form("column (#eta direction)(SuperMod No. %d)",i));
			fHistCellEtaPhi[i]->GetYaxis()->SetTitle("row (#phi direction)");
			fOutputList_Cell->Add(fHistCellEtaPhi[i]);
		}
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Other Stuff
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Int_t dim = 0;
	TString axistitle[40];
	Int_t nbins[40] = {0};
	Double_t min[40] = {0};
	Double_t max[40] = {0};

	if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp)
	{
		axistitle[dim] = "Centrality %";
		nbins[dim] = 101;
		min[dim] = 0;
		max[dim] = 101;
		dim++;

		if (!fCentMethod2.IsNull())
		{
			axistitle[dim] = Form("Centrality %s %%", fCentMethod2.Data());
			nbins[dim] = 101;
			min[dim] = 0;
			max[dim] = 101;
			dim++;
		}

		if (!fCentMethod3.IsNull())
		{
			axistitle[dim] = Form("Centrality %s %%", fCentMethod3.Data());
			nbins[dim] = 101;
			min[dim] = 0;
			max[dim] = 101;
			dim++;
		}

		if (fDoV0QA==1)
		{
			axistitle[dim] = "V0A total multiplicity";
			nbins[dim] = 200;
			min[dim] = 0;
			max[dim] = 20000;
			dim++;

			axistitle[dim] = "V0C total multiplicity";
			nbins[dim] = 200;
			min[dim] = 0;
			max[dim] = 20000;
			dim++;
		}
		else if (fDoV0QA==2)
		{
			axistitle[dim] = "V0A+V0C total multiplicity";
			nbins[dim] = 300;
			min[dim] = 0;
			max[dim] = 30000;
			dim++;
		}

		if (fDoEPQA)
		{
			axistitle[dim] = "#psi_{EP}";
			nbins[dim] = 200;
			min[dim] = -TMath::Pi();
			max[dim] = TMath::Pi();
			dim++;
		}
	}

	if (fParticleCollArray.GetEntriesFast()>0)
	{
		axistitle[dim] = "No. of tracks";
		if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp)
		{
			nbins[dim] = 3000;
			min[dim] = -0.5;
			max[dim] = 6000-0.5;
		}
		else {
			nbins[dim] = 200;
			min[dim] = 0;
			max[dim] = 200;
		}
		dim++;

		axistitle[dim] = "#it{p}_{T,track}^{leading} (GeV/c)";
		nbins[dim] = nPtBins;
		min[dim] = 0;
		max[dim] = fMaxPt;
		dim++;

		if (fDoLeadingObjectPosition)
		{
			axistitle[dim] = "#eta_{track}^{leading}";
			nbins[dim] = 100;
			min[dim] = -1;
			max[dim] = 1;
			dim++;

			axistitle[dim] = "#phi_{track}^{leading}";
			nbins[dim] = 100;
			min[dim] = 0;
			max[dim] = TMath::TwoPi();
			dim++;
		}
	}

	if (fClusterCollArray.GetEntriesFast()>0)
	{
		axistitle[dim] = "No. of clusters";

		if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp)
		{
			nbins[dim] = 2000;
			min[dim] = -0.5;
			max[dim] = 4000-0.5;
		}
		else
		{
			nbins[dim] = 200;
			min[dim] = 0;
			max[dim] = 200;
		}
		dim++;

		if (fSeparateEMCalDCal)
		{
			axistitle[dim] = "#it{E}_{EMCal cluster}^{leading} (GeV)";
			nbins[dim] = nPtBins;
			min[dim] = 0;
			max[dim] = fMaxPt;
			dim++;

			axistitle[dim] = "#it{E}_{DCal cluster}^{leading} (GeV)";
			nbins[dim] = nPtBins;
			min[dim] = 0;
			max[dim] = fMaxPt;
			dim++;

			if (fDoLeadingObjectPosition)
			{
				axistitle[dim] = "#eta_{EMCal cluster}^{leading}";
				nbins[dim] = 100;
				min[dim] = -1;
				max[dim] = 1;
				dim++;

				axistitle[dim] = "#phi_{EMCal cluster}^{leading}";
				nbins[dim] = 100;
				min[dim] = 0;
				max[dim] = TMath::TwoPi();
				dim++;

				axistitle[dim] = "#eta_{DCal cluster}^{leading}";
				nbins[dim] = 100;
				min[dim] = -1;
				max[dim] = 1;
				dim++;

				axistitle[dim] = "#phi_{DCal cluster}^{leading}";
				nbins[dim] = 100;
				min[dim] = 0;
				max[dim] = TMath::TwoPi();
				dim++;
			}
		}
		else
		{
			axistitle[dim] = "#it{E}_{cluster}^{leading} (GeV)";
			nbins[dim] = nPtBins;
			min[dim] = 0;
			max[dim] = fMaxPt;
			dim++;

			if (fDoLeadingObjectPosition)
			{
				axistitle[dim] = "#eta_{cluster}^{leading}";
				nbins[dim] = 100;
				min[dim] = -1;
				max[dim] = 1;
				dim++;

				axistitle[dim] = "#phi_{cluster}^{leading}";
				nbins[dim] = 100;
				min[dim] = 0;
				max[dim] = TMath::TwoPi();
				dim++;
			}
		}
	}

	if (!fCaloCellsName.IsNull())
	{
		axistitle[dim] = "No. of cells";

		if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp)
		{
			nbins[dim] = 5000;
			min[dim] = -0.5;
			max[dim] = 10000-0.5;
		}
		else
		{
			nbins[dim] = 500;
			min[dim] = -0.5;
			max[dim] = 500-0.5;
		}

		dim++;
	}

	THnSparse* hn = fHistManager.CreateTHnSparse("fHistEventQA","fHistEventQA",dim,nbins,min,max);
	for (Int_t i = 0; i < dim; i++)
		hn->GetAxis(i)->SetTitle(axistitle[i]);

	fOutput->Add(fHistManager.GetListOfHistograms());

	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskEmcalRun2QA::ExecOnce()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskEmcalRun2QA::ExecOnce()"<<endl;

	if (fClusterCollArray.GetEntriesFast() == 0  && fCaloCellsName.IsNull())
	{
		fNeedEmcalGeom = kFALSE;
	}
	else
	{
		fNeedEmcalGeom = kTRUE;
	}

	//produces the fGeom Object, ...
	AliAnalysisTaskEmcalLight::ExecOnce();

	if (fDoV0QA)
	{
		fVZERO = InputEvent()->GetVZEROData();
		if (!fVZERO)
		{
			AliError("AliVVZERO not available");
		}
	}
	//Initialize the calo utils object with the current settings of the event
	if(fCaloUtils)
	{
		cout<<"- - - - - fCaloUtils initialized - - - "<<endl;
		// Set geometry matrices before filling arrays, in case recalibration/position calculation etc is needed
		fCaloUtils->AccessGeometry(InputEvent()); // InputEvent()->GetRunNumber()
		// Set the AODB calibration, bad channels etc. parameters at least once
		fCaloUtils->AccessOADB(InputEvent());
		//apparently not initialized correctly like eg in AliEMCALGeometry!
		fCaloUtils->SetNumberOfSuperModulesUsed(fGeom->GetNumberOfSuperModules());
		cout<<"- - - - - fCaloUtils initialized end- - - "<<endl;
	}
	else
	{
		cout<<"- - - - - no fCaloUtils initialized - - - "<<endl;
	}

	if(fNumOfSuperMod != fGeom->GetNumberOfSuperModules())
	{
		cout<<"Will cause major error in histogram arrays - change hard coded value of fNumOfSuperMod!"<<endl;
		cout<<"Hard coded value is: "<<fNumOfSuperMod<<", value from fGeom is: "<<fGeom->GetNumberOfSuperModules()<<endl;
	}

	//Modules - super modules max rows max cells
	cout<<"Info: "<<endl;
	//number of rows
	Int_t NRows = AliEMCALGeoParams::fgkEMCALRows;
	Int_t NCol  = AliEMCALGeoParams::fgkEMCALCols;
	Int_t NSMod = AliEMCALGeoParams::fgkEMCALModules;
    cout<<" Number of super modules :"<<fGeom->GetNumberOfSuperModules()<<", "<<fCaloUtils->GetNumberOfSuperModulesUsed()<<"(fGeom, fCaloUtils)"<<endl;
    cout<<"number of mo from geoparams:  "<<NSMod<<endl;
    cout<<"number of col from geoparams: "<<NCol<<endl;
    cout<<"number of row from geoparams: "<<NRows<<endl;
	cout<<"Number of rows: "<<NRows<<", rowsMax: "<<fNMaxRowsAbs<<", Number of collumns:"<<NCol<<", colMax: "<<fNMaxColsAbs<<endl;
	//Number of rows: 24, rowsMax: 120, Number of collumns:48, colMax: 96

	//check that! (also for DCal!!!)
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalRun2QA::RetrieveEventObjects()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskEmcalRun2QA::RetrieveEventObjects()"<<endl;

	if (!AliAnalysisTaskEmcalLight::RetrieveEventObjects()) return kFALSE;

	if (!fCentMethod2.IsNull() || !fCentMethod3.IsNull())
	{
		if (fBeamType == kAA || fBeamType == kpA )
		{
			AliCentrality *aliCent = InputEvent()->GetCentrality();
			if (aliCent)
			{
				if (!fCentMethod2.IsNull())
					fCent2 = aliCent->GetCentralityPercentile(fCentMethod2);
				if (!fCentMethod3.IsNull())
					fCent3 = aliCent->GetCentralityPercentile(fCentMethod3);
			}
		}
	}

	if (fVZERO)
	{
		fV0ATotMult = AliESDUtils::GetCorrV0A(fVZERO->GetMTotV0A(),fVertex[2]);
		fV0CTotMult = AliESDUtils::GetCorrV0C(fVZERO->GetMTotV0C(),fVertex[2]);
	}

	return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalRun2QA::FillHistograms()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskEmcalRun2QA::FillHistograms()"<<endl;

	// Fill histograms.
	EventQA_t eventQA;

	DoClusterLoop();
	AliDebug(2,Form("%d clusters found in EMCal and %d in DCal", fNTotClusters[0], fNTotClusters[1]));
	for (Int_t i = 0; i < 2; i++)
	{
		eventQA.fMaxCluster[i] = fLeadingCluster[i];
	}

	if (fCaloCells)
	{
		eventQA.fNCells = DoCellLoop();
		AliDebug(2,Form("%d cells found in the event", eventQA.fNCells));
	}

	eventQA.fCent = fCent;
	eventQA.fCent2 = fCent2;
	eventQA.fCent3 = fCent3;
	eventQA.fV0A = fV0ATotMult;
	eventQA.fV0C = fV0CTotMult;
	eventQA.fEP = fEPV0;
	eventQA.fNClusters[0] = fNTotClusters[0];
	eventQA.fNClusters[1] = fNTotClusters[1];

	FillEventQAHisto(eventQA);

	return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalRun2QA::FillEventQAHisto(const EventQA_t& eventQA)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskEmcalRun2QA::FillEventQAHisto()"<<endl;
	Double_t contents[40]={0};

	Int_t globalNclusters = eventQA.fNClusters[0] + eventQA.fNClusters[1];

	AliTLorentzVector globalMaxCluster = eventQA.fMaxCluster[0].E() > eventQA.fMaxCluster[1].E() ? eventQA.fMaxCluster[0] : eventQA.fMaxCluster[1];

	THnSparse* histEventQA = static_cast<THnSparse*>(fHistManager.FindObject("fHistEventQA"));

	for (Int_t i = 0; i < histEventQA->GetNdimensions(); i++)
	{
		TString title(histEventQA->GetAxis(i)->GetTitle());
		if (title=="Centrality %")
			contents[i] = eventQA.fCent;
		else if (title==Form("Centrality %s %%", fCentMethod2.Data()))
			contents[i] = eventQA.fCent2;
		else if (title==Form("Centrality %s %%", fCentMethod3.Data()))
			contents[i] = eventQA.fCent3;
		else if (title=="V0A total multiplicity")
			contents[i] = eventQA.fV0A;
		else if (title=="V0C total multiplicity")
			contents[i] = eventQA.fV0C;
		else if (title=="V0A+V0C total multiplicity")
			contents[i] = eventQA.fV0A+eventQA.fV0C;
		else if (title=="#psi_{RP}")
			contents[i] = eventQA.fEP;
		else if (title=="No. of clusters")
			contents[i] = globalNclusters;
		else if (title=="No. of cells")
			contents[i] = eventQA.fNCells;
		else if (title=="#it{p}_{T,track}^{leading} (GeV/c)")
			contents[i] = eventQA.fMaxTrack.Pt();
		else if (title=="#eta_{track}^{leading}")
			contents[i] = eventQA.fMaxTrack.Eta();
		else if (title=="#phi_{track}^{leading}")
			contents[i] = eventQA.fMaxTrack.Phi_0_2pi();
		else if (title=="#it{E}_{cluster}^{leading} (GeV)")
			contents[i] = globalMaxCluster.E();
		else if (title=="#eta_{cluster}^{leading}")
			contents[i] = globalMaxCluster.Eta();
		else if (title=="#phi_{cluster}^{leading}")
			contents[i] = globalMaxCluster.Phi();
		else if (title=="#it{E}_{EMCal cluster}^{leading} (GeV)")
			contents[i] = eventQA.fMaxCluster[0].E();
		else if (title=="#eta_{EMCal cluster}^{leading}")
			contents[i] = eventQA.fMaxCluster[0].Phi_0_2pi();
		else if (title=="#phi_{EMCal cluster}^{leading}")
			contents[i] = eventQA.fMaxCluster[0].Eta();
		else if (title=="#it{E}_{DCal cluster}^{leading} (GeV)")
			contents[i] = eventQA.fMaxCluster[1].E();
		else if (title=="#phi_{DCal cluster}^{leading}")
			contents[i] = eventQA.fMaxCluster[1].Phi_0_2pi();
		else if (title=="#eta_{DCal cluster}^{leading}")
			contents[i] = eventQA.fMaxCluster[1].Eta();
		else
			AliWarning(Form("Unable to fill dimension %s!",title.Data()));
	}

	histEventQA->Fill(contents);
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalRun2QA::DoCellLoop()
{
	// Do cell loop.
	if (!fCaloCells) return 0;

	if(fNumOfSuperMod != fGeom->GetNumberOfSuperModules())
	{
		cout<<"Will cause major error in histogram arrays - change hard coded value of fNumOfSuperMod!"<<endl;
		return 0;
	}

	TString histname = TString::Format("%s/fHistCellsAbsIdEnergy_%d", fCaloCellsName.Data(), fCentBin);

	const Int_t ncells = fCaloCells->GetNumberOfCells();
	Int_t nAccCells = 0;
	Double_t CellEta,CellPhi;
	Float_t amp   = -1;
	Int_t absId   = -1;

	//will be set by GetModuleNumberCellIndexesAbsCaloMap
	Int_t icol    = -1;
	Int_t icolAbs = -1;
    Int_t irow    = -1;
    Int_t irowAbs = -1;
    Int_t iRCU    = -1;

	for (Int_t pos = 0; pos < ncells; pos++)
	{
		amp   = fCaloCells->GetAmplitude(pos);
		if (amp < fCellEnergyCut) continue;
        absId = fCaloCells->GetCellNumber(pos);
		fGeom->EtaPhiFromIndex(absId,CellEta,CellPhi);
		Int_t smNo = fGeom->GetSuperModuleNumber(absId);

		//get the cell row and collumn
		//enum detector { kEMCAL = 0, kPHOS = 1, kCTS = 2, kDCAL = 3, kDCALPHOS = 4 };
		fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(absId,0,icol,irow,iRCU,icolAbs,irowAbs);
		if(icolAbs> fNMaxColsAbs || irowAbs>fNMaxRowsAbs)
		{
			cout<<" Problem! wrong calculated number of max col and max rows"<<endl;
			cout<<"current col: "<<icolAbs<<", max col"<<fNMaxColsAbs<<endl;
			cout<<"current row: "<<irowAbs<<", max row"<<fNMaxRowsAbs<<endl;
		}

		fHistCellIDvsE[smNo]   ->Fill(absId,amp);
		fHistCellIDvsELow[smNo]->Fill(absId,amp);
		fHistCellEtaPhi[smNo]  ->Fill(icolAbs,irowAbs); //icolAbs,irowAbs

		//all centralities
		fHistCellIDvsE[fNumOfSuperMod]   ->Fill(absId,amp);
		fHistCellIDvsELow[fNumOfSuperMod]->Fill(absId,amp);
		fHistCellEtaPhi[fNumOfSuperMod]  ->Fill(icolAbs,irowAbs);


		//IsInDCAL
		//IsInEMCAL

		nAccCells++;
	}
	return nAccCells;
}
//________________________________________________________________________
void AliAnalysisTaskEmcalRun2QA::DoClusterLoop()
{
	// Do cluster loop.
	TString histname;

	memset(fNTotClusters, 0, sizeof(Int_t)*2);
	for (Int_t i = 0; i < 2; i++) fLeadingCluster[i].SetPxPyPzE(0,0,0,0);

	AliClusterContainer* clusters = 0;
	TIter nextClusColl(&fClusterCollArray);
	while ((clusters = static_cast<AliClusterContainer*>(nextClusColl())))
	{
		// Cluster loop
		AliClusterIterableMomentumContainer itcont = clusters->all_momentum();
		//will be set by GetModuleNumberCellIndexesAbsCaloMap
		Int_t icol    = -1;
		Int_t icolAbs = -1;
	    Int_t irow    = -1;
	    Int_t irowAbs = -1;
	    Int_t iRCU    = -1;

		for (AliClusterIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++)
		{
			UInt_t rejectionReason = 0;
//doesn't work			if (!clusters->AcceptCluster(clusters->GetCurrentID(), rejectionReason))
//doesn't work			if (!clusters->AcceptCluster(it->second->GetID(), rejectionReason))
			if (!clusters->AcceptCluster(it->second, rejectionReason))
			{
				histname = TString::Format("%s/fHistRejectionReason_%d", clusters->GetArrayName().Data(), fCentBin);
				fHistManager.FillTH2(histname, clusters->GetRejectionReasonBitPosition(rejectionReason), it->first.E());
				continue;
			}

			Float_t pos[3]={0};
			it->second->GetPosition(pos);

			//- - - - - - - - - - - - - - - - - - - - - - - - - -
			Int_t    NCells = it->second->GetNCells();
			Double_t FirstCellEnergy=0;
			Double_t FirstCellID=0;
			Double_t CellEta,CellPhi;
			Double_t ClusterRawEnergy = it->second->E();

			for (Int_t iCell = 0; iCell < NCells; iCell++)
			{
				Int_t    cellId   = it->second->GetCellAbsId(iCell);
				Double_t cellFrac = it->second->GetCellAmplitudeFraction(iCell);
				Double_t Eseed    = fCaloCells->GetCellAmplitude(cellId);
				if(iCell==0)
				{
					FirstCellEnergy=Eseed;
					FirstCellID    =cellId;
				}
				if(Eseed>FirstCellEnergy)cout<<"Problem: Energy bigger than first."<<iCell<<endl;
			}
			//use the id of the first cell (should be the leading cell)
			Int_t smNo = fGeom->GetSuperModuleNumber(FirstCellID);
			fGeom->EtaPhiFromIndex(FirstCellID,CellEta,CellPhi);
			//get the cell row and collumn
			//enum detector { kEMCAL = 0, kPHOS = 1, kCTS = 2, kDCAL = 3, kDCALPHOS = 4 };
			fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(FirstCellID,0,icol,irow,iRCU,icolAbs,irowAbs);

			//number od cells in cluster vs module (Maybe eta-phi raender weniger etc)
			fHistClusterIDvsE[smNo]   ->Fill(FirstCellID,ClusterRawEnergy);
			fHistClusterIDvsELow[smNo]->Fill(FirstCellID,ClusterRawEnergy);
		    fHistClusterEtaPhi[smNo]  ->Fill(icolAbs,irowAbs);

			//all centralities/modules
			fHistClusterIDvsE[fNumOfSuperMod]   ->Fill(FirstCellID,ClusterRawEnergy);
			fHistClusterIDvsELow[fNumOfSuperMod]->Fill(FirstCellID,ClusterRawEnergy);
			fHistClusterEtaPhi[fNumOfSuperMod]  ->Fill(icolAbs,irowAbs);

			//one dimensional // would be nice to have two dimensional.q
			fHistNrCellsInCluster[smNo]  ->Fill(NCells);
			//- - - - - - - - - - - - - - - - - - - - - - - - - -


			histname = TString::Format("%s/fHistClusPosition_%d", clusters->GetArrayName().Data(), fCentBin);
			fHistManager.FillTH3(histname, pos[0], pos[1], pos[2]);

			Double_t phi = it->first.Phi_0_2pi();

			Int_t isDcal = Int_t(phi > fgkEMCalDCalPhiDivide);

			histname = TString::Format("%s/fHistClusPhiEtaEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
			fHistManager.FillTH3(histname, it->first.Eta(), it->first.Phi_0_2pi(), it->first.E());

			if (fLeadingCluster[isDcal].E() < it->first.E()) fLeadingCluster[isDcal] = it->first;


			histname = TString::Format("%s/fHistClusMCEnergyFraction_%d", clusters->GetArrayName().Data(), fCentBin);
			if (fHistManager.FindObject(histname))
			{
				fHistManager.FillTH1(histname, it->second->GetMCEnergyFraction());
			}

			histname = TString::Format("%s/fHistClusEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
			fHistManager.FillTH1(histname, it->second->E());

			if (it->second->GetNonLinCorrEnergy() > 0.)
			{
				histname = TString::Format("%s/fHistClusNonLinCorrEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
				fHistManager.FillTH1(histname, it->second->GetNonLinCorrEnergy());
			}

			if (it->second->GetHadCorrEnergy() > 0.)
			{
				histname = TString::Format("%s/fHistClusHadCorrEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
				fHistManager.FillTH1(histname, it->second->GetHadCorrEnergy());
			}

			fNTotClusters[isDcal]++;
		}
	}
}
