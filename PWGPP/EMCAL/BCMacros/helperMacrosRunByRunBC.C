/// \file helperMacrosRunByRunBC.C
/// \ingroup EMCALOfflineMacros
/// \brief This macro is used for a run-by-run evaluation of the bad channels and finding the best period splitting for long runlists
///
///
/// See https://twiki.cern.ch/twiki/bin/view/ALICE/BadChannelAnalysis for general documentation
///
/// ---------------------
///  Running the macro
/// ---------------------
/// use root -b to speed up (no canvas drawn)                              <br>
/// root [0] .L helperMacrosRunByRunBC.C++                                 <br>
/// root [2] SummarizeRunByRun("LHC15o","Train_771","INT7","runList45",45) <br>
/// root [2] SummarizeRunByRun("LHC15o","Train_771","INT7","GloballyGood") <br>
/// root [2] GetBestPeriodSplitting("LHC15o",771,105,3)                    <br>
/// root [2] CompareTwoBCstrategies(TString period="LHC15n",Int_t trainNo=603,Int_t version=5)  <br>
///
/// \author Eliane Epple <eliane.epple@yale.edu>, Yale University
/// \date June 29, 2017

// --- ROOT system ---
#include <Riostream.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TEnv.h>
#include <TSystem.h>


// --- ANALYSIS system ---
#include "AliEMCALGeometry.h"          //include when compile
#include "AliCalorimeterUtils.h"       //include when compile
#include "AliAODEvent.h"               //include when compile

//colors
const Int_t RainbowColors[]= {kRed, kRed-4, kRed-7, kRed-9, kRed-10, kYellow, kYellow-4, kYellow-7, kYellow-9, kYellow-10, kGreen, kGreen-4 , kGreen-7, kGreen-9, kGreen-10, kCyan, kCyan-4, kCyan-7, kCyan-9, kCyan-10, kBlue, kBlue-4, kBlue-7, kBlue-9, kBlue-10, kMagenta, kMagenta-4, kMagenta-7, kMagenta-9, kMagenta-10};

//definition of methods
TH2F* CompressHistogram(TH2 *Histo,Int_t totalCells, Int_t badCells,std::vector<Int_t> runIdVec);

void BuildMaxMinHisto(TH1D* inHisto, TH1D* minHist,TH1D* maxHist);
void PlotLowFractionCells(TString pdfName, std::vector<Int_t> cellVector,TH2F* badVsCell[],Int_t nRuns,TH2F* ampID[],TH1D* hCellGoodMean[]);
Bool_t IsItReallyBadRatio(TH1D* minHistoRatio,TH1D* maxHistoRatio,TH1D* meanHistoRatior,TString& crit);

void PlotHorLineRange(Double_t y_val, Double_t xLow, Double_t xHigh, Int_t Line_Col);
void SetHisto(TH2 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto);
void SetHisto(TH1 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto);
Bool_t IsCellMaskedByHand(Int_t cell, std::vector<Int_t> cellVector);
void CreateCellCompPDF(TH2F* hAmpIDMasked, std::vector<Int_t> cellVector, TH1* goodCellsMerged, TH1* goodCellsRbR, TString pdfName);
void Plot2DCells(TString Block, Int_t runNo, std::vector<Int_t> cellVectorRbR, std::vector<Int_t> cellVectorMerge);

/// Draw the good, bad, dead channel maps, and the amplitude distribution per each run and save them in pdf files in analysisOutput/train/RunByRunSummary
/// -- can improve: write different pages in one pdf
/// -- decide how to treat this info
//________________________________________________________________________
void SummarizeRunByRun(TString period = "LHC15o", TString train = "Train_641", TString trigger= "AnyINTnoBC", TString listName="runList",Int_t runsUsed=-1)
{
	gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //..to supress a lot of standard output

	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	//gStyle->SetPalette(53);  //standard is 1
	gStyle->SetCanvasColor(10);
	TGaxis::SetMaxDigits(4);
	gStyle->SetPadTopMargin(0.07);//0.05
	gStyle->SetPadBottomMargin(0.18);//0.15
	gStyle->SetPadRightMargin(0.10);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetFrameFillColor(10);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetTitleSize(5.0,"X");
	gStyle->SetTitleSize(5.0,"Y");
	gEnv->SetValue("Canvas.ShowEventStatus",1);  //shows the status bar in the canvas

	//..............................................
	//..manually disable cells
	std::vector<Int_t> badcellsBlock1;
	std::vector<Int_t> badcellsBlock2;
	std::vector<Int_t> badcellsBlock3;
	std::vector<Int_t> badcellsBlock4;
	//badcellsBlock1.push_back(13483);
	//badcellsBlock2.push_back(13483);
	//badcellsBlock3.push_back(13483);
	//badcellsBlock4.push_back(13483);

	//..select runs after which a new bad map is built
	//..you get these numbers after you run this function and then the GetBestPeriodSplitting function
	Int_t splitRuns1=34;  //run bock is inclusive of this run
	Int_t splitRuns2=66;  //run bock is inclusive of this run
	Int_t splitRuns3=74;  //run bock is inclusive of this run
	//..............................................
	TString analysisInput  = Form("AnalysisInput/%s",period.Data());
	TString analysisOutput = Form("AnalysisOutput/%s/%s",period.Data(),train.Data());
	TString runList        = Form("./%s/%s/%s.txt", analysisInput.Data(), train.Data(),listName.Data());
	gSystem->mkdir(TString::Format("%s/RunByRunSummary%i/", analysisOutput.Data(),runsUsed));

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..open the text file and save the run IDs into the RunId[] array
	cout<<"o o o Open .txt file with run indices. Name = " << runList << endl;
	FILE *pFile = fopen(runList.Data(), "r");
	if(!pFile)
	{
		cout<<"couldn't open file "<<runList<<"!"<<endl;
		return;
	}
	Int_t q;
	Int_t ncols;
	Int_t nlines = 0 ;
	Int_t RunId[500] ;
	std::vector<Int_t> RunIdVec;
	while (1)
	{
		ncols = fscanf(pFile,"  %d ",&q);
		if (ncols< 0) break;
		RunId[nlines]=q;
		RunIdVec.push_back(q);
		nlines++;
	}
	fclose(pFile);
	//..sort the vector by size to be shure to use the right order
	std::sort (RunIdVec.begin(), RunIdVec.end());
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Create different canvases for run-by-runcomparision
	cout<<"o o o Found " << RunIdVec.size() <<" files in list"<< endl;
	Int_t intRun;

	intRun= RunIdVec.size();
	if(runsUsed>0 && runsUsed<intRun)intRun = runsUsed; //for test purposes
	const Int_t nRun = intRun;
	Int_t nRunsUsed = nRun;
	//ELI for MartinInt_t totalperCv = 4;
	Int_t totalperCv = 16;
	Int_t nPad = TMath::Sqrt(totalperCv);
	Int_t nCv = nRun/totalperCv+1;


	if(nCv<1)nCv=1;

	//..canvases per run
	TCanvas **cBad  = new TCanvas*[nCv];
	TCanvas **cGood = new TCanvas*[nCv];
	TCanvas **cDead = new TCanvas*[nCv];
	TCanvas **cAmp  = new TCanvas*[nCv];

	for(Int_t ic = 0; ic<nCv; ic++)
	{
		cBad [ic] = new TCanvas(TString::Format("badcells%d", ic), TString::Format("I) badcells  (%d/%d)", ic+1, nCv), 1000,750);
		cGood[ic] = new TCanvas(TString::Format("goodcells%d", ic),TString::Format("I) goodcells (%d/%d)", ic+1, nCv),1000,750);
		cDead[ic] = new TCanvas(TString::Format("deadcells%d", ic),TString::Format("I) deadcells (%d/%d)", ic+1, nCv),1000,750);
		cAmp [ic] = new TCanvas(TString::Format("Amplitide%d", ic),TString::Format("I) Amplitide (%d/%d)", ic+1, nCv),1000,750);

		cBad [ic] ->Divide(nPad,nPad,0.001,0.001);
		cGood[ic] ->Divide(nPad,nPad,0.001,0.001);
		cDead[ic] ->Divide(nPad,nPad,0.001,0.001);
		cAmp [ic] ->Divide(nPad,nPad,0.001,0.001);
	}

	//..summary figures for all runs
	Int_t nFlags = 3;
	TH2F** hFlagvsRun = new TH2F*[nFlags]; 
	TH2F* hFlagNew  = 0x0;
	TH2F* hFlagNewClean  = 0x0;
	TH2F** ampID         = new TH2F*[nRun];
	TH2F** ampIDCl       = new TH2F*[nRun];
	TH2F** ampIDCl3Block = new TH2F*[nRun];
	TH2F** ampIDCl1Block = new TH2F*[nRun];
	TH2F** ampIDDelete   = new TH2F*[nRun];
	TH1D** hCellGoodMean = new TH1D*[nRun];
	TH1D** hNEvent       = new TH1D*[nRun];
	TH2F* hBadVsEvent      = new TH2F("hBadVsEvent","hBadVsEvent",100,100000,25000000,60,700,1800);
	TH2F* hDeadBadVsEvent  = new TH2F("hDeadBadVsEvent","hDeadBadVsEvent",100,100000,25000000,60,700,1800);
	TH1D* deadbadCellsVsRun;
	TH1D* deadCellsVsRun;
	TH1D* badCellsVsRun;
	TH1D* deadCellsVsRunC;
	TH1D* badCellsVsRunC;
	TH1D* projSum;
	TH1D* projSumC;
	TH1D* projSumC3Blocks;
	TH1D* projSumC3BlocksA;
	TH1D* projSumC3BlocksB;
	TH1D* projSumC3BlocksC;
	TH1D* projSumC3BlocksD;
	TH1D* projSumC1Block;
	TH1D* nEventsVsRuns;
	TH2F* Sum2DSingleMask;
	TH2F* Sum2D3BlockMask;
	TH2F* Sum2D3BlockMaskA;
	TH2F* Sum2D3BlockMaskB;
	TH2F* Sum2D3BlockMaskC;
	TH2F* Sum2D3BlockMaskD;
	TH2F* Sum2DOrig;
	TH2F* Sum2DIdeal;
	TH1D* hgoodMean;

	for(Int_t i = 0; i<nFlags; i++)
	{
		hFlagvsRun[i] = 0x0;
	}
	TCanvas *cFlagDeadBadI      = new TCanvas("cFlagDeadBadI", "II) Flag dead or bad a", 1600, 1000);
	TCanvas *cFlagDeadBadII     = new TCanvas("cFlagDeadBadII", "II) Flag dead or bad b", 1600, 1000);
	TCanvas *cFlagSumI          = new TCanvas("cFlagSumI", "II) Flag dead&bad a", 1600, 1000);
	TCanvas *cFlagSumII         = new TCanvas("cFlagSumII", "II) Flag dead&bad b", 1600, 1000);
	TCanvas *cFlagSumCleanedI   = new TCanvas("cFlagSumCleanI", "III) cleanded Flag dead&bad a", 1600, 1000);
	TCanvas *cFlagSumCleanedII  = new TCanvas("cFlagSumCleanII", "III) cleanded Flag dead&bad b", 1600, 1000);
	TCanvas *cFlagSumCompAllI   = new TCanvas("cFlagSumCompAllI", "III) compressed Flag dead&bad a", 1600, 1000);
	TCanvas *cFlagSumCompAllII  = new TCanvas("cFlagSumCompAllII", "III) compressed Flag dead&bad b", 1600, 1000);
	TCanvas *cFlagSumCompCleanI = new TCanvas("cFlagSumCompI", "III) compressed&cleaned Flag dead&bad a", 1600, 1000);
	TCanvas *cFlagSumCompCleanII= new TCanvas("cFlagSumCompII", "III) compressed&cleaned Flag dead&bad b", 1600, 1000);

	TCanvas *cFlagNew         = new TCanvas("cFlagNew", "IV) Frac dead&bad 2D", 1600, 800);
	TCanvas *cellSummaryCan   = new TCanvas("cellSummaryCan", "I) run overview", 1600, 800);
	TCanvas *cellSummaryCan2  = new TCanvas("cellSummaryCan2","I) run overview II",1600,800);
	TCanvas *cAmpSum          = new TCanvas("SumOfAmplitudes","I) Sum of Amplitides over runs",1500,750);
	TCanvas *cAmpSum2D        = new TCanvas("SumOf2DAmplitudes","I) Sum of 2D Amplitides over runs",1500,750);
	TCanvas *cAmpSum2D4Blocks = new TCanvas("cAmpSum2D4Blocks","I) Sum of 2D Amplitides in 4 run blocks",1500,750);
	TCanvas *cAmpSum2D4BlocksRatio = new TCanvas("cAmpSum2D4BlocksRatio","I) Ratio of 2D Amplitides in 4 run blocks",1500,750);

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..defining variables and histograms
	TString rootFileName;
    TString badChannelOutput;
	/*TString badChannelOutput[4];
	badChannelOutput[0]="AnalysisOutput/LHC16h/Version2/Train_622AnyINTnoBC_Histograms_V2.root";
	badChannelOutput[1]="AnalysisOutput/LHC16i/Version2/Train_623AnyINTnoBC_Histograms_V2.root";
	badChannelOutput[2]="AnalysisOutput/LHC16k/Version0/Train_658AnyINTnoBC_Histograms_V0.root";
	badChannelOutput[3]="AnalysisOutput/LHC16o/Version3/Train_663AnyINTnoBC_Histograms_V3.root";
    */
	Int_t   noOfCells=0;
	Int_t usedRuns=0;
	Bool_t iFirstRun=0;

	AliCalorimeterUtils *fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(RunIdVec.at(0));
	fCaloUtils->AccessGeometry(aod);
    noOfCells=fCaloUtils->GetEMCALGeometry()->GetNCells(); //..Very important number, never change after that point!

	hFlagvsRun[0] = new TH2F("hFlag1vsRun", "hFlag1vsRun (?); cell ID; Run number", noOfCells+10, 0, noOfCells+10, nRun, 0, nRun); // update this axis, need to have the run number
	hFlagvsRun[1] = new TH2F("hFlag2vsRun", "hFlag2vsRun (?); cell ID; Run number", noOfCells+10, 0, noOfCells+10, nRun, 0, nRun);
	hFlagvsRun[2] = new TH2F("hFlag3vsRun", "hFlag3vsRun (?); cell ID; Run number", noOfCells+10, 0, noOfCells+10, nRun, 0, nRun);
	nEventsVsRuns = new TH1D("nEventVsRun", "number of events in run", nRun, 0, nRun);

	for(Int_t i=0;i<nRun;i++)
	{
		hFlagvsRun[0]->GetYaxis()->SetBinLabel(i+1,Form("%i",RunIdVec.at(i)));
		hFlagvsRun[1]->GetYaxis()->SetBinLabel(i+1,Form("%i",RunIdVec.at(i)));
		hFlagvsRun[2]->GetYaxis()->SetBinLabel(i+1,Form("%i",RunIdVec.at(i)));
	}
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Open the different .root files with help of the run numbers from the text file
	//..loop over the amount of run numbers found in the previous text file.
	for(Int_t i = 0 ; i < nRun ; i++)  //Version%i" LHC16i_muon_caloLego_Histograms_V255539
	{
		rootFileName      = Form("%s_Histograms_V%i.root", trigger.Data(), RunIdVec.at(i));
		badChannelOutput  = Form("%s/Version%i/%s", analysisOutput.Data(), RunIdVec.at(i),rootFileName.Data());

		//Martin cout<<"Open root file: "<<badChannelOutput[i]<<endl;
		cout<<"Open root file No: "<<i+1<<" - "<<badChannelOutput<<" - "<<flush;
		TFile *f = TFile::Open(badChannelOutput);
		if(!f)
		{
			cout<<"Couldn't open/find .root file: "<<badChannelOutput<<endl;
			cout<<endl;
			break;
		}

		//..you may introduce a cut here, selecting only
		//..runs with a certain number of events
		hNEvent[usedRuns]      = (TH1D*)f->Get("hNEvents");
		cout<<hNEvent[usedRuns]->Integral()<<" evt."<<endl;
//		if(hNEvent[usedRuns]->Integral()>1000000)continue;

		TH2F *goodCells        = (TH2F*)f->Get("2DChannelMap_Flag0");
		TH2F *deadCells        = (TH2F*)f->Get("2DChannelMap_Flag1");
		TH2F *badCells         = (TH2F*)f->Get("2DChannelMap_Flag2");
		TH1F *hCellFlag        = (TH1F*)f->Get("fhCellFlag");
		ampID[usedRuns]        = (TH2F*)f->Get("hCellAmplitude");
		hCellGoodMean[usedRuns]= (TH1D*)f->Get("hgoodMean");

		if(!badCells || !goodCells || !deadCells || !ampID[usedRuns] || !hCellFlag)
		{
			if(!badCells) Printf("2DChannelMap_Flag2 not found");
			if(!goodCells)Printf("2DChannelMap_Flag0 not found");
			if(!deadCells)Printf("2DChannelMap_Flag1 not found");
			if(!ampID[usedRuns]) Printf("hCellAmplitude not found");
			if(!hCellFlag)Printf("fhCellFlag not found");
			cout<<endl;
			continue;
		}

		nEventsVsRuns->SetBinContent(usedRuns+1,hNEvent[usedRuns]->Integral());
		ampID[usedRuns]        ->SetName(Form("hCellAmplitudeRun%i",usedRuns));
		ampIDCl[usedRuns]      = (TH2F*)ampID[usedRuns]->Clone(Form("hCellAmplitudeClRun%i",usedRuns));
		ampIDCl3Block[usedRuns]= (TH2F*)ampID[usedRuns]->Clone(Form("hCellAmplitudeCl3RunBlock%i",usedRuns));
		ampIDCl1Block[usedRuns]= (TH2F*)ampID[usedRuns]->Clone(Form("hCellAmplitudeCl1RunBlock%i",usedRuns));
		ampIDDelete[usedRuns]  = (TH2F*)ampID[usedRuns]->Clone(Form("hCellAmplitudeTest%i",usedRuns));

		hCellGoodMean[usedRuns]->SetLineColor(kGreen);
		if(iFirstRun==0)
		{
			Sum2DOrig      = (TH2F*)ampID[usedRuns]->Clone("Sum2DOrig");
			Sum2DSingleMask= (TH2F*)ampID[usedRuns]->Clone("Sum2DSingleMask");
			Sum2D3BlockMask= (TH2F*)ampID[usedRuns]->Clone("Sum2D3BlockMask");
			Sum2DIdeal     = (TH2F*)ampID[usedRuns]->Clone("Sum2DIdealDistr");
			Sum2DOrig      ->Reset();
			Sum2DSingleMask->Reset();
			Sum2D3BlockMask->Reset();
			Sum2DIdeal     ->Reset();
			hgoodMean      = (TH1D*)hCellGoodMean[usedRuns]->Clone("MeanSpectrumAllRuns");
			hgoodMean      ->Reset();
			Sum2D3BlockMaskA =(TH2F*)Sum2D3BlockMask->Clone("Sum2D3BlockMaskA");
			Sum2D3BlockMaskB =(TH2F*)Sum2D3BlockMask->Clone("Sum2D3BlockMaskB");
			Sum2D3BlockMaskC =(TH2F*)Sum2D3BlockMask->Clone("Sum2D3BlockMaskC");
			Sum2D3BlockMaskD =(TH2F*)Sum2D3BlockMask->Clone("Sum2D3BlockMaskD");
		}
		//if(i<30)hCellGoodMean[i]->SetLineColor(RainbowColors[i]);
		hgoodMean->Add(hCellGoodMean[usedRuns]); //..add all good distributions to build a mean of all runs
		hCellGoodMean[usedRuns]->SetLineWidth(3);

		//..fill the histo bad cell vs. run number
		Int_t percBad=0;
		for(Int_t icell = 0; icell < noOfCells; icell++)
		{
			Int_t flag =  hCellFlag->GetBinContent(icell+1);
			//..dead
			if(flag == 1) hFlagvsRun[0]->Fill(icell, usedRuns, 1);
			//..bad or warm
			if(flag>1)    hFlagvsRun[1]->Fill(icell, usedRuns, 1); //fill, use the x, y values
			//..dead+bad
			if(flag>0)
			{
				hFlagvsRun[2]->Fill(icell, usedRuns, 1);
				percBad++;
			}
		}
		if(1.0*percBad/noOfCells>0.3)cout<<"Problem in this run detected. Large number of bad+dead cells (>30%) - please double check!"<<endl;

		if(!hFlagNew)
		{
			hFlagNew = (TH2F*)goodCells->Clone(TString::Format("h2DChannelMapNew_FlagAll"));
			hFlagNew->Reset();
			hFlagNew->SetTitle("Selected flag greater than 0; cell column (#eta direction); cell raw (#phi direction)");
			hFlagNewClean = (TH2F*)goodCells->Clone(TString::Format("h2DChannelMapNew_FlagAllClean"));
			hFlagNewClean->Reset();
			hFlagNewClean->SetTitle("Selected flag greater than 0; cell column (#eta direction); cell raw (#phi direction)");
		}

		// Drawing histograms for each run
		//....................................
		cBad[usedRuns/totalperCv]->cd(usedRuns%totalperCv+1);
		SetHisto(badCells,"","",0);
		badCells->Draw("colz");
		TLatex* text = new TLatex(0.2,0.85,Form("Bad Cells - Run %i",RunIdVec.at(i)));
		text->SetTextSize(0.06);
		text->SetNDC();
		text->SetTextColor(1);
		text->Draw();
		//....................................
		cGood[usedRuns/totalperCv]->cd(usedRuns%totalperCv+1);
		SetHisto(goodCells,"","",0);
		goodCells->Draw("colz");
		TLatex* text1 = new TLatex(0.2,0.85,Form("Good Cells - Run %i",RunIdVec.at(i)));
		text1->SetTextSize(0.06);
		text1->SetNDC();
		text1->SetTextColor(1);
		text1->Draw();
		//....................................
		cDead[usedRuns/totalperCv]->cd(usedRuns%totalperCv+1);
		SetHisto(deadCells,"","",0);
		deadCells->Draw("colz");
		TLatex* text2 = new TLatex(0.2,0.85,Form("Dead Cells - Run %i",RunIdVec.at(i)));
		text2->SetTextSize(0.06);
		text2->SetNDC();
		text2->SetTextColor(1);
		text2->Draw();
		//....................................
		cAmp[usedRuns/totalperCv]->cd(usedRuns%totalperCv+1)->SetLogy();
		TH1D* proj = ampID[usedRuns]->ProjectionX(TString::Format("hampIDProj_Run%d",RunIdVec.at(i)));
		SetHisto(proj,"","",0);
		proj->SetLineColor(kCyan+2);
		proj->GetYaxis()->SetRangeUser(0.0000001,100);
		proj->Draw("hist");
		TLatex* text3 = new TLatex(0.2,0.85,Form("Amplitudes - Run %i",RunIdVec.at(i)));
		text3->SetTextSize(0.06);
		text3->SetNDC();
		text3->SetTextColor(1);
		text3->Draw();
		//..create a summ version
		if(iFirstRun==0)projSum = ampID[usedRuns]->ProjectionX("hampIDProj_Sum");
		if(iFirstRun>0) projSum->Add(proj);
		Sum2DOrig->Add(ampID[usedRuns]);

		iFirstRun=1;
		usedRuns++;
	}
	nRunsUsed=usedRuns;  //why that?? I forgot
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	// Count number of cells that are bad in at least one run
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	Int_t nBadCells=0;
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		TH1D *htmpCell = hFlagvsRun[2]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		Double_t sumRun = 0;
		for(Int_t ir = 0; ir < nRunsUsed ; ir++)
		{
			sumRun += htmpCell->GetBinContent(ir+1);
			if(htmpCell->GetBinContent(ir+1)==1)
			{
				//..mask the bad cells for this run
				//..Direction of amplitude (Checks energies from 0-nBins GeV)
				for (Int_t amp = 1; amp <= ampIDCl[ir]->GetNbinsX(); amp++)
				{
					ampIDCl[ir]->SetBinContent(amp,ic+1,0);
				}
			}
		}
		if(sumRun!=0)nBadCells++; //only count for the dead+bad case
	}

	hgoodMean->Scale(1.0/nRunsUsed);

    //..create an ideal 2D energy distribution for a later division
	//..helps to identify where cells have been unmasked and whether
	//..this was a good or bad unmasking desicion (e.g. creating spikes)
	for(Int_t eBin=0;eBin<Sum2DIdeal->GetNbinsX();eBin++)
	{
		Double_t binVal=hgoodMean->GetBinContent(eBin+1);
		for(Int_t icell=0;icell<Sum2DIdeal->GetNbinsY();icell++)
		{
			Sum2DIdeal->SetBinContent(eBin+1,icell+1,binVal);
		}
	}
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	// Draw masked cell amplitude by masking cells that were identified bad or dead in this specific run
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	for(Int_t ir = 0; ir < nRunsUsed ; ir++)
	{
		cAmp[ir/totalperCv]->cd(ir%totalperCv+1)->SetLogy();
		TH1D* proj = ampIDCl[ir]->ProjectionX(TString::Format("hampIDMaskedProj_Run%d",RunIdVec.at(ir)));
		SetHisto(proj,"","",0);
		proj->SetLineColor(kSpring-2);
		proj->Draw("hist same");
		//..create a sum version
		if(ir==0)projSumC = ampIDCl[ir]->ProjectionX("hampIDProj_SumMasked");
		if(ir>0) projSumC->Add(proj);
		Sum2DSingleMask->Add(ampIDCl[ir]);
	}
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	// Draw summary histogram with dead and bad cells vs run
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	deadCellsVsRun    =hFlagvsRun[0]->ProjectionY("deadCellsVsRun");
	badCellsVsRun     =hFlagvsRun[1]->ProjectionY("badCellsVsRun");
	deadbadCellsVsRun =hFlagvsRun[2]->ProjectionY("badDeadCellsVsRun");
	cellSummaryCan->Divide(2);
	cellSummaryCan->cd(1);
	SetHisto(badCellsVsRun,"Run","No. of cells",0);
	badCellsVsRun->GetYaxis()->SetTitleOffset(1.7);
	badCellsVsRun->GetYaxis()->SetRangeUser(0,1300);
	badCellsVsRun->SetLineColor(kCyan+2);
	badCellsVsRun->SetLineWidth(2);
	badCellsVsRun->DrawCopy("hist");
	deadCellsVsRun->SetLineColor(kMagenta-2);
	deadCellsVsRun->SetLineWidth(2);
  	deadCellsVsRun->DrawCopy("same");

  	cellSummaryCan->cd(2);
  	SetHisto(nEventsVsRuns,"Run","No. of Events",0);
  	nEventsVsRuns->DrawCopy("hist");

  	cellSummaryCan2->Divide(2);
  	for(Int_t iRun=0;iRun<nRunsUsed;iRun++)
  	{
  		hBadVsEvent    ->Fill(nEventsVsRuns->GetBinContent(iRun+1),badCellsVsRun->GetBinContent(iRun+1));
  		hDeadBadVsEvent->Fill(nEventsVsRuns->GetBinContent(iRun+1),deadbadCellsVsRun->GetBinContent(iRun+1),1);
  	}
  	cellSummaryCan2->cd(1);
  	SetHisto(hBadVsEvent,"events in run","bad cells in run",0);
	hBadVsEvent->DrawCopy("colz");
 	cellSummaryCan2->cd(2);
	SetHisto(hDeadBadVsEvent,"events in run","bad+dead cells in run",0);
	hDeadBadVsEvent->DrawCopy("colz");
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//Draw bad & dead cells vs. ID
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	Color_t histCol=kCyan-8;
	//..Draw summary for all runs
	cFlagDeadBadI->cd()->SetLeftMargin(0.05);
	cFlagDeadBadI->cd()->SetRightMargin(0.05);
	cFlagDeadBadI->cd()->SetBottomMargin(0.06);
	cFlagDeadBadI->cd()->SetTopMargin(0.02);
	SetHisto(hFlagvsRun[0],"dead cell ID","Run No.",1); //ELI TT
	hFlagvsRun[0]->SetFillColor(histCol);
	hFlagvsRun[0]->Draw("BOX");
	cFlagDeadBadII->cd()->SetLeftMargin(0.05);
	cFlagDeadBadII->cd()->SetRightMargin(0.05);
	cFlagDeadBadII->cd()->SetBottomMargin(0.05);
	cFlagDeadBadII->cd()->SetTopMargin(0.02);
	SetHisto(hFlagvsRun[1],"bad cell ID","Run No.",1);
	hFlagvsRun[1]->SetFillColor(histCol);
	hFlagvsRun[1]->Draw("BOX");

	cFlagSumI->cd()->SetLeftMargin(0.05);
	cFlagSumI->cd()->SetRightMargin(0.05);
	cFlagSumI->cd()->SetBottomMargin(0.05);
	cFlagSumI->cd()->SetTopMargin(0.02);
	SetHisto(hFlagvsRun[2],"dead+bad cell ID pt.1","Run No.",1);
	hFlagvsRun[2]->SetFillColor(histCol);
	hFlagvsRun[2]->GetXaxis()->SetRangeUser(0,8837);
	hFlagvsRun[2]->DrawCopy("BOX");
	cFlagSumII->cd()->SetLeftMargin(0.05);
	cFlagSumII->cd()->SetRightMargin(0.05);
	cFlagSumII->cd()->SetBottomMargin(0.05);
	cFlagSumII->cd()->SetTopMargin(0.02);
	SetHisto(hFlagvsRun[2],"dead+bad cell ID pt.2","Run No.",1);
	hFlagvsRun[2]->GetXaxis()->SetRangeUser(8838,17674-2); //ELI why does that not work?
	hFlagvsRun[2]->DrawCopy("BOX");

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..compress the histogram for visibility since
	//..90% of the space is filled by empty good cells
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TH2F* CompressedAll = CompressHistogram(hFlagvsRun[2],noOfCells , nBadCells,RunIdVec);
	cFlagSumCompAllI->cd()->SetLeftMargin(0.05);
	cFlagSumCompAllI->cd()->SetRightMargin(0.05);
	cFlagSumCompAllI->cd()->SetBottomMargin(0.05);
	cFlagSumCompAllI->cd()->SetTopMargin(0.02);
	SetHisto(CompressedAll,"certain dead+bad cell pt.1","Run No.",1);
	CompressedAll->GetXaxis()->SetRangeUser(0,nBadCells/2);
	CompressedAll->SetFillColor(histCol);
	CompressedAll->DrawCopy("BOX");
	cFlagSumCompAllII->cd()->SetLeftMargin(0.05);
	cFlagSumCompAllII->cd()->SetRightMargin(0.05);
	cFlagSumCompAllII->cd()->SetBottomMargin(0.05);
	cFlagSumCompAllII->cd()->SetTopMargin(0.02);
	SetHisto(CompressedAll,"certain dead+bad cell pt.2","Run No.",1);
	CompressedAll->GetXaxis()->SetRangeUser(nBadCells/2,nBadCells+2);
	CompressedAll->DrawCopy("BOX");

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..................................................................
	// Find cells that are bad in a low fraction of runs
	//..................................................................
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"    o Summary: "<<nBadCells<<" bad cells of "<<noOfCells<<" total cells and "<<nRunsUsed<<" runs"<<endl;
	cout<<"    o 1 bad run out of "<<nRunsUsed<<" is "<<1.0/nRunsUsed<<endl;
	Double_t percbad = 0.20;       //..if a cell is only bad at 20% of the runs, test it again
	Double_t setBadCompletley=0.8; //..If a cell is bad in >80% of the runs then set it bad completley
	cout<<"    o Cells with "<<percbad<<"% of bad runs are double checked. These are "<<nRunsUsed*percbad<<" runs"<<endl;
	cout<<"    o Cell id with low fraction of bad runs:"<<endl;

	std::vector<Int_t> cellVector; //Filled with cells that onl have a small fraction of bad runs
	Double_t fracRun = 0, sumRun = 0;
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		TH1D *htmpCell = hFlagvsRun[2]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		fracRun = 0, sumRun = 0;
		for(Int_t ir = 0; ir < nRunsUsed ; ir++)
		{
			sumRun += htmpCell->GetBinContent(ir+1);
		}
		fracRun = sumRun/(Double_t)(nRunsUsed);

		//..loose selection criteria to remove the runs with zero entries
		if(fracRun>0)
		{
			Int_t cellColumn=0, cellRow=0;
			Int_t cellColumnAbs=0, cellRowAbs=0;
			Int_t trash = 0 ;
			fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(ic,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
			hFlagNew->Fill(cellColumnAbs, cellRowAbs, fracRun*100);
			//..If a cell is bad in >80% of the runs then set it bad completley
			if(fracRun>setBadCompletley)
			{
				for(Int_t j = 0 ; j < nRunsUsed; j++)
				{
					hFlagvsRun[2]->SetBinContent(ic+1,j+1,1);
				}
			}
			//..If a cell is bad in a low fraction of runs double check it
			if(fracRun<percbad)
			{
				cout<<ic<<", "<<flush;
				cellVector.push_back(ic);
			}
		}
	}
	cout<<endl;
	cout<<"    o In total "<<cellVector.size()<<" cells fall under this category"<<endl;
	//..................................................................
	// Plot cells with low fraction if bad runs
	// Double checks if they are really bad and re-includes them
	//..................................................................
	TString pdfName  = Form("%s/RunByRunSummary/%s_LowFractionCells",analysisOutput.Data(),listName.Data());
	PlotLowFractionCells(pdfName,cellVector,hFlagvsRun,nRunsUsed,ampID,hCellGoodMean);

	//..................................................................
	// Draw masked cell amplitude by masking cells that were identified bad or dead in a certain runblock
	//..................................................................
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		TH1D *htmpCellAllRuns     =hFlagvsRun[2]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		Double_t integralBlock1   =htmpCellAllRuns->Integral(0,splitRuns1);
		Double_t integralBlock2   =htmpCellAllRuns->Integral(splitRuns1+1,splitRuns2);
		Double_t integralBlock3   =htmpCellAllRuns->Integral(splitRuns2+1,splitRuns3);
		Double_t integralBlock4   =htmpCellAllRuns->Integral(splitRuns3+1,nRunsUsed);

		//..manually mask cells
		if(badcellsBlock1.size()>0 && ic==badcellsBlock1.at(0))integralBlock1=1;
		if(badcellsBlock2.size()>0 && ic==badcellsBlock2.at(0))integralBlock2=1;
		if(badcellsBlock3.size()>0 && ic==badcellsBlock3.at(0))integralBlock3=1;
		if(badcellsBlock4.size()>0 && ic==badcellsBlock4.at(0))integralBlock4=1;

		Double_t integralBlock3Sum=integralBlock1+integralBlock2+integralBlock3;
		//..only if the cell is bad in 1 run we will start the
		//..masking proceedure
		if(integralBlock3Sum>0)
		{
			for(Int_t ir = 0; ir < nRunsUsed ; ir++)
			{
				if((integralBlock1>0 && ir<=splitRuns1) ||
				   (integralBlock2>0 && ir>splitRuns1 && ir<=splitRuns2) ||
				   (integralBlock3>0 && ir>splitRuns2 && ir<=splitRuns3) ||
				   (integralBlock4>0 && ir>splitRuns3))
				{
					for (Int_t amp = 1; amp <= ampIDCl3Block[ir]->GetNbinsX(); amp++)
					{
						//..mask the cell, if it is once bad in the specific block
						ampIDCl3Block[ir]->SetBinContent(amp,ic+1,0);
					}
				}
                //..mask the cell if it is bad in any run in the period
				for (Int_t amp = 1; amp <= ampIDCl1Block[ir]->GetNbinsX(); amp++)
				{
					ampIDCl1Block[ir]->SetBinContent(amp,ic+1,0);
				}
			}
		}
	}
	//..build projections of amplitudes
	for(Int_t ir = 0; ir < nRunsUsed ; ir++)
	{
		cAmp[ir/totalperCv]->cd(ir%totalperCv+1)->SetLogy();
		TH1D* proj      = ampIDCl3Block[ir]->ProjectionX(TString::Format("hampIDMaskedCleanedProj_Run%d",RunIdVec.at(ir)));
		TH1D* projBlock = ampIDCl1Block[ir]->ProjectionX(TString::Format("hampIDMaskedCleaned1blockProj_Run%d",RunIdVec.at(ir)));
		SetHisto(proj,"","",0);
		proj->SetLineColor(kBlue-1);
		proj->Draw("hist same");
		//..create a summ version
		if(ir==0)
		{
			projSumC3Blocks = ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlock");
			projSumC3BlocksA= ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlockA");
			projSumC3BlocksB= ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlockB");
			projSumC3BlocksC= ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlockC");
			projSumC3BlocksD= ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlockD");
			projSumC3BlocksB->Reset();
			projSumC3BlocksC->Reset();
			projSumC3BlocksD->Reset();
			projSumC1Block = ampIDCl1Block[ir]->ProjectionX("hampIDProj_SumMaskedOneBlock");
		}
		else
		{
			                                    projSumC3Blocks->Add(proj);
			if(ir<=splitRuns1)                  projSumC3BlocksA->Add(proj);
			if(ir>splitRuns1 && ir<=splitRuns2) projSumC3BlocksB->Add(proj);
			if(ir>splitRuns2 && ir<=splitRuns3) projSumC3BlocksC->Add(proj);
			if(ir>splitRuns3)                   projSumC3BlocksD->Add(proj);
			projSumC1Block ->Add(projBlock);
		}
		Sum2D3BlockMask->Add(ampIDCl3Block[ir]);
		if(ir<=splitRuns1)                 Sum2D3BlockMaskA->Add(ampIDCl3Block[ir]);
		if(ir>splitRuns1 && ir<=splitRuns2)Sum2D3BlockMaskB->Add(ampIDCl3Block[ir]);
		if(ir>splitRuns2 && ir<=splitRuns3)Sum2D3BlockMaskC->Add(ampIDCl3Block[ir]);
		if(ir>splitRuns3)                  Sum2D3BlockMaskD->Add(ampIDCl3Block[ir]);
	}
	//..................................................................
	// Draw the cells masked in each respective run and summed amplitudes
	//..................................................................
	cAmpSum->Divide(2);
	cAmpSum->cd(1)->SetLogy();
	SetHisto(projSum,"","No. of hits/event",0);
	projSum->GetYaxis()->SetRangeUser(0.0000001,1000);
	projSum->SetLineColor(kCyan+3);
	projSum->SetLineWidth(2);
	projSum->Draw("hist");

	projSumC->SetLineColor(kSpring-2);
	projSumC->SetLineWidth(2);
	projSumC->DrawCopy("hist same");

	projSumC3Blocks->SetLineColor(kCyan-8);
	projSumC3Blocks->SetLineWidth(2);
	projSumC3Blocks->DrawCopy("hist same");
	projSumC1Block->SetLineColor(2);
	projSumC1Block->DrawCopy("hist same");
	projSumC3BlocksA->SetLineColor(4);
	projSumC3BlocksA->DrawCopy("hist same");
	projSumC3BlocksB->SetLineColor(5);
	projSumC3BlocksB->DrawCopy("hist same");
	projSumC3BlocksC->SetLineColor(6);
	projSumC3BlocksC->DrawCopy("hist same");
	projSumC3BlocksD->SetLineColor(kRed-7);
	projSumC3BlocksD->DrawCopy("hist same");
	TLegend *legSum = new TLegend(0.35,0.70,0.55,0.85);
	legSum->AddEntry(projSum,"Original engery distr.","l");
	legSum->AddEntry(projSumC,"Cells masked in each run","l");
	legSum->AddEntry(projSumC3Blocks,"Cells reincluded and masked in 3 blocks","l");
	legSum->SetBorderSize(0);
	legSum->SetTextSize(0.03);
	legSum->Draw("same");

	cAmpSum->cd(2);
	projSumC3Blocks->Divide(projSumC);
	SetHisto(projSumC3Blocks,"","block masked/single masked",0);
	projSumC3Blocks->SetLineColor(30);
	projSumC3Blocks->DrawCopy("hist");
	projSumC1Block->SetLineColor(4);
	projSumC1Block->Divide(projSumC);
	projSumC1Block->DrawCopy("same");

	projSumC1Block->SetLineStyle(3);
	projSumC1Block->SetLineColor(6);
	projSumC1Block->Divide(projSumC3Blocks);
	projSumC1Block->DrawCopy("same");

	TLegend *legRatio = new TLegend(0.25,0.30,0.4,0.45);
	legRatio->AddEntry(projSumC1Block,"re-incl & masked as one big block","l");
	legRatio->AddEntry(projSumC3Blocks,"re-incl & masked in three blocks","l");
	legRatio->AddEntry(projSumC1Block,"ratio 1Block/3Blocks","l");
	legRatio->SetBorderSize(0);
	legRatio->SetTextSize(0.03);
	legRatio->Draw("same");

	//..................................................................
	// Draw the Amp vs E histogram and the ratio to the masked 2D one
	//..................................................................
	cAmpSum2D->Divide(2);
	cAmpSum2D->cd(1)->SetLogz();
	SetHisto(Sum2D3BlockMask,"","",0);
	Sum2D3BlockMask->GetZaxis()->SetRangeUser(10e-8,10);
	Sum2D3BlockMask->DrawCopy("colz");

	cAmpSum2D->cd(2);
	Sum2D3BlockMask->Divide(Sum2DIdeal);
	SetHisto(Sum2D3BlockMask,"","",0);
	Sum2D3BlockMask->SetTitle("Amplitude Vs. run / av. amp. vs. run");
	//Sum2D3BlockMask->GetZaxis()->SetRangeUser(2,20);
	Sum2D3BlockMask->DrawCopy("colz");

	//..................................................................
	TLatex* textA = new TLatex(0.5,0.8,Form("Block over run %i-%i",0,splitRuns1));
	textA->SetTextSize(0.04);
	textA->SetTextColor(1);
	textA->SetNDC();
	TLatex* textB = new TLatex(0.5,0.8,Form("Block over run %i-%i",splitRuns1,splitRuns2));
	textB->SetTextSize(0.04);
	textB->SetTextColor(1);
	textB->SetNDC();
	TLatex* textC = new TLatex(0.5,0.8,Form("Block over run %i-%i",splitRuns2,splitRuns3));
	textC->SetTextSize(0.04);
	textC->SetTextColor(1);
	textC->SetNDC();
	TLatex* textD = new TLatex(0.5,0.8,Form("Block over run %i-%i",splitRuns3,nRunsUsed));
	textD->SetTextSize(0.04);
	textD->SetTextColor(1);
	textD->SetNDC();

	cAmpSum2D4Blocks->Divide(2,2);
	cAmpSum2D4Blocks->cd(1)->SetLogz();
	SetHisto(Sum2D3BlockMaskA,"","",0);
	Sum2D3BlockMaskA->GetZaxis()->SetRangeUser(10e-8,10);
	Sum2D3BlockMaskA->DrawCopy("colz");
	textA->Draw();
	cAmpSum2D4Blocks->cd(2)->SetLogz();
	SetHisto(Sum2D3BlockMaskB,"","",0);
	Sum2D3BlockMaskB->GetZaxis()->SetRangeUser(10e-8,10);
	Sum2D3BlockMaskB->DrawCopy("colz");
	textB->Draw();
	cAmpSum2D4Blocks->cd(3)->SetLogz();
	SetHisto(Sum2D3BlockMaskC,"","",0);
	Sum2D3BlockMaskC->GetZaxis()->SetRangeUser(10e-8,10);
	Sum2D3BlockMaskC->DrawCopy("colz");
	textC->Draw();
	cAmpSum2D4Blocks->cd(4)->SetLogz();
	SetHisto(Sum2D3BlockMaskD,"","",0);
	Sum2D3BlockMaskD->GetZaxis()->SetRangeUser(10e-8,10);
	Sum2D3BlockMaskD->DrawCopy("colz");
	textD->Draw();

	cAmpSum2D4BlocksRatio->Divide(2,2);
	cAmpSum2D4BlocksRatio->cd(1)->SetLogz();
	Sum2D3BlockMaskA->Divide(Sum2DIdeal);
	SetHisto(Sum2D3BlockMaskA,"","",0);
	Sum2D3BlockMaskA->SetTitle("Amplitude Vs. run / av. amp. vs. run");
	Sum2D3BlockMaskA->DrawCopy("colz");
	textA->Draw();
	cAmpSum2D4BlocksRatio->cd(2)->SetLogz();
	Sum2D3BlockMaskB->Divide(Sum2DIdeal);
	SetHisto(Sum2D3BlockMaskB,"","",0);
	Sum2D3BlockMaskB->SetTitle("Amplitude Vs. run / av. amp. vs. run");
	Sum2D3BlockMaskB->DrawCopy("colz");
	textB->Draw();
	cAmpSum2D4BlocksRatio->cd(3)->SetLogz();
	Sum2D3BlockMaskC->Divide(Sum2DIdeal);
	SetHisto(Sum2D3BlockMaskC,"","",0);
	Sum2D3BlockMaskC->SetTitle("Amplitude Vs. run / av. amp. vs. run");
	Sum2D3BlockMaskC->DrawCopy("colz");
	textC->Draw();
	cAmpSum2D4BlocksRatio->cd(4)->SetLogz();
	Sum2D3BlockMaskD->Divide(Sum2DIdeal);
	SetHisto(Sum2D3BlockMaskD,"","",0);
	Sum2D3BlockMaskD->SetTitle("Amplitude Vs. run / av. amp. vs. run");
	Sum2D3BlockMaskD->DrawCopy("colz");
	textD->Draw();

	//..................................................................
	//..Plot the cleaned up histogram again
	//..................................................................
	cFlagSumCleanedI->cd()->SetLeftMargin(0.05);
	cFlagSumCleanedI->cd()->SetRightMargin(0.05);
	cFlagSumCleanedI->cd()->SetBottomMargin(0.05);
	cFlagSumCleanedI->cd()->SetTopMargin(0.02);
	SetHisto(hFlagvsRun[2],"dead+bad cell ID (cleaned) pt.1","Run No.",1);
	hFlagvsRun[2]->GetXaxis()->SetRangeUser(0,8837);
	hFlagvsRun[2]->SetFillColor(histCol);
	hFlagvsRun[2]->DrawCopy("BOX");
	cFlagSumCleanedII->cd()->SetLeftMargin(0.05);
	cFlagSumCleanedII->cd()->SetRightMargin(0.05);
	cFlagSumCleanedII->cd()->SetBottomMargin(0.05);
	cFlagSumCleanedII->cd()->SetTopMargin(0.02);
	SetHisto(hFlagvsRun[2],"dead+bad cell ID (cleaned) pt.2","Run No.",1);
	hFlagvsRun[2]->GetXaxis()->SetRangeUser(8837,17674);
	hFlagvsRun[2]->DrawCopy("BOX");

	//..................................................................
	// Draw summary histogram with dead and bad cells vs run
	//..................................................................
	badCellsVsRunC  =hFlagvsRun[1]->ProjectionY("badCellsVsRunC");
	cellSummaryCan->cd(1);
	badCellsVsRunC->SetLineColor(kGreen-9);
	badCellsVsRunC->SetLineWidth(2);
	badCellsVsRunC->SetFillColorAlpha(10,0);
	badCellsVsRunC->DrawCopy("same");

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..compress the histogram for visibility since
	//..90% of the space is filled by empty good cells
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<endl;
	TH2F* CompressedClean = CompressHistogram(hFlagvsRun[2],noOfCells , nBadCells,RunIdVec);

	cFlagSumCompCleanI->cd()->SetLeftMargin(0.05);
	cFlagSumCompCleanI->cd()->SetRightMargin(0.05);
	cFlagSumCompCleanI->cd()->SetBottomMargin(0.05);
	cFlagSumCompCleanI->cd()->SetTopMargin(0.02);
	SetHisto(CompressedClean,"certain dead+bad cell (cleaned) pt.1","Run No.",1);
	CompressedClean->GetXaxis()->SetRangeUser(0,nBadCells/2);
	CompressedClean->SetFillColor(histCol);
	CompressedClean->DrawCopy("BOX");
	cFlagSumCompCleanII->cd()->SetLeftMargin(0.05);
	cFlagSumCompCleanII->cd()->SetRightMargin(0.05);
	cFlagSumCompCleanII->cd()->SetBottomMargin(0.05);
	cFlagSumCompCleanII->cd()->SetTopMargin(0.02);
	SetHisto(CompressedClean,"certain dead+bad cell (cleaned) pt.2","Run No.",1);
	CompressedClean->GetXaxis()->SetRangeUser(nBadCells/2,nBadCells+2);
	CompressedClean->DrawCopy("BOX");

	//..................................................................
	//..build and draw a new 2D histogram after all the cleaning/resetting was done
	//..................................................................
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		TH1D *htmpCell = hFlagvsRun[2]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		fracRun = 0, sumRun = 0;
		for(Int_t ir = 0; ir < nRunsUsed ; ir++)
		{
			sumRun += htmpCell->GetBinContent(ir+1);
		}
		fracRun = sumRun/(Double_t)(nRunsUsed);
		if(fracRun>0)
		{
			Int_t cellColumn=0, cellRow=0;
			Int_t cellColumnAbs=0, cellRowAbs=0;
			Int_t trash = 0 ;
			fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(ic,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
			hFlagNewClean->Fill(cellColumnAbs, cellRowAbs, fracRun*100);
		}
	}
    //..2D fraction
    cFlagNew->Divide(2);
    cFlagNew->cd(1);
	SetHisto(hFlagNew,"","",0);
	hFlagNew->Draw("colz");
    cFlagNew->cd(2);
	SetHisto(hFlagNewClean,"","",0);
	hFlagNewClean->Draw("colz");

	//..................................................................
	// Print out an overview of this run-by-run analysis
	//..................................................................
	Int_t bcBolckSum=0;

	TH1D *htmpCellSum = hFlagvsRun[2]->ProjectionX("countSum");
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		if(htmpCellSum->GetBinContent(ic+1)>0)bcBolckSum++;
	}

    cout<<"...................................."<<endl;
	cout<<"Summary: "<<nBadCells<<" bad cells of "<<noOfCells<<" total cells and "<<nRunsUsed<<" runs"<<endl;
    cout<<"o All bad cells: "<<nBadCells<<endl;
    cout<<"o low fraction bad cells: "<<cellVector.size()<<endl;
    cout<<"o bad cells after reinclusion: "<<bcBolckSum<<endl;
    cout<<"o rescued cells: "<<nBadCells-bcBolckSum<<endl;
    if(cellVector.size()!=0)cout<<"o These are: "<<100*(nBadCells-bcBolckSum)/cellVector.size()<<"% of low frac cells"<<endl;
    cout<<"...................................."<<endl;

	//_______________________________________________________________________
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . Save histograms to file
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//_______________________________________________________________________
	TString fileName       = Form("%s/RunByRunSummary%i/%s_Results.root",analysisOutput.Data(),runsUsed,listName.Data());
	TFile* rootFile        = new TFile(fileName,"recreate");

	for(Int_t ic = 0; ic<nCv; ic++)
	{
		rootFile->WriteObject(cBad [ic],cBad [ic]->GetName());
		rootFile->WriteObject(cGood[ic],cGood[ic]->GetName());
		rootFile->WriteObject(cDead[ic],cDead[ic]->GetName());
		rootFile->WriteObject(cAmp [ic],cAmp [ic]->GetName());
	}
	rootFile->WriteObject(cellSummaryCan,cellSummaryCan->GetName());
	rootFile->WriteObject(cAmpSum,cAmpSum->GetName());
	rootFile->WriteObject(cFlagDeadBadI,cFlagDeadBadI->GetName());
	rootFile->WriteObject(cFlagDeadBadII,cFlagDeadBadII->GetName());
	rootFile->WriteObject(cFlagSumI,cFlagSumII->GetName());
	rootFile->WriteObject(cFlagSumII,cFlagSumII->GetName());
	rootFile->WriteObject(cFlagSumCleanedI,cFlagSumCleanedI->GetName());
	rootFile->WriteObject(cFlagSumCleanedII,cFlagSumCleanedII->GetName());
	rootFile->WriteObject(cFlagSumCompAllI,cFlagSumCompAllI->GetName());
	rootFile->WriteObject(cFlagSumCompAllII,cFlagSumCompAllII->GetName());
	rootFile->WriteObject(cFlagSumCompCleanI,cFlagSumCompCleanI->GetName());
	rootFile->WriteObject(cFlagSumCompCleanII,cFlagSumCompCleanII->GetName());
	rootFile->WriteObject(cFlagNew,cFlagNew->GetName());
	rootFile->WriteObject(cAmpSum2D,cAmpSum2D->GetName());
	rootFile->WriteObject(cAmpSum2D4Blocks,cAmpSum2D4Blocks->GetName());
	rootFile->WriteObject(cAmpSum2D4BlocksRatio,cAmpSum2D4BlocksRatio->GetName());
	//..histograms
	rootFile->WriteObject(hFlagNew,hFlagNew->GetName());
	rootFile->WriteObject(CompressedAll,CompressedAll->GetName());
	rootFile->WriteObject(hFlagvsRun[2],hFlagvsRun[2]->GetName());      //..bad/dead vs. run number cleaned for low frac runs
	rootFile->WriteObject(CompressedClean,CompressedClean->GetName());  //..bad/dead vs. run number cleaned&compressed
	rootFile->WriteObject(projSum,projSum->GetName());
	rootFile->WriteObject(projSumC,projSumC->GetName());
	rootFile->WriteObject(projSumC3Blocks,projSumC3Blocks->GetName());
	rootFile->WriteObject(nEventsVsRuns,nEventsVsRuns->GetName());
	for(Int_t ir = 0; ir < nRunsUsed ; ir++)
	{
		rootFile->WriteObject(ampIDCl[ir],ampIDCl[ir]->GetName()); //..single run masked (original bad ch. analysis)
		rootFile->WriteObject(ampID[ir],ampID[ir]->GetName()); //..original distribution - unmasked
		rootFile->WriteObject(ampIDDelete[ir],ampIDDelete[ir]->GetName());
	}


	//..plot the canvases of cells into a .pdf file
	TString badCanvasName =Form("%s/RunByRunSummary%i/%s_BadCells.pdf", analysisOutput.Data(),nRunsUsed,listName.Data());
	TString goodCanvasName=Form("%s/RunByRunSummary%i/%s_GoodCells.pdf", analysisOutput.Data(),nRunsUsed,listName.Data());
	TString deadCanvasName=Form("%s/RunByRunSummary%i/%s_DeadCells.pdf", analysisOutput.Data(),nRunsUsed,listName.Data());
	TString ampCanvasName =Form("%s/RunByRunSummary%i/%s_Amplitudes.pdf", analysisOutput.Data(),nRunsUsed,listName.Data());
	TString summaryPDF    =Form("%s/RunByRunSummary%i/%s_FigureCollection.pdf", analysisOutput.Data(),nRunsUsed,listName.Data());

	for(Int_t ic = 0; ic<nCv; ic++)
	{
		if(ic==0 && nCv>1)
		{
			//..first pad
			cBad [ic] ->Print(Form("%s(",badCanvasName.Data()));
			cGood[ic] ->Print(Form("%s(",goodCanvasName.Data()));
			cDead[ic] ->Print(Form("%s(",deadCanvasName.Data()));
			cAmp [ic] ->Print(Form("%s(",ampCanvasName.Data()));
		}
		else if(ic==(nCv-1) && nCv>1)
		{
			//..last pad
			cBad [ic] ->Print(Form("%s)",badCanvasName.Data()));
			cGood[ic] ->Print(Form("%s)",goodCanvasName.Data()));
			cDead[ic] ->Print(Form("%s)",deadCanvasName.Data()));
			cAmp [ic] ->Print(Form("%s)",ampCanvasName.Data()));
		}
		else
		{
			//..all pads in between
			cBad [ic] ->Print(Form("%s",badCanvasName.Data()));
			cGood[ic] ->Print(Form("%s",goodCanvasName.Data()));
			cDead[ic] ->Print(Form("%s",deadCanvasName.Data()));
			cAmp [ic] ->Print(Form("%s",ampCanvasName.Data()));
		}
	}

	//..Add figures to the summary canvas
	cellSummaryCan   ->Print(Form("%s(",summaryPDF.Data()));
	cellSummaryCan2  ->Print(Form("%s",summaryPDF.Data()));
	cAmpSum          ->Print(Form("%s",summaryPDF.Data()));
	cAmpSum2D        ->Print(Form("%s",summaryPDF.Data()));
	cAmpSum2D4Blocks ->Print(Form("%s",summaryPDF.Data()));
	cFlagNew         ->Print(Form("%s",summaryPDF.Data()));
	cFlagDeadBadI    ->Print(Form("%s",summaryPDF.Data()));
	cFlagDeadBadII   ->Print(Form("%s",summaryPDF.Data()));
	cFlagSumI        ->Print(Form("%s",summaryPDF.Data()));
	cFlagSumII       ->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCleanedI ->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCleanedII->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCompAllI ->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCompAllII->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCompCleanI->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCompCleanII->Print(Form("%s)",summaryPDF.Data()));

}
///
/// Calculate how to best split periods into run blocks
/// So that the number of masked cells is minimized
//________________________________________________________________________
void GetBestPeriodSplitting(TString period = "LHC15o", Int_t train = 771,Int_t Nruns=105, Int_t noOfSplits=4)
{
	if(noOfSplits>4)cout<<"Error: so far only implemented for 1-4 splits"<<endl;

	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	//gStyle->SetPalette(53);  //standard is 1
	gStyle->SetCanvasColor(10);
	TGaxis::SetMaxDigits(4);
	gStyle->SetPadTopMargin(0.07);//0.05
	gStyle->SetPadBottomMargin(0.18);//0.15
	gStyle->SetPadRightMargin(0.04);
	gStyle->SetPadLeftMargin(0.06);
	gStyle->SetFrameFillColor(10);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetTitleSize(5.0,"X");
	gStyle->SetTitleSize(5.0,"Y");
	gEnv->SetValue("Canvas.ShowEventStatus",1);  //shows the status bar in the canvas

	//.......................................
	//..Settings
	//.......................................
	//TString runList = "GloballyGood";
	TString runList = "GloballyGood10";
	Int_t runNumber=245145;
	Bool_t weightByNevt=1;  //..weight the run period splitting by the number of events in the run


	//.......................................
	//..Get the Number of cells
	//.......................................
	AliCalorimeterUtils *fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(runNumber);
	fCaloUtils->AccessGeometry(aod);
    Int_t noOfCells=fCaloUtils->GetEMCALGeometry()->GetNCells(); //..Very important number, never change after that point!

	//.......................................
	//..Get run-by-run analysis output
	//.......................................
	TString fileName = Form("AnalysisOutput/%s/Train_%i/RunByRunSummary%i/%s_Results.root",period.Data(),train,Nruns,runList.Data());

	cout<<"Open root file: "<<fileName<<endl;
	TFile *f = TFile::Open(fileName);
	if(!f)
	{
		cout<<"Couldn't open/find .root file: "<<fileName<<endl;
	}
	//..Get the histograms
	TH2F *hbadAndDeadvsRun = (TH2F*)f->Get("hFlag3vsRun_Comp");
	TH2F *hbadAndDeadvsRunw = (TH2F*)hbadAndDeadvsRun->Clone("hbadAndDeadvsRunWeight");
	if(!hbadAndDeadvsRun)cout<<"couldn't find histogram - hFlag3vsRun_Comp"<<endl;
	hbadAndDeadvsRun->GetXaxis()->UnZoom();
	hbadAndDeadvsRunw->GetXaxis()->UnZoom();

	TH1D* hEventsPerRun = (TH1D*)f->Get("nEventVsRun");
	if(!hEventsPerRun)cout<<"couldn't find histogram - nEventVsRun"<<endl;

	//TH1D *hEventsPerRunHI = (TH1D)hEventsPerRun->Clone("hEventsPerRunHI");
    // 2DPlot
	// flag vs run

	//.......................................
	//..Prepare the histogram for the splitting analysis
	//.......................................
	if(weightByNevt==1)
	{
		//..weight bad cells by numbers of events
		for(Int_t icell = 0; icell < noOfCells ; icell++)
		{
			for(Int_t iRun = 0; iRun < Nruns ; iRun++)
			{
				if(hbadAndDeadvsRun->GetBinContent(icell+1,iRun+1)==1)
				{
					hbadAndDeadvsRunw->SetBinContent(icell+1,iRun+1,hEventsPerRun->GetBinContent(iRun+1));
				}
			}
		}
	}
	//..Draw the histogram
	TCanvas *canEvt= new TCanvas("canEvt", "canEvt", 500, 500);
	canEvt->cd()->SetLeftMargin(0.08);
	SetHisto(hEventsPerRun,"Run","No of events",0);
	//hEventsPerRun->GetYaxis()->SetTitleOffset(0.35);
	hEventsPerRun->DrawCopy("hist"); //box

	TCanvas *can1= new TCanvas("compressedIDs1", "compressed cell ID's A)", 1600, 1000);
	can1->cd()->SetLeftMargin(0.05);
	can1->cd()->SetRightMargin(0.05);
	can1->cd()->SetBottomMargin(0.06);
	can1->cd()->SetTopMargin(0.02);
	SetHisto(hbadAndDeadvsRun,"","Run in List",1);
	hbadAndDeadvsRun->GetYaxis()->SetTitleOffset(0.35);
	hbadAndDeadvsRun->GetXaxis()->SetRangeUser(0,2300);
	hbadAndDeadvsRun->DrawCopy("colz"); //box
	TCanvas *can2= new TCanvas("compressedIDs2", "compressed cell ID's B)", 1600, 1000);
	can2->cd()->SetLeftMargin(0.05);
	can2->cd()->SetRightMargin(0.05);
	can2->cd()->SetBottomMargin(0.06);
	can2->cd()->SetTopMargin(0.02);
	SetHisto(hbadAndDeadvsRunw,"","Run in List",1);
	hbadAndDeadvsRunw->GetYaxis()->SetTitleOffset(0.35);
	hbadAndDeadvsRunw->GetXaxis()->SetRangeUser(0,2300);
	hbadAndDeadvsRunw->DrawCopy("colz");

	//.......................................
	//..Start analysis for getting the best possible splitting
	//.......................................
	//..Analyze the histogram
	//..split into three blocks and see what performs better
	Int_t splitRun1=0;
	Int_t splitRun2=0;
	Int_t splitRun3=0;
	Int_t splitRun1w=0;
	Int_t splitRun2w=0;
	Int_t splitRun3w=0;
	Int_t endBlock2;
	Int_t endBlock3;
	Double_t totalCellsBadRun=0;
	Double_t totalCellsBadEvt=0;
	Double_t nCellRunsBlock1=0;
	Double_t nCellRunsBlock2=0;
	Double_t nCellRunsBlock3=0;
	Double_t nCellRunsBlock4=0;
	Double_t nCellEvtBlock1=0;
	Double_t nCellEvtBlock2=0;
	Double_t nCellEvtBlock3=0;
	Double_t nCellEvtBlock4=0;

	Int_t sumRunBlockTotal=0;
	Int_t nBlockTotal=0;
	TH1D *hOneBigBlock = hbadAndDeadvsRunw->ProjectionX(TString::Format("%s_projSum",hbadAndDeadvsRunw->GetName()),0,Nruns);

	for(Int_t iRun=1; iRun<=Nruns; iRun++)
	{
		cout<<"Round "<<iRun<<" of "<<Nruns<<endl;

		if(noOfSplits<3)
		{
			endBlock2=Nruns-1;
		}
		else
		{
			endBlock2=iRun+1;
		}
		for(Int_t iRun2=endBlock2; iRun2<=Nruns; iRun2++)
		{
			if(noOfSplits<4)
			{
				endBlock3=Nruns-1;
			}
			else
			{
				endBlock3=iRun2+1;
			}

			for(Int_t iRun3=endBlock3; iRun3<=Nruns; iRun3++)
			{
				hbadAndDeadvsRun ->GetXaxis()->UnZoom();
				hbadAndDeadvsRunw->GetXaxis()->UnZoom();
				TH1D *htmpCell1;
				TH1D *htmpCell2;
				TH1D *htmpCell3;
				TH1D *htmpCell4;

				if(noOfSplits<3)iRun2=Nruns;
				if(noOfSplits<4)iRun3=Nruns;

				if(weightByNevt==1)
				{
					htmpCell1 = hbadAndDeadvsRunw->ProjectionX(TString::Format("%s_proj1",hbadAndDeadvsRunw->GetName()),0,iRun);
					htmpCell2 = hbadAndDeadvsRunw->ProjectionX(TString::Format("%s_proj2",hbadAndDeadvsRunw->GetName()),iRun+1,iRun2);
					if(noOfSplits>2)htmpCell3 = hbadAndDeadvsRunw->ProjectionX(TString::Format("%s_proj3",hbadAndDeadvsRunw->GetName()),iRun2+1,iRun3);
					if(noOfSplits>3)htmpCell4 = hbadAndDeadvsRunw->ProjectionX(TString::Format("%s_proj4",hbadAndDeadvsRunw->GetName()),iRun3+1,Nruns);
				}
				else
				{
					htmpCell1 = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj1",hbadAndDeadvsRun->GetName()),0,iRun);
					htmpCell2 = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj2",hbadAndDeadvsRun->GetName()),iRun+1,iRun2);
					if(noOfSplits>2)htmpCell3 = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj3",hbadAndDeadvsRun->GetName()),iRun2+1,iRun3);
					if(noOfSplits>3)htmpCell4 = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj4",hbadAndDeadvsRun->GetName()),iRun3+1,Nruns);
				}

				Double_t nEvtBlock1 = hEventsPerRun->Integral(0,iRun);
				Double_t nEvtBlock2 = hEventsPerRun->Integral(iRun+1,iRun2);
				Double_t nEvtBlock3 = hEventsPerRun->Integral(iRun2+1,iRun3);
				Double_t nEvtBlock4 = hEventsPerRun->Integral(iRun3+1,Nruns);

				sumRunBlockTotal=0;
				Int_t sumRunBlock1=0;
				Int_t sumRunBlock2=0;
				Int_t sumRunBlock3=0;
				Int_t sumRunBlock4=0;
				Int_t nBlock1=0;
				Int_t nBlock2=0;
				Int_t nBlock3=0;
				Int_t nBlock4=0;
				nBlockTotal=0;

				for(Int_t icell = 0; icell < htmpCell1->GetNbinsX(); icell++)
				{
					sumRunBlockTotal             = hOneBigBlock->GetBinContent(icell+1);
					sumRunBlock1                 = htmpCell1->GetBinContent(icell+1);
					sumRunBlock2                 = htmpCell2->GetBinContent(icell+1);
					if(noOfSplits>2)sumRunBlock3 = htmpCell3->GetBinContent(icell+1);
					if(noOfSplits>3)sumRunBlock4 = htmpCell4->GetBinContent(icell+1);

					if(sumRunBlock1>0)nBlock1++;
					if(sumRunBlock2>0)nBlock2++;
					if(sumRunBlock3>0)nBlock3++;
					if(sumRunBlock4>0)nBlock4++;
					if(sumRunBlockTotal>0)nBlockTotal++;
				}
				//..bad cells in block
				nCellRunsBlock1=nBlock1*iRun;
				nCellRunsBlock2=nBlock2*(iRun2-iRun+1);
				nCellRunsBlock3=nBlock3*(iRun3-iRun2+1);  //..is 0 for 2block splitting
				nCellRunsBlock4=nBlock4*(Nruns-iRun3+1);  //..is 0 for 3block splitting

				//..bad cells in block weightet by nEvents in block
				nCellEvtBlock1 =nBlock1*nEvtBlock1;
				nCellEvtBlock2 =nBlock2*nEvtBlock2;
				nCellEvtBlock3 =nBlock3*nEvtBlock3;
				nCellEvtBlock4 =nBlock4*nEvtBlock4;

				//..not weighted by nuber of events in run
				if(totalCellsBadRun==0 || (nCellRunsBlock1+nCellRunsBlock2+nCellRunsBlock3+nCellRunsBlock4)<totalCellsBadRun)
				{
					totalCellsBadRun=nCellRunsBlock1+nCellRunsBlock2+nCellRunsBlock3+nCellRunsBlock4;
					splitRun1=iRun;
					splitRun2=iRun2;
					splitRun3=iRun3;
				}
				//..weighted by nuber of events in run
				if(totalCellsBadEvt==0 || (nCellEvtBlock1+nCellEvtBlock2+nCellEvtBlock3+nCellEvtBlock4)<totalCellsBadEvt)
				{
					totalCellsBadEvt=nCellEvtBlock1+nCellEvtBlock2+nCellEvtBlock3+nCellEvtBlock4;
					splitRun1w=iRun;
					splitRun2w=iRun2;
					splitRun3w=iRun3;
				}
			}
		}
	}
	hbadAndDeadvsRun->GetXaxis()->UnZoom();
	TH1D *htmpCell1p = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj1",hbadAndDeadvsRun->GetName()),0,splitRun1);
	TH1D *htmpCell2p = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj2",hbadAndDeadvsRun->GetName()),splitRun1+1,splitRun2);
	TH1D *htmpCell3p = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj3",hbadAndDeadvsRun->GetName()),splitRun2+1,splitRun3);
	TH1D *htmpCell4p = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj4",hbadAndDeadvsRun->GetName()),splitRun3+1,Nruns);

	TCanvas *canSplit= new TCanvas("canSplit", "Split compressed cell ID's", 1600, 500);
	canSplit->Divide(2,2);
	canSplit->cd(1);
	SetHisto(htmpCell1p,"","nruns bad",1);
	htmpCell1p->GetYaxis()->SetTitleOffset(0.35);
	htmpCell1p->DrawCopy("hist");
	canSplit->cd(2);
	SetHisto(htmpCell2p,"","nruns bad",1);
	htmpCell2p->GetYaxis()->SetTitleOffset(0.35);
	htmpCell2p->DrawCopy("hist");
	canSplit->cd(3);
	SetHisto(htmpCell3p,"","nruns bad",1);
	htmpCell3p->GetYaxis()->SetTitleOffset(0.35);
	htmpCell3p->DrawCopy("hist");
	canSplit->cd(4);
	SetHisto(htmpCell4p,"","nruns bad",1);
	htmpCell4p->GetYaxis()->SetTitleOffset(0.35);
	htmpCell4p->DrawCopy("hist");

	//..Draw the split lines into the canvas
	can1->cd();
	PlotHorLineRange(splitRun1w,0,2000,9);
	PlotHorLineRange(splitRun2w,0,2000,9);
	PlotHorLineRange(splitRun3w,0,2000,9);
	can2->cd();
	PlotHorLineRange(splitRun1w,0,2000,1);
	PlotHorLineRange(splitRun2w,0,2000,1);
	PlotHorLineRange(splitRun3w,0,2000,1);


	cout<<"Best results are achieved by splitting into:"<<endl;
	cout<<"- - - - - - - - GENERAL - - - - - - - - -"<<endl;
	cout<<"Number of runs: "<<Nruns<<endl;
	cout<<"Number events: "<<hEventsPerRun->Integral(0,Nruns)<<endl;
	cout<<"Bad cells if masked as 1 block: "<<nBlockTotal<<endl;
	cout<<"Bad cells * runs: "<<nBlockTotal*Nruns<<endl;
	cout<<"Bad cells * evt: "<<nBlockTotal*hEventsPerRun->Integral(0,Nruns)<<endl;
	cout<<"- - - - - - - - UNWEIGHTED Splitting - - - - - - - - -"<<endl;
	cout<<"Run: 0-"<<splitRun1<<endl;
	cout<<"Run: "<<splitRun1+1<<"-"<<splitRun2<<endl;
	if(noOfSplits>2)cout<<"Run: "<<splitRun2+1<<"-"<<splitRun3<<endl;
	if(noOfSplits>3)cout<<"Run: "<<splitRun3+1<<"-"<<Nruns<<endl;
	cout<<"Number of bad cells*runs   ="<<totalCellsBadRun<<endl;
	cout<<"- - - - - - - - WEIGHTED Splitting - - - - - - - - -"<<endl;
	cout<<"Run: 0-"<<splitRun1w<<endl;
	cout<<"Run: "<<splitRun1w+1<<"-"<<splitRun2w<<endl;
	if(noOfSplits>2)cout<<"Run: "<<splitRun2w+1<<"-"<<splitRun3w<<endl;
	if(noOfSplits>3)cout<<"Run: "<<splitRun3w+1<<"-"<<Nruns<<endl;
	cout<<"Number of bad cells*events ="<<totalCellsBadEvt<<endl;

	cout<<" events block1 : "<<hEventsPerRun->Integral(0,splitRun1)<<endl;
	cout<<" events block2 : "<<hEventsPerRun->Integral(splitRun1+1,splitRun2)<<endl;
	cout<<" events block3 : "<<hEventsPerRun->Integral(splitRun2+1,splitRun3)<<endl;
	cout<<" events block4 : "<<hEventsPerRun->Integral(splitRun3+1,Nruns)<<endl;

	//..With this splitting produce bad map lists and show how the masked blocks look like
	//
}
//
// Compares masked amplidudes from two different versions of masking cells
// Version A sum all runs in a certrain range up and analyse the channels all together
// Version B do a run-by-run analsis and check wheather channels are masked in almost all runs and maks them then entirely
//________________________________________________________________________
void CompareTwoBCstrategies(TString period="LHC15n",Int_t trainNo=603,Int_t version=5)
{
	//..this was originally used or LHC15o
	gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //..to supress a lot of standard output
	gStyle->SetOptStat(0); //..Do not plot stat boxes
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadTopMargin(0.02);

    Int_t noOfCells=17674; //
    const Int_t nRunsUsed=105;
    //const Int_t nRunsUsed=10;
    //..select runs after which a new bad map is built
    Int_t splitRuns1=34;  //run bock is inclusive of this run
    Int_t splitRuns2=66;  //run bock is inclusive of this run
    Int_t splitRuns3=74;  //run bock is inclusive of this run
    Int_t nBadCellMerged[4]={noOfCells,noOfCells,noOfCells,noOfCells};
    Int_t nBadCellRbR[4]   ={noOfCells,noOfCells,noOfCells,noOfCells};

    //..............................................
    //..manually disable cells
    std::vector<Int_t> badcellsBlock1;
    std::vector<Int_t> badcellsBlock2;
    std::vector<Int_t> badcellsBlock3;
    std::vector<Int_t> badcellsBlock4;

    badcellsBlock1.insert(badcellsBlock1.end(),{14655,14622,14640,14728,14726});

    badcellsBlock2.insert(badcellsBlock2.end(),{6644,6655,10140,12036,12037,12038,12039,12040,12041,12926,13067,13066,13125});
    badcellsBlock2.insert(badcellsBlock2.end(),{13133,13483,13971,13978,14116,14118,14122,14411,14593,14599,14600,14606,14699});
    badcellsBlock2.insert(badcellsBlock2.end(),{14728,15158,15462,16309});

    badcellsBlock3.insert(badcellsBlock3.end(),{292,294,297,301,13483, 13975, 14116, 14320, 14326});

    badcellsBlock4.insert(badcellsBlock4.end(),{3472,3473,3474,3475,3476,3477,3478,3479,3480,3481,3482,3483,3484,3485,3486,3487});
    badcellsBlock4.insert(badcellsBlock4.end(),{3520,3521,3522,3523,3524,3525,3526,3527,3528,3529,3530,3531,3532,3533,3534,3535});
    badcellsBlock4.insert(badcellsBlock4.end(),{3665,3666,3667,3668,3669,3670,3671,3672,3673,3674,3675,3676,3677,3678,3679});
    badcellsBlock4.insert(badcellsBlock4.end(),{3712,3713,3714,3715,3716,3717,3718,3719,3720,3721,3722,3723,3724,3725,3726,3727});
    badcellsBlock4.insert(badcellsBlock4.end(),{8848,8849,8850,8851,8852,8853,8854,8855,8856,8857,8858,8859,8860,8861,8862,8863});
    badcellsBlock4.insert(badcellsBlock4.end(),{8896,8897,8898,8899,8900,8901,8902,8903,8904,8905,89106,8907,8908,8909,8910});
    badcellsBlock4.insert(badcellsBlock4.end(),{11906,11907,11908,11909,11910,11911,11912,11913,11914,11915,11916,11917,11918,11919});
    badcellsBlock4.insert(badcellsBlock4.end(),{12034,12035,12036,12037,12038,12039,12040,12041,12042,12043,12044,12045,13469,13483,16427,16430});

    std::vector<Int_t> vOnlyMaskedInMergedA;
	std::vector<Int_t> vOnlyMaskedInRunbRunA;
	std::vector<Int_t> vOnlyMaskedInMergedB;
	std::vector<Int_t> vOnlyMaskedInRunbRunB;
	std::vector<Int_t> vOnlyMaskedInMergedC;
	std::vector<Int_t> vOnlyMaskedInRunbRunC;
	std::vector<Int_t> vOnlyMaskedInMergedD;
	std::vector<Int_t> vOnlyMaskedInRunbRunD;

	cout<<"run splitting: "<<endl;
	cout<<0<<"-"<<splitRuns1<<" -> "<<splitRuns1<<" runs"<<endl;
	cout<<splitRuns1+1<<"-"<<splitRuns2<<" -> "<<splitRuns2-splitRuns1<<" runs"<<endl;
	cout<<splitRuns2+1<<"-"<<splitRuns3<<" -> "<<splitRuns3-splitRuns2<<" runs"<<endl;
	cout<<splitRuns3+1<<"-"<<nRunsUsed <<" -> "<<nRunsUsed-splitRuns3 <<" runs"<<endl;

	//......................................................
	//..Get the .root file from analyzing runs as 1 block together
	//......................................................
	cout<<"**Open files with merged runlist analysis: "<<endl;
	TString pathA        = Form("./AnalysisOutput/%s/Train_%i/Version%i",period.Data(),trainNo,version);
	TString rootFileNameA= Form("INT7_Histograms_V%i.root",version);
	TFile* outputRootA   = TFile::Open(Form("%s/%s",pathA.Data(),rootFileNameA.Data()));
	if(!outputRootA)cout<<"File "<<outputRootA->GetName()<<" does not exist"<<endl;
	else cout<<"file A: "<<outputRootA->GetName()<<endl;
	TString pathB        = Form("./AnalysisOutput/%s/Train_%i/Version%i",period.Data(),trainNo,version+1);
	TString rootFileNameB= Form("INT7_Histograms_V%i.root",version+1);
	TFile* outputRootB   = TFile::Open(Form("%s/%s",pathB.Data(),rootFileNameB.Data()));
	if(!outputRootB)cout<<"File "<<outputRootB->GetName()<<" does not exist"<<endl;
	else cout<<"file B: "<<outputRootB->GetName()<<endl;
	TString pathC        = Form("./AnalysisOutput/%s/Train_%i/Version%i",period.Data(),trainNo,version+2);
	TString rootFileNameC= Form("INT7_Histograms_V%i.root",version+2);
	TFile* outputRootC   = TFile::Open(Form("%s/%s",pathC.Data(),rootFileNameC.Data()));
	if(!outputRootC)cout<<"File "<<outputRootC->GetName()<<" does not exist"<<endl;
	else cout<<"file C: "<<outputRootC->GetName()<<endl;
	TString pathD        = Form("./AnalysisOutput/%s/Train_%i/Version%i",period.Data(),trainNo,version+3);
	TString rootFileNameD= Form("INT7_Histograms_V%i.root",version+3);
	TFile* outputRootD   = TFile::Open(Form("%s/%s",pathD.Data(),rootFileNameD.Data()));
	if(!outputRootD)cout<<"File "<<outputRootD->GetName()<<" does not exist"<<endl;
	else cout<<"file D: "<<outputRootD->GetName()<<endl;

	//..get the necessary histograms
	TH2F* hCellAmplitudeA    =(TH2F*)outputRootA->Get("hCellAmplitude");
	TH1F* hnEventsA          =(TH1F*)outputRootA->Get("hNEvents");
	hCellAmplitudeA->Scale(hnEventsA->GetBinContent(1));
	TH2F* hCellAmplitudeBlockMaskA=(TH2F*)hCellAmplitudeA->Clone("hCellAmplitudeMask");
	TH1F* hCellFlagA         =(TH1F*)outputRootA->Get("fhCellFlag");

	TH2F* hCellAmplitudeB    =(TH2F*)outputRootB->Get("hCellAmplitude");
	TH1F* hnEventsB          =(TH1F*)outputRootB->Get("hNEvents");
	hCellAmplitudeB->Scale(hnEventsB->GetBinContent(1));
	TH2F* hCellAmplitudeBlockMaskB=(TH2F*)hCellAmplitudeB->Clone("hCellAmplitudeMask");
	TH1F* hCellFlagB         =(TH1F*)outputRootB->Get("fhCellFlag");

	TH2F* hCellAmplitudeC    =(TH2F*)outputRootC->Get("hCellAmplitude");
	TH1F* hnEventsC          =(TH1F*)outputRootC->Get("hNEvents");
	hCellAmplitudeC->Scale(hnEventsC->GetBinContent(1));
	TH2F* hCellAmplitudeBlockMaskC=(TH2F*)hCellAmplitudeC->Clone("hCellAmplitudeMask");
	TH1F* hCellFlagC         =(TH1F*)outputRootC->Get("fhCellFlag");

	TH2F* hCellAmplitudeD    =(TH2F*)outputRootD->Get("hCellAmplitude");
	TH1F* hnEventsD          =(TH1F*)outputRootD->Get("hNEvents");
	hCellAmplitudeD->Scale(hnEventsD->GetBinContent(1));
	TH2F* hCellAmplitudeBlockMaskD=(TH2F*)hCellAmplitudeD->Clone("hCellAmplitudeMask");
	TH1F* hCellFlagD         =(TH1F*)outputRootD->Get("fhCellFlag");

	//......................................................
	//..Get the .root file with run-by-run bad channel mask
	//......................................................
	cout<<"**Open file with masking based on run-by-run info within a block: "<<endl;
	TString pathRbR        = Form("./AnalysisOutput/%s/Train_%i/RunByRunSummary105",period.Data(),trainNo);
	TString rootFileNameRbR = Form("GloballyGood_Results.root");
	TFile* outputRootRbR    = TFile::Open(Form("%s/%s",pathRbR .Data(),rootFileNameRbR .Data()));
	if(!outputRootRbR )cout<<"File "<<outputRootRbR ->GetName()<<" does not exist"<<endl;
	else cout<<"file R-b-R: "<<outputRootRbR ->GetName()<<endl;

	TH2F* BadCellsVsRun =(TH2F*)outputRootRbR ->Get("hFlag3vsRun"); //..that is the original bad map - cleaned and filled up for channels with>80% bad runs
	TH1F* hEvtVsRun     =(TH1F*)outputRootRbR ->Get("nEventVsRun");
	TH2F* AmpVsID[nRunsUsed];
    for(Int_t nRun=0;nRun<nRunsUsed;nRun++)
    {
		AmpVsID[nRun] = (TH2F*)outputRootRbR ->Get(Form("hCellAmplitudeRun%i",nRun));
		AmpVsID[nRun]->Scale(hEvtVsRun->GetBinContent(nRun+1));
    }
    //..build the 3/4 blocks
	TH2F* AmpVsIDrBrBlockA = (TH2F*)AmpVsID[0]->Clone("AmpVsIDrBrBlockA");
	TH2F* AmpVsIDrBrBlockB = (TH2F*)AmpVsID[0]->Clone("AmpVsIDrBrBlockB");
	TH2F* AmpVsIDrBrBlockC = (TH2F*)AmpVsID[0]->Clone("AmpVsIDrBrBlockC");
	TH2F* AmpVsIDrBrBlockD = (TH2F*)AmpVsID[0]->Clone("AmpVsIDrBrBlockD");

	//..use only the histogram properties not the content
	AmpVsIDrBrBlockA->Reset();
	AmpVsIDrBrBlockB->Reset();
	AmpVsIDrBrBlockC->Reset();
	AmpVsIDrBrBlockD->Reset();

	for(Int_t nRun=0;nRun<nRunsUsed;nRun++)
    {
		if(nRun<=splitRuns1)                   AmpVsIDrBrBlockA->Add(AmpVsID[nRun]);
		if(nRun>splitRuns1 && nRun<=splitRuns2)AmpVsIDrBrBlockB->Add(AmpVsID[nRun]);
		if(nRun>splitRuns2 && nRun<=splitRuns3)AmpVsIDrBrBlockC->Add(AmpVsID[nRun]);
		if(nRun>splitRuns3)                    AmpVsIDrBrBlockD->Add(AmpVsID[nRun]);
    }

	//--------------------------------------------------
	//..unfortunatley this part is necessary. The amplidude vs. ID figures
	//..are scaled at the end of the bad channel analysis and alathough
	//..we can undo the scaling here by multiplying with nEvents it does
	//..not look the same. It does when one comments the scaling in the bad
	//..channel analysis. So to avaoid intoducing artificial differences
	//..we use the exact same amplidute vs. ID figure.
	AmpVsIDrBrBlockA->Reset();
	AmpVsIDrBrBlockA=(TH2F*)hCellAmplitudeA->Clone("AmpVsIDrBrBlockA");
	AmpVsIDrBrBlockB->Reset();
	AmpVsIDrBrBlockB=(TH2F*)hCellAmplitudeB->Clone("AmpVsIDrBrBlockB");
	AmpVsIDrBrBlockC->Reset();
	AmpVsIDrBrBlockC=(TH2F*)hCellAmplitudeC->Clone("AmpVsIDrBrBlockC");
	AmpVsIDrBrBlockD->Reset();
	AmpVsIDrBrBlockD=(TH2F*)hCellAmplitudeD->Clone("AmpVsIDrBrBlockD");
	//--------------------------------------------------

	//......................................................
	//..mask the bad cells according to both versions
	//......................................................
	Int_t maskA=0;
	Int_t maskB=0;
	Int_t maskC=0;
	Int_t maskD=0;

	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		//......................................................
		//..Part A - analyzed as one merged runblock
		//......................................................
		maskA=0;
		maskB=0;
		maskC=0;
		maskD=0;
		maskA=IsCellMaskedByHand(ic,badcellsBlock1);
		maskB=IsCellMaskedByHand(ic,badcellsBlock2);
		maskC=IsCellMaskedByHand(ic,badcellsBlock3);
		maskD=IsCellMaskedByHand(ic,badcellsBlock4);

		//..mask the bad cells
		for (Int_t amp = 1; amp <= hCellAmplitudeBlockMaskA->GetNbinsX(); amp++)
		{
			if(hCellFlagA->GetBinContent(ic+1)>0 || maskA==1)
			{
				hCellAmplitudeBlockMaskA->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellMerged[0]--;
			}
			if(hCellFlagB->GetBinContent(ic+1)>0 || maskB==1)
			{
				hCellAmplitudeBlockMaskB->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellMerged[1]--;
			}
			if(hCellFlagC->GetBinContent(ic+1)>0 || maskC==1)
			{
				hCellAmplitudeBlockMaskC->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellMerged[2]--;
			}
			if(hCellFlagD->GetBinContent(ic+1)>0 || maskD==1)
			{
				hCellAmplitudeBlockMaskD->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellMerged[3]--;
			}
		}

		//......................................................
		//..Part B - analyzed run-by-run, corrected and masked in block
		//......................................................

		TH1D *htmpCellAllRuns     =BadCellsVsRun->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		Double_t integralBlock1   =htmpCellAllRuns->Integral(0,splitRuns1);
		Double_t integralBlock2   =htmpCellAllRuns->Integral(splitRuns1+1,splitRuns2);
		Double_t integralBlock3   =htmpCellAllRuns->Integral(splitRuns2+1,splitRuns3);
		Double_t integralBlock4   =htmpCellAllRuns->Integral(splitRuns3+1,nRunsUsed);

		//..manually mask cells
		integralBlock1+=IsCellMaskedByHand(ic,badcellsBlock1);
		integralBlock2+=IsCellMaskedByHand(ic,badcellsBlock2);
		integralBlock3+=IsCellMaskedByHand(ic,badcellsBlock3);
		integralBlock4+=IsCellMaskedByHand(ic,badcellsBlock4);

		for(Int_t amp = 1; amp <= AmpVsIDrBrBlockA->GetNbinsX(); amp++)
		{
			if(integralBlock1>0)
			{
				AmpVsIDrBrBlockA->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellRbR[0]--;
			}
			if(integralBlock2>0)
			{
				AmpVsIDrBrBlockB->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellRbR[1]--;
			}
			if(integralBlock3>0)
			{
				AmpVsIDrBrBlockC->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellRbR[2]--;
			}
			if(integralBlock4>0)
			{
				AmpVsIDrBrBlockD->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellRbR[3]--;
			}
		}

		//......................................................
		//..Compare the different channels that are marked in the two versions
		//......................................................
		if(!IsCellMaskedByHand(ic,badcellsBlock1))
		{
			if(hCellFlagA->GetBinContent(ic+1)>0  && integralBlock1==0)vOnlyMaskedInMergedA.push_back(ic);
			if(hCellFlagA->GetBinContent(ic+1)==0 && integralBlock1>0) vOnlyMaskedInRunbRunA.push_back(ic);
		}
		if(!IsCellMaskedByHand(ic,badcellsBlock2))
		{
			if(hCellFlagB->GetBinContent(ic+1)>0  && integralBlock2==0)vOnlyMaskedInMergedB.push_back(ic);
			if(hCellFlagB->GetBinContent(ic+1)==0 && integralBlock2>0) vOnlyMaskedInRunbRunB.push_back(ic);
		}
		if(!IsCellMaskedByHand(ic,badcellsBlock3))
		{
			if(hCellFlagC->GetBinContent(ic+1)>0  && integralBlock3==0)vOnlyMaskedInMergedC.push_back(ic);
			if(hCellFlagC->GetBinContent(ic+1)==0 && integralBlock3>0) vOnlyMaskedInRunbRunC.push_back(ic);
		}
		if(!IsCellMaskedByHand(ic,badcellsBlock4))
		{
			if(hCellFlagD->GetBinContent(ic+1)>0  && integralBlock4==0)vOnlyMaskedInMergedD.push_back(ic);
			if(hCellFlagD->GetBinContent(ic+1)==0 && integralBlock4>0) vOnlyMaskedInRunbRunD.push_back(ic);
		}

	}
	//..merged runblock
	TH1* projMaskedCellsA = hCellAmplitudeBlockMaskA->ProjectionX("MaskedCellsMergedBlockA");
	TH1* projMaskedCellsB = hCellAmplitudeBlockMaskB->ProjectionX("MaskedCellsMergedBlockB");
	TH1* projMaskedCellsC = hCellAmplitudeBlockMaskC->ProjectionX("MaskedCellsMergedBlockC");
	TH1* projMaskedCellsD = hCellAmplitudeBlockMaskD->ProjectionX("MaskedCellsMergedBlockD");

	//..RbR mask
	TH1* projMaskedCellsBlockA = AmpVsIDrBrBlockA->ProjectionX("MaskedCellsRbRBlockA");
	TH1* projMaskedCellsBlockB = AmpVsIDrBrBlockB->ProjectionX("MaskedCellsRbRBlockB");
	TH1* projMaskedCellsBlockC = AmpVsIDrBrBlockC->ProjectionX("MaskedCellsRbRBlockC");
	TH1* projMaskedCellsBlockD = AmpVsIDrBrBlockD->ProjectionX("MaskedCellsRbRBlockD");

	//......................................................
	//..Plot results
	//......................................................
	TCanvas* C1 = new TCanvas("-1-","Bock A and B",900,900);
	C1->Divide(2,2);
	C1->cd(1)->SetLogy();
  	SetHisto(projMaskedCellsA,"","hits/event",0);
	projMaskedCellsA->DrawCopy("hist");
	projMaskedCellsBlockA->SetLineColor(8);
	projMaskedCellsBlockA->DrawCopy("same hist");
	C1->cd(2);
	projMaskedCellsA->Divide(projMaskedCellsBlockA);
  	SetHisto(projMaskedCellsA,"","Merged block/Run-by-Run block",0);
	projMaskedCellsA->DrawCopy("hist");
	projMaskedCellsA->Multiply(projMaskedCellsBlockA);

	C1->cd(3)->SetLogy();
  	SetHisto(projMaskedCellsB,"","hits/event",0);
	projMaskedCellsB->DrawCopy("hist");
	projMaskedCellsBlockB->SetLineColor(8);
	projMaskedCellsBlockB->DrawCopy("same hist");
	C1->cd(4);
	projMaskedCellsB->Divide(projMaskedCellsBlockB);
  	SetHisto(projMaskedCellsB,"","Merged block/Run-by-Run block",0);
	projMaskedCellsB->DrawCopy("hist");
	projMaskedCellsB->Multiply(projMaskedCellsBlockB);

	TCanvas* C2 = new TCanvas("-2-","Block C and D",900,900);
	C2->Divide(2,2);
	C2->cd(1)->SetLogy();
  	SetHisto(projMaskedCellsC,"","hits/event",0);
	projMaskedCellsC->DrawCopy("hist");
	projMaskedCellsBlockC->SetLineColor(8);
	projMaskedCellsBlockC->DrawCopy("same hist");
	C2->cd(2);
	projMaskedCellsC->Divide(projMaskedCellsBlockC);
  	SetHisto(projMaskedCellsC,"","Merged block/Run-by-Run block",0);
	projMaskedCellsC->DrawCopy("hist");
	projMaskedCellsC->Multiply(projMaskedCellsBlockC);

	C2->cd(3)->SetLogy();
  	SetHisto(projMaskedCellsD,"","hits/event",0);
	projMaskedCellsD->DrawCopy("hist");
	projMaskedCellsBlockD->SetLineColor(8);
	projMaskedCellsBlockD->DrawCopy("same hist");
	C2->cd(4);
	projMaskedCellsD->Divide(projMaskedCellsBlockD);
  	SetHisto(projMaskedCellsD,"","Merged block/Run-by-Run block",0);
	projMaskedCellsD->DrawCopy("hist");
	projMaskedCellsD->Multiply(projMaskedCellsBlockD);

	TCanvas* C3 = new TCanvas("-3-","Merge A, B, C, and D",900,900);
	C3->Divide(2,2);
	C3->cd(1)->SetLogz();
  	SetHisto(hCellAmplitudeBlockMaskA,"","cell ID",0);
  	hCellAmplitudeBlockMaskA->DrawCopy("colz");
	C3->cd(2)->SetLogz();
  	SetHisto(hCellAmplitudeBlockMaskB,"","cell ID",0);
  	hCellAmplitudeBlockMaskB->DrawCopy("colz");
	C3->cd(3)->SetLogz();
  	SetHisto(hCellAmplitudeBlockMaskC,"","cell ID",0);
  	hCellAmplitudeBlockMaskC->DrawCopy("colz");
	C3->cd(4)->SetLogz();
  	SetHisto(hCellAmplitudeBlockMaskD,"","cell ID",0);
  	hCellAmplitudeBlockMaskD->DrawCopy("colz");

	TCanvas* C4 = new TCanvas("-4-","Run-by-Run A, B, C, and D",900,900);
	C4->Divide(2,2);
	C4->cd(1)->SetLogz();
  	SetHisto(AmpVsIDrBrBlockA,"","cell ID",0);
  	AmpVsIDrBrBlockA->DrawCopy("colz");
	C4->cd(2)->SetLogz();
  	SetHisto(AmpVsIDrBrBlockB,"","cell ID",0);
  	AmpVsIDrBrBlockB->DrawCopy("colz");
	C4->cd(3)->SetLogz();
  	SetHisto(AmpVsIDrBrBlockC,"","cell ID",0);
  	AmpVsIDrBrBlockC->DrawCopy("colz");
	C4->cd(4)->SetLogz();
  	SetHisto(AmpVsIDrBrBlockD,"","cell ID",0);
  	AmpVsIDrBrBlockD->DrawCopy("colz");

 	//......................................................
	//..Print out compared cells and plot the spectra
	//......................................................
	projMaskedCellsA->Scale(1.0/nBadCellMerged[0]);
	projMaskedCellsB->Scale(1.0/nBadCellMerged[1]);
	projMaskedCellsC->Scale(1.0/nBadCellMerged[2]);
	projMaskedCellsD->Scale(1.0/nBadCellMerged[3]);
	projMaskedCellsBlockA->Scale(1.0/nBadCellRbR[0]);
	projMaskedCellsBlockB->Scale(1.0/nBadCellRbR[1]);
	projMaskedCellsBlockC->Scale(1.0/nBadCellRbR[2]);
	projMaskedCellsBlockD->Scale(1.0/nBadCellRbR[3]);

	cout<<0<<"-"<<splitRuns1<<" -> "<<splitRuns1<<" runs"<<endl;
	cout<<splitRuns1+1<<"-"<<splitRuns2<<" -> "<<splitRuns2-splitRuns1<<" runs"<<endl;
	cout<<splitRuns2+1<<"-"<<splitRuns3<<" -> "<<splitRuns3-splitRuns2<<" runs"<<endl;
	cout<<splitRuns3+1<<"-"<<nRunsUsed<<" -> "<<nRunsUsed-splitRuns3<<" runs"<<endl;

	cout<<"o Run block A ("<<"0-"<<splitRuns1<<")"<<endl;
	cout<<"  Cells masked in Merged version and not in Run-b-run version ("<<vOnlyMaskedInMergedA.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInMergedA.size();i++)
	{
		cout<<vOnlyMaskedInMergedA.at(i)<<","<<flush;
	}
	cout<<endl;
	CreateCellCompPDF(AmpVsIDrBrBlockA,vOnlyMaskedInMergedA,projMaskedCellsA,projMaskedCellsBlockA,"./cOnlyMergedBlockA.pdf");
	cout<<"  Cells masked in Run-b-run version and not in Merged version ("<<vOnlyMaskedInRunbRunA.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInRunbRunA.size();i++)
	{
		cout<<vOnlyMaskedInRunbRunA.at(i)<<","<<flush;
	}
	cout<<endl;
	CreateCellCompPDF(hCellAmplitudeBlockMaskA,vOnlyMaskedInRunbRunA,projMaskedCellsA,projMaskedCellsBlockA,"./cOnlyRbRBlockA.pdf");

	cout<<"o Run block B ("<<splitRuns1+1<<"-"<<splitRuns2<<")"<<endl;
	cout<<"  Cells masked in Merged version and not in Run-b-run version ("<<vOnlyMaskedInMergedB.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInMergedB.size();i++)
	{
		cout<<vOnlyMaskedInMergedB.at(i)<<","<<flush;
	}
	cout<<endl;
	CreateCellCompPDF(AmpVsIDrBrBlockB,vOnlyMaskedInMergedB,projMaskedCellsB,projMaskedCellsBlockB,"./cOnlyMergedBlockB.pdf");
	cout<<"  Cells masked in Run-b-run version and not in Merged version ("<<vOnlyMaskedInRunbRunB.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInRunbRunB.size();i++)
	{
		cout<<vOnlyMaskedInRunbRunB.at(i)<<","<<flush;
	}
	cout<<endl;
	CreateCellCompPDF(hCellAmplitudeBlockMaskB,vOnlyMaskedInRunbRunB,projMaskedCellsB,projMaskedCellsBlockB,"./cOnlyRbRBlockB.pdf");

	cout<<"o Run block C ("<<splitRuns2+1<<"-"<<splitRuns3<<")"<<endl;
	cout<<"  Cells masked in Merged version and not in Run-b-run version ("<<vOnlyMaskedInMergedC.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInMergedC.size();i++)
	{
		cout<<vOnlyMaskedInMergedC.at(i)<<","<<flush;
	}
	cout<<endl;
	CreateCellCompPDF(AmpVsIDrBrBlockC,vOnlyMaskedInMergedC,projMaskedCellsC,projMaskedCellsBlockC,"./cOnlyMergedBlockC.pdf");
	cout<<"  Cells masked in Run-b-run version and not in Merged version ("<<vOnlyMaskedInRunbRunC.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInRunbRunC.size();i++)
	{
		cout<<vOnlyMaskedInRunbRunC.at(i)<<","<<flush;
	}
	cout<<endl;
	CreateCellCompPDF(hCellAmplitudeBlockMaskC,vOnlyMaskedInRunbRunC,projMaskedCellsC,projMaskedCellsBlockC,"./cOnlyRbRBlockC.pdf");

	cout<<"o Run block D ("<<splitRuns3+1<<"-"<<nRunsUsed<<")"<<endl;
	cout<<"  Cells masked in Merged version and not in Run-b-run version ("<<vOnlyMaskedInMergedD.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInMergedD.size();i++)
	{
		cout<<vOnlyMaskedInMergedD.at(i)<<","<<flush;
	}
	cout<<endl;
	CreateCellCompPDF(AmpVsIDrBrBlockD,vOnlyMaskedInMergedD,projMaskedCellsD,projMaskedCellsBlockD,"./cOnlyMergedBlockD.pdf");
	cout<<"  Cells masked in Run-b-run version and not in Merged version ("<<vOnlyMaskedInRunbRunD.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInRunbRunD.size();i++)
	{
		cout<<vOnlyMaskedInRunbRunD.at(i)<<","<<flush;
	}
	cout<<endl;
	CreateCellCompPDF(hCellAmplitudeBlockMaskD,vOnlyMaskedInRunbRunD,projMaskedCellsD,projMaskedCellsBlockD,"./cOnlyRbRBlockD.pdf");

 	//......................................................
	//..build two dimensional histogram with cells rejected from
	//..the one or the other method
	//......................................................
	Plot2DCells("A",244917,vOnlyMaskedInRunbRunA,vOnlyMaskedInMergedA);
	Plot2DCells("B",244917,vOnlyMaskedInRunbRunB,vOnlyMaskedInMergedB);
	Plot2DCells("C",244917,vOnlyMaskedInRunbRunC,vOnlyMaskedInMergedC);
	Plot2DCells("D",244917,vOnlyMaskedInRunbRunD,vOnlyMaskedInMergedD);

	vOnlyMaskedInRunbRunA.insert( vOnlyMaskedInRunbRunA.end(), vOnlyMaskedInRunbRunB.begin(), vOnlyMaskedInRunbRunB.end() );
	vOnlyMaskedInRunbRunA.insert( vOnlyMaskedInRunbRunA.end(), vOnlyMaskedInRunbRunC.begin(), vOnlyMaskedInRunbRunC.end() );
	vOnlyMaskedInRunbRunA.insert( vOnlyMaskedInRunbRunA.end(), vOnlyMaskedInRunbRunD.begin(), vOnlyMaskedInRunbRunD.end() );
	vOnlyMaskedInMergedA.insert( vOnlyMaskedInMergedA.end(), vOnlyMaskedInMergedB.begin(), vOnlyMaskedInMergedB.end() );
	vOnlyMaskedInMergedA.insert( vOnlyMaskedInMergedA.end(), vOnlyMaskedInMergedC.begin(), vOnlyMaskedInMergedC.end() );
	vOnlyMaskedInMergedA.insert( vOnlyMaskedInMergedA.end(), vOnlyMaskedInMergedD.begin(), vOnlyMaskedInMergedD.end() );
	Plot2DCells("sum",244917,vOnlyMaskedInRunbRunA,vOnlyMaskedInMergedA);

}
//
//
//
//
void Plot2DCells(TString Block, Int_t runNo, std::vector<Int_t> cellVectorRbR, std::vector<Int_t> cellVectorMerge)
{
	//......................................................
	//..Initialize EMCal/DCal geometry
	 AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(runNo);
	fCaloUtils->AccessGeometry(aod);
	//......................................................
	//..setings for the 2D histogram
	Int_t nModules=fCaloUtils->GetEMCALGeometry()->GetNumberOfSuperModules();
	Int_t fNMaxColsAbs = 2*48;
	Int_t fNMaxRowsAbs = Int_t (nModules/2)*24; //multiply by number of supermodules

	TH2F *plot2D_RbR   = new TH2F(Form("Block%s_MaskedRbR",Block.Data()),Form("Block%s_MaskedRbR",Block.Data()),fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	plot2D_RbR->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_RbR->GetYaxis()->SetTitle("cell row (#phi direction)");
	TH2F *plot2D_Merge = new TH2F(Form("Block%s_MaskedMerge",Block.Data()),Form("Block%s_MaskedMerge",Block.Data()),fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	plot2D_Merge->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_Merge->GetYaxis()->SetTitle("cell row (#phi direction)");

	Int_t cellColumn=0,cellRow=0;
	Int_t cellColumnAbs=0,cellRowAbs=0;
	Int_t trash;

	for(Int_t i = 0; i < (Int_t)cellVectorRbR.size(); i++)
	{
		Int_t cell=cellVectorRbR.at(i);
		//..Do that only for cell ids also accepted by the code
		if(!fCaloUtils->GetEMCALGeometry()->CheckAbsCellId(cell))continue;

		//..Get Row and Collumn for cell ID c
		fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(cell,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
		if(cellColumnAbs> fNMaxColsAbs || cellRowAbs>fNMaxRowsAbs)
		{
			cout<<"Problem! wrong calculated number of max col and max rows"<<endl;
			cout<<"current col: "<<cellColumnAbs<<", max col"<<fNMaxColsAbs<<endl;
			cout<<"current row: "<<cellRowAbs<<", max row"<<fNMaxRowsAbs<<endl;
		}
		plot2D_RbR->Fill(cellColumnAbs,cellRowAbs);
	}
	for(Int_t i = 0; i < (Int_t)cellVectorMerge.size(); i++)
	{
		Int_t cell=cellVectorMerge.at(i);
		//..Do that only for cell ids also accepted by the code
		if(!fCaloUtils->GetEMCALGeometry()->CheckAbsCellId(cell))continue;

		//..Get Row and Collumn for cell ID c
		fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(cell,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
		if(cellColumnAbs> fNMaxColsAbs || cellRowAbs>fNMaxRowsAbs)
		{
			cout<<"Problem! wrong calculated number of max col and max rows"<<endl;
			cout<<"current col: "<<cellColumnAbs<<", max col"<<fNMaxColsAbs<<endl;
			cout<<"current row: "<<cellRowAbs<<", max row"<<fNMaxRowsAbs<<endl;
		}
		plot2D_Merge->Fill(cellColumnAbs,cellRowAbs,1);
	}
	//. . . . . . . . . . . . . . . . . . . .
	TCanvas *c1 = new TCanvas(Form("2DMapForBlock%s",Block.Data()),Form("2DMapForBlock%s",Block.Data()),900,500);
	c1->ToggleEventStatus();
	c1->Divide(2);
	c1->cd(1);
	plot2D_RbR->Draw("colz");
	c1->cd(2);
	plot2D_Merge->Draw("colz");
}
//
//________________________________________________________________________
void BuildMaxMinHisto(TH1D* inHisto, TH1D* minHist,TH1D* maxHist)
{
	Double_t ref;
	Double_t min;
	Double_t max;
	for(Int_t bin=0;bin<50;bin++)
	{
		ref = inHisto->GetBinContent(bin);
	    max = maxHist->GetBinContent(bin);
	    min = minHist->GetBinContent(bin);
	    if((ref!=0 && ref<min) || min==0)minHist->SetBinContent(bin,ref);
	    if(ref>max)maxHist->SetBinContent(bin,ref);
	}
}
//
//________________________________________________________________________
Bool_t IsItReallyBadRatio(TH1D* minHistoRatio,TH1D* maxHistoRatio,TH1D* meanHistoRatio, TString& crit)
{
	//minHistoRatio  ->Smooth();
	//maxHistoRatio  ->Smooth();
	//meanHistoRatio ->Smooth();

	//..do the check only until 1 GeV because
	//..then we are running out of statistic
	/*Int_t OneGeVBin  =minHistoRatio->FindBin(1.47);  //=30bins
	Int_t HalfGeVBin =minHistoRatio->FindBin(0.73);  //=15bins
	Int_t ThirdGeVBinA =minHistoRatio->FindBin(0.35);
	Int_t ThirdGeVBinB =minHistoRatio->FindBin(0.7);
	Int_t ThirdGeVBinC =minHistoRatio->FindBin(1.05);
	*/
    Double_t mBlock1=0,mBlock2=0,mBlock3=0,mBlock4=0,mBlock5=0;
    Double_t zBlock1=0,zBlock2=0,zBlock3=0,zBlock4=0,zBlock5=0;

	Double_t minMean=0,    maxMean=0,    meanMean=0;
	Double_t zminMean=0,   zmaxMean=0,   zmeanMean=0;
/*	Double_t meanMean1=0,  meanMean2=0;
	Double_t meanMeanA3=0, meanMeanB3=0, meanMeanC3=0;
	Int_t zeroBinsHalf1=0, zeroBinsHalf2=0;
	Int_t zeroBinsThird1=0, zeroBinsThird2=0, zeroBinsThird3=0;
*/
	//for(Int_t bin=0;bin<30;bin++)
	for(Int_t bin=3;bin<53;bin++)
	{
		if(bin<33)
		{
			minMean  += minHistoRatio ->GetBinContent(bin);
			maxMean  += maxHistoRatio ->GetBinContent(bin);
			meanMean += meanHistoRatio->GetBinContent(bin);
			if(minHistoRatio->GetBinContent(bin)==0)zminMean++;
			if(maxHistoRatio->GetBinContent(bin)==0)zmaxMean++;
			if(meanHistoRatio->GetBinContent(bin)==0)zmeanMean++;
		}

/*		if(bin<18)      meanMean1+=meanHistoRatio->GetBinContent(bin);
		else if(bin<33) meanMean2+=meanHistoRatio->GetBinContent(bin);

*/
		/*		if(bin<10)     meanMeanA3+=meanHistoRatio->GetBinContent(bin+1);
		else if(bin<20)meanMeanB3+=meanHistoRatio->GetBinContent(bin+1);
		else           meanMeanC3+=meanHistoRatio->GetBinContent(bin+1);
		 */
		if(bin<13)      mBlock1+=meanHistoRatio->GetBinContent(bin);
		else if(bin<23) mBlock2+=meanHistoRatio->GetBinContent(bin);
		else if(bin<33) mBlock3+=meanHistoRatio->GetBinContent(bin);
		else if(bin<43) mBlock4+=meanHistoRatio->GetBinContent(bin);
		else if(bin<53) mBlock5+=meanHistoRatio->GetBinContent(bin);
		//..count zero bins
		if(meanHistoRatio->GetBinContent(bin)==0)
		{
			if(bin<13)      zBlock1++;
			else if(bin<23) zBlock2++;
			else if(bin<33) zBlock3++;
			else if(bin<43) zBlock4++;
			else if(bin<53) zBlock5++;
		}
	/*	//..correct for zero bins
		if(meanHistoRatio->GetBinContent(bin+1)==0 && bin<15) zeroBinsHalf1++;
		else if(meanHistoRatio->GetBinContent(bin+1)==0)      zeroBinsHalf2++;

		if(meanHistoRatio->GetBinContent(bin+1)==0 && bin<10)     zeroBinsThird1++;
		else if(meanHistoRatio->GetBinContent(bin+1)==0 && bin<20)zeroBinsThird2++;
		else if(meanHistoRatio->GetBinContent(bin+1)==0)          zeroBinsThird3++;
	*/}
	//.......................................
	//..correct for zero bins
	//.......................................
/*
	if(zeroBinsHalf1!=0)meanMean1=meanMean1/(1.0-1.0*zeroBinsHalf1/15);
	if(zeroBinsHalf2!=0)meanMean2=meanMean2/(1.0-1.0*zeroBinsHalf2/15);

	if(zeroBinsThird1!=0)meanMeanA3=meanMeanA3/(1.0-1.0*zeroBinsThird1/10);
	if(zeroBinsThird2!=0)meanMeanB3=meanMeanB3/(1.0-1.0*zeroBinsThird2/10);
	if(zeroBinsThird3!=0)meanMeanC3=meanMeanC3/(1.0-1.0*zeroBinsThird3/10);
*/
	if(zBlock1<10 && zBlock1!=0)mBlock1=mBlock1/(1.0-1.0*zBlock1/10);
	if(zBlock2<10 && zBlock2!=0)mBlock2=mBlock2/(1.0-1.0*zBlock2/10);
	if(zBlock3<10 && zBlock3!=0)mBlock3=mBlock3/(1.0-1.0*zBlock3/10);
	if(zBlock4<10 && zBlock4!=0)mBlock4=mBlock4/(1.0-1.0*zBlock4/10);
	if(zBlock5<10 && zBlock5!=0)mBlock5=mBlock5/(1.0-1.0*zBlock5/10);

	if(zminMean!=0)minMean  =minMean/(1.0-1.0*zminMean/30);
	if(zmaxMean!=0)maxMean  =maxMean/(1.0-1.0*zmaxMean/30);
	if(zmeanMean!=0)meanMean=meanMean/(1.0-1.0*zmeanMean/30);

	//..if more than half of the bins in the block were 0 exclude block
	if(zBlock1>5)mBlock1=0;
	if(zBlock2>5)mBlock2=0;
	if(zBlock3>5)mBlock3=0;
	if(zBlock4>5)mBlock4=0;
	if(zBlock5>5)mBlock5=0;
	//.......................................
	//..check criteria
	//.......................................

	//..bad channel is 5times lower than max distr.
	crit = "spectr. too low";
	if(meanMean/30.0>5) return 1;
	if(minMean/30.0 >5) return 1;   //5 times lower than the lowest run
	//..bad channel is 5times higher than max distr.
	crit = "spectr. too high";
	if(meanMean/30.0<0.2) return 1;
	if(maxMean/30.0 <0.2) return 1; //5 times higher than the highest run

	//..if there is a slope down (should be more than 10% decrease)
	crit = "slope down";
	Int_t down=0;
	if(mBlock1!=0 && mBlock2!=0 && mBlock1>mBlock2 && mBlock1>1.4*mBlock2)down++;
	if(mBlock2!=0 && mBlock3!=0 && mBlock2>mBlock3 && mBlock2>1.4*mBlock3)down++;
	if(mBlock3!=0 && mBlock4!=0 && mBlock3>mBlock4 && mBlock3>1.4*mBlock4)down++;
	if(mBlock4!=0 && mBlock5!=0 && mBlock4>mBlock5 && mBlock4>1.4*mBlock5)down++;
	if(down>=2)return 1; //..means at least three blocks have to be staggered

	crit = "slope up";
	Int_t up=0;
	if(mBlock1!=0 && mBlock2 !=0 && mBlock1<mBlock2 && 1.4*mBlock1<mBlock2)up++;
	if(mBlock2!=0 && mBlock3 !=0 && mBlock2<mBlock3 && 1.4*mBlock2<mBlock3)up++;
	if(mBlock3!=0 && mBlock4 !=0 && mBlock3<mBlock4 && 1.4*mBlock3<mBlock4)up++;
	if(mBlock4!=0 && mBlock5 !=0 && mBlock4<mBlock5 && 1.4*mBlock4<mBlock5)up++;
	if(up>=2)return 1; //..means at least three blocks have to be staggered

	//..if there is a step at 2.1 GeV
    crit = "step 2.1 GeV";
    if(mBlock4!=0 && mBlock5!=0 && mBlock4>mBlock5 && mBlock4>20*mBlock5) return 1;

    //..if there is a steep step at 1.1 GeV (can only be performed if this is not dominated by "zero" bins)
//    crit = "step 1.1";
//    if(zeroBinsThird3<5 && meanMeanB3>meanMeanC3 && meanMeanB3>3.5*meanMeanC3) return 1;
	//..is good
    crit = "";
	return 0;
}
//
//________________________________________________________________________
void PlotLowFractionCells(TString pdfName, std::vector<Int_t> cellVector,TH2F* badVsCell[],Int_t nRuns,TH2F* ampID[],TH1D* hCellGoodMean[])
{
	gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //..to supress a lot of standard output

	Int_t nRescuableChannels=cellVector.size();
	Int_t totalperCv = 16;
	Int_t nPad = TMath::Sqrt(totalperCv);
	Int_t nCv  = nRescuableChannels/totalperCv+1;
	Int_t lastGood=0;

	TLatex* text = new TLatex(0.45,0.6,"*Indeed bad*");//..
	text->SetTextSize(0.07);
	text->SetTextColor(kAzure-8);
	text->SetNDC();

	TH1D *maxHisto  = ampID[0]->ProjectionX("hMaxCells", 1, 1);
	TH1D *minHisto  = ampID[0]->ProjectionX("hMinCells", 1, 1);
	TH1D* hgoodMean = ampID[0]->ProjectionX("hMeanofRuns", 1, 1);
	if(nCv<1)nCv=1;

	cout<<"    o create: "<<nCv<<" Canvases with "<<nPad*nPad<<" pads"<<endl;
	//..to compare specific cells over the runs
	TCanvas **cComp     = new TCanvas*[nCv];
	TCanvas **cCompDiv  = new TCanvas*[nCv];
	for(Int_t i=0;i<nCv;i++)
	{
		cComp[i]    = new TCanvas(TString::Format("CompareGood%d", i), TString::Format("V) Candidates (%d/%d)", i+1, nCv), 1000,750);
		cComp[i]    ->Divide(nPad,nPad,0.001,0.001);
		cCompDiv[i] = new TCanvas(TString::Format("CompareGood Ratio%d", i), TString::Format("V) Candidates Ratio (%d/%d)", i+1, nCv), 1000,750);
		cCompDiv[i] ->Divide(nPad,nPad,0.001,0.001);
	}

	Int_t notBadCounter=0;
	//TH1F** hCellSpectr = new TH1F*[nRescuableChannels];
	for(Int_t cell=0;cell< (Int_t)cellVector.size();cell++)
	{
		if(cell%400==0)cout<<"cell No."<<cell<<endl;
		if(cell%20==0) cout<<"."<<flush;
		maxHisto->Reset();
		minHisto->Reset();
		hgoodMean->Reset();
		TH1D* declaredBad;
		std::vector<TH1D*> badHistVector;
		Int_t badRun=-1;
		for(Int_t i = 0 ; i < nRuns ; i++)
		{
			TH1D *htmpCell = ampID[i]->ProjectionX(TString::Format("hIDProj_cell%dRun%i", cellVector.at(cell),i), cellVector.at(cell)+1, cellVector.at(cell)+1);
			htmpCell->SetLineColor(kGray+1);
			if(badVsCell[2]->GetBinContent(cellVector.at(cell)+1,i+1)==1)
			{
				htmpCell->SetLineColor(2);
				htmpCell->SetFillColor(2);
				htmpCell->SetFillStyle(3002);
				declaredBad = (TH1D*)htmpCell->Clone("saveForLater");
				badHistVector.push_back(htmpCell);
				badRun=i;
				if(htmpCell->Integral()==0)cout<<"cell "<<cell<<" is dead for run: "<<i<<endl;
			}
			else
			{
				BuildMaxMinHisto(htmpCell,minHisto,maxHisto);
				hgoodMean->Add(htmpCell);
			}
			//..go to the last pad and draw the mean of all good cell distribution
			//..for all the tested runs
			if(i==0)
			{
				cComp[nCv-1]   ->cd(nPad*nPad)->SetLogy();
				SetHisto(hCellGoodMean[i],"","Number of Hits",0);
				hCellGoodMean[i]->GetXaxis()->SetRangeUser(0,3);
				hCellGoodMean[i]->Draw("hist");
			}
			else
			{
				cComp[nCv-1]   ->cd(nPad*nPad)->SetLogy();
				hCellGoodMean[i]->Draw("same hist");
			}
			//..go to pads and draw cell in all runs
			lastGood=(cell-notBadCounter)/totalperCv;//..you can overwrite good canvases to lessen the amount of canvases
			if(i==0)
			{
				cComp[(cell-notBadCounter)/totalperCv]->cd(((cell-notBadCounter)%totalperCv)+1)->SetLogy();
				SetHisto(htmpCell,Form("Energy of cell %i",cellVector.at(cell)),"Number of Hits",0);
				htmpCell->GetXaxis()->SetRangeUser(0,3);
				htmpCell->Draw("hist");
			}
			else
			{
				cComp[(cell-notBadCounter)/totalperCv]->cd(((cell-notBadCounter)%totalperCv)+1)->SetLogy();
				htmpCell->DrawCopy("same hist");
			}
		}//..end of loop over runs
		hgoodMean->Scale(1.0/nRuns);

		hCellGoodMean[badRun]->DrawCopy("same hist"); //..draw the mean of all good cells for the run where the cell was bad
		minHisto->SetLineColor(1);
		maxHisto->SetLineColor(1);
		minHisto->DrawCopy("same");     //..draw the combined minimum of that cell for all the runs
		maxHisto->DrawCopy("same");     //..draw the combined maximum of that cell for all the runs

		//..Draw bad again
		for(Int_t j=0;j< (Int_t)badHistVector.size();j++)
		{
			badHistVector[j]->SetLineColor(2);
			badHistVector[j]->DrawCopy("same");
		}

		TLegend *leg = new TLegend(0.35,0.65,0.75,0.85);
		leg->AddEntry(hCellGoodMean[badRun],"mean of good cells","l");
		leg->AddEntry(declaredBad,"Cell in the ''bad'' run","l");
		leg->AddEntry(maxHisto,"max and min values","l");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.07);
		leg->Draw("same");

		//. . . . . . . . . . . . . . . . . . . . . . . . .
		//..Fill the ratio canvas
		cCompDiv[(cell-notBadCounter)/totalperCv]->cd(((cell-notBadCounter)%totalperCv)+1);

		//..Draw bad again
		Int_t runTry=0;
		Bool_t badRatio=0;
		TString failCrit;
		for(Int_t j=0;j< (Int_t)badHistVector.size();j++)
		{
			TH1D* minHistoCopy =(TH1D*)minHisto->Clone("minHistoCopy");
			TH1D* maxHistoCopy =(TH1D*)maxHisto->Clone("maxHistoCopy");
			TH1D* meanHistoCopy=(TH1D*)hgoodMean->Clone("meanHistoCopy");
			minHistoCopy ->Divide(badHistVector[j]);
			maxHistoCopy ->Divide(badHistVector[j]);
			meanHistoCopy->Divide(badHistVector[j]);

			SetHisto(maxHistoCopy,Form("Energy of cell %i",cellVector.at(cell)),"max,min,mean / bad cell",0);
			maxHistoCopy->GetXaxis()->SetRangeUser(0,3);
			maxHistoCopy->DrawCopy("hist");
			minHistoCopy->DrawCopy("same hist");
			meanHistoCopy->SetLineColor(kBlue-8);
			meanHistoCopy->DrawCopy("same hist");
			//cout<<"cell: "<<cell<<endl;
			badRatio = IsItReallyBadRatio(minHistoCopy,maxHistoCopy,meanHistoCopy,failCrit);
			runTry=j;
			if(badRatio==1)break; //if its bad for one of the runs its enough
		}
		if(badRatio==1)
		{
			gPad->SetFrameFillColor(kRed-10);
			text->DrawLatexNDC(0.45,0.8,Form("#splitline{(%i)Excluded: }{%s}",runTry,(const char*)failCrit));;
		}
		else
		{
			//de-mask cells not declared as bad. (re-inclusion)
			for(Int_t j = 0 ; j < nRuns ; j++)
			{
				badVsCell[0]->SetBinContent(cellVector.at(cell)+1,j+1,0);
				badVsCell[1]->SetBinContent(cellVector.at(cell)+1,j+1,0);
				badVsCell[2]->SetBinContent(cellVector.at(cell)+1,j+1,0);
			}
		}
		if(cell==(Int_t)cellVector.size()-1)text->SetText(0.45,0.8,"test");

	}
	cout<<endl;

	//..plot the canvases of cells into a .pdf file
	for(Int_t can=0;can<nCv;can++)
	{
		TString internal_pdfName1=pdfName+"Low.pdf";
		TString internal_pdfName2=pdfName+"High.pdf";
		TString internal_pdfName3=pdfName+"Ratio.pdf";
		if(can==0)
		{
			//..first pad
			internal_pdfName1 +="(";
			internal_pdfName2 +="(";
			internal_pdfName3 +="(";
			cComp[can]    ->Print(internal_pdfName1.Data());
			cCompDiv[can] ->Print(internal_pdfName3.Data());
		}
		else if(can==(nCv-1))
		{
			//..last pad
			internal_pdfName1 +=")";
			internal_pdfName2 +=")";
			internal_pdfName3 +=")";
			cComp[can]    ->Print(internal_pdfName1.Data());
			cCompDiv[can] ->Print(internal_pdfName3.Data());
		}
		else
		{
			//..all pads in between
			cComp[can]    ->Print(internal_pdfName1.Data());
			cCompDiv[can] ->Print(internal_pdfName3.Data());
		}
	}
	for(Int_t i=lastGood+1;i<nCv-1;i++)
	{
		cout<<"last good "<<lastGood<<endl;
		cout<<"nCv "<<nCv<<endl;
		cout<<"round "<<i<<endl;
		delete cComp[i];
	}
	cout<<endl;

}
///
/// compresses the bad cell histogram
/// to delete entries of good cells (cells that are good in EVERY run) (only 10-15% of cells are bad)
//________________________________________________________________________
TH2F* CompressHistogram(TH2 *Histo,Int_t totalCells, Int_t badCells,std::vector<Int_t> runIdVec)
{
	TH2F* cpmpressed = new TH2F(Form("%s_Comp",Histo->GetName()),Form("%s_Comp",Histo->GetName()),badCells+2, 0,badCells+2,runIdVec.size(),0,runIdVec.size());

	Histo->GetXaxis()->UnZoom();
	TH1D *htmpCell = Histo->ProjectionX(TString::Format("%s_proj",Histo->GetName()),0,runIdVec.size());
	Int_t sumRun=0,newHistoBin=0;
	for(Int_t icell = 0; icell < totalCells ; icell++)
	{
		sumRun = htmpCell->GetBinContent(icell+1);
		//cout<<"enties cell("<<icell<<"): "<<sumRun<<endl;
		//..Fill non zero entries into the new histogram
		if(sumRun>0)
		{
			newHistoBin++;
			if(newHistoBin>badCells)cout<<"PROBLEM"<<endl;
			for(Int_t iRun = 0; iRun < (Int_t)runIdVec.size() ; iRun++)
			{
				cpmpressed->Fill(newHistoBin,iRun,Histo->GetBinContent(icell+1,iRun+1));
				//cout<<"fill bin "<<newHistoBin<<" with value: 2"<<endl;
			}
		}
	}
	for(Int_t i=0;i<(Int_t)runIdVec.size();i++)
	{
		cpmpressed->GetYaxis()->SetBinLabel(i+1,Form("%i",runIdVec.at(i)));
	}
	return cpmpressed;
}
///
/// A vertical line that can be plotted
///
//________________________________________________________________________
void PlotHorLineRange(Double_t y_val, Double_t xLow, Double_t xHigh, Int_t Line_Col)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetY1(y_val);
    Zero_line -> SetY2(y_val);
    Zero_line -> SetX1(xLow);
    Zero_line -> SetX2(xHigh);
    //cout << "x_val = " << x_val << ", Bin = " << Histo->FindBin(x_val) << ", Y2 = " << Histo->GetBinContent(Histo->FindBin(x_val)) << endl;
    Zero_line -> SetLineWidth(2);
    Zero_line -> SetLineStyle(9);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
///
/// Function to set TH1 histograms to a similar style
///
//________________________________________________________________________
void SetHisto(TH2 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto)
{
  //ELI
	Histo->SetStats(0);
	Histo->SetTitle("");
	if(longhisto==0)
	{
		Histo->GetYaxis()->SetTitleOffset(1.4);
		Histo->GetXaxis()->SetTitleOffset(1.4);
		Histo->GetXaxis()->SetLabelSize(0.05);
		Histo->GetYaxis()->SetLabelSize(0.05);
		Histo->GetXaxis()->SetTitleSize(0.045);
		Histo->GetYaxis()->SetTitleSize(0.045);
	}
	//..these are the run number vs. ID
	if(longhisto==1)
	{
		Histo->GetYaxis()->SetTitleOffset(0.6);
		Histo->GetXaxis()->SetTitleOffset(0.8);
		Histo->GetYaxis()->SetLabelOffset(0.002);
//		Histo->GetXaxis()->SetLabelSize(0.07);
//		Histo->GetYaxis()->SetLabelSize(0.07);
		Histo->GetXaxis()->SetLabelSize(0.025);
		Histo->GetYaxis()->SetLabelSize(0.015);
		Histo->GetXaxis()->SetTitleSize(0.03);
		Histo->GetYaxis()->SetTitleSize(0.04);
		Histo->GetYaxis()->SetTickLength(0.02);
	}
	//Histo->GetXaxis()->CenterTitle();
	//Histo->GetYaxis()->CenterTitle();

	if(longhisto==1)
	{
		Histo->GetXaxis()->SetNdivisions(520);
		//Histo->GetYaxis()->SetNdivisions(10);
	}
	else
	{
		Histo->GetXaxis()->SetNdivisions(505);
		Histo->GetYaxis()->SetNdivisions(505);
	}

	//..make nice font
	Histo->GetXaxis()->SetLabelFont(42);
	Histo->GetYaxis()->SetLabelFont(42);
	Histo->GetXaxis()->SetTitleFont(42);
	Histo->GetYaxis()->SetTitleFont(42);
	if(Xtitel!="")Histo->GetXaxis()->SetTitle(Xtitel);
	if(Ytitel!="")Histo->GetYaxis()->SetTitle(Ytitel);
	Histo->SetLineColor(1);
	Histo->SetMarkerColor(1);
	Histo->SetMarkerStyle(20);
	Histo->SetMarkerSize(0.5);
}
///
/// Funtion to set TH1 histograms to a similar style
///
//________________________________________________________________________
void SetHisto(TH1 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto)
{
	Histo->SetStats(0);
	Histo->SetTitle("");
	if(longhisto==0)
	{
		Histo->GetYaxis()->SetTitleOffset(1.4);
		Histo->GetXaxis()->SetTitleOffset(1.4);
		Histo->GetXaxis()->SetLabelSize(0.05);
		Histo->GetYaxis()->SetLabelSize(0.05);
		Histo->GetXaxis()->SetTitleSize(0.045);
		Histo->GetYaxis()->SetTitleSize(0.045);
		Histo->GetXaxis()->SetNdivisions(505);
		Histo->GetYaxis()->SetNdivisions(505);
	}

	if(longhisto==1)
	{
		Histo->GetYaxis()->SetTitleOffset(0.2);
		Histo->GetXaxis()->SetTitleOffset(1.0);
		//if(big==1)	Histo->GetYaxis()->SetLabelOffset(0.015);
		//if(big==1)	Histo->GetXaxis()->SetLabelOffset(0.015);
		Histo->GetXaxis()->SetLabelSize(0.07);
		Histo->GetYaxis()->SetLabelSize(0.07);
		Histo->GetXaxis()->SetTitleSize(0.08);
		Histo->GetYaxis()->SetTitleSize(0.08);
		//Histo->GetXaxis()->CenterTitle();
		//Histo->GetYaxis()->CenterTitle();
		Histo->GetXaxis()->SetNdivisions(520);
		Histo->GetYaxis()->SetNdivisions(10);
	}

	Histo->GetXaxis()->SetLabelFont(42);
	Histo->GetYaxis()->SetLabelFont(42);
	Histo->GetXaxis()->SetTitleFont(42);
	Histo->GetYaxis()->SetTitleFont(42);
	if(Xtitel!="")Histo->GetXaxis()->SetTitle(Xtitel);
	if(Ytitel!="")Histo->GetYaxis()->SetTitle(Ytitel);

	Histo->SetLineColor(1);
	Histo->SetMarkerColor(1);
	Histo->SetMarkerStyle(20);
	Histo->SetMarkerSize(0.5);
}
//
//
// checks if the cell is part of manually masked cells
//
Bool_t IsCellMaskedByHand(Int_t cell, std::vector<Int_t> cellVector)
{
	Bool_t bad=0;
	for(Int_t i=0; i<(Int_t)cellVector.size();i++)
	{
		if(cell==cellVector.at(i))bad=1;
	}

	return bad;
}
void CreateCellCompPDF(TH2F* hAmpIDMasked, std::vector<Int_t> cellVector,TH1* goodCellsMerged, TH1* goodCellsRbR, TString pdfName)
{
	Int_t NoOfCells=cellVector.size();
	Bool_t firstCanvas=0;
	TString name;
	/*TLatex* textA = new TLatex(0.2,0.8,"*test*");
	textA->SetTextSize(0.08);
	textA->SetTextColor(1);
	textA->SetNDC();
    */
	for(Int_t cell=0;cell<NoOfCells;cell++)
	{
		TString internal_pdfName=pdfName;
		TCanvas *c1;
		if(cell%9==0)
		{
			c1 = new TCanvas(Form("badcells%i",cell),"badcells",1000,750);
			if(cellVector.size() > 6)        c1->Divide(3,3);
			else if (cellVector.size() > 3)  c1->Divide(3,2);
			else                             c1->Divide(3,1);
		}
		TH1 *hCell  = hAmpIDMasked->ProjectionX(Form("Cell %d",cellVector.at(cell)),cellVector.at(cell)+1,cellVector.at(cell)+1);
		TH1 *hCell2 = (TH1*)hCell->Clone("hCell2");

		c1->cd(cell%9 + 1);
		hCell->Divide(goodCellsRbR);
		hCell2->Divide(goodCellsMerged);

		hCell->SetLineColor(kBlue-8);
		hCell2->SetLineColor(kRed-9);
		hCell->GetXaxis()->SetTitle("E (GeV)");
		hCell->GetYaxis()->SetTitle("cell/mean of good");
		hCell->GetXaxis()->SetRangeUser(0.,10.);
		hCell->SetLineWidth(1) ;
		hCell->SetTitle(Form("Cell No. %d",cellVector.at(cell)));
		hCell->Draw("hist");
		hCell2->DrawCopy("same hist");

		//textA->SetTitle(Form("Cell No. %d",cellVector.at(cell)));
		//textA->Draw();

		//..page is full or end of loop
		if(cell%9==8 ||cell == NoOfCells-1)
		{
			if(cell == NoOfCells-1)
			{
				//internal_pdfName +=")";
				c1->Print(Form("%s)",pdfName.Data()));
			}
			else if(firstCanvas==0)
			{
				internal_pdfName +="(";
				c1->Print(internal_pdfName.Data());
				firstCanvas=1;
			}
			else
			{
				c1->Print(internal_pdfName.Data());
			}
			delete c1;
		}
	}

}
