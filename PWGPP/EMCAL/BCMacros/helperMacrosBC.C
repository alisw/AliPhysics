/// \file helperMacrosBC.C
/// \ingroup EMCALOfflineMacros
/// \brief small collection of helper functions
///
///
/// See https://twiki.cern.ch/twiki/bin/view/ALICE/BadChannelAnalysis for general documentation
///
/// ---------------------
///  Running the macro
/// ---------------------
/// use root -b to speed up (no canvas drawn)           <br>
/// root [0] .L helperMacrosBC.C++                      <br>
/// root [0] Get_RowCollumnID(244917,-1,4,16,4)         <br>
/// root [0] Get_RowCollumnID(244917,2001,-1,-1,-1)     <br>
/// root [0] Compare2Blocks("LHC16l",803,3,5)           <br>
///
/// \author Eliane Epple <eliane.epple@yale.edu>, Yale University
///
/// \date June 29, 2017

// --- ROOT system ---
#include <Riostream.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>


// --- ANALYSIS system ---
#include "AliEMCALGeometry.h"          //include when compile
#include "AliCalorimeterUtils.h"       //include when compile
#include "AliAODEvent.h"               //include when compile


//definition of methods
void SetHisto(TH1 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto);
Bool_t IsCellMaskedByHand(Int_t cell, std::vector<Int_t> cellVector);
void CreateCellCompPDF(TH2F* hAmpIDMasked, std::vector<Int_t> cellVector, TH1* goodCellsMerged, TH1* goodCellsRbR, TString pdfName);
void Plot2DCells(TString Block, Int_t runNo, std::vector<Int_t> cellVectorMergeB, std::vector<Int_t> cellVectorMergeA);

/// get all runnumbers from different groups
/// and sort them to give a min and max range for the bad map.
/// Run numbers are in your runList file.
///________________________________________________________________________
void Compare_RunNumbers(TString period="LHC15n",Int_t trainNo=603,TString runList1="",TString runList2="")
{
	//......................................................
	//..PT. II open the 1st text file and save the run IDs into the RunId[] array
	//......................................................
	cout<<"o o o Open .txt file with run indices. Name = " << runList1 << endl;
	TString RunPath1        = Form("./AnalysisInput/%s/Train_%i/%s",period.Data(),trainNo,runList1.Data());
	cout<<"o o o Open .txt file with run indices = " << RunPath1 << endl;
	FILE *pFile1 = fopen(RunPath1.Data(), "r");
	if(!pFile1)
	{
		cout<<"couldn't open file "<<RunPath1<<"!"<<endl;
		return;
	}
	Int_t q1;
	Int_t ncols1;
	Int_t nlines1 = 0 ;
	std::vector<Int_t> RunIdVec1;
	while (1)
	{
		ncols1 = fscanf(pFile1,"  %d ",&q1);
		if (ncols1< 0) break;
		RunIdVec1.push_back(q1);
		nlines1++;
	}
	fclose(pFile1);
	std::sort (RunIdVec1.begin(), RunIdVec1.end());
	//......................................................
	//..PT. II open the 2nd text file and save the run IDs into the RunId[] array
	//......................................................
	cout<<"o o o Open .txt file with run indices. Name = " << runList2 << endl;
	TString RunPath2       = Form("./AnalysisInput/%s/Train_%i/%s",period.Data(),trainNo,runList2.Data());
	cout<<"o o o Open .txt file with run indices = " << RunPath2 << endl;
	FILE *pFile2 = fopen(RunPath2.Data(), "r");
	if(!pFile2)
	{
		cout<<"couldn't open file "<<RunPath2<<"!"<<endl;
		return;
	}
	Int_t q2;
	Int_t ncols2;
	Int_t nlines2 = 0 ;
	std::vector<Int_t> RunIdVec2;
	while (1)
	{
		ncols2 = fscanf(pFile2,"  %d ",&q2);
		if (ncols2< 0) break;
		RunIdVec2.push_back(q2);
		nlines2++;
	}
	fclose(pFile2);
	std::sort (RunIdVec2.begin(), RunIdVec2.end());

	//......................................................
	//..Find different elements in the runlists
	//......................................................
	//find runnumbers present in both lists

	//this part is not yet working
	//what it should do is find unique elements.
	//means I compare list 1 and list 2 and I am interested in elements
	//that are exclusivley in list1 or exclusivley in list2

	std::vector<Int_t> RunIdVecSum;
	cout<<"Elements in Vector Sum a ("<<RunIdVecSum.size()<<")"<<endl;
	RunIdVecSum.insert( RunIdVecSum.end(),RunIdVec1.begin(), RunIdVec1.end());
	cout<<"Elements in Vector Sum b("<<RunIdVecSum.size()<<")"<<endl;
	RunIdVecSum.insert( RunIdVecSum.end(),RunIdVec2.begin(), RunIdVec2.end());
	cout<<"Elements in Vector Sum c("<<RunIdVecSum.size()<<")"<<endl;
	std::sort (RunIdVecSum.begin(), RunIdVecSum.end());

	// using default comparison:
	std::vector<int>::iterator it;
	it = std::unique (RunIdVecSum.begin(), RunIdVecSum.end());
	RunIdVecSum.resize( std::distance(RunIdVecSum.begin(),it) );
	cout<<"Elements in Vector Sum d("<<RunIdVecSum.size()<<")"<<endl;


	cout<<"Elements in Vector 1 ("<<RunIdVec1.size()<<")"<<endl;
	cout<<"Elements in Vector 2 ("<<RunIdVec2.size()<<")"<<endl;

	cout<<"##Unique elements in vector(not really) ("<<RunIdVecSum.size()<<"):"<<endl;
	for(Int_t i=0;i<(Int_t)RunIdVecSum.size();i++)
	{
		cout<<RunIdVecSum.at(i)<<", "<<flush;
	}
	cout<<endl;
}
///
/// check where the bad cell is (row collumn), if you have only its ID
/// or get it's ID if you have its row and collumn
///________________________________________________________________________
void Get_RowCollumnID(Int_t runNum= 244411,Int_t inputCellID=-1,Int_t inputRow=-1,Int_t inputCollumn=-1,Int_t inputSM=-1)
{
	//......................................................
	//..Initialize EMCal/DCal geometry
	AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(runNum);
	fCaloUtils->AccessGeometry(aod);
	AliEMCALGeometry * geom = fCaloUtils->GetEMCALGeometry();

	//..get row collumn from cell ID
	Int_t cellColumn=0,cellRow=0;
	Int_t cellColumnAbs=0,cellRowAbs=0;
	Int_t cellID=0;
	Int_t trash;

	cout<<"...................................................."<<endl;
	cout<<""<<endl;
	if(inputCellID!=-1)
	{
		fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(inputCellID,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
		cout<<"Cell Id provided: "<<inputCellID<<endl;
		cout<<"This corresponds to absolute row: "<<cellRowAbs<<" and absolute collumn: "<<cellColumnAbs<<endl;
	}
	else if(inputRow!=-1 && inputCollumn!=-1 && inputSM!=-1)
	{
		cout<<"Supermodule provided: "<<inputSM<<endl;
		cout<<"Absolute row provided: "<< inputRow<<" and absolute collumn provided : "<<inputCollumn<<endl;
		cellID=geom->GetAbsCellIdFromCellIndexes(inputSM,inputRow,inputCollumn);
		cout<<"This corresponds to Cell Id: "<<cellID<<endl;
	}
	else cout<<"need more information!"<<endl;
	cout<<""<<endl;
	cout<<"...................................................."<<endl;
}

///
/// Compares masked amplidudes from 2 different merged blocks
/// or from 2 different versions of the same block so that one can test
/// the effectivness of added periods
///
//________________________________________________________________________
void Compare2Blocks(TString period="LHC15n",Int_t trainNo=603,Int_t versionA=0, Int_t versionB=1)
{
	gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //..to supress a lot of standard output
	gStyle->SetOptStat(0); //..Do not plot stat boxes
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin(0.06);

    Int_t noOfCells=17674;
    Int_t nBadCellMerged =noOfCells;
    Int_t nBadCellRbR    =noOfCells;
	Int_t array_StartCellSM_Value[21]   ={0,1152,2304,3456,4608,5760,6912,8064,9216,10368,11520,11904,12288,13056,13824,14592,15360,16128,16896,17280,17664};

    //..............................................
    //..manually disable cells
    std::vector<Int_t> badcellsBlock1;
    //badcellsBlock1.insert(badcellsBlock1.end(),{6644,6655,10140,12036,12037,12038,12039,12040,12041,12926,13067,13066,13125});
    //badcellsBlock1.insert(badcellsBlock1.end(),{13133,13483,13971,13978,14116,14118,14122,14411,14593,14599,14600,14606,14699});

    std::vector<Int_t> vOnlyMaskedInMergedA;
	std::vector<Int_t> vOnlyMaskedInMergedB;

	//......................................................
	//..Get the .root file from masking runs in Version1
	//......................................................
	cout<<"** Open file A with merged runlist analysis: "<<endl;
	TString pathA        = Form("./AnalysisOutput/%s/Train_%i/Version%i",period.Data(),trainNo,versionA);
	TString rootFileNameA= Form("%s_INT7_Histograms_V%i.root",period.Data(),versionA);
	TFile* outputRootA   = TFile::Open(Form("%s/%s",pathA.Data(),rootFileNameA.Data()));
	if(!outputRootA)cout<<"File "<<outputRootA->GetName()<<" does not exist"<<endl;
	else cout<<"file A: "<<outputRootA->GetName()<<endl;

	//..get the necessary histograms
	TH2F* hCellAmplitudeA    =(TH2F*)outputRootA->Get("hCellAmplitude");
	TH1F* hnEventsA          =(TH1F*)outputRootA->Get("hNEvents");
	hCellAmplitudeA->Scale(hnEventsA->GetBinContent(1));
	TH2F* hCellAmplitudeBlockMaskA=(TH2F*)hCellAmplitudeA->Clone("hCellAmplitudeMask");
	TH1F* hCellFlagA         =(TH1F*)outputRootA->Get("fhCellFlag");

	//......................................................
	//..Get the .root file from masking runs in Version2
	//......................................................
	cout<<endl;
	cout<<"**Open file B with merged runlist analysis: "<<endl;
	TString pathB        = Form("./AnalysisOutput/%s/Train_%i/Version%i",period.Data(),trainNo,versionB);
	TString rootFileNameB= Form("%s_INT7_Histograms_V%i.root",period.Data(),versionB);
	TFile* outputRootB   = TFile::Open(Form("%s/%s",pathB.Data(),rootFileNameB.Data()));
	if(!outputRootB)cout<<"File "<<outputRootB->GetName()<<" does not exist"<<endl;
	else cout<<"file B: "<<outputRootB->GetName()<<endl;

	//..get the necessary histograms
	TH2F* hCellAmplitudeB    =(TH2F*)outputRootB->Get("hCellAmplitude");
	TH1F* hnEventsB          =(TH1F*)outputRootB->Get("hNEvents");
	hCellAmplitudeB->Scale(hnEventsB->GetBinContent(1));
	TH2F* hCellAmplitudeBlockMaskB=(TH2F*)hCellAmplitudeB->Clone("hCellAmplitudeMask");
	TH1F* hCellFlagB         =(TH1F*)outputRootB->Get("fhCellFlag");

	//......................................................
	//..mask the bad cells according to both versions
	//......................................................
	Int_t maskA=0;
	Int_t maskB=0;

	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		//......................................................
		//..Part A - analyzed as one merged runblock
		//......................................................
		maskA=0;
		maskA=IsCellMaskedByHand(ic,badcellsBlock1);

		//..mask the bad cells
		for (Int_t amp = 1; amp <= hCellAmplitudeBlockMaskA->GetNbinsX(); amp++)
		{
			if(hCellFlagA->GetBinContent(ic+1)>0 || maskA==1)
			{
				hCellAmplitudeBlockMaskA->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellMerged--;
			}
		}

		//......................................................
		//..Part B - analyzed run-by-run, corrected and masked in block
		//......................................................
		maskB=0;
		maskB=IsCellMaskedByHand(ic,badcellsBlock1);

		//..mask the bad cells
		for (Int_t amp = 1; amp <= hCellAmplitudeBlockMaskB->GetNbinsX(); amp++)
		{
			if(hCellFlagB->GetBinContent(ic+1)>0 || maskB==1)
			{
				hCellAmplitudeBlockMaskB->SetBinContent(amp,ic+1,0);
				if(amp==1)nBadCellRbR--;
			}
		}

		//......................................................
		//..Compare the different channels that are marked in the two versions
		//......................................................
		if(!IsCellMaskedByHand(ic,badcellsBlock1))
		{
			if(hCellFlagA->GetBinContent(ic+1)>0  && hCellFlagB->GetBinContent(ic+1)==0)vOnlyMaskedInMergedA.push_back(ic);
			if(hCellFlagA->GetBinContent(ic+1)==0 && hCellFlagB->GetBinContent(ic+1)>0) vOnlyMaskedInMergedB.push_back(ic);
		}
	}
	//..merged runblock
	TH1* projMaskedCellsA = hCellAmplitudeBlockMaskA->ProjectionX("MaskedCellsMergedBlockA");
	TH1* projMaskedCellsB = hCellAmplitudeBlockMaskB->ProjectionX("MaskedCellsMergedBlockB");

	//......................................................
	//..Plot results
	//......................................................
	TCanvas* C1 = new TCanvas("-1-","Projections of A and B",900,900);
	C1->Divide(2,2);
	C1->cd(1)->SetLogy();
  	SetHisto(projMaskedCellsA,"","hits/event",0);
	projMaskedCellsA->DrawCopy("hist");
	projMaskedCellsB->SetLineColor(8);
	projMaskedCellsB->DrawCopy("same hist");
	C1->cd(2);
	projMaskedCellsA->Divide(projMaskedCellsB);
  	SetHisto(projMaskedCellsA,"","Merged Version A/Merged Version B",0);
	projMaskedCellsA->DrawCopy("hist");
	projMaskedCellsA->Multiply(projMaskedCellsB);

	TLatex* textSM = new TLatex(0.1,0.1,"*test*");
	textSM->SetTextSize(0.06);
	textSM->SetTextColor(1);
	textSM->SetNDC();

	TCanvas *c1_projSM = new TCanvas("-1a-","Projections of A and B per SM",1200,900);
	c1_projSM->Divide(5,4,0.001,0.001);
	TH1* projEnergyMaskASM[20];
	TH1* projEnergyMaskBSM[20];
	for(Int_t iSM=0;iSM<20;iSM++)
	{
		c1_projSM->cd(iSM+1)->SetLogy();
		gPad->SetTopMargin(0.03);
		gPad->SetBottomMargin(0.11);
		projEnergyMaskASM[iSM] = hCellAmplitudeBlockMaskA->ProjectionX(Form("%sMask_ProjSM%i",hCellAmplitudeBlockMaskA->GetName(),iSM),array_StartCellSM_Value[iSM],array_StartCellSM_Value[iSM+1]-1);
		projEnergyMaskASM[iSM]->SetTitle("");
		projEnergyMaskASM[iSM]->SetXTitle(Form("Cell Energy [GeV], SM%i",iSM));
		projEnergyMaskASM[iSM]->GetYaxis()->SetTitleOffset(1.6);
		projEnergyMaskASM[iSM]->GetYaxis()->SetLabelSize(0.06);
		projEnergyMaskASM[iSM]->GetXaxis()->SetLabelSize(0.06);
		projEnergyMaskASM[iSM]->GetXaxis()->SetRangeUser(0,20);
		projEnergyMaskASM[iSM]->GetXaxis()->SetTitleSize(0.06);
		projEnergyMaskASM[iSM]->SetLineColor(kGreen+1);
		projEnergyMaskASM[iSM]->DrawCopy(" hist");

		projEnergyMaskBSM[iSM] = hCellAmplitudeBlockMaskB->ProjectionX(Form("%s_ProjSM%i",hCellAmplitudeBlockMaskB->GetName(),iSM),array_StartCellSM_Value[iSM],array_StartCellSM_Value[iSM+1]-1);
		projEnergyMaskBSM[iSM]->DrawCopy("same hist");

		textSM->SetTitle(Form("Includes cell IDs %d-%d",array_StartCellSM_Value[iSM],array_StartCellSM_Value[iSM+1]-1));
		textSM->DrawLatex(0.2,0.8,Form("Includes cell IDs %d-%d",array_StartCellSM_Value[iSM],array_StartCellSM_Value[iSM+1]-1));
	}
	TCanvas *c1_projSMR = new TCanvas("-1b-","Ratio of A/B per SM",1200,900);
	c1_projSMR->Divide(5,4,0.001,0.001);
	TH1* projEnergyMaskSM[20];
	TH1* projEnergySM[20];
	for(Int_t iSM=0;iSM<20;iSM++)
	{
		c1_projSMR->cd(iSM+1);
		gPad->SetTopMargin(0.03);
		gPad->SetBottomMargin(0.11);
		projEnergyMaskASM[iSM]->Divide(projEnergyMaskBSM[iSM]);
		projEnergyMaskASM[iSM]->SetLineColor(kGray+1);
		projEnergyMaskASM[iSM]->DrawCopy(" hist");

		textSM->SetTitle(Form("Includes cell IDs %d-%d",array_StartCellSM_Value[iSM],array_StartCellSM_Value[iSM+1]-1));
		textSM->DrawLatex(0.2,0.8,Form("Includes cell IDs %d-%d",array_StartCellSM_Value[iSM],array_StartCellSM_Value[iSM+1]-1));
	}

	TCanvas* C3 = new TCanvas("-3-","2D of A and B",900,900);
	C3->Divide(2,2);
	C3->cd(1)->SetLogz();
  	SetHisto(hCellAmplitudeBlockMaskA,"","cell ID",0);
  	hCellAmplitudeBlockMaskA->DrawCopy("colz");
	C3->cd(2)->SetLogz();
  	SetHisto(hCellAmplitudeBlockMaskB,"","cell ID",0);
  	hCellAmplitudeBlockMaskB->DrawCopy("colz");


 	//......................................................
	//..Print out compared cells and plot the spectra
	//......................................................
	projMaskedCellsA->Scale(1.0/nBadCellMerged);
	projMaskedCellsB->Scale(1.0/nBadCellRbR);


	cout<<"  Cells masked in version A and not in version B ("<<vOnlyMaskedInMergedA.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInMergedA.size();i++)
	{
		cout<<vOnlyMaskedInMergedA.at(i)<<","<<flush;
	}
	cout<<endl;
	TString outNameA = Form("./AnalysisOutput/%s/Train_%i/cOnlyMergedBlockA.pdf",period.Data(),trainNo);
	CreateCellCompPDF(hCellAmplitudeBlockMaskB,vOnlyMaskedInMergedA,projMaskedCellsA,projMaskedCellsB,outNameA);
	cout<<"  Cells masked in version B and not in version A ("<<vOnlyMaskedInMergedB.size()<<"):"<<endl;
	for(Int_t i=0; i<(Int_t)vOnlyMaskedInMergedB.size();i++)
	{
		cout<<vOnlyMaskedInMergedB.at(i)<<","<<flush;
	}
	cout<<endl;
	TString outNameB = Form("./AnalysisOutput/%s/Train_%i/cOnlyMergedBlockB.pdf",period.Data(),trainNo);
	CreateCellCompPDF(hCellAmplitudeBlockMaskA,vOnlyMaskedInMergedB,projMaskedCellsA,projMaskedCellsB,outNameB);

 	//......................................................
	//..build two dimensional histogram with cells rejected from
	//..the one or the other method
	//......................................................
	Plot2DCells("-",244917,vOnlyMaskedInMergedB,vOnlyMaskedInMergedA);
}
//-----------------------------------------------------------
// All functions below this point are helping the
// Compare2Blocks function
//-----------------------------------------------------------

///
///
///
void Plot2DCells(TString Block, Int_t runNo, std::vector<Int_t> cellVectorMergeB, std::vector<Int_t> cellVectorMergeA)
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

	TH2F *plot2D_VersionB   = new TH2F(Form("Block%s_VersionB",Block.Data()),Form("Block%s_VersionB",Block.Data()),fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	plot2D_VersionB->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_VersionB->GetYaxis()->SetTitle("cell row (#phi direction)");
	TH2F *plot2D_VersionA = new TH2F(Form("Block%s_VersionA",Block.Data()),Form("Block%s_VersionA",Block.Data()),fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	plot2D_VersionA->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_VersionA->GetYaxis()->SetTitle("cell row (#phi direction)");

	Int_t cellColumn=0,cellRow=0;
	Int_t cellColumnAbs=0,cellRowAbs=0;
	Int_t trash;

	for(Int_t i = 0; i < (Int_t)cellVectorMergeB.size(); i++)
	{
		Int_t cell=cellVectorMergeB.at(i);
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
		plot2D_VersionB->Fill(cellColumnAbs,cellRowAbs);
	}
	for(Int_t i = 0; i < (Int_t)cellVectorMergeA.size(); i++)
	{
		Int_t cell=cellVectorMergeA.at(i);
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
		plot2D_VersionA->Fill(cellColumnAbs,cellRowAbs,1);
	}
	//. . . . . . . . . . . . . . . . . . . .
	TCanvas *c1 = new TCanvas(Form("2DMapForBlock%s",Block.Data()),Form("2DMapForBlock%s",Block.Data()),900,500);
	c1->ToggleEventStatus();
	c1->Divide(2);
	c1->cd(1);
	plot2D_VersionB->Draw("colz");
	c1->cd(2);
	plot2D_VersionA->Draw("colz");
}
///
/// checks if the cell is part of manually masked cells
///
Bool_t IsCellMaskedByHand(Int_t cell, std::vector<Int_t> cellVector)
{
	Bool_t bad=0;
	for(Int_t i=0; i<(Int_t)cellVector.size();i++)
	{
		if(cell==cellVector.at(i))bad=1;
	}

	return bad;
}
///
///
///
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

		c1->cd(cell%9 + 1)->SetLogy();
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
///
/// Funtion to set TH1 histograms to a similar style
///
///________________________________________________________________________
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
