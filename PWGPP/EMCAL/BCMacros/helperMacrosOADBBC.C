/// \file helperMacrosOADBBC.C
/// \ingroup EMCALOfflineMacros
/// \brief Macro to test the files on OADB and provide a run range to be committed to OADB
///
///
/// See https://twiki.cern.ch/twiki/bin/view/ALICE/BadChannelAnalysis for general documentation
///
/// ---------------------
///  Running the macro
/// ---------------------
/// use root -b to speed up (no canvas drawn)                   <br>
/// root [0] .L helperMacrosOADBBC.C++                          <br>
/// root [1] Sort_RunNumbers("LHC16o",663,"runList.txt")        <br>
/// root [2] Test_OADB("LHC16k",804,0,"runList16k.txt")         <br>
/// root [3] Plot_CellList("LHC15o",771,"1","List115o.txt")         <br>
///
/// \author Eliane Epple <eliane.epple@yale.edu>, Yale University
/// \date June 29, 2017



// --- ROOT system ---
#include <Riostream.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>

// --- ANALYSIS system ---
#include "AliEMCALGeometry.h"          //include when compile
#include "AliCalorimeterUtils.h"       //include when compile
#include "AliAODEvent.h"               //include when compile
#include "AliOADBContainer.h"          //include when compile


/// get all runnumbers from different groups
/// and sort them to give a min and max range for the bad map.
/// Run numbers are in your runList file.
///________________________________________________________________________
void Sort_RunNumbers(TString period="LHC15n",Int_t trainNo=603,TString runList="")
{
	//......................................................
	//..open the text file and save the run IDs into the RunId[] array
	cout<<"o o o Open .txt file with run indices. Name = " << runList << endl;
	TString RunPath1        = Form("./AnalysisInput/%s/Train_%i/%s",period.Data(),trainNo,runList.Data());
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
	//.......................................................
	cout<<"Runs in Order: "<<endl;
	for(Int_t iRun = 0; iRun < (Int_t)RunIdVec1.size(); iRun++)
	{
		cout<<RunIdVec1.at(iRun)<<", "<<flush;
	}
	cout<<endl;
}
///
/// Test if the file committed to OADB is correct
/// look at your local .root file with the bad maps inside and see what
/// the bad map looks at a certain runNumber. If everything is committed
/// correctly they should show the same
///________________________________________________________________________
void Test_OADB(TString period="LHC15n",Int_t trainNo=603,TString version="INT7Emc",TString runList="")
{
    gStyle->SetOptStat(0); //..Do not plot stat boxes
	//......................................................
	// Test if the file committed to OADB is correct
	//......................................................
	Int_t nSM = 20;
	TH1C *h[20];
	Int_t cellColumnAbs,cellRowAbs,trash,cellID;
	TString summaryPDF=Form("./AnalysisOutput/%s/Train_%i/OADB_Summary.pdf",period.Data(),trainNo);

	//......................................................
	//..open the text file and save the run IDs into the RunId[] array
	cout<<"o o o Open .txt file with run indices. Name = " << runList << endl;
	TString RunPath        = Form("./AnalysisInput/%s/Train_%i/%s",period.Data(),trainNo,runList.Data());
	cout<<"o o o Open .txt file with run indices = " << RunPath << endl;
	FILE *pFile = fopen(RunPath.Data(), "r");
	if(!pFile)
	{
		cout<<"couldn't open file "<<RunPath<<"!"<<endl;
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
    Int_t nRuns=RunIdVec.size();
	//......................................................
	//..Get the OADB information
    //TString fBasePath="/Users/Eliane/Software/alice/sw/osx_x86-64/AliPhysics/latest-ali-master/OADB/EMCAL"; //..from AliPhysics
    TString fBasePath="/Users/Eliane/Software/BadChannelAnalysis"; //..locally from OADB commit/e-mail to test quickly

    AliOADBContainer *cont=new AliOADBContainer("");
	cont->InitFromFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"AliEMCALBadChannels");
	//......................................................
	//..Get the .root file with the original histogram to compare if they coincide
	//TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Train_%i/VersionINT7Glob",period.Data(),trainNo);
	//TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Train_%i/Version4ManMasked",period.Data(),trainNo);
	TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Train_%i/Version%s",period.Data(),trainNo,version.Data());
	TString rootFileName= Form("%s_INT7_Histograms_V%s.root",period.Data(),version.Data());
	TFile* outputRoot   = TFile::Open(Form("%s/%s",path.Data(),rootFileName.Data()));

	if(!outputRoot)cout<<"File "<<outputRoot->GetName()<<" does not exist"<<endl;
	TH2F* h2DChannelMap_FlagBad =(TH2F*)outputRoot->Get("2DChannelMap_Flag2");
	TH2F* h2DChannelMap_FlagDead=(TH2F*)outputRoot->Get("2DChannelMap_Flag1");

	//......................................................
	//..Initialize EMCal/DCal geometry
	AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(RunIdVec.at(0));
	fCaloUtils->AccessGeometry(aod);
	AliEMCALGeometry * geom = fCaloUtils->GetEMCALGeometry();

	//.......................................................
	//..build two dimensional histogram with values row vs. column
	//..with info from OADB
	Int_t fNMaxCols    = 48;  //eta direction
	Int_t fNMaxRows    = 24;  //phi direction

	Int_t fNMaxColsAbs = 2*fNMaxCols;
	Int_t fNMaxRowsAbs = Int_t (nSM/2)*fNMaxRows; //multiply by number of supermodules
	TString histoName;
	histoName = Form("2DChannelMap_Flag_Bad");
	TH2F *plot2D_Bad_OADB = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	//*OADB  looks like Marcels figures*/	TH2F *plot2D_Bad_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
	plot2D_Bad_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_Bad_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");
	histoName = Form("2DChannelMap_Flag_Dead");
	TH2F *plot2D_Dead_OADB = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	//*OADB looks like Marcels figures*/	TH2F *plot2D_Dead_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
	plot2D_Dead_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_Dead_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");

	//.......................................................
	//.. Read the Bad Channel map from OADB
	TCanvas* C2 = new TCanvas("CompareCanvas","Compare OADB to LocalFile",900,900);
	C2->Divide(2,2);
	TCanvas* C4 = new TCanvas("Subtraction of OADB-Orig.","Subtraction of OADB-Orig.",900,900);
	C4->Divide(2,2);
	TLatex* textA = new TLatex(0.5,0.8,"If empty -> good!");
	textA->SetTextSize(0.04);
	textA->SetTextColor(1);
	textA->SetNDC();

	cout<<"Checking "<<nRuns<<" runs: "<<endl;
	std::vector<Int_t> RunsWithoutMap;
	std::vector<Int_t> RunsWithMap;

	for(Int_t iRun = 0; iRun < nRuns; iRun++)
	{
		//cout<<"------ run "<<RunIdVec.at(iRun)<<endl;
		if(iRun%5==0)cout<<"."<<flush;
		if(iRun%20==0)cout<<"Run No."<<iRun<<endl;
		TObjArray *recal=(TObjArray*)cont->GetObject(RunIdVec.at(iRun));
		if(!recal)
		{
			cout<<"Error - No bad map for run Number "<<RunIdVec.at(iRun)<<" online!!"<<endl;
			RunsWithoutMap.push_back(RunIdVec.at(iRun));
			continue;
		}
		RunsWithMap.push_back(RunIdVec.at(iRun));

		plot2D_Bad_OADB ->Reset();
		plot2D_Dead_OADB->Reset();
		for(Int_t iSM = 0; iSM < nSM; iSM++)
		{
			//..Bad map fore each module
			h[iSM]=(TH1C *)recal->FindObject(Form("EMCALBadChannelMap_Mod%d",iSM));
			//..Loop though the SM to set which cells are bad
			for(Int_t column=0;column<48;column++)
			{
				for(Int_t row=0;row<24;row++)
				{
					Int_t inRow=row;
					Int_t inCol=column;
					cellID=geom->GetAbsCellIdFromCellIndexes(iSM,inRow,inCol);
					fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(cellID,0,inCol,inRow,trash,cellColumnAbs,cellRowAbs);
					if(h[iSM]->GetBinContent(column,row)>1)//..bad and warm
					{
						plot2D_Bad_OADB->SetBinContent(cellColumnAbs,cellRowAbs,1);
					}
					if(h[iSM]->GetBinContent(column,row)==1)//..dead
					{
						plot2D_Dead_OADB->SetBinContent(cellColumnAbs,cellRowAbs,1);
					}
				}
			}
		}
		//..................................................................
		TLatex* text=new TLatex(0.15,0.85,Form("Run Number %i",RunIdVec.at(iRun)));
		text->SetTextSize(0.05);
		text->SetNDC();
		text->SetTextColor(1);
		text->SetTextAngle(0);

		C2->cd(1);
		plot2D_Bad_OADB->SetTitle("Bad Cells OADB");
		plot2D_Bad_OADB->DrawCopy("colz");
		text->Draw();
		C2->cd(2);
		plot2D_Dead_OADB->SetTitle("Dead Cells OADB");
		plot2D_Dead_OADB->DrawCopy("colz");
		C2->cd(3);
		h2DChannelMap_FlagBad->SetTitle("Bad Cells LocalFile");
		h2DChannelMap_FlagBad->DrawCopy("colz");
		C2->cd(4);
		h2DChannelMap_FlagDead->SetTitle("Dead Cells LocalFile");
		h2DChannelMap_FlagDead->DrawCopy("colz");

		//..................................................................
		C4->cd(1);
		plot2D_Bad_OADB->SetTitle("OADB-Local File (Bad)");
		plot2D_Bad_OADB->Add(h2DChannelMap_FlagBad,-1);
		plot2D_Bad_OADB->DrawCopy("colz");
		plot2D_Bad_OADB->Add(h2DChannelMap_FlagBad,+1);
		textA->Draw();
		text->Draw();
		C4->cd(2);
		plot2D_Dead_OADB->SetTitle("OADB-Local File (Dead)");
		plot2D_Dead_OADB->Add(h2DChannelMap_FlagDead,-1);
		plot2D_Dead_OADB->DrawCopy("colz");
		plot2D_Dead_OADB->Add(h2DChannelMap_FlagDead,+1);
		textA->Draw();
		C4->cd(3);
		plot2D_Bad_OADB->SetTitle("OADB/Local File (Bad)");
		plot2D_Bad_OADB->Divide(h2DChannelMap_FlagBad);
		plot2D_Bad_OADB->DrawCopy("colz");
		C4->cd(4);
		plot2D_Dead_OADB->SetTitle("OADB/Local File (Dead)");
		plot2D_Dead_OADB->Divide(h2DChannelMap_FlagDead);
		plot2D_Dead_OADB->DrawCopy("colz");

		//..................................................................
        //..Save to PDF
		//..Add figures to the summary canvas
		if(iRun==0)C2   ->Print(Form("%s(",summaryPDF.Data()));
		else       C2   ->Print(Form("%s(",summaryPDF.Data()));
		if(iRun==nRuns-1)C4   ->Print(Form("%s)",summaryPDF.Data()));
		else             C4   ->Print(Form("%s",summaryPDF.Data()));
	}//end of run loop
	cout<<"==Total Summary=="<<endl;
	//..print runs without a correct map.
	cout<<"Runs with a bad map ("<<RunsWithMap.size()<<"): "<<flush;
	for(Int_t i=0;i<(Int_t)RunsWithMap.size();i++)
	{
		cout<<RunsWithMap.at(i)<<","<<flush;
	}
	cout<<endl;
	cout<<"Runs without a bad map ("<<RunsWithoutMap.size()<<"): "<<flush;
	for(Int_t i=0;i<(Int_t)RunsWithoutMap.size();i++)
	{
		cout<<RunsWithoutMap.at(i)<<","<<flush;
	}
	cout<<endl;
}
///
/// After OADB maps are online you will get
/// suggestions from other users that think certain
/// channels should be masked. You can look at them here
///________________________________________________________________________
void Plot_CellList(TString period="LHC15n",Int_t trainNo=603,TString version="5",TString cellList="")
{
	//......................................................
	//..open the text file and save the run IDs into the RunId[] array
	cout<<"o o o Open .txt file with run indices. Name = " << cellList << endl;
	FILE *pFile = fopen(cellList.Data(), "r");
	if(!pFile)
	{
		cout<<"couldn't open file "<<cellList<<"!"<<endl;
		return;
	}
	Int_t q;
	Int_t ncols;
	Int_t nlines = 0 ;
	Int_t cellId[500] ;
	std::vector<Int_t> cellIdVec;
	while (1)
	{
		ncols = fscanf(pFile,"  %d ",&q);
		if (ncols< 0) break;
		cellId[nlines]=q;
		cellIdVec.push_back(q);
		nlines++;
	}
	fclose(pFile);
	//..sort the vector by size to be shure to use the right order
	std::sort (cellIdVec.begin(), cellIdVec.end());
    Int_t nCells=cellIdVec.size();

	//......................................................
	//..Get the .root file with the original histogram to compare if they coincide
	TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Train_%i/Version%s",period.Data(),trainNo,version.Data());
	TString rootFileName= Form("%s_INT7_Histograms_V%s.root",period.Data(),version.Data());
	TFile* outputRoot   = TFile::Open(Form("%s/%s",path.Data(),rootFileName.Data()));
	if(!outputRoot)cout<<"File "<<outputRoot->GetName()<<" does not exist"<<endl;
	TH2F* h2DAmp   =(TH2F*)outputRoot->Get("hCellAmplitude");
	TH2F* h2DRatio =(TH2F*)outputRoot->Get("ratio2DAmp");

	TCanvas *c1 = new TCanvas(1);
	c1->Divide(2);
	c1->cd(1);
	h2DAmp->Draw("colz");
	c1->cd(2);
	h2DRatio->Draw("colz");



	Int_t totalperCv = 16;
	Int_t nPad = TMath::Sqrt(totalperCv);
	Int_t nCv  = nCells/totalperCv+1;
	if(nCv<1)nCv=1;

	cout<<"    o create: "<<nCv<<" Canvases with "<<nPad*nPad<<" pads"<<endl;
	//..to compare specific cells over the runs
	TCanvas **cComp     = new TCanvas*[nCv];
	TCanvas **cCompR     = new TCanvas*[nCv];
	for(Int_t i=0;i<nCv;i++)
	{
		cComp[i]    = new TCanvas(TString::Format("CompareGood%d", i), TString::Format("V) Candidates (%d/%d)", i+1, nCv), 1000,750);
		cComp[i]    ->Divide(nPad,nPad,0.001,0.001);
		cCompR[i]    = new TCanvas(TString::Format("CompareGoodR%d", i), TString::Format("V) Ratio (%d/%d)", i+1, nCv), 1000,750);
		cCompR[i]    ->Divide(nPad,nPad,0.001,0.001);
	}
	cout<<"    o Fill Canvases with bad cells canvases"<<endl;

	for(Int_t icell = 0; icell < nCells; icell++)
	{
		TH1D *htmpCellAllRuns     =h2DAmp->ProjectionX(TString::Format("hIDProj_cell%d", icell), icell+1, icell+1);
		cout<<"cell number:"<<cellIdVec.at(icell)<<endl;
		cComp[(icell)/totalperCv]->cd(((icell)%totalperCv)+1)->SetLogy();
		//SetHisto(htmpCellAllRuns,Form("Energy of cell %i",cellIdVec.at(icell)),"Number of Hits",0);
		//htmpCellAllRuns->GetXaxis()->SetRangeUser(0,3);
		htmpCellAllRuns->Draw("hist");


		TH1D *htmpCellRatioAllRuns     =h2DRatio->ProjectionX(TString::Format("hIDRatioProj_cell%d", icell), icell+1, icell+1);
		cout<<"cell number:"<<cellIdVec.at(icell)<<endl;
		cCompR[(icell)/totalperCv]->cd(((icell)%totalperCv)+1)->SetLogy();
		//SetHisto(htmpCellAllRuns,Form("Energy of cell %i",cellIdVec.at(icell)),"Number of Hits",0);
		//htmpCellAllRuns->GetXaxis()->SetRangeUser(0,3);
		htmpCellRatioAllRuns->Draw("hist");
}

	TString pdfName="MoreBadCellsCandidates.pdf";
	//..plot the canvases of cells into a .pdf file
	for(Int_t can=0;can<nCv;can++)
	{
		//TString internal_pdfName1=pdfName+"Low.pdf";
		if(can==0)
		{
			//..first pad
			//internal_pdfName1 +="(";
			cComp[can]    ->Print(Form("%s(",pdfName.Data()));
		}
		else if(can==(nCv-1))
		{
			//..last pad
			//internal_pdfName1 +=")";
			cComp[can]    ->Print(Form("%s)",pdfName.Data()));
		}
		else
		{
			//..all pads in between
			cComp[can]    ->Print(Form("%s",pdfName.Data()));
		}
	}



}


