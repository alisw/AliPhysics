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
#include <TLegend.h>
#include <TStyle.h>

// --- ANALYSIS system ---
#include "AliEMCALGeometry.h"          //include when compile
#include "AliCalorimeterUtils.h"       //include when compile
#include "AliAODEvent.h"               //include when compile
#include "AliOADBContainer.h"          //include when compile

void SetHisto(TH1 *Histo,TString Xtitel,TString Ytitel);
void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2, Float_t lMargin = 0.15, Float_t rMargin = 0.05, Float_t bMargin = 0.15, Float_t tMargin = 0.05);

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
    TString fBasePath="."; //..locally from OADB commit/e-mail to test quickly

    AliOADBContainer *cont=new AliOADBContainer("");
	cont->InitFromFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"AliEMCALBadChannels");
	//......................................................
	//..Get the .root file with the original histogram to compare if they coincide
	//TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Train_%i/VersionINT7Glob",period.Data(),trainNo);
	//TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Train_%i/Version4ManMasked",period.Data(),trainNo);
	TString path        = Form("./AnalysisOutput/%s/Train_%i/Version%s",period.Data(),trainNo,version.Data());
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
	TCanvas* C2 = new TCanvas("CompareCanvas","Compare OADB to LocalFile",1400,800);
	C2->Divide(4,2);
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
		plot2D_Dead_OADB->SetTitle("Dead Cells OADB");
		plot2D_Dead_OADB->DrawCopy("colz");
		text->Draw();
		C2->cd(5);
		plot2D_Bad_OADB->SetTitle("Bad Cells OADB");
		plot2D_Bad_OADB->DrawCopy("colz");
		C2->cd(2);
		h2DChannelMap_FlagDead->SetTitle("Dead Cells LocalFile");
		h2DChannelMap_FlagDead->DrawCopy("colz");
		C2->cd(6);
		h2DChannelMap_FlagBad->SetTitle("Bad Cells LocalFile");
		h2DChannelMap_FlagBad->DrawCopy("colz");

		//..................................................................
		C2->cd(4);
		plot2D_Dead_OADB->SetTitle("OADB-Local File (Dead)");
		plot2D_Dead_OADB->Add(h2DChannelMap_FlagDead,-1);
		plot2D_Dead_OADB->DrawCopy("colz");
		plot2D_Dead_OADB->Add(h2DChannelMap_FlagDead,+1);
		textA->Draw();
		C2->cd(8);
		plot2D_Bad_OADB->SetTitle("OADB-Local File (Bad)");
		plot2D_Bad_OADB->Add(h2DChannelMap_FlagBad,-1);
		plot2D_Bad_OADB->DrawCopy("colz");
		plot2D_Bad_OADB->Add(h2DChannelMap_FlagBad,+1);
		textA->Draw();
		C2->cd(3);
		plot2D_Dead_OADB->SetTitle("OADB/Local File (Dead)");
		plot2D_Dead_OADB->Divide(h2DChannelMap_FlagDead);
		plot2D_Dead_OADB->DrawCopy("colz");
		C2->cd(7);
		plot2D_Bad_OADB->SetTitle("OADB/Local File (Bad)");
		plot2D_Bad_OADB->Divide(h2DChannelMap_FlagBad);
		plot2D_Bad_OADB->DrawCopy("colz");

		//..................................................................
        //..Save to PDF
		//..Add figures to the summary canvas
		if(iRun==0)             C2   ->Print(Form("%s(",summaryPDF.Data()));
		else if (iRun==nRuns-1) C2   ->Print(Form("%s)",summaryPDF.Data()));
		else                    C2   ->Print(Form("%s",summaryPDF.Data()));
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
	gStyle->SetPadTopMargin(0.05);//0.05
	gStyle->SetPadBottomMargin(0.18);//0.15
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadLeftMargin(0.17);
	gStyle->SetFrameFillColor(10);
	//......................................................
	//..open the text file and save the run IDs into the RunId[] array
	cout<<"o o o Open .txt suggested cell IDs. Name = " << cellList << endl;
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
	TString path        = Form("./AnalysisOutput/%s/Train_%i",period.Data(),trainNo);
	//TString rootFileName= Form("Version%s/%s_INT7_Histograms_V%s.root",period.Data(),version.Data());
	TString rootFileName1= Form("Version1ManMasked/%s_INT7_Histograms_V1.root",period.Data());
	TString rootFileName2= Form("Version2ManMasked/%s_INT7_Histograms_V2.root",period.Data());
	TString rootFileName3= Form("Version3ManMasked/%s_INT7_Histograms_V3.root",period.Data());
	TString rootFileName4= Form("Version4ManMasked/%s_INT7_Histograms_V4.root",period.Data());
	TFile* outputRoot1   = TFile::Open(Form("%s/%s",path.Data(),rootFileName1.Data()));
	if(!outputRoot1)cout<<"File "<<outputRoot1->GetName()<<" does not exist"<<endl;
	TH2F* h2DAmp1   =(TH2F*)outputRoot1->Get("hCellAmplitude");
	TH2F* h2DRatio1 =(TH2F*)outputRoot1->Get("ratio2DAmp");
	TFile* outputRoot2   = TFile::Open(Form("%s/%s",path.Data(),rootFileName2.Data()));
	if(!outputRoot2)cout<<"File "<<outputRoot2->GetName()<<" does not exist"<<endl;
	TH2F* h2DAmp2   =(TH2F*)outputRoot2->Get("hCellAmplitude");
	TH2F* h2DRatio2 =(TH2F*)outputRoot2->Get("ratio2DAmp");
	TFile* outputRoot3   = TFile::Open(Form("%s/%s",path.Data(),rootFileName3.Data()));
	if(!outputRoot3)cout<<"File "<<outputRoot3->GetName()<<" does not exist"<<endl;
	TH2F* h2DAmp3   =(TH2F*)outputRoot3->Get("hCellAmplitude");
	TH2F* h2DRatio3 =(TH2F*)outputRoot3->Get("ratio2DAmp");
	TFile* outputRoot4   = TFile::Open(Form("%s/%s",path.Data(),rootFileName4.Data()));
	if(!outputRoot3)cout<<"File "<<outputRoot4->GetName()<<" does not exist"<<endl;
	TH2F* h2DAmp4   =(TH2F*)outputRoot4->Get("hCellAmplitude");
	TH2F* h2DRatio4 =(TH2F*)outputRoot4->Get("ratio2DAmp");

	TCanvas *c1 = new TCanvas(1);
	c1->Divide(2);
	c1->cd(1);
	SetHisto(h2DAmp1,"","");
	h2DAmp1->Draw("colz");
	c1->cd(2);
	SetHisto(h2DRatio1,"","");
	h2DRatio1->Draw("colz");

	//.. be aware of the special sturcture of the canvas
	//.. canvas has totalperCv*2 pads
//	Int_t totalperCv = 4;
	Int_t totalperCv = 8;
	Int_t nPad = TMath::Sqrt(totalperCv);
	Int_t nCv  = nCells/totalperCv+1;
	if(nCv<1)nCv=1;

	cout<<"    o create: "<<nCv<<" Canvases with "<<nPad*nPad<<" pads"<<endl;
	//..to compare specific cells over the runs
	TCanvas **cCompAll  = new TCanvas*[nCv];
	for(Int_t i=0;i<nCv;i++)
	{
		cCompAll[i] = new TCanvas(TString::Format("CompareGoodAll%d", i), TString::Format("V) Both (%d/%d)", i+1, nCv), 1000,750);
		CanvasPartition(cCompAll[i],4,totalperCv/2,0.15,0.02,0.13,0.05);
	}
	TLegend *leg2 = new TLegend(0.60,0.60,0.9,0.85);
	cout<<"    o Fill Canvases with bad cells histograms"<<endl;

	Int_t zoomRange=20;
	Double_t range=10;
	Int_t shift=0;
	cout<<"cell number: "<<endl;
	for(Int_t icell = 0; icell < nCells; icell++)
	{
		Int_t cellID=cellIdVec.at(icell);
		cout<<cellID<<", "<<flush;

		TH1D *htmpCellAllRuns1     =h2DAmp1->ProjectionX(TString::Format("hIDProj1_cell%d", cellID), cellID+1, cellID+1);
		SetHisto(htmpCellAllRuns1,Form("Energy [cell %i]",cellID),"Number of Hits/Events");

		TH1D *htmpCellRatioAllRuns1=h2DRatio1->ProjectionX(TString::Format("hIDRProj1_cell%d", cellID), cellID+1, cellID+1);
		SetHisto(htmpCellRatioAllRuns1,Form("Energy [cell %i]",cellID),"No. Hits/av. No. Hits");

		TH1D *htmpCellAllRuns2     =h2DAmp2->ProjectionX(TString::Format("hIDProj2_cell%d", cellID), cellID+1, cellID+1);
		TH1D *htmpCellRatioAllRuns2=h2DRatio2->ProjectionX(TString::Format("hIDRProj2_cell%d", cellID), cellID+1, cellID+1);
		TH1D *htmpCellAllRuns3     =h2DAmp3->ProjectionX(TString::Format("hIDProj3_cell%d", cellID), cellID+1, cellID+1);
		TH1D *htmpCellRatioAllRuns3=h2DRatio3->ProjectionX(TString::Format("hIDRProj3_cell%d", cellID), cellID+1, cellID+1);
		TH1D *htmpCellAllRuns4     =h2DAmp4->ProjectionX(TString::Format("hIDProj4_cell%d", cellID), cellID+1, cellID+1);
		TH1D *htmpCellRatioAllRuns4=h2DRatio4->ProjectionX(TString::Format("hIDRProj4_cell%d", cellID), cellID+1, cellID+1);

		if(((icell)%8)<4) shift=0;
		else              shift=8;
		cCompAll[(icell)/totalperCv]->cd(((icell)%4)+1+shift)->SetLogy();
		htmpCellAllRuns1->GetXaxis()->SetRangeUser(0,range);
		htmpCellAllRuns1->Draw("hist");
		htmpCellAllRuns2->SetLineColor(kBlue-7);
		htmpCellAllRuns2->DrawCopy("same hist");
		htmpCellAllRuns3->SetLineColor(kGreen-2);
		htmpCellAllRuns3->DrawCopy("same hist");
		htmpCellAllRuns4->SetLineColor(kViolet-1);
		htmpCellAllRuns4->DrawCopy("same hist");

		if(icell==0)
		{
			leg2->AddEntry(htmpCellAllRuns1,"Block1","l");
			leg2->AddEntry(htmpCellAllRuns2,"Block2","l");
			leg2->AddEntry(htmpCellAllRuns3,"Block3","l");
			leg2->AddEntry(htmpCellAllRuns4,"Block4","l");
			leg2->SetTextSize(0.07);
			leg2->SetBorderSize(0);
			leg2->SetFillColorAlpha(10, 0);
			leg2->Draw("same");
		}
		else if((icell)%4==0)
		{
			leg2->Draw("same");
		}

		cCompAll[(icell)/totalperCv]->cd(((icell)%4)+5+shift)->SetLogy();
		htmpCellRatioAllRuns1->GetXaxis()->SetRangeUser(0,range);
		//htmpCellRatioAllRuns1->GetYaxis()->SetRangeUser(0,5);
		htmpCellRatioAllRuns1->Draw("hist");
		htmpCellRatioAllRuns2->SetLineColor(kBlue-7);
		htmpCellRatioAllRuns2->Draw("same hist");
		htmpCellRatioAllRuns3->SetLineColor(kGreen-2);
		htmpCellRatioAllRuns3->Draw("same hist");
		htmpCellRatioAllRuns4->SetLineColor(kViolet-1);
		htmpCellRatioAllRuns4->Draw("same hist");

	}
	cout<<endl;
	TString pdfName= Form("%s_MoreBadCellsCandidates.pdf",period.Data());
	//..plot the canvases of cells into a .pdf file
	for(Int_t can=0;can<nCv;can++)
	{
		//TString internal_pdfName1=pdfName+"Low.pdf";
		if(can==0)
		{
			//..first pad
			cCompAll[can]    ->Print(Form("%s(",pdfName.Data()));
		}
		else if(can==(nCv-1))//..last canvas
		{
			//..last pad
			cCompAll[can]    ->Print(Form("%s)",pdfName.Data()));
		}
		else
		{
			//..all pads in between
			cCompAll[can]    ->Print(Form("%s",pdfName.Data()));
		}
	}
}
///
/// Funtion to set TH1 histograms to a similar style
///
//________________________________________________________________________
void SetHisto(TH1 *Histo,TString Xtitel,TString Ytitel)
{
	Histo->SetStats(0);
	Histo->SetTitle("");

	Histo->GetYaxis()->SetTitleOffset(1.1);
	Histo->GetXaxis()->SetTitleOffset(1.1);
	Histo->GetXaxis()->SetLabelSize(0.05);
	Histo->GetYaxis()->SetLabelSize(0.05);
	Histo->GetXaxis()->SetTitleSize(0.06);
	Histo->GetYaxis()->SetTitleSize(0.06);
	Histo->GetXaxis()->SetNdivisions(505);
	Histo->GetYaxis()->SetNdivisions(505);

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
/// Function to set up canvas pads such that they share a common axis
///
//___________________________________________
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

   C->Divide(Nx,Ny,0.000,0.000);
   //..Top row
   for (Int_t i=0;i<Nx;i++)
   {
	   C->cd(i+1)->SetLogy();
	   gPad->SetLeftMargin(lMargin);
	   gPad->SetRightMargin(rMargin);
	   gPad->SetBottomMargin(0);
	   gPad->SetTopMargin(tMargin);
       //pad->Draw();??
  }
   //..Bottom row
   for (Int_t i=0;i<Nx;i++)
   {
	   C->cd(i+1+Nx);
	   gPad->SetLeftMargin(lMargin);
	   gPad->SetRightMargin(rMargin);
	   gPad->SetBottomMargin(bMargin);
	   gPad->SetTopMargin(0);
   }
   //..Top row
   for (Int_t i=0;i<Nx;i++)
   {
	   C->cd(i+1+2*Nx)->SetLogy();
	   gPad->SetLeftMargin(lMargin);
	   gPad->SetRightMargin(rMargin);
	   gPad->SetBottomMargin(0);
	   gPad->SetTopMargin(tMargin);
   }
   //..Bottom row
   for (Int_t i=0;i<Nx;i++)
   {
	   C->cd(i+1+3*Nx);
	   gPad->SetLeftMargin(lMargin);
	   gPad->SetRightMargin(rMargin);
	   gPad->SetBottomMargin(bMargin);
	   gPad->SetTopMargin(0);
   }
}
