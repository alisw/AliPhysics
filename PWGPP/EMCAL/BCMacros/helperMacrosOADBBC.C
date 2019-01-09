/// \file helperMacrosOADBBC.C
/// \ingroup EMCALOfflineMacros
/// \brief Macro to test the files on OADB
///
///
/// See https://twiki.cern.ch/twiki/bin/view/ALICE/BadChannelAnalysis for general documentation
///
/// ---------------------
///  Running the macro
/// ---------------------
/// use root -b to speed up (no canvas drawn)              <br>
/// IMPORTANT dont forget to update your local OADB copy!  <br>
/// rsync -av --delete yourCERNuserName@lxplus.cern.ch:/eos/experiment/alice/analysis-data/ $ALICE_DATA <br>
/// To run the macro you have to compile it first          <br>
/// root [0] .L $ALICE_WORK_DIR/../ali-master/AliPhysics/PWGPP/EMCAL/BCMacros/helperMacrosOADBBC.C++ <br>
/// or                                                     <br>
/// root [0] .L $ALICE_WORK_DIR/../AliPhysics/PWGPP/EMCAL/BCMacros/helperMacrosOADBBC.C++ <br>
/// Then you can execute the function in it                <br>
/// This is a function for everyone interested:            <br>
/// root [1] Plot_BCMap("runList.txt")   In case your local copy is in $ALICE_DATA/OADB/EMCAL/..   <br>
/// root [1] Plot_BCMap("runList.txt","pathToYourOtherOADBCopy")  If you have an OADB copy other than in $ALICE_DATA/ you can add an explicit path here <br>
/// These are functions for BC exerts:                       <br>
/// root [1] Test_OADB("LHC16k",804,0,"runList16k.txt")     <br>
/// root [1] Plot_CellList("LHC15o",771,"List115o.txt") <br>
///
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
#include <TROOT.h>


// --- ANALYSIS system ---
#include "AliEMCALGeometry.h"          //include when compile
#include "AliCalorimeterUtils.h"       //include when compile
#include "AliAODEvent.h"               //include when compile
#include "AliOADBContainer.h"          //include when compile
#include "AliDataFile.h"               //include when compile

void SetHisto(TH1 *Histo,TString Xtitel,TString Ytitel);
void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2, Float_t lMargin = 0.15, Float_t rMargin = 0.05, Float_t bMargin = 0.15, Float_t tMargin = 0.05);

///
/// -Function for everyone who wants to plot bad maps-
/// Test how the bad channel map looks like for a given run
///________________________________________________________________________
Bool_t Plot_BCMap(TString runList="", TString pathOADB ="")
{
    gStyle->SetOptStat(0); //..Do not plot stat boxes
	//......................................................
	// Test if the file committed to OADB is correct
	//......................................................
	Int_t nSM = 20;
	TH1C *h[20];
	Int_t cellColumnAbs,cellRowAbs,trash,cellID;
	TString summaryPDF=Form("./OADB_Summary.pdf");

	//......................................................
	//..open the text file and save the run IDs into the RunId vector
	cout<<"  o Open .txt file with run indices. Name = " << runList << endl;
	TString RunPath        = Form("./%s",runList.Data());
	FILE *pFile = fopen(RunPath.Data(), "r");
	if(!pFile)
	{
		cout<<"couldn't open file "<<RunPath<<"!"<<endl;
		return 0;
	}
	Int_t q;
	Int_t ncols;
	Int_t nlines = 0 ;
	std::vector<Int_t> RunIdVec;
	while (1)
	{
		ncols = fscanf(pFile,"  %d ",&q);
		if (ncols< 0) break;
		RunIdVec.push_back(q);
		nlines++;
	}
	fclose(pFile);
	//..sort the vector by size to be shure to use the right order
	std::sort (RunIdVec.begin(), RunIdVec.end());
    Int_t nRuns=RunIdVec.size();

	//......................................................
	//..Get the OADB information
    //$ALICE_DATA
    cout<<"  o Alice path: "<<pathOADB<<endl;
    cout<<"  o Alice path+: "<<pathOADB<<"/OADB/EMCAL/EMCALBadChannels.root"<<endl;
    //pathOADB="."; //..locally from OADB commit/e-mail to test quickly

    AliOADBContainer *cont=new AliOADBContainer("");
    if(pathOADB!="")
    {
    		//..open the .root file directly to see if one can access it
    		TFile *fbad=new TFile(Form("%s/OADB/EMCAL/EMCALBadChannels.root",pathOADB.Data()),"read");
    		if (!fbad || fbad->IsZombie())
    		{
    			cout<<"couldn't find OADB container with help of pathOADB !"<<endl;
    			return 0;
    		}
    		else delete fbad;
    		cont->InitFromFile(Form("%s/OADB/EMCAL/EMCALBadChannels.root",pathOADB.Data()),"AliEMCALBadChannels");
    	}
    	else
    	{
    		//..open the .root file directly to see if one can access it
    		TFile *fbad=new TFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALBadChannels.root").data(),"read");
    		if (!fbad || fbad->IsZombie())
    		{
    			cout<<"OADB/EMCAL/EMCALBadChannels.root was not found !"<<endl;
    			return 0;
    		}
    		else delete fbad;
    		cont->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALBadChannels.root").data(),"AliEMCALBadChannels");
    }
    //......................................................
	//..Initialize EMCal/DCal geometry
	AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(RunIdVec.at(0));
	fCaloUtils->AccessGeometry(aod);
	AliEMCALGeometry * geom = fCaloUtils->GetEMCALGeometry();
	cout<<"  o Geometry loaded for Run "<<RunIdVec.at(0)<<endl;

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
	//*OADB figure that  looks like Marcels figures*/	TH2F *plot2D_Bad_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
	plot2D_Bad_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_Bad_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");
	histoName = Form("2DChannelMap_Flag_Dead");
	TH2F *plot2D_Dead_OADB = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	//*OADB figure that looks like Marcels figures*/	TH2F *plot2D_Dead_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
	plot2D_Dead_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_Dead_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");
	histoName = Form("2DChannelMap_Flag_Good");
	TH2F *plot2D_Good_OADB = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	//*OADB figure that looks like Marcels figures*/	TH2F *plot2D_Good_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
	plot2D_Good_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_Good_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");

	//.......................................................
	//.. Read the Bad Channel map from OADB
	TCanvas* C2 = new TCanvas("CompareCanvas","Compare OADB to LocalFile",1400,500);
	C2->Divide(3);

	cout<<"  o Checking "<<nRuns<<" runs: "<<endl;
	std::vector<Int_t> RunsWithoutMap;
	std::vector<Int_t> RunsWithMap;
	//Int_t dummyRun = 000276462;
	Int_t dummyRun = 282367;
	for(Int_t iRun = 0; iRun < nRuns; iRun++)
	{
		//cout<<"run number: "<<RunIdVec.at(iRun)<<endl;
		if(iRun%5==0 && iRun!=0)cout<<"."<<flush;
		if(iRun%20==0)cout<<"  o Run No."<<RunIdVec.at(iRun)<<endl;
		TObjArray *recal=(TObjArray*)cont->GetObject(RunIdVec.at(iRun));
		if(!recal)
		{
			cout<<"Error - No bad map for run Number "<<RunIdVec.at(iRun)<<" online!!"<<endl;
			cout<<"Check if you have the latest EOS version and/or check the BadChannelTwiki"<<endl;
			RunsWithoutMap.push_back(RunIdVec.at(iRun));
			continue;
		}
		RunsWithMap.push_back(RunIdVec.at(iRun));

		//..clear content from previous loop
		plot2D_Bad_OADB ->Reset();
		plot2D_Dead_OADB->Reset();
		plot2D_Good_OADB->Reset();
		for(Int_t iSM = 0; iSM < nSM; iSM++)
		{
			//..Bad map for each module
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
					if(h[iSM]->GetBinContent(column,row)==0)//..Good
					{
						plot2D_Good_OADB->SetBinContent(cellColumnAbs,cellRowAbs,1);
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
		plot2D_Good_OADB->SetTitle("Good Cells OADB");
		plot2D_Good_OADB->DrawCopy("colz");
		text->Draw();
		C2->cd(2);
		plot2D_Dead_OADB->SetTitle("Dead Cells OADB");
		plot2D_Dead_OADB->DrawCopy("colz");
		C2->cd(3);
		plot2D_Bad_OADB->SetTitle("Bad Cells OADB");
		plot2D_Bad_OADB->DrawCopy("colz");

		//..................................................................
        //..Save to PDF
		//..Add figures to the summary canvas
		if(iRun==0)             C2   ->Print(Form("%s(",summaryPDF.Data()));
		else if (iRun==nRuns-1) C2   ->Print(Form("%s)",summaryPDF.Data()));
		else                    C2   ->Print(Form("%s",summaryPDF.Data()));
	}//..end of run loop

	//..................................................................
	//..Summarize how many bad maps were found
	cout<<"==Total Summary=="<<endl;
	//..print run numbers without a correct map.
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

    return 1;
}

///
/// -Function for Bad Channel experts-
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
	//..open the text file and save the run IDs into the RunId vector
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
	std::vector<Int_t> RunIdVec;
	while (1)
	{
		ncols = fscanf(pFile,"  %d ",&q);
		if (ncols< 0) break;
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
/// -Function for Bad Channel experts-
/// After OADB maps are online you will get
/// suggestions from other users that think certain
/// channels should be masked. You can look at them here
///________________________________________________________________________
void Plot_CellList(TString period="LHC15n",Int_t trainNo=603,TString cellList="")
{
	gROOT->SetBatch(1); //..Prevent ROOT from stealing focus when plotting
	gStyle->SetPadTopMargin(0.05);//0.05
	gStyle->SetPadBottomMargin(0.18);//0.15
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadLeftMargin(0.17);
	gStyle->SetFrameFillColor(10);


	//......................................................
	//..Settings
	const Int_t nBlocks=3; //..number of different runblocks of the period

	Int_t zoomRange=20;
	Double_t range=5; //10 GeV
	Double_t tRange1,tRange2;
	tRange1=-200;
    tRange2=400;
	//..uncalibrated:
	//tRange1 =400;
	//tRange2 =900;
	Int_t shift=0;

	//......................................................
	//..open the text file and save the run IDs into the RunId vector
	TString ListName  = Form("./AnalysisOutput/%s/Train_%i/%s",period.Data(),trainNo,cellList.Data());
	cout<<"o o o Open .txt suggested cell IDs. Name = " << ListName << endl;
	FILE *pFile = fopen(ListName.Data(), "r");
	if(!pFile)
	{
		cout<<"couldn't open file "<<ListName<<"!"<<endl;
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
    cout<<"o found "<<cellIdVec.size()<<" cells in list"<<endl;
	//......................................................
	//..Get the .root file with the original histogram to compare if they coincide
	TString path        = Form("./AnalysisOutput/%s/Train_%i",period.Data(),trainNo);
	TString rootFileName[nBlocks];
	TFile* outputRoot[nBlocks];
	TH2F* h2DAmp[nBlocks];
	TH2F* h2DRatio[nBlocks];
	TH2F* h2DCellTime[nBlocks];
	//TString rootFileName= Form("Version%s/%s_INT7_Histograms_V%s.root",period.Data(),version.Data());
	for(Int_t iBlock=0;iBlock<nBlocks;iBlock++)
	{
		//..mostly set by hand
		//..15l
		/*
		if(iBlock==0)rootFileName[iBlock]= Form("Version1OADB/%s_INT7_Histograms_V1.root",period.Data());
		if(iBlock==1)rootFileName[iBlock]= Form("Version2OADB/%s_INT7_Histograms_V2.root",period.Data());
		if(iBlock==2)rootFileName[iBlock]= Form("Version3OADB/%s_INT7_Histograms_V3.root",period.Data());
		if(iBlock==3)rootFileName[iBlock]= Form("Version4OADB/%s_INT7_Histograms_V4.root",period.Data());
		 */
		//..15o
		/*if(iBlock==0)rootFileName[iBlock]= Form("Version1ManMasked/%s_INT7_Histograms_V1.root",period.Data());
		if(iBlock==1)rootFileName[iBlock]= Form("Version2ManMasked/%s_INT7_Histograms_V2.root",period.Data());
		if(iBlock==2)rootFileName[iBlock]= Form("Version3ManMasked/%s_INT7_Histograms_V3.root",period.Data());
		if(iBlock==3)rootFileName[iBlock]= Form("Version4ManMasked/%s_INT7_Histograms_V4.root",period.Data());
        */
		/*
		//..16s
		if(iBlock==0)rootFileName[iBlock]= Form("Version1OADB/%s_INT7_Histograms_V1.root",period.Data());
		if(iBlock==1)rootFileName[iBlock]= Form("Version2OADB/%s_INT7_Histograms_V2.root",period.Data());
		if(iBlock==2)rootFileName[iBlock]= Form("Version4OADB/%s_INT7_Histograms_V4.root",period.Data());
		 */
		//..16k
		if(period=="LHC16i")
		{
		if(iBlock==0)rootFileName[iBlock]= Form("Version0/%s_INT7_Histograms_V0.root",period.Data());
		}
		//..17n
		if(period=="LHC17n")
		{
		if(iBlock==0)rootFileName[iBlock]= Form("Version2/%s_INT7_Histograms_V2.root",period.Data());
		}
		if(period=="LHC16i")
		{
			if(iBlock==0)rootFileName[iBlock]= "Version2/LHC16i_INT7_Histograms_V1.root";
			if(iBlock==1)rootFileName[iBlock]= "Version2/LHC16i_INT7_Histograms_V2.root";
		}
		if(period=="LHC16g")
		{
			if(iBlock==0)rootFileName[iBlock]= "LHC16g_Block1_INT7_Histograms_V0.root";
			if(iBlock==1)rootFileName[iBlock]= "LHC16g_Block2_INT7_Histograms_V0.root";
			if(iBlock==2)rootFileName[iBlock]= "LHC16g_Block3_INT7_Histograms_V0.root";
		}
		cout<<"Look for filename: "<<rootFileName[iBlock]<<endl;
		cout<<"in path: "<<path<<endl;
		outputRoot[iBlock]   = TFile::Open(Form("%s/%s",path.Data(),rootFileName[iBlock].Data()));
		if(!outputRoot[iBlock])cout<<"File "<<outputRoot[iBlock]->GetName()<<" does not exist"<<endl;
		else cout<<"open File "<<outputRoot[iBlock]->GetName()<<endl;
		h2DAmp[iBlock]     =(TH2F*)outputRoot[iBlock]->Get("hCellAmplitude");
		h2DRatio[iBlock]   =(TH2F*)outputRoot[iBlock]->Get("ratio2DAmp");
		h2DCellTime[iBlock]=(TH2F*)outputRoot[iBlock]->Get("hCellTime");
		if(!h2DAmp[iBlock])cout<<"Problem 1"<<endl;
		if(!h2DRatio[iBlock])cout<<"Problem 2"<<endl;
		if(!h2DCellTime[iBlock])cout<<"Problem 3"<<endl;
	}

	TCanvas *c1 = new TCanvas(1);
	c1->Divide(2);
	c1->cd(1);
	SetHisto(h2DAmp[0],"","");
	h2DAmp[0]->Draw("colz");
	c1->cd(2);
	SetHisto(h2DRatio[0],"","");
	h2DRatio[0]->Draw("colz");

	//.. be aware of the special sturcture of the canvas
	//.. canvas has totalperCv*2 pads
	Int_t totalperCv = 5;
//	Int_t totalperCv = 8;
	Int_t nPad = TMath::Sqrt(totalperCv);
	Int_t nCv  = nCells/totalperCv+1;
	if(nCv<1)nCv=1;

	cout<<"    o create: "<<nCv<<" Canvases with "<<nPad*nPad<<" pads"<<endl;
	//..to compare specific cells over the runs
	TCanvas **cCompAll      = new TCanvas*[nCv];
	for(Int_t i=0;i<nCv;i++)
	{
		cCompAll[i] = new TCanvas(TString::Format("CompareGoodAll%d", i), TString::Format("V) Both (%d/%d)", i+1, nCv), 1200,600);
//		CanvasPartition(cCompAll[i],4,totalperCv/2,0.15,0.02,0.13,0.05);
		CanvasPartition(cCompAll[i],5,3,0.15,0.02,0.13,0.05);
	}
	TLegend *leg2 = new TLegend(0.60,0.60,0.9,0.85);
	cout<<"    o Fill Canvases with bad cells histograms"<<endl;

	TH1D *htmpCellAllRuns[nBlocks];
	TH1D *htmpCellTimeRuns[nBlocks];
	TH1D *htmpCellRatioAllRuns[nBlocks];
	cout<<"cell number: "<<endl;
	for(Int_t icell = 0; icell < nCells; icell++)
	{
		Int_t cellID=cellIdVec.at(icell);
		cout<<cellID<<", "<<flush;

		for(Int_t iBlock=0;iBlock<nBlocks;iBlock++)
		{
			htmpCellAllRuns[iBlock]       =h2DAmp[iBlock]  ->ProjectionX(TString::Format("hIDProj%i_cell%d",iBlock, cellID), cellID+1, cellID+1);
			SetHisto(htmpCellAllRuns[iBlock],Form("Energy [cell %i]",cellID),"Number of Hits/Events");
			htmpCellRatioAllRuns[iBlock]  =h2DRatio[iBlock]->ProjectionX(TString::Format("hIDRProj%i_cell%d", iBlock, cellID), cellID+1, cellID+1);
			SetHisto(htmpCellRatioAllRuns[iBlock],Form("Energy [cell %i]",cellID),"No. Hits/av. No. Hits");

			htmpCellTimeRuns[iBlock]  =h2DCellTime[iBlock]->ProjectionX(TString::Format("hIDTimeProj%i_cell%d", iBlock, cellID), cellID+1, cellID+1);
			SetHisto(htmpCellTimeRuns[iBlock],Form("Time, ns [cell %i]",cellID),"Entries/Events");

			//..Amplitude
			cCompAll[(icell)/totalperCv]->cd(((icell)%5)+1)->SetLogy();
			if(iBlock==0)htmpCellAllRuns[iBlock]->GetXaxis()->SetRangeUser(0,range);
			if(iBlock==0)htmpCellAllRuns[iBlock]->Draw("hist");
			if(iBlock==1)htmpCellAllRuns[iBlock]->SetLineColor(kBlue-7);
			if(iBlock==2)htmpCellAllRuns[iBlock]->SetLineColor(kGreen-2);
			if(iBlock==3)htmpCellAllRuns[iBlock]->SetLineColor(kViolet-1);
			if(iBlock>0)htmpCellAllRuns[iBlock]->DrawCopy("same hist");

			if(icell==0)
			{
				leg2->AddEntry(htmpCellAllRuns[iBlock],Form("Block%i",iBlock),"l");
				leg2->SetTextSize(0.07);
				leg2->SetBorderSize(0);
				leg2->SetFillColorAlpha(10, 0);
				leg2->Draw("same");
			}
			else if((icell)%5==0)
			{
				leg2->Draw("same");
			}

			//..Ratio
			cCompAll[(icell)/totalperCv]->cd(((icell)%5)+6)->SetLogy();
			if(iBlock==0)htmpCellRatioAllRuns[iBlock]->GetXaxis()->SetRangeUser(0,range);
			if(iBlock==0)htmpCellRatioAllRuns[iBlock]->Draw("hist");
			if(iBlock==1)htmpCellRatioAllRuns[iBlock]->SetLineColor(kBlue-7);
			if(iBlock==2)htmpCellRatioAllRuns[iBlock]->SetLineColor(kGreen-2);
			if(iBlock==3)htmpCellRatioAllRuns[iBlock]->SetLineColor(kViolet-1);
			if(iBlock>0)htmpCellRatioAllRuns[iBlock]->DrawCopy("same hist");

			//..Time
			cCompAll[(icell)/totalperCv]->cd(((icell)%5)+11)->SetLogy();
			if(iBlock==0)htmpCellTimeRuns[iBlock]->GetXaxis()->SetRangeUser(tRange1,tRange2);
			if(iBlock==0)htmpCellTimeRuns[iBlock]->Draw("hist");
			if(iBlock==1)htmpCellTimeRuns[iBlock]->SetLineColor(kBlue-7);
			if(iBlock==2)htmpCellTimeRuns[iBlock]->SetLineColor(kGreen-2);
			if(iBlock==3)htmpCellTimeRuns[iBlock]->SetLineColor(kViolet-1);
			if(iBlock>0)htmpCellTimeRuns[iBlock]->DrawCopy("same hist");
		}
	}
	cout<<endl;
//	TString pdfName= Form("%s_MoreBadCellsCandidates.pdf",period.Data());
	TString pdfName= Form("./AnalysisOutput/%s/Train_%i/%s_MoreBadCellsCandidates.pdf",period.Data(),trainNo,period.Data());
	//..plot the canvases of cells into a .pdf file
	for(Int_t can=0;can<nCv;can++)
	{
		//TString internal_pdfName1=pdfName+"Low.pdf";
		if(can==0)
		{
			//..first pad
			if(nCv>1)
			{
				cCompAll[can]    ->Print(Form("%s(",pdfName.Data()));
			}
			else
			{
				cCompAll[can]    ->Print(Form("%s",pdfName.Data()));
			}
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
/// Function to set TH1 histograms to a similar style
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
	   //gPad->SetBottomMargin(0);
	   gPad->SetBottomMargin(bMargin);
	   gPad->SetTopMargin(tMargin);
   }
 /*

   //..Bottom row
   for (Int_t i=0;i<Nx;i++)
   {
	   C->cd(i+1+3*Nx);
	   gPad->SetLeftMargin(lMargin);
	   gPad->SetRightMargin(rMargin);
	   gPad->SetBottomMargin(bMargin);
	   gPad->SetTopMargin(0);
   }
*/
}
