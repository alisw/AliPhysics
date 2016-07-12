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

// --- ROOT system ---
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TString.h>
#include <TLatex.h>

#include <Riostream.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <alloca.h>
#include <string>
#include <cstring>

// --- ANALYSIS system ---
#include <AliCalorimeterUtils.h>
#include <AliEMCALGeometry.h>
#include <AliAnaCaloChannelAnalysis.h>
#include <AliAODEvent.h>

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnaCaloChannelAnalysis);
/// \endcond

///
/// Default constructor
///
//________________________________________________________________________
AliAnaCaloChannelAnalysis::AliAnaCaloChannelAnalysis():
		TObject(),
    fCurrentRunNumber(-1),fPeriod(),fPass(),fTrigger(),fNoOfCells(),
    fMergeOutput(),fAnalysisOutput(),fAnalysisInput(),fRunList(),
    fQADirect(), fMergedFileName(), fAnalysisVector(),
    fRunListFileName(),fWorkdir(),fTrial(),fExternalFileName(),
    fCaloUtils()
{
	fCurrentRunNumber = 254381;
	fPeriod           = "LHC16h";
	fPass             = "muon_caloLego";
	fTrigger          = "AnyINT";

	Init();
}

///
/// Constructor
///
//________________________________________________________________________
AliAnaCaloChannelAnalysis::AliAnaCaloChannelAnalysis(TString period, TString pass, TString trigger, Int_t RunNumber):
		TObject(),
    fCurrentRunNumber(-1),fPeriod(),fPass(),fTrigger(),fNoOfCells(),
    fMergeOutput(),fAnalysisOutput(),fAnalysisInput(),fRunList(),
    fQADirect(), fMergedFileName(), fAnalysisVector(),
    fRunListFileName(),fWorkdir(),fTrial(),fExternalFileName(),
    fCaloUtils()
{
	fCurrentRunNumber = RunNumber;
	fPeriod           = period;
	fPass             = pass;
	fTrigger          = trigger;

	Init();
}

///
/// Init default parameters
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::Init()
{
	//......................................................
	//..Default values - can be set by functions
	fWorkdir="./";
	fRunListFileName="runList.txt";
	fTrial = 0;
	fExternalFileName="";

	//..Settings for the input/output structure (hard coded)
	fAnalysisInput  ="AnalysisInput";
	fMergeOutput    ="ConvertOutput";
	fAnalysisOutput ="AnalysisOutput";
	//..Stuff for the convert function
	gSystem->mkdir(fMergeOutput);
	fMergedFileName= Form("%s/%s_%s_Merged.root", fMergeOutput.Data(), fPeriod.Data(),fPass.Data());
	fQADirect      = Form("CaloQA_%s",fTrigger.Data());
	fRunList       = Form("%s/%s/%s/%s", fAnalysisInput.Data(), fPeriod.Data(), fPass.Data(), fRunListFileName.Data());

	//.. make sure the vector is empty
	fAnalysisVector.clear();

	//......................................................
	//..Initialize EMCal/DCal geometry
	fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(fCurrentRunNumber);
	fCaloUtils->AccessGeometry(aod);
	//..Set the AODB calibration, bad channels etc. parameters at least once
	//fCaloUtils->AccessOADB(aod);

    fNoOfCells    =fCaloUtils->GetEMCALGeometry()->GetNCells(); //..Very important number, never change after that point!
    Int_t NModules=fCaloUtils->GetEMCALGeometry()->GetNumberOfSuperModules();
    fGoodCellID   =2377;  //..This is the ID of a good cell ELI where does this number come from ->write a setter function
    fCellStartDCal=12288; //..ELI this should be automatized from the geometry information!!

    //..This is how the calorimeter looks like in the current period (defined by example run ID fCurrentRunNumber)
	cout<<"Called geometry for run number: "<<fCurrentRunNumber<<endl;
	cout<<"Number of supermod: "<<NModules<<endl;
	cout<<"Number of cells: "<<fNoOfCells<<endl;
	cout<<"Cell ID of first DCal cell: "<<fCellStartDCal<<endl;
	//cout<<"Number of supermod utils: "<<fCaloUtils->GetNumberOfSuperModulesUsed()<<endl; //..will always be 22 unless set by hand

	//......................................................
	//..Initialize flag array to store how the cell is categorized
	//..In the histogram: bin 1= cellID 0, bin 2= cellID 1 etc
	//..In the array: fFlag[cellID]= some information
	fFlag   = new Int_t[fNoOfCells];
	fFlag[fNoOfCells]={0};  //..flagged as good by default

	//......................................................
	//..setings for the 2D histogram
	fNMaxCols    = 48;  //eta direction
	fNMaxRows    = 24;  //phi direction
	fNMaxColsAbs = 2*fNMaxCols;
	fNMaxRowsAbs = Int_t (NModules/2)*fNMaxRows; //multiply by number of supermodules
}

///
///	Main execution method.
///
/// First use Convert() to merge historgrams from a runlist .txt file.
/// The merged outputfile contains 3 different histograms.
/// In a second step analyse these merged histograms by calling BCAnalysis().
///
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::Run()
{
	TString inputfile = "";

	if(fExternalFileName=="")
	{
		cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;
		cout<<". . .Start process by converting files. . . . . . . . . . . ."<<endl;
		cout<<endl;
		fMergedFileName = Convert();
		cout<<endl;
	}
	else
	{
		cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;
		cout<<". . .Start process by loading external file. . . . . . . . . . ."<<endl;
		fMergedFileName = Form("%s/%s", fMergeOutput.Data(), fExternalFileName.Data());
	}

	cout<<". . .Load inputfile with name: "<<inputfile<<" . . . . . . . ."<<endl;
	cout<<". . .Continue process by . . . . . . . . . . . ."<<endl;
	cout<<endl;
	cout<<"o o o Bad channel analysis o o o"<<endl;


	BCAnalysis();


	cout<<endl;
	cout<<". . .End of process . . . . . . . . . . . . . . . . . . . . ."<<endl;
	cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;

}

///
/// Creates one file for the analysis from several QA output files listed in runlist.txt.
///
//________________________________________________________________________
TString AliAnaCaloChannelAnalysis::Convert()
{
	cout<<"o o o Start conversion process o o o"<<endl;
	cout<<"o o o period: " << fPeriod << ", pass: " << fPass << ",  trigger: "<<fTrigger<< endl;

	//..Create histograms needed for adding all the files together
	TH1D *hNEventsProcessedPerRun = new TH1D("hNEventsProcessedPerRun","Number of processed events vs run number",200000,100000,300000);
	//ELI a little problematic to hard code properties of histograms??
	TH2F *hCellAmplitude          = new TH2F("hCellAmplitude","Cell Amplitude",200,0,10,23040,0,23040);
	TH2F *hCellTime               = new TH2F("hCellTime","Cell Time",250,-275,975,23040,0,23040);

	//..Open the text file with the run list numbers and run index
	cout<<"o o o Open .txt file with run indices. Name = " << fRunList << endl;
	FILE *pFile = fopen(fRunList.Data(), "r");
	if(!pFile)cout<<"count't open file!"<<endl;
	Int_t Nentr;
	Int_t q;
	Int_t ncols;
	Int_t nlines = 0 ;
	Int_t RunId[500] ;
	while (1)
	{
		ncols = fscanf(pFile,"  %d ",&q);
		if (ncols< 0) break;
		RunId[nlines]=q;
		nlines++;
	}
	fclose(pFile);


	//..Open the different .root files with help of the run numbers from the text file
	const Int_t nRun = nlines ;
	TString base;
	TString infile;

	cout<<"o o o Start merging process of " << nRun <<" files"<< endl;
	//..loop over the amount of run numbers found in the previous text file.
	for(Int_t i = 0 ; i < nRun ; i++)
	{
		base  = Form("%s/%s/%s/%d", fAnalysisInput.Data(), fPeriod.Data(), fPass.Data(), RunId[i]);
		if ((fPass=="cpass1_pass2")||(fPass=="cfPass1-2"))
		{
			if (fTrigger=="default")
			{
				infile = Form("%s_barrel.root",base.Data());
			}
			else
			{
				infile = Form("%s_outer.root",base.Data());
			}
		}
		else
		{   //..This is a run2 case
			infile = Form("%s.root",base.Data()) ;
		}

		cout<<"    o Open .root file with name: "<<infile<<endl;
		TFile *f = TFile::Open(infile);

		//..Do some basic checks
		if(!f)
		{
			cout<<"Couldn't open/find .root file: "<<infile<<endl;
			continue;
		}
		TDirectoryFile *dir = (TDirectoryFile *)f->Get(fQADirect);
		if(!dir)
		{
			cout<<"Couln't open directory file in .root file: "<<infile<<", directory: "<<fQADirect<<endl;
			continue;
		}
		TList *outputList = (TList*)dir->Get(fQADirect);
		if(!outputList)
		{
			cout << "Couln't get TList from directory file: "<<fQADirect<<endl;
			continue;
		}

		//ELI should one maybe clone the hAmpId histos eg to hCellAmplitude, then one does't need to hard code them.
		TH2F *hAmpId;
		TH2F *hTimeId;
		TH2F *hNEvents;

		hAmpId =(TH2F *)outputList->FindObject("EMCAL_hAmpId");
		if(!hAmpId)
		{
			Printf("hAmpId not found");
			outputList->ls();
			continue;
		}
		hTimeId =(TH2F *)outputList->FindObject("EMCAL_hTimeId");
		if(!hTimeId)
		{
			Printf("hTimeId not found");
			outputList->ls();
			continue;
		}
		hNEvents =(TH2F *)outputList->FindObject("hNEvents");
		if(!hNEvents)
		{
			Printf("hNEvents not found");
			outputList->ls();
			continue;
		}
		Nentr =  (Int_t)hNEvents->GetEntries();

		//..does that mean do not merge small files?
		if (Nentr<100)
		{
			cout <<"    o File to small to be merged. Only N entries " << Nentr << endl;
			continue ;
		}
		cout <<"    o File with N entries " << Nentr<<" will be merged"<< endl;

		hNEventsProcessedPerRun->SetBinContent(RunId[i]-100000,(Double_t)Nentr);
		hCellAmplitude->Add(hAmpId);
		hCellTime->Add(hTimeId);

		outputList->Delete();
		dir->Delete();
		f->Close();
		delete f;
	}

	//.. Save the merged histograms
	cout<<"o o o Save the merged histogramms to .root file with name: "<<fMergedFileName<<endl;
	TFile *BCF = TFile::Open(fMergedFileName,"recreate");
	hNEventsProcessedPerRun->Write();
	hCellAmplitude->Write();
	hCellTime->Write();
	BCF->Close();
	cout<<"o o o End conversion process o o o"<<endl;
	return fMergedFileName;
}


/// Configure a complete analysis with different criteria, it provides bad+dead cells lists
/// You can manage criteria used and their order, the first criteria will use the original
/// output file from AliAnalysisTaskCaloCellsQA task,
//..Run over the list of stored period analysed settings (in fAnalysisVector) and
//..pass them to the function PeriodAnalysis()
//..The last PeriodAnalysis() is always executed with criteria 7
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::BCAnalysis()
{
	cout<<"o o o Bad channel analysis o o o"<<endl;

    TArrayD PeriodArray;
    for(Int_t i=0;i<fAnalysisVector.size();i++)
    {
		TFile::Open(fMergedFileName);
		PeriodArray=fAnalysisVector.at(i);
    		PeriodAnalysis(PeriodArray.At(0),PeriodArray.At(1),PeriodArray.At(2),PeriodArray.At(3));
    }

    //..In the end summarize results
    //..in a .pdf and a .txt file
	SummarizeResults();
	cout<<"o o o End of bad channel analysis o o o"<<endl;
}

///
/// Set period analysis parameters
///
/// \param criteria -- selection criteria (?)
/// \param nsigma   -- n sigma cut
/// \param emin     -- minimum energy
/// \param emax     -- maximum energy
///
///
/// This function does perform different checks depending on the given criterium variable
/// different possibilities for criterium are:
/// 1 : average E for E>emin
/// 2 : entries for E>emin
/// 3 : ki²/ndf  (from fit of each cell Amplitude between emin and emax)
/// 4 : A parameter (from fit of each cell Amplitude between emin and emax)
/// 5 : B parameter (from fit of each cell Amplitude between emin and emax)
/// 6 :
/// 7 : give bad + dead list
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::AddPeriodAnalysis(Int_t criteria, Double_t nsigma, Double_t emin, Double_t emax)
{
  TArrayD periodArray;
  periodArray.Set(4);
  periodArray.AddAt((Double_t)criteria,0);
  periodArray.AddAt(nsigma,1);
  periodArray.AddAt(emin,2);
  periodArray.AddAt(emax,3);
  fAnalysisVector.push_back(periodArray);
}

//..This function does perform different checks depending on the given criterium variable
//..diffrent possibilities for criterium are:
// 1 : average E for E>Emin and E<Emax
// 2 : entries for E>Emin and E<Emax

// 3 : ki²/ndf  (from fit of each cell Amplitude between Emin and Emax)
// 4 : A parameter (from fit of each cell Amplitude between Emin and Emax)
// 5 : B parameter (from fit of each cell Amplitude between Emin and Emax)
// 6 :
// 7 : give bad + dead list
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::PeriodAnalysis(Int_t criterum, Double_t nsigma, Double_t emin, Double_t emax)
{
	//ELI criterum should be between 1-4

	cout<<""<<endl;
	cout<<""<<endl;
	cout<<""<<endl;
	cout<<"o o o o o o o o o o o o o o o o o o o o o o  o o o"<<endl;
	cout<<"o o o PeriodAnalysis for flag "<<criterum<<" o o o"<<endl;
	cout<<"o o o Done in the energy range E "<<emin<<"-"<<emax<<endl;

	Int_t CellID, nb1=0, nb2=0;
	TString output;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. DEAD CELLS
	//.. Flage Dead cells with fFlag=1
	//.. this excludes cells from analysis (will not appear in results)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"o o o Flag Dead Cells o o o"<<endl;
	FlagAsDead();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. ANALYSIS OF CELLS WITH ENTRIES
	//.. Build average distributions and fit them
	//.. Three tests for bad cells:
	//.. 1) Average energy per hit
	//.. 2) Average hit per event
	//.. 3) ...
	//.. 4) ...
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TH1F* hisogram;
	if(criterum < 6)cout<<"o o o Analyze average cell distributions o o o"<<endl;
	//..For case 1 or 2
	if(criterum < 3)   hisogram = TestCellEandN(criterum, emin, emax,nsigma);
	//..For case 3, 4 or 5
	else if (criterum < 6) TestCellShapes(criterum, emin, emax, nsigma);

	Int_t dnbins = 200;
	if(criterum==1)              FlagAsBad(criterum, hisogram, nsigma, dnbins,-1);
	if(criterum==2 && emin==0.5) FlagAsBad(criterum, hisogram, nsigma, dnbins*9000,-1); //ELI I did massivley increase the binning now but it helps a lot
	if(criterum==2 && emin>0.5)  FlagAsBad(criterum, hisogram, nsigma, dnbins*17,-1);

	/*
	if(criterum==3)              FlagAsBad(criterum, hisogram, nsigma, dnbins, maxval3);
	if(criterum==4)              FlagAsBad(criterum, hisogram, nsigma, dnbins, maxval1);
	if(criterum==5)              FlagAsBad(criterum, hisogram, nsigma, dnbins, maxval2);
	 */

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. RESULTS
	//.. 1) Print the bad cells
	//..    and write the results to a file
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	//..Print the results on the screen and
	//..write the results in a file
	output.Form("%s/Criterion%d_Emin-%.2f_Emax-%.2f.txt", fAnalysisOutput.Data(), criterum,emin,emax);
	ofstream file(output, ios::out | ios::trunc);
	if(!file)
	{
		cout<<"#### Major Error. Check the textfile!"<<endl;
	}
	file<<"Criterion : "<<criterum<<", emin = "<<emin<<" GeV"<<", emax = "<<emax<<" GeV"<<endl;
	file<<"Bad by lower value : "<<endl;
	cout<<"    o bad cells by lower value (for cell E between "<<emin<<"-"<<emax<<")"<<endl;
	cout<<"      ";
	nb1=0;
	for(CellID=0;CellID<fNoOfCells;CellID++)
	{
		if(fFlag[CellID]==2)
		{
			nb1++;
			file<<CellID<<", ";
			cout<<CellID<<",";
		}
	}
	file<<"("<<nb1<<")"<<endl;
	cout<<"("<<nb1<<")"<<endl;
	file<<"Bad by higher value : "<<endl;
	cout<<"    o bad cells by higher value (for cell E between "<<emin<<"-"<<emax<<")"<<endl;
	cout<<"      ";
	nb2=0;
	for(CellID=0;CellID<fNoOfCells;CellID++)
	{
		if(fFlag[CellID]==3)
		{
			nb2++;
			file<<CellID<<", ";
			cout<<CellID<<",";
		}
	}
	file<<"("<<nb2<<")"<<endl;
	cout<<"("<<nb2<<")"<<endl;

	file<<"Total number of bad cells"<<endl;
	file<<"("<<nb1+nb2<<")"<<endl;
	file.close();
	cout<<"    o Total number of bad cells "<<endl;
	cout<<"      ("<<nb1+nb2<<")"<<endl;

}
//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
//.. 1) summarize all dead and bad cells in a text file
//.. 2) plot all bad cell E distributions in a .pdf file
//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::SummarizeResults()
{
	Int_t CellID, DCalCells=0, EMCalCells=0;
	TString CellSummaryFile, DeadPdfName, BadPdfName;

	DeadPdfName     = Form("%s/%s%sDC_SummaryResults_V%i.pdf", fAnalysisOutput.Data(), fPeriod.Data(), fPass.Data(), fTrial);
	BadPdfName      = Form("%s/%s%sBC_SummaryResults_V%i.pdf", fAnalysisOutput.Data(), fPeriod.Data(), fPass.Data(), fTrial);
	CellSummaryFile = Form("%s/%s%sBC_SummaryResults_V%i.txt", fAnalysisOutput.Data(), fPeriod.Data(), fPass.Data(), fTrial); ;
	cout<<"    o Final results o "<<endl;
	cout<<"    o write results into .txt file: "<<CellSummaryFile<<endl;
	cout<<"    o write results into .pdf file: "<<BadPdfName<<endl;
	ofstream file(CellSummaryFile, ios::out | ios::trunc);
	if(file)
	{
		file<<"Dead cells : "<<endl;
		cout<<"    o Dead cells : "<<endl;
		EMCalCells =0;
		DCalCells  =0;
		for(CellID=0; CellID<fNoOfCells; CellID++)
		{
			if(CellID==0)
			{
				file<<"In EMCal : "<<endl;
				cout<<"    o In EMCal : "<<endl;
			}
			if(CellID==fCellStartDCal)
			{
				file<<"In DCal : "<<endl;
				cout<<endl;
				cout<<"    o In DCal : "<<endl;
			}
			if(fFlag[CellID]==1)
			{
				file<<CellID<<"\n" ;
				cout<<CellID<<"," ;
				if(CellID<fCellStartDCal)EMCalCells++;
				else                     DCalCells++;
			}
		}
		cout<<endl;
		file<<"EMCal ("<<EMCalCells<<" ="<<100*EMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<DCalCells<<" ="<<100*DCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;
		cout<<"    o EMCal ("<<EMCalCells<<" ="<<100*EMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<DCalCells<<" ="<<100*DCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;

		file<<"Bad cells: "<<endl;
		cout<<"    o Bad cells: "<<endl;
		EMCalCells =0;
		DCalCells  =0;
		for(CellID=0;CellID<fNoOfCells;CellID++)
		{
			if(CellID==0)
			{
				file<<"In EMCal : "<<endl;
				cout<<"    o In EMCal : "<<endl;
			}
			if(CellID==fCellStartDCal)
			{
				file<<"In DCal : "<<endl;
				cout<<endl;
				cout<<"    o In DCal : "<<endl;
			}
			if(fFlag[CellID]>1)
			{
				file<<CellID<<"\n" ;
				cout<<CellID<<"," ;
				if(CellID<fCellStartDCal)EMCalCells++;
				else                     DCalCells++;
			}
		}
		cout<<endl;
		file<<"EMCal ("<<EMCalCells<<" ="<<100*EMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<DCalCells<<" ="<<100*DCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;
		cout<<"    o EMCal ("<<EMCalCells<<" ="<<100*EMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<DCalCells<<" ="<<100*DCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;
	}
	file.close();

	TFile::Open(fMergedFileName);
	//cout<<"    o Save the Dead channel spectra to a .pdf file"<<endl;
	//SaveBadCellsToPDF(0,DeadPdfName);
	cout<<"    o Save the bad channel spectra to a .pdf file"<<endl;
	SaveBadCellsToPDF(1,BadPdfName) ;
}

///
/// Allow to produce a pdf file with badcells candidates (red) compared to a refence cell (black).
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::SaveBadCellsToPDF(Int_t version, TString pdfName)
{
	//..version=0 ->Print dead cells
	//..version=1 ->print bad cells
	//..ELI can be refined with more versions!
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetFillColor(kWhite);
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetPalette(1);

	char title[100];
	char name[100];

	TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");
	TH1 *hCellref = hCellAmplitude->ProjectionX("badcells",fGoodCellID+1,fGoodCellID+1);
    Int_t FirstCanvas=0;
	//..collect cells in an internal vector.
	//..when the vector is of size=9 or at the end of being filled
	//..plot the channels into a canvas
	std::vector<Int_t> ChannelVector;
	ChannelVector.clear();
	for(Int_t cell=0;cell<fNoOfCells;cell++)
	{
		if(fFlag[cell]==1 && version==0)ChannelVector.push_back(cell);
		if(fFlag[cell]>1  && version==1)ChannelVector.push_back(cell);

		//..when 9 bad cells are collected or we are at the end of the list, fill the canvas
		if(ChannelVector.size()==9 || cell == fNoOfCells-1)
		{
			TString Internal_pdfName=pdfName;
			TCanvas *c1 = new TCanvas("badcells","badcells",1000,750);
			if(ChannelVector.size() > 6)        c1->Divide(3,3);
			else if (ChannelVector.size() > 3)  c1->Divide(3,2);
			else                                c1->Divide(3,1);

			TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
			for(Int_t i=0; i<ChannelVector.size() ; i++)
			{
				sprintf(name, "Cell %d",ChannelVector.at(i)) ;
				TH1 *hCell = hCellAmplitude->ProjectionX(name,ChannelVector.at(i)+1,ChannelVector.at(i)+1);
				sprintf(title,"Cell %d      Entries : %d  Ref : %d",ChannelVector.at(i), (Int_t)hCell->GetEntries(), (Int_t)hCellref->GetEntries() ) ;
				TString reflegend = Form("reference Cell %i",fGoodCellID);

				c1->cd(i%9 + 1);
				c1->cd(i%9 + 1)->SetLogy();
				hCell->SetLineColor(2);
				hCell->SetMaximum(1e6);
				hCell->SetAxisRange(0.,10.);
				hCell->GetXaxis()->SetTitle("E (GeV)");
				hCell->GetYaxis()->SetTitle("N Entries");
				hCell->SetLineWidth(1) ;
				hCell->SetTitle(title);
				hCellref->SetAxisRange(0.,8.);
				hCellref->SetLineWidth(1);
				hCellref->SetLineColor(1);

				if(i==0)
				{
					leg->AddEntry(hCellref,reflegend,"l");
					leg->AddEntry(hCell,"current","l");
				}
				hCell->Draw() ;
				hCellref->Draw("same") ;
				leg->Draw();
			}

			if(ChannelVector.size()<9 || cell == fNoOfCells-1)
			{
				Internal_pdfName +=")";
				//cout<<"Print canvas to file: "<<Internal_pdfName.Data()<<endl;
				c1->Print(Internal_pdfName.Data());
			}
			else if(FirstCanvas==0)
			{
				Internal_pdfName +="(";
				//cout<<"Print canvas to file: "<<Internal_pdfName.Data()<<endl;
				c1->Print(Internal_pdfName.Data());
				FirstCanvas=1;
			}
			else
			{
				//cout<<"Print canvas to file: "<<Internal_pdfName.Data()<<endl;
				c1->Print(Internal_pdfName.Data());
			}
			delete c1;
			delete leg;
			ChannelVector.clear();
		}
	}
	delete hCellref;
}

///
///  1) create a distribution for the input histogram;
///  2) fit the distribution with a gaussian;
///  3) define good area within +-nsigma to identfy badcells.
///
/// \param pflag   -- flag with channel criteria to be filled
/// \param inhisto -- input histogram;
/// \param dnbins  -- number of bins in distribution;
/// \param dmaxval -- maximum value on distribution histogram.
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::FlagAsBad(Int_t crit, TH1* inhisto, Double_t nsigma, Int_t dnbins, Double_t dmaxval)
{  
	gStyle->SetOptStat(1); // MG modif
	gStyle->SetOptFit(1);  // MG modif

	if(crit==1)cout<<"    o Fit average energy per hit distribution"<<endl;
	if(crit==2)cout<<"    o Fit average hit per event distribution"<<endl;

	Int_t CellColumn=0,CellRow=0;
	Int_t CellColumnAbs=0,CellRowAbs=0;
	Int_t Trash;

	TString HistoName=inhisto->GetName();
	Double_t goodmax= 0. ;
	Double_t goodmin= 0. ;
	if (dmaxval < 0.)
	{
		dmaxval = inhisto->GetMaximum()*1.01;  // 1.01 - to see the last bin
		if(crit==2 && dmaxval > 1) dmaxval =1. ;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .build the distribution of average values
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TH1 *distrib = new TH1F(Form("%sDistr",(const char*)HistoName), "", dnbins, inhisto->GetMinimum(), dmaxval);
	distrib->SetXTitle(inhisto->GetYaxis()->GetTitle());
	distrib->SetYTitle("Entries");

	//..build two dimensional histogram with values row vs. column
	TH2F *Plot2D = new TH2F(Form("%s_HitRowColumn",(const char*)HistoName),Form("%s_HitRowColumn",(const char*)HistoName),fNMaxColsAbs+2,-1.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+2,-1.5,fNMaxRowsAbs+0.5);
	Plot2D->GetXaxis()->SetTitle("cell column (#eta direction)");
	Plot2D->GetYaxis()->SetTitle("cell row (#phi direction)");

	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		//..Do that only for cell ids also accepted by the code
		if(!fCaloUtils->GetEMCALGeometry()->CheckAbsCellId(cell))continue;
		//ELI throw away the zeros if(inhisto->GetBinContent(c+1)!=0)
		//..fill the distribution of avarge cell values
		distrib->Fill(inhisto->GetBinContent(cell+1));

		//..Get Row and Collumn for cell ID c
		fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(cell,0,CellColumn,CellRow,Trash,CellColumnAbs,CellRowAbs);
		if(CellColumnAbs> fNMaxColsAbs || CellRowAbs>fNMaxRowsAbs)
		{
			cout<<"Problem! wrong calculated number of max col and max rows"<<endl;
			cout<<"current col: "<<CellColumnAbs<<", max col"<<fNMaxColsAbs<<endl;
			cout<<"current row: "<<CellRowAbs<<", max row"<<fNMaxRowsAbs<<endl;
		}
		Plot2D->SetBinContent(CellColumnAbs,CellRowAbs,inhisto->GetBinContent(cell+1));
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .draw histogram + distribution
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	//Produ
	TCanvas *c1 = new TCanvas(HistoName,HistoName,900,900);
	c1->ToggleEventStatus();
	TPad*    upperPad    = new TPad("upperPad", "upperPad",.005, .5, .995, .995);
	TPad*    lowerPadLeft = new TPad("lowerPadL", "lowerPadL",.005, .005, .5, .5);
	TPad*    lowerPadRight = new TPad("lowerPadR", "lowerPadR",.5, .005, .995, .5);
	upperPad->Draw();
	lowerPadLeft->Draw();
	lowerPadRight->Draw();

	upperPad->cd();
	upperPad->SetLeftMargin(0.045);
	upperPad->SetRightMargin(0.03);
	upperPad->SetLogy();
	inhisto->SetTitleOffset(0.6,"Y");
	inhisto->GetXaxis()->SetRangeUser(0,17000);

	inhisto->SetLineColor(kBlue+1);
	inhisto->Draw();

	lowerPadRight->cd();
	lowerPadRight->SetLeftMargin(0.09);
	lowerPadRight->SetRightMargin(0.06);
	Plot2D->Draw("colz");

	lowerPadLeft->cd();
	lowerPadLeft->SetLeftMargin(0.09);
	lowerPadLeft->SetRightMargin(0.06);
	lowerPadLeft->SetLogy();
	distrib->SetLineColor(kBlue+1);
	distrib->Draw();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .fit histogram
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	Int_t higherbin=0,i;
	for(i = 2; i <= dnbins; i++)
	{
		if(distrib->GetBinContent(higherbin) < distrib->GetBinContent(i))  higherbin = i ;
	}
	//..good range is around the max value as long as the
	//..bin content is larger than 2 entries
	for(i = higherbin ; i<=dnbins ; i++)
	{
		if(distrib->GetBinContent(i)<2) break ;
		goodmax = distrib->GetBinCenter(i);
	}
	for(i = higherbin ; i>1 ; i--)
	{
		if(distrib->GetBinContent(i)<2) break ;
		goodmin = distrib->GetBinLowEdge(i);
	}
	//cout<<"higherbin : "<<higherbin<<endl;
	//cout<<"good range : "<<goodmin<<" - "<<goodmax<<endl;

	TF1 *fit2 = new TF1("fit2", "gaus");
	//..start the fit with a mean of the highest value
	fit2->SetParameter(1,higherbin);

	distrib->Fit(fit2, "0LQEM", "", goodmin, goodmax);
	Double_t sig, mean, chi2ndf;
	// Marie midif to take into account very non gaussian distrig
	mean    = fit2->GetParameter(1);
	sig     = fit2->GetParameter(2);
	chi2ndf = fit2->GetChisquare()/fit2->GetNDF();

	if (mean <0.) mean=0.; //ELI is this not a highly problematic case??

	goodmin = mean - nsigma*sig ;
	goodmax = mean + nsigma*sig ;

	if (goodmin<0) goodmin=0.;

	cout<<"    o Result of fit: "<<endl;
	cout<<"    o  "<<endl;
	cout<<"    o Mean: "<<mean <<" sigma: "<<sig<<endl;
	cout<<"    o good range : "<<goodmin <<" - "<<goodmax<<endl;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .Add info to histogram
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TLine *lline = new TLine(goodmin, 0, goodmin, distrib->GetMaximum());
	lline->SetLineColor(kGreen+2);
	lline->SetLineStyle(7);
	lline->Draw();

	TLine *rline = new TLine(goodmax, 0, goodmax, distrib->GetMaximum());
	rline->SetLineColor(kGreen+2);
	rline->SetLineStyle(7);
	rline->Draw();

	TLegend *leg = new TLegend(0.60,0.82,0.9,0.88);
	leg->AddEntry(lline, "Good region boundary","l");
	leg->Draw("same");

	fit2->SetLineColor(kOrange-3);
	fit2->SetLineStyle(1);//7
	fit2->Draw("same");

	TLatex* text = 0x0;
	if(crit==1) text = new TLatex(0.2,0.8,Form("Good range: %.2f-%.2f",goodmin,goodmax));
	if(crit==2) text = new TLatex(0.2,0.8,Form("Good range: %.2f-%.2fx10^-5",goodmin*100000,goodmax*100000));
	text->SetTextSize(0.06);
	text->SetNDC();
	text->SetTextColor(1);
	//text->SetTextAngle(angle);
	text->Draw();
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .Save histogram
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	c1->Update();
	gSystem->mkdir(fAnalysisOutput);
	TString name   =Form("%s/criteria-_%d.gif", fAnalysisOutput.Data(), crit);
	if(crit==1)name=Form("%s/AverageEperHit_%s.gif", fAnalysisOutput.Data(), (const char*)HistoName);
	if(crit==2)name=Form("%s/AverageHitperEvent_%s.gif", fAnalysisOutput.Data(), (const char*)HistoName);
	c1->SaveAs(name);


	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . . Mark the bad cells in the pflag array
	//. . .(0= bad because cell average value lower than min allowed)
	//. . .(2= bad because cell average value higher than max allowed)
	//. . .(1 by default)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"    o Flag bad cells that are outside the good range "<<endl;
	for(Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		//cel=0 and bin=1, cel=1 and bin=2
		if (inhisto->GetBinContent(cell+1) <= goodmin && fFlag[cell]!=1) //ELI <= but not >=
		{
			 fFlag[cell]=2;
		}
		if (inhisto->GetBinContent(cell+1) > goodmax && fFlag[cell]!=1)
		{
			 fFlag[cell]=3;
		}
	}
	cout<<"    o "<<endl;

}

///
/// Average hit per event and the average energy per hit is caluclated for each cell.
///
//_________________________________________________________________________
TH1F* AliAnaCaloChannelAnalysis::TestCellEandN(Int_t crit, Double_t emin, Double_t emax, Double_t nsigma)
{
	TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");

	cout<<"    o Calculate average cell hit per event and average cell energy per hit "<<endl;
	//..binning parameters
	Int_t ncells  = hCellAmplitude->GetNbinsY();
	Double_t amin = hCellAmplitude->GetYaxis()->GetXmin();
	Double_t amax = hCellAmplitude->GetYaxis()->GetXmax();
	TH1F *Histogram;
	if(crit==1)Histogram = new TH1F(Form("hCellEtoNtotal_E%.2f-%.2f",emin,emax),Form("Average energy per hit, %.2f < E < %.2f GeV",emin,emax), ncells,amin,amax);
	if(crit==2)Histogram = new TH1F(Form("hCellNtotal_E%.2f-%.2f",emin,emax),Form("Number of hits per events, %.2f < E < %.2f GeV",emin,emax), ncells,amin,amax);
	Histogram->SetXTitle("AbsId");
	Histogram->SetYTitle("Av. hits per events");
	Histogram->GetXaxis()->SetNdivisions(505);

	TH1* hNEventsProcessedPerRun = (TH1*) gFile->Get("hNEventsProcessedPerRun");
	Double_t totalevents = hNEventsProcessedPerRun->Integral(1, hNEventsProcessedPerRun->GetNbinsX());

	//..here the average hit per event and the average energy per hit is caluclated for each cell.
	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		Double_t Esum = 0;
		Double_t Nsum = 0;

		for (Int_t j = 1; j <= hCellAmplitude->GetNbinsX(); j++)
		{
			Double_t E = hCellAmplitude->GetXaxis()->GetBinCenter(j);
			Double_t N = hCellAmplitude->GetBinContent(j, cell+1);
			//..investigate only cells that were not flagged as dead ore bad
			if (E < emin || E > emax || fFlag[cell]!=0) continue;
			Esum += E*N;
			Nsum += N;
		}
		if(totalevents> 0. && crit==2)Histogram->SetBinContent(cell+1, Nsum/totalevents);  //..number of hits per event
		if(Nsum > 0.       && crit==1)Histogram->SetBinContent(cell+1, Esum/Nsum);         //..average energy per hit
	}
	delete hCellAmplitude;

	return Histogram;
}

/// ELI this method is currently not used
/// Test cells shape using fit function f(x)=A*exp(-B*x)/x^2.
/// Produce values per cell + distributions for A,B and chi2/ndf parameters.
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::TestCellShapes(Int_t crit, Double_t fitemin, Double_t fitemax, Double_t nsigma)
{
	TString hname= "hCellAmplitude";
	Int_t dnbins = 1000;
	TH2 *hCellAmplitude = (TH2*) gFile->Get(Form("%s",(const char*)hname));

	// binning parameters
	Int_t  ncells = hCellAmplitude->GetNbinsY();
	Double_t amin = hCellAmplitude->GetYaxis()->GetXmin();
	Double_t amax = hCellAmplitude->GetYaxis()->GetXmax();
	cout << "ncells " << ncells << " amin = " << amin << "amax = " << amax<< endl;

	// initialize histograms
	TH1 *hFitA = new TH1F(Form("hFitA_%s",(const char*)hname),"Fit A value", ncells,amin,amax);
	hFitA->SetXTitle("AbsId");
	hFitA->SetYTitle("A");

	TH1 *hFitB = new TH1F(Form("hFitB_%s",(const char*)hname),"Fit B value", ncells,amin,amax);
	hFitB->SetXTitle("AbsId");
	hFitB->SetYTitle("B");

	TH1 *hFitChi2Ndf = new TH1F(Form("hFitChi2Ndf_%s",(const char*)hname),"Fit #chi^{2}/ndf value", ncells,amin,amax);
	hFitChi2Ndf->SetXTitle("AbsId");
	hFitChi2Ndf->SetYTitle("#chi^{2}/ndf");

	Double_t maxval1=0., maxval2=0., maxval3=0.;
	Double_t prev=0., MSA=0., AvA = 0. ; //those param are used to automaticaly determined a reasonable maxval1
	Double_t prev2=0., MSB=0., AvB = 0.  ; //those param are used to automaticaly determined a reasonable maxval2
	Double_t prev3=0., MSki2=0., Avki2 = 0. ; //those param are used to automaticaly determined a reasonable maxval3
	Double_t ki2=0.0 ;
	for (Int_t k = 1; k <= ncells; k++)
	{
		TF1 *fit = new TF1("fit", "[0]*exp(-[1]*x)/x^2");
		TH1 *hCell = hCellAmplitude->ProjectionX("",k,k);
		if (hCell->GetEntries() == 0) continue;
		// hCell->Rebin(3);
		hCell->Fit(fit, "0QEM", "", fitemin, fitemax);
		delete hCell;

		if(fit->GetParameter(0) < 5000.)
		{
			hFitA->SetBinContent(k, fit->GetParameter(0));
			if(k<3000)
			{
				AvA +=  fit->GetParameter(0);
				if(k==2999)  maxval1  = AvA/3000. ;
				if (prev < fit->GetParameter(0)) MSA += fit->GetParameter(0) - prev;
				else MSA -= (fit->GetParameter(0) - prev) ;
				prev = fit->GetParameter(0);
			}
			else
			{
				if((fit->GetParameter(0) - maxval1) > 0. && (fit->GetParameter(0) - maxval1) < (MSA/1000.))
				{
					maxval1 = fit->GetParameter(0);
				}
			}
		}
		else hFitA->SetBinContent(k, 5000.);

		if(fit->GetParameter(1) < 5000.)
		{
			hFitB->SetBinContent(k, fit->GetParameter(1));
			if(k<3000)
			{
				AvB +=  fit->GetParameter(1);
				if(k==2999)  maxval2  = AvB/3000. ;
				if (prev2 < fit->GetParameter(1)) MSB += fit->GetParameter(1) - prev2;
				else MSB -= (fit->GetParameter(1) - prev2) ;
				prev2 = fit->GetParameter(1);
			}
			else
			{
				if((fit->GetParameter(1) - maxval2) > 0. && (fit->GetParameter(1) - maxval2) < (MSB/1000.))
				{
					maxval2 = fit->GetParameter(1);
				}
			}
		}
		else hFitB->SetBinContent(k, 5000.);


		if (fit->GetNDF() != 0 ) ki2 =  fit->GetChisquare()/fit->GetNDF();
		else ki2 = 1000.;

		if(ki2 < 1000.)
		{
			hFitChi2Ndf->SetBinContent(k, ki2);
			if(k<3000)
			{
				Avki2 +=  ki2;
				if(k==2999)  maxval3  = Avki2/3000. ;
				if (prev3 < ki2) MSki2 += ki2 - prev3;
				else MSki2 -= (ki2 - prev3) ;
				prev3 = ki2;
			}
			else
			{
				if((ki2 - maxval3) > 0. && (ki2 - maxval3) < (MSki2/1000.))
				{
					maxval3 = ki2;
				}
			}
		}
		else hFitChi2Ndf->SetBinContent(k, 1000.);

		delete fit ;
	}

	delete hCellAmplitude;

	// if you have problem with automatic parameter :
	//  maxval1 =
	//  maxval2 =
	//  maxval3 =
	/*
	if(crit==3)
		Process(crit, hFitChi2Ndf, nsigma, dnbins, maxval3);
	if(crit==4)
		Process(crit, hFitA, nsigma, dnbins,  maxval1);
	if(crit==5)
		Process(crit, hFitB, nsigma, dnbins, maxval2);
	*/

	//ELI something like thie in the future: return histogram;
}


/// --
/// This function finds cells with zero entries
/// It flags them by setting the fFlag content of the cell to 1.
///
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::FlagAsDead()
{
	TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");
	Int_t SumOfExcl=0;
	//..Direction of cell ID
	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		Double_t Nsum = 0;
		//..Direction of amplitude
		for (Int_t amp = 1; amp <= hCellAmplitude->GetNbinsX(); amp++)
		{
			//..cellID+1 = histogram bin
			Double_t N = hCellAmplitude->GetBinContent(amp,cell+1);
			Nsum += N;
		}
		//..If the amplitude in one cell is basically 0
		//..mark the cell as excluded
		//ELI I just wonder how you can have less than one but more than 0.5
		//shouldnt everything below 1 be excluded?
		if(Nsum >= 0.5 && Nsum < 1)cout<<"-----------------------small but non zero!!!!"<<endl;
		if(Nsum < 0.5 && Nsum != 0)cout<<"-----------------------non zero!!!!"<<endl;

		if(Nsum < 0.5 && fFlag[cell]==0)
		{
			//..Cell flagged as dead.
			//..Flag only if it hasn't been flagged before
			fFlag[cell]   =1;
			SumOfExcl++;
		}
	}
	delete hCellAmplitude;
	cout<<"    o Number of dead cells: "<<SumOfExcl<<endl;
	cout<<"     ("<<SumOfExcl<<")"<<endl;
}
