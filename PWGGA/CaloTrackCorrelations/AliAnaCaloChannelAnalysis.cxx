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
#include <TSpectrum.h>
#include <TROOT.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TList.h>
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

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::flush;
using std::ios;

/// \cond CLASSIMP
ClassImp(AliAnaCaloChannelAnalysis);
/// \endcond

///
/// Default constructor
///
//________________________________________________________________________
AliAnaCaloChannelAnalysis::AliAnaCaloChannelAnalysis():
TObject(),
	fCurrentRunNumber(-1),
	fPeriod(),
	fPass(),
	fTrigger(),
	fNoOfCells(),
	fCellStartDCal(12288),
	fAnalysisOutput(),
	fAnalysisInput(),
	fRunList(),
	fRunListFileName(),
	fWorkdir(),
	fQADirect(),
	fMergedFileName(),
	fAnalysisVector(),
	fTrial(),
	fExternalFileName(),
	fTestRoutine(),
	fNMaxCols(48),
	fNMaxRows(24),
	fNMaxColsAbs(),
	fNMaxRowsAbs(),
	fFlag(),
	fCriterionCounter(),
	fWarmCell(),
	fCaloUtils(),
	fRootFile(),
	fCellAmplitude(),
	fCellTime(),
	fProcessedEvents(),
	fhCellFlag(),
	fhCellWarm()
{
	fCurrentRunNumber = 254381;
	fPeriod           = "LHC16h";
	fPass             = "muon_caloLego";
	fTrigger          = "AnyINT";
	fWorkdir          = ".";
	fRunListFileName  = "runList.txt";
	fTrial            = 0;

	Init();
}

///
/// Constructor
///
//________________________________________________________________________
AliAnaCaloChannelAnalysis::AliAnaCaloChannelAnalysis(TString period, TString pass, TString trigger, Int_t runNumber,Int_t trial, TString workDir, TString listName):
	TObject(),
	fCurrentRunNumber(-1),
	fPeriod(),
	fPass(),
	fTrigger(),
	fNoOfCells(),
	fCellStartDCal(12288),
	fAnalysisOutput(),
	fAnalysisInput(),
	fRunList(),
	fRunListFileName(),
	fWorkdir(),
	fQADirect(),
	fMergedFileName(),
	fAnalysisVector(),
	fTrial(),
	fExternalFileName(),
	fTestRoutine(),
	fNMaxCols(48),
	fNMaxRows(24),
	fNMaxColsAbs(),
	fNMaxRowsAbs(),
	fFlag(),
	fCriterionCounter(),
	fWarmCell(),
	fCaloUtils(),
	fRootFile(),
	fCellAmplitude(),
	fCellTime(),
	fProcessedEvents(),
	fhCellFlag(),
	fhCellWarm()
{
	fCurrentRunNumber = runNumber;
	fPeriod           = period;
	fPass             = pass;    //ELI only for folder structure
	fTrigger          = trigger; //ELI important to select trigger in output file == different wagons in lego train
	fWorkdir          = workDir;
	fRunListFileName  = listName;
	fTrial            = trial;

	Init();
}

///
/// Initialize default parameters
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::Init()
{
	//......................................................
	//..Default values - can be set by functions
	fExternalFileName="";
	fTestRoutine=0;

	//..Settings for the input/output structure (hard coded)
	// TO BE CHANGED
	fAnalysisInput  =Form("AnalysisInput/%s",fPeriod.Data());
	fAnalysisOutput =Form("AnalysisOutput/%s/Version%i",fPeriod.Data(),fTrial);

	//..Make output directory if it doesn't exist
	//..first the period folder
	gSystem->mkdir(Form("%s/AnalysisOutput/%s",fWorkdir.Data(),fPeriod.Data()));
	//..then the version folder
	gSystem->mkdir(Form("%s/%s",fWorkdir.Data(),fAnalysisOutput.Data()));

	fMergedFileName= Form("%s/%s/%s/MergedRuns_%s.root",fWorkdir.Data(),fAnalysisInput.Data(),fPass.Data(),fTrigger.Data());
	fRunList       = Form("%s/%s/%s/%s",fWorkdir.Data(), fAnalysisInput.Data(), fPass.Data(), fRunListFileName.Data());
	fQADirect      = Form("CaloQA_%s",fTrigger.Data());

	TString fileName = Form("%s/%s/%s%s_Histograms_V%i.root",fWorkdir.Data(),fAnalysisOutput.Data(), fPass.Data(),fTrigger.Data() ,fTrial);
	fRootFile = new TFile(fileName,"recreate");
	//.. make sure the vector is empty
	fAnalysisVector.clear();

	//......................................................
	//..Initialize EMCal/DCal geometry
	fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(fCurrentRunNumber);
	fCaloUtils->AccessGeometry(aod);

    fNoOfCells    =fCaloUtils->GetEMCALGeometry()->GetNCells(); //..Very important number, never change after that point!
    Int_t nModules=fCaloUtils->GetEMCALGeometry()->GetNumberOfSuperModules();
    //..This is how the calorimeter looks like in the current period (defined by example run ID fCurrentRunNumber)
	cout<<"Called geometry for run number: "<<fCurrentRunNumber<<endl;
	cout<<"Number of supermod: "<<nModules<<endl;
	cout<<"Number of cells: "<<fNoOfCells<<endl;
	cout<<"Cell ID of first DCal cell: "<<fCellStartDCal<<endl;
	//cout<<"Number of supermod utils: "<<fCaloUtils->GetNumberOfSuperModulesUsed()<<endl; //..will always be 22 unless set by hand

	//......................................................
	//..Initialize flag array to store how the cell is categorized
	//..In the histogram: bin 1= cellID 0, bin 2= cellID 1 etc
	//..In the array: fFlag[cellID]= some information
	fFlag     = new Int_t[fNoOfCells];
	fWarmCell = new Bool_t[fNoOfCells];
	fFlag[fNoOfCells] = {0};      //..flagged as good by default
	fWarmCell[fNoOfCells] = {0};  //..flagged as not warm by default
	fCriterionCounter=2; //This value will be written in fflag and updates after each PeriodAnalysis
	//......................................................
	//..setings for the 2D histogram
	//fNMaxCols    = 48;  //eta direction
	//fNMaxRows    = 24;  //phi direction
	fNMaxColsAbs = 2*fNMaxCols;
	fNMaxRowsAbs = Int_t (nModules/2)*fNMaxRows; //multiply by number of supermodules

	//......................................................
	//..Create TLists to organize the outputfile
	fOutputListBad       = new TList();
	fOutputListGood      = new TList();
	fOutputListBadRatio  = new TList();
	fOutputListGoodRatio = new TList();

	fOutputListBad      ->SetName("BadCell_Amplitudes");
	fOutputListGood     ->SetName("GoodCell_Amplitudes");
	fOutputListBadRatio ->SetName("BadCell_AmplitudeRatios");
	fOutputListGoodRatio->SetName("GoodCell_AmplitudeRatios");

	//fOutputListGood    ->SetOwner();//ELI instead of delete in destructor??
	//fOutputListBadRatio    ->SetOwner();
	//fOutputListGoodRatio    ->SetOwner();

	//......................................................
	//..Create Histograms to store the flag in a root file
	fhCellFlag = new TH1F("fhCellFlag","fhCellFlag",fNoOfCells+10,0,fNoOfCells+10); //..cellID+1 = histogram bin
	fhCellWarm = new TH1F("fhCellWarm","fhCellWarm",fNoOfCells,0,fNoOfCells); //..cellID+1 = histogram bin
}

///
///	Main execution method.
///
/// 1)
/// If no external file is provided use MergeRuns() to merge historgrams from a runlist .txt file.
/// The merged outputfile contains 3 different histograms (hCellAmplitude, hCellTime and hNEventsProcessedPerRun).
/// 2)
/// Flags dead cells
/// 3)
/// Analyses merged histograms by calling BCAnalysis() and flags bad cells.
/// 4)
/// It calls SummarizeResults() to store all information in output files (.gif .txt .pdf)
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::Run()
{
//	cout<<"fired trigger class"<<AliAODEvent::GetFiredTriggerClasses()<<endl;

	if(fExternalFileName=="")
	{
		//..If no extrenal file is provided merge different runs together
		cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;
		cout<<". . .Start process by converting files. . . . . . . . . . . ."<<endl;
		cout<<endl;
		fMergedFileName = MergeRuns();
		if(fMergedFileName.IsNull())
		{
			Printf("File not produced, exit");
			return;
		}
		cout<<endl;
	}
	else
	{
		//..If extrenal file is provided load it
		cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;
		cout<<". . .Start process by loading external file. . . . . . . . . . ."<<endl;
		fMergedFileName= Form("%s/%s/%s/%s",fWorkdir.Data(),fAnalysisInput.Data(),fPass.Data(),fExternalFileName.Data());
	}
	cout<<". . .Load inputfile with name: "<<fMergedFileName<<" . . . . . . . ."<<endl;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..	Read all the needed input for the Bad/Dead Channel analysis
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TFile *mergedFileInput = new TFile(fMergedFileName);
	if(!mergedFileInput->IsOpen())
	{
		Printf("Error! Input file not found, abort");
		return;
	}
	fCellAmplitude   = (TH2F*) mergedFileInput->Get("hCellAmplitude");
	fCellTime        = (TH2F*) mergedFileInput->Get("hCellTime");
	fProcessedEvents = (TH1F*) mergedFileInput->Get("hNEventsProcessedPerRun");

	cout<<". . .Continue process by . . . . . . . . . . . ."<<endl;
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. DEAD CELLS
	//.. Flag dead cells with fFlag=1
	//.. this excludes cells from analysis (will not appear in results)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"o o o Flag dead cells o o o"<<endl;
	FlagAsDead();
	cout<<endl;
	cout<<endl;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. BAD CELLS
	//.. Flag dead cells with fFlag=2 and bigger
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"o o o Flag bad cells o o o"<<endl;
	BCAnalysis();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..In the end summarize results
	//..in a .pdf and a .txt file
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"o o o Write .txt for each period analyis with bad cells  o o o"<<endl;
	SummarizeResultsByFlag();

	cout<<"o o o Create summary documents for the entire analysis o o o"<<endl;
	SummarizeResults();
	fRootFile->WriteObject(fFlag,"FlagArray");
	fRootFile->Close();
	cout<<endl;

	//..make a reccomendation about the used energy range to be investigated
	//..and the binning
	TH1D *hRefDistr = BuildMeanFromGood();
	//..Find bin where reference has value 1, and the corresponding x-value
	Int_t binHeihgtOne            = hRefDistr->FindLastBinAbove(1);
	Double_t binCentreHeightOne   = hRefDistr->GetBinCenter(binHeihgtOne);
	cout<<". . .Recomendation:"<<endl;
	cout<<". . .With the current statistic on average a cell has 1 hit at "<<binCentreHeightOne<<" GeV"<<endl;
	cout<<". . .so it makes no sense to select energy ranges >"<<binCentreHeightOne<<" as cells will be maked bad just due to the lack of statistic"<<endl;
	cout<<". . .End of process . . . . . . . . . . . . . . . . . . . . ."<<endl;
	cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;
}

///
/// This function takes the QA outputs from several runs and merges the histograms
/// together, in case there were more than 100 events in this run.
/// The list of runs to be merged is given by fRunList
/// The output file contains three histograms hCellAmplitude, hCellTime and hNEventsProcessedPerRun
///
//________________________________________________________________________
TString AliAnaCaloChannelAnalysis::MergeRuns()
{
	cout<<"o o o Start conversion process o o o"<<endl;
	cout<<"o o o period: " << fPeriod << ", pass: " << fPass << ",  trigger: "<<fTrigger<< endl;

	//..Create histograms needed for adding all the files together
	//ELI a little problematic to hard code properties of histograms??
	TH2F *hCellAmplitude          = new TH2F("hCellAmplitude","Cell Amplitude",200,0,10,23040,0,23040);
	TH2F *hCellTime               = new TH2F("hCellTime","Cell Time",250,-275,975,23040,0,23040);
	TH1F *hNEventsProcessedPerRun = new TH1F("hNEventsProcessedPerRun","Number of processed events in analyzed runs",1,0,1);

	//..Open the text file with the run list numbers and run index
	cout<<"o o o Open .txt file with run indices. Name = " << fRunList << endl;
	FILE *pFile = fopen(fRunList.Data(), "r");
	if(!pFile)
	{
		cout<<"couldn't open file!"<<endl;
		return "";
	}
	Int_t nEntr;
	Int_t nEntrTot=0;
	Int_t q;
	Int_t ncols;
	Int_t nlines = 0 ;
	Int_t runId[500] ;
	while (1)
	{
		ncols = fscanf(pFile,"  %d ",&q);
		if (ncols< 0) break;
		runId[nlines]=q;
		nlines++;
	}
	fclose(pFile);


	//..Open the different .root files with help of the run numbers from the text file
	const Int_t nRun = nlines ;
	TString base;
	TString infile;
	TString singleRunFileName;

	cout<<"o o o Start merging process of " << nRun <<" files"<< endl;
	//..loop over the amount of run numbers found in the previous text file.
	for(Int_t i = 0 ; i < nRun ; i++)
	{
		base  = Form("%s/%s/%s/%d", fWorkdir.Data(), fAnalysisInput.Data(), fPass.Data(), runId[i]);
		/*if ((fPass=="cpass1_pass2")||(fPass=="cfPass1-2"))
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
		{*/
		//..This is a run2 case
		infile = Form("%s.root",base.Data()) ;
		//}

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
		TH1F *hNEvents;

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
		hNEvents =(TH1F *)outputList->FindObject("hNEvents");
		if(!hNEvents)
		{
			Printf("hNEvents not found");
			outputList->ls();
			continue;
		}
		nEntr =  (Int_t)hNEvents->GetEntries();

		//..does that mean do not merge small files?
		if (nEntr<100)
		{
			cout <<"    o File to small to be merged. Only N entries " << nEntr << endl;
			continue ;
		}
		cout <<"    o File with N entries " << nEntr<<" will be merged"<< endl;
		nEntrTot+=nEntr;
		hCellAmplitude->Add(hAmpId);
		hCellTime->Add(hTimeId);
		hNEventsProcessedPerRun->Add(hNEvents);

		//..Create copies of the original root files just with the bad channel QA
		//..So that a run by run bad channel analysis can be performed more easily
		singleRunFileName= Form("%s/%s/%s/%d_Filtered.root",fWorkdir.Data(),fAnalysisInput.Data(),fPass.Data(),runId[i]);
		TFile *singleRunFile = TFile::Open(singleRunFileName,"recreate");
		hNEventsProcessedPerRun->Write();
		hCellAmplitude->Write();
		hCellTime->Write();
		singleRunFile->Close();

		outputList->Delete();
		dir->Delete();
		f->Close();
		delete f;
	}

	//.. Save the merged histograms
	cout<<"o o o Save the merged histogramms to .root file with name: "<<fMergedFileName<<endl;
	cout<<"o o o "<<nEntrTot<<" events were merged"<<endl;
	TFile *BCF = TFile::Open(fMergedFileName,"recreate");
	hNEventsProcessedPerRun->Write();
	hCellAmplitude->Write();
	hCellTime->Write();
	BCF->Close();
	cout<<"o o o End conversion process o o o"<<endl;
	return fMergedFileName;
}


///
/// This function checks how many different criteria should be analysed.
/// It checks how many period analyses criteria are stored in the fAnalysisVector
/// and passes it to the period analyis for execution.
///
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::BCAnalysis()
{
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. BAD CELLS
	//.. Flag bad cells with fFlag= 2,3,4,5.. etc
	//.. this excludes cells from subsequent analysis
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TArrayD periodArray;
	for(Int_t i=0;i<(Int_t)fAnalysisVector.size();i++)
	{
		periodArray=fAnalysisVector.at(i);
		PeriodAnalysis(periodArray.At(0),periodArray.At(1),periodArray.At(2),periodArray.At(3));
		fCriterionCounter++;
		cout<<""<<endl;
		cout<<""<<endl;
	}
	cout<<"o o o End of bad channel analysis o o o"<<endl;
}

///
/// This function adds period analyses to the Bad Channel analysis.
/// Each period analysis needs to be specified with four parameters.
///
/// \param criteria -- selection criteria. Determines whether to check hits/event, energy/hit, time-, energy distr. etc
/// \param nsigma   -- n sigma cut. Select cells as bad outside this sigma range away from the mean of the distribution of cell values
/// \param emin     -- minimum energy. Perform this test only for a given energy range. This sets the minimum E
/// \param emax     -- maximum energy. Perform this test only for a given energy range. This sets the maximum E
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
///
/// This function builds the distribution of cell properties like eg. hit/event in a cell by calling BuildHitAndEnergyMean()
/// Then the cells are flagged as bad if the given cell property is far off (nsigma) of the mean distribution of all cell values.
///
/// There are different checks done depending on the given criteria:
/// 1 : deposited E/hit for E>Emin and E<Emax
/// 2 : recorded hits/event for E>Emin and E<Emax
/// 3 and 4 are currently not used - needs carefull checking
/// for the future checks on the time distribution can be implemented
///
/// \param crit -- criterium that distinguishs the type of distribution (1= E/hit, 2= hit/event)
/// \param nsigma -- range that defines good cells
/// \param emin -- min. energy for cell amplitudes
/// \param emax -- max. energy for cell amplitudes
///
//____________________________________________________________________
void AliAnaCaloChannelAnalysis::PeriodAnalysis(Int_t criterion, Double_t nsigma, Double_t emin, Double_t emax)
{
	//ELI criterion should be between 1-4
	cout<<"o o o o o o o o o o o o o o o o o o o o o o o o o"<<endl;
	cout<<"o o o PeriodAnalysis for flag "<<criterion<<" o o o"<<endl;
	cout<<"o o o Done in the energy range E "<<emin<<"-"<<emax<<endl;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. ANALYSIS OF CELLS WITH ENTRIES
	//.. Build average distributions and fit them
	//.. Three tests for bad cells:
	//.. 1) Average energy per hit
	//.. 2) Average hit per event
	//.. 3) ...
	//.. 4) ...
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TH1F* histogram;
	if(criterion < 6)cout<<"o o o Analyze average cell distributions o o o"<<endl;
	//..For case 1 or 2
	if(criterion < 3)   histogram = BuildHitAndEnergyMean(criterion, emin, emax,nsigma);
	//..For case 3, 4 or 5
	else if (criterion < 6) TestCellShapes(criterion, emin, emax, nsigma);

	Double_t range=0.001;
	if(emin>0.49)range=0.0005;
	if(emin>0.99)range=0.0001;
	if(emin>1.99)range=0.00005;

	if(criterion==1) FlagAsBad(criterion, histogram, nsigma, 200,-1);
	if(criterion==2) FlagAsBad(criterion, histogram, nsigma, 600,range);

	/*
	if(criterion==3)              FlagAsBad(criterion, histogram, nsigma, dnbins, maxval3);
	if(criterion==4)              FlagAsBad(criterion, histogram, nsigma, dnbins, maxval1);
	if(criterion==5)              FlagAsBad(criterion, histogram, nsigma, dnbins, maxval2);
	 */
}

///
/// Builds average hit per event and the average energy per hit is caluclated for each cell.
/// The output is a histogram with either of these two values as a function of cell ID.
///
/// \param crit -- criterium that distinguishs the type of distribution (1= E/hit, 2= hit/event)
/// \param emin -- min. energy for cell amplitudes
/// \param emax -- max. energy for cell amplitudes
/// \param nsigma -- range that defines good cells
///
//_________________________________________________________________________
TH1F* AliAnaCaloChannelAnalysis::BuildHitAndEnergyMean(Int_t crit, Double_t emin, Double_t emax, Double_t nsigma)
{
	cout<<"    o Calculate average cell hit per event and average cell energy per hit "<<endl;
	TH1F *histogram;
	if(crit==1)histogram = new TH1F(Form("hCellEtoNtotal_E%.2f-%.2f",emin,emax),Form("Energy per hit, %.2f < E < %.2f GeV",emin,emax), fNoOfCells,0,fNoOfCells);
	if(crit==2)histogram = new TH1F(Form("hCellNtotal_E%.2f-%.2f",emin,emax),Form("Number of hits per event, %.2f < E < %.2f GeV",emin,emax), fNoOfCells,0,fNoOfCells);
	histogram->SetXTitle("Abs. Cell Id");
	if(crit==1)histogram->SetYTitle("Energy per hit");
	if(crit==2)histogram->SetYTitle("Number of hits per event");
	histogram->GetXaxis()->SetNdivisions(510);
	Double_t totalevents = fProcessedEvents->Integral(1, fProcessedEvents->GetNbinsX());

	//..here the average hit per event and the average energy per hit is caluclated for each cell.
	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		Double_t Esum = 0;
		Double_t Nsum = 0;

		for (Int_t j = 1; j <= fCellAmplitude->GetNbinsX(); j++)
		{
			Double_t E = fCellAmplitude->GetXaxis()->GetBinCenter(j);
			Double_t N = fCellAmplitude->GetBinContent(j, cell+1);
			//..investigate only cells that were not flagged as dead ore bad
			if (E < emin || E > emax || fFlag[cell]!=0) continue;
			Esum += E*N;
			Nsum += N;
		}
		//..Set the values only for cells that are not yet marked as bad
		if(fFlag[cell]==0)
		{
			if(totalevents> 0. && crit==2)histogram->SetBinContent(cell+1, Nsum/totalevents);  //..number of hits per event
			if(Nsum > 0.       && crit==1)histogram->SetBinContent(cell+1, Esum/Nsum);         //..average energy per hit
		}
	}
	return histogram;
}
///
/// Empty function - to be implemented!
/// Possibility to add there a check on the cell time too, if the time is calibrated for the period
///
//_________________________________________________________________________
TH1F* AliAnaCaloChannelAnalysis::BuildTimeMean(Int_t crit, Double_t emin, Double_t emax, Double_t nsigma)
{
	TH1F *Histogram = new TH1F(Form("hSomethingWithTime_E%.2f-%.2f",emin,emax),Form("I don't know, %.2f < E < %.2f GeV",emin,emax), 2,0,2);

	return Histogram;
}
///
/// This method is currently not used
/// Test cells shape using fit function f(x)=A*exp(-B*x)/x^2.
/// Produce values per cell + distributions for A,B and chi2/ndf parameters.
///
/// \param crit -- criterium that distinguishs the type of distribution
/// \param fitemin --
/// \param fitemax --
/// \param nsigma --
///
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::TestCellShapes(Int_t crit, Double_t fitemin, Double_t fitemax, Double_t nsigma)
{
	//..binning parameters
	Int_t  ncells = fCellAmplitude->GetNbinsY();
	Double_t amin = fCellAmplitude->GetYaxis()->GetXmin();
	Double_t amax = fCellAmplitude->GetYaxis()->GetXmax();
	cout << "ncells " << ncells << " amin = " << amin << "amax = " << amax<< endl;

	//..initialize histograms
	TH1 *hFitA = new TH1F("hFitA_hCellAmplitude","Fit A value", ncells,amin,amax);
	hFitA->SetXTitle("AbsId");
	hFitA->SetYTitle("A");

	TH1 *hFitB = new TH1F("hFitB_hCellAmplitude","Fit B value", ncells,amin,amax);
	hFitB->SetXTitle("AbsIdhname");
	hFitB->SetYTitle("B");

	TH1 *hFitChi2Ndf = new TH1F("hFitChi2Ndf_hCellAmplitude","Fit #chi^{2}/ndf value", ncells,amin,amax);
	hFitChi2Ndf->SetXTitle("AbsId");
	hFitChi2Ndf->SetYTitle("#chi^{2}/ndf");

	Double_t maxval1=0., maxval2=0., maxval3=0.;
	Double_t prev=0., MSA=0., AvA = 0. ; //those param are used to automaticaly determined a reasonable maxval1
	Double_t prev2=0., MSB=0., AvB = 0.  ; //those param are used to automaticaly determined a reasonable maxval2
	Double_t prev3=0., MSki2=0., Avki2 = 0. ; //those param are used to automaticaly determined a reasonable maxval3
	Double_t ki2=0.0 ;
	for (Int_t k = 1; k <= fNoOfCells; k++)
	{
		TF1 *fit = new TF1("fit", "[0]*exp(-[1]*x)/x^2");
		TH1 *hCell = fCellAmplitude->ProjectionX("",k,k);
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

	//ELI something like thie in the future:
	//return histogram;
}


///
/// This function finds cells with zero entries
/// It flags them by setting the fFlag[CellID] to 1.
///
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::FlagAsDead()
{
	Int_t sumOfExcl=0;
	//..Direction of cell ID
	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		Double_t nSum = 0;
		//..Direction of amplitude (Checks energies from 0-10 GeV)
		for (Int_t amp = 1; amp <= fCellAmplitude->GetNbinsX(); amp++)
		{
			//..cellID+1 = histogram bin
			Double_t N = fCellAmplitude->GetBinContent(amp,cell+1);
			nSum += N;
		}
		//..If the amplitude in one cell is basically 0
		//..mark the cell as excluded
		if(nSum < 0.5 && fFlag[cell]==0)
		{
			//..Cell flagged as dead.
			//..Flag only if it hasn't been flagged before
			fFlag[cell] =1;
			sumOfExcl++;
		}
	}
	cout<<"    o Number of dead cells: "<<sumOfExcl<<endl;
	cout<<"     ("<<sumOfExcl<<")"<<endl;
}

///
///  This function flags bad cells
///  It uses an histogram set in BuildHitAndEnergyMean()
///  The value for each cell is collected and the abundancy
///  of this values is plotted in a summery histogram (distrib)
///
///  The summary histogram (distribution of cell values) is fitted with
///  a gaussion histogram. The good area is -+ nsigma.
///  Cells with values beyond that are flagged as bad
///
/// \param crit    -- flag with channel criteria to be filled
/// \param inhisto -- input histogram;
/// \param dnbins  -- number of bins in distribution;
/// \param dmaxval -- maximum value on distribution histogram.
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::FlagAsBad(Int_t crit, TH1F* inhisto, Double_t nsigma, Int_t dnbins, Double_t dmaxval)
{  
	gStyle->SetOptStat(0); // MG modif
	gStyle->SetOptFit(0);  // MG modif

	if(crit==1)cout<<"    o Fit average energy per hit distribution"<<endl;
	if(crit==2)cout<<"    o Fit average hit per event distribution"<<endl;

	Int_t cellColumn=0,cellRow=0;
	Int_t cellColumnAbs=0,cellRowAbs=0;
	Int_t trash;

	TString histoName=inhisto->GetName();
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
	Double_t totalevents = fProcessedEvents->Integral(1, fProcessedEvents->GetNbinsX());
    //..Be aware that the bin width should never be larger than 1/totalevents because this
	//..is the minimum differce between the cells. One hit more/event otherwise you will see ugly empty bins in the histogram
	if(((dmaxval-inhisto->GetMinimum())/(dnbins*1.0))<1.0/totalevents)
	{
	    dnbins=(dmaxval-inhisto->GetMinimum())*totalevents;
	    cout<<"Problem - Reset dnbins to new value:"<<dnbins<<endl;
	}

	//..build histos
	TH1F *distrib = new TH1F(Form("%sDistr",(const char*)histoName), "", dnbins, inhisto->GetMinimum(), dmaxval);
	distrib->SetXTitle(inhisto->GetYaxis()->GetTitle());
	distrib->SetYTitle("Entries");
	TH1F *distrib_wTRDStruc = new TH1F(Form("%sDistr_wTRD",(const char*)histoName), "", dnbins, inhisto->GetMinimum(), dmaxval);
	TH1F *distrib_woTRDStruc= new TH1F(Form("%sDistr_woTRD",(const char*)histoName), "", dnbins, inhisto->GetMinimum(), dmaxval);


	//..build two dimensional histogram with values row vs. column
	TH2F *plot2D = new TH2F(Form("%s_HitRowColumn",(const char*)histoName),Form("%s_HitRowColumn",(const char*)histoName),fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	plot2D->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D->GetYaxis()->SetTitle("cell row (#phi direction)");

	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		//..Do that only for cell ids also accepted by the code
		if(!fCaloUtils->GetEMCALGeometry()->CheckAbsCellId(cell))continue;

		//..Fill histograms only for cells that are not yet flagged as bad
		if(fFlag[cell]==0)
		{
			//..fill the distribution of avarge cell values
			distrib->Fill(inhisto->GetBinContent(cell+1));
			//if(cell<fCellStartDCal)distrib_wTRDStruc->Fill(inhisto->GetBinContent(cell+1));
			//else                   distrib_woTRDStruc ->Fill(inhisto->GetBinContent(cell+1));
			//..Get Row and Collumn for cell ID
			fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(cell,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
			if(cellColumnAbs> fNMaxColsAbs || cellRowAbs>fNMaxRowsAbs)
			{
				cout<<"Problem! wrong calculated number of max col and max rows"<<endl;
				cout<<"current col: "<<cellColumnAbs<<", max col"<<fNMaxColsAbs<<endl;
				cout<<"current row: "<<cellRowAbs<<", max row"<<fNMaxRowsAbs<<endl;
			}
			plot2D->SetBinContent(cellColumnAbs,cellRowAbs,inhisto->GetBinContent(cell+1));
			//..check TRD support structure
			if(IsCoveredByTRD(cellRowAbs,cellColumnAbs)==1)
			{
				distrib_wTRDStruc->Fill(inhisto->GetBinContent(cell+1));
			}
			else
			{
				distrib_woTRDStruc ->Fill(inhisto->GetBinContent(cell+1));
			}
		}
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .draw histogram + distribution
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	TCanvas *c1 = new TCanvas(histoName,histoName,900,900);
	c1->ToggleEventStatus();
	TPad*    upperPad    = new TPad("upperPad", "upperPad",.005, .5, .995, .995);
	TPad*    lowerPadLeft = new TPad("lowerPadL", "lowerPadL",.005, .005, .5, .5);
	TPad*    lowerPadRight = new TPad("lowerPadR", "lowerPadR",.5, .005, .995, .5);
	upperPad->Draw();
	lowerPadLeft->Draw();
	lowerPadRight->Draw();

	upperPad->cd();
	upperPad->SetLeftMargin(0.045);
	upperPad->SetRightMargin(0.05);
	upperPad->SetLogy();
	inhisto->SetTitleOffset(0.6,"Y");
	inhisto->GetXaxis()->SetRangeUser(0,fNoOfCells+1);

	inhisto->SetLineColor(kBlue+1);
	inhisto->Draw();

	lowerPadRight->cd();
	lowerPadRight->SetLeftMargin(0.09);
	lowerPadRight->SetRightMargin(0.1);
	plot2D->Draw("colz");

	lowerPadLeft->cd();
	lowerPadLeft->SetLeftMargin(0.09);
	lowerPadLeft->SetRightMargin(0.07);
	lowerPadLeft->SetLogy();
	distrib->SetLineColor(kBlue+1);
	distrib->Draw();
	distrib_wTRDStruc->SetLineColor(kGreen+1);
	distrib_wTRDStruc->DrawCopy("same");
	distrib_woTRDStruc->SetLineColor(kMagenta+1);
	distrib_woTRDStruc->DrawCopy("same");

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

	TF1 *fit2 = new TF1("fit2", "gaus",0,10);
	//..start the fit with a mean of the highest value
	fit2->SetParameter(1,higherbin);

	distrib->Fit(fit2, "0LQEM", "", goodmin, goodmax);
	Double_t sig, mean;// chi2ndf;
	// Marie midif to take into account very non gaussian distrig
	mean    = fit2->GetParameter(1);
	sig     = fit2->GetParameter(2);
	//chi2ndf = fit2->GetChisquare()/fit2->GetNDF();

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

	TLegend *leg = new TLegend(0.60,0.70,0.9,0.85);
	leg->AddEntry(lline,"Good region boundary","l");
	leg->AddEntry(distrib_wTRDStruc,"Covered by TRD","l");
	leg->AddEntry(distrib_woTRDStruc,"wo TRD structure","l");
	leg->SetBorderSize(0);
	leg->Draw("same");

	fit2->SetLineColor(kOrange-3);
	fit2->SetLineStyle(1);//7
	fit2->Draw("same");

	TLatex* text = 0x0;
	if(crit==1) text = new TLatex(0.12,0.85,Form("Good range: %.2f-%.2f",goodmin,goodmax));
	if(crit==2) text = new TLatex(0.12,0.85,Form("Good range: %.2f-%.2fx10^-5",goodmin*100000,goodmax*100000));
	text->SetTextSize(0.06);
	text->SetNDC();
	text->SetTextColor(1);
	text->Draw();


	upperPad->cd();
	TLine *uline = new TLine(0, goodmax,fNoOfCells,goodmax);
	uline->SetLineColor(kGreen+2);
	uline->SetLineStyle(7);
	uline->Draw();
	TLine *lowline = new TLine(0, goodmin,fNoOfCells,goodmin);
	lowline->SetLineColor(kGreen+2);
	lowline->SetLineStyle(7);
	lowline->Draw();
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .Save histogram
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	c1->Update();
	TString name   =Form("%s/%s/criteria-_%d.gif",fWorkdir.Data(), fAnalysisOutput.Data(), crit);
	if(crit==1)name=Form("%s/%s/AverageEperHit_%s.gif",fWorkdir.Data(), fAnalysisOutput.Data(), (const char*)histoName);
	if(crit==2)name=Form("%s/%s/AverageHitperEvent_%s.gif",fWorkdir.Data(), fAnalysisOutput.Data(), (const char*)histoName);
	c1->SaveAs(name);

	fRootFile->WriteObject(c1,c1->GetName());
	fRootFile->WriteObject(plot2D,plot2D->GetName());
	fRootFile->WriteObject(distrib,distrib->GetName());
	fRootFile->WriteObject(inhisto,inhisto->GetName());
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . . Mark the bad cells in the fFlag array
	//. . .(fCriterionCounter= bad because cell average value lower than min allowed)
	//. . .(fCriterionCounter= bad because cell average value higher than max allowed)
	//. . .(0 by default - good cell)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"    o Flag bad cells that are outside the good range "<<endl;
	for(Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		//..cell=0 and bin=1, cell=1 and bin=2
		//.. <= throws out zeros, might not be a dead cell but still have zero entries in a given energy range
		if (inhisto->GetBinContent(cell+1) <= goodmin && fFlag[cell]==0)//ELI
		{
			fFlag[cell]=fCriterionCounter;
		}
		if (inhisto->GetBinContent(cell+1) > goodmax && fFlag[cell]==0)
		{
			fFlag[cell]=fCriterionCounter;
		}
	}
	cout<<"    o o o o o o o o o o o o o o o o o o o o o o o"<<endl;
}




///
/// In this function the final status of the analysis is summarized for each flag/period analysis.
/// A .txt file with dead and bad channel IDs is created for each check
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::SummarizeResultsByFlag()
{
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. RESULTS
	//.. 1) Print the bad cells
	//..    and write the results to a file
	//..    for each added period analysis
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TArrayD periodArray;
	Double_t emin,emax;
	Int_t criterion;
	TString output;
	Int_t nb1=0;

	for(Int_t i=0;i<(Int_t)fAnalysisVector.size();i++)
	{
		periodArray=fAnalysisVector.at(i);
		criterion  =periodArray.At(0);
		emin       =periodArray.At(2);
		emax       =periodArray.At(3);

		//..Print the results on the screen and
		//..write the results in a file
		output.Form("%s/%s/Criterion%d_Emin-%.2f_Emax-%.2f.txt",fWorkdir.Data(), fAnalysisOutput.Data(), criterion,emin,emax);
		ofstream file(output, ios::out | ios::trunc);
		if(!file)
		{
			cout<<"#### Major Error. Check the textfile!"<<endl;
		}
		file<<"fFlag="<<i+2<<"means Criterion : "<<criterion<<", emin = "<<emin<<" GeV"<<", emax = "<<emax<<" GeV"<<endl;
		cout<<"    o Criterion : "<<criterion<<", emin = "<<emin<<" GeV"<<", emax = "<<emax<<" GeV"<<" (Method "<<i<<")"<<endl;

		nb1=0;
		for(Int_t cellID=0;cellID<fNoOfCells;cellID++)
		{
			if(fFlag[cellID]==(i+2))
			{
				nb1++;
				file<<cellID<<", ";
			}
		}
		file<<"Total number of bad cells with fFlag=="<<i+2<<endl;
		file<<"("<<nb1<<")"<<endl;
		file.close();
		cout<<"    o Total number of bad cells ("<<nb1<<")"<<endl;
		cout<<endl;
	}
}
///
/// In this function the final status of the analysis is summarized.
/// .txt file with dead and bad channel IDs.
/// .pdf file with ampltidues and amplitude ratios of bad cells
/// .gif files with a 2D map of bad, dead and good channels
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::SummarizeResults()
{
	Int_t cellID, nDCalCells = 0, nEMCalCells = 0;
	TString cellSummaryFile, deadPdfName, badPdfName, ratioOfBad,goodCells;
	TH2F* cellAmp_masked= (TH2F*)fCellAmplitude->Clone("cellAmp_masked");
	TH2F* cellTime_masked= (TH2F*)fCellTime->Clone("fCellTime");

	deadPdfName     = Form("%s/%s/%s%s_Dead_Amplitudes_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fPass.Data(), fTrigger.Data() ,fTrial);
	badPdfName      = Form("%s/%s/%s%s_Bad_Amplitudes_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fPass.Data(), fTrigger.Data() ,fTrial);
	cellSummaryFile = Form("%s/%s/%s%s_Bad_Amplitudes_V%i.txt",fWorkdir.Data(), fAnalysisOutput.Data(), fPass.Data(), fTrigger.Data() ,fTrial); ;
	ratioOfBad      = Form("%s/%s/%s%s_BCRatio_Amplitudes_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fPass.Data(), fTrigger.Data() ,fTrial);
	goodCells       = Form("%s/%s/%s%s_Good_Amplitudes_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fPass.Data(), fTrigger.Data() ,fTrial);

	cout<<"    o Final results o "<<endl;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Write the final results of dead and bad cells in a file and on screen
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	ofstream file(cellSummaryFile, ios::out | ios::trunc);
	if(file)
	{
		file<<"Dead cells : "<<endl;
		cout<<"    o Dead cells : "<<endl;
		nEMCalCells =0;
		nDCalCells  =0;
		for(cellID=0; cellID<fNoOfCells; cellID++)
		{
			if(cellID==0)
			{
				file<<"In EMCal: "<<endl;
				//cout<<"    o In EMCal : "<<endl;
			}
			if(cellID==fCellStartDCal)
			{
				file<<"\n"<<endl;
				file<<"In DCal: "<<endl;
				//cout<<endl;
				//cout<<"    o In DCal : "<<endl;
			}
			if(fFlag[cellID]==1)
			{
				//file<<cellID<<"\n" ;
				file<<cellID<<", ";
				//cout<<cellID<<"," ;
				if(cellID<fCellStartDCal)nEMCalCells++;
				else                     nDCalCells++;
			}
		}
		//cout<<endl;
		file<<"\n"<<endl;
		file<<"EMCal ("<<nEMCalCells<<" ="<<100*nEMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<nDCalCells<<" ="<<100*nDCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;
		cout<<"    o EMCal ("<<nEMCalCells<<" ="<<100*nEMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<nDCalCells<<" ="<<100*nDCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;

		file<<"Bad cells: "<<endl;
		cout<<"    o Bad cells: "<<endl;
		nEMCalCells =0;
		nDCalCells  =0;
		for(cellID=0;cellID<fNoOfCells;cellID++)
		{
			if(cellID==0)
			{
				file<<"In EMCal: "<<endl;
				//cout<<"    o In EMCal : "<<endl;
			}
			if(cellID==fCellStartDCal)
			{
				file<<"\n"<<endl;
				file<<"In DCal: "<<endl;
				//cout<<endl;
				//cout<<"    o In DCal : "<<endl;
			}
			if(fFlag[cellID]>1)
			{
				//file<<cellID<<"\n" ;
				file<<cellID<<", ";
				//cout<<cellID<<"," ;
				if(cellID<fCellStartDCal)nEMCalCells++;
				else                     nDCalCells++;
			}
		}
		//cout<<endl;
		file<<"\n"<<endl;
		file<<"EMCal ("<<nEMCalCells<<" ="<<100*nEMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<nDCalCells<<" ="<<100*nDCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;
		cout<<"    o EMCal ("<<nEMCalCells<<" ="<<100*nEMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<nDCalCells<<" ="<<100*nDCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;
	}
	file.close();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Save the flagged cells to .pdf files
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"    o Save the bad channel spectra to a .pdf file"<<endl;
	SaveBadCellsToPDF(1,badPdfName) ;
	SaveBadCellsToPDF(10,ratioOfBad) ; //..Special case
	if(fTestRoutine==1)SaveBadCellsToPDF(2,goodCells) ;   //..Special case all good cells to check, should all have a flag naming them *Candidate*


	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Create a masked version of the Amp vs. ID and Time vs. ID histograms
	//..And Fill the histograms with the flag information
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		//..Direction of amplitude (Checks energies from 0-10 GeV)
		for (Int_t amp = 1; amp <= fCellAmplitude->GetNbinsX(); amp++)
		{
			if(fFlag[cell]!=0)
			{
				//..cellID+1 = histogram bin
				cellAmp_masked->SetBinContent(amp,cell+1,0);
				cellTime_masked->SetBinContent(amp,cell+1,0);
			}
		}
		fhCellFlag->SetBinContent(cell+1,fFlag[cell]);
		fhCellWarm->SetBinContent(cell+1,fWarmCell[cell]);
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Plot the 2D distribution of cells by flag
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	PlotFlaggedCells2D(0);    //..all good cells
	PlotFlaggedCells2D(1);    //..all dead cells
	PlotFlaggedCells2D(2,fCriterionCounter);  //..all bad cells
	PlotFlaggedCells2D(0,0);  //..Special case - Warm cells


	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Plot some summary canvases
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TCanvas *c1 = new TCanvas("CellProp","summary of cell properties",1000,1000);
	c1->ToggleEventStatus();
	c1->Divide(2,2);
	c1->cd(1)->SetLogz();
	//lowerPadRight->SetLeftMargin(0.09);
	//lowerPadRight->SetRightMargin(0.06);
	fCellAmplitude->SetXTitle("Cell Energy [GeV]");
	fCellAmplitude->SetYTitle("Abs. Cell Id");
	fCellAmplitude->Draw("colz");
	c1->cd(2)->SetLogz();
	fCellTime->SetXTitle("Cell Time [ns]");
	fCellTime->SetYTitle("Abs. Cell Id");
	fCellTime->Draw("colz");
	c1->cd(3)->SetLogz();
	//lowerPadRight->SetLeftMargin(0.09);
	//lowerPadRight->SetRightMargin(0.06);
	cellAmp_masked->SetTitle("Masked Cell Amplitude");
	cellAmp_masked->SetXTitle("Cell Energy [GeV]");
	cellAmp_masked->SetYTitle("Abs. Cell Id");
	cellAmp_masked->Draw("colz");
	c1->cd(4)->SetLogz();
	cellTime_masked->SetTitle("Masked Cell Time");
	cellTime_masked->SetXTitle("Cell Time [ns]");
	cellTime_masked->SetYTitle("Abs. Cell Id");
	cellTime_masked->Draw("colz");
	c1->Update();

	TCanvas *c2 = new TCanvas("CellFlag","summary of cell flags",1200,800);
	c2->ToggleEventStatus();
	c2->Divide(1,2);
	c2->cd(1);
	fhCellFlag->SetTitle("cell flag");
	fhCellFlag->SetXTitle("Abs. Cell Id");
	fhCellFlag->SetYTitle("flag");
	fhCellFlag->DrawCopy("hist");
	c2->cd(2);
	fhCellWarm->SetTitle("Warm cells");
	fhCellWarm->SetXTitle("Abs. Cell Id");
	fhCellWarm->SetYTitle("warm=1");
	fhCellWarm->DrawCopy("hist");
	c2->Update();
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Add different histograms/canvases to the output root file
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TString name   =Form("%s/%s/CellProperties.gif", fWorkdir.Data(),fAnalysisOutput.Data());
	c1->SaveAs(name);
	fRootFile->WriteObject(c1,c1->GetName());
	fRootFile->WriteObject(c2,c2->GetName());
	fRootFile->WriteObject(fCellAmplitude,fCellAmplitude->GetName());
	fRootFile->WriteObject(fCellTime,fCellTime->GetName());
	fRootFile->WriteObject(fhCellFlag,fhCellFlag->GetName());
	fRootFile->WriteObject(fhCellWarm,fhCellWarm->GetName());
	//..Save all amplitudes to the root file
	SaveHistoToFile();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Save also the identified warm channels into a text file.
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	file.open(cellSummaryFile, ios::out | ios::app);
	if(file)
	{
		file<<"Warm cells : "<<endl;
		cout<<"    o Warm cells : "<<endl;
		nEMCalCells =0;
		nDCalCells  =0;
		for(cellID=0; cellID<fNoOfCells; cellID++)
		{
			if(cellID==0)
			{
				file<<"In EMCal: "<<endl;
			}
			if(cellID==fCellStartDCal)
			{
				file<<"\n"<<endl;
				file<<"In DCal: "<<endl;
			}
			if(fWarmCell[cellID]==1)
			{
				file<<cellID<<", ";
				if(cellID<fCellStartDCal)nEMCalCells++;
				else                     nDCalCells++;
			}
		}
		file<<"\n"<<endl;
		file<<"EMCal ("<<nEMCalCells<<" ="<<100*nEMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<nDCalCells<<" ="<<100*nDCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;
		cout<<"    o EMCal ("<<nEMCalCells<<" ="<<100*nEMCalCells/(1.0*fCellStartDCal)<<"%), DCal ("<<nDCalCells<<" ="<<100*nDCalCells/(1.0*fNoOfCells-fCellStartDCal)<<"%)"<<endl;
	}
	file.close();

	//cout<<"    o Results can be found in : "<<endl;
	//cout<<"    o "<<cellSummaryFile<<endl;
	//cout<<"    o "<<badPdfName<<endl;
	//cout<<"    o "<<badPdfName<<endl;
}


///
/// Allow to produce a .pdf file with 9 histograms per page
/// They contain the energy distribution of bad cells (blue)
/// and compare them to the mean of all good cells (gray).
/// Different options are possible. To be selected with \param version
/// version=0 ->Print dead cells
/// version=1 ->Print bad cells
/// version=2 ->Print ratio of good cells to mean of all good cells
/// version=10->Print the ratio of badcell distr. and mean good cell distr.
///
/// \param version  -- flag that selects how and which cells are plotted into the pdf
/// \param pdfName  -- name of the .pdf file
///
///
//________________________________________________________________________
void AliAnaCaloChannelAnalysis::SaveBadCellsToPDF(Int_t version, TString pdfName)
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetFillColor(kWhite);
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetPalette(1);

	char title[100];
	char name[100];

	TH1D *hRefDistr = BuildMeanFromGood();
	fRootFile->WriteObject(hRefDistr,hRefDistr->GetName());
	Int_t firstCanvas=0;
	Bool_t candidate;
	TLatex* text = new TLatex(0.2,0.8,"*Candidate*");
	text->SetTextSize(0.07);
	text->SetTextColor(kOrange-3);
	text->SetNDC();

	TLatex* text2 = new TLatex(0.2,0.8,"*Not a Candidate*");
	text2->SetTextSize(0.08);
	text2->SetTextColor(8);

	TLatex* textA = new TLatex(0.65,0.62,"*test*");
	textA->SetTextSize(0.04);
	textA->SetTextColor(1);
	textA->SetNDC();

	//..collect cells in an internal vector.
	//..when the vector is of size=9 or at the end of being filled
	//..plot the channels into a canvas
	std::vector<Int_t> channelVector;
	channelVector.clear();
	cout<<"Start printing into .pdf for version: "<<version<<endl;
	for(Int_t cell=0;cell<fNoOfCells;cell++)
	{
		if(fFlag[cell]==1 && version==0)channelVector.push_back(cell);
		if(fFlag[cell]>1  && version==1)channelVector.push_back(cell);
		if(fFlag[cell]==0 && version==2)channelVector.push_back(cell);
		if(fFlag[cell]>1  && version==10)channelVector.push_back(cell);

		if(cell%2000==0)cout<<"...cell No."<<cell<<endl;
		//..when 9 bad cells are collected or we are at the end of the list, fill the canvas
		if(channelVector.size()==9 || cell == fNoOfCells-1)
		{
			cout<<"."<<flush;

			TString internal_pdfName=pdfName;
			TCanvas *c1 = new TCanvas("badcells","badcells",1000,750);
			if(channelVector.size() > 6)        c1->Divide(3,3);
			else if (channelVector.size() > 3)  c1->Divide(3,2);
			else                                c1->Divide(3,1);

			TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
			for(Int_t i=0; i< (Int_t)channelVector.size() ; i++)
			{
				sprintf(name, "Cell %d",channelVector.at(i)) ;
				TH1 *hCell = fCellAmplitude->ProjectionX(name,channelVector.at(i)+1,channelVector.at(i)+1);
				sprintf(title,"Cell No: %d    Entries: %d",channelVector.at(i), (Int_t)hCell->GetEntries()) ;

				c1->cd(i%9 + 1);
				c1->cd(i%9 + 1)->SetLogy();
				if(i==0)
				{
					leg->AddEntry(hRefDistr,"mean of good","l");
					leg->AddEntry(hCell,"current","l");
				}
				//..Check the distribution whether it looks like a *Candidate* for a miscalibrated warm cell
				candidate = CheckDistribution(hCell,hRefDistr);
                if(candidate==1)fWarmCell[channelVector.at(i)]=candidate;
				if(version>1)//..These are ratio plots of energy distr. of cell and mean of all good cells
				{
					hCell->Divide(hRefDistr);
				}
				//.. save histograms to file
				if(version==1) fOutputListBad->Add(hCell);
				if(version==10)fOutputListBadRatio->Add(hCell);
				if(version==2) fOutputListGoodRatio->Add(hCell);

				hCell->SetLineColor(kBlue+1);
				hCell->GetXaxis()->SetTitle("E (GeV)");
				hCell->GetYaxis()->SetTitle("N Entries");
				hCell->GetXaxis()->SetRangeUser(0.,10.);
				hCell->SetLineWidth(1) ;
				hCell->SetTitle(title);
				hRefDistr->SetLineColor(kGray+2);
				hRefDistr->SetLineWidth(1);

				hCell->Draw("hist");

				if(version==1)hRefDistr->Draw("same") ;

				//..Mark the histogram that could be miscalibrated and labelled as warm
				if(candidate==1 && (version==1 || version==10))
				{
					gPad->SetFrameFillColor(kYellow-10);
					text->Draw();
				}
				if(version==1)
				{
					textA->SetTitle(Form("Excluded by No. %d",fFlag[channelVector.at(i)]));
					textA->Draw();
				}
				if(version==2  && candidate==0)
				{
					gPad->SetFrameFillColor(kYellow-10);
					text2->Draw(); //..Draw a marker in the histogram that was falsley missed as a good candidate
					leg->Draw();
				}
				if(version<2)leg->Draw();
			}

			if(channelVector.size()<9 || cell == fNoOfCells-1)
			{
				internal_pdfName +=")";
				c1->Print(internal_pdfName.Data());
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
			delete leg;
			channelVector.clear();
		}
	}
	delete hRefDistr;
	//..Add the subdirectories to the file
	if(version==1) fRootFile->WriteObject(fOutputListBad,fOutputListBad->GetName());
	if(version==10)fRootFile->WriteObject(fOutputListBadRatio,fOutputListBadRatio->GetName());
	if(version==2) fRootFile->WriteObject(fOutputListGoodRatio,fOutputListGoodRatio->GetName());

	cout<<endl;
}
////
//// Build the mean cell amplitude distribution of all good cells
////
//_________________________________________________________________________
TH1D* AliAnaCaloChannelAnalysis::BuildMeanFromGood()
{
	TH1D* hGoodAmp;
	TH1D* hgoodMean;
	Int_t NrGood=0;
	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
		if(fFlag[cell]!=0)continue;
		NrGood++;
		if(NrGood==1)hgoodMean = (TH1D*)fCellAmplitude->ProjectionX("hgoodMean",cell+1,cell+1);
		else
		{
			hGoodAmp = (TH1D*)fCellAmplitude->ProjectionX("hGoodCells",cell+1,cell+1);
			hgoodMean->Add(hGoodAmp);
		}
	}
	hgoodMean->Scale(1.0/NrGood);

	return hgoodMean;
}
///
/// This is an automatic check of the amplitude ratio of cell/(mean of good cells)
/// It should help identifying cells that are candidates for recalibration
/// By default all cells are candidates. These checks identify obviously odd looking
/// cells and remove the candidate status (candidate=0).
/// These cells might have spikes, cliffs very steep slopes etc.
///
/// \param ratio  -- histogram that should be checked if it is bad or warm
/// \param reference  -- good reference histogram
///
//_________________________________________________________________________
Bool_t AliAnaCaloChannelAnalysis::CheckDistribution(TH1* histogram, TH1* reference)
{
	TString Name = Form("%s_ratio",histogram->GetName());
	TH1* ratio=(TH1*)histogram->Clone(Name);
	ratio->Divide(reference);

	Double_t percentageOfLast=0.6;
	Double_t higherThanMean=5;
	Double_t highestRatio=1000;
	Double_t maxMean=10;
	Double_t minMean=0.1;
	Double_t cliffSize=2;       //height before cliff shouldn't be double as high than after cliff

	//..By default each cell is a candidate for recalibration
	Bool_t candidate=1;
	//..Find bin where reference has value 1, and the corresponding x-value
	Int_t binHeihgtOne            = reference->FindLastBinAbove(1);
	Double_t binCentreHeightOne   = reference->GetBinCenter(binHeihgtOne);
	Double_t thirdBinCentre      = reference->GetBinCenter(3);

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Check the histogram
	//..Different checks to see whether the
	//..cell is really bad. Set candidate to 0.

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..check end of spectrum, should be larger than "percentageOfLast"% of the end of the mean histogram
	if(ratio->FindLastBinAbove(0)<ratio->FindBin(binCentreHeightOne*percentageOfLast))
	{
		candidate=0;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Check maximum of ratio. Cell should not have "highestRatio" times more entries than reference in any bin
	//ELI check that crieteria carfully - seems to work but not shure about it
	ratio->GetXaxis()->SetRangeUser(thirdBinCentre,10);//..zoom in to find the  maximum between "not first 2 bins" - 10 GeV
	if(ratio->GetMaximum()>highestRatio)//
	{
		candidate=0;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..check whether the ratio is much larger than 1
	//..calculate the mean in the relevant energy range
	Double_t mean=0;
	for(Int_t i=2;i<binHeihgtOne;i++)
	{
		mean+=ratio->GetBinContent(i);
	}
	mean*=1.0/(binHeihgtOne-1);//..divide by number of bins
	if(mean>maxMean || mean<minMean)
	{
		candidate=0;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..check whether there are large spikes in the histogram
	//..compare bin values to mean of the ratio. If there is a bin value with
	//..content "higherThanMean" times lareger than mean it's losing it candidate status
	mean=0;
	//..Find the maximum in the mean range (0-binHeihgtOne)
	ratio->GetXaxis()->SetRangeUser(0,binCentreHeightOne);
	Double_t localMaxBin=ratio->GetMaximumBin();

	for(Int_t i=2;i<binHeihgtOne;i++)
	{
		//..Exclude 0 bins and exclude bins near the maximum
		if(ratio->GetBinContent(i)<=0)        continue;
		if(i>localMaxBin-3 && i<localMaxBin+3)continue;
		mean+=ratio->GetBinContent(i);
	}
	mean*=1.0/(binHeihgtOne-1);//..divide by number of bins
	ratio->GetXaxis()->SetRangeUser(thirdBinCentre,binCentreHeightOne);//..zoom in to find the  maximum between 0-BinOne
	//cout<<"mean: "<<mean<<", max: "<<ratio->GetMaximum()<<endl;
	if(ratio->GetMaximum()>mean*higherThanMean)
	{
		candidate=0;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    //..Check for a cliff at 4GeV, happens in some cases
	Double_t beforeCliff=0;
	Double_t afterCliff=0;
	Int_t binsBefore=0;
	Int_t binsAfter=0;
	Int_t cliffBin = ratio->FindBin(4);
	for(Int_t i=cliffBin-10;i<cliffBin+11;i++)
	{
		//..Exclude 0 bins
		if(ratio->GetBinContent(i)<=0)continue;
		if(i<=cliffBin) binsBefore++;
		if(i>cliffBin)  binsAfter++;
		if(i<=cliffBin) beforeCliff+=ratio->GetBinContent(i);
		if(i>cliffBin)   afterCliff+=ratio->GetBinContent(i);
		beforeCliff*=1.0/binsBefore;
		afterCliff*=1.0/binsAfter;
	}
	if((beforeCliff-afterCliff)>cliffSize*afterCliff)
	{
		if(beforeCliff!=0 && afterCliff!=0)candidate=0;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Employ peak finder
/*	if(candidate==1)
	{
		TSpectrum *spec = new TSpectrum(2,1);  //..Nr peaks, resol. 1=3sigma distance
		Int_t nfound = spec->Search(ratio,4.3,"nobackground new",0.70);
		//GetNPeaks();
		//cout<<"found N peaks: "<<nfound<<endl;
	}
*/
	return candidate;
}

///
/// This function checks wether the cell is covered by the TRD
/// support structure which causes the cell hit/event number
/// to be shifted to lower values
///
/// \param row -- absolute row No. of the cell
/// \param collumn -- absolute cullumn No. of the cell
///
//_________________________________________________________________________
Bool_t AliAnaCaloChannelAnalysis::IsCoveredByTRD(Int_t row, Int_t collumn)
{
	//..TRD support structure
	//..(determined by eye, could be improved, but is already very acurate):
	//..collumn 4,5,6,7,8   33,34,35,36     58,59,60   85,86,87,88,89
	//..row     1    (21),22,23,24   45,46,47,(48)     69,70,71,(72)  (92),93,94,95   117,118,(119)  127   149,150,151   (173),174,175,(176)    198,199,200
	Bool_t coveredByTRDSupportStruc=0;

	if((collumn>3 && collumn<9) || (collumn>32 && collumn<37) || (collumn>57 && collumn<61) || (collumn>84 && collumn<90) ||
	   (row==1) ||(row>20 && row<25) || (row>44 && row<49) || (row>68 && row<73) || (row>91 && row<96) ||
	   (row>116 && row<120)|| row==127 || (row>148 && row<152) || (row>172 && row<177) || (row>197 && row<201)
	)
	{
		coveredByTRDSupportStruc=1;
	}
	return coveredByTRDSupportStruc;
}
///
/// Plots a 2D map of flagged cells, dead, bad, good
/// depending on the selected value of fFlag[]
///
/// \param flag1 -- plot the cells that have fFlag[cell]==flag1
/// \param flag2 -- plot the cells that have fFlag[cell]==flag2
/// \param flag3 -- plot the cells that have fFlag[cell]==flag3
///
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::PlotFlaggedCells2D(Int_t flagBegin,Int_t flagEnd)
{
	//..build two dimensional histogram with values row vs. column
	TString histoName;
	histoName = Form("2DChannelMap_Flag%d",flagBegin);
	if(flagBegin==0 && flagEnd==0)histoName = Form("2DChannelMap_Flag100");

	TH2F *plot2D = new TH2F(histoName,histoName,fNMaxColsAbs+2,-1.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+2,-1.5,fNMaxRowsAbs+0.5);
	plot2D->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D->GetYaxis()->SetTitle("cell row (#phi direction)");

	Int_t cellColumn=0,cellRow=0;
	Int_t cellColumnAbs=0,cellRowAbs=0;
	Int_t trash;

	for (Int_t cell = 0; cell < fNoOfCells; cell++)
	{
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
		if(flagEnd==-1 && fFlag[cell]==flagBegin)                        plot2D->SetBinContent(cellColumnAbs,cellRowAbs,1);
		if(flagEnd!=0 && flagEnd!=-1 && fFlag[cell]>=flagBegin && fFlag[cell]<=flagEnd)plot2D->SetBinContent(cellColumnAbs,cellRowAbs,1);
		if(flagBegin==0 && flagEnd==0 && fWarmCell[cell]==1)             plot2D->SetBinContent(cellColumnAbs,cellRowAbs,1); //warm cells


	}
	TCanvas *c1 = new TCanvas(histoName,histoName,500,500);
	c1->ToggleEventStatus();
	c1->cd();
	//lowerPadRight->SetLeftMargin(0.09);
	//lowerPadRight->SetRightMargin(0.06);
	plot2D->Draw("colz");

	TLatex* text = 0x0;
	if(flagBegin==0 && flagEnd==-1) text = new TLatex(0.2,0.8,"Good Cells");
	if(flagBegin==1) text = new TLatex(0.2,0.8,"Dead Cells");
	if(flagBegin>1)  text = new TLatex(0.2,0.8,"Bad Cells");
	if(flagBegin==0 && flagEnd==0) text = new TLatex(0.2,0.8,"Warm Cells");
	text->SetTextSize(0.06);
	text->SetNDC();
	text->SetTextColor(1);
	text->Draw();

	c1->Update();
	TString name =Form("%s/%s/%s.gif", fWorkdir.Data(),fAnalysisOutput.Data(), histoName.Data());
	c1->SaveAs(name);

	fRootFile->WriteObject(plot2D,plot2D->GetName());

}
///
/// This function saves all good cells amplitudes to a root file
///
//_________________________________________________________________________
void AliAnaCaloChannelAnalysis::SaveHistoToFile()
{
	char name[100];
	for(Int_t cell=0;cell<fNoOfCells;cell++)
	{
		sprintf(name, "Cell %d",cell) ;
		TH1 *hCell = fCellAmplitude->ProjectionX(name,cell+1,cell+1);
		if(fFlag[cell]==0)fOutputListGood->Add(hCell);
	}
	fRootFile->WriteObject(fOutputListGood,fOutputListGood->GetName());
}
