
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
#include <TError.h>

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
#include <BadChannelAna.h>
#include <AliAODEvent.h>

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::flush;
using std::ios;

/// \cond CLASSIMP
ClassImp(BadChannelAna);
/// \endcond

///
/// Default constructor
///
//________________________________________________________________________
BadChannelAna::BadChannelAna():
TObject(),
	fCurrentRunNumber(-1),
	fPeriod(),
	fTrainNo(),
	fTrigger(),
	fNoOfCells(),
	fCellStartDCal(12288),
	fStartCell(0),
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
	fPrint(0),
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
	fTrainNo          = "muon_caloLego";
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
BadChannelAna::BadChannelAna(TString period, TString train, TString trigger, Int_t runNumber,Int_t trial, TString workDir, TString listName):
	TObject(),
	fCurrentRunNumber(-1),
	fPeriod(),
	fTrainNo(),
	fTrigger(),
	fNoOfCells(),
	fCellStartDCal(12288),
	fStartCell(0),
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
	fPrint(0),
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
	fTrainNo          = train;   //only for folder structure
	fTrigger          = trigger; //important to select trigger in output file == different wagons in lego train
	fWorkdir          = workDir;
	fRunListFileName  = listName;
	fTrial            = trial;

	Init();
}

///
/// Initialize default parameters
///
//________________________________________________________________________
void BadChannelAna::Init()
{
	gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //Does not work -
	//......................................................
	//..Default values - can be set by functions
	fExternalFileName="";
	fTestRoutine=0;

	//..Settings for the input/output structure (hard coded)
	// TO BE CHANGED
	fAnalysisInput  =Form("AnalysisInput/%s",fPeriod.Data());
	fAnalysisOutput =Form("AnalysisOutput/%s/%s/Version%i",fPeriod.Data(),fTrainNo.Data(),fTrial);

	//..Make output directory if it doesn't exist
	//..first the period folder
	gSystem->mkdir(Form("%s/AnalysisOutput/%s",fWorkdir.Data(),fPeriod.Data()));
	//..then the Train folder
	gSystem->mkdir(Form("%s/AnalysisOutput/%s/%s",fWorkdir.Data(),fPeriod.Data(),fTrainNo.Data()));
	//..then the version folder
	gSystem->mkdir(Form("%s/%s",fWorkdir.Data(),fAnalysisOutput.Data()));

	fMergedFileName= Form("%s/%s/%s/MergedRuns_%s.root",fWorkdir.Data(),fAnalysisInput.Data(),fTrainNo.Data(),fTrigger.Data());
	fRunList       = Form("%s/%s/%s/%s",fWorkdir.Data(), fAnalysisInput.Data(), fTrainNo.Data(), fRunListFileName.Data());
	fQADirect      = Form("CaloQA_%s",fTrigger.Data());

	TString fileName = Form("%s/%s/%s_%s_Histograms_V%i.root",fWorkdir.Data(),fAnalysisOutput.Data(),fPeriod.Data(),fTrigger.Data() ,fTrial);
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
	fhCellWarm = new TH1F("fhCellWarm","fhCellWarm",fNoOfCells+10,0,fNoOfCells+10); //..cellID+1 = histogram bin
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
void BadChannelAna::Run(Bool_t mergeOnly)
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
		fMergedFileName= Form("%s/%s/%s/%s",fWorkdir.Data(),fAnalysisInput.Data(),fTrainNo.Data(),fExternalFileName.Data());
	}
	//..if ==1 only produce filtered and merged files and do not perform a BC analysis
	if(mergeOnly==0)
	{
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
		fProcessedEvents = (TH1F*) mergedFileInput->Get("hNEvents");

		cout<<". . .Continue process by . . . . . . . . . . . ."<<endl;
		//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
		//.. DEAD CELLS
		//.. Flag dead cells with fFlag=1
		//.. this excludes cells from analysis (will not appear in results)
		//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
		cout<<"o o o Flag dead cells o o o"<<endl;
		FlagAsDead();
		if(fPrint==1)cout<<endl;
		if(fPrint==1)cout<<endl;

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
		if(fPrint==1)cout<<"o o o Write .txt for each period analyis with bad cells  o o o"<<endl;
		SummarizeResultsByFlag();

		if(fPrint==1)cout<<"o o o Create summary documents for the entire analysis o o o"<<endl;
		SummarizeResults();
		fRootFile->WriteObject(fFlag,"FlagArray");
		fRootFile->Close();
		cout<<endl;

		//..make a reccomendation about the used energy range to be investigated
		//..and the binning
		TH1D *hRefDistr = BuildMeanFromGood();
		Double_t totalevents = fProcessedEvents->Integral();
		//..Find bin where reference has value "totalevents" (means 1 hit/event), and the corresponding x-value
		Int_t binHeightOne            = hRefDistr->FindLastBinAbove(1.0/totalevents);
		Double_t binCentreHeightOne   = hRefDistr->GetBinCenter(binHeightOne);
		cout<<". . .Recomendation:"<<endl;
		cout<<". . .With the current statistic on average a cell has 1 hit at "<<binCentreHeightOne<<" GeV"<<endl;
		cout<<". . .so it makes no sense to select energy ranges >"<<binCentreHeightOne<<" as cells will be"<<endl;
		cout<<". . .marked bad just due to the lack of statistic"<<endl;
		cout<<". . .End of process . . . . . . . . . . . . . . . . . . . . ."<<endl;
		cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;
	}
}

///
/// This function takes the QA outputs from several runs and merges the histograms
/// together, in case there were more than 100 events in this run.
/// The list of runs to be merged is given by fRunList
/// The output file contains three histograms hCellAmplitude, hCellTime and hNEventsProcessedPerRun
///
//________________________________________________________________________
TString BadChannelAna::MergeRuns()
{
	cout<<"o o o Start conversion process o o o"<<endl;
	cout<<"o o o period: " << fPeriod << ", pass: " << fTrainNo << ",  trigger: "<<fTrigger<< endl;

	//..Create histograms needed for adding all the files together
	TH1F *hNEventsProcessedPerRun= new TH1F("hNEvents","Number of processed events in analyzed runs",1,0,1);
	TH2F *hCellAmplitude;
	TH2F *hCellTime;

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
	Int_t usedFiles=0;
	//..loop over the amount of run numbers found in the previous text file.
	for(Int_t i = 0 ; i < nRun ; i++)
	{
		base  = Form("%s/%s/%s/%d", fWorkdir.Data(), fAnalysisInput.Data(), fTrainNo.Data(), runId[i]);
		//..This is a run2 case
		infile = Form("%s.root",base.Data()) ;

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
        usedFiles++;
		TH2F *hAmpId;
		TH2F *hTimeId;
		TH1F *hNEvents;

		hAmpId =(TH2F*)outputList->FindObject("EMCAL_hAmpId");
		if(!hAmpId)
		{
			Printf("hAmpId not found");
			outputList->ls();
			continue;
		}
		else
		{
			hAmpId->SetName("hCellAmplitude");
			hAmpId->SetTitle("Cell Amplitude");
		}

		hTimeId =(TH2F*)outputList->FindObject("EMCAL_hTimeId");
		if(!hTimeId)
		{
			Printf("hTimeId not found");
			outputList->ls();
			continue;
		}
		else
		{
			hTimeId->SetName("hCellTime");
			hTimeId->SetTitle("Cell Time");
		}

		if(i==0)
		{
			//..copy the properties of the mother histogram for adding them all up
			hCellAmplitude=(TH2F*)hAmpId->Clone("DummyName1");
			hCellAmplitude->Reset();
			hCellAmplitude->SetDirectory(0); //otherwise histogram will dissapear when file is closed
			//..copy the properties of the mother histogram for adding them all up
			hCellTime=(TH2F*)hTimeId->Clone("DummyName2");
			hCellTime->Reset();
			hCellTime->SetDirectory(0); //otherwise histogram will dissapear when file is closed
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
		singleRunFileName= Form("%s/%s/%s/%d_%sFiltered.root",fWorkdir.Data(),fAnalysisInput.Data(),fTrainNo.Data(),runId[i],fTrigger.Data());
		TFile *singleRunFile = TFile::Open(singleRunFileName,"recreate");
		//..do not yet normalize by number of events
		//..due to binning issues in later histograms
		//..but rather do it at the very end of the analysis
		hAmpId ->Write();
		hTimeId->Write();
		hNEvents->Write();
		singleRunFile->Close();

		outputList->Delete();
		dir->Delete();
		f->Close();
		delete f;
	}
	//..avoid creating empty files
	if(usedFiles>0)
	{
		//.. Save the merged histograms
		cout<<"o o o Save the merged histogramms to .root file with name: "<<fMergedFileName<<endl;
		cout<<"o o o "<<nEntrTot<<" events were merged"<<endl;
		TFile *BCF = TFile::Open(fMergedFileName,"recreate");
		//hNEventsProcessedPerRun->SetName("hNEvents");
		//hNEventsProcessedPerRun->SetTitle("Number of processed events");
		hNEventsProcessedPerRun->Write();
		hCellAmplitude->SetName("hCellAmplitude");
		hCellAmplitude->SetTitle("Cell Amplitude");
		hCellAmplitude->Write();
		hCellTime->SetName("hCellTime");
		hCellTime->SetTitle("Cell Time");
		hCellTime->Write();
		BCF->Close();
		cout<<"o o o End conversion process o o o"<<endl;
		return fMergedFileName;
	}
	else
	{
		return "";
	}
}


///
/// This function checks how many different criteria should be analysed.
/// It checks how many period analyses criteria are stored in the fAnalysisVector
/// and passes it to the period analyis for execution.
///
//_________________________________________________________________________
void BadChannelAna::BCAnalysis()
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
		if(fPrint==1)cout<<endl;
		if(fPrint==1)cout<<endl;
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
void BadChannelAna::AddPeriodAnalysis(Int_t criteria, Double_t nsigma, Double_t emin, Double_t emax)
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
void BadChannelAna::PeriodAnalysis(Int_t criterion, Double_t nsigma, Double_t emin, Double_t emax)
{
	//.. criterion should be between 1-4
	if(fPrint==1)cout<<"o o o o o o o o o o o o o o o o o o o o o o o o o"<<endl;
	if(fPrint==1)cout<<"o o o PeriodAnalysis for flag "<<criterion<<" o o o"<<endl;
	if(fPrint==1 && criterion < 3)cout<<"o o o Done in the energy range E "<<emin<<" - "<<emax<<endl;
	if(fPrint==1 && criterion == 3)cout<<"o o o Done in the time range t "<<emin<<" - "<<emax<<endl;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. ANALYSIS OF CELLS WITH ENTRIES
	//.. Build average distributions and fit them
	//.. Three tests for bad cells:
	//.. 1) Average energy per hit
	//.. 2) Average hit per event
	//.. 3) Average max of cell time distribution
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TH1F* histogram;
	if(fPrint==1)cout<<"o o o Analyze average cell distributions o o o"<<endl;
	//..For case 1 or 2
	if(criterion < 3)   histogram = BuildHitAndEnergyMean(criterion, emin, emax);
	if(criterion == 3)  histogram = BuildTimeMean(criterion, emin, emax); //.. in case of crit 3 emin is tmin and emax is tmax

	if(criterion==1)
	{
		if(emin>2)FlagAsBad(criterion, histogram, nsigma, -1);//..do not apply a lower boundary
		else      FlagAsBad(criterion, histogram, nsigma, 400);
	}
	if(criterion==2)
	{
		if(emin>2)FlagAsBad(criterion, histogram, nsigma, -1);//..do not narrow the integration window
		else      FlagAsBad(criterion, histogram, nsigma, 601);
		//FlagAsBad(criterion, histogram, nsigma,601);// ELIANE ELIANE
		//FlagAsBad(criterion, histogram, nsigma, -1);//..do not narrow the integration window
	}
	if(criterion==3) FlagAsBad(criterion, histogram, nsigma, 602);
}

///
/// Builds average hit per event and the average energy per hit is caluclated for each cell.
/// The output is a histogram with either of these two values as a function of cell ID.
///
/// \param crit -- criterium that distinguishs the type of distribution (1= E/hit, 2= hit/event)
/// \param emin -- min. energy for cell amplitudes
/// \param emax -- max. energy for cell amplitudes
///
//_________________________________________________________________________
TH1F* BadChannelAna::BuildHitAndEnergyMean(Int_t crit, Double_t emin, Double_t emax)
{
	if(fPrint==1)cout<<"    o Calculate average cell hit per event and average cell energy per hit "<<endl;
	TH1F *histogram;
	if(crit==1)histogram = new TH1F(Form("hCellEtoN_E%.2f-%.2f",emin,emax),Form("Energy per hit, %.2f < E < %.2f GeV",emin,emax), fNoOfCells,0,fNoOfCells);
	if(crit==2)histogram = new TH1F(Form("hCellNtoEvt_E%.2f-%.2f",emin,emax),Form("Number of hits in cell, %.2f < E < %.2f GeV",emin,emax), fNoOfCells,0,fNoOfCells);
	histogram->SetXTitle("Abs. Cell Id");
	if(crit==1)histogram->SetYTitle("Energy per hit");
	if(crit==2)histogram->SetYTitle("Number of hits in cell");
	histogram->GetXaxis()->SetNdivisions(510);

	//..count Events in Energy Range
	TH1F* pojection  = (TH1F*)fCellAmplitude->ProjectionX("Intermediate");
	fnEventsInRange  = pojection->Integral(pojection->GetBinContent(emin),pojection->GetBinContent(emax));

	//..here the average hit per event and the average energy per hit is caluclated for each cell.
	for (Int_t cell = fStartCell; cell < fNoOfCells; cell++)
	{
		Double_t Esum = 0;
		Double_t Nsum = 0;

		for (Int_t j = 1; j <= fCellAmplitude->GetNbinsX(); j++)
		{
			//To Do: I think one should also take the different bin width into account
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
			if(crit==2)             histogram->SetBinContent(cell+1, Nsum);  //..number of hits (per event -> not anymore)
			if(Nsum > 0. && crit==1)histogram->SetBinContent(cell+1, Esum/(Nsum)); //..average energy per hit
		}
	}
	return histogram;
}
///
/// Empty function - to be implemented!
/// Possibility to add there a check on the cell time too, if the time is calibrated for the period
///
//_________________________________________________________________________
TH1F* BadChannelAna::BuildTimeMean(Int_t crit, Double_t tmin, Double_t tmax)
{
	if(fPrint==1)cout<<"    o Calculate maximum of cell time distribution "<<endl;
	TString name;
	TH1F *histogram;
	histogram = new TH1F(Form("timeMax_t%.2f-%.2f",tmin,tmax),Form("Maximum of time distr., %.2f < t < %.2f GeV",tmin,tmax), fNoOfCells,0,fNoOfCells);
	histogram->SetXTitle("Abs. Cell Id");
	histogram->SetYTitle("time max");
	histogram->GetXaxis()->SetNdivisions(510);

	//..Time information
	for (Int_t cell = fStartCell; cell < fNoOfCells; cell++)
	{
		Double_t timeMax = -100;
		//..search for the maximum bin between tmin and tmax
		name=Form("Cell %d",cell);
		//sprintf(name, "Cell %d",channelVector.at(i)) ;
		TH1 *hCell = fCellTime->ProjectionX(name,cell+1,cell+1);
		hCell->GetXaxis()->SetRangeUser(tmin,tmax);

		timeMax = hCell->GetXaxis()->GetBinCenter(hCell->GetMaximumBin());
		if(cell>55 &&cell<65)cout<<"Integral: "<<hCell->Integral()<<endl;
		if(hCell->Integral()<=0)timeMax=-40;
		//..Set the values only for cells that are not yet marked as bad
		if(fFlag[cell]==0)
		{
			histogram->SetBinContent(cell+1, timeMax);  //..number of hits per event
		}
		else
		{
			histogram->SetBinContent(cell+1, -50);  //..number of hits per event
		}

	}
	return histogram;
}

///
/// This function finds cells with zero entries
/// It flags them by setting the fFlag[CellID] to 1.
///
//_________________________________________________________________________
void BadChannelAna::FlagAsDead()
{
	Int_t sumOfExcl=0;
	Int_t manualMaskCounter=0;
	//..sort the elements in manual mask list
	std::sort (fManualMask.begin(), fManualMask.end());

	//..Direction of cell ID
	for (Int_t cell = fStartCell; cell < fNoOfCells; cell++)
	{
		Double_t nSum = 0;
		//..Direction of amplitude (Checks energies from 0-nBins GeV)
		for (Int_t amp = 1; amp <= fCellAmplitude->GetNbinsX(); amp++)
		{
			//..cellID+1 = histogram bin
			Double_t N = fCellAmplitude->GetBinContent(amp,cell+1);
			nSum += N;
		}
		//..If the amplitude in one cell is basically 0
		//..mark the cell as excluded
		if(nSum == 0 && fFlag[cell]==0)
		{
			//..Cell flagged as dead.
			//..Flag only if it hasn't been flagged before
			fFlag[cell] =1;
			sumOfExcl++;
		}
		//..add here the manual masking
		if(manualMaskCounter<(Int_t)fManualMask.size() && fManualMask.at(manualMaskCounter)==cell)
		{
			fFlag[cell] = 2;
			manualMaskCounter++;
		}
	}
	if(fPrint==1)cout<<"    o Number of dead cells: "<<sumOfExcl<<endl;
	if(fPrint==1)cout<<"     ("<<sumOfExcl<<")"<<endl;
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
/// \param dmaxVal -- maximum value on distribution histogram.
//_________________________________________________________________________
void BadChannelAna::FlagAsBad(Int_t crit, TH1F* inhisto, Double_t nsigma, Double_t dnbins)
{  
    gStyle->SetOptStat(0); //..Do not plot stat boxes
	gStyle->SetOptFit(0);  //
	if(fPrint==1 && crit==1)cout<<"    o Fit average energy per hit distribution"<<endl;
	if(fPrint==1 && crit==2)cout<<"    o Fit average hit per event distribution"<<endl;
	if(fPrint==1 && crit==3)cout<<"    o Fit average hit maximum distribution"<<endl;

	Int_t cellColumn=0,cellRow=0;
	Int_t cellColumnAbs=0,cellRowAbs=0;
	Int_t trash;

	TString histoName=inhisto->GetName();
	Double_t goodmax = 0. ;
	Double_t goodmin = 0. ;
	//..set a user range so that the min and max is only found in a certain range
	inhisto->GetXaxis()->SetRangeUser(fStartCell,fNoOfCells);//..
	Double_t dminVal = inhisto->GetMinimum(0);
	Double_t dmaxVal = inhisto->GetMaximum();
	Double_t inputBins=dnbins;
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .determine settings for the histograms (range and binning)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//cout<<"max value: "<<dmaxVal<<", min value: "<<dminVal<<endl;
	if(crit==2 && inputBins==-1)	dnbins=dmaxVal-dminVal;
	if(crit==1 && inputBins==-1)	dnbins=400;

	if(crit==2 && inputBins!=-1)
	{
		//..calculate and print the "median"
		Int_t numBins = inhisto->GetXaxis()->GetNbins();

		numBins-=fStartCell;
		Double_t *x = new Double_t[numBins];
		Double_t* y = new Double_t[numBins];
		for (int i = 0; i < numBins; i++)
		{
			x[i] = inhisto->GetBinContent(i+fStartCell);
			y[i] = 1; //each value with weight 1
			if(x[i]==0)y[i] = 0;
		}
		Double_t medianOfHisto = TMath::Median(numBins,x,y);

		//..if dmaxVal is too far away from medianOfHisto the histogram
		//..range will be too large -> reduce the range
		//cout<<"max value: "<<dmaxVal<<" median of histogram: "<<medianOfHisto<<endl;
		if(medianOfHisto*10<dmaxVal)
		{
			//cout<<"- - - median too far away from max range"<<endl;
			dmaxVal=medianOfHisto+0.2*(dmaxVal-medianOfHisto);  //..reduce the distance between max and mean drastically to cut out the outliers
		}
		dnbins=dmaxVal-dminVal;

		if(dmaxVal-dminVal>100)
		{
			if(dnbins>2000)dnbins=0.01*(dmaxVal-dminVal); //..maximum 5000 bins. changed to 3000 .. lets see..
			if(dnbins>2000)dnbins=0.001*(dmaxVal-dminVal);//..maximum 5000 bins.
			if(dnbins<100) dnbins=0.02*(dmaxVal-dminVal); //..minimum 100 bins.
		}
	}

	if(crit==3)
	{
		//..obtain the binwidth for
		Int_t binWidth = fCellTime->GetXaxis()->GetBinWidth(1);
		dnbins = fCellTime->GetXaxis()->GetNbins();
		dminVal= fCellTime->GetXaxis()->GetBinCenter(1)-(binWidth*1.0)/2.0;
		dmaxVal= fCellTime->GetXaxis()->GetBinCenter(dnbins)+(binWidth*1.0)/2.0;
		cout<<"set the new histo with following settings- bins: "<<dnbins<<", min val = "<<dminVal<<", max val:"<<dmaxVal<<endl;
	}
	//..build histos
	TH1F *distrib = new TH1F(Form("%sDistr",(const char*)histoName), "", dnbins, dminVal, dmaxVal);
	distrib->SetXTitle(inhisto->GetYaxis()->GetTitle());
	distrib->SetYTitle("Entries");
	TH1F *distrib_wTRDStruc = new TH1F(Form("%sDistr_wTRD",(const char*)histoName), "", dnbins, dminVal, dmaxVal);
	TH1F *distrib_woTRDStruc= new TH1F(Form("%sDistr_woTRD",(const char*)histoName), "", dnbins, dminVal, dmaxVal);

	//..build two dimensional histogram with values row vs. column
	TH2F *plot2D = new TH2F(Form("%s_HitRowColumn",(const char*)histoName),Form("%s_HitRowColumn",(const char*)histoName),fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	plot2D->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D->GetYaxis()->SetTitle("cell row (#phi direction)");

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .build the distribution of average values
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	for (Int_t cell = fStartCell; cell < fNoOfCells; cell++)
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
	TPad*    upperPad      = new TPad("upperPad", "upperPad",.005, .5, .995, .995);
	TPad*    lowerPadLeft  = new TPad("lowerPadL", "lowerPadL",.005, .005, .5, .5);
	TPad*    lowerPadRight = new TPad("lowerPadR", "lowerPadR",.5, .005, .995, .5);
	upperPad->Draw();
	lowerPadLeft->Draw();
	lowerPadRight->Draw();

	upperPad->cd();
	upperPad->SetLeftMargin(0.05);
	upperPad->SetRightMargin(0.05);
	if(crit<3)upperPad->SetLogy();
	inhisto->SetTitleOffset(0.6,"Y");
	inhisto->GetXaxis()->SetRangeUser(0,fNoOfCells+1);
	inhisto->GetYaxis()->SetTitleOffset(0.7);
	inhisto->SetLineColor(kBlue+1);
	inhisto->Draw();

	lowerPadRight->cd();
	lowerPadRight->SetLeftMargin(0.09);
	lowerPadRight->SetRightMargin(0.12);
	plot2D->GetYaxis()->SetTitleOffset(1.3);
	plot2D->Draw("colz");

	lowerPadLeft->cd();
	lowerPadLeft->SetLeftMargin(0.09);
	lowerPadLeft->SetRightMargin(0.07);
	lowerPadLeft->SetLogy();
	distrib->SetLineColor(kBlue+1);
	distrib->GetYaxis()->SetTitleOffset(1.3);
	distrib->Draw();
	distrib_wTRDStruc->SetLineColor(kGreen+1);
	distrib_wTRDStruc->DrawCopy("same");
	distrib_woTRDStruc->SetLineColor(kMagenta+1);
	distrib_woTRDStruc->DrawCopy("same");

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .fit histogram
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	//..to exclude first 2 bins from max finding
	distrib->GetXaxis()->SetRangeUser(distrib->GetBinCenter(2),distrib->GetBinCenter(dnbins));
	Double_t maxDistr  = distrib->GetMaximum();
	Int_t    maxBin    = distrib->GetMaximumBin();
	Double_t maxCenter = distrib->GetBinCenter(maxBin);

	//..good range is around the max value as long as the
	//..bin content is larger than 2 entries
	for(Int_t i = maxBin ; i<=dnbins ; i++)
	{
		if(distrib->GetBinContent(i)<2) break ;
		goodmax = distrib->GetBinCenter(i);
	}
	for(Int_t i = maxBin ; i>2 ; i--)
	{
		if(distrib->GetBinContent(i)<2) break ;
		goodmin = distrib->GetBinLowEdge(i);
	}
	//if(maxBin<2)goodmin = distrib->GetBinLowEdge(2);

	//..Define min/max range of the fit function
	Double_t   minFitRange=goodmin;
    Double_t   maxFitRange=goodmax;
	//if(crit==2)minFitRange=distrib->GetBinLowEdge(2);//..exclude first 2 bins
	//if(crit==2)maxFitRange=dmaxVal;
    if(crit==3)minFitRange=-20;
	if(crit==3)maxFitRange=20;

	//cout<<"max bin:    "<<maxBin<<", max value: "<<maxDistr<<endl;
	//cout<<"start mean: "<<maxCenter<<", mean range: 0-"<<dmaxVal<<endl;
	//cout<<"fit range:  "<<minFitRange<<" - "<<maxFitRange<<endl;
	TF1 *fit2 = new TF1("fit2", "gaus",minFitRange,maxFitRange);
	//..start the fit with a mean of the highest value
	fit2->SetParLimits(0,0,maxDistr); //..do not allow negative ampl values
	fit2->SetParameter(1,maxCenter);  //..set the start value for the mean
	fit2->SetParLimits(1,0,dmaxVal);  //..do not allow negative mean values

	//..ELI - TO DO: eventually fit with TRD and without TRD seperatley
	distrib->Fit(fit2, "0LQEM", "", minFitRange, maxFitRange);
	Double_t sig, mean;// chi2ndf;
	mean    = fit2->GetParameter(1);
	sig     = fit2->GetParameter(2);

	//..for case 1 and 2 lower than 0 is an unphysical value
	if(crit<3 && mean <0.) mean=0.;

	goodmin = mean - nsigma*sig ;
	goodmax = mean + nsigma*sig ;
	//..for case 1 and 2 lower than 0 is an unphysical value
	if(crit<3 && goodmin <0.) goodmin=0.;
	if(inputBins==-1)         goodmin=-1; //..this is a special case for the very last histogram 3-40 GeV

	if(fPrint==1)cout<<"    o Result of fit: "<<endl;
	if(fPrint==1)cout<<"    o  "<<endl;
	if(fPrint==1)cout<<"    o Mean: "<<mean <<" sigma: "<<sig<<endl;
	if(fPrint==1)cout<<"    o good range : "<<goodmin <<" - "<<goodmax<<endl;

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
	if(crit==2) text = new TLatex(0.12,0.85,Form("Good range: %.2f-%.2f",goodmin,goodmax));
	if(crit==3) text = new TLatex(0.12,0.85,Form("Good range: %.2f-%.2f",goodmin,goodmax));
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
	if(crit==3)name=Form("%s/%s/AverageTimeMax_%s.gif",fWorkdir.Data(), fAnalysisOutput.Data(), (const char*)histoName);
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
	if(fPrint==1)cout<<"    o Flag bad cells that are outside the good range "<<endl;
	for(Int_t cell = fStartCell; cell < fNoOfCells; cell++)
	{
		//..cell=0 and bin=1, cell=1 and bin=2
		//.. <= throws out zeros, might not be a dead cell but still have zero entries in a given energy range
		if(inhisto->GetBinContent(cell+1) <= goodmin && fFlag[cell]==0)
		{
			//cout<<"smaller than min value: "<<inhisto->GetBinContent(cell+1)<<endl;
			fFlag[cell]=fCriterionCounter;
		}
		if(inhisto->GetBinContent(cell+1) > goodmax && fFlag[cell]==0)
		{
			//cout<<"bigger than max value: "<<inhisto->GetBinContent(cell+1)<<endl;
			fFlag[cell]=fCriterionCounter;
		}
	}
	if(fPrint==1)cout<<"    o o o o o o o o o o o o o o o o o o o o o o o"<<endl;
}

///
/// In this function the final status of the analysis is summarized for each flag/period analysis.
/// A .txt file with dead and bad channel IDs is created for each check
///
//________________________________________________________________________
void BadChannelAna::SummarizeResultsByFlag()
{
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. RESULTS
	//.. 1) Print the bad cells
	//..    and write the results to a file
	//..    for each added period analysis
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TArrayD periodArray;
	Double_t emin,emax,sig;
	Int_t criterion;
	TString output;
	Int_t nb1=0;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Save the results in a tWiki format for the webpage (https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalQABadChannels2)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TString aliceTwikiTable = Form("%s/%s/%s_TwikiTable_V%i.txt",fWorkdir.Data(), fAnalysisOutput.Data(), fTrigger.Data() ,fTrial); ;
	ofstream file2(aliceTwikiTable, ios::out | ios::trunc);
	if(file2)
	{
		file2<<"|*Criterium* | *E_min !GeV* | *E_max !GeV* | *sigma* | *Excluded Cells* |"<<endl;
	}
	file2.close();

	for(Int_t i=0;i<(Int_t)fAnalysisVector.size();i++)
	{
		periodArray=fAnalysisVector.at(i);
		criterion  =periodArray.At(0);
		emin       =periodArray.At(2);
		emax       =periodArray.At(3);
		sig        =periodArray.At(1);

		//..Print the results on the screen and
		//..write the results in a file
		output.Form("%s/%s/Criterion%d_Emin-%.2f_Emax-%.2f.txt",fWorkdir.Data(), fAnalysisOutput.Data(), criterion,emin,emax);
		ofstream file(output, ios::out | ios::trunc);
		if(!file)
		{
			cout<<"#### Major Error. Check the textfile!"<<endl;
		}
		file<<"fFlag="<<i+2<<"means Criterion : "<<criterion<<", emin = "<<emin<<" GeV"<<", emax = "<<emax<<" GeV"<<endl;
		if(fPrint==1)cout<<"    o Criterion : "<<criterion<<", emin = "<<emin<<" GeV"<<", emax = "<<emax<<" GeV"<<" (Method "<<i<<")"<<endl;

		nb1=0;
		for(Int_t cellID=fStartCell ;cellID<fNoOfCells;cellID++)
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
		if(fPrint==1)cout<<"    o Total number of bad cells ("<<nb1<<")"<<endl;
		if(fPrint==1)cout<<endl;
		//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
		//..Save the results in a tWiki format for the webpage (https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalQABadChannels2)
		//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
		ofstream file2(aliceTwikiTable, ios::out | ios::app);
		if(file2)
		{
			file2<<"|     "<<criterion<<"      |    "<<emin<<"      |    "<<emax<<"      |    "<<sig<<"   |       "<<nb1<<"        |"<<endl;
		}
		file2.close();
	}
}
///
/// In this function the final status of the analysis is summarized.
/// .txt file with dead and bad channel IDs.
/// .pdf file with ampltidues and amplitude ratios of bad cells
/// .gif files with a 2D map of bad, dead and good channels
///
//________________________________________________________________________
void BadChannelAna::SummarizeResults()
{
	Int_t cellID, nDeadDCalCells = 0, nDeadEMCalCells = 0, nDCalCells = 0, nEMCalCells = 0;
	Double_t perDeadEMCal,perDeadDCal,perBadEMCal,perBadDCal,perWarmEMCal,perWarmDCal;
	TString aliceTwikiTable, cellSummaryFile, deadPdfName, badPdfName, ratioOfBad,goodCells,goodCellsRatio,cellProp;
	TH2F* cellAmp_masked = (TH2F*)fCellAmplitude->Clone("hcellAmp_masked");
	TH2F* cellTime_masked= (TH2F*)fCellTime->Clone("fCellTime");

	deadPdfName     = Form("%s/%s/%s_Dead_Ampl_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fTrigger.Data() ,fTrial);
	badPdfName      = Form("%s/%s/%s_Bad_Ampl_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fTrigger.Data() ,fTrial);
	ratioOfBad      = Form("%s/%s/%s_Bad_Ampl_Ratio_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fTrigger.Data() ,fTrial);
	goodCells       = Form("%s/%s/%s_Good_Ampl_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fTrigger.Data() ,fTrial);
	goodCellsRatio  = Form("%s/%s/%s_Good_Ampl_Ratio_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fTrigger.Data() ,fTrial);
	cellSummaryFile = Form("%s/%s/%s_%s_Bad_Ampl_V%i.txt",fWorkdir.Data(), fAnalysisOutput.Data(),fPeriod.Data(), fTrigger.Data() ,fTrial); ;
	aliceTwikiTable = Form("%s/%s/%s_TwikiTable_V%i.txt",fWorkdir.Data(), fAnalysisOutput.Data(), fTrigger.Data() ,fTrial); ;
	cellProp        = Form("%s/%s/%s_CellProp_V%i.pdf",fWorkdir.Data(), fAnalysisOutput.Data(), fTrigger.Data() ,fTrial);

	cout<<"    o Final results o "<<endl;
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Create a masked version of the Amp vs. ID and Time vs. ID histograms
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	for (Int_t cell = fStartCell; cell < fNoOfCells; cell++)
	{
		//..Direction of amplitude (Checks energies from 0-nBins GeV)
		for (Int_t amp = 1; amp <= fCellAmplitude->GetNbinsX(); amp++)
		{
			if(fFlag[cell]!=0)
			{
				//..cellID+1 = histogram bin
				cellAmp_masked->SetBinContent(amp,cell+1,0);
			}
		}
		//..Direction of time (Checks times from -275-975 GeV)
		for (Int_t time = 1; time <= fCellTime->GetNbinsX(); time++)
		{
			if(fFlag[cell]!=0)
			{
				//..cellID+1 = histogram bin
				cellTime_masked->SetBinContent(time,cell+1,0);
			}
		}
	}
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Plot some summary canvases
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TCanvas *c1 = new TCanvas("CellProp","I summary of cell properties",1000,1000);
	c1->ToggleEventStatus();
	c1->Divide(2,2);
	c1->cd(1)->SetLogz();
	//lowerPadRight->SetLeftMargin(0.09);
	//lowerPadRight->SetRightMargin(0.06);
	fCellAmplitude->SetXTitle("Cell Energy [GeV]");
	fCellAmplitude->SetYTitle("Abs. Cell Id");
	fCellAmplitude->GetYaxis()->SetTitleOffset(1.6);
	fCellAmplitude->Draw("colz");
	c1->cd(2)->SetLogz();
	fCellTime->SetXTitle("Cell Time [ns]");
	fCellTime->SetYTitle("Abs. Cell Id");
	fCellTime->GetYaxis()->SetTitleOffset(1.6);
	fCellTime->Draw("colz");
	c1->cd(3)->SetLogz();
	//lowerPadRight->SetLeftMargin(0.09);
	//lowerPadRight->SetRightMargin(0.06);
	cellAmp_masked->SetTitle("Masked Cell Amplitude");
	cellAmp_masked->SetXTitle("Cell Energy [GeV]");
	cellAmp_masked->SetYTitle("Abs. Cell Id");
	cellAmp_masked->GetYaxis()->SetTitleOffset(1.6);
	cellAmp_masked->Draw("colz");
	c1->cd(4)->SetLogz();
	cellTime_masked->SetTitle("Masked Cell Time");
	cellTime_masked->SetXTitle("Cell Time [ns]");
	cellTime_masked->SetYTitle("Abs. Cell Id");
	cellTime_masked->GetYaxis()->SetTitleOffset(1.6);
	cellTime_masked->Draw("colz");
	c1->Update();

	TCanvas *c1_ratio = new TCanvas("CellPropRatio","II summary of cell properties ratio",1000,500);
	c1_ratio->ToggleEventStatus();
	c1_ratio->Divide(2);
	c1_ratio->cd(1)->SetLogz();
	cellAmp_masked->SetTitle("Masked Cell Amplitude");
	cellAmp_masked->GetZaxis()->SetRangeUser(0.0001,10e7);
	cellAmp_masked->Draw("colz");
	c1_ratio->cd(2)->SetLogz();

	TH1D *hRefDistr  = BuildMeanFromGood();
	TH2F* ratio2DAmp =(TH2F*)cellAmp_masked->Clone("ratio2DAmp");
	TH2F* Sum2DIdeal =(TH2F*)cellAmp_masked->Clone("Sum2DIdeal");
	Sum2DIdeal->Reset();
    //..Create an ideal 2D energy distribution for division.
	//..Helps to identify whether there are some cells that still
	//..need to be masked by hand
	for(Int_t eBin=0;eBin<Sum2DIdeal->GetNbinsX();eBin++)
	{
		Double_t binVal=hRefDistr->GetBinContent(eBin+1);
		for(Int_t icell=0;icell<Sum2DIdeal->GetNbinsY();icell++)
		{
			Sum2DIdeal->SetBinContent(eBin+1,icell+1,binVal);
		}
	}
	ratio2DAmp->SetTitle("Ratio of cell Amplitude to mean cell ampl.");
	ratio2DAmp->Divide(Sum2DIdeal);
	ratio2DAmp->GetZaxis()->UnZoom();
	ratio2DAmp->Draw("colz");

	TCanvas *c1_proj = new TCanvas("CellPropPProj","III summary of cell properties",1000,500);
	c1_proj->ToggleEventStatus();
	c1_proj->Divide(2);
	c1_proj->cd(1)->SetLogy();
	TH1* projEnergyMask = cellAmp_masked->ProjectionX(Form("%sMask_Proj",cellAmp_masked->GetName()),fStartCell,fNoOfCells);
	projEnergyMask->SetXTitle("Cell Energy [GeV]");
	projEnergyMask->GetYaxis()->SetTitleOffset(1.6);
	projEnergyMask->SetLineColor(kTeal+3);
	projEnergyMask->DrawCopy(" hist");

	TH1* projEnergy = fCellAmplitude->ProjectionX(Form("%s_Proj",fCellAmplitude->GetName()),fStartCell,fNoOfCells);
	projEnergy->DrawCopy("same hist");

	c1_proj->cd(2)->SetLogy();
	TH1* projTimeMask = cellTime_masked->ProjectionX(Form("%s_Proj",cellTime_masked->GetName()),fStartCell,fNoOfCells);
	projTimeMask->SetXTitle("Cell Time [ns]");
	projTimeMask->GetYaxis()->SetTitleOffset(1.6);
	projTimeMask->SetLineColor(kGreen+3);
	projTimeMask->DrawCopy("hist");
	TH1* projTime = fCellTime->ProjectionX(Form("%s_Proj",fCellTime->GetName()),fStartCell,fNoOfCells);
	projTime->DrawCopy("same hist");
	c1_proj->Update();

	//..save to a PDF
	c1       ->Print(Form("%s(",cellProp.Data()));
	c1_ratio ->Print(Form("%s",cellProp.Data()));
	c1_proj  ->Print(Form("%s)",cellProp.Data()));
	//..Scale the histogtams by the number of events
	//..so that they are more comparable for a run-by-run
	//..analysis
	Double_t totalevents = fProcessedEvents->Integral();
	fCellAmplitude ->Scale(1.0/totalevents);
	cellAmp_masked ->Scale(1.0/totalevents);
	fCellTime      ->Scale(1.0/totalevents);

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Write the final results of dead and bad cells in a file and on screen
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	ofstream file(cellSummaryFile, ios::out | ios::trunc);
	if(file)
	{
		file<<"Dead cells : "<<endl;
		cout<<"    o Dead cells : "<<endl;
		for(cellID=fStartCell; cellID<fNoOfCells; cellID++)
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
			if(fFlag[cellID]==1)
			{
				file<<cellID<<", ";
				if(cellID<fCellStartDCal)nDeadEMCalCells++;
				else                     nDeadDCalCells++;
			}
		}
		file<<"\n"<<endl;
		perDeadEMCal=100*nDeadEMCalCells/(1.0*fCellStartDCal);
		perDeadDCal=100*nDeadDCalCells/(1.0*fNoOfCells-fCellStartDCal);
		file<<"EMCal ("<<nDeadEMCalCells<<" ="<<perDeadEMCal<<"%), DCal ("<<nDeadDCalCells<<" ="<<perDeadDCal<<"%)"<<endl;
		cout<<"    o EMCal ("<<nDeadEMCalCells<<" ="<<perDeadEMCal<<"%), DCal ("<<nDeadDCalCells<<" ="<<perDeadDCal<<"%)"<<endl;

		file<<"Bad cells: "<<endl;
		cout<<"    o Bad cells: "<<endl;
		for(cellID=fStartCell;cellID<fNoOfCells;cellID++)
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
			if(fFlag[cellID]>1)
			{
				file<<cellID<<", ";
				if(cellID<fCellStartDCal)nEMCalCells++;
				else                     nDCalCells++;
			}
		}
		file<<"\n"<<endl;
		perBadEMCal=100*nEMCalCells/(1.0*fCellStartDCal);
		perBadDCal =100*nDCalCells/(1.0*fNoOfCells-fCellStartDCal);
		file<<"EMCal ("<<nEMCalCells<<" ="<<perBadEMCal<<"%), DCal ("<<nDCalCells<<" ="<<perBadDCal<<"%)"<<endl;
		cout<<"    o EMCal ("<<nEMCalCells<<" ="<<perBadEMCal<<"%), DCal ("<<nDCalCells<<" ="<<perBadDCal<<"%)"<<endl;
	}
	file.close();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Determine the number of warm cells and save the flagged cells to .pdf files
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	if(fPrint==1)cout<<"    o Save the bad channel spectra to a .pdf file"<<endl;
	SaveBadCellsToPDF(1,badPdfName) ;
	SaveBadCellsToPDF(10,ratioOfBad) ; //..Special case
	if(fTestRoutine==1)SaveBadCellsToPDF(2,goodCells) ;   //..Special case all good cells to check, should all have a flag naming them *Candidate*
	if(fTestRoutine==1)SaveBadCellsToPDF(20,goodCellsRatio) ;   //..Special case all good cells to check, should all have a flag naming them *Candidate*

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Fill the histograms with the flag information
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	for (Int_t cell = fStartCell; cell < fNoOfCells; cell++)
	{
		fhCellFlag->SetBinContent(cell+1,fFlag[cell]);
		fhCellWarm->SetBinContent(cell+1,fWarmCell[cell]);
	}
	TCanvas *c2 = new TCanvas("CellFlag","summary of cell flags",1200,800);
	c2->ToggleEventStatus();
	c2->Divide(1,2);
	c2->cd(1);
	fhCellFlag->SetTitle("cell flag");
	fhCellFlag->SetXTitle("Abs. Cell Id");
	fhCellFlag->SetYTitle("flag by which cell was excluded");
	fhCellFlag->DrawCopy("hist");
	c2->cd(2);
	fhCellWarm->SetTitle("Warm cells");
	fhCellWarm->SetXTitle("Abs. Cell Id");
	fhCellWarm->SetYTitle("warm=1");
	fhCellWarm->DrawCopy("hist");
	c2->Update();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Plot the 2D distribution of cells by flag
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	PlotFlaggedCells2D(0);    //..all good cells
	PlotFlaggedCells2D(1);    //..all dead cells
	PlotFlaggedCells2D(2,fCriterionCounter);  //..all bad cells
	PlotFlaggedCells2D(0,0);  //..Special case - Warm cells

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Add different histograms/canvases to the output root file
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TString name1,name2,name3;
	name1   = Form("%s/%s/CellPropertiesRatio.gif", fWorkdir.Data(),fAnalysisOutput.Data());
	c1_ratio->SaveAs(name1);
	name2   = Form("%s/%s/CellProperties.gif", fWorkdir.Data(),fAnalysisOutput.Data());
	c1->SaveAs(name2);
	name3   = Form("%s/%s/CellPropertiesProj.gif", fWorkdir.Data(),fAnalysisOutput.Data());
	c1_proj->SaveAs(name3);

	fRootFile->WriteObject(c1,c1->GetName());
	fRootFile->WriteObject(c1_ratio,c1_ratio->GetName());
	fRootFile->WriteObject(c1_proj,c1_proj->GetName());
	fRootFile->WriteObject(c2,c2->GetName());
	fRootFile->WriteObject(fCellAmplitude,fCellAmplitude->GetName());
	fRootFile->WriteObject(cellAmp_masked,cellAmp_masked->GetName());
	fRootFile->WriteObject(ratio2DAmp,ratio2DAmp->GetName());
	fRootFile->WriteObject(fCellTime,fCellTime->GetName());
	fRootFile->WriteObject(fProcessedEvents,fProcessedEvents->GetName());
	fRootFile->WriteObject(fhCellFlag,fhCellFlag->GetName());
	fRootFile->WriteObject(fhCellWarm,fhCellWarm->GetName());
	fRootFile->WriteObject(projEnergyMask,projEnergyMask->GetName());
	fRootFile->WriteObject(projEnergy,projEnergy->GetName());
	//..Save all amplitudes to the root file
	SaveHistoToFile();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Save also the identified warm channels into a text file.
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	Int_t nWEMCalCells =0;
	Int_t nWDCalCells  =0;
	file.open(cellSummaryFile, ios::out | ios::app);
	if(file)
	{
		file<<"Warm cells : "<<endl;
		if(fPrint==1)cout<<"    o Warm cells : "<<endl;
		for(cellID=fStartCell; cellID<fNoOfCells; cellID++)
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
				if(cellID<fCellStartDCal)nWEMCalCells++;
				else                     nWDCalCells++;
			}
		}
		file<<"\n"<<endl;
		perWarmEMCal= 100*nWEMCalCells/(1.0*fCellStartDCal);
		perWarmDCal = 100*nWDCalCells/(1.0*fNoOfCells-fCellStartDCal);
		file<<"EMCal ("<<nWEMCalCells<<" ="<<perWarmEMCal<<"%), DCal ("<<nWDCalCells<<" ="<<perWarmDCal<<"%)"<<endl;
		if(fPrint==1)cout<<"    o EMCal ("<<nWEMCalCells<<" ="<<perWarmEMCal<<"%), DCal ("<<nWDCalCells<<" ="<<perWarmDCal<<"%)"<<endl;
	}
	file.close();
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Save the results in a tWiki format for the webpage (https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalQABadChannels2)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	ofstream file2(aliceTwikiTable, ios::out | ios::app);
	if(file2)
	{
		file2<<"1=energy/hit, 2= hit/event"<<endl;
		file2<<"\n"<<endl;
		file2<<"| *Detector* |    *No of cells*    |  *percentage*  |"<<endl;
		file2<<"| Dead EMCal |    "<<nDeadEMCalCells<<"    |  "<<perDeadEMCal<<"%  |"<<endl;
		file2<<"| Bad EMCal |    "<<nEMCalCells<<"    |  "<<perBadEMCal<<"%  |"<<endl;
		file2<<"| - Warm EMCal |    "<<nWEMCalCells<<"    |  "<<perWarmEMCal<<"%  |"<<endl;
		file2<<"| Dead DCal |    "<<nDeadDCalCells<<"    |  "<<perDeadDCal<<"%  |"<<endl;
		file2<<"| Bad DCal |    "<<nDCalCells<<"    |  "<<perBadDCal<<"%  |"<<endl;
		file2<<"|  - Warm DCal |    "<<nWDCalCells<<"    |  "<<perWarmDCal<<"%  |"<<endl;
		file2<<"\n"<<endl;
	}
	file2.close();

}

///
/// Allow to produce a .pdf file with 9 histograms per page
/// They contain the energy distribution of bad cells (blue)
/// and compare them to the mean of all good cells (gray).
/// Different options are possible. To be selected with \param version
/// version=0 ->Print dead cells
/// version=1 ->Print bad cells
/// version=2 ->Print ratio of good cells to mean of all good cells
/// version=10->Print the ratio of BadCell distr. and mean good cell distr.
/// version=20->Print the ratio of GoodCell distr. and mean good cell distr.
///
/// \param version  -- flag that selects how and which cells are plotted into the pdf
/// \param pdfName  -- name of the .pdf file
///
///
//________________________________________________________________________
void BadChannelAna::SaveBadCellsToPDF(Int_t version, TString pdfName)
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
	text2->SetNDC();

	TLatex* textA = new TLatex(0.65,0.62,"*test*");
	textA->SetTextSize(0.04);
	textA->SetTextColor(1);
	textA->SetNDC();

	//..collect cells in an internal vector.
	//..when the vector is of size=9 or at the end of being filled
	//..plot the channels into a canvas
	std::vector<Int_t> channelVector;
	channelVector.clear();
	cout<<"    o Start printing into .pdf for version: "<<version<<endl;
	for(Int_t cell=fStartCell;cell<fNoOfCells;cell++)
	{
		if(fFlag[cell]==1 && version==0)channelVector.push_back(cell);
		if(fFlag[cell]>1  && version==1)channelVector.push_back(cell);
		if(fFlag[cell]==0 && version==2)channelVector.push_back(cell);
		if(fFlag[cell]>1  && version==10)channelVector.push_back(cell);
		if(fFlag[cell]==0 && version==20)channelVector.push_back(cell);

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
				if(version>2)//..These are ratio plots of energy distr. of cell and mean of all good cells
				{
					hCell->Divide(hRefDistr);
				}
				//.. save histograms to file
				if(version==1) fOutputListBad->Add(hCell);
				if(version==10)fOutputListBadRatio->Add(hCell);
				if(version==2) fOutputListGood->Add(hCell);
				if(version==20)fOutputListGoodRatio->Add(hCell);

				hCell->SetLineColor(kBlue+1);
				hCell->GetXaxis()->SetTitle("E (GeV)");
				hCell->GetYaxis()->SetTitle("N Entries");
				hCell->GetXaxis()->SetRangeUser(0.,10.);
				hCell->SetLineWidth(1) ;
				hCell->SetTitle(title);
				hRefDistr->SetLineColor(kGray+2);
				hRefDistr->SetLineWidth(1);

				hCell->Draw("hist");

				if(version==1 || version==2)hRefDistr->Draw("same") ;

				//..Mark the histogram that could be miscalibrated and labelled as warm
				if(candidate==1 && (version==1 || version==10))
				{
					gPad->SetFrameFillColor(kYellow-10);
					text->Draw();
				}
				if(version==1)
				{
					textA->SetTitle(Form("Excluded by No. %d",fFlag[channelVector.at(i)]));
					textA->DrawLatex(0.65,0.62,Form("Excluded by No. %d",fFlag[channelVector.at(i)]));
				}
				if(candidate==0 && (version==2 || version==20))
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
	cout<<endl;
	delete hRefDistr;
	//..Add the subdirectories to the file
	if(version==1) fRootFile->WriteObject(fOutputListBad,fOutputListBad->GetName());
	if(version==10)fRootFile->WriteObject(fOutputListBadRatio,fOutputListBadRatio->GetName());
	if(version==2) fRootFile->WriteObject(fOutputListGoodRatio,fOutputListGoodRatio->GetName());
	if(version==20)fRootFile->WriteObject(fOutputListGood,fOutputListGood->GetName());

	if(fPrint==1)cout<<endl;
}
////
//// Build the mean cell amplitude distribution of all good cells
////
//_________________________________________________________________________
TH1D* BadChannelAna::BuildMeanFromGood(Int_t warmIn)
{
	TH1D* hGoodAmp;
	TH1D* hgoodMean = (TH1D*)fCellAmplitude->ProjectionX("hgoodMean");
	hgoodMean->Reset();
	Int_t NrGood=0;

	for (Int_t cell = fStartCell; cell < fNoOfCells; cell++)
	{
		if(warmIn==0 && fFlag[cell]!=0 )continue;
		if(warmIn==1 && fFlag[cell]!=0 && fWarmCell[cell]==0)continue;
		if(warmIn==2 && fWarmCell[cell]==0)continue;
		NrGood++;

		hGoodAmp = (TH1D*)fCellAmplitude->ProjectionX("hGoodCells",cell+1,cell+1);
		hgoodMean->Add(hGoodAmp);
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
Bool_t BadChannelAna::CheckDistribution(TH1* histogram, TH1* reference)
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
	Double_t totalevents        = fProcessedEvents->Integral();
	Int_t binHeightOne          = reference->FindLastBinAbove(1.0/totalevents);//..look at the spectrum until there is 1hit/event GeV
	Double_t binCentreHeightOne = reference->GetBinCenter(binHeightOne);
	Double_t thirdBinCentre     = reference->GetBinCenter(3);

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Check the histogram
	//..Different checks to see whether the
	//..cell is really bad. Set candidate to 0.

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..check end of spectrum, should be larger than "percentageOfLast"% of the end of the mean histogram
	if(ratio->FindLastBinAbove(0)<ratio->FindBin(binCentreHeightOne*percentageOfLast))
	{
		candidate=0;
		//cout<<"1"<<endl;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Check maximum of ratio. Cell should not have "highestRatio" times more entries than reference in any bin
	//ELI - TO DO: check that crieteria carfully - seems to work but not shure about it
	ratio->GetXaxis()->SetRangeUser(thirdBinCentre,10);//..zoom in to find the  maximum between "not first 2 bins" - 10 GeV
	if(ratio->GetMaximum()>highestRatio)//
	{
		candidate=0;
		//cout<<"2"<<endl;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..check whether the ratio is much larger than 1
	//..calculate the mean in the relevant energy range
	Double_t mean=0;
	Int_t nullEntries=0;
	for(Int_t i=2;i<binHeightOne;i++)
	{
		if(ratio->GetBinContent(i)!=0)mean+=ratio->GetBinContent(i);
		else nullEntries++;
	}
	mean*=1.0/(binHeightOne-1-nullEntries);//..divide by number of bins (excluding bins without entries)
	if(mean>maxMean || mean<minMean)
	{
		candidate=0;
		//cout<<"3"<<endl;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..check whether there are large spikes in the histogram
	//..compare bin values to mean of the ratio. If there is a bin value with
	//..content "higherThanMean" times lareger than mean it's losing it candidate status
	mean=0;
	//..Find the maximum in the mean range (0-binHeightOne)
	ratio->GetXaxis()->SetRangeUser(0,binCentreHeightOne);
	Double_t localMaxBin=ratio->GetMaximumBin();

	for(Int_t i=2;i<binHeightOne;i++)
	{
		//..Exclude 0 bins and exclude bins near the maximum
		if(ratio->GetBinContent(i)<=0)        continue;
		if(i>localMaxBin-3 && i<localMaxBin+3)continue;
		mean+=ratio->GetBinContent(i);
	}
	mean*=1.0/(binHeightOne-1);//..divide by number of bins
	ratio->GetXaxis()->SetRangeUser(thirdBinCentre,binCentreHeightOne);//..zoom in to find the  maximum between 0-BinOne
	//cout<<"mean: "<<mean<<", max: "<<ratio->GetMaximum()<<endl;
	if(ratio->GetMaximum()>mean*higherThanMean)
	{
		candidate=0;
		//cout<<"4"<<endl;
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
		if(beforeCliff!=0 && afterCliff!=0)cout<<"5"<<endl;
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
Bool_t BadChannelAna::IsCoveredByTRD(Int_t row, Int_t collumn)
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
void BadChannelAna::PlotFlaggedCells2D(Int_t flagBegin,Int_t flagEnd)
{
	//..build two dimensional histogram with values row vs. column
	TString histoName;
	histoName = Form("2DChannelMap_Flag%d",flagBegin);
	if(flagBegin==0 && flagEnd==0)histoName = Form("2DChannelMap_Flag100");

	TH2F *plot2D = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	plot2D->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D->GetYaxis()->SetTitle("cell row (#phi direction)");

	Int_t cellColumn=0,cellRow=0;
	Int_t cellColumnAbs=0,cellRowAbs=0;
	Int_t trash;

	for (Int_t cell = fStartCell; cell < fNoOfCells; cell++)
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
		if(flagEnd==-1  && fFlag[cell]==flagBegin)                        plot2D->SetBinContent(cellColumnAbs,cellRowAbs,1);
		if(flagEnd!=0   && flagEnd!=-1 && fFlag[cell]>=flagBegin && fFlag[cell]<=flagEnd)plot2D->SetBinContent(cellColumnAbs,cellRowAbs,1);
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
	TString name =Form("%s/%s/%s_%s.gif", fWorkdir.Data(),fAnalysisOutput.Data(),fPeriod.Data() , histoName.Data());
	c1->SaveAs(name);
    ///cout<<"gErrorIgnoreLevel: "<<gErrorIgnoreLevel<<endl;
	fRootFile->WriteObject(plot2D,plot2D->GetName());

}
///
/// This function saves all good cells amplitudes to a root file
///
//_________________________________________________________________________
void BadChannelAna::SaveHistoToFile()
{
	char name[100];
	for(Int_t cell=fStartCell;cell<fNoOfCells;cell++)
	{
		sprintf(name, "Cell %d",cell) ;
		TH1 *hCell = fCellAmplitude->ProjectionX(name,cell+1,cell+1);
		if(fFlag[cell]==0)fOutputListGood->Add(hCell);
	}
	fRootFile->WriteObject(fOutputListGood,fOutputListGood->GetName());
}
