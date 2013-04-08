// -------------------------------------------------------------------------
// Object to handle AOD Track cuts for the DiHadronPID analysis. Should
// inherit from: ?TObject or TNamed?, so that it can be written to a file.
//
// Still open question: should this object hold multiple track cuts, or
// should the analysis hold a TObjArray, with all the track cuts?
// -------------------------------------------------------------------------
// Author: Misha Veldhoen (misha.veldhoen@cern.ch)

#include <iostream>

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNamed.h"
#include "TFormula.h"
#include "TIterator.h"

// AOD includes.
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"

// PID includes.
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"

#include "AliTrackDiHadronPID.h"

#include "AliAODTrackCutsDiHadronPID.h"

using namespace std;

ClassImp(AliAODTrackCutsDiHadronPID);

// -------------------------------------------------------------------------
AliAODTrackCutsDiHadronPID::AliAODTrackCutsDiHadronPID():
	TNamed(),
	fMinPt(0.),
	fMaxPt(10.),
	fFilterMask(0),
	fMaxEta(-999.),
	fMaxRapidity(-999.),
	fDemandedFlags(0),
	fMinSPDHitsForPtDeptDCAcut(0),
	fPtDeptDCAxyCutFormula(0x0),
	fDCAzCut(999.),
	fIsMC(kFALSE),
	fTestFilterMask(kFALSE),
	fTestMaxEta(kFALSE),
	fTestMaxRapidity(kFALSE),
	fTestFlags(kFALSE),
	fTestTOFmismatch(kFALSE),
	fTestPtDeptDCAcut(kFALSE),
	fDataTrackQAHistos(0x0),
	fPrimRecMCTrackQAHistos(0x0),
	fPrimGenMCTrackQAHistos(0x0),
	fSecRecMCTrackQAHistos(0x0),
	fSecGenMCTrackQAHistos(0x0),
	fDebug(0)
	
 {

 	// 
 	// Default constructor
 	//

 	// Q: When an object is loaded from a TFile, then the default constructor is called.
 	//    I am assuming that the valiables initialized here are overwritten by the values
 	//    which were stored in the file, except for the data members which have an
 	//    exclamation mark in the header file. -> Check this! 

 	cout<<"AliAODTrackCutsDiHadronPID Default Constructor Called."<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	for (Int_t iHistoName = 0; iHistoName < 12; iHistoName++) {
 		
 		// Initialize the data histograms.
 		if (iHistoName < 3) {
 			fHistDataPt[iHistoName] = 0x0;
 		 	fHistDataNTracks[iHistoName] = 0x0;
  			fNTracks[iHistoName] = 0;

  			// DCA histograms. 
  			fHistDataDCAxy[iHistoName] = 0x0;
  			fHistDataDCAz[iHistoName] = 0x0;

  			// Initialize PID histograms.
  			for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
  				for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
  					fHistDataPID[iHistoName][iSpecies][iPtClass] = 0x0;
   					fHistTOFMismatch[iHistoName][iSpecies][iPtClass] = 0x0; 				  					
  				}
  			}
 		}
 		// Initialize data DCA for one sigma.
 		fHistDataDCAxyOneSigma[iHistoName] = 0x0;

 		// Initialize MC histograms.
 		fHistPrimGenMCPt[iHistoName] = 0x0;
 		fHistPrimRecMCPt[iHistoName] = 0x0;
 		fHistPrimRecNTracks[iHistoName] = 0x0;

 		fHistSecGenMCPt[iHistoName] = 0x0;
 		fHistSecRecMCPt[iHistoName] = 0x0;

 		// Initialze MC DCA histograms
 		fHistPrimRecMCDCA[iHistoName] = 0x0;
 		fHistSecRecMCDCAMat[iHistoName] = 0x0;
 		fHistSecRecMCDCAWeak[iHistoName] = 0x0;

 		// Initialize other stuff.
 		fHistRequested[iHistoName] = kFALSE;
		f3DSpectraEnabeled[iHistoName] = kFALSE;
		fPIDHistosEnabeled[iHistoName] = kFALSE;
 	}

 	// TODO: It would be a bit neater to initialize all values to zero
 	// instead of the default values...
 	InitializeDefaultHistoNamesAndAxes();

}

// -------------------------------------------------------------------------
AliAODTrackCutsDiHadronPID::AliAODTrackCutsDiHadronPID(const char* name):
	TNamed(name,"AOD Track Cuts"),
	fMinPt(0.),
	fMaxPt(10.),
	fFilterMask(0),
	fMaxEta(-999.),
	fMaxRapidity(-999.),
	fDemandedFlags(0),
	fMinSPDHitsForPtDeptDCAcut(0),
	fPtDeptDCAxyCutFormula(0x0),
	fDCAzCut(999.),
	fIsMC(kFALSE),
	fTestFilterMask(kFALSE),
	fTestMaxEta(kFALSE),
	fTestMaxRapidity(kFALSE),
	fTestFlags(kFALSE),
	fTestTOFmismatch(kFALSE),
	fTestPtDeptDCAcut(kFALSE),
	fDataTrackQAHistos(0x0),
	fPrimRecMCTrackQAHistos(0x0),
	fPrimGenMCTrackQAHistos(0x0),
	fSecRecMCTrackQAHistos(0x0),
	fSecGenMCTrackQAHistos(0x0),
	fDebug(0)

{

	//
	// Named constructor
	//

	cout<<"AliAODTrackCutsDiHadronPID Named Constructor Called."<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}
		
	for (Int_t iHistoName = 0; iHistoName < 12; iHistoName++) {
 		
 		// Initialize the data histograms.
 		if (iHistoName < 3) {
 			fHistDataPt[iHistoName] = 0x0;
 		 	fHistDataNTracks[iHistoName] = 0x0;
  			fNTracks[iHistoName] = 0;

  			// DCA histograms. 
  			fHistDataDCAxy[iHistoName] = 0x0;
  			fHistDataDCAz[iHistoName] = 0x0;

 			// Initialize PID histograms.
  			for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
  				for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
  					fHistDataPID[iHistoName][iSpecies][iPtClass] = 0x0;
   					fHistTOFMismatch[iHistoName][iSpecies][iPtClass] = 0x0; 				  					
  				}
  			}
 		}
		// Initialize data DCA for one sigma.
 		fHistDataDCAxyOneSigma[iHistoName] = 0x0;

 		// Initialize MC histograms.
 		fHistPrimGenMCPt[iHistoName] = 0x0;
 		fHistPrimRecMCPt[iHistoName] = 0x0;
 		fHistPrimRecNTracks[iHistoName] = 0x0;

 		fHistSecGenMCPt[iHistoName] = 0x0;
 		fHistSecRecMCPt[iHistoName] = 0x0;

 		// Initialze MC DCA histograms
 		fHistPrimRecMCDCA[iHistoName] = 0x0;
		fHistSecRecMCDCAMat[iHistoName] = 0x0;
 		fHistSecRecMCDCAWeak[iHistoName] = 0x0;

 		// Initialize other stuff.
 		fHistRequested[iHistoName] = kFALSE;
		f3DSpectraEnabeled[iHistoName] = kFALSE;
		fPIDHistosEnabeled[iHistoName] = kFALSE;
 	}

 	InitializeDefaultHistoNamesAndAxes();

}

// -------------------------------------------------------------------------
void AliAODTrackCutsDiHadronPID::InitializeDefaultHistoNamesAndAxes() {

	// Initializes the histogram name conventions and the ranges of all the histograms
	// by their default values. This method should only be called by the
	// (default) constructor.

	// TODO: User should be able to change these standard values at initialization.
	// TODO: User should be able to retrieve all these variables with appropriate Getters.

	cout<<"AliAODTrackCutsDiHadronPID - Initializing Default Histogram Names and axes..."<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

 	// Setting the Pt axis for all histograms except the PID and Mismatch histograms.
 	Double_t ptaxis[47] = {0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,
							1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,
							3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
	for (Int_t iPtBins = 0; iPtBins < 47; iPtBins++) {fPtAxis[iPtBins] = ptaxis[iPtBins];}
 	fNPtBins = 46;

	// Setting the pt range of the five PID histogram pt classes. 	
	Double_t ptboundarypid[6] = {0.5,0.7,1.0,1.7,3.0,5.0};
	for (Int_t iPtBoundary = 0; iPtBoundary < 6; iPtBoundary++) {fPtBoundaryPID[iPtBoundary] = ptboundarypid[iPtBoundary];}

	// Setting the number of bins in pt for the five PID histogram pt classes.
	Int_t nptbinspid[5] = {8,6,14,13,10};
	for (Int_t iPtBins = 0; iPtBins < 5; iPtBins++) {fNPtBinsPID[iPtBins] = nptbinspid[iPtBins];}

	// Setting the TOF axes for the PID histograms.
	Double_t toflowerbound[5][3] = {{-2000.,-6000.,-10000.},{-2000.,-4000.,-10000.},{-1000.,-2000.,-5000.},{-1000.,-1000.,-2500.},{-500.,-500.,-1000.}};
	Double_t tofupperbound[5][3] = {{10000.,10000.,6000.},{10000.,8000.,6000.},{6000.,6000.,6000.},{6000.,6000.,6000.},{4000.,4000.,6000.}};
	Int_t tofbins[5][3] = {{120,160,160},{120,120,160},{140,140,165},{140,140,170},{90,90,140}};
	for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
		for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
			fTOFLowerBound[iPtClass][iSpecies] = toflowerbound[iPtClass][iSpecies];
			fTOFUpperBound[iPtClass][iSpecies] = tofupperbound[iPtClass][iSpecies];
			fTOFbins[iPtClass][iSpecies] = tofbins[iPtClass][iSpecies];
		}
	}

	// Setting the TPC axes for the PID histograms.
	Double_t tpclowerbound[5][3] = {{-20.,-50.,-100.},{-20.,-30.,-80.},{-25.,-25.,-45.},{-25.,-25.,-45.},{-25.,-20.,-20.}};
	Double_t tpcupperbound[5][3] = {{60.,30.,20.},{60.,40.,20.},{50.,50.,25.},{45.,45.,25.},{25.,30.,30.}}; // Check highest pT bin boundaries for K,p
	Int_t tpcbins[5][3] = {{80,80,120},{80,70,100},{75,75,70},{70,70,70},{50,50,50}};
	for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
		for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
			fTPCLowerBound[iPtClass][iSpecies] = tpclowerbound[iPtClass][iSpecies];
			fTPCUpperBound[iPtClass][iSpecies] = tpcupperbound[iPtClass][iSpecies];
			fTPCbins[iPtClass][iSpecies] = tpcbins[iPtClass][iSpecies];
		}
	}

	// Names for the 12 (species,charge) combinations considered.
	const char* histoname[12] = {"AllCharged","Pos","Neg","AllPion","PosPion","NegPion","AllKaon","PosKaon","NegKaon","AllProton","PosProton","NegProton"};
	const char* histolatex[12] = {"ch","ch^{+}","ch^{-}","#pi^{+/-}","#pi^{+}","#pi^{-}","K^{+/-}","K^{+}","K^{-}","p/#bar{p}","p","#bar{p}"};
	for (Int_t iHistoName = 0; iHistoName < 12; iHistoName++) {
 		fHistoName[iHistoName] = histoname[iHistoName];
 		fHistoLatex[iHistoName] = histolatex[iHistoName];
 	}

 	// Names of the 3 species considered.
	const char* particlename[3] = {"Pion","Kaon","Proton"};
	for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {fParticleName[iSpecies] = particlename[iSpecies];}

	// Names of the 5 pt classes considered for the PID histograms.
	const char* ptclassname[5] = {"VeryLowPt","LowPt","MedPt","MedHighPt","HighPt"};
	for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {fPtClassName[iPtClass] = ptclassname[iPtClass];}

}

// -------------------------------------------------------------------------
AliAODTrackCutsDiHadronPID::~AliAODTrackCutsDiHadronPID() {

	//
	// Destructor
	//

	cout<<"AliAODTrackCutsDiHadronPID Destructor Called."<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

 	if (fDataTrackQAHistos) delete fDataTrackQAHistos;
 	fDataTrackQAHistos = 0x0;
 	if (fPrimRecMCTrackQAHistos) delete fPrimRecMCTrackQAHistos;
 	fPrimRecMCTrackQAHistos = 0x0;
 	if (fPrimGenMCTrackQAHistos) delete fPrimGenMCTrackQAHistos;
 	fPrimGenMCTrackQAHistos = 0x0;
 	if (fSecRecMCTrackQAHistos) delete fSecRecMCTrackQAHistos;
 	fSecRecMCTrackQAHistos = 0x0;
 	if (fSecGenMCTrackQAHistos) delete fSecGenMCTrackQAHistos;
	fSecGenMCTrackQAHistos = 0x0;

}

// -------------------------------------------------------------------------
Long64_t AliAODTrackCutsDiHadronPID::Merge(TCollection* list) {

	//
	// Merger. 
	// 

	cout<<"AliAODTrackCutsDiHadronPID Merger Called."<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	Bool_t HistosOK = kTRUE;

	if (fIsMC) {
		if (!fPrimGenMCTrackQAHistos||!fSecGenMCTrackQAHistos||!fPrimRecMCTrackQAHistos||!fSecRecMCTrackQAHistos) {HistosOK = kFALSE;}
	} else {
		if (!fDataTrackQAHistos) {HistosOK = kFALSE;}
	}

	if (!HistosOK) {
		cout<<"AliAODTrackCutsDiHadronPID::Merge() - Warning, current object's histograms are missing... Generating."<<endl;
		CreateHistos();
	}

	if (!list) {
		//cout<<"No list available..."<<endl;
		return 0;
	}
	if (list->IsEmpty()) {
		//cout<<"List is empty..."<<endl;
		return 1;
	}

	//Int_t NEntries = list->GetEntries();
	//cout<<"Supplied list has "<<NEntries<<" entries."<<endl;

	TIterator* iter = list->MakeIterator();
	TObject* obj;

	// List of collections
	TList collection_fDataTrackQAHistos;
	TList collection_fPrimRecMCTrackQAHistos;
	TList collection_fPrimGenMCTrackQAHistos;
	TList collection_fSecRecMCTrackQAHistos;
	TList collection_fSecGenMCTrackQAHistos;

	Int_t count = 0;

  	while ((obj = iter->Next())) {
    	AliAODTrackCutsDiHadronPID* entry = dynamic_cast<AliAODTrackCutsDiHadronPID*> (obj);
    	if (entry == 0) continue;

    	// Check if the object to be merged really has the same name! (FIXME!)

    	// Getting the lists from obj.
    	TList* list_fDataTrackQAHistos = entry->GetListOfDataQAHistos();
    	TList* list_fPrimRecMCTrackQAHistos = entry->GetListOfPrimRecMCTrackQAHistos();
    	TList* list_fPrimGenMCTrackQAHistos = entry->GetListOfPrimGenMCTrackQAHistos();
    	TList* list_fSecRecMCTrackQAHistos = entry->GetListOfSecRecMCTrackQAHistos();
    	TList* list_fSecGenMCTrackQAHistos = entry->GetListOfSecGenMCTrackQAHistos();

    	// Adding the retrieved lists to the collection.
    	if (list_fDataTrackQAHistos) collection_fDataTrackQAHistos.Add(list_fDataTrackQAHistos);
    	if (list_fPrimRecMCTrackQAHistos) collection_fPrimRecMCTrackQAHistos.Add(list_fPrimRecMCTrackQAHistos);
    	if (list_fPrimGenMCTrackQAHistos) collection_fPrimGenMCTrackQAHistos.Add(list_fPrimGenMCTrackQAHistos);
    	if (list_fSecRecMCTrackQAHistos) collection_fSecRecMCTrackQAHistos.Add(list_fSecRecMCTrackQAHistos);
    	if (list_fSecGenMCTrackQAHistos) collection_fSecGenMCTrackQAHistos.Add(list_fSecGenMCTrackQAHistos);

    	//cout<<"Entries intermediate list Data: "<<collection_fDataTrackQAHistos.GetEntries()<<endl;
    	//cout<<"Entries intermediate list RecMC: "<<collection_fPrimRecMCTrackQAHistos.GetEntries()<<endl;

    	count++;
    }

    // Merging. Note that we require the original list to exist.
    //  * Assume that if the collection happens to be empty, then nothing will happen.
    //  |-> This of course leads to some problems, since if the first file does not
    //  |   have such a list, then the merged file will have no results...
    //  |   IDEA: if fDataTrackQAHistos does not exist, then create a new one. (TO DO)
    //  * All other variables are taken from the original object.
    if (fDataTrackQAHistos) fDataTrackQAHistos->Merge(&collection_fDataTrackQAHistos);
    if (fPrimRecMCTrackQAHistos) fPrimRecMCTrackQAHistos->Merge(&collection_fPrimRecMCTrackQAHistos);
    if (fPrimGenMCTrackQAHistos) fPrimGenMCTrackQAHistos->Merge(&collection_fPrimGenMCTrackQAHistos);    
    if (fSecRecMCTrackQAHistos) fSecRecMCTrackQAHistos->Merge(&collection_fSecRecMCTrackQAHistos);
    if (fSecGenMCTrackQAHistos) fSecGenMCTrackQAHistos->Merge(&collection_fSecGenMCTrackQAHistos);

    delete iter;

	return count+1;
	
}

// -------------------------------------------------------------------------
void AliAODTrackCutsDiHadronPID::PrintCuts() {}

// -------------------------------------------------------------------------
void AliAODTrackCutsDiHadronPID::CreateHistos() {

	// Function should be called by the UserCreateOutput() function of the
	// analysis task. This function then generates all the histograms that
	// were requested locally on the workernode. Even if case that the
	// histograms are not filled, it is imperative that they are still
	// created, and in the same order, otherwise problems with merging will
	// arise.

	// TODO: In principle it should never happen that if this method is called,
	//       that the lists of QA histograms already exist, so this may become a fatal error
	//       instead of a warning.

	cout<<"AliAODTrackCutsDiHadronPID - Creating histograms for cuts: "<<fName<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (fIsMC) {
	
		if (!fPrimGenMCTrackQAHistos) {
			cout<<"AliAODTrackCutsDiHadronPID - Creating Prim. Gen. MC Track QA TList..."<<endl;	
			fPrimGenMCTrackQAHistos = new TList();
			fPrimGenMCTrackQAHistos->SetName("PrimGenMCTrackQAHistos");
			fPrimGenMCTrackQAHistos->SetOwner(kTRUE);			
		} else {cout<<"AliAODTrackCutsDiHadronPID - Warning, Prim. Gen. MC Track QA TList was already created..."<<endl;}
		
		if (!fSecGenMCTrackQAHistos) {
			cout<<"AliAODTrackCutsDiHadronPID - Creating Sec. Gen. MC Track QA TList..."<<endl;				
			fSecGenMCTrackQAHistos = new TList();
			fSecGenMCTrackQAHistos->SetName("SecGenMCTrackQAHistos");
			fSecGenMCTrackQAHistos->SetOwner(kTRUE);
			for (Int_t iHistoClass = 0; iHistoClass < 12; iHistoClass++) {if (fHistRequested[iHistoClass]) {InitializeGenMCHistos(iHistoClass);}}				
		} else {cout<<"AliAODTrackCutsDiHadronPID - Warning, Sec. Gen. MC Track QA TList was already created..."<<endl;}

		if (!fPrimRecMCTrackQAHistos) {
			cout<<"AliAODTrackCutsDiHadronPID - Creating Prim. Rec. MC Track QA TList..."<<endl;				
			fPrimRecMCTrackQAHistos = new TList();
			fPrimRecMCTrackQAHistos->SetName("PrimRecMCTrackQAHistos");
			fPrimRecMCTrackQAHistos->SetOwner(kTRUE);
		} else {cout<<"AliAODTrackCutsDiHadronPID - Warning, Prim. Rec. MC Track QA TList was already created..."<<endl;}

		if (!fSecRecMCTrackQAHistos) {
			cout<<"AliAODTrackCutsDiHadronPID - Creating Sec. Rec. MC Track QA TList..."<<endl;				
			fSecRecMCTrackQAHistos = new TList();
			fSecRecMCTrackQAHistos->SetName("SecRecMCTrackQAHistos");
			fSecRecMCTrackQAHistos->SetOwner(kTRUE);
			for (Int_t iHistoClass = 0; iHistoClass < 12; iHistoClass++) {
				if (fHistRequested[iHistoClass]) {
					InitializeRecMCHistos(iHistoClass);
				}
			}	
		} else {cout<<"AliAODTrackCutsDiHadronPID - Warning, Sec. Rec. MC Track QA TList was already created..."<<endl;}

	} else {
		
		if (!fDataTrackQAHistos) {
			cout<<"AliAODTrackCutsDiHadronPID - Creating Data Track QA TList..."<<endl;
			fDataTrackQAHistos = new TList();
			fDataTrackQAHistos->SetName("DataTrackQAHistos");
			fDataTrackQAHistos->SetOwner(kTRUE);
			for (Int_t iHistoClass = 0; iHistoClass < 12; iHistoClass++) {
				fHistDataDCAxyOneSigma[iHistoClass] = InitializeDCASpectrum("fHistDataDCAxyOneSigma",iHistoClass);
				fDataTrackQAHistos->Add(fHistDataDCAxyOneSigma[iHistoClass]);
				if (fHistRequested[iHistoClass]) { if (iHistoClass < 3) {InitializeDataHistos(iHistoClass);}}
			}	
		} else {cout<<"AliAODTrackCutsDiHadronPID - Warning, Data QA TList was already created..."<<endl;}
	
	}
}

// -------------------------------------------------------------------------
void AliAODTrackCutsDiHadronPID::StartNewEvent() {

	// FIXME: This method is now only suited for Data (3 histo classes only.) 
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Resetting the counters.
	for (Int_t iHistoClass = 0; iHistoClass < 3; iHistoClass++) {
		fNTracks[iHistoClass] = 0;
 	}
}

// -------------------------------------------------------------------------
void AliAODTrackCutsDiHadronPID::EventIsDone(Bool_t IsMC) {

	// FIXME: This method is now only suited for Data (3 histo classes only.) 
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Fill NTracks histos.

	for (Int_t iHistoClass = 0; iHistoClass < 3; iHistoClass++) {

		if (fHistRequested[iHistoClass]) {
			
			if (IsMC) {
				// THE FOLLOWING SHOULD NEVER HAPPEN.
				//if (!fHistPrimRecMCPt[iHistoClass]) InitializeRecMCHistos(iHistoClass);
				fHistPrimRecNTracks[iHistoClass]->Fill(fNTracks[iHistoClass]);
			} else {
				// THE FOLLOWING SHOULD NEVER HAPPEN.
				//if (!fHistDataPt[iHistoClass]) InitializeDataHistos(iHistoClass);
				fHistDataNTracks[iHistoClass]->Fill(fNTracks[iHistoClass]);
			}
		}
	}
}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::IsSelectedData(AliTrackDiHadronPID* track, Double_t randomhittime) {
	
	//
	// Checks performed on a data track.
	//

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (!track) return kFALSE;

	if (!fDataTrackQAHistos) {cout<<"AliAODTrackCutsDiHadronPID - Histograms were not created, you should have called CreateHistos()..."<<endl;}

	// Check the track cuts (NB. check function return kTRUE if track is accepted.)
	if (!CheckPt(track->Pt())) return kFALSE;
	if (!CheckMaxEta(track->Eta())) return kFALSE;
	if (!CheckFilterMask(track->GetFilterMap())) return kFALSE;
	if (!CheckFlags(track->GetFlags())) return kFALSE;
	if (!CheckTOFmismatch(track->IsTOFmismatch())) return kFALSE;

	Int_t NSPDhits = 0;
	if (track->HasPointOnITSLayer(0)) NSPDhits++;
	if (track->HasPointOnITSLayer(1)) NSPDhits++;
	if (!CheckPtDeptDCACut(track->GetZAtDCA(),track->GetXYAtDCA(),track->Pt(),NSPDhits)) return kFALSE;

	// Track has passed the cuts, fill QA histograms.
	for (Int_t iHistoClass = 0; iHistoClass < 3; iHistoClass++) {

		if (fHistRequested[iHistoClass]) {

			// Check the charge (could be neater).
			if ((iHistoClass == 0) && (track->Charge() == 0)) continue;
			if ((iHistoClass == 1) && (track->Charge() <= 0)) continue;
			if ((iHistoClass == 2) && (track->Charge() >= 0)) continue; 

			FillDataHistos(iHistoClass, track);

			// Ignore if random hit is < -1.e20.
			if (randomhittime > -1.e20) FillTOFMismatchHistos(iHistoClass,track,randomhittime);

			fNTracks[iHistoClass]++;

		}
	}

	// Track Has Passed.
	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::IsSelectedGeneratedMC(AliAODMCParticle* particle) {
	
	//
	// Checks performed on a generated MC particle.
	//

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// PDG codes for particles:
	//  pi+ 211, K+ 321, p 2212

	if (!particle) return kFALSE;

	if (!fPrimGenMCTrackQAHistos||!fSecGenMCTrackQAHistos) {cout<<"AliAODTrackCutsDiHadronPID - Histograms were not created, you should have called CreateHistos()..."<<endl;}

	// Check the track cuts.
	if (!CheckPt(particle->Pt())) return kFALSE;
	if (!CheckMaxEta(particle->Eta())) return kFALSE;
	if (!CheckRapidity(particle->Y())) return kFALSE;		// NEW: rapidity cut. (not done on data)

	// NOT YET IMPLEMENTED! - Histoclasses for specific particles.
	for (Int_t iHistoClass = 0; iHistoClass < 12; iHistoClass++) {

		if (fHistRequested[iHistoClass]) {

			// Check the charge (could be neater).
			if ((iHistoClass == 0) && (particle->Charge() == 0)) continue;
			if ((iHistoClass == 1) && (particle->Charge() <= 0)) continue;
			if ((iHistoClass == 2) && (particle->Charge() >= 0)) continue;
			if ((iHistoClass == 3) && (TMath::Abs(particle->GetPdgCode()) != 211)) continue;
			if ((iHistoClass == 4) && (particle->GetPdgCode()) != 211) continue;
			if ((iHistoClass == 5) && (particle->GetPdgCode()) != -211) continue;
			if ((iHistoClass == 6) && (TMath::Abs(particle->GetPdgCode()) != 321)) continue;
			if ((iHistoClass == 7) && (particle->GetPdgCode()) != 321) continue;
			if ((iHistoClass == 8) && (particle->GetPdgCode()) != -321) continue;
			if ((iHistoClass == 9) && (TMath::Abs(particle->GetPdgCode()) != 2212)) continue;
			if ((iHistoClass == 10) && (particle->GetPdgCode()) != 2212) continue;
			if ((iHistoClass == 11) && (particle->GetPdgCode()) != -2212) continue;

			// Secondary specification not set.
			//cout<<"Charge: "<<particle->Charge()<<" PDG code: "<<particle->GetPdgCode()<<endl;
			//cout<<particle->IsPhysicalPrimary()<<" "<<particle->IsSecondaryFromWeakDecay()<<" "<<particle->IsSecondaryFromMaterial()<<endl;

			// These two functions are not implemented...
			if (particle->IsSecondaryFromWeakDecay()) cout<<"Secondary From Weak Decay!"<<endl;
			if (particle->IsSecondaryFromMaterial()) cout<<"Secondary From Material!"<<endl;

			FillGenMCHistos(iHistoClass, particle);	

		}
	}

	// Particle has Passed.
	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::IsSelectedReconstructedMC(AliTrackDiHadronPID* track) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (!track) return kFALSE;

	if (!fPrimRecMCTrackQAHistos||!fSecRecMCTrackQAHistos) {cout<<"AliAODTrackCutsDiHadronPID - Histograms were not created, you should have called CreateHistos()..."<<endl;}

	// Check the track cuts.
	if (!CheckPt(track->MCPt())) return kFALSE;					// Kinematic cuts are done on the MC particles
	if (!CheckMaxEta(track->MCEta())) return kFALSE;			// and not on the data.
	if (!CheckRapidity(track->MCY())) return kFALSE;			// NEW: rapidity cut.
	if (!CheckFilterMask(track->GetFilterMap())) return kFALSE;
	if (!CheckFlags(track->GetFlags())) return kFALSE;

	Int_t NSPDhits = 0;
	if (track->HasPointOnITSLayer(0)) NSPDhits++;
	if (track->HasPointOnITSLayer(1)) NSPDhits++;
	if (!CheckPtDeptDCACut(track->GetZAtDCA(),track->GetXYAtDCA(),track->Pt(),NSPDhits)) return kFALSE;

	// Track has passed the cuts, fill QA histograms.
	for (Int_t iHistoClass = 0; iHistoClass < 12; iHistoClass++) {

		if (fHistRequested[iHistoClass]) {

			// Check the charge (could be neater).
			if ((iHistoClass == 0) && (track->Charge() == 0)) continue;
			if ((iHistoClass == 1) && (track->Charge() <= 0)) continue;
			if ((iHistoClass == 2) && (track->Charge() >= 0)) continue;
			if ((iHistoClass == 3) && (TMath::Abs(track->GetPdgCode()) != 211)) continue;
			if ((iHistoClass == 4) && (track->GetPdgCode()) != 211) continue;
			if ((iHistoClass == 5) && (track->GetPdgCode()) != -211) continue;
			if ((iHistoClass == 6) && (TMath::Abs(track->GetPdgCode()) != 321)) continue;
			if ((iHistoClass == 7) && (track->GetPdgCode()) != 321) continue;
			if ((iHistoClass == 8) && (track->GetPdgCode()) != -321) continue;
			if ((iHistoClass == 9) && (TMath::Abs(track->GetPdgCode()) != 2212)) continue;
			if ((iHistoClass == 10) && (track->GetPdgCode()) != 2212) continue;
			if ((iHistoClass == 11) && (track->GetPdgCode()) != -2212) continue;

			FillRecMCHistos(iHistoClass, track);

			fNTracks[iHistoClass]++;
		}
	}

	// Track Has Passed.
	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::FillDataHistos(Int_t histoclass, AliTrackDiHadronPID* track) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Fill the histograms.
	fHistDataPt[histoclass]->Fill(track->Pt());

	// Fill DCA histograms.
	fHistDataDCAxy[histoclass]->Fill(track->GetXYAtDCA());
	fHistDataDCAz[histoclass]->Fill(track->GetZAtDCA());

	// Fill DCA_{xy} one sigma histograms.
	Int_t checkSum = 0;
	// histoclass = 0: all charges; histoclass = 1: positive; histoclass = 2: negative;
	if (TMath::Sqrt(track->GetNumberOfSigmasTOF(0) * track->GetNumberOfSigmasTOF(0) + 
					track->GetNumberOfSigmasTPC(0) * track->GetNumberOfSigmasTPC(0)) < 1.) {
		fHistDataDCAxyOneSigma[0 + histoclass]->Fill(track->Pt(),track->GetXYAtDCA());	// All species.
		fHistDataDCAxyOneSigma[3 + histoclass]->Fill(track->Pt(),track->GetXYAtDCA());	// Pions.
	//	checkSum++;	cout<<"Pion found: nSigTOF: "<<track->GetNumberOfSigmasTOF(0)<<"; nSigTPC: "<<track->GetNumberOfSigmasTPC(0)<<endl;
	}
	if (TMath::Sqrt(track->GetNumberOfSigmasTOF(1) * track->GetNumberOfSigmasTOF(1) + 
					track->GetNumberOfSigmasTPC(1) * track->GetNumberOfSigmasTPC(1)) < 1.) {
		fHistDataDCAxyOneSigma[0 + histoclass]->Fill(track->Pt(),track->GetXYAtDCA());	// All species.
		fHistDataDCAxyOneSigma[6 + histoclass]->Fill(track->Pt(),track->GetXYAtDCA());	// Kaons.
	//	checkSum++;	cout<<"Kaon found: nSigTOF: "<<track->GetNumberOfSigmasTOF(1)<<"; nSigTPC: "<<track->GetNumberOfSigmasTPC(1)<<endl;
	}
	if (TMath::Sqrt(track->GetNumberOfSigmasTOF(2) * track->GetNumberOfSigmasTOF(2) + 
					track->GetNumberOfSigmasTPC(2) * track->GetNumberOfSigmasTPC(2)) < 1.) {		
		fHistDataDCAxyOneSigma[0 + histoclass]->Fill(track->Pt(),track->GetXYAtDCA());	// All species.
		fHistDataDCAxyOneSigma[9 + histoclass]->Fill(track->Pt(),track->GetXYAtDCA());	// Protons.
	//	checkSum++;	cout<<"Proton found: nSigTOF: "<<track->GetNumberOfSigmasTOF(2)<<"; nSigTPC: "<<track->GetNumberOfSigmasTPC(2)<<endl;
	}

	// check for double identification (or triple..)
	if(checkSum > 1) {AliError("More than one particle identified for the same track!"); }

	// Fill PID histos.
	for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

		// Note that a possible Y cut is only done on the PID cuts!
		if (!CheckRapidity(track->Y(iSpecies))) continue;

		for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
			fHistDataPID[histoclass][iSpecies][iPtClass]->Fill(track->Pt(),track->GetTOFsignalMinusExpected(iSpecies),track->GetTPCsignalMinusExpected(iSpecies));		
		}
	}

	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::FillTOFMismatchHistos(Int_t histoclass, AliTrackDiHadronPID* track, Double_t randomhittime) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Fill TOF mismatch.
	for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

		// Note that a possible Y cut is only done on the PID cuts!
		if (!CheckRapidity(track->Y(iSpecies))) continue;

		for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
			fHistTOFMismatch[histoclass][iSpecies][iPtClass]->Fill(track->Pt(),randomhittime - track->GetTOFsignalExpected(iSpecies));		
		}
	}

	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::FillGenMCHistos(Int_t histoclass, AliAODMCParticle* particle) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Fill the histograms.
	if (particle->IsPhysicalPrimary()) {
		fHistPrimGenMCPt[histoclass]->Fill(particle->Pt());
	} else {
		fHistSecGenMCPt[histoclass]->Fill(particle->Pt());
	}

	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::FillRecMCHistos(Int_t histoclass, AliTrackDiHadronPID* track) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Fill the Pt histograms.
	if (track->IsPhysicalPrimary()) {
		fHistPrimRecMCPt[histoclass]->Fill(track->Pt());
	} else {
		fHistSecRecMCPt[histoclass]->Fill(track->Pt());
	}

	// Fill the DCA histograms.
	if (track->IsPhysicalPrimary()) {
		fHistPrimRecMCDCA[histoclass]->Fill(track->Pt(), track->GetXYAtDCA());
	} 
	if (track->IsSecondaryFromMaterial()) {
		fHistSecRecMCDCAMat[histoclass]->Fill(track->Pt(), track->GetXYAtDCA());
	}
	if (track->IsSecondaryFromWeakDecay()) {
		fHistSecRecMCDCAWeak[histoclass]->Fill(track->Pt(), track->GetXYAtDCA());
	}

	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::InitializeDataHistos(Int_t histoclass) {

	cout<<"AliAODTrackCutsDiHadronPID - Creating Data Histograms of Class "<<histoclass<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Add the Pt spectra.
	fHistDataPt[histoclass] = InitializePtSpectrum("fHistDataPt",histoclass);
	fDataTrackQAHistos->Add(fHistDataPt[histoclass]);

	// Add the NTrack histograms.
	fHistDataNTracks[histoclass] = InitializeNTracksHisto("fHistDataNTracks",histoclass);
	fDataTrackQAHistos->Add(fHistDataNTracks[histoclass]);

	// Add the DCA histograms.
	fHistDataDCAxy[histoclass] = InitializeDCAxyHisto("fHistDCAxy",histoclass);
	fDataTrackQAHistos->Add(fHistDataDCAxy[histoclass]);
	fHistDataDCAz[histoclass] = InitializeDCAzHisto("fHistDCAz",histoclass);
	fDataTrackQAHistos->Add(fHistDataDCAz[histoclass]);

	// Add the PID histograms. (FIXME shoudl be able to turn the creation of these histos off.)
	for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
		for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
			fHistDataPID[histoclass][iSpecies][iPtClass] = InitializePIDHisto("fHistDataPID",histoclass,iSpecies,iPtClass);
			fDataTrackQAHistos->Add(fHistDataPID[histoclass][iSpecies][iPtClass]);
			fHistTOFMismatch[histoclass][iSpecies][iPtClass] = InitializeTOFMismatchHisto("fHistTOFMismatch",histoclass,iSpecies,iPtClass);
			fDataTrackQAHistos->Add(fHistTOFMismatch[histoclass][iSpecies][iPtClass]);
		}
	}

	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::InitializeGenMCHistos(Int_t histoclass) {

	cout<<"AliAODTrackCutsDiHadronPID - Creating Generated MC Histograms of Class "<<histoclass<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Primary Particles.
	fHistPrimGenMCPt[histoclass] = InitializePtSpectrum("fHistPrimGenMCPt",histoclass);
	fPrimGenMCTrackQAHistos->Add(fHistPrimGenMCPt[histoclass]);	

	// Secondary Particles.
	fHistSecGenMCPt[histoclass] = InitializePtSpectrum("fHistSecGenMCPt",histoclass);
	fSecGenMCTrackQAHistos->Add(fHistSecGenMCPt[histoclass]);	

	return kTRUE;

}

// -------------------------------------------------------------------------
Bool_t AliAODTrackCutsDiHadronPID::InitializeRecMCHistos(Int_t histoclass) {

	cout<<"AliAODTrackCutsDiHadronPID - Creating Reconstructed MC Histograms of Class "<<histoclass<<endl;
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Primary Particles.
	fHistPrimRecMCPt[histoclass] = InitializePtSpectrum("fHistPrimRecMCPt",histoclass);
	fPrimRecMCTrackQAHistos->Add(fHistPrimRecMCPt[histoclass]);

	fHistPrimRecNTracks[histoclass] = InitializeNTracksHisto("fHistPrimRecNTracks",histoclass);
	fPrimRecMCTrackQAHistos->Add(fHistPrimRecNTracks[histoclass]);

	fHistPrimRecMCDCA[histoclass] = InitializeDCASpectrum("fHistPrimRecDCA",histoclass);
	fPrimRecMCTrackQAHistos->Add(fHistPrimRecMCDCA[histoclass]);

	// Secondary Particles.
	fHistSecRecMCPt[histoclass] = InitializePtSpectrum("fHistSecRecMCPt",histoclass);
	fSecRecMCTrackQAHistos->Add(fHistSecRecMCPt[histoclass]);

	fHistSecRecMCDCAMat[histoclass] = InitializeDCASpectrum("fHistSecRecDCAMat",histoclass);
	fSecRecMCTrackQAHistos->Add(fHistSecRecMCDCAMat[histoclass]);

	fHistSecRecMCDCAWeak[histoclass] = InitializeDCASpectrum("fHistSecRecDCAWeak",histoclass);
	fSecRecMCTrackQAHistos->Add(fHistSecRecMCDCAWeak[histoclass]);

	return kTRUE;

}

// -------------------------------------------------------------------------
TH1F* AliAODTrackCutsDiHadronPID::InitializePtSpectrum(const char* name, Int_t histoclass) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	TH1F* hout = new TH1F(Form("%s%s",name,fHistoName[histoclass].Data()),
		Form("p_{T} Spectrum (%s);p_{T} (GeV/c);Count",fHistoName[histoclass].Data()),fNPtBins,fPtAxis);

	hout->SetDirectory(0);

	return hout;

}

// -------------------------------------------------------------------------
TH2F* AliAODTrackCutsDiHadronPID::InitializeDCASpectrum(const char* name, Int_t histoclass){

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	TH2F* hout = new TH2F(Form("%s%s",name,fHistoName[histoclass].Data()),
		Form("DCA_{xy} (%s); p_{T} (GeV/c); DCA_{xy} (cm)",fHistoName[histoclass].Data()),
		fNPtBins,fPtAxis,
		300,-3.,3.); 

	hout->SetDirectory(0);

	return hout;

}

// -------------------------------------------------------------------------
TH1F* AliAODTrackCutsDiHadronPID::InitializeNTracksHisto(const char* name, Int_t histoclass) {

	TH1F* hout = new TH1F(Form("%s%s",name,fHistoName[histoclass].Data()),
		Form("Number of Accepted Tracks (%s);N%s;Count",fHistoName[histoclass].Data(),fHistoLatex[histoclass].Data()),
		100,0,4000);

	hout->SetDirectory(0);

	return hout;

}

// -------------------------------------------------------------------------
TH1F* AliAODTrackCutsDiHadronPID::InitializeDCAxyHisto(const char* name, Int_t histoclass) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	TH1F* hout = new TH1F(Form("%s%s",name,fHistoName[histoclass].Data()),
		Form("DCAxy (%s);DCAxy (cm);Count",fHistoName[histoclass].Data()),300,-15.,15.);

	hout->SetDirectory(0);

	return hout;

}

// -------------------------------------------------------------------------
TH1F* AliAODTrackCutsDiHadronPID::InitializeDCAzHisto(const char* name, Int_t histoclass) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	TH1F* hout = new TH1F(Form("%s%s",name,fHistoName[histoclass].Data()),
	Form("DCAz (%s);DCAz (cm);Count",fHistoName[histoclass].Data()),300,-15.,15.);

	hout->SetDirectory(0);

	return hout;

}

// -------------------------------------------------------------------------
TH3F* AliAODTrackCutsDiHadronPID::InitializeAcceptanceHisto(const char* /*name*/, Int_t /*histoclass*/) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}
	return 0x0;

}

// -------------------------------------------------------------------------
TH3F* AliAODTrackCutsDiHadronPID::InitializePIDHisto(const char* name, Int_t histoclass, Int_t expspecies, Int_t ptclass) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	TH3F* hout = new TH3F(Form("%s%s%s%s",name,fHistoName[histoclass].Data(),fParticleName[expspecies].Data(),fPtClassName[ptclass].Data()),
		Form("PID %s (Exp: %s);p_{T} (GeV/c);#Delta t (ms);dE/dx (a.u.)",fHistoName[histoclass].Data(),fParticleName[expspecies].Data()),
		fNPtBinsPID[ptclass],fPtBoundaryPID[ptclass],fPtBoundaryPID[ptclass+1],
		fTOFbins[ptclass][expspecies],fTOFLowerBound[ptclass][expspecies],fTOFUpperBound[ptclass][expspecies],
		fTPCbins[ptclass][expspecies],fTPCLowerBound[ptclass][expspecies],fTPCUpperBound[ptclass][expspecies]);

	hout->SetDirectory(0);

	return hout;

}

// -------------------------------------------------------------------------
TH2F* AliAODTrackCutsDiHadronPID::InitializeTOFMismatchHisto(const char* name, Int_t histoclass, Int_t expspecies, Int_t ptclass) {

	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	TH2F* hout = new TH2F(Form("%s%s%s%s",name,fHistoName[histoclass].Data(),fParticleName[expspecies].Data(),fPtClassName[ptclass].Data()),
		Form("TOF Mismatch %s (Exp: %s);p_{T} (GeV/c);#Delta t (ms)",fHistoName[histoclass].Data(),fParticleName[expspecies].Data()),
		fNPtBinsPID[ptclass],fPtBoundaryPID[ptclass],fPtBoundaryPID[ptclass+1],
		fTOFbins[ptclass][expspecies],fTOFLowerBound[ptclass][expspecies],fTOFUpperBound[ptclass][expspecies]);

	hout->SetDirectory(0);

	return hout;

}

// -------------------------------------------------------------------------
TH1F* AliAODTrackCutsDiHadronPID::GetHistDataTOFProjection(Int_t charge, Int_t species, Int_t ptbin) {

	// Returns a projection in TOF of the data histogram.
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	Int_t ptclass = GetPtClass(ptbin);
	if (ptclass == -1) return 0x0;
	Int_t bininptclass = GetBinInPtClass(ptbin);

	// Retrieve original the 3D histogram.
	TH3F* htmp = (TH3F*)GetHistData(Form("fHistDataPID%s%s%s",fHistoName[charge].Data(),fParticleName[species].Data(),fPtClassName[ptclass].Data()));

	// Make the projection.
	TH1F* htmp_proj = (TH1F*)htmp->ProjectionY(Form("TOFprojection_%i_%i_%i",charge,species,ptbin),bininptclass,bininptclass);

	// Some settings of the output histogram.
	htmp_proj->SetDirectory(0);
	htmp_proj->SetTitle(Form("%5.3f < p_{T} < %5.3f",GetPtMinPID(ptbin),GetPtMaxPID(ptbin)));
	htmp_proj->Sumw2();

	return htmp_proj;

}

// -------------------------------------------------------------------------
TH1F* AliAODTrackCutsDiHadronPID::GetHistDataTOFMismatch(Int_t charge, Int_t species, Int_t ptbin) {

	// Returns a projection in TOF of the mismatch histogram.
	if (fDebug) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	Int_t ptclass = GetPtClass(ptbin);
	if (ptclass == -1) return 0x0;
	Int_t bininptclass = GetBinInPtClass(ptbin);

	// Retrieve the original 2D histogram.
	TH2F* htmp = (TH2F*)GetHistData(Form("fHistTOFMismatch%s%s%s",fHistoName[charge].Data(),fParticleName[species].Data(),fPtClassName[ptclass].Data()));

	// Make the projection.
	TH1F* htmp_proj = (TH1F*)htmp->ProjectionY(Form("TOFprojectionsMismatch_%i_%i_%i",charge,species,ptbin),bininptclass,bininptclass);

	// Some settings of the output histogram.
	htmp_proj->SetDirectory(0);
	htmp_proj->SetTitle(Form("%5.3f < p_{T} < %5.3f",GetPtMinPID(ptbin),GetPtMaxPID(ptbin)));
	htmp_proj->Sumw2();

	return htmp_proj;

}
