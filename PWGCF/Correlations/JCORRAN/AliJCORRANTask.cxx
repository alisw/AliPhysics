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

//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects  
//////////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TList.h"
#include "TTree.h"
#include "TFile.h"


#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"

#include "AliJCORRANTask.h" 
#include "AliAnalysisManager.h"

#include "AliJTrack.h"
#include "AliJMCTrack.h"
#include "AliJPhoton.h"
//#include "AliJCaloCell.h"
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"

//______________________________________________________________________________
AliJCORRANTask::AliJCORRANTask() :   
    AliAnalysisTaskSE("PWG4JCORRAN"),
    fFilter(0x0),
    fAODName("jcorran.root"),
    fJODTree(0x0),
    fAliJRunHeader(0x0),
    fDoStoreJOD(kFALSE)

{

  DefineInput (0, TChain::Class());
  DefineOutput (1, TTree::Class());
  DefineOutput (2, TList::Class());
  
   fFilter = new AliJFilter();
}

//______________________________________________________________________________
AliJCORRANTask::AliJCORRANTask(const char *name, TString inputformat):
    AliAnalysisTaskSE(name), 
    fFilter(0x0),
    fAODName("jcorran.root"),
    fJODTree(0x0),
    fAliJRunHeader(0x0),
    fDoStoreJOD(kFALSE)
{
  // Constructor
  AliInfo("---- AliJCORRANTask Constructor ----");

  JUNUSED(inputformat);

  DefineInput (0, TChain::Class());
  DefineOutput (1, TTree::Class());
  DefineOutput (2, TList::Class());

   fFilter = new AliJFilter( Form("%sFilter",name), this );
}

//____________________________________________________________________________
AliJCORRANTask::AliJCORRANTask(const AliJCORRANTask& ap) :
    AliAnalysisTaskSE(ap.GetName()), 
    fFilter(ap.fFilter),
    fAODName(ap.fAODName),
    fJODTree(ap.fJODTree),
    fAliJRunHeader(ap.fAliJRunHeader),
    fDoStoreJOD(ap.fDoStoreJOD)
{ 

  AliInfo("----DEBUG AliJCORRANTask COPY ----");

}

//_____________________________________________________________________________
AliJCORRANTask& AliJCORRANTask::operator = (const AliJCORRANTask& ap)
{
  // assignment operator

  AliInfo("----DEBUG AliJCORRANTask operator= ----");
  this->~AliJCORRANTask();
  new(this) AliJCORRANTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJCORRANTask::~AliJCORRANTask()
{
  // destructor 

   delete fFilter;
   delete fJODTree;
   delete fAliJRunHeader;

}

//________________________________________________________________________

void AliJCORRANTask::UserCreateOutputObjects()
{  
  //=== create the jcorran outputs objects
  if(fDebug > 1) printf("AliJCORRANTask::UserCreateOutPutData() \n");
  
  //=== Get AnalysisManager
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if(!man->GetOutputEventHandler()) {
    Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
    return;
  }

   // run the filter class
   fFilter->SetMyTask( this );
   fFilter->SetAliJRunHeader( fAliJRunHeader );
   fFilter->UserCreateOutputObjects();

  // register ouput branches

  if(fDoStoreJOD){
		TFile * file1 = OpenFile(1);
		file1-> SetCompressionLevel(9);
		TTree * tree = new TTree("JODTree","JYFL Object Data");
		int split = 2;
		int basketsize = 32000;
		tree->Branch("TrackList", fFilter->GetTrackList(),basketsize, split);
		if( fFilter->IsMC() ) 
				tree->Branch("MCTrackList", fFilter->GetMCTrackList(),basketsize, split );
		//== Event Header
		tree->Branch("HeaderList", fFilter->GetHeaderList(),basketsize, split );
		//== EventPlane SRC
		if( fFilter->GetStoreEventPlaneSource() ){
				tree->Branch("AliESDVZERO", fFilter->GetESDVZERO());
				tree->Branch("AliESDTZERO", fFilter->GetESDTZERO());
				tree->Branch("AliESDZDC",   fFilter->GetESDZDC());
		}
		fJODTree = tree;
  }

  PostData( 1, fJODTree );
  //OpenFile(2);
  PostData( 2,fFilter->GetRunInfoList());


  cout << "Add(fAliRunHeader) in UserCreateObject() ======= " << endl;

}

//______________________________________________________________________________
void AliJCORRANTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) cout << "------- AliJCORRANTask Exec-------"<<endl;
	if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry())); 

	fFilter->UserExec("");

	if(  1 || fFilter->GetEventSuccess() ){ // TODO
		if( fDoStoreJOD ){
			fJODTree->Fill();
		}
	}
	PostData(1,fJODTree);
	PostData(2,fFilter->GetRunInfoList());


	if(fDebug > 5) cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJCORRANTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 

	fFilter->Init();

	//   TString formula(fEsdTrackCuts->GetMaxDCAToVertexXYPtDep());
	//   if(formula.Length()>0){ // momentum dep DCA cut for AOD
	//     formula.ReplaceAll("pt","x");
	//   }
}

//______________________________________________________________________________
void AliJCORRANTask::Terminate(Option_t * option)
{
	fFilter->Terminate();

	// Processing when the event loop is ended
	fAliJRunHeader->PrintOut();
	cout<<"AliJCORRANTask Analysis DONE !!"<<endl; 
	// Printout fRunInfoList here
	TList* fRunInfoList = dynamic_cast<TList*> (GetOutputData(1));
	if(fRunInfoList)
	{
		AliJRunHeader *fAliRunHeader = dynamic_cast<AliJRunHeader*> (fRunInfoList->FindObject("AliJRunHeader"));
		if(fAliRunHeader) {fAliRunHeader->Print();}
	}
	else
	{
		cout << "WARNING : Run Information List is empty" << endl;
	}

}
