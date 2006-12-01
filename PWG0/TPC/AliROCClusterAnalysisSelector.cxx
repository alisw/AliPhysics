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

/* $Id$ */

// 
//  This class analyses TPC cosmics data from clusters
//
//  Authors: Jan.Fiete.Grosse-Oetringhaus@cern.ch, Claus.Jorgensen@cern.ch
//

#include "AliROCClusterAnalysisSelector.h"

#include <AliLog.h>
#include <AliTPCclusterMI.h>
#include <AliRunLoader.h>

#include <AliTPCClustersRow.h>
#include <AliESD.h>


#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TTimeStamp.h>
#include <TRandom.h>

#include "TPC/AliTPCClusterHistograms.h"

extern TSystem* gSystem;

ClassImp(AliROCClusterAnalysisSelector)

AliROCClusterAnalysisSelector::AliROCClusterAnalysisSelector() :
  AliSelectorRL(),
  fObjectsToSave(0)
{
  //
  // Constructor. Initialization of pointers
  //

  fNMaxObjectsToSave = 50;
  fObjectsToSave     = new TObjArray();
  
  for (Int_t i=0; i<kTPCHists; i++)
    fClusterHistograms[i] = 0;
}

AliROCClusterAnalysisSelector::~AliROCClusterAnalysisSelector()
{
  //
  // Destructor
  //
}

void AliROCClusterAnalysisSelector::SlaveBegin(TTree* tree)
{
  //
  
  AliSelectorRL::SlaveBegin(tree);
} 

void AliROCClusterAnalysisSelector::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  AliSelectorRL::Init(tree);
  
  // Set branch address
  if (tree) {  
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fTimeStamp", 1);
  }
}

Bool_t AliROCClusterAnalysisSelector::Process(Long64_t entry)
{
  //
  // Implement your analysis here. Do not forget to call the parent class Process by
  // if (AliSelectorRL::Process(entry) == kFALSE)
  //   return kFALSE;
  //

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;


  // reset counters
  for (Int_t i=0; i<kTPCHists; i++)
    if (fClusterHistograms[i])
      fClusterHistograms[i]->StartEvent();
  
  //  runLoader->Dump();

  Int_t flag = ProcessEvent(entry, kFALSE);
  if (flag != 0)
    ProcessEvent(entry, kTRUE, "");


  Int_t time = 0;  
  if (fESD) 
    if (fESD->GetTimeStamp()>1160000000) {
      time = fESD->GetTimeStamp();      
    }  
  
  // finish event
  for (Int_t i=0; i<kTPCHists; i++)
    if (fClusterHistograms[i])
      fClusterHistograms[i]->FinishEvent(time);

  
  // TODO This should not be needed, the TTree::GetEntry() should take care of this, maybe because it has a reference member, to be analyzed
  // if the ESDfriend is not deleted we get a major memory leak
  // here the esdfriend seems to be also deleted, very weird behaviour....

  delete fESD;
  fESD = 0;    

  return kTRUE;
}

Int_t AliROCClusterAnalysisSelector::ProcessEvent(Long64_t entry, Bool_t detailedHistogram, const Char_t* label)
{
  //
  // Looping over clusters in event and filling histograms 
  //
  // - if detailedHistogram = kTRUE special histograms are saved (in fObjectsToSave)
  //   

  // save maximum N objects
  if (detailedHistogram) 
    if (fObjectsToSave->GetEntries() > fNMaxObjectsToSave) 
      return 0;
      
  // TODO: find a clever way to handle the time      
  Int_t time = 0;  
  if (fESD) 
    if (fESD->GetTimeStamp()>1160000000) {
      time = fESD->GetTimeStamp();      
    }  

  // for saving single events
  AliTPCClusterHistograms* clusterHistograms[kTPCSectors];
  Bool_t keepEvent[kTPCSectors];
  for (Int_t i=0; i<kTPCSectors; i++) { 
    clusterHistograms[i] = 0;  
    keepEvent[i] = kFALSE;

    if (fClusterHistograms[i]) {
      TString why;
      keepEvent[i] = fClusterHistograms[i]->KeepThisEvent(why);
    }
  }
  
  Bool_t intToReturn = 0;

  // --------------

  AliRunLoader* runLoader = GetRunLoader();
  runLoader->LoadRecPoints("TPC");
  
  if (!runLoader) {
    AliDebug(AliLog::kError, " Run loader not found");
    return kFALSE;
  }
  
  TTree* tree = runLoader->GetTreeR("TPC", kFALSE);
  
  // load clusters to the memory
  AliTPCClustersRow* clrow = new AliTPCClustersRow;
  clrow->SetClass("AliTPCclusterMI");
  clrow->SetArray(0);
  clrow->GetArray()->ExpandCreateFast(10000);
  //

  if (!tree) {
    AliDebug(AliLog::kError, " TPC cluster tree not found");
    return kFALSE;
  }

  TBranch* br = tree->GetBranch("Segment");
  br->SetAddress(&clrow);
  //

  Int_t j = Int_t(tree->GetEntries());
  for (Int_t i=0; i<j; i++) {
    br->GetEntry(i);
    //  
    for (Int_t icl=0; icl<clrow->GetArray()->GetEntriesFast(); icl++){
      
      AliTPCclusterMI* cluster = (AliTPCclusterMI*)clrow->GetArray()->At(icl);

      if (!cluster) 
	continue;
      
     Int_t detector = cluster->GetDetector();
     
     if (detector < 0 || detector >= kTPCSectors) 
     {
       AliDebug(AliLog::kDebug, Form("We found a cluster from invalid sector %d", detector));
       continue;
     }

     if (!detailedHistogram) {
       
       if (!fClusterHistograms[detector])
	 fClusterHistograms[detector] = new AliTPCClusterHistograms(detector,"",time,time+5*60*60);
       
   	if (!fClusterHistograms[detector+kTPCSectors])
   	  fClusterHistograms[detector+kTPCSectors] = new AliTPCClusterHistograms(detector,"",time,time+5*60*60, kTRUE);
	
   	fClusterHistograms[detector]->FillCluster(cluster, time);
   	fClusterHistograms[detector+kTPCSectors]->FillCluster(cluster, time);
			
     } // end of if !detailedHistograms
     else {
       // if we need the detailed histograms for this event
       
       if (keepEvent[detector]) {	 

	 TString why(fClusterHistograms[detector]->WhyKeepEvent());
	
	 why.Append(Form("_entry_%d",entry));
	 why.Append(label);

	 // if clusterHistograms for this event is not there, construct it
	 if (!clusterHistograms[detector]) {
	   clusterHistograms[detector] = new AliTPCClusterHistograms(detector, why.Data());
	   
	   // adding file and time as comment
	   TString comment = TString(Form("%s",fTree->GetCurrentFile()->GetName()));
	   comment.Append(Form(" entry %d", entry));
	   if (time!=0) {
	     TString timeStr(TTimeStamp(time).AsString());
	     timeStr.Remove(26);
	     
	     comment.Append(Form(" (%s)",timeStr.Data()));
	   }  
	   clusterHistograms[detector]->SetCommentToHistograms(comment.Data());
	   
	 }  
	 else 
	   clusterHistograms[detector]->FillCluster(cluster);
       } // end of (keep this event)
     } // 
    }
  }
  
  if (!detailedHistogram) {
    for (Int_t i=0; i<kTPCSectors; i++)
      if (fClusterHistograms[i]) {
	TString why;
	if (fClusterHistograms[i]->KeepThisEvent(why))
 	  intToReturn = 1;
      } 
  }
  else {
    for (Int_t i=0; i< kTPCSectors; i++) {
      if (clusterHistograms[i]) {
	fObjectsToSave->Add(clusterHistograms[i]);
      }
    }    
  }  
  
  delete clrow;
   
  return intToReturn;
}


void AliROCClusterAnalysisSelector::SlaveTerminate()
{
  //
  
  if (fOutput)
  {
    for (Int_t i=0; i<kTPCHists; i++)
      if (fClusterHistograms[i])
        fOutput->Add(fClusterHistograms[i]);
  }
} 

void AliROCClusterAnalysisSelector::Terminate()
{
  // 
  // read the objects from the output list and write them to a file
  // the filename is modified by the object comment passed in the tree info or input list
  //

  if (fOutput)
  {  
    fOutput->Print();
        
    for (Int_t i=0; i<kTPCSectors; i++)
      fClusterHistograms[i] = dynamic_cast<AliTPCClusterHistograms*> (fOutput->FindObject(AliTPCClusterHistograms::FormDetectorName(i, kFALSE)));
    for (Int_t i=0; i<kTPCSectors; i++)
      fClusterHistograms[kTPCSectors+i] = dynamic_cast<AliTPCClusterHistograms*> (fOutput->FindObject(AliTPCClusterHistograms::FormDetectorName(i, kTRUE)));
  }
  
  TNamed* comment = 0;
  if (fTree && fTree->GetUserInfo())
    comment = dynamic_cast<TNamed*>(fTree->GetUserInfo()->FindObject("comment"));
  if (!comment && fInput)
    comment = dynamic_cast<TNamed*>(fInput->FindObject("comment"));

  if (comment)
  {
    AliDebug(AliLog::kInfo, Form("INFO: Found comment in input list: %s", comment->GetTitle()));
  }
  else
    return;

  TFile* file = TFile::Open(Form("rocCluster_%s.root",comment->GetTitle()), "RECREATE");

  for (Int_t i=0; i<kTPCHists; i++)
     if (fClusterHistograms[i]) {
       fClusterHistograms[i]->SaveHistograms();

//        TCanvas* c = fClusterHistograms[i]->DrawHistograms();
//        TString dir;
//        dir.Form("WWW/%s/%s", comment->GetTitle(), c->GetName());
//        gSystem->mkdir(dir, kTRUE);
//        c->SaveAs(Form("%s/plots_%s_%s.eps",dir.Data(),comment->GetTitle(),c->GetName()));
//        c->SaveAs(Form("%s/plots_%s_%s.gif",dir.Data(),comment->GetTitle(),c->GetName()));
    
//        c->Close();
//        delete c;
     }
  
  gDirectory->mkdir("saved_objects");
  gDirectory->cd("saved_objects");

  for (Int_t i=0; i<fObjectsToSave->GetSize(); i++) {
    if (fObjectsToSave->At(i)) {
      AliTPCClusterHistograms* clusterHistograms = dynamic_cast<AliTPCClusterHistograms*> (fObjectsToSave->At(i));
      if (clusterHistograms)
	clusterHistograms->SaveHistograms();
      else
	fObjectsToSave->At(i)->Write();
    }
  }

  gDirectory->cd("../");


  file->Close();
}
