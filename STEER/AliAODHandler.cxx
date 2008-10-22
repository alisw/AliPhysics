/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     Implementation of the Virtual Event Handler Interface for AOD
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------


#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TList.h>

#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODTracklets.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"



ClassImp(AliAODHandler)

//______________________________________________________________________________
AliAODHandler::AliAODHandler() :
    AliVEventHandler(),
    fIsStandard(kTRUE),
    fFillAOD(kTRUE),
    fNeedsHeaderReplication(kFALSE),
    fNeedsTracksBranchReplication(kFALSE),
    fNeedsVerticesBranchReplication(kFALSE),
    fNeedsV0sBranchReplication(kFALSE),
    fNeedsTrackletsBranchReplication(kFALSE),
    fNeedsPMDClustersBranchReplication(kFALSE),
    fNeedsJetsBranchReplication(kFALSE),
    fNeedsFMDClustersBranchReplication(kFALSE),
    fNeedsCaloClustersBranchReplication(kFALSE),
    fAODIsReplicated(kFALSE),
    fAODEvent(NULL),
    fMCEventH(NULL),
    fTreeA(NULL),
    fFileA(NULL),
    fFileName("")
{
  // default constructor
}

//______________________________________________________________________________
AliAODHandler::AliAODHandler(const char* name, const char* title):
    AliVEventHandler(name, title),
    fIsStandard(kTRUE),
    fFillAOD(kTRUE),
    fNeedsHeaderReplication(kFALSE),
    fNeedsTracksBranchReplication(kFALSE),
    fNeedsVerticesBranchReplication(kFALSE),
    fNeedsV0sBranchReplication(kFALSE),
    fNeedsTrackletsBranchReplication(kFALSE),
    fNeedsPMDClustersBranchReplication(kFALSE),
    fNeedsJetsBranchReplication(kFALSE),
    fNeedsFMDClustersBranchReplication(kFALSE),
    fNeedsCaloClustersBranchReplication(kFALSE),
    fAODIsReplicated(kFALSE),
    fAODEvent(NULL),
    fMCEventH(NULL),
    fTreeA(NULL),
    fFileA(NULL),
    fFileName("")
{
}

//______________________________________________________________________________
AliAODHandler::~AliAODHandler() 
{
  delete fAODEvent;
  if(fFileA){
    // is already handled in TerminateIO
    fFileA->Close();
    delete fFileA;
  }
  delete fTreeA;
 // destructor
}

//______________________________________________________________________________
Bool_t AliAODHandler::Init(Option_t* opt)
{
  // Initialize IO
  //
  // Create the AODevent object
  if(!fAODEvent){
    fAODEvent = new AliAODEvent();
    if (fIsStandard) fAODEvent->CreateStdContent();
  }
  //
  // File opening according to execution mode
  TString option(opt);
  option.ToLower();
  if (option.Contains("proof")) {
    // proof
    if (option.Contains("special")) {
       // File for tree already opened on slave -> merging via files
       fFileA = gFile;
       CreateTree(1);
    } else {   
       // Merging in memory
       CreateTree(0);
    }   
  } else {
    // local and grid
    TDirectory *owd = gDirectory;
    fFileA = new TFile(fFileName.Data(), "RECREATE");
    CreateTree(1);
    owd->cd();
  }
  return kTRUE;
}


void AliAODHandler::StoreMCParticles(){

  // 
  // Remap the labels from ESD stack and store
  // the AODMCParticles, makes only sense if we have
  // the mcparticles branch
  // has to be done here since we cannot know in advance 
  // which particles are needed (e.g. by the tracks etc.)
  //
  // Particles have been selected by AliMCEventhanlder->SelectParticle()
  // To use the MCEventhandler here we need to set it from the outside
  // can vanish when Handler go to the ANALYSISalice library
  //
  // The Branch booking for mcParticles and mcHeader has to happen 
  // in an external task for now since the AODHandler does not have access
  // the AnalysisManager. For the same reason the pointer t o the MCEventH
  // has to passed to the AOD Handler by this task 
  // (doing this in the steering macro would not work on PROOF)

  TClonesArray *mcarray = (TClonesArray*)fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()); 
  if(!mcarray)return;
  mcarray->Delete();

  AliAODMCHeader *mcHeader = (AliAODMCHeader*)fAODEvent->FindListObject(AliAODMCHeader::StdBranchName()); 
  if(!mcHeader)return;

  mcHeader->Reset();

  // Get the MC Infos.. Handler needs to be set before 
  // while adding the branch
  // This needs to be done, not to depend on the AnalysisManager

  if(!fMCEventH)return;
  if(!fMCEventH->MCEvent())return;
  AliStack *pStack = fMCEventH->MCEvent()->Stack();
  if(!pStack)return;

  fMCEventH->CreateLabelMap();

  //
  // Get the Event Header 
  // 

  AliHeader* header = fMCEventH->MCEvent()->Header();
  if (!header)return;

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);
  mcHeader->SetVertex(vtxMC[0],vtxMC[1],vtxMC[2]);

  // we search the MCEventHeaders first 
  // Two cases, cocktail or not...
  AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
  if(genCocktailHeader){
    // we have a coktail header
    mcHeader->AddGeneratorName(genHeader->GetName());
    // Loop from the back so that the first one sets the process type
    TList* headerList = genCocktailHeader->GetHeaders();
    for(int i = headerList->GetEntries()-1;i>=0;--i){
      AliGenEventHeader *headerEntry = dynamic_cast<AliGenEventHeader*>(headerList->At(i));
      SetMCHeaderInfo(mcHeader,headerEntry);
    }
  }
  else{
    // No Cocktail just take the first one
    SetMCHeaderInfo(mcHeader,genHeader);
  }




  // Store the AliAODParticlesMC

  Int_t np    = pStack->GetNtrack();
  Int_t nprim = pStack->GetNprimary();


  Int_t j = 0;
  TClonesArray& l = *mcarray;

  for(int i = 0;i < np;++i){
    if(fMCEventH->IsParticleSelected(i)){

      Int_t flag = 0;
      TParticle *part = pStack->Particle(i);
      if(i<nprim)flag |= AliAODMCParticle::kPrimary;
      if(pStack->IsPhysicalPrimary(i))flag |= AliAODMCParticle::kPhysicalPrim;

      if(fMCEventH->GetNewLabel(i)!=j){
	AliError(Form("MISMATCH New label %d j: %d",fMCEventH->GetNewLabel(i),j));
      }
      AliAODMCParticle mcpart_tmp(part,i,flag);

      // 
      Int_t d0 =  mcpart_tmp.GetDaughter(0);
      Int_t d1 =  mcpart_tmp.GetDaughter(1);
      Int_t m =  mcpart_tmp.GetMother();

      // other than for the track labels, negative values mean
      // no daughter/mother so preserve it
 
      if(d0<0 && d1<0){
	// no first daughter -> no second daughter
	// nothing to be done
	// second condition not needed just for sanity check at the end
	mcpart_tmp.SetDaughter(0,d0);
	mcpart_tmp.SetDaughter(1,d1);
      }
      else if(d1 < 0 && d0 >= 0){
	// Only one daughter
	// second condition not needed just for sanity check at the end
	if(fMCEventH->IsParticleSelected(d0)){
	  mcpart_tmp.SetDaughter(0,fMCEventH->GetNewLabel(d0));
	}
	else{
	  mcpart_tmp.SetDaughter(0,-1);
	}
	mcpart_tmp.SetDaughter(1,d1);
      }
      else if (d0 > 0 && d1 > 0 ){
	// we have two or more daughters loop on the stack to see if they are
	// selected
	Int_t d0_tmp = -1;
	Int_t d1_tmp = -1;
	for(int id = d0; id<=d1;++id){
	  if(fMCEventH->IsParticleSelected(id)){
	    if(d0_tmp==-1){
	      // first time
	      d0_tmp = fMCEventH->GetNewLabel(id);
	      d1_tmp = d0_tmp; // this is to have the same schema as on the stack i.e. with one daugther d0 and d1 are the same 
	    }
	    else d1_tmp = fMCEventH->GetNewLabel(id);
	  }
	}
	mcpart_tmp.SetDaughter(0,d0_tmp);
	mcpart_tmp.SetDaughter(1,d1_tmp);
      }
      else{
	AliError(Form("Unxpected indices %d %d",d0,d1));
      }

      if(m<0){
	mcpart_tmp.SetMother(m);
      }
      else{
	if(fMCEventH->IsParticleSelected(m))mcpart_tmp.SetMother(fMCEventH->GetNewLabel(m));
	else AliError("PROBLEM Mother not selected");
      }

      new (l[j++]) AliAODMCParticle(mcpart_tmp);
      
    }
  }
  AliInfo(Form("AliAODHandler::StoreMCParticles: Selected %d (Primaries %d / total %d) after validation",
	       j,nprim,np));
  
  // Set the labels in the AOD output...
  // Remapping

  // AODTracks
  TClonesArray* tracks = fAODEvent->GetTracks();
  if(tracks){
    for(int it = 0; it < fAODEvent->GetNTracks();++it){
      AliAODTrack *track = fAODEvent->GetTrack(it);
      
      if(TMath::Abs(track->GetLabel())>np||track->GetLabel()==0){
	AliWarning(Form("Wrong ESD track label %d",track->GetLabel()));
      }
      if(fMCEventH->GetNewLabel(track->GetLabel())==0){
	AliWarning(Form("New label not found for %d",track->GetLabel()));
      }
      track->SetLabel(fMCEventH->GetNewLabel(track->GetLabel()));
    }
  }
  
  // AOD calo cluster
  TClonesArray *clusters = fAODEvent->GetCaloClusters();
  if(clusters){
    for (Int_t iClust = 0;iClust < fAODEvent->GetNCaloClusters(); ++iClust) {
      AliAODCaloCluster * cluster = fAODEvent->GetCaloCluster(iClust);
      UInt_t nLabel    = cluster->GetNLabel();
      // Ugly but do not want to fragment memory by creating 
      // new Int_t (nLabel)
      Int_t* labels    = const_cast<Int_t*>(cluster->GetLabels());
      if (labels){
	for(UInt_t i = 0;i < nLabel;++i){
	  labels[i] = fMCEventH->GetNewLabel(cluster->GetLabel(i));
	}
      }
      //      cluster->SetLabels(labels,nLabel);
    }// iClust
  }// clusters

  // AOD tracklets
  AliAODTracklets *tracklets = fAODEvent->GetTracklets();
  if(tracklets){
    for(int it = 0;it < tracklets->GetNumberOfTracklets();++it){
      int label0 = tracklets->GetLabel(it,0);
      int label1 = tracklets->GetLabel(it,1);
      if(label0>=0)label0 = fMCEventH->GetNewLabel(label0);      
      if(label1>=0)label1 = fMCEventH->GetNewLabel(label1);
      tracklets->SetLabel(it,0,label0);
      tracklets->SetLabel(it,1,label1);
    }
  }

}

Bool_t AliAODHandler::FinishEvent()
{
  // Fill data structures
  if(fFillAOD){
    fAODEvent->MakeEntriesReferencable();
    StoreMCParticles();
    FillTree();
  }

  if (fIsStandard) fAODEvent->ResetStd();
  // Reset AOD replication flag
  fAODIsReplicated = kFALSE;
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODHandler::Terminate()
{
    // Terminate 
    AddAODtoTreeUserInfo();
    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODHandler::TerminateIO()
{
    // Terminate IO
    if (fFileA) {
	fFileA->Close();
	delete fFileA;
    }
    return kTRUE;
}

//______________________________________________________________________________
void AliAODHandler::CreateTree(Int_t flag)
{
    // Creates the AOD Tree
    fTreeA = new TTree("aodTree", "AliAOD tree");
    fTreeA->Branch(fAODEvent->GetList());
    if (flag == 0) fTreeA->SetDirectory(0);
}

//______________________________________________________________________________
void AliAODHandler::FillTree()
{
    // Fill the AOD Tree
  fTreeA->Fill();
}

//______________________________________________________________________________
void AliAODHandler::AddAODtoTreeUserInfo()
{
  // Add aod event to tree user info
  fTreeA->GetUserInfo()->Add(fAODEvent);
}

//______________________________________________________________________________
void AliAODHandler::AddBranch(const char* cname, void* addobj)
{
    // Add a new branch to the aod 
    TDirectory *owd = gDirectory;
    if (fFileA) {
      fFileA->cd();
    }
    char** apointer = (char**) addobj;
    TObject* obj = (TObject*) *apointer;

    fAODEvent->AddObject(obj);
 
    const Int_t kSplitlevel = 99; // default value in TTree::Branch()
    const Int_t kBufsize = 32000; // default value in TTree::Branch()

    if (!fTreeA->FindBranch(obj->GetName())) {
      // Do the same as if we book via 
      // TTree::Branch(TCollection*)
      
      fTreeA->Bronch(obj->GetName(), cname, fAODEvent->GetList()->GetObjectRef(obj),
		     kBufsize, kSplitlevel - 1);
      //    fTreeA->Branch(obj->GetName(), cname, addobj);
    }
    owd->cd();
}

//______________________________________________________________________________
void AliAODHandler::SetOutputFileName(const char* fname)
{
// Set file name.
   fFileName = fname;
}

//______________________________________________________________________________
const char *AliAODHandler::GetOutputFileName()
{
// Get file name.
   return fFileName.Data();
}

void  AliAODHandler::SetMCHeaderInfo(AliAODMCHeader *mcHeader,AliGenEventHeader *genHeader){


  // Utility function to cover different cases for the AliGenEventHeader
  // Needed since different ProcessType and ImpactParamter are not 
  // in the base class...
  // We don't encode process types for event cocktails yet
  // coul be done e.g. by adding offsets depnding on the generator


  mcHeader->AddGeneratorName(genHeader->GetName());



  if(!genHeader)return;
  AliGenPythiaEventHeader *pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
    if (pythiaGenHeader) {
      mcHeader->SetEventType(pythiaGenHeader->ProcessType());
      return;
    }
    
    AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);

  if (dpmJetGenHeader){
    mcHeader->SetEventType(dpmJetGenHeader->ProcessType());
    return;
  } 

  AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
  if(hijingGenHeader){
    mcHeader->SetImpactParameter(hijingGenHeader->ImpactParameter());
    return;
  }
  
  AliWarning(Form("MC Eventheader not known: %s",genHeader->GetName()));

}
