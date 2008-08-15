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
 

/*
  Realisation of an AliVEventPool based on THnSparse
  Author Peter Hristov
  Peter.Hristov@cern.ch

  Possible usage: three steps
  1) Creation of a XML tag collection in aliensh

  aliensh:> find -x charm /alice/sim/PDC_08/LHC08x/180001 tag.root > 180001.xml

  2) Merging the tag files

  TGrid::Connect("alien://");
  TAlienCollection *collection = TAlienCollection::Open(xmlfile);
  TGridResult* result = collection->GetGridResult("",0);
  AliTagCreator *t = new AliTagCreator();
  t->MergeTags("ESD",result);

  3) Chain the merged tag files and test the event pool

void testpool(const char * dirname = ".", const char * pattern = "Run180001") {

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");

  // Create a chain
  TChain * fChain = new TChain("T");


  // Chain the tag files in the working directory
  TString fTagFilename;
  
  // Open the working directory
  void * dirp = gSystem->OpenDirectory(dirname);
  const char * name = 0x0;
  // Add all files matching *pattern* to the chain
  while((name = gSystem->GetDirEntry(dirp))) {
    cout << name << endl;
    if (strstr(name,pattern)) { 
      fTagFilename = dirname;
      fTagFilename += "/";
      fTagFilename += name;
	  	
      fChain->Add(fTagFilename);  
    }//pattern check
  }//directory loop


  Int_t nruns = fChain->GetEntries();

  cout << nruns << " run(s) found in the tag chain." << endl;

  Int_t dim = 3;
  const char * vars[] = {"fNumberOfPositiveTracks","fNumberOfNegativeTracks","fPrimaryVertexZ"};
  Int_t nbins[] = {10,10,10};
  Double_t xmin[] ={-0.5,-0.5,-20};
  Double_t xmax[] ={49/5,49.5,20};
  Int_t chunksize = 100;

  AliEventPoolSparse * pool = 
    new AliEventPoolSparse("test", "test", fChain, dim, vars, nbins, xmin, xmax, chunksize);

  pool->Init();

  TChain * esdchain = 0x0;
  Int_t ichain = 0;
  while (esdchain=pool->GetNextChain()) {
    cout << "Chain: "<< ichain <<" Events: " << esdchain->GetEntries() << endl;
    ichain++;
  }

  delete fChain;

}

*/

#include "AliEventPoolSparse.h"
#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliLog.h"

#include <TObjArray.h>
#include <TAxis.h>
#include <TTreeFormula.h>
#include <TChain.h>
#include <TFile.h>
#include <Riostream.h>
#include <cstring>

ClassImp(AliEventPoolSparse)

// _________________________________________________________________________
AliEventPoolSparse::AliEventPoolSparse() :
  fHnSparseI(),
  fChunkSize(1024 * 16),
  fN(0),
  fPool(0x0),
  fCurrentBin(-1),
  fTagChain(0x0),
  fVars(0x0),
  fRunCut(0x0),
  fLHCCut(0x0),
  fDetCut(0x0),
  fEvCut(0x0),
  fBinNumber(0)
{
  // Default constructor. Initializes the THnSparseI,
  // the initial size of the array and the array itself
  fN = fChunkSize;
  fPool = new TEntryList * [fN];
  memset(fPool,0x0,fN*sizeof(TEntryList*));
}

// _________________________________________________________________________
AliEventPoolSparse::AliEventPoolSparse(const char* name, const char* title, TChain * tagchain, Int_t dim,
				       const char ** vars, const Int_t* nbins, const Double_t* xmin,
				       const Double_t* xmax, Int_t chunksize):
  fHnSparseI(name, title, dim, nbins, xmin, xmax, chunksize),
  fChunkSize(chunksize),
  fN(0),
  fPool(0x0),
  fCurrentBin(-1),
  fTagChain(tagchain),
  fVars(0x0),
  fRunCut(0x0),
  fLHCCut(0x0),
  fDetCut(0x0),
  fEvCut(0x0),
  fBinNumber(0){
  // Constructor. Initializes the THnSparseI,
  // the initial size of the pool array and the array itself
  // It uses the provided array of variables to create TTreeFormulas
  // that are used when the pools are filled. This is the reason to require the input
  // tag chain in the constructor.

  fN = fChunkSize;
  fPool = new TEntryList * [fN];
  memset(fPool,0x0,fN*sizeof(TArrayI*));

  // Pool variables
  fVars = new TTreeFormula*[dim];

  for (Int_t ivar=0; ivar<dim; ++ivar) {
    fVars[ivar] = new TTreeFormula(vars[ivar],vars[ivar],fTagChain);
  }


}

// _________________________________________________________________________
AliEventPoolSparse::~AliEventPoolSparse() {
  // Destructor. Delete the pool, the array of TTreeFormula
  // and the pointers to cuts
  for (Int_t i=0; i<fN; ++i) delete fPool[i];
  if (fN>0) delete [] fPool;

  Int_t ndim = fHnSparseI.GetNdimensions();
  for (Int_t i=0; i<ndim; ++i) delete fVars[i];
  delete [] fVars;

  delete fRunCut;
  delete fLHCCut;
  delete fDetCut;
  delete fEvCut;
}


// Implementation of the interface functions
// _________________________________________________________________________
TChain* AliEventPoolSparse::GetNextChain(){
  // Return the chains one by one. The output is 0x0 if the pool is not initialized
  // or the last chain is already reached
  if (fCurrentBin<0) {
    AliError("The event pool is not initialized");
    return 0x0;
  }

  if (fCurrentBin>=fHnSparseI.GetNbins()) { // Check if >= or >
    AliInfo("No more chains");
    return 0x0;
  }

  fBinNumber++;
  
  fChain->SetEntryList(fPool[fCurrentBin++],"ne"); 
  return fChain;
}

// _________________________________________________________________________
void  AliEventPoolSparse::GetCurrentBin(Float_t* xbin) {
  // This method fills the center of the current bin in xbin

  if (fCurrentBin<0) {
    AliError("The event pool is not initialized");
    return;
  }

  Int_t ndim = fHnSparseI.GetNdimensions();
  Int_t * coord = new Int_t[ndim]; 
  fHnSparseI.GetBinContent(fCurrentBin,coord);

  TObjArray * axes = fHnSparseI.GetListOfAxes();
  for (Int_t i=0; i<ndim; ++i) 
    xbin[i]=((TAxis*)axes->At(i+1))->GetBinCenter(coord[i]);

  delete [] coord;
}

// _________________________________________________________________________
void  AliEventPoolSparse::Init(){
  // Loop on the tag chain and select the events according
  // to the Run, LHC, detector, and event cuts.
  // Fill the THnSparse bin and add the event to the corresponding pool
  // Taken and modified from AliAnalysisTag

  if (!fTagChain) {
    AliError("Please provide a tag chain!");
    return;
  }

  Int_t ndim = fHnSparseI.GetNdimensions();
  if (ndim<=0) return;

  Double_t * x = new Double_t[ndim];
  
  // Tag objects.
  AliRunTag *tag = new AliRunTag; 	 
  AliEventTag *evTag = new AliEventTag;  
  fTagChain->SetBranchAddress("AliTAG",&tag); 	 
  
  TString guid(""); 	 
  TString turl(""); 	 
  TString path(""); 	 
  
  Int_t current = -1; 	 // Current tree number
  for(Int_t iTagFiles = 0; iTagFiles < fTagChain->GetEntries(); iTagFiles++) { 	 
    fTagChain->GetEntry(iTagFiles); 	 

    if (current != fTagChain->GetTreeNumber()) { 	 
      // Update the formula leaves if a new file is processed by the chain
      if (fRunCut) fRunCut->UpdateFormulaLeaves(); 	 
      if (fLHCCut) fLHCCut->UpdateFormulaLeaves(); 	 
      if (fDetCut) fDetCut->UpdateFormulaLeaves(); 	 
      if (fEvCut)  fEvCut->UpdateFormulaLeaves();

      for (Int_t ivar=0; ivar<fHnSparseI.GetNdimensions(); ++ivar)
	if (fVars[ivar]) fVars[ivar]->UpdateFormulaLeaves();

      // Create the ESD/AOD chain if not done
      if (!fChain) {
	// Decide if we have ESD or AOD
	TFile * tagfile = fTagChain->GetFile();
	if (strstr(tagfile->GetName(),"ESD")) fChain = new TChain("esdTree");
	else if (strstr(tagfile->GetName(),"AOD")) fChain = new TChain("aodTree");
	else {
	  AliError("Only ESD and AOD type is implemented!!!");
	  delete [] x;
	  return;
	}
      }
      
      // Update the tree number
      current = fTagChain->GetTreeNumber();
    }

    // Apply Run, LHC, and detector cuts if they exist
    if(!fRunCut || fRunCut->EvalInstance(iTagFiles) == 1) { 	 
      if(!fLHCCut || fLHCCut->EvalInstance(iTagFiles) == 1) { 	 
	if(!fDetCut || fDetCut->EvalInstance(iTagFiles) == 1) {
 	 

	  // Get access to the event data in the TTreeFormula
	  if (fEvCut) fEvCut->GetNdata();
	  for (Int_t ivar=0; ivar<fHnSparseI.GetNdimensions(); ++ivar)
	    if (fVars[ivar]) fVars[ivar]->GetNdata();
 	 
	  // Loop on events
	  const TClonesArray *tagList = tag->GetEventTags();
	  Int_t iEvents = tagList->GetEntries();
	  for(Int_t i = 0; i < iEvents; i++) { 	 
	    evTag = (AliEventTag *) tagList->At(i); 	 

	    guid = evTag->GetGUID(); 	 
	    turl = evTag->GetTURL(); 	 
	    path = evTag->GetPath(); 	 


	    if(!fEvCut || fEvCut->EvalInstance(i) == 1) {
	      TEntryList *fLocalList = new TEntryList();
	      fLocalList->SetTreeName(fChain->GetName());
	      fLocalList->SetFileName(turl.Data());
	      fLocalList->Enter(i);


	      // Add this event to the corresponding pool
	      {
		// Increment the bin content corrresponding to the vector "x" by "w",
		// and store the event index iev to the array associated with the bin,
		// then return the bin index.

		for (Int_t ivar=0; ivar<ndim; ++ivar) x[ivar] = fVars[ivar]->EvalInstance(i);
		
		Int_t bin =  fHnSparseI.Fill(x);
		// Check if we have to enlarge the array of pointers
		if (bin>=fN) Set(bin+fChunkSize);
		// Allocate the TEntryList if this is the first use of it
		if (!fPool[bin]) fPool[bin] = new TEntryList();
		// Add the event iev to the corresponding bin
		fPool[bin]->Add(fLocalList);
	      }
	    }
	  }//event loop 	 
	  
	  for (Int_t ipool=0; ipool<fHnSparseI.GetNbins(); ++ipool) 
	    fPool[ipool]->OptimizeStorage();

	  // Add the current file to the ESD/AOD chain
	  if(!path.IsNull()) fChain->AddFile(path); 	 
	  else if(!turl.IsNull()) fChain->AddFile(turl);
 	 
	}//detector tag cuts
      }//lhc tag cuts
    }//run tag cut 	 
  }//tag file loop 	 

  delete [] x;
  fCurrentBin = 0; // Initialize the current bin
}

// _________________________________________________________________________
void AliEventPoolSparse::SetRunCut(const char * cut){
  // Run selection cuts
  if (fRunCut) delete fRunCut;
  fRunCut = new TTreeFormula("fRun",cut,fTagChain);
}

// _________________________________________________________________________
void AliEventPoolSparse::SetLHCCut(const char * cut){
  // LHC selection cuts
  if (fLHCCut) delete fLHCCut;
  fLHCCut = new TTreeFormula("fLHC",cut,fTagChain);
}

// _________________________________________________________________________
void AliEventPoolSparse::SetDetCut(const char * cut){
  // Detector selection cuts
  if (fDetCut) delete fDetCut;
  fDetCut = new TTreeFormula("fDet",cut,fTagChain);
}

// _________________________________________________________________________
void AliEventPoolSparse::SetEventCut(const char * cut){
  // Event selection cuts
  if (fEvCut) delete fEvCut;
  fEvCut = new TTreeFormula("fEv",cut,fTagChain);
}

// _________________________________________________________________________
void AliEventPoolSparse::Set(Int_t n){
  // Set size of the array of pointers to n.
  // A new array is created, the old contents copied to the new array,
  // then the old array is deleted.
  // This function is taken from TArrayI

  if (n < 0) return;
  if (n != fN) {
    TEntryList **temp = fPool;
    if (n != 0) {
      fPool = new TEntryList*[n];
      if (n < fN) memcpy(fPool,temp, n*sizeof(TEntryList*));
      else {
	memcpy(fPool,temp,fN*sizeof(TEntryList*));
	memset(&fPool[fN],0x0,(n-fN)*sizeof(TEntryList*));
      }
    } else {
      fPool = 0x0;
    }
    if (fN) delete [] temp;
    fN = n;
  }
}
