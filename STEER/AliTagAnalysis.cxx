/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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

//-----------------------------------------------------------------
//           AliTagAnalysis class
//   This is the class to deal with the tag analysis
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

//ROOT
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TEventList.h>

//ROOT-AliEn
#include <TGridResult.h>

#include "AliLog.h"

#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliTagAnalysis.h"
#include "AliEventTagCuts.h"

class TTree;

ClassImp(AliTagAnalysis)

TChain *AliTagAnalysis::fgChain = 0;

//______________________________________________________________________________
AliTagAnalysis::AliTagAnalysis(): TObject()//local mode
{
  //==============Default constructor for a AliTagAnalysis==================
  ftagresult = 0;
}

//______________________________________________________________________________
AliTagAnalysis::~AliTagAnalysis()
{
//================Default destructor for a AliTagAnalysis=======================
}

//______________________________________________________________________________
void AliTagAnalysis::ChainLocalTags(const char *dirname) //local version
{
  //Searches the entries of the provided direcory
  //Chains the tags that are stored locally
  fTagDirName = dirname;
  TString fTagFilename;
  
  TChain *fgChain = new TChain("T");
  fChain = fgChain;
  
  const char * tagPattern = "tag";
  // Open the working directory
  void * dirp = gSystem->OpenDirectory(fTagDirName);
  const char * name = 0x0;
  // Add all files matching *pattern* to the chain
  while((name = gSystem->GetDirEntry(dirp)))
    {
      if (strstr(name,tagPattern))
	{ 
	  fTagFilename = fTagDirName;
	  fTagFilename += "/";
	  fTagFilename += name;
	  	
	  TFile * fTag = TFile::Open(fTagFilename);
	  if((!fTag) || (!fTag->IsOpen()))
	    {
	      AliError(Form("Tag file not opened!!!"));
	      continue;
	    } 
	  fChain->Add(fTagFilename);  
	  fTag->Close();
	  delete fTag;
	}//pattern check
    }//directory loop
  AliInfo(Form("Chained tag files: %d ",fChain->GetEntries()));
}


//______________________________________________________________________________
void AliTagAnalysis::ChainGridTags(TGridResult *res)
{
  //Loops overs the entries of the TGridResult
  //Chains the tags that are stored in the GRID
  ftagresult = res;
  Int_t nEntries = ftagresult->GetEntries();
 
  TChain *fgChain = new TChain("T");
  fChain = fgChain;

  TString gridname = "alien://";
  TString alienUrl;
 
  for(Int_t i = 0; i < nEntries; i++)
    {
      alienUrl = ftagresult->GetKey(i,"turl");
      TFile *f = TFile::Open(alienUrl,"READ");
      fChain->Add(alienUrl);
      //f->Close();
      delete f;	 
    }//grid result loop  
}


//______________________________________________________________________________
TList *AliTagAnalysis::QueryTags(AliEventTagCuts *EvTagCuts)
{
  //Queries the tag chain using the defined 
  //event tag cuts from the AliEventTagCuts object
  AliInfo(Form("Querying the tags........"));
  
  //file info list
  TList *list = new TList();
  
  Int_t iAccepted = 0, evCounter = 0;
  
  //Defining tag objects
  AliRunTag *tag = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  Long64_t size = -1;
  const char* md5 = 0;
  const char* guid = 0;
  const char* turl = 0;

  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++)
    {
      TEventList *fEventList = new TEventList();
      evCounter = 0;
      fChain->GetEntry(iTagFiles);
      Int_t iEvents = tag->GetNEvents();
      const TClonesArray *tagList = tag->GetEventTags();
      for(Int_t i = 0; i < iEvents; i++)
	{
	  evTag = (AliEventTag *) tagList->At(i);
	  size = evTag->GetSize();
	  md5 = evTag->GetMD5();
	  guid = evTag->GetGUID(); 
	  turl = evTag->GetTURL(); 
	  if(EvTagCuts->IsAccepted(evTag))
	    {
	      fEventList->Enter(i);
	      evCounter++;
	    
	      iAccepted++;
	    }
	}//event loop

      //adding a TFileInfo object to the list
      if(evCounter != 0)
	list->Add(new TFileInfo(turl,size,guid,md5,-1,-1,-1,fEventList));
      fEventList->Clear("");
      delete fEventList;
    }//tag file loop
  AliInfo(Form("Accepted events: %d",iAccepted));
   
  return list;
}


