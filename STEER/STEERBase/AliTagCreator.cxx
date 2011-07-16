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
//           AliTagCreator class
//   This is the class to deal with the tag creation (post process)
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

//ROOT
#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>

//ROOT-AliEn
#include <TGrid.h>
#include <TGridResult.h>

//AliRoot
#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliLog.h"

#include "AliTagCreator.h"


ClassImp(AliTagCreator)


//______________________________________________________________________________
  AliTagCreator::AliTagCreator() :
    TObject(),
    fSE("ALICE::CERN::se"),
    fgridpath(""),
    fStorage(0)
{
  //==============Default constructor for a AliTagCreator==================
}

//______________________________________________________________________________
AliTagCreator::~AliTagCreator() {
//================Default destructor for a AliTagCreator=======================
}

//______________________________________________________________________________
void AliTagCreator::SetStorage(Int_t storage) {
  // Sets correctly the storage: 0 for local, 1 for GRID
  fStorage = storage;
  if(fStorage == 0)
    AliInfo(Form("Tags will be stored locally...."));
  if(fStorage == 1)
    AliInfo(Form("Tags will be stored in the grid...."));
  if((fStorage != 0)&&(fStorage != 1))
    {
      AliInfo(Form("Storage was not properly set!!!"));
      abort();
    }  
}

//__________________________________________________________________________
Bool_t AliTagCreator::MergeTags(const char *type) {
  //Merges the tags and stores the merged tag file 
  //locally if fStorage=0 or in the grid if fStorage=1
  AliInfo(Form("Merging tags....."));
  TChain *fgChain = new TChain("T");
  TString tagPattern = type; tagPattern += ".tag.root"; 

  if(fStorage == 0) {
    // Open the working directory
    void * dirp = gSystem->OpenDirectory(gSystem->pwd());
    const char * name = 0x0;
    // Add all files matching *pattern* to the chain
    while((name = gSystem->GetDirEntry(dirp))) {
      if (strstr(name,tagPattern)) fgChain->Add(name);  
    }//directory loop
    AliInfo(Form("Chained tag files: %lld",fgChain->GetEntries()));
  }//local mode

  else if(fStorage == 1) {
    TString alienLocation = gGrid->Pwd();
    alienLocation += fgridpath.Data();
    alienLocation += "/";
    TString queryPattern = "*."; queryPattern += tagPattern;
    TGridResult *tagresult = gGrid->Query(alienLocation,queryPattern.Data(),"","");
    Int_t nEntries = tagresult->GetEntries();
    for(Int_t i = 0; i < nEntries; i++) {
      TString alienUrl = tagresult->GetKey(i,"turl");
      fgChain->Add(alienUrl);
    }//grid result loop      
    AliInfo(Form("Chained tag files: %lld",fgChain->GetEntries()));
  }//grid mode
 
  AliRunTag *tag = new AliRunTag;
  fgChain->SetBranchAddress("AliTAG",&tag);
  fgChain->GetEntry(0);
  TString localFileName = "Run"; localFileName += tag->GetRunId(); 
  localFileName += ".Merged."; localFileName += tagPattern.Data();
     
  TString filename;
  
  if(fStorage == 0) {
    filename = localFileName.Data();      
    AliInfo(Form("Writing merged tags to local file: %s",filename.Data()));
  } 
  else if(fStorage == 1) {
    TString alienFileName = "/alien";
    alienFileName += gGrid->Pwd();
    alienFileName += fgridpath.Data();
    alienFileName += "/";
    alienFileName +=  localFileName;
    alienFileName += "?se=";
    alienFileName += fSE.Data();
    filename = alienFileName.Data();
    AliInfo(Form("Writing merged tags to grid file: %s",filename.Data()));     
  }

  fgChain->Merge(filename);
  gSystem->Exec("rm Run*.Event*");

  return kTRUE;
}


//__________________________________________________________________________
Bool_t AliTagCreator::MergeTags(const char *type, const char *inflist) {
  //Merges the tags and stores the merged tag file 
  //locally if fStorage=0 or in the grid if fStorage=1
  AliInfo(Form("Merging tags....."));
  TChain *fgChain = new TChain("T");
  TString tagPattern = type; tagPattern += ".tag.root"; 

//   if(fStorage == 0) {
//     // Open the working directory
//     void * dirp = gSystem->OpenDirectory(gSystem->pwd());
//     const char * name = 0x0;
//     // Add all files matching *pattern* to the chain
//     while((name = gSystem->GetDirEntry(dirp))) {
//       if (strstr(name,tagPattern)) fgChain->Add(name);  
//     }//directory loop
//     AliInfo(Form("Chained tag files: %d",fgChain->GetEntries()));
//   }//local mode

//   else if(fStorage == 1) {
//     TString alienLocation = gGrid->Pwd();
//     alienLocation += fgridpath.Data();
//     alienLocation += "/";
//     TString queryPattern = "*."; queryPattern += tagPattern;
//     TGridResult *tagresult = gGrid->Query(alienLocation,queryPattern.Data(),"","");
//     Int_t nEntries = tagresult->GetEntries();
//     for(Int_t i = 0; i < nEntries; i++) {
//       TString alienUrl = tagresult->GetKey(i,"turl");
//       fgChain->Add(alienUrl);
//     }//grid result loop      
//     AliInfo(Form("Chained tag files: %d",fgChain->GetEntries()));
//   }//grid mode

  if (fStorage == 0) {
    ifstream *istr = new ifstream(inflist);
    char fname[200];
    while (!(istr->eof())) {
      (*istr) >> fname;
      if (strstr(fname, tagPattern)) fgChain->Add(fname);
    }
    AliInfo(Form("Chained tag files: %lld",fgChain->GetEntries()));
  }
 
  AliRunTag *tag = new AliRunTag;
  fgChain->SetBranchAddress("AliTAG",&tag);
  fgChain->GetEntry(0);
  TString localFileName = "Run"; localFileName += tag->GetRunId(); 
  localFileName += ".Merged."; localFileName += tagPattern.Data();
     
  TString filename;
  
  if(fStorage == 0) {
    filename = localFileName.Data();      
    AliInfo(Form("Writing merged tags to local file: %s",filename.Data()));
  } 
  else if(fStorage == 1) {
    TString alienFileName = "/alien";
    alienFileName += gGrid->Pwd();
    alienFileName += fgridpath.Data();
    alienFileName += "/";
    alienFileName +=  localFileName;
    alienFileName += "?se=";
    alienFileName += fSE.Data();
    filename = alienFileName.Data();
    AliInfo(Form("Writing merged tags to grid file: %s",filename.Data()));     
  }

  fgChain->Merge(filename);
  gSystem->Exec("rm Run*.Event*");

  return kTRUE;
}

//__________________________________________________________________________
Bool_t AliTagCreator::MergeTags(const char *type, TGridResult *result) {
  //Merges the tags that are listed in the TGridResult 
  AliInfo(Form("Merging tags....."));
  TChain *fgChain = new TChain("T");
  TString tagPattern = "."; tagPattern += type; tagPattern += ".tag.root";

  Int_t nEntries = result->GetEntries();

  TString alienUrl;
  for(Int_t i = 0; i < nEntries; i++) {
    alienUrl = result->GetKey(i,"turl");
    fgChain->Add(alienUrl);  
  }
  AliInfo(Form("Chained tag files: %lld",fgChain->GetEntries()));
  AliRunTag *tag = new AliRunTag;
  fgChain->SetBranchAddress("AliTAG",&tag);
  fgChain->GetEntry(0);
    
  TString localFileName = "Run"; localFileName += tag->GetRunId(); 
  localFileName += ".Merged"; localFileName += tagPattern.Data();
     
  TString filename;
  
  if(fStorage == 0) {
    filename = localFileName.Data();      
    AliInfo(Form("Writing merged tags to local file: %s",filename.Data()));
  } 
  else if(fStorage == 1) {
    TString alienFileName = "/alien";
    alienFileName += gGrid->Pwd();
    alienFileName += fgridpath.Data();
    alienFileName += "/";
    alienFileName +=  localFileName;
    alienFileName += "?se=";
    alienFileName += fSE.Data();
    filename = alienFileName.Data();
    AliInfo(Form("Writing merged tags to grid file: %s",filename.Data()));     
  }
  
  fgChain->Merge(filename);

  return kTRUE;
}

Bool_t AliTagCreator::MergeTagsForRun(const char* type) {
  //Merges the tags and stores the merged tag file 
  //locally if fStorage=0 or in the grid if fStorage=1
  AliInfo(Form("Merging tags....."));
  TChain *fgChain = new TChain("T");
  TString tagPattern = type; tagPattern += ".tag.root"; 

  if(fStorage == 0) {
    // Open the working directory
    void * dirp = gSystem->OpenDirectory(gSystem->pwd());
    const char * name = 0x0;
    // Add all files matching *pattern* to the chain
    while((name = gSystem->GetDirEntry(dirp))) {
      if (strstr(name,tagPattern)) fgChain->Add(name);  
    }//directory loop
    AliInfo(Form("Chained tag files: %lld",fgChain->GetEntries()));
  }//local mode

  else if(fStorage == 1) {
    TString alienLocation = gGrid->Pwd();
    alienLocation += fgridpath.Data();
    alienLocation += "/";
    TString queryPattern = "*."; queryPattern += tagPattern;
    TGridResult *tagresult = gGrid->Query(alienLocation,queryPattern.Data(),"","");
    Int_t nEntries = tagresult->GetEntries();
    for(Int_t i = 0; i < nEntries; i++) {
      TString alienUrl = tagresult->GetKey(i,"turl");
      fgChain->Add(alienUrl);
    }//grid result loop      
    AliInfo(Form("Chained tag files: %lld",fgChain->GetEntries()));
  }//grid mode
 
  AliRunTag *tag = new AliRunTag;
  fgChain->SetBranchAddress("AliTAG",&tag);
  fgChain->GetEntry(0);
  TString localFileName = "Run"; localFileName += tag->GetRunId(); 
  localFileName += ".Merged."; localFileName += tagPattern.Data();
     
  TString filename;
  
  if(fStorage == 0) {
    filename = localFileName.Data();      
    AliInfo(Form("Writing merged tags to local file: %s",filename.Data()));
  } 
  else if(fStorage == 1) {
    TString alienFileName = "/alien";
    alienFileName += gGrid->Pwd();
    alienFileName += fgridpath.Data();
    alienFileName += "/";
    alienFileName +=  localFileName;
    alienFileName += "?se=";
    alienFileName += fSE.Data();
    filename = alienFileName.Data();
    AliInfo(Form("Writing merged tags to grid file: %s",filename.Data()));     
  }

  //  fgChain->Merge(filename);
  MergeToSingleRunTag(fgChain, filename);

  gSystem->Exec("rm Run*.Event*");

  return kTRUE;
}

Bool_t AliTagCreator::MergeTagsForRun(const char* type, TGridResult *result) {
  //Merges the tags that are listed in the TGridResult 
  AliInfo(Form("Merging tags....."));
  TChain *fgChain = new TChain("T");
  TString tagPattern = "."; tagPattern += type; tagPattern += ".tag.root";

  Int_t nEntries = result->GetEntries();

  TString alienUrl;
  for(Int_t i = 0; i < nEntries; i++) {
    alienUrl = result->GetKey(i,"turl");
    fgChain->Add(alienUrl);  
  }
  AliInfo(Form("Chained tag files: %lld",fgChain->GetEntries()));
  AliRunTag *tag = new AliRunTag;
  fgChain->SetBranchAddress("AliTAG",&tag);
  fgChain->GetEntry(0);
    
  TString localFileName = "Run"; localFileName += tag->GetRunId(); 
  localFileName += ".Merged"; localFileName += tagPattern.Data();
     
  TString filename;
  
  if(fStorage == 0) {
    filename = localFileName.Data();      
    AliInfo(Form("Writing merged tags to local file: %s",filename.Data()));
  } 
  else if(fStorage == 1) {
    TString alienFileName = "/alien";
    alienFileName += gGrid->Pwd();
    alienFileName += fgridpath.Data();
    alienFileName += "/";
    alienFileName +=  localFileName;
    alienFileName += "?se=";
    alienFileName += fSE.Data();
    filename = alienFileName.Data();
    AliInfo(Form("Writing merged tags to grid file: %s",filename.Data()));     
  }
  
  //  fgChain->Merge(filename);
  MergeToSingleRunTag(fgChain, filename);

  return kTRUE;
}

Bool_t AliTagCreator::MergeTagsForRun(const char* type, const char *inflist) {
  //Merges the tags and stores the merged tag file 
  //locally if fStorage=0 or in the grid if fStorage=1
  AliInfo(Form("Merging tags....."));
  TChain *fgChain = new TChain("T");
  TString tagPattern = type; tagPattern += ".tag.root"; 

//   if(fStorage == 0) {
//     // Open the working directory
//     void * dirp = gSystem->OpenDirectory(gSystem->pwd());
//     const char * name = 0x0;
//     // Add all files matching *pattern* to the chain
//     while((name = gSystem->GetDirEntry(dirp))) {
//       if (strstr(name,tagPattern)) fgChain->Add(name);  
//     }//directory loop
//     AliInfo(Form("Chained tag files: %d",fgChain->GetEntries()));
//   }//local mode

//   else if(fStorage == 1) {
//     TString alienLocation = gGrid->Pwd();
//     alienLocation += fgridpath.Data();
//     alienLocation += "/";
//     TString queryPattern = "*."; queryPattern += tagPattern;
//     TGridResult *tagresult = gGrid->Query(alienLocation,queryPattern.Data(),"","");
//     Int_t nEntries = tagresult->GetEntries();
//     for(Int_t i = 0; i < nEntries; i++) {
//       TString alienUrl = tagresult->GetKey(i,"turl");
//       fgChain->Add(alienUrl);
//     }//grid result loop      
//     AliInfo(Form("Chained tag files: %d",fgChain->GetEntries()));
//   }//grid mode

  if (fStorage == 0) {
    ifstream *istr = new ifstream(inflist);
    char fname[200];
    while (!(istr->eof())) {
      (*istr) >> fname;
      if (strstr(fname, tagPattern)) fgChain->Add(fname);
    }
    AliInfo(Form("Chained tag files: %lld",fgChain->GetEntries()));
  }
 
  AliRunTag *tag = new AliRunTag;
  fgChain->SetBranchAddress("AliTAG",&tag);
  fgChain->GetEntry(0);
  TString localFileName = "Run"; localFileName += tag->GetRunId(); 
  localFileName += ".Merged."; localFileName += tagPattern.Data();
     
  TString filename;
  
  if(fStorage == 0) {
    filename = localFileName.Data();      
    AliInfo(Form("Writing merged tags to local file: %s",filename.Data()));
  } 
  else if(fStorage == 1) {
    TString alienFileName = "/alien";
    alienFileName += gGrid->Pwd();
    alienFileName += fgridpath.Data();
    alienFileName += "/";
    alienFileName +=  localFileName;
    alienFileName += "?se=";
    alienFileName += fSE.Data();
    filename = alienFileName.Data();
    AliInfo(Form("Writing merged tags to grid file: %s",filename.Data()));     
  }

  //  fgChain->Merge(filename);
  MergeToSingleRunTag(fgChain, filename);

  gSystem->Exec("rm Run*.Event*");

  return kTRUE;
}


Bool_t AliTagCreator::MergeToSingleRunTag(TChain *chain, const char *filename)
{
  // Merge all tags for a given run into a single RunTag
  TFile* ftag = TFile::Open(filename, "recreate");

  AliRunTag *tag = new AliRunTag;
  TTree * ttag = new TTree("T","A Tree with event tags");
  TBranch * btag = ttag->Branch("AliTAG", &tag, 1000000);
  btag->SetCompressionLevel(9);
  ttag->AutoSave("10000");

  AliRunTag *rtag = new AliRunTag();
  chain->SetBranchAddress("AliTAG", &rtag);

  AliFileTag *evt;

  if (chain->GetEntries()) {
    chain->GetEntry(0);
    tag->CopyStandardContent(rtag);
    tag->Clear();

    int runno = rtag->GetRunId();

    for (int iter=0; iter<chain->GetEntries(); iter++) {
      chain->GetEntry(iter);
      if (runno != rtag->GetRunId()) {
	AliInfo(Form("Run tag ID %i is different from the Run ID for the merged run: %i\n", rtag->GetRunId(), runno));
	continue;
      }

      for (int iev=0; iev<rtag->GetNFiles(); iev++) {
	evt = (AliFileTag *) rtag->GetFileTag(iev);
	tag->AddFileTag(new AliFileTag(*evt));
      }
    }
  }
  else {
    AliInfo("Found no tag files to merge.");
    return kFALSE;
  }

  ttag->Fill();

  ftag->cd();
  tag->Clear();
  ttag->Write();
  ftag->Close();
  
  return kTRUE;
}


