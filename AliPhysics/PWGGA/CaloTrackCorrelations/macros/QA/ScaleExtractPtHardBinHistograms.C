/// \file ScaleExtractPtHardBinHistograms.C
/// \ingroup CaloTrackCorrMacrosQAPtHard
/// \brief Scale and extract pT hard dependent histograms of analysis EMCal PWG-GA QA wagon
///
/// Macro to extract and scale QA MC productions  histograms 
/// done with pT hard bins (pythia jet-jet)./
/// Extract the histograms scaled and not scaled
///
/// To execute: root -q -b -l ScaleExtractPtHardBinHistograms.C'(1,50,kFALSE)'
///
/// The input files Scaled.root and NotScaled.root are obtained executing the script:
/// * DownloadExtractScaleMergePtHardAnalysisFiles.sh
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

//---------------------------------------------------------
// Set includes and declare methods for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TList.h"
#include "TDirectoryFile.h"

#endif // CINT

void ForBin
(
 Int_t bin = -1, 
 TString listName = "Pi0IM_GammaTrackCorr_EMCAL", 
 TString subListName = "default" 
);

///
/// Main method, loop on pT hard bins if requested
///
/// Input:
/// \param bin: bin integer to extract. If -1 loop on all bins from 1 to 20
/// \param listName:  string name of directory list 
/// \param subListName: string name of list 
//_______________________________________________________________________
void ScaleExtractPtHardBinHistograms
(Int_t bin = -1, 
 TString listName = "Pi0IM_GammaTrackCorr_EMCAL", 
 TString subListName = "default" )
{ 
  if(bin < 0)
  {
    for(Int_t i = 1; i < 21; i++)
    {
      ForBin(i,listName,subListName);
    }
  }
  else           
    ForBin(bin,listName,subListName);
}

///
/// Per pT hard bin scaling and extraction
///
/// Input:
/// \param bin: bin integer to extract. 
/// \param listName:  string name of directory list 
/// \param subListName: string name of list 
//_______________________________________________________________________
void ForBin ( Int_t bin, TString listName, TString subListName )
{ 
  printf("Scale pT hard bin %d\n",bin);
  
  TFile * f= TFile::Open(Form("%d/Merged.root",bin),"read");
  
  if(!f)
  {
    printf("No file available \n");
    return;
  }

  TDirectoryFile* dir = (TDirectoryFile*) f->Get(Form("%s",listName.Data()));
  
  if(!dir)
  {
    printf("No directory <%s> available \n",listName.Data());
    return;
  }
  
  TList* list = (TList*) dir->Get(Form("%s",subListName.Data()));
    
  if(!list)
  {
    printf("No list <%s> available \n", subListName.Data());
    return;
  }
       
  //
  // Extract into separate file without scaling
  //
  TFile * foutNo = new TFile(Form("%d/NotScaledMerged.root",bin),"recreate");
  
  foutNo->Print("");
  
  list->Write();
  
  foutNo->Close();
  delete foutNo;
  
  //
  // Calculate scaling factor
  //
  Double_t nEve   = ((TH1F*)list->FindObject("hNEvents"))->GetBinContent(1);
  Double_t xsec   = 1;
  Double_t trials = 1;
  Int_t    nfiles = 1;
  if((TH1F*)list->FindObject("hXsec"   ))
  {
    xsec   = ((TH1F*)list->FindObject("hXsec"   ))->GetBinContent(1);
    trials = ((TH1F*)list->FindObject("hTrials" ))->GetBinContent(1);
    nfiles = ((TH1F*)list->FindObject("hTrials" ))->GetEntries();
  }
  else printf("XSec histogram not found!!!!\n");
  printf(" nevents %e, xsec %2.2e, nTrials %e, files %d; ",nEve,xsec,trials,nfiles);

  trials/=nfiles;
  xsec  /=nfiles;
  
  Double_t scale = xsec/trials/nEve;

  printf(" average: xsec %2.2e, trials %2.2f, scale factor %2.3e\n",xsec,trials,scale);

  //
  // Scale each histogram, apply Sumw2 before.
  //
  
  printf("\t scaling ...\n ");
  
  TObject * h = 0 ; 

  for(Int_t i = 0; i < list->GetEntries(); i++) 
  { 
    h = list->At(i);
    if(h)
    {
      if ( !strncmp(h->ClassName(),"TH",2) ) 
      {
        //snprintf(name, buffersize, "%sScaled", h->GetName()) ; 
        
        TH1 * hout = dynamic_cast<TH1*> (h);//(h->Clone(name)) ; 
        
        if(hout)
        {
          hout->Sumw2();
          hout->Scale(scale) ;  
         // fOutputList->Add(hout) ;
        }// casting not null
      } 
      //else  fOutputList->Add(h) ; 
    }
  }
          
  printf("\t ... end scaling, write ...\n");
  
  TFile * fout = new TFile(Form("%d/ScaledMerged.root",bin),"recreate");
  
  fout->Print("");
  
  list->Write();
  
  fout->Close();
  delete fout;
  f->Close();
  delete f;
  delete h;

  printf("\t ... done.\n");
  
}
