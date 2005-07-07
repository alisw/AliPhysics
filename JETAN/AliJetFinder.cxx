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
  
 
//---------------------------------------------------------------------
// Jet finder base class
// manages the search for jets 
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TFile.h>
#include "AliGenPythiaEventHeader.h"
#include "AliJetFinder.h"
#include "AliJet.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliJetControlPlots.h"
#include "AliLeading.h"
#include "AliHeader.h"


ClassImp(AliJetFinder)

////////////////////////////////////////////////////////////////////////

AliJetFinder::AliJetFinder()
{
  //
  // Constructor
  //
  fOut = 0;
  fJets = new AliJet();
  fGenJets = new AliJet();
  fLeading = new AliLeading();
  fReader = 0;
  fPlots = 0;
  SetPlotMode(kFALSE);
}


////////////////////////////////////////////////////////////////////////

AliJetFinder::~AliJetFinder()
{
  //
  // destructor
  //

  // here reset and delete jets
  fJets->ClearJets();
  delete fJets;
  fGenJets->ClearJets();
  delete fGenJets;
  // close file
  if (fOut) {
    fOut->Close();
    fOut->Delete();
  }
  delete fOut;
  // reset and delete control plots
  if (fPlots) delete fPlots;
  // delete fLeading;
}


////////////////////////////////////////////////////////////////////////

void AliJetFinder::SetOutputFile(const char *name)
{
  // opens output file 
  fOut = new TFile(name,"recreate");
}


////////////////////////////////////////////////////////////////////////

void AliJetFinder::PrintJets()
{
//
// Print jet information
  cout << " Jets found with jet algorithm:" << endl;
  fJets->PrintJets();
  cout << " Jets found by pythia:" << endl;
  fGenJets->PrintJets();
}


////////////////////////////////////////////////////////////////////////

void AliJetFinder::SetPlotMode(Bool_t b)
{
// Sets the plotting mode
  fPlotMode=b;
  if (b && !fPlots) fPlots = new AliJetControlPlots(); 
}

////////////////////////////////////////////////////////////////////////

void AliJetFinder::WriteJetsToFile(Int_t i)
{
// Writes the jets to file
  fOut->cd();
  char hname[30];
  sprintf(hname,"TreeJ%d",i);
  TTree* jetT = new TTree(hname,"AliJet");
  jetT->Branch("FoundJet",&fJets,1000);
  jetT->Branch("GenJet",&fGenJets,1000);
  jetT->Branch("LeadingPart",&fLeading,1000);
  jetT->Fill();
  jetT->Write(hname);
  delete jetT;
}

////////////////////////////////////////////////////////////////////////

void AliJetFinder::WriteRHeaderToFile()
{
  // write reader header
  fOut->cd();
  AliJetReaderHeader *rh = fReader->GetReaderHeader();
  rh->Write();
}


////////////////////////////////////////////////////////////////////////

void AliJetFinder::GetGenJets()
{
// Get the generated jet information from mc header
  fGenJets->SetNinput(1);
  AliHeader* alih = fReader->GetAliHeader(); 
  if (alih == 0) return;
  AliGenEventHeader * genh = alih->GenEventHeader();
  if (genh == 0) return;
  Int_t nj =((AliGenPythiaEventHeader*)genh)->NTriggerJets(); 
  Int_t* m = new Int_t[nj];
  Int_t* k = new Int_t[nj];
  for (Int_t i=0; i< nj; i++) {
    Float_t p[4];
    ((AliGenPythiaEventHeader*)genh)->TriggerJet(i,p);
    fGenJets->AddJet(p[0],p[1],p[2],p[3]);
    m[i]=1;
    k[i]=i;
  }
  fGenJets->SetMultiplicities(m);
  fGenJets->SetInJet(k);
}

////////////////////////////////////////////////////////////////////////

void AliJetFinder::Run()
{
  // do some initialization
  Init();

  // connect files
  fReader->OpenInputFiles();

  // write headers
  if (fOut) {
      fOut->cd();
      WriteRHeaderToFile();
      WriteJHeaderToFile();
  }

  // loop over events
  Int_t nFirst,nLast;
  nFirst = fReader->GetReaderHeader()->GetFirstEvent();
  nLast = fReader->GetReaderHeader()->GetLastEvent();
  // loop over events
  for (Int_t i=nFirst;i<nLast;i++) {
      fReader->FillMomentumArray(i);
      fLeading->FindLeading(fReader);
      GetGenJets();
      FindJets();
      if (fOut) WriteJetsToFile(i);
      if (fPlots) fPlots->FillHistos(fJets,fReader);
      fLeading->Reset();
      fGenJets->ClearJets();
      Reset();
  } 
  // write out
  if (fPlots) {
      fPlots->Normalize();
      fPlots->PlotHistos();
  }
  if (fOut) {
      fOut->cd();
      fPlots->Write();
      fOut->Close();
  }
}

