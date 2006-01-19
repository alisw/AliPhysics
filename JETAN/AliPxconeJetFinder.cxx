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
// Pxcone Jet finder 
// manages the search for jets 
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include "AliPxconeJetFinder.h"
#include "AliPxconeJetHeader.h"
#include "AliJetReader.h"
#include "AliJet.h"

ClassImp(AliPxconeJetFinder)

////////////////////////////////////////////////////////////////////////

AliPxconeJetFinder::AliPxconeJetFinder()

{
  //
  // Constructor
  //
  fHeader = 0;
}

////////////////////////////////////////////////////////////////////////

AliPxconeJetFinder::~AliPxconeJetFinder()

{
  //
  // destructor
  //

  // reset and delete header
}

////////////////////////////////////////////////////////////////////////

#ifndef WIN32
# define pxcone pxcone_
# define type_of_call

#else
# define pxcone PXCONE
# define type_of_call _stdcall
#endif

extern "C" void type_of_call
pxcone_(int* mode, int* ntrak, int* dim, double* ptrak,
	double* coner, double* epslon, double* ovlim, 
	int * maxjet, int* njet, double* pjet, int* ipass, 
	int* ijmul, int* ierr);

void AliPxconeJetFinder::FindJets()

{
  // get number of entries
  // Int_t ne=fReader->GetNEvents();
  // loops over entries (check user request for number of events)
  // (this info should be stored in reader-header)

  // test with one event
  TClonesArray *lvArray = fReader->GetMomentumArray();
  Int_t nIn = lvArray->GetEntries();
  
  // local arrays for input
  const Int_t kNmaxVec = 30000;
  if (nIn > kNmaxVec) {
    cout << " AliPxconeJetFinder::FindJets: Too many input vectors."
	 << endl;
    cout << " Using only the first " << kNmaxVec << endl;
    nIn = kNmaxVec;
  }
  
  Double_t pIn[kNmaxVec][4];
  Double_t pJet[kNmaxVec][5];
  int ipass[kNmaxVec];
  int ijmul[kNmaxVec];
  Int_t ierr;

  // load input vectors
  for (Int_t i=0; i<nIn;i++){
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    pIn[i][0]= lv->Px();
    pIn[i][1]= lv->Py();
    pIn[i][2]= lv->Pz();
    pIn[i][3]= lv->E();
  }
  fJets->SetNinput(nIn);

  // run the algorithm. Get parameters from header
  Int_t dim=4;
  Int_t mode=fHeader->GetMode();
  Double_t radius=fHeader->GetRadius();
  Double_t minpt=fHeader->GetMinPt();
  Double_t ov=fHeader->GetOverlap();
  Int_t nj;
  pxcone(&mode,&nIn,&dim,&pIn[0][0],&radius,&minpt,&ov,&nIn,&nj,
	  &pJet[0][0],&ipass[0],&ijmul[0],&ierr);

  // download jets
  fJets->SetInJet(ipass);
  for(Int_t i=0; i<nj; i++) 
    fJets->AddJet(pJet[i][0],pJet[i][1],
		  pJet[i][2],pJet[i][3]);
  fJets->SetMultiplicities(ijmul);
}

////////////////////////////////////////////////////////////////////////

void AliPxconeJetFinder::WriteJHeaderToFile()
{
// Write Header to file

  fOut->cd();
  fHeader->Write();
}

////////////////////////////////////////////////////////////////////////

void AliPxconeJetFinder::Reset()
{
// Reset jet list
  fJets->ClearJets();
}

