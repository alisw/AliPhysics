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

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*--
//*-- Author: Yves Schutz  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TGeometry.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"


// --- Standard library ---
#include <iostream.h>
#include <stdlib.h>   

// --- AliRoot header files ---
#include "AliRun.h" 
#include "AliPHOSClusterizer.h"
#include "AliHeader.h" 
#include "AliPHOSGetter.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSDigitizer.h"

ClassImp(AliPHOSClusterizer)

//____________________________________________________________________________
  AliPHOSClusterizer::AliPHOSClusterizer():TTask("","")
{
  // ctor
  fSplitFile= 0 ; 
  fHitsFileName = "" ; 
  fSDigitsFileName = "" ; 
  fDigitsFileName = "" ; 

}

//____________________________________________________________________________
AliPHOSClusterizer::AliPHOSClusterizer(const char* headerFile, const char* name):
TTask(name, headerFile)
{
  // ctor
  fSplitFile= 0 ; 
  fDigitsFileName = headerFile ; 
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(headerFile, name) ; 
  gime->Event(0,"D") ; 
  fSDigitsFileName = gime->Digitizer()->GetTitle() ; 
  gime = AliPHOSGetter::GetInstance(fSDigitsFileName, name) ; 
  gime->Event(0,"S") ; 
  fHitsFileName = gime->SDigitizer()->GetTitle() ; 
}

//____________________________________________________________________________
AliPHOSClusterizer::~AliPHOSClusterizer()
{
  // dtor
         
  fSplitFile = 0 ; ;
}

//____________________________________________________________________________
void AliPHOSClusterizer::SetSplitFile(const TString splitFileName) 
{
  // Diverts the RecPoints in a file separate from the Digits file
  

  TDirectory * cwd = gDirectory ;
  fSplitFile = gAlice->InitTreeFile("R",splitFileName.Data());
  fSplitFile->cd() ; 
  gAlice->Write(0, TObject::kOverwrite);

  TTree *treeE  = gAlice->TreeE();
  if (!treeE) {
    cerr << "ERROR: AliPHOSClusterizer::SetSplitFile -> No TreeE found "<<endl;
    abort() ;
  }      
  
  // copy TreeE
  AliHeader *header = new AliHeader();
  treeE->SetBranchAddress("Header", &header);
  treeE->SetBranchStatus("*",1);
  TTree *treeENew =  treeE->CloneTree();
  treeENew->Write(0, TObject::kOverwrite);
    
  // copy AliceGeom
  TGeometry *AliceGeom = static_cast<TGeometry*>(cwd->Get("AliceGeom"));
  if (!AliceGeom) {
    cerr << "ERROR: AliPHOSClusterizer::SetSplitFile -> AliceGeom was not found in the input file "<<endl;
    abort() ;
  }
  AliceGeom->Write(0, TObject::kOverwrite);
  
  gAlice->MakeTree("R", fSplitFile);
  cwd->cd() ; 
  cout << "INFO: AliPHOSClusterizer::SetSPlitMode -> RecPoints will be stored in " << splitFileName.Data() << endl ;   
}
