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
//  Algorithm class for the identification of particles detected in PHOS        
//  base  class  of identified particle  
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                
                         
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko


// --- ROOT system ---
#include "TGeometry.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
 
// --- Standard library ---
#include <stdlib.h>


// --- AliRoot header files ---
#include "AliRun.h" 
#include "AliPHOSPID.h"
#include "AliHeader.h" 

ClassImp(AliPHOSPID)

//____________________________________________________________________________
  AliPHOSPID::AliPHOSPID():TTask("","")
{
  // ctor
  fSplitFile= 0 ; 

}


//____________________________________________________________________________
AliPHOSPID::AliPHOSPID(const char* headerFile, const char * name, const Bool_t toSplit):TTask(name, headerFile)
{
  // ctor

  fToSplit = toSplit ;
  fSplitFile= 0 ; 
}

//____________________________________________________________________________
AliPHOSPID::~AliPHOSPID()
{
  // dtor
        
  fSplitFile = 0 ;
}

// //____________________________________________________________________________
// void AliPHOSPID::SetSplitFile(const TString splitFileName) const
// {
//   // Diverts the Digits in a file separate from the hits file
  

//   TDirectory * cwd = gDirectory ;
//   TFile * splitFile = gAlice->InitTreeFile("R",splitFileName.Data());
//   splitFile->cd() ; 
//   gAlice->Write(0, TObject::kOverwrite);
  
//   TTree *treeE  = gAlice->TreeE();
//   if (!treeE) {
//     cerr << "ERROR: AliPHOSPID::SetSplitFile -> No TreeE found "<<endl;
//     abort() ;
//   }      
  
//   // copy TreeE
//   AliHeader *header = new AliHeader();
//   treeE->SetBranchAddress("Header", &header);
//   treeE->SetBranchStatus("*",1);
//   TTree *treeENew =  treeE->CloneTree();
//   treeENew->Write(0, TObject::kOverwrite);
  
//   // copy AliceGeom
//   TGeometry *AliceGeom = static_cast<TGeometry*>(cwd->Get("AliceGeom"));
//   if (!AliceGeom) {
//     cerr << "ERROR: AliPHOSPID::SetSplitFile -> AliceGeom was not found in the input file "<<endl;
//     abort() ;
//     }
//   AliceGeom->Write(0, TObject::kOverwrite) ;
  
//   gAlice->MakeTree("R",splitFile);
//   cwd->cd() ; 
//   cout << "INFO: AliPHOSPID::SetSPlitMode -> RecParticles will be stored in " << splitFileName.Data() << endl ;   
// }
