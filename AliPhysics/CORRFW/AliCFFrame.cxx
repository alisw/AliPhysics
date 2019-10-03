/* $Id$ */
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
//--------------------------------------------------------------------//
//                                                                    //
// AliCFFrame Class                                                 //
// Class to accumulate data on an N-dimensional grid, to be used      //
// as input to get corrections for Reconstruction & Trigger efficiency// 
//                                                                    //
// -- Author : S.Arcelli                                              //
// Still to be done:                                                  //
// --Implement methods to merge cells                                 //
// --Interpolate among bins in a range                                // 
//--------------------------------------------------------------------//
//
//

#include "TSystem.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliCFFrame.h"

//____________________________________________________________________
ClassImp(AliCFFrame)

//____________________________________________________________________
AliCFFrame::AliCFFrame() : 
  TNamed()
{
  // default constructor
}

//____________________________________________________________________
AliCFFrame::AliCFFrame(const Char_t* name, const Char_t* title) : 
  TNamed(name,title)
{
  // named constructor
}

//____________________________________________________________________
void AliCFFrame::Save(const Char_t *outfile) const
{
  //
  // Save 'this' to a root file
  //

  const char *dirname = "./";
  TString filename = outfile;
  TFile *file=0x0;
  if((gSystem->FindFile(dirname,filename))!=NULL){
    file = new TFile( outfile,"UPDATE");
  }
  else{
    file = new TFile( outfile,"RECREATE");
  } 
  file->cd();
  //write the object to a file
  this->Write(GetName(),TObject::kSingleKey);
  file->Close();
  delete file;
}
