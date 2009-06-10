/************************************************************************* 
* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. * 
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

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// parent class of all analyzers
// keeps the obj array of histograms filled by the daughter
// takes care of storing them on file at the end
//=============================================================================

#include <TROOT.h>
#include <TFile.h>
#include "AliUnicorAnal.h"

ClassImp(AliUnicorAnal)
  
TDatabasePDG AliUnicorAnal::fgPDG;

//=============================================================================
AliUnicorAnal::AliUnicorAnal(char *nam) : TNamed(nam,nam), fHistos() 
{
  // constructor

  fHistos.SetOwner(1);
  TDirectory *dir = gROOT->mkdir(GetName());
  dir->cd();

  printf("%s object named %s created\n",ClassName(),GetName());
}
//=============================================================================
void AliUnicorAnal::Save(const char *outfil, const char *mode) 
{
  // store histograms on file in a directory named after the object
  // mode should be "update" (default) or "new"

  printf("%s saving  histograms on %s (%s)\n",GetName(),outfil,mode);  
  TFile * f = TFile::Open(outfil, mode);
  TDirectory *dest = f->mkdir(GetName());
  dest->cd();
  fHistos.Write();
  gROOT->cd();
  f->Close();
}
//=============================================================================
