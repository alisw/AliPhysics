/**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
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


#include <TTree.h> 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>


#include "AliTOFMerger.h"
#include "AliTOF.h"
#include "AliTOFSDigitizer.h"
#include "AliTOFhit.h"
#include "AliTOFdigit.h"

#include "AliRun.h"
#include "AliPDG.h"

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

ClassImp(AliTOFMerger)

//___________________________________________
  AliTOFMerger::AliTOFMerger() 
{
// Default constructor    
    fEvNrSig = 0;
    fEvNrBgr = 0;
    fMerge =kDigitize;
    fFnBgr   = 0;
}

//------------------------------------------------------------------------
AliTOFMerger::~AliTOFMerger()
{
// Destructor
  if(fSDigits)  {
    fSDigits->Delete();
    delete fSDigits ;
    fSDigits = 0;
  }
}


//------------------------------------------------------------------------
void AliTOFMerger::Init()
{
// Initialisation
    if (fMerge) fBgrFile = InitBgr();
    
}



//------------------------------------------------------------------------
TFile* AliTOFMerger::InitBgr()
{
// Initialise background event
    TFile *file = new TFile(fFnBgr);
// add error checking later
    printf("\n AliTOFMerger has opened %s file with background event \n", fFnBgr);
    return file;
}

//------------------------------------------------------------------------
void AliTOFMerger::Digitise()
{

#ifdef DEBUG
  cout<<"ALiTOFMerger::>SDigits2Digits start...\n";
#endif
}
