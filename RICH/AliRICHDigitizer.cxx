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


#include <Riostream.h> 

#include <TTree.h> 
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TParticle.h>

#include <AliRunLoader.h>
#include <AliLoader.h>

#include "AliRICHDigitizer.h"
#include "AliRICH.h"
#include "AliRun.h"
#include "AliPDG.h"
#include "AliRunDigitizer.h"

ClassImp(AliRICHDigitizer)

//__________________________________________________________________________________________________
AliRICHDigitizer::AliRICHDigitizer() 
{//default constructor
}//default ctor
//__________________________________________________________________________________________________
AliRICHDigitizer::AliRICHDigitizer(AliRunDigitizer *pManager) 
                 :AliDigitizer(pManager)
{
//main ctor which should be used
  if(GetDebug())Info("main ctor","Start.");
  fRich=(AliRICH*)gAlice->GetDetector("RICH");
}//main ctor
//__________________________________________________________________________________________________
AliRICHDigitizer::~AliRICHDigitizer()
{
//dtor
  if(GetDebug())Info("dtor","Start.");
}//dtor
//__________________________________________________________________________________________________
void AliRICHDigitizer::Exec(Option_t*)
{
  if(GetDebug())Info("Exec","Start with %i input(s) and %i",fManager->GetNinputs(),fManager->GetOutputEventNr());
  
  AliRunLoader *pInAL, *pOutAL;//in and out Run loaders
  AliLoader    *pInRL, *pOutRL;//in and out RICH loaders
 
  pOutAL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  pOutRL = pOutAL->GetLoader("RICHLoader");

  
     
  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
  
  if(!pOutRL->TreeD()) pOutRL->MakeTree("D");  pRICH->MakeBranch("D");
  
  
  for(Int_t inputFile=0;inputFile<fManager->GetNinputs();inputFile++){//files loop
    pInAL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    pInRL = pInAL->GetLoader("RICHLoader");
    pInRL->LoadSDigits();
    pInRL->TreeD()->GetEntries();
  }//files loop

      
      
  pOutRL->TreeD()->Fill();

    
  pRICH->ResetDigits(); /// ??? should it be here???
  
  pOutRL->WriteDigits("OVERWRITE");

  pOutRL->UnloadHits();
  pOutRL->UnloadDigits();
  if(GetDebug())Info("Exec","Stop");
}//Exec()
//__________________________________________________________________________________________________
