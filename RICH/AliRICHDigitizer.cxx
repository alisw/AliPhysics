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


//Piotr.Skowronski@cern.ch :
//Corrections applied in order to compile (only) with new I/O and folder structure
//To be implemented correctly by responsible

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
AliRICHDigitizer::AliRICHDigitizer(AliRunDigitizer *pManager) 
                 :AliDigitizer(pManager)
{//main ctor which should be used
}//main ctor
//__________________________________________________________________________________________________
AliRICHDigitizer::~AliRICHDigitizer()
{//dtor
  if(fManager->GetDebug())Info("dtor","Start.");
}//dtor
//__________________________________________________________________________________________________
void AliRICHDigitizer::Exec(Option_t*)
{

  AliRunLoader *pInAL, *pOutAL;//in and out Run loaders
  AliLoader    *pInRL, *pOutRL;//in and out RICH loaders
 
  pOutAL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  pOutRL = pOutAL->GetLoader("RICHLoader");

  
     
  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
  
  if(!pOutRL->TreeD()) pOutRL->MakeTree("D");  pRICH->MakeBranch("D");
  
  
  for(Int_t inputFile=0;inputFile<fManager->GetNinputs();inputFile++){//files loop
    pInAL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    pInRL = pInAL->GetLoader("RICHLoader");
    
    
  }//files loop

      
      
  pOutRL->TreeD()->Fill();

    
  pRICH->ResetDigits(); /// ??? should it be here???
  
  pOutRL->WriteDigits("OVERWRITE");

  pOutRL->UnloadHits();
  pOutRL->UnloadDigits();
}//Exec()
