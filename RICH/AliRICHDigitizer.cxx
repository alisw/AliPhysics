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
  Info("main ctor","Start.");
  fRich=(AliRICH*)gAlice->GetDetector("RICH");
  Rich()->Param()->GenSigmaThMap();
}//main ctor
//__________________________________________________________________________________________________
AliRICHDigitizer::~AliRICHDigitizer()
{
//dtor
  Info("dtor","Start.");
}//dtor
//__________________________________________________________________________________________________
void AliRICHDigitizer::Exec(Option_t*)
{
  Info("Exec","\n\n\n\n");
  Info("Exec","Start with %i input(s) for event %i",fManager->GetNinputs(),fManager->GetOutputEventNr());
  
  AliRunLoader *pInAL=0, *pOutAL;//in and out Run loaders
  AliLoader    *pInRL=0, *pOutRL;//in and out RICH loaders
 
  pOutAL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  pOutRL = pOutAL->GetLoader("RICHLoader");
  pOutRL->MakeTree("D");   Rich()->MakeBranch("D"); //create TreeD with RICH branches in output stream

  for(Int_t inFileN=0;inFileN<fManager->GetNinputs();inFileN++){//files loop
    pInAL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inFileN)); pInRL = pInAL->GetLoader("RICHLoader");
    pInAL->GetEvent(fManager->GetOutputEventNr());    pInRL->LoadSDigits();    pInRL->TreeS()->GetEntry(0);
    Info("Exec","Run Loader %p RICH Loader %p RICH %p",pInAL,pInRL,Rich());
    Rich()->SDigits()->Print();
    pInRL->UnloadSDigits();
  }//files loop
  
//     Rich()->SDigits()->Sort();                     //sort them according to Id() methode
//       
//     Int_t combiPid=0,chamber=0,x=0,y=0,tid[3],id=0; Double_t q=0;
//     Int_t iNdigitsPerPad=0;//how many sdigits for a given pad
//     for(Int_t i=0;i<Rich()->SDigits()->GetEntries();i++){//sdigits loop (sorted)
//       AliRICHdigit *pSdig=(AliRICHdigit*)Rich()->SDigits()->At(i);
//       if(pSdig->Id()==id){//still the same pad
//         iNdigitsPerPad++;
//         q+=pSdig->Q();
//         combiPid+=pSdig->CombiPid();
//         if(iNdigitsPerPad<=3)
//           tid[iNdigitsPerPad-1]=pSdig->Tid(0);
//         else
//           Warning("SDigits2Digits","More then 3 sdigits for the given pad");
//       }else{//new pad, add the pevious one
//         if(id!=kBad&&Rich()->Param()->IsOverTh(chamber,x,y,q)) {
//            Rich()->AddDigit(chamber,x,y,(Int_t)q,combiPid,tid); //add newly created digit to the list of digits
//          }
//         combiPid=pSdig->CombiPid();chamber=pSdig->C();id=pSdig->Id();
//         x=pSdig->X();y=pSdig->Y();
//         q=pSdig->Q();
//         tid[0]=pSdig->Tid(0);
//         iNdigitsPerPad=1;tid[1]=tid[2]=kBad;
//       }
//     }//sdigits loop (sorted)
//   
//     if(Rich()->SDigits()->GetEntries()&&Rich()->Param()->IsOverTh(chamber,x,y,q))
//       Rich()->AddDigit(chamber,x,y,(Int_t)q,combiPid,tid);//add the last digit
//         
  pOutRL->TreeD()->Fill();              //fill the tree with the list of digits
  pOutRL->WriteDigits("OVERWRITE");     //serialize them to file
            
  Info("Exec","Stop\n\n\n\n");
}//Exec()
//__________________________________________________________________________________________________
