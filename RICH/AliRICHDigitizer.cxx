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

#include "AliRICHDigitizer.h"
#include "AliRICH.h"
#include "AliRICHDigit.h"
#include <AliRun.h>
#include <AliRunLoader.h>
#include "AliRunDigitizer.h"
#include <AliLoader.h>
#include <AliLog.h>


ClassImp(AliRICHDigitizer)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHDigitizer::Exec(Option_t*)
{
// This methode is responsible for merging sdigits to a list of digits
//Disintegration leeds to the fact that one hit affects several neighbouring pads, which means that the same pad might be affected by few hits.     
  AliDebug(1,Form("Start with %i input(s) for event %i",fManager->GetNinputs(),fManager->GetOutputEventNr()));
//First we read all sdigits from all inputs  
  AliRunLoader *pInRunLoader=0;//in and out Run loaders
  AliLoader    *pInRichLoader=0;//in and out RICH loaders  
  TClonesArray sdigs("AliRICHDigit");//tmp storage for sdigits sum up from all input files
  Int_t total=0;
  for(Int_t inFileN=0;inFileN<fManager->GetNinputs();inFileN++){//files loop
    pInRunLoader  = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inFileN));          //get run loader from current input 
    pInRichLoader = pInRunLoader->GetLoader("RICHLoader"); if(pInRichLoader==0) continue;       //no RICH in this input, check the next input
    if (!pInRunLoader->GetAliRun()) pInRunLoader->LoadgAlice();
    AliRICH* pInRich=(AliRICH*)pInRunLoader->GetAliRun()->GetDetector("RICH");                  //take RICH from current input
    pInRichLoader->LoadSDigits(); pInRichLoader->TreeS()->GetEntry(0);                          //take list of RICH sdigits from current input 
    AliDebug(1,Form("input %i has %i sdigits",inFileN,pInRich->SdiLst()->GetEntries()));
    for(Int_t i=0;i<pInRich->SdiLst()->GetEntries();i++){                                        //collect sdigits from current input
      AliRICHDigit *pSDig=(AliRICHDigit*)pInRich->SdiLst()->At(i);
      pSDig->AddTidOffset(fManager->GetMask(inFileN));                                          //apply TID shift since all inputs count tracks independently starting from 0
      new(sdigs[total++]) AliRICHDigit(*pSDig);       
    }
    pInRichLoader->UnloadSDigits();   pInRich->SdiReset(); //close current input and reset 
  }//files loop

  //PH  if(sdigs.GetEntries()==0) return;                                                              //no sdigits collected, nothing to convert  
  
  AliRunLoader *pOutRunLoader  = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());    //open output stream (only 1 possible)
  AliLoader    *pOutRichLoader = pOutRunLoader->GetLoader("RICHLoader");                         //take output RICH loader
  AliRICH      *pOutRich       = (AliRICH*)pOutRunLoader->GetAliRun()->GetDetector("RICH");      //take output RICH
  pOutRichLoader->MakeTree("D");   pOutRich->MakeBranch("D");                                    //create TreeD in output stream

  Sdi2Dig(&sdigs,pOutRich->DigLst());
  
  pOutRichLoader->TreeD()->Fill();              //fill the output tree with the list of digits
  pOutRichLoader->WriteDigits("OVERWRITE");     //serialize them to file
  
  sdigs.Clear();                      //remove all tmp sdigits
  pOutRichLoader->UnloadDigits();   pOutRich->DigReset();
}//Exec()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHDigitizer::Sdi2Dig(TClonesArray *pSdiLst,TObjArray *pDigLst)
{
// Converts list of sdigits to 7 lists of digits, one per each chamber
// Arguments: pSDigLst - list of all sdigits
//            pDigLst  - list of 7 lists of digits        
//   Returns: none  
  
  TClonesArray *pLst[7]; Int_t iCnt[7];
  
  for(Int_t i=0;i<7;i++){
    pLst[i]=(TClonesArray*)(*pDigLst)[i];
    iCnt[i]=pLst[i]->GetEntries();         //in principle those lists should be empty                                                                       
  }
  
  pSdiLst->Sort();  
                     
  Int_t iPad=-1,iCh=-1,iNdigPad=-1,aTids[3]={-1,-1,-1}; Float_t q=-1;
  for(Int_t i=0;i<pSdiLst->GetEntries();i++){                                                                  //sdigits loop (sorted)
    AliRICHDigit *pSdig=(AliRICHDigit*)pSdiLst->At(i);                                                         //take current sdigit
    if(pSdig->Pad()==iPad){                                                                                    //if the same pad 
      q+=pSdig->Q();                                                                                           //sum up charge
      iNdigPad++; if(iNdigPad<=3) aTids[iNdigPad-1]=pSdig->GetTrack(0);                                        //collect TID 
      continue;
    }
    if(i!=0 && AliRICHDigit::IsOverTh(q))  new((*pLst[iCh])[iCnt[iCh]++]) AliRICHDigit(iPad,(Int_t)q,aTids);   //do not create digit for the very first sdigit 
    iPad=pSdig->Pad(); iCh=AliRICHDigit::A2C(iPad);                                                            //new sdigit comes, reset collectors
    iNdigPad=1;
    aTids[0]=pSdig->GetTrack(0);aTids[1]=aTids[2]=-1; 
    q=pSdig->Q();    
  }//sdigits loop (sorted)
  
  if(AliRICHDigit::IsOverTh(q))  new((*pLst[iCh])[iCnt[iCh]++]) AliRICHDigit(iPad,(Int_t)q,aTids);             //add the last one, in case of empty sdigits list q=-1
                                                                                                               //so digit is not created    
}//Sdi2Dig()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
