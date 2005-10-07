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
#include <AliRun.h>
#include <AliRunLoader.h>
#include "AliRunDigitizer.h"
#include <AliLoader.h>
#include <AliLog.h>


ClassImp(AliRICHDigitizer)

//__________________________________________________________________________________________________
void AliRICHDigitizer::Exec(Option_t*)
{
//This methode is responsible for merging sdigits to a list of digits
//Disintegration leeds to the fact that one hit affected several neighbouring pads, which means that the same pad might be
//affected by few hits.     
  AliDebug(1,Form("Start with %i input(s) for event %i",fManager->GetNinputs(),fManager->GetOutputEventNr()));
//First we read all sdigits from all inputs  
  AliRunLoader *pInRunLoader=0;//in and out Run loaders
  AliLoader    *pInRichLoader=0;//in and out RICH loaders  
  TClonesArray tmpCA("AliRICHDigit");//tmp storage for sdigits sum up from all input files
  Int_t total=0;
  for(Int_t inFileN=0;inFileN<fManager->GetNinputs();inFileN++){//files loop
    pInRunLoader  = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inFileN));          //get run loader from current input 
    pInRichLoader = pInRunLoader->GetLoader("RICHLoader"); if(pInRichLoader==0) continue;       //no RICH in this input, check the next input
    if (!pInRunLoader->GetAliRun()) pInRunLoader->LoadgAlice();
    AliRICH* pInRich=(AliRICH*)pInRunLoader->GetAliRun()->GetDetector("RICH");                  //take RICH from current input
    pInRichLoader->LoadSDigits(); pInRichLoader->TreeS()->GetEntry(0);                          //take list of RICH sdigits from current input 
    AliDebug(1,Form("input %i has %i sdigits",inFileN,pInRich->SDigits()->GetEntries()));
    for(Int_t i=0;i<pInRich->SDigits()->GetEntries();i++){//collect sdigits from current input to tmpCA
      new(tmpCA[total++]) AliRICHDigit(*(AliRICHDigit*)pInRich->SDigits()->At(i)); 
      ((AliRICHDigit*)tmpCA[total-1])->AddTidOffset(fManager->GetMask(inFileN));//apply TID shift since all inputs count tracks independently starting from 0
    }
    pInRichLoader->UnloadSDigits();   pInRich->SDigitsReset(); //close current input and reset 
  }//files loop
  
  tmpCA.Sort();                     //at this point we have a list of all sdigits from all inputs, now sort them according to fPad field
  
  AliRunLoader *pOutRunLoader  = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());    //open output stream (only 1 possible)
  AliLoader    *pOutRichLoader = pOutRunLoader->GetLoader("RICHLoader");                         //take output RICH loader
  AliRICH      *pOutRich       = (AliRICH*)pOutRunLoader->GetAliRun()->GetDetector("RICH");      //take output RICH
  pOutRichLoader->MakeTree("D");   pOutRich->MakeBranch("D");                                    //create TreeD in output stream
  
  TVector pad(2); pad[0]=0; pad[1]=0; Int_t iChamber=0,iCfm=0,aTids[3]={0,0,0},iId=0; Double_t dQdc=0;//current pad info   
  Int_t iNdigsPerPad=0;                   //how many sdigits for a given pad
  for(Int_t i=0;i<tmpCA.GetEntries();i++){//sdigits loop (sorted)
    AliRICHDigit *pSdig=(AliRICHDigit*)tmpCA.At(i);//get new sdigit
    if(pSdig->PadAbs()==iId){//still the same pad
      iNdigsPerPad++;         dQdc+=pSdig->Qdc();      iCfm+=pSdig->Cfm();//sum up charge and cfm
      if(pSdig->Cfm()==1) aTids[0] = pSdig->GetTrack(0); // force the first tid to be mip's tid if it exists in the current pad
      if(iNdigsPerPad<=3)        aTids[iNdigsPerPad-1]=pSdig->GetTrack(0);
      else                         AliDebug(1,Form("More then 3 sdigits for the given pad X:%3.0f Y:%3.0f",pSdig->Pad()(0),pSdig->Pad()(1)));
    }else{//new pad, add the pevious one
        if(iId!=-1 && AliRICHParam::IsOverTh(iChamber,pad,dQdc)) pOutRich->DigitAdd(iChamber,pad,(Int_t)dQdc,iCfm,aTids); //add newly created dig
        iChamber=pSdig->Chamber(); pad=pSdig->Pad(); iCfm=pSdig->Cfm(); dQdc=pSdig->Qdc();  iId=pSdig->PadAbs();                    //init all values by current sdig
        iNdigsPerPad=1; aTids[0]=pSdig->GetTrack(0); aTids[1]=aTids[2]=-1;
      }
  }//sdigits loop (sorted)
  if(tmpCA.GetEntries() && AliRICHParam::IsOverTh(iChamber,pad,dQdc)) pOutRich->DigitAdd(iChamber,pad,(Int_t)dQdc,iCfm,aTids);//add the last dig
  
  pOutRichLoader->TreeD()->Fill();              //fill the output tree with the list of digits
  pOutRichLoader->WriteDigits("OVERWRITE");     //serialize them to file
  
  tmpCA.Clear();                      //remove all tmp sdigits
  pOutRichLoader->UnloadDigits();   pOutRich->DigitsReset();
}//Exec()
