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
Bool_t AliRICHDigitizer::Init()
{
//This methode is called from AliRunDigitizer after the corresponding file is open
  if(GetDebug())Info("Init","Start.");
  AliRunLoader *pOutAL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  if (!pOutAL->GetAliRun()) pOutAL->LoadgAlice();
  fRich=(AliRICH*)pOutAL->GetAliRun()->GetDetector("RICH");//retrive RICH pointer from OUTPUT stream
  return kTRUE;
}//Init()
//__________________________________________________________________________________________________
void AliRICHDigitizer::Exec(Option_t*)
{
//this method invoked   
  if(GetDebug())Info("Exec","Start with %i input(s) for event %i",fManager->GetNinputs(),fManager->GetOutputEventNr());
  
  AliRunLoader *pInAL=0, *pOutAL;//in and out Run loaders
  AliLoader    *pInRL=0, *pOutRL;//in and out RICH loaders
 
  pOutAL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  pOutRL = pOutAL->GetLoader("RICHLoader");
  pOutRL->MakeTree("D");   R()->MakeBranch("D"); //create TreeD with RICH branches in output stream
  
  TClonesArray tmpCA("AliRICHdigit");//tmp storage for sdigits sum up from all input files
  Int_t total=0;
  for(Int_t inFileN=0;inFileN<fManager->GetNinputs();inFileN++){//files loop
    pInAL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inFileN));//get run loader from current input 
    pInRL = pInAL->GetLoader("RICHLoader"); if(pInRL==0) continue;//no RICH in this input, check the next input
    if (!pInAL->GetAliRun()) pInAL->LoadgAlice();
    AliRICH* rich=(AliRICH*)pInAL->GetAliRun()->GetDetector("RICH");
    pInRL->LoadSDigits(); pInRL->TreeS()->GetEntry(0);
    if(GetDebug())Info("Exec","input %i has %i sdigits",inFileN,rich->SDigits()->GetEntries());
    for(Int_t i=0;i<rich->SDigits()->GetEntries();i++){//collect sdigits from current input to tmpCA
      new(tmpCA[total++]) AliRICHdigit(*(AliRICHdigit*)rich->SDigits()->At(i)); 
      ((AliRICHdigit*)tmpCA[total-1])->AddTidOffset(fManager->GetMask(inFileN));
    }
    pInRL->UnloadSDigits();   rich->ResetSDigits();
  }//files loop
  
  tmpCA.Sort();                     //sort them according to Id() method
  
  if(GetDebug()) {tmpCA.Print();Info("Exec","Totally %i sdigits in %i inputs",tmpCA.GetEntries(),fManager->GetNinputs());}
  Int_t cfm=0,chamber=0,id=0; //cfm is cerenkov feedback mip mixture
  TVector pad(2); pad[0]=0;pad[1]=0;
  Double_t q=0;
  Int_t tid[3]={0,0,0};
  Int_t iNdigitsPerPad=0;//how many sdigits for a given pad
  for(Int_t i=0;i<tmpCA.GetEntries();i++){//sdigits loop (sorted)
    AliRICHdigit *pSdig=(AliRICHdigit*)tmpCA.At(i);//get new sdigit
    if(pSdig->Id()==id){//still the same pad
      iNdigitsPerPad++;         q+=pSdig->Q();       cfm+=pSdig->ChFbMi();//sum up charge and cfm
      if(pSdig->ChFbMi()==1) tid[0] = pSdig->GetTrack(0); // force the first tid to be mip's tid if it exists in the current pad
      if(iNdigitsPerPad<=3)        tid[iNdigitsPerPad-1]=pSdig->GetTrack(0);
      else                         if(GetDebug())Warning("Exec","More then 3 sdigits for the given pad");
    }else{//new pad, add the pevious one
        if(id!=kBad&&R()->P()->IsOverTh(chamber,pad,q)) R()->AddDigit(chamber,pad,(Int_t)q,cfm,tid); //add newly created dig
        cfm=pSdig->ChFbMi(); chamber=pSdig->C(); id=pSdig->Id();  pad=pSdig->Pad(); q=pSdig->Q();  //init all values by current sdig
        iNdigitsPerPad=1; tid[0]=pSdig->GetTrack(0); tid[1]=tid[2]=kBad;
      }
  }//sdigits loop (sorted)
  if(tmpCA.GetEntries()&&R()->P()->IsOverTh(chamber,pad,q)) R()->AddDigit(chamber,pad,(Int_t)q,cfm,tid);//add the last dig
  
  pOutRL->TreeD()->Fill();              //fill the tree with the list of digits
  pOutRL->WriteDigits("OVERWRITE");     //serialize them to file
  
  tmpCA.Clear();
  pOutRL->UnloadDigits();   R()->ResetDigits();
            
  if(GetDebug())Info("Exec","Stop.");
}//Exec()
//__________________________________________________________________________________________________
