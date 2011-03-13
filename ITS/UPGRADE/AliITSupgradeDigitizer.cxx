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
//.

/* $Id$ */

#include <AliRun.h>
#include <AliRunLoader.h>
#include "AliRunDigitizer.h"
#include <AliLoader.h>
#include <AliLog.h>
#include "AliITSupgradeDigitizer.h"
#include "AliITSDigitUpgrade.h"
#include "AliITSupgrade.h"
#include "AliITSsegmentationUpgrade.h"
#include <TObjArray.h>
#include <TClonesArray.h>

extern TRandom *gRandom;

ClassImp(AliITSupgradeDigitizer)
    
  void AliITSupgradeDigitizer::SetConfiguration(TArrayD xcell, TArrayD zcell)
{
  
  if(xcell.GetSize()!=zcell.GetSize()) AliError(" !! The # of X cells and Z cells differ !!");
  
  fNlayers = xcell.GetSize();
  
  if(fNlayers > 9) {
    AliError("*  Only 9 layers can be be filled  ...Exiting!!! *");
    return;
  }
  
  fNxCells.Set(fNlayers);
  fNzCells.Set(fNlayers);
  for(Int_t i=0; i<fNlayers; i++){
    fNxCells.AddAt(xcell.At(i),i);
    fNzCells.AddAt(zcell.At(i),i); 
  }  
}   
//______________________________________________________________________________      
void AliITSupgradeDigitizer::Exec(Option_t*)
{
  // This method is responsible for merging sdigits to a list of digits
  //  Disintegration leeds to the fact that one hit affects several neighbouring pads, 
  // which means that the same pad might be affected by few hits.     
  
  AliDebug(1,Form("Start with %i input(s) for event %i",fManager->GetNinputs(),fManager->GetOutputEventNr()));
  
  
  AliITSsegmentationUpgrade *s = new AliITSsegmentationUpgrade();
  SetConfiguration(s->GetFullCellSizeX(),s->GetFullCellSizeZ());
  delete s;
  //First we read all sdigits from all inputs  
  AliRunLoader *pInRunLoader=0;//in and out Run loaders
  AliLoader    *pITSLoader=0;//in and out ITS loaders  
  
  TClonesArray sdigits[10];
  for(Int_t i=0; i<fNlayers; i++) sdigits[i].SetClass("AliITSDigitUpgrade");//tmp storage for sdigits sum up from all input files
  
  
  Int_t total[10]={0,0,0,0,0,0,0,0,0,0};
  for(Int_t inFileN=0;inFileN<fManager->GetNinputs();inFileN++){//files loop
    
    pInRunLoader  = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inFileN));          //get run loader from current input 
    pITSLoader = pInRunLoader->GetLoader("ITSLoader"); 
    if(pITSLoader==0) {
      continue;       //no ITS in this input, check the next input
      AliDebug(1,"no ITS lodader, checking in the other input \n"); 
    }
    
    if (!pInRunLoader->GetAliRun()) pInRunLoader->LoadgAlice();
    AliITSupgrade* pITS=(AliITSupgrade*)pInRunLoader->GetAliRun()->GetDetector("ITS"); 
    
    pITSLoader->LoadSDigits();  
    
    pITSLoader->TreeS()->GetEntry(0);                          //take list of ITS sdigits from current input 
    
    for(Int_t is=0;is<pITS->SDigitsList()->GetEntries();is++){      
      
      //collect sdigits from current input
      for(Int_t ientr =0; ientr < ((TClonesArray*)pITS->SDigitsList()->At(is))->GetEntries(); ientr++){
	AliITSDigitUpgrade *pSDig=(AliITSDigitUpgrade*)((TClonesArray*)pITS->SDigitsList()->At(is))->At(ientr);
	pSDig->AddTidOffset(fManager->GetMask(inFileN)); // -> To be introduced for merging (apply TID shift since all inputs count tracks independently starting from 0)
	new((sdigits[is])[total[is]++]) AliITSDigitUpgrade(*pSDig);  
      }
    }
    
    pITSLoader->UnloadSDigits();   
    pITS->SDigitsReset(); //close current input and reset 
  }//files loop  
  
  AliRunLoader *pOutRunLoader  = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());  //open output stream (only 1 possible)
  AliLoader    *pOutITSLoader = pOutRunLoader->GetLoader("ITSLoader");                        
  AliRun *pArun = pOutRunLoader->GetAliRun();
  AliITSupgrade      *pOutITS       = (AliITSupgrade*)pArun->GetDetector("ITS");      
  pOutITSLoader->MakeTree("D");   pOutITS->MakeBranch("D");                                    //create TreeD in output stream
  pOutITS->SetTreeAddress();
  
  Sdigits2Digits(sdigits,pOutITS->DigitsList());
  
  pOutITSLoader->TreeD()->Fill();              //fill the output tree with the list of digits
  pOutITSLoader->WriteDigits("OVERWRITE");     //serialize them to file
  
  for(Int_t i=0; i< fNlayers; i++) sdigits[i].Clear();                      //remove all tmp sdigits
  pOutITSLoader->UnloadDigits();   
  pOutITS->DigitsReset(); 
}//Exec()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliITSupgradeDigitizer::Sdigits2Digits(TClonesArray *pSDigitList,TObjArray *pDigitList)
{   
  TClonesArray *pLst[100]; Int_t iCnt[100];
 
  for(Int_t i=0;i<fNlayers;i++){ 
    pLst[i]=(TClonesArray*)(*pDigitList)[i];
    iCnt[i]=0; if(pLst[i]->GetEntries()!=0) AliErrorClass("Some of digits lists is not empty");  //in principle those lists should be empty 
  }
   
  //AliInfo("starting loop over gli sdigits to create the digits");
  Double_t eloss =0.;
  Double_t nele = 0.;
  ULong_t  pixid = 999;
  Int_t tids[3]={-1,-1,-1};
  
  Float_t elossID[3]={-1.,-1.,-1.};
  //AliInfo("starting layers");
  AliDebug(1,"starting loop over layers");
   
  for(Int_t ilay=0;ilay<fNlayers;ilay++){ 
     
    AliITSDigitUpgrade *tmpdig=0x0;
    pSDigitList[ilay].Sort();
    Int_t module=999; 
    Int_t iNdigPart=0; 
    AliDebug(1,"starting loop over sdigits to create digits");
    for(Int_t isdigentr=0; isdigentr<pSDigitList[ilay].GetEntries(); isdigentr++){
      tmpdig = (AliITSDigitUpgrade*)(pSDigitList[ilay].At(isdigentr) )  ;
     if(tmpdig->GetPixId()==pixid && tmpdig->GetModule()==module) {
	iNdigPart++; 
	if(iNdigPart<=3) {
          tids[iNdigPart-1] = tmpdig->GetTrack(0);
          elossID[iNdigPart-1] = tmpdig->GetSignal();
	}     
	eloss+=tmpdig->GetSignal();
	nele+=tmpdig->GetNelectrons();
	continue;
      }
      AliITSDigitUpgrade digit(pixid,eloss);
      
      digit.SetNelectrons(nele); 
      digit.SetLayer(ilay);
      digit.SetModule(module);
      digit.SetTids(tids);
      digit.SetSignalID(elossID);       
      if(isdigentr!=0) new((*pLst[ilay])[iCnt[ilay]++]) AliITSDigitUpgrade(digit);
      eloss = tmpdig->GetSignal();
      nele = tmpdig->GetNelectrons();
      pixid=tmpdig->GetPixId(); 
      tids[0]=tmpdig->GetTrack(0);
      tids[1]=tids[2]=-1;      
      elossID[0]=tmpdig->GetSignal();
      elossID[1]=elossID[2]=-1;
      module=tmpdig->GetModule();   
    }
     
    if(!tmpdig) AliDebug(1,"\n \n---------> tmpdig is null...break is expected ! \n");
    else AliDebug(1," tmpdig exists \n");
     
    if(tmpdig){
      tmpdig->SetSignal(eloss);  
      tmpdig->SetPixId(pixid); 
      tmpdig->SetTids(tids); 
      tmpdig->SetSignalID(elossID);
      tmpdig->SetNelectrons(nele);  
      tmpdig->SetLayer(ilay); 
      tmpdig->SetModule(module);
      //cout<<" tmpdigit : pixid "<< pixid<< "  tids "<< tids << " nele " << nele << " ilay "<<ilay<<endl;     
      new((*pLst[ilay])[iCnt[ilay]++]) AliITSDigitUpgrade(*tmpdig);
    }     
    AliDebug(1,"ending loop over sdigits to create digits");

  }
  AliDebug(1,"ending loop over layers");  

}

