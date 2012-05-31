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
//.
//.
//.
//.
#include <AliRun.h>
#include <AliRunLoader.h>
#include "AliDigitizationInput.h"
#include <AliLoader.h>
#include "AliLog.h"
#include <AliCDBEntry.h> 
#include <AliCDBManager.h>
#include "AliHMPIDDigitizer.h"
#include "AliHMPIDReconstructor.h"
#include "AliHMPIDDigit.h"
#include "AliHMPID.h"
#include "AliHMPIDParam.h"
#include <TRandom.h>
#include <TMath.h>
#include <TTree.h>
#include <TObjArray.h>

ClassImp(AliHMPIDDigitizer)

Bool_t AliHMPIDDigitizer::fgDoNoise=kTRUE;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigitizer::Digitize(Option_t*)
{
// This methode is responsible for merging sdigits to a list of digits
//Disintegration leeds to the fact that one hit affects several neighbouring pads, which means that the same pad might be affected by few hits.     
  AliDebug(1,Form("Start with %i input(s) for event %i",fDigInput->GetNinputs(),fDigInput->GetOutputEventNr()));
//First we read all sdigits from all inputs  
  AliRunLoader *pInRunLoader=0;//in and out Run loaders
  AliLoader    *pInRichLoader=0;//in and out HMPID loaders  
  static TClonesArray sdigs("AliHMPIDDigit");//tmp storage for sdigits sum up from all input files
  Int_t total=0;
  for(Int_t inFileN=0;inFileN<fDigInput->GetNinputs();inFileN++){//files loop
    pInRunLoader  = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(inFileN));          //get run loader from current input 
    pInRichLoader = pInRunLoader->GetLoader("HMPIDLoader"); if(pInRichLoader==0) continue;       //no HMPID in this input, check the next input
    if (!pInRunLoader->GetAliRun()) pInRunLoader->LoadgAlice();
    AliHMPID* pInRich=(AliHMPID*)pInRunLoader->GetAliRun()->GetDetector("HMPID");                  //take HMPID from current input
    pInRichLoader->LoadSDigits(); pInRichLoader->TreeS()->GetEntry(0);                          //take list of HMPID sdigits from current input 
    AliDebug(1,Form("input %i has %i sdigits",inFileN,pInRich->SdiLst()->GetEntries()));
    for(Int_t i=0;i<pInRich->SdiLst()->GetEntries();i++){                                        //collect sdigits from current input
      AliHMPIDDigit *pSDig=(AliHMPIDDigit*)pInRich->SdiLst()->At(i);
      pSDig->AddTidOffset(fDigInput->GetMask(inFileN));                                          //apply TID shift since all inputs count tracks independently starting from 0
      new(sdigs[total++]) AliHMPIDDigit(*pSDig);       
    }
    pInRichLoader->UnloadSDigits();   pInRich->SdiReset(); //close current input and reset 
  }//files loop

  //PH  if(sdigs.GetEntries()==0) return;                                                              //no sdigits collected, nothing to convert  
  
  AliRunLoader *pOutRunLoader  = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());    //open output stream (only 1 possible)
  AliLoader    *pOutRichLoader = pOutRunLoader->GetLoader("HMPIDLoader");                         //take output HMPID loader
  AliRun *pArun = pOutRunLoader->GetAliRun();
  AliHMPID      *pOutRich       = (AliHMPID*)pArun->GetDetector("HMPID");      //take output HMPID
  pOutRichLoader->MakeTree("D");   pOutRich->MakeBranch("D");                                    //create TreeD in output stream

  Sdi2Dig(&sdigs,pOutRich->DigLst());
  
  pOutRichLoader->TreeD()->Fill();              //fill the output tree with the list of digits
  pOutRichLoader->WriteDigits("OVERWRITE");     //serialize them to file
  
  sdigs.Clear();                      //remove all tmp sdigits
  pOutRichLoader->UnloadDigits();   pOutRich->DigReset();
}//Exec()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigitizer::Sdi2Dig(TClonesArray *pSdiLst,TObjArray *pDigLst)
{
// Converts list of sdigits to 7 lists of digits, one per each chamber
// Arguments: pSDigLst - list of all sdigits
//            pDigLst  - list of 7 lists of digits        
//   Returns: none  
  
  TClonesArray *pLst[7]; Int_t iCnt[7];

  for(Int_t i=0;i<7;i++){
    pLst[i]=(TClonesArray*)(*pDigLst)[i];
    iCnt[i]=0; if(pLst[i]->GetEntries()!=0) AliErrorClass("Some of digits lists is not empty");         //in principle those lists should be empty                                                                       
  }
  
  // make noise array
  Float_t arrNoise[7][6][80][48], arrSigmaPed[7][6][80][48];
  if(fgDoNoise) {
    for (Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++)
      for (Int_t iPc=AliHMPIDParam::kMinPc;iPc<=AliHMPIDParam::kMaxPc;iPc++)
        for(Int_t iPx=AliHMPIDParam::kMinPx;iPx<=AliHMPIDParam::kMaxPx;iPx++)
          for(Int_t iPy=AliHMPIDParam::kMinPy;iPy<=AliHMPIDParam::kMaxPy;iPy++){
            arrNoise[iCh][iPc][iPx][iPy] = gRandom->Gaus(0,1.);
            arrSigmaPed[iCh][iPc][iPx][iPy] = 1.;
          }
          
    AliCDBEntry *pDaqSigEnt = AliCDBManager::Instance()->Get("HMPID/Calib/DaqSig");  //contains TObjArray of TObjArray 14 TMatrixF sigmas values for pads 
   
    if(pDaqSigEnt){
      TObjArray *pDaqSig = (TObjArray*)pDaqSigEnt->GetObject();
      for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){                  //chambers loop    
	TMatrixF *pM = (TMatrixF*)pDaqSig->At(iCh);     
	for (Int_t iPc=AliHMPIDParam::kMinPc;iPc<=AliHMPIDParam::kMaxPc;iPc++)
	  for(Int_t iPx=AliHMPIDParam::kMinPx;iPx<=AliHMPIDParam::kMaxPx;iPx++)
	    for(Int_t iPy=AliHMPIDParam::kMinPy;iPy<=AliHMPIDParam::kMaxPy;iPy++){
	      Int_t padX = (iPc%2)*AliHMPIDParam::kPadPcX+iPx; 
	      Int_t padY = (iPc/2)*AliHMPIDParam::kPadPcY+iPy;
	      if((*pM)(padX,padY)>0.){
		arrNoise[iCh][iPc][iPx][iPy] = gRandom->Gaus(0,(*pM)(padX,padY));
		arrSigmaPed[iCh][iPc][iPx][iPy] = (*pM)(padX,padY);}
	      else{
		arrNoise[iCh][iPc][iPx][iPy] = gRandom->Gaus(0,1.);
		arrSigmaPed[iCh][iPc][iPx][iPy] = 1.;}
	    } 
      }
    }
  }  
  
  pSdiLst->Sort();  
                     
  Int_t iPad=-1,iCh=-1,iNdigPad=-1,aTids[3]={-1,-1,-1}; Float_t q=-1;
  for(Int_t i=0;i<pSdiLst->GetEntries();i++){                                                                  //sdigits loop (sorted)
    AliHMPIDDigit *pSdig=(AliHMPIDDigit*)pSdiLst->At(i);                                                       //take current sdigit
    if(pSdig->Pad()==iPad){                                                                                    //if the same pad 
      q+=pSdig->Q();                                                                                           //sum up charge
      iNdigPad++; if(iNdigPad<=3) aTids[iNdigPad-1]=pSdig->GetTrack(0);                                        //collect TID 
      continue;
    }
    if(i!=0 && iCh>=AliHMPIDParam::kMinCh && iCh<=AliHMPIDParam::kMaxCh){
      AliHMPIDParam::Instance()->SetThreshold((TMath::Nint(arrSigmaPed[iCh][pSdig->Pc()][pSdig->PadPcX()][pSdig->PadPcY()])*AliHMPIDParam::Nsig()));
      if(AliHMPIDParam::IsOverTh(q)) new((*pLst[iCh])[iCnt[iCh]++]) AliHMPIDDigit(iPad,(Int_t)q,aTids);}  //do not create digit for the very first sdigit 
    
    iPad=pSdig->Pad(); iCh=AliHMPIDParam::A2C(iPad);                                                            //new sdigit comes, reset collectors
    iNdigPad=1;
    aTids[0]=pSdig->GetTrack(0);aTids[1]=aTids[2]=-1; 
    q=pSdig->Q();
    if(fgDoNoise) q+=arrNoise[iCh][pSdig->Pc()][pSdig->PadPcX()][pSdig->PadPcY()];
    arrNoise[iCh][pSdig->Pc()][pSdig->PadPcX()][pSdig->PadPcY()]=0;
  }//sdigits loop (sorted)
  
  if(iCh>=AliHMPIDParam::kMinCh && iCh<=AliHMPIDParam::kMaxCh){
    Int_t pc = AliHMPIDParam::A2P(iPad);
    Int_t px = AliHMPIDParam::A2X(iPad); 
    Int_t py = AliHMPIDParam::A2Y(iPad);
    AliHMPIDParam::Instance()->SetThreshold((TMath::Nint(arrSigmaPed[iCh][pc][px][py])*AliHMPIDParam::Nsig()));
    if(AliHMPIDParam::IsOverTh(q)) new((*pLst[iCh])[iCnt[iCh]++]) AliHMPIDDigit(iPad,(Int_t)q,aTids);
  }  //add the last one, in case of empty sdigits list q=-1 
  
// add noise pad above threshold with no signal merged...if any
  if(!fgDoNoise) return;
  
  aTids[0]=aTids[1]=aTids[2]=-1;
  for (Int_t iChCurr=AliHMPIDParam::kMinCh;iChCurr<=AliHMPIDParam::kMaxCh;iChCurr++){
    for (Int_t iPc=AliHMPIDParam::kMinPc;iPc<=AliHMPIDParam::kMaxPc;iPc++)
      for(Int_t iPx=AliHMPIDParam::kMinPx;iPx<=AliHMPIDParam::kMaxPx;iPx++)
        for(Int_t iPy=AliHMPIDParam::kMinPy;iPy<=AliHMPIDParam::kMaxPy;iPy++) {
          Float_t qNoise = arrNoise[iChCurr][iPc][iPx][iPy];
          AliHMPIDParam::Instance()->SetThreshold((TMath::Nint(arrSigmaPed[iChCurr][iPc][iPx][iPy])*AliHMPIDParam::Nsig()));
          if(AliHMPIDParam::IsOverTh(qNoise)) new((*pLst[iChCurr])[iCnt[iChCurr]++]) AliHMPIDDigit(AliHMPIDParam::Abs(iChCurr,iPc,iPx,iPy),(Int_t)qNoise,aTids);
        }
  }        
}//Sdi2Dig()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
