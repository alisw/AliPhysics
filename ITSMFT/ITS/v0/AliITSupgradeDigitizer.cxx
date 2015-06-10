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
#include "AliDigitizationInput.h"
#include <AliLoader.h>
#include <AliLog.h>
#include "AliITSupgradeDigitizer.h"
#include "AliITSDigitUpgrade.h"
#include "AliITSupgrade.h"
#include "AliITSsegmentationUpgrade.h"
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TTree.h>

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
void AliITSupgradeDigitizer::Digitize(Option_t*)
{
  // This method is responsible for merging sdigits to a list of digits
  //  Disintegration leeds to the fact that one hit affects several neighbouring pads, 
  // which means that the same pad might be affected by few hits.     
  
  AliDebug(1,Form("Start with %i input(s) for event %i",fDigInput->GetNinputs(),fDigInput->GetOutputEventNr()));
  
  
  AliITSsegmentationUpgrade *s = new AliITSsegmentationUpgrade();
  SetConfiguration(s->GetFullCellSizeX(),s->GetFullCellSizeZ());
  delete s;
  //First we read all sdigits from all inputs  
  AliRunLoader *pInRunLoader=0;//in and out Run loaders
  AliLoader    *pITSLoader=0;//in and out ITS loaders  
  
  TClonesArray sdigits[10];
  for(Int_t i=0; i<fNlayers; i++) sdigits[i].SetClass("AliITSDigitUpgrade");//tmp storage for sdigits sum up from all input files
  
  
  Int_t total[10]={0,0,0,0,0,0,0,0,0,0};
  for(Int_t inFileN=0;inFileN<fDigInput->GetNinputs();inFileN++){//files loop
    
    pInRunLoader  = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(inFileN));          //get run loader from current input 
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
	pSDig->AddTidOffset(fDigInput->GetMask(inFileN)); // -> To be introduced for merging (apply TID shift since all inputs count tracks independently starting from 0)
	new((sdigits[is])[total[is]++]) AliITSDigitUpgrade(*pSDig);  
      }
    }
    
    pITSLoader->UnloadSDigits();   
    pITS->SDigitsReset(); //close current input and reset 
  }//files loop  
  
  AliRunLoader *pOutRunLoader  = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());  //open output stream (only 1 possible)
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
   
  Double_t eloss= 0.;
  Double_t nele = 0.;
  
  Int_t tids[maxLab];           // track with id#
  Float_t elossPart[maxLab];    // eloss produced by track with id#
  for(Int_t i=0; i<maxLab ; i++) {
    tids[i]=-1;
    elossPart[i]=-1;
  }
  TArrayI labelsMC(0);          // all track with id#   

  AliDebug(1,"starting loop over layers");
   
  for(Int_t ilay=0;ilay<fNlayers;ilay++){ 

    AliITSDigitUpgrade *tmpdig=0x0;
    pSDigitList[ilay].Sort();

    Int_t prevModule=999; 
    ULong_t  prevPixId = 999;
    Int_t iNdigPart=0;  
    Int_t iNtrackPart = 0;

    AliDebug(1,"starting loop over sdigits to create digits");

    Int_t nSDigits = pSDigitList[ilay].GetEntries();
    Int_t nDigits = 0;

    for(Int_t isdigentr=0; isdigentr<nSDigits; isdigentr++){

      tmpdig = (AliITSDigitUpgrade*)(pSDigitList[ilay].At(isdigentr) )  ;
      
      Int_t    module  = tmpdig->GetModule(); 
      ULong_t  pixId   = tmpdig->GetPixId();
      Int_t    trackId = tmpdig->GetTrackID(0);	 
      
      AliDebug(3,Form("  #tracks %d; #summed digits = %d, TrackIds (%d,%d,%d) ... ",iNtrackPart,iNdigPart,tids[0],tids[1],tids[2]));
      AliDebug(3,Form("\t adding lay:%d  mod:%d  pixId:%lu   trackID:%d \n",tmpdig->GetLayer(),tmpdig->GetModule(),tmpdig->GetPixId(),
		      tmpdig->GetTrackID(0)));
      
      if (trackId<0)
	AliError("Screw you! A track with label<0 produced a SDigit? Something is wrong with Geant?\n");

      if(pixId==prevPixId && module==prevModule) { 
	// This sdigit belongs to the same pixel as before

	// check if trackId is already in the list, if yes just add eloss ...
	Int_t iid =0;
	while (iid<labelsMC.GetSize()) {
	  if ( trackId==labelsMC.At(iid)) {
	    if (iid<maxLab) elossPart[iid] += tmpdig->GetSignal(); // hardcoded limit for elossPart
	    break;
	  }
	  iid++;
	}
	if (iid==labelsMC.GetSize()) { // the trackId is NEW!!
	  if (iid<maxLab) { // hardcoded limits for trackID 
	    tids[iNtrackPart] = tmpdig->GetTrackID(0);
	    elossPart[iNtrackPart] = tmpdig->GetSignal();
	  }
	  iNtrackPart++; 
	  // complete list of trackIDs
	  labelsMC.Set(iNtrackPart); 
	  labelsMC.SetAt(trackId,iNtrackPart-1);
	} else {
	  AliDebug(3,"  -> Track ID already in the list\n");
	}

	if(iNtrackPart>maxLab) {	
	  AliWarning(Form(" Event %d: Number of summable digits for this pixel (lay=%d,mod=%d,pixId=%lu) is too large (%d<%d). Sum is ok but skipping track Id %d ... ",
			  fDigInput->GetOutputEventNr(),ilay,module,pixId,maxLab,iNtrackPart,trackId));
	}
	eloss+=tmpdig->GetSignal();
	nele +=tmpdig->GetNelectrons();
	iNdigPart++;

      } else { // new Pixel
	
	// write "previous" Digit
	if(isdigentr!=0) {
	  AliITSDigitUpgrade digit(prevPixId,eloss);
	  digit.SetNelectrons(nele); 
	  digit.SetLayer(ilay);
	  digit.SetModule(prevModule);
	  digit.SetTids(tids);
	  digit.SetSignalID(elossPart);  
	  digit.SetNTracksIdMC(iNtrackPart);
	  //	  for (Int_t i=0; i<12; i++) { if (i<12)printf("%d ",tids[i]); if (i==11) printf("| ");}
	  //	  for (Int_t i=0; i<digit.GetNTracksIdMC(); i++) { printf("%d ",digit.GetTrackID(i)); if (i==labelsMC.GetSize()-1) printf("\n"); };
	  new((*pLst[ilay])[iCnt[ilay]++]) AliITSDigitUpgrade(digit);
	  nDigits++;
	  AliDebug(3,Form(" -> Wrote NEW digit in layer %d (%d)\n",ilay,nDigits));
	}
	
	// Prepare newly found Pixel
	eloss = tmpdig->GetSignal();
	nele = tmpdig->GetNelectrons();
	
	iNtrackPart=1;  iNdigPart=1;   
	labelsMC.Set(iNtrackPart); 
	labelsMC.SetAt(tmpdig->GetTrackID(0),iNtrackPart-1);
	tids[0]=tmpdig->GetTrackID(0);
	elossPart[0]=tmpdig->GetSignal();
	for(Int_t i=1; i<maxLab ; i++) {
	  tids[i]=-1;
	  elossPart[i]=-1;
	}
      }
      
      prevPixId=tmpdig->GetPixId(); 
      prevModule=tmpdig->GetModule(); 
      
      
    }  
    // write "last" Digit
    AliITSDigitUpgrade digit(prevPixId,eloss);
    digit.SetNelectrons(nele); 
    digit.SetLayer(ilay);
    digit.SetModule(prevModule);
    digit.SetTids(tids);
    digit.SetSignalID(elossPart);      
    digit.SetNTracksIdMC(iNtrackPart);
    //   for (Int_t i=0; i<12; i++) { if (i<12)printf("%d ",tids[i]); if (i==11) printf("| ");}
    //    for (Int_t i=0; i<digit.GetNTracksIdMC(); i++) { printf("%d ",digit.GetTrackID(i)); if (i==labelsMC.GetSize()-1) printf("\n"); };
    new((*pLst[ilay])[iCnt[ilay]++]) AliITSDigitUpgrade(digit);
    nDigits++;
    AliDebug(3,Form(" -> Wrote LAST digit in layer %d (%d)\n",ilay,nDigits));

    AliDebug(3,Form(" -> layer %d: Number of created digits %d",ilay,nDigits));	
    AliDebug(1,"ending loop over sdigits to create digits");

  }
  AliDebug(1,"ending loop over layers");  

}

