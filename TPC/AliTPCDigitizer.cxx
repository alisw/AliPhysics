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

/* $Id$ */

#include <stdlib.h>
#include <TTree.h> 
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <Riostream.h>

#include "AliTPCDigitizer.h"

#include "AliTPC.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h" 
#include "AliRun.h"
#include "AliPDG.h"
#include "AliRunDigitizer.h"
#include "AliSimDigits.h"
#include "AliLog.h"

#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"

ClassImp(AliTPCDigitizer)

//___________________________________________
  AliTPCDigitizer::AliTPCDigitizer() :AliDigitizer(),fDebug(0)
{
  //
// Default ctor - don't use it
//
  
}

//___________________________________________
AliTPCDigitizer::AliTPCDigitizer(AliRunDigitizer* manager) 
  :AliDigitizer(manager),fDebug(0)
{
  //
// ctor which should be used
//  
  AliDebug(2,"(AliRunDigitizer* manager) was processed");
}

//------------------------------------------------------------------------
AliTPCDigitizer::~AliTPCDigitizer()
{
// Destructor
}



//------------------------------------------------------------------------
Bool_t AliTPCDigitizer::Init()
{
// Initialization 
    
 return kTRUE;
}


//------------------------------------------------------------------------
void AliTPCDigitizer::Exec(Option_t* option)
{
  ExecFast(option);  
}
//------------------------------------------------------------------------
void AliTPCDigitizer::ExecFast(Option_t* option)
{
  
  // merge input tree's with summable digits
  //output stored in TreeTPCD
  char s[100]; 
  char ss[100];
  TString optionString = option;
  if (optionString.Data() == "deb") {
    cout<<"AliTPCDigitizer::Exec: called with option deb "<<endl;
    fDebug = 3;
  }
  //get detector and geometry


  AliRunLoader *rl, *orl;
  AliLoader *gime, *ogime;
  
  if (gAlice == 0x0)
   {
     Warning("ExecFast","gAlice is NULL. Loading from input 0");
     rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
     if (rl == 0x0)
      {
        Error("ExecFast","Can not find Run Loader for input 0. Can not proceed.");
        return;
      }
     rl->LoadgAlice();
     rl->GetAliRun();
   }
  AliTPC *pTPC  = (AliTPC *) gAlice->GetModule("TPC");
  AliTPCParam * param = pTPC->GetParam();
  
  sprintf(s,param->GetTitle());
  sprintf(ss,"75x40_100x60");
  if(strcmp(s,ss)==0){
    printf("2 pad-length geom hits with 3 pad-lenght geom digits...\n");
    delete param;
    param=new AliTPCParamSR();
  }
  else{
   sprintf(ss,"75x40_100x60_150x60");
   if(strcmp(s,ss)!=0) {
     printf("No TPC parameters found...\n");
     exit(2); 
   }
  }
  
  pTPC->GenerNoise(500000); //create teble with noise
  //
  Int_t nInputs = fManager->GetNinputs();
  Int_t * masks = new Int_t[nInputs];
  for (Int_t i=0; i<nInputs;i++)
    masks[i]= fManager->GetMask(i);
  Short_t **pdig= new Short_t*[nInputs];   //pointers to the expanded digits array
  Int_t **ptr=  new Int_t*[nInputs];       //pointers to teh expanded tracks array
  Bool_t *active=  new Bool_t[nInputs];    //flag for active input segments

  
  //create digits array for given sectors
  // make indexes
  
  AliSimDigits ** digarr = new AliSimDigits*[nInputs]; 
  for (Int_t i1=0;i1<nInputs; i1++)
    {
      digarr[i1]=0;
     //    intree[i1]
      rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(i1));
      gime = rl->GetLoader("TPCLoader");
      gime->LoadSDigits("read");
      TTree * treear =  gime->TreeS();
     
      if (!treear) 
       {
        cerr<<"AliTPCDigitizer: Input tree with SDigits not found in"
            <<" input "<< i1<<endl;
        return;
       }
    
      if (treear->GetIndex()==0)  
	treear->BuildIndex("fSegmentID","fSegmentID");
      treear->GetBranch("Segment")->SetAddress(&digarr[i1]);
    }


  //create branch's in TPC treeD
  AliSimDigits * digrow = new AliSimDigits;

  orl = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  ogime = orl->GetLoader("TPCLoader");
  
  TTree * tree  = ogime->TreeD();
  if (tree == 0x0)
   {
     ogime->MakeTree("D");
     tree  = ogime->TreeD();
   }
  tree->Branch("Segment","AliSimDigits",&digrow);
  //

  param->SetZeroSup(2);

  Int_t zerosup = param->GetZeroSup(); 
  AliTPCCalPad * gainTPC = AliTPCcalibDB::Instance()->GetPadGainFactor(); 
  //
  //Loop over segments of the TPC
    
  for (Int_t segmentID=0; segmentID<param->GetNRowsTotal(); segmentID++) 
   {
    Int_t sec, row;
    if (!param->AdjustSectorRow(segmentID,sec,row)) 
     {
      cerr<<"AliTPC warning: invalid segment ID ! "<<segmentID<<endl;
      continue;
     }
    AliTPCCalROC * gainROC = gainTPC->GetCalROC(sec);  // pad gains per given sector
    digrow->SetID(segmentID);

    Int_t nrows = 0;
    Int_t ncols = 0;

    Bool_t digitize = kFALSE;
    for (Int_t i=0;i<nInputs; i++) 
     { 

      rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(i));
      gime = rl->GetLoader("TPCLoader");
      
      if (gime->TreeS()->GetEntryWithIndex(segmentID,segmentID) >= 0) {
	digarr[i]->ExpandBuffer();
	digarr[i]->ExpandTrackBuffer();
        nrows = digarr[i]->GetNRows();
        ncols = digarr[i]->GetNCols();
	active[i] = kTRUE;
	if (!fRegionOfInterest || (i == 0)) digitize = kTRUE;
      } else {
	active[i] = kFALSE;
      }
      if (fRegionOfInterest && !digitize) break;
     }   
    if (!digitize) continue;

    digrow->Allocate(nrows,ncols);
    digrow->AllocateTrack(3);

    Float_t q=0;
    Int_t label[1000]; //stack for 300 events 
    Int_t labptr = 0;

    Int_t nElems = nrows*ncols;     
 
    for (Int_t i=0;i<nInputs; i++)
     if (active[i]) { 
       pdig[i] = digarr[i]->GetDigits();
       ptr[i]  = digarr[i]->GetTracks();
      }
     
    Short_t *pdig1= digrow->GetDigits();
    Int_t   *ptr1= digrow->GetTracks() ;

    

    for (Int_t elem=0;elem<nElems; elem++)
     {    

       q=0;
       labptr=0;
       // looop over digits 
        for (Int_t i=0;i<nInputs; i++) if (active[i]) 
         { 
          //          q  += digarr[i]->GetDigitFast(rows,col);
            q  += *(pdig[i]);
         
           for (Int_t tr=0;tr<3;tr++) 
            {
             //             Int_t lab = digarr[i]->GetTrackIDFast(rows,col,tr);
             Int_t lab = ptr[i][tr*nElems];
             if ( (lab > 1) && *(pdig[i])>zerosup) 
              {
                label[labptr]=lab+masks[i];
                labptr++;
              }          
            }
           pdig[i]++;
           ptr[i]++;
         }
        q/=16.;  //conversion factor
	Float_t gain = gainROC->GetValue(row,elem/nrows);  // get gain for given - pad-row pad
	if (gain<0.5){
	  printf("problem\n");
	}
	q*= gain;
        //       Float_t noise  = gRandom->Gaus(0,param->GetNoise()*param->GetNoiseNormFac());  
        Float_t noise  = pTPC->GetNoise();
        q+=noise;
        q=TMath::Nint(q);
        if (q > zerosup)
         { 
          if(q >= param->GetADCSat()) q = (Short_t)(param->GetADCSat() - 1);
          //digrow->SetDigitFast((Short_t)q,rows,col);  
          *pdig1 =Short_t(q);
          for (Int_t tr=0;tr<3;tr++)
           {
            if (tr<labptr) 
             ptr1[tr*nElems] = label[tr];
           }
          }
        pdig1++;
        ptr1++;
     }
    
    digrow->CompresBuffer(1,zerosup);
    digrow->CompresTrackBuffer(1);
    tree->Fill();
    if (fDebug>0) cerr<<sec<<"\t"<<row<<"\n";  
   } //for (Int_t n=0; n<param->GetNRowsTotal(); n++) 
  

  orl = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  ogime = orl->GetLoader("TPCLoader");
  ogime->WriteDigits("OVERWRITE");
  
  //fManager->GetTreeDTPC()->Write(0,TObject::kOverwrite);
  
  delete digrow;     
  for (Int_t i1=0;i1<nInputs; i1++) delete digarr[i1];
  delete []masks;
  delete []pdig;
  delete []ptr;
  delete []active;
  delete []digarr;  
}



//------------------------------------------------------------------------
void AliTPCDigitizer::ExecSave(Option_t* option)
{
  
  // merge input tree's with summable digits
  //output stored in TreeTPCD

  TString optionString = option;
  if (optionString.Data() == "deb") {
    cout<<"AliTPCDigitizer::Exec: called with option deb "<<endl;
    fDebug = 3;
  }
  //get detector and geometry 
  AliRunLoader *rl, *orl;
  AliLoader *gime, *ogime;

  
  orl = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  ogime = orl->GetLoader("TPCLoader");
  
  rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
  gime = rl->GetLoader("TPCLoader");
  
  rl->LoadgAlice();
  AliRun* alirun = rl->GetAliRun();
  
  AliTPC *pTPC  = (AliTPC *) alirun->GetModule("TPC");
  AliTPCParam * param = pTPC->GetParam();
  pTPC->GenerNoise(500000); //create teble with noise
  printf("noise %f \n",  param->GetNoise()*param->GetNoiseNormFac());
  //
  Int_t nInputs = fManager->GetNinputs();
  Int_t * masks = new Int_t[nInputs];
  for (Int_t i=0; i<nInputs;i++)
    masks[i]= fManager->GetMask(i);

  AliSimDigits ** digarr = new AliSimDigits*[nInputs]; 
  for (Int_t i1=0;i1<nInputs; i1++)
   {
    digarr[i1]=0;
    //    intree[i1]
    rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(i1));
    gime = rl->GetLoader("TPCLoader");

    TTree * treear =  gime->TreeS();
    TBranch * br = treear->GetBranch("fSegmentID");
    if (br) br->GetFile()->cd();
    if (!treear) {      
      cerr<<" TPC -  not existing input = \n"<<i1<<" ";      
    } 
    treear->GetBranch("Segment")->SetAddress(&digarr[i1]);
  }
  
  rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
  gime = rl->GetLoader("TPCLoader");
  Stat_t nentries = gime->TreeS()->GetEntries();
  

  //create branch's in TPC treeD
  AliSimDigits * digrow = new AliSimDigits;
  TTree * tree  = ogime->TreeD();

  tree->Branch("Segment","AliSimDigits",&digrow);
  param->SetZeroSup(2);

  Int_t zerosup = param->GetZeroSup();
  //Loop over segments of the TPC
    
  AliTPCCalPad * gainTPC = AliTPCcalibDB::Instance()->GetPadGainFactor();
  for (Int_t n=0; n<nentries; n++) {
    rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
    gime = rl->GetLoader("TPCLoader");
    gime->TreeS()->GetEvent(n);

    digarr[0]->ExpandBuffer();
    digarr[0]->ExpandTrackBuffer();


    for (Int_t i=1;i<nInputs; i++){ 
//      fManager->GetInputTreeTPCS(i)->GetEntryWithIndex(digarr[0]->GetID(),digarr[0]->GetID());      
      rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(i));
      gime = rl->GetLoader("TPCLoader");
      gime->TreeS()->GetEntryWithIndex(digarr[0]->GetID(),digarr[0]->GetID());  
      digarr[i]->ExpandBuffer();
      digarr[i]->ExpandTrackBuffer();
      if ((digarr[0]->GetID()-digarr[i]->GetID())>0) 
       printf("problem\n");
    
    }   
    
    Int_t sec, row;
    if (!param->AdjustSectorRow(digarr[0]->GetID(),sec,row)) {
      cerr<<"AliTPC warning: invalid segment ID ! "<<digarr[0]->GetID()<<endl;
      continue;
    }

    AliTPCCalROC * gainROC = gainTPC->GetCalROC(sec);  // pad gains per given sector
    digrow->SetID(digarr[0]->GetID());

    Int_t nrows = digarr[0]->GetNRows();
    Int_t ncols = digarr[0]->GetNCols();
    digrow->Allocate(nrows,ncols);
    digrow->AllocateTrack(3);

    Float_t q=0;
    Int_t label[1000]; //stack for 300 events 
    Int_t labptr = 0;

    

    for (Int_t rows=0;rows<nrows; rows++){
      for (Int_t col=0;col<ncols; col++){
    
       q=0;
       labptr=0;
       // looop over digits 
        for (Int_t i=0;i<nInputs; i++){ 
         q  += digarr[i]->GetDigitFast(rows,col);
          //q  += *(pdig[i]);
         
          for (Int_t tr=0;tr<3;tr++) {
           Int_t lab = digarr[i]->GetTrackIDFast(rows,col,tr);
           //Int_t lab = ptr[i][tr*nElems];
            if ( (lab > 1) ) {
              label[labptr]=lab+masks[i];
              labptr++;
            }          
          }
         // pdig[i]++;
         //ptr[i]++;
         
        }
       q/=16.;  //conversion factor
       //       Float_t noise  = gRandom->Gaus(0,param->GetNoise()*param->GetNoiseNormFac());  
       Float_t gain = gainROC->GetValue(row,col);
       q*= gain;
       Float_t noise  = pTPC->GetNoise();
       q+=noise;
        q=TMath::Nint(q);
        if (q > zerosup){ 
         
         if(q >= param->GetADCSat()) q = (Short_t)(param->GetADCSat() - 1);
         digrow->SetDigitFast((Short_t)q,rows,col);  
         // *pdig1 =Short_t(q);
         for (Int_t tr=0;tr<3;tr++){
           if (tr<labptr) 
             ((AliSimDigits*)digrow)->SetTrackIDFast(label[tr],rows,col,tr);
           //ptr1[tr*nElems] = label[tr];
           //else
             //           ((AliSimDigits*)digrow)->SetTrackIDFast(-1,rows,col,tr);          
           //  ptr1[tr*nElems] = 1;
         }
       }
       //pdig1++;
       //ptr1++;
    }
    }
    
    digrow->CompresBuffer(1,zerosup);
    digrow->CompresTrackBuffer(1);
    tree->Fill();
    if (fDebug>0) cerr<<sec<<"\t"<<row<<"\n";  
  } 
//  printf("end TPC merging - end -Tree %s\t%p\n",fManager->GetInputTreeH(0)->GetName(),fManager->GetInputTreeH(0)->GetListOfBranches()->At(3));
  //fManager->GetTreeDTPC()->Write(0,TObject::kOverwrite);
    ogime->WriteDigits("OVERWRITE");

    for (Int_t i=1;i<nInputs; i++) 
     { 
      rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(i));
      gime = rl->GetLoader("TPCLoader");
      gime->UnloadSDigits();
     }
    ogime->UnloadDigits();
    
  delete digrow;     
  for (Int_t i1=0;i1<nInputs; i1++) delete digarr[i1];
  delete []masks;
  delete digarr;  
}
