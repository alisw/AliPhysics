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

/*
$Log$
*/

#include <TTree.h> 
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <iostream.h>

#include "AliTPCDigitizer.h"

#include "AliTPC.h"
#include "AliTPCParam.h"
#include "AliRun.h"
#include "AliPDG.h"
#include "AliRunDigitizer.h"
#include "AliSimDigits.h"


ClassImp(AliTPCDigitizer)

//___________________________________________
AliTPCDigitizer::AliTPCDigitizer() :AliDigitizer()
{
// Default ctor - don't use it
}

//___________________________________________
AliTPCDigitizer::AliTPCDigitizer(AliRunDigitizer* manager) 
    :AliDigitizer(manager)
{
// ctor which should be used
  if (GetDebug()>2) 
    cerr<<"AliTPCDigitizer::AliTPCDigitizer"
	<<"(AliRunDigitizer* manager) was processed"<<endl;
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
  
  // merge input tree's with summable digits
  //output stored in TreeTPCD

  TString optionString = option;
  if (optionString.Data() == "deb") {
    cout<<"AliTPCDigitizer::Exec: called with option deb "<<endl;
    fDebug = 3;
  }
  //get detector and geometry
  AliTPC *pTPC  = (AliTPC *) gAlice->GetModule("TPC");
  AliTPCParam * param = pTPC->GetParam();

  printf("noise %f \n",  param->GetNoise()*param->GetNoiseNormFac());
  //create digits array for given sectors
  AliSimDigits ** digarr = new AliSimDigits*[fManager->GetNinputs()]; 
  for (Int_t i1=0;i1<fManager->GetNinputs(); i1++){
    digarr[i1]=0;
    //    intree[i1]
    TTree * treear =  fManager->GetInputTreeTPCS(i1);
    if (!treear) {      
      cerr<<" TPC -  not existing input = \n"<<i1<<" ";      
    }
    treear->GetBranch("Segment")->SetAddress(&digarr[i1]);
  }
  Stat_t nentries = fManager->GetInputTreeTPCS(0)->GetEntries();
  

  //create branch's in TPC treeD
  AliSimDigits * digrow = new AliSimDigits;
  TTree * tree  = fManager->GetTreeDTPC();
  if (tree->GetBranch("Segment") ) tree->GetBranch("Segment")->SetAddress(&digrow);
  else
    tree->Branch("Segment","AliSimDigits",&digrow);
  //

  Int_t zerosup = param->GetZeroSup();
  //Loop over segments of the TPC
    
  for (Int_t n=0; n<nentries; n++) {
    
    for (Int_t i=0;i<fManager->GetNinputs(); i++){ 
      fManager->GetInputTreeTPCS(i)->GetEvent(n);      
      digarr[i]->ExpandBuffer();
      digarr[i]->ExpandTrackBuffer();
    }   
    if ((digarr[0]->GetID()-digarr[1]->GetID())>0) 
      printf("problem\n");
    Int_t sec, row;
    if (!param->AdjustSectorRow(digarr[0]->GetID(),sec,row)) {
      cerr<<"AliTPC warning: invalid segment ID ! "<<digarr[0]->GetID()<<endl;
      continue;
    }

    digrow->SetID(digarr[0]->GetID());

    Int_t nrows = digarr[0]->GetNRows();
    Int_t ncols = digarr[0]->GetNCols();
    digrow->Allocate(nrows,ncols);
    digrow->AllocateTrack(3);

    for (Int_t rows=0;rows<nrows; rows++){
      for (Int_t col=0;col<ncols; col++){
        Float_t q=0;
        Int_t label[1000]; //stack for 300 events 
        Int_t labptr = 0;
        // looop over digits
        for (Int_t i=0;i<fManager->GetNinputs(); i++){ 
          q  += digarr[i]->GetDigitFast(rows,col);
          for (Int_t tr=0;tr<3;tr++) {
            Int_t lab = digarr[i]->GetTrackIDFast(rows,col,tr);
            if ( lab > 1) {
              label[labptr]=lab+fManager->GetMask(i)-2;
              labptr++;
            }
          }
        }
        //add noise
	q/=16.;  //conversion factor
	Float_t noise  = gRandom->Gaus(0,param->GetNoise()*param->GetNoiseNormFac());     
	//printf("%f\n",noise);
	q+=noise;
        q=TMath::Nint(q);
	printf("%f\n",q);

        if (q <= zerosup)continue;
	  
	if(q > param->GetADCSat()) q = (Short_t)(param->GetADCSat());
	digrow->SetDigitFast((Short_t)q,rows,col);  
	
	for (Int_t tr=0;tr<3;tr++){
	  if (tr<labptr)
	    ((AliSimDigits*)digrow)->SetTrackIDFast(label[tr],rows,col,tr);
	  else
	    ((AliSimDigits*)digrow)->SetTrackIDFast(-1,rows,col,tr);          
        }
      } 
    }
    digrow->CompresBuffer(1,zerosup);
    digrow->CompresTrackBuffer(1);
    tree->Fill();
    if (fDebug>0) cerr<<sec<<"\t"<<row<<"\n";  
  }
  delete digrow;     
  for (Int_t ii1=0;ii1<fManager->GetNinputs(); ii1++) delete digarr[ii1];
  delete digarr;  
}

