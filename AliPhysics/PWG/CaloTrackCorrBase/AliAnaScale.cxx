/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// Root system
#include <TH1.h>
#include <TH1F.h>

// Analysis system
#include "AliAnaScale.h" 

/// \cond CLASSIMP
ClassImp(AliAnaScale) ;
/// \endcond

//__________________________
/// Default constructor.
//__________________________
AliAnaScale::AliAnaScale() : 
  fDebug(0),
  fScale(1.0),
  fInputList(0x0), 
  fOutputList(0x0),
  fSumw2(0),
  fhCount() 
{
}

//__________________________________________
/// Default constructor.
//__________________________________________
AliAnaScale::AliAnaScale(const char *name) : 
  AliAnalysisTask(name,""),  
  fDebug(0),
  fScale(1.0), 
  fInputList(0x0), 
  fOutputList(0x0), 
  fSumw2(0),
  fhCount(0) 
{  
  // Called only after the event loop
  SetPostEventLoop(kTRUE);
  
  // Input slot #0 
  DefineInput(0,  TList::Class()) ; 
  
  // Output slot 
  DefineOutput(0,  TList::Class()) ;
}

//_________________________________________________
/// Initialisation of branch container with histograms. 
//_________________________________________________
void AliAnaScale::ConnectInputData(const Option_t*)
{
  if(fDebug > 1) printf("*** Initialization of %s \n", GetName()) ; 
  
  fInputList = dynamic_cast<TList*>(GetInputData(0)) ; 
}

//_____________________________________
/// Create the outputs containers.
//_____________________________________
void AliAnaScale::CreateOutputObjects()
{  
  fOutputList = new TList() ; 
  fOutputList->SetName(GetName()) ; 

  fhCount =new TH1F("hCount","count files",1,0,1);  
  fOutputList->Add(fhCount);

  fOutputList->SetOwner(kTRUE);
  
  PostData(0, fOutputList);
}

//________________________________
/// Do the histogram scaling.
//________________________________
void AliAnaScale::Exec(Option_t *) 
{
  if(fDebug > 0 ) printf(">>>>> Scaling factor %e, do Sumw2 %d <<<<< \n",fScale,fSumw2) ;
  
  const Int_t buffersize = 255;
  char name[buffersize] ; 
  
  TIter next(fInputList) ; 	
  TObject * h ; 
  while ( (h = next()) ) 
  { 
    if(h)
    {
      if ( !strncmp(h->ClassName(),"TH",2) ) 
      {
        snprintf(name, buffersize, "%sScaled", h->GetName()) ; 
        
        TH1 * hout = dynamic_cast<TH1*> (h->Clone(name)) ; 
        
        if(hout)
        {
          if(fSumw2) hout->Sumw2();
          hout->Scale(fScale) ;  
          fOutputList->Add(hout) ;
        }// casting not null
      } 
      else  fOutputList->Add(h) ; 
    }
  }
  // number of files
  
  //File scaled, needed for file merging on grid
  fhCount->Fill(0);
  
  PostData(0, fOutputList);
}

//______________________
/// Intialization of parameters.
//______________________
void AliAnaScale::Init()
{  
  if(fDebug > 0 )printf("No initialization in scale class \n") ;
}

