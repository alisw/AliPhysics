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

#include "Riostream.h"                 //needed as include
#include "AliFlowCommonConstants.h"    //needed as include
#include "AliFlowCommonHistResults.h"

#include "TString.h" 
#include "TH1D.h"   //needed as include
#include "TMath.h"  //needed as include
#include "TList.h"
#include "TBrowser.h"

class TH1F;
class AliFlowVector;
class AliFlowCommonHist;

// AliFlowCommonHistResults:
// Class to organize the common histograms for Flow Analysis
// Holds v2(pt), integrated v2 and chi (resolution)
//
// authors: N.K A.B. R.S

ClassImp(AliFlowCommonHistResults)

//-----------------------------------------------------------------------

  AliFlowCommonHistResults::AliFlowCommonHistResults(): 
    TNamed(),
    fHistIntFlow(0),
    fHistDiffFlow(0),
    fHistChi(0),
    fHistList(NULL)
{
  //default constructor
} 

//-----------------------------------------------------------------------

  AliFlowCommonHistResults::AliFlowCommonHistResults(const char *anInput,const char *title): 
    TNamed(anInput,title),
    fHistIntFlow(0),
    fHistDiffFlow(0),
    fHistChi(0),
    fHistList(NULL)
{
  //constructor creating histograms 
  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  TString name;

  Double_t  dPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  dPtMax = AliFlowCommonConstants::GetPtMax();
  
  //integrated flow
  name = "Flow_Integrated_";
  name += anInput;
  fHistIntFlow = new TH1D(name.Data(), name.Data(),1,0.5,1.5);
  fHistIntFlow ->SetXTitle("");
  fHistIntFlow ->SetYTitle("V_{2}");

  //differential flow
  name = "Flow_Differential_Pt_";
  name += anInput;
  fHistDiffFlow = new TH1D(name.Data(), name.Data(),iNbinsPt,dPtMin,dPtMax);
  fHistDiffFlow ->SetXTitle("P_{t}");
  fHistDiffFlow ->SetYTitle("v_{2}");
  
  //Chi (needed for rebinning later on)
  name = "Flow_Chi_";
  name += anInput;
  fHistChi = new TH1D(name.Data(), name.Data(),1,0.5,1.5);
  fHistChi ->SetXTitle("");
  fHistChi ->SetYTitle("#Chi");

  //list of histograms
  fHistList = new TList();
  fHistList-> Add(fHistIntFlow);
  fHistList-> Add(fHistDiffFlow);
  fHistList-> Add(fHistChi);

  }

//----------------------------------------------------------------------- 

AliFlowCommonHistResults::~AliFlowCommonHistResults()
{
  //deletes histograms
  delete fHistIntFlow;
  delete fHistDiffFlow;
  delete fHistChi;
  delete fHistList;
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillIntegratedFlow(Double_t aV, Double_t anError)
{
  //Fill fHistIntFlow
  fHistIntFlow -> SetBinContent(1,aV);
  fHistIntFlow -> SetBinError(1,anError);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillDifferentialFlow(Int_t aBin, Double_t av, Double_t anError)
{
  //Fill fHistDiffFlow
  fHistDiffFlow ->SetBinContent(aBin,av);
  fHistDiffFlow ->SetBinError(aBin,anError);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillChi(Double_t aChi)
{
  //Fill fHistChi
  fHistChi -> SetBinContent(1,aChi);
  
  return kTRUE; 
}

//----------------------------------------------------------------------- 
 Double_t AliFlowCommonHistResults::Merge(TCollection *aList)
{
  //merge fuction
  cout<<"entering merge function"<<endl;
  if (!aList) return 0;
  if (aList->IsEmpty()) return 0; //no merging is needed

  Int_t iCount = 0;
  TIter next(aList); // list is supposed to contain only objects of the same type as this
  AliFlowCommonHistResults *toMerge;
  // make a temporary list
  TList *pTemp = new TList();
  while ((toMerge=(AliFlowCommonHistResults*)next())) {
    pTemp->Add(toMerge->GetHistList()); 
    iCount++;
  }
  // Now call merge for fHistList providing temp list
  fHistList->Merge(pTemp);
  // Cleanup
  delete pTemp;
    
  cout<<"Merged"<<endl;
  return (double)iCount;
    
}

//----------------------------------------------------------------------- 
void AliFlowCommonHistResults::Print(Option_t *option) const
{
  //   -*-*-*-*-*Print some global quantities for this histogram collection class *-*-*-*-*-*-*-*
  //             ===============================================
  //   printf( "TH1.Print Name  = %s, Entries= %d, Total sum= %g\n",GetName(),Int_t(fEntries),GetSumOfWeights());
  printf( "Class.Print Name = %s, Histogram list:\n",GetName());

  if (fHistList) {  
    fHistList->Print(option);
  }
  else
    {
      printf( "Empty histogram list \n");
    }
}

//----------------------------------------------------------------------- 
 void AliFlowCommonHistResults::Browse(TBrowser *b)
{

  if (!b) return;
  if (fHistList) b->Add(fHistList,"AliFlowCommonHistResultsList");
}




