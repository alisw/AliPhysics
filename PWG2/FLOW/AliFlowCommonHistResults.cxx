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
 
 /*************************************
 *   AliFlowCommonHistResults:       *
 *   class to organize the common    *
 *   histograms for Flow Analysis    * 
 *                                   * 
 * authors: Naomi van der Kolk       *
 *           (kolk@nikhef.nl)        *  
 *          Raimond Snellings        *
 *           (snelling@nikhef.nl)    * 
 *          Ante Bilandzic           *
 *           (anteb@nikhef.nl)       * 
 * **********************************/

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

ClassImp(AliFlowCommonHistResults)

//-----------------------------------------------------------------------

  AliFlowCommonHistResults::AliFlowCommonHistResults(): 
    TNamed(),
    fHistIntFlow(NULL),//to be removed
    fHistDiffFlow(NULL),//to be removed
    fHistChi(NULL),//to be removed
    fHistIntFlowRP(NULL),
    fHistChiRP(NULL),
    fHistDiffFlowPtRP(NULL),
    fHistDiffFlowEtaRP(NULL),
    fHistIntFlowPOI(NULL),
    fHistChiPOI(NULL),
    fHistDiffFlowPtPOI(NULL),
    fHistDiffFlowEtaPOI(NULL), 
    fHistList(NULL)
{
  //default constructor
} 

//-----------------------------------------------------------------------

  AliFlowCommonHistResults::AliFlowCommonHistResults(const char *anInput, const char *title): 
    TNamed(anInput,title),
    fHistIntFlow(NULL),//to be removed
    fHistDiffFlow(NULL),//to be removed
    fHistChi(NULL),//to be removed
    fHistIntFlowRP(NULL),
    fHistChiRP(NULL),
    fHistDiffFlowPtRP(NULL),
    fHistDiffFlowEtaRP(NULL),
    fHistIntFlowPOI(NULL),
    fHistChiPOI(NULL),
    fHistDiffFlowPtPOI(NULL),
    fHistDiffFlowEtaPOI(NULL),  
    fHistList(NULL)
{
  //constructor creating histograms 
  //Pt:
  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Double_t dPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t dPtMax = AliFlowCommonConstants::GetPtMax();
  //eta:
  Int_t iNbinsEta = AliFlowCommonConstants::GetNbinsEta();
  Double_t dEtaMin = AliFlowCommonConstants::GetEtaMin();	     
  Double_t dEtaMax = AliFlowCommonConstants::GetEtaMax();
  
  TString name;


  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
  //                                 !!!     to be removed    !!!
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
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  //integrated flow
  name = "Flow_Integrated_RP_";
  name += anInput;
  fHistIntFlowRP = new TH1D(name.Data(), name.Data(),1,0.5,1.5);
  fHistIntFlowRP->SetLabelSize(0.06);
  fHistIntFlowRP->SetLabelOffset(0.01);
  (fHistIntFlowRP->GetXaxis())->SetBinLabel(1,"V_{2}");
  
  //Chi (needed for rebinning later on)
  name = "Flow_Chi_RP_";
  name += anInput;
  fHistChiRP = new TH1D(name.Data(), name.Data(),1,0.5,1.5);
  fHistChiRP->SetLabelSize(0.06);
  fHistChiRP->SetLabelOffset(0.01);
  (fHistChiRP->GetXaxis())->SetBinLabel(1,"#chi");
  
  //differential flow (Pt)
  name = "Flow_Differential_Pt_RP_";
  name += anInput;
  fHistDiffFlowPtRP = new TH1D(name.Data(), name.Data(),iNbinsPt,dPtMin,dPtMax);
  fHistDiffFlowPtRP->SetXTitle("P_{t}");
  fHistDiffFlowPtRP->SetYTitle("v_{2}");
  
  //differential flow (eta)
  name = "Flow_Differential_Eta_RP_";
  name += anInput;
  fHistDiffFlowEtaRP = new TH1D(name.Data(), name.Data(),iNbinsEta,dEtaMin,dEtaMax);
  fHistDiffFlowEtaRP->SetXTitle("#eta");
  fHistDiffFlowEtaRP->SetYTitle("v_{2}");
  
  //integrated flow
  name = "Flow_Integrated_POI_";
  name += anInput;
  fHistIntFlowPOI = new TH1D(name.Data(), name.Data(),1,0.5,1.5);
  fHistIntFlowPOI->SetLabelSize(0.06);
  fHistIntFlowPOI->SetLabelOffset(0.01);
  (fHistIntFlowPOI->GetXaxis())->SetBinLabel(1,"V_{2}");
  
  //Chi (needed for rebinning later on)
  name = "Flow_Chi_POI_";
  name += anInput;
  fHistChiPOI = new TH1D(name.Data(), name.Data(),1,0.5,1.5);
  fHistChiPOI->SetLabelSize(0.06);
  fHistChiPOI->SetLabelOffset(0.01);
  (fHistChiPOI->GetXaxis())->SetBinLabel(1,"#chi");
  
  //differential flow (Pt)
  name = "Flow_Differential_Pt_POI_";
  name += anInput;
  fHistDiffFlowPtPOI = new TH1D(name.Data(), name.Data(),iNbinsPt,dPtMin,dPtMax);
  fHistDiffFlowPtPOI->SetXTitle("P_{t}");
  fHistDiffFlowPtPOI->SetYTitle("v_{2}");
  
  //differential flow (eta)
  name = "Flow_Differential_Eta_POI_";
  name += anInput;
  fHistDiffFlowEtaPOI = new TH1D(name.Data(), name.Data(),iNbinsEta,dEtaMin,dEtaMax);
  fHistDiffFlowEtaPOI->SetXTitle("#eta");
  fHistDiffFlowEtaPOI->SetYTitle("v_{2}");
  
  //list of histograms
  fHistList = new TList();
  
  
  fHistList-> Add(fHistIntFlow);//to be removed
  fHistList-> Add(fHistDiffFlow);//to be removed
  fHistList-> Add(fHistChi);//to be removed
  
  
  fHistList->Add(fHistIntFlowRP);
  fHistList->Add(fHistChiRP);
  fHistList->Add(fHistDiffFlowPtRP);
  fHistList->Add(fHistDiffFlowEtaRP);
  fHistList->Add(fHistIntFlowPOI);
  fHistList->Add(fHistChiPOI);
  fHistList->Add(fHistDiffFlowPtPOI);
  fHistList->Add(fHistDiffFlowEtaPOI);  
  
 }

//----------------------------------------------------------------------- 

AliFlowCommonHistResults::~AliFlowCommonHistResults()
{
  //deletes histograms
  delete fHistIntFlow;//to be removed
  delete fHistDiffFlow;//to be removed
  delete fHistChi;//to be removed
  delete fHistIntFlowRP;
  delete fHistChiRP;
  delete fHistDiffFlowPtRP;
  delete fHistDiffFlowEtaRP;
  delete fHistIntFlowPOI;
  delete fHistChiPOI;
  delete fHistDiffFlowPtPOI;
  delete fHistDiffFlowEtaPOI;
  delete fHistList;
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
//                                 !!!     to be removed    !!!
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
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillIntegratedFlowRP(Double_t aV, Double_t anError)
{
  //Fill fHistIntFlowRP
  fHistIntFlowRP->SetBinContent(1,aV);
  fHistIntFlowRP->SetBinError(1,anError);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillChiRP(Double_t aChi)
{
  //Fill fHistChiRP
  fHistChiRP->SetBinContent(1,aChi);
  
  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillDifferentialFlowPtRP(Int_t aBin, Double_t av, Double_t anError)
{
  //Fill fHistDiffFlowPtRP
  fHistDiffFlowPtRP->SetBinContent(aBin,av);
  fHistDiffFlowPtRP->SetBinError(aBin,anError);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillDifferentialFlowEtaRP(Int_t aBin, Double_t av, Double_t anError)
{
  //Fill fHistDiffFlowEtaRP
  fHistDiffFlowEtaRP->SetBinContent(aBin,av);
  fHistDiffFlowEtaRP->SetBinError(aBin,anError);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillIntegratedFlowPOI(Double_t aV, Double_t anError)
{
  //Fill fHistIntFlowPOI
  fHistIntFlowPOI->SetBinContent(1,aV);
  fHistIntFlowPOI->SetBinError(1,anError);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillChiPOI(Double_t aChi)
{
  //Fill fHistChiPOI
  fHistChiPOI->SetBinContent(1,aChi);
  
  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillDifferentialFlowPtPOI(Int_t aBin, Double_t av, Double_t anError)
{
  //Fill fHistDiffFlowPtPOI
  fHistDiffFlowPtPOI->SetBinContent(aBin,av);
  fHistDiffFlowPtPOI->SetBinError(aBin,anError);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillDifferentialFlowEtaPOI(Int_t aBin, Double_t av, Double_t anError)
{
  //Fill fHistDiffFlowEtaPOI
  fHistDiffFlowEtaPOI->SetBinContent(aBin,av);
  fHistDiffFlowEtaPOI->SetBinError(aBin,anError);

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




