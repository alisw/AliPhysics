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
    fHistIntFlow(NULL),
    fHistChi(NULL),
    fHistIntFlowRP(NULL),
    fHistDiffFlowPtRP(NULL),
    fHistDiffFlowEtaRP(NULL),
    fHistIntFlowPOI(NULL),
    fHistDiffFlowPtPOI(NULL),
    fHistDiffFlowEtaPOI(NULL), 
    fHistList(NULL)
{
  //default constructor
} 

//-----------------------------------------------------------------------

  AliFlowCommonHistResults::AliFlowCommonHistResults(const char *anInput, const char *title, Int_t harmonic): 
    TNamed(anInput,title),
    fHistIntFlow(NULL),
    fHistChi(NULL),
    fHistIntFlowRP(NULL),
    fHistDiffFlowPtRP(NULL),
    fHistDiffFlowEtaRP(NULL),
    fHistIntFlowPOI(NULL),
    fHistDiffFlowPtPOI(NULL),
    fHistDiffFlowEtaPOI(NULL),  
    fHistList(NULL)
{
  //constructor creating histograms 
  //Pt:
  Int_t iNbinsPt = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  Double_t dPtMin = AliFlowCommonConstants::GetMaster()->GetPtMin();	     
  Double_t dPtMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
  //eta:
  Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
  Double_t dEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();	     
  Double_t dEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();
  
  TString name;

  // Reference flow: (TBI: rename eventually integrated flow => reference flow)
  name = "Flow_Integrated_";
  name += anInput;
  fHistIntFlow = new TH1D(name.Data(),"Reference Flow",1,0.5,1.5);
  fHistIntFlow->SetStats(kFALSE);
  fHistIntFlow->SetMarkerStyle(kOpenSquare);
  fHistIntFlow->SetLabelSize(0.06,"X");
  fHistIntFlow->SetLabelOffset(0.015,"X");
  fHistIntFlow->GetXaxis()->SetBinLabel(1,Form("v_{%d}",harmonic));
 
  // chi (resolution):
  name = "Flow_Chi_";
  name += anInput;
  fHistChi = new TH1D(name.Data(),"Resolution",1,0.5,1.5);
  fHistChi->SetStats(kFALSE);
  fHistChi->SetMarkerStyle(kOpenSquare);
  fHistChi->SetLabelSize(0.06,"X");
  fHistChi->SetLabelOffset(0.015,"X");
  fHistChi->GetXaxis()->SetBinLabel(1,"#chi");
  
  // Integrated flow of RPs:
  name = "Flow_Integrated_RP_";
  name += anInput;
  fHistIntFlowRP = new TH1D(name.Data(),"Integrated Flow (RP)",1,0.5,1.5);
  fHistIntFlowRP->SetStats(kFALSE);
  fHistIntFlowRP->SetMarkerStyle(kOpenSquare);
  fHistIntFlowRP->SetLabelSize(0.06,"X");
  fHistIntFlowRP->SetLabelOffset(0.015,"X");
  fHistIntFlowRP->GetXaxis()->SetBinLabel(1,Form("v_{%d}",harmonic));

  // Differential flow (Pt) of RPs:
  name = "Flow_Differential_Pt_RP_";
  name += anInput;
  fHistDiffFlowPtRP = new TH1D(name.Data(),"Differential Flow vs p_{t} (RP)",iNbinsPt,dPtMin,dPtMax);
  fHistDiffFlowPtRP->SetStats(kFALSE);
  fHistDiffFlowPtRP->SetXTitle("p_{t}");
  fHistDiffFlowPtRP->SetYTitle(Form("v_{%d}",harmonic));
  
  // Differential flow (eta) of RPs:
  name = "Flow_Differential_Eta_RP_";
  name += anInput;
  fHistDiffFlowEtaRP = new TH1D(name.Data(),"Differential Flow vs #eta (RP)",iNbinsEta,dEtaMin,dEtaMax);
  fHistDiffFlowEtaRP->SetStats(kFALSE);
  fHistDiffFlowEtaRP->SetXTitle("#eta");
  fHistDiffFlowEtaRP->SetYTitle(Form("v_{%d}",harmonic));
  
  // Integrated flow of POIs:
  name = "Flow_Integrated_POI_";
  name += anInput;
  fHistIntFlowPOI = new TH1D(name.Data(),"Integrated Flow (POI)",1,0.5,1.5);
  fHistIntFlowPOI->SetStats(kFALSE);
  fHistIntFlowPOI->SetMarkerStyle(kOpenSquare);
  fHistIntFlowPOI->SetLabelSize(0.06,"X");
  fHistIntFlowPOI->SetLabelOffset(0.015,"X");
  fHistIntFlowPOI->GetXaxis()->SetBinLabel(1,Form("v_{%d}",harmonic));

  // Differential flow (Pt) of POIs:
  name = "Flow_Differential_Pt_POI_";
  name += anInput;
  fHistDiffFlowPtPOI = new TH1D(name.Data(),"Differential Flow vs p_{t} (POI)",iNbinsPt,dPtMin,dPtMax);
  fHistDiffFlowPtPOI->SetXTitle("p_{t}");
  fHistDiffFlowPtPOI->SetYTitle(Form("v_{%d}",harmonic));
  
  // Differential flow (eta) of POIs:
  name = "Flow_Differential_Eta_POI_";
  name += anInput;
  fHistDiffFlowEtaPOI = new TH1D(name.Data(),"Differential Flow vs #eta (POI)",iNbinsEta,dEtaMin,dEtaMax);
  fHistDiffFlowEtaPOI->SetStats(kFALSE);
  fHistDiffFlowEtaPOI->SetXTitle("#eta");
  fHistDiffFlowEtaPOI->SetYTitle(Form("v_{%d}",harmonic));
  
  // List of histograms:
  fHistList = new TList();
  fHistList-> Add(fHistIntFlow);
  fHistList-> Add(fHistChi);
  fHistList->Add(fHistIntFlowRP);
  fHistList->Add(fHistDiffFlowPtRP);
  fHistList->Add(fHistDiffFlowEtaRP);
  fHistList->Add(fHistIntFlowPOI);
  fHistList->Add(fHistDiffFlowPtPOI);
  fHistList->Add(fHistDiffFlowEtaPOI);    
 }

//----------------------------------------------------------------------- 

AliFlowCommonHistResults::~AliFlowCommonHistResults()
{
 // Deletes histograms:
 delete fHistIntFlow;
 delete fHistChi;
 delete fHistIntFlowRP;
 delete fHistDiffFlowPtRP;
 delete fHistDiffFlowEtaRP;
 delete fHistIntFlowPOI;
 delete fHistDiffFlowPtPOI;
 delete fHistDiffFlowEtaPOI;
 delete fHistList;
}

//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHistResults::FillIntegratedFlow(Double_t aV, Double_t anError)
{
 // Fill fHistIntFlow:
 fHistIntFlow -> SetBinContent(1,aV);
 fHistIntFlow -> SetBinError(1,anError);

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

Bool_t AliFlowCommonHistResults::FillIntegratedFlowRP(Double_t aV, Double_t anError)
{
  //Fill fHistIntFlowRP
  fHistIntFlowRP->SetBinContent(1,aV);
  fHistIntFlowRP->SetBinError(1,anError);

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
  //cout<<"entering merge function"<<endl;
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
    
  //cout<<"Merged"<<endl;
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




