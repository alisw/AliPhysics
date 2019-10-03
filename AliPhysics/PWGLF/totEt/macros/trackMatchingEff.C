#ifndef __CINT__
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include <iostream>
#endif

int trackMatchingEff(TString file = "Et.ESD.simPbPb.PHOS.root")
{
  TFile *f = TFile::Open(file, "READ");
  if(!f)
  {
    std::cout << "Could not open file: " << file << " !" << std::endl;
    return -1;
  }
  
  TList *l = (TList*)(f->Get("out1"));
  
  if(!l)
  {
    std::cout << "Could not find list!" << std::endl;
    return -1;
  }
  
  TTree *primTree = (TTree*)(l->FindObject("fPrimaryTreePhosMC"));
  
  if(!primTree)
  {
    std::cout << "Could not find tree!" << std::endl;
    return -1;
  }
  
  
  TString emSelect = "(fPrimaryCode==22||fPrimaryCode==221||TMath::Abs(fPrimaryCode)==11)";
  TString chargeSelect = "(fPrimaryCharge!=0 && TMath::Abs(fPrimaryCode)!=11)";
  TString neutralSelect = "(!"+emSelect+")&&fPrimaryCharge==0&&(!fSecondary)"; 
  TString secondarySelect = "(fSecondary)";
  emSelect += "&&(!fSecondary)";
  chargeSelect += "&&(!fSecondary)";
  
  TString matchedSelect = "fPrimaryMatched==1&&";
  TString notMatchedSelect = "fPrimaryMatched==0&&";
  
  
  int n = primTree->Draw("fDepositedEt", notMatchedSelect+ emSelect);
  int nRemoved = primTree->Draw("fDepositedEt", matchedSelect + emSelect);
  std::cout << "EM: " << float(n)/(n+nRemoved) << std::endl;
  
  n = primTree->Draw("fDepositedEt", notMatchedSelect+ chargeSelect);
  nRemoved = primTree->Draw("fDepositedEt", matchedSelect + chargeSelect);
  std::cout << "Charged: " << float(n)/(n+nRemoved) << std::endl;

  n = primTree->Draw("fDepositedEt", notMatchedSelect+ neutralSelect);
  nRemoved = primTree->Draw("fDepositedEt", matchedSelect + neutralSelect);
  std::cout << "Neutral: " << float(n)/(n+nRemoved) << std::endl;

  n = primTree->Draw("fDepositedEt", notMatchedSelect+ secondarySelect);
  nRemoved = primTree->Draw("fDepositedEt", matchedSelect + secondarySelect);
  if(n+nRemoved) std::cout << "Secondary: " << float(n)/(n+nRemoved) << std::endl;
  else std::cout << "No secondaries" << std::endl;
  return 0;
  
  
}