/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Class AliAnalysisBtoJPSItoEle
//                  Unbinned log-likelihood fit analysis class
//
//                             Origin: C.Di Giglio
//        Contact: Carmelo.Digiglio@ba.infn.it , giuseppe.bruno@ba.infn.it
//-------------------------------------------------------------------------
class TH1F ;
#include "TNtuple.h"
#include "TMath.h"

#include "AliBtoJPSItoEleCDFfitFCN.h"
#include "AliBtoJPSItoEleCDFfitHandler.h"
#include "AliAnalysisBtoJPSItoEle.h"
#include "AliLog.h"

ClassImp(AliAnalysisBtoJPSItoEle)

//_______________________________________________________________________________ 
AliAnalysisBtoJPSItoEle::AliAnalysisBtoJPSItoEle() :
fFCNfunction(0),
fPtBin(0),
fMCtemplate(0)
{
  //
  // default constructor
  //
}
//___________________________________________________________________________________
AliAnalysisBtoJPSItoEle::AliAnalysisBtoJPSItoEle(const AliAnalysisBtoJPSItoEle& source) :
TNamed(source),
fFCNfunction(source.fFCNfunction),
fPtBin(source.fPtBin),
fMCtemplate(source.fMCtemplate)
{
  //
  // copy constructor
  //
}
//_________________________________________________________________________________________________

AliAnalysisBtoJPSItoEle &AliAnalysisBtoJPSItoEle::operator=(const AliAnalysisBtoJPSItoEle& source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;
  fFCNfunction = source.fFCNfunction;
  fPtBin = source.fPtBin;
  fMCtemplate = source.fMCtemplate;

  return *this;
}
//_________________________________________________________________________________________________
AliAnalysisBtoJPSItoEle::~AliAnalysisBtoJPSItoEle()
{
  //
  // destructor
  //
  delete fFCNfunction;
  delete fMCtemplate;
}
//_________________________________________________________________________________________________
Int_t AliAnalysisBtoJPSItoEle::DoMinimization(Double_t* x,
                                              Double_t* m, Int_t ncand)
{
  //
  // performs the minimization
  //
  AliInfo(Form("Number of candidates used for the minimisation is %d",ncand));
  fFCNfunction = new AliBtoJPSItoEleCDFfitHandler(x,m,ncand);
  SetResolutionConstants(fPtBin);
  SetCsiMC(fMCtemplate);
  fFCNfunction->SetErrorDef(0.5); // tells Minuit that the error interval is the one in which
                                  // the function differs from the minimum for less than setted value
  Int_t iret=fFCNfunction->DoMinimization();

  return iret;
}
//_________________________________________________________________________________________________
void AliAnalysisBtoJPSItoEle::SetResolutionConstants(Int_t BinNum)
{
  //
  // sets constants for parametrized resolution function
  //
  if(!fFCNfunction) {
    AliInfo("fFCNfunction not istanziated  ---> nothing done");
    return;
  }
  AliInfo("Call likelihood SetResolutionConstants method ---> OK");
  AliBtoJPSItoEleCDFfitFCN* loglikePnt = fFCNfunction->LikelihoodPointer();
  if(!loglikePnt) {
     AliWarning("Pointer to AliBtoJPSItoEleCDFfitFCN class not found!");
     return;
    }
  loglikePnt->SetResolutionConstants(BinNum);
}
//_________________________________________________________________________________________________
void AliAnalysisBtoJPSItoEle::ReadCandidates(TNtuple* nt, Double_t* &pseudoproper,Double_t* &invmass, Int_t& ncand)
{
  //
  // Read N-tuple with X and M values
  //
  Float_t mJPSI = 0; Float_t x = 0;
  Int_t nentries = 0;
  ncand=0;
  nt->SetBranchAddress("Mass",&mJPSI);
  nt->SetBranchAddress("Xdecaytime",&x);
  nentries = (Int_t)nt->GetEntries();
  pseudoproper = new Double_t[nentries];
  invmass      = new Double_t[nentries];
  for(Int_t i = 0; i < nentries; i++) {
      nt->GetEntry(i);
      ncand++;
      pseudoproper[i]=(Double_t)(10000*x);
      invmass[i]=(Double_t)mJPSI;
    }

 return; 
}
//_________________________________________________________________________________________________
void AliAnalysisBtoJPSItoEle::SetCsiMC(TH1F* MCtemplates)
{
  //
  // Sets X distribution used as MC template for JPSI from B
  //
  fFCNfunction->LikelihoodPointer()->SetCsiMC(MCtemplates);

  return;
}
//_________________________________________________________________________________________________
