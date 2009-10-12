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
Int_t AliAnalysisBtoJPSItoEle::DoMinimization()
{
  //
  // performs the minimization
  //
  Int_t iret=fFCNfunction->DoMinimization();

  return iret;
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
void AliAnalysisBtoJPSItoEle::SetCsiMC()
{
  //
  // Sets X distribution used as MC template for JPSI from B
  //
  fFCNfunction->LikelihoodPointer()->SetCsiMC(fMCtemplate);

  return;
}
//_________________________________________________________________________________________________
void AliAnalysisBtoJPSItoEle::SetFitHandler(Double_t* x /*pseudoproper*/, Double_t* m /*inv mass*/, Int_t ncand /*candidates*/) 
{
  //
  // Create the fit handler object to play with different params of the fitting function
  //

  fFCNfunction = new AliBtoJPSItoEleCDFfitHandler(x,m,ncand);
  if(!fFCNfunction) {

     AliInfo("fFCNfunction not istanziated  ---> nothing done");
     return;

     } 
}
