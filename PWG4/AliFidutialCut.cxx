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
/* $Id: AliFidutialCut.cxx 21839 2007-10-29 13:49:42Z gustavo $ */

//_________________________________________________________________________
// Class for track/cluster acceptance selection
// Selection in Central barrel, EMCAL and PHOS
//                
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TMath.h>
#include <TLorentzVector.h>
#include <TString.h>
//#include <TArrayF.h>
  
//---- ANALYSIS system ----
#include "AliLog.h"
#include "AliFidutialCut.h"

ClassImp(AliFidutialCut)


//____________________________________________________________________________
AliFidutialCut::AliFidutialCut() : 
  TObject(),
  fEMCALFidutialCut(0),  fPHOSFidutialCut(0),  fCTSFidutialCut(0),
  fCTSFidCutMinEta(0x0),fCTSFidCutMinPhi(0x0),fCTSFidCutMaxEta(0x0), fCTSFidCutMaxPhi(0x0),
  fEMCALFidCutMinEta(0x0),fEMCALFidCutMinPhi(0x0),fEMCALFidCutMaxEta(0x0), fEMCALFidCutMaxPhi(0x0),
  fPHOSFidCutMinEta(0x0),fPHOSFidCutMinPhi(0x0),fPHOSFidCutMaxEta(0x0), fPHOSFidCutMaxPhi(0x0)

{
  //Ctor

  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliFidutialCut::AliFidutialCut(const AliFidutialCut & g) :   
  TObject(g), 
  fEMCALFidutialCut(g.fEMCALFidutialCut),  fPHOSFidutialCut(g.fPHOSFidutialCut),  fCTSFidutialCut(g. fCTSFidutialCut),
  fCTSFidCutMinEta(g.fCTSFidCutMinEta?new TArrayF(*g.fCTSFidCutMinEta):0x0),
  fCTSFidCutMinPhi(g.fCTSFidCutMinPhi?new TArrayF(*g.fCTSFidCutMinPhi):0x0),
  fCTSFidCutMaxEta(g.fCTSFidCutMaxEta?new TArrayF(*g.fCTSFidCutMaxEta):0x0),
  fCTSFidCutMaxPhi(g.fCTSFidCutMaxPhi?new TArrayF(*g.fCTSFidCutMaxPhi):0x0),
  fEMCALFidCutMinEta(g.fEMCALFidCutMinEta?new TArrayF(*g.fEMCALFidCutMinEta):0x0),
  fEMCALFidCutMinPhi(g.fEMCALFidCutMinPhi?new TArrayF(*g.fEMCALFidCutMinPhi):0x0),
  fEMCALFidCutMaxEta(g.fEMCALFidCutMaxEta?new TArrayF(*g.fEMCALFidCutMaxEta):0x0),
  fEMCALFidCutMaxPhi(g.fEMCALFidCutMaxPhi?new TArrayF(*g.fEMCALFidCutMaxPhi):0x0),
  fPHOSFidCutMinEta(g.fPHOSFidCutMinEta?new TArrayF(*g.fPHOSFidCutMinEta):0x0),
  fPHOSFidCutMinPhi(g.fPHOSFidCutMinPhi?new TArrayF(*g.fPHOSFidCutMinPhi):0x0),
  fPHOSFidCutMaxEta(g.fPHOSFidCutMaxEta?new TArrayF(*g.fPHOSFidCutMaxEta):0x0),
  fPHOSFidCutMaxPhi(g.fPHOSFidCutMaxPhi?new TArrayF(*g.fPHOSFidCutMaxPhi):0x0)

{
  // cpy ctor

}

//_________________________________________________________________________
AliFidutialCut & AliFidutialCut::operator = (const AliFidutialCut & source)
{
  // assignment operator
  
  if(&source == this) return *this;
  
  fEMCALFidutialCut = source.fEMCALFidutialCut;  
  fPHOSFidutialCut = source.fPHOSFidutialCut;
  fCTSFidutialCut = source.fCTSFidutialCut;
  
  fCTSFidCutMinEta = source.fCTSFidCutMinEta?new TArrayF(*source.fCTSFidCutMinEta):0x0;
  fCTSFidCutMinPhi = source.fCTSFidCutMinPhi?new TArrayF(*source.fCTSFidCutMinPhi):0x0;
  fCTSFidCutMaxEta = source.fCTSFidCutMaxEta?new TArrayF(*source.fCTSFidCutMaxEta):0x0;
  fCTSFidCutMaxPhi = source.fCTSFidCutMaxPhi?new TArrayF(*source.fCTSFidCutMaxPhi):0x0;
  fEMCALFidCutMinEta = source.fEMCALFidCutMinEta?new TArrayF(*source.fEMCALFidCutMinEta):0x0;
  fEMCALFidCutMinPhi = source.fEMCALFidCutMinPhi?new TArrayF(*source.fEMCALFidCutMinPhi):0x0;
  fEMCALFidCutMaxEta = source.fEMCALFidCutMaxEta?new TArrayF(*source.fEMCALFidCutMaxEta):0x0;
  fEMCALFidCutMaxPhi = source.fEMCALFidCutMaxPhi?new TArrayF(*source.fEMCALFidCutMaxPhi):0x0;
  fPHOSFidCutMinEta = source.fPHOSFidCutMinEta?new TArrayF(*source.fPHOSFidCutMinEta):0x0;
  fPHOSFidCutMinPhi = source.fPHOSFidCutMinPhi?new TArrayF(*source.fPHOSFidCutMinPhi):0x0;
  fPHOSFidCutMaxEta = source.fPHOSFidCutMaxEta?new TArrayF(*source.fPHOSFidCutMaxEta):0x0;
  fPHOSFidCutMaxPhi = source.fPHOSFidCutMaxPhi?new TArrayF(*source.fPHOSFidCutMaxPhi):0x0;

  return *this;

}

//_________________________________
AliFidutialCut::~AliFidutialCut() {
  //Dtor

  delete fCTSFidCutMinEta ;
  delete fCTSFidCutMinPhi ;
  delete fCTSFidCutMaxEta ;
  delete fCTSFidCutMaxPhi ;
  
  delete fEMCALFidCutMinEta ;
  delete fEMCALFidCutMinPhi ;
  delete fEMCALFidCutMaxEta ; 
  delete fEMCALFidCutMaxPhi ;
  
  delete fPHOSFidCutMinEta ; 
  delete fPHOSFidCutMinPhi ; 
  delete fPHOSFidCutMaxEta ;
  delete fPHOSFidCutMaxPhi ;

}


//_______________________________________________________________
Bool_t AliFidutialCut::IsInFidutialCut(TLorentzVector momentum, TString det) const
{
  //Selects EMCAL or PHOS cluster or CTS track if it is inside eta-phi defined regions

  Double_t phi = momentum.Phi();
  if(phi<0) phi+=TMath::TwoPi() ;
  Double_t eta =  momentum.Eta();
  //printf("IsInFidutialCut::Det: %s, phi = %f, eta = %f\n", det.Data(),phi*TMath::RadToDeg(), eta);

  Float_t * maxeta = new Float_t;
  Float_t * maxphi = new Float_t;
  Float_t * mineta = new Float_t;
  Float_t * minphi = new Float_t;
  Int_t nphiregions = 0;
  Int_t netaregions = 0;
  Bool_t selection = kFALSE;

  if(det == "CTS"){
    netaregions =  fCTSFidCutMaxEta->GetSize();
    nphiregions =  fCTSFidCutMaxPhi->GetSize();
    if(netaregions !=  fCTSFidCutMinEta->GetSize() || nphiregions !=  fCTSFidCutMinPhi->GetSize())
      AliFatal(Form("Wrong number of CTS fidutial cut regions: nmaxeta %d != nmineta %d; nmaxphi %d != nminphi %d",
		    netaregions, fCTSFidCutMinEta->GetSize(),  nphiregions, fCTSFidCutMinPhi->GetSize()));
    
    maxeta = fCTSFidCutMaxEta->GetArray();
    maxphi = fCTSFidCutMaxPhi->GetArray();
    mineta = fCTSFidCutMinEta->GetArray();
    minphi = fCTSFidCutMinPhi->GetArray();
    selection =  fCTSFidutialCut ; 
  }
  else   if(det == "EMCAL"){
    netaregions =  fEMCALFidCutMaxEta->GetSize();
    nphiregions =  fEMCALFidCutMaxPhi->GetSize();
    if(netaregions !=  fEMCALFidCutMinEta->GetSize() || nphiregions !=  fEMCALFidCutMinPhi->GetSize())
      AliFatal(Form("Wrong number of EMCAL fidutial cut regions: nmaxeta %d != nmineta %d; nmaxphi %d != nminphi %d",
		    netaregions, fEMCALFidCutMinEta->GetSize(),  nphiregions, fEMCALFidCutMinPhi->GetSize()));
    
    maxeta = fEMCALFidCutMaxEta->GetArray();
    maxphi = fEMCALFidCutMaxPhi->GetArray();
    mineta = fEMCALFidCutMinEta->GetArray();
    minphi = fEMCALFidCutMinPhi->GetArray();
    selection =  fEMCALFidutialCut ; 
  }
  else   if(det == "PHOS"){
    netaregions =  fPHOSFidCutMaxEta->GetSize();
    nphiregions =  fPHOSFidCutMaxPhi->GetSize();
    if(netaregions !=  fPHOSFidCutMinEta->GetSize() || nphiregions !=  fPHOSFidCutMinPhi->GetSize())
      AliFatal(Form("Wrong number of PHOS fidutial cut regions: nmaxeta %d != nmineta %d; nmaxphi %d != nminphi %d",
		    netaregions, fPHOSFidCutMinEta->GetSize(),  nphiregions, fPHOSFidCutMinPhi->GetSize()));
    
    maxeta = fPHOSFidCutMaxEta->GetArray();
    maxphi = fPHOSFidCutMaxPhi->GetArray();
    mineta = fPHOSFidCutMinEta->GetArray();
    minphi = fPHOSFidCutMinPhi->GetArray();
    selection =  fPHOSFidutialCut ; 
  }
  else
    AliFatal(Form("Wrong detector name = %s", det.Data()));

  //printf("IsInFidutialCut::nphiregions = %d, netaregions = %d\n",nphiregions, netaregions);

  if(!selection) return kTRUE; //No cuts applied, all tracks/clusters used

  //Eta fidutial cut
  Bool_t bInEtaFidCut = kFALSE;
  for(Int_t ieta = 0; ieta < netaregions; ieta++)
    if(eta > mineta[ieta] && eta < maxeta[ieta]) bInEtaFidCut = kTRUE;

  if(bInEtaFidCut){
    //printf("Eta cut passed\n");
    //Phi fidutial cut
    Bool_t bInPhiFidCut = kFALSE;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      if(phi > minphi[iphi] *TMath::DegToRad()&& phi < maxphi[iphi]*TMath::DegToRad()) bInPhiFidCut = kTRUE ;
    
    if(bInPhiFidCut) {
      //printf("IsInFidutialCut:: %s cluster/track accepted\n",det.Data());
      return kTRUE;
    }
    else return kFALSE;
    
  }//In eta fid cut
  else
    return kFALSE;
  
}


//_______________________________________________________________
void AliFidutialCut::InitParameters()
{
 
  //Initialize the parameters of the analysis.

  fEMCALFidutialCut = kTRUE ;  
  fPHOSFidutialCut = kTRUE ;
  fCTSFidutialCut = kTRUE ;

  fCTSFidCutMinEta = new TArrayF(1);
  fCTSFidCutMinEta->SetAt(-0.9,0); 
  fCTSFidCutMaxEta = new TArrayF(1);
  fCTSFidCutMaxEta->SetAt(0.9,0); 

  fCTSFidCutMinPhi = new TArrayF(1);
  fCTSFidCutMinPhi->SetAt(0,0); 
  fCTSFidCutMaxPhi = new TArrayF(1);
  fCTSFidCutMaxPhi->SetAt(360.,0); 

  fEMCALFidCutMinEta = new TArrayF(1);
  fEMCALFidCutMinEta->SetAt(-0.7,0); 
  fEMCALFidCutMaxEta = new TArrayF(1);
  fEMCALFidCutMaxEta->SetAt(0.7,0); 

  fEMCALFidCutMinPhi = new TArrayF(1);
  fEMCALFidCutMinPhi->SetAt(80.,0); 
  fEMCALFidCutMaxPhi = new TArrayF(1);
  fEMCALFidCutMaxPhi->SetAt(190.,0); 

  fPHOSFidCutMinEta = new TArrayF(1);
  fPHOSFidCutMinEta->SetAt(-0.13,0); 
  fPHOSFidCutMaxEta = new TArrayF(1);
  fPHOSFidCutMaxEta->SetAt(0.13,0); 

  fPHOSFidCutMinPhi = new TArrayF(1);
  fPHOSFidCutMinPhi->SetAt(220.,0); 
  fPHOSFidCutMaxPhi = new TArrayF(1);
  fPHOSFidCutMaxPhi->SetAt(320.,0); 

}


//________________________________________________________________
void AliFidutialCut::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;

  if(fCTSFidutialCut){
    Int_t netaregions =  fCTSFidCutMaxEta->GetSize();
    Int_t nphiregions =  fCTSFidCutMaxPhi->GetSize();
    printf(">> CTS Fidutial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fCTSFidCutMinEta->GetAt(ieta), fCTSFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fCTSFidCutMinPhi->GetAt(iphi), fCTSFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fidutial cuts in CTS\n");

  if(fEMCALFidutialCut){
    Int_t netaregions =  fEMCALFidCutMaxEta->GetSize();
    Int_t nphiregions =  fEMCALFidCutMaxPhi->GetSize();
    printf(">>EMCAL Fidutial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fEMCALFidCutMinEta->GetAt(ieta), fEMCALFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fEMCALFidCutMinPhi->GetAt(iphi), fEMCALFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fidutial cuts in EMCAL\n");

  if(fPHOSFidutialCut){
    Int_t netaregions =  fPHOSFidCutMaxEta->GetSize();
    Int_t nphiregions =  fPHOSFidCutMaxPhi->GetSize();
    printf(">>PHOS Fidutial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fPHOSFidCutMinEta->GetAt(ieta), fPHOSFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fPHOSFidCutMinPhi->GetAt(iphi), fPHOSFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fidutial cuts in PHOS\n");
  printf("    \n") ;
} 


//_______________________________________________________________
void AliFidutialCut::SetSimpleCTSFidutialCut(const Float_t eta, const Float_t minphi, const Float_t maxphi){

  //Method to set simple acceptance cut to CTS
  fCTSFidCutMinEta->Set(1);
  fCTSFidCutMaxEta->Set(1);
  fCTSFidCutMinPhi->Set(1);
  fCTSFidCutMaxPhi->Set(1);

  fCTSFidCutMinEta->SetAt(-eta,0);
  fCTSFidCutMaxEta->SetAt(eta,0);
  fCTSFidCutMinPhi->SetAt(minphi,0);
  fCTSFidCutMaxPhi->SetAt(maxphi,0);


}


//_______________________________________________________________
void AliFidutialCut::SetSimpleEMCALFidutialCut(const Float_t eta, const Float_t minphi, const Float_t maxphi){
  //Method to set simple acceptance cut to EMCAL

  fEMCALFidCutMinEta->Set(1);
  fEMCALFidCutMaxEta->Set(1);
  fEMCALFidCutMinPhi->Set(1);
  fEMCALFidCutMaxPhi->Set(1);

  fEMCALFidCutMinEta->SetAt(-eta,0);
  fEMCALFidCutMaxEta->SetAt(eta,0);
  fEMCALFidCutMinPhi->SetAt(minphi,0);
  fEMCALFidCutMaxPhi->SetAt(maxphi,0);


}

//_______________________________________________________________
void AliFidutialCut::SetSimplePHOSFidutialCut(const Float_t eta, const Float_t minphi, const Float_t maxphi){

  //Method to set simple acceptance cut to PHOS
  fPHOSFidCutMinEta->Set(1);
  fPHOSFidCutMaxEta->Set(1);
  fPHOSFidCutMinPhi->Set(1);
  fPHOSFidCutMaxPhi->Set(1);

  fPHOSFidCutMinEta->SetAt(-eta,0);
  fPHOSFidCutMaxEta->SetAt(eta,0);
  fPHOSFidCutMinPhi->SetAt(minphi,0);
  fPHOSFidCutMaxPhi->SetAt(maxphi,0);


}
