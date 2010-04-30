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
/* $Id: AliFiducialCut.cxx 21839 2007-10-29 13:49:42Z gustavo $ */

//_________________________________________________________________________
// Class for track/cluster acceptance selection
// Selection in Central barrel, EMCAL and PHOS
//  
// Several selection regions possible for the different
// detectors
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TMath.h>
#include <TLorentzVector.h>
#include <TString.h>
//#include <TArrayF.h>

//---- ANALYSIS system ----
#include "AliFiducialCut.h"

ClassImp(AliFiducialCut)


//____________________________________________________________________________
AliFiducialCut::AliFiducialCut() : 
  TObject(),
  fEMCALFiducialCut(0),  fPHOSFiducialCut(0),  fCTSFiducialCut(0),
  fCTSFidCutMinEta(0x0),fCTSFidCutMinPhi(0x0),fCTSFidCutMaxEta(0x0), fCTSFidCutMaxPhi(0x0),
  fEMCALFidCutMinEta(0x0),fEMCALFidCutMinPhi(0x0),fEMCALFidCutMaxEta(0x0), fEMCALFidCutMaxPhi(0x0),
  fPHOSFidCutMinEta(0x0),fPHOSFidCutMinPhi(0x0),fPHOSFidCutMaxEta(0x0), fPHOSFidCutMaxPhi(0x0)

{
  //Ctor

  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliFiducialCut::AliFiducialCut(const AliFiducialCut & g) :   
  TObject(g), 
  fEMCALFiducialCut(g.fEMCALFiducialCut),  fPHOSFiducialCut(g.fPHOSFiducialCut),  fCTSFiducialCut(g. fCTSFiducialCut),
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
AliFiducialCut & AliFiducialCut::operator = (const AliFiducialCut & source)
{
  // assignment operator
  
  if(&source == this) return *this;
  
  fEMCALFiducialCut = source.fEMCALFiducialCut;  
  fPHOSFiducialCut = source.fPHOSFiducialCut;
  fCTSFiducialCut = source.fCTSFiducialCut;
  
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
AliFiducialCut::~AliFiducialCut() {
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
Bool_t AliFiducialCut::IsInFiducialCut(const TLorentzVector momentum, const TString det) const
{
  //Selects EMCAL or PHOS cluster or CTS track if it is inside eta-phi defined regions

  if(det == "CTS"){
	  if(!fCTSFiducialCut)  return kTRUE; //Fiducial cut not requested, accept all tracks  
	  else return CheckFiducialRegion(momentum, fCTSFidCutMinPhi  , fCTSFidCutMaxPhi , fCTSFidCutMinEta  , fCTSFidCutMaxEta  );
  }
  else   if(det == "EMCAL") {
	  if(!fEMCALFiducialCut) return kTRUE; //Fiducial cut not requested, accept all clusters  
	  else return CheckFiducialRegion(momentum, fEMCALFidCutMinPhi, fEMCALFidCutMaxPhi, fEMCALFidCutMinEta, fEMCALFidCutMaxEta);
  }
  else   if(det == "PHOS")  {
	  if(!fPHOSFiducialCut) return kTRUE; //Fiducial cut not requested, accept all clusters 
	  else return CheckFiducialRegion(momentum, fPHOSFidCutMinPhi , fPHOSFidCutMaxPhi , fPHOSFidCutMinEta , fPHOSFidCutMaxEta );
  }
  else{
	  printf("AliFiducialCut::IsInFiducialCut() - Wrong detector name = %s\n", det.Data());
	  return kFALSE;
  }

}

//_______________________________________________________________
Bool_t AliFiducialCut::CheckFiducialRegion(const TLorentzVector momentum, const TArrayF* minphi, const TArrayF* maxphi, const TArrayF* mineta, const TArrayF* maxeta) const {
  //Given the selection regions in Eta and Phi, check if particle is in this region.

    Double_t phi = momentum.Phi();
	if(phi < 0) phi+=TMath::TwoPi() ;
	Double_t eta =  momentum.Eta();
	//printf("IsInFiducialCut::Det: %s, phi = %f, eta = %f\n", det.Data(),phi*TMath::RadToDeg(), eta);

    Int_t netaregions = maxeta->GetSize();
    Int_t nphiregions = maxphi->GetSize();
    if(netaregions !=  mineta->GetSize() || nphiregions !=  minphi->GetSize())
		printf("AliFiducialCut::IsInFiducialCut() - Wrong number of fiducial cut regions: nmaxeta %d != nmineta %d; nmaxphi %d != nminphi %d\n",
			   netaregions, mineta->GetSize(),  nphiregions, minphi->GetSize());
	
	//Eta fiducial cut
	Bool_t bInEtaFidCut = kFALSE;
	for(Int_t ieta = 0; ieta < netaregions; ieta++)
		if(eta > mineta->GetAt(ieta) && eta < maxeta->GetAt(ieta)) bInEtaFidCut = kTRUE;

	if(bInEtaFidCut){
		//printf("Eta cut passed\n");
		//Phi fiducial cut
		Bool_t bInPhiFidCut = kFALSE;
		for(Int_t iphi = 0; iphi < nphiregions; iphi++)
			if(phi > minphi->GetAt(iphi) *TMath::DegToRad()&& phi < maxphi->GetAt(iphi)*TMath::DegToRad()) bInPhiFidCut = kTRUE ;
	  
		if(bInPhiFidCut) {
			//printf("IsInFiducialCut:: %s cluster/track accepted\n",det.Data());
			return kTRUE;
		}
		else return kFALSE;
    
	}//In eta fid cut
	else
    return kFALSE;
}


//_______________________________________________________________
void AliFiducialCut::InitParameters()
{
 
  //Initialize the parameters of the analysis.

  fEMCALFiducialCut = kTRUE ;  
  fPHOSFiducialCut = kTRUE ;
  fCTSFiducialCut = kTRUE ;

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
void AliFiducialCut::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;

  if(fCTSFiducialCut){
    Int_t netaregions =  fCTSFidCutMaxEta->GetSize();
    Int_t nphiregions =  fCTSFidCutMaxPhi->GetSize();
    printf(">> CTS Fiducial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fCTSFidCutMinEta->GetAt(ieta), fCTSFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fCTSFidCutMinPhi->GetAt(iphi), fCTSFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fiducial cuts in CTS\n");

  if(fEMCALFiducialCut){
    Int_t netaregions =  fEMCALFidCutMaxEta->GetSize();
    Int_t nphiregions =  fEMCALFidCutMaxPhi->GetSize();
    printf(">>EMCAL Fiducial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fEMCALFidCutMinEta->GetAt(ieta), fEMCALFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fEMCALFidCutMinPhi->GetAt(iphi), fEMCALFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fiducial cuts in EMCAL\n");

  if(fPHOSFiducialCut){
    Int_t netaregions =  fPHOSFidCutMaxEta->GetSize();
    Int_t nphiregions =  fPHOSFidCutMaxPhi->GetSize();
    printf(">>PHOS Fiducial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fPHOSFidCutMinEta->GetAt(ieta), fPHOSFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fPHOSFidCutMinPhi->GetAt(iphi), fPHOSFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fiducial cuts in PHOS\n");
  printf("    \n") ;
} 


//_______________________________________________________________
void AliFiducialCut::SetSimpleCTSFiducialCut(const Float_t eta, const Float_t minphi, const Float_t maxphi){

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
void AliFiducialCut::SetSimpleEMCALFiducialCut(const Float_t eta, const Float_t minphi, const Float_t maxphi){
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
void AliFiducialCut::SetSimplePHOSFiducialCut(const Float_t eta, const Float_t minphi, const Float_t maxphi){

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
