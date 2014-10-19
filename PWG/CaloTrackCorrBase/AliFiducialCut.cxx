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
//#include <TLorentzVector.h>
#include <TString.h>

//---- ANALYSIS system ----
#include "AliFiducialCut.h"
#include <AliLog.h>

ClassImp(AliFiducialCut)


//________________________________
AliFiducialCut::AliFiducialCut() : 
TObject(),
fEMCALFiducialCut(0),   fPHOSFiducialCut(0),    fCTSFiducialCut(0),
fCTSFidCutMinEta(0x0),  fCTSFidCutMinPhi(0x0),  fCTSFidCutMaxEta(0x0),   fCTSFidCutMaxPhi(0x0),
fEMCALFidCutMinEta(0x0),fEMCALFidCutMinPhi(0x0),fEMCALFidCutMaxEta(0x0), fEMCALFidCutMaxPhi(0x0),
fPHOSFidCutMinEta(0x0), fPHOSFidCutMinPhi(0x0), fPHOSFidCutMaxEta(0x0),  fPHOSFidCutMaxPhi(0x0),
fDCALFidCutMinEta(0x0), fDCALFidCutMinPhi(0x0), fDCALFidCutMaxEta(0x0),  fDCALFidCutMaxPhi(0x0)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
  
}

//_______________________________
AliFiducialCut::~AliFiducialCut()
{
  //Dtor
  
  if(fCTSFidCutMinEta)   delete fCTSFidCutMinEta ;
  if(fCTSFidCutMinPhi)   delete fCTSFidCutMinPhi ;
  if(fCTSFidCutMaxEta)   delete fCTSFidCutMaxEta ;
  if(fCTSFidCutMaxPhi)   delete fCTSFidCutMaxPhi ;
  
  if(fEMCALFidCutMinEta) delete fEMCALFidCutMinEta ;
  if(fEMCALFidCutMinPhi) delete fEMCALFidCutMinPhi ;
  if(fEMCALFidCutMaxEta) delete fEMCALFidCutMaxEta ; 
  if(fEMCALFidCutMaxPhi) delete fEMCALFidCutMaxPhi ;
  
  if(fPHOSFidCutMinEta)  delete fPHOSFidCutMinEta ; 
  if(fPHOSFidCutMinPhi)  delete fPHOSFidCutMinPhi ; 
  if(fPHOSFidCutMaxEta)  delete fPHOSFidCutMaxEta ;
  if(fPHOSFidCutMaxPhi)  delete fPHOSFidCutMaxPhi ;
  
  if(fDCALFidCutMinEta)  delete fDCALFidCutMinEta ;
  if(fDCALFidCutMinPhi)  delete fDCALFidCutMinPhi ;
  if(fDCALFidCutMaxEta)  delete fDCALFidCutMaxEta ;
  if(fDCALFidCutMaxPhi)  delete fDCALFidCutMaxPhi ;
}


////________________________________________________________________________________
//Bool_t AliFiducialCut::IsInFiducialCut(TLorentzVector momentum, TString det) const
//{
//  // Selects EMCAL or PHOS cluster or CTS track if it is inside eta-phi defined regions
//  Int_t idet = -1;
//  if     (det=="EMCAL") idet = kEMCAL;
//  else if(det=="PHOS" ) idet = kPHOS;
//  else if(det=="CTS")   idet = kCTS;
//  else if(det=="DCAL")  idet = kDCAL;
//  else if(det.Contains("DCAL") && det.Contains("PHOS")) idet = kDCALPHOS;
//  else
//  {
//    AliFatal(Form("Detector < %s > not known!", det.Data()));
//    return kFALSE;
//  }
//  
//  return IsInFiducialCut(momentum.Eta(), momentum.Phi(), idet);
//}

//________________________________________________________________________________
Bool_t AliFiducialCut::IsInFiducialCut(Float_t eta, Float_t phi, Int_t det) const
{
  // Selects EMCAL or PHOS cluster or CTS track if it is inside eta-phi defined regions
  
  if(det == kCTS)
  {
	  if(!fCTSFiducialCut)  
      return kTRUE; //Fiducial cut not requested, accept all tracks  
	  else 
      return CheckFiducialRegion(eta,phi, fCTSFidCutMinPhi  , fCTSFidCutMaxPhi , fCTSFidCutMinEta  , fCTSFidCutMaxEta  );
  }
  else   if(det == kEMCAL)
  {
	  if(!fEMCALFiducialCut) 
      return kTRUE; //Fiducial cut not requested, accept all clusters  
	  else                   
      return CheckFiducialRegion(eta,phi, fEMCALFidCutMinPhi, fEMCALFidCutMaxPhi, fEMCALFidCutMinEta, fEMCALFidCutMaxEta);
  }
  else   if(det == kPHOS)
  {
    if(!fPHOSFiducialCut)
      return kTRUE; //Fiducial cut not requested, accept all clusters
    else
      return CheckFiducialRegion(eta,phi, fPHOSFidCutMinPhi , fPHOSFidCutMaxPhi , fPHOSFidCutMinEta , fPHOSFidCutMaxEta );
  }
  else   if(det == kDCAL || det == kDCALPHOS)
  {
    if(!fDCALFiducialCut)
      return kTRUE; //Fiducial cut not requested, accept all clusters
    else
      return CheckFiducialRegion(eta,phi, fDCALFidCutMinPhi , fDCALFidCutMaxPhi , fDCALFidCutMinEta , fDCALFidCutMaxEta );
  }
  else
  {
    return kFALSE;
    AliFatal(Form("Detector < %d > not known!", det));
  }
  
}

//___________________________________________________________________________________________
Bool_t AliFiducialCut::CheckFiducialRegion(Float_t eta, Float_t phiOrg,
                                           const TArrayF* minphi, const TArrayF* maxphi, 
                                           const TArrayF* mineta, const TArrayF* maxeta) const 
{
  //Given the selection regions in Eta and Phi, check if particle is in this region.
  
  Float_t phi = phiOrg;
	if(phi < 0) phi+=TMath::TwoPi() ;
  
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
  fPHOSFiducialCut  = kTRUE ;
  fCTSFiducialCut   = kTRUE ;
  fDCALFiducialCut  = kTRUE ;
  
  fCTSFidCutMinEta = new TArrayF(1);
  fCTSFidCutMinEta->SetAt(-0.9,0); 
  fCTSFidCutMaxEta = new TArrayF(1);
  fCTSFidCutMaxEta->SetAt( 0.9,0); 
  
  fCTSFidCutMinPhi = new TArrayF(1);
  fCTSFidCutMinPhi->SetAt(0.  ,0); 
  fCTSFidCutMaxPhi = new TArrayF(1);
  fCTSFidCutMaxPhi->SetAt(360.,0); 
  
  fEMCALFidCutMinEta = new TArrayF(1);
  fEMCALFidCutMinEta->SetAt(-0.7,0); 
  fEMCALFidCutMaxEta = new TArrayF(1);
  fEMCALFidCutMaxEta->SetAt( 0.7,0); 
  
  fEMCALFidCutMinPhi = new TArrayF(1);
  fEMCALFidCutMinPhi->SetAt(80.,0); 
  fEMCALFidCutMaxPhi = new TArrayF(1);
  fEMCALFidCutMaxPhi->SetAt(187.,0); 
  
  fPHOSFidCutMinEta = new TArrayF(1);
  fPHOSFidCutMinEta->SetAt(-0.13,0); 
  fPHOSFidCutMaxEta = new TArrayF(1);
  fPHOSFidCutMaxEta->SetAt( 0.13,0); 
  
  fPHOSFidCutMinPhi = new TArrayF(1);
  fPHOSFidCutMinPhi->SetAt(260.,0); 
  fPHOSFidCutMaxPhi = new TArrayF(1);
  fPHOSFidCutMaxPhi->SetAt(320.,0); 
  
//  fDCALFidCutMinEta = new TArrayF(1);
//  fDCALFidCutMinEta->SetAt(-0.7,0);
//  fDCALFidCutMaxEta = new TArrayF(1);
//  fDCALFidCutMaxEta->SetAt( 0.7,0);
//  
//  fDCALFidCutMinPhi = new TArrayF(1);
//  fDCALFidCutMinPhi->SetAt(260.,0);
//  fDCALFidCutMaxPhi = new TArrayF(1);
//  fDCALFidCutMaxPhi->SetAt(327.,0);
  
  // Divide DCal in 3 regions:
  // A (C?) side : -0.70<eta<-0.15, 260<phi<320
  // C (A?) side :  0.15<eta< 0.70, 260<phi<320
  // 1/3 SM      : -0.70<eta< 0.70, 320<phi<327
  
  fDCALFidCutMinEta = new TArrayF(3);
  fDCALFidCutMinEta->SetAt(-0.7 ,0);
  fDCALFidCutMinEta->SetAt( 0.15,1);
  fDCALFidCutMinEta->SetAt(-0.7 ,2);
  fDCALFidCutMaxEta = new TArrayF(3);
  fDCALFidCutMaxEta->SetAt(-0.15,0);
  fDCALFidCutMaxEta->SetAt( 0.7 ,1);
  fDCALFidCutMaxEta->SetAt( 0.7 ,2);
  
  fDCALFidCutMinPhi = new TArrayF(3);
  fDCALFidCutMinPhi->SetAt(260.,0);
  fDCALFidCutMinPhi->SetAt(260.,1);
  fDCALFidCutMinPhi->SetAt(320.,2);
  fDCALFidCutMaxPhi = new TArrayF(3);
  fDCALFidCutMaxPhi->SetAt(320.,0);
  fDCALFidCutMaxPhi->SetAt(320.,0);
  fDCALFidCutMaxPhi->SetAt(327.,0);

}


//________________________________________________________________
void AliFiducialCut::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  
  if(fCTSFiducialCut)
  {
    Int_t netaregions =  fCTSFidCutMaxEta->GetSize();
    Int_t nphiregions =  fCTSFidCutMaxPhi->GetSize();
    printf(">> CTS Fiducial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fCTSFidCutMinEta->GetAt(ieta), fCTSFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fCTSFidCutMinPhi->GetAt(iphi), fCTSFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fiducial cuts in CTS\n");
  
  if(fEMCALFiducialCut)
  {
    Int_t netaregions =  fEMCALFidCutMaxEta->GetSize();
    Int_t nphiregions =  fEMCALFidCutMaxPhi->GetSize();
    printf(">>EMCAL Fiducial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fEMCALFidCutMinEta->GetAt(ieta), fEMCALFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fEMCALFidCutMinPhi->GetAt(iphi), fEMCALFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fiducial cuts in EMCAL\n");
  
  if(fPHOSFiducialCut)
  {
    Int_t netaregions =  fPHOSFidCutMaxEta->GetSize();
    Int_t nphiregions =  fPHOSFidCutMaxPhi->GetSize();
    printf(">>PHOS Fiducial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fPHOSFidCutMinEta->GetAt(ieta), fPHOSFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fPHOSFidCutMinPhi->GetAt(iphi), fPHOSFidCutMaxPhi->GetAt(iphi)) ; 
  }
  else printf(">>No fiducial cuts in PHOS\n");
  
  if(fDCALFiducialCut)
  {
    Int_t netaregions =  fDCALFidCutMaxEta->GetSize();
    Int_t nphiregions =  fDCALFidCutMaxPhi->GetSize();
    printf(">>DCAL Fiducial regions : phi %d eta %d\n", netaregions, nphiregions) ;
    for(Int_t ieta = 0; ieta < netaregions; ieta++)
      printf(" region %d : %3.2f < eta < %3.2f\n", ieta, fDCALFidCutMinEta->GetAt(ieta), fDCALFidCutMaxEta->GetAt(ieta)) ;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
      printf(" region %d : %3.1f < phi < %3.1f\n", iphi, fDCALFidCutMinPhi->GetAt(iphi), fDCALFidCutMaxPhi->GetAt(iphi)) ;
  }
  else printf(">>No fiducial cuts in DCAL\n");
  
  printf("    \n") ;
  
} 

//_______________________________________________________________________________________
void AliFiducialCut::SetSimpleCTSFiducialCut(Float_t eta, Float_t minphi, Float_t maxphi)
{
  
  //Method to set simple acceptance cut to CTS
  
  fCTSFidCutMinEta->Set(1);
  fCTSFidCutMaxEta->Set(1);
  fCTSFidCutMinPhi->Set(1);
  fCTSFidCutMaxPhi->Set(1);
  
  fCTSFidCutMinEta->SetAt(-eta,0);
  fCTSFidCutMaxEta->SetAt( eta,0);
  fCTSFidCutMinPhi->SetAt(minphi,0);
  fCTSFidCutMaxPhi->SetAt(maxphi,0);
  
}

//_________________________________________________________________________________________
void AliFiducialCut::SetSimpleEMCALFiducialCut(Float_t eta, Float_t minphi, Float_t maxphi)
{
  //Method to set simple acceptance cut to EMCAL
  
  fEMCALFidCutMinEta->Set(1);
  fEMCALFidCutMaxEta->Set(1);
  fEMCALFidCutMinPhi->Set(1);
  fEMCALFidCutMaxPhi->Set(1);
  
  fEMCALFidCutMinEta->SetAt(-eta,0);
  fEMCALFidCutMaxEta->SetAt( eta,0);
  fEMCALFidCutMinPhi->SetAt(minphi,0);
  fEMCALFidCutMaxPhi->SetAt(maxphi,0);
  
}

//________________________________________________________________________________________
void AliFiducialCut::SetSimplePHOSFiducialCut(Float_t eta, Float_t minphi, Float_t maxphi)
{
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

//_________________________________________________________________________________________
void AliFiducialCut::SetSimpleDCALFiducialCut(Float_t eta, Float_t minphi, Float_t maxphi)
{
  //Method to set simple acceptance cut to DCAL
  
  fDCALFidCutMinEta->Set(1);
  fDCALFidCutMaxEta->Set(1);
  fDCALFidCutMinPhi->Set(1);
  fDCALFidCutMaxPhi->Set(1);
  
  fDCALFidCutMinEta->SetAt(-eta,0);
  fDCALFidCutMaxEta->SetAt( eta,0);
  fDCALFidCutMinPhi->SetAt(minphi,0);
  fDCALFidCutMaxPhi->SetAt(maxphi,0);
  
}