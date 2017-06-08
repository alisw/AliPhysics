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

// --- ROOT system ---
#include <TMath.h>
#include <TString.h>

//---- ANALYSIS system ----
#include "AliFiducialCut.h"
#include <AliLog.h>

/// \cond CLASSIMP
ClassImp(AliFiducialCut) ;
/// \endcond

//________________________________
/// Default constructor. 
/// Initialize parameters
//________________________________
AliFiducialCut::AliFiducialCut() : 
TObject(),
fEMCALFiducialCut(0),   fDCALFiducialCut(0),    fPHOSFiducialCut(0),     fCTSFiducialCut(0),
fCTSFidCutMinEta(0x0),  fCTSFidCutMinPhi(0x0),  fCTSFidCutMaxEta(0x0),   fCTSFidCutMaxPhi(0x0),
fEMCALFidCutMinEta(0x0),fEMCALFidCutMinPhi(0x0),fEMCALFidCutMaxEta(0x0), fEMCALFidCutMaxPhi(0x0),
fPHOSFidCutMinEta(0x0), fPHOSFidCutMinPhi(0x0), fPHOSFidCutMaxEta(0x0),  fPHOSFidCutMaxPhi(0x0),
fDCALFidCutMinEta(0x0), fDCALFidCutMinPhi(0x0), fDCALFidCutMaxEta(0x0),  fDCALFidCutMaxPhi(0x0)
{
  InitParameters();
}

//_______________________________
/// Destructor
//_______________________________
AliFiducialCut::~AliFiducialCut()
{  
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

//________________________________________________________________________________
/// Select EMCAL or PHOS cluster or CTS track or particle 
/// if it is inside eta-phi defined regions
/// \param eta: track/cluster/particle pseudorapidity
/// \param phi: track/cluster/particle azimuthal angle
/// \param det: detector tag where region is checked
/// \return kTRUE if cluster/track is in une of the defined regions
//________________________________________________________________________________
Bool_t AliFiducialCut::IsInFiducialCut(Float_t eta, Float_t phi, Int_t det) const
{  
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
/// Given the selection regions in Eta and Phi, check if particle is in the region 
/// defined by the TArray.
/// \param eta: track/cluster/particle pseudorapidity.
/// \param phiOrg: track/cluster/particle azimuthal angle.
/// \param phimin: array with list of minimum azimuthal angle regions.
/// \param phimax: array with list of maximum azimuthal angle regions.
/// \param etamin: array with list of minimum pseudorapidity regions.
/// \param etamax: array with list of maximum pseudorapidity regions.
/// \return kTRUE if track/cluster/particle is in one of the regions
//___________________________________________________________________________________________
Bool_t AliFiducialCut::CheckFiducialRegion(Float_t eta, Float_t phiOrg,
                                           const TArrayF* phimin, const TArrayF* phimax, 
                                           const TArrayF* etamin, const TArrayF* etamax) const 
{  
  Float_t phi = phiOrg;
  if(phi < 0) phi+=TMath::TwoPi() ;
  phi*=TMath::RadToDeg();
  
  //printf("IsInFiducialCut::Det: %s, phi = %f, eta = %f\n", det.Data(),phi*TMath::RadToDeg(), eta);
  
  Int_t netaregions = etamax->GetSize();
  Int_t nphiregions = phimax->GetSize();
  
  if ( netaregions !=  etamin->GetSize() || nphiregions !=  phimin->GetSize() )
  AliWarning(Form("Wrong number of fiducial cut regions: netamax %d != netamin %d; nphimax %d != nphimin %d\n",
                  netaregions, etamin->GetSize(),  nphiregions, phimin->GetSize()));
  
  //printf("n eta %d, nphi %d\n",netaregions,nphiregions);
  //printf("eta %2.2f, phi %2.2f\n",eta,phi);
  
  //
  // Same number of regions in eta and phi
  //
  if(netaregions == nphiregions)
  {
    for(Int_t iRegion = 0; iRegion < netaregions; iRegion++)
    {
      //      printf("region %d, min %2.2f < eta %2.2f < max %2.2f;min %2.2f < phi %2.2f < max %2.2f\n",
      //             iRegion,
      //             etamin->GetAt(iRegion),eta,etamax->GetAt(iRegion),
      //             phimax->GetAt(iRegion),phi,phimax->GetAt(iRegion));
      
      if ( eta > etamin->GetAt(iRegion) && eta < etamax->GetAt(iRegion) && 
          phi > phimin->GetAt(iRegion) && phi < phimax->GetAt(iRegion)    )
      return kTRUE ;
    }
    
    return kFALSE;
  }
  
  //
  // Different number of regions in eta and phi
  // Careful, better use same number.
  //
  
  // Eta fiducial cut
  Bool_t bInEtaFidCut = kFALSE;
  for(Int_t ieta = 0; ieta < netaregions; ieta++)
  {
    if(eta > etamin->GetAt(ieta) && eta < etamax->GetAt(ieta)) bInEtaFidCut = kTRUE;
  }
  
  //printf("Eta ok? %d\n",bInEtaFidCut);
  if(bInEtaFidCut)
  {
    //printf("Eta cut passed\n");
    //Phi fiducial cut
    Bool_t bInPhiFidCut = kFALSE;
    for(Int_t iphi = 0; iphi < nphiregions; iphi++)
    {
      if ( phi > phimin->GetAt(iphi) && phi < phimax->GetAt(iphi) ) bInPhiFidCut = kTRUE ;
    }
    //printf("Phi ok? %d\n",bInPhiFidCut);
    
    if(bInPhiFidCut) 
    {
      //printf("IsInFiducialCut:: %s cluster/track accepted\n",det.Data());
      return kTRUE;
    }
    else return kFALSE;
    
  } // In eta fid cut
  else
  {
    return kFALSE;
  }
}


//____________________________________________
/// Initialize the parameters.
//____________________________________________
void AliFiducialCut::InitParameters()
{
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
    
  // Divide DCal in 3 regions:
  // A (C?) side : -0.70<eta<-0.22, 260<phi<320
  // C (A?) side :  0.22<eta< 0.70, 260<phi<320
  // 1/3 SM      : -0.70<eta< 0.70, 320<phi<327
  
  fDCALFidCutMinEta = new TArrayF(3);
  fDCALFidCutMinEta->SetAt(-0.7 ,0);
  fDCALFidCutMinEta->SetAt( 0.22,1);
  fDCALFidCutMinEta->SetAt(-0.7 ,2);
  
  fDCALFidCutMaxEta = new TArrayF(3);
  fDCALFidCutMaxEta->SetAt(-0.22,0);
  fDCALFidCutMaxEta->SetAt( 0.7 ,1);
  fDCALFidCutMaxEta->SetAt( 0.7 ,2);
  
  fDCALFidCutMinPhi = new TArrayF(3);
  fDCALFidCutMinPhi->SetAt(260.,0);
  fDCALFidCutMinPhi->SetAt(260.,1);
  fDCALFidCutMinPhi->SetAt(320.,2);
  
  fDCALFidCutMaxPhi = new TArrayF(3);
  fDCALFidCutMaxPhi->SetAt(320.,0);
  fDCALFidCutMaxPhi->SetAt(320.,1);
  fDCALFidCutMaxPhi->SetAt(327.,2);
}


//________________________________________________________________
/// Print some relevant parameters set.
//________________________________________________________________
void AliFiducialCut::Print(const Option_t * opt) const
{  
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
/// Define simple acceptance cut to tracks.
/// \param eta: absolute maximum value of track pseudorapidity.
/// \param phimin track minimum azimuthal angle.
/// \param phimax track maximum azimuthal angle.
//_______________________________________________________________________________________
void AliFiducialCut::SetSimpleCTSFiducialCut(Float_t eta, Float_t phimin, Float_t phimax)
{
  fCTSFidCutMinEta->Set(1);
  fCTSFidCutMaxEta->Set(1);
  fCTSFidCutMinPhi->Set(1);
  fCTSFidCutMaxPhi->Set(1);
  
  fCTSFidCutMinEta->SetAt(-eta,0);
  fCTSFidCutMaxEta->SetAt( eta,0);
  fCTSFidCutMinPhi->SetAt(phimin,0);
  fCTSFidCutMaxPhi->SetAt(phimax,0);
}

//_________________________________________________________________________________________
/// Define simple acceptance cut to EMCal clusters.
/// \param eta: absolute maximum value of cluster pseudorapidity.
/// \param phimin cluster minimum azimuthal angle.
/// \param phimax cluster maximum azimuthal angle.
//_________________________________________________________________________________________
void AliFiducialCut::SetSimpleEMCALFiducialCut(Float_t eta, Float_t phimin, Float_t phimax)
{  
  fEMCALFidCutMinEta->Set(1);
  fEMCALFidCutMaxEta->Set(1);
  fEMCALFidCutMinPhi->Set(1);
  fEMCALFidCutMaxPhi->Set(1);
  
  fEMCALFidCutMinEta->SetAt(-eta,0);
  fEMCALFidCutMaxEta->SetAt( eta,0);
  fEMCALFidCutMinPhi->SetAt(phimin,0);
  fEMCALFidCutMaxPhi->SetAt(phimax,0);
}

//_______________________________________________________________________________________
/// Define simple acceptance cut to PHOS clusters.
/// \param eta: absolute maximum value of cluster pseudorapidity.
/// \param phimin cluster minimum azimuthal angle.
/// \param phimax cluster maximum azimuthal angle.
//_______________________________________________________________________________________
void AliFiducialCut::SetSimplePHOSFiducialCut(Float_t eta, Float_t phimin, Float_t phimax)
{  
  fPHOSFidCutMinEta->Set(1);
  fPHOSFidCutMaxEta->Set(1);
  fPHOSFidCutMinPhi->Set(1);
  fPHOSFidCutMaxPhi->Set(1);
  
  fPHOSFidCutMinEta->SetAt(-eta,0);
  fPHOSFidCutMaxEta->SetAt(eta,0);
  fPHOSFidCutMinPhi->SetAt(phimin,0);
  fPHOSFidCutMaxPhi->SetAt(phimax,0);
}

//_________________________________________________________________________________________
/// Define simple acceptance cut to DCal clusters.
/// \param eta: absolute maximum value of cluster pseudorapidity.
/// \param phimin cluster minimum azimuthal angle.
/// \param phimax cluster maximum azimuthal angle.
//_________________________________________________________________________________________
void AliFiducialCut::SetSimpleDCALFiducialCut(Float_t eta, Float_t phimin, Float_t phimax)
{  
  fDCALFidCutMinEta->Set(1);
  fDCALFidCutMaxEta->Set(1);
  fDCALFidCutMinPhi->Set(1);
  fDCALFidCutMaxPhi->Set(1);
  
  fDCALFidCutMinEta->SetAt(-eta,0);
  fDCALFidCutMaxEta->SetAt( eta,0);
  fDCALFidCutMinPhi->SetAt(phimin,0);
  fDCALFidCutMaxPhi->SetAt(phimax,0);
}

//_________________________________________________________________________________________
/// Define acceptance cut to DCal clusters 
/// with more accurate description than in SetSimpleDCALFiducialCut()
/// \param etaminFull  minimum value of cluster pseudorapidity in Full SMs.
/// \param etamaxFull  maximum value of cluster pseudorapidity in Full SMs.
/// \param phiminFull  cluster minimum azimuthal angle in Full SMs.
/// \param phimaxFull  cluster maximum azimuthal angle in Full SMs.
/// \param etaminThird minimum value of cluster pseudorapidity in 1/3 SMs.
/// \param etamaxThird maximum value of cluster pseudorapidity in 1/3 SMs.
/// \param phiminThird cluster minimum azimuthal angle in 1/3 SMs.
/// \param phimaxThird cluster maximum azimuthal angle in 1/3 SMs.
//_________________________________________________________________________________________
void AliFiducialCut::SetDCALFiducialCut(Float_t etaminFull , Float_t etamaxFull,
                                        Float_t phiminFull , Float_t phimaxFull,
                                        Float_t etaminThird, Float_t etamaxThird,
                                        Float_t phiminThird, Float_t phimaxThird)
{  
  fDCALFidCutMinEta->Set(3);
  fDCALFidCutMaxEta->Set(3);
  fDCALFidCutMinPhi->Set(3);
  fDCALFidCutMaxPhi->Set(3);
  
  // 1/3 SM
  fDCALFidCutMinEta->SetAt(etaminThird,2);
  fDCALFidCutMaxEta->SetAt(etamaxThird,2);
  fDCALFidCutMinPhi->SetAt(phiminThird,2);
  fDCALFidCutMaxPhi->SetAt(phimaxThird,2);
  
  // Full SM positive eta
  fDCALFidCutMinEta->SetAt(etaminFull,1);
  fDCALFidCutMaxEta->SetAt(etamaxFull,1);
  fDCALFidCutMinPhi->SetAt(phiminFull,1);
  fDCALFidCutMaxPhi->SetAt(phimaxFull,1); 

  // Full SM negative eta
  fDCALFidCutMinEta->SetAt(-etamaxFull,0);
  fDCALFidCutMaxEta->SetAt(-etaminFull,0);
  fDCALFidCutMinPhi->SetAt( phiminFull,0);
  fDCALFidCutMaxPhi->SetAt( phimaxFull,0); 
}


