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

/*
$Log$
Revision 1.6  2003/11/13 14:21:57  morsch
Coding Rule violation corrections.

Revision 1.5  2003/08/13 17:37:29  hristov
Bug fix (Alpha)

Revision 1.4  2003/08/05 16:14:20  morsch
Some problems with too big fluctuations corrected. (A. de Falco)

Revision 1.1  2003/01/06 10:13:09  morsch
First commit.

*/

// Implementation of AliFastResponse for the Muon Spectrometer resolution.
// The response depends on the charge of the muon and
// the background level.
// The class uses the instance of an object of type AliMUONFastTracking to 
// obtain the smearing parameters.
// Author: andreas.morsch@cern.ch

#include "AliFastMuonTrackingRes.h"
#include "AliMUONFastTracking.h"
#include <TRandom.h>
#include <TF1.h>

ClassImp(AliFastMuonTrackingRes)


AliFastMuonTrackingRes::AliFastMuonTrackingRes() :
    AliFastResponse("Resolution", "Muon Tracking Resolution")
{
// Deafault constructor
    SetBackground();
}

void AliFastMuonTrackingRes::Init()
{
// Initialisation
    fFastTracking = AliMUONFastTracking::Instance();
    fFastTracking->Init(fBackground);
}



void AliFastMuonTrackingRes::Evaluate(Float_t   p,  Float_t  theta , Float_t   phi,
					 Float_t& pS,  Float_t& thetaS, Float_t&  phiS)
{
//
// Evaluate Gaussian smearing from given kinematics 
//

    Double_t meanp    = fFastTracking->MeanP  (p, theta, phi, Int_t(fCharge));
    Double_t sigmap   = fFastTracking->SigmaP (p, theta, phi, Int_t(fCharge));
    Double_t sigma1p  = fFastTracking->Sigma1P(p, theta, phi, Int_t(fCharge));
    Double_t normg2   = fFastTracking->NormG2 (p, theta, phi, Int_t(fCharge));
    Double_t meang2   = fFastTracking->MeanG2 (p, theta, phi, Int_t(fCharge));
    Double_t sigmag2  = fFastTracking->SigmaG2(p, theta, phi, Int_t(fCharge));
    
    Int_t ip,itheta,iphi;
    fFastTracking->GetIpIthetaIphi(p, theta, phi, Int_t(fCharge), ip, itheta, iphi);
    TF1* fitp = fFastTracking->GetFitP(ip,itheta,iphi);
    
    Float_t curmeanp   = fitp->GetParameter(0); 
    Float_t cursigmap  = fitp->GetParameter(1); 
    Float_t cursigma1p = fitp->GetParameter(2); 
    Float_t curnormg2  = fitp->GetParameter(3); 
    Float_t curmeang2  = fitp->GetParameter(4); 
    Float_t cursigmag2 = fitp->GetParameter(5); 
    if (curmeanp != meanp || cursigmap != sigmap || cursigma1p != sigma1p || 
	curnormg2 != normg2 || curmeang2 != meang2 || cursigmag2 != sigmag2){ 
      printf ("Setting new parameters for ip=%d itheta=%d iphi=%d\n",ip,itheta,iphi); 
      fitp->SetParameters(meanp,sigmap,sigma1p,normg2,meang2,sigmag2);
    }
  
    Double_t meantheta  = fFastTracking->MeanTheta (p, theta, phi, Int_t(fCharge));
    Double_t sigmatheta = fFastTracking->SigmaTheta(p, theta, phi, Int_t(fCharge));
    Double_t meanphi    = fFastTracking->MeanPhi   (p, theta, phi, Int_t(fCharge));
    Double_t sigmaphi   = fFastTracking->SigmaPhi  (p, theta, phi, Int_t(fCharge));
    
    if (sigmatheta<0 || sigmaphi<0) 
      printf ("bin %d %d %d sigmatheta = %f, sigmaphi = %f\n",
	      ip,itheta,iphi,sigmatheta,sigmaphi);
    // Components different from ip=0 have the RMS bigger than mean
    Float_t ptp[3]  =  { 1.219576,-0.354764,-0.690117 };
    Float_t ptph[3] =  { 0.977522, 0.016269, 0.023158 }; 
    Float_t pphp[3] =  { 1.303256,-0.464847,-0.869322 };

    // Smeared momentum
    pS = -1.;
    //    Float_t dpmax = 5. + ip * 2.5; 
    //    Float_t dpmax = 5. + ip * 2; 
    Float_t dpmax;
    if (sigmag2<999.) dpmax = 5. * TMath::Abs(sigmap + sigmag2); 
    else dpmax = 5. * TMath::Abs(sigmap);
    Float_t dp = 100;
    while (pS<0 || TMath::Abs(dp)>dpmax) { 
      pS = p + fitp->GetRandom();
      dp = pS - p; 
    }
    // Smeared phi
    Float_t sigmaphiold=sigmaphi;
    if (ip==0) sigmaphi *= pphp[0] + pphp[1] * dp + pphp[2] * dp*dp; 
    if (sigmaphi<0.5 * sigmaphiold) sigmaphi = 0.5 * sigmaphiold;
    if (sigmaphi>2.  * sigmaphiold) sigmaphi = 2.  * sigmaphiold;
    phiS = phi + gRandom->Gaus(meanphi, sigmaphi);
    Float_t dphi = phiS - phi; 
    // Smeared theta
    Float_t sigmathetaold=sigmatheta;
    if (ip==0) sigmatheta *= ptp[0] + ptp[1] * dp + ptp[2] * dp*dp; 
    if (ip==0) sigmatheta *= ptph[0] + ptph[1] * dphi + ptph[2] * dphi*dphi; 
    if (sigmatheta<0.5 * sigmathetaold) sigmatheta = 0.5 * sigmathetaold;
    if (sigmatheta>2.  * sigmathetaold) sigmatheta = 2.  * sigmathetaold;
    thetaS = theta + gRandom->Gaus(meantheta,sigmatheta);
}






