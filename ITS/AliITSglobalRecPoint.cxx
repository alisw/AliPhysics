#include <stdlib.h>
#include <Riostream.h>

#include <TROOT.h>
#include <TMath.h>
#include <TString.h>

#include "AliITSglobalRecPoint.h"


ClassImp(AliITSglobalRecPoint)
//
//
//
AliITSglobalRecPoint::AliITSglobalRecPoint()
{
	// Default constructor:
	// sets to zero all members
	fLayer = 0;
	
	fGX = fGSX = 0.0;
	fGY = fGSY = 0.0;
	fGZ = fGSZ = 0.0;
	
	fR2 = fR3 = fPhi = fTheta = 0.0;
	
	fLabel[0] = fLabel[1] = fLabel[2] = -1;
	fKalmanLabel = 0;
	fSign = 0;
	fModule = 0;
	fPosInModule = 0;
	
	fUsed = 0;
}
//
//
//
AliITSglobalRecPoint::AliITSglobalRecPoint(Double_t gx, Double_t gy, 
	Double_t gz, Double_t gsx, Double_t gsy, Double_t gsz, Int_t l) : 
	fGX(gx), fGY(gy), fGZ(gz), fGSX(gsx), fGSY(gsy), fGSZ(gsz), fLayer(l)
{
	fR2 = TMath::Sqrt(fGX*fGX + fGY*fGY);
	fR3 = TMath::Sqrt(fGX*fGX + fGY*fGY + fGZ*fGZ); 
	fPhi = TMath::ATan2(fGY, fGX);
	if (fPhi < 0.0) fPhi = 2.0 * TMath::Pi() + fPhi;
	fTheta = TMath::ATan2(fR2, fGZ);
	
	fLabel[0] = fLabel[1] = fLabel[2] = -1;
	fUsed = 0;
	fKalmanLabel = 0;
	fSign = 0;
	fModule = 0;
	fPosInModule = 0;
}
//
//
/*
Double_t AliITSglobalRecPoint::GetGlobalSigmaR2()
{	
	// Square sigma of the 2-D radius:
	// = X^2 * Sx^2 + Y^2 * Sy^2
	
	Double_t answer = fGX*fGX*fGSigmaX + fGY*fGY*fGSigmaY;
	return answer / (fR2 * fR2);
}
//
//
//
Double_t AliITSglobalRecPoint::GetGlobalSigmaR3()
{
	// Square sigma of the 3-D radius:
	// = (X^2 * Sx^2 + Y^2 * Sy^2 + Z^2 * Sz^2) / R^2
	// R in 3-D, r in 2-D
	
	Double_t answer = fGX*fGX*fGSigmaX + fGY*fGY*fGSigmaY + fGZ*fGZ*fGSigmaZ;
	return answer / (fR3 * fR3);
}
//
//
//
Double_t AliITSglobalRecPoint::GetGlobalSigmaTheta()
{
	// Square sigma of theta:
	// = (Z^2 * (X^2 * Sx^2 + Y^2 * Sy^2) + r^4 * Sz^2) / (R^4 * r^2)
	// R in 3-D, r in 2-D
	
	Double_t answer = fGZ*fGZ*(fGX*fGX*fGSigmaX + fGY*fGY*fGSigmaY) + fR2*fR2*fR2*fR2*fGSigmaZ;
	return answer / (fR3*fR3*fR3*fR3*fR2*fR2);
}
//
//
//
Double_t AliITSglobalRecPoint::GetGlobalSigmaPhi()
{
	// Square sigma of phi:
	// = (Y^2 * Sx^2 + X^2 * Sy^2) / r^4
	// R in 3-D, r in 2-D
	
	Double_t answer = fGY*fGY*fGSigmaX + fGX*fGX*fGSigmaY;
	return answer / (fR2 * fR2 * fR2 * fR2);
}
//
//
//
Bool_t AliITSglobalRecPoint::SharesID(AliITSRecPoint* pt)
{
	// Controls if there is a track index shared by two points
	
	Bool_t ok = kFALSE;
	Int_t i, j;
	for (i = 0; i < 3; i++) {
		if (fTracks[i] < 0) continue;
		for (j = 0; j < 3; j++) {
			if (pt->fTracks[j] < 0) continue;
			ok = ok || (fTracks[i] == pt->fTracks[j]);
		}
	}
	return ok;
}
*/
//
//
Double_t AliITSglobalRecPoint::DPhi(AliITSglobalRecPoint *p)
{
	// Absolute value of the difference
	// between 'phi' coordinates of two points
	// controlled in order to avoid that, for
	// reasons of initial values, it come > 180. degrees
	
	Double_t phi = TMath::Abs(fPhi - p->fPhi);
	if (phi > TMath::Pi())
		phi = 2.0 * TMath::Pi() - phi;
	return phi;
}
//
//
//
Double_t AliITSglobalRecPoint::DTheta(AliITSglobalRecPoint *p)
{
	// Absolute value of the difference
	// between 'theta' coordinates of two points
	
	return TMath::Abs(fTheta - p->fTheta);
}
//
//
//
Int_t AliITSglobalRecPoint::Compare(const TObject *O) const
{
	// Criterion for sorting:
	// p1 < p2 if p1.phi < p2.phi;
	
	// Casting to convert the given argument to rec-point
	AliITSglobalRecPoint *you = (AliITSglobalRecPoint*)O;
	
	// Comparation
	if (fLayer < you->fLayer)
		return -1;
	else if (fLayer > you->fLayer)
		return 1;
	else {
		if (fPhi < you->fPhi)
			return -1;
		else if (fPhi > you->fPhi)
			return 1;
		else
			return 0;
	}
}
//
//
//
