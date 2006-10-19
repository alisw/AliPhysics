///////////////////////////////////////////////////////////////
// AliITSneuralPoint                                         //
//                                                           //
// A class which resumes the information of ITS clusters     //
// in the global reference frame.                            //
// Author: A. Pulvirenti                                     //
///////////////////////////////////////////////////////////////
//#include <stdlib.h>
//#include <Riostream.h>

#include <TString.h>

#include "AliITSRecPoint.h"
#include "AliITSRecPoint.h"
#include "AliITSgeom.h"
#include "AliITSgeomMatrix.h"

#include "AliITSNeuralPoint.h"


ClassImp(AliITSNeuralPoint)
//
//------------------------------------------------------------------------------------------------------
//
AliITSNeuralPoint::AliITSNeuralPoint():
fX(0),
fY(0),
fZ(0),
fConfX(0),
fConfY(0),
fEX(0),
fEY(0),
fEZ(0),
fCharge(0),
fModule(0),
fIndex(0),
fLayer(0),
fUser(0),
fZSort(0){
	// Default constructor.
	// Defines the point as a noise point in the origin.
	
	fLabel[0] = fLabel[1] = fLabel[2] = -1;

}
//
//------------------------------------------------------------------------------------------------------
//
AliITSNeuralPoint::AliITSNeuralPoint(AliITSNeuralPoint *p) :
fX(p->fX), 
fY(p->fY), 
fZ(p->fZ),
fConfX(p->fConfX),
fConfY(p->fConfY), 
fEX(p->fEX), 
fEY(p->fEY), 
fEZ(p->fEZ),
fCharge(p->fCharge),
fModule(p->fModule),
fIndex(p->fIndex),
fLayer(p->fLayer),
fUser(p->fUser),
fZSort(p->fZSort)
{
	// Modified copy constructor.
	// Accepts a pointer to a like object and copies its datamembers.
	

	for (Int_t i = 0; i < 3; i++) fLabel[i] = p->fLabel[i];

}
//
//------------------------------------------------------------------------------------------------------
//
AliITSNeuralPoint::AliITSNeuralPoint(AliITSRecPoint *rp, AliITSgeomMatrix *gm):
fX(0),
fY(0),
fZ(0),
fConfX(0),
fConfY(0),
fEX(0),
fEY(0),
fEZ(0),
fCharge(0),
fModule(0),
fIndex(0),
fLayer(0),
fUser(0),
fZSort(0)
{
	// Conversion constructor.
	// Accepts a AliITSRecPoint and a AliITSgeomMatrix,
	// and converts the local coord of the AliITSRecPoint object into global
	
	Int_t i, k;
	Double_t locPos[3], globPos[3], locErr[3][3], globErr[3][3];
	for (i = 0; i < 3; i++) {
		globPos[i] = 0.0;
		for (k = 0; k < 3; k++) {
			locErr[i][k] = 0.0;
			globErr[i][k] = 0.0;
		}
	}
	
	// local to global conversions of coords
	locPos[0] = rp->GetDetLocalX();
	locPos[1] = 0.0;
	locPos[2] = rp->GetDetLocalZ();
	gm->LtoGPosition(locPos, globPos);
	fX = globPos[0];
	fY = globPos[1];
	fZ = globPos[2];

	// local to global conversions of sigmas
	locErr[0][0] = rp->GetSigmaDetLocX2();
	locErr[2][2] = rp->GetSigmaZ2();
	gm->LtoGPositionError(locErr, globErr);
	for (i = 0; i < 3; i++) fLabel[i] = rp->GetLabel(i);
	fEX = TMath::Sqrt(globErr[0][0]);
	fEY = TMath::Sqrt(globErr[1][1]);
	fEZ = TMath::Sqrt(globErr[2][2]);

	// copy of other data-members
	fCharge = rp->GetQ();
	fLayer = 0;
	fIndex = 0;
	fModule = 0;
	fUser = 0;
}
//
//-------------------------------------------------------------------------------------------------
//
AliITSNeuralPoint::AliITSNeuralPoint(AliITSRecPoint *rp, AliITSgeom *geom, Short_t module, Short_t index):
fX(0),
fY(0),
fZ(0),
fConfX(0),
fConfY(0),
fEX(0),
fEY(0),
fEZ(0),
fCharge(0),
fModule(0),
fIndex(0),
fLayer(0),
fUser(0),
fZSort(0)
{
	// Conversion constructor.
	// Accepts a AliITSRecPoint and an AliITSgeom,
	// and converts the local coord of the AliITSRecPoint object into global
	
	Int_t mod = (Int_t)module, lay, lad, det;
	fModule = module;
	fIndex = index;
	geom->GetModuleId(mod, lay, lad, det);
	fLayer = (Short_t)lay;
	
	Double_t rot[9];  
	Float_t tx, ty, tz;
	geom->GetRotMatrix(fModule, rot);
	geom->GetTrans(fLayer, lad, det, tx, ty, tz);
	
	Double_t r, phi, cosPhi, sinPhi;
	r = -tx*rot[1] + ty*rot[0];
	if (lay == 1) r = -r;
	phi = TMath::ATan2(rot[1], rot[0]);
	if (lay==1) phi -= 3.1415927;
	cosPhi = TMath::Cos(phi);
	sinPhi = TMath::Sin(phi);
	fX =  r*cosPhi + rp->GetY()*sinPhi;
	fY = -r*sinPhi + rp->GetY()*cosPhi;
	fZ = rp->GetZ();
	fEX = TMath::Sqrt(rp->GetSigmaY2())*sinPhi;
	fEY = TMath::Sqrt(rp->GetSigmaY2())*cosPhi;
	fEZ = TMath::Sqrt(rp->GetSigmaZ2());
	fLayer--;
}
//
//-------------------------------------------------------------------------------------------------
//
Double_t AliITSNeuralPoint::GetPhi() const
{
	// Returns the azimuthal coordinate in the range 0-2pi
	Double_t q;
	q = TMath::ATan2(fY,fX); 
	if (q >= 0.) 
		return q;
	else 
		return q + 2.0*TMath::Pi();
}
//
//------------------------------------------------------------------------------------------------------
//
Double_t AliITSNeuralPoint::GetError(Option_t *option)
{
// Returns the error or the square error of
// values related to the coordinates in different systems.
// The option argument specifies the coordinate error desired:
//
// "R2"     --> error in transverse radius
// "R3"     --> error in spherical radius
// "PHI"    --> error in azimuthal angle
// "THETA"  --> error in polar angle
// "SQ"     --> get the square of error
//
// In order to get the error on the cartesian coordinates
// reference to the inline ErrX(), ErrY() adn ErrZ() methods.

	TString opt(option);
	Double_t errorSq = 0.0;
	opt.ToUpper();

	if (opt.Contains("R2")) {
		errorSq  = fX*fX*fEX*fEX + fY*fY*fEY*fEY;
		errorSq /= GetR2sq();
	}
	else if (opt.Contains("R3")) {
		errorSq  = fX*fX*fEX*fEX + fY*fY*fEY*fEY + fZ*fZ*fEZ*fEZ;
		errorSq /= GetR3sq();
	}
	else if (opt.Contains("PHI")) {
		errorSq  = fY*fY*fEX*fEX;
		errorSq += fX*fX*fEY*fEY;
		errorSq /= GetR2sq() * GetR2sq();
	}
	else if (opt.Contains("THETA")) {
		errorSq = fZ*fZ * (fX*fX*fEX*fEX + fY*fY*fEY*fEY);
		errorSq += GetR2sq() * GetR2sq() * fEZ*fEZ;
		errorSq /= GetR3sq() * GetR3sq() * GetR2() * GetR2();
	}
	
	if (!opt.Contains("SQ")) 
		return TMath::Sqrt(errorSq);
	else 
		return errorSq;
}
//
//------------------------------------------------------------------------------------------------------
//
Bool_t AliITSNeuralPoint::HasID(Int_t ID) const
{
// Checks if the recpoint belongs to the GEANT track
// whose label is specified in the argument
	if (ID<0) 
		return kFALSE; 
	else 
		return (fLabel[0]==ID || fLabel[1]==ID || fLabel[2]==ID);
}
//
//------------------------------------------------------------------------------------------------------
//
Int_t* AliITSNeuralPoint::SharedID(AliITSNeuralPoint *p) const
{
// Checks if there is a GEANT track owning both
// <this> and the recpoint in the argument
// The return value is an array of 4 integers.
// The firs integer returns the count of matches between labels of 
// <this> and labels of the argument (0 to 3)
// The other three return the matched labels.
// If a NULL pointer is passed, the array will be returned as:
// {0, -1, -1, -1}

	Int_t i, *shared = new Int_t[4];
	for (i = 0; i < 4; i++) shared[i] = -1;
	shared[0] = 0;
	if (!p) return shared;
	for (i = 0; i < 3; i++) {
		if (HasID(p->fLabel[i])) shared[i + 1] = p->fLabel[i];
		shared[0]++;
	}
	return shared;
}
//
//
//
void AliITSNeuralPoint::ConfMap(Double_t vx, Double_t vy)
{
// Performs conformal mapping vertex-constrained

	Double_t dx = fX - vx;
	Double_t dy = vy - fY;
	Double_t r2 = dx*dx + dy*dy;
	fConfX = dx / r2;
	fConfY = dy / r2;
}
