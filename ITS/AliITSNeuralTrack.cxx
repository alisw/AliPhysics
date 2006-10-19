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

/* $Id$ */

// AliITSNeuralTrack
//
// The format of output data from Neural Tracker
// It can export data in the format of AliITSIOTrack
// (V1) tracking.
// Compatibility adaptation to V2 tracking is on the way.
// Author: A. Pulvirenti

#include <Riostream.h>
//#include <cstdlib>
//#include <cstring>

//#include <TObject.h>
//#include <TROOT.h>
#include <TMath.h>
#include <TString.h>
//#include <TObjArray.h>
//#include <TH1.h>
#include <TMatrixD.h>
#if ROOT_VERSION_CODE >= 262146
#include <TMatrixDEigen.h>
#endif

//#include "AliITSVertex.h"
#include "AliITSIOTrack.h"
#include "AliITSNeuralPoint.h"

#include "AliITSNeuralTrack.h"



ClassImp(AliITSNeuralTrack)
//
//
//
AliITSNeuralTrack::AliITSNeuralTrack() :  
fXC(0.0),
fYC(0.0),
fR(0.0),
fC(0.0),
fTanL(0.0),
fG0(0.0),
fDt(0.0),
fDz(0.0),
fStateR(0.0),
fStatePhi(0.0),
fStateZ(0.0),
fMatrix(5,5),
fChi2(0.0),
fNSteps(0.0),
fMass(0.1396),// default assumption: pion
fField(2.0),// default assumption: B = 0.4 Tesla
fPDG(0),
fLabel(0),
fCount(0),
fVertex(){
	// Default constructor
	
	Int_t i;

	for (i = 0; i < 6; i++) fPoint[i] = 0;

	fVertex.X() = 0.0;
	fVertex.Y() = 0.0;
 	fVertex.Z() = 0.0;
	fVertex.ErrX() = 0.0;
	fVertex.ErrY() = 0.0;
	fVertex.ErrZ() = 0.0;
}
//
//
//
AliITSNeuralTrack::AliITSNeuralTrack(const AliITSNeuralTrack &track) 
: TObject((TObject&)track), 
fXC(0.0),
fYC(0.0),
fR(0.0),
fC(0.0),
fTanL(0.0),
fG0(0.0),
fDt(0.0),
fDz(0.0),
fStateR(0.0),
fStatePhi(0.0),
fStateZ(0.0),
fMatrix(5,5),
fChi2(0.0),
fNSteps(0.0),
fMass(0.1396),// default assumption: pion
fField(2.0),// default assumption: B = 0.4 Tesla
fPDG(0),
fLabel(0),
fCount(0),
fVertex(){
// copy constructor

	Int_t i;

	fMass  = 0.1396;   // default assumption: pion
	fField = 2.0;      // default assumption: B = 0.4 Tesla

	fXC = fYC = fR = fC = 0.0;
	fTanL = fG0 = fDt = fDz = 0.0;
	fStateR = fStatePhi = fStateZ = fChi2 = fNSteps = 0.0;

	fLabel = 0;
	fCount = 0;
	for (i = 0; i < 6; i++) fPoint[i] = track.fPoint[i];

	fVertex.X() = 0.0;
	fVertex.Y() = 0.0;
 	fVertex.Z() = 0.0;
	fVertex.ErrX() = 0.0;
	fVertex.ErrY() = 0.0;
	fVertex.ErrZ() = 0.0;	
}
//
//
//

AliITSNeuralTrack& AliITSNeuralTrack::operator=(const AliITSNeuralTrack& track){
  //assignment operator
  this->~AliITSNeuralTrack();
  new(this) AliITSNeuralTrack(track);
  return *this;
}
AliITSNeuralTrack::~AliITSNeuralTrack()
{
	Int_t l;
	for (l = 0; l < 6; l++) fPoint[l] = 0;
}
//
//
//
void AliITSNeuralTrack::AssignLabel()
{	
// Assigns a GEANT label to the found track. 
// Every cluster has up to three labels (it can have less). Then each label is 
// recorded for each point. Then, counts are made to check if some of the labels
// appear more than once. Finally, the label which appears most times is assigned
// to the track in the field fLabel.
// The number of points containing that label is counted in the fCount data-member. 

	Bool_t found;
	Int_t i, j, l, lab, max = 0;
	
	// We have up to 6 points for 3 labels each => up to 18 possible different values
	Int_t idx[18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Int_t count[18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	for (l = 0; l < 6; l++) {
		if (!fPoint[l]) continue;
		// Sometimes the same label appears two times in the same recpoint.
		// With these if statements, such problem is solved by turning
		// one of them to -1.
		if (fPoint[l]->GetLabel(1) >= 0 && fPoint[l]->GetLabel(1) == fPoint[l]->GetLabel(0))
			fPoint[l]->SetLabel(1, -1);
		if (fPoint[l]->GetLabel(2) >= 0 && fPoint[l]->GetLabel(2) == fPoint[l]->GetLabel(0))
			fPoint[l]->SetLabel(2, -1);
		if (fPoint[l]->GetLabel(2) >= 0 && fPoint[l]->GetLabel(2) == fPoint[l]->GetLabel(1))
			fPoint[l]->SetLabel(2, -1);
		for (i = 0; i < 3; i++) {
			lab = fPoint[l]->GetLabel(i);
			if (lab < 0) continue;
			found = kFALSE;
			for (j = 0; j < max; j++) {
				if (idx[j] == lab) {
					count[j]++;
					found = kTRUE;
				}
			}
			if(!found) {
				max++;
				idx[max - 1] = lab;
				count[max - 1] = 1;
			}
		}
	}
	
	j = 0, max = count[0];
	for (i = 0; i < 18; i++) {
		if (count[i] > max) {
			j = i;
			max = count[i];
		}
	}
	fLabel = idx[j];
	fCount = count[j];
}
//
//
//
void AliITSNeuralTrack::CleanSlot(Int_t i, Bool_t del)
{
// Removes a point from the corresponding layer slot in the found track.
// If the argument is TRUE, the point object is also deleted from heap.

	if (i >= 0 && i < 6) {
		if (del) delete fPoint[i];
		fPoint[i] = 0;
	}
}
//
//
//
void AliITSNeuralTrack::GetModuleData(Int_t layer, Int_t &mod, Int_t &pos)
{
// Returns the point coordinates according to the TreeR philosophy in galice.root files
// that consist in the module number (mod) and the position in the TClonesArray of
// the points reconstructed in that module for the run being examined.

	if (layer < 0 || layer > 5) {
		Error("GetModuleData", "Layer out of range: %d", layer);
		return;
	}
	mod = fPoint[layer]->GetModule();
	pos = fPoint[layer]->GetIndex();
}
//
//
//
void AliITSNeuralTrack::Insert(AliITSNeuralPoint *point)
{
// A trivial method to insert a point in the tracks;
// the point is inserted to the slot corresponding to its ITS layer.

	if (!point) return;
	
	Int_t layer = point->GetLayer();
	if (layer < 0 || layer > 6) {
		Error("Insert", "Layer index %d out of range", layer);
		return;
	} 
	
	fPoint[layer] = point;
}
//
//
//
Int_t AliITSNeuralTrack::OccupationMask() const
{
// Returns a byte which maps the occupied slots. 
// Each bit represents a layer going from the less significant on.

	Int_t i, check, mask = 0;
	for (i = 0; i < 6; i++) {
		check = 1 << i;
		if (fPoint[i]) mask |= check;
	}
	return mask;
}
//
//
//
void AliITSNeuralTrack::PrintLabels()
{
// Prints the results of the AssignLabel() method, together with
// the GEANT labels assigned to each point, in order to evaluate 
// how the assigned label is distributed among points.

	cout << "Assigned label = " << fLabel << " -- counted " << fCount << " times: " << endl << endl;
	for (Int_t i = 0; i < 6; i++) {
		cout << "Point #" << i + 1 << " --> ";
		if (fPoint[i]) {
			cout << "labels = " << fPoint[i]->GetLabel(0) << ", ";
			cout << fPoint[i]->GetLabel(1) << ", ";
			cout << fPoint[i]->GetLabel(2) << endl;
		}
		else {
			cout << "not assigned" << endl;
		}
	}
	cout << endl;
}
//
//
//
Bool_t AliITSNeuralTrack::AddEL(Int_t layer, Double_t sign)
{
// Calculates the correction for energy loss

	Double_t width = 0.0;
	switch (layer) {
		case 0: width = 0.00260 + 0.00283; break;
		case 1: width = 0.0180; break;
		case 2: width = 0.0094; break;
		case 3: width = 0.0095; break;
		case 4: width = 0.0091; break;
		case 5: width = 0.0087; break;
		default:
			Error("AddEL", "Layer value %d out of range!", layer);
			return kFALSE;
	}
	width *= 1.7;

	if((layer == 5) && (fStatePhi < 0.174 || fStatePhi > 6.100 || (fStatePhi > 2.960 && fStatePhi < 3.31))) {
		width += 0.012;
	}

	Double_t invSqCosL = 1. + fTanL * fTanL;            // = 1 / (cos(lambda)^2) = 1 + tan(lambda)^2
	Double_t invCosL   = TMath::Sqrt(invSqCosL);        // = 1  / cos(lambda)
	Double_t pt = GetPt();                              // = transverse momentum
	Double_t p2 = pt *pt * invSqCosL;                   // = square modulus of momentum
	Double_t energy = TMath::Sqrt(p2 + fMass * fMass);  // = energy
	Double_t beta2 = p2 / (p2 + fMass * fMass);         // = (v / c) ^ 2
	if (beta2 == 0.0) {
		printf("Anomaly in AddEL: pt=%8.6f invSqCosL=%8.6f fMass=%8.7f --> beta2 = %8.7f\n", pt, invSqCosL, fMass, beta2);
		return kFALSE;
	}

	Double_t dE = 0.153 / beta2 * (log(5940. * beta2 / (1. - beta2)) - beta2) * width * 21.82 * invCosL;
	dE = sign * dE * 0.001;

	energy += dE;
	p2 = energy * energy - fMass * fMass;
	pt = TMath::Sqrt(p2) / invCosL;
	if (fC < 0.) pt = -pt;
	fC = (0.299792458 * 0.2 * fField) / (pt * 100.);

	return kTRUE;
}
//
//
//
Bool_t AliITSNeuralTrack::AddMS(Int_t layer)
{
// Calculates the noise perturbation due to multiple scattering

	Double_t width = 0.0;
	switch (layer) {
		case 0: width = 0.00260 + 0.00283; break;
		case 1: width = 0.0180; break;
		case 2: width = 0.0094; break;
		case 3: width = 0.0095; break;
		case 4: width = 0.0091; break;
		case 5: width = 0.0087; break;
		default:
			Error("AddEL", "Layer value %d out of range!", layer);
			return kFALSE;
	}
	width *= 1.7;

	if((layer == 5) && (fStatePhi < 0.174 || fStatePhi > 6.100 || (fStatePhi > 2.960 && fStatePhi < 3.31))) {
		width += 0.012;
	}

	Double_t cosL = TMath::Cos(TMath::ATan(fTanL));
	Double_t halfC = fC / 2.;
	Double_t q20 = 1. / (cosL * cosL);
	Double_t q30 = fC * fTanL;

	Double_t q40 = halfC * (fStateR * fStateR - fDt * fDt) / (1. + 2. * halfC * fDt);
	Double_t dd  = fDt + halfC * fDt * fDt - halfC * fStateR * fStateR;
	Double_t dprova = fStateR * fStateR - dd * dd;
	Double_t q41 = 0.;
	if(dprova > 0.) q41 = -1. / cosL * TMath::Sqrt(dprova) / (1. + 2. * halfC *fDt);

	Double_t p2 = (GetPt()*GetPt()) / (cosL * cosL);
	Double_t beta2 = p2 / (p2 + fMass * fMass);
	Double_t theta2 = 14.1 * 14.1 / (beta2 * p2 * 1.e6) * (width / TMath::Abs(cosL));

	fMatrix(2,2) += theta2 * (q40 * q40 + q41 * q41);
	fMatrix(3,2) += theta2 * q20 * q40;
	fMatrix(2,3) += theta2 * q20 * q40;
	fMatrix(3,3) += theta2 * q20 * q20;
	fMatrix(4,2) += theta2 * q30 * q40;
	fMatrix(2,4) += theta2 * q30 * q40;
	fMatrix(4,3) += theta2 * q30 * q20;
	fMatrix(3,4) += theta2 * q30 * q20;
	fMatrix(4,4) += theta2 * q30 * q30;
	
	return kTRUE;
}
//
//
//
Int_t AliITSNeuralTrack::PropagateTo(Double_t rk)
{
	// Propagation method.
	// Changes the state vector according to a new radial position
	// which is specified by the passed 'r' value (in cylindircal coordinates).
	// The covariance matrix is also propagated (and enlarged) according to
	// the FCFt technique, where F is the jacobian of the new parameters
	// w.r.t. their old values.
	// The option argument forces the method to add also the energy loss
	// and the multiple scattering effects, which respectively have the effect
	// of changing the curvature and widening the covariance matrix.

	if (rk < TMath::Abs(fDt)) {
		Error("PropagateTo", Form("Impossible propagation to r (=%17.15g) < Dt (=%17.15g)", rk, fDt));
		return 0;
	}
	
	Double_t duepi = 2. * TMath::Pi();
	Double_t rkm1 = fStateR;
	Double_t aAk = ArgPhi(rk), aAkm1 = ArgPhi(rkm1);
	Double_t ak = ArgZ(rk), akm1 = ArgZ(rkm1);

	fStatePhi += TMath::ASin(aAk) - TMath::ASin(aAkm1);
	if(fStatePhi > duepi) fStatePhi -= duepi;
	if(fStatePhi < 0.) fStatePhi += duepi;

	Double_t halfC = 0.5 * fC;
	fStateZ += fTanL / halfC * (TMath::ASin(ak)-TMath::ASin(akm1));
	
	Double_t bk = ArgB(rk), bkm1 = ArgB(rkm1);
	Double_t ck = ArgC(rk), ckm1 = ArgC(rkm1);
	
	Double_t f02 = ck / TMath::Sqrt(1. - aAk * aAk) - ckm1 / TMath::Sqrt(1. - aAkm1 * aAkm1);
	Double_t f04 = bk / TMath::Sqrt(1. - aAk * aAk) - bkm1 / TMath::Sqrt(1. - aAkm1 * aAkm1);
	Double_t f12 = fTanL * fDt * (1. / rk - 1. / rkm1);
	Double_t f13 = rk - rkm1;
	
	Double_t c00 = fMatrix(0,0);
	Double_t c10 = fMatrix(1,0);
	Double_t c11 = fMatrix(1,1);
	Double_t c20 = fMatrix(2,0);
	Double_t c21 = fMatrix(2,1);
	Double_t c22 = fMatrix(2,2);
	Double_t c30 = fMatrix(3,0);
	Double_t c31 = fMatrix(3,1);
	Double_t c32 = fMatrix(3,2);
	Double_t c33 = fMatrix(3,3);
	Double_t c40 = fMatrix(4,0);
	Double_t c41 = fMatrix(4,1);
	Double_t c42 = fMatrix(4,2);
	Double_t c43 = fMatrix(4,3);
	Double_t c44 = fMatrix(4,4);

	Double_t r10 = c10 + c21*f02 + c41*f04;
	Double_t r20 = c20 + c22*f02 + c42*f04;
	Double_t r30 = c30 + c32*f02 + c43*f04;
	Double_t r40 = c40 + c42*f02 + c44*f04;
	Double_t r21 = c21 + c22*f12 + c32*f13;
	Double_t r31 = c31 + c32*f12 + c33*f13;
	Double_t r41 = c41 + c42*f12 + c43*f13;

	fMatrix(0,0) = c00 + c20*f02 + c40*f04 + f02*r20 + f04*r40;
	fMatrix(1,0) = fMatrix(0,1) = r10 + f12*r20 + f13*r30;
	fMatrix(1,1) = c11 + c21*f12 + c31*f13 + f12*r21 + f13*r31;
	fMatrix(2,0) = fMatrix(0,2) = r20;
	fMatrix(2,1) = fMatrix(1,2) = r21;
	fMatrix(3,0) = fMatrix(0,3) = r30;
	fMatrix(3,1) = fMatrix(1,3) = r31;
	fMatrix(4,0) = fMatrix(0,4) = r40;
	fMatrix(4,1) = fMatrix(1,4) = r41;

	fStateR = rk;

	if (rkm1 < fStateR)  // going to greater R --> energy LOSS
		return -1;
	else                 // going to smaller R --> energy GAIN
		return 1;
}
//
//
//
Bool_t AliITSNeuralTrack::SeedCovariance()
{
	// generate a covariance matrix depending on the results obtained from
	// the preliminary seeding fit procedure.
	// It calculates the variances for C, D ans TanL, according to the
	// differences of the fitted values from the requested ones necessary
	// to make the curve exactly pass through each point.

	/*
	Int_t i, j;
	AliITSNeuralPoint *p = 0;
	Double_t r, argPhi, phiC, phiD, argZ, zL;
	Double_t sumC = 0.0, sumD = 0.0, sumphi = 0., sumz = 0., sumL = 0.;
	for (i = 0; i < fNum; i++) {
		p = At(i);
		if (!p) continue;
		r = p->GetR2();
		// weight and derivatives of phi and zeta w.r.t. various params
		sumphi += 1./ p->ErrorGetPhi();
		argPhi = ArgPhi(r);
		argZ = ArgZ(r);
		if (argPhi > 100.0 || argZ > 100.0) {
			Error("InitCovariance", "Argument error");
			return kFALSE;
		}
		phiC = DerArgPhiC(r) / TMath::Sqrt(1.0 - argPhi * argPhi);
		phiD = DerArgPhiD(r) / TMath::Sqrt(1.0 - argPhi * argPhi);
		if (phiC > 100.0 || phiD > 100.0) {
			Error("InitCovariance", "Argument error");
			return kFALSE;
		}
		zL = asin(argZ) / fC;
		sumL += zL * zL;
		sumC += phiC * phiC;
		sumD += phiD * phiD;
		sumz += 1.0 / (p->fError[2] * p->fError[2]);
	}

	for (i = 0; i < 5; i++) for (j = 0; j < 5; j++) fMatrix(i,j) = 0.;
	fMatrix(0,0) = 1. / sumphi;
	fMatrix(1,1) = 1. / sumz;
	fMatrix(2,2) = 1. / sumD;
	fMatrix(3,3) = 1. / sumL;
	fMatrix(4,4) = 1. / sumC;
	fMatrix.Print();
	*/
	
	AliITSNeuralPoint *p = 0;
	Double_t delta, cs, sn, r, argz;
	Double_t diffC, diffD, diffL, calcC, calcD, calcL;

	Int_t l;
	for (l = 0; l < 6; l++) {
		p = fPoint[l];
		if (!p) break;
		sn = TMath::Sin(p->GetPhi() - fG0);
		cs = TMath::Cos(p->GetPhi() - fG0);
		r  = p->GetR2();
		calcC = (fDt/r - sn) / (2.*fDt*sn - r - fDt*fDt/r);
		argz = ArgZ(r);
		if (argz > 1000.0) {
			Error("Covariance", "Value too high");
			return kFALSE;
		}
		calcL = (p->Z() - fDz) * fC / asin(argz);
		delta = fR*fR + r*r + 2.0*fR*r*sin(p->GetPhi() - fG0);
		if (delta < 0.E0) {
			if (delta >= -0.5)
				delta = 0.;
			else {
				Error("Covariance", Form("Discriminant = %g --- Dt = %g", delta, fDt));
				return kFALSE;
			}
		}
		delta = sqrt(delta);
		if (fC >= 0)
			calcD = delta - fR;
		else
			calcD = fR - delta;
		diffD = calcD - fDt;
		diffL = calcL - fTanL;
		diffC = fC - calcC;
		fMatrix(0,0) += 100000000.0 * p->GetError("phi") * p->GetError("phi");
		fMatrix(1,1) += 10000.0 * p->ErrZ() * p->ErrZ();
		fMatrix(2,2) += 100000.0 * diffD * diffD;
		fMatrix(3,3) += diffL * diffL;
		fMatrix(4,4) += 100000000.0 * diffC * diffC;
	}
	Double_t n = 0.;
	for (l = 0; l < 6; l++) if (fPoint[l]) n++;
	fMatrix *= 1./(n * (n+1));
	return kTRUE;
}
//
//
//
Bool_t AliITSNeuralTrack::Filter(AliITSNeuralPoint *test)
{
	// Makes all calculations which apply the Kalman filter to the
	// stored guess of the state vector, after propagation to a new layer

	if (!test) {
		Error("Filter", "Null pointer passed");
		return kFALSE;
	}
	
	Double_t m[2];
	Double_t rk, phik, zk;
	rk = test->GetR2();
	phik = test->GetPhi();
	zk = test->Z();
	m[0]=phik;
	m[1]=zk;
	
	//////////////////////// Evaluation of the error matrix V  /////////
	Double_t v00 = test->GetError("phi") * rk;
	Double_t v11 = test->ErrZ();
	////////////////////////////////////////////////////////////////////  
	
	// Get the covariance matrix
	Double_t cin00, cin10, cin20, cin30, cin40;
	Double_t cin11, cin21, cin31, cin41, cin22;
	Double_t cin32, cin42, cin33, cin43, cin44;
	cin00 = fMatrix(0,0);
	cin10 = fMatrix(1,0);
	cin20 = fMatrix(2,0);
	cin30 = fMatrix(3,0);
	cin40 = fMatrix(4,0);
	cin11 = fMatrix(1,1);
	cin21 = fMatrix(2,1);
	cin31 = fMatrix(3,1);
	cin41 = fMatrix(4,1);
	cin22 = fMatrix(2,2);
	cin32 = fMatrix(3,2);
	cin42 = fMatrix(4,2);
	cin33 = fMatrix(3,3);
	cin43 = fMatrix(4,3);
	cin44 = fMatrix(4,4);
	
	// Calculate R matrix
	Double_t rold00 = cin00 + v00;
	Double_t rold10 = cin10;
	Double_t rold11 = cin11 + v11;
	
	////////////////////// R matrix inversion  /////////////////////////
	Double_t det = rold00*rold11 - rold10*rold10;
	Double_t r00 = rold11/det;
	Double_t r10 = -rold10/det;
	Double_t r11 = rold00/det;
	////////////////////////////////////////////////////////////////////
	
	// Calculate Kalman matrix
	Double_t k00 = cin00*r00 + cin10*r10;
	Double_t k01 = cin00*r10 + cin10*r11;
	Double_t k10 = cin10*r00 + cin11*r10;  
	Double_t k11 = cin10*r10 + cin11*r11;
	Double_t k20 = cin20*r00 + cin21*r10;  
	Double_t k21 = cin20*r10 + cin21*r11;  
	Double_t k30 = cin30*r00 + cin31*r10;  
	Double_t k31 = cin30*r10 + cin31*r11;  
	Double_t k40 = cin40*r00 + cin41*r10;
	Double_t k41 = cin40*r10 + cin41*r11;
	
	// Get state vector (will keep the old values for phi and z)
	Double_t x0, x1, x2, x3, x4, savex0, savex1;
	x0 = savex0 = fStatePhi;
	x1 = savex1 = fStateZ;
	x2 = fDt;
	x3 = fTanL;
	x4 = fC;

	// Update the state vector
	x0 += k00*(m[0]-savex0) + k01*(m[1]-savex1);
	x1 += k10*(m[0]-savex0) + k11*(m[1]-savex1);
	x2 += k20*(m[0]-savex0) + k21*(m[1]-savex1);
	x3 += k30*(m[0]-savex0) + k31*(m[1]-savex1);
	x4 += k40*(m[0]-savex0) + k41*(m[1]-savex1);
	
	// Update the covariance matrix
	Double_t cout00, cout10, cout20, cout30, cout40;
	Double_t cout11, cout21, cout31, cout41, cout22;
	Double_t cout32, cout42, cout33, cout43, cout44;
	
	cout00 = cin00 - k00*cin00 - k01*cin10;
	cout10 = cin10 - k00*cin10 - k01*cin11;
	cout20 = cin20 - k00*cin20 - k01*cin21;
	cout30 = cin30 - k00*cin30 - k01*cin31;
	cout40 = cin40 - k00*cin40 - k01*cin41;
	cout11 = cin11 - k10*cin10 - k11*cin11;
	cout21 = cin21 - k10*cin20 - k11*cin21;
	cout31 = cin31 - k10*cin30 - k11*cin31;
	cout41 = cin41 - k10*cin40 - k11*cin41;
	cout22 = cin22 - k20*cin20 - k21*cin21;
	cout32 = cin32 - k20*cin30 - k21*cin31;
	cout42 = cin42 - k20*cin40 - k21*cin41;
	cout33 = cin33 - k30*cin30 - k31*cin31;
	cout43 = cin43 - k30*cin40 - k31*cin41;
	cout44 = cin44 - k40*cin40 - k41*cin41;
	
	// Store the new covariance matrix
	fMatrix(0,0) = cout00;
	fMatrix(1,0) = fMatrix(0,1) = cout10;
	fMatrix(2,0) = fMatrix(0,2) = cout20;
	fMatrix(3,0) = fMatrix(0,3) = cout30;
	fMatrix(4,0) = fMatrix(0,4) = cout40;
	fMatrix(1,1) = cout11;
	fMatrix(2,1) = fMatrix(1,2) = cout21;
	fMatrix(3,1) = fMatrix(1,3) = cout31;
	fMatrix(4,1) = fMatrix(1,4) = cout41;
	fMatrix(2,2) = cout22;
	fMatrix(3,2) = fMatrix(2,3) = cout32;
	fMatrix(4,2) = fMatrix(2,4) = cout42;
	fMatrix(3,3) = cout33;
	fMatrix(4,3) = fMatrix(3,4) = cout43;
	fMatrix(4,4) = cout44;
	
	// Calculation of the chi2 increment
	Double_t vmcold00 = v00 - cout00;
	Double_t vmcold10 = -cout10;
	Double_t vmcold11 = v11 - cout11;
	////////////////////// Matrix vmc inversion  ///////////////////////
	det = vmcold00*vmcold11 - vmcold10*vmcold10;
	Double_t vmc00=vmcold11/det;
	Double_t vmc10 = -vmcold10/det;
	Double_t vmc11 = vmcold00/det;
	////////////////////////////////////////////////////////////////////
	Double_t chi2 = (m[0] - x0)*( vmc00*(m[0] - x0) + 2.*vmc10*(m[1] - x1) ) + (m[1] - x1)*vmc11*(m[1] - x1);
	fChi2 += chi2;
	fNSteps++;
	
	return kTRUE;
}
//
//
//
Bool_t AliITSNeuralTrack::KalmanFit()
{
// Applies the Kalman Filter to improve the track parameters resolution.
// First, thre point which lies closer to the estimated helix is chosen.
// Then, a fit is performed towards the 6th layer
// Finally, the track is refitted to the 1st layer

	Double_t rho;
	Int_t l, layer, sign;
	
	fStateR = fPoint[0]->GetR2();
	fStatePhi = fPoint[0]->GetPhi();
	fStateZ = fPoint[0]->Z();
	
	if (!PropagateTo(3.0)) {
		Error("KalmanFit", "Unsuccessful initialization");
		return kFALSE;
	}
	l=0;

	// Performs a Kalman filter going from the actual state position
	// towards layer 6 position
	// Now, the propagation + filtering operations can be performed
	Double_t argPhi = 0.0, argZ = 0.0;
	while (l <= 5) {
		if (!fPoint[l]) {
			Error("KalmanFit", "Not six points!");
			return kFALSE;
		}
		layer = fPoint[l]->GetLayer();
		rho = fPoint[l]->GetR2();
		sign = PropagateTo(rho);
		if (!sign) return kFALSE;
		AddEL(layer, -1.0);
		AddMS(layer);
		if (!Filter(fPoint[l])) return kFALSE;
		// these two parameters are update according to the filtered values
		argPhi = ArgPhi(fStateR);
		argZ = ArgZ(fStateR);
		if (argPhi > 1.0 || argPhi < -1.0 || argZ > 1.0 || argZ < -1.0) {
			Error("Filter", Form("Filtering returns too large values: %g, %g", argPhi, argZ));
			return kFALSE;
		}
		fG0 = fStatePhi - asin(argPhi);
		fDz = fStateZ - (2.0 * fTanL / fC) * asin(argZ);
		l++;
	}

	// Now a Kalman filter i performed going from the actual state position
	// towards layer 1 position and then propagates to vertex
	if (l >= 5) l = 5;
	while (l >= 1) {
		layer = fPoint[l]->GetLayer();
		rho = fPoint[l]->GetR2();
		AddEL(layer, 1.0);
		sign = PropagateTo(rho);
		if (!sign) return kFALSE;
		AddMS(layer);
		if (!Filter(fPoint[l])) return kFALSE;
		// these two parameters are update according to the filtered values
		argPhi = ArgPhi(fStateR);
		argZ = ArgZ(fStateR);
		if (argPhi > 1.0 || argPhi < -1.0 || argZ > 1.0 || argZ < -1.0) {
			Error("Filter", Form("Filtering returns too large values: %g, %g", argPhi, argZ));
			return kFALSE;
		}
		fG0 = fStatePhi - asin(argPhi);
		fDz = fStateZ - (2.0 * fTanL / fC) * asin(argZ);
		l--;
	}
	return kTRUE;
}
//
//
//
Bool_t AliITSNeuralTrack::RiemannFit()
{
	// Method which executes the circle fit via a Riemann Sphere projection
	// with the only improvement of a weighted mean, due to different errors
	// over different point measurements.
	// As an output, it returns kTRUE or kFALSE respectively if the fit succeeded or not
	// in fact, if some variables assume strange values, the fit is aborted,
	// in order to prevent the class from raising a floating point error;

	Int_t i, j;

	// M1 - matrix of ones
	TMatrixD m1(6,1);
	for (i = 0; i < 6; i++) m1(i,0) = 1.0;

	// X - matrix of Rieman projection coordinates
	TMatrixD mX(6,3);
	for (i = 0; i < 6; i++) {
		mX(i,0) = fPoint[i]->X();
		mX(i,1) = fPoint[i]->Y();
		mX(i,2) = fPoint[i]->GetR2sq();
	}

	// W - matrix of weights
	Double_t xterm, yterm, ex, ey;
	TMatrixD mW(6,6);
	for (i = 0; i < 6; i++) {
		xterm = fPoint[i]->X() * fPoint[i]->GetPhi() - fPoint[i]->Y() / fPoint[i]->GetR2();
		ex = fPoint[i]->ErrX();
		yterm = fPoint[i]->Y() * fPoint[i]->GetPhi() + fPoint[i]->X() / fPoint[i]->GetR2();
		ey = fPoint[i]->ErrY();
		mW(i,i) = fPoint[i]->GetR2sq() / (xterm * xterm * ex * ex + yterm * yterm * ey * ey);
	}

	// Xm - weighted sample mean
	Double_t meanX = 0.0, meanY = 0.0, meanW = 0.0, sw = 0.0;
	for (i = 0; i < 6; i++) {
		meanX += mW(i,i) * mX(i,0);
		meanY += mW(i,i) * mX(i,1);
		meanW += mW(i,i) * mX(i,2);
		sw += mW(i,i);
	}
	meanX /= sw;
	meanY /= sw;
	meanW /= sw;

	// V - sample covariance matrix
	for (i = 0; i < 6; i++) {
		mX(i,0) -= meanX;
		mX(i,1) -= meanY;
		mX(i,2) -= meanW;
	}
	TMatrixD mXt(TMatrixD::kTransposed, mX);
	TMatrixD mWX(mW, TMatrixD::kMult, mX);
	TMatrixD mV(mXt, TMatrixD::kMult, mWX);
	for (i = 0; i < 3; i++) {
		for (j = i + 1; j < 3; j++) {
			mV(i,j)  = mV(j,i)  = (mV(i,j) + mV(j,i)) * 0.5;
		}
	}

	// Eigenvalue problem solving for V matrix
	Int_t ileast = 0;
	TVectorD eval(3), n(3);
	//	TMatrixD evec = mV.EigenVectors(eval);
#if ROOT_VERSION_CODE >= 262146
	TMatrixDEigen ei(mV);
	TMatrixD evec = ei.GetEigenVectors();
	eval = ei.GetEigenValues();
#else
	TMatrixD evec = mV.EigenVectors(eval);
#endif

	if (eval(1) < eval(ileast)) ileast = 1;
	if (eval(2) < eval(ileast)) ileast = 2;
	n(0) = evec(0, ileast);
	n(1) = evec(1, ileast);
	n(2) = evec(2, ileast);

	// c - known term in the plane intersection with Riemann axes
	Double_t c = -(meanX * n(0) + meanY * n(1) + meanW * n(2));

	fXC = -n(0) / (2. * n(2));
	fYC = -n(1) / (2. * n(2));
	fR  = (1. - n(2)*n(2) - 4.*c*n(2)) / (4. * n(2) * n(2));

	if (fR <= 0.E0) {
		Error("RiemannFit", "Radius comed less than zero!!!");
		return kFALSE;
	}
	fR = TMath::Sqrt(fR);
	fC = 1.0 / fR;

	// evaluating signs for curvature and others
	Double_t phi1 = 0.0, phi2, temp1, temp2, sumdphi = 0.0, ref = TMath::Pi();
	AliITSNeuralPoint *p = fPoint[0];
	phi1 = p->GetPhi();
	for (i = 1; i < 6; i++) {
		p = (AliITSNeuralPoint*)fPoint[i];
		if (!p) break;
		phi2 = p->GetPhi();
		temp1 = phi1;
		temp2 = phi2;
		if (temp1 > ref && temp2 < ref)
			temp2 += 2.0 * ref;
		else if (temp1 < ref && temp2 > ref)
			temp1 += 2.0 * ref;
		sumdphi += temp2 - temp1;
		phi1 = phi2;
	}
	if (sumdphi < 0.E0) fC = -fC;
	Double_t diff, angle = TMath::ATan2(fYC, fXC);
	if (fC < 0.E0)
		fG0 = angle + 0.5 * TMath::Pi();
	else
		fG0 = angle - 0.5 * TMath::Pi();
	diff = angle - fG0;

	Double_t d = TMath::Sqrt(fXC*fXC + fYC*fYC) - fR;
	if (fC >= 0.E0)
		fDt = d;
	else
		fDt = -d;

	Int_t nn = 6;
	Double_t halfC = 0.5 * fC;
	Double_t *s = new Double_t[nn], *z = new Double_t[nn], *ws = new Double_t[nn];
	for (j = 0; j < 6; j++) {
		p = fPoint[j];
		if (!p) break;
		s[j] = ArgZ(p->GetR2());
		if (s[j] > 100.0) return kFALSE;
		z[j] = p->Z();
		s[j] = asin(s[j]) / halfC;
		ws[j] = 1.0 / (p->ErrZ()*p->ErrZ());
	}

	// second tep final fit
	Double_t sums2 = 0.0, sumz = 0.0, sumsz = 0.0, sums = 0.0, sumw = 0.0;
	for (i = 0; i < nn; i++) {
		sums2 += ws[i] * s[i] * s[i];
		sumz  += ws[i] * z[i];
		sums  += ws[i] * s[i];
		sumsz += ws[i] * s[i] * z[i];
		sumw += ws[i];
	}
	sums2 /= sumw;
	sumz /= sumw;
	sums /= sumw;
	sumsz /= sumw;
	d = sums2 - sums*sums;

	fDz = (sums2*sumz - sums*sumsz) / d;
	fTanL = (sumsz - sums*sumz) / d;

	delete [] s;
	delete [] z;
	delete [] ws;

	return kTRUE;
}
//
//
//
void AliITSNeuralTrack::PrintState(Bool_t matrix)
{
// Prints the state vector values.
// The argument switches on or off the printing of the covariance matrix.

	cout << "\nState vector: " << endl;
	cout << " Rho = " << fStateR << "\n";
	cout << " Phi = " << fStatePhi << "\n";
	cout << "   Z = " << fStateZ << "\n";
	cout << "  Dt = " << fDt << "\n";
	cout << "  Dz = " << fDz << "\n";
	cout << "TanL = " << fTanL << "\n";
	cout << "   C = " << fC << "\n";
	cout << "  G0 = " << fG0 << "\n";
	cout << "  XC = " << fXC << "\n";
	cout << "  YC = " << fYC << "\n";
	if (matrix) {
		cout << "\nCovariance Matrix: " << endl;
		fMatrix.Print();
	}
	cout << "Actual square chi = " << fChi2;
}
//
//
//
Double_t AliITSNeuralTrack::GetDz() const
{
//	Double_t argZ = ArgZ(fStateR);
//	if (argZ > 9.9) {
//		Error("GetDz", "Too large value: %g", argZ);
//		return 0.0;
//	}
//	fDz = fStateZ - (2.0 * fTanL / fC) * asin(argZ);
	return fDz;
}
//
//
//
Double_t AliITSNeuralTrack::GetGamma() const
{
// these two parameters are update according to the filtered values
//	Double_t argPhi = ArgPhi(fStateR);
//	if (argPhi > 9.9) {
//		Error("Filter", "Too large value: %g", argPhi);
//		return kFALSE;
//	}
//	fG0 = fStatePhi - asin(argPhi);
	return fG0;
}
//
//
//
Double_t AliITSNeuralTrack::GetPhi(Double_t r) const
{
// Gives the value of azymuthal coordinate in the helix
// as a function of cylindric radius

	Double_t arg = ArgPhi(r);
	if (arg > 0.9) return 0.0;
	arg = fG0 + asin(arg);
	while (arg >= 2.0 * TMath::Pi()) { arg -= 2.0 * TMath::Pi(); }
	while (arg < 0.0) { arg += 2.0 * TMath::Pi(); }
	return arg;
}
//
//
//
Double_t AliITSNeuralTrack::GetZ(Double_t r) const
{
// gives the value of Z in the helix
// as a function of cylindric radius

	Double_t arg = ArgZ(r);
	if (arg > 0.9) return 0.0;
	return fDz + fTanL * asin(arg) / fC;
}
//
//
//
Double_t AliITSNeuralTrack::GetdEdX()
{
// total energy loss of the track

	Double_t q[4] = {0., 0., 0., 0.}, dedx = 0.0;
        Int_t i = 0, swap = 0;
        for (i = 2; i < 6; i++) {
                if (!fPoint[i]) continue;
                q[i - 2] = (Double_t)fPoint[i]->GetCharge();
              	q[i - 2] /= (1 + fTanL*fTanL);
        }
	q[0] /= 280.;
	q[1] /= 280.;
	q[2] /= 38.;
	q[3] /= 38.;
	do {
		swap = 0;
		for (i = 0; i < 3; i++) {
			if (q[i] <= q[i + 1]) continue;
			Double_t tmp = q[i];
			q[i] = q[i + 1]; 
			q[i+1] = tmp;
			swap++;
		}
	} while(swap); 
	if(q[0] < 0.) {
		q[0] = q[1];
		q[1] = q[2];
		q[2] = q[3];
		q[3] = -1.;
	} 
	dedx = (q[0] + q[1]) / 2.;
        return dedx;
}
//
//
//
void AliITSNeuralTrack::SetVertex(Double_t *pos, Double_t *err)
{
	// Stores vertex data

	if (!pos || !err) return;
	fVertex.ErrX() = err[0];
	fVertex.ErrY() = err[1];
	fVertex.ErrZ() = err[2];
	fVertex.SetLayer(0);
	fVertex.SetModule(0);
	fVertex.SetIndex(0);
	fVertex.SetLabel(0, -1);
	fVertex.SetLabel(1, -1);
	fVertex.SetLabel(2, -1);
	fVertex.SetUser(1);
}
//
//
//
AliITSIOTrack* AliITSNeuralTrack::ExportIOtrack(Int_t min)
{
// Exports an object in the standard format for reconstructed tracks

	Int_t layer = 0;
	AliITSIOTrack *track = new AliITSIOTrack;

	// covariance matrix
	track->SetCovMatrix(fMatrix(0,0), fMatrix(1,0), fMatrix(1,1),
	                    fMatrix(2,0), fMatrix(2,1), fMatrix(2,2),
			    fMatrix(3,0), fMatrix(3,1), fMatrix(3,2),
			    fMatrix(3,3), fMatrix(4,0), fMatrix(4,1),
			    fMatrix(4,2), fMatrix(4,3), fMatrix(4,4));
	
	// labels
	track->SetLabel(IsGood(min) ? fLabel : -fLabel);
	track->SetTPCLabel(-1);

	// points characteristics
	for (layer = 0; layer < 6; layer++) {
		if (fPoint[layer]) {
			track->SetIdModule(layer, fPoint[layer]->GetModule());
			track->SetIdPoint(layer, fPoint[layer]->GetIndex());
		}
	}
	
	// state vector
	track->SetStatePhi(fStatePhi);
	track->SetStateZ(fStateZ);
	track->SetStateD(fDt);
	track->SetStateTgl(fTanL);
	track->SetStateC(fC);
	track->SetRadius(fStateR);
	track->SetCharge((fC > 0.0) ? -1 : 1);
	track->SetDz(fDz);

	// track parameters in the closest point
	track->SetX(fStateR * cos(fStatePhi));
	track->SetY(fStateR * cos(fStatePhi));
	track->SetZ(fStateZ);
	track->SetPx(GetPt() * cos(fG0));
	track->SetPy(GetPt() * sin(fG0));
	track->SetPz(GetPt() * fTanL);

	// PID
	track->SetPid(fPDG);
	track->SetMass(fMass);

	return track;
}
//
//
//====================================================================================
//============================       PRIVATE METHODS      ============================
//====================================================================================
//
//
Double_t AliITSNeuralTrack::ArgPhi(Double_t r) const
{
	// calculates the expression ((1/2)Cr + (1 + (1/2)CD) D/r) / (1 + CD)

	Double_t arg, num, den;
	num = (0.5 * fC * r) + (1. + (0.5 * fC * fDt)) * (fDt / r);
	den = 1. + fC * fDt;
	if (den == 0.) {
		Error("ArgPhi", "Denominator = 0!");
		return 10.0;
	}
	arg = num / den;
	if (TMath::Abs(arg) < 1.) return arg;
	if (TMath::Abs(arg) <= 1.00001) return (arg > 0.) ? 0.99999999999 : -0.9999999999;
	Error("ArgPhi", "Value too large: %17.15g", arg);
	return 10.0;
}
//
//
//
Double_t AliITSNeuralTrack::ArgZ(Double_t r) const
{
	// calculates the expression (1/2)C * sqrt( (r^2 - Dt^2) / (1 + CD) )

	Double_t arg;
	arg = (r * r - fDt * fDt) / (1. + fC * fDt);
	if (arg < 0.) {
		if (TMath::Abs(arg) < 1.E-6) arg = 0.;
		else {
			Error("ArgZ", "Square root argument error: %17.15g < 0", arg);
			return 10.;
		}
	}
	arg = 0.5 * fC * TMath::Sqrt(arg);
	if (TMath::Abs(arg) < 1.) return arg;
	if (TMath::Abs(arg) <= 1.00001) return (arg > 0.) ? 0.99999999999 : -0.9999999999;
	Error("ArgZ", "Value too large: %17.15g", arg);
	return 10.0;
}
//
//
//
Double_t AliITSNeuralTrack::ArgB(Double_t r) const 
{
// UTILITY FUNCTION 

	Double_t arg;
	arg = (r*r - fDt*fDt);
	arg /= (r*(1.+ fC*fDt)*(1.+ fC*fDt));
	return arg;
}
//
//
//
Double_t AliITSNeuralTrack::ArgC(Double_t r) const 
{
// UTILITY FUNCTION

	Double_t arg;
	arg = (1./r - fC * ArgPhi(r));
	arg /= 1.+ fC*fDt;
	return arg;
}
