#include <Riostream.h>
#include <Riostream.h>
//#include <stdlib.h>

//#include <TObject.h>
//#include <TROOT.h>
#include <TMath.h>
#include <TString.h>

#include "AliITSglobalRecPoint.h"
//#include "AliITSneuralTracker.h"

#include "AliITSneuralTrack.h"

///////////////////////////////////////////////////////////////////////
//
// AliITSneuralTrack: 
// track object produced with a
// neural network MFT algorithm 
// for track finding in ITS stand alone
// according to the Denby-Peterson model with adaptments to the 
// ALICE multiplicity
// 
///////////////////////////////////////////////////////////////////////

ClassImp(AliITSneuralTrack)
//
//
//
AliITSneuralTrack::AliITSneuralTrack()
{
  // default constructor
	Int_t i;
	for (i = 0; i < 6; i++) fPoint[i] = 0;
}
//__________________________________________________________
AliITSneuralTrack::AliITSneuralTrack(const AliITSneuralTrack &src):TObject(src){
  // copy constructor
  fFitXC = src.fFitXC;
  fFitYC = src.fFitYC;
  fFitRadius = src.fFitRadius;
  fFitTanL = src.fFitTanL;
  fSqChi = src.fSqChi;
  for(Int_t i=0;i<6;i++)fPoint[i] = src.fPoint[i];
}

//__________________________________________________________
AliITSneuralTrack& AliITSneuralTrack::operator=(const AliITSneuralTrack &src){
  // assignment operator
  if(this != &src){
    TObject::operator=(src);
    fFitXC = src.fFitXC;
    fFitYC = src.fFitYC;
    fFitRadius = src.fFitRadius;
    fFitTanL = src.fFitTanL;
    fSqChi = src.fSqChi;
    for(Int_t i=0;i<6;i++)fPoint[i] = src.fPoint[i];
  }
  return *this;
}

//
//
//
AliITSneuralTrack::~AliITSneuralTrack() 
{
	Int_t i;
	for (i = 0; i < 6; i++) delete fPoint[i];
}
//
//
//
Int_t AliITSneuralTrack::CheckMe(Bool_t verbose) const {
  // returns the number of points (it should be 6)
	Int_t l, stored = 0;
	TString empty("Not filled slots: ");
	for (l = 0; l < 6; l++) {
		if (fPoint[l]) 
			stored++;
		else {
			empty += l;
			empty += ' ';
		}
	}
		
	if (stored < 6 && verbose) Warning("", empty);
	return stored;
}
//
//
//
Int_t AliITSneuralTrack::EvaluateTrack(Bool_t verbose, Int_t min, Int_t* &good){
  //Used for comparison (good and fake tracks)
	Int_t i, j, k = 0, count[18], id[18];
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 3; j++) {
			if (fPoint[i])
				id[k] = fPoint[i]->fLabel[j];
			else
				id[k] = -1;
			count[k] = 0;
			k++;
		}
	}
	for (i = 0; i < 18; i++) {
		for (j = i+1; j < 18; j++) {
			if (id[i] == id[j] || id[j] < 0) id[j] = -1;
		}
	}
	for (i = 0; i < 18; i++) {
		if (id[i] < 0) continue;
		for (j = 0; j < 6; j++) {
			if (fPoint[j] && fPoint[j]->HasID(id[i])) count[i]++;
		}
	}
	Int_t index[18], best;
	TMath::Sort(18, count, index);
	best = id[index[0]];
	if (verbose) {
		if (count[index[0]]<min) cout << "\t";
		cout << best << " " << count[index[0]] << endl;
	}
		
	if (good) delete [] good;
	good = new Int_t[6];
	for (i = 0; i < 6; i++) {
		good[i] = 0;
		if (!fPoint[i])
			good[i] = -1;
		else if (fPoint[i] && fPoint[i]->HasID(best))
			good[i] = 1;
		else if (fPoint[i]->fLabel[0] < 0 && fPoint[i]->fLabel[1] < 0 && fPoint[i]->fLabel[2] < 0)
			good[i] = -1;
	}
	
	if (count[index[0]] < min) best = -best;
	
	return best;
	
	
	
	/*
	if (good) delete [] good;
	good = new Int_t[6];
	Int_t count = CheckMe(verbose);
	
	if (count < min) return -9999999;
	
	Int_t i, l, max = 0, best = 0;
	for (l = 0; l < 6; l++) {
		for (i = 0; i < 3; i++) {
			if (fPoint[l].fLabel[i] > max) max = fPoint[l].fLabel[i];
		}
	}
	
	count = max + 1;
	Int_t *counter = new Int_t[count];
	for (i = 0; i < count; i++) counter[i] = 0;
	
	for (l = 0; l < 6; l++) {
		for (i = 0; i < 3; i++) {
			if (fPoint[l].fLabel[i] >= 0) 
				counter[fPoint[l].fLabel[i]]++;
		}
	}
	
	for (i = 0; i < count; i++) {
		if (counter[i] > counter[best]) best = i;
	}
	
	for (l = 0; l < 6; l++) {
		good[l] = 0;
		if (fPoint[l].fLabel[0] == best || fPoint[l].fLabel[1] == best || fPoint[l].fLabel[2] == best)
			good[l] = 1;
		else if (fPoint[l].fLabel[0] < 0 && fPoint[l].fLabel[1] < 0 && fPoint[l].fLabel[2] < 0)
			good[l] = -1;
	}
	
	if (counter[best] < min) best = -best;
	delete counter;
	return best;
	*/
}
//
//
//
void AliITSneuralTrack::GetCoords(Double_t* &x, Double_t* &y, Double_t* &z) const
{
  // returns the coordinates of a point
	if (x) delete [] x; x = new Double_t[6];
	if (y) delete [] y; y = new Double_t[6];
	if (z) delete [] z; z = new Double_t[6];
	
	Int_t i;
	
	for (i = 0; i < 6; i++) {
		x[i] = y[i] = z[i] = 0.0;
		if (!fPoint[i]) continue;
		x[i] = fPoint[i]->fGX;
		y[i] = fPoint[i]->fGY;
		z[i] = fPoint[i]->fGZ;
	}
}
//
//
//
void AliITSneuralTrack::CopyPoint(AliITSglobalRecPoint *p)
{
  // copy p in a new globalrecpoint obj
	Int_t layer = p->fLayer;
	if (layer < 0 || layer > 6) {
		Error("", Form("Wrong layer [%d]", layer));
		return;
	}
	
	fPoint[layer] = new AliITSglobalRecPoint;
	fPoint[layer]->fGX = p->fGX;
	fPoint[layer]->fGY = p->fGY;
	fPoint[layer]->fGZ = p->fGZ;
	fPoint[layer]->fGSX = p->fGSX;
	fPoint[layer]->fGSY = p->fGSY;
	fPoint[layer]->fGSZ = p->fGSZ;
	fPoint[layer]->fLayer = layer;
	fPoint[layer]->fLabel[0] = p->fLabel[0];
	fPoint[layer]->fLabel[1] = p->fLabel[1];
	fPoint[layer]->fLabel[2] = p->fLabel[2];
}
//
//
//
void AliITSneuralTrack::Print(Option_t *option, Int_t min)
{
  // print this
	Int_t *vuoto = 0;
	TString opt(option);
	opt.ToUpper();
	Int_t id = EvaluateTrack(0, min, vuoto);
	if (opt.Contains("A")) {
		cout << "\nEvaluated ID for this track: ";
		cout.width(8);
		if (id >= 0) {
			cout << id;
			cout << " [good]";
		}
		else {
			cout << -id;
			cout << " [fake]";
		}
	}
	cout << endl << endl;
}
//
//
//
void AliITSneuralTrack::Kinks(Int_t &pos, Int_t &neg, Int_t &incr, Int_t &decr)
{
  // Kinks
	Int_t i;
	Double_t dphi, dphiold = 0.0;
	pos = neg = incr = decr = 0;
	for (i = 1; i < 6; i++) {
		dphi = fPoint[i]->fPhi - fPoint[i-1]->fPhi;
		if (dphi > 0.0) pos++; else neg++;
		if (TMath::Abs(dphi) > dphiold) incr++; else decr++;
		dphiold = TMath::Abs(dphi);
	}
}
//
//
//
Double_t AliITSneuralTrack::FitXY(Double_t VX, Double_t VY)
{
  // determines XC, YC, radius and chi squared
	Int_t i;
	Double_t xX, yY, dD, rR;
	Double_t rx(0.0), ry(0.0), x2(0.0), y2(0.0), xy(0.0);
	for (i = 0; i < 6; i++) {
		xX = fPoint[i]->fGX - VX;
		yY = fPoint[i]->fGY - VY;
		rR = xX * xX + yY * yY;
		rx += rR * xX;
		ry += rR * yY;
		x2 += xX * xX;
		y2 += yY * yY;
		xy += xX * yY;
	}
  
	dD = 2 * (x2 * y2 - xy * xy);
	if (dD == 0.0) 
		return 1000.0;
	else {
		xX = (rx * y2 - ry * xy) / dD;
		yY = (ry * x2 - rx * xy) / dD;
		fFitRadius  = TMath::Sqrt(xX * xX + yY * yY);
		fFitXC = xX + VX;
		fFitYC = yY + VY;
	}
  
	fSqChi = 0.0;
	for (i = 0; i < 6; i++) {
		xX = fPoint[i]->fGX - fFitXC;
		yY = fPoint[i]->fGY - fFitYC;
		fSqChi += ((xX * xX + yY * yY) / (fFitRadius * fFitRadius)) * ((xX * xX + yY * yY) / (fFitRadius * fFitRadius));
	}
	fSqChi /= 6.0;
	fSqChi = TMath::Sqrt(fSqChi);
	return fSqChi;
}
