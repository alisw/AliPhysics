/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Alberto Pulvirenti.                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 * track finder for ITS stand-alone with neural network algorithm         *
 * This class defines a Hopfield MFT neural network simulation            *
 * which reads all recpoints from an event and produces a tree with       *
 * the points belonging to recognized tracks                              *
 * the TTree obtained as the output file will be saved in a .root file    * 
 * and it must be read by groups of six entries at a time                 *
 **************************************************************************/

#include <Riostream.h>
#include <Riostream.h>
#include <stdlib.h>

#include <TROOT.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSgeomMatrix.h"
#include "AliITSRecPoint.h"
#include "AliITSglobalRecPoint.h"
#include "AliITSneuralTrack.h"
#include "AliITSneuralTracker.h"

ClassImp(AliITSneuralTracker)


// ====================================================================================================
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AliITSneuralTracker <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// ====================================================================================================

AliITSneuralTracker::AliITSneuralTracker() 
{
	//
	// Default constructor (not for use)
	//
	Int_t i = 0;
	fCurvCut = 0;
	for (i = 0; i < 6; i++) {
		fPoints[i] = 0;
		if (i < 5) fNeurons[i] = 0;
	}
	fTracks = 0;
}
//---------------------------------------------------------------------------------------------------------
AliITSneuralTracker::AliITSneuralTracker (Int_t nsecs, Int_t ncurv, Double_t *curv, Double_t theta)
{ 
	// Argumented constructor
	// ----------------------
	// gets the number of azymuthal sectors,
	// the number of curvature cut steps and their cut values, 
	// and then the polar angle cut.
	// Be careful not to put 'ncurv' greater than the dimension of the 
	// 'curv' array (--> segmentation fault)
	
	fSectorNum   = nsecs;
	fSectorWidth = 2.0 * TMath::Pi()/(Double_t)nsecs;
		
	fCurvNum = ncurv;
	fCurvCut = new Double_t[ncurv];
	Int_t i;
	for (i = 0; i < ncurv; i++) fCurvCut[i] = curv[i];
	
	fThetaCut = theta;
	fThetaNum = (Int_t)(TMath::Pi() / theta) + 1;
	if (fThetaNum < 1) fThetaNum = 1;
		
	cout << "\nClass AliITSneuralTracker invoked with " << fSectorNum << " phi secs.";
	cout << " and " << fThetaNum << " theta secs.\n" << endl;
	
	Int_t k, j;
	for (k = 0; k < 6; k++) {
		fPoints[k] = new TObjArray**[fSectorNum];
		for (i = 0; i < fSectorNum; i++) {
			fPoints[k][i] = new TObjArray*[fThetaNum];
			for (j = 0; j < fThetaNum; j++) fPoints[k][i][j] = new TObjArray(0);
		}
		if (k < 6) fNeurons[k] = new TObjArray(0);
	}
	
	fTracks = new TObjArray(1);
}
//---------------------------------------------------------------------------------------------------------
AliITSneuralTracker::~AliITSneuralTracker() 
{
	// Destructor
	// ----------
	// clears the TObjArrays and all heap objects
	
	Int_t lay, sec, k;
	delete [] fCurvCut;
	
	for (lay = 0; lay < 6; lay++) {
		for (sec = 0; sec < fSectorNum; sec++) {
			for (k = 0; k < fThetaNum; k++) {
				fPoints[lay][sec][k]->Delete();
			}
		}
		delete [] fPoints[lay];
		if (lay < 5) {
			fNeurons[lay]->Delete();
			delete fNeurons[lay];
		}
	}
	
	fTracks->Delete();
}
//---------------------------------------------------------------------------------------------------------
Int_t AliITSneuralTracker::ReadFile(const Text_t *fname, Int_t evnum)
{
	// File reader
	// -----------
	// Reads a ROOT file in order to retrieve recpoints.
	// the 'evnum' argument can have two meanings:
	// if it is >= 0 it represents the event number to retrieve the correct gAlice
	// if it is < 0, instead, this means that the open file contains 
	// a 'TreeP' tree, with all points which remained unused 
	// after the Kalman tracking procedure (produced via the ITSafterKalman.C macro)
	// the return value is the number of collected points 
	// (then, if errors occur, the method returns 0)
	
	TFile *file = new TFile(fname, "read");
	TTree *tree = 0;
	Int_t i, nentries, total, sector, mesh;
	
	if (evnum < 0) {
		tree = (TTree*)file->Get("TreeP0");
		if (!tree) {
			Error("", "Specifying evnum < 0, I expect to find a 'TreeP' tree in the file, but there isn't any");
			return 0;
		}
		AliITSglobalRecPoint *p = 0;
		tree->SetBranchAddress("Points", &p);
		nentries = (Int_t)tree->GetEntries();
		total = 0; 
		for (i = 0; i < nentries; i++) {
			tree->GetEntry(i);
			AliITSglobalRecPoint *q = new AliITSglobalRecPoint(*p);
			sector = (Int_t)(q->fPhi / fSectorWidth);
			mesh = (Int_t)(q->fTheta / fThetaCut);
			q->fUsed = 0;
			fPoints[q->fLayer][sector][mesh]->AddLast(q);
			total++;
		}
		return total;
	}
	
	if (gAlice) { delete gAlice; gAlice = 0; }
	gAlice = (AliRun*)file->Get("gAlice");
	if (!gAlice) {
		Error("", "gAlice is NULL, maybe wrong filename...");
		return 0;
	}

	Int_t nparticles = (Int_t)gAlice->GetEvent(evnum);
	if (!nparticles) {
		Error("", "No particles found");
		return 0;
	}

	AliITS *its = (AliITS*)gAlice->GetModule("ITS");
	if (!its) {
		Error("", "No ITS found");
		return 0;
	}

	AliITSgeom *geom = (AliITSgeom*)its->GetITSgeom();
	if (!geom) {
		Error("", "AliITSgeom is NULL");
		return 0;
	}

	tree = gAlice->TreeR();
	if (!tree) {
		Error("", "TreeR() returned NULL");
		return 0;
	}

	nentries = (Int_t)tree->GetEntries();
	if (!nentries) {
		Error("", "TreeR is empty");
		return 0;
	}

	// Reading objects initialization
	Int_t k, entry, lay, lad, det;
	Double_t l[3], g[3], ll[3][3], gg[3][3];
	
	TObjArray *recsArray = its->RecPoints();
	AliITSRecPoint *recp;
	AliITSglobalRecPoint *pnt = 0;
	
	total = 0;
	for(entry = 0; entry < nentries; entry++) {
		if (entry > geom->GetLastSSD()) continue;
		its->ResetRecPoints();
		tree->GetEvent(entry);
		TObjArrayIter recs(recsArray);
		AliITSgeomMatrix *gm = geom->GetGeomMatrix(entry);
		geom->GetModuleId(entry, lay, lad, det);
		lay--;
		while ((recp = (AliITSRecPoint*)recs.Next())) {
			// conversion to global coordinate system
			for (i = 0; i < 3; i++) {
				l[i] = g[i] = 0.0;
				for (k = 0; k < 3; k++) {
					ll[i][k] = 0.0;
					gg[i][k] = 0.0;
				}
			}
			// local to global conversions of coords
			l[0] = recp->fX;
			l[1] = 0.0;
			l[2] = recp->fZ;
			gm->LtoGPosition(l, g);
			// local to global conversions of sigmas
			ll[0][0] = recp->fSigmaX2;
			ll[2][2] = recp->fSigmaZ2;
			gm->LtoGPositionError(ll, gg);
			// instantiation of a new global rec-point
			pnt = new AliITSglobalRecPoint(g[0], g[1], g[2], gg[0][0], gg[1][1], gg[2][2], lay);
			// copy of owner track labels
			for (i = 0; i < 3; i++) pnt->fLabel[i] = recp->fTracks[i];
			sector = (Int_t)(pnt->fPhi / fSectorWidth);
			mesh = (Int_t)(pnt->fTheta / fThetaCut);
			pnt->fUsed = 0;
			fPoints[lay][sector][mesh]->AddLast(pnt);
			total++;
		}
	}
	
	return total;
}
//---------------------------------------------------------------------------------------------------------
void AliITSneuralTracker::Go(const char* filesave, Bool_t flag)
{
	// Global working method
	// ---------------------
	// It's the method which calls all other ones,
	// in order to analyze the whole ITS points ensemble
	// It's the only method to use, besides the parameter setters,
	// so, all other methods are defined as private
	// (see header file)
	// the flag meaning is explained in the 'Initialize' method
	
	Int_t i; 
	for (i = 0; i < fSectorNum; i++) DoSector(i, flag);
	cout << endl << "Saving results in file " << filesave << "..." << flush;
	TFile *f = new TFile(filesave, "recreate");
	fTracks->Write();
	cout << "done" << endl;
	f->Close();
}

//=========================================================================================================
// * * *   P R I V A T E   S E C T I O N   * * *
//=========================================================================================================

Int_t AliITSneuralTracker::Initialize(Int_t secnum, Int_t curvnum, Bool_t flag)
{
	// Network initializer
	// -------------------
	// Fills the neuron arrays, following the cut criteria for the selected step
	// (secnum = sector to analyze, curvnum = curvature cut step to use)
	// It also sets the initial random activation.
	// In order to avoid a large number of 'new' operations, all existing neurons
	// are reset and initialized with the new values, and are created new unit only if
	// it turns out to be necessary
	// the 'flag' argument is used to decide if the lower edge in the intevral
	// of the accepted curvatures is given by zero (kFALSE) or by the preceding used cut (kTRUE)
	// (if this is the first step, anyway, the minimum is always zero)
	
	Int_t l, m, k, neurons = 0;
	
	for (l = 0; l < 5; l++) fNeurons[l]->Delete();
	
	AliITSneuron *unit = 0;
	Double_t abscurv, max = fCurvCut[curvnum], min = 0.0;
	if (flag && curvnum > 0) min = fCurvCut[curvnum - 1];
	AliITSglobalRecPoint *inner = 0, *outer = 0;
	for (l = 0; l < 5; l++) {
		for (m = 0; m < fThetaNum; m++) {
			TObjArrayIter inners(fPoints[l][secnum][m]);
			while ((inner = (AliITSglobalRecPoint*)inners.Next())) {
				if (inner->fUsed > 0) continue; // points can't be recycled
				for (k = m-1; k <= m+1; k++) {
					if (k < 0 || k >= fThetaNum) continue; // to avoid seg faults 
					TObjArrayIter outers(fPoints[l+1][secnum][k]);
					while ((outer = (AliITSglobalRecPoint*)outers.Next())) {
						if (outer->fUsed > 0) continue;  // points can't be recycled
						if (inner->DTheta(outer) > fThetaCut) continue;
						unit = new AliITSneuron;
						unit->fInner = inner;
						unit->fOuter = outer;
						CalcParams(unit);
						abscurv = TMath::Abs(unit->fCurv);
						if (unit->fDiff > fDiff || abscurv < min || abscurv > max) {
							delete unit;
							continue;
						}
						unit->fActivation = gRandom->Rndm() * (fMax - fMin) + fMin;
						fNeurons[l]->AddLast(unit);
						neurons++;
					} // end loop on candidates for outer point for neurons
				}
			} // end loop on candidates inner points for neuron
		}
	} // end loop on layers
	return neurons;
}
//---------------------------------------------------------------------------------------------------------
Bool_t AliITSneuralTracker::Update()
{
	// Updating cycle
	// --------------
	// Performs an updating cycle, by summing all the
	// gain and cost contribution for each neuron and
	// updating its activation.
	// An asynchronous method is followed, with the order
	// according that the neurons have been created
	// Time wasting is avoided by dividing the neurons 
	// into different arrays, depending on the layer of their inner point.
	// The return value answers the question: "The network has stabilized?"
	
	Int_t j, l;
	Double_t actOld, actNew, argNew, totDiff = 0.0;
	Double_t sumInh, sumExc, units = 0.0;
	AliITSneuron *me = 0, *it = 0;
	for (l = 0; l < 5; l++) {
		TObjArrayIter meIter(fNeurons[l]);
		while ((me = (AliITSneuron*)meIter.Next())) {
			units++;
			TObjArrayIter inhIter(fNeurons[l]);
			
			// inhibition (search only for neurons starting in the same layer)
			sumInh = 0.0;
			while((it = (AliITSneuron*)inhIter.Next())) {
				if (it->fOuter == me->fOuter || it->fInner == me->fInner)	
					sumInh += it->fActivation;
			}
			
			// excitation (search only for neurons 
			// which start in the previous or next layers)
			sumExc = 0.0;
			for (j = l - 1; j <= l + 1; j += 2) {
				if (j < 0 || j >= 5) continue;
				TObjArrayIter itIter(fNeurons[j]);
				while ((it = (AliITSneuron*)itIter.Next())) {
					if (it->fInner == me->fOuter || it->fOuter == me->fInner && it->fCurv * me->fCurv > 0.0) {
						sumExc += Weight(me, it) * it->fActivation;
					}
				}
			} // end search for excitories
			
			actOld = me->fActivation;
			argNew = fRatio * sumExc - sumInh;
			actNew = 1.0 / (1.0 + TMath::Exp(-argNew / fTemp));
			me->fActivation = actNew;
			// calculation the relative activation variation 
			// (for stabilization check)
			totDiff += TMath::Abs((actNew - actOld) / actOld);
		} // end loop over 'me' (updated neuron)
	} // end loop over layers
	totDiff /= units;
	
	return (totDiff < fVar);
}
//---------------------------------------------------------------------------------------------------------
Int_t AliITSneuralTracker::Save()
{
	// Tracki saving method
	// --------------------
	// Stores the recognized tracks into the TObjArray 'fTracks'
	// but doesn't save it into a file until the whole work has ended
	
	Int_t l;
	Double_t test = 0.5;
	
	// Before all, for each neuron is found the best active sequences
	// among the incoming and outgoing other units
	AliITSneuron *main, *fwd, *back, *start;
	for (l = 0; l < 5; l++) {
		TObjArrayIter mainIter(fNeurons[l]);
		while((main = (AliITSneuron*)mainIter.Next())) {
			if (l < 4) {
				TObjArrayIter fwdIter(fNeurons[l+1]);
				test = 0.5;
				while ((fwd = (AliITSneuron*)fwdIter.Next())) {
					if (main->fOuter == fwd->fInner && fwd->fActivation > test) {
						test = fwd->fActivation;
						main->fNext = fwd;
					}
				};
			}
			if (l > 0) {
				TObjArrayIter backIter(fNeurons[l-1]);
				test = 0.5;
				while ((back = (AliITSneuron*)backIter.Next())) {
					if (main->fInner == back->fOuter && back->fActivation > test) {
						test = back->fActivation;
						main->fPrev = back;
					}
				};
			}
		}
	}
	
	// Then, are resolved the mismatches in the chains found above 
	// (when unit->next->prev != unit or unit->prev->next != unit)
	for (l = 0; l < 5; l++) {
		TObjArrayIter mainIter(fNeurons[l]);
		while((main = (AliITSneuron*)mainIter.Next())) {
			fwd  = main->fNext;
			back = main->fPrev;
			if (fwd  != 0 && fwd->fPrev  != main) main->fNext = 0; 
			if (back != 0 && back->fNext != main) main->fPrev = 0; 
		}
	}
	
	Int_t pos, neg, incr, decr, ntracks;
	AliITSneuralTrack *trk = 0;
		
	TObjArrayIter startIter(fNeurons[0]);
	for (;;) {
		start = (AliITSneuron*)startIter.Next();
		if (!start) break;
		Int_t pts = 0;
		// with the chain above defined, a track can be followed
		// like a linked list structure
		for (main = start; main; main = main->fNext) pts++;
		// the count will be  > 5 only if the track contains 6 points
		// (what we want)
		if (pts < 5) continue;
		// track storage
		trk = new AliITSneuralTrack;
		trk->CopyPoint(start->fInner);
		for (main = start; main; main = main->fNext) {
			//main->fInner->fUsed = kTRUE;
			//main->fOuter->fUsed = kTRUE;
			trk->CopyPoint(main->fOuter);
		}
		trk->Kinks(pos, neg, incr, decr);
		if (pos != 0 && neg != 0) {
			cout << "pos, neg kinks = " << pos << " " << neg << endl;
			continue;
		}
		if (incr != 0 && decr != 0) {
			cout << "increment, decrements in delta phi = " << incr << " " << decr << endl;
			continue;
		}
		for (main = start; main; main = main->fNext) {
			main->fInner->fUsed++;
			main->fOuter->fUsed++;
		}
		fTracks->AddLast(trk);
	}
	ntracks = (Int_t)fTracks->GetEntriesFast();
	return ntracks;
}
//---------------------------------------------------------------------------------------------------------
void AliITSneuralTracker::DoSector(Int_t sect, Bool_t flag)
{
	// Sector recognition
	// ------------------
	// This is a private method which works on a single sector
	// just for an easier readability of the code.
	// The 'flag' argument is needed to be given to the 'Initialize' method (see it)
	
	Int_t curvnum, nTracks = 0, nUnits;
	Bool_t isStable = kFALSE;
	for (curvnum = 0; curvnum < fCurvNum; curvnum++) {
		cout << "\rSector " << sect << ":" << ((curvnum+1)*100)/fCurvNum << "%" << flush;
		nUnits = Initialize(sect, curvnum, flag);
		if (!nUnits) continue;
		do {
			isStable = Update();
		} while (!isStable);
		nTracks = Save();
	}
	cout << "\rTracks stored after sector " << sect << ": " << nTracks << endl;
}
//---------------------------------------------------------------------------------------------------------
Int_t AliITSneuralTracker::CountExits(AliITSneuron *n)
{
	// Counter for neurons which enter in 'n' tail
	
	Int_t count = 0, l = n->fOuter->fLayer;
	if (l == 5) return 0;
	
	AliITSneuron *test = 0;
	TObjArrayIter iter(fNeurons[l]);
	while ((test = (AliITSneuron*)iter.Next())) {
		if (test->fInner == n->fOuter) count++;
	}
	return count;
}
//---------------------------------------------------------------------------------------------------------
Int_t AliITSneuralTracker::CountEnters(AliITSneuron *n)
{
	// Counter for neurons which exit from 'n' head
	
	Int_t count = 0, l = n->fInner->fLayer;
	if (l == 0) return 0;
	
	AliITSneuron *test = 0;
	TObjArrayIter iter(fNeurons[l-1]);
	while ((test = (AliITSneuron*)iter.Next())) {
		if (test->fOuter == n->fInner) count++;
	}
	return count;
}
//---------------------------------------------------------------------------------------------------------
Bool_t AliITSneuralTracker::ResetNeuron(AliITSneuron *n)
{
	// A method which resets the neuron's pointers to zero
	
	if (!n) return kFALSE;
	n->fNext = 0;
	n->fPrev = 0;
	n->fInner = 0;
	n->fOuter = 0;
	// the other datamembers don't need to be reinitialized
	return kTRUE;
}
//---------------------------------------------------------------------------------------------------------
Bool_t AliITSneuralTracker::SetEdge(AliITSneuron *n, AliITSglobalRecPoint *p, const char what)
{
	// Sets the indicated edge point to the neuron
	// (if the arguments are suitable ones...)
	
	if (!n || !p || (what != 'i' && what != 'I' && what != 'o' && what != 'O')) return kFALSE;
	if (what == 'i' || what == 'I') 
		n->fInner = p;
	else
		n->fOuter = p;
	return kTRUE;
}
//---------------------------------------------------------------------------------------------------------
Bool_t AliITSneuralTracker::CalcParams(AliITSneuron *n)
{
	// This method evaluates the helix parameters from the edged of a neuron
	// (curvature, phase parameter, tangent of lambda angle)
	// if the p[arameters assume too strange (and unphysical) values, the 
	// method return kFALSE, else it return kTRUE;
	
	if (!n) return kFALSE;
	
	Double_t sign(1.0), det, rc = 0.0;
	Double_t tanl1, tanl2, l1, l2;
	// changing the reference frame into the one centered in the vertex
	Double_t x1 = n->fInner->fGX - fVPos[0];
	Double_t y1 = n->fInner->fGY - fVPos[1];
	Double_t z1 = n->fInner->fGZ - fVPos[2];
	Double_t r1 = x1 * x1 + y1 * y1;
	Double_t x2 = n->fOuter->fGX - fVPos[0];
	Double_t y2 = n->fOuter->fGY - fVPos[1];
	Double_t z2 = n->fOuter->fGZ - fVPos[2];
	Double_t r2 = x2 * x2 + y2 * y2;
	// initialization of the XY-plane data (curvature and centre)
	// (will remain these values if we encounter errors)
	n->fCurv = n->fCX = n->fCY = n->fTanL = n->fPhase = 0.0;
	n->fLength = (n->fOuter->fGX - n->fInner->fGX) * (n->fOuter->fGX - n->fInner->fGX);
	n->fLength += (n->fOuter->fGY - n->fInner->fGY) * (n->fOuter->fGY - n->fInner->fGY);
	n->fLength += (n->fOuter->fGZ - n->fInner->fGZ) * (n->fOuter->fGZ - n->fInner->fGZ);
	n->fLength = TMath::Sqrt(n->fLength);
	// this determinant which is zero when the points are in a line
	det = 2.0 * (x1 * y2 - x2 * y1);
	if (det == 0.0)
		return kFALSE; // it is difficult having perfectly aligned points --- it's more likely to be an error
	else {
		n->fCX = (r1*y2 - r2*y1) / det;
		n->fCY = (r2*x1 - r1*x2) / det;
		rc = TMath::Sqrt(n->fCX * n->fCX + n->fCY * n->fCY);
		r1 = TMath::Sqrt(r1);
		r2 = TMath::Sqrt(r2);
		sign = (n->fOuter->fPhi >= n->fInner->fPhi) ? 1.0 : -1.0;
		if (rc > 0.0) {
			n->fCurv = sign / rc;
			n->fPhase = TMath::ATan2(-n->fCX, -n->fCY);
			if (n->fPhase < 0.0) n->fPhase += 2.0 * TMath::Pi(); // angles in 0 --> 2pi
			if (r1 <= 2.0 * rc && r2 <= 2.0 * rc) {
				l1 = 2.0 * rc * TMath::ASin(r1 / (2.0 * rc));
				l2 = 2.0 * rc * TMath::ASin(r2 / (2.0 * rc));
				tanl1 = z1 / l1;
				tanl2 = z2 / l2;
				n->fDiff = TMath::Abs(tanl1 - tanl2);
				n->fTanL = 0.5 * (tanl1 + tanl2);
			}
			else {
				n->fTanL = 0.0;
				n->fDiff = 100.0;
				return kFALSE; // it' even more difficult that the arguments above aren't acceptable for arcsine
			}
		}
		else
			return kFALSE;
	}
	return kTRUE;
}
//---------------------------------------------------------------------------------------------------------
Double_t AliITSneuralTracker::Angle(AliITSneuron *n, AliITSneuron *m)
{
	// calculates the angle between two segments
	// but doesn't check if they have a common point.
	// The calculation is made by the ratio of their scalar product
	// to the product of their magnitudes
	
	Double_t x = (m->fOuter->fGX - m->fInner->fGX) * (n->fOuter->fGX - n->fInner->fGX);
	Double_t y = (m->fOuter->fGY - m->fInner->fGY) * (n->fOuter->fGY - n->fInner->fGY);
	Double_t z = (m->fOuter->fGZ - m->fInner->fGZ) * (n->fOuter->fGZ - n->fInner->fGZ);
		
	Double_t cosine = (x + y + z) / (n->fLength * m->fLength);
	
	if (cosine > 1.0 || cosine < -1.0) {
		Warning("AliITSneuronV1::Angle", "Strange value of cosine");
		return 0.0;
	}
	
	return TMath::ACos(cosine);
}
//---------------------------------------------------------------------------------------------------------
Double_t AliITSneuralTracker::Weight(AliITSneuron *n, AliITSneuron *m)
{
	// calculation of the excitoy weight
	
//	Double_t v1 = TMath::Abs(fTanL - n->fTanL);
//	Double_t v2 = TMath::Abs(TMath::Abs(fCurv) - TMath::Abs(n->fCurv));
//	Double_t v3 = TMath::Abs(fPhase - n->fPhase);
//	Double_t q = (v1 + v2 + v3) / 3.0;
//	Double_t a = 1.0; // 0.02 / (fDiff + n->fDiff);
	Double_t b = 1.0 - TMath::Sin(Angle(n, m));
	return TMath::Power(b, fExp);
//	return TMath::Exp(q / par1);

}
//---------------------------------------------------------------------------------------------------------
// ====================================================================================================
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AliITSneuron <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// ====================================================================================================

AliITSneuron::AliITSneuron()
{
	// Default constructor
	// -------------------
	// sets all parameters to default (but unuseful) 
	// (in particular all pointers are set to zero)
	// There is no alternative constructor, because 
	// this class should not be used by itself, but 
	// it is only necessary to allow the neural tracker to work
	// In fact, all methods which operate on neuron's datamembers
	// are stored in the AliITSneuralTracker class anyway
	
	fActivation = 0.0;
	
	fCurv = 0.0;
	fCX = 0.0;
	fCY = 0.0;
	fTanL = 0.0;
	fPhase = 0.0;
	fDiff = 0.0;
	fLength = 0.0;
	
	fNext = fPrev = 0;
	fInner = fOuter = 0;
}
//---------------------------------------------------------------------------------------------------------
AliITSneuron::~AliITSneuron() 
{
	// Destructor
	// ----------
	// does nothing, because there is no need to affect parameters,
	// while the AliITSglobalRecPoint pointers can't be deleted,
	// and the AliITSneuron pointers belong to a TObjArray which will
	// be cleared when necessary...
}

