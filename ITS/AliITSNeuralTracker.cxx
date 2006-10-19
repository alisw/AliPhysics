///**************************************************************************
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: Alberto Pulvirenti.                                            *
// *                                                                        *
// * Permission to use, copy, modify and distribute this software and its   *
// * documentation strictly for non-commercial purposes is hereby granted   *
// * without fee, provided that the above copyright notice appears in all   *
// * copies and that both the copyright notice and this permission notice   *
// * appear in the supporting documentation. The authors make no claims     *
// * about the suitability of this software for any purpose.                *
// * It is provided "as is" without express or implied warranty.            *
// *                                                                        *
// * AN ITS STAND-ALONE "NEURAL" TRACK FINDER                               *
// * ----------------------------------------                               *
// * This class implements the Denby-Peterson algorithm for track finding   *
// * in the ITS stand-alone, by means of a neural network simulation.       *
// * Many parameters have to be set for the neural network to operate       *
// * correctly and with a good efficiency.                                  *
// * The neural tracker must be feeded with a TTree filled with objects     *
// * of the class "AliITSNeuralPoint", into a single branch called       *
// * "Points".                                                              *
// **************************************************************************/

//#include <fstream>
#include <Riostream.h>
#include <stdlib.h>

//#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLine.h>
//#include <TMarker.h>
#include <TRandom.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
//#include <TParticle.h>
#include <TObjArray.h>
//#include <TList.h>

#include "AliITSNeuralTracker.h"

//using namespace std;
ClassImp(AliITSNeuralTracker)

//--------------------------------------------------------------------------------------------

AliITSNeuralTracker::AliITSNeuralTracker():
fSectorNum(1),
fSectorWidth(0),
fPolarInterval(45.0),
fCurvNum(0),
fCurvCut(0),
fStructureOK(kFALSE),
fVX(0.0),
fVY(0.0),
fVZ(0.0),
fActMinimum(0.5),
fEdge1(0.3),
fEdge2(0.7),
fTemperature(1.0),
fStabThreshold(0.001),
fGain2CostRatio(1.0),
fAlignExponent(1.0),
fChains(0),
fNeurons(0){
	// CONSTRUCTOR
	//
	// Initializes some data-members:
	// - all pointers to NULL
	// - theta cut to 180 deg. (= no theta cut)
	// - initial choice of only 1 azimuthal sector
	// - initial values for the neural network parameters.
	//
	// With these settings the neural tracker can't work
	// because it has not any curvature cut set.
	

	fSectorWidth = TMath::Pi() * 2.0;

	Int_t ilayer, itheta;
	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) fPoints[ilayer][itheta] = 0;
		fThetaCut2DMin[ilayer] = 0.0;
		fThetaCut2DMax[ilayer] = TMath::Pi();
		fThetaCut3DMin[ilayer] = 0.0;
		fThetaCut3DMax[ilayer] = TMath::Pi();
		fHelixMatchCutMin[ilayer] = 1.0;
		fHelixMatchCutMax[ilayer] = 1.0;
	}




	fChains = new TTree("TreeC", "Sextines of points");
	fChains->Branch("l0", &fPoint[0], "l0/I");
	fChains->Branch("l1", &fPoint[1], "l1/I");
	fChains->Branch("l2", &fPoint[2], "l2/I");
	fChains->Branch("l3", &fPoint[3], "l3/I");
	fChains->Branch("l4", &fPoint[4], "l4/I");
	fChains->Branch("l5", &fPoint[5], "l5/I");
}

AliITSNeuralTracker::AliITSNeuralTracker(const AliITSNeuralTracker &n):TObject(n),
fSectorNum(n.fSectorNum),
fSectorWidth(n.fSectorWidth),
fPolarInterval(n.fPolarInterval),
fCurvNum(n.fCurvNum),
fCurvCut(n.fCurvCut),
fStructureOK(n.fStructureOK),
fVX(n.fVX),
fVY(n.fVY),
fVZ(n.fVZ),
fActMinimum(n.fActMinimum),
fEdge1(n.fEdge1),
fEdge2(n.fEdge2),
fTemperature(n.fTemperature),
fStabThreshold(n.fStabThreshold),
fGain2CostRatio(n.fGain2CostRatio),
fAlignExponent(n.fAlignExponent),
fChains(n.fChains),
fNeurons(n.fNeurons){
  //copy contructor
}

AliITSNeuralTracker& AliITSNeuralTracker::operator=(const AliITSNeuralTracker& t){
  //assignment operator
  this->~AliITSNeuralTracker();
  new(this) AliITSNeuralTracker(t);
  return *this;
}
//
//--------------------------------------------------------------------------------------------
//
AliITSNeuralTracker::~AliITSNeuralTracker()
{
	// DESTRUCTOR
	//
	// It Destroys all the dynamic arrays and
	// clears the TCollections and the points tree
	
	delete [] fCurvCut;
	
	Int_t ilayer, itheta;
	if (fStructureOK) {
		for (ilayer = 0; ilayer < 6; ilayer++) {
			for (itheta = 0; itheta < 180; itheta++) {
				fPoints[ilayer][itheta]->SetOwner();
				delete fPoints[ilayer][itheta];
			}
		}
		fNeurons->SetOwner();
		delete fNeurons;
	}
	
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::Display(TCanvas*& canv) const
{
	// Displays the neural segments
	
	Double_t x1, y1, x2, y2;
	canv->Clear();
	TObjArrayIter iter(fNeurons);
	for (;;) {
		AliITSneuron *unit = (AliITSneuron*)iter.Next();
		if (!unit) break;
		if (unit->Activation() < fActMinimum) continue;
		x1 = unit->Inner()->X();
		x2 = unit->Outer()->X();
		y1 = unit->Inner()->Y();
		y2 = unit->Outer()->Y();
		TLine *line = new TLine(x1, y1, x2, y2);
		canv->cd();
		line->Draw();
	}
	canv->Update();
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::SetThetaCuts2D(Double_t *min, Double_t *max)
{
	// Cut setter
	
	Int_t i;
	Double_t temp;
	for (i = 0; i < 5; i++) {
		if (min[i] > max[i]) {
			temp = min[i];
			min[i] = max[i];
			max[i] = temp;
		}
		fThetaCut2DMin[i] = min[i];
		fThetaCut2DMax[i] = max[i];
	}
	for (i = 0; i < 5; i++) {
		cout << "Theta 2D cut for layer " << i << " = " << fThetaCut2DMin[i] << " to " << fThetaCut2DMax[i] << endl;
	}
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::SetThetaCuts3D(Double_t *min, Double_t *max)
{
	// Cut setter
	
	Int_t i;
	Double_t temp;
	for (i = 0; i < 5; i++) {
		if (min[i] > max[i]) {
			temp = min[i];
			min[i] = max[i];
			max[i] = temp;
		}
		fThetaCut3DMin[i] = min[i];
		fThetaCut3DMax[i] = max[i];
	}
	for (i = 0; i < 5; i++) {
		cout << "Theta 3D cut for layer " << i << " = " << fThetaCut3DMin[i] << " to " << fThetaCut3DMax[i] << endl;
	}
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::SetHelixMatchCuts(Double_t *min, Double_t *max)
{
	// Cut setter
	
	Int_t i;
	Double_t temp;
	for (i = 0; i < 5; i++) {
		if (min[i] > max[i]) {
			temp = min[i];
			min[i] = max[i];
			max[i] = temp;
		}
		fHelixMatchCutMin[i] = min[i];
		fHelixMatchCutMax[i] = max[i];
	}
	for (i = 0; i < 5; i++) {
		cout << "Helix-match cut for layer " << i << " = " << fHelixMatchCutMin[i] << " to " << fHelixMatchCutMax[i] << endl;
	}
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::SetCurvatureCuts(Int_t ncuts, Double_t *curv)
{
	// CURVATURE CUTS SETTER
	//
	// Requires an array of double values and its dimension
	// After sorting it in increasing order, the array of curvature cuts
	// is dinamically allocated, and filled with the sorted cuts array.
	// A message is shown which lists all the curvature cuts.

	Int_t i, *ind = new Int_t[ncuts];
	TMath::Sort(ncuts, curv, ind, kFALSE);
	fCurvCut = new Double_t[ncuts];
	cout << "\n" << ncuts << " curvature cuts defined" << endl << "-----" << endl;
	for (i = 0; i < ncuts; i++) {
		fCurvCut[i] = curv[ind[i]];
		cout << "Cut #" << i + 1 << " : " << fCurvCut[i] << endl;
	}
	cout << "-----" << endl;
	fCurvNum = ncuts;
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::CreateArrayStructure(Int_t nsectors)
{
	// ARRAY CREATOR
	//
	// Organizes the array structure to store all points in.
	//
	// The array is organized into a "multi-level" TCollection:
	// - 6 fPoints[] TObjArray containing a TObjArray for each layer
	// - each TObject contains a TObjArray for each sector.

	// sets the number of sectors and their width.
	fSectorNum   = nsectors;
	fSectorWidth = TMath::Pi() * 2.0 / (Double_t)fSectorNum;

	// creates the TCollection structure
	Int_t ilayer, isector, itheta;
	TObjArray *sector = 0;
	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			if (fPoints[ilayer][itheta]) {
				fPoints[ilayer][itheta]->SetOwner();
				delete fPoints[ilayer][itheta];
			}
			fPoints[ilayer][itheta] = new TObjArray(nsectors);
			for (isector = 0; isector < nsectors; isector++) {
				sector = new TObjArray;
				sector->SetOwner();
				fPoints[ilayer][itheta]->AddAt(sector, isector);
			}
		}
	}

	// Sets a checking flag to TRUE. 
	// This flag is checked before filling up the arrays with the points.
	fStructureOK = kTRUE;
}
//
//--------------------------------------------------------------------------------------------
//
Int_t AliITSNeuralTracker::ArrangePoints(TTree* ptstree)
{
	// POINTS STORAGE INTO ARRAY
	//
	// Reads points from the tree and creates AliITSNode objects for each one,
	// storing them into the array structure defined above.
	// Returns the number of points collected (if successful) or 0 (otherwise)

	// check: if the points tree is NULL or empty, there is nothing to do...
	if ( !ptstree || (ptstree && !(Int_t)ptstree->GetEntries()) ) {
		Error("ArrangePoints", "Points tree is NULL or empty: no points to arrange");
		return 0;
	}

	if (!fStructureOK) {
		Error("ArrangePoints", "Structure NOT defined. Call CreateArrayStructure() first");
		return 0;
	}

	Int_t isector, itheta, ientry, ilayer, nentries, pos;
	TObjArray *sector = 0;
	AliITSNode *created = 0;
	AliITSNeuralPoint *cursor = 0;

	ptstree->SetBranchAddress("pos", &pos);
	ptstree->SetBranchAddress("Points", &cursor);
	nentries = (Int_t)ptstree->GetEntries();

	for (ientry = 0; ientry < nentries; ientry++) {
		ptstree->GetEntry(ientry);
		// creates the object
		created = new AliITSNode(cursor, kTRUE);
		created->SetUser(-1);
		created->PosInTree() = pos;
		// finds the sector in phi
		isector = created->GetSector(fSectorWidth);
		itheta  = created->GetThetaCell();
		ilayer  = created->GetLayer();
		if (ilayer < 0 || ilayer > 5) {
			Error("ArrangePoints", "Layer value %d not allowed. Aborting...", ilayer);
			return 0;
		}
		// selects the right TObjArray to store object into
		sector = (TObjArray*)fPoints[ilayer][itheta]->At(isector);
		sector->AddLast(created);
	}

	// returns the number of points processed
	return ientry;
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::PrintPoints()
{
	// creates the TCollection structure
	TObjArray *sector = 0;
	Int_t ilayer, isector, itheta;
	fstream file("test.log", ios::out);
	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (isector = 0; isector < fSectorNum; isector++) {
			for (itheta = 0; itheta < 180; itheta++) {
				sector = (TObjArray*)fPoints[ilayer][itheta]->At(isector);
				file << ilayer << " " << isector << " " << itheta;
				file << " " << sector->GetSize() << " points" << endl;
			}
		}
	}
	file.close();
}
//
//--------------------------------------------------------------------------------------------
//
Bool_t AliITSNeuralTracker::PassCurvCut
(AliITSNode *p1, AliITSNode *p2, Int_t curvindex, Double_t vx, Double_t vy, Double_t vz)
{
	// CURVATURE CUT EVALUATOR
	//
	// Checks the passsed point pair w.r. to the current curvature cut
	// Returns the result of the check.

	if (curvindex < 0 || curvindex >= fCurvNum) {
		Error("PassCurvCut", "Curv index %d out of range", curvindex);
		return kFALSE;
	}
	
	// Find the reference layer
	Int_t lay1 = p1->GetLayer();
	Int_t lay2 = p2->GetLayer();
	Int_t reflayer = (lay1 < lay2) ? lay1 : lay2;

	Double_t x1 = p1->X() - vx;
	Double_t x2 = p2->X() - vx;
	Double_t y1 = p1->Y() - vy;
	Double_t y2 = p2->Y() - vy;
	Double_t z1 = p1->Z() - vz;
	Double_t z2 = p2->Z() - vz;
	Double_t r1 = sqrt(x1*x1 + y1*y1);
	Double_t r2 = sqrt(x2*x2 + y2*y2);

	// calculation of curvature
	Double_t dx = p1->X() - p2->X(), dy = p1->Y() - p2->Y();
	Double_t num = 2 * (x1*y2 - x2*y1);
	Double_t den = r1*r2*sqrt(dx*dx + dy*dy);
	Double_t curv = 0.;
	/* FOR OLD VERSION
	if (den != 0.) {
		curv = fabs(num / den);
		if (curv > fCurvCut[curvindex]) return kFALSE;
		return kTRUE;
	}
	else
		return kFALSE;
	*/
	// NEW VERSION
	if (den != 0.) {
		curv = TMath::Abs(num / den);
		if (curv > fCurvCut[curvindex]) return kFALSE;
	}
	else
		return kFALSE;
	// calculation of helix matching
	Double_t arc1 = 2.0 * r1 * curv;
	Double_t arc2 = 2.0 * r2 * curv;
	Double_t helmatch = 0.0;
	if (arc1 > -1.0 && arc1 < 1.0) arc1 = asin(arc1);
	else arc1 = ((arc1 > 0.0) ? 0.5 : 1.5) * TMath::Pi();
	if (arc2 > -1.0 && arc2 < 1.0) arc2 = asin(arc2);
	else arc2 = ((arc2 > 0.0) ? 0.5 : 1.5) * TMath::Pi();
	arc1 /= 2.0 * curv;
	arc2 /= 2.0 * curv;
	if (arc1 == 0.0 || arc2 == 0.0) return kFALSE;
	helmatch = TMath::Abs(z1 / arc1 - z2 / arc2);
	return (helmatch >= fHelixMatchCutMin[reflayer] && helmatch <= fHelixMatchCutMax[reflayer]);
	// END NEW VERSION
}
//
//--------------------------------------------------------------------------------------------
//
Int_t AliITSNeuralTracker::PassAllCuts
(AliITSNode *p1, AliITSNode *p2, Int_t curvindex, Double_t vx, Double_t vy, Double_t vz)
{
	// GLOBAL CUT EVALUATOR
	//
	// Checks all cuts for the passed point pair.
	// Return values:
	// 
	// 0 - All cuts passed
	// 1 - theta 2D cut not passed
	// 2 - theta 3D cut not passed
	// 3 - curvature calculated but cut not passed
	// 4 - curvature not calculated (division by zero)
	// 5 - helix cut not passed
	// 6 - curvature inxed out of range

	if (curvindex < 0 || curvindex >= fCurvNum) return 6;

	// Find the reference layer
	Int_t lay1 = p1->GetLayer();
	Int_t lay2 = p2->GetLayer();
	Int_t reflayer = (lay1 < lay2) ? lay1 : lay2;
	
	// Swap points in order that r1 < r2
	AliITSNode *temp = 0;
	if (p2->GetLayer() < p1->GetLayer()) {
		temp = p2;
		p2 = p1;
		p1 = temp;
	}

	// shift XY coords to the reference to the vertex position,
	// for easier calculus of quantities.
	Double_t x1 = p1->X() - vx;
	Double_t x2 = p2->X() - vx;
	Double_t y1 = p1->Y() - vy;
	Double_t y2 = p2->Y() - vy;
	Double_t z1 = p1->Z() - vz;
	Double_t z2 = p2->Z() - vz;
	Double_t r1 = sqrt(x1*x1 + y1*y1);
	Double_t r2 = sqrt(x2*x2 + y2*y2);
	
	// Check for theta cut
	Double_t dtheta, dtheta3;
	TVector3 v01(z1, r1, 0.0);
	TVector3 v12(z2 - z1, r2 - r1, 0.0);
	dtheta = v01.Angle(v12) * 180.0 / TMath::Pi();
	TVector3 vv01(x1, y1, z1);
	TVector3 vv12(x2 - x1, y2 - y1, z2 - z1);
	dtheta3 = vv01.Angle(vv12) * 180.0 / TMath::Pi();
	if (dtheta < fThetaCut2DMin[reflayer] || dtheta > fThetaCut2DMax[reflayer]) return 1;
	if (dtheta3 < fThetaCut3DMin[reflayer] || dtheta3 > fThetaCut3DMax[reflayer]) return 2;

	// calculation of curvature
	Double_t dx = x1 - x2, dy = y1 - y2;
	Double_t num = 2 * (x1*y2 - x2*y1);
	Double_t den = r1*r2*sqrt(dx*dx + dy*dy);
	Double_t curv = 0.;
	if (den != 0.) {
		curv = TMath::Abs(num / den);
		if (curv > fCurvCut[curvindex]) return 3;
	}
	else
		return 4;

	// calculation of helix matching
	Double_t arc1 = 2.0 * r1 * curv;
	Double_t arc2 = 2.0 * r2 * curv;
	Double_t helmatch = 0.0;
	if (arc1 > -1.0 && arc1 < 1.0) arc1 = asin(arc1);
	else arc1 = ((arc1 > 0.0) ? 0.5 : 1.5) * TMath::Pi();
	if (arc2 > -1.0 && arc2 < 1.0) arc2 = asin(arc2);
	else arc2 = ((arc2 > 0.0) ? 0.5 : 1.5) * TMath::Pi();
	arc1 /= 2.0 * curv;
	arc2 /= 2.0 * curv;
	if (arc1 == 0.0 || arc2 == 0.0) return kFALSE;
	helmatch = TMath::Abs(z1 / arc1 - z2 / arc2);
	if (helmatch < fHelixMatchCutMin[reflayer] || helmatch > fHelixMatchCutMax[reflayer]) return 5;
	
	return 0;
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::StoreAbsoluteMatches()
{
	// Stores in the 'fMatches' array of each node all the points in the
	// adjacent layers which allow to create neurons accordin to the
	// helix and theta cut, and also to the largest curvature cut

	Int_t ilayer, isector, itheta1, itheta2, check;
	TObjArray *list1 = 0, *list2 = 0;
	AliITSNode *node1 = 0, *node2 = 0;
	Double_t thetamin, thetamax;
	Int_t imin, imax;

	for (isector = 0; isector < fSectorNum; isector++) {
		// sector is chosen once for both lists
		for (ilayer = 0; ilayer < 5; ilayer++) {
			for (itheta1 = 0; itheta1 < 180; itheta1++) {
				list1 = (TObjArray*)fPoints[ilayer][itheta1]->At(isector);
				TObjArrayIter iter1(list1);
				while ( (node1 = (AliITSNode*)iter1.Next()) ) {
					if (node1->GetUser() >= 0) continue;
					node1->Matches()->Clear();
					thetamin = node1->ThetaDeg() - fPolarInterval;
					thetamax = node1->ThetaDeg() + fPolarInterval;
					imin = (Int_t)thetamin;
					imax = (Int_t)thetamax;
					if (imin < 0) imin = 0;
					if (imax > 179) imax = 179;
					for (itheta2 = imin; itheta2 <= imax; itheta2++) {
						list2 = (TObjArray*)fPoints[ilayer + 1][itheta2]->At(isector);
						TObjArrayIter iter2(list2);
						while ( (node2 = (AliITSNode*)iter2.Next()) ) {
							check = PassAllCuts(node1, node2, fCurvNum - 1, fVX, fVY, fVZ);
							switch (check) {
								case 0:
									node1->Matches()->AddLast(node2);
									break;
								case 1:
									//Info("StoreAbsoluteMatches", "Layer %d: THETA 2D cut not passed", ilayer);
									break;
								case 2:
									//Info("StoreAbsoluteMatches", "Layer %d: THETA 3D cut not passed", ilayer);
									break;
								case 3:
									//Info("StoreAbsoluteMatches", "Layer %d: CURVATURE cut not passed", ilayer);
									break;
								case 4:
									//Info("StoreAbsoluteMatches", "Layer %d: Division by ZERO in curvature evaluation", ilayer);
									break;
								case 5:
									//Info("StoreAbsoluteMatches", "Layer %d: HELIX-MATCH cut not passed", ilayer);
									break;
								case 6:
									//Error("PassAllCuts", "Layer %d: Curv index out of range", ilayer);
									break;
								default:
									Warning("StoreAbsoluteMatches", "Layer %d: %d: unrecognized return value", ilayer, check);
							}
						} // while (node2...)
					} // for (itheta2...)
				} // while (node1...)
			} // for (itheta...)
		} // for (ilayer...)
	} // for (isector...)
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::PrintMatches(Bool_t stop)
{
	// Prints the matches for each point
	
	TFile *ft = new TFile("its_findables_v1.root");
	TTree *tt = (TTree*)ft->Get("Tracks");
	Int_t it, nP, nU, lab, nF = (Int_t)tt->GetEntries();
	tt->SetBranchAddress("nhitsP", &nP);
	tt->SetBranchAddress("nhitsU", &nU);
	tt->SetBranchAddress("label", &lab);
	TString strP("|"), strU("|");
	for (it = 0; it < nF; it++) {
		tt->GetEntry(it);
		if (nP >= 5) strP.Append(Form("%d|", lab));
		if (nU >= 5) strU.Append(Form("%d|", lab));
	}

	TObjArray *sector = 0;
	Int_t ilayer, isector, itheta, id[3];
	AliITSNode *node1 = 0, *node2 = 0;

	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (isector = 0; isector < fSectorNum; isector++) {
			for (itheta = 0; itheta < 180; itheta++) {
				sector = (TObjArray*)fPoints[ilayer][itheta]->At(isector);
				TObjArrayIter points(sector);
				while ( (node1 = (AliITSNode*)points.Next()) ) {
					for (it = 0; it < 3; it++) id[it] = node1->GetLabel(it);
					nF = (Int_t)node1->Matches()->GetSize();
					cout << "Node layer: " << node1->GetLayer() << ", labels: ";
					cout << id[0] << " " << id[1] << " " << id[2] << " --> ";
					if (!nF) {
						cout << "NO MatchES!!!" << endl;
						continue;
					}
					cout << nF << " Matches" << endl;
					for (it = 0; it < 3; it++) {
						if (strP.Contains(Form("|%d|", id[it])))
							cout << "Belongs to findable (originary) track " << id[it] << endl;
						if (strU.Contains(Form("|%d|", id[it])))
							cout << "Belongs to findable (post-Kalman) track " << id[it] << endl;
					}
					TObjArrayIter matches(node1->Matches());
					while ( (node2 = (AliITSNode*)matches.Next()) ) {
						cout << "Match with " << node2;
						Int_t *sh = node1->SharedID(node2);
						for (Int_t k = 0; k < 3; k++)
							if (sh[k] > -1) cout << " " << sh[k];
						cout << endl;
      			}
					if (stop) {
						cout << "Press a key" << endl;
						cin.get();
					}
				}
			}
		}
	}
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::ResetNodes(Int_t isector)
{
	// Clears some utility arrays for each node
	
	Int_t ilayer, itheta;
	TObjArray *sector = 0;
	AliITSNode *node = 0;
	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			sector = (TObjArray*)fPoints[ilayer][itheta]->At(isector);
			TObjArrayIter iter(sector);
			for (;;) {
				node = (AliITSNode*)iter.Next();
				if (!node) break;
				node->InnerOf()->Clear();
				node->OuterOf()->Clear();
				delete node->InnerOf();
				delete node->OuterOf();
				node->InnerOf() = new TObjArray;
				node->OuterOf() = new TObjArray;
			}
		}
	}
}
//
//--------------------------------------------------------------------------------------------
//
void AliITSNeuralTracker::NeuralTracking(const char* filesave, TCanvas*& display)
{
// This is the public method that operates the tracking.
// It works sector by sector, and at the end  saves the found tracks.
// Other methods are privare because they have no meaning id used alone,
// and sometimes they could get segmentation faults due to uninitialized
// datamembert they require to work on.
// The argument is a file where the final results have to be stored.

	Bool_t isStable = kFALSE;
	Int_t i, nUnits = 0, nLinks = 0, nTracks = 0, sectTracks = 0, totTracks = 0;

	// tracking through sectors
	cout << endl;
	Int_t sector, curv;
	for (sector = 0; sector < fSectorNum; sector++) {
		cout << "\rSector " << sector << ": " << endl;
		sectTracks = 0;
		for(curv = 0; curv < fCurvNum; curv++) {
			cout << "- curvature step " << curv + 1;
			cout << " (" << fCurvCut[curv] << "): ";
			// units creation
			nUnits = CreateNeurons(sector, curv);
			if (!nUnits) {
				cout << "no neurons --> skipped" << endl;
				continue;
			}
			cout << nUnits << " units, " << flush;
			// units connection
			nLinks = LinkNeurons();
			if (!nLinks) {
				cout << "no connections --> skipped" << endl;
				continue;
			}
			cout << nLinks << " connections, " << flush;
			// neural updating
			for (i = 0;; i++) {
				isStable = Update();
				if (display) Display(display);
				TObjArrayIter iter(fNeurons);
				for (;;) {
					AliITSneuron *n = (AliITSneuron*)iter.Next();
					if (!n) break;
				}
				if (isStable) break;
			}
			cout << i << " cycles --> " << flush;
			// tracks saving
			CleanNetwork();
			nTracks = Save(sector);
			cout << nTracks << " tracks" << endl;
			sectTracks += nTracks;
		}
		totTracks += sectTracks;
		//cout << sectTracks << " tracks found (total = " << totTracks << ")      " << flush;
	}

	cout << endl << totTracks << " tracks found!!!" << endl;
	cout << endl << "Saving results in file " << filesave << "..." << flush;
	TFile *f = new TFile(filesave, "recreate");
	fChains->Write("TreeC");
	f->Close();
}
//
//--------------------------------------------------------------------------------------------
//
Int_t AliITSNeuralTracker::CreateNeurons(Int_t sectoridx, Int_t curvidx)
{
// Fills the neuron arrays, following the cut criteria for the selected step
// (secnum = sector to analyze, curvnum = curvature cut step to use)
// It also sets the initial random activation.
// In order to avoid a large number of 'new' operations, all existing neurons
// are reset and initialized with the new values, and are created new unit only if
// it turns out to be necessary
// the 'flag' argument is used to decide if the lower edge in the intevral
// of the accepted curvatures is given by zero (kFALSE) or by the preceding used cut (kTRUE)
// (if this is the first step, anyway, the minimum is always zero)

	ResetNodes(sectoridx);

	if (fNeurons) delete fNeurons;
	fNeurons = new TObjArray;

	AliITSneuron *unit = 0;
	Int_t itheta, neurons = 0;
	TObjArray *lstsector = 0;
	
	// NEW VERSION
	Double_t vx[6], vy[6], vz[6];
	AliITSNode *p[6] = {0, 0, 0, 0, 0, 0};
	for (itheta = 0; itheta < 180; itheta++) {
		lstsector = (TObjArray*)fPoints[0][itheta]->At(sectoridx);
		TObjArrayIter lay0(lstsector);
		while ( (p[0] = (AliITSNode*)lay0.Next()) ) {
			if (p[0]->GetUser() >= 0) continue;
			vx[0] = fVX;
			vy[0] = fVY;
			vz[0] = fVZ;
			TObjArrayIter lay1(p[0]->Matches());
			while ( (p[1] = (AliITSNode*)lay1.Next()) ) {
				if (p[1]->GetUser() >= 0) continue;
				if (!PassCurvCut(p[0], p[1], curvidx, fVX, fVY, fVZ)) continue;
				unit = new AliITSneuron;
				unit->Inner() = p[0];
				unit->Outer() = p[1];
				unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
				unit->Gain() = new TObjArray;
				fNeurons->AddLast(unit);
				p[0]->InnerOf()->AddLast(unit);
				p[1]->OuterOf()->AddLast(unit);
				neurons++;
				vx[1] = p[0]->X();
				vy[1] = p[0]->Y();
				vz[1] = p[0]->Z();
				TObjArrayIter lay2(p[1]->Matches());
				while ( (p[2] = (AliITSNode*)lay2.Next()) ) {
					if (p[2]->GetUser() >= 0) continue;
					if (!PassCurvCut(p[1], p[2], curvidx, vx[1], vy[1], vz[1])) continue;
					unit = new AliITSneuron;
					unit->Inner() = p[1];
					unit->Outer() = p[2];
					unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
					unit->Gain() = new TObjArray;
					fNeurons->AddLast(unit);
					p[1]->InnerOf()->AddLast(unit);
					p[2]->OuterOf()->AddLast(unit);
					neurons++;
					vx[2] = p[1]->X();
					vy[2] = p[1]->Y();
					vz[2] = p[1]->Z();
					TObjArrayIter lay3(p[2]->Matches());
					while ( (p[3] = (AliITSNode*)lay3.Next()) ) {
						if (p[3]->GetUser() >= 0) continue;
						if (!PassCurvCut(p[2], p[3], curvidx, vx[2], vy[2], vz[2])) continue;
						unit = new AliITSneuron;
						unit->Inner() = p[2];
						unit->Outer() = p[3];
						unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
						unit->Gain() = new TObjArray;
						fNeurons->AddLast(unit);
						p[2]->InnerOf()->AddLast(unit);
						p[3]->OuterOf()->AddLast(unit);
						neurons++;
						vx[3] = p[2]->X();
						vy[3] = p[2]->Y();
						vz[3] = p[2]->Z();
						TObjArrayIter lay4(p[3]->Matches());
						while ( (p[4] = (AliITSNode*)lay4.Next()) ) {
							if (p[4]->GetUser() >= 0) continue;
							if (!PassCurvCut(p[3], p[4], curvidx, vx[3], vy[3], vz[3])) continue;
							unit = new AliITSneuron;
							unit->Inner() = p[3];
							unit->Outer() = p[4];
							unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
							unit->Gain() = new TObjArray;
							fNeurons->AddLast(unit);
							p[3]->InnerOf()->AddLast(unit);
							p[4]->OuterOf()->AddLast(unit);
							neurons++;
							vx[4] = p[3]->X();
							vy[4] = p[3]->Y();
							vz[4] = p[3]->Z();
							TObjArrayIter lay5(p[4]->Matches());
							while ( (p[5] = (AliITSNode*)lay5.Next()) ) {
								if (p[5]->GetUser() >= 0) continue;
								if (!PassCurvCut(p[4], p[5], curvidx, vx[4], vy[4], vz[4])) continue;
								unit = new AliITSneuron;
								unit->Inner() = p[4];
								unit->Outer() = p[5];
								unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
								unit->Gain() = new TObjArray;
								fNeurons->AddLast(unit);
								p[4]->InnerOf()->AddLast(unit);
								p[5]->OuterOf()->AddLast(unit);
								neurons++;
							} // while (p[5])
						} // while (p[4])
					} // while (p[3])
				} // while (p[2])
			} // while (p[1])
		} // while (p[0])
	} // for (itheta...)
	// END OF NEW VERSION

	/* OLD VERSION
	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			lstsector = (TObjArray*)fPoints[ilayer][itheta]->At(sectoridx);
			TObjArrayIter inners(lstsector);
			while ( (inner = (AliITSNode*)inners.Next()) ) {
				if (inner->GetUser() >= 0) continue;
				TObjArrayIter outers(inner->Matches());
				while ( (outer = (AliITSNode*)outers.Next()) ) {
					if (outer->GetUser() >= 0) continue;
					if (!PassCurvCut(inner, outer, curvidx, fVX, fVY, fVZ)) continue;
					unit = new AliITSneuron;
					unit->Inner() = inner;
					unit->Outer() = outer;
					unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
					unit->Gain() = new TObjArray;
					fNeurons->AddLast(unit);
					inner->InnerOf()->AddLast(unit);
					outer->OuterOf()->AddLast(unit);
					neurons++;
				} // for (;;)
			} // for (;;)
		} // for (itheta...)
	} // for (ilayer...)
	*/
	
	fNeurons->SetOwner();
	return neurons;
}
//
//
//
Int_t AliITSNeuralTracker::LinkNeurons() const
{
// Creates the neural synapses among all neurons
// which share a point.

	Int_t total = 0;
	TObjArrayIter iter(fNeurons), *testiter;
	AliITSneuron *neuron = 0, *test = 0;
	for (;;) {
		neuron = (AliITSneuron*)iter.Next();
		if (!neuron) break;
		// gain contributors
		testiter = (TObjArrayIter*)neuron->Inner()->OuterOf()->MakeIterator();
		for (;;) {
			test = (AliITSneuron*)testiter->Next();
			if (!test) break;
			neuron->Add2Gain(test, fGain2CostRatio, fAlignExponent);
			total++;
		}
		delete testiter;
		testiter = (TObjArrayIter*)neuron->Outer()->InnerOf()->MakeIterator();
		for (;;) {
			test = (AliITSneuron*)testiter->Next();
			if (!test) break;
			neuron->Add2Gain(test, fGain2CostRatio, fAlignExponent);
			total++;
		}
		delete testiter;
	}
	return total;
}
//
//
//
Bool_t AliITSNeuralTracker::Update()
{
// A complete updating cycle with the asynchronous method
// when the neural network is stable, kTRUE is returned.

	Double_t actVar = 0.0, totDiff = 0.0;
	TObjArrayIter iter(fNeurons);
	AliITSneuron *unit;
	for (;;) {
		unit = (AliITSneuron*)iter.Next();
		if (!unit) break;
		actVar = Activate(unit);
		// calculation the relative activation variation
		totDiff += actVar;
	}
	totDiff /= fNeurons->GetSize();
	return (totDiff < fStabThreshold);
}
//
//
//
void AliITSNeuralTracker::CleanNetwork()
{
// Removes all deactivated neurons, and resolves the competitions
// to the advantage of the most active unit

	AliITSneuron *unit = 0;
	TObjArrayIter neurons(fNeurons);
	while ( (unit = (AliITSneuron*)neurons.Next()) ) {
		if (unit->Activation() < fActMinimum) {
			unit->Inner()->InnerOf()->Remove(unit);
			unit->Outer()->OuterOf()->Remove(unit);
			delete fNeurons->Remove(unit);
		}
	}
	return;
	Bool_t removed;
	Int_t nIn, nOut;
	AliITSneuron *enemy = 0;
	neurons.Reset();
	while ( (unit = (AliITSneuron*)neurons.Next()) ) {
		nIn = (Int_t)unit->Inner()->InnerOf()->GetSize();
		nOut = (Int_t)unit->Outer()->OuterOf()->GetSize();
		if (nIn < 2 && nOut < 2) continue;
		removed = kFALSE;
		if (nIn > 1) {
			TObjArrayIter competing(unit->Inner()->InnerOf());
			while ( (enemy = (AliITSneuron*)competing.Next()) ) {
				if (unit->Activation() > enemy->Activation()) {
					enemy->Inner()->InnerOf()->Remove(enemy);
					enemy->Outer()->OuterOf()->Remove(enemy);
					delete fNeurons->Remove(enemy);
				}
				else {
					unit->Inner()->InnerOf()->Remove(unit);
					unit->Outer()->OuterOf()->Remove(unit);
					delete fNeurons->Remove(unit);
					removed = kTRUE;
					break;
				}
			}
			if (removed) continue;
		}
		if (nOut > 1) {
			TObjArrayIter competing(unit->Outer()->OuterOf());
			while ( (enemy = (AliITSneuron*)competing.Next()) ) {
				if (unit->Activation() > enemy->Activation()) {
					enemy->Inner()->InnerOf()->Remove(enemy);
					enemy->Outer()->OuterOf()->Remove(enemy);
					delete fNeurons->Remove(enemy);
				}
				else {
					unit->Inner()->InnerOf()->Remove(unit);
					unit->Outer()->OuterOf()->Remove(unit);
					delete fNeurons->Remove(unit);
					removed = kTRUE;
					break;
				}
			}
		}
	}
}
//
//
//
Bool_t AliITSNeuralTracker::CheckOccupation() const
{
// checks if a track recognized has a point for each layer
	Int_t i;
	for (i = 0; i < 6; i++) if (fPoint[i] < 0) return kFALSE;
	return kTRUE;
}
//
//
//
Int_t AliITSNeuralTracker::Save(Int_t sectorid)
{
// Creates chains of active neurons, in order to
// find the tracks obtained as the neural network output.

	// every chain is started from the neurons in the first 2 layers
	Int_t i, check, stored = 0;
	Int_t a = sectorid; a = 0;
	Double_t testact = 0;
	AliITSneuron *unit = 0, *cursor = 0, *fwd = 0;
	AliITSNode *node = 0;
	TObjArrayIter iter(fNeurons), *fwditer;
	TObjArray *removedUnits = new TObjArray(0);
	TObjArray *removedPoints = new TObjArray(6);
	removedUnits->SetOwner(kFALSE);
	removedPoints->SetOwner(kFALSE);
	for (;;) {
		for (i = 0; i < 6; i++) fPoint[i] = -1;
		unit = (AliITSneuron*)iter.Next();
		if (!unit) break;
		if (unit->Inner()->GetLayer() > 0) continue;
		fPoint[unit->Inner()->GetLayer()] = unit->Inner()->PosInTree();
		fPoint[unit->Outer()->GetLayer()] = unit->Outer()->PosInTree();
		node = unit->Outer();
		removedUnits->AddLast(unit);
		removedPoints->AddAt(unit->Inner(), unit->Inner()->GetLayer());
		removedPoints->AddAt(unit->Outer(), unit->Outer()->GetLayer());
		while (node) {
			testact = fActMinimum;
			fwditer = (TObjArrayIter*)node->InnerOf()->MakeIterator();
			fwd = 0;
			for (;;) {
				cursor = (AliITSneuron*)fwditer->Next();
				if (!cursor) break;
				if (cursor->Used()) continue;
				if (cursor->Activation() >= testact) {
					testact = cursor->Activation();
					fwd = cursor;
				}
			}
			if (!fwd) break;
			removedUnits->AddLast(fwd);
			node = fwd->Outer();
			fPoint[node->GetLayer()] = node->PosInTree();
			removedPoints->AddAt(node, node->GetLayer());
		}
		check = 0;
		for (i = 0; i < 6; i++) if (fPoint[i] >= 0) check++;
		if (check >= 6) {
			stored++;
			fChains->Fill();
			for (i = 0; i < 6; i++) {
				node = (AliITSNode*)removedPoints->At(i);
				node->SetUser(1);
			}
			fwditer = (TObjArrayIter*)removedUnits->MakeIterator();
			for (;;) {
				cursor = (AliITSneuron*)fwditer->Next();
				if(!cursor) break;
				cursor->Used() = 1;
			}
		}
	}

	return stored;
}
/*
Int_t AliITSNeuralTracker::Save(Int_t sectoridx)
// Creates chains of active neurons, in order to
// find the tracks obtained as the neural network output.
{
	// Reads the final map of neurons and removes all
	// units with too small activation
	cout << "saving..." << flush;
	
	// every chain is started from the neurons in the first 2 layers
	//                     00111111  00011111  00101111  00110111
	Int_t allowmask[] = {     0x3F,     0x1F,     0x2F,     0x37,
										0x3B,     0x3D,     0x3E };
	//                     00111011  00111101  00111110

	Double_t test_fwd = 0., testback = 0.;
	Int_t ilayer, itheta;
	AliITSNode *node = 0;
	AliITSneuron *fwd = 0, *back = 0;
	TList *listsector = 0;

	cout << "A -" << fActMinimum << "-" << flush;
	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			listsector = (TList*)fPoints[ilayer][itheta]->At(sectoridx);
			TListIter iter(listsector);
			while ( (node = (AliITSNode*)iter.Next()) ) {
				TListIter fwditer(node->InnerOf());
				TListIter backiter(node->OuterOf());
				testfwd = testback = fActMinimum;
				while ( (fwd = (AliITSneuron*)fwditer.Next()) ) {
					if (fwd->Activation() > testfwd) {
						testfwd = fwd->Activation();
						node->fNext = fwd->Outer();
					}
				}
				while ( (back = (AliITSneuron*)backiter.Next()) ) {
					if (back->Activation() > testback) {
						testback = back->Activation();
						node->fPrev = back->Inner();
					}
				}
			}
		}
	}

	cout << "B" << flush;
	for (ilayer = 0; ilayer < 5; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			listsector = (TList*)fPoints[ilayer][itheta]->At(sectoridx);
			TListIter iter(listsector);
			while ( (node = (AliITSNode*)iter.Next()) ) {
				if (node->fNext) {
					if (node->fNext->fPrev != node) node->fNext = 0;
				}
			}
		}
	}

	cout << "C" << flush;
	Bool_t ok;
	Int_t stored = 0, l1, l2, i, mask, im;
	AliITSNode *temp = 0;
	TList *list[2];
	for (itheta = 0; itheta < 180; itheta++) {
		list[0] = (TList*)fPoints[0][itheta]->At(sectoridx);
		list[1] = (TList*)fPoints[1][itheta]->At(sectoridx);
		for (i = 0; i < 2; i++) {
			TListIter iter(list[i]);
			while ( (node = (AliITSNode*)iter.Next()) ) {
				//cout << endl << node << flush;
				AliITSneuralChain *chain = new AliITSneuralChain;
				for (temp = node; temp; temp = temp->fNext) {
					chain->Insert((AliITSNeuralPoint*)temp);
				}
				ok = kFALSE;
				mask = chain->OccupationMask();
				for (im = 0; im < 7; im++) ok = ok || (mask == allowmask[im]);
				if (!ok) {
					delete chain;
					continue;
				}
				//cout << " lab " << flush;
				chain->AssignLabel();
				//cout << " add " << flush;
				fTracks->AddLast(chain);
				stored++;
				//cout << " " << stored << flush;
			}
		}
	}

	cout << "D" << flush;
	return stored;
}
*/
//
//
//
Double_t AliITSNeuralTracker::Activate(AliITSneuron* &unit)
{
// Computes the new neural activation and
// returns the variation w.r.t the previous one

	if (!unit) return 0.0;
	Double_t sumgain = 0.0, sumcost = 0.0, input, actOld, actNew;

	// sum gain contributions
	TObjArrayIter *iter = (TObjArrayIter*)unit->Gain()->MakeIterator();
	AliITSlink *link;
	for(;;) {
		link = (AliITSlink*)iter->Next();
		if (!link) break;
		sumgain += link->Linked()->Activation() * link->Weight();
	}

	// sum cost contributions
	TObjArrayIter *testiter = (TObjArrayIter*)unit->Inner()->InnerOf()->MakeIterator();
	AliITSneuron *linked = 0;
	for (;;) {
		linked = (AliITSneuron*)testiter->Next();
		if (!linked) break;
		if (linked == unit) continue;
		sumcost += linked->Activation();
	}
	delete testiter;
	testiter = (TObjArrayIter*)unit->Outer()->OuterOf()->MakeIterator();
	for (;;) {
		linked = (AliITSneuron*)testiter->Next();
		if (!linked) break;
		if (linked == unit) continue;
		sumcost += linked->Activation();
	}

	//cout << "gain = " << sumgain << ", cost = " << sumcost << endl;

	input = (sumgain - sumcost) / fTemperature;
	actOld = unit->Activation();
#ifdef NEURAL_LINEAR
	if (input <= -2.0 * fTemperature)
		actNew = 0.0;
	else if (input >= 2.0 * fTemperature)
		actNew = 1.0;
	else
		actNew = input / (4.0 * fTemperature) + 0.5;
#else
	actNew = 1.0 / (1.0 + TMath::Exp(-input));
#endif
	unit->Activation() = actNew;
	return TMath::Abs(actNew - actOld);
}
//
//
//
void AliITSNeuralTracker::AliITSneuron::Add2Gain(AliITSneuron *n, Double_t multconst, Double_t exponent)
{
// Adds a neuron to the collection of sequenced ones of another neuron			

	AliITSlink *link = new AliITSlink;
	link->Linked() = n;
	Double_t weight = Weight(n);
	link->Weight() = multconst * TMath::Power(weight, exponent);
	fGain->AddLast(link);
}
//
//
//
Double_t AliITSNeuralTracker::AliITSneuron::Weight(AliITSneuron *n)
{
// computes synaptic weight
	TVector3 me(Outer()->X() - Inner()->X(), Outer()->Y() - Inner()->Y(), Outer()->Z() - Inner()->Z());
	TVector3 it(n->Outer()->X() - n->Inner()->X(), n->Outer()->Y() - n->Inner()->Y(), n->Outer()->Z() - n->Inner()->Z());

	Double_t angle = me.Angle(it);
	Double_t weight = 1.0 - sin(angle);
	return weight;
}
//
//
//
//


AliITSNeuralTracker::AliITSNode::AliITSNode(const AliITSNode &t): AliITSNeuralPoint((AliITSNeuralPoint&)t),
fPosInTree(t.fPosInTree),
fInnerOf(t.fInnerOf),
fOuterOf(t.fOuterOf),
fMatches(t.fMatches),
fNext(t.fNext),
fPrev(t.fPrev)
{ 
	//copy constructor
}

AliITSNeuralTracker::AliITSNode& AliITSNeuralTracker::AliITSNode::operator=(const AliITSNode& t)
{ 
	// assignment operator
  this->~AliITSNode(); 
  new(this) AliITSNode(t);
  return *this;
}


AliITSNeuralTracker::AliITSneuron::AliITSneuron(const AliITSneuron &t)
  : TObject((TObject&)t),
fUsed(t.fUsed),
fActivation(t.fActivation),
fInner(t.fInner),
fOuter(t.fOuter),
fGain(t.fGain)
{ 
	//copy constructor
}

AliITSNeuralTracker::AliITSneuron& AliITSNeuralTracker::AliITSneuron::operator=(const AliITSneuron& t)
{ 
  //assignment operator
  this->~AliITSneuron(); 
  new(this) AliITSneuron(t);
  return *this;
}


AliITSNeuralTracker::AliITSlink::AliITSlink(const AliITSlink &t)
  : TObject((TObject&)t),
fWeight(t.fWeight), 
fLinked(t.fLinked)
{ 
  //copy constructor
}

AliITSNeuralTracker::AliITSlink& AliITSNeuralTracker::AliITSlink::operator=(const AliITSlink& t)
{ 
  // assignment operator
  this->~AliITSlink(); 
  new(this) AliITSlink(t);
  return *this;
}



