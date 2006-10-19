// ---------------------------------------------------------------------------------
// AliITStrackerANN
// ---------------------------------------------------------------------------------
// Neural Network-based algorithm for track reconstruction in the ALICE ITS.
// The class contains all that is necessary in order to perform the tracking
// using the ITS clusters in the V2 format.
// The class is organized with some embedded objects which serve for 
// data management, and all setters and getters for all parameters.
// The usage of the class from a tracking macro should be based essentially
// in the initialization (constructor) and the call of a method which
// performs in the right sequence all the operations.
// Temporarily (and maybe definitively) all the intermediate steps are 
// public functions which could be called independently, but they do not
// produce any result if all the preventively required steps have not been 
// performed.
// ---------------------------------------------------------------------------------
// Author: Alberto Pulvirenti (university of Catania)
// Email : alberto.pulvirenti@ct.infn.it
// ---------------------------------------------------------------------------------

#include <TString.h>
#include <TObjArray.h>
#include <TVector3.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMatrixD.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,2)
  #include <TMatrixDEigen.h>
#endif

#include "AliITSgeom.h"
#include "AliITStrackSA.h"
#include "AliITSRecPoint.h"
#include "AliITStrackV2.h"

#include "AliITStrackerANN.h"

const Double_t AliITStrackerANN::fgkPi     = 3.141592653; // pi
const Double_t AliITStrackerANN::fgkHalfPi = 1.570796327; // pi / 2
const Double_t AliITStrackerANN::fgkTwoPi  = 6.283185307; // 2 * pi

ClassImp(AliITStrackerANN)


AliITStrackerANN::AliITStrackerANN() : AliITStrackerV2(), 
fVertexX(0.0),
fVertexY(0.0),
fVertexZ(0.0),
fSectorNum(0),
fSectorWidth(0),
fPolarInterval(0),
fCurvNum(0),
fCurvCut(0),
fActMinimum(0),
fEdge1(0),
fEdge2(0),
fStabThreshold(0),
fTemperature(0),
fGain2CostRatio(0),
fExponent(0),
fNLayers(0),
fFirstModInLayer(0),
fIndexMap(0),
fFoundTracks(0),
fNodes(0),
fNeurons(0),
fMsgLevel(0),
fStructureOK(0),
fGeom(0){
  //Default constructor
}
//__________________________________________________________________________________
AliITStrackerANN::AliITStrackerANN(const AliITSgeom *geom, Int_t msglev) 
: AliITStrackerV2(geom), 
fVertexX(0.0),
fVertexY(0.0),
fVertexZ(0.0),
fSectorNum(0),
fSectorWidth(0),
fPolarInterval(0),
fCurvNum(0),
fCurvCut(0),
fActMinimum(0),
fEdge1(0),
fEdge2(0),
fStabThreshold(0),
fTemperature(0),
fGain2CostRatio(0),
fExponent(0),
fNLayers(0),
fFirstModInLayer(0),
fIndexMap(0),
fFoundTracks(0),
fNodes(0),
fNeurons(0),
fMsgLevel(msglev),
fStructureOK(0),
fGeom(0)
{
/**************************************************************************

		CONSTRUCTOR (standard-to-use version)
		
		Arguments:
			1) the ITS geometry used in the generated event 
			2) the flag for log-messages writing
			
		The AliITSgeometry is used along the class, 
		in order to translate the local AliITSRecPoint coordinates
		into the Global reference frame, which is necessary for the
		Neural Network algorithm to operate.
		In case of serialized use, the log messages should be excluded, 
		in order to save real execution time.
		
		Operations:
			- all pointer data members are initialized 
			  (if possible, otherwise are set to NULL)
			- all numeric data members are initialized to
			  values which allow the Neural Network to operate
			  even if nothing is externally set.
			  
		 NOTE: it is possible that tracking an event 
		       with these default values results in a non-sense.

 **************************************************************************/

	Int_t i;
	
	// Get ITS geometry
	fGeom = (AliITSgeom*)geom;
	
	// Initialize the array of first module indexes per layer
	fNLayers = geom->GetNlayers();
	fFirstModInLayer = new Int_t[fNLayers + 1];
	for (i = 0; i < fNLayers; i++) {
		fFirstModInLayer[i] = fGeom->GetModuleIndex(i + 1, 1, 1);
	}
	fFirstModInLayer[fNLayers] = geom->GetIndexMax();
	
	// initialization: no curvature cut steps
	fCurvNum = 0;
	fCurvCut = 0;

	// initialization: 4 sectors (one for each quadrant)
	fSectorNum   = 4;
	fSectorWidth = fgkHalfPi;
	
	// initialization: theta offset of 20 degrees
	fPolarInterval = 20.0;

	// initialization: array structure not defined
	fStructureOK = kFALSE;

	// initialization: vertex in the origin
	fVertexX = 0.0;
	fVertexY = 0.0;
	fVertexZ = 0.0;
	
	// initialization: uninitialized point array
	fNodes = 0;

	// initialization: very large (unuseful) cut values
	Int_t ilayer;
	for (ilayer = 0; ilayer < 6; ilayer++) {
		fThetaCut2D[ilayer] = TMath::Pi();
		fThetaCut3D[ilayer] = TMath::Pi();
		fHelixMatchCut[ilayer] = 1.0;
	}

	// initialization: inictial activation range between 0.3 and 0.7
	fEdge1 = 0.3;
	fEdge2 = 0.7;

	// initialization: neural network operation & weights
	fTemperature = 1.0;
	fStabThreshold = 0.001;
	fGain2CostRatio = 1.0;
	fExponent = 1.0;
	fActMinimum = 0.5;

	// initialization: uninitialized array of neurons
	fNeurons = 0;
	
	// initialization: uninitialized array of tracks
	fFoundTracks = 0;
}


AliITStrackerANN::AliITStrackerANN(const AliITStrackerANN &n) : AliITStrackerV2((AliITStrackerV2&)n),  
fVertexX(n.fVertexX),
fVertexY(n.fVertexY),
fVertexZ(n.fVertexZ),
fSectorNum(n.fSectorNum),
fSectorWidth(n.fSectorWidth),
fPolarInterval(n.fPolarInterval),
fCurvNum(n.fCurvNum),
fCurvCut(n.fCurvCut),
fActMinimum(n.fActMinimum),
fEdge1(n.fEdge1),
fEdge2(n.fEdge2),
fStabThreshold(n.fStabThreshold),
fTemperature(n.fTemperature),
fGain2CostRatio(n.fGain2CostRatio),
fExponent(n.fExponent),
fNLayers(n.fNLayers),
fFirstModInLayer(n.fFirstModInLayer),
fIndexMap(n.fIndexMap),
fFoundTracks(n.fFoundTracks),
fNodes(n.fNodes),
fNeurons(n.fNeurons),
fMsgLevel(n.fMsgLevel),
fStructureOK(n.fStructureOK),
fGeom(n.fGeom){
  //Copy constructor
}

AliITStrackerANN& AliITStrackerANN::operator=(const AliITStrackerANN& arg){
  //Assignment operator
  this->~AliITStrackerANN();
  new(this) AliITStrackerANN(arg);
  return *this;
}
	
//__________________________________________________________________________________
void AliITStrackerANN::SetCuts
(Int_t ncurv, Double_t *curv, Double_t *theta2D, Double_t *theta3D, Double_t *helix)

/**************************************************************************
 
	CUT SETTER
	
	Allows for the definition of all kind of geometric cuts
	which have been studied in for the creation of a neuron
	from a pair of clusters C1 and C2 on consecutive layers.
	Neuron will be created only if the pair passes ALL of these cuts.
	
	At the moment, we define 4 kinds of geometrical cuts:
		a) cut on the difference of the polar 'theta' angle;
		b) cut on the angle between origin->C1 and C1->C2 in space;
		c) cut on the curvature of the circle passing 
		   through C1, C2 and the primary vertex;
		d) cut on heli-matching of the same three points.
	
	Arguments:
		1) the number of curvature cut steps
		2) the array of curvature cuts for each step 
		   (its dimension is given by the argument 1)
		3) array of 5 cut values (one for each consecutive lauer pair)
		   related to cut a)
		4) array of 5 cut values (one for each consecutive lauer pair)
		   related to cut b)
		5) array of 5 cut values (one for each consecutive lauer pair)
		     related to cut c)
			  
	Operations:
		- gets the values for each cut and stores them in data members
		- in the case of curvature cuts, the cut array 
		  (whose size is not fixed) is allocated
	  
	NOTE: in case that the user wants to set onyl ONE of the 4 cuts array,
	      he can simply pass NULL arguments for other cuts, and (eventually)
			a ZERO as the first argument (if curvature cuts have not to be set).
	
	Anyway, all the cuts have to be set at least once.
	
 **************************************************************************/
{
	// counter
	Int_t i;
	
	/*** Curvature cut setting ***/
	
	// first of all, the curvature cuts are sorted in increasing order
	// (from the smallest to the largest one)
	Int_t *ind = new Int_t[ncurv];
	TMath::Sort(ncurv, curv, ind, kFALSE);
	// then, the curvature cut array is allocated and filled
	// (a message with the list of defined cuts can be optionally shown)
	fCurvCut = new Double_t[ncurv];
	if (fMsgLevel >= 1) cout << "Number of curvature cuts: " << ncurv << endl;
	for (i = 0; i < ncurv; i++) {
		fCurvCut[i] = curv[ind[i]];
		if (fMsgLevel >= 1) cout << " - " << fCurvCut[i] << endl;
	}
	fCurvNum = ncurv;
	
	/*** 'Fixed' cuts setting ***/
	
	// checks what cuts have to be set
	Bool_t doTheta2D = (theta2D != 0);
	Bool_t doTheta3D = (theta3D != 0);
	Bool_t doHelix = (helix != 0);
	// sets the cuts for all layer pairs
	for (i = 0; i < fNLayers - 1; i++) {
		if (doTheta2D) fThetaCut2D[i] = theta2D[i];
		if (doTheta3D) fThetaCut3D[i] = theta3D[i];
		if (doHelix) fHelixMatchCut[i] = helix[i];
	}
	// if required, lists the cuts
	if (!fMsgLevel < 2) return;
	cout << "Theta 2D cuts: ";
	if (!doTheta2D) {
		cout << "<not set>" << endl;
	}
	else {
		cout << endl;
		for (i = 0; i < fNLayers - 1; i++) {
			cout << "For layers " << i << " --> " << i + 1;
			cout << " cut = " << fThetaCut2D[i] << endl;
		}
	}
	cout << "---" << endl;
	cout << "Theta 3D cuts: ";
	if (!doTheta3D) {
		cout << "<not set>" << endl;
	}
	else {
		cout << endl;
		for (i = 0; i < fNLayers - 1; i++) {
			cout << "For layers " << i << " --> " << i + 1;
			cout << " cut = " << fThetaCut3D[i] << endl;
		}
	}
	cout << "---" << endl;
	cout << "Helix-match cuts: ";
	if (!doHelix) {
		cout << "<not set>" << endl;
	}
	else {
		cout << endl;
		for (i = 0; i < fNLayers - 1; i++) {
			cout << "For layers " << i << " --> " << i + 1;
			cout << " cut = " << fHelixMatchCut[i] << endl;
		}
	}
	cout << "---" << endl;
}

//__________________________________________________________________________________ 
Bool_t AliITStrackerANN::GetGlobalXYZ
(Int_t refIndex, 
 Double_t &x, Double_t &y, Double_t &z, Double_t &ex2, Double_t &ey2, Double_t &ez2)
 
/**************************************************************************

	LOCAL TO GLOBAL TRANSLATOR

	Taking information from the ITS geometry stored in the class,
	gets a stored AliITScluster and calculates it coordinates
	and errors in the global reference frame.
	These values are stored in the variables, 
	which are passed by reference.
		
	Arguments:
		1) reference index for the cluster to use
		2) (by reference) place to store the global X-coord into
		3) (by reference) place to store the global Y-coord into
		4) (by reference) place to store the global Z-coord into
		5) (by reference) place to store the global X-coord error into
		6) (by reference) place to store the global Y-coord error into
		7) (by reference) place to store the global Z-coord error into
		
	Operations:
		essentially, determines the ITS module index from the 
		detector index of the AliITSRecPoint object, and extracts
		the roto-translation from the ITS geometry, to convert
		the local module coordinates into the global ones.
		
	Return value:
		- kFALSE if the given cluster index points to a non-existing
		cluster, or if the layer number makes no sense (< 0 or > 6).
		- otherwise, kTRUE (meaning a successful operation).

 **************************************************************************/ 
{
	// checks if the layer is correct
	Int_t ilayer = (refIndex & 0xf0000000) >> 28;
	if (ilayer < 0 || ilayer >= fNLayers) {
		Error("GetGlobalXYZ", "Wrong layer number: %d [range: %d - %d]", ilayer, 0, fNLayers);
		return kFALSE;
	}
	// checks if the referenced cluster exists and corresponds to the passed reference
	AliITSRecPoint *refCluster = (AliITSRecPoint*) GetCluster(refIndex);
	if (!refCluster) {
		Error("GetGlobalXYZ", "Cluster not found for index %d", refIndex);
		return kFALSE;
	}
	
	// determine the detector number
	Int_t detID = refCluster->GetDetectorIndex() + fFirstModInLayer[ilayer];
	
	// get rotation matrix
	Double_t rot[9];
	fGeom->GetRotMatrix(detID, rot);
	
	// get translation vector
	Float_t tx,ty,tz;
	fGeom->GetTrans(detID, tx, ty, tz);
	
	// determine r and phi for the reference conversion
	Double_t r = -(Double_t)tx * rot[1] + (Double_t)ty * rot[0];
	if (ilayer == 0) r = -r;
	Double_t phi = TMath::ATan2(rot[1],rot[0]);
	if (ilayer == 0) phi -= fgkPi;
	
	// sets values for X, Y, Z in global coordinates and their errors
	Double_t cosPhi = TMath::Cos(phi);
	Double_t sinPhi = TMath::Sin(phi);
	x =  r*cosPhi + refCluster->GetY()*sinPhi;
	y = -r*sinPhi + refCluster->GetY()*cosPhi;
	z =  refCluster->GetZ();
	ex2 = refCluster->GetSigmaY2()*sinPhi*sinPhi;
	ey2 = refCluster->GetSigmaY2()*cosPhi*cosPhi;
	ez2 = refCluster->GetSigmaZ2();
	
	return kTRUE;
}

//__________________________________________________________________________________ 
AliITStrackerANN::AliITSnode* AliITStrackerANN::AddNode(Int_t refIndex)

/**************************************************************************

	GENERATOR OF NEURAL NODES

	Fills the array of neural 'nodes', which are the ITS clusters
	translated in the global reference frame.
	Given that the global coordinates are used many times, they are 
	stored in a well-defined structure, in the form of an embedded class.
	Moreover, this class allows a faster navigation among points
	and neurons, by means of some object arrays, storing only the 
	neurons which start from, or end to, the given node.
	Finally, each node contains all the other nodes which match it
	with respect to the fixed walues, in order to perform a faster
	neuron-creation phase.
	
	Arguments:
		1) reference index of the correlated AliITSRecPoint object
		
	Operations:
		- allocates the new AliITSnode objects
		- initializes its object arrays
		- from the global coordinates, calculates the
		  'phi' and 'theta' coordinates, in order to store it
		  into the correct theta-slot and azimutal sector.
		  
	REturn values:
		- the pointer of the creater AliITSnode object
		- in case of errors, a waring is given and a NULL is returned

 **************************************************************************/
{
	// create object and set the reference
	AliITSnode *node = new AliITSnode;
	if (!node) {
		Warning("AddNode", "Error occurred when allocating AliITSnode");
		return 0;
	}
	node->ClusterRef() = refIndex;
	
	// calls the conversion function, which makes also some checks 
	// (layer number within range, existence of referenced cluster)
	if ( !GetGlobalXYZ (
				refIndex, 
				node->X(), node->Y(), node->Z(), 
				node->ErrX2(), node->ErrY2(), node->ErrZ2()
			) ) {return 0;}
	
	// initializes the object arrays
	node->Matches() = new TObjArray;
	node->InnerOf() = new TObjArray;
	node->OuterOf() = new TObjArray;
	
	// finds azimutal and polar sector (in degrees)
	Double_t phi = node->GetPhi() * 180.0 / fgkPi;
	Double_t theta = node->GetTheta() * 180.0 / fgkPi;
	Int_t isector = (Int_t)(phi / fSectorWidth);
	Int_t itheta = (Int_t)theta;
	Int_t ilayer = (refIndex & 0xf0000000) >> 28;
	
	// selects the right TObjArray to store object into
	TObjArray *sector = (TObjArray*)fNodes[ilayer][itheta]->At(isector);
	sector->AddLast(node);
	
	return node;
}

//__________________________________________________________________________________
void AliITStrackerANN::CreateArrayStructure(Int_t nsectors)

/**************************************************************************

	ARRAY STRUCTURE CREATOR
	
	Creates a structure of nested TObjArray's where the AliITSnode's
	have to be stored:
	- the first level is made by 6 arrays (one for each layer)
	- the second level is made by 180 arrays (one for each int part of theta)
	- the third level is made by a variable number of arrays 
	  (one for each azimutal sector)
	  
	Arguments:
		1) the number of azimutal sectors
		
	Operations:
		- calculates the width of each sector, from the argument
		- allocates and initializes all array levels
		- sets a flag which tells the user if this NECESSARY operation
		  has been performed (it is needed BEFORE performing tracking)

 **************************************************************************/
{
	// Set the number of sectors and their width.
	fSectorNum   = nsectors;
	fSectorWidth = 360.0 / (Double_t)fSectorNum;
	if (fMsgLevel >= 2) {
		cout << fSectorNum << " sectors --> sector width (degrees) = " << fSectorWidth << endl;
	}
		
	// Meaningful indexes
	Int_t ilayer, isector, itheta;
	
	// Mark for the created objects
	TObjArray *sector = 0;
	
	// First index: layer
	fNodes = new TObjArray**[fNLayers];
	for (ilayer = 0; ilayer < fNLayers; ilayer++) {
		fNodes[ilayer] = new TObjArray*[180];
		for (itheta = 0; itheta < 180; itheta++) fNodes[ilayer][itheta] = 0;
		for (itheta = 0; itheta < 180; itheta++) {
			fNodes[ilayer][itheta] = new TObjArray(nsectors);
			for (isector = 0; isector < nsectors; isector++) {
				sector = new TObjArray;
				sector->SetOwner();
				fNodes[ilayer][itheta]->AddAt(sector, isector);
			}
		}
	}

	// Sets a checking flag to TRUE. 
	// This flag is checked before filling up the arrays with the points.
	fStructureOK = kTRUE;
}

//__________________________________________________________________________________
Int_t AliITStrackerANN::ArrangePoints(char *exportFile)

/**************************************************************************

	POINTS LOCATOR
	
	This function assembles the operation from the other above methods, 
	and fills the arrays with the clusters already stored in the 
	layers of the tracker.
	Then, in order to use this method, the user MUSTs call LoadClusters()
	before.
	
	Arguments: 
		1) string for a file name where the global ccordinates
		   of all points can be exported (optional).
			If this file must not be created, simply pass a NULL argument
			
	Operations:
		- for each AliITSRecPoint in each AliITSlayer, a ne AliITSnode
		  is created and stored in the correct location.
		  
	Return values:
		- the number of stored points
		- when errors occur, or no points are found, 0 is returned
		
 **************************************************************************/
{
	// Check if the array structure has been created
	if (!fStructureOK) {
		Error("ArrangePoints", "Structure NOT defined. Call CreateArrayStructure() first");
		return 0;
	}

	// meaningful indexes
	Int_t ientry, ilayer, nentries = 0, index;
	Int_t nPtsLayer = 0;
	
	// if the argument is not NULL, a file is opened
	fstream file(exportFile, ios::out);
	if (!exportFile || file.fail()) {
		file.close();
		exportFile = 0;
	}
	
	// scan all layers for node creation
	for (ilayer = 0; ilayer < fNLayers; ilayer++) {
		nPtsLayer = GetNumberOfClustersLayer(ilayer);
		if (fMsgLevel >= 1) {
			cout << "Layer " << ilayer << " --> " << nPtsLayer << " clusters" << endl;
		}
		for (ientry = 0; ientry < nPtsLayer; ientry++) {
			// calculation of cluster index : (Bit mask LLLLIIIIIIIIIIII)
			// (L = bits used for layer)
			// (I = bits used for position in layer)
			index = ilayer << 28;
			index += ientry;
			// add new AliITSnode object
			AliITSnode *n = AddNode(index);
			if ( (n != NULL) && exportFile ) {
				file << index << ' ' << n->X() << ' ' << n->Y() << ' ' << n->Z() << endl;
			}
		}
		nentries += nPtsLayer;
	}
	
	// a conventional final message is put at the end of file
	if (exportFile) {
		file << "-1 0.0 0.0 0.0" << endl;
		file.close();
	}

	// returns the number of points processed
	return nentries;
}

//__________________________________________________________________________________
void AliITStrackerANN::StoreOverallMatches()

/**************************************************************************

	NODE-MATCH ANALYSIS
	
	Once the nodes have been created, a firs analysis is to check
	what pairs will satisfy at least the 'fixed' cuts (theta, helix-match)
	and the most permissive (= larger) curvature cut.
	All these node pairs are suitable for neuron creation. 
	In fact, when performing a Neural Tracking step, the only further check 
	will be a check against the current curvature step, while the other 
	are always the same.
	For thi purpose, each AliITSnode has a data member, named 'fMatches'
	which contains references to all other AliITSnodes in the successive layer
	that form, with it, a 'good' pair, with respect to the above cited cuts.
	Then, in each step for neuron creation, the possible neurons starting from
	each node will be searched ONLY within the nodes referenced in fMatches.
	This, of course, speeds up a lot the neuron creation procedure, at the 
	cost of some memory occupation, which results not to be critical.
	
	Operations:
		- for each AliITSnode, matches are found according to the criteria
		  expressed above, and stored in the node->fMatches array

 **************************************************************************/
{
	// meaningful counters
	Int_t ilayer, isector, itheta1, itheta2, check;
	TObjArray *list1 = 0, *list2 = 0;
	AliITSnode *node1 = 0, *node2 = 0;
	Double_t thetaMin, thetaMax;
	Int_t imin, imax;

	// Scan for each sector
	for (isector = 0; isector < fSectorNum; isector++) {
		// sector is chosen once for both lists
		for (ilayer = 0; ilayer < fNLayers - 1; ilayer++) {
			for (itheta1 = 0; itheta1 < 180; itheta1++) {
				list1 = (TObjArray*)fNodes[ilayer][itheta1]->At(isector);
				TObjArrayIter iter1(list1);
				while ( (node1 = (AliITSnode*)iter1.Next()) ) {
					if (node1->IsUsed()) continue;
					// clear an eventually already present array
					// node1->Matches()->Clear();
					// get the global coordinates and defines the theta interval from cut
					thetaMin = (node1->GetTheta() * 180.0 / fgkPi) - fPolarInterval;
					thetaMax = (node1->GetTheta() * 180.0 / fgkPi) + fPolarInterval;
					imin = (Int_t)thetaMin;
					imax = (Int_t)thetaMax;
					if (imin < 0) imin = 0;
					if (imax > 179) imax = 179;
					// loop on the selected theta slots
					for (itheta2 = imin; itheta2 <= imax; itheta2++) {
						list2 = (TObjArray*)fNodes[ilayer + 1][itheta2]->At(isector);
						TObjArrayIter iter2(list2);
						while ( (node2 = (AliITSnode*)iter2.Next()) ) {
							check = PassAllCuts(node1, node2, fCurvNum - 1, fVertexX, fVertexY, fVertexZ);
							if (check == 0) {
								node1->Matches()->AddLast(node2);
							}
						} // while (node2...)
					} // for (itheta2...)
				} // while (node1...)
			} // for (itheta...)
		} // for (ilayer...)
	} // for (isector...)
}

//__________________________________________________________________________________
Int_t AliITStrackerANN::PassAllCuts
(AliITSnode *inner, AliITSnode *outer, Int_t curvStep,                     
 Double_t vx, Double_t vy, Double_t vz)  
{
// ***********************************************************************************
//
//	This check is called in the above method for finding the matches of each node
//	It check the passed point pair against all the fixed cuts and a specified
//	curvature cut, among all the ones which have been defined.
//	The cuts need a vertex-constraint, which is not absolute, but it is passed
//	as argument.
//	
//	Arguments:
//		1) the point in the inner layer
//		2) the point in the outer layer
//		3) curvature step for the curvature cut check (preferably the last)
//		4) X of the used vertex
//		5) Y of the used vertex
//		6) Z of the used vertex
//		
//	Operations:
//		- if necessary, swaps the two points 
//		  (the first must be in the innermost of the two layers)
//		- checks for the theta cuts
//		- calculates the circle passing through the vertex
//		  and the given points and checks for the curvature cut
//		- using the radius calculated there, checks for the helix-math cut
//	
//	Return values:
//		0 - All cuts passed
//		1 - theta 2D cut not passed
//		2 - theta 3D cut not passed
//		3 - curvature calculated but cut not passed
//		4 - curvature not calculated (division by zero)
//		5 - helix cut not passed
//		6 - curvature inxed out of range
// 
// ***********************************************************************************
	
	// Check for curvature index
	if (curvStep < 0 || curvStep >= fCurvNum) return 6;
	
	// Swap points in order that r1 < r2
	AliITSnode *temp = 0;
	if (outer->GetLayer() < inner->GetLayer()) {
		temp = outer;
		outer = inner;
		inner = temp;
	}
	
	// All cuts are variable according to the layer of the 
	// innermost point (the other point will surely be 
	// in the further one, because we don't check poin pairs
	// which are not in adjacent layers)
	// The reference is given by the innermost point.
	Int_t layRef = inner->GetLayer();
	
	// The calculations in the transverse plane are made in 
	// a shifted reference frame, whose origin corresponds to
	// the reference point passed in the argument.
	Double_t xIn = inner->X() - vx;
	Double_t xOut = outer->X() - vx;
	Double_t yIn = inner->Y() - vy;
	Double_t yOut = outer->Y() - vy;
	Double_t zIn = inner->Z() - vz;
	Double_t zOut = outer->Z() - vz;
	Double_t rIn = TMath::Sqrt(xIn*xIn + yIn*yIn);
	Double_t rOut = TMath::Sqrt(xOut*xOut + yOut*yOut);
	
	// Check for theta cut.
	// There are two different angular cuts:
	// one w.r. to the angle in the 2-dimensional r-z plane...
	Double_t dthetaRZ;
	TVector3 origin2innerRZ(zIn, rIn, 0.0);
	TVector3 inner2outerRZ(zOut - zIn, rOut - rIn, 0.0);
	dthetaRZ = origin2innerRZ.Angle(inner2outerRZ) * 180.0 / fgkPi;
	if (dthetaRZ > fThetaCut2D[layRef]) {
		return 1;
		// r-z theta cut not passed ---> 1
	}
	// ...and another w.r. to the angle in the 3-dimensional x-y-z space
	Double_t dthetaXYZ;
	TVector3 origin2innerXYZ(xIn, yIn, zIn);
	TVector3 inner2outerXYZ(xOut - xIn, yOut - yIn, zOut - zIn);
	dthetaXYZ = origin2innerXYZ.Angle(inner2outerXYZ) * 180.0 / fgkPi;
	if (dthetaXYZ > fThetaCut3D[layRef]) {
		return 2;
		// x-y-z theta cut not passed ---> 2
	}
	
	// Calculation & check of curvature
	Double_t dx = xIn - xOut;
	Double_t dy = yIn - yOut;
	Double_t num = 2.0 * (xIn*yOut - xOut*yIn);
	Double_t den = rIn*rOut*sqrt(dx*dx + dy*dy);
	Double_t curv = 0.;
	if (den != 0.) {
		curv = TMath::Abs(num / den);
		if (curv > fCurvCut[curvStep]) {
			return 3;
			// curvature too large for cut ---> 3
		}
	}
	else {
		Error("PassAllCuts", "Curvature calculations gives zero denominator");
		return 4;
		// error: denominator = 0 ---> 4
	}
		
	// Calculation & check of helix matching
	Double_t helMatch = 0.0;
	Double_t arcIn = 2.0 * rIn * curv;
	Double_t arcOut = 2.0 * rOut * curv;
	if (arcIn > -1.0 && arcIn < 1.0) 
		arcIn = TMath::ASin(arcIn);
	else 
		arcIn = ((arcIn > 0.0) ? 0.5 : 1.5) * TMath::Pi();
	if (arcOut > -1.0 && arcOut < 1.0) 
		arcOut = TMath::ASin(arcOut);
	else 
		arcOut = ((arcOut > 0.0) ? 0.5 : 1.5) * TMath::Pi();
	arcIn /= 2.0 * curv;
	arcOut /= 2.0 * curv;
	if (arcIn == 0.0 || arcOut == 0.0) {
		Error("PassAllCuts", "Calculation returns zero-length arcs: l1=%f, l2=%f", arcIn, arcOut);
		return 4;
		// error: circumference arcs seem to equal zero ---> 4
	}
	helMatch = TMath::Abs(zIn / arcIn - zOut / arcOut);
	if (helMatch > fHelixMatchCut[layRef]) {
		return 5;
		// helix match cut not passed ---> 5
	}
	
	// ALL CUTS PASSED ---> 0
	return 0;
}

//__________________________________________________________________________________
Bool_t AliITStrackerANN::PassCurvCut
(AliITSnode *inner, AliITSnode *outer, 
 Int_t curvStep, 
 Double_t vx, Double_t vy, Double_t vz)
{
//***********************************************************************************
//
//	This method operates essentially like the above one, but it is used
//	during a single step of Neural Tracking, where the curvature cut
//	changes.
//	Then, not necessaryly all the nodes stored in the fMatches array
//	will be suitable for neuron creation in an intermediate step.
//	
//	It has the same arguments of the PassAllCuts() method, but 
//  the theta cut is not checked.
//  Moreover, it has a boolean return value.
//	
//***********************************************************************************
	
	// Check for curvature index
	if (curvStep < 0 || curvStep >= fCurvNum) return 6;
	
	// Find the reference layer
	Int_t layIn = inner->GetLayer();
	Int_t layOut = outer->GetLayer();
	Int_t layRef = (layIn < layOut) ? layIn : layOut;
	
	// The calculations in the transverse plane are made in 
	// a shifted reference frame, whose origin corresponds to
	// the reference point passed in the argument.
	Double_t xIn = inner->X() - vx;
	Double_t xOut = outer->X() - vx;
	Double_t yIn = inner->Y() - vy;
	Double_t yOut = outer->Y() - vy;
	Double_t zIn = inner->Z() - vz;
	Double_t zOut = outer->Z() - vz;
	Double_t rIn = TMath::Sqrt(xIn*xIn + yIn*yIn);
	Double_t rOut = TMath::Sqrt(xOut*xOut + yOut*yOut);
	
	// Calculation & check of curvature
	Double_t dx = xIn - xOut;
	Double_t dy = yIn - yOut;
	Double_t num = 2.0 * (xIn*yOut - xOut*yIn);
	Double_t den = rIn*rOut*sqrt(dx*dx + dy*dy);
	Double_t curv = 0.;
	/* OLD VERSION
	if (den != 0.) {
		curv = TMath::Abs(num / den);
		if (curv > fCurvCut[curvStep]) return kFALSE;
		return kTRUE;
	}
	else
		return kFALSE;
	*/
	// NEW VERSION
	if (den != 0.) {
		curv = TMath::Abs(num / den);
		if (curv > fCurvCut[curvStep]) {
			return kFALSE;
		}
	}
	else {
		Error("PassAllCuts", "Curvature calculations gives zero denominator");
		return kFALSE;
	}
		
	// Calculation & check of helix matching
	Double_t helMatch = 0.0;
	Double_t arcIn = 2.0 * rIn * curv;
	Double_t arcOut = 2.0 * rOut * curv;
	if (arcIn > -1.0 && arcIn < 1.0) 
		arcIn = TMath::ASin(arcIn);
	else 
		arcIn = ((arcIn > 0.0) ? 0.5 : 1.5) * TMath::Pi();
	if (arcOut > -1.0 && arcOut < 1.0) 
		arcOut = TMath::ASin(arcOut);
	else 
		arcOut = ((arcOut > 0.0) ? 0.5 : 1.5) * TMath::Pi();
	arcIn /= 2.0 * curv;
	arcOut /= 2.0 * curv;
	if (arcIn == 0.0 || arcOut == 0.0) {
		Error("PassAllCuts", "Calculation returns zero-length arcs: l1=%f, l2=%f", arcIn, arcOut);
		return 4;
		// error: circumference arcs seem to equal zero ---> 4
	}
	helMatch = TMath::Abs(zIn / arcIn - zOut / arcOut);
	return (helMatch <= fHelixMatchCut[layRef]);
}

//__________________________________________________________________________________
void AliITStrackerANN::PrintMatches(Bool_t stop)
{
//	Prints the list of points which appear to match
//	each one of them, according to the preliminary 
//	overall cuts.
//	The arguments states if a pause is required after printing
//	the matches for each one. In this case, a keypress is needed.

	TObjArray *sector = 0;
	Int_t ilayer, isector, itheta, nF;
	AliITSnode *node1 = 0, *node2 = 0;
	//AliITSRecPoint *cluster1 = 0, *cluster2 = 0;

	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (isector = 0; isector < fSectorNum; isector++) {
			for (itheta = 0; itheta < 180; itheta++) {
				sector = (TObjArray*)fNodes[ilayer][itheta]->At(isector);
				TObjArrayIter points(sector);
				while ( (node1 = (AliITSnode*)points.Next()) ) {
					nF = (Int_t)node1->Matches()->GetEntries();
					cout << "Node layer: " << node1->GetLayer() << " --> ";
					if (!nF) {
						cout << "NO Matches!!!" << endl;
						continue;
					}
					cout << nF << " Matches" << endl;
					cout << "Reference cluster: " << hex << node1->ClusterRef() << endl;
					TObjArrayIter matches(node1->Matches());
					while ( (node2 = (AliITSnode*)matches.Next()) ) {
						cout << "Match with " << hex << node2->ClusterRef() << endl;
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

//__________________________________________________________________________________
void AliITStrackerANN::ResetNodes(Int_t isector)
{
/***********************************************************************************

	NODE NEURON ARRAY CLEANER
	
	After a neural tracking step, this method
	clears the arrays 'fInnerOf' and 'fOuterOf' of each AliITSnode
	
	Arguments:
		- the sector where the operation is being executed
		
 ***********************************************************************************/

	Int_t ilayer, itheta;
	TObjArray *sector = 0;
	AliITSnode *node = 0;
	for (ilayer = 0; ilayer < fNLayers; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			sector = (TObjArray*)fNodes[ilayer][itheta]->At(isector);
			TObjArrayIter iter(sector);
			for (;;) {
				node = (AliITSnode*)iter.Next();
				if (!node) break;
				node->InnerOf()->Clear();
				node->OuterOf()->Clear();
				/*
				delete node->InnerOf();
				delete node->OuterOf();
				node->InnerOf() = new TObjArray;
				node->OuterOf() = new TObjArray;
				*/
			}
		}
	}
}

//__________________________________________________________________________________
Int_t AliITStrackerANN::CreateNeurons
(AliITSnode *node, Int_t curvStep, Double_t vx, Double_t vy, Double_t vz)
{
	// This method is used to create alle suitable neurons starting from
	// a single AliITSnode. Each unit is also stored in the fInnerOf array
	// of the passed node, and in the fOuterOf array of the other neuron edge.
	// In the new implementation of the intermediate check steps, a further one
	// is made, which chechs how well a helix is matched by three points
	// in three consecutive layers.
	// Then, a vertex-constrained check is made with vertex located
	// in a layer L, for points in layers L+1 and L+2.
	// 	
	// In order to do this, the creator works recursively, in a tree-visit like operation.
	// The neurons are effectively created only if the node argument passed is in 
	// the 5th layer (they are created between point of 5th and 6th layer).
	// If the node is in an inner layer, its coordinates are passet as vertex for a nested
	// call of the same function in the next two layers.
	// 
	// Arguments:
	// 	1) reference node
	// 	2) current curvature cut step
	// 	3) X of referenced temporary (not primary) vertex
	// 	4) Y of referenced temporary (not primary) vertex
	// 	5) Z of referenced temporary (not primary) vertex
	// 	
	// Operations:
	// 	- if the layer is the 5th, neurons are created with nodes
	// 	  in the fMatches array of the passed node
	// 	- otherwise, the X, Y, Z of the passed node are given as 
	// 	  vertex and the same procedure is recursively called for all
	// 	  nodes in the fMatches array of the passed one.
	// 	  
	// Return values:
	// 	- the total number of neurons created from the passed one
	// 	  summed with all neurons created from all nodes well matched with it
	// 	  (assumes a meaning only for nodes in the first layer)
	
	// local variables
	Int_t made = 0;         // counter
	Bool_t found = 0;       // flag
	AliITSnode *match = 0;  // cursor for a AliITSnode array
	AliITSneuron *unit = 0; // cursor for a AliITSneuron array
	
	// --> Case 0: the passed node has already been used 
	// as member of a track found in a previous step. 
	// In this case, of course, the function exits.
	if (node->IsUsed()) return 0;
	
	// --> Case 1: there are NO well-matched points.
	// This can happen in all ITS layers, but it happens **for sure**
	// for a node in the 'last' layer.
	// Even in this case, the function exits.
	if (node->Matches()->IsEmpty()) return 0;
	
	// --> Case 2: there are well-matched points.
	// In this case, the function creates a neuron for each
	// well-matched pair (according to the cuts for the current step)
	// Moreover, before storing the neuron, a check is necessary
	// to avoid the duplicate creation of the same neuron twice.
	// (This could happen if the 3 last arguments of the function
	// are close enough to cause a good match for the current step
	// between two points, independently of their difference).
	// Finally, a node is skipped if it has already been used.
	// For each matched point for which a neuron is created, the procedure is 
	// recursively called.
	TObjArrayIter matches(node->Matches());
	while ( (match = (AliITSnode*)matches.Next()) ) {
		if (match->IsUsed()) continue;
		if (!PassCurvCut(node, match, curvStep, vx, vy, vz)) continue;
		found = kFALSE;
		if (!node->InnerOf()->IsEmpty()) {
			TObjArrayIter presentUnits(node->InnerOf());
			while ( (unit = (AliITSneuron*)presentUnits.Next()) ) {
				if (unit->Inner() == node && unit->Outer() == match) {
					found = kTRUE;
					break;
				}
			}
		}
		if (found) continue;
		AliITSneuron *unit = new AliITSneuron(node, match, fEdge2, fEdge1);
		fNeurons->AddLast(unit);
		node->InnerOf()->AddLast(unit);
		match->OuterOf()->AddLast(unit);
		made += CreateNeurons(match, curvStep, node->X(), node->Y(), node->Z());
		made++;
	}
	
	// Of course, the return value contains the number of neurons
	// counting in also the oned created in all levels of recursive calls.
	return made;
}
 
//__________________________________________________________________________________
Int_t AliITStrackerANN::CreateNetwork(Int_t sector, Int_t curvStep)
{
	// This function simply recalls the CreateNeurons() method for each node
	// in the first layer, for the current sector. 
	// This generates the whole network, thanks to the recursive calls.
	//
	// Arguments:
	//	 1) current sector
	//	 2) current curvature step
	//	
	// Operations:
	//	- scans the nodes array for all theta's in the current sector
	//	  and layer 0, and calls the CreateNeurons() function.
	
	// removes all eventually present neurons
	if (fNeurons) delete fNeurons;
	fNeurons = new TObjArray;
	fNeurons->SetOwner(kTRUE);
	
	// calls the ResetNodes() function to free the AliITSnode arrays
	if (fMsgLevel >= 2) {
		cout << "Sector " << sector << " PHI = ";
		cout << fSectorWidth * (Double_t)sector << " --> ";
		cout << fSectorWidth * (Double_t)(sector + 1) << endl;
		cout << "Curvature step " << curvStep << " [cut = " << fCurvCut[curvStep] << "]" << endl;
	}
	ResetNodes(sector);

	// meaningful counters
	Int_t itheta, neurons = 0;
	TObjArray *lstSector = 0;
	
	// NEW VERSION
	Double_t vx[6], vy[6], vz[6];
	AliITSnode *p[6] = {0, 0, 0, 0, 0, 0};
	for (itheta = 0; itheta < 180; itheta++) {
		lstSector = (TObjArray*)fNodes[0][itheta]->At(sector);
		TObjArrayIter lay0(lstSector);
		while ( (p[0] = (AliITSnode*)lay0.Next()) ) {
			if (p[0]->IsUsed()) continue;
			vx[0] = fVertexX;
			vy[0] = fVertexY;
			vz[0] = fVertexZ;
			neurons += CreateNeurons(p[0], curvStep, fVertexX, fVertexY, fVertexZ);
			/*
			TObjArrayIter lay1(p[0]->Matches());
			while ( (p[1] = (AliITSnode*)lay1.Next()) ) {
				if (p[1]->IsUsed()) continue;
				if (!PassCurvCut(p[0], p[1], curvStep, vx[0], vy[0], vz[0])) continue;
				unit = new AliITSneuron;
				unit->Inner() = p[0];
				unit->Outer() = p[1];
				unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
				unit->fGain = new TObjArray;
				fNeurons->AddLast(unit);
				p[0]->InnerOf()->AddLast(unit);
				p[1]->OuterOf()->AddLast(unit);
				neurons++;
				vx[1] = p[0]->X();
				vy[1] = p[0]->Y();
				vz[1] = p[0]->Z();
				TObjArrayIter lay2(p[1]->Matches());
				while ( (p[2] = (AliITSnode*)lay2.Next()) ) {
					if (p[2]->IsUsed()) continue;
					if (!PassCurvCut(p[1], p[2], curvStep, vx[1], vy[1], vz[1])) continue;
					unit = new AliITSneuron;
					unit->Inner() = p[1];
					unit->Outer() = p[2];
					unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
					unit->fGain = new TObjArray;
					fNeurons->AddLast(unit);
					p[1]->InnerOf()->AddLast(unit);
					p[2]->OuterOf()->AddLast(unit);
					neurons++;
					vx[2] = p[1]->X();
					vy[2] = p[1]->Y();
					vz[2] = p[1]->Z();
					TObjArrayIter lay3(p[2]->Matches());
					while ( (p[3] = (AliITSnode*)lay3.Next()) ) {
						if (p[3]->IsUsed()) continue;
						if (!PassCurvCut(p[2], p[3], curvStep, vx[2], vy[2], vz[2])) continue;
						unit = new AliITSneuron;
						unit->Inner() = p[2];
						unit->Outer() = p[3];
						unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
						unit->fGain = new TObjArray;
						fNeurons->AddLast(unit);
						p[2]->InnerOf()->AddLast(unit);
						p[3]->OuterOf()->AddLast(unit);
						neurons++;
						vx[3] = p[2]->X();
						vy[3] = p[2]->Y();
						vz[3] = p[2]->Z();
						TObjArrayIter lay4(p[3]->Matches());
						while ( (p[4] = (AliITSnode*)lay4.Next()) ) {
							if (p[4]->IsUsed()) continue;
							if (!PassCurvCut(p[3], p[4], curvStep, vx[3], vy[3], vz[3])) continue;
							unit = new AliITSneuron;
							unit->Inner() = p[3];
							unit->Outer() = p[4];
							unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
							unit->fGain = new TObjArray;
							fNeurons->AddLast(unit);
							p[3]->InnerOf()->AddLast(unit);
							p[4]->OuterOf()->AddLast(unit);
							neurons++;
							vx[4] = p[3]->X();
							vy[4] = p[3]->Y();
							vz[4] = p[3]->Z();
							TObjArrayIter lay5(p[4]->Matches());
							while ( (p[5] = (AliITSnode*)lay5.Next()) ) {
								if (p[5]->IsUsed()) continue;
								if (!PassCurvCut(p[4], p[5], curvStep, vx[4], vy[4], vz[4])) continue;
								unit = new AliITSneuron;
								unit->Inner() = p[4];
								unit->Outer() = p[5];
								unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
								unit->fGain = new TObjArray;
								fNeurons->AddLast(unit);
								p[4]->InnerOf()->AddLast(unit);
								p[5]->OuterOf()->AddLast(unit);
								neurons++;
							} // while (p[5])
						} // while (p[4])
					} // while (p[3])
				} // while (p[2])
			} // while (p[1])
			*/
		} // while (p[0])
	} // for (itheta...)
	// END OF NEW VERSION

	/* OLD VERSION
	for (ilayer = 0; ilayer < 6; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			lstSector = (TObjArray*)fNodes[ilayer][itheta]->At(sector_idx);
			TObjArrayIter inners(lstSector);
			while ( (inner = (AliITSnode*)inners.Next()) ) {
				if (inner->GetUser() >= 0) continue;
				TObjArrayIter outers(inner->Matches());
				while ( (outer = (AliITSnode*)outers.Next()) ) {
					if (outer->GetUser() >= 0) continue;
					if (!PassCurvCut(inner, outer, curvStep, fVX, fVY, fVZ)) continue;
					unit = new AliITSneuron;
					unit->Inner() = inner;
					unit->Outer() = outer;
					unit->Activation() = gRandom->Rndm() * (fEdge1 - fEdge2) + fEdge2;
					unit->fGain = new TObjArray;
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
 
 //__________________________________________________________________________________
Int_t AliITStrackerANN::LinkNeurons()
{
/***********************************************************************************
	
	SYNAPSIS GENERATOR
	
	Scans the whole neuron array, in order to find all neuron pairs
	which are connected in sequence and share a positive weight.
	For each of them, an AliITSlink is created, which stores 
	the weight value, and will allow for a faster calculation
	of the total neural input for each updating cycle.
	
	Every neuron contains an object array which stores all AliITSlink
	objects which point to sequenced units, with the respective weights.
	
	Return value:
		- the number of link created for the neural network. 
		  If they are 0, no updating can be done and the step is skipped.
	
 ***********************************************************************************/

	// meaningful indexes
	Int_t total = 0;
	Double_t weight = 0.0;
	TObjArrayIter neurons(fNeurons), *iter;
	AliITSneuron *neuron = 0, *test = 0;
	
	// scan in the neuron array
	for (;;) {
		neuron = (AliITSneuron*)neurons.Next();
		if (!neuron) break;
		// checks for neurons startin from its head ( -> )
 		iter = (TObjArrayIter*)neuron->Inner()->OuterOf()->MakeIterator();
		for (;;) {
			test = (AliITSneuron*)iter->Next();
			if (!test) break;
			weight = Weight(test, neuron);
			if (weight > 0.0) neuron->Gain()->AddLast(new AliITSlink(weight, test));
			total++;
		}
		delete iter;
		// checks for neurons ending in its tail ( >- )
		iter = (TObjArrayIter*)neuron->Outer()->InnerOf()->MakeIterator();
		for (;;) {
			test = (AliITSneuron*)iter->Next();
			if (!test) break;
			weight = Weight(neuron, test);
			if (weight > 0.0) neuron->Gain()->AddLast(new AliITSlink(weight, test));
			total++;
		}
		delete iter;
	}
	return total;
}
 
//__________________________________________________________________________________
Bool_t AliITStrackerANN::Update()
{
/***********************************************************************************
	
	Performs a single updating cycle.
	
	Operations:
		- for each neuron, gets the activation with the neuron Activate() method
		- checks if stability has been reached (compare mean activation variation
		  with the stability threshold data member)
	
	Return values:
		- kTRUE means that the neural network has stabilized
		- kFALSE means that another updating cycle is needed
		
 ***********************************************************************************/

	Double_t actVar = 0.0, totDiff = 0.0;
	TObjArrayIter iter(fNeurons);
	AliITSneuron *unit;
	for (;;) {
		unit = (AliITSneuron*)iter.Next();
		if (!unit) break;
		actVar = unit->Activate(fTemperature);
		// calculation the relative activation variation
		totDiff += actVar;
	}
	totDiff /= fNeurons->GetSize();
	return (totDiff < fStabThreshold);
}

//__________________________________________________________________________________
void AliITStrackerANN::FollowChains(Int_t sector)
{
/***********************************************************************************
	
	CHAINS CREATION
	
	After that the neural network has stabilized, 
	the final step is to create polygonal chains 
	of clusters, one in each layer, which represent 
	the tracks recognized by the neural algorithm.
	This is made by means of a choice of the best 
	neuron among the ones starting from each point.
	
	Once that such neuron is selected, its inner point
	will set the 'fNext' field to its outer point, and 
	similarly, its outer point will set the 'fPrev' field
	to its inner point. 
	This defines a bi-directional sequence. 
	
	In this procedure, it can happen that many neurons
	which have the head of the arrow in a given node, will
	all select as best following the neuron with the largest 
	activation starting in that point. 
	This results in MANY nodes which have the same 'fNext'.
	But, this field will be set to NULL for all these points, 
	but the only one which is pointed by the 'fPrev' field
	of this shared node.
 
 ***********************************************************************************/

	// meaningful counters
	Int_t itheta, ilayer;
	TObjArray *lstSector = 0;
	Double_t test = fActMinimum;
	AliITSnode *p = 0;
	AliITSneuron *n = 0;
	
	// scan the whole collection of nodes
	for (ilayer = 0; ilayer < fNLayers; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			// get the array of point in a given layer/theta-slot/sector
			lstSector = (TObjArray*)fNodes[ilayer][itheta]->At(sector);
			TObjArrayIter nodes(lstSector);
			while ( (p = (AliITSnode*)nodes.Next()) ) {
				// if the point is used, it is skipped
				if (p->IsUsed()) continue;
				// initially, fNext points to nothing, and
				// the comparison value is set to the minimum activation
				// which allows to say that a neuron is turned 'on'
				// a node from which only 'off' neurons start is probably
				// a noise point, which will be excluded from the reconstruction.
				test = fActMinimum;
				p->Next() = 0;
				TObjArrayIter innerof(p->InnerOf());
				while ( (n = (AliITSneuron*)innerof.Next()) ) {
					// if the examined neuron has not the largest activation
					// it is skipped and removed from array of all neurons
					// and of its outer point (its inner is the cursor p)
					if (n->Activation() < test) {
						p->InnerOf()->Remove(n);
						n->Outer()->OuterOf()->Remove(n);
						delete fNeurons->Remove(n);
						continue;
					}
					// otherwise, its activation becomes the maximum reference
					p->Next() = n->Outer();
					// at the exit of the while(), the fNext will point
					// to the outer node of the neuron starting in p, whose
					// activation is the largest.
				}
				// the same procedure is made now for all neurons
				// for which p is the outer point
				test = fActMinimum;
				p->Prev() = 0;
				TObjArrayIter outerof(p->OuterOf());
				while ( (n = (AliITSneuron*)outerof.Next()) ) {
					// if the examined neuron has not the largest activation
					// it is skipped and removed from array of all neurons
					// and of its inner point (its outer is the cursor p)
					if (n->Activation() < test) {
						p->OuterOf()->Remove(n);
						n->Inner()->InnerOf()->Remove(n);
						delete fNeurons->Remove(n);
						continue;
					}
					// otherwise, its activation becomes the maximum reference
					p->Prev() = n->Inner();
					// at the exit of the while(), the fPrev will point
					// to the inner node of the neuron ending in p, whose
					// activation is the largest.
				}
			} // end while (p ...)
		} // end for (itheta ...)
	} // end for (ilayer ...)
	
	// now the mismatches are solved
	Bool_t matchPrev, matchNext;
	for (ilayer = 0; ilayer < fNLayers; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			// get the array of point in a given layer/theta-slot/sector
			lstSector = (TObjArray*)fNodes[ilayer][itheta]->At(sector);
			TObjArrayIter nodes(lstSector);
			while ( (p = (AliITSnode*)nodes.Next()) ) {
				// now p will point to a fPrev and a fNext node.
				// Ideally they are placed this way: fPrev --> P --> fNext
				// A mismatch happens if the point addressed as fPrev does NOT
				// point to p as its fNext. And the same for the point addressed
				// as fNext.
				// In this case, the fNext and fPrev pointers are set to NULL
				// and p is excluded from the reconstruction
				matchPrev = matchNext= kFALSE;
				if (ilayer > 0 && p->Prev() != NULL) 
					if (p->Prev()->Next() == p) matchPrev = kTRUE;
				if (ilayer < 5 && p->Next() != NULL) 
					if (p->Next()->Prev() == p) matchNext = kTRUE;
				if (ilayer == 0) 
					matchPrev = kTRUE;
				else if (ilayer == 5)
					matchNext = kTRUE;
				if (!matchNext || !matchPrev) {
					p->Prev() = p->Next() = 0;
				}
			} // end while (p ...)
		} // end for (itheta ...)
	} // end for (ilayer ...)
}

//__________________________________________________________________________________
Int_t AliITStrackerANN::SaveTracks(Int_t sector)
{
/********************************************************************************

	TRACK SAVING
	------------	
	Using the fNext and fPrev pointers, the chain is followed 
	and the track is fitted and saved.
	Of course, the track is followed as a chain with a point
	for each layer, then the track following starts always
	from the clusters in layer 0.
	
***********************************************************************************/
 
	// if not initialized, the tracks TobjArray is initialized
	if (!fFoundTracks) fFoundTracks = new TObjArray;

	// meaningful counters
	Int_t itheta, ilayer, l;
	TObjArray *lstSector = 0;
	AliITSnode *p = 0, *q = 0, **node = new AliITSnode*[fNLayers];
	for (l = 0; l < fNLayers; l++) node[l] = 0;
	
	/*
	array = new AliITSnode*[fNLayers + 1];
	for (l = 0; l <= fNLayers; l++) array[l] = 0;
	array[0] = new AliITSnode();
	array[0]->X() = fVertexX;
	array[0]->Y() = fVertexY;
	array[0]->Z() = fVertexZ;
	array[0]->ErrX2() = fVertexErrorX;
	array[0]->ErrY2() = fVertexErrorY;
	array[0]->ErrZ2() = fVertexErrorZ;
	*/
	Double_t *param = new Double_t[8];
	
	// scan the whole collection of nodes
	for (ilayer = 0; ilayer < 1; ilayer++) {
		for (itheta = 0; itheta < 180; itheta++) {
			// get the array of point in a given layer/theta-slot/sector
			lstSector = (TObjArray*)fNodes[ilayer][itheta]->At(sector);
			TObjArrayIter nodes(lstSector);
			while ( (p = (AliITSnode*)nodes.Next()) ) {
				for (q = p; q; q = q->Next()) {
					l = q->GetLayer();
					node[l] = q;
				}
				//if (!RiemannFit(fNLayers, node, param)) continue;
				// initialization of Kalman Filter Tracking
				AliITSRecPoint *cluster = (AliITSRecPoint*)GetCluster(node[0]->ClusterRef());
				Int_t mod = cluster->GetDetectorIndex();
				Int_t lay, lad, det;
				fGeom->GetModuleId(mod, lay, lad, det);
				Float_t y0 = cluster->GetY();
				Float_t z0 = cluster->GetZ();
				AliITStrackSA* trac = new AliITStrackSA(fGeom,lay, lad, det, 
				                                        y0, z0, 
																	 param[4], param[7], param[3], 1);
				for (l = 0; l < fNLayers; l++) {
					cluster = (AliITSRecPoint*)GetCluster(node[l]->ClusterRef());
					if (cluster) trac->AddClusterV2(l, (node[l]->ClusterRef() & 0x0fffffff)>>0);
				}
            AliITStrackV2* ot = new AliITStrackV2(*trac);
				ot->ResetCovariance(10.);
				ot->ResetClusters();
				if (RefitAt(49.,ot,trac)) { //fit from layer 1 to layer 6
					AliITStrackV2 *otrack2 = new AliITStrackV2(*ot);
					otrack2->ResetCovariance(10.);
					otrack2->ResetClusters();
					//fit from layer 6 to layer 1
					if (RefitAt(3.7,otrack2,ot)) fFoundTracks->AddLast(otrack2);
				}       
				// end of Kalman Filter fit
			}
		}
	}
	
	return 1;
}
 
//__________________________________________________________________________________
void AliITStrackerANN::ExportTracks(const char *filename) const
{
// Exports found tracks into a TTree of AliITStrackV2 objects
	TFile *file = new TFile(filename, "RECREATE");
	TTree *tree = new TTree("TreeT-ANN", "Tracks found in ITS stand-alone with Neural Tracking");
	AliITStrackV2 *track = 0;
	tree->Branch("Tracks", &track, "AliITStrackV2");
	TObjArrayIter tracks(fFoundTracks);
	while ( (track = (AliITStrackV2*)tracks.Next()) ) {
		tree->Fill();
	}
	file->cd();
	tree->Write();
	file->Close();
}


//__________________________________________________________________________________
void AliITStrackerANN::CleanNetwork()
{
	// Removes deactivated units from the network

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
 
//__________________________________________________________________________________
Int_t AliITStrackerANN::StoreTracks()
{
	// Stores the tracks found in a single neural tracking step.
	// In order to do this, it sects each neuron which has a point
	// in the first layer.
	// Then 
	
	// if not initialized, the tracks TobjArray is initialized
	if (!fFoundTracks) fFoundTracks = new TObjArray;
	
	Int_t i, check, stored = 0;
	Double_t testAct = 0;
	AliITSneuron *unit = 0, *cursor = 0, *fwd = 0;
	AliITSnode *node = 0;
	TObjArrayIter iter(fNeurons), *fwdIter;
	TObjArray *removedUnits = new TObjArray(0);
	removedUnits->SetOwner(kFALSE);
	AliITStrackANN annTrack(fNLayers);
	
	for (;;) {
		unit = (AliITSneuron*)iter.Next();
		if (!unit) break;
		if (unit->Inner()->GetLayer() > 0) continue;
		annTrack.SetNode(unit->Inner()->GetLayer(), unit->Inner());
		annTrack.SetNode(unit->Outer()->GetLayer(), unit->Outer());
		node = unit->Outer();
		removedUnits->AddLast(unit);
		while (node) {
			testAct = fActMinimum;
			fwdIter = (TObjArrayIter*)node->InnerOf()->MakeIterator();
			fwd = 0;
			for (;;) {
				cursor = (AliITSneuron*)fwdIter->Next();
				if (!cursor) break;
				if (cursor->Used()) continue;
				if (cursor->Activation() >= testAct) {
					testAct = cursor->Activation();
					fwd = cursor;
				}
			}
			if (!fwd) break;
			removedUnits->AddLast(fwd);
			node = fwd->Outer();
			annTrack.SetNode(node->GetLayer(), node);
		}
		check = annTrack.CheckOccupation();
		if (check >= 6) {
			stored++;
			// FIT
			//if (!RiemannFit(fNLayers, trackitem, param)) continue;
			if (!annTrack.RiemannFit()) continue;
			// initialization of Kalman Filter Tracking
			AliITSRecPoint *cluster = (AliITSRecPoint*)GetCluster(annTrack[0]->ClusterRef());
			Int_t mod = cluster->GetDetectorIndex();
			Int_t lay, lad, det;
			fGeom->GetModuleId(mod, lay, lad, det);
			Float_t y0 = cluster->GetY();
			Float_t z0 = cluster->GetZ();
			AliITStrackSA* trac = new AliITStrackSA(fGeom,lay, lad, det, y0, z0, 
																 annTrack.Phi(), annTrack.TanLambda(), 
																 annTrack.Curv(), 1);
			for (Int_t l = 0; l < fNLayers; l++) {
				if (!annTrack[l]) continue;
				cluster = (AliITSRecPoint*)GetCluster(annTrack[l]->ClusterRef());
				if (cluster) trac->AddClusterV2(l, (annTrack[l]->ClusterRef() & 0x0fffffff)>>0);
			}
			AliITStrackV2* ot = new AliITStrackV2(*trac);
			ot->ResetCovariance(10.);
			ot->ResetClusters();
			if (RefitAt(49.,ot,trac)) { //fit from layer 1 to layer 6
				AliITStrackV2 *otrack2 = new AliITStrackV2(*ot);
				otrack2->ResetCovariance(10.);
				otrack2->ResetClusters();
				//fit from layer 6 to layer 1
				if (RefitAt(3.7,otrack2,ot)) fFoundTracks->AddLast(otrack2);
			}       
			// end of Kalman Filter fit
			// END FIT
			for (i = 0; i < fNLayers; i++) {
				//node = (AliITSnode*)removedPoints->At(i);
				//node->Use();
				annTrack[i]->Use();
			}
			fwdIter = (TObjArrayIter*)removedUnits->MakeIterator();
			for (;;) {
				cursor = (AliITSneuron*)fwdIter->Next();
				if(!cursor) break;
				cursor->Used() = 1;
			}
		}
	}

	return stored;
}

Double_t AliITStrackerANN::Weight(AliITSneuron *nAB, AliITSneuron *nBC)
{
/***********************************************************************************
 * Calculation of neural weight.
 * The implementation of positive neural weight is set only in the case
 * of connected units (e.g.: A->B with B->C).
 * Given that B is the **common** point. care should be taken to pass 
 * as FIRST argument the neuron going "to" B, and
 * as SECOND argument the neuron starting "from" B
 * anyway, a check is put in order to return 0.0 when arguments are not well given.
 ***********************************************************************************/
 
	if (nAB->Outer() != nBC->Inner()) {
		if (nBC->Outer() == nAB->Inner()) {
			AliITSneuron *temp = nAB;
			nAB = nBC;
			nBC = temp;
			temp = 0;
			if (fMsgLevel >= 3) {
				Info("Weight", "Switching wrongly ordered arguments.");
			}
		}
		Warning("Weight", "Not connected segments. Returning 0.0");
		return 0.0;
	}
	
	AliITSnode *pA = nAB->Inner();
	AliITSnode *pB = nAB->Outer();
	AliITSnode *pC = nBC->Outer();
	
	TVector3 vAB(pB->X() - pA->X(), pB->Y() - pA->Y(), pB->Z() - pA->Z());
	TVector3 vBC(pC->X() - pB->X(), pC->Y() - pB->Y(), pC->Z() - pB->Z());

	Double_t weight = 1.0 - sin(vAB.Angle(vBC));
	return fGain2CostRatio * TMath::Power(weight, fExponent);
}



/******************************************
 ******************************************
 *** AliITStrackerANN::AliITSnode class ***
 ******************************************
 ******************************************/
 
//__________________________________________________________________________________ 
AliITStrackerANN::AliITSnode::AliITSnode(): 
fX(0.0),
fY(0.0),
fZ(0.0),
fEX2(0.0),
fEY2(0.0),
fEZ2(0.0),
fUsed(kFALSE), 
fClusterRef(-1), 
fMatches(NULL), 
fInnerOf(NULL), 
fOuterOf(NULL), 
fNext(NULL), 
fPrev(NULL)
{
	// Constructor for the embedded 'AliITSnode' class.
	// It initializes all pointer-like objects.
	
}

AliITStrackerANN::AliITSnode::AliITSnode(const AliITSnode &n) : TObject((TObject&)n),
fX(n.fX),
fY(n.fY),
fZ(n.fZ),
fEX2(n.fEX2),
fEY2(n.fEY2),
fEZ2(n.fEZ2),
fUsed(n.fUsed),
fClusterRef(n.fClusterRef),
fMatches(n.fMatches),
fInnerOf(n.fInnerOf),
fOuterOf(n.fOuterOf),
fNext(n.fNext),
fPrev(n.fPrev){
  //copy constructor
} 
//__________________________________________________________________________________ 
AliITStrackerANN::AliITSnode::~AliITSnode()
{
	// Destructor for the embedded 'AliITSnode' class.
	// It should clear the object arrays, but it is possible
	// that some objects still are useful after the point deletion
	// then the arrays are cleared but their objects are owed by
	// another TCollection object, and not deleted.
	// For safety reasons, all the pointers are set to zero.
	
	fInnerOf->SetOwner(kFALSE); 
	fInnerOf->Clear(); 
	delete fInnerOf;
	fInnerOf = 0;
	fOuterOf->SetOwner(kFALSE); 
	fOuterOf->Clear(); 
	delete fOuterOf;
	fOuterOf = 0;
	fMatches->SetOwner(kFALSE); 
	fMatches->Clear();
	delete fMatches;
	fMatches = 0;
	fNext = 0;
	fPrev = 0;
}

//__________________________________________________________________________________ 
Double_t AliITStrackerANN::AliITSnode::GetPhi() const
{
	// Calculates the 'phi' (azimutal) angle, and returns it
	// in the range between 0 and 2Pi radians.
	
	Double_t q;
	q = TMath::ATan2(fY,fX); 
	if (q >= 0.) 
		return q;
	else 
		return q + TMath::TwoPi();
}
 
//__________________________________________________________________________________ 
Double_t AliITStrackerANN::AliITSnode::GetError(Option_t *option)
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
	// reference to the inline ErrX2(), ErrY2() adn ErrZ2() methods.
	
	TString opt(option);
	Double_t errorSq = 0.0;
	opt.ToUpper();

	if (opt.Contains("R2")) {
		errorSq  = fX*fX*fEX2 + fY*fY*fEY2;
		errorSq /= GetR2sq();
	}
	else if (opt.Contains("R3")) {
		errorSq  = fX*fX*fEX2 + fY*fY*fEY2 + fZ*fZ*fEZ2;
		errorSq /= GetR3sq();
	}
	else if (opt.Contains("PHI")) {
		errorSq  = fY*fY*fEX2;
		errorSq += fX*fX*fEY2;
		errorSq /= GetR2sq() * GetR2sq();
	}
	else if (opt.Contains("THETA")) {
		errorSq = fZ*fZ * (fX*fX*fEX2 + fY*fY*fEY2);
		errorSq += GetR2sq() * GetR2sq() * fEZ2;
		errorSq /= GetR3sq() * GetR3sq() * GetR2() * GetR2();
	}
	
	if (!opt.Contains("SQ")) 
		return TMath::Sqrt(errorSq);
	else 
		return errorSq;
}



/********************************************
 ********************************************
 *** AliITStrackerANN::AliITSneuron class ***
 ********************************************
 ********************************************/
 
//__________________________________________________________________________________
AliITStrackerANN::AliITSneuron::AliITSneuron
(AliITSnode *inner, AliITSnode *outer, Double_t minAct, Double_t maxAct): 
fUsed(0), 
fActivation(0),
fInner(inner), 
fOuter(outer),
fGain(0) 
{
	// Default neuron constructor
	fActivation = gRandom->Rndm() * (maxAct-minAct) + minAct;
	fGain = new TObjArray;
}

AliITStrackerANN::AliITSneuron::AliITSneuron(const AliITSneuron &n) : TObject((TObject&)n),
fUsed(n.fUsed),
fActivation(n.fActivation),
fInner(n.fInner),
fOuter(n.fOuter),
fGain(n.fGain){
  //copy constructor
} 
//__________________________________________________________________________________
Double_t AliITStrackerANN::AliITSneuron::Activate(Double_t temperature)
{
	// This computes the new activation of a neuron, and returns
	// its activation variation as a consequence of the updating.
	// 
	// Arguments:
	// - the 'temperature' parameter for the neural activation logistic function
	// 
	// Operations:
	// - collects the total gain, by summing the products
	//   of the activation of each sequenced unit by the relative weight.
	// - collects the total cost, by summing the activations of 
	//   all competing units
	// - passes the sum of gain - cost to the activation function and
	//   calculates the new activation
	//   
	// Return value:
	// - the difference between the old activation and the new one
	//   (absolute value)

	// local variables
	Double_t sumGain = 0.0;      // total contribution from chained neurons
	Double_t sumCost = 0.0;      // total contribution from crossing neurons
	Double_t input;              // total input
	Double_t actOld, actNew;     // old and new values for the activation
	AliITSneuron *linked = 0;    // cursor for scanning the neuron arrays (for link check)
	AliITSlink *link;            // cursor for scanning the synapses arrays (for link check)
	TObjArrayIter *iterator = 0; // pointer to the iterator along the neuron arrays

	// sum contributions from the correlated units
	iterator = (TObjArrayIter*)fGain->MakeIterator();
	for(;;) {
		link = (AliITSlink*)iterator->Next();
		if (!link) break;
		sumGain += link->Contribution();
	}
	delete iterator;

	// sum contributions from the competing units:
	// the ones which have the same starting point...
	iterator = (TObjArrayIter*)fInner->InnerOf()->MakeIterator();
	for (;;) {
		linked = (AliITSneuron*)iterator->Next();
		if (!linked) break;
		if (linked == this) continue;
		sumCost += linked->fActivation;
	}
	delete iterator;
	// ...and the ones which have the same ending point
	iterator = (TObjArrayIter*)fOuter->OuterOf()->MakeIterator();
	for (;;) {
		linked = (AliITSneuron*)iterator->Next();
		if (!linked) break;
		if (linked == this) continue;
		sumCost += linked->fActivation;
	}

	// calculate the total input as the difference between gain and cost
	input = (sumGain - sumCost) / temperature;
	actOld = fActivation;
	// calculate the final output
#ifdef NEURAL_LINEAR
	if (input <= -2.0 * temperature)
		actNew = 0.0;
	else if (input >= 2.0 * temperature)
		actNew = 1.0;
	else
		actNew = input / (4.0 * temperature) + 0.5;
#else
	actNew = 1.0 / (1.0 + TMath::Exp(-input));
#endif
	fActivation = actNew;
	
	// return the activation variation
	return TMath::Abs(actNew - actOld);
}



/******************************************
 ******************************************
 *** AliITStrackerANN::AliITSlink class ***
 ******************************************
 ******************************************/
 
 // No methods defined non-inline
 
 
 
 /**********************************************
  **********************************************
  *** AliITStrackerANN::AliITStrackANN class ***
  **********************************************
  **********************************************/
 
//__________________________________________________________________________________
AliITStrackerANN::AliITStrackANN::AliITStrackANN(Int_t dim) : 
fNPoints(dim),
fXCenter(0.0),
fYCenter(0.0),
fRadius(0.0),
fCurv(0.0),
fDTrans(0.0),
fDLong(0.0),
fTanLambda(0.0),
fPhi(0.0),
fNode(0)
{
	// Default constructor for the AliITStrackANN class
	
	if (! dim) {
		fNode = 0;
	}
	else{
		Int_t i = 0;
		fNode = new AliITSnode*[dim];
		for (i = 0; i < dim; i++) fNode[i] = 0;
	}
}

AliITStrackerANN::AliITStrackANN::AliITStrackANN(const AliITStrackANN &n) : TObject((TObject&)n),
fNPoints(n.fNPoints),
fXCenter(n.fXCenter),
fYCenter(n.fYCenter),
fRadius(n.fRadius),
fCurv(n.fCurv),
fDTrans(n.fDTrans),
fDLong(n.fDLong),
fTanLambda(n.fTanLambda),
fPhi(n.fPhi),
fNode(n.fNode)
{
  //copy constructor
}
//__________________________________________________________________________________
Int_t AliITStrackerANN::AliITStrackANN::CheckOccupation() const
{
	// Returns the number of pointers fNode which are not NULL	

	Int_t i;         // cursor
	Int_t count = 0; // counter for how many points are stored in the track
	
	for (i = 0; i < fNPoints; i++) {
		if (fNode[i] != NULL) count++;
	}
	
	return count;
}

//__________________________________________________________________________________
Bool_t AliITStrackerANN::AliITStrackANN::RiemannFit()
{ 
	// Performs the Riemann Sphere fit for the given points to a circle
	// and then uses the fit parameters to fit a helix in space.
	//
	// Return values:
	// - kTRUE if all operations have been performed
	// - kFALSE if the numbers risk to lead to an arithmetic violation

	Int_t i, j, count, dim = fNPoints;
	
	// First check for all points
	count = CheckOccupation();
	if (count != fNPoints) {
		Error ("AliITStrackANN::RiemannFit", "CheckOccupations returns %d, fNPoints = %d ==> MISMATCH", count, fNPoints);
		return kFALSE;
	}

	// matrix of ones
	TMatrixD m1(dim,1);
	for (i = 0; i < dim; i++) m1(i,0) = 1.0;

	// matrix of Rieman projection coordinates
	TMatrixD coords(dim,3);
	for (i = 0; i < dim; i++) {
		coords(i,0) = fNode[i]->X();
		coords(i,1) = fNode[i]->Y();
		coords(i,2) = fNode[i]->GetR2sq();
	}

	// matrix of weights
	Double_t xterm, yterm, ex, ey;
	TMatrixD weights(dim,dim);
	for (i = 0; i < dim; i++) {
		xterm = fNode[i]->X() * fNode[i]->GetPhi() - fNode[i]->Y() / fNode[i]->GetR2();
		ex = fNode[i]->ErrX2();
		yterm = fNode[i]->Y() * fNode[i]->GetPhi() + fNode[i]->X() / fNode[i]->GetR2();
		ey = fNode[i]->ErrY2();
		weights(i,i) = fNode[i]->GetR2sq() / (xterm * xterm * ex + yterm * yterm * ey );
	}

	// weighted sample mean
	Double_t meanX = 0.0, meanY = 0.0, meanW = 0.0, sw = 0.0;
	for (i = 0; i < dim; i++) {
		meanX += weights(i,i) * coords(i,0);
		meanY += weights(i,i) * coords(i,1);
		meanW += weights(i,i) * coords(i,2);
		sw += weights(i,i);
	}
	meanX /= sw;
	meanY /= sw;
	meanW /= sw;

	// sample covariance matrix
	for (i = 0; i < dim; i++) {
		coords(i,0) -= meanX;
		coords(i,1) -= meanY;
		coords(i,2) -= meanW;
	}
	TMatrixD coordsT(TMatrixD::kTransposed, coords);
	TMatrixD weights4coords(weights, TMatrixD::kMult, coords);
	TMatrixD sampleCov(coordsT, TMatrixD::kMult, weights4coords);
	for (i = 0; i < 3; i++) {
		for (j = i + 1; j < 3; j++) {
			sampleCov(i,j)  = sampleCov(j,i)  = (sampleCov(i,j) + sampleCov(j,i)) * 0.5;
		}
	}

	// Eigenvalue problem solving for V matrix
	Int_t ileast = 0;
	TVectorD eval(3), n(3);
#if ROOT_VERSION_CODE < ROOT_VERSION(4,0,2)
	TMatrixD evec = sampleCov.EigenVectors(eval);
#else
	TMatrixDEigen ei(sampleCov);
	TMatrixD evec = ei.GetEigenVectors();
	eval = ei.GetEigenValuesRe();
#endif
	if (eval(1) < eval(ileast)) ileast = 1;
	if (eval(2) < eval(ileast)) ileast = 2;
	n(0) = evec(0, ileast);
	n(1) = evec(1, ileast);
	n(2) = evec(2, ileast);

	// c - known term in the plane intersection with Riemann axes
	Double_t c = -(meanX * n(0) + meanY * n(1) + meanW * n(2));

	// center and radius of fitted circle
	Double_t xc, yc, radius, curv;
	xc = -n(0) / (2. * n(2));
	yc = -n(1) / (2. * n(2));
	radius = (1. - n(2)*n(2) - 4.*c*n(2)) / (4. * n(2) * n(2));
	
	if (radius <= 0.E0) {
		Error("RiemannFit", "Radius = %f less than zero!!!", radius);
		return kFALSE;
	}
	radius = TMath::Sqrt(radius);
	curv = 1.0 / radius;
		
	// evaluating signs for curvature and others
	Double_t phi1 = 0.0, phi2, temp1, temp2, phi0, sumdphi = 0.0;
	AliITSnode *p = fNode[0];
	phi1 = p->GetPhi();
	for (i = 1; i < dim; i++) {
		p = (AliITSnode*)fNode[i];
		if (!p) break;
		phi2 = p->GetPhi();
		temp1 = phi1;
		temp2 = phi2;
		if (temp1 > fgkPi && temp2 < fgkPi)
			temp2 += fgkTwoPi;
		else if (temp1 < fgkPi && temp2 > fgkPi)
			temp1 += fgkTwoPi;
		sumdphi += temp2 - temp1;
		phi1 = phi2;
	}
	if (sumdphi < 0.E0) curv = -curv;
	Double_t diff, angle = TMath::ATan2(yc, xc);
	if (curv < 0.E0)
		phi0 = angle + 0.5 * TMath::Pi();
	else
		phi0 = angle - 0.5 * TMath::Pi();
	diff = angle - phi0;

	Double_t dt, temp = TMath::Sqrt(xc*xc + yc*yc) - radius;
	if (curv >= 0.E0)
		dt = temp;
	else
		dt = -temp;
	//cout << "Dt = " << dt << endl;
	
	Double_t halfC = 0.5 * curv, test;
	Double_t *s = new Double_t[dim], *zz = new Double_t[dim], *ws = new Double_t[dim];
	for (j = 0; j < 6; j++) {
		p = fNode[j];
		if (!p) break;
		//----
		s[j] = (p->GetR2sq() - dt * dt) / (1. + curv * dt);
		if (s[j] < 0.) {
			if (TMath::Abs(s[j]) < 1.E-6) s[j] = 0.;
			else {
				Error("RiemannFit", "Square root argument error: %17.15g < 0", s[j]);
				return kFALSE;
			}
		}
		s[j] = TMath::Sqrt(s[j]);
		//cout << "Curv = " << halfC << " --- s[" << j << "] = " << s[j] << endl;
		s[j] *= halfC;
		test = TMath::Abs(s[j]);
		if (test > 1.) {
			if (test <= 1.1) 
				s[j] = ((s[j] > 0.) ? 0.99999999999 : -0.9999999999);
			else {
				Error("RiemannFit", "Value too large: %17.15g", s[j]);
				return kFALSE;
			}
		}
		//----
		zz[j] = p->Z();
		s[j] = TMath::ASin(s[j]) / halfC;
		ws[j] = 1.0 / (p->ErrZ2());
	}

	// second tep final fit
	Double_t s2Sum = 0.0, zSum = 0.0, szSum = 0.0, sSum = 0.0, sumw = 0.0;
	for (i = 0; i < dim; i++) {
		s2Sum += ws[i] * s[i] * s[i];
		zSum  += ws[i] * zz[i];
		sSum  += ws[i] * s[i];
		szSum += ws[i] * s[i] * zz[i];
		sumw += ws[i];
	}
	s2Sum /= sumw;
	zSum /= sumw;
	sSum /= sumw;
	szSum /= sumw;
	temp = s2Sum - sSum*sSum;

	Double_t dz, tanL;
	dz = (s2Sum*zSum - sSum*szSum) / temp;
	tanL = (szSum - sSum*zSum) / temp;
	
	fXCenter = xc;
	fYCenter = yc;
	fRadius = radius;
	fCurv = curv;
	fPhi = phi0;
	fDTrans = dt;
	fDLong = dz;
	fTanLambda = tanL;

	delete [] s;
	delete [] zz;
	delete [] ws;
	
	return kTRUE;
}
