#ifndef AliITSneuralTracker_h
#define AliITSneuralTracker_h

#include <TObject.h>

#include "AliITSglobalRecPoint.h"
class TObjArray;

class AliITSneuron;
class AliITSneuralTrack;
//class AliITSglobalRecPoint;

///////////////////////////////////////////////////////////////////////
//
// AliITSneuralTracker: 
//
// neural network MFT algorithm 
// for track finding in ITS stand alone
// according to the Denby-Peterson model with adaptments to the 
// ALICE multiplicity
// 
///////////////////////////////////////////////////////////////////////

class AliITSneuralTracker : public TObject {

public:
	
	// constructors & destructor
	AliITSneuralTracker();
	AliITSneuralTracker(Int_t nsecs, Int_t nc, Double_t *c, Double_t theta);
	virtual ~AliITSneuralTracker();
	
	// file reading and points array population
	Int_t ReadFile(const Text_t *fname, Int_t evnum);
	
	// setters for working parameters (NOT cuts)
	void SetDiff(Double_t a) {fDiff=a;}                 // helix adaptment parameter
	void SetExponent(Double_t a) {fExp = a;}            // exponent which selects the high excitory values
	void SetGainToCostRatio(Double_t a) {fRatio = a;}   // ratio between the gain and cost contributions
	void SetTemperature(Double_t a) {fTemp = a;}        // temperature parameter in activation function
	void SetVariationLimit(Double_t a) {fVar = a;}      // stability threshold
	void SetInitLimits(Double_t a, Double_t b)          // intervals for initial random activations
		{fMin = (a<=b) ? a : b; fMax = (a<=b)? b : a;}
	void SetVertex(Double_t x, Double_t y, Double_t z)  // estimated vertex position
		{ fVPos[0] = x; fVPos[1] = y; fVPos[2] = z;}
	
	// neuron managin methods
	Bool_t   SetEdge(AliITSneuron *n, AliITSglobalRecPoint *p, const char what); // sets the neuron's edges
	Bool_t   CalcParams(AliITSneuron *n);              // calculation of helix parameters for the neuron
	
	// neural tracking (the argument is the ROOT file to store results)
	void Go(const char* filesave, Bool_t flag = kFALSE); 
			
private:
		
	// These methods are private to prevent the users to use them in incorrect sequences
		
	// neuron management methods (not to be used publically)
	Bool_t   ResetNeuron(AliITSneuron *n);             // resets all neuron's fields like in the constructor
	Int_t    CountExits (AliITSneuron *n);             // counter for # of neurons going into n's tail
	Int_t    CountEnters(AliITSneuron *n);             // counter for # of neurons startins from n's head
	Double_t Angle(AliITSneuron *n, AliITSneuron *m);  // angle between two neural segments
	Double_t Weight(AliITSneuron *n, AliITSneuron *m); // synaptic weight between two neural segments
	
	// neural network work-flow
	Int_t  Initialize(Int_t snum, Int_t cnum, Bool_t flag = kFALSE); 
			                                      // resets (or creates) neurons array for a new neural tracking
	Bool_t Update();                            // updates the neurons and checks if stabilization has occurred
	Int_t  Save();                              // saves the found tracks after stabilization
	void   DoSector(Int_t sect, Bool_t flag);   // complete work on a single sector
		
	Int_t     fCurvNum;  //  # of curvature cut steps
	Double_t *fCurvCut;  //! value of all curvature cuts
	
	Int_t     fSectorNum;    // no. of azymuthal sectors
	Double_t  fSectorWidth;  // width of an azymuthal sector
	
	Double_t fVPos[3]; // estimated vertex coords
	
	Double_t fMin;      // min initial random activations
	Double_t fMax;      // max initial random activations
	Double_t fThetaCut; // polar angle cut
	Int_t    fThetaNum; // size of theta sectionement
	Double_t fTemp;     // logistic function parameter
	Double_t fVar;      // stability threshold (for mean rel. activations)
	Double_t fRatio;    // ratio between inhibitory and excitory contributions
	Double_t fExp;      // alignment weight
	Double_t fDiff;     // max allowed difference between TanL exstimations from the two neuron edges
	
	TObjArray ***fPoints[6];  //! Collection of recpoints (sectioned in azym. secs)
	TObjArray   *fNeurons[6]; //! Collection of neurons
	TObjArray   *fTracks;	  //! Collection of tracks
	
	ClassDef(AliITSneuralTracker, 1)
};
	

////////////////////////////////////////////////////////////////////////////////


class AliITSneuron : public TObject {
	
	friend class AliITSneuralTracker;
	
public:
		
	         AliITSneuron();
	virtual ~AliITSneuron();
	
	Bool_t   ContainsID(Int_t ID) 
		{ return (fInner->HasID(ID) && fOuter->HasID(ID)); }
	
private:
		
	Double_t fActivation; // Activation value

	Double_t fCurv;    // curvature [= 1 / R]
	Double_t fCX;      // curvature centre (X) in changed RF
	Double_t fCY;      // curvature centre (X) in changed RF
	Double_t fTanL;    // tan(dip angle) = C / freq
	Double_t fPhase;   // 'phase' parameter
	Double_t fDiff;    // difference between tan(lambda) estimated by the two different points
	Double_t fLength;  // segment length
	
	AliITSneuron *fNext; // best outgoing unit
	AliITSneuron *fPrev; // best incoming unit
	
	AliITSglobalRecPoint *fInner; // inner point 
	AliITSglobalRecPoint *fOuter; // outer point

	ClassDef(AliITSneuron,1) // neural tracker helper class
};

	
#endif
