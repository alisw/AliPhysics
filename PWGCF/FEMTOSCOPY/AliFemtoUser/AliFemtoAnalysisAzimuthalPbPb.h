////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFemtoAnalysisReactionPlane - Femtoscopic analysis which mixes event //
// with respect to the z position of the primary vertex and event total   //
// multiplicity and uses only events in certain reaction plane angle bin  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * Author: Johanna Gramling, University of Heidelberg, jgramlin@cern.ch
 *         Jorge Mercado, University of Heidelberg, jmercado@cern.ch
 *         Vera Loggins, Wayne State University, veraloggins@wayne.edu
 *
 **************************************************************************/

#ifndef ALIFEMTOANALYSISAZIMUTHALPB_H
#define ALIFEMTOANALYSISAZIMUTHALPB_H

#include "AliFemtoSimpleAnalysis.h"        // base analysis class
#include "TH1.h"
#include "TVector2.h"
#include "AliAODEvent.h"
#include "AliFemtoPicoEventRP.h"
#include "AliFemtoPairCutRadialDistanceLM.h"

class TVector2;
class AliFemtoPicoEventRP;

class AliFemtoAnalysisAzimuthalPbPb : public AliFemtoSimpleAnalysis {

public:

  AliFemtoAnalysisAzimuthalPbPb(unsigned int binsVertex=10, double minVertex=-100., double maxVertex=+100., unsigned int binsMult=10, double minMult=-1.e9, double maxMult=+1.e9, unsigned short binsRP=10);
  AliFemtoAnalysisAzimuthalPbPb(const AliFemtoAnalysisAzimuthalPbPb& TheOriginalAnalysis);  // copy constructor

  AliFemtoAnalysisAzimuthalPbPb& operator=(const AliFemtoAnalysisAzimuthalPbPb& aAna);

//   virtual void FillHBTParticleCollectionRP(AliFemtoParticleCut* partCut, AliFemtoEvent* hbtEvent, AliFemtoPicoEventRP* picoevent);
  virtual void ProcessEvent(const AliFemtoEvent* ProcessThisEvent);
  virtual ~AliFemtoAnalysisAzimuthalPbPb();
  virtual unsigned int OverflowVertexZ() const { return fOverFlowVertexZ;}
  virtual unsigned int UnderflowVertexZ() const { return fUnderFlowVertexZ;}
  virtual unsigned int OverflowMult() const { return fOverFlowMult;}
  virtual unsigned int UnderflowMult() const { return fUnderFlowMult;}
  double GetCurrentReactionPlane();
  TVector2 GetQVector(AliFemtoParticleCollection* particlecollection);
  virtual void MakePairs(const char* typeIn, AliFemtoPicoEventRP *coll1, AliFemtoPicoEventRP *coll2=0);
  virtual TList* GetOutputList();
  virtual void Finish() {;}
	
 // Get the particle cuts
  virtual AliFemtoParticleCut*   FirstParticleCut() {return fFirstParticleCut;}
  virtual AliFemtoParticleCut*   SecondParticleCut() {return fSecondParticleCut;}
 // Set the cuts
  void SetFirstParticleCut(AliFemtoParticleCut* x) {fFirstParticleCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
  void SetSecondParticleCut(AliFemtoParticleCut* x) {fSecondParticleCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
  void SetEventCut(AliFemtoEventCut* x) {fEventCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
  void SetPairCut(AliFemtoPairCut* x) {fPairCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
  void SetPairCutRD(AliFemtoPairCutRadialDistanceLM* x) {fPairCutRD = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
  void SetEPhistname(char* histname);

protected:

  AliFemtoParticleCut*         	fFirstParticleCut;    //  select particles of type #1 
  AliFemtoParticleCut*  	fSecondParticleCut;   //  select particles of type #2 
  AliFemtoPairCutRadialDistanceLM* fPairCutRD;
  AliFemtoPicoEventRP*		fPicoEventRP;

  double fVertexZ[2];                 /* min/max z-vertex position allowed to be processed */
  unsigned int fVertexZBins;          /* number of VERTEX mixing bins in z-vertex in EventMixing Buffer */
  unsigned int fOverFlowVertexZ;      /* number of events encountered which had too large z-vertex */
  unsigned int fUnderFlowVertexZ;     /* number of events encountered which had too small z-vertex */
  double fMult[2];                    /* min/max multiplicity allowed for event to be processed */
  unsigned int fMultBins;             /* number of MULTIPLICITY mixing bins in z-vertex in EventMixing Buffer */
  unsigned int fOverFlowMult;         /* number of events encountered which had too large multiplicity */
  unsigned int fUnderFlowMult;        /* number of events encountered which had too small multiplicity */
  unsigned short fRPBins;             // Number of reaction plane angle orientation bins
  double fRP;                		// Reaction plane angle of the current event
  TH1F* fphidist;			// Phi distribution as control histogram
  TH1F* fpairphi;			// Phi distribution as control histogram
  TH1F* fRPdist;			// RP distribution as control histogram
  TH1F* fsubRPdist;
  TH1F* frealpsi;			// Psi distribution for real pairs as control histogram
  TH1F* fmixedpsi;			// Psi distribution for mixed pairs as control histogram

#ifdef __ROOT__
  ClassDef(AliFemtoAnalysisAzimuthalPbPb, 0)
#endif
    
};

#endif
