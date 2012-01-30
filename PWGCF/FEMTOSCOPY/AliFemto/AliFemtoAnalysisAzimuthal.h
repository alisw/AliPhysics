////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFemtoAnalysisReactionPlane - Femtoscopic analysis which mixes event //
// with respect to the z position of the primary vertex and event total   //
// multiplicity and uses only events in certain reaction plane angle bin  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOANALYSISAZIMUTHAL_H
#define ALIFEMTOANALYSISAZIMUTHAL_H

#include "AliFemtoSimpleAnalysis.h"        // base analysis class

class TVector2;

class AliFemtoAnalysisAzimuthal : public AliFemtoSimpleAnalysis {

public:

  AliFemtoAnalysisAzimuthal(unsigned int binsVertex=10, double minVertex=-100., double maxVertex=+100., unsigned int binsMult=10, double minMult=-1.e9, double maxMult=+1.e9, unsigned short binsRP=10);
  AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& TheOriginalAnalysis);  // copy constructor

  AliFemtoAnalysisAzimuthal& operator=(const AliFemtoAnalysisAzimuthal& aAna);

  virtual void ProcessEvent(const AliFemtoEvent* ProcessThisEvent);
  virtual ~AliFemtoAnalysisAzimuthal();
  virtual unsigned int OverflowVertexZ() const { return fOverFlowVertexZ;}
  virtual unsigned int UnderflowVertexZ() const { return fUnderFlowVertexZ;}
  virtual unsigned int OverflowMult() const { return fOverFlowMult;}
  virtual unsigned int UnderflowMult() const { return fUnderFlowMult;}
  double GetCurrentReactionPlane();
  TVector2 GetQVector(AliFemtoParticleCollection* particlecollection);
  virtual void MakePairs(const char* typeIn, AliFemtoParticleCollection *partCollection1, AliFemtoParticleCollection *partCollection2=0);
  virtual TList* GetOutputList();

 // Get the particle cuts
  virtual AliFemtoParticleCut*   FemtoParticleCut() {return fFemtoParticleCut;}
  virtual AliFemtoParticleCut*   FlowParticleCut() {return fFlowParticleCut;}
 // Set the cuts
  void SetFemtoParticleCut(AliFemtoParticleCut* x) {fFemtoParticleCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
  void SetFlowParticleCut(AliFemtoParticleCut* x) {fFlowParticleCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
  void SetEventCut(AliFemtoEventCut* x) {fEventCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
  void SetPairCut(AliFemtoPairCut* x) {fPairCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}

protected:

  AliFemtoParticleCut*         	fFemtoParticleCut;    //  select particles of type #1 
  AliFemtoParticleCut*  	fFlowParticleCut;   //  select particles of type #2 

  double fVertexZ[2];                 /* min/max z-vertex position allowed to be processed */
  unsigned int fVertexZBins;          /* number of VERTEX mixing bins in z-vertex in EventMixing Buffer */
  unsigned int fOverFlowVertexZ;      /* number of events encountered which had too large z-vertex */
  unsigned int fUnderFlowVertexZ;     /* number of events encountered which had too small z-vertex */
  double fMult[2];                    /* min/max multiplicity allowed for event to be processed */
  unsigned int fMultBins;             /* number of MULTIPLICITY mixing bins in z-vertex in EventMixing Buffer */
  unsigned int fOverFlowMult;         /* number of events encountered which had too large multiplicity */
  unsigned int fUnderFlowMult;        /* number of events encountered which had too small multiplicity */
  unsigned short fRPBins;             // Number of reaction plane angle orientation bins
  double fPsi;                		// Reaction plane angle of the current event

#ifdef __ROOT__
  ClassDef(AliFemtoAnalysisAzimuthal, 0)
#endif
    
};

#endif
