/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoEventReaderKinematicsChain - the reader class for the Alice ESD and     //
// the model Kinematics information tailored for the Task framework and the        //
// Reads in AliESDfriend to create shared hit/quality information                  //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//  	       				                                           //
/////////////////////////////////////////////////////////////////////////////////////


#ifndef ALIFEMTOEVENTREADERKINEMATICSCHAIN_H
#define ALIFEMTOEVENTREADERKINEMATICSCHAIN_H

#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"
#include "AliFemtoV0.h"

#include <string>
#include <vector>
#include <TTree.h>
#include <AliStack.h>
#include <list>
#include <AliGenEventHeader.h>

class AliFemtoEvent;

class AliFemtoEventReaderKinematicsChain : public AliFemtoEventReader 
{
 public:
  enum EventMult {kGlobalCount=0, kVZERO=1};
  typedef enum EventMult EstEventMult;

  AliFemtoEventReaderKinematicsChain();
  AliFemtoEventReaderKinematicsChain(const AliFemtoEventReaderKinematicsChain& aReader);
  ~AliFemtoEventReaderKinematicsChain();

  AliFemtoEventReaderKinematicsChain& operator=(const AliFemtoEventReaderKinematicsChain& aReader);

  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;

  //void SetESDSource(AliESDEvent *aESD);
  void SetStackSource(AliStack *aStack);
  void SetGenEventHeader(AliGenEventHeader *aGenHeader);
  void SetRotateToEventPlane(short dorotate);
  void SetUseMultiplicity(EstEventMult aType);
  void CopyAODtoFemtoV0(TParticle *tv0, AliFemtoV0 *tFemtoV0 );

  void ReadOnlyPrimaries(bool primaries);
  void ReadOnlyPrimariesV0(bool primaries);
  void ReadPrimariesSecWeakMaterial(bool primaries);
  void ReadPrimariesSecWeakMaterialV0(bool primaries);

 protected:

 private:
  string         fFileName;      // name of current ESD file
  bool           fConstrained;   // flag to set which momentum from ESD file will be use
  int            fNumberofEvent; // number of Events in ESD file
  int            fCurEvent;      // number of current event
  unsigned int   fCurFile;       // number of current file
  AliStack       *fStack;         // Kinematics stack pointer
  AliGenEventHeader *fGenHeader; // Link to the generator event header
  EstEventMult   fEstEventMult;  // Type of the event multiplicity estimator
  short          fRotateToEventPlane; // Rotate the event so that event plane is at x=0

  bool           fReadOnlyPrimaries; // read only primaries
  bool           fReadOnlyPrimariesV0; // read only primaries, V0 analysis
  bool           fReadPrimariesSecWeakMaterial; //read only primaries, secondaries from weak decays and secondaries from material
  bool           fReadPrimariesSecWeakMaterialV0; //read only primaries, secondaries from weak decays and secondaries from material, V0 analysis

  Float_t GetSigmaToVertex(double *impact, double *covar);

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventReaderKinematicsChain, 1);
  /// \endcond
#endif

    };
  
#endif


