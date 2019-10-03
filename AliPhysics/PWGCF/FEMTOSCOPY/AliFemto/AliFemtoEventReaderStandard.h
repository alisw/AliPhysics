////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderStandard - the reader class for the Alice ESD, AOD      //
// the model Kinematics information tailored for the Task framework           //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOEVENTREADERSTANDARD_H
#define ALIFEMTOEVENTREADERSTANDARD_H

#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include <TTree.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDfriend.h>
#include <AliStack.h>
#include <list>
#include <AliGenEventHeader.h>
#include <AliPhysicsSelection.h>
#include <AliESDtrackCuts.h>
#include <AliAODMCParticle.h>

class AliFemtoEvent;
class AliFemtoTrack;

class AliFemtoEventReaderStandard : public AliFemtoEventReader 
{
 public:
  enum InputType {kUnknown = -1, kESD=0, kAOD=1, kESDKine=2, kAODKine=3};
  typedef enum InputType AliFemtoInputType;

  AliFemtoEventReaderStandard();
  AliFemtoEventReaderStandard(const AliFemtoEventReaderStandard& aReader);
  ~AliFemtoEventReaderStandard();

  AliFemtoEventReaderStandard& operator=(const AliFemtoEventReaderStandard& aReader);

  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();

  void SetESDSource(AliESDEvent *aESD);
  void SetAODSource(AliAODEvent *aAOD);
  void SetStackSource(AliStack *aStack);
  void SetInputType(AliFemtoInputType aInput);
  void SetGenEventHeader(AliGenEventHeader *aGenHeader);
  void SetUsePhysicsSelection(const bool usephysics);

  void SetESDTrackCuts(AliESDtrackCuts *esdcuts);
  void SetUseTPCOnly(const bool usetpconly);

 protected:

  void CopyAODtoFemtoTrack(const AliAODTrack *tAodTrack, AliFemtoTrack *tFemtoTrack);
  AliAODMCParticle* GetParticleWithLabel(TClonesArray *mcP, Int_t aLabel);

 private:
  string             fFileName;      //  name of current ESD file
  int                fNumberofEvent; //  number of Events in ESD file
  int                fCurEvent;      //  number of current event
  unsigned int       fCurFile;       //  number of current file
  AliESDEvent       *fESDEvent;      //! ESD event
  AliAODEvent       *fAODEvent;      //! AOD event
  AliStack          *fStack;         //! Kinematics stack pointer
  AliGenEventHeader *fGenHeader;     //! Link to the generator event header
  AliFemtoInputType  fInputType;     //  Select data input type

  bool                 fUsePhysicsSel; // if true the physics selection class will be used
  AliPhysicsSelection *fSelect;        // Class to select only physics events

  AliESDtrackCuts     *fTrackCuts;     // Link to external ESD track cut
  bool                 fUseTPCOnly;    // if true the TPC only parameters will be used

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventReaderStandard, 1);
  /// \endcond
#endif

};
  
#endif


