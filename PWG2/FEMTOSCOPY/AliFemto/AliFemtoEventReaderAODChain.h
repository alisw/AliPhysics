////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderAODChain - the reader class for the Alice AOD in Chain  //
// Reads in AOD information and converts it into internal AliFemtoEvent       //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOEVENTREADERAODCHAIN_H
#define ALIFEMTOEVENTREADERAODCHAIN_H
#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "TChain.h"
#include "TBits.h"
#include "AliAODEvent.h"
#include <list>
#include "AliPWG2AODTrack.h"
#include "AliFemtoEventReaderAOD.h"

class AliFemtoEvent;
class AliFemtoTrack;

class AliFemtoEventReaderAODChain : public AliFemtoEventReaderAOD 
{
 public:
  AliFemtoEventReaderAODChain();
  AliFemtoEventReaderAODChain(const AliFemtoEventReaderAODChain &aReader);
  virtual ~AliFemtoEventReaderAODChain();

  AliFemtoEventReaderAODChain& operator=(const AliFemtoEventReaderAODChain& aReader);

  virtual AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetAODSource(AliAODEvent *aAOD);

 protected:

 private:

#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderAODChain, 1)
#endif

    };
  
#endif


