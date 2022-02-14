////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderNanoAODChain - the reader class for the Alice AOD in Chain  //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOEVENTREADERNANOAODCHAIN_H
#define ALIFEMTOEVENTREADERNANOAODCHAIN_H
#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "TChain.h"
#include "TBits.h"
#include "AliAODEvent.h"
#include <list>  
/* #include "AliPWG2AODTrack.h" */
#include "AliFemtoEventReaderNanoAOD.h"

 

class AliFemtoEvent;
class AliFemtoTrack;

class AliFemtoEventReaderNanoAODChain : public AliFemtoEventReaderNanoAOD 
{
 public:
  AliFemtoEventReaderNanoAODChain();
  AliFemtoEventReaderNanoAODChain(const AliFemtoEventReaderNanoAODChain &aReader);
  virtual ~AliFemtoEventReaderNanoAODChain();

  AliFemtoEventReaderNanoAODChain& operator=(const AliFemtoEventReaderNanoAODChain& aReader);

  virtual AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetAODSource(AliAODEvent *aAOD);


 protected:


 private:

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventReaderNanoAODChain, 1);
  /// \endcond
#endif

    };
  
#endif


