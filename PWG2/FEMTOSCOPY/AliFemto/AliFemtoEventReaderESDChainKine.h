////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderESDChainKine - the reader class for the Alice ESD and   //
// the model Kinematics information tailored for the Task framework and the   //
// Reads in AliESDfriend to create shared hit/quality information             //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOEVENTREADERESDCHAINKINE_H
#define ALIFEMTOEVENTREADERESDCHAINKINE_H

#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include <TTree.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliStack.h>
#include <list>
#include <AliGenEventHeader.h>

class AliFemtoEvent;

class AliFemtoEventReaderESDChainKine : public AliFemtoEventReader 
{
 public:
  AliFemtoEventReaderESDChainKine();
  AliFemtoEventReaderESDChainKine(const AliFemtoEventReaderESDChainKine& aReader);
  ~AliFemtoEventReaderESDChainKine();

  AliFemtoEventReaderESDChainKine& operator=(const AliFemtoEventReaderESDChainKine& aReader);

  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;
  void SetUseTPCOnly(const bool usetpconly);
  bool GetUseTPCOnly() const;

  void SetESDSource(AliESDEvent *aESD);
  void SetStackSource(AliStack *aStack);
  void SetGenEventHeader(AliGenEventHeader *aGenHeader);

 protected:

 private:
  string         fFileName;      // name of current ESD file
  bool           fConstrained;   // flag to set which momentum from ESD file will be use
  bool           fUseTPCOnly;    // flog to set to read TPC only momentum instead of the full
  int            fNumberofEvent; // number of Events in ESD file
  int            fCurEvent;      // number of current event
  unsigned int   fCurFile;       // number of current file
  AliESDEvent*   fEvent;         // ESD event
  AliStack      *fStack;         // Kinematics stack pointer
  AliGenEventHeader *fGenHeader; // Link to the generator event header

  Float_t GetSigmaToVertex(double *impact, double *covar);

#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESDChainKine, 1)
#endif

    };
  
#endif


