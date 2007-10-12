////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderESDChain - the reader class for the Alice ESD           //
// tailored for the Task framework                                            //
// Reads in AliESDfriend to create shared hit/quality information             //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOEVENTREADERESDCHAIN_H
#define ALIFEMTOEVENTREADERESDCHAIN_H

#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include <list>

class AliFemtoEvent;

class AliFemtoEventReaderESDChain : public AliFemtoEventReader 
{
 public:
  AliFemtoEventReaderESDChain();
  AliFemtoEventReaderESDChain(const AliFemtoEventReaderESDChain& aReader);
  ~AliFemtoEventReaderESDChain();

  AliFemtoEventReaderESDChain& operator=(const AliFemtoEventReaderESDChain& aReader);

  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;

  void SetESDSource(AliESDEvent *aESD);
  //  void SetESDfriendSource(AliESDfriend *aFriend);

 protected:

 private:
  string         fFileName;      //name of current ESD file
  bool           fConstrained;   //flag to set which momentum from ESD file will be use
  int            fNumberofEvent; //number of Events in ESD file
  int            fCurEvent;      //number of current event
  unsigned int   fCurFile;       //number of current file
  AliESDEvent*   fEvent;         //ESD event
  //  AliESDfriend*  fEventFriend;

/*   list<Int_t>  **fSharedList;       //! Table (one list per padrow) of clusters which are shared */
/*   list<Int_t>  **fClusterPerPadrow; //! Table (one list per padrow) of clusters in each padrow */
		
#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESDChain, 1)
#endif

    };
  
#endif


