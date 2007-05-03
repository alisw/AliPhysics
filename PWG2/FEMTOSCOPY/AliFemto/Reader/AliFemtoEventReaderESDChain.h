////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoEventReaderESDChain - the reader class for the Alice ESD         ///
/// tailored for the Analysis/Task framework                                 ///
/// Reads in AliESDfriend to create shared hit/quality information           ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoEventReaderESDChain_hh
#define AliFemtoEventReaderESDChain_hh
#include "Base/AliFemtoEventReader.h"
#include "Infrastructure/AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "AliESD.h"
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

  void SetESDSource(AliESD *);
  void SetESDfriendSource(AliESDfriend *);

 protected:

 private:
  string         fFileName;      //name of current ESD file
  bool           fConstrained;   //flag to set which momentum from ESD file will be use
  int            fNumberofEvent; //number of Events in ESD file
  int            fCurEvent;      //number of current event
  unsigned int   fCurFile;       //number of current file
  AliESD*        fEvent;         //ESD event
  AliESDfriend*  fEventFriend;

  list<Int_t>  **fSharedList;       //! Table (one list per padrow) of clusters which are shared
  list<Int_t>  **fClusterPerPadrow; //! Table (one list per padrow) of clusters in each padrow
		
#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESDChain, 1)
#endif

    };
  
#endif


