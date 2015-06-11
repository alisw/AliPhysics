//
//  AliEveDataSource.h
//
//  blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//

#ifndef __AliEveDataSource__
#define __AliEveDataSource__

#include "TFile.h"
#include "TTree.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliAODEvent.h"
#include "TNamed.h"

#include "TQObject.h"

struct AliEveData {
  TFile        *fESDFile;		// ESD file.
  TTree        *fESDTree;		// ESD tree.
  TTree        *fHLTESDTree;	// HLT ESD tree.
  AliESDEvent  *fESD;			// ESDEvent object.
  AliESDfriend *fESDfriend;		// ESDfriend object.
  TFile        *fAODFile;		// AOD file.
  TTree        *fAODTree;		// AOD tree.
  AliAODEvent  *fAOD;			// AODEvent object.

  AliEveData()
      : fESDFile(NULL)
      , fESDTree(NULL)
      , fHLTESDTree(NULL)
      , fESD(NULL)
      , fESDfriend(NULL)
      , fAODFile(NULL)
      , fAODTree(NULL)
      , fAOD(NULL)
  {}
};

class AliEveDataSource : public TNamed
{
public:
    AliEveDataSource(bool storageManager=false);
    ~AliEveDataSource();
    
    virtual void Init();
    virtual void GotoEvent(Int_t event);
    virtual AliEveData* GetData() const {return fCurrentData;}

    void StorageManagerOk();     // *SIGNAL*
    void StorageManagerDown();   // *SIGNAL*
    void EventServerOk();        // *SIGNAL*
    void EventServerDown();      // *SIGNAL*
    
protected:
    AliEveData* fCurrentData;
    
private:
    ClassDef(AliEveDataSource, 0); // Interface for getting all event components in a uniform way.
};

#endif
