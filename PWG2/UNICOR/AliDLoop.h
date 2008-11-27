// Author: Dariusz Miskowiec 2007

//=============================================================================
// simple event loop manager
//=============================================================================

#ifndef ALIDLOOP_H
#define ALIDLOOP_H

#include <TObject.h>

class TTree;
class AliDEvent;

//=============================================================================
class AliDLoop : public TObject {
   
 public:
                                        // constructor
  AliDLoop(TTree *tr, AliDEvent *ev0, char *outfil="result.root"); 
  AliDLoop(const AliDLoop &loop);             // copy constructor
  AliDLoop &operator=(const AliDLoop &loop);  // substitution operator
  virtual ~AliDLoop() {}                   // destructor
  void Run(int n=-1) const;             // process n events; n=-1 means all

 protected:
  TTree *fTree0;                        //! input event tree
  TTree *fTree1;                        //! clone of fTree0 for event mixing
  AliDEvent *fEv0;                         //! event
  AliDEvent *fEv1;                         //! clone of fEv0 for event mixing
  TString fOutputFilename;              //! output filename
  Int_t Mem() const;                    // virtual memory in MB

  ClassDef(AliDLoop,1)
};
//=============================================================================

#endif
