// $Header$

// Tools for import of kinematics. 
// Preliminary/minimal solution.

#ifndef ALIEVE_KineTools_H
#define ALIEVE_KineTools_H

#include <TEveUtil.h>
#include <TObject.h>

class TTree;
class AliStack;

class TEveTrackList;

namespace Alieve {

class KineTools
{
private:
  KineTools(const KineTools&);            // Not implemented
  KineTools& operator=(const KineTools&); // Not implemented

protected:
  // data from TreeK
public:
  KineTools();
  virtual ~KineTools(){}
 
  // data from TreeTR
  void SetDaughterPathMarks(TEveElement* cont, AliStack* stack, Bool_t recurse=kFALSE);
  void SetTrackReferences  (TEveElement* cont, TTree* treeTR=0, Bool_t recurse=kFALSE);
  void SortPathMarks       (TEveElement* cont, Bool_t recurse=kFALSE);

  ClassDef(KineTools, 1);
}; // endclass KineTools

} // end namespace Alieve

#endif
