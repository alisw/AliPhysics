// $Header$

// Tools for import of kinematics. 
// Preliminary/minimal solution.

#ifndef ALIEVE_KineTools_H
#define ALIEVE_KineTools_H

#include <Reve/Reve.h>
#include <TObject.h>

class TTree;
class AliStack;

namespace Reve {
  class TrackList;
}

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
  void SetDaughterPathMarks(Reve::TrackList* cont,  AliStack* stack);
  void SetTrackReferences(Reve::TrackList* cont, TTree* treeTR = 0);

  ClassDef(KineTools, 1);
}; // endclass KineTools

} // end namespace Alieve

#endif
