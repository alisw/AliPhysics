#ifndef ITSTRACK_H
#define ITSTRACK_H


////////////////////////////////////////////////////////////////////////
//           Track class for set: ITS                                 //
////////////////////////////////////////////////////////////////////////

#include <TObject.h>



//_____________________________________________________________________________
class AliITStrack : public TObject {

public:

  AliITStrack() {}
  AliITStrack(const AliITStrack& t) {}
  virtual ~AliITStrack() {}


  Bool_t IsSortable() const {return kTRUE;}

  ClassDef(AliITStrack,1)  // ITS reconstructed tracks
};

#endif
