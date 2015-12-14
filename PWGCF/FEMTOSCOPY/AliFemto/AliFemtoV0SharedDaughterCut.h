#ifndef AliFemtoV0SharedDaughterCut_H
#define AliFemtoV0SharedDaughterCut_H

#include "AliFemtoTrackCut.h"
#include "AliFemtoV0.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoV0Collection.h"

class AliFemtoV0SharedDaughterCut {
public:

  AliFemtoV0SharedDaughterCut(); ///< default constructor, does nothing
  virtual AliFemtoV0Collection AliFemtoV0SharedDaughterCutCollection(AliFemtoV0Collection*, AliFemtoV0Cut*);
  virtual ~AliFemtoV0SharedDaughterCut();

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0SharedDaughterCut, 0);
  /// \endcond
#endif
};

#endif
