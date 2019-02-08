///
/// AliFemtoEventReaderAlt.h
///

#include "AliFemtoEventReaderAODMultSelection.h"


class AliFemtoEventReaderAlt : public AliFemtoEventReaderAODMultSelection {
public:
  AliFemtoEventReaderAlt();
  ~AliFemtoEventReaderAlt();

protected:
  virtual AliFemtoTrack* CopyAODtoFemtoTrack(AliAODTrack *src);
};

