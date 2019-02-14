///
/// AliFemtoEventReaderAlt.h
///

#include "AliFemtoEventReaderAODMultSelection.h"

class TRandom3;

/// \class AliFemtoEventReaderAlt
/// \brief Alternative Event Reader
///
/// AliFemtoEventReaderAODMultSelection
///
class AliFemtoEventReaderAlt : public AliFemtoEventReaderAODMultSelection {
public:

  AliFemtoEventReaderAlt();
  ~AliFemtoEventReaderAlt();

  void SetEnhanceSmearing(double);

protected:
  virtual AliFemtoTrack* CopyAODtoFemtoTrack(AliAODTrack *src);

  TRandom3 *fRng;
  double fEnhanceSmearing;
};
