#ifndef ALIVERTEXERHYPERTRITON3BODY_H
#define ALIVERTEXERHYPERTRITON3BODY_H

class TClonesArray; /// This will be removed as soon as alisw/AliRoot#898 is merged and a new tag is available

#include <AliVertexerTracks.h>

class AliESDVertex;
class AliESDtrack;
class AliExternalTrackParam;

class AliVertexerHyperTriton3Body {
public:
  AliVertexerHyperTriton3Body();
  ~AliVertexerHyperTriton3Body();

  AliESDVertex *GetCurrentVertex() { return mCurrentVertex; }
  int GetGuessCompatibility() { return mCurrentGuessCompatibility; }

  bool FindDecayVertex(AliExternalTrackParam *deuteronTrack, AliExternalTrackParam *protonTrack,
                       AliExternalTrackParam *pionTrack, float b);
  static void Find2ProngClosestPoint(AliExternalTrackParam *track1, AliExternalTrackParam *track2, float b, float *pos);

  void SetMaxDinstanceInit(float maxD) { mMaxDistanceInitialGuesses = maxD; }
  void SetToleranceGuessCompatibility(int tol) { mToleranceGuessCompatibility = tol; }

  AliVertexerTracks mVertexerTracks;

private:
  AliESDVertex *mCurrentVertex;
  int mCurrentGuessCompatibility;

  float mMaxDistanceInitialGuesses;
  int mToleranceGuessCompatibility;
};

#endif