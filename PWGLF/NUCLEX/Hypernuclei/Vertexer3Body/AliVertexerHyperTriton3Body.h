#ifndef ALIVERTEXERHYPERTRITON3BODY_H
#define ALIVERTEXERHYPERTRITON3BODY_H

#include <AliVertexerTracks.h>

class AliESDVertex;
class AliESDtrack;
class AliExternalTrackParam;

class AliVertexerHyperTriton3Body
{
public:
  AliVertexerHyperTriton3Body();

  AliESDVertex* GetCurrentVertex() { return mCurrentVertex; }

  bool FindDecayVertex(AliESDtrack *track1, AliESDtrack *track2, AliESDtrack* track3, float b);
  static void Find2ProngClosestPoint(AliExternalTrackParam *track1, AliExternalTrackParam *track2, float b, float* pos);

  void SetMaxDinstanceInit(float maxD) { mMaxDistanceInitialGuesses = maxD; }
  void SetToleranceGuessCompatibility(int tol) { mToleranceGuessCompatibility = tol; }

  AliVertexerTracks mVertexerTracks;

private:
  AliESDVertex* mCurrentVertex;

  float mPosition[3];
  float mCovariance[6];

  float mMaxDistanceInitialGuesses;
  int mToleranceGuessCompatibility;
};

#endif