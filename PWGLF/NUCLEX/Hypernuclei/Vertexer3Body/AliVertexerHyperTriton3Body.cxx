#include "AliVertexerHyperTriton3Body.h"

#include <cmath>

#include <AliESDVertex.h>
#include <AliESDtrack.h>
#include <AliExternalTrackParam.h>
#include <TObjArray.h>

namespace
{
float Point2PointDistance(float *p0, float *p1)
{
  float d2 = 0.0;
  for (int iDim = 0; iDim < 3; ++iDim)
  {
    d2 += (p0[iDim] - p1[iDim]) * (p0[iDim] - p1[iDim]);
  }
  return std::sqrt(d2);
}
} // namespace

AliVertexerHyperTriton3Body::AliVertexerHyperTriton3Body() : mVertexerTracks{},
                                                             mCurrentVertex{nullptr},
                                                             mPosition{0.f, 0.f, 0.f},
                                                             mCovariance{999.f, 999.f, 999.f, 999.f, 999.f, 999.f},
                                                             mMaxDistanceInitialGuesses{4.f},
                                                             mToleranceGuessCompatibility{1}
{
}

void AliVertexerHyperTriton3Body::Find2ProngClosestPoint(AliExternalTrackParam *track1, AliExternalTrackParam *track2, float b, float *pos)
{
  /// Extracted from the constructor of AliESDv0

  track1->PropagateToDCA(track2, b);
  track2->PropagateToDCA(track1, b);

  constexpr double ss = 0.0005 * 0.0005; //a kind of a residual misalignment precision
  double tmp[3]{0.0, 0.0, 0.0};

  //Trivial estimation of the vertex parameters
  double alpha = track1->GetAlpha(), cs = std::cos(alpha), sn = std::sin(alpha);
  track1->GetXYZ(tmp);
  double x1 = tmp[0], y1 = tmp[1], z1 = tmp[2];
  double sx1 = sn * sn * track1->GetSigmaY2() + ss, sy1 = cs * cs * track1->GetSigmaY2() + ss;

  alpha = track2->GetAlpha();
  cs = std::cos(alpha);
  sn = std::sin(alpha);
  track2->GetXYZ(tmp);
  double x2 = tmp[0], y2 = tmp[1], z2 = tmp[2];
  double sx2 = sn * sn * track2->GetSigmaY2() + ss, sy2 = cs * cs * track2->GetSigmaY2() + ss;

  double sz1 = track1->GetSigmaZ2(), sz2 = track2->GetSigmaZ2();
  double wx1 = sx2 / (sx1 + sx2), wx2 = 1. - wx1;
  double wy1 = sy2 / (sy1 + sy2), wy2 = 1. - wy1;
  double wz1 = sz2 / (sz1 + sz2), wz2 = 1. - wz1;

  pos[0] = wx1 * x1 + wx2 * x2;
  pos[1] = wy1 * y1 + wy2 * y2;
  pos[2] = wz1 * z1 + wz2 * z2;
}

bool AliVertexerHyperTriton3Body::FindDecayVertex(AliESDtrack *track1, AliESDtrack *track2, AliESDtrack *track3, float b)
{
  float initialGuesses[3][3];
  AliExternalTrackParam *tracks[3]{
      static_cast<AliExternalTrackParam *>(track1),
      static_cast<AliExternalTrackParam *>(track2),
      static_cast<AliExternalTrackParam *>(track3)};

  TObjArray trArray(3);
  trArray.Add(track1);
  trArray.Add(track2);
  trArray.Add(track3);

  for (int iTrack = 0; iTrack < 3; ++iTrack)
  {
    Find2ProngClosestPoint(tracks[iTrack], tracks[(iTrack + 1) % 3], b, initialGuesses[iTrack]);
  }

  int agreement = 0;
  for (int iPoint = 0; iPoint < 3; ++iPoint)
  {
    float distance = Point2PointDistance(initialGuesses[iPoint], initialGuesses[(iPoint + 1) % 3]);
    if (distance < mMaxDistanceInitialGuesses)
    {
      agreement++;
    }
  }

  if (agreement < mToleranceGuessCompatibility)
  {
    return false;
  }

  for (int iDim = 0; iDim < 3; ++iDim)
  {
    for (int iPoint = 0; iPoint < 3; ++iPoint)
    {
      mPosition[iDim] += initialGuesses[iPoint][iDim];
    }
    mPosition[iDim] /= 3.f;
  }

  mVertexerTracks.SetFieldkG(b);
  mVertexerTracks.SetVtxStart(mPosition[0], mPosition[1], mPosition[2]);
  mCurrentVertex = mVertexerTracks.VertexForSelectedESDTracks(&trArray);
  return mCurrentVertex ? true : false;
}