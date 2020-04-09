/// \file Track.h
/// \brief O2 interface with AliRoot under the hood.

#ifndef ALICEO2_BASE_TRACK
#define ALICEO2_BASE_TRACK

#include "AliExternalTrackParam.h"
#include "MathUtils.h"
#include "Primitive2D.h"

namespace o2
{

namespace track
{
// aliases for track elements

class TrackParCov : public AliExternalTrackParam
{ // track+error parameterization
 public:

  TrackParCov() : AliExternalTrackParam() {}
  TrackParCov(const TrackParCov &t) : AliExternalTrackParam(static_cast<AliExternalTrackParam>(t)) {}
  virtual ~TrackParCov() {}

  const double* getParams() const { return this->GetParameter(); }
  float getParam(int i) const { return getParams()[i]; }
  float getX() const { return this->GetX(); }
  float getAlpha() const { return this->GetAlpha(); }
  float getY() const { return this->GetY(); }
  float getZ() const { return this->GetZ(); }
  float getSnp() const { return this->GetSnp(); }
  float getTgl() const { return this->GetTgl(); }
  float getQ2Pt() const { return this->GetSigned1Pt(); }


  // derived getters
  float getCurvature(float b) const { return this->GetC(b); }
  float getSign() const { return this->GetSign(); }
  float getPhi() const { return this->Phi(); }
  float getPhiPos() const { return this->PhiPos(); }

  float getP() const { return this->P(); }
  float getPt() const { return this->Pt(); }

  float getTheta() const { return this->Theta(); }
  float getEta() const { return -std::log(std::tan(0.5 * getTheta())); }
  // void getXYZGlo(std::array<float, 3>& xyz) const;
  // bool getPxPyPzGlo(std::array<float, 3>& pxyz) const;
  // bool getPosDirGlo(std::array<float, 9>& posdirp) const;

  // // methods for track params estimate at other point
  // bool getYZAt(float xk, float b, float& y, float& z) const;
  // float getZAt(float xk, float b) const;
  // float getYAt(float xk, float b) const;

  // parameters manipulation
  bool rotateParam(float alpha) { return this->RotateParamOnly(alpha); }
  bool propagateParamTo(float xk, float b) { return this->PropagateParamOnlyTo(xk, b); }


  const double* getCov() const { return this->GetCovariance(); }
  float getSigmaY2() const { return this->GetSigmaY2(); }
  float getSigmaZY() const { return this->GetSigmaZY(); }
  float getSigmaZ2() const { return this->GetSigmaZ2(); }
  float getSigmaSnpY() const { return this->GetSigmaSnpY(); }
  float getSigmaSnpZ() const { return this->GetSigmaSnpZ(); }
  float getSigmaSnp2() const { return this->GetSigmaSnp2(); }
  float getSigmaTglY() const { return this->GetSigmaTglY(); }
  float getSigmaTglZ() const { return this->GetSigmaTglZ(); }
  float getSigmaTglSnp() const { return this->GetSigmaTglSnp(); }
  float getSigmaTgl2() const { return this->GetSigmaTgl2(); }
  float getSigma1PtY() const { return this->GetSigma1PtY(); }
  float getSigma1PtZ() const { return this->GetSigma1PtZ(); }
  float getSigma1PtSnp() const { return this->GetSigma1PtSnp(); }
  float getSigma1PtTgl() const { return this->GetSigma1PtTgl(); }
  float getSigma1Pt2() const { return this->GetSigma1Pt2(); }

  // parameters + covmat manipulation
  bool rotate(float alpha) { return this->Rotate(alpha); }
  bool propagateTo(float xk, float b) { return this->PropagateTo(xk, b); }
  void invert() { this->Invert(); }

  void getCircleParamsLoc(float bz, o2::utils::CircleXY& c) const;
  void getCircleParams(float bz, o2::utils::CircleXY& c, float& sna, float& csa) const;

};

//_______________________________________________________
inline void TrackParCov::getCircleParams(float bz, o2::utils::CircleXY& c, float& sna, float& csa) const
{
  // get circle params in loc and lab frame
  getCircleParamsLoc(bz, c);
  o2::utils::sincosf(getAlpha(), sna, csa);
  o2::utils::rotateZ(c.xC, c.yC, c.xC, c.yC, sna, csa); // center in global frame
}

//_______________________________________________________
inline void TrackParCov::getCircleParamsLoc(float bz, o2::utils::CircleXY& c) const
{
  // get circle params in track local frame
  c.rC = 1.f / getCurvature(bz);
  float sn = getSnp(), cs = sqrtf((1. - sn) * (1. + sn));
  c.xC = getX() - sn * c.rC; // center in tracking
  c.yC = getY() + cs * c.rC; // frame. Note: r is signed!!!
  c.rC = fabs(c.rC);
}

} // namespace track
} // namespace o2

#endif