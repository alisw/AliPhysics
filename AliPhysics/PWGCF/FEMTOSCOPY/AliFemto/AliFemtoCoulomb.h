/// \class AliFemtoCoulomb
/// \brief Coulomb correction class
///
/// This is a Coulomb correction class which:
/// 1. Reads in the dat from a file
/// 2. Performs a linear interpolation in R and creates any array of
///    interpolations
/// 3. Interpolates in eta and returns the Coulomb correction to user

#ifndef ALIFEMTOCOULOMB_H
#define ALIFEMTOCOULOMB_H

class TH1D;
class TH3D;
class AliFemtoPair;

class AliFemtoCoulomb {
public:
  AliFemtoCoulomb();
  AliFemtoCoulomb(const char *readFile, const double& radius, const double& charge);
  AliFemtoCoulomb(const AliFemtoCoulomb& aCoul);
  virtual ~AliFemtoCoulomb();

  AliFemtoCoulomb& operator=(const AliFemtoCoulomb& aCoul);

  void SetRadius(const double& radius);
  double GetRadius() const;
  void SetFile(const char *readFile);
  void SetChargeProduct(const double& charge);

  // These have different names so eta/Qinv don't confuse the compiler
  double CoulombCorrect(const double& eta);
  double CoulombCorrect(const double& eta, const double& radius);
  double CoulombCorrect(const AliFemtoPair* pair);
  double CoulombCorrect(const AliFemtoPair* pair, const double& radius);
  double CoulombCorrect(const double& mass,
                        const double& charge,
                        const double& radius,
                        const double& qInv);

  TH1D* CorrectionHistogram(const double& mass1,
                            const double& mass2,
                            const int& nBins,
                            const double& low,
                            const double& high);

#ifdef __ROOT__
  TH1D* CorrectionHistogram(const TH1D*, const double);
  TH3D* CorrectionHistogram(const TH3D*, const double);
#endif
private:
  double Eta(const AliFemtoPair* pair);          ///< Calculates eta
  void CreateLookupTable(const double& radius);  ///< Creates look-up table
  const char* fFile;                             ///< File to interpolate corrections from
  double fRadius;                                ///< Radius from previous iteration
  double fZ1Z2;                                  ///< Charge product of particles
  double fEta[1000];                             ///< interpolated Coulomb correction table
  double fCoulomb[1000];                         ///< interpolated Coulomb correction table
  int fNLines;                                   ///< Number of Eta's in lookup-table

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCoulomb, 0);
  /// \endcond
#endif
};


#endif
