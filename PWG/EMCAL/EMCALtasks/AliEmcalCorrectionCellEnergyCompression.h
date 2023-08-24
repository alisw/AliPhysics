#ifndef ALIEMCALCORRECTIONCELLENERGYCOMPRESSION_H
#define ALIEMCALCORRECTIONCELLENERGYCOMPRESSION_H

#include "AliEmcalCorrectionComponent.h"

#include "TF1.h"

/**
 * @class AliEmcalCorrectionCellEnergyCompression
 * @ingroup EMCALCORRECTIONFW
 * @brief Cell energy Compression component in the EMCal correction framework.
 *
 * This component is to test different compression algorithms as planned to be used in O2.
 * The cell energy is multiplied with a certain Compression factor (typially 1/0.0153) and is then stored as int16_t (Short_t)
 * To retrieve the energy the int16_t (Short_t) is devided by the Compression factor
 * All
 * @author Joshua Koenig, joshua.konig@cern.ch
 * @date Apr 25 2022
 */

class AliEmcalCorrectionCellEnergyCompression : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellEnergyCompression();
  virtual ~AliEmcalCorrectionCellEnergyCompression();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  void ExecOnce();
  Bool_t Run();

protected:
  void SetEnergy(Double_t e);
  Double_t GetEnergy() const;
  void SetTime(Double_t t);
  Double_t GetTime() const;

  TH1F                  *fCellEnergyCompressionLoss;            //!<! difference between energy with no compression and with compression
  TH1F                  *fCellTimeCompressionLoss;              //!<! difference between cell time with no compression and with compression

  Bool_t                 fDoEnergyCompression;                  ///< switch to turn on energy Compression
  Short_t                fEnergy;                               ///< current energy
  Double_t               fTruncationMaxE;                       ///< Maximum cell energy
  Double_t               fEnergyCompression;                    ///< Energy Compression
  Double_t               fCompressionBias;                      ///< bias added to energy before converting it to int16_t (Short_t)
  Bool_t                 fDoTimeCompression;                    ///< switch to turn on time Compression
  Short_t                fTime;                                 ///< current time
  Double_t               fTimeCompression;                      ///< Time Compression
  Double_t               fTimeShift;                            ///< Time shift (~600ns for Run2 before time calib)

 private:
  AliEmcalCorrectionCellEnergyCompression(const AliEmcalCorrectionCellEnergyCompression &);               // Not implemented
  AliEmcalCorrectionCellEnergyCompression &operator=(const AliEmcalCorrectionCellEnergyCompression &);    // Not implemented

  // Allows the registration of the class so that it is available to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellEnergyCompression> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellEnergyCompression, 1); // EMCal cell energy variation component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCELLENERGYCOMPRESSION_H */
