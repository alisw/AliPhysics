#ifndef ALIEMCALCORRECTIONCELLEMULATECROSSTALK_H
#define ALIEMCALCORRECTIONCELLEMULATECROSSTALK_H

#include <TRandom3.h>

#include "AliEmcalCorrectionComponent.h"

/**
 * @class AliEmcalCorrectionCellEmulateCrosstalk
 * @ingroup EMCALCOREFW
 * @brief Correction component to emulate cell-level crosstalk in the EMCal correction framework.
 *
 * Performs energy smearing of cells to mimic T-card crosstalk. The original cell information in the event **will be overwritten**.
 *
 * Based on code in AliAnalysisTaskEMCALClusterize.
 *
 * @author Gustavo Conesa Balbastre, LPSC-Grenoble, AliAnalysisTaskEMCALClusterize
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Feb 11 2018
 */

class AliEmcalCorrectionCellEmulateCrosstalk : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellEmulateCrosstalk();
  virtual ~AliEmcalCorrectionCellEmulateCrosstalk();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  
  /// Constant energy lost by max energy cell in one of T-Card cells, same for all SM
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossConstant(Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    for(Int_t ism = 0; ism < 22; ism++) {
      fTCardCorrInduceEner[0][ism] = ud; fTCardCorrInduceEner[1][ism] = udlr;
      fTCardCorrInduceEner[2][ism] = lr; fTCardCorrInduceEner[3][ism] = sec ; } }
  
  /// Fraction of energy lost by max energy cell in one of T-Card cells, same for all SM
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFraction(Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    for(Int_t ism = 0; ism < 22; ism++) {
      fTCardCorrInduceEnerFrac[0][ism] = ud; fTCardCorrInduceEnerFrac[1][ism] = udlr;
      fTCardCorrInduceEnerFrac[2][ism] = lr; fTCardCorrInduceEnerFrac[3][ism] = sec ; } }
  
  /// Slope parameter of fraction of energy lost by max energy cell in one of T-Card cells, same for all SM
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionP1(Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    for(Int_t ism = 0; ism < 22; ism++) {
      fTCardCorrInduceEnerFracP1[0][ism] = ud; fTCardCorrInduceEnerFracP1[1][ism] = udlr;
      fTCardCorrInduceEnerFracP1[2][ism] = lr; fTCardCorrInduceEnerFracP1[3][ism] = sec ; } }
  
  /// Constant energy lost by max energy cell in one of T-Card cells, per SM
  /// \param sm super module index
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossConstantPerSM(Int_t sm, Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    if ( sm < 22 && sm >= 0 ) {
      fTCardCorrInduceEner[0][sm] = ud; fTCardCorrInduceEner[1][sm] = udlr;
      fTCardCorrInduceEner[2][sm] = lr; fTCardCorrInduceEner[3][sm] = sec ; } }
  
  /// Fraction of energy lost by max energy cell in one of T-Card cells, per SM
  /// \param sm super module index
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionPerSM(Int_t sm, Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    if ( sm < 22 && sm >= 0 ) {
      fTCardCorrInduceEnerFrac[0][sm] = ud; fTCardCorrInduceEnerFrac[1][sm] = udlr;
      fTCardCorrInduceEnerFrac[2][sm] = lr; fTCardCorrInduceEnerFrac[3][sm] = sec ; } }
  
  /// Slope parameter of fraction of energy lost by max energy cell in one of T-Card cells, per SM
  /// \param sm super module index
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionP1PerSM(Int_t sm, Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    if ( sm < 22 && sm >= 0 ) {
      fTCardCorrInduceEnerFracP1[0][sm] = ud; fTCardCorrInduceEnerFracP1[1][sm] = udlr;
      fTCardCorrInduceEnerFracP1[2][sm] = lr; fTCardCorrInduceEnerFracP1[3][sm] = sec ; } }
  
  /// Fraction of energy lost by max energy cell in one of T-Card cells, width of random gaussian, same for all SM
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionWidth(Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    for(Int_t ism = 0; ism < 22; ism++) {
      fTCardCorrInduceEnerFracWidth[0][ism] = ud; fTCardCorrInduceEnerFracWidth[1][ism] = udlr;
      fTCardCorrInduceEnerFracWidth[2][ism] = lr; fTCardCorrInduceEnerFracWidth[3][ism] = sec ; } }
  
  /// Fraction of energy lost by max energy cell in one of T-Card cells, width of random gaussian, per SM
  /// \param sm super module index
  /// \param ud energy lost in upper/lower cell, same column
  /// \param udlr energy lost in upper/lower cell, left or right
  /// \param lr   energy lost in left or right cell, same row
  void           SetInducedEnergyLossFractionWidthPerSM(Int_t sm, Float_t ud, Float_t udlr, Float_t lr, Float_t sec) {
    if ( sm < 22 && sm >= 0 ) {
      fTCardCorrInduceEnerFracWidth[0][sm] = ud; fTCardCorrInduceEnerFracWidth[1][sm] = udlr;
      fTCardCorrInduceEnerFracWidth[2][sm] = lr; fTCardCorrInduceEnerFracWidth[3][sm] = sec ; } }
  
  /// Maximum induced energy fraction when linear dependency is set, per SM number
  /// \param max maximum fraction
  /// \param sm  super-module number
  void           SetInducedEnergyLossMaximumFractionPerSM(Float_t max, Int_t sm) {
    if ( sm < 22 && sm >= 0 ) fTCardCorrInduceEnerFracMax[sm] = max ; }
  
  /// Minimum induced energy fraction when linear dependency is set, per SM number
  /// \param min minimum fraction
  /// \param sm  super-module number
  void           SetInducedEnergyLossMinimumFractionPerSM(Float_t min, Int_t sm) {
    if ( sm < 22 && sm >= 0 ) fTCardCorrInduceEnerFracMin[sm] = min ; }
  
  /// Maximum induced energy fraction when linear dependency is set, same for all SM
  /// \param max maximum fraction
  void           SetInducedEnergyLossMaximumFraction(Float_t max) {
    for(Int_t ism = 0; ism < 22; ism++) fTCardCorrInduceEnerFracMax[ism] = max ; }
  
  /// Minimum induced energy fraction when linear dependency is set, same for all SM
  /// \param min minimum fraction
  void           SetInducedEnergyLossMinimumFraction(Float_t min) {
    for(Int_t ism = 0; ism < 22; ism++) fTCardCorrInduceEnerFracMin[ism] = min ; }
  
  /// fraction of times max cell energy correlates with cross cells, different for each super-module
  /// \param prob probability per event, from 0 to 1
  /// \param sm   probability assigned to this super-module number
  void           SetInducedEnergyLossProbabilityPerSM(Float_t prob, Int_t sm) {
    if ( sm < 22 && sm >= 0 ) fTCardCorrInduceEnerProb[sm] = prob ; }
  
  void           SwitchOnRandomizeTCardInducedEnergy()          { fRandomizeTCard = kTRUE   ; }
  void           SwitchOffRandomizeTCardInducedEnergy()         { fRandomizeTCard = kFALSE  ; }
  
  void           SetInducedTCardMinimumCellEnergy(Float_t mi)   { fTCardCorrMinAmp     = mi ; }
  void           SeInducedTCardMaximum(Float_t ma)              { fTCardCorrMaxInduced = ma ; }
  
  void           PrintTCardParam();
  
protected:
  
  // T-Card correlation emulation, do on MC
  void           MakeCellTCardCorrelation();
  void           AddInducedEnergiesToExistingCells();
  void           AddInducedEnergiesToNewCells();
  
  virtual void   ResetArrays();
  Bool_t         AcceptCell(Int_t absID);
  
  TH1F* fCellEnergyDistBefore;        //!<! cell energy distribution, before energy smearing
  TH1F* fCellEnergyDistAfter;         //!<! cell energy distribution, after energy smearing
  
  static const Int_t fgkNEMCalCells = 17664;       ///< Total number of cells in the calorimeter, 10*48*24 (EMCal) + 4*48*8 (EMCal/DCal 1/3) + 6*32*24 (DCal)
  
  // T-Card correlation emulation, do on MC
  Bool_t                fTCardCorrClusEnerConserv; ///< When making correlation, subtract from the reference cell the induced energy on the neighbour cells
  Float_t               fTCardCorrCellsEner[fgkNEMCalCells]; ///<  Array with induced cell energy in T-Card neighbour cells
  Bool_t                fTCardCorrCellsNew [fgkNEMCalCells]; ///<  Array with induced cell energy in T-Card neighbour cells, that before had no signal
  
  Float_t               fTCardCorrInduceEner         [4 ][22]; ///< Induced energy loss gauss constant on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param 0
  Float_t               fTCardCorrInduceEnerFrac     [4 ][22]; ///< Induced energy loss gauss fraction param0 on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param 0
  Float_t               fTCardCorrInduceEnerFracP1   [4 ][22]; ///< Induced energy loss gauss fraction param1 on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param1
  Float_t               fTCardCorrInduceEnerFracWidth[4 ][22]; ///< Induced energy loss gauss witdth on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells
  Float_t               fTCardCorrInduceEnerFracMax[22];   ///< In case fTCardCorrInduceEnerFracP1  is non null, restrict the maximum fraction of induced energy per SM
  Float_t               fTCardCorrInduceEnerFracMin[22];   ///< In case fTCardCorrInduceEnerFracP1  is non null, restrict the minimum fraction of induced energy per SM
  Float_t               fTCardCorrInduceEnerProb[22];      ///< Probability to induce energy loss per SM
  
  TRandom3              fRandom   ;                ///<  Random generator
  Bool_t                fRandomizeTCard ;          ///<  Use random induced energy
  
  Float_t               fTCardCorrMinAmp;          ///<  Minimum cell energy to induce signal on adjacent cells
  Float_t               fTCardCorrMaxInduced;      ///<  Maximum induced energy signal on adjacent cells
  
  Bool_t                fPrintOnce;                ///< Print once analysis parameters
  
private:
  
  AliEmcalCorrectionCellEmulateCrosstalk(const AliEmcalCorrectionCellEmulateCrosstalk &);               // Not implemented
  AliEmcalCorrectionCellEmulateCrosstalk &operator=(const AliEmcalCorrectionCellEmulateCrosstalk &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellEmulateCrosstalk> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellEmulateCrosstalk, 1); // EMCal cell crosstalk emulation correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCellEmulateCrosstalk_H */
