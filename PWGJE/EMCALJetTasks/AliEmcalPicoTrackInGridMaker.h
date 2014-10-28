#ifndef ALIEMCALPICOTRACKINGRIDMAKER_H
#define ALIEMCALPICOTRACKINGRIDMAKER_H

class TClonesArray;
class AliVTrack;
class AliVParticle;
class TProfile;
class TH3F;

#include "AliAnalysisTaskSE.h"

class AliEmcalPicoTrackInGridMaker : public AliAnalysisTaskSE {
 public:
  AliEmcalPicoTrackInGridMaker();
  AliEmcalPicoTrackInGridMaker(const char *name);
  virtual ~AliEmcalPicoTrackInGridMaker();

  void               SetTracksInName(const char *name)                 { fTracksInName      = name; }
  void               SetTracksOutName(const char *name)                { fTracksOutName     = name; }

  void               SetCellSize(Double_t a)                           { fCellSize          = a;    }
  void               SetMinCellE(Double_t e)                           { fMinCellE          = e;    }

  void               SetEMCalAcceptance(Double_t phiMin, Double_t phiMax, Double_t etaMin, Double_t etaMax) { fPhiMin[0]=phiMin; fPhiMax[0]=phiMax; fEtaMin[0]=etaMin; fEtaMax[0]=etaMax; }
  void               SetDCalAcceptance(Double_t phiMin, Double_t phiMax, Double_t etaMin, Double_t etaMax) { fPhiMin[1]=phiMin; fPhiMax[1]=phiMax; fEtaMin[1]=etaMin; fEtaMax[1]=etaMax; }

  void               SetExcludeLeadingPatch(Int_t i)                   { fExclLeadingPatch = i; }

  void               SetPatchTypeForSubtraction(Int_t i)               { fPatchSub = i; }
  void               SetPatchTypeForSubtraction(Int_t dim, Int_t lev)  { fPatchSub = GetPatchType(dim,lev); }
  void               SetMeanRho(Double_t r)                            { fRhoMean = r; }

 protected:
  void               UserCreateOutputObjects();
  void               UserExec(Option_t *option);

  Bool_t             InitCells();
  Bool_t             CreateGridCells();
  Bool_t             CheckEdges();
  void               PrintAcceptance();

  Bool_t             InitMiniPatches();
  Bool_t             CreateGridMiniPatches();
  Bool_t             CreateGridPatches(Int_t dim, Int_t level);

  Bool_t             InitPatches(Int_t dim, Int_t level); //give dimension in cell units

  Double_t           CalculateMedian(Int_t patchType, Int_t type, Int_t areaType = 0);
  Double_t           CalculateSum(Int_t patchType);

  //Getters
  Int_t              GetGridID(AliVParticle *vp) {return GetGridID(vp->Eta(),vp->Phi());}
  Int_t              GetGridID(Double_t eta, Double_t phi);
  Int_t              GetGridID(Int_t row, Int_t col, Int_t type);
  Int_t              GetCellType(Double_t eta, Double_t phi);
  Int_t              GetNCellsRow(Int_t type);
  Int_t              GetNCellsCol(Int_t type);

  Int_t              GetNRowMiniPatches(Int_t type);
  Int_t              GetNColMiniPatches(Int_t type);
  Int_t              GetMiniPatchID(Int_t row, Int_t col, Int_t type);

  Int_t              GetPatchType(Int_t dim, Int_t level);
  Int_t              GetSlidingStepSizeCells(Int_t dim);
  Int_t              GetSlidingStepSizeMiniPatches(Int_t dim);

  Double_t           GetPatchArea(Int_t ipatch);
  Double_t           GetPatchAreaActive(Int_t id, Int_t ipatch, Int_t type);

  TString            fTracksOutName;        // name of output track array
  TString            fTracksInName;         // name of input jet array
  TClonesArray      *fTracksIn;             //!jet array in
  TClonesArray      *fTracksOut;            //!track array out
  Bool_t             fL1Slide;              // sliding window on

  Double_t           fPhiMin[2];            // min phi of EMCal an DCal
  Double_t           fPhiMax[2];            // max phi of EMCal an DCa
  Double_t           fEtaMin[2];            // min eta of EMCal an DCal
  Double_t           fEtaMax[2];            // max eta of EMCal an DCal
  Double_t           fCellSize;             // size of cell (equal in eta and phi)
  Double_t           fMinCellE;             // minimum cell energy
  Int_t              fExclLeadingPatch;     // exclude leading patch from median calculation

  Int_t              fPatchSub;             // patch type to use for subtraction
  Double_t           fRhoMean;              // mean rho

  Int_t              fNCells;               // total number of cells
  Int_t              fNCellsEMCal;          // total number of EMCal cells
  Int_t              fNCellsDCal;           // total number of DCal cells
  TArrayD            fCellGrid;             // grid of cells in EMCal and DCal
  TArrayD            fMiniPatchGrid;        // grid of mini patches in EMCal and DCal
  TArrayI            fActiveAreaMP;         // active area for each mini patch
  TArrayD            fPatchGrid[5];         // grid of trigger patches: 4x4 L0, 4x4 L1, 8x8 L1, 16x16 L1, 32x32 L1
  TArrayI            fActiveAreaMPP[5];     // active area in mini patches for each trigger patch
  TArrayI            fActiveAreaCP[5];      // active area in cells for each trigger patch
  Int_t              fNPatchesEMCal[5];     // number of patches in EMCal

 private:
  AliEmcalPicoTrackInGridMaker(const AliEmcalPicoTrackInGridMaker&);            // not implemented
  AliEmcalPicoTrackInGridMaker &operator=(const AliEmcalPicoTrackInGridMaker&); // not implemented

  TH2F              *fh2MedianTypeEmcal[3]; //! median vs patch type for 3 types of area calculation
  TH2F              *fh2MedianTypeDcal[3];  //! median vs patch type for 3 types of area calculation
  TProfile          *fpMedianTypeEmcal[3];  //! median vs patch type for 3 types of area calculation
  TProfile          *fpMedianTypeDcal[3];   //! median vs patch type for 3 types of area calculation
  TH1F              *fh1RhoEmcal[5];        //! rho distributions for passive area
  TH1F              *fh1RhoDcal[5];         //! rho distributions for passive area
  TH2F              *fPatchEnVsActivityEmcal[5]; //! patch energy vs active cells
  TH2F              *fPatchEnVsActivityDcal[5];  //! patch energy vs active cells

  TH1F              *fPatchECorr[2][5];     //! corrected patch energy for EMCal and DCal
  TH1F              *fPatchECorrPar[2][5];  //! corrected patch energy with inclusive mean rho for EMCal and DCal
  TH2F              *fPatchECorrRho[2][5];  //! corrected patch energy vs rho opposite side
  TH2F              *fPatchECorrRhoDijet[2][5]; //! corrected patch energy vs rho opposite side
  TH3F              *fPatchECorrECorrRho[2][5]; //! Ecorr,det1 vs Ecorr,det2 vs rho,det2 opposite side for dijet in acceptance like events

  TH2F              *fMultVsRho;            //! track multiplicity vs rho from EMCal

  TList             *fHistList;             //! List of Histograms

  ClassDef(AliEmcalPicoTrackInGridMaker, 1); // Task to make PicoTracks in a grid corresponding to EMCAL/DCAL acceptance
};
#endif
