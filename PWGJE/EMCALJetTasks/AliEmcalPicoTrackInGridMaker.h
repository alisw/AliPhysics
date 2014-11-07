#ifndef ALIEMCALPICOTRACKINGRIDMAKER_H
#define ALIEMCALPICOTRACKINGRIDMAKER_H

class TClonesArray;
class AliVTrack;
class AliVParticle;
class TProfile;
class TH3F;
class AliEmcalJet;

#include "AliAnalysisTaskEmcalJet.h"

class AliEmcalPicoTrackInGridMaker : public AliAnalysisTaskEmcalJet {
 public:
  AliEmcalPicoTrackInGridMaker();
  AliEmcalPicoTrackInGridMaker(const char *name);
  virtual ~AliEmcalPicoTrackInGridMaker();

  void               SetTracksOutName(const char *name)                { fTracksOutName     = name; }

  void               SetCellSize(Double_t a)                           { fCellSize          = a;    }
  void               SetMinCellE(Double_t e)                           { fMinCellE          = e;    }
  void               SetL1Slide(Bool_t b)                              { fL1Slide           = b;    }

  void               SetEMCalAcceptance(Double_t phiMin, Double_t phiMax, Double_t etaMin, Double_t etaMax) { fPhiMin[0]=phiMin; fPhiMax[0]=phiMax; fEtaMin[0]=etaMin; fEtaMax[0]=etaMax; }
  void               SetDCalAcceptance(Double_t phiMin, Double_t phiMax, Double_t etaMin, Double_t etaMax)  { fPhiMin[1]=phiMin; fPhiMax[1]=phiMax; fEtaMin[1]=etaMin; fEtaMax[1]=etaMax; }

  void               SetExcludeLeadingPatch(Int_t i)                   { fExclLeadingPatch = i; }

  void               SetPatchTypeForSubtraction(Int_t i)               { fPatchSub = i; }
  void               SetPatchTypeForSubtraction(Int_t dim, Int_t lev)  { fPatchSub = GetPatchType(dim,lev); }
  void               SetMeanRho(Double_t r)                            { fRhoMean = r; }

 protected:
  void               UserCreateOutputObjects();
  Bool_t             Run();

  Bool_t             InitCells();
  Bool_t             CreateGridCells();
  Bool_t             CheckEdges();
  void               PrintAcceptance() const;

  Bool_t             InitMiniPatches();
  Bool_t             CreateGridMiniPatches();
  Bool_t             CreateGridPatches(const Int_t dim, const Int_t level);

  Bool_t             InitPatches(const Int_t dim, const Int_t level); //give dimension in cell units

  AliEmcalJet*       GetClosestJet(const Double_t eta, const Double_t phi, const Int_t icont = 0) const;

  Double_t           CalculateMedian(const Int_t patchType, const Int_t type, const Int_t areaType = 0);
  Double_t           CalculateSum(const Int_t patchType) const;

  //Getters
  Int_t              GetCellType(const Double_t eta, const Double_t phi) const;
  Int_t              GetCellType(const AliVParticle *vp) const {return GetCellType(vp->Eta(),vp->Phi());}
 
  Int_t              GetGridID(const AliVParticle *vp) const {return GetGridID(vp->Eta(),vp->Phi());}
  Int_t              GetGridID(const Double_t eta, const Double_t phi) const;
  Int_t              GetGridID(const Int_t row, const Int_t col, const Int_t type) const;
  void               GetEtaPhiFromGridID(const Int_t id, const Int_t type, Double_t &eta, Double_t &phi) const;
  Int_t              GetNCellsRow(const Int_t type) const;
  Int_t              GetNCellsCol(const Int_t type) const;

  Int_t              GetNRowMiniPatches(const Int_t type) const;
  Int_t              GetNColMiniPatches(const Int_t type) const;
  Int_t              GetMiniPatchID(const Int_t row, const Int_t col, const Int_t type) const;
  void               GetEtaPhiFromMiniPatchID(const Int_t id, const Int_t type, Double_t &eta, Double_t &phi) const;

  Int_t              GetPatchType(const Int_t dim, const Int_t level) const;
  Int_t              GetPatchDim(const Int_t ipatch) const;
  Int_t              GetSlidingStepSizeCells(const Int_t dim, const Int_t level = 1) const;
  Int_t              GetSlidingStepSizeMiniPatches(const Int_t dim, const Int_t level = 1) const;
  Int_t              GetTriggerPatchIdStepSizeNoOverlap(const Int_t dim, const Int_t level = 1) const;
  Int_t              GetNTriggerPatches(const Int_t type, const Int_t dim, const Int_t level) const;
  Int_t              GetNColTriggerPatches(const Int_t type, const Int_t dim, const Int_t patchType) const ;
  Int_t              GetNRowTriggerPatches(const Int_t type, const Int_t dim, const Int_t patchType) const;
  Int_t              GetTriggerPatchID(const Int_t row, const Int_t col, const Int_t type, const Int_t dim, const Int_t patchType) const;
  void               GetEtaPhiFromTriggerPatchID(const Int_t id, const Int_t type, const Int_t dim, const Int_t level, Double_t &eta, Double_t &phi) const;

  Double_t           GetPatchArea(const Int_t ipatch) const;
  Double_t           GetPatchAreaActive(const Int_t id, const Int_t type, const Int_t ipatch, const Int_t atype) const;

  TString            fTracksOutName;        // name of output track array
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
  TArrayD            fCellGrid[2];          // grid of cells in EMCal and DCal
  TArrayD            fMiniPatchGrid[2];     // grid of mini patches in EMCal and DCal
  TArrayI            fActiveAreaMP[2];      // active area for each mini patch
  TArrayD            fPatchGrid[2][5];      // grid of trigger patches: 4x4 L0, 4x4 L1, 8x8 L1, 16x16 L1, 32x32 L1
  TArrayI            fActiveAreaMPP[2][5];  // active area in mini patches for each trigger patch
  TArrayI            fActiveAreaCP[2][5];   // active area in cells for each trigger patch
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

  TH1F              *fPatchERaw[2][5];      //! uncorrected patch energy for EMCal and DCal
  TH2F              *fPatchECorrRho[2][5];  //! corrected patch energy vs rho opposite side
  TH3F              *fPatchECorrECorrRho[2][5]; //! Ecorr,det1 vs Ecorr,det2 vs rho,det2 opposite side for dijet in acceptance like events
  TH2F              *fh2PatchEtaPhiEmcal[5];    //! patch positions in EMCal
  TH2F              *fh2PatchEtaPhiDcal[5];     //! patch positions in DCal
  
  //jet histos
  TH2F              *fh2JetPtPatchECorr[2][5];  //! jet pt vs leading patch energy

  TH2F              *fMultVsRho;            //! track multiplicity vs rho from EMCal

  ClassDef(AliEmcalPicoTrackInGridMaker, 3); // Task to make PicoTracks in a grid corresponding to EMCAL/DCAL acceptance
};
#endif
