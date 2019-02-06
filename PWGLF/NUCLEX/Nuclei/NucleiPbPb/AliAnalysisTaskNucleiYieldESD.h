#ifndef __AliAnalysisTaskNucleiYieldESD__
#define __AliAnalysisTaskNucleiYieldESD__

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>
#include <TString.h>
#include <TArrayF.h>
#include "AliPIDResponse.h"
#include "AliPID.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"

class AliPIDResponse;
class AliVTrack;

class TH1I;
class TH2F;
class TH3F;
class TList;
class AliESDtrackCuts;
class AliEventCuts;
class TDirectory;

class AliAnalysisTaskNucleiYieldESD : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskNucleiYieldESD(TString taskname = "NucleiYieldESDTask");
  virtual ~AliAnalysisTaskNucleiYieldESD();

  void SetParticleType(AliPID::EParticleType part);
  void SetIsMC (bool isMc) { fIsMC = isMc; }

  void SetCentBins (int nbins, float *bins);
  void SetDCABins (int nbins, float min, float max);
  void SetDCABins (int nbins, float* bins);
  void SetPtBins (int nbins, float *bins);
  void SetTOFBins (int nbins, float min, float max);
  void SetDCAzBins (int nbins, float limit);
  void SetYcut (float y) { fYregion = y; }
  void SetTPCpidNsigmaCut (float nSig) { fTPCnSigmaCut = nSig; }

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  AliESDtrackCuts    fTrackCuts;
  AliEventCuts fEventCuts;
private:
  AliAnalysisTaskNucleiYieldESD (const AliAnalysisTaskNucleiYieldESD &source);
  AliAnalysisTaskNucleiYieldESD &operator=(const AliAnalysisTaskNucleiYieldESD &source);

  float HasTOF(AliVTrack *t);

  bool Flatten (float cent);
  void PtCorrection (float &pt, bool pos);

  TF1                *fTOFfunction;           //!<! TOF signal function
  float               fTOFtail;               ///< tail of the TOF function

  TList              *fList;                  ///<  Output list
  float               fYregion;               ///<  Symmetric rapidity cut (|y| < fYregion)
  float               fTPCnSigmaCut;          ///<  TPC pid cut
  AliPID::EParticleType fParticle;            ///<  Particle specie
  int                 fPDG;                   ///<  PDG code of the particle of interest
  float               fPDGMass;               ///<  PDG mass
  float               fPDGMassOverZ;          ///<  PDG mass over z
  bool                fIsMC;                  ///<  Switch between MC and data

  AliPIDResponse     *fPID;                   //!<! PID response class
  float               fMagField;              ///<  Magnetic field value for the current event

  float               fDCAzLimit;             ///<  Limits of the \f$DCA_{z}\f$ histograms
  int                 fDCAzNbins;             ///<  Number of bins used for \f$DCA_{z}\f$ distributions

  TArrayF             fPtCorrectionA;         ///<  Array containing the parametrisation of the \f$p_{T}\$ correction for anti-matter
  TArrayF             fPtCorrectionM;         ///<  Array containing the parametrisation of the \f$p_{T}\$ correction for matter

  float               fTOFlowBoundary;        ///<  Lower limit for the TOF mass spectra histograms
  float               fTOFhighBoundary;       ///<  Upper limit for the TOF mass spectra histograms
  int                 fTOFnBins;              ///<  Number of bins used for the TOF mass spectra

  TArrayF             fCentBins;              ///<  Centrality bins
  TArrayF             fDCABins;               ///<  DCA bins
  TArrayF             fPtBins;                ///<  Transverse momentum bins
  TArrayF             fFlatteningProbs;       ///<  Flattening probabilities

  // Event related histograms
  TH1F               *fCentralityClasses;     //!<! Events statistics per centrality classes

  // MC only histograms
  TH1F               *fProduction[2];         //!<! *(MC only)* Total number of produced particles
  TH2F               *fITS_TPC[2];            //!<! *(MC only)* Tracks reconstructed in ITS-TPC acceptance
  TH2F               *fITS_TPC_TOF[2];        //!<! *(MC only)* Tracks reconstructed in ITS-TPC-TOF acceptance
  TH2F               *fTotal[2];              //!<! *(MC only)* Particles in acceptance (anti-matter)
  TH2F               *fPtCorrection[2];       //!<! *(MC only)* \f$p_{T}^{rec}-p_{T}^{MC}\f$ as a function of \f$p_{T}^{rec}\f$
  TH3F               *fDCAPrimaryTPC[2];      //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of primaries for ITS+TPC tracks
  TH3F               *fDCASecondaryTPC[2];    //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of secondaries for ITS+TPC tracks
  TH3F               *fDCAPrimaryTOF[2];      //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of primaries for ITS+TPC+TOF tracks
  TH3F               *fDCASecondaryTOF[2];    //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of secondaries for ITS+TPC+TOF tracks

  // Data histograms
  TH3F               *fTOFsignal[2];          //!<!
  TH3F               *fTPCcounts[2];          //!<!
  TH2F               *fTPCdEdx[2];            //!<!
  TH2F               *fTPCdEdxTpcCut[2];      //!<!
  TH2F               *fTPCdEdxTofCut[2];      //!<!
  TH3F               *fDCAxyTPC[2];           //!<! *(Data only)* \f$DCA_{xy}\f$ distribution for ITS+TPC tracks
  TH3F               *fDCAzTPC[2];            //!<! *(Data only)* \f$DCA_{z}\f$ distribution for ITS+TPC tracks
  TH3F               *fDCAxyTOF[2];           //!<! *(Data only)* \f$DCA_{xy}\f$ distribution for ITS+TPC+TOF tracks
  TH3F               *fDCAzTOF[2];            //!<! *(Data only)* \f$DCA_{z}\f$ distribution for ITS+TPC+TOF tracks
  TH3F               *fTOFtemplates[5];       //!<! *(Data only)* TOF signal templates for pi/k/p/d/t

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskNucleiYieldESD, 2);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskNucleiYieldESD__) */
