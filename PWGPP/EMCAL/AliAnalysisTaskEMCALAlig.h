#ifndef AliAnalysisTaskEMCALAlig_H
#define AliAnalysisTaskEMCALAlig_H

//_________________________________________________________________________
/// \class AliAnalysisTaskEMCALAlig
/// \ingroup EMCALPerformance 
/// \brief Alignment checks
///
// Matching residual for electrons on EMCal.
// Both default residual (until the surface) and a improved one, using
// electron mass hypotehsis and propagation until the cluster, are available
///
/// \author Henrique Zanoli <Henrique.Zanoli@cern.ch>, University of Sao Paulo and Utrecht University
//_________________________________________________________________________


#include "AliAnalysisTaskEmcal.h"

class AliEMCALRecoUtils;
class AliPIDResponse;
class TH1F;
class TH2F;
class AliEMCALGeometry;

class AliAnalysisTaskEMCALAlig : public AliAnalysisTaskEmcal {
public:
    
    AliAnalysisTaskEMCALAlig();
    AliAnalysisTaskEMCALAlig(const char *name);
    virtual ~AliAnalysisTaskEMCALAlig();
    
    void                        UserCreateOutputObjects();
    void                        Terminate(Option_t *option);
    
protected:
    void                        ExecOnce();
    Bool_t                      FillHistograms();
    Bool_t                      Run();
    
    void                        DoTrackLoop();
    
    
    AliEMCALRecoUtils *fEMCALRecoUtils; //! EMCAL Reco utils used to recalculate the matching
    AliEMCALGeometry *fEMCALGeo; //! EMCAL geometry class
    AliPIDResponse *fPIDResponse; //! PID response task used to perform electron identification
    
    //Negative Particle Histograms;
    TH1F *fNPartPt; //! Electron transverse momemtum distribution
    TH1F *fNPartPhi; //! Electron azimuthal angle distribution
    TH1F *fNPartEta; //! Electron pseudorapidity distribution
    
    //Positive Particle Histograms;
    TH1F *fPPartPt; //! Positron transverse momemtum distribution
    TH1F *fPPartPhi; //! Positron azimuthal angle distribution
    TH1F *fPPartEta; //! Positron pseudorapidity distribution
    
    //PID Plots
    TH2F* fTPCnSgima; //! TPC Nsigma as function of pT
    TH2F* fEOverP; //! E/p as function of pT
    
    //Matching Residual Plots
    TH2F        **fElectronPhiRes; //! Electron matching residual in phi as function of EMCAL SuperModule using propagation to EMCAL surface
    TH2F        **fElectronEtaRes; //! Electron matching residual in eta as function of EMCAL SuperModule using propagation to EMCAL surface
    TH2F        **fElectronXRes; //! Electron matching residual in x as function of EMCAL SuperModule using propagation to EMCAL surface
    TH2F        **fElectronYRes; //! Electron matching residual in y as function of EMCAL SuperModule using propagation to EMCAL surface
    TH2F        **fElectronPositronZRes; //! Electron and positron matching residual in z as function of EMCAL SuperModule using propagation to EMCAL surface
    
    TH2F        **fPositronPhiRes; //! Positron matching residual in phi as function of EMCAL SuperModule using propagation to EMCAL surface
    TH2F        **fPositronEtaRes; //! Positron matching residual in eta as function of EMCAL SuperModule using propagation to EMCAL surface
    TH2F        **fPositronXRes; //! Positron matching residual in x as function of EMCAL SuperModule using propagation to EMCAL surface
    TH2F        **fPositronYRes; //! Positron  matching residual in y as function of EMCAL SuperModule using propagation to EMCAL surface
    
    TH2F        **fAllMatchedTracksPhiRes; //! All tracks mathing residual in Phi
    TH2F        **fAllMatchedTracksEtaRes; //! All tracks mathing residual in Eta
    
    TH2F        **fElectronPhiResRU; //! Electron matching residual in phi as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)
    TH2F        **fElectronEtaResRU; //! Electron matching residual in eta as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)
    TH2F        **fElectronXResRU; //! Electron matching residual in x as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)
    TH2F        **fElectronYResRU; //! Electron matching residual in y as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)
    TH2F        **fElectronPositronZResRU; //! Electron and positron matching residual in z as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)
    
    TH2F        **fPositronPhiResRU; //! Positron matching residual in phi as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)
    TH2F        **fPositronEtaResRU; //! Positron matching residual in eta as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)
    TH2F        **fPositronXResRU; //! Positron matching residual in x as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)
    TH2F        **fPositronYResRU; //! Positron  matching residual in y as function of EMCAL SuperModule using propagation to cluster (done by EMCALRecoUtils)

    
private:
    AliAnalysisTaskEMCALAlig(const AliAnalysisTaskEMCALAlig&)           ; // not implemented
    AliAnalysisTaskEMCALAlig &operator=(const AliAnalysisTaskEMCALAlig&); // not implemented
    
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALAlig, 1);
  /// \endcond

};
#endif
