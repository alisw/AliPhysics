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
#include "TObject.h"

class AliEMCALRecoUtils;
class AliPIDResponse;
class TH1F;
class TH2F;
class AliEMCALGeometry;
class TTree;

class ElectronForAlignment : public TObject {
public:
    
    //General Track properties
    Short_t charge;
    Float_t pt;
    Float_t pz;
    Float_t eta_track;
    Float_t phi_track;
    Float_t x_track;
    Float_t y_track;
    Float_t z_track;
    Float_t zvtx;

    
    //cluster properties
    Float_t energy;
    Float_t M20;
    Float_t M02;
    Float_t eta_cluster;
    Float_t phi_cluster;
    Float_t x_cluster;
    Float_t y_cluster;
    Float_t z_cluster;

    //mathing properties using electron mass
    UShort_t super_module_number;
    Int_t distance_bad_channel;
    Bool_t is_in_fid_region;

    //PID properties
    Float_t n_sigma_electron_TPC;
    
    ElectronForAlignment()
    {
        charge = -9;
        pt = -999;
        pz = -999;
        eta_track = -999;
        phi_track = -999;
        x_track = -999;
        y_track = -999;
        z_track = -999;
        zvtx = -999;

        energy = -999;
        M20 = -999;
        M02 = -999;
        eta_cluster = -999;
        phi_cluster = -999;
        x_cluster = -999;
        y_cluster = -999;
        z_cluster = -999;
        

        super_module_number = 99;
        distance_bad_channel = -99;
        is_in_fid_region = kFALSE;

        n_sigma_electron_TPC = -999;
    }
    
    void Reset()
    {
               charge = -9;
        pt = -999;
        pz = -999;
        eta_track = -999;
        phi_track = -999;
        x_track = -999;
        y_track = -999;
        z_track = -999;
        zvtx = -999;


        energy = -999;
        M20 = -999;
        M02 = -999;
        eta_cluster = -999;
        phi_cluster = -999;
        x_cluster = -999;
        y_cluster = -999;
        z_cluster = -999;
        

        super_module_number = 99;
        distance_bad_channel = -99;
        is_in_fid_region = kFALSE;
        n_sigma_electron_TPC = -999;
    }
    
    //Default Initializer
    
    ClassDef(ElectronForAlignment, 3);
    
};


class AliAnalysisTaskEMCALAlig : public AliAnalysisTaskEmcal {
public:
    
    AliAnalysisTaskEMCALAlig();
    AliAnalysisTaskEMCALAlig(const char *name);
    virtual ~AliAnalysisTaskEMCALAlig();
    
    void UserCreateOutputObjects();
    void Terminate(Option_t *option);
    void SetSuffix(TString suff) { fTreeSuffix = suff;};
    
protected:
    void                        ExecOnce();
    Bool_t                      FillHistograms();
    Bool_t                      Run();
    
    void                        DoTrackLoop();
    
    
    AliEMCALRecoUtils *fEMCALRecoUtils; //! EMCAL Reco utils used to recalculate the matching
    AliEMCALGeometry *fEMCALGeo; //! EMCAL geometry class
    AliPIDResponse *fPIDResponse; //! PID response task used to perform electron identification
    
    ElectronForAlignment fElectronInformation; //! Object to hold the electron information
    TTree* fElectronTree; //! Electron tree output
    TString fTreeSuffix; // Suffix for tree name

    
private:
    AliAnalysisTaskEMCALAlig(const AliAnalysisTaskEMCALAlig&)           ; // not implemented
    AliAnalysisTaskEMCALAlig &operator=(const AliAnalysisTaskEMCALAlig&); // not implemented
    
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALAlig, 3);
  /// \endcond

};
#endif
