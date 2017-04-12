#ifndef AliAnalysisTaskEmcalAlig_H
#define AliAnalysisTaskEmcalAlig_H

#include "AliAnalysisTaskEmcal.h"

class AliEMCALRecoUtils;
class AliPIDResponse;
class TH1F;
class TH2F;
class AliEMCALGeometry;

class AliAnalysisTaskEmcalAlig : public AliAnalysisTaskEmcal {
public:
    
    AliAnalysisTaskEmcalAlig();
    AliAnalysisTaskEmcalAlig(const char *name);
    virtual ~AliAnalysisTaskEmcalAlig();
    
    void                        UserCreateOutputObjects();
    void                        Terminate(Option_t *option);
    
protected:
    void                        ExecOnce();
    Bool_t                      FillHistograms();
    Bool_t                      Run();
    
    void                        DoTrackLoop();
    
    
    AliEMCALRecoUtils *fEMCALRecoUtils;
    AliEMCALGeometry *fEMCALGeo;
    AliPIDResponse *fPIDResponse;
    
    //Negative Particle Histograms;
    TH1F *fNPartPt; //!
    TH1F *fNPartPhi; //!
    TH1F *fNPartEta; //!
    
    //Positive Particle Histograms;
    TH1F *fPPartPt; //!
    TH1F *fPPartPhi; //!
    TH1F *fPPartEta; //!
    
    //PID Plots
    TH2F* fTPCnSgima; //! TPC Nsigma as function of pT
    TH2F* fEOverP; //! E/p as function of pT
    
    //Matching Residual Plots
    TH2F        **fElectronPhiRes; //!
    TH2F        **fElectronEtaRes; //!
    TH2F        **fElectronXRes; //!
    TH2F        **fElectronYRes; //!
    TH2F        **fElectronPositronZRes; //!
    
    TH2F        **fPositronPhiRes; //!
    TH2F        **fPositronEtaRes; //!
    TH2F        **fPositronXRes; //!
    TH2F        **fPositronYRes; //!
    
    TH2F        **fAllMatchedTracksPhiRes; //!
    TH2F        **fAllMatchedTracksEtaRes; //!
    
    TH2F        **fElectronPhiResRU; //!
    TH2F        **fElectronEtaResRU; //!
    TH2F        **fElectronXResRU; //!
    TH2F        **fElectronYResRU; //!
    TH2F        **fElectronPositronZResRU; //!
    
    TH2F        **fPositronPhiResRU; //!
    TH2F        **fPositronEtaResRU; //!
    TH2F        **fPositronXResRU; //!
    TH2F        **fPositronYResRU; //!

    
private:
    AliAnalysisTaskEmcalAlig(const AliAnalysisTaskEmcalAlig&)           ; // not implemented
    AliAnalysisTaskEmcalAlig &operator=(const AliAnalysisTaskEmcalAlig&); // not implemented
    
    ClassDef(AliAnalysisTaskEmcalAlig, 1);
};
#endif
