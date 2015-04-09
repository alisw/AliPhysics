#ifndef AliAnalysisTaskHFEDiJets_cxx
#define AliAnalysisTaskHFEDiJets_cxx

//QA task for EMCAL electron analysis

class TH1F;
class THnSparse;
class AliESDEvent;
class AliAODEvent;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class TLorentzVector;
class TGraphErrors;
class AliMagF;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"

#include <fastjet/config.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
# include <fastjet/ClusterSequenceArea.hh>
# include <fastjet/AreaDefinition.hh>

class AliAnalysisTaskHFEDiJets : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskHFEDiJets();
    AliAnalysisTaskHFEDiJets(const char *name);
    virtual ~AliAnalysisTaskHFEDiJets();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    virtual void  Clear(const Option_t* /*opt*/ = "");
    
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    
    void SetJetRcut  	(Double_t JetRcut)		{ fJetRcut = JetRcut; };

    void EMCalClusterInfo();
    void InputTracksForJetReco(Double_t &fMaxpT);
    void DoJetReconstruction();
    void DoBkgJetReconstruction();
    void PlotRecoBkgJetProperties(TH1F* HisBkgArea, TH1F* HisBkgPhi, TH1F* HisBkgEta, TH1F* HisBkgPt, TH1F* HisBkgP);
    void PlotRecoJetProperties(TH1F *HisSigArea,TH1F *HisSigPhi, TH1F *HisSigEta, TH1F *HisSigPt_BfrSub, TH1F *HisSigPt_Sub, TH1F *HisSigP, TH1F *HisSigR, TH2F *HisSigEtaPhi, TH1F *HisSigEtainEMC);
    void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);
    Bool_t TagElectronJet(vector<fastjet::PseudoJet> constJet, AliVTrack *Eletrack, Double_t &maxpT, Double_t &maxP, Double_t &maxPx, Double_t &maxPy, Double_t &maxPz);
    void PlotJetswithElec(AliVTrack *Eletrack,TH1F *HisNthJetwEle,TH1F *HisPtJetwElec,TH1F *HisPtElecinJet,TH1F *HisPtElecLeadInJet,TH2F *HisEtaPhiJetwElec,TH1F *HisEtaJetwElecinEMC);
    void ElecJetDeltaPhi(unsigned iTrigJet, AliVTrack *Eletrack, THnSparse *fSprDphiDeta);
    
private:
    enum{
        kAODanalysis = BIT(20),
    };
    
    AliVEvent   *fVevent;  //!event object
    AliESDEvent *fESD;    //!ESD object
    AliAODEvent *fAOD;    //!AOD object
    AliPIDResponse *fpidResponse; //!pid response
    
    Bool_t      fFlagSparse;// switch to THnspare
    
    Double_t	fJetRcut;// Jet R cut
    vector<fastjet::PseudoJet> fInputParticlesJet;//!vector of input tracks for jet reco
    Double_t    fBkgMedian;//!
    Double_t    fBkgSigma;//!
    Double_t    fBkgMeanArea;//!
    Double_t    fInvmassCut;//!
    Double_t    fMinJetPt;//!
    
    fastjet::JetDefinition *fJetDef;//!
    fastjet::GhostedAreaSpec *fGhostArea;//!
    fastjet::AreaDefinition *fAreaDef;//!
    fastjet::ClusterSequenceArea *fCSA;//!
    vector<fastjet::PseudoJet> fJets;//!
    
    fastjet::JetDefinition *fBkgJetDef;//!
    fastjet::GhostedAreaSpec *fBkgGhostArea;//!
    fastjet::AreaDefinition *fBkgAreaDef;//!
    fastjet::ClusterSequenceArea *fBkgCSA;//!
    fastjet::RangeDefinition *fRange; //!
    vector<fastjet::PseudoJet> fBkgJets;//!
    
    //---General histos
    TList       *fOutputList; //!Output list
    TH1F        *fNevents;//!no of events
    TH1F        *fHistClustE;//!cluster energy
    TH2F        *fEMCClsEtaPhi;//!EMC cluster eta and phi
    TH1F        *fTrkPt;//!track pt
    TH1F        *fTrketa;//!track eta
    TH1F        *fTrkphi;//!track phi
    TH2F        *fdEdx;//!dedx vs pt
    TH2F        *fTPCnsig;//!TPC Nsigma
    TH1F        *fHistPtMatch;//!Pt tracks matched to EMCAL
    TH2F        *fEMCTrkMatch;//!Distance of EMC cluster to closest track in phi and z
    TH1F        *fEMCTrkPt;//!Pt tracks with EMCAL cluster
    TH1F        *fEMCTrketa;//!EMC trk eta
    TH1F        *fEMCTrkphi;//!EMC trk phi
    TH2F        *fEMCTPCnsig;//! EMC trk nsig
    TH2F        *fHistNsigEop;//!E/p vs nsig
    TH2F        *fHistEop;//!pt vs E/p
    TH2F        *fM20;//!pt vs M20
    TH2F        *fM02;//!pt vs M20
    TH2F        *fM20EovP;//!M20 vs E/p
    TH2F        *fM02EovP;//!M20 vs E/p
    TH1F        *fInvmassLS;//!Invmass of LS pairs
    TH1F        *fInvmassULS;//!Invmass of ULS pairs
    
    //---Jet histos
    TH1F        *fJetInputTrPhi;//!phi distribution of input tracks to reco jets
    TH1F        *fJetInputTrEta;//!eta distribution
    TH1F        *fBgMedian_Brf;//!median of Bkg before removing leading 2 jets
    TH1F        *fBgArea_Bfr;//!area of Bkg before removing leading 2 jets
    TH1F        *fBgPhi_Bfr;//!phi of Bkg before removing leading 2 jets
    TH1F        *fBgEta_Bfr;//!eta of Bkg before removing leading 2 jets
    TH1F        *fBgPt_Bfr;//!pt of Bkg before removing leading 2 jets
    TH1F        *fBgP_Bfr;//!P of Bkg before removing leading 2 jets
    TH1F        *fBgMedian_Afr;//!median of Bkg after removing leading jets
    TH1F        *fSigArea;//!area of jets
    TH1F        *fSigPt_BfrSub;//!Pt of jets before bkg pt removal
    TH1F        *fSigPt_Sub;//!Pt of jets after bkg pt removal
    TH1F        *fSigPhi;//!phi of jets
    TH1F        *fSigEta;//!eta of jets
    TH1F        *fSigP;//!P of jets
    TH1F        *fSigR;//!Jet R = sqrt(phi2 + eta2)
    TH2F        *fSigEtaPhi;//!Jet eta and phi
    TH1F        *fSigEtainEMC;//!jet eta if in EMC phi acceptance
    TH1F        *fBgArea_Afr;//!area of Bkg after removing leading jets
    TH1F        *fBgPhi_Afr;//!phi of Bkg after removing leading jets
    TH1F        *fBgEta_Afr;//!eta of Bkg after removing leading jets
    TH1F        *fBgPt_Afr;//!Pt of Bkg after removing leading jets
    TH1F        *fBgP_Afr;//!P of Bkg after removing leading jets
    TH1F        *fPtAllElec;//!Pt of electrons
    TH1F        *fNthJetwEle;//!Nth jet with electron
    TH1F        *fPtJetwElec;//!Pt of jet with electron
    TH1F        *fPtElecinJet;//!Pt of electron in jet
    TH1F        *fPtElecLeadInJet;//!Pt of ele, if ele is the lead trk in jet
    TH2F        *fEtaPhiJetwElec;//!Eta and Phi of jets with elec
    TH1F        *fEtaJetwElecinEMC;//!Eta of jets with elec, if in EMC phi acceptance
    TH1F        *fNthJetwPEle;//!Nth jet with photonic electron
    TH1F        *fPtJetwPElec;//!Pt of jet with photonic electron
    TH1F        *fPtPElecinJet;//!Pt of photonic electron in jet
    TH1F        *fPtPElecLeadInJet;//!Pt of photonic ele, if ele is the lead trk in jet
    TH2F        *fEtaPhiJetwPElec;//!Eta and Phi of jets with photonic elec
    TH1F        *fEtaJetwPElecinEMC;//!Eta of jets with photonic elec, if in EMC phi acceptance
    TH1F        *fPtLeadJetwElec;//!Pt of lead jet with electron
    TH1F        *fPtElecinLeadJet;//!Pt of electron if in lead jet
    TH1F        *fPtElecLeadInLeadJet;//!Pt of electron, if elec is lead track in lead jet
    TH1F        *fPtTrigJetwInclE;//!Pt of selected trig jet containing inclusive elec
    TH1F        *fPtTrigJetwSemiInclE;//!Pt of selected trig jet containing photonic elec
    TH1F        *fPtTrigJetwPhotoE;//!Pt of selected trig jet containing semi-inclusive elec

    THnSparse   *fSprsInclEleDphiDeta;//!Sparse of Dphi and Deta for inclusive elec
    THnSparse   *fSprsPhotEleDphiDeta;//!Sparse of Dphi and Deta for photonic elec
    THnSparse   *fSprsSemiIncEleDphiDeta;//!Sparse of Dphi and Deta for semi-inclusive elec
    
    
    AliAnalysisTaskHFEDiJets(const AliAnalysisTaskHFEDiJets&); //
    AliAnalysisTaskHFEDiJets& operator=(const AliAnalysisTaskHFEDiJets&); //
    
    ClassDef(AliAnalysisTaskHFEDiJets, 1); //
};

#endif
