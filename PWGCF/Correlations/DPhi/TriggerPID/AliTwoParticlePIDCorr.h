#ifndef ALITWOPARTICLEPIDCORR_H
#define ALITWOPARTICLEPIDCORR_H

#include "THn.h" // in cxx file causes .../THn.h:257: error: conflicting declaration ‘typedef class THnT<float> THnF’


class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TString;
class TList;
//class AliESDtrackCuts;
class TSeqCollection;
class AliPIDResponse;
class AliAODEvent;
class AliAODTrack;
class AliAODVertex;
class AliEventPoolManager;
class TFormula;
//class AliAnalysisUtils;
class LRCParticlePID;
class AliVParticle;


#include <TObject.h> //LRCParticlePID is a derived class from"TObject"
#include "TMath.h"
#include "TNamed.h"
#include "AliUEHist.h"
#include "AliPID.h"
#include "AliAnalysisTask.h"
#include "AliUEHist.h"
#include "TString.h"
#include "AliVParticle.h"
#include "TParticle.h"
#include "AliLog.h"


#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

namespace AliPIDNameSpace {
  
  enum PIDType
  {
    NSigmaTPC = 0,
    NSigmaTOF,
    NSigmaTPCTOF,  // squared sum
    NSigmaPIDType=NSigmaTPCTOF
  };
    
  enum AliDetectorType
  {
    TPC = 0,
    TOF,
    NDetectors
  };
  
  
  enum AliParticleSpecies
  {
    SpPion = 0,
    SpKaon,
    SpProton,
    NSpecies,
    SpUndefined=999
  }; // Particle species used in plotting
  
  
  enum AliCharge
  {
    Posch = 0,
    Negch,
    NCharge
  };
}


using namespace AliPIDNameSpace;

class AliTwoParticlePIDCorr : public AliAnalysisTaskSE {
 public:
    AliTwoParticlePIDCorr();
    AliTwoParticlePIDCorr(const char *name);
    virtual ~AliTwoParticlePIDCorr();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void    doAODevent();
    virtual void    doMCAODevent();
    virtual void     Terminate(Option_t *);

 private:
 //histograms
    TList *fOutput;        //! Output list
    TH1F *fhistcentrality;//!
    TH1F *fEventCounter; //!
    TH2F *fEtaSpectrasso;//!
    TH2F *fphiSpectraasso;//!
    TH1F *fEtaSpectraTrigall;//!
    TH2F* fCentralityCorrelation;  //! centrality vs multiplicity

    TH2F *fHistoTPCdEdx;//!
    TH2F *fHistoTOFbeta;//!
    // TH3F *fHistocentNSigmaTPC;//! nsigma TPC
    // TH3F *fHistocentNSigmaTOF;//! nsigma TOF 
    TH2F *fsame;
    TH2F *fmix;
    TH1F *fdeletasame;//!
    TH1F *fdelphisame;//!
    TH1F *fdeletamixed;//!
    TH1F *fdelphimixed;//!
    TH1F *fdeletamixedproton;//!
    TH1F *fdelphimixedproton;//!
    TH1F *fdeletamixedkaonpion;//!
    TH1F *fdelphimixedkaonpion;//!
        

    TH1F *fEventno[4][10]; //!
    TH1F *fEventnobaryon[4][10]; //!
    TH1F *fEventnomeson[4][10]; //!


    TH1F *fHistQA[16]; //!

    TH2F *fTPCTOFpion2d[9];//!

    TH1F *fhistoassopioncont[5];//!
    TH1F *fhistoassokaoncont[5];//!
    TH1F *fhistoassoprotoncont[5];//!
    TH1F *fhistotrigbaryoncont[4];//!
    TH1F *fhistotrigmesoncont[4];//!


    TH1F *fHistoNSigmaTPCpion[4][9];//! 
    TH1F *fHistoNSigmaTOFpion[4][9];//!
    TH1F *fHistoNSigmaTPCTOFpion[4][9];//!
    TH1F *fhistopionnsigmaTPCMC[4][9];//! 
    TH1F *fhistopionnsigmaTOFMC[4][9];//!
    TH1F *fhistopionnsigmaTPCTOFMC[4][9];//! 
    TH1F *fhistokaonnsigmaTPCMC[4][9];//! 
    TH1F *fhistokaonnsigmaTOFMC[4][9];//!
    TH1F *fhistokaonnsigmaTPCTOFMC[4][9];//! 
    TH1F *fhistoprotonnsigmaTPCMC[4][9];//! 
    TH1F *fhistoprotonnsigmaTOFMC[4][9];//!
    TH1F *fhistoprotonnsigmaTPCTOFMC[4][9];//! 
    TH1F *fhistoelectronnsigmaTPCMC[4][9];//! 
    TH1F *fhistoelectronnsigmaTOFMC[4][9];//!
    TH1F *fhistoelectronnsigmaTPCTOFMC[4][9];//!  

    //after electron removal cut
    // TH1F *fHistoNSigmaTPCpioncut[4][9];//! 
    //TH1F *fHistoNSigmaTOFpioncut[4][9];//!
    
    TH1F *fEtaSpectraTrig[4][10];//!
    TH1F *fEtaSpectraTrigbaryon[4][10];//!
    TH1F *fEtaSpectraTrigmeson[4][10];//! 

    TH2F *falltrigallasso[4][2][10]; //!
    TH2F *falltrigpionasso[4][2][10];//!
    TH2F *falltrigkaonasso[4][2][10];//! 
    TH2F *falltrigprotonasso[4][2][10];//!   
    TH2F *fbaryontrigallasso[4][2][10]; //!
    TH2F *fbaryontrigpionasso[4][2][10]; //!
    TH2F *fbaryontrigkaonasso[4][2][10]; //!
    TH2F *fbaryontrigprotonasso[4][2][10]; //!
    TH2F *fmesontrigallasso[4][2][10]; //!
    TH2F *fmesontrigpionasso[4][2][10]; //!
    TH2F *fmesontrigkaonasso[4][2][10]; //!
    TH2F *fmesontrigprotonasso[4][2][10]; //!


    TH2F *falltrigallassomix[4][2][10]; //!
    TH2F *falltrigpionassomix[4][2][10];//!
    TH2F *falltrigkaonassomix[4][2][10];//! 
    TH2F *falltrigprotonassomix[4][2][10];//! 
    TH2F *fbaryontrigallassomix[4][2][10]; //!
    TH2F *fbaryontrigpionassomix[4][2][10]; //!
    TH2F *fbaryontrigkaonassomix[4][2][10]; //!
    TH2F *fbaryontrigprotonassomix[4][2][10]; //!
    TH2F *fmesontrigallassomix[4][2][10]; //!
    TH2F *fmesontrigpionassomix[4][2][10]; //!
    TH2F *fmesontrigkaonassomix[4][2][10]; //!
    TH2F *fmesontrigprotonassomix[4][2][10]; //!


    THnF *fTHnrecoallPid[6];//0 pion, 1 kaon,2 proton,3 others,4 mesons,5 all
    THnF *fTHngenprimPidTruth[6];
    THnF *effcorection[6];
    // THnF *effmap[6];





    TH1F *recoallpt[4];//!
    TH1F *recoalleta[4];//!
    TH1F *alltrigeta[4];//!
    TH1F *allassoeta[4];//!
    TH1F *baryontrigeta[4];//!
    TH1F *mesontrigeta[4];//!
    TH1F *pionassoeta[4];//!
    TH1F *kaonassoeta[4];//!
    TH1F *protonassoeta[4];//!
    TH1F *recoallphi[4];//!
    TH1F *MCrecomatchedprimpt[4];//!
    TH1F *MCrecomatchedprimeta[4];//!
    TH1F *MCrecomatchedprimphi[4];//!
    TH1F *MCtruthpt[4];//! 
    TH1F *MCtrutheta[4];//! 
    TH1F *MCtruthphi[4];//!

    TH1F *MCrecomatchedprimpionpt[4];//!
    TH1F *MCrecomatchedprimpioneta[4];//!
    TH1F *MCrecomatchedprimpionphi[4];//!

    TH1F *MCrecomatchedprimkaonpt[4];//!
    TH1F *MCrecomatchedprimkaoneta[4];//!
    TH1F *MCrecomatchedprimkaonphi[4];//!

    TH1F *MCrecomatchedprimprotonpt[4];//!
    TH1F *MCrecomatchedprimprotoneta[4];//!
    TH1F *MCrecomatchedprimprotonphi[4];//!

    TH1F *MCtruthpionpt[4];//!
    TH1F *MCtruthpioneta[4];//!
    TH1F *MCtruthpionphi[4];//!

    TH1F *MCtruthkaonpt[4];//!
    TH1F *MCtruthkaoneta[4];//!
    TH1F *MCtruthkaonphi[4];//!

    TH1F *MCtruthprotonpt[4];//!
    TH1F *MCtruthprotoneta[4];//!
    TH1F *MCtruthprotonphi[4];//! 
   
      
    //TObjArray* trackstrig;
    // TObjArray* tracksasso;
    //TH2F* fControlConvResoncances; //! control histograms for cuts on conversions and resonances

    Int_t ClassifyTrack(AliAODTrack* track,Int_t centbin,AliAODVertex* vertex,Float_t magfield,Bool_t dcacut);
 Int_t Getzbin(Float_t z);
 Int_t Getcentbin(Float_t cent_v0m);
 Int_t Getptbin(Float_t pt);

 

 void Fillcorrelation(TObjArray *trackstrig,TObjArray *tracksasso,Int_t centbin,Int_t vtx,Float_t bSign,Bool_t twoTrackEfficiencyCut,Bool_t mixcase);//mixcase=kTRUE in case of mixing
 Float_t GetTrackbyTrackeffvalue(AliAODTrack* track,Float_t cent,Float_t evzvtx, Int_t parpid);

//Mixing functions
  void DefineEventPool();
  // AliAnalysisUtils *fUtils;//!
  AliEventPoolManager    *fPoolMgr;//! 
  TClonesArray          *fArrayMC;//!
  TString          fAnalysisType;//!          // "MC", "ESD", "AOD"


    //PID part histograms

  //PID functions
    Bool_t HasTPCPID(AliAODTrack *track) const; // has TPC PID
    Bool_t HasTOFPID(AliAODTrack *track) const; // has TOF PID
    Float_t GetBeta(AliAODTrack *track);
    void CalculateNSigmas(AliAODTrack *track,Int_t centbin);
    Int_t FindMinNSigma(AliAODTrack *track);
    Bool_t* GetDoubleCounting(AliAODTrack * trk);
    Int_t GetParticle(AliAODTrack * trk,Int_t centbin);  
   
   
     	   
   Float_t twoTrackEfficiencyCutValue;
  //Pid objects
  AliPIDResponse *fPID; //! PID
  Int_t eventno;
  Float_t fPtTOFPID; //lower pt bound for the TOF pid
  Bool_t fRequestTOFPID;//if true returns kSpUndefined if the TOF signal is missing
  PIDType fPIDType; // PID type  Double_t fNSigmaPID; // number of sigma for PID cut
  Double_t fNSigmaPID; // number of sigma for PID cut
  Double_t fNSigmaPIDtrig1; // number of sigma for PID cut
  Bool_t fUseExclusiveNSigma;//if true returns the identity only if no double counting(i.e not in the overlap area)
  Bool_t fRemoveTracksT0Fill;//if true remove tracks for which only StartTime from To-Fill is available (worst resolution)
 Int_t fSelectCharge;           // (un)like sign selection when building correlations: 0: no selection; 1: unlike sign; 2: like sign
    Int_t fTriggerSelectCharge;    // select charge of trigger particle: 1: positive; -1 negative
    Int_t fAssociatedSelectCharge; // select charge of associated particle: 1: positive; -1 negative
    Float_t fTriggerRestrictEta;   // restrict eta range for trigger particle (default: -1 [off])
    Bool_t fEtaOrdering;           // eta ordering, see AliUEHistograms.h for documentation
    Bool_t fCutConversions;        // cut on conversions (inv mass)
    Bool_t fCutResonances;         // cut on resonances (inv mass)
    Int_t fRejectResonanceDaughters; // reject all daughters of all resonance candidates (1: test method (cut at m_inv=0.9); 2: k0; 3: lambda)
 Int_t fOnlyOneEtaSide;       // decides that only trigger particle from one eta side are considered (0 = all; -1 = negative, 1 = positive)
 Bool_t fPtOrder;		   // apply pT,a < pt,t condition; default: kTRUE
    Bool_t		fInjectedSignals;	  // check header to skip injected signals in MC
    Bool_t fRemoveWeakDecays;	   // remove secondaries from weak decays from tracks and particles
    Bool_t fRemoveDuplicates;// remove particles with the same label (double reconstruction)
    Bool_t applyefficiency;//if kTRUE then eff correction calculation starts
  TFormula*      fDCAXYCut;          // additional pt dependent cut on DCA XY (only for AOD)


  Float_t fnsigmas[NSpecies][NSigmaPIDType+1]; //nsigma values
  Bool_t fHasDoubleCounting[NSpecies];//array with compatible identities

  //Int_t fPIDMethod; // PID method

 //functions
  Float_t PhiRange(Float_t DPhi);
  Float_t GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
Float_t GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);
	   
    
    AliTwoParticlePIDCorr(const AliTwoParticlePIDCorr&); // not implemented
    AliTwoParticlePIDCorr& operator=(const AliTwoParticlePIDCorr&); // not implemented
    
    ClassDef(AliTwoParticlePIDCorr, 1); // example of analysis
};
class LRCParticlePID : public TObject {
public:
 LRCParticlePID(Int_t par,Short_t icharge,Float_t pt,Float_t eta, Float_t phi, Int_t cent, Int_t zvtx,Float_t effcorrectionval)
   :fparticle(par),fcharge(icharge),fPt(pt), fEta(eta), fPhi(phi),fcent(cent),fzvtx(zvtx),feffcorrectionval(effcorrectionval)  {}
  virtual ~LRCParticlePID() {}

  
    virtual Float_t Eta()        const { return fEta; }
    virtual Float_t Phi()        const { return fPhi; }
    virtual  Int_t getcent() const {return fcent;}
    virtual Int_t getzvtx() const {return fzvtx;}
    virtual Float_t Pt() const { return fPt; }
    Int_t getparticle() const {return fparticle;}
    virtual Short_t Charge()      const { return fcharge; }
    Float_t geteffcorrectionval() const {return feffcorrectionval;}
    virtual Bool_t IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }


private:
  LRCParticlePID(const LRCParticlePID&);  // not implemented
   LRCParticlePID& operator=(const LRCParticlePID&);  // not implemented
  //Double_t fcent;
  //Double_t fzvtx;
  //Int_t feventno;
  Int_t fparticle;
  Short_t fcharge;
  Float_t fPt;
  Float_t fEta;
  Float_t fPhi;
  Int_t fcent;
  Int_t fzvtx;
  Float_t feffcorrectionval;
  ClassDef(LRCParticlePID, 1);
} ;

#endif

