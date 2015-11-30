/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
// QA class of Heavy Flavor quark and fragmeted/decayed particles
// -Check kinematics of Heavy Quarks/hadrons, and decayed leptons
//    pT, rapidity
//    decay lepton kinematics w/wo acceptance
//    heavy hadron decay length, electron pT fraction carried from decay 
// -Check yield of Heavy Quarks/hadrons
//    Number of produced heavy quark
//    Number of produced hadron of given pdg code
//

#ifndef ALIHFEMCQA_H
#define ALIHFEMCQA_H

#ifndef ROOT_TObject
//#include <TObject.h>
#endif

class TH1F;
class TH2F;
class TList;
class TParticle;
class TString;
class AliMCEvent;
class AliGenEventHeader;
class AliMCParticle;
class AliAODMCParticle;
class AliHFEcollection;
class TTreeSRedirector;
class TClonesArray;

//________________________________________________________________
class AliHFEmcQA: public TObject {

  public: 
    enum heavyType {kCharm=4, kBeauty=5, kOthers=6, kElectronPDG=11};
    enum qType {kQuark, kantiQuark, kHadron, keHadron, kDeHadron, kElectron, kElectron2nd};
    enum SourceType {
        kDirectCharm=1,  	// electrons from primary charmed hadrons and primary resonance charmed hadrons(primary charmed hadrons, charmed hadrons decaying from the charmed hadron resonances(ex. D*): D->e ) 
        kDirectBeauty=2,  	// electrons from primary beauty hadrons and primary resonance beauty hadrons (primary beauty hadrons, beauty hadrons decaying from the beauty hadron resonances: B->e)
        kBeautyCharm=3, 	// electrons from charmed hadrons decaying from the beauty hadrons (B->D->e)
        kGamma=4,        	// should be obsolete -> please let me know if you see something! 
        kPi0=5, 	        // electrons from p0 Dalitz
        kElse=6, 	        // all the other sources which was not in this enumeration
        kMisID=7, 	        // not the electrons (hadrons)
        kEta=8, 	        // electrons from eta Dalitz
        kOmega=9, 	        // electrons from omega decay (Dalitz and di-electrons)
        kPhi=10, 	        // electrons from phi decay (di-electron)
        kEtaPrime=11, 	        // electrons from eta prime decay (Dalitz and 2charged-pions&di-electrons) 
        kRho0=12,  	        // electrons from rho decay (di-electron)
        kGammaPi0=13, 	        // electrons from photon conversion where the photon originated from pi0
        kGammaEta=14, 	        // electrons from photon conversion where the photon originated from eta
        kGammaOmega=15, 	// electrons from photon conversion where the photon originated from omega
        kGammaPhi=16, 	        // electrons from photon conversion where the photon originated from phi
        kGammaEtaPrime=17, 	// electrons from photon conversion where the photon originated from eta prime
        kGammaRho0=18, 	        // electrons from photon conversion where the photon originated from rho
        kJpsi=19, 	        // electrons from primary J/psi decay
        kB2Jpsi=20, 	        // electrons from J/psi decay where the J/psi originated from the beauty hadrons
        kKe3=21, 	        // Ke3 electrons
        kGammaB2M=22, 	        // electrons from photon conversion from meson where the meson originated from the beauty hadrons
        kGammaD2M=23, 	        // electrons from photon conversion from meson where the meson originated from the charm hadrons
        kGammaM2M=24, 	        // electrons from photon conversion from the light meson decay where the light meson originated from other light meson 
        kB2M=25, 	        // electrons from the meson where the meson originated from the beauty hadrons 
        kD2M=26, 	        // electrons from the meson where the meson originated from the charm hadrons
        kM2M=27, 	        // electrons from the light meson decay where the light meson originated from other light meson 
        kK0L=28, 	        // K0L -> e +X
    	kK0s2P=31,              // electrons from secondary pions from K0s
        kK0l2P=32,              // electrons from secondary pions from K0l
        kLamda2P=33,            // electrons from secondary pions from Lamda
        kSigma2P=34,            // electrons from secondary pions from Sigma
    	kGammaK0s2P=35,         // electrons from photon conversion where the photon originated from secondary pion from K0s
        kGammaK0l2P=36,         // electrons from photon conversion where the photon originated from secondary pion from K0l
        kGammaLamda2P=37,       // electrons from photon conversion where the photon originated from secondary pion from lamda
        kGammaSigma2P=38,       // electrons from photon conversion where the photon originated from secondary pion from sigma
        kK2P=39,                // electrons from secondary pions from K  
        kGammaK2P=40            // electrons from photon conversion where the photon originated from secondary pion from K
    };
    enum ProcessType {
      kPairCreationFromq,  kPairCreationFromg,  kFlavourExitation,  kGluonSplitting, kInitialPartonShower, kLightQuarkShower
    };
    enum{
      kBgPtBins = 44,
	kElecBgSpecies = 9,
	kCentBins = 11,
	kBgLevels = 3
    };


    AliHFEmcQA();
    AliHFEmcQA(const AliHFEmcQA &p); // copy constructor
    AliHFEmcQA &operator=(const AliHFEmcQA &); // assignment operator

    virtual ~AliHFEmcQA();

    TList *GetList() const { return fQAhistos; };
    void PostAnalyze() const;
    void CreatDefaultHistograms(TList * const qaList); // create default histograms  
    void CreateHistograms(const Int_t kquark); // create histograms for mc qa analysis
    void SetMCEvent(AliMCEvent* const mcEvent){fMCEvent = mcEvent;} 
    void SetGenEventHeader(AliGenEventHeader* const mcHeader){fMCHeader=mcHeader;} // set stack pointer
    void SetMCArray(TClonesArray* const mcarry){fMCArray=mcarry;} // set mcarray pointer
    void Init();

    void GetQuarkKine(TParticle *part, Int_t iTrack, const Int_t kquark); // get heavy quark kinematics distribution
    void GetHadronKine(TParticle *part, const Int_t kquark); // get heavy hadron kinematics distribution
    void GetDecayedKine(TParticle *part, const Int_t kquark, const Int_t kdecayed); // get decay electron kinematics distribution
    void GetDecayedKine(AliAODMCParticle *mcpart, const Int_t kquark, Int_t kdecayed); // get decay electron kinematics for AOD 
    void GetMesonKine(); // get meson and its decay electron pt spectra
    void EndOfEventAna(const Int_t kquark); // run analysis which should be done at the end of the event loop
    Int_t GetSource(const TParticle * const mcpart) const; // return source id 
    Int_t GetElecSource(const AliVParticle * const mctrack, Bool_t isElec) const;
    Int_t GetElecSource(const AliVParticle * const mctrack, Bool_t isElec, Double_t &mpt) const;
    Int_t GetElecSource(TParticle * const mcpart, Bool_t isElec) const; // return electron source id 
    Int_t GetElecSource(TParticle * const mcpart, Bool_t isElec, Double_t &mpt) const; // return electron source id 
    Int_t GetElecSource(const AliAODMCParticle * const mcpart, Bool_t isElec) const;
    Int_t GetSource(const AliVParticle * const mcpart) const; // return electron source id for AOD
    Double_t GetWeightFactor(AliMCParticle *mctrack, const Int_t iBgLevel); // return best/lower/upper weighting factor for electron's mother meson
    Double_t GetWeightFactor(const AliAODMCParticle * const mcpart, const Int_t iBgLevel);
    Double_t GetWeightFactorForPrimaries(const AliAODMCParticle * const mcpart, const Int_t iBgLevel);
    Int_t GetWeightCentralityBin(const Float_t percentile) const; //translate the centrality percentile into the centrality bin of the reference weighting histograms for electron background
    void EnableDebugStreamer() { fIsDebugStreamerON = kTRUE;};
    void EnableGetWeightHist() { fGetWeightHist = kTRUE;};

    void SetBackgroundWeightFactor(Double_t *elecBackgroundFactor, Double_t *binLimit);
    void SetContainerStep(Int_t containerStep) { fContainerStep = containerStep;};
    void SetHFEImpactParameters(Double_t hfeimpactR, Double_t hfeimpactnsigmaR) {fHfeImpactR = hfeimpactR; fHfeImpactnsigmaR = hfeimpactnsigmaR; };
    void SetTrkKine(Double_t pt, Double_t eta, Double_t phi) {fRecPt = pt; fRecEta = eta; fRecPhi = phi;};
    void SetITSInfo(Double_t ilyrhit, Double_t ilyrstat) { fLyrhit = ilyrhit; fLyrstat = ilyrstat;};

    void SetCentrality(Int_t centrality) { fCentrality = centrality; };
    void SetPercentrality(Int_t percentrality) { fPerCentrality = percentrality; };//centrality percentile
    void SetPbPb() { fIsPbPb = kTRUE; };
    void SetPP() { fIsPbPb = kFALSE; };
    void SetPPMultiBin() { fIsppMultiBin = kFALSE; };
    Bool_t IsPbPb() const { return fIsPbPb; };
    Bool_t IsPP() const { return !fIsPbPb; };
    Bool_t IsPPMultiBin() const { return fIsppMultiBin; };

  protected:
    Int_t GetMother(const AliVParticle * const track) const;
    void IdentifyMother(Int_t motherlabel, Int_t &motherpdg, Int_t &grandmotherlabel); // 
    void HardScattering(const Int_t kquark, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from hard scattering
    void ReportStrangeness(Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // report if the quark production process is unknown
    Bool_t IsFromInitialShower(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from initial parton shower 
    Bool_t IsFromFinalParton(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from final parton shower

    AliMCEvent* fMCEvent; // mcevent pointer
    AliGenEventHeader* fMCHeader; // mcheader pointer
    TClonesArray *fMCArray; // mc array pointer

    static const Int_t fgkGluon=21; // gluon pdg code
    static const Int_t fgkMaxGener=10; // ancester level wanted to be checked 
    static const Int_t fgkMaxIter=100; // number of iteration to find out matching particle 
    static const Int_t fgkqType=7; // number of particle type to be checked
    static const Int_t fgkEtaRanges=3; // cuts for different eta ranges

    struct AliHists{
      TH1F *fPdgCode; // histogram to store particle pdg code
      TH1F *fPt; // histogram to store pt
      TH1F *fY; // histogram to store rapidity
      TH1F *fEta; // histogram to store eta

			AliHists()
			  : fPdgCode()
			  , fPt()
			  , fY()
			  , fEta()
      {
			  // default constructor
			};
			AliHists(const AliHists & p)
			  : fPdgCode(p.fPdgCode)
			  , fPt(p.fPt)
			  , fY(p.fY)
			  , fEta(p.fEta)
      {
			  // copy constructor
			};
			AliHists &operator=(const AliHists &)
			{
			  // assignment operator, not yet implemented 
			  return *this;
			}
      void FillList(TList *l) const;
    };
    struct AliHistsComm {
      TH1F *fNq; // histogram to store number of quark
      TH1F *fProcessID; // histogram to store process id 
      TH2F *fePtRatio; // fraction of electron pT from D or B hadron
      TH2F *fPtCorr; // pt correlation between e and direct D or B 
      TH2F *fPtCorrDp; // pt correlation between e and direct D+
      TH2F *fPtCorrD0; // pt correlation between e and direct D0
      TH2F *fPtCorrDrest; // pt correlation between e and direct D rest
      TH2F *fDePtRatio; // fraction of D electron pT from B hadron 
      TH2F *feDistance; // distance between electron production point to mother particle 
      TH2F *fDeDistance; // distance between D electron production point to mother particle

      TH2F *fPtCorrDinein; // pt correlation between e and direct D or B 
      TH2F *fPtCorrDineout; // pt correlation between e and direct D or B 
      TH2F *fPtCorrDoutein; // pt correlation between e and direct D or B 
      TH2F *fPtCorrDouteout; // pt correlation between e and direct D or B 
      TH2F *fPtCorrDpDinein; // pt correlation between e and direct D+
      TH2F *fPtCorrDpDineout; // pt correlation between e and direct D+
      TH2F *fPtCorrDpDoutein; // pt correlation between e and direct D+
      TH2F *fPtCorrDpDouteout; // pt correlation between e and direct D+
      TH2F *fPtCorrD0Dinein; // pt correlation between e and direct D0
      TH2F *fPtCorrD0Dineout; // pt correlation between e and direct D0
      TH2F *fPtCorrD0Doutein; // pt correlation between e and direct D0
      TH2F *fPtCorrD0Douteout; // pt correlation between e and direct D0
      TH2F *fPtCorrDrestDinein; // pt correlation between e and direct D rest
      TH2F *fPtCorrDrestDineout; // pt correlation between e and direct D rest
      TH2F *fPtCorrDrestDoutein; // pt correlation between e and direct D rest
      TH2F *fPtCorrDrestDouteout; // pt correlation between e and direct D rest

      TH2F *fEtaCorrD; // eta correlation between e and direct D 
      TH2F *fEtaCorrDp; // eta correlation between e and direct D 
      TH2F *fEtaCorrD0; // eta correlation between e and direct D 
      TH2F *fEtaCorrDrest; // eta correlation between e and direct D 
      TH2F *fEtaCorrGD; // eta correlation between e and direct D 
      TH2F *fEtaCorrGDp; // eta correlation between e and direct D 
      TH2F *fEtaCorrGD0; // eta correlation between e and direct D 
      TH2F *fEtaCorrGDrest; // eta correlation between e and direct D 

      TH2F *fEtaCorrB; // eta correlation between e and B
      TH2F *fEtaCorrGB; // eta correlation between e and B
      TH2F *fPtCorrBinein; // eta correlation between e and B
      TH2F *fPtCorrBineout; // eta correlation between e and B
      TH2F *fPtCorrBoutein; // eta correlation between e and B
      TH2F *fPtCorrBouteout; // eta correlation between e and B

			AliHistsComm()
			  : fNq()
			  , fProcessID()
			  , fePtRatio()
			  , fPtCorr()
			  , fPtCorrDp()
			  , fPtCorrD0()
			  , fPtCorrDrest()
			  , fDePtRatio()
			  , feDistance()
			  , fDeDistance()

                          , fPtCorrDinein()
                          , fPtCorrDineout()
                          , fPtCorrDoutein()
                          , fPtCorrDouteout()
                          , fPtCorrDpDinein()
                          , fPtCorrDpDineout()
                          , fPtCorrDpDoutein()
                          , fPtCorrDpDouteout()
                          , fPtCorrD0Dinein()
                          , fPtCorrD0Dineout()
                          , fPtCorrD0Doutein()
                          , fPtCorrD0Douteout()
                          , fPtCorrDrestDinein()
                          , fPtCorrDrestDineout()
                          , fPtCorrDrestDoutein()
                          , fPtCorrDrestDouteout()

                          , fEtaCorrD()
                          , fEtaCorrDp()
                          , fEtaCorrD0()
                          , fEtaCorrDrest()
                          , fEtaCorrGD()
                          , fEtaCorrGDp()
                          , fEtaCorrGD0()
                          , fEtaCorrGDrest()

                          , fEtaCorrB()
                          , fEtaCorrGB()
                          , fPtCorrBinein()
                          , fPtCorrBineout()
                          , fPtCorrBoutein()
                          , fPtCorrBouteout()
      {  
			  // default constructor
			};
			AliHistsComm(const AliHistsComm & p)
			  : fNq(p.fNq)
			  , fProcessID(p.fProcessID)
			  , fePtRatio(p.fePtRatio)
			  , fPtCorr(p.fPtCorr)
			  , fPtCorrDp(p.fPtCorrDp)
			  , fPtCorrD0(p.fPtCorrD0)
			  , fPtCorrDrest(p.fPtCorrDrest)
			  , fDePtRatio(p.fDePtRatio)
			  , feDistance(p.feDistance)
			  , fDeDistance(p.fDeDistance)

                          , fPtCorrDinein(p.fPtCorrDinein)
                          , fPtCorrDineout(p.fPtCorrDineout)
                          , fPtCorrDoutein(p.fPtCorrDoutein)
                          , fPtCorrDouteout(p.fPtCorrDouteout)
                          , fPtCorrDpDinein(p.fPtCorrDpDinein)
                          , fPtCorrDpDineout(p.fPtCorrDpDineout)
                          , fPtCorrDpDoutein(p.fPtCorrDpDoutein)
                          , fPtCorrDpDouteout(p.fPtCorrDpDouteout)
                          , fPtCorrD0Dinein(p.fPtCorrD0Dinein)
                          , fPtCorrD0Dineout(p.fPtCorrD0Dineout)
                          , fPtCorrD0Doutein(p.fPtCorrD0Doutein)
                          , fPtCorrD0Douteout(p.fPtCorrD0Douteout)
                          , fPtCorrDrestDinein(p.fPtCorrDrestDinein)
                          , fPtCorrDrestDineout(p.fPtCorrDrestDineout)
                          , fPtCorrDrestDoutein(p.fPtCorrDrestDoutein)
                          , fPtCorrDrestDouteout(p.fPtCorrDrestDouteout)

                          , fEtaCorrD(p.fEtaCorrD)
                          , fEtaCorrDp(p.fEtaCorrDp)
                          , fEtaCorrD0(p.fEtaCorrD0)
                          , fEtaCorrDrest(p.fEtaCorrDrest)
                          , fEtaCorrGD(p.fEtaCorrGD)
                          , fEtaCorrGDp(p.fEtaCorrGDp)
                          , fEtaCorrGD0(p.fEtaCorrGD0)
                          , fEtaCorrGDrest(p.fEtaCorrGDrest)

                          , fEtaCorrB(p.fEtaCorrB)
                          , fEtaCorrGB(p.fEtaCorrGB)
                          , fPtCorrBinein(p.fPtCorrBinein)
                          , fPtCorrBineout(p.fPtCorrBineout)
                          , fPtCorrBoutein(p.fPtCorrBoutein)
                          , fPtCorrBouteout(p.fPtCorrBouteout)
      {
			  // copy constructor
			};
			AliHistsComm &operator=(const AliHistsComm &)
			{
			  // assignment operator, not yet implemented 
			  return *this;
			}
      void FillList(TList *l) const;
    };

    AliHists fHist[3][7][6][11]; // struct of histograms to store kinematics of given particles
    AliHistsComm fHistComm[2][6]; // struct of additional histograms of given particles
    TH2F *fhD[9]; // D meson pt,Y spectra

    TList *fQAhistos;           // Container for QA histos
    AliHFEcollection *fMCQACollection;      //! Tasks own QA collection
    TParticle *fHeavyQuark[50]; //! store pointer of heavy flavour quark 
    Int_t fIsHeavy[2]; // count of heavy flavour
    Int_t fNparents; // number of heavy hadrons to be considered
    Int_t fParentSelect[2][7]; // heavy hadron species
    
    /*
     -------------------------------------------------------------------------------------
    fParentSelect[0][0] =  411; //D+  
    fParentSelect[0][1] =  421; //D0
    fParentSelect[0][2] =  431; //Ds+
    fParentSelect[0][3] = 4122; //Lambdac+
    fParentSelect[0][4] = 4132; //Ksic0
    fParentSelect[0][5] = 4232; //Ksic+
    fParentSelect[0][6] = 4332; //OmegaC0

        -------------------------------------------------------------------------------------
    fParentSelect[1][0] =  511; //B0
    fParentSelect[1][1] =  521; //B+
    fParentSelect[1][2] =  531; //Bs0
    fParentSelect[1][3] = 5122; //Lambdab0
    fParentSelect[1][4] = 5132; //Ksib-
    fParentSelect[1][5] = 5232; //Ksib0
    fParentSelect[1][6] = 5332; //Omegab-
    */


    Double_t fElecBackgroundFactor[kBgLevels][kCentBins][kElecBgSpecies][kBgPtBins];     // Electron background factors
    Double_t fBinLimit[kBgPtBins+1];       // Electron background bins

private:
    Int_t              fCentrality;        // Centrality
    Int_t              fPerCentrality;     // Centrality percentile
    Bool_t             fIsPbPb;            // Analysis Type: pp or PbPb
    Bool_t             fIsppMultiBin;      // pp multiplicity bin analysis
    Int_t              fContainerStep;     // step the weighting factor called
    Bool_t             fIsDebugStreamerON; // check if the debugstreamer is on

    Double_t           fRecPt;             //reconstructed pt
    Double_t           fRecEta;            //reconstructed eta
    Double_t           fRecPhi;            //reconstructed phi
    Double_t           fLyrhit;            //its layer hit
    Double_t           fLyrstat;           //its layer status

    Double_t fHfeImpactR;                  //absolute impact parameter R
    Double_t fHfeImpactnsigmaR;            //absolute impact parameter sigma R 
 
    TTreeSRedirector *fTreeStream;         // TreeStream
    Bool_t             fGetWeightHist;     // Write histogram to calculate weights

  ClassDef(AliHFEmcQA,1);
};

#endif
