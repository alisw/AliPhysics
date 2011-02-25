#ifndef ALIHFEMCQA_H
#define ALIHFEMCQA_H

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

/* $Id$ */ 

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
class AliAODMCParticle;

//________________________________________________________________
class AliHFEmcQA: public TObject {

  public: 
    enum heavyType {kCharm=4, kBeauty=5, kOthers=6, kElectronPDG=11};
    enum qType {kQuark, kantiQuark, kHadron, keHadron, kDeHadron, kElectron, kElectron2nd};
    enum SourceType {kDirectCharm=1, kDirectBeauty=2, kBeautyCharm=3, kGamma=4, kPi0=5, kElse=6, kMisID=7, kEta=8, kOmega=9, kPhi=10, kEtaPrime=11, kRho0=12};
    enum ProcessType {
      kPairCreationFromq,  kPairCreationFromg,  kFlavourExitation,  kGluonSplitting, kInitialPartonShower, kLightQuarkShower
    };

    AliHFEmcQA();
    AliHFEmcQA(const AliHFEmcQA &p); // copy constructor
    AliHFEmcQA &operator=(const AliHFEmcQA &); // assignment operator

    virtual ~AliHFEmcQA();

    TList *GetList() const { return fQAhistos; };
    void PostAnalyze() const;
    void CreatDefaultHistograms(TList * const qaList); // create default histograms  
    void CreateHistograms(const Int_t kquark, Int_t icut, TString hnopt=""); // create histograms for mc qa analysis
    void SetMCEvent(AliMCEvent* const mcEvent){fMCEvent = mcEvent;} 
    void SetGenEventHeader(AliGenEventHeader* const mcHeader){fMCHeader=mcHeader;} // set stack pointer
    void SetMCArray(TClonesArray* const mcarry){fMCArray=mcarry;} // set mcarray pointer
    void Init();

    void GetQuarkKine(TParticle *part, Int_t iTrack, const Int_t kquark); // get heavy quark kinematics distribution
    void GetHadronKine(TParticle *part, const Int_t kquark); // get heavy hadron kinematics distribution
    void GetDecayedKine(TParticle *part, const Int_t kquark, const Int_t kdecayed, Int_t icut); // get decay electron kinematics distribution
		void GetDecayedKine(AliAODMCParticle *mcpart, const Int_t kquark, Int_t kdecayed, Int_t icut); // get decay electron kinematics for AOD 
    void EndOfEventAna(const Int_t kquark); // run analysis which should be done at the end of the event loop
		Int_t GetSource(TParticle * const mcpart); // return source id 
		Int_t GetElecSource(TParticle * const mcpart); // return electron source id 
		Int_t GetSource(AliAODMCParticle * const mcpart); // return electron source id for AOD

  protected:
    void IdentifyMother(Int_t motherlabel, Int_t &motherpdg, Int_t &grandmotherlabel); // 
    void HardScattering(const Int_t kquark, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from hard scattering
    void ReportStrangeness(Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // report if the quark production process is unknown
    Bool_t IsFromInitialShower(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from initial parton shower 
    Bool_t IsFromFinalParton(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from final parton shower

    AliMCEvent* fMCEvent; // mcevent pointer
    AliGenEventHeader* fMCHeader; // mcheader pointer
    TClonesArray *fMCArray; // mc array pointer

    static const Int_t fgkGluon; // gluon pdg code
    static const Int_t fgkMaxGener; // ancester level wanted to be checked 
    static const Int_t fgkMaxIter; // number of iteration to find out matching particle 
    static const Int_t fgkqType; // number of particle type to be checked

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
      TH2F *fDePtRatio; // fraction of D electron pT from B hadron 
      TH2F *feDistance; // distance between electron production point to mother particle 
      TH2F *fDeDistance; // distance between D electron production point to mother particle

			AliHistsComm()
			  : fNq()
			  , fProcessID()
			  , fePtRatio()
			  , fPtCorr()
			  , fDePtRatio()
			  , feDistance()
			  , fDeDistance()
      {  
			  // default constructor
			};
			AliHistsComm(const AliHistsComm & p)
			  : fNq(p.fNq)
			  , fProcessID(p.fProcessID)
			  , fePtRatio(p.fePtRatio)
			  , fPtCorr(p.fPtCorr)
			  , fDePtRatio(p.fDePtRatio)
			  , feDistance(p.feDistance)
			  , fDeDistance(p.fDeDistance)
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

    AliHists fHist[3][7][6]; // struct of histograms to store kinematics of given particles
    AliHistsComm fHistComm[2][6]; // struct of additional histograms of given particles
    TH2F *fhD[9]; // D meson pt,Y spectra
    TH1F *fhDLogbin[9]; // D meson pt w diff binning

    TList *fQAhistos;           // Container for QA histos
    TParticle *fHeavyQuark[50]; //! store pointer of heavy flavour quark 
    Int_t fIsHeavy[2]; // count of heavy flavour
    Int_t fNparents; // number of heavy hadrons to be considered
    Int_t fParentSelect[2][7]; // heavy hadron species


  ClassDef(AliHFEmcQA,1);
};

#endif
