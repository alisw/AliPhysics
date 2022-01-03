/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Sigma Plus Class                                                      //
//                                                                        //
//  This class is used for the reconstruction of Sigma+ baryons in the    //
//  Sigma+ -> p+ + pi0 decay channel.                                     //
//                                                                        //
//  Author: B.Heybeck (b.heybeck@cern.ch)                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskSigmaPlus_H
#define AliAnalysisTaskSigmaPlus_H

#include "AliAnalysisTaskSE.h"

// Includes needed for KFParticle
#define HomogeneousField ///Homogenous Magnetic Field in z-Direction
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPTrack.h"

class AliMCEvent;
class AliMCEventHandler;
class AliAODVZERO;

class AliPIDResponse;    //forward declaration of PID Response Object

class KFParticleCD : public KFParticle 
{
    //This class inherits from KFParticle but adds the "CheckDaughter" function 
    //to avoid floating exceptions in KFParticle. Checks "GetMeasurement" which
    //is private in KF. Many Thanks to M. Hartung 
    public:
    Bool_t CheckDaughter(KFParticle daughter) {
        Float_t m[8], mV[36], D[3][3];
        Bool_t CheckMeasurement = KFParticleBase::GetMeasurement(daughter, m, mV, D);
        if(!CheckMeasurement) {if(fWarnings) AliWarning("Checking KF Daughter failed: GetMeasurement() returned false"); return kFALSE;}
        //Require the Covariance Matrix between Daughters to have a reasonable value. 10^6 exceedingly large? Maybe, but it works fine.
        for(Int_t i=0;i<3;i++){for(Int_t j=0;j<3;j++) {if(TMath::Abs(D[i][j])>1000000) {if(fWarnings) AliWarning("Checking KF Daughter failed: Covariance limit exceeded"); return kFALSE;}}}
        return kTRUE;  
    }
    void ActivateWarnings() {fWarnings = kTRUE;} //Activate Warnings to check how much is discarded. Caution: Very anoying!
    protected:
    Bool_t fWarnings = kFALSE;
};

class AliAnalysisTaskSigmaPlus : public AliAnalysisTaskSE  
{
    public:
                            AliAnalysisTaskSigmaPlus();
                            AliAnalysisTaskSigmaPlus(const char *name);
        virtual             ~AliAnalysisTaskSigmaPlus();
                            
        virtual void        UserCreateOutputObjects();
        virtual void        UserExec(Option_t* option);
        virtual void        Terminate(Option_t* option);
        
    private:

        //Functions for Filling Histograms
        void FillHistogram(const char* key,Double_t x) const;                         // Fill Histogram with the name "key"
        void FillHistogram(const char* key,Double_t x, Double_t y) const;             // Histogram must be in TList "fOutputList"
        void FillHistogram(const char* key,Double_t x, Double_t y, Double_t z) const; // 1, 2 and 3 Dimensions is supported (overload)

        //Functions for Event Processing
        void FillProtonArray();                           // Select (Anti-)Protons and fill array (fProtonArray)
        void ProcessMCParticles() const;                  // Fill MC Histograms
        void FillV0PhotonArray();                         // Select V0 Photons and fill arrays (fConvPhotonArray)
        void FindAddPhotons();                            // Find track pairs which could be Photons but which were not found by the V0 Finders 
        void FillElectronArray();                         // Select Electrons/Positrons and fill array (fElectronArray)
        void ReconstructParticles();                      // Reconstruct Pi0 and Sigma Plus from selected Protons and Photons  
        void ReconstructParticlesOff();                   // Reconstruct Pi0 and Sigma Plus from one Finder Photon and one additional one  
        void ReconstructParticlesOff2();                  // Reconstruct Pi0 and Sigma Plus only from additional Photons  
        void ReconstructParticlesOneGamma();              // Reconstruct Pi0 and Sigma Plus in an "exotic" decay channel with only one photon!  
        AliESDv0 Tracks2V0vertex(Int_t v0index) const;    // Build an ESDv0 from track pairs found by FindAddPhotons(). Uses ESDv0::refit() 

        TList*                fOutputList;                //!<! Output list which contains all Histograms

        AliAODEvent*          aodEvent;                   //!<! AOD Event to be processed
        AliMCEvent*           mcEvent;                    //!<! MC Event to be processed   
        TClonesArray*         AODMCTrackArray;            //!<! TClonesArray containing the MC Particles  

        AliPIDResponse*       fPIDResponse;               //!<! PID response object 

        TTree*                fSigmaCandTree;             //! Tree with Sigma candidates
        TTree*                fSigmaCandTreeExtra;        //! Tree with Extra Sigma candidates

        //V0 Maps
        std::map<std::pair<int, int>, int> fONMap;        //!<! Map the track pairs found by the V0
        std::map<std::pair<int, int>, int> fOFFMap;       //!<! Finders to avoid double counting. Uses ESD IDs!

        //Particle Arrays
        std::vector<int>      fProtonArray;               //!<! Proton candidates
        std::vector<int>      fElectronArray;             //!<! Electron candidates
        std::vector<int>      fConvPhotonArray;           //!<! V0 Photon candidates
        std::vector<std::pair<int,int>>  PairIndexArray;  //!<! V0 Photon candidates
        std::vector<int>      fV0ParticleIDArray;         //!<! Save IDs of Particles found by any of the V0 Finders

        Bool_t                isMonteCarlo;               //True if MC information is available  

        Double_t              cElectronMass;              //Physical constants
        Double_t              cProtonMass;                //Physical constants
        Double_t              cPionMass;                  //Physical constants
        Double_t              cPi0Mass;                   //Physical constants
        Double_t              c;                          //Physical constants

        Double_t              Bz;                         //Magnetic Field in the Event
        Double_t              primaryVtxPosX;             //Position of the primary Vertex of the Event
        Double_t              primaryVtxPosY;             //Position of the primary Vertex of the Event
        Double_t              primaryVtxPosZ;             //Position of the primary Vertex of the Event
        Int_t                 nTracks;                    //Number of Tracks in the Event
        Int_t                 Centrality;                 //Centrality of the Event

        //Activate some extra output. Caution: Very anoying!
        Bool_t                fDebug = kFALSE;
        void SetActivateDebug(Bool_t debug) {fDebug = debug;}

        //Activate and deactivate parts of the Analysis 
        Bool_t                fProcessProtons = kTRUE; 
        Bool_t                fProcessMCParticles = kTRUE; 
        Bool_t                fProcessV0s = kTRUE; 
        Bool_t                fProcessAddPhoton = kTRUE; 
        Bool_t                fProcessElectrons = kTRUE; 
        Bool_t                fProcessReco = kTRUE; 
        Bool_t                fProcessRecoOff = kTRUE; 
        Bool_t                fProcessRecoOff2 = kTRUE; 
        Bool_t                fProcessOneGamma = kFALSE; 
        void SetActivateProton(Bool_t processprotons) {fProcessProtons = processprotons;}
        void SetActivateMCParticles(Bool_t processmcpart) {fProcessMCParticles = processmcpart;}
        void SetActivateV0s(Bool_t processv0s) {fProcessV0s = processv0s;}
        void SetActivateAddPhoton(Bool_t processaddphoton) {fProcessAddPhoton = processaddphoton;}
        void SetActivateElectron(Bool_t processelectrons) {fProcessElectrons = processelectrons;}
        void SetActivateReco(Bool_t processreco) {fProcessReco = processreco;}
        void SetActivateRecoOff(Bool_t processrecooff) {fProcessRecoOff = processrecooff;}
        void SetActivateRecoOff2(Bool_t processrecooff2) {fProcessRecoOff2 = processrecooff2;}
        void SetActivateOneGamma(Bool_t processonegamma) {fProcessOneGamma = processonegamma;}

        //FillProtonArray Cuts
        Double_t              fMaxProtEta = 0.9;    
        Double_t              fMinTPCClustProt = 70;
        Double_t              fMaxNsigProtTPC = 3;
        Double_t              fMaxNsigProtTOF = 3;  
        Double_t              fMaxpOnlyTPCPID = 0.75;
        Double_t              fMinProtpt = 0; 
        Double_t              fMaxProtpt = 15;  
        void SetProtonMaxEta(Double_t maxeta) {fMaxProtEta = maxeta;}
        void SetProtonMinTPCCluster(Double_t minclst) {fMinTPCClustProt = minclst;}
        void SetProtonMaxNSigmaTPC(Double_t nsigtpc) {fMaxNsigProtTPC = nsigtpc;}
        void SetProtonMaxNSigmaTOF(Double_t nsigtof) {fMaxNsigProtTOF = nsigtof;}
        void SetProtonMaxpOnlyTPCPID(Double_t pmaxonlytpc) {fMaxpOnlyTPCPID = pmaxonlytpc;}
        void SetProtonMinpt(Double_t minpt) {fMinProtpt = minpt;}
        void SetProtonMaxpt(Double_t maxpt) {fMaxProtpt = maxpt;}

        //ProcessMCParticles Cuts
        Double_t              fMaxMCEta = 0.9;    
        void SetMCMaxEta(Double_t maxeta) {fMaxMCEta = maxeta;}

        //FillV0PhotonArray Cuts
        Double_t              fMaxDaughtEta = 0.9;    
        Double_t              fMinTPCClustDaught = 40;
        Double_t              fMaxNsigDaughtTPC = 3;
        Double_t              fMaxalpha = 0.9;
        Double_t              fMaxqt = 0.04;
        Double_t              fMaxopenangle = 0.2;
        Double_t              fMaxdeltatheta = 0.1;
        Double_t              fMinV0CPA = 0.85;
        Double_t              fMinV0Radius = 3; 
        Double_t              fMaxV0Radius = 200;
        Double_t              fMaxphotonmass = 0.05;
        void SetV0DaughtMaxEta(Double_t maxeta) {fMaxDaughtEta = maxeta;}
        void SetV0DaughtMinTPCCluster(Double_t minclst) {fMinTPCClustDaught = minclst;}
        void SetV0DaughtMaxNSigmaTPC(Double_t nsigtpc) {fMaxNsigDaughtTPC = nsigtpc;}
        void SetV0DaughtMaxAlpha(Double_t maxalpha) {fMaxalpha = maxalpha;}
        void SetV0DaughtMaxQt(Double_t maxqt) {fMaxqt = maxqt;}
        void SetV0DaughtMaxAngle(Double_t maxangle) {fMaxopenangle = maxangle;}
        void SetV0DaughtMaxTheta(Double_t maxtheta) {fMaxdeltatheta = maxtheta;}
        void SetAddPhotonMinCPA(Double_t mincpa) {fMinV0CPA = mincpa;}
        void SetAddPhotonMinConRad(Double_t minrad) {fMinV0Radius = minrad;}
        void SetAddPhotonMaxConRad(Double_t maxrad) {fMaxV0Radius = maxrad;}
        void SetV0DaughtMaxMass(Double_t maxmass) {fMaxphotonmass = maxmass;}

        //FindAddPhotons Cuts. CAUTION: Uses same Cuts for Daughters as used for V0s 
        Double_t              fMinDCADaughtPV = 0.1;
        Double_t              fMaxDCADaught = 5;
        void SetAddPhotonMinDCAPV(Double_t mindca) {fMinDCADaughtPV = mindca;}
        void SetAddPhotonMaxDaughtDCA(Double_t maxdca) {fMaxDCADaught = maxdca;}

        //FillElectronArray Cuts
        Double_t              fMaxElecEta = 0.9;    
        Double_t              fMinTPCClustElec = 40;
        Double_t              fMaxNsigElecTPC = 3;
        Double_t              fMinNsigHadronTPC = 0; //If 0, no Hadron rejection or TOF Cut is made
        Double_t              fMaxNsigElecTOF = 3; 
        Double_t              fMaxElecpt = 5; 
        void SetElectronMaxEta(Double_t maxeta) {fMaxElecEta = maxeta;}
        void SetElectronMinTPCCluster(Double_t minclst) {fMinTPCClustElec = minclst;}
        void SetElectronMaxNSigmaTPC(Double_t nsigtpc) {fMaxNsigElecTPC = nsigtpc;}
        void SetElecMinNSigHadronTPC(Double_t nsighadtpc) {fMinNsigHadronTPC = nsighadtpc;}
        void SetElectronMaxNSigmaTOF(Double_t nsigtof) {fMaxNsigElecTOF = nsigtof;}
        void SetElectronMaxpt(Double_t maxpt) {fMaxElecpt = maxpt;}

        //Fill Sigma Tree Cuts     (Coarse Cuts to reduce Tree Size)
        Double_t              fRejectGammaAutoCorr = 0.02;
        Double_t              fMinPi0Mass = 0.06;    
        Double_t              fMaxPi0Mass = 0.19;    
        Double_t              fMaxSigmaPA = 0.2;     
        Double_t              fMaxSigmaMass = 1.7;    
        Double_t              fMinProtonDCAxy = 0;    
        Double_t              fMinProtonDCAz = 0;    
        Double_t              flowkstar = 0.6;    
        Double_t              fverylowkstar = 0.3;    
        Double_t              fveryverylowkstar = 0.1;    
        void SetRejectGammaAutoCorr(Double_t minrejectmass) {fRejectGammaAutoCorr = minrejectmass;}
        void SetPi0MinMass(Double_t minpi0mass) {fMinPi0Mass = minpi0mass;}
        void SetPi0MaxMass(Double_t maxpi0mass) {fMaxPi0Mass = maxpi0mass;}
        void SetSigmaMaxPA(Double_t maxpa) {fMaxSigmaPA = maxpa;}
        void SetSigmaMaxMass(Double_t maxsigmass) {fMaxSigmaMass = maxsigmass;}
        void SetProtonMinDCAxy(Double_t minprotondcaxy) {fMinProtonDCAxy = minprotondcaxy;}
        void SetProtonMinDCAz(Double_t minprotondcaz) {fMinProtonDCAz = minprotondcaz;}
        void SetLowkstar(Double_t lowkstar) {flowkstar = lowkstar;}
        void SetVeryLowkstar(Double_t verylowkstar) {fverylowkstar = verylowkstar;}
        void SetVeryVeryLowkstar(Double_t veryverylowkstar) {fveryverylowkstar = veryverylowkstar;}

        //Branches of TTree "fSigmaCandTree"
        Bool_t                fIsMCSigma;
        Bool_t                fIsMCPrimary;
        Bool_t                fIsV01fromFinder;
        Bool_t                fIsV02fromFinder;
        Bool_t                fIsV01Onthefly;
        Bool_t                fIsV02Onthefly;
        Int_t                 fSigRunnumber;               
        ULong_t               fSigTriggerMask;               
        Float_t               fSigCentrality;               
        Float_t               fSigBField;               
        Float_t               fInvSigMass;               
        Float_t               fSigPA;               
        Float_t               fSigCharge;        
        Float_t               fSigPx;        
        Float_t               fSigPy;        
        Float_t               fSigPz;        
        Float_t               fPrimVertX;        
        Float_t               fPrimVertY;        
        Float_t               fPrimVertZ;        
        Float_t               fSigDecayVertX;        
        Float_t               fSigDecayVertY;        
        Float_t               fSigDecayVertZ;        
        Float_t               fPhoton1Radius;        
        Float_t               fPhoton2Radius;        
        Float_t               fPhoton1DCAPV;        
        Float_t               fPhoton2DCAPV;        
        Float_t               fInvPi0Mass;        
        Float_t               fPi0Px;        
        Float_t               fPi0Py;        
        Float_t               fPi0Pz;
        Float_t               fPi0DecayVertX;        
        Float_t               fPi0DecayVertY;        
        Float_t               fPi0DecayVertZ;
        Float_t               fPi0PhotPhotDCA;
        Float_t               fProtonPx;        
        Float_t               fProtonPy;        
        Float_t               fProtonPz;        
        Float_t               fProtonDCAtoPVxy;        
        Float_t               fProtonDCAtoPVz;        
        Float_t               fProtonPi0DCA;        
        Short_t               fnPair;        
        Short_t               fnPairlowkstar;        
        Short_t               fnPairverylowkstar;        
        Short_t               fnPairveryverylowkstar;        
        //End of Branches of "fSigmaCandTree"

        //Extra Branches of TTree "fSigmaCandTreeExtra"
        Bool_t                fIsV0fromFinder;
        Bool_t                fIsV0Onthefly;
        Float_t               fPhotonPx;        
        Float_t               fPhotonPy;        
        Float_t               fPhotonPz;
        Float_t               fPhotonRadius;        
        Float_t               fPhotonDCAPV;        
        Float_t               fExtPhotProtDCA;
        //End of Extra Branches of "fSigmaCandTreeExtra"

        AliAnalysisTaskSigmaPlus(const AliAnalysisTaskSigmaPlus&); // not implemented
        AliAnalysisTaskSigmaPlus& operator=(const AliAnalysisTaskSigmaPlus&); // not implemented

        ClassDef(AliAnalysisTaskSigmaPlus, 1);
};

#endif
