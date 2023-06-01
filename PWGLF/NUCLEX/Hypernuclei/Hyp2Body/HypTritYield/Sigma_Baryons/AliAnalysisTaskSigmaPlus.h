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

class AliEventPoolManager;

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
        void FillProtonTree();                            // Fill Tree of Protons (fProtonTree)
        void ProcessMCParticles() const;                  // Fill MC Histograms
        void FillV0PhotonArray();                         // Select V0 Photons and fill arrays (fConvPhotonArray)
        void FillCaloClusterArray();                      // Select Calo Photons and fill arrays (fCaloPhotonArray)
        void ReconstructParticles();                      // Reconstruct Pi0 and Sigma Plus from selected Protons and Photons  
        void ReconstructParticlesPHOS();                  // Reconstruct Pi0 and Sigma Plus from selected Protons and Photons where 1 Photons comes from PHOS  
        void ReconstructParticlesCalc();                  // Reconstruct Pi0 and Sigma Plus from selected Protons and Photons where 1 Photon is calculated  

        TList*                fOutputList;                //!<! Output list which contains all Histograms

        AliAODEvent*          aodEvent;                   //!<! AOD Event to be processed
        AliMCEvent*           mcEvent;                    //!<! MC Event to be processed   

        TClonesArray*         AODMCTrackArray;            //!<! TClonesArray containing the MC Particles  
        AliPIDResponse*       fPIDResponse;               //!<! PID response object 

        AliEventPoolManager*  fEvPoolMgr;                 //!<! Event Pool Manager for Event-Mixing - Contains Protons for PCM-PCM
        AliEventPoolManager*  fEvPoolMgr2;                //!<! Event Pool Manager for Event-Mixing - Contains primary Protons for Femptoscopy
        AliEventPoolManager*  fEvPoolMgr3;                //!<! Event Pool Manager for Event-Mixing - Contains Clusters for PCM-Calo

        //PCM PCM Trees
        TTree*                fSigmaCandTree;             //!<! Tree with Sigma candidates
        TTree*                fSigmaPairTreeSE;           //!<! Tree with Sigma Proton Pairs Same Event
        TTree*                fSigmaPairTreeME;           //!<! Tree with Sigma Proton Pairs Mixed Event
        TTree*                fSigmaMEBackgroundTree;     //!<! Tree with Sigma candidate Mixed Event Background
        
        //PCM Calo Trees                                     
        TTree*                fSigmaPHOSCandTree;         //!<! Tree with Sigma candidates from PHOS
        TTree*                fSigmaPairTreePHOSSE;       //!<! Tree with Sigma Proton Pairs Same Event
        TTree*                fSigmaPairTreePHOSME;       //!<! Tree with Sigma Proton Pairs Mixed Event
        TTree*                fSigmaPHOSMEBkgTree;        //!<! Tree with Sigma candidates Mixed Event Background

        //Test
        TTree*                fSigmaCalcCandTree;         //!<! Tree with Sigma candidates calculated

        //Proton Tree                                     
        TTree*                fProtonTree;                //!<! Tree with Protons

        //Resolution functions
        TF1*                  fThetaFunc;                 //!<!   
        TF1*                  fPhiFunc;                   //!<!           

        //V0 Arrays
        std::vector<int>      fOnFlyVector;               //!<! Save the track IDs used by the V0 Finders 
        std::vector<int>      fFinderVector;              //!<! to avoid double counting. Uses ESD IDs!
        std::vector<int>      fV0ParticleIDArray;         //!<! Save IDs of Particles found by any of the V0 Finders

        //Particle Arrays
        std::vector<int>      fProtonArray;               //!<! Proton candidates
        std::vector<int>      fProtonArray2;              //!<! Stricter Selection for Mixing
        std::vector<int>      fConvPhotonArray;           //!<! V0 Photon candidates
        std::vector<int>      fCaloPhotonArray;           //!<! Calo Photon candidates

        Bool_t                isMonteCarlo;               //True if MC information is available  

        Double_t              cElectronMass;              //Physical constants
        Double_t              cProtonMass;                //Physical constants
        Double_t              cSigmaMass;                 //Physical constants
        Double_t              cPi0Mass;                   //Physical constants
        Double_t              c;                          //Physical constants

        Double_t              Bz;                         //Magnetic Field in the Event
        Double_t              primaryVtxPosX;             //Position of the primary Vertex of the Event
        Double_t              primaryVtxPosY;             //Position of the primary Vertex of the Event
        Double_t              primaryVtxPosZ;             //Position of the primary Vertex of the Event
        Double_t              primaryVtxPosXMC;           //Position of the MC primary Vertex of the Event
        Double_t              primaryVtxPosYMC;           //Position of the MC primary Vertex of the Event
        Double_t              primaryVtxPosZMC;           //Position of the MC primary Vertex of the Event
        Int_t                 nTracks;                    //Number of Tracks in the Event

        Double_t              Centrality;                 //Centrality of the Event
        UInt_t                EventTriggers;              //Triggers of the Event
        Short_t               fRefMultComb05;             //combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.5
        Short_t               fRefMultComb08;             //combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8
        Short_t               fRefMultComb10;             //combined reference multiplicity (tracklets + ITSTPC) in |eta|<1.0
        ULong64_t             fGlobalEventID;             //Global ID of the Collision  

        Double_t              fMeanProtonpt;              //Mean pt of protons for correlations
        Bool_t                fEventhasSigma;             //Event has Sigma Candidate(s)
        Bool_t                fEventhasProton;            //Event has Proton Candidate(s)

        Bool_t                fEventhasSigmaCand;             //Event has Sigma Candidate(s)

        //Remove generated pile-up events (included in pass2/pass3)
        Bool_t                fRemoveGenPileup = kTRUE;   

        //Activate some extra output. Caution: Very anoying!
        Bool_t                fDebug = kFALSE;

        //Activate and deactivate parts of the Analysis 
        Bool_t                fProcessMCParticles = kTRUE; 
        Bool_t                fProcessProtons = kTRUE;
        Bool_t                fSaveProtons = kTRUE;
        Bool_t                fProcessV0s = kTRUE; 
        Bool_t                fProcessReco = kTRUE; 
        Bool_t                fSavePartCand = kTRUE; 
        Bool_t                fFillPairTreeSE = kTRUE; 
        Bool_t                fFillPairTreeME = kTRUE; 
        Bool_t                fSaveMixedBackground = kTRUE; 
        Bool_t                fProcessClusters = kTRUE;
        Bool_t                fMapClusterstoTracks = kTRUE;
        Bool_t                fProcessRecoPHOS = kTRUE;
        Bool_t                fProcessRecoCalc = kTRUE;
        Bool_t                fSavePartCandPHOS = kTRUE;
        Bool_t                fSavePartCandCalc = kTRUE;
        Bool_t                fFillPHOSPairTreeSE = kTRUE;
        Bool_t                fFillPHOSPairTreeME = kTRUE;
        Bool_t                fSavePHOSMixedBackground = kTRUE;
        Bool_t                fSaveAdditionalBranches = kFALSE;
        Bool_t                fSaveMCBranches = kTRUE;
        Bool_t                fSaveAddMCBranches = kTRUE;

        //Event Cuts
        Double_t              fMaxVertexZ = 10;    

        //Event Mixing Settings
        Int_t                 fEvPoolSize = 10;        //Number of Events used for Event Mixing for Invariant Mass
        Int_t                 fEvTrackSize = 10000;    //Targeted Track Number
        Int_t                 fCentralityBins = 29;
        Double_t              fMinCentBin = -5;
        Double_t              fMaxCentBin = 140;
        Int_t                 fZvtxBins = 10;
        Double_t              fMinZBin = 0;
        Double_t              fMaxZBin = 10;
    	
        //Event Mixing Settings
        Int_t                 fEvPoolSize2 = 20;        //Number of Events used for Event Mixing for Correlation
        Int_t                 fEvTrackSize2 = 10000;    //Targeted Track Number
        Int_t                 fCentralityBins2 = 29;
        Double_t              fMinCentBin2 = -5;
        Double_t              fMaxCentBin2 = 140;
        Int_t                 fZvtxBins2 = 20;
        Double_t              fMinZBin2 = -10;
        Double_t              fMaxZBin2 = 10;

        //Event Mixing Settings
        Bool_t                fRequireSigma = kFALSE;
        Bool_t                fRequireProton = kFALSE;
        Bool_t                fRequireSigmaCand = kTRUE;
        Bool_t                fUseAbsZ = kTRUE;
        Bool_t                fUseAbsZCorr = kTRUE;

        //FillProtonArray Cuts
        Bool_t                fRejectNegIDs = kTRUE;
        Bool_t                fRejectZeroFilterBit = kFALSE;
        Double_t              fMaxProtEta = 1;                   // 0.9  
        Double_t              fMinTPCClustProt = 40;             // 60
        Double_t              fMaxNsigProtTPC = 4;               // 3
        Bool_t                fRequireProtonTPC = kTRUE;  
        Bool_t                fRequireProtonTOF = kFALSE;  
        Bool_t                fRequireProtonTOFforPairs = kTRUE;  
        Double_t              fMaxNsigProtTOF = 6;                // 5
        Double_t              fMaxpOnlyTPCPID = 0.9;              // - 
        Double_t              fMinProtpt = 0;                     // -
        Double_t              fMaxProtpt = 15;                    // - 

        //FillProtonArray2 Cuts
        Double_t              fStrictMaxProtEta = 0.9;               
        Double_t              fStrictMinTPCClustProt = 60;        
        Double_t              fStrictMaxNsigProtTPC = 3;         
        Double_t              fStrictMaxNsigProtTOF = 5;         
        Double_t              fStrictMaxpOnlyTPCPID = 0.8;        
        Double_t              fStrictMinProtpt = 0;              
        Double_t              fStrictMaxProtpt = 5;              

        //ProcessMCParticles Cuts
        Double_t              fMaxMCEta = 0.9;                    // 0.9

        //FillV0PhotonArray Cuts
        Double_t              fMaxDaughtEta = 1.2;                // 0.9    
        Double_t              fMinTPCClustDaught = 20;            // 30
        Double_t              fMaxNsigDaughtTPC = 6;              // 5
        Double_t              fMaxalpha = 1.1;                    // 0.9
        Double_t              fMaxqt = 0.04;                      // 0.03
        Double_t              fMaxopenangle = 0.4;                // 0.3
        Double_t              fMaxdeltatheta = 0.15;              // 0.1
        Double_t              fMinV0CPA = 0.8;                    // 0.8
        Double_t              fMinV0Radius = 1;                   // 3
        Double_t              fMaxV0Radius = 250;                 // 220
        Double_t              fMaxphotonmass = 0.08;               // 0.06

        //FillCaloClusterArray Cuts
        Double_t              fRequirePHOS = kFALSE;
        Double_t              fMinClusterBeta = 0;
        Double_t              fMinClusterDy = 8;    
        Double_t              fMaxClusterM02 = 4;

        //Fill Sigma Tree Cuts     (Coarse Cuts to reduce Tree Size)
        Bool_t                fCleanAutoCorr = kTRUE;
        Double_t              fMinPi0Mass = 0.06;                   // 0.1
        Double_t              fMaxPi0Mass = 0.19;                   // 0.16
        Double_t              fMaxSigmaPA = 0.05;                   // 0.04
        Double_t              fMaxSigmaY  = 0.9;                    // 0.8
        Double_t              fMaxSigmaMass = 1.4;    
        Double_t              fMinProtonDCAxy = 0.005;              // 0.01
        Double_t              fMinProtonDCAz = -1; //Dummy Value
        Double_t              fMaxProtonDCAxy = 5;                  // 2
        Double_t              fMaxProtonDCAz = 9999; //Dummy Value
        Bool_t                fRequireDCACut = kFALSE;

        //Fill Sigma PHOS Tree Cuts     (Coarse Cuts to reduce Tree Size)
        Double_t              fMinPi0MassPHOS = 0.09;                   // 0.1
        Double_t              fMaxPi0MassPHOS = 0.16;                   // 0.15
        Double_t              fMaxSigmaPAPHOS = 0.02;                   // 0.02
        Double_t              fMaxSigmaPAPHOSHM = 0.02;                 // 0.012
        Double_t              fMinSigmaAntiPAPHOS = 0.005;               // 0.012
        Double_t              fMaxProtPhotDCA = 2;                     // 0.05
        Double_t              fMinSigmaDCAtoPVPHOS = 0.5;               // 0.01
        Double_t              fMaxSigmaDCAtoPVPHOS = 60;               // 0.01
        Double_t              fMaxSigmaYPHOS  = 0.9;                    // 0.8
        Double_t              fMaxSigmaMassPHOS = 1.35;    
        Double_t              fMinProtonDCAxyPHOS = 0.005;              // 0.01
        Double_t              fMinProtonDCAzPHOS = -1; //Dummy Value
        Double_t              fMaxProtonDCAxyPHOS = 5;                  // 2
        Double_t              fMaxProtonDCAzPHOS = 9999; //Dummy Value
        Bool_t                fRequireDCACutPHOS = kFALSE;

        //Fill Pair Tree Cuts     (Coarse Cuts to reduce Tree Size)
        Double_t              fMinCorrPi0Mass = 0.1;    
        Double_t              fMaxCorrPi0Mass = 0.16;    
        Double_t              fMaxCorrSigmaPA = 0.06;     
        Double_t              fMinCorrSigmaMass = 1.13;    
        Double_t              fMaxCorrSigmaMass = 1.25;    
        Double_t              fMinCorrProtonDCAxy = 0.005;    
        Double_t              fMaxCorrPairProtonDCAxy = 0.2;    // Pair only with
        Double_t              fMaxCorrPairProtonDCAz = 0.4;     // primary Protons
        Double_t              fMaxCorrkstar = 600;   //Cut in MeV/c!!!

        //Fill Pair Tree Cuts PHOS    (Coarse Cuts to reduce Tree Size)
        Double_t              fMinCorrPi0MassPHOS = 0.09;    
        Double_t              fMaxCorrPi0MassPHOS = 0.16;    
        Double_t              fMaxCorrSigmaPAPHOS = 0.03;     

        //Calculation parameters and cuts
        Bool_t                fCheckProtonV0IDs = kTRUE;
        Float_t               fMaxSigmaMassCalc = 1.35;
        Double_t              fMinSigmaAntiPACalc = 0.005;
        Double_t              fMaxSigmaAntiPACalc = 0.6;    
        Int_t                 fMaxIter=10;
        Bool_t                fRejectNegPar1 = kTRUE;
        Double_t              fMaxSigmaYCalc = 0.9;
        Double_t              fThetaPar0 = 0.000133299;
        Double_t              fThetaPar1 = 0.0016761;
        Double_t              fThetaRange = 0.02;
        Double_t              fPhiPar0 = 0.000129845;
        Double_t              fPhiPar1 = 0.00199688;
        Double_t              fPhiRange = 0.02;

        //Fill Proton Tree Cuts    (Coarse Cuts to reduce Tree Size)
        Float_t              fMaxNsigPionTPC = 2;
        Float_t              fMaxNsigKaonTPC = 2;
        Float_t              fMaxNsigElecTPC = 2;
        Float_t              fMaxNsigPionTOF = 2;
        Float_t              fMaxNsigKaonTOF = 2;
        Float_t              fMaxNsigElecTOF = 2;
        Bool_t               fRejOtherTOF    = kTRUE;

    public: //Setter-Functions

        //Remove generated pile-up
        void SetRemovePileup(Bool_t removepileup) {fRemoveGenPileup = removepileup;}

        //Activate some extra output. Caution: Very anoying!
        void SetActivateDebug(Bool_t debug) {fDebug = debug;}

        //Activate and deactivate parts of the Analysis 
        void SetActivateMCParticles(Bool_t processmcpart) {fProcessMCParticles = processmcpart;}
        void SetActivateProton(Bool_t processprotons) {fProcessProtons = processprotons;}
        void SetActivateSaveProton(Bool_t saveprotons) {fSaveProtons = saveprotons;}
        void SetActivateV0s(Bool_t processv0s) {fProcessV0s = processv0s;}
        void SetActivateReco(Bool_t processreco) {fProcessReco = processreco;}
        void SetActivateCandidates(Bool_t savecand) {fSavePartCand = savecand;}
        void SetActivateSETree(Bool_t saveSE) {fFillPairTreeSE = saveSE;}
        void SetActivateMETree(Bool_t saveME) {fFillPairTreeME = saveME;}
        void SetActivateMixedBackground(Bool_t savemebkg) {fSaveMixedBackground = savemebkg;}
        void SetActivateClusters(Bool_t processclusters) {fProcessClusters = processclusters;}
        void SetActivateClusterMapping(Bool_t mapclusters) {fMapClusterstoTracks = mapclusters;}
        void SetActivateRecoPhos(Bool_t recophos) {fProcessRecoPHOS = recophos;}
        void SetActivateRecoCalc(Bool_t recocalc) {fProcessRecoCalc = recocalc;}
        void SetActivatePhosCand(Bool_t savephos) {fSavePartCandPHOS = savephos;}
        void SetActivateCalcCand(Bool_t savecalc) {fSavePartCandCalc = savecalc;}
        void SetActivatePHOSSETree(Bool_t savePhosSE) {fFillPHOSPairTreeSE = savePhosSE;}
        void SetActivatePHOSMETree(Bool_t savePhosME) {fFillPHOSPairTreeME = savePhosME;}
        void SetActivatePHOSMEBkg(Bool_t savephosmebkg) {fSavePHOSMixedBackground = savephosmebkg;}
        void SetActivateExtraBranches(Bool_t savextrabr) {fSaveAdditionalBranches = savextrabr;}
        void SetActivateMCBranches(Bool_t savemcbr) {fSaveMCBranches = savemcbr;}
        void SetActivateAddMCBranches(Bool_t saveaddmcbr) {fSaveAddMCBranches = saveaddmcbr;}

        //Event Cuts
        void SetMaxVertexZ(Double_t maxvertexz) {fMaxVertexZ = maxvertexz;}
        
        //Event Mixing Settings
        void SetEventPoolSize(Int_t poolsize) {fEvPoolSize = poolsize;}
        void SetTrackBufferSize(Int_t tracksize) {fEvTrackSize = tracksize;}
        void SetNMultBins(Int_t nmult) {fCentralityBins = nmult;}
        void SetMinMultBin(Double_t minmult) {fMinCentBin = minmult;}
        void SetMaxMultBin(Double_t maxmult) {fMaxCentBin = maxmult;}
        void SetNZVertBins(Int_t nzvert) {fZvtxBins = nzvert;}
        void SetMinZVertBin(Double_t minz) {fMinZBin = minz;}
        void SetMaxZVertBin(Double_t maxz) {fMaxZBin = maxz;}

        //Event Mixing Settings
        void SetEventPoolSize2(Int_t poolsize) {fEvPoolSize2 = poolsize;}
        void SetTrackBufferSize2(Int_t tracksize) {fEvTrackSize2 = tracksize;}
        void SetNMultBins2(Int_t nmult) {fCentralityBins2 = nmult;}
        void SetMinMultBin2(Double_t minmult) {fMinCentBin2 = minmult;}
        void SetMaxMultBin2(Double_t maxmult) {fMaxCentBin2 = maxmult;}
        void SetNZVertBins2(Int_t nzvert) {fZvtxBins2 = nzvert;}
        void SetMinZVertBin2(Double_t minz) {fMinZBin2 = minz;}
        void SetMaxZVertBin2(Double_t maxz) {fMaxZBin2 = maxz;}

        //Event Mixing Settings
        void SetRequireSigma(Bool_t reqsig) {fRequireSigma = reqsig;}
        void SetRequireProton(Bool_t reqprt) {fRequireProton = reqprt;}
        void SetRequireSigmaCand(Bool_t reqsigc) {fRequireSigmaCand = reqsigc;}
        void SetUseAbsZ(Bool_t useabsz) {fUseAbsZ = useabsz;}
        void SetUseAbsZCorr(Bool_t useabszcorr) {fUseAbsZCorr = useabszcorr;}

        //FillProtonArray Cuts
        void SetRejectNegIDs(Bool_t rejnegid) {fRejectNegIDs = rejnegid;}
        void SetRejectZeroFilterBit(Bool_t rejzerofbit) {fRejectZeroFilterBit = rejzerofbit;}
        void SetProtonMaxEta(Double_t maxeta) {fMaxProtEta = maxeta;}
        void SetProtonMinTPCCluster(Double_t minclst) {fMinTPCClustProt = minclst;}
        void SetProtonMaxNSigmaTPC(Double_t nsigtpc) {fMaxNsigProtTPC = nsigtpc;}
        void SetProtonRequireTPC(Bool_t requiretpc) {fRequireProtonTPC = requiretpc;}
        void SetProtonRequireTOF(Bool_t requiretof) {fRequireProtonTOF = requiretof;}
        void SetProtonRequireTOFforPairs(Bool_t requirepairtof) {fRequireProtonTOFforPairs = requirepairtof;}
        void SetProtonMaxNSigmaTOF(Double_t nsigtof) {fMaxNsigProtTOF = nsigtof;}
        void SetProtonMaxpOnlyTPCPID(Double_t pmaxonlytpc) {fMaxpOnlyTPCPID = pmaxonlytpc;}
        void SetProtonMinpt(Double_t minpt) {fMinProtpt = minpt;}
        void SetProtonMaxpt(Double_t maxpt) {fMaxProtpt = maxpt;}

        //FillProtonArray2 Cuts
        void SetStrictProtonMaxEta(Double_t maxeta) {fStrictMaxProtEta = maxeta;}
        void SetStrictProtonMinTPCCluster(Double_t minclst) {fStrictMinTPCClustProt = minclst;}
        void SetStrictProtonMaxNSigmaTPC(Double_t nsigtpc) {fStrictMaxNsigProtTPC = nsigtpc;}
        void SetStrictProtonMaxNSigmaTOF(Double_t nsigtof) {fStrictMaxNsigProtTOF = nsigtof;}
        void SetStrictProtonMaxpOnlyTPCPID(Double_t pmaxonlytpc) {fStrictMaxpOnlyTPCPID = pmaxonlytpc;}
        void SetStrictProtonMinpt(Double_t minpt) {fStrictMinProtpt = minpt;}
        void SetStrictProtonMaxpt(Double_t maxpt) {fStrictMaxProtpt = maxpt;}

        //ProcessMCParticles Cuts
        void SetMCMaxEta(Double_t maxeta) {fMaxMCEta = maxeta;}

        //FillV0PhotonArray Cuts
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

        //FillCaloClusterArray Cuts
        void SetRequirePHOS(Bool_t requirephos) {fRequirePHOS = requirephos;}
        void SetMinClusterBeta(Double_t minclustbeta) {fMinClusterBeta = minclustbeta;}
        void SetMinClusterDy(Double_t minclustdy) {fMinClusterDy = minclustdy;}
        void SetMaxClusterM02(Double_t maxclustm02) {fMaxClusterM02 = maxclustm02;}

        //Fill Sigma Tree Cuts     (Coarse Cuts to reduce Tree Size)
        void SetCleanAutoCorrelations(Bool_t cleanpairs) {fCleanAutoCorr = cleanpairs;}        
        void SetPi0MinMass(Double_t minpi0mass) {fMinPi0Mass = minpi0mass;}
        void SetPi0MaxMass(Double_t maxpi0mass) {fMaxPi0Mass = maxpi0mass;}
        void SetSigmaMaxPA(Double_t maxpa) {fMaxSigmaPA = maxpa;}
        void SetSigmaMaxY(Double_t maxy) {fMaxSigmaY = maxy;}
        void SetSigmaMaxMass(Double_t maxsigmass) {fMaxSigmaMass = maxsigmass;}
        void SetProtonMinDCAxy(Double_t minprotondcaxy) {fMinProtonDCAxy = minprotondcaxy;}
        void SetProtonMinDCAz(Double_t minprotondcaz) {fMinProtonDCAz = minprotondcaz;}
        void SetProtonMaxDCAxy(Double_t maxprotondcaxy) {fMaxProtonDCAxy = maxprotondcaxy;}
        void SetProtonMaxDCAz(Double_t maxprotondcaz) {fMaxProtonDCAz = maxprotondcaz;}
        void SetRequireDCA(Bool_t requiredca) {fRequireDCACut = requiredca;}

        //Fill Sigma PHOS Tree Cuts     (Coarse Cuts to reduce Tree Size)
        void SetPi0MinMassPHOS(Double_t minpi0mass) {fMinPi0MassPHOS = minpi0mass;}
        void SetPi0MaxMassPHOS(Double_t maxpi0mass) {fMaxPi0MassPHOS = maxpi0mass;}
        void SetSigmaMaxPAPHOS(Double_t maxpa) {fMaxSigmaPAPHOS = maxpa;}
        void SetSigmaMaxPAPHOSHM(Double_t maxhmpa) {fMaxSigmaPAPHOSHM = maxhmpa;}
        void SetSigmaMinAntiPAPHOS(Double_t minantipa) {fMinSigmaAntiPAPHOS = minantipa;}
        void SetMaxProtPhotDCA(Double_t maxppdca) {fMaxProtPhotDCA = maxppdca;}
        void SetSigmaMinDCAtoPVPHOS(Double_t mindcatopv) {fMinSigmaDCAtoPVPHOS = mindcatopv;}
        void SetSigmaMaxDCAtoPVPHOS(Double_t maxdcatopv) {fMaxSigmaDCAtoPVPHOS = maxdcatopv;}
        void SetSigmaMaxYPHOS(Double_t maxy) {fMaxSigmaYPHOS = maxy;}
        void SetSigmaMaxMassPHOS(Double_t maxsigmass) {fMaxSigmaMassPHOS = maxsigmass;}
        void SetProtonMinDCAxyPHOS(Double_t minprotondcaxy) {fMinProtonDCAxyPHOS = minprotondcaxy;}
        void SetProtonMinDCAzPHOS(Double_t minprotondcaz) {fMinProtonDCAzPHOS = minprotondcaz;}
        void SetProtonMaxDCAxyPHOS(Double_t maxprotondcaxy) {fMaxProtonDCAxyPHOS = maxprotondcaxy;}
        void SetProtonMaxDCAzPHOS(Double_t maxprotondcaz) {fMaxProtonDCAzPHOS = maxprotondcaz;}
        void SetRequireDCAPHOS(Bool_t requiredca) {fRequireDCACutPHOS = requiredca;}

        //Fill Pair Tree Cuts     (Coarse Cuts to reduce Tree Size)
        void SetCorrPi0MinMass(Double_t minpi0mass) {fMinCorrPi0Mass = minpi0mass;}
        void SetCorrPi0MaxMass(Double_t maxpi0mass) {fMaxCorrPi0Mass = maxpi0mass;}
        void SetCorrSigmaMaxPA(Double_t maxpa) {fMaxCorrSigmaPA = maxpa;}
        void SetCorrSigmaMinMass(Double_t minsigmass) {fMinCorrSigmaMass = minsigmass;}
        void SetCorrSigmaMaxMass(Double_t maxsigmass) {fMaxCorrSigmaMass = maxsigmass;}
        void SetCorrProtonMinDCAxy(Double_t minprotondcaxy) {fMinCorrProtonDCAxy = minprotondcaxy;}
        void SetCorrPairProtonMaxDCAxy(Double_t maxpairprotondcaxy) {fMaxCorrPairProtonDCAxy = maxpairprotondcaxy;}
        void SetCorrPairProtonMaxDCAz(Double_t maxpairprotondcaz) {fMaxCorrPairProtonDCAz = maxpairprotondcaz;}
        void SetPairMaxkstar(Double_t maxkstar) {fMaxCorrkstar = maxkstar;}

        //Fill Pair Tree Cuts PHOS    (Coarse Cuts to reduce Tree Size)
        void SetCorrPi0MinMassPHOS(Double_t minpi0mass) {fMinCorrPi0MassPHOS = minpi0mass;}
        void SetCorrPi0MaxMassPHOS(Double_t maxpi0mass) {fMaxCorrPi0MassPHOS = maxpi0mass;}
        void SetCorrSigmaMaxPAPHOS(Double_t maxpa) {fMaxCorrSigmaPAPHOS = maxpa;}

        //Calculation parameters and cuts
        void SetCheckProtonV0IDs(Bool_t checkid) {fCheckProtonV0IDs = checkid;}
        void SetMaxSigmaMassCalc(Float_t maxmass) {fMaxSigmaMassCalc = maxmass;}
        void SetSigmaMinAntiPACalc(Double_t minantipa) {fMinSigmaAntiPACalc = minantipa;}
        void SetSigmaMaxAntiPACalc(Double_t maxantipa) {fMaxSigmaAntiPACalc = maxantipa;}
        void SetMaxIter(Int_t maxiter) {fMaxIter = maxiter;}
        void SetRejectNegPar1(Bool_t rejnegpar1) { fRejectNegPar1 = rejnegpar1;}
        void SetSigmaMaxYCalc(Double_t maxy) {fMaxSigmaYCalc = maxy;}
        void SetThetaPar0(Double_t thetapar0) {fThetaPar0 = thetapar0;}
        void SetThetaPar1(Double_t thetapar1) {fThetaPar1 = thetapar1;}
        void SetThetaRange(Double_t thetarange) {fThetaRange = thetarange;}
        void SetPhiPar0(Double_t phipar0) {fPhiPar0 = phipar0;}
        void SetPhiPar1(Double_t phipar1) {fPhiPar1 = phipar1;}
        void SetPhiRange(Double_t phirange) {fPhiRange = phirange;}

        //Fill Proton Tree Cuts    (Coarse Cuts to reduce Tree Size)
        void SetPionMaxNSigmaTPC(Float_t nsigtpcpion) {fMaxNsigPionTPC = nsigtpcpion;}
        void SetKaonMaxNSigmaTPC(Float_t nsigtpckaon) {fMaxNsigKaonTPC = nsigtpckaon;}
        void SetElecMaxNSigmaTPC(Float_t nsigtpcelec) {fMaxNsigElecTPC = nsigtpcelec;}
        void SetPionMaxNSigmaTOF(Float_t nsigtofpion) {fMaxNsigPionTOF = nsigtofpion;}
        void SetKaonMaxNSigmaTOF(Float_t nsigtofkaon) {fMaxNsigKaonTOF = nsigtofkaon;}
        void SetElecMaxNSigmaTOF(Float_t nsigtofelec) {fMaxNsigElecTOF = nsigtofelec;}
        void SetRejectOtherTOF(Bool_t rejothertof) {fRejOtherTOF = rejothertof;}

    private:

        //Branches of TTree "fSigmaCandTree"
        Bool_t                fIsMCSigma;
        Bool_t                fIsMCDalitz;
        Bool_t                fIsMCPrimary;
        Bool_t                fIsV01Onthefly;
        Bool_t                fIsV02Onthefly;
        Bool_t                fHas4DiffIDs;
        Int_t                 fSigRunnumber;
        UInt_t                fSigTriggerMask;               
        Int_t                 fSigMCLabel;
        Int_t                 fSigProtonID;
        ULong64_t             fSigProtonStatus;
        UInt_t                fSigProtonFilterMap;
        ULong64_t             fSigEventID;
        Float_t               fSigCentrality;               
        Short_t               fSigRefMultComb05;        
        Short_t               fSigRefMultComb08;        
        Short_t               fSigRefMultComb10;        
        Float_t               fSigBField;               
        Float_t               fInvSigMass;               
        Float_t               fInvSigpropMass;
        Float_t               fInvSigMassUncorr;
        Float_t               fSigY;               
        Float_t               fSigYprop;               
        Float_t               fSigPA;               
        Float_t               fSigPAprop;               
        Short_t               fSigCharge;        
        Float_t               fSigAntiPA;               
        Float_t               fSigPx;        
        Float_t               fSigPy;        
        Float_t               fSigPz;        
        Float_t               fSigPt;        
        Float_t               fSigPxprop;        
        Float_t               fSigPyprop;        
        Float_t               fSigPzprop;        
        Float_t               fPrimVertX;        
        Float_t               fPrimVertY;        
        Float_t               fPrimVertZ;        
        Float_t               fPrimVertXMC;        
        Float_t               fPrimVertYMC;        
        Float_t               fPrimVertZMC;        
        Float_t               fSigDecayVertX;        
        Float_t               fSigDecayVertY;        
        Float_t               fSigDecayVertZ;        
        Float_t               fSigFlightDist;        
        Float_t               fSigDecayVertXMC;        
        Float_t               fSigDecayVertYMC;        
        Float_t               fSigDecayVertZMC;        
        Float_t               fSigPxMC;        
        Float_t               fSigPyMC;        
        Float_t               fSigPzMC;        
        Float_t               fPhoton1Px;        
        Float_t               fPhoton1Py;        
        Float_t               fPhoton1Pz;        
        Float_t               fPhoton2Px;        
        Float_t               fPhoton2Py;        
        Float_t               fPhoton2Pz;        
        Float_t               fPhotonDaughtMaxEta;        
        Float_t               fPhotonsMaxDeltaTheta;        
        Float_t               fPhoton1CPA;        
        Float_t               fPhoton2CPA;        
        Float_t               fPhoton1Radius;        
        Float_t               fPhoton2Radius;        
        Float_t               fPhoton1DCAPV;        
        Float_t               fPhoton2DCAPV;
        Float_t               fPhoton1DCASV;        
        Float_t               fPhoton2DCASV;
        Float_t               fTrackDCASV;  
        Float_t               fTrackDCASVKF;  
        Float_t               fKFChi2;  
        Short_t               fPhotonsMinCluster;
        Short_t               fPhotonsMinITSCluster;
        Float_t               fPhotonsMaxalpha;
        Float_t               fPhotonsMaxqt;
        Float_t               fPhotonsMaxOpenAngle;
        Float_t               fPhotonsMaxinvmass;        
        Float_t               fPhotonsMaxNSigTPC;
        Float_t               fPhotonsMaxChi2;
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
        Float_t               fProtonX;        
        Float_t               fProtonY;        
        Float_t               fProtonZ;        
        Float_t               fProtonEta;        
        Float_t               fProtonpropPx;        
        Float_t               fProtonpropPy;        
        Float_t               fProtonpropPz;        
        Float_t               fProtonDCAtoPVxy;        
        Float_t               fProtonDCAtoPVz;        
        Float_t               fProtonPi0DCA;
        Float_t               fProtonNSigTPC;
        Float_t               fProtonNSigTOF;
        Short_t               fProtonNCluster;  
        Short_t               fProtonNITSCluster;  
        Float_t               fProtonChi2;  
        Float_t               fProtonNSigTPCPion;
        Float_t               fProtonNSigTPCKaon;
        Float_t               fProtonNSigTPCElec;
        Float_t               fProtonNSigTOFPion;
        Float_t               fProtonNSigTOFKaon;
        Float_t               fProtonNSigTOFElec;
        Int_t                 fTrackLabel;
        Int_t                 fTrackPDGCode;
        Int_t                 fTrackMotherID;       // -1 if Particle is primary
        Int_t                 fTrackMotherPDGCode;
        Float_t               fTrackMotherMCPx;
        Float_t               fTrackMotherMCPy;
        Float_t               fTrackMotherMCPz;
        Float_t               fTrackMCPx;
        Float_t               fTrackMCPy;
        Float_t               fTrackMCPz;
        Int_t                 fPhoton1Label;       // -1 if fake V0
        Int_t                 fPhoton1PDGCode;
        Int_t                 fPhoton1MotherID;
        Int_t                 fPhoton1MotherPDGCode;
        Int_t                 fPhoton1GMotherID;
        Int_t                 fPhoton1GMotherPDGCode;
        Float_t               fPhoton1MCPx;
        Float_t               fPhoton1MCPy;
        Float_t               fPhoton1MCPz;
        Int_t                 fPhoton2Label;       // -1 if fake V0
        Int_t                 fPhoton2PDGCode;
        Int_t                 fPhoton2MotherID;
        Int_t                 fPhoton2MotherPDGCode;
        Int_t                 fPhoton2GMotherID;
        Int_t                 fPhoton2GMotherPDGCode;
        Float_t               fPhoton2MCPx;
        Float_t               fPhoton2MCPy;
        Float_t               fPhoton2MCPz;
        //End of Branches of "fSigmaCandTree"

        //Extra Branches of TTree "fSigmaPHOSCandTree"
        Bool_t                fIsClusterEMCAL;
        Bool_t                fIsClusterPHOS; 
        Float_t               fCaloPhotonX;
        Float_t               fCaloPhotonY;
        Float_t               fCaloPhotonZ;
        Float_t               fConvPhotonX;
        Float_t               fConvPhotonY;
        Float_t               fConvPhotonZ;
        Float_t               fConvPhotonSecPA;
        Float_t               fCaloPhotonPxMC;
        Float_t               fCaloPhotonPyMC;
        Float_t               fCaloPhotonPzMC;
        Float_t               fCaloPhotonE;
        Float_t               fCaloPhotonELead;
        Float_t               fCaloPhotonEcorr;
        Float_t               fCaloPhotonEMC;
        Int_t                 fClustNLabels;
        Int_t                 fClustPDGCode;
        Float_t               fClustDispersion;
        Float_t               fClustM20;
        Float_t               fClustM02;
        Int_t                 fClustNTracksMatched;
        Float_t               fClustTrackDx;
        Float_t               fClustTrackDz;
        Float_t               fClustTrackD;
        Float_t               fClustTOF;
        Float_t               fClustBeta;
        Int_t                 fClustNCells;
        Float_t               fClustDisttoBC;
        Short_t               fCellsAbsId[20];
        //End of Extra Branches of "fSigmaPHOSCandTree"

        //Extra Branches of TTree "fSigmaPairTree"
        Bool_t                fPairProtonIsMC;        
        Int_t                 fPairProtonPDG;        
        Bool_t                fPairProtonIsPrimary;        
        Float_t               fPairProtonPx;        
        Float_t               fPairProtonPy;        
        Float_t               fPairProtonPz;        
        Float_t               fPairProtonPxMC;        
        Float_t               fPairProtonPyMC;        
        Float_t               fPairProtonPzMC;        
        Float_t               fPairProtonP;        
        Float_t               fPairProtonEta;        
        Short_t               fPairProtonCharge;        
        Float_t               fPairProtonDCAtoPVxy;
        Float_t               fPairProtonDCAtoPVz;        
        Float_t               fPairProtonNSigTPC;        
        Float_t               fPairProtonNSigTOF;
        Float_t               fPairProtNSigTPCPion;        
        Float_t               fPairProtNSigTPCKaon;
        Float_t               fPairProtNSigTPCElec;
        Float_t               fPairProtNSigTOFPion;        
        Float_t               fPairProtNSigTOFKaon;
        Float_t               fPairProtNSigTOFElec;
        Float_t               fPairProtonChi2;
        Short_t               fPairProtonCluster;
        Short_t               fPairProtonITSCluster;
        Int_t                 fPairProtonID;
        ULong64_t             fPairProtonStatus;
        UInt_t                fPairProtonFilterMap;
        //End of Extra Branches of "fSigmaPairTree"

        //End of Extra Branches of "fSigmaCalcCandTree"        
        Short_t               fNIter;
        Float_t               fDeltaPhi;
        Float_t               fDeltaTheta;
        //End of Extra Branches of "fSigmaCalcCandTree"

        //Extra Branches of TTree "fSigmaPairTreeSE/ME"
        Float_t               fSigmaProtonkstar;  //Saved in MeV/c
        Float_t               fSigmaProtonpropkstar; 
        //End of Extra Branches of "fSigmaPairTreeSE/ME"

        AliAnalysisTaskSigmaPlus(const AliAnalysisTaskSigmaPlus&); // not implemented
        AliAnalysisTaskSigmaPlus& operator=(const AliAnalysisTaskSigmaPlus&); // not implemented

        ClassDef(AliAnalysisTaskSigmaPlus, 1);
};

class AliAODTrackreduced : public TObject //Version with Cov Matrix for Reconstruction
{
    //This class stores only the info of the AliAODTrack needed for  
    //the Event Mixing. Inherits from TObject so it can be stored in
    //a TObjArray needed for AliEventPool. Uses (mostly) same getters
    //as in AliAODTrack for simplicity
    public:
        //AliAODTrackreduced(){}
        AliAODTrackreduced() :
        charge(-999), tpcncls(-999), id(-999), status(0), filtermap(0), tpcchi2(-999), dcaxy(-999), dcaz(-999),
        nsigmatpcproton(-999), nsigmatpcpion(-999), nsigmatpckaon(-999), nsigmatpcelectron(-999), 
        nsigmatofproton(-999), nsigmatofpion(-999), nsigmatofkaon(-999), nsigmatofelectron(-999)
        {
            for(Int_t i=0; i<3; i++){x[i]=-999; p[i]=-999;}            
            for(Int_t i=0; i<21; i++){covMatrix[i]=-999;}            
        }
        virtual ~AliAODTrackreduced() {}

        //Setter
        void InitfromTrack(const AliAODTrack *aodTrack, const AliPIDResponse* pidresp=0, Double_t primaryVertexX = 0, Double_t primaryVertexY = 0,  Double_t primaryVertexZ = 0) { 
            if(!aodTrack) {std::cout << "WARNING: Input source 'AOD Track' does not exit!\n"; return;}            
            aodTrack->GetXYZ(x);
            aodTrack->GetPxPyPz(p);
            aodTrack->GetCovarianceXYZPxPyPz(covMatrix);
            aodTrack->GetImpactParameters(dcaxy,dcaz);
            //Substract PV and move Track to PV of the Pion      
            x[0] -= primaryVertexX;
            x[1] -= primaryVertexY;
            x[2] -= primaryVertexZ;
            charge = aodTrack->Charge();
            tpcncls = aodTrack->GetTPCNcls();
            itsncls = aodTrack->GetITSNcls();
            id = aodTrack->GetID();
            status = aodTrack->GetStatus();
            filtermap = aodTrack->GetFilterMap();
            tpcchi2 = aodTrack->GetTPCchi2();
            if(!pidresp) {std::cout << "WARNING: Input source 'PID Response' does not exit!\n"; return;}            
            nsigmatpcproton = pidresp->NumberOfSigmasTPC(aodTrack,AliPID::kProton);
            nsigmatpcpion = pidresp->NumberOfSigmasTPC(aodTrack,AliPID::kPion);
            nsigmatpckaon = pidresp->NumberOfSigmasTPC(aodTrack,AliPID::kKaon);
            nsigmatpcelectron = pidresp->NumberOfSigmasTPC(aodTrack,AliPID::kElectron);
            nsigmatofproton = pidresp->NumberOfSigmasTOF(aodTrack,AliPID::kProton);
            nsigmatofpion = pidresp->NumberOfSigmasTOF(aodTrack,AliPID::kPion);
            nsigmatofkaon = pidresp->NumberOfSigmasTOF(aodTrack,AliPID::kKaon);
            nsigmatofelectron = pidresp->NumberOfSigmasTOF(aodTrack,AliPID::kElectron);
            return;
        }

        //Getters
        Double_t Px() const { return p[0]; } 
        Double_t Py() const { return p[1]; } 
        Double_t Pz() const { return p[2]; } 
        Short_t Charge() const { return charge; } 
        Int_t GetTPCNcls() const { return tpcncls; } 
        Int_t GetITSNcls() const { return itsncls; } 
        Int_t GetID() const { return id; } 
        ULong64_t GetStatus() const { return status; } 
        UInt_t GetFilterMap() const { return filtermap; } 
        Double_t GetTPCchi2() const { return tpcchi2; } 

        void GetXYZ(Double_t xx[3]) const { for(Int_t i=0; i<3; i++){xx[i]=x[i];} return; }
        void GetPxPyPz(Double_t pp[3]) const { for(Int_t i=0; i<3; i++){pp[i]=p[i];} return; }
        void GetCovarianceXYZPxPyPz(Double_t cv[21]) const { for(Int_t i=0; i<21; i++){cv[i]=covMatrix[i];} return; }
        void GetImpactParameters(Float_t &xy, Float_t &z) const { xy = dcaxy; z = dcaz; return; }

        Float_t NumberOfSigmasTPCProton() const { return nsigmatpcproton; }
        Float_t NumberOfSigmasTPCPion() const { return nsigmatpcpion; }
        Float_t NumberOfSigmasTPCKaon() const { return nsigmatpckaon; }
        Float_t NumberOfSigmasTPCElectron() const { return nsigmatpcelectron; }
        Float_t NumberOfSigmasTOFProton() const { return nsigmatofproton; }
        Float_t NumberOfSigmasTOFPion() const { return nsigmatofpion; }
        Float_t NumberOfSigmasTOFKaon() const { return nsigmatofkaon; }
        Float_t NumberOfSigmasTOFElectron() const { return nsigmatofelectron; }        

    protected:
        Double_t  p[3];
        Double_t  x[3];
        Double_t  covMatrix[21];
        Short_t   charge;
        Int_t     tpcncls;
        Int_t     itsncls;
        Int_t     id;
        ULong64_t status;
        UInt_t    filtermap;
        Double_t  tpcchi2;
        Float_t   dcaxy;
        Float_t   dcaz;
        Float_t   nsigmatpcproton; 
        Float_t   nsigmatpcpion;
        Float_t   nsigmatpckaon; 
        Float_t   nsigmatpcelectron;
        Float_t   nsigmatofproton; 
        Float_t   nsigmatofpion;
        Float_t   nsigmatofkaon; 
        Float_t   nsigmatofelectron;

    ClassDef(AliAODTrackreduced, 1); // class required for event mixing
};

class AliAODTrackcorrelation : public TObject //Version without Cov Matrix for Correlations
{
    //This class stores only the info of the AliAODTrack needed for  
    //the Event Mixing. Inherits from TObject so it can be stored in
    //a TObjArray needed for AliEventPool. Uses (mostly) same getters
    //as in AliAODTrack for simplicity
    public:
        //AliAODTrackcorrelation(){}
        AliAODTrackcorrelation() :
        charge(-999), tpcncls(-999), id(-999), status(0), filtermap(0), tpcchi2(-999), dcaxy(-999), dcaz(-999), 
        nsigmatpcproton(-999), nsigmatpckaon(-999), nsigmatpcpion(-999), nsigmatofproton(-999), nsigmatofkaon(-999), nsigmatofpion(-999)
        {
            for(Int_t i=0; i<3; i++){p[i]=-999;}            
        }
        virtual ~AliAODTrackcorrelation() {}

        //Setter
        void InitfromTrack(const AliAODTrack *aodTrack, const AliPIDResponse* pidresp = 0) { 
            if(!aodTrack) {std::cout << "WARNING: Input source 'AOD Track' does not exit!\n"; return;}            
            aodTrack->GetPxPyPz(p);
            aodTrack->GetImpactParameters(dcaxy,dcaz);
            charge = aodTrack->Charge();
            tpcncls = aodTrack->GetTPCNcls();
            itsncls = aodTrack->GetITSNcls();
            id = aodTrack->GetID();
            status = aodTrack->GetStatus();
            filtermap = aodTrack->GetFilterMap();
            tpcchi2 = aodTrack->GetTPCchi2();
            if(!pidresp) {std::cout << "WARNING: Input source 'PID Response' does not exit!\n"; return;}            
            nsigmatpcproton = pidresp->NumberOfSigmasTPC(aodTrack,AliPID::kProton);
            nsigmatpckaon = pidresp->NumberOfSigmasTPC(aodTrack,AliPID::kKaon);
            nsigmatpcpion = pidresp->NumberOfSigmasTPC(aodTrack,AliPID::kPion);
            nsigmatofproton = pidresp->NumberOfSigmasTOF(aodTrack,AliPID::kProton);
            nsigmatofkaon = pidresp->NumberOfSigmasTOF(aodTrack,AliPID::kKaon);
            nsigmatofpion = pidresp->NumberOfSigmasTOF(aodTrack,AliPID::kPion);
            return;
        }

        //Getters
        Double_t Px() const { return p[0]; } 
        Double_t Py() const { return p[1]; } 
        Double_t Pz() const { return p[2]; } 
        Short_t Charge() const { return charge; } 
        Int_t GetTPCNcls() const { return tpcncls; } 
        Int_t GetITSNcls() const { return itsncls; } 
        Int_t GetID() const { return id; } 
        ULong64_t GetStatus() const { return status; } 
        UInt_t GetFilterMap() const { return filtermap; } 
        Float_t GetTPCchi2() const { return tpcchi2; } 

        void GetPxPyPz(Double_t pp[3]) const { for(Int_t i=0; i<3; i++){pp[i]=p[i];} return; }
        void GetImpactParameters(Float_t &xy, Float_t &z) const { xy = dcaxy; z = dcaz; return; }

        Float_t NumberOfSigmasTPCProton() const { return nsigmatpcproton; }
        Float_t NumberOfSigmasTPCKaon() const { return nsigmatpckaon; }
        Float_t NumberOfSigmasTPCPion() const { return nsigmatpcpion; }
        Float_t NumberOfSigmasTOFProton() const { return nsigmatofproton; }
        Float_t NumberOfSigmasTOFKaon() const { return nsigmatofkaon; }
        Float_t NumberOfSigmasTOFPion() const { return nsigmatofpion; }

    protected:
        Double_t   p[3];
        Short_t   charge;
        Int_t     tpcncls;
        Int_t     itsncls;
        Int_t     id;
        ULong64_t status;
        UInt_t    filtermap;
        Float_t   tpcchi2;
        Float_t   dcaxy;
        Float_t   dcaz;
        Float_t   nsigmatpcproton; 
        Float_t   nsigmatpckaon; 
        Float_t   nsigmatpcpion; 
        Float_t   nsigmatofproton; 
        Float_t   nsigmatofkaon; 
        Float_t   nsigmatofpion; 

    ClassDef(AliAODTrackcorrelation, 1); // class required for event mixing
};

class AliAODClusterreduced : public TObject 
{
    //This class stores only the info of the AliAODCaloCluster needed for  
    //the Event Mixing. Inherits from TObject so it can be stored in
    //a TObjArray needed for AliEventPool. Uses (mostly) same getters
    //as in AliAODCaloCluster for simplicity
    public:
        //AliAODClusterreduced(){}
        AliAODClusterreduced() :
        isEmcal(-999), isPhos(-999), Energy(-999), EnergyLead(-999), Dispersion(-999), M20(-999), M02(-999), 
        NTracksMatched(-999), TrackDx(-999), TrackDz(-999), TOF(-999), NCells(-999), DisttoBC(-999)
        {
            for(Int_t i=0; i<3; i++){x[i]=-999;}            
            for(Int_t i=0; i<20; i++){CellIds[i]=-999;}            
        }
        virtual ~AliAODClusterreduced() {}

        //Setter
        void InitfromCluster(const AliAODCaloCluster *aodCluster, AliAODCaloCells *caloCells) { 
            if(!aodCluster) {std::cout << "WARNING: Input source 'AOD Cluster' does not exit!\n"; return;}            
            isEmcal        = aodCluster->IsEMCAL();
            isPhos         = aodCluster->IsPHOS();
            Energy         = aodCluster->E();
            Dispersion     = aodCluster->GetDispersion();
            M20            = aodCluster->GetM20();
            M02            = aodCluster->GetM02();
            NTracksMatched = aodCluster->GetNTracksMatched();
            TrackDx        = aodCluster->GetTrackDx();
            TrackDz        = aodCluster->GetTrackDz();
            TOF            = aodCluster->GetTOF();
            NCells         = aodCluster->GetNCells();
            DisttoBC       = aodCluster->GetDistanceToBadChannel();
            aodCluster->GetPosition(x);
            for(Int_t i=0; i<20; i++){CellIds[i] = (Short_t)aodCluster->GetCellAbsId(i);}
            EnergyLead = -999;
            if(caloCells){
              for(Int_t iCell=0; iCell<NCells; iCell++){double cellamp = caloCells->GetCellAmplitude(CellIds[iCell]); if(cellamp>EnergyLead) EnergyLead = cellamp;}
            }
            return;
        }

        //Getters
        Bool_t   IsEMCAL() const { return isEmcal; }
        Bool_t   IsPHOS() const { return isPhos; }
        Float_t  E() const { return Energy; }
        Float_t  ELead() const { return EnergyLead; }
        Float_t  GetDispersion() const { return Dispersion; }
        Float_t  GetM20() const { return M20; }
        Float_t  GetM02() const { return M02; }
        Int_t    GetNTracksMatched() const { return NTracksMatched; }
        Float_t  GetTrackDx() const { return TrackDx; }
        Float_t  GetTrackDz() const { return TrackDz; }
        Double_t GetTOF() const { return TOF; }
        Int_t    GetNCells() const { return NCells; }
        Float_t  GetDistanceToBadChannel() const { return DisttoBC; }

        void GetPosition(Float_t xx[3]) const { for(Int_t i=0; i<3; i++){xx[i]=x[i];} return; }

        Short_t GetCellAbsId(Int_t i) const {  
            if(i>=0&&i<20){ return CellIds[i]; }    
            return -1;
        }

    protected:
        Bool_t   isEmcal;
        Bool_t   isPhos;
        Float_t  x[3];
        Float_t  Energy;
        Float_t  EnergyLead;
        Float_t  Dispersion;
        Float_t  M20;
        Float_t  M02;
        Int_t    NTracksMatched;
        Float_t  TrackDx;
        Float_t  TrackDz;
        Double_t TOF;
        Int_t    NCells;
        Float_t  DisttoBC;
        Short_t  CellIds[20];

    ClassDef(AliAODClusterreduced, 1); // class required for event mixing
};

#endif