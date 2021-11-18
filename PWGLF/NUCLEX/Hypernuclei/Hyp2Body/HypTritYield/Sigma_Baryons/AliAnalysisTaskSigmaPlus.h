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

class AliAnalysisTaskSigmaPlus : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskSigmaPlus();
                                AliAnalysisTaskSigmaPlus(const char *name);
        virtual                 ~AliAnalysisTaskSigmaPlus();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        
    private:

        //Functions for Filling Histograms
        void FillHistogram(const char* key,Double_t x) const;                         // Fill Histogram with the name "key"
        void FillHistogram(const char* key,Double_t x, Double_t y) const;             // Histogram must be in TList "fOutputList"
        void FillHistogram(const char* key,Double_t x, Double_t y, Double_t z) const; // 1, 2 and 3 Dimensions is supported (overload)

        //Functions for Event Processing
        void FillProtonArray() const;                                                 // Select (Anti-)Protons and fill arrays (fProtonArray)
        void FillElectronArray() const;                                               // Select Electrons and fill arrays (ElectronCandID)
        void ProcessMCParticles() const;                                              // Fill MC Histograms
        void FillV0PhotonArray() const;                                               // Select V0 Photons and fill arrays (fConvPhotonArray)
        void PairElectrons() const;                                                   // Pair selected electrons using KFParticle  
        void PairPhotonandElectron();                                                 // Pair V0 Photon with a single electron  
        void ReconstructParticles();                                                  // Reconstruct Pi0 and Sigma Plus from selected Protons and Photons  
        void ReconstructParticlesfromee() const;                                      // Reconstruct Pi0 (and Sigma Plus) from selected e+ e- pairs with KF  

        TList*                fOutputList;                //!<! Output list which contains all Histograms

        AliAODEvent*          aodEvent;                   //!<! AOD Event to be processed
        AliMCEvent*           mcEvent;                    //!<! MC Event to be processed   
        TClonesArray*         AODMCTrackArray;            //!<! TClonesArray containing the MC Particles  

        AliPIDResponse*       fPIDResponse;               //!<! PID response object 

        TTree*                fSigmaCandTree;             //!  Tree with Sigma candidates
        TTree*                f4PartSigmaCandTree;        //!  Tree with Sigma candidates
        TTree*                f3PartPi0CandTree;          //!  Tree with Sigma candidates

        //Particle Arrays
        std::vector<int>*     fProtonArray;               //!<! Proton candidates
        std::vector<int>*     fAntiProtonArray;           //!<! Anti-Proton candidates
        std::vector<int>*     fConvPhotonArray;           //!<! V0 Photon candidates
        std::vector<int>*     fElectronArray;             //!<! Electron candidates
        std::vector<int>*     fPositronArray;             //!<! Positron candidates
        std::vector<int>*     fElectronPairArray;         //!<! Electron from photon candidates
        std::vector<int>*     fPositronPairArray;         //!<! Positron from photon candidates

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

        //FillProtonArray Cuts
        Double_t              fMaxProtEta = 0.9;    
        Double_t              fMinTPCClustProt = 70;
        Double_t              fMaxNsigProtTPC = 3;
        Double_t              fMaxNsigProtTOF = 3;  
        Double_t              fMaxpOnlyTPCPID = 0.75;
        Double_t              fMinProtpt = 0.5; 
        Double_t              fMaxProtpt = 10; //was 4  

        //FillElectronArray Cuts
        Double_t              fMaxElecEta = 0.9;    
        Double_t              fMinTPCClustElec = 60;
        Double_t              fMaxNsigElecTPC = 3;
        Double_t              fMaxNsigElecTOF = 3;  
        Double_t              fNsigHadronreject = 2;   

        //ProcessMCParticles Cuts
        Double_t              fMaxMCEta = 0.9;    

        //FillV0PhotonArray Cuts
        Double_t              fMaxDaughtEta = 0.9;    
        Double_t              fMinTPCClustDaught = 60;
        Double_t              fMaxNsigDaughtTPC = 3;
        Double_t              fMaxalpha = 0.9;
        Double_t              fMaxqt = 0.04;
        Double_t              fMaxopenangle = 0.2;
        Double_t              fMaxdeltatheta = 0.1;
        Double_t              fMaxphotonmass = 0.05;

        //PairElectrons Cuts
        Double_t              fMaxPairalpha = 0.9;
        Double_t              fMaxPairqt = 0.04;
        Double_t              fMaxPairopenangle = 0.2;
        Double_t              fMaxPairDCA = 2;
        Double_t              fMaxPairphotonmass = 0.06;

        //PairPhotonandElectron Cuts
        Double_t              fMinSingleElectronPt = 0.25;
        Double_t              fMinPi0mass3Part = 0.12;
        Double_t              fMaxPi0mass3Part = 0.15;
        Double_t              fMaxElecPt3Part = 0.9;
        Double_t              fMaxSigmaPA3Part = 0.175;
        Double_t              fMinPi0DCAtoPV3Part = 1;
        Double_t              fMaxProtPi0DCA3Part = 2;
        Double_t              fMaxPhotElecDCA = 2;

        //ReconstructParticles Cuts
        Double_t              fMinPi0mass = 0.12;
        Double_t              fMaxPi0mass = 0.15;
        Double_t              fMaxElecPt = 0.9;
        Double_t              fMaxSigmaPA = 0.175;
        Double_t              fMinPi0DCAtoPV = 1;
        Double_t              fMaxProtPi0DCA = 2;
        Double_t              fMaxPhotPhotDCA = 2;
        Double_t              fMindcazProt = 0;
        Double_t              fMindcaxyProt = 0; 

        //Branches of TTree "fSigmaCandTree"
        Bool_t                fIsMCSigma;
        Float_t               fInvSigMass;               
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
        Float_t               fInvPi0Mass;        
        Float_t               fPi0Px;        
        Float_t               fPi0Py;        
        Float_t               fPi0Pz;
        Float_t               fPi0DaughtPt1;
        Float_t               fPi0DaughtPt2;
        Float_t               fPi0DaughtPt3;
        Float_t               fPi0DaughtPt4;
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
        //End of Branches of "fSigmaCandTree"

        //Branches of TTree "f4PartSigmaCandTree"
        Bool_t                f4PartIsMCSigma;
        Float_t               f4PartInvSigMass;               
        Float_t               f4PartSigCharge;        
        Float_t               f4PartSigPx;        
        Float_t               f4PartSigPy;        
        Float_t               f4PartSigPz;        
        Float_t               f4PartPrimVertX;        
        Float_t               f4PartPrimVertY;        
        Float_t               f4PartPrimVertZ;        
        Float_t               f4PartSigDecayVertX;        
        Float_t               f4PartSigDecayVertY;        
        Float_t               f4PartSigDecayVertZ;        
        Float_t               f4PartInvPi0Mass;        
        Float_t               f4PartPi0Px;        
        Float_t               f4PartPi0Py;        
        Float_t               f4PartPi0Pz;
        Float_t               f4PartPi0DaughtPt1;
        Float_t               f4PartPi0DaughtPt2;
        Float_t               f4PartPi0DaughtPt3;
        Float_t               f4PartPi0DecayVertX;        
        Float_t               f4PartPi0DecayVertY;        
        Float_t               f4PartPi0DecayVertZ;
        Float_t               f4PartPi0PhotPhotDCA;
        Float_t               f4PartProtonPx;        
        Float_t               f4PartProtonPy;        
        Float_t               f4PartProtonPz;        
        Float_t               f4PartProtonDCAtoPVxy;        
        Float_t               f4PartProtonDCAtoPVz;        
        Float_t               f4PartProtonPi0DCA;        
        //End of Branches of "f4PartSigmaCandTree"

        //Branches of TTree "f3PartPi0CandTree"
        Float_t               fInv3PartPi0Mass;
        Float_t               f3PartPi0Px;        
        Float_t               f3PartPi0Py;        
        Float_t               f3PartPi0Pz;
        Float_t               f3PartPi0DaughtPt1;
        Float_t               f3PartPi0DaughtPt2;
        Float_t               f3PartPi0DaughtPt3;
        Float_t               f3PartPi0DecayVertX;        
        Float_t               f3PartPi0DecayVertY;        
        Float_t               f3PartPi0DecayVertZ;
        Float_t               f3PartPi0PhotPhotDCA;
        Float_t               f3PartPi0PhotConvx;        
        Float_t               f3PartPi0PhotConvy;        
        Float_t               f3PartPi0PhotConvz;
        Float_t               f3PartPrimVertX;        
        Float_t               f3PartPrimVertY;        
        Float_t               f3PartPrimVertZ;        
        Float_t               f3PartPi0ElecPosx;        
        Float_t               f3PartPi0ElecPosy;        
        Float_t               f3PartPi0ElecPosz;
        //End of Branches of "f3PartPi0CandTree"

        AliAnalysisTaskSigmaPlus(const AliAnalysisTaskSigmaPlus&); // not implemented
        AliAnalysisTaskSigmaPlus& operator=(const AliAnalysisTaskSigmaPlus&); // not implemented

        ClassDef(AliAnalysisTaskSigmaPlus, 1);
};

#endif
