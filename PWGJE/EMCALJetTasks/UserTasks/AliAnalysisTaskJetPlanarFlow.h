/// \class AliAnalysisTaskJetPlanarFlow
/// \brief Analysis task for planar flow and wake analysis
///
/// The main output is stored in a TTree.
///
/// \author Nima Zardoshti

#ifndef ALIANALYSISTASKJETPLANARFLOW_H
#define ALIANALYSISTASKJETPLANARFLOW_H

/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

class TClonesArray;
class AliAODEvent;
class AliVParticle;
class AliAODMCParticle;
class AliParticleContainer;
class AliClusterContainer;
class THnSparse;
class THashList;
class TTree;
class AliEMCALGeometry;
class TRandom;
class AliRhoParameter;
class TFile;
class AliJetContainer;
class AliEmcalJetFinder;
class AliFJWrapper;
class TH1;
class TH2;
class TH3;
class TArrayI;
class AliAnalysisManager;

//C++
#include <exception>
#include <list>
#include <vector>
#include <map>

#include "AliTLorentzVector.h"
#include "THistManager.h"

#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliFJWrapper.h"
#include "AliJetContainer.h"
#include "TFile.h"
#include "AliAnalysisTaskEmcalJet.h"



class AliAnalysisTaskJetPlanarFlow : public AliAnalysisTaskEmcalJet
{
 public:

  class AliEventNotFound : public std::exception
   {
   public:
     AliEventNotFound(const std::string& class_name, const std::string& method_name);
#if !(defined(__CINT__) || defined(__MAKECINT__))
     const char* what() const noexcept;
#endif


   private:
     std::string  fClassName        ; 
     std::string  fAccessMethodName ; 
     std::string  fWhat             ;
   };


 enum JetAnalysisType {
   kData   = 0,  
   kDet   = 1,   
   kTrue   = 2,
   kTrueDet   = 3, 
   kEmb   = 4
 };

 enum JetSubType {
   kNoSub   = 0,  
   kAreaSub   = 1
 };

 enum TreeSize {
   nVar = 12,
   nVar_Particles =16
  };


 AliAnalysisTaskJetPlanarFlow();
 AliAnalysisTaskJetPlanarFlow(const char *name);
 virtual ~AliAnalysisTaskJetPlanarFlow();



/*
AliAnalysisTaskJetPlanarFlow* AddTaskJetPlanarFlow(
                    const char * njet1,
                    const char * ntracks1,
								    const char * njets2,
                    const char * ntracks2,
								    const char * njets3,
                    const char * ntracks3, 
								    const Double_t R,
								    const char * nrhoBase, 
								    const char *type,				      
								    const char *CentEst,
								    Int_t       pSel,
								    AliAnalysisTaskJetPlanarFlow::JetAnalysisType jetAnalysisType = AliAnalysisTaskJetPlanarFlow::kData,
								    AliAnalysisTaskJetPlanarFlow::JetSubType jetSubType = AliAnalysisTaskJetPlanarFlow::kNoSub);
                    */

Float_t RelativePhi(Float_t mphi, Float_t vphi);
void SetTree(AliEmcalJet *jet, AliJetContainer *jetContainer, AliTrackContainer *trackContainer, Float_t jetPt, Int_t level);

void UserCreateOutputObjects();
void Terminate(Option_t *option);

void SetJetAnalysisType(JetAnalysisType tJetAnalysisType) { fJetAnalysisType = tJetAnalysisType; }
void SetJetSubType(JetSubType tJetSubType) { fJetSubType = tJetSubType; }

void SetJetRadius(Float_t JetRadius) { fJetRadius = JetRadius; }
Float_t GetJetRadius() { return fJetRadius; }
void SetMinJetPt(Float_t MinJetPt) { fMinJetPt = MinJetPt; }
Float_t GetMinJetPt() { return fMinJetPt; }
void SetMaxJetPt(Float_t MaxJetPt) { fMaxJetPt = MaxJetPt; }
Float_t GetMaxJetPt() { return fMaxJetPt; }
void SetMinJetConstiteuntPAPt(Float_t MinJetConstiteuntPAPt) { fMinJetConstiteuntPAPt = MinJetConstiteuntPAPt; }
Float_t GetMinJetConstiteuntPAPt() { return fMinJetConstiteuntPAPt; }
void SetMaxJetConstiteuntPAPt(Float_t MaxJetConstiteuntPAPt) { fMaxJetConstiteuntPAPt = MaxJetConstiteuntPAPt; }
Float_t GetMaxJetConstiteuntPAPt() { return fMaxJetConstiteuntPAPt; }
void SetTrackingEfficiency(Bool_t TrackingEfficiency)                              {fTrackingEfficiency = TrackingEfficiency;}
Bool_t GetTrackingEfficiency()                                           {return fTrackingEfficiency;}
void SetSharedFractionPtMin(Float_t SharedFractionPtMin) { fSharedFractionPtMin = SharedFractionPtMin; }
Float_t GetSharedFractionPtMin() { return fSharedFractionPtMin; }
void SetCentralitySelection(Bool_t CentralitySelection)                              {fCentralitySelection = CentralitySelection;}
Bool_t GetCentralitySelection()                                           {return fCentralitySelection;}
void SetCentralityMin(Float_t CentralityMin)                              {fCentralityMin = CentralityMin;}
Float_t GetCentralityMin()                                           {return fCentralityMin;}
void SetCentralityMax(Float_t CentralityMax)                              {fCentralityMax = CentralityMax;}
Float_t GetCentralityMax()                                           {return fCentralityMax;}
void SetFillNsubjettiness(Bool_t FillNsubjettiness)                              {fFillNsubjettiness = FillNsubjettiness;}
void SetFillDeltaR(Bool_t FillDeltaR)                              {fFillDeltaR = FillDeltaR;}
void SetWithoutRotations(Bool_t DoRotations)                              {fDoRotations = DoRotations;}

 
 


 
 protected:

 Bool_t                              RetrieveEventObjects();
 Bool_t                              Run()                 ;
 Bool_t                              FillHistograms()      ;

 JetAnalysisType fJetAnalysisType;
 JetSubType fJetSubType;
 



 Float_t                           fJetRadius             ;
 Float_t                           fMinJetPt              ;
 Float_t                           fMaxJetPt              ;
 Float_t                           fMinJetConstiteuntPAPt ;
 Float_t                           fMaxJetConstiteuntPAPt ;
 Float_t                          fSharedFractionPtMin    ;
 Float_t                           fTrackingEfficiency    ;
 Bool_t                           fCentralitySelection    ;
 Float_t                           fCentralityMin         ;
 Float_t                           fCentralityMax         ;
 Bool_t                            fFillNsubjettiness     ;
 Bool_t                            fFillDeltaR            ;
 Bool_t                            fDoRotations;
 Float_t                           fShapesVar[nVar];

 TRandom3                          fRandom                ;

 std::vector<Int_t>               fJetConstituentLabels    ;
 std::vector<Float_t>             fShapesVar_Particles_E;
 std::vector<Float_t>             fShapesVar_Particles_E_Truth;
 std::vector<Float_t>             fShapesVar_Particles_pT;
 std::vector<Float_t>             fShapesVar_Particles_pT_Truth;
 std::vector<Float_t>             fShapesVar_Particles_Phi;
 std::vector<Float_t>             fShapesVar_Particles_Phi_Truth;
 std::vector<Float_t>             fShapesVar_Particles_Theta;
 std::vector<Float_t>             fShapesVar_Particles_Theta_Truth;
 std::vector<Float_t>             fShapesVar_Particles_InJet;
 std::vector<Float_t>             fShapesVar_Particles_InJet_Truth;
 std::vector<Float_t>             fShapesVar_Particles_DeltaR;
 std::vector<Float_t>             fShapesVar_Particles_DeltaR_Truth;
 std::vector<Float_t>             fShapesVar_Particles_NRPhi;
 std::vector<Float_t>             fShapesVar_Particles_NRPhi_Truth;
 std::vector<Float_t>             fShapesVar_Particles_Eta;
 std::vector<Float_t>             fShapesVar_Particles_Eta_Truth;
 TTree                            *fTreeJet;
 TTree                            *fTreeParticles;




 TH1D                             *fhEvent;
   


 friend class OutputHandler;

 private:

 AliAnalysisTaskJetPlanarFlow(const AliAnalysisTaskJetPlanarFlow&);            
 AliAnalysisTaskJetPlanarFlow &operator=(const AliAnalysisTaskJetPlanarFlow&); 

 ClassDef(AliAnalysisTaskJetPlanarFlow, 2)
    
   };

#endif
