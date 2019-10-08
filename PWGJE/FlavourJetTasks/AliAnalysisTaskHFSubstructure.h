/// \class AliAnalysisTaskDmesonJets
/// \brief Analysis task for D meson jets
///
/// This task selects D meson candidates according to predefined cuts,
/// then runs a jet finder to reconstruct the jets that contain
/// the D meson candidates.
///
/// The main output is stored in a THnSparse histogram or in a TTree.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Oct 13, 2017

#ifndef ALIANALYSISTASKHFSUBSTRUCTURE_H
#define ALIANALYSISTASKHFSUBSTRUCTURE_H

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
class AliRDHFCuts;
class AliAODEvent;
class AliAODRecoDecay;
class AliAODRecoDecayHF2Prong;
class AliAODRecoCascadeHF;
class AliVParticle;
class AliAODMCParticle;
class AliHFAODMCParticleContainer;
class AliHFTrackContainer;
class AliParticleContainer;
class AliClusterContainer;
class THnSparse;
class AliFJWrapper;
class THashList;
class TTree;
class AliEMCALGeometry;
class TRandom;
class AliRhoParameter;
class TFile;

//C++
#include <exception>
#include <list>
#include <vector>
#include <map>

#include "AliTLorentzVector.h"
#include "THistManager.h"

#include "AliAnalysisTaskEmcal.h"
#include "AliJetContainer.h"
#include "TFile.h"
#include "AliRDHFCutsD0toKpi.h"


class AliAnalysisTaskHFSubstructure : public AliAnalysisTaskEmcal
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

 typedef AliJetContainer::EJetType_t EJetType_t;
 typedef AliJetContainer::EJetAlgo_t EJetAlgo_t;
 typedef AliJetContainer::ERecoScheme_t ERecoScheme_t;

 enum JetShapeType {
   kData   = 0,  
   kDetSignal   = 1,  
   kDetBackground   = 2,  
   kDetReflection   = 3, 
   kTrueDet   = 4,  
   kTrue   = 5,
   kDataInclusive   = 6,
   kDet   = 7, 
 };
 enum EMesonOrigin_t {
   kUnknownQuark = BIT(0),
   kFromDown     = BIT(1),
   kFromUp       = BIT(2),
   kFromStrange  = BIT(3),
   kFromCharm    = BIT(4),
   kFromBottom   = BIT(5),
   kFromTop      = BIT(6),
   kFromGluon    = BIT(7),
   kAnyOrigin    = kUnknownQuark | kFromDown | kFromUp | kFromStrange | kFromCharm | kFromBottom | kFromTop | kFromGluon
 };

 enum EMesonDecayChannel_t {
   kAnyDecay            = 0,
   kUnknownDecay        = BIT(0),
   kDecayD0toKpi        = BIT(1),
   kDecayDStartoKpipi   = BIT(2)
 };

 enum TreeSize {
   nVar = 10,
   nVar_Splittings =10
  };

 //enum ECandidateType_t  { kD0toKpi, kDstartoKpipi, kD0toKpiLikeSign };
 enum ECandidateType_t  {kD0toKpi};

 AliAnalysisTaskHFSubstructure();
 AliAnalysisTaskHFSubstructure(const char *name);
 virtual ~AliAnalysisTaskHFSubstructure();




 AliAnalysisTaskHFSubstructure* AddTaskAliAnalysisTaskHFSubstructure(const char * ntracksData,
								     const char * ntracksTrue,
								     const char * ntracksDet,
								     const Double_t R,
								     AliAnalysisTaskHFSubstructure::ECandidateType_t ECandidateType = AliAnalysisTaskHFSubstructure::kD0toKpi,
								     AliAnalysisTaskHFSubstructure::JetShapeType jetShapeType = AliAnalysisTaskHFSubstructure::kData);




 

 void                                UserCreateOutputObjects();
 void                                Terminate(Option_t *option);

 JetShapeType                        fJetShapeType;               
 ECandidateType_t                    fECandidateType;


 void SetJetShapeType(JetShapeType tJetShapeType)                 {fJetShapeType = tJetShapeType;}
 void SetECandidateType_t (ECandidateType_t ECandidateType)       {fECandidateType = ECandidateType;}

 void SetIncludeInclusive(Bool_t IncludeInclusive)                {fIncludeInclusive = IncludeInclusive;}
 Bool_t GetIncludeInclusive()                                     {return fIncludeInclusive;}
 void SetIsBDecay(Bool_t IsBDecay)                                {fIsBDecay = IsBDecay;}
 Bool_t GetIsBDecay()                                             {return fIsBDecay;} 
 void SetBranchName(TString BranchName)                           {fBranchName = BranchName;}
 TString GetBranchName()                                          {return fBranchName;}
 void SetCutsType(TString CutsType)                               {fCutsType = CutsType;}
 TString GetCutsType()                                            {return fCutsType;}
 void SetCandidatePDG(Int_t CandidatePDG)                         {fCandidatePDG = CandidatePDG;}
 TString GetCandidatePDG()                                        {return fCandidatePDG;}
 void SetJetRadius(Double_t JetRadius)                            {fJetRadius = JetRadius;}
 Double_t GetJetRadius()                                          {return fJetRadius;}
 void SetJetMinPt(Double_t JetMinPt)                              {fJetMinPt = JetMinPt;}
 Double_t GetJetMinPt()                                           {return fJetMinPt;}
 void SetTrackingEfficiency(Double_t TrackingEfficiency)         {fTrackingEfficiency = TrackingEfficiency;}
 Double_t GetTrackingEfficiency()                                 {return fTrackingEfficiency;}
 void SetRejectedOrigin(UInt_t RejectedOrigin)                    {fRejectedOrigin = RejectedOrigin;}
 void SetAcceptedDecay(UInt_t AcceptedDecay)                      {fAcceptedDecay = AcceptedDecay;}
 void SetRejectISR(Bool_t RejectISR)                              {fRejectISR = RejectISR;}
 Bool_t GetRejectISR()                                            {return fRejectISR;}
 void SetPromptReject(Bool_t PromptReject)                        {fPromptReject = PromptReject;}
 Bool_t GetPromptRejectR()                                        {return fPromptReject;}
 void SetAlienConnect(Bool_t AlienConnect)                        {fAlienConnect = AlienConnect;}
 void SetCuts(TString CutsFileName)                               {TFile *fCutsFile = TFile::Open(CutsFileName);
                                                                   TString cutsname="D0toKpiCuts";
                                                                   if (fCutsType!="") cutsname += TString::Format("_%s", fCutsType.Data()); 
                                                                   fRDHFCuts = dynamic_cast<AliRDHFCuts*>(fCutsFile->Get(cutsname));}
 


 
 protected:

 Bool_t                              RetrieveEventObjects();
 Bool_t                              Run()                 ;
 Bool_t                              FillHistograms()      ;
 Bool_t                              fIncludeInclusive     ;
 Bool_t                              fIsBDecay             ;
 Bool_t                              fRejectISR            ;
 Bool_t                              fPromptReject         ;
 Bool_t                              fAlienConnect         ;
 
 TString                            fBranchName            ; 
 TString                            fCutsType;
 Int_t                              fCandidatePDG          ; 
 UInt_t                             fRejectedOrigin        ; 
 UInt_t                             fAcceptedDecay         ; 

 Double_t                           fJetRadius             ;
 Double_t                           fJetMinPt              ;
 Double_t                           fTrackingEfficiency    ;
 Double_t                           fShapesVar[nVar]       ;

   

 TClonesArray                      *fCandidateArray        ; 
   
 AliAODEvent                       *fAodEvent              ; 

 AliRDHFCuts                       *fRDHFCuts              ; 
 AliFJWrapper                      *fFastJetWrapper        ; //!<! 
 AliFJWrapper                      *fFastJetWrapper_Truth  ; //!<! 

 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_DeltaR;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_DeltaR_Truth;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_Zg;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_Zg_Truth;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_LeadingSubJetpT;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_LeadingSubJetpT_Truth;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_HardestSubJetD0;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_HardestSubJetD0_Truth;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_RadiatorE;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_RadiatorE_Truth;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_RadiatorpT;
 std::vector<std::vector<Double_t>>            fShapesVar_Splittings_RadiatorpT_Truth;
 TTree                               *fTreeResponseMatrixAxis;
 TTree                               *fTreeSplittings;




 TH1D                                *fhEvent;
   

 friend class AliAnalysisTaskDmesonJetsSub;
 friend class OutputHandler;

 private:

 AliAnalysisTaskHFSubstructure(const AliAnalysisTaskHFSubstructure&);            
 AliAnalysisTaskHFSubstructure &operator=(const AliAnalysisTaskHFSubstructure&); 

 ClassDef(AliAnalysisTaskHFSubstructure, 3)
    
   };

#endif
