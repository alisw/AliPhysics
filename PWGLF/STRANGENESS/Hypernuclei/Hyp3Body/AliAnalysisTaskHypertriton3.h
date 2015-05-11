#ifndef ALIANALYSISTASKHYPERTRITON3_H
#define ALIANALYSISTASKHYPERTRITON3_H


/**************************************************************************
 *                                                                        *                                                                         
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



///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskHypertriton3 class
// analysis task for the study of the production of hypertriton
// which decays in 3 prongs: d+p+pi^-
// This task is optimized for ESDs.root
//
// Author: 
// S. Trogolo, trogolo@to.infn.it
///////////////////////////////////////////////////////////////////////////

#include <TROOT.h>

#include "AliAnalysisTaskSE.h"

class TChain;
class TH1F;
class TH2F;
class TList;
class TObjArray;
class TTree;

class AliESDEvent;
class AliESDtrackCuts;
class AliESDVertex;
class AliPIDResponse;
class AliVertexerTracks; 

class AliAnalysisTaskHypertriton3 : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskHypertriton3();
  virtual ~AliAnalysisTaskHypertriton3();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t*);
  virtual void   Terminate(Option_t*);

  void SetReadMC(Bool_t flag = kTRUE) {fMC = flag;}
  void SetFillTree(Bool_t outTree = kFALSE) {fFillTree = outTree;}

  Double_t GetDCAcut(Int_t part, Double_t dca) const;

 private:

  AliESDEvent        *fESDevent;                   ///< ESD event
  AliESDtrackCuts    *fESDtrackCuts;               ///< First set of ESD track cuts
  AliESDtrackCuts    *fESDtrackCutsV0;             ///< Track cuts applied only to \f$\pi^{-}\f$ and \f$\pi^{+}\f$ candidate
  AliESDVertex       *fPrimaryVertex;              //!<! Primary vertex of the current event
  AliPIDResponse     *fPIDResponse;                //!<! PID response class
  AliVertexerTracks  *fVertexer;                   //!<! Secondary vertex reconstructed with three candidate tracks

  TObjArray          *fTrkArray;                   //!<! Array containing the three tracks candidated to the secondary vertex reconstruction
  
  //Variables
  Bool_t             fMC;                          ///< variables for MC selection
  Bool_t             fFillTree;                    ///< variables to fill the Tree
  Float_t            fCentrality;                  ///< Centrality class
  Float_t            fCentralityPercentile;        ///< Centrality percentile

  //Cut variables
  Double_t           fDCAPiPVmin;                  ///< Cut on Min DCA of \f$\pi\f$ from primary vertex
  Double_t           fCosPointingAngle;            ///< Cut on Cosine of the pointing angle
  Double_t           fDecayLength;                 ///< Cut on Decay length
  Double_t           fPtMother;                    ///< Cut on mother reconstructed \f$p_{T}\f$
  Double_t           fDCAPiSVxymax;                ///< Cut on \f$\pi DCA_{xy}\f$ from reconstructed secondary vertex
  Double_t           fDCAPiSVzmax;                 ///< Cut on \f$\pi DCA_{z}\f$ from reconstructed secondary vertex
  Double_t           fDCAProSVmax;                 ///< Cut on proton DCA from reconstructed secondary vertex
  Double_t           fDCADeuSVmax;                 ///< Cut on deuteron DCA from reconstructed secondary vertex
  Double_t           fDCAdp;                       ///< Cut DCA deuteron-proton
  Double_t           fDCApip;                      ///< Cut DCA pion-proton
  Double_t           fDCAdpi;                      ///< Cut DCA deuteron-pion
  
  //Output list
  TList              *fOutput;                     ///< Output list

  //Histograms
  TH1F               *fHistCount;                  //!<! Total number of events
  TH1F               *fHistCentralityClass;        //!<! Event statistic per centrality class
  TH1F               *fHistCentralityPercentile;   //!<! Event statistic per centrality percentile
  TH1F               *fHistTrigger;                //!<! Event trigger statistic
  TH1F               *fHistMultiplicity;           //!<! Event multiplicity
  TH1F               *fHistZPrimaryVtx;            //!<! Primary vertex Z coordinate
  TH1F               *fHistXPrimaryVtx;            //!<! Primary vertex X coordinate
  TH1F               *fHistYPrimaryVtx;            //!<! Primary vertex Y coordinate
  TH1F               *fHistChi2perTPCcluster;      //!<! TPC \f$\chi^{2}/NDF\f$ tracks distribution


  //PID
  //--> TPC
  TH2F               *fHistTPCpid;                 //!<! TPC dE/dx vs \f$p_{TPC}\f$
  TH2F               *fHistTPCdeusignal;           //!<! TPC PID: deuteron candidates
  TH2F               *fHistTPCprosignal;           //!<! TPC PID: proton candidates
  TH2F               *fHistTPCpionsignal;          //!<! TPC PID: \f$\pi^{-}\f$ candidates
  TH2F               *fHistTPCantideusignal;       //!<! TPC PID: anti-deuteron candidates
  TH2F               *fHistTPCantiprosignal;       //!<! TPC PID: anti-proton candidates
  TH2F               *fHistTPCpionplussignal;      //!<! TPC PID: \f$\pi^{+}\f$ candidates

  //--> TOF
  TH2F               *fHistTOFsignal;              //!<! TOF \f$\beta\f$ vs \f$p_{TPC}\f$
  //TH2F               *fHistTOFdeusignal;           //!<! TOF PID: deuteron candidates
  //TH2F               *fHistTOFprosignal;           //!<! TOF PID: proton candidates
  //TH2F               *fHistTOFantideusignal;       //!<! TOF PID: anti-deuteron candidates
  //TH2F               *fHistTOFantiprosignal;       //!<! TOF PID: anti-proton candidates
  TH1F               *fHistTOFdeumass;             //!<! TOF mass of deuteron identified with TPC
  TH1F               *fHistTOFpromass;             //!<! TOF mass of proton identified with TPC

  //Candidate combination
  // Data and MC histograms
  TH1F               *fHistpionTPCcls;             //!<! TPC clusters distribution of candidate \f$\pi\f$
  TH1F               *fHistpTpion;                 //!<! \f$\p^{T}\f$ distribution of candidate \f$\pi\f$
  TH2F               *fHistCorrDCAdprimary;        //!<! Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ deuteron-primary vertex 
  TH2F               *fHistCorrDCApprimary;        //!<! Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ proton-primary vertex
  TH2F               *fHistCorrDCApiprimary;       //!<! Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ pion-primary vertex
  TH1F               *fHistDCApiprimary;           //!<! DCA pion-primary vertex distribution
  TH1F               *fHistDCAdeupro;              //!<! DCA deuteron-proton distribution
  TH1F               *fHistDCApiondeu;	           //!<! DCA pion-deuteron distribution
  TH1F               *fHistDCApionpro;             //!<! DCA pion-proton distribution
  TH1F               *fHistZDecayVtx;              //!<! Reco secondary vertex Z coordinate
  TH1F               *fHistXDecayVtx;              //!<! Reco secondary vertex X coordinate
  TH1F               *fHistYDecayVtx;              //!<! Reco secondary vertex Y coordinate
  TH1F               *fHistDCAXYdeuvtx;            //!<! \f$DCA_{xy}\f$ candidate deuteron-secondary vertex
  TH1F               *fHistDCAZdeuvtx;             //!<! \f$DCA_{z}\f$ candidate deuteron-secondary vertex
  TH1F               *fHistDCAXYprovtx;            //!<! \f$DCA_{xy}\f$ candidate proton-secondary vertex
  TH1F               *fHistDCAZprovtx;             //!<! \f$DCA_{z}\f$ candidate proton-secondary vertex
  TH1F               *fHistDCAXYpionvtx;           //!<! \f$DCA_{xy}\f$ candidate pion-secondary vertex
  TH1F               *fHistDCAZpionvtx;            //!<! \f$DCA_{z}\f$ candidate pion-secondary vertex
  TH1F               *fHistDecayLengthH3L;         //!<! Decay length distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistCosPointingAngle;       //!<! Cosine of pointing angle distribution of candidate mother particle
  TH1F               *fHistMassHypertriton;        //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$
  // MC only histograms
  TH1F               *fHistParticle;               //!<! *(MC only)* Reconstructed particles distribution per species through PDGCode cross-check
  TH1F               *fHistpionTPCclsMCt;          //!<! *(MC only)* TPC clusters distribution of candidate \f$\pi\f$ through PDGCode cross-check
  TH1F               *fHistpTpionMCt;              //!<! *(MC only)* \f$\p^{T}\f$ distribution of \f$\pi\f$ identified with PDGCode
  TH1F               *fHistpTproMCt;               //!<! *(MC only)* \f$\p^{T}\f$ distribution of proton identified with PDGCode
  TH1F               *fHistpTdeuMCt;               //!<! *(MC only)* \f$\p^{T}\f$ distribution of deuteron identified with PDGCode
  TH2F               *fHistCorrDCAdprimaryMCt;     //!<! *(MC only)* Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ deuteron(PDGCode identified)-primary vertex 
  TH2F               *fHistCorrDCApprimaryMCt;     //!<! *(MC only)* Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ proton(PDGCode identified)-primary vertex
  TH2F               *fHistCorrDCApiprimaryMCt;    //!<! *(MC only)* Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ pion(PDGCode identified)-primary vertex
  TH1F               *fHistDCApiprimaryMCt;        //!<! *(MC only)* DCA pion(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCApprimaryMCt;         //!<! *(MC only)* DCA proton(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAdprimaryMCt;         //!<! *(MC only)* DCA deuteron(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAdeuproMCt;           //!<! *(MC only)* DCA deuteron-proton distribution identified with PdgCode
  TH1F               *fHistDCApiondeuMCt;          //!<! *(MC only)* DCA pion-deuteron distribution identified with PdgCode
  TH1F               *fHistDCApionproMCt;          //!<! *(MC only)* DCA pion-proton distribution identified with PdgCode
  TH1F               *fHistZDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex Z coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistXDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex X coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistYDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex Y coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistDCAXYdeuvtxMCt;         //!<! *(MC only)* \f$DCA_{xy}\f$ candidate deuteron(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAZdeuvtxMCt;          //!<! *(MC only)* \f$DCA_{z}\f$ candidate deuteron(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAXYprovtxMCt;         //!<! *(MC only)* \f$DCA_{xy}\f$ candidate proton(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAZprovtxMCt;          //!<! *(MC only)* \f$DCA_{z}\f$ candidate proton(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAXYpionvtxMCt;        //!<! *(MC only)* \f$DCA_{xy}\f$ candidate pion(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAZpionvtxMCt;         //!<! *(MC only)* \f$DCA_{z}\f$ candidate pion(PDGCode identified)-secondary vertex
  TH1F               *fHistMassHypertritonMCt;     //!<! *(MC only)* Invariant mass distribution of reconstructed \f$H^{3}_{\Lambda}\f$ - daughters particles identified with PDGCode

  //TTree
  TTree              *fTTree;                      //!<! Tree used for local tests and cross-check
  // Deuteron
  Float_t            fTchi2NDFdeu;
  UShort_t           fTPCclsdeu;
  UShort_t           fTPCclsPIDdeu;
  Float_t            fTpTPCdeu;
  Float_t            fTpTdeu;
  Float_t            fTpdeu;
  Float_t            fTTPCnsigmadeu;
  Float_t            fTTOFmassdeu;
  Float_t            fTDCAXYdeuprvtx;
  Float_t            fTDCAZdeuprvtx;
  // Proton
  Float_t            fTchi2NDFpro;
  UShort_t           fTPCclspro;
  UShort_t           fTPCclsPIDpro;
  Float_t            fTpTPCpro;
  Float_t            fTpTpro;
  Float_t            fTppro;
  Float_t            fTTPCnsigmapro;
  Float_t            fTTOFmasspro;
  Float_t            fTDCAXYproprvtx;
  Float_t            fTDCAZproprvtx;
  // Pion
  Float_t            fTchi2NDFpion;
  UShort_t           fTPCclspion;
  UShort_t           fTPCclsPIDpion;
  Float_t            fTpTPCpion;
  Float_t            fTpTpion;
  Float_t            fTppion;
  Float_t            fTTPCnsigmapion;
  Float_t            fTDCAXYpioprvtx;
  Float_t            fTDCAZpioprvtx;
  // Triplet
  Float_t            fTDCAdp;
  Float_t            fTDCAdpi;
  Float_t            fTDCAppi;
   
  Float_t            fTDCAXYdvtx;
  Float_t            fTDCAZdvtx;
  Float_t            fTDCAXYpvtx;
  Float_t            fTDCAZpvtx;
  Float_t            fTDCAXYpivtx;
  Float_t            fTDCAZpivtx;
   
  Float_t            fTDecayLength;
  Float_t            fTCosPA;
  Float_t            fTInvariantMass;
 
  
  AliAnalysisTaskHypertriton3(const AliAnalysisTaskHypertriton3&); // not implemented
  AliAnalysisTaskHypertriton3& operator=(const AliAnalysisTaskHypertriton3&); // not implemented
  
  ClassDef(AliAnalysisTaskHypertriton3, 1); // analysisclass
  
};

#endif
