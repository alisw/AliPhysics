/**************************************************************************
 * Author : Benjamin Dönigus (benjamin.doenigus@cern.ch)                  *
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskHdibaryonLPpi class
//          used to search for the H-Dibaryon in weak 
//          (Lambda Proton Pion) and strong (Lambda Lambda) decays
//-----------------------------------------------------------------


#ifndef ALIANALYSISTASKHDIBARYONLPPI_H
#define ALIANALYSISTASKHDIBARYONLPPI_H

// analysis task for the H-Dibaryon analysis
// Author: Benjamin Doenigus

#include "AliAnalysisTaskSE.h"
//#include "AliESDv0Cuts.h"

class TH1F;
class AliPIDResponse;
class AliESDEvent;
class TFile;
class TH2F;
class TH3F;
class THnSparse;
class AliESDtrackCuts;
class AliESDv0Cuts;
class AliMCEvent;
class TList;
class TParticle;
class TLegend;
class AliESDpid;
class AliHFEpid;
class AliHFEtools;

class AliMCEventHandler;
class TTree;
class TChain;
class TParallelCoord;

class AliVertexerTracks;

class AliAnalysisTaskHdibaryonLPpi : public AliAnalysisTaskSE {
 
public:
  AliAnalysisTaskHdibaryonLPpi();
  AliAnalysisTaskHdibaryonLPpi(const char *name);
  //  AliAnalysisTaskHdibaryonLPpi();
  virtual ~AliAnalysisTaskHdibaryonLPpi();
  
  void SetHasMC(Bool_t hasMC=kTRUE){SetBit(kHasMC,hasMC);}
  Bool_t HasMC(){return TestBit(kHasMC);}

  //  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  //virtual void   Exec(Option_t *option);
  virtual void   UserExec(Option_t* /*option*/);

  virtual void   Terminate(Option_t *);

  /*
   *Get the MC stack 
   */
  AliMCEvent* GetMCEvent() const{return fMCEvent;}   // for CF

private:
  AliESDEvent *fESD; //AliESD event	
  //  AliESDtrackCuts* fEsdTrackCuts; // track cuts
  AliESDtrackCuts   *fESDtrackCutsV0;     // basic cut variables for v0's
  AliESDv0Cuts *fESDCutsV0; // V0 track cuts
  AliESDtrackCuts* fEsdTrackCuts; // track cuts
  
  Int_t       fBin;        //bin for the angular plot

  AliVEvent *fEvent;                //! current event

  TList       *fHistList;  //TList
  TH1F        *fHistMassDPi; // V0 online finder 
  TH1F        *fHistMassLPi; // V0 offline finder
  TH1F        *fHistMassLambdaPi; // inv mass lambda pi
  TH1F        *fHistMassLambda; // inv mass lambda from tlorentz
  TH1F        *fHistMassLambdaPPi; // inv mass Lambda p pi before cuts
  TH1F        *fHistMassLambdaP; // inv mass Lambda p 
  TH1F        *fHistMassLambdaK; // inv mass Lambda K

  TH1F        *fHistMassK0onl; // K0 mass
  TH1F        *fHistMassK0offl; // K0 mass
  TH1F        *fHistMassK0onlC; // K0 mass
  TH1F        *fHistMassK0offlC;// K0 mass
  TH1F        *fHistMassPQonl;  // K0 mass
  TH1F        *fHistMassPQoffl; // K0 mass
  TH1F        *fHistDC; // DC
  TH2F        *fHistArmenterosPodolanski; // AP distribution
  TH2F        *fHistArmenterosPodolanskiCut; // AP distributian after cuts
  TH1F        *fHistHDibaryonInvaMassGen;                //! for MC
  TH1F        *fHistHDibaryonInvaMassGenRes;                //! for MC
  TH1F        *fHistAntiHDibaryonInvaMassGen;            //! for MC
  TH1F        *fHistHDibaryonInvaMassAso;                //! for MC
  TH1F        *fHistHDibaryonInvaMassAsoReso;   //! for MC
  TH1F        *fHistAntiHDibaryonInvaMassAso;            //! for MC
  TH2F        *fHistCheck; // Online vs.  Offline V0 finder
  TH1F        *fHistHPointingAngle; // Poniting angle dist
  TH1F        *fHistMassH; // Rec. H mass
  TH1F        *fHistMassLambdaFromH; // True Lambdas from H
  TH1F        *fHistMassLambdaFromHtLorentz; // via TLorentz
  TH1F        *fHistMassPpi; //inv mass ppi
  TH1F        *fHistMassPpiReso; // inv mass ppi reso
  TH1F        *fHistMassLpi; // inv mass lambda pi
  TH1F        *fHistMassLP; // inv mass lambda p
  TH2F        *fHistProtonPIDBb;//TPC BB for Proton PID
  TH2F        *fHistPionPIDBb;//TPC BB for Pion PID
  TH2F        *fHistProtonPIDLambda;//TPC BB for Proton PID from Lambda
  TH2F        *fHistPionPIDLambda;//TPC BB for Pion PID from Lambda
  TH1F        *fHistMCdcaPvtxDvtx; //mc dca
  TH1F        *fHistMCdcaPvtxLvtx; //mc dca
  TH1F        *fHistMCdcaDvtxLvtx; // mc dca
  TH1F        *fHistMCangleLH; // angle between H and Lambda
  TH1F        *fHistMCdecayAngle; // decay angle
  TH1F        *fHistMCpointingAngle; // true PA
  TH2F        *fHistMCap; // true AP
  TH1F        *fHistMCdcaPvtxDvtxReso; //mc dca
  TH1F        *fHistMCdcaPvtxLvtxReso; // mc dca
  TH1F        *fHistMCdcaDvtxLvtxReso; // mc dca
  TH1F        *fHistMCangleLHReso; // mc anle
  TH1F        *fHistMCdecayAngleReso; // mc decay angle
  TH1F        *fHistMCpointingAngleReso; // mc PA
  TH2F        *fHistMCapReso; // MC ap for reso
  TH1F        *fHistCentrality; // Centrality histogram
  TH1F        *fHistCentralityAC; // Centrality histogram after selection cuts
  TH1F        *fHistMultiplicity; // Multiplicity histogram
 
  TH1F        *fHistHilf1; // Support
  TH1F        *fHistHilf2; // support
  TH1F        *fHistHilf3; // support
  TH1F        *fHistHilf4; // support
  TH1F        *fHistHilf5; // support
  TH1F        *fHistHilf6; // support
  TH2F        *fHistPtvsEtaGen; // eta vs pt generated
  TH2F        *fHistPtvsEtaAso; // eta vs pt associated
  TH2F        *fHistPtvsYGen; // y vs pt generated
  TH2F        *fHistPtvsYAso; // y vs pt associated
  TH1F        *fHistRap; // rapidity
  
  TH1F        *fHistCount; //Counting of events and different decay parameters
  
  AliPIDResponse *fPIDtpcESD;     //! PID response object 
  TH1F        *fHistTriggerStat;                         //! Trigger statistics
  TH1F        *fHistTriggerStatAfterEventSelection;      //! Trigger statistics
  TH3F        *fHistMassHcentMult; // Inv. mass vs. centrality vs. multiplicity
  THnSparse   *fHistNdim; //THnSparse

  enum{
        kHasMC=BIT(18)
  };

    AliAnalysisTaskHdibaryonLPpi(const AliAnalysisTaskHdibaryonLPpi&); // not implemented
    AliAnalysisTaskHdibaryonLPpi& operator=(const AliAnalysisTaskHdibaryonLPpi&); // not implemented


  ClassDef(AliAnalysisTaskHdibaryonLPpi, 1); // analysisclass


};

#endif
