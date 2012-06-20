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

//=================================================================================
// AliAnalysysTaskJetsHMPID - Class performing PID analysis in jets with the HMPID
// A set of histograms is created.
//=================================================================================

#ifndef AliAnalysisTaskJetsHMPID_H
#define AliAnalysisTaskJetsHMPID_H

#include <TList.h>
#include <TH2F.h>
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskJetsHMPID : public AliAnalysisTaskSE
{
  public:
    AliAnalysisTaskJetsHMPID();
    AliAnalysisTaskJetsHMPID(const Char_t* name);
    AliAnalysisTaskJetsHMPID& operator= (const AliAnalysisTaskJetsHMPID& c);
    AliAnalysisTaskJetsHMPID(const AliAnalysisTaskJetsHMPID& c);
    virtual ~AliAnalysisTaskJetsHMPID();

    virtual void   ConnectInputData(Option_t *);
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetJetBranch(const Char_t *br){ fJetBranch=br; }
    void SetBkgBranch(const Char_t *br){ fBkgBranch=br; }
    void SetJetPtCut(Float_t ptcut)    { fJetPtCut=ptcut; }

 protected:

 private:
  TString       fJetBranch;          // jet branch to read
  TString       fBkgBranch;          // background branch to read
  Float_t       fJetPtCut;           // jet pT threshold

  AliAODEvent   *fAOD;               // AOD object
  TList         *fHistList;          // list of histograms

  TH2F          *fThetaChJet;        // theta Cherenkov distribution in the jets
  TH2F          *fThetaChBkg;        // theta Cherenkov distribution out of the jets
  TH2F          *fThetaChRndCone;    // theta Cherenkov distribution in random cone from background
  TH1F          *fEvSelDPhi;         // Delta phi jet-HMPID in different events selected
  TH1F          *fJetsPt;            // pT of jets into the HMPID
  TH1F          *fRndConePt;         // pT of random cones into the HMPID
  TH1F          *fAwayJetPt;         // pT of jets on the away side
  TH2F          *fJetsEtaPhi;        // eta and phi of jets into the HMPID
  TH2F          *fTrksEtaPhiJet;     // eta and phi of jet tracks into the HMPID
  TH2F          *fTrksEtaPhiBkg;     // eta and phi of bkg tracks into the HMPID

  TTree         *fTree;              // tree with useful data for subsequent analysis
  Float_t        fTrackPt;           // track pt
  Float_t        fJetPt;             // jet pt
  Float_t        fPionBkg;           // pions probability out of the jets
  Float_t        fKaonBkg;           // kaons probability out of the jets
  Float_t        fProtBkg;           // prots probability out of the jets
  Float_t        fPionJet;           // pions probability in the jets
  Float_t        fKaonJet;           // kaons probability in the jets
  Float_t        fProtJet;           // prots probability in the jets

  ClassDef(AliAnalysisTaskJetsHMPID, 1);
};

#endif
