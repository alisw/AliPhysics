/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskEx02.h 45956 2010-12-10 12:55:37Z agheata $ */
/* AliAnalysisTaskEx02.h
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#ifndef ALIANALYSISTASKEX02_H
#define ALIANALYSISTASKEX02_H

class TH1F;
class TList;
class AliESDtrackCuts;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskEx02 : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskEx02();
   AliAnalysisTaskEx02(const char *name);
   virtual ~AliAnalysisTaskEx02();

   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option);
   virtual void     UserExecMix(Option_t *);
   virtual void     Terminate(Option_t *);

   void Loop(AliESDEvent *esd);
   void LoopESDMC();
   void Loop(AliAODEvent *aod);
   void LoopAODMC();

private:
   TList           *fOutput;        // Output list
   TH1F            *fHistPt;        // Pt spectrum
   TH1F            *fHistEta;       // pseudorapidity spectrum

   TH1F            *fHistMultiDiff; // multiplicity difference
   TH1F            *fHistZVertexDiff; // multiplicity difference

   AliAnalysisTaskEx02(const AliAnalysisTaskEx02 &); // not implemented
   AliAnalysisTaskEx02 &operator=(const AliAnalysisTaskEx02 &); // not implemented

   ClassDef(AliAnalysisTaskEx02, 1); // example of analysis
};

#endif

