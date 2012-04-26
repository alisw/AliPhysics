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
class AliAnalysisManager;
class AliESDtrackCuts;
class AliMultiInputEventHandler;
class AliMixInputEventHandler;

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

   // function to loop over ESD event
   void Loop(AliESDEvent *esd);
   void LoopV0(AliESDEvent *esd);
   void LoopESDMC();

   // function to loop over AOD event
   void Loop(AliAODEvent *aod);
   void LoopV0(AliAODEvent *aod);
   void LoopAODMC();

   void SetLoopInUserExecMix(Bool_t useMix=kTRUE) { fUseLoopInUserExecMix = useMix;}
   void SetUseLoopMixedEvent(Bool_t useMixEvent=kTRUE) { fUseLoopMixedEvent = useMixEvent;}
   void SetLoopV0(Bool_t useV0=kTRUE) { fUseLoopV0 = useV0;}

   AliVEvent* GetMainEvent();
   AliVEvent* GetMixedEvent(Int_t buffId=0);

   AliMultiInputEventHandler *SetMainInputHandler(AliAnalysisManager *mgr);
   AliMixInputEventHandler *SetMixingInputHandler(AliMultiInputEventHandler *mainIH);

private:
   TList           *fOutput;               // Output list
   TH1F            *fHistPt;               // Pt spectrum
   TH1F            *fHistEta;              // pseudorapidity spectrum

   TH1F            *fHistMultiDiff;        // multiplicity difference
   TH1F            *fHistZVertexDiff;      // Vz difference

   // Following options are only for testing purposes (user can ignore)
   Bool_t           fUseLoopInUserExecMix; //
   Bool_t           fUseLoopMixedEvent; //
   Bool_t           fUseLoopV0; //
   
   AliMultiInputEventHandler *fMainInputHandler;    //! tmp pointer to mixing input handler
   AliMixInputEventHandler   *fMixingInputHandler;  //! tmp pointer to mixing input handler

   AliAnalysisTaskEx02(const AliAnalysisTaskEx02 &); // not implemented
   AliAnalysisTaskEx02 &operator=(const AliAnalysisTaskEx02 &); // not implemented

   ClassDef(AliAnalysisTaskEx02, 1); // example of analysis
};

#endif
