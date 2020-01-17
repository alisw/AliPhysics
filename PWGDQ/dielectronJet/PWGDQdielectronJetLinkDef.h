#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#ifdef HAVE_FASTJET
// Classes which need direct access only to Fastjet objects (not
// needed if wrapped into ALICE objects)
#pragma link C++ class AliAnalysisTaskJpsiJet+;
#pragma link C++ class AliAnalysisTaskJpsiJetFilter+;
#endif
#endif
