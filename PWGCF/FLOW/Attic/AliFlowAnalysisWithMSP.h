/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

#ifndef ALIFLOWANALYSISWITHMSP_H
#define ALIFLOWANALYSISWITHMSP_H

// The original implementation of the scalar produc method was done by Naomi van der Kolk and Ante Bilandzic
// This file contains a completely rewritten implementation that extends the SP method to allow for subevents 
// with different v2.
//
// Author: Paul Kuijer <mailto:Paul.Kuijer@nikhef.nl>
// 

#include <TNamed.h>

class TH1D;
class TProfile;
class TDirectoryFile;

class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowMSPHistograms;

class AliFlowAnalysisWithMSP : public TNamed
{
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// AliFlowAnalysisWithMSP                                                    //
//                                                                           //
// Use the modified scalar produc method to analyze flow                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
public:
   AliFlowAnalysisWithMSP();                             // Default constructor for root I/O
   AliFlowAnalysisWithMSP(const unsigned int harmonic, const bool commonConst=true, const bool commonHist=false);   // Construct and define params
   AliFlowAnalysisWithMSP(TDirectory *file);         // Construct and read from file
   AliFlowAnalysisWithMSP(TList *list);                  // Construct and load histograms from a TList. This object becomes owner.
   AliFlowAnalysisWithMSP(const AliFlowAnalysisWithMSP &x); // Copy constructor
   ~AliFlowAnalysisWithMSP();                            // Destructor assumes everything is owned

   void Init();                                          // (re)define all output objects, clears results
   void Make(AliFlowEventSimple *event);                 // Analyze one event
   void Finish();                                        // Calculate final results without redefining NUA usage
   void Finish(bool nua){fNUA=nua;Finish();};            // Calculate after redefining if NUA is applied
   void WriteHistograms(TDirectoryFile *file) const;     // Write all histograms to a file
   void WriteCommonResults(TDirectoryFile *file) const;  // Write an AliFlowCommonHistResults() to file
   TList *ListHistograms();                              // Create a TList from all histograms, ownership is kept in this class

   void Print(const Option_t *opt=0)const;               // Summarize results on std::cout

   void SetHarmonic(const unsigned int h){fHarmonic=h;}; // Define harmonic for the flow analysis, should be followed by Init()
   void UseCommonConstants(bool b=true){fUseCommonConstants=b;};        // Get parameters from AliFlowCommonConstants. MUST be called before Init() 
   void EnableCommonHistograms(bool c=true){fBookCommonHistograms=c;};  // Create and fill AliFlowCommonHistograms() MUST be called before Init()
   void EnableNUA(bool c=true){fNUA=c;};                 // Switch on/off NUA corrections. MUST be called before Finish()
   // TODO: Allow NUA off, NUAu integrated, NUAu per bin instead of On(per bin)/OFF. Default should be backward compatible
   // TODO: Allow setting of pt and eta binning in a simpler way, i.e. independent of AliFlowCommonConstants.

   const TH1D* GetIntegratedFlow(){return fIntegratedFlow;};   // Access result
   const TH1D* GetDiffFlowPt(){return fDiffFlowPt;};           // Access to result
   const TH1D* GetDiffFlowEta(){return fDiffFlowEta;};         // Access to result

   AliFlowAnalysisWithMSP &operator=(const AliFlowAnalysisWithMSP &x);  // Assignment 

private:
   UInt_t   fHarmonic;                       // Harmonic (n) in v_n analysis, used in Init() and Make()
   Bool_t   fNUA;                            // Correct for Non Uniform Acceptance during Finish()
   Bool_t   fUseCommonConstants;             //! Get parameters from AliFlowCommonConstants during Init()
   Bool_t   fBookCommonHistograms;           //! Book standard histograms in Init() 

   AliFlowCommonHist    *fCommonHist;        // Standard control histograms, if enabled

   // Correlations
   AliFlowMSPHistograms *fQaComponents;      // Averages of Qa components per event for NUA
   AliFlowMSPHistograms *fQbComponents;      // Averages of Qb components per event for NUA
   AliFlowMSPHistograms *fQaQb;              // Average of QaQb per event
   AliFlowMSPHistograms *fPtUComponents;     // ux and uy per pt bin for NUA
   AliFlowMSPHistograms *fEtaUComponents;    // ux and uy per eta bin for NUA
   AliFlowMSPHistograms *fAllStatistics;     // Correlations for uQa uQb and QaQb (integrated)
   AliFlowMSPHistograms *fPtStatistics;      // Correlations for uQa uQb and QaQb per pt bin
   AliFlowMSPHistograms *fEtaStatistics;     // Correlations for uQa uQb and QaQb per eta bin

   // Result histograms (if calculated)
   TH1D        *fIntegratedFlow;    // vn for POI and subevents
   TH1D        *fDiffFlowPt;        // vn as function of pt 
   TH1D        *fDiffFlowEta;       // vn as function of eta
   TProfile    *fFlags;             // Stores fHarmonic and fNUA

   TList       *fHistList;          // List of all histograms if requested by ListHistograms()

   // Helper functions 
   bool Calculate(double &vn, double &vnerror, const AliFlowMSPHistograms *hist, const AliFlowMSPHistograms *comp, const int bin=1, const int poi=0)const;  // Calculates flow from histograms with NUA corrections
   bool Calculate(double &vn, double &vnerror, const AliFlowMSPHistograms *hist, const int bin=1, const int poi=0)const{return Calculate(vn,vnerror,hist,0,bin,poi);}; // No NUA version
   void ReadHistograms(TDirectory *file);    // Restore histograms from a file
   void ReadHistograms(TList *list);         // Restore histograms from a TList

   ClassDef(AliFlowAnalysisWithMSP,0);       // Class version
};

#endif //ALIFLOWANALYSISWITHMSP_H
