/************************************************************************* 
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
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

#include "AliFlowAnalysisWithMSP.h"

#include "AliFlowVector.h"
#include "AliFlowMSPHistograms.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

#include <TList.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TVector2.h>
#include <TDirectoryFile.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>

AliFlowAnalysisWithMSP::AliFlowAnalysisWithMSP() 
   : TNamed(), 
   fHarmonic(2), fNUA(kFALSE), fUseCommonConstants(kFALSE), fBookCommonHistograms(kFALSE), fCommonHist(0), fQaComponents(0), 
   fQbComponents(0), fQaQb(0), fPtUComponents(0), fEtaUComponents(0), fAllStatistics(0), 
   fPtStatistics(0), fEtaStatistics(0),
   fIntegratedFlow(0), fDiffFlowPt(0), fDiffFlowEta(0), fFlags(0), fHistList(0)
{
   // Default constructor. Intended for root IO purposes only
}     

AliFlowAnalysisWithMSP::AliFlowAnalysisWithMSP(TDirectory *file)
   : TNamed(), 
   fHarmonic(2), fNUA(kFALSE), fUseCommonConstants(kFALSE), fBookCommonHistograms(kFALSE), fCommonHist(0), fQaComponents(0), 
   fQbComponents(0), fQaQb(0), fPtUComponents(0), fEtaUComponents(0), fAllStatistics(0), 
   fPtStatistics(0), fEtaStatistics(0),
   fIntegratedFlow(0), fDiffFlowPt(0), fDiffFlowEta(0), fFlags(0), fHistList(0)
{
   // Constructor reads the internal state from a root file.
   // No check for consistency is done. If flags or histograms are not found they are left at default (0).
   ReadHistograms(file);
}

AliFlowAnalysisWithMSP::AliFlowAnalysisWithMSP(TList *list)
   : TNamed(), 
   fHarmonic(2), fNUA(kFALSE), fUseCommonConstants(kFALSE), fBookCommonHistograms(kFALSE), fCommonHist(0), fQaComponents(0), 
   fQbComponents(0), fQaQb(0), fPtUComponents(0), fEtaUComponents(0), fAllStatistics(0), 
   fPtStatistics(0), fEtaStatistics(0),
   fIntegratedFlow(0), fDiffFlowPt(0), fDiffFlowEta(0), fFlags(0), fHistList(0)
{
   // Constructor reads the internal state from a root file.
   // No check for consistency is done. If flags or histograms are not found they are left at default (0).
   ReadHistograms(list);
}


AliFlowAnalysisWithMSP::AliFlowAnalysisWithMSP(const unsigned int harmonic, const bool commonConst, const bool commonHist) 
   : TNamed(), 
   fHarmonic(harmonic), fNUA(kFALSE), fUseCommonConstants(commonConst), fBookCommonHistograms(commonHist), fCommonHist(0), fQaComponents(0), 
   fQbComponents(0), fQaQb(0), fPtUComponents(0), fEtaUComponents(0), fAllStatistics(0), 
   fPtStatistics(0), fEtaStatistics(0),
   fIntegratedFlow(0), fDiffFlowPt(0), fDiffFlowEta(0), fFlags(0), fHistList(0)
{
   // Constructor defining the harmonic, usage of non uniform acceptance corrections and the AliFlowCommonHist() histograms
   // Calls Init() to set-up the histograms. This is equivalent to:
   // AliFlowAnalysisWithMSP *ana= new AliFlowAnalysisWithMSP);
   // ana->SetHarmonic(harmonic);
   // ana->UseCommonConstants(commonConst);
   // ana->EnableCommonHistograms(commonHist);
   // ana->Init();
   // This is the constructor intended for the user
   SetNameTitle("MSP","Flow analysis with the Modified Scalar Product method");
   Init();
}

AliFlowAnalysisWithMSP::AliFlowAnalysisWithMSP(const AliFlowAnalysisWithMSP &x)
   : TNamed(),
   fHarmonic(x.fHarmonic), fNUA(x.fNUA), fUseCommonConstants(x.fUseCommonConstants), fBookCommonHistograms(x.fBookCommonHistograms), fCommonHist(0), 
   fQaComponents(0), fQbComponents(0), fQaQb(0), fPtUComponents(0), fEtaUComponents(0), fAllStatistics(0), 
   fPtStatistics(0), fEtaStatistics(0),
   fIntegratedFlow(0), fDiffFlowPt(0), fDiffFlowEta(0), fFlags(0), fHistList(0)
{
   // Copy constructor
   SetNameTitle("MSP","Flow analysis with the Modified Scalar Product method");
   if( x.fQaComponents )   fQaComponents=(AliFlowMSPHistograms *)(x.fQaComponents)->Clone();
   if( x.fQbComponents )   fQbComponents=(AliFlowMSPHistograms *)(x.fQbComponents)->Clone();
   if( x.fQaQb )           fQaQb=(AliFlowMSPHistograms *)(x.fQaQb)->Clone();
   if( x.fPtUComponents)     fPtUComponents=(AliFlowMSPHistograms *)(x.fPtUComponents)->Clone();
   if( x.fEtaUComponents )   fEtaUComponents=(AliFlowMSPHistograms *)(x.fEtaUComponents)->Clone();
   if( x.fAllStatistics )    fAllStatistics=(AliFlowMSPHistograms *)(x.fAllStatistics)->Clone();
   if( x.fPtStatistics )     fPtStatistics=(AliFlowMSPHistograms *)(x.fPtStatistics)->Clone();
   if( x.fEtaStatistics )    fEtaStatistics=(AliFlowMSPHistograms *)(x.fEtaStatistics)->Clone();
   if( x.fIntegratedFlow )   fIntegratedFlow=(TH1D *)(x.fIntegratedFlow)->Clone();
   if( x.fDiffFlowPt )       fDiffFlowPt=(TH1D *)(x.fDiffFlowPt)->Clone();
   if( x.fDiffFlowEta )      fDiffFlowEta=(TH1D *)(x.fDiffFlowEta)->Clone();
   if( x.fFlags )            fFlags=(TProfile *)(x.fFlags)->Clone();
   if( x.fCommonHist )       fCommonHist=new AliFlowCommonHist(*(x.fCommonHist));
   // The histogram list fHistList is not cloned because it is regenerated when requested
}

AliFlowAnalysisWithMSP::~AliFlowAnalysisWithMSP() 
{
   // Destructor. All internal objects are owned by this class and deleted here.
   delete fCommonHist;
   delete fQaComponents;
   delete fQbComponents;
   delete fQaQb;
   delete fPtUComponents;
   delete fEtaUComponents;
   delete fAllStatistics;
   delete fPtStatistics;
   delete fEtaStatistics;
   delete fIntegratedFlow;
   delete fDiffFlowPt;
   delete fDiffFlowEta;
   delete fFlags;
   if(fHistList) {
      fHistList->SetOwner(kFALSE);   // Histograms were already deleted manually
      delete fHistList;              // Delete the TList itself
   }
}

void AliFlowAnalysisWithMSP::Init()
{
   // Create all output objects. Memory consumption can be reduced by switching off some of the control histograms.
   // pt and eta binning are set to the values defined in the static class AliFlowCommonConstants if enabled. Otherwise defaults are set.
   delete fCommonHist;     fCommonHist=0;    // Delete existing histograms
   delete fQaComponents;   fQaComponents=0;
   delete fQbComponents;   fQbComponents=0;
   delete fQaQb;           fQaQb=0;
   delete fPtUComponents;  fPtUComponents=0;
   delete fEtaUComponents; fEtaUComponents=0;
   delete fAllStatistics;  fAllStatistics=0;
   delete fPtStatistics;   fPtStatistics=0;
   delete fEtaStatistics;  fEtaStatistics=0;
   delete fIntegratedFlow; fIntegratedFlow=0;
   delete fDiffFlowPt;     fDiffFlowPt=0;
   delete fDiffFlowEta;    fDiffFlowEta=0;
   delete fFlags;          fFlags=0;
   if(fHistList) {
      fHistList->SetOwner(kFALSE);   // Histograms were already deleted manually
      delete fHistList;              // Delete the TList itself
      fHistList=0;                   // Clear pointer which is invalid afcter delete
   }

   // Default binning for histograms
   // TODO: allow variable binning in a simpler way than by AliFlowCommonConstants() only

   // Defaults
   int nBinsPt=100;
   double ptMin=0;
   double ptMax=10.0;
   int nBinsEta=100;
   double etaMin=-0.8;
   double etaMax=+0.8;

   if( fUseCommonConstants || fBookCommonHistograms ) {
      AliFlowCommonConstants *c=AliFlowCommonConstants::GetMaster();
      nBinsPt=c->GetNbinsPt();
      ptMin=c->GetPtMin();
      ptMax=c->GetPtMax();
      nBinsEta=c->GetNbinsEta();
      etaMin=c->GetEtaMin();
      etaMax=c->GetEtaMax();
   }

   if( fBookCommonHistograms ) {           // Use the common constants from the flow package for histogram definitions
      std::cerr << "AliFlowAnalysisWithMSP::Init() Creating common histograms" << std::endl;
      fCommonHist=new AliFlowCommonHist("AliFlowCommonHist_MSP","AliFlowCommonHist",kTRUE);
   }

   std::cerr << "AliFlowAnalysisWithMSP::Init() creating MSPHistograms" << std::endl;
   // Averages for NUA corrections:
   fQaComponents=new AliFlowMSPHistograms(2,"QaComponents",1,0.5,1.5);           // QaX, QaY
   fQaComponents->SetXName("all");
   fQaComponents->SetVarName("QaX",0);
   fQaComponents->SetVarName("QaY",1);
   fQbComponents=new AliFlowMSPHistograms(2,"QbComponents",1,0.5,1.5);           // QbX, QbY
   fQbComponents->SetXName("all");
   fQbComponents->SetVarName("QbX",0);
   fQbComponents->SetVarName("QbY",1);
   fQaQb=new AliFlowMSPHistograms(1,"QaQb",1,0.5,1.5);                           // QaQb
   fQaQb->SetXName("all");
   fQaQb->SetVarName("QaQb",0);
   fPtUComponents=new AliFlowMSPHistograms(2,"PtUComponents",nBinsPt,ptMin,ptMax);      // ux(pt), uy(pt)
   fPtUComponents->SetXName("pt");
   fPtUComponents->SetVarName("ux",0);
   fPtUComponents->SetVarName("uy",1);
   fEtaUComponents=new AliFlowMSPHistograms(2,"EtaUComponents",nBinsEta,etaMin,etaMax);    // ux(eta), uy(eta)
   fEtaUComponents->SetXName("eta");
   fEtaUComponents->SetVarName("ux",0);
   fEtaUComponents->SetVarName("uy",1);

   //  Correlation terms
   fAllStatistics=new AliFlowMSPHistograms(3,"AllStatistics",1,0.5,1.5);          // terms integrated over pt and eta
   fAllStatistics->SetXName("all");
   fAllStatistics->SetVarName("uQa",0);
   fAllStatistics->SetVarName("uQb",1);
   fAllStatistics->SetVarName("QaQb",2);
   fPtStatistics=new AliFlowMSPHistograms(3,"PtStatistics",nBinsPt,ptMin,ptMax);         // terms per pt bin
   fPtStatistics->SetXName("pt");
   fPtStatistics->SetVarName("uQa",0);
   fPtStatistics->SetVarName("uQb",1);
   fPtStatistics->SetVarName("QaQb",2);
   fEtaStatistics=new AliFlowMSPHistograms(3,"EtaStatistics",nBinsEta,etaMin,etaMax);      // terms per eta bin
   fEtaStatistics->SetXName("eta");
   fEtaStatistics->SetVarName("uQa",0);
   fEtaStatistics->SetVarName("uQb",1);
   fEtaStatistics->SetVarName("QaQb",2);

   fIntegratedFlow=0;   // Created in Finish()
   fDiffFlowPt=0;       // Created in Finish
   fDiffFlowEta=0;      // Created in Finish

   fFlags=new TProfile("Flags","Flags for AliFlowAnalysisWithMSP",10,0.5,10.5,"s");
   fFlags->Fill("Harmonic",fHarmonic); // bin 1
}

void AliFlowAnalysisWithMSP::Make(AliFlowEventSimple *event)
{
   // Analyze one event. The modified scalar product method estimates flow using the formula:
   // BEGIN_LATEX v_2(MSP) = \sqrt( uQ^{a} uQ^{b} / (Q^{a}Q^{b}) ) END_LATEX
   // The Q vectors are calculated for the harmonic set previously: fHarmonic.
   // Make(event) does not use the new operator, thus avoiding memory leaks in a simple way
   // Depending on the compiler about 200 bytes may be pushed/popped on the stack for local variables. 
   // 

   if( !event ) return;          // Protect against running without events. Can this happen???

   AliFlowVector flowVectors[2];                   // Calculate the two subevent Q vectors:
   event->Get2Qsub(flowVectors, fHarmonic);        // TODO: Check how do the phi, pt, eta weights work in the flow event?

   AliFlowVector &Qa=flowVectors[0];               // Define some mnemonics for the subevent flow vectors:
   AliFlowVector &Qb=flowVectors[1];

   const double QaW=Qa.GetMult();      // Weight for Qa and combinations of Qa
   const double QbW=Qb.GetMult();      // Weight for Qb and combinations of Qb

   if( QaW<2 || QbW<2 ) return;  // Require at least 2 particles in each subevent

   if(fCommonHist)fCommonHist->FillControlHistograms(event);   // Standard for all flow analysis

   // Fill NUA correction histograms for Qa, Qb and QaQb per event:
   const double qaxy[]={Qa.X()/QaW,Qa.Y()/QaW};                // Two variables: Qa components
   const double wqaxy[]={QaW,QaW};
   fQaComponents->Fill(1, qaxy, wqaxy);                        // only one bin (all pt)

   const double qbxy[]={Qb.X()/QbW,Qb.Y()/QbW};                // Two variables: Qb components
   const double wqbxy[]={QbW,QbW};
   fQbComponents->Fill(1, qbxy, wqbxy);                        // only one bin (all pt)

   const double QaQbW[]={QaW*QbW};
   const double weightedQaQb[] = {(Qa*Qb)/QaQbW[0]};                  // Scalar product of subevent Q vectors with weight, only 1 variable
   fQaQb->Fill(1,weightedQaQb,QaQbW);                        // Average of QaQb per event
   
   int iTrack=0;
   while( AliFlowTrackSimple *track=event->GetTrack(iTrack++) ) {// Loop over the tracks in the event 
      if(! track->InPOISelection() ) continue;  // Ignore if not a POI
      const double trackWeight=track->Weight(); // Get the track vector
      const double phi=track->Phi();
      const double pt=track->Pt();
      const double eta=track->Eta();

      AliFlowVector u;
      u.SetMagPhi(trackWeight, fHarmonic*phi, trackWeight);    // Length = track weight = multiplicity for a single track

      // Remove track from subevent a 
      AliFlowVector mQa(Qa);                 // Initialize with Qa flow vector of this event
      if( track->InSubevent(0) ) {           // Correct for autocorrelation with POI for this track
         mQa-=u;                             
      }

      // Remove track from subevent b 
      AliFlowVector mQb(Qb);                 // Initialize with Qb flow vector of this event
      if( track->InSubevent(1) ) {           // Correct for autocorrelation with POI for this track
         mQb-=u;                             
      }

      const double uQaW = mQa.GetMult()*trackWeight;     // Weight is multiplicity of Q vector times weight of POI
      const double uQbW = mQb.GetMult()*trackWeight;
      const double uQa=u*mQa;                            // Correlate POI with subevent 
      const double uQb=u*mQb;                

      const double uxy[]={u.X(),u.Y()};
      const double wxy[]={1.0,1.0};
      fPtUComponents->Fill(pt, uxy, wxy);    // NUA for POI vs pt
      fEtaUComponents->Fill(eta, uxy, wxy);  // NUA for POI vs eta

      const double par[]={uQa/uQaW, uQb/uQbW, weightedQaQb[0]};   // Three variables to correlate: uQa, uQb and QaQb
      const double wgt[]={uQaW, uQbW, QaQbW[0]};
      fAllStatistics->Fill(1,   par, wgt  );    // only 1 bin, integrated over pt and eta
      fPtStatistics->Fill(pt,   par, wgt );     // pt differential correlations
      fEtaStatistics->Fill(eta, par, wgt );     // eta differential correlations
   }
}

void AliFlowAnalysisWithMSP::Finish()
{
   // Calculate the final result from the stored correlations.
   // The NUA corrections are applied if the flag fNUA was set before the call to Finish()
   // If the output histograms already exist then they are replaced by the newly calculated result

   Print();       // Print a summary of the NUA terms and integrated flow

   // Create result histograms
   fFlags->Fill("NUA",fNUA);
   delete fIntegratedFlow; fIntegratedFlow=0;   // First delete existing results (if any)
   delete fDiffFlowPt;     fDiffFlowPt=0;
   delete fDiffFlowEta;    fDiffFlowEta=0;

   bool oldAddStatus=TH1::AddDirectoryStatus();       // Do not store the next histograms automatically
   TH1::AddDirectory(kFALSE);                         // We need full control over the writing of the hostograms

   fIntegratedFlow=new TH1D("IntegratedFlow","Integrated flow results",10,0.5,10.5);
   fIntegratedFlow->SetDirectory(0);
   double vn, vnerror;
   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 0);    
   fIntegratedFlow->SetBinContent(1,vn);
   fIntegratedFlow->SetBinError(1,vnerror);
   fIntegratedFlow->GetXaxis()->SetBinLabel(1,"POI");
   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 1);
   fIntegratedFlow->SetBinContent(2,vn);
   fIntegratedFlow->SetBinError(2,vnerror);
   fIntegratedFlow->GetXaxis()->SetBinLabel(2,"A");
   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 2);
   fIntegratedFlow->SetBinContent(3,vn);
   fIntegratedFlow->SetBinError(3,vnerror);
   fIntegratedFlow->GetXaxis()->SetBinLabel(3,"B");

   int nbinspt=fPtStatistics->Nbins();
   double ptlow=fPtStatistics->XLow();
   double pthigh=fPtStatistics->XHigh();
   fDiffFlowPt = new TH1D("DiffFlowPt","flow of POI vs pt",nbinspt,ptlow,pthigh);
   fDiffFlowPt->SetDirectory(0);          // Do not automatically store in file
   fDiffFlowPt->SetStats(kFALSE);
   fDiffFlowPt->SetXTitle("p_{t}");
   fDiffFlowPt->SetYTitle(Form("v_{%d}",fHarmonic));

   for(int i=1; i<=nbinspt; ++i) {     
      double vnpt=0;
      double evnpt=0;
      if( Calculate(vnpt, evnpt, fPtStatistics, fPtUComponents, i) && evnpt<1 ) {
         fDiffFlowPt->SetBinContent(i,vnpt);
         fDiffFlowPt->SetBinError(i,evnpt);
      }
   }

   int nbinseta=fEtaStatistics->Nbins();
   double etalow=fEtaStatistics->XLow();
   double etahigh=fEtaStatistics->XHigh();
   fDiffFlowEta= new TH1D("DiffFlowEta","flow of POI vs #eta",nbinseta,etalow,etahigh);
   fDiffFlowEta->SetDirectory(0);          // Do not automatically store in file
   fDiffFlowEta->SetStats(kFALSE);
   fDiffFlowEta->SetXTitle("#eta");
   fDiffFlowEta->SetYTitle(Form("v_{%d}",fHarmonic));

   for(int i=1; i<=nbinseta; ++i) {     
      double vneta=0;
      double evneta=0;
      if( Calculate(vneta, evneta, fEtaStatistics, fEtaUComponents, i) && evneta<1 ) {
         fDiffFlowEta->SetBinContent(i,vneta);
         fDiffFlowEta->SetBinError(i,evneta);
      }
   }

   TH1::AddDirectory(oldAddStatus);
}

void AliFlowAnalysisWithMSP::WriteHistograms(TDirectoryFile *file) const
{
   // Write the internal status to file. If the AliFlowCommonHistograms were enabled they are written also.
   // Results are written only if they were calculated by Finish().
   // From these histograms the internal state of *this can be reconstructed with ReadHistograms()

   //file->Write(file->GetName(), TObject::kSingleKey);     // Make sure the directory itself is written

   if(fCommonHist) file->WriteTObject(fCommonHist);

   TList *t=new TList();                                    // This object is written as a marker for redoFinish()
   t->SetName("cobjMSP");
   file->WriteTObject(t,0,"Overwrite");
   file->WriteTObject(fQaComponents,0,"Overwrite");         // Averages of Qa components per event
   file->WriteTObject(fQbComponents,0,"Overwrite");         // Averages of Qb components per event
   file->WriteTObject(fQaQb,0,"Overwrite");                 // Average of QaQb per event
   file->WriteTObject(fPtUComponents,0,"Overwrite");        // u components vs pt
   file->WriteTObject(fEtaUComponents,0,"Overwrite");       // u components vs eta
   file->WriteTObject(fAllStatistics,0,"Overwrite");        // Integrated uQa, uQb and QaQa
   file->WriteTObject(fPtStatistics,0,"Overwrite");         // uQa, uQb and QaQb vs pt
   file->WriteTObject(fEtaStatistics,0,"Overwrite");        // uQa, uQb and QaQb vs eta

   if( fIntegratedFlow ) file->WriteTObject(fIntegratedFlow,0,"Overwrite");  // Integrated flow for POI and subevents
   if( fDiffFlowPt ) file->WriteTObject(fDiffFlowPt,0,"Overwrite");          // Differential flow vs pt if calculated
   if( fDiffFlowEta ) file->WriteTObject(fDiffFlowEta,0,"Overwrite");        // Differential flow vs eta if calculated

   file->WriteTObject(fFlags,0,"Overwrite");               // fHarmonic, fNUA

   file->WriteKeys();   // Make sure it happens now
}

void AliFlowAnalysisWithMSP::WriteCommonResults(TDirectoryFile *file) const
{
   // Export the results to a AliFlowCommonHistResults() class and write to file
   // If the results were not calculated then Finish() is called to generate them
   
   int nBinsPt=fPtStatistics->Nbins();       // Get the actually used binning from the Statistics histograms
   double ptMin=fPtStatistics->XLow();
   double ptMax=fPtStatistics->XHigh();

   int nBinsEta=fEtaStatistics->Nbins();
   double etaMin=fEtaStatistics->XLow();
   double etaMax=fEtaStatistics->XHigh();

   AliFlowCommonConstants *c=AliFlowCommonConstants::GetMaster(); // Get the static common constants object
   // Save the old AliFlowCommonConstants status
   int oldNbinsPt=c->GetNbinsPt();     
   double oldPtMin=c->GetPtMin();      
   double oldPtMax=c->GetPtMax();      
   int oldNbinsEta=c->GetNbinsEta();   
   double oldEtaMin=c->GetEtaMin();    
   double oldEtaMax=c->GetEtaMax(); 
   // Modify AliFlowCommonConstants to make sure that AliFlowCommonResults is generated with the correct binning
   c->SetNbinsPt(nBinsPt);
   c->SetPtMin(ptMin);
   c->SetPtMax(ptMax);
   c->SetNbinsEta(nBinsEta);
   c->SetEtaMin(etaMin);
   c->SetEtaMax(etaMax);

   bool oldAddStatus=TH1::AddDirectoryStatus();       // Do not store the next histograms automatically
   TH1::AddDirectory(kFALSE);                         // We need full control over the writing of the hostograms
   AliFlowCommonHistResults *h=new AliFlowCommonHistResults("AliFlowCommonHistResults_MSP","AliFlowCommonHistResults from the MSP method",fHarmonic);

   double ivn, ivnerror;
   Calculate(ivn, ivnerror, fAllStatistics, fPtUComponents, 0, 0);    
   h->FillIntegratedFlowPOI(ivn, ivnerror);

   for(int bin=1; bin<=nBinsPt; ++bin) {
      double vn=0;
      double evn=0;
      if( Calculate(vn, evn, fPtStatistics, fPtUComponents, bin) && evn>0 ) {
         h->FillDifferentialFlowPtPOI(bin, vn, evn);
      }
   }

   for(int bin=1; bin<=nBinsEta; ++bin) {
      double vn=0;
      double evn=0;
      if( Calculate(vn, evn, fEtaStatistics, fEtaUComponents, bin) && evn>0 ) {
         h->FillDifferentialFlowEtaPOI(bin, vn, evn);
      }
   }

   file->WriteTObject(h,0,"Overwrite");

   TH1::AddDirectory(oldAddStatus);    // Restore the automatic storage of histograms to its original status

   // Restore AliFlowCommonConstants to make sure that no other analysis are affected
   c->SetNbinsPt(oldNbinsPt);
   c->SetPtMin(oldPtMin);
   c->SetPtMax(oldPtMax);
   c->SetNbinsEta(oldNbinsEta);
   c->SetEtaMin(oldEtaMin);
   c->SetEtaMax(oldEtaMax);
}

TList *AliFlowAnalysisWithMSP::ListHistograms()
{
   if( fHistList ) {
      fHistList->SetOwner(kFALSE);
      delete fHistList;
   }
   fHistList = new TList();
   fHistList->SetOwner(kFALSE);

   if(fCommonHist)      fHistList->Add(fCommonHist);        // Standard control histograms, if enabled

   //Correlations
   if(fQaComponents)    fHistList->Add(fQaComponents);      // Averages of Qa components per event for NUA
   if(fQbComponents)    fHistList->Add(fQbComponents);      // Averages of Qb components per event for NUA
   if(fQaQb)            fHistList->Add(fQaQb);              // Average of QaQb per event
   if(fPtUComponents)   fHistList->Add(fPtUComponents);     // ux and uy per pt bin for NUA
   if(fEtaUComponents)  fHistList->Add(fEtaUComponents);    // ux and uy per eta bin for NUA
   if(fAllStatistics)   fHistList->Add(fAllStatistics);     // Correlations for uQa uQb and QaQb (integrated)
   if(fPtStatistics)    fHistList->Add(fPtStatistics);      // Correlations for uQa uQb and QaQb per pt bin
   if(fEtaStatistics)   fHistList->Add(fEtaStatistics);     // Correlations for uQa uQb and QaQb per eta bin

   // Result histograms (if calculated)
   if(fIntegratedFlow)  fHistList->Add(fIntegratedFlow);    // vn for POI and subevents
   if(fDiffFlowPt)      fHistList->Add(fDiffFlowPt);        // vn as function of pt 
   if(fDiffFlowEta)     fHistList->Add(fDiffFlowEta);       // vn as function of eta
   if(fFlags)           fHistList->Add(fFlags);             // Stores fHarmonic and fNUA

   return fHistList;
}


void AliFlowAnalysisWithMSP::Print(const Option_t *opt)const
{
   if( opt ) std::cout << std::endl;

   std::cout << "****************************************************" << std::endl;
   std::cout << "    Integrated flow from Modified Scalar Product    " << std::endl;
   std::cout << "                                                    " << std::endl;

   double vn=0;
   double vnerror=0;

   std::cout << setprecision(4);
   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 0);      
   std::cout << "v" << fHarmonic << " for POI       : " << setw(11) << vn << " +- " << setw(9) << vnerror << std::endl;
   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 1);
   std::cout << "v" << fHarmonic << " for subevent A: " << setw(11) << vn << " +- " << setw(9) << vnerror << std::endl;
   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 2);
   std::cout << "v" << fHarmonic << " for subevent B: " << setw(11) << vn << " +- " << setw(9) << vnerror << std::endl;
   std::cout << std::endl;

   std::cout << "NUA terms: " << (fNUA?"(applied)":"(NOT applied)") << std::endl;
   std::cout << setprecision(3);
   const double ux=fPtUComponents->Average(0);       // Average over all bins
   const double eux=TMath::Sqrt(fPtUComponents->Variance(0));
   std::cout << "<ux>       " << setw(12) << ux << " +- " << setw(12) << eux << (TMath::Abs(ux)<2*eux?" NOT significant ":" ") << std::endl;
   const double ux0eta=fEtaUComponents->Average(fEtaUComponents->FindBin(0.0),0);
   const double eux0eta=TMath::Sqrt(fEtaUComponents->Variance(fEtaUComponents->FindBin(0.0),0));
   std::cout << "<ux> eta=0 " << setw(12) << ux0eta << " +- " <<  setw(12) << eux0eta << (TMath::Abs(ux0eta)<2*eux0eta?" NOT significant":" ") << std::endl;
   const double ux0pt=fPtUComponents->Average(fPtUComponents->FindBin(1.0),0);
   const double eux0pt=TMath::Sqrt(fPtUComponents->Variance(fPtUComponents->FindBin(1.0),0));
   std::cout << "<ux> pt=1  " << setw(12) << ux0pt << " +- " <<  setw(12) << eux0pt << (TMath::Abs(ux0pt)<2*eux0pt?" NOT significant":" ") << std::endl;

   const double uy=fPtUComponents->Average(0);        // Average over all bins 
   const double euy=TMath::Sqrt(fPtUComponents->Variance(0));
   std::cout << "<uy>       " << setw(12) << uy << " +- " << setw(12) << euy << (TMath::Abs(uy)<2*euy?" NOT significant ":" ") << std::endl;;
   const double uy0eta=fEtaUComponents->Average(fEtaUComponents->FindBin(0.0),1);
   const double euy0eta=TMath::Sqrt(fEtaUComponents->Variance(fEtaUComponents->FindBin(0.0),1));
   std::cout << "<uy> eta=0 " << setw(12) << uy0eta << " +- " << setw(12) << euy0eta << (TMath::Abs(uy0eta)<2*euy0eta?" NOT significant ":" ") << std::endl;
   const double uy0pt=fPtUComponents->Average(fPtUComponents->FindBin(1.0),1);
   const double euy0pt=TMath::Sqrt(fPtUComponents->Variance(fPtUComponents->FindBin(1.0),1));
   std::cout << "<uy> pt=1  " << setw(12) << uy0pt << " +- " << setw(12) << euy0pt << (TMath::Abs(uy0pt)<2*euy0pt?" NOT significant ":" ") << std::endl;

   const double ax=fQaComponents->Average(0);
   const double eax=TMath::Sqrt(fQaComponents->Variance(0));
   std::cout << "<QaX>      " << setw(12) << ax << " +- " << setw(12) <<eax << (TMath::Abs(ax)<2*eax?" NOT significant ":" ") << std::endl;
   const double ay=fQaComponents->Average(1);
   const double eay=TMath::Sqrt(fQaComponents->Variance(1));
   std::cout << "<QaY>      " << setw(12) << ay << " +- " << setw(12) << eay << (TMath::Abs(ay)<2*eay?" NOT significant ":" ") << std::endl;
   const double bx=fQbComponents->Average(0);
   const double ebx=TMath::Sqrt(fQbComponents->Variance(0));
   std::cout << "<QbX>      " << setw(12) << bx << " +- " << setw(12) << ebx << (TMath::Abs(bx)<2*ebx?" NOT significant ":" ") << std::endl;
   const double by=fQbComponents->Average(1);
   const double eby=TMath::Sqrt(fQbComponents->Variance(1));
   std::cout << "<QbY>      " << setw(12) << by << " +- " << setw(12) << eby << (TMath::Abs(by)<2*eby?" NOT significant ":" ") << std::endl;
   const double ab=fQaQb->Average(0);
   const double eab=TMath::Sqrt(fQbComponents->Variance(0));
   std::cout << "<QaQb>     " << setw(12) << ab << " +- " << setw(12) << eab << (TMath::Abs(ab)<2*eab?" NOT significant ":" ") << std::endl;
   std::cout << std::endl;

   std::cout << "Covariance matrix: " << std::endl;
   std::cout << "     " << setw(12) << "uQa"                           << setw(12) << "uQb"                           << setw(12) << "QaQb" << std::endl;
   std::cout << "uQa  " << setw(12) << fAllStatistics->Covariance(0,0) << std::endl;
   std::cout << "uQb  " << setw(12) << fAllStatistics->Covariance(1,0) << setw(12) << fAllStatistics->Covariance(1,1) << std::endl;
   std::cout << "QaQb " << setw(12) << fAllStatistics->Covariance(2,0) << setw(12) << fAllStatistics->Covariance(2,1) << setw(12) << fQaQb->Variance(0) << std::endl;
   std::cout << "****************************************************" << std::endl;
   std::cout << std::endl;
}


AliFlowAnalysisWithMSP &AliFlowAnalysisWithMSP::operator=(const AliFlowAnalysisWithMSP &x)
{
   if (&x==this) return *this; //handle self assignmnet
   SetNameTitle("MSP","Flow analysis with the Modified Scalar Product method");
   delete fQaComponents; fQaComponents=0;
   if( x.fQaComponents )   fQaComponents=(AliFlowMSPHistograms *)(x.fQaComponents)->Clone();
   delete fQbComponents; fQbComponents=0;
   if( x.fQbComponents )   fQbComponents=(AliFlowMSPHistograms *)(x.fQbComponents)->Clone();
   delete fQaQb; fQaQb=0;
   if( x.fQaQb )           fQaQb=(AliFlowMSPHistograms *)(x.fQaQb)->Clone();
   delete fPtUComponents; fPtUComponents=0;
   if( fPtUComponents)     fPtUComponents=(AliFlowMSPHistograms *)(x.fPtUComponents)->Clone();
   delete fEtaUComponents; fEtaUComponents=0;
   if( fEtaUComponents )   fEtaUComponents=(AliFlowMSPHistograms *)(x.fEtaUComponents)->Clone();
   delete fAllStatistics; fAllStatistics=0;
   if( fAllStatistics )    fAllStatistics=(AliFlowMSPHistograms *)(x.fAllStatistics)->Clone();
   delete fPtStatistics; fPtStatistics=0;
   if( fPtStatistics )     fPtStatistics=(AliFlowMSPHistograms *)(x.fPtStatistics)->Clone();
   delete fEtaStatistics; fEtaStatistics=0;
   if( fEtaStatistics )    fEtaStatistics=(AliFlowMSPHistograms *)(x.fEtaStatistics)->Clone();
   delete fIntegratedFlow; fIntegratedFlow=0;
   if( fIntegratedFlow )   fIntegratedFlow=(TH1D *)(x.fIntegratedFlow)->Clone();
   delete fDiffFlowPt; fDiffFlowPt=0;
   if( fDiffFlowPt )       fDiffFlowPt=(TH1D *)(x.fDiffFlowPt)->Clone();
   delete fDiffFlowEta; fDiffFlowEta=0;
   if( fDiffFlowEta )      fDiffFlowEta=(TH1D *)(x.fDiffFlowEta)->Clone();
   delete fFlags; fFlags=0;
   if( fFlags )            fFlags=(TProfile *)(x.fFlags)->Clone();
   delete fCommonHist; fCommonHist=0;
   if( fCommonHist ) fCommonHist=new AliFlowCommonHist(*(x.fCommonHist));
   return *this;
}

// private functions --------------------------------------------------------------------------------------
bool AliFlowAnalysisWithMSP::Calculate(double &vn, double &vnerror, const AliFlowMSPHistograms *hist, const AliFlowMSPHistograms *components, const int bin, const int poi) const
{
   // Get all averages and correlations need for the flow calculation
   double uQa=hist->Average(bin,0);                // <<uQa>>
   double VuQa=hist->Variance(bin,0);              // Var(<<uQa>>)
   double uQb=hist->Average(bin,1);                // <<uQb>>
   double VuQb=hist->Variance(bin,1);              // Var(<<uQb>>)
   double QaQb=fQaQb->Average(1,0);                // <QaQb> Should not be taken from hist(bin) because there QaQb is entered multiple times: <<QaQb>>!!
   double VQaQb=fQaQb->Variance(1,0);              // V(<QaQb>)
   double CuQauQb=hist->Covariance(bin,0,1);       // Cov(<<uQa>>,<<uQb>>)
   double CuQaQaQb=hist->Covariance(bin,0,2);      // Cov(<<uQa>>,<QaQb>)
   double CuQbQaQb=hist->Covariance(bin,1,2);      // Cov(<<uQb>>,<QaQb>)

   if( fNUA && components ) {      
      // Apply NUA correction to QaQb.
      double QaX=fQaComponents->Average(0);
      double QaY=fQaComponents->Average(1);
      double QbX=fQbComponents->Average(0);
      double QbY=fQbComponents->Average(1);

      QaQb=QaQb-QaX*QbX-QaY*QbY;

      // Apply NUA to uQa and uQb (per bin)
      double uX=components->Average(bin,0);  // bin 0 is integrated over all bins
      double uY=components->Average(bin,1);

      uQa = uQa - uX*QaX - uY*QaY;
      uQb = uQb - uX*QbX - uY*QbY;
      // Error calculation not fully modified: only above terms, the spread in <<u>> and <Qa> and <Qb> are not used
      // therefore this should only be applied if significant!
      // Check if not fully NUA correcting the error calculation is justified (but this is compatible with the original SP method)!
      // In general it is not justified but this can only be checked by splitting the event sample in many subsamples and looking at
      // the variance of the result. This may not be feasible if statistics is low
   }

   // Some sanity checks:
   if( uQa*uQb*QaQb <= 0 ) {        // Catch imaginary results
      vn=-99;
      vnerror=-99;
      return false;
   }

   // Sanity checks passed, calculate, print and store
   switch (poi) {
   case 1:                                        // Subevent A reference flow
      {
         if( TMath::Abs(uQb) < 1e-30*TMath::Abs(uQa*QaQb) ) {  // Protect against infinity
            vn=0;
            vnerror=-1;
            return false;
         }
         double vnA = TMath::Sqrt( uQa*QaQb / uQb ); // vn
         double VvnA = QaQb*VuQa/(4*uQa*uQb)         // Variance of vn
            + uQa*VQaQb/(4*QaQb*uQb) 
            + uQa*QaQb*VuQb/(4*TMath::Power(uQb,3)) 
            + CuQaQaQb/(2*uQb)
            - QaQb*CuQauQb/(2*TMath::Power(uQb,2))
            - uQa*CuQbQaQb/(2*TMath::Power(uQb,2));
         vn=vnA;
         if( VvnA<0 ) {
            vnerror=VvnA;
            return false;
         }
         vnerror=TMath::Sqrt(VvnA);
      }
      break;
   case 2:                                        // Subevent B reference flow
      {
         if( TMath::Abs(uQa) < 1e-30*TMath::Abs(uQb*QaQb) ) {  // Protect against infinity
            vn=0;
            vnerror=-1;
            return false;
         }
         double vnB = TMath::Sqrt( uQb*QaQb / uQa );        // vn
         double VvnB = uQb*VQaQb/(4*QaQb*uQa)               // Variance of vn
            + QaQb*VuQb/(4*uQb*uQa) 
            + QaQb*uQb*VuQa/(4*TMath::Power(uQa,3)) 
            + CuQbQaQb/(2*uQa)
            - uQb*CuQaQaQb/(2*TMath::Power(uQa,2))
            - QaQb*CuQauQb/(2*TMath::Power(uQa,2));
         vn=vnB;
         if( VvnB<0 ) {
            vnerror=VvnB;
            return false;
         }
         vnerror=TMath::Sqrt(VvnB);
      }
      break;
   default:                                       // POI flow
      {
         if( TMath::Abs(QaQb) < 1e-30*TMath::Abs(uQa*uQb) ) {  // Catch infinity
            vn=0;
            vnerror=-1;
            return false;
         }
         double vnP = TMath::Sqrt( uQa*uQb / QaQb );        // vn
         if(   TMath::Abs(uQa*QaQb) < 1e-30*TMath::Abs(uQb*VuQa) 
            || TMath::Abs(uQb*QaQb) < 1e-30*TMath::Abs(uQa*VuQb) 
            || TMath::Abs(QaQb) < 1e-30*TMath::Abs(uQa*uQb*VQaQb)
            || TMath::Abs(QaQb) < 1e-30*TMath::Abs(CuQauQb)
            || TMath::Abs(QaQb) < 1e-30*TMath::Abs(uQb*CuQaQaQb)
            || TMath::Abs(QaQb) < 1e-30*TMath::Abs(uQa*CuQbQaQb)
            ){
               vnerror=-98;
               return false;
         }
         double VvnP = uQb*VuQa/(4*uQa*QaQb)                // Variance of vn
            + uQa*VuQb/(4*uQb*QaQb) 
            + uQa*uQb*VQaQb/(4*TMath::Power(QaQb,3)) 
            + CuQauQb/(2*QaQb)
            - uQb*CuQaQaQb/(2*TMath::Power(QaQb,2))
            - uQa*CuQbQaQb/(2*TMath::Power(QaQb,2));
         vn=TMath::Sign(vnP,uQb);
         if( VvnP<0 ) {
            vnerror=VvnP;
            return false;
         }
         vnerror=TMath::Sqrt(VvnP);
      }
   }  // Switch between POI and subevents

   return (vnerror>=0);
}

void AliFlowAnalysisWithMSP::ReadHistograms(TDirectory *file)
{
   delete fCommonHist;     fCommonHist=0;    // Delete existing histograms
   delete fQaComponents;   fQaComponents=0;
   delete fQbComponents;   fQbComponents=0;
   delete fQaQb;           fQaQb=0;
   delete fPtUComponents;  fPtUComponents=0;
   delete fEtaUComponents; fEtaUComponents=0;   
   delete fAllStatistics;  fAllStatistics=0;
   delete fPtStatistics;   fPtStatistics=0;
   delete fEtaStatistics;  fEtaStatistics=0;
   delete fIntegratedFlow; fIntegratedFlow=0;
   delete fDiffFlowPt;     fDiffFlowPt=0;
   delete fDiffFlowEta;    fDiffFlowEta=0;
   delete fFlags;          fFlags=0;
   if( fHistList ) {
      fHistList->SetOwner(kFALSE);
      delete fHistList;
      fHistList=0;
   }

   file->GetObject("QaComponents",fQaComponents);
   file->GetObject("QbComponents",fQbComponents);
   file->GetObject("QaQb",fQaQb);
   file->GetObject("PtUComponents",fPtUComponents);
   file->GetObject("EtaUComponents",fEtaUComponents);
   file->GetObject("AllStatistics",fAllStatistics);
   file->GetObject("PtStatistics",fPtStatistics);
   file->GetObject("EtaStatistics",fEtaStatistics);
   if( !fQaComponents || !fQbComponents || !fQaQb || !fPtUComponents || !fEtaUComponents || !fAllStatistics || !fPtStatistics || !fEtaStatistics ) {
      std::cerr << "AliFlowAnalysisWithMSP::ReadHistograms() : One or more histograms were not read correctly from " << file->GetPath() << std::endl;
   }

   file->GetObject("Flags",fFlags);          // Flags are required

   if( !fFlags ){
      std::cerr << "AliFlowAnalysisWithMSP::ReadHistograms(TDirectoryFile *) : Flags histogram not found, using defaults" << std::endl;
      fHarmonic=2;
      fNUA=false;
   }else{
      fHarmonic=(UInt_t)(fFlags->GetBinContent(1));
      double harmonicSpread=fFlags->GetBinError(1);
      if( harmonicSpread!=0 ) {
         std::cerr << "AliFlowAnalysisWithMSP::ReadHistograms(TDirectoryFile *) :These histograms seem to be merged from analysis with different Harmonics. Results removed!" << std::endl;
         delete fIntegratedFlow; fIntegratedFlow=0;
         delete fDiffFlowPt;     fDiffFlowPt=0;
         delete fDiffFlowEta;    fDiffFlowEta=0;
      }
      fNUA=fFlags->GetBinContent(2);   // Mixing NUA does not matter since it needs to be recalculated anyway
      double nuaSpread=fFlags->GetBinError(2);
      if( nuaSpread!=0 ) {
         std::cerr << "AliFlowAnalysisWithMSP::ReadHistograms(TDirectoryFile *) :These histograms seem to be merged from analysis with different NUA corrections. Results removed" << std::endl;
         delete fIntegratedFlow; fIntegratedFlow=0;
         delete fDiffFlowPt;     fDiffFlowPt=0;
         delete fDiffFlowEta;    fDiffFlowEta=0;
      }
   }

   // Optional histograms, may return a zero pointer
   file->GetObject("AliFlowCommonHist_MSP",fCommonHist);             // The AliFlowCommonHist is optional
   file->GetObject("IntegratedFlow",fIntegratedFlow);                // Results are optional
   file->GetObject("DiffFlowPt",fDiffFlowPt);
   file->GetObject("DiffFlowEta",fDiffFlowEta);

   fBookCommonHistograms=(fCommonHist!=0);
}


void AliFlowAnalysisWithMSP::ReadHistograms(TList *list)
{
   if( !list ) return;
   delete fCommonHist;     fCommonHist=0;    // Delete existing histograms if any
   delete fQaComponents;   fQaComponents=0;
   delete fQbComponents;   fQbComponents=0;
   delete fQaQb;           fQaQb=0;
   delete fPtUComponents;  fPtUComponents=0;
   delete fEtaUComponents; fEtaUComponents=0;   
   delete fAllStatistics;  fAllStatistics=0;
   delete fPtStatistics;   fPtStatistics=0;
   delete fEtaStatistics;  fEtaStatistics=0;
   delete fIntegratedFlow; fIntegratedFlow=0;
   delete fDiffFlowPt;     fDiffFlowPt=0;
   delete fDiffFlowEta;    fDiffFlowEta=0;
   delete fFlags;          fFlags=0;
   if( fHistList ) {
      fHistList->SetOwner(kFALSE);
      delete fHistList;
      fHistList=0;
   }

   fQaComponents = static_cast<AliFlowMSPHistograms *>(list->FindObject("QaComponents"));
   fQbComponents = static_cast<AliFlowMSPHistograms *>(list->FindObject("QbComponents"));
   fQaQb = static_cast<AliFlowMSPHistograms *>(list->FindObject("QaQb"));
   fPtUComponents = static_cast<AliFlowMSPHistograms *>(list->FindObject("PtUComponents"));
   fEtaUComponents = static_cast<AliFlowMSPHistograms *>(list->FindObject("EtaUComponents"));
   fAllStatistics = static_cast<AliFlowMSPHistograms *>(list->FindObject("AllStatistics"));
   fPtStatistics = static_cast<AliFlowMSPHistograms *>(list->FindObject("PtStatistics"));
   fEtaStatistics = static_cast<AliFlowMSPHistograms *>(list->FindObject("EtaStatistics"));
   if( !fQaComponents || !fQbComponents || !fQaQb || !fPtUComponents || !fEtaUComponents || !fAllStatistics || !fPtStatistics || !fEtaStatistics ) {
      std::cerr << "AliFlowAnalysisWithMSP::ReadHistograms(Tlist *) : One or more histograms were not read correctly from TList" << std::endl;
   }

   fFlags = static_cast<TProfile *>(list->FindObject("Flags"));          // Flags are required

   if( !fFlags ){
      std::cerr << "AliFlowAnalysisWithMSP::ReadHistograms(TList *) : Flags histogram not found, using defaults" << std::endl;
      fHarmonic=2;
      fNUA=false;
   }else{
      fHarmonic=(UInt_t)(fFlags->GetBinContent(1));
      double harmonicSpread=fFlags->GetBinError(1);
      if( harmonicSpread!=0 ) {
         std::cerr << "These histograms seem to be merged from analysis with different Harmonics. Results removed!" << std::endl;
         delete fIntegratedFlow; fIntegratedFlow=0;
         delete fDiffFlowPt;     fDiffFlowPt=0;
         delete fDiffFlowEta;    fDiffFlowEta=0;
      }
      fNUA=fFlags->GetBinContent(2);   // Mixing NUA does not matter since it needs to be recalculated anyway
      double nuaSpread=fFlags->GetBinError(2);
      if( nuaSpread!=0 ) {
         std::cerr << "These histograms seem to be merged from analysis with different NUA corrections. Results removed" << std::endl;
         delete fIntegratedFlow; fIntegratedFlow=0;
         delete fDiffFlowPt;     fDiffFlowPt=0;
         delete fDiffFlowEta;    fDiffFlowEta=0;
      }
   }

   // Optional histograms, may return a zero pointer
   fCommonHist = static_cast<AliFlowCommonHist *>(list->FindObject("AliFlowCommonHist_MSP"));             // The AliFlowCommonHist is optional
   fIntegratedFlow = static_cast<TH1D *>(list->FindObject("IntegratedFlow"));                // Results are optional
   fDiffFlowPt = static_cast<TH1D *>(list->FindObject("DiffFlowPt"));
   fDiffFlowEta = static_cast<TH1D *>(list->FindObject("DiffFlowEta"));

   fBookCommonHistograms=(fCommonHist!=0);
}
