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
   : TNamed(), fHarmonic(2), fNUA(kFALSE), fBookCommonHistograms(kFALSE), fCommonHist(0), fQaComponents(0), 
   fQbComponents(0), fQaQb(0), fPtUComponents(0), fEtaUComponents(0), fAllStatistics(0), 
   fPtStatistics(0), fEtaStatistics(0),
   fIntegratedFlow(0), fDiffFlowPt(0), fDiffFlowEta(0), fFlags(0)
{
   // Default constructor. Intended for root IO purposes only
}     

AliFlowAnalysisWithMSP::AliFlowAnalysisWithMSP(TDirectoryFile *file)
   : TNamed(), fHarmonic(2), fNUA(kFALSE), fUseCommonConstants(kFALSE), fBookCommonHistograms(kFALSE), fCommonHist(0), fQaComponents(0), 
   fQbComponents(0), fQaQb(0), fPtUComponents(0), fEtaUComponents(0), fAllStatistics(0), 
   fPtStatistics(0), fEtaStatistics(0),
   fIntegratedFlow(0), fDiffFlowPt(0), fDiffFlowEta(0), fFlags(0)
{
   ReadHistograms(file);
}

AliFlowAnalysisWithMSP::AliFlowAnalysisWithMSP(const unsigned int harmonic, const bool commonConst, const bool commonHist) 
   : TNamed(), fHarmonic(harmonic), fNUA(kFALSE), fUseCommonConstants(commonConst), fBookCommonHistograms(commonHist), fCommonHist(0), fQaComponents(0), 
   fQbComponents(0), fQaQb(0), fPtUComponents(0), fEtaUComponents(0), fAllStatistics(0), 
   fPtStatistics(0), fEtaStatistics(0),
   fIntegratedFlow(0), fDiffFlowPt(0), fDiffFlowEta(0), fFlags(0)
{
   // Constructor defining the harmonic, usage of non uniform acceptance corrections and the AliFlowCommonHist() histograms
   // This is the constructor intended for the user
   SetNameTitle("MSP","Flow analysis with the Modified Scalar Product method");
}

AliFlowAnalysisWithMSP::~AliFlowAnalysisWithMSP() 
{
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
}

void AliFlowAnalysisWithMSP::Init()
{
   // Create all output objects. Memory consumption can be reduced by switching off some of the control histograms.
   // 
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

   // Default binning for histograms
   // TODO: allow variable binning

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
   event->Get2Qsub(flowVectors, fHarmonic);        // No phi, pt or eta weights implemented yet

   AliFlowVector &Qa=flowVectors[0];               // Define some mnemonics for the subevent flow vectors:
   AliFlowVector &Qb=flowVectors[1];

   const double QaW=Qa.GetMult();      // Weight for Qa and combinations of Qa
   const double QbW=Qb.GetMult();      // Weight for Qb and combinations of Qb

   if( QaW<2 || QbW<2 ) return;  // Require at least 2 particles in each subevent

   if(fCommonHist)fCommonHist->FillControlHistograms(event);   // Standard for all flow analysis


   const double qaxy[]={Qa.X()/QaW,Qa.Y()/QaW};                // Two variables expected
   const double wqaxy[]={QaW,QaW};
   fQaComponents->Fill(1, qaxy, wqaxy);                        // only one bin (all pt)

   const double qbxy[]={Qb.X()/QbW,Qb.Y()/QbW};                // Two variables expected
   const double wqbxy[]={QbW,QbW};
   fQbComponents->Fill(1, qbxy, wqbxy);                        // only one bin (all pt)

   const double QaQbW=QaW*QbW;
   const double weightedQaQb = (Qa*Qb)/QaQbW;                  // Scalar product of subevent Q vectors with weight
   fQaQb->Fill(1,&weightedQaQb,&QaQbW);                        // Average of QaQb per event

   
   int iTrack=0;
   while( AliFlowTrackSimple *track=event->GetTrack(iTrack++) ) {// Loop over the tracks in the event 
      // Get the track vector
      if(! track->InPOISelection() ) continue;
      const double trackWeight=track->Weight();
      const double phi=track->Phi();
      const double pt=track->Pt();
      const double eta=track->Eta();

      AliFlowVector u;
      u.SetMagPhi(1, fHarmonic*phi, trackWeight);
      u.SetMult(1);

      // Remove track from subevent a 
      AliFlowVector mQa(Qa);                 // Initialize with Qa flow vector
      if( track->InSubevent(0) ) {
         mQa-=u;                             // Should introduce phi weights here
      }

      // Remove track from subevent b 
      AliFlowVector mQb(Qb);                 // Initialize with Qb flow vector
      if( track->InSubevent(1) ) {
         mQb-=u;                             // Should introduce phi weights here
      }

      const double uQaW = mQa.GetMult();     // Weight is multiplicity of Q vector
      const double uQbW = mQb.GetMult();
      const double uQa=u*mQa;                // Correlate POI with subevent 
      const double uQb=u*mQb;                

      const double uxy[]={u.X(),u.Y()};
      const double wxy[]={1.0,1.0};
      fPtUComponents->Fill(pt, uxy, wxy);    // vs pt
      fEtaUComponents->Fill(eta, uxy, wxy);  // vs eta

      const double par[]={uQa/uQaW, uQb/uQbW, weightedQaQb};
      const double wgt[]={uQaW, uQbW, QaQbW};
      fAllStatistics->Fill(1,   par, wgt  );
      fPtStatistics->Fill(pt,   par, wgt );
      fEtaStatistics->Fill(eta, par, wgt );
   }
}

void AliFlowAnalysisWithMSP::Finish()
{
   // Calculate the final result from the stored correlations.
   // The NUA corrections are applied if the flag fNUA was set before the call to Finish()
   // If the output histograms already exist then they are replaced by the newly calculated result

   Print();       // Print a summary of the NUA terms and integrated flow
   // TODO: Create result histograms for integrated flow and store the flags to be able to restore the initial state from histograms only

   // Create result histograms
   fFlags->Fill("NUA",fNUA);
   delete fIntegratedFlow; fIntegratedFlow=0;   // First delete existing results (if any)
   delete fDiffFlowPt;     fDiffFlowPt=0;
   delete fDiffFlowEta;    fDiffFlowEta=0;

   fIntegratedFlow=new TH1D("IntegratedFlow","Integrated flow results",10,0.5,10.5);
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
      double vn=0;
      double evn=0;
      if( Calculate(vn, evn, fPtStatistics, fPtUComponents, i) && evn<1 ) {
         fDiffFlowPt->SetBinContent(i,vn);
         fDiffFlowPt->SetBinError(i,evn);
      }
   }

   int nbinseta=fEtaStatistics->Nbins();
   double etalow=fEtaStatistics->XLow();
   double etahigh=fEtaStatistics->XHigh();
   fDiffFlowEta= new TH1D("DiffFlowEta","flow of POI vs #eta",nbinseta,etalow,etahigh);
   fDiffFlowEta->SetDirectory(0);          // Do not automatically store in file
   fDiffFlowEta->SetStats(kFALSE);
   fDiffFlowEta->SetXTitle("p_{t}");
   fDiffFlowEta->SetYTitle(Form("v_{%d}",fHarmonic));

   for(int i=1; i<=nbinseta; ++i) {     
      double vn=0;
      double evn=0;
      if( Calculate(vn, evn, fEtaStatistics, fEtaUComponents, i) && evn<1 ) {
         fDiffFlowEta->SetBinContent(i,vn);
         fDiffFlowEta->SetBinError(i,evn);
      }
   }
}

void AliFlowAnalysisWithMSP::WriteHistograms(TDirectoryFile *file) const
{
   //file->Write(file->GetName(), TObject::kSingleKey);      // Make sure the directory itself is written

   if(fCommonHist) file->WriteTObject(fCommonHist);

   file->WriteTObject(fQaComponents,0,"Overwrite");        // Averages of Qa components per event
   file->WriteTObject(fQbComponents,0,"Overwrite");        // Averages of Qb components per event
   file->WriteTObject(fQaQb,0,"Overwrite");                // Average of QaQb per event
   file->WriteTObject(fPtUComponents,0,"Overwrite");       // u components vs pt
   file->WriteTObject(fEtaUComponents,0,"Overwrite");      // u components vs eta
   file->WriteTObject(fAllStatistics,0,"Overwrite");       // Integrated uQa, uQb and QaQa
   file->WriteTObject(fPtStatistics,0,"Overwrite");        // uQa, uQb and QaQb vs pt
   file->WriteTObject(fEtaStatistics,0,"Overwrite");       // uQa, uQb and QaQb vs eta

   if( fIntegratedFlow ) file->WriteTObject(fIntegratedFlow,0,"Overwrite");  // Integrated flow for POI and subevents
   if( fDiffFlowPt ) file->WriteTObject(fDiffFlowPt,0,"Overwrite");          // Differential flow vs pt if calculated
   if( fDiffFlowEta ) file->WriteTObject(fDiffFlowEta,0,"Overwrite");        // Differential flow vs eta if calculated

   file->WriteTObject(fFlags,0,"Overwrite");               // fHarmonic, fNUA

   file->WriteKeys();   // Make sure it happens now
}

void AliFlowAnalysisWithMSP::WriteCommonResults(TDirectoryFile *file) const
{
   // Copy the results to a AliFlowCommonHistResults() class and write to file
   // If the results were not calculated again then Finish() is called to generate them
   // The global AliFlowCommonConstants() is modified to adapt to the binning used in this analysis

   AliFlowCommonConstants *c=AliFlowCommonConstants::GetMaster();
   int nBinsPt=fPtStatistics->Nbins();
   double ptMin=fPtStatistics->XLow();
   double ptMax=fPtStatistics->XHigh();

   int nBinsEta=fEtaStatistics->Nbins();
   double etaMin=fEtaStatistics->XLow();
   double etaMax=fEtaStatistics->XHigh();

   c->SetNbinsPt(nBinsPt);
   c->SetPtMin(ptMin);
   c->SetPtMax(ptMax);
   c->SetNbinsEta(nBinsEta);
   c->SetEtaMin(etaMin);
   c->SetEtaMax(etaMax);

   bool oldAddStatus=TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
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

   TH1::AddDirectory(oldAddStatus);
}



void AliFlowAnalysisWithMSP::Print(const Option_t *opt)const
{
   std::cout << setprecision(3);   
   
   std::cout << "****************************************************" << std::endl;
   std::cout << "    Integrated flow from Modified Scalar Product    " << std::endl;
   std::cout << "                                                    " << std::endl;

   double vn=0;
   double vnerror=0;

   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 0);      
   std::cout << "v" << fHarmonic << " for POI       : " << setw(10) << vn << " +- " << setw(8) << vnerror << std::endl;
   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 1);
   std::cout << "v" << fHarmonic << " for subevent A: " << setw(10) << vn << " +- " << setw(8) << vnerror << std::endl;
   Calculate(vn, vnerror, fAllStatistics, fPtUComponents, 0, 2);
   std::cout << "v" << fHarmonic << " for subevent B: " << setw(10) << vn << " +- " << setw(8) << vnerror << std::endl;
   std::cout << std::endl;

   std::cout << "NUA terms: " << (fNUA?"(applied)":"(NOT applied)") << std::endl;

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

   std::cout << "<QaX>      " << setw(12) << fQaComponents->Average(0) << " +- " << setw(12) << TMath::Sqrt(fQaComponents->Variance(0)) << std::endl;
   std::cout << "<QaY>      " << setw(12) << fQaComponents->Average(1) << " +- " << setw(12) << TMath::Sqrt(fQaComponents->Variance(1)) << std::endl;
   std::cout << "<QbX>      " << setw(12) << fQbComponents->Average(0) << " +- " << setw(12) << TMath::Sqrt(fQbComponents->Variance(0)) << std::endl;
   std::cout << "<QbY>      " << setw(12) << fQbComponents->Average(1) << " +- " << setw(12) << TMath::Sqrt(fQbComponents->Variance(1)) << std::endl;
   std::cout << "<QaQb>     " << setw(12) << fQaQb->Average(0) << " +- " << setw(12) << TMath::Sqrt(fQbComponents->Variance(0)) << std::endl;
   std::cout << std::endl;

   std::cout << "Covariance matrix: " << std::endl;
   std::cout << "     " << setw(12) << "uQa"                           << setw(12) << "uQb"                           << setw(12) << "QaQb" << std::endl;
   std::cout << "uQa  " << setw(12) << fAllStatistics->Covariance(0,0) << std::endl;
   std::cout << "uQb  " << setw(12) << fAllStatistics->Covariance(1,0) << setw(12) << fAllStatistics->Covariance(1,1) << std::endl;
   std::cout << "QaQb " << setw(12) << fAllStatistics->Covariance(2,0) << setw(12) << fAllStatistics->Covariance(2,1) << setw(12) << fQaQb->Variance(0) << std::endl;
   std::cout << "****************************************************" << std::endl;
   std::cout << std::endl;
}


// private functions --------------------------------------------------------------------------------------
bool AliFlowAnalysisWithMSP::Calculate(double &vn, double &vnerror, const AliFlowMSPHistograms *hist, const AliFlowMSPHistograms *components, const int bin, const int poi) const
{
   // Get all averages and correlations need for the flow calculation
   double uQa=hist->Average(bin,0);                // <<uQa>>
   double VuQa=hist->Variance(bin,0);              // Var(<<uQa>>)
   double uQb=hist->Average(bin,1);                // <<uQb>>
   double VuQb=hist->Variance(bin,1);              // Var(<<uQb>>)
   double QaQb=fQaQb->Average(1,0);                // <QaQb> Should not be taken from hist(bin) buit from fQaQa because there QaQb is entered multiple times: <<QaQb>>!!
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
      // TODO: Check if not fully NUA correcting the error calculation is justified (but this is compatible with the original SP method)!
   }

   // Some sanity checks:
   if( uQa*uQb*QaQb <= 0 ) {                                    // Catch imaginary results
      vn=-99;
      vnerror=-99;
      return false;
   }

   // Decent results, calculate, print and store
   // TODO: The three cases are cyclic permutations of u, Qa, Qb so there should be a simpler way to do this.
   // However the difficulty is that we deal here with uQa, uQb and QaQb
   switch (poi) {
   case 1:                                        // Subevent A reference flow
      {
         if( TMath::Abs(uQb) < 1e-30*TMath::Abs(uQa*QaQb) ) {  // Protect against infinity
            vn=0;
            vnerror=-1;
            return false;
         }
         double vnB = TMath::Sqrt( uQa*QaQb / uQb ); // vn
         double VvnB = QaQb*VuQa/(4*uQa*uQb)         // Variance of vn
            + uQa*VQaQb/(4*QaQb*uQb) 
            + uQa*QaQb*VuQb/(4*TMath::Power(uQb,3)) 
            + CuQaQaQb/(2*uQb)
            - QaQb*CuQauQb/(2*TMath::Power(uQb,2))
            - uQa*CuQbQaQb/(2*TMath::Power(uQb,2));
         vn=vnB;
         if( VvnB<0 ) {
            vnerror=VvnB;
            return false;
         }
         vnerror=TMath::Sqrt(VvnB);
      }
      break;
   case 2:                                        // Subevent B reference flow
      {
         if( TMath::Abs(uQa) < 1e-30*TMath::Abs(uQb*QaQb) ) {  // Protect against infinity
            vn=0;
            vnerror=-1;
            return false;
         }
         double vnA = TMath::Sqrt( uQb*QaQb / uQa );        // vn
         double VvnA = uQb*VQaQb/(4*QaQb*uQa)               // Variance of vn
            + QaQb*VuQb/(4*uQb*uQa) 
            + QaQb*uQb*VuQa/(4*TMath::Power(uQa,3)) 
            + CuQbQaQb/(2*uQa)
            - uQb*CuQaQaQb/(2*TMath::Power(uQa,2))
            - uQa*CuQauQb/(2*TMath::Power(uQa,2));
         vn=vnA;
         if( VvnA<0 ) {
            vnerror=VvnA;
            return false;
         }
         vnerror=TMath::Sqrt(VvnA);
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
         vn=vnP;
         if( VvnP<0 ) {
            vnerror=VvnP;
            return false;
         }
         vnerror=TMath::Sqrt(VvnP);
      }
   }  // Switch between POI and subevents

   return (vnerror>=0);
}

void AliFlowAnalysisWithMSP::ReadHistograms(TDirectoryFile *file)
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

   // Optional histograms, may return a zero pointer
   file->GetObject("AliFlowCommonHist_MSP",fCommonHist);              // The AliFlowCommonHist is optional
   file->GetObject("IntegratedFlow",fIntegratedFlow);
   file->GetObject("DiffFlowPt",fDiffFlowPt);
   file->GetObject("DiffFlowEta",fDiffFlowEta);

   file->GetObject("Flags",fFlags);

   if( !fFlags ){
      std::cerr << "AliFlowAnalysisWithMSP::ReadHistograms() : Flags histogram not found, using defaults" << std::endl;
      fHarmonic=2;
      fNUA=false;
   }else{
      fHarmonic=fFlags->GetBinContent(1);
      double harmonicSpread=fFlags->GetBinError(1);
      if( harmonicSpread!=0 ) {
         std::cerr << "These histograms seem to be merged from analysis with different Harmonics. Results are invalid!" << std::endl;
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
   fBookCommonHistograms=(fCommonHist!=0);
}