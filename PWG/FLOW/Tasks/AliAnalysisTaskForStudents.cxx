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

/************************************** 
* template class for student projects * 
**************************************/ 
  
#include "Riostream.h"
#include "AliAnalysisTaskForStudents.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "TLegend.h"

#include "TCanvas.h" // TBI
#include "TFile.h" // TBI
#include "TRandom.h" 
#include "TGraphErrors.h"
#include "TComplex.h"
#include "TProfile.h"
#include "TStyle.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskForStudents)

//================================================================================================================

AliAnalysisTaskForStudents::AliAnalysisTaskForStudents(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(10.),
 fCentralityHist(NULL),
 fNCentralityBins(10),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 fMultHist(NULL),
 fPhiHist(NULL),
 fEtaHist(NULL),
 fhf(NULL),
 fhf2(NULL),
 leg_hist(NULL),
 // Final results:
 fFinalResultsList(NULL),
 num(0),
 fhr2bin(0), // me
 fsize(9),
 fbinnum(0),
 fcounter1(0),
 fcounter2(0),
 fcounter3(0),
 fmult(0),
 c2(0.),
 fhr2min(0.),
 fhr2max(0.),
 fper2(0.),
 fper3(0.),
 fper4(0.),
 fvarb(0.),
 fmu(0.),
 fymin(0.), // me
 fymax(0.), // me
 fc4(0.),
 fc3(0.),
 fc22(0.),
 fc24(0.),
 fc23(0.),
 fhr2(NULL),
 fhr3(NULL),
 fhrc22(NULL),
 fhrc23(NULL),
 fhrc24(NULL),
 fhrc3(NULL),
 fhrc4(NULL),
 fgr(NULL),
 fq6(0.,0.),
 fq4(0.,0.),
 fq(0.,0.),
 fq2(0.,0.),
 fqq(0.,0.),
 fq5(0.,0.),
 fq3(0.,0.),
 fq1(0.,0.),
 fq3x(0.,0.),
 fqq3(0.,0.)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskForStudents::AliAnalysisTaskForStudents(const char *name, Bool_t useParticleWeights)");

  // Base list:
  fHistList = new TList();
  fHistList->SetName("outputStudentAnalysis");
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays:
  this->InitializeArrays();

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  //DefineInput(0, AliFlowEventSimple::Class());  
  // Input slot #1 is needed for the weights input file:
  //if(useParticleWeights)
  //{
  // DefineInput(1, TList::Class());   
  //}  
  // Output slot #0 is reserved              
  // Output slot #1 writes into a TList container

  //num =1;

  DefineOutput(1, TList::Class());  

  if(useParticleWeights)
  {
   // not needed for the time being
  }

} // AliAnalysisTaskForStudents::AliAnalysisTaskForStudents(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskForStudents::AliAnalysisTaskForStudents(): 
 AliAnalysisTaskSE(),
 fHistList(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(10.),
 fCentralityHist(NULL),
 fNCentralityBins(10),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 fMultHist(NULL),
 fPhiHist(NULL),
 fEtaHist(NULL),
 fhf(NULL),
 fhf2(NULL),
 leg_hist(NULL),
 // Final results:
 fFinalResultsList(NULL),
 num(0),
 fhr2bin(0),
 fsize(9),
 fbinnum(0),
 fcounter1(0),
 fcounter2(0),
 fcounter3(0),
 fmult(0),
 c2(0.),
 fhr2min(0.),
 fhr2max(0.),
 fper2(0.),
 fper3(0.),
 fper4(0.),
 fvarb(0.),
 fmu(0.),
 fymin(0.), // me
 fymax(0.), // me
 fc4(0.),
 fc3(0.),
 fc22(0.),
 fc24(0.),
 fc23(0.),
 fhr2(NULL),
 fhr3(NULL),
 fhrc22(NULL),
 fhrc23(NULL),
 fhrc24(NULL),
 fhrc3(NULL),
 fhrc4(NULL),
 fgr(NULL),
 fq6(0.,0.),
 fq4(0.,0.),
 fq(0.,0.),
 fq2(0.,0.),
 fqq(0.,0.),
 fq5(0.,0.),
 fq3(0.,0.),
 fq1(0.,0.),
 fq3x(0.,0.),
 fqq3(0.,0.)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskForStudents::AliAnalysisTaskForStudents()");

} // AliAnalysisTaskForStudents::AliAnalysisTaskForStudents():

//================================================================================================================

AliAnalysisTaskForStudents::~AliAnalysisTaskForStudents()
{
 // Destructor.

 if(fHistList) delete fHistList;
  
} // AliAnalysisTaskForStudents::~AliAnalysisTaskForStudents()

//================================================================================================================

void AliAnalysisTaskForStudents::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Trick to avoid name clashes, part 1;
 // b) Book and nest all lists;
 // c) Book all objects;
 // *) Trick to avoid name clashes, part 2.
  
 // a) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);

 // b) Book and nest all lists:
 this->BookAndNestAllLists();

 // c) Book all objects:
 this->BookControlHistograms();
 this->BookFinalResultsHistograms();

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskForStudents::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskForStudents::UserExec(Option_t *) 
{
 Int_t im=0, jm=0;

  fq6 = TComplex(0.,0.);
  fq4 = TComplex(0.,0.);
  fq2 = TComplex(0.,0.);
  fq5 = TComplex(0.,0.);
  fq3 = TComplex(0.,0.);
  fq1 = TComplex(0.,0.);
  
 // Main loop (called for each event).

 // a) Get pointer to AOD event:
 // b) Start analysis over AODs;
 // c) Reset event-by-event objects;
 // d) PostData.


 // a) Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){return;}

 AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
 if(!ams){return;}
 if(ams->GetMultiplicityPercentile("V0M") >= fMinCentrality && ams->GetMultiplicityPercentile("V0M") < fMaxCentrality)
 {
  fCentralityHist->Fill(ams->GetMultiplicityPercentile("V0M"));
 }
 else
 {
  return; // this event do not belong to the centrality class specified for this particular analysis
 }


 // b) Start analysis over AODs:
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {
  AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to a track
  if(!aTrack){continue;} // protection against NULL pointers
  if(!aTrack->TestFilterBit(128)){continue;} // filter bit 128 denotes TPC-only tracks, use only them for the analysis

  // example variables for each track: (for more options, please see class /home/deniz/alicesw/aliroot/v5-07-20/src/STEER/AOD/AliAODTrack.h )  
  Double_t px = aTrack->Px(); // x-component of momenta
  Double_t py = aTrack->Py(); // y-component of momenta
  Double_t pz = aTrack->Pz(); // z-component of momenta
  Double_t e = aTrack->E();  // energy
  Double_t phi = aTrack->Phi(); // azimuthal angle
  Double_t eta = aTrack->Eta(); // pseudorapidity
  Double_t charge = aTrack->Charge(); // charge
  Double_t pt = aTrack->Pt(); // Pt
 
  // apply some cuts: e.g. take for the analysis only particles in -0.8 < eta < 0.8, and 0.2 < pT < 5.0
  // ... implementation of particle cuts ...
  if ( (-0.8 < eta) && (eta < 0.8) && (0.2 < pt) && (pt < 5.0)  ) {

        fq2 += TComplex(cos(2.*phi),sin(2.*phi));      
        fq4 += TComplex(cos(4.*phi),sin(4.*phi));
        fq6 += TComplex(cos(6.*phi),sin(6.*phi));
        fq1 += TComplex(cos(phi),sin(phi));
        fq3 += TComplex(cos(3.*phi),sin(3.*phi));
        fq5 += TComplex(cos(5.*phi),sin(5.*phi));

        im++;
    
        fPtHist->Fill(pt);
        fPhiHist->Fill(phi); 
        fEtaHist->Fill(eta);

  } // if ( (-0.8 < eta) && (eta < 0.8) && (0.2 < pT) && (pT < 5.0)  )


  // do some analysis only with the particles which passed the cuts
  // ... your analysis code ... 

 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
  
    fmult = im;

    fMultHist->Fill(fmult);

    if ( fmult > 3) {
      fper2 = fmult*(fmult-1.);
      fper4 = fmult*(fmult-1.)*(fmult-2.)*(fmult-3.);
      //fper3 = fmult*(fmult-1)*(fmult-2);

      fc22=(fq2.Rho2()-fmult)/fper2;
      fhr2->Fill(0.1,fc22,fper2);
   
      fc24=(fq4.Rho2()-fmult)/fper2;
      fhr2->Fill(1.1,fc24,fper2);

      fq=fq6*TComplex::Conjugate(fq2)*TComplex::Conjugate(fq4);
      fqq=fq4*TComplex::Conjugate(fq2)*TComplex::Conjugate(fq2);

      fc4=(fq2.Rho2()*fq4.Rho2() - 2.*fq.Re() - 2.*fqq.Re() + fq6.Rho2()+fq2.Rho2()-(fmult-4.)*(fq2.Rho2()+fq4.Rho2())+fmult*(fmult-6.))/fper4;
      fhr2->Fill(2.1,fc4,fper4);

      fc23=(fq3.Rho2()-fmult)/fper2;
      fhr3->Fill(0.1,fc23,fper2);

      fq3x=fq5*TComplex::Conjugate(fq2)*TComplex::Conjugate(fq3);
      fqq3=fq3*TComplex::Conjugate(fq1)*TComplex::Conjugate(fq2);
   
      fc3=(fq2.Rho2()*fq3.Rho2() - 2.*fq3x.Re() - 2.*fqq3.Re() + fq5.Rho2()+fq2.Rho2()-(fmult-4.)*(fq2.Rho2()+fq3.Rho2())+fmult*(fmult-6.))/fper4;
      fhr3->Fill(0.6,fc3,fper4);

       //cout<<endl<<"Aod disinda "<<fcounter1++<<"  "<<fc3<<"  "<<fc23<<"  "<<fc22<<"  "<< fc23*fc22<<endl;
       //cout<<endl<<"Aod disinda "<<fcounter1++<<"  "<<fc3<<"  "<<fc23<<"  "<<fc22<<"  "<< fc23*fc22<<endl;
    
      fbinnum = gRandom->Uniform(0,10);
      fhrc22->Fill(fbinnum,fc22,fper2);
      fhrc24->Fill(fbinnum,fc24,fper2);
      fhrc4->Fill(fbinnum,fc4,fper4);

      fhrc23->Fill(fbinnum,fc23,fper2);
      fhrc3->Fill(fbinnum,fc3,fper4);
    } // if ( fmult > 3)
 

 // c) Reset event-by-event objects:
 // ...

 // d) PostData:
 PostData(1,fHistList);

} // void AliAnalysisTaskForStudents::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskForStudents::Terminate(Option_t *)
{

  Double_t varb4=0., mu4=0., varb3=0., mu3=0., sc42=0., sc32=0.;

  sc42 = fhr2->GetBinContent(3)-fhr2->GetBinContent(1)*fhr2->GetBinContent(2);
  sc32 = fhr3->GetBinContent(2)-fhr2->GetBinContent(1)*fhr3->GetBinContent(1);
 // Accessing the merged output list.


 // access average two particle correlation from TProfile, that you have filled in UserExec
 // raise that two particle correlation to power 7
 // store that information in the histogram holding final results

 cout<<endl<<" fhr2(3) "<<fhr2->GetBinContent(3)<<" "<<fhr2->GetBinContent(2)<<"  "<<fhr2->GetBinContent(1)<<" sc(4,2) "<<sc42<<endl;
 fhf->SetBinContent(1, sc32);
 fhf->SetBinContent(2,sc42);

  for (Int_t i=1;i<11;i++) { mu4 +=fhrc4->GetBinContent(i)-fhrc24->GetBinContent(i)*fhrc22->GetBinContent(i); }
  for (Int_t s=1;s<11;s++) {
    varb4+= pow(fhrc4->GetBinContent(s)-fhrc24->GetBinContent(s)*fhrc22->GetBinContent(s) - mu4/10.,2);
  }
 fhf->SetBinError(2, sqrt(varb4/9.));

 for (Int_t i=1;i<11;i++) { mu3 +=fhrc3->GetBinContent(i)-fhrc23->GetBinContent(i)*fhrc22->GetBinContent(i); }
  for (Int_t s=1;s<11;s++) {
    varb3+= pow(fhrc3->GetBinContent(s)-fhrc23->GetBinContent(s)*fhrc22->GetBinContent(s) - mu3/10.,2);
  }
 fhf->SetBinError(1,sqrt(varb3/9.));
 
 cout<<endl<<" sc32 "<<fhf->GetBinContent(1)<<" "<<" sc42 "<<fhf->GetBinContent(2)<<endl;

 fHistList = (TList*)GetOutputData(1);
 if(!fHistList){exit(1);}

 // Do some calculation in offline mode here:
 // ...

 TFile *f = new TFile("AnalysisResultsSC.root","RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete f;

} // end of void AliAnalysisTaskForStudents::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskForStudents::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.

 for(Int_t i=0;i<9;i++)
 {
  fsc4[i] = 0.;
  fcentral[i] = 0.;
  fyerr[i] = 0.;
 }

 
 // end

} // void AliAnalysisTaskForStudents::InitializeArrays()

//=======================================================================================================================

void AliAnalysisTaskForStudents::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskForStudents::BookAndNestAllLists()";
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is NULL");}

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("ControlHistograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

 // b) Book and nest lists for final results:
 fFinalResultsList = new TList();
 fFinalResultsList->SetName("FinalResults");
 fFinalResultsList->SetOwner(kTRUE);
 fHistList->Add(fFinalResultsList);

} // void AliAnalysisTaskForStudents::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskForStudents::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt spectra;
 // b) ...

 // a) Book histogram to hold pt spectra:
 fPtHist = new TH1F("fPtHist","atrack->Pt()",fNbins,fMinBin,fMaxBin);
 fPtHist->SetStats(kFALSE);
 fPtHist->SetFillColor(kBlue-10);
 fPtHist->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHist);

 fCentralityHist = new TH1F("fCentralityHist","ams->GetMultiplicityPercentile(\"V0M\")",fNCentralityBins,fMinCentrality,fMaxCentrality);
 fCentralityHist->SetFillColor(kBlue-10);
 fCentralityHist->GetXaxis()->SetTitle("centrality percentile");
 fControlHistogramsList->Add(fCentralityHist);

 fMultHist = new TH1F("fMultHist","Multiplicity Distribution",1000,0,3000);
 fMultHist->GetXaxis()->SetTitle("m");
 fMultHist->SetLineColor(4);
 fControlHistogramsList->Add(fMultHist);

 fPhiHist = new TH1F("fPhiHist","Phi Distribution",1000,0.,6.3);
 fPhiHist->GetXaxis()->SetTitle("Phi");
 fPhiHist->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHist);

 fEtaHist = new TH1F("fEtaHist","Eta Distribution",1000,-1.,1.);
 fEtaHist->GetXaxis()->SetTitle("Eta");
 fEtaHist->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHist);
 //fControlHistogramsList->Add(fhf);
 

 fgr = new TGraphErrors(fsize,fcentral,fsc4,0,fyerr);
 fgr->SetName("gr");
 fgr->SetMarkerColor(2);
 fgr->SetMarkerStyle(21);
 fgr->SetMinimum(fymin);
 fgr->SetMaximum(fymax);
 fgr->SetMarkerSize(1.5);
 fgr->GetYaxis()->SetTitleOffset(1.2);
 fgr->GetXaxis()->SetTitle("centrality");
 fgr->GetYaxis()->SetTitle("SC(4,2)");

 fhr2 = new TProfile("fhr2","fhr2 ",3,0.,3.);
  fhr2->Sumw2();
  fControlHistogramsList->Add(fhr2);
 
 fhr3 = new TProfile("fhr3","fhr3 ",2,0.,1.);
  fhr3->Sumw2();
  fControlHistogramsList->Add(fhr3);

 fhrc22 = new TProfile("fhrc22","fhrc22 ",10,0.,10);
  fhrc22->Sumw2();
 fControlHistogramsList->Add(fhrc22);
 fhrc23 = new TProfile("fhrc23","fhrc23 ",10,0.,10);
  fhrc23->Sumw2();
  fControlHistogramsList->Add(fhrc23);
 fhrc24 = new TProfile("fhrc24","fhrc24 ",10,0.,10);
  fhrc24->Sumw2();
  fControlHistogramsList->Add(fhrc24);
 fhrc4 = new TProfile("fhrc4","fhrc4 ",10,0.,10);
  fhrc4->Sumw2();
  fControlHistogramsList->Add(fhrc4);
 fhrc3 = new TProfile("fhrc3","fhrc3 ",10,0.,10);
  fhrc3->Sumw2();
  fControlHistogramsList->Add(fhrc3);

 
 // b) ...

} // void AliAnalysisTaskForStudents::BookControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskForStudents::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

 fhf = new TH1F("fhf","Profil ",2,0.,100);
 fhf->GetXaxis()->SetTitle(" Centrality percentile");
 fhf->GetYaxis()->SetTitle("SC(m,n)");
 fhf->SetTitle(" ");
 fhf->SetMarkerColor(4);
 fhf->SetLineColor(1);
 fhf->SetLineWidth(1.8);
 fhf->SetMarkerSize(1.3);
 fhf->SetMarkerStyle(8);
 fhf->SetStats(0);
 //gStyle->SetErrorX(0.0);
 fhf->SetOption("pex0");

 fFinalResultsList->Add(fhf);

} // void AliAnalysisTaskForStudents::BookFinalResultsHistograms()

//=======================================================================================================================

