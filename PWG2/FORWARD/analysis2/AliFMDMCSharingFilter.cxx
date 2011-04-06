//
// Class to do the sharing correction for MC data.
//
// Input: 
//    - AliESDFMD object  - from reconstruction
//    - Kinematics
//    - Track-References
//
// Output: 
//    - AliESDFMD object  - copy of input, but with signals merged 
//
// Corrections used: 
//    - ELoss fits
//
// Histograms: 
//    - For each ring (FMD1i, FMD2i, FMD2o, FMD3i, FMD3o) the distribution of 
//      signals before and after the filter.  
//    - For each ring (see above), an array of distributions of number of 
//      hit strips for each vertex bin (if enabled - see Init method)
// 
#include "AliFMDMCSharingFilter.h"
#include <AliESDFMD.h>
#include <AliMCEvent.h>
#include <AliTrackReference.h>
#include <AliStack.h>
#include <TAxis.h>
#include <TList.h>
#include <TH1.h>
#include <TMath.h>
#include "AliFMDStripIndex.h"
#include "AliFMDFloatMap.h"
#include <AliLog.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>

ClassImp(AliFMDMCSharingFilter)
#if 0
; // This is for Emacs
#endif 


//____________________________________________________________________
AliFMDMCSharingFilter::AliFMDMCSharingFilter(const char* title)
  : AliFMDSharingFilter(title), 
    fFMD1i(0),
    fFMD2i(0),
    fFMD2o(0),
    fFMD3i(0),
    fFMD3o(0),
    fSumEta(0),
    fOperComp(0),
    fThetaVsNr(0), 
    fOnlyPrimary(false)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    title Title of object  - not significant 
  //
  fFMD1i = new TH2D("FMD1i_corr", "Merged vs MC", 21, -.5, 20.5, 300, 0, 15);
  fFMD2i = new TH2D("FMD2i_corr", "Merged vs MC", 21, -.5, 20.5, 300, 0, 15);
  fFMD2o = new TH2D("FMD2o_corr", "Merged vs MC", 21, -.5, 20.5, 300, 0, 15);
  fFMD3i = new TH2D("FMD3i_corr", "Merged vs MC", 21, -.5, 20.5, 300, 0, 15);
  fFMD3o = new TH2D("FMD3o_corr", "Merged vs MC", 21, -.5, 20.5, 300, 0, 15);
  fFMD1i->SetYTitle("#Delta E/#Delta_{mip} (ESD)");
  fFMD1i->SetXTitle("Hits (MC)");
  fFMD2i->SetYTitle("#Delta E/#Delta_{mip} (ESD)");
  fFMD2i->SetXTitle("Hits (MC)");
  fFMD2o->SetYTitle("#Delta E/#Delta_{mip} (ESD)");
  fFMD2o->SetXTitle("Hits (MC)");
  fFMD3i->SetYTitle("#Delta E/#Delta_{mip} (ESD)");
  fFMD3i->SetXTitle("Hits (MC)");
  fFMD3o->SetYTitle("#Delta E/#Delta_{mip} (ESD)");
  fFMD3o->SetXTitle("Hits (MC)");
  fFMD1i->SetDirectory(0);
  fFMD2i->SetDirectory(0);
  fFMD2o->SetDirectory(0);
  fFMD3i->SetDirectory(0);
  fFMD3o->SetDirectory(0);
  fSumEta = new TH1D("mcSumEta", "MC INEL Truth", 200, -4, 6);
  fSumEta->SetXTitle("#eta");
  fSumEta->SetYTitle("dN_{ch}/d#eta");
  fSumEta->SetDirectory(0);
  fSumEta->Sumw2();
  fSumEta->SetMarkerColor(kOrange+2);
  fSumEta->SetMarkerStyle(22);
  fSumEta->SetFillColor(0);
  fSumEta->SetFillStyle(0);

  fOper     = new AliFMDFloatMap(0,0,0,0);
  fOperComp = new TH2I("operComp", "Operation vs # track refs", 
		       kMergedInto, kNone-.5, kMergedInto+.5, 
		       20, -.5, 19.5);
  fOperComp->SetXTitle("Operation");
  fOperComp->SetYTitle("# of track refs in sector");
  fOperComp->SetZTitle("Observations");
  fOperComp->GetXaxis()->SetBinLabel(kNone,            "None");
  fOperComp->GetXaxis()->SetBinLabel(kCandidate,       "Candidate");
  fOperComp->GetXaxis()->SetBinLabel(kMergedWithOther, "Merged w/other");
  fOperComp->GetXaxis()->SetBinLabel(kMergedInto,      "Merged into");
  fOperComp->SetDirectory(0);
  
  fThetaVsNr = new TH2D("thetaVsNr", "#theta of track vs # track references",
			360, 0, 360, 20, -.5, 19.5);
  fThetaVsNr->SetXTitle("#theta [degrees]");
  fThetaVsNr->SetYTitle("# of track references");
  fThetaVsNr->SetDirectory(0);
}

//____________________________________________________________________
AliFMDMCSharingFilter::AliFMDMCSharingFilter(const AliFMDMCSharingFilter& o)
  : AliFMDSharingFilter(o), 
    fFMD1i(o.fFMD1i),
    fFMD2i(o.fFMD2i),
    fFMD2o(o.fFMD2o),
    fFMD3i(o.fFMD3i),
    fFMD3o(o.fFMD3o),
    fSumEta(o.fSumEta),
    fOperComp(o.fOperComp),
    fThetaVsNr(o.fThetaVsNr),
    fOnlyPrimary(o.fOnlyPrimary)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}

//____________________________________________________________________
AliFMDMCSharingFilter::~AliFMDMCSharingFilter()
{
  // 
  // Destructor
  //
  if (fFMD1i)  delete fFMD1i;
  if (fFMD2i)  delete fFMD2i;
  if (fFMD2o)  delete fFMD2o;
  if (fFMD3i)  delete fFMD3i;
  if (fFMD3o)  delete fFMD3o;
  if (fSumEta) delete fSumEta;
}

//____________________________________________________________________
AliFMDMCSharingFilter&
AliFMDMCSharingFilter::operator=(const AliFMDMCSharingFilter& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this 
  //
  AliFMDSharingFilter::operator=(o);
  fOnlyPrimary = o.fOnlyPrimary;
  return *this;
}

//____________________________________________________________________
void
AliFMDMCSharingFilter::StoreParticle(UShort_t   d, 
				     Char_t     r,
				     UShort_t   s, 
				     UShort_t   t, 
				     UShort_t   nR,
				     Double_t   theta,
				     AliESDFMD& output) const
{
  // 
  // Store a particle hit in FMD<i>dr</i>[<i>s,t</i>] in @a output
  // 
  // Parameters:
  //    d       Detector
  //    r       Ring
  //    s       Sector
  //    t       Strip
  //    nR      Number of references to this particle in this sector
  //    output  Output ESD object
  //
  Double_t old = output.Multiplicity(d,r,s,t);
  if (old == AliESDFMD::kInvalidMult) old = 0;
  if (fOper) fOperComp->Fill(fOper->operator()(d,r,s,t), nR);
  if (theta < 0) theta += 2*TMath::Pi();
  theta *= 180. / TMath::Pi();
  fThetaVsNr->Fill(theta, nR);
  output.SetMultiplicity(d,r,s,t,old+1);
}

//____________________________________________________________________
Bool_t
AliFMDMCSharingFilter::FilterMC(const AliESDFMD&  input, 
				const AliMCEvent& event,
				Double_t          vz,
				AliESDFMD&        output, 
				TH2D*             primary)
{
  // 
  // Filter the input kinematics and track references, using 
  // some of the ESD information
  // 
  // Parameters:
  //    input   Input ESD event
  //    event   Input MC event
  //    vz      Vertex position 
  //    output  Output ESD-like object
  //    primary Per-event histogram of primaries 
  //
  // Return:
  //    True on succes, false otherwise 
  //
  output.Clear();

  // Increase event count - stored in 
  // underflow bin 
  fSumEta->AddBinContent(0); 

  // Copy eta values to output 
  for (UShort_t ed = 1; ed <= 3; ed++) { 
    UShort_t nq = (ed == 1 ? 1 : 2);
    for (UShort_t eq = 0; eq < nq; eq++) {
      Char_t   er = (eq == 0 ? 'I' : 'O');
      UShort_t ns = (eq == 0 ?  20 :  40);
      UShort_t nt = (eq == 0 ? 512 : 256);
      for (UShort_t es = 0; es < ns; es++) 
	for (UShort_t et = 0; et < nt; et++) 
	  output.SetEta(ed, er, es, et, input.Eta(ed, er, es, et));
    }
  }
  AliStack* stack = const_cast<AliMCEvent&>(event).Stack();
  Int_t nTracks   = stack->GetNtrack();//event.GetNumberOfTracks();
  Int_t nPrim     = stack->GetNprimary();//event.GetNumberOfPrimary();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { 
    AliMCParticle* particle = 
      static_cast<AliMCParticle*>(event.GetTrack(iTr));
    
    // Check the returned particle 
    if (!particle) continue;
    
    // Check if this charged and a primary 
    Bool_t isCharged = particle->Charge() != 0;
    if (!isCharged) continue;
    Bool_t isPrimary = stack->IsPhysicalPrimary(iTr);

    // Pseudo rapidity and azimuthal angle 
    Double_t eta = particle->Eta();
    Double_t phi = particle->Phi();
    
    // Fill 'dn/deta' histogram 
    if (isPrimary && iTr < nPrim) { 
      // Avoid under count - used to store event count
      if (eta >= fSumEta->GetXaxis()->GetXmin()) fSumEta->Fill(eta);
      primary->Fill(eta, phi);
    }

    // Bail out if we're only processing primaries - perhaps we should
    // track back to the original primary?
    if (fOnlyPrimary && !isPrimary) continue;

    Int_t    nTrRef  = particle->GetNumberOfTrackReferences();
    Int_t    longest = -1;
    Double_t angle   = 0;
    UShort_t oD = 0, oS = 1024, oT = 1024;
    Char_t   oR = '\0';
    UShort_t nC = 0;
    Double_t oTheta = 0;
    for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) { 
      AliTrackReference* ref = particle->GetTrackReference(iTrRef);
      
      // Check existence 
      if (!ref) continue;

      // Check that we hit an FMD element 
      if (ref->DetectorId() != AliTrackReference::kFMD) 
	continue;

      // Count number of track refs in this sector 
      nC++;

      // Get the detector coordinates 
      UShort_t d, s, t;
      Char_t r;
      AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);
      // If this is a new detector/ring, then reset the other one 
      if (oD > 0 && oR != '\0' && oS != 1024 && 
	  (d != oD || r != oR || s != oS)) {
	longest = -1;
	angle   = 0;
	StoreParticle(oD, oR, oS, oT, nC, oTheta, output);
	nC = 0;
	oD = 0;
	oR = '\0';
	oS = 1024;
	oT = 1024;
      }

      oD = d;
      oR = r;
      oS = s;
      oT = t;

      // The longest passage is determined through the angle 
      Double_t x    = ref->X();
      Double_t y    = ref->Y();
      Double_t z    = ref->Z()-vz;
      Double_t rr   = TMath::Sqrt(x*x+y*y);
      Double_t theta= TMath::ATan2(rr,z);
      Double_t ang  = TMath::Abs(TMath::Pi()-theta);
      if (ang > angle) {
	longest = iTrRef;
	angle   = ang;
      }
      oTheta = theta;
    } // Loop over track references
    if (longest < 0) continue;

    // Get the reference corresponding to the longest path through the detector
    AliTrackReference* ref = particle->GetTrackReference(longest);

    // Get the detector coordinates 
    UShort_t d, s, t;
    Char_t r;
    AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);
    
    StoreParticle(d,r,s,t,nC,particle->Theta(),output);
  } // Loop over tracks
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDMCSharingFilter::CompareResults(const AliESDFMD&  esd, 
				      const AliESDFMD&  mc)
{
  // 
  // Compare the result of merging to the monte-carlo truth.  This
  // fills the correlation histograms
  // 
  // Parameters:
  //    esd  ESD after sharing correction
  //    mc   MC ESD 
  //

  // Copy eta values to output 
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nq = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nq; q++) {
      Char_t   r  = (q == 0 ? 'I' : 'O');
      UShort_t ns = (q == 0 ?  20 :  40);
      UShort_t nt = (q == 0 ? 512 : 256);
      TH2*     co = 0;
      switch (d) { 
      case 1: co = fFMD1i; break;
      case 2: co = (q == 0 ? fFMD2i : fFMD2o); break;
      case 3: co = (q == 0 ? fFMD3i : fFMD3o); break;
      }

      for (UShort_t s = 0; s < ns; s++) {
	for (UShort_t t = 0; t < nt; t++) { 
	  Float_t mEsd = esd.Multiplicity(d, r, s, t);
	  Float_t mMc  = mc.Multiplicity(d, r, s, t);

	  co->Fill(mMc, mEsd);
	} 
      }
    }
  }
}
  
//____________________________________________________________________
void
AliFMDMCSharingFilter::DefineOutput(TList* dir)
{
  // 
  // Define the output histograms.  These are put in a sub list of the
  // passed list.   The histograms are merged before the parent task calls 
  // AliAnalysisTaskSE::Terminate 
  // 
  // Parameters:
  //    dir Directory to add to 
  //
  AliFMDSharingFilter::DefineOutput(dir);
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  TList* cd = new TList;
  cd->SetName("esd_mc_comparion");
  d->Add(cd);
  cd->Add(fFMD1i);
  cd->Add(fFMD2i);
  cd->Add(fFMD2o);
  cd->Add(fFMD3i);
  cd->Add(fFMD3o);
  dir->Add(fSumEta);
  cd->Add(fOperComp);
  cd->Add(fThetaVsNr);
}

//____________________________________________________________________
void
AliFMDMCSharingFilter::ScaleHistograms(const TList* dir, Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // 
  // Parameters:
  //    dir     Where the output is 
  //    nEvents Number of events 
  //
  AliFMDSharingFilter::ScaleHistograms(dir, nEvents);
  TH1D* sumEta = static_cast<TH1D*>(dir->FindObject("mcSumEta"));
  if (!sumEta) { 
    AliWarning(Form("No mcSumEta histogram found in output list"));
    return;
  }
  Double_t n = nEvents; // sumEta->GetBinContent(0);
  sumEta->Scale(1. / n, "width");
}

//____________________________________________________________________
void
AliFMDMCSharingFilter::Print(Option_t* option) const
{
  // 
  // Print information
  // 
  // Parameters:
  //    option Not used 
  //
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  AliFMDSharingFilter::Print(option);
  std::cout << std::boolalpha 
	    << ind << " Only primary tracks:    " << fOnlyPrimary 
	    << std::noboolalpha << std::endl;
}

//____________________________________________________________________
//
// EOF
//
