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
    fTrackDensity(title),
    fFMD1i(0),
    fFMD2i(0),
    fFMD2o(0),
    fFMD3i(0),
    fFMD3o(0),
    fOperComp(0)
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
}

//____________________________________________________________________
AliFMDMCSharingFilter::AliFMDMCSharingFilter(const AliFMDMCSharingFilter& o)
  : AliFMDSharingFilter(o), 
    fTrackDensity(o.fTrackDensity),
    fFMD1i(o.fFMD1i),
    fFMD2i(o.fFMD2i),
    fFMD2o(o.fFMD2o),
    fFMD3i(o.fFMD3i),
    fFMD3o(o.fFMD3o),
    fOperComp(o.fOperComp)
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
  fTrackDensity = o.fTrackDensity;
  return *this;
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


  fTrackDensity.Calculate(input, event, vz, output, primary);

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
  cd->SetOwner();
  cd->SetName("esd_mc_comparion");
  d->Add(cd);
  cd->Add(fFMD1i);
  cd->Add(fFMD2i);
  cd->Add(fFMD2o);
  cd->Add(fFMD3i);
  cd->Add(fFMD3o);
  cd->Add(fOperComp);
  fTrackDensity.DefineOutput(d);
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
  AliFMDSharingFilter::Print(option);
  gROOT->IncreaseDirLevel();
  fTrackDensity.Print(option);
  gROOT->DecreaseDirLevel();

}

//____________________________________________________________________
//
// EOF
//
