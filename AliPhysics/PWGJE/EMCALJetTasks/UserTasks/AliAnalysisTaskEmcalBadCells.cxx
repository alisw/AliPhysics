//
// Bad cells analysis task.
//
// Author: M. Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliVCaloCells.h"

#include "AliAnalysisTaskEmcalBadCells.h"

ClassImp(AliAnalysisTaskEmcalBadCells)

//________________________________________________________________________
AliAnalysisTaskEmcalBadCells::AliAnalysisTaskEmcalBadCells() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskEmcalBadCells"),
  fh2AmplitudeCellNumber(0x0)
{
  // Default constructor.


  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalBadCells::AliAnalysisTaskEmcalBadCells(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fh2AmplitudeCellNumber(0x0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(histo);
}

//________________________________________________________________________
AliAnalysisTaskEmcalBadCells::~AliAnalysisTaskEmcalBadCells()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalBadCells::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fh2AmplitudeCellNumber = new TH2F("fh2AmplitudeCellNumber","fh2AmplitudeCellNumber",11520,0,11520,100,0,10);
  fOutput->Add(fh2AmplitudeCellNumber);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalBadCells::FillHistograms()
{
  // Fill histograms.

  if(!fGeom)
    fGeom = AliEMCALGeometry::GetInstance();

  Short_t absId = -1;
  Int_t nCells  =  0;

  Int_t absIdMin = 1000000;
  Int_t absIdMax = -1;

  if (fCaloCells) {
    if(fCaloCells->IsEMCAL()){
      for (Int_t icell = 0; icell <  fCaloCells->GetNumberOfCells(); icell++) {

	nCells++;

	Double_t amp =0., time = 0., efrac = 0;
	Int_t mclabel = -1;

	fCaloCells->GetCell(icell, absId, amp, time,mclabel,efrac);
	if(absId<absIdMin) absIdMin=absId;
	if(absId>absIdMax) absIdMax=absId;

	fh2AmplitudeCellNumber->Fill(absId,amp);

      }
      //      AliInfo(Form("%s: absId min: %d  max: %d  nCells: %d",GetName(),absIdMin,absIdMax,nCells));
    }

  }
  
 

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalBadCells::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalBadCells::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
