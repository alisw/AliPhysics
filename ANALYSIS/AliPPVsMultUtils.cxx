/*****************************************************************
 
 General utility class for doing pp vs multiplicity studies
 
 --- Functionality being added on demand.
 
 --- Please report any bugs, comments, suggestions to:
 david.dobrigkeit.chinellato@cern.ch
 
 *****************************************************************/

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include <TMatrixDSym.h>
#include <TMath.h>
#include "AliESDUtils.h"
#include "AliPPVsMultUtils.h"
#include <TFile.h>


ClassImp(AliPPVsMultUtils)

//______________________________________________________________________
AliPPVsMultUtils::AliPPVsMultUtils():TObject(),
fRunNumber(0),
fCalibrationLoaded(0),
fBoundaryHisto_V0M(0),
fBoundaryHisto_V0A(0),
fBoundaryHisto_V0C(0),
fBoundaryHisto_V0MEq(0),
fBoundaryHisto_V0AEq(0),
fBoundaryHisto_V0CEq(0)
{
    // Default contructor
}


//_____________________________________________________________________________
AliPPVsMultUtils::AliPPVsMultUtils(const AliPPVsMultUtils &c) : TObject(c),
fRunNumber(0),
fCalibrationLoaded(0),
fBoundaryHisto_V0M(0),
fBoundaryHisto_V0A(0),
fBoundaryHisto_V0C(0),
fBoundaryHisto_V0MEq(0),
fBoundaryHisto_V0AEq(0),
fBoundaryHisto_V0CEq(0)
{
    //
    // copy constructor - untested
    //
    ((AliPPVsMultUtils &) c).Copy(*this);
}

//_____________________________________________________________________________
AliPPVsMultUtils &AliPPVsMultUtils::operator=(const AliPPVsMultUtils &c)
{
    //
    // Assignment operator - untested
    //
    
    if (this != &c) ((AliPPVsMultUtils &) c).Copy(*this);
    return *this;
}

//______________________________________________________________________
Float_t AliPPVsMultUtils::GetMultiplicityPercentile(AliVEvent *event, TString lMethod)
// Function to get multiplicity quantiles 
{
    Int_t lRequestedRunNumber = event->GetRunNumber();
    if( lRequestedRunNumber != fRunNumber ) fCalibrationLoaded = LoadCalibration( lRequestedRunNumber );
    
    if ( !fCalibrationLoaded ){ 
      return -1000; //Return absurd value (hopefully noone will use this...) 
    }
    
    Float_t lreturnval = -1;
    
    //Get VZERO Information for multiplicity later
    AliVVZERO* esdV0 = event->GetVZEROData();
    if (!esdV0) {
        AliError("AliVVZERO not available");
        return -1;
    }
    
    // VZERO PART
    Float_t  multV0A  = 0;            //  multiplicity from V0 reco side A
    Float_t  multV0C  = 0;            //  multiplicity from V0 reco side C
    Float_t  multV0AEq  = 0;          //  multiplicity from V0 reco side A
    Float_t  multV0CEq  = 0;          //  multiplicity from V0 reco side C
    Float_t  multV0ACorr  = 0;            //  multiplicity from V0 reco side A
    Float_t  multV0CCorr  = 0;            //  multiplicity from V0 reco side C
    
    //Non-Equalized Signal: copy of multV0ACorr and multV0CCorr from AliCentralitySelectionTask
    //Getters for uncorrected multiplicity
    multV0A=esdV0->GetMTotV0A();
    multV0C=esdV0->GetMTotV0C();
   
    //Getting around to the SPD vertex -> typecast to ESD/AOD
    const AliVVertex *lPrimarySPDVtx = NULL;
    /* get ESD vertex SPD */
    if (event->InheritsFrom("AliESDEvent")) {
      AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
      if (!esdevent) return kFALSE;
      lPrimarySPDVtx = esdevent->GetPrimaryVertexSPD();
    }
    /* get AOD vertex SPD */
    else if (event->InheritsFrom("AliAODEvent")) {
      AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
      if (!aodevent) return kFALSE;
      lPrimarySPDVtx = aodevent->GetPrimaryVertexSPD();
    }
   
    //Get Z vertex position of SPD vertex (why not Tracking if available?)
    Float_t zvtx = lPrimarySPDVtx->GetZ();
    
    //Acquire Corrected multV0A
    multV0ACorr = AliESDUtils::GetCorrV0A(multV0A,zvtx);
    multV0CCorr = AliESDUtils::GetCorrV0C(multV0C,zvtx);
    
    // Equalized signals // From AliCentralitySelectionTask // Updated
    for(Int_t iCh = 32; iCh < 64; ++iCh) {
        Double_t mult = event->GetVZEROEqMultiplicity(iCh);
        multV0AEq += mult;
    }
    for(Int_t iCh = 0; iCh < 32; ++iCh) {
        Double_t mult = event->GetVZEROEqMultiplicity(iCh);
        multV0CEq += mult;
    }
    
    if ( lMethod == "V0M" ) lreturnval =
        fBoundaryHisto_V0M -> GetBinContent( fBoundaryHisto_V0M->FindBin(multV0ACorr+multV0CCorr) );
    if ( lMethod == "V0A" ) lreturnval =
        fBoundaryHisto_V0A -> GetBinContent( fBoundaryHisto_V0A->FindBin(multV0ACorr) );
    if ( lMethod == "V0C" ) lreturnval =
        fBoundaryHisto_V0C -> GetBinContent( fBoundaryHisto_V0C->FindBin(multV0CCorr) );
    //equalized
    if ( lMethod == "V0MEq" ) lreturnval =
        fBoundaryHisto_V0MEq -> GetBinContent( fBoundaryHisto_V0MEq->FindBin(multV0AEq+multV0CEq) );
    if ( lMethod == "V0AEq" ) lreturnval =
        fBoundaryHisto_V0AEq -> GetBinContent( fBoundaryHisto_V0AEq->FindBin(multV0AEq) );
    if ( lMethod == "V0CEq" ) lreturnval =
        fBoundaryHisto_V0CEq -> GetBinContent( fBoundaryHisto_V0CEq->FindBin(multV0CEq) );
    
    return lreturnval;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::LoadCalibration(Int_t lLoadThisCalibration)
//To be called if starting analysis on a new run
{
    //If Histograms exist, de-allocate to prevent memory leakage
    if( fBoundaryHisto_V0M ) {fBoundaryHisto_V0M->Delete(); fBoundaryHisto_V0M = 0x0; }
    if( fBoundaryHisto_V0A ) {fBoundaryHisto_V0A->Delete(); fBoundaryHisto_V0A = 0x0; }
    if( fBoundaryHisto_V0C ) {fBoundaryHisto_V0C->Delete(); fBoundaryHisto_V0C = 0x0; }
    if( fBoundaryHisto_V0MEq ) {fBoundaryHisto_V0MEq->Delete(); fBoundaryHisto_V0MEq = 0x0; }
    if( fBoundaryHisto_V0AEq ) {fBoundaryHisto_V0AEq->Delete(); fBoundaryHisto_V0AEq = 0x0; }
    if( fBoundaryHisto_V0CEq ) {fBoundaryHisto_V0CEq->Delete(); fBoundaryHisto_V0CEq = 0x0; }
    
    AliInfo(Form( "Loading calibration file for run %i",lLoadThisCalibration) );
    TFile *lCalibFile_V0M = 0x0;
    TFile *lCalibFile_V0A = 0x0;
    TFile *lCalibFile_V0C = 0x0;
    TFile *lCalibFile_V0MEq = 0x0;
    TFile *lCalibFile_V0AEq = 0x0;
    TFile *lCalibFile_V0CEq = 0x0;
    
    //AliInfo("Calling TFile::Open");
    lCalibFile_V0M = TFile::Open("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/corrections/calibration_V0M.root");
    lCalibFile_V0A = TFile::Open("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/corrections/calibration_V0A.root");
    lCalibFile_V0C = TFile::Open("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/corrections/calibration_V0C.root");
    lCalibFile_V0MEq = TFile::Open("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/corrections/calibration_V0MEq.root");
    lCalibFile_V0AEq = TFile::Open("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/corrections/calibration_V0AEq.root");
    lCalibFile_V0CEq = TFile::Open("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/corrections/calibration_V0CEq.root");

    //AliInfo("Casting");
    //check memory consumption later... 
    fBoundaryHisto_V0M   = dynamic_cast<TH1F *>(lCalibFile_V0M  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0A   = dynamic_cast<TH1F *>(lCalibFile_V0A  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0C   = dynamic_cast<TH1F *>(lCalibFile_V0C  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0MEq   = dynamic_cast<TH1F *>(lCalibFile_V0MEq  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0AEq   = dynamic_cast<TH1F *>(lCalibFile_V0AEq  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0CEq   = dynamic_cast<TH1F *>(lCalibFile_V0CEq  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    
    if ( !fBoundaryHisto_V0M   || !fBoundaryHisto_V0A   || !fBoundaryHisto_V0C ||
        !fBoundaryHisto_V0MEq || !fBoundaryHisto_V0AEq || !fBoundaryHisto_V0CEq ){
        AliInfo(Form("No calibration for run %i exists at the moment!",lLoadThisCalibration));
        return kFALSE; //return denial
    }
    
    //Careful with manual cleanup if needed: to be implemented
    fBoundaryHisto_V0M->SetDirectory(0);
    fBoundaryHisto_V0A->SetDirectory(0);
    fBoundaryHisto_V0C->SetDirectory(0);
    fBoundaryHisto_V0MEq->SetDirectory(0);
    fBoundaryHisto_V0AEq->SetDirectory(0);
    fBoundaryHisto_V0CEq->SetDirectory(0);
    
    //AliInfo("Closing");
    
    if( lCalibFile_V0M ) lCalibFile_V0M->Close();
    if( lCalibFile_V0A ) lCalibFile_V0A->Close();
    if( lCalibFile_V0C ) lCalibFile_V0C->Close();
    if( lCalibFile_V0MEq ) lCalibFile_V0MEq->Close();
    if( lCalibFile_V0AEq ) lCalibFile_V0AEq->Close();
    if( lCalibFile_V0CEq ) lCalibFile_V0CEq->Close();
    
    fRunNumber = lLoadThisCalibration; //Loaded!
    AliInfo(Form("Finished loading calibration for run %i",lLoadThisCalibration));
    return kTRUE;
}
