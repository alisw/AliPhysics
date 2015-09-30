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
#include "AliESDtrackCuts.h"
#include "AliPPVsMultUtils.h"
#include <TFile.h>
#include "AliAODHeader.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"


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
    fBoundaryHisto_V0CEq(0),
    fBoundaryHisto_V0B(0),
    fBoundaryHisto_V0Apartial(0),
    fBoundaryHisto_V0Cpartial(0),
    fBoundaryHisto_V0S(0),
    fBoundaryHisto_V0SB(0),
    fAverageAmplitudes(0)
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
    fBoundaryHisto_V0CEq(0),
    fBoundaryHisto_V0B(0),
    fBoundaryHisto_V0Apartial(0),
    fBoundaryHisto_V0Cpartial(0),
    fBoundaryHisto_V0S(0),
    fBoundaryHisto_V0SB(0),
    fAverageAmplitudes(0)
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

//_____________________________________________________________________________
Float_t AliPPVsMultUtils::MinVal( Float_t A, Float_t B ) {
    if( A < B ) {
        return A;
    }
    else {
        return B;
    }
}

//______________________________________________________________________
Float_t AliPPVsMultUtils::GetMultiplicityPercentile(AliVEvent *event, TString lMethod, Bool_t lEmbedEventSelection)
// Function to get multiplicity quantiles
//
// Estimators available (use as strings in "lMethod"), e.g.
//  ::GetMultiplicityPercentile( [event object] , "V0M")
//
//  --- V0M: Sum of amplitudes in V0A and V0C
//  --- V0A: VZERO amplitudes (A side)
//  --- V0C: VZERO amplitudes (C side)
//  --- V0MEq: Sum of amplitudes in V0A and V0C, equalized (experimental)
//  --- V0AEq: VZERO amplitudes (A side), equalized (experimental)
//  --- V0CEq: VZERO amplitudes (C side)
//  --- V0B: Simultaneous selection in V0A and V0C
//           ( implemented via (V0A/<V0A>) > x and (V0C/<V0C>) > x )
//  --- V0Apartial: 2 rings selected such that  2.8 < eta <  3.9
//  --- V0Cpartial: 2 rings selected such that -3.7 < eta < -2.7
//  --- V0S: Symmetrized selection in V0A and V0C
//           ( implemented via (V0Apartial/<V0Apartial> + V0Cpartial/<V0Cpartial>) > x )
//  --- V0SB: Symmetrized simultaneous selection in V0A and V0C
//           ( implemented via (V0Apartial/<V0A>) > x and (V0Cpartial/<V0Cpartial>) > x )
//
//  This getter automatically includes event selection by default and will return negative
//  values for the following types of events:
//
//  --- Events that don't have at least one tracklet
//  --- Events without reconstructed SPD vertex
//  --- Events with a PV falling outside |z|<10cm
//  --- Events that are tagged as pileup with IsPileupFromSPDInMultBins
//
//  For more info, please consult AliESDtrackCuts::GetReferenceMultiplicity
//  (and, more specifically, the kTracklets option)
{
    Int_t lRequestedRunNumber = event->GetRunNumber();
    if( lRequestedRunNumber != fRunNumber ) fCalibrationLoaded = LoadCalibration( lRequestedRunNumber );

    if ( !fCalibrationLoaded ) {
        return -1000; //Return absurd value (hopefully noone will use this...)
    }
    if ( !fBoundaryHisto_V0M   || !fBoundaryHisto_V0A   || !fBoundaryHisto_V0C ||
            !fBoundaryHisto_V0MEq || !fBoundaryHisto_V0AEq || !fBoundaryHisto_V0CEq || !fBoundaryHisto_V0B || !fBoundaryHisto_V0Apartial || !fBoundaryHisto_V0Cpartial ||
            !fBoundaryHisto_V0S || !fBoundaryHisto_V0SB || !fAverageAmplitudes ) {
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
    Float_t multV0Apartial = 0;
    Float_t multV0Cpartial = 0;


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

    for(Int_t iCh = 32; iCh < 48; iCh++) {
        Double_t mult = esdV0->GetMultiplicity(iCh);
        multV0Apartial += mult;
    }
    for(Int_t iCh = 0; iCh < 16; iCh++) {
        Double_t mult = esdV0->GetMultiplicity(iCh);
        multV0Cpartial += mult;
    }

    if ( lMethod == "V0M" ) lreturnval =
            fBoundaryHisto_V0M -> GetBinContent( fBoundaryHisto_V0M->FindBin(multV0A+multV0C) );
    if ( lMethod == "V0A" ) lreturnval =
            fBoundaryHisto_V0A -> GetBinContent( fBoundaryHisto_V0A->FindBin(multV0A) );
    if ( lMethod == "V0C" ) lreturnval =
            fBoundaryHisto_V0C -> GetBinContent( fBoundaryHisto_V0C->FindBin(multV0C) );
    //equalized
    if ( lMethod == "V0MEq" ) lreturnval =
            fBoundaryHisto_V0MEq -> GetBinContent( fBoundaryHisto_V0MEq->FindBin(multV0AEq+multV0CEq) );
    if ( lMethod == "V0AEq" ) lreturnval =
            fBoundaryHisto_V0AEq -> GetBinContent( fBoundaryHisto_V0AEq->FindBin(multV0AEq) );
    if ( lMethod == "V0CEq" ) lreturnval =
            fBoundaryHisto_V0CEq -> GetBinContent( fBoundaryHisto_V0CEq->FindBin(multV0CEq) );
    //extra stuff
    if ( lMethod == "V0B" ) lreturnval =
            fBoundaryHisto_V0B -> GetBinContent( fBoundaryHisto_V0B->FindBin( MinVal( multV0A / fAverageAmplitudes->GetBinContent(1) , multV0C / fAverageAmplitudes->GetBinContent(2) ) ) );
    if ( lMethod == "V0Apartial" ) lreturnval =
            fBoundaryHisto_V0Apartial -> GetBinContent( fBoundaryHisto_V0Apartial->FindBin( multV0Apartial ) );
    if ( lMethod == "V0Cpartial" ) lreturnval =
            fBoundaryHisto_V0Cpartial -> GetBinContent( fBoundaryHisto_V0Cpartial->FindBin( multV0Cpartial ) );
    if ( lMethod == "V0S" ) lreturnval =
            fBoundaryHisto_V0S -> GetBinContent( fBoundaryHisto_V0S->FindBin( (multV0Apartial/fAverageAmplitudes->GetBinContent(3)) + (multV0Cpartial/fAverageAmplitudes->GetBinContent(4)) ) );
    if ( lMethod == "V0SB" ) lreturnval =
            fBoundaryHisto_V0SB -> GetBinContent( fBoundaryHisto_V0SB->FindBin( MinVal( multV0Apartial / fAverageAmplitudes->GetBinContent(3) , multV0Cpartial / fAverageAmplitudes->GetBinContent(4) ) ) );

    if ( lEmbedEventSelection ) {
        if(IsSelectedTrigger                        ( event ) == kFALSE ) lreturnval = -200;
        if(IsINELgtZERO                         ( event ) == kFALSE ) lreturnval = -201;
        if(IsAcceptedVertexPosition             ( event ) == kFALSE ) lreturnval = -202;
        if(IsNotPileupSPDInMultBins             ( event ) == kFALSE ) lreturnval = -203;
        if(HasNoInconsistentSPDandTrackVertices ( event ) == kFALSE ) lreturnval = -204;
    }

    return lreturnval;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::LoadCalibration(Int_t lLoadThisCalibration)
//To be called if starting analysis on a new run
{
    //If Histograms exist, de-allocate to prevent memory leakage
    if( fBoundaryHisto_V0M ) {
        fBoundaryHisto_V0M->Delete();
        fBoundaryHisto_V0M = 0x0;
    }
    if( fBoundaryHisto_V0A ) {
        fBoundaryHisto_V0A->Delete();
        fBoundaryHisto_V0A = 0x0;
    }
    if( fBoundaryHisto_V0C ) {
        fBoundaryHisto_V0C->Delete();
        fBoundaryHisto_V0C = 0x0;
    }
    if( fBoundaryHisto_V0MEq ) {
        fBoundaryHisto_V0MEq->Delete();
        fBoundaryHisto_V0MEq = 0x0;
    }
    if( fBoundaryHisto_V0AEq ) {
        fBoundaryHisto_V0AEq->Delete();
        fBoundaryHisto_V0AEq = 0x0;
    }
    if( fBoundaryHisto_V0CEq ) {
        fBoundaryHisto_V0CEq->Delete();
        fBoundaryHisto_V0CEq = 0x0;
    }
    if( fBoundaryHisto_V0B ) {
        fBoundaryHisto_V0B->Delete();
        fBoundaryHisto_V0B = 0x0;
    }
    if( fBoundaryHisto_V0Apartial ) {
        fBoundaryHisto_V0Apartial->Delete();
        fBoundaryHisto_V0Apartial = 0x0;
    }
    if( fBoundaryHisto_V0Cpartial ) {
        fBoundaryHisto_V0Cpartial->Delete();
        fBoundaryHisto_V0Cpartial = 0x0;
    }
    if( fBoundaryHisto_V0S ) {
        fBoundaryHisto_V0S->Delete();
        fBoundaryHisto_V0S = 0x0;
    }
    if( fBoundaryHisto_V0SB ) {
        fBoundaryHisto_V0SB->Delete();
        fBoundaryHisto_V0SB = 0x0;
    }
    if( fAverageAmplitudes ) {
        fAverageAmplitudes->Delete();
        fAverageAmplitudes = 0x0;
    }

    AliInfo(Form( "Loading calibration file for run %i",lLoadThisCalibration) );
    TFile *lCalibFile_V0M = 0x0;
    TFile *lCalibFile_V0A = 0x0;
    TFile *lCalibFile_V0C = 0x0;
    TFile *lCalibFile_V0MEq = 0x0;
    TFile *lCalibFile_V0AEq = 0x0;
    TFile *lCalibFile_V0CEq = 0x0;
    TFile *lCalibFile_V0B = 0x0;
    TFile *lCalibFile_V0Apartial = 0x0;
    TFile *lCalibFile_V0Cpartial = 0x0;
    TFile *lCalibFile_V0S = 0x0;
    TFile *lCalibFile_V0SB = 0x0;
    TFile *lCalibFile_Averages = 0x0;

    //AliInfo("Calling TFile::Open");
    lCalibFile_V0M = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0M.root");
    lCalibFile_V0A = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0A.root");
    lCalibFile_V0C = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0C.root");
    lCalibFile_V0MEq = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0MEq.root");
    lCalibFile_V0AEq = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0AEq.root");
    lCalibFile_V0CEq = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0CEq.root");
    lCalibFile_V0B = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0B.root");
    lCalibFile_V0Apartial = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0Apartial.root");
    lCalibFile_V0Cpartial = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0Cpartial.root");
    lCalibFile_V0S = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0S.root");
    lCalibFile_V0SB = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calibration_adaptive_V0SB.root");
    lCalibFile_Averages = TFile::Open("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/corrections/calib-averages.root");

    //AliInfo("Casting");
    fBoundaryHisto_V0M        = dynamic_cast<TH1F *>(lCalibFile_V0M  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0A        = dynamic_cast<TH1F *>(lCalibFile_V0A  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0C        = dynamic_cast<TH1F *>(lCalibFile_V0C  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0MEq      = dynamic_cast<TH1F *>(lCalibFile_V0MEq  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0AEq      = dynamic_cast<TH1F *>(lCalibFile_V0AEq  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0CEq      = dynamic_cast<TH1F *>(lCalibFile_V0CEq  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0B        = dynamic_cast<TH1F *>(lCalibFile_V0B  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0Apartial = dynamic_cast<TH1F *>(lCalibFile_V0Apartial  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0Cpartial = dynamic_cast<TH1F *>(lCalibFile_V0Cpartial  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0S        = dynamic_cast<TH1F *>(lCalibFile_V0S  -> Get(Form("histocalib%i",lLoadThisCalibration)) );
    fBoundaryHisto_V0SB       = dynamic_cast<TH1F *>(lCalibFile_V0SB  -> Get(Form("histocalib%i",lLoadThisCalibration)) );

    //Average Amplitudes for weighting
    fAverageAmplitudes       = dynamic_cast<TH1D *>(lCalibFile_Averages  -> Get(Form("hcalib_averages_%i",lLoadThisCalibration)) );

    if ( !fBoundaryHisto_V0M   || !fBoundaryHisto_V0A   || !fBoundaryHisto_V0C ||
            !fBoundaryHisto_V0MEq || !fBoundaryHisto_V0AEq || !fBoundaryHisto_V0CEq || !fBoundaryHisto_V0B || !fBoundaryHisto_V0Apartial || !fBoundaryHisto_V0Cpartial ||
            !fBoundaryHisto_V0S || !fBoundaryHisto_V0SB || !fAverageAmplitudes ) {
        AliInfo(Form("No calibration for run %i exists at the moment!",lLoadThisCalibration));
        if( lCalibFile_V0M ) {
            lCalibFile_V0M->Close();
            lCalibFile_V0M->Delete();
            delete lCalibFile_V0M;
        }
        if( lCalibFile_V0A ) {
            lCalibFile_V0A->Close();
            lCalibFile_V0A->Delete();
            delete lCalibFile_V0A;
        }
        if( lCalibFile_V0C ) {
            lCalibFile_V0C->Close();
            lCalibFile_V0C->Delete();
            delete lCalibFile_V0C;
        }
        if( lCalibFile_V0MEq ) {
            lCalibFile_V0MEq->Close();
            lCalibFile_V0MEq->Delete();
            delete lCalibFile_V0MEq;
        }
        if( lCalibFile_V0AEq ) {
            lCalibFile_V0AEq->Close();
            lCalibFile_V0AEq->Delete();
            delete lCalibFile_V0AEq;
        }
        if( lCalibFile_V0CEq ) {
            lCalibFile_V0CEq->Close();
            lCalibFile_V0CEq->Delete();
            delete lCalibFile_V0CEq;
        }
        if( lCalibFile_V0B ) {
            lCalibFile_V0B->Close();
            lCalibFile_V0B->Delete();
            delete lCalibFile_V0B;
        }
        if( lCalibFile_V0Apartial ) {
            lCalibFile_V0Apartial->Close();
            lCalibFile_V0Apartial->Delete();
            delete lCalibFile_V0Apartial;
        }
        if( lCalibFile_V0Cpartial ) {
            lCalibFile_V0Cpartial->Close();
            lCalibFile_V0Cpartial->Delete();
            delete lCalibFile_V0Cpartial;
        }
        if( lCalibFile_V0S ) {
            lCalibFile_V0S->Close();
            lCalibFile_V0S->Delete();
            delete lCalibFile_V0S;
        }
        if( lCalibFile_V0SB ) {
            lCalibFile_V0SB->Close();
            lCalibFile_V0SB->Delete();
            delete lCalibFile_V0SB;
        }
        if( lCalibFile_Averages ) {
            lCalibFile_Averages->Close();
            lCalibFile_Averages->Delete();
            delete lCalibFile_Averages;
        }
        fRunNumber = lLoadThisCalibration;
        return kFALSE; //return denial
    }

    fBoundaryHisto_V0M->SetName("fBoundaryHisto_V0M");
    fBoundaryHisto_V0A->SetName("fBoundaryHisto_V0A");
    fBoundaryHisto_V0C->SetName("fBoundaryHisto_V0C");
    fBoundaryHisto_V0MEq->SetName("fBoundaryHisto_V0MEq");
    fBoundaryHisto_V0AEq->SetName("fBoundaryHisto_V0AEq");
    fBoundaryHisto_V0CEq->SetName("fBoundaryHisto_V0CEq");
    fBoundaryHisto_V0B->SetName("fBoundaryHisto_V0B");
    fBoundaryHisto_V0Apartial->SetName("fBoundaryHisto_V0Apartial");
    fBoundaryHisto_V0Cpartial->SetName("fBoundaryHisto_V0Cpartial");
    fBoundaryHisto_V0S->SetName("fBoundaryHisto_V0S");
    fBoundaryHisto_V0SB->SetName("fBoundaryHisto_V0SB");
    fAverageAmplitudes->SetName("fBoundaryHisto_V0SB");

    //Careful with manual cleanup if needed: to be implemented
    fBoundaryHisto_V0M->SetDirectory(0);
    fBoundaryHisto_V0A->SetDirectory(0);
    fBoundaryHisto_V0C->SetDirectory(0);
    fBoundaryHisto_V0MEq->SetDirectory(0);
    fBoundaryHisto_V0AEq->SetDirectory(0);
    fBoundaryHisto_V0CEq->SetDirectory(0);
    fBoundaryHisto_V0B->SetDirectory(0);
    fBoundaryHisto_V0Apartial->SetDirectory(0);
    fBoundaryHisto_V0Cpartial->SetDirectory(0);
    fBoundaryHisto_V0S->SetDirectory(0);
    fBoundaryHisto_V0SB->SetDirectory(0);
    fAverageAmplitudes->SetDirectory(0);

    //AliInfo("Closing");
    if( lCalibFile_V0M ) {
        lCalibFile_V0M->Close();
        lCalibFile_V0M->Delete();
        delete lCalibFile_V0M;
    }
    if( lCalibFile_V0A ) {
        lCalibFile_V0A->Close();
        lCalibFile_V0A->Delete();
        delete lCalibFile_V0A;
    }
    if( lCalibFile_V0C ) {
        lCalibFile_V0C->Close();
        lCalibFile_V0C->Delete();
        delete lCalibFile_V0C;
    }
    if( lCalibFile_V0MEq ) {
        lCalibFile_V0MEq->Close();
        lCalibFile_V0MEq->Delete();
        delete lCalibFile_V0MEq;
    }
    if( lCalibFile_V0AEq ) {
        lCalibFile_V0AEq->Close();
        lCalibFile_V0AEq->Delete();
        delete lCalibFile_V0AEq;
    }
    if( lCalibFile_V0CEq ) {
        lCalibFile_V0CEq->Close();
        lCalibFile_V0CEq->Delete();
        delete lCalibFile_V0CEq;
    }
    if( lCalibFile_V0B ) {
        lCalibFile_V0B->Close();
        lCalibFile_V0B->Delete();
        delete lCalibFile_V0B;
    }
    if( lCalibFile_V0Apartial ) {
        lCalibFile_V0Apartial->Close();
        lCalibFile_V0Apartial->Delete();
        delete lCalibFile_V0Apartial;
    }
    if( lCalibFile_V0Cpartial ) {
        lCalibFile_V0Cpartial->Close();
        lCalibFile_V0Cpartial->Delete();
        delete lCalibFile_V0Cpartial;
    }
    if( lCalibFile_V0S ) {
        lCalibFile_V0S->Close();
        lCalibFile_V0S->Delete();
        delete lCalibFile_V0S;
    }
    if( lCalibFile_V0SB ) {
        lCalibFile_V0SB->Close();
        lCalibFile_V0SB->Delete();
        delete lCalibFile_V0SB;
    }
    if( lCalibFile_Averages ) {
        lCalibFile_Averages->Close();
        lCalibFile_Averages->Delete();
        delete lCalibFile_Averages;
    }

    fRunNumber = lLoadThisCalibration; //Loaded!
    AliInfo(Form("Finished loading calibration for run %i",lLoadThisCalibration));
    return kTRUE;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::IsMinimumBias(AliVEvent* event)
// Function to check for minimum-bias trigger (AliVEvent::kMB)
{
    //Code to reject events that aren't kMB
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
    return isSelected;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::IsSelectedTrigger(AliVEvent* event, AliVEvent::EOfflineTriggerTypes trigType)
// Function to check for a specific trigger class available in AliVEvent (default AliVEvent::kMB)
{
    //Code to reject events that aren't trigType
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & trigType) == trigType;
    return isSelected;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::IsSelectedTrigger(AliVEvent* event, TString trigName)
// Function to check for a specific trigger class with name "trigName"
{
    //Code to reject events that aren't kMB
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    TString classTrigg = event->GetFiredTriggerClasses(); 
    isSelected = maskIsSelected && classTrigg.Contains(trigName.Data());
    return isSelected;
}


//______________________________________________________________________
Bool_t AliPPVsMultUtils::IsINELgtZERO(AliVEvent *event)
// Function to check for INEL > 0 condition
// Makes use of tracklets and requires at least and SPD vertex
{
    Bool_t lReturnValue = kFALSE;
    //Use Ref.Mult. code...
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        if ( AliESDtrackCuts::GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets, 1.0) >= 1 ) lReturnValue = kTRUE;
    }
    //Redo equivalent test
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;

        //FIXME --- Actually, here we can come up with a workaround.
        // We can check for the reference multiplicity stored and look for error codes!

        AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
        Int_t lStoredRefMult = header->GetRefMultiplicityComb08();

        //Get Multiplicity object
        AliAODTracklets *spdmult = aodevent->GetMultiplicity();
        for (Int_t i=0; i<spdmult->GetNumberOfTracklets(); ++i)
        {
            if ( lStoredRefMult != -1 && lStoredRefMult != -2 && TMath::Abs(spdmult->GetEta(i)) < 1.0 ) lReturnValue = kTRUE;
        }
    }
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::IsAcceptedVertexPosition(AliVEvent *event)
// Simple check for the best primary vertex Z position:
// Will accept events only if |z| < 10cm
{
    Bool_t lReturnValue = kFALSE;
    //Getting around to the best vertex -> typecast to ESD/AOD
    const AliVVertex *lPrimaryVtx = NULL;
    /* get ESD vertex */
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        lPrimaryVtx = esdevent->GetPrimaryVertex();
    }
    /* get AOD vertex */
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;
        lPrimaryVtx = aodevent->GetPrimaryVertex();
    }
    if ( TMath::Abs( lPrimaryVtx->GetZ() ) <= 10.0 ) lReturnValue = kTRUE;
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices(AliVEvent *event)
// This function checks if track and SPD vertices are consistent.
// N.B.: It is rigorously a "Not Inconsistent" function which will
// let events with only SPD vertex go through without troubles.
{
    //It's consistent until proven otherwise...
    Bool_t lReturnValue = kTRUE;

    //Getting around to the best vertex -> typecast to ESD/AOD


    /* get ESD vertex */
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        const AliESDVertex *lPrimaryVtxSPD    = NULL;
        const AliESDVertex *lPrimaryVtxTracks = NULL;
        
        lPrimaryVtxSPD    = esdevent->GetPrimaryVertexSPD   ();
        lPrimaryVtxTracks = esdevent->GetPrimaryVertexTracks();

        //Only continue if track vertex defined
        if( lPrimaryVtxTracks->GetStatus() && lPrimaryVtxSPD->GetStatus() ){
            //Copy-paste from refmult estimator
            // TODO value of displacement to be studied
            const Float_t maxDisplacement = 0.5;
            //check for displaced vertices
            Double_t displacement = TMath::Abs(lPrimaryVtxSPD->GetZ() - lPrimaryVtxTracks->GetZ());
            if (displacement > maxDisplacement) lReturnValue = kFALSE;
        }
    }
    /* get AOD vertex */
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;

        //FIXME - Hack to deal with the fact that no
        //        AliAODEvent::GetPrimaryVertexTracks() exists...
        AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
        Int_t lStoredRefMult = header->GetRefMultiplicityComb08();
        if( lStoredRefMult == -4 ) lReturnValue = kFALSE;
    }
    return lReturnValue;
}

//______________________________________________________________________
Int_t AliPPVsMultUtils::GetStandardReferenceMultiplicity(AliVEvent *event, Bool_t lEmbedEventSelection)
// Wrapper for AliESDtrackCuts::GetReferenceMultiplicity
//
// By default, uses kTrackletsITSTPC if tracking vertex exists
// but will fall back to kTracklets if only SPD vertex exists
{
    //It's consistent until proven otherwise...
    Long_t lReturnValue = -10; //Kill this event, please
    
    /* get ESD vertex */
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        //It should always be this easy...
        const AliESDVertex *lPrimaryVtxTracks = esdevent->GetPrimaryVertexTracks();
        if (!lPrimaryVtxTracks->GetStatus()) {
            //If no track vertex available, fall back on kTracklets
            lReturnValue = AliESDtrackCuts::GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets, 0.8);
        } else {
            //If track vertex available, use combined estimator
            lReturnValue = AliESDtrackCuts::GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8);
        }
    }
    /* get AOD vertex */
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;

        AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
        Long_t lStoredRefMult = header->GetRefMultiplicityComb08();
        // (-1) -> kept, user might discard, but will be killed anyhow by ev. sel.
        // (-2) -> kept, user might discard, but will be killed anyhow by ev. sel.
        // (-3) -> only if old AOD, but then -> recompute
        // (-4) -> kept, user might discard, but will be killed anyhow by ev. sel.

        //Hack: if this is -3, this is an old AOD filtering and we need to fall back on tracklets.
        //Let's count them...
        if( lStoredRefMult == -3 ) { //(then -1 and -2 are NOT the case, -4 unchecked but impossible since no track vtx)
            //Get Multiplicity object
            AliAODTracklets *spdmult = aodevent->GetMultiplicity();
            Long_t lNTracklets = 0;
            for (Int_t i=0; i<spdmult->GetNumberOfTracklets(); ++i)
            {
                if ( TMath::Abs(spdmult->GetEta(i)) < 0.8 ) lNTracklets++;
            }
            lReturnValue = lNTracklets;
        } else {
            lReturnValue = lStoredRefMult;
        }
    }

    if ( lEmbedEventSelection ) {
        if(IsSelectedTrigger                        ( event ) == kFALSE ) lReturnValue = -200;
        if(IsINELgtZERO                         ( event ) == kFALSE ) lReturnValue = -201;
        if(IsAcceptedVertexPosition             ( event ) == kFALSE ) lReturnValue = -202;
        if(IsNotPileupSPDInMultBins             ( event ) == kFALSE ) lReturnValue = -203;
        if(HasNoInconsistentSPDandTrackVertices ( event ) == kFALSE ) lReturnValue = -204;
    }
    
    return lReturnValue;
}


//______________________________________________________________________
Bool_t AliPPVsMultUtils::IsNotPileupSPDInMultBins(AliVEvent *event)
// Checks if not pileup from SPD (via IsPileupFromSPDInMultBins)
{
    Bool_t lReturnValue = kTRUE;
    //Getting around to the SPD vertex -> typecast to ESD/AOD
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        if ( esdevent->IsPileupFromSPDInMultBins() == kTRUE ) lReturnValue = kFALSE;
    }
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;
        if ( aodevent->IsPileupFromSPDInMultBins() == kTRUE ) lReturnValue = kFALSE;
    }
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::IsEventSelected(AliVEvent *event, AliVEvent::EOfflineTriggerTypes trigType)
// Single function which does a series of selections:
//  1) kMB trigger selection
//  2) INEL > 0 (with tracklets)
//  3) Checks for accepted vertex position (|eta|<10cm)
//  4) Checks for consistent SPD and track vertex (if track vertex exists)
//  5) Rejects events tagged with IsPileupFromSPDInMultBins()
{
    Bool_t lReturnValue = kFALSE;
    if ( IsNotPileupSPDInMultBins               ( event ) == kTRUE &&
            IsINELgtZERO                        ( event ) == kTRUE &&
            IsAcceptedVertexPosition            ( event ) == kTRUE &&
            HasNoInconsistentSPDandTrackVertices( event ) == kTRUE &&
            IsSelectedTrigger                       ( event , trigType) == kTRUE
       ) lReturnValue = kTRUE;
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliPPVsMultUtils::IsEventSelected(AliVEvent *event, TString trigName)
// Single function which does a series of selections:
//  1) kMB trigger selection
//  2) INEL > 0 (with tracklets)
//  3) Checks for accepted vertex position (|eta|<10cm)
//  4) Checks for consistent SPD and track vertex (if track vertex exists)
//  5) Rejects events tagged with IsPileupFromSPDInMultBins()
{
    Bool_t lReturnValue = kFALSE;
    if ( IsNotPileupSPDInMultBins               ( event ) == kTRUE &&
            IsINELgtZERO                        ( event ) == kTRUE &&
            IsAcceptedVertexPosition            ( event ) == kTRUE &&
            HasNoInconsistentSPDandTrackVertices( event ) == kTRUE &&
            IsSelectedTrigger                       ( event , trigName) == kTRUE
       ) lReturnValue = kTRUE;
    return lReturnValue;
}

