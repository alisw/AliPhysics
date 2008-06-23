/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//                     Class AliRsnPairSimple
//-------------------------------------------------------------------------
// This class computes the invariant mass spectrum of a specified pair of
// particles, throughout a list of AliRsnEvents, and returns it as a TH1D.
// This object is not supposed to be used directly: an AliRsnAnalysis
// should be initialized in a macro and filled with one or more AliRsnPairSimple's
// which are then processed with a given sample of events.
//
// author: A. Pulvirenti
// email : alberto.pulvirenti@ct.infn.it
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TH1.h>
#include <TString.h>
#include <TArrayI.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRsnMCInfo.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnCut.h"
#include "AliRsnCutMgr.h"
#include "AliRsnPairDef.h"
#include "AliRsnPID.h"

#include "AliRsnPairSimple.h"

ClassImp(AliRsnPairSimple)

//--------------------------------------------------------------------------------------------------------
AliRsnPairSimple::AliRsnPairSimple(AliRsnPairDef *pd, const char *name, const char *title) :
  TNamed(name, title),
  fPIDMethod(AliRsnPID::kNone),
  fForMixing(kFALSE),
  fStoreOnlyTrue(kFALSE),
  fCuts(0x0),
  fPair(),
  fPairDef(pd),
  fHistogram(0x0),
  fHistogramMC(0x0)
{
//
// Constructor.
// This constructor allows to define some of the initialization values:
//  - name and title of the object
//  - histogram binning and edges
// The other parameters are initialized as in the default constructor.
//
}
//--------------------------------------------------------------------------------------------------------
AliRsnPairSimple::AliRsnPairSimple(const AliRsnPairSimple &copy) :
  TNamed(copy),
  fPIDMethod(copy.fPIDMethod),
  fForMixing(copy.fForMixing),
  fStoreOnlyTrue(copy.fStoreOnlyTrue),
  fCuts(0x0),
  fPair(),
  fPairDef(copy.fPairDef),
  fHistogram(0x0),
  fHistogramMC(0x0)
{
//
// Copy constructor.
// Default behavior as a copy constructor for what concerns non-array data-members.
// The arrays are cloned if they are not NULL.
//
}
//--------------------------------------------------------------------------------------------------------
const AliRsnPairSimple& AliRsnPairSimple::operator=(const AliRsnPairSimple &copy)
{
//
// Assignment operator.
// Default behavior like copy constructor.
//
    SetName(copy.GetName());
    SetTitle(copy.GetTitle());

    fHistogram = 0x0;
    fHistogramMC = 0x0;
    fPairDef = copy.fPairDef;

    fPIDMethod = copy.fPIDMethod;
    fForMixing = copy.fForMixing;
    fStoreOnlyTrue = copy.fStoreOnlyTrue;
    if (copy.fHistogram) fHistogram = (TH1D*)(copy.fHistogram->Clone());

    fCuts = 0x0;

    return (*this);
}
//--------------------------------------------------------------------------------------------------------
void AliRsnPairSimple::Clear(Option_t* /*option*/)
{
//
// Clear arrays and histogram.
// For the sake of security, all pointers are also set explicitly to NULL.
//
    delete fHistogram;
    fHistogram = 0x0;
    delete fHistogramMC;
    fHistogramMC = 0x0;
}

//--------------------------------------------------------------------------------------------------------
void AliRsnPairSimple::InitHistograms()
{
//
// Initialize histograms
//

    Int_t nbins = fPairDef->GetNBins();
    Double_t min = fPairDef->GetMin(), max = fPairDef->GetMax();

    if (strlen(GetName()) > 0) {
        fHistogram = new TH1D(GetName(), GetTitle(), nbins, min, max);
        fHistogramMC = new TH1D(Form("MC_%s", GetName()), Form("%s (MC)", GetTitle()), nbins, min, max);
    }
    else {
        Char_t name[200], title[200];
        strcpy(name, GetHistName());
        strcpy(title, GetHistTitle());
        fHistogram = new TH1D(name, title, nbins, min, max);
        fHistogramMC = new TH1D(Form("MC_%s", name), Form("%s (MC)", title), nbins, min, max);
    }
    fHistogram->Sumw2();
    fHistogramMC->Sumw2();
}

const char* AliRsnPairSimple::GetHistName ()
{
//
// Creates the histogram name, given a cut manager
//

    TString strName("");
    strName.Append(AliRsnPID::ParticleName(fPairDef->GetType(0)));
    strName += '(';
    strName += fPairDef->GetCharge(0);
    strName += ')';
    strName += '_';
    strName.Append(AliRsnPID::ParticleName(fPairDef->GetType(1)));
    strName += '(';
    strName += fPairDef->GetCharge(1);
    strName += ')';
    strName += '_';
    if (fCuts) {
        strName.Append("cuts:");
        strName.Append(fCuts->GetName());
    }
    else {
        strName.Append("NoCuts");
    }
    if (fStoreOnlyTrue) strName.Append("_true");
    if (fForMixing) strName.Append("_mix");

    return strName.Data();
}

const char* AliRsnPairSimple::GetHistTitle()
{
//
// Creates the histogram title, given a cut manager
//

    TString strName("Inv. mass of ");
    strName.Append(AliRsnPID::ParticleName(fPairDef->GetType(0), kFALSE));
    strName += 's';
    strName += ' ';
    strName += '(';
    strName += fPairDef->GetCharge(0);
    strName += ')';
    strName.Append(" and ");
    strName.Append(AliRsnPID::ParticleName(fPairDef->GetType(1), kFALSE));
    strName += 's';
    strName += ' ';
    strName += '(';
    strName += fPairDef->GetCharge(1);
    strName += ')';
    if (fCuts) {
        strName.Append(" [cuts: ");
        strName.Append(fCuts->GetTitle());
        strName.Append("] ");
    }
    else {
        strName.Append(" [No cuts] ");
    }
    if (fStoreOnlyTrue) strName.Append(" [true pairs]");
    if (fForMixing) strName.Append(" [event mixing]");

    return strName.Data();
}

//--------------------------------------------------------------------------------------------------------
Stat_t AliRsnPairSimple::Process(AliRsnEvent *event1, AliRsnEvent *event2)
{
//
// Scans the two events specified in argument to fill the histogram.
// This method essentially calls the AliRsnPairSimple::Fill() method one or many times.
// When the "noPID" argument is kFALSE, the analysis is done with identified particles
// and this causes the Fill() method to be called only once, for the two lists of
// identified particles of the two kinds specified in AliRsnPairSimple datamembers.
// When the "noPID" argument is kTRUE, the analysis is done with all collections
// of particles of the same sign as specified in the two arguments of the pair.
// ---
// Particles of type #1 are taken in 'event1', and particles of type #2 are taken in 'event2'.
// When doing single-event analysis (for resonance signal or like-sign background),
// the second argument can be simply skipped.
// When doing event mixing, the two arguments must be not null and different.
// If argument #1 is NULL, an error is raised, while if argument #2 is NULL, no error is raised,
// and 'event2' argument is set equal to 'event1' (= single event processing).
// ---
// Return value is the total number of pairs processed.
//

    if (!event1) {
        // argument #1 cannot be NULL
        AliError("Argument #1 cannot be NULL.");
        return 0.0;
    }
    if (!event2) {
        // if argument #2 is NULL, it is put equal to argument #1
        event2 = event1;
    }
    if (!fPairDef) {
        AliError("No pairdef defined");
        return 0.0;
    }
    if (!fHistogram) {
        AliError("Histograms not initialized");
        return 0.0;
    }

    // assign pointers to the list of indexes to be used
    TArrayI *listCharged1 = 0x0, *listCharged2 = 0x0;
    //if (pidMethod == AliRsnPID::kNone) {
        listCharged1 = event1->GetCharged(fPairDef->GetCharge(0));
        listCharged2 = event2->GetCharged(fPairDef->GetCharge(1));
    //}
    //else {
        //listCharged1 = event1->GetTracksArray(pidMethod, fPairDef->GetCharge(0), fPairDef->GetType(0));
        //listCharged2 = event2->GetTracksArray(pidMethod, fPairDef->GetCharge(1), fPairDef->GetType(1));
    //}
    if (!listCharged1 || !listCharged2) return 0.0;
    TArrayI &list1 = *listCharged1;
    TArrayI &list2 = *listCharged2;

    Int_t   i1, i2, start2;
    Stat_t  nPairs = 0;
    AliRsnDaughter *track1 = 0x0, *track2 = 0x0;

    // loop on particle of type 1 (in first event)
    for (i1 = 0; i1 < list1.GetSize(); i1++) {
        track1 = event1->GetTrack(list1[i1]);
        //if (!track1) continue;
        if (!fCuts->IsSelected(AliRsnCut::kParticle, track1)) continue;
        // loop on particle of type 2 (in second event)
        // in case we are building a like-sign histogram with particles
        // of the same type in the same event, we must avoid that
        // each pair is computed twice
        start2 = 0;
        if (listCharged1 == listCharged2) start2 = i1 + 1;
        for (i2 = start2; i2 < list2.GetSize(); i2++) {
            track2 = event2->GetTrack(list2[i2]);
            //if (!track2) continue;
            if (!fCuts->IsSelected(AliRsnCut::kParticle, track2)) continue;
            nPairs += Process(track1, track2);
        }
    }

    return nPairs;
}

//_____________________________________________________________________________
Stat_t AliRsnPairSimple::Process
(AliRsnDaughter *track1, AliRsnDaughter *track2)
{
//
// Checks the single tracks and the pair against track cuts and,
// if the cuts are passed, fill the histograms with a weight
// given by the product of the appropriate PID probabilities of both.
// The method returns a boolean success value for eventually counting.
//

    // setup pair and check pair cuts
    fPair.SetPair(track1, track2);
    if (!fCuts->IsSelected(AliRsnCut::kPair, &fPair)) return 0.0;

    // if there is a request to process only the pairs born from a true resonance
    // this is checked here
    if (fStoreOnlyTrue && !fPair.IsTruePair(fPairDef->GetMotherPDG())) return 0.0;

    // computation variables
    Double_t mass1 = fPairDef->GetMass(0);
    Double_t mass2 = fPairDef->GetMass(1);
    Double_t invmass, invmassMC, weight;

    // assign nominal masses to daughters
    track1->SetM(mass1);
    track2->SetM(mass2);

    // if we are here, all cuts are passed - fill histograms
    invmass = fPair.GetInvMass(mass1, mass2);
    invmassMC = fPair.GetInvMassMC(mass1, mass2);
    if (fPIDMethod == AliRsnPID::kNone) {
        fHistogram->Fill(invmass);
        fHistogramMC->Fill(invmassMC);
        if (fPairDef->IsLikeSign() && !fPairDef->HasEqualTypes()) {
            track1->SetM(mass1);
            track2->SetM(mass2);
            if (fCuts->IsSelected(AliRsnCut::kPair, &fPair)) {
                invmass = fPair.GetInvMass(mass2, mass1);
                invmassMC = fPair.GetInvMassMC(mass2, mass1);
                fHistogram->Fill(invmass);
                fHistogramMC->Fill(invmassMC);
            }
        }
    }
    else {
        weight = fPairDef->ComputeWeight(track1, track2);
        if (weight > 0.0) {
            fHistogram->Fill(invmass, weight);
            fHistogramMC->Fill(invmassMC, weight);
        }
        // if we are treating a like-sign pair with different types for track #1 and #2,
        // we must also exchange the tracks and fill again the histogram
        if (fPairDef->IsLikeSign() && !fPairDef->HasEqualTypes()) {
            track1->SetM(mass1);
            track2->SetM(mass2);
            if (fCuts->IsSelected(AliRsnCut::kPair, &fPair)) {
                weight = fPairDef->ComputeWeight(track2, track1);
                invmass = fPair.GetInvMass(mass2, mass1);
                invmassMC = fPair.GetInvMassMC(mass2, mass1);
                if (weight > 0.0) {
                    fHistogram->Fill(invmass, weight);
                    fHistogramMC->Fill(invmassMC, weight);
                }
            }
        }
    }

    return 1.0;
}
