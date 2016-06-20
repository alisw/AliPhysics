//
// Created by mgrochow on 9/3/15.
//

#ifndef ALIROOT_CONVERSIONCONSTANTS_H
#define ALIROOT_CONVERSIONCONSTANTS_H


enum TrackType {
    kStandard,

    kKinkMother,
    kKinkDaughter,

    kV0NegativeDaughter,
    kV0PositiveDaughter,
    kV0Mother,

    kCascadePrimaryMother,
    kCascadePrimaryDaughter,
    kCascadeSecondaryMother,
    kCascadeNegativeDaughter,
    kCascadePositiveDaughter,

    kMuonMatched,
    kMuonNotMatched,
    kMuonGhost
};

const TString fgkDetector[23] = {
        "Invalid Layer",
        "First Layer",

        "SPD1",
        "SPD2",
        "SDD1",
        "SDD2",
        "SSD1",
        "SSD2",

        "TPC1",
        "TPC2",

        "TRD1",
        "TRD2",
        "TRD3",
        "TRD4",
        "TRD5",
        "TRD6",

        "TOF",

        "PHOS1",
        "PHOS2",

        "HMPID",
        "MUON",
        "EMCAL",
        "LastLayer"
};

const TString fgkTrackTypes[14] = {
        "standard",

        "kink_mother",
        "kink_daughter",

        "V0_negative_daughter",
        "V0_positive_daughter",
        "V0_mother",

        "cascade_primary_mother",
        "cascade_primary_daughter",
        "cascade_secondary_mother",
        "cascade_negative_daughter",
        "cascade_positive_daughter",

        "muon_matched",
        "muon_not_matched",
        "muon_ghost"
};

#endif //ALIROOT_CONVERSIONCONSTANTS_H
