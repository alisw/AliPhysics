#ifndef ALICONVEVENTCUTS_H
#define ALICONVEVENTCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Friederike Bock, Daniel Muehlheim
#include <TObjString.h>
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVTrack.h"
#include "AliAnalysisCuts.h"
#include "AliEMCALGeometry.h"
#include "AliDataFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TObjArray.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"
#include "TRandom3.h"
#include "AliVCaloTrigger.h"
#include "AliTimeRangeCut.h"

class AliESDEvent;
class AliAODEvent;
class TH1F;
class TH2F;
class TF1;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;
class AliAODMCParticle;
class AliEMCALTriggerPatchInfo;
class AliCaloTriggerMimicHelper;
class AliV0ReaderV1;
class AliAODConversionPhoton;

/**
 * @class AliConvEventCuts
 * @brief Class handling all kinds of selection cuts for Gamma Conversion analysis
 * @author Friederike Bock
 * @author Daniel Muehlheim
 * @ingroup GammaConv
 *
 * The cut configuration is set as a string with an 8 digit number.
 * Each digit in the string corresponds to a certain cut type, while
 * its values represent the cut values. The cut configuration is listed here:
 *
 * | Position in the cut string (from the end) | Cut type                     |
 * |-------------------------------------------|------------------------------|
 * |                  0                        | HeavyIon                     |
 * |                  1                        | CentralityMin                |
 * |                  2                        | CentralityMax                |
 * |                  3                        | SelectSpecialTrigger         |
 * |                  4                        | SelectSpecialSubTriggerClass |
 * |                  5                        | RemovePileUp                 |
 * |                  6                        | RejectExtraSignals           |
 * |                  7                        | VertexCut                    |
 */
class AliConvEventCuts : public AliAnalysisCuts {

    public:
      enum cutIds {
        kisHeavyIon,
        kCentralityMin,
        kCentralityMax,
        kSelectSpecialTriggerAlias,
        kSelectSubTriggerClass,
        kremovePileUp,
        kExtraSignals,
        kVertex,
        kNCuts
      };

      enum TriggerTypeEMCAL {
        kND       = -1,  //not defined
        kJ1       = 1,
        kJ2       = 2,
        kG1       = 3,
        kG2       = 4,
        kL0       = 5,
      };

      /**
       * @enum PeriodVar
       * @brief Collection of supported periods
       */
      enum PeriodVar {
        // data periods
        kNoPeriod=0,  //!< kNoPeriod
        // 2010
        kLHC10bg,         //!< pp 7 TeV (LHC10c incl 900 GeV)
        kLHC10h,          //!< PbPb 2.76TeV
        // MC's corresponding to 2010 data
        kLHC10d1,         //!< anchored LHC10b pass 2
        kLHC10d2,         //!< anchored LHC10b pass 2
        kLHC10d4a,        //! anchored LHC10c pass 2
        kLHC10d4,         //! anchored LHC10c pass 2
        kLHC10e12,        //! anchored LHC10c pass 2
        kLHC10e13,        //! anchored LHC10c pass 2
        kLHC10f6a,        //! anchored LHC10d pass 2
        kLHC10f6,         //!< anchored LHC10d pass 2
        kLHC10e20,        //!< anchored LHC10e pass 2
        kLHC10e21,        //!< anchored LHC10e pass 2
        kLHC14j4,         //!< anchored LHC10[b-g] pass 4
        kLHC14b7,         //!< anchored LHC11 pass 1
        kLHC13d2,         //!< anchored LHC10h pass 2
        kLHC13d2b,        //!< anchored LHC10h pass 2
        kLHC12a11a,       //!< anchored LHC10h pass 2
        kLHC12a11b,       //!< anchored LHC10h pass 2
        kLHC12a11c,       //!< anchored LHC10h pass 2
        kLHC12a11d,       //!< anchored LHC10h pass 2
        kLHC12a11e,       //!< anchored LHC10h pass 2
        kLHC12a11f,       //!< anchored LHC10h pass 2

        // 2011
        kLHC11a,          //!< pp 2.76TeV (part 7TeV)
        kLHC11b,          //!< pp 7TeV
        kLHC11cg,         //!< pp 7TeV
        kLHC11h,          //!< PbPb 2.76TeV
        // MC's corresponding to 2011 data
        kLHC12a15c,       //!< anchored LHC11a pass 2 - JJ
        kLHC12f1a,        //!< anchored LHC11a pass 4
        kLHC12f1b,        //!< anchored LHC11a pass 4
        kLHC12i3,         //!< anchored LHC11a pass 4
        kLHC15g1a,        //!< anchored LHC11a pass 4 - JJ
        kLHC15g1b,        //!< anchored LHC11a pass 4 - JJ
        kLHC13e4,         //!< anchored LHC11c pass 1 - GJ
        kLHC13e5,         //!< anchored LHC11c pass 1 - JJ
        kLHC14k1a,        //!< anchored LHC11[c-d] pass 1 - JJ
        kLHC14k1b,        //!< anchored LHC11[c-d] pass 1 - JJ
        kLHC12a15f,       //!< anchored LHC11d pass 1 - JJ
        kLHC12a15g,       //!< anchored LHC11d pass 1 - GJ
        kLHC12f2a,        //!< anchored LHC11d pass 1 - JJ
        kLHC14a1a,        //!< anchored LHC11h pass 2
        kLHC14a1b,        //!< anchored LHC11h pass 2
        kLHC14a1c,        //!< anchored LHC11h pass 2

        // 2012
        kLHC12,           //!< pp 8TeV
        // MC's corresponding to 2012 data
        kLHC14e2b,        //!< anchored LHC12[a-h] pass 1
        kLHC15h1,         //!< anchored LHC12[a-h] pass 2
        kLHC15h2,         //!< anchored LHC12[a-h] pass 2
        kLHC12P2JJ,       //!< anchored LHC12[a-h] pass 2 - JJ
        kLHC17g5b,        //!< anchored LHC12[a-h] pass 2 - dec gamma JJ
        kLHC17g5c,        //!< anchored LHC12[a-h] pass 2 - dec gamma JJ
        kLHC17g5a1,        //!< anchored LHC12[a-h] pass 2 - GJ Geant3
        kLHC17g5a2,        //!< anchored LHC12[a-h] pass 2 - GJ Geant3

        // 2013
        kLHC13bc,         //!< pPb 5.023TeV
        kLHC13de,         //!< pPb 5.023TeV
        kLHC13f,          //!< Pbp 5.023TeV
        kLHC13g,          //!< pp 2.76TeV
        // MC's corresponding to 2013 data
        kLHC13b2_efix,    //!< anchored LHC13[b-c] pass 2
        kLHC13e7,         //!< anchored LHC13[b-c] pass 2
        kLHC14b2,         //!< anchored LHC13[b-c] pass 2
        kLHC13b4_fix,     //!< anchored LHC13[b-c] pass 2 - JJ
        kLHC13b4_plus,    //!< anchored LHC13[b-c] pass 2 - JJ
        kLHC16c3a,        //!< anchored LHC13[d-e] pass 2 - JJ
        kLHC16c3b,        //!< anchored LHC13[d-e] pass 2 - JJ
        kLHC16c3c,        //!< anchored LHC13[d-e] pass 2 - GJ
        kLHC17g6a1,       //!< anchored LHC13[d-f] pass 2 - GJ
        kLHC17g6a2,       //!< anchored LHC13[d-f] pass 2 - JJ low
        kLHC17g6a3,       //!< anchored LHC13[d-f] pass 2 - JJ high
        kLHC18j5,         //!< anchored LHC13[b-c] pass 4 - General Purpose
        kLHC19a4,         //!< anchored LHC13[b-f] pass 4 - jj
        kLHC15g2,         //!< anchored LHC13g pass 1
        kLHC15a3a,        //!< anchored LHC13g pass 1 - JJ
        kLHC15a3a_plus,   //!< anchored LHC13g pass 1 - JJ
        kLHC15a3b,        //!< anchored LHC13g pass 1 - JJ
        kLHC15d3a,        //!< anchored LHC13g pass 1
        kLHC15d3b,        //!< anchored LHC13g pass 1
        // 2015
        kLHC15fm,         //!< pp 13 TeV
        kLHC15n,          //!< pp 5 TeV
        kLHC15o,          //!< PbPb 5 TeV
        kLHC18qr,         //!< PbPb 5 TeV
        // MC's corresponding to 2015 data
        kLHC15g3a3,       //!< anchored LHC15f pass 1
        kLHC15g3a,        //!< anchored LHC15f pass 1
        kLHC15g3c2,       //!< anchored LHC15f pass 1
        kLHC15g3c3,       //!< anchored LHC15f pass 1
        kLHC15g3,         //!< anchored LHC15f pass 1
        kLHC16a2a,        //!< anchored LHC15h pass 1
        kLHC16a2b,        //!< anchored LHC15h pass 1
        kLHC16a2c,        //!< anchored LHC15h pass 1
        kLHC15P2EPos,     //!< anchored LHC15f pass 2
        kLHC15P2Pyt8,     //!< anchored LHC15[h,i] pass 2
        kLHC15l1a2,       //!< anchored LHC15n pass 1
        kLHC15l1b2,       //!< anchored LHC15n pass 1
        kLHC15k1a1,       //!< LHC15o low IR firstPhysics
        kLHC15k1a2,       //!< LHC15o low IR firstPhysics
        kLHC15k1a3,       //!< LHC15o low IR firstPhysics
        kLHC16j7,         //!< LHC15o low IR pass4
        kLHC16g2,         //!< anchored LHC15o pass1 - general purpose EPOS-LHC
        kLHC16g3,         //!< anchored LHC15o pass1 - general purpose DPMJET
        kLHC16h4,         //!< anchored LHC15o pass1 - injected signals 0-100%
        kLHC16i1a,        //!< anchored LHC15o pass1 - LF added (multi-)strange 0-10%
        kLHC16i1b,        //!<                                                  10-50%
        kLHC16i1c,        //!<                                                  50-90%
        kLHC16i2a,        //!< anchored LHC15o pass1 - HF added hadronic decays 0-10%
        kLHC16i2b,        //!<                                                  10-50%
        kLHC16i2c,        //!<                                                  50-90%
        kLHC16i3a,        //!< anchored LHC15o pass1 - HF added electron decays 0-10%
        kLHC16i3b,        //!<                                                  10-50%
        kLHC16i3c,        //!<                                                  50-90%
        kLHC16h2a,        //!< anchored LHC15o pass1 - jet-jet 0-10%
        kLHC16h2b,        //!< anchored LHC15o pass1 - jet-jet 10-50%
        kLHC16h2c,        //!< anchored LHC15o pass1 - jet-jet 50-90%
        kLHC16h3,         //!< anchored LHC15n pass4 - jet-jet MC Pythia8 reproduction
        kLHC16h8a,        //!< anchored LHC15n pass2 - general purpose Pythia8
        kLHC16h8b,        //!< anchored LHC15n pass2 - general purpose Pythia6
        kLHC16k3a,        //!< anchored LHC15n pass2 - gen. purpose Pyt6wpileup
        kLHC16k3b,        //!< anchored LHC15o pass3 - gen. purpose Pyt6wpileup
        kLHC16k3a2,       //!< anchored LHC15n pass2 - gen. purpose Pyt6wopileup
        kLHC16k3b2,       //!< anchored LHC15o pass3 - gen. purpose Pyt6wopileup
        kLHC16k5a,        //!< anchored LHC15n pass3 - general purpose Pythia8
        kLHC16k5b,        //!< anchored LHC15n pass3 - general purpose Pythia6
        kLHC17e2,         //!< anchored LHC15n pass4 - general purpose Pythia8
        kLHC18j3,         //!< anchored LHC15n pass4 - general purpose Pythia8
        kLHC15k5a,        //!< anchored LHC15f pass2 - HF-forced MC for D2H analyses
        kLHC15k5b,        //!< anchored LHC15f pass2 - HF-forced MC for HFE analyses
        kLHC15k5c,        //!< anchored LHC15f pass2 - HF-forced MC for HFCJ analyses
        kLHC18b11a,       //!< anchored to LHC15o    - gamma-jets Pythia events embedded in HI MC events
        kLHC18b11b,       //!< anchored to LHC15o    - gamma-jets Pythia events embedded in HI MC events
        kLHC18b11c,       //!< anchored to LHC15o    - gamma-jets Pythia events embedded in HI MC events
        kLHC18e1,         //!< anchored to LHC15o    - general purpose - fixed MC
        kLHC18e1a,        //!< anchored LHC15o pass1 - general purpose - 0-10%
        kLHC18e1b,        //!< anchored LHC15o pass1 - general purpose - 10-50%
        kLHC18e1c,        //!< anchored LHC15o pass1 - general purpose - 50-90%
        kLHC18l8a,        //!< anchored to LHC18qr    - general purpose Pythia8
        kLHC18l8b,        //!< anchored to LHC18qr    - general purpose Pythia8
        kLHC18l8c,        //!< anchored to LHC18qr    - general purpose Pythia8
        kLHC19h2a,        //!< anchored to LHC18qr    - general purpose Pythia8
        kLHC19h2b,        //!< anchored to LHC18qr    - general purpose Pythia8
        kLHC19h2c,        //!< anchored to LHC18qr    - general purpose Pythia8
        kLHC19h3,         //!< anchored to LHC18qr    - general purpose Pythia8 with added GA signals
        kLHC20e3a,        //!< anchored to LHC18qr pass3 - general purpose Pythia8
        kLHC20e3b,        //!< anchored to LHC18qr pass3 - general purpose Pythia8
        kLHC20e3c,        //!< anchored to LHC18qr pass3 - general purpose Pythia8
        kLHC20g10,        //!< anchored to LHC18qr pass3 - general purpose Pythia8, with added GA signals

        // MC upgrade
        kLHC13d19,        //!< upgrade 5.5TeV PbPb

        // 2016
        kLHC16NomB,         //!< pp 13 TeV nominal B field
        kLHC16LowB,         //!< pp 13 TeV low B field
        kLHC16qt,           //!< pPb 5 TeV
        kLHC16r,            //!< pPb 8 TeV
        kLHC16s,            //!< pPb 8 TeV
        // MC's corresponding to 2016 data
        kLHC16P1Pyt8,       //!< anchored LHC16x pass 1 nom B-field - general purpose Pythia8
        kLHC16P1Pyt8LowB,   //!< anchored LHC16f pass 1 low B-field - general purpose Pythia8
        kLHC16P1EPOS,       //!< anchored LHC16x pass 1 nom B-field - general purpose EPOS
        kLHC16P1PHO,        //!< anchored LHC16d pass 1 nom B- field - for MBW Phojet
        kLHC16P1JJ,         //!< anchored LHC16x pass 1 nom B-field - Pythia8 JJ
        kLHC16P1JJLowB,     //!< anchored LHC16f pass 1 low B-field - Pythia8 JJ
        kLHC17h8a,          //!< anchored LHC16d,e,g,h,j,o,p pass 1 - heavy flavour MC Pythia6
        kLHC17h8b,          //!< anchored LHC16d,e,g,h,j,o,p pass 1 - heavy flavour MC Pythia6
        kLHC17h8c,          //!< anchored LHC16i,j,o,p pass 1 - heavy flavour MC Pythia6
        kLHC17c3b1,         //!< anchored LHC16k pass 1 - heavy flavour MC Pythia6
        kLHC17c3a1,         //!< anchored LHC16k pass 1 - heavy flavour MC Pythia6
        kLHC17c3b2,         //!< anchored LHC16l pass 1 - heavy flavour MC Pythia6
        kLHC17c3a2,         //!< anchored LHC16l pass 1 - heavy flavour MC Pythia6
        kLHC17i3a1,         //!< anchored LHC16i,j,k,l,o,p GammaJet - EMCal triggered
        kLHC17i3b1,         //!< anchored LHC16i,j,k,l,o,p JetJet - 3.5 GeV in EMCal acc.
        kLHC17i3b2,         //!< anchored LHC16i,j,k,l,o,p JetJet - 3.5 GeV in DCal/PHOS acc.
        kLHC17i3c1,         //!< anchored LHC16i,j,k,l,o,p JetJet - 7 GeV in EMCal acc.
        kLHC17i3c2,         //!< anchored LHC16i,j,k,l,o,p JetJet - 7 GeV in DCal/PHOS acc.
        kLHC20b1b1,         //!< anchored LHC16i,j,k,l,o,p JetJet - 3.5 GeV in EMCal acc. new prod.
        kLHC20b1b2,         //!< anchored LHC16i,j,k,l,o,p JetJet - 3.5 GeV in DCal/PHOS acc. new prod.
        kLHC20b1c1,         //!< anchored LHC16i,j,k,l,o,p JetJet - 7 GeV in EMCal acc. new prod.
        kLHC20b1c2,         //!< anchored LHC16i,j,k,l,o,p JetJet - 7 GeV in DCal/PHOS acc. new prod.

        //General purpose- pPb
        kLHC17a3a,            //!< anchored LHC16r pass 1 - general purpose EPOSLHC
        kLHC17a3b,            //!< anchored LHC16r pass 1 - general purpose DPMJET
        kLHC17a4a,            //!< anchored LHC16s pass 1 - general purpose EPOSLHC
        kLHC17a4b,            //!< anchored LHC16s pass 1 - general purpose DPMJET
        kLHC17g6b2a,          //!< anchored LHC16rs pass 1 - decay gamma 3.5 GeV EMCal GeV JJ
        kLHC17g6b2b,          //!< anchored LHC16rs pass 1 - decay gamma 3.5 GeV DCal/PHOS GeV JJ
        kLHC17g6b3a,          //!< anchored LHC16rs pass 1 - decay gamma 7 GeV EMCal GeV JJ
        kLHC17g6b3b,          //!< anchored LHC16rs pass 1 - decay gamma 7 GeV DCal/PHOS GeV JJ
        kLHC18f3bc,           //!< anchored LHC16rs pass 1 - general purpose DPMJET
        kLHC17f2a,            //!< anchored LHC16qt pass 1 - general purpose EPOSLHC
        kLHC17f2b,            //!< anchored LHC16qt pass 1 - general purpose DPMJET
        kLHC18f3,             //!< anchored LHC16qt pass 1 - general purpose DPMJET
        kLHC17g8a,            //!< anchored LHC16qt pass 1 - jet-jet MC in EPOSLHC
        kLHC17f3,             //!< anchored LHC16r pass 1 - general purpose
        kLHC17f3a,            //!< anchored LHC16r pass 1 - general purpose EPOSLHC
        kLHC17f3b,            //!< anchored LHC16r pass 1 - general purpose DPMJET
        kLHC17f4,             //!< anchored LHC16s pass 1 - general purpose
        kLHC17f4a,            //!< anchored LHC16s pass 1 - general purpose EPOSLHC
        kLHC17f4b,            //!< anchored LHC16s pass 1 - general purpose DPMJET
        kLHC16rP1JJ,          //!< anchored LHC16r pass 1 - jet-jet MC in EPOSLHC
        kLHC16sP1JJ,          //!< anchored LHC16s pass 1 - jet-jet MC in EPOSLHC
        kLHC16rsGJ,           //!< anchored LHC16rs pass 1 - Gamma-jet MC in EMCal acc

        //heavy flavour MC pPb k17d2a_fast,
        kLHC17d2a,          //!< anchored LHC16q,t pass 1 - heavy flavour MC Hijing, fast only
        kLHC17d2b,          //!< anchored LHC16q,t pass 1 - heavy flavour MC Hijing, fast only

        // 2017
        kLHC17NomB,           //!< pp 13 TeV nominal B field
        kLHC17LowB,           //!< pp 13 TeV low B field
        kLHC17n,              //!< Xe-Xe 5.44 TeV
        kLHC17pq,             //!< pp 5 TeV
        // MC Xe-Xe
        kLHC17XeXeHi,         //!< MC for Xe-Xe 5.44 TeV HIJING
        // 5 TeV MC 2017
        kLHC17l3b,            //!< anchored LHC17p/q pass 1 - general purpose w/GEANT3,
        kLHC18j2,             //!< anchored LHC17p/q pass 1 - general purpose w/GEANT3,
        kLHC17l4b,            //!< anchored LHC17p/q pass 1 - general purpose w/GEANT4,
        kLHC18b8,             //!< anchored LHC17p/q pass 1 - jet-jet MC w/GEANT3,
        kLHC18b10,            //!< anchored LHC17p/q pass 1 - gamma-jet MC w/GEANT3,
        kLHC18l2,             //!< anchored LHC17p/q pass 1 - gamma-jet MC w/GEANT3,
        kLHC17P1PHO,          //!< anchored LHC17p only low Intensity Phojet 5 TeV
        //13 TeV MC 2017
        kLHC17P1Pyt8NomB,     //!LHC17x Pythia8 MB productions nom B anchored to LHC17x
        kLHC17P1Pyt6NomB,     //!LHC17x Pythia8 MB productions nom B anchored to LHC17x
        kLHC17P1PHONomB13TeV, //!LHC17x Phojet MB productions nom B anchored to LHC17x
        kLHC17P1Pyt8LowB,     //!LHC17x Pythia8 MB productions low B anchored to LHC17g
        kLHC17j5a,            //!LHC17k Strangeness enhanced
        kLHC17j5b,            //!LHC17l Strangeness enhanced
        kLHC17j5c,            //!LHC17o Strangeness enhanced
        //13 TeV LHC2017 JJ
        kLHC17P1JJ,           //!LHC17k JJ
        kLHC17P1JJLowB,       //!LHC17k JJ
        kLHC18l6b1,           //!JJ MC anchored to LHC17 with decay photon > 3.5 GeV in EMCal acc.
        kLHC18l6b2,           //!JJ MC anchored to LHC17 with decay photon > 3.5 GeV in DCal/PHOS acc.
        kLHC18l6c1,           //!JJ MC anchored to LHC17 with decay photon > 7 GeV in EMCal acc.
        kLHC18l6c2,           //!JJ MC anchored to LHC17 with decay photon > 7 GeV in DCal/PHOS acc.
        // 2018
        kLHC18NomB,           //!< pp 13 TeV nominal B field
        kLHC18LowB,           //!< pp 13 TeV low B field
        kLHC18P1JJ,           //!< pp 13 TeV JJ MCs
        kLHC19i3b1,           //!JJ MC anchored to LHC17 with decay photon > 3.5 GeV in EMCal acc.
        kLHC19i3b2,           //!JJ MC anchored to LHC17 with decay photon > 3.5 GeV in DCal/PHOS acc.
        kLHC19i3c1,           //!JJ MC anchored to LHC17 with decay photon > 7 GeV in EMCal acc.
        kLHC19i3c2,           //!JJ MC anchored to LHC17 with decay photon > 7 GeV in DCal/PHOS acc.

        //13 TeV LHC2018
        kLHC18P1Pyt8NomB,     //!LHC18x Pythia8 MB productions nom B anchored to LHC18x
        kLHC18P1Pyt8LowB,     //!LHC18x Pythia8 MB productions low B anchored to LHC18c

        kUnknownPeriod//!< kUnknownPeriod
      };

      /**
       * @enum EnergyVar
       * @brief Supported collision systems
       */
      enum EnergyVar {
        kUnset        = 0,   //!< not defined
        k900GeV       = 1,   //!< pp 900 GeV
        k2760GeV      = 2,   //!< pp 2.76TeV
        k5TeV         = 3,   //!< pp 5 TeV
        k7TeV         = 4,   //!< pp 7 TeV
        k8TeV         = 5,   //!< pp 8 TeV
        k13TeV        = 6,   //!< pp 13 TeV
        k13TeVLowB    = 7,   //!< pp 13 TeV low B
        kpPb5TeV      = 8,   //!< pPb 5 TeV
        kpPb8TeV      = 9,   //!< pPb 8 TeV
        kPbPb2760GeV  = 10,  //!< PbPb 2.76TeV
        kPbPb5TeV     = 11,  //!< PbPb 5 TeV
        kXeXe5440GeV  = 12,  //!< XeXe 5.44 TeV
        kpPb5TeVR2    = 13   //!< pPb 5 TeV run 2

      };

      enum phosTriggerType{kPHOSAny,kPHOSL0,kPHOSL1low,kPHOSL1med,kPHOSL1high} ;


      AliConvEventCuts(const char *name="EventCuts", const char * title="Event Cuts");
      AliConvEventCuts(const AliConvEventCuts&);
      AliConvEventCuts& operator=(const AliConvEventCuts&);

      virtual ~AliConvEventCuts();                            //virtual destructor

  //    static AliConvEventCuts * GetStandardCuts2010PbPb();
  //    static AliConvEventCuts * GetStandardCuts2010pp();

      Int_t     fCuts[kNCuts];
      Bool_t    UpdateCutString();
      static const char * fgkCutNames[kNCuts];

      // Setters
      Bool_t    SetCutIds (TString cutString);
      Bool_t    SetCut (cutIds cutID, Int_t cut);
      Bool_t    SetIsHeavyIon (Int_t isHeavyIon);
      Bool_t    SetCentralityMax (Int_t centralityBin);
      Bool_t    SetCentralityMin (Int_t centralityBin);
      Bool_t    SetRemovePileUp (Int_t removePileUp);
      Bool_t    SetMultiplicityMethod (Int_t multiplicityMethod);
      Bool_t    SetSelectSpecialTrigger (Int_t selectSpecialTrigger);
      Bool_t    SetSelectSubTriggerClass (Int_t selectSpecialSubTriggerClass);
      Bool_t    SetRejectExtraSignalsCut (Int_t extraSignal);
      Bool_t    SetVertexCut(Int_t vertexCut);

      void    SetPeriodEnum (TString periodName);
      void    SetPeriodEnumExplicit ( PeriodVar periodEnum )                        { fPeriodEnum = periodEnum                                  ; }
      void    SetCorrectionTaskSetting(TString setting)                             { fCorrTaskSetting = setting                                ; }
      void    SetTriggerMimicking(Int_t value)                                      { fMimicTrigger = value                                     ;
                                                                                      if(value)AliInfo("enabled trigger mimicking")             ; }
      void    SetTriggerOverlapRejecion (Bool_t value)                              { fRejectTriggerOverlap = value                             ;
                                                                                      if(value)AliInfo("enabled trigger overlap rejection")     ; }
      void    SetPHOSTrigger(phosTriggerType t=kPHOSL0)                             { fPHOSTrigger=t                                            ; }

      void    SetV0ReaderName (TString name)                                        { fV0ReaderName = name                                      ; }
      void    SetCaloTriggerHelperName (TString name)                               { CaloTriggerHelperName = name                                      ; }

      void    SetAddedSignalPDGCode (Int_t addedSignalPDGcode)                      { fAddedSignalPDGCode = addedSignalPDGcode                  ; }
      void    SetPreSelectionCutFlag (Bool_t preSelFlag)                            { fPreSelCut = preSelFlag                                   ; }
      void    SetCaloTriggerPatchInfoName(const char *n)                            { fCaloTriggerPatchInfoName = n                             ; }
      void    SetCaloTriggersName(const char *n)                                    { fCaloTriggersName  = n                                    ; }
      void    SetAcceptedHeader(TList *HeaderList)                                  { fHeaderList = HeaderList                                  ; }
      void    SetFillCutHistograms( TString name="",
                                    Bool_t preCut = kTRUE)                          { if(!fHistograms){ InitCutHistograms(name,preCut);}        ; }
      void    SetEtaShift(Double_t etaShift)                                        { fEtaShift = etaShift                                      ; } // Eta shift Setting
      void    SetUseJetFinderForOutliers(Bool_t useJetFinder)                       { fUseJetFinderForOutlier = useJetFinder                    ; } // Eta shift Setting
      void    SetUsePtHardBinFromFile(Bool_t useFilePathPth)                        { fUseFilePathForPthard = useFilePathPth                    ; } // Eta shift Setting
      void    SetUseAdditionalOutlierRejection(Bool_t useOutlierRej)                { fUseAdditionalOutlierRejection = useOutlierRej            ; } // Eta shift Setting
      void    SetEtaShift(TString pPbOrPbp)                                         { Double_t etaShift = 0.0                                   ;
                                                                                      if(!pPbOrPbp.CompareTo("pPb"))      etaShift = -0.465     ;
                                                                                      else if(!pPbOrPbp.CompareTo("Pbp")) etaShift =  0.465     ;
                                                                                      fEtaShift = etaShift                                      ; }
  //    void GetHistoCentralityFlattening
      void    SetUseWeightFlatCentralityFromFile( Int_t doFlattening = 1,
                              TString pathC="$ALICE_PHYSICS/PWGGA/GammaConv/InterpValuesAndFlattening.root",
                              TString histoCentNotFlat="")
                                                                                    {
                                                                                      AliInfo(Form("enabled centrality flattening with weights from file: %s",pathC.Data()));
                                                                                      fDoCentralityFlat = doFlattening                          ;
                                                                                      fPathWeightsFlatCent=pathC                                ;
                                                                                      fNameHistoNotFlatCentrality = histoCentNotFlat            ;
                                                                                    }
      void    SetCustomTriggerMimicOADBFile(TString pathOADB="")
                                                                                    {
                                                                                      AliInfo(Form("setting custom trigger mimic OADB from file: %s",pathOADB.Data()));
                                                                                      fPathTriggerMimicSpecialInput=pathOADB                                ;
                                                                                    }
      void    SetUseReweightingWithHistogramFromFile( Bool_t pi0reweight=kTRUE,
                                Bool_t etareweight=kFALSE,
                                Bool_t k0sreweight=kFALSE,
                                                              TString path="$ALICE_PHYSICS/PWGGA/GammaConv/MCSpectraInput.root",
                                TString histoNamePi0 = "",
                                TString histoNameEta = "",
                                TString histoNameK0s = "",
                                TString fitNamePi0 = "",
                                TString fitNameEta = "",
                                TString fitNameK0s ="" )
                                                                                    {
                                                                                      AliInfo(Form("enabled reweighting for: pi0 : %i, eta: %i, K0s: %i",pi0reweight, etareweight, k0sreweight));
                                                                                      fDoReweightHistoMCPi0 = pi0reweight                       ;
                                                                                      fDoReweightHistoMCEta = etareweight                       ;
                                                                                      fDoReweightHistoMCK0s = k0sreweight                       ;
                                                                                      fPathTrFReweighting=path                                  ;
                                                                                      fNameHistoReweightingPi0 =histoNamePi0                    ;
                                                                                      fNameHistoReweightingEta =histoNameEta                    ;
                                                                                      fNameHistoReweightingK0s =histoNameK0s                    ;
                                                                                      fNameFitDataPi0 =fitNamePi0                               ;
                                                                                      fNameFitDataEta =fitNameEta                               ;
                                                                                      fNameFitDataK0s =fitNameK0s                               ;
                                                                                    }
      void    SetUseWeightMultiplicityFromFile( Int_t doWeighting = 0,
                                                TString pathC="$ALICE_PHYSICS/PWGGA/GammaConv/MultiplicityInput.root",
                                                TString nameHistoMultData="",
                                                TString nameHistoMultMC=""
                                              )
                                                                                    {
                                                                                      AliInfo(Form("enabled multiplicity weights from file: %s",pathC.Data()));
                                                                                      fDoMultiplicityWeighting = doWeighting                    ;
                                                                                      fPathReweightingMult=pathC                                ;
                                                                                      fNameHistoReweightingMultData = nameHistoMultData         ;
                                                                                      fNameHistoReweightingMultMC = nameHistoMultMC             ;
                                                                                    }
      void    SetUseGammaPtReweightingWithHistogramFromFile(  Bool_t gammareweight=kTRUE,
                                                              TString path="$ALICE_PHYSICS/PWGGA/GammaConv/MCGammaSpectraInput.root",
                                                              TString histoNameGamma = "",
                                                              TString histoDataNameGamma = "")
                                                                                    {
                                                                                      AliInfo(Form("enabled pT reweighting for: Gamma : %i ", gammareweight));
                                                                                      fDoReweightHistoMCGamma = gammareweight                   ;
                                                                                      fPathTrFGammaReweighting = path                           ;
                                                                                      fNameHistoReweightingGamma = histoNameGamma               ;
                                                                                      fNameDataHistoReweightingGamma = histoDataNameGamma       ;
                                                                                    }

      void    SetMinFacPtHard(Float_t value)                                        { fMinFacPtHard = value                                     ;
                                                                                      AliInfo(Form("minimum factor between pt hard and jet put to: %2.2f",fMinFacPtHard));
                                                                                    }
      void    SetMaxFacPtHard(Float_t value)                                        { fMaxFacPtHard = value                                     ;
                                                                                      AliInfo(Form("maximum factor between pt hard and jet put to: %2.2f",fMaxFacPtHard));
                                                                                    }
      void    SetMaxFacPtHardSingleParticle(Float_t value)                          { fMaxFacPtHardSingleParticle = value                       ;
                                                                                      AliInfo(Form("maximum factor between pt hard and pt of pi0 or eta put to: %2.2f",fMaxFacPtHardSingleParticle));
                                                                                    }
      void    SetDebugLevel( Int_t value)                                           { fDebugLevel = value                                       ; }

      // Geters
      AliV0ReaderV1* GetV0Reader();
      TString   GetCutNumber();
      TString*  GetFoundHeader()                                                    { return fGeneratorNames                                    ; }
      Int_t     GetEventQuality()                                                   { return fEventQuality                                      ; }
      Bool_t    GetIsFromPileupSPD()                                                { return fRemovePileUpSPD                                   ; }
      Int_t     GetUseSphericity()                                                  { return fUseSphericity                                     ; }
      Bool_t    GetUseSphericityTrue()                                              { return fUseSphericityTrue                                 ; }
      Int_t     GetPastFutureLowBC()                                                { return fPastFutureRejectionLow                            ; }
      Int_t     GetPastFutureHighBC()                                               { return fPastFutureRejectionHigh                           ; }
      Bool_t    GetDoPileUpRejectV0MTPCout()                                        { return fDoPileUpRejectV0MTPCout                           ; }
      void      GetCentralityRange(Double_t range[2])                               { range[0]=10*fCentralityMin                                ;
                                                                                      range[1]=10*fCentralityMax                                ; }
      TList*    GetCutHistograms()                                                  { return fHistograms                                        ; }
      Int_t     GetMultiplicityMethod()                                             { return fMultiplicityMethod                                ; }
      Int_t     GetSignalRejection()                                                { return fRejectExtraSignals                                ; }
      Int_t     GetNAcceptedHeaders()                                               { return fnHeaders                                          ; }
      TString * GetAcceptedHeaderNames()                                            { return fGeneratorNames                                    ; }
      Int_t *   GetAcceptedHeaderStart()                                            { return fNotRejectedStart                                  ; }
      Int_t *   GetAcceptedHeaderEnd()                                              { return fNotRejectedEnd                                    ; }
      Int_t     GetAcceptedHeaderStart(Int_t headernumber)                          { if (headernumber < fnHeaders)
                                                                                        return fNotRejectedStart[headernumber]                  ;
                                                                                      else
                                                                                        return -1                                               ;
                                                                                    }
      Int_t     GetAcceptedHeaderEnd(Int_t headernumber)                            { if (headernumber < fnHeaders)
                                                                                        return fNotRejectedEnd[headernumber]                    ;
                                                                                      else
                                                                                        return -1                                               ;
                                                                                    }
      TList*    GetAcceptedHeader()                                                 { return fHeaderList                                        ; }
      Int_t     GetNumberOfContributorsVtx(AliVEvent *event);
      Double_t  GetEtaShift()                                                       { return fEtaShift                                          ; }
      Bool_t    GetDoEtaShift()                                                     { return fDoEtaShift                                        ; }
      Bool_t    GetUseJetFinderForOutliers()                                        { return fUseJetFinderForOutlier                            ; }
      Bool_t    GetUsePtHardBinFromFile()                                           { return fUseFilePathForPthard                              ; }
      TString   GetSpecialTriggerName()                                             { return fSpecialTriggerName                                ; }
      const TString& GetLabelNamePileupCutTPC() const                               { return fLabelNamePileupCutTPC                             ; }
      AliEMCALTriggerPatchInfo   *GetMainTriggerPatch();
      ULong_t   GetTriggerList();
      phosTriggerType GetPHOSTrigger()                                              { return fPHOSTrigger                                       ; }
      Int_t    GetTriggerMimicking()                                                { return fMimicTrigger                                      ; }
      Float_t   GetWeightForCentralityFlattening(AliVEvent *event = 0x0);
      Float_t   GetWeightForMultiplicity(Int_t mult);
      Float_t   GetWeightForMeson( Int_t index, AliMCEvent *mcEvent, AliVEvent *event = 0x0);
      Float_t   GetWeightForGamma( Int_t index, Double_t gammaPTrec, AliMCEvent *mcEvent, AliVEvent *event = 0x0);
      Float_t   GetCentrality(AliVEvent *event);
      Bool_t    GetUseNewMultiplicityFramework();
      void      GetCorrectEtaShiftFromPeriod();
      void      GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliVEvent *event);
      Double_t  GetV0Multiplicity(AliVEvent *event) const;
      Int_t     GetNumberOfTPCClusters(AliVEvent *event) const;
      TClonesArray*     GetArrayFromEvent(AliVEvent* event, const char *name, const char *clname=0);
      AliEMCALGeometry* GetGeomEMCAL()                                              { return fGeomEMCAL;}

      Bool_t    InitializeCutsFromCutString(const TString analysisCutSelection);
      void      SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kAny) {
                                                                                      fOfflineTriggerMask = offlineTriggerMask                  ;
                                                                                      fTriggerSelectedManually = kTRUE                          ;
                                                                                    }
      void    SelectSpecialTrigger( UInt_t offlineTriggerMask = AliVEvent::kAny,
                                    TString TriggerClassName = "AliVEvent::kAny" ) {
                                                                                      fOfflineTriggerMask = offlineTriggerMask                  ;
                                                                                      fSpecialTriggerName = TriggerClassName                    ;
                                                                                      AliInfo(fSpecialTriggerName)                              ;

                                                                                    }

      PeriodVar GetPeriodEnum ()                                                    { return fPeriodEnum                                        ; }
      EnergyVar GetEnergyEnum ()                                                    { return fEnergyEnum                                        ; }
      virtual   Bool_t IsSelected(TObject* /*obj*/)                                 { return kTRUE                                              ; }
      virtual   Bool_t IsSelected(TList* /*list*/)                                  { return kTRUE                                              ; }


      // Cut Selection
      Bool_t    EventIsSelected(AliVEvent *fInputEvent,
                                AliMCEvent *fMCEvent);
      Int_t     IsEventAcceptedByCut( AliConvEventCuts *ReaderCuts,
                                      AliVEvent *event,
                                      AliMCEvent *mcEvent,
                                      Int_t isHeavyIon,
                                      Bool_t isEMCALAnalysis);

      void    PrintCuts();
      void    PrintCutsWithValues();
      void    InitCutHistograms(  TString name="",
                                  Bool_t preCut = kTRUE);
      void    SetLightOutput( Int_t flag ){fDoLightOutput = flag; return;}
      void    SetUseSphericityTrue( Bool_t flag ){fUseSphericityTrue = flag;}
      void    FillTPCPileUpHistograms(AliVEvent *event);

      ///Cut functions
      Int_t   IsParticleFromBGEvent(  Int_t index,
                                      AliMCEvent *mcEvent,
                                      AliVEvent *event = 0x0,
                                      Int_t debug = 0
                                   );
      TString   GetParticleHeaderName(  Int_t index,
                                      AliMCEvent *mcEvent,
                                      AliVEvent *event = 0x0,
                                      Int_t debug = 0
                                   );

      Bool_t PhotonPassesAddedParticlesCriterion(AliMCEvent             *theMCEvent,
                                                 AliVEvent              *theInputEvent,
                                                 AliAODConversionPhoton &thePhoton,
                                                 Bool_t                 &theIsFromSelectedHeader); // future todo: make this const

      void    LoadWeightingFlatCentralityFromFile ();
      void    LoadWeightingMultiplicityFromFile ();
      void    LoadReweightingHistosMCFromFile ();
      void    LoadGammaPtReweightingHistosMCFromFile ();

      // Event Cuts
      Bool_t    IsCentralitySelected(AliVEvent *event, AliMCEvent *mcEvent);
      Bool_t    IsOutOfBunchPileupPastFuture(AliVEvent *event);
      Bool_t    IsPileUpV0MTPCout(AliVEvent *event);
      Bool_t    IsPileUpSDDSSDTPC(AliVEvent *event);
      Bool_t    VertexZCut(AliVEvent *event);
      Bool_t    IsJetJetMCEventAccepted(AliMCEvent *mcEvent, Double_t& weight, Float_t& pthard, AliVEvent* event = 0x0, Double_t maxJetPt = -1);
      Float_t   GetPtHard(AliMCEvent *mcEvent, AliVEvent* event = 0x0);
      Int_t     GetPtHardBinFromPath(const char* currFile, AliVEvent *event);
      void      GetXSectionAndNTrials(AliMCEvent *mcEvent, Float_t &XSection, Float_t &NTrials, AliVEvent* event = 0x0 );
      Float_t   GetMaxPtJet()                                                       { return fMaxPtJetMC                                        ; }
      Bool_t    MimicTrigger( AliVEvent *event,
                              Bool_t isMC );
      Bool_t    IsTriggerSelected(  AliVEvent *event,
                                    Bool_t isMC);
      Bool_t    HasV0AND()                                                          { return fHasV0AND                                          ; }
      Bool_t    IsSDDFired()                                                        { return fIsSDDFired                                        ; }
      Int_t     IsSpecialTrigger()                                                  { return fSpecialTrigger                                    ; }
      Int_t     IsSpecialSubTrigger()                                               { return fSpecialSubTrigger                                 ; }
      void      InitializeEMCALTrigger( AliVEvent *event);
      Bool_t    HasTriggerType(TriggerTypeEMCAL t);

      // Request Flags
      Int_t     IsHeavyIon()                                                        { return fIsHeavyIon                                        ; }
      void      DoEtaShift(Bool_t doEtaShift)                                       { fDoEtaShift = doEtaShift                                  ; }

      //MC particle flags - determine whether particle is primary or secondary
      Bool_t    IsConversionPrimaryESD( AliMCEvent *mcEvent,
                                        Long_t eventpos,
                                        Double_t prodVtxX,
                                        Double_t prodVtxY,
                                        Double_t prodVtxZ);
      Bool_t    IsConversionPrimaryAOD( AliVEvent *event,
                                        AliAODMCParticle* AODMCParticle,
                                        Double_t prodVtxX,
                                        Double_t prodVtxY,
                                        Double_t prodVtxZ);

      Int_t     SecondaryClassificationPhoton(  TParticle *particle,
                                                AliMCEvent *mcEvent,
                                                Bool_t isConversion );
      Int_t     SecondaryClassificationPhotonAOD( AliAODMCParticle *particle,
                                                  TClonesArray *aodmcArray,
                                                  Bool_t isConversion );

    protected:
      TList*                      fHistograms;                            ///<
      TList*                      fHeaderList;                            ///<

      Int_t                       fDoLightOutput;                         ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
      Int_t                       fEventQuality;                          ///< EventQuality
      AliEMCALGeometry*           fGeomEMCAL;                             ///< pointer to EMCal geometry
      TClonesArray*               fAODMCTrackArray;                       ///< pointer to track array
      AliV0ReaderV1*              fV0Reader;                              //!
      //cuts
      Int_t                       fIsHeavyIon;                            ///< flag for heavy ion
      Int_t                       fDetectorCentrality;                    ///< centrality detecotor V0M or CL1
      Int_t                       fModCentralityClass;                    ///< allows to select smaller centrality classes
      Bool_t                      fEnableVertexCut;                       ///< enable vertex cut
      Double_t                    fMaxVertexZ;                            ///< max z offset of vertex
      Int_t                       fCentralityMin;                         ///< centrality selection lower bin value
      Int_t                       fCentralityMax;                         ///< centrality selection upper bin value
      Int_t                       fMultiplicityMethod;                    ///< selected multiplicity method
      Int_t                       fSpecialTrigger;                        ///< flag
      Int_t                       fSpecialSubTrigger;                     ///< flag
      Bool_t                      fRemovePileUp;                          ///< flag specifies if any pileup cut is applied
      Bool_t                      fRemovePileUpSPD;                       ///< flag specifies if SPD pileup cuts are applied
      Int_t                       fUseSphericity;                         ///< flag that specifies the sphericityCut
      Bool_t                      fUseSphericityTrue;                     ///< switch for true sphericity cuts
      Int_t                       fPastFutureRejectionLow;                ///< sets bunch crossing event rejection in past
      Int_t                       fPastFutureRejectionHigh;               ///< sets bunch crossing event rejection in future. If both are 0, the cut is not applied
      Int_t                       fDoPileUpRejectV0MTPCout;               ///< reject event if # TPCout tracks does not follow expected V0M mult
      TF1 *                       fFPileUpRejectV0MTPCout;                ///< Pol1 function to compute the cut
      Bool_t                      fRemovePileUpSDDSSDTPC;                 //<  reject event if too many TPC clusters with respect to SDD+SSD clusters
      TF1 *                       fFPileUpRejectSDDSSDTPC;                //<  Pol2 function to compute cut
      Int_t                       fRejectExtraSignals;                    ///<
      UInt_t                      fOfflineTriggerMask;                    ///< Task processes collision candidates only
      Bool_t                      fHasV0AND;                              ///< V0AND Offline Trigger
      Bool_t                      fIsSDDFired;                            ///< SDD FIRED to select with SDD events
      TRandom3                    fRandom;                                ///<
      Int_t                       fnHeaders;                              ///< Number of Headers
      Int_t*                      fNotRejectedStart;                      //[fnHeaders]
      Int_t*                      fNotRejectedEnd;                        //[fnHeaders]
      TString*                    fGeneratorNames;                        //[fnHeaders]
      PeriodVar                   fPeriodEnum;                            ///< period selector
      EnergyVar                   fEnergyEnum;                            ///< energy selector
      AliTimeRangeCut             fTimeRangeCut;                          //!

      TObjString*                 fCutString;                             ///< cut number used for analysis
      TString                     fCutStringRead;                         ///<
      AliAnalysisUtils*           fUtils;                                 ///<
      Double_t                    fEtaShift;                              ///<
      Bool_t                      fDoEtaShift;                            ///< Flag for Etashift
      Bool_t                      fUseJetFinderForOutlier;                ///< Flag for Etashift
      Bool_t                      fUseFilePathForPthard;                  ///< Flag for Etashift
      Bool_t                      fUseAdditionalOutlierRejection;         ///< Flag for Etashift
      Int_t                       fDoCentralityFlat;                      ///<
      TString                     fPathWeightsFlatCent;                   ///<
      TString                     fNameHistoNotFlatCentrality;            ///<
      Bool_t                      fDoReweightHistoMCPi0;                  ///< Flag for reweighting Pi0 input with histogram
      Bool_t                      fDoReweightHistoMCEta;                  ///< Flag for reweighting Eta input with histogram
      Bool_t                      fDoReweightHistoMCK0s;                  ///< Flag for reweighting K0s input with histogram
      TString                     fPathTrFReweighting;                    ///< Path for file used in reweighting
      TString                     fNameHistoReweightingPi0;               ///< Histogram name for reweighting Pi0
      TString                     fNameHistoReweightingEta;               ///< Histogram name for reweighting Eta
      TString                     fNameHistoReweightingK0s;               ///< Histogram name for reweighting K0s
      TString                     fNameFitDataPi0;                        ///< Fit name for fit to spectrum of pi0s in Data
      TString                     fNameFitDataEta;                        ///< Fit name for fit to spectrum of etas in Data
      TString                     fNameFitDataK0s;                        ///< Fit name for fit to spectrum of k0s in Data
      Bool_t                      fDoReweightHistoMCGamma;                ///< Flag for reweighting Gamma input with histogram
      TString                     fPathTrFGammaReweighting;               ///< Path for file used in gamma reweighting
      TString                     fNameHistoReweightingGamma;             ///< Histogram name for reweighting Gamma
      TString                     fNameDataHistoReweightingGamma;         ///< Histogram Data name for reweighting Gamma
      TString                     fLabelNamePileupCutTPC;                 //<  Label for NEvents histograms depending on pileup cut used
      // Histograms
      TH1F*                       fHistoEventCuts;                        ///< bookkeeping for event selection cuts
      TH1F*                       fHistoPastFutureBits;                   ///< bookkeeping for event selection cuts
      TH1F*                       hCentrality;                            ///< centrality distribution for selected events
      TH1D*                       hCentralityNotFlat;                     ///< centrality distribution loaded for cent. flattening
      //TH2F*                      hCentralityVsNumberOfPrimaryTracks;    ///< centrality distribution for selected events
      TH1F*                       hVertexZ;                               ///< vertex z distribution for selected events
      TH1F*                       hNPileupVertices;                       ///< number of SPD pileup vertices
      TH1F*                       hPileupVertexToPrimZ;                   ///< distance of SPD pileup vertex to prim vertex in z
      TH1F*                       hPileupVertexToPrimZSPDPileup;          ///< distance of SPD pileup vertex to prim vertex in z for SPD pileup flagged events
      TH1F*                       hPileupVertexToPrimZTrackletvsHits;     ///< distance of SPD pileup vertex to prim vertex in z for Tracklet vs Hits flagged events
      TH1F*                       hEventPlaneAngle;                       ///<
      Double_t                    fEventPlaneAngle;                       ///< EventPlaneAngle
      TH1F*                       hTriggerClass;                          ///< fired offline trigger class
      TH1F*                       hTriggerClassSelected;                  ///< selected fired offline trigger class
      TH1F*                       hTriggerClassesCorrelated;              ///< selected trigger class correlation with others
      TH1D*                       hReweightMCHistPi0;                     ///< histogram input for reweighting Pi0
      TH1D*                       hReweightMCHistEta;                     ///< histogram input for reweighting Eta
      TH1D*                       hReweightMCHistK0s;                     ///< histogram input for reweighting K0s
      TF1*                        fFitDataPi0;                            ///< fit to pi0 spectrum in Data
      TF1*                        fFitDataEta;                            ///< fit to eta spectrum in Data
      TF1*                        fFitDataK0s;                            ///< fit to K0s spectrum in Data
      TH1D*                       hReweightMCHistGamma;                   ///< histogram MC   input for reweighting Gamma
      TH1D*                       hReweightDataHistGamma;                 ///< histogram data input for reweighting Gamma
      Int_t                       fAddedSignalPDGCode;
      Bool_t                      fPreSelCut;                             // Flag for preselection cut used in V0Reader
      Bool_t                      fTriggerSelectedManually;               // Flag for manual trigger selection
      TString                     fSpecialTriggerName;                    // Name of the Special Triggers
      TString                     fSpecialSubTriggerName;                 // Name of the Special Triggers
      TString                     fSpecialSubTriggerNameAdditional;       // Name of an additional Special Trigger
      Int_t                       fNSpecialSubTriggerOptions;
      TH2F*                       hSPDClusterTrackletBackgroundBefore;    ///< SPD tracklets vs SPD clusters for background-correction before cut
      TH2F*                       hSPDClusterTrackletBackground;          ///< SPD tracklets vs SPD clusters for background-correction
      TH2F*                       hV0MultVsNumberTPCoutTracks;            ///< correlation V=Mult vs number TPC out Tracks
      TH2F*                       hTPCSDDSSDClusters;                     ///< x: TPC clusters, y: SDD+SSD clusters
      // trigger information
      TString                     fV0ReaderName;                          ///< Name of V0Reader
      TString                     CaloTriggerHelperName;                  ///< Name of CaloTriggerHelper 
      TString                     fCorrTaskSetting;                       ///< Name of Corr Task Setting
      AliVCaloTrigger*            fCaloTriggers;                          //!<! calo triggers
      TClonesArray*               fTriggerPatchInfo;                      //!<! trigger patch info array
      AliEMCALTriggerPatchInfo *  fMainTriggerPatchEMCAL;                 ///< main trigger patch, will be cached after first call
      TString                     fCaloTriggersName;                      ///< name of calo triggers collection
      TString                     fCaloTriggerPatchInfoName;              ///< trigger patch info array name
      ULong_t                     fTriggersEMCAL;                         ///< list of fired EMCAL triggers
      ULong_t                     fTriggersEMCALSelected;                 ///< list of accepted triggers
      Bool_t                      fEMCALTrigInitialized;                  ///< EMCAL triggers initialized
      TH1S*                       fHistoTriggThresh;                      ///< EMCal trigger thresholds
      Int_t                       fRunNumberTriggerOADB;                  ///< last used runnumber of OADB trigger object
      // Primary secondary distinction
      Double_t                    fSecProdBoundary;                       ///< 3D radius of production (cm) for primary-secodary distinction
      Float_t                     fMaxPtJetMC;                            ///< maximum jet pt in event
      Float_t                     fMinFacPtHard;                          ///< minimum factor between maximum jet pt and pt hard generated
      Float_t                     fMaxFacPtHard;                          ///< maximum factor between maximum jet pt and pt hard generated
      Float_t                     fMaxFacPtHardSingleParticle;            ///< maximum factor between maximum single particle pt (pi0/eta) and pt hard generated
      Int_t                       fMimicTrigger;                          ///< enable trigger mimiking
      TString                     fPathTriggerMimicSpecialInput;          ///< set special trigger mimiking OADB file
      Bool_t                      fRejectTriggerOverlap;                  ///< enable trigger overlap rejections
      //
      Bool_t                      fDoMultiplicityWeighting;               ///< Flag for multiplicity weighting
      TString                     fPathReweightingMult;                   ///< Path for file used in multiplicity reweighting
      TString                     fNameHistoReweightingMultData;          ///< Histogram name for reweighting Pi0
      TString                     fNameHistoReweightingMultMC;            ///< Histogram name for reweighting Eta
      TH1D*                       hReweightMultData;                      ///< histogram input for reweighting Eta
      TH1D*                       hReweightMultMC;                        ///< histogram input for reweighting Pi0
      phosTriggerType             fPHOSTrigger;                           // Kind of PHOS trigger: L0,L1
      Int_t                       fDebugLevel;                            ///< debug level for interactive debugging
  private:

      /// \cond CLASSIMP
      ClassDef(AliConvEventCuts,84)
      /// \endcond
};


#endif
