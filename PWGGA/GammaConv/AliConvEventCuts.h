#ifndef ALICONVEVENTCUTS_H
#define ALICONVEVENTCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Friederike Bock, Daniel Muehlheim

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVTrack.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "TF1.h"
#include "TObjArray.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"
#include "TRandom3.h"
#include "AliVCaloTrigger.h"

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
        kLHC14e2a,        //!< anchored LHC12[a-h] pass 1
        kLHC14e2b,        //!< anchored LHC12[a-h] pass 1
        kLHC14e2c,        //!< anchored LHC12[a-h] pass 1
        kLHC15h1,         //!< anchored LHC12[a-h] pass 2
        kLHC15h2,         //!< anchored LHC12[a-h] pass 2
        kLHC16c2,         //!< anchored LHC12[a-h] pass 2 - JJ
        kLHC16c2_plus,    //!< anchored LHC12[a-h] pass 2 - JJ - additional stat
        
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
        // MC's corresponding to 2015 data
        kLHC15g3a3,       //!< anchored LHC15f pass 1
        kLHC15g3a,        //!< anchored LHC15f pass 1
        kLHC15g3c2,       //!< anchored LHC15f pass 1
        kLHC15g3c3,       //!< anchored LHC15f pass 1
        kLHC15g3,         //!< anchored LHC15f pass 1
        kLHC16a2a,        //!< anchored LHC15h pass 1
        kLHC16a2b,        //!< anchored LHC15h pass 1
        kLHC16a2c,        //!< anchored LHC15h pass 1
        kLHC15l1a2,       //!< anchored LHC15n pass 1
        kLHC15l1b2,       //!< anchored LHC15n pass 1
        kLHC15k1a1,       //!< LHC15o low IR firstPhysics
        kLHC15k1a2,       //!< LHC15o low IR firstPhysics
        kLHC15k1a3,       //!< LHC15o low IR firstPhysics
        kLHC16j7,         //!< LHC15o low IR pass4
        kLHC16g1,         //!< anchored LHC15o pass1 - general purpose
        kLHC16g1a,        //!< anchored LHC15o pass1 - general purpose 0-10%
        kLHC16g1b,        //!< anchored LHC15o pass1 - general purpose 10-50%
        kLHC16g1c,        //!< anchored LHC15o pass1 - general purpose 50-90%
        kLHC16g2,         //!< anchored LHC15o pass1 - general purpose EPOS-LHC
        kLHC16g3,         //!< anchored LHC15o pass1 - general purpose DPMJET
        kLHC16h4,         //!< anchored LHC15o pass1 - injected signals 0-100%
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
        // MC upgrade
        kLHC13d19,        //!< upgrade 5.5TeV PbPb
        
        // 2016
        kLHC16kl,         //!< pp 13 TeV
	kLHC16d,          //!< pp 13 TeV
	kLHC16e,          //!< pp 13 TeV
	kLHC16f,          //!< pp 13 TeV
	kLHC16g,          //!< pp 13 TeV
	kLHC16h,          //!< pp 13 TeV
	kLHC16i,          //!< pp 13 TeV
	kLHC16j,          //!< pp 13 TeV
	kLHC16o,          //!< pp 13 TeV
	kLHC16p,          //!< pp 13 TeV        
        kLHC16q,          //!< pPb 5 TeV
	kLHC16r,          //!< pPb 8 TeV
	kLHC16s,          //!< pPb 8 TeV
	kLHC16t,          //!< pPb 5 TeV
        // MC's corresponding to 2016 data
        kLHC16j2a1,       //!< anchored LHC16k pass 1 - general purpose Pythia8
        kLHC16j2b1,       //!< anchored LHC16k pass 1 - general purpose EPOSLHC
        kLHC16j2a2,       //!< anchored LHC16l pass 1 - general purpose Pythia8
        kLHC16j2b2,       //!< anchored LHC16l pass 1 - general purpose EPOSLHC
	// General purpose pp-13TeV
	kLHC17f6,         //!< anchored LHC16d pass 1 - general purpose Pythia8  
	kLHC17f9,         //!< anchored LHC16e pass 1 - general purpose Pythia8  
	kLHC17d1,         //!< anchored LHC16f pass 1 - general purpose Pythia8 Nominal/LowB field 
	kLHC17d17,         //!< anchored LHC16g pass 1 - general purpose Pythia8  
	kLHC17f5,         //!< anchored LHC16h pass 1 - general purpose Pythia8  
	kLHC17d3,         //!< anchored LHC16i pass 1 - general purpose Pythia8  
	kLHC17e5,         //!< anchored LHC16j pass 1 - general purpose Pythia8  
	kLHC17d20a1,         //!< anchored LHC16k pass 1 - general purpose Pythia8  
	kLHC17d20a1_extra,   //!< anchored LHC16k pass 1 - general purpose Pythia8  
	kLHC17d20a2,         //!< anchored LHC16l pass 1 - general purpose Pythia8  
	kLHC17d20a2_extra,   //!< anchored LHC16l pass 1 - general purpose Pythia8  
	kLHC17d16,         //!< anchored LHC16o pass 1 - general purpose Pythia8  
	kLHC17d18,         //!< anchored LHC16p pass 1 - general purpose Pythia8  
	// Pythia +JJ pp-13TeV
	kLHC17f8a,         //!< anchored LHC16k,l pass 1 - Pythia8+JJ  
	kLHC17f8b,         //!< anchored LHC16f pass 1 - Pythia8+JJ Nominal field 
	kLHC17f8c,         //!< anchored LHC16g pass 1 - Pythia8+JJ  
	kLHC17f8d,         //!< anchored LHC16j pass 1 - Pythia8+JJ 
	kLHC17f8e,         //!< anchored LHC16o pass 1 - Pythia8+JJ   
	//General purpose- pPb
        kLHC17a2a,            //!< anchored LHC16qt pass 1 - general purpose EPOSLHC
        kLHC17a2a_fast,       //!< anchored LHC16qt pass 1 - general purpose EPOSLHC, fast only
        kLHC17a2a_cent,       //!< anchored LHC16qt pass 1 - general purpose EPOSLHC, CENT
        kLHC17a2a_cent_woSDD, //!< anchored LHC16qt pass 1 - general purpose EPOSLHC, CENT woSDD
        kLHC17a2b,            //!< anchored LHC16qt pass 1 - general purpose DPMJET
        kLHC17a2b_fast,       //!< anchored LHC16qt pass 1 - general purpose DPMJET,  fast only
        kLHC17a2b_cent,       //!< anchored LHC16qt pass 1 - general purpose DPMJET,  CENT
        kLHC17a2b_cent_woSDD, //!< anchored LHC16qt pass 1 - general purpose DPMJET,  CENT woSDD
        kLHC17a3a,            //!< anchored LHC16r pass 1 - general purpose EPOSLHC
        kLHC17a3a_fast,       //!< anchored LHC16r pass 1 - general purpose EPOSLHC, fast only
        kLHC17a3a_cent,       //!< anchored LHC16r pass 1 - general purpose EPOSLHC, CENT
        kLHC17a3a_cent_woSDD, //!< anchored LHC16r pass 1 - general purpose EPOSLHC, CENT woSDD
        kLHC17a3b,            //!< anchored LHC16r pass 1 - general purpose DPMJET
        kLHC17a3b_fast,       //!< anchored LHC16r pass 1 - general purpose DPMJET,  fast only
        kLHC17a3b_cent,       //!< anchored LHC16r pass 1 - general purpose DPMJET,  CENT
        kLHC17a3b_cent_woSDD, //!< anchored LHC16r pass 1 - general purpose DPMJET,  CENT woSDD
        kLHC17a4a,            //!< anchored LHC16s pass 1 - general purpose EPOSLHC
        kLHC17a4a_fast,       //!< anchored LHC16s pass 1 - general purpose EPOSLHC, fast only
        kLHC17a4a_cent,       //!< anchored LHC16s pass 1 - general purpose EPOSLHC, CENT
        kLHC17a4a_cent_woSDD, //!< anchored LHC16s pass 1 - general purpose EPOSLHC, CENT woSDD
        kLHC17a4b,            //!< anchored LHC16s pass 1 - general purpose DPMJET
        kLHC17a4b_fast,       //!< anchored LHC16s pass 1 - general purpose DPMJET,  fast only
        kLHC17a4b_cent,       //!< anchored LHC16s pass 1 - general purpose DPMJET,  CENT
        kLHC17a4b_cent_woSDD, //!< anchored LHC16s pass 1 - general purpose DPMJET,  CENT woSDD

        kLHC17f2a,            //!< anchored LHC16qt pass 1 - general purpose EPOSLHC
        kLHC17f2a_fast,       //!< anchored LHC16qt pass 1 - general purpose EPOSLHC, fast only
        kLHC17f2a_cent,       //!< anchored LHC16qt pass 1 - general purpose EPOSLHC, CENT
        kLHC17f2a_cent_woSDD, //!< anchored LHC16qt pass 1 - general purpose EPOSLHC, CENT woSDD
        kLHC17f2b,            //!< anchored LHC16qt pass 1 - general purpose DPMJET
        kLHC17f2b_fast,       //!< anchored LHC16qt pass 1 - general purpose DPMJET,  fast only
        kLHC17f2b_cent,       //!< anchored LHC16qt pass 1 - general purpose DPMJET,  CENT
        kLHC17f2b_cent_woSDD, //!< anchored LHC16qt pass 1 - general purpose DPMJET,  CENT woSDD
        kLHC17f3a,            //!< anchored LHC16r pass 1 - general purpose EPOSLHC
        kLHC17f3a_fast,       //!< anchored LHC16r pass 1 - general purpose EPOSLHC, fast only
        kLHC17f3a_cent,       //!< anchored LHC16r pass 1 - general purpose EPOSLHC, CENT
        kLHC17f3a_cent_woSDD, //!< anchored LHC16r pass 1 - general purpose EPOSLHC, CENT woSDD
        kLHC17f3b,            //!< anchored LHC16r pass 1 - general purpose DPMJET
        kLHC17f3b_fast,       //!< anchored LHC16r pass 1 - general purpose DPMJET,  fast only
        kLHC17f3b_cent,       //!< anchored LHC16r pass 1 - general purpose DPMJET,  CENT
        kLHC17f3b_cent_woSDD, //!< anchored LHC16r pass 1 - general purpose DPMJET,  CENT woSDD
        kLHC17f4a,            //!< anchored LHC16s pass 1 - general purpose EPOSLHC
        kLHC17f4a_fast,       //!< anchored LHC16s pass 1 - general purpose EPOSLHC, fast only
        kLHC17f4a_cent,       //!< anchored LHC16s pass 1 - general purpose EPOSLHC, CENT
        kLHC17f4a_cent_woSDD, //!< anchored LHC16s pass 1 - general purpose EPOSLHC, CENT woSDD
        kLHC17f4b,            //!< anchored LHC16s pass 1 - general purpose DPMJET
        kLHC17f4b_fast,       //!< anchored LHC16s pass 1 - general purpose DPMJET,  fast only
        kLHC17f4b_cent,       //!< anchored LHC16s pass 1 - general purpose DPMJET,  CENT
        kLHC17f4b_cent_woSDD, //!< anchored LHC16s pass 1 - general purpose DPMJET,  CENT woSDD

        // 
        kUnknownPeriod//!< kUnknownPeriod
      };

      /**
       * @enum EnergyVar
       * @brief Supported collision systems
       */
      enum EnergyVar {
        kUnset        = 0,  //!< not defined
        k900GeV       = 1,  //!< pp 900 GeV
        k2760GeV      = 2,  //!< pp 2.76TeV
        k5TeV         = 3,  //!< pp 5 TeV
        k7TeV         = 4,  //!< pp 7 TeV
        k8TeV         = 5,  //!< pp 8 TeV
        k13TeV        = 6,  //!< pp 13 TeV
        kpPb5TeV      = 7,  //!< pPb 5 TeV
        kpPb8TeV      = 8,  //!< pPb 8 TeV
        kPbPb2760GeV  = 9,  //!< PbPb 2.76TeV
        kPbPb5TeV     = 10, //!< PbPb 5 TeV
      };

      
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
      void    SetTriggerMimicking(Bool_t value)                                     { fMimicTrigger = value                                     ; 
                                                                                      if(value)AliInfo("enabled trigger mimicking")             ; }
      void    SetTriggerOverlapRejecion (Bool_t value)                              { fRejectTriggerOverlap = value                             ; 
                                                                                      if(value)AliInfo("enabled trigger overlap rejection")     ; }

      void    SetV0ReaderName (TString name)                                        { fV0ReaderName = name                                      ; }
      
      void    SetAddedSignalPDGCode (Int_t addedSignalPDGcode)                      { fAddedSignalPDGCode = addedSignalPDGcode                  ; }
      void    SetPreSelectionCutFlag (Bool_t preSelFlag)                            { fPreSelCut = preSelFlag                                   ; }   
      void    SetCaloTriggerPatchInfoName(const char *n)                            { fCaloTriggerPatchInfoName = n                             ; }
      void    SetCaloTriggersName(const char *n)                                    { fCaloTriggersName  = n                                    ; }
      void    SetAcceptedHeader(TList *HeaderList)                                  { fHeaderList = HeaderList                                  ; }   
      void    SetFillCutHistograms( TString name="",
                                    Bool_t preCut = kTRUE)                          { if(!fHistograms){ InitCutHistograms(name,preCut);}        ; }
      void    SetEtaShift(Double_t etaShift)                                        { fEtaShift = etaShift                                      ; } // Eta shift Setting
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

      void    SetMaxFacPtHard(Float_t value)                                        { fMaxFacPtHard = value                                     ; 
                                                                                      AliInfo(Form("maximum factor between pt hard and jet put to: %2.2f",fMaxFacPtHard));
                                                                                    }  
      void    SetDebugLevel( Int_t value)                                           { fDebugLevel = value                                       ; }
      
      // Geters
      TString   GetCutNumber();
      TString*  GetFoundHeader()                                                    { return fGeneratorNames                                    ; }
      Int_t     GetEventQuality()                                                   { return fEventQuality                                      ; }
      Bool_t    GetIsFromPileup()                                                   { return fRemovePileUp                                      ; }
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
      TString   GetSpecialTriggerName()                                             { return fSpecialTriggerName                                ; }
      AliEMCALTriggerPatchInfo   *GetMainTriggerPatch();
      ULong_t   GetTriggerList();
      Float_t   GetWeightForCentralityFlattening(AliVEvent *event = 0x0);
      Float_t   GetWeightForMultiplicity(Int_t mult);
      Float_t   GetWeightForMeson( Int_t index, AliMCEvent *mcEvent, AliVEvent *event = 0x0);
      Float_t   GetCentrality(AliVEvent *event);
      Bool_t    GetUseNewMultiplicityFramework(); 
      void      GetCorrectEtaShiftFromPeriod();
      void      GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliVEvent *event);
      TClonesArray*     GetArrayFromEvent(AliVEvent* event, const char *name, const char *clname=0);
      
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
      void    SetLightOutput( Bool_t flag ){fDoLightOutput = flag; return;}
      
      ///Cut functions
      Int_t   IsParticleFromBGEvent(  Int_t index, 
                                      AliMCEvent *mcEvent,
                                      AliVEvent *event = 0x0);
      
      void    LoadWeightingFlatCentralityFromFile ();
      void    LoadWeightingMultiplicityFromFile ();
      void    LoadReweightingHistosMCFromFile ();

      // Event Cuts
      Bool_t    IsCentralitySelected(AliVEvent *event, AliMCEvent *mcEvent);
      Bool_t    IsOutOfBunchPileupPastFuture(AliVEvent *event);
      Bool_t    IsPileUpV0MTPCout(AliVEvent *event);
      Bool_t    VertexZCut(AliVEvent *event);
      Bool_t    IsJetJetMCEventAccepted(AliMCEvent *event, Double_t& weight);
      Float_t   GetPtHard(AliMCEvent *event);
      void      GetXSectionAndNTrials(AliMCEvent *event, Float_t &XSection, Float_t &NTrials);
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

      Bool_t                      fDoLightOutput;                         ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
      Int_t                       fEventQuality;                          ///< EventQuality
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
      Bool_t                      fRemovePileUp;                          ///< flag
      Int_t                       fPastFutureRejectionLow;                ///< sets bunch crossing event rejection in past
      Int_t                       fPastFutureRejectionHigh;               ///< sets bunch crossing event rejection in future
      Int_t                       fDoPileUpRejectV0MTPCout;               ///< reject event if # TPCout tracks does not follow expected V=M mult
      TF1 *                       fFPileUpRejectV0MTPCout;                ///< Pol1 function to compute the cut
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

      TObjString*                 fCutString;                             ///< cut number used for analysis
      TString                     fCutStringRead;                         ///<
      AliAnalysisUtils*           fUtils;                                 ///<
      Double_t                    fEtaShift;                              ///<
      Bool_t                      fDoEtaShift;                            ///< Flag for Etashift
      Int_t                       fDoCentralityFlat;                      ///<
      TString                     fNameHistoNotFlatCentrality;            ///<
      TString                     fPathWeightsFlatCent;                   ///<
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
      Int_t                       fAddedSignalPDGCode;
      Bool_t                      fPreSelCut;                             // Flag for preselection cut used in V0Reader
      Bool_t                      fTriggerSelectedManually;               // Flag for manual trigger selection
      TString                     fSpecialTriggerName;                    // Name of the Special Triggers
      TString                     fSpecialSubTriggerName;                 // Name of the Special Triggers
      Int_t                       fNSpecialSubTriggerOptions;
      TH2F*                       hSPDClusterTrackletBackgroundBefore;    ///< SPD tracklets vs SPD clusters for background-correction before cut
      TH2F*                       hSPDClusterTrackletBackground;          ///< SPD tracklets vs SPD clusters for background-correction
      // trigger information
      TString                     fV0ReaderName;                          ///< Name of V0Reader
      AliVCaloTrigger*            fCaloTriggers;                          //!<! calo triggers
      TClonesArray*               fTriggerPatchInfo;                      //!<! trigger patch info array
      AliEMCALTriggerPatchInfo *  fMainTriggerPatchEMCAL;                 ///< main trigger patch, will be cached after first call
      TString                     fCaloTriggersName;                      ///< name of calo triggers collection
      TString                     fCaloTriggerPatchInfoName;              ///< trigger patch info array name
      ULong_t                     fTriggersEMCAL;                         ///< list of fired EMCAL triggers
      ULong_t                     fTriggersEMCALSelected;                 ///< list of accepted triggers
      Bool_t                      fEMCALTrigInitialized;                  ///< EMCAL triggers initialized
      // Primary secondary distinction
      Double_t                    fSecProdBoundary;                       ///< 3D radius of production (cm) for primary-secodary distinction
      Float_t                     fMaxPtJetMC;                            ///< maximum jet pt in event
      Float_t                     fMaxFacPtHard;                          ///< maximum factor between maximum jet pt and pt hard generated
      Float_t                     fMaxFacPtHardSingleParticle;            ///< maximum factor between maximum single particle pt (pi0/eta) and pt hard generated
      Bool_t                      fMimicTrigger;                          ///< enable trigger mimiking
      Bool_t                      fRejectTriggerOverlap;                  ///< enable trigger overlap rejections
      // 
      Bool_t                      fDoMultiplicityWeighting;               ///< Flag for multiplicity weighting
      TString                     fPathReweightingMult;                   ///< Path for file used in multiplicity reweighting
      TString                     fNameHistoReweightingMultData;          ///< Histogram name for reweighting Pi0
      TString                     fNameHistoReweightingMultMC;            ///< Histogram name for reweighting Eta
      TH1D*                       hReweightMultData;                      ///< histogram input for reweighting Eta
      TH1D*                       hReweightMultMC;                        ///< histogram input for reweighting Pi0
      Int_t                       fDebugLevel;                            ///< debug level for interactive debugging
  private:

      /// \cond CLASSIMP
      ClassDef(AliConvEventCuts,34)
      /// \endcond
};


#endif
