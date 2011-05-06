#ifndef ALIANALYSISTASKDIELECTRONREADAODBRANCH_H
#define ALIANALYSISTASKDIELECTRONREADAODBRANCH_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//  Class to retrieve the branch of the dielectron candidates stored in   //
//  filtered AOD files (AliAOD.Dielectron.root). It is possible to        //
//  apply tighter cuts to candidates stored in the branch and do the      //
//  matching with MC truth                                                //
//                                                                        //
//  Contacts : G. Bruno / Giuseppe.Bruno@ba.infn.it,                      //
//             F. Fionda / Fiorella.Fionda@ba.infn.it,                    //
//             C. Di Giglio / Carmelo.Digiglio@ba.infn.it                 // 
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"

class TNtuple;
class TH1F;
class TH2F;

class AliAnalysisTaskDielectronReadAODBranch : public AliAnalysisTaskSE
{
	public:
		AliAnalysisTaskDielectronReadAODBranch();
		AliAnalysisTaskDielectronReadAODBranch(const char *name);
		virtual ~AliAnalysisTaskDielectronReadAODBranch();

		// Implementation of interface methods
		virtual void UserCreateOutputObjects();
		virtual void Init();
		virtual void LocalInit() {Init();}
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *option);
 		Bool_t MatchToMC(Int_t pdgabs,const TObjArray *mcArray,const Int_t *dgLabels,Int_t ndg,Int_t ndgCk,const Int_t *pdgDg) const;
              
                //setters
                void SetPtLeg(Double_t cutPt){ fPtCut = cutPt;} 
                void SetSpdFirstRequired(Bool_t spdfirst){fSpdFirstRequired = spdfirst;}
                void SetNclsTPC(Int_t nCls){fClsTPC = nCls;}
                void SetPairType(Int_t type){fPairType = type;}
                void SetPtJpsi(Double_t ptjpsi) {fPtJpsi = ptjpsi;}
                void SetInvMassSignalRegion(Double_t lowLimit, Double_t upLimit){fInvMassSignalLimits[0]=lowLimit; fInvMassSignalLimits[1]=upLimit;}
                void SetInvMassSidebandRegion(Double_t lowLimit, Double_t upLimit){fInvMassSideBandsLimits[0]=lowLimit; fInvMassSideBandsLimits[1]=upLimit;}
                void SetHasMC(Bool_t mcFlag) {fHasMC = mcFlag;}


     enum {
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
      kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
      kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
      kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
      kHMPIDout=0x10000,kHMPIDpid=0x20000
  };

	private:
		AliAnalysisTaskDielectronReadAODBranch(const AliAnalysisTaskDielectronReadAODBranch &c);
		AliAnalysisTaskDielectronReadAODBranch& operator= (const AliAnalysisTaskDielectronReadAODBranch &c);

		TList      *fOutput;                 // output list with histo
	        TNtuple    *fNtupleJPSI;             // ntupla to store invariant mass and psproper of candidates
                TH1F       *fNentries;               // number of filtered events
		TH1F       *fInvMass;                // invariant mass of candidates in branch after tighter cuts application
		TH1F       *fInvMassNoCuts;          // invariant mass of candidates in branch without cuts
                TH1F       *fpsproperSignal;         // psproper distribution in the mass signal region           
		TH1F       *fpsproperSidebands;      // psproper distribution in the mass sidebands region           
		TH1F       *fpsproperAll;            // whole psproper distribution          
                TH1F       *fpsproperUnder;          // psproper distribution in the lower mass sidebands region 
                TH1F       *fpsproperUpper;          // psproper distribution in the upper mass sidebands region
		TH2F       *fLxyVsPtleg1;            // psproper vs pt for leg1
		TH2F       *fLxyVsPtleg2;            // psproper vs pt for leg2          
		TH2F       *fLxyVsPt;                // psproper vs pt jpsi 
		TH2F       *fMeeVsPt;                // invariant mass vs pt jpsi
		TH2F       *fMeeVsLxy;               // invariant mass vs psproper

		// QA plots
		TH1F       *fprimvtxZ;               // z coord of primary vertex        
		TH1F       *fsecvtxZ;                // z coord of secondary vertex        
		TH1F       *fprimvtxX;               // x coord of primary vertex          
		TH1F       *fsecvtxX;                // x coord of secondary vertex           
		TH1F       *fprimvtxY;               // y coord of primary vertex          
		TH1F       *fsecvtxY;                // y coord of secondary vertex          
		TH1F       *fPt;                     // pT(J/psi) distribution               
		TH1F       *fPtLeg1;                 // pT(leg1)  distribution. Warning: mixture of pos and neg tracks               
		TH1F       *fPtLeg2;                 // pT(leg2)  distribution. Warning: mixture of pos and neg tracks               
		TH2F       *fdEdxP;                  // dEdx vs P for legs
		Bool_t     fHasMC;                   // flag for read MC branch
                TObjArray  *fobj;                    // TObjArray with stored reconstructed candidates
		TObjArray  *fobjMC;                  // TObjArray with MC tracks
 
                // cuts on candidates
                Double_t  fPtCut;                    // ptLeg 
                Bool_t    fSpdFirstRequired;         // spd any/first 
                Int_t     fClsTPC;                   // #clsTPC
                Int_t     fPairType;                 // pair Type
                Double_t  fPtJpsi;                   // pt jpsi
                Double_t  *fInvMassSignalLimits;     // invariant mass signal region to extract psproper distribution
                Double_t  *fInvMassSideBandsLimits;  // invariant mass sideband region to extract psproper distribution 

		ClassDef(AliAnalysisTaskDielectronReadAODBranch,2); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};


#endif
