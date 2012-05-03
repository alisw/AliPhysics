#ifndef ALIANALYSISTASKPARTONDISC_H
#define ALIANALYSISTASKPARTONDISC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 

//////////////////////////////////////////////////////
//                                                  //
// Analysis task for parton discrimination studies  //
//                                                  //
//////////////////////////////////////////////////////

class TH1F;
class TH1I;
class TH2F;
class TH2I;
class TH3F;
class TProfile;
class AliAODEvent;
class AliMCEvent;
class AliAODJet;
class AliAODTrack;
class AliAODMCParticle;
class AliAODVertex;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPartonDisc : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPartonDisc();
  AliAnalysisTaskPartonDisc(const char *name);
  virtual ~AliAnalysisTaskPartonDisc() {}
  
  virtual Bool_t UserNotify();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t *);
  
  virtual void SetAODwithMC(Bool_t flag)     {fUseAODMC = flag;}
  virtual void SetMCBranch(const char* mc)   {fBranchMC = mc;}
  virtual void SetRecBranch(const char* rec) {fBranchRec = rec;}
  virtual void SetSecondRecBranch(const char* secrec) {fBranchSecRec = secrec;}
  virtual void SetSqrtS(const Double_t sqrts) {fSqrts = sqrts;}
  virtual void SetXNtX(const Int_t x) {fNtX = x;}
  virtual void SetJetRadius(const Double_t jetradius) {fJetRadius = jetradius;}
  virtual void SetFilterBitTracks(const UInt_t bval) {fFilterBit = bval;}
  virtual void SetFlavorRadius(const Double_t fradius) {fFlavorRadius = fradius;}
  virtual void SetPhojetMC(Bool_t flagmc) {fPhojetMC = flagmc;}
  Int_t        GetMCEventType(AliMCEvent *mcEvent);
  Int_t        GetPhojetEventType(AliMCEvent *mcEvent);
  Bool_t       IsInsideAcceptance(AliAODJet *jet);
  Int_t        GetJetFlavour(AliAODJet *jet, Int_t ntracks,TClonesArray *mcarray);
  Double_t     GetDeltaR(Double_t eta1, Double_t phi1,Double_t eta2, Double_t phi2);
  Int_t        GetNumberOfMcChargedTracks(Int_t percentage,AliAODJet *Jet, Int_t ntracks, TClonesArray *mcarray, Double_t jr);
  Int_t        GetNumberOfChargedTracks(Int_t percentage,AliAODJet *jet, Int_t ntracks, AliAODEvent *aode, Double_t jr);
  virtual void AllocateStaticContainer(Int_t size);
  virtual void InitializeStaticContainer(Int_t size);
  virtual void SortArray(Double_t *pointer, Int_t arraySize);
  Int_t        TracksForPercentage(Double_t *array, Int_t arraysize, Int_t percentage, Double_t jetenergy);
  Bool_t       IsMCTrackInsideThisJet(AliAODMCParticle *MCParticle, AliAODJet *Jet, Double_t jr);
  Bool_t       IsTrackInsideThisJet(AliAODTrack *aodT, AliAODJet *Jet, Double_t jr);
  Bool_t       VertexInJet(AliAODVertex *pvtx, AliAODVertex *vtx, AliAODJet *jet, Double_t jr);
  Double_t     GetEtaValue(Double_t theta) const;
  Double_t     GetThetaAngle(Double_t xval, Double_t yval, Double_t zval);
  Double_t     GetPhiAngle(Double_t xval, Double_t yval);
  virtual void SetAODMCInput(Bool_t b){fUseAODJetInput = b;}
  Double_t     DeltaPhiMC(AliAODJet *jet, AliAODMCParticle *particle);
  Double_t     DeltaEtaMC(AliAODJet *jet, AliAODMCParticle *particle);
  Double_t     DeltaPhiSqMC(AliAODJet *jet, AliAODMCParticle *particle);
  Double_t     DeltaEtaSqMC(AliAODJet *jet, AliAODMCParticle *particle);
  Double_t     DeltaPhiTrack(AliAODJet *jet, AliAODTrack *track);
  Double_t     DeltaEtaTrack(AliAODJet *jet, AliAODTrack *track);
  Double_t     DeltaPhiSqTrack(AliAODJet *jet, AliAODTrack *track);
  Double_t     DeltaEtaSqTrack(AliAODJet *jet, AliAODTrack *track);
  virtual void SetMinPtTrackCut(const Double_t minptval) {fMinpTVal = minptval;}
  static Bool_t NumberOfReadEventsAOD(const char* currFile, Int_t &fNEvents);
  virtual void  SetOnlyMC(Bool_t flagOnlyMC)     {fUseOnlyMC = flagOnlyMC;}
  virtual void  SetCheckMCStatus(Bool_t flagMCStatus)     {fCheckMCStatus = flagMCStatus;}
  virtual void  HasOverlapedCones(TClonesArray *JetArray);
  virtual void  ResetJetFlags();
  virtual void  SetEnablePrints(Bool_t flagEnablePrints) {fEnablePrints = flagEnablePrints;}
  Bool_t        HasPerpendicularCone() const {return fHasPerpCone;}
  virtual void  SetHIEvent(Bool_t flagHIEvent)     {fIsHIevent = flagHIEvent;}
  Int_t         GetNMcChargedTracksAboveThreshold(AliAODJet *jet, Int_t ntracks, TClonesArray *mcarray, Double_t jr);
  Int_t         GetRecalcNTXMc(Int_t percentage, AliAODJet *originaljet, Int_t ntracks, TClonesArray *mcarray, Double_t jr);
  Int_t         GetRecalcNMcChTrUpThr(AliAODJet *jet, Int_t ntracks, TClonesArray *mcarray, Double_t jr);
  Int_t         GetNRecChargedTracksAboveThreshold(AliAODJet *jet, Int_t ntracks, AliAODEvent *aode, Double_t jr);
  Int_t         GetRecalcNTXRec(Int_t percentage,AliAODJet *originaljet, Int_t ntracks, AliAODEvent *aode, Double_t jr);
  Int_t         GetRecalcNRecChTrUpThr(AliAODJet *jet, Int_t ntracks, AliAODEvent *aode, Double_t jr);
  Int_t         TracksForPercentageRecalc(Double_t *array, Int_t arraysize, Int_t percentage, Double_t jetenergy);
  Bool_t        IsTrackInsideExcludedArea(Double_t tracketa, Double_t trackphi, TClonesArray *recojets);
  Double_t      GetV0ExcludedMultiplicity(TClonesArray *recojets);
  virtual void  SetMinPtUE(const Double_t minptvalUE) {fMinpTValUE = minptvalUE;}
  virtual void  SetMaxPtUE(const Double_t maxptvalUE) {fMaxpTValUE = maxptvalUE;}
  virtual void  SetMinPtMC(const Double_t minptvalMC) {fMinpTValMC = minptvalMC;}
  Int_t         GetV0LikeExcludedMultMC(TClonesArray *mcjets, TClonesArray *mcparticles);
  virtual void  SetIncreaseOfExclusionR(const Double_t increaseExclR) {fIncExcR = increaseExclR;}
  virtual void  ForceNotUseTrackRefs(const Bool_t flagForce) {fForceNotTR = flagForce;};
  virtual void  NotExtendDiJetExclusion(const Bool_t flagNotDiJ) {fNotExtDiJEx = flagNotDiJ;};
  virtual void  FillPerpConeHisto(TH3F *currenthisto, Int_t ntracks, AliAODEvent *aode, Int_t CentralityBin, Int_t pTBin);
  virtual void  ForceSkipSingleTrackJets(const Bool_t flagForceSJ) {fForceSkipSJ = flagForceSJ;};
  Bool_t        IsEqualRel(Double_t vA, Double_t vB);
  virtual void  SetEnableJetEtaRestriction(Bool_t flagEnableJetEtaRes) {fIncreasingExcl = flagEnableJetEtaRes;}
  virtual void  SetTrackRandomRejectionPerc(const Double_t perctrackrr) {fTTrackRandomRejection = perctrackrr;}
  virtual void  SetTrackInJetRandomRejectionPerc(const Double_t perctrackijrr) {fJTrackRandomRejection = perctrackijrr;}
  virtual void  SetMinPtCutGlobMult(const Double_t minptglobmult) {fMinPtInGlobMult = minptglobmult;}

 private:
  AliAODEvent *fAOD;         //! AOD object
  Bool_t      fUseAODMC;     // Flag for MC info in the AOD
  Bool_t      fPhojetMC;     // Flag for Phojet MC
  TString     fBranchMC;     // AOD branch name for MC jets
  TString     fBranchRec;    // AOD branch name for reconstructed jets
  TString     fBranchSecRec; // AOD branch name for secondary reconstructed jets 
  Double_t    fSqrts;        // Value of sqrt{s}  default 0 to spot errors
  Int_t       fNtX;          // X Value of NTX  default 0 to spot errors
  Double_t    fJetRadius;    // Radius used in jet finding default 0 to spot errors
  Double_t    fFlavorRadius; // Radius used in flavor asignment
  UInt_t      fFilterBit;    // Filterbit value: 16= tracks with standard cuts, 128 = tracks with hit in SDD, 144 (16+128) = all the previous
  TList       *fOutputList;  //! Output list
  TH1F        *fJetPt;      //! Pt spectrum of reco jets
  TH1F        *fJetPtSec;   //! Pt spectrum of secondary reco jets
  TH1F        *fJetPtMC;    //! Pt spectrum of mc jets
  TH2F        *fJetEta;     //! Eta of reco jets
  TH2F        *fJetEtaSec;  //! Eta of secondary reco jets
  TH2F        *fJetPhi;     //! Phi of reco jets
  TH2F        *fJetPhiSec;  //! Phi of secondary reco jets
  TH2F        *fJetEtaMC;   //! Eta of MC jets
  TH2F        *fJetPhiMC;   //! Phi of MC jets
  TH2F        *fPtAODMC;    //! Pt spectrum of MC tracks in AOD 
  TH2F        *fPtAOD;      //! Pt spectrum of tracks in AOD (reco MC or real)
  TH2F        *fEtaAODMC;   //! Eta distribution of MC AOD tracks 
  TH2F        *fPhiAODMC;   //! Phi distribution of MC AOD tracks 
  TH2F        *fEtaAOD;     //! Eta distribution of AOD tracks
  TH2F        *fPhiAOD;     //! Phi distribution of AOD tracks
  TH2F        *fFlavor;     //! Flavor distribution of jets
  TH2F        *fNJetsMC;    //! Number of jets per event in MC
  TH2F        *fNJetsRD;    //! Number of jets per event in real data or reco MC
  TH2F        *fNJetsRDSeco;     //! Number of jets per event in real data or reco MC secondary branch
  TH2F        *fJetsMultPtMC;    //! Jet multiplicity in jet pT in MC
  TH2F        *fJetsMultPtRD;    //! Jet multiplicity in jet pT in real data or reco MC
  static Double_t *fgContainer;     //! static container for track counting
  TH2F        *fNChTr[12];          //! Number of charged tracks in the jets as a function of jet pt (MC)
  TH2F        *fNChTrRD;            //! Number of charged tracks in the jets as a function of jet pt (Real Data)
  TProfile    *fProfNChTrRD;        //! Number of charged tracks in the jets as a function of jet pt (Real Data)
  TH2I        *fProcessPDG[6];      //! Pythia process and pT of the jet
  TH1F        *fHistPtParton[12];   //! Pt distribution of jets per flavor, mc and reco
  TH2F        *fFragPion[6];        //! Fragmentation of jet in pions, jet energy
  TH2F        *fFragKaon[6];        //! Fragmentation of jet in kaons, jet energy
  TH2F        *fFragProton[6];      //! Fragmentation of jet in protons, jet energy
  TH2F        *fHistContainerR4[6];    //! Temporary containers for fragmentation of reco jets R4
  TH2F        *fHistContainerR3[6];    //! Temporary containers for fragmentation of reco jets R3
  TH2F        *fHistContainerR2[6];    //! Temporary containers for fragmentation of reco jets R2
  TH2F        *fFragChargedR4[6];   //! Fragmentation of jet in charged part, jet energy R4
  TH2F        *fFragChargedR3[6];   //! Fragmentation of jet in charged part, jet energy R3
  TH2F        *fFragChargedR2[6];   //! Fragmentation of jet in charged part, jet energy R2
  TH2F        *fFragCandidates[2];  //! Tagged candidates for FF
  TH1F        *fFracQQ;             //! process qq as a function of x_{T} 
  TH1F        *fFracGQ;             //! process gq as a function of x_{T} 
  TH1F        *fFracGG;             //! process gg as a function of x_{T}
  TH1F        *fFracOutGoingQQ;     //! process outgoing qq as a function of x_{T} 
  TH1F        *fFracOutGoingGQ;     //! process outgoing gq as a function of x_{T} 
  TH1F        *fFracOutGoingGG;     //! process outgoing gg as a function of x_{T}
  TProfile    *fh1Xsec;             //! xsection from pyxsec.root
  TH1F        *fh1Trials;           //! ntrials from pyxsec.root            
  Short_t      fMpdg;               //! PDG code of mother of parton
  TH2F        *fProcessJetPt;       //! Pythia Process and jet pT 
  TH2F        *fFlavorLead;         //! Flavor distribution of the leading jet
  TH2F        *fProcessLeadJetPt;   //! Pythia Process and jet pT of the leading jet
  TH3F        *fPDGMothLPart;       //! PDG code of the mother of the leading parton, leading parton, jet pT
  TH2F        *fFlavProc;           //! Flavor, Flavor status code
  Float_t      fAvgTrials;          //  Average number of trials
  Bool_t       fUseAODJetInput;     //  take MC from input AOD not from output AOD
  Double_t     fMinTrackPtInNTX;    //  Minimum track pT taken into the NTX calculation
  Double_t     fMaxTrackPtInNTX;    //  Maximum track pT taken into the NTX calculation
  TH3F        *fMinTrackPtInNTXh[2];//! Histo to save fMinTrackPtInNTX as a function of jet pT
  TH2F        *fMaxTrackPtInNTXh[2];//! Histo to save fMaxTrackPtInNTX as a function of jet pT
  TH2F        *fSCM[12];            //! Second central moment as a function of jet pt (MC)
  TH2F        *fSCMRD;              //! Second central moment as a function of jet pt (Real Data)
  Double_t     fMinpTVal;           //  Minimum pT track cut for SCM analysis
  TH2F        *fZVertex;            //! Z coordinate vertex position, number of reco jets
  TH1F        *fh1Events;           //! nevents read out from PWG4_JetTasksOutput.root
  Bool_t       fUseOnlyMC;          //  Flag to signal only MC input
  Bool_t       fCheckMCStatus;      //  Flag to check the status of MC, not working for old aliroot
  Bool_t       fJetFlags[16];       //  Flag to mark if this jet is ok (acceptance and no overlap), kTRUE if usable
  Int_t        fEvtCount;           //  Event counter for debugging
  TH2F        *fNAccJetsMC;         //! Number of accepted jets per event in MC
  TH2F        *fNAccJetsRD;         //! Number of accepted jets per event in real data or reco MC
  TH2F        *fNAccJetsRDSeco;     //! Number of jets accepted per event in real data or reco MC secondary branch                    
  Bool_t       fEnablePrints;       //  Flag to enable print outs
  TH1F        *fRecJetPtInclusive;  //! Pt spectrum of inclusive reco jets
  TH1F        *fMCJetPtInclusive;   //! Pt spectrum of inclusive MC jets
  TH1F        *fRecJetPtLeading;    //! Pt spectrum of leading reco jets
  TH1F        *fMCJetPtLeading;     //! Pt spectrum of leading MC jets
  TH1F        *fSecRecJetPtInclusive;  //! Pt spectrum of inclusive reco jets (2nd branch)
  TH1F        *fSecRecJetPtLeading;    //! Pt spectrum of leading reco jets (2nd branch)
  Bool_t       fHasPerpCone;           // Flag to indicate if possible to use a perpendicular jet for bckg
  Double_t     fEtaPerpCoord;          // Eta of perpendicular cone
  Double_t     fPhiPerpCoord;          // Phi of perpendicular cone
  Double_t     fPtPerpCoord;           // pT of perpendicular cone
  Bool_t       fJetEvent;              // Flag to indicate a jet event(in acceptance)
  AliAODJet   *fPerpCone;              // Perpendicular Cone
  TH2F        *fNChTrMCPerp;           //! Number of charged tracks in the perpendicular cone MC
  TH2F        *fNChTrRecPerp;          //! Number of charged tracks in the perpendicular cone Reco or Real
  TH2F        *fSCMMCPerp;             //! Second central moment as a function of jet pt for the perp cone MC
  TH2F        *fSCMRecPerp;            //! Second central moment as a function of jet pt for the perp cone Rec
  Bool_t       fIsHIevent;             // Flag to indicate that is reading a HI event
  Double_t     fCurrentJetMinPtNT90;   // Min pT used in the NT90 calculation of the current jet
  Double_t    *fBckgSbsJet;            //! Current jet, background substracted
  Double_t     fCurrentJetMinPtNT90Recalc;  // Min pT used in the NT90 re-calculation of the current jet 
                                            // after energy correction of the jet
  TH2F        *fNChTrCorrMCQuark;           //! Number of charged tracks after jet energy correction, as a function of corrected pT, MC quarks
  TH2F        *fNChTrCorrMCGluon;           //! Number of charged tracks after jet energy correction, as a function of corrected pT, MC gluons
  TH2F        *fNChTrCorrMCPerp;            //! Number of charged tracks in the perpendicular cone MC, after correction
  Bool_t       fIsPossibleToSubstBckg;      //  Flag to signal that there was a perpendicular cone, as is possible to substract background
  TH3F        *fNChTrRecECorr;              //! Number of charged tracks in the energy corrected jet, as a function of corrected jet pt, centrality
  TH3F        *fNChTrRecPerpECorr;          //! Number of charged tracks in the perpendicular after correction, as a func. of corr. jet pt, centrality
  TH1F        *fRefMult;                    //! Reference multiplicity in the AOD
  TH2F        *fNChTrRDMult[8];             //! Number of charged tracks in the jets as a function of jet pt (Real Data), for reference multiplicities in pp
  TH1F        *fNAccJetsRDMult[8];          //! Number of accepted jets per event in real data or reco MC, for reference multiplicities in pp
  Int_t        fCurrentJetCharge;           //  Charge in the current jet
  TH1F        *fTotalJetCharge[8];          //! Charge of this jet (dependent on event multiplicity)
  TH2F        *fRefMultWOJet;               //! Reference multiplicity in the AOD and multiplicity without jets
  TH2F        *fSCMRDMult[8];               //! Second central moment as a function of jet pt (Real Data), for reference multiplicities in pp
  TH2F        *fNChTrRDMultMC[8];           //! Number of charged tracks in the jets as a function of jet pt (MC Data), for reference multiplicities in pp
  TH2F        *fSCMRDMultMC[8];             //! Second central moment as a function of jet pt (MC Data), for reference multiplicities in pp
  TH2F        *fVZEROMult;                  //! Multiplicity in V0A and V0C
  TH2F        *fMultWOJetVZero;             //! Multiplicity without jets, and VZERO multiplicity
  AliAODVZERO *fVZero;                      //! AOD VZERO object
  TH2F        *fNChTrRDMultSE[8];           //! Number of charged tracks in the jets as a function of jet pt (Real Data), for 2nd reference multiplicities in pp
  TH1F        *fNAccJetsRDMultSE[8];        //! Number of accepted jets per event in real data or reco MC, for 2nd reference multiplicities in pp
  TH1F        *fTotalJetChargeSE[8];        //! Charge of this jet (dependent on 2nd event multiplicity)
  TH2F        *fSCMRDMultSE[8];             //! Second central moment as a function of jet pt (Real Data), for 2nd reference multiplicities in pp
  TH2F        *fRefMultFullV0;              //! Reference multiplicity in the AOD and multiplicity from the full V0
  TH2F        *fRefMultV0Corr;              //! Reference multiplicity in the AOD and multiplicity from the V0 sectors with no jets (2 leading)
  TH2F        *fFullV0V0Corr;               //! Multiplicity from the full V0, Multiplicity from the V0 sectors with no jets (2 leading)
  TH3F        *fNTXV0MultPt;                //! NTX, V0 corrected multiplicity, and jet pT
  TH3F        *fNTXCBMultPt;                //! NTX, Central Barrel corrected multiplicity, and jet pT
  TH2F        *fNChTrRDMultOJ[8];           //! Number of charged tracks in the jets as a function of jet pt (Real Data), for reference multiplicities in pp, 1 Jet
  TH2F        *fSCMRDMultOJ[8];             //! Second central moment as a function of jet pt (Real Data), for reference multiplicities in pp, 1 Jet
  TH2F        *fNChTrRDMultSEOJ[8];         //! Number of charged tracks in the jets as a function of jet pt (Real Data), for 2nd reference multiplicities in pp, 1 Jet
  TH2F        *fSCMRDMultSEOJ[8];           //! Second central moment as a function of jet pt (Real Data), for 2nd reference multiplicities in pp, 1 Jet
  Double_t     fMinpTValUE;                 //  Minimum pT track cut for the UE multiplicity, default 2 GeV for debugging
  TH2F        *fRefMultFullV0UJ;            //! Reference multiplicity in the AOD and multiplicity from the full V0, un jet
  TH2F        *fRefMultV0CorrUJ;            //! Reference multiplicity in the AOD and multiplicity from the V0 sectors with no jets (2 leading), un jet
  TH2F        *fFullV0V0CorrUJ;             //! Multiplicity from the full V0, Multiplicity from the V0 sectors with no jets (2 leading), un jet
  TH2F        *fMultWOJetVZeroUJ;           //! Multiplicity without jets, and VZERO multiplicity, un jet
  TH2F        *fRefMultWOJetUJ;             //! Reference multiplicity in the AOD and multiplicity without jets, un jet
  Double_t     fMaxpTValUE;                 //  Maximum pT track cut for the UE multiplicity, default 2 GeV for debugging
  TH2F        *fRefAODTrackCount;           //! Correlation between ref aod mult. and my own counting
  TH2F        *fRefAODTrackCountUJ;         //! Correlation between ref aod mult. and my own counting, single jet event
  TH2F        *fTrackCountWOJet;            //! Correlation between my own counting TPC & soft TPC
  TH2F        *fTrackCountWOJetUJ;          //! Correlation between my own counting TPC & soft TPC UJ
  TH2F        *fTrackCountWOJetUJMC;        //! Correlation between my own counting TPC & soft TPC UJ solo MC
  TH2F        *fFullV0V0CorrUJMC;           //! Multiplicity from the full V0, Multiplicity from the V0 sectors with no jets, un jet MC no real V0
  TH2F        *fNChTrRDMultOJMC[8];         //! Number of charged tracks in the jets as a function of jet pt (MC), for reference multiplicities in pp, 1 Jet
  TH2F        *fSCMRDMultOJMC[8];           //! Second central moment as a function of jet pt (MC), for reference multiplicities in pp, 1 Jet
  TH2F        *fNChTrRDMultSEOJMC[8];       //! Number of charged tracks in the jets as a function of jet pt (MC), for 2nd reference multiplicities in pp, 1 Jet
  TH2F        *fSCMRDMultSEOJMC[8];         //! Second central moment as a function of jet pt (MC), for 2nd reference multiplicities in pp, 1 Jet
  Double_t     fMinpTValMC;                 //  Minimum pT track cut for the MC multiplicity, default 2 GeV for debugging
  Double_t     fIncExcR;                    //  Increase in the exclusion radius value
  Bool_t       fForceNotTR;                 //  Force NOT to use track references
  Bool_t       fNotExtDiJEx;                //  Not extend the exclusion in the dijet area, old behaviour
  TH3F        *fMinTrackPtInNTXRecalc;      //! Histo to save fMinTrackPtInNTX after recalculation as a function of jet pT
  TH2F        *fMaxTrackPtInNTXRecalc;      //! Histo to save fMaxTrackPtInNTX after recalculation as a function of jet pT
  TH3F        *fPtDistInJetConeRaw;         //! pT distributions of tracks inside the cone for jet pT (raw) ranges and centralities
  TH3F        *fPtDistInPerpConeRaw;        //! pT distributions of tracks inside the perpendicular cone for jet pT (raw) ranges and centralities
  TH3F        *fPtInPerpCon;                //! summed pT from the perpendicular cone for jet pT (raw) ranges and centralities
  Double_t     fMinTrackPtInNTXR;           //  Minimum track pT taken into the NTX re-calculation
  Double_t     fMaxTrackPtInNTXR;           //  Maximum track pT taken into the NTX re-calculation
  Double_t     fEventCent;                  //  event centrality
  TH2F        *fNChTrRecPerpMultSEOJ[8];    //! Number of charged tracks in the perpendicular cone reco or real, with multiplicities
  TH1F        *fJetEtaAll;                  //! Eta distribution of all the found jets, no cuts
  TH1F        *fJetEtaOnlyTPCcut;           //! Eta distribution of all the found jets, only with eta acceptance cut
  TH1F        *fJetEtaJetPt[3];             //! Eta distribution of analyzed jets, 3 ranges of pT
  TH3F        *fNChTrRecECorrPPMult;        //! Number of charged tracks in the energy corrected jet, as a function of corrected jet pt, pp mult
  TH3F        *fNChTrRecPerpECorrPPMult;    //! Number of charged tracks in the perpendicular after correction, as a func. of corr. jet pt, pp mult
  Bool_t       fForceSkipSJ;                //  Force to skip single track jets
  TH2F        *fJetPtCentPbPbRaw;           //! Raw pT spectrum of reco jets, centrality in PbPb
  TH2F        *fJetPtCentPbPbCorr;          //! Corrected pT spectrum of reco jets, centrality in PbPb
  Double_t     fJetAcceptance;              //  Acceptance cut on jets, for multiplicity in PbPb
  Bool_t       fIncreasingExcl;             //  Flag to indicate that the analyis increases exclusion beyond jet radius
  TH3F        *fTotTracksCone;              //! total number of tracks in the jet cone, for jet pT (raw) ranges and centralities
  Int_t        fTotTracksInCone;            //  Total number of tracks in the jet cone

  Double_t     fTTrackRandomRejection;      //  Percentage of tracks from the event randomly rejected
  Double_t     fJTrackRandomRejection;      //  Percentage of tracks from the jet randomly rejected
  TH1F        *fJEtaMCMultOJ[8];            //! Eta distribution of jets as a function of jet pt (MC), for V0-like multiplicities in pp, 1 Jet
  TH1F        *fJEtaMCMultSEOJ[8];          //! Eta distribution of jets as a function of jet pt (MC), for TPC-like multiplicities in pp, 1 Jet
  TH1F        *fJEtaRDMultOJ[8];            //! Eta distribution of jets as a function of jet pt (Reco Data), for V0 multiplicities in pp, 1 Jet
  TH1F        *fJEtaRDMultSEOJ[8];          //! Eta distribution of jets as a function of jet pt (Reco Data), for TPC-like multiplicities in pp, 1 Jet
  TH1F        *fJetPtMCMultOJ[8];           //! Pt spectrum jets (MC), for V0-like multiplicities in pp, 1 jet
  TH1F        *fJetPtMCMultSEOJ[8];         //! Pt spectrum jets (MC), for TPC-like multiplicities in pp, 1 jet
  TH1F        *fJetPtRDMultOJ[8];           //! Pt spectrum jets (Reco Data), for V0 multiplicities in pp, 1 jet
  TH1F        *fJetPtRDMultSEOJ[8];         //! Pt spectrum jets (Reco Data), for TPC multiplicities in pp, 1 jet
  TH2F        *fEntriesQuark[8];            //! Quark NT90 in MC in the multiplicity bins
  TH2F        *fEntriesGluon[8];            //! Gluon NT90 in MC in the multiplicity bins
  Double_t     fMinPtInGlobMult;            //  Min pT used in the global multiplicity calculation

  AliAnalysisTaskPartonDisc(const AliAnalysisTaskPartonDisc&); // not implemented
  AliAnalysisTaskPartonDisc& operator=(const AliAnalysisTaskPartonDisc&); // not implemented
  
  ClassDef(AliAnalysisTaskPartonDisc, 4); 
};

#endif
