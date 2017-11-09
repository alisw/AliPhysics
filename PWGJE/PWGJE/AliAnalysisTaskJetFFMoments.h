#ifndef ALIANALYSISTASKJETFFMOMENTS_H
#define ALIANALYSISTASKJETFFMOMENTS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// task used for comparing different jets D parmaters from fastjet 
// *******************************************

#include  "AliAnalysisTaskSE.h"
#include  "AliAnalysisManager.h"
#include <cmath>

#ifndef __CINT__
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/CompositeJetStructure.hh"
#include "fastjet/Selector.hh"
#else
namespace fastjet {
  class PseudoJet;
  class ClusterSequenceArea;
  class GhostedAreaSpec;
  class AreaType;
  class JetDefinition;
  class JetAlgorithm;
  class CompositeJetStructure;
  class Selector;
  class SelectorWorker;
  class Strategy;
  class RecombinationScheme;
}
#endif

using std::vector;

/////////////////
class TList;
class TClonesArray;
class TRefArray;
class TH2F;
class THnSparse;
class TProfile;
class TProfile2D;
class TF1;

class AliAODEvent;
class AliAODExtension;
class AliAODJet;
class AliAODTrack;
class AliAODMCParticle;
class AliAODJetEventBackground;
class AliAnalysisHelperJetTasks;
#include "AliGenPythiaEventHeader.h"
#include "AliGenHerwigEventHeader.h"

class AliAnalysisTaskJetFFMoments : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskJetFFMoments();
  AliAnalysisTaskJetFFMoments(const char* name);
  virtual ~AliAnalysisTaskJetFFMoments();

  //
  // Implementation of analysis framework
  //

  void UserCreateOutputObjects();
  void LocalInit();
  void UserExec(Option_t *option);
  void Terminate(Option_t *option);
  Bool_t Notify();

  //
  // Setters / Getters
  //

  // Event
  void SetAODTrackInput(Bool_t b)                      {fkUseAODTrackInput = b;}  			   // Use input (kTrue) AOD for tracks  
  void SetAODMCInput(Bool_t b)                         {fkUseAODMCInput = b;}        			   // Use input AOD (kTrue) for MC tracks
  void SetFilterPt(Float_t f)                          {fFilterPt = f;}                                    // Min jet pt in output AOD
  void SetEventSelection(Bool_t b)                     {fkEventSelection = b;}   			   // Bypass event selection?
  void SetIsPbPb()                                     {fkIsPbPb=kTRUE;} 			 	   // PbPb case
  void SetVtxCuts(Float_t z,Float_t r = 1)             {fVtxZMax = z; fVtxR2Max = r*r;} 		   // Set vertex cuts in z and x*y    
  void SetCentralityCut(Float_t xLo,Float_t xUp)       {fCentCutLo = xLo; fCentCutUp = xUp;} 		   // Global centrality range
  void SetCentralityClasses(Int_t bin0 = 10, Int_t bin1 = 30, Int_t bin2 = 50, Int_t bin3 = 80) 
  {fCentClass[0] = bin0; fCentClass[1] = bin1; fCentClass[2] = bin2; fCentClass[3] = bin3;} 	           // Set centrality upper bin limits
  // Setters for detector level effects
  void SetRequireT0vtx(Bool_t b = true)                {fkRequireTZEROvtx = b;} 			   // Set to require T0 vtx 
  void SetRequireV0AC(Bool_t b = true)                 {fkRequireVZEROAC = b;} 			           // Set to require V0 AC 
  void SetRejectPileup(Bool_t b = true)                {fkRejectPileup = b;}                               // Set to reject pileup    
  void SetRejectFastCluster(Bool_t b = false)          {fkRejectFastOnly = b;}                             // Reject fast cluster
  void SetPtHardCuts(Double_t jetpt = 0, Double_t trackpt = 0) {fPtHardAndPythiaJetPtFactor = jetpt; fPtHardAndTrackPtFactor = trackpt;}

  // Tracks
  void SetFilterMask(UInt_t i,Int_t iType = 0)         {fFilterMask = i; fFilterType = iType;} 	            // Set filter mask and types
  void SetFilterMaskBestPt(UInt_t i)                   {fFilterMaskBestPt = i;} 			   //  
  void SetTrackTypeGen(Int_t i)                        {fTrackType[0] = i;} 				   // Set generated track type 
  void SetTrackTypeRec(Int_t i)                        {fTrackType[1] = i;} 				   // Set reconstructed track type
  void SetTrackPtMin(Float_t x)                        {fTrackPtMin = x;} 				   // Set track min pt cut
  void SetTrackPtMax(Float_t x)                        {fTrackPtMax = x;} 				   // Set selected track max pt cut
  void SetTrackEtaWindow(Float_t f)                    {fTrackEtaMin = -1.*f; fTrackEtaMax = f;} 	   // Set track eta window
  void SetTrackEtaWindow(Float_t fmin, Float_t fmax)   {fTrackEtaMin = fmin; fTrackEtaMax = fmax;} 	   // Set track eta window
  Float_t GetTrackEtaMin()                             {return fTrackEtaMin;}
  Float_t GetTrackEtaMax()                             {return fTrackEtaMax;}
  void SetExtTrackCutType(Int_t x)                     {fExtTrackCutType = x;} 			           // Temporary setter to select track qa cuts 
  void SetRequireITSRefit(Int_t i)                     {fkRequireITSRefit=i;}
  void SetSharedClusterCut(Int_t docut)                {fkApplySharedClusterCut=docut;}
  void SetTrackEfficiency(TObject* effi)               {fEffi=effi;}
  void SetTrackEfficiencyVarPercent(Double_t effivar) {fEffivar=effivar;}
  void SetTrackResolution(Double_t resol)              {fResol=resol;}
  void SetTrackResolutionVarPercent(Double_t resolvar) {fResolvar=resolvar;}
  void SetTrackFastSimParams(TObject* effi,Double_t resol, Int_t resolmeth = 0 ,Double_t effivar = 0, Double_t resolvar = 0) {fEffi=effi; fResol=resol; fResolMeth=resolmeth; fEffivar=effivar; fResolvar=resolvar;}

  // Jets
  void SetJetBranches(Int_t i, TString name)           {fJetBranch[i] = name;} 			           // Set the name of the jet branches  
  void SetBackgroundBranch(TString name)               {for(Int_t i=0; i<2; i++) fBkgJetBranch[i] = name;} // Set the name of the background jet branch -
  void SetJetOutputBranch(const char *c)               {fNonStdBranch = c;}        			   // Set the jet output branch name
  const char* GetJetOutputBranch()                     {return fNonStdBranch.Data();} 		           // Get the jet output branch name
  void SetJetOutputFile(const char *c)                 {fNonStdFile = c;} 				   // Set the name of the file used to store the output AliAODJet 
  const char* GetJetOutputFile()                       {return fNonStdFile.Data();} 			   // Get the name of the file used to store the output AliAODJet
  void SetAnaJetType( const char* c )                  {fAnaJetType = c;} 			           // Type of Jets Analysed
  void SetJetPtMin(Double_t f)                         {fJetPtMin = f;} 				   // Minimum pt of selected jets 
  void SetJetEtaWindow(Float_t f)                      {fJetEtaMin = -1.*f; fJetEtaMax = f;}   	           // Set jet eta window
  void SetJetEtaWindow(Float_t fmin, Float_t fmax)     {fJetEtaMin = fmin; fJetEtaMax = fmax;}             // Set jet eta window
  void SetJetDeltaPhiCut(Float_t dphi)                 {fJetDeltaPhiCut = dphi;}                           // Set delta phi limit for dijet selection. NB: pi+/-fJetDeltaPhiCut
  void SetJetDoMatching(Bool_t b, Bool_t hm = kFALSE)  {fkDoJetMatching=b; fkFillMismatchHisto=hm;}        // Enable jet matching and mismatched histograms
  void SetJetUseClosestJetMatching(Bool_t b)           {fkUseClosestJetsforMatching = b;}                  // Enable use of Closest jet for matching (fkDoJetMatching hast to be on)
  void SetJetMatchingFractionMin(Float_t x)            {fJetMatchingFractionMin = x;}                      // Set the minimum energy fraction for matching
  void SetJetMatchedDistMax(Double_t f)                {fJetMatchedDistMax = f;}                           // Set maximum distance between matched jets
  void SetJetMatchingParams(Float_t x, Double_t f, Bool_t hm = kFALSE)  
      {fkDoJetMatching=kTRUE; fkFillMismatchHisto=hm; fJetMatchingFractionMin = x;fJetMatchedDistMax = f;} // Set jet matching parameters
  void SetUseAODInputJets(Bool_t b)                    {fkUseJetFromInput = b;}                            // Set Read jets from Input
  void SetNUsedJets(Int_t f)                           {fNUsedJets = f;}                                   // Set number of used jet
  void SetUseTrackPtSumAsJetPt(Bool_t b)               {fkUseTrackPtSumAsJetPt =b;}                        // Use Track Pt Sum As Jet Pt
  void SetGenJetType(Bool_t f=0)                       {fkGenJetType = f;}                                 // Set Gen Jet type (0 = gen , 1 = rec)
  void SetEffJetType(Bool_t f=0)                       {fkEffJetType = f;}                                 // Set Eff Jet type (0 = gen , 1 = rec)
  void SetJetMinLTrackPt(Float_t pt = -1) { fJetMinLTrackPt = pt; }
  void SetJetMaxTrackPt(Float_t pt = -1) { fJetMaxTrackPt = pt; }
  void SetJetMinNTracks(Int_t nTracks = 0) { fJetMinnTracks = nTracks; }
  void SetTracksInJetMethod(Int_t method = 0) {fTracksInJetMethod = method; }
  void SetEffParams(Int_t method = 0, Bool_t genType = 0, Bool_t effType = 0) {fTracksInJetMethod = method; fkGenJetType = genType; fkEffJetType = effType;}
  void SetFFJtMode(Int_t mode =0, Int_t njT = 50, Double_t minjT= 4.-2*TMath::Pi(), Double_t maxjT = 4 )            {fFFJtValue = mode;  fnBinsAxis[6] = njT; fBinMinAxis[6] = minjT; fBinMaxAxis[6] = maxjT;}
  void SetFFBckgMode(Int_t mode = 0)         {fFFBckgMode = mode;}
  Float_t  GetJetMinLTrackPt() const { return fJetMinLTrackPt; }
  Float_t  GetJetMaxTrackPt() const { return fJetMaxTrackPt; }
  Float_t  GetJetMinNTracks() const { return fJetMinnTracks; }
  // For FFM analysis
  void SetFFMNParam(Int_t n, Float_t min, Float_t max) {fFFMMomMax = n; fFFMNMin = min; fFFMNMax = max;}   // Set FFM Nmin and Nmax
  void SetFFMScalePower(Double_t f)                       {fFFMScalePower = f;}                               // Set FFM scale power
  void SetFFMBckgTypeAndBounds(Int_t bType, Double_t bpar1=0, Double_t bpar2=0, double mu = 25)
  {fFFMBckgType = bType; fFFMBckgPar1 = bpar1; 
    fFFMBckgPar2 = bpar2; fFFMBckgMu = mu;}
  void SetHistosLevel(Int_t i)                         {fHistosLevel = i;}                               // Set output verbosity level 
  void SetHighResolution(Bool_t b)                     {fkHighResolution = b;}                           // Set parameter to set binning of histograms (*5)
  // Axis parameters
  void SetAxisForTracks(Int_t nPt = 40,  Double_t minPt = 0.,   Double_t maxPt = 200.,
			Int_t nEta = 100, Double_t minEta = -1., Double_t maxEta = 1.,
			Int_t nPhi = 90, Double_t minPhi = 0.,  Double_t maxPhi = 2*TMath::Pi())
  { fnBinsAxis[20] = nPt;  fBinMinAxis[20] = minPt;  fBinMaxAxis[20] = maxPt;
    fnBinsAxis[21] = nEta; fBinMinAxis[21] = minEta; fBinMaxAxis[21] = maxEta;
    fnBinsAxis[22] = nPhi; fBinMinAxis[22] = minPhi; fBinMaxAxis[22] = maxPhi;
  }
  void SetAxisForJets(Int_t nPt = 40,            Double_t minPt = 0.,             Double_t maxPt = 200.,
		      Int_t nEta = 100,           Double_t minEta = -1.,           Double_t maxEta = 1.,
		      Int_t nPhi = 90,           Double_t minPhi = 0.,            Double_t maxPhi = 2*TMath::Pi(),
		      Int_t nArea = 50,          Double_t minArea = 0.,           Double_t maxArea = 1.,
		      Int_t nNConstituents = 50, Double_t minNConstituents = 0.5, Double_t maxNConstituents = 50.5)
  { for (Int_t i = 0; i < 2; i++) {
      fnBinsAxis[10+i*5] = nPt;            fBinMinAxis[10+5*i] = minPt;            fBinMaxAxis[10+5*i] = maxPt;
      fnBinsAxis[11+i*5] = nEta;           fBinMinAxis[11+5*i] = minEta;           fBinMaxAxis[11+5*i] = maxEta;
      fnBinsAxis[12+i*5] = nPhi;           fBinMinAxis[12+5*i] = minPhi;           fBinMaxAxis[12+5*i] = maxPhi;
      fnBinsAxis[13+i*5] = nArea;          fBinMinAxis[13+5*i] = minArea;          fBinMaxAxis[13+5*i] = maxArea;
      fnBinsAxis[14+i*5] = nNConstituents; fBinMinAxis[14+5*i] = minNConstituents; fBinMaxAxis[14+5*i] = maxNConstituents;
    }
  }
  void SetFFMAxis(Int_t nFFM = 675,  Double_t minFFM = 0., Double_t maxFFM = 30.)
  {
    fnBinsAxis[8] = nFFM; fBinMinAxis[8] = minFFM; fBinMaxAxis[8] = maxFFM;
    fnBinsAxis[9] = nFFM; fBinMinAxis[9] = minFFM; fBinMaxAxis[9] = maxFFM;
  }
  void SetJetStructureAxis( Int_t nz = 44,   Double_t minz = 0.,                 Double_t maxz = 1.1,            //z
			    Int_t nxi = 70,  Double_t minxi = 0.,                Double_t maxxi = 7.,            //xi
                      //    Int_t njT = 50,  Double_t minjT = 4.-2*TMath::Pi(),  Double_t maxjT = 4.,            //lnjT
                            Int_t njT = 500,  Double_t minjT = 0.,                Double_t maxjT = 50.,          //1/jT
			    Int_t nDT = 50,  Double_t minDT = 0,                 Double_t maxDT = 0.5 )          //DeltaTheta
  { fnBinsAxis[4] = nz;  fBinMinAxis[4] = minz;  fBinMaxAxis[4] = maxz;
    fnBinsAxis[5] = nxi; fBinMinAxis[5] = minxi; fBinMaxAxis[5] = maxxi;
    fnBinsAxis[6] = njT; fBinMinAxis[6] = minjT; fBinMaxAxis[6] = maxjT;
    fnBinsAxis[7] = nDT; fBinMinAxis[7] = minDT; fBinMaxAxis[7] = maxDT;
  } 
  // Needed for PbPb
  void SetCentralityAxis( Int_t nv0 = 20,    Double_t minv0 = 0.,    Double_t maxv0 = 100.,          //v0 centrality
			  Int_t nntr = 50,   Double_t minntr = 0.,   Double_t maxntr = 50.,          //ntr
			  Int_t nep = 40,    Double_t minep = -20.,  Double_t maxep = 20.,           //event plane
			  Int_t nepb = 30,   Double_t minepb = 0.,   Double_t maxepb = TMath::Pi())  //event plane bin
  { fnBinsAxis[0] = nv0;  fBinMinAxis[0] = minv0;  fBinMaxAxis[0] = maxv0;
    fnBinsAxis[1] = nntr; fBinMinAxis[1] = minntr; fBinMaxAxis[1] = maxntr;
    fnBinsAxis[2] = nep;  fBinMinAxis[2] = minep;  fBinMaxAxis[2] = maxep;
    fnBinsAxis[3] = nepb; fBinMinAxis[3] = minepb; fBinMaxAxis[3] = maxepb;}
  // For jet reco
  void SetDoJetReco( Bool_t c = kFALSE)                {fkDoJetReco = c;}                                // Do jet reconstruction in this code (kTRUE) or read existing jet branch 
  void SetBackgroundCalc(Bool_t b)                     {fkUseBackgroundCalc = b;}                        // - 
  void SetRparam(Double_t f)                           {fRparam = f;}                                    // Set jet radius
  void SetAlgorithm(int f)                             {fAlgorithm = (fastjet::JetAlgorithm) f;}         // fj
  void SetAlgorithm(fastjet::JetAlgorithm f)           {fAlgorithm = f;}                                 // fj
  void SetBkgAlgorithm(int f)                          {fBkgAlgorithm = (fastjet::JetAlgorithm) f;}      // fj
  void SetBkgAlgorithm(fastjet::JetAlgorithm f)        {fBkgAlgorithm = f;}                              // fj
  void SetStrategy(int f)                              {fStrategy = (fastjet::Strategy) f;}              // fj
  void SetStrategy(fastjet::Strategy f)                {fStrategy = f;}                                  // fj
  void SetRecombScheme(int f) {                        fRecombScheme = (fastjet::RecombinationScheme) f;}// fj 
  void SetRecombScheme(fastjet::RecombinationScheme f) {fRecombScheme = f;}                              // fj
  void SetAreaType(int f)                              {fAreaType = (fastjet::AreaType) f;}              // fj
  void SetAreaType(fastjet::AreaType f)                {fAreaType = f;}                                  // fj 
  void SetGhostArea(Double_t f)                        {fGhostArea = f;}                                 // fj
  void SetActiveAreaRepeats(Int_t f)                   {fActiveAreaRepeats = f;}                         // fj
  void SetGhostEtaMax(Double_t f)                      {fGhostEtaMax = f;}                               // fj

  // for Fast Jet

  fastjet::JetAlgorithm        GetAlgorithm()    const {return fAlgorithm;}                              // fj
  fastjet::JetAlgorithm        GetBkgAlgorithm() const {return fBkgAlgorithm;}                           // fj
  fastjet::Strategy            GetStrategy()     const {return fStrategy;}                               // fj
  fastjet::RecombinationScheme GetRecombScheme() const {return fRecombScheme;}                           // fj
  fastjet::AreaType            GetAreaType()     const {return fAreaType;}                               // fj

  // we have different cases
  // AOD reading -> MC from AOD
  // ESD reading -> MC from Kinematics
  // this has to match with our selection of input events
  enum {kTrackUndef = 0, kTrackKineAll, kTrackKineCharged, kTrackKineAcceptance, kTrackKineChargedAcceptance,kTrackKineAcceptanceDet, kTrackKineChargedAcceptanceDet, kTrackAOD, kTrackAODQualityCuts, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance,  kTrackAODMCChargedSecS, kTrackAODMCChargedSecNS, kTrackAODextra, kTrackAODextraonly, kTrackAODMCextra, kTrackAODMCextraonly, kTrackAODMCHF};
  enum {kUndef = -1, kDoughnut = 0, kRectangle, kEtaRange, kPerp , kPerp2, kRapPhiRange};
  enum {kGenJet = 0, kRecJet};

 private:


  enum { fgkFFMNJetBranches = 2 };

  AliAnalysisTaskJetFFMoments(const AliAnalysisTaskJetFFMoments&);                                      // Not implemented 
  AliAnalysisTaskJetFFMoments& operator=(const AliAnalysisTaskJetFFMoments&);                           // Not implemented

  // Helper
  // 
  void       AssociateGenRec(TList* tracksAODMCCharged, TList* tracksRec, TArrayI& indexAODTr, TArrayI& indexMCTr, TArrayS& isRefGen, TH2D ** fh2RecVsGen);
  void       CalcTrackListFFM(AliAODJet * jet, TList * tracks, TArrayI * indexTr, TH2D** TrackProp, TProfile2D* fp2FFM);
  void       CalcFFAndFFMFromTracksInJet( AliAODJet* jet,  TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, TArrayI indexAODTr, const TArrayS& isRefGen, TList* jetTrackListTR, Bool_t scaleStrangeness,TH2D** hist, TProfile2D* prof);
  void       CalcSingleTrackEff(TH1D* trackEffPtGen, TH2D* trackEffEtaPhiGen, TH1D* trackEffPtRec, TH2D* trackEffEtaPhiRec, TList* tracksGen, const TArrayI& indexAODTr, const TArrayS& isRefGen, Bool_t scaleStrangeness=kFALSE);
 void FillTrackAssoList(AliAODJet* jet,TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, const TArrayI& indexAODTr, const TArrayS& isRefGen, TList* jetTrackListTR, Bool_t scaleStrangeness,  TArrayD* sWeight, TList* listRecTracks);
 void  FillFF(TVector3 TrackV, TVector3 JetV, Double_t sWeight, TH2D** histo);
  Float_t ComputeFF(Bool_t SelectTrackAsso, AliAODJet * jet, TArrayI indexTr, TList* tracks, TArrayD* sWeight, TH2D** histo);
  void ComputeFFM(Bool_t SelectTrackAsso, TList* tracks,  TArrayI indexTr, Float_t JetPt, TArrayD* sWeight, TProfile2D * p2FFM);
  void	     GetJetTracksTrackrefs(TList* l, const AliAODJet* j, Double_t minPtL, Double_t maxPt, Bool_t& isBadPt);
  void	     GetJetTracksPointing(TList* in, TList* out, const AliAODJet* j, Double_t r, Double_t& sumPt,  Double_t minPtL, Double_t maxPt, Bool_t& isBadPt);
  void       GetTracksInJet(TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t& sumPt, Bool_t& isBadPt);
  void       GetTracksTiltedwrpJetAxis(Float_t alpha, TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius,Double_t& sumPt);
  Int_t      GetListOfTracks(TList *list, Int_t type); 
  Int_t      GetListOfMatchedJets(TList** listJets, Int_t* nUsedJets, TList* listMatchedJets1, TList* listMatchedJets2, Int_t ifirstBr);
  Int_t      GetListOfClosestJets(TList** listJets, Int_t* nUsedJets, TList* listMatchedJets1, TList* listMatchedJets2, Int_t ifirstBr);
  Bool_t     TrackQAFilter(TList* list, Int_t type);
  Int_t      GetListOfJets(TList *list, Int_t js);     
  void       CreateHistos();                           
  Bool_t     ClassifyJetEvent(Int_t njets, Int_t ijet); // To be checked
  Bool_t     SelectJet(AliAODJet * jet, TList* tracksAfterCut); // Jet selection: eta, pt, track pt, etc. cuts
  TList*     GetTrackInJetAssoList(AliAODJet* jet,  TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, TArrayI indexAODTr, const TArrayS& isRefGen, TList* jetTrackListTR, Bool_t scaleStrangeness);
  void       PseudoJetsToAODJets( const vector<fastjet::PseudoJet> & fPseudojet, fastjet::ClusterSequenceArea & CluSeq, TList & tracks, TList *jets); // Convert a pseudo jet in AOD jet
  int        AliAODJetToPseudoJet(AliAODJet * jet, fastjet::PseudoJet & fPseudojet); // Convert a AOD jet in a pseudo jet
  int        AliAODJetToPseudoJet(AliAODJet * jet, TList* list, fastjet::PseudoJet & fCurrentPseudojet);
  Bool_t     PropertiesInJet(TVector3 & jetV, TVector3 & trackV, Double_t & Z, Double_t & xi, Double_t & lnjT, Double_t & deltaTheta, Double_t & pt, Double_t & eta, Double_t & phi);
  THnSparse* CreateTHnSparseF(const char* name, UInt_t entries, UInt_t res);
  TH2D*      CreateTH2D(const char* name, Int_t iXAxis, Int_t iYAxis, Bool_t res);
  TProfile2D*      CreateTProfile2D(const char* name, Int_t iXAxis, Int_t iYAxis, Int_t iZAxis, Bool_t res);
  TH1D*      CreateTH1D(const char* name, Int_t iXAxis, Bool_t res);
  void       GetDimParams(Int_t iEntry, Bool_t hr, const char* &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  Int_t      AddDaughters(TList * list, AliAODMCParticle *part, TClonesArray * tca);
  Bool_t     AvoidDoubleCountingHF(AliAODEvent *aod, AliAODTrack *tr1);
  bool       IsBMeson(int pc);
  bool       IsDMeson(int pc);
  AliGenPythiaEventHeader* GetPythiaHeader();
  AliGenHerwigEventHeader* GetHerwigHeader();
  Bool_t IsOutlier(AliGenPythiaEventHeader * const header);
  bool       IsProof() {return AliAnalysisManager::GetAnalysisManager()->IsProofMode();}  
  const char* ProofClearOpt() { if(IsProof()) {return "nodelete";} else { return "";}} 
  fastjet::PseudoJet           join_with_area(const vector <fastjet::PseudoJet> & pieces,
                                              const fastjet::PseudoJet & area_4vector,
                                              const double area,
                                              const double area_error = 0,
                                              const bool   is_pure_ghost = false);                       // fj
 Double_t GetMCStrangenessFactorCMS(AliAODMCParticle* daughter);

  AliAODEvent*       fAOD;			    //! Input AOD
  AliAODEvent*       fAODJets;			    // Input jet AOD
  AliAODExtension*   fAODExtension;		    //! AOD extension in case we write a non-std branch to a separate file
  Float_t            fFilterPt;                     //  Use this as a switch for writing the AOD, minium p_T of leading jet
  TRefArray*         fRef;			    //! TRefArray for track references within the jet
  Bool_t             fkIsPbPb;                      // PbPb case
  Bool_t 	     fkEventSelection;		    // Use to bypass event selection?
  Bool_t             fkRejectFastOnly;              // Switch to reject FastCluster
  Bool_t 	     fkRequireVZEROAC;		    // Switch to require V0 AC
  Bool_t 	     fkRequireTZEROvtx;		    // Switch to require T0 vtx
  Bool_t             fkRejectPileup;                // Switch to reject pileup
  Float_t 	     fCentCutUp;		    // Upper limit on centrality
  Float_t 	     fCentCutLo;		    // Lower limit on centrality
  Int_t              fCentClass[4];                 // Centrality classes
  Float_t 	     fVtxZMax;			    // Maximum vtx Z
  Float_t 	     fVtxR2Max;		            // Maximum vtx R^2 

  // Tracks
  UInt_t 	     fFilterMask;		    // Filter mask for track selection
  UInt_t 	     fFilterMaskBestPt;		    // Filter bit to mark jets with high quality leading tracks
  UInt_t 	     fFilterType;		    // Filter type: 0 = all, 1 = ITSTPC, 2 = TPC
  Bool_t 	     fkUseAODTrackInput;	    // Use input (kTrue) AOD for tracks
  Bool_t 	     fkUseAODMCInput;		    // Use input (kTrue) AOD for MC tracks
  Int_t 	     fTrackType[fgkFFMNJetBranches]; // Track type, 0: generated, 1: reconstructed 
  Float_t 	     fTrackPtMin;		    // Minimum track pt
  Float_t 	     fTrackPtMax;	            // Maximum track pt (Get rid of bad tracks)
  Int_t 	     fExtTrackCutType;		    // Temporary setter to select track qa cuts: 1: pt_Chi2, 2: eta_Chi2, etc. 
  Bool_t             fkUseHFcuts;                   // HF cuts for track selection
  Int_t              fkRequireITSRefit;             // to select hybrids with ITS refit only
  Int_t              fkApplySharedClusterCut;       // flag to apply shared cluster cut (needed for some AODs where this cut was not applied in the filtering)
  Float_t 	     fTrackEtaMin;		    // Eta track min
  Float_t 	     fTrackEtaMax;		    // Eta track max
  Float_t          fTrackPhiMin;   // track phi cut
  Float_t          fTrackPhiMax;   // track phi cut

  // Jets
  TClonesArray*      fTCAJetsOut;	                  // TCA of output jets
  AliAODJetEventBackground* fAODJetBackgroundOut;         //! Jet background to be written out
  TList*             fListJets[fgkFFMNJetBranches];        // Jet list, 0: generated, 1: reconstructed
  TList*             fListMatchedJets[fgkFFMNJetBranches]; // Matched jet list, 0: generated, 1: reconstructed
  TList* fTracksAODMCCharged; //! AOD MC tracks
  TList* fTracksAODMCChargedSecNS; //! AOD MC tracks - secondaries (non-strangeness)
  TList* fTracksAODMCChargedSecS;  //! AOD MC tracks - secondaries (from strangeness)
  TList* fTracksRecQualityCuts;    //! reconstructed tracks after quality cuts, no acceptance/pt cut
  TString   	     fJetBranch[fgkFFMNJetBranches];       // Jet branch name, 0: generated, 1: reconstructed
  TString            fBkgJetBranch[fgkFFMNJetBranches];    // Bkg jet branch name, 0: generated, 1: reconstructed
  // output configuration
  TString 	     fNonStdBranch;		    // Non-std branch name
  TString 	     fNonStdFile;		    // The optional name of the output file the non-std branch is written to
  TString            fAnaJetType;	            // Leadingjet, Dijet, AllJet
  Double_t 	     fJetPtMin;			    // Min Pt of jets
  Float_t 	     fJetEtaMin;		    // Eta jet min
  Float_t 	     fJetEtaMax;		    // Eta jet max
  Float_t            fJetDeltaPhiCut;               // Delta phi limit for dijet selection (pi+/-fJetDeltaPhiCut)
  Bool_t             fkDoJetMatching;               // Switch on jet matching (kTRUE)
  Bool_t             fkUseClosestJetsforMatching;   // Switch on CLosest jet matching (kTRUE);
  Bool_t             fkFillMismatchHisto;           // Fill matched histos (kFALSE)
  Float_t 	     fJetMatchingFractionMin;       // Minimum energy fraction for matching
  Double_t           fJetMatchedDistMax;     	    // Maximum distance in eta-phi space for matching
  Bool_t             fkUseJetFromInput;             // Read jets from input
  Int_t              fNUsedJets;                    // Number of jets used
  Float_t        fJetMinLTrackPt;   // Reject jets with leading track with pt smaller than this value
  Float_t        fJetMaxTrackPt;    // Reject jets containing any track with pt larger than this value
  Int_t            fJetMinnTracks;    // Reject jets with less tracks than this value
  Bool_t           fkEffJetType;       // Eff jet type
  Bool_t           fkGenJetType;       // Gen Jet type
  Int_t            fTracksInJetMethod; // Method of track selection in jet
  Int_t           fFFBckgMode;       // FF background mode (kFALSE=off)
  Int_t           fFFJtValue;
  // For FFM analysis 
  Float_t  	     fFFMNMin;			    // Minimum N of M_N
  Float_t  	     fFFMNMax;			    // Maximum N of M_N
  Int_t	             fFFMMomMax;		    // Number of bins in M_N
  Double_t           fFFMScalePower;		    // ((N+1)/2)^(fFFMScalePower)*M_{N}
  Int_t              fFFMBckgType;                  // FFM background type: 0 = doughnut, 1 = rectangle, 2 = eta range
  Double_t           fFFMBckgPar1;                  // First parameter for FFM bckg bound definition
  Double_t           fFFMBckgPar2;                  // Second parameter for FFM bckg bound definition
  double             fFFMBckgMu;                    // Mu parameter for improved background subtraction in FFM
  Int_t 	     fHistosLevel;		    // Output verbosity level
  Bool_t             fkHighResolution;        	    // Parameter to change binning resolution of histograms (*5)
  Bool_t            fkUseTrackPtSumAsJetPt;            // Use MC jet to get the associated track FFM, otherwise, use sum pt
  Int_t	             fnBinsAxis[32];		    // Possible nBins for axis used in this class for all plots
  Double_t	     fBinMinAxis[32];		    // Possible Bin Min for axis used in this class for all plots
  Double_t	     fBinMaxAxis[32];		    // Possible Bin Max for axis used in this class for all plots
  //  0 - 9
  //  0 |  1  |  2  |  3   |  4  |  5 |   6   |      7      |    8    |    9    |
  // vz | ntr | ep  | epb  |  z  | xi | lnjT  |  DeltaTheta | FFM_gen | FFM_rec |
  // 10 - 24
  //  0 |  1  |  2  |    3     |      4        |||      n
  // pt | eta | phi |   area   | nconstituents |||     gen  2x5+n
  // pt | eta | phi |   area   | nconstituents |||     rec  3x5+n
  // pt | eta | phi |          |    N(order)   |||    track 4x5+n
  // 25 - 32
 //     25    |    26    |
  //  fraction | RDeltaPt | R for Relative
  // For jet reco (fastjet)
  Bool_t                        fkDoJetReco;                   // Do jet reco in this code (kTRUE) or read jet input branch (kFALSE)
  Bool_t 	                  fkUseBackgroundCalc;	       // jet background subtraction 
  Double_t 			  fRparam;		       // fastjet distance parameter
  fastjet::JetAlgorithm   	  fAlgorithm;		       // fastjet::antikt_algorithm
  fastjet::JetAlgorithm   	  fBkgAlgorithm;	       // fastjet::kt_algorithm
  fastjet::Strategy 		  fStrategy;		       // fastjet::Best;
  fastjet::RecombinationScheme  fRecombScheme;	               // fastjet::BIpt_scheme;
  fastjet::AreaType 		  fAreaType;		       // fastjet area type
  Double_t 			  fGhostArea;		       // fasjet ghost area
  Int_t 			  fActiveAreaRepeats;	       // fast jet active area repeats
  Double_t 			  fGhostEtaMax;		       // fast jet ghost area

  // Output lists
  TList* 	    fHistListJets[fgkFFMNJetBranches];	       // To be checked (how they are deleted) 
  TList* 	    fHistList;				       //! Final list of histos

  // Output histos (verbosity 2)
  // Parameters, pythia, event
  TProfile*       fPrParameters;		               //! Summary of the analysis parameters
  TProfile*       fh1Xsec;			               //! Pythia cross section profile
  TH1F*           fh1Trials;			               //! Pythia number of trials
  TH1F*           fh1AvgTrials;		                       //! Average number of trials from pythia
  TH1F*           fh1CentralityPhySel;	                       //! Centrality after physics selection
  TH1F*           fh1CentralitySelect;	                       //! Centrality with physics selection bypassed
  TH1F*           fh1vZPhySel;		                       //! Vertex z after physics selection
  TH1F*           fh1vZSelect;		                       //! Vertex z with physics selection bypassed
  // Tracks qa
  TH1D*	    fh1GenTracks[3];		                       //! Histos in 0: pt, 1: eta and 2: phi of generated tracks
  TH1D*	    fh1RecTracks[3];		                       //! Histos in 0: pt, 1: eta and 2: phi of reconstructed tracks
  TH2D*         fh2PtRecVsGenPrim[3];                             //! Histos in 0: pt, 1: eta and 2: phi of reconstructed tracks associated to Physics Primary
  TH2D*         fh2PtRecVsGenSec[3];                             //! Histos in 0: pt, 1: eta and 2: phi of reconstructed tracks associated to Physics Primary
  TH1D*         fh1TrackEffPtGen;                                //! Single track efficiency (pT Gen)
  TH2D*         fh2TrackEffEtaPhiGen;                            //! Single track efficiency (eta,phi Gen)
  TH1D*         fh1TrackEffPtRec;                                //! Single track efficiency (pT Rec)
  TH2D*         fh2TrackEffEtaPhiRec;                            //! Single track efficiency (eta,phi Rec)
  TObject*          fEffi;                                           // Efficiency parametrisation (fast mc)
  Double_t      fResol;                                          // Track resolution
  Int_t         fResolMeth;                                      // Track resolution smearing method
  Double_t      fEffivar;                                        // Efficiency shift (%)
  Double_t      fResolvar;                                       // Track resolution syst (%)
  TH2D*         fh2TrackResPt;                                   //! Track resolution Histo
  TH1D*         fh1TrackResPtInv;                                //! Track resolution Histo (1/pT)
  Double_t      fPtHardAndPythiaJetPtFactor;                     //  Outliers cut : PtHard vs Pythia jets
  Double_t      fPtHardAndTrackPtFactor;                         //  Outliers cut : PtHard vs Track pt
  TH2D*         fh2PtHardVsPt[2];                                //! Outliers Histo (jet,track)
  TH2D*         fh2PtHardVsPtCut[2];                             //! Cut Outliers Histo (jet,track)
  // Jet analysis plots: Aj, jet structure and matching
  TH1F* 	    fh1Njets[fgkFFMNJetBranches];	       //! Number of jets, 0: generated, 1: reconstructed
  TH1F* 	    fh1Asy_DiJets[fgkFFMNJetBranches];  	       //! Aj, 0: generated, 1: reconstructed
  TH2D*         fh2TracksInJets[fgkFFMNJetBranches][7];                  //! Distributions of (0: generated, 1: reconstructed)
                                                                        //  tracks in jets: 0: Z, 1: log(1/Z), 2: ln(jT), 3: Delta Theta,4: pt, 5: eta, 6: phi
  TH2D*         fh2AssociatedTracksInJets[fgkFFMNJetBranches][7]; //! Distributions of (0: generated, 1: reconstructed)
  TH2D*         fh2AssociatedTracksInJetsSecNS[fgkFFMNJetBranches][7]; //! Distributions of (0: generated, 1: reconstructed)
  TH2D*         fh2AssociatedTracksInJetsSecS[fgkFFMNJetBranches][7]; //! Distributions of (0: generated, 1: reconstructed)
  TH2D*         fh2AssociatedTracksInJetsSecSsc[fgkFFMNJetBranches][7]; //! Distributions of (0: generated, 1: reconstructed)
                                                                        //  associated tracks in matched jets: 0: Z, 1: log(1/Z), 2: ln(jT), 3: Delta Theta,4: pt, 5: eta, 6: phi
  TH1D*             fh1JetPr_Mismatched[fgkFFMNJetBranches][5]; //! Mismatched (0: generated, 1: reconstructed) jet property distributions: 0: pt, 1: eta, 2: phi, 3: area, 4: constituents
  TH2D*         fh2MismatchedJetsAreaVSPt[fgkFFMNJetBranches];           //! Histos of area vs Pt
  TH2D*         fh2MatchedJetsRDPtVSPt[fgkFFMNJetBranches];              //! Histos of DeltaPt vs Pt
  TH2D*         fh2MatchedJetsUE[fgkFFMNJetBranches];            //! Histos of UEptSum vs Pt
  TH2D*         fh2MatchedJetsAreaVSPt[fgkFFMNJetBranches];              //! Histos of area vs Pt
  TH2D*	            fh2MatchedJets[5];		               //! Histos in 0: pt, 1: eta, 2: phi, 3: area, 4: constituents of jets after matching

  // FFM analysis output (verbosity 2)
  THnSparse*	    fhnJetFFM_Raw[fgkFFMNJetBranches];	       //! Raw FFM, 0: generated, 1: reconstructed
  THnSparse*	    fhnJetFFM_Sub[fgkFFMNJetBranches];	       //! Background subtracted FFM, 0: generated, 1: reconstructed
  THnSparse*	    fhnJetFFM_Imp[fgkFFMNJetBranches];	       //! Improved background subtracted FFM, 0: generated, 1: reconstructed

  TProfile2D*  fp2TracksInJetFFM[fgkFFMNJetBranches];             //! FFM of tracks in Jets (computed outside fj), 0: generated, 1: reconstructed
  TProfile2D*   fp2AssociatedTracksJetFFM[fgkFFMNJetBranches];       //! associated tracks FFM, 0: generated, 1: reconstructed
  TProfile2D*   fp2AssociatedTracksJetFFMSecNS[fgkFFMNJetBranches];       //! associated tracks FFM, 0: generated, 1: reconstructed
  TProfile2D*   fp2AssociatedTracksJetFFMSecS[fgkFFMNJetBranches];       //! associated tracks FFM, 0: generated, 1: reconstructed
  TProfile2D*   fp2AssociatedTracksJetFFMSecSsc[fgkFFMNJetBranches];       //! associated tracks FFM, 0: generated, 1: reconstructed
  TProfile2D*   fp2JetFFM_Raw[fgkFFMNJetBranches];                   //! Raw FFM, 0: generated, 1: reconstructed
  TProfile2D*   fp2JetFFM_Sub[fgkFFMNJetBranches];                   //! Background subtracted FFM, 0: generated, 1: reconstructed
  TProfile2D*   fp2JetFFM_Imp[fgkFFMNJetBranches];                   //! Improved background subtracted FFM, 0: generated, 1: reconstructed
  // (verbosity 9)
  THnSparse**     fhnJetMomN_Raw;			       //! Raw M_N distributions (fFFMMomMax distributions)
  THnSparse**     fhnJetMomN_Sub;			       //! M_N distributions after background subtraction (Cacciari's et al. method)
  THnSparse**     fhnJetMomN_Imp;			       //! M_N distributions after improved background subtraction (Cacciari's et al. method)

  //----------------------------------------------------------------------
  /// a class that behaves like CompositeJetStructure, except that the user
  /// must fully specify the area properties.

  class AliCompositeJetStructureUserArea : public fastjet::CompositeJetStructure 
    {

    public:
    AliCompositeJetStructureUserArea(const vector<fastjet::PseudoJet> & particles,
				     const fastjet::PseudoJet & area_4vector,
				     const double area,
				     const double area_error = 0,
				     const bool   is_pure_ghost = false,
				     const fastjet::JetDefinition::Recombiner * recombiner = 0
				     ) : fastjet::CompositeJetStructure(particles, recombiner),
	_area_4vector(area_4vector),
	_area(area),
	_area_error(area_error),
	_is_pure_ghost(is_pure_ghost)
	{
	}

      /// check if it has a well-defined area
      bool has_area() const {return true;}

      /// return the jet (scalar) area.
      double area(const fastjet::PseudoJet &reference) const {return _area;}

      /// return the error (uncertainty) associated with the determination
      /// of the area of this jet.
      ///
      /// Be conservative: return the sum of the errors
      double area_error(const fastjet::PseudoJet &reference) const {return _area_error;}

      /// return the jet 4-vector area.
      fastjet::PseudoJet area_4vector(const fastjet::PseudoJet &reference) const {return _area_4vector;}

      /// true if this jet is made exclusively of ghosts.
      ///
      /// In this case, it will be true if all pieces are pure ghost
      bool is_pure_ghost(const fastjet::PseudoJet &reference) const {return _is_pure_ghost;}


    private:
      fastjet::PseudoJet _area_4vector;
      double _area, _area_error;
      bool   _is_pure_ghost;

    };

     /// helper for selecting on objects within a distance 'radius' of a eta-tiled reference (by an angle alpha)
     class AliFJetSW_Perp : public  fastjet::SelectorWorker {
     public:
      AliFJetSW_Perp(const double &radius, const double &alpha) : fastjet::SelectorWorker(),_is_initialised(false), _radius2(radius*radius), _alpha(alpha) {}
      /// returns true if the worker takes a reference jet
      virtual bool takes_reference() const { return true;}

     /// sets the reference jet
    virtual void set_reference(const fastjet::PseudoJet &centre){
     _is_initialised = true;
     _reference = centre;
    // cout << "ref eta: " << _reference.eta() << ", phi: " << _reference.phi() << ", pt: " << _reference.pt() << ", area: " << _reference.area() << endl;
   }

    /// return a copy of the current object
    virtual fastjet::SelectorWorker* copy(){ return new AliFJetSW_Perp(*this);}

    /// returns true if a given object passes the selection criterium
    /// this has to be overloaded by derived workers
    virtual bool pass(const fastjet::PseudoJet & jet) const {
    // make sure the centre is initialised
    if (! _is_initialised)
    throw fastjet::Error("To use a SelectorPerp (or any selector that requires a reference), you first have to call set_reference(...)");

    // tilt reference jet phi by angle _alpha
    Double_t etaTilted = _reference.eta();
    Double_t phiTilted = _reference.phi_02pi() + _alpha;
    if(phiTilted > 2*TMath::Pi()) phiTilted = phiTilted - 2*TMath::Pi();

    double dphi = std::abs(jet.phi() - phiTilted);
    if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
    double drap = jet.eta() - etaTilted;

    return (dphi*dphi + drap*drap) <= _radius2;

   }

    /// returns a description of the worker
    virtual std::string description() const {
     std::ostringstream ostr;
     ostr << "distance from the centre <= " << sqrt(_radius2);
     ostr << " angle: " << _alpha;
     return ostr.str();
     }

   /// returns the rapidity range for which it may return "true"
   virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
     // make sure the centre is initialised
     if (! _is_initialised)
       throw fastjet::Error("To use a SelectorPerp (or any selector that requires a reference), you first have to call set_reference(...)");

     rapmax = _reference.rap()+sqrt(_radius2);
     rapmin = _reference.rap()-sqrt(_radius2);
   }

    virtual bool is_geometric() const { return true;}    ///< implies a finite area
    virtual bool has_finite_area() const { return true;} ///< regardless of the reference
    virtual bool has_known_area() const { return true;}  ///< the area is analytically known
    virtual double known_area() const {
    return TMath::Pi() * _radius2;
   }

      protected:
      fastjet::PseudoJet _reference;
      bool _is_initialised;
      double _radius2;
      double _alpha;
 };

   // select on objets within a distance 'radius' of a variable location
   fastjet::Selector SelectorPerp(const double & radius, const double &alpha) {
   return fastjet::Selector(new AliFJetSW_Perp(radius,alpha));
  }


  ClassDef(AliAnalysisTaskJetFFMoments, 2);

};

#endif
