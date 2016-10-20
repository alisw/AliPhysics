/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
//===========================================================

#ifndef ALIJJETJTANALYSIS_H
#define ALIJJETJTANALYSIS_H

#include <vector>
#include <TObjArray.h>
#include <TVector.h>
#include "AliJCard.h"
#include "AliJJet.h"
#include "AliJHistManager.h"
#include "AliJJetAnalysis.h"
#include <TRandom3.h>  //FK//

class AliJEfficiency;
class TString;

class AliJJetJtAnalysis{
  public:

    enum { kJUndefined = -1 , kJRecoTrack, kJMCParticle };

    AliJJetJtAnalysis();
    AliJJetJtAnalysis( AliJCard * card );
    AliJJetJtAnalysis(const AliJJetJtAnalysis& ap);
    AliJJetJtAnalysis& operator = (const AliJJetJtAnalysis& ap);
    ~AliJJetJtAnalysis();

    //type = charged jet or full jet, js = cone size 

    void FillHistosJets();

    void AddJets(TObjArray * jets, int type ){ 
      if( !jets ) {
        //return;
      }
      fJetListOfList.Add( (TObject*)jets ); 
      fTrackOrMCParticle.push_back(type);
      //if( !jets ) return;
      for( int i=0;i<jets->GetEntriesFast();i++ ){
        //((AliJJet*)jets->At(i))->ReSum();
      }
    } // TODO clean before event

    void AddMCJets(TObjArray * jets ){ 
      if( !jets ) {
        //return;
      }
      fMCJetListOfList.Add( (TObject*)jets ); 
      //if( !jets ) return;
      for( int i=0;i<jets->GetEntriesFast();i++ ){
        //((AliJJet*)jets->At(i))->ReSum();
      }
    } // TODO clean before event
    void SetJTracks(TClonesArray *tracks ){fTracks = tracks ;}
    void SetMCJTracks(TClonesArray *tracks ){fMCTracks = tracks ;}


    int GetNJets(){ return GetJetList()->GetEntriesFast(); }
    TObjArray* GetJetList(){ return fJetList; }
    //Double_t GetJetEtaRange(){ return fJetEtaRange; }
    //void SetJetEtaRange(double eta){ fJetEtaRange=eta; }
    void SetJetList(TObjArray* jetlist){ fJetList=jetlist; }
    void SetMCJetList(TObjArray* jetlist){ fMCJetList=jetlist; }
    void SetInputList(TObjArray * ilist){ fInputList = ilist;}
    void SetTrackOrMCParticle( UInt_t i, int v ){ fTrackOrMCParticle[i] = v; }
    int  GetTrackOrMCParticle( UInt_t i ){ return fTrackOrMCParticle.at( i ); }
    //void SetTrackJetMap(std::vector<int> * v){ fTrackJetMap=v;}
    //void SetJetPtBin( TVector * b){ fJetPtBins=b; }
    void SetCard(AliJCard * card) {fCard = card;}

    AliJJet & Jet( int i ){ return *(AliJJet*) fJetList->At(i); }
    AliJJet & MCJet( int i ){ return *(AliJJet*) fMCJetList->At(i); }

    void UserCreateOutputObjects();
    void UserExec();

    void ClearBeforeEvent();
    void CreateMCHistograms();
    void FillJtHistogram( TObjArray *Jets, int iContainer , int mc);
    void FillRandomBackground(double jetpT, double jetE, TObjArray *Jets , int iContainer, int mc);

    void FillCorrelation( TObjArray *Jets, TObjArray *mcJets, int iContainer, int iContainerParticle);
    void FillPythia(TObjArray *Jets, int iContainer); 
    int FindPythiaJet(int iContainer, int itrack);

    void FillBgJtWithSmallerR(const TClonesArray &Jets, 
        double nR, int iHist);

    void FillBgJtWithDiffAxes (int iao, int ia, int iHist);

    void WriteHistograms();

    double getDiffR(double phi1, double phi2, double eta1, double eta2);
    int GetBin(TVector *array, double val){

      int iBin=-1;
      for(int i=1; i< array->GetNoElements(); i++){
        if((*array)[i] <= val && val<(*array)[i+1]){
          iBin=i-1;
          break;
        }
      }

      return iBin;
    }
    void SetJetFinderName(vector<TString> JetFinderName){ fJetFinderName = JetFinderName; }
    void SetMCJetFinderName(vector<TString> JetFinderName){ fMCJetFinderName = JetFinderName; }

    // Need for event loop
    void SetCentralityBin( int cbin) { cBin = cbin;}
    void SetCentrality( float cent) { fcent = cent;}
    void SetZVertex( float zvtx) { zVert = zvtx;}
    void SetZVertexBin( int zbin) { zBin = zbin;}
    void SetNrandom(int nRndm) { Nrandom = nRndm;}
    void SetMoveJet(int move) { moveJet = move;}
    void SetMC(int mc) {fDoMC = mc;};
    void SetNumberOfJetFinders( int njfinder ) { nJetContainer = njfinder;}
    AliJEfficiency* GetAliJEfficiency() { return fEfficiency;}




  private:
    TObjArray * fInputList; // comment needed
    TObjArray * fJetList; // comment needed
    TObjArray * fMCJetList; // comment needed
    TObjArray fJetListOfList; // !comment needed
    TObjArray fMCJetListOfList; // !comment needed
    vector<int>               fTrackOrMCParticle;
    vector<TClonesArray>      fJetBgListOfList;
    double fJetEtaCut;
    TRandom3 *frandom; // comment me

    TVector  *fJetTriggPtBorders;
    TVector  *fJetConstPtLowLimits;
    TVector  *fJetAssocPtBorders;
    TVector  *fDeltaRBorders;
    TVector *fJetLeadPtBorders;
    TVector *fJetMultBorders;
    int nJetContainer;

    AliJCard * fCard; // comment needed
    AliJJetAnalysis *fJJetAnalysis;
    vector<TString> fJetFinderName;
    vector<TString> fMCJetFinderName;
    vector<double>  fConeSizes;
    // Need for events
    AliJEfficiency *fEfficiency;
    int cBin;
    float fcent;
    int zBin;
    float zVert;
    TClonesArray *fTracks;
    TClonesArray *fMCTracks;
    int iConst;

    TVector *fTrackJt;
    TVector *fConstJt;
    TVector *fTrackPt;
    TVector *fConstPt;
    TVector *fConstLabels;
    TVector *fJetPt;
    TVector *fJetPt2;
    TVector *fTrackFound;
    TVector *fConstFound;
    TVector *fBin2;
    TVector *fBin3;
    TVector *fpta;
    TVector *fptt;
    int Nrandom;
    int moveJet;
    int fDoMC;
    //Histograms
    AliJHistManager * fHMG;
    AliJHistManager * fHMGMC;

    AliJBin fJetFinderBin; 
    AliJBin fJetTriggerBin; 
    AliJBin fTrkPtBin; 
    AliJBin fTrkLimPtBin; 
    AliJBin fJetLeadPtBin;
    AliJBin fJetMultBin;
    AliJBin fdRBin;
    AliJBin fiHist;
    AliJBin fJetFinderBinMC; 
    AliJBin fJetTriggerBinMC; 
    AliJBin fTrkPtBinMC; 
    AliJBin fTrkLimPtBinMC; 
    AliJBin fJetLeadPtBinMC;
    AliJBin fJetMultBinMC;
    AliJBin fdRBinMC;
    AliJBin fiHistMC;
    AliJBin fiHist2MC;
    AliJTH1D fhNumber;
    AliJTH1D fhKNumber;
    AliJTH1D fhJetPt ;
    AliJTH1D fhJetPtBin;
    AliJTH1D fhJetPtMultiplicityCutBin;
    AliJTH1D fhJetPtTrackCutBin;
    AliJTH1D fhJetPtWeight ;
    AliJTH1D fhJetPtWeightBin;
    AliJTH1D fhJetMultiplicityBin;
    AliJTH1D fhZ ;
    AliJTH1D fhZBin;
    AliJTH1D fhJt ;
    AliJTH1D fhJtBin;
    AliJTH1D fhJtWeightBin;
    AliJTH1D fhJtWeightBinTest;
    AliJTH1D fhJtWeightBinTest2;
    AliJTH1D fhLogJtWeightBin;
    AliJTH1D fhLogJtWeight2Bin;
    AliJTH1D fhJtWithPtCutWeightBinBin;
    AliJTH1D fhLogJtWithPtCutWeightBinBin;
    AliJTH1D fhLogJtWithPtCutWeight2BinBin;
    AliJTH1D fhJtBinLimBin;
    AliJTH1D fhJtWeightBinLimBin;
    AliJTH1D fhLogJtWeightBinLimBin;
    AliJTH1D fhLogJtWeight2BinLimBin;

    //Histograms for jt in cone
    AliJTH1D fhJetConeTrkPt; //
    AliJTH1D fhJetConeTrkPtBin; //
    AliJTH1D fhJetConeTrkPtWeightBin; //
    AliJTH1D fhJetConeZ;
    AliJTH1D fhJetConeZBin;
    AliJTH1D fhJetConeJt;
    AliJTH1D fhJetConeJtBin;
    AliJTH1D fhJetConeJtWeightBin;
    AliJTH1D fhJetConeJtWeightBinTest;
    AliJTH1D fhJetConeJtWeightBinTest2;
    AliJTH1D fhJetConeJtWeightWithTrackCutBinBin;
    AliJTH1D fhJetConeJtWeightWithMultiplicityCutBinBin;
    AliJTH1D fhJetConeLogJtWeightBin;
    AliJTH1D fhJetConeLogJtWeight2Bin;
    AliJTH1D fhJetConeJtWithPtCutWeightBinBin;
    AliJTH1D fhJetConeLogJtWithPtCutWeightBinBin;
    AliJTH1D fhJetConeLogJtWithPtCutWeight2BinBin;

    //Unfolding Background histograms
    AliJTH1D fhJetConeJtUnfBg;
    AliJTH1D fhJetConeJtBinUnfBg;
    AliJTH1D fhJetConeJtWeightBinUnfBg;
    AliJTH1D fhJetConeLogJtWeightBinUnfBg;
    AliJTH1D fhJetConeLogJtWeight2BinUnfBg;

    AliJTH1D fhJetBgPt ;
    AliJTH1D fhJetBgPtBin;
    AliJTH1D fhBgZ ;
    AliJTH1D fhBgZBin;
    AliJTH1D fhBgJt ;
    AliJTH1D fhBgJtBin;
    AliJTH1D fhBgJtWeightBin;
    AliJTH1D fhBgLogJtWeightBin;
    AliJTH1D fhBgLogJtWeight2Bin;
    AliJTH1D fhBgJtWithPtCutWeightBinBin;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBin;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBin;
    AliJTH1D fhBgJtWithPtCutWeightBinBinSmallerR;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinSmallerR;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinSmallerR;
    AliJTH1D fhBgJtWithPtCutWeightBinBinDiffR;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinDiffR;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinDiffR;
    AliJTH1D fhBgJtBinLimBin;
    AliJTH1D fhBgJtWeightBinLimBin;
    AliJTH1D fhBgLogJtWeightBinLimBin;
    AliJTH1D fhBgLogJtWeight2BinLimBin;
    AliJTH1D fhTrkPt;
    AliJTH1D fhTrkPtBin;
    AliJTH1D fhTrkPtWeightBin;
    AliJTH1D fhLeadingTrkPtBin;
    AliJTH1D fhBgTrkPt;
    AliJTH1D fhBgTrkPtBin;
    AliJTH1D fhBgTrkPtWeightBin;
    AliJTH1D fhBgTrkNumber;
    AliJTH1D fhBgTrkNumberBin;


    //Randomized background histograms
    AliJTH1D fhBgRndmTrkPt;
    AliJTH1D fhBgRndmZ;
    AliJTH1D fhBgRndmJt;
    AliJTH1D fhBgRndmLogJt;
    AliJTH1D fhBgRndmJtWithPtCutWeightBin;
    AliJTH1D fhBgRndmLogJtWithPtCutWeight2Bin;
    AliJTH1D fhBgRndmJtWithPtCutWeightBinBin;
    AliJTH1D fhBgRndmLogJtWithPtCutWeight2BinBin;
    AliJTH1D fhBgRndmTrkNumber;

    AliJTH1D fhdeltaE;
    AliJTH1D fhdeltaN;
    AliJTH1D fhFullJetEChJetBin;
    AliJTH1D fhFullChdRChJetBin;
    AliJTH2D fh2DFullEvsChEdN0;
    AliJTH2D fh2DFullEvsChEdNnot0;
    AliJTH2D fhJetEtaPhi;
    AliJTH2D fhTrackEtaPhi;
    AliJTH1D fhJetdPt;

    //Jet Correlation Histograms
    AliJTH2D fhTrackJtCorrBin;	
    //AliJTH3D fhTrackJtCorr2D; //FIXME
    AliJTH2D fhTrackJtCorrBinTest;	
    AliJTH2D fhTrackJtCorrBinTest2;
    AliJTH2D fhTrackPtCorr;	
    AliJTH2D fhConstPtCorr;	
    AliJTH2D fhJetPtCorr;	
    AliJTH2D fhJetPtCorr2;	
    AliJTH2D fhJetPtCorr3;	
    AliJTH2D fhJetPtCorrCoarse;	
    AliJTH1D fhJetdR;
    AliJTH1D fhTrackMatchSuccess;
    AliJTH2D fhConstJtCorrBin;	
    AliJTH1D fhConstMatchSuccess;

    //PYthia correlation histograms
    AliJTH2D fhTrackJtCorrBinPythia;	
    AliJTH2D fhTrackPtCorrPythia;	
    AliJTH2D fhJetPtCorrPythia;	
    AliJTH2D fhJetPtCorr2Pythia;	
    AliJTH2D fhJetPtCorrPythiaCoarse;	
    AliJTH1D fhJetdRPythia;
    AliJTH1D fhTrackMatchSuccessPythia;
    AliJTH1D fhJetdPtPythia;

    //Pythia jets histograms
    AliJTH1D fhJetPtPythia;
    AliJTH1D fhJetPtBinPythia;
    AliJTH1D fhJtPythia;
    AliJTH1D fhJtBinPythia;
    AliJTH1D fhJtWeightBinPythia;
    AliJTH1D fhLogJtWeightBinPythia;
    AliJTH1D fhLogJtWeight2BinPythia;
    AliJTH1D fhJtWithPtCutWeightBinBinPythia;
    AliJTH1D fhLogJtWithPtCutWeightBinBinPythia;
    AliJTH1D fhLogJtWithPtCutWeight2BinBinPythia;

    //Pythia background
    AliJTH1D fhBgJtPythia;
    AliJTH1D fhBgJtBinPythia;
    AliJTH1D fhBgJtWeightBinPythia;
    AliJTH1D fhBgLogJtWeightBinPythia;
    AliJTH1D fhBgLogJtWeight2BinPythia;
    AliJTH1D fhBgJtWithPtCutWeightBinBinPythia;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinPythia;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinPythia;
};

#endif

