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
    void FillJtHistogramMC( TObjArray *Jets, int iContainer );
    void FillRandomBackground(TObjArray *Jets , int iContainer, int MC);
    void FillRandomBackground(double jetpT, double jetE, TObjArray *Jets , int iContainer, int MC);

    void FillCorrelation( TObjArray *Jets, TObjArray *mcJets, int iContainer);

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

    TVector *fTrackJt;
    TVector *fTrackPt;
    TVector *fJetPt;
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
    AliJTH1D fhJetConeJtWeightWithTrackCutBinBin;
    AliJTH1D fhJetConeJtWeightWithMultiplicityCutBinBin;
    AliJTH1D fhJetConeLogJtWeightBin;
    AliJTH1D fhJetConeLogJtWeight2Bin;
    AliJTH1D fhJetConeJtWithPtCutWeightBinBin;
    AliJTH1D fhJetConeLogJtWithPtCutWeightBinBin;
    AliJTH1D fhJetConeLogJtWithPtCutWeight2BinBin;

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


    //double   fJetPtMinCut;
    //Monte Carlo Truth
    AliJTH1D fhNumberMC;
    AliJTH1D fhKNumberMC;
    AliJTH1D fhJetPtMC ;

    AliJTH1D fhJetPtBinMC;
    AliJTH1D fhJetPtMultiplicityCutBinMC;
    AliJTH1D fhJetPtTrackCutBinMC;
    AliJTH1D fhJetPtWeightMC;
    AliJTH1D fhJetPtWeightBinMC;
    AliJTH1D fhJetMultiplicityBinMC;
    AliJTH1D fhZMC;
    AliJTH1D fhZBinMC;
    AliJTH1D fhJtMC;
    AliJTH1D fhJtBinMC;
    AliJTH1D fhJtWeightBinMC;
    AliJTH1D fhLogJtWeightBinMC;
    AliJTH1D fhLogJtWeight2BinMC;
    AliJTH1D fhJtWithPtCutWeightBinBinMC;
    AliJTH1D fhLogJtWithPtCutWeightBinBinMC;
    AliJTH1D fhLogJtWithPtCutWeight2BinBinMC;
    AliJTH1D fhJtBinLimBinMC;
    AliJTH1D fhJtWeightBinLimBinMC;
    AliJTH1D fhLogJtWeightBinLimBinMC;
    AliJTH1D fhLogJtWeight2BinLimBinMC;

    //Histograms for jt in cone
    AliJTH1D fhJetConeTrkPtMC; //
    AliJTH1D fhJetConeTrkPtBinMC; //
    AliJTH1D fhJetConeTrkPtWeightBinMC; //
    AliJTH1D fhJetConeZMC;
    AliJTH1D fhJetConeZBinMC;
    AliJTH1D fhJetConeJtMC;
    AliJTH1D fhJetConeJtBinMC;
    AliJTH1D fhJetConeJtWeightWithTrackCutBinBinMC;
    AliJTH1D fhJetConeJtWeightWithMultiplicityCutBinBinMC;
    AliJTH1D fhJetConeJtWeightBinMC;
    AliJTH1D fhJetConeLogJtWeightBinMC;
    AliJTH1D fhJetConeLogJtWeight2BinMC;
    AliJTH1D fhJetConeJtWithPtCutWeightBinBinMC;
    AliJTH1D fhJetConeLogJtWithPtCutWeightBinBinMC;
    AliJTH1D fhJetConeLogJtWithPtCutWeight2BinBinMC;

    AliJTH1D fhJetBgPtMC;
    AliJTH1D fhJetBgPtBinMC;
    AliJTH1D fhBgZMC;
    AliJTH1D fhBgZBinMC;
    AliJTH1D fhBgJtMC;
    AliJTH1D fhBgJtBinMC;
    AliJTH1D fhBgJtWeightBinMC;
    AliJTH1D fhBgLogJtWeightBinMC;
    AliJTH1D fhBgLogJtWeight2BinMC;
    AliJTH1D fhBgJtWithPtCutWeightBinBinMC;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinMC;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinMC;
    AliJTH1D fhBgJtWithPtCutWeightBinBinSmallerRMC;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinSmallerRMC;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinSmallerRMC;
    AliJTH1D fhBgJtWithPtCutWeightBinBinDiffRMC;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinDiffRMC;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinDiffRMC;
    AliJTH1D fhBgJtBinLimBinMC;
    AliJTH1D fhBgJtWeightBinLimBinMC;
    AliJTH1D fhBgLogJtWeightBinLimBinMC;
    AliJTH1D fhBgLogJtWeight2BinLimBinMC;
    AliJTH1D fhTrkPtMC;
    AliJTH1D fhTrkPtBinMC;
    AliJTH1D fhTrkPtWeightBinMC;
    AliJTH1D fhLeadingTrkPtBinMC;
    AliJTH1D fhBgTrkPtMC;
    AliJTH1D fhBgTrkPtBinMC;
    AliJTH1D fhBgTrkPtWeightBinMC;
    AliJTH1D fhBgTrkNumberMC;
    AliJTH1D fhBgTrkNumberBinMC;


    //Randomized background histograms
    AliJTH1D fhBgRndmTrkPtMC;
    AliJTH1D fhBgRndmZMC;
    AliJTH1D fhBgRndmJtMC;
    AliJTH1D fhBgRndmLogJtMC;
    AliJTH1D fhBgRndmJtWithPtCutWeightBinMC;
    AliJTH1D fhBgRndmLogJtWithPtCutWeight2BinMC;
    AliJTH1D fhBgRndmJtWithPtCutWeightBinBinMC;
    AliJTH1D fhBgRndmLogJtWithPtCutWeight2BinBinMC;
    AliJTH1D fhBgRndmTrkNumberMC;

    AliJTH1D fhdeltaEMC;
    AliJTH1D fhdeltaNMC;
    AliJTH1D fhFullJetEChJetBinMC;
    AliJTH1D fhFullChdRChJetBinMC;
    AliJTH2D fh2DFullEvsChEdN0MC;
    AliJTH2D fh2DFullEvsChEdNnot0MC;
    AliJTH2D fhJetEtaPhiMC;
    AliJTH2D fhTrackEtaPhiMC;


    //Jet Correlation Histograms
    AliJTH2D fhTrackJtCorrBin;	
    AliJTH2D fhTrackPtCorr;	
    AliJTH2D fhJetPtCorr;	
    AliJTH2D fhJetPtCorr2;	
    AliJTH1D fhJetdR;
    AliJTH1D fhTrackMatchSuccess;
};

#endif

