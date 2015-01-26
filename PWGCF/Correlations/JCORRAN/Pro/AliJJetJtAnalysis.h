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

class AliJEfficiency;
class TString;

class AliJJetJtAnalysis{
    public:


        AliJJetJtAnalysis();
        AliJJetJtAnalysis( AliJCard * card );
        AliJJetJtAnalysis(const AliJJetJtAnalysis& ap);
        AliJJetJtAnalysis& operator = (const AliJJetJtAnalysis& ap);
        ~AliJJetJtAnalysis();

        //type = charged jet or full jet, js = cone size 

        void FillHistosJets();

        void AddJets(TObjArray * jets ){ 
          if( !jets ) {
            //return;
          }
          fJetListOfList.Add( (TObject*)jets ); 
          //if( !jets ) return;
          for( int i=0;i<jets->GetEntriesFast();i++ ){
              //((AliJJet*)jets->At(i))->ReSum();
          }
        } // TODO clean before event

        void SetJTracks(TClonesArray *tracks ){fTracks = tracks ;}


        int GetNJets(){ return GetJetList()->GetEntriesFast(); }
        TObjArray* GetJetList(){ return fJetList; }
        //Double_t GetJetEtaRange(){ return fJetEtaRange; }
        //void SetJetEtaRange(double eta){ fJetEtaRange=eta; }
        void SetJetList(TObjArray* jetlist){ fJetList=jetlist; }
        void SetInputList(TObjArray * ilist){ fInputList = ilist;}
        //void SetTrackJetMap(std::vector<int> * v){ fTrackJetMap=v;}
        //void SetJetPtBin( TVector * b){ fJetPtBins=b; }
        void SetCard(AliJCard * card) {fCard = card;}

        AliJJet & Jet( int i ){ return *(AliJJet*) fJetList->At(i); }

        void UserCreateOutputObjects();
        void UserExec();

        void ClearBeforeEvent();
        void FillJtHistogram( TObjArray *Jets, int iContainer );

        void FillBgJtWithSmallerR(const TClonesArray &Jets, 
            double nR, int iHist);

        void FillBgJtWithDiffAxes (int iao, int ia, int iHist);

        void WriteHistograms();
        
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

		// Need for event loop
        void SetCentralityBin( int cbin) { cBin = cbin;}
        void SetCentrality( float cent) { fcent = cent;}
        void SetZVertex( float zvtx) { zVert = zvtx;}
        void SetZVertexBin( int zbin) { zBin = zbin;}
		void SetNumberOfJetFinders( int njfinder ) { nJetContainer = njfinder;}
        AliJEfficiency* GetAliJEfficiency() { return fEfficiency;}



    private:
        TObjArray * fInputList; // comment needed
        TObjArray * fJetList; // comment needed
        TObjArray fJetListOfList; // !comment needed
        vector<TClonesArray>      fJetBgListOfList;
        double fJetEtaCut;

        TVector  *fJetTriggPtBorders;
        TVector  *fJetConstPtLowLimits;
        TVector  *fJetAssocPtBorders;
        TVector  *fDeltaRBorders;
		int nJetContainer;

        AliJCard * fCard; // comment needed
        AliJJetAnalysis *fJJetAnalysis;
		vector<TString> fJetFinderName;
        vector<double>  fConeSizes;
		// Need for events
		AliJEfficiency *fEfficiency;
		int cBin;
		float fcent;
		int zBin;
		float zVert;
        TClonesArray *fTracks;

        //Histograms
        AliJHistManager * fHMG;

        AliJBin fJetFinderBin; 
        AliJBin fJetTriggerBin; 
        AliJBin fTrkPtBin; 
        AliJBin fTrkLimPtBin; 
        AliJBin fdRBin;
        AliJBin fiHist;
        AliJTH1D fhNumber;
        AliJTH1D fhKNumber;
        AliJTH1D fhJetPt ;
        AliJTH1D fhJetPtBin;
        AliJTH1D fhZ ;
        AliJTH1D fhZBin;
        AliJTH1D fhJt ;
        AliJTH1D fhJtBin;
        AliJTH1D fhJtWeightBin;
        AliJTH1D fhLogJtWeightBin;
        AliJTH1D fhJtWithPtCutWeightBinBin;
        AliJTH1D fhLogJtWithPtCutWeightBinBin;
        AliJTH1D fhJtBinLimBin;
        AliJTH1D fhJtWeightBinLimBin;
        AliJTH1D fhLogJtWeightBinLimBin;
        
        
        AliJTH1D fhJetBgPt ;
        AliJTH1D fhJetBgPtBin;
        AliJTH1D fhBgZ ;
        AliJTH1D fhBgZBin;
        AliJTH1D fhBgJt ;
        AliJTH1D fhBgJtBin;
        AliJTH1D fhBgJtWeightBin;
        AliJTH1D fhBgLogJtWeightBin;
        AliJTH1D fhBgJtWithPtCutWeightBinBin;
        AliJTH1D fhBgLogJtWithPtCutWeightBinBin;
        AliJTH1D fhBgJtWithPtCutWeightBinBinSmallerR;
        AliJTH1D fhBgLogJtWithPtCutWeightBinBinSmallerR;
        AliJTH1D fhBgJtWithPtCutWeightBinBinDiffR;
        AliJTH1D fhBgLogJtWithPtCutWeightBinBinDiffR;
        AliJTH1D fhBgJtBinLimBin;
        AliJTH1D fhBgJtWeightBinLimBin;
        AliJTH1D fhBgLogJtWeightBinLimBin;
        
        AliJTH1D fhdeltaE;
        AliJTH1D fhdeltaN;
        AliJTH1D fhFullJetEChJetBin;
        AliJTH1D fhFullChdRChJetBin;
        AliJTH2D fh2DFullEvsChEdN0;
        AliJTH2D fh2DFullEvsChEdNnot0;


		//double   fJetPtMinCut;
};

#endif

