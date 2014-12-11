/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
//===========================================================

#ifndef ALIJJETCORRANALYSIS_H
#define ALIJJETCORRANALYSIS_H

#include <vector>
#include <TObjArray.h>
#include <TClonesArray.h>
#include "AliJConst.h"
#include "AliJCard.h"
#include "AliJJet.h"
#include "AliJHistManager.h"

class AliJCard;
class AliJHistos;
class AliJCorrelations;
class AliJEfficiency;
class AliJEventPool;

class AliJJetCORRAnalysis{
    public:
        AliJJetCORRAnalysis();
        AliJJetCORRAnalysis( AliJCard * card );
        AliJJetCORRAnalysis( const AliJJetCORRAnalysis & obj );
		AliJJetCORRAnalysis& operator=(const AliJJetCORRAnalysis & obj);
		virtual ~AliJJetCORRAnalysis();
		void SetDebugMode( int debug) { fDebugMode = debug; };

		void CreateHistos();
		void FillHistosJets(TObjArray *Jets);

		void AddJets(TObjArray * jets ){ 
			if( !jets ) {
				cout<<"JWARN_C1 in AddJets jets="<<jets<<endl;
				//return;
			}
			fJetListOfList.Add( (TObject*)jets ); 
			if( !jets ) return;
			for( int i=0;i<jets->GetEntriesFast();i++ ){
				//((AliJJet*)jets->At(i))->ReSum();
			}
		} // TODO clean before event

		int GetNJets(){ return GetJetList()->GetEntriesFast(); }
		TObjArray* GetJetList(){ return fJetList; }
		Double_t GetJetEtaRange(){ return fJetEtaRange; }
		void SetJetEtaRange(double eta){ fJetEtaRange=eta; }
		void SetJetList(TObjArray* jetlist){ fJetList=jetlist; }
		void SetInputList(TObjArray * ilist){ fInputList = ilist;}
		void SetTrackJetMap(std::vector<int> * v){ fTrackJetMap=v;}
		void SetJetPtBin( TVector * b){ fJetPtBins=b; }
		void SetTargetJetIndex( int jfindex ) { fTargetJetIndex = jfindex; } // PlayCorrelation only for selected jetfinder

		AliJJet & Jet( int i ){ return *(AliJJet*) fJetList->At(i); }
		void UserCreateOutputObjects();
		void UserExec();
		void Terminate() const;

		void ClearBeforeEvent();

		// Related JCORRAN classes
		void PlayCorrelation(TClonesArray *triggList, TClonesArray *assocList);
		void SetCentralityBin( int cbin) { cBin = cbin;}
		void SetCentrality( float cent) { fcent = cent;}
		void SetZVertex( float zvtx) { zVert = zvtx;}
		void SetZVertexBin( int zbin) { zBin = zbin;}
		AliJEfficiency* GetAliJEfficiency() { return fEfficiency;}
		AliJHistos *GetAliJHistos() { return fHistos;}

	private:
		TObjArray * fInputList; // comment needed
		TObjArray * fJetList; // comment needed
		TObjArray fJetListOfList; // !comment needed
		bool   fIsFullAcceptance; // comment needed
		double fJetEtaRange; // comment needed
		std::vector<int> *fTrackJetMap; // comment needed

		TVector *fJetPtBins;
		double   fJetPtMinCut;

		AliJCard * fCard; // comment needed
		AliJHistos * fHistos;               //!
		AliJEfficiency *fEfficiency;
		TClonesArray * ftriggList;
		TClonesArray * fassocList;
		AliJCorrelations *fcorrelations;    //!
		AliJEventPool *fassocPool;          //!
		int cBin;
		float fcent;
		int zBin;
		float zVert;
		Int_t fevt; // event number
		int fTargetJetIndex;
		int fDebugMode;
};

#endif

