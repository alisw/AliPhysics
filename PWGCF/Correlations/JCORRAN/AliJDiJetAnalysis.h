/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
//===========================================================

#ifndef ALIJDIJETANALYSIS_H
#define ALIJDIJETANALYSIS_H

#include <vector>
#include <TObjArray.h>
#include "AliJBaseCard.h"
#include "AliJConst.h"
#include "AliJJet.h"
//#include "AliJHistosDiJetAnalysis.h"
#include "AliJHistManager.h"

class AliJDiJet {
    public:
        AliJDiJet():fJet(){ ; }
        AliJDiJet(AliJJet*i,AliJJet*j):fJet(){SetJets(i,j);}
        void SetJet( int i, AliJJet* j ){ Safe(i);fJet[i]=j; }
        void SetJets( AliJJet* i, AliJJet*j){ Safe(2);fJet[0]=i;fJet[1]=j; }
        double DPt( AliJDiJet & dijet, int i, int j );
        double PtPair(){ return (jet(0)+jet(1)).Pt(); }
        double kTA(){ return LeadingPt()*TMath::Sin( DPhi() ); }
        double DPhi(){ return jet(0).DeltaPhi(jet(1)); }
        double DR(){ return jet(0).DeltaR(jet(1)); }
        double DEta(){ return jet(0).Eta() - jet(1).Eta(); }
        double DPt(){ return jet(0).Pt() - jet(1).Pt();}
        double PtAsymm(){ return  (jet(0).Pt() -jet(1).Pt())/(jet(0).Pt()+jet(1).Pt()); }
        double EAsymm(){ return  (jet(0).E() -jet(1).E())/(jet(0).E()+jet(1).E()); }
        double InvM(){ return (jet(0)+jet(1)).M(); }
        void Clear(){ fJet.clear(); }
        void Add(AliJJet*j){ fJet.push_back(j); } 
        int GetEntries(){ return fJet.size(); }
        void Safe( UInt_t x ){ if( fJet.size()<x) fJet.resize(x); }
        double LeadingPt(){ double pt0=jet(0).Pt();double pt1=jet(1).Pt();return pt0>pt1?pt0:pt1;}

        AliJJet& jet(int i){ Safe(i);return (*this)(i); }
        AliJJet& operator()(int i ){
            if( i > 1 || i< 0) { cout<<"wrong index"<<endl;exit(2);  }
            if( !fJet[i] ){ cout<<"Empty jet"<<endl;exit(3); }
            return* fJet[i];
        }

        vector<AliJJet*> fJet;
};


class AliJDiJetAnalysis{
    public:
        enum { 
            kJDiJetLeadingSubLeading,
            kLeadingSubLeadingOpposite,
            kJDiJetMarta,
            //kJDiJetAtlas,
            kJNDiJetSelection
        };
        enum { 
            kJIncomingJet, kJOutgoingJet, 
            kJFullJet, kJFullJet08, 
            kJChargedJet, kJChargedJet08, 
            //kJIdealJet, 
            kJMultiPartonAll, 
            kJNJetType 
        };
        enum{ kJEtaAll, kJEtaAlice, kJNJetSelection};

        AliJDiJetAnalysis();
        AliJDiJetAnalysis( AliJBaseCard * card );
        AliJDiJetAnalysis( const AliJDiJetAnalysis & obj );
        AliJDiJetAnalysis& operator=(const AliJDiJetAnalysis & obj);

        AliJDiJet * GetDiJets(int type, int js, int djs);
        TObjArray * GetJets(int type, int js);

        void FillHistosJets();
        void FillHistosDiJet();
        void CreateHistos();


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

        AliJJet & Jet( int i ){ return *(AliJJet*) fJetList->At(i); }
        TObjArray * GetDiJetList( int i ){ return (TObjArray*)fJetListOfList[i]; }
        void FillJetPythiaComperisonHistos(int ref, int typ, int js);

        void UserCreateOutputObjects();
        void UserExec();
        void Terminate() const;

        void ClearBeforeEvent();

        void FillPythiaDiJet(int js );

    private:
        TObjArray * fInputList; // comment needed
        TObjArray * fJetList; // comment needed
        TObjArray fJetListOfList; // !comment needed
        AliJDiJet *fDiJets[kJNJetType][kJNJetSelection][kJNDiJetSelection]; // comment needed
        TObjArray *fJets[kJNJetType][kJNJetSelection]; // comment needed
        bool   fIsFullAcceptance; // comment needed
        double fJetEtaRange; // comment needed
        int   fIsMartasDiJet; // comment needed
        std::vector<int> *fTrackJetMap; // comment needed

        TVector *fJetPtBins;
        TVector *fJetPtPairBins;
        double   fJetPtMinCut;



        double fMeanDiJetPt;

        // PYTHIA8
        double fParton23Pt;

        AliJBaseCard * fCard; // comment needed


        //== Histogram
        AliJHistManager * fHMG;
        AliJBin fJetPtBin;
        AliJBin fJetPtPairBin;
        AliJBin fInvMBin;
        AliJBin fJetTypeBin;
        AliJBin fJetSelectionBin;
        AliJBin fDiJetSelectionBin;
        AliJBin fJetRBin;
        AliJBin fJetDRBin;
        AliJBin fDiJetBinTypeBin;
        AliJBin fPYJetTypeBin;
        AliJBin fBin2;

        //== Jets
        AliJTH1D fhJetPt ;
        AliJTH1D fhJetPhi;
        AliJTH1D fhJetEta;
        AliJTH1D fhJetNConst;
        AliJTH1D fhJetEnergy;

        AliJTProfile fhJetEnergyComp;
        AliJTProfile fhJetPtComp;
        AliJTH1D fhJetDRToRef;
        AliJTH1D fhJetDPhiToRef;
        //== DiJets;
        AliJTH1D fhDiJetPtPair;
        AliJTH1D fhDiJetInvM;
        AliJTH1D fhDiJetKtA;
        AliJTH1D fhDiJetDeltaR;     //!
        AliJTH1D fhDiJetDeltaPhi;
        AliJTH1D fhDiJetDeltaEta;
        AliJTH1D fhDiJetPtAsymm;
        AliJTH1D fhDiJetEAsymm;
        AliJTH1D fhDiJetMultiplicity;

        //== PYTHIA
        AliJTH1D fhPythiaJetPtPair;
        AliJTH1D fhPythiaJetSum;

        int fDiJetBin[kJNJetType][kJNJetSelection][kJNDiJetSelection][10];
};

#endif

