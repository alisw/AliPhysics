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
        double Mt(){ return (jet(0)+jet(1)).Mt(); }
        double Mt0(){ return sqrt( pow(jet(0).Pt()+jet(1).Pt(),2) - pow(PtPair(),2)); }
        double Mt1(){ return sqrt( pow(sqrt(jet(0).M2()+jet(0).Perp2())+sqrt(jet(1).M2()+jet(1).Perp2()),2) - pow(PtPair(),2)); }
        double Mt2(){ return sqrt( pow(sqrt(jet(0).E2nd2()-jet(0).Pz()*jet(0).Pz())+sqrt(jet(1).E2nd2()-jet(1).Pz()*jet(1).Pz()),2) - pow(PtPair(),2)); }


        AliJJet& jet(int i){ Safe(i);return (*this)(i); }
        AliJJet& operator()(int i ){
            if( i > 1 || i< 0) { cout<<"wrong index"<<endl;exit(2);  }
            if( !fJet[i] ){ cout<<"Empty jet"<<endl;exit(3); }
            return* fJet[i];
        }
        void Print(int type){
            static int count = 0;
            if( GetEntries() < 2 ){
                cout<<"number of dijet is "<<GetEntries()<<endl;
                return;
            }
            TLorentzVector lvsum = jet(0)+jet(1);
            TLorentzVector * lv[3];
            lv[0] = &jet(0);
            lv[1] = &jet(1);
            lv[2] = &lvsum;

            cout<<"DEBUG_DiJet_Mass ==== START"<<endl;
            for( int i=0;i<3;i++ ){
                TLorentzVector &l = *lv[i];
                cout<<Form("DEBUG_DiJet_Mass : %d %d %d : %10.3f %10.3f %10.3f %10.3f %10.3f",count, type, i,l.Px(),l.Py(),l.Pz(),l.E(), l.M())<<endl;
            }
            lv[2]->Print();
            cout<<"DEBUG_DiJet_Mass === EMD"<<endl;
            count++;
        }

        vector<AliJJet*> fJet;
};


class AliJDiJetAnalysis{
    public:
        enum { 
            kJDiJetLeadingSubLeading5,
            kJDiJetLeadingSubLeading10,
            kJDiJetLeadingSubLeading20,
            kJDiJetLeadingSubLeading50,
            kLeadingSubLeadingOpposite,
            kLeadingSubLeadingEtaWindow,
            //kJDiJetAtlas,
            kJNDiJetSelection,
            kJDiJetMarta,
            kJDiJetCMS,
        };
        enum { 
            kJIncomingJet, kJOutgoingJet, 
            kJFullJetR03, kJFullJetR04, kJFullJetR05, kJFullJetR06, 
            kJChargedJetR03, kJChargedJetR04, kJChargedJetR05, kJChargedJetR06,
            //kJIdealJet, 
            kJMultiPartonAll, 
            kJNJetType 
        };
        enum{ kJEtaAll, kJEtaAlice, kJNJetSelection};
        enum{ kJNJetTypeMax=12 }; // FIXME:

        AliJDiJetAnalysis();
        AliJDiJetAnalysis( AliJBaseCard * card );
        AliJDiJetAnalysis( const AliJDiJetAnalysis & obj );
        AliJDiJetAnalysis& operator=(const AliJDiJetAnalysis & obj);

        AliJDiJet * GetDiJets(int type, int js, int djs);
        TObjArray * GetJets(int type, int js);

        void FillHistosJets();
        void FillHistosDiJet();
        void CreateHistos();

        static bool CompareTrackByPt( AliJBaseTrack * a, AliJBaseTrack* b );
        vector<AliJBaseTrack*> SortTrackByPt( TObjArray * os );


        void AddJets(TObjArray * jets, int isTrackOrParticle ){ 
            if( !jets ) {
                cout<<"JWARN_C1 in AddJets jets="<<jets<<endl;
                //return;
            }
            fJetListOfList.Add( (TObject*)jets ); 
            fTrackOrMCParticle.push_back( isTrackOrParticle );
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
        std::vector<int> fTrackOrMCParticle;
        std::vector<int> fChargedOrFull;
        AliJDiJet *fDiJets[kJNJetTypeMax][kJNJetSelection][kJNDiJetSelection]; // comment needed
        TObjArray *fJets[kJNJetTypeMax][kJNJetSelection]; // comment needed
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
        AliJTH1D fhJetInvM;
        AliJTH1D fhJetMjjEbin;
        AliJTH1D fhJetMtEbin;

        AliJTProfile fhJetEnergyComp;
        AliJTProfile fhJetPtComp;
        AliJTH1D fhJetDRToRef;
        AliJTH1D fhJetDPhiToRef;
        //== DiJets;
        AliJTH1D fhDiJetPtPair;
        AliJTH1D fhDiJetPtPairRaw;
        AliJTH1D fhDiJetPt1;
        AliJTH1D fhDiJetPt2;
        AliJTH1D fhDiJetInvM;
        AliJTH1D fhDiJetInvMInclu;
        AliJTH1D fhDiJetMtInclu;
        AliJTH1D fhDiJetKtA;

        AliJTH1D fhDiJetSingleJetMass;
        AliJTH1D fhDiJetSingleJetArea;
        AliJTH1D fhDiJetSingleJetActivity;
        AliJTH1D fhDiJetSingleJetNCont;

        AliJTH1D fhDiJetDeltaR;     //!
        AliJTH1D fhDiJetDeltaPhi;
        AliJTH1D fhDiJetDeltaEta;
        AliJTH1D fhDiJetPtAsymm;
        AliJTH1D fhDiJetEAsymm;
        AliJTH1D fhDiJetMultiplicity;

        AliJTH2D fhDiJetTypeCor;
        AliJTH2D fhInvMPttCor;
        AliJTH2D fhMtMjjCor;

        //== PYTHIA
        AliJTH1D fhPythiaJetPtPair;
        AliJTH1D fhPythiaJetSum;


        int fDiJetBin[kJNJetTypeMax][kJNJetSelection][kJNDiJetSelection][10];
        double fDiJetMass[kJNJetTypeMax][kJNJetSelection][kJNDiJetSelection][10];
};

#endif

