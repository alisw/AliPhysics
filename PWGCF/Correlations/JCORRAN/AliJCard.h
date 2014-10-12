//=================================================
// AliJCard.h
// last modified FK 6.NOV 2009
//=================================================

#ifndef ALIJCARD_H
#define ALIJCARD_H

#include <TObject.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include <string.h>
#include <TString.h>
#include <TVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <TVector3.h>
#include <THashList.h>
#include <TNamed.h>

#include "AliJConst.h"
#include "AliJBaseCard.h"

class AliJBaseTrack;
#include "AliJPhoton.h"
class AliJTrack;

using namespace std;


class AliJCard : public AliJBaseCard {

    public:

        AliJCard(); // constructor
        AliJCard(const char *filename); // constructor
        AliJCard(const AliJCard& obj);
        AliJCard& operator=(const AliJCard& obj);

        virtual ~AliJCard();

        void MakeFastCorrTypeIndex();

        void   PrintOut(); 

        float  Get(TString keyword, int VectorComponent=0){ return AliJBaseCard::Get(keyword, VectorComponent); }
        float  Get(corrType ctype, int VectorComponent =0);  //get TVector component
        float  GetFast(corrType ctype, int VectorComponent=0);  //get TVector component
        int    GetN(TString keyword){ return AliJBaseCard::GetN(keyword); }       //get TVector dimension
        int    GetN(corrType ctype);        //get TVector dimension
        int    GetNFast(corrType ctype);    //get TVector dimension

        //---- Collision Species --
        float  GetCollisionType()  { return Get("CollisionType");}

        //------  v e r t e x -----

        float  VertInZRange(float Z) {
            //cout<< "zv " << Z <<" "<< ( GetBinBorder(kZVertType,0) < Z && Z < GetBinBorder(kZVertType, GetNoOfBins(kZVertType)) ) <<endl; 
            return ( GetBinBorder(kZVertType,0) < Z && Z < GetBinBorder(kZVertType, GetNoOfBins(kZVertType)) ); }

        //Alice vertex  cuts
        bool CheckEventZVetrex(double fZVertex, double ZVertErr){
            if(Get("ZVertexRange",0)<fZVertex && fZVertex <Get("ZVertexRange",1) && (Get("MaxZVertexError") > ZVertErr))        return true; 
            else    return false;
        }

        int IsLessThanUpperPairPtCut(double inPairPt);


        //--- c o r r e l a t i o n  bins & borders --
        int    GetNoOfBins (corrType ctype){ return GetNFast(ctype)-1; }
        float  GetBinBorder(corrType ctype, int ii){ return GetFast(ctype,ii); }
        int    GetBin(corrType ctype, float val);
        int    GetBinFast(corrType ctype, float val);
        double GetBinCenter(corrType ctype, int ii){ return (GetFast(ctype,ii)+GetFast(ctype,ii+1))/2.;}

        //-----  m i x i n g ----
        int    GetEventPoolDepth(int cBin){ return (int) Get("EventPoolDepth",cBin);}
        bool   SimilarVertZ(float Z1, float Z2);
        bool   SimilarMultiplicity(float mult1, float mult2);
        bool   SimilarCentrality(float c1, float c2, int cbin);

        //run characteristics
        bool   IsGoodRun(int runID);

        //trigger
        bool   MbTrigger(int triggin) const;

        //photon
        bool   IsPhoton(AliJPhoton *g){ return g->GetProbPhot() > Get("probPhot"); }
        float  GetPhotEnergyCut(){ return Get("minEnergy"); }

        bool   InPhiRange(float Phi);
        bool   IsInZEDandThetaRange(float zedDC, float theta){
            return (fabs(zedDC)<Get("zedRange") && fabs(kJPi/2-theta)<Get("thetaRange"));
        }

        bool   IsInEtaRange(float eta){return fabs(eta)<Get("EtaRange"); }
        //  bool   likeSgnCheck(PhJCgl *cgl1, PhJCgl *cgl2);
        bool   DeltaEtaCheck(const AliJBaseTrack *ftk1, const AliJBaseTrack *ftk2);
        bool   NotGhost(float zedDC, float PhiDC){ 
            return (zedDC > Get("deltaZEDPhiDC",0) && fabs(PhiDC) > Get("deltaZEDPhiDC",1));
        }

        /* bool CheckCovDiagonalElements(double *element){
           if(Get("MaxCovDiagonalElements",0)<0) return true;
           bool isGoodTrack = true;
           for(Int_t i=0;i<5;i++){
           if(Get("MaxCovDiagonalElements",i) < element[i]) isGoodTrack = false;
           }
           return isGoodTrack; 
           }*/

        bool CheckTrackParamsInTPC(int NClustersTPC,float Chi2PerClusterTPC);

        bool CheckMinNumTPCClustPt(int NClustersTPC, float fpt);

        bool CheckTrackImpact(float xyIm, float zIm, float fpt);


        bool AcceptKinkDaughters(int kinkIndex){
            if(( ! (bool) Get("AcceptKinkDaughters")) && kinkIndex > 0 ) return false;//we do not want kinks but kink findex>0
            else return true;
        }


        double GetCutOnBkgActivity(){
            return (double) Get("CutOnBkgActivity");
        }

        void   SetEventV3kv(double inV3kv) {feventV3kv = inV3kv;}
        double GetEventV3kv() const {return feventV3kv;}

        //Alice CALO
        bool ReadEmcalSm(int sm){ return (bool) Get("EMCAL_SM",sm);}
        bool ReadPhosSm(int sm){ return (bool) Get("PHOS_SM",sm);}

        bool  MixPhotonForPi0Mass() {return ((int) Get("massMix")==1);}//FK//
        bool  MixMBForPi0Mass() {return ((int) Get("massMixMB")>0);}
        bool  MixMBMBForPi0Mass() {return ((int) Get("massMixMB")==2);}

        virtual void InitCard(); //TODO
        void FinishCard(); // TODO
        void ReCompile();


    protected:


        TString GetKeyWord(corrType ctype);
        corrType GetCorrType( TString inStr );

        //====   D a t a    M e m b e r s  ========

        //   double effPar[16];
        //   double corrCent[10]; //FK// additional scaling factor to effPar to correct on hi fcentrality
        TH2D *fhCorr;  // comment me

        double feventV3kv;  // comment me

        vector< int >       fIndexVector;     //array of float number confg parameter vectors 

        Double_t **fpi0massbin; //! faster access to pi0 mass bins

        //ClassDef(AliJCard, 1); // EMCAL for jcorran
};

#endif






















