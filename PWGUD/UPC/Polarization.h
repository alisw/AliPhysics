/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef Polarization_H
#define Polarization_H

#include "AliAnalysisTaskSE.h"

class Polarization : public AliAnalysisTaskSE  
{
    public:
                                Polarization();
                                Polarization(const char *name);
        virtual                 ~Polarization();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;  //! input event
        AliAODZDC*              fZDCdata; 
        
        AliPIDResponse*         fPIDResponse;            
        TList*                  fOutputList;    //! output list
        //One Dimensional Counters
        TH1I*                   fHistITS;
        TH1I*                   fHistEtaCut;
        TH1I*                   fHistPtCut;
        TH1I*                   fHistPairCut;
        TH1I*                   fHistChargeCut;
        TH1I*                   fHistPIDCut;
        TH1I*                   fHistFilterbitCut;
        TH1I*                   fHistTriggerCut;
        
        TH1I*                   fHistCounter; 
        TH1F*                   fHistPt;        //! dummy histogram
        TH2F*                   fHistP_TPC;
                          
        Float_t                 fPt;
        Float_t                 fM;
        Float_t                 fPt0;
        Float_t                 fPt1;
       // Float_t                 fPt1;
        Float_t                 fPiSigma0;
        Float_t                 fPiSigma1;
        Float_t                 fMuSigma0;
        Float_t                 fMuSigma1;
        Float_t                 fRabs1;
        Float_t                 fRabs2;
        Float_t                 fKaonSigma0;
        Float_t                 fKaonSigma1;
        Float_t                 fTheta;
        Float_t                 fHelicityTheta;
        Float_t                 fCollinTheta;
        Float_t                 fPhi;
        Float_t                 fHelicityPhi;
        Float_t                 fCollinPhi;
       
        
        Float_t                 fTPCcluster1;
        Float_t                 fEta1;
        Float_t                 fTPCcluster2;
        Float_t                 fEta2;
        Float_t                 fDCAxy1;
        Float_t                 fDCAz1;
        Float_t                 fDCAxy2;
        Float_t                 fDCAz2;
        TTree*                  fTree;
        Int_t                   fTriggerClass;
        Float_t                 fPd;
        Float_t                 fPp;
        Float_t                 fdEdX;
        Float_t                 fPtd;
         
        Float_t                 fZNAenergy; 
        Float_t                 fZNCenergy;     
         
        Float_t                 fZDCAtime;
        Float_t                 fZDCCtime;
        
        
        
         TLorentzVector daughter1;
         TLorentzVector daughter2; 
         TLorentzVector parent;
   
        Polarization(const Polarization&); // not implemented
        Polarization& operator=(const Polarization&); // not implemented

        ClassDef(Polarization, 1);
};

#endif
