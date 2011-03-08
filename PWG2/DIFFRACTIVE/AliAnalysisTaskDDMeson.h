/*************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Select events accroding to L0 trigger input
// Reconstruct 2pi mass
// save charge, Armenteros' x and y,  pid, momentum
//

#ifndef ALIANALYSISTASKDDMESON_H
#define ALIANALYSISTASKDDMESON_H

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

class AliESDEvent; 
class AliESDtrack;

class TH1I;
class TH2I;
class TH2D;
class TList; 
class THnSparse;
class TLorentzVector;

class AliAnalysisTaskDDMeson : public AliAnalysisTaskSE{
  public:

    AliAnalysisTaskDDMeson(const TString opt);
    AliAnalysisTaskDDMeson(const AliAnalysisTaskDDMeson  &p);
    AliAnalysisTaskDDMeson& operator=(const AliAnalysisTaskDDMeson  &p);
    virtual ~AliAnalysisTaskDDMeson();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    void IniTask();

    Bool_t CheckESD();
    Bool_t CheckBit();
    Int_t CutESD(const AliESDtrack *outtrk[]);
    //Bool_t CutTrack(const AliESDtrack * esdtrack) const;

    //------
    void SwapTrack(const AliESDtrack * trks[]) const;

    Int_t GetV0() const;
    Int_t GetCombCh(const Double_t s1, const Double_t s2) const;

    TLorentzVector GetKinematics(const Double_t *pa, const Double_t *pb, const Double_t ma, const Double_t mb, Double_t & cts) const;
    Double_t GetCtlab(const Double_t *pa, const Double_t *pb) const;

    void CheckRange(Double_t &ptv, Double_t &pta, Double_t &etaa
                    , Double_t &mpi
                    ) const;

    void FillBit();
    void CalcBit(TH1I *hc, Double_t tot[]) const;

    //---
    void SPDLoadGeom() const;
    Bool_t SPDLoc2Glo(const Int_t id, const Double_t *loc, Double_t *glo) const;
    Bool_t CheckChipEta(const Int_t chipKey) const;
    void GetNFO(Int_t &ni, Int_t &no) const;

 private:
    TString fOpt;                            //option
    AliESDEvent *fESD;                          //esd event                    
    //------

    Int_t fnmass;                          //nbins for mass
    Double_t fmass1;                       //upper edge of axis   

    Int_t fnptv;                             //nbins for p
    Double_t fptv1;                         //upper edge of axis  

    Int_t fnpta;                             //nbins for p
    Double_t fpta1;                         //upper edge of axis  

    Int_t fneta;                           //nbins for eta
    Double_t feta;                        //upper edge of axis

    Int_t fnnsel;                           //nbins for nsel
    Double_t fnsel1;                        //upper edge of axis

    Int_t fncts;                            //nbins for cts
    Int_t fnctlab;                           //nbins for ctlab

    //------
    Int_t fCHECKVBA;                       //V0A bit
    Int_t fCHECKVBC;                       //V0C bit

    TH1I *fHBIT;                                //histogram of bits
    Int_t fBitcg;                                 //trigger bit configuration
 
    Int_t fRun;                                   //run number
    TH1I *fat;                                  //V0A-only
    TH1I *fct;                                  //V0C-only
    TH1I *fbt;                                  //V0A & V0C
    TH1I *fnt;                                  //!V0A & !V0C
    TH1I *ftt;                                  //TOTAL

    //------
    TH2D *fv0ntrk;                               //v0bit vs. nch
    TH2D *frsntrk;                               //raw vs. sel

    TH1I *fhps;                                 //V0 BG
    TH2I *fhfo; //SPD fastor hardware, in eta acceptance
    TH1I *fhspd;//SPD fastor offline
    TH2I *fhv0fmd; //v0 vs fmd
    TH1I *fhpriv;                               //primary vertex cut effect
    TH1I *fhntrk;                               //n-trk-after-cut effect

    //------
    TList *fHist;                             //output list
    THnSparse *fThnMass;                    //ThnSparse for mother pt and Mass
    THnSparse *fThnDPt;                    //ThnSparse for pt and eta of daughter
    THnSparse *fThnDEta;                    //ThnSparse for pt and eta of daughter
    THnSparse *fThnKF; //kf
    ClassDef(AliAnalysisTaskDDMeson, 1);
};


#endif

