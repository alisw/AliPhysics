/**************************************************************************
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
// Get improved dca info. by ITS upgrade implemented by the ALICE Heavy Flavour Electron Group
//
#ifndef ALIHFESMEARDCA_H
#define ALIHFESMEARDCA_H

class TList;
class TGraph;

class AliVEvent;
class AliVParticle;
class AliVTrack;
class AliMCEvent;

class AliHFEsmearDCA : public TObject{
  public:
  AliHFEsmearDCA(const Char_t */*name*/, const char *resfileCurURI, const char *resfileUpgURI, const Char_t */*title*/);
    AliHFEsmearDCA(const AliHFEsmearDCA &c);
    virtual ~AliHFEsmearDCA();
    
    virtual void SetRecEventInfo(const TObject *event);

    void SetMCEvent(AliMCEvent* const mcEvent){fMCEvent=mcEvent;};  // set stack pointer
    void GetImproveITSImpactParameters(AliVTrack *track, Double_t &dcaxyn, Double_t &dcaxyo, Double_t &dcaxysign, Double_t &dcaxysigo, Double_t &dcazn, Double_t &dcazo, Double_t &dcazsign, Double_t &dcazsigo); // to check improvement by the ITS upgrade
    Double_t EvalGraph(const TGraph *graph,Double_t x) const;
    

  private:
    AliHFEsmearDCA &operator=(const AliHFEsmearDCA &);
 
    AliVEvent *fEvent;                //! working event
    AliMCEvent *fMCEvent;             //! MCEvent pointer

    TGraph *fD0ZResPCur  ; // old pt dep. d0 res. in z for protons
    TGraph *fD0ZResKCur  ; // old pt dep. d0 res. in z for kaons
    TGraph *fD0ZResPiCur ; // old pt dep. d0 res. in z for pions
    TGraph *fD0ZResECur  ; // old pt dep. d0 res. in z for electrons 
    TGraph *fD0RPResPCur ; // old pt dep. d0 res. in rphi for protons
    TGraph *fD0RPResKCur ; // old pt dep. d0 res. in rphi for kaons
    TGraph *fD0RPResPiCur; // old pt dep. d0 res. in rphi for pions
    TGraph *fD0RPResECur ; // old pt dep. d0 res. in rphi for electrons
    TGraph *fPt1ResPCur  ; // old pt dep. 1/pt res. for protons
    TGraph *fPt1ResKCur  ; // old pt dep. 1/pt res. for kaons
    TGraph *fPt1ResPiCur ; // old pt dep. 1/pt res. for pions
    TGraph *fPt1ResECur  ; // old pt dep. 1/pt res. for electrons
    TGraph *fD0ZResPUpg  ; // new pt dep. d0 res. in z for protons
    TGraph *fD0ZResKUpg  ; // new pt dep. d0 res. in z for kaons
    TGraph *fD0ZResPiUpg ; // new pt dep. d0 res. in z for pions
    TGraph *fD0ZResEUpg  ; // new pt dep. d0 res. in z for electrons
    TGraph *fD0RPResPUpg ; // new pt dep. d0 res. in rphi for protons
    TGraph *fD0RPResKUpg ; // new pt dep. d0 res. in rphi for kaons
    TGraph *fD0RPResPiUpg; // new pt dep. d0 res. in rphi for pions
    TGraph *fD0RPResEUpg ; // new pt dep. d0 res. in rphi for electrons
    TGraph *fPt1ResPUpg  ; // new pt dep. 1/pt res. for protons
    TGraph *fPt1ResKUpg  ; // new pt dep. 1/pt res. for kaons
    TGraph *fPt1ResPiUpg ; // new pt dep. 1/pt res. for pions
    TGraph *fPt1ResEUpg  ; // new pt dep. 1/pt res. for electrons
  
  ClassDef(AliHFEsmearDCA, 1)      // Additional cuts implemented by the ALICE HFE group
};

#endif
