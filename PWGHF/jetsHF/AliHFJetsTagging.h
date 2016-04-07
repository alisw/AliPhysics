#ifndef ALIHFJETSTAGGING_H
#define ALIHFJETSTAGGING_H

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

/* Manager class for HF jet analysis   */

/* Mailto: andrea.rossi@cern.ch,       *
 *         elena.bruna@cern.ch,        *
 *         min.jung.kweon@cern.ch,     *
 *         linus.feldkamp@cern.ch,     *
 *         svallero@to.infn.it         *
 *         s.lapointe@cern.ch          */

class AliAODEvent;
class AliVEvent;
class AliVTrack;
class AliAODTrack;
class AliAODMCParticle;
class TClonesArray;
class AliEmcalJet;
class TVector2;
class AliMCEvent;
class AliMCParticle;

#include <TNamed.h>
#include <vector>

#include "AliHFJetsUtils.h"

class AliHFJetsTagging : public TNamed {
public:
  AliHFJetsTagging();
  AliHFJetsTagging(const char* name);
  ~AliHFJetsTagging();

  AliHFJetsTagging& operator=(const AliHFJetsTagging& corr);  // Assignment operator
  AliAODMCParticle* IsMCJetParton(const TClonesArray *arrayMC,const AliEmcalJet *jet,Double_t radius=0.7);
  AliMCParticle* IsMCJetParton(const AliMCEvent *mcevent,const AliEmcalJet *jet,Double_t radius=0.7);

  AliAODMCParticle* IsMCJetMeson(const TClonesArray *arrayMC,const AliEmcalJet *jet ,Double_t radius=0.7);
  AliMCParticle* IsMCJetMeson(const AliMCEvent *mcEvent,const AliEmcalJet *jet ,Double_t radius=0.7);

  AliAODMCParticle* GetAODMCParticleFromAodTrack(const TClonesArray *arrayMC,const AliAODTrack* track);

  Bool_t IsBMeson(Int_t pc);   
  Bool_t IsDMeson(Int_t pc); 

  Bool_t IsEventSelectedTrackCounting(const AliAODEvent*event);

  Bool_t GetSignedRPhiImpactParameter(const AliVEvent*event, const AliVTrack* track, const AliEmcalJet* jet,Double_t &signedimpactparameter, Double_t d[2],Double_t cov[3]);

  std::vector<std::pair<Double_t,Int_t> > GetSortedListOfSignedRPhiImpactParameters(const AliVEvent*event, const AliEmcalJet* jet,const TClonesArray * tracksAccepted);


  
protected:  
  Double_t RelativePhi(Double_t mphi,Double_t vphi);
private:
  struct sort_descend
  { // sort in decreasing order
    bool operator () (const std::pair<Double_t, Int_t>& p1, const std::pair<Double_t, Int_t>& p2)  { return p1.first > p2.first ; }
  };
  
  TClonesArray *fSelectedTracks;   //! array with selected tracks to be used for tagging
  
  ClassDef(AliHFJetsTagging,1);

};
#endif
