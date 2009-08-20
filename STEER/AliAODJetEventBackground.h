#ifndef AliAODJETBACKGROUND_H
#define AliAODJETBACKGROUND_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     AOD jet background class
//     Stores Different background calculations on an event by event level
//     Author: Christian Klein-Boesing, IKP Muenster
//-------------------------------------------------------------------------

#include "TNamed.h"
#include "TString.h"


class AliAODJetEventBackground : public TNamed {

 public:
    AliAODJetEventBackground();
    virtual ~AliAODJetEventBackground();
    AliAODJetEventBackground(const AliAODJetEventBackground& jet); 
    AliAODJetEventBackground& operator=(const AliAODJetEventBackground& jet);


    virtual void SetBackground(Int_t i,Double_t back){
      fBackground[i] = back;
    }

    virtual Double_t GetBackground(Int_t i){
      if(i>=0&&i<kMaxBackground)return fBackground[i];
      return 0;
    }

    static const char* StdBranchName(){return fgkStdBranchName.Data();}
    virtual void       Print(Option_t* /*option*/) const;
    virtual void       Reset();

    enum { kSmallR = 0,
	   kOnlyCharged,
	   kOutOfCone,
	   kStatistical,
	   kMaxBackground};

 private:
    static TString fgkStdBranchName;                    // Standard branch name
    Double32_t      fBackground[kMaxBackground];        // Background from different schemes, normalized to area
    ClassDef(AliAODJetEventBackground,1);

};


#endif
