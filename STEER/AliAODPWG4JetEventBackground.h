#ifndef AliPWG4AODJETBACKGROUND_H
#define AliPWG4AODJETBACKGROUND_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     AOD jet background class
//     Stores Different background calculations on an event by event level
//     Author: Christian Klein-Boesing, IKP Muenster
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TString.h"


class AliAODPWG4JetEventBackground : public TObject {

 public:
    AliAODPWG4JetEventBackground();
    virtual ~AliAODPWG4JetEventBackground();
    AliAODPWG4JetEventBackground(const AliAODPWG4JetEventBackground& jet); 
    AliAODPWG4JetEventBackground& operator=(const AliAODPWG4JetEventBackground& jet);


    virtual void SetBackground(Int_t i,Double_t back){
      fBackground[i] = back;
    }

    virtual Double_t GetBackground(Int_t i){
      if(i>=0&&i<kMaxBackground)return fBackground[i];
      return 0;
    }

    static const char* StdBranchName(){return fgkStdBranchName.Data();}
    virtual void       Print(Option_t* /*option*/) const;

    enum { kSmallR = 0,
	   kOnlyCharged,
	   kOutOfCone,
	   kStatistical,
	   kMaxBackground};

 private:
    static TString fgkStdBranchName;                    // Standard branch name
    Double32_t      fBackground[kMaxBackground];        // Background from different schemes, normalized to area
    ClassDef(AliAODPWG4JetEventBackground,1);

};


#endif
