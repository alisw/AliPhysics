#ifndef ALICMEANALYSISCUTS_H
#define ALICMEANALYSISCUTS_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2015                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/



#include <TNamed.h>
#include <TClonesArray.h>
#include <TBits.h>
#include <TMath.h>
#include <TArrayS.h>
#include <TString.h>
#include "iostream"



//_____________________________________________________________________
class AliCMEAnalysisCuts : public TNamed {

  enum Constants {
    kNCuts=100
  };

 public:
  AliCMEAnalysisCuts();
  AliCMEAnalysisCuts(const Char_t* name);
  AliCMEAnalysisCuts(const AliCMEAnalysisCuts &c);
  ~AliCMEAnalysisCuts();

  // setters
  void AddCut(Int_t type, Float_t min, Float_t max, Bool_t excludeRange=kFALSE) {fCuts[fNcuts][0]= type; fCuts[fNcuts][1]=min; fCuts[fNcuts][2]=max; fCuts[fNcuts][3]= (excludeRange ? 1.0 : 0.0); fNcuts++;}
  void AddFlag(Int_t type, Int_t flag, Bool_t accept) {fCuts[fNcuts][0]= type; fCuts[fNcuts][1]=-999.0; fCuts[fNcuts][2]=flag; fCuts[fNcuts][3]= (accept ? 1.0 : 0.0); fNcuts++;}
  void CopyCuts(AliCMEAnalysisCuts* cuts);
  //void SetName(TString name) {fCutsName=name;}

  // getters
  Int_t Type(Int_t cut) const {return ((Int_t) (fCuts[cut][0]+10e-6));}
  Float_t Min(Int_t cut) const {return fCuts[cut][1];}
  Float_t Max(Int_t cut) const {return fCuts[cut][2];}
  Bool_t ExcludeRange(Int_t cut) const {return (fCuts[cut][3]<0.5 ? kFALSE : kTRUE);}
  Int_t Ncuts() const {return fNcuts;}
  Bool_t  IsSelected(Float_t* values) const;
  //TString GetCutsName() const {return fCutsName;}


 private:

  Float_t fCuts[kNCuts][4];
  Int_t fNcuts;
  //TString fCutsName;


  ClassDef(AliCMEAnalysisCuts, 1);
};



//_____________________________________________________________________
inline Bool_t AliCMEAnalysisCuts::IsSelected(Float_t* values) const {
  //
  // track selection for TPC event plane
  //

  for(Int_t icut=0; icut<Ncuts(); icut++){


    if(Min(icut)>-999.5&&Min(icut)<-998.5){
      Int_t flag=(Int_t) (values[(Int_t) (Type(icut)+Max(icut))]);
      if(!ExcludeRange(icut)&&flag) return kFALSE;
      if(ExcludeRange(icut)&&flag==0) return kFALSE;
    }
    else if(!ExcludeRange(icut)){
      if(values[Type(icut)]<=Min(icut)) return kFALSE;
      if(values[Type(icut)]>=Max(icut)) return kFALSE;
    }
    else {
      if(values[Type(icut)]>Min(icut)&&values[Type(icut)]<Max(icut)) return kFALSE;
    }

  }
  return kTRUE;
}




#endif
