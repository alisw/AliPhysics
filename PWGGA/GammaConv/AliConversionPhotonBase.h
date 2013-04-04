#ifndef ALICONVERSIONPHOTONBASE_H
#define ALICONVERSIONPHOTONBASE_H

#include "TMath.h"
#include "TParticle.h"
#include "AliStack.h"
#include "AliLog.h"
#include "TObject.h"
#include "AliMCEvent.h"   
#include "AliESDEvent.h"
#include "AliKFParticle.h"
#include "TParticle.h"
#include <vector>
#include "AliESDpid.h"
#include "TF1.h"
#include "TRandom3.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliKFVertex.h"
#include "AliMCEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "TList.h"


using namespace std;

class AliConversionPhotonBase {

 public: 

  //Constructors
  AliConversionPhotonBase();    

  //Copy Constructor
  AliConversionPhotonBase(const AliConversionPhotonBase & g);           
  //assignment operator
  AliConversionPhotonBase & operator = (const AliConversionPhotonBase & g);

  //Destructor
  virtual ~AliConversionPhotonBase();

  ///Set the tag for decay meson
  void SetTag( Bool_t tagged ) { fTagged = tagged; }
  Bool_t IsTagged(){return fTagged;}

  //Get the Chi2 of particle
  void SetChi2perNDF(Float_t chi2) {fChi2perNDF = chi2;}
  Float_t GetChi2perNDF() const {return fChi2perNDF;}

  ///Track labels
  void SetLabelPositive(Int_t label){fLabel[0] = label;}
  void SetLabelNegative(Int_t label){fLabel[1] = label;}
  void SetTrackLabels(Int_t label1, Int_t label2){fLabel[0] = label1; fLabel[1] = label2;}
  Int_t GetTrackLabelPositive() const{return fLabel[0];}
  Int_t GetTrackLabelNegative() const {return fLabel[1];}
  Int_t GetTrackLabel(Int_t i) const {return fLabel[i];}
  virtual Int_t GetLabel(Int_t i) const { return GetTrackLabel(i); };
  virtual Int_t GetLabel1() const { return GetTrackLabelPositive(); };
  virtual Int_t GetLabel2() const { return GetTrackLabelNegative(); };

  // MC Label

  void SetMCLabel(Int_t* label){fMCLabel[0]=label[0];fMCLabel[1]=label[1];}
  void SetMCLabelPositive(Int_t label){fMCLabel[0]=label;}
  void SetMCLabelNegative(Int_t label){fMCLabel[1]=label;}
  Int_t GetMCLabel(Int_t i) const{return fMCLabel[i];}
  Int_t GetMCLabelPositive() const{return fMCLabel[0];}
  Int_t GetMCLabelNegative() const{return fMCLabel[1];}
  Int_t GetMCParticleLabel(AliStack *fMCStack);

  // GetMCParticle

  Bool_t IsTruePhoton(AliStack *fMCStack);
  TParticle *GetMCParticle(AliStack *fMCStack);
  TParticle *GetPositiveMCDaughter(AliStack *fMCStack){return GetMCDaughter(fMCStack,0);};
  TParticle *GetNegativeMCDaughter(AliStack *fMCStack){return GetMCDaughter(fMCStack,1);};
  TParticle *GetMCDaughter(AliStack *fMCStack,Int_t label);

  // V0Index
  Int_t GetV0Index() const {return fV0Index;}
  void SetV0Index(Int_t index) {fV0Index=index;}

  // Conversion Point
   void SetConversionPoint(Double_t convpoint[3]){fConversionPoint[0]=convpoint[0];fConversionPoint[1]=convpoint[1];fConversionPoint[2]=convpoint[2];}
  void GetConversionPoint(Double_t convpoint[3]){convpoint[0]=fConversionPoint[0];convpoint[1]=fConversionPoint[1];convpoint[2]=fConversionPoint[2];}
  Double_t GetConversionRadius() const {return TMath::Sqrt(fConversionPoint[0]*fConversionPoint[0]+fConversionPoint[1]*fConversionPoint[1]);}
  Double_t GetConversionX() const {return fConversionPoint[0];}
  Double_t GetConversionY() const {return fConversionPoint[1];}
  Double_t GetConversionZ() const {return fConversionPoint[2];}

  // Armenteros Qt Alpha
  void GetArmenterosQtAlpha(Double_t qtalpha[2]){qtalpha[0]=fArmenteros[0];qtalpha[1]=fArmenteros[1];}
  Double_t GetArmenterosQt() const {return fArmenteros[0];}
  Double_t GetArmenterosAlpha() const {return fArmenteros[1];}

  // virtual functions to be implemented in KF/AOD classes

  virtual Double_t GetPhotonMass() const = 0;
  virtual Double_t GetPhotonPt()const = 0;
  virtual Double_t GetPhotonP() const = 0;
  virtual Double_t GetPhotonEta() const = 0;
  virtual Double_t GetPhotonPhi() const =0;
//  virtual Double_t GetPhotonTheta() const =0;
  virtual Double_t GetPx() const = 0;
  virtual Double_t GetPy() const = 0;
  virtual Double_t GetPz() const = 0;


  Float_t GetMass() const { return fIMass; }
  void SetMass( Float_t mass) { fIMass = mass; }

  Float_t GetPsiPair() const {return fPsiPair;}
  void SetPsiPair(Float_t PsiPair){fPsiPair=PsiPair;}

 protected:

  Int_t fLabel[2]; // Electron/Positron Track Label
  Int_t fV0Index; // Index of the V0
  Int_t fMCLabel[2]; // Electron/Positron MC Label
  Float_t fChi2perNDF; // Chi2perNDF
  Double_t fArmenteros[2]; // Armenteros Paramters
  Double_t fConversionPoint[3]; // Conversion Point
  Bool_t fTagged; // Is it tagged as decay pion (only for gammas)
  Float_t fIMass; // Invariant Mass of dilepton pair
  Float_t fPsiPair; // Psi Pair Value

  ClassDef(AliConversionPhotonBase,3);
};


#endif



