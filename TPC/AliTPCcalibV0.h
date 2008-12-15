#ifndef AliTPCCALIBV0_H
#define AliTPCCALIBV0_H


#include <AliTPCcalibBase.h>


class TTreeSRedirector;
class AliTPCROC;
class AliTPCseed;
class AliESDtrack;
class AliESDEvent;
class TH3F;
class TH1F;
class TH2F;
class TH1I;
class TDatabasePDG;
class AliKFParticle;
class AliKFVertex;
class AliESDv0;
class TArrayI;
class TTree;
class AliStack;

class AliTPCcalibV0 : public AliTPCcalibBase {
public :

   // List of branches

  AliTPCcalibV0();
  virtual ~AliTPCcalibV0();
  virtual void     Process(AliESDEvent *event) {return ProcessESD(event,0);}

  //
  //
  //
  void ProcessESD(AliESDEvent *esd, AliStack *stack=0);
  void MakeMC();
  void MakeV0s();
  void ProcessV0(Int_t ftype);
  void ProcessPI0();
  TH2F * GetHistograms();
  void BinLogX(TH2F * h);
  //
  //
  //  
  static AliKFParticle * Fit(AliKFVertex & primVtx, AliESDv0 *v0, Int_t PDG1, Int_t PDG2);
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}

protected:
private:

   AliTPCcalibV0(const AliTPCcalibV0&); // Not implemented
   AliTPCcalibV0& operator=(const AliTPCcalibV0&); // Not implemented


   AliStack       *fStack;        // pointer to kinematic tree        
   AliESDEvent    *fESD;              //! current ED to proccess - NOT OWNER
   TDatabasePDG   *fPdg;              // particle database
   TObjArray      *fParticles;         // array of selected MC particles
   TObjArray      *fV0s;               // array of V0s
   TObjArray      *fGammas;           // gamma conversion candidates
   //
   TArrayI        *fV0type;            // array of types for V0s       
   TH2F           *fTPCdEdx;              // dEdx spectra
   TH2F           *fTPCdEdxPi;            // dEdx spectra - pion anti-pion
   TH2F           *fTPCdEdxEl;            // dEdx spectra - electroen -positrons from gamma
   TH2F           *fTPCdEdxP;             // dEdx spectra - proton antiproton - lambda -  antilambda
   //       
   ClassDef(AliTPCcalibV0,1);
};


#endif
