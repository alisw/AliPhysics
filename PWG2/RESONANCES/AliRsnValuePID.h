#ifndef ALIRSNVALUEPID_H
#define ALIRSNVALUEPID_H

#include "AliPID.h"
#include "AliRsnValue.h"

class AliVTrack;
class AliPIDResponse;

class AliRsnValuePID : public AliRsnValue {

public:

   enum EValuePID {
      kITSsignal,
      kITSnsigma,
      kTPCsignal,
      kTPCnsigma,
      kTOFsignal,
      kTOFnsigma,
      kTOFtime,
      kTOFsigma,
      kValues
   };

   AliRsnValuePID();
   AliRsnValuePID(const char *name, EValuePID type, AliPID::EParticleType species, Int_t nbins = 0, Double_t min = 0.0, Double_t max = 0.0);
   AliRsnValuePID(const char *name, EValuePID type, AliPID::EParticleType species, Double_t min, Double_t max, Double_t step);
   AliRsnValuePID(const char *name, EValuePID type, AliPID::EParticleType species, Int_t nbins, Double_t *array);
   AliRsnValuePID(const AliRsnValuePID& copy);
   AliRsnValuePID& operator=(const AliRsnValuePID& copy);

   virtual ~AliRsnValuePID() { }
   
   void            SetValuePID(EValuePID type) {fValuePID = type;}
   EValuePID       GetValuePID()               {return fValuePID;}
   
   virtual Bool_t  Eval(TObject *object, Bool_t useMC = kFALSE);
   virtual void    Print(Option_t *option = "") const;

protected:

   void   InitializePID();
   Bool_t TOFComputations(AliVTrack *track);

   AliPID::EParticleType  fSpecies;                    //  particle species
   EValuePID              fValuePID;                   //  output object
   AliPIDResponse        *fPID;                        //  PID response object
   Double_t               fTOFtimes[AliPID::kSPECIES]; //! TOF times
   Double_t               fTOFsigma[AliPID::kSPECIES]; //! TOF sigma

   ClassDef(AliRsnValuePID,1)                          //   AliRsnValuePID class

};

#endif
