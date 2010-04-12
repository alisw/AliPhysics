#ifndef ALIHFEPIDESD_H
#define ALIHFEPIDESD_H

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

class AliVParticle;

class AliHFEpidESD : public AliHFEpidBase{
  public:
    AliHFEpidESD();
    AliHFEpidESD(const AliHFEpidESD &ref);
    AliHFEpidESD &operator=(const AliHFEpidESD &ref);

    virtual void InitializePID();
    virtual Int_t IsSelected(const AliVParticle *track);
    void SetRequireTOFRange(Double_t pmin, Double_t pmax);
    void SetRequireTRDRange(Double_t pmin, Double_t pmax);
    void SetRequireMinTRDtracklets();
  private:
    Double_t fPminTRD;		// Min. Momentum where TRD PID is required
    Double_t fPmaxTRD;		// Max. Momentum where TRD PID is required
    Double_t fPminTOF;		// Min. Momentum where TOF PID is required
    Double_t fPmaxTOF;		// Max. Momentum where TOF PID is required
    Int_t fMinTrackletsTRD;	// Min. Number of TRD tracklets in region where TRD PID is required
  ClassDef(AliHFEpidESD) // ESD PID class
};
#endif
