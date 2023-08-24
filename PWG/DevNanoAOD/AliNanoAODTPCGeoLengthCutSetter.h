#ifndef AliNanoAODTPCGeoLengthCutSetter_h
#define AliNanoAODTPCGeoLengthCutSetter_h

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliNanoAODCustomSetter.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"

class AliNanoAODTPCGeoLengthCutSetter : public AliNanoAODCustomSetter {
 public:
  AliNanoAODTPCGeoLengthCutSetter(
      const char* name = "AliNanoAODTPCGeoLengthCutSetter");
  virtual ~AliNanoAODTPCGeoLengthCutSetter();
  virtual void SetNanoAODHeader(const AliAODEvent* event,
                                AliNanoAODHeader* head,
                                TString varListHeader);
  virtual void SetNanoAODTrack(const AliAODTrack* aodTrack,
                               AliNanoAODTrack* spTrack);

  int fMode;
  Double_t fDeltaY;
  Double_t fDeltaZ;
  Double_t fMagField;
  Double_t fRequireCutGeoNcrNclLength;
  Double_t fRequireCutGeoNcrNclGeom1Pt;
  Double_t fCutGeoNcrNclFractionNcr;
  Double_t fCutGeoNcrNclFractionNcl;

 protected:
  bool fGoodToGo;
  int fIndex;

  ClassDef(AliNanoAODTPCGeoLengthCutSetter, 1)
};

#endif /* AliNanoAODTPCGeoLengthCutSetter_h */
