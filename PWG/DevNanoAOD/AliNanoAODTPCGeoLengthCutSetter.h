#ifndef AliNanoAODTPCGeoLengthCutSetter_h
#define AliNanoAODTPCGeoLengthCutSetter_h

#include "AliNanoAODCustomSetter.h"
#include "AliNanoAODHeader.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"


class AliNanoAODTPCGeoLengthCutSetter : public AliNanoAODCustomSetter
{
public:
  AliNanoAODTPCGeoLengthCutSetter(const char * name = "AliNanoAODTPCGeoLengthCutSetter") : 
    AliNanoAODCustomSetter(name), fMode(0), fDeltaY(3.0), fDeltaZ(220.0), fMagField(-1),
    fRequireCutGeoNcrNclLength(130), fRequireCutGeoNcrNclGeom1Pt(1.5),
    fCutGeoNcrNclFractionNcr(0.85), fCutGeoNcrNclFractionNcl(0.7), fIndex(-1) {;}
  virtual ~AliNanoAODTPCGeoLengthCutSetter() {;}
  inline virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head , TString varListHeader  ) {
      fMagField = event->GetMagneticField();
      return;
  }
  inline virtual void SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * spTrack) {
      if (fIndex == -1)
          fIndex = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstTPCGeoLength");
      if ((fIndex >= 0) && (fMagField != -1) ){
        auto checkResult = true;
        AliESDtrack fESDTrack(aodTrack);
        fESDTrack.SetTPCClusterMap(aodTrack->GetTPCClusterMap());
        fESDTrack.SetTPCSharedMap(aodTrack->GetTPCSharedMap());
        fESDTrack.SetTPCPointsF(aodTrack->GetTPCNclsF());

        auto nCrossedRowsTPC = fESDTrack.GetTPCCrossedRows();
        auto lengthInActiveZoneTPC=fESDTrack.GetLengthInActiveZone(fMode,fDeltaY,fDeltaZ,fMagField);
        auto cutGeoNcrNclLength=fRequireCutGeoNcrNclLength-TMath::Power(TMath::Abs(fESDTrack.GetSigned1Pt()),fRequireCutGeoNcrNclGeom1Pt);

        if (lengthInActiveZoneTPC < cutGeoNcrNclLength) 
            checkResult = false;
        if (nCrossedRowsTPC<fCutGeoNcrNclFractionNcr*cutGeoNcrNclLength) 
            checkResult = false;
        if (fESDTrack.GetTPCncls()<fCutGeoNcrNclFractionNcl*cutGeoNcrNclLength) 
            checkResult = false;

        spTrack->SetVar(fIndex, (checkResult) ? 1. : 0.);
      }
  };

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
