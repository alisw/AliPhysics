#ifndef ALIRSNCUTPHIRAA_H
#define ALIRSNCUTPHIRAA_H

//
// This cut implements all the checks done to accept a track as a Kaon
// for the PbPb analysis using 2010 runs.
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"
#include "AliESDtrackCuts.h"

class AliRsnCutPhiRAA : public AliRsnCut {

public:

   enum ECutMode {
      k2010 = 0,
      k2011_0,
      k2011_1,
      k2011_1_05,
      k2011_1_075
   };

   AliRsnCutPhiRAA(const char *name = "");
   AliRsnCutPhiRAA(const AliRsnCutPhiRAA &copy);
   AliRsnCutPhiRAA &operator=(const AliRsnCutPhiRAA &copy);
   virtual ~AliRsnCutPhiRAA() { }

   virtual Bool_t IsSelected(TObject *obj);

   AliRsnCutTrackQuality *CutQuality()            {return &fCutQuality;}

   void   SetMode(ECutMode mode)            {fMode = mode;}

private:

   ECutMode              fMode;          // how the cut is applied
   AliRsnCutTrackQuality fCutQuality;    // track quality cut
   AliESDtrackCuts *cut1;
   AliESDtrackCuts *cut2;
   AliESDtrackCuts *cut3;


   ClassDef(AliRsnCutPhiRAA,1)

};


#endif
