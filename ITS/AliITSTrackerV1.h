#ifndef ALIITSTRACKERV1_H
#define ALIITSTRACKERV1_H


class AliITS;
class AliITSRad;
 
//   ITS Tracker V1 Class
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it

class AliITSTrackerV1 : public TObject {

  public:
    
	 AliITSTrackerV1(AliITS* IITTSS);

    AliITStrack Tracking(AliITStrack &track, AliITStrack *reference, TObjArray *fpoints, Int_t **vettid,
	 Bool_t flagvert,  AliITSRad *rl, AliITSgeoinfo *geoinfo);  

    void DoTracking(Int_t evNumber, Int_t min_t, Int_t max_t, TFile *file, Bool_t flagvert);
	 
  private:

    AliITS* ITS;
	 

    ClassDef(AliITSTrackerV1,1)
};

#endif
