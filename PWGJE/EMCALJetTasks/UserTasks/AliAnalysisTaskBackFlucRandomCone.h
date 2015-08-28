/// \class AliAnalysisTaskBackFlucRandomCone
/// \brief Study background fluctuations with random cones
///
///
/// \author Chiara Bianchin
/// \date Aug, 2015

#ifndef ALIANALYSISTASKBACKFLUCRANDOMCONE_H
#define ALIANALYSISTASKBACKFLUCRANDOMCONE_H

class TH1F;
class TH2F;
class AliRhoParameter;
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskBackFlucRandomCone : public AliAnalysisTaskEmcalJet {

public:
   AliAnalysisTaskBackFlucRandomCone();
   AliAnalysisTaskBackFlucRandomCone(const char* name);
   virtual ~AliAnalysisTaskBackFlucRandomCone();

   void        UserCreateOutputObjects();
   void        Terminate(Option_t *option);
   
   //setters
   void        SetEtaRange(Double_t min, Double_t max) {fEtaMin = min; fEtaMax = max;}
   void        SetConeRadius(Double_t R) {fRCone = R;}
   void        SetRhoName();
   
protected:
   Bool_t      Run();
   Bool_t      FillHistograms();

private:
   Double_t    fEtaMin;                       ///< Minimum eta limit (consider the R!), default -0.5
   Double_t    fEtaMax;                       ///< Maximum eta limit (consider the R!), default 0.5
   Double_t    fEta;                          ///!<! Eta coordinate of the cone centre
   Double_t    fPhi;                          ///!<! Phi coordinate of the cone centre
   Double_t    fRCone;                        ///< Radius of the cone, default 0.3
   
   TRandom3    fRnd;                          ///< Random number generator
   TString     fRhoName;                      ///< Name of rho run previously in the event
   
   TH2F       *fhEtaPhiConeCentre;            ///!<! Cone centre
   TH1F       *fNInOut;                       ///!<! Count track in and out of cone per event
   TH2F       *fhMasspTInCone;                ///!<! Mass and pT of particles in the cone
   TH1F       *fhNConstituents;               ///!<! Number of constituents in random cone
   TH1F       *fhDeltapT;                     ///!<! pT cone minus rho
   
   AliAnalysisTaskBackFlucRandomCone(const AliAnalysisTaskBackFlucRandomCone&);
   AliAnalysisTaskBackFlucRandomCone &operator=(const AliAnalysisTaskBackFlucRandomCone&);
   
   ClassDef(AliAnalysisTaskBackFlucRandomCone, 1)
};
#endif


