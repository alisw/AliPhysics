/// \class AliAnalysisTaskPsEfficiency
/// \brief This task fills histograms required to perform the analysis on the light nuclei yield.
///
/// The histograms filled in this tasks are used by the analysis macro
/// ## Monte Carlo
/// There are mainly three items studied here:
/// * the acceptance x efficiency;
/// * the difference of the reconstructed \f$p_{T}\f$ with respect to the MC truth \f$p_{T}\f$
///
/// \author Maximiliano Puccio <maximiliano.puccio@unito.it> and Luca Barioglio <luca.barioglio@unito.it>, University and INFN Torino
/// \date Apr 14, 2017

#ifndef __AliAnalysisTaskPsEfficiency__
#define __AliAnalysisTaskPsEfficiency__

#include <Rtypes.h>
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliPID.h"

class AliVTrack;
class TH2F;
class TList;

class AliAnalysisTaskPsEfficiency : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPsEfficiency(const char* taskname = "PsEfficiency");
  virtual ~AliAnalysisTaskPsEfficiency();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *) {}

  AliEventCuts fEventCut;               /// Standard cuts for event selection
  unsigned int fFilterBit;              /// AOD filter bit for the tracks used in this analysis

private:
  AliAnalysisTaskPsEfficiency (const AliAnalysisTaskPsEfficiency&);
  AliAnalysisTaskPsEfficiency& operator=(const AliAnalysisTaskPsEfficiency&) {return *this;}

  bool HasTOF(AliVTrack *track);

  TList                *fList;                      ///<  Output list

  float                fRequireYmin;                ///<  Cut on tracks: mimimum y for the track (using PDG mass)
  float                fRequireYmax;                ///<  Cut on tracks: maximum y for the track (using PDG mass)

  AliPIDResponse       *fPID;                       //!<! PID response class

  // MC only histograms
  TH2F                 *fProduction[2][2];             //!<! *(MC only)* Total number of produced particles [Ps state][Matter-Antimatter]
  TH2F                 *fReconstructed[2][2][3];    //!<! *(MC only)* Positive and negative tracks reconstructed in the acceptance (ITS-TPC,ITS-TPC-TOF,ITS-TPC-(TOF)) [Ps state][Matter-Antimatter][Detector]
  TH2F                 *fTotal[2][2];               //!<! *(MC only)* Positively and negatively charged particles in acceptance : [Ps state][Matter-Antimatter]

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskPsEfficiency,1)
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskPsEfficiency__) */
