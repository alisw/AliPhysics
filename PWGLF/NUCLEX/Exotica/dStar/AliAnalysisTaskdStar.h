/// \class AliAnalysisdStar
/// \brief This task fills histograms required to perform the analysis on the light nuclei yield.
///
/// The histograms filled in this tasks are used by the analysis macro
/// ## Monte Carlo
/// There are mainly three items studied here:
/// * the acceptance x efficiency;
/// * the difference of the reconstructed \f$p_{T}\f$ with respect to the MC truth \f$p_{T}\f$
///
/// \author Maximiliano Puccio <maximiliano.puccio@unito.it> and Pietro Fecchio <pietro.fecchio@edu.unito.it>, University and INFN Torino
/// \date Apr 14, 2017

#ifndef __AliAnalysisTaskdStar__
#define __AliAnalysisTaskdStar__

#include <Rtypes.h>
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliPID.h"

#include <TMath.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> FourVector_t ;
using std::vector;

class AliVTrack;
class TH2F;
class TList;


// Define bit flags for datacompression
const unsigned char c = 0x1;   // on if charge is +1
const unsigned char p = 0x2;   // on if is Physical Primary
const unsigned char s = 0x4;   // on if is Seconbdary from Material
const unsigned char t = 0x8;   // on if has TOF
const unsigned char i = 0x10;  // on if has ITS standalone PID

struct mother_struct{
  int id;
  bool pi_tof;
  bool pi_tpc;
  bool deuteron_tof;
  int n_daughters;
  FourVector_t vec;
  bool operator==(const int &id_comp) const {return id==id_comp;}
};

struct daughter_struct{
  int mother_pdg;
  int mother_id;
  int mc_truth;     // pdg code of the generated particle
  unsigned char properties;
  FourVector_t vec;
};


class AliAnalysisTaskdStar : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskdStar(const char* taskname = "dStar");
  virtual ~AliAnalysisTaskdStar();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *) {}

  AliEventCuts fEventCut;               /// Standard cuts for event selection
  unsigned int fFilterBit;              /// AOD filter bit for the tracks used in this analysis

private:
  AliAnalysisTaskdStar (const AliAnalysisTaskdStar&);
  AliAnalysisTaskdStar& operator=(const AliAnalysisTaskdStar&) {return *this;}

  bool HasTOF(AliVTrack *track);

  TList                *fList;                      //!<  Output list

  float                fRequireYmin;                ///<  Cut on tracks: mimimum y for the track (using PDG mass)
  float                fRequireYmax;                ///<  Cut on tracks: maximum y for the track (using PDG mass)
  float                fDalitPlotMassCutMin;        ///<  Cut on min d* candidate's invariant mass for Dalitz Plot
  float                fDalitPlotMassCutMax;        ///<  Cut on max d* candidate's invariant mass for Dalitz Plot

  AliPIDResponse       *fPID;                       //!<! PID response class

  // MC only histograms
  TH2F                 *fProduction[2];             //!<! *(MC only)* Total number of produced particles dStar state][Matter-Antimatter]
  TH2F                 *fReconstructed[2][4];       //!<! *(MC only)* Positive and negative tracks reconstructed in the acceptance (ITS-TPC,ITS-TPC-TOF,ITS-TPC-(TOF),ITSsa) [Matter-Antimatter][Detector]
  TH2F                 *fTotal[2];                  //!<! *(MC only)* Positively and negatively charged particles in acceptance : [dStar state][Matter-Antimatter]
  TH2F                 *fMCDalitzPlot;              //!<!

  TTree                *fTree;                      //!<!

  Short_t               fZvtx;                      //<

  vector<daughter_struct>   fDeuteronVector;        //<
  vector<daughter_struct>   fPiPlusVector;          //<
  vector<daughter_struct>   fPiMinusVector;         //<
  vector<daughter_struct>   fMCDeuteronVector;      //<
  vector<daughter_struct>   fMCPiPlusVector;        //<
  vector<daughter_struct>   fMCPiMinusVector;       //<


  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskdStar,1)
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskdStar__) */
