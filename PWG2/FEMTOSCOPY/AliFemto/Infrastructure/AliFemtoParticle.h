/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   Particle objects are part of the PicoEvent, which is what is
 *   stored in the EventMixingBuffers
 *   A Track object gets converted to a Particle object if it
 *   passes the ParticleCut of an Analysis
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.2  2007-04-03 16:00:08  mchojnacki
 * Changes to iprove memory managing
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.19  2003/01/14 09:41:16  renault
 * changes on average separation calculation, hit shared finder and memory optimisation
 * for Z,U and Sectors variables.
 *
 * Revision 1.18  2002/12/12 17:01:50  kisiel
 * Hidden Information handling and purity calculation
 *
 * Revision 1.17  2002/11/19 23:35:52  renault
 * Enable calculation of exit/entrance separation for V0 daughters
 *
 * Revision 1.16  2001/12/14 23:11:30  fretiere
 * Add class HitMergingCut. Add class fabricesPairCut = HitMerginCut + pair purity cuts. Add TpcLocalTransform function which convert to local tpc coord (not pretty). Modify AliFemtoTrack, AliFemtoParticle, AliFemtoHiddenInfo, AliFemtoPair to handle the hit information and cope with my code
 *
 * Revision 1.15  2001/05/25 23:23:59  lisa
 * Added in AliFemtoKink stuff
 *
 * Revision 1.14  2001/05/23 00:19:05  lisa
 * Add in Smearing classes and methods needed for momentum resolution studies and correction
 *
 * Revision 1.13  2001/04/03 21:04:36  kisiel
 *
 *
 *   Changes needed to make the Theoretical code
 *   work. The main code is the ThCorrFctn directory.
 *   The most visible change is the addition of the
 *   HiddenInfo to AliFemtoPair.
 *
 * Revision 1.12  2000/10/05 23:09:05  lisa
 * Added kT-dependent radii to mixed-event simulator AND implemented AverageSeparation Cut and CorrFctn
 *
 * Revision 1.11  2000/07/17 20:03:17  lisa
 * Implemented tools for addressing and assessing trackmerging
 *
 * Revision 1.10  2000/07/16 21:38:23  laue
 * AliFemtoCoulomb.cxx AliFemtoSectoredAnalysis.cxx : updated for standalone version
 * AliFemtoV0.cc AliFemtoV0.h : some cast to prevent compiling warnings
 * AliFemtoParticle.cc AliFemtoParticle.h : pointers fTrack,fV0 initialized to 0
 * AliFemtoIOBinary.cc : some printouts in #ifdef STHBTDEBUG
 * AliFemtoEvent.cc : B-Field set to 0.25Tesla, we have to think about a better
 *                 solution
 *
 * Revision 1.9  2000/05/03 17:44:43  laue
 * AliFemtoEvent, AliFemtoTrack & AliFemtoV0 declared friend to AliFemtoIOBinary
 * AliFemtoParticle updated for V0 pos,neg track Id
 *
 * Revision 1.8  2000/04/03 16:21:51  laue
 * some include files changed
 * Multi track cut added
 *
 * Revision 1.6  2000/02/26 19:04:52  laue
 * Some unnecessary includes removed.
 * StThreeVectorD replace by AliFemtoThreeVector.
 * AliFemtoCoulomb modified to compile without Root (ClassDef embraced into
 *   #ifdef __ROOT__  ..... #endif)
 * AliFemtoParticle now returns references (FourMomentum(),Helix(),
 *   DecayVertexPosiion())
 *
 * Revision 1.5  1999/12/11 15:58:29  lisa
 * Add vertex decay position datum and accessor to AliFemtoParticle to allow pairwise cuts on seperation of V0s
 *
 * Revision 1.4  1999/09/17 22:38:02  lisa
 * first full integration of V0s into AliFemto framework
 *
 * Revision 1.3  1999/09/01 19:04:54  lisa
 * update Particle class AND add parity cf and Randys Coulomb correction
 *
 * Revision 1.2  1999/07/06 22:33:23  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#ifndef AliFemtoParticle_hh
#define AliFemtoParticle_hh

//#include "math.h"
#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoTrack.h"
#include "Infrastructure/AliFemtoV0.h"
#include "Infrastructure/AliFemtoKink.h"
#include "Infrastructure/AliFemtoXi.h"
#include "AliFmPhysicalHelixD.h"
// ***
class AliFemtoHiddenInfo;
// ***
class AliFemtoParticle{
public:
  AliFemtoParticle();
  AliFemtoParticle(const AliFemtoTrack* const hbtTrack, const double& mass);
  AliFemtoParticle(const AliFemtoV0* const hbtV0, const double& mass);
  AliFemtoParticle(const AliFemtoKink* const hbtKink, const double& mass);
  AliFemtoParticle(const AliFemtoXi* const hbtXi, const double& mass);
  ~AliFemtoParticle();

  const AliFemtoLorentzVector& FourMomentum() const;

  AliFmPhysicalHelixD& Helix();

  const AliFemtoThreeVector DecayVertexPosition() const;
  unsigned long TopologyMap(const int word) const;
  int NumberOfHits() const;

  unsigned long TrackId() const;         // only for particles from tracks 
  unsigned short   NegTrackId() const;   // only for particles from v0 
  unsigned short   PosTrackId() const;   // only for particles from v0 

  AliFemtoTrack* Track() const;
  AliFemtoV0* V0() const;
  AliFemtoKink* Kink() const;

  const AliFemtoThreeVector& NominalTpcExitPoint() const;     // position track exits TPC assuming start at (0,0,0)
  const AliFemtoThreeVector& NominalTpcEntrancePoint() const; // position track crosses IFC assuming start at (0,0,0)
  const AliFemtoThreeVector& TpcV0PosExitPoint() const;  
  const AliFemtoThreeVector& TpcV0PosEntrancePoint() const;
  const AliFemtoThreeVector& TpcV0NegExitPoint() const;  
  const AliFemtoThreeVector& TpcV0NegEntrancePoint() const;

  // the following method is for explicit internal calculation to fill datamembers.
  // It is invoked automatically if AliFemtoParticle constructed from AliFemtoTrack
  //void CalculateNominalTpcExitAndEntrancePoints(); // NOTE - this requires the fHelix, so be sure this is filled


  AliFemtoThreeVector fNominalPosSample[11];  // I make this public for convenience and speed of AliFemtoPair()
  float fZ[45];
  float fU[45];
  int fSect[45];

  void ResetFourMomentum(const AliFemtoLorentzVector& fourMomentum);

  const AliFemtoHiddenInfo*  HiddenInfo() const;
  // Fab private
  AliFemtoHiddenInfo*  getHiddenInfo() const;
  void SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo);
  void CalculatePurity();
  double GetPionPurity();
  double GetKaonPurity();
  double GetProtonPurity();
  void CalculateTpcExitAndEntrancePoints( AliFmPhysicalHelixD* tHelix,
					  AliFemtoThreeVector* PrimVert,
					  AliFemtoThreeVector* SecVert,
					  AliFemtoThreeVector* tmpTpcEntrancePoint,
					  AliFemtoThreeVector* tmpTpcExitPoint,
					  AliFemtoThreeVector* tmpPosSample,
					  float* tmpZ,float* tmpU,int* tmpSect);

  // For V0 Neg Daugthers TpcEntrance/ExitPoints
  AliFemtoThreeVector* fTpcV0NegPosSample;
  float* fV0NegZ;
  float* fV0NegU;
  int* fV0NegSect;
 
private:
  AliFemtoTrack* fTrack;  // copy of the track the particle was formed of, else Null
  AliFemtoV0* fV0;        // copy of the v0 the particle was formed of, else Null
  AliFemtoKink* fKink;        // copy of the v0 the particle was formed of, else Null
  AliFemtoXi* fXi;

  AliFemtoLorentzVector fFourMomentum;
  AliFmPhysicalHelixD fHelix;
  //unsigned long  fMap[2]; 
  //int fNhits;
  AliFemtoThreeVector fNominalTpcExitPoint;
  AliFemtoThreeVector fNominalTpcEntrancePoint;
  AliFemtoHiddenInfo* fHiddenInfo;  // Fab private

  double fPurity[6];

  static double fPrimPimPar0;
  static double fPrimPimPar1;
  static double fPrimPimPar2;
  static double fPrimPipPar0;
  static double fPrimPipPar1;
  static double fPrimPipPar2;
  static double fPrimPmPar0;
  static double fPrimPmPar1;
  static double fPrimPmPar2;
  static double fPrimPpPar0;
  static double fPrimPpPar1;
  static double fPrimPpPar2;

   // For V0 Daugthers TpcEntrance/ExitPoints
  AliFemtoThreeVector fPrimaryVertex;
  AliFemtoThreeVector fSecondaryVertex;

  AliFmPhysicalHelixD fHelixV0Pos;
  AliFemtoThreeVector fTpcV0PosEntrancePoint;
  AliFemtoThreeVector fTpcV0PosExitPoint;

  AliFmPhysicalHelixD fHelixV0Neg;
  AliFemtoThreeVector fTpcV0NegEntrancePoint;
  AliFemtoThreeVector fTpcV0NegExitPoint;
};

inline AliFemtoTrack* AliFemtoParticle::Track() const { return fTrack; }
inline unsigned long  AliFemtoParticle::TrackId() const { return fTrack->TrackId(); }; 
inline const AliFemtoLorentzVector& AliFemtoParticle::FourMomentum() const {return fFourMomentum;}
inline AliFmPhysicalHelixD& AliFemtoParticle::Helix() {return fHelix;}
//inline unsigned long AliFemtoParticle::TopologyMap(const int word) const {return fMap[word];}
//inline int AliFemtoParticle::NumberOfHits() const {return fNhits;}
//by marek chojnacki to could compile 

inline unsigned long AliFemtoParticle::TopologyMap(const int word) const {return 1;}
inline int AliFemtoParticle::NumberOfHits() const {return 1;}

inline AliFemtoV0* AliFemtoParticle::V0() const { return fV0; }
inline unsigned short AliFemtoParticle::NegTrackId() const { return fV0->idNeg(); }
inline unsigned short AliFemtoParticle::PosTrackId() const { return fV0->idPos(); }
inline const AliFemtoThreeVector AliFemtoParticle::DecayVertexPosition() const {return fV0->decayVertexV0(); }
// ***
inline AliFemtoHiddenInfo* AliFemtoParticle::getHiddenInfo() const
{return fHiddenInfo;}
inline const AliFemtoHiddenInfo* AliFemtoParticle::HiddenInfo() const
{return fHiddenInfo;}
inline void AliFemtoParticle::SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo)
{ fHiddenInfo = aHiddenInfo->clone();}
// ***

inline void AliFemtoParticle::ResetFourMomentum(const AliFemtoLorentzVector& vec){fFourMomentum = vec;}

inline AliFemtoKink* AliFemtoParticle::Kink() const { return fKink; }

#endif
