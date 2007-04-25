/***********************************************************************
 *
 * $Id: AliFemtoV0.h,v 1.0 1999/09/07
 *
 * Authors: Helen Caines, Tom Humanic 07-Sep-1999
 *
 ***********************************************************************
 *
 * Description: V0 class with members copied from StV0hbtV0.h
 *
 ***********************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.13  2003/09/02 17:58:32  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.12  2003/01/14 09:45:14  renault
 * enable this class to be used for theorical calculations.
 *
 * Revision 1.11  2002/11/19 23:37:56  renault
 * Get V0 daughters helix
 *
 * Revision 1.10  2002/02/09 19:25:36  laue
 * updates (dedx length)
 *
 * Revision 1.9  2001/09/05 20:41:43  laue
 * Updates of the hbtMuDstTree microDSTs
 *
 * Revision 1.8  2000/10/09 21:54:23  laue
 * Helens changes to the V0s
 *
 * Revision 1.7  2000/08/09 14:50:22  laue
 * 'const' removed to compile on solaris
 *
 * Revision 1.6  2000/07/16 21:38:23  laue
 * AliFemtoCoulomb.cxx AliFemtoSectoredAnalysis.cxx : updated for standalone version
 * AliFemtoV0.cc AliFemtoV0.h : some cast to prevent compiling warnings
 * AliFemtoParticle.cc AliFemtoParticle.h : pointers mTrack,mV0 initialized to 0
 * AliFemtoIOBinary.cc : some printouts in #ifdef STHBTDEBUG
 * AliFemtoEvent.cc : B-Field set to 0.25Tesla, we have to think about a better
 *                 solution
 *
 * Revision 1.5  2000/05/03 17:44:43  laue
 * AliFemtoEvent, AliFemtoTrack & AliFemtoV0 declared friend to AliFemtoIOBinary
 * AliFemtoParticle updated for V0 pos,neg track Id
 *
 * Revision 1.4  2000/02/18 21:32:24  laue
 * franksTrackCut changed. If mCharge is set to '0' there will be no cut
 * on charge. This is important for front-loaded cuts.
 *
 * copy constructor implemented for AliFemtoEvent, AliFemtoTrack and AliFemtoV0.
 *
 * franks1HistoD.cxx franks1HistoD.h franks2HistoD.cxx franks2HistoD.h
 * removed. We can now (CC5 on Solaris) use the versions (no D)
 *
 * Revision 1.3  2000/02/01 00:33:32  laue
 * namespaces changed to run on the new Solaris Compiler CC5
 * since we can use member templates in franks1Histo.h we are doing it
 *
 * Revision 1.2  1999/10/15 02:32:37  lisa
 * Helen has added method AliFemtoV0::UpdateV0() to fill derived information - this leads to smaller microDSTs
 *
 * Revision 1.1  1999/09/16 18:47:59  lisa
 * replace placeholder HbtV0Track stuff with Helens AliFemtoV0 classes
 *
 * 
 *
 ***********************************************************************/
#ifndef AliFemtoV0_hh
#define AliFemtoV0_hh

#include "Infrastructure/AliFemtoTypes.h" //same as in AliFemtoTrack.h
#include "AliFmPhysicalHelixD.h" // Gael 12 Sept 02
#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#include "StStrangeMuDstMaker/StV0MuDst.h"
#endif
#endif
/* Th stuff */
#include "Base/AliFemtoHiddenInfo.h"
/***/

class AliFemtoV0 {
public:
  AliFemtoV0(){/* no-op */}
  AliFemtoV0( const AliFemtoV0&); // copy constructor
#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
  AliFemtoV0( StV0MuDst&); // from strangeness V0 micro dst structure
#endif
#endif
  ~AliFemtoV0(){if(fHiddenInfo) delete fHiddenInfo;}


  float decayLengthV0() const;       // 3-d decay distance
  AliFemtoThreeVector decayVertexV0() const; // Coordinates of decay vertex
  AliFemtoThreeVector primaryVertex() const; // Coordinates of primary vertex
  float decayVertexV0X() const; // Coordinates of decay vertex
  float decayVertexV0Y() const; // Coordinates of decay vertex
  float decayVertexV0Z() const; // Coordinates of decay vertex
  float dcaV0Daughters() const;      // DCA of v0 daughters at decay vertex
  float dcaV0ToPrimVertex() const;   // DCA of v0 to primary vertex
  float dcaPosToPrimVertex() const;  // DCA of pos v0 daughter to pri vertex
  float dcaNegToPrimVertex() const;  // DCA of neg v0 daughter to pri vertex
  AliFemtoThreeVector momPos() const;   // Momentum components of pos. daughter
  float momPosX() const;   // Momentum components of pos. daughter
  float momPosY() const;   // Momentum components of pos. daughter
  float momPosZ() const;   // Momentum components of pos. daughter
  AliFemtoThreeVector momNeg() const;   // Momentum components of neg. daughter
  float momNegX() const;   // Momentum components of neg. daughter
  float momNegY() const;   // Momentum components of neg. daughter
  float momNegZ() const;   // Momentum components of neg. daughter

  int   tpcHitsPos() const;          // Number of TPC hits on pos. daughter
  int   tpcHitsNeg() const;          // Number of TPC hits on neg. daughter
  unsigned long trackTopologyMapPos(unsigned int) const;
  unsigned long trackTopologyMapNeg(unsigned int) const;

  AliFemtoThreeVector momV0() const ;    // Momentum components of V0
  float momV0X() const ;    // Momentum components of V0
  float momV0Y() const ;    // Momentum components of V0
  float momV0Z() const ;    // Momentum components of V0
  float alphaV0() const ;             // Armenteros-Podolanski variable
  float ptArmV0() const ;            // Armenteros-Podolanski variable
  float eLambda() const ;             // Energy assuming lambda hypothesis
  float eK0Short() const ;            // Energy assuming k-short hypothesis
  float ePosProton() const ;          // Energy of pos. daughter assuming proton
  float ePosPion() const ;            // Energy of pos. daughter assuming pion
  float eNegProton() const ;          // Energy of neg. daughter assuming antiproton
  float eNegPion() const ;            // Energy of neg. daughter assuming pion
  float massLambda() const ;          // Mass assuming lambda hypothesis
  float massAntiLambda() const ;      // Mass assuming antilambda hypothesis
  float massK0Short() const ;         // Mass assuming k-short hypothesis
  float rapLambda() const ;           // Rapidity assuming (anti) constlambda
  float rapK0Short() const ;          // Rapidity assuming k-short
  float cTauLambda() const ;          // Lifetime (ctau) const assuming (anti) constlambda
  float cTauK0Short() const ;         // Lifetime (ctau) const assuming k-short
  float ptV0() const ;                // Transverse momentum
  float ptotV0() const ;              // Total momentum
  float ptPos() const ;               // Transverse momentum of pos. daughter
  float ptotPos() const ;             // Total momentum of pos. daughter
  float dedxPos() const;              // dedx of Positive track
  float numdedxPos() const;
// number of hits in dE/dX track of pos. daughter--Gael04Fev 2002
  float errdedxPos() const;              
// error on dedx of Positive track--Gael04Fev 2002
  float lendedxPos() const;              
// Length of dE/dX track of pos. daughter--Gael04Fev 2002
  float pseudoRapPos() const;              
// Length of dE/dX track of neg. daughter--Gael04Fev 2002

  float ptNeg() const ;               // Transverse momentum of neg. daughter
  float ptotNeg() const ;             // Total momentum of neg. daughter
  float dedxNeg() const;              // dedx of Negative track
  float numdedxNeg() const;
// number of hits in dE/dX track of neg. daughter--Gael04Fev 2002
  float errdedxNeg() const;              
// error on dedx of Negative track--Gael04Fev 2002
  float lendedxNeg() const;              
// Length of dE/dX track of neg. daughter--Gael04Fev 2002
  float pseudoRapNeg() const;              
// Length of dE/dX track of neg. daughter--Gael04Fev 2002

  unsigned short     idNeg() const;               // Id of negative track
  unsigned short   idPos() const;               // Id of positive track
  unsigned short   keyNeg() const;               // Id of negative track
  unsigned short   keyPos() const;               // Id of positive track

  const AliFmPhysicalHelixD& HelixPos() const; // Gael 12 Sept 02
  const AliFmPhysicalHelixD& HelixNeg() const; // Gael 12 Sept 02

  void UpdateV0(); // Fills derived info
  void SetdecayLengthV0(const float);  
  void SetdecayVertexV0(const AliFemtoThreeVector);  
  void SetdecayVertexV0X(const float);
  void SetdecayVertexV0Y(const float);
  void SetdecayVertexV0Z(const float);
  void SetdcaV0Daughters(const float); 
  void SetdcaV0ToPrimVertex(const float);  
  void SetdcaPosToPrimVertex(const float); 
  void SetdcaNegToPrimVertex(const float); 
  void SetmomPos(const AliFemtoThreeVector);  
  void SetmomPosX(const float);  
  void SetmomPosY(const float);  
  void SetmomPosZ(const float);  
  void SetmomNeg(const AliFemtoThreeVector);  
  void SetmomNegX(const float);  
  void SetmomNegY(const float);  
  void SetmomNegZ(const float);  

  void SettpcHitsPos(const int&);      
  void SettpcHitsNeg(const int&);      

  void SetTrackTopologyMapPos(unsigned int, const unsigned long&);
  void SetTrackTopologyMapNeg(unsigned int, const unsigned long&);      

  void SetmomV0( AliFemtoThreeVector);
  void SetmomV0X( float);
  void SetmomV0Y( float);
  void SetmomV0Z( float);
  void SetalphaV0( float);       
  void SetptArmV0( float);       
  void SeteLambda( float);     
  void SeteK0Short( float);    
  void SetePosProton( float);  
  void SetePosPion( float);    
  void SeteNegProton( float);  
  void SeteNegPion( float);    
  void SetmassLambda( float);  
  void SetmassAntiLambda( float);
  void SetmassK0Short( float);  
  void SetrapLambda( float);    
  void SetrapK0Short( float);   
  void SetcTauLambda( float);   
  void SetcTauK0Short( float);  
  void SetptV0( float);         
  void SetptotV0( float);       
  void SetptPos( float);        
  void SetptotPos( float);      
  void SetptNeg( float);        
  void SetptotNeg( float);
  void SetidNeg(const unsigned short&);
  void SetidPos(const unsigned short&);
  void SetdedxNeg(float);
  void SeterrdedxNeg(float x);//Gael 04Fev2002
  void SetlendedxNeg(float x);//Gael 04Fev2002
  void SetpseudoRapNeg(float x);//Gael 04Fev2002
  void SetdedxPos(float);
  void SeterrdedxPos(float x);//Gael 04Fev2002
  void SetlendedxPos(float x);//Gael 04Fev2002
  void SetpseudoRapPos(float x);//Gael 04Fev2002
  void SetkeyNeg(const unsigned short&);
  void SetkeyPos(const unsigned short&);
     
  void SetHelixPos(const AliFmPhysicalHelixD&); // Gael 12 Sept 02
  void SetHelixNeg(const AliFmPhysicalHelixD&); // Gael 12 Sept 02

  void SetprimaryVertex(const AliFemtoThreeVector v);//Gael 24 Sept 02
  /* Th stuff */
  void SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo);
  bool ValidHiddenInfo() const;
  // Fab private : (official : const AliFemtoHiddenInfo* HiddenInfo() const;
  AliFemtoHiddenInfo* getHiddenInfo() const;
  /***/

protected:
 
  float fDecayLengthV0;
  AliFemtoThreeVector fDecayVertexV0;
  AliFemtoThreeVector fPrimaryVertex;
  float fDcaV0Daughters;
  float fDcaV0ToPrimVertex;
  float fDcaPosToPrimVertex;
  float fDcaNegToPrimVertex;
  AliFemtoThreeVector fMomPos;
  AliFemtoThreeVector fMomNeg;

  unsigned long  fTrackTopologyMapPos[2];
  unsigned long  fTrackTopologyMapNeg[2];
       
  int   fTpcHitsPos;
  int   fTpcHitsNeg;

  float fChi2V0;
  float fClV0;
  float fChi2Pos;
  float fClPos;
  float fChi2Neg;
  float fClNeg;

  float fDedxPos;
  float fErrDedxPos;
  float fLenDedxPos;
  
  float fDedxNeg;
  float fErrDedxNeg;
  float fLenDedxNeg;

  unsigned short fNufDedxPos;
  unsigned short fNufDedxNeg;

  AliFmPhysicalHelixD fHelixPos; // Gael 12 Sept 02
  AliFmPhysicalHelixD fHelixNeg; // Gael 12 Sept 02

  // the following variables are not in the persistent version and can be calculated via UpdateV0();
  AliFemtoThreeVector fMomV0;
  float fAlphaV0;
  float fPtArmV0;
  float fELambda;
  float fEK0Short;
  float fEPosProton;
  float fEPosPion;
  float fENegProton;
  float fENegPion;
  float fMassLambda;
  float fMassAntiLambda;
  float fMassK0Short;
  float fRapLambda;
  float fRapK0Short;
  float fCTauLambda;
  float fCTauK0Short;
  float fPtV0;
  float fPtotV0;
  float fPtPos;
  float fPtotPos;
  float fPtNeg;
  float fPtotNeg;

  unsigned short   fKeyNeg;
  unsigned short   fKeyPos;
  /* Th stuff */
  // Fab private : add mutable
  mutable AliFemtoHiddenInfo* fHiddenInfo; //!
  /***/


};



inline float AliFemtoV0::decayLengthV0() const { return fDecayLengthV0; }
inline AliFemtoThreeVector AliFemtoV0::decayVertexV0() const { return fDecayVertexV0; } 
inline AliFemtoThreeVector AliFemtoV0::primaryVertex() const { return fPrimaryVertex; }
inline float AliFemtoV0::decayVertexV0X() const { return fDecayVertexV0.x(); } 
inline float AliFemtoV0::decayVertexV0Y() const { return fDecayVertexV0.y(); } 
inline float AliFemtoV0::decayVertexV0Z() const { return fDecayVertexV0.z(); } 
inline float AliFemtoV0::dcaV0Daughters() const { return fDcaV0Daughters; }
inline float AliFemtoV0::dcaV0ToPrimVertex() const { return fDcaV0ToPrimVertex; }
inline float AliFemtoV0::dcaPosToPrimVertex() const { return fDcaPosToPrimVertex; }
inline float AliFemtoV0::dcaNegToPrimVertex() const { return fDcaNegToPrimVertex; }
inline AliFemtoThreeVector AliFemtoV0::momPos() const { return fMomPos; }
inline float AliFemtoV0::momPosX() const { return fMomPos.x(); }
inline float AliFemtoV0::momPosY() const { return fMomPos.y(); }
inline float AliFemtoV0::momPosZ() const { return fMomPos.z(); }
inline AliFemtoThreeVector AliFemtoV0::momNeg() const { return fMomNeg; }
inline float AliFemtoV0::momNegX() const { return fMomNeg.x(); }
inline float AliFemtoV0::momNegY() const { return fMomNeg.y(); }
inline float AliFemtoV0::momNegZ() const { return fMomNeg.z(); }
inline AliFemtoThreeVector AliFemtoV0::momV0() const { return fMomV0; }
inline float AliFemtoV0::momV0X() const { return fMomV0.x(); }
inline float AliFemtoV0::momV0Y() const { return fMomV0.y(); }
inline float AliFemtoV0::momV0Z() const { return fMomV0.z(); }
inline float AliFemtoV0::alphaV0() const { return fAlphaV0; }
inline float AliFemtoV0::ptArmV0() const {return fPtArmV0;}
inline float AliFemtoV0::eLambda() const {return fELambda;}
inline float AliFemtoV0::eK0Short() const {return fEK0Short;}
inline float AliFemtoV0::ePosProton() const {return fEPosProton;}
inline float AliFemtoV0::ePosPion() const {return fEPosPion;}
inline float AliFemtoV0::eNegProton() const {return fENegProton;}
inline float AliFemtoV0::eNegPion() const {return fENegPion;}
inline float AliFemtoV0::massLambda() const {return fMassLambda;}
inline float AliFemtoV0::massAntiLambda() const {return fMassAntiLambda;}
inline float AliFemtoV0::massK0Short() const {return fMassK0Short;}
inline float AliFemtoV0::rapLambda() const {return fRapLambda;}
inline float AliFemtoV0::rapK0Short() const {return fRapK0Short;}
inline float AliFemtoV0::cTauLambda() const {return fCTauLambda;}
inline float AliFemtoV0::cTauK0Short() const {return fCTauK0Short;}
inline float AliFemtoV0::ptV0() const {return fPtV0;}
inline float AliFemtoV0::ptotV0() const {return fPtotV0;}
inline float AliFemtoV0::ptPos() const {return fPtPos;}
inline float AliFemtoV0::ptotPos() const {return fPtotPos;}
inline float AliFemtoV0::ptNeg() const {return fPtNeg;}
inline float AliFemtoV0::ptotNeg() const {return fPtotNeg;}
inline int   AliFemtoV0::tpcHitsPos() const { return fTpcHitsPos; }
inline int   AliFemtoV0::tpcHitsNeg() const { return fTpcHitsNeg; }
inline float AliFemtoV0::dedxNeg() const {return fDedxNeg;}
inline float AliFemtoV0::numdedxNeg() const {return fNufDedxNeg;} //Gael 04Fev2002
inline float AliFemtoV0::errdedxNeg() const {return fErrDedxNeg;} //Gael 04Fev2002
inline float AliFemtoV0::lendedxNeg() const {return fLenDedxNeg;} //Gael 04Fev2002
inline float AliFemtoV0::pseudoRapNeg() const {return fMomNeg.pseudoRapidity();} //Gael 04Fev2002
inline float AliFemtoV0::dedxPos() const {return fDedxPos;}
inline float AliFemtoV0::numdedxPos() const {return fNufDedxPos;} //Gael 04Fev2002
inline float AliFemtoV0::errdedxPos() const {return fErrDedxPos;} //Gael 04Fev2002
inline float AliFemtoV0::lendedxPos() const {return fLenDedxPos;} //Gael 04Fev2002
inline float AliFemtoV0::pseudoRapPos() const {return fMomPos.pseudoRapidity();} //Gael 04Fev2002


inline unsigned long   AliFemtoV0::trackTopologyMapPos(unsigned int word) const { return fTrackTopologyMapPos[word]; }
inline unsigned long   AliFemtoV0::trackTopologyMapNeg(unsigned int word) const { return fTrackTopologyMapNeg[word]; }
inline unsigned short   AliFemtoV0::idNeg() const { return fKeyNeg; } 
inline unsigned short   AliFemtoV0::keyNeg() const { return fKeyNeg; }
inline unsigned short   AliFemtoV0::idPos() const { return fKeyPos; } 
inline unsigned short   AliFemtoV0::keyPos() const { return fKeyPos; }

inline void AliFemtoV0::SetdecayLengthV0(const float x){ fDecayLengthV0= x;}   
inline void AliFemtoV0::SetdecayVertexV0X(const float x){ fDecayVertexV0.setX(x);}
inline void AliFemtoV0::SetdecayVertexV0Y(const float x){ fDecayVertexV0.setY(x);}
inline void AliFemtoV0::SetdecayVertexV0Z(const float x){ fDecayVertexV0.setZ(x);}
inline void AliFemtoV0::SetdecayVertexV0(const AliFemtoThreeVector v){ fDecayVertexV0 = v; }
inline void AliFemtoV0::SetdcaV0Daughters(const float x){fDcaV0Daughters= x;} 
inline void AliFemtoV0::SetdcaV0ToPrimVertex(const float x){fDcaV0ToPrimVertex= x;}   
inline void AliFemtoV0::SetdcaPosToPrimVertex(const float x){fDcaPosToPrimVertex = x;} 
inline void AliFemtoV0::SetdcaNegToPrimVertex(const float x){fDcaNegToPrimVertex = x;} 
inline void AliFemtoV0::SetmomPos(const AliFemtoThreeVector v){fMomPos = v; }
inline void AliFemtoV0::SetmomPosX(const float x){fMomPos.setX(x);}
inline void AliFemtoV0::SetmomPosY(const float x){fMomPos.setY(x);}
inline void AliFemtoV0::SetmomPosZ(const float x){fMomPos.setZ(x);}
inline void AliFemtoV0::SetmomNeg(const AliFemtoThreeVector v){fMomNeg = v; }
inline void AliFemtoV0::SetmomNegX(const float x){fMomNeg.setX(x);}
inline void AliFemtoV0::SetmomNegY(const float x){fMomNeg.setY(x);}
inline void AliFemtoV0::SetmomNegZ(const float x){fMomNeg.setZ(x);}
inline void AliFemtoV0::SetTrackTopologyMapPos(unsigned int word, const unsigned long& m){fTrackTopologyMapPos[word]=m;} 
inline void AliFemtoV0::SetTrackTopologyMapNeg(unsigned int word, const unsigned long& m){fTrackTopologyMapNeg[word]=m;} 
inline void AliFemtoV0::SetmomV0(AliFemtoThreeVector v){fMomV0= v; }
inline void AliFemtoV0::SetmomV0X(const float x){fMomV0.setX(x);}
inline void AliFemtoV0::SetmomV0Y(const float x){fMomV0.setY(x);}
inline void AliFemtoV0::SetmomV0Z(const float x){fMomV0.setZ(x);}

inline void AliFemtoV0::SetalphaV0( float x){fAlphaV0= x;}
inline void AliFemtoV0::SetptArmV0( float x){fPtArmV0 = x;}
inline void AliFemtoV0::SeteLambda( float x){fELambda= x;}       
inline void AliFemtoV0::SeteK0Short( float x){fEK0Short= x;}
inline void AliFemtoV0::SetePosProton( float x){fEPosProton= x;}      
inline void AliFemtoV0::SetePosPion( float x){fEPosPion= x;}      
inline void AliFemtoV0::SeteNegProton( float x){fENegProton= x;} 
inline void AliFemtoV0::SeteNegPion( float x){fENegPion= x;}       
inline void AliFemtoV0::SetmassLambda( float x){fMassLambda = x;} 
inline void AliFemtoV0::SetmassAntiLambda( float x){fMassAntiLambda= x;} 
inline void AliFemtoV0::SetmassK0Short( float x){fMassK0Short= x;}  
inline void AliFemtoV0::SetrapLambda( float x){fRapLambda= x;}
inline void AliFemtoV0::SetrapK0Short( float x){fRapK0Short = x;}   
inline void AliFemtoV0::SetcTauLambda( float x){fCTauLambda = x;}   
inline void AliFemtoV0::SetcTauK0Short( float x){fCTauK0Short = x;}   
inline void AliFemtoV0::SetptV0( float x){fPtV0 = x;}          
inline void AliFemtoV0::SetptotV0( float x){fPtotV0 = x;}
inline void AliFemtoV0::SetptPos( float x){fPtPos = x;}
inline void AliFemtoV0::SetptotPos( float x){fPtotPos = x;}    
inline void AliFemtoV0::SetptNeg( float x){ fPtNeg= x;}    
inline void AliFemtoV0::SetptotNeg( float x){ fPtotNeg= x;}
inline void AliFemtoV0::SetidNeg(const unsigned short& s){ fKeyNeg= s;}
inline void AliFemtoV0::SetidPos(const unsigned short& s){ fKeyPos= s;}
inline void AliFemtoV0::SetkeyNeg(const unsigned short& s){ fKeyNeg= s;}
inline void AliFemtoV0::SetkeyPos(const unsigned short& s){ fKeyPos= s;}
inline void AliFemtoV0::SettpcHitsPos(const int& i){fTpcHitsPos=i;} 
inline void AliFemtoV0::SettpcHitsNeg(const int& i){fTpcHitsNeg=i;}
inline void AliFemtoV0::SetdedxNeg(float x){fDedxNeg=x;}
inline void AliFemtoV0::SeterrdedxNeg(float x){fErrDedxNeg=x;}//Gael 04Fev2002
inline void AliFemtoV0::SetlendedxNeg(float x){fLenDedxNeg=x;}//Gael 04Fev2002
inline void AliFemtoV0::SetdedxPos(float x){fDedxPos=x;}
inline void AliFemtoV0::SeterrdedxPos(float x){fErrDedxPos=x;}//Gael 04Fev2002
inline void AliFemtoV0::SetlendedxPos(float x){fLenDedxPos=x;}//Gael 04Fev2002
inline void AliFemtoV0::SetprimaryVertex(const AliFemtoThreeVector v) { fPrimaryVertex = v; }//Gael 24 Sept 02
#endif


















