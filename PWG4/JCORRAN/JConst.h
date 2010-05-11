// $Id: JConst.h,v 1.5 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file JConst.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.5 $
  \date $Date: 2008/05/08 15:19:52 $
*/
////////////////////////////////////////////////////

#ifndef JCONST_H
#define JCONST_H

//
//  Constants
//
// for JHisto and PhxJHisto
#define kMaxNoCentrBin 10   // Maximum no of centrality bins defined in JCard.h
#define kPtDim 20           // 
#define kMaxNoRuns 100      // Maximum no of runs in a nanoDST file

const double kJPi           = 3.14159265358979; //TMath::Pi();
const double kJTwoPi        = 2*kJPi;
const double kJToRadian     = kJPi/180.;
const double kJToDegree     = 180./kJPi;
const double kElectronMass = .51099906e-3;
const double kPionMass     = .1395675;

enum expName {kPHENIX, kALICE};
const int kNumberOfExperiments = 2; //numberOfExperiments

//====================== particle types ============================================
const int kNumberOfParticleTypes = 9;
enum particleType   {kHadron, kPion, kKaon, kProton, kPhoton, kDecayphoton, kPizero, kEta, kNone};

const char* const kParticleTypeStrName[kNumberOfParticleTypes] =
                    {"hadron", "pion", "kaon", "proton", "photon", "decayphoton", "pizero", "eta", "none"};
const char* const kParticleProtoType[kNumberOfExperiments][kNumberOfParticleTypes] =
                    {{"PhJCgl",    "PhJCgl",    "PhJCgl",    "PhJCgl",    "PhJPhoton",  "PhJPhoton",  "AliPhJPiZero", "AliPhJPiZero", "None"},
                     {"AliJTrack", "AliJTrack", "AliJTrack", "AliJTrack", "AliJPhoton", "AliJPhoton", "AliPhJPiZero", "AliPhJPiZero", "None"}};


//=======================JCorran trigger table definition===========================
enum TriggerBitJCorran   {kMinBiasTriggerBitJCorran};  //internal JCorran trigger mask  TBit=0 is MinBias
const int kRangeTriggerTableAlice   = 50;
const int kRangeTriggerTableJCorran = 16;

//==================================================================================
enum fillType { kReal, kMixed, kRotated };
enum corrType { kTriggType, kAssocType, kXeType, kCentrType, kMassType, kNoType }; 


// PHENIX  constants
enum TEMC {kPbSc, kPbGl}; 

#endif
