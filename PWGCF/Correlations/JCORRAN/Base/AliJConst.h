/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Comment describing what this class does needed!

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

#define JUNUSED(expr) do { (void)(expr); } while (0)

//
//  Constants
//
// for JHisto and PhxJHisto
#define kMaxNoCentrBin 10   // Maximum no of centrality bins defined in JCard.h
#define kMaxTriggParticleType 35           // 
#define kMaxJetType  32          // 
#define kPtDim 15           // 
#define kMaxNoRuns 100      // Maximum no of runs in a nanoDST file
#define kMaxNoEPBin 24

const double kJPi           = 3.14159265358979; //TMath::Pi();
const double kJTwoPi        = 2*kJPi;
const double kJToRadian     = kJPi/180.;
const double kJToDegree     = 180./kJPi;
const double kElectronMass = .51099906e-3;
const double kPionMass     = .1395675;


//====================== particle types ============================================
const int kNumberOfParticleTypes = 10;
enum particleType   {kJHadron, kJPion, kJKaon, kJProton, kJPhoton, kJDecayphoton, kJPizero, kJEta,kJHadronMC, kJJet, kJNone};

const char* const kParticleTypeStrName[kNumberOfParticleTypes] =
                    {"hadron", "pion", "kaon", "proton", "photon", "decayphoton", "pizero", "eta", "hadronMC", "none"};
const char* const kParticleProtoType[kNumberOfParticleTypes] =
                     {"AliJTrack", "AliJTrack", "AliJTrack", "AliJTrack", "AliJPhoton", "AliJPhoton", "AliJPiZero", "AliJPiZero", "None"};

//=======================JCorran trigger table definition===========================
//internal JCorran trigger mask  TBit=0 is MinBias, TBit=1 HighMultiplicityTrigger
enum TriggerBitJCorran   {kMinBiasTriggerBitJCorran, kHighMultTriggerBitJCorran,
                          kEmc0TriggerBitJCorran,
                          kEmc1GammaTriggerBitJCorran, kEmc1JetTriggerBitJCorran,
                          kCentralTriggerBitJCorran,kSemiCentralTriggerBitJCorran,
                          kFastOnlyBitJCorran,
						  kINT7TriggerBitJCorran, kJNTriggerBit
                          };  
const int kRangeTriggerTableAlice   = 50;
const int kRangeTriggerTableJCorran = 16;

//==================================================================================
enum fillType { kReal, kMixed, kEtaGap };
enum corrFillType {kAzimuthFill=0,kPionFill=1};
enum corrType { kTriggType, kAssocType, kXeType, kLongType, kCentrType, kZVertType, kMassType, kEtaGapType, kDiJetType, kRGapType, kNoType, kNcorrType };
enum TriggerParticleType { kTriggParticles, kLeadingParticle, kIsolatedParticle }; 

const char* const kTriggerParticleTypeName[] = 
      {"TriggParticles","LeadingParticle","IsolatedParticle"};

// JETs =====
const int kNJetAlg = 10;
enum JetAlg   {kkt,kantikt,ksiscone,krecomE,krecomB,kcdfmidpoint,kjade,kd0run2cone,kGF,kSimpleCone};
const char* const kJetAlgStrName[kNJetAlg] =
                    {"kt","antikt","siscone","recomE","recomB","cdfmidpoint","jade","d0run2cone","GF","SimpleCone"};

enum EPType { kEPV0A, kEPV0C, kEPV0AC, kNEPType };
const int kNHarmonics = 5;
// PHENIX  constants
enum TEMC {kPbSc, kPbGl}; 

#endif
