#ifndef ALIJCORRANSETUP_H
#define ALIJCORRANSETUP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//______________________________________________________________________________
//                 Implementation Class AliJCorranSetup.h
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////


#include <string>

#define kMaxNumOfBadRuns 100
#define kMaxNumOfFiles 500
#define kMaxNumOfQual  5

using namespace std;

class AliJCORRANSetup {

public:
  AliJCORRANSetup(); // default constructor
  AliJCORRANSetup(const char *SelectionsFile);
  virtual ~AliJCORRANSetup(){;}

  void PrintOut();
  void WriteToFile(const char *fout);
  
  // GENERAL SELECTIONS:
  string  GetSpecies()  const     {return fSpecies;}
  float GetVtxRange()   const     {return fVtxRange;}
  // Trigger SELECTIONS:
  float GetOffLineTriggThreshold() const {return fOffLineTriggThreshold;}
  float GetOffLinePi0TriggThreshold() const  {return fOffLinePi0TriggThreshold;}
  float GetOffLineHadronTriggThreshold() const {return fOffLineHadronTriggThreshold;}
  float GetOffLineAccHadronTriggThreshold() const {return fOffLineAccHadronTriggThreshold;}
  int   GetMinBiasScaleFact() const  {return fMinBiasScaleFact;}
  // CLUSTER SELECTIONS:
  float GetClusEnerMin() const    {return fClusEnerMin;}       
  float GetClusEnerMax() const    {return fClusEnerMax;}       
  float GetChi2Cut()     const    {return fChi2Cut;}      
  float GetProbPhotCut() const    {return fProbPhotCut;}
  int   GetClusNTowMin() const    {return fClusNTowMin;}
  int   GetClusNTowMax() const    {return fClusNTowMax;}
  // TRACK SELECTIONS:
  float GetTrkPtMin()    const    {return fTrkPtMin;}         
  float GetTrkPtMax()    const    {return fTrkPtMax;}
  // External file
  string  GetWarnMapFile()  const {return fWarnmapfile;}

protected:

  // GENERAL SELECTIONS:
  string fSpecies;  //name
  float fVtxRange;  //vertex range
  // Trigger SELECTIONS: 
  float fOffLineTriggThreshold;     //offline trigger threshold
  float fOffLinePi0TriggThreshold;  //offline pi0 trigger threshold
  float fOffLineHadronTriggThreshold;  //offline hadron trigger threshold
  float fOffLineAccHadronTriggThreshold; //offline acc hadron trigg thr
  int   fMinBiasScaleFact; //min bias scale factor
  // CLUSTER SELECTIONS:
  float fClusEnerMin; //minimal cluster energy
  float fClusEnerMax; //maximal cluster energy
  float fChi2Cut;     //chi2 cut
  float fProbPhotCut; //probability photon cut
  int   fClusNTowMin; //minimum of N towers
  int   fClusNTowMax;  //maximum of N towers
  // TRACK SELECTIONS:
  float fTrkPtMin;   //minimum pt track
  float fTrkPtMax;   //maximum pt track       
  // External file
  string fWarnmapfile;  //warn map file
  string fComment;    //comment
};
#endif // AliJCORRANSetup_H
