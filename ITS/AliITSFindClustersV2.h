#ifndef ALIITSFINDCLUSTERSV2_H
#define ALIITSFINDCLUSTERSV2_H
/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */
 
/*
  $Id$
 */
///////////////////////////////////////////////////////////////////
//Class for reconstruction of clusters V2                        //
///////////////////////////////////////////////////////////////////
#include <TTask.h>

class TString;
class TFile;
class AliRun;
class AliITSgeom;

class AliITSFindClustersV2 : public TTask{
 public:
    AliITSFindClustersV2(); // default constructor
    // standard constructor files not opened, by file name
    AliITSFindClustersV2(const TString infile,const TString outfile = "");
    // standard constructor for files already opened.
    AliITSFindClustersV2(TFile *in,TFile *out=0);
    // Standard constructor for AliRun already read in.
    AliITSFindClustersV2(AliRun *ar, const TString outfile = "");
    AliITSFindClustersV2(const AliITSFindClustersV2& rec);
    AliITSFindClustersV2& operator=(const AliITSFindClustersV2 &source);

    virtual ~AliITSFindClustersV2();//Destructor
    virtual Bool_t FastSimulation() const {return fSlowFast;}
    virtual void SetSlowSimulation(){fSlowFast = kFALSE;}
    virtual void SetFastSimulation(){fSlowFast = kTRUE;}
    virtual void Exec(const Option_t *opt="");
 private:
    AliRun *fAr;           //! Pointer of type AliRun
    Bool_t fDeletfAr;      //! Logical to indecate if fAr should be deleted.
    AliITSgeom *fGeom;     //! Pointer to ITS geometry
    TString *fInFileName;  //! Pointer to input file name string.
    TString *fOutFileName; //! Pointer to output file name string.
    TFile   *fIn;          //! pointer to input file
    TFile   *fOut;         //! pointer to output file
    Bool_t fInit;          //! true if Init was successfull
    Bool_t fSlowFast;      //! if true then using fast ITS simulation.

    ClassDef(AliITSFindClustersV2,1) // Task to Reconstruct ITS from Digits.

};
#endif
