// -*- mode: c++ -*-
#ifndef ALIFMDANARING_H
#define ALIFMDANARING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
//____________________________________________________________________
//
// Utility class for analysing ESD data. 
// This class does sharing and background correction 
//
#include <TObject.h>
#include <TH2.h>
#include <TH1.h>
class TBrowser;

/** Base class for analysing FMD ESD data */
class AliFMDAnaRing : public TObject
{
public:
  /** Constructor */
  AliFMDAnaRing();
  /** Constructor 
      @param det   Detector 
      @param ring  Ring 
      @param bg    Background - not owned
      @param c0    Lower cut 
      @param c1    higher cut */
  AliFMDAnaRing(UShort_t det, Char_t ring, TH2* bg, Float_t c0, Float_t c1);
  /** Destructor */
  virtual ~AliFMDAnaRing() {}
  /** Called at beginning of run */
  virtual void Init() { fNEvents = 0; }
  /** Called at beginning of event */
  virtual void Begin() { fNEvents++; }
  /** Process ESD 
      @param phi Azimuthal angle @f$ \varphi @f$ of the hit (radians)
      @param eta Psuedo-rapidity @f$ \eta@f$ of hit
      @param m1  Multiplicity of this strip
      @param m2  Multiplicity of neighbor strip 
      @return @c true if hits are merged */
  virtual Bool_t ProcessESD(Float_t phi, Float_t eta, Float_t& m1, Float_t m2);
  /** User defined member function. 
      @param phi Azimuthal angle @f$ \varphi @f$ of the hit  (radians) 
      @param eta Psuedo-rapidity @f$ \eta@f$ of hit
      @param mult Corrected multiplicity */
  virtual void Fill(Float_t phi, Float_t eta, Float_t mult) = 0;
  /** Called at end of event */
  virtual void End() {}
  /** Called at end of run */
  virtual void Finish();
  /** Get the detector identifier 
      @return Detector identifier */
  UShort_t Detector() const { return fDet; }
  /** Get the ring identifier 
      @return Ring identifier */
  Char_t Ring() const { return fRing; }
  /** Get the number of seqments */ 
  UShort_t NSeq() const { return fNSeq; }
  /** Get the name 
      @return static string */ 
  const Char_t* Name() const { return fName; }
  /** Get the ring color 
      @return color */
  Int_t  Color() const;
  /** Browse this object */
  virtual void Browse(TBrowser* b);
  /** This is a folder */ 
  virtual Bool_t IsFolder() const { return kTRUE; }
protected: 
  /** Hidden copy constructor */
  AliFMDAnaRing(const AliFMDAnaRing&);
  /** Hidden Assignement operator */
  AliFMDAnaRing& operator=(const AliFMDAnaRing&);  
  /** Detector number */
  UShort_t fDet; 	// Detector number 
  /** Ring identifier */
  Char_t   fRing;	// Ring identifier 
  /** Name */
  Char_t   fName[6];	// Name 
  /** Background correction */
  TH2*     fBg;	// Background correction 
  /** Lower cut */
  Float_t  fCut0;	// Lower cut 
  /** Higher cut */
  Float_t  fCut1;	// Higher cut 
  /** Whether to use bacground correction */
  Bool_t   fUseBgCor;	// Whether to use bacground correction 
  /** # of segments */
  Int_t    fNSeq;	// # of segments 
  /** Histogram of multiplicity before any cuts */ 
  TH2D     fBareMult;
  /** Histogram of multiplicity after merging */ 
  TH2D     fMergedMult;
  /** Histogram of removed multiplicity */ 
  TH2D     fRemovedMult;
  /** Histogram of background corrected */ 
  TH2D     fMult;
  /** Step histogram 1 */
  TH2D     fStep1;
  /** Step histogram 1 */
  TH2D     fStep2;
  /** Step histogram 1 */
  TH2D     fStep3;

  /** Event counter */ 
  Int_t    fNEvents;
  ClassDef(AliFMDAnaRing,1) // Analysis of ESD in a ring 
};


#endif 
//
// EOF
//
