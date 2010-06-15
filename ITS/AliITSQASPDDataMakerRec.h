#ifndef ALIITSQASPDDATAMAKERREC_H
#define ALIITSQASPDDATAMAKERREC_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//
//
//  W. Ferrarese + P. Cerello Feb 2008
//  M. Nicassio D. Elia INFN Bari March 2008
//  maria.nicassio@ba.infn.it

/* $Id$  */

class TObjArray;
class AliRawReader;
class AliITSRawStreamSPDErrorLog;
class AliITSQADataMakerRec;
class AliQAv1;

class AliITSQASPDDataMakerRec : public TObject {

public:
  AliITSQASPDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode = kFALSE, Short_t ldc = 0,
                          AliITSRawStreamSPDErrorLog *aliITSRawStreamSPDErrorLog = NULL); //ctor
  AliITSQASPDDataMakerRec(const AliITSQASPDDataMakerRec& qadm);
  AliITSQASPDDataMakerRec& operator = (const AliITSQASPDDataMakerRec& qac);
  virtual Int_t InitRaws();
  virtual Int_t InitDigits();
  virtual Int_t InitRecPoints();
  virtual Int_t MakeRaws(AliRawReader *rawReader);
  virtual Int_t MakeRecPoints(TTree *clustersTree);
  virtual Int_t MakeDigits()  {return 0;}
  virtual Int_t MakeDigits(TTree *clustersTree);
  virtual void  StartOfDetectorCycle();
  virtual void  EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASPDDataMakerRec();   // dtor
  Int_t GetOffset(AliQAv1::TASKINDEX_t task,Int_t specie=0) const;
  void  SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset, Int_t specie = 0);
  Int_t GetTaskHisto(AliQAv1::TASKINDEX_t task) const;
  virtual void ResetDetector(AliQAv1::TASKINDEX_t){;};

  enum {kAmoreFoOffset=10, kAmoreErrorsOffset=21};

private: 

  static const Short_t fgknSPDmodules = 240;    //number of SPD modules
  static const Short_t fgkLADDonLay1  = 80;     //number of modules on layer 1
  static const Short_t fgkLADDonLay2  = 160;    //number of modules on layer 2
  static const Short_t fgkSPDchips    = 1200;   //number of chips

  AliITSQADataMakerRec *fAliITSQADataMakerRec;//pointer to the main ctor
  Bool_t  fkOnline;                           //online (1) or offline (0) use
  Int_t   fLDC;                               //LDC number (0 for offline, 1 to 4 for online) 
  Int_t   fSPDhRawsTask;                      // number of booked SPD histograms for the Raws Task
  Int_t   fSPDhDigitsTask;                    // number of booked SPD histograms for the RecPoints Task
  Int_t   fSPDhRecPointsTask;                 // number of booked SPD histograms for the RecPoints Task
  Int_t   *fGenRawsOffset;                     // QAchecking Raws offset
  Int_t   *fGenDigitsOffset;                   // QAchecking Digits offset
  Int_t   *fGenRecPointsOffset;                // QAchecking RecPoints offset
  AliITSRawStreamSPDErrorLog *fAdvLogger;  // pointer to special error logger object

  ClassDef(AliITSQASPDDataMakerRec,6)      // description 

};

#endif


