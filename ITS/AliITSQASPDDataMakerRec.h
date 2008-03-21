#ifndef AliITSQASPDDataMakerRec_H
#define AliITSQASPDDataMakerRec_H
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

/* $Id:$  */

class TObjArray;
class AliRawReader;
class AliITSQADataMakerRec;
class AliQA;

class AliITSQASPDDataMakerRec : public TObject {

public:
  AliITSQASPDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode = kFALSE, Short_t ldc = 0); //ctor
  AliITSQASPDDataMakerRec(const AliITSQASPDDataMakerRec& qadm);
  AliITSQASPDDataMakerRec& operator = (const AliITSQASPDDataMakerRec& qac);
  virtual void InitRaws();
  virtual void InitRecPoints();
  virtual void MakeRaws(AliRawReader *rawReader);
  virtual void MakeRecPoints(TTree *clustersTree);
  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASPDDataMakerRec();   // dtor
  Int_t Raws() const { return fSPDhRaws; }
  Int_t Recs() const { return fSPDhRecs; }

private: 

  static const Int_t fgknSPDmodules = 240;   //number of SPD modules
  static const Int_t fgkLADDonLay1 = 80;     //number of modules on layer 1
  static const Int_t fgkLADDonLay2 = 160;    //number of modules on layer 2


  AliITSQADataMakerRec *fAliITSQADataMakerRec;//pointer to the main ctor
  Bool_t  fkOnline;                           //online (1) or offline (0) use
  Int_t   fLDC;                               //LDC number (0 for offline, 1 to 4 for online) 
  Int_t   fSPDhRaws;                          //number of booked SPD Raws histograms;
  Int_t   fSPDhRecs;                          //number of booked SPD Recs histograms;
  Int_t   fRawsOffset;                        // number of histo booked when SPD start 
  Int_t   fRecsOffset;                        // number of histo booked when SPD start
  

  ClassDef(AliITSQASPDDataMakerRec,2)      // description 

};

#endif

