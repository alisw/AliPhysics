#ifndef ALIQADATAMAKERREC_H
#define ALIQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */
//
//  Base Class:
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  Y. Schutz CERN July 2007
//

// --- ROOT system ---
class TNtupleD ; 
// --- Standard library ---

// --- AliRoot header files ---
class AliDetectorRecoParam ;
#include "AliQADataMaker.h"
#include "AliQAv1.h"

class AliQADataMakerRec: public AliQADataMaker {
  
public:
	
	AliQADataMakerRec(const char * name="", const char * title="") ;          // ctor
	AliQADataMakerRec(const AliQADataMakerRec& qadm) ;   
	AliQADataMakerRec& operator = (const AliQADataMakerRec& qadm) ;
	virtual ~AliQADataMakerRec() ; // dtor
  
 	virtual Int_t Add2DigitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)    
  { return Add2List(hist, index, fDigitsQAList, expert, image) ; }
  virtual Int_t Add2ESDsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)                      
    { return Add2List(hist, index, fESDsQAList, expert, image) ; }
  virtual Int_t Add2HitsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/, const Bool_t /*image = kFALSE*/)      
    { return -1 ; }  
  virtual Int_t Add2RecPointsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)                 
    { return Add2List(hist, index, fRecPointsQAList, expert, image) ; }
  virtual Int_t Add2RawsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE, const Bool_t saveForCorr = kFALSE)  
  { return Add2List(hist, index, fRawsQAList, expert, image, saveForCorr) ; }
  virtual Int_t Add2SDigitsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/, const Bool_t /*image = kFALSE*/)   { return -1 ; } 
	
  virtual void        Exec(AliQAv1::TASKINDEX_t task, TObject * data) ;
  virtual void        EndOfCycle() ;
  virtual void        EndOfCycle(AliQAv1::TASKINDEX_t task) ;
  virtual void        EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** ) {AliInfo("To be implemented by detectors");} 
  virtual const AliDetectorRecoParam * GetRecoParam() { return fRecoParam ; }

  virtual TObject*    GetDigitsData(const Int_t index   )   { return GetData(fDigitsQAList, index); } 
  virtual TObject*    GetESDsData(const Int_t index)        { return GetData(fESDsQAList, index); }
  virtual TObject*    GetHitsData(const Int_t /*index*/)    { return NULL ; }
  virtual TObject*    GetRecPointsData(const Int_t index)   { return GetData(fRecPointsQAList, index); }
  virtual TObject*    GetRawsData(const Int_t index)        { return GetData(fRawsQAList, index); }
  virtual TObject*    GetSDigitsData(const Int_t /*index*/) { return NULL; }  
  //
  virtual TH1*        GetDigitsData(const Int_t index, int cloneID)      { return GetData(fDigitsQAList, index, cloneID); } 
  virtual TH1*        GetESDsData(const Int_t index, int cloneID)        { return GetData(fESDsQAList, index, cloneID); }
  virtual TH1*        GetHitsData(const Int_t /*index*/, int)            { return NULL ; }
  virtual TH1*        GetRecPointsData(const Int_t index, int cloneID)   { return GetData(fRecPointsQAList, index, cloneID); }
  virtual TH1*        GetRawsData(const Int_t index, int cloneID)        { return GetData(fRawsQAList, index, cloneID); }
  virtual TH1*        GetSDigitsData(const Int_t /*index*/, int)         { return NULL; }  
  //
  virtual TObjArray*  GetDigitsDataOfTrigClass(int cloneID, TObjArray *dest=0)           {return GetDataOfTrigClass(fDigitsQAList,cloneID,dest);}
  virtual TObjArray*  GetSDigitsDataOfTrigClass(int /*cloneID*/, TObjArray */*dest*/=0)  {return NULL;}
  virtual TObjArray*  GetESDsDataOfTrigClass(int cloneID, TObjArray *dest=0)             {return GetDataOfTrigClass(fESDsQAList,cloneID,dest);}
  virtual TObjArray*  GetHitsDataOfTrigClass(int /*cloneID*/, TObjArray */*dest*/=0)     {return NULL;}
  virtual TObjArray*  GetRecPointsDataOfTrigClass(int cloneID, TObjArray *dest=0)        {return GetDataOfTrigClass(fRecPointsQAList,cloneID,dest);}
  virtual TObjArray*  GetRawsDataOfTrigClass(int cloneID, TObjArray *dest=0)             {return GetDataOfTrigClass(fRawsQAList,cloneID,dest);}
  //
  virtual TObjArray** Init(AliQAv1::TASKINDEX_t task, Int_t cycles = -1) ;
  virtual void        Init(AliQAv1::TASKINDEX_t task, TObjArray ** list, Int_t run, Int_t cycles = -1) ;
  virtual void        InitRaws()                          {AliInfo("To be implemented by detectors");}
  virtual void        InitRecPoints()                     {AliInfo("To be implemented by detectors");}
  virtual void        InitDigits()                        {AliInfo("To be implemented by detectors");}
  virtual void        InitESDs()                          {AliInfo("To be implemented by detectors");}
  virtual void        ResetDetector(AliQAv1::TASKINDEX_t task) ;
  virtual void        StartOfCycle(Int_t run = -1) ;
  virtual void        StartOfCycle(AliQAv1::TASKINDEX_t task, Int_t run, const Bool_t sameCycle = kFALSE) ;
  virtual void        SetRecoParam(const AliDetectorRecoParam *param) { fRecoParam = param; }
  //
  virtual TObjArray* GetMatchingHitsData(const Int_t, TObjArray*)                      {return 0;}
  virtual TObjArray* GetMatchingSDigitsData(const Int_t, TObjArray*)                   {return 0;}
  virtual TObjArray* GetMatchingDigitsData(const Int_t index, TObjArray* optDest=0)    {return GetMatchingHistos(fDigitsQAList,index,optDest);}
  virtual TObjArray* GetMatchingRawsData(const Int_t index, TObjArray* optDest=0)      {return GetMatchingHistos(fRawsQAList,index,optDest);}
  virtual TObjArray* GetMatchingRecPointsData(const Int_t index, TObjArray* optDest=0) {return GetMatchingHistos(fRecPointsQAList,index,optDest);}
  virtual TObjArray* GetMatchingESDsData(const Int_t index, TObjArray* optDest=0)      {return GetMatchingHistos(fESDsQAList,index,optDest);}
  //
  virtual TH1*       GetMatchingHitsHisto(Int_t , Int_t )                  {return 0;}
  virtual TH1*       GetMatchingSDigitsHisto(Int_t , Int_t )               {return 0;}
  virtual TH1*       GetMatchingDigitsHisto(Int_t index, Int_t trigId)     {return GetMatchingHisto(fDigitsQAList,index,trigId);}
  virtual TH1*       GetMatchingRawsHisto(Int_t index, Int_t trigId)       {return GetMatchingHisto(fRawsQAList,index,trigId);}
  virtual TH1*       GetMatchingRecPointsHisto(Int_t index, Int_t trigId)  {return GetMatchingHisto(fRecPointsQAList,index,trigId);}
  virtual TH1*       GetMatchingESDsHisto(Int_t index, Int_t trigId)       {return GetMatchingHisto(fESDsQAList,index,trigId);}
  //
  virtual TObjArray* GetMatchingHitsHistosSet(const Int_t*, Int_t ,Int_t)                            {return 0;}
  virtual TObjArray* GetMatchingSDigitsHistosSet(const Int_t* , Int_t ,Int_t)                        {return 0;}
  virtual TObjArray* GetMatchingDigitsHistosSet(const Int_t* indexList, Int_t nHist,Int_t trigId)    {return GetMatchingHistosSet(fDigitsQAList,indexList,nHist,trigId);}
  virtual TObjArray* GetMatchingRawsHistosSet(const Int_t* indexList, Int_t nHist,Int_t trigId)      {return GetMatchingHistosSet(fRawsQAList,indexList,nHist,trigId);}
  virtual TObjArray* GetMatchingRecPointsHistosSet(const Int_t* indexList, Int_t nHist,Int_t trigId) {return GetMatchingHistosSet(fRecPointsQAList,indexList,nHist,trigId);}
  virtual TObjArray* GetMatchingESDsHistosSet(const Int_t* indexList, Int_t nHist,Int_t trigId)      {return GetMatchingHistosSet(fESDsQAList,indexList,nHist,trigId);}
  //
  virtual Int_t  FillHitsData(Int_t, double )                          {return -1;}
  virtual Int_t  FillSDigitsData(Int_t, double)                        {return -1;}
  virtual Int_t  FillDigitsData(Int_t index, double x)                 {return FillData(fDigitsQAList, index, x);}
  virtual Int_t  FillRawsData(Int_t index, double x)                   {return FillData(fRawsQAList, index, x);}
  virtual Int_t  FillRecPointsData(Int_t index, double x)              {return FillData(fRecPointsQAList, index, x);}
  virtual Int_t  FillESDsData(Int_t index, double x)                   {return FillData(fESDsQAList, index, x);}
  //
  virtual Int_t  FillHitsData(Int_t, double, double)                   {return -1;}
  virtual Int_t  FillSDigitsData(Int_t, double, double)                {return -1;}
  virtual Int_t  FillDigitsData(Int_t index, double x, double y)       {return FillData(fDigitsQAList, index, x, y);}
  virtual Int_t  FillRawsData(Int_t index, double x, double y)         {return FillData(fRawsQAList, index, x, y);}
  virtual Int_t  FillRecPointsData(Int_t index, double x, double y)    {return FillData(fRecPointsQAList, index, x, y);}
  virtual Int_t  FillESDsData(Int_t index, double x, double y)         {return FillData(fESDsQAList, index, x, y);}
  //
  virtual Int_t  FillHitsData(Int_t, double, double, double)                     {return -1;}
  virtual Int_t  FillSDigitsData(Int_t, double, double, double)                  {return -1;}
  virtual Int_t  FillDigitsData(Int_t index, double x, double y, double z)       {return FillData(fDigitsQAList, index, x,y,z);}
  virtual Int_t  FillRawsData(Int_t index, double x, double y, double z)         {return FillData(fRawsQAList, index, x,y,z);}
  virtual Int_t  FillRecPointsData(Int_t index, double x, double y, double z)    {return FillData(fRecPointsQAList, index, x,y,z);}
  virtual Int_t  FillESDsData(Int_t index, double x, double y, double z)         {return FillData(fESDsQAList, index, x,y,z);}
  //
  virtual Int_t  SetHitsDataBinContent(Int_t, int, double)                       {return -1;}
  virtual Int_t  SetSDigitsDataBinContent(Int_t, int, double)                    {return -1;}
  virtual Int_t  SetDigitsDataBinContent(Int_t index, int bin, double w)         {return SetDataBinContent(fDigitsQAList, index,bin,w);}
  virtual Int_t  SetRawsDataBinContent(Int_t index, int bin, double w)           {return SetDataBinContent(fRawsQAList, index,bin,w);}
  virtual Int_t  SetRecPointsDataBinContent(Int_t index, int bin, double w)      {return SetDataBinContent(fRecPointsQAList, index,bin,w);}
  virtual Int_t  SetESDsDataBinContent(Int_t index, int bin, double w)           {return SetDataBinContent(fESDsQAList, index,bin,w);}
  //
  virtual Int_t  SetHitsDataBinContent(Int_t, int, int, double)                         {return -1;}
  virtual Int_t  SetSDigitsDataBinContent(Int_t, int, int, double)                      {return -1;}
  virtual Int_t  SetDigitsDataBinContent(Int_t index, int binX, int binY, double w)     {return SetDataBinContent(fDigitsQAList, index,binX,binY,w);}
  virtual Int_t  SetRawsDataBinContent(Int_t index, int binX, int binY, double w)       {return SetDataBinContent(fRawsQAList, index,binX,binY,w);}
  virtual Int_t  SetRecPointsDataBinContent(Int_t index, int binX, int binY, double w)  {return SetDataBinContent(fRecPointsQAList, index,binX,binY,w);}
  virtual Int_t  SetESDsDataBinContent(Int_t index, int binX, int binY, double w)       {return SetDataBinContent(fESDsQAList, index,binX,binY,w);}
  //
  virtual Int_t  SetHitsDataBinError(Int_t, int, double)                        {return -1;}
  virtual Int_t  SetSDigitsDataBinError(Int_t, int, double)                     {return -1;}
  virtual Int_t  SetDigitsDataBinError(Int_t index, int bin, double err)        {return SetDataBinError(fDigitsQAList, index,bin,err);}
  virtual Int_t  SetRawsDataBinError(Int_t index, int bin, double err)          {return SetDataBinError(fRawsQAList, index,bin,err);}
  virtual Int_t  SetRecPointsDataBinError(Int_t index, int bin, double err)     {return SetDataBinError(fRecPointsQAList, index,bin,err);}
  virtual Int_t  SetESDsDataBinError(Int_t index, int bin, double err)          {return SetDataBinError(fESDsQAList, index,bin,err);}
  //
  virtual Int_t  SetHitsDataBinError(Int_t, int, int, double)                              {return -1;}
  virtual Int_t  SetSDigitsDataBinError(Int_t, int, int, double)                           {return -1;}
  virtual Int_t  SetDigitsDataBinError(Int_t index, int binX, int binY, double err)        {return SetDataBinError(fDigitsQAList, index,binX,binY,err);}
  virtual Int_t  SetRawsDataBinError(Int_t index, int binX, int binY, double err)          {return SetDataBinError(fRawsQAList, index,binX,binY,err);}
  virtual Int_t  SetRecPointsDataBinError(Int_t index, int binX, int binY, double err)     {return SetDataBinError(fRecPointsQAList, index,binX,binY,err);}
  virtual Int_t  SetESDsDataBinError(Int_t index, int binX, int binY, double err)          {return SetDataBinError(fESDsQAList, index,binX,binY,err);}
  //
  virtual Int_t  ResetHitsData(Int_t, Option_t*)                      {return -1;}
  virtual Int_t  ResetSDigitsData(Int_t, Option_t*)                   {return -1;}
  virtual Int_t  ResetDigitsData(Int_t index, Option_t* opt="")       {return ResetData(fDigitsQAList, index, opt);}
  virtual Int_t  ResetRawsData(Int_t index, Option_t* opt="")         {return ResetData(fRawsQAList, index, opt);}
  virtual Int_t  ResetRecPointsData(Int_t index, Option_t* opt="")    {return ResetData(fRecPointsQAList, index, opt);}
  virtual Int_t  ResetESDsData(Int_t index, Option_t* opt="")         {return ResetData(fESDsQAList, index, opt);}
  //
  virtual Int_t  ResetStatsHitsData(Int_t)                            {return -1;}
  virtual Int_t  ResetStatsSDigitsData(Int_t)                         {return -1;}
  virtual Int_t  ResetStatsDigitsData(Int_t index)                    {return ResetStatsData(fDigitsQAList, index);}
  virtual Int_t  ResetStatsRawsData(Int_t index)                      {return ResetStatsData(fRawsQAList, index);}
  virtual Int_t  ResetStatsRecPointsData(Int_t index)                 {return ResetStatsData(fRecPointsQAList, index);}
  virtual Int_t  ResetStatsESDsData(Int_t index)                      {return ResetStatsData(fESDsQAList, index);}
  //
  virtual void   ClonePerTrigClass(AliQAv1::TASKINDEX_t task);
  //
protected: 

  virtual void   InitRecoParams() ; 
  virtual void   InitHits()                          {AliWarning("Call not valid") ; }
  //virtual void   InitRecParticles()                {AliInfo("To be implemented by detectors");}
  virtual void   InitSDigits()                       {AliWarning("Call not valid") ; }
  //virtual void   InitTrackSegments()               {AliInfo("To ne implemented by detectors");}
  virtual void   MakeESDs(AliESDEvent * )            {AliInfo("To be implemented by detectors");} 
  virtual void   MakeHits()                          {AliWarning("Call not valid") ; }
  virtual void   MakeHits(TTree * )                  {AliWarning("Call not valid") ; }  
  virtual void   MakeDigits()                        {AliInfo("To be implemented by detectors");}   
  virtual void   MakeDigits(TTree * )                {AliInfo("To be implemented by detectors");}   
  //virtual void   MakeRecParticles()                {AliInfo("To be implemented by detectors");} 
  virtual void   MakeRaws(AliRawReader *)            {AliInfo("To be implemented by detectors");} 
  virtual void   MakeRecPoints(TTree * )             {AliInfo("To be implemented by detectors");} 
  virtual void   MakeSDigits()                       {AliWarning("Call not valid") ; }     
  virtual void   MakeSDigits(TTree * )               {AliWarning("Call not valid") ; }    
  virtual void   StartOfDetectorCycle()              {AliInfo("To be implemented by detectors");} 
  
  TObjArray * *               fDigitsQAList ;     //! list of the digits QA data objects
  TObjArray * *               fESDsQAList ;       //! list of the ESDs QA data objects
  TObjArray * *               fRawsQAList ;       //! list of the raws QA data objects
  TObjArray * *               fRecPointsQAList ;  //! list of the RecPoints QA data objects
  TNtupleD  **                fCorrNt ;           //! This is used by Corr only to hold its Ntuple. 
  const AliDetectorRecoParam *fRecoParam;         //! const pointer to the reco parameters to be used in the reco QA
  TClonesArray *              fRecPointsArray;    //! Array that contains the RecPoints    
  
 ClassDef(AliQADataMakerRec,4)  // description 

};

#endif // ALIQADATAMAKERREC_H
