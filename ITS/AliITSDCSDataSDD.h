#ifndef ALIITSDCSDATASDD_H
#define ALIITSDCSDATASDD_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to define object containing SDD DCS data                //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//         V.Pospisil, CTU Prague, gdermog@seznam.cz             //
///////////////////////////////////////////////////////////////////

#include<TObject.h>
#include<TArrayF.h>
#include<TArrayI.h>
#include<TArrayC.h>

class AliITSDCSDataSDD : public TObject 
{ 

 public:
  AliITSDCSDataSDD( void );
                        // Default constructor
  ~AliITSDCSDataSDD( void ){};
                        // Destructor is void

  void SetNPointsTempLeft( Int_t npts );
  void SetNPointsTempRight( Int_t npts );
  void SetNPointsHV( Int_t npts );
  void SetNPointsMV( Int_t npts );
  void SetNPointsStatus( Int_t npts );
                        // Sets sizes of the DCS variable arrays

  void SetValueTempLeft(Int_t time, Float_t temperature );
  void SetValueTempRight(Int_t time, Float_t temperature );
  void SetValueHV(Int_t time, Float_t voltage );
  void SetValueMV(Int_t time, Float_t voltage );
  void SetValueStatus(Int_t time, Char_t status );
                        // Inserts value of a DCS variable into the appropriate array.
                        //  Resizes and sorts array if necessary.

  void Compress();      // Tries to minimize array sizes

  Int_t GetTempLeftRecords()  const {return fTempLeftSetPoints;}
  Int_t GetTempRightRecords() const {return fTempRightSetPoints;}
  Int_t GetHVRecords()        const {return fHVSetPoints;}
  Int_t GetMVRecords()        const {return fMVSetPoints;}
  Int_t GetStatusRecords()    const {return fStatusSetPoints;}
                       // Returns number of stored values of specific DCS variable

  Int_t GetTempLeftSize()  const {return fTempLeftMaxPoints;}
  Int_t GetTempRightSize() const {return fTempRightMaxPoints;}
  Int_t GetHVSize()        const {return fHVMaxPoints;}
  Int_t GetMVSize()        const {return fMVMaxPoints;}
  Int_t GetStatusSize()    const {return fStatusMaxPoints;}
                       // Returns size of selected array

  Float_t GetTempLeftIdx( Int_t index )  const {return fTempLeft.At(index);};
  Float_t GetTempRightIdx( Int_t index ) const {return fTempRight.At(index);};
  Float_t GetHVIdx( Int_t index )        const {return fHV.At(index);};
  Float_t GetMVIdx( Int_t index )        const {return fMV.At(index);};
  Char_t  GetStatusIdx( Int_t index )    const {return fStatus.At(index);};
                       // Returns value of specific DCS variable by index in the array

  Bool_t GetOKStatIdx( Int_t index )        const { return (Bool_t)(fStatus.At(index) & 1 ); };
  Bool_t GetTempLeftStatIdx( Int_t index )  const { return (Bool_t)(fStatus.At(index) & 2 ); };
  Bool_t GetTempRightStatIdx( Int_t index ) const { return (Bool_t)(fStatus.At(index) & 4 ); };
                       // Return status of a readout device by index in the array

  Int_t   GetTempLeftTimeIdx( Int_t index )  const {return fTempLeftTimeStamp.At(index);};
  Int_t   GetTempRightTimeIdx( Int_t index ) const {return fTempRightTimeStamp.At(index);};
  Int_t   GetHVTimeIdx( Int_t index )        const {return fHVTimeStamp.At(index);};
  Int_t   GetMVTimeIdx( Int_t index )        const {return fMVTimeStamp.At(index);};
  Int_t   GetStatusTimeIdx( Int_t index )    const {return fStatusTimeStamp.At(index);};
                       // Returns time stamp of specific DCS variable by index in the array

  Float_t GetTempLeft( Int_t time );
  Float_t GetTempRight( Int_t time );
  Float_t GetHV( Int_t time );
  Float_t GetMV( Int_t time );
  Char_t  GetStatus( Int_t time );
                       // Returns value of specific DCS variable by time stamp

  Bool_t GetOKStat( Int_t time )        { Char_t stat = GetStatus( time ); return (Bool_t)( stat & 1 ); };
  Bool_t GetTempLeftStat( Int_t time )  { Char_t stat = GetStatus( time ); return (Bool_t)( stat & 2 ); };
  Bool_t GetTempRightStat( Int_t time ) { Char_t stat = GetStatus( time ); return (Bool_t)( stat & 4 ); };
                       // Return status of a readout device in given time

  Float_t GetDriftField( Int_t timeStamp );
                       // Returns drift field counted for specific time

  Float_t GetDriftSpeed( Int_t timeStamp ) const;
                       // Returns drift speed counted for specific time. This metod is not dedicated
                       //  for normal usage - it should be used only in cases that the injectors for
                       //  given module fails

  void PrintValues( FILE *output = stdout ) const;
                       // Displays stored DCS varable values or writes it into a text file
 private:

  TArrayF fTempLeft;              //|| "don't split" automatic streamer instruction
                        // Temperature on left side. If there is stored a negative value
                        //  something wrong happend with the temperature chip.
  TArrayI fTempLeftTimeStamp;     //||
                        // Time stamps of the temperatures on left side
  Int_t   fTempLeftMaxPoints;    // Size of the arrays
  Int_t   fTempLeftSetPoints;    // Number of filled array cells (number of set values)


  TArrayF fTempRight;             //||
                        // Temperature on right side. If there is stored a negative value
                        //  something wrong happend with the temperature chip.
  TArrayI fTempRightTimeStamp;    //||
                        // Time stamps of temperatures on right side
  Int_t   fTempRightMaxPoints;    // Size of the arrays
  Int_t   fTempRightSetPoints;    // Number of filled array cells (number of set values)


  TArrayF fHV;                    //||
                        // High voltage on SDD
  TArrayI fHVTimeStamp;           //||
                        // Time stamps of HV
  Int_t   fHVMaxPoints; // Size of the arrays
  Int_t   fHVSetPoints; // Number of filled array cells (number of set values)


  TArrayF fMV;                    //||
                        // Medium voltage on SDD
  TArrayI fMVTimeStamp;           //||
                        // Time stamps of MV
  Int_t   fMVMaxPoints; // Size of the arrays
  Int_t   fMVSetPoints; // Number of filled array cells (number of set values)


  TArrayC fStatus;                //||
                        // Status of temperature and voltage readout
                        //
                        // bit 0 - _OK
                        // bit 1 - _TEMP_L_STATE
                        // bit 2 - _TEMP_R_STATE

  TArrayI fStatusTimeStamp;       //||
                        // Time stamps of MV
  Int_t   fStatusMaxPoints; // Size of the arrays
  Int_t   fStatusSetPoints; // Number of filled array cells (number of set values)


  Int_t   FindIndex( Int_t timeStamp, TArrayI &timeStampArray, Int_t n ) const;
                        // Returns index in the array of time stamps by selected value. Returns -1
                        //  if the time is less than time stamp in the timeArray[0]

  ClassDef(AliITSDCSDataSDD, 2)

}; /*class AliITSDCSDataSDD*/

#endif
