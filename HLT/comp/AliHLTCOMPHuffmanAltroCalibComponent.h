//-*- Mode: C++ -*-
#ifndef ALIHLTCOMPHUFFMANALTROCALIBCOMPONENT_H
#define ALIHLTCOMPHUFFMANALTROCALIBCOMPONENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTCOMPHuffmanAltroCalibComponent.h
    @author Jochen Thaeder
    @author extended by Jenny Wagner
    @date   20-11-2007
    @brief  A calibration component for the Huffman code creator 
*/

#include "AliHLTCalibrationProcessor.h"

class AliHLTCOMPHuffmanAltro;
class AliHLTCOMPHuffmanData;

/**
 * @class AliHLTCOMPHuffmanAltroCalibComponent
 * Component ID: \b HuffmanAltroCalibComponent <br>
 * Library: \b libAliHLTComp
 * This class is the calibration component for the AliTPCCalibHuffmanAltro class 
 * used for calibration of the Huffman code table (which is created here). 
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -origin <i> detector <\i> <br>
 *  set origin of data for code creation to specify output table (parameter transient)
 * \li -runnumber <i> decimal number <\i> <br>
 *  set runnumber to specify output table (parameter transient)
 * \li -dataspec <i> 0xYYXXaabb <\i> <br>
 *  set usual HLT dataspec (last slice, first slice, last patch, first patch)_Hexadezimal to specify output table
 * \li -trailerwords <i> decimal number <\i> <br>
 *  set number of trailerwords of incoming data (ranging from 1 to 3)
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -tablepath <i> path to Huffman code table <\i> <br>
 *  set path to out put Huffman code table as root file, if no path is given, output path is set to current path (parameter transient)
 *
 * It inherits from the @ref AliHLTCalibrationProcessor and uses the high-level 
 * interface. The output is the class @ref HuffmanData as a TObject which is written to the data base
 *
 * @ingroup alihlt_comp
 */
class AliHLTCOMPHuffmanAltroCalibComponent : public AliHLTCalibrationProcessor
    {
    public:
      /** constructor */
      AliHLTCOMPHuffmanAltroCalibComponent();
      
      /** destructor */
      virtual ~AliHLTCOMPHuffmanAltroCalibComponent();
      
      // Public functions to implement AliHLTComponent's interface.
      // These functions are required for the registration process

      const char* GetComponentID();
      void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
      AliHLTComponentDataType GetOutputDataType();
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
      AliHLTComponent* Spawn();

    protected:
      
      // Protected functions to implement AliHLTComponent's interface.
      // These functions provide initialization as well as the actual processing
      // capabilities of the component. 
      
      /** Initialize the calibration component. */
      Int_t InitCalibration();

      /** Scan commandline arguments of the calibration component. */
      Int_t ScanArgument( Int_t argc, const char** argv );

      /** DeInitialize the calibration component. */
      Int_t DeinitCalibration();

      /** Process the data in the calibration component. */
      Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

      /** Ship the data to the FXS at end of run or eventmodulo. */
      Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

    private:

      /** copy constructor prohibited */
      AliHLTCOMPHuffmanAltroCalibComponent(const AliHLTCOMPHuffmanAltroCalibComponent&);

      /** assignment operator prohibited */
      AliHLTCOMPHuffmanAltroCalibComponent& operator=(const AliHLTCOMPHuffmanAltroCalibComponent&);

      /** Huffman compressor class */
      AliHLTCOMPHuffmanAltro * fHuffmanCompressor; //! instance of Huffman Compressor in this component
      
      /** pointer to output Huffman code table and occurrence table (togehter in this class) */
      AliHLTCOMPHuffmanData * fHuffmanData;            //! instance of output (Huffman Data, containing code table)

      /** The Specification for this component */
      /** explicit specification of the origin of the data (transient) */
      TString fOrigin;                      // input line argument to determine origin for Huffman table
                                           // -> no input --> default taken from incoming data
      /** explicit specification of the run number */
      AliHLTUInt64_t fRunNumber;           // input line argument to determine run number for Huffman table

      /** specifications of the data */
      AliHLTUInt64_t fSpecification;      // see above

      /** explicit path to Huffman code table which will be put out */
      TString fTablePath;                 // input line argument to determine path for Huffman table
                                          // -> no input --> default set to current path name
      //AliHLTUInt8_t fSlice;             // slice 
      //AliHLTUInt8_t fPatch;            // patch   

      /** number of NRCU trailer words of input data */
      Int_t fNRCUTrailerWords;          // see above

      ClassDef(AliHLTCOMPHuffmanAltroCalibComponent, 0)

    };
#endif
