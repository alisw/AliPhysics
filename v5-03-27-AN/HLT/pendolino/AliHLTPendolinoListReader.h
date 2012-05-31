//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPENDOLINOLISTREADER_H
#define ALIHLTPENDOLINOLISTREADER_H

/************************************************************************
**
**
** This file is property of and copyright by the Department of Physics
** Institute for Physic and Technology, University of Bergen,
** Bergen, Norway, 2007
** This file has been written by Sebastian Bablok,
** sebastian.bablok@ift.uib.no
**
** Important: This file is provided without any warranty, including
** fitness for any particular purpose.
**
**
*************************************************************************/

//  @file   AliHLTPendolinoListReader.h
//  @author Sepastian Bablok
//  @date   
//  @brief  Helper class for pendolino list handling
//  @note   maintained by matthias.richter@cern.ch

//#include <vector>
//#include <map>
#include <string>

#include <TString.h>
#include <TObject.h>
#include <TMap.h>


//namespace alice { namespace hlt { namespace pendolino {

//#define MAX_LINE_LENGTH 256

/**
 * This class reads the file containing the list of calibration objects the
 * Pendolino shall fetch.
 *
 * @author Sebastian Bablok
 *
 * @date 2007-02-23
 */
class AliHLTPendolinoListReader : public TObject {
	public:

		/**
		 * Static that defines the maximum line length inside the list file
		 */
		static int kMAX_LINE_LENGTH;
			
		/**
		 * Constructor for the AliHLTPendolinoListReader.
		 */
		AliHLTPendolinoListReader();

		/**
		 * Constaructor for AliHLTPendolinoListReader, which also loads the 
		 * list from file
		 *
		 * @param filename of the list file to read given as char*
		 */
		AliHLTPendolinoListReader(const char* filename);

		/**
		 * Destructor for the AliHLTPendolinoListReader.
		 */
		virtual ~AliHLTPendolinoListReader();

		/**
		 * Reads the list from the file given by filename
		 * and stores the calibration object names in this class/object.
		 *
		 * @param filename of the list file to read given as char*
		 *
		 * @return false, if file does not exist or file is empty, else true
		 */
		bool ReadListFromFile(const char* filename);

		/**
		 * Reads the list from the file given by filename
		 * and stores the calibration object names in this class/object.
         *
         * @param filename of the list file to read given as string
		 *
		 * @return false, if file does not exist or file is empty, else true
		 */
		bool ReadListFromFile(const std::string filename);

		/**
		 * Function to get the list of calibration object names.
		 *
		 * @return pointer to a TMap containing the the lists grouped by 
		 * 			detector. The Alias names are already included as the key
		 * 			in a map, so they can be directly used for contacting
		 * 			the AliDCSClient. Structure of the returned TMap:
		 * <pre>
		 * TMap( detector name [key], lists of alias names [value] )
		 *       {TObjString}         {TSeqCollection( Alias names )}
		 *                                             {TObjString}       
		 * e.g.:
		 * TMap( TRD [key], TRD alias list [value]
		 *                  (    [Alias names]             
		 *                    -> TRD_low_volt_xyz, 
		 *                    -> TRD_low_volt_123, 
		 *                    -> TRD_high_volt,    
		 *                    -> TRD_temp_1234,    )
		 *       TPC [key], TPC alias list [value]
		 *                  (    [Alias names]             
		 *                    -> TPC_low_volt_abc, 
		 *                    -> TPC_temp_1234,    
		 *                    -> Pressure_sens_01, )
		 *       ...)
		 * </pre>
		 */
		TMap* GetCalibObjList();

		/**
		 * Function to retrieve the latest run number
		 *
		 * @param path path to the directory containing the lastRunNumber file,
		 *			NOTE: the file name must not be included to the path
		 *
		 * @return latest run number
		 */
		static int RetrieveLastRunNumber(std::string path);

		/**
         * Function to retrieve the latest run number
         *
         * @param path path to the directory containing the lastRunNumber file,
         *          NOTE: the file name must not be included to the path
         *
         * @return latest run number
         */
		static int RetrieveLastRunNumber(char* path);

		/**
		 * Function to print out the content of the Map containing the lists.
		 */
		void Print() const;

// check here !!!
		/**
		 * Retrieves and returns the value of a given property name.
		 *
		 * @param [in] propertyName the name of the requested property (as char)
		 * @param [out] propertyValue pointer to the string representing the
		 *			value of the requested property
		 *
		 * @return true, if requested property has been found in the class, else
		 * 			false (false also when PropertyReader is not valid or
		 *			propertyName is NULL)
		 */
//		bool retrievePropertyValue(const char* propertyName,
//					std::string** propertyValue);

		/**
		 * Retrieves and returns the value of a given property name.
		 *
		 * @param [in] propertyName pointer to the string containing the name of
		 *			the requested property
		 * @param [out] propertyValue pointer to the string representing the
		 *			value of the requested property
		 *
		 * @return true, if requested property has been found in the class, else
		 * 			false (false also when PropertyReader is not valid or
		 *			propertyName is NULL)
		 */
//		bool retrievePropertyValue(const std::string* propertyName,
//			std::string** propertyValue);


	private:

        /**
         * Flag, indicating if list has already been read from file.
         */
        bool fValid;  //! see above

        /**
         * Vector containing the list of calibration objects
         */
//		std::vector<std::string> mCalibObjList;



		/**
		 * Map containing all PropertyName PropertyValue pairs read from the
		 * property file.
		 */
		TMap fCalibObjList;//! mProperties;

        /**
         * AliRoot required stuff
         */
        ClassDef(AliHLTPendolinoListReader, 5);
}; //end of class


inline TMap* AliHLTPendolinoListReader::GetCalibObjList() {
	return &fCalibObjList;
}


#endif // ALIHLTPENDOLINOLISTREADER_H

