#ifndef ALIEMCALTRIGGERTRUDCSCONFIG_H
#define ALIEMCALTRIGGERTRUDCSCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"
#include <iosfwd>

/**
 * @class AliEMCALTriggerTRUDCSConfig
 * \brief TRU configuration in OCDB
 * @ingroup EMCALbase
 * @author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
 * @author: Jiri Kral, JYU
 * 
 * This object is the OCDB object containing the 
 * TRU configuration for a single TRU provided
 * by the Preprocessor. Configuration values are
 * - Peak finder
 * - L0 algorithm version
 * - 2x2 threshold
 * - 4x4 threshold
 * - Masking (96 channels)
 * - Rollback
 * - Firmware version
 * Note that the firmare version was not set
 * for run2 data from the preprocessor, instead
 * always the default firmware version is provided.
 * 
 * The OCDB object provides also streamer versions
 * to ASCII file and to JSON strings.
 */
class AliEMCALTriggerTRUDCSConfig : public TObject 
{

public:

	/**
	 * @brief Default constructor
	 */
	AliEMCALTriggerTRUDCSConfig();

	/**
	 * @brief Destructor
	 */
	virtual ~AliEMCALTriggerTRUDCSConfig() {}

	/**
 	 * @brief equalty operator
   * 
	 * Checking if the two TRU DCS configurations are equal. For equalty
	 * all settings must be the same.
	 */ 
	bool operator==(const AliEMCALTriggerTRUDCSConfig &other) const;

	/**
	 * @brief Streaming operator
	 * 
	 * Printing all settings of the given TRU on the output stream
	 */
	friend std::ostream &operator<<(std::ostream &stream, const AliEMCALTriggerTRUDCSConfig &other);

	/**
	 * @brief Serialize object to JSON format
	 * 
	 * @return JSON-serialized TRU DCS config object 
	 */
	std::string ToJSON() const;

	/**
	 * @brief Set peak finder algorithm
	 * @param pf Peak finder algorithm
	 * 
	 * The peak finder algorithm consists of two parts:
	 * - Pattern: bits 8-14
	 * - Mask:  bits 0-6
	 * 
	 * Both are used in the L0 algorithm to find the number of samples
	 * before / after the peak
	 */
	void    SetSELPF(  UInt_t pf)              { fSELPF  = pf;        }

	/**
	 * @brief Set the L0 algorithm version
	 * @param la L0 algorithm version
	 * 
	 * Following versions are supported:
	 * - 0x0001: L0v1
	 * - 0x1: L0v0
	 */
	void    SetL0SEL(  UInt_t la)              { fL0SEL  = la;        }

	/**
	 * @brief Set L0 cosmic threshold (2x2 patches)
	 * @param lc L0 cosmic threshold
	 */
	void    SetL0COSM( UInt_t lc)              { fL0COSM = lc;        }

	/**
	 * @brief Set the L0 threshold (4x4 patches)
	 * @param lg L0 threshold
	 */
	void    SetGTHRL0( UInt_t lg)              { fGTHRL0 = lg;        }

	/**
	 * @brief Set the L0 mask object for a given TRU
	 * @param msk 16 bit mask
	 * @param pos Position in the array
	 * 
	 * The mask consists of 96 bits split to 6 array entries of each 16 bit. 
	 * Each bit relates to a given FastOR within the TRU, where the online 
	 * masking is applied. In case the bit for a given FastOR is set to 1
	 * the FastOR is enabled, otherwise it is disabled.
	 */
	void    SetMaskReg(UInt_t msk, Int_t pos)  { fMaskReg[pos] = msk; }

	/**
	 * @brief Set the TRU rollback
	 * @param rb TRU rollback
	 * 
	 * TRU rollback is defined as the number of samples in the 
	 * time sliding window with respect to the L0 alignment
	 */
	void    SetRLBKSTU(UInt_t rb)              { fRLBKSTU = rb;       }

	/**
	 * @brief Set the firmware version
	 * @param fw Firmware version
	 */
	void    SetFw(     UInt_t fw)              { fFw = fw;            }
			
	/**
	 * @brief Get the peak finder algorithm
	 * @return peak finder algorithm (encoding pattern and mask) 
	 * 
	 * See SetSELPF for explanation of the peak finder algorithm
	 */
	UInt_t  GetSELPF()                   const { return fSELPF;       }

	/**
	 * @brief Get the L0 algorithm version
	 * @return L0 algorithm version
	 * 
	 * Following versions are supported:
	 * - 0x0001: L0v1
	 * - 0x1: L0v0
	 */
	UInt_t  GetL0SEL()                   const { return fL0SEL;       }

	/**
	 * @brief Get L0 cosmic threshold (2x2 patches)
	 * @return L0 cosmic threshold 
	 */
	UInt_t  GetL0COSM()                  const { return fL0COSM;      }

	/**
	 * @brief 
	 * 
	 * @return UInt_t 
	 */
	UInt_t  GetGTHRL0()                  const { return fGTHRL0;      }

	/**
	 * @brief Get FastOR mask index for a given position
	 * @param pos Position in the array
	 * @return Mask bits for the given position 
	 * 
	 * See SetMaskReg for the defintion of the L0 mask
	 */
	UInt_t  GetMaskReg(Int_t pos)        const { return fMaskReg[pos];}

	/**
	 * @brief Get the TRU rollback
	 * @return TRU rollback
	 * 
	 * Set 
	 */
	UInt_t  GetRLBKSTU()                 const { return fRLBKSTU;     }

	/**
	 * @brief Get the firmware version
	 * @return Firmware version 
	 */
	UInt_t  GetFw()                      const { return fFw;          }
	
	/**
	 * @brief Get size of the L0 patch
	 * @return Size of the L0 patch (in number of FastORs per direction) 
	 */
	Int_t   GetSegmentation();
	
protected:

	AliEMCALTriggerTRUDCSConfig           (const AliEMCALTriggerTRUDCSConfig &cd);
	AliEMCALTriggerTRUDCSConfig &operator=(const AliEMCALTriggerTRUDCSConfig &cd);

private:
	
  UInt_t   fSELPF;                         ///< PeakFinder setup
  UInt_t   fL0SEL;                         ///< L0 Algo selection
  UInt_t   fL0COSM;                        ///< 2x2 threshold
  UInt_t   fGTHRL0;                        ///< 4x4 threshold
  UInt_t   fMaskReg[6];                    ///< 6*16 = 96 mask bits per TRU
  UInt_t   fRLBKSTU;                       ///< TRU circular buffer rollback
  UInt_t   fFw;                            ///< TRU firmware version
	
  ClassDef(AliEMCALTriggerTRUDCSConfig,4) ;
};

#endif // ALIEMCALTRIGGERTRUDCSCONFIG_H
