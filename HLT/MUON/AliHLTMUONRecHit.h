#ifndef ALIHLTMUONRECHIT_H
#define ALIHLTMUONRECHIT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONRecHit.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Declaration of a reconstructed hit ROOT object to store 3D hit coordinates.
///

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include <ostream>

/**
 * A 3D hit object used to store hits reconstructed on the tracking chambers by
 * the dHLT. These objects store information translated into ROOT format from
 * dHLT raw data. Reconstructed hit values of (0, 0, 0) indicate an invalid or
 * nil hit.
 * This class is mainly for testing or as a helper object for dHLT specific analysis,
 * since it is sometimes easier to store and handle ROOT objects.
 */
class AliHLTMUONRecHit : public TObject
{
	/**
	 * Stream operator for usage with std::ostream classes.
	 * Allows usage such as:
	 *   AliHLTMUONRecHit h; std::cout << h;
	 */
	friend std::ostream& operator << (std::ostream& stream, const AliHLTMUONRecHit& hit);

public:

	/**
	 * The AliChannel class stores extra debugging information about the channels
	 * and raw data words that were considered during reconstruction of a hit
	 * by the dHLT hit reconstructor component.
	 */
	class AliChannel : public TObject
	{
		/**
		 * Stream operator for usage with std::ostream classes.
		 */
		friend std::ostream& operator << (std::ostream& stream, const AliChannel& c);
	
	public:
		
		/**
		 * Constructor.
		 * \param manu  The MANU ID of the channel as found in the raw data word.
		 * \param channel  The MANU channel ID as found in the raw data word.
		 * \param signal  The ADC signal value as found in the raw data word.
		 * \param rawDataWord  The actual raw data word.
		 */
		AliChannel(
				Short_t manu = -1,
				Short_t channel = -1,
				Short_t signal = -1,
				UInt_t rawDataWord = 0
			) :
			TObject(),
			fManu(manu), fAddress(channel), fSignal(signal),
			fRawDataWord(rawDataWord)
		{}
		
		/**
		 * Default destructor.
		 */
		virtual ~AliChannel() {}
		
		/**
		 * Returns the MANU address.
		 */
		Short_t Manu() const { return fManu; }
		
		/**
		 * Returns the channel address of the MANU.
		 */
		Short_t Address() const { return fAddress; }
		
		/**
		 * Returns the ADC signal measured on the channel.
		 */
		Short_t Signal() const { return fSignal; }
		
		/**
		 * Returns the raw data word as found in the tracking DDL payload.
		 */
		UInt_t RawDataWord() const { return fRawDataWord; }
		
		virtual void Print(Option_t* option = NULL) const;
	
		// Methods inherited from TObject
		virtual Bool_t IsSortable() const { return kTRUE; }
		Int_t Compare(const TObject* obj) const;

		// Implement comparison operators.
		bool operator == (const AliChannel& c) const
		{
			return fManu == c.fManu and fAddress == c.fAddress 
				and fSignal == c.fSignal
				and fRawDataWord == c.fRawDataWord;
		}

		bool operator != (const AliChannel& c) const
		{
			return not this->operator == (c);
		}
	
	private:
	
		Short_t fManu;       ///< The MANU address on the electronics.
		Short_t fAddress;    ///< The channel address on the electronics.
		Short_t fSignal;     ///< ADC value of signal.
		UInt_t fRawDataWord; ///< The raw data word as found in the DDL stream.
		
		ClassDef(AliHLTMUONRecHit::AliChannel, 3); // A MANU channel forming part of a cluster that was considered during hit reconstruction in dHLT.
	};

	/**
	 * Construct a new AliHLTMUONRecHit object with coordinate (x, y, z).
	 * @param x           X coordinate of hit
	 * @param y           Y coordinate of hit
	 * @param z           Z coordinate of hit
	 * @param sourceDDL   The DDL from which this hit originates.
	 * @param detectorId  The ID number for the AliRoot detector element on
	 *                    which this hit resides.
	 * @param clusterId   The cluster ID number assigned to the hit's cluster.
	 * @param nChExp      The expected number of channels that form the cluster.
	 * @param charge      The total charge of the cluster.
	 */
	AliHLTMUONRecHit(
			Float_t x = 0,
			Float_t y = 0,
			Float_t z = 0,
			Int_t sourceDDL = -1,
			Int_t detElemId = -1,
			Int_t clusterId = -1,
			Int_t nChExp = 0,
			Float_t charge = 0
		) :
		TObject(), fCoordinate(x, y, z), fSourceDDL(sourceDDL),
		fDetElemId(detElemId), fClusterId(clusterId), fNchExp(nChExp),
		fChannels("AliHLTMUONRecHit::AliChannel", 6), fCharge(charge)
	{}
	
	/**
	 * Default destructor.
	 */
	virtual ~AliHLTMUONRecHit() {}

	/**
	 * Returns the 3D hit coordinate.
	 */
	const TVector3& Coordinate() const { return fCoordinate; }

	/**
	 * Returns the X coordinate of the reconstructed hit in centimetres.
	 */
	Double_t X() const { return fCoordinate.X(); }

	/**
	 * Returns the Y coordinate of the reconstructed hit in centimetres.
	 */
	Double_t Y() const { return fCoordinate.Y(); }

	/**
	 * Returns the Z coordinate of the reconstructed hit in centimetres.
	 */
	Double_t Z() const { return fCoordinate.Z(); }
	
	/**
	 * Returns the source DDL from which this hit originates.
	 * -1 is returned if this was not set.
	 */
	Int_t SourceDDL() const { return fSourceDDL; }
	
	/**
	 * Returns the detector element ID on which this reconstructed hit resides.
	 * -1 is returned if this was not set.
	 */
	Int_t DetElemId() const { return fDetElemId; }
	
	/**
	 * Returns the chamber number of this hit in the range [1..14].
	 * If -1 is returned then the chamber number is not known because the
	 * detector element ID was not set or is invalid.
	 * @param warn  Indicates if any warning should be printed in case of problems.
	 */
	Int_t Chamber(bool warn = true) const;
	
	/**
	 * Returns the ID number given to the hit's cluster.
	 */
	Int_t ClusterId() const { return fClusterId; }
	
	/**
	 * Returns the expected number of channels that are to be added to this hit.
	 * If the number of calls to AddChannel does not correspond to this value
	 * then we lost some debugging information along the way.
	 */
	UInt_t ExpectedNchannels() const { return fNchExp; }
	
	/**
	 * Returns the total charge of the cluster.
	 */
	Float_t TotalCharge() const { return fCharge; }
	
	/**
	 * Sets the extra debugging information for this hit.
	 * @param detElemId  The detector element ID.
	 * @param clusterId  Cluster ID of the hit's cluster.
	 * @param nChExp     Number of expected channels forming the cluster.
	 * @param charge     The total charge of the cluster
	 * @param sourceDDL  The source DDL of this hit.
	 */
	void SetDebugInfo(
			Int_t detElemId, Int_t clusterId, UInt_t nChExp,
			Float_t charge, Int_t sourceDDL = -1
		);
	
	/**
	 * Sets the hit coordinate.
	 */
	void SetHit(Float_t x, Float_t y, Float_t z)
	{
		fCoordinate.SetXYZ(x, y, z);
	}
	
	/**
	 * Returns the number of channels associated with this hit.
	 */
	Int_t Nchannels() const { return fChannels.GetEntriesFast(); }
	
	/**
	 * Returns the i'th channel associated with this hit.
	 * @param i  Should be a number in the range [0..n), where n = Nchannels().
	 */
	const AliChannel* GetChannel(Int_t i) const
	{
		return static_cast<const AliChannel*>(fChannels[i]);
	}
	
	/**
	 * Adds a new channel to this hit if it is on a tracking chamber.
	 * @param manu    The MANU number
	 * @param channel The MANU channel address.
	 * @param signal  The ADC signal value measured on the channel.
	 * @param rawDataWord This is the raw data word as read from the DDL.
	 */
	void AddChannel(
			Short_t manu, Short_t channel, Short_t signal,
			UInt_t rawDataWord
		);

	/**
	 * Prints the details of the reconstructed hit.
	 * @param option  A case sensitive string that can contain one of the
	 *     following strings:
	 *       "compact" - Prints just the coordinates of the hit in a terse format.
	 *       "detail" - Prints the coordinates and detector element ID.
	 *       "all" - Prints all known information about this hit including
	 *               channel information forming the cluster that was reconstructed.
	 *     If the string contains an empty option or NULL then the default is
	 *     to print compactly.
	 */
	virtual void Print(Option_t* option = NULL) const;
	
	// Methods inherited from TObject
	virtual Bool_t IsSortable() const { return kTRUE; }
	Int_t Compare(const TObject* obj) const;

	// Implement comparison operators.
	bool operator == (const AliHLTMUONRecHit& hit) const
	{
		return X() == hit.X() and Y() == hit.Y() and Z() == hit.Z();
	}

	bool operator != (const AliHLTMUONRecHit& hit) const
	{
		return not this->operator == (hit);
	}

private:

	TVector3 fCoordinate; ///< The 3D coordinate of the hit in AliRoot global coordinates (cm).
	
	// The following is debugging information and may not be filled if the
	// dHLT components were not set to produce this information.
	Int_t fSourceDDL;  ///< The DDL from which this hit originates.
	Int_t fDetElemId;  ///< Detector element ID number.
	Int_t fClusterId;  ///< The cluster ID number used to relate all the channels to each other.
	UInt_t fNchExp;    ///< The number of channels that were supposed to be found.
	TClonesArray fChannels; ///< The channels forming part of the cluster from which this hit was reconstructed.
	Float_t fCharge;   ///< The total charge of the cluster.

	ClassDef(AliHLTMUONRecHit, 4); // A reconstructed hit translated from dHLT raw data into ROOT format.
};

#endif // ALIHLTMUONRECHIT_H
