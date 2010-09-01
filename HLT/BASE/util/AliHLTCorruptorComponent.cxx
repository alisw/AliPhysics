// $Id: $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTCorruptorComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   5 Aug 2010
/// @brief  Implementation of the AliHLTCorruptorComponent class.
///
/// The AliHLTCorruptorComponent is used to generate corruption of various kind
/// in the input data blocks it sees. The blocks will be copied to the output
/// but with various bit flips and garbage inserted.
/// This kind of component is used for testing purposes.

#include "AliHLTCorruptorComponent.h"
#include "TRandom3.h"
#include <cstdlib>
#include <cstring>
#include <cerrno>

ClassImp(AliHLTCorruptorComponent)


AliHLTCorruptorComponent::AliHLTCorruptorComponent() :
	AliHLTProcessor(),
	fBufferMultiplier(2.),
	fBlockIdList(),
	fCorruptionInfoList(),
	fPatterns()
{
	// Default constructor.
}


AliHLTCorruptorComponent::~AliHLTCorruptorComponent()
{
	// Default destructor.
}


const char* AliHLTCorruptorComponent::GetComponentID()
{
	// Returns the component ID.
	return "CorruptorComponent";
}


void AliHLTCorruptorComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	// Returns the list of input data types that are handled.
	list.push_back(kAliHLTAnyDataType);
}


AliHLTComponentDataType AliHLTCorruptorComponent::GetOutputDataType()
{
	// Returns kAliHLTMultipleDataType.
	return kAliHLTMultipleDataType;
}


int AliHLTCorruptorComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	// Returns the list of output data blocks handled.
	list.push_back(kAliHLTAnyDataType);
	return int(list.size());
}


void AliHLTCorruptorComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
	// Returns the buffer size requirements.
	
	unsigned long total = 0;
	for (size_t i = 0; i < fCorruptionInfoList.size(); ++i)
	{
		AliCorruptionInfo& info = fCorruptionInfoList[i];
		if (info.fType == kInsert)
		{
			total += (info.fLastPattern - info.fFirstPattern) / sizeof(AliPattern) * sizeof(AliHLTUInt64_t);
		}
		else if (info.fType == kInsertRandom)
		{
			total += info.fMaxBits / 8 + 1;
		}
	}
	
	constBase = total;
	inputMultiplier = fBufferMultiplier;
}


AliHLTComponent* AliHLTCorruptorComponent::Spawn()
{
	// Creates a new instance of the component.
	return new AliHLTCorruptorComponent;
}


bool AliHLTCorruptorComponent::ConvertToPositiveInt(const char* value, AliHLTUInt64_t& num, bool printErrors) const
{
	// Converts a string value to a positive integer.
	
	char* err = NULL;
	errno = 0;
	unsigned long long tmpnum = strtoull(value, &err, 0);
	if (err == NULL or *err != '\0')
	{
		if (printErrors) HLTError("Cannot convert '%s' to a positive integer.", value);
		return false;
	}
	if (errno == ERANGE)
	{
		if (printErrors) HLTError("The specified value '%s' is out of range.", value);
		return false;
	}
	num = tmpnum;
	return true;
}


bool AliHLTCorruptorComponent::ConvertToBitPosition(
		const char* value, AliHLTUInt64_t& pos, bool& relative, bool printErrors
	) const
{
	// Converts a string value to a bit position.
	
	if (strcmp(value, "min") == 0)
	{
		// Special case where the value is the special symbol "min".
		pos = 0;
		relative = false;
		return true;
	}
	if (strcmp(value, "max") == 0)
	{
		// Special case where the value is the special symbol "max".
		pos = 0xFFFFFFFFFFFFFFFFull;
		relative = false;
		return true;
	}
	
	int valuelength = strlen(value);
	const char* valstr = value;
	relative = false;
	// Check if the position is relative.
	if (valuelength >= 3 and value[0] == 'm' and value[1] == 'i' and value[2] == 'n')
	{
		if (valuelength >= 4 and value[3] == '-')
		{
			if (printErrors)
			{
				HLTError("Cannot use 'min-' in option '%s' since that would be outside the buffer."
					" Should use 'min+' instead.",
					value
				);
				return false;
			}
		}
		else if (valuelength >= 4 and value[3] == '+')
		{
			// Value starting with min+ does not need anything special because
			// min is always == 0. Just adjust the value string pointer to
			// skip this part of the string.
			valstr = value + 4;
		}
	}
	else if (valuelength >= 3 and value[0] == 'm' and value[1] == 'a' and value[2] == 'x')
	{
		if (valuelength >= 4 and value[3] == '+')
		{
			if (printErrors)
			{
				HLTError("Cannot use 'max+' in option '%s' since that would be outside the buffer."
					" Should use 'max-' instead.",
					value
				);
				return false;
			}
		}
		else if (valuelength >= 4 and value[3] == '-')
		{
			// Need to mark the result as a relative value and skip the
			// first part of the string.
			relative = true;
			valstr = value + 4;
		}
	}
	
	// Find the location of ':' if any, to mark the start of the byte and bit strings.
	int bytestrlen = strlen(valstr);
	const char* bitstr = valstr + bytestrlen;
	for (int i = 0; valstr[i] != 0x0; ++i)
	{
		if (valstr[i] == ':')
		{
			bytestrlen = i;
			bitstr = valstr + (i+1);
			break;
		}
	}
	TString bytestr(valstr, bytestrlen);
	
	char* err = NULL;
	errno = 0;
	unsigned long long bytenum = strtoull(bytestr.Data(), &err, 0);
	if (err == NULL or *err != '\0')
	{
		if (printErrors)
		{
			HLTError("Cannot convert '%s' given in option '%s' to a positive integer.",
				bytestr.Data(), value
			);
		}
		return false;
	}
	if (errno == ERANGE or bytenum > 0x1FFFFFFFFFFFFFFFull)
	{
		if (printErrors)
		{
			HLTError("The specified byte number '%s' given in option '%s' is out of range."
				" The value should be in the range [0..%lld].",
				bytestr.Data(), value, 0x1FFFFFFFFFFFFFFFull
			);
		}
		return false;
	}
	err = NULL;
	unsigned long bitnum = strtoul(bitstr, &err, 0);
	if (err == NULL or *err != '\0')
	{
		if (printErrors)
		{
			HLTError("Cannot convert '%s' given in option '%s' to a positive integer.",
				bitstr, value
			);
		}
		return false;
	}
	if (errno == ERANGE or bitnum > 7)
	{
		if (printErrors)
		{
			HLTError("The specified bit number '%s' given in option '%s' is out of range."
				" The value should be in the range [0..7].",
				bitstr, value
			);
		}
		return false;
	}
	pos = (bytenum << 3) | (bitnum & 0x7);
	return true;
}


bool AliHLTCorruptorComponent::ConvertToPercentage(const char* value, double& num, bool printErrors) const
{
	// Converts a string value as a percentage to a floating point number in the range [0..1].
	
	TString str = value;
	bool percent = false;
	if (str.Length() > 0 and str[str.Length()-1] == '%')
	{
		percent = true;
		str[str.Length()-1] = '\0';
	}
	char* err = NULL;
	errno = 0;
	double tmpnum = strtod(str.Data(), &err);
	if (err == NULL or *err != '\0')
	{
		if (printErrors) HLTError("Cannot convert '%s' to a valid floating point value.", value);
		return false;
	}
	if (errno == ERANGE or tmpnum < 0. or (tmpnum > 100. and percent) or (tmpnum > 1. and not percent))
	{
		if (printErrors)
		{
			int range = percent ? 100 : 1;
			HLTError("The specified value '%s' is outside the range [0..%d]", value, range);
		}
		return false;
	}
	num = percent ? tmpnum / 100. : tmpnum;
	return true;
}


bool AliHLTCorruptorComponent::ConvertToPattern(const char* value, AliPattern& pattern, bool printErrors) const
{
	// Coverts a string value into a pattern structure.
	
	if (strcmp(value, "removed") == 0)
	{
		// Special case where the pattern is the special symbol to used remove bits.
		pattern.fPattern = 0;
		pattern.fWidth = 0;
		pattern.fUseRemoved = true;
		return true;
	}
	
	// Find the location of '/' if any, to mark the start of the pattern and width fields.
	int patstrlen = strlen(value);
	const char* widthstr = value + patstrlen;
	for (int i = 0; value[i] != 0x0; ++i)
	{
		if (value[i] == '/')
		{
			patstrlen = i;
			widthstr = value + (i+1);
			break;
		}
	}
	TString patstr(value, patstrlen);
	
	char* err = NULL;
	errno = 0;
	AliHLTUInt64_t patnum = AliHLTUInt64_t( strtoull(patstr.Data(), &err, 0) );
	if (err == NULL or *err != '\0')
	{
		if (printErrors)
		{
			HLTError("Cannot convert '%s' given in pattern '%s' to a positive integer.",
				patstr.Data(), value
			);
		}
		return false;
	}
	if (errno == ERANGE)
	{
		if (printErrors)
		{
			HLTError("The specified pattern '%s' is out of range."
				"Cannot be bigger than a 64 bit integer.",
				 value
			);
		}
		return false;
	}
	unsigned long widthnum = 1;
	if (*widthstr != '\0')  // check if widthstr string is not empty
	{
		err = NULL;
		widthnum = strtoul(widthstr, &err, 0);
		if (err == NULL or *err != '\0')
		{
			if (printErrors)
			{
				HLTError("Cannot convert the width '%s' given in pattern '%s'"
					" to a positive integer.",
					widthstr, value
				);
			}
			return false;
		}
		if (errno == ERANGE or widthnum < 1 or 64 < widthnum)
		{
			if (printErrors)
			{
				HLTError("The specified width '%s' given in pattern '%s' is out of range."
					" The value should be in the range [1..64].",
					widthstr, value
				);
			}
			return false;
		}
	}
	else
	{
		// Find the width of the pattern if not explicitly given.
		// This is done by finding the most significant set bit.
		for (int i = 63; i >= 0; --i)
		{
			if (((patnum >> i) & 0x1) == 0x1)
			{
				widthnum = i+1;
				break;
			}
		}
	}
	pattern.fPattern = patnum;
	pattern.fWidth = AliHLTUInt8_t(widthnum);
	pattern.fUseRemoved = false;
	return true;
}


bool AliHLTCorruptorComponent::AddBlockTypeId(const char* type)
{
	// Sets the block type ID to filter on.
	
	if (strlen(type) > (unsigned)kAliHLTComponentDataTypefIDsize)
	{
		HLTError("The specified block type '%s' must not be longer than %d characters.",
			type, kAliHLTComponentDataTypefIDsize
		);
		return false;
	}
	if (fBlockIdList.size() > 0 and not fBlockIdList[fBlockIdList.size()-1].fTypeSet)
	{
		fBlockIdList[fBlockIdList.size()-1].fType =
			AliHLTComponentDataTypeInitializer(
					type,
					fBlockIdList[fBlockIdList.size()-1].fType.fOrigin
				);
		fBlockIdList[fBlockIdList.size()-1].fTypeSet = true;
	}
	else
	{
		AliBlockId id;
		id.fType = AliHLTComponentDataTypeInitializer(type, kAliHLTDataOriginAny);
		id.fSpec = kAliHLTVoidDataSpec;
		id.fTypeSet = true;
		id.fOriginSet = false;
		id.fSpecSet = false;
		fBlockIdList.push_back(id);
	}
	return true;
}


bool AliHLTCorruptorComponent::AddBlockOrigin(const char* origin)
{
	// Sets the block origin to filter on.
	
	if (strlen(origin) > (unsigned)kAliHLTComponentDataTypefOriginSize)
	{
		HLTError("The specified origin '%s' must not be longer than %d characters.",
			origin, kAliHLTComponentDataTypefOriginSize
		);
		return false;
	}
	if (fBlockIdList.size() > 0 and not fBlockIdList[fBlockIdList.size()-1].fOriginSet)
	{
		fBlockIdList[fBlockIdList.size()-1].fType =
			AliHLTComponentDataTypeInitializer(
					fBlockIdList[fBlockIdList.size()-1].fType.fID,
					origin
				);
		fBlockIdList[fBlockIdList.size()-1].fOriginSet = true;
	}
	else
	{
		AliBlockId id;
		id.fType = AliHLTComponentDataTypeInitializer(kAliHLTAnyDataTypeID, origin);
		id.fSpec = kAliHLTVoidDataSpec;
		id.fTypeSet = false;
		id.fOriginSet = true;
		id.fSpecSet = false;
		fBlockIdList.push_back(id);
	}
	return true;
}


bool AliHLTCorruptorComponent::AddBlockSpec(const char* spec)
{
	// Sets the block specification to filter on.
	
	char* err = NULL;
	errno = 0;
	unsigned long long num = strtoull(spec, &err, 0);
	if (err == NULL or *err != '\0')
	{
		HLTError("Cannot convert '%s' to a specification number.", spec);
		return false;
	}
	if (num > 0xFFFFFFFF or errno == ERANGE)
	{
		HLTError("The specification number is not an unsigned 32-bit value or out of range.");
		return false;
	}
	AliHLTUInt32_t specVal = AliHLTUInt32_t(num);
	
	if (fBlockIdList.size() > 0 and not fBlockIdList[fBlockIdList.size()-1].fSpecSet)
	{
		fBlockIdList[fBlockIdList.size()-1].fSpec = specVal;
		fBlockIdList[fBlockIdList.size()-1].fSpecSet = true;
	}
	else
	{
		AliBlockId id;
		id.fType = AliHLTComponentDataTypeInitializer(kAliHLTAnyDataTypeID, kAliHLTDataOriginAny);
		id.fSpec = specVal;
		id.fTypeSet = false;
		id.fOriginSet = false;
		id.fSpecSet = true;
		fBlockIdList.push_back(id);
	}
	
	return true;
}


void AliHLTCorruptorComponent::AddCorruptionCommand(
		AliCorruptionType type,
		AliHLTUInt64_t minrange,
		AliHLTUInt64_t maxrange,
		AliHLTUInt64_t alignment,
		bool minRangeRelative,
		bool maxRangeRelative,
		bool userate,
		double rate,
		AliHLTUInt64_t minerrors,
		AliHLTUInt64_t maxerrors,
		AliHLTUInt64_t burstlength,
		AliHLTUInt64_t burstcount,
		AliHLTUInt64_t minbits,
		AliHLTUInt64_t maxbits,
		size_t firstpattern,
		size_t lastpattern
	)
{
	// Adds a new entry into the corruptions list to apply.
	
	AliCorruptionInfo info;
	memset(&info, 0x0, sizeof(info));
	info.fType = type;
	info.fMinRange = minrange;
	info.fMaxRange = maxrange;
	info.fAlignment = alignment;
	info.fMinRangeRelative = minRangeRelative;
	info.fMaxRangeRelative = maxRangeRelative;
	info.fUseRate = userate;
	if (userate)
	{
		info.fRate = rate;
	}
	else
	{
		info.fMinErrors = minerrors;
		info.fMaxErrors = maxerrors;
	}
	switch (type)
	{
	case kBurstErrors:
		info.fBurstErrorLength = burstlength;
		info.fBurstErrorCount = burstcount;
		break;
	case kReplace:
		info.fFirstPattern = firstpattern;
		info.fLastPattern = lastpattern;
		break;
	case kReplaceRandom:
		info.fMinBits = minbits;
		info.fMaxBits = maxbits;
		break;
	case kInsert:
		info.fFirstPattern = firstpattern;
		info.fLastPattern = lastpattern;
		break;
	case kInsertRandom:
		info.fMinBits = minbits;
		info.fMaxBits = maxbits;
		break;
	case kRemove:
		info.fMinBits = minbits;
		info.fMaxBits = maxbits;
		break;
	default:
		// Nothing to do.
		break;
	}
	fCorruptionInfoList.push_back(info);
}


int AliHLTCorruptorComponent::CheckForCommandWithPatterns(
		const char* name, AliCorruptionType type,
		int& i, int argc, const char** argv,
		AliHLTUInt64_t minrange,
		AliHLTUInt64_t maxrange,
		AliHLTUInt64_t alignment,
		bool minRangeRelative,
		bool maxRangeRelative,
		bool userate,
		double errorrate,
		AliHLTUInt64_t minerrors,
		AliHLTUInt64_t maxerrors
	)
{
	// Check for a particular command line option that has pattern parameters.
	
	AliPattern pattern = {0, 0, false};
	if (strcmp(argv[i], name) == 0)
	{
		// Check for at least 1 pattern.
		if (argc <= i+1)
		{
			HLTError("At least one pattern must be specified for %s.", name);
			return -EINVAL;
		}
		if (not ConvertToPattern(argv[i+1], pattern)) return -EPROTO;
		// Mark the pattern found as the first pattern and add it.
		size_t firstpattern = fPatterns.size();
		fPatterns.push_back(pattern);
		++i;
		// Now check for more patterns. Do not print any error messages,
		// since here we assume that if the string is not a pattern string
		// then it might be another command.
		while (i+1 < argc and ConvertToPattern(argv[i+1], pattern, false))
		{
			fPatterns.push_back(pattern);
			++i;
		}
		// Mark the last pattern and add corruption command.
		size_t lastpattern = fPatterns.size();
		AddCorruptionCommand(
			type, minrange, maxrange, alignment,
			minRangeRelative, maxRangeRelative, userate,
			errorrate, minerrors, maxerrors, 0, 0, 0, 0,
			firstpattern, lastpattern
		);
		return 1;
	}
	return 0;
}


int AliHLTCorruptorComponent::CheckForCommandWithMinMax(
		const char* name, AliCorruptionType type,
		int& i, int argc, const char** argv,
		AliHLTUInt64_t minrange,
		AliHLTUInt64_t maxrange,
		AliHLTUInt64_t alignment,
		bool minRangeRelative,
		bool maxRangeRelative,
		bool userate,
		double errorrate,
		AliHLTUInt64_t minerrors,
		AliHLTUInt64_t maxerrors
	)
{
	// Check for a particular command line option that has 2 (min/max) parameters.
	
	if (strcmp(argv[i], name) == 0)
	{
		if (argc <= i+1)
		{
			HLTError("Minimum number of bits to replace not specified for %s.", name);
			return -EINVAL;
		}
		AliHLTUInt64_t minbits = 0;
		if (strcmp(argv[i+1], "min") == 0)
		{
			minbits = 0;
		}
		else if (strcmp(argv[i+1], "max") == 0)
		{
			minbits = 0xFFFFFFFFFFFFFFFFull;
		}
		else
		{
			if (not ConvertToPositiveInt(argv[i+1], minbits)) return -EPROTO;
		}
		if (argc <= i+2)
		{
			HLTError("Maximum number of bits to replace not specified for %s.", name);
			return -EINVAL;
		}
		AliHLTUInt64_t maxbits = 0;
		if (strcmp(argv[i+2], "min") == 0)
		{
			maxbits = 0;
		}
		else if (strcmp(argv[i+2], "max") == 0)
		{
			maxbits = 0xFFFFFFFFFFFFFFFFull;
		}
		else
		{
			if (not ConvertToPositiveInt(argv[i+2], maxbits)) return -EPROTO;
		}
		if (maxbits < minbits)
		{
			HLTError("Maximum value (%s) is smaller than minimum (%s).", argv[i+2], argv[i+1]);
			return -EPROTO;
		}
		AddCorruptionCommand(
			type, minrange, maxrange, alignment,
			minRangeRelative, maxRangeRelative, userate,
			errorrate, minerrors, maxerrors, 0, 0, minbits, maxbits
		);
		i += 2;
		return 1;
	}
	return 0;
}


Int_t AliHLTCorruptorComponent::DoInit(int argc, const char** argv)
{
	// Initialises the corruptor component from the command line.
	
	int result = 0;
	fBlockIdList.clear();
	fCorruptionInfoList.clear();
	fPatterns.clear();
	AliHLTUInt64_t minpos = 0;
	AliHLTUInt64_t maxpos = 0xFFFFFFFFFFFFFFFFull;
	AliHLTUInt64_t alignment = 1;
	bool minposrel = false;
	bool maxposrel = false;
	bool useErrorRate = true;
	double errorrate = 0.001;
	AliHLTUInt64_t minerrors = 1;
	AliHLTUInt64_t maxerrors = 1;
	bool seedSet = false;
	UInt_t seed = 0;
	
	for (int i = 0; i < argc; ++i)
	{
		if (strcmp(argv[i], "-seed") == 0)
		{
			if (seedSet)
			{
				HLTError("The option \"-seed\" has already been used with"
					" the value %u. Can only use this option once",
					seed
				);
				return -EINVAL;
			}
			if (argc <= i+1)
			{
				HLTError("The random number generator seed was not specified for -seed.");
				return -EINVAL;
			}
			AliHLTUInt64_t tmpseed = 0;
			if (not ConvertToPositiveInt(argv[i+1], tmpseed)) return -EPROTO;
			if (tmpseed > 0xFFFFFFFF)
			{
				HLTError("The random number generator seed value cannot be larger than %u.", 0xFFFFFFFF);
				return -EPROTO;
			}
			seed = UInt_t(tmpseed);
			seedSet = true;
			++i;
			continue;
		}
		
		if (strcmp(argv[i], "-datatype") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The data type identifier was not specified for -datatype.");
				return -EINVAL;
			}
			if (not AddBlockTypeId(argv[i+1])) return -EPROTO;
			if (argc <= i+2)
			{
				HLTError("The origin identifier was not specified for -datatype.");
				return -EINVAL;
			}
			if (not AddBlockOrigin(argv[i+2])) return -EPROTO;
			i += 2;
			continue;
		}
		
		if (strcmp(argv[i], "-origin") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The origin identifier was not specified for -origin.");
				return -EINVAL;
			}
			if (not AddBlockOrigin(argv[i+1])) return -EPROTO;
			++i;
			continue;
		}
		
		if (strcmp(argv[i], "-typeid") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The data type identifier was not specified for -typeid.");
				return -EINVAL;
			}
			if (not AddBlockTypeId(argv[i+1])) return -EPROTO;
			++i;
			continue;
		}
		
		if (strcmp(argv[i], "-dataspec") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The data specification was not specified for -dataspec.");
				return -EINVAL;
			}
			if (not AddBlockSpec(argv[i+1])) return -EPROTO;
			++i;
			continue;
		}
		
		if (strcmp(argv[i], "-range") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("Minimum byte/bit position not given for -range.");
				return -EINVAL;
			}
			if (not ConvertToBitPosition(argv[i+1], minpos, minposrel)) return -EPROTO;
			if (argc <= i+2)
			{
				HLTError("Maximum byte/bit position not given for -range.");
				return -EINVAL;
			}
			if (not ConvertToBitPosition(argv[i+2], maxpos, maxposrel)) return -EPROTO;
			if (minposrel == maxposrel and maxpos < minpos)
			{
				HLTError("Maximum position (%s) is smaller than minimum position (%s).",
					argv[i+2], argv[i+1]
				);
				return -EPROTO;
			}
			i += 2;
			continue;
		}
		
		if (strcmp(argv[i], "-alignment") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The alignment value was not specified for -errorrate.");
				return -EINVAL;
			}
			if (not ConvertToPositiveInt(argv[i+1], alignment)) return -EPROTO;
			if (alignment == 0)
			{
				HLTError("The alignment value cannot be zero.");
				return -EPROTO;
			}
			++i;
			continue;
		}
		
		if (strcmp(argv[i], "-errorrate") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The current error rate to use was not specified for -errorrate.");
				return -EINVAL;
			}
			if (not ConvertToPercentage(argv[i+1], errorrate)) return -EPROTO;
			useErrorRate = true;
			++i;
			continue;
		}
		
		if (strcmp(argv[i], "-errorcount") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("Minimum number of errors not specified for -errorcount.");
				return -EINVAL;
			}
			if (not ConvertToPositiveInt(argv[i+1], minerrors)) return -EPROTO;
			if (argc <= i+2)
			{
				HLTError("Maximum number of errors not specified for -errorcount.");
				return -EINVAL;
			}
			if (not ConvertToPositiveInt(argv[i+2], maxerrors)) return -EPROTO;
			if (maxerrors < minerrors)
			{
				HLTError("Maximum error count (%s) is smaller than minimum (%s).",
					argv[i+2], argv[i+1]
				);
				return -EPROTO;
			}
			useErrorRate = false;
			i += 2;
			continue;
		}
		
		if (strcmp(argv[i], "-singleflips") == 0)
		{
			AddCorruptionCommand(
				kSingleFlips, minpos, maxpos, alignment,
				minposrel, maxposrel, useErrorRate,
				errorrate, minerrors, maxerrors
			);
			continue;
		}
		
		if (strcmp(argv[i], "-bursterrors") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("Maximum burst error length not specified for -bursterrors.");
				return -EINVAL;
			}
			AliHLTUInt64_t burstLength;
			if (not ConvertToPositiveInt(argv[i+1], burstLength)) return -EPROTO;
			if (argc <= i+2)
			{
				HLTError("Number of bit flips in a burst error not specified for -bursterrors.");
				return -EINVAL;
			}
			AliHLTUInt64_t burstCount;
			if (not ConvertToPositiveInt(argv[i+2], burstCount)) return -EPROTO;
			AddCorruptionCommand(
				kBurstErrors, minpos, maxpos, alignment,
				minposrel, maxposrel, useErrorRate, errorrate,
				minerrors, maxerrors, burstLength, burstCount
			);
			i += 2;
			continue;
		}
		
		result = CheckForCommandWithPatterns("-replace", kReplace,
				i, argc, argv, minpos, maxpos, alignment,
				minposrel, maxposrel, useErrorRate,
				errorrate, minerrors, maxerrors
			);
		if (result == 1) continue;
		if (result < 0) return result;
		
		result = CheckForCommandWithPatterns("-insert", kInsert,
				i, argc, argv, minpos, maxpos, alignment,
				minposrel, maxposrel, useErrorRate,
				errorrate, minerrors, maxerrors
			);
		if (result == 1) continue;
		if (result < 0) return result;
		
		result = CheckForCommandWithMinMax("-replace-random", kReplaceRandom,
				i, argc, argv, minpos, maxpos, alignment,
				minposrel, maxposrel, useErrorRate,
				errorrate, minerrors, maxerrors
			);
		if (result == 1) continue;
		if (result < 0) return result;
		
		result = CheckForCommandWithMinMax("-insert-random", kInsertRandom,
				i, argc, argv, minpos, maxpos, alignment,
				minposrel, maxposrel, useErrorRate,
				errorrate, minerrors, maxerrors
			);
		if (result == 1) continue;
		if (result < 0) return result;
		
		result = CheckForCommandWithMinMax("-remove", kRemove,
				i, argc, argv, minpos, maxpos, alignment,
				minposrel, maxposrel, useErrorRate,
				errorrate, minerrors, maxerrors
			);
		if (result == 1) continue;
		if (result < 0) return result;
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	} // for loop
	
	// If no errors were specified then set to 1% single bit corruption.
	if (fCorruptionInfoList.size() == 0)
	{
		minpos = 0;
		maxpos = 0xFFFFFFFFFFFFFFFFull;
		alignment = 1;
		useErrorRate = true;
		errorrate = 0.001;
		AddCorruptionCommand(
			kSingleFlips, minpos, maxpos, alignment, useErrorRate, errorrate
		);
	}
	
	gRandom->SetSeed(0); // Use current time.
	
	HLTInfo("Starting Corruptor Component.");
	return 0;
}


Int_t AliHLTCorruptorComponent::DoDeinit()
{
	// Cleans up the corruptor component.
	fBlockIdList.clear();
	return 0;
}


bool AliHLTCorruptorComponent::ShouldProcess(const AliHLTComponentBlockData* block) const
{
	if (fBlockIdList.size() == 0) return true;
	for (size_t i = 0; i < fBlockIdList.size(); ++i)
	{
		AliBlockId id = fBlockIdList[i];
		AliHLTComponentDataType type = AliHLTComponentDataTypeInitializer(
				id.fTypeSet ? id.fType : kAliHLTAnyDataType,
				id.fOriginSet ? id.fType.fOrigin : kAliHLTDataOriginAny
			);
		if (type != block->fDataType) continue;
		if (id.fSpecSet and id.fSpec != block->fSpecification) continue;
		return true;
	}
	return false;
}


void AliHLTCorruptorComponent::ApplyCorruption(std::vector<bool>& bits) const
{
	// Applies corruption to the bit string.
	
	std::vector<bool> removed;
	std::vector<bool> pattern;
	for (size_t i = 0; i < fCorruptionInfoList.size(); ++i)
	{
		AliCorruptionInfo info = fCorruptionInfoList[i];
#ifdef DEBUG
		const char* cmdtype = NULL;
		switch (info.fType)
		{
		case kSingleFlips:    cmdtype = "kSingleFlips"; break;
		case kBurstErrors:    cmdtype = "kBurstErrors"; break;
		case kReplace:        cmdtype = "kReplace"; break;
		case kReplaceRandom:  cmdtype = "kReplaceRandom"; break;
		case kInsert:         cmdtype = "kInsert"; break;
		case kInsertRandom:   cmdtype = "kInsertRandom"; break;
		case kRemove:         cmdtype = "kRemove"; break;
		default:              cmdtype = "UNKNOWN"; break;
		}
		HLTDebug("Processing corruption command %u of %u, fType = %s,"
			" fMinRange = %llu, fMaxRange = %llu, fAlignment = %llu,"
			" fUseRate = %s, fRate = %.16f, fMinErrors = %llu, fMaxErrors = %llu"
			" fBurstErrorLength = %llu, fBurstErrorCount = %llu,"
			" fMinBits = %llu, fMaxBits = %llu, fFirstPattern = %llu, fLastPattern = %llu.",
			i, fCorruptionInfoList.size(), cmdtype,
			info.fMinRange, info.fMaxRange, info.fAlignment,
			(info.fUseRate ? "true" : "false"), info.fRate,
			info.fMinErrors, info.fMaxErrors,
			info.fBurstErrorLength, info.fBurstErrorCount,
			info.fMinBits, info.fMaxBits,
			AliHLTUInt64_t(info.fFirstPattern), AliHLTUInt64_t(info.fLastPattern)
		);
#endif // DEBUG
		
		// Calculate the range (and adjust for buffer size)
		AliHLTUInt64_t minrange = info.fMinRangeRelative ? bits.size() - info.fMinRange : info.fMinRange;
		AliHLTUInt64_t maxrange = info.fMaxRangeRelative ? bits.size() - info.fMaxRange : info.fMaxRange;
		if (minrange == 0xFFFFFFFFFFFFFFFFull)
		{
			if (info.fType == kInsert or info.fType == kInsertRandom or bits.size() == 0)
			{
				minrange = bits.size();
			}
			else
			{
				minrange = bits.size() - 1;
			}
		}
		if (maxrange == 0xFFFFFFFFFFFFFFFFull)
		{
			if (info.fType == kInsert or info.fType == kInsertRandom or bits.size() == 0)
			{
				maxrange = bits.size();
			}
			else
			{
				maxrange = bits.size() - 1;
			}
		}
		if (maxrange < minrange) maxrange = minrange;
		
		// Select the correct number of errors to generate.
		AliHLTUInt64_t errorcount = 0;
		if (info.fUseRate)
		{
			errorcount = AliHLTUInt64_t(bits.size() * info.fRate);
		}
		else
		{
			AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
						| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
			errorcount = rnum % (info.fMaxErrors - info.fMinErrors + 1) + info.fMinErrors;
		}
		
		// Build the pattern bits for a replace or insert command.
		if (info.fType == kReplace or info.fType == kInsert)
		{
			pattern.clear();
			for (size_t k = info.fFirstPattern; k < info.fLastPattern; ++k)
			{
				const AliPattern& pat = fPatterns[k];
				if (pat.fUseRemoved)
				{
					pattern.insert(pattern.end(), removed.begin(), removed.end());
				}
				else
				{
					for (AliHLTUInt8_t j = 0; j < pat.fWidth; ++j)
					{
						pattern.push_back( ((pat.fPattern >> j) & 0x1) == 0x1 );
					}
				}
			}
		}
		
		// Find a random value for the number of bits to use for the corresponding
		// replace-random, insert-random and remove commands.
		AliHLTUInt64_t bitcount = 0;
		if (info.fType == kReplaceRandom or info.fType == kInsertRandom or info.fType == kRemove)
		{
			AliHLTUInt64_t minbits = info.fMinBits;
			AliHLTUInt64_t maxbits = info.fMaxBits;
			if (minbits == 0xFFFFFFFFFFFFFFFFull) minbits = bits.size();
			if (maxbits == 0xFFFFFFFFFFFFFFFFull) maxbits = bits.size();
			if (maxbits < minbits) maxbits = minbits;
			AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
					      | AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
			bitcount = rnum % (maxbits - minbits + 1) + minbits;
		}
		
		switch (info.fType)
		{
		case kSingleFlips:
			for (AliHLTUInt64_t n = 0; n < errorcount; ++n)
			{
				AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
							| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
				AliHLTUInt64_t pos = rnum % (maxrange - minrange + 1) + minrange;
				pos = ((pos - minrange) / info.fAlignment) * info.fAlignment + minrange;
				if (not (pos < bits.size())) continue;
				bits[pos] = not bits[pos];
			}
			break;
		case kBurstErrors:
			for (AliHLTUInt64_t n = 0; n < errorcount; ++n)
			{
				AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
							| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
				AliHLTUInt64_t firstpos = rnum % (maxrange - minrange + 1) + minrange;
				firstpos = ((firstpos - minrange) / info.fAlignment) * info.fAlignment + minrange;
				AliHLTUInt64_t lastpos = firstpos + info.fBurstErrorLength;
				if (lastpos < firstpos)
				{
					HLTWarning("Burst error length %llu is too large and caused an overflow.",
						info.fBurstErrorLength
					);
					lastpos = firstpos;
				}
				if (lastpos == firstpos) continue;
				for (AliHLTUInt64_t j = 0; j < info.fBurstErrorCount; ++j)
				{
					rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
							| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
					AliHLTUInt64_t pos = rnum % (lastpos - firstpos) + firstpos;
					if (not (pos < bits.size())) continue;
					bits[pos] = not bits[pos];
				}
			}
			break;
		case kReplace:
			for (AliHLTUInt64_t n = 0; n < errorcount; ++n)
			{
				AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
							| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
				AliHLTUInt64_t pos = rnum % (maxrange - minrange + 1) + minrange;
				pos = ((pos - minrange) / info.fAlignment) * info.fAlignment + minrange;
				for (UInt_t j = pos; j < pos + pattern.size() and j < bits.size(); ++j)
				{
					bits[j] = pattern[j-pos];
				}
			}
			break;
		case kReplaceRandom:
			for (AliHLTUInt64_t n = 0; n < errorcount; ++n)
			{
				AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
							| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
				AliHLTUInt64_t pos = rnum % (maxrange - minrange + 1) + minrange;
				pos = ((pos - minrange) / info.fAlignment) * info.fAlignment + minrange;
				for (UInt_t j = pos; j < pos + bitcount and j < bits.size(); ++j)
				{
					bits[j] = gRandom->Integer(2) == 1;
				}
			}
			break;
		case kInsert:
			for (AliHLTUInt64_t n = 0; n < errorcount; ++n)
			{
				AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
							| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
				AliHLTUInt64_t pos = rnum % (maxrange - minrange + 1) + minrange;
				pos = ((pos - minrange) / info.fAlignment) * info.fAlignment + minrange;
				if (not (pos <= bits.size())) continue;
				bits.insert(bits.begin() + pos, pattern.begin(), pattern.end());
			}
			break;
		case kInsertRandom:
			for (AliHLTUInt64_t n = 0; n < errorcount; ++n)
			{
				AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
							| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
				AliHLTUInt64_t pos = rnum % (maxrange - minrange + 1) + minrange;
				pos = ((pos - minrange) / info.fAlignment) * info.fAlignment + minrange;
				if (not (pos <= bits.size())) continue;
				std::vector<bool> newbits;
				for (UInt_t j = 0; j < bitcount; ++j)
				{
					newbits.push_back(gRandom->Integer(2) == 1);
				}
				bits.insert(bits.begin() + pos, newbits.begin(), newbits.end());
			}
			break;
		case kRemove:
			removed.clear();
			for (AliHLTUInt64_t n = 0; n < errorcount; ++n)
			{
				AliHLTUInt64_t rnum = (AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF)) << 32)
							| AliHLTUInt64_t(gRandom->Integer(0xFFFFFFFF));
				AliHLTUInt64_t pos = rnum % (maxrange - minrange + 1) + minrange;
				pos = ((pos - minrange) / info.fAlignment) * info.fAlignment + minrange;
				if (not (pos < bits.size())) continue;
				AliHLTUInt64_t end = pos + bitcount;
				if (end > bits.size()) end = bits.size();
				removed.insert(removed.begin(), bits.begin() + pos, bits.begin() + end);
				bits.erase(bits.begin() + pos, bits.begin() + end);
			}
			break;
		default:
			HLTWarning("Received an unknown corruption type. There must be a program bug!");
			continue;
		}
	}
}


int AliHLTCorruptorComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks, 
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	)
{
	// Processes the input data blocks.
	// The data blocks will be copied to the output and corrupted in some way.
	
	AliHLTUInt8_t* output = outputPtr;
	AliHLTUInt32_t totalSize = 0;
	
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; ++n)
	{
		HLTDebug("Processing block %d, of type '%s' with specification 0x%8.8X.",
			n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
		);
		if (not ShouldProcess(blocks + n))
		{
			Forward();
			continue;
		}
		
		// Convert data block buffer to a bit string.
		std::vector<bool> bits;
		const AliHLTUInt8_t* buffer = reinterpret_cast<const AliHLTUInt8_t*>(blocks[n].fPtr);
		for (AliHLTUInt32_t i = 0; i < blocks[n].fSize; ++i)
		{
			for (AliHLTUInt8_t j = 0; j < 8; ++j)
			{
				bits.push_back( ((buffer[i] >> j) & 0x1) == 0x1 );
			}
		}
		
		ApplyCorruption(bits);
		
		// Convert bit string back to a data buffer.
		if (bits.size() > 0xFFFFFFFFull * 8)
		{
			HLTWarning("The internal bit string representation is too large for the output buffer. Skipping block.");
			continue;
		}
		AliHLTUInt32_t blockSize = bits.size() / 8;
		if (bits.size() % 8 > 0) ++blockSize;
		
		if (totalSize + blockSize > size)
		{
			HLTError("Out of buffer space. Require %d bytes but have %d bytes of buffer space.",
				totalSize + blockSize, size
			);
			fBufferMultiplier *= 8;
			return -ENOSPC;
		}
		
		memset(output, 0x0, blockSize);
		for (size_t i = 0; i < bits.size(); ++i)
		{
			AliHLTUInt32_t byte = i / 8;
			AliHLTUInt32_t bit = i % 8;
			output[byte] = output[byte] | (bits[i] << bit);
		}
		
		AliHLTComponentBlockData bd;
		FillBlockData(bd);
		bd.fPtr = outputPtr;
		bd.fOffset = totalSize;
		bd.fSize = blockSize;
		bd.fDataType = blocks[n].fDataType;
		bd.fSpecification = blocks[n].fSpecification;
		outputBlocks.push_back(bd);
		
		totalSize += blockSize;
		output += blockSize;
	}
	
	size = totalSize; // Set the correct output size.
	return 0;
}
