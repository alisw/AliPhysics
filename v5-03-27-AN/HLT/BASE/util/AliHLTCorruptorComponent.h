//-*- Mode: C++ -*-
// $Id: $
#ifndef ALIHLTCORRUPTORCOMPONENT_H
#define ALIHLTCORRUPTORCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// \file   AliHLTCorruptorComponent.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   5 Aug 2010
/// \brief  Declaration of the AliHLTCorruptorComponent component class.

#include "AliHLTProcessor.h"
#include <vector>

/**
 * \class AliHLTCorruptorComponent
 * The corruptor component is used for testing purposes. It should be used to corrupt
 * input data meant for a normal processing component to check the stable behaviour
 * of such a component to data corruption.
 * The component will copy all input data blocks that it will corrupt to the output
 * and modify the data by corrupting it in some way. This includes removing parts of
 * the data, adding garbage to it or flipping single bits.
 * One can select which data block types to modify; then unmodified data blocks are
 * just forwarded.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b CorruptorComponent <br>
 * Library: \b libAliHLTUtil.so   <br>
 * Input Data Types:  ::kAliHLTAnyDataType <br>
 * Output Data Types: ::kAliHLTAnyDataType <br>
 *
 * <h2>Mandatory arguments:</h2>
 * None.
 *
 * <h2>Optional arguments:</h2>
 * \li -seed <i>number</i> <br>
 *      Sets the random number generator's seed value. The default value of zero
 *      forces the generator to use the current time for the seed.
 * \li -datatype <i>id</i> <i>origin</i> <br>
 *      Data block type ID and origin to filter on. Only blocks with this type are
 *      modified.
 * \li -origin <i>origin</i> <br>
 *      Data block origin to filter on. Only blocks with this origin are modified.
 * \li -typeid <i>id</i> <br>
 *      Data block type ID to filter on. Only blocks with this type ID are modified.
 * \li -dataspec <i>specification</i> <br>
 *      Data block specification to filter on. Only blocks with this specification
 *      are modified.
 * \li -range <i>min</i> <i>max</i> <br>
 *      Sets the range that the current corruption settings will be applied to.
 *      The current corruption settings are selected with "-singleflips",
 *      "-bursterrors", "-replace", "-insert" and "-remove".
 *      The min/max values specify the range in bytes. One can also specify it on
 *      the bit level with the following syntax: "\<byte\>:\<bit\>",
 *      where \<byte\> indicates the byte position counted from zero, \<bit\>
 *      indicates the bit position inside the byte, which can be one of [0..7].
 *      An example of using the "-range" option: "-range 3:2 4:7". This indicates
 *      a range of 14 bits from the 4th byte, bit 3 to the 5th byte up to and
 *      including bit 8.
 *      The special tags "min" and "max" can be used as placeholders for <i>min</i>
 *      and <i>max</i> respectively. These will automatically use the minimum and
 *      maximum size for the data block.
 *      One can also use the syntax: "max-\<byte\>:\<bit\>", which allows to select
 *      positions relative to the end of the block.
 * \li -alignment <i>value</i> <br>
 *      Sets the alignment in bits to use for the applied errors.
 *      For example if one wants to apply errors to locations aligned to 16 bit
 *      words, then <i>value</i> should be 16. With the default value equal to 1.
 *      The alignment is relative to the start of the range given by "-range".
 *      Thus, the location <i>l</i> is calculated with:
 *        <i>l</i> = floor((<i>r</i> - <i>m</i>) / <i>value</i>) * <i>value</i> + <i>m</i>
 *      where <i>r</i> is a random number chosen in the current range, <i>m</i> is
 *      the minimum bit position of the range and <i>value</i> is the alignment value.
 * \li -errorrate <i>rate</i> <br>
 *      Sets the rate of errors to apply to the current range as a percentage or
 *      fraction of the current block size. If no range has been selected then the
 *      rate is for the whole data block by default.
 * \li -errorcount <i>min</i> <i>max</i> <br>
 *      This is an alternative to "-errorrate", where a random number of errors is
 *      generated within the range [<i>min</i>..<i>max</i>]. If <i>min</i> is set
 *      to equal <i>max</i> then a fixed number of errors is always applied to the
 *      current range.
 * \li -singleflips <br>
 *      Generates <i>n</i> single bit errors distributed evenly within the current
 *      data range given by "-range". <i>n</i> is determined either from the rate
 *      given by "-errorrate" or from the range given by "-errorcount".
 * \li -bursterrors <i>length</i> <i>count</i> <br>
 *      Generates burst errors in the current data range. Burst errors have <i>count</i>
 *      number of bit flips within a <i>length</i> number of consecutive bits
 *      closely spaced together. The total number of burst errors generated is
 *      controlled by the "-errorrate" or "-errorcount" options.
 * \li -replace <i>pattern</i> [<i>pattern</i> ...] <br>
 *      Allows one to replace parts of the data with the given pattern or consecutive
 *      patterns. The format of <i>pattern</i> should be "\<value\>/\<width\>",
 *      where \<value\> is a number or hexadecimal value (interpreted in little
 *      endian format) to be used as the replacement. \<width\> is an optional
 *      field to indicate the width of the pattern in bits. The maximum valid value
 *      for \<width\> is 64.
 *      A special symbol "removed" can be used for the pattern, which will cause
 *      the command to use the bits last removed by the "-remove" command. If no
 *      bits were removed or no previous "-remove" command was issued then nothing
 *      is replaced.
 *      Note that the location of the replacement is controlled by the "-range"
 *      option. Normally the position is randomly selected within this range.
 *      For specific placement one can use "-range x x" where x is the exact
 *      position of the replacement. The total number of replacements is
 *      controlled by the "-errorrate" or "-errorcount" options.
 * \li -replace-random <i>min</i> <i>max</i> <br>
 *      Replaces a random number of bits within the range [<i>min</i>..<i>max</i>]
 *      from locations selected within the current range and for a total number of
 *      times controlled by the "-errorrate" or "-errorcount" options.
 *      The placeholder symbols "min" and "max" can be used.
 * \li -insert <i>pattern</i> [<i>pattern</i> ...] <br>
 *      Allows one to insert a pattern or patterns into the data within the current
 *      data range. The format of <i>pattern</i> is the same as for "-replace".
 *      The total number of insertions is controlled by the "-errorrate" or
 *      "-errorcount" options and the locations of the insertions by "-range".
 * \li -insert-random <i>min</i> <i>max</i> <br>
 *      Inserts a random number of bits within the range [<i>min</i>..<i>max</i>]
 *      into locations selected within the current range and for a total number of
 *      times controlled by the "-errorrate" or "-errorcount" options.
 * \li -remove <i>min</i> <i>max</i> <br>
 *      Removes a random number of bits within the range [<i>min</i>..<i>max</i>]
 *      at locations selected within the current range and for a total number of
 *      times controlled by the "-errorrate" or "-errorcount" options.
 *      For removing bytes one needs to specify a multiple of 8. Also, to remove
 *      a fixed number of bits one must set <i>min</i> equal <i>max</i>.
 *
 * \note All corruption operations are applied in the order they as specified on
 *       the command line. So the first argument will be applied first, and the
 *       last one last.
 * \note The default corruption settings is to introduce only 1% single bit flips
 *       if no other settings are set and for all data block types if none specified.
 *
 * <h2>Configuration:</h2>
 * Can only be configured with the command line arguments.
 *
 * <h2>Default CDB entries:</h2>
 * None.
 *
 * <h2>Performance:</h2>
 * Depends on complexity of corruption, but typically should easily run with a
 * 3kHz event rate.
 *
 * <h2>Memory consumption:</h2>
 * The amount of memory is used as up to 8 times the size of the largest input
 * data block.
 *
 * <h2>Output size:</h2>
 * The same as the input data size, unless garbage insertion options are used.
 * In that case the output size depends on the exact settings.
 *
 * \ingroup alihlt_util_components
 */
class AliHLTCorruptorComponent : public AliHLTProcessor
{
public:
	
	AliHLTCorruptorComponent();
	virtual ~AliHLTCorruptorComponent();
	
	// Methods inherited from AliHLTComponent:
	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();
	virtual Int_t DoInit(int argc, const char** argv);
	virtual Int_t DoDeinit();
	
protected:
	
	// Method inherited from AliHLTProcessor:
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks, 
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr, 
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		);
	
	using AliHLTProcessor::DoEvent;
	
private:
	
	enum AliCorruptionType
	{
		kUnknownErrorType = 0,  // Unknown corruption type.
		kSingleFlips,           // Single bit flip errors.
		kBurstErrors,           // Burst errors.
		kReplace,               // Replacement of data.
		kReplaceRandom,         // Replacement of data with random bits.
		kInsert,                // Insertion of data.
		kInsertRandom,          // Insertion of random data.
		kRemove,                // Removal of data.
	};
	
	struct AliBlockId
	{
		AliHLTComponentDataType fType;  // Data type and origin.
		AliHLTUInt32_t fSpec;           // Data block specification.
		bool fTypeSet;    // Flag indicating if the fType.fID field was set.
		bool fOriginSet;  // Flag indicating if the fType.fOrigin field was set.
		bool fSpecSet;    // Flag indicating if the fSpec field was set.
	};
	
	typedef std::vector<AliBlockId> AliBlockIdList;
	
	struct AliPattern
	{
		AliHLTUInt64_t fPattern;  // The bit pattern.
		AliHLTUInt8_t fWidth;     // The width of the bit pattern in bits.
		bool fUseRemoved;         // Indicates if the removed bits should be used.
	};
	
	typedef std::vector<AliPattern> AliPatternList;
	
	struct AliCorruptionInfo
	{
		AliCorruptionType fType;  // The type of error/modification being applied.
		AliHLTUInt64_t fMinRange;  // The minimum range to apply the errors in bits.
		AliHLTUInt64_t fMaxRange;  // The maximum range to apply the errors in bits.
		AliHLTUInt64_t fAlignment;  // The alignment to use for error location selection in bits.
		bool fMinRangeRelative; // Indicates that the fMinRange value is relative to the end of the buffer.
		bool fMaxRangeRelative; // Indicates that the fMaxRange value is relative to the end of the buffer.
		bool fUseRate;  // Indicates if the rate or the range should be used.
		union
		{
			double fRate;  // The rate that should be used if fUseRate == true.
			struct
			{
				AliHLTUInt64_t fMinErrors;  // The minimum error count that should be used if fUseRate == false.
				AliHLTUInt64_t fMaxErrors;  // The maximum error count that should be used if fUseRate == false.
			};
		};
		union
		{
			struct
			{
				AliHLTUInt64_t fBurstErrorLength;  // Maximum burst error length in bits.
				AliHLTUInt64_t fBurstErrorCount;   // Number of bit flips in a burst error.
			};
			struct
			{
				AliHLTUInt64_t fMinBits;  // Minimum number of bits to replace/insert/remove.
				AliHLTUInt64_t fMaxBits;  // Maximum number of bits to replace/insert/remove.
			};
			struct
			{
				size_t fFirstPattern;  // First pattern to use for replacement or insertion.
				size_t fLastPattern;   // Last pattern to use for replacement or insertion.
			};
		};
	};
	
	typedef std::vector<AliCorruptionInfo> AliCorruptionInfoList;
	
	/**
	 * Converts the value to a positive 64 bit integer.
	 * \param [in]  value  The command line parameter to convert.
	 * \param [out] num  The result is stored in this variable.
	 * \param [in]  printErrors Indicates if error messages should be generated.
	 * \returns true if the convertion was successful.
	 */
	bool ConvertToPositiveInt(const char* value, AliHLTUInt64_t& num, bool printErrors = true) const;
	
	/**
	 * Converts the value to a bit position.
	 * \param [in]  value  The command line parameter to convert.
	 * \param [out] pos  The result is stored in this variable.
	 * \param [out] relative  Flag filled to indicate if the resulting position
	 *                        is relative to the end of the data block buffer.
	 * \param [in]  printErrors Indicates if error messages should be generated.
	 * \returns true if the convertion was successful.
	 */
	bool ConvertToBitPosition(const char* value, AliHLTUInt64_t& pos, bool& relative, bool printErrors = true) const;
	
	/**
	 * Converts the value to a percentage fraction in the range [0..1].
	 * \param [in]  value  The command line parameter to convert.
	 * \param [out] num  The result is stored in this variable.
	 * \param [in]  printErrors Indicates if error messages should be generated.
	 * \returns true if the convertion was successful.
	 */
	bool ConvertToPercentage(const char* value, double& num, bool printErrors = true) const;
	
	/**
	 * Converts a string to a pattern for the -replace and -insert command line options.
	 * \param [in]  value  The command line parameter to convert.
	 * \param [out] pattern  The resulting pattern structure.
	 * \param [in]  printErrors Indicates if error messages should be generated.
	 * \returns true if the convertion was successful.
	 */
	bool ConvertToPattern(const char* value, AliPattern& pattern, bool printErrors = true) const;
	
	/**
	 * Adds a data block type to the current data block filter rule.
	 * A new filter rule is created if the type is already set.
	 * \param type  The type of the data block.
	 * \returns true if the type string was valid.
	 */
	bool AddBlockTypeId(const char* type);
	
	/**
	 * Adds a data block origin to the current data block filter rule.
	 * A new filter rule is created if the origin is already set.
	 * \param origin  The origin of the data block.
	 * \returns true if the origin string was valid.
	 */
	bool AddBlockOrigin(const char* origin);
	
	/**
	 * Adds a data block specification to the current data block filter rule.
	 * A new filter rule is created if the specification is already set.
	 * \param spec  The specification of the data block.
	 * \returns true if the specification string was a valid string.
	 */
	bool AddBlockSpec(const char* spec);
	
	/**
	 * Adds a corruption command as a AliCorruptionInfo structure to the list
	 * of commands to apply to the input data blocks.
	 */
	void AddCorruptionCommand(
			AliCorruptionType type,
			AliHLTUInt64_t minrange,
			AliHLTUInt64_t maxrange,
			AliHLTUInt64_t alignment,
			bool minRangeRelative,
			bool maxRangeRelative,
			bool userate = true,
			double rate = 0,
			AliHLTUInt64_t minerrors = 0,
			AliHLTUInt64_t maxerrors = 0,
			AliHLTUInt64_t burstlength = 0,
			AliHLTUInt64_t burstcount = 0,
			AliHLTUInt64_t minbits = 0,
			AliHLTUInt64_t maxbits = 0,
			size_t firstpattern = 0,
			size_t lastpattern = 0
		);
	
	/// Parses command line parameters for a command with pattern fields.
	int CheckForCommandWithPatterns(
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
		);
	
	/// Parses command line parameters for a command with min and max fields.
	int CheckForCommandWithMinMax(
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
		);
	
	/// Checks to see if a data block should be processed or not.
	bool ShouldProcess(const AliHLTComponentBlockData* block) const;
	
	/// Applies the corruption commands to a bit string.
	void ApplyCorruption(std::vector<bool>& bits) const;
	
	double fBufferMultiplier;  // The multiplier used by GetOutputDataSize.
	AliBlockIdList fBlockIdList;  // The list of data block filters.
	AliCorruptionInfoList fCorruptionInfoList;  // The list of corruptions to be applied.
	AliPatternList fPatterns;  // Patterns to use for replacements or insertions.
	
	ClassDef(AliHLTCorruptorComponent, 0)  // Data corruption component for testing.
};

#endif // ALIHLTCORRUPTORCOMPONENT_H
