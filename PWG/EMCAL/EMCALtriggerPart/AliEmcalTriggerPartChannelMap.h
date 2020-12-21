/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef AliEMCALTRIGGERPARTCHANNELMAP_H
#define AliEMCALTRIGGERPARTCHANNELMAP_H

#include <exception>
#include <sstream>
#include <string>
#include <TObject.h>

namespace PWG {

namespace EMCAL {

namespace TriggerPart {

class AliEmcalTriggerPartChannelMap : public TObject {
public:
	class BoundaryException : public std::exception{
	public:
		BoundaryException():
			exception(),
			fMessage(""),
			fRow(),
			fCol(),
			fNRow(),
			fNCol()
		{
		}
		BoundaryException(int row, int col, int nrow, int ncol):
			exception(),
			fMessage(""),
			fRow(row),
			fCol(col),
			fNRow(nrow),
			fNCol(ncol)
		{
			std::stringstream msgstream;
			msgstream << "Boundary exception: row(" << fRow << ", max " << fNRow -1 << "), col(" << fCol << ", max " << fNCol -1 << ")";
			fMessage = msgstream.str().c_str();
		}
		virtual ~BoundaryException() throw() {}

		int GetRow() const throw() { return fRow; }
		int GetCol() const throw() { return fCol; }
		int GetNRow() const throw() { return fNRow; }
		int GetNCol() const throw() { return fNCol; }

		const char *what() const throw(){
			return fMessage.c_str();
		}
	private:
		std::string 			fMessage;			///< Error message
		int						fRow;				///< Row of the position
		int						fCol;				///< Col of the position
		int						fNRow;				///< Number of rows in the channel map
		int						fNCol;				///< Number of cols in the channel map
	};

	AliEmcalTriggerPartChannelMap(int cols, int rows);
	virtual ~AliEmcalTriggerPartChannelMap();

	void Reset();

	void SetADC(int col, int row, double adc);
	void AddADC(int col, int row, double adc);
	double GetADC(int col, int row) const;
	/**
	 * Get the number of columns in the map
	 * @return The number of colums
	 */
	int GetNumberOfCols() const { return fNADCCols; }
	/**
	 * Get the number of rows in the map
	 * @return The number of cols
	 */
	int GetNumberOfRows() const { return fNADCRows; }

protected:

	inline int GetIndexInArray(int col, int row) const;
	int                     fNADCCols;      ///< Number of columns
	int                     fNADCRows;      ///< Number of rows
	double                  *fADC;          ///< Array of Trigger ADC values

	ClassDef(AliEmcalTriggerPartChannelMap, 1);
};

int AliEmcalTriggerPartChannelMap::GetIndexInArray(int col, int row) const {
	return fNADCCols * row + col;
}

}
}
}
#endif /* AliEMCALTRIGGERPARTCHANNELMAP_H */
