//-*- Mode: C++ -*-
// $Id$

//
// C++ Interface: AliHLTPHOSTRUMapper
//
// Description: 
//
//
// Author: Øystein Djuvsland <oystein.djuvsland@gmail.com>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
 /**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 

#ifndef ALIHLTPHOSTRUMAPPER_H
#define ALIHLTPHOSTRUMAPPER_H

/**
	@author Øystein Djuvsland <oystein.djuvsland@gmail.com>
*/
class AliHLTPHOSTRUMapper : public AliHLTPHOSBase
{
public:
    AliHLTPHOSTRUMapper();

    ~AliHLTPHOSTRUMapper();
    
    Int_t GetChannelMask(Int_t x, Int_t z, Short_t *maskArray);
    Int_t GetTRUMask(Int_t x, Int_t z, Short_t *maskArray);

};

#endif
