///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 213                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2015-08-15 23:08:02 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H


class randomGenerator
{
	public:
	void SetSeed(unsigned int seed);
	double Rndom(int i=0);
	
	private:
	unsigned int _Mt[624];
	int _count624;
	
};


#endif  // RANDOMGENERATOR_H
	
