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
// $Rev:: 28                          $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-12-10 19:30:01 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef VECTOR3_H
#define VECTOR3_H


#include <iostream>
#include <cmath>


class vector3
{
   public:
      vector3();
      vector3(double *vec);
      vector3(double x, double y, double z);
      virtual ~vector3();
      
      const double* GetVector() const { return _vec; }
      
      void SetVector(double x, double y, double z);
      void SetVector(double *vec);

	    vector3& operator =(const vector3& vec)
	    {
		    if (this != &vec)
			    for (unsigned int i = 0; i < 3; ++i)
				    _vec[i] = vec._vec[i];
		    return *this;
	    }

	    vector3& operator +=(const vector3& vec)
	    {
		    for (unsigned int i = 0; i < 3; ++i)
			    _vec[i] += vec._vec[i];
		    return *this;
	    }
	    vector3& operator -=(const vector3& vec)
	    {
		    for (unsigned int i = 0; i < 3; ++i)
			    _vec[i] -= vec._vec[i];
		    return *this;
	    }

	    double X() const { return _vec[0]; }
	    double Y() const { return _vec[1]; }
	    double Z() const { return _vec[2]; }

	    double Mag2() const { return _vec[0] * _vec[0] + _vec[1] * _vec[1] + _vec[2] * _vec[2]; }
	    double Mag () const { return sqrt(Mag2()); }
      
	    friend std::ostream& operator << (std::ostream&  out,
	                                      const vector3& vec)
	    {
		    out << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << ")";
		    return out;
	    }
	
   private:
      
      double _vec[3];
   
};


#endif  // VECTOR3_H
