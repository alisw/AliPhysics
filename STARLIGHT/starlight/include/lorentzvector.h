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


#ifndef LORENTZVECTOR_H
#define LORENTZVECTOR_H


#include "vector3.h"
#include <vector>


class lorentzVector
{
   public:
      
      lorentzVector();
      virtual ~lorentzVector();
      
      lorentzVector(double x, double y, double z, double t);
      
      void SetXYZT(double x, double y, double z, double t);
      void SetPxPyPzE(double px, double py, double pz, double e) { SetXYZT(px, py, pz, e); };
      
      double GetPx() const { return fSpaceVec.GetVector()[0]; }
      double GetPy() const { return fSpaceVec.GetVector()[1]; }
      double GetPz() const { return fSpaceVec.GetVector()[2]; }
      double GetE() const { return fTime; }

	    lorentzVector& operator +=(const lorentzVector& vec)
	    {
		    fSpaceVec += vec.fSpaceVec;
		    fTime     += vec.fTime;
		    return *this;
	    }
	    lorentzVector& operator -=(const lorentzVector& vec)
	    {
		    fSpaceVec -= vec.fSpaceVec;
		    fTime     -= vec.fTime;
		    return *this;
	    }

	    double M2() const { return fTime * fTime - fSpaceVec.Mag2(); }
      double M () const
	    {
	      const double mag2 = M2();
	      return (mag2 < 0) ? -sqrt(-mag2) : sqrt(mag2);
      }

	    vector3 BoostVector() const
	    { return vector3(fSpaceVec.X() / fTime, fSpaceVec.Y() / fTime, fSpaceVec.Z() / fTime); }
	    void Boost(const vector3& beta)
	    {
		    const double beta2        = beta.Mag2();
		    const double gamma        = 1 / sqrt(1 - beta2);
		    const double betaTimesMom = beta.X() * fSpaceVec.X() + beta.Y() * fSpaceVec.Y() + beta.Z() * fSpaceVec.Z();
		    const double gamma2       = (beta2 > 0) ? (gamma - 1) / beta2 : 0;
		    SetXYZT(fSpaceVec.X() + gamma2 * betaTimesMom * beta.X() + gamma * beta.X() * fTime,
		            fSpaceVec.Y() + gamma2 * betaTimesMom * beta.Y() + gamma * beta.Y() * fTime,
		            fSpaceVec.Z() + gamma2 * betaTimesMom * beta.Z() + gamma * beta.Z() * fTime,
		            gamma * (fTime + betaTimesMom));
	    }
      
	    friend std::ostream& operator << (std::ostream&        out,
	                                      const lorentzVector& vec)
	    {
		    out << "(" << vec.GetPx() << ", " << vec.GetPy() << ", " << vec.GetPz()
		        << "; " << vec.GetE() << ")";
		    return out;
	    }

   private:
      
      vector3 fSpaceVec;
      double fTime;
      
};


#endif  // LORENTZVECTOR_H
