 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// --- Standard library ---
#include <Riostream.h>

// --- AliRoot header files ---
#include "AliEMCALHit.h"

using std::endl;

/// \cond CLASSIMP
ClassImp(AliEMCALHit) ;
/// \endcond

///
/// Default constructor
//______________________________________________________________________
AliEMCALHit::AliEMCALHit()
  : fId(0),
    fELOS(0.),
    fPrimary(0),
    fPx(0.),
    fPy(0.),
    fPz(0.),
    fPe(0.),
    fIparent(0),
    fIenergy(0.),
    fTime(0.)
{ }

///
/// Copy constructor
//______________________________________________________________________
AliEMCALHit::AliEMCALHit(const AliEMCALHit & hit) 
  : AliHit(hit),
    fId(hit.fId),
    fELOS(hit.fELOS),
    fPrimary(hit.fPrimary),
    fPx(hit.fPx),
    fPy(hit.fPy),
    fPz(hit.fPz),
    fPe(hit.fPe),
    fIparent(hit.fIparent),
    fIenergy(hit.fIenergy),
    fTime(hit.fTime)
{ }

///
/// Assignment operator; use copy constructor
//_____________________________________________________________________
AliEMCALHit& AliEMCALHit::operator = (const AliEMCALHit &source)
{ 
  if (&source == this) return *this;
  
  new (this) AliEMCALHit(source);
  return *this;
}

///
/// Create an EMCal hit object
///
/// \param shunt: level of primary selection to store
/// \param primary: index label
/// \param track: index
/// \param iparent: index
/// \param ienergy: deposited energy?
/// \param id: cell id?
/// \param hits: hit position time energy
/// \param p: hit momentum
//______________________________________________________________________
AliEMCALHit::AliEMCALHit(Int_t shunt, Int_t primary, Int_t track,
                         Int_t iparent, Float_t ienergy, Int_t id,
                         Float_t *hits,Float_t *p)
  : AliHit(shunt, track),
    fId(id),
    fELOS(0.),
    fPrimary(primary),
    fPx(0.),
    fPy(0.),
    fPz(0.),
    fPe(0.),
    fIparent(iparent),
    fIenergy(ienergy),
    fTime(0.)
{
  fX          = hits[0];
  fY          = hits[1];
  fZ          = hits[2];
  fTime       = hits[3] ;
  fELOS       = hits[4];
  fPx         = p[0];
  fPy         = p[1];
  fPz         = p[2];
  fPe         = p[3];
}

///
/// Two hits are identical if they have the same Id and originat
/// from the same enterring Particle 
//______________________________________________________________________
Bool_t AliEMCALHit::operator==(AliEMCALHit const &rValue) const
{ 
  Bool_t rv = kFALSE;
  
  if ( (fId == rValue.GetId()) && ( fIparent == rValue.GetIparent()) )
    rv = kTRUE;
  
  return rv;
}

///
/// Add the energy of the hit
//______________________________________________________________________
AliEMCALHit AliEMCALHit::operator+(const AliEMCALHit &rValue)
{
  fELOS += rValue.GetEnergy() ;
  
  if(rValue.GetTime() < fTime)
    fTime = rValue.GetTime() ;
  
  return *this;
}

///
/// Dump hit info
//______________________________________________________________________
ostream& operator << (ostream& out,AliEMCALHit& hit)
{
  out << "AliEMCALHit:";
  out << "id=" <<  hit.GetId();
  out << ", Eloss=" <<  hit.GetEnergy();
  out << ", Time=" << hit.GetTime();
  out << "GeV , Track no.=" << hit.GetPrimary();
  out << ", (xyz)=(" << hit.X()<< ","<< hit.Y()<< ","<<hit.Z()<<") cm";
  out << ", fTrack=" << hit.GetTrack();
  out << ", P=(" << hit.GetPx() << "," << hit.GetPy() << "," << hit.GetPz()
  << "," <<hit.GetPe() << ") GeV"  ;
  out << ", Enterring particle ID" << hit.GetIparent();
  out << ", Enterring particle initial energy = " << hit.GetIenergy() << " GeV" ;
  out << endl;
  
  return out;
}
