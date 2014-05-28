#ifndef HEPMC_GEN_CROSS_SECTION_H
#define HEPMC_GEN_CROSS_SECTION_H

//--------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////
// garren@fnal.gov, May 2009
// 
//////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------

#include <iostream>

namespace HepMC {

//! The GenCrossSection class stores the generated cross section

///
/// \class  GenCrossSection
/// HepMC::GenCrossSection is used to store the generated cross section.
/// This class is meant to be used to pass, on an event by event basis,
/// the current best guess of the total cross section. 
/// It is expected that the final cross section will be stored elsewhere.
///
///   - double cross_section;  		// cross section in pb
///   - double cross_section_error;  	// error associated with this cross section   
/// 
/// The units of cross_section and cross_section_error are expected to be pb.
/// 
/// GenCrossSection information will be written if GenEvent contains a pointer 
/// to a valid GenCrossSection object.
///
class GenCrossSection {

public:
  GenCrossSection()
    : m_cross_section(0),
      m_cross_section_error(0),
      m_is_set(false)
    {}
  ~GenCrossSection() {}

  GenCrossSection( GenCrossSection const & orig ); //!< copy

  void swap( GenCrossSection & other); //!< swap
  GenCrossSection &  operator = ( GenCrossSection const & rhs ); //!< shallow
  /// check for equality
  bool         operator==( const GenCrossSection& ) const;
  /// check for inequality
  bool         operator!=( const GenCrossSection& ) const;


  // ---  accessors:

  /// cross section in pb
  double cross_section()          const { return m_cross_section; }
  /// error associated with this cross section in pb
  double cross_section_error()    const { return m_cross_section_error; }

  /// True if the cross section has been set.  False by default.
  bool   is_set()                 const { return m_is_set; }

  // ---  mutators:
  /// Set cross section and error in pb
  void   set_cross_section( double xs, double xs_err );
  /// set cross section in pb
  void   set_cross_section( double );
  /// set error associated with this cross section in pb
  void   set_cross_section_error( double );
  /// Clear all GenCrossSection info 
  /// (disables output of GenCrossSection until the cross section is set again)
  void   clear();
 
  // ---  I/O:
  /// write to an output stream
  std::ostream &  write( std::ostream & ) const;
  /// read from an input stream
  std::istream &  read( std::istream & );

private: // data members
    double m_cross_section;
    double m_cross_section_error;
    bool   m_is_set;

};

//
// streaming I/O

inline std::ostream & operator << ( std::ostream & os, GenCrossSection & xs )
{ return xs.write(os); }

inline std::istream & operator >> ( std::istream & is, GenCrossSection & xs )
{ return xs.read(is); }

//
// inline methods

inline void GenCrossSection::set_cross_section( double xs, double xserr ) { 
  set_cross_section(xs);
  set_cross_section_error(xserr); 
}

inline void GenCrossSection::set_cross_section( double xs )        
{
  m_cross_section = xs;
  m_is_set = true;
}

inline void GenCrossSection::set_cross_section_error( double xserr )  
{
  m_cross_section_error = xserr;
}

} // HepMC

#endif  // HEPMC_GEN_CROSS_SECTION_H
