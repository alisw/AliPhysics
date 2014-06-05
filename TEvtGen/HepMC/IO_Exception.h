#ifndef HEPMC_IO_EXCEPTION_H
#define HEPMC_IO_EXCEPTION_H
// ----------------------------------------------------------------------
//
// IO_Exception.h
// Author:  Lynn Garren
//
// IO exception handling
// IO_GenEvent, etc. catch the throw and set data members with the error type and message 
// Some of the messages are constructed with transient information 
//      (e.g., contents of a bad GenParticle)
//
// ----------------------------------------------------------------------



#include <stdexcept>
 
namespace HepMC {

//! IO exception handling

///
/// \class  IO_Exception
/// IO_GenEvent, etc. catch the throw and set data members with the error type and message 
/// Some of the messages are constructed with transient information 
///      (e.g., contents of a bad GenParticle)
class IO_Exception : public std::runtime_error {
public:
  IO_Exception(const std::string & msg) 
  : std::runtime_error(msg) { }

  /// IO error types
  enum ErrorType{ OK,
                  NullEvent, 
                  WrongFileType, 
                  MissingStartKey, 
		  EndOfStream, 
		  EndKeyMismatch, 
		  MissingEndKey, 
		  InvalidData,
                  InputAndOutput,
		  BadOutputStream,
		  BadInputStream };

};

}	// namespace HepMC

#endif // HEPMC_IO_EXCEPTION_H
