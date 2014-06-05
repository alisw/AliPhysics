//--------------------------------------------------------------------------
//
// PdfInfo.cc
// Author:  Lynn Garren
//
// Implement operator >> and operator <<
//
// ----------------------------------------------------------------------

#include <iostream>
#include <ostream>
#include <istream>
#include <sstream>

#include "HepMC/PdfInfo.h"
#include "HepMC/StreamHelpers.h"
#include "HepMC/IO_Exception.h"

namespace HepMC {

std::ostream & operator << ( std::ostream & os, PdfInfo const * pdf)
{
    if ( !os ) {
	std::cerr << "operator << for PdfInfo: !os, "
		  << " setting badbit" << std::endl;
	os.clear(std::ios::badbit); 
	return os;
    }
    os << 'F';
    // PdfInfo* is set to 0 by default
    if ( !pdf ) {
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0. );
	detail::output( os, 0. );
	detail::output( os, 0. );
	detail::output( os, 0. );
	detail::output( os, 0. );
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os,'\n');
	return os;
    }
    //
    detail::output( os, pdf->id1() );
    detail::output( os, pdf->id2() );
    detail::output( os, pdf->x1() );
    detail::output( os, pdf->x2() );
    detail::output( os, pdf->scalePDF() );
    detail::output( os, pdf->pdf1() );
    detail::output( os, pdf->pdf2() );
    detail::output( os, pdf->pdf_id1() );
    detail::output( os, pdf->pdf_id2() );
    detail::output( os,'\n');

    return os;
}

std::istream & operator >> (std::istream & is, PdfInfo * pdf)
{
    // make sure the stream is valid
    if ( !is ) {
	std::cerr << "PdfInfo input stream setting badbit." << std::endl;
	is.clear(std::ios::badbit); 
	return is;
    } 
    //
    // get the PdfInfo line
    std::string line;
    std::getline(is,line);
    std::istringstream iline(line);
    std::string firstc;
    iline >> firstc;
    // test to be sure the next entry is of type "F" then ignore it
    if ( firstc != "F" ) {
	std::cerr << "PdfInfo input stream invalid line type: " 
	          << firstc << std::endl;
	// this is non-recoverable, so throw here 
	throw IO_Exception("PdfInfo input stream encounterd invalid data");
    } 
    // read values into temp variables, then create a new PdfInfo object
    int id1 =0, id2 =0, pdf_id1=0, pdf_id2=0;
    double  x1 = 0., x2 = 0., scale = 0., pdf1 = 0., pdf2 = 0.; 
    iline >> id1 ;
    if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
    // check now for empty PdfInfo line
    if( id1 == 0 ) return is;
    // continue reading
    iline >> id2 ;
    if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
    iline >> x1 ;
    if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
    iline >> x2 ;
    if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
    iline >> scale ;
    if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
    iline >> pdf1 ;
    if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
    iline >> pdf2;
    if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
    // check to see if we are at the end of the line
    if( !iline.eof() ) {
        iline >> pdf_id1 ;
        if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
	iline >> pdf_id2;
        if(!iline) throw IO_Exception("PdfInfo input stream encounterd invalid data");
    }
    pdf->set_id1( id1 );
    pdf->set_id2( id2 );
    pdf->set_pdf_id1( pdf_id1 );
    pdf->set_pdf_id2( pdf_id2 );
    pdf->set_x1( x1 );
    pdf->set_x2( x2 );
    pdf->set_scalePDF( scale );
    pdf->set_pdf1( pdf1 );
    pdf->set_pdf2( pdf2 );

    return is;
}

} // HepMC
