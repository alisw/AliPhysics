/////////////////////////////////////////////////////////////////////////
// CompareGenEvent.cc
//
// garren@fnal.gov, January 2008
// Free functions used to compare two copies of GenEvent
//////////////////////////////////////////////////////////////////////////
//

#include <iostream>

#include "HepMC/CompareGenEvent.h"
#include "HepMC/GenEvent.h"

namespace HepMC {

bool compareGenEvent( GenEvent* e1, GenEvent* e2)
{
   //std::cout << "compareGenEvent: comparing event " << e1->event_number() << " to event " 
   //          << e2->event_number() << std::endl;
   if( e1->event_number() != e2->event_number() ) { 
       std::cerr << "compareGenEvent: event numbers differ " << std::endl;
       return false; 
   }
   if( e1->signal_process_id() != e2->signal_process_id() ) { 
       std::cerr << "compareGenEvent: signal process ids differ " << std::endl;
       return false; 
   }
   if( e1->event_scale() != e2->event_scale() ) { 
       std::cerr << "compareGenEvent: event scales differ " << std::endl;
       return false; 
   }
   if( e1->alphaQCD() != e2->alphaQCD() ) { 
       std::cerr << "compareGenEvent: alphaQCD differs " << std::endl;
       return false; 
   }
   if( e1->alphaQED() != e2->alphaQED() ) { 
       std::cerr << "alphaQED differs " << std::endl;
       return false; 
   }
   if( e1->mpi() != e2->mpi() ) { 
       std::cerr << "compareGenEvent: mpi differs " << std::endl;
       return false; 
   }
   if ( !compareSignalProcessVertex( e1, e2 ) ) { return false; }
   if ( !compareBeamParticles( e1, e2 ) ) { return false; }
   if ( !compareWeights( e1, e2 ) ) { return false; }
   if( e1->random_states() != e2->random_states() ) { 
       std::cerr << "compareGenEvent: random states differ " << std::endl;
       return false; 
   }
   if( e1->heavy_ion() != e2->heavy_ion() ) { 
       std::cerr << "compareGenEvent: heavy ions differ " << std::endl;
       return false; 
   }
   if( e1->pdf_info() != e2->pdf_info() ) { 
       std::cerr << "compareGenEvent: pdf info differs " << std::endl;
       return false; 
   }
   if ( !compareParticles( e1, e2 ) ) { return false; }
   if ( !compareVertices( e1, e2 ) ) { return false; }
   return true;
}

bool compareSignalProcessVertex( GenEvent* e1, GenEvent* e2 ) {
   // compare signal process vertex
   GenVertex* s1 = e1->signal_process_vertex();
   GenVertex* s2 = e2->signal_process_vertex();
   if( s1 && s2 ) {
       if( (*s1) != (*s2) ) { 
	   std::cerr << "compareSignalProcessVertex: signal process vertices differ " << std::endl;
	   return false; 
       }
   }
   return true;
}

bool compareBeamParticles( GenEvent* e1, GenEvent* e2 ) {
   GenParticle* e1b1 = e1->beam_particles().first;
   GenParticle* e1b2 = e1->beam_particles().second;
   GenParticle* e2b1 = e2->beam_particles().first;
   GenParticle* e2b2 = e2->beam_particles().second;
   if( e1b1 && e1b2 && e2b1 && e2b2 ) {
       if( (*e1b1) == (*e2b1)  && (*e1b2) == (*e2b2) ) {
       } else { 
	   std::cerr << "compareBeamParticles: beam particles differ " << std::endl;
	   return false; 
       }
   }
   return true;
}

bool compareWeights( GenEvent* e1, GenEvent* e2 ) {
   if( e1->weights() == e2->weights() ) return true;
   std::cerr << "compareWeights: weight containers differ " << std::endl;
   return false;
}

bool compareParticles( GenEvent* e1, GenEvent* e2 ) {
   if( e1->particles_size() != e2->particles_size() ) { 
       std::cerr << "compareParticles: number of particles differs " << std::endl;
       return false; 
   }
   if( e1->particles_size() == 0 ) { return true; }
   for ( GenEvent::particle_const_iterator p1 =  e1->particles_begin(),
         p2 = e2->particles_begin();
         p1 !=  e1->particles_end(); ++p1, ++p2 ) {
       /* std::cout << "compareParticles: particle " 
		 << (*p1)->barcode() << " " << (*p2)->barcode()
		 << std::endl; */
       if ( **p1 != **p2 ) {
	   std::cerr << "compareParticles: particle " 
	             << (*p1)->barcode() << " differs from " 
		     << (*p2)->barcode() << std::endl;
	   return false; 
       }
   }
   return true;
}

bool compareVertices( GenEvent* e1, GenEvent* e2 ) {
   if( e1->vertices_size() != e2->vertices_size() ) { 
       std::cerr << "compareVertices: number of vertices differs " << std::endl;
       return false; 
   }
   for ( GenEvent::vertex_const_iterator v =  e1->vertices_begin();
         v !=  e1->vertices_end(); ++v ) {
       //std::cout << "compareVertices:  comparing vertex " 
	//         << (*v)->barcode() << std::endl;
       GenVertex* v1 = (*v);
       GenVertex* v2 = e2->barcode_to_vertex((*v)->barcode());
       compareVertex( (*v), e2->barcode_to_vertex((*v)->barcode()));
       if ( (*v1) != (*v2) ) {
	   std::cerr << "compareVertices: vertex " 
	             << (*v)->barcode() << " differs" << std::endl;
	   return false; 
       }
   }
   return true;
}

bool compareVertex( GenVertex* v1, GenVertex* v2 ) {
       if ( v1->position() !=  v2->position() ) {
	  std::cerr << "compareVertex: position " 
	            << v1->barcode() << " differs" << std::endl;
	  return false; 
       }
       // if the size of the inlist differs, return false.
       if ( v1->particles_in_size() !=  v2->particles_in_size() ) {
	  std::cerr << "compareVertex: particles_in_size " 
	            << v1->barcode() << " differs" << std::endl;
	  return false; 
       }
       // loop over the inlist and ensure particles are identical
       if ( v1->particles_in_const_begin() != v1->particles_in_const_end() ) {
	   for ( GenVertex::particles_in_const_iterator 
		    ia = v1->particles_in_const_begin(),
		    ib = v2->particles_in_const_begin();
		ia != v1->particles_in_const_end(); ia++, ib++ ){
	      if ( **ia != **ib ) {
		 std::cerr << "compareVertex: incoming particle " 
	        	   << v1->barcode() << " differs: " 
			   << (*ia)->barcode() << " " << (*ib)->barcode()
			   << std::endl;
		  //return false; 
	      }
	   }
       }
       // if the size of the outlist differs, return false.
       if ( v1->particles_out_size() !=  v2->particles_out_size() ) {
	  std::cerr << "compareVertex: particles_out_size " 
	            << v1->barcode() << " differs" << std::endl;
	  return false; 
       }
       // loop over the outlist and ensure particles are identical
       if ( v1->particles_out_const_begin() != v1->particles_out_const_end() ) {
	   for ( GenVertex::particles_out_const_iterator 
		     ia = v1->particles_out_const_begin(),
		     ib = v2->particles_out_const_begin();
		 ia != v1->particles_out_const_end(); ia++, ib++ ){
	       if ( **ia != **ib ) {
		   std::cerr << "compareVertex: outgoing particle " 
	        	     << v1->barcode() << " differs: " 
			     << (*ia)->barcode() << " " << (*ib)->barcode()
			     << std::endl;
		   //return false; 
	       }
	   }
       }
   return true;
}

} // HepMC
