// from Andy Buckley

#include "HepMC/GenEvent.h"

  void filterEvent(HepMC::GenEvent* ge) {
    // We treat beam particles a bit specially
    const std::pair<HepMC::GenParticle*, HepMC::GenParticle*> beams = ge->beam_particles();

    // Make list of non-physical particle entries
    std::vector<HepMC::GenParticle*> unphys_particles;
    for (HepMC::GenEvent::particle_const_iterator pi = ge->particles_begin(); 
         pi != ge->particles_end(); ++pi) {
      // Beam particles might not have status = 4, but we want them anyway
      if (beams.first == *pi || beams.second == *pi) continue;
      // Filter by status
      const int status = (*pi)->status();
      if (status != 1 && status != 2 && status != 4) {
        unphys_particles.push_back(*pi);
      }      
    }

    // Remove each unphysical particle from the list
    while (unphys_particles.size()) {
      HepMC::GenParticle* gp = unphys_particles.back();

      // Get start and end vertices
      HepMC::GenVertex* vstart = gp->production_vertex();
      HepMC::GenVertex* vend = gp->end_vertex();

      if (vend == vstart) {
        // Deal with loops
        vstart->remove_particle(gp);
      } else {

        // Connect decay particles from end vertex to start vertex
        /// @todo Have to build a list, since the GV::add_particle_out method modifies the end vertex!
        if (vend && vend->particles_out_size()) {
          std::vector<HepMC::GenParticle*> end_particles;
          for (HepMC::GenVertex::particles_out_const_iterator gpe = vend->particles_out_const_begin(); 
               gpe != vend->particles_out_const_end(); ++gpe) {
            end_particles.push_back(*gpe);
          }
          // Reset production vertices of child particles to bypass unphysical particle
          for (std::vector<HepMC::GenParticle*>::const_iterator gpe = end_particles.begin();
               gpe != end_particles.end(); ++gpe) {
            //std::cout << vstart << ", " << vend << std::endl;
            if (vstart) vstart->add_particle_out(*gpe);
          }
        } else {
          // If null end_vertex... stable unphysical particle?
        }

        // Delete unphysical particle and its orphaned end vertex
        delete vend;
        if (vstart) {
          delete vstart->remove_particle(gp);
        }// else {
        /// @todo Why does this cause an error?
        //  delete gp;
        //}
      }

      // Remove deleted particle from list
      unphys_particles.pop_back();
      //std::cout << unphys_particles.size() << std::endl;
    }

    // Delete any orphaned vertices
    std::vector<HepMC::GenVertex*> orphaned_vtxs;
    for (HepMC::GenEvent::vertex_const_iterator vi = ge->vertices_begin(); 
         vi != ge->vertices_end(); ++vi) {
      if ((*vi)->particles_in_size() == 0 && (*vi)->particles_out_size() == 0) {
        orphaned_vtxs.push_back(*vi);
      }
    }
    while (orphaned_vtxs.size()) {
      delete orphaned_vtxs.back();
      orphaned_vtxs.pop_back();
    }
  }


