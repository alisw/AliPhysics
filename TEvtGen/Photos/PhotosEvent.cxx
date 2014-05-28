#include <vector>
#include "PhotosParticle.h"
#include "PhotosBranch.h"
#include "PhotosEvent.h"
#include "Log.h"
using std::vector;

namespace Photospp
{

PhotosEvent::~PhotosEvent()
{
	while(m_branch_points.size()!=0)
	{
		PhotosBranch *temp = m_branch_points.back();
		m_branch_points.pop_back();
		delete temp;
	}
}

void PhotosEvent::process()
{
	//print();
	vector<PhotosParticle*> particles = filterParticles( getParticleList() );
	m_branch_points = PhotosBranch::createBranches(particles);

	for(int i=0;i<(int)m_branch_points.size();i++)
		m_branch_points.at(i)->process();
	//print();
}

vector<PhotosParticle *> PhotosEvent::filterParticles(vector<PhotosParticle *> particles)
{
	vector<PhotosParticle *> filtered;
	for(int i=0;i<(int)particles.size();i++)
	{
		PhotosParticle *p = particles.at(i);
		if(!p) continue;

		//check that the particle decays
		if(p->getStatus()==PhotosParticle::STABLE) continue;

		//check for self decays
		vector<PhotosParticle *> daughters = p->getDaughters();
		int j=0;
		for(j=0;j<(int)daughters.size();j++)
			if(daughters.at(j)->getPdgID()==p->getPdgID()) break;
		if(j!=(int)daughters.size()) continue;

		Log::Debug(2)<<"Passed particle filter"<<endl;
		filtered.push_back(p);
	}
	return filtered;
}

} // namespace Photospp
