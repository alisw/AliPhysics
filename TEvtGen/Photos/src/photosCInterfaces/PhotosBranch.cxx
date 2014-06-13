#include <vector>
#include <list>
#include "PH_HEPEVT_Interface.h"
#include "PhotosParticle.h"
#include "PhotosBranch.h"
#include "Photos.h"
#include "Log.h"
using std::vector;
using std::list;
using std::endl;

namespace Photospp
{

PhotosBranch::PhotosBranch(PhotosParticle* p)
{
	daughters = p->getDaughters();

	//Suppress if somehow got stable particle
	if(daughters.size()==0)
	{
		Log::Debug(1)<<"Stable particle."<<endl;
		suppression = 1;
		forcing     = 0;
		particle    = NULL;
		return;
	}
	else if(daughters.at(0)->getMothers().size()==1)
	{
		// Regular case - one mother
		Log::Debug(1)<<"Regular case."<<endl;
		particle  = p;
		mothers   = p->findProductionMothers();
	}
	else
	{
		// Advanced case - branch with multiple mothers - no mid-particle
		Log::Debug(1)<<"Advanced case."<<endl;
		particle  = NULL;
		mothers   = daughters.at(0)->getMothers();
	}

  //--------------------------------------------------
  // Finalize suppression/forcing checks
  // NOTE: if user forces decay of specific particle,
  //       this overrides any suppresion
  //--------------------------------------------------

	forcing = checkForcingLevel();
	if(!forcing) suppression = checkSuppressionLevel();
	else         suppression = 0;

	// Even if forced or passed suppression check, we still have to check few things
	if(!suppression)
	{
		// Check momentum conservation
		suppression=!checkMomentumConservation();
		if(suppression) Log::Warning()<<"Branching ignored due to 4-momentum non conservation"<<endl;

		// Check if advanced case has only one daughter
		if(!particle && daughters.size()==1) suppression=-1;

		// If any of special cases is true, we're not forcing this branch
		if(suppression) forcing=0;
	}
}

void PhotosBranch::process()
{
	Log::Debug(703)<<"   Processing barcode: "<<( (particle) ? particle->getBarcode() : ( (mothers.size()) ? mothers.at(0)->getBarcode() : -1) )<<endl;
	/*
	cout<<"Particles send to photos (with barcodes in brackets):"<<endl;
	vector<PhotosParticle *> get = getParticles();
	for(int i=0;i<(int)get.size();i++) cout<<"ID: "<<get.at(i)->getPdgID()<<" ("<<get.at(i)->getBarcode()<<"), "; cout<<endl;
	*/
	int index = PH_HEPEVT_Interface::set(this);
	PH_HEPEVT_Interface::prepare();
	PHOTOS_MAKE_C(index);
	PH_HEPEVT_Interface::complete();
	PH_HEPEVT_Interface::get();
	checkMomentumConservation();
}

vector<PhotosParticle *> PhotosBranch::getParticles()
{
	vector<PhotosParticle *> ret = mothers;
	if(particle) ret.push_back(particle);
	ret.insert(ret.end(),daughters.begin(),daughters.end());
	return ret;
}

bool PhotosBranch::checkMomentumConservation()
{
	if(particle)           return particle->checkMomentumConservation();
	if(mothers.size()>0)   return mothers.at(0)->checkMomentumConservation();
	return true;
}

vector<PhotosBranch *> PhotosBranch::createBranches(vector<PhotosParticle *> particles)
{
	Log::Debug(700)<<"PhotosBranch::createBranches - filtering started"<<endl;
	list<PhotosParticle *> list(particles.begin(),particles.end());
	vector<PhotosBranch *> branches;

	// First - add all forced decays
	if(Photos::forceBremList)
	{
		std::list<PhotosParticle *>::iterator it;
		for(it=list.begin();it!=list.end();it++)
		{
			PhotosBranch *branch = new PhotosBranch(*it);
			int forcing = branch->getForcingStatus();
			if(forcing)
			{
				Log::Debug(701)<<" Forced: "<<(*it)->getPdgID()<<" (barcode: "<<(*it)->getBarcode()<<") with forcing status= "<<forcing<<endl;
				branches.push_back(branch);
				it = list.erase(it);
				--it;
				// If forcing consecutive decays
				if(forcing==2)
				{
					PhotosParticle *p = branch->getDecayingParticle();
					if(!p)
					{
						if(branch->getMothers().size()>0) p = branch->getMothers().at(0);
						else continue;
					}
					vector<PhotosParticle *> tree = p->getDecayTree();
					//Add branches for all particles from the list - max O(n*m)
					std::list<PhotosParticle *>::iterator it2;
					for(it2=list.begin();it2!=list.end();it2++)
					{
						for(int i=0;i<(int)tree.size();i++)
						{
							if(tree.at(i)->getBarcode()==(*it2)->getBarcode())
							{
								PhotosBranch *b = new PhotosBranch(*it2);
								branches.push_back(b);
								// If we were to delete our next particle in line
								if(it==it2) --it;
								it2 = list.erase(it2);
								--it2;
								break;
							}
						}
					}
				}
			}
			else delete branch;
		}
	}
	// Quit if we're suppressing everything
	if(Photos::isSuppressed) return branches;
	// Now - check if remaining decays are suppressed
	while(!list.empty())
	{
		PhotosParticle *particle = list.front();
		list.pop_front();
		if(!particle) continue;

		PhotosBranch *branch = new PhotosBranch(particle);
		int suppression = branch->getSuppressionStatus();
		if(!suppression) branches.push_back(branch);
		else
		{
			Log::Debug(702)<<"  Suppressed: "<<particle->getPdgID()<<" (barcode: "<<particle->getBarcode()<<") with suppression status= "<<suppression<<endl;
			//If suppressing consecutive decays
			if(suppression==2)
			{
				PhotosParticle *p = branch->getDecayingParticle();
				if(!p)
				{
					if(branch->getMothers().size()>0) p = branch->getMothers().at(0);
					else continue;
				}
				vector<PhotosParticle *> tree = p->getDecayTree();
				//Remove all particles from the list - max O(n*m)
				std::list<PhotosParticle *>::iterator it;
				for(it=list.begin();it!=list.end();it++)
				{
					for(int i=0;i<(int)tree.size();i++)
					{
						if(tree.at(i)->getBarcode()==(*it)->getBarcode())
						{
							it = list.erase(it);
							--it;
							break;
						}
					}
				}
			}
			delete branch;
			continue;
		}

		//In case we don't have mid-particle erase rest of the mothers from list
		if(!branch->getDecayingParticle())
		{
			vector<PhotosParticle *> mothers = branch->getMothers();
			for(int i=0;i<(int)mothers.size();i++)
			{
				PhotosParticle *m = mothers.at(i);
				if(m->getBarcode()==particle->getBarcode()) continue;
				std::list<PhotosParticle *>::iterator it;
				for(it=list.begin();it!=list.end();it++)
					if(m->getBarcode()==(*it)->getBarcode())
					{
						it = list.erase(it);
						break;
					}
			}
		}
	}
	return branches;
}

int PhotosBranch::checkList(bool forceOrSuppress)
{
	vector< vector<int>* > *list = (forceOrSuppress) ? Photos::forceBremList : Photos::supBremList;
	if(!list) return 0;

	// Can't check without pdgid
	int motherID;
	if(particle) motherID = particle->getPdgID();
	else
	{
		if(mothers.size()==0) return 0;
		motherID = mothers.at(0)->getPdgID();
	}

	// Create list of daughters
	vector<int> dID;
	for(int j=0;j<(int)daughters.size();j++) dID.push_back(daughters[j]->getPdgID());

	vector< vector<int> *> &patternList = *list;

	// Check if the mother and list of daughters matches any of the declared patterns
	for(int j=0; j<(int)patternList.size();j++)
	{
		// Skip patterns that don't have our mother
		if(motherID!=(*patternList[j])[0]) continue;

		// Compare decay daughters with pattern - max O(n*m)
		vector<int> &pattern = *patternList[j];
		bool fullMatch=true;
		for(int k = 1; k<(int)pattern.size()-1; k++)
		{
			bool oneMatch=false;
			for(int l=0;l<(int)dID.size(); l++)
				if(pattern[k]==dID[l]) { oneMatch=true; break; }
			if(!oneMatch) { fullMatch=false; break; }
		}
		// Check if the matching pattern is set for consecutive suppression
		/*
		  Currently minimal matching is used.
		  If e.g. 25 -> -15 is suppressed, then 25 -> 15,-15 and 25 -> 15,-15,22 etc. is suppressed too
		  For exact matching (then suppress 25 -> 15,-15 ; 25 -> 15,-15,22 etc. must be done separately) uncoment line ...:
		*/
		// if(pattern.size()<=2 || (fullMatch && dID.size()==pattern.size()-2) )
		// ...and comment out line:
		if(pattern.size()<=2 || fullMatch)
			return (pattern.back()==1) ? 2 : 1;
	}
	return 0;
}

} // namespace Photospp
