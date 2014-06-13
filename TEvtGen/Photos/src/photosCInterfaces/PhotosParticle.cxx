#include <vector>
#include <math.h>
#include "PhotosParticle.h"
#include "Log.h"
using std::vector;

namespace Photospp
{

bool PhotosParticle::hasDaughters()
{
	if(getDaughters().size()==0) return false;
	else                         return true;
}

PhotosParticle * PhotosParticle::findLastSelf()
{
	vector<PhotosParticle*> daughters = getDaughters();
	vector<PhotosParticle*>::iterator pcl_itr = daughters.begin();

	//get all daughters and look for stable with same pgd id
	for(;pcl_itr != daughters.end();pcl_itr++)
	{
		if((*pcl_itr)->getPdgID()==this->getPdgID())
		return (*pcl_itr)->findLastSelf();
	}

	return this;
}

vector<PhotosParticle*> PhotosParticle::findProductionMothers()
{
	vector<PhotosParticle*> mothers = getMothers();
	vector<PhotosParticle*>::iterator pcl_itr = mothers.begin();

	//get all mothers and check none have pdg id of this one
	for(;pcl_itr != mothers.end();pcl_itr++)
	{
		if((*pcl_itr)->getPdgID()==this->getPdgID())
		return (*pcl_itr)->findProductionMothers();
	}
	return mothers;
}

vector<PhotosParticle *> PhotosParticle::getDecayTree()
{
	vector<PhotosParticle *> particles;
	particles.push_back(this);
	vector<PhotosParticle *> daughters = getDaughters();
	for(int i=0;i<(int)daughters.size();i++)
	{
		// Check if we are the first mother of each daughters
		// If not - skip this daughter
		PhotosParticle *p = daughters.at(i);
		vector<PhotosParticle *> mothers = p->getMothers();
		if(mothers.size()>1 && mothers.at(0)->getBarcode()!=getBarcode()) continue;
		vector<PhotosParticle *> tree = p->getDecayTree();
		particles.insert(particles.end(),tree.begin(),tree.end());
	}
	return particles;
}

void PhotosParticle::boostDaughtersFromRestFrame(PhotosParticle * tau_momentum)
{
	if(!hasDaughters()) //if there are no daughters
	return;

	// get all daughters, granddaughters, etc. then rotate and boost them
	vector<PhotosParticle*> list = getAllDecayProducts();
	vector<PhotosParticle*>::iterator pcl_itr = list.begin();

	for(;pcl_itr != list.end();pcl_itr++)
	{
		(*pcl_itr)->boostFromRestFrame(tau_momentum);
	}

	//checkMomentumConservation();
}

void PhotosParticle::boostDaughtersToRestFrame(PhotosParticle * tau_momentum)
{
	if(!hasDaughters()) //if there are no daughters
	return;

	// get all daughters, granddaughters, etc. then rotate and boost them
	vector<PhotosParticle*> list = getAllDecayProducts();
	vector<PhotosParticle*>::iterator pcl_itr = list.begin();

	for(;pcl_itr != list.end();pcl_itr++)
	{
		(*pcl_itr)->boostToRestFrame(tau_momentum);
	}

	//checkMomentumConservation();
}


void PhotosParticle::boostToRestFrame(PhotosParticle * tau_momentum)
{
	double theta = tau_momentum->getRotationAngle(Y_AXIS);
	tau_momentum->rotate(Y_AXIS,theta);

	double phi = tau_momentum->getRotationAngle(X_AXIS);
	tau_momentum->rotate(Y_AXIS,-theta);

	//Now rotate coordinates to get boost in Z direction.
	rotate(Y_AXIS,theta);
	rotate(X_AXIS,phi);
	boostAlongZ(-1*tau_momentum->getP(),tau_momentum->getE());
	rotate(X_AXIS,-phi);
	rotate(Y_AXIS,-theta);
}

void PhotosParticle::boostFromRestFrame(PhotosParticle * tau_momentum)
{
	//get the rotation angles
	//and boost z

	double theta = tau_momentum->getRotationAngle(Y_AXIS);
	tau_momentum->rotate(Y_AXIS,theta);

	double phi = tau_momentum->getRotationAngle(X_AXIS);
	tau_momentum->rotate(Y_AXIS,-theta);

	//Now rotate coordinates to get boost in Z direction.
	rotate(Y_AXIS,theta);
	rotate(X_AXIS,phi);
	boostAlongZ(tau_momentum->getP(),tau_momentum->getE());
	rotate(X_AXIS,-phi);
	rotate(Y_AXIS,-theta);
}

/** Get the angle needed to rotate the 4 momentum vector so that
    the x (y) component disapears. (and the Z component is > 0) */
double PhotosParticle::getRotationAngle(int axis, int second_axis)
{
	/**if(getP(axis)==0){
	if(getPz()>0)
	return 0; //no rotaion required
	else
	return M_PI;
	}**/
	if(getP(second_axis)==0)
	{
		if(getP(axis)>0) return -M_PI/2.0;
		else             return  M_PI/2.0;
	}
	if(getP(second_axis)>0) return     -atan(getP(axis)/getP(second_axis));
	else                    return M_PI-atan(getP(axis)/getP(second_axis));

}

/** Boost this vector along the Z direction.
    Assume no momentum components in the X or Y directions. */
void PhotosParticle::boostAlongZ(double boost_pz, double boost_e)
{
	// Boost along the Z axis
	double m_tau=sqrt(boost_e*boost_e-boost_pz*boost_pz);

	double p=getPz();
	double e=getE();

	setPz((boost_e*p + boost_pz*e)/m_tau);
	setE((boost_pz*p + boost_e*e )/m_tau);
}

/** Rotation around an axis X or Y */
void PhotosParticle::rotate(int axis,double theta, int second_axis)
{
	double temp_px=getP(axis);
	double temp_pz=getP(second_axis);
	setP(axis,cos(theta)*temp_px + sin(theta)*temp_pz);
	setP(second_axis,-sin(theta)*temp_px + cos(theta)*temp_pz);
}

void PhotosParticle::rotateDaughters(int axis,double theta, int second_axis)
{
	if(!hasDaughters()) //if there are no daughters
	return;

	vector<PhotosParticle*> daughters = getDaughters();
	vector<PhotosParticle*>::iterator pcl_itr = daughters.begin();

	//get all daughters then rotate and boost them.
	for(;pcl_itr != daughters.end();pcl_itr++)
	{
		(*pcl_itr)->rotate(axis,theta,second_axis);
		(*pcl_itr)->rotateDaughters(axis,theta,second_axis);
	}
	//checkMomentumConservation();
}

double PhotosParticle::getVirtuality()
{
	double e_sq=getE()*getE();
	double p_sq=getP()*getP();

	if(e_sq>p_sq) return    sqrt(e_sq-p_sq);
	else          return -1*sqrt(p_sq-e_sq); //if it's negative
}

double PhotosParticle::getP()
{
	return sqrt(getPx()*getPx()+getPy()*getPy()+getPz()*getPz());
}

double PhotosParticle::getP(int axis)
{
	if(axis==X_AXIS) return getPx();
	if(axis==Y_AXIS) return getPy();
	if(axis==Z_AXIS) return getPz();
	return 0;
}

void PhotosParticle::setP(int axis, double p_component)
{
	if(axis==X_AXIS) setPx(p_component);
	if(axis==Y_AXIS) setPy(p_component);
	if(axis==Z_AXIS) setPz(p_component);
}

} // namespace Photospp
