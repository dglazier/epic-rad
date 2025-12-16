#pragma once
/**
 * Derived class of ParticleModifier to add ePIC specific utilities
 */

#include "ParticleModifier.h"
#include <TMath.h>

namespace rad{
  namespace epic{


    ///////////////////////////////////////////////////////////////
    ///Class Definition
    ///////////////////////////////////////////////////////////////
    class ePICParticleModifier : public rad::config::ParticleModifier {
      
    public:
      
    ePICParticleModifier(rad::config::ConfigReaction& cr): ParticleModifier{cr}{};
      
      void UndoCrossAngle(const std::string& particle){
	std::vector<std::string> vars = {particle, Rec()+"px", Rec()+"py", Rec()+"pz", Rec()+"m"};
	auto func_call_str = rad::utils::createFunctionCallStringFromVec("rad::epic::UndoCrossAnglePart", vars);
	//_modifications.push_back(func_call_str);
	addModification(func_call_str);
	
      }

    };
     
      template <typename Tp, typename Tm>
	bool UndoCrossAnglePart(const int idx, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m){
	    
	auto cross_angle = 25e-3;
	    auto vec = PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]);
	    RotationY UndoAB(cross_angle);
	    auto new_vec = UndoAB(vec);
	    px[idx] = new_vec.X();
	    py[idx] = new_vec.Y();
	    pz[idx] = new_vec.Z();
	    
	    /* if(!std::isnan(vec.X())){ */
	    /*   cout << "Before: " << vec << endl; */
	    /*   cout << "After: " << new_vec << endl; */
	    /* } */
	    
	    return true;
      }
    



    
  }//epic
}//rad
