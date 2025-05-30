#pragma once

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <string>

namespace rad{
  namespace epic{
    using RVecS = ROOT::RVec<std::string>;
    using RVecU  = ROOT::RVecU;
    
    class PodioMetadata {

    public:
      PodioMetadata() = default;
      PodioMetadata(const std::string& filename, const std::string& treename="podio_metadata"){

	std::cout<<"PodioMetadata file : "  <<filename<<" "<<treename<<std::endl;
	ROOT::RDataFrame df(treename,filename);
	//Get the collectionId and names column data
	auto ids_ptr = df.Take<RVecU>("m_collectionIDs");
	auto names_ptr = df.Take<RVecS>("m_names");

	//Take the vectors of the first event
	//and keep them
	_collectionIDs = (*(ids_ptr ))[0];
	_names = (*(names_ptr ))[0];

	/* for(uint i=0;i<_names.size();++i){ */
	/*   cout<<_names[i]<<"\t collection id = "<<_collectionIDs[i]<<endl; */
	/* } */
	
      }

      bool Exists(const std::string name){
	if(std::find(_names.begin(),_names.end(),name)==_names.end()) return false;
	return true;
      }

      UInt_t CollectionIDFor(const std::string& name){
	//no checks made, use Exists first
	auto it = std::find(_names.begin(),_names.end(),name);
	UInt_t index =  it - _names.begin();
	return _collectionIDs[index];
      }
    private:
      
      RVecU _collectionIDs;
      RVecS _names;
      
    };

  }
}
