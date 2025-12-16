#pragma once

#include "StringUtilities.h"
#include <ROOT/RVec.hxx>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

namespace podio{
  namespace root_utils{
    // Struct definition from podio::root_utils
    struct CollectionWriteInfo {
      uint32_t collectionID{static_cast<uint32_t>(-1)};
      std::string dataType{};
      bool isSubset{false};
      unsigned int schemaVersion{0};
      std::string name{};
      std::string storageType{};
    };
    
  }
}

#pragma link C++ class podio::root_utils::CollectionWriteInfo+;
#pragma link C++ class vector<podio::root_utils::CollectionWriteInfo>+;

namespace rad{
  namespace podio{
    using RVecS = ROOT::RVec<std::string>;
    using RVecU  = ROOT::RVecU;
    
    class PodioMetadata {

    public:
      PodioMetadata() = default;
      PodioMetadata(const std::string& filename, const std::string& treename="podio_metadata"){
	
	std::cout<<"PodioMetadata file : "  <<filename<<" "<<treename<<std::endl;
	/* ROOT::RDataFrame df(treename,filename); */
	/* //Get the collectionId and names column data */
	/* auto ids_ptr = df.Take<RVecU>("m_collectionIDs"); */
	/* auto names_ptr = df.Take<RVecS>("m_names"); */

	/* //Take the vectors of the first event */
	/* //and keep them */
	/* _collectionIDs = (*(ids_ptr ))[0]; */
	/* _names = (*(names_ptr ))[0]; */
	
	/* for(uint i=0;i<_names.size();++i){ */
	/*   cout<<_names[i]<<"\t collection id = "<<_collectionIDs[i]<<endl; */
	/* } */
	
	// g. penman 17.11.25
	//new podio version fix
        
	// Open file
        auto file = std::unique_ptr<TFile>(TFile::Open(filename.c_str()));
        if (!file || file->IsZombie()) {
          throw std::runtime_error("Failed to open file: " + filename);
        }

        // Get metadata tree
        TTree *tree = file->Get<TTree>(treename.c_str());
        if (!tree) {
          throw std::runtime_error("Metadata tree not found: " + treename);
        }
	
        if (auto branch = tree->GetBranch("m_collectionIDs")){
	  // Old schema fallback
	  std::cout << "Branches exist, using old schema." << std::endl;
	  ROOT::RDataFrame df(treename, filename);
	  auto ids_ptr = df.Take<RVecU>("m_collectionIDs");
	  auto names_ptr = df.Take<RVecS>("m_names");
	  _collectionIDs = (*ids_ptr)[0];
	  _names = (*names_ptr)[0];
	} 
	else if (auto branch = tree->GetBranch("events___CollectionTypeInfo")) {
	  std::cout << "Branches dont exist, forwarding to new schema" << std::endl;
	  // e.g., vector<podio::root_utils::CollectionWriteInfo> or vector<tuple<...>>
	  std::string typeName = branch->GetClassName(); 
	  
	  if (typeName.find("CollectionWriteInfo") != std::string::npos) {
	    std::vector<::podio::root_utils::CollectionWriteInfo> *info = nullptr;
	    tree->SetBranchAddress("events___CollectionTypeInfo", &info);
	    tree->GetEntry(0);
	    for (auto &entry : *info) {
	      _collectionIDs.push_back(entry.collectionID);
	      _names.push_back(entry.name);
	    }
	  } 
	  else if (typeName.find("tuple") != std::string::npos) {
	    // Intermediate schema
	    std::cout << "Using old schema, new branches" << std::endl;
	    using TupleType = std::tuple<unsigned int, std::string, bool, unsigned int>;
	    std::vector<TupleType> *info = nullptr;
	    tree->SetBranchAddress("events___CollectionTypeInfo", &info);
	    tree->GetEntry(0);
	    for (auto &entry : *info) {
	      std::cout << std::get<0>(entry)  << std::endl;
	      std::cout << std::get<1>(entry) << std::endl;
	      _collectionIDs.push_back(std::get<0>(entry));
	      _names.push_back(std::get<1>(entry));
	    }
	  }else {
	    throw std::runtime_error("Unknown CollectionTypeInfo schema: " + typeName);
	  }
	}
      }
	
      bool Exists(const std::string &name) const {
	return std::find(_names.begin(), _names.end(), name) != _names.end();
      }

      UInt_t CollectionIDFor(const std::string& name){
	//no checks made, use Exists first
	auto it = std::find(_names.begin(),_names.end(),name);
	UInt_t index =  it - _names.begin();
	return _collectionIDs[index];
      }

      RVecS FilterNames(std::string sub_string){
	return rad::utils::filterStrings(_names,sub_string);
      }
    private:
      
      RVecU _collectionIDs;
      RVecS _names;
      
    };

  }
}
  
