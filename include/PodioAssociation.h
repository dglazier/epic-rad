/**
 * @file PodioAssociation.h
 * @brief Templates for handling One-To-Many PODIO associations in RDataFrame.
 */

#pragma once

#include <ROOT/RVec.hxx>
#include <vector>
#include <unordered_map>
#include "Constants.h"

namespace rad {
  namespace podio {

    using ROOT::RVecI;
    using ROOT::RVecU;

    /**
     * @class IdToLocalIndexMap
     * @brief Helper to map Global Collection IDs (from PODIO) to Local Vector Indices.
     * @details Used when merging multiple collections (e.g. EcalBarrel, EcalEndcap) 
     * to determine which sub-collection a specific hit belongs to.
     */
    class IdToLocalIndexMap {
    public:
      IdToLocalIndexMap(const RVecU& validCollIDs) {
        for (size_t i = 0; i < validCollIDs.size(); ++i) {
          _map[validCollIDs[i]] = static_cast<int>(i);
        }
      }

      RVecI operator()(const RVecU& inputIDs) const {
        RVecI result(inputIDs.size(), rad::consts::InvalidIndex());
        for (size_t i = 0; i < inputIDs.size(); ++i) {
          auto it = _map.find(inputIDs[i]);
          if (it != _map.end()) result[i] = it->second;
        }
        return result;
      }

    private:
      std::unordered_map<unsigned int, int> _map;
    };

    /**
     * @brief Merges data from multiple sub-detector collections into a unified flat vector.
     * @details Uses `local_collIds` to select which collection vector to read from,
     * and `indices` to select the element within that collection.
     */
    template <typename Tval>
    ROOT::RVec<Tval> MergeCollections(const ROOT::RVec<ROOT::RVec<Tval>>& collections,
                                      const RVecI& local_collIds,
                                      const RVecI& indices) {
      const size_t n = local_collIds.size();
      ROOT::RVec<Tval> result(n, rad::consts::InvalidEntry<Tval>());
      
      for (size_t i = 0; i < n; ++i) {
        int collIdx = local_collIds[i];
        if (collIdx != rad::consts::InvalidIndex() && collIdx < (int)collections.size()) {
           if(indices[i] < (int)collections[collIdx].size()) {
               result[i] = collections[collIdx][indices[i]];
           }
        }
      }
      return result;
    }

    /**
     * @brief Slices a flat unified vector into per-particle vectors (One-To-Many).
     * @param begin The start index for each particle's group.
     * @param end The end index (exclusive) for each particle's group.
     * @param flat_values The unified vector of associated object data.
     */
    template <typename Tval>
    ROOT::RVec<ROOT::RVec<Tval>> GroupOneToMany(const RVecI& begin,
                                                const RVecI& end,
                                                const ROOT::RVec<Tval>& flat_values) {
      const size_t n_particles = begin.size();
      ROOT::RVec<ROOT::RVec<Tval>> result;
      result.reserve(n_particles);

      for (size_t i = 0; i < n_particles; ++i) {
        // Construct RVec view/copy from the range [begin, end)
        if (end[i] >= begin[i] && end[i] <= (int)flat_values.size()) {
             auto start_it = flat_values.begin() + begin[i];
             auto end_it   = flat_values.begin() + end[i];
             result.emplace_back(start_it, end_it);
        } else {
             result.emplace_back(); // Empty vector for this particle
        }
      }
      return result;
    }

    /**
     * @brief High-level Wrapper: Merges collections AND groups them by particle.
     * @details This is the function called by RDataFrame to produce the final "rec_cluster_energy" column.
     */
    template <typename Tval>
    ROOT::RVec<ROOT::RVec<Tval>> CreateAssociation(const ROOT::RVec<ROOT::RVec<Tval>>& collections,
                                                   const RVecI& local_collIds,
                                                   const RVecI& indices,
                                                   const RVecI& rec_begin,
                                                   const RVecI& rec_end) {
      auto merged = MergeCollections(collections, local_collIds, indices);
      return GroupOneToMany(rec_begin, rec_end, merged);
    }

  } // namespace podio
} // namespace rad
