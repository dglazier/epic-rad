/**
 * @file ePICSource.h
 * @brief Configuration factory for ePIC detector collections.
 */

#pragma once

#include <string>
#include <vector>
#include <ROOT/RDataFrame.hxx>
#include "ParticleInjector.h"
#include "ePICUtilities.h" // For UndoAfterBurn

namespace rad {
  namespace epic {

    template <typename T> class ePICSource;

    /**
     * @class ePICSource
     * @brief A templated configuration factory for ePIC detector collections.
     * @tparam TReaction The reaction class type (e.g., ePICReaction).
     */
    template <typename TReaction>
    class ePICSource {
    public:
      /**
       * @brief Constructor
       * @param branch Raw branch name in the tree (e.g., "ReconstructedParticles")
       * @param prefix Internal prefix for this detector (e.g., "Central_")
       * @param detID  Integer ID assigned to this detector source
       */
      ePICSource(const std::string& branch, const std::string& prefix, int detID);

      /** @brief Set the PDG ID for truth matching (e.g. 2212 for Protons) */
      void SetTargetPID(int pid) { _targetPID = pid; }
      
      /** @brief Set a minimum momentum cut in GeV/c */
      void SetMinP(double p) { _minP = p; }
      
      /** @brief If true, use 'corr_px/y/z' columns instead of raw tree branches */
      void SetIsCorrected(bool b) { _isCorrected = b; }

      /**
       * @brief Orchestrates the column definitions and adds the source to the injector.
       */
      void Process(TReaction* reaction, rad::ParticleInjector& injector, bool truthMatched);

    private:
      /** @brief Internal helper to dispatch the crossing-angle correction. */
      void ApplyCorrection(TReaction* reaction);

      /** @brief Implementation of correction for specific types. */
      template <typename Tp, typename Tm>
      void ApplyCorrectionImpl(TReaction* reaction);

      std::string _branch;    ///< Raw EDM4hep branch name
      std::string _prefix;    ///< Internal RAD prefix
      int _detID=-1;             ///< Numerical ID for detector identification
      int _targetPID = 0;     ///< PID used for forward matching logic
      double _minP = 0.0;     ///< Minimum momentum magnitude for filtering
      bool _isCorrected = false; ///< True if AfterBurner corrections are applied
    };

    // =================================================================================
    // IMPLEMENTATION
    // =================================================================================

    template <typename TReaction>
    inline ePICSource<TReaction>::ePICSource(const std::string& branch, const std::string& prefix, int detID)
      : _branch(branch), _prefix(prefix), _detID(detID) {}

    template <typename TReaction>
    inline void ePICSource<TReaction>::Process(TReaction* reaction, rad::ParticleInjector& injector, bool truthMatched) {

      if (!reaction->OriginalColumnExists(_branch + ".momentum.x")) return;

      std::string dnw = DoNotWriteTag();
      
      // 1. Calculate crossing-angle corrected columns if requested
      if (_isCorrected) {
          ApplyCorrection(reaction);
      }

      // 2. Define Detector ID
      reaction->Define(_prefix + "det_id" + dnw, 
          Form("return ROOT::RVecI(%s.momentum.x.size(), %d);", _branch.c_str(), _detID));

      // 3. Setup Filter
      std::string filter = "";
      if (_minP > 0) {
          filter = Form("sqrt(%s.momentum.x*%s.momentum.x + %s.momentum.y*%s.momentum.y + %s.momentum.z*%s.momentum.z) > %f", 
                        _branch.c_str(), _branch.c_str(), _branch.c_str(), _branch.c_str(), _branch.c_str(), _branch.c_str(), _minP);
      }

      // 4. Collect Columns
      std::vector<std::string> cols;
      if (_isCorrected) {
          cols = {_prefix + "corr_px" + dnw, _prefix + "corr_py" + dnw, _prefix + "corr_pz" + dnw};
      } else {
          cols = {_branch + ".momentum.x", _branch + ".momentum.y", _branch + ".momentum.z"};
      }
      cols.push_back(_branch + ".mass");
      cols.push_back(_branch + ".PDG");

      if (reaction->OriginalColumnExists(_branch + ".charge")) {
          cols.push_back(_branch + ".charge");
      } else {
          reaction->Define(_prefix + "charge" + dnw, Form("return ROOT::RVec<short>(%s.momentum.x.size(), 0);", _branch.c_str()));
          cols.push_back(_prefix + "charge" + dnw);
      }
      cols.push_back(_prefix + "det_id" + dnw);

      if (truthMatched) {
          if (!reaction->ColumnExists(_prefix + "match_id" + dnw)) {
              reaction->DefineForwardMatching(_prefix, _targetPID);
          }
          cols.push_back(_prefix + "match_id" + dnw);
      }

      injector.AddSource(rad::consts::data_type::Rec(), cols, filter);
    }

    template <typename TReaction>
    inline void ePICSource<TReaction>::ApplyCorrection(TReaction* reaction) {
        // Use rad::DeduceColumnVectorType utility directly
        auto typeP = rad::DeduceColumnVectorType(reaction, _branch + ".momentum.x");
        auto typeM = rad::DeduceColumnVectorType(reaction, _branch + ".mass");

        if (typeP == ColType::Float && typeM == ColType::Float) ApplyCorrectionImpl<float, float>(reaction);
        else if (typeP == ColType::Double && typeM == ColType::Double) ApplyCorrectionImpl<double, double>(reaction);
        else if (typeP == ColType::Float && typeM == ColType::Double) ApplyCorrectionImpl<float, double>(reaction);
        else if (typeP == ColType::Double && typeM == ColType::Float) ApplyCorrectionImpl<double, float>(reaction);
    }

    template <typename TReaction>
    template <typename Tp, typename Tm>
    inline void ePICSource<TReaction>::ApplyCorrectionImpl(TReaction* reaction) {
        // Get beams (now returning PxPyPzMVector)
        rad::epic::UndoAfterBurn<Tp, Tm> transformer(reaction->GetBeamIonP4(), reaction->GetBeamElectronP4());
        
        std::string dnw = DoNotWriteTag();
        std::string vecName = _prefix + "corr_vec" + dnw;
        
        reaction->Define(vecName, 
            [transformer](const ROOT::RVec<Tp>& px, const ROOT::RVec<Tp>& py, const ROOT::RVec<Tp>& pz, const ROOT::RVec<Tm>& m){
                return transformer(px, py, pz, m);
            }, 
            {_branch + ".momentum.x", _branch + ".momentum.y", _branch + ".momentum.z", _branch + ".mass"});

        reaction->Define(_prefix + "corr_px" + dnw, [](const ROOT::RVec<rad::epic::PxPyPzEVector>& v){ 
            ROOT::RVecD out(v.size()); for(size_t i=0;i<v.size();++i) out[i] = v[i].Px(); return out; 
        }, {vecName});
        reaction->Define(_prefix + "corr_py" + dnw, [](const ROOT::RVec<rad::epic::PxPyPzEVector>& v){ 
            ROOT::RVecD out(v.size()); for(size_t i=0;i<v.size();++i) out[i] = v[i].Py(); return out; 
        }, {vecName});
        reaction->Define(_prefix + "corr_pz" + dnw, [](const ROOT::RVec<rad::epic::PxPyPzEVector>& v){ 
            ROOT::RVecD out(v.size()); for(size_t i=0;i<v.size();++i) out[i] = v[i].Pz(); return out; 
        }, {vecName});
    }

  } // namespace epic
} // namespace rad
