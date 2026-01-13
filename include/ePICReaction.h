/**
 * @file ePICReaction.h
 * @brief Specialized Reaction class for ePIC RAD data analysis.
 */

#pragma once

#include "ElectroIonReaction.h"
#include "ParticleInjector.h"
#include "ePICSource.h" 
#include "ReactionUtilities.h"

namespace rad {
  namespace epic {
    
    using rad::consts::data_type::Rec;
    using rad::consts::data_type::Truth;
    using rad::Indices_t; 

    /**
     * @class ePICReaction
     * @brief High-level management of ePIC-specific data processing.
     * @details
     * Orchestrates the initialization of beams, the creation of unified 
     * reconstructed and truth vectors, and truth flagging metadata.
     * It merges multiple detector sources (Central, RP, ZDC) into a coherent event structure.
     */
    class ePICReaction : public ElectroIonReaction {
    public:
      /** * @brief Constructor for globbed filenames. 
       * @param treeName Name of the input TTree (e.g. "events").
       * @param fileNameGlob File pattern (e.g. "data/XXX.root").
       * @param columns Optional list of columns to read.
       */
      ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns ={} );
      
      /** * @brief Constructor for an existing RDataFrame. 
       * @param rdf Existing RDataFrame object.
       */
      ePICReaction(ROOT::RDataFrame rdf);

      // =================================================================================
      // Setup Methods
      // =================================================================================

      /** * @brief Sets up unified reconstructed vectors from detector sources.
       * @param isEnd If true, finalizes definitions (e.g. Beams). Default kTRUE.
       */
      void SetupReconstructed(Bool_t isEnd = kTRUE);

      /** * @brief Sets up truth vectors from MCParticles.
       * @param isEnd If true, finalizes definitions. Default kTRUE.
       */
      void SetupTruth(Bool_t isEnd = kTRUE);

      /** * @brief Sets up matching metadata between Rec and Truth.
       * @details Creates association maps and triggers reconstruction setup with matching enabled.
       * @param isEnd If true, finalizes definitions. Default kTRUE.
       */
      void SetupMatching(Bool_t isEnd = kTRUE);
      
      // =================================================================================
      // Functional Methods
      // =================================================================================

      /** * @brief Determines MC beam energies by scanning the header of the input file.
       * @details Automatically registers beam indices for both Rec and Truth streams.
       * @param ielIdx Index of the electron beam in MCParticles.
       * @param iionIdx Index of the ion beam in MCParticles.
       * @param nRows Number of rows to scan for energy averaging. Default 100.
       */
      void SetBeamsFromMC(UInt_t ielIdx, UInt_t iionIdx, Long64_t nRows = 100);

      /** * @brief Flags a candidate based on Truth ID.
       * @details Creates a column `[candidateName]_is_true` returning 1 if the truth PID matches.
       * @param candidateName Name of the particle candidate (e.g. "ele").
       * @param targetPid The PDG code to match against.
       */
      void MatchCandidateToTruth(const std::string& candidateName, int targetPid);

      // =================================================================================
      // Helper API (Public for ePICSource Policy)
      // =================================================================================

      /** * @brief Heuristic matching helper for forward detectors. 
       * @param prefix Column prefix for the detector (e.g. "rp_").
       * @param pidToMatch PDG ID to look for in MCParticles (heuristic).
       */
      void DefineForwardMatching(const std::string& prefix, int pidToMatch);

      /** * @brief Direct PID tagging helper from Truth data. 
       * @details Defines `[prefix]true_pid` by looking up the match ID in MCParticles.
       * @param prefix Column prefix for the detector.
       */
      void DefineTruePID(const std::string& prefix);

      // =================================================================================
      // Accessors
      // =================================================================================

      /** @return The beam electron 4-vector for crossing angle corrections. */
      const rad::epic::PxPyPzMVector& GetBeamElectronP4() const { return _p4el_beam; }
      
      /** @return The beam ion 4-vector for crossing angle corrections. */
      const rad::epic::PxPyPzMVector& GetBeamIonP4() const { return _p4ion_beam; }

    private:
      enum DetID { BEAM=0, CENTRAL=1, RP=2, ZDC=3, B0=4 };
      bool _truthMatched = false;   ///< Flag for match_id availability
      bool _beamsCorrected = false; ///< Flag for AfterBurner application
      Int_t _idxBeamEle = -1;       ///< Stored index for Beam Electron
      Int_t _idxBeamIon = -1;       ///< Stored index for Beam Ion
    }; 

    // =================================================================================
    // IMPLEMENTATION: ePICReaction
    // =================================================================================

    inline ePICReaction::ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns) 
      : ElectroIonReaction{treeName, fileNameGlob, columns} {}

    inline ePICReaction::ePICReaction(ROOT::RDataFrame rdf) 
      : ElectroIonReaction{rdf} {}

    inline void ePICReaction::SetupReconstructed(Bool_t isEnd) {
        AddType(Rec()); 
        rad::ParticleInjector injector(this);
        
        // Define unified vector structure
        std::vector<std::string> suffixes = {"double px", "double py", "double pz", "double m", "int pid", "short charge", "int det_id"};
        if(_truthMatched) {
            suffixes.push_back("int match_id");
            suffixes.push_back("int true_pid"); // Requires [Source]_true_pid to exist
        }
        injector.DefineParticleInfo(suffixes);

        // 1. Inject Beams (Scalar based, handled directly)
        // Automatically register these indices for analysis
        SetBeamElectronIndex(_idxBeamEle, Rec()); 
        SetBeamIonIndex(_idxBeamIon, Rec());
        
        DefineBeamComponents(); 
        // [Internal Beam Logic as discussed previously...]

        // 2. Define Detector Sources
        using Source = ePICSource<ePICReaction>;

        // Central Tracker: Apply crossing angle corrections
        Source central("ReconstructedParticles", "Central_", CENTRAL);
        central.SetIsCorrected(true);
        central.SetMinP(0.3); 
        central.Process(this, injector, _truthMatched);

        // Roman Pots: Forward matching to protons
        Source rp("ForwardRomanPotRecParticles", "rp_", RP);
        rp.SetTargetPID(2212);
        rp.SetMinP(1);
        rp.Process(this, injector, _truthMatched);

        // Zero Degree Calorimeter
        Source zdc("ReconstructedFarForwardZDCNeutrons", "ZDC_", ZDC);
        zdc.SetTargetPID(2112);
        zdc.SetIsCorrected(true);
        zdc.SetMinP(1);
        zdc.Process(this, injector, _truthMatched);

        // 3. Finalize Injection
        injector.CreateUnifiedVectors();
        util::CountParticles(this, Rec());
        if (isEnd) { DefineBeamElectron(); DefineBeamIon(); }
    }

    inline void ePICReaction::SetupTruth(Bool_t isEnd) {
        AddType(Truth());

        // Automatically register these indices for analysis
        SetBeamElectronIndex(_idxBeamEle, Truth()); 
        SetBeamIonIndex(_idxBeamIon, Truth());

        rad::ParticleInjector injector(this);
        injector.DefineParticleInfo({"double px", "double py", "double pz", "double m", "int pid", "int genStat", "int charge"});

        // Filter for final state particles (status 1) and beams (status 4)
        std::string filterGen = "MCParticles.generatorStatus==1||MCParticles.generatorStatus==4";
        injector.AddSource(Truth(), {
          "MCParticles.momentum.x", "MCParticles.momentum.y", "MCParticles.momentum.z",
          "MCParticles.mass", "MCParticles.PDG", "MCParticles.generatorStatus", "MCParticles.charge"
        }, filterGen); 

        injector.CreateUnifiedVectors();
        util::CountParticles(this, Truth());
        if (isEnd) { DefineBeamElectron(); DefineBeamIon(); }
    }

    inline void ePICReaction::SetupMatching(Bool_t isEnd) {
        _truthMatched = true;

        // Association Map lookup for Central Tracking
        Define("Central_match_id" + DoNotWriteTag(), 
            [](const ROOT::RVecU& recID, const ROOT::RVecU& simID, const Indices_t& rec_ind) {
                Indices_t match_id(rec_ind.size(), -1);
                for(size_t i = 0; i < recID.size(); ++i) {
                    if(recID[i] < (uint)match_id.size()) match_id[recID[i]] = (int)simID[i];
                }
                return match_id;
            }, 
            {"ReconstructedParticleAssociations.recID", "ReconstructedParticleAssociations.simID", "ReconstructedParticles.PDG"}
        );

        // First type setup call becomes GetDefaultType
        SetupReconstructed(kFALSE);
        SetupTruth(kFALSE);

        if (isEnd) {}
    }

    inline void ePICReaction::SetBeamsFromMC(UInt_t iel, UInt_t iion, Long64_t nRows) {
        _useBeamsFromMC = true;
        _idxBeamEle = iel;
        _idxBeamIon = iion;

        auto nthreads = ROOT::GetThreadPoolSize();
        if (nthreads) ROOT::DisableImplicitMT();

        auto tempframe = GetFileNames().empty() ? ROOT::RDataFrame{GetTreeName(), GetFileName()} : ROOT::RDataFrame{GetTreeName(), GetFileNames()};
        auto beamdf = tempframe.Range(nRows)
            .Define("emean", Form("MCParticles.momentum.z[%d]", iel))
            .Define("pzmean", Form("MCParticles.momentum.z[%d]", iion))
            .Define("pxmean", Form("MCParticles.momentum.x[%d]", iion));
          
        auto pze = beamdf.Mean("emean");
        auto pzp = beamdf.Mean("pzmean");
        auto pxp = beamdf.Mean("pxmean");

        // Set members directly to avoid premature Define() calls
        _p4el_beam.SetPxPyPzE(0, 0, *pze, std::abs(*pze));
        _p4ion_beam.SetPxPyPzE(*pxp, 0, *pzp, std::sqrt(*pxp * *pxp + *pzp * *pzp + 0.938*0.938));

        std::cout << " [ePICReaction] Beams Initialized: " << std::endl;
        std::cout << "    Electron: " << _p4el_beam.Pz() << " GeV" << std::endl;
        std::cout << "    Ion:       " << _p4ion_beam.Pz() << " GeV" << std::endl;
    
        if (nthreads) ROOT::EnableImplicitMT(nthreads);
    }

    inline void ePICReaction::MatchCandidateToTruth(const std::string& candidateName, int targetPid) {
        std::string idxCol = candidateName + "_idx";
        Define(candidateName + "_is_true", 
               [targetPid](const Indices_t& candIdx, const Indices_t& recMatchId, const Indices_t& rawPdg) {
                   if(candIdx.empty() || candIdx[0] == -1) return 0;
                   int rIdx = candIdx[0];
                   if(rIdx >= (int)recMatchId.size()) return 0;
                   int tIdx = recMatchId[rIdx];
                   if(tIdx == -1 || tIdx >= (int)rawPdg.size()) return 0;
                   return (std::abs(rawPdg[tIdx]) == targetPid) ? 1 : 0;
               }, 
               {idxCol, Rec() + "match_id", "MCParticles.PDG"});
    }

    inline void ePICReaction::DefineForwardMatching(const std::string& prefix, int pidToMatch) {
        Define(prefix + "match_id" + DoNotWriteTag(), 
            [pidToMatch](const Indices_t& pdg, const Indices_t& stat, const Indices_t& id_vec) {
                int best_idx = -1;
                for(size_t i=0; i<pdg.size(); ++i) {
                    if(std::abs(pdg[i]) == pidToMatch && stat[i] == 1) {
                        best_idx = i; break; 
                    }
                }
                auto result = Indices_t(id_vec.size(), best_idx);
                return result;
            }, {"MCParticles.PDG", "MCParticles.generatorStatus", prefix + "det_id" + DoNotWriteTag()});
    }

    inline void ePICReaction::DefineTruePID(const std::string& prefix) {
        Define(prefix + "true_pid" + DoNotWriteTag(), 
               [](const Indices_t& match_id, const Indices_t& tru_pdg) {
                   Indices_t tpid(match_id.size(), 0);
                   for(size_t i=0; i<(int)match_id.size(); ++i) {
                       if(match_id[i] != -1 && match_id[i] < (int)tru_pdg.size()) tpid[i] = tru_pdg[match_id[i]];
                   }
                   return tpid;
               }, {prefix + "match_id" + DoNotWriteTag(), "MCParticles.PDG"});
    }

  } // namespace epic
} // namespace rad
