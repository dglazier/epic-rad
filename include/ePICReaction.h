/**
 * @file ePICReaction.h
 * @brief Unified Reaction class for ePIC analysis (Physics + Detector).
 */

#pragma once

#include "ElectroIonReaction.h"
#include "ParticleInjector.h"
#include "ePICSource.h" 
#include "ReactionUtilities.h"
#include "PodioMetadata.h"
#include "PodioAssociation.h"
#include <TChain.h>
#include <memory>
#include <regex> 

namespace rad {
  namespace epic {
    
    using rad::consts::data_type::Rec;
    using rad::consts::data_type::Truth;
    using rad::Indices_t; 
    using ROOT::RVecU;
    enum DetID { BEAM=0, CENTRAL=1, RP=2, ZDC=3, B0=4 };
 
    /**
     * @class ePICReaction
     * @brief High-level management of ePIC-specific data processing.
     * @details
     * This class acts as the central coordinator for ePIC analysis. It bridges two worlds:
     * 1. **Physics Reconstruction:** Inherits from `ElectroIonReaction` to handle Beams, Q2, x, and Kinematics.
     * 2. **Detector Associations:** Manages PODIO metadata to link high-level Reconstructed Particles
     * back to their low-level detector objects (Clusters, Tracks, Hits).
     */
    class ePICReaction : public ElectroIonReaction {
    public:
      /** * @brief Constructor for globbed filenames. 
       * @details Automatically loads PODIO metadata from the first file in the glob pattern.
       * @param treeName Name of the input TTree (e.g. "events").
       * @param fileNameGlob File pattern (e.g. "data/XXX.root").
       * @param columns Optional list of columns to read (optimization).
       */
      ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns ={} );
      
      /** * @brief Constructor for a vector of filenames. 
       * @param treeName Name of the input TTree.
       * @param filenames Vector of file paths.
       * @param columns Optional list of columns to read.
       */
      ePICReaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns ={} );
      
      /** * @brief Constructor for existing RDataFrame. 
       * @note Metadata loading is skipped here; Detector Associations will throw an error if used.
       * @param rdf Existing RDataFrame object.
       */
      ePICReaction(ROOT::RDataFrame rdf);

      /** * @brief Constructor for existing RDataFrame. 
       * @note Metadata loading is skipped here; Detector Associations will throw an error if used.
       * @param rdf Existing RDataFrame object.
       */
      ePICReaction(ROOT::RDF::RNode rdf);

      // =================================================================================
      // Setup Methods
      // =================================================================================

      /** * @brief Sets up unified Reconstructed vectors (px, py, pz, m, pid...). 
       * @param isEnd If true, finalize definitions (Beam Particles). Default kTRUE.
       */
      void SetupReconstructed(Bool_t isEnd = kTRUE);

      /** * @brief Sets up unified Truth vectors from MCParticles. 
       * @param isEnd If true, finalize definitions. Default kTRUE.
       */
      void SetupTruth(Bool_t isEnd = kTRUE);

      /** * @brief Sets up matching metadata between Rec and Truth.
       * @details Creates the `rec_match_id` column linking Reconstructed tracks to MCParticles.
       * @param isEnd If true, finalize definitions. Default kTRUE.
       */
      void SetupMatching(Bool_t isEnd = kTRUE);
      
      // =================================================================================
      // Detector Association API
      // =================================================================================
      /**
       * @brief Creates a 1/0 flag aligned with ReconstructedParticles indicating 
       * if the particle was also reconstructed in a specific subdetector.
       * @param outName The name of the new column (e.g., "rec_has_tagger")
       * @param subdetAssocName The association collection of the subdetector.
       */
      void DefineDetectorFlag(const std::string& outName, const std::string& subdetAssocName);
      
      // =================================================================================
      // Functional Methods
      // =================================================================================

      /** * @brief Scans the input file header to determine Beam Energies. 
       * @param ielIdx Index of the electron beam in MCParticles.
       * @param iionIdx Index of the ion beam in MCParticles.
       * @param nRows Number of events to average over (default 100).
       */
      void SetBeamsFromMC(UInt_t ielIdx, UInt_t iionIdx, Long64_t nRows = 100);

      /** * @brief Flags a candidate based on Truth Matching.
       * @details Creates `[candidateName]_is_true` column. Returns 1 if the candidate's
       * matched truth particle has the specified PDG ID.
       * @param candidateName Name of the particle (e.g. "scat_ele").
       * @param targetPid PDG ID to match against (e.g. 11).
       */
      void MatchCandidateToTruth(const std::string& candidateName, int targetPid);

      /** @brief Heuristic matching for Forward Detectors (RP, ZDC). */
      void DefineForwardMatching(const std::string& prefix, int pidToMatch);
      
      /** @brief Defines [prefix]true_pid column based on [prefix]match_id. */
      void DefineTruePID(const std::string& prefix);
      
      const rad::epic::PxPyPzMVector& GetBeamElectronP4() const { return _p4el_beam; }
      const rad::epic::PxPyPzMVector& GetBeamIonP4() const { return _p4ion_beam; }
      
      /** @return Raw pointer to metadata manager (nullable). */
      rad::podio::PodioMetadata* GetMetadata() { return _podioMetadata.get(); }

      void SetIonPdg(const int pdg){_ionPDG=pdg;}

      struct DetectorDef {
        std::string branch;
        std::string prefix;
        int id;
    };

    // New method required by ePICAssociationManager
    DetectorDef GetDetector(const std::string& name) const {
        if (_detectors.find(name) == _detectors.end()) {
            throw std::invalid_argument("Detector " + name + " not registered in ePICReaction.");
        }
        return _detectors.at(name);
    }

    // Optional: Let users register custom detectors dynamically!
    void RegisterDetector(const std::string& name, const std::string& branch, const std::string& prefix, int id) {
        _detectors[name] = {branch, prefix, id};
    }
      /** @brief Extracts inner type string from "RVec<T>". Needed for template metaprogramming. */
      std::string GetElementTypeName(const std::string& colName);

    private:
      /** @brief internal helper to load metadata safely. */
      void InitMetadata(const std::string& filename);
      
  
      bool _truthMatched = false;   
      bool _beamsCorrected = false; 
      Int_t _idxBeamEle = 0;       
      Int_t _idxBeamIon = 1;
      Int_t _ionPDG=2212;
      
      std::map<std::string, DetectorDef> _detectors;
      std::shared_ptr<rad::podio::PodioMetadata> _podioMetadata;
    }; 

    // =================================================================================
    // IMPLEMENTATION
    // =================================================================================

    inline ePICReaction::ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns) 
      : ElectroIonReaction{treeName, fileNameGlob, columns} 
    {
        // Resolve Glob to single file for metadata loading
        TChain chain("podio_metadata"); 
        chain.Add(fileNameGlob.data());
        // Only load metadata if we found files
        if(chain.GetListOfFiles() && chain.GetListOfFiles()->GetEntries() > 0) {
            InitMetadata(chain.GetListOfFiles()->At(0)->GetTitle());
        }
    }

    inline ePICReaction::ePICReaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns) 
      : ElectroIonReaction{treeName, filenames, columns} 
    {
        if(!filenames.empty()) InitMetadata(filenames[0]);
    }

    inline ePICReaction::ePICReaction(ROOT::RDataFrame rdf) : ElectroIonReaction{rdf} {}

    inline ePICReaction::ePICReaction(ROOT::RDF::RNode rdf) : ElectroIonReaction{rdf} {}

    inline void ePICReaction::InitMetadata(const std::string& filename) {
        try { _podioMetadata = std::make_shared<rad::podio::PodioMetadata>(filename); } 
        catch (const std::exception& e) { std::cerr << "[ePICReaction] Warning: Failed to load metadata: " << e.what() << std::endl; }
    }

    inline std::string ePICReaction::GetElementTypeName(const std::string& colName) {
         // Get the full type string (e.g. "ROOT::VecOps::RVec<float>")
         std::string fullType = ColObjTypeString(colName);
         
         // Extract content between < and >
         auto start = fullType.find('<');
         auto end   = fullType.rfind('>');
         if(start != std::string::npos && end != std::string::npos) {
             return fullType.substr(start + 1, end - start - 1);
         }
         return fullType; // Fallback (shouldn't happen for PODIO vectors)
    }

   
    // --- Standard Methods ---
 /** * @brief Sets up unified Reconstructed vectors and populates the detector registry.
     * @details 
     * This method initializes the `ParticleInjector` to build the unified kinematics arrays 
     * (e.g., `rec_px`, `rec_det_id`). It registers the active detectors (Central, RP, ZDC) 
     * into the internal `_detectors` map. The `ePICAssociationManager` later relies on this 
     * map to correctly pad auxiliary detector data (like ECal clusters) into the unified coordinate space.
     * * @param isEnd If true, finalize definitions. Default kTRUE.
     */
    inline void ePICReaction::SetupReconstructed(Bool_t isEnd) {
        // 1. Initialize data type and base kinematic columns
        AddType(Rec());

        DefineBeamComponents(Rec());
        SetBeamElectronIndex(_idxBeamEle, Rec());
        SetBeamIonIndex(_idxBeamIon, Rec());

        // 2. Initialize Particle Injector
        rad::ParticleInjector injector(this);
        
        // Define standard particle columns (Structure of Arrays)
        ROOT::RVec<std::string> suffixes = {
            "double px", "double py", "double pz", "double m", 
            "int pid", "short charge", "int det_id"
        };
        if (_truthMatched) { suffixes.push_back("int match_id"); }
        injector.DefineParticleInfo(suffixes);

        // 3. Inject Beams First (Indices 0 and 1, det_id = BEAM = 0)
        std::string bEle = Rec() + consts::BeamEle() + "_src_";
        std::string bIon = Rec() + consts::BeamIon() + "_src_";

        ROOT::RVec<std::string> ele_cols = {
            bEle + consts::NamePx(), bEle + consts::NamePy(), bEle + consts::NamePz(), 
            bEle + consts::NameM(),  bEle + consts::NamePid(), 
            "ROOT::RVecI{-1}", "ROOT::RVecI{0}" // Charge=-1 (Electron), DetID=0
        };

        ROOT::RVec<std::string> ion_cols = {
            bIon + consts::NamePx(), bIon + consts::NamePy(), bIon + consts::NamePz(), 
            bIon + consts::NameM(),  bIon + consts::NamePid(), 
            "ROOT::RVecI{1}", "ROOT::RVecI{0}" // Charge=1 (Proton/Ion), DetID=0
        };

        // Handle Truth Matching Columns for Beams
        if (_truthMatched) {
            ele_cols.push_back("rad::Indices_t{0}"); // Electron Beam: Matches Truth Index 0
            ion_cols.push_back("rad::Indices_t{1}"); // Ion Beam: Matches Truth Index 1
        }

        // Add Beams to Injector
        injector.AddSource(Rec(), ele_cols);
        injector.AddSource(Rec(), ion_cols);
 
        // 4. Register Detectors to the Internal Map
        // This decouples PODIO associations from the Reaction class.
        RegisterDetector("Central", "ReconstructedParticles", "Central_", CENTRAL);
        RegisterDetector("RP", "ForwardRomanPotRecParticles", "rp_", RP);
        RegisterDetector("ZDC", "ReconstructedFarForwardZDCNeutrons", "ZDC_", ZDC);

        // 5. Inject Detectors via ePICSource Factory
        using Source = ePICSource<ePICReaction>;
        
        // --- Central Tracker ---
        Source central(_detectors.at("Central").branch, _detectors.at("Central").prefix, _detectors.at("Central").id);
        central.SetIsCorrected(true); 
        central.SetMinP(0.3); 
        central.Process(this, injector, _truthMatched);
        
        // --- Roman Pots ---
        Source rp(_detectors.at("RP").branch, _detectors.at("RP").prefix, _detectors.at("RP").id);
        rp.SetTargetPID(2212); 
        rp.SetMinP(1); 
        rp.Process(this, injector, _truthMatched);
        
        // --- ZDC ---
        Source zdc(_detectors.at("ZDC").branch, _detectors.at("ZDC").prefix, _detectors.at("ZDC").id);
        zdc.SetTargetPID(2112); 
        zdc.SetIsCorrected(true); 
        zdc.SetMinP(1); 
        zdc.Process(this, injector, _truthMatched);

        // 6. Finalize Structure of Arrays
        injector.CreateUnifiedVectors();
        
        util::CountParticles(this, Rec());
        
        if (isEnd) {
            // Placeholder for end-of-setup logic
        }
    }
    
inline void ePICReaction::SetupTruth(Bool_t isEnd) {
    AddType(Truth());
    
    DefineBeamComponents(Truth());
    SetBeamElectronIndex(_idxBeamEle, Truth());
    SetBeamIonIndex(_idxBeamIon, Truth());

    rad::ParticleInjector injector(this);
    injector.DefineParticleInfo({"double px", "double py", "double pz", "double m", "int pid", "int genStat", "int charge"});
    
    std::string bEle = Truth() + consts::BeamEle() + "_src_";
    std::string bIon = Truth() + consts::BeamIon() + "_src_";
    
    injector.AddSource(Truth(), {
        bEle + consts::NamePx(), bEle + consts::NamePy(), bEle + consts::NamePz(), 
        bEle + consts::NameM(),  bEle + consts::NamePid(), 
        "ROOT::RVecI{4}", "ROOT::RVecI{-1}" 
    });

    injector.AddSource(Truth(), {
        bIon + consts::NamePx(), bIon + consts::NamePy(), bIon + consts::NamePz(), 
        bIon + consts::NameM(),  bIon + consts::NamePid(), 
        "ROOT::RVecI{4}", "ROOT::RVecI{1}" 
    });
    
    rad::epic::ApplyCrossingAngleCorrection(
        this, Truth(),
        "MCParticles.momentum.x", "MCParticles.momentum.y", "MCParticles.momentum.z", "MCParticles.mass",
        GetBeamIonP4(), GetBeamElectronP4()
    );

    std::string dnw = DoNotWriteTag();
    injector.AddSource(Truth(), {
      Truth() + "corr_px" + dnw, Truth() + "corr_py" + dnw, Truth() + "corr_pz" + dnw,
      "MCParticles.mass", "MCParticles.PDG", "MCParticles.generatorStatus", "MCParticles.charge"
    }, "MCParticles.generatorStatus>0 && MCParticles.generatorStatus!=4");
    
    injector.CreateUnifiedVectors();
    util::CountParticles(this, Truth());
    
    if (isEnd) { }
} 

    // inline void ePICReaction::SetupTruth(Bool_t isEnd) {
    //     AddType(Truth());
	
    // 	// 1. DEFINE BEAM COLUMNS (Truth Version)
    //     DefineBeamComponents(Truth());
    // 	SetBeamElectronIndex(_idxBeamEle,Truth());
    // 	SetBeamIonIndex(_idxBeamIon,Truth());


    //     rad::ParticleInjector injector(this);
    //     injector.DefineParticleInfo({"double px", "double py", "double pz", "double m", "int pid", "int genStat", "int charge"});
	
    // 	// INJECT BEAMS FIRST
    // 	std::string bEle = Truth() + consts::BeamEle() + "_src_";
    //     std::string bIon = Truth() + consts::BeamIon() + "_src_";
	
    //      // Electron Beam (Index 0)
    //     injector.AddSource(Truth(), {
    //         bEle + consts::NamePx(), bEle + consts::NamePy(), bEle + consts::NamePz(), 
    //         bEle + consts::NameM(),  bEle + consts::NamePid(), 
    //         "ROOT::RVecI{4}", "ROOT::RVecI{-1}" // Status 4 (Beam), Charge -1
    //     });

    //     // Ion Beam (Index 1)
    //     injector.AddSource(Truth(), {
    //         bIon + consts::NamePx(), bIon + consts::NamePy(), bIon + consts::NamePz(), 
    //         bIon + consts::NameM(),  bIon + consts::NamePid(), 
    //         "ROOT::RVecI{4}", "ROOT::RVecI{1}" // Status 4 (Beam), Charge +1
    //     });
	
    //     // Filter: Status 1 (Stable) 
    //     injector.AddSource(Truth(), {
    //       "MCParticles.momentum.x", "MCParticles.momentum.y", "MCParticles.momentum.z",
    //       "MCParticles.mass", "MCParticles.PDG", "MCParticles.generatorStatus", "MCParticles.charge"
    //     }, "MCParticles.generatorStatus>0 && MCParticles.generatorStatus!=4"); 
       
    //     injector.CreateUnifiedVectors();
    //     util::CountParticles(this, Truth());
        
    // 	if (isEnd) {

    // 	}
    // }

    inline void ePICReaction::SetupMatching(Bool_t isEnd) {
      _truthMatched = true;

     
      // Define central matching ID column
      // Uses the ReconstructedParticleAssociations table to map SimID -> filtered SimID
      Define("Central_match_id" + DoNotWriteTag(),
      [ionPDG=_ionPDG](const Indices_t& recID, const Indices_t& simID, const Indices_t& rec_ind, const Indices_t& genStat, const ROOT::RVecI& pid){
        
        const auto nTru0 = genStat.size(); 
        const auto nTru = ROOT::VecOps::Sum(genStat > 0); 
        Indices_t map(nTru0, -1); 

        const auto beam_e_sim = ROOT::VecOps::ArgMax(genStat == 4 && pid == 11);
        const auto beam_ion_sim = ROOT::VecOps::ArgMax(genStat == 4 && pid == ionPDG);

        map[beam_e_sim]  = 0;
        map[beam_ion_sim] = 1;

        uint pos = 2; 
        for (size_t i = 0; i < nTru0; ++i) {
          if(genStat[i] > 0 && genStat[i] != 4){ 
            map[i] = pos++;
          }
        }
        
        const int num_reco_particles = static_cast<int>(rec_ind.size()); 
        Indices_t match_id(num_reco_particles, -1);
        
        const int num_associations = static_cast<int>(recID.size()); 
        const int nMap = static_cast<int>(map.size()); 
        
        for(Int_t i = 0; i < num_associations; ++i) {
          const int r = recID[i];
          const int s = simID[i];
          if (r >= 0 && r < num_reco_particles && s >= 0 && s < nMap) {
            match_id[r] = map[s]; 
          }
        }
        return match_id;
      }, 
      {"_ReconstructedParticleAssociations_rec.index", "_ReconstructedParticleAssociations_sim.index", "ReconstructedParticles.PDG", "MCParticles.generatorStatus", "MCParticles.PDG"});
      
    
      
      SetupReconstructed(kFALSE); SetupTruth(kFALSE); DefineTruePID(Rec());
    }
  

    inline void ePICReaction::SetBeamsFromMC(UInt_t iel, UInt_t iion, Long64_t nRows) {
        
        // Heuristic: Scan first nRows to find average beam energy
        auto nthreads = ROOT::GetThreadPoolSize();
        if (nthreads) ROOT::DisableImplicitMT(); // Disable MT for simple range scan
        
        auto tempframe = GetFileNames().empty() ? ROOT::RDataFrame{GetTreeName(), GetFileName()} : ROOT::RDataFrame{GetTreeName(), utils::as_stdvector(GetFileNames())};
        auto beamdf = tempframe.Range(nRows).Define("emean", Form("MCParticles.momentum.z[%d]", iel)).Define("pzmean", Form("MCParticles.momentum.z[%d]", iion)).Define("pxmean", Form("MCParticles.momentum.x[%d]", iion));
        
        auto pze = beamdf.Mean("emean"); 
        auto pzp = beamdf.Mean("pzmean"); 
        auto pxp = beamdf.Mean("pxmean");

        // Set internal beam vectors (used for crossing angle correction in ePICSource)
        _p4el_beam.SetPxPyPzE(0, 0, *pze, std::abs(*pze));
        _p4ion_beam.SetPxPyPzE(*pxp, 0, *pzp, std::sqrt(*pxp * *pxp + *pzp * *pzp + 0.938*0.938));
        
        std::cout << " [ePICReaction] Beams: Ele=" << _p4el_beam.Pz() << " GeV, Ion=" << _p4ion_beam.Pz() << " GeV" << std::endl;
        if (nthreads) ROOT::EnableImplicitMT(nthreads);
    }

    inline void ePICReaction::MatchCandidateToTruth(const std::string& candidateName, int targetPid) {
        std::string idxCol = candidateName + "_idx";
        Define(candidateName + "_is_true", 
               [targetPid](const Indices_t& candIdx, const Indices_t& recMatchId, const Indices_t& rawPdg) {
                   if(candIdx.empty() || candIdx[0] == -1) return 0;
                   int rIdx = candIdx[0]; if(rIdx >= (int)recMatchId.size()) return 0;
                   int tIdx = recMatchId[rIdx]; if(tIdx == -1 || tIdx >= (int)rawPdg.size()) return 0;
                   return (std::abs(rawPdg[tIdx]) == targetPid) ? 1 : 0;
               }, {idxCol, Rec() + "match_id", "MCParticles.PDG"});
    }

    inline void ePICReaction::DefineForwardMatching(const std::string& prefix, int pidToMatch) {
        // Heuristic: Match simply by PID for forward detectors where associations might be missing
        Define(prefix + "match_id" + DoNotWriteTag(), 
            [pidToMatch](const Indices_t& pdg, const Indices_t& stat, const Indices_t& id_vec) {
                int best_idx = -1;
                for(size_t i=0; i<pdg.size(); ++i) { if(std::abs(pdg[i]) == pidToMatch && stat[i] == 1) { best_idx = i; break; } }
                auto result = Indices_t(id_vec.size(), best_idx); return result;
            }, {"MCParticles.PDG", "MCParticles.generatorStatus", prefix + "det_id" + DoNotWriteTag()});
    }

    inline void ePICReaction::DefineTruePID(const std::string& prefix) {
      rad::PrintDefinedColumnNames(CurrFrame());

        // Helper to pull Truth PID into the Reconstructed stream for convenience
        Define(prefix + "true_pid", 
               [](const Indices_t& match_id, const Indices_t& tru_pdg) {
                   Indices_t tpid(match_id.size(), 0);
                   for(UInt_t i=0; i<match_id.size(); ++i) { if(match_id[i] != -1 && match_id[i] < (int)tru_pdg.size()) tpid[i] = tru_pdg[match_id[i]]; }
                   return tpid;
              }, {prefix + "match_id", "tru_pid"});
    }
    inline void ePICReaction::DefineDetectorFlag(const std::string& outName, const std::string& subdetAssocName) {
        
        // The raw simulated indices seen by the subdetector
        std::string subdetSimCol = "_" + subdetAssocName + "_sim.index";
        
        Define(outName, 
            [ionPDG=_ionPDG](const rad::Indices_t& rec_match_id, 
                             const rad::Indices_t& subdet_sim_idx,
                             const rad::Indices_t& genStat,
                             const ROOT::RVecI& pid) {
                   
                // 1. Initialize output vector perfectly aligned with rec_px / rec_match_id
                ROOT::RVecI result(rec_match_id.size(), 0);
                if (subdet_sim_idx.empty() || rec_match_id.empty()) return result;

                // 2. Recreate the Truth Map (Raw MC Index -> Filtered MC Index)
                // This guarantees we are speaking the same "language" as rec_match_id
                const auto nTru0 = genStat.size();
                rad::Indices_t map(nTru0, -1); 
                
                const auto beam_e_sim = ROOT::VecOps::ArgMax(genStat == 4 && pid == 11);
                const auto beam_ion_sim = ROOT::VecOps::ArgMax(genStat == 4 && pid == ionPDG);
                
                map[beam_e_sim] = 0;
                map[beam_ion_sim] = 1;
                
                uint pos = 2;
                for (size_t i = 0; i < nTru0; ++i) {
                    if(genStat[i] > 0 && genStat[i] != 4){ 
                        map[i] = pos++;
                    }
                }

                // 3. Map the subdetector's raw truth hits to the filtered truth IDs
                int max_mapped_sim = -1;
                for (size_t raw_idx : subdet_sim_idx) {
                    if (raw_idx < nTru0) {
                        int mapped_idx = map[raw_idx];
                        if (mapped_idx > max_mapped_sim) max_mapped_sim = mapped_idx;
                    }
                }

                // If the subdetector saw no valid truth particles, return 0s
                if (max_mapped_sim < 0) return result;

                // Create the fast O(1) lookup table for mapped truth indices
                std::vector<bool> subdet_saw_truth(max_mapped_sim + 1, false);
                for (size_t raw_idx : subdet_sim_idx) {
                    if (raw_idx < nTru0) {
                        int mapped_idx = map[raw_idx];
                        if (mapped_idx >= 0) subdet_saw_truth[mapped_idx] = true;
                    }
                }

                // 4. Loop over the unified rec_match_id and flag matches
                for(size_t i = 0; i < rec_match_id.size(); ++i) {
                    int m_idx = (int)rec_match_id[i];
                    
                    // If this unified particle has a valid truth match, 
                    // and the subdetector ALSO saw that same truth particle, flag it!
                    if (m_idx >= 0 && m_idx <= max_mapped_sim) {
                        if (subdet_saw_truth[m_idx]) {
                            result[i] = 1; 
                        }
                    }
                }
                
                return result;
            }, 
            {Rec() + "match_id", // The unified match ID synced with rec_px
             subdetSimCol, 
             "MCParticles.generatorStatus", 
             "MCParticles.PDG"}
        );
    }
    
  } 
}

// /**
//  * @file ePICReaction.h
//  * @brief Unified Reaction class for ePIC analysis (Physics + Detector).
//  */

// #pragma once

// #include "ElectroIonReaction.h"
// #include "ParticleInjector.h"
// #include "ePICSource.h" 
// #include "ReactionUtilities.h"
// #include "PodioMetadata.h"
// #include "PodioAssociation.h"
// #include <TChain.h>
// #include <memory>

// namespace rad {
//   namespace epic {
    
//     using rad::consts::data_type::Rec;
//     using rad::consts::data_type::Truth;
//     using rad::Indices_t; 
//     using ROOT::RVecU;

//     /**
//      * @class ePICReaction
//      * @brief High-level management of ePIC-specific data processing.
//      * @details
//      * Merges physics object reconstruction (Tracks/Beams) with detector-level 
//      * association logic (Clusters/Hits).
//      */
//     class ePICReaction : public ElectroIonReaction {
//     public:
//       /** * @brief Constructor for globbed filenames. */
//       ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns ={} );
      
//       /** * @brief Constructor for a vector of filenames. */
//       ePICReaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns ={} );
      
//       /** * @brief Constructor for existing RDataFrame (No Metadata support). */
//       ePICReaction(ROOT::RDataFrame rdf);

//       // =================================================================================
//       // Setup Methods
//       // =================================================================================

//       void SetupReconstructed(Bool_t isEnd = kTRUE);
//       void SetupTruth(Bool_t isEnd = kTRUE);
//       void SetupMatching(Bool_t isEnd = kTRUE);
      
//       // =================================================================================
//       // Detector Association API
//       // =================================================================================

//       /**
//        * @brief Creates a unified column linking Rec Particles to Detector Objects.
//        * @details 
//        * Creates a column named `[prefix]_[objName]_[member]` (e.g. "rec_cluster_energy").
//        * The output is of type `RVec<RVec<T>>` (One-to-Many).
//        * * @param objName The object name (e.g. "clusters").
//        * @param collectionNames List of valid PODIO collection names (e.g. "EcalBarrelClusters").
//        * @param memberName The leaf to extract (e.g. "energy").
//        */
//       void DefineAssociation(const std::string& objName, 
//                              const ROOT::RVec<std::string>& collectionNames, 
//                              const std::string& memberName);

//       /**
//        * @brief Projects a One-To-Many association into a One-To-One column.
//        * @details Creates a flat column aligned with the track list (e.g. "rec_cal_energy").
//        * * @param outName The name of the new column.
//        * @param inputAssociation The existing One-To-Many column (e.g. "rec_cluster_energy").
//        * @param funcName The reducer function (e.g. "rad::helpers::First").
//        */
//       void DefineProjection(const std::string& outName, 
//                             const std::string& inputAssociation, 
//                             const std::string& funcName);

//       // =================================================================================
//       // Functional Methods
//       // =================================================================================

//       void SetBeamsFromMC(UInt_t ielIdx, UInt_t iionIdx, Long64_t nRows = 100);
//       void MatchCandidateToTruth(const std::string& candidateName, int targetPid);

//       // --- Helpers ---
//       void DefineForwardMatching(const std::string& prefix, int pidToMatch);
//       void DefineTruePID(const std::string& prefix);
      
//       const rad::epic::PxPyPzMVector& GetBeamElectronP4() const { return _p4el_beam; }
//       const rad::epic::PxPyPzMVector& GetBeamIonP4() const { return _p4ion_beam; }
      
//       rad::podio::PodioMetadata* GetMetadata() { return _podioMetadata.get(); }

//     private:
//       void InitMetadata(const std::string& filename);

//       enum DetID { BEAM=0, CENTRAL=1, RP=2, ZDC=3, B0=4 };
//       bool _truthMatched = false;   
//       bool _beamsCorrected = false; 
//       Int_t _idxBeamEle = -1;       
//       Int_t _idxBeamIon = -1;       
      
//       std::shared_ptr<rad::podio::PodioMetadata> _podioMetadata;
//     }; 

//     // =================================================================================
//     // IMPLEMENTATION: ePICReaction
//     // =================================================================================

//     inline ePICReaction::ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns) 
//       : ElectroIonReaction{treeName, fileNameGlob, columns} 
//     {
//         // Resolve Glob to single file for metadata
//         TChain chain("podio_metadata"); 
//         chain.Add(fileNameGlob.data());
//         if(chain.GetListOfFiles() && chain.GetListOfFiles()->GetEntries() > 0) {
//             InitMetadata(chain.GetListOfFiles()->At(0)->GetTitle());
//         }
//     }

//     inline ePICReaction::ePICReaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns) 
//       : ElectroIonReaction{treeName, filenames, columns} 
//     {
//         if(!filenames.empty()) InitMetadata(filenames[0]);
//     }

//     inline ePICReaction::ePICReaction(ROOT::RDataFrame rdf) 
//       : ElectroIonReaction{rdf} {} // Metadata remains null

//     inline void ePICReaction::InitMetadata(const std::string& filename) {
//         try {
//             _podioMetadata = std::make_shared<rad::podio::PodioMetadata>(filename);
//         } catch (const std::exception& e) {
//             std::cerr << "[ePICReaction] Warning: Failed to load metadata: " << e.what() << std::endl;
//         }
//     }

//     // --- Association Logic ---
//     inline void ePICReaction::DefineAssociation(const std::string& objName, 
//                                                 const ROOT::RVec<std::string>& collectionNames, 
//                                                 const std::string& memberName) 
//     {
//       if(!_podioMetadata) {
//           throw std::runtime_error("DefineAssociation requires PodioMetadata. Construct ePICReaction with filenames.");
//       }
//       if (collectionNames.empty()) return;

//       // 1. Resolve Collection IDs
//       ROOT::RVec<std::string> validNames;
//       ROOT::RVecU validIDs;
      
//       for(const auto& name : collectionNames) {
//           if(_podioMetadata->Exists(name)) {
//               validNames.push_back(name);
//               validIDs.push_back(_podioMetadata->CollectionIDFor(name));
//           }
//       }
      
//       if(validNames.empty()) {
//           std::cerr << "[ePICReaction] Warning: No valid collections found for " << objName << std::endl;
//           return;
//       }

//       // 2. Define ID->Index Map (Unified Index)
//       std::string unifiedIdxName = Rec() + objName + "_unified_idx" + DoNotWriteTag();
//       std::string assocCollID = "_ReconstructedParticles_" + objName + ".collectionID";
      
//       if(!ColumnExists(unifiedIdxName)) {
//            rad::podio::IdToLocalIndexMap mapper(validIDs);
//            Define(unifiedIdxName, [mapper](const RVecU& ids){ return mapper(ids); }, {assocCollID});
//       }
      
//       // 3. Prepare Inputs
//       ROOT::RVec<std::string> leafCols;
//       for(const auto& name : validNames) leafCols.push_back(name + "." + memberName);
      
//       std::string packName = Rec() + objName + "_" + memberName + "_pack" + DoNotWriteTag();
//       std::string typeName = ColObjTypeString(validNames[0] + "." + memberName); // e.g. "float" or "double"
      
//       // Define the "Vector of Vectors" packing column
//       Define(packName, Form("ROOT::RVec<ROOT::RVec<%s>>{%s}", 
//              typeName.c_str(), rad::util::combineVectorToString(leafCols).c_str()));
      
//       // 4. Define Final Output (One-To-Many)
//       std::string outputCol = Rec() + objName + "_" + memberName; 
      
//       // FIX: Ensure all string concatenations are converted to C-Strings for Form()
//       Define(outputCol, 
//              Form("rad::podio::CreateAssociation<%s>(%s, %s, %s, %s, %s)", 
//                   typeName.c_str(), 
//                   packName.c_str(), 
//                   unifiedIdxName.c_str(), 
//                   ("_ReconstructedParticles_" + objName + ".index").c_str(), 
//                   ("ReconstructedParticles." + objName + "_begin").c_str(), 
//                   ("ReconstructedParticles." + objName + "_end").c_str())
//       );
//     }

//     inline void ePICReaction::DefineProjection(const std::string& outName, 
//                                                const std::string& inputAssociation, 
//                                                const std::string& funcName) 
//     {
//       // Generates a JIT loop to apply the reduction defined in funcName
//       //i.e. converting RVec<RVec<>> to RVec<RVecResultType>
//       std::string expr = Form(
//             "rad::RVecResultType result(%s.size()); "
//             "for(const auto& v : %s) { result.push_back(static_cast<rad::ResultType_t>( %s(v) ) ); } "
//             "return result;",
//             inputAssociation.c_str(), 
//             inputAssociation.c_str(),
// 	    funcName.c_str()
//         );
//         Define(outName, expr);
//     }

//     // --- Standard Reconstructed/Truth Setup ---

//     inline void ePICReaction::SetupReconstructed(Bool_t isEnd) {
//         AddType(Rec()); 
//         rad::ParticleInjector injector(this);
        
//         ROOT::RVec<std::string> suffixes = {"double px", "double py", "double pz", "double m", "int pid", "short charge", "int det_id"};
//         if(_truthMatched) { suffixes.push_back("int match_id"); suffixes.push_back("int true_pid"); }
//         injector.DefineParticleInfo(suffixes);

//         SetBeamElectronIndex(_idxBeamEle, Rec()); 
//         SetBeamIonIndex(_idxBeamIon, Rec());
//         DefineBeamComponents(); 

//         using Source = ePICSource<ePICReaction>;
//         Source central("ReconstructedParticles", "Central_", CENTRAL);
//         central.SetIsCorrected(true); central.SetMinP(0.3); 
//         central.Process(this, injector, _truthMatched);

//         Source rp("ForwardRomanPotRecParticles", "rp_", RP);
//         rp.SetTargetPID(2212); rp.SetMinP(1);
//         rp.Process(this, injector, _truthMatched);

//         Source zdc("ReconstructedFarForwardZDCNeutrons", "ZDC_", ZDC);
//         zdc.SetTargetPID(2112); zdc.SetIsCorrected(true); zdc.SetMinP(1);
//         zdc.Process(this, injector, _truthMatched);

//         injector.CreateUnifiedVectors();
//         util::CountParticles(this, Rec());
//         if (isEnd) { DefineBeamElectron(); DefineBeamIon(); }
//     }

//     inline void ePICReaction::SetupTruth(Bool_t isEnd) {
//         AddType(Truth());
//         SetBeamElectronIndex(_idxBeamEle, Truth()); 
//         SetBeamIonIndex(_idxBeamIon, Truth());

//         rad::ParticleInjector injector(this);
//         injector.DefineParticleInfo({"double px", "double py", "double pz", "double m", "int pid", "int genStat", "int charge"});

//         injector.AddSource(Truth(), {
//           "MCParticles.momentum.x", "MCParticles.momentum.y", "MCParticles.momentum.z",
//           "MCParticles.mass", "MCParticles.PDG", "MCParticles.generatorStatus", "MCParticles.charge"
//         }, "MCParticles.generatorStatus==1||MCParticles.generatorStatus==4"); 

//         injector.CreateUnifiedVectors();
//         util::CountParticles(this, Truth());
//         if (isEnd) { DefineBeamElectron(); DefineBeamIon(); }
//     }

//     inline void ePICReaction::SetupMatching(Bool_t isEnd) {
//         _truthMatched = true;
//         Define("Central_match_id" + DoNotWriteTag(), 
//             [](const ROOT::RVecU& recID, const ROOT::RVecU& simID, const Indices_t& rec_ind) {
//                 Indices_t match_id(rec_ind.size(), -1);
//                 for(size_t i = 0; i < recID.size(); ++i) {
//                     if(recID[i] < (uint)match_id.size()) match_id[recID[i]] = (int)simID[i];
//                 }
//                 return match_id;
//             }, 
//             {"ReconstructedParticleAssociations.recID", "ReconstructedParticleAssociations.simID", "ReconstructedParticles.PDG"}
//         );
//         SetupReconstructed(kFALSE);
//         SetupTruth(kFALSE);
//     }

//     inline void ePICReaction::SetBeamsFromMC(UInt_t iel, UInt_t iion, Long64_t nRows) {
//         _useBeamsFromMC = true;
//         _idxBeamEle = iel; _idxBeamIon = iion;
//         auto nthreads = ROOT::GetThreadPoolSize();
//         if (nthreads) ROOT::DisableImplicitMT();

//         auto tempframe = GetFileNames().empty() ? ROOT::RDataFrame{GetTreeName(), GetFileName()} : ROOT::RDataFrame{GetTreeName(), GetFileNames()};
//         auto beamdf = tempframe.Range(nRows)
//             .Define("emean", Form("MCParticles.momentum.z[%d]", iel))
//             .Define("pzmean", Form("MCParticles.momentum.z[%d]", iion))
//             .Define("pxmean", Form("MCParticles.momentum.x[%d]", iion));
          
//         auto pze = beamdf.Mean("emean");
//         auto pzp = beamdf.Mean("pzmean");
//         auto pxp = beamdf.Mean("pxmean");

//         _p4el_beam.SetPxPyPzE(0, 0, *pze, std::abs(*pze));
//         _p4ion_beam.SetPxPyPzE(*pxp, 0, *pzp, std::sqrt(*pxp * *pxp + *pzp * *pzp + 0.938*0.938));
    
//         std::cout << " [ePICReaction] Beams Initialized: Ele=" << _p4el_beam.Pz() << " GeV, Ion=" << _p4ion_beam.Pz() << " GeV" << std::endl;
//         if (nthreads) ROOT::EnableImplicitMT(nthreads);
//     }

//     inline void ePICReaction::MatchCandidateToTruth(const std::string& candidateName, int targetPid) {
//         std::string idxCol = candidateName + "_idx";
//         Define(candidateName + "_is_true", 
//                [targetPid](const Indices_t& candIdx, const Indices_t& recMatchId, const Indices_t& rawPdg) {
//                    if(candIdx.empty() || candIdx[0] == -1) return 0;
//                    int rIdx = candIdx[0];
//                    if(rIdx >= (int)recMatchId.size()) return 0;
//                    int tIdx = recMatchId[rIdx];
//                    if(tIdx == -1 || tIdx >= (int)rawPdg.size()) return 0;
//                    return (std::abs(rawPdg[tIdx]) == targetPid) ? 1 : 0;
//                }, 
//                {idxCol, Rec() + "match_id", "MCParticles.PDG"});
//     }

//     inline void ePICReaction::DefineForwardMatching(const std::string& prefix, int pidToMatch) {
//         Define(prefix + "match_id" + DoNotWriteTag(), 
//             [pidToMatch](const Indices_t& pdg, const Indices_t& stat, const Indices_t& id_vec) {
//                 int best_idx = -1;
//                 for(size_t i=0; i<pdg.size(); ++i) {
//                     if(std::abs(pdg[i]) == pidToMatch && stat[i] == 1) {
//                         best_idx = i; break; 
//                     }
//                 }
//                 auto result = Indices_t(id_vec.size(), best_idx);
//                 return result;
//             }, {"MCParticles.PDG", "MCParticles.generatorStatus", prefix + "det_id" + DoNotWriteTag()});
//     }

//     inline void ePICReaction::DefineTruePID(const std::string& prefix) {
//         Define(prefix + "true_pid" + DoNotWriteTag(), 
//                [](const Indices_t& match_id, const Indices_t& tru_pdg) {
//                    Indices_t tpid(match_id.size(), 0);
//                    for(size_t i=0; i<(int)match_id.size(); ++i) {
//                        if(match_id[i] != -1 && match_id[i] < (int)tru_pdg.size()) tpid[i] = tru_pdg[match_id[i]];
//                    }
//                    return tpid;
//                }, {prefix + "match_id" + DoNotWriteTag(), "MCParticles.PDG"});
//     }

//   } // namespace epic
// } // namespace rad
