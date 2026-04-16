#include "AnalysisManager.h"
#include "ePICReaction.h"
#include "ePICAssociationsManager.h" // <-- NEW: Include the manager
#include "KinematicsProcessor.h"
#include "BasicKinematicsRDF.h"
#include <TBenchmark.h>
#include <vector>
#include <memory>
#include <string>

const std::string my_out_dir = "TCS_18x275_hplus_test/";

void ProcessRecTruthTCSSingleTrack(
    std::vector<std::string> infiles={"/home/dglazier/EIC/EpIC1.1.6-1.0_DDVCS_18x275_q2_0_10_edecay_hplus_abconv_run0.0000.eicrecon.edm4eic.root"}, 
    std::string outdir = my_out_dir,
    const int Role_ScatEle  = 2,
    const int Role_Recoil   = 5,
    const int Role_DecayEle = 6, 
    const int Role_DecayPos = 7) 
{
  using namespace rad;
  using namespace rad::consts::data_type; 
  using Reaction = epic::ePICReaction;
  using Processor = KinematicsProcessor; // Use base processor to avoid Q2/W overhead

  gBenchmark->Start("df");

  if(infiles[0] == ""){
    std::cout << "Infiles empty, grabbing test file!" << std::endl;
    // Assuming GetXRootDFiles is defined elsewhere in your macros
    infiles = {"root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/26.02.0/epic_craterlake/EXCLUSIVE/DDVCS_ABCONV/EpIC1.1.6-1.0/18x275/q2_0_10/edecay/hplus/edm4eic.root"};
  }

  // =================================================================================
  // 1. BASE REACTION SETUP (Shared by all managers)
  // =================================================================================
  Reaction base_df{"events", infiles};
  base_df.SetBeamsFromMC(0, 3); 
  base_df.SetupMatching();
  
  // --- NEW: Auxiliary Data Handling via Fluent Builder ---
  epic::ePICAssociationManager assoc(base_df);

  // Extract Central Calorimeter Energies
  assoc.For("Central")
       .From({"EcalBarrelClusters", "EcalEndcapPClusters"}) // Optionally add "EcalEndcapNClusters" if needed!
       .Extract("energy")
       .As("cal_energy"); // Maps to "rec_cal_energy"

  // Extract Tagger Tracker PIDs
  assoc.For("Central")
       .Relation("tracks") // Override default "clusters" relation
       .From({"TaggerTrackerReconstructedParticles"})
       .Extract("PDG")
       .As("tracks_pid", "int"); // Maps to "rec_tracks_pid", explicitly set to integer

  // Execute padding and build the unified arrays for the base dataframe
  assoc.Build();
  // -------------------------------------------------------
  
  // The detector flag logic remains securely in ePICReaction
  // requires podio association linked with detector you are flagging
  base_df.DefineDetectorFlag("rec_from_tagger", "TaggerTrackerReconstructedParticleAssociations");

  // =================================================================================
  // 2. RECIPE GENERATORS
  // =================================================================================
  
  // Creates a Kinematics recipe for a specific particle
  auto make_kine_recipe = [](const std::string& pName) {
      return [pName](Processor& p) {
          p.ParticleTheta({pName});
          p.ParticlePhi({pName});
          p.ParticleP({pName});    
          p.ParticleEta({pName});
      };
  };

  // Creates a Correction recipe for a specific particle
  auto make_corr_recipe = [](const std::string& pName, double mass) {
      return [pName, mass](Processor& p) {
          p.PreModifier().FixMass(pName, mass);
      };
  };

  // Creates a Histogram recipe for a specific particle
  auto make_hist_recipe = [](const std::string& pName) {
      return [pName](histo::Histogrammer& h) {
          h.Create(pName + "_pmag", pName + " Momentum; p [GeV/c]", 100, 0, 275, pName + "_pmag");
          h.Create(pName + "_eta", pName + " Pseudorapidity; #eta", 100, -10, 10, pName + "_eta");
          h.Create(pName + "_theta", pName + " Polar Angle; #theta [rad]", 100, 0.0, TMath::Pi(), pName + "_theta");
          h.Create(pName + "_phi", pName + " Azimuthal Angle; #phi [rad]", 100, -TMath::Pi(), TMath::Pi(), pName + "_phi");
      };
  };

  // Truth match requirement
  auto match_recipe = [](PhysicsSelection& s) {
      s.AddCutBool("match_cut", rad::consts::TruthMatchedCombi()); 
  };

  // =================================================================================
  // 3. SPAWN PARTICLE MANAGERS (Using a Builder Lambda)
  // =================================================================================

  // Vector to securely hold our managers in memory
  std::vector<std::unique_ptr<AnalysisManager<Reaction, Processor>>> managers;

  // The Builder Lambda: 'auto filter_func' allows both FilterIndices and FilterIndicesWithFlag
  auto build_stream = [&](const std::string& pName, int role, auto filter_func, 
                          const std::vector<std::string>& filterCols, double mass) 
  {
      // Clones the base_df so it automatically inherits the built auxiliary columns!
      auto mgr = std::make_unique<AnalysisManager<Reaction, Processor>>("TCS_" + pName, base_df);
      mgr->SetOutputDir(outdir);
      
      // 1. Setup Candidates
      mgr->Reaction().SetParticleCandidates(pName, role, filter_func, filterCols); 

      // 2. The Hack: Register EVERYTHING as a "Meson" just to get it into the ReactionMap
      mgr->Reaction().SetMesonParticles({pName});

      // 3. Make Combinations
      mgr->Reaction().MakeCombinations();
      
      // 4. Add Streams
      mgr->AddStream(Truth(), "");
      mgr->AddStream(Rec(), "");
      
      // 5. Apply Recipes
      mgr->ConfigureKinematics(make_kine_recipe(pName));
      mgr->ConfigureKinematics(Rec(), make_corr_recipe(pName, mass));
      ///////mgr->ConfigureSelection(Rec(), match_recipe);
      mgr->ConfigureHistograms(make_hist_recipe(pName));
      
      // 6. Store safely
      managers.push_back(std::move(mgr));
  };

  // --- Build all particle streams ---
  // Now they all use the exact same signature!
  build_stream(consts::ScatEle(), Role_ScatEle, rad::index::FilterIndicesWithFlag(11), {"rec_true_pid", "rec_from_tagger"}, consts::M_ele());
  build_stream("ele", Role_DecayEle, rad::index::FilterIndices(11), {"rec_true_pid"}, consts::M_ele());
  build_stream("pos", Role_DecayPos, rad::index::FilterIndices(-11), {"rec_true_pid"}, consts::M_ele());
  build_stream("pprime", Role_Recoil, rad::index::FilterIndices(2212), {"rec_true_pid"}, consts::M_pro());
  
  // =================================================================================
  // 4. BOOK TREES & RUN
  // =================================================================================
  
  // Book all snapshots (Lazy Execution)
  for (auto& mgr : managers) {
      mgr->Snapshot();
  }

  gBenchmark->Start("analysis");

  // The first Run() triggers the single event loop across the entire underlying RDataFrame.
  // The subsequent loops just finalize the TFile writes.
  for (auto& mgr : managers) {
      mgr->Run();
  }

  gBenchmark->Stop("analysis");
  gBenchmark->Print("analysis");
}
