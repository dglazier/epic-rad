#include "CommonDefines.h"
#include "ePICReaction.h"
#include "KinematicsProcElectro.h"
#include "KineCalculation.h"
#include "Indicing.h"
#include "ElectronScatterKinematics.h"
#include "BasicKinematicsRDF.h"
#include <TBenchmark.h>

/**
 * @brief Analysis Example: Y(4260) -> J/psi pi+ pi-
 * * Strategy:
 * 1. Define Master Processor for the "Detected Recoil" topology (most inclusive).
 * 2. Clone it for the "Missing Recoil" hypothesis.
 * 3. Clone both for Truth.
 */
void ProcessRecTruthYCombi() {
  
  using namespace rad;
  using namespace rad::consts::data_type; 

  gBenchmark->Start("df");

  // =================================================================================
  // 1. SETUP & MATCHING
  // =================================================================================
  epic::ePICReaction rad_df{"events", "/home/dglazier/EIC/data/Y4260/jpac_y4260_18_275_10day_*_recon.root"};
  rad_df.SetBeamsFromMC(0, 1); 
  rad_df.SetupMatching(); 

  const int Role_ScatEle  = 2; 
  const int Role_Recoil   = 3; // Proton
  const int Role_PiM      = 4; 
  const int Role_PiP      = 5; 
  const int Role_DecayEle = 6; 
  const int Role_DecayPos = 7; 
  
  rad_df.SetParticleCandidates(consts::ScatEle(), Role_ScatEle, rad::index::FilterIndices(11), {"rec_true_pid"});
  rad_df.SetParticleCandidates("ele", Role_DecayEle, rad::index::FilterIndices(11),  {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pos", Role_DecayPos, rad::index::FilterIndices(-11), {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pip", Role_PiP, rad::index::FilterIndices(211),  {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pim", Role_PiM, rad::index::FilterIndices(-211), {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("recoil", Role_Recoil, rad::index::FilterIndices(2212), {"rec_true_pid"}); 

  
  rad_df.MakeCombinations();
  rad_df.DefineTrueMatchedCombi("rec_match_id", Rec());

  // =================================================================================
  // 2. KINEMATICS CONFIGURATION (Master: DETECTED Recoil)
  // =================================================================================
  // We start with the topology that explicitly uses the "recoil" proton.
  // This ensures "recoil" is mapped naturally without needing RequireParticle.
  
  KinematicsProcElectro kine_det_rec(&rad_df, Rec());
 
  // A. Topological Construction
  kine_det_rec.Creator().Sum("Jpsi", {{"ele", "pos"}});       
  kine_det_rec.Creator().Sum("Y",    {{"Jpsi", "pip", "pim"}});
  
  //Also define the missing mass variable here so it's available for plotting/cuts
  kine_det_rec.Creator().Diff("recoil_miss", 
         {{consts::BeamIon(), consts::BeamEle()},  
          {"Y", consts::ScatEle()}}             
  );

  // B. Physics Grouping (Master = Detected)
  kine_det_rec.SetMesonParticles({"Y"});
  kine_det_rec.SetBaryonParticles({"recoil"}); // Uses detected proton

  // C. Calculations
  kine_det_rec.Mass("MassMeson", {"Y"});          
  kine_det_rec.Mass("MassMiss", {"recoil_miss"});
  kine_det_rec.Q2();
  
  // Register 't' calculation (using detected 'recoil')
  kine_det_rec.RegisterCalc("tb", rad::physics::TBot); 

  // =================================================================================
  // 3. CLONING PHYSICS STREAMS
  // =================================================================================

  // Stream 2: Detected Recoil (Truth)
  // Inherits "recoil" grouping and "t_det" calculation.
  auto kine_det_tru = kine_det_rec.CloneForType(Truth());

  // Stream 3: Missing Recoil (Rec) - [Linked Clone]
  // We CloneLinked to create a new hypothesis on the same data.
  // Suffix "_miss" ensures we get unique columns like "rec_t_miss".
   auto kine_miss_rec = kine_det_rec.CloneLinked("_miss");
  
  // Override: Use the missing mass particle instead of detected recoil
  kine_miss_rec->SetBaryonParticles({"recoil_miss"}); 
  
  // Stream 4: Missing Recoil (Truth)
  // Clone from the Missing Rec stream to inherit the correct grouping.
  auto kine_miss_tru = kine_miss_rec->CloneForType(Truth());

  // =================================================================================
  // 4. EXECUTION
  // =================================================================================
  
  // Initializing the Master (kine_det_rec) triggers the entire chain.
  kine_det_rec.Init(); 

  gBenchmark->Start("snapshot");
  
  // File will contain:
  // - rec_t_det, tru_t_det (from Master & Stream 2)
  // - rec_t, tru_t (from Stream 3 & 4)
  rad_df.Snapshot("ePIC_RecTruth_Y4260_Full.root");
  
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");

  // Debug Maps
  std::cout << "\n--- Master Map (Rec Detected) ---" << std::endl;
  kine_det_rec.PrintReactionMap();
  
 }
