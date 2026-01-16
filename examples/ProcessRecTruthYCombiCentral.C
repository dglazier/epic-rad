#include "CommonDefines.h"
#include "ePICReaction.h"
#include "KinematicsProcElectro.h"
#include "KineCalculation.h"
#include "Indicing.h"
#include "ElectronScatterKinematics.h"
#include "BasicKinematicsRDF.h"
#include "PhysicsSelection.h"
#include "Histogrammer.h" 
#include <TBenchmark.h>

/**
 * @brief Analysis Example: Y(4260) -> J/psi pi+ pi-
 * @details
 * This macro demonstrates the complete RDataFrame analysis chain:
 * 1. Data Loading & Particle Matching
 * 2. Kinematic Reconstruction (Combinatorics) for Rec AND Truth
 * 3. Physics Selection (Lazy Masking)
 * 4. Snapshotting (Single Tree with Rec + Scalarized Truth)
 */
void ProcessRecTruthYCombiCentral() {
    
    // Enable multi-threading for parallel processing
    ROOT::EnableImplicitMT(2);

    using namespace rad;
    using namespace rad::consts::data_type; 

    gBenchmark->Start("df");

    // =================================================================================
    // 1. SETUP & MATCHING
    // =================================================================================
    // Initialize the Reaction framework 
    epic::ePICReaction rad_df{"events", "/home/dglazier/EIC/data/Y4260/jpac_y4260_18_275_10day_*_recon.root"};
    
    // Configure Beam Energies (Electron=0, Proton=1 from MC truth)
    rad_df.SetBeamsFromMC(0, 1); 
    rad_df.SetupMatching(); 

    // Define Particle Roles (MC indices)
    const int Role_ScatEle  = 2; 
    const int Role_Recoil   = 3; // Proton
    const int Role_PiM      = 4; 
    const int Role_PiP      = 5; 
    const int Role_DecayEle = 6; 
    const int Role_DecayPos = 7; 
    
    // Define Particle Candidates (Rec ID -> MC Role mapping)
    rad_df.SetParticleCandidates("ele", Role_DecayEle, rad::index::FilterIndices(11),  {"rec_true_pid"}); 
    rad_df.SetParticleCandidates("pos", Role_DecayPos, rad::index::FilterIndices(-11), {"rec_true_pid"}); 
    rad_df.SetParticleCandidates("pip", Role_PiP, rad::index::FilterIndices(211),  {"rec_true_pid"}); 
    rad_df.SetParticleCandidates("pim", Role_PiM, rad::index::FilterIndices(-211), {"rec_true_pid"}); 
    
    // Generate all valid combinations of candidates
    rad_df.MakeCombinations();
    
    // Flag the "True" combination using MC truth matching
    // This creates "rec_is_signal_combi" based on the Roles defined above
    rad_df.DefineTrueMatchedCombi("rec_match_id", Rec());

    // --- Auxiliary Data Handling ---
    rad_df.DefineAssociation("clusters", {"EcalBarrelClusters", "EcalEndcapPClusters"}, "energy");
    rad_df.DefineProjection("rec_cal_energy", "rec_clusters_energy", "rad::util::First");


    // =================================================================================
    // 2. KINEMATICS CONFIGURATION (REC)
    // =================================================================================
    // Create the Processor for Reconstructed Data
    KinematicsProcessor kine_det_rec(&rad_df, Rec());
    
    // Modifiers
    kine_det_rec.PreModifier().SetMomentumFrom("ele", "rec_cal_energy");
    kine_det_rec.PreModifier().SetMomentumFrom("pos", "rec_cal_energy");

    // Topological Construction
    kine_det_rec.Creator().Sum("Jpsi", {{"ele", "pos"}});         
    kine_det_rec.Creator().Sum("Y",    {{"Jpsi", "pip", "pim"}});
    
    // Physics Grouping
    kine_det_rec.SetMesonParticles({"Y"});
    kine_det_rec.SetBaryonParticles({}); 

    // Calculations
    kine_det_rec.Mass("MassMeson", {"Y"});            
    kine_det_rec.Mass("MassJ", {"Jpsi"});            
    kine_det_rec.RegisterCalc("tb", rad::physics::TBot); 


    // =================================================================================
    // 3. KINEMATICS CONFIGURATION (TRUTH)
    // =================================================================================
    // Clone configuration for Truth (MC) particles. 
    // This automatically maps "rec_ele" -> "tru_ele", etc.
    auto kine_det_tru = kine_det_rec.CloneForType(Truth());

    // Initialize BOTH chains (Note: DefineNewComponentVecs is called internally by Init)
    kine_det_rec.Init(); 
    kine_det_tru->Init();


    // =================================================================================
    // 4. SELECTION (The "Funnel")
    // =================================================================================
    // Bind Selection to the Processor.
    rad::PhysicsSelection sel(kine_det_rec);
    
    // Define Cuts (Boolean Masks)
    sel.AddCutRange("Y_mass_cut", "MassMeson", 3.5, 4.5); 
    sel.AddCutMin("Jpsi_mass_cut", "MassJ", 2.8);         
    
    sel.Compile(); 

    // =================================================================================
    // 5. HISTOGRAMMING
    // =================================================================================
    // Bind Histogrammer to Processor (names) and Selection (mask)
    rad::histo::Histogrammer hist(kine_det_rec, &sel);

    // Optional: Add a Split (e.g., bin in 't'). All histograms become 2D (Var vs t).
    // hist.AddSplit("t-channel", "tb", 10, 0, 2.0); 

    // Define Histograms: (Name, Title, Bins, Min, Max, VariableName)
    // Only combinations passing 'sel' are filled.
    hist.Create("MassMeson", "Invariant Mass Y; Mass [GeV]", 100, 0, 5, "MassMeson");
    hist.Create("MassJ",     "Invariant Mass J/psi; Mass [GeV]", 100, 0, 4, "MassJ");

    // =================================================================================
    // 5. SCALARIZE TRUTH VARIABLES
    // =================================================================================
    // We want to write the Truth Mass to the Rec Tree. 
    // Since Rec is combinatorial (N entries) and Truth is per-event (1 entry), 
    // we extract the 0-th element of the Truth vector to create a SCALAR.
    // SnapshotCombi will automatically broadcast this scalar to all Rec combinations.
    
    rad_df.Define("tru_MassMeson_scalar", "tru_MassMeson.at(0)");
    rad_df.Define("tru_MassJ_scalar",     "tru_MassJ.at(0)");


    // =================================================================================
    // 6. SNAPSHOT (SINGLE FLAT TREE)
    // =================================================================================
    
    std::vector<std::string> tree_vars = {
        // --- Reconstructed Kinematics (Vectors -> Flattened) ---
        "rec_MassMeson", "rec_MassJ",
        "rec_ele_px", "rec_ele_py", "rec_ele_pz", "rec_ele_m",
        "rec_pos_px", "rec_pos_py", "rec_pos_pz", "rec_pos_m",
        "rec_pip_px", "rec_pip_py", "rec_pip_pz", "rec_pip_m",
        "rec_pim_px", "rec_pim_py", "rec_pim_pz", "rec_pim_m",
         "is_signal_combi",
	 "tru_MassMeson_scalar", 
        "tru_MassJ_scalar"
    };

    // Lazy Booking
    rad_df.BookSnapshotCombi("ePIC_RecTruth_Y4260_Tree.root", "rad_tree", tree_vars, sel.GetMaskColumn());


    // =================================================================================
    // 7. EXECUTION
    // =================================================================================
    hist.File("ePIC_RecTruth_Y4260_Hist.root");
    // Explicitly trigger the snapshots to ensure data is written
    rad_df.TriggerSnapshots();

    gBenchmark->Stop("df");
    gBenchmark->Print("df");

    // Debug Maps
    std::cout << "\n--- Master Map (Rec Detected) ---" << std::endl;
    kine_det_rec.PrintReactionMap();
}
