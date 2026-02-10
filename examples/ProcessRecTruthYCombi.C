#include "AnalysisManager.h"
#include "ePICReaction.h"
#include "KinematicsProcElectro.h"
#include "ElectronScatterKinematics.h"
#include <TBenchmark.h>

/**
 * @brief Analysis Example: Y(4260) -> J/psi pi+ pi-
 * * Strategy:
 * 1. Define base processor for given ReconstructedParticles.
 * 2. Define additoinal processor for using cluster energy for e+,e-
 * 3. Produce and analyse all valid combinations
 */
void ProcessRecTruthYCombi() {
  //ROOT::EnableImplicitMT();

  using namespace rad;
  using namespace rad::consts::data_type; 

  gBenchmark->Start("df");

  //Define template ConfigReaction and KinematicsProcessor for ePIC analysis
  using Reaction = epic::ePICReaction;
  using Processor = KinematicsProcElectro;
  
  // =================================================================================
  // 1. SETUP & MATCHING
  // =================================================================================
  AnalysisManager<Reaction,Processor>  mgr{
    "Y4260",
    "events",
    "/home/dglazier/EIC/data/Y4260/jpac_y4260_18_275_10day_*_recon.root"};

  mgr.SetOutputDir("output");
  auto& rad_df = mgr.Reaction();

  //Indices of electron,ion in MCParticles
  rad_df.SetBeamsFromMC(0, 1); 
  rad_df.SetupMatching(); 

  //Truth Ordering in MCParticles vector
  const int Role_ScatEle  = 2; 
  const int Role_Recoil   = 3; // Proton
  const int Role_PiM      = 4; 
  const int Role_PiP      = 5; 
  const int Role_DecayEle = 6; 
  const int Role_DecayPos = 7; 

  //arguments (name, truth_idx, candidate list function, arguments for list function)
  rad_df.SetParticleCandidates(consts::ScatEle(), Role_ScatEle, rad::index::FilterIndices(11), {"rec_true_pid"});
  rad_df.SetParticleCandidates("ele", Role_DecayEle, rad::index::FilterIndices(11),  {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pos", Role_DecayPos, rad::index::FilterIndices(-11), {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pip", Role_PiP, rad::index::FilterIndices(211),  {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pim", Role_PiM, rad::index::FilterIndices(-211), {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("p", Role_Recoil, rad::index::FilterIndices(2212), {"rec_true_pid"}); 

  //create indices for each combination, based in given candidates
  rad_df.MakeCombinations();

  //Add different processing streams for these combi events
  mgr.AddStream(Rec(),"base"); //base rec stream
  mgr.AddStream(Truth(),"base"); //base tru stream
  mgr.AddStream(Rec(),"calE"); //rec stream using calorimeter energies for e+e-

  // --- Auxiliary Data Handling ---
  //collect and associate cluster energies with particles
  //rec_cal_energy will be synched with momentum ordering
  rad_df.DefineAssociation("clusters", {"EcalBarrelClusters", "EcalEndcapPClusters"}, "energy");
  rad_df.DefineProjection("rec_cal_energy", "rec_clusters_energy", "rad::util::First");


  // =================================================================================
  // 2. ANALYSIS CONFIGURATION 
  // =================================================================================
   
  // [A] SHARED KINEMATICS (Topology)
  // Applied to both Rec and Truth streams. Defines the decay chain.
  auto topology_recipe = [](Processor& p) {
    using namespace consts;
    // 1. Reconstruct J/psi -> e+ e-
    // A. Topological Construction
    p.Creator().Sum("Jpsi", {{"ele", "pos"}});       
    p.Creator().Sum("TwoPi", {{"pip", "pim"}});       
    p.Creator().Sum("Y",    {{"Jpsi", "TwoPi"}});
    //p.Creator().Diff("Miss",    {{BeamEle(),BeamIon()},{"Jpsi", "TwoPi"}});
    p.Creator().Diff("Miss",    {{BeamEle(),BeamIon()},{"Jpsi", "TwoPi","p"}});
  
    p.SetMesonParticles({"Jpsi","TwoPi"});
    p.SetBaryonParticles({"p"});
     //p.SetBaryonParticles({"Miss"});
    
    // 3. Calculate Invariant Masses
    p.Mass("MassJ",     {"ele","pos"});             
    //    p.Mass("MassJ",     {"Jpsi"});//fix mass in postmodifer             
    p.Mass("MassTwoPi",     {"TwoPi"});             
    p.Mass("MassY",     {"Y"});             
    p.Mass2("MissMass2",     {"Miss"});             

    p.Q2();
    
    //4. Calculate Mandelstam t (requires beam definition)
    p.RegisterCalc("tb", rad::physics::TBot);
    p.RegisterCalc("DeltaPhiYxP", rad::DeltaPhi,{{"Y","p"}});
    
    p.ParticleTheta({"scat_ele","Y"});
    p.ParticlePhi({"scat_ele","Y"});
    p.ParticleP({"scat_ele","Y"});
    
  };

  // Apply Topology to ALL streams
  mgr.ConfigureKinematics(topology_recipe);

  
  // [B] (i) Apply some pre and post momentum organising modifiers
  // Fixes all masses to PDG values
  auto mass_recipe = [](KinematicsProcessor& p) {
    using namespace consts;
    //use PDG mas values for all particles in case
    //they were PIDed wrong
    p.PreModifier().FixMass(ScatEle(), M_ele());
    p.PreModifier().FixMass("ele", M_ele());
    p.PreModifier().FixMass("pos", M_ele());
    p.PreModifier().FixMass("pip", M_pi());
    p.PreModifier().FixMass("pim", M_pi());
    p.PreModifier().FixMass("p", M_pro());
  
    //Fix reconstructed Jpsi mass after it is calculted
    p.PostModifier().FixMass("Jpsi", M_Jpsi());

    
  };
  // Apply Mass Corrections to REC stream ONLY
  mgr.ConfigureKinematics(Rec(), mass_recipe);

  // (ii) Corrects electron momentum using ECal energy.
  auto correction_recipe = [](KinematicsProcessor& p) {
    p.PreModifier().SetMomentumFrom("ele", "rec_cal_energy");
    p.PreModifier().SetMomentumFrom("pos", "rec_cal_energy");
  };
  // Apply Corrections to REC calE stream ONLY
  mgr.ConfigureKinematics(Rec()+"calE", correction_recipe);

    
  // [C] SELECTION CUTS
  // Applied to Rec stream (Primary). Determines which combinations are saved.
  // Can use variables defined in kinematics recipe
  auto selection_recipe = [](PhysicsSelection& s) {
    // Wide window for Y(4260)
    s.AddCutRange("Y_mass_cut",    "MassY", 3.5, 4.5); 
    // Loose cut for J/psi
    s.AddCutMin("Jpsi_mass_cut", "MassJ",     2.8);          
  };

  // Apply Cuts to Rec_base only
  //mgr.ConfigureSelection(Rec()+"base", selection_recipe);
  // Apply Cuts to ALL
  mgr.ConfigureSelection(selection_recipe);

  // [D] HISTOGRAMS
  // Shared definitions for Rec and Truth.
  auto histogram_recipe = [](histo::Histogrammer& h) {
    h.Create("hQ2",     "Q^{2}; [GeV]", 100, 0, 1.0, "Q2");
    h.Create("hMassJ",     "Invariant Mass e+e-;  [GeV]", 100, 2.0, 4.0, "MassJ");
    h.Create("hMassTwoPi",     "Invariant Mass 2#pi; [GeV]", 100, 0, 2, "MassTwoPi");
    h.Create("hMassY",     "Invariant Mass J/#psi 2#pi; [GeV]", 100, 2.0, 6.0, "MassY");
    h.Create("hMissMass2",     "Missing Mass squared; [GeV]", 1000, -50, 50, "MissMass2");
  };
 
  // Apply Histograms to ALL streams
  mgr.ConfigureHistograms(histogram_recipe);
  
  // [D] TREES
  mgr.Snapshot({consts::TruthMatchedCombi()});//currently need to add isTruth branch

  
  // Print diagnostics BEFORE running expensive event loop
  std::cout << "\n=== CHECKING ANALYSIS SETUP ===\n" << std::endl;
  mgr.PrintDiagnostics();
  //  PrintDefinedColumnNames(mgr.Reaction().CurrFrame());
  
  // =================================================================================
  // 3. RUN IT ALL 
  // =================================================================================
  gBenchmark->Start("analysis");
  mgr.Run();
  gBenchmark->Stop("analysis");
  gBenchmark->Print("analysis");
 }
