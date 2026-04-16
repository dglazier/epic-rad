# ePIC-RAD: Reaction Analysis & Design for ePIC

**A Dedicated ePIC/EDM4hep Extension for the RAD Framework**

ePIC-RAD is an extension of the base RAD framework specifically tailored for processing electron scattering reactions using PODIO/EDM4hep style data. It abstracts the complexities of detector hits, truth matching, and combinatorics into a high-level declarative API.

### Project Goals
1. Simplifying analysis code for general final states.
2. Maintaining the same code structure for different types of data (e.g., HepMC, ePIC).
3. Requiring no additional dependencies beyond ROOT, with no build system required (header-only, parsed by ROOT at runtime).
4. Defining final state particles and standardizing branch outputs with minimal lines of code.
5. Hiding boilerplate and low-level C++ mechanics from the user.
6. Automating MC matching, the calculation of equivalent truth variables, and full combinatorial analysis.

---

## ⚡ Installation

To install, simply download the code from Git and add the paths to your `ROOT_INCLUDE_PATH`. 

If you do not have the base `rad` code already installed, you can clone it simultaneously as a submodule:
```bash
git clone --recurse-submodules [https://github.com/dglazier/epic-rad.git](https://github.com/dglazier/epic-rad.git)
export EPICRAD=/path/to/epic-rad
export RAD=${EPICRAD}/rad
```

If you already have `rad` installed elsewhere on your system, you can clone just the extension:
```bash
git clone [https://github.com/dglazier/epic-rad.git](https://github.com/dglazier/epic-rad.git)
export EPICRAD=/path/to/epic-rad
export RAD=/path/to/existing/rad
```

Finally, make the headers visible to your interactive ROOT session by appending them to your include path:
```bash
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${RAD}/include:${EPICRAD}/include
```

---

## 🚀 Quick Start Example

Here is a complete analysis steering macro reconstructing $Y(4260) \to J/\psi(\to e^+e^-) \pi^+ \pi^-$. This demonstrates the core workflow: setting up the data, extracting detector hits, generating combinations, and defining kinematics.

```cpp
#include "AnalysisManager.h"
#include "ePICReaction.h"
#include "ePICAssociationsManager.h" 
#include "KinematicsProcElectro.h"

void Analysis() {
  // 1. Initialize Reaction & Data Source
  // ePICReaction handles Podio file reading and Rec-to-Truth matching automatically
  rad::AnalysisManager<epic::ePICReaction, rad::KinematicsProcElectro> mgr{
    "Y4260", "events", "input_data/*_recon.root"
  };
  
  auto& rad_df = mgr.Reaction();
  rad_df.SetBeamsFromMC(0, 1); 
  rad_df.SetupMatching(); 

  // 2. Auxiliary Data (Detector Associations)
  // Extract calorimeter energies and dynamically pad them to sync with combinatorial tracks
  epic::ePICAssociationManager assoc(rad_df);
  assoc.For("Central")
       .From({"EcalBarrelClusters", "EcalEndcapPClusters", "EcalEndcapNClusters"})
       .Extract("energy")
       .As("cal_energy"); 
  assoc.Build(); 

  // 3. Define Candidates (Name, MC_Role_ID, Filter, Columns)
  rad_df.SetParticleCandidates(rad::consts::ScatEle(), 2, rad::index::FilterIndices(11), {"rec_true_pid"});
  rad_df.SetParticleCandidates("ele", 6, rad::index::FilterIndices(11),  {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pos", 7, rad::index::FilterIndices(-11), {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pip", 5, rad::index::FilterIndices(211), {"rec_true_pid"}); 
  rad_df.SetParticleCandidates("pim", 4, rad::index::FilterIndices(-211),{"rec_true_pid"}); 
  rad_df.SetParticleCandidates("p",   3, rad::index::FilterIndices(2212),{"rec_true_pid"}); 

  // 4. Generate All Unique Combinations
  rad_df.MakeCombinations();

  // 5. Setup Processing Streams
  mgr.AddStream(rad::consts::data_type::Rec(), "base");
  mgr.AddStream(rad::consts::data_type::Truth(), "base");

  // 6. Configure Kinematics (Topology & Variables)
  auto topology_recipe = [](rad::KinematicsProcElectro& p) {
    using namespace rad::consts;
    
    // Topological Construction
    p.Creator().Sum("Jpsi", {{"ele", "pos"}});       
    p.Creator().Sum("TwoPi", {{"pip", "pim"}});       
    p.Creator().Sum("Y",    {{"Jpsi", "TwoPi"}});
    p.Creator().Diff("Miss", {{BeamEle(), BeamIon()}, {"Jpsi", "TwoPi", "p", ScatEle()}});
    
    // Invariant Masses & Physics
    p.Mass("MassJ", {"ele","pos"});             
    p.Mass("MassY", {"Y"});             
    p.Q2();
    p.RegisterCalc("tb", rad::physics::TBot, {{BeamIon(), "p"}});

    //  Calculate Mandelstam t and other injected calculations
    p.RegisterCalc("tb", rad::physics::TBot);
    p.RegisterCalc("DeltaPhiYxP", rad::DeltaPhi, {{"Y","p"}});
    
    p.ParticleTheta({"scat_ele","Y","ele","pos"});
    p.ParticlePhi({"scat_ele","Y","ele","pos"});
    p.ParticleP({"scat_ele","Y","ele","pos"});

    // Auxiliary Passthrough
    // Pulls the perfectly unified auxiliary column directly into the Flat output tree
    p.PassThrough("ele", "rec_cal_energy", "_cal_energy");
    p.PassThrough("pos", "rec_cal_energy", "_cal_energy"); 
 
  };
  mgr.ConfigureKinematics(topology_recipe);

  // 7. Save Flat Tree
  mgr.Snapshot({rad::consts::TruthMatchedCombi()});
  
  // 8. Run Event Loop
  mgr.Run();
}
```

---

## 🧩 ePIC & Podio Integration

### 1. Truth Matching
In EDM4hep, Rec-to-Truth links are stored in separate association arrays. `ePICReaction` traverses these arrays automatically. By providing a "Truth Role ID" when defining your particle candidates (e.g., `SetParticleCandidates("ele", 6, ...)`), the framework will automatically generate truth-matching flags like `TruthMatchedCombi()`.

### 2. Auxiliary Detector Data
ePIC-RAD introduces a robust Fluent Builder API (`ePICAssociationManager`) to seamlessly unpack nested PODIO One-To-Many relations.

Instead of dealing with raw 2D arrays, you can extract specific detector hits (like cluster energies or track times) and project them directly onto your physics tracks:
```cpp
epic::ePICAssociationManager assoc(rad_df);

assoc.For("Central")
     .From({"EcalBarrelClusters", "EcalEndcapPClusters"})
     .Extract("energy")
     .As("cal_energy");

assoc.Build(); 
```
This automatically produces a padded, flat RDataFrame column (`rec_cal_energy`) that perfectly synchronizes with the combinatorial engine, using `NaN` safely where tracks lack detector hits. You can then pass this into your kinematics processor to correct momenta on the fly.

---

## 📁 Output Structure

The framework utilizes `SnapshotCombi` to safely output complex multi-particle data from multi-threaded execution environments using `ROOT::TBufferMerger`.

The resulting `TTree` is **flat** and strictly analysis-friendly:
* **Scalars:** Event-level variables (e.g., Beam Energy, Q2) are broadcasted and repeated for every valid combination in the event.
* **Vectors:** Combinatorial variables (e.g., specific particle momenta) are flattened so that 1 Tree Entry = 1 Combination hypothesis.