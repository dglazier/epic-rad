//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "ePICDetectorReaction.h"
#include "ePICParticleCreator.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>


void PythiaRho(){
  
  using namespace rad::names::data_type; //for Rec(), Truth()
  gBenchmark->Start("df");

  //create reaction dataframe
  rad::config::ePICDetectorReaction epic{"events","/home/dglazier/eic_shell/benchmarks/reco/epic_craterlake_10x100/lowQ2_rho_10x100_1_0001.edm4eic.root"};
  
  //Note the following contains approx 75k rho events and takes 2hours
 // rad::config::ePICDetectorReaction epic{"events","root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.05.0/epic_craterlake/SIDIS/pythia6-eic/1.0.0/10x100/q2_0to1/pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run7.ab.00*.eicrecon.edm4eic.root"}; //Take the beam energy and angle from the MCParticles branch

 
//Here we actually take the mean over all events,
  //as this is likely to be what we have in the experiment
  epic.SetBeamsFromMC();

  epic.AliasColumnsAndMatchWithMC();
  //epic.AliasColumnsMC();

 
  epic.setScatElectronIndex(rad::indice::useNthOccurance(1,11),{"tru_pid"});
  //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  epic.setParticleIndex("pip",rad::indice::useNthOccurance(1,211),{"tru_pid"});
  epic.setParticleIndex("pim",rad::indice::useNthOccurance(1,-211),{"tru_pid"});
  epic.setParticleIndex("prot",rad::indice::useNthOccurance(1,2212),{"tru_pid"});

  //particle creator
  rad::epic::ePICParticleCreator p_creator{epic};
  //Get Farforward protons into stream, note name must match particle
  //p_creator.MCMatchedFarForwardProton("prot");
  p_creator.MCMatchedRomanPotProton("prot");
  //p_creator.MCMatchedB0Proton("prot");
 
  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  epic.setBaryonParticles({"prot"});
  //epic.setBaryonParticles({});

  epic.Particles().Sum("rho",{"pip","pim"});
  epic.setMesonParticles({"pip","pim"});

  //Also calculate proton as a missing particle
  p_creator.Miss("calc_prot",{rad::names::ScatEle().data(),"rho"});

  //must call this after all particles are configured
  epic.makeParticleMap();
  
  //////////////////////////////////////////////////////////////////
  ///Add some detector associations
  //////////////////////////////////////////////////////////////////
  
  epic.AssociateClusters({"EcalBarrelClusters","EcalBarrelImagingClusters",
      "EcalBarrelScFiClusters",
      "EcalEndcapNClusters","EcalEndcapPClusters","EcalEndcapPInsertClusters",
      "HcalBarrelClusters","HcalEndcapNClusters","LFHCALClusters",
      "EcalFarForwardZDCClusters","HcalFarForwardZDCClusters"},
    {"energy"}); //just going to associate the cluster energy from the list of detectors
  epic.AssociateTracks({"TaggerTrackerTracks"},
   		       {"momentum.z"});

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  rad::rdf::Mass(epic,"RhoMass","{pip,pim}");
  rad::rdf::Mass(epic,"DppMass","{pip,prot}");
  rad::rdf::Mass(epic,"D0Mass","{pim,prot}");
  
  //t distribution, column name
  rad::rdf::TTop(epic,"t_top");
  rad::rdf::TBot(epic,"t_bot");
  rad::rdf::TPrimeTop(epic,"tp_top");
  rad::rdf::TPrimeBot(epic,"tp_bot");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");
  rad::rdf::Q2(epic,"Q2");
  
  //decay angles
  rad::rdf::gn2s0s0s12::HelicityAngles(epic,"Heli");
  
  //exlusivity
  rad::rdf::MissMass(epic,"MissMass_Meson","{scat_ele,rho}");
  rad::rdf::MissMass2(epic,"MissMass2_Meson","{scat_ele,rho}");
  rad::rdf::MissP(epic,"MissP_Meson","{scat_ele,rho}");
  rad::rdf::MissPt(epic,"MissPt_Meson","{scat_ele,rho}");
  rad::rdf::MissPz(epic,"MissPz_Meson","{scat_ele,rho}");
  rad::rdf::MissTheta(epic,"MissTheta_Meson","{scat_ele,rho}");
  
  rad::rdf::MissMass(epic,"MissMass_All","{scat_ele,pip,pim,prot}");
  rad::rdf::MissMass2(epic,"MissMass2_All","{scat_ele,pip,pim,prot}");
  rad::rdf::MissP(epic,"MissP_All","{scat_ele,prot,rho}");
  rad::rdf::MissPt(epic,"MissPt_All","{scat_ele,prot,rho}");
  rad::rdf::MissPz(epic,"MissPz_All","{scat_ele,prot,rho}");
  rad::rdf::MissTheta(epic,"MissTheta_All","{scat_ele,prot,rho}");

  rad::rdf::DeltaPhi(epic,"DeltaPhi_Prot","{prot,calc_prot}");
  rad::rdf::DeltaTheta(epic,"DeltaTheta_Prot","{prot,calc_prot}");
  rad::rdf::DeltaP(epic,"DeltaP_Prot","{prot,calc_prot}");
  epic.Define("rec_DeltaMagP_Prot","rec_pmag[prot]-rec_pmag[calc_prot]");
  epic.Define("tru_DeltaMagP_Prot","tru_pmag[prot]-tru_pmag[calc_prot]");

 /*
  */
    //Create a column tagger_pz for particles "scat_ele" based on "rec_tracks_momentum_z" (which is a vector over all particles) and define the variable as the Mean of the entries associated with scat_ele
  epic.DefineForParticles("tagger_pz",{"scat_ele"},{"-rec_tracks_momentum_z"},"rad::helpers::First");
  epic.DefineForParticles("cal_energy_sum",{"scat_ele"},{"rec_clusters_energy"},"rad::helpers::Sum");

  //compare scat_ele calorimeter energy and track momenta to truth values
  epic.Define("Delta_cal_energy_sum","(cal_energy_sum_scat_ele-tru_pmag[scat_ele])/tru_pmag[scat_ele]");
  epic.Define("Delta_momentum_scat_ele","res_pmag[scat_ele]");
 
  //exclusivity cut on truth variables
  //i.e. make sure we are analysing exclusive 2pi final state
  auto epic_exc= epic;
  epic_exc.Filter("(tru_Npip==1)*(tru_Npim==1)","reaction_filter");
  epic_exc.Filter("tru_MissMass_All<0&&tru_MissMass_Meson<1.1","exclusivity_filter");
  epic_exc.Filter("(rec_Npip==1)*(rec_Npim==1)*(scat_ele>-1)*(rec_pmag[scat_ele]>2)","rec_filter");
  // ///////////////////////////////////////////////////////////
  // //Define histograms
  // ///////////////////////////////////////////////////////////
  std::vector<rad::config::ePICDetectorReaction*> epics {&epic,&epic_exc};
  std::vector<rad::histo::Histogrammer> histos;
  histos.reserve(epics.size());
  
  for(auto* epi:epics){
    ///////////////////////////////////////////////////////////
    //Draw histograms
    ///////////////////////////////////////////////////////////
    rad::histo::Histogrammer histo{"set1",*epi};
    histo.SetVerbose(1);
    // we can create many histograms by splitting events into
    // bins, where the same histogram is produced for the given bin
    // e.g. create 10 bins in tru_W between 4 and 54 GeV 
    //histo.Splitter().AddRegularDimension(Truth()+"W", rad::histo::RegularSplits(10,4,14) );
    //can add as many split dimensions as we like
    //histo.Splitter().AddRegularDimension("xxx", rad::histo::RegularSplits(nbins,low,high) );
    histo.Init({Rec(),Truth()});//will create histograms for rec and truth
  
    histo.Create<TH1D,double>({"hQ2",";Q^{2} [GeV^{2}]",100,0,0.2},{"Q2"});
    histo.Create<TH1D,double>({"hWelec",";W (electro miss mass) [GeV/c^{2}]",100,0,200.},{"W"});
    histo.Create<TH1D,double>({"hRhoMass",";M(#pi-,#pi+) [GeV/c^{2}]",100,.2,3.},{"RhoMass"});
    histo.Create<TH1D,double>({"hD0Mass",";M(#pi-p) [GeV/c^{2}]",100,1,3.},{"D0Mass"});
    histo.Create<TH1D,double>({"hDppMass",";M(#pi+p) [GeV/c^{2}]",100,1,3.},{"DppMass"});
  
    histo.Create<TH1D,double>({"httop",";t(top vertex) [GeV^{2}]",100,-1,5},{"t_top"});
    histo.Create<TH1D,double>({"htbot",";t(bottom vertex) [GeV^{2}]",100,-1,5},{"t_bot"});
    histo.Create<TH1D,double>({"htptop",";t'(top vertex) [GeV^{2}]",100,-1,5},{"tp_top"});
    histo.Create<TH1D,double>({"htpbot",";t'(bottom vertex) [GeV^{2}]",100,-1,5},{"tp_bot"});
  
    histo.Create<TH1D,double>({"hcthCM",";cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
    histo.Create<TH1D,double>({"hphCM",";#phi_{CM}",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  
    histo.Create<TH1D,double>({"hMissMass_Meson",";M^{2}_{miss} [GeV/c^{2}]",1000,-50,50},{"MissMass2_Meson"});
    histo.Create<TH1D,double>({"hmissP_Meson",";p_{miss}(e',Y)",1000,80,120},{"MissP_Meson"});
    histo.Create<TH1D,double>({"hmissPt_Meson",";p_{t,miss}(e',Y)",100,0,2},{"MissPt_Meson"});
    histo.Create<TH1D,double>({"hmissPz_Meson",";p_{z,miss}(e',Y)",1000,80,120},{"MissPz_Meson"});
    histo.Create<TH1D,double>({"hmissTheta_Meson",";#theta_{miss}(e',Y)",100,-0.1,TMath::Pi()+0.1},{"MissTheta_Meson"});

    histo.Create<TH1D,double>({"hMissMass_All",";M^{2}_{miss} [GeV/c^{2}]",1000,-10,10},{"MissMass_All"});
    histo.Create<TH1D,double>({"hmissP_All",";p_{miss}(e',Y)",1000,-10,10},{"MissP_All"});
    histo.Create<TH1D,double>({"hmissPt_All",";p_{t,miss}(e',Y)",100,0,2},{"MissPt_All"});
    histo.Create<TH1D,double>({"hmissPz_All",";p_{z,miss}(e',Y)",1000,-10,10},{"MissPz_All"});
    histo.Create<TH1D,double>({"hmissTheta_All",";#theta_{miss}(e',Y)",100,-0.1,TMath::Pi()+0.1},{"MissTheta_All"});

    histo.Create<TH1D,double>({"hDeltaPhi_Prot",";#Delta#phi_{p} [rad]",200,-3.2,3.2},{"DeltaPhi_Prot"});
    histo.Create<TH2D,double,double>({"hDeltaPhi_ProtVPhi",";#Delta#phi_{p} [rad] v #phi",100,-1.2,1.2,100,-3.2,3.2},{"DeltaPhi_Prot","phi[calc_prot]"});
    histo.Create<TH1D,double>({"hDeltaTheta_Prot",";#Delta#theta_{p} []",200,-0.001,0.02},{"DeltaTheta_Prot"});
    histo.Create<TH1D,double>({"hDeltaP_Prot",";#Delta P_{p} [GeV/c]",200,-1,10},{"DeltaP_Prot"});
    histo.Create<TH1D,double>({"hDeltaMagP_Prot",";#Delta P_{p} [GeV/c]",200,-10,10},{"DeltaMagP_Prot"});

  
    //particle momenta
    histo.Create<TH1D,double>({"hpmag_elec",";p_{e'} [GeV/c]",100,-1,20},{"pmag[scat_ele]"});
    histo.Create<TH1D,double>({"heta_elec",";#eta_{e'} ",100,-10,10},{"eta[scat_ele]"});
    histo.Create<TH1D,double>({"htheta_elec",";#theta_{e'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"theta[scat_ele]"});
    histo.Create<TH1D,double>({"hphi_elec",";#phi_{e'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"phi[scat_ele]"});
    histo.Create<TH2D,double,ROOT::RVecF>({"EleP_v_Ecal","Electron momentum v Cluster Energy",100,0,10,100,0,20},{"pmag[scat_ele]","clusters_energy[scat_ele]"});

    histo.Create<TH1D,double>({"hpmag_pprime",";p_{p'} [GeV/c]",100,90,110},{"pmag[prot]"});
    histo.Create<TH1D,double>({"heta_pprime",";#eta_{p'} ",100,4,10},{"eta[prot]"});
    histo.Create<TH1D,double>({"htheta_pprime",";#theta_{p'} [rad]",200,0,0.02},{"theta[prot]"});
    histo.Create<TH1D,double>({"hphi_pprime",";#phi_{p'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"phi[prot]"});

    histo.Create<TH1D,double>({"hpmag_pcalc",";p_{p'} [GeV/c]",100,90,110},{"pmag[calc_prot]"});
    histo.Create<TH1D,double>({"heta_pcalc",";#eta_{p'} ",100,4,10},{"eta[calc_prot]"});
    histo.Create<TH1D,double>({"htheta_pcalc",";#theta_{p'} [rad]",200,0,0.02},{"theta[calc_prot]"});
    histo.Create<TH1D,double>({"hphi_pcalc",";#phi_{p'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"phi[calc_prot]"});
   
    histo.Create<TH2D,double,float>({"EleP_v_TrackP","Electron momentum v track pz",100,0,10,100,0,20},{"pmag[scat_ele]","tagger_pz_scat_ele"});
    histo.Create<TH1D,float>({"TrackP","track pz",100,0,20},{"tagger_pz_scat_ele"});

    histo.Create<TH2D,double,double>({"truEleP_v_DeltaCalEsum","Electron momentum v (Cluster Energysum - momentum)",100,0,10,100,-1,1},{"pmag[scat_ele]","Delta_cal_energy_sum"});

    histo.Create<TH1D,double>({"DeltaCalEsum","(Cluster Energysum - TruMomentum)",200,-1,1},{"Delta_cal_energy_sum"});
    histo.Create<TH1D,double>({"DeltaTrackP","(RecMomentum - TruMomentum)",200,-1,1},{"Delta_momentum_scat_ele"});

    histos.emplace_back( std::move(histo) );
  }
 
    // rad::config::PrintDefinedColumnNames(epic.CurrFrame());
    //Save all calculations to a root tree 
    //epic.BookLazySnapshot("PythiaRho_10x100.root");
    //save all histograms to file
    //histo.DrawAll("tagger_pz");
  // rad::histo::Histogrammer histo{"set1",epic};
  // histo.SetVerbose(1);
  // histo.Init({Rec(),Truth()});//will create histograms for rec and truth
  // histo.Create<TH1D,double>({"hQ2",";Q^{2} [GeV^{2}]",100,0,0.2},{"Q2"});
  // histo.File("test.root");
  
    gBenchmark->Start("processing");
    histos[0].File("PythiaRho_hists_Inc.root");
    histos[1].File("PythiaRho_hists_Exc.root");

    gBenchmark->Stop("processing");
    gBenchmark->Print("processing");
    gBenchmark->Stop("df");
    gBenchmark->Print("df");

 
  }
