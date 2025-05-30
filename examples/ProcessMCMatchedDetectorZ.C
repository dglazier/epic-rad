//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "ePICDetectorReaction.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"

#include <TBenchmark.h>
#include <TCanvas.h>

void ProcessMCMatchedDetectorZ(){
   using namespace rad::names::data_type; //for Rec(), Truth()

  auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

  // ROOT::EnableImplicitMT(4);
  gBenchmark->Start("df");
  
  rad::config::ePICDetectorReaction epic{"events", "/home/dglazier/EIC/data/sim/jpac_z3900_10x100/AB_jpac_z3900_10x100_*_.recon.root"};
  epic.SetBeamsFromMC();

 // rad::config::ePICDetectorReaction epic{"events", "/home/dglazier/EIC/data/sim/jpac_z3900_10x100/AB_jpac_z3900_10x100_0_.recon.root"};
  epic.AliasColumnsAndMatchWithMC();
  //epic.AliasColumns();
   
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  //epic.setBeamIonIndex(rad::beams::InitBotFix());
  //epic.setBeamElectronIndex(rad::beams::InitTopFix());
  epic.setScatElectronIndex(4);

 
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  epic.setParticleIndex("ele",0,11);
  epic.setParticleIndex("pos",1,-11);
  epic.setParticleIndex("pip",2,211);
  epic.setParticleIndex("n",3,2112);
  
  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  
  epic.Particles().Sum("Jpsi",{"ele","pos"});
  epic.Particles().Sum("Zc",{"pip","Jpsi"});
  
  epic.setMesonParticles({"Jpsi","pip"});
  
  //can also add missing particles
  //And use those in calculations
  //Miss(name,{other final state particles})
  epic.Particles().Miss("calc_n",{rad::names::ScatEle().data(),"Zc"});
  epic.setBaryonParticles({"calc_n"});

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

 //////////////////////////////////////////////////////////
  // Now define calculated variables
  // Note reconstructed variables will have rec_ prepended
  // truth variables will have tru_ prepended
  //////////////////////////////////////////////////////////

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  rad::rdf::Mass(epic,"Whad","{Zc,n}");
  rad::rdf::Mass(epic,"JMass","{Jpsi}");
  rad::rdf::Mass(epic,"ZMass","{Zc}");
  rad::rdf::Mass(epic,"MissNMass","{calc_n}");
  rad::rdf::Q2(epic,"Q2");

  //t distribution, column name
  rad::rdf::TTop(epic,"t_gZ");
  rad::rdf::TBot(epic,"t_pn");
  rad::rdf::TPrimeBot(epic,"tp_pn");
  rad::rdf::TPrimeTop(epic,"tp_gZ");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");

  //exlusivity
  rad::rdf::MissMass(epic,"MissMass","{scat_ele,n,Zc}");
  rad::rdf::MissP(epic,"MissP_Meson","{scat_ele,Zc}");
  rad::rdf::MissPt(epic,"MissPt_Meson","{scat_ele,Zc}");
  rad::rdf::MissPz(epic,"MissPz_Meson","{scat_ele,Zc}");
  rad::rdf::MissTheta(epic,"MissTheta_Meson","{scat_ele,Zc}");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  rad::histo::Histogrammer histo{"set1",epic};
  // we can create many histograms by splitting events into
  // bins, where the same histogram is produced for the given bin
  // e.g. create 10 bins in tru_W between 4 and 54 GeV 
  histo.Splitter().AddRegularDimension(Truth()+"W", rad::histo::RegularSplits(10,4,14) );
  //can add as many split dimensions as we like
  //histo.Splitter().AddRegularDimension("xxx", rad::histo::RegularSplits(nbins,low,high) );
  histo.Init({Rec(),Truth()});//will create histograms for rec and truth

  histo.Create<TH1D,double>({"Q2","Q2",500,0,2.},{"Q2"});
  histo.Create<TH1D,double>({"W","W",100,0,20.},{"W"});
  histo.Create<TH1D,double>({"MesonMass","M(e-,e+, #pi) [GeV]",100,.3,5.},{"ZMass"});
  histo.Create<TH1D,double>({"MissMass","Mmiss [GeV]",1000,-10,10},{"MissMass"});
  histo.Create<TH1D,double>({"JMass","M(e-,e+) [GeV]",100,.3,5.},{"JMass"});
  histo.Create<TH1D,double>({"tpn","t(p,n) [GeV^{2}]",100,-2,5},{"t_pn"});
  histo.Create<TH1D,double>({"tgZ","t(g,Z) [GeV^{2}]",100,-2,5},{"t_gZ"});
  histo.Create<TH1D,double>({"cthCM","cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"phCM","#phi_{CM}",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  histo.Create<TH1D,double>({"missP","p_{miss}(e',Z)",105,0,105},{"MissP_Meson"});
  histo.Create<TH1D,double>({"missPt","p_{t,miss}(e',Z)",100,0,10},{"MissPt_Meson"});
  histo.Create<TH1D,double>({"missPz","p_{z,miss}(e',Z)",105,0,105},{"MissPz_Meson"});
  histo.Create<TH1D,double>({"missTheta","#theta_{miss}(e',Z)",100,0,1},{"MissTheta_Meson"});

  histo.Create<TH1D,double>({"EleP","p_{e-})",100,0,20},{"pmag[ele]"});
  histo.Create<TH1D,float>({"Ecal","Cluster Energy",100,0,10},{"clustersenergy[ele]"});
  histo.Create<TH2D,double,float>({"EleP_v_Ecal","Electron momentum v Cluster Energy",100,0,10,100,0,20},{"pmag[ele]","clustersenergy[ele]"});

  gBenchmark->Start("processing");
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
 
  //Draw all meson masss histograms on 1 canvas
  histo.DrawAll("MesonMass");
  histo.DrawAll("Ecal");
   
  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");

  //save all histograms to file
  histo.File("MCMatchedDetectorZ_hists.root");

  
  gBenchmark->Stop("df");
  gBenchmark->Print("df");
  
  gBenchmark->Start("snapshot");
  //  epic.Snapshot("MCMatchedDetectorZ.root");
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");

}
