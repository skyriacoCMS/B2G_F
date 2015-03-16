#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

// dR and dPhi                                                                                                        

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// some basic root definitions
#include <TLorentzVector.h>
#include <vector>
#include <TROOT.h>


// Photon dependencies
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" // gives access to the (release cycle dependent) trigger object codes

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

// Fastjet (for creating subjets)                                          
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>



// Jets
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"





#include<vector>

using namespace reco;
using namespace edm;
using namespace std;
using namespace fastjet;
using namespace trigger;

typedef std::vector<pat::Jet> PatJetCollection;

class  PhotonJets : public edm::EDProducer {
public:
  PhotonJets( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  bool isMatchedWithTrigger();
  bool passIDWP();
  float IsoCalc();

  InputTag phoLabel_;
  InputTag pvLabel_, convLabel_;
  InputTag rhoLabel_;
  InputTag pckPFCdsLabel_;
  InputTag jLabel_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
  // edm::EDGetTokenT<std::vector<pat::Jet> > jetToken_;
  
  InputTag triggerResultsLabel_, triggerSummaryLabel_;
  InputTag hltElectronFilterLabel_;
  TString hltPath_;
  HLTConfigProvider hltConfig;
  int triggerBit;
  int debug_; 
  
};



//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}



//Necessary for the  Isolation Rho Correction
namespace EffectiveAreas {
  const int nEtaBins = 7;
  const float etaBinLimits[nEtaBins+1] = {
    0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};

  const float areaPhotons[nEtaBins] = {
    0.0894, 0.0750, 0.0423, 0.0561, 0.0882, 0.1144, 0.1684
  };
  const float areaNeutralHadrons[nEtaBins] = {
    0.049, 0.0108, 0.0019, 0.0037, 0.0062, 0.0130, 0.1699
  };
  const float areaChargedHadrons[nEtaBins] = {
    0.0089, 0.0062, 0.0086, 0.0041, 0.0113, 0.0085, 0.0039
  };
}

// The Cut based Photon ID's - the current is PHYS 14 first Iteration

namespace CPID_B{

  // this ID is the PHYS14 Iteration 1 for Barrel
  const float H_o_E[3]={0.032,0.020,0.012};
  const float s_IEIE[3]={0.01000,0.0099,0.0098};
  const float isoC[3]={2.94,2.62,1.91};
  const float isoN[3]={3.16,2.69,2.55};
  const float isoP[3]={4.43,1.35,1.29};
  const float slope_n = 0.0023;
  const float slope_p = 0.0004;
  

}


namespace CPID_E{

  // this ID is the PHYS14 Iteration 1 for EndCap
  const float H_o_E[3]={0.023,0.011,0.011};
  const float s_IEIE[3]={0.0270,0.0269,0.0264};
  const float isoC[3]={3.07,1.40,1.26};
  const float isoN[3]={17.16,4.92,2.71};
  const float isoP[3]={2.11,2.11,1.91};
  const float slope_n = 0.0116;
  const float slope_p = 0.0037;
  
}



PhotonJets::PhotonJets(const edm::ParameterSet& iConfig):
   phoLabel_(iConfig.getParameter<edm::InputTag>("phoLabel")),
   pvLabel_(iConfig.getParameter<edm::InputTag>("pv")),   // "offlinePrimaryVertex"
   rhoLabel_(iConfig.getParameter<edm::InputTag>("rho")), //rhofixedgridRhoFastjet All"
   pckPFCdsLabel_(iConfig.getParameter<edm::InputTag>("packedPFCands")),
   jLabel_       (iConfig.getParameter<edm::InputTag>("jetLabel")),
   ebReducedRecHitCollection_(consumes <EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))), //Lazy tool additions
   eeReducedRecHitCollection_(consumes <EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection")))  // Lazy tool additions  


{
  debug_ = iConfig.getUntrackedParameter<int>("debugLevel",int(0));
  produces<std::vector<pat::Jet> >();
}

void PhotonJets::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  
  //PV
  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByLabel(pvLabel_, vertices);


  //Photons
  edm::Handle<std::vector<pat::Photon> > phoHandle;
  iEvent.getByLabel(phoLabel_, phoHandle);
  auto_ptr<vector<pat::Photon> > phoColl( new vector<pat::Photon> (*phoHandle) );
  
  
  //Packed PF Cands
  edm::Handle<std::vector<pat::PackedCandidate>> pfCnd1Handle;
  iEvent.getByLabel(pckPFCdsLabel_,pfCnd1Handle); 
  auto_ptr<vector<pat::PackedCandidate> > CandColl( new vector<pat::PackedCandidate> (*pfCnd1Handle) );


  edm::Handle< edm::View<reco::Candidate>> pfCndHandle;
  iEvent.getByLabel(pckPFCdsLabel_,pfCndHandle);

  //Jet collection
  edm::Handle<std::vector<pat::Jet> > jetHandle, packedjetHandle;
  iEvent.getByLabel(jLabel_, jetHandle);
  auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );





  // Ecal Rec Hits

  //rho
  //  float rho_;
  edm::Handle< double > rhoH;
  iEvent.getByLabel(rhoLabel_,rhoH);
  //rho_ = *rhoH;


  if(debug_>=1) cout<<"vtx size " << vertices->size()<<endl; 

  reco::TrackBase::Point vtxPoint(0,0, 0);
  if(  vertices->size() >= 1 ) {
    vtxPoint = vertices->at(0).position();
  }

  if(debug_>=1) cout<<"vtxPoint " <<vtxPoint.X()<<" "<< vtxPoint.Y()<<" "<< vtxPoint.Z()<<endl; 
  
  // noZS::EcalClusterLazyTools *lazyToolnoZS;
  //lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_);

  std::cout<<"Eentered Event ok"<<endl;


  std::cout<<"ABout to enter the collection"<<endl;
   
  int ijet = -1; 
  for (size_t i = 0; i< jetColl->size(); i++){
    pat::Jet & jet = (*jetColl)[i];
    ijet ++;
    double jet_pt  = jet.pt(); 
    double jet_eta = jet.eta(); 
    double jet_phi = jet.phi(); 
    

    if(jet_pt < 200 ) continue;
       

    int pho_max_pt_indx = -99; 
    double pho_max_pt = 0; 
    //  int pho_sec_pt_indx = -99; 
    //    double pho_sec_pt = 0; 
    int ipho = -1; 
    for(size_t j = 0; j < phoColl->size();j++){
      pat::Photon & pho = (*phoColl)[j];
       double pho_pt  = pho.pt(); 
       double pho_eta = pho.eta(); 
       double pho_phi = pho.phi(); 
       ipho++;
          
       if(pho_pt < 15 || fabs(pho_eta) > 2.5 ) continue;
       
       
       double dr_pho_jet = deltaR(pho_eta,pho_phi,jet_eta,jet_phi);
       
       if(dr_pho_jet > 0.8 ) continue; 
       
       if(pho_pt > pho_max_pt){
	 pho_max_pt_indx = j;
	 pho_max_pt  = pho_pt; 
       }
    }//EOF Loop over gammas
    
    // hardwired nee to change!!!!!!!!1 for debuggin now only     
    // pho_max_pt_indx = 0;
    
    // select the highest pt photon in the jet for now
    float spt0 = -99;    float spt1 = -99;    float spt2 = -99;
    float sphi0 = -99;    float sphi1 = -99;    float sphi2 = -99;
    float seta0 = -99;    float seta1 = -99;    float seta2 = -99;
    float sen0 = -99;    float sen1 = -99;    float sen2 = -99;

    float photon_subjet_frac = -99;
   

    if(pho_max_pt_indx != -99  && phoColl->size() > 0 ){    
    pat::Photon & pho = (*phoColl)[pho_max_pt_indx];
    std::vector<fastjet::PseudoJet> FJparticles;
    TLorentzVector Pho_Sums; 
   

    Pho_Sums.SetPtEtaPhiE(0,0,0,0);
    for (unsigned int k = 0; k < jet.numberOfDaughters(); k++){
	const edm::Ptr<reco::Candidate> & this_constituent = jet.daughterPtr(k);

	

	bool inFootprint = isInFootprint(pho.associatedPackedPFCandidates(), this_constituent);
       	if(inFootprint)	  continue;
	FJparticles.push_back( fastjet::PseudoJet(this_constituent->px(), this_constituent->py(), this_constituent->pz(), this_constituent->energy()) );
	

    }
    
    fastjet::PseudoJet PhoC(pho.px(),pho.py(),pho.pz(),pho.energy());
    
    PhoC.set_user_index(0456);
    FJparticles.push_back( PhoC );

    
    fastjet::JetDefinition jet_def_3(fastjet::antikt_algorithm,1.5708); 
    fastjet::ClusterSequence clust_seq_3(FJparticles, jet_def_3);
   
    int nSubJets = 3;
    std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(clust_seq_3.exclusive_jets_up_to(nSubJets));


     if(subjets.size() == 3){ 
      int sub_wgamma = -99; 
      spt0 = subjets[0].perp();  spt1 = subjets[1].perp();  spt2 = subjets[2].perp();
      sphi0 = subjets[0].phi();  sphi1 = subjets[1].phi();  sphi2 = subjets[2].phi();
      seta0 = subjets[0].eta();  seta1 = subjets[1].eta();  seta2 = subjets[2].eta();
      sen0 = subjets[0].e();     sen1 = subjets[1].e();     sen2 = subjets[2].e();
      for(int k = 0 ; k < nSubJets; k++){
       	
	vector<fastjet::PseudoJet> Sconst = subjets[k].constituents();
	for(unsigned int jconst = 0; jconst < Sconst.size(); jconst++){
	  if(Sconst[jconst].user_index() == 0456){
	    sub_wgamma = k;
	    break;
	  }
	}//eof subjet const loop 
      }//eof subjets loop  
      
      photon_subjet_frac = pho_max_pt/subjets[sub_wgamma].pt();
      
     }
     
     
    }
    
    

   
    jet.addUserInt("jetIndex",ijet);
    jet.addUserInt("phoIndex",pho_max_pt_indx);
    jet.addUserFloat("phoSubjetPtFrac",photon_subjet_frac);

    
    jet.addUserFloat("SubPt0", spt0);
    jet.addUserFloat("SubPt1", spt1);
    jet.addUserFloat("SubPt2", spt2);
    jet.addUserFloat("SubEta0",seta0);
    jet.addUserFloat("SubEta1",seta1);
    jet.addUserFloat("SubEta2",seta2); 
    jet.addUserFloat("SubPhi0",sphi0);
    jet.addUserFloat("SubPhi1",sphi1);
    jet.addUserFloat("SubPhi2",sphi2);
    jet.addUserFloat("SubEne0",sen0);
    jet.addUserFloat("SubEne1",sen1);
    jet.addUserFloat("SubEne2",sen2);
   
    
    
  }//EOF Loop over Jets
  
  /*
  for(size_t i = 0; i < phoColl->size();i++){
    pat::Photon & pho = (*phoColl)[i];
    
    float pho_eta = pho.superCluster()->eta();
    float pho_phi = pho.superCluster()->phi();
    float pho_pt  = pho.pt();

    if(pho_pt < 15 || fabs(pho_eta) < 2.5) continue;

    reco::Vertex pv = vertices->at(0);
    //setting the photon directions with respect to the primary vertex
    math::XYZVector photon_directionWrtVtxs(pho.superCluster()->x() - pv.x(),
					    pho.superCluster()->y() - pv.y(),
					    pho.superCluster()->z() - pv.z());
    

    // shower shape variables
    float r_9  = pho.r9();
    float hoe = pho.hadTowOverEm();
  
    //extracting sigma I eta I eat 5X5 using the Lazy tool
    const auto& theseed = *(pho.superCluster()->seed());
    float see = -999;
    std::vector<float> vCov = lazyToolnoZS->localCovariances( theseed );
    see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    float sigmaIetaIeta = see;

    //isolation variables 

    //Calculating the Isolation from the pf candidates
    //CandColl

    float iso_ch = 0; 
    float iso_n = 0; 
    float iso_p = 0; 
   
    const float coneSizeDR = 0.3; 
    const float dxyMax     = 0.1; 
    const float dzMax      = 0.2; 


    for(size_t k = 0; k< CandColl->size();i++){
      pat::PackedCandidate & pfC = (*CandColl)[i];
      
      const auto& iCand = pfCndHandle->ptrAt(k);                                                          


      float DR = deltaR(photon_directionWrtVtxs.eta(),photon_directionWrtVtxs.phi(),pfC.eta(),pfC.phi());
      if(DR > coneSizeDR) continue;

	bool inFootprint = isInFootprint(pho.associatedPackedPFCandidates(), iCand);
	if(inFootprint) continue;
	
	// Now  we need to check the type of the Candidate and also associate the PF::h with the primary vertex
	
	reco::PFCandidate::ParticleType thisCandidateType = reco::PFCandidate::X;
	const int pdgId = pfC.pdgId();
        if( pdgId == 22 )
          thisCandidateType = reco::PFCandidate::gamma;
        else if( abs(pdgId) == 130) // PDG ID for K0L                                                                  
          thisCandidateType = reco::PFCandidate::h0;
        else if( abs(pdgId) == 211) // PDG ID for pi+-                                                                 
          thisCandidateType = reco::PFCandidate::h;

	if(thisCandidateType == reco::PFCandidate::h ){
	 float dz = pfC.pseudoTrack().dz(pv.position());
	 float dxy =pfC.pseudoTrack().dxy(pv.position());
	  if( dxyMax < dxy ) continue;
	  if( dzMax < dz ) continue;
	  iso_ch += pfC.pt();
	}
	
	if(thisCandidateType == reco::PFCandidate::h0    ) iso_n += pfC.pt();	
	if(thisCandidateType == reco::PFCandidate::gamma ) iso_p += pfC.pt();
    }



    



    // setting effective areas region
    int ieBin = 0; 
    while(ieBin < (EffectiveAreas::nEtaBins -1) && pho_eta > EffectiveAreas::etaBinLimits[ieBin+1]) ieBin++;

    float isoC_withEA = std::max(float(0.0),iso_ch - rho_ * EffectiveAreas::areaChargedHadrons[ieBin]);
    float isoN_withEA = std::max(float(0.0),iso_n - rho_ * EffectiveAreas::areaNeutralHadrons[ieBin]);
    float isoP_withEA = std::max(float(0.0),iso_p - rho_ * EffectiveAreas::areaPhotons[ieBin]);


    

    // other variables
    int hasPixelSeed    = pho.hasPixelSeed(); 


    // ID bools
    int isLoose  = 0; 
    int isMedium = 0; 
    int isTight  = 0; 
    
    if(fabs(pho_eta) < 1.479){
      if(hoe < CPID_B::H_o_E[0]  && sigmaIetaIeta < CPID_B::s_IEIE[0]  && isoC_withEA < CPID_B::isoC[0] &&  isoN_withEA  < CPID_B::isoN[0] + CPID_B::slope_n*pho_pt && isoP_withEA < CPID_B::isoP[0] + CPID_B::slope_p*pho_pt ) isLoose = 1;
      
      if(hoe < CPID_B::H_o_E[1]  && sigmaIetaIeta < CPID_B::s_IEIE[1]  && isoC_withEA < CPID_B::isoC[1] &&  isoN_withEA  < CPID_B::isoN[1] + CPID_B::slope_n*pho_pt && isoP_withEA < CPID_B::isoP[1] + CPID_B::slope_p*pho_pt ) isMedium = 1;
      if(hoe < CPID_B::H_o_E[2]  && sigmaIetaIeta < CPID_B::s_IEIE[2]  && isoC_withEA < CPID_B::isoC[2] &&  isoN_withEA  < CPID_B::isoN[2] + CPID_B::slope_n*pho_pt && isoP_withEA < CPID_B::isoP[2] + CPID_B::slope_p*pho_pt ) isTight = 1;
    }else{
      if(hoe < CPID_E::H_o_E[0]  && sigmaIetaIeta < CPID_E::s_IEIE[0]  && isoC_withEA < CPID_E::isoC[0] &&  isoN_withEA  < CPID_E::isoN[0] + CPID_E::slope_n*pho_pt && isoP_withEA < CPID_E::isoP[0] + CPID_E::slope_p*pho_pt ) isLoose = 1;
      if(hoe < CPID_E::H_o_E[1]  && sigmaIetaIeta < CPID_E::s_IEIE[1]  && isoC_withEA < CPID_E::isoC[1] &&  isoN_withEA  < CPID_E::isoN[1] + CPID_E::slope_n*pho_pt && isoP_withEA < CPID_E::isoP[1] + CPID_E::slope_p*pho_pt ) isMedium = 1;
      if(hoe < CPID_E::H_o_E[2]  && sigmaIetaIeta < CPID_E::s_IEIE[2]  && isoC_withEA < CPID_E::isoC[2] &&  isoN_withEA  < CPID_E::isoN[2] + CPID_E::slope_n*pho_pt && isoP_withEA < CPID_E::isoP[2] + CPID_E::slope_p*pho_pt ) isTight = 1;
    }
    


    
    pho.addUserInt("phoSceta",pho_eta);
    pho.addUserInt("phoScphi",pho_phi);
    pho.addUserInt("phopt",pho_pt);

    pho.addUserInt("hasPixelSeed",   hasPixelSeed);
    pho.addUserFloat("sigmaIetaIeta",    sigmaIetaIeta);
    pho.addUserFloat("hoe",     hoe);
    pho.addUserFloat("r9",     r_9);

    pho.addUserFloat("isoCwithEA",isoC_withEA);
    pho.addUserFloat("isoPwithEA",isoP_withEA);
    pho.addUserFloat("isoNwithEA",isoN_withEA);

    pho.addUserFloat("isLoose",    isLoose);
    pho.addUserFloat("isMedium",    isMedium);
    pho.addUserFloat("isTight",    isTight);

  }

  */
  //  iEvent.put( phoColl );
  iEvent.put( jetColl );

}

float
PhotonJets::IsoCalc(){

  return 2.0;

}


bool
PhotonJets::isMatchedWithTrigger(){
  return true;
}


bool PhotonJets::passIDWP(){
 
  return true;
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhotonJets);
