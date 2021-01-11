#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

using std::vector;

TH1D *h_dijet_jet_pT;
TH1D *h_dijet_jet_E;
TH1D *h_dijet_jet_eta;
TH1D *h_dijet_jet_m0;
TH1D *h_dijet_jet_nJet;
vector<float> vDijet_jet_pT_;
vector<float> vDijet_jet_m0_;
vector<float> vDijet_jet_eta_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet( TTree* tree, edm::Service<TFileService> &fs ) {

  h_dijet_jet_pT    = fs->make<TH1D>("h_jet_pT"  , "p_{T};p_{T};Particles", 300,  0., 1500.);
  h_dijet_jet_E     = fs->make<TH1D>("h_jet_E"   , "E;E;Particles"        , 100,  0., 800.);
  h_dijet_jet_eta   = fs->make<TH1D>("h_jet_eta" , "#eta;#eta;Particles"  , 100, -5., 5.);
  h_dijet_jet_nJet  = fs->make<TH1D>("h_jet_nJet", "nJet;nJet;Events"     ,  10,  0., 10.);
  h_dijet_jet_m0    = fs->make<TH1D>("h_jet_m0"  , "m0;m0;Events"         , 100,  0., 100.);

  tree->Branch("jetPt",  &vDijet_jet_pT_);
  tree->Branch("jetM",   &vDijet_jet_m0_);
  tree->Branch("jetEta", &vDijet_jet_eta_);

} // branchesEvtSel_jet_dijet()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(jetCollectionT_, jets);
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  vJetIdxs.clear();
  vDijet_jet_pT_.clear();
  vDijet_jet_m0_.clear();
  vDijet_jet_eta_.clear();

  int nJet = 0;
  if (isTTbar_) { //is a ttbar sample
    std::vector<TLorentzVector> had_tops,bdau,wdau;
    for (const auto & p : *genParticles.product()) {
      int id = p.pdgId();
      if ( std::abs(p.pt()) < minTopPt_ ) continue;
      if ( std::abs(p.eta()) > maxTopEta_) continue;
      if(abs(id) != 6 || p.numberOfDaughters()!=2) continue;
      int iw=-1;
      int ib=-1;
      if (abs(p.daughter(0)->pdgId())==24 && abs(p.daughter(1)->pdgId())==5) {
	iw=0;ib=1;
      }
      else if (abs(p.daughter(1)->pdgId())==24 && abs(p.daughter(0)->pdgId())==5) {
	iw=1;ib=0;
      }
      else continue;  
    
      if (debug) std::cout << " >> topPt: " << p.pt() << " topEta: " << p.eta() << " topE: " << p.energy() << std::endl; 
    
      const reco::Candidate *d = p.daughter(iw);
      const reco::Candidate *b = p.daughter(ib);
      while(d->numberOfDaughters() == 1) d = d->daughter(0);
      if(!(abs(d->daughter(0)->pdgId()) < 10 && abs(d->daughter(1)->pdgId()) < 10)) continue;
      TLorentzVector the_top,the_w,the_b;
      the_top.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
      the_w.SetPtEtaPhiE(d->pt(),d->eta(),d->phi(),d->energy());
      the_b.SetPtEtaPhiE(b->pt(),b->eta(),b->phi(),b->energy());
      had_tops.push_back(the_top);
      wdau.push_back(the_w);
      bdau.push_back(the_b);
    }

    // Loop over jets
    for ( unsigned ihad=0;ihad<had_tops.size();ihad++) {
      float drMin = 999;
      reco::PFJetRef theJet;
      int indJet = -1;
      for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
	reco::PFJetRef iJet( jets, iJ );
	TLorentzVector vjet;
	vjet.SetPtEtaPhiE(iJet->pt(),iJet->eta(),iJet->phi(),iJet->energy());
	//if ( std::abs(iJet->pt()) < minJetPt_ ) continue;
	//if ( std::abs(iJet->eta()) > maxJetEta_) continue;
	if (had_tops[ihad].DeltaR(vjet)<drMin) {
          drMin = had_tops[ihad].DeltaR(vjet);
          indJet = iJ;
          theJet = iJet;
        }
      } // Jets
      if (drMin>0.6) continue;
      //if (wdau[ihad].DeltaR(vjet)>0.8) continue;
      //if (bdau[ihad].DeltaR(vjet)>0.8) continue;
      
      if ( debug ) std::cout << " >> top[" << ihad << "] E: " << had_tops[ihad].E() << std::endl;
      if ( debug ) std::cout << " >> jet[" << indJet << "] Pt: " << theJet->pt() << " jetEta: " << theJet->eta() << " jetE: " << theJet->energy() << " jetM: " << theJet->mass() << std::endl;

      vJetIdxs.push_back(indJet);
      nJet++;
      //break; This should allow two hadronic tops
      if ( (nJets_ > 0) && (nJet >= nJets_) ) break;
    } // hadronic tops
  } // isTTbar

  else { //is QCD
    std::vector<TLorentzVector> qcds;
    for (const auto & p : *genParticles.product()) { 
      int id = p.pdgId();
      if ( std::abs(p.pt()) < minTopPt_ ) continue;
      if ( std::abs(p.eta()) > maxTopEta_) continue;
      if(abs(id)==1 || abs(id)==2 || abs(id)==3 || abs(id)==4 || abs(id)==5 || abs(id)==21) {
        if ( !p.isHardProcess() ) continue;
        TLorentzVector the_qcd;
	the_qcd.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
        qcds.push_back(the_qcd);
      }
    }
    for ( unsigned ihad=0;ihad<qcds.size();ihad++) {
      float drMin = 999;
      reco::PFJetRef theJet;
      int indJet = -1;
      for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
	reco::PFJetRef iJet( jets, iJ );
        TLorentzVector vjet;
        vjet.SetPtEtaPhiE(iJet->pt(),iJet->eta(),iJet->phi(),iJet->energy());
        //if ( std::abs(iJet->pt()) < minJetPt_ ) continue;                                                
        //if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;                        
        if (qcds[ihad].DeltaR(vjet)<drMin) {
	  drMin = qcds[ihad].DeltaR(vjet);
	  indJet = iJ;
	  theJet = iJet;
	}
      } //Jets
      if ( drMin > 0.6 ) continue;
      if ( debug ) std::cout << " >> jet[" << indJet << "]Pt:" << theJet->pt() << " jetE:" << theJet->energy() << " jetM:" << theJet->mass() << std::endl;
      vJetIdxs.push_back(indJet);
      nJet++;
      if ( (nJets_ > 0) && (nJet > nJets_) ) break;  //  changed greater than or equal to to greater than                           
    } //QCD loop

  } // is QCD

  if ( debug ) {
    for(int thisJetIdx : vJetIdxs)
      std::cout << " >> vJetIdxs:" << thisJetIdx << std::endl;
  }

  if ( (nJets_ > 0) && (nJet > nJets_) ){ // changed not equal to greater than
    if ( debug ) std::cout << " Fail jet multiplicity:  " << nJet << " < " << nJets_ << std::endl;
    return false;
  }

  if ( vJetIdxs.size() == 0){
    if ( debug ) std::cout << " No passing jets...  " << std::endl;
    return false;
  }

  if ( debug ) std::cout << " >> has_jet_dijet: passed" << std::endl;
  return true;

} // runEvtSel_jet_dijet() 

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_dijet_jet_nJet->Fill( vJetIdxs.size() );
  // Fill branches and histograms 
  for(int thisJetIdx : vJetIdxs){
    reco::PFJetRef thisJet( jets, thisJetIdx );
    if ( debug ) std::cout << " >> Jet[" << thisJetIdx << "] Pt:" << thisJet->pt() << std::endl;
    h_dijet_jet_pT->Fill( std::abs(thisJet->pt()) );
    h_dijet_jet_E->Fill( thisJet->energy() );
    h_dijet_jet_m0->Fill( thisJet->mass() );
    h_dijet_jet_eta->Fill( thisJet->eta() );
    vDijet_jet_pT_.push_back( std::abs(thisJet->pt()) );
    vDijet_jet_m0_.push_back( thisJet->mass() );
    vDijet_jet_eta_.push_back( thisJet->eta() );
  }

} // fillEvtSel_jet_dijet()
