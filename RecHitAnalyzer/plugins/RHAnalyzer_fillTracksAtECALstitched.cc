#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill Tracks into stitched EEm_EB_EEp image //////////////////////
// Store all Track positions into a stitched EEm_EB_EEp image 

TH2F *hEvt_EE_tracksPt[nEE];
TH2F *hEvt_EE_tracksQPt[nEE];
TH2F *hEvt_EE_tracksPtAtECAL[nEE];
TH2F *hEvt_EE_tracksQPtAtECAL[nEE];
TH2F *hEvt_EE_tracksPtAtHCAL[nEE];
TH2F *hEvt_EE_tracksQPtAtHCAL[nEE];

TProfile2D *hECAL_tracks;
TProfile2D *hECAL_tracksPt;
TProfile2D *hECAL_tracksQPt;
TProfile2D *hECAL_tracksAtECAL;
TProfile2D *hECAL_tracksPtAtECAL;
TProfile2D *hECAL_tracksQPtAtECAL;
TProfile2D *hECAL_tracksAtHCAL;
TProfile2D *hECAL_tracksPtAtHCAL;
TProfile2D *hECAL_tracksQPtAtHCAL;

std::vector<float> vECAL_tracksPt_;
std::vector<float> vECAL_tracksQPt_;
std::vector<float> vECAL_tracksPtAtECAL_;
std::vector<float> vECAL_tracksQPtAtECAL_;
std::vector<float> vECAL_tracksPtAtHCAL_;
std::vector<float> vECAL_tracksQPtAtHCAL_;

// Initialize branches _______________________________________________________________//
void RecHitAnalyzer::branchesTracksAtECALstitched ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("ECAL_tracksPt",    &vECAL_tracksPt_);
  tree->Branch("ECAL_tracksQPt",    &vECAL_tracksQPt_);

  tree->Branch("ECAL_tracksPtAtECAL",    &vECAL_tracksPtAtECAL_);
  tree->Branch("ECAL_tracksQPtAtECAL",    &vECAL_tracksQPtAtECAL_);

  tree->Branch("ECAL_tracksPtAtHCAL",    &vECAL_tracksPtAtHCAL_);
  tree->Branch("ECAL_tracksQPtAtHCAL",    &vECAL_tracksQPtAtHCAL_);

  // Intermediate helper histogram (single event only)
  hEvt_EE_tracksPt[0] = new TH2F("evt_EEm_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksPt[1] = new TH2F("evt_EEp_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  hEvt_EE_tracksQPt[0] = new TH2F("evt_EEm_tracksQPt", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksQPt[1] = new TH2F("evt_EEp_tracksQPt", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );


  hEvt_EE_tracksPtAtECAL[0] = new TH2F("evt_EEm_tracksPtAtECAL", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksPtAtECAL[1] = new TH2F("evt_EEp_tracksPtAtECAL", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  hEvt_EE_tracksQPtAtECAL[0] = new TH2F("evt_EEm_tracksQPtAtECAL", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksQPtAtECAL[1] = new TH2F("evt_EEp_tracksQPtAtECAL", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );


  hEvt_EE_tracksPtAtHCAL[0] = new TH2F("evt_EEm_tracksPtAtHCAL", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksPtAtHCAL[1] = new TH2F("evt_EEp_tracksPtAtHCAL", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  hEvt_EE_tracksQPtAtHCAL[0] = new TH2F("evt_EEm_tracksQPtAtHCAL", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksQPtAtHCAL[1] = new TH2F("evt_EEp_tracksQPtAtHCAL", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );


  // Histograms for monitoring
  hECAL_tracks = fs->make<TProfile2D>("ECAL_tracks", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

  hECAL_tracksPt = fs->make<TProfile2D>("ECAL_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

  hECAL_tracksQPt = fs->make<TProfile2D>("ECAL_tracksQPt", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );


  hECAL_tracksAtECAL = fs->make<TProfile2D>("ECAL_tracksAtECAL", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

  hECAL_tracksPtAtECAL = fs->make<TProfile2D>("ECAL_tracksPtAtECAL", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

  hECAL_tracksQPtAtECAL = fs->make<TProfile2D>("ECAL_tracksQPtAtECAL", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );


  hECAL_tracksAtHCAL = fs->make<TProfile2D>("ECAL_tracksAtHCAL", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

  hECAL_tracksPtAtHCAL = fs->make<TProfile2D>("ECAL_tracksPtAtHCAL", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

  hECAL_tracksQPtAtHCAL = fs->make<TProfile2D>("ECAL_tracksQPtAtHCAL", "qxPt(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

} // branchesTracksAtECALstitched()

// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
//void fillTracksAtECAL_with_EEproj ( TH2F *hEvt_EE_tracksPt_, TH2F *hEvt_EE_tracksQPt_, int ieta_global_offset, int ieta_signed_offset ) {
void fillTracksAtECAL_with_EEproj ( int side, int ieta_global_offset, int ieta_signed_offset ) {

  int ieta_global_, ieta_signed_;
  int ieta_, iphi_, idx_;
  float trackPt_,trackPtAtECAL_,trackPtAtHCAL_;
  float trackQPt_,trackQPtAtECAL_,trackQPtAtHCAL_;

  for (int ieta = 1; ieta < hEvt_EE_tracksPt[side]->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    ieta_global_ = ieta_ + ieta_global_offset;
    ieta_signed_ = ieta_ + ieta_signed_offset;
    for (int iphi = 1; iphi < hEvt_EE_tracksPt[side]->GetNbinsX()+1; iphi++) {

      trackPt_ = hEvt_EE_tracksPt[side]->GetBinContent( iphi, ieta );
      trackQPt_ = hEvt_EE_tracksQPt[side]->GetBinContent( iphi, ieta );
      trackPtAtECAL_ = hEvt_EE_tracksPtAtECAL[side]->GetBinContent( iphi, ieta );
      trackQPtAtECAL_ = hEvt_EE_tracksQPtAtECAL[side]->GetBinContent( iphi, ieta );
      trackPtAtHCAL_ = hEvt_EE_tracksPtAtHCAL[side]->GetBinContent( iphi, ieta );
      trackQPtAtHCAL_ = hEvt_EE_tracksQPtAtHCAL[side]->GetBinContent( iphi, ieta );

      //if ( (trackPt_ == 0.) ) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 5*38; // shift
      iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_global_*EB_IPHI_MAX + iphi_;

      if ( (trackPt_ != 0.) ){ 
        // Fill vector for image
        vECAL_tracksPt_[idx_] = trackPt_;
        vECAL_tracksQPt_[idx_] = trackQPt_;
      
        // Fill histogram for monitoring
        hECAL_tracks->Fill( iphi_, ieta_signed_, 1. );
        hECAL_tracksPt->Fill( iphi_, ieta_signed_, trackPt_ );
        hECAL_tracksQPt->Fill( iphi_, ieta_signed_, trackQPt_ );
      }

      if ( (trackPtAtECAL_ != 0.) ){
        vECAL_tracksPtAtECAL_[idx_] = trackPtAtECAL_;
        vECAL_tracksQPtAtECAL_[idx_] = trackQPtAtECAL_;
        hECAL_tracksAtECAL->Fill( iphi_, ieta_signed_, 1. );
        hECAL_tracksPtAtECAL->Fill( iphi_, ieta_signed_, trackPtAtECAL_ );
        hECAL_tracksQPtAtECAL->Fill( iphi_, ieta_signed_, trackQPtAtECAL_ );
      }

      if ( (trackPtAtHCAL_ != 0.) ){
        vECAL_tracksPtAtHCAL_[idx_] = trackPtAtHCAL_;
        vECAL_tracksQPtAtHCAL_[idx_] = trackQPtAtHCAL_;
        hECAL_tracksAtHCAL->Fill( iphi_, ieta_signed_, 1. );
        hECAL_tracksPtAtHCAL->Fill( iphi_, ieta_signed_, trackPtAtHCAL_ );
        hECAL_tracksQPtAtHCAL->Fill( iphi_, ieta_signed_, trackQPtAtHCAL_ );
      }

    } // iphi_
  } // ieta_

} // fillTracksAtECAL_with_EEproj

// Fill stitched EE-, EB, EE+ rechits ________________________________________________________//
void RecHitAnalyzer::fillTracksAtECALstitched ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int iphi_, ieta_, iz_, idx_;
  int ieta_global, ieta_signed;
  int ieta_global_offset, ieta_signed_offset;
  float eta, phi, etaAtECAL, phiAtECAL=0., etaAtHCAL, phiAtHCAL=0., trackPt_, trackQPt_;
  GlobalPoint pos;

  edm::ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);

  vECAL_tracksPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksQPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_tracksPt[iz]->Reset();
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_tracksQPt[iz]->Reset();

  vECAL_tracksPtAtECAL_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksQPtAtECAL_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_tracksPtAtECAL[iz]->Reset();
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_tracksQPtAtECAL[iz]->Reset();

  vECAL_tracksPtAtHCAL_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksQPtAtHCAL_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_tracksPtAtHCAL[iz]->Reset();
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_tracksQPtAtHCAL[iz]->Reset();

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByToken( EBRecHitCollectionT_, EBRecHitsH_ );
  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByToken( EERecHitCollectionT_, EERecHitsH_ );
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  edm::Handle<reco::TrackCollection> tracksH_;
  iEvent.getByToken( trackCollectionT_, tracksH_ );

  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) {
    if ( !(iTk->quality(tkQt_)) ) continue;
    eta = iTk->eta();
    phi = iTk->phi();
    //if ( std::abs(eta) > 3. ) continue;
    auto propagatedECALTrack = spr::propagateTrackToECAL(&*iTk, magfield.product());
    auto propagatedHCALTrack = spr::propagateTrackToHCAL(&*iTk, magfield.product());
    if (propagatedECALTrack.ok)
    {
      etaAtECAL = propagatedECALTrack.direction.eta();
      phiAtECAL = propagatedECALTrack.direction.phi();
    }
    if (propagatedHCALTrack.ok)
    {
      etaAtHCAL = propagatedHCALTrack.direction.eta();
      phiAtHCAL = propagatedHCALTrack.direction.phi();
    }
    if ( std::abs(eta) <= 3. ){
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    //if ( id.subdetId() == EcalBarrel ) continue;
    if ( id.subdetId() == EcalEndcap ) {
      iz_ = (eta > 0.) ? 1 : 0;
      // Fill intermediate helper histogram by eta,phi
      hEvt_EE_tracksPt[iz_]->Fill( phi, eta, iTk->pt() );
      hEvt_EE_tracksQPt[iz_]->Fill( phi, eta, iTk->charge()*iTk->pt() );
    }}
    if ( std::abs(etaAtECAL) <= 3. && propagatedECALTrack.ok){
    DetId idAtECAL( spr::findDetIdECAL( caloGeom, etaAtECAL, phiAtECAL, false ) );
    if ( idAtECAL.subdetId() == EcalEndcap ) {
      iz_ = (etaAtECAL > 0.) ? 1 : 0;
      hEvt_EE_tracksPtAtECAL[iz_]->Fill( phiAtECAL, etaAtECAL, iTk->pt() );
      hEvt_EE_tracksQPtAtECAL[iz_]->Fill( phiAtECAL, etaAtECAL, iTk->charge()*iTk->pt() );
    }}
    if ( std::abs(etaAtHCAL) <= 3. && propagatedHCALTrack.ok){
    DetId idAtHCAL( spr::findDetIdECAL( caloGeom, etaAtHCAL, phiAtHCAL, false ) );
    if ( idAtHCAL.subdetId() == EcalEndcap ) {
      iz_ = (etaAtHCAL > 0.) ? 1 : 0;
      hEvt_EE_tracksPtAtHCAL[iz_]->Fill( phiAtHCAL, etaAtHCAL, iTk->pt() );
      hEvt_EE_tracksQPtAtHCAL[iz_]->Fill( phiAtHCAL, etaAtHCAL, iTk->charge()*iTk->pt() );
    }}
  } // tracks


  // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
  ieta_global_offset = 0;
  ieta_signed_offset = -ECAL_IETA_MAX_EXT;
  fillTracksAtECAL_with_EEproj( 0, ieta_global_offset, ieta_signed_offset );

  // Fill middle part of ECAL(iphi,ieta) with the EB rechits.
  ieta_global_offset = 55;

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) { 
    if ( !(iTk->quality(tkQt_)) ) continue;
    eta = iTk->eta();
    phi = iTk->phi();
    trackPt_ = iTk->pt();
    trackQPt_ = (iTk->charge()*iTk->pt());
    auto propagatedECALTrack = spr::propagateTrackToECAL(&*iTk, magfield.product());
    auto propagatedHCALTrack = spr::propagateTrackToHCAL(&*iTk, magfield.product());
    if (propagatedECALTrack.ok)
    {
      etaAtECAL = propagatedECALTrack.direction.eta();
      phiAtECAL = propagatedECALTrack.direction.phi();
    }
    if (propagatedHCALTrack.ok)
    {
      etaAtHCAL = propagatedHCALTrack.direction.eta();
      phiAtHCAL = propagatedHCALTrack.direction.phi();
    }
    //if ( std::abs(eta) > 3. ) continue;
    if ( std::abs(eta) <= 3. ){
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    //if ( id.subdetId() == EcalEndcap ) continue;
    if ( id.subdetId() == EcalBarrel ) { 
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      //if ( trackPt_ == 0. ) continue;
      // Fill vector for image
      ieta_signed = ieta_;
      ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
      idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
      vECAL_tracksPt_[idx_] += trackPt_;
      vECAL_tracksQPt_[idx_] += trackQPt_;
      // Fill histogram for monitoring
      hECAL_tracks->Fill( iphi_, ieta_signed, 1. );
      hECAL_tracksPt->Fill( iphi_, ieta_signed, trackPt_ );
      hECAL_tracksQPt->Fill( iphi_, ieta_signed, trackQPt_ );
    }}

    if ( std::abs(etaAtECAL) <= 3. && propagatedECALTrack.ok){
    DetId idAtECAL( spr::findDetIdECAL( caloGeom, etaAtECAL, phiAtECAL, false ) );
    if ( idAtECAL.subdetId() == EcalBarrel ) { 
      EBDetId ebId( idAtECAL );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      ieta_signed = ieta_;
      ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
      idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
      vECAL_tracksPtAtECAL_[idx_] += trackPt_;
      vECAL_tracksQPtAtECAL_[idx_] += trackQPt_;
      hECAL_tracksAtECAL->Fill( iphi_, ieta_signed, 1. );
      hECAL_tracksPtAtECAL->Fill( iphi_, ieta_signed, trackPt_ );
      hECAL_tracksQPtAtECAL->Fill( iphi_, ieta_signed, trackQPt_ );
    }}

    if ( std::abs(etaAtHCAL) <= 3. && propagatedHCALTrack.ok){
    DetId idAtHCAL( spr::findDetIdECAL( caloGeom, etaAtHCAL, phiAtHCAL, false ) );
    if ( idAtHCAL.subdetId() == EcalBarrel ) { 
      EBDetId ebId( idAtHCAL );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      ieta_signed = ieta_;
      ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
      idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
      vECAL_tracksPtAtHCAL_[idx_] += trackPt_;
      vECAL_tracksQPtAtHCAL_[idx_] += trackQPt_;
      hECAL_tracksAtHCAL->Fill( iphi_, ieta_signed, 1. );
      hECAL_tracksPtAtHCAL->Fill( iphi_, ieta_signed, trackPt_ );
      hECAL_tracksQPtAtHCAL->Fill( iphi_, ieta_signed, trackQPt_ );
    }}

  } // EB Tracks



  // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
  ieta_global_offset = ECAL_IETA_MAX_EXT + EB_IETA_MAX;
  ieta_signed_offset = EB_IETA_MAX;
  fillTracksAtECAL_with_EEproj( 1, ieta_global_offset, ieta_signed_offset );

} // fillTracksAtECALstitched()
