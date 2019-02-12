#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill Tracks into stitched EEm_EB_EEp image //////////////////////
// Store all Track positions into a stitched EEm_EB_EEp image 

TH2F *hEvt_Adj_tracksPt;
TH2F *hEvt_Adj_tracksQPt;
TH1F *hEvt_Adj;
TProfile2D *hECALadj_tracks;
TProfile2D *hECALadj_tracksPt;
TProfile2D *hECALadj_tracksQPt;

TProfile2D *hECALadj_tracks_id;
TProfile2D *hECALadj_tracksPt_id;
TProfile2D *hECALadj_tracksQPt_id;

std::vector<float> vECALadj_tracksPt_;
std::vector<float> vECALadj_tracksQPt_;



// Initialize branches _______________________________________________________________//
void RecHitAnalyzer::branchesTracksAtECALadjustable ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("ECALadj_tracksPt",    &vECALadj_tracksPt_);
  tree->Branch("ECALadj_tracksQPt",    &vECALadj_tracksQPt_);

    // int granularityMultiPhi;
    // int granularityMultiEta;
    std::vector<double> adjEtaBins;
    std::vector<double> adjPhiBins;

  // Intermediate helper histogram (single event only)
  // hEvt_EE_tracksPt[0] = new TH2F("evt_EEm_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
  //     EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
  //     5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  // hEvt_EE_tracksPt[1] = new TH2F("evt_EEp_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
  //     EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
  //     5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  // hEvt_EE_tracksQPt[0] = new TH2F("evt_EEm_tracksQPt", "qxPt(i#phi,i#eta);i#phi;i#eta",
  //     EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
  //     5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  // hEvt_EE_tracksQPt[1] = new TH2F("evt_EEp_tracksQPt", "qxPt(i#phi,i#eta);i#phi;i#eta",
  //     EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
  //     5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
// static const int HBHE_IPHI_NUM = hcaldqm::constants::IPHI_NUM;//72;
// static const int HBHE_IPHI_MIN = hcaldqm::constants::IPHI_MIN;//1;
// static const int HBHE_IPHI_MAX = hcaldqm::constants::IPHI_MAX;//72;
// static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
// static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;

  int eta_nbins_HBHE = 2*(hcaldqm::constants::IETA_MAX_HE-1);
  //int eta_bins_HBHE_size = 2*(hcaldqm::constants::IETA_MAX_HE-1)+1;
  int totalMultiEta = granularityMultiEta * granularityMultiECAL;
  // totalMultiEta * ( eta_bins_HBHE_size - 1 ) + 1
  // granularityMultiEta * granularityMultiECAL * ( eta_bins_HBHE_size - 1 ) + 1
  // granularityMultiEta * granularityMultiECAL * 2 * (hcaldqm::constants::IETA_MAX_HE - 1 ) + 1
  for (int i=0; i<eta_nbins_HBHE; i++)
  {
    //adjEtaBins.push_back(eta_bins_HBHE[i]);
    double step=(eta_bins_HBHE[i+1]-eta_bins_HBHE[i])/totalMultiEta;
    for (int j=0; j<totalMultiEta; j++)
    {
      adjEtaBins.push_back(eta_bins_HBHE[i]+step*j);
    }
  }
  adjEtaBins.push_back(eta_bins_HBHE[eta_nbins_HBHE]);

  //int totalMultiPhi = granularityMultiPhi * granularityMultiECAL;
  totalEtaBins = totalMultiEta*(eta_nbins_HBHE);
  totalPhiBins = granularityMultiPhi * granularityMultiECAL*HBHE_IPHI_NUM;
  //double phi_step=2*TMath::Pi()/(totalMultiPhi*HBHE_IPHI_NUM);

  //intermediate histograms
  hEvt_Adj_tracksPt = new TH2F("evt_Adj_tracksPt", "Pt(#phi,#eta);#phi;#eta",
       totalPhiBins, -TMath::Pi(), TMath::Pi(),
       adjEtaBins.size()-1, &adjEtaBins[0] );
  hEvt_Adj_tracksQPt = new TH2F("evt_Adj_tracksQPt", "qxPt(#phi,#eta);#phi;#eta",
       totalPhiBins, -TMath::Pi(), TMath::Pi(),
       adjEtaBins.size()-1, &adjEtaBins[0] );
  hEvt_Adj = new TH1F("evt_Adj", "",
       adjEtaBins.size()-1, &adjEtaBins[0] );

  // static const double eta_bins_HBHE[2*(hcaldqm::constants::IETA_MAX_HE-1)+1] =
  //                 {-3.000, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
  //                  -1.218, -1.131, -1.044, -0.957, -0.870, -0.783, -0.695, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000,
  //                   0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,  0.695,  0.783,  0.870,  0.957,  1.044,  1.131,  1.218,
  //                   1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,  1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  3.000}; // 57
  // hEvt_EE_energy[0] = new TH2F("evt_EEm_energy", "E(i#phi,i#eta);i#phi;i#eta",
  //     EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
  //     5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  // hEvt_EE_energy[1] = new TH2F("evt_EEp_energy", "E(i#phi,i#eta);i#phi;i#eta",
  //     EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
  //     5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  // // Histograms for monitoring
  // hECAL_energy = fs->make<TProfile2D>("ECAL_energy", "E(i#phi,i#eta);i#phi;i#eta",
  //     EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
  //     2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );




  // Histograms for monitoring
  // hECALadj_tracks = fs->make<TProfile2D>("ECALadj_tracks", "E(#phi,#eta);#phi;#eta",
  //      totalPhiBins, -TMath::Pi(), TMath::Pi(),
  //      adjEtaBins.size()-1, &adjEtaBins[0] );

  // hECALadj_tracksPt = fs->make<TProfile2D>("ECALadj_tracksPt", "E(#phi,#eta);#phi;#eta",
  //      totalPhiBins, -TMath::Pi(), TMath::Pi(),
  //      adjEtaBins.size()-1, &adjEtaBins[0] );

  // hECALadj_tracksQPt = fs->make<TProfile2D>("ECALadj_tracksQPt", "qxPt(#phi,#eta);#phi;#eta",
  //      totalPhiBins, -TMath::Pi(), TMath::Pi(),
  //      adjEtaBins.size()-1, &adjEtaBins[0] );




  // hECALadj_tracks_id = fs->make<TProfile2D>("ECALadj_tracks_id", "E(i#phi,i#eta);i#phi;i#eta",
  //      totalPhiBins, 0, totalPhiBins,
  //      totalEtaBins, 1, totalEtaBins+1 );

  hECALadj_tracksPt_id = fs->make<TProfile2D>("ECALadj_tracksPt_id", "E(i#phi,i#eta);i#phi;i#eta",
       totalPhiBins, 0, totalPhiBins,
       totalEtaBins, -140, totalEtaBins-140 );

  // hECALadj_tracksQPt_id = fs->make<TProfile2D>("ECALadj_tracksQPt_id", "qxPt(i#phi,i#eta);i#phi;i#eta",
  //      totalPhiBins, 0, totalPhiBins,
  //      totalEtaBins, -140, totalEtaBins-140 );

} // branchesTracksAtECALadjustable()

// // Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
// void fillTracksAtECAL_with_EEproj ( TH2F *hEvt_EE_tracksPt_, TH2F *hEvt_EE_tracksQPt_, int ieta_global_offset, int ieta_signed_offset ) {

//   int ieta_global_, ieta_signed_;
//   int ieta_, iphi_, idx_;
//   float trackPt_;
//   float trackQPt_;

//   for (int ieta = 1; ieta < hEvt_EE_tracksPt_->GetNbinsY()+1; ieta++) {
//     ieta_ = ieta - 1;
//     ieta_global_ = ieta_ + ieta_global_offset;
//     ieta_signed_ = ieta_ + ieta_signed_offset;
//     for (int iphi = 1; iphi < hEvt_EE_tracksPt_->GetNbinsX()+1; iphi++) {

//       trackPt_ = hEvt_EE_tracksPt_->GetBinContent( iphi, ieta );
//       trackQPt_ = hEvt_EE_tracksQPt_->GetBinContent( iphi, ieta );
//       if ( (trackPt_ == 0.) ) continue;
//       // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
//       iphi_ = iphi  + 5*38; // shift
//       iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
//       iphi_ = iphi_ - 1;
//       idx_  = ieta_global_*EB_IPHI_MAX + iphi_;
//       // Fill vector for image
//       vECAL_tracksPt_[idx_] = trackPt_;
//       vECAL_tracksQPt_[idx_] = trackQPt_;
//       // Fill histogram for monitoring
//       hECAL_tracks->Fill( iphi_, ieta_signed_, 1. );
//       hECAL_tracksPt->Fill( iphi_, ieta_signed_, trackPt_ );
//       hECAL_tracksQPt->Fill( iphi_, ieta_signed_, trackQPt_ );

//     } // iphi_
//   } // ieta_

// } // fillTracksAtECAL_with_EEproj

// Fill stitched EE-, EB, EE+ rechits ________________________________________________________//
void RecHitAnalyzer::fillTracksAtECALadjustable ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  vECALadj_tracksPt_.assign( totalEtaBins*totalPhiBins, 0. );
  vECALadj_tracksQPt_.assign( totalEtaBins*totalPhiBins, 0. );
  hEvt_Adj_tracksPt->Reset();
  hEvt_Adj_tracksQPt->Reset();

  // int iphi_, ieta_, iz_, idx_;
  // int ieta_global, ieta_signed;
  // int ieta_global_offset, ieta_signed_offset;
  float eta=0, phi=0, trackPt_=0, trackQPt_=0;
  // GlobalPoint pos;

  // vECAL_tracksPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  // vECAL_tracksQPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  // for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_tracksPt[iz]->Reset();
  // for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_tracksQPt[iz]->Reset();

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

    eta =       iTk->eta();
    phi =       iTk->phi();
    trackPt_ =  iTk->pt();
    trackQPt_ = iTk->charge()*iTk->pt();

    if ( std::abs(eta) <= 3. ){
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel )
    { 
      auto subDetGeometry = caloGeom->getSubdetectorGeometry(id);
      auto caloCellGeometry = subDetGeometry->getGeometry(id);
      auto corners = caloCellGeometry->getCorners();
      auto reference= caloCellGeometry->getPosition();
      auto reference_phi = reference.phi();
      auto reference_eta = reference.eta();
      float kappa= 4*TMath::Pi()/HBHE_IPHI_NUM;
      if (reference_phi>-kappa)
        reference_phi=reference_phi+kappa-TMath::Pi();
      else  
        reference_phi=reference_phi+kappa+TMath::Pi();

      EBDetId ebId( id );
      hEvt_Adj_tracksPt->SetBinContent(  ebId.iphi() - 1 +1, (ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta()) +141,
        hEvt_Adj_tracksPt->GetBinContent(ebId.iphi() - 1 +1, (ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta()) +141)+trackPt_ );
      hEvt_Adj_tracksQPt->Fill( reference_phi, reference_eta, trackQPt_ );


      Int_t bin_number;
      bin_number=hEvt_Adj->Fill( reference_eta, trackPt_ );
      
      
      
    if(ebId.ieta()>30 && ebId.ieta()<70)
      std::cout<<" eta:"<<eta<<
                 " reference_eta:"<<reference_eta<<
                 " bin_number:"<<bin_number<<
                 " predicted ieta:"<<bin_number-141<<
                 " real ieta:"<<ebId.ieta()<<
                 " bin width:"<<hEvt_Adj->GetBinWidth(bin_number)<<
                 " bin center:"<<hEvt_Adj->GetBinCenter(bin_number)<<
                 " bin low edge:"<<hEvt_Adj->GetBinLowEdge(bin_number)<<std::endl;
      ///Evt_Adj_tracksPt->GetBinXYZ


    }
    else
    {
      //mess with phi in order to match original index arithmetic. Ugh.
      float kappa= 4*TMath::Pi()/HBHE_IPHI_NUM;
      if (phi>-kappa)
        phi=phi+kappa-TMath::Pi();
      else  
        phi=phi+kappa+TMath::Pi();

      hEvt_Adj_tracksPt->Fill(  phi, eta, trackPt_ );
      hEvt_Adj_tracksQPt->Fill( phi, eta, trackQPt_ );

    }
    }

  } // tracks

  int index1d=0;
  for ( int ieta=1; ieta<=totalEtaBins; ieta++ )
  {
    for ( int iphi=1; iphi<=totalPhiBins; iphi++ )
    {
      index1d= (ieta-1)*totalPhiBins+iphi-1;//ieta_global*EB_IPHI_MAX + iphi_; 
      vECALadj_tracksPt_[index1d] += hEvt_Adj_tracksPt->GetBinContent(iphi,ieta);
      vECALadj_tracksQPt_[index1d] += hEvt_Adj_tracksQPt->GetBinContent(iphi,ieta);
      if(hEvt_Adj_tracksPt->GetBinContent(iphi,ieta)>0.0001)
      {
        //hECALadj_tracks_id->Fill(    iphi-1, ieta-141, 1. );
        hECALadj_tracksPt_id->Fill(  iphi-1, ieta-141, hEvt_Adj_tracksPt->GetBinContent(iphi,ieta ));
        //hECALadj_tracksQPt_id->Fill( iphi-1, ieta-141, hEvt_Adj_tracksQPt->GetBinContent(iphi,ieta) );
      }
    }
  }




  // for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
  //       iTk != tracksH_->end(); ++iTk ) {
  //   if ( !(iTk->quality(tkQt_)) ) continue;
  //   eta = iTk->eta();
  //   phi = iTk->phi();
  //   if ( std::abs(eta) > 3. ) continue;
  //   DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
  //   if ( id.subdetId() == EcalBarrel ) continue;
  //   if ( id.subdetId() == EcalEndcap ) {
  //     iz_ = (eta > 0.) ? 1 : 0;
  //     // Fill intermediate helper histogram by eta,phi
  //     hEvt_EE_tracksPt[iz_]->Fill( phi, eta, iTk->pt() );
  //     hEvt_EE_tracksQPt[iz_]->Fill( phi, eta, iTk->charge()*iTk->pt() );
  //   }
  // } // tracks


  // // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
  // ieta_global_offset = 0;
  // ieta_signed_offset = -ECAL_IETA_MAX_EXT;
  // fillTracksAtECAL_with_EEproj( hEvt_EE_tracksPt[0], hEvt_EE_tracksQPt[0], ieta_global_offset, ieta_signed_offset );

  // // Fill middle part of ECAL(iphi,ieta) with the EB rechits.
  // ieta_global_offset = 55;

  // for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
  //       iTk != tracksH_->end(); ++iTk ) { 
  //   if ( !(iTk->quality(tkQt_)) ) continue;
  //   eta = iTk->eta();
  //   phi = iTk->phi();
  //   trackPt_ = iTk->pt();
  //   trackQPt_ = (iTk->charge()*iTk->pt());
  //   if ( std::abs(eta) > 3. ) continue;
  //   DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
  //   if ( id.subdetId() == EcalEndcap ) continue;
  //   if ( id.subdetId() == EcalBarrel ) { 
  //     EBDetId ebId( id );
  //     iphi_ = ebId.iphi() - 1;
  //     ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
  //     if ( trackPt_ == 0. ) continue;
  //     // Fill vector for image
  //     ieta_signed = ieta_;
  //     ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
  //     idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
  //     vECAL_tracksPt_[idx_] += trackPt_;
  //     vECAL_tracksQPt_[idx_] += trackQPt_;
  //     // Fill histogram for monitoring
  //     hECAL_tracks->Fill( iphi_, ieta_signed, 1. );
  //     hECAL_tracksPt->Fill( iphi_, ieta_signed, trackPt_ );
  //     hECAL_tracksQPt->Fill( iphi_, ieta_signed, trackQPt_ );
  //   }

  // } // EB Tracks



  // // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
  // ieta_global_offset = ECAL_IETA_MAX_EXT + EB_IETA_MAX;
  // ieta_signed_offset = EB_IETA_MAX;
  // fillTracksAtECAL_with_EEproj( hEvt_EE_tracksPt[1], hEvt_EE_tracksQPt[1], ieta_global_offset, ieta_signed_offset );

} // fillTracksAtECALadjustable()
