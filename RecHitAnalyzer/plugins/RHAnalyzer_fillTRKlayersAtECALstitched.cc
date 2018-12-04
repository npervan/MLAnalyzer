#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill TRK rec hits ////////////////////////////////
// by layer at ECAL stitched

TH2F *hTOB_ECAL[nTOB];
std::vector<float> vTOB_ECAL_[nTOB];
TH2F *hEvt_EE_TOB[nTOB][nEE];

TH2F *hTEC_ECAL[nTEC];
std::vector<float> vTEC_ECAL_[nTEC];
TH2F *hEvt_EE_TEC[nTEC][nEE];

TH2F *hTIB_ECAL[nTIB];
std::vector<float> vTIB_ECAL_[nTIB];
TH2F *hEvt_EE_TIB[nTIB][nEE];

TH2F *hTID_ECAL[nTID];
std::vector<float> vTID_ECAL_[nTID];
TH2F *hEvt_EE_TID[nTID][nEE];

TH2F *hBPIX_ECAL[nBPIX];
std::vector<float> vBPIX_ECAL_[nBPIX];
TH2F *hEvt_EE_BPIX[nBPIX][nEE];

TH2F *hFPIX_ECAL[nFPIX];
std::vector<float> vFPIX_ECAL_[nFPIX];
TH2F *hEvt_EE_FPIX[nFPIX][nEE];

// Initialize branches ____________________________________________________________//
void RecHitAnalyzer::branchesTRKlayersAtECALstitched ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images

  int layer;
  char hname[50], htitle[50];
  const float * eta_bins_EE[2] = {eta_bins_EEm,eta_bins_EEp};

  //TOB
  for ( int iL(0); iL < nTOB; iL++ ) {
    // Branches for images
    layer = iL + 1;
    sprintf(hname, "TOB_layer%d_ECAL",layer);
    tree->Branch(hname,        &vTOB_ECAL_[iL]);

    // Histograms for monitoring
    sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
    hTOB_ECAL[iL] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

    for ( int iz(0); iz < nEE; iz++ ) {
      const char *zside = (iz > 0) ? "p" : "m";
      sprintf(hname, "evt_TOB_layer%d_EE%s",layer,zside);
      sprintf(htitle,"N(ix,iy);ix;iy");
      hEvt_EE_TOB[iL][iz] = new TH2F(hname, htitle,
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
    } // iz

  } // iL

  //TEC
  for ( int iL(0); iL < nTEC; iL++ ) {
    // Branches for images
    layer = iL + 1;
    sprintf(hname, "TEC_layer%d_ECAL",layer);
    tree->Branch(hname,        &vTEC_ECAL_[iL]);

    // Histograms for monitoring
    sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
    hTEC_ECAL[iL] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

    for ( int iz(0); iz < nEE; iz++ ) {
      const char *zside = (iz > 0) ? "p" : "m";
      sprintf(hname, "evt_TEC_layer%d_EE%s",layer,zside);
      sprintf(htitle,"N(ix,iy);ix;iy");
      hEvt_EE_TEC[iL][iz] = new TH2F(hname, htitle,
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
    } // iz

  } // iL

  //TIB
  for ( int iL(0); iL < nTIB; iL++ ) {
    // Branches for images
    layer = iL + 1;
    sprintf(hname, "TIB_layer%d_ECAL",layer);
    tree->Branch(hname,        &vTIB_ECAL_[iL]);

    // Histograms for monitoring
    sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
    hTIB_ECAL[iL] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

    for ( int iz(0); iz < nEE; iz++ ) {
      const char *zside = (iz > 0) ? "p" : "m";
      sprintf(hname, "evt_TIB_layer%d_EE%s",layer,zside);
      sprintf(htitle,"N(ix,iy);ix;iy");
      hEvt_EE_TIB[iL][iz] = new TH2F(hname, htitle,
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
    } // iz

  } // iL

  //TID
  for ( int iL(0); iL < nTID; iL++ ) {
    // Branches for images
    layer = iL + 1;
    sprintf(hname, "TID_layer%d_ECAL",layer);
    tree->Branch(hname,        &vTID_ECAL_[iL]);

    // Histograms for monitoring
    sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
    hTID_ECAL[iL] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

    for ( int iz(0); iz < nEE; iz++ ) {
      const char *zside = (iz > 0) ? "p" : "m";
      sprintf(hname, "evt_TID_layer%d_EE%s",layer,zside);
      sprintf(htitle,"N(ix,iy);ix;iy");
      hEvt_EE_TID[iL][iz] = new TH2F(hname, htitle,
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
    } // iz

  } // iL

  //BPIX
  for ( int iL(0); iL < nBPIX; iL++ ) {
    // Branches for images
    layer = iL + 1;
    sprintf(hname, "BPIX_layer%d_ECAL",layer);
    tree->Branch(hname,        &vBPIX_ECAL_[iL]);

    // Histograms for monitoring
    sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
    hBPIX_ECAL[iL] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

    for ( int iz(0); iz < nEE; iz++ ) {
      const char *zside = (iz > 0) ? "p" : "m";
      sprintf(hname, "evt_BPIX_layer%d_EE%s",layer,zside);
      sprintf(htitle,"N(ix,iy);ix;iy");
      hEvt_EE_BPIX[iL][iz] = new TH2F(hname, htitle,
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
    } // iz

  } // iL

  //FPIX
  for ( int iL(0); iL < nFPIX; iL++ ) {
    // Branches for images
    layer = iL + 1;
    sprintf(hname, "FPIX_layer%d_ECAL",layer);
    tree->Branch(hname,        &vFPIX_ECAL_[iL]);

    // Histograms for monitoring
    sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
    hFPIX_ECAL[iL] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

    for ( int iz(0); iz < nEE; iz++ ) {
      const char *zside = (iz > 0) ? "p" : "m";
      sprintf(hname, "evt_FPIX_layer%d_EE%s",layer,zside);
      sprintf(htitle,"N(ix,iy);ix;iy");
      hEvt_EE_FPIX[iL][iz] = new TH2F(hname, htitle,
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
    } // iz

  } // iL

} // branchesEB()


// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
void fillTRKatECAL_with_EEproj ( TH2F *hEvt_EE_tracksPt_, int ieta_global_offset, int ieta_signed_offset
TH2F *hSUBDET_ECAL[nTOB];
std::vector<float> vSUBDET_ECAL_[nTOB];
TH2F *hEvt_EE_SUBDET[nTOB][nEE];


  ) {

  int ieta_global_, ieta_signed_;
  int ieta_, iphi_, idx_;
  float trackPt_;
  float trackQPt_;

  for (int ieta = 1; ieta < hEvt_EE_tracksPt_->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    ieta_global_ = ieta_ + ieta_global_offset;
    ieta_signed_ = ieta_ + ieta_signed_offset;
    for (int iphi = 1; iphi < hEvt_EE_tracksPt_->GetNbinsX()+1; iphi++) {

      trackPt_ = hEvt_EE_tracksPt_->GetBinContent( iphi, ieta );
      trackQPt_ = hEvt_EE_tracksQPt_->GetBinContent( iphi, ieta );
      if ( (trackPt_ == 0.) ) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 5*38; // shift
      iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_global_*EB_IPHI_MAX + iphi_;
      // Fill vector for image
      vECAL_tracksPt_[idx_] = trackPt_;
      vECAL_tracksQPt_[idx_] = trackQPt_;
      // Fill histogram for monitoring
      hECAL_tracks->Fill( iphi_, ieta_signed_, 1. );
      hECAL_tracksPt->Fill( iphi_, ieta_signed_, trackPt_ );
      hECAL_tracksQPt->Fill( iphi_, ieta_signed_, trackQPt_ );

    } // iphi_
  } // ieta_

} // fillTracksAtECAL_with_EEproj


void fillTRKatEB ( EBDetId ebId, int iL, TH2F *hTRK_EB[], std::vector<float> vTRK_EB_[] ) {
  int iphi_, ieta_, idx_;
  iphi_ = ebId.iphi() - 1;
  ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
  // Fill histograms for monitoring
  hTRK_EB[iL]->Fill( iphi_, ieta_ );
  idx_ = ebId.hashedIndex(); // (ieta_+ECAL_IETA_MAX_EXT)*EB_IPHI_MAX + iphi_
  // Fill vectors for images
  vTRK_EB_[iL][idx_] += 1.;
}

//template <std::size_t N, std::size_t M> void fillTRKatEE ( EEDetId eeId, int iL, TH2F (*hTRK_EE)[N][M], std::vector<float> vTRK_layers_EE[][nEE] ) {
void fillTRKatEE ( EEDetId eeId, int iL, TH2F *hTRK_EE[][nEE], std::vector<float> vTRK_EE_[][nEE] ) {
  int ix_, iy_, iz_, idx_;
  ix_ = eeId.ix() - 1;
  iy_ = eeId.iy() - 1;
  iz_ = (eeId.zside() > 0) ? 1 : 0;
  // Fill histograms for monitoring
  hTRK_EE[iL][iz_]->Fill( ix_, iy_ );
  // Create hashed Index: maps from [iy][ix] -> [idx_]
  idx_ = iy_*EEDetId::IX_MAX + ix_;
  // Fill vectors for images
  vTRK_EE_[iL][iz_][idx_] += 1.;
}


unsigned int getLayer(const DetId& detid)
{

  unsigned int subid=detid.subdetId();

          switch(subid)
          {
            case 1://BPIX
            {
              PXBDetId pdetId = PXBDetId(detid);
              return pdetId.layer();
            }

            case 2://FPIX
            {
              PXFDetId pdetId = PXFDetId(detid.rawId());
              return pdetId.disk();
            }

            case 3://TIB
            {
              TIBDetId pdetId = TIBDetId(detid);
              return pdetId.layer();
            }
            break;

            case 4://TID
            {
              TIDDetId pdetId = TIDDetId(detid);
              return pdetId.wheel();
            }
            break;

            case 5://TOB
            {
              TOBDetId pdetId = TOBDetId(detid);
              return pdetId.layer();
            }
            break;

            case 6://TEC
            {
              TECDetId pdetId = TECDetId(detid);
              return pdetId.wheel();
            }
            break;
          }
          return 999;

}


// Fill TRK rechits at EB/EE ______________________________________________________________//
void RecHitAnalyzer::fillTRKlayersAtECALstitched ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  //int ix_, iy_, iz_;
  //int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  //int layer;
  float eta, phi;//, rho;
  GlobalPoint pos;

  for ( int iL(0); iL < nTOB; iL++ ) {
    vTOB_ECAL_[iL].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_TOB[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nTEC; iL++ ) {
    vTEC_ECAL_[iL].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_TEC[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nTIB; iL++ ) {
    vTIB_ECAL_[iL].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_TIB[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nTID; iL++ ) {
    vTID_ECAL_[iL].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_TID[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nBPIX; iL++ ) {
    vBPIX_ECAL_[iL].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_BPIX[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nFPIX; iL++ ) {
    vFPIX_ECAL_[iL].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_FPIX[iL][iz]->Reset();
  }

  edm::Handle<TrackingRecHitCollection> TRKRecHitsH_;
  iEvent.getByToken( TRKRecHitCollectionT_, TRKRecHitsH_ );
  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();


//sipixel
  edm::ESHandle<TrackerGeometry> geom;
  iSetup.get<TrackerDigiGeometryRecord>().get( geom );
  const TrackerGeometry& theTracker(*geom);


  edm::Handle<SiPixelRecHitCollection> recHitColl;
  iEvent.getByToken(siPixelRecHitCollectionT_, recHitColl);

  SiPixelRecHitCollection::const_iterator recHitIdIterator      = (recHitColl.product())->begin();
  SiPixelRecHitCollection::const_iterator recHitIdIteratorEnd   = (recHitColl.product())->end();

  for ( ; recHitIdIterator != recHitIdIteratorEnd; recHitIdIterator++)
  {
    SiPixelRecHitCollection::DetSet detset = *recHitIdIterator;
    DetId detId = DetId(detset.detId()); // Get the Detid object
    unsigned int subid=detId.subdetId(); //subdetector type, barrel=1, fpix=2
    unsigned int layer = getLayer(detId);
    const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detId) );

    SiPixelRecHitCollection::DetSet::const_iterator pixeliter=detset.begin();
    SiPixelRecHitCollection::DetSet::const_iterator rechitRangeIteratorEnd   = detset.end();
    for(;pixeliter!=rechitRangeIteratorEnd;++pixeliter)
    {//loop on the rechit

      if (pixeliter->isValid())
      {

        LocalPoint lp = pixeliter->localPosition();
        GlobalPoint GP = theGeomDet->surface().toGlobal(Local3DPoint(lp));
        phi = GP.phi();
        eta = GP.eta();
        if ( std::abs(eta) > 3. ) continue;
        DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
        if ( subid == PixelSubdetector::PixelBarrel )
        {
          hBPIX_layers->Fill(layer);
          if ( ecalId.subdetId() == EcalBarrel )
            fillTRKatEB( EBDetId(ecalId), layer-1, hBPIX_EB, vBPIX_EB_ );
          else
            if ( ecalId.subdetId() == EcalEndcap )
              fillTRKatEE( EEDetId(ecalId), layer-1, hBPIX_EE, vBPIX_EE_ );

        }
        else 
          if ( subid == PixelSubdetector::PixelEndcap )
          {
            hFPIX_layers->Fill(layer);
            if ( ecalId.subdetId() == EcalBarrel )
              fillTRKatEB( EBDetId(ecalId), layer-1, hFPIX_EB, vFPIX_EB_ );
            else
              if ( ecalId.subdetId() == EcalEndcap )
                fillTRKatEE( EEDetId(ecalId), layer-1, hFPIX_EE, vFPIX_EE_ );

          }
      }
    }
  }


  for (const auto & itoken: siStripRecHitCollectionT_)
  {
    edm::Handle<SiStripRecHit2DCollection> stripRecHitColl;
    iEvent.getByToken( itoken , stripRecHitColl);

    SiStripRecHit2DCollection::const_iterator stripRecHitIdIterator      = (stripRecHitColl.product())->begin();
    SiStripRecHit2DCollection::const_iterator stripRecHitIdIteratorEnd   = (stripRecHitColl.product())->end();

    for (; stripRecHitIdIterator != stripRecHitIdIteratorEnd; ++stripRecHitIdIterator)
    { 

     SiStripRecHit2DCollection::DetSet detset = *stripRecHitIdIterator;
     DetId detId = DetId(detset.detId()); // Get the Detid object
     unsigned int subid=detId.subdetId(); //subdetector type, barrel=1, fpix=2
     unsigned int layer = getLayer(detId);
     const StripGeomDetUnit* theGeomDet = dynamic_cast<const StripGeomDetUnit*>( theTracker.idToDet( detId ) );

     SiStripRecHit2DCollection::DetSet::const_iterator stripiter=detset.begin();
     SiStripRecHit2DCollection::DetSet::const_iterator stripRechitRangeIteratorEnd   = detset.end();
     for(;stripiter!=stripRechitRangeIteratorEnd;++stripiter)
      {
        if (stripiter->isValid())
        {

          LocalPoint lp = stripiter->localPosition();
          GlobalPoint GP = theGeomDet->surface().toGlobal(Local3DPoint(lp));
          phi = GP.phi();
          eta = GP.eta();
          if ( std::abs(eta) > 3. ) continue;
          DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );

          if ( subid == StripSubdetector::TOB ) {

            hTOB_layers->Fill(layer);

            if ( ecalId.subdetId() == EcalBarrel )
              fillTRKatEB( EBDetId(ecalId), layer-1, hTOB_EB, vTOB_EB_ );
            else if ( ecalId.subdetId() == EcalEndcap )
              fillTRKatEE( EEDetId(ecalId), layer-1, hTOB_EE, vTOB_EE_ );

          } else if ( subid == StripSubdetector::TEC ) {
    
            hTEC_layers->Fill(layer);
            if ( ecalId.subdetId() == EcalBarrel )
              fillTRKatEB( EBDetId(ecalId), layer-1, hTEC_EB, vTEC_EB_ );
            else if ( ecalId.subdetId() == EcalEndcap )
              fillTRKatEE( EEDetId(ecalId), layer-1, hTEC_EE, vTEC_EE_ );

          } else if ( subid == StripSubdetector::TIB ) {
    
            hTIB_layers->Fill(layer);
            if ( ecalId.subdetId() == EcalBarrel )
              fillTRKatEB( EBDetId(ecalId), layer-1, hTIB_EB, vTIB_EB_ );
            else if ( ecalId.subdetId() == EcalEndcap )
              fillTRKatEE( EEDetId(ecalId), layer-1, hTIB_EE, vTIB_EE_ );

          } else if ( subid == StripSubdetector::TID ) {
    
            hTID_layers->Fill(layer);
            if ( ecalId.subdetId() == EcalBarrel )
              fillTRKatEB( EBDetId(ecalId), layer-1, hTID_EB, vTID_EB_ );
            else if ( ecalId.subdetId() == EcalEndcap )
              fillTRKatEE( EEDetId(ecalId), layer-1, hTID_EE, vTID_EE_ );

          }
        }
      }
    }
  }

} // fillEB()
