// -*- C++ -*-
//
// Package:    SiStripRecHitsAnalyzer
// Class:      SiStripRecHitsAnalyzer
// 
/**\class SiStripRecHitsAnalyzer SiStripRecHitsAnalyzer.cc Validation/SiStripRecHitsAnalyzer/src/SiStripRecHitsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ivan Amos Cali,32 4-A08,+41227673039,
//         Created:  Sat Nov 13 21:24:20 CET 2010
// $Id$
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
//
// class declaration
//

class SiStripRecHitsAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SiStripRecHitsAnalyzer(const edm::ParameterSet&);
      ~SiStripRecHitsAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
	  std::vector<edm::InputTag> inputTags_;
	  const TrackerGeometry* trGeo_;
	  CentralityProvider * centrality_;
	  uint32_t CentralityBin_;
	  
	  uint32_t detId_;
	  uint16_t Nstrips_;
	  double rhX_;
	  double rhY_;
	  double rhZ_;
      double rhR_;
      double rhEta_;
	  double rhPhi_;
	  
	  uint32_t nClusters_;
	  uint32_t clSize_;
	  uint32_t clFistStrip_;
	  //const std::vector<uint8_t>& vclAdcs_;
	  uint32_t clTotCharge_;
	  uint16_t APVwithSatStip_;
	  
	  TH1F* clWidth_[11];
	  TH1F* clCharge_[11];
      TH2F* clWidthVsEta_[11];
	  TH2F* clWidthVsOcc_[11];
	  TH2F* clChargeVsOcc_[11];
	  
      TH2F* clChargeVsEta_[11];
	  TH2F* clChargeVsWidth_[11];
	  
	  
	  TH2F* clChargeVsCMShift_;
	  TH2F* clWidthVsCMShift_;
	  TH2F* nClusVsCentrality_;
      


};


SiStripRecHitsAnalyzer::SiStripRecHitsAnalyzer(const edm::ParameterSet& conf){


   inputTags_ = conf.getParameter<std::vector<edm::InputTag> >("RecHitCollections");
   centrality_ = 0;
      
   edm::Service<TFileService> fs;
   for( int i = 0; i<10; i++)
  {
     clWidth_[i]= fs->make<TH1F>(Form("clusterWidth%d",i),Form("clusterWidth%d",i), 250, 1, 250 );
	 clCharge_[i]= fs->make<TH1F>(Form("clusterCharge%d",i),Form("clusterCharge%d",i), 30001, -0.5, 30000.5 );
     clWidthVsEta_[i]= fs->make<TH2F>(Form("clusterWidthVsEta%d",i),Form("clusterWidthVsEta%d",i), 1000, -5, 5, 250, 1, 250  );
	 clWidthVsOcc_[i]= fs->make<TH2F>(Form("clusterWidthVsOcc%d",i),Form("clusterWidthVsOcc%d",i),  128, 0,100, 250, 1, 250 );
	 clChargeVsOcc_[i]= fs->make<TH2F>(Form("clusterChargeVsOcc%d",i),Form("clusterChargeVsOcc%d",i),  100, 0,100, 250, 1, 250 );
     clChargeVsEta_[i]= fs->make<TH2F>(Form("clusterChargeVsEta%d",i),Form("clusterChargeVsEta%d",i),  1000, -5, 5 , 30001, -0.5, 30000.5);
	 clChargeVsWidth_[i]= fs->make<TH2F>(Form("clusterChargeVsWidth%d",i),Form("clusterChargeVsWidth%d",i), 250, 1, 250 , 30001, -0.5, 30000.5 );
     
  }
  
  
   clWidth_[10]= fs->make<TH1F>("clusterWidth","clusterWidth", 250, 1, 250 );
   clCharge_[10]= fs->make<TH1F>("clusterCharge","clusterCharge", 30001, -0.5, 30000.5 );
   clWidthVsEta_[10]= fs->make<TH2F>("clusterWidthVsEta","clusterWidthVsEta", 1000, -5, 5, 250, 1, 250  );
   clWidthVsOcc_[10]= fs->make<TH2F>("clusterWidthVsOcc","clusterWidthVsOcc", 128, 0,100,  250, 1, 250 );
   clChargeVsOcc_[10]= fs->make<TH2F>("clusterChargeVsOcc","clusterChargeVsOcc",  100, 0,100, 250, 1, 250 );
   clChargeVsEta_[10]= fs->make<TH2F>("clusterChargeVsEta","clusterChargeVsEta",  1000, -5, 5, 30001, -0.5, 30000.5 );
   clChargeVsWidth_[10]= fs->make<TH2F>("clusterChargeVsWidth","clusterChargeVsWidth", 250, 1, 250, 30001, -0.5, 30000.5 );
   
   //clWidthVsCMShift_=fs->make<TH2F>("clusterWidthVsCMShift", "clusterWidthVsCMShift", 500, -200, 300, 250, 1, 250);
   //clChargeVsCMShift_=fs->make<TH2F>("clusterChargeVsCMShift", "clusterChargeVsCMShift",500, -200, 300, 300001, -0.5, 30000.5);
   
   nClusVsCentrality_= fs->make<TH2F>("nClusVsCentrality","nClusVsCentrality", 11, -0.5, 10.5, 100000, 0, 100000 );
   
}


SiStripRecHitsAnalyzer::~SiStripRecHitsAnalyzer()
{
 

}


void
SiStripRecHitsAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
   using namespace edm;
   
   if(!centrality_) centrality_ = new CentralityProvider(es);
   centrality_->newEvent(e,es); // make sure you do this first in every event
   CentralityBin_ = centrality_->getBin();
   uint32_t hCbin = CentralityBin_/4;
   
   
    
   edm::ESHandle<TrackerGeometry> tracker;
   es.get<TrackerDigiGeometryRecord>().get( tracker );    
   trGeo_ = tracker.product(); 
   
   
  for(std::vector<edm::InputTag>::const_iterator inputTag = inputTags_.begin(); inputTag != inputTags_.end(); ++inputTag) {
   edm::Handle<SiStripRecHit2DCollection> siStripRecCollection;
   e.getByLabel(*inputTag,siStripRecCollection);
   const SiStripRecHit2DCollection RecHits = *(siStripRecCollection.product());
   
   
   SiStripRecHit2DCollection::DataContainer::const_iterator itRecHits = RecHits.data().begin();
   uint32_t olddetId=(*itRecHits).geographicalId();
   uint16_t nstripwithSignal =0; 
   for(; itRecHits != RecHits.data().end();++itRecHits) {
	  nClusters_ = RecHits.data().size();
	  detId_= (*itRecHits).geographicalId();
	  
	  const StripGeomDetUnit* stripGedmDetUnit = dynamic_cast<const StripGeomDetUnit*>(trGeo_->idToDetUnit(detId_));
      Nstrips_ = stripGedmDetUnit->specificTopology().nstrips();
	  
     //Position of the cluster
     GlobalPoint globalPosition = stripGedmDetUnit->toGlobal((*itRecHits).localPosition());
     rhX_= globalPosition.x();
	 rhY_= globalPosition.y();
	 rhZ_ = globalPosition.z();
     rhR_ = sqrt(rhX_*rhX_+rhY_*rhY_);
     rhEta_ = -log(tan(atan2(rhR_,rhZ_)/2.));
	 rhPhi_= atan2(rhY_,rhX_);
          
     
	 //SiPixelCluster* cluster1 = (*IPRH_1).cluster();
	 clSize_= (*itRecHits).cluster()->amplitudes().size();
	 clFistStrip_ = (*itRecHits).cluster()->firstStrip();
	 const std::vector<uint8_t>& vclAdcs = (*itRecHits).cluster()->amplitudes(); 
	 
	 
	 uint16_t strip = clFistStrip_;
	 clTotCharge_=0;
	 for(size_t i =0; i < vclAdcs.size(); ++i){
		clTotCharge_ += vclAdcs[i];
	    //if(vclAdcs[i] ==255) APVwithSatStip_ = 
		++strip;
	 }
	 
	 
	 //if(detId_!=olddetId){
	  // float occ = (float)nstripwithSignal / (float)Nstrips_;
	  // clWidthVsOcc_[hCbin]->Fill(occ, clSize_);
	  // clWidthVsOcc_[10]->Fill(occ, clSize_);
	  // clChargeVsOcc_[hCbin]->Fill(occ, clTotCharge_);
	  // clChargeVsOcc_[hCbin]->Fill(occ, clTotCharge_);
       //olddetId= detId_;
	   //nstripwithSignal =0;
	 //}  
	 // nstripwithSignal += clSize_;
	 
     clWidth_[hCbin]->Fill(clSize_);
	 clWidth_[10]->Fill(clSize_);
	 clCharge_[hCbin]->Fill(clTotCharge_);
	 clCharge_[10]->Fill(clTotCharge_);
     clWidthVsEta_[hCbin]->Fill(rhEta_, clSize_);
	 clWidthVsEta_[10]->Fill(rhEta_, clSize_);
	 clChargeVsEta_[hCbin]->Fill(rhEta_, clTotCharge_);
	 clChargeVsEta_[10]->Fill(rhEta_, clTotCharge_);
	 
	 clChargeVsWidth_[hCbin]->Fill(clSize_, clTotCharge_);
	 clChargeVsWidth_[10]->Fill(clSize_, clTotCharge_);
     nClusVsCentrality_->Fill(hCbin, nClusters_); 		
	
	// clWidthVsCMShift_->Fill(CMShift_, clSize_);	  
	// clChargeVsCMShift_->Fill(CMShift_, clTotCharge_);	
	
	
  }
 
 }	
}



void 
SiStripRecHitsAnalyzer::beginJob()
{
}


void 
SiStripRecHitsAnalyzer::endJob() {
}


DEFINE_FWK_MODULE(SiStripRecHitsAnalyzer);
