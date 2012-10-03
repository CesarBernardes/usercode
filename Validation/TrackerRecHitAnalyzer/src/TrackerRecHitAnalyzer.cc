// -*- C++ -*-
//
// Package:    TrackerRecHitAnalyzer
// Class:      TrackerRecHitAnalyzer
// 
/**\class TrackerRecHitAnalyzer TrackerRecHitAnalyzer.cc Validation/TrackerRecHitAnalyzer/src/TrackerRecHitAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ivan Amos Cali,32 4-A08,+41227673039,
//         Created:  Tue Oct  2 17:13:49 CEST 2012
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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"


#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"



#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"

#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//ROOT inclusion
#include "TFile.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

//
// class declaration
//

class TrackerRecHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrackerRecHitAnalyzer(const edm::ParameterSet&);
      ~TrackerRecHitAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      void GetGlobalPosition(GlobalPoint const &);
      void BarrelFillHisto(int const &);
      void ForwardFillHisto(int const & , int const &);


      // ----------member data ---------------------------
      std::vector<edm::InputTag> StripInputTags_;
      edm::InputTag PixelInputTag_; 
      std::string outputFile_;
      edm::InputTag src_;
      const TrackerGeometry* trGeo_;
      TFile* oFile_;
      uint32_t detId_;
      double recHitX_;
      double recHitY_;
      double recHitZ_;
      double recHitR_;
      double recHitEta_; 
      double recHitPhi_;
      
      TH2F* TrackerBarrelEtaPhi_[13];
      TH2F* TrackerBarrelZPhi_[13];
      TH2F* TrackerForwardYXPlus_[14];
      TH2F* TrackerForwardYXMinus_[14];
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackerRecHitAnalyzer::TrackerRecHitAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    StripInputTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("StripRecHitCollections");
    PixelInputTag_ = iConfig.getParameter<edm::InputTag >("PixelRecHitCollections");
    double pi=3.14159265359;	
    edm::Service<TFileService> fs;
    for( int i = 0, j=1; i<13; ++i, ++j) TrackerBarrelEtaPhi_[i]= fs->make<TH2F>(Form("EtaPhiL%d",j),Form("EtaPhiL%d",j), 500, -2.5, 2.5, 2*1000*pi, -pi, pi);
    for( int i = 0, j=1; i<13; ++i, ++j) TrackerBarrelZPhi_[i]= fs->make<TH2F>(Form("ZPhiL%d",j),Form("ZPhiL%d",j), 20200, -101., 101., 2*1000*pi, -pi, pi);
    for( int i = 0, j=1; i<14; ++i, ++j) TrackerForwardYXPlus_[i]= fs->make<TH2F>(Form("XYPlusD%d",j),Form("XYPlusD%d",j), 2020, -101., 101., 2020, -101., 101.);
    for( int i = 0, j=1; i<14; ++i, ++j) TrackerForwardYXMinus_[i]= fs->make<TH2F>(Form("XYMinusD%d",j),Form("XYMinusD%d",j), 2020, -101., 101., 2020, -101., 101.);
}


TrackerRecHitAnalyzer::~TrackerRecHitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackerRecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



 edm::Handle<SiPixelRecHitCollection> siPixelRecCollection;
   iEvent.getByLabel(PixelInputTag_,siPixelRecCollection);
   const SiPixelRecHitCollection pixelRecHits = *(siPixelRecCollection.product());
    
   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get( tracker );    
   trGeo_ = tracker.product(); 
   
    std::cout << "Starting Pixels" << std::endl;
   //looping over first rechit collection  
   for(SiPixelRecHitCollection::DataContainer::const_iterator IPRH = pixelRecHits.data().begin(); IPRH != pixelRecHits.data().end(); ++IPRH) {
          detId_=(*IPRH).geographicalId();
	  const PixelGeomDetUnit* pgdu = dynamic_cast<const PixelGeomDetUnit*>(trGeo_->idToDetUnit(detId_));             
	  
	  //Position of the cluster
	  GlobalPoint globalPosition = pgdu->toGlobal((*IPRH).localPosition());
          this->GetGlobalPosition(globalPosition);
	  
     
	  //barrel
          if(pgdu->subDetector() == GeomDetEnumerators::PixelBarrel){
	    //layer
	    PXBDetId detId(detId_);
	    BarrelFillHisto(detId.layer() -1);
	  } else if(pgdu->subDetector() == GeomDetEnumerators::PixelEndcap){
            PXFDetId detId(detId_); 
            ForwardFillHisto(detId.side(), detId.disk()-1);         
          }  
	} 

   //----------------------------------------------------

 for(std::vector<edm::InputTag>::const_iterator inputTag = StripInputTags_.begin(); inputTag != StripInputTags_.end(); ++inputTag) {
   edm::Handle<SiStripRecHit2DCollection> siStripRecCollection;
   iEvent.getByLabel(*inputTag,siStripRecCollection);
   const SiStripRecHit2DCollection RecHits = *(siStripRecCollection.product());
   
   
   SiStripRecHit2DCollection::DataContainer::const_iterator itRecHits = RecHits.data().begin();
   for(; itRecHits != RecHits.data().end();++itRecHits) {
	  detId_= (*itRecHits).geographicalId();
	  
	  const StripGeomDetUnit* stripGedmDetUnit = dynamic_cast<const StripGeomDetUnit*>(trGeo_->idToDetUnit(detId_));
	  
	  //Position of the cluster
	  GlobalPoint globalPosition = stripGedmDetUnit->toGlobal((*itRecHits).localPosition());
	  this->GetGlobalPosition(globalPosition);
	  if(stripGedmDetUnit->subDetector() == GeomDetEnumerators::TIB){
            TIBDetId detId(detId_);
            BarrelFillHisto(detId.layer() + 2);	  
	  } else if(stripGedmDetUnit->subDetector() ==  GeomDetEnumerators::TOB ){ 								
	    TOBDetId detId(detId_);
            BarrelFillHisto(detId.layer() + 6); 
          } else if(stripGedmDetUnit->subDetector() ==  GeomDetEnumerators::TID ){ 								
	    TIDDetId detId(detId_);
            ForwardFillHisto(detId.side(), detId.wheel()+2);
	  } else if(stripGedmDetUnit->subDetector() ==  GeomDetEnumerators::TEC ){ 								
	    TECDetId detId(detId_);
            ForwardFillHisto(detId.side(), detId.wheel()+4);
	  }
   }

 }
  
}

void TrackerRecHitAnalyzer::BarrelFillHisto(int const &layer){
  TrackerBarrelEtaPhi_[layer]->Fill(recHitEta_, recHitPhi_);
  TrackerBarrelZPhi_[layer]->Fill(recHitZ_, recHitPhi_);
}

void TrackerRecHitAnalyzer::ForwardFillHisto(int  const & side, int  const & wheel){
   if(side==2) TrackerForwardYXPlus_[wheel]->Fill(recHitX_, recHitY_);
  else TrackerForwardYXMinus_[wheel]->Fill(recHitX_, recHitY_);
}
void TrackerRecHitAnalyzer::GetGlobalPosition(GlobalPoint const & globalPosition){
          recHitX_= globalPosition.x();
	  recHitY_= globalPosition.y();
	  recHitZ_ = globalPosition.z();
	  recHitR_ = sqrt(recHitX_*recHitX_+recHitY_*recHitY_);
	  recHitEta_ = -log(tan(atan2(recHitR_,recHitZ_)/2.));
	  recHitPhi_= atan2(recHitY_,recHitX_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrackerRecHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackerRecHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TrackerRecHitAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TrackerRecHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TrackerRecHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TrackerRecHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackerRecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackerRecHitAnalyzer);

/*
unsigned int layer, stringorrod, layerg, module, side, external,TOB, moduleg;
     //TIB-------------------------------------------
     if(detId.subdetId()==StripSubdetector::TIB){  //SubDetector { TIB=3,TID=4,TOB=5,TEC=6 }
       TIBDetId tibid(id);
        layer = tibid.layer();            //1,2,3,4
        layerg = layer;                   
        module = tibid.module();          //1,2,3
        side = tibid.string()[0] -1;      //string()[0] = 1 for TIB- and = 2 for TIB+
        external = tibid.string()[1] -1;  //string()[1] = 1 for internal and = 2 for external 
        stringorrod = tibid.string()[2];
        TOB = 0;
        moduleg = (module-1)*2 + side * 6 + external + 1; //1-12: even external and odd internal  
       
     }      
        //TOB---------------------------------------------
     else if(detId.subdetId()==StripSubdetector::TOB){
       TOBDetId tobid(id);
        layer = tobid.layer();            //1,2,3,4
        layerg = layer +4;
        module = tobid.module();          //1,2,3
        side = tobid.rod()[0] -1;         //rod()[0] = 1 for TOB- and = 2 for TOB+
        external = 0;                     //all the module are internal
        stringorrod = tobid.rod()[1];
        TOB =1;
        moduleg = module + side * 6 ; //1-12 
     }
 */ 
