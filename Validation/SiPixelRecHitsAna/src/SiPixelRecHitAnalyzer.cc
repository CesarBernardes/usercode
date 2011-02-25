// -*- C++ -*-
//
// Package:    SiPixelRecHitsAna
// Class:      SiPixelRecHitsAna
// 
/**\class SiPixelRecHitsAna SiPixelRecHitsAna.cc Validation/SiPixelRecHitsAna/src/SiPixelRecHitsAna.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ivan Amos CALI
//         Created:  Fri Dec 11 16:22:55 CET 2009
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
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

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
//ROOT inclusion
#include "TFile.h"
#include "TTree.h"


//
// class decleration
//

class SiPixelRecHitsAna : public edm::EDAnalyzer {
   public:
      explicit SiPixelRecHitsAna(const edm::ParameterSet&);
      ~SiPixelRecHitsAna();
       

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      // ----------member data ---------------------------
	  std::string outputFile_;
      edm::InputTag src_;
      const TrackerGeometry* trGeo_;
      TFile* oFile_;
      TTree* Clusters_;
      TTree* ClustersDist_;
	  unsigned int nClusters_;
      unsigned int EventN_;
	  unsigned int SubEventN_;
	  int nClustersL1_;
	  int nClustersL2_;
	  int nClustersL3_;
	  
	  unsigned int ClusterNumber1_;
	  double recHitX1_;
      double recHitY1_;
	  double recHitZ1_;
	  double recHitR1_;
	  double recHitEta1_;
	  double recHitPhi1_;
      int sizeX1_;
      int sizeY1_;
	  int clSize1_;
	  int pixelsX1_[500];
	  int pixelsY1_[500];
	  int pixelsADC1_[500];
	  float charge1_;
      int minPixelRow1_;
      int minPixelCol1_;
      int maxPixelRow1_;
      int maxPixelCol1_;
      
      
      unsigned int ClusterNumber2_;
      double recHitX2_;
      double recHitY2_;
	  double recHitZ2_;
	  double recHitR2_;
	  double recHitEta2_;
	  double recHitPhi2_;
      int sizeX2_;
      int sizeY2_;
	  int clSize2_;
      int pixelsX2_[500];
	  int pixelsY2_[500];
	  int pixelsADC2_[500];
	  float charge2_;
      int minPixelRow2_;
      int minPixelCol2_;
      int maxPixelRow2_;
      int maxPixelCol2_;
      
      float dX_;
      float dY_;
	  float dR_;
      double minDist_;
  
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
SiPixelRecHitsAna::SiPixelRecHitsAna(const edm::ParameterSet& iConfig)

{
    outputFile_ = iConfig.getUntrackedParameter<std::string>("outputFile", "clusters.root");
    src_ =  iConfig.getParameter<edm::InputTag>( "src" );


}


SiPixelRecHitsAna::~SiPixelRecHitsAna()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiPixelRecHitsAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   // get trigger decisions
  bool accept = false;

  bool is_bit36 = false;
  bool is_bit37 = false;
  bool is_bit38 = false;
  bool is_bit39 = false;
  bool is_bit40 = false;
  bool is_bit41 = false;

  edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  iEvent.getByLabel( edm::InputTag("gtDigis"), gtReadoutRecord);
  const TechnicalTriggerWord&  technicalTriggerWord = gtReadoutRecord->technicalTriggerWord();
  is_bit36 = technicalTriggerWord.at(36);
  is_bit37 = technicalTriggerWord.at(37);
  is_bit38 = technicalTriggerWord.at(38);
  is_bit39 = technicalTriggerWord.at(39);
  is_bit40 = technicalTriggerWord.at(40);
  is_bit41 = technicalTriggerWord.at(41);


  bool is_vtx = false;
  edm::Handle<reco::VertexCollection> vertexCollection;
  iEvent.getByLabel("pixel3Vertices",vertexCollection);
  const reco::VertexCollection * vertices = vertexCollection.product();
  if(vertices->size() > 0)
    is_vtx = true;

  float rver_z = 0.;
  unsigned int vertex_tracks = 0;
  for(reco::VertexCollection::const_iterator
      vertex = vertices->begin(); vertex!= vertices->end(); vertex++)
  {
     if(vertex->tracksSize()>vertex_tracks) {
       vertex_tracks = vertex->tracksSize();
       rver_z = vertex->z();
     }
  }

  // define trigger
 if(!is_bit36 && !is_bit37 && !is_bit38 && !is_bit39 && is_bit40 && is_vtx) accept = true;
 //if(is_bit40 && is_vtx) accept = true;
  
if(accept){   //starting the process 
   EventN_ = iEvent.id().event();		  
   ++SubEventN_;
   
     
   // Get recHits
   edm::Handle<SiPixelRecHitCollection> siPixelRecCollection;
   iEvent.getByLabel("siPixelRecHits",siPixelRecCollection);
   const SiPixelRecHitCollection pixelRecHits = *(siPixelRecCollection.product());
    
   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get( tracker );    
   trGeo_ = tracker.product(); 
   
   //looping over first rechit collection  
   ClusterNumber1_=0; 
  
   for(SiPixelRecHitCollection::DataContainer::const_iterator IPRH_1 = pixelRecHits.data().begin(); IPRH_1 != pixelRecHits.data().end(); ++IPRH_1) {
	  unsigned int nClusters = pixelRecHits.data().size();
	  if(nClusters > 300) continue;
	  
	  const PixelGeomDetUnit* pgdu_1 = dynamic_cast<const PixelGeomDetUnit*>(trGeo_->idToDetUnit((*IPRH_1).geographicalId()));
      const RectangularPixelTopology* topol_1 = dynamic_cast<const RectangularPixelTopology*>(&(pgdu_1->specificTopology()));

     //barrel
     if(pgdu_1->subDetector() != GeomDetEnumerators::PixelBarrel)
       continue;
      
	  ++ClusterNumber1_; //need to be after selction on barrel
	  
     //Position of the cluster
     GlobalPoint globalPosition1 = pgdu_1->toGlobal((*IPRH_1).localPosition());
     recHitX1_= globalPosition1.x();
	 recHitY1_= globalPosition1.y();
	 recHitZ1_ = globalPosition1.z()-rver_z;
     recHitR1_ = sqrt(recHitX1_*recHitX1_+recHitY1_*recHitY1_);
     recHitEta1_ = -log(tan(atan2(recHitR1_,recHitZ1_)/2.));
	 recHitPhi1_= atan2(recHitY1_,recHitX1_);
     
     //layer
     PXBDetId pdetId_1 = PXBDetId(IPRH_1->geographicalId());
     int layer_1=pdetId_1.layer();
	
	 //SiPixelCluster* cluster1 = (*IPRH_1).cluster();
	 sizeX1_= (*IPRH_1).cluster()->sizeX();
	 sizeY1_= (*IPRH_1).cluster()->sizeY();
	 std::vector<SiPixelCluster::Pixel> pixels1 = (*IPRH_1).cluster()->pixels();
	 clSize1_ = pixels1.size();
	 if(clSize1_ > 500) clSize1_=500;
	 for(int i =0; i < clSize1_; ++i){
		pixelsX1_[i] = pixels1[i].x;
		pixelsY1_[i] = pixels1[i].y;
	    pixelsADC1_[i] =  pixels1[i].adc;
	 }
	charge1_ = (*IPRH_1).cluster()->charge();
	minPixelRow1_= (*IPRH_1).cluster()->minPixelRow();
	minPixelCol1_= (*IPRH_1).cluster()->minPixelCol();
	maxPixelRow1_= (*IPRH_1).cluster()->maxPixelRow();
	maxPixelCol1_= (*IPRH_1).cluster()->maxPixelCol();
	Clusters_->Fill();
			
	ClusterNumber2_=0;
	//Loop over again the same cluster collection o make the matching
	for(SiPixelRecHitCollection::DataContainer::const_iterator IPRH_2 = pixelRecHits.data().begin(); IPRH_2 != pixelRecHits.data().end();++IPRH_2) {
        
		
		const PixelGeomDetUnit* pgdu_2 = dynamic_cast<const PixelGeomDetUnit*>(trGeo_->idToDetUnit((*IPRH_2).geographicalId()));
        const RectangularPixelTopology* topol_2 = dynamic_cast<const RectangularPixelTopology*>(&(pgdu_2->specificTopology()));

        //barrel
        if(pgdu_2->subDetector() != GeomDetEnumerators::PixelBarrel) continue;

        // no double counting (should be after "barrel")
        ++ClusterNumber2_; 
	    if(ClusterNumber2_ <= ClusterNumber1_) continue; 
        
		//layer
        PXBDetId pdetId_2 = PXBDetId(IPRH_2->geographicalId());
        int layer_2=pdetId_2.layer();
        
		if(layer_2 != layer_1)  continue;  //check that the layer is the same 
		
		GlobalPoint globalPosition2 = pgdu_2->toGlobal((*IPRH_1).localPosition());
        recHitX2_= globalPosition2.x();
	    recHitY2_= globalPosition2.y();
	    recHitZ2_ = globalPosition2.z()-rver_z;
        recHitR2_ = sqrt(recHitX2_*recHitX2_+recHitY2_*recHitY2_);
        recHitEta2_ = -log(tan(atan2(recHitR2_,recHitZ2_)/2.));
	    recHitPhi2_= atan2(recHitY2_,recHitX2_);
		
		//SiPixelCluster* cluster2 = (*IPRH_2).cluster();
		sizeX2_= (*IPRH_2).cluster()->sizeX();
	    sizeY2_= (*IPRH_2).cluster()->sizeY();
	    std::vector<SiPixelCluster::Pixel> pixels2 = (*IPRH_2).cluster()->pixels();
	    clSize2_ = pixels2.size();
	    if(clSize2_ > 500) clSize2_=500;
	 	for(int i =0; i < clSize2_; ++i){
			pixelsX2_[i] = pixels2[i].x;
			pixelsY2_[i] = pixels2[i].y;
	    	pixelsADC2_[i] =  pixels2[i].adc;
	 	}
		charge2_ = (*IPRH_2).cluster()->charge();
		minPixelRow2_= (*IPRH_2).cluster()->minPixelRow();
		minPixelCol2_= (*IPRH_2).cluster()->minPixelCol();
		maxPixelRow2_= (*IPRH_2).cluster()->maxPixelRow();
		maxPixelCol2_= (*IPRH_2).cluster()->maxPixelCol();		
		dX_ = abs(recHitX2_ - recHitX1_);
		dY_ = abs(recHitY2_ - recHitY1_);
        dR_ = sqrt(recHitR2_- recHitR1_);
		
		
		//calculate the min distance--------------						
		minDist_=2000.;
		for(std::vector<SiPixelCluster::Pixel>::const_iterator pixel_1 = pixels1.begin(); pixel_1!= pixels1.end(); pixel_1++) {
			for(std::vector<SiPixelCluster::Pixel>::const_iterator pixel_2 = pixels2.begin(); pixel_2!= pixels2.end(); pixel_2++) {

           	 LocalPoint lp_1 = topol_1->localPosition(MeasurementPoint(pixel_1->x, pixel_1->y));
             GlobalPoint gp_1 = pgdu_1->surface().toGlobal(Local3DPoint(lp_1.x(),lp_1.y(),lp_1.z()));

			 LocalPoint lp_2 = topol_2->localPosition(MeasurementPoint(pixel_2->x, pixel_2->y));
             GlobalPoint gp_2 = pgdu_2->surface().toGlobal(Local3DPoint(lp_2.x(),lp_2.y(),lp_2.z()));

             float distance = (gp_1-gp_2).mag();
             if(distance < minDist_)minDist_ = distance;
		    }//first pixel
		}//second pixel
		
	  	if(minDist_ < 0.1) printf("%f\n", minDist_);
		ClustersDist_->Fill();
	}//rechit2 
  }//reehit1          
}//accept


	
}


// ------------ method called once each job just before starting event loop  ------------
void 
SiPixelRecHitsAna::beginJob()
{
    SubEventN_=0;
    oFile_ = new TFile((const char*)outputFile_.c_str(), "RECREATE");
     
	  
    Clusters_ = new TTree("PixelClusters", "Pixel Clusters");
   	Clusters_->Branch("EventN", &EventN_);
	Clusters_->Branch("SubEventN", &SubEventN_);
	Clusters_->Branch("recHitX", &recHitX1_);
	Clusters_->Branch("recHitY", &recHitY1_);
	Clusters_->Branch("recHitZ", &recHitZ1_);
	Clusters_->Branch("recHitR", &recHitR1_);
	Clusters_->Branch("recHitEta", &recHitEta1_);
	Clusters_->Branch("recHitPhi", &recHitPhi1_);
	Clusters_->Branch("ClusterNumber", &ClusterNumber1_);
    Clusters_->Branch("sizeX", &sizeX1_);
    Clusters_->Branch("sizeY", &sizeY1_);
	Clusters_->Branch("clSize", &clSize1_);
	Clusters_->Branch("pixelsX", pixelsX1_);
	Clusters_->Branch("pixelsY", pixelsY1_);
	Clusters_->Branch("pixelsADC", pixelsADC1_);
    Clusters_->Branch("charge", &charge1_);
    Clusters_->Branch("minPixelRow", &minPixelRow1_);
	Clusters_->Branch("minPixelCol", &minPixelCol1_);
	Clusters_->Branch("maxPixelRow", &maxPixelRow1_);
	Clusters_->Branch("maxPixelCol", &maxPixelCol1_);
       
    ClustersDist_ = new TTree("PixelClustersDistance", "Pixel Clusters Distance (n pixels)");
   	ClustersDist_->Branch("EventN", &EventN_);
	ClustersDist_->Branch("SubEventN", &SubEventN_);
		
    ClustersDist_->Branch("ClusterNumber1", &ClusterNumber1_);
	ClustersDist_->Branch("recHitX1", &recHitX1_);
	ClustersDist_->Branch("recHitY1", &recHitY1_);
    ClustersDist_->Branch("recHitZ1", &recHitZ1_);
	ClustersDist_->Branch("recHitR1", &recHitR1_);
	ClustersDist_->Branch("recHitEta1", &recHitEta1_);
    ClustersDist_->Branch("recHitPhi1", &recHitPhi1_);
    ClustersDist_->Branch("sizeX1", &sizeX1_);
    ClustersDist_->Branch("sizeY1", &sizeY1_);
	ClustersDist_->Branch("clSize1", &clSize1_);
    ClustersDist_->Branch("pixelsX1", pixelsX1_);
	ClustersDist_->Branch("pixelsY1", pixelsY1_);
	ClustersDist_->Branch("pixelsADC1", pixelsADC1_);
    ClustersDist_->Branch("charge1", &charge1_);
    ClustersDist_->Branch("minPixelRow1", &minPixelRow1_);
	ClustersDist_->Branch("minPixelCol1", &minPixelCol1_);
	ClustersDist_->Branch("maxPixelRow1", &maxPixelRow1_);
	ClustersDist_->Branch("maxPixelCol1", &maxPixelCol1_);
    	
    ClustersDist_->Branch("ClusterNumber2", &ClusterNumber2_);
    ClustersDist_->Branch("recHitX2", &recHitX2_);
	ClustersDist_->Branch("recHitY2", &recHitY2_);
    ClustersDist_->Branch("recHitZ2", &recHitZ2_);
	ClustersDist_->Branch("recHitR2", &recHitR2_);
	ClustersDist_->Branch("recHitEta2", &recHitEta2_);
    ClustersDist_->Branch("recHitPhi2", &recHitPhi2_);
    ClustersDist_->Branch("sizeX2", &sizeX2_);
    ClustersDist_->Branch("sizeY2", &sizeY2_);
	ClustersDist_->Branch("clSize2", &clSize2_);
    ClustersDist_->Branch("pixelsX2", pixelsX2_);
	ClustersDist_->Branch("pixelsY2", pixelsY2_);
	ClustersDist_->Branch("pixelsADC2", pixelsADC2_);
    ClustersDist_->Branch("charge2", &charge2_);
    ClustersDist_->Branch("minPixelRow2", &minPixelRow2_);
	ClustersDist_->Branch("minPixelCol2", &minPixelCol2_);
	ClustersDist_->Branch("maxPixelRow2", &maxPixelRow2_);
	ClustersDist_->Branch("maxPixelCol2", &maxPixelCol2_);
    
	ClustersDist_->Branch("dX", &dX_);
    ClustersDist_->Branch("dY", &dY_);
	ClustersDist_->Branch("dR", &dR_);
    ClustersDist_->Branch("minDist", &minDist_);
	
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiPixelRecHitsAna::endJob() {
   oFile_->Write();
   oFile_->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelRecHitsAna);
