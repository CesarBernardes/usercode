// -*- C++ -*-
//
// Package:    SiStripRecHitToDetId
// Class:      SiStripRecHitToDetId
// 
/**\class SiStripRecHitToDetId SiStripRecHitToDetId.cc Validation/SiStripRecHitToDetId/src/SiStripRecHitToDetId.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ivan Amos Cali,32 4-A08,+41227673039,
//         Created:  Sat Nov 13 21:24:20 CET 2010
// $Id: SiStripRecHitToDetId.cc,v 1.1 2011/02/28 09:27:40 icali Exp $
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

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"



//
// class declaration
//

class SiStripRecHitToDetId : public edm::EDAnalyzer {
   public:
      explicit SiStripRecHitToDetId(const edm::ParameterSet&);
      ~SiStripRecHitToDetId();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
	  std::vector<edm::InputTag> inputTagsRecHits_;
	  edm::InputTag inputTagRestoredAPVs_;
	  
	 
    


};


SiStripRecHitToDetId::SiStripRecHitToDetId(const edm::ParameterSet& conf){


   inputTagsRecHits_ = conf.getParameter<std::vector<edm::InputTag> >("RecHitCollections");
  
   inputTagRestoredAPVs_ = conf.getParameter<edm::InputTag>("RestoredAPVCollections");
      
  
  
   
}


SiStripRecHitToDetId::~SiStripRecHitToDetId()
{
 

}


void
SiStripRecHitToDetId::analyze(const edm::Event& e, const edm::EventSetup& es)
{
   using namespace edm;
   
    
  //create bad APVs map
  edm::Handle< edm::DetSetVector<SiStripRawDigi> > inputraw;
  e.getByLabel(inputTagRestoredAPVs_,inputraw);
	
  //create the map with the bad APVs  
	std::map< uint32_t, std::vector < uint16_t> > BadAPVsMap;
	BadAPVsMap.clear();
  	
	if (inputraw->size()){
		for ( edm::DetSetVector<SiStripRawDigi>::const_iterator rawDigis = inputraw->begin(); rawDigis != inputraw->end(); ++rawDigis) {
    	   
			edm::DetSet<SiStripRawDigi>::const_iterator itRawDigis = rawDigis->begin();
			uint16_t nAPV = rawDigis->size()/128;
			uint32_t rawDetId = rawDigis->id;
          
			
			std::vector<uint16_t> restoredAPV;
			restoredAPV.clear();
			
			std::vector<bool> tmprestoredAPV;
			tmprestoredAPV.clear();
			tmprestoredAPV.insert(tmprestoredAPV.begin(), nAPV, false); 
          	
          	bool foundRestoredAPV = false;
			for( uint16_t strip =0; strip < rawDigis->size();++strip){
				if(itRawDigis[strip].adc()!=0){
				  tmprestoredAPV[strip/128] = true;       
				  foundRestoredAPV = true;
				}
			}
 
            for(size_t itAPV=0; itAPV < tmprestoredAPV.size(); ++itAPV){
                if(tmprestoredAPV[itAPV])  restoredAPV.push_back(itAPV);
                
            }
   		  	
   		  	BadAPVsMap.insert(BadAPVsMap.end(), std::pair< uint32_t, std::vector < uint16_t> >(rawDetId, restoredAPV));
		}//loop over raw data collection
						
    }//if inputraw.size
 
  
  //0000000000000000000000000000000000000000000000000000000000000000000000000000000
  //loop over rec-hits and find bad APVs    
  for(std::vector<edm::InputTag>::const_iterator inputTagRH = inputTagsRecHits_.begin(); inputTagRH != inputTagsRecHits_.end(); ++inputTagRH) {
   edm::Handle<SiStripRecHit2DCollection> siStripRecCollection;
   e.getByLabel(*inputTagRH,siStripRecCollection);
   const SiStripRecHit2DCollection RecHits = *(siStripRecCollection.product());
   
     
   SiStripRecHit2DCollection::DataContainer::const_iterator itRecHits = RecHits.data().begin();
   for(; itRecHits != RecHits.data().end();++itRecHits) {
		uint32_t RecHitDetId= (*itRecHits).geographicalId();
	    uint32_t clSize= (*itRecHits).cluster()->amplitudes().size();
	 	uint32_t clFirstStrip = (*itRecHits).cluster()->firstStrip();
		    
		
		
		//is RecHit APV bad?
		if(BadAPVsMap.find(RecHitDetId) != BadAPVsMap.end()){
			std::vector<uint16_t> ClusterAPV;
		    ClusterAPV.clear();
		
		    uint16_t FirstAPV = clFirstStrip/128;
		    uint16_t SecondAPV = (clFirstStrip+clSize)/128;
		    ClusterAPV.push_back(FirstAPV);
		    if(SecondAPV != FirstAPV) ClusterAPV.push_back(SecondAPV);
		    
		    //here I can put the code that I want........
		    
		
		}
	 	
   }
 
 }	
}



void 
SiStripRecHitToDetId::beginJob()
{
}


void 
SiStripRecHitToDetId::endJob() {
}


DEFINE_FWK_MODULE(SiStripRecHitToDetId);
