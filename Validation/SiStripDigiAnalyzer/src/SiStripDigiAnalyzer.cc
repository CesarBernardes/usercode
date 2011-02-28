// -*- C++ -*-
//
// Package:    SiStripDigiAnalyzer
// Class:      SiStripDigiAnalyzer
// 
/**\class SiStripDigiAnalyzer SiStripDigiAnalyzer.cc Validation/SiStripAnalyzer/src/SiStripDigiAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ivan Amos Cali
//         Created:  Mon Jul 28 14:10:52 CEST 2008
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
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h" 
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "DataFormats/SiStripDigi/interface/SiStripProcessedRawDigi.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripPedestalsSubtractor.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripCommonModeNoiseSubtractor.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripRawProcessingFactory.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"


//
// class decleration
//

class SiStripDigiAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SiStripDigiAnalyzer(const edm::ParameterSet&);
      ~SiStripDigiAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
   edm::InputTag src_;
   edm::InputTag srcProcessedRawDigi_;
   edm::InputTag srcCM_;
   CentralityProvider * centrality_;
   
   std::auto_ptr<SiStripPedestalsSubtractor>   subtractorPed_;
   
  //const TrackerGeometry& trGeo_;
   
   TTree* Digis_;
   TTree* AnaInfo_;   
   
   bool doComparisonRaw_;
   bool PassedOccTh_;
   uint32_t subDetId_;
   uint32_t detId_;
   uint32_t EventN_;
   uint32_t RunN_;
   uint32_t CentralityBin_;
   uint32_t MarkedAsBadInZs_;
   double detZ_;
   double detR_;
   double detEta_;
   double detPhi_;
   uint32_t Nstrips_;
   uint32_t NAPVs_;   
   double APVoccupancy_[6];
   double APVcharge_[6];
   uint16_t stripDigis_[768];
   int16_t RawDigis_[768];
   float APVCM_[6];
   
  
   uint32_t TotEvProcessed_;
   uint32_t TotAPVProcessed_;
   uint32_t TotAPVSelected_;

};


SiStripDigiAnalyzer::SiStripDigiAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  
  src_ =  iConfig.getParameter<edm::InputTag>( "src" );
  subtractorPed_ = SiStripRawProcessingFactory::create_SubtractorPed(iConfig.getParameter<edm::ParameterSet>("Algorithms"));
  srcProcessedRawDigi_ =  iConfig.getParameter<edm::InputTag>( "srcProcessedRawDigi" );
  srcCM_ =  iConfig.getParameter<edm::InputTag>( "srcAPVCM" );
  doComparisonRaw_ =  iConfig.getParameter<bool>( "doComparisonRaw" );
   
   centrality_ = 0;
   
   edm::Service<TFileService> fs;
  
   Digis_ = fs->make<TTree>("Digis", "Digis");
   Digis_->Branch("EventN", &EventN_);
   Digis_->Branch("RunN", &RunN_);
   Digis_->Branch("CentralityBin", &CentralityBin_);
  
   Digis_->Branch("subDetId",&subDetId_);
   Digis_->Branch("detId", &detId_);
   Digis_->Branch("PassedOccTh",&PassedOccTh_);
   Digis_->Branch("MarkedAsBadInZs",&MarkedAsBadInZs_);
   Digis_->Branch("detZ", &detZ_);
   Digis_->Branch("detR", &detR_);
   Digis_->Branch("detEta", &detEta_);
   Digis_->Branch("detPhi", &detPhi_);
   Digis_->Branch("nstrips", &Nstrips_);
   Digis_->Branch("NAPVs", &NAPVs_);
   Digis_->Branch("APVCM", APVCM_);
   Digis_->Branch("APVoccupancy", APVoccupancy_);
   Digis_->Branch("APVcharge", APVcharge_);
   Digis_->Branch("Digis", stripDigis_);
   Digis_->Branch("RawDigis", RawDigis_);
  

   AnaInfo_ = fs->make<TTree>("AnaInfo", "AnaInfo");
   AnaInfo_->Branch("TotEvProcessed", &TotEvProcessed_);
   AnaInfo_->Branch("TotAPVProcessed", &TotAPVProcessed_);
   AnaInfo_->Branch("TotAPVSelected", &TotAPVSelected_);
}


SiStripDigiAnalyzer::~SiStripDigiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiStripDigiAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
   using namespace edm;
   ++TotEvProcessed_;
   //if(!centrality_) centrality_ = new CentralityProvider(es);

   //centrality_->newEvent(e,es); // make sure you do this first in every event
   //CentralityBin_ = centrality_->getBin();
   ////bin = bin/4;
   
   
   
   RunN_  = e.id().run();
   EventN_ = e.id().event();

  //int lumiBlock = e.luminosityBlock();
  //uint32_t bx        = e.bunchCrossing();
  //int orbit     = e.orbitNumber();
   
   
   edm::ESHandle<TrackerGeometry> tracker;
   es.get<TrackerDigiGeometryRecord>().get( tracker );
   const TrackerGeometry &trGeo(*tracker);
   

  // std::string digiProducer = "siStripDigis";
   edm::Handle<edm::DetSetVector<SiStripDigi> > stripDigis;
   e.getByLabel(src_, stripDigis);
  
   subtractorPed_->init(es);
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > moduleRawDigi;
   e.getByLabel(srcProcessedRawDigi_,moduleRawDigi);
   
   edm::Handle<edm::DetSetVector<SiStripProcessedRawDigi> > moduleCM;
   e.getByLabel(srcCM_, moduleCM);
   
   //Processing Digis
   //-------------------------------------------------------------------------------
   edm::DetSetVector<SiStripDigi>::const_iterator DSViter = stripDigis->begin();
   for( ; DSViter != stripDigis->end(); DSViter++) {
	 uint32_t detId_ = DSViter->id;
	 DetId  detId(detId_);
	 subDetId_ = detId.subdetId();
     
    
     //detector Geometry-- the coordinate are of the module center                                                           
	 const StripGeomDetUnit* StripModuleGeom =(const StripGeomDetUnit*)trGeo.idToDetUnit(detId);   //detector geometry -> it returns the center of the module                                                                                      
	 detZ_ = StripModuleGeom->surface().position().z();        //module z                                              
	 detR_ = StripModuleGeom->surface().position().perp();        //module R                                           
	 detEta_ = StripModuleGeom->surface().position().eta();    //module eta                                            
	 detPhi_ = StripModuleGeom->surface().position().phi();    //module phi                                            
	 Nstrips_ = StripModuleGeom->specificTopology().nstrips(); //n of strips    
	 NAPVs_ = Nstrips_/128;
	 TotAPVProcessed_+=NAPVs_;
  
   
     uint16_t i;
     for(i=0; i<6; ++i){
		APVoccupancy_[i]=0;
        APVcharge_[i]=0;
	 }
     for(i=0; i<768; ++i){
	   stripDigis_[i]=0;
	   RawDigis_[i]=0;
     }
	
	 //Loop over Digis------------------------------------------
	 edm::DetSet<SiStripDigi>::const_iterator  iter = DSViter->data.begin();
     for( ;iter != DSViter->data.end();++iter ) {  
		uint16_t strip = (*iter).strip();
	    uint16_t adc = (*iter).adc();
	    uint16_t APVn = strip/ 128;
	    APVoccupancy_[APVn] += 1./128.;
        APVcharge_[APVn] +=adc;
		stripDigis_[strip] = adc;
	 }
      
     PassedOccTh_=false;
     for(i=0; i< NAPVs_; ++i){
		if(APVoccupancy_[i]>0.2){
    		PassedOccTh_=true;
			++TotAPVSelected_;
		}
     }
     
	 
	 //looking for raw data
	 MarkedAsBadInZs_ = false;
     if(PassedOccTh_==true&&doComparisonRaw_){
				
		edm::DetSetVector<SiStripRawDigi>::const_iterator itRawDigis = moduleRawDigi->begin();
		for (; itRawDigis != moduleRawDigi->end(); ++itRawDigis) {
           uint32_t detIdraw = itRawDigis->id;
		   
		   if(detIdraw == detId_){
			MarkedAsBadInZs_ = true;
			
			std::vector<int16_t> ProcessedRawDigis(itRawDigis->size());
			subtractorPed_->subtract( *itRawDigis, ProcessedRawDigis);
		
	        std::vector<int16_t>::const_iterator itProcessedRawDigis;

			int strip =0;
			for(itProcessedRawDigis = ProcessedRawDigis.begin();itProcessedRawDigis != ProcessedRawDigis.end(); ++itProcessedRawDigis){
				RawDigis_[strip] = *itProcessedRawDigis;
			    ++strip;
			}	  
			itRawDigis = moduleRawDigi->end() -1;
			
		   }
		}
   
	  }
	  
	  //Adding the calculated CM
	  for(uint16_t i =0; i < 6; ++i) APVCM_[i]=0;
	  edm::DetSetVector<SiStripProcessedRawDigi>::const_iterator itAPVCM = moduleCM->begin();
	  for (; itAPVCM != moduleCM->end(); ++itAPVCM) {
	    uint32_t detIdCM = itAPVCM->id;
		if(detIdCM == detId_){
			edm::DetSet<SiStripProcessedRawDigi>::const_iterator itCM = itAPVCM->begin();
			uint16_t APVn =0;
			for (; itCM != itAPVCM->end(); ++itCM){
				APVCM_[APVn] = itCM->adc();
				++APVn;
			}
			itAPVCM = moduleCM->end() -1;
		}
	  }
	  
    Digis_->Fill();		
  }
}


// ------------ method called once each job just before starting event loop  ------------
void SiStripDigiAnalyzer::beginJob()
{
  TotEvProcessed_=0;
  TotAPVProcessed_=0;
  TotAPVSelected_=0;
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripDigiAnalyzer::endJob() {
  AnaInfo_->Fill();
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripDigiAnalyzer);
