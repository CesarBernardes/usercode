// -*- C++ -*-
//
// Package:    SiStripAnalyzerRawLight
// Class:      SiStripAnalyzerRawLight
// 
/**\class SiStripAnalyzerRawLight SiStripAnalyzerRawLight.cc Validation/SiStripAnalyzer/src/SiStripAnalyzerRawLight.cc

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
//here
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
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


//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"


//
// class decleration
//

class SiStripAnalyzerRawLight : public edm::EDAnalyzer {
   public:
      explicit SiStripAnalyzerRawLight(const edm::ParameterSet&);
      ~SiStripAnalyzerRawLight();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
   std::string outputFile_;
   edm::InputTag src_;
   edm::InputTag srcOrig_;
   edm::InputTag srcCMN_;
   
  //const TrackerGeometry& trGeo_;
   TFile* oFile_;
   TNtuple* ModuleNt_;
   TNtuple* DigisNt_;
   TNtuple* DigisOrigNt_;
   TNtuple* DigisPerocAndOrigNt_;
   TNtuple* CMNNt_;

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
SiStripAnalyzerRawLight::SiStripAnalyzerRawLight(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outputFile", "StripHistos.root");
  src_ =  iConfig.getParameter<edm::InputTag>( "src" );

}


SiStripAnalyzerRawLight::~SiStripAnalyzerRawLight()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiStripAnalyzerRawLight::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get( tracker );
   const TrackerGeometry &trGeo(*tracker);
   //   trGeo_ = tracker.product();


  // std::string digiProducer = "siStripDigis";
   //here
   edm::Handle<edm::DetSetVector<SiStripRawDigi> > stripDigis;
   iEvent.getByLabel(src_, stripDigis);
  
   
   //Processing Digis
   //-------------------------------------------------------------------------------
   //can also be raw
   //here
   edm::DetSetVector<SiStripRawDigi>::const_iterator DSViter = stripDigis->begin();
   
   std::map<uint32_t, std::vector<uint16_t> > APVMultMap;
   APVMultMap.clear();
        
   std::map<uint32_t, std::vector<uint16_t> > ModuleDigisMap;
   ModuleDigisMap.clear();
   
   for( ; DSViter != stripDigis->end(); DSViter++) {
     uint32_t id = DSViter->id;
     float idL = id & 0xFFFF;
     float idH = (id & 0xFFFF0000) >> 16;
     DetId  detId(id);
     //here
     edm::DetSet<SiStripRawDigi>::const_iterator  begin = DSViter->data.begin();
     edm::DetSet<SiStripRawDigi>::const_iterator  end   = DSViter->data.end();
     edm::DetSet<SiStripRawDigi>::const_iterator  iter;
    
     //detector Geometry-- the coordinate are of the module center                                                           
     const StripGeomDetUnit* StripModuleGeom =(const StripGeomDetUnit*)trGeo.idToDetUnit(detId);   //detector geometry -> it returns the center of the module                                                                                      
     double detZ = StripModuleGeom->surface().position().z();        //module z                                              
     double detR = StripModuleGeom->surface().position().perp();        //module R                                           
     double detEta = StripModuleGeom->surface().position().eta();    //module eta                                            
     double detPhi = StripModuleGeom->surface().position().phi();    //module phi                                            
     unsigned int Nstrips = StripModuleGeom->specificTopology().nstrips(); //n of strips    
     unsigned int NAPVs = Nstrips/128;
     
     uint32_t SubDet = detId.subdetId();        
	 unsigned int ModuleDigisNz =0, ModuleDigis =0;
	
	 float* ModuleDigisVect= new float[Nstrips];
	 for(unsigned int i=0; i<Nstrips; ++i) ModuleDigisVect[i]=0;
	 
	 std::vector<float> TotADCCharge;
	 TotADCCharge.insert(TotADCCharge.begin(), NAPVs, 0);
	 
	 std::vector<float> APVMultVect;
	 APVMultVect.insert(APVMultVect.begin(), NAPVs, 0);     
	 //Loop over Digis------------------------------------------
	 int strip =0;
   	 for ( iter = begin ; iter != end; iter++ ) {  
	      int16_t adc = (*iter).adc();
	      uint16_t APVn = strip/ 128;
	      ModuleDigisVect[strip] = adc;
		  TotADCCharge[APVn] += adc;
	      if(adc >0){
             ++ModuleDigisNz; 
             APVMultVect[APVn]++;
          }
	      ++ModuleDigis;
	      ++strip;
	      DigisNt_->Fill(idH,idL, SubDet, detZ, detR, detPhi, detEta,strip,adc, APVn);

     }
	        
      ModuleNt_->Fill(idH,idL, SubDet,  detZ, detR, detPhi, detEta, ModuleDigisNz, ModuleDigis, Nstrips);
   
     for(unsigned int APVn =0; APVn< NAPVs; ++APVn){
   	    Double_t APVMedian = TMath::Median(128, ModuleDigisVect+128*APVn);
        CMNNt_->Fill(idH, idL, SubDet, NAPVs, APVn,APVMedian ,TotADCCharge[APVn], APVMultVect[APVn] );
     }  
   
   
   }

 
     
}


// ------------ method called once each job just before starting event loop  ------------
void SiStripAnalyzerRawLight::beginJob()
{
  
  oFile_ = new TFile((const char*)outputFile_.c_str(), "RECREATE");
  
  ModuleNt_ = new TNtuple("TrackerModules", "TrackerModules","idH:idL:SubDet:z:R:phi:eta:MhitsNz:MHits:nstrips",1000);
  DigisNt_ = new TNtuple("TrackerDigis", "TrackerDigis","idH:idL:SubDet:z:R:phi:eta:strip:adc:APVn",1000);
  CMNNt_ = new TNtuple("TrackerCMN", "TrackerCMN", "idH:idL:SubDet:ModAPV:APVn:CMN:FullADC:APVMult", 1000);
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripAnalyzerRawLight::endJob() {
  oFile_->Write();
  oFile_->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripAnalyzerRawLight);

//  SiStripCMN CNoise(APVn, CM);   Median
//  SiStripCMN CNoise(APVn, offset, slope);   fast Linear
//  SiStripCMN CNoise(APVn, CM, FixedBias);   TT6 

//SubDetector { TIB=3,TID=4,TOB=5,TEC=6 }
/* TOBDetId(uint32_t layer,
            uint32_t rod_fw_bw,
            uint32_t rod,
            uint32_t module,
            uint32_t ster) 
            
    TIBDetId(uint32_t layer,
            uint32_t str_fw_bw,
            uint32_t str_int_ext,
            uint32_t str,
            uint32_t module,
            uint32_t ster)
 
    TIDDetId(uint32_t side,
           uint32_t wheel,
           uint32_t ring,
           uint32_t module_fw_bw,
           uint32_t module,
           uint32_t ster)   
           
    TECDetId(uint32_t side,
           uint32_t wheel,
           uint32_t petal_fw_bw,
           uint32_t petal,
           uint32_t ring,
           uint32_t module,
           uint32_t ster)
*/

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
