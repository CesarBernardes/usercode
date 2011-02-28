// -*- C++ -*-
//
// Package:    SiStripBadAPVAnalyzer
// Class:      SiStripBadAPVAnalyzer
// 
/**\class SiStripBadAPVAnalyzer SiStripBadAPVAnalyzer.cc Validation/SiStripAnalyzer/src/SiStripBadAPVAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ivan Amos Cali
//         Created:  Mon Jul 28 14:10:52 CEST 2008
// $Id: SiStripBadAPVAnalyzer.cc,v 1.1 2010/11/04 15:29:18 edwenger Exp $
//
//
 

// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"

#include "DataFormats/SiStripDigi/interface/SiStripProcessedRawDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"

#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripPedestalsSubtractor.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripCommonModeNoiseSubtractor.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripRawProcessingFactory.h"

#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileDirectory.h"


//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "THStack.h"


//
// class decleration
//

class SiStripBadAPVAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SiStripBadAPVAnalyzer(const edm::ParameterSet&);
      ~SiStripBadAPVAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
	  void initHistos();
	  
	  std::auto_ptr<SiStripPedestalsSubtractor>   subtractorPed_;
      CentralityProvider * centrality_;
	  uint32_t CentralityBin_;
	  
	  
	  edm::InputTag srcBaseline_;
	  edm::InputTag srcProcessedRawDigi_;
      
	  edm::Service<TFileService> fs_;
	  
	  
	  TH1F* h1BadAPVperEvent_;
	  TH1F* h1BadAPVsPerModule_;
	  TH1F* h1ClusterSizeAllAPVs_;
	  TH1F* h1ClusterChargeAllAPVs_;
	  TH1F* h1ClusterSizeBadAPVs_;
	  TH1F* h1ClusterChargeBadAPVs_;
	  TH1F* h1ClusterSizeGoodAPVs_;
	  TH1F* h1ClusterChargeGoodAPVs_;
	  TH1F* h1NHIPinStereo_;
	  TH1F* h1SatStripsPerAPVAll_;
	  TH1F* h1SatStripsPerAPVBad_;
	  TH1F* h1SatStripsPerAPVGood_;
	       	  
	  TH2F* h2BadAPVVsEtaVsR_;	  
	  TH2F* h2BadAPVperEventVsCent_;
	    		
	  TCanvas* Canvas_;
	  TH1F* h1ProcessedRawDigis_;
	  TH1F* h1Baseline_;
	  TH1F* h1Clusters_;
	   	  
	  uint16_t nModuletoDisplay_;
	  uint16_t actualModule_;
};


SiStripBadAPVAnalyzer::SiStripBadAPVAnalyzer(const edm::ParameterSet& conf){
   
   centrality_ = 0;
   actualModule_ =0;
     
  srcBaseline_ =  conf.getParameter<edm::InputTag>( "srcBaseline" );
  srcProcessedRawDigi_ =  conf.getParameter<edm::InputTag>( "srcProcessedRawDigi" );
  subtractorPed_ = SiStripRawProcessingFactory::create_SubtractorPed(conf.getParameter<edm::ParameterSet>("Algorithms"));
  nModuletoDisplay_ = conf.getParameter<uint32_t>( "nModuletoDisplay" );
 
  
  
  this->initHistos();
  
}


SiStripBadAPVAnalyzer::~SiStripBadAPVAnalyzer()
{
 
   

}

void
SiStripBadAPVAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
   using namespace edm;
   
   
   subtractorPed_->init(es);
   
   if(!centrality_) centrality_ = new CentralityProvider(es);
   centrality_->newEvent(e,es); // make sure you do this first in every event
   CentralityBin_ = centrality_->getBin();
   std::cout << "cen " << CentralityBin_ << std::endl;
   
   edm::ESHandle<TrackerGeometry> tracker;
   es.get<TrackerDigiGeometryRecord>().get( tracker );
   const TrackerGeometry &trGeo(*tracker);
   
   
   edm::Handle<edm::DetSetVector<SiStripProcessedRawDigi> > moduleBaseline;
   e.getByLabel(srcBaseline_, moduleBaseline);
   
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > moduleRawDigi;
   e.getByLabel(srcProcessedRawDigi_,moduleRawDigi);

   edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusters;
   edm::InputTag clusLabel("siStripClusters");
   e.getByLabel(clusLabel, clusters);

   
   char detIds[20];
   char evs[20];
   char runs[20];    
    
   edm::DetSetVector<SiStripProcessedRawDigi>::const_iterator itDSBaseline = moduleBaseline->begin();
   edm::DetSetVector<SiStripRawDigi>::const_iterator itRawDigis = moduleRawDigi->begin();
   
   //uint32_t NBabAPVs = moduleRawDigi->size();   
   uint32_t NBabAPVs =0;   
    
   TFileDirectory sdProcessedRawDigis_= fs_->mkdir("ProcessedRawDigis");
   TFileDirectory sdBaseline_= fs_->mkdir("Baseline");
   TFileDirectory sdClusters_= fs_->mkdir("Clusters");   
   for (; itRawDigis != moduleRawDigi->end(); ++itRawDigis, ++itDSBaseline) {
     // if(actualModule_ > nModuletoDisplay_) return;
      uint32_t detId = itRawDigis->id;
	  
	 const StripGeomDetUnit* StripModuleGeom =(const StripGeomDetUnit*)trGeo.idToDetUnit(detId);   //detector geometry -> it returns the center of the module                                                                                      
	// double detZ = StripModuleGeom->surface().position().z();        //module z                                              
	 double detR = StripModuleGeom->surface().position().perp();        //module R                                           
	 double detEta = StripModuleGeom->surface().position().eta();    //module eta                                            
	 //double detPhi = StripModuleGeom->surface().position().phi();    //module phi                                            
	 //uint32_t Nstrips = StripModuleGeom->specificTopology().nstrips(); //n of strips    
	 //NAPVs = Nstrips_/128;
	  
	  h2BadAPVVsEtaVsR_->Fill(detEta, detR);
	  
      if(itDSBaseline->id != detId){
		std::cout << "Collections out of Synch. Something of fishy is going on ;-)" << std::endl;
		return;
      }	  
      
	  actualModule_++;
	  if(actualModule_ < nModuletoDisplay_){
	  uint32_t event = e.id().event();
	  uint32_t run = e.id().run();
	  //std::cout << "processing module N: " << actualModule_<< " detId: " << detId << " event: "<< event << std::endl; 

      //Creating histograms for baseline, clusters and raw Digis
      //--------------------------------------------------------------------------------------------------------------------------	   
	  sprintf(detIds,"%ul", detId);
	  sprintf(evs,"%ul", event);
	  sprintf(runs,"%ul", run);
	  char* dHistoName = Form("Id:%s_run:%s_ev:%s",detIds, runs, evs);
	  h1ProcessedRawDigis_ = sdProcessedRawDigis_.make<TH1F>(dHistoName,dHistoName, 768, -0.5, 767.5); 
	  h1Baseline_ = sdBaseline_.make<TH1F>(dHistoName,dHistoName, 768, -0.5, 767.5); 
      h1Clusters_ = sdClusters_.make<TH1F>(dHistoName,dHistoName, 768, -0.5, 767.5);
	  
	  
	  h1ProcessedRawDigis_->SetXTitle("strip#"); h1ProcessedRawDigis_->SetYTitle("ADC");  h1ProcessedRawDigis_->SetMaximum(1024.);   h1ProcessedRawDigis_->SetMinimum(-300.); h1ProcessedRawDigis_->SetLineWidth(2);
      h1Baseline_->SetXTitle("strip#"); h1Baseline_->SetYTitle("ADC"); h1Baseline_->SetMaximum(1024.); h1Baseline_->SetMinimum(-300.);  h1Baseline_->SetLineWidth(2);  h1Baseline_->SetLineStyle(2); h1Baseline_->SetLineColor(2);
	  h1Clusters_->SetXTitle("strip#"); h1Clusters_->SetYTitle("ADC"); h1Clusters_->SetMaximum(1024.); h1Clusters_->SetMinimum(-300.);  h1Clusters_->SetLineWidth(2); h1Clusters_->SetLineStyle(2); h1Clusters_->SetLineColor(3);
	  }
	  
	
      //Filling up the baseline and Raw Digis Histos
	  //-----------------------------------------------------------------------------------------------------------------------------------------------
	  edm::DetSet<SiStripProcessedRawDigi>::const_iterator  itBaseline; 
	  std::vector<int16_t>::const_iterator itProcessedRawDigis;
	  
	  
	  std::vector<int16_t> ProcessedRawDigis(itRawDigis->size());
	  subtractorPed_->subtract( *itRawDigis, ProcessedRawDigis);
	  
	  
	  int strip =0, satstrip=0;
	  std::vector<uint8_t> vRestoredAPV;
	  vRestoredAPV.clear();
	  bool Restored = false; 
	  float oldbasADC= itDSBaseline->begin()->adc();
	  for(itProcessedRawDigis = ProcessedRawDigis.begin(), itBaseline = itDSBaseline->begin();itProcessedRawDigis != ProcessedRawDigis.end(); ++itProcessedRawDigis, ++itBaseline){
		float rawDigADC = *itProcessedRawDigis;
		if(rawDigADC>512) ++satstrip;
		h1ProcessedRawDigis_->Fill(strip, rawDigADC);
		float basADC = itBaseline->adc();
		h1Baseline_->Fill(strip,basADC ); 
		if(strip %128 ==0) oldbasADC = basADC;
		if(basADC != oldbasADC) Restored = true;
		if(strip % 127 ==0){
		   h1SatStripsPerAPVAll_->Fill(satstrip);
		   if(Restored == true){
		     h1SatStripsPerAPVBad_->Fill(satstrip);
			 vRestoredAPV.push_back((int)(strip/128));
			 Restored = false;
		   } else{
		      h1SatStripsPerAPVGood_->Fill(satstrip);
		   }
		   satstrip =0;
		}
		++strip;
      }	  
	  
	  NBabAPVs += vRestoredAPV.size();
	 // std::cout << NBabAPVs << std::endl;
	  h1BadAPVsPerModule_->Fill(vRestoredAPV.size());
	  
	  //looping over clusters
	  //------------------------------------------------------------------------------------------------------------------------------------------------------------
	  edmNew::DetSetVector<SiStripCluster>::const_iterator itClusters = clusters->begin();
	  for ( ; itClusters != clusters->end(); ++itClusters ){
	  	for ( edmNew::DetSet<SiStripCluster>::const_iterator clus =	itClusters->begin(); clus != itClusters->end(); ++clus){
            if(clus->geographicalId() == detId){
				      	
     			strip=clus->firstStrip();
				uint16_t stAPVn = strip/128, APVn;
				uint32_t clSize=0, clCharge=0; 
				for( std::vector<uint8_t>::const_iterator itAmpl = clus->amplitudes().begin(); itAmpl != clus->amplitudes().end(); ++itAmpl){
                   APVn = strip/128;
				   if(APVn != stAPVn){
				     bool Restored =false;
				     for(size_t i=0; i< vRestoredAPV.size(); ++i)if(stAPVn==vRestoredAPV[i]) Restored=true;
					 if(Restored){
					   h1ClusterSizeBadAPVs_->Fill(clSize);
					   h1ClusterChargeBadAPVs_->Fill(clCharge);
					 }else{
					   h1ClusterSizeGoodAPVs_->Fill(clSize);
					   h1ClusterChargeGoodAPVs_->Fill(clCharge);
					 }
					 h1ClusterSizeAllAPVs_->Fill(clSize);
					 h1ClusterChargeAllAPVs_->Fill(clCharge);
				   
         		     clSize=0;
				     clCharge=0;
				     stAPVn=APVn;
				   }
				   h1Clusters_->Fill(strip, *itAmpl);
				   ++clSize;
				   clCharge += *itAmpl;
				   
				   ++strip;
				}
				
				if(clSize>0){  //chCharge not needed because already implemented
				  bool Restored =false;
				  for(size_t i=0; i< vRestoredAPV.size(); ++i) if(APVn==vRestoredAPV[i])Restored = true;
					
					if(Restored){
					   h1ClusterSizeBadAPVs_->Fill(clSize);
					   h1ClusterChargeBadAPVs_->Fill(clCharge);
				   	}else{
					   h1ClusterSizeGoodAPVs_->Fill(clSize);
					   h1ClusterChargeGoodAPVs_->Fill(clCharge);
					}
					
				  h1ClusterSizeAllAPVs_->Fill(clSize);
				  h1ClusterChargeAllAPVs_->Fill(clCharge);
				}
         		    
				   
				
			}              
    	}
	  }
	  
	 	 	 
	 
	}
	
	std::cout << "Number of module with HIP in this event: " << NBabAPVs << std::endl;
    h1BadAPVperEvent_->Fill(NBabAPVs);	
    h2BadAPVperEventVsCent_->Fill(CentralityBin_, NBabAPVs);
}


// ------------ method called once each job just before starting event loop  ------------
void SiStripBadAPVAnalyzer::beginJob()
{  
   
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripBadAPVAnalyzer::endJob() {
   
	
	
	
  
}

void SiStripBadAPVAnalyzer::initHistos(){


   
  h1BadAPVperEvent_ = fs_->make<TH1F>("BadModulesPerEvent","BadAPVPerEvent", 20001, -0.5, 20000.5); //done
  h1BadAPVperEvent_->SetXTitle("# Modules with Bad APVs");
  h1BadAPVperEvent_->SetYTitle("Entries");
  h1BadAPVperEvent_->SetLineWidth(2);
  h1BadAPVperEvent_->SetLineStyle(2);
  
  h1BadAPVsPerModule_ = fs_->make<TH1F>("BadAPVsPerModule", "BadAPVsPerModule", 7, -0.5, 6.5) ;  //done
  h1BadAPVsPerModule_->SetXTitle("# Bad APVs per module");
  h1BadAPVsPerModule_->SetYTitle("Entries");
  h1BadAPVsPerModule_->SetLineWidth(2);
  h1BadAPVsPerModule_->SetLineStyle(2);
	  
  h1NHIPinStereo_= fs_->make<TH1F>("NHIPinStereo", "NHIPinStereo", 11, -0.5, 10.5);
  h1NHIPinStereo_->SetXTitle("# HIPs in Stereo");
  h1NHIPinStereo_->SetYTitle("Entries");
  h1NHIPinStereo_->SetLineWidth(2);
  h1NHIPinStereo_->SetLineStyle(2);
  
  h1SatStripsPerAPVAll_= fs_->make<TH1F>("SatStripPerAPVAll", "SatStripPerAPV All", 41, -0.5, 40.5);     //done
  h1SatStripsPerAPVAll_->SetXTitle("# Saturatd Strip / APV");
  h1SatStripsPerAPVAll_->SetYTitle("Entries");
  h1SatStripsPerAPVAll_->SetLineWidth(2);
  h1SatStripsPerAPVAll_->SetLineStyle(2);
  
  h1SatStripsPerAPVBad_= fs_->make<TH1F>("SatStripPerAPVBad", "SatStripPerAPV Bad", 41, -0.5, 40.5);  //done
  h1SatStripsPerAPVBad_->SetXTitle("# Saturatd Strip / APV");
  h1SatStripsPerAPVBad_->SetYTitle("Entries");
  h1SatStripsPerAPVBad_->SetLineWidth(2);
  h1SatStripsPerAPVBad_->SetLineStyle(2);
  
  h1SatStripsPerAPVGood_= fs_->make<TH1F>("SatStripPerAPVGood", "SatStripPerAPV Good", 41, -0.5, 40.5); //done
  h1SatStripsPerAPVGood_->SetXTitle("# Saturatd Strip / APV");
  h1SatStripsPerAPVGood_->SetYTitle("Entries");
  h1SatStripsPerAPVGood_->SetLineWidth(2);
  h1SatStripsPerAPVGood_->SetLineStyle(2);
  
  h1ClusterSizeAllAPVs_= fs_->make<TH1F>("ClusterSizeAllAPV", "Cluster Size All APVs", 128, 1, 128); //done
  h1ClusterSizeAllAPVs_->SetXTitle("Cluster Size [strips]");
  h1ClusterSizeAllAPVs_->SetYTitle("Entries");
  h1ClusterSizeAllAPVs_->SetLineWidth(2);
  h1ClusterSizeAllAPVs_->SetLineStyle(2);
  
  h1ClusterChargeAllAPVs_= fs_->make<TH1F>("ClusterChargeAllAPV", "Cluster Charge All APVs" , 30001, 1., 30000.); //done
  h1ClusterChargeAllAPVs_->SetXTitle("Cluster charge [adc]");
  h1ClusterChargeAllAPVs_->SetYTitle("Entries");
  h1ClusterChargeAllAPVs_->SetLineWidth(2);
  h1ClusterChargeAllAPVs_->SetLineStyle(2);
  
  h1ClusterSizeBadAPVs_= fs_->make<TH1F>("ClusterSizeBadAPV", "Cluster Size Bad APVs", 128, 1, 128); //done
  h1ClusterSizeBadAPVs_->SetXTitle("Cluster Size [strips]");
  h1ClusterSizeBadAPVs_->SetYTitle("Entries");
  h1ClusterSizeBadAPVs_->SetLineWidth(2);
  h1ClusterSizeBadAPVs_->SetLineStyle(2);
  
  h1ClusterChargeBadAPVs_= fs_->make<TH1F>("ClusterChargeBadAPV", "Cluster Charge Bad APVs", 30001, 1., 30000.); //done
  h1ClusterChargeBadAPVs_->SetXTitle("Cluster charge [adc]");
  h1ClusterChargeBadAPVs_->SetYTitle("Entries");
  h1ClusterChargeBadAPVs_->SetLineWidth(2);
  h1ClusterChargeBadAPVs_->SetLineStyle(2);
  
  h1ClusterSizeGoodAPVs_= fs_->make<TH1F>("ClusterSizeGoodAPV", "Cluster Size Good APVs", 128, 1, 128);  //done
  h1ClusterSizeGoodAPVs_->SetXTitle("Cluster Size [strips]");
  h1ClusterSizeGoodAPVs_->SetYTitle("Entries");
  h1ClusterSizeGoodAPVs_->SetLineWidth(2);
  h1ClusterSizeGoodAPVs_->SetLineStyle(2);
  
  h1ClusterChargeGoodAPVs_= fs_->make<TH1F>("ClusterChargeGoodAPV", "Cluster Charge Good APVs", 30001, 1., 30000.); //done
  h1ClusterChargeGoodAPVs_->SetXTitle("Cluster charge [adc]");
  h1ClusterChargeGoodAPVs_->SetYTitle("Entries");
  h1ClusterChargeGoodAPVs_->SetLineWidth(2);
  h1ClusterChargeGoodAPVs_->SetLineStyle(2);
  
  h2BadAPVperEventVsCent_ = fs_->make<TH2F>("BadAPVperEventVsCent","BadAPVperEventVsCent", 40, -0.5, 39.5,  50001, -0.5, 50000.5); //done
  h2BadAPVperEventVsCent_->SetXTitle(" Centrality");
  h2BadAPVperEventVsCent_->SetYTitle("# Bad Modules/ev.");
  
  h2BadAPVVsEtaVsR_= fs_->make<TH2F>("BadAPVVsEtaVsR","Bad APV Vs Eta Vs R",1000, -5, 5,  1300, 0, 130);;
  h2BadAPVVsEtaVsR_->SetXTitle("Eta");
  h2BadAPVVsEtaVsR_->SetYTitle("R [cm]");
 
  
}
DEFINE_FWK_MODULE(SiStripBadAPVAnalyzer);

