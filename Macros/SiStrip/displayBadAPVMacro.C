#include "TH1F.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TFile.h"
#include "TPaveStats.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TF1.h"
#include "TLegend.h"
#include "TKey.h"

#include "iostream"
#include "vector"
#include "math.h"
#include "map"

class displayBadAPVMacro {

public:
  TFile *f;
  TString BaseDir;
  TString dir[3];
  TString fullPath, title, subDet, genSubDet;
  TCanvas *C;
  void loop(){
    
    int NEntries =0; 
   // f->cd();
    f->cd(dir[0]);
       
    TIter nextkey(gDirectory->GetListOfKeys());
    TKey *key ;
    while ((key = (TKey*)nextkey()) && NEntries<1000) {
      NEntries++;
      TObject *obj = key->ReadObj();
      std::cout << " object " << obj->GetName() << std::endl;

      if ( obj->IsA()->InheritsFrom( "TH1" ) ) {

     	C->Clear();
    	
		TH1F* h_raw = (TH1F*)key->ReadObj();
        h_raw->Draw();

    	
		TH1F* h_bas = (TH1F*) f->Get(dir[1] +"/"+ obj->GetName());
		if(h_bas!=0){
      		h_bas->SetLineWidth(2);
	  		h_bas->SetLineStyle(2);
	  		h_bas->SetLineColor(2);
	  		h_bas->Draw("same");
		}
    	

    	/*   	
    	TH1F* h_clus = (TH1F*) f->Get(dir[2] +"/"+ obj->GetName());
		if(h_clus!=0){
       		h_clus->SetLineWidth(2);
	  		h_clus->SetLineStyle(2);
	  		h_clus->SetLineColor(4);
	  		h_clus->Draw("same");
		}
       */
		C->Update();
		C->SaveAs(TString("img/")+obj->GetName()+TString(".gif"));
      }
    }

  };

  displayBadAPVMacro(TString file){
    C = new TCanvas();
  
    f = new TFile(file);
  
    BaseDir="SiStripBadAPVAnalyzer/";
    dir[0]=BaseDir+"ProcessedRawDigis";
    dir[1]=BaseDir+"Baseline";
    dir[2]=BaseDir+"Clusters";

    f->cd(BaseDir);
    loop();
  }
  
  ~displayBadAPVMacro(){};
};
