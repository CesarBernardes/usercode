/**
  USAGE: 
    DATA-only:    root.exe -b -l -q TnP_ICHEP_Trigger_data_all.root "plotTrigger.cxx(\"data_all\")"
    MC data-like: root.exe -b -l -q TnP_ICHEP_Trigger_datalike_mc.root  "plotTrigger.cxx(\"datalike_mc\")"
    DATA+MC:      root.exe -b -l -q TnP_ICHEP_Trigger_data_all.root TnP_ICHEP_Trigger_datalike_mc.root  "plotTrigger.cxx(\"data_vs_mc\")"

  REQUIRES:
   1) mkdir -p plots_ichep/muonid/
   2) provide a suitable "tdrStyle.cc" macro or similar
      (by default, it's taken from ~/cpp/tdrstyle.cc;
       if you need one you might want to grab ~gpetrucc/cpp/tdrstyle.cc)
*/
#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"
TString prefix = "plots_ichep_dev/trigger/";
TString basedir = "histoTrigger";

TFile *ref = 0;

TCanvas *c1 = 0;
void plotTrigger_ICHEP(TString scenario) {
    prefix = prefix+scenario+"/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (gROOT->GetListOfFiles()->GetEntries() == 2) {
        ref = (TFile *) gROOT->GetListOfFiles()->At(1);
        ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();
    }
    plotTriggerData();
}

void plotTriggerData() {
    char *ids[2] = { "Mu3", "L1DoubleMuOpen" };
    for (size_t i = 0; i < 2; ++i) {
        TString idname(ids[i]);
        plotTriggerData(idname);

    }

    for (size_t i = 0; i < 1; ++i) {
        TString idname = TString(ids[i])+"overL1";
        plotTriggerData(idname);
    }
    
}
void plotTriggerData(TString idname) {
    TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
    if (ref != 0) {
        TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_eta-0.0-0.8",  "pt_PLOT_abseta_bin0_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_eta-0.8-1.2",  "pt_PLOT_abseta_bin1_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_eta-1.2-1.6",  "pt_PLOT_abseta_bin2_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_eta-1.6-2.4",  "pt_PLOT_abseta_bin3_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_eta_pt-3-5",      "abseta_PLOT_pt_bin1_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_eta_pt-5-12",     "abseta_PLOT_pt_bin2_");
    } else {
        TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_mcTrue/");
        if (mc_pt_eta) {
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_eta-0.0-0.8", "pt_PLOT_abseta_bin0_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_eta-0.8-1.2", "pt_PLOT_abseta_bin1_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_eta-1.2-1.6", "pt_PLOT_abseta_bin2_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_eta-1.6-2.4", "pt_PLOT_abseta_bin3_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_eta_pt-3-5",     "abseta_PLOT_pt_bin1_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_eta_pt-5-12",    "abseta_PLOT_pt_bin2_");
        } else {
            single(fit_pt_eta, idname+"_pt_eta-0.0-0.8", "pt_PLOT_abseta_bin0_");
            single(fit_pt_eta, idname+"_pt_eta-0.8-1.2", "pt_PLOT_abseta_bin1_");
            single(fit_pt_eta, idname+"_pt_eta-1.2-1.6", "pt_PLOT_abseta_bin2_");
            single(fit_pt_eta, idname+"_pt_eta-1.6-2.4", "pt_PLOT_abseta_bin3_");
            single(fit_pt_eta, idname+"_eta_pt-3-5",     "abseta_PLOT_pt_bin1_");
            single(fit_pt_eta, idname+"_eta_pt-5-12",    "abseta_PLOT_pt_bin2_");
        }
    }

    doCanvas(fit_pt_eta, 4, 3, idname+"_eta_%d_pt_%d",   "abseta_bin%d__pt_bin%d_");
}


