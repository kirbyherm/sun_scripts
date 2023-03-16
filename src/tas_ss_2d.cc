// Kirby Hermansen 
// 2022-10-19
// Purpose: To create a histogram of SuN's TAS

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>
#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

extern const ConfigParams p;
int ss_bin_width = p.Get_SS_Bin_Width();
int ss_lower_limit = p.Get_SS_Window()[0];
int ss_upper_limit = p.Get_SS_Window()[1];
int bin_width = p.Get_TAS_Bin_Width();
int tas_lower_limit = p.Get_TAS_Window()[0];
int tas_upper_limit = p.Get_TAS_Window()[1];
std::string rootfiles_path = p.Get_Rootfiles_Path();

TH2F* get_histogram_tas_ss(std::string isotope = "ti56", int time_window = 10000, std::string decay_type = "true", bool isDaughter=false){

    std::string str_time_window = std::to_string(time_window);
    if (p.Use_Daughter() && isDaughter)
        isotope = p.Get_Isotope_Daughter();
    else
        isotope = p.Get_Isotope();
   
    std::string dssd_window = std::to_string(p.Get_DSSD_Window());
    TChain* chain = new TChain("t");
    for (int i=2052; i<2169; i++)
    {
        std::string run_str = std::to_string(i);
        std::string filepath = rootfiles_path+str_time_window+"ms/true/Run" + run_str + "_SuN_" + str_time_window + "ms_" + isotope + "_" + dssd_window + "px_bound.root";
        if (decay_type == "rand")
            filepath = rootfiles_path+str_time_window+"ms/bg/BWRun" + run_str + "_SuN_" + str_time_window + "ms_" + isotope  + "_" + dssd_window + "px_bound.root";
        const char * c = filepath.c_str();
        chain->Add(c);
    } 
    std::string cuts = p.Get_Cuts(isDaughter);
    chain->Draw(Form("SuN_TAS_high_gain_cal:SuN_Seg_high_gain_cal>>%s(%d,%d,%d,%d,%d,%d)",decay_type.c_str(),int((ss_upper_limit - ss_lower_limit)/ss_bin_width),ss_lower_limit,ss_upper_limit,int((tas_upper_limit - tas_lower_limit)/bin_width),tas_lower_limit,tas_upper_limit)
                ,Form("%s",cuts.c_str()),"");
    TH2F* tas = ((TH2F*) gDirectory->Get(Form("%s",decay_type.c_str())));
    return tas;
}


void single_tas_ss_rand_true(std::string isotope = "v56", int time_window = 1000, TFile * fOut = new TFile("tas.root","RECREATE"), bool isDaughter=false){

    std::string str_time_window = std::to_string(time_window);
    TH2F* bg = get_histogram_tas_ss(isotope, time_window, "rand", isDaughter);
    TH2F* tas = get_histogram_tas_ss(isotope, time_window, "true", isDaughter);
    tas->Draw();
    tas->SetName(Form("%s_tas_ss",isotope.c_str()));
    tas->Write();
    bg->Draw("same");
    bg->SetName(Form("%s_tas_ss_bg",isotope.c_str()));
    bg->Write();
    tas->Add(bg,-1);
    tas->SetStats(0);
    tas->Draw("hist");
    tas->SetName(Form("%s_tas_ss_bgs",isotope.c_str()));
    tas->Write();
    return;
}


void tas_ss_parent(){

    int time_window = get_correlator_max_time_window();
    TFile *fOut = new TFile("output_tas_ss_parent.root","RECREATE");
    std::string isotope = p.Get_Isotope();
    single_tas_ss_rand_true(isotope, time_window, fOut);
    fOut->Close();
    std::cout << p.Get_Cuts(false) << std::endl;
    return;

}

void tas_ss_daughter(){

    int time_window = get_correlator_max_time_window();
    TFile *fOut = new TFile("output_tas_ss_daughter.root","RECREATE");
    std::string isotope = p.Get_Isotope_Daughter();
    single_tas_ss_rand_true(isotope, time_window, fOut, true);
    fOut->Close();
    std::cout << p.Get_Cuts(true) << std::endl;
    return;

}

int main(){
    tas_ss_parent();
    tas_ss_daughter();
    return 0;
}
