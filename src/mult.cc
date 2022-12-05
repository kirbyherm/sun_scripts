// Kirby Hermansen 
// 2022-10-19
// Purpose: To create a histogram of SuN's event multiplicities

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include <iostream>
#include <string>
#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

extern const ConfigParams p;
std::string rootfiles_path = p.Get_Rootfiles_Path();

TH1F* get_histogram_mult(std::string isotope = "ti56", int time_window = 10000, std::string decay_type = "true", bool isDaughter=false){

    std::string str_time_window = std::to_string(time_window);
    if (p.Use_Daughter() && isDaughter)
        isotope = p.Get_Isotope_Daughter();
    else
        isotope = p.Get_Isotope();
   
    TChain* chain = new TChain("t");
    for (int i=2052; i<2169; i++)
    {
        std::string run_str = std::to_string(i);
        std::string filepath = rootfiles_path+str_time_window+"ms/true/Run" + run_str + "_SuN_" + str_time_window + "ms_" + isotope + "_9px_bound.root";
        if (decay_type == "rand")
            filepath = rootfiles_path+str_time_window+"ms/bg/BWRun" + run_str + "_SuN_" + str_time_window + "ms_" + isotope  + "_9px_bound.root";
        const char * c = filepath.c_str();
        chain->Add(c);
    } 
    std::string cuts = p.Get_Cuts(isDaughter);
    chain->Draw(Form("SuN_m_high_gain_cal>>%s(20,0,20)",decay_type.c_str())
               ,Form("%s",cuts.c_str())
               ,"");
    TH1F* mult = ((TH1F*) gDirectory->Get(Form("%s",decay_type.c_str())));
    return mult;
}


void single_mult_rand_true(std::string isotope = "v56", int time_window = 1000, TFile * fOut = new TFile("mult.root","RECREATE"), bool isDaughter=false){

    std::string str_time_window = std::to_string(time_window);
    TH1F* bg = get_histogram_mult(isotope, time_window, "rand", isDaughter);
    TH1F* mult = get_histogram_mult(isotope, time_window, "true", isDaughter);
    mult->Draw();
    mult->SetName(Form("%s_mult",isotope.c_str()));
    mult->Write();
    bg->Draw("same");
    bg->SetName(Form("%s_mult_bg",isotope.c_str()));
    bg->Write();
//    add_bg_error(mult, bg);
    mult->Add( bg,-1);
    mult->Draw();
    mult->SetName(Form("%s_mult_bgs",isotope.c_str()));
    mult->Write();
    return;
}


void mult_parent(){

    int time_window = get_correlator_max_time_window();
    TFile *fOut = new TFile("output_mult_parent.root","RECREATE");
    std::string isotope = p.Get_Isotope();
    single_mult_rand_true(isotope, time_window, fOut);
    fOut->Close();
    return;

}

void mult_daughter(){

    int time_window = get_correlator_max_time_window();
    TFile *fOut = new TFile("output_mult_daughter.root","RECREATE");
    std::string isotope = p.Get_Isotope_Daughter();
    single_mult_rand_true(isotope, time_window, fOut, true);
    fOut->Close();
    return;

}

int main(){
    mult_parent();
    mult_daughter();
    return 0;
}
