// Kirby Hermansen 
// 2022-10-19
// Purpose: To create a histogram of the sum of SuN's individual segments

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
int bin_width = p.Get_SS_Bin_Width();
int ss_lower_limit = p.Get_SS_Window()[0];
int ss_upper_limit = p.Get_SS_Window()[1];
std::string rootfiles_path = p.Get_Rootfiles_Path();

TH1F* get_histogram_ss(std::string isotope = "ti56", int time_window = 10000, std::string decay_type = "true", bool isDaughter=false, int seg=0){

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
    chain->Draw(Form("SuN_Seg_high_gain_cal[%d]>>%s%d(%d,%d,%d)",seg,decay_type.c_str(),seg,int((ss_upper_limit-ss_lower_limit)/bin_width),ss_lower_limit,ss_upper_limit)
               ,Form("%s",cuts.c_str())
               ,"");
    TH1F* ss = ((TH1F*) gDirectory->Get(Form("%s%d",decay_type.c_str(),seg)));
    delete chain;
    return ss;
}

void write_hists(TH1F* ss, TH1F* bg, std::string isotope = "v56", std::string filename="test.root"){

    TFile *fOut = new TFile(Form("%s",filename.c_str()),"RECREATE");
    ss->SetName(Form("%s_ss",isotope.c_str()));
    ss->Write();
    bg->SetName(Form("%s_ss_bg",isotope.c_str()));
    bg->Write();
//    add_bg_error(ss, bg, -1);
    ss->Add(bg, -1);
    ss->SetName(Form("%s_ss_bgs",isotope.c_str()));
    ss->Write();
    fOut->Close();
    return;
}

void single_ss_rand_true(std::string isotope = "v56", int time_window = 1000, std::string filename="test.root", bool isDaughter=false){

    std::string str_time_window = std::to_string(time_window);
    TFile *fOut = new TFile(Form("%s",filename.c_str()),"RECREATE");
    if (p.Use_Middle_Segments_Only()){
        TH1F* bg = get_histogram_ss(isotope, time_window, "rand",isDaughter,1);
        TH1F* ss = get_histogram_ss(isotope, time_window, "true",isDaughter,1);
        ss->Write();
        for (int i=2; i<7; i++){
            if( i==3 || i==4)
                continue;
            bg->Add(get_histogram_ss(isotope, time_window, "rand",isDaughter,i));
            ss->Add(get_histogram_ss(isotope, time_window, "true",isDaughter,i));
        }
        write_hists(ss, bg, isotope, filename);
    } else {     
        TH1F* bg = get_histogram_ss(isotope, time_window, "rand",isDaughter,0);
        TH1F* ss = get_histogram_ss(isotope, time_window, "true",isDaughter,0);
        for (int i=1; i<8; i++){
            bg->Add(get_histogram_ss(isotope, time_window, "rand",isDaughter,i));
            ss->Add(get_histogram_ss(isotope, time_window, "true",isDaughter,i));
        }
        write_hists(ss, bg, isotope, filename);
    }
    return;
}


void ss_parent(){

    int time_window = get_correlator_max_time_window();
    std::string filename = "output_ss_parent.root";
    std::string isotope = p.Get_Isotope();
    single_ss_rand_true(isotope, time_window, filename);
    return;

}

void ss_daughter(){

    int time_window = get_correlator_max_time_window();
    std::string filename = "output_ss_daughter.root";
    std::string isotope = p.Get_Isotope_Daughter();
    single_ss_rand_true(isotope, time_window, filename, true);
    return;

}

int main(){
    ss_parent();
    ss_daughter();
    return 0;
}
