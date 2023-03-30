// Kirby Hermansen 
// 2022-10-19
// Purpose: build a set of configs to be used in all the histograms

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TH1F.h"
#include "TSystem.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

template <typename T>
T readJSON(boost::property_tree::ptree loadPtreeRoot, std::string var="foo"){
    return loadPtreeRoot.get_child(var).get_value<T>();
}

/*********************************************
****CONFIG PARAMS FOR BOTH ROOT AND PYTHON****
**********************************************/
class ConfigParams{

public:

    void LoadParams(std::string filename);
    ConfigParams(std::string filename){
        LoadParams(filename);
    }
    std::string Get_Cuts(bool, int) const;
    std::string Get_simCuts(int) const;
    std::vector<std::string> Get_Levels() const;
    std::string Get_Isotope() const;
    std::string Get_Rootfiles_Path() const;
    std::string Get_Isotope_Daughter() const;
    std::string Get_Isotope_Granddaughter() const;
    std::string Get_Isotope_Neutron_Daughter() const;
    std::string Get_Isotope_Neutron_Granddaughter() const;
    double Get_Pn() const;
    double Get_Daughter_HL() const;
    double Get_Granddaughter_HL() const;
    double Get_Neutron_Daughter_HL() const;
    double Get_Neutron_Granddaughter_HL() const;
    bool Use_Daughter() const;
    bool Use_Middle_Segments_Only() const;
    int Get_DSSD_Window() const;
    int Get_SS_Bin_Width() const;
    int Get_Decay_Bin_Width() const;
    int Get_TAS_Bin_Width() const;
    double Get_Scale_Factor() const;
    std::vector<int> Get_TAS_Window() const;
    std::vector<int> Get_SS_Window() const;
    std::vector<int> Get_simTAS_Window() const;
    std::vector<int> Get_simSS_Window() const;
    std::vector<int> Get_Parent_Decay_Window() const;
    std::vector<int> Get_Daughter_Decay_Window() const;

private:

    std::string rootfiles_path;
    std::string isotope;
    std::string isotope_daughter;
    std::string isotope_gdaughter;
    std::string isotope_ndaughter;
    std::string isotope_ngdaughter;
    bool use_daughter;
    bool middle_segments_only;
    int dssd_window;
    int ss_bin_width;
    int tas_bin_width;
    int decay_bin_width;
    std::vector<int> tas_window;
    std::vector<int> ss_window;
    bool cut_on_time;
    int time_start;
    int time_end;
    int daughter_time_start;
    int daughter_time_end;
    bool cut_on_veto;
    int veto_start;
    int veto_end;
    bool cut_on_tas;
    int tas_start;
    int tas_end;
    bool cut_on_ss;
    bool not_single;
    int ss_start;
    int ss_end;
    bool cut_on_mult;
    int mult_start;
    int mult_end;
    double pn;
    double scale_factor;
    double daughterhl;
    double gdaughterhl;
    double ndaughterhl;
    double ngdaughterhl;
    std::vector<std::string> levels;


};

extern const ConfigParams p("config.json");

void ConfigParams::LoadParams(std::string filename="config.json"){

    namespace pt = boost::property_tree;
    pt::ptree loadPtreeRoot;
    pt::read_json(filename.c_str(),loadPtreeRoot);
    rootfiles_path = readJSON<std::string>(loadPtreeRoot,"rootfiles_path"); 
    isotope = readJSON<std::string>(loadPtreeRoot,"isotope"); 
    isotope_daughter = readJSON<std::string>(loadPtreeRoot,"isotope_daughter"); 
    isotope_gdaughter = readJSON<std::string>(loadPtreeRoot,"isotope_gdaughter"); 
    isotope_ndaughter = readJSON<std::string>(loadPtreeRoot,"isotope_ndaughter"); 
    isotope_ngdaughter = readJSON<std::string>(loadPtreeRoot,"isotope_ngdaughter"); 
    use_daughter = readJSON<bool>(loadPtreeRoot,"use_daughter"); 
    middle_segments_only = readJSON<bool>(loadPtreeRoot,"middle_segments_only"); 
    tas_bin_width = readJSON<int>(loadPtreeRoot,"tas_bin_width"); 
    tas_window.push_back(readJSON<int>(loadPtreeRoot,"tas_lower_limit")); 
    tas_window.push_back(readJSON<int>(loadPtreeRoot,"tas_upper_limit")); 
    dssd_window = readJSON<int>(loadPtreeRoot,"dssd_window"); 
    ss_bin_width = readJSON<int>(loadPtreeRoot,"ss_bin_width"); 
    ss_window.push_back(readJSON<int>(loadPtreeRoot,"ss_lower_limit")); 
    ss_window.push_back(readJSON<int>(loadPtreeRoot,"ss_upper_limit")); 

    cut_on_time = readJSON<bool>(loadPtreeRoot,"cut_on_time");
    time_start = readJSON<int>(loadPtreeRoot,"time_start");
    time_end = readJSON<int>(loadPtreeRoot,"time_end");
    daughter_time_start = readJSON<int>(loadPtreeRoot,"daughter_time_start");
    daughter_time_end = readJSON<int>(loadPtreeRoot,"daughter_time_end");
    cut_on_veto = readJSON<bool>(loadPtreeRoot,"cut_on_veto");
    veto_start = readJSON<int>(loadPtreeRoot,"veto_start");
    veto_end = readJSON<int>(loadPtreeRoot,"veto_end");
    cut_on_tas = readJSON<bool>(loadPtreeRoot,"cut_on_tas");
    tas_start = readJSON<int>(loadPtreeRoot,"tas_start");
    tas_end = readJSON<int>(loadPtreeRoot,"tas_end");
    cut_on_ss = readJSON<bool>(loadPtreeRoot,"cut_on_ss");
    not_single = readJSON<bool>(loadPtreeRoot,"not_single");
    ss_start = readJSON<int>(loadPtreeRoot,"ss_start");
    ss_end = readJSON<int>(loadPtreeRoot,"ss_end");
    cut_on_mult = readJSON<bool>(loadPtreeRoot,"cut_on_mult");
    mult_start = readJSON<int>(loadPtreeRoot,"mult_start");
    mult_end = readJSON<int>(loadPtreeRoot,"mult_end");
    
    scale_factor = readJSON<double>(loadPtreeRoot,"scale_factor");
    pn = readJSON<double>(loadPtreeRoot,"pn");

    daughterhl = readJSON<double>(loadPtreeRoot,"daughterhl");
    gdaughterhl = readJSON<double>(loadPtreeRoot,"gdaughterhl");
    ndaughterhl = readJSON<double>(loadPtreeRoot,"ndaughterhl");
    ngdaughterhl = readJSON<double>(loadPtreeRoot,"ngdaughterhl");

    decay_bin_width = readJSON<int>(loadPtreeRoot,"decay_bin_width"); 

    for (pt::ptree::value_type &level : loadPtreeRoot.get_child("levels"))
    {
        levels.push_back(level.second.data());
    }

}

std::string ConfigParams::Get_Rootfiles_Path() const{
    return rootfiles_path;
}

std::string ConfigParams::Get_Isotope() const{
    return isotope;
}

std::string ConfigParams::Get_Isotope_Daughter() const{
    return isotope_daughter;
}

std::string ConfigParams::Get_Isotope_Granddaughter() const{
    return isotope_gdaughter;
}

std::string ConfigParams::Get_Isotope_Neutron_Daughter() const{
    return isotope_ndaughter;
}

std::string ConfigParams::Get_Isotope_Neutron_Granddaughter() const{
    return isotope_ngdaughter;
}

int ConfigParams::Get_Decay_Bin_Width() const{
    return decay_bin_width;
}

int ConfigParams::Get_TAS_Bin_Width() const{
    return tas_bin_width;
}

int ConfigParams::Get_SS_Bin_Width() const{
    return ss_bin_width;
}

int ConfigParams::Get_DSSD_Window() const{
    return dssd_window;
}

std::vector<int> ConfigParams::Get_TAS_Window() const{
    return tas_window; 
}

std::vector<int> ConfigParams::Get_simTAS_Window() const{
    return {tas_start, tas_end}; 
}

std::vector<int> ConfigParams::Get_simSS_Window() const{
    return {ss_start, ss_end}; 
}

std::vector<int> ConfigParams::Get_SS_Window() const{
    return ss_window; 
}

std::vector<int> ConfigParams::Get_Parent_Decay_Window() const{
    return {int(time_start/1e6), int(time_end/1e6)}; 
}

std::vector<int> ConfigParams::Get_Daughter_Decay_Window() const{
    return {int(daughter_time_start/1e6), int(daughter_time_end/1e6)}; 
}

std::vector<std::string> ConfigParams::Get_Levels() const{
    return levels; 
}

std::string ConfigParams::Get_simCuts(int seg=-1) const{

    std::string cuts = "(eneDSSD>0)";
    if(cut_on_veto){
        if (cuts.size() != 0)
            cuts += "&&";
        cuts += Form("(eneVeto<%d)",veto_end);
    }

    if(cut_on_ss){
        std::string seg_cut="";
        if (cuts.size() != 0)
            cuts += "&&";
        seg_cut += "(";
        if (not_single){
        for (int i_seg=0; i_seg<8; i_seg++){
            if (i_seg == seg) 
                continue;
            if (seg_cut.size() > 1 && i_seg < 8)
                seg_cut += "&&";
            seg_cut += Form("(SuN_Seg_high_gain_cal[%d]<%d||SuN_Seg_high_gain_cal[%d]>%d)",i_seg,i_seg,ss_start,ss_end);
            }
        } else {
            for (int i_seg=0; i_seg<8; i_seg++){
                if (i_seg == seg) 
                    continue;
                if (seg_cut.size() > 1 && i_seg < 8)
                    seg_cut += "||";
                seg_cut += Form("(SuN_Seg_high_gain_cal[%d]>%d&&SuN_Seg_high_gain_cal[%d]<%d)",i_seg,i_seg,ss_start,ss_end);
            }
        }
        seg_cut += ")";
        cuts += seg_cut;
    }

    if(cut_on_tas){
        if (cuts.size() != 0)
            cuts += "&&";
        cuts += Form("(eneAll>=%d&&eneAll<%d)",tas_start,tas_end);
    }

    if(cut_on_mult){
        if (cuts.size() != 0)
            cuts += "&&";
        cuts += Form("(multi>%d&&multi<%d)",mult_start,mult_end);
    }

    return cuts;
}
std::string ConfigParams::Get_Cuts(bool isDaughter=false, int seg=-1) const{

    std::string cuts = "";
    if(cut_on_veto){
        if (cuts.size() != 0)
            cuts += "&&";
        cuts += Form("(E_veto_raw<%d)",veto_end);
    }

    if(cut_on_time){
        if (cuts.size() != 0)
            cuts += "&&";
        std::vector<int> corr_times;
        if (isDaughter){
            corr_times.push_back(daughter_time_start);
            corr_times.push_back(daughter_time_end);
        } else {
            corr_times.push_back(time_start);
            corr_times.push_back(time_end);
        }
        cuts += Form("(correlation_implant_decay_time>%d&&correlation_implant_decay_time<%d)",corr_times[0],corr_times[1]);
    }

    if(cut_on_ss){
        std::string seg_cut="";
        if (cuts.size() != 0)
            cuts += "&&";
        seg_cut += "(";
        if (not_single){
        for (int i_seg=0; i_seg<8; i_seg++){
            if (i_seg == seg) 
                continue;
            if (use_daughter && (i_seg==0 || i_seg==3 || i_seg==4 || i_seg==7))
                continue;
            if (seg_cut.size() > 1 && i_seg < 8)
                seg_cut += "&&";
            seg_cut += Form("(SuN_Seg_high_gain_cal[%d]<%d||SuN_Seg_high_gain_cal[%d]>%d)",i_seg,ss_start,i_seg,ss_end);
            }
        } else {
            for (int i_seg=0; i_seg<8; i_seg++){
                if (i_seg == seg) 
                    continue;
                if (use_daughter && (i_seg==0 || i_seg==3 || i_seg==4 || i_seg==7))
                    continue;
                if (seg_cut.size() > 1 && i_seg < 8)
                    seg_cut += "||";
                seg_cut += Form("(SuN_Seg_high_gain_cal[%d]>%d&&SuN_Seg_high_gain_cal[%d]<%d)",i_seg,ss_start,i_seg,ss_end);
            }
        }
        seg_cut += ")";
        cuts += seg_cut;
    }

    if(cut_on_tas){
        if (cuts.size() != 0)
            cuts += "&&";
        cuts += Form("(SuN_TAS_high_gain_cal>%d&&SuN_TAS_high_gain_cal<%d)",tas_start,tas_end);
    }

    if(cut_on_mult){
        if (cuts.size() != 0)
            cuts += "&&";
        cuts += Form("(SuN_m_high_gain_cal>%d&&SuN_m_high_gain_cal<%d)",mult_start,mult_end);
    }

    return cuts;
}

double ConfigParams::Get_Pn() const{
    return pn;
}

double ConfigParams::Get_Daughter_HL() const{
    return daughterhl;
}

double ConfigParams::Get_Granddaughter_HL() const{
    return gdaughterhl;
}

double ConfigParams::Get_Neutron_Daughter_HL() const{
    return ndaughterhl;
}

double ConfigParams::Get_Neutron_Granddaughter_HL() const{
    return ngdaughterhl;
}

double ConfigParams::Get_Scale_Factor() const{
    return scale_factor;

}

bool ConfigParams::Use_Daughter() const{
    return use_daughter;
}

bool ConfigParams::Use_Middle_Segments_Only() const{
    return middle_segments_only;
}

/********************************************
****FUNCTIONS AND VARIABLES ONLY FOR ROOT****
*********************************************/
void save_to_png( TCanvas *c1, std::string name ){

    TImage *img = TImage::Create();
    img->FromPad(c1);
    img->WriteImage(Form("%s.png",name.c_str()));
    return;
}

void add_bg_error(TH1F* tas, TH1F* bg, double scale=1.0){
//    tas->Add(bg,-scale);
    for (int i = 1; i < tas->GetNbinsX()+1; i++){
        if (tas->GetBinContent(i) ==0)
            tas->SetBinContent(i,1);
        if (bg->GetBinContent(i) ==0)
            bg->SetBinContent(i,1);
        tas->SetBinContent(i, tas->GetBinContent(i) - bg->GetBinContent(i)*scale);
//        tas->SetBinError(i, tas->GetBinContent(i)*sqrt(tas->GetBinError(i)/tas->GetBinContent(i)*tas->GetBinError(i)/tas->GetBinContent(i)+bg->GetBinError(i)/bg->GetBinContent(i)*bg->GetBinError(i)/bg->GetBinContent(i)));
        tas->SetBinError(i, sqrt(tas->GetBinContent(i)));
    }
}

int get_correlator_max_time_window(){

    int corr_ms = 2000; 
    return corr_ms;

}

