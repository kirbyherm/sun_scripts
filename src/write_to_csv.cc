#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

extern const ConfigParams p;

void write_mult(){

    std::string rootFile = Form("%s_mult_final.root",p.Get_Isotope().c_str());
    std::string histName = Form("%s_mult_final",p.Get_Isotope().c_str());
    std::string outFile = Form("%s_mult_obs.csv",p.Get_Isotope().c_str());
    TFile * f1 = TFile::Open(rootFile.c_str());
    TH1F * h1;
    h1 = (TH1F*)(f1->Get(histName.c_str()));
    std::ofstream myFile(outFile.c_str());
    for ( int i = 1; i<h1->GetNbinsX()+1; i++){
        int binout = TMath::Max((Double_t)(h1->GetBinContent(i)),(Double_t)0);
        int errorout = TMath::Max((Double_t)(h1->GetBinError(i)),(Double_t)0);
        myFile << i << "," << binout << "," << errorout << "\n";
    }
    myFile.close();
    return;
}

void write_ss(){

    std::string rootFile = Form("%s_ss_final.root",p.Get_Isotope().c_str());
    std::string histName = Form("%s_ss_final",p.Get_Isotope().c_str());
    std::string outFile = Form("%s_ss_obs.csv",p.Get_Isotope().c_str());
    TFile * f1 = TFile::Open(rootFile.c_str());
    TH1F * h1;
    h1 = (TH1F*)(f1->Get(histName.c_str()));
    std::ofstream myFile(outFile.c_str());
    
    for ( int i = 1; i<h1->GetNbinsX()+1; i++){
        int binout = TMath::Max((Double_t)(h1->GetBinContent(i)),(Double_t)0);
        int errorout = TMath::Max((Double_t)(h1->GetBinError(i)),(Double_t)0);
        myFile << i << "," << binout << "," << errorout << "\n";
    }
    myFile.close();
    return;
}
 
void write_tas(){

    std::cout << p.Get_Isotope() << std::endl;
    std::string rootFile = Form("%s_tas_final.root",p.Get_Isotope().c_str());
    std::string histName = Form("%s_tas_final",p.Get_Isotope().c_str());
    std::string outFile = Form("%s_tas_obs.csv",p.Get_Isotope().c_str());
    std::cout << rootFile << std::endl;
    TFile * f1 = TFile::Open(rootFile.c_str());
    TH1F * h1;
    h1 = (TH1F*)(f1->Get(histName.c_str()));
    std::ofstream myFile(outFile.c_str());
    for ( int i = 1; i<h1->GetNbinsX()+1; i++){
        int binout = TMath::Max((Double_t)(h1->GetBinContent(i)),(Double_t)0);
        int errorout = TMath::Max((Double_t)(h1->GetBinError(i)),(Double_t)0);
        myFile << i << "," << binout << "," << errorout << "\n";
        if (i > 15 && i < 20)
        std::cout << i << "," << binout << " : " << h1->GetBinContent(i) << "," << errorout  << " : " << h1->GetBinError(i)<< "\n";
    }
    myFile.close();
    
    return;
}

void write_2d(){

    std::string rootFile = Form("%s_tas_ss_final.root",p.Get_Isotope().c_str());
    std::string histNames[4] = {Form("%s_tas",p.Get_Isotope().c_str()),Form("%s_tas_bg",p.Get_Isotope().c_str()),Form("%s_dtas",p.Get_Isotope().c_str()),Form("%s_dtas_bg",p.Get_Isotope().c_str())};
    std::string outFile = Form("%s_2d_obs.csv",p.Get_Isotope().c_str());
    std::string outErrorFile = Form("%s_2d_error.csv",p.Get_Isotope().c_str());
    TFile * f1 = TFile::Open(rootFile.c_str());
    TH2F * hError;
    hError = (TH2F*)(f1->Get(histNames[0].c_str()));
    for (int i = 1; i < 4; i++){
        TH2F * h2;
        std::cout << "hErrorBinContent: " << hError->GetBinContent(33,13) << " hErrorErrorContent: " << hError->GetBinError(33,13) << std::endl;
        hError->Sumw2();
        std::cout << "hErrorBinContent: " << hError->GetBinContent(33,13) << " hErrorErrorContent: " << hError->GetBinError(33,13) << std::endl;
        h2 = (TH2F*)(f1->Get(histNames[i].c_str()));
        h2->Sumw2();
        std::cout << "h2BinContent: " << h2->GetBinContent(33,13) << " h2ErrorContent: " << h2->GetBinError(33,13) << std::endl;
        hError->Add(h2);
        std::cout << "hErrorBinContent: " << hError->GetBinContent(33,13) << " hErrorErrorContent: " << hError->GetBinError(33,13) << std::endl;
        std::cout << "h2BinContent: " << h2->GetBinContent(33,13) << " h2ErrorContent: " << h2->GetBinError(33,13) << std::endl;
        delete h2;
    }
 
    std::ofstream myErrorFile(outErrorFile.c_str());
//    h1->Rebin(2);
    for ( int i = 1; i<hError->GetNbinsY()+1; i++){
        for ( int j = 1; j<hError->GetNbinsX()+1; j++){
            int binout = (Double_t)(hError->GetBinContent(j,i));
            myErrorFile << binout; 
            if (j < hError->GetNbinsX())
                myErrorFile << " "; 
            else 
                myErrorFile << "\n";
        }
    }
    myErrorFile.close();


    TH2F * hData;
    std::string histName = Form("%s_tas_bgs_dbgs",p.Get_Isotope().c_str());
//    std::string histName = "sc57_tas_bgs";
    hData = (TH2F*)(f1->Get(histName.c_str()));
    std::ofstream myDataFile(outFile.c_str());
    for ( int i = 1; i<hData->GetNbinsY()+1; i++){
        for ( int j = 1; j<hData->GetNbinsX()+1; j++){
            int binout = (Double_t)(hData->GetBinContent(j,i));
            myDataFile << binout; 
            if (j < hData->GetNbinsX())
                myDataFile << " "; 
            else 
                myDataFile << "\n";
        }
    }
    myDataFile.close();
    
    return;

}

void write_to_csv(){
//    write_2d();
    write_tas();
    write_ss();
    write_mult();
    return;
}

int main(){
    write_to_csv();
    return 0;
}
