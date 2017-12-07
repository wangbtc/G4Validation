//*******************************************************************************
/* G4_103_Validation_macro.cc

This code serve the purpose of validation of new GEANT4 version. One should first simulate
desired events in BACCARAT and convert to root file. Subsequently using Bacc2AnalysisTree to
extract the results into analysis tree. This code then read in the output file of Bacc2AnalysisTree
and compare the energy spectrum between two different GEANT versions. -- Ryan Wang (12/2/2017)
****************************************************************************************************
chnge log
    12/2/2017 file created. --Ryan Wang

*/
//*******************************************************************************




#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TObjString.h"
#include "TKey.h"
#include "TCollection.h"

//#include "/home/vagrant/BACCARAT/tools/BaccRootConverterEvent.hh"
TFile* fin;
TFile* fin2;
TTree* tree;
TTree* tree2;
TString parName;
float pE,pX,pY,pZ,XeX,XeY,XeZ,XeER,XeNR,skinX,skinY,skinZ,skinER,skinNR;
float pE2,pX2,pY2,pZ2,XeX2,XeY2,XeZ2,XeER2,XeNR2,skinX2,skinY2,skinZ2,skinER2,skinNR2;
using namespace std;
void NormalizedSum(TH1F* hsum){

    double norm = hsum->Integral(1,hsum->GetNbinsX(),"width");
    hsum->Scale(1 / norm , "width");

}
bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
    infile.clear();
}
int dirExists(const char *path)
{
    struct stat info;

    if(stat( path, &info ) != 0)
        return 0;
    else if(info.st_mode & S_IFDIR)
        return 1;
    else
        return 0;
}
void GetTObjString(TFile *temf){
    TKey *key;
    TIter next(temf->GetListOfKeys());
    TObjString *instring;

    while (key = (TKey*)next()){
        instring = (TObjString*)key->ReadObj();
        cout<<" Primary particle is : "<<instring->GetString()<<endl;
    }
}
static void show_usage(string name){
    cout<<" Usage : ./Extract_info_analysis_tree [-co] file1 (file2)"<<name<<" Options:\n"
    <<" -c : compare two root files, if specified then enter the comparison mode.\n"
    <<" -o : Set the name of output file.\n"
    <<" -i : Read the input file.\n"
    <<" -wd : Working directory\n"
    <<" -debug : Get in the debugging mode.\n"
    <<" -h or --help : Show the usage\n"
    <<" Enjoy ! -Ryan Wang"<<endl;
}

Int_t main(int argc, char *argv[]){
    //Define control variable here
    bool compare = false ;// compare two root file
    bool debug_mode = false ;// Flag of debugging mode
    string infile1 ; // File name of first file
    string infile2 ; // File name of second file
    string ofilename ; //output filename
    string working_dir ;// working directory
    TString full_inputfile1 ;// full path to input file 1
    TString full_inputfile2 ;// full path to input file 2
    TString full_outfile ;// full path to out put file
    TObjString *ParticleNamee; // Name of parent particle
    TObjString *ParticleName; // Name of parent particle

    // Read in the info from command line
    if (argc<2){
        show_usage(argv[0]);
        return 1;
    }
    for (int i=1;i<argc;++i){
        string arg = argv[i];
        if ((arg=="-h") || (arg=="--help")){
            show_usage(argv[0]);
            return 0;
        }
        else if (arg=="-c"){
            compare = true;
        }
        else if (arg=="-wd"){
            working_dir = argv[i+1] ;
        }
        else if (arg=="-i") {
            if (!compare){
                infile1 = argv[i+1];
            }
            else {
                infile1 = argv[i+1];
                infile2 = argv[i+2];
            }
        }
        else if (arg=="-o"){
            ofilename = argv[i+1];
        }
        else if (arg=="-debug"){
            debug_mode = true;
        }
    }

    // Debugging code
    if (debug_mode){
        cout<<" The comparison mode is : "<<compare<<" file name is : ";
        if (!compare){
            cout<<infile1;
        }
        else {
            cout<<infile1<<" and "<<infile2;
        }
        cout<<" Output filename is : "<<ofilename<<endl;
    }

    // Set up file name and check if it is exist
    if (!compare){
        full_inputfile1 = working_dir + infile1 + TString(".root") ;
        full_outfile = working_dir + ofilename + TString(".root") ;
        if (!is_file_exist(full_inputfile1.Data())){
            cout<<" File is not exist ! "<<endl;
            return 1;
        }
    }
    else {
        full_inputfile1 = working_dir + infile1 + TString(".root") ;
        full_inputfile2 = working_dir + infile2 + TString(".root") ;
        if (!is_file_exist(full_inputfile1.Data()) || !is_file_exist(full_inputfile1.Data())){
            cout<<" Files"<<full_inputfile1<<" and "<<full_inputfile2<<" are not exist ! "<<endl;
            return 1;
        }
        full_outfile = working_dir + ofilename + TString(".root") ;
    }
    // Read root file

    TFile *fin = TFile::Open(full_inputfile1.Data(),"READ");
    if (fin == NULL){
        cout<<" File "<<full_inputfile1<<" is corrupted ! "<<endl;
        return 1;
    }
    else{
        tree = (TTree*) fin->Get("lzsim_analysis_tree");
        if (tree != NULL){
            cout<<" Successfully reading the tree ! "<<endl;
        }
        else{
            cout<<" No tree in this "<<full_inputfile1<<" file ! "<<endl;
        }
    }
    // Set Branch address for each parameter in analysis root file
    tree->SetBranchAddress("fPrimaryParE_keV",&pE) ;
    tree->SetBranchAddress("fPrimaryParX_cm",&pX) ;
    tree->SetBranchAddress("fPrimaryParY_cm",&pY) ;
    tree->SetBranchAddress("fPrimaryParZ_cm",&pZ) ;
    tree->SetBranchAddress("fLXeEDepER_keV",&XeER) ;
    tree->SetBranchAddress("fLXeX_cm",&XeX) ;
    tree->SetBranchAddress("fLXeY_cm",&XeY) ;
    tree->SetBranchAddress("fLXeZ_cm",&XeZ) ;
    tree->SetBranchAddress("fSkinEDepER_keV",&skinER) ;
    tree->SetBranchAddress("fSkinX_cm",&skinX) ;
    tree->SetBranchAddress("fSkinY_cm",&skinY) ;
    tree->SetBranchAddress("fSkinZ_cm",&skinZ) ;

    // Create some histogram
    TH1F* LXe_target_ER = new TH1F("GEANT4.10.3_target","LXe_target_ER",500,0,500);
    LXe_target_ER->SetLineColor(2);
    LXe_target_ER->SetXTitle(" keV ");
    LXe_target_ER->SetYTitle(" Normalized count ");
    LXe_target_ER->Sumw2();
    TH1F* LXe_skin_ER = new TH1F("GEANT4.10.3_Skin","LXe_skin_ER",500,0,500);
    LXe_skin_ER->SetLineColor(2);
    LXe_skin_ER->SetXTitle(" keV ");
    LXe_skin_ER->SetYTitle(" Normalized count ");
    LXe_skin_ER->Sumw2();
    TH1F* LXe_target_ER2 = new TH1F("GEANT4.9.5_target","LXe_target_ER2",500,0,500);
    LXe_target_ER2->SetLineColor(4);
    LXe_target_ER2->SetXTitle(" keV ");
    LXe_target_ER2->SetYTitle(" Normalized count ");
    LXe_target_ER2->Sumw2();
    TH1F* LXe_skin_ER2 = new TH1F("GEANT4.9.5_skin","LXe_skin_ER2",500,0,500);
    LXe_skin_ER2->SetLineColor(4);
    LXe_skin_ER2->SetXTitle(" keV ");
    LXe_skin_ER2->SetYTitle(" Normalized count ");
    LXe_skin_ER2->Sumw2();


    // Loop the ntuple to get the value
    for (int i=0;i<tree->GetEntries();i++){

        if (i==0){
            cout<<" Start Processing "<<full_inputfile1<<endl;
        }
        else if (i>0){
            cout<<" Finish processing "<<i<<" events. "<<endl;
        }
        //tree->Print();
        tree->GetEntry(i);

        LXe_target_ER->Fill(XeER);
        LXe_skin_ER->Fill(skinER);
    }


    if (compare){

        TFile *fin2 = TFile::Open(full_inputfile2.Data(),"READ");
        if (fin2 == NULL){
            cout<<" File is corrupted ! "<<endl;
                return 1;
        }
        else{
            tree2 = (TTree*) fin2->Get("lzsim_analysis_tree");
            if (tree2 != NULL){
                cout<<" Successfully reading the tree ! "<<endl;
            }
            else{
                cout<<" No tree in this "<<full_inputfile2<<" file ! "<<endl;
            }
        }
        //TBranch* temp2 = (TBranch*)tree2->GetBranch("cPrimaryParName");
        //temp2->SetAddress(&ParticleName2);
        tree2->SetBranchAddress("fPrimaryParE_keV",&pE2) ;
        tree2->SetBranchAddress("fPrimaryParX_cm",&pX2) ;
        tree2->SetBranchAddress("fPrimaryParY_cm",&pY2) ;
        tree2->SetBranchAddress("fPrimaryParZ_cm",&pZ2) ;
        tree2->SetBranchAddress("fLXeEDepER_keV",&XeER2) ;
        tree2->SetBranchAddress("fLXeX_cm",&XeX2) ;
        tree2->SetBranchAddress("fLXeY_cm",&XeY2) ;
        tree2->SetBranchAddress("fLXeZ_cm",&XeZ2) ;
        tree2->SetBranchAddress("fSkinEDepER_keV",&skinER2) ;
        tree2->SetBranchAddress("fSkinX_cm",&skinX2) ;
        tree2->SetBranchAddress("fSkinY_cm",&skinY2) ;
        tree2->SetBranchAddress("fSkinZ_cm",&skinZ2) ;
        //tree2->SetBranchAddress("cPrimaryParName",&ParticleNamee);

        for (int i=0;i<10000;i++){
            if (i==0){
                cout<<" Start Processing "<<full_inputfile2<<endl;
            }
            else if (i%100000==0){
                cout<<" Finish processing "<<i<<" events. "<<endl;
            }
            tree2->GetEntry(i);

            /*if (ParName[0].compare(ParticleNamee->String())==0){
                LXe_skin_ER2_Ac228->Fill(skinER2);
                cout<<" The Parname is : "<<ParticleNamee->String()<<" ParName is  : "<<ParName[0]<<" EDep in Skin is : "<<skinER2<<endl;
            }*/

            LXe_target_ER2->Fill(XeER2);
            LXe_skin_ER2->Fill(skinER2);
        }
    }

    //Normalize histogram w.r.t. Integral
    NormalizedSum(LXe_target_ER);
    NormalizedSum(LXe_skin_ER);

    if (compare){
        NormalizedSum(LXe_target_ER2);
        NormalizedSum(LXe_skin_ER2);
        NormalizedSum(LXe_skin_ER2_Ac228);
        NormalizedSum(LXe_skin_ER_Ac228);

    }

    TFile* outfile_root = new TFile(full_outfile.Data(),"RECREATE") ;
    TH1F* LXe_target_ER_clone = (TH1F*)LXe_target_ER->Clone("LXe_target_ER_clone");
    LXe_target_ER_clone->SetLineColor(kBlack);
    LXe_target_ER_clone->SetMinimum(0);  // Define Y ..
    LXe_target_ER_clone->SetMaximum(2); // .. range
    LXe_target_ER_clone->Sumw2();
    LXe_target_ER_clone->SetStats(0);      // No statistics on lower plot
    LXe_target_ER_clone->Divide(LXe_target_ER2);
    LXe_target_ER_clone->SetMarkerStyle(21);

    TH1F* LXe_skin_ER_clone = (TH1F*)LXe_skin_ER->Clone("LXe_skin_ER_clone");
    LXe_skin_ER_clone->SetLineColor(kBlack);
    LXe_skin_ER_clone->SetMinimum(0);  // Define Y ..
    LXe_skin_ER_clone->SetMaximum(2); // .. range
    LXe_skin_ER_clone->Sumw2();
    LXe_skin_ER_clone->SetStats(0);      // No statistics on lower plot
    LXe_skin_ER_clone->Divide(LXe_skin_ER2);
    LXe_skin_ER_clone->SetMarkerStyle(21);

    TCanvas* ER_com = new TCanvas("ER_com","ER_com");
    ER_com->SetLogy();
    LXe_target_ER->Draw();
    LXe_skin_ER->Draw("sames");
    ER_com->SaveAs("LXe_target_skin_ER_com.pdf");

    TCanvas* ER_com2 = new TCanvas("ER_com2","ER_com2");
    TPad *p1 = new TPad("p1","pad",0.0,0.25,1.0,1.0);
    TPad *p2 = new TPad("p2","pad",0.0,0.0,1.0,0.25);
    ER_com2->SetLogy();
    p1->Draw();
    p2->Draw();
    p1->cd();
    LXe_target_ER->Draw();
    LXe_target_ER2->Draw("sames");
    p2->cd();
    TLine *line = new TLine(0,1,500,1);
    LXe_target_ER_clone->Draw("p");
    line->SetLineColor(2);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw("");
    ER_com2->SaveAs("LXe_target_skin_ER_com2.pdf");

    TCanvas* ER_com3 = new TCanvas("ER_com3","ER_com3");
    ER_com3->SetLogy();
    TPad *p3 = new TPad("p3","pad",0.0,0.25,1.0,1.0);
    TPad *p4 = new TPad("p4","pad",0.0,0.0,1.0,0.25);
    p3->Draw();
    p4->Draw();
    p3->cd();
    LXe_skin_ER->Draw();
    LXe_skin_ER2->Draw("sames");
    p4->cd();
    TLine *line2 = new TLine(0,1,500,1);
    LXe_target_ER_clone->Draw();
    line2->SetLineColor(2);
    line2->SetLineStyle(2);
    line2->SetLineWidth(3);

    LXe_skin_ER_clone->Draw("p");
    line->Draw("");
    ER_com3->SaveAs("LXe_target_skin_ER_com3.pdf");


    ER_com->Write();
    ER_com2->Write();
    ER_com3->Write();
    LXe_skin_ER->Write();
    LXe_skin_ER2->Write();
    outfile_root->Write();
    outfile_root->Close();
    return 0 ;
}
