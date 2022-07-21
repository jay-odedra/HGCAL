#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <list>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <typeinfo>
#include <TGraph.h>
#include <TCanvas.h>
#include <TNtuple.h> 
#include <TH1.h>
#include <string>
#include <TString.h>
#include <unordered_map>
#include <array>
#include <optional>
#include <THStack.h>
#include <TLegend.h>

float gethistvalue(std::string a,std::string b ,std::string c ){
    typedef std::unordered_map<std::string, float> histproperties;
    typedef std::unordered_map<std::string ,histproperties> variablename;
    std::unordered_map<std::string,variablename> treeleveldict {
        {"tracksters",
        {{"raw_energy",{{"binl", -0.1},{"binr",50},{"nbins",250}}},
        {"raw_em_energy",{{"binl", -0.1},{"binr",40},{"nbins",250}}},
        {"barycenter_x",{{"binl", -220},{"binr",220},{"nbins",250}}},
        {"barycenter_y",{{"binl", -220},{"binr",220},{"nbins",250}}},
        {"barycenter_z",{{"binl", -550},{"binr",550},{"nbins",250}}},
        {"trackster_barycenter_eta",{{"binl", -4.5},{"binr",4.5},{"nbins",250}}},
        {"trackster_barycenter_phi",{{"binl", -4},{"binr",4},{"nbins",250}}},
        {"eVector0_x",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"eVector0_y",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"eVector0_z",{{"binl", 1.2},{"binr",1.2},{"nbins",250}}},
        {"id_probabilities",{{"binl", -0.2},{"binr",1.2},{"nbins",250}}}}
        },
        {"candidates",
        {{"candidate_px",{{"binl", -20},{"binr",20},{"nbins",250}}},
        {"candidate_py",{{"binl", -20},{"binr",20},{"nbins",250}}},
        {"candidate_pz",{{"binl", -40},{"binr",40},{"nbins",250}}},
        {"candidate_charge",{{"binl", -2},{"binr",2},{"nbins",2500}}},
        {"candidate_id_probabilities",{{"binl", -0.2},{"binr",1.2},{"nbins",250}}},
        {"tracksters_in_candidate",{{"binl", 0},{"binr",80},{"nbins",80}}}}
        },
        {"trackstersMerged",
        {{"barycenter_x",{{"binl", -220},{"binr",220},{"nbins",250}}},
        {"barycenter_y",{{"binl", -220},{"binr",220},{"nbins",250}}},
        {"barycenter_z",{{"binl", -550},{"binr",550},{"nbins",250}}},
        {"eVector0_x",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"eVector0_y",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"eVector0_z",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"id_probabilities",{{"binl", -0.2},{"binr",1.2},{"nbins",250}}}}
        },
        {"associations",
        {{"tsCLUE3D_recoToSim_SC",{{"binl", -1},{"binr",120},{"nbins",250}}},
        {"tsCLUE3D_recoToSim_SC_score",{{"binl", -0.2},{"binr",1.2},{"nbins",250}}}}
        },
        {"simtrackstersSC",
        {{"stsSC_raw_energy",{{"binl", -0.1},{"binr",50},{"nbins",250}}},
        {"stsSC_raw_em_energy",{{"binl", -0.1},{"binr",50},{"nbins",250}}},
        {"stsSC_barycenter_x",{{"binl", -220},{"binr",220},{"nbins",250}}},
        {"stsSC_barycenter_y",{{"binl", -220},{"binr",220},{"nbins",250}}},
        {"stsSC_barycenter_z",{{"binl", -550},{"binr",550},{"nbins",250}}},
        {"stsSC_trackster_barycenter_eta",{{"binl", -4.5},{"binr",4.5},{"nbins",250}}},
        {"stsSC_trackster_barycenter_phi",{{"binl", -4},{"binr",4},{"nbins",250}}},
        {"stsSC_eVector0_x",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"stsSC_eVector0_y",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"stsSC_eVector0_z",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"stsSC_id_probabilities",{{"binl", -0.2},{"binr",1.2},{"nbins",250}}}}
        },
        {"simtrackstersCP",
        {{"stsCP_raw_energy",{{"binl", -0.1},{"binr",700},{"nbins",250}}},
        {"stsCP_raw_em_energy",{{"binl", -0.1},{"binr",700},{"nbins",250}}},
        {"stsCP_barycenter_x",{{"binl", -220},{"binr",220},{"nbins",250}}},
        {"stsCP_barycenter_y",{{"binl", -220},{"binr",220},{"nbins",250}}},
        {"stsCP_barycenter_z",{{"binl", -550},{"binr",550},{"nbins",250}}},
        {"stsCP_trackster_barycenter_eta",{{"binl", -4.5},{"binr",4.5},{"nbins",250}}},
        {"stsCP_trackster_barycenter_phi",{{"binl", -4},{"binr",4},{"nbins",250}}},
        {"stsCP_eVector0_x",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"stsCP_eVector0_y",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"stsCP_eVector0_z",{{"binl", -1.2},{"binr",1.2},{"nbins",250}}},
        {"stsCP_id_probabilities",{{"binl", -0.2},{"binr",1.2},{"nbins",250}}}}
        },
        };
    if(b =="SIZE" && c == "SIZE"){
        return (float)treeleveldict[a].size();
    }
    else if(c =="SIZE"){
        return (float)treeleveldict[a][b].size();
    }
    else {
        return treeleveldict[a][b][c];
    }
}

std::vector<std::string> getvarnamevalue(std::string name){
  std::map<std::string, std::vector<std::string>>  dict {
    {"tracksters",{"raw_energy","raw_em_energy","barycenter_x","barycenter_y","barycenter_z","trackster_barycenter_eta",
    "trackster_barycenter_phi","eVector0_x","eVector0_y","eVector0_z","id_probabilities"}},

    {"candidates",{"candidate_px","candidate_py","candidate_pz","candidate_charge","candidate_id_probabilities","tracksters_in_candidate"}},

    {"trackstersMerged",{"barycenter_x","barycenter_y","barycenter_z","eVector0_x","eVector0_y","eVector0_z","id_probabilities"}},

    {"associations",{"tsCLUE3D_recoToSim_SC","tsCLUE3D_recoToSim_SC_score"}},

    {"simtrackstersSC",{"stsSC_raw_energy","stsSC_raw_em_energy","stsSC_barycenter_x","stsSC_barycenter_y","stsSC_barycenter_z","stsSC_trackster_barycenter_eta",
    "stsSC_trackster_barycenter_phi","stsSC_eVector0_x","stsSC_eVector0_y","stsSC_eVector0_z","stsSC_id_probabilities"}},

    {"simtrackstersCP",{"stsCP_raw_energy","stsCP_raw_em_energy","stsCP_barycenter_x","stsCP_barycenter_y","stsCP_barycenter_z","stsCP_trackster_barycenter_eta",
    "stsCP_trackster_barycenter_phi","stsCP_eVector0_x","stsCP_eVector0_y","stsCP_eVector0_z","stsCP_id_probabilities"}},
    };

    return dict[name];
}









void hgcalntupleplotternotfancy(){
    TString dir("SinglePionPhotonElectronCompNoPU"); //foldername
    TString filetag("combinedPionPhotonElectronNoPU");
    TFile *pionfile = TFile::Open("D:/Downloads/50filespionNopileup.root"); //file1
    TFile *electronfile = TFile::Open("D:/Downloads/50fileselectronNopileup.root"); //file2
    TFile * photonfile = TFile::Open("D:/Downloads/50photonfilesNoPu.root");
    std::array< TString ,6> treenames = {"tracksters","candidates","trackstersMerged","associations","simtrackstersSC","simtrackstersCP"};
    //std::cout<<gethistvalue("tracksters","raw_energy","SIZE");
    for(int i=0; i<treenames.size(); i++){
        TCanvas c (  "can1"  ,  "Canvas t i t l e" , 1200 , 900 );
        
        TTree *pionfiletree = (TTree*)pionfile->Get("ticlNtuplizer/"+TString(treenames[i]));
        TTree *electronfiletree = (TTree*)electronfile->Get("ticlNtuplizer/"+TString(treenames[i]));
        TTree *photonfiletree = (TTree*)photonfile->Get("ticlNtuplizer/"+TString(treenames[i]));
        //std::cout<<gethistvalue((std::string)treenames[i],"SIZE","SIZE");


        for(int j=0; j<gethistvalue((std::string)treenames[i],"SIZE","SIZE");j++){
            std::vector<std::string> variablenamesvec = getvarnamevalue((std::string)treenames[i]);
            TString varnametstring(variablenamesvec[j]);
            TH1D *h1 = new TH1D("h1",filetag+"_"+treenames[i]+"_"+varnametstring,gethistvalue((std::string)treenames[i],(std::string)varnametstring,"nbins"),gethistvalue((std::string)treenames[i],(std::string)varnametstring,"binl"),gethistvalue((std::string)treenames[i],(std::string)varnametstring,"binr"));
            TH1D *h2 = new TH1D("h2",filetag+"_"+treenames[i]+"_"+varnametstring,gethistvalue((std::string)treenames[i],(std::string)varnametstring,"nbins"),gethistvalue((std::string)treenames[i],(std::string)varnametstring,"binl"),gethistvalue((std::string)treenames[i],(std::string)varnametstring,"binr"));
            TH1D *h3 = new TH1D("h3",filetag+"_"+treenames[i]+"_"+varnametstring,gethistvalue((std::string)treenames[i],(std::string)varnametstring,"nbins"),gethistvalue((std::string)treenames[i],(std::string)varnametstring,"binl"),gethistvalue((std::string)treenames[i],(std::string)varnametstring,"binr"));

            h1->SetLineColor(2);
            h2->SetLineColor(4);
            h3->SetLineColor(8);
            h1->SetStats(FALSE);
            double factor = 1.;
            //h1->Scale(factor/(h1->Integral()),"nosw2");
            //h2->Scale(factor/(h2->Integral()),"nosw2");
            
            pionfiletree->Draw(varnametstring +">>h1");
            TH1* htemp1 = (TH1*)gPad->GetPrimitive("h1");
            htemp1->Scale((factor/htemp1->Integral()),"nosw2");
            
            //htemp1->GetYaxis()->SetRange(0,binmax*3.2);
            electronfiletree->Draw(varnametstring +">>h2","","same");
            TH1* htemp2 = (TH1*)gPad->GetPrimitive("h2");
            htemp2->Scale((factor/htemp2->Integral()),"nosw2");

            photonfiletree->Draw(varnametstring +">>h3","","same");
            TH1* htemp3 = (TH1*)gPad->GetPrimitive("h3");
            htemp3->Scale((factor/htemp3->Integral()),"nosw2");




            double binmax2 = htemp2->GetMaximum();
            double binmax = htemp1->GetMaximum();
            double binmax3 = htemp3->GetMaximum();
            std::cout<<binmax<<endl;
            std::cout<<binmax2<<endl;
            std::cout<<binmax3<<endl;
            if(binmax > binmax2 && binmax3){
                //std::cout<<"binmax1 wins"<<endl;
                htemp1->GetYaxis()->SetRangeUser(0.0001, binmax*1.9);
                htemp2->GetYaxis()->SetRangeUser(0.0001, binmax*1.9);
                htemp3->GetYaxis()->SetRangeUser(0.0001, binmax*1.9);

            }else if(binmax2 > binmax && binmax3) {
                htemp1->GetYaxis()->SetRangeUser(0.0001, binmax2*1.3);
                htemp2->GetYaxis()->SetRangeUser(0.0001, binmax2*1.3);
                htemp3->GetYaxis()->SetRangeUser(0.0001, binmax2*1.3);
            }else {
                htemp1->GetYaxis()->SetRangeUser(0.0001, binmax3*1.3);
                htemp2->GetYaxis()->SetRangeUser(0.0001, binmax3*1.3);
                htemp3->GetYaxis()->SetRangeUser(0.0001, binmax3*1.3);
                
            }



            TLegend *legend = new TLegend(0.15, 0.8, 0.25, 0.9);
            legend->SetMargin(1.1);
            legend->SetTextFont(61);
            legend->AddEntry(h1,"Pion", "l");
            legend->SetBorderSize(0);
            legend->SetTextSize(0.04);
            legend->SetFillStyle(0);
            legend->AddEntry(h2,"Electron", "l");
            legend->AddEntry(h3,"Photon", "l");

            legend->Draw();
            c.SetLogy();

            c.SaveAs(dir+"/"+dir+"_"+treenames[i]+"_"+varnametstring+"_LogPlots.pdf");
            delete h1;
            delete h3;
            delete h2;
        }
    }
    
}