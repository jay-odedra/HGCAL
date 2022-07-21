#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TMath.h>
#include <TTree.h>
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TLeaf.h>
#include <list>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <typeinfo>
#include <TGraph.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TNtuple.h> 
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include "ROOT/RVec.hxx"
#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include <iostream>
#include <typeinfo>
#include <algorithm>
#include <typeinfo>
#include <iostream>
#include <fstream>

class matchedevent {
    public:
        float particleidindex;
        unsigned short matchedtrackster;
        float score; 
        int indexes[3];

};

float getIndex(std::vector<float> v, float K){
    auto it = find(v.begin(), v.end(), K);
    // If element was found
    if (it != v.end()){
        int index = it - v.begin();
        return index;
    }
    else{
        std::cout<<"ERROR NO INDEX FOUND"<<endl;
        return 600000;
    }
}

std::vector<matchedevent> minimizedvector(std::vector<matchedevent> original, std::set<float> uniq){
    std::vector<matchedevent> copyd=original;
    std::set<int> indexestoremove;
    std::set<float>::iterator matchedeventval;
    std::vector<int> indexestoremovevec;
    for( matchedeventval = uniq.begin(); matchedeventval != uniq.end(); matchedeventval++){
        bool flag= false;
        float minscoreval= 200;
        int previndex = -1;
        for(int j=0;j<original.size();j++){
            if(original[j].matchedtrackster==(*matchedeventval)){
                if(original[j].score<minscoreval){
                    minscoreval=original[j].score;
                    std::cout<<"minscoreval  ===="<<minscoreval<<"  index  =="<<j<<" , "<<endl;
                    if(flag){
                        indexestoremove.insert(previndex);
                    }
                    flag = true;
                    previndex = j;

                }
                else{
                    std::cout<<"index to remov ======  "<<j<<endl;
                    indexestoremove.insert(j);
                }
            }
            else{
                continue;
            }

            //std::sort(indexestoremove.rbegin(), indexestoremove.rend());

        }
    }
    indexestoremovevec.insert(indexestoremovevec.end(), indexestoremove.begin(), indexestoremove.end());
    sort(indexestoremovevec.begin(),indexestoremovevec.end(),std::greater<int>());
    for( int k:indexestoremovevec){
        copyd.erase(copyd.begin() + k);
    }

    return copyd;
}











void Confusionmatrixofficialclean(){
    //gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
    int counternumberofelectrons = 0;
    int confmatrix[8][8] = {};
    std::string line;
    std::ifstream filelist;
    filelist.open("electrontest.txt");
    if(!filelist.is_open()) {
      perror("Error open");
      exit(EXIT_FAILURE);
    }
    while(getline(filelist, line)) {
        TString filename(line);
        std::cout<<line<<endl;
    //TFile *infile = TFile::Open("/eos/user/j/jodedra/HGCALHACKATHON/photonsample/50photonfilesNoPu.root", "READ"); //open file
        TFile *infile = TFile::Open(filename, "READ"); //open file

        //TTree * s2reco=0; //sim to reco matched event
        //TTree * s2recoscore=0; //sim to reco matched event score
        //TTree * PIDSCtrackster=0; //SIM PID
        //TTree * trackstermergePID=0;

        TTreeReader S2RReader;
        TTreeReader S2RScoreReader;
        TTreeReader PIDReader;
        TTreeReader trackstermergedPIDReader;



        TTreeReaderArray<std::vector<unsigned short>> stsSCidprob2Dvec ={S2RReader,"Mergetstracksters_simToReco_SC"};
        TTreeReaderArray<std::vector<float>> ssts2dscorevec ={S2RScoreReader,"Mergetstracksters_simToReco_SC_score"};
        TTreeReaderArray<std::vector<float>> probvec2d {PIDReader,"stsSC_id_probabilities"};
        TTreeReaderArray<std::vector<float>> probmergedvec2d {trackstermergedPIDReader,"id_probabilities"};

        TTree *trackstermergePID = (TTree*)infile->Get("ticlNtuplizer/trackstersMerged");
        trackstermergedPIDReader.SetTree(trackstermergePID);

        TTree *s2reco = (TTree*)infile->Get("ticlNtuplizer/associations");
        S2RReader.SetTree(s2reco);

        TTree *s2recoscore = (TTree*)infile->Get("ticlNtuplizer/associations");
        S2RScoreReader.SetTree(s2recoscore);

        TTree *PIDSCtrackster = (TTree*)infile->Get("ticlNtuplizer/simtrackstersSC");
        PIDReader.SetTree(PIDSCtrackster);
        int entriesnum = s2reco->GetEntries();
        std::cout<<entriesnum<<endl;
        int counternumbpredele =0;
        int countermatchedele =0;
        for(int kjl=0;kjl<entriesnum;kjl++){

            int evtnum=kjl;
            
            S2RReader.SetLocalEntry(evtnum);
            S2RScoreReader.SetLocalEntry(evtnum);
            PIDReader.SetLocalEntry(evtnum);
            trackstermergedPIDReader.SetLocalEntry(evtnum);

            std::vector<matchedevent> matchedeventvector;
            std::set<float> uniquetracksters;
            std:cout<<"simTracksterSC EVENT NUMBER ======================================================" <<evtnum<<endl;

            //for loop over the number of maps from event x to merged tracksters
            for( int j =0; j<stsSCidprob2Dvec.GetSize();j++){ 
                matchedevent loopobject;
                for (int z =0;z<stsSCidprob2Dvec[j].size();z++){
                loopobject.particleidindex = getIndex((std::vector<float>)probvec2d[j],1); //gives index of particle id 1 ie index 0 = photon index 1 = electron
                float matchedeventtrackster = stsSCidprob2Dvec[j][z]; 
                loopobject.matchedtrackster = matchedeventtrackster;


                if( (float)ssts2dscorevec[j][z]<1.0 && (float)ssts2dscorevec[j][z]>=0.0){
                        std::cout<<"duringloop scoreval  "<< ssts2dscorevec[j][z]<<"  trackmatched  ="<<matchedeventtrackster<<" particle ID == "<<loopobject.particleidindex<<endl;
                        uniquetracksters.insert(matchedeventtrackster); //puts matchedevents to set
                        loopobject.score = ssts2dscorevec[j][z];
                        loopobject.indexes[0] = evtnum;
                        loopobject.indexes[1] = j;
                        loopobject.indexes[2] = z;
                        if(matchedeventtrackster<60000 && matchedeventtrackster>-1){
                            matchedeventvector.push_back(loopobject);//pushes matchedevent object to vector to sort
                        }
                    }else{
                        continue;
                    }
                }

            }
            std::cout<<"score ==== "<< matchedeventvector.size()<<endl;

            std::vector<matchedevent> minvector = minimizedvector(matchedeventvector,uniquetracksters);
            std::cout<<" final minimized vector size =="<< minvector.size()<<endl;
            for( int v =0;v<minvector.size();v++){
                int jindex = getIndex((std::vector<float>)probmergedvec2d[minvector[v].matchedtrackster],1);
                int iindex = minvector[v].particleidindex;
                confmatrix[iindex][jindex] = confmatrix[iindex][jindex] + 1;
                std::cout<<"simEvent ==" <<evtnum<<"matchedtracksetermerge event ==" <<minvector[v].matchedtrackster<<endl;
                std::cout<<" PID OF SIMEVENT =="<<iindex<<" PID OF MATCHEDTRCKSTERMERG == "<<jindex<<endl;
            }
                //std::cout<<minvector[0].matchedtrackster;
                //std::cout<<"size of probs ========" <<probmergedvec2d.GetSize()<<endl;
                //std::cout<<"PID value matched=====" <<getIndex((std::vector<float>)probmergedvec2d[minvector[0].matchedtrackster],1)<<endl;
                //std::cout<<" simtrackster PID VALUE MATCHED ====== "<< minvector[0].particleidindex<<endl;
        }
        //infile->Close();
        //delete infile;
    }
    
    std::ofstream myfile ("outputs/MultiParticlesNoPUMATRIX.txt");
    if(myfile.is_open()){
        for(int i =0;i<8;i++){
            for(int j =0;j<8;j++){
                myfile<<confmatrix[i][j]<<"  ";
            }
        myfile<<"\n";
    }
    }
    myfile.close();
    
    

    //getIndex((std::vector<float>)probmergedvec2d[j],1)
}