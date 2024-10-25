

#include "include/InitTree.h"
#include "include/coordinateTools.h"
//#include "include/high_vn.h"













void Gosling::Loop(){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    //define C/A algorithm , the frist step ( reclustering  ) of softdrop
    fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm ,fastjet::JetDefinition::max_allowable_R);
    
    //define softdrop
    std::vector< fastjet::contrib::SoftDrop > sd;
    for(int ibeta=0;ibeta<beta_bin;ibeta++){
        sd.push_back( fastjet::contrib::SoftDrop( beta_SD[ibeta],z_cut,R)  );
        sd[ibeta].set_reclustering(false,0);
    }

    //define histogram
    TH1D *h_SD_JetMass_dis[beta_bin][trackbin];
    TH1D *h_SD_Zg_dis[beta_bin][trackbin];
    for(int ibeta=0;ibeta<beta_bin;ibeta++){
        for(int itrk=0;itrk<trackbin;itrk++){
            h_SD_JetMass_dis[ibeta][itrk]=new TH1D(Form("h_SD_JetMass_dis_beta_%d_Nch_%d_%d",(int)beta_SD[ibeta],trackbinbounds[itrk],trackbinboundsUpper[itrk] ),Form("h_SD_JetMass_dis_beta_%d_Nch_%d_%d",(int)beta_SD[ibeta],trackbinbounds[itrk],trackbinboundsUpper[itrk] ),2000,0,2000 );
            h_SD_Zg_dis[ibeta][itrk]=new TH1D(Form("h_SD_Zg_dis_beta_%d_Nch_%d_%d",(int)beta_SD[ibeta],trackbinbounds[itrk],trackbinboundsUpper[itrk]),Form("h_SD_Zg_dis_beta_%d_Nch_%d_%d",(int)beta_SD[ibeta],trackbinbounds[itrk],trackbinboundsUpper[itrk]),50,0.0,0.5 );
        }
    }

   



   


    //***********************************ENTE FILE LOOP********************************
    for(int ifile=0;ifile<fileList.size();ifile++){
        fFile=TFile::Open(fileList[ifile].c_str(),"READ");
        std::cout<<"current file:"<<fileList[ifile]<<std::endl;
        std::cout<<(ifile+1)<<"/"<<fileList.size()<<std::endl;

        //open tree
        TTree *tree=(TTree*)fFile->Get("trackTree");
        Init(tree);
        Long64_t nentries = fChain->GetEntriesFast();
        std::cout<<"Total Entries ="<<nentries<<std::endl;
       

        //********************************ENETR EVENT LOOP******************************
        for(int ievent=0;ievent<nentries;ievent++){
            fChain->GetEntry(ievent);
            int JetCounter=(*genJetPt).size();

            
            //*******************************ENTER JET LOOP*****************************
            for(int ijet=0;ijet<JetCounter;ijet++){
                //cut on jet 
                if(   fabs( (*genJetEta)[ijet] ) > jetEtaCut   ) continue;
                if(         (*genJetPt)[ijet]    < jetPtCut    ) continue;

               
                //determine the Nch
                int Nch=0;
                for(int idau=0;idau<(*genDau_pid)[ijet].size();idau++){
                    if(  particleData.charge( (*genDau_pid)[ijet][idau] )   == 0   ) continue;//charge
                    if(    fabs(  (*genDau_pt)[ijet][idau]     )                   < 0.3  ) continue;//lab pt
                    if(    fabs(  (*genDau_eta)[ijet][idau]    )                   > 2.4  ) continue;//lab eta
                    Nch++;
                }

                //fill the daughter particles of this jet into particles 
                std::vector< fastjet::PseudoJet > particles;
                for(int idau=0;idau<(*genDau_pid)[ijet].size();idau++){
                    TLorentzVector v;
                    int tmp_pid=(*genDau_pid)[ijet][idau];
                    v.SetPtEtaPhiM((*genDau_pt)[ijet][idau],(*genDau_eta)[ijet][idau],(*genDau_phi)[ijet][idau],particleData.m0(tmp_pid));
                    fastjet::PseudoJet particle(v.Px(),v.Py(),v.Pz(),v.E() );
                    particle.set_user_index(tmp_pid);
                    particles.emplace_back(particle);
                }

                //apply C/A algorithm 
                fastjet::ClusterSequence cs(particles,jet_def);
                vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
                //jets should contain only one jet 




                for(int itrk=0;itrk<trackbin;itrk++){
                    if(Nch>=trackbinbounds[itrk] && Nch<trackbinboundsUpper[itrk]){
                        for(int ibeta=0;ibeta<beta_bin;ibeta++){
                            fastjet::PseudoJet sd_jet=sd[ibeta]( jets[0] );
                            double temp_sdJetMass=sd_jet.m();
                            double temp_zg=sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();

                            h_SD_Zg_dis[ibeta][itrk]->Fill(temp_zg);
                            h_SD_JetMass_dis[ibeta][itrk]->Fill(temp_sdJetMass);



                        }
                    }



                }
            }//**********************END JET LOOP********************************
               



            



        }//********************************END EVENT LOOP ******************************


        fFile->Close();




    }//***********************************END FILE LOOP*********************************


    TFile *fout=new TFile("HighNch.root","RECREATE");
    for(int ibeta=0;ibeta<beta_bin;ibeta++){
        for(int itrk=0;itrk<trackbin;itrk++){
            h_SD_Zg_dis[ibeta][itrk]->Write();
            h_SD_JetMass_dis[ibeta][itrk]->Write();
        }
    }

    fout->Close();
   

                        
            
}


//Code enters execution here
int main(int argc, const char* argv[])
{   

    
    std::vector<std::string> listOfFiles;
    for(int i=1;i<131+1;i++){
        listOfFiles.push_back(Form("/home/fl52/ppJetAnalysis/Pythia_CP5_pp/Pythia_CP5_pp_rootfile/pp_highMultGen_nChGT60_%d.root",i) );
    }
    
    Gosling m = Gosling(listOfFiles);
    m.Loop();
    

    

    return 0;
}

