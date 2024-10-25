

#include "include/InitTree.h"
#include "include/coordinateTools.h"
//#include "include/high_vn.h"













void Gosling::Loop(){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);
    //double R=0.8;
    //construct N-subjettiness

    //define Nsubjettiness functions for beta=1.0 using WTA KT axes
    double beta1=1.0;
    
    fastjet::contrib::Nsubjettiness       nSub1_beta1(1  ,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::UnnormalizedMeasure(beta1) );
    fastjet::contrib::Nsubjettiness       nSub2_beta1(2  ,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::UnnormalizedMeasure(beta1) );
    fastjet::contrib::Nsubjettiness       nSub3_beta1(3  ,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::UnnormalizedMeasure(beta1) ); 
    fastjet::contrib::NsubjettinessRatio nSub21_beta1(2,1,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::UnnormalizedMeasure(beta1) );
    fastjet::contrib::NsubjettinessRatio nSub32_beta1(3,2,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::UnnormalizedMeasure(beta1) );
    
    /*
    fastjet::contrib::Nsubjettiness       nSub1_beta1(1  ,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::NormalizedMeasure(beta1,R) );
    fastjet::contrib::Nsubjettiness       nSub2_beta1(2  ,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::NormalizedMeasure(beta1,R) );
    fastjet::contrib::Nsubjettiness       nSub3_beta1(3  ,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::NormalizedMeasure(beta1,R) ); 
    fastjet::contrib::NsubjettinessRatio nSub21_beta1(2,1,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::NormalizedMeasure(beta1,R) );
    fastjet::contrib::NsubjettinessRatio nSub32_beta1(3,2,fastjet::contrib::WTA_KT_Axes(),fastjet::contrib::NormalizedMeasure(beta1,R) );  
    */                                                                        
    double beta2=2.0;
    
     // Define Nsubjettiness functions for beta = 2.0, using Half-KT axes

    fastjet::contrib::Nsubjettiness         nSub1_beta2(1,   fastjet::contrib::HalfKT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta2));
    fastjet::contrib::Nsubjettiness         nSub2_beta2(2,   fastjet::contrib::HalfKT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta2));
    fastjet::contrib::Nsubjettiness         nSub3_beta2(3,   fastjet::contrib::HalfKT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta2));
    fastjet::contrib::NsubjettinessRatio   nSub21_beta2(2,1, fastjet::contrib::HalfKT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta2));
    fastjet::contrib::NsubjettinessRatio   nSub32_beta2(3,2, fastjet::contrib::HalfKT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta2));
    /*
    fastjet::contrib::Nsubjettiness         nSub1_beta2(1,   fastjet::contrib::HalfKT_Axes(), fastjet::contrib::NormalizedMeasure(beta2,R));
    fastjet::contrib::Nsubjettiness         nSub2_beta2(2,   fastjet::contrib::HalfKT_Axes(), fastjet::contrib::NormalizedMeasure(beta2,R));
    fastjet::contrib::Nsubjettiness         nSub3_beta2(3,   fastjet::contrib::HalfKT_Axes(), fastjet::contrib::NormalizedMeasure(beta2,R));
    fastjet::contrib::NsubjettinessRatio   nSub21_beta2(2,1, fastjet::contrib::HalfKT_Axes(), fastjet::contrib::NormalizedMeasure(beta2,R));
    fastjet::contrib::NsubjettinessRatio   nSub32_beta2(3,2, fastjet::contrib::HalfKT_Axes(), fastjet::contrib::NormalizedMeasure(beta2,R));
    */

    //define jet algorithm 
    
    
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm,R);


    double jetEtaCut=1.6;
    double jetPtCut=550;


    //
    TH1D *h_Sub1_beta1[trackbin];
    TH1D *h_Sub2_beta1[trackbin];
    TH1D *h_Sub3_beta1[trackbin];
    TH1D *h_Sub21_beta1[trackbin];
    TH1D *h_Sub32_beta1[trackbin];


    TH1D *h_Sub1_beta2[trackbin];
    TH1D *h_Sub2_beta2[trackbin];
    TH1D *h_Sub3_beta2[trackbin];
    TH1D *h_Sub21_beta2[trackbin];
    TH1D *h_Sub32_beta2[trackbin];


    //only for double check
    TH1D *h_Jet_pt_stored=new TH1D("h_Jet_pt_stored","h_Jet_pt_stored",200,400,3000);
    TH1D *h_Jet_pt_claculated=new TH1D("h_Jet_pt_calculaetd","h_Jet_pt_calculated",200,400,3000);

    TH1D *h_Jet_eta_stored=new TH1D("h_Jet_eta_stored","h_Jet_eta_stored",100,-2,2);
    TH1D *h_Jet_eta_claculated=new TH1D("h_Jet_eta_calculaetd","h_Jet_eta_calculated",100,-2,2);



    for(int i=0;i<trackbin;i++){
        h_Sub1_beta1[i]=new TH1D(Form("h_Sub1_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub1_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,200);
        h_Sub2_beta1[i]=new TH1D(Form("h_Sub2_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub2_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,200);
        h_Sub3_beta1[i]=new TH1D(Form("h_Sub3_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub3_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,200);
        h_Sub21_beta1[i]=new TH1D(Form("h_Sub21_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub21_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,1);
        h_Sub32_beta1[i]=new TH1D(Form("h_Sub32_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub32_beta1_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,1);


        h_Sub1_beta2[i]=new TH1D(Form("h_Sub1_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub1_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,200);
        h_Sub2_beta2[i]=new TH1D(Form("h_Sub2_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub2_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,200);
        h_Sub3_beta2[i]=new TH1D(Form("h_Sub3_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub3_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,200);
        h_Sub21_beta2[i]=new TH1D(Form("h_Sub21_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub21_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,1);
        h_Sub32_beta2[i]=new TH1D(Form("h_Sub32_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),Form("h_Sub32_beta2_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i] ),100,0,1);
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

            std::vector< fastjet::PseudoJet > particles;
            //*******************************ENTER JET LOOP*****************************
            for(int ijet=0;ijet<JetCounter;ijet++){
                for(int idau=0;idau<(*genDau_pid)[ijet].size();idau++){
                    TLorentzVector v;
                    int tmp_pid=(*genDau_pid)[ijet][idau];
                    v.SetPtEtaPhiM((*genDau_pt)[ijet][idau],(*genDau_eta)[ijet][idau],(*genDau_phi)[ijet][idau],particleData.m0(tmp_pid));
                    fastjet::PseudoJet particle(v.Px(),v.Py(),v.Pz(),v.E() );
                    particle.set_user_index(tmp_pid);
                    particles.emplace_back(particle);
                }
            }//*******************************END JET LOOP**************************


            for(int ijet=0;ijet<JetCounter;ijet++){
                if(   fabs( (*genJetEta)[ijet] ) > jetEtaCut   ) continue;
                if(         (*genJetPt)[ijet]    < jetPtCut    ) continue;

                h_Jet_pt_stored->Fill(  (*genJetPt)[ijet] );
                h_Jet_eta_stored->Fill( (*genJetEta)[ijet]  );
            }


            //recluster 
            fastjet::ClusterSequence cs(particles,jet_def);
            std::vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
            
            //*****************************************ENTER JET LOOP****************************
            for(int ijet=0;ijet<jets.size();ijet++){
                if(   fabs( jets[ijet].eta() ) > jetEtaCut   ) continue;
                if(         jets[ijet].pt()    < jetPtCut    ) continue;

                double  tau1_beta1 =  nSub1_beta1(jets[ijet]);
                double  tau2_beta1 =  nSub2_beta1(jets[ijet]);
                double  tau3_beta1 =  nSub3_beta1(jets[ijet]);
                double tau21_beta1 = nSub21_beta1(jets[ijet]);
                double tau32_beta1 = nSub32_beta1(jets[ijet]);


                double  tau1_beta2 =  nSub1_beta2(jets[ijet]);
                double  tau2_beta2 =  nSub2_beta2(jets[ijet]);
                double  tau3_beta2 =  nSub3_beta2(jets[ijet]);
                double tau21_beta2 = nSub21_beta2(jets[ijet]);
                double tau32_beta2 = nSub32_beta2(jets[ijet]);


                h_Jet_pt_claculated->Fill( jets[ijet].pt()  );
                h_Jet_eta_claculated->Fill(jets[ijet].eta() );

                std::vector<fastjet::PseudoJet> constituents=jets[ijet].constituents();
                int Nch=0;
                for(int icon=0;icon<constituents.size();icon++){
                    if(  particleData.charge( constituents[icon].user_index() )   == 0   ) continue;//charge
                    if(    fabs(  constituents[icon].pt()     )                   < 0.3  ) continue;//lab pt
                    if(    fabs(  constituents[icon].eta()    )                   > 2.4  ) continue;//lab eta
                    Nch++;
                }

               for(int itrk=0;itrk<trackbin;itrk++){
                    if( Nch>=trackbinbounds[itrk] && Nch<trackbinboundsUpper[itrk] ){
                        h_Sub1_beta1[itrk]->Fill(tau1_beta1);
                        h_Sub2_beta1[itrk]->Fill(tau2_beta1);
                        h_Sub3_beta1[itrk]->Fill(tau3_beta1);
                        h_Sub21_beta1[itrk]->Fill(tau21_beta1);
                        h_Sub32_beta1[itrk]->Fill(tau32_beta1);


                        h_Sub1_beta2[itrk]->Fill(tau1_beta2);
                        h_Sub2_beta2[itrk]->Fill(tau2_beta2);
                        h_Sub3_beta2[itrk]->Fill(tau3_beta2);
                        h_Sub21_beta2[itrk]->Fill(tau21_beta2);
                        h_Sub32_beta2[itrk]->Fill(tau32_beta2);
                    }
                }



               


            }//*****************************************END JET LOOP********************************



        }//********************************END EVENT LOOP ******************************


        fFile->Close();




    }//***********************************END FILE LOOP*********************************


    TFile *fout=new TFile("test.root","RECREATE");
    for(int i=0;i<trackbin;i++){
        h_Sub1_beta1[i]->Write();
        h_Sub2_beta1[i]->Write();
        h_Sub3_beta1[i]->Write();
        h_Sub21_beta1[i]->Write();
        h_Sub32_beta1[i]->Write();


        h_Sub1_beta2[i]->Write();
        h_Sub2_beta2[i]->Write();
        h_Sub3_beta2[i]->Write();
        h_Sub21_beta2[i]->Write();
        h_Sub32_beta2[i]->Write();
       


    }    
    h_Jet_pt_claculated->Write();
    h_Jet_eta_claculated->Write();
    h_Jet_pt_stored->Write();
    h_Jet_eta_stored->Write();
    fout->Close();


                        
            
}


//Code enters execution here
int main(int argc, const char* argv[])
{   

    
    std::vector<std::string> listOfFiles;
    for(int i=1;i<131+1;i++){
        listOfFiles.push_back(Form("/home/fl52/ppJetAnalysis/Pythia_CP5_pp/Pythia_CP5_pp_rootfile/pp_highMultGen_nChGT60_%d.root ",i) );
    }
    
    Gosling m = Gosling(listOfFiles);
    m.Loop();
    

    

    return 0;
}

