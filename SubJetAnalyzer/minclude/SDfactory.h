#include "/home/fl52/fastjet/ClusterSequence.hh"
#include "/home/fl52/fastjet/contrib.055/RecursiveTools/SoftDrop.cc"
#include "1d2d_constants.h"
#include "Pythia8/Pythia.h"
using namespace fastjet;
Pythia8::Pythia pythia;
Pythia8::ParticleData &particleData = pythia.particleData;
#ifndef SDTree_h
#define SDTree_h



class MyUserInfo : public PseudoJet::UserInfoBase{
public:
   MyUserInfo(const int & pid_in) :
      _pid(pid_in){}

   int pid() const{return _pid;}
   


protected:
   int _pid;
   

};


class SDfactory{
public :
   TFile          *fFile;
   TTree          *fChain;
   long int            TotalEntry;
   std::vector<std::string> fileList;
   std::vector<std::string> SDfileList;

   // Fixed size dimensions of array or collections stored in the TTree if any.
   //trackTree
   vector<float>   *genJetEta;
   vector<float>   *genJetPt;
   vector<float>   *genJetPhi;
   vector<int>     *genJetChargedMultiplicity;
   vector<vector<int> > *genDau_chg;
   vector<vector<int> > *genDau_pid;
   vector<vector<float> > *genDau_pt;
   vector<vector<float> > *genDau_eta;
   vector<vector<float> > *genDau_phi;
   //SDtrackTree
   std::vector< float > *jetPt;
   std::vector< float > *jetEta;
   std::vector< float > *jetPhi;
   std::vector< int   > *jetChargedMultiplicity;
   std::vector< std::vector<int>  > *has_parents;
   std::vector< std::vector<float> > *jetZg;
   std::vector< std::vector<float> > *jetRg;
   std::vector< std::vector<float> > *subjet1_pt;
   std::vector< std::vector<float> > *subjet1_eta;
   std::vector< std::vector<float> > *subjet1_phi;
   std::vector< std::vector<float> > *subjet1_m;
   std::vector< std::vector<float> > *subjet2_pt;
   std::vector< std::vector<float> > *subjet2_eta;
   std::vector< std::vector<float> > *subjet2_phi;
   std::vector< std::vector<float> > *subjet2_m;
   std::vector< std::vector<float> > *dau_pt;
   std::vector< std::vector<int> >   *dau_chg;
   std::vector< std::vector<float> > *dau_eta;
   std::vector< std::vector<float> > *dau_phi;
   std::vector< std::vector<int> >   *dau_pid;
   std::vector< std::vector< std::vector<int> > * > dau_Flag;   // if flag=0 , groomed
                                                                // if flag=1 , this particle belongs to subjet1
                                                                // if flga=2 , this particle belongs to subjet2   

   TBranch        *b_genJetEta;  
   TBranch        *b_genJetPt;   
   TBranch        *b_genJetPhi;   
   TBranch        *b_genJetChargedMultiplicity;   
   TBranch        *b_genDau_chg;   
   TBranch        *b_genDau_pid;   
   TBranch        *b_genDau_pt;   
   TBranch        *b_genDau_eta;   
   TBranch        *b_genDau_phi; 



   SDfactory(std::vector<std::string> _fileList,std::vector<std::string>_SDfileList);
   virtual ~SDfactory();
   //virtual Int_t    GetEntry(Long64_t entry);
   virtual void     LoadData(TTree *tree);
   virtual void     GetMyEntry(int entry);
   virtual void     SDtreePlanter();
   
};	

SDfactory::SDfactory(std::vector<std::string> _fileList,std::vector<std::string>_SDfileList) 
   : fChain(0), fileList(_fileList) ,SDfileList(_SDfileList)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("total_output_ak8.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("total_output_ak8.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("total_output_ak8.root:/analyzer");
//      dir->GetObject("trackTree",tree);
//
//   }
   
  
}




SDfactory::~SDfactory()
{
   //if (!fChain) return;
   //delete fChain->GetCurrentFile();
}


void SDfactory::LoadData(TTree *tree){
   //set object pointer
   genJetEta=0;
   genJetPt=0;
   genJetPhi=0;
   genJetChargedMultiplicity=0;
   genDau_chg=0;
   genDau_pid=0;
   genDau_pt=0;
   genDau_eta=0;
   genDau_phi=0;
   //
   jetPt=new std::vector< float >;
   jetEta=new std::vector< float >;
   jetPhi=new std::vector< float >;
   jetChargedMultiplicity=new std::vector< int   >;
   has_parents=new std::vector< std::vector<int>   >;
   jetZg=new std::vector< std::vector<float> >;
   jetRg=new std::vector< std::vector<float> >;
   subjet1_pt=new std::vector< std::vector<float> >;
   subjet1_eta=new  std::vector< std::vector<float> >;
   subjet1_phi=new  std::vector< std::vector<float> >;
   subjet1_m=new  std::vector< std::vector<float> >;
   subjet2_pt=new  std::vector< std::vector<float> >;
   subjet2_eta=new  std::vector< std::vector<float> >;
   subjet2_phi=new  std::vector< std::vector<float> >;
   subjet2_m=new  std::vector< std::vector<float> >;
   dau_pt=new  std::vector< std::vector<float> >;
   dau_chg=new std::vector< std::vector<int> >;
   dau_eta=new  std::vector< std::vector<float> >;
   dau_phi=new  std::vector< std::vector<float> >;
   dau_pid=new std::vector< std::vector<int> >;

   for(int ibeta=0;ibeta<beta_bin;ibeta++){
      dau_Flag.push_back(new std::vector< std::vector<int> >);
   }

   if(! tree) {
      std::cout<<"Input tree does not exist!"<<std::endl;
      return;
   }

   fChain=tree;
   TotalEntry=fChain->GetEntriesFast();
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("genJetEta"                ,&genJetEta                ,&b_genJetEta                );
   fChain->SetBranchAddress("genJetPt"                 ,&genJetPt                 ,&b_genJetPt                 );
   fChain->SetBranchAddress("genJetPhi"                ,&genJetPhi                ,&b_genJetPhi                );
   fChain->SetBranchAddress("genJetChargedMultiplicity",&genJetChargedMultiplicity,&b_genJetChargedMultiplicity);
   fChain->SetBranchAddress("genDau_chg"               ,&genDau_chg               ,&b_genDau_chg               );
   fChain->SetBranchAddress("genDau_pid"               ,&genDau_pid               ,&b_genDau_pid               );
   fChain->SetBranchAddress("genDau_pt"                ,&genDau_pt                ,&b_genDau_pt                );
   fChain->SetBranchAddress("genDau_eta"               ,&genDau_eta               ,&b_genDau_eta               );
   fChain->SetBranchAddress("genDau_phi"               ,&genDau_phi               ,&b_genDau_phi               );

}


void SDfactory::GetMyEntry(int entry){
   if(!fChain) return ;
   fChain->GetEntry(entry);
};


void SDfactory::SDtreePlanter(){
   //waring
   int warnp1=0;
   int warnp2=0;


   //create jet algorithm
   double R=0.8;
   JetDefinition jet_def(antikt_algorithm,R);
   //create softdrop algorithm
   double z_cut=0.1;
   std::vector<contrib::SoftDrop> sd;
   for(int ibeta=0;ibeta<beta_bin;ibeta++){
      sd.push_back( contrib::SoftDrop(beta_SD[ibeta],z_cut,R) );
   }

   //print SoftDroper information 

   std::cout<<"There are "<<beta_bin<< " SoftDroper in total "<<std::endl;
   std::cout<<"*********************************************"<<std::endl;
   for(int ibeta=0;ibeta<beta_bin;ibeta++){
     std::cout<< sd[ibeta].description()<<std::endl;
   }
   std::cout<<"*********************************************"<<std::endl;

   long int n_jet_total=0;
   long int n_jet_nonsd_total[beta_bin]={0};
   long int n_jet_nonp_total[beta_bin]={0};
   long int n_jet_normal_total[beta_bin]={0};


   

  
   std::cout<<"Start plant SoftDrop tree"<<std::endl;
   std::cout<<"Total input file number is "<< fileList.size() << std::endl;
    //*******************************ENTER FILE LOOP**************************************
   for(int ifile=0;ifile<fileList.size();ifile++){

      long int n_jet_file=0;
      long int n_jet_nonsd_file[beta_bin]={0};
      long int n_jet_nonp_file[beta_bin]={0};
      long int n_jet_normal_file[beta_bin]={0};

      std::cout<<"Current file :"<< ifile+1 <<" / " << fileList.size() <<std::endl;
      fFile=TFile::Open( fileList[ifile].c_str() ,"READ" );
      TTree *tree=(TTree *)fFile->Get("trackTree");
      LoadData(tree);
      std::cout<<"Total entries in this file is "<< TotalEntry << std::endl;
      
      TFile *SDfile=new TFile(SDfileList[ifile].c_str(),"RECREATE");
      TTree *SoftDroptree=new TTree("SDtrackTree","SDtrackTree");
      SoftDroptree->Branch("dau_pt" ,               dau_pt );
      SoftDroptree->Branch("dau_chg",               dau_chg );
      SoftDroptree->Branch("dau_eta",               dau_eta);
      SoftDroptree->Branch("dau_phi",               dau_phi );
      SoftDroptree->Branch("dau_pid",               dau_pid );
      SoftDroptree->Branch("jetPt",                jetPt  );
      SoftDroptree->Branch("jetEta",               jetEta );
      SoftDroptree->Branch("jetPhi",               jetPhi   );
      SoftDroptree->Branch("jetChargedMultiplicity", jetChargedMultiplicity);
      SoftDroptree->Branch("has_parents",            has_parents);
      SoftDroptree->Branch("jetZg",                  jetZg);
      SoftDroptree->Branch("jetRg",                  jetRg);
      SoftDroptree->Branch("subjet1_pt",             subjet1_pt);
      SoftDroptree->Branch("subjet1_eta",            subjet1_eta);
      SoftDroptree->Branch("subjet1_phi",            subjet1_phi);
      SoftDroptree->Branch("subjet1_m",              subjet1_m);
      SoftDroptree->Branch("subjet2_pt",             subjet2_pt);
      SoftDroptree->Branch("subjet2_eta",            subjet2_eta);
      SoftDroptree->Branch("subjet2_phi",            subjet2_phi);
      SoftDroptree->Branch("subjet2_m",              subjet2_m);
      for(int ibeta=0;ibeta<beta_bin;ibeta++){
         SoftDroptree->Branch( Form("dau_flag_beta_%d",(int) beta_SD[ibeta] ) ,   dau_Flag[ibeta]  );
      }


      //**************************ENTER EVENT LOOP**********************************
      for(int ievent=0;ievent<TotalEntry;ievent++){
         if(ievent%10000==0) std::cout<<"ievent="<<ievent<<std::endl;
         GetMyEntry(ievent);
         std::vector<PseudoJet> particles;


         jetPt->clear();
         jetEta->clear();
         jetPhi->clear();
         jetChargedMultiplicity->clear();
         jetZg->clear();
         jetRg->clear();
         has_parents->clear();
         subjet1_pt->clear();
         subjet1_eta->clear();
         subjet1_phi->clear();
         subjet1_m->clear();
         subjet2_pt->clear();
         subjet2_eta->clear();
         subjet2_phi->clear();
         subjet2_m->clear();
         dau_pt->clear();
         dau_chg->clear();
         dau_eta->clear();
         dau_pid->clear();
         dau_phi->clear();
         for(int ibeta=0;ibeta<beta_bin;ibeta++){
            dau_Flag[ibeta]->clear();
         }

         int n_par=0;

         //fill particles will all the daughter particles , this should be modified for more general purposes
         for(int itmp=0;itmp<(*genJetPt).size();itmp++){
            for(int jtmp=0;jtmp<(*genDau_pid)[itmp].size();jtmp++){
               TLorentzVector v;
               int tmp_pid=(*genDau_pid)[itmp][jtmp];
               v.SetPtEtaPhiM((*genDau_pt)[itmp][jtmp],(*genDau_eta)[itmp][jtmp],(*genDau_phi)[itmp][jtmp],particleData.m0(tmp_pid) );
               PseudoJet particle( v.Px(),v.Py(),v.Pz(),v.E() );
               particle.set_user_info( new MyUserInfo(tmp_pid) );
               particle.set_user_index(n_par);
               n_par++;
               particles.emplace_back(particle);
            }
         }

         ClusterSequence cs(particles,jet_def);
         std::vector<PseudoJet> jets=cs.inclusive_jets();

         //************************ENTER JET LOOP**********************************
         for(int ijet=0;ijet<jets.size();ijet++){
            n_jet_file++;
            jetPt->emplace_back ( ( jets[ijet] ).pt()  );
            jetEta->emplace_back( ( jets[ijet] ).eta() );
            jetPhi->emplace_back( ( jets[ijet] ).phi() );
            std::vector<PseudoJet> constituents=jets[ijet].constituents();//All the daughter particles of the current jet 
            //dauguter particles information
            std::vector<float> pt;
            std::vector<float> eta;
            std::vector<float> phi;
            std::vector<int>   chg;
            std::vector<int>   pid;
            //subjet information 
            std::vector<float> z_g;
            std::vector<float> r_g;
            std::vector<float> subj1_pt;
            std::vector<float> subj1_eta;
            std::vector<float> subj1_phi;
            std::vector<float> subj1_m;
            std::vector<float> subj2_pt;
            std::vector<float> subj2_eta;
            std::vector<float> subj2_phi;
            std::vector<float> subj2_m;
            std::vector<int>   has_p; 
            std::vector<std::vector<int> > d_flag;
            int Nch=0;



            for(int iconstituent=0;iconstituent<constituents.size();iconstituent++){
               pt.emplace_back(constituents[iconstituent].pt()  );
               eta.emplace_back(constituents[iconstituent].eta());
               phi.emplace_back(constituents[iconstituent].phi() );
               int tmp_pid=constituents[iconstituent].user_info<MyUserInfo>().pid();
               if( particleData.charge(tmp_pid) ){Nch++;}
               chg.emplace_back( particleData.charge(tmp_pid) );
               pid.emplace_back(tmp_pid);
               
            }

            std::vector<int> flag(constituents.size());
            for(int ibeta=0;ibeta<beta_bin;ibeta++){
               d_flag.emplace_back(  flag  );
            }


            jetChargedMultiplicity->emplace_back(Nch);

            dau_pt->emplace_back(pt);
            dau_eta->emplace_back(eta);
            dau_phi->emplace_back(phi);
            dau_chg->emplace_back(chg);
            dau_pid->emplace_back(pid);


            //*************************ENTER SoftDrop LOOP**************************************
            
            for(int ibeta=0;ibeta<beta_bin;ibeta++){
               PseudoJet sd_jet=sd[ibeta]( jets[ijet] );//softdrop jet
               
               PseudoJet parent1;
               PseudoJet parent2;
               
               if(sd_jet.has_structure_of<contrib::SoftDrop>()){ 
                  if( sd_jet.has_parents(parent1,parent2) ){
                     // This sd_jet has sd stucture and parents
                     has_p.emplace_back(+1);
                     subj1_pt.emplace_back(parent1.pt());
                     subj1_eta.emplace_back(parent1.eta());
                     subj1_phi.emplace_back(parent1.phi());
                     subj1_m.emplace_back(parent1.m());
                     subj2_pt.emplace_back(parent2.pt());
                     subj2_eta.emplace_back(parent2.eta());
                     subj2_phi.emplace_back(parent2.phi());
                     subj2_m.emplace_back(parent2.m());



                     std::vector<PseudoJet> p1_con=parent1.constituents();
                     std::vector<PseudoJet> p2_con=parent2.constituents();
                     int isucc=0;
                     for(int ip1=0;ip1<p1_con.size();ip1++){
                       int index=p1_con[ip1].user_index();
                       for(int iconstituent=0;iconstituent<constituents.size();iconstituent++){
                           if(constituents[iconstituent].user_index()==index){
                              d_flag[ibeta][iconstituent]=1;
                              isucc++;
                              break;
                           } 
                       }
                       

                     }

                     if( isucc!=p1_con.size() ) {std::cout<<"waring p1"<<std::endl;warnp1++;}

                     isucc=0;
                     for(int ip2=0;ip2<p2_con.size();ip2++){
                        int index=p2_con[ip2].user_index();
                        for(int iconstituent=0;iconstituent<constituents.size();iconstituent++){
                           if(constituents[iconstituent].user_index()==index){
                              d_flag[ibeta][iconstituent]=2;
                              isucc++;
                              break;
                           }
                        } 
                        
                     }

                     if( isucc!=p2_con.size() ) {std::cout<<"waring p2"<<std::endl;warnp2++;}


                     
                     z_g.emplace_back(sd_jet.structure_of<contrib::SoftDrop>().symmetry());
                     r_g.emplace_back(sd_jet.structure_of<contrib::SoftDrop>().delta_R() );

                     if( (sd_jet.structure_of<contrib::SoftDrop>().symmetry() > 0.) && (sd_jet.structure_of<contrib::SoftDrop>().symmetry() < 0.5) ) n_jet_normal_file[ibeta]++;





                  }

                  else {// no parents
                     has_p.emplace_back(-1);
                     subj1_pt.emplace_back(0);
                     subj1_eta.emplace_back(0);
                     subj1_phi.emplace_back(0);
                     subj1_m.emplace_back(0);
                     subj2_pt.emplace_back(0);
                     subj2_eta.emplace_back(0);
                     subj2_phi.emplace_back(0);
                     subj2_m.emplace_back(0);
                     z_g.emplace_back(-1);
                     r_g.emplace_back( -1 );
                     n_jet_nonp_file[ibeta]++;
                  }
               }



               else{// no structure of soft drop 

                  has_p.emplace_back(-1);
                  subj1_pt.emplace_back(0);
                  subj1_eta.emplace_back(0);
                  subj1_phi.emplace_back(0);
                  subj1_m.emplace_back(0);
                  subj2_pt.emplace_back(0);
                  subj2_eta.emplace_back(0);
                  subj2_phi.emplace_back(0);
                  subj2_m.emplace_back(0);
                  z_g.emplace_back(-1);
                  r_g.emplace_back( -1 );
                  n_jet_nonsd_file[ibeta]++;

               }


            }//*********************************END SoftDrop LOOP****************************************

            
            jetZg->emplace_back(z_g);
            jetRg->emplace_back(r_g);
            has_parents->emplace_back(has_p);
            subjet1_pt->emplace_back(subj1_pt);
            subjet1_eta->emplace_back(subj1_eta);
            subjet1_phi->emplace_back(subj1_phi);
            subjet1_m->emplace_back(subj1_m);
            subjet2_pt->emplace_back(subj2_pt);
            subjet2_eta->emplace_back(subj2_eta);
            subjet2_phi->emplace_back(subj2_phi);
            subjet2_m->emplace_back(subj2_m);
            for(int ibeta=0;ibeta<beta_bin;ibeta++){
               dau_Flag[ibeta]->emplace_back(d_flag[ibeta]);
            }


         } //**********************END JET LOOP************************************************
         SoftDroptree->Fill();
      }//************************END EVENT LOOP***********************************************
      SDfile->Write();
      SDfile->Close();
      fFile->Close();
      std::cout<<"***************************************"<<std::endl;
      std::cout<<"number of jets in current file is " << n_jet_file <<std::endl;
      n_jet_total+=n_jet_file;
      for(int ibeta=0;ibeta<beta_bin;ibeta++){
         std::cout<<"beta="<<beta_SD[ibeta]<<std::endl;
         std::cout<<"number of jets without subjets is " <<n_jet_nonp_file[ibeta]<<std::endl;
         std::cout<<"number of jets without sd structure is "<<n_jet_nonsd_file[ibeta]<<std::endl;
         std::cout<<"number of jets with normal zg is " <<n_jet_normal_file[ibeta]<<std::endl;

         n_jet_nonp_total[ibeta]+=n_jet_nonp_file[ibeta];
         n_jet_nonsd_total[ibeta]+=n_jet_nonsd_file[ibeta];
         n_jet_normal_total[ibeta]+=n_jet_normal_file[ibeta];


      }
      std::cout<<"***************************************"<<std::endl;


   }//*******************************END FILE LOOP *******************************************

   std::cout<<"***************************************"<<std::endl;

   std::cout<<"number of jets in total  is " << n_jet_total <<std::endl;
   for(int ibeta=0;ibeta<beta_bin;ibeta++){
         std::cout<<"beta="<<beta_SD[ibeta]<<std::endl;
         std::cout<<"number of jets without subjets is " <<n_jet_nonp_total[ibeta]<<std::endl;
         std::cout<<"number of jets without sd structure is "<<n_jet_nonsd_total[ibeta]<<std::endl;
         std::cout<<"number of jets with normal zg is " <<n_jet_normal_total[ibeta]<<std::endl;
   }
   std::cout<<"***************************************"<<std::endl;
   if(!(warnp1+warnp2) ) std::cout<<"GOOD LUCK !"<<std::endl;

}


#endif










