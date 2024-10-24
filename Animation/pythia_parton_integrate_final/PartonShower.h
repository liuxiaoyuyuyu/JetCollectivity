#include <iostream>    //....C++ headers
#include <fstream>
#include <string>
#include <ctime>
#include <assert.h>
#include <unistd.h>   // access() funciton
#include <time.h>
#include <sys/stat.h>  // mkdir() function 
#include <sstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include "coordinateTools.h"
#include "Pythia8/Pythia.h"    //....PYTHIA8 headers
#include "fastjet/ClusterSequence.hh"


#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLorentzVector.h"





using namespace Pythia8;
using namespace std;
using namespace fastjet;


void rotate(double px,double py,double pz,double pr[4],int icc);



class PartonShower{


public:
	std::vector< int         > ID;
    std::vector< double      > Px;
    std::vector< double      > Py;
    std::vector< double      > Pz;
    std::vector< double      > E;
    std::vector< double      > X;
    std::vector< double      > Y;
    std::vector< double      > Z;
    std::vector< double      > T;
    std::vector< int         > IDD1;
    std::vector< int         > IDD2;
    std::vector< int    	 > IDM1;
    std::vector< int    	 > IDM2;
    std::vector< int    	 > STATUS;
    std::vector< std::string > NAME;
    std::vector< bool        > isFinal;
    std::vector< bool        > isParton;
    std::vector< double      > Dau1_appear_T;
    std::vector< double      > Dau2_appear_T;
    std::vector< double      > self_appear_T;
    std::vector< int         > cal_flag;

    double current_time;
    double atime_step;
    double JetPt;
    double JetEta;
    double JetPhi;
    double JetRapidity;

    Pythia pythia;
    Event *event;
    

    //
    //TCanvas *c;
    //TH2D *hxy;
    //TH2D *hxz;
    //TH2D *hyz;
    //TH3D *hxyz;


   
    
    PartonShower();
    virtual ~PartonShower();
    virtual int initPythia();
    virtual void initCanvas();
    virtual void initAnimation();
    virtual void getEvent();
    virtual void Replace1(int j, int iddd1);
    virtual void Replace2(int j, int idd1, int idd2);
    virtual void Erase1(int j);
    virtual void calculate_daughter_time(int j);
    virtual void make_animation();
    virtual void FreeStreaming();
};

PartonShower::PartonShower()
    : current_time(0),atime_step(0.01),event(0),JetPt(0),JetEta(0),JetPhi(0),JetRapidity(0)
{
    ID=std::vector<int > {};
    Px=std::vector<double>{};
    Py=std::vector<double>{};
    Pz=std::vector<double>{};
    E=std::vector<double>{};
    X=std::vector<double>{};
    Y=std::vector<double>{};
    Z=std::vector<double>{};
    T=std::vector<double>{};
    IDD1=std::vector<int>{};
    IDD2=std::vector<int>{};
    IDM1=std::vector<int>{};
    IDM2=std::vector<int>{};
    STATUS=std::vector<int>{};
    NAME=std::vector<std::string>{};
    isFinal=std::vector<bool>{};
    isParton=std::vector<bool>{};
    Dau1_appear_T=std::vector<double>{};
    Dau2_appear_T=std::vector<double>{};
    self_appear_T=std::vector<double>{};
    cal_flag=std::vector<int>{};

}


PartonShower::~PartonShower(){

}


int PartonShower::initPythia(){
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");  //92549 792 458 871
    //pythia.readString("Random:seed = " + random_str);
    pythia.readString("Beams:eCM = 13000.0");//5020,2760 GeV
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212"); //2212 is proton; 1000791970 for 197Au; 1000822080 for 208Pb;
    // Standard settings
    pythia.readString("HardQCD:all = on");
    //pythia.readString("Tune:pp=18");
    pythia.readString("PhaseSpace:pTHatMin = 500.0");
    pythia.readString("PhaseSpace:pTHatMax = -1.");

    pythia.readString("Tune:pp=14");
    //pythia.readString("PDF:pSet = 20");
    pythia.readString("Tune:ee=7");
    pythia.readString("MultipartonInteractions:ecmPow=0.03");
    pythia.readString("MultipartonInteractions:bProfile=2");
    pythia.readString("MultipartonInteractions:pT0Ref=1.41");
    pythia.readString("MultipartonInteractions:coreRadius=0.76");
    pythia.readString("MultipartonInteractions:coreFraction=0.63");
    pythia.readString("ColourReconnection:range=5.18");
    pythia.readString("SigmaTotal:zeroAXB=off");
    pythia.readString("SpaceShower:rapidityOrder=on");
    pythia.readString("SpaceShower:alphaSorder=2");
    pythia.readString("SpaceShower:alphaSvalue=0.118");
    pythia.readString("SigmaProcess:alphaSvalue=0.118");
    pythia.readString("SigmaProcess:alphaSorder=2");
    pythia.readString("MultipartonInteractions:alphaSvalue=0.118");
    pythia.readString("MultipartonInteractions:alphaSorder=2");
    pythia.readString("TimeShower:alphaSorder=2");
    pythia.readString("TimeShower:alphaSvalue=0.118");
    pythia.readString("SigmaTotal:mode = 0");
    pythia.readString("SigmaTotal:sigmaEl = 22.08");
    pythia.readString("SigmaTotal:sigmaTot = 101.037");
    pythia.readString("PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118");
    pythia.readString("HadronLevel:Hadronize = off");

    if (!pythia.init()) {
        std::cerr << "Pythia initialization failed!" << std::endl;
        return 0;
    }

    else{return 1;}


}

void PartonShower::initCanvas(){
    //c=new TCanvas("","",800,800);
    //c->SetFrameLineWidth(3);

    //hxy=new TH2D("T1=","T1=",100,-5,5,100,-5,5);
    //hxy->SetStats(0);
    //hxy->GetXaxis()->SetTitle("X [GeV^{-1}]");
    //hxy->GetXaxis()->CenterTitle();
    //hxy->GetYaxis()->SetTitle("Y [GeV^{-1}]");
    //hxy->GetYaxis()->CenterTitle();

    //hxz=new TH2D("T2=","T2=",100,-5,5,100,-5,5);
    //hxz->SetStats(0);
    //hxz->GetXaxis()->SetTitle("X [GeV^{-1}]");
    //hxz->GetXaxis()->CenterTitle();
    //hxz->GetYaxis()->SetTitle("Z [GeV^{-1}]");
    //hxz->GetYaxis()->CenterTitle();

    //hyz=new TH2D("T3=","T3=",100,-5,5,100,-5,5);
    //hyz->SetStats(0);
    //hyz->GetXaxis()->SetTitle("Y [GeV^{-1}]");
    //hyz->GetXaxis()->CenterTitle();
    //hyz->GetYaxis()->SetTitle("Z [GeV^{-1}]");
    //hyz->GetYaxis()->CenterTitle();

    //hxyz=new TH3D("T3=","T3=",100,-5,5,100,-5,5,100,-5,5);
    //hxyz->SetStats(0);
    //hxyz->GetXaxis()->SetTitle("Y [ GeV^{-1} ]");
    //hxyz->GetXaxis()->SetTitleOffset(1.95);
    //hxyz->GetXaxis()->CenterTitle();
    //hxyz->GetYaxis()->SetTitle("Z [ GeV^{-1} ]");
    //hxyz->GetYaxis()->SetTitleOffset(1.5);
    //hxyz->GetYaxis()->CenterTitle();
    //hxyz->GetZaxis()->SetTitle("X [ GeV^{-1} ]");
    //hxyz->GetZaxis()->CenterTitle();




}

void PartonShower::initAnimation(){
    //c->Clear();
    //hxy->Reset();
    //hxz->Reset();
    //hyz->Reset();
    //hxyz->Reset();
    current_time=0;
    ID.clear();
    Px.clear();
    Py.clear();
    Pz.clear();
    E.clear();
    X.clear();
    Y.clear();
    Z.clear();
    T.clear();
    IDD1.clear();
    IDD2.clear();
    IDM1.clear();
    IDM2.clear();
    STATUS.clear();
    NAME.clear();
    isFinal.clear();
    isParton.clear();
    Dau1_appear_T.clear();
    Dau2_appear_T.clear();
    self_appear_T.clear();
    cal_flag.clear();

    // the production of  hardest  interaction 
    for( int init=1;init<(*event).size();init++){
        if( (*event)[init].mother1()==1 || (*event)[init].mother1()==2   ){
            ID.push_back( init                        );
            Px.push_back( (*event)[init].px()            );
            Py.push_back( (*event)[init].py()            );
            Pz.push_back( (*event)[init].pz()            );
            E.push_back(  (*event)[init].e()             );
            X.push_back(  0.                          );
            Y.push_back(  0.                          );
            Z.push_back(  0.                          );
            T.push_back(  0.                          );
            IDD1.push_back(  (*event)[init].daughter1()  );
            IDD2.push_back(  (*event)[init].daughter2()  );
            IDM1.push_back(  (*event)[init].mother1()    );
            IDM2.push_back(  (*event)[init].mother2()    );
            STATUS.push_back((*event)[init].status()     );
            NAME.push_back(  (*event)[init].name()       ); 
            isFinal.push_back( (*event)[init].isFinal()  );
            isParton.push_back( (*event)[init].isParton());
            Dau1_appear_T.push_back(0.);
            Dau2_appear_T.push_back(0.);
            self_appear_T.push_back(0.);
            cal_flag.push_back(0.     );
        }
    
    }

}





void PartonShower::getEvent(){
    int Ncounter=0;
    while(1){
        if(!pythia.next()){continue;}
        else {
            
            Ncounter=0;
            for( int j=0;j<pythia.event.size();j++){
                if( pythia.event[j].isFinal() && pythia.event[j].isParton() ) Ncounter++;
            }
            
            //if( Ncounter < 300 ) continue;
            ///else break;
            break;
        }
    }
    std::cout<<"Ncounter="<<Ncounter<<std::endl;
    event=&(pythia.event);
    std::cout<<"getEvent succ"<<std::endl;
}




void PartonShower::Replace1(int j, int iddd1){
    ID[j]=iddd1;
    Px[j]=(*event)[iddd1].px();
    Py[j]=(*event)[iddd1].py();
    Pz[j]=(*event)[iddd1].pz();
    E[j]=(*event)[iddd1].e();
    IDD1[j]=(*event)[iddd1].daughter1();
    IDD2[j]=(*event)[iddd1].daughter2();
    IDM1[j]=(*event)[iddd1].mother1();
    IDM2[j]=(*event)[iddd1].mother2();
    STATUS[j]=(*event)[iddd1].status();
    NAME[j]=(*event)[iddd1].name();
    isFinal[j]=(*event)[iddd1].isFinal();
    isParton[j]=(*event)[iddd1].isParton();
    self_appear_T[j]=Dau1_appear_T[j];
    cal_flag[j]=1;
}


void PartonShower::Replace2(int j,int idd1,int idd2){
    double temp_x=X[j];
    double temp_y=Y[j];
    double temp_z=Z[j];
    double temp_t=T[j];
    double temp_dau1_appear_T=Dau1_appear_T[j];
    double temp_dau2_appear_T=Dau2_appear_T[j];
    ID.erase(ID.begin()+j);
    Px.erase(Px.begin()+j);
    Py.erase(Py.begin()+j);
    Pz.erase(Pz.begin()+j);
    E.erase( E.begin()+j );
    X.erase(X.begin()+j);
    Y.erase(Y.begin()+j);
    Z.erase(Z.begin()+j);
    T.erase(T.begin()+j);
    IDD1.erase(IDD1.begin()+j);
    IDD2.erase(IDD2.begin()+j);
    IDM1.erase(IDM1.begin()+j);
    IDM2.erase(IDM2.begin()+j);
    STATUS.erase(STATUS.begin()+j);
    NAME.erase(NAME.begin()+j);
    isFinal.erase(isFinal.begin()+j);
    isParton.erase(isParton.begin()+j);


    Dau1_appear_T.erase( Dau1_appear_T.begin()+j );
    Dau2_appear_T.erase( Dau2_appear_T.begin()+j );
    self_appear_T.erase( self_appear_T.begin()+j );
    cal_flag.erase(cal_flag.begin()+j);

    //do replecing
    ID.insert( ID.begin()+j, idd1  );
    Px.insert( Px.begin()+j, (*event)[idd1].px()  );
    Py.insert( Py.begin()+j, (*event)[idd1].px()  );
    Pz.insert( Pz.begin()+j, (*event)[idd1].pz()  );
    E.insert(  E.begin()+j,  (*event)[idd1].e()   );
    X.insert(  X.begin()+j,  temp_x            );
    Y.insert(  Y.begin()+j,  temp_y            );
    Z.insert(  Z.begin()+j,  temp_z            );
    T.insert(  T.begin()+j,  temp_t            );
    IDD1.insert( IDD1.begin()+j, (*event)[idd1].daughter1() );
    IDD2.insert( IDD2.begin()+j, (*event)[idd1].daughter2() );
    IDM1.insert( IDM1.begin()+j, (*event)[idd1].mother1()   );
    IDM2.insert( IDM2.begin()+j, (*event)[idd1].mother2()   );
    STATUS.insert(STATUS.begin()+j,(*event)[idd1].status()  );
    NAME.insert(NAME.begin()+j,  (*event)[idd1].name()      );
    isFinal.insert(isFinal.begin()+j,(*event)[idd1].isFinal());
    isParton.insert(isParton.begin()+j,(*event)[idd1].isParton());
    Dau1_appear_T.insert(Dau1_appear_T.begin()+j,temp_dau1_appear_T);
    Dau2_appear_T.insert(Dau2_appear_T.begin()+j,temp_dau1_appear_T);
    self_appear_T.insert(self_appear_T.begin()+j,temp_dau1_appear_T);
    cal_flag.insert(cal_flag.begin()+j,1);


    ID.insert( ID.begin()+j, idd2  );
    Px.insert( Px.begin()+j, (*event)[idd2].px()  );
    Py.insert( Py.begin()+j, (*event)[idd2].px()  );
    Pz.insert( Pz.begin()+j, (*event)[idd2].pz()  );
    E.insert(  E.begin()+j,  (*event)[idd2].e()   );
    X.insert(  X.begin()+j,  temp_x            );
    Y.insert(  Y.begin()+j,  temp_y            );
    Z.insert(  Z.begin()+j,  temp_z            );
    T.insert(  T.begin()+j,  temp_t            );
    IDD1.insert( IDD1.begin()+j, (*event)[idd2].daughter1() );
    IDD2.insert( IDD2.begin()+j, (*event)[idd2].daughter2() );
    IDM1.insert( IDM1.begin()+j, (*event)[idd2].mother1()   );
    IDM2.insert( IDM2.begin()+j, (*event)[idd2].mother2()   );
    STATUS.insert(STATUS.begin()+j,(*event)[idd2].status()  );
    NAME.insert(NAME.begin()+j,  (*event)[idd2].name()      );
    isFinal.insert(isFinal.begin()+j,(*event)[idd2].isFinal());
    isParton.insert(isParton.begin()+j,(*event)[idd2].isParton());
    Dau1_appear_T.insert(Dau1_appear_T.begin()+j,temp_dau2_appear_T);
    Dau2_appear_T.insert(Dau2_appear_T.begin()+j,temp_dau2_appear_T);
    self_appear_T.insert(self_appear_T.begin()+j,temp_dau2_appear_T);
    cal_flag.insert(cal_flag.begin()+j,1);
}


void PartonShower::Erase1(int j ){
    ID.erase(ID.begin()+j);
    Px.erase(Px.begin()+j);
    Py.erase(Py.begin()+j);
    Pz.erase(Pz.begin()+j);
    E.erase( E.begin()+j );
    X.erase(X.begin()+j);
    Y.erase(Y.begin()+j);
    Z.erase(Z.begin()+j);
    T.erase(T.begin()+j);
    IDD1.erase(IDD1.begin()+j);
    IDD2.erase(IDD2.begin()+j);
    IDM1.erase(IDM1.begin()+j);
    IDM2.erase(IDM2.begin()+j);
    STATUS.erase(STATUS.begin()+j);
    NAME.erase(NAME.begin()+j);
    isFinal.erase(isFinal.begin()+j);
    isParton.erase(isParton.begin()+j);
    Dau1_appear_T.erase( Dau1_appear_T.begin()+j );
    Dau2_appear_T.erase( Dau2_appear_T.begin()+j );
    self_appear_T.erase( self_appear_T.begin()+j );
    cal_flag.erase(cal_flag.begin()+j);
}


void PartonShower::calculate_daughter_time(int j){
    int idd1=IDD1[j];
    int idd2=IDD2[j];
    double p0[4] = {0.0}; // particle's four momentum
    double p4[4] = {0.0}; // mother's four momentum
    double qt,time_step;
    //start calculating daughter appearing time
    p4[0] = (*event)[ idd1 ].e() +(*event)[ idd2 ].e();
    p4[1] = (*event)[ idd1 ].px()+(*event)[ idd2 ].px();
    p4[2] = (*event)[ idd1 ].py()+(*event)[ idd2 ].py();
    p4[3] = (*event)[ idd1 ].pz()+(*event)[ idd2 ].pz();
    double x_split1 = (*event)[ idd1 ].e()/p4[0];
    double x_split2 = (*event)[ idd2 ].e()/p4[0];
    if (x_split1>1) x_split1=1.0/x_split1; // revise the mother and daughter
    if (x_split2>1) x_split2=1.0/x_split2; // revise the mother and daughter
    p0[0] = (*event)[ idd1 ].e();
    p0[1] = (*event)[ idd1 ].px();
    p0[2] = (*event)[ idd1 ].py();
    p0[3] = (*event)[ idd1 ].pz();     
    rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
    qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
    rotate(p4[1],p4[2],p4[3],p0,-1);
    double kt_daughter = qt;
    if(x_split1<0.5){
            if (kt_daughter > 0.0001) {                     
                time_step = 2.0*p4[0]*x_split1*(1-x_split1)/pow(kt_daughter,2)*0.197;
                Dau1_appear_T[j]+=time_step;
                                                                 
            }
    }  else {
            if(kt_daughter > 0.0001) {
                time_step = 1/p4[0];
                Dau1_appear_T[j]+=time_step;
                                         
                        
            }   
    }

    p0[0] = (*event)[ idd2 ].e();
    p0[1] = (*event)[ idd2 ].px();
    p0[2] = (*event)[ idd2 ].py();
    p0[3] = (*event)[ idd2 ].pz();     
    rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
    qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
    rotate(p4[1],p4[2],p4[3],p0,-1);
    kt_daughter = qt;
    if(x_split2<0.5){
        if (kt_daughter > 0.0001) {                     
            time_step = 2.0*p4[0]*x_split2*(1-x_split2)/pow(kt_daughter,2)*0.197;
                                        
            Dau2_appear_T[j]+=time_step;                         
        }
    }  else {
        if(kt_daughter > 0.0001) {
            time_step = 1/p4[0];
                                        
            Dau2_appear_T[j]+=time_step;   
                        
        }   
    }


    //end calculating daughter appearing time
    cal_flag[j]=0; 

}


void PartonShower::make_animation(){
    initAnimation();

    int breakflag=0;
    int count=0;
    TLatex *lx=new TLatex();

    std::string outfile_name="animation_xyz.dat";
    std::string outfile_name2="animation_Jetxy.dat";
   
    ofstream outfile(outfile_name.c_str());
    ofstream outfile2(outfile_name2.c_str());
    if(!outfile.is_open() ){
        std::cout<<" can not open outfile:"<<outfile_name<<std::endl;
        return ;

    }
    if(!outfile2.is_open() ){
        std::cout<<" can not open outfile:"<<outfile_name2<<std::endl;
        return ;

    }
    

    
    while(1){
        breakflag=1;     
        if(count>500) break;
        std::cout<<"current_time:"<<current_time<<std::endl;
        std::cout<<"Px.size()="<<Px.size()<<std::endl;
        int appear_partons=0;

        for( int j=0; j<Px.size(); j++ ){                
            std::cout<<"Px="<<Px[j]<<" "<<"Py="<<Py[j]<<" "<<"Pz="<<Pz[j]<<" "<<"E="<<E[j]<<" "<<"IDD1="<<IDD1[j]<<" "<<"IDD2="<<IDD2[j]<<" "<<"IDM1="<<IDM1[j]<<" "<<"IDM2="<<IDM2[j]<<" "<<"STATUS="<<STATUS[j]<<" "<<"NAME:"<<NAME[j]<<" "<<"isFinal:"<<isFinal[j]<<" "<<"isParton:"<<isParton[j]<<std::endl; 
            if( (isFinal[j]==0) || (isParton[j]==0) ) breakflag=0;
            if(self_appear_T[j]>current_time ) continue;
            else appear_partons++; 
        }



        if(breakflag==1) break;
        outfile<<appear_partons<<std::endl;
        outfile2<<appear_partons<<std::endl;
        for(int j=0;j<Px.size();j++){
            if(self_appear_T[j]>current_time ) continue;
           
            TLorentzVector v;
            v.SetPxPyPzE(Px[j],Py[j],Pz[j],E[j]);

            double current_dau_y=v.Rapidity(); // momentum rapidity
            double current_dau_phi=v.Phi();    // momentum phi

            double delta_phi=TMath::ACos( TMath::Cos(current_dau_phi - JetPhi) );      // momentum phi
            double delta_y=current_dau_y- JetRapidity;                                 //momentum rapidity         
            double delta_eta=v.Eta()-JetEta;                                           // momentum eta
            double delte_R1=TMath::Sqrt(     delta_phi*delta_phi+delta_y*delta_y  );   // the distance in y-phi plane 
            double delta_R2=TMath::Sqrt(     delta_phi*delta_phi+delta_eta*delta_eta  );// the distance in eta_phi plane 


            double lab_dau_eta_s=TMath::ATanH(  Z[j]/current_time  );                 // spacetime rapidity

            outfile<<X[j]<<" "<<Y[j]<<" "<<Z[j]<<" "<<Px[j]<<" "<<Py[j]<<" "<<Pz[j]<<" "<<E[j]<<" "<<lab_dau_eta_s<<" "<<delte_R1<<" "<<NAME[j]<<std::endl;
            //spatial information in jet frame 
            fastjet::PseudoJet p = PseudoJet(X[j],Y[j], Z[j], 0);
            double jet_dau_Sxy    =  ptWRTJet(JetPt  ,JetEta, JetPhi     , p.pt()  , p.eta(), p.phi() );// sqrt(x**2,y**2) in jet frame
            double jet_dau_Seta   = etaWRTJet(JetPt  ,JetEta, JetPhi     , p.pt()  , p.eta(), p.phi() );// spatial eta in jet frame 
            double jet_dau_Sphi   = phiWRTJet(JetPt  ,JetEta, JetPhi     , p.pt()  , p.eta(), p.phi() );// phi in jet frame 
            double jet_dau_Stheta = thetaWRTJet(JetPt,JetEta, JetPhi     , p.pt()  , p.eta(), p.phi() );// theta in jet frame 
            
          
            //momentum information in jet frame 
            double jet_dau_Pt    =  ptWRTJet(JetPt, JetEta     , JetPhi     ,  v.Pt()  , v.Eta(), v.Phi() );
            double jet_dau_Peta   = etaWRTJet(JetPt, JetEta    , JetPhi     ,  v.Pt()  , v.Eta(), v.Phi() );
            double jet_dau_Pphi   = phiWRTJet(JetPt, JetEta    , JetPhi     ,  v.Pt()  , v.Eta(), v.Phi() );


            double x_star=jet_dau_Sxy*TMath::Cos(jet_dau_Sphi); // x* in jet frame 
            double y_star=jet_dau_Sxy*TMath::Sin(jet_dau_Sphi); // y* in jet frame 
            double z_star=jet_dau_Sxy/TMath::Tan(jet_dau_Stheta);//z* in jet frame 
            double jet_dau_eta_s=TMath::ATanH( z_star/current_time ); //eta*_s in jet frame 

            outfile2<<x_star<<" "<<y_star<<" "<<jet_dau_Seta<<" "<<jet_dau_Pt<<" "<<jet_dau_Peta<<" "<<jet_dau_Pphi<<" "<<E[j]<<" "<<delte_R1<<" "<<jet_dau_eta_s<<" "<<NAME[j]<<std::endl;


        }
        

        std::cout<<std::endl;
        std::cout<<std::endl;
        std::cout<<std::endl;

        //!!!!!BE CAREFUL WITH THE CHANGING  INDEX!!!!!!

        current_time=current_time+atime_step;
        count++;
        //****************BEGIN J LOOP****************
        for( int j=0;j<Px.size(); j++){
            int idd1=IDD1[j];
            int idd2=IDD2[j];
            

            //****************SATRT case 1*****************
            if(idd1==0 && idd2==0 ){
                X[j]=X[j]+Px[j]/E[j]*atime_step;
                Y[j]=Y[j]+Py[j]/E[j]*atime_step;
                Z[j]=Z[j]+Pz[j]/E[j]*atime_step;
                T[j]=T[j]+atime_step;
                continue;
            }
            //****************END   case 1*****************

            //****************SATRT case 2***************** (idd1==idd2 && idd2>0 ) || (idd1>0 && idd2)  
            if( (idd1==idd2 && idd2>0) || (idd1>0 && idd2==0)   ){
                int idm1=(*event)[idd1].mother1();
                int idm2=(*event)[idd1].mother2();
                //if( idm1!=idm2 ) {warning1++;std::cout<<std::endl;std::cout<<"warning1_1 occurs at "<<j<<std::endl;std::cout<<std::endl;return 1; }
                //if( idm1!=ID[j] || idm2!=ID[j] ) {warning1++;std::cout<<std::endl;std::cout<<"warning1_2 occurs at "<<j<<std::endl;std::cout<<std::endl;return 1;}
                if( idm1!=idm2 && idm1!=0 && idm2!=0 ) { 
                    //std::cout<<"two mothers, exit"<<std::endl;return ;
                    if( idm2==ID[j]){
                        Erase1(j);
                        j=j-1;
                        continue;
                    }

                }

                       
                while(1){

                    int iddd1=IDD1[j];
                    int iddd2=IDD2[j];

                    if(iddd1==0 && iddd2==0 ){break;}
                    if( iddd1!=0 && iddd2 !=0 && iddd1!=iddd2  ){ idd1=iddd1;idd2=iddd2;break; }

                    if(iddd1==iddd2 && iddd2>0 || (iddd1>0 &&iddd2==0 ) ){
                        Replace1(j,iddd1);

                    }   

                    else{
                        std::cout<<"while warning: iddd1="<<iddd1<<" iddd2="<<iddd2<<std::endl;
                        return ;
                    }

                }

            }
            //****************END   case 2*****************


            //****************SATRT case 4*****************
            if( (idd1<idd2 && idd1>0) || (idd2<idd1 && idd2>0) ){
                int idm1=(*event)[idd1].mother1();
                int idm2=(*event)[idd1].mother2();


                if(idd1<idd2 && (idd2-idd1)>1 ) {std::cout<<"more than two daughters,exit"<<std::endl; return ; }
                if( idm1!=idm2 && idm1!=0 && idm2!=0 ){
                    
                    //std::cout<<"two mothers,exit"<<std::endl;return ;
                    if( ID[j]==idm2 ){
                        Erase1(j);
                        j=j-1;
                        continue;
                    }
                    if( ID[j]==idm1){
                        Replace2(j,idd1,idd2);
                        j=j-1;
                        continue;
                    }

                }
                //if( idm1!=ID[j] ) { 
                //    std::cout<<"waring2_1 occurs at "<<j<<std::endl; return ;
                //}
                 


                if(cal_flag[j]>0){
                   
                    calculate_daughter_time(j);
                }

                if(Dau1_appear_T[j]>current_time && Dau2_appear_T[j]>current_time ) {
                    X[j]=X[j]+Px[j]/E[j]*atime_step;
                    Y[j]=Y[j]+Py[j]/E[j]*atime_step;
                    Z[j]=Z[j]+Pz[j]/E[j]*atime_step;
                    T[j]=T[j]+atime_step;
                }


                else{
                    Replace2(j,idd1,idd2);
                    j=j+1;
                }

            }
            //****************END   case 4*****************
        }//******************END J LOOP**************** 

        
    }//end while(1)
    outfile.close();
    outfile2.close();

}

void rotate(double px,double py,double pz,double pr[4], int icc) {
    //     input:  (px,py,pz), (wx,wy,wz), argument (i)
    //     output: new (wx,wy,wz)
    //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
    //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)

       
    double wx,wy,wz,E,pt, cosa,sina,cosb,sinb;
    double wx1,wy1,wz1;        

    wx=pr[1];
    wy=pr[2];
    wz=pr[3];

    E=sqrt(px*px+py*py+pz*pz);
    pt=sqrt(px*px+py*py);

    //w = sqrt(wx*wx+wy*wy+wz*wz);

    //  if(pt==0)
    if (pt<1e-6) {
        cosa=1;
        sina=0;
    } 
    else {
        cosa=px/pt;
        sina=py/pt;
    }

    if (E>1e-6) {
        cosb=pz/E;
        sinb=pt/E;
        if (icc==1) {
            wx1=wx*cosb*cosa+wy*cosb*sina-wz*sinb;
              wy1=-wx*sina+wy*cosa;
              wz1=wx*sinb*cosa+wy*sinb*sina+wz*cosb;
          } else {
              wx1=wx*cosa*cosb-wy*sina+wz*cosa*sinb;
              wy1=wx*sina*cosb+wy*cosa+wz*sina*sinb;
              wz1=-wx*sinb+wz*cosb;
          }
        wx=wx1;
        wy=wy1;
        wz=wz1;
    } else {
        cout << "warning: small E in rotation" << endl;
    }

    pr[1]=wx;
    pr[2]=wy;
    pr[3]=wz;      
}
 









