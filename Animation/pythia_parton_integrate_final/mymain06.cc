#include "PartonShower.h"
#include "Free.h"
#include "frag.h"
int main()
{
    PartonShower ps=PartonShower();
    stringFragmentation sf=stringFragmentation(); 
    if(! ps.initPythia() ){
        return 1;
    }
    if(! sf.initPythia() ){
        return 1;
    }



    while(1){
        if(sf.Pass_flag==1) break;
        ps.getEvent();
        ps.FreeStreaming();
        sf.frag();
        sf.JetFinder();

    }
    sf.PrintResult();
    ps.JetPt=sf.JetPt;
    ps.JetEta=sf.JetEta;
    ps.JetPhi=sf.JetPhi;
    ps.JetRapidity=sf.JetRapidity;
    ps.initCanvas();
    ps.initAnimation();
    ps.make_animation();    
    return 0;
}


  
 

