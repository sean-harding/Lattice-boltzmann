#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "lBoltz.h"

typedef std::vector<double> vec;
typedef std::vector<int> ivec;
typedef std::vector<latt2D> lattice;

//Map a 2-parameter function across multiple argument pairs given in arg1 and arg2
template<class FUNC, class T>
auto map(FUNC func, const std::vector<T> &arg1, const std::vector<T> &arg2)
{
    int nvalues = std::min( arg1.size(), arg2.size());
    auto result = std::vector<T>(nvalues);
    for (int i=0; i<nvalues; ++i)
    {result[i] = func(arg1[i], arg2[i]);}
    return result;
}

//Constructor function for lattice site objects. See header file for member definitions
latt2D::latt2D(int index,double dens, vec iVel,double om, vec &wx, vec &wy,char type,int lx,int ly){
    //Initialize the current state of the system
    this->vDist = vec {0,0,0,0,0,0,0,0,0};
    this->vDist_new = vec {0,0,0,0,0,0,0,0,0};
    this->index = index;           
    this->type = type;
    this->om = om;
    //Initialize macroscopic variables
    this->dens = dens;
    this->fVel = iVel;

    //This works out nearest-neighbors and next-nearest neighbors. index = y*ly + x
    int x = index%lx;
    int y = index/lx;
        
    //ivec new_x = {x,(x+1)%lx,x,(x-1)%lx,x,(x+1)%lx,(x-1)%lx,(x-1)%lx,(x+1)%lx};
    //ivec new_y = {y,y,(y-1)%ly,y,(y+1)%ly,(y-1)%ly,(y-1)%ly,(y+1)%ly,(y+1)%ly};
    
    ivec new_x = {x,(x-1)%lx,x,(x+1)%lx,x,(x-1)%lx,(x+1)%lx,(x+1)%lx,(x-1)%lx};
    ivec new_y = {y,y,(y+1)%ly,y,(y-1)%ly,(y-1)%ly,(y-1)%ly,(y+1)%ly,(y+1)%ly};
    for(int i=0; i<new_x.size();i++){
        if(new_x[i]<0){new_x[i] = lx-1;}
        if(new_y[i]<0){new_y[i]=ly-1;}
        int neigh =new_y[i]*lx + new_x[i];
        this->streamChannels.push_back(neigh);
        } 
}   
/*  The velocity vector in position [i] in the array points from zero to adjacent sites. The order of these sites
    is:
        6 2 5
        3 0 1
        7 4 8
    for example, the magnitude of the velocity vector (-1,0) is stored in array index [3].
vec reflect = {0,3,4,1,2,7,8,5,6};
*/

//This function takes a lattice site l and streams particles from it to neighbors
void latt2D::stream(lattice &l, vec &refl){
    for(int i=0;i<(this->vDist).size();i++){
        int whichSite = this->streamChannels[i];  //Which site to stream from. Work this out when objects
        if(l[whichSite].type != 'w'){                         //Check if there is a solid wall at that pixel
        this->vDist_new[i] = l[whichSite].vDist[i];
        }
        else{                                                 //Reflect particles off of the wall
        this->vDist_new[i] = this->vDist[refl[i]];            
        }
    }
}
/*
This function relaxes the state of the fluid at a given lattice site depending on how close the state is
  to thermodynamic equilibrium
*/
void latt2D::collide(int iter, vec &wx, vec &wy){
    /*
    Relax particles at a given lattice site
    STEP 1.) Calculate the new macroscopic quanities from vDist_new
    STEP 2.) Calculate equilibrium distributions from macroscopic quantities, (F*_eq)
    STEP 3.) Update distribution as F_new = (1-t/T)*F_old +(t/T)*F*_eq
    
    If t/T is kept fixed, then we can pre-compute t/T to avoid having to perform a bunch of similar divisions
    at each iteration. This is done in the current implementation but should be changed if we want to be able to
    change the timestep 
    */

    //STEP 1
    if(iter!=0){
    this->fVel = vec {0.1,0};
    this->dens = 0;
    for(int i=0;i<vDist_new.size();i++){
        this->dens+=this->vDist_new[i];
        this->fVel[0] += wx[i]*this->vDist_new[i];
        this->fVel[1] += wy[i]*this->vDist_new[i];
    }}
    //STEP 2
    //mag_u is |velocity|^2
    double mag_u = 0.0;
    for(int i =0;i<(this->fVel).size();i++){
        mag_u += pow(this->fVel[i],2.0);        
    }
    //Define lambda expression for the f1->f4
    auto expr1 = [this, &mag_u](double a,double b){
        return this->dens*(2+6*a+9*b+3*mag_u)/18;};
    
    //Define lambda expression for the f5->f8
    auto expr2 = [this, &mag_u](double a,double b){
        return this->dens*(1+3*(a+b)+9*a*b+3*mag_u)/36;};
    
    //These are the function arguments
    std::vector<double> args1 = {this->fVel[0],this->fVel[1],-1.0*this->fVel[0],-1.0*this->fVel[1]};
    std::vector<double> args2 = {pow(this->fVel[0],2.0),pow(this->fVel[1],2.0),pow(this->fVel[0],2.0),pow(this->fVel[1],2.0)};
    std::vector<double> args3 = {this->fVel[1],-1.0*this->fVel[0],-1.0*this->fVel[1],this->fVel[0]};
    
    //STEP 3
    //The following could be written more cleanly as each line of the if/else has the same structure
    //Problem in the arguments??
    for(int i=0;i<this->vDist.size();i++){
        if(i==0){
            this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*(2-3*mag_u)*this->dens;}
        else if(i<5){
            this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*expr1(args1[i-1],args2[i-1]);}
        else{
            this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*expr2(args1[i-5],args3[i-5]);}     //f2 through f8 inclusive 
    }
}
int main(int argc, char *argv[]){
    if(argc==1){
        std::cout<<"Program terminated: please specify #iterations, lx, ly"<<std::endl;
        return 0;
    }
    vec wx = {0.0,1,0.0,-1.0,0.0,1.0/sqrt(2),-1.0/sqrt(2),-1.0/sqrt(2),1.0/sqrt(2)};
    vec wy = {0.0,0.0,1.0,0.0,-1.0,1.0/sqrt(2),1.0/sqrt(2),-1.0/sqrt(2),-1.0/sqrt(2)};
    vec iVel = {0.1,0};
    vec reflect = {0,3,4,1,2,7,8,5,6};
    int nIter = std::atoi(argv[1]);
    int lx = std::atoi(argv[2]);
    int ly = std::atoi(argv[3]);
    int nSites = lx*ly;
    lattice l;
    std::ofstream outputFile;
    outputFile.open("velocities.txt");
    outputFile<<lx<<","<<ly<<"\n";
    l.reserve(nSites);
    //Put walls at the extreme left and right boundaries of the system
    for(int i=0;i<nSites;i++){
        if(i%lx==0||i==lx){
            l.push_back(latt2D(i,1.0,iVel,pow(10,-2),wx,wy,'f',lx,ly));}   //Set wall
        else{
            l.push_back(latt2D(i,1.0,iVel,pow(10,-2),wx,wy,'f',lx,ly));}
        }   

    //Main body of program. Performs a number of streaming and collision operations if the site has fluid in it
    for(int n=0;n<nIter;n++){
        for(int i=0;i<nSites;i++){
            if(l[i].type!='w'){
            l[i].collide(n,wx,wy);}}
        for(int i=0;i<nSites;i++){
            if(l[i].type!='w'){
            l[i].stream(l,reflect);}}
        //Write particle velocities to a file
        for(int i=0;i<nSites;i++){
            outputFile<<l[i].fVel[0]<<","<<l[i].fVel[1]<<"\n";
            //outputFile<<l[i].dens<<",";
        };
        outputFile<<"\n";}
    outputFile.close();
    return 0;
}
    /*
    For parallelization using tbb, we wish to pass:
    [&](tbb::blocked_range<int> r){
        for(int i = r.start();i<r.end();i++){
            l[i].stream(l,reflect);
        }}
    Followed by:
    [&](tbb::blocked_range<int> r){
        for(int i = r.start();i<r.end();i++){
            l[i].collide(wx,wy);
        }}
    on each iteration to a tbb map function
    */