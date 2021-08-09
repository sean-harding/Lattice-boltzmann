#include <vector>
#include <iostream>
#include "lBoltz.h"
#include <cmath>

typedef std::vector<double> vec;
typedef std::vector<int> ivec;
typedef std::vector<latt2D> lattice; //Use pointers to allow pass by reference

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
latt2D::latt2D(int index,vec initState,double om, vec &wx, vec &wy,char type,int lx,int ly){
    //Initialize the current state of the system
    this-> vDist = initState;
    this->vDist_new = initState;   
    this->index = index;           
    this->type = type;
    this->om = om;
    //Initialize macroscopic variables
    this->dens = 0;
    this->fVel = vec {0,0};
    for(int i=0;i<initState.size();i++){
        this->dens+=initState[i];
        this->fVel[0] += wx[i]*initState[i];
        this->fVel[1] += wy[i]*initState[i];}
    
    //This works out nearest-neighbors and next-nearest neighbors
    int y = index%lx;
    int x = (index - y)/lx;
    ivec args1 = {x,x-1,x,x+1,x,x-1,x+1,x+1,x-1};
    ivec args2 = {y,y,y-1,y,y+1,y-1,y-1,y+1,y+1};
    for(int i = 0;i<args1.size();i++){
        if(args1[i]<0){
            args1[i] = lx-1;}
        if(args2[i]<0){
            args2[i] = ly-1;}
        if(args1[i]>lx){
            args1[i] = 0;}
        if(args2[i]>ly){
            args2[i] = 0;}
        this->streamChannels.push_back(args2[i]*lx+args1[i]);
        if(index==4){
        std::cout<<args2[i]*lx + args1[i]<<std::endl;}
    }
}
/*
625
301
748
*/

//Takes a lattice site l and streams particles from it to neighbors
void latt2D::stream(lattice &l, vec &refl){
    for(int i=0;i<(this->vDist).size();i++){
        int whichSite = this->streamChannels[i];  //Which site to stream from. Work this out when objects
        if(l[whichSite].type != 'w'){                         //Check if there is a solid wall at that pixel
        this->vDist_new[i] = l[whichSite].vDist[i];}
        else{                                                 //Reflect particles off of the wall
        this->vDist_new[i] = this->vDist[refl[i]];            
        }
    }
}
void latt2D::collide(vec &wx, vec &wy){
    double o1;          //1- dt/T
    double o2;          //dt/T
    /*Relax particles at a given lattice site
    STEP 1.) Calculate the new macroscopic quanities from vDist_new
    STEP 2.) Calculate equilibrium distributions from macroscopic quantities, (F*_eq)
    STEP 3.) Update distribution as F_new = (1-t/T)*F_old +(t/T)*F*_eq
    
    If t/T is kept fixed, then we can pre-compute t/T to avoid having to perform a bunch of similar divisions
    at each iteration. This is done in the current implementation but should be changed if we want to be able to
    change the timestep 
    */

    //STEP 1
    for(int i=0;i<vDist_new.size();i++){
        this->dens+=this->vDist_new[i];
        this->fVel[0] += wx[i]*this->vDist_new[i];
        this->fVel[1] += wy[i]*this->vDist_new[i];
    }
    
    //STEP 2
    //mag_u is |velocity|^2
    double mag_u = 0;
    for(int i =0;i<this->fVel.size();i++){
        mag_u += pow(this->fVel[i],2);
    }
    
    //Define lambda expression for the f1->f4
    auto expr1 = [this, &mag_u](double a,double b){
        return this->dens*(2+6*a+9*b+3*mag_u)/18;};
    
    //Define lambda expression for the f5->f8
    auto expr2 = [this, &mag_u](double a,double b){
        return this->dens*(1+3*(a+b)+9*a*b+3*mag_u)/36;};
    
    //These are the function arguments
    std::vector<double> args1 = {this->fVel[0],this->fVel[1],-1*this->fVel[0],-1*this->fVel[1]};
    std::vector<double> args2 = {pow(this->fVel[0],2),pow(this->fVel[1],2),pow(this->fVel[0],2),pow(this->fVel[1],2)};
    std::vector<double> args3 = {this->fVel[1],-1*this->fVel[0],-1*this->fVel[1],this->fVel[0]};
    
    //STEP 3
    //The following could be written more cleanly as each line of the if/else has the same structure
    for(int i=0;i<this->vDist.size();i++){
        if(i==0){
            this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*(2-3*mag_u)*this->dens;}
        else if(i<5){
            this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*expr1(args1[i-1],args2[i-1]);}
        else{
            this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*expr2(args1[i-1],args3[i-1]);}     //f2 through f8 inclusive 
    }
}
int main(){
    vec initState = {1/sqrt(9),1/sqrt(9),1,1/sqrt(9),1/sqrt(9),1/sqrt(9),1/sqrt(9),1/sqrt(9),1/sqrt(9)};
    vec wx = {0,1,0,-1,0,1/sqrt(2),-1/sqrt(2),-1/sqrt(2),1/sqrt(2)};
    vec wy = {0,0,1,0,-1,1/sqrt(2),1/sqrt(2),-1/sqrt(2),-1/sqrt(2)};
    vec reflect = {0,3,4,1,2,7,8,5,6};
    int nIter = 10;
    lattice l;
    std::vector<vec> densities; //Store the particle densities at each timestep
    l.reserve(25);
    densities.reserve(nIter);
    for(int i=0;i<25;i++){
        l.push_back(latt2D(i,initState,0.1,wx,wy,'f',5,5));}
    
    //Main body of program. Performs a number of streaming and collision operations
     for(int i=0;i<25;i++){
        l[i].stream(l,reflect);
        std::cout<<i<<std::endl;}


    //l[0].stream(l,reflect);
    //l[1].stream(l,reflect);

    return 0;
    
    /*For parallelization, we wish to pass:
    [&](tbb::blocked_range<int> r){
        for(int i = r.start();i<r.end();i++){
            L[i].stream(L);
        }}
    Followed by:
    [&](tbb::blocked_range<int> r){
        for(int i = r.start();i<r.end();i++){
            L[i].collide(vx,vy);
        }}
int index,vec initState,double dT, vec &wx, vec &wy,char type,int lx,int ly
*/
}