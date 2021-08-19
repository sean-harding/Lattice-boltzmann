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
    /*Start by initializing the distribution to its equilibrium value.
    The method of Mei. et.al will be used to find a correct initialization
    following this*/
    this->collide(wx,wy,[](double arg1, double arg2){return arg2;});
    this->vDist_new = this->vDist;  //Check to make sure these are not pointing to the same thing
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
        this->vDist_new[i] = this->vDist[refl[i]];            //Reflect particles off of the wall
        }
    }
}
/*
This function relaxes the state of the fluid at a given lattice site depending on how close the state is
  to thermodynamic equilibrium
*/
template<class FUNC> void latt2D::collide(vec &wx, vec &wy,FUNC f){
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
    //this->update('d',wx,wy);
    //this->update('v',wx,wy);
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
            this->vDist[i] = f(this->vDist_new[i],(2*3*mag_u)*this->dens);}
            //this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*(2-3*mag_u)*this->dens;}
        else if(i<5){
            this->vDist[i] = f(this->vDist_new[i],expr1(args1[i-1],args2[i-1]));}
            //this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*expr1(args1[i-1],args2[i-1]);}
        else{
            this->vDist[i] = f(this->vDist_new[i],expr2(args1[i-5],args3[i-5]));}            
            //this->vDist[i] = (1-this->om)*this->vDist_new[i] + this->om*expr2(args1[i-5],args3[i-5]);}     //f2 through f8 inclusive 
    }
}
void latt2D::update(char which, vec &wx, vec &wy){
    if(which=='d'){
        this->dens = 0;
        for(int i=0;i<vDist_new.size();i++){
            this->dens+=this->vDist_new[i];}
            }    
    else if(which=='v'){
        this->fVel = vec {0.0,0.0};
        for(int i=0;i<vDist_new.size();i++){
            this->fVel[0] += wx[i]*this->vDist_new[i];
            this->fVel[1] += wy[i]*this->vDist_new[i];}}
    }
int main(int argc, char *argv[]){
    if(argc==1){
        std::cout<<"Program terminated: please specify #iterations, lx, ly"<<std::endl;
        return 0;
    }

    /*Simulate decaying Taylor-Green flow as a benchmark test. The following functions describe a particular flow
    profile in 2D which solves the Navier-Stokes equation. The equations specified give the solution at time t=0
    with the NS solution being an exponential decay of the velocity and density profiles. tg_x, tg_y and tg_d
    are the values of the x-component of velocity, y-component of velocity and density of the TG solution at t=0
    */
    double u0 = 0.1;
    double kx = M_PI/10;
    double ky = M_PI/10;
    double p0 = 1.0;        //This I suppose is somewhat arbitrary
    double d0 = 1.0; 
    auto tg_x = [&kx,&ky,&u0](double x,double y){return -u0*sqrt(kx/ky)*cos(kx*x)*sin(ky*y);};
    auto tg_y = [&kx,&ky,&u0](double x,double y){return u0*sqrt(kx/ky)*sin(kx*x)*cos(ky*y);};
    auto tg_d = [&](double x,double y){return p0 + 3*0.25*d0*pow(u0,2)*(ky*cos(2*kx*x)/kx+kx*cos(2*ky*y)/ky);};
    //vec iVel = {0.1,0.1};
    vec wx = {0.0,1,0.0,-1.0,0.0,1.0/sqrt(2),-1.0/sqrt(2),-1.0/sqrt(2),1.0/sqrt(2)};
    vec wy = {0.0,0.0,1.0,0.0,-1.0,1.0/sqrt(2),1.0/sqrt(2),-1.0/sqrt(2),-1.0/sqrt(2)};
    vec reflect = {0,3,4,1,2,7,8,5,6};
    double mix = 0.9*pow(10.0,-1.0);
    
    int nIter = std::atoi(argv[1]);
    int lx = std::atoi(argv[2]);
    int ly = std::atoi(argv[3]);
    
    int nSites = lx*ly;
    lattice l;
    
    std::ofstream outputFile;
    outputFile.open("velocities.txt");
    outputFile<<lx<<","<<ly<<"\n";
    l.reserve(nSites);
    /*Put walls at the extreme left and right boundaries of the system if we wish. For the TG vortex
    simulation I have turned walls off.
    */
    
    for(int i=0;i<nSites;i++){
        vec iVel = {tg_x(i%lx,i/lx),tg_y(i%lx,i/lx)};
        if(i%lx==0||i%lx==lx-1){
            l.push_back(latt2D(i,tg_d(i%lx,i/lx),iVel,mix,wx,wy,'f',lx,ly));}   //Allows us to set solid walls in the simulation
        else{                                                       //with no-slip boundary conditions implemented
            l.push_back(latt2D(i,tg_d(i%lx,i/lx),iVel,mix,wx,wy,'f',lx,ly));}
        }   
    
    double tol = pow(10.0,-10)/nSites;
    double delta = 10;
    int iter = 0;
    
    while(delta>tol&&iter<400){
        iter+=1;
        std::cout<<"Delta = "<<delta<<" Tolerance = "<<tol<<std::endl;
        delta=0;
        for(int i=0;i<nSites;i++){
            if(l[i].type!='w'){
                l[i].update('d',wx,wy);    
                l[i].collide(wx,wy,[&mix](double arg1,double arg2){return (1-mix)*arg1 + mix*arg2;});}
        }
        for(int i=0;i<nSites;i++){
            if(l[i].type!='w'){
                l[i].stream(l,reflect);}
                for(int k=0;k<l[i].vDist.size();k++){
                    delta+=pow(l[i].vDist_new[k]-l[i].vDist[k],2.0);}
        }
    }
    if(delta>tol){
        std::cout<<"Initialization failed in "<<iter<<" iterations. Convergence threshold reached: " << delta<<std::endl;
        return 0;
    }
    //Main body of program. Performs a number of streaming and collision operations if the site has fluid in it
    
    for(int n=0;n<nIter;n++){
        for(int i=0;i<nSites;i++){
            if(l[i].type!='w'){
            l[i].update('d',wx,wy);
            l[i].update('v',wx,wy);
            l[i].collide(wx,wy,[&mix](double arg1,double arg2){return (1-mix)*arg1 + mix*arg2;});}}
        for(int i=0;i<nSites;i++){
            if(l[i].type!='w'){
            l[i].stream(l,reflect);}}
        //Write particle velocities to a file
        for(int i=0;i<nSites;i++){
            outputFile<<l[i].fVel[0]<<","<<l[i].fVel[1]<<"\n";
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