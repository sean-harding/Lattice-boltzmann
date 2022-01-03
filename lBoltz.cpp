#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "lBoltz.h"

typedef std::vector<double> vec;
typedef std::vector<int> ivec;
typedef std::vector<latt2D> lattice;

//Constructor function for lattice site objects. See header file for member definitions
latt2D::latt2D(int index,double dens, vec iVel,double om, vec &wx, vec &wy,char type,int lx,int ly){
    //Initialize the current state of the system
    this->vDist = vec {0,0,0,0,0,0,0,0,0};
    this->vDist_new = vec {0,0,0,0,0,0,0,0,0};
    this->index = index;           
    this->type = type;
    //Initialize macroscopic variables
    this->dens = dens;
    this->fVel = iVel;
    //This works out nearest-neighbors and next-nearest neighbors. index = y*ly + x
    int x = index%lx;
    int y = index/lx;
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

//This function takes a lattice site l and streams particles to it from neighbors
void latt2D::stream(lattice &l, vec &refl){
    for(int i=0;i<(this->vDist).size();i++){
        int whichSite = this->streamChannels[i];              //Which site to stream from
        if(l[whichSite].type != 'w'){                         //Check if there is a solid wall at that pixel
        this->vDist_new[i] = l[whichSite].vDist[i];
        }
        else{                                                 //Reflect particles off of any walls
        this->vDist_new[i] = this->vDist[refl[i]];            
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
    STEP 1.) Calculate the new macroscopic quanities from vDist_new (this is done externally to the function to reuse the same code for initialization)
    STEP 2.) Calculate equilibrium distributions from macroscopic quantities, (F*_eq)
    STEP 3.) Update distribution as F_new = (1-t/T)*F_old +(t/T)*F*_eq
    */

    //STEP 2
    //Calculate |velocity|^2
    double mag_u = 0.0;
    for(int i =0;i<(this->fVel).size();i++){
        mag_u += pow(this->fVel[i],2.0);        
    }
    //These are expressions for the equilibrium distribution functions for f_1->f_4
    auto expr1 = [this, &mag_u](double a,double b){
        return this->dens*(2+6*a+9*b+3*mag_u)/18;};
    
    //These are expressions for the equilibrium distribution functions for f_5->f_8
    auto expr2 = [this, &mag_u](double a,double b){
        return this->dens*(1+3*(a+b)+9*a*b+3*mag_u)/36;};
    
    //These are the function arguments to calculate each discretant of the distribution function
    std::vector<double> args1 = {this->fVel[0],this->fVel[1],-1.0*this->fVel[0],-1.0*this->fVel[1]};
    std::vector<double> args2 = {pow(this->fVel[0],2.0),pow(this->fVel[1],2.0),pow(this->fVel[0],2.0),pow(this->fVel[1],2.0)};
    std::vector<double> args3 = {this->fVel[1],-1.0*this->fVel[0],-1.0*this->fVel[1],this->fVel[0]};
    
    //STEP 3. The function f is used to perform the update
    for(int i=0;i<this->vDist.size();i++){
        if(i==0){
            this->vDist[i] = f(this->vDist_new[i],(2*3*mag_u)*this->dens);}
        else if(i<5){
            this->vDist[i] = f(this->vDist_new[i],expr1(args1[i-1],args2[i-1]));}
        else{
            this->vDist[i] = f(this->vDist_new[i],expr2(args1[i-5],args3[i-5]));}            
    }
}
//Is used to update the density (with flag 'd') or velocity (with flag 'v')
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
    //------------------|TEST CASE PARAMETERS|------------------
    double u0 = 0.3;
    double kx = 0.8*M_PI/std::atof(argv[2]);
    double ky = 0.8*M_PI/std::atof(argv[3];
    double p0 = 1.0;        //This choice I suppose is somewhat arbitrary, but works fine
    //double d0 = 1.0; 
    double x0 = std::atof(argv[2])*0.5;
    double y0 = std::atof(argv[2])*0.5;
    auto tg_x = [&kx,&ky,&u0,&x0,&y0](double x,double y){return -u0*sqrt(kx/ky)*cos(kx*(x-x0))*sin(ky*(y-y0));};
    auto tg_y = [&kx,&ky,&u0,&x0,&y0](double x,double y){return u0*sqrt(kx/ky)*sin(kx*x)*cos(ky*y);};
    auto tg_d = [&](double x,double y){return p0 + 0.25*pow(u0,2)*(ky*cos(2*kx*(x-x0))/kx+kx*cos(2*ky*(y-y0))/ky);};
    
    //------------------|SIMULATION VARIABLES|------------------
    /* wx and wy are weights used in the calculation of the equilibrium distibution function
       store once and then access by reference in the relevant functions
    */
    vec wx = {0.0,1,0.0,-1.0,0.0,1.0/sqrt(2),-1.0/sqrt(2),-1.0/sqrt(2),1.0/sqrt(2)};
    vec wy = {0.0,0.0,1.0,0.0,-1.0,1.0/sqrt(2),1.0/sqrt(2),-1.0/sqrt(2),-1.0/sqrt(2)};
    vec reflect = {0,3,4,1,2,7,8,5,6};  //When a particle of velocity v[k] is reflected off a wall, the outgoing particle's velocity v[reflec[k]] is relplaced by v[k].
    /*In simulation the dist. function is updated by mixing the old distribution function with the equilibrium dist. function
    calculated with the density and velocity obtained after a streaming step. This determines what proportion of the equilibrium
    function to mix with the old non equilibrium one. It is proportional to dT/T where dT is the simulation timestep and T is the
    characteristic relaxation time of the fluid (the fluid relaxes to local equilibrium on a timescale T*/
    double mix = 0.9*pow(10.0,-1.0);
    
    //Command line inputs for the number of iterations and the lattice dimensions
    int nIter = std::atoi(argv[1]);
    int lx = std::atoi(argv[2]);
    int ly = std::atoi(argv[3]);
    int nSites = lx*ly;
    lattice l;
    
    //Data is output to a txt file which is then read into python for animating the simulation
    std::ofstream outputFile;
    outputFile.open("velocities.txt");
    outputFile<<lx<<","<<ly<<"\n";
    
     //------------------|INITIALIZATION|------------------
    l.reserve(nSites);
    /*Put walls at the extreme left and right boundaries of the system if we wish. For the TG vortex
    simulation I have turned walls off.
    */
    
    for(int i=0;i<nSites;i++){
        vec iVel = {tg_x(i%lx,i/lx),tg_y(i%lx,i/lx)};
        if(i%lx==0||i%lx==lx-1){
            l.push_back(latt2D(i,tg_d(i%lx,i/lx),iVel,mix,wx,wy,'f',lx,ly));}   //Allows us to set solid walls in the simulation
        else{
            l.push_back(latt2D(i,tg_d(i%lx,i/lx),iVel,mix,wx,wy,'f',lx,ly));}
        }   
    
    /* It is not correct to use any old density, velocity pair to initialize the simulation. This is because the Poisson equation which
    relates the velocity to the density is not generally satisfied. The approach of Mei. et al. performs a series of collide/stream steps
    in which the particle density is allowed to relax, keeping the initial velocity fixed. This remedies the issue and removes instabilities
    */
    double tol = pow(10.0,-10)/nSites;  //Simulation is started when the change in the distribution function is small. This is the numerical tolerance
    double delta = 10;                  //The mod square difference in the dist. function between iterations when initializing
    int iter = 0;                       //To prevent the program from never finishing, we stop after a set number of iterations if the initialization fails
    
    while(delta>tol&&iter<2000){
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
            l[i].update('d',wx,wy); //Update the macroscopic variables each time
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
