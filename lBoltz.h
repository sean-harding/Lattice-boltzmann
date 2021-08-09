#ifndef _lboltz
#define _lboltz

#include <vector>

typedef std::vector<double> vec;
typedef std::vector<int> ivec;

/* NOTES:
1.) We do not need to distinguish vDist from vDist_inc if the code is paralellized. We would simply pass
    vDist to each process by value before reducing the results back to the main thread. The main thread then
    performs the state update once all new states are calcuated, meaning the data in vDist can be overwritten
2.)
*/

class latt2D{
    public:  
        int index;                       
        vec vDist;                  //Probability distribution for velocities: for reading old state from 
        vec vDist_new;              //Probability distribution for velocities: for writing new state to
        vec fVel;                 //Fluid velocity in 2D
        double dens;                //Fluid density
        latt2D(int index,vec initState,double dT, vec &wx, vec &wy,char type,int lx,int ly);    //Constructor
        void stream(std::vector<latt2D> &l,vec &refl);
        void collide(vec &wx,vec &wy);
    private:
        /*wx and wy are the weight factors used to calculate velocities. All 2D lattice sites have the same
        eqWeight is used in the collision method to evaluate the equilibrium distribution*/
        double om;                       //Relaxation time divided by timestep
        ivec streamChannels;
        char type;
};
#endif