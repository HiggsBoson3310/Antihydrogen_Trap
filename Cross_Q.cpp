#include <iostream> 
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <complex>
#include <random> //Random number library
using namespace std;
const double PI=3.14159265359, h = 6.626e-34, mue = 928.5e-26, z=2*PI*(mue/h);
//Matrix-Vector Multiplication
vector<complex<double>> MM(vector<complex<double>> M, vector<complex<double>> a){
    vector<complex<double>> r{complex<double> (0,0),complex<double> (0,0),complex<double> (0,0)};
    r.at(0)=M.at(0)*a.at(0)+M.at(1)*a.at(1)+M.at(2)*a.at(2);
    r.at(1)=M.at(3)*a.at(0)+M.at(4)*a.at(1)+M.at(5)*a.at(2);
    r.at(2)=M.at(6)*a.at(0)+M.at(7)*a.at(1)+M.at(8)*a.at(2);
    return r;
}
//Inner product
complex<double> inner(vector<complex<double>> b, vector<complex<double>> a){
    complex<double> r=0;
    for(int j=0; j<b.size(); j++){
        r+=conj(b.at(j))*a.at(j);
    }
    return r;
}
double inner(vector<double> b, vector<double> a){
    double r=0;
    for(int j=0; j<a.size(); j++){
        r+=b.at(j)*a.at(j);
    }
    return r;
}
//Matrix-Matrix multiplication
vector<complex<double>> Msq(vector<complex<double>> M, vector<complex<double>> A){
    vector<complex<double>> r(A.size(), 0);
    r.at(0)=M.at(0)*A.at(0)+M.at(1)*A.at(3)+M.at(2)*A.at(6);
    r.at(1)=M.at(0)*A.at(1)+M.at(1)*A.at(4)+M.at(2)*A.at(7);
    r.at(2)=M.at(0)*A.at(2)+M.at(1)*A.at(5)+M.at(2)*A.at(8);
    r.at(3)=M.at(3)*A.at(0)+M.at(4)*A.at(3)+M.at(5)*A.at(6);
    r.at(4)=M.at(3)*A.at(1)+M.at(4)*A.at(4)+M.at(5)*A.at(7);
    r.at(5)=M.at(3)*A.at(2)+M.at(4)*A.at(5)+M.at(5)*A.at(8);
    r.at(6)=M.at(6)*A.at(0)+M.at(7)*A.at(3)+M.at(8)*A.at(6);
    r.at(7)=M.at(6)*A.at(1)+M.at(7)*A.at(4)+M.at(8)*A.at(7);
    r.at(8)=M.at(6)*A.at(2)+M.at(7)*A.at(5)+M.at(8)*A.at(8);
    return r;
}
//Vector addition
vector<complex<double>> add_M(vector<complex<double>> b, vector<complex<double>> a){
    vector<complex<double>> r(a.size(), 0);
    for(int i=0; i<a.size();i++){
        r.at(i)=b.at(i)+a.at(i);
    }
    return r;
}
vector<double> add_M(vector<double> b, vector<double> a){
    vector<double> r(a.size(), 0);
    for(int i=0; i<a.size();i++){
        r.at(i)=b.at(i)+a.at(i);
    }
    return r;
}
//Scalar multiplication
vector<complex<double>> MS(complex<double> b, vector<complex<double>> a){
    vector<complex<double>> r(a.size(), 0);
    for(int i=0; i<a.size();i++){
        r.at(i)=b*a.at(i);
    }
    return r;
}
vector<double> MS(double b, vector<double> a){
    vector<double> r(a.size(), 0);
    for(int i=0; i<a.size();i++){
        r.at(i)=b*a.at(i);
    }
    return r;
}
//Cross product
vector<double> cross(vector<double> a, vector<double> b){
    vector<double> r(3,0);
    r.at(0)=a.at(1)*b.at(2)-b.at(1)*a.at(2);
    r.at(1)=a.at(2)*b.at(0)-b.at(2)*a.at(0);
    r.at(2)=a.at(0)*b.at(1)-b.at(0)*a.at(1);
    return r;
}
//function that obtains an orthogonal unitary vector to B
vector<double> ortho(vector<double> B){
    vector<double> b_unit, L, T, a={0,0,1};
    b_unit = MS(1/sqrt(inner(B,B)),B);
    double dot = inner(a,b_unit);
    L = MS(-dot,b_unit);
    T = add_M(a,L);
    return MS(1/sqrt(inner(T,T)),T);
}
//Function that generates a hamiltonian given a field.
vector<complex<double>> Hamil(vector<complex<double>> B){
    //Spin matrices of S=1
    complex<double> i(0,1);
    vector<complex<double>> Sx={0,1/sqrt(2),0,1/sqrt(2),0,1/sqrt(2),0,1/sqrt(2),0},Sy={0,-i/sqrt(2),0,i/sqrt(2),0,-i/sqrt(2),0,i/sqrt(2),0},Sz={1,0,0,0,0,0,0,0,-1},H;
    H=add_M(MS(mue*B.at(0),Sx),MS(mue*B.at(1),Sy)); H=add_M(H,MS(mue*B.at(2),Sz));
    return H;
}
//Expectation value
complex<double> obs(vector<complex<double>> O, vector<complex<double>> v){
    return inner(v,MM(O,v));
}
//Spin projection in the direction of the field
vector<complex<double>> Sp1(vector<complex<double>> b){
    complex<double> i(0,1);
    vector<complex<double>> Sx={0,1/sqrt(2),0,1/sqrt(2),0,1/sqrt(2),0,1/sqrt(2),0},Sy={0,-i/sqrt(2),0,i/sqrt(2),0,-i/sqrt(2),0,i/sqrt(2),0},Sz={1,0,0,0,0,0,0,0,-1}, Sp;
    Sp = add_M(MS(b.at(0),Sx),MS(b.at(1),Sy)); Sp=add_M(Sp,MS(b.at(2),Sz));
    return Sp;
    
}
//Time varyiing field
vector<complex<double>> B_field(vector<double> x, double Bc){
    vector<complex<double>> r={0,0,0};
    r.at(0)=Bc*(-80*pow(x.at(2),3)*x.at(0)+60*pow(x.at(0),3)*x.at(2)+60*x.at(0)*pow(x.at(1),2)*x.at(2));
    r.at(1)=Bc*(-80*pow(x.at(2),3)*x.at(1)+60*pow(x.at(1),3)*x.at(2)+60*x.at(1)*pow(x.at(0),2)*x.at(2));
    r.at(2)=Bc*(40*pow(x.at(2),4)-120*pow(x.at(0)*x.at(2),2)-120*pow(x.at(1)*x.at(2),2)+15*pow(x.at(0),4)+15*pow(x.at(1),4)+30*pow(x.at(0)*x.at(1),2));
    return r;
}
//-------------------------------------------Time Steps-----------------------------------------------------
vector<complex<double>> RK(vector<complex<double>> v, vector<complex<double>> Ht, vector<complex<double>> Hdt, double dt){
    complex<double> D(0,-dt);   
    vector<complex<double>> Tdt=MS(D,Hdt), Tt=MS(D,Ht), T2=MS(0.5,Msq(Tdt,Tt));
    return  add_M(add_M(v,MM(Tdt,v)),MM(T2,v));
}
//Crank-Nicolson
vector<complex<double>> CN(vector<complex<double>> v, vector<complex<double>> H, double dt){
    complex<double> D(0,0.5*dt*2*PI*(1/h));
    vector<complex<double>> I={1,0,0,0,1,0,0,0,1}, A=add_M(I,MS(D,H)), phi={0,0,0}, psi=MS(-1,v);
    phi.at(1)=v.at(1)-(A.at(5)/A.at(8))*v.at(2)-(A.at(3)/A.at(0))*v.at(0);
    phi.at(1)=phi.at(1)/(A.at(4)-(A.at(3)/A.at(0))*A.at(1)-(A.at(5)/A.at(8))*A.at(7));
    phi.at(0)=(v.at(0)-A.at(1)*phi.at(1))/(A.at(0));
    phi.at(2)=(v.at(2)-A.at(7)*phi.at(1))/(A.at(8));
    return add_M(MS(2,phi),psi);
}

int cross(double Vo, ofstream &out){
    out << setprecision(9) << scientific;
    //Uniformly distributed random numbers in the interval [0,1].
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0,1.0);
    //Variables in which the random number is stored.
    //For initial position (e.i. the velocity direction)
    double cthr=0, phr=0;
    //For the impact parameter
    double alpha=0, imp_N=0;
    //Maximum Montecarlo tries:
    int N_Mc=500;
    //Integrations: one if for the montecarlo and the otherone for the approximate expression b_mx^2*PI*<P>
    double mont=0, app_i=0, mont1=0;
    //Imaginary unit
    complex<double> i(0,1);
    //Octopole and linear amplitudes
    double Bc=100;
    //Position, Velocity and impact prameter vectors
    vector<double> v, xo, imp;
    //Initial B-fields
    vector<complex<double>> B, b, Sph, HCNh;
    //Initial projections
    //Sph=SpH(b);
    //B field direction and magnitude
    double BM, th, ph;
    //Landau-Zener factors
    double wN, A, C, K;
    vector<double> Bv,Bimp, uBv, CB;
    //Definition of the initial states as the up direction eigenstate 
    vector<complex<double>> psihCN, psiCN2, psiCN3, psihf, psihf2, psihf3;
    //Time constants
    double Rp=1e-3, Tf, dt, b_max=0.15e-3*Vo, l, uu=0, ud=0, uz=0, dd=0, dz=0, zz=0;
    vector<double> e1, e2, temp;
    int k; 
    for(int m=0; m<N_Mc; m++){
        //Generation of the random position parameters
        //Random polar angle
        cthr = 2*dis(gen)-1;
        //Random azimuthal angle
        phr = 2*PI*dis(gen);
        //Randon cylindrical angle for the impact parameter
        alpha = 2*PI*dis(gen);
        //Random magnitude of the impact parameter
        imp_N = b_max*sqrt(dis(gen));
        //Direction oposed to the velocity.n);
        temp = {-sqrt(1-cthr*cthr)*cos(phr),-sqrt(1-cthr*cthr)*sin(phr),-cthr};
        //Orthogonal vectors that form a basis with the direction of movement as the z direction.
        e1 = ortho(temp);
        e2 = cross(e1,temp);
        //Initial position
        xo = MS(1.5*Rp,temp);
        imp = add_M(MS(imp_N*cos(alpha),e1),MS(imp_N*sin(alpha),e2));
        xo = add_M(xo,imp);
        //Velocity
        v = MS(-Vo,temp);
        //Final time and time-step.
        Tf=Rp/Vo;
        //Initial field and unit vector in the field direction
        B=B_field(xo,Bc);
        b=MS(1/sqrt(real(inner(B,B))),B);
        dt=Tf/3e4;/*(0.1)*(h/(sqrt(real(inner(B,B)))*mue*2*PI))*/;
        k= -2*Tf/dt;
        //Instantaneous spin projections
        //Sph=Sp1(b);
        //Spherical coordinates of the field.
        BM=real(sqrt(inner(B,B)));
        th=acos(real(B.at(2))/BM);
        if(real(B.at(1))>=0){ph=acos(real(B.at(0))/(BM*sin(th)));}
        else{ph=2*PI-acos(real(B.at(0))/(BM*sin(th)));}
        //Initial states. We start the system in the up state of the local initial field. 
        psihCN={0.5*(1+cos(th))*exp(-i*ph),(1/sqrt(2))*sin(th),0.5*(1-cos(th))*exp(i*ph)};
        psiCN2 = {-(1/sqrt(2))*sin(th)*exp(-i*ph),cos(th),(1/sqrt(2))*sin(th)*exp(i*ph)};
        psiCN3 = {0.5*(1-cos(th))*exp(-i*ph),-(1/sqrt(2))*sin(th),0.5*(1+cos(th))*exp(i*ph)};
        //Time evolution:
        while(k<=2*Tf/dt){
            //Position update
            //if((k%1000)%150==0){pos <<k*dt<<" "<< xo.at(0)*1e3 <<" "<< xo.at(1)*1e3<<" "<< xo.at(2)*1e3<< endl;}
            xo.at(0)+=v.at(0)*dt;
            xo.at(1)+=v.at(1)*dt;
            xo.at(2)+=v.at(2)*dt;
            //Speed update
            //if((k%1000)%150==0){speed <<k*dt<<" "<< sqrt(real(inner(v,v))) << endl;}
            l=Vo*0.5*(1+erf(5e5*(k*dt+1.5*Tf)))-Vo*0.5*(1+erf(5e5*(k*dt-1.5*Tf)));
            v=MS(-l , temp);
            //Field Update
            //if((k%1000)%150==0){field <<k*dt<<" "<< real(B.at(0)) <<" "<< real(B.at(1))<<" "<< real(B.at(2))<< endl;}
            B=B_field(xo,Bc);
            //b=MS(1/sqrt(real(inner(B,B))),B);
            //Update of the projections
            //Sph=SpH(b);
            //Update of the Hamiltonians
            HCNh=Hamil(B);
            //Time evolution of the states
            //if((k%1000)%150==0){temp_f <<k*dt<<" "<< real(obs(Sph,psihCN))<< endl;}
            psihCN=CN(psihCN,HCNh,dt);
            psiCN2=CN(psiCN2,HCNh,dt);
            psiCN3=CN(psiCN3,HCNh,dt);
            //Counter update.
            k+=1;
        }
        //Here we obatin the spherical coordinates for the final magnetic field 
        BM=real(sqrt(inner(B,B)));
        th=acos(real(B.at(2))/BM);
        if(real(B.at(1))>=0){ph=acos(real(B.at(0))/(BM*sin(th)));}
        else{ph=2*PI-acos(real(B.at(0))/(BM*sin(th)));}
        //Here we find the final down eigenstates
        psihf={0.5*(1-cos(th))*exp(-i*ph),-(1/sqrt(2))*sin(th),0.5*(1+cos(th))*exp(i*ph)}; 
        psihf2 = {-(1/sqrt(2))*sin(th)*exp(-i*ph),cos(th),(1/sqrt(2))*sin(th)*exp(i*ph)};
        psihf3 = {0.5*(1+cos(th))*exp(-i*ph),(1/sqrt(2))*sin(th),0.5*(1-cos(th))*exp(i*ph)};
        ud+=b_max*b_max*PI*norm(inner(psihCN,psihf));
        uz+=b_max*b_max*PI*norm(inner(psihCN,psihf2));
        dz+=b_max*b_max*PI*norm(inner(psiCN2,psihf));

    }
    
    out <<Vo<<" "<< ud/N_Mc <<" "<< uz/N_Mc <<" "<< dz/N_Mc << endl;
    return 0;
}
int main(){
    double dVo=0.7;
    ofstream out;
    out.open("out.txt"); out << setprecision(9) << scientific;
    for(int i=1; i<2; i++){
        cout << i << endl;
        cross(i*dVo,out);
    }
    return 0;
}