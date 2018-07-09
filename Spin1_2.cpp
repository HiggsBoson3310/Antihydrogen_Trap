#include <iostream> 
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <complex>
#include <random> //Random number library
using namespace std;
const double PI=3.14159265359, h = 6.626e-34, mue = 928.5e-26, z=2*PI*(mue/h);
//Definition of matrix multiplication in two dimensions for both complex vectors
vector<complex<double>> MMh(vector<complex<double>> M, vector<complex<double>> a){
    complex<double> r0,r1;
    r0 = M.at(0)*a.at(0)+M.at(1)*a.at(1);
    r1 = M.at(2)*a.at(0)+M.at(3)*a.at(1);
    vector <complex<double>> r={r0,r1};
    return r;
}
//Definition of inner product in the complex space and real space
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
//We define matrix multiplication between 2x2 matrices
vector<complex<double>> Msqh(vector<complex<double>> M, vector<complex<double>> A){
    vector<complex<double>> r(A.size(), 0);
    r.at(0)=M.at(0)*A.at(0)+M.at(1)*A.at(2);
    r.at(1)=M.at(0)*A.at(1)+M.at(1)*A.at(3);
    r.at(2)=M.at(2)*A.at(0)+M.at(3)*A.at(2);
    r.at(3)=M.at(2)*A.at(1)+M.at(3)*A.at(3);
    return r;
}
//Addition of vectors
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
//Multiply a vector and a scalar complex or real.
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
//Function that defines the Hamiltonian given a field
vector<complex<double>> HamilH(vector<complex<double>> B){
    //Spin matrices of S=1/2
    complex<double> i(0,1);
    vector<complex<double>> Sx={0,1,1,0}, Sy={0,-i,i,0}, Sz={1,0,0,-1}, H;
    H=add_M(MS(-mue*B.at(0),Sx),MS(-mue*B.at(1),Sy)); H=add_M(H,MS(-mue*B.at(2),Sz));
    return H;
}
//Expectation value of operator
complex<double> obs(vector<complex<double>> O, vector<complex<double>> v){
    return inner(v,MMh(O,v));
}
//Projections
vector<complex<double>> SpH(vector<complex<double>> b){
    complex<double> i(0,1);
    vector<complex<double>>  Sx={0,1,1,0}, Sy={0,-i,i,0}, Sz={1,0,0,-1}, Sp;
    Sp = add_M(MS(b.at(0),Sx),MS(b.at(1),Sy)); Sp=add_M(Sp,MS(b.at(2),Sz));
    return Sp;
    
}
//---------------------------------------------------------------Numerical Steps----------------------------------------------------------------------------------------

//Second order RK method.
vector<complex<double>> RK(vector<complex<double>> v, vector<complex<double>> Ht, vector<complex<double>> Hdt, double dt){
    complex<double> D(0,-dt);
    vector<complex<double>> Tdt=MS(D,Hdt), Tt=MS(D,Ht), T2=MS(0.5,Msqh(Tdt,Tt));
    return  add_M(add_M(v,MMh(Tdt,v)),MMh(T2,v));
}
//Crank-Nicolson
vector<complex<double>> CNh(vector<complex<double>> v, vector<complex<double>> H, double dt){
    complex<double> D(0,0.5*dt*2*PI*(1/h));
    vector<complex<double>> I={1,0,0,1}, A=add_M(I,MS(D,H)), phi={0,0}, psi=MS(-1,v);
    phi.at(1) = v.at(1)-(A.at(2)/A.at(0))*v.at(0);
    phi.at(1) = phi.at(1)/(A.at(3)-((A.at(2)/A.at(0))*A.at(1)));
    phi.at(0) = (v.at(0)-A.at(1)*phi.at(1))/A.at(0);
    return add_M(MS(2,phi),psi);

}
//New Varying field.
vector<complex<double>> B_field(vector<double> x, double Bo, double Bl){
    vector<complex<double>> r={0,0,0};
    r.at(0)=Bo*(4*pow(x.at(0),3)-12*pow(x.at(1),2)*x.at(0))-Bl*(0.5*x.at(0));
    r.at(1)=Bo*(-12*pow(x.at(0),2)*x.at(1)+4*pow(x.at(1),3))-Bl*(0.5*x.at(1));
    r.at(2)=Bl*(x.at(2));
    return r;
}


int main(){
    cout << setprecision(9) << scientific;
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
    int N_Mc=25000;
    //Integrations: one if for the montecarlo and the otherone for the approximate expression b_mx^2*PI*<P>
    double mont=0, app_i=0, mont1=0;
    //Output files
    /*ofstream  temp_f, field, pos;
    pos.open("pos.txt"), field.open("field.txt"), temp_f.open("temp.txt");
    pos<< setprecision(9) << scientific;
    field<< setprecision(9) << scientific;
    temp_f << setprecision(9) << scientific;
    ofstream prob, speed;
    prob.open("prob.txt");
    speed.open("speed.txt");
    prob << setprecision(9) << scientific; 
    speed << setprecision(9) << scientific;*/

    ofstream out_f;
    out_f.open("out.txt"); out_f << setprecision(9) << scientific;
    //Imaginary unit
    complex<double> i(0,1);
    //Octopole and linear amplitudes
    double Bo=0, Bl=1.139/2;
    //Velocity and initial position directions and magnitudes
    double Vo=70;
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
    vector<complex<double>> psihCN, psihf;
    //Time constants
    double Rp=2e-3, Tf, dt, b_max=0.1e-3, l, bm=1e6;
    vector<double> e1, e2, temp;
    int k; 
    cout << "fine" << endl;
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
        xo = MS(Rp,temp);
        imp = add_M(MS(imp_N*cos(alpha),e1),MS(imp_N*sin(alpha),e2));
        xo = add_M(xo,imp);
        //Velocity
        v = MS(-Vo,temp);
        //Final time and time-step.
        Tf=sqrt(inner(xo,xo))/Vo;
        //Initial field and unit vector in the field direction
        B=B_field(xo,Bo,Bl);  
        //b=MS(1/sqrt(real(inner(B,B))),B);
        dt=(0.1)*(h/(sqrt(real(inner(B,B)))*mue*2*PI));
        k= -2*Tf/dt;
        //Instantaneous spin projections
        //Sph=SpH(b);
        //Spherical coordinates of the field.
        BM=real(sqrt(inner(B,B)));
        th=acos(real(B.at(2))/BM);
        if(real(B.at(1))>=0){ph=acos(real(B.at(0))/(BM*sin(th)));}
        else{ph=2*PI-acos(real(B.at(0))/(BM*sin(th)));}
        //Initial states. We start the system in the up state of the local initial field. 
        psihCN={cos(th/2),exp(i*ph)*sin(th/2)};
        Bv = {-0.5*v.at(0),-0.5*v.at(1),v.at(2)};
        Bimp = {-0.5*imp.at(0),-0.5*imp.at(1), imp.at(2)};
        wN=sqrt(inner(Bv,Bv));
        uBv = MS(1/wN,Bv);
        CB = {real(B_field(imp,Bo,Bl).at(0)),real(B_field(imp,Bo,Bl).at(1)),real(B_field(imp,Bo,Bl).at(2))};
        CB = add_M(CB,MS(-inner(uBv,CB),uBv));
        A=Bl*wN;
        C=sqrt(inner(CB,CB));
        K=2*PI*mue*(1/h)*C*C*(1/A);
        //Time evolution:
        while(k<=2*Tf/dt){
            //Position update
            //if((k%1000)%150==0){pos <<k*dt<<" "<< xo.at(0)*1e3 <<" "<< xo.at(1)*1e3<<" "<< xo.at(2)*1e3<< endl;}
            xo.at(0)+=v.at(0)*dt;
            xo.at(1)+=v.at(1)*dt;
            xo.at(2)+=v.at(2)*dt;
            //Speed update
            //if((k%1000)%150==0){speed <<k*dt<<" "<< sqrt(real(inner(v,v))) << endl;}
            l=Vo*0.5*(1+erf(2e5*(k*dt+1.5*Tf)))-Vo*0.5*(1+erf(2e5*(k*dt-1.5*Tf)));
            v=MS(-l , temp);
            //Field Update
            //if((k%1000)%150==0){field <<k*dt<<" "<< real(B.at(0)) <<" "<< real(B.at(1))<<" "<< real(B.at(2))<< endl;}
            B=B_field(xo,Bo,Bl);
            //b=MS(1/sqrt(real(inner(B,B))),B);
            //Update of the projections
            //Sph=SpH(b);
            //Update of the Hamiltonians
            HCNh=HamilH(B);
            //Time evolution of the states
            //if((k%1000)%150==0){temp_f <<k*dt<<" "<< real(obs(Sph,psihCN))<< endl;}
            psihCN=CNh(psihCN,HCNh,dt);
            //Counter update.
            k+=1;
        }
        //Here we obtain the spherical coordinates for the final magnetic field 
        BM=real(sqrt(inner(B,B)));
        th=acos(real(B.at(2))/BM);
        if(real(B.at(1))>=0){ph=acos(real(B.at(0))/(BM*sin(th)));}
        else{ph=2*PI-acos(real(B.at(0))/(BM*sin(th)));}
        //Here we find the final down eigenstates
        psihf={-exp(-i*ph)*sin(th/2),cos(th/2)};
        /* 
        //Here we print the probability of a spin flip as the norm of the projection between the final state and the final down eigenstate.
        prob << imp_N*cos(alpha)*1e3 <<" "<< imp_N*sin(alpha)*1e3 << " " << norm(inner(psihCN,psihf)) << endl;
        //Finally we print the probability given by the landau-Zener formula.
        temp_f << imp_N*cos(alpha)*1e3 <<" "<< imp_N*sin(alpha)*1e3 << " " << exp(-K*PI)<< endl;*/
        
        //Output print:
        if(norm(inner(psihCN,psihf))<=exp(-1) && bm>imp_N){bm=imp_N;}
        out_f <<norm(inner(psihCN,psihf))<<" "<<exp(-K*PI)<<endl;

        
    }
    return 0;
}