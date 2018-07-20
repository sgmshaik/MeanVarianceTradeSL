/* 
 * File:   newmain.cpp
 * Author: mousa
 *
 * Created on 22 March 2011, 00:13
 */

#include "matrices.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <sstream>
#include <ctime>
//#include <boost/timer.hpp>
//#include <boost/progress.hpp>
#include <string>
#include<cstring>
#include<time.h>
#include <stdlib.h>

/*
 * 
 */



using namespace matrices;
using namespace std;

double point(double istar, int maxsize, int minsize)
{
    int intistar = floor(istar);
    double value;
    value= max(minsize,min(intistar,maxsize-1));
    return value;
}

inline double f(double v,double kt,double ks ,double beta)
{
    //cout << v << endl;
    double tv =0;
    if (v == 0 )
        tv = 0;
    else if(v < 0)
        tv = -1;
    else if (v > 0)
        tv = 1;

     double tempf = (1 + ks*(tv))*exp(kt*tv*pow(abs(v), beta));
   
    return tempf;
   
}

double g(double v, double kp)
{
    return (kp*v);

}


inline double interpolate2d(double hS, double halpha, double ds, double da, const matrix &vec, int SMAX, int AMAX)
{

    int s1;
    int s2;
    int a1;
    int a2;
    int spoint1;
    int spoint2;

    double ds1;
    double da1;


    ds1 = hS / ds;
    da1 = halpha / da;


    if(ds1 < -SMAX/2.)
        cout << "error  "  << ds1 << "  "<< endl;
    else if(ds1 > SMAX/2.)
        cout << "error  "  << ds1 << "  "<< endl;

    s1 = point(ds1,SMAX/2.,-SMAX/2.);

    if(s1 >= 0)
    {
    s2 = s1 +1 ;

    spoint1 = point(ds1,SMAX/2.,-SMAX/2.)+SMAX/2.;
    spoint2 = spoint1 + 1;
    }
    else if(s1 < 0)
    {

        s2 = s1 + 1;

        spoint1 = point(ds1,SMAX/2.,-SMAX/2.)+SMAX/2.;

        spoint2 = spoint1 + 1;

    }



    a1 = point(da1,AMAX,0);

    a2 = a1 + 1;


    double step1 = ((s2 * ds - hS) / double(s2 * ds - s1 * ds)) * vec[a1].at(spoint1) + ((hS - s1 * ds) / double(s2 * ds - s1 * ds)) * vec[a1].at(spoint2);
    double step2 = ((s2 * ds - hS) / double(s2 * ds - s1 * ds)) * vec[a2].at(spoint1) + ((hS - s1 * ds) / double(s2 * ds - s1 * ds)) * vec[a2].at(spoint2);
    double Val   =  ((a2 * da - halpha) / double(a2 * da - a1 * da)) * step1 + ((halpha - a1 * da) / double(a2 * da - a1 * da)) * step2;
    return Val;

}

array triSolver(matrix &mat)
{
    double a1 = 0;
    double b1 = 0; // these are current row coeffs
    double c1 = 0;
    double d1 = 0;

    double a2 = 0; // these are next row coeffs
    double b2 = 0;
    double c2 = 0;
    double d2 = 0;

    array solution(mat.size()); // since the matrix is from i = 1 to n-1

    for (int i = 0; i != (mat.size() - 1)/2.; i++)
    {

        a1 = mat[i][0];
        b1 = mat[i][1];
        c1 = mat[i][2];
        d1 = mat[i][3];
        a2 = mat[i+1][0];
        b2 = mat[i+1][1];
        c2 = mat[i+1][2];
        d2 = mat[i+1][3];

        if (i == 0)
        {
            mat[i+1][0] = 0;
            mat[i+1][1] = b2; //b1 = 0
            mat[i+1][2] = c2; //c1 = 0
            mat[i+1][3] = d2 - (a2*d1/a1);

        }
        else if (i > 0 && i < (mat.size() - 1)/2.  -1)
        {

            mat[i+1][0] = 0;
            mat[i+1][1] = (b2) - (a2*c1)/b1;
            mat[i+1][2] = (c2);
            mat[i+1][3] = d2 - (a2*d1/b1);
        }


    }
    solution[(mat.size() - 1)/2.] = mat[(mat.size()-1)/2.][3]/mat[(mat.size()-1)/2.][2];

    for (int m = (mat.size() - 1)/2. -1. ; m >= 1; m--)
    {


        solution[m] = (1./mat[m][1])*(mat[m][3] - mat[m][2]*solution[m+1]);
        //cout << " solution " << solution[m] << endl ;

    }

    solution[0] = (mat[0][3])/(mat[0][0]);


    for ( int i = (mat.size()-1)/2.; i != mat.size() - 1; i++ )
    {

        a1 = mat[i][0];
        b1 = mat[i][1];
        c1 = mat[i][2];
        d1 = mat[i][3];
        a2 = mat[i+1][0];
        b2 = mat[i+1][1];
        c2 = mat[i+1][2];
        d2 = mat[i+1][3];

        if (i == (mat.size()-1)/2. )
        {
            mat[i+1][0] = 0;
            mat[i+1][1] = b2; //b1 = 0
            mat[i+1][2] = c2; //c1 = 0
            mat[i+1][3] = d2 -a2*solution[(mat.size()-1)/2.];
        }

        else if (i > (mat.size()-1)/2. && i < (mat.size() - 2) )
        {
            mat[i+1][0] = 0;
            mat[i+1][1] = (b2) - (a2*c1)/b1;
            mat[i+1][2] = (c2);
            mat[i+1][3] = d2 - (a2*d1/b1);
        }

        }

    solution[mat.size()-1] = mat[(mat.size()-1)][3]/mat[(mat.size()-1)][2];
    for (int m = (mat.size() - 2)  ; m >= (mat.size()-1)/2. + 1; m--)
    {
        solution[m] = (1./mat[m][1])*(mat[m][3] - mat[m][2]*solution[m+1]);
        //cout << " solution " << solution[m] << endl ;
    }

    return (solution);
}

inline array hjbDiscretisation(int i, double dP, double dtow, double vol, double drift)
{
    array coeff(3);
    double alpha; // i-1 pivot location
    double beta; // i
    double gamma; //i+1
    double alphacentral;
    double betacentral;
    double gammacentral;
    double alphaforward;
    double betaforward;
    double gammaforward;
    double gammaback;
    double alphaback;



    alphacentral = (vol * vol)*(i * i) / (2.) - drift * i / ((2.));
    gammacentral = (vol*vol)*((i)*(i)) / (2.) + drift * i / (2.);
    alphaforward = (vol*vol)*((i)*(i)/2.);
    gammaforward = (vol*vol) * i * i / (2.) + drift*i;
    alphaback = (vol*vol)*i*i/2. - drift*i;
    gammaback = (vol * vol)*((i)*(i)) / 2.;

    if (alphacentral >= 0 && gammacentral >= 0) {
        alpha = alphacentral;
        gamma = gammacentral;

        beta = (alpha + gamma);

        coeff[0] = -alpha*dtow; //;// -dtow*(vol*vol*i*i/2. - drift*i/2.);//

        coeff[1] = 1 + dtow*beta; ////1 + dtow*(vol*vol*i*i);//

        coeff[2] = -gamma*dtow; //-dtow*(vol*vol*i*i/2. + drift*i/2. );//-dtow*(vol*vol*i*i/2. + drift*i/2. );//

        return coeff;
    }
    else if (gammaforward >= 0)
    {

        alpha = alphaforward;
        gamma = gammaforward;

        beta = (alpha + gamma);

        coeff[0] = -alpha*dtow; //;// -dtow*(vol*vol*i*i/2. - drift*i/2.);//

        coeff[1] = +1 + dtow*beta; ////1 + dtow*(vol*vol*i*i);//

        coeff[2] = -gamma*dtow; //-dtow*(vol*vol*i*i/2. + drift*i/2. );//-dtow*(vol*vol*i*i/2. + drift*i/2. );//


        return coeff;
    }
    else if (alphaback >= 0)
    {
        alpha = alphaback;
        gamma = gammaback;

        beta = (alpha + gamma);

        coeff[0] = -alpha*dtow; //;// -dtow*(vol*vol*i*i/2. - drift*i/2.);//

        coeff[1] = 1 + dtow*beta; ////1 + dtow*(vol*vol*i*i);//

        coeff[2] = -gamma*dtow; //-dtow*(vol*vol*i*i/2. + drift*i/2. );//-dtow*(vol*vol*i*i/2. + drift*i/2. );//


        return coeff;
    }


}

inline double limitv(double v,double dtow, double r,double kp)
{
    if(r-g(v,kp)!=0)
    {
    return (exp(r*dtow)-exp(g(v,kp)*dtow))/(r-g(v,kp));

    }
    else
        return dtow;
}

inline double characteristicsS(double HATS,double v,double dtow, double kp, double kt,double ks,double beta,double r)
{

    if (r-g(v,kp)!= 0)
    {
     return HATS*exp(g(v,kp)*dtow)/(exp(r*dtow)-(v*HATS*f(v,kt,ks,beta))*(exp(r*dtow)-exp(g(v,kp)*dtow))/(r-g(v,kp))) ;
    }
    else
    {
     return HATS*exp(g(v, kp)*dtow)/(exp(r*dtow)-(v*HATS*f(v,kt,ks,beta))*(dtow)) ;
    }

}
inline double coeff(double HATS,double v, double r, double dtow,double kp,double kt ,double ks ,double beta)
{

        double val = pow(exp(r*dtow) - (v*HATS*f(v,kt,ks,beta)*(limitv(v,dtow,r,kp))),2.);

        return(val);
}

inline double minControl(double valpha, double usermin)
{
    double val;
    if(usermin < valpha )
    {
       val =  valpha  ;
    }else
    {
        val = usermin;
    }
    return val;
}

matrix semiLagrangian1(int gS, int gA, int gT, int gC,double MAXS,double MAXA,double T,double dT,double vmin,matrix &optControl,matrix &expReturn, matrix &initial, matrix3d &objectivefunction, double kp, double kt ,double ks , double r,double vol,double mu,double beta)
{

    double dHATS = 2*(MAXS)/(double)gS;

    double da = MAXA/(double)gA;

    double lvmin = vmin;

    double vmax = 0;

    array unknown(3);

    array maxim(0);

    array control(0);

    matrix tri(gS+1,array(4));

    matrix testTriError(gS+1,array(4));

    matrix solution(gA+1,array((gS+1)));

    array solutionerror(gS+1);

    matrix preSolution(gA+1,array((gS+1))); //solution at previous timestep

    matrix opt_control(gA+1,array((gS+1)));

    matrix localobjective;

    array co_ord(2);

     double dtow = T/double(gT);
    for(int k = 0; k <= gA; k++)
    {
        double vT = -da*k/dT;
        int spoint;
        for(int i = 0; i <= gS; i++)
        {
            spoint = (i-gS/2.);
            solution[k][i] = 1.; // + 2*k*da*f(vT)*((dHATS*spoint)) +  pow(k*da*f(vT)*((dHATS*spoint)),2.);
        }
    // cout << "f  " << f(vT)<<endl ;
    }

    initial = solution;
    preSolution = solution;

    for(int n = 0; n !=4; n++)
    {
       //dtow = deltat[n];
        
            int spoint;
                for(int i = 0; i<=gS; i++)
                {
                    spoint = i - gS/2.;
                 
                    if(i == 0) //s=SMIN
                    {
                        tri[i][0] = 1.;
                        tri[i][1] = 0;
                        tri[i][2] = 0;

                        tri[i][3] = pow(exp(r*dtow),2.)*preSolution[0][0];

                        opt_control[0][i] = 0;

                        continue;
                    }

                    if(i!=0 && i < gS/2.)
                    {
                        unknown = hjbDiscretisation(spoint,dHATS,dtow, vol, mu);
                        tri[i][0]=  unknown[0];
                        tri[i][1] = unknown[1];
                        tri[i][2] = unknown[2];
                        tri[i][3]= pow(exp(r*dtow),2.)*preSolution[0][i]; //double hS,double halpha,double ds,double da,double smax,double alphamax,const matrix &vec

                        opt_control[0][i] = 0;
                  
                        continue;
                    }

                     if(i == gS/2.) //s=smax
                     {
                        tri[i][0] = 0;
                        tri[i][1] = 0;
                        tri[i][2] = 1.;
                        tri[i][3] =  pow(exp(r*dtow),2.)*preSolution[0][i];

                        opt_control[0][i]=0;
                        continue;
                      }
                    if(i!=gS)
                    {
                        unknown = hjbDiscretisation(spoint,dHATS,dtow, vol, mu);
                        tri[i][0]=  unknown[0];
                        tri[i][1] = unknown[1];
                        tri[i][2] = unknown[2];
                        tri[i][3]= pow(exp(r*dtow),2.)*preSolution[0][i]; //double hS,double halpha,double ds,double da,double smax,double alphamax,const matrix &vec
                        
                        opt_control[0][i] = 0;

                        continue;
                    }
                    if(i==gS)
                    {
                        tri[i][0] = 0;
                        tri[i][1] = 0;
                        tri[i][2] = 1.;
                        tri[i][3] =  pow(exp(r*dtow),2.)*preSolution[0][i];

                        opt_control[0][i]=0;
                    }

                }


                solution[0] = triSolver(tri);


//     #pragma omp parallel default(none) firstprivate(tri,maxim,control,unknown,lvmin,da,dHATS,gA,gC,MAXS,gS,dtow,preSolution, vol, r, mu,vmin,MAXA,beta,kp,ks,kt) shared(solution,opt_control)
       {
//     #pragma omp for nowait schedule(static)
            for(int k = 1; k<=gA; k++)
            {

              // cout << " k  "  << k << endl;
                maxim.clear();
                control.clear();
             //   cout << n  << " " << 1/dtow << " MAX "<<endl;
                int spoint;

                lvmin =  minControl (-da*k/dtow, vmin) ;
                double  dc= (lvmin-0)/(double)gC;
                for(int i = 0; i<=gS; i++)
                {
                     
                   
                    spoint = i - (gS)/2.;
                    maxim.clear();
                    control.clear();
                    if(i == 0) //s=0
                    {
                        tri[i][0] = 1.;
                        tri[i][1] = 0;
                        tri[i][2]=  0;
                        for (int c=0;c<=gC;c++)
                        {

                           double v = c*dc ; //c*dc;

                           double ha = da*k + v*dtow;

                           double hs =spoint*dHATS;

                           if(ha < 0. )
                           {
                              dc = (c+1)*dc/double(gC);
                               continue;
                           }

                           //#pragma omp flush(preSolution)
                            maxim.push_back((coeff(spoint*dHATS,v,r,dtow,kp,kt,ks,beta))*interpolate2d(hs,ha,dHATS,da,preSolution,gS,gA));
                           //#pragma omp flush(preSolution)


                            control.push_back(v);

                        }


                        tri[i][3] = *(min_element(maxim.begin(),maxim.end()));
                        int pos = (min_element(maxim.begin(),maxim.end())) -maxim.begin();

                        opt_control[k][i] =control[pos];

                        double opt_ha = da*k + opt_control[k][i]*dtow;
                       // cout << opt_ha << " opt_ha " << endl;
                        continue;
                    }

                    if(i!=0 && i <gS/2.)
                    {
                        unknown = hjbDiscretisation(spoint,dHATS,dtow,vol,mu);
                        tri[i][0] = unknown[0];
                        tri[i][1] = unknown[1];
                        tri[i][2] = unknown[2];

                   //     controlgen =  controlGenerator(MAXS, -MAXS,spoint*dHATS, dtow,da*k, gC,vmin);

                        for (int c=0;c<=gC;c++)
                        {
                       //   cout << "c "  <<controlgen[i][k].size()<< endl;

                             //#pragma omp flush(preSolution)
                            double v = c*dc;
                           // cout << v;
//#pragma omp flush(preSolution)
                            double ha = da*k + v*dtow;
                            double hs = characteristicsS(spoint*dHATS,v,dtow,  kp, kt , ks,beta,r);

                              if(ha < 0. ||hs < -MAXS || hs > MAXS)
                            {

                              
                                continue;
                            }

                            //#pragma omp flush(preSolution)
                            maxim.push_back((coeff(spoint*dHATS,v,r,dtow,kp,kt,ks,beta))*interpolate2d(hs,ha,dHATS,da,preSolution,gS,gA));
                            //#pragma omp flush(preSolution)
                            control.push_back(v);

                        }


                        tri[i][3] = *(min_element(maxim.begin(),maxim.end()));
                        int pos = min_element(maxim.begin(),maxim.end()) - maxim.begin();

                          opt_control[k][i] = control[pos];
           
                          double opt_v = control[pos];

                          double opt_hs = characteristicsS(spoint*dHATS,opt_v,dtow,  kp,  kt,ks,beta,r);

                          double opt_ha = da*k + opt_v*dtow;

                          continue;
                    }

                   if(i == gS/2.) //s=smax
                   {
                        tri[i][0] = 0;
                        tri[i][1] = 0;
                        tri[i][2] = 1.;

                     
                        for(int c = 0; c<=gC; c++)
                        {
                           double v = c*dc;

                            double ha = da*k +v*dtow;

                              if(ha < 0.)
                              {
                                   dc = (c+1)*dc/double(gC);

                                continue;
                              }
                               //#pragma omp flush(preSolution)
                            maxim.push_back(pow(exp(r*dtow),2.)*interpolate2d(0,ha,dHATS,da,preSolution,gS,gA));
                             //#pragma omp flush(preSolution)
                            control.push_back(v);



                        }

                               // maxim.push_back(pow((exp(r*dtow) - (v*spoint*dHATS*f(v,kt,ks,beta)))*((exp(r*dtow) - exp(g(v,kp)*dtow)))/(r-g(v,kp)),2.)*interpolate2d(spoint*dHATS,ha,dHATS,da,MAXS,MAXA,preSolution));



                        tri[i][3] = *(min_element(maxim.begin(),maxim.end()));

                        int pos = (min_element(maxim.begin(),maxim.end())) -maxim.begin();
                        opt_control[k][i] = control[pos];

                        double opt_v = control[pos];
                        double opt_ha = da*k + opt_v*dtow;

                        continue;
                   }

                    if(i!=gS)
                    {
                        unknown = hjbDiscretisation(spoint,dHATS,dtow,vol,mu);
                        tri[i][0] = unknown[0];
                        tri[i][1] = unknown[1];
                        tri[i][2] = unknown[2];
                        //controlgen =  controlGenerator(MAXS, -MAXS,spoint*dHATS, dtow,da*k, gC, vmin);
                        for (int c=0;c<=gC;c++)
                        {
                              //#pragma omp flush(preSolution)
                             double v =c*dc ;
                               //#pragma omp flush(preSolution)

                            double ha = da*k + v*dtow;
                            double hs =  characteristicsS(spoint*dHATS, v,dtow, kp, kt , ks,beta,r);
//(r-g(v,kp))/((-1*v*f(v,kt,ks,beta)+((r-g(v,kp))/(spoint*dHATS)))*exp(dtow*(r-g(v,kp))) + v*f(v,kt,ks,beta));

                            
                            if(ha < 0. ||hs < -MAXS || hs > MAXS)
                            {

                                dc = (c+1)*dc/double(gC);
                                continue;
                            }

                            maxim.push_back((coeff(spoint*dHATS,v,r,dtow,kp,kt,ks,beta))*interpolate2d(hs,ha,dHATS,da,preSolution,gS,gA));
                               //#pragma omp flush(preSolution)
                            control.push_back(v);


                        }

                        tri[i][3] = *(min_element(maxim.begin(),maxim.end()));

                        int pos = min_element(maxim.begin(),maxim.end()) - maxim.begin();

                          opt_control[k][i] = control[pos];
                //          cout << "S != 0" <<  opt_control[k][i];
                          double opt_v = control[pos];

                          double opt_hs =characteristicsS(spoint*dHATS,opt_v,dtow,  kp, kt, ks,beta,r);
                          double opt_ha = da*k + opt_v*dtow;

                          continue;
                    }

                    if(i==gS) //s=smax
                    {
                        tri[i][0] = 0;
                        tri[i][1] = 0;
                        tri[i][2] = 1. ;
                        for(int c = 0; c<=gC; c++)
                        {
                            double v = c*dc;

                            double ha = da*k +v*dtow;
                             double hs = spoint*dHATS;
                            
                            if(ha < 0.)
                            {
                         
                                dc = (c+1)*dc/double(gC);
                                continue;
                            }

                                //#pragma omp flush(preSolution)
                            maxim.push_back((coeff(spoint*dHATS,v,r,dtow,kp,kt,ks,beta))*interpolate2d(hs,ha,dHATS,da,preSolution,gS,gA));

                            control.push_back(v);


                        }


                        tri[i][3] = *(min_element(maxim.begin(),maxim.end()));

                        int pos = (min_element(maxim.begin(),maxim.end())) - maxim.begin();

                        opt_control[k][i] = control[pos];

                        double opt_v = control[pos];

                        double opt_ha = da*k + opt_v*dtow;
                   }

                }
                    solution[k] = triSolver(tri);

            }

        }//pragma


        preSolution  = solution;

    }// n


    optControl = opt_control;

     cout <<  " time   " << gT*dtow  <<  "value   " <<interpolate2d(-1.,1.,dHATS,da,solution,gS,gA) << endl;

    cout <<  " time   " << gT*dtow  <<  "control value   " <<interpolate2d(-1,1.,dHATS,da,opt_control,gS,gA) << endl;

        return solution;
}



void grid(int i ,int &gS, int &gA, int &gT, int &gC ,int s, int a ,int t,int c)
{
            gS = (i==0)? (gS) : gS*s ;

            gA= (i==0) ? (gA):  gA*a ; // why wont it work with 161;

            gT= (i==0) ? gT : gT*t;
            gC= (i==0) ? gC : gC*c;
}

void printFileVhatSA(const matrix &results, double ds, double da, char refine[], double MAXS, double vol, double time, double vmin,char name [],double mu,double r, double kt, double ks ,double kp, double beta)
{
//const matrix &results, double ds, double da,double lambda,char refine[],double MAXS,double vol,double time,double vmin,char name []
    char rname [400];
    sprintf(rname,"./resultVSA/VSArefine%sSmax%fvol%ftime%fvmin%fmu%fr%fkt%fks%fkp%fb%f%s.dat" , refine, MAXS,vol,time,vmin,mu, r,kt,ks, kp,beta,name);
    ofstream dataFile(rname, ios::out);
    dataFile.precision(10);
    cout.precision(10);

    if (!dataFile) {
        cout << " Error " << endl;
    }

    dataFile << "# " << setw(25)<<  "Alpha" << setw(25) << "S" << setw(20) << "V" << '\n' << endl;
    //cout <<"# " << setw(4) << "Alpha"<< setw(15) << "S" << setw(15) <<"V"  <<'\n' <<endl;

    for (int k = 0; k <= results.size() - 1; k++)
    {
        int spoint ;
        for (int i = 0; i <= results[k].size() - 1; i++)
        {
            spoint = i - (results[k].size() - 1)/2.;
            dataFile << showpoint << setw(20)  << k * da << setw(20) << ds * spoint << setw(20) << results[k][i] << endl;
            // cout <<showpoint <<setw(6) << k*da <<setw(15) << ds*(i)<< setw(15)  << results[k][i] << endl;
        }
        dataFile << endl;
        //cout << endl;
    }
    dataFile << endl;
    dataFile << endl;
    //dataFile << "end" << endl;
    //dataFile << "pause -1" << endl;
    dataFile << endl;
    dataFile.close();
}
void printConvergence(const matrix &results,double T,double vmin,double SMAX, double SMIN, char  name [],double vol, double r, double mu, double beta,double kp,double ks, double kt)
{
    char rname [800];
    sprintf(rname, "./convergence/%sT%fvmin%fSMAX%fSMIN%fvol%fr%fmu%fbeta%fkp%fks%fkt%f.dat",name,T,vmin, SMAX, SMIN ,vol, r, mu,beta, kp, ks , kt);
    ofstream dataFile(rname, ios::out);
    dataFile.precision(10);
    cout.precision(10);

    if (!dataFile)
    {
        cout << " Error " << endl;
    }

    dataFile << " #Refine  " << setw(6)  << " [gC]  " << setw(6)<< "  [gS]  " << setw(6) << "  [gA]  " << setw(6) << "  [gT]  "<< setw(15) << " E[(B_L)^2] " << setw(20) << " abs[Refine_i - Refine_(i-1)] "  << setw(20)<< "opt_v" << setw(20) << "duration" <<endl;

            for (int i = 0; i <= results.size() - 1; i++)
            {
                dataFile <<setw(9)  <<i << setw(6)  << results[i][0] << setw(6)<< results[i][1] << setw(6) << results[i][2] << setw(6) << results[i][3]<< setw(15) << results[i][4] << setw(20) <<results[i][5] << setw(20) <<results[i][6]<< setw(20) << results[i][7]  << endl;
            // cout <<showpoint <<setw(6) << k*da <<setw(15) << ds*(i)<< setw(15)  << results[k][i] << endl;
            }

            dataFile << endl;

    dataFile << endl;
    dataFile << endl;
    //dataFile << "end" << endl;
    //dataFile << "pause -1" << endl;
    dataFile << endl;
    dataFile.close();
}

void resultsrefinedcontrol(void  (*grid)(int i, int & , int &, int &,int &, int ,int a,int t, int c),int gS, int gA,int gT, int gC,double vol, double r, double mu,double ks, double kt,double kp,double beta ,double MAXS ,double T,double dT,int irefine,int s,int a, int t,int c,double maxvmin)
{

    matrix solution;
    matrix optControl;
    matrix initial;
    matrix expReturn;
    matrix obj_func;
    matrix characteristichats;
    matrix solutionprerefine;
    matrix optControlprerefine;
    matrix value(irefine+2,array(8));
    matrix differenceValue;
    matrix differenceControl;
    matrix3d objectivefunction(gS+1);

    int lgS,lgA, lgT, lgC;

    lgS = gS;
    lgA = gA;
    lgT = gT;
    lgC = gC;

    double vmin;

    for(int i =0; i <=irefine; i++)
    {
      
         (*grid)(i, lgS,lgA,lgT,lgC,s,a,t,c);

         objectivefunction.resize(lgS+1);

         double  MAXA= 1.;

         double dtow = T/(double(lgT));

         double vmin = -1/dtow;

         double ds = 2*MAXS/(double)lgS;

         double da = 1./(double)lgA;

         stringstream refinestream (stringstream::out);

         refinestream <<"refined" <<"GS" << lgS << "GA" << lgA << "GT" << lgT<< "GC"<<lgC ;

         string strrefine = refinestream.str();
          //cout << refinestream << endl;g

         char refine [strrefine.size()+1];

         strcpy(refine,strrefine.c_str());
         time_t start;
         time_t end;

         time(&start);// semiLagrangian1(int gS, int gA, int gT, int gC,double MAXS,double MAXA,double T,double dT,double vmin,matrix &optControl,matrix &expReturn, matrix &initial, matrix3d &objectivefunction, double kp, double kt ,double ks , double r,double vol,double mu,double beta)
         solution = semiLagrangian1(lgS, lgA, lgT,lgC, MAXS , MAXA , T, dT ,maxvmin ,optControl , expReturn , initial,objectivefunction, kp,  kt , ks ,  r, vol, mu,beta);
         time(&end);

         double duration = difftime(end,start);

         if(maxvmin <= vmin)
            {
                vmin=  vmin  ;
            }else
            {
               vmin = maxvmin;
            }

         value[0][3] = 0;
         value[0][0] = 0;
         value[0][1] = 0;
         value[0][2] = 0;
         value[0][3] = 0;
         value[0][4] = 0;
         value[0][5] = 0;
         value[i+1][0] = lgC;
         value[i+1][1] = lgS;
         value[i+1][2] = lgA;
         value[i+1][3] = lgT;
         value[i+1][4] = interpolate2d(-1.,1.,ds,da,solution,lgS,lgA);
         value[i+1][5] = abs(value[i+1][4] - value[i][4]);
         value[i+1][6] = interpolate2d(-1.,1.,ds,da,optControl,lgS,lgA);
         value[i+1][7] = duration;
         cout <<" [gC]  "<< lgC << "  [gS]  " << lgS << "  [gA]  " << lgA << "  [gT]   " << lgT << " [Solution] " << value[i+1][4]  << "  Convergence  " <<  value[i+1][5]  << "opt_v"<< value[i+1][6] << "duration" << value[i+1][7]  <<endl;

        char name[] = "resultreducedrefined";


        matrix vsaresult(lgA+1,array(lgS+1));

       

         char name1[] = "initrefined";

     
         char name2[]  = "origresultrefined";

        printFileVhatSA(solution, ds, da, refine, MAXS, vol, T, vmin,name2, mu,r, kt, ks ,kp,  beta);

        char name3[] = "controlrefined";

   
        char name4[] = "objectivefunctionrefined";


        objectivefunction.clear();


  char name6[] = "differenceValue";
  if(i > 0)
      {

//      gridCheck(solution,solutionprerefine, differenceValue);
   //   printFilegridCheck(differenceValue,refine, da ,ds ,MAXS,  vol, T,  vmin , name6);
      }
  solutionprerefine = solution;



  char name7[] = "differenceControl";
  if(i >0)
  {
 // gridCheck(optControl,optControlprerefine, differenceControl);
 // printFilegridCheck(differenceControl, refine, da ,ds,MAXS, vol, T, vmin, name7);
   }
//optControlprerefine = optControl;

   }


   char name5[] = "refinedcontroltimerParallel";
   printConvergence(value, T,maxvmin,MAXS, -MAXS,name5,vol,r,mu,beta,kp,ks,kt);

}



int resultsMain(int onoff)
{
    if (onoff == 1)
    {

     resultsrefinedcontrol(grid,250,80,50, 50, 1. , 0, 0., 0. ,2E-6, 0, 1. ,25, 1/250. , 10E-9 ,3, 2,2 , 2,2,-250*25);
     //resultsrefinedcontrol(grid,1000,200,50, 250, 1 , 0, 0., 0. ,2E-6, 0, 1 ,25, 1/250., 10E-9 ,5, 2,2 , 2,2,-250*25);
    }

    else
        return 1;

    return 1;
}

int main()
{
   // int gS, int gA,int gT, int gC,

    resultsMain(1);

     return (0);
}
