#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <chrono>

double DxV(double x, double y);
double DyV(double x, double y);
int toIndex(double a, double b, double x, int N);
double test(double a, double b, double beta, std::vector<double> bias, int N, int t=0);
double IntyV(double x, double beta, double pas);
void tpsSortie();
void SchemaStoNaif();
void Schema122();
void Schema250();
int ABF_SchemaSto();
int ABF_SchemaSto_Rep();
int ABF_Schema250();
int ABF_Schema250_Rep();
double E(double qx, double qy, double px, double py);
double V(double x, double y);

int main(int argc, char *argv[])
{
    tpsSortie();


    return 0;
}

void tpsSortie(){
    std::mt19937_64 generator;
    std::normal_distribution<double> d(0.0,1.0);
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration t = myclock::now() - beginning;
    unsigned seed = t.count();
    generator.seed(seed);
/*
    std::ofstream resultat;
    resultat.open("resultatStoNaif.txt");
*/
    int N = 100;
    double Dt = 0.01;
    double a = -2;
    double b = 2;
    double beta = 5.5;
    double varminT = 0.000001;
    double varmin = 0.01;
    double var = 1;
    double Tmoyen = 0.0;
    int j = 0;
    double carresT = 0.0;
    while((((j<=2)||(var>varminT)))&&(j<10000)){
        j++;
        double tmp = Tmoyen;
        double X = -1.5;
        double Y = 0;
        /*
        std::vector<double> den(N+2);
        std::vector<double> num(N+2);
        std::vector<double> bias(N+2);
        std::vector<double> carres(N+2);

        for(int i=0;i<N+2;i++){
            num[i] = 0;
            den[i] = 0;
            bias[i] = 0;
            carres[i] = 0;
        }
*/
        while(X<0.5){
            Tmoyen++;/*
            int ind = toIndex(a,b,X,N);
            den[ind]++;
            num[ind] += DxV(X,Y);
            carres[ind] += DxV(X,Y)*DxV(X,Y);
            double m2 = (num[ind]/den[ind])*(num[ind]/den[ind]);
            double var = (carres[ind]/den[ind]-m2)/den[ind];
            if( ((var/m2)<varmin) && (ind<N) )
                bias[ind] = num[ind]/den[ind];

*/
            double x = X;
            double y = Y;

            X = x-Dt*DxV(x,y)+sqrt(2*Dt/beta)*d(generator);
            Y = y-Dt*DyV(x,y)+sqrt(2*Dt/beta)*d(generator);
        }

        double delta = Tmoyen-tmp;
        std::cout<<"Pour le tour "<<j<<" : "<<delta<<std::endl;
        carresT += delta*delta;
        var = (carresT/j-(Tmoyen/j)*(Tmoyen/j))/(j*(Tmoyen/j)*(Tmoyen/j));
    }
    std::cout<<"Temps moyen obtenu pour "<<j<<" tours : "<<Tmoyen/j<<std::endl;
    std::cout<<"Ecart type : "<<sqrt(var)<<std::endl;
}

double test(double a, double b, double beta, std::vector<double> bias, int N, int t){
/*
    std::ofstream resultat;
    resultat.open((std::string("estimationV")+std::string(std::to_string(t))+std::string(".txt")).c_str());
*/
    std::vector<double> Int(N);
    double h = (b-a)/((double) N);

    Int[0] = IntyV(a,beta,0.01) - bias[0]*a;
    for(int i=1; i<N; i++){
        Int[i] = (a+i*h)*(bias[i-1]-bias[i]) + Int[i-1];
    }

    double diffL1 = 0.0;
    for(int i=0; i<N; i++){
        diffL1 += std::abs(IntyV(a+i*h,beta,0.01)-bias[i]*(a+i*h)-Int[i]);
  //      resultat<<bias[i]*(a+i*h)+Int[i]<<" ";
    }

   // std::cout<<"Erreur L1 = "<<diffL1<<std::endl;

    //resultat.close();

    return diffL1;
}

double DxV(double x, double y){
    return (-16.0*x*(1-x*x-y*y)+8*x*(x*x-2)+4*(x+y)*((x+y)*(x+y)-1)+4*(x-y)*((x-y)*(x-y)-1))/6.0;
}

double DyV(double x, double y){
    return (-16.0*y*(1-x*x-y*y)+4*(x+y)*((x+y)*(x+y)-1)+4*(y-x)*((x-y)*(x-y)-1))/6.0;
}

void SchemaStoNaif(){

    std::mt19937_64 generator;
    std::normal_distribution<double> d(0.0,1.0);
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration t = myclock::now() - beginning;
    unsigned seed = t.count();
    generator.seed(seed);

    std::ofstream resultat;
    resultat.open("resultatStoNaif_beta6.txt");


    int T = 10000000;
    double Dt = 0.01;
    double beta = 6.0;


    double X=0;
    double Y=0;

    resultat<<X<<" ";

    for(int i=1;i<T;i++){

        double x = X;
        double y = Y;

        X = x-Dt*DxV(x,y)+sqrt(2*Dt/beta)*d(generator);
        Y = y-Dt*DyV(x,y)+sqrt(2*Dt/beta)*d(generator);

        resultat<<X<<" ";
    }

    resultat.close();

}

double E(double qx, double qy, double px, double py){

    return 0.5*(px*px+py*py)+V(qx,qy);
}

double V(double x, double y){

    return (4*(1-x*x-y*y)*(1-x*x-y*y)+2*(x*x-2)*(x*x-2)+((x+y)*(x+y)-1)*((x+y)*(x+y)-1)+((x-y)*(x-y)-1)*((x-y)*(x-y)-1))/6.0;
}

void Schema122(){

    std::mt19937_64 generator;
    std::normal_distribution<double> d(0.0,1.0);
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration t = myclock::now() - beginning;
    unsigned seed = t.count();
    generator.seed(seed);

    std::ofstream resultat;
    resultat.open("resultat122.txt");


    int T = 1000000;
    double Dt = 0.01;
    std::vector<double> PX(T);
    std::vector<double> PY(T);
    std::vector<double> QX(T);
    std::vector<double> QY(T);

    PX[0]=d(generator);
    PY[0]=d(generator);
    QX[0]=10.0;
    QY[0]=2.0;

    for(int i=1;i<T;i++){
        double PX12 = PX[i-1]-0.5*Dt*DxV(QX[i-1],QY[i-1]);
        double PY12 = PY[i-1]-0.5*Dt*DyV(QX[i-1],QY[i-1]);

        QX[i] = QX[i-1]+Dt*PX12;
        QY[i] = QY[i-1]+Dt*PY12;

        PX[i] = PX12-0.5*Dt*DxV(QX[i],QY[i]);
        PY[i] = PY12-0.5*Dt*DyV(QX[i],QY[i]);


        resultat<<QX[i]<<" ";
    }

    resultat.close();

    std::cout<< QX[std::distance(std::begin(QX), std::min_element(std::begin(QX), std::end(QX)))]<<std::endl<<QX[std::distance(std::begin(QX), std::max_element(std::begin(QX), std::end(QX)))]<<std::endl;
}

void Schema250(){
    std::mt19937_64 generator;
    std::normal_distribution<double> d(0.0,1.0);
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration t = myclock::now() - beginning;
    unsigned seed = t.count();
    generator.seed(seed);

    std::ofstream resultat;
    resultat.open("resultat250_beta6.txt");


    int T = 10000000;
    double Dt = 0.01;
    double beta = 6.0;
    double sigma = 2.0;
    double gamma = beta*sigma*sigma*0.5;
    double PX;
    double PY;
    double QX;
    double QY;

    PX=d(generator);
    PY=d(generator);
    QX=0.0;
    QY=0.0;

    for(int i=1;i<T;i++){

        double PX14 = (PX-0.25*Dt*gamma*PX+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);
        double PY14 = (PY-0.25*Dt*gamma*PY+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);

        double PX12 = PX14-0.5*Dt*DxV(QX,QY);
        double PY12 = PY14-0.5*Dt*DyV(QX,QY);

        QX = QX+Dt*PX12;
        QY = QY+Dt*PY12;

        double PX34 = PX12-0.5*Dt*DxV(QX,QY);
        double PY34 = PY12-0.5*Dt*DyV(QX,QY);

        PX = (PX34-0.25*gamma*Dt*PX34+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);
        PY = (PY34-0.25*gamma*Dt*PY34+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);


        resultat<<QX<<" ";
    }

    resultat.close();

}

int toIndex(double a, double b, double x, int N){

    if(x<a)
        return N;
    else if(x>=b)
        return N+1;
    else{
        return (int) floor(((double) N)*(x-a)/(b-a));
    }

}

double IntyV(double x, double beta, double pas){

    double y = -100.0;

    double Int = 0.0;

    while(y<=100.0){
        Int += exp(-beta*V(x,y))*pas;
        y += pas;
    }

    return -log(Int)/beta;
}

int ABF_SchemaSto(){
    std::mt19937_64 generator;
    std::normal_distribution<double> d(0.0,1.0);
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration t = myclock::now() - beginning;
    unsigned seed = t.count();
    generator.seed(seed);
/*
    std::ofstream resultat;
    resultat.open("resultatABFOverdamped_10traj.txt");
*/

    int N = 120;
    double a = -2;
    double b = 2;
    int L = 1;
    double Dt = 0.01;
    double beta = 6.0;
    double varmin = 0.01;
    int T = 0;

    std::vector<double> X(L);
    std::vector<double> Y(L);
    for(int l=0;l<L;l++){
        X[l] = 0;
        Y[l] = 0;
    }

    std::vector<double> den(N+2);
    std::vector<double> num(N+2);
    std::vector<double> bias(N+2);
    std::vector<double> carres(N+2);
    for(int i=0;i<N+2;i++){
        num[i] = 0;
        den[i] = 0;
        bias[i] = 0;
        carres[i] = 0;
    }

    while(1){
        T++;

        if(T%5000==0){
            //std::cout<<"Tour :"<<T<<std::endl;
            double tst = test(a,b,beta,bias,N);
            //std::cout<<tst<<std::endl;
            if(tst<2.5)
                return T;
        }

        for(int l=0;l<L;l++){
            int ind = toIndex(a,b,X[l],N);
            den[ind]++;
            num[ind] += DxV(X[l],Y[l]);
            carres[ind] += DxV(X[l],Y[l])*DxV(X[l],Y[l]);
            double m2 = (num[ind]/den[ind])*(num[ind]/den[ind]);
            double var = (carres[ind]/den[ind]-m2)/den[ind];
            if( ((var/m2)<varmin) && (ind<N) )
                bias[ind] = num[ind]/den[ind];
        }
        for(int l=0;l<L;l++){
            int ind = toIndex(a,b,X[l],N);

            double x = X[l];
            double y = Y[l];

            X[l] = x-Dt*DxV(x,y)+bias[ind]*Dt+sqrt(2*Dt/beta)*d(generator);
            Y[l] = y-Dt*DyV(x,y)+sqrt(2*Dt/beta)*d(generator);

//            resultat<<X[l];


/*
            if(l!=(L-1))
                resultat<<" ";
            else
                resultat<<std::endl;
  */      }
    }
    /*
    std::cout<<"Fin"<<std::endl;

    test(a,b,beta,bias,N);

  //  resultat.close();
*/}

int ABF_Schema250(){
    std::mt19937_64 generator;
    std::normal_distribution<double> d(0.0,1.0);
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration t = myclock::now() - beginning;
    unsigned seed = t.count();
    generator.seed(seed);
/*
    std::ofstream resultat;
    resultat.open("resultatABFLangevin_10traj_beta6.txt");
*/
    int T = 0;
    double Dt = 0.01;
    double beta = 6.0;
    double sigma = 2.0;
    double gamma = beta*sigma*sigma*0.5;
    int N = 120;
    double a = -2;
    double b = 2;
    int L = 1;
    double varmin = 0.01;

    std::vector<double> QX(L);
    std::vector<double> QY(L);
    std::vector<double> PX(L);
    std::vector<double> PY(L);
    for(int l=0;l<L;l++){
        QX[l] = 0;
        QY[l] = 0;
        PX[l] = d(generator);
        PY[l] = d(generator);
    }

    std::vector<double> den(N+2);
    std::vector<double> num(N+2);
    std::vector<double> bias(N+2);
    std::vector<double> carres(N+2);
    for(int i=0;i<N+2;i++){
        num[i] = 0;
        den[i] = 0;
        bias[i] = 0;
        carres[i] = 0;
    }

    while(1){
        /*if(i%100==0){
            std::cout<<"Tour :"<<i<<std::endl;
            std::cout<<"QX = "<<QX[0]<<std::endl;
        }*/
        T++;
        if(T%5000==0){
            //std::cout<<"Tour :"<<T<<std::endl;
            double tst = test(a,b,beta,bias,N);
            //std::cout<<tst<<std::endl;
            if(tst<2.5)
                return T;
        }
        for(int l=0;l<L;l++){
            int ind = toIndex(a,b,QX[l],N);
            den[ind]++;
            num[ind] += DxV(QX[l],QY[l]);
            carres[ind] += DxV(QX[l],QY[l])*DxV(QX[l],QY[l]);
            double m2 = (num[ind]/den[ind])*(num[ind]/den[ind]);
            double var = (carres[ind]/den[ind]-m2)/den[ind];
            if( ((var/m2)<varmin) && (ind<N) )
                bias[ind] = num[ind]/den[ind];
        }
        for(int l=0;l<L;l++){
            int ind = toIndex(a,b,QX[l],N);

            double qx = QX[l];
            double qy = QY[l];
            double px = PX[l];
            double py = PY[l];


            double px14 = (px-0.25*Dt*gamma*px+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);
            double py14 = (py-0.25*Dt*gamma*py+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);

            double px12 = px14-0.5*Dt*(DxV(qx,qy)-bias[ind]);
            double py12 = py14-0.5*Dt*DyV(qx,qy);

            QX[l] = qx+Dt*px12;
            QY[l] = qy+Dt*py12;

            double px34 = px12-0.5*Dt*(DxV(qx,qy)-bias[ind]);
            double py34 = py12-0.5*Dt*DyV(QX[l],QY[l]);

            PX[l] = (px34-0.25*gamma*Dt*px34+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);
            PY[l] = (py34-0.25*gamma*Dt*py34+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);


  /*          resultat<<QX[l];


            if(l!=(L-1))
                resultat<<" ";
            else
                resultat<<std::endl;
    */    }
    }
    std::cout<<"Fin";

    test(a,b,beta,bias,N,T);

//    resultat.close();
}

int ABF_SchemaSto_Rep(){
    std::mt19937_64 generator;
    std::normal_distribution<double> d(0.0,1.0);
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration t = myclock::now() - beginning;
    unsigned seed = t.count();
    generator.seed(seed);
/*
    std::ofstream resultat;
    resultat.open("resultatABFOverdamped_10traj.txt");
*/

    int T = 0;
    int N = 120;
    double a = -2;
    double b = 2;
    int L = 1000;
    double Dt = 0.01;
    double beta = 6.0;
    int nmin = 10;

    std::vector<double> X(L);
    std::vector<double> Y(L);
    for(int l=0;l<L;l++){
        X[l] = 0;
        Y[l] = 0;
    }

    std::vector<double> den(N+2);
    std::vector<double> num(N+2);
    std::vector<double> bias(N+2);
    for(int i=0;i<N+2;i++){
        num[i] = 0;
        den[i] = 0;
        bias[i] = 0;
    }

    while(1){
        /*if(i%1000==0){
            std::cout<<"Tour :"<<i<<std::endl;
            test(a,b,bias,N);
        }*/
        T++;
        if(T%1000==0){
            //std::cout<<"Tour :"<<T<<std::endl;
            double tst = test(a,b,beta,bias,N);
            std::cout<<tst<<std::endl;
            if(tst<2.5)
                return T;
        }
        for(int i=0;i<N;i++){
            num[i] = 0;
            den[i] = 0;
            bias[i] = 0;
        }
        for(int l=0;l<L;l++){
            int ind = toIndex(a,b,X[l],N);
            den[ind]++;
            num[ind] += DxV(X[l],Y[l]);
            if((den[ind]>nmin)&&(ind<N))
                bias[ind] = num[ind]/den[ind];
        }
        for(int l=0;l<L;l++){
            int ind = toIndex(a,b,X[l],N);

            double x = X[l];
            double y = Y[l];

            X[l] = x-Dt*DxV(x,y)+bias[ind]*Dt+sqrt(2*Dt/beta)*d(generator);
            Y[l] = y-Dt*DyV(x,y)+sqrt(2*Dt/beta)*d(generator);

//            resultat<<X[l];


/*
            if(l!=(L-1))
                resultat<<" ";
            else
                resultat<<std::endl;
  */      }
    }
    std::cout<<"Fin";

    test(a,b,beta,bias,N);

  //  resultat.close();
}

int ABF_Schema250_Rep(){
    std::mt19937_64 generator;
    std::normal_distribution<double> d(0.0,1.0);
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration t = myclock::now() - beginning;
    unsigned seed = t.count();
    generator.seed(seed);
/*
    std::ofstream resultat;
    resultat.open("resultatABFLangevin_10traj.txt");
*/
    int T = 0;
    double Dt = 0.01;
    double beta = 6.0;
    double sigma = 2.0;
    double gamma = beta*sigma*sigma*0.5;
    int N = 120;
    double a = -2;
    double b = 2;
    int L = 100000;
    int nmin = 10;

    std::vector<double> QX(L);
    std::vector<double> QY(L);
    std::vector<double> PX(L);
    std::vector<double> PY(L);
    for(int l=0;l<L;l++){
        QX[l] = 0;
        QY[l] = 0;
        PX[l] = d(generator);
        PY[l] = d(generator);
    }

    std::vector<double> den(N+2);
    std::vector<double> num(N+2);
    std::vector<double> bias(N+2);
    for(int i=0;i<N+2;i++){
        num[i] = 0;
        den[i] = 0;
        bias[i] = 0;
    }

    while(1){
        /*if(i%100==0){
            std::cout<<"Tour :"<<i<<std::endl;
            std::cout<<"QX = "<<QX[0]<<std::endl;
            test(a,b,beta,bias,N);
        }*/
        T++;
        if(T%1000==0){
            //std::cout<<"Tour :"<<T<<std::endl;
            double tst = test(a,b,beta,bias,N);
            std::cout<<tst<<std::endl;
            if(tst<2.5)
                return T;
        }
        for(int i=0;i<N;i++){
            num[i] = 0;
            den[i] = 0;
            bias[i] = 0;
        }
        for(int l=0;l<L;l++){
            int ind = toIndex(a,b,QX[l],N);
            den[ind]++;
            num[ind] += DxV(QX[l],QY[l]);
            if((den[ind]>nmin)&&(ind<N))
                bias[ind] = num[ind]/den[ind];
        }
        for(int l=0;l<L;l++){
            int ind = toIndex(a,b,QX[l],N);

            double qx = QX[l];
            double qy = QY[l];
            double px = PX[l];
            double py = PY[l];


            double px14 = (px-0.25*Dt*gamma*px+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);
            double py14 = (py-0.25*Dt*gamma*py+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);

            double px12 = px14-0.5*Dt*(DxV(qx,qy)-bias[ind]);
            double py12 = py14-0.5*Dt*DyV(qx,qy);

            QX[l] = qx+Dt*px12;
            QY[l] = qy+Dt*py12;

            double px34 = px12-0.5*Dt*DxV(QX[l],QY[l]);
            double py34 = py12-0.5*Dt*DyV(QX[l],QY[l]);

            PX[l] = (px34-0.25*gamma*Dt*px34+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);
            PY[l] = (py34-0.25*gamma*Dt*py34+sqrt(0.5*Dt)*sigma*d(generator))/(1.0+0.25*Dt*gamma);


   //         resultat<<QX[l];

/*
            if(l!=(L-1))
                resultat<<" ";
            else
                resultat<<std::endl;
 */       }
    }
    std::cout<<"Fin";

    test(a,b,beta,bias,N);

//    resultat.close();
}

