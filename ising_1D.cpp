// 2021/11/10 ising one-dimensional
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <time.h>

using namespace std;

const int ne = 1000; /* number of objects */
double kt; // = 10.0; /* temperature */
const double J = -1.0; /* exchange energy */


//---- maximam of random value
#ifndef RAND_MAX
#define RAND_MAX     (0x7fffffff)
#endif

//---- sets seed of random from time
void Srand(void)
{
    srand(time(NULL));
}

//---- uniformal random within [0,1]
double Urand(void)
{
    return (double)rand() / RAND_MAX;
}

int dice() {
    //0から(ne-2)までの乱数
    return rand() % (ne - 1); 
}


/* function returns energy of the system */
double energy(vector<int> array, double J)
{
    double sum = 0.0;
    int max = array.size();
    for (int i = 0; i < (max - 1); i++)		/* loop through elements */
    {
        sum += array[i] * array[i + 1];
    }
    return (J * sum);
}


void write_data(vector<vector<double>> &seqX) {
    string fname = "D:/VCPP/VC++2019/data01.txt";
    int dataSizeMax = 100000;
    int dataSize = seqX.size();
    int interval = 1;
    if (dataSize > dataSizeMax) interval = dataSize / dataSizeMax; 	
    ofstream outputfile(fname);
    outputfile.precision(6); // 小数点以下の出力桁を6桁に指定 	
    for (int i = 0; i < seqX.size(); i++) {
        if (i % interval == 0) {
            outputfile << seqX[i][1] << " " << seqX[i][2] << " " << seqX[i][3] << endl; //データは空白区切り 		} 	
        }
        outputfile.close();
    }
}

void gnuplot() {
    FILE* gp = _popen("C:/PROGRA~1/gnuplot/bin/gnuplot.exe -persist", "w");
    if (gp == NULL) return;
    fputs("load 'D:/VCPP/VC++2019/computationPhysics/ising_1D/ising_1D/gChart_2D.txt'\n", gp);
    fflush(gp);
    _pclose(gp);
}


void calc_physical()
{
    //2012 Jun 16, 熱力学量算出
    //2012 Jun 11, 1次元ising, Metropolis algorithm


    int element;
    vector<int> array(ne);
    double olden, newen;
    int Nmps = 50;//熱力学量算出回数

    vector<double> en(Nmps);
    double C;

    newen = 0.0;//dummy
    int ms = ne * 4;//1 Metropolis step数

    double dkt = 0.2;
    int nkt = (int)(10.0 / dkt);
    vector<double> u(nkt);
    vector<double> sh(nkt);

    //string fname = "D:/VCPP/VC++2019/computationPhysics/ising_1D/ising_1D/gChart_2D.txt";
    string fname1 = "D:/VCPP/VC++2019/computationPhysics/ising_1D/ising_1D/energy.dat";
    string fname2 = "D:/VCPP/VC++2019/computationPhysics/ising_1D/ising_1D/energy_exact.dat";
    ofstream outputfile(fname1);
    outputfile.precision(6); // 小数点以下の出力桁を6桁に指定 	
    ofstream outputfile2(fname2);
    outputfile2.precision(6);
    //for (int i = 0; i < seqX.size(); i++) {
    //    outputfile << seqX[i][1] << " " << seqX[i][2] << " " << seqX[i][3] << endl; //データは空白区切り 		
    //}
    //outputfile.close();
    

    //random_device rd;
    //mt19937 mt(rd()); //seedにrandom_deviceで生成した乱数
    //uniform_int_distribution<int> dice(0, ne - 2);
    //uniform_real_distribution<double> real_rnd(0.0, 1.0);
    srand(time(NULL));

    for (int i = 0; i < ne; i++) array[i] = 1;          /* uniform start */

    for (int j = 0; j < nkt; j++)
    {
        kt = 0.1 + dkt * (double)j;
        //kt<=1で熱平衡に達するのが緩やかにつき
        if (kt <= 0.5) ms = ne * 10;
        else  if (kt <= 4.0) ms = ne * 2; //ne * 5; //ne * 10;
        else ms = ne * 4;  ne * 2;
        for (int k = 0; k < Nmps; k++)
        {

            //mt19937 mt(rd()); //seedにrandom_deviceで生成した乱数
            //uniform_int_distribution<int> dice(0, ne - 2);
            //uniform_real_distribution<double> real_rnd(0.0, 1.0);
            //1 Metropolis step
            //初期状態設定(周期境界条件)
            //uniform startの方が理論値に近い、ランダムstartではkt<1で大きく乖離
            //for (int i = 0; i < ne; i++) array[i] = 1;          /* uniform start */

            for (int i = 0; i <= ms; i++)
            {
                olden = energy(array, J); /* initial energy  */
                //element = rnd.Next(ne - 1); //0から(ne-2)までの乱数
                //element = dice(mt);
                element = dice();
                array[element] *= -1; /* change spin */
                if (element == 0) array[ne - 1] = array[0];

                newen = energy(array, J); /* calculate new energy */
                //if ((newen > olden) && (exp((-newen + olden) / kt) <= rnd.NextDouble())) //0から1.0までの乱数
                //if ((newen > olden) && (exp((-newen + olden) / kt) <= real_rnd(mt)))
                if ((newen > olden) && (exp((-newen + olden) / kt) <= Urand()))
                {
                    array[element] = array[element] * (-1);	    /* reject change */
                    if (element == 0) array[ne - 1] = array[0];
                }
            }
            newen = energy(array, J);
            en[k] = newen;
        }

        double s1 = 0.0, s2 = 0.0;
        for (int i = 0; i < Nmps; i++)
        {
            s1 += en[i];
            s2 += en[i] * en[i];
        }
        s1 /= (Nmps);
        s2 /= (Nmps);
        double U = s1 / ne;
        u[j] = s1 / ne;
        C = (s2 - s1 * s1) / (kt * kt * ne);
        sh[j] = (s2 - s1 * s1) / (kt * kt * ne);
        cout << "kT= " << kt << endl;
        cout << "エネルギー= " << U << ", " << "比熱= " << C << endl;
        outputfile << kt << " " << U << endl;
        outputfile2 << kt << " " << -J * tanh(J / kt) << endl;

    }
    outputfile.close();
    outputfile2.close();
    gnuplot();



    //chartSetup();
    //addChart2();
    //text("1格子当りの比熱");
    ////xyLimit(0.0, 10.0, 0.0, 0.001);//x軸y軸上下限設定
    //for (i = 0; i < nkt; i++)
    //{
    //    kt = 0.1 + dkt * (double)i;
    //    setLine(kt, sh[i]);
    //    double s = J * J / (kt * kt * Math.Cosh(J / kt) * Math.Cosh(J / kt));
    //    setLine2(kt, s);
    //}
    //showChart();

}



int main()
{
    calc_physical();

}

