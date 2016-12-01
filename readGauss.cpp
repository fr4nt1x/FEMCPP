#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>


using namespace std;

int main()
{
    int n = 3;
    ifstream leggauss;
    leggauss.open("leggauss",ios::in);
    string line;
    string deg = "Degree: ";
    
    cout<<setprecision(64);
    if (leggauss.is_open())
    {
        bool startread= false;
        while(getline(leggauss,line))
        {

            if(line ==deg +to_string(n+1))
            {
                startread = false;
            }

            if (startread)
            {
                istringstream ss(line);
                while(!ss.eof())
                {
                    string x;
                    getline(ss,x,' ');
                    double test = stod(x);
                    cout<<test<<endl;
                }
                cout<<endl;
            }

            if (line ==deg + to_string(n))
            {
                startread = true;

            }
        }

    }

    leggauss.close();
    return 0;
}

