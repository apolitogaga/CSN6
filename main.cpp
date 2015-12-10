#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <sstream>

using namespace std;

const int N0 = 100;
const int M0 = 10;
const int TMAX = 10000;

void create_degree_sequence_from_stubs(vector<int> stubs, int t)
{
    const int N = N0 + t;

    vector<int> degree;
    degree.resize(N);

    //Creates degree sequence by counting each node appearance in stubs
    for(int i=0; i<stubs.size(); i++)
    {
        degree[stubs[i]-1]++;
    }

    //If it should be final degree sequence, sort degrees descending
    if(t==TMAX) sort(degree.begin(),degree.end(), greater<int>());

    stringstream ss;

    //Prepare filename
    ss << "C:\\users\\Pawel\\Desktop\\degree_";
    if(t<TMAX) ss  << t;
    else ss << "sequence";
    ss << ".txt";

    string filename = ss.str();

    ofstream out(filename);

    for(int i=0;i<degree.size();i++)
    {
        out << degree[i] << endl;
    }

    out.close();
}

int main()
{
    default_random_engine generator;

    vector<int> stubs;

    //Adds n0 initial nodes with every node connected to next one
    for(int i=1; i<N0; i++)
    {
        stubs.push_back(i);
        stubs.push_back(i+1);
    }

    int node = N0;
    int index = 0;
    vector<int> chosen;

    //Adds new node and connects it randomly with existing nodes TMAX times
    for(int t=1; t<=TMAX; t++)
    {
        node++;
        chosen.clear();

        //Creates random number distribution
        uniform_int_distribution<int> distribution(0, stubs.size()-1);

        for(int m=0; m<M0; m++)
        {
            //Looks for random stub which was not chosen before
            while(true)
            {
                index = distribution(generator);
                //If it is the first stub break
                if(chosen.empty()) break;
                //If it was not chosen break
                if(find(chosen.begin(), chosen.end(), stubs[index]) == chosen.end()) break;
                //Else try again
            }

            //Creates new edge with new node
            stubs.push_back(node);
            stubs.push_back(stubs[index]);

            //Saves chosen stub to not use it again
            chosen.push_back(stubs[index]);
        }

        if(t==1||t==10||t==100||t==1000)
        {
            create_degree_sequence_from_stubs(stubs,t);
        }
    }

    create_degree_sequence_from_stubs(stubs,TMAX);

    return 0;
}
