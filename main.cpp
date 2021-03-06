#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <sstream>
#include <map>

using namespace std;

const int N0 = 1000;
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
    ss << "C:\\users\\Pawel\\Desktop\\BA_degree_";
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

void save_degrees(vector<int> degrees, int t)
{
    //If it should be final degree sequence, sort degrees descending
    if(t==TMAX) sort(degrees.begin(),degrees.end(), greater<int>());

    stringstream ss;

    //Prepare filename
    ss << "C:\\users\\Pawel\\Desktop\\RG_degree_";
    if(t<TMAX) ss  << t;
    else ss << "sequence";
    ss << ".txt";

    string filename = ss.str();

    ofstream out(filename);

    for(int i=0;i<degrees.size();i++)
    {
        out << degrees[i] << endl;
    }

    out.close();
}

void create_degree_from_adjacency_list(map<int,vector<int>> adjacency_list, int t)
{
    vector<int> degrees;

    //Creates degree sequence by counting each node connections
    typedef map<int, vector<int>>::iterator it_type;
    for(it_type iterator = adjacency_list.begin(); iterator != adjacency_list.end(); iterator++)
    {
       degrees.push_back(iterator->second.size());
    }

    //If it should be final degree sequence, sort degrees descending
    if(t==TMAX) sort(degrees.begin(),degrees.end(), greater<int>());

    stringstream ss;

    //Prepare filename
    ss << "C:\\users\\Pawel\\Desktop\\NG_degree_";
    if(t<TMAX) ss  << t;
    else ss << "sequence";
    ss << ".txt";

    string filename = ss.str();

    ofstream out(filename);

    for(int i=0;i<degrees.size();i++)
    {
        out << degrees[i] << endl;
    }

    out.close();
}

int main()
{
    default_random_engine generator;

    //##################################
    //##   PART 1 - BARABASI-ALBERT   ##
    //##################################

    vector<int> stubs;

    //Adds n0 initial nodes with every node connected to next one
    for(int i=1; i<N0; i++)
    {
        stubs.push_back(i);
        stubs.push_back(i+1);
    }

    //Closes the circle
    stubs.push_back(N0);
    stubs.push_back(1);

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
                //If it is the first stub chosen break
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



    //##################################
    //##    PART 2 - RANDOM GROWTH    ##
    //##################################

    //Creates initial graph in form of degrees set
    vector<int> degrees;
    for(int i=0; i<N0; i++)
    {
        degrees.push_back(2);
    }

    for(int t=1; t<=TMAX; t++)
    {
        chosen.clear();

        //Creates random number distribution
        uniform_int_distribution<int> distribution(0, degrees.size()-1);

        //Looks for M0 different vertices to connect
        for(int m=0; m<M0; m++)
        {
            //Looks for random vertex which was not chosen before
            while(true)
            {
                index = distribution(generator);
                //If it is the first vertex chosen break
                if(chosen.empty()) break;
                //If it was not chosen break
                if(find(chosen.begin(), chosen.end(), index) == chosen.end()) break;
                //Else try again
            }

            //Increase degree of found vertex
            degrees[index]++;
            //Save found vertex
            chosen.push_back(index);
        }

        //Add new node with initial degree M0
        degrees.push_back(M0);

        if(t==1||t==10||t==100||t==1000)
        {
            save_degrees(degrees, t);
        }
    }

    save_degrees(degrees,TMAX);



    //##################################
    //##      PART 3 - NO GROWTH      ##
    //##################################

    //Creates initial graph in form of adjacency list
    int before;
    int after;
    map<int, vector<int>> adjacency_list;
    for(int i=0;i<N0;i++)
    {
        if(i==0) before = N0-1;
        else before = i-1;
        if(i==N0-1) after = 0;
        else after = i+1;

        vector<int> connections;
        connections.push_back(before);
        connections.push_back(after);

        adjacency_list.insert(pair<int, vector<int>>(i,connections));
    }

    int vertex;
    for(int t=1; t<=TMAX; t++)
    {
        chosen.clear();

        //Creates random number distribution
        uniform_int_distribution<int> distribution(0, adjacency_list.size()-1);
        //Choose randomly one of vertices
        vertex = distribution(generator);

        vector<int> conn = adjacency_list.find(vertex)->second;

        //Looks for M0 different vertices to connect
        for(int m=0; m<M0; m++)
        {

            while (true)
            {
                //Look for another vertex
                index = distribution(generator);
                //If they were not connected break
                if (find(conn.begin(), conn.end(), index) == conn.end()) break;
                //Else try again
            }

            //Add connections
            adjacency_list.find(vertex)->second.push_back(index);
            adjacency_list.find(index)->second.push_back(vertex);
        }

        if(t==1||t==10||t==100||t==1000)
        {
            create_degree_from_adjacency_list(adjacency_list,t);
        }
    }

    create_degree_from_adjacency_list(adjacency_list,TMAX);

    return 0;
}
