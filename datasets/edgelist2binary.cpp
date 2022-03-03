//
//  main.cpp
//  gen_bin_bigraph
//
//  Created by kai on 14/9/2021.
//

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <stack>
#include <set>
#include <map>
//#define NDEBUG // must precede cassert to disable assert.
#include <cassert>
using namespace std;
using ui = unsigned int;

static FILE *open_file(const char *file_name, const char *mode) {
    FILE *f = fopen(file_name, mode);
    if(f == nullptr) {
        printf("Can not open file: %s\n", file_name);
        exit(1);
    }
    return f;
}

static std::string integer_to_string(long long number) {
    std::vector<ui> sequence;
    if(number == 0) sequence.push_back(0);
    while(number > 0) {
        sequence.push_back(number%1000);
        number /= 1000;
    }

    char buf[5];
    std::string res;
    for(unsigned int i = sequence.size();i > 0;i --) {
        if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
        else sprintf(buf, ",%03u", sequence[i-1]);
        res += std::string(buf);
    }
    return res;
}

int main(int argc, const char * argv[]) {
    if(argc < 2) {
        printf("Usage: [1].exe [2].original_graph \n");
        return 0;
    }
    cout<<"Graph name = "<<argv[1]<<endl;
    
    vector<ui> nodes;
    vector<pair<ui,ui> > edges;
    string original_graph = string(argv[1]);
    
    FILE *f = open_file(original_graph.c_str(), "r");
    
    char buf[1024];
    ui t_n1, t_n2, t_m;
    if(fgets(buf, 1024, f) != NULL) {
        for(ui j = 0;buf[j] != '\0';j ++) if(buf[j] < '0'||buf[j] > '9') buf[j] = ' ';
        sscanf(buf, "%u%u%u", &t_n1, &t_n2, &t_m);
    }
    
    ui a, b;
    
    while(fgets(buf, 1024, f)) {
        char comment = 1;
        for(ui j = 0;buf[j] != '\0';j ++) if(buf[j] != ' '&&buf[j] != '\t') {
            if(buf[j] >= '0'&&buf[j] <= '9') comment = 0;
            break;
        }
        if(comment) continue;

        for(ui j = 0;buf[j] != '\0';j ++) if(buf[j] < '0'||buf[j] > '9') buf[j] = ' ';
        sscanf(buf, "%u%u", &a, &b);
        if(a == b) continue;
        nodes.push_back(a);
        nodes.push_back(b);
        edges.push_back(make_pair(a,b));
        edges.push_back(make_pair(b,a));
    }

    fclose(f);
    
    sort(nodes.begin(), nodes.end());
    nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
    
    if(nodes.size() != (t_n1+t_n2)) {cout<<"Nodes.size() != (t_n1+t_n2)"<<endl; exit(1);}
    
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());

    map<ui,ui> M;
    for(ui i = 0;i < nodes.size();i ++) M[nodes[i]] = i;

    char preserved = 1;
    for(ui i = 0;i < nodes.size();i ++) if(nodes[i] != i) preserved = 0;
    if(!preserved) printf("Node ids are not preserved!\n");

    ui n = (ui)nodes.size();
    ui m = (ui)edges.size();
    printf("n = %s, m = %s\n", integer_to_string(n).c_str(), integer_to_string(m/2).c_str());
    
    if(m/2 != t_m) {cout<<"m/2 != t_m"<<endl; exit(1);}
    
    original_graph.erase(original_graph.end() - 4, original_graph.end());
    f = open_file((original_graph + ".bin").c_str(), "wb");
    ui tt = sizeof(ui);
    fwrite(&tt, sizeof(ui), 1, f);
    fwrite(&t_n1, sizeof(ui), 1, f);
    fwrite(&t_n2, sizeof(ui), 1, f);
    fwrite(&m, sizeof(ui), 1, f);
    ui j = 0;
    for(ui i = 0;i < n;i ++) {
        while(j < m&&edges[j].first == nodes[i]) {
            fwrite(&i, sizeof(ui), 1, f);
            fwrite(&M[edges[j].second], sizeof(ui), 1, f);
            ++ j;
        }
    }
    fclose(f);
    
    cout<<"Finish generating binary version!"<<endl;
    return 0;
}
