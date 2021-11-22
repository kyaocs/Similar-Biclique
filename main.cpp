//
//  main.cpp
//  MSBE_GitHub
//
//  Created by kai on 20/4/2021.
//

#include "Timer.h"
#include "Utility.h"
#include "LinearHeap.h"

int maxTime = 3600;  //seconds
bool over_time_flag = false;
double startTime;
double durTime;
int binary_version;

ui n, n1, n2, m;
ui * pstart;
ui * edges;
int * degree; //vertex degree
int * TMPdeg; //used in vertex reduction
int * Sdeg; //used in vertex reduction
int * os;

//the remaining graph
ui r_n, r_m;
ui * r_pstart;
ui * r_edges;
int * r_degree;
ui * peel_s;
ui * oid;
ui * nid;

ui max_CS_deg;
int ** Matrix;
ui * trans;
int * inPQ;
int * deg_inP;
int * pdeg_inP;
int * ndeg_inP;
double seg_times;
double rg_limit;
bool * del_ver;
int vr_method;
double epsilon;
int tau;
long long time_on_vr = 0;
long long time_on_enum = 0;


vector<pair<vector<ui>, vector<ui>>> results;
long long MSBC_num = 0;
int max_MSBC_size = 0;
int min_MSBC_size = INF;
int thre_seg = 12;
double triscr = 0.05;
long long cnt_prune_by_ub = 0;
long long cnt_cmp_by_js = 0;
long long total_js = 0;
double min_density;
double max_density;
double ave_density;
long long ave_density_cnt;

bool cmp_degP (ui u, ui v)
{
    return deg_inP[u] < deg_inP[v];
}

class Itval
{
public:
        ui s_idx; //true vertex id
        ui e_idx; //true vertex id
        double min_score;
        double max_score;
        int c;
    Itval() {
        s_idx = 0;
        e_idx = 0;
        min_score = 0;
        max_score = 0;
        c = 0;
    }
    Itval(ui _s, ui _e, double _mins, double _maxs, int _c){
        s_idx = _s;
        e_idx = _e;
        min_score = _mins;
        max_score = _maxs;
        c = _c;
    }
};

class Range
{
public:
    ui coreV;  //vertex id
    int rgC;  //L vid -> R vid
    ui Lidx;  //L position
    ui Ridx;  //R position
    Range() {
        coreV = 0;
        rgC = 0;
        Lidx = 0;
        Ridx = 0;
    }
    Range(ui _coreV, int _rgC, ui _Lidx, ui _Ridx){
        coreV = _coreV;
        rgC = _rgC;
        Lidx = _Lidx;
        Ridx = _Ridx;
    }
};

class node
{
public:
    node * L;
    node * R;
    node * P;
    double mins;
    double maxs;
    int mark;
    int isv;
    int idx;
    node() {
        L = nullptr;
        R = nullptr;
        P = nullptr;
        mins = 0.0;
        maxs = 0.0;
        mark = 0;
        isv = 0;
        idx = 0;
    }
    node(node * _L, node * _R, node * _P, double _mins, double _maxs, int _mark, int _isv, int _idx) {
        L = _L;
        R = _R;
        P = _P;
        mins = _mins;
        maxs = _maxs;
        mark = _mark;
        isv = _isv;
        idx = _idx;
    }
};

vector<vector<Itval>> vsn;
vector<vector<pair<ui, double>>> index_storeall;

void load_graph_binary(string graph_name)
{
    FILE *f = Utility::open_file(graph_name.c_str(), "rb");

    ui tt;
    fread(&tt, sizeof(ui), 1, f);
    if(tt != sizeof(ui)) {
        printf("sizeof unsigned int is different: .bin(%u), machine(%lu)\n", tt, sizeof(ui));
        return ;
    }
    fread(&n1, sizeof(ui), 1, f);
    fread(&n2, sizeof(ui), 1, f);
    n = n1 + n2;
    fread(&m, sizeof(ui), 1, f);
    cout<<"Graph name = "<<graph_name<<" ( n1 = "<<n1<<", n2 = "<<n2<<", m = "<<m<<" )."<<endl;
    
    map<ui, set<ui>> biG;
    ui read_edges = 0;
    ui tu, tv;
    while (read_edges ++ < m) {
        fread(&tu, sizeof(ui), 1, f);
        fread(&tv, sizeof(ui), 1, f);
        biG[tu].insert(tv);
        biG[tv].insert(tu);
    }
    fclose(f);
    
    pstart = new ui[n+1];
    edges = new ui[m];
    degree = new int[n];
    TMPdeg = new int[n];
    Sdeg = new int[n];
    os = new int[n];
    del_ver = new bool[n];
    memset(del_ver, 0, sizeof(bool)*n);
        
    pstart[0] = 0;
    for(ui i = 0; i < n; i++){
        const set<ui> & nei = biG[i];
        ui s_idx = pstart[i];
        for(auto e : nei) edges[s_idx++] = e;
        pstart[i+1] = s_idx;
        degree[i] = nei.size();
        TMPdeg[i] = nei.size();
        Sdeg[i] = 0;
    }
    for(ui i = 0; i < n1; i++) os[i] = 1;
    for(ui i = n1; i < n; i++) os[i] = 2;
}

void load_graph_txt(string graph_name)
{
    ifstream input_file(graph_name, ios::in);
    map<ui, set<ui>> biG;
    if (!input_file.is_open()){
        cout << "cannot open file : "<<graph_name<<endl;exit(1);
    }
    else{
        input_file >> n1 >> n2 >> m;
        n = n1 + n2;
        cout<<graph_name<<" : n1 = "<<n1<<", n2 = "<<n2<<", m = "<<m<<endl;
        ui tu, tv;
        while (input_file >> tu >> tv) {
            assert(tu != tv);
            assert(tu >= 0 && tu < n);
            assert(tv >= 0 && tv < n);
            biG[tu].insert(tv);
            biG[tv].insert(tu);
        }
        assert(biG.size() == n);
        m = 0;
        for(auto e : biG) m += e.second.size();
        assert(m%2 == 0); m /= 2;
        cout<<"after loading, n = "<<n<<", m = "<<m<<", ";
        input_file.close();
    }
    
    pstart = new ui[n+1];
    edges = new ui[2*m];
    degree = new int[n];
    TMPdeg = new int[n];
    Sdeg = new int[n];
        
    pstart[0] = 0;
    for(ui i = 0; i < n; i++){
        const set<ui> & nei = biG[i];
        ui s_idx = pstart[i];
        for(auto e : nei) edges[s_idx++] = e;
        pstart[i+1] = s_idx;
        degree[i] = nei.size();
        TMPdeg[i] = nei.size();
        Sdeg[i] = 0;
    }
    assert(pstart[n] == 2*m);
}

void load_index_LG(string graph_name)
{
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_" + to_string((int)(seg_times*100)));
    graph_name.append("_" + to_string(0));
    graph_name.append("_LG.bin");
    cout<<"Start loading index : "<<graph_name<<endl;
    
    FILE * f = Utility::open_file(graph_name.c_str(), "rb");
    vsn.resize(n);
    
    int vid, numofI;
    ui s_vid, e_vid;
    char minS, maxS;
    int cnt;
    unsigned long result;
    
    while(fread(&vid, sizeof(int), 1, f) == 1) {
        fread(&numofI, sizeof(int), 1, f);
        vector<Itval> & tmp_vec = vsn[vid];
        if(numofI > 0){
            while (1) {
                fread(&s_vid, sizeof(ui), 1, f);
                fread(&e_vid, sizeof(ui), 1, f);
                fread(&minS, sizeof(char), 1, f);
                fread(&maxS, sizeof(char), 1, f);
                fread(&cnt, sizeof(int), 1, f);
                tmp_vec.push_back(Itval(s_vid, e_vid, (double)minS/100, (double)maxS/100, cnt));
                -- numofI;
                if(numofI == 0) break;
            }
        }
        else { //indi
            numofI = -numofI;
            while (1) {
                fread(&s_vid, sizeof(ui), 1, f);
                fread(&minS, sizeof(char), 1, f);
                tmp_vec.push_back(Itval(s_vid, s_vid, (double)minS/100, (double)minS/100, 1));
                -- numofI;
                if(numofI == 0) break;
            }
        }
    }
    fclose(f);
}

void load_index_GR(string graph_name)
{
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_GR.txt");
    ifstream fin;
    fin.open(graph_name);
    if(!fin.is_open()){
        cout<<"cannot open index GR."<<endl; exit(1);
    }
    vsn.resize(n);
    string line;
    ui vid, s_vid, e_vid, cnt;
    double minS, maxS;
    while (getline(fin, line)) {
        stringstream ss(line);
        ss >> vid;
        assert(vid >=0 && vid < n);
        vector<Itval> & tmp_vec = vsn[vid];
        while (ss >> s_vid >> e_vid >> minS >> maxS >> cnt) {
            tmp_vec.push_back(Itval(s_vid, e_vid, minS, maxS, cnt));
        }
    }
    fin.close();
    cout<<"*** finish load_index_GR ***"<<endl;
}

void load_index_GRL(string graph_name)
{
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_GRL.txt");
    ifstream fin;
    fin.open(graph_name);
    if(!fin.is_open()){
        cout<<"cannot open index GRL."<<endl; exit(1);
    }
    
    vsn.resize(n);
    
    string line;
    ui vid, s_vid, e_vid, cnt;
    double minS, maxS;
    while (getline(fin, line)) {
        stringstream ss(line);
        ss >> vid;
        assert(vid >=0 && vid < n);
        vector<Itval> & tmp_vec = vsn[vid];
        while (ss >> s_vid >> e_vid >> minS >> maxS >> cnt) {
            tmp_vec.push_back(Itval(s_vid, e_vid, minS, maxS, cnt));
        }
    }
    fin.close();
    cout<<"*** finish load_index_GRL ***"<<endl;
}

void load_index_SS(string graph_name)
{
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_" + to_string((int)(seg_times*100)));
    graph_name.append("_" + to_string((int)(rg_limit*100)));
    graph_name.append("_SS.bin");

    cout<<"Start loading index : "<<graph_name<<endl;
    
    FILE * f = Utility::open_file(graph_name.c_str(), "rb");
    vsn.resize(n);
    
    int vid, numofI;
    ui s_vid, e_vid;
    char minS, maxS;
    int cnt;
    unsigned long result;
    
    while(fread(&vid, sizeof(int), 1, f) == 1) {
        fread(&numofI, sizeof(int), 1, f);
        vector<Itval> & tmp_vec = vsn[vid];
        if(numofI > 0){
            while (1) {
                fread(&s_vid, sizeof(ui), 1, f);
                fread(&e_vid, sizeof(ui), 1, f);
                fread(&minS, sizeof(char), 1, f);
                fread(&maxS, sizeof(char), 1, f);
                fread(&cnt, sizeof(int), 1, f);
                tmp_vec.push_back(Itval(s_vid, e_vid, (double)minS/100, (double)maxS/100, cnt));
                -- numofI;
                if(numofI == 0) break;
            }
        }
        else { //indi
            numofI = -numofI;
            while (1) {
                fread(&s_vid, sizeof(ui), 1, f);
                fread(&minS, sizeof(char), 1, f);
                tmp_vec.push_back(Itval(s_vid, s_vid, (double)minS/100, (double)minS/100, 1));
                -- numofI;
                if(numofI == 0) break;
            }
        }
    }
    fclose(f);
}

void load_index_SA(string graph_name)
{
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_SA.txt");
    ifstream fin;
    fin.open(graph_name);
    if(!fin.is_open()){
        cout<<"cannot open index SA."<<endl; exit(1);
    }
    
    index_storeall.resize(n);
    
    string line;
    ui u, v;
    double simscore;
    while (getline(fin, line)) {
        stringstream ss(line);
        ss >> u;
        assert(u >=0 && u < n);
        vector<pair<ui, double>> & tmp_vec = index_storeall[u];
        while (ss >> v >> simscore ) {
            assert(v >= 0 && v < n);
            assert(simscore >= 0 && simscore <= 1);
            tmp_vec.push_back(make_pair(v, simscore));
        }
    }
    fin.close();
    cout<<"*** finish load_index_SA ***"<<endl;
//    cout<<"index_storeall information : "<<endl;
//    for(ui i = 0; i < n; i++){
//        cout<<i<<" : ";
//        for(auto e : index_storeall[i]){
//            cout<<e.first<<","<<e.second<<" ";
//        }
//        cout<<endl;
//    }
}

void graph_property()
{
    int u_max_span = 0, u_min_span = INF, v_max_span = 0, v_min_span = INF;
    int u_ave_span = 0, v_ave_span = 0;
        
    for(ui i = 0; i < n1; i++){
        int first_id = -1;
        int last_id;
        for(ui j = pstart[i]; j < pstart[i+1]; j++){
            ui v = edges[j];
            if(first_id == -1) first_id = v;
            last_id = v;
        }
        assert(last_id >= first_id);
        int tmp_span = last_id - first_id + 1;
        u_ave_span += tmp_span;
        if(tmp_span > u_max_span) u_max_span = tmp_span;
        if(tmp_span < u_min_span) u_min_span = tmp_span;
        assert(tmp_span >= degree[i]);
    }
    for(ui i = n1; i < n; i++){
        int first_id = -1;
        int last_id;
        for(ui j = pstart[i]; j < pstart[i+1]; j++){
            ui v = edges[j];
            if(first_id == -1) first_id = v;
            last_id = v;
        }
        assert(last_id >= first_id);
        int tmp_span = last_id - first_id + 1;
        v_ave_span += tmp_span;
        if(tmp_span > v_max_span) v_max_span = tmp_span;
        if(tmp_span < v_min_span) v_min_span = tmp_span;
        assert(tmp_span >= degree[i]);
    }
    
    cout<<"U: ave spa = "<<u_ave_span/n1<<" ["<<u_min_span<<", "<<u_max_span<<"]"<<endl;
    cout<<"V: ave spa = "<<v_ave_span/n2<<" ["<<v_min_span<<", "<<v_max_span<<"]"<<endl;
        
    long long eps00 = 0, eps01 = 0, eps12 = 0, eps23 = 0, eps34 = 0, eps45 = 0, eps56 = 0, eps67 = 0, eps78 = 0, eps89 = 0, eps910 = 0, eps111 = 0;
    long long Ueps00 = 0, Ueps01 = 0, Ueps12 = 0, Ueps23 = 0, Ueps34 = 0, Ueps45 = 0, Ueps56 = 0, Ueps67 = 0, Ueps78 = 0, Ueps89 = 0, Ueps910 = 0, Ueps111 = 0;
    long long Veps00 = 0, Veps01 = 0, Veps12 = 0, Veps23 = 0, Veps34 = 0, Veps45 = 0, Veps56 = 0, Veps67 = 0, Veps78 = 0, Veps89 = 0, Veps910 = 0, Veps111 = 0;
    
    long long U_2hop_nei_cnt = 0, U_2hop_nei_span = 0;
    long long U_min_2hop_deg = INF, U_max_2hop_deg = 0, U_min_2hop_span = INF, U_max_2hop_span = 0;
    
    ui * c = new ui[n];
    memset(c, 0, sizeof(ui)*n);
    
    for(ui u = 0; u < n1; u++){
        
        int spanL = INF, spanR = 0;
        
        vector<ui> two_hop_nei;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(c[w] == 0) {
                    two_hop_nei.push_back(w);
                    if(w > spanR) spanR = w;
                    if(w < spanL) spanL = w;
                }
                ++ c[w];
            }
        }
        if(two_hop_nei.empty()) continue;

        U_2hop_nei_cnt += two_hop_nei.size();
        assert(spanR >= spanL);
        int range_span = spanR - spanL + 1;
        U_2hop_nei_span += range_span;
        
        if(two_hop_nei.size() > U_max_2hop_deg) U_max_2hop_deg = two_hop_nei.size();
        if(two_hop_nei.size() < U_min_2hop_deg) U_min_2hop_deg = two_hop_nei.size();
        
        if(range_span > U_max_2hop_span) U_max_2hop_span = range_span;
        if(range_span < U_min_2hop_span) U_min_2hop_span = range_span;
        
        for(auto e : two_hop_nei){
            assert(degree[u] >= c[e]);
            assert(degree[e] >= c[e]);
            double simscore = (double) c[e] / (degree[u] + degree[e] - c[e]);
            
            c[e] = 0;
            
            if(simscore == 0)
                {eps00 ++; Ueps00 ++;}
            else if(simscore > 0 && simscore < 0.1)
                {eps01 ++; Ueps01 ++;}
            else if (simscore >= 0.1 && simscore < 0.2)
                {eps12 ++; Ueps12 ++;}
            else if (simscore >= 0.2 && simscore < 0.3)
                {eps23 ++; Ueps23 ++;}
            else if (simscore >= 0.3 && simscore < 0.4)
                {eps34 ++; Ueps34 ++;}
            else if (simscore >= 0.4 && simscore < 0.5)
                {eps45 ++; Ueps45 ++;}
            else if (simscore >= 0.5 && simscore < 0.6)
                {eps56 ++; Ueps56 ++;}
            else if (simscore >= 0.6 && simscore < 0.7)
                {eps67 ++; Ueps67 ++;}
            else if (simscore >= 0.7 && simscore < 0.8)
                {eps78 ++; Ueps78 ++;}
            else if (simscore >= 0.8 && simscore < 0.9)
                {eps89 ++; Ueps89 ++;}
            else if (simscore >= 0.9 && simscore < 1)
                {eps910 ++; Ueps910 ++;}
            else{
                assert(simscore == 1);
                eps111 ++; Ueps111 ++;
            }
        }
    }
        
    cout<<endl<<"U: 2hop deg = "<<U_2hop_nei_cnt/n1<<" ["<<U_min_2hop_deg<<", "<<U_max_2hop_deg<<"]"<<endl;
    cout<<"U: 2hop spa = "<<U_2hop_nei_span/n1<<" ["<<U_min_2hop_span<<", "<<U_max_2hop_span<<"]"<<endl;
    
    long long V_2hop_nei_cnt = 0, V_2hop_nei_span = 0;
    long long V_min_2hop_deg = INF, V_max_2hop_deg = 0, V_min_2hop_span = INF, V_max_2hop_span = 0;
    
    for(ui u = n1; u < n; u++){
        if(u%100000 == 0) cout<<"v"<<u<<" ";
        
        int spanL = INF, spanR = 0;
        
        vector<ui> two_hop_nei;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(c[w] == 0) {
                    two_hop_nei.push_back(w);
                    if(w > spanR) spanR = w;
                    if(w < spanL) spanL = w;
                }
                ++ c[w];
            }
        }
        if(two_hop_nei.empty()) continue;
        
        V_2hop_nei_cnt += two_hop_nei.size();
        assert(spanR >= spanL);
        int range_span = spanR - spanL + 1;
        V_2hop_nei_span += range_span;
                
        if(two_hop_nei.size() > V_max_2hop_deg) V_max_2hop_deg = two_hop_nei.size();
        if(two_hop_nei.size() < V_min_2hop_deg) V_min_2hop_deg = two_hop_nei.size();
        
        if(range_span > V_max_2hop_span) V_max_2hop_span = range_span;
        if(range_span < V_min_2hop_span) V_min_2hop_span = range_span;
        
        for(auto e : two_hop_nei){
            assert(degree[u] >= c[e]);
            assert(degree[e] >= c[e]);
            double simscore = (double) c[e] / (degree[u] + degree[e] - c[e]);
            
            c[e] = 0;
            
            if(simscore == 0)
                {eps00 ++; Veps00 ++;}
            else if(simscore > 0 && simscore < 0.1)
                {eps01 ++; Veps01 ++;}
            else if (simscore >= 0.1 && simscore < 0.2)
                {eps12 ++; Veps12 ++;}
            else if (simscore >= 0.2 && simscore < 0.3)
                {eps23 ++; Veps23 ++;}
            else if (simscore >= 0.3 && simscore < 0.4)
                {eps34 ++; Veps34 ++;}
            else if (simscore >= 0.4 && simscore < 0.5)
                {eps45 ++; Veps45 ++;}
            else if (simscore >= 0.5 && simscore < 0.6)
                {eps56 ++; Veps56 ++;}
            else if (simscore >= 0.6 && simscore < 0.7)
                {eps67 ++; Veps67 ++;}
            else if (simscore >= 0.7 && simscore < 0.8)
                {eps78 ++; Veps78 ++;}
            else if (simscore >= 0.8 && simscore < 0.9)
                {eps89 ++; Veps89 ++;}
            else if (simscore >= 0.9 && simscore < 1)
                {eps910 ++; Veps910 ++;}
            else{
                assert(simscore == 1);
                eps111 ++; Veps111 ++;
            }
        }
    }
        
    cout<<endl<<"V: 2hop deg = "<<V_2hop_nei_cnt/n2<<" ["<<V_min_2hop_deg<<", "<<V_max_2hop_deg<<"]"<<endl;
    cout<<"V: 2hop spa = "<<V_2hop_nei_span/n2<<" ["<<V_min_2hop_span<<", "<<V_max_2hop_span<<"]"<<endl;
    cout<<endl;
    
    long long totalcnt = U_2hop_nei_cnt + V_2hop_nei_cnt;
    
    cout<<"total = "<<totalcnt<<endl;
    
    cout<<"0.0        : "<<((double) eps00/totalcnt )*100<<"%"<<endl;
    cout<<"(0.0, 0.1) : "<<((double) eps01/totalcnt )*100<<"%"<<endl;
    cout<<"[0.1, 0.2) : "<<((double) eps12/totalcnt )*100<<"%"<<endl;
    cout<<"[0.2, 0.3) : "<<((double) eps23/totalcnt )*100<<"%"<<endl;
    cout<<"[0.3, 0.4) : "<<((double) eps34/totalcnt )*100<<"%"<<endl;
    cout<<"[0.4, 0.5) : "<<((double) eps45/totalcnt )*100<<"%"<<endl;
    cout<<"[0.5, 0.6) : "<<((double) eps56/totalcnt )*100<<"%"<<endl;
    cout<<"[0.6, 0.7) : "<<((double) eps67/totalcnt )*100<<"%"<<endl;
    cout<<"[0.7, 0.8) : "<<((double) eps78/totalcnt )*100<<"%"<<endl;
    cout<<"[0.8, 0.9) : "<<((double) eps89/totalcnt )*100<<"%"<<endl;
    cout<<"[0.9, 1.0) : "<<((double) eps910/totalcnt )*100<<"%"<<endl;
    cout<<"1.0        : "<<((double) eps111/totalcnt )*100<<"%"<<endl;
    
    cout<<endl;

    cout<<"U_2hop_nei_cnt = "<<U_2hop_nei_cnt<<endl;
    cout<<"0.0        : "<<((double) Ueps00/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"(0.0, 0.1) : "<<((double) Ueps01/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.1, 0.2) : "<<((double) Ueps12/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.2, 0.3) : "<<((double) Ueps23/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.3, 0.4) : "<<((double) Ueps34/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.4, 0.5) : "<<((double) Ueps45/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.5, 0.6) : "<<((double) Ueps56/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.6, 0.7) : "<<((double) Ueps67/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.7, 0.8) : "<<((double) Ueps78/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.8, 0.9) : "<<((double) Ueps89/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.9, 1.0) : "<<((double) Ueps910/U_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"1.0        : "<<((double) Ueps111/U_2hop_nei_cnt )*100<<"%"<<endl;

    cout<<endl;
    cout<<"V_2hop_nei_cnt = "<<V_2hop_nei_cnt<<endl;
    cout<<"0.0        : "<<((double) Veps00/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"(0.0, 0.1) : "<<((double) Veps01/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.1, 0.2) : "<<((double) Veps12/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.2, 0.3) : "<<((double) Veps23/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.3, 0.4) : "<<((double) Veps34/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.4, 0.5) : "<<((double) Veps45/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.5, 0.6) : "<<((double) Veps56/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.6, 0.7) : "<<((double) Veps67/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.7, 0.8) : "<<((double) Veps78/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.8, 0.9) : "<<((double) Veps89/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"[0.9, 1.0) : "<<((double) Veps910/V_2hop_nei_cnt )*100<<"%"<<endl;
    cout<<"1.0        : "<<((double) Veps111/V_2hop_nei_cnt )*100<<"%"<<endl;
    
    delete [] c;
}

inline bool jsub(ui & u, ui & v)
{
    for(auto & e : vsn[v]){
        if(u >= e.s_idx && u <= e.e_idx && e.max_score >= epsilon) return true;
        if(e.s_idx > u) break;
    }
    return false;
}

inline double js(ui & u, ui & v)
{
    if(u == v) return 0;
    int du = degree[u], dv = degree[v];
    if( (double) min(du, dv) / max(du, dv) < epsilon) return 0;
    
    ui idx1 = pstart[u], idx2 = pstart[v];
    ui idx1_end = idx1 + du, idx2_end = idx2 + dv;
    double common = 0;
    
    while (idx1 < idx1_end && idx2 < idx2_end) {
        if(edges[idx1] == edges[idx2]){
            ++ common; ++ idx1; ++ idx2;
        }
        else{
            if(edges[idx1] > edges[idx2]) {
                ++ idx2;
            }
            else {
                ++ idx1;
            }
        }
    }
    return common/(du + dv - common);
}

vector<ui>  get_sim_nei(ui u)
{
    assert(u >= 0 && u < n);
    vector<ui> sim_nei;
    if(u < n1){ //u is from U
        for(ui i = 0; i < n1; i++){
            if(i == u) continue;
            if(js(u,i) >= epsilon) sim_nei.push_back(i);
        }
    }
    else{ // u is from V
        for(ui i = n1; i < n; i++){
            if(i == u) continue;
            if(js(u,i) >= epsilon) sim_nei.push_back(i);
        }
    }
    return sim_nei;
}

vector<ui> get_sim_nei2(ui u)
{
    assert(u >= 0 && u < n);
    
    vector<ui> sim_nei;
    ui * UV_cnt = new ui[n];
    memset(UV_cnt, 0, sizeof(ui)*n);
    vector<ui> may_sim_list;
    
    for(ui i = pstart[u]; i < pstart[u+1]; i++){
        ui v = edges[i];
        for(ui j = pstart[v]; j < pstart[v+1]; j++){
            ui w = edges[j];
            if(w == u) continue;
            if(UV_cnt[w] == 0) {
                may_sim_list.push_back(w);
            }
            ++ UV_cnt[w];
        }
    }
    
    for(auto e : may_sim_list){
        
        if(del_ver[e]) continue;
        
        assert(degree[e] >= UV_cnt[e]);
        assert(degree[u] >= UV_cnt[e]);
        double sim_score = (double) UV_cnt[e] / (degree[u] + degree[e] - UV_cnt[e]);
        if(sim_score >= epsilon) sim_nei.push_back(e);
        
    }
    delete [] UV_cnt;
    return sim_nei;
}

void naive_enum_core(vector<ui> CL, vector<ui> CR, vector<ui> PL, vector<ui> PR, vector<ui> QL, vector<ui> QR)
{
    if(PL.size() == 0 && PR.size() == 0){
        if(QL.size() == 0 && QR.size() == 0){
            if(CL.size() >= tau && CR.size() >= tau){
                results.push_back(make_pair(CL, CR));
                if(CL.size() + CR.size() > max_MSBC_size) max_MSBC_size = CL.size() + CR.size();
                if(CL.size() + CR.size() < min_MSBC_size) min_MSBC_size = CL.size() + CR.size();
            }
        }
        return;
    }
    
    for(ui idx = 0; idx < PL.size(); idx++){
        ui next_u = PL[idx];
        
        vector<ui> new_CL = CL;
        new_CL.push_back(next_u);
        
        vector<ui> new_CR = CR;
        
        vector<ui> new_PL;
        for(ui idx2 = idx + 1; idx2 < PL.size(); idx2++){
            ui ver = PL[idx2];
            if(js(next_u, ver) >= epsilon) new_PL.push_back(ver);
        }
        
        vector<ui> new_PR;
        for(ui i = 1; i < PR.size(); i++) assert(PR[i-1] < PR[i]);
        vector<ui> L1;
        for(ui i = pstart[next_u]; i < pstart[next_u+1]; i++) if(!del_ver[edges[i]]) L1.push_back(edges[i]);
        ui p1 = 0, p2 = 0;
        while (p1 < L1.size() && p2 < PR.size()) {
            if(L1[p1] == PR[p2]){
                new_PR.push_back(L1[p1]);
                ++ p1; ++ p2;
            }
            else{
                if(L1[p1] > PR[p2]) ++ p2;
                else ++ p1;
            }
        }
        
        vector<ui> new_QL;
        for(auto e : QL) if(js(e, next_u) >= epsilon) new_QL.push_back(e);
        
        vector<ui> new_QR;
        for(ui i = 1; i < QR.size(); i++) assert(QR[i-1] < QR[i]);
        p1 = 0; p2 = 0;
        while (p1 < L1.size() && p2 < QR.size()) {
            if(L1[p1] == QR[p2]){
                new_QR.push_back(L1[p1]);
                ++ p1; ++ p2;
            }
            else{
                if(L1[p1] > QR[p2]) ++ p2;
                else ++ p1;
            }
        }
        if(new_CL.size() + new_PL.size() >= tau && new_CR.size() + new_PR.size() >= tau)
            naive_enum_core(new_CR, new_CL, new_PR, new_PL, new_QR, new_QL);
        QL.push_back(next_u);
    }
    for(ui idx = 0; idx < PR.size(); idx++){
        ui next_u = PR[idx];
        
        vector<ui> new_CL = CL;
        
        vector<ui> new_CR = CR;
        new_CR.push_back(next_u);
        
        vector<ui> new_PL;
        
        vector<ui> new_PR;
        for(ui idx2 = idx + 1; idx2 < PR.size(); idx2++){
            ui ver = PR[idx2];
            if(js(ver, next_u) >= epsilon) new_PR.push_back(ver);
        }
        
        vector<ui> new_QL;
        for(ui i = 1; i < QL.size(); i++) assert(QL[i-1] < QL[i]);
        vector<ui> L1;
        for(ui i = pstart[next_u]; i < pstart[next_u+1]; i++) if(!del_ver[edges[i]]) L1.push_back(edges[i]);
        ui p1 = 0, p2 = 0;
        while (p1 < L1.size() && p2 < QL.size()) {
            if(L1[p1] == QL[p2]){
                new_QL.push_back(L1[p1]);
                ++ p1; ++ p2;
            }
            else{
                if(L1[p1] > QL[p2]) ++ p2;
                else ++ p1;
            }
        }
        
        vector<ui> new_QR;
        for(auto e : QR) if(js(e, next_u) >= epsilon) new_QR.push_back(e);
        if(new_CL.size() + new_PL.size() >= tau && new_CR.size() + new_PR.size() >= tau)
            naive_enum_core(new_CR, new_CL, new_PR, new_PL, new_QR, new_QL);
        QR.push_back(next_u);
    }
}

void vr()
{
    cout<<"Start VR."<<endl;
    Timer tt;
    assert(tau >= 2);
    
    queue<ui> Q;
    for(ui i = 0; i < n; i++) if(TMPdeg[i] < tau || get_sim_nei2(i).size() < tau - 1) {
        Q.push(i);
    }
    int delcnt = 0;
    while (!Q.empty()) {
        
        ui u = Q.front();
        Q.pop();
        del_ver[u] = 1;
        
        delcnt++;
        
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            -- TMPdeg[v];
            if(TMPdeg[v] == tau - 1 && get_sim_nei2(v).size() >= tau - 1) {
                Q.push(v);
            }
        }

        vector<ui> sim_nei = get_sim_nei2(u);
        for(ui w : sim_nei){
            assert(!del_ver[w]);
            if(get_sim_nei2(w).size() == tau - 2 && TMPdeg[w] >= tau) {
                Q.push(w);
            }
        }
        
    }
    int cnt = 0;
    for(ui i = 0; i < n; i++) if(!del_ver[i]) ++ cnt;
    for(ui i = 0; i < n; i++) if(!del_ver[i]) Sdeg[i] = get_sim_nei2(i).size();
    cout<<"Delete "<<delcnt<<" vertices. (remain "<<cnt<<")"<<endl;
}

void index_vr()
{
    if(vsn.empty()) {cout<<"vsn is empty."<<endl; exit(1);}
    assert(tau >= 2);
    cout<<"Start index based VR."<<endl;
    
    Timer tt;
    queue<ui> Q;
    for(ui i = 0; i < n; i++) if(TMPdeg[i] < tau) Q.push(i);
    
    int delcnt = 0;
    
    while (!Q.empty()) {
        ui u = Q.front(); Q.pop();
        del_ver[u] = 1;
        ++ delcnt;
        for(ui i = pstart[u]; i < pstart[u+1]; i++) if(TMPdeg[edges[i]]-- == tau) Q.push(edges[i]);
    }
    
    int sd;
    for(ui u = 0; u < n; u++) if(!del_ver[u]) {
        const vector<Itval> & cand_list = vsn[u];
        sd = 0;
        for(auto & e : cand_list) if(e.max_score >= epsilon) sd += e.c;
        if(sd < tau - 1) {
            Q.push(u);
            del_ver[u] = 1;
        }
    }
        
    while (!Q.empty()) {
        ui u = Q.front(); Q.pop();
        ++ delcnt;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            if(TMPdeg[v]-- == tau && !del_ver[v]) {
                Q.push(v);
                del_ver[v] = 1;
            }
        }
    }
    cout<<"     Phase-I, deleted "<<delcnt<<", (remain "<<n-delcnt<<"), time cost = "<<integer_to_string(tt.elapsed())<<endl;
    tt.restart();
            
    for(ui u = 0; u < n; u++) if(!del_ver[u]) {
        assert(TMPdeg[u] >= tau);
        vector<Itval> & cand_list = vsn[u];
        for(auto & e : cand_list) if(e.max_score >= epsilon) {
#ifdef _SHRINK_
            bool f = false;
#endif
            for(ui i = e.s_idx; i <= e.e_idx; i++) if(!del_ver[i] && js(u, i) >= epsilon && i != u) {
                ++ Sdeg[u];
#ifdef _SHRINK_
                f = true;
#endif
            }
#ifdef _SHRINK_
            if(!f) {
                e.max_score = 0;
                e.min_score = 0;
                e.c = 0;
            }
#endif
        }
        if(Sdeg[u] < tau - 1) Q.push(u);
    }
    
    int phase2_total_del = 0;

    while (!Q.empty()) {
        ui u = Q.front(); Q.pop();
        del_ver[u] = 1; ++ delcnt; ++ phase2_total_del;
        
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            if(del_ver[v]) continue;
            if(TMPdeg[v] -- == tau && Sdeg[v] >= tau - 1) {
                Q.push(v);
            }
        }
        
        vector<Itval> & cand_list = vsn[u];
        for(auto & e : cand_list){
            if(e.max_score < epsilon) continue;
#ifdef _SHRINK_
            bool f = false;
#endif
            for(ui i = e.s_idx; i <= e.e_idx; i++){
                if(del_ver[i]) continue;
                if(js(u, i) >= epsilon && i != u) {
                    if(Sdeg[i] -- == tau - 1 && TMPdeg[i] >= tau) {
                        Q.push(i);
                    }
#ifdef _SHRINK_
                    f = true;
#endif
                }
            }
#ifdef _SHRINK_
            if(!f) {
                e.max_score = 0;
                e.min_score = 0;
                e.c = 0;
            }
#endif
        }
    }
    cout<<"     Phase-II, deleted "<<phase2_total_del<<", (remain "<<n-delcnt<<"), time cost = "<<integer_to_string(tt.elapsed())<<endl;
}

void materialize_the_remaining_signed_graph()
{
    cout<<"** in materialize_the_remaining_signed_graph() **"<<endl;
    unordered_map<ui, ui> newid;
    ui idx = 0;
    for(ui i = 0; i < n; i++) if(del_ver[i] == 0) newid[i] = idx ++;
    
    map<ui, map<ui, int>> sG;
    
    for(ui u = 0; u < n; u++) if(del_ver[u] == 0) {
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            if(del_ver[v] == 0){
                sG[u][v] = -1;
                sG[v][u] = -1;
            }
        }
    }
    
    ui * UV_cnt = new ui[n];
    memset(UV_cnt, 0, sizeof(ui)*n);
    
    for(ui u = 0; u < n; u++) if(del_ver[u] == 0) {
        vector<ui> may_sim_list;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(del_ver[w]) continue;
                if(w == u) continue;
                if(UV_cnt[w] == 0){
                    may_sim_list.push_back(w);
                }
                ++ UV_cnt[w];
            }
        }
                
        for(auto e : may_sim_list){
            assert(degree[e] >= UV_cnt[e]);
            assert(degree[u] >= UV_cnt[e]);
            double sim_score = (double) UV_cnt[e] / (degree[e] + degree[u] - UV_cnt[e]);
            if(sim_score >= epsilon){
                sG[u][e] = 1;
                sG[e][u] = 1;
            }
            UV_cnt[e] = 0;
        }
    }
    delete [] UV_cnt;
    
    cout<<"the materialized signed graph : ";
    ui sg_n = sG.size();
    ui sg_m = 0;
    for(auto e : sG) sg_m += e.second.size();
    assert(sg_m%2 == 0);
    sg_m /= 2;
    cout<<"sg_n = "<<sg_n<<", sg_m = "<<sg_m<<endl;
    
    ofstream fout;
    string graph_name = "/Users/yaokai/datasets/MSBC/mat-datasets/mat_remain_";
    int int_epsilon = epsilon * 100;
    graph_name.append(to_string(int_epsilon));
    graph_name.append(".txt");
    //cout<<graph_name<<endl;
    fout.open(graph_name);
    assert(fout.is_open());
    
    fout<<sg_n<<"\t"<<sg_m;
    for(auto e : sG){
        for(auto x : e.second){
            fout<<endl<<newid[e.first]<<"\t"<<newid[x.first]<<"\t"<<x.second;
        }
    }
    fout.close();
    cout<<"finish writing the materialized remaining signed graph to disk!"<<endl;
}

void naive_enum()
{
    cout<<"     ** in naive_enum() **"<<endl;
    for(ui u = 0; u < n; u++){
        if(del_ver[u]) continue;
        vector<ui> CL;
        CL.push_back(u);
        
        vector<ui> CR;
        
        vector<ui> PL, QL;
        vector<ui> tmp_vec;
        tmp_vec = get_sim_nei2(u);
        sort(tmp_vec.begin(), tmp_vec.end(), less<>());
        for(auto e : tmp_vec){
            if(e > u) PL.push_back(e);
            else QL.push_back(e);
        }
        
        vector<ui> PR, QR;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            
            if(del_ver[v]) continue;
            
            if(v > u) PR.push_back(v);
            else QR.push_back(v);
        }
        if(PL.size() < tau - 1 || PR.size() < tau) {
            continue;
        }
        naive_enum_core(CL, CR, PL, PR, QL, QR);
    }
}

void peel_remaining_graph_by_connection_degree()
{
    r_n = 0;
    r_m = 0;
    for(ui u = 0; u < n; u++) if(del_ver[u] == 0) {
        ++ r_n;
        r_m += TMPdeg[u];
    }
    assert(r_m%2 == 0); r_m /= 2;
    
    if(r_n == 0) return;
    
    peel_s = new ui[r_n];
    ui idx = 0;
    for(ui u = 0; u < n; u++) if(del_ver[u] == 0) peel_s[idx ++] = u;
    assert(idx == r_n);
    
    ui * deg = new ui[n];
    for(ui u = 0; u < n; u++) deg[u] = TMPdeg[u];
    
    ui * core = new ui[n];
    memset(core, 0, sizeof(ui)*n);
    ui max_core = 0;
    
    ListLinearHeap *linear_heap = new ListLinearHeap(n, n-1);
    linear_heap->init(idx, idx-1, peel_s, deg);
    for(ui i = 0; i < idx; i ++) {
        ui u, key;
        linear_heap->pop_min(u, key);
        if(key > max_core) max_core = key;
        peel_s[i] = u;
        core[u] = max_core;
        for(ui j = pstart[u];j < pstart[u+1];j ++) if(del_ver[edges[j]] == 0 && core[edges[j]] == 0) {
            linear_heap->decrement(edges[j], 1);
        }
    }
    delete linear_heap;
    delete [] core;
    delete [] deg;
}

void rebuild_remaining_graph()
{
    oid = new ui[r_n];
    nid = new ui[n];
    
    for(ui u = 0; u < r_n; u++){
        nid[peel_s[u]] = u;
        oid[u] = peel_s[u];
    }
    
    r_pstart = new ui[r_n+1];
    r_edges = new ui[2*r_m];
    r_degree = new int[r_n];
    
    r_pstart[0] = 0;
    ui r_idx = 0;
    for(ui u = 0; u < r_n; u++){
        ui start_pos = r_pstart[r_idx];
        ui tdeg = 0;
        ui orivid = oid[u];
        for(ui i = pstart[orivid]; i < pstart[orivid+1]; i++) if(del_ver[edges[i]] == 0) {
            r_edges[start_pos++] = nid[edges[i]];
            ++ tdeg;
        }
        r_degree[u] = tdeg;
        r_pstart[++r_idx] = start_pos;
    }
    
    assert(r_idx == r_n);
}

void adv_enum_core(vector<ui> CL, vector<ui> CR, vector<ui> PL, vector<ui> PR, vector<ui> QL, vector<ui> QR)
{
    if(PL.empty() && PR.empty()){
        if(QL.empty() && QR.empty() && CL.size() >= tau && CR.size() >= tau){
            results.push_back(make_pair(CL, CR));
            if(CL.size() + CR.size() > max_MSBC_size) max_MSBC_size = CL.size() + CR.size();
            if(CL.size() + CR.size() < min_MSBC_size) min_MSBC_size = CL.size() + CR.size();
        }
        return;
    }
    
    if(CL.size() + PL.size() < tau || CR.size() + PR.size() < tau) return;
    
    
    for(ui idx = 0; idx < PL.size(); idx++){
        ui next_u = PL[idx];
        vector<ui> new_CL, new_CR, new_PL, new_PR, new_QL, new_QR;
        //new_CL
        new_CL = CL;
        new_CL.push_back(next_u);
        //new_CR
        new_CR = CR;
        //new_PL
        for(ui idx2 = idx + 1; idx2 < PL.size(); idx2++) if(js(oid[next_u], oid[PL[idx2]]) >= epsilon) new_PL.push_back(PL[idx2]);
        //new_PR
        ui p1s = 0, p1e = PR.size(), p2s = r_pstart[next_u], p2e = r_pstart[next_u+1];
        while (p1s < p1e && p2s < p2e) {
            if(PR[p1s] == r_edges[p2s]){
                new_PR.push_back(PR[p1s]);
                ++ p1s; ++ p2s;
            }
            else{
                if(PR[p1s] < r_edges[p2s]) ++ p1s;
                else ++ p2s;
            }
        }
        //new_QL
        for(auto e : QL) if(js(oid[next_u], oid[e]) >= epsilon) new_QL.push_back(e);
        //new_QR
        p1s = 0, p1e = QR.size(), p2s = r_pstart[next_u];
        while (p1s < p1e && p2s < p2e) {
            if(QR[p1s] == r_edges[p2s]){
                new_QR.push_back(QR[p1s]);
                ++ p1s; ++ p2s;
            }
            else{
                if(QR[p1s] < r_edges[p2s]) ++ p1s;
                else ++ p2s;
            }
        }
        
        if(new_CL.size() + new_PL.size() >= tau && new_CR.size() + new_PR.size() >= tau){
            adv_enum_core(new_CR, new_CL, new_PR, new_PL, new_QR, new_QL);
        }
        QL.push_back(next_u);
    }
    
    for(ui idx = 0; idx < PR.size(); idx++){
        ui next_u = PR[idx];
        vector<ui> new_CL, new_CR, new_PL, new_PR, new_QL, new_QR;
        //new_CL
        new_CL = CL;
        //new_CR
        new_CR = CR;
        new_CR.push_back(next_u);
        //new_PR
        for(ui idx2 = idx + 1; idx2 < PR.size(); idx2 ++) if(js(oid[next_u], oid[PR[idx2]]) >= epsilon) new_PR.push_back(PR[idx2]);
        //new_QL
        ui p1s = 0, p1e = QL.size(), p2s = r_pstart[next_u], p2e = r_pstart[next_u+1];
        while (p1s < p1e && p2s < p2e) {
            if(QL[p1s] == r_edges[p2s]){
                new_QL.push_back(QL[p1s]);
                ++ p1s; ++ p2s;
            }
            else{
                if(QL[p1s] < r_edges[p2s]) ++ p1s;
                else ++ p2s;
            }
        }
        //new_QR
        for(auto e : QR) if(js(oid[e], oid[next_u]) >= epsilon) new_QR.push_back(e);
        if(new_CL.size() + new_PL.size() >= tau && new_CR.size() + new_PR.size() >= tau){
            adv_enum_core(new_CR, new_CL, new_PR, new_PL, new_QR, new_QL);
        }
        QR.push_back(next_u);
    }
}

void build_matrix_for_PQ(vector<ui> &QL, vector<ui> &PL, vector<ui> &QR, vector<ui> &PR)
{
    for(ui i = 0; i < max_CS_deg; i++) memset(Matrix[i], 0, sizeof(int)*max_CS_deg);
    
    int idx = 0;
    for(auto e : QL) trans[e] = idx ++;
    for(auto e : PL) trans[e] = idx ++;
    for(auto e : QR) trans[e] = idx ++;
    for(auto e : PR) trans[e] = idx ++;
    
    long long num_of_edges = 0;
    
    //L
    for(ui i = 0; i < QL.size(); i++){
        for(ui j = i + 1; j < QL.size(); j++){
            ui u = QL[i], v = QL[j];
            ++ total_js;
#ifdef _JSUB_
            if(!jsub(oid[u], oid[v])) {++ cnt_prune_by_ub; continue;}
#endif
            ++ cnt_cmp_by_js;
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
                ++ num_of_edges;
            }
        }
    }
    for(auto u : QL){
        for(auto v : PL){
            ++ total_js;
#ifdef _JSUB_
            if(!jsub(oid[u], oid[v])) {++ cnt_prune_by_ub; continue;}
#endif
            ++ cnt_cmp_by_js;
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
                ++ num_of_edges;
            }
        }
    }
    for(ui i = 0; i < PL.size(); i++){
        for(ui j = i + 1; j < PL.size(); j++){
            ui u = PL[i], v = PL[j];
            ++ total_js;
#ifdef _JSUB_
            if(!jsub(oid[u], oid[v])) {++ cnt_prune_by_ub; continue;}
#endif
            ++ cnt_cmp_by_js;
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
                ++ num_of_edges;
            }
        }
    }
    
    //R
    for(ui i = 0; i < QR.size(); i++){
        for(ui j = i + 1; j < QR.size(); j++){
            ui u = QR[i], v = QR[j];
            ++ total_js;
#ifdef _JSUB_
            if(!jsub(oid[u], oid[v])) {++ cnt_prune_by_ub; continue;}
#endif
            ++ cnt_cmp_by_js;
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
                ++ num_of_edges;
            }
        }
    }
    for(auto u : QR){
        for(auto v : PR){
            ++ total_js;
#ifdef _JSUB_
            if(!jsub(oid[u], oid[v])) {++ cnt_prune_by_ub; continue;}
#endif
            ++ cnt_cmp_by_js;
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
                ++ num_of_edges;
            }
        }
    }
    for(ui i = 0; i < PR.size(); i++){
        for(ui j = i + 1; j < PR.size(); j++){
            ui u = PR[i], v = PR[j];
            ++ total_js;
#ifdef _JSUB_
            if(!jsub(oid[u], oid[v])) {++ cnt_prune_by_ub; continue;}
#endif
            ++ cnt_cmp_by_js;
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
                ++ num_of_edges;
            }
        }
    }
    for(auto u : QL){
        for(ui i = r_pstart[u]; i < r_pstart[u+1]; i++){
            ui v = r_edges[i];
            if(inPQ[v] != 0) {
                Matrix[trans[u]][trans[v]] = -1;
                Matrix[trans[v]][trans[u]] = -1;
                ++ num_of_edges;
            }
        }
    }
    for(auto u : PL){
        for(ui i = r_pstart[u]; i < r_pstart[u+1]; i++){
            ui v = r_edges[i];
            if(inPQ[v] != 0) {
                Matrix[trans[u]][trans[v]] = -1;
                Matrix[trans[v]][trans[u]] = -1;
                ++ num_of_edges;
            }
        }
    }

    num_of_edges *= 2;
    
    if(idx > 1) {
        double den = ((double)num_of_edges/((idx)*(idx - 1)));
        if(den > max_density) max_density = den;
        if(den < min_density) min_density = den;
        ave_density += den;
        ++ ave_density_cnt;
    }
}

void obtain_degree_inP(vector<ui> &PL, vector<ui> &PR)
{
    for(auto e : PL) {
        pdeg_inP[e] = 0;
        ndeg_inP[e] = 0;
    }
    for(auto e : PR) {
        pdeg_inP[e] = 0;
        ndeg_inP[e] = 0;
    }
    
    for(ui i = 0; i < PL.size(); i++){
        for(ui j = i + 1; j < PL.size(); j++){
            ui u = PL[i], v = PL[j];
            if(Matrix[trans[u]][trans[v]] != 0){
                ++ pdeg_inP[u]; ++ pdeg_inP[v];
            }
        }
    }
    for(ui i = 0; i < PR.size(); i++){
        for(ui j = i + 1; j < PR.size(); j++){
            ui u = PR[i], v = PR[j];
            if(Matrix[trans[u]][trans[v]] != 0){
                ++ pdeg_inP[u]; ++ pdeg_inP[v];
            }
        }
    }
    for(auto u : PL){
        for(auto v : PR){
            if(Matrix[trans[u]][trans[v]] != 0){
                ++ ndeg_inP[u]; ++ ndeg_inP[v];
            }
        }
    }
}

void pruneP_by_degree(vector<ui> &CL, vector<ui> &CR, vector<ui> &PL, vector<ui> &PR)
{
    int L_pt = tau - (int)CL.size() - 1;
    int L_nt = tau - (int)CR.size();
    int R_pt = tau - (int)CR.size() - 1;
    int R_nt = tau - (int)CL.size();
    if(L_pt <= 0 && L_nt <= 0 && R_pt <= 0 && R_nt <= 0) return;
    queue<ui> Q;
    for(auto e : PL) if(pdeg_inP[e] < L_pt || ndeg_inP[e] < L_nt) {
        Q.push(e);
    }
    for(auto e : PR) if(pdeg_inP[e] < R_pt || ndeg_inP[e] < R_nt) {
        Q.push(e);
    }
    
    ui * x = new ui[r_n];
    memset(x, 0, sizeof(ui) * r_n);
    
    while (!Q.empty()) {
        ui u = Q.front();
        Q.pop();
        x[u] = 1;
        for(auto v : PL){
            if(Matrix[trans[u]][trans[v]] == 1){
                assert(pdeg_inP[v] > 0);
                if(pdeg_inP[v] -- == L_pt && ndeg_inP[v] >= L_nt) {
                    Q.push(v);
                }
            }
            else if (Matrix[trans[u]][trans[v]] == -1){
                assert(ndeg_inP[v] > 0);
                if(ndeg_inP[v] -- == L_nt && pdeg_inP[v] >= L_pt) {
                    Q.push(v);
                }
            }
        }
        for(auto v : PR){
            if(Matrix[trans[u]][trans[v]] == 1){
                assert(pdeg_inP[v] > 0);
                if(pdeg_inP[v] -- == R_pt && ndeg_inP[v] >= R_nt) {
                    Q.push(v);
                }
            }
            else if (Matrix[trans[u]][trans[v]] == -1){
                assert(ndeg_inP[v] > 0);
                if(ndeg_inP[v] -- == R_nt && pdeg_inP[v] >= R_pt) {
                    Q.push(v);
                }
            }
        }
    }
    vector<ui> newPL, newPR;
    for(auto e : PL) if(x[e] == 0) newPL.push_back(e);
    for(auto e : PR) if(x[e] == 0) newPR.push_back(e);
    PL = newPL;
    PR = newPR;
    delete [] x;
    x = nullptr;
}

ui pivot_selection(vector<ui> &PL, vector<ui> &PR, vector<ui> &QL, vector<ui> &QR)
{
    for(auto e : PL) deg_inP[e] = pdeg_inP[e] + ndeg_inP[e];
    for(auto e : PR) deg_inP[e] = pdeg_inP[e] + ndeg_inP[e];
    
    for(auto e : QL) deg_inP[e] = 0;
    for(auto e : QR) deg_inP[e] = 0;
    for(auto u : QL){
        for(auto v : PL) if(Matrix[trans[u]][trans[v]] != 0) ++ deg_inP[u];
        for(auto v : PR) if(Matrix[trans[u]][trans[v]] != 0) ++ deg_inP[u];
    }
    for(auto u : QR){
        for(auto v : PL) if(Matrix[trans[u]][trans[v]] != 0) ++ deg_inP[u];
        for(auto v : PR) if(Matrix[trans[u]][trans[v]] != 0) ++ deg_inP[u];
    }
    
    ui pivot = 0;
    int max_local_deg = -1;
    for(auto e : PL) if(deg_inP[e] > max_local_deg) {max_local_deg = deg_inP[e]; pivot = e; }
    for(auto e : PR) if(deg_inP[e] > max_local_deg) {max_local_deg = deg_inP[e]; pivot = e; }
    for(auto e : QL) if(deg_inP[e] > max_local_deg) {max_local_deg = deg_inP[e]; pivot = e; }
    for(auto e : QR) if(deg_inP[e] > max_local_deg) {max_local_deg = deg_inP[e]; pivot = e; }
    
    assert(max_local_deg >= 0);
    
    return pivot;
}

void matrix_based_enum_core(vector<ui> CL, vector<ui> CR, vector<ui> PL, vector<ui> PR, vector<ui> QL, vector<ui> QR)
{
    if(over_time_flag) return;
    if(((double)clock() / CLOCKS_PER_SEC - startTime) > maxTime) over_time_flag = true;
    
    if(PL.empty() && PR.empty()){
        if(QL.empty() && QR.empty()){
            if(CL.size() >= tau && CR.size() >= tau){
                
#ifdef _CheckResults_
                vector<ui> C1, C2;
                for(auto e : CL) C1.push_back(oid[e]);
                for(auto e : CR) C2.push_back(oid[e]);
                sort(C1.begin(), C1.end(), less<>());
                sort(C2.begin(), C2.end(), less<>());
                if(C1[0] < C2[0]) results.push_back(make_pair(C1, C2));
                else results.push_back(make_pair(C2, C1));
#endif
                ++ MSBC_num;
                if(CL.size() + CR.size() > max_MSBC_size) max_MSBC_size = CL.size() + CR.size();
                if(CL.size() + CR.size() < min_MSBC_size) min_MSBC_size = CL.size() + CR.size();
            }
        }
        return;
    }
    if(CL.size() + PL.size() < tau || CR.size() + PR.size() < tau) return;
    obtain_degree_inP(PL, PR);  //validate pdeg_inP[] and ndeg_inP[], which will be use for P pruning.
    pruneP_by_degree(CL, CR, PL, PR);
    if(PL.empty() && PR.empty()) return;
    if(CL.size() + PL.size() < tau || CR.size() + PR.size() < tau) return;
    
    ui pivot = pivot_selection(PL, PR, QL, QR);
    
    sort(PL.begin(), PL.end(), cmp_degP);
    sort(PR.begin(), PR.end(), cmp_degP);
    
    ui PLsize = (ui)PL.size();
    ui PRsize = (ui)PR.size();
    
    vector<int> PL_pivot(PLsize, 1);
    vector<int> PR_pivot(PRsize, 1);
    for(ui i = 0; i < PL.size(); i++) if(Matrix[trans[PL[i]]][trans[pivot]] != 0) PL_pivot[i] = 0;
    for(ui i = 0; i < PR.size(); i++) if(Matrix[trans[PR[i]]][trans[pivot]] != 0) PR_pivot[i] = 0;
    
    vector<int> PL_exp(PLsize, 1);
    vector<int> PR_exp(PRsize, 1);
    
    for(ui i = 0; i < PL.size(); i++){
        
        if(PL_pivot[i] == 0) continue;
        
        ui u = PL[i];
        
        vector<ui> newCL = CL;
        newCL.push_back(u);
        
        vector<ui> newCR = CR;
        
        vector<ui> newPL;
        for(ui j = 0; j < PL.size(); j++) if(Matrix[trans[PL[j]]][trans[u]] != 0 && PL_exp[j] != 0)
            newPL.push_back(PL[j]);
        
        vector<ui> newPR;
        for(ui j = 0; j < PR.size(); j++) if(Matrix[trans[PR[j]]][trans[u]] != 0)
            newPR.push_back(PR[j]);
        
        vector<ui> newQL;
        for(ui j = 0; j < QL.size(); j++) if(Matrix[trans[QL[j]]][trans[u]] != 0)
            newQL.push_back(QL[j]);
        
        vector<ui> newQR;
        for(ui j = 0; j < QR.size(); j++) if(Matrix[trans[QR[j]]][trans[u]] != 0)
            newQR.push_back(QR[j]);
        
        matrix_based_enum_core(newCR, newCL, newPR, newPL, newQR, newQL);  //swap L and R
        
        PL_exp[i] = 0;
        QL.push_back(u);
    }
    
    for(ui i = 0; i < PR.size(); i++){
        
        if(PR_pivot[i] == 0) continue;
        
        ui u = PR[i];
        
        vector<ui> newCL = CL;
        
        vector<ui> newCR = CR;
        newCR.push_back(u);
        
        vector<ui> newPL;
        for(ui j = 0; j < PL.size(); j++) if(Matrix[trans[PL[j]]][trans[u]] != 0 && PL_exp[j] != 0)
            newPL.push_back(PL[j]);
        
        vector<ui> newPR;
        for(ui j = 0; j < PR.size(); j++) if(Matrix[trans[PR[j]]][trans[u]] != 0 && PR_exp[j] != 0)
            newPR.push_back(PR[j]);
        
        vector<ui> newQL;
        for(ui j = 0; j < QL.size(); j++) if(Matrix[trans[QL[j]]][trans[u]] != 0)
            newQL.push_back(QL[j]);
        
        vector<ui> newQR;
        for(ui j = 0; j < QR.size(); j++) if(Matrix[trans[QR[j]]][trans[u]] != 0)
            newQR.push_back(QR[j]);
        
        matrix_based_enum_core(newCR, newCL, newPR, newPL, newQR, newQL);  //swap L and R
        
        PR_exp[i] = 0;
        QR.push_back(u);
    }
    
}

void matrix_based_enum()
{
    cout<<"Start enumeration procedure."<<endl;
    Timer t;
    peel_remaining_graph_by_connection_degree();  //get an approximate core decomposition
    t.restart();
    
    if(r_n == 0) return;
    
    t.restart();
    rebuild_remaining_graph();  //rebuild the remaining graph to r_pstart[] and r_edges[]
    t.restart();

    max_CS_deg = 0;
    for(ui u = 0; u < r_n; u++) if(TMPdeg[oid[u]] + Sdeg[oid[u]] > max_CS_deg) {
         max_CS_deg = TMPdeg[oid[u]] + Sdeg[oid[u]];
    }
    
    t.restart();
    
    Matrix = new int * [max_CS_deg];
    for(ui i = 0; i < max_CS_deg; i++) Matrix[i] = new int [max_CS_deg];
    trans = new ui[r_n];
    inPQ = new int[r_n];  //in PL/PR : 1, in QL/QR : -1
    memset(inPQ, 0, sizeof(int)*r_n);
    
    deg_inP = new int[r_n];
    pdeg_inP = new int[r_n];
    ndeg_inP = new int[r_n];
    
    long long enum_findCPQ_T = 0;
    long long enum_build_sg_T = 0;
    long long enum_core_T = 0;
    long long last_MSBC_num = 0;
    
    long long num_of_I = 0;
    
    Timer tt;
    startTime = (double)clock() / CLOCKS_PER_SEC;
    
    min_density = INF;
    max_density = 0;
    ave_density = 0;
    ave_density_cnt = 0;
    
    vector<ui> all_new_vertices_inR;
    for(ui u = 0; u < r_n; u++) {
        if(os[oid[u]] == 2) {
            all_new_vertices_inR.push_back(u);
        }
    }
    
    for(ui u = 0; u < r_n; u++){
            
        tt.restart();
        
        vector<ui> CL;
        CL.push_back(u);
        
        vector<ui> CR;
        
        vector<ui> PL, QL;
        
        if (vr_method == 1) { //vr_method == 1, no index is loaded.
            vector<ui> vec = get_sim_nei2(oid[u]);
            for(auto & e : vec) {
                ui v = nid[e];
                if(v > u) {
                    PL.push_back(v);
                    inPQ[v] = 1;
                }
                else if (v < u) {
                    QL.push_back(v);
                    inPQ[v] = -1;
                }
            }
        }
        else {
            vector<Itval> & cand_list = vsn[oid[u]];
            for(auto & e : cand_list) if(e.max_score >= epsilon) {
                for(ui i = e.s_idx; i <= e.e_idx; i++) if(del_ver[i] == 0) {
                    if(js(oid[u], i) >= epsilon){
                        ui v = nid[i];
                        if(v > u) {
                            PL.push_back(v);
                            inPQ[v] = 1;
                        }
                        else if(v < u) {
                            QL.push_back(v);
                            inPQ[v] = -1;
                        }
                    }
                }
                ++ num_of_I;
            }
        }
        
        vector<ui> PR, QR;
        for(ui i = r_pstart[u]; i < r_pstart[u+1]; i++){
            ui v = r_edges[i];
            if(v > u) {
                PR.push_back(v);
                inPQ[v] = 1;
            }
            else if (v < u) {
                QR.push_back(v);
                inPQ[v] = -1;
            }
        }
        
        enum_findCPQ_T += tt.elapsed();
    
        if(PL.size() < tau - 1 || PR.size() < tau) {
            for(auto e : QL) inPQ[e] = 0;
            for(auto e : PL) inPQ[e] = 0;
            for(auto e : QR) inPQ[e] = 0;
            for(auto e : PR) inPQ[e] = 0;
            continue;
        }
    
        tt.restart();
        build_matrix_for_PQ(QL, PL, QR, PR);  //construct matrix
        enum_build_sg_T += tt.elapsed();
        
        last_MSBC_num = MSBC_num;
        tt.restart();
        matrix_based_enum_core(CL, CR, PL, PR, QL, QR);
        enum_core_T += tt.elapsed();
                
        for(auto e : QL) inPQ[e] = 0;
        for(auto e : PL) inPQ[e] = 0;
        for(auto e : QR) inPQ[e] = 0;
        for(auto e : PR) inPQ[e] = 0;
    
        if (over_time_flag) break;
    }
    long long total_T = enum_findCPQ_T + enum_build_sg_T + enum_core_T;
    cout<<"     -- Time on find CPQ   = "<<integer_to_string(enum_findCPQ_T)<<"    ("<<((double)enum_findCPQ_T/total_T)*100<<"%)."<<endl;
    cout<<"     -- Time on build sg   = "<<integer_to_string(enum_build_sg_T)<<"    ("<<((double)enum_build_sg_T/total_T)*100<<"%)."<<endl;
    cout<<"     -- Time on enum core = "<<integer_to_string(enum_core_T)<<"    ("<<((double)enum_core_T/total_T)*100<<"%)."<<endl;
    cout<<"     As for all subgraphs, ave density = "<<(double)ave_density/ave_density_cnt<<" ("<<min_density<<" -> "<<max_density<<")."<<endl;
}

void MSBC_Enum()
{
    Timer t;
    
    if(vr_method == 1) vr();  //vertex reduction
    else if (vr_method == 2) index_vr();  //index based vertex reduction
    else cout<<"No matching vr."<<endl;
    
    time_on_vr += t.elapsed();
    t.restart();
    
    matrix_based_enum();  //enumeration procedure
    
    time_on_enum += t.elapsed();
}

void materialize_graph_only_count()
{
    cout<<"*** in mat only count ***"<<endl;

    ui * UV_cnt = new ui[n];
    
    long long eps0 = 0;
    long long eps1 = 0;
    long long eps2 = 0;
    long long eps3 = 0;
    long long eps5 = 0;
    long long eps6 = 0;
    long long eps8 = 0;
    long long eps10 = 0;
    
    for(ui u = 0; u < n; u++){
        memset(UV_cnt, 0, sizeof(ui)*n);
        vector<ui> may_sim_list;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(UV_cnt[w] == 0){
                    may_sim_list.push_back(w);
                }
                ++ UV_cnt[w];
            }
        }
        
        for(auto e : may_sim_list){
            assert(degree[e] >= UV_cnt[e]);
            assert(degree[u] >= UV_cnt[e]);
            double sim_score = (double) UV_cnt[e] / (degree[e] + degree[u] - UV_cnt[e]);
            ++ eps0;
            if(sim_score >= 0.1) ++ eps1;
            if(sim_score >= 0.2) ++ eps2;
            if(sim_score >= 0.3) ++ eps3;
            if(sim_score >= 0.5) ++ eps5;
            if(sim_score >= 0.6) ++ eps6;
            if(sim_score >= 0.8) ++ eps8;
            if(sim_score >= 1) ++ eps10;
        }
    }
    delete [] UV_cnt;
    assert(eps0%2 == 0);
    assert(eps1%2 == 0);
    assert(eps2%2 == 0);
    assert(eps3%2 == 0);
    assert(eps5%2 == 0);
    assert(eps6%2 == 0);
    assert(eps8%2 == 0);
    assert(eps10%2 == 0);
    cout<<endl;
    cout<<"     eps0 = "<<eps0<<", total edges = "<<eps0+m<<endl;
    cout<<"     eps1 = "<<eps1<<", total edges = "<<eps1+m<<endl;
    cout<<"     eps2 = "<<eps2<<", total edges = "<<eps2+m<<endl;
    cout<<"     eps3 = "<<eps3<<", total edges = "<<eps3+m<<endl;
    cout<<"     eps5 = "<<eps5<<", total edges = "<<eps5+m<<endl;
    cout<<"     eps6 = "<<eps6<<", total edges = "<<eps6+m<<endl;
    cout<<"     eps8 = "<<eps8<<", total edges = "<<eps8+m<<endl;
    cout<<"     eps10 = "<<eps10<<", total edges = "<<eps10+m<<endl;
}

void materialize_the_signed_graph(string graph_name, bool write_the_materialized_signed_graph_on_disk)
{
    Timer t;
    
    cout<<"*** in mat ***"<<endl;
    map<ui, map<ui, int>> sG;
    for(ui i = 0; i < n; i++){
        for(ui j = pstart[i]; j < pstart[i+1]; j++){
            sG[i][edges[j]] = -1;
            sG[edges[j]][i] = -1;
        }
    }

    ui * UV_cnt = new ui[n];
    
    for(ui u = 0; u < n; u++){
        memset(UV_cnt, 0, sizeof(ui)*n);
        vector<ui> may_sim_list;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(UV_cnt[w] == 0){
                    may_sim_list.push_back(w);
                }
                ++ UV_cnt[w];
            }
        }
        
        for(auto e : may_sim_list){
            assert(degree[e] >= UV_cnt[e]);
            assert(degree[u] >= UV_cnt[e]);
            double sim_score = (double) UV_cnt[e] / (degree[e] + degree[u] - UV_cnt[e]);
            if(sim_score >= epsilon){
                sG[u][e] = 1;
                sG[e][u] = 1;
            }
        }
    }
    delete [] UV_cnt;
    ui sg_n = sG.size();
    assert(sg_n == n);
    ui sg_m = 0;
    for(auto e : sG) sg_m += e.second.size();
    assert(sg_m%2 == 0);
    sg_m /= 2;
    cout<<endl<<"sg_n = "<<sg_n<<", sg_m = "<<sg_m<<endl;
    assert(sg_m >= m);
    int U_new_edges = 0;
    for(ui i = 0; i < n1; i++){
        assert(sG.find(i) != sG.end());
        for(auto e : sG[i]) if(e.second == 1) ++ U_new_edges;
    }
    assert(U_new_edges%2 == 0); U_new_edges /= 2;
    int V_new_edges = 0;
    for(ui i = n1; i < n; i++){
        assert(sG.find(i) != sG.end());
        for(auto e : sG[i]) if(e.second == 1) ++ V_new_edges;
    }
    assert(V_new_edges%2 == 0); V_new_edges /= 2;
    cout<<"adding "<<sg_m - m<<" ("<<U_new_edges<<" + "<<V_new_edges<<") more edges."<<endl;
    cout<<"Time cost to materialize the similar graph (without I/O) : "<<integer_to_string(t.elapsed())<<endl;
    if(!write_the_materialized_signed_graph_on_disk) return;
    cout<<"start write."<<endl;
    ofstream fout;
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_m_");
    int int_epsilon = epsilon * 100;
    graph_name.append(to_string(int_epsilon));
    graph_name.append(".txt");
    size_t found = graph_name.find_last_of("/");
    if(found != string::npos) {
        graph_name.insert(found, "/materialized_graphs");
    }
    else {
        cout<<"graph_name cannot find / "<<endl;
        exit(1);
    }
    fout.open(graph_name);
    assert(fout.is_open());
    fout<<sg_n<<"\t"<<sg_m;
    for(auto e : sG){
        for(auto x : e.second){
            fout<<endl<<e.first<<"\t"<<x.first<<"\t"<<x.second;
        }
    }
    fout.close();
    cout<<"finish write."<<endl;
}

void dele_memo()
{
    if(pstart != nullptr){
        delete [] pstart;
        pstart = nullptr;
    }
    if(edges != nullptr){
        delete [] edges;
        edges = nullptr;
    }
    if(degree != nullptr){
        delete [] degree;
        degree = nullptr;
    }
    if(del_ver != nullptr){
        delete [] del_ver;
        del_ver = nullptr;
    }
    if(TMPdeg != nullptr){
        delete [] TMPdeg;
        TMPdeg = nullptr;
    }
    if(Sdeg != nullptr){
        delete [] Sdeg;
        Sdeg = nullptr;
    }
    if(r_pstart != nullptr){
        delete [] r_pstart;
        r_pstart = nullptr;
    }
    if(r_edges != nullptr){
        delete [] r_edges;
        r_edges = nullptr;
    }
    if(r_degree != nullptr){
        delete [] r_degree;
        r_degree = nullptr;
    }
    if(peel_s != nullptr){
        delete [] peel_s;
        peel_s = nullptr;
    }
    if(oid != nullptr){
        delete [] oid;
        oid = nullptr;
    }
    if(nid != nullptr){
        delete [] nid;
        nid = nullptr;
    }
    for(int i = 0; i < max_CS_deg; i++) if(Matrix[i] != nullptr){
        delete [] Matrix[i];
        Matrix[i] = nullptr;
    }
    if(trans != nullptr){
        delete [] trans;
        trans = nullptr;
    }
    if(inPQ != nullptr){
        delete [] inPQ;
        inPQ = nullptr;
    }
    if(deg_inP != nullptr){
        delete [] deg_inP;
        deg_inP = nullptr;
    }
    if(pdeg_inP != nullptr){
        delete [] pdeg_inP;
        pdeg_inP = nullptr;
    }
    if(ndeg_inP != nullptr){
        delete [] ndeg_inP;
        ndeg_inP = nullptr;
    }
    if(os != nullptr){
        delete [] os;
        os = nullptr;
    }
}

void build_index_SA(string graph_name)
{
    ofstream fout;
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_SA.txt");

    fout.open(graph_name);
    assert(fout.is_open());

    //for each vertex, we compute its 2-hop neighbors and store them.
    ui * c = new ui[n];
    memset(c, 0, sizeof(ui)*n);
    for(ui u = 0; u < n; u++){
        vector<ui> two_hop_nei;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(c[w] == 0) {
                    two_hop_nei.push_back(w);
                }
                ++ c[w];
            }
        }
        if(two_hop_nei.empty()) continue;
        sort(two_hop_nei.begin(), two_hop_nei.end(), less<>()); //increasing order
        vector<pair<ui, double>> ordered_2hop_neis;
        for(auto e : two_hop_nei){
            assert(degree[u] >= c[e]);
            assert(degree[e] >= c[e]);
            double simscore = (double) c[e] / (degree[u] + degree[e] - c[e]);
            ordered_2hop_neis.push_back(make_pair(e, simscore));
            c[e] = 0;
        }
        fout<<u<<" ";
        for(auto e : ordered_2hop_neis){
            fout<<e.first<<" "<<setprecision(4)<<e.second<<" ";
        }
        fout<<endl;
    }
    fout.close();
    delete [] c;
    cout<<"*** finish build_index_SA ***"<<endl;
}

void build_index_LG(string graph_name)
{
    cout<<"Start building index LG ( seg_times = "<<seg_times<<", rg_limit = "<<0<<" )."<<endl;
    vector<int> INDEX_vid;
    vector<ui> INDEX_flag;
    vector<vector<Itval>> INDEX_list;
    
    Timer tt;
    long long T_find_2hopneis = 0;
    long long T_cal_sim_and_sort_for_each = 0;
    long long T_build_ranges = 0;
    long long total_2hop_nei_size = 0;
    long long total_shrinked_2hop_nei_size = 0;
    
    long long Phi_total = 0;
    long long Phi_exist = 0;
    
    long long make_idvidual = 0;
    long long make_seg = 0;
    
    //for each vertex, we compute its 2-hop neighbors and store them.
    ui * c = new ui[n];
    memset(c, 0, sizeof(ui)*n);
    for(ui u = 0; u < n; u++){
        tt.restart();
        
        vector<ui> two_hop_nei;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(c[w] == 0) {
                    two_hop_nei.push_back(w);
                }
                ++ c[w];
            }
        }
        
        T_find_2hopneis += tt.elapsed();
        tt.restart();
        
        total_2hop_nei_size += two_hop_nei.size();
        
        vector<pair<ui, double>> ordered_2hop_neis;
        for(auto e : two_hop_nei){
            assert(degree[u] >= c[e]);
            assert(degree[e] >= c[e]);
            double simscore = (double) c[e] / (degree[u] + degree[e] - c[e]);
            if(simscore >= triscr)
                ordered_2hop_neis.push_back(make_pair(e, simscore));
            c[e] = 0;
        }

        total_shrinked_2hop_nei_size += ordered_2hop_neis.size();
        
        ++ Phi_total;
        
        if(ordered_2hop_neis.empty()) continue;
        
        ++ Phi_exist;
        
        sort(ordered_2hop_neis.begin(), ordered_2hop_neis.end(), less<>()); //increasing order
        
        T_cal_sim_and_sort_for_each += tt.elapsed();
        tt.restart();
        
        //start to generate summary ranges
        assert(ordered_2hop_neis.size() >= 1);
        thre_seg = 2;
        INDEX_vid.push_back(u);
        
        if(ordered_2hop_neis.size() < thre_seg){
            vector<Itval> tmpV;
            for(auto e : ordered_2hop_neis){
                tmpV.push_back(Itval(e.first, e.first, e.second, e.second, 1));
            }
            INDEX_list.push_back(tmpV);
            INDEX_flag.push_back(1);
            ++ make_idvidual;
        }
        else{
            vector<Itval> tmpV;
            assert(ordered_2hop_neis.size() >= thre_seg);
            assert(ordered_2hop_neis.size() >= 2);
            int gc = (int)log2(ordered_2hop_neis.size());
            assert(gc < ordered_2hop_neis.size());
            
            gc = gc * (seg_times/2);
            if(gc < 1) gc = 1;
            if(gc > ordered_2hop_neis.size()-1) gc = ordered_2hop_neis.size()-1;
            priority_queue<pair<int, ui>, vector<pair<int, ui>>, less<pair<int, ui>>> kset;
            for(ui i = 0; i < gc; i++){
                int gap_value = ordered_2hop_neis[i+1].first - ordered_2hop_neis[i].first;
                assert(gap_value >= 1);
                kset.push(make_pair(gap_value, i));
            }
            for(ui i = gc; i < ordered_2hop_neis.size() - 1; i++){
                int gap_value = ordered_2hop_neis[i+1].first - ordered_2hop_neis[i].first;
                assert(gap_value >= 1);
                if(gap_value > kset.top().first){
                    kset.pop();
                    kset.push(make_pair(gap_value, i));
                }
            }
            assert(kset.size() == gc);
            vector<ui> positions; //store all index in the ordered_2hop_neis.
            while (!kset.empty()) {
                ui idx = kset.top().second;
                positions.push_back(idx);
                kset.pop();
            }
            assert(positions.size() == gc);
            sort(positions.begin(), positions.end(), less<>()); //increasing order
                        
            double minS, maxS;
            ui start_idx = 0;
            for(ui i = 0; i < positions.size(); i++){
                ui end_idx = positions[i];
                assert(end_idx >= start_idx);
                                
                minS = INF;
                maxS = 0;
                for(ui j = start_idx; j <= end_idx; j++){
                    if(ordered_2hop_neis[j].second < minS) minS = ordered_2hop_neis[j].second;
                    if(ordered_2hop_neis[j].second > maxS) maxS = ordered_2hop_neis[j].second;
                }

                tmpV.push_back(Itval(ordered_2hop_neis[start_idx].first, ordered_2hop_neis[end_idx].first, minS, maxS, end_idx-start_idx+1));
                start_idx = end_idx + 1;
            }
            minS = INF;
            maxS = 0;
            for(ui i = start_idx; i < ordered_2hop_neis.size(); i ++){
                if(ordered_2hop_neis[i].second < minS) minS = ordered_2hop_neis[i].second;
                if(ordered_2hop_neis[i].second > maxS) maxS = ordered_2hop_neis[i].second;
            }

            tmpV.push_back(Itval(ordered_2hop_neis[start_idx].first, ordered_2hop_neis[ordered_2hop_neis.size()-1].first, minS, maxS, ordered_2hop_neis.size() - start_idx));
            INDEX_list.push_back(tmpV);
            INDEX_flag.push_back(2);
            ++ make_seg;
        }

        T_build_ranges += tt.elapsed();
        tt.restart();
        
    }

    long long totalT = T_find_2hopneis + T_cal_sim_and_sort_for_each + T_build_ranges;
    
    cout<<"     -- Time on find 2-hop neis = "<<integer_to_string(T_find_2hopneis)<<" ( "<<((double)T_find_2hopneis/(totalT) )*100<<"% )"<<endl;
    cout<<"     -- Time on cal sim and sort = "<<integer_to_string(T_cal_sim_and_sort_for_each)<<" ( "<<((double)T_cal_sim_and_sort_for_each/(totalT) )*100<<"% )"<<endl;
    cout<<"     -- Time on build segments = "<<integer_to_string(T_build_ranges)<<" ( "<<((double)T_build_ranges/(totalT) )*100<<"% )"<<endl;
    
    assert(Phi_total == n);
        
    assert(INDEX_vid.size() == INDEX_list.size());
    assert(INDEX_vid.size() == INDEX_flag.size());
    
    //binary version
    
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_" + to_string((int)(seg_times*100)));
    graph_name.append("_" + to_string(0));
    graph_name.append("_LG.bin");
    
    FILE * f = Utility::open_file(graph_name.c_str(), "wb");
    
    for(ui i = 0; i < INDEX_vid.size(); i++) {
        fwrite(&INDEX_vid[i], sizeof(int), 1, f);
        int num = (int)INDEX_list[i].size();
        assert(num > 0);
        if(INDEX_flag[i] == 1) {
            num = -num;
            fwrite(&num, sizeof(int), 1, f);
            for(auto &e : INDEX_list[i]) {
                assert(e.e_idx == e.s_idx);
                assert(e.min_score == e.max_score);
                assert(e.c == 1);
                fwrite(&e.s_idx, sizeof(ui), 1, f);
                char aaa = (char)(e.min_score*100);
                fwrite(&aaa, sizeof(char), 1, f);
            }
        }
        else {
            fwrite(&num, sizeof(int), 1, f);
            
            for(auto &e : INDEX_list[i]) {
                char a = (char)(e.min_score*100);
                char b = (char)(e.max_score*100);
                fwrite(&e.s_idx, sizeof(ui), 1, f);
                fwrite(&e.e_idx, sizeof(ui), 1, f);
                fwrite(&a, sizeof(char), 1, f);
                fwrite(&b, sizeof(char), 1, f);
                fwrite(&e.c, sizeof(int), 1, f);
            }
        }
    }
    fclose(f);
    delete [] c;
}

void build_index_rglmt(string graph_name)
{
    cout<<"start to build index rglmt."<<endl;
    
    vector<ui> INDEX_vid;
    vector<vector<Itval>> INDEX_list;
    double RGLMT = 0.1;
    
    Timer tt;
    long long T_find_2hopneis = 0;
    long long T_cal_sim_and_sort_for_each = 0;
    long long T_build_ranges = 0;
    
    ui * c = new ui[n];
    memset(c, 0, sizeof(ui)*n);
    for(ui u = 0; u < n; u++){
        tt.restart();
        vector<ui> two_hop_nei;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(c[w] == 0) {
                    two_hop_nei.push_back(w);
                }
                ++ c[w];
            }
        }
        if(two_hop_nei.empty()) continue;
        
        T_find_2hopneis += tt.elapsed();
        tt.restart();
        
        sort(two_hop_nei.begin(), two_hop_nei.end(), less<>()); //increasing order

        vector<pair<ui, double>> ordered_2hop_neis;
        for(auto e : two_hop_nei){
            assert(degree[u] >= c[e]);
            assert(degree[e] >= c[e]);
            double simscore = (double) c[e] / (degree[u] + degree[e] - c[e]);
            ordered_2hop_neis.push_back(make_pair(e, simscore));
            c[e] = 0;
        }
        T_cal_sim_and_sort_for_each += tt.elapsed();
        tt.restart();
        
        //start to generate summary ranges
        assert(ordered_2hop_neis.size() >= 1);
        INDEX_vid.push_back(u);
        
        vector<Itval> tmpVec;
        
        ui t_size = ordered_2hop_neis.size();
        ui s_idx = 0;
        while (1) {
            pair<ui, double> & s_node = ordered_2hop_neis[s_idx];
            double std_scr = s_node.second;
            Itval tmpI(s_node.first, s_node.first, s_node.second, s_node.second, 1);
            ++ s_idx;
            while (s_idx < t_size) {
                pair<ui, double> & nt_node = ordered_2hop_neis[s_idx];
                assert(nt_node.first > tmpI.e_idx);
                if(fabs(nt_node.second - std_scr) > RGLMT) break;
                tmpI.e_idx = nt_node.first;
                if(nt_node.second < tmpI.min_score) tmpI.min_score = nt_node.second;
                if(nt_node.second > tmpI.max_score) tmpI.max_score = nt_node.second;
                ++ tmpI.c;
                ++ s_idx;
            }
            tmpVec.push_back(tmpI);
            if(s_idx >= t_size) break;
        }
        INDEX_list.push_back(tmpVec);
        
        T_build_ranges += tt.elapsed();
        tt.restart();
    }

    long long totalT = T_find_2hopneis + T_cal_sim_and_sort_for_each + T_build_ranges;
    
    cout<<"     T_find_2hopneis = "<<integer_to_string(T_find_2hopneis)<<" ( "<<((double)T_find_2hopneis/(totalT) )*100<<"% )"<<endl;
    cout<<"     T_cal_sim_and_sort_for_each = "<<integer_to_string(T_cal_sim_and_sort_for_each)<<" ( "<<((double)T_cal_sim_and_sort_for_each/(totalT) )*100<<"% )"<<endl;
    cout<<"     T_build_ranges = "<<integer_to_string(T_build_ranges)<<" ( "<<((double)T_build_ranges/(totalT) )*100<<"% )"<<endl;

    
    assert(INDEX_vid.size() == INDEX_list.size());
    
    ofstream fout;
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_rglmt.txt");

    //cout<<graph_name<<endl;
    fout.open(graph_name);
    assert(fout.is_open());
    
    for(ui i = 0; i < INDEX_vid.size(); i++){
        fout<<INDEX_vid[i]<<" ";
        for(auto e : INDEX_list[i]){
            fout<<e.s_idx<<" "<<e.e_idx<<" "<<setprecision(4)<<e.min_score<<" "<<setprecision(4)<<e.max_score<<" "<<e.c<<" ";
        }
        fout<<endl;
    }
    fout.close();
    
    delete [] c;
    cout<<"finish building index rglmt."<<endl;
}

void build_index_GR(string graph_name)
{
    cout<<"*** build_index_GR ***"<<endl;
    
    vector<ui> INDEX_vid;
    vector<vector<Itval>> INDEX_list;
    double RGLMT = 0.1;
    
    Timer tt;
    long long T_find_2hopneis = 0;
    long long T_cal_sim_and_sort_for_each = 0;
    long long T_build_ranges = 0;
    
    ui * c = new ui[n];
    memset(c, 0, sizeof(ui)*n);
    
    Timer t1;
    long long t1_allranges = 0, t1_selectrange = 0, t1_partition = 0;
    
    for(ui u = 0; u < n; u++){
        tt.restart();
        vector<ui> two_hop_nei;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(c[w] == 0) {
                    two_hop_nei.push_back(w);
                }
                ++ c[w];
            }
        }
        if(two_hop_nei.empty()) continue;
        
        T_find_2hopneis += tt.elapsed();
        tt.restart();
        
        sort(two_hop_nei.begin(), two_hop_nei.end(), less<>()); //increasing order

        vector<pair<ui, double>> ordered_2hop_neis;
        for(auto e : two_hop_nei){
            assert(degree[u] >= c[e]);
            assert(degree[e] >= c[e]);
            double simscore = (double) c[e] / (degree[u] + degree[e] - c[e]);
            ordered_2hop_neis.push_back(make_pair(e, simscore));
            c[e] = 0;
        }
        
        T_cal_sim_and_sort_for_each += tt.elapsed();
        tt.restart();

        ui t_size = ordered_2hop_neis.size();
        INDEX_vid.push_back(u);

        if(t_size < 10){
            vector<Itval> tmpV;
            for(auto e : ordered_2hop_neis){
                tmpV.push_back(Itval(e.first, e.first, e.second, e.second, 1));
            }
            INDEX_list.push_back(tmpV);
        }
        else{
            assert(t_size >= 10);
            
            t1.restart();
            
            vector<Range> allRanges;
            for(ui i = 0; i < t_size; i++){
                ui v = ordered_2hop_neis[i].first;
                double std_scr = ordered_2hop_neis[i].second;
                
                Range tmpR(v, 1, i, i); //Range
                
                if(i > 0){ // need go left
                    ui j = i - 1;
                    assert(j >= 0);
                    while (1) {
                        ui tmpID = ordered_2hop_neis[j].first;
                        double tmpScore = ordered_2hop_neis[j].second;
                        if(fabs(tmpScore - std_scr) > RGLMT) {break;}
                        assert(v > tmpID);
                        tmpR.Lidx = j;
                        if(j == 0) break;
                        -- j;
                    }
                }
                //go right
                ui j = i + 1;
                while (j < t_size) {
                    ui tmpID = ordered_2hop_neis[j].first;
                    double tmpScore = ordered_2hop_neis[j].second;
                    if(fabs(tmpScore - std_scr) > RGLMT) {break;}
                    assert(v < tmpID);
                    tmpR.Ridx = j;
                    ++ j;
                }
                
                assert(tmpR.Ridx >= tmpR.Lidx);
                assert(tmpR.Lidx >= 0 && tmpR.Lidx < t_size);
                assert(tmpR.Ridx >= 0 && tmpR.Ridx < t_size);
                
                tmpR.rgC = ordered_2hop_neis[tmpR.Ridx].first - ordered_2hop_neis[tmpR.Lidx].first + 1;
                allRanges.push_back(tmpR);
            }
            assert(allRanges.size() == t_size);
            
            t1_allranges += t1.elapsed();
            t1.restart();

            int rc = (int)log2(ordered_2hop_neis.size()) - 1;
            assert(rc < t_size);
            
            vector<pair<ui, ui>> rgS;
            for(ui i = 0; i < rc; i++){
                int max_rgC = -1;
                ui idx = 0;
                for(ui j = 0; j < allRanges.size(); j++){
                    if(allRanges[j].rgC > max_rgC && c[j] == 0) {
                        max_rgC = allRanges[j].rgC;
                        idx = j;
                    }
                }
                if(max_rgC == -1) break;
                
                Range selected_range = allRanges[idx];
                rgS.push_back(make_pair(selected_range.Lidx, selected_range.Ridx));
                
                for(ui j = selected_range.Lidx; j <= selected_range.Ridx; j++) c[j] = 1;
                
                //update the influenced ranges
                for(ui j = 0; j < allRanges.size(); j++){
                    if(c[j]) continue;
                    Range & influenced_range = allRanges[j];
                    if(j < selected_range.Lidx) {
                        assert(influenced_range.Lidx < selected_range.Lidx);
                        assert(influenced_range.Ridx < selected_range.Ridx);
                        if(influenced_range.Ridx >= selected_range.Lidx){
                            assert(selected_range.Lidx >= 1);
                            influenced_range.Ridx = selected_range.Lidx - 1;
                            assert(influenced_range.Lidx <= influenced_range.Ridx);
                            influenced_range.rgC = ordered_2hop_neis[influenced_range.Ridx].first - ordered_2hop_neis[influenced_range.Lidx].first + 1;
                            assert(influenced_range.rgC >= 1);
                        }
                    }
                    else{
                        assert(j > selected_range.Ridx);
                        assert(influenced_range.Ridx > selected_range.Ridx);
                        assert(influenced_range.Lidx > selected_range.Lidx);
                        if(influenced_range.Lidx <= selected_range.Ridx){
                            influenced_range.Lidx = selected_range.Ridx + 1;
                            assert(influenced_range.Lidx < allRanges.size());
                            assert(influenced_range.Lidx <= influenced_range.Ridx);
                            influenced_range.rgC = ordered_2hop_neis[influenced_range.Ridx].first - ordered_2hop_neis[influenced_range.Lidx].first + 1;
                            assert(influenced_range.rgC >= 1);
                        }
                    }
                }
                
            }
            assert(rgS.size() <= rc);
            for(auto & e : rgS){
                for(ui i = e.first; i <= e.second; i++) c[i] = 0;
            }

            t1_selectrange += t1.elapsed();
            t1.restart();
            
            sort(rgS.begin(), rgS.end(), less<>());
            
            vector<Itval> tmpVec;
            
            ui start_idx = 0;
            for(ui i = 0; i < rgS.size(); i++){
                ui rg_s = rgS[i].first;
                ui rg_e = rgS[i].second;
                assert(rg_e >= rg_s);
                if(rg_s > start_idx){
                    pair<ui, double> & tmpv = ordered_2hop_neis[start_idx];
                    Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                    ++ start_idx;
                    while (start_idx < rg_s) {
                        pair<ui, double> & nxtv = ordered_2hop_neis[start_idx];
                        tmpI.e_idx = nxtv.first;
                        if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                        if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                        ++ tmpI.c;
                        ++ start_idx;
                    }
                    assert(start_idx == rg_s);
                    tmpVec.push_back(tmpI);
                }
                //construct rg_s -> rg_e
                pair<ui, double> & tmpv = ordered_2hop_neis[rg_s];
                Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                for(ui j = rg_s + 1; j <= rg_e; j++){
                    pair<ui, double> & nxtv = ordered_2hop_neis[j];
                    tmpI.e_idx = nxtv.first;
                    if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                    if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                    ++ tmpI.c;
                }
                tmpVec.push_back(tmpI);
                start_idx = rg_e + 1;
            }
            
            if(start_idx < t_size){
                pair<ui, double> & tmpv = ordered_2hop_neis[start_idx];
                Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                ++ start_idx;
                while (start_idx < t_size) {
                    pair<ui, double> & nxtv = ordered_2hop_neis[start_idx];
                    tmpI.e_idx = nxtv.first;
                    if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                    if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                    ++ tmpI.c;
                    ++ start_idx;
                }
                tmpVec.push_back(tmpI);
            }
            INDEX_list.push_back(tmpVec);
            
            t1_partition += t1.elapsed();
            t1.restart();
            
        }
        
        T_build_ranges += tt.elapsed();
        tt.restart();
    }

    long long totalT = T_find_2hopneis + T_cal_sim_and_sort_for_each + T_build_ranges;
    
    cout<<"     T_find_2hopneis = "<<integer_to_string(T_find_2hopneis)<<" ( "<<((double)T_find_2hopneis/(totalT) )*100<<"% )"<<endl;
    cout<<"     T_cal_sim_and_sort_for_each = "<<integer_to_string(T_cal_sim_and_sort_for_each)<<" ( "<<((double)T_cal_sim_and_sort_for_each/(totalT) )*100<<"% )"<<endl;
    cout<<"     T_build_ranges = "<<integer_to_string(T_build_ranges)<<" ( "<<((double)T_build_ranges/(totalT) )*100<<"% )"<<endl;
        
    assert(INDEX_vid.size() == INDEX_list.size());
    
    ofstream fout;
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_GR.txt");

    fout.open(graph_name);
    assert(fout.is_open());
    
    for(ui i = 0; i < INDEX_vid.size(); i++){
        fout<<INDEX_vid[i]<<" ";
        for(auto e : INDEX_list[i]){
            fout<<e.s_idx<<" "<<e.e_idx<<" "<<setprecision(4)<<e.min_score<<" "<<setprecision(4)<<e.max_score<<" "<<e.c<<" ";
        }
        fout<<endl;
    }
    fout.close();
    
    delete [] c;
    cout<<"*** finish build_index_GR ***"<<endl;

}

void build_index_GRL(string graph_name)
{
    cout<<"*** build_index_GRL ***"<<endl;
    
    vector<ui> INDEX_vid;
    vector<vector<Itval>> INDEX_list;
    double RGLMT = 0.1;
    
    Timer tt;
    long long T_find_2hopneis = 0;
    long long T_cal_sim_and_sort_for_each = 0;
    long long T_build_ranges = 0;
    long long total_2hop_nei_size = 0;
    long long total_shrinked_2hop_nei_size = 0;
    
    ui * c = new ui[n];
    memset(c, 0, sizeof(ui)*n);
    
    Timer t1;
    long long t1_allranges = 0, t1_selectrange = 0, t1_partition = 0;
    
    for(ui u = 0; u < n; u++){
        tt.restart();
        
        vector<ui> two_hop_nei;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(c[w] == 0) {
                    two_hop_nei.push_back(w);
                }
                ++ c[w];
            }
        }
        
        T_find_2hopneis += tt.elapsed();
        tt.restart();
        
        total_2hop_nei_size += two_hop_nei.size();

        vector<pair<ui, double>> ordered_2hop_neis;
        for(auto e : two_hop_nei){
            assert(degree[u] >= c[e]);
            assert(degree[e] >= c[e]);
            double simscore = (double) c[e] / (degree[u] + degree[e] - c[e]);
            if(simscore >= 0.05)
                ordered_2hop_neis.push_back(make_pair(e, simscore));
            c[e] = 0;
        }
        
        total_shrinked_2hop_nei_size += ordered_2hop_neis.size();
        
        if(ordered_2hop_neis.empty()) continue;
        sort(ordered_2hop_neis.begin(), ordered_2hop_neis.end(), less<>()); //increasing order
        
        T_cal_sim_and_sort_for_each += tt.elapsed();
        tt.restart();

        ui t_size = ordered_2hop_neis.size();
        INDEX_vid.push_back(u);
        
        if(t_size < 8){
            vector<Itval> tmpV;
            for(auto e : ordered_2hop_neis){
                tmpV.push_back(Itval(e.first, e.first, e.second, e.second, 1));
            }
            INDEX_list.push_back(tmpV);
        }
        else{
            assert(t_size >= 8);
            t1.restart();
            vector<pair<int, pair<ui, ui>>> allRanges;
            ui s_idx = 0;
            
            while (1) { //linear scan of ordered_2hop_neis, at the same time, record all ranges
                pair<int, pair<ui, ui>> tmpR (1, make_pair(s_idx, s_idx));
                double std_scr = ordered_2hop_neis[s_idx].second;
                
                ++ s_idx;
                while (s_idx < t_size) {
                    pair<ui, double> & tnode = ordered_2hop_neis[s_idx];
                    if(fabs(tnode.second - std_scr) > RGLMT) break;
                    assert(tmpR.second.second < s_idx);
                    tmpR.second.second = s_idx;
                    ++ s_idx;
                }
                tmpR.first = ordered_2hop_neis[tmpR.second.second].first - ordered_2hop_neis[tmpR.second.first].first + 1;
                allRanges.push_back(tmpR);
                if(s_idx >= t_size) break;
            }
            
            t1_allranges += t1.elapsed();
            t1.restart();

            int rc = (int)log2(ordered_2hop_neis.size()) - 1;
            assert(rc < t_size);
            
            vector<pair<ui, ui>> rgS;
            if(rc >= allRanges.size()){
                for(auto & e : allRanges)
                    rgS.push_back(e.second);
            }
            else{
                assert(rc < allRanges.size());
                priority_queue<pair<int, pair<ui, ui>>, vector<pair<int, pair<ui, ui>>>, greater<pair<int, pair<ui, ui>>>> kset;
                for(ui i = 0 ; i < rc; i++){
                    pair<int, pair<ui, ui>> & x = allRanges[i];
                    kset.push(x);
                }
                for(ui i = rc; i < allRanges.size(); i++){
                    pair<int, pair<ui, ui>> & x = allRanges[i];
                    if(x.first > kset.top().first){
                        kset.pop();
                        kset.push(x);
                    }
                }
                assert(kset.size() == rc);
                while (!kset.empty()) {
                    rgS.push_back(kset.top().second);
                    kset.pop();
                }
            }
            
            t1_selectrange += t1.elapsed();
            t1.restart();
            
            sort(rgS.begin(), rgS.end(), less<>());
            
            vector<Itval> tmpVec;
            
            ui start_idx = 0;
            for(ui i = 0; i < rgS.size(); i++){
                ui rg_s = rgS[i].first;
                ui rg_e = rgS[i].second;
                assert(rg_e >= rg_s);
                if(start_idx < rg_s){
                    pair<ui, double> & tmpv = ordered_2hop_neis[start_idx];
                    Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                    ++ start_idx;
                    while (start_idx < rg_s) {
                        pair<ui, double> & nxtv = ordered_2hop_neis[start_idx];
                        tmpI.e_idx = nxtv.first;
                        if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                        if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                        ++ tmpI.c;
                        ++ start_idx;
                    }
                    assert(start_idx == rg_s);
                    tmpVec.push_back(tmpI);
                }
                pair<ui, double> & tmpv = ordered_2hop_neis[rg_s];
                Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                for(ui j = rg_s + 1; j <= rg_e; j++){
                    pair<ui, double> & nxtv = ordered_2hop_neis[j];
                    tmpI.e_idx = nxtv.first;
                    if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                    if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                    ++ tmpI.c;
                }
                tmpVec.push_back(tmpI);
                start_idx = rg_e + 1;
            }
            
            if(start_idx < t_size){
                pair<ui, double> & tmpv = ordered_2hop_neis[start_idx];
                Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                ++ start_idx;
                while (start_idx < t_size) {
                    pair<ui, double> & nxtv = ordered_2hop_neis[start_idx];
                    tmpI.e_idx = nxtv.first;
                    if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                    if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                    ++ tmpI.c;
                    ++ start_idx;
                }
                tmpVec.push_back(tmpI);
            }
            INDEX_list.push_back(tmpVec);
            
            t1_partition += t1.elapsed();
            t1.restart();
            
        }
        
        T_build_ranges += tt.elapsed();
        tt.restart();
    }

    long long totalT = T_find_2hopneis + T_cal_sim_and_sort_for_each + T_build_ranges;
    
    cout<<"     T_find_2hopneis = "<<integer_to_string(T_find_2hopneis)<<" ( "<<((double)T_find_2hopneis/(totalT) )*100<<"% )"<<endl;
    cout<<"     T_cal_sim_and_sort_for_each = "<<integer_to_string(T_cal_sim_and_sort_for_each)<<" ( "<<((double)T_cal_sim_and_sort_for_each/(totalT) )*100<<"% )"<<endl;
    cout<<"     T_build_ranges = "<<integer_to_string(T_build_ranges)<<" ( "<<((double)T_build_ranges/(totalT) )*100<<"% )"<<endl;

    assert(INDEX_vid.size() == INDEX_list.size());
    
    ofstream fout;
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_GRL.txt");

    //cout<<graph_name<<endl;
    fout.open(graph_name);
    assert(fout.is_open());
    
    for(ui i = 0; i < INDEX_vid.size(); i++){
        fout<<INDEX_vid[i]<<" ";
        for(auto e : INDEX_list[i]){
            fout<<e.s_idx<<" "<<e.e_idx<<" "<<setprecision(3)<<e.min_score<<" "<<setprecision(3)<<e.max_score<<" "<<e.c<<" ";
        }
        fout<<endl;
    }
    fout.close();
    
    delete [] c;
    cout<<"*** finish build_index_GRL ***"<<endl;

}

node * build_tree(vector<node *> A)
{
    if(A.size() == 1) {
        return A[0];
    }

    int pos = 0;
    vector<node *> tmp_A;
    while (pos < A.size()) {
        if(pos + 1 < A.size()) {
            node * p1 = A[pos];
            node * p2 = A[pos+1];
            node * new_node = new node;
            new_node->L = p1;
            new_node->R = p2;
            new_node->P = nullptr;
            p1->P = new_node;
            p2->P = new_node;
            new_node->mins = min(p1->mins, p2->mins);
            new_node->maxs = max(p1->maxs, p2->maxs);
            new_node->mark = 0;
            new_node->isv = 0;
            new_node->idx = 0;
            tmp_A.push_back(new_node);
            pos += 2;
        }
        else { //remain 1 node in A
            node * p1 = A[pos];
            node * new_node = new node;
            new_node->L = p1;
            new_node->R = nullptr;
            new_node->P = nullptr;
            p1->P = new_node;
            new_node->mins = p1->mins;
            new_node->maxs = p1->maxs;
            new_node->mark = 0;
            new_node->isv = 0;
            new_node->idx = 0;
            tmp_A.push_back(new_node);
            pos += 2;
        }
    }
    return build_tree(tmp_A);
}

void show_tree(node * ptr)
{
    if(ptr == nullptr) return;
    if(ptr->isv == 1) {
        cout<<"v"<<ptr->idx<<endl;
    }
    else {
        cout<<"["<<ptr->mins<<","<<ptr->maxs<<"]"<<endl;
    }
    show_tree(ptr->L);
    show_tree(ptr->R);
}

void cp_tree(vector<pair<ui, double>> & ordered_2hop_neis, vector<node *> & A, node * root, vector<pair<ui, int>> & B, double RGLMT)
{
    assert(ordered_2hop_neis.size() == A.size());
    assert(A.size() == B.size());
    assert(root != nullptr);
    for(ui i = 0; i < ordered_2hop_neis.size(); i++) {
        vector<node *> broadway;
        node * unode = A[i];
        assert(unode != nullptr);
        while (unode != nullptr) { //find broadway: from unode -> root
            unode->mark = 1;
            broadway.push_back(unode);
            unode = unode->P;
        }
        assert(unode == nullptr);
        double minscore = ordered_2hop_neis[i].second;
        double maxscore = ordered_2hop_neis[i].second;
        
        //broadway is just the path from unode to root, which must exist
        //next, we need to find the downnode pointed by a rightward branch along this broadway, which may be empty
        //downnode is the first node that makes current segment non-steady!
        node * downnode = nullptr;
        for(auto e : broadway) {
            if(e->R != nullptr && e->R->mark == 0) { //e->R may be a downnode
                downnode = e->R;
                double tmp_minscore = min(downnode->mins, minscore);
                double tmp_maxscore = max(downnode->maxs, maxscore);
                if((tmp_maxscore - tmp_minscore) > RGLMT) {
                    break; //find the node where we should go down
                }
                minscore = min(downnode->mins, minscore);
                maxscore = max(downnode->maxs, maxscore);
            }
        }
        
        if(downnode == nullptr) { // we are processing the last vertex
            assert(i == (ordered_2hop_neis.size() - 1) );
            B[i].first = i;
            B[i].second = 1;
            for (auto e : broadway) e->mark = 0;
        }
        else {
            assert(downnode != nullptr);
            assert((maxscore - minscore) <= RGLMT);
            while (downnode != nullptr) { //go down till a vertex (leaf)
                if(downnode->isv == 1) { //reach a vertex (leaf)
                    ui e_idx = downnode->idx;
                    assert(e_idx > 0);
                    double tmp_minscore = min(downnode->mins, minscore);
                    double tmp_maxscore = max(downnode->maxs, maxscore);
                    if((tmp_maxscore - tmp_minscore) > RGLMT) {
                        -- e_idx;
                    }
                    assert(e_idx >= i);
                    B[i].first = e_idx;
                    B[i].second = ordered_2hop_neis[e_idx].first - ordered_2hop_neis[i].first + 1;
                    for (auto e : broadway) e->mark = 0;
                    break;
                }
                //firstly check its left child
                assert(downnode->L != nullptr);
                double tmp_minscore = min(downnode->L->mins, minscore);
                double tmp_maxscore = max(downnode->L->maxs, maxscore);
                if((tmp_maxscore - tmp_minscore) > RGLMT) {
                    downnode = downnode->L;
                    continue;
                }
                minscore = min(downnode->L->mins, minscore);
                maxscore = max(downnode->L->maxs, maxscore);
                
                //already involve the left branch
                downnode = downnode->R;
                if(downnode == nullptr) {
                    B[i].first = B.size() - 1;
                    B[i].second = ordered_2hop_neis[B.size() - 1].first - ordered_2hop_neis[i].first + 1;
                    for (auto e : broadway) e->mark = 0;
                    break;
                }
            }
        }
    }
}

void destroy_tree(node * root)
{
    if(root == nullptr) return;
    destroy_tree(root->L);
    destroy_tree(root->R);
    delete root;
    root = nullptr;
}

void build_index_SS(string graph_name)
{
    cout<<"Start building index SS ( seg_times = "<<seg_times<<", rg_limit = "<<rg_limit<<" )."<<endl;
    
    vector<ui> INDEX_vid;
    vector<ui> INDEX_flag;
    vector<vector<Itval>> INDEX_list;
    double RGLMT = rg_limit;
    
    Timer tt;
    long long T_find_2hopneis = 0;
    long long T_cal_sim_and_sort_for_each = 0;
    long long T_build_ranges = 0;
    long long total_2hop_nei_size = 0;
    long long total_shrinked_2hop_nei_size = 0;
    
    Timer ttt;
    long long T_build_tree = 0;
    long long T_cp_tree = 0;
    long long T_destroy_tree = 0;
    long long T_range_tree = 0;
    
    long long Phi_total = 0;
    long long Phi_exist = 0;
    
    long long make_idvidual = 0;
    long long make_seg = 0;
    
    //for each vertex, we compute its 2-hop neighbors and store them.
    ui * c = new ui[n];
    memset(c, 0, sizeof(ui)*n);
    
    Timer t1;
    long long t1_allranges = 0, t1_selectrange = 0, t1_partition = 0;
    
    for(ui u = 0; u < n; u++){
        tt.restart();
        
        vector<ui> two_hop_nei;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            for(ui j = pstart[v]; j < pstart[v+1]; j++){
                ui w = edges[j];
                if(w == u) continue;
                if(c[w] == 0) {
                    two_hop_nei.push_back(w);
                }
                ++ c[w];
            }
        }
        
        T_find_2hopneis += tt.elapsed();
        tt.restart();
        
        total_2hop_nei_size += two_hop_nei.size();

        vector<pair<ui, double>> ordered_2hop_neis;
        for(auto e : two_hop_nei){
            assert(degree[u] >= c[e]);
            assert(degree[e] >= c[e]);
            double simscore = (double) c[e] / (degree[u] + degree[e] - c[e]);
            if(simscore >= triscr)
                ordered_2hop_neis.push_back(make_pair(e, simscore));
            c[e] = 0;
        }
        
        total_shrinked_2hop_nei_size += ordered_2hop_neis.size();
        
        ++ Phi_total;
        
        if(ordered_2hop_neis.empty()) continue;
        
        ++ Phi_exist;
        
        sort(ordered_2hop_neis.begin(), ordered_2hop_neis.end(), less<>()); //increasing order
        
        T_cal_sim_and_sort_for_each += tt.elapsed();
        tt.restart();
        
        INDEX_vid.push_back(u);
        
        if(ordered_2hop_neis.size() < thre_seg){
            vector<Itval> tmpV;
            for(auto e : ordered_2hop_neis){
                tmpV.push_back(Itval(e.first, e.first, e.second, e.second, 1));
            }
            INDEX_list.push_back(tmpV);
            INDEX_flag.push_back(1);
            ++ make_idvidual;
        }
        else{
            
            assert(ordered_2hop_neis.size() >= thre_seg);
            assert(ordered_2hop_neis.size() >= 2);
            
            //tree manner
            ttt.restart();
            vector<node *> A;
            vector<pair<ui, int>> B;
            int cnt = 0;
            for(auto &e : ordered_2hop_neis) {
                node * tmp_node = new node;
                tmp_node->L = nullptr;
                tmp_node->R = nullptr;
                tmp_node->P = nullptr;
                tmp_node->mins = e.second;
                tmp_node->maxs = e.second;
                tmp_node->mark = 0;
                tmp_node->isv = 1;
                tmp_node->idx = cnt++;
                A.push_back(tmp_node);
                B.push_back(make_pair(0, 0)); //unexpanded
            }
            assert(cnt == A.size());
            node * tree_root = build_tree(A);
            
            T_build_tree += ttt.elapsed();
            ttt.restart();
            
            //find temporal steady segment for each vertex in ordered_2hop_neis
            cp_tree(ordered_2hop_neis, A, tree_root, B, RGLMT);
            
            T_cp_tree += ttt.elapsed();
            ttt.restart();
                
            destroy_tree(tree_root);
            
            T_destroy_tree += ttt.elapsed();
            ttt.restart();
            //tree manner
            
            //B has been set
            //select rc segments
            int rc = (int)log2(ordered_2hop_neis.size());
            assert(rc < ordered_2hop_neis.size());
            
            rc = rc * (seg_times/2);
            
            if(rc > ordered_2hop_neis.size()) rc = ordered_2hop_neis.size();
            
            assert(rc >= 0 && rc <= ordered_2hop_neis.size());
            
            assert(A.size() == B.size());
            assert(A.size() == ordered_2hop_neis.size());
                        
            vector<int> C;
            C.resize(ordered_2hop_neis.size(), 1);
            
            vector<pair<ui, ui>> rgS;
            if(rc == 0) {
                rgS.push_back(make_pair(0, ordered_2hop_neis.size() - 1));
            }
            while (rc > 0) {
                pair<ui, ui> tmp_p;
                int cur_largest_cp = 0;
                for(ui i = 0; i < ordered_2hop_neis.size(); i ++) if(C[i] == 1) {
                    if(B[i].second > cur_largest_cp) {
                        tmp_p = make_pair(i, B[i].first);
                        cur_largest_cp = B[i].second;
                    }
                }
                if(cur_largest_cp == 0) break;
                rgS.push_back(tmp_p);
                assert(tmp_p.first <= tmp_p.second);
                for(ui j = tmp_p.first; j <= tmp_p.second; j++) C[j] = 0;
                
                //update other influenced temporal segments
                for(ui i = 0; i < tmp_p.first; i++) if (C[i] == 1) {
                    if(B[i].first >= tmp_p.first) {
                        assert(tmp_p.first >= 1);
                        ui new_eidx = tmp_p.first - 1;
                        int new_cp = ordered_2hop_neis[new_eidx].first - ordered_2hop_neis[i].first + 1;
                        B[i].first = new_eidx;
                        B[i].second = new_cp;
                    }
                }
                -- rc;
            }
            
            sort(rgS.begin(), rgS.end(), less<>());
            
            //check rgS
            assert(rgS.size() > 0);
            
            vector<Itval> tmpVec;
            
            ui start_idx = 0;
            for(ui i = 0; i < rgS.size(); i++){
                ui rg_s = rgS[i].first;
                ui rg_e = rgS[i].second;
                assert(rg_e >= rg_s);
                if(start_idx < rg_s){
                    pair<ui, double> & tmpv = ordered_2hop_neis[start_idx];
                    Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                    ++ start_idx;
                    while (start_idx < rg_s) {
                        pair<ui, double> & nxtv = ordered_2hop_neis[start_idx];
                        tmpI.e_idx = nxtv.first;
                        if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                        if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                        ++ tmpI.c;
                        ++ start_idx;
                    }
                    assert(start_idx == rg_s);
                    tmpVec.push_back(tmpI);
                }
                //construct rg_s -> rg_e
                pair<ui, double> & tmpv = ordered_2hop_neis[rg_s];
                Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                for(ui j = rg_s + 1; j <= rg_e; j++){
                    pair<ui, double> & nxtv = ordered_2hop_neis[j];
                    tmpI.e_idx = nxtv.first;
                    if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                    if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                    ++ tmpI.c;
                }
                tmpVec.push_back(tmpI);
                start_idx = rg_e + 1;
            }
            
            if(start_idx < ordered_2hop_neis.size()){
                pair<ui, double> & tmpv = ordered_2hop_neis[start_idx];
                Itval tmpI(tmpv.first, tmpv.first, tmpv.second, tmpv.second, 1);
                ++ start_idx;
                while (start_idx < ordered_2hop_neis.size()) {
                    pair<ui, double> & nxtv = ordered_2hop_neis[start_idx];
                    tmpI.e_idx = nxtv.first;
                    if(nxtv.second < tmpI.min_score) tmpI.min_score = nxtv.second;
                    if(nxtv.second > tmpI.max_score) tmpI.max_score = nxtv.second;
                    ++ tmpI.c;
                    ++ start_idx;
                }
                tmpVec.push_back(tmpI);
            }
            
            INDEX_list.push_back(tmpVec);
            INDEX_flag.push_back(2);
            ++ make_seg;
            
            T_range_tree += ttt.elapsed();
            ttt.restart();
            
            t1_partition += t1.elapsed();
            t1.restart();
            
        }
        
        T_build_ranges += tt.elapsed();
        tt.restart();
    }

    long long totalT = T_find_2hopneis + T_cal_sim_and_sort_for_each + T_build_ranges;
    
    cout<<"     -- Time on find 2-hop neis = "<<integer_to_string(T_find_2hopneis)<<" ( "<<((double)T_find_2hopneis/(totalT) )*100<<"% )"<<endl;
    cout<<"     -- Time on cal sim and sort = "<<integer_to_string(T_cal_sim_and_sort_for_each)<<" ( "<<((double)T_cal_sim_and_sort_for_each/(totalT) )*100<<"% )"<<endl;
    cout<<"     -- Time on build segments = "<<integer_to_string(T_build_ranges)<<" ( "<<((double)T_build_ranges/(totalT) )*100<<"% )"<<endl;
    
    assert(INDEX_vid.size() == INDEX_list.size());
    assert(INDEX_vid.size() == INDEX_flag.size());
        
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_" + to_string((int)(seg_times*100)));
    graph_name.append("_" + to_string((int)(rg_limit*100)));
    graph_name.append("_SS.bin");
    
    FILE * f = Utility::open_file(graph_name.c_str(), "wb");
    
    for(ui i = 0; i < INDEX_vid.size(); i++) {
        fwrite(&INDEX_vid[i], sizeof(int), 1, f);
        int num = (int)INDEX_list[i].size();
        assert(num > 0);
        if(INDEX_flag[i] == 1) {
            num = -num;
            fwrite(&num, sizeof(int), 1, f);
            for(auto &e : INDEX_list[i]) {
                assert(e.e_idx == e.s_idx);
                assert(e.min_score == e.max_score);
                assert(e.c == 1);
                fwrite(&e.s_idx, sizeof(ui), 1, f);
                char aaa = (char)(e.min_score*100);
                fwrite(&aaa, sizeof(char), 1, f);
            }
        }
        else {
            fwrite(&num, sizeof(int), 1, f);
            
            for(auto &e : INDEX_list[i]) {
                char a = (char)(e.min_score*100);
                char b = (char)(e.max_score*100);
                fwrite(&e.s_idx, sizeof(ui), 1, f);
                fwrite(&e.e_idx, sizeof(ui), 1, f);
                fwrite(&a, sizeof(char), 1, f);
                fwrite(&b, sizeof(char), 1, f);
                fwrite(&e.c, sizeof(int), 1, f);
            }
        }
    }
    fclose(f);
    delete [] c;
}

void check_results(vector<pair<vector<ui>, vector<ui>>> r)
{
    cout<<"\t = = = = checking results = = = = \t"<<endl;
    if(r.empty()) {
        cout<<"Note that results is empty !!!"<<endl;
    }
    bool flag = true;
    for(auto clique : r) {
        vector<ui> C1 = clique.first;
        vector<ui> C2 = clique.second;

        //check connection links
        for(auto u : C1) {
            unordered_set<ui> C1nei;
            for(ui i = pstart[u]; i < pstart[u+1]; i++ ) {
                ui v = edges[i];
                C1nei.insert(v);
            }
            for(auto x : C2) {
                if(C1nei.find(x) == C1nei.end()) {
                    cout<<"find a vertex in C2 that is not neighbors of u in C1."<<endl;
                    flag = false;
                    exit(1);
                }
            }
        }

        //check similarity
        for(ui i = 0; i < C1.size(); i++) {
            ui u = C1[i];
            for(ui j = i + 1; j < C1.size(); j++) {
                ui v = C1[j];
                if(js(u, v) < epsilon) {
                    cout<<"find dissimilar pairs in C1."<<endl;
                    flag = false;
                    exit(1);
                }
            }
        }

        for(ui i = 0; i < C2.size(); i++) {
            ui u = C2[i];
            for(ui j = i + 1; j < C2.size(); j++) {
                ui v = C2[j];
                if(js(u, v) < epsilon) {
                    cout<<"find dissimilar pairs in C2."<<endl;
                    flag = false;
                    exit(1);
                }
            }
        }

    }
    if(flag == false) cout<<"\t = = = = = = WRONG! = = = = = = \t"<<endl;
    else cout<<"\t = = = = = = CORRECT! = = = = = = \t"<<endl;
    
    double total_js = 0;
    double total_C1_js = 0;
    double total_C2_js = 0;
    for(auto clique : r) {
        vector<ui> C1 = clique.first;
        vector<ui> C2 = clique.second;
        
        double tmp_C1_js = 0;
        for(ui i = 0; i < C1.size(); i++) {
            ui u = C1[i];
            for(ui j = i + 1; j < C1.size(); j++) {
                ui v = C1[j];
                tmp_C1_js += js(u, v);
            }
        }
        double allpair = (C1.size()*(C1.size()-1))/2;
        tmp_C1_js = (double) tmp_C1_js / allpair;
        
        total_C1_js += tmp_C1_js;
        
        double tmp_C2_js = 0;
        for(ui i = 0; i < C2.size(); i++) {
            ui u = C2[i];
            for(ui j = i + 1; j < C2.size(); j++) {
                ui v = C2[j];
                tmp_C2_js += js(u, v);
                
            }
        }
        allpair = (C2.size()*(C2.size()-1))/2;
        tmp_C2_js = (double) tmp_C2_js / allpair;
        
        total_C2_js += tmp_C2_js;
    }
    cout<<"Average Jaccard similarity = "<<((double)total_C1_js/(double)(r.size()) + (double)total_C2_js/(double)(r.size()) ) / 2 <<endl;
}

int main(int argc, const char * argv[]) {
    
    string graph_name = argv[1];
    load_graph_binary(graph_name);
    
    cout<<"*** ";
#ifdef _SHRINK_
    cout<<" _SHRINK_ ,";
#else
    cout<<" NO _SHRINK_ ,";
#endif
#ifdef _JSUB_
    cout<<" _JSUB_ ,";
#else
    cout<<" NO _JSUB_ ,";
#endif
#ifdef _CheckResults_
    cout<<" _CheckResults_ ";
#else
    cout<<" NO _CheckResults_ ";
#endif
    cout<<" ***"<<endl;

    seg_times = atof(argv[3]);  //control the number of segments
    assert(seg_times >= 0.01 && seg_times <= 100);
    
    rg_limit = atof(argv[4]);  //steady threshold
    assert(rg_limit >= 0.01 && rg_limit <= 1);

    int build_idx = atoi(argv[2]);
    if(build_idx == 1){
        assert(argc >= 6);
        Timer t;

        string way = argv[5];
        
        if (way.compare("LG") == 0) build_index_LG(graph_name);  //build index LG
        else if (way.compare("SS") == 0) build_index_SS(graph_name);  //build index SS
        else cout<<"No matching build index way."<<endl;
        
        cout<<"Finish building index, time cost = "<<integer_to_string(t.elapsed())<<endl;
        dele_memo();
        return 0;
    }

    int load_idx = atoi(argv[6]);
    if(load_idx == 1){
        assert(argc >= 8);
        Timer t;
        
        string index_name = argv[7];
        
        if (index_name.compare("LG") == 0) load_index_LG(graph_name);  //build index LG
        else if (index_name.compare("SS") == 0) load_index_SS(graph_name);  //build index SS
        else cout<<"No matching index name."<<endl;
        
        cout<<"Finish loading index, time cost = "<<integer_to_string(t.elapsed())<<endl;
    }

    assert(argc == 10);

    vr_method = atoi(argv[8]);  //1: vertex reduction, 2: index based vertex reduction

    epsilon = atof(argv[9]);  //epsilon: similarity threshold
    assert(epsilon >= 0 && epsilon <= 1);
    cout<<"epsilon = "<<epsilon<<", ";

    tau = atoi(argv[10]);  //tau: size constraint
    assert(tau >= 2);
    cout<<"tau = "<<tau<<endl;
    
    Timer t;

    MSBC_Enum();

    cout<<" - - - - - - - - - - - - - - - - - - - - - - - - -  "<<endl;
    if(over_time_flag) cout<<"| ###### OVER TIME ######"<<endl;
    cout<<"| Number of similar-bicliques = "<<MSBC_num<<endl;
    if(MSBC_num != 0){
        cout<<"|    min size = "<<min_MSBC_size<<endl;
        cout<<"|    max size = "<<max_MSBC_size<<endl;
    }
    cout<<"| VR time cost = "<<integer_to_string(time_on_vr)<<endl;
    cout<<"| ENUM time cost = "<<integer_to_string(time_on_enum)<<endl;
    cout<<"| TOTAL time cost = "<<integer_to_string(t.elapsed())<<endl;
    cout<<" - - - - - - - - - - - - - - - - - - - - - - - - -  "<<endl;

#ifdef _CheckResults_
    check_results(results);
#endif

        
#ifdef _PrintResults_
    cout<<"Results : "<<endl;
    for(auto e : results) {
        cout<<"\tCL : "; for(auto x : e.first) cout<<x<<", "; cout<<endl;
        cout<<"\tCR : "; for(auto x : e.second) cout<<x<<", "; cout<<endl;
        cout<<"\t----------------------------------------"<<endl;
    }
#endif
    
    dele_memo();
    
    return 0;
}
