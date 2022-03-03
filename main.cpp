//
//  main.cpp
//  MSBE
//
//  Created by kai on 20/12/2021.
//

#include "Timer.h"
#include "Utility.h"
#include "LinearHeap.h"

//the original graph
ui n, n1, n2, m; //# of vertices, # of L-side vertices, # of R-side vertices, # of edges
ui * pstart; //adjacent array
ui * edges; //adjacent array
int * degree; //vertex degree
int * TMPdeg; //used in vertex reduction
int * Sdeg; //used in vertex reduction
int * os; //vertex side

//the remaining graph
ui r_n, r_n1, r_n2, r_m; //# of vertices, # of L-side vertices, # of R-side vertices, # of edges
ui * r_pstart; //adjacent array
ui * r_edges; //adjacent array
int * r_degree; //vertex degree
ui * oid; //original vertex id
ui * nid; //new vertex id

//the subgraph in enumeration procedure
ui max_CS_deg; //used to construct Matrix
int ** Matrix; //adjacent Matrix
ui * trans; //mapping vertex id in Matrix
bool * del_ver; //record vertex status
int * inCR; //in CR or not
int * domCR; //domination set
int * deg_inCR; //degree in CR

//input parameters
double epsilon; //similarity constraint
int tau; //size constraint
double seg_num_times; //control the number of segments, i.e., \alpha
double rg_limit; //steady threshold, i.e., \gamma

int build_index_flag;
int load_index_flag;

vector<pair<vector<ui>, vector<ui>>> results;
long long MDBC_num = 0;
int max_MDBC_size = 0;
int min_MDBC_size = INF;

bool cmp_degCR (ui u, ui v)
{
    return deg_inCR[u] < deg_inCR[v];
}

class Itval
{
public:
        ui s_idx; //true vertex id
        ui e_idx; //true vertex id
        double min_score;
        double max_score;
        int c; //count
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
    ui coreV; //vertex id
    int rgC; //L vertex id -> R vertex id
    ui Lidx; //L position
    ui Ridx; //R position
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

void load_graph_bin(string graph_name)
{
    cout<<"\tStart reading graph: "<<graph_name<<endl;
    FILE *f = Utility::open_file(graph_name.c_str(), "rb");

    ui tt;
    fread(&tt, sizeof(ui), 1, f);
    if(tt != sizeof(ui)) {
        printf("sizeof unsigned int is different: .bin(%u), machine(%lu)\n", tt, sizeof(ui));
        exit(1);
    }
    fread(&n1, sizeof(ui), 1, f);
    fread(&n2, sizeof(ui), 1, f);
    n = n1 + n2;
    fread(&m, sizeof(ui), 1, f);
    cout<<"\tn1 = "<<n1<<", n2 = "<<n2<<", m = "<<m<<endl;
    
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
    cout<<"\tStart reading graph: "<<graph_name<<endl;
    ifstream input_file(graph_name, ios::in);
    map<ui, set<ui>> biG;
    if (!input_file.is_open()){
        cout << "\tCannot open file : "<<graph_name<<endl;exit(1);
    }
    else{
        input_file >> n1 >> n2 >> m;
        n = n1 + n2;
        cout<<"\tn1 = "<<n1<<", n2 = "<<n2<<", m = "<<m<<endl;
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
        input_file.close();
    }
    
    pstart = new ui[n+1];
    edges = new ui[2*m];
    degree = new int[n];
    TMPdeg = new int[n];
    Sdeg = new int[n];
    os = new int[n];
        
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
    for(ui i = 0; i < n1; i++) os[i] = 1;
    for(ui i = n1; i < n; i++) os[i] = 2;
}

void load_index_LG(string graph_name)
{
    rg_limit = 0;
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_" + to_string((int)(seg_num_times*100)));
    graph_name.append("_" + to_string((int)(rg_limit*100)));
    graph_name.append("_LG.bin");

    cout<<"\tLoad index: "<<graph_name<<endl;
    
    FILE * f = Utility::open_file(graph_name.c_str(), "rb");
    vsn.resize(n);
    
    int vid, numofI;
    ui s_vid, e_vid;
    char minS, maxS;
    int cnt;
    
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
        else {
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

void load_index_SS(string graph_name)
{
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_" + to_string((int)(seg_num_times*100)));
    graph_name.append("_" + to_string((int)(rg_limit*100)));
    graph_name.append("_SS.bin");

    cout<<"\tLoad index: "<<graph_name<<endl;
    
    FILE * f = Utility::open_file(graph_name.c_str(), "rb");
    vsn.resize(n);
    
    int vid, numofI;
    ui s_vid, e_vid;
    char minS, maxS;
    int cnt;
    
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
        else {
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

void vr()
{
    cout<<"\tStart vertex reduction: "<<endl;
    Timer tt;
    assert(tau >= 2);
    
    ui * inQ = new ui[n];
    memset(inQ, 0, sizeof(ui)*n);
    
    queue<ui> Q;
    for(ui i = 0; i < n; i++) if((i < n1 && (get_sim_nei2(i).size() < tau - 1)) || TMPdeg[i] < tau ) {
        Q.push(i);
        inQ[i] = 1;
    }
    int delcnt = 0;
    while (!Q.empty()) {
        ui u = Q.front();
        Q.pop();
        del_ver[u] = 1;
        delcnt++;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            if(del_ver[v] == 1 || inQ[v] == 1) continue;
            -- TMPdeg[v];
            if(TMPdeg[v] == tau - 1) {
                Q.push(v);
                inQ[v] = 1;
            }
        }
        if(os[u] == 2) continue;
        vector<ui> sim_nei = get_sim_nei2(u);
        for(ui w : sim_nei){
            assert(!del_ver[w]);
            if(inQ[w] == 1) continue;
            if(get_sim_nei2(w).size() < tau - 1) {
                Q.push(w);
                inQ[w] = 1;
            }
        }
    }
    int cnt = 0;
    for(ui i = 0; i < n; i++){
        if(!del_ver[i]) ++ cnt;
    }
    for(ui i = 0; i < n1; i++){
        if(!del_ver[i]) {
            Sdeg[i] = get_sim_nei2(i).size();
        }
    }
    int num_of_vertice_inR = 0;
    for(ui u = n1; u < n; u++) if(!del_ver[u]) ++ num_of_vertice_inR;
    -- num_of_vertice_inR;
    for(ui u = n1; u < n; u++) if(!del_ver[u]) Sdeg[u] = num_of_vertice_inR;
    cout<<"\t    Delete "<<delcnt<<" vertices."<<endl;
    delete [] inQ;
}

void index_vr()
{
    cout<<"\tStart index based vertex reduction: "<<endl;
    if(vsn.empty()) {
        cout<<"\tvsn is empty."<<endl; exit(1);
    }
    assert(tau >= 2);
    
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
    for(ui u = 0; u < n1; u++) if(!del_ver[u]) {
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
    cout<<"\t    (Phase I). Delete "<<delcnt<<" vertices, time cost: "<<integer_to_string(tt.elapsed())<<endl;
    tt.restart();
        
    for(ui u = 0; u < n1; u++) if(!del_ver[u]) {
        assert(TMPdeg[u] >= tau);
        vector<Itval> & cand_list = vsn[u];
        for(auto & e : cand_list) if(e.max_score >= epsilon) {
            for(ui i = e.s_idx; i <= e.e_idx; i++) if(!del_ver[i] && js(u, i) >= epsilon && i != u) {
                ++ Sdeg[u];
            }
        }
        if(Sdeg[u] < tau - 1) Q.push(u);
    }
    
    int num_of_vertice_inR = 0;
    for(ui u = n1; u < n; u++) if(!del_ver[u]) ++ num_of_vertice_inR;
    -- num_of_vertice_inR;
    for(ui u = n1; u < n; u++) if(!del_ver[u]) Sdeg[u] = num_of_vertice_inR;
        
    int phaseII_delv = 0;
    while (!Q.empty()) {
        ui u = Q.front(); Q.pop();
        del_ver[u] = 1; ++ delcnt; ++ phaseII_delv;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i];
            if(del_ver[v]) continue;
            if(TMPdeg[v] -- == tau && Sdeg[v] >= tau - 1) {
                Q.push(v);
            }
        }
        if(os[u]==2) continue;
        vector<Itval> & cand_list = vsn[u];
        for(auto & e : cand_list){
            if(e.max_score < epsilon) continue;
            for(ui i = e.s_idx; i <= e.e_idx; i++){
                if(del_ver[i]) continue;
                if(js(u, i) >= epsilon && i != u) {
                    if(Sdeg[i] -- == tau - 1 && TMPdeg[i] >= tau) {
                        Q.push(i);
                    }
                }
            }
        }
    }
    num_of_vertice_inR = 0;
    for(ui u = n1; u < n; u++) if(!del_ver[u]) ++ num_of_vertice_inR;
    -- num_of_vertice_inR;
    for(ui u = n1; u < n; u++) if(!del_ver[u]) Sdeg[u] = num_of_vertice_inR;
    
    cout<<"\t    (Phase II). Delete "<<phaseII_delv<<" vertices, time cost: "<<integer_to_string(tt.elapsed())<<endl;
}

void Enum_noRSim_adv_core(vector<ui> CL, vector<ui> CR, vector<ui> PL, vector<ui> QL)
{
    if(CL.size() + PL.size() < tau || CR.size() < tau) return;
    bool maximality = true;
    //for each vertex in PL and QL, obtain its degree in CR
    for(auto u : QL) {
        int d = 0;
        for(auto v : CR) {
            if(Matrix[trans[u]][trans[v]] == 1) ++ d;
        }
        if (d == CR.size()) {
            maximality = false;
#ifdef _ET_
            bool sim2PL = true;
            for(auto e : PL) {
                if(Matrix[trans[u]][trans[e]] != 1) {
                    sim2PL = false;
                    break;
                }
            }
            if(sim2PL == true) return;
#endif
        }
        deg_inCR[u] = d;
    }
    for(auto u : PL) {
        int d = 0;
        for(auto v : CR) {
            if(Matrix[trans[u]][trans[v]] == 1) ++ d;
        }
        if (d == CR.size()) maximality = false;
        deg_inCR[u] = d;
    }
    if(maximality == true) {
        if(CL.size() >= tau && CR.size() >= tau) {
#ifdef _CheckResults_
                vector<ui> C1, C2;
                for(auto e : CL) C1.push_back(oid[e]);
                for(auto e : CR) C2.push_back(oid[e]);
                sort(C1.begin(), C1.end(), less<>());
                sort(C2.begin(), C2.end(), less<>());
                if(C1[0] < C2[0]) results.push_back(make_pair(C1, C2));
                else results.push_back(make_pair(C2, C1));
#endif
                ++ MDBC_num;
                if(CL.size() + CR.size() > max_MDBC_size) max_MDBC_size = CL.size() + CR.size();
                if(CL.size() + CR.size() < min_MDBC_size) min_MDBC_size = CL.size() + CR.size();
        }
    }
    if(PL.empty()) return;
    
    sort(PL.begin(), PL.end(), cmp_degCR);
    
    vector<ui> dom_cand;
#ifdef _DOM_
    ui ustar;
    int dom_num = -1;
    int tmpd = -1;
    ui tmpv;
    for(auto e : QL) if(deg_inCR[e] > tmpd) {
        tmpd = deg_inCR[e]; tmpv = e;
    }
    for(auto e : PL) if(deg_inCR[e] > tmpd) {
        tmpd = deg_inCR[e]; tmpv = e;
    }
    assert(tmpd >= 0);
    vector<ui> mixvec;
    mixvec.push_back(tmpv);
    for(auto u : mixvec) {
        //u's similar neighbors in QL and PL
        vector<ui> simnei;
        for(auto v : PL) if(Matrix[trans[u]][trans[v]] == 1) simnei.push_back(v);
        //u's connection neighbors in CR
        vector<ui> connei;
        for(auto v : CR) if(Matrix[trans[u]][trans[v]] == 1) {
            connei.push_back(v); domCR[v] = 1;
        }
        vector<ui> dom_collection;
        for(auto v : simnei) { //i.e., select dominated vertex from u's simnei
            vector<ui> neiconnei;
            for(auto w : CR) if(Matrix[trans[v]][trans[w]] == 1) neiconnei.push_back(w);
            bool dom_flag = true;
            for(auto w : neiconnei) if(domCR[w] == 0) {
                dom_flag = false;
                break;
            }
            if(dom_flag == true) dom_collection.push_back(v);
        }
        if((int)dom_collection.size() > dom_num) {
            ustar = u;
            dom_num = (int) dom_collection.size();
            dom_cand = dom_collection;
        }
        for(auto e : connei) domCR[e] = 0;
        if(dom_num >= 0) break;
    }
    if(dom_num >= 0){
        unordered_set<ui> dom_set;
        for(auto e : dom_cand) dom_set.insert(e);
        vector<ui> tmpPL;
        for(auto e : PL) if(dom_set.find(e) == dom_set.end()) tmpPL.push_back(e);
        for(auto e : dom_cand) tmpPL.push_back(e);
        PL = tmpPL;
    }
#endif
        
    for(ui i = 0; i < PL.size() - dom_cand.size(); i++) {
        ui u = PL[i];
        vector<ui> new_CL = CL;
        new_CL.push_back(u);
        vector<ui> new_CR;
        for(auto v : CR) if(Matrix[trans[u]][trans[v]] == 1) new_CR.push_back(v);
        vector<ui> new_PL;
        for(ui j = i + 1; j < PL.size(); j++) {
            ui v = PL[j];
            if(Matrix[trans[u]][trans[v]] == 1) new_PL.push_back(v);
        }
        vector<ui> new_QL;
        for(auto e : QL) if(Matrix[trans[u]][trans[e]] == 1) new_QL.push_back(e);
        Enum_noRSim_adv_core(new_CL, new_CR, new_PL, new_QL);
        QL.push_back(u);
    }
}

void build_matrix_for_CRPLQL(vector<ui> &CR, vector<ui> &PL, vector<ui> &QL)
{
    for(ui i = 0; i < max_CS_deg; i++) memset(Matrix[i], 0, sizeof(int)*max_CS_deg);
    int idx = 0;
    for(auto e : CR) trans[e] = idx ++;
    for(auto e : PL) trans[e] = idx ++;
    for(auto e : QL) trans[e] = idx ++;
    //CR
    for(auto u : PL) {
        for(ui i = r_pstart[u]; i < r_pstart[u+1]; i++) {
            ui v = r_edges[i];
            if(inCR[v] == 1) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
            }
        }
    }
    for(auto u : QL) {
        for(ui i = r_pstart[u]; i < r_pstart[u+1]; i++) {
            ui v = r_edges[i];
            if(inCR[v] == 1) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
            }
        }
    }
    //L
    for(ui i = 0; i < PL.size(); i++){
        for(ui j = i + 1; j < PL.size(); j++){
            ui u = PL[i], v = PL[j];
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
            }
        }
    }
    for(ui i = 0; i < QL.size(); i++){
        for(ui j = i + 1; j < QL.size(); j++){
            ui u = QL[i], v = QL[j];
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
            }
        }
    }
    for(auto u : QL){
        for(auto v : PL){
            if(js(oid[u], oid[v]) >= epsilon) {
                Matrix[trans[u]][trans[v]] = 1;
                Matrix[trans[v]][trans[u]] = 1;
            }
        }
    }
    
}

void MDBC_Enum_noRSim_adv()
{
    //vertex reduction
    del_ver = new bool[n];
    memset(del_ver, 0, sizeof(bool)*n);
    
    if(load_index_flag == 1) index_vr();
    else vr();

    //rebuild remaining graph
    r_n = 0; r_n1 = 0; r_n2 = 0;
    r_m = 0;
    for(ui u = 0; u < n1; u++) if(del_ver[u] == 0) {
        ++ r_n; ++r_n1;
        r_m += TMPdeg[u];
    }
    for(ui u = n1; u < n; u++) if(del_ver[u] == 0) {
        ++ r_n; ++r_n2;
        r_m += TMPdeg[u];
    }
    assert(r_m%2 == 0); r_m /= 2;
    
    if(r_n == 0) return;
    
    oid = new ui[r_n];
    nid = new ui[n];
    
    ui vididx = 0;
    for(ui u = 0; u < n; u++) if(del_ver[u] == 0) {
        nid[u] = vididx;
        oid[vididx] = u;
        ++ vididx;
    }
    assert(r_n == vididx);
    
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
    
    max_CS_deg = 0;
    for(ui i = 0; i < r_n1; i++) if(Sdeg[oid[i]] + r_degree[i] > max_CS_deg)
        max_CS_deg = Sdeg[oid[i]] + r_degree[i];
    
    Matrix = new int * [max_CS_deg];
    for(ui i = 0; i < max_CS_deg; i++) Matrix[i] = new int [max_CS_deg];
    trans = new ui[r_n];
    
    inCR = new int[r_n];
    memset(inCR, 0, sizeof(int)*r_n);
    deg_inCR = new int[r_n];
    domCR = new int[r_n];
    memset(domCR, 0, sizeof(int)*r_n);
    
    cout<<"\tStart enumeration procedure: "<<endl;
    for(ui u = 0; u < r_n1; u++) {
        //find maximal similar-bicliques containing u
        vector<ui> CL;
        CL.push_back(u);

        vector<ui> CR;
        for(ui i = r_pstart[u]; i < r_pstart[u+1]; i++) {
            ui v = r_edges[i];
            CR.push_back(v);
            inCR[v] = 1;
        }
        
        vector<ui> PL, QL;
        
        if(load_index_flag == 1) {
            vector<Itval> & cand_list = vsn[oid[u]];
            for(auto & e : cand_list) if(e.max_score >= epsilon) {
                for(ui i = e.s_idx; i <= e.e_idx; i++) if(del_ver[i] == 0) {
                    if(js(oid[u], i) >= epsilon){
                        ui v = nid[i];
                        if(v > u) { PL.push_back(v); }
                        else if(v < u) { QL.push_back(v); }
                    }
                }
            }
        }
        else {
            vector<ui> vec = get_sim_nei2(oid[u]);
            for(auto & e : vec) {
                ui v = nid[e];
                if(v > u) { PL.push_back(v); }
                else if(v < u) { QL.push_back(v); }
            }
        }
        if((CL.size() + PL.size()) < tau || CR.size() < tau) {
            for(auto e : CR) inCR[e] = 0;
            continue;
        }
        build_matrix_for_CRPLQL(CR, PL, QL);
        Enum_noRSim_adv_core(CL, CR, PL, QL);
        for(auto e : CR) inCR[e] = 0;
    }
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
    if(os != nullptr){
        delete [] os;
        os = nullptr;
    }
    if(inCR != nullptr){
        delete [] inCR;
        inCR = nullptr;
    }
    if(domCR != nullptr){
        delete [] domCR;
        domCR = nullptr;
    }
    if(deg_inCR != nullptr){
        delete [] deg_inCR;
        deg_inCR = nullptr;
    }
}

void build_index_LG(string graph_name)
{
    cout<<"\tIn build_index_LG()"<<endl;
    rg_limit = 0;
    vector<int> INDEX_vid;
    vector<ui> INDEX_flag;
    vector<vector<Itval>> INDEX_list;
    
    Timer tt;
    long long T_find_2hopneis = 0;
    long long T_cal_sim_and_sort_for_each = 0;
    long long T_build_ranges = 0;
    long long total_2hop_nei_size = 0;
    long long total_shrinked_2hop_nei_size = 0;
    
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
            if(simscore >= TS) ordered_2hop_neis.push_back(make_pair(e, simscore));
            c[e] = 0;
        }
        total_shrinked_2hop_nei_size += ordered_2hop_neis.size();
        if(ordered_2hop_neis.empty()) continue;
        sort(ordered_2hop_neis.begin(), ordered_2hop_neis.end(), less<>()); //increasing order
        T_cal_sim_and_sort_for_each += tt.elapsed();
        tt.restart();
        //start to generate summary ranges
        assert(ordered_2hop_neis.size() >= 1);
        INDEX_vid.push_back(u);
        if(ordered_2hop_neis.size() < 2){
            vector<Itval> tmpV;
            for(auto e : ordered_2hop_neis){
                tmpV.push_back(Itval(e.first, e.first, e.second, e.second, 1));
            }
            INDEX_list.push_back(tmpV);
            INDEX_flag.push_back(1);
        }
        else{
            vector<Itval> tmpV;
            assert(ordered_2hop_neis.size() >= 2);
            int gc = (int)log2(ordered_2hop_neis.size());
            assert(gc < ordered_2hop_neis.size());
            gc = gc * seg_num_times;
            gc /= 2;
            if(gc < 1) gc = 1;
            if(gc > ordered_2hop_neis.size()-1) gc = ordered_2hop_neis.size()-1;
            priority_queue<pair<int, ui>, vector<pair<int, ui>>, greater<pair<int, ui>>> kset;  //min heap
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
        }
        T_build_ranges += tt.elapsed();
        tt.restart();
    }
    long long totalT = T_find_2hopneis + T_cal_sim_and_sort_for_each + T_build_ranges;
    cout<<"\t  (1). Find 2hop neis, time cost: "<<integer_to_string(T_find_2hopneis)<<" ( "<<((double)T_find_2hopneis/(totalT) )*100<<"% )"<<endl;
    cout<<"\t  (2). Comp similarity, time cost: "<<integer_to_string(T_cal_sim_and_sort_for_each)<<" ( "<<((double)T_cal_sim_and_sort_for_each/(totalT) )*100<<"% )"<<endl;
    cout<<"\t  (3). Build segments, time cost: "<<integer_to_string(T_build_ranges)<<" ( "<<((double)T_build_ranges/(totalT) )*100<<"% )"<<endl;
    assert(INDEX_vid.size() == INDEX_list.size());
    assert(INDEX_vid.size() == INDEX_flag.size());
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_" + to_string((int)(seg_num_times*100)));
    graph_name.append("_" + to_string((int)(rg_limit*100)));
    graph_name.append("_LG.bin");
    FILE * f = Utility::open_file(graph_name.c_str(), "wb");
    for(ui i = 0; i < INDEX_vid.size(); i++) {
        fwrite(&INDEX_vid[i], sizeof(int), 1, f);
        int num = INDEX_list[i].size();
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
    cout<<"\tIn build_index_SS()"<<endl;
    
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
            if(simscore >= TS) ordered_2hop_neis.push_back(make_pair(e, simscore));
            c[e] = 0;
        }
        total_shrinked_2hop_nei_size += ordered_2hop_neis.size();
        if(ordered_2hop_neis.empty()) continue;
        sort(ordered_2hop_neis.begin(), ordered_2hop_neis.end(), less<>()); //increasing order
        T_cal_sim_and_sort_for_each += tt.elapsed();
        tt.restart();
        INDEX_vid.push_back(u);
        if(ordered_2hop_neis.size() < 12){
            vector<Itval> tmpV;
            for(auto e : ordered_2hop_neis){
                tmpV.push_back(Itval(e.first, e.first, e.second, e.second, 1));
            }
            INDEX_list.push_back(tmpV);
            INDEX_flag.push_back(1);
        }
        else{
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

            //select rc segments
            int rc = (int)log2(ordered_2hop_neis.size());
            assert(rc < ordered_2hop_neis.size());
            rc = rc * seg_num_times;
            rc /= 2;
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
            T_range_tree += ttt.elapsed();
            ttt.restart();
            t1_partition += t1.elapsed();
            t1.restart();
        }
        T_build_ranges += tt.elapsed();
        tt.restart();
    }
    long long totalT = T_find_2hopneis + T_cal_sim_and_sort_for_each + T_build_ranges;
    cout<<"\t  (1). Find 2hop neis, time cost: "<<integer_to_string(T_find_2hopneis)<<" ( "<<((double)T_find_2hopneis/(totalT) )*100<<"% )"<<endl;
    cout<<"\t  (2). Comp similarity, time cost: "<<integer_to_string(T_cal_sim_and_sort_for_each)<<" ( "<<((double)T_cal_sim_and_sort_for_each/(totalT) )*100<<"% )"<<endl;
    cout<<"\t  (3). Build segments, time cost: "<<integer_to_string(T_build_ranges)<<" ( "<<((double)T_build_ranges/(totalT) )*100<<"% )"<<endl;
    long long total_tree_T = T_build_tree + T_cp_tree + +T_destroy_tree + T_range_tree;
    cout<<"\t    (3.1). Build tree, time cost: "<<integer_to_string(T_build_tree)<<" ( "<<((double)T_build_tree/(total_tree_T) )*100<<"% )"<<endl;
    cout<<"\t    (3.2). Steady seg, time cost: "<<integer_to_string(T_cp_tree)<<" ( "<<((double)T_cp_tree/(total_tree_T) )*100<<"% )"<<endl;
    cout<<"\t    (3.3). Destroy tree, time cost: "<<integer_to_string(T_destroy_tree)<<" ( "<<((double)T_destroy_tree/(total_tree_T) )*100<<"% )"<<endl;
    cout<<"\t    (3.4). Build segments, time cost: "<<integer_to_string(T_range_tree)<<" ( "<<((double)T_range_tree/(total_tree_T) )*100<<"% )"<<endl;
    assert(INDEX_vid.size() == INDEX_list.size());
    assert(INDEX_vid.size() == INDEX_flag.size());
    
    //binary version
    
    graph_name.erase(graph_name.end() - 4, graph_name.end());
    graph_name.append("_" + to_string((int)(seg_num_times*100)));
    graph_name.append("_" + to_string((int)(rg_limit*100)));
    graph_name.append("_SS.bin");
    
    FILE * f = Utility::open_file(graph_name.c_str(), "wb");
    for(ui i = 0; i < INDEX_vid.size(); i++) {
        fwrite(&INDEX_vid[i], sizeof(int), 1, f);
        int num = INDEX_list[i].size();
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
    cout<<"\t= = = = checking results = = = = \t"<<endl;
    if(r.empty()) {
        cout<<"\tNote that results is empty !!!"<<endl;
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
                    cout<<"\tFind a vertex in C2 that is not neighbors of u in C1."<<endl;
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
                    cout<<"\tFind dissimilar pairs in C1."<<endl;
                    flag = false;
                    exit(1);
                }
            }
        }
    }
    if(flag == false) cout<<"\t= = = = = = WRONG! = = = = = = \t"<<endl;
    else cout<<"\t= = = = = = CORRECT! = = = = = = \t"<<endl;
}

void Case_for_DBLP()
{
    ifstream inputf;
//    string user_info = "/Users/kai/pyWorkspace/extract_dblp/bi_author_map.txt";
    string user_info = "/data/MDBC/txt_graphs_and_old_index/mydblp/bi_author_map.txt";
    inputf.open(user_info);
    assert(inputf.is_open());
    unordered_map<int, string> user_profile;
    string buffer;
    while (!inputf.eof()) {
        getline(inputf, buffer);
        if(buffer.empty()) continue;
        vector<string> splits = split(buffer, "\t");
        assert(splits.size() == 2);
        ui uid = stoi(splits[1]);
        user_profile[uid] = splits[0];
    }
    inputf.close();
    
//    string mv_info = "/Users/kai/datasets/MDBC/case-study/ml-1m/movies_map_over3.txt";
    string mv_info = "/data/MDBC/txt_graphs_and_old_index/mydblp/bi_paper_map.txt";
    inputf.open(mv_info);
    assert(inputf.is_open());
    unordered_map<int, string> mv_profile;
    while (!inputf.eof()) {
        getline(inputf, buffer);
        if(buffer.empty()) continue;
        vector<string> splits = split(buffer, "\t");
        assert(splits.size() == 2);
        ui mvid = stoi(splits[1]);
        mv_profile[mvid] = splits[0];
    }
    inputf.close();
    
    cout<<"Results : "<<endl;
    for(auto e : results) {
        cout<<"\tCL : ";
        for(auto x : e.first) {
            cout<<user_profile[x]<<", ";
        }
        cout<<endl;
        
        cout<<"\tCR : ";
        for(auto x : e.second) {
            cout<<"\t "<<mv_profile[x]<<endl;
        }
        cout<<"\t----------------------------------------"<<endl;
    }
}

void Case_for_ul()
{
    ifstream inputf;
    string user_info = "/Users/kai/datasets/MDBC/source-datasets/unicodelang/country_map.txt";
    inputf.open(user_info);
    assert(inputf.is_open());
    unordered_map<int, string> user_profile;
    string buffer;
    while (!inputf.eof()) {
        getline(inputf, buffer);
        if(buffer.empty()) continue;
        vector<string> splits = split(buffer, "\t");
        assert(splits.size() == 2);
        ui uid = stoi(splits[0]);
        user_profile[uid] = splits[1];
    }
    inputf.close();
    
    string mv_info = "/Users/kai/datasets/MDBC/source-datasets/unicodelang/lan_map.txt";
    inputf.open(mv_info);
    assert(inputf.is_open());
    unordered_map<int, string> mv_profile;
    while (!inputf.eof()) {
        getline(inputf, buffer);
        if(buffer.empty()) continue;
        vector<string> splits = split(buffer, "\t");
        assert(splits.size() == 2);
        ui mvid = stoi(splits[0]);
        mv_profile[mvid] = splits[1];
    }
    inputf.close();
    cout<<"---------------------case study---------------------"<<endl;
    for(auto C : results) {
        cout<<"C : "<<endl;
        cout<<"\tCL ("<<C.first.size()<<") : "; for(auto e : C.first) cout<<user_profile[e]<<","; cout<<endl;
        cout<<"\tCR ("<<C.second.size()<<") : "; for(auto e : C.second) cout<<mv_profile[e]<<","; cout<<endl;
    }
}

void print_results()
{
    cout<<"\tResults : "<<endl;
    for(auto e : results) {
        cout<<"\tCL : "; for(auto x : e.first) cout<<x<<", "; cout<<endl;
        cout<<"\tCR : "; for(auto x : e.second) cout<<x<<", "; cout<<endl;
        cout<<"\t----------------------------------------"<<endl;
    }
}

int main(int argc, const char * argv[]) {
    
    string graph_name = argv[1];
    load_graph_bin(graph_name); //load graph in binary format
//    load_graph_txt(graph_name); //load graph in txt format
    
    build_index_flag = atoi(argv[2]);
    if(build_index_flag == 1){
        assert(argc >= 6);
        
        string index_name = argv[3];
        
        seg_num_times = atof(argv[4]);
        assert(seg_num_times >= 0);
        
        rg_limit = atof(argv[5]);
        assert(rg_limit >= 0 && rg_limit <= 1);
        
        Timer t;
        
        if (index_name.compare("LG") == 0) build_index_LG(graph_name);
        else if (index_name.compare("SS") == 0) build_index_SS(graph_name);
        else cout<<"\tNo matched algorithm."<<endl;
        
        cout<<"\tBuilding index time cost: "<<integer_to_string(t.elapsed())<<endl;
        
        dele_memo();
        
        return 0;
    }

    load_index_flag = atoi(argv[6]);
    if(load_index_flag == 1){
        assert(argc >= 10);
        
        string index_name = argv[7];
        
        seg_num_times = atof(argv[8]);
        assert(seg_num_times >= 0);
        
        rg_limit = atof(argv[9]);
        assert(rg_limit >= 0 && rg_limit <= 1);
        
        if (index_name.compare("LG") == 0) load_index_LG(graph_name);
        else if (index_name.compare("SS") == 0) load_index_SS(graph_name);
        else cout<<"\tNo matched algorithm."<<endl;
    }

    assert(argc == 12);

    epsilon = atof(argv[10]);
    assert(epsilon >= 0 && epsilon <= 1);
    
    tau = atoi(argv[11]);
    assert(tau >= 1);
    
    cout<<"\tepsilon="<<epsilon<<", "<<"tau="<<tau<<endl;
    
    Timer t;
    MDBC_Enum_noRSim_adv();
    cout<<"\t------------------------------------------------"<<endl;
    cout<<"\tNumber of Similar-Bicliques: "<<MDBC_num<<" (min size: "<<min_MDBC_size<<", max size: "<<max_MDBC_size<<")."<<endl;
    cout<<"\tTime cost (without I/O): "<<integer_to_string(t.elapsed())<<endl;
    cout<<"\t------------------------------------------------"<<endl;
    
#ifdef _CheckResults_
    check_results(results);
//    print_results();
#endif
    
//    Case_for_DBLP();
//    Case_for_ul();
    
    dele_memo();
    
    return 0;
}
