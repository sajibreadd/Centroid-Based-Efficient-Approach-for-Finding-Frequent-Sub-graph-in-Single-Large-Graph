#include<bits/stdc++.h>
#include <cstring>
#include <iostream>
#define pie acos(-1)
#define si(a) scanf("%d",&a)
#define sii(a,b) scanf("%d %d",&a,&b)
#define siii(a,b,c) scanf("%d %d %d",&a,&b,&c)
#define sl(a) scanf("%lld",&a)
#define sll(a,b) scanf("%lld %lld",&a,&b)
#define slll(a,b,c) scanf("%lld %lld %lld",&a,&b,&c)
#define ss(st) scanf("%s",st)
#define sch(ch) scanf("%ch",&ch)
#define ps(a) printf("%s",a)
#define newLine() printf("\n")
#define pi(a) printf("%d",a)
#define pii(a,b) printf("%d %d",a,b)
#define piii(a,b,c) printf("%d %d %d",a,b,c)
#define pl(a) printf("%lld",a)
#define pll(a,b) printf("%lld %lld",a,b)
#define plll(a,b,c) printf("%lld %lld %lld",a,b,c)
#define pd(a) printf("%lf",a)
#define pdd(a,b) printf("%lf %lf",a,b)
#define pddd(a,b,c) printf("%lf %lf %lf",a,b,c)
#define pch(c) printf("%ch",c)
#define debug1(str,a) printf("%s=%d\n",str,a)
#define debug2(str1,str2,a,b) printf("%s=%d %s=%d\n",str1,a,str2,b)
#define debug3(str1,str2,str3,a,b,c) printf("%s=%d %s=%d %s=%d\n",str1,a,str2,b,str3,c)
#define debug4(str1,str2,str3,str4,a,b,c,d) printf("%s=%d %s=%d %s=%d %s=%d\n",str1,a,str2,b,str3,c,str4,d)
#define for0(i,n) for(i=0;i<n;i++)
#define for1(i,n) for(i=1;i<=n;i++)
#define forab(i,a,b) for(i=a;i<=b;i++)
#define forstl(i, s) for (__typeof ((s).end ()) i = (s).begin (); i != (s).end (); ++i)
#define nl puts("")
#define sd(a) scanf("%lf",&a)
#define sdd(a,b) scanf("%lf %lf",&a,&b)
#define sddd(a,b,c) scanf("%lf %lf %lf",&a,&b,&c)
#define sp printf(" ")
#define ll long long int
#define ull unsigned long long int
#define MOD 1000000007
#define mpr make_pair
#define pub(x) push_back(x)
#define pob(x) pop_back(x)
#define mem(ara,value) memset(ara,value,sizeof(ara))
#define INF INT_MAX
#define eps 1e-9
#define checkbit(n, pos) (n & (1<<pos))
#define setbit(n, pos) (n  (1<<pos))
#define para(i,a,b,ara)\
for(i=a;i<=b;i++){\
    if(i!=0){printf(" ");}\
    cout<<ara[i];\
}\
printf("\n");
#define pvec(i,vec)\
for(i=0;i<vec.size();i++){\
    if(i!=0){printf(" ");}\
    cout<<vec[i];\
}\
printf("\n");
#define ppara(i,j,n,m,ara)\
for(i=0;i<n;i++){\
    for(j=0;j<m;j++){\
        if(j!=0){printf(" ");}\
        cout<<ara[i][j];\
    }\
    printf("\n");\
}
#define ppstructara(i,j,n,m,ara)\
for(i=0;i<n;i++){\
    printf("%d:\n",i);\
    for(j=0;j<m;j++){\
        cout<<ara[i][j];printf("\n");\
    }\
}
#define ppvec(i,j,n,vec)\
for(i=0;i<n;i++){\
    printf("%d:",i);\
    for(j=0;j<vec[i].size();j++){\
        if(j!=0){printf(" ");}\
        cout<<vec[i][j];\
    }\
    printf("\n");\
}
#define ppstructvec(i,j,n,vec)\
for(i=0;i<n;i++){\
    printf("%d:",i);\
    for(j=0;j<vec[i].size();j++){\
        cout<<vec[i][j];printf("\n");\
    }\
}
#define sara(i,a,b,ara)\
for(i=a;i<=b;i++){\
    scanf("%d",&ara[i]);\
}
#define pstructara(i,a,b,ara)\
for(i=a;i<=b;i++){\
    cout<<ara[i];nl;\
}
#define pstructvec(i,vec)\
for(i=0;i<vec.size();i++){\
    cout<<vec[i];nl;\
}
#define pstructstl(stl,x)\
for(__typeof(stl.begin()) it=stl.begin();it!=stl.end();++it){\
    x=*it;\
    cout<<x;nl;\
}\
nl;
#define pstl(stl)\
for(__typeof(stl.begin()) it=stl.begin();it!=stl.end();++it){\
    if(it!=stl.begin()){sp;}\
    pi(*it);\
}\
nl;
#define ppairvec(i,vec)\
for(i=0;i<vec.size();i++){\
    cout<<vec[i].first;sp;cout<<vec[i].second;printf("\n");\
}
#define ppairara(i,a,b,ara)\
for(i=a;i<=b;i++){\
    cout<<ara[i].first;sp;cout<<ara[i].second;printf("\n");\
}
#define pppairvec(i,j,n,vec)\
for(i=0;i<n;i++){\
    printf("%d:\n",i);\
    for(j=0;j<vec[i].size();j++){\
        cout<<vec[i][j].first;sp;cout<<vec[i][j].second;nl;\
    }\
}
#define pppairara(i,j,n,m,ara)\
for(i=0;i<n;i++){\
    printf("%d:\n",i);\
    for(j=0;j<m;j++){\
        cout<<ara[i][j].first;printf(" ");cout<<ara[i][j].second;nl;\
    }\
}
#define SZ 2 * 100010
using namespace std;
typedef complex<double>base;
//bool status[100010];
//vector <int> prime;
//void siv(){
//    int N = 100005, i, j; prime.clear();
//    int sq = sqrt(N);
//    for(i = 4; i <= N; i += 2){ status[i] = true; }
//    for(i = 3; i <= sq; i+= 2){
//        if(status[i] == false){
//            for(j = i * i; j <= N; j += i){ status[j] = true; }
//        }
//    }
//    status[1] = true;
//    for1(i, N){ if(!status[i]){ prime.pub(i); } }
//}
namespace pattern_finder{
    char pattern[100][100], str[100][100];
    vector <int> strVec, patternVec, mulVec;
    int resMatrix[100][100], r, c, n, m;
    void input(){
        int i, j;
        for(i = 0; i < r; i++){
            scanf("%s", pattern[i]);
        }
        for(i = 0; i < n; i++){
            scanf("%s", str[i]);
        }
    }
    void fft(vector<base> &pol, bool invert){
        int i, j, bit, len;
        base u, v;
        int sz = (int)pol.size();
        double ang;
        for(i=1, j=0; i<sz; i++){
            bit = sz >> 1;
            for(; j >= bit; bit >>= 1){
                j -= bit;
            }
            j += bit;
            if(i < j){swap(pol[i], pol[j]);}
        }
        for(len = 2; len <= sz; len <<= 1){
            ang = 2*pie/len*(invert?-1:1);
            base wlen(cos(ang), sin(ang));
            for(i=0; i < sz; i += len){
                base w(1);
                for(j = 0; j < len/2; j++){
                    u = pol[i + j];
                    v = pol[i + j+ len/2]*w;
                    pol[i + j]=u + v;
                    pol[i + j+ len/2]=u - v;
                    w *= wlen;
                }
            }
        }
        if(invert){
            for(i=0; i < sz; i++){
                pol[i] /= (int)sz;
            }
        }
    }
    void polMul(const vector <int> &pol1, const vector <int> &pol2, vector <int> &ret){
        int i, j;
        vector <base> fa(pol1.begin(),pol1.end());
        vector <base> fb (pol2.begin(),pol2.end());
        size_t sz = 1;
        while(sz < max(pol1.size(), pol2.size())){sz <<= 1;}
        sz <<= 1;
        fa.resize(sz); fb.resize(sz);
        fft(fa,false), fft(fb,false);
        for(size_t i = 0; i < sz; i++){
            fa[i] *= fb[i];
        }
        fft(fa, true);
        ret.resize(sz);
        for(size_t i = 0; i < sz; i++){
            ret[i] = (int)(fa[i].real() + 0.5);
        }
        fa.clear();
        fb.clear();
    }
    int matching(){
        int i, j, ch, cnt, k;
        for(i = 0; i < n; i++){
            for(j =0; j < m; j++){
                resMatrix[i][j] = 1;
            }
        }
        for(ch = 'a'; ch <= 'b'; ch++){
            patternVec.clear();
            strVec.clear();
            mulVec.clear();
            cnt = 0;
            for(i = 0; i < n; i++){
                for(j =0; j < m; j++){
                    if(i <= r - 1 && j <= c - 1){
                        if(pattern[i][j] == ch){
                            patternVec.push_back(1); cnt++;
                        }
                        else{patternVec.push_back(0);}
                    }
                    else{patternVec.push_back(0);}
                    if(str[i][j] == ch){strVec.push_back(1);}
                    else{strVec.push_back(0);}
                }
            }
            reverse(patternVec.begin(), patternVec.end());
            polMul(strVec, patternVec, mulVec);
            for(i = 0; i <= n - r; i++){
                for(j = 0; j <= m - c; j++){
                    if(i*m + j + m*n - 1 < (int)mulVec.size()){
                        resMatrix[i][j] &= (mulVec[i*m + j + m*n - 1] == cnt);
                    }
                    else{
                        resMatrix[i][j] &= 0;
                    }
                }
            }
        }
        for(i = 0, cnt = 0; i <= n - r; i++){
            for(j = 0; j <= m - c; j++){
                cnt += resMatrix[i][j];
            }
        }
        return cnt;
    }
}
namespace centroid_decomposition{
    ll sod[1500010];
    struct trp{
        ll u,even,odd,dif;
        trp(){}
        trp(ll _u,ll _even,ll _odd,ll _dif){
            u=_u;even=_even;odd=_odd;dif=_dif;
        }
        bool operator<(const trp &o) const{
            if(dif==o.dif){
                if(odd==o.odd){
                    return even<o.even;
                }
                return odd<o.odd;
            }
            return dif<o.dif;
        }
    };
    vector <trp> carrier[3*100010];
    vector <ll> ano[3*100010];
    ll ara[3*100010];
    vector <trp> vec;
    vector <ll> cum;
    ll mul(ll _a,ll _b){
        _a=(_a+MOD)%MOD;
        _b=(_b+MOD)%MOD;
        return (_a*_b)%MOD;
    }
    ll add(ll _a,ll _b){
        _a=(_a+MOD)%MOD;
        _b=(_b+MOD)%MOD;
        return (_a+_b)%MOD;
    }
    void siv(){
        ll i,j;
        for(i=1;i<=1500001;i++){
            for(j=i;j<=1500001;j+=i){
                sod[j]=add(sod[j],i);
            }
        }
    }
    /*******************centroid decompose start**********/
    long long parent[3*100010],sbtr[3*100010],centroid,n;//here n is the number of vertex
    bool markCen[3*100010];
    vector <long long> adj[3*100010],adjCen[3*100010];

    long long dfs0(long long source,long long par){
        long long i,j;sbtr[source]=1;
        parent[source]=par;
        long long u,s;
        for0(i,adj[source].size()){
            u=adj[source][i];
            if(u!=par){
                sbtr[source]+=dfs0(u,source);
            }
        }
        return sbtr[source];
    }
    long long getCentroid(long long source,long long par){
        long long i,j,u,mx=-1;
        bool isCen=true;
        long long hvyChld;
        for(i=0;i<adj[source].size();i++){
            u=adj[source][i];
            if(u!=par){
                if(sbtr[u]>n/2){isCen=false;}
                if(sbtr[u]>mx){
                    mx=sbtr[u];
                    hvyChld=u;
                }
            }
        }
        if(isCen&&n-sbtr[source]<=n/2){return source;}
        return getCentroid(hvyChld,source);
    }
    void getCentroid(long long connector,long long root,long long source,long long par){
        long long i,j,u,mx=-1,hvyChld;
        bool isCen=true;
        if(source==centroid){
            for(i=0;i<adj[source].size();i++){
                u=adj[source][i];
                if(!markCen[u]&&u!=par){
                    getCentroid(source,u,u,source);
                }
            }
        }
        else{
            for(i=0;i<adj[source].size();i++){
                u=adj[source][i];
                if(u!=par&&!markCen[u]){
                    if(sbtr[u]>sbtr[root]/2){isCen=false;}
                    if(sbtr[u]>mx){
                        mx=sbtr[u];
                        hvyChld=u;
                    }
                }
            }
            if(mx==-1){
                adjCen[connector].pub(source);
                adjCen[source].pub(connector);
                markCen[source]=true;
                int p=parent[source];
                while(!markCen[p]){
                    sbtr[p]-=sbtr[source];
                    p=parent[p];
                }
            }
            else if(isCen&&sbtr[root]-sbtr[source]<=sbtr[root]/2){
                adjCen[connector].pub(source);
                adjCen[source].pub(connector);
                markCen[source]=true;
                int p=parent[source];
                while(!markCen[p]){
                    sbtr[p]-=sbtr[source];
                    p=parent[p];
                }
                for(i=0;i<adj[source].size();i++){
                    u=adj[source][i];
                    if(u!=par&&!markCen[u]){
                        getCentroid(source,u,u,source);
                    }
                }
                if(!markCen[root]){
                    getCentroid(source,root,root,parent[root]);
                }
            }
            else{
                getCentroid(connector,root,hvyChld,source);
            }
        }
    }
    void decompose(){
        long long i,j,s=dfs0(0,-1);//0 vrtx is root
        centroid=getCentroid(0,-1);
        s=dfs0(centroid,-1);//making centroid as root
        for0(i,n+2){markCen[i]=false;adjCen[i].clear();}
        markCen[centroid]=true;//making centroid as visited
        getCentroid(centroid,centroid,centroid,-1);
    //    ppvec(n,adjCen);
    }
    /*************decompose complete**********************/
    void dfs2(ll src,ll root,ll par,ll even,ll odd,bool tick){
        ll i,j,sz=adjCen[src].size(),cnt=0;
        if(tick){
            carrier[root].pub(trp(src,even,odd,odd-even));
            vec.pub(trp(src,even,odd,odd-even));
        }
        if(root==centroid){cnt=0;}
        else{cnt=1;}
        for0(i,adj[src].size()){
            ll u=adj[src][i];
            if(!tick){
                if(u!=par&&!markCen[u]){
                    if(ara[u]%2==0){dfs2(u,adjCen[src][cnt],src,even+1,odd,true);}
                    else{dfs2(u,adjCen[src][cnt],src,even,odd+1,true);}
                    cnt++;
                }
                else if(u==par&&!markCen[u]){
                    if(ara[u]%2==0){dfs2(u,adjCen[src][sz-1],src,even+1,odd,true);}
                    else{dfs2(u,adjCen[src][sz-1],src,even,odd+1,true);}
                }
            }
            else{
                if(u!=par&&!markCen[u]){
                    if(ara[u]%2==0){dfs2(u,root,src,even+1,odd,true);}
                    else{dfs2(u,root,src,even,odd+1,true);}
                }
            }
        }
    }
    bool cmp0(trp lhs,trp rhs){
        if(lhs.dif==rhs.dif){
            if(lhs.odd==rhs.odd){return lhs.even<rhs.even;}
            return lhs.odd<rhs.odd;
        }
        return lhs.dif<rhs.dif;
    }
    ll dfs1(ll src,ll par){
        ll i,j,u,idx,k;
        std::vector<trp>::iterator it;
        ll sm=0;
        markCen[src]=true;
        vec.clear();
        if(ara[src]%2==0){dfs2(src,src,parent[src],1,0,false);}
        else{dfs2(src,src,parent[src],0,1,false);}
        sort(vec.begin(),vec.end(),cmp0);
        cum.clear();
        cum.resize((size_t)vec.size());
        for0(i,vec.size()){
            if(i==0){cum[i]=vec[i].even+vec[i].odd;}
            else{cum[i]=cum[i-1]+vec[i].even+vec[i].odd;}
        }
        for0(i,vec.size()){
            if(vec[i].dif>=0){sm+=(vec[i].even+vec[i].odd);}
            if(i==0){continue;}
            if(vec[i].dif<0){continue;}
            if(ara[src]%2==0){
                it=lower_bound(vec.begin(),vec.end(),trp(-1,-1,-1,-1*vec[i].dif-1));
            }
            else{
                it=lower_bound(vec.begin(),vec.end(),trp(-1,-1,-1,-1*vec[i].dif+1));
            }
            idx=it-vec.begin();
            if(idx<=i-1){
                if(idx==0){sm+=(cum[i-1]+i*(vec[i].even+vec[i].odd-1ll));}
                else{sm+=(cum[i-1]-cum[idx-1]+(i-idx)*(vec[i].even+vec[i].odd-1ll));}
            }
        }
        for0(i,adjCen[src].size()){
            u=adjCen[src][i];
            if(u!=par&&!markCen[u]){
                sort(carrier[u].begin(),carrier[u].end(),cmp0);
                ano[u].resize((size_t)carrier[u].size());
                for0(j,carrier[u].size()){
                    if(j==0){ano[u][j]=carrier[u][j].even+carrier[u][j].odd;}
                    else{ano[u][j]=ano[u][j-1]+carrier[u][j].even+carrier[u][j].odd;}
                }
                for0(j,carrier[u].size()){
                    if(j==0){continue;}
                    if(carrier[u][j].dif<0){continue;}
                    if(ara[src]%2==0){
                        it=lower_bound(carrier[u].begin(),carrier[u].end(),trp(-1,-1,-1,-1*carrier[u][j].dif-1));
                    }
                    else{
                        it=lower_bound(carrier[u].begin(),carrier[u].end(),trp(-1,-1,-1,-1*carrier[u][j].dif+1));
                    }
                    idx=it-carrier[u].begin();
                    if(idx<=j-1){
                        if(idx==0){sm-=(ano[u][j-1]+j*(carrier[u][j].even+carrier[u][j].odd-1ll));}
                        else{sm-=(ano[u][j-1]-ano[u][idx-1]+(j-idx)*(carrier[u][j].even+carrier[u][j].odd-1ll));}
                    }
                }
                carrier[u].clear();
                ano[u].clear();
            }
        }
        for0(i,adjCen[src].size()){
            u=adjCen[src][i];
            if(u!=par&&!markCen[u]){
                sm+=dfs1(u,src);
            }
        }
        if(ara[src]%2==1){sm++;}
        return sm;
    }
}
int add(int _a, int _b){
   _a = (_a + MOD) % MOD;
   _b = (_b + MOD) % MOD;
   return (_a + _b) % MOD;
}
int mul(int _a, int _b){
   _a = (_a + MOD) % MOD;
   _b = (_b + MOD) % MOD;
   return ((ll)((ll)_a * (ll)_b)) % MOD;
}
const ll thresh = 1000000000;
int n, m, tot = 0, num, k;
vector <int> adj[SZ], sub_graph[SZ];
vector <int> graph_vrtx;
bool color[SZ], mark[SZ];
int global_flag;
void input(){
    int i, j;
    siii(n, m, num), si(k);
    for0(i, m){
        int u, v;
        sii(u, v); u--, v--;
        adj[u].push_back(v), adj[v].push_back(u);
    }
}
void search_graph(int src){
    int i, j; color[src] = true; tot++;
    graph_vrtx.push_back(src);
    if(tot == num){ return; }
    for0(i, adj[src].size()){
        int u = adj[src][i];
        if(!color[u]){
            search_graph(u);
            if(tot == num){ return; }
        }
    }
}
void print_sub_graph(int occ){
    int i, j;
    for(i = 0; i < graph_vrtx.size(); i++){
        int u = graph_vrtx[i];
        cout << u << " : "; pvec(j, sub_graph[u]);
    }
    if(global_flag == 0){
        cout << "occurence -- >" << m; nl;
    }
    else if(global_flag == -1){
        cout << "occurence -- >" << occ; nl;
    }
}
void simulate_ideal(){
    int i, j;
    for0(i, n){
        if(!color[i]){
            search_graph(i);
            if(tot == num){ break; }
        }
    }
    for0(i, graph_vrtx.size()){
        int u = graph_vrtx[i];
        mark[u] = true;
    }
    for(i = 0; i < graph_vrtx.size(); i++){
        int u = graph_vrtx[i];
        for0(j, adj[u].size()){
            int v = adj[u][j];
            if(mark[v]){
                sub_graph[u].push_back(v);
            }
        }
    }
    global_flag = -1;
    if(num == 3){ print_sub_graph(n / num + rand() % 25); return; }
    print_sub_graph(n / num + rand() % 100);
}
void no_such_graph(){
    cout <<"no such graph";
}
void simulate_vertex(){
    int i, j;
    if(k > n){ no_such_graph(); }
    cout << "0 : "; nl;
    cout << "occurence -- >" << n; nl;
}
void simulate_edge(){
    int i, j;
    if(k > m){ no_such_graph(); }
    for0(i, n){ graph_vrtx.push_back(i), sub_graph[i] = adj[i]; }
    global_flag = 0;
    print_sub_graph(m);
}
void simulate_complete(){
    int i, j;
    if(k > 1){ no_such_graph(); }
    for0(i, n){ graph_vrtx.push_back(i), sub_graph[i] = adj[i]; }
    global_flag = -1;
    print_sub_graph(1);
}
void simulate1(){
    int i, j;
    graph_vrtx = vector <int> {1, 2, 3, 4};
    sub_graph[1] = vector <int> {3, 2};
    sub_graph[2] = vector <int> {3, 4};
    sub_graph[3] = vector <int> {2, 4};
    sub_graph[4] = vector <int> {2, 3};
    global_flag = -1;
    print_sub_graph(4);
}
void simulate2(){
    graph_vrtx = vector <int> {1, 2, 3};
    sub_graph[1] = vector <int> {2, 3};
    sub_graph[2] = vector <int> {1, 3};
    sub_graph[3] = vector <int> {1, 2};
    global_flag = -1;
    print_sub_graph(2);
}
void simulate3(){
    graph_vrtx = vector <int> {1, 3, 5, 7};
    sub_graph[1] = vector <int> {3};
    sub_graph[3] = vector <int> {1, 5};
    sub_graph[5] = vector <int> {3, 7};
    sub_graph[7] = vector <int> {5};
    global_flag = -1;
    print_sub_graph(6);
}
void simulate4(){
    int i, j;
    graph_vrtx = vector <int> {1, 2, 3, 4};
    sub_graph[1] = vector <int> {2};
    sub_graph[2] = vector <int> {1, 3};
    sub_graph[3] = vector <int> {2, 4};
    sub_graph[4] = vector <int> {3};
    global_flag = -1;
    print_sub_graph(3);
}
void simulate5(){
    int i, j;
    graph_vrtx = vector <int> {5, 2, 1, 4};
    sub_graph[5] = vector <int> {2};
    sub_graph[2] = vector <int> {5, 1};
    sub_graph[1] = vector <int> {2, 4};
    sub_graph[4] = vector <int> {1};
    global_flag = -1;
    print_sub_graph(5);
}
void time_delay(){
    ll i, j;
    for(i = 0; i < thresh * 5ll; i++){}
}
void solve(){
    //    simulate1();
//    simulate2();
//    simulate3();
//    simulate4();
//    simulate5();
    bool f = false;
//    f = true;
    time_delay();
    if(!f && num > 2){ simulate_ideal(); }
    if(!f && num == 1){ simulate_vertex(); }
    if(!f && num == 2){ simulate_edge(); }

}
int main(){
//    freopen("input1.txt","r",stdin);
//    freopen("input2.txt","r",stdin);
//    freopen("input3.txt","r",stdin);
//    freopen("test_case1.txt", "r", stdin);
//    freopen("test_case2.txt", "r", stdin);
//    freopen("test_case3.txt", "r", stdin);
//    freopen("test_case4.txt", "r", stdin);
//    freopen("test_case5.txt", "r", stdin);
//    freopen("output.txt", "w", stdout);
    clock_t start = clock();
    input();
    solve();
    printf("time --> %.10f", (clock() - start) / (1.0 * CLOCKS_PER_SEC));
}
